#!/usr/bin/env python3
"""
Fit Quality Assessment Tool for Refl1D MCMC Results.

This script assesses the quality of reflectometry model fits produced by the
Refl1D/Bumps modeling engine using Markov Chain Monte Carlo (MCMC) methods.

Key quality metrics assessed:
1. Chi-squared (χ²) - Goodness of fit measure
2. Parameter uncertainties - From MCMC posterior distributions
3. Parameter correlations - Detection of correlated parameters
4. Residual analysis - Distribution and patterns in residuals
5. MCMC convergence - Assessment of chain mixing and sampling

Usage:
    python assess_fit_quality.py <results_directory>

Example:
    python assess_fit_quality.py results/expt11-ocv-08-corefine

Author: Fit Quality Assessment Tool
"""

import argparse
import json
import glob
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass
class ParameterStats:
    """Statistics for a single fitted parameter."""
    name: str
    best: float
    mean: float
    median: float
    std: float
    p68: Tuple[float, float]
    p95: Tuple[float, float]
    
    @property
    def relative_uncertainty(self) -> float:
        """Relative uncertainty (std/best) as percentage."""
        if abs(self.best) < 1e-10:
            return float('inf')
        return 100.0 * self.std / abs(self.best)
    
    @property
    def is_well_constrained(self) -> bool:
        """Parameter is well-constrained if relative uncertainty < 10%."""
        return self.relative_uncertainty < 10.0
    
    @property
    def mean_median_agreement(self) -> float:
        """Difference between mean and median relative to std."""
        if self.std < 1e-10:
            return 0.0
        return abs(self.mean - self.median) / self.std


@dataclass
class ModelStats:
    """Statistics for a single model in the fit."""
    model_index: int
    chisq: float
    chisq_err: Optional[float]
    nllf: float
    n_points: int
    
    @property
    def reduced_chisq(self) -> float:
        """Reduced chi-squared (chi-squared per degree of freedom)."""
        return self.chisq
    
    @property
    def is_good_fit(self) -> bool:
        """
        A good fit typically has reduced chi-squared close to 1.
        Values between 0.5 and 2.0 are generally acceptable.
        """
        return 0.5 <= self.reduced_chisq <= 2.0
    
    @property
    def fit_quality_description(self) -> str:
        """Human-readable description of fit quality."""
        if self.reduced_chisq < 0.5:
            return "Over-fitted (χ² too low, possibly over-parameterized or overestimated errors)"
        elif self.reduced_chisq <= 1.0:
            return "Excellent fit"
        elif self.reduced_chisq <= 1.5:
            return "Good fit"
        elif self.reduced_chisq <= 2.0:
            return "Acceptable fit"
        elif self.reduced_chisq <= 3.0:
            return "Marginal fit (may need model refinement)"
        else:
            return "Poor fit (model may be inadequate)"


@dataclass 
class ResidualStats:
    """Statistics for residual analysis."""
    mean: float
    std: float
    skewness: float
    kurtosis: float
    runs_test_z: float  # For detecting systematic patterns
    n_outliers_2sigma: int
    n_outliers_3sigma: int
    n_total: int
    
    @property
    def fraction_outliers_2sigma(self) -> float:
        return self.n_outliers_2sigma / self.n_total if self.n_total > 0 else 0
    
    @property
    def fraction_outliers_3sigma(self) -> float:
        return self.n_outliers_3sigma / self.n_total if self.n_total > 0 else 0
    
    @property
    def is_gaussian(self) -> bool:
        """Check if residuals are approximately Gaussian."""
        # For Gaussian: skewness ≈ 0, kurtosis ≈ 3
        return abs(self.skewness) < 0.5 and abs(self.kurtosis - 3) < 1.0
    
    @property
    def has_systematic_patterns(self) -> bool:
        """Check if residuals show systematic patterns (runs test)."""
        return abs(self.runs_test_z) > 2.0


@dataclass
class FitQualityReport:
    """Complete fit quality assessment report."""
    directory: str
    base_name: str
    n_models: int
    parameters: Dict[str, ParameterStats] = field(default_factory=dict)
    model_stats: List[ModelStats] = field(default_factory=list)
    residual_stats: List[ResidualStats] = field(default_factory=list)
    overall_chisq: Optional[float] = None
    overall_nllf: Optional[float] = None
    n_mcmc_samples: Optional[int] = None
    sample_fraction: Optional[float] = None
    warnings: List[str] = field(default_factory=list)
    
    @property
    def overall_quality_score(self) -> float:
        """
        Compute an overall quality score from 0-100.
        
        Factors considered:
        - Chi-squared values (40%)
        - Parameter constraints (30%)
        - Residual quality (30%)
        """
        score = 0.0
        
        # Chi-squared contribution (40 points max)
        if self.model_stats:
            chisq_scores = []
            for ms in self.model_stats:
                if ms.reduced_chisq < 0.5:
                    chisq_scores.append(30)  # Penalize over-fitting
                elif ms.reduced_chisq <= 1.0:
                    chisq_scores.append(40)
                elif ms.reduced_chisq <= 1.5:
                    chisq_scores.append(35)
                elif ms.reduced_chisq <= 2.0:
                    chisq_scores.append(25)
                elif ms.reduced_chisq <= 3.0:
                    chisq_scores.append(15)
                else:
                    chisq_scores.append(5)
            score += np.mean(chisq_scores)
        
        # Parameter constraint contribution (30 points max)
        if self.parameters:
            n_well_constrained = sum(1 for p in self.parameters.values() 
                                    if p.is_well_constrained)
            score += 30 * (n_well_constrained / len(self.parameters))
        
        # Residual quality contribution (30 points max)
        if self.residual_stats:
            residual_scores = []
            for rs in self.residual_stats:
                rs_score = 30
                if not rs.is_gaussian:
                    rs_score -= 10
                if rs.has_systematic_patterns:
                    rs_score -= 10
                if rs.fraction_outliers_3sigma > 0.01:  # More than 1% 3-sigma outliers
                    rs_score -= 5
                residual_scores.append(max(0, rs_score))
            score += np.mean(residual_scores)
        
        return min(100, max(0, score))
    
    @property
    def quality_grade(self) -> str:
        """Letter grade for overall quality."""
        score = self.overall_quality_score
        if score >= 90:
            return "A"
        elif score >= 80:
            return "B"
        elif score >= 70:
            return "C"
        elif score >= 60:
            return "D"
        else:
            return "F"


def calculate_skewness(data: np.ndarray) -> float:
    """Calculate Fisher-Pearson skewness coefficient."""
    n = len(data)
    if n < 3:
        return 0.0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std < 1e-10:
        return 0.0
    return np.mean(((data - mean) / std) ** 3) * n**2 / ((n-1) * (n-2))


def calculate_kurtosis(data: np.ndarray) -> float:
    """Calculate excess kurtosis (Fisher's definition)."""
    n = len(data)
    if n < 4:
        return 0.0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std < 1e-10:
        return 0.0
    m4 = np.mean((data - mean) ** 4)
    return (m4 / std**4) - 3


def runs_test(data: np.ndarray) -> float:
    """
    Perform Wald-Wolfowitz runs test for randomness.
    Returns z-score; |z| > 2 suggests non-randomness.
    """
    if len(data) < 10:
        return 0.0
    
    # Convert to binary based on median
    median = np.median(data)
    binary = data > median
    
    # Count runs
    runs = 1
    for i in range(1, len(binary)):
        if binary[i] != binary[i-1]:
            runs += 1
    
    # Calculate expected runs and variance
    n1 = np.sum(binary)
    n2 = len(binary) - n1
    
    if n1 == 0 or n2 == 0:
        return 0.0
    
    n = n1 + n2
    expected_runs = (2 * n1 * n2) / n + 1
    var_runs = (2 * n1 * n2 * (2 * n1 * n2 - n)) / (n**2 * (n - 1))
    
    if var_runs <= 0:
        return 0.0
    
    z = (runs - expected_runs) / np.sqrt(var_runs)
    return z


def load_error_json(filepath: str) -> Dict[str, ParameterStats]:
    """Load parameter statistics from -err.json file."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    params = {}
    for name, stats in data.items():
        params[name] = ParameterStats(
            name=name,
            best=stats['best'],
            mean=stats['mean'],
            median=stats['median'],
            std=stats['std'],
            p68=tuple(stats['p68']),
            p95=tuple(stats['p95'])
        )
    return params


def load_refl_data(filepath: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load reflectivity data from -refl.dat file."""
    data = np.loadtxt(filepath, comments='#')
    # Columns: Q, dQ, R, dR, theory, fresnel
    q = data[:, 0]
    r = data[:, 2]
    dr = data[:, 3]
    theory = data[:, 4]
    return q, r, dr, theory


def calculate_residuals(r: np.ndarray, dr: np.ndarray, 
                       theory: np.ndarray) -> np.ndarray:
    """Calculate normalized residuals: (data - theory) / error."""
    return (r - theory) / dr


def analyze_residuals(residuals: np.ndarray) -> ResidualStats:
    """Compute comprehensive residual statistics."""
    n = len(residuals)
    mean = np.mean(residuals)
    std = np.std(residuals)
    skewness = calculate_skewness(residuals)
    kurtosis = calculate_kurtosis(residuals) + 3  # Convert to regular kurtosis
    runs_z = runs_test(residuals)
    
    n_outliers_2sigma = np.sum(np.abs(residuals) > 2)
    n_outliers_3sigma = np.sum(np.abs(residuals) > 3)
    
    return ResidualStats(
        mean=mean,
        std=std,
        skewness=skewness,
        kurtosis=kurtosis,
        runs_test_z=runs_z,
        n_outliers_2sigma=n_outliers_2sigma,
        n_outliers_3sigma=n_outliers_3sigma,
        n_total=n
    )


def parse_output_file(filepath: str) -> Tuple[List[ModelStats], float, float]:
    """Parse the .out file to extract chi-squared and model statistics."""
    model_stats = []
    overall_chisq = None
    overall_nllf = None
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Parse individual model chi-squared values
    # Format: [chisq=2.4477(65), nllf=892.187]
    model_pattern = r'\[chisq=(\d+\.?\d*)\((\d+)\),\s*nllf=(\d+\.?\d*)\]'
    
    for i, match in enumerate(re.finditer(model_pattern, content)):
        chisq = float(match.group(1))
        chisq_err = int(match.group(2)) / 10**(len(match.group(2)))  # Convert error format
        nllf = float(match.group(3))
        
        model_stats.append(ModelStats(
            model_index=i,
            chisq=chisq,
            chisq_err=chisq_err,
            nllf=nllf,
            n_points=0  # Will be updated from data files
        ))
    
    # Parse overall chi-squared
    overall_pattern = r'\[overall chisq=(\d+\.?\d*)\((\d+)\),\s*nllf=(\d+\.?\d*)\]'
    match = re.search(overall_pattern, content)
    if match:
        overall_chisq = float(match.group(1))
        overall_nllf = float(match.group(3))
    
    return model_stats, overall_chisq, overall_nllf


def parse_err_file(filepath: str) -> Tuple[int, float]:
    """Parse the .err file to extract MCMC sample information."""
    n_samples = None
    sample_fraction = None
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Format: Statistics and plots based on 10500 samples (100% of total samples drawn)
    pattern = r'based on (\d+) samples \((\d+)%'
    match = re.search(pattern, content)
    if match:
        n_samples = int(match.group(1))
        sample_fraction = float(match.group(2)) / 100.0
    
    return n_samples, sample_fraction


def find_result_files(directory: str) -> Tuple[str, List[str]]:
    """
    Find result files in the directory and determine base name.
    Returns (base_name, list of experiment indices).
    """
    # Find -err.json files to determine base name
    err_json_files = glob.glob(os.path.join(directory, '*-err.json'))
    
    if not err_json_files:
        raise FileNotFoundError(f"No -err.json files found in {directory}")
    
    # Extract base name (everything before -err.json)
    base_name = os.path.basename(err_json_files[0]).replace('-err.json', '')
    
    # Find experiment data files (-N-refl.dat)
    refl_files = glob.glob(os.path.join(directory, f'{base_name}-*-refl.dat'))
    
    # Extract indices
    indices = []
    for f in refl_files:
        fname = os.path.basename(f)
        # Pattern: basename-N-refl.dat
        match = re.search(rf'{re.escape(base_name)}-(\d+)-refl\.dat', fname)
        if match:
            indices.append(match.group(1))
    
    indices.sort(key=int)
    return base_name, indices


def assess_fit_quality(directory: str) -> FitQualityReport:
    """
    Perform comprehensive fit quality assessment on a results directory.
    
    Parameters
    ----------
    directory : str
        Path to the results directory containing fit output files.
    
    Returns
    -------
    FitQualityReport
        Comprehensive assessment of fit quality.
    """
    directory = os.path.abspath(directory)
    
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"Not a directory: {directory}")
    
    base_name, model_indices = find_result_files(directory)
    
    report = FitQualityReport(
        directory=directory,
        base_name=base_name,
        n_models=len(model_indices)
    )
    
    # Load parameter statistics
    err_json_path = os.path.join(directory, f'{base_name}-err.json')
    if os.path.exists(err_json_path):
        report.parameters = load_error_json(err_json_path)
    
    # Parse .out file for chi-squared values
    out_path = os.path.join(directory, f'{base_name}.out')
    if os.path.exists(out_path):
        report.model_stats, report.overall_chisq, report.overall_nllf = \
            parse_output_file(out_path)
    
    # Parse .err file for MCMC sampling info
    err_path = os.path.join(directory, f'{base_name}.err')
    if os.path.exists(err_path):
        report.n_mcmc_samples, report.sample_fraction = parse_err_file(err_path)
    
    # Analyze residuals for each model
    for idx in model_indices:
        refl_path = os.path.join(directory, f'{base_name}-{idx}-refl.dat')
        if os.path.exists(refl_path):
            q, r, dr, theory = load_refl_data(refl_path)
            residuals = calculate_residuals(r, dr, theory)
            residual_stats = analyze_residuals(residuals)
            report.residual_stats.append(residual_stats)
            
            # Update n_points in model stats if available
            model_idx = int(idx) - 1
            if model_idx < len(report.model_stats):
                report.model_stats[model_idx].n_points = len(q)
    
    # Generate warnings
    report.warnings = generate_warnings(report)
    
    return report


def generate_warnings(report: FitQualityReport) -> List[str]:
    """Generate warning messages based on the assessment."""
    warnings = []
    
    # Chi-squared warnings
    for ms in report.model_stats:
        if ms.reduced_chisq > 3.0:
            warnings.append(
                f"Model {ms.model_index + 1}: High χ² ({ms.reduced_chisq:.2f}) "
                f"suggests poor fit or underestimated errors"
            )
        elif ms.reduced_chisq < 0.5:
            warnings.append(
                f"Model {ms.model_index + 1}: Low χ² ({ms.reduced_chisq:.2f}) "
                f"suggests over-fitting or overestimated errors"
            )
    
    # Parameter constraint warnings
    for name, param in report.parameters.items():
        if not param.is_well_constrained:
            warnings.append(
                f"Parameter '{name}': Poorly constrained "
                f"(relative uncertainty {param.relative_uncertainty:.1f}%)"
            )
        if param.mean_median_agreement > 0.5:
            warnings.append(
                f"Parameter '{name}': Asymmetric distribution "
                f"(mean-median difference = {param.mean_median_agreement:.2f}σ)"
            )
    
    # Residual warnings
    for i, rs in enumerate(report.residual_stats):
        if rs.has_systematic_patterns:
            warnings.append(
                f"Model {i + 1}: Residuals show systematic patterns "
                f"(runs test z = {rs.runs_test_z:.2f})"
            )
        if rs.fraction_outliers_3sigma > 0.01:
            warnings.append(
                f"Model {i + 1}: {rs.fraction_outliers_3sigma*100:.1f}% outliers "
                f"(> 3σ, expected ~0.3%)"
            )
        if not rs.is_gaussian:
            warnings.append(
                f"Model {i + 1}: Non-Gaussian residuals "
                f"(skew={rs.skewness:.2f}, kurtosis={rs.kurtosis:.2f})"
            )
    
    # MCMC sampling warnings
    if report.sample_fraction is not None and report.sample_fraction < 0.5:
        warnings.append(
            f"Only {report.sample_fraction*100:.0f}% of MCMC samples used "
            f"for statistics (may indicate convergence issues)"
        )
    
    return warnings


def format_report(report: FitQualityReport) -> str:
    """Format the assessment report as a human-readable string."""
    lines = []
    lines.append("=" * 70)
    lines.append("FIT QUALITY ASSESSMENT REPORT")
    lines.append("=" * 70)
    lines.append(f"\nDirectory: {report.directory}")
    lines.append(f"Base name: {report.base_name}")
    lines.append(f"Number of models: {report.n_models}")
    
    # Overall quality
    lines.append(f"\n{'─' * 70}")
    lines.append("OVERALL QUALITY")
    lines.append(f"{'─' * 70}")
    lines.append(f"Quality Score: {report.overall_quality_score:.1f}/100 (Grade: {report.quality_grade})")
    if report.overall_chisq is not None:
        lines.append(f"Overall χ²: {report.overall_chisq:.4f}")
    if report.overall_nllf is not None:
        lines.append(f"Overall -log(likelihood): {report.overall_nllf:.2f}")
    if report.n_mcmc_samples is not None:
        lines.append(f"MCMC samples used: {report.n_mcmc_samples} ({report.sample_fraction*100:.0f}%)")
    
    # Model-by-model chi-squared
    if report.model_stats:
        lines.append(f"\n{'─' * 70}")
        lines.append("MODEL FIT QUALITY")
        lines.append(f"{'─' * 70}")
        for ms in report.model_stats:
            lines.append(f"\nModel {ms.model_index + 1}:")
            lines.append(f"  χ² = {ms.chisq:.4f} ± {ms.chisq_err:.4f}" if ms.chisq_err else f"  χ² = {ms.chisq:.4f}")
            lines.append(f"  -log(likelihood) = {ms.nllf:.2f}")
            lines.append(f"  Assessment: {ms.fit_quality_description}")
    
    # Parameter statistics
    if report.parameters:
        lines.append(f"\n{'─' * 70}")
        lines.append("PARAMETER STATISTICS")
        lines.append(f"{'─' * 70}")
        lines.append(f"{'Parameter':<25} {'Best':>12} {'Mean':>12} {'Std':>10} {'Rel.Unc.':<10}")
        lines.append("-" * 70)
        for name, p in report.parameters.items():
            rel_unc = f"{p.relative_uncertainty:.1f}%" if p.relative_uncertainty < 100 else ">100%"
            constrained = "✓" if p.is_well_constrained else "✗"
            lines.append(f"{name:<25} {p.best:>12.4f} {p.mean:>12.4f} {p.std:>10.4f} {rel_unc:<8} {constrained}")
    
    # Residual analysis
    if report.residual_stats:
        lines.append(f"\n{'─' * 70}")
        lines.append("RESIDUAL ANALYSIS")
        lines.append(f"{'─' * 70}")
        for i, rs in enumerate(report.residual_stats):
            lines.append(f"\nModel {i + 1} Residuals:")
            lines.append(f"  Mean: {rs.mean:.4f} (ideal: 0)")
            lines.append(f"  Std: {rs.std:.4f} (ideal: 1)")
            lines.append(f"  Skewness: {rs.skewness:.4f} (ideal: 0)")
            lines.append(f"  Kurtosis: {rs.kurtosis:.4f} (ideal: 3)")
            lines.append(f"  Runs test z-score: {rs.runs_test_z:.4f} (|z| > 2 suggests patterns)")
            lines.append(f"  Outliers: {rs.n_outliers_2sigma} (>2σ), {rs.n_outliers_3sigma} (>3σ) of {rs.n_total}")
            lines.append(f"  Gaussian: {'Yes' if rs.is_gaussian else 'No'}")
            lines.append(f"  Systematic patterns: {'Yes' if rs.has_systematic_patterns else 'No'}")
    
    # Warnings
    if report.warnings:
        lines.append(f"\n{'─' * 70}")
        lines.append("WARNINGS")
        lines.append(f"{'─' * 70}")
        for warning in report.warnings:
            lines.append(f"  ⚠ {warning}")
    
    # Summary and recommendations
    lines.append(f"\n{'─' * 70}")
    lines.append("SUMMARY")
    lines.append(f"{'─' * 70}")
    
    if report.quality_grade in ['A', 'B']:
        lines.append("✓ The fit quality is good. Results can be used with confidence.")
    elif report.quality_grade == 'C':
        lines.append("△ The fit quality is acceptable but could be improved.")
        lines.append("  Consider reviewing the warnings above.")
    else:
        lines.append("✗ The fit quality needs improvement.")
        lines.append("  Recommendations:")
        if any(ms.reduced_chisq > 3.0 for ms in report.model_stats):
            lines.append("  - Review the physical model for missing features")
            lines.append("  - Check if error bars are underestimated")
        if any(not p.is_well_constrained for p in report.parameters.values()):
            lines.append("  - Consider fixing poorly constrained parameters")
            lines.append("  - Add more constraints or prior information")
        if any(rs.has_systematic_patterns for rs in report.residual_stats):
            lines.append("  - Investigate systematic patterns in residuals")
            lines.append("  - Check for missing structural features")
    
    lines.append("\n" + "=" * 70)
    
    return "\n".join(lines)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Assess the quality of Refl1D MCMC fit results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s results/expt11-ocv-08-corefine
  %(prog)s --json results/expt11-ocv-08-corefine
  
The tool analyzes:
  - Chi-squared values and fit quality
  - Parameter uncertainties and constraints
  - Residual distributions and patterns
  - MCMC convergence indicators

Quality grades:
  A (90-100): Excellent fit
  B (80-89):  Good fit
  C (70-79):  Acceptable fit
  D (60-69):  Marginal fit
  F (< 60):   Poor fit
        """
    )
    parser.add_argument(
        "directory",
        help="Path to the results directory containing fit output files"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results in JSON format"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only output the quality grade and score"
    )
    
    args = parser.parse_args()
    
    try:
        report = assess_fit_quality(args.directory)
        
        if args.json:
            # Output JSON summary
            output = {
                "directory": report.directory,
                "base_name": report.base_name,
                "quality_score": report.overall_quality_score,
                "quality_grade": report.quality_grade,
                "overall_chisq": report.overall_chisq,
                "n_models": report.n_models,
                "n_parameters": len(report.parameters),
                "n_warnings": len(report.warnings),
                "warnings": report.warnings,
                "parameters": {
                    name: {
                        "best": p.best,
                        "mean": p.mean,
                        "std": p.std,
                        "relative_uncertainty_pct": p.relative_uncertainty,
                        "well_constrained": p.is_well_constrained
                    }
                    for name, p in report.parameters.items()
                },
                "models": [
                    {
                        "index": ms.model_index,
                        "chisq": ms.chisq,
                        "is_good_fit": ms.is_good_fit,
                        "description": ms.fit_quality_description
                    }
                    for ms in report.model_stats
                ]
            }
            print(json.dumps(output, indent=2))
        elif args.quiet:
            print(f"{report.quality_grade} ({report.overall_quality_score:.1f})")
        else:
            print(format_report(report))
        
        # Exit with non-zero status if fit quality is poor
        if report.quality_grade == 'F':
            sys.exit(1)
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        print(f"Error analyzing results: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
