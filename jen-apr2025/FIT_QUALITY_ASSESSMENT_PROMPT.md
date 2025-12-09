# Fit Quality Assessment - AI Prompt Template

Use the following prompt to assess reflectometry fit quality using the `assess_fit_quality.py` script:

---

## Prompt for AI Assessment

```
I need you to assess the quality of a reflectometry MCMC fit. Please run the fit quality assessment tool and interpret the results.

1. Run the assessment:
   python assess_fit_quality.py <results_directory>

2. Analyze the output and provide:
   - A summary of the overall fit quality (grade and score)
   - Key issues identified in the fit
   - Specific recommendations for improvement
   - Whether the results are reliable enough for publication

3. For each warning, explain:
   - What the issue means physically
   - How it might affect the interpretation
   - Potential remedies

Focus on these key quality indicators:
- Chi-squared values (ideal: 1.0, acceptable: 0.5-2.0)
- Parameter uncertainties (well-constrained: <10% relative uncertainty)
- Residual patterns (systematic patterns indicate missing physics)
- Outlier fraction (>3σ outliers should be <1%)
```

---

## Example Usage

```bash
# Basic assessment
python assess_fit_quality.py results/expt11-ocv-08-corefine

# JSON output for programmatic use
python assess_fit_quality.py results/expt11-ocv-08-corefine --json

# Quick grade only
python assess_fit_quality.py results/expt11-ocv-08-corefine --quiet
```

---

## Interpreting Results

### Quality Grades
- **A (90-100)**: Excellent fit - results can be published with high confidence
- **B (80-89)**: Good fit - results are reliable, minor improvements possible
- **C (70-79)**: Acceptable fit - use with caution, consider improvements
- **D (60-69)**: Marginal fit - significant issues, investigate before using
- **F (< 60)**: Poor fit - results are unreliable, major revision needed

### Key Metrics

| Metric | Ideal | Acceptable | Problematic |
|--------|-------|------------|-------------|
| χ² (reduced) | 1.0 | 0.5 - 2.0 | < 0.5 or > 3.0 |
| Relative uncertainty | < 5% | < 10% | > 10% |
| Residual std | 1.0 | 0.8 - 1.2 | > 1.5 |
| 3σ outliers | 0.3% | < 1% | > 2% |
| Runs test |z| | < 1 | < 2 | > 2 |

### Common Issues and Solutions

1. **High χ²** (> 2.0)
   - Model missing structural features
   - Error bars underestimated
   - Data quality issues
   
2. **Low χ²** (< 0.5)
   - Over-parameterized model
   - Error bars overestimated
   - Correlated parameters

3. **Poorly constrained parameters** (> 10% uncertainty)
   - Parameter not well-determined by data
   - Consider fixing to known value
   - May need different Q-range or contrast

4. **Systematic residual patterns**
   - Missing layer or feature in model
   - Incorrect interfacial roughness model
   - Background or resolution issues

5. **Non-Gaussian residuals**
   - May indicate data quality issues
   - Could suggest model inadequacy
   - Check for systematic errors

---

## File Requirements

The script expects a results directory containing:
- `*-err.json` - Parameter statistics from MCMC
- `*-N-refl.dat` - Reflectivity data and theory for each model
- `*.out` - Summary output with chi-squared values
- `*.err` - MCMC sampling information

These files are standard outputs from Refl1D/Bumps MCMC fitting.
