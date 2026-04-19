# Project Guidelines

## Overview

This repository develops **neutron reflectometry (NR) models** for fitting experimental data from the **REF_L instrument** at the Spallation Neutron Source (SNS). Models are written in Python using the [refl1d](https://refl1d.readthedocs.io/) library and [bumps](https://bumps.readthedocs.io/) fitting framework.

## Repository Structure

Each top-level folder is an **experiment**, identified by researcher name, month, and year (e.g., `jen-apr2025/`, `jen-mar2026/`). Within each experiment:

- `data/` — Reduced reflectometry data files (`REFL_*_combined_data_auto.txt`, `REFL_*_partial.txt`)
- `models/` — Reusable model templates
- `Sample<N>/` — Per-sample directories containing sample-specific models and raw data
- `results/` — Fit outputs (JSON, CSV, plots)
- `notebooks/` — Jupyter notebooks for analysis and visualization

## Code Style

- All model files are standalone Python scripts that define a `problem = FitProblem(...)` at module scope (required by bumps/refl1d).
- Use `from refl1d.names import *` for refl1d API access.
- Prefer descriptive `SLD(name, rho=...)` material names that match the physical material.
- Use the pipe `|` operator for layer stacks, ordered **surface → substrate** (left to right).
- Silicon substrate (`Si`, rho=2.07) is always the rightmost (last) element; never vary its SLD.

## Key Domain Knowledge

See [docs/neutron-reflectometry.md](docs/neutron-reflectometry.md) for:
- Data file structure and header format
- Probe creation methods (QProbe vs make_probe)
- Common SLD values and range guidelines
- Chi-squared interpretation and BIC model complexity rules
- Roughness constraints and refinement strategies

## Build and Test

Models are run with bumps:
```bash
# Preview a model
bumps --preview model.py

# Fit a model using DREAM
bumps --fit=dream model.py --store=results/ --burn=1000 --steps=500
```

Requires Python with `refl1d`, `bumps`, and `numpy`. Activate the project venv:
```bash
source venv/bin/activate
```

## Conventions

- Data column order in combined files: `Q, R, dR, dQ` (columns 0, 1, 2, 3)
- `dQ` in data files is **FWHM**; convert to sigma with `dq /= 2.355` when using `QProbe`
- Multi-segment partial files follow the naming pattern: `REFL_{run}_{segment}_{subrun}_partial.txt`
- Standard angles for multi-angle measurements: `[0.45, 1.2, 3.5]` degrees (theta, not two-theta)
