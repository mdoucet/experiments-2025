---
description: >
  Use when editing or creating refl1d model Python files (.py) in experiment directories.
  Enforces layer stack conventions, SLD API rules, QProbe/make_probe patterns,
  parameter range guidelines, and FitProblem requirements for bumps/refl1d models.
applyTo: "**/jen-*/**/*.py"
---

# Refl1d Model File Conventions

## File Structure Requirements

Every model file MUST define `problem = FitProblem(...)` at module scope.
This is required by the bumps/refl1d fitting framework.

## Layer Stack Order

Stacks are ordered **surface → substrate** using the pipe `|` operator.
Silicon substrate is always the rightmost element.

```python
sample = Solvent(0, roughness) | ... | Cu(thickness, roughness) | Ti(thickness, roughness) | Si
```

## API Safety

Set parameter bounds on `sample['LayerName']` (Slab objects), NEVER on `SLD(...)` objects:
```python
sample['Cu'].material.rho.range(5.0, 7.0)   # Correct
# Cu.material.rho.range(5.0, 7.0)           # WRONG — crashes
```

## Data Loading

- Combined files: use `QProbe`, convert dQ from FWHM to sigma with `dq /= 2.355`
- Partial files: use `make_probe` with explicit theta angle
- Column order: `Q, R, dR, dQ` (indices 0, 1, 2, 3)

## Substrate

- Silicon SLD is always 2.07 — never add `.range()` or `.pm()` to it
- Ti adhesion layer SLD range: at least -4.0 to 0.0
