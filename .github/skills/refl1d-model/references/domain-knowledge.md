## Data Structure for the REF_L Instrument

The Liquids Reflectometer (REF_L) at the Spallation Neutron Source (SNS) produces data files
with segments corresponding to different angle/wavelength configurations.

### Combined Data File

Single file with all segments merged. Columns: `Q, R, dR, dQ` (dQ is FWHM).

Header example:
```
# Experiment IPTS-36897 Run 226613
# DataRun   NormRun   TwoTheta(deg)  LambdaMin(A)  LambdaMax(A)
# 226613    226559    0.739463       2.74975       9.4987
# 226614    226560    2.39957        2.74978       9.49867
# 226615    226561    7.00027        2.74977       9.49868
# Q [1/Angstrom]  R  dR  dQ [FWHM]
```

### Partial Files (Multi-Segment)

Individual angle files: `REFL_{run}_{segment}_{subrun}_partial.txt`
- Segment 1: `REFL_226613_1_226613_partial.txt` → theta = TwoTheta/2
- Segment 2: `REFL_226613_2_226614_partial.txt`
- Segment 3: `REFL_226613_3_226615_partial.txt`

Standard angles: `[0.45, 1.2, 3.5]` degrees (theta, NOT two-theta).

## Common SLD Values (×10⁻⁶ Å⁻²)

| Material | SLD | Notes |
|----------|-----|-------|
| Silicon | 2.07 | Substrate — never vary |
| SiO₂ | 3.47 | Native oxide — omit by default |
| Air | 0.0 | |
| Gold | 4.5 | |
| Copper | 6.55 | |
| Titanium | -1.95 | Adhesion layer |
| Platinum | 6.288 | |
| D₂O | 6.19–6.3 | |
| THF | 5.8–6.0 | |

### SLD Range Guidelines

- Set ranges to at least ±2.0 around nominal value
- Never use ranges narrower than ±1.0
- For adhesion layers like Ti that can intermix: use ±3.0 or wider (e.g., -5.0 to 1.0)

## Chi-Squared Interpretation

| χ² | Interpretation |
|----|----------------|
| ≈ 1 | Ideal fit |
| < 0.5 | Possible overfitting |
| 1–2 | Excellent fit |
| 2–5 | Good fit |
| 5–10 | Marginal fit |
| > 10 | Poor fit |

## Model Complexity (BIC)

- BIC = n·ln(χ²) + k·ln(n), where n = data points, k = free parameters
- Each layer adds 3 free parameters (thickness, SLD, roughness)
- Do NOT add layers unless BIC would clearly improve
- Do NOT split layers (e.g., CuO + Cu₂O) unless χ² > 10 with residual evidence
- Never make multiple structural changes at once

## Roughness Constraints

- Minimum roughness: 5 Å
- Must be < half the thickness of either adjacent layer
- Typical range: 5–30 Å

## Refl1d API Rules

**CRITICAL**: `SLD(...)` objects do NOT have `.material`, `.thickness`, or `.interface`
attributes. Those only exist on `Slab` objects inside the sample stack.

```python
# CORRECT
sample['Cu'].material.rho.range(5.5, 7.0)

# WRONG — crashes
Cu.material.rho.range(5.5, 7.0)
```

## General Constraints

- Never change the fitting engine/method
- Never reverse the layer order or change back-reflection geometry
- Never change error bars, resolution, or Q-range
- Never vary the substrate SLD (unless explicitly requested)
- Minimum layer thickness: 5 Å

## Refinement Strategy

Priority order when χ² is above threshold:

1. Constrain unphysical parameters (values far from nominal)
2. Widen bounds on parameters hitting limits
3. Adjust starting values to best-fit from previous iteration
4. Check ambient SLD (common source of high χ²)
5. Enable sample_broadening for multi-segment data
6. Structural changes as last resort (only if χ² > 10 with residual evidence)
