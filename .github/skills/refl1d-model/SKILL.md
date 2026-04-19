---
name: refl1d-model
description: >
  Build a single-measurement refl1d neutron reflectometry model script.
  Use when: creating a new model file, writing a bumps FitProblem, defining SLD materials
  and layer stacks, setting parameter ranges, loading combined or partial data files,
  creating QProbe or make_probe probes. Covers Cu, Pt, Ti, ionomer, and oxide layer stacks.
---

# Single-Measurement Refl1d Model

## When to Use

- Creating a new model `.py` file for a single reflectometry measurement
- Setting up materials, layer stacks, and parameter bounds for fitting
- Loading data from combined or partial files
- Choosing between QProbe (combined data) and make_probe (multi-angle partial data)

## Procedure

### 1. Choose the Data Loading Strategy

There are two probe types. The data file determines which to use.

**QProbe — for combined data files** (`REFL_{run}_combined_data_auto.txt`):
```python
from refl1d.names import *
import numpy as np

_refl = np.loadtxt(data_file).T
q, R, dR, dq = _refl[0], _refl[1], _refl[2], _refl[3]

# CRITICAL: dQ in data files is FWHM — convert to sigma
dq /= 2.355

probe = QProbe(q, dq, data=(R, dR))
probe.intensity = Parameter(value=1, name='intensity')
probe.intensity.pm(0.05)
```

**make_probe — for multi-angle partial files** (`REFL_{run}_{seg}_{subrun}_partial.txt`):
```python
from refl1d.names import *
from refl1d.probe import make_probe
import numpy as np

def create_probe(data_file, theta):
    """Create a probe from a single-angle partial file.
    
    Args:
        data_file: Path to partial data file
        theta: Incident angle in degrees (NOT two-theta)
    """
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0.0 * q  # Wavelength resolution (set to 0 or use moderator model)
    
    probe = make_probe(
        T=theta, dT=dT, L=wl, dL=dL,
        data=(data, errors),
        radiation="neutron",
        resolution="uniform",
    )
    return probe
```

Standard angles for REF_L: `thetas = [0.45, 1.2, 3.5]` degrees.

The TwoTheta value from the data file header is **2× theta** — divide by 2 before using.

### 2. Define Materials

Use `SLD(name, rho=value)` for each material. Common materials and typical SLD values (×10⁻⁶ Å⁻²):

| Material | Typical SLD | Notes |
|----------|-------------|-------|
| Si | 2.07 | Substrate — **never vary** |
| SiOx | 2.8–3.47 | Native oxide — often omitted, fix thickness < 20 Å if present |
| Ti | -1.95 | Adhesion layer — can intermix, use wide range |
| Cu | 6.3–6.55 | Working electrode |
| Pt | 6.288 | Working electrode — constrain tightly (pm 0.07) |
| CuOx | 4.9–5.6 | Copper oxide |
| PtOx | 5.0–6.0 | Platinum oxide |
| D2O | 6.19–6.3 | Deuterated solvent |
| H2O | 5.5–6.2 | Hydrogenated water (in D2O contrast) |
| THF | 5.8–6.0 | Tetrahydrofuran solvent |
| Air | 0.0 | For air measurements |
| ionomer | 0.6–4.8 | Polymer coating — wide range, poorly constrained |

```python
Si = SLD('Si', rho=2.07)
Ti = SLD('Ti', rho=-1.95)
Cu = SLD('Cu', rho=6.4)
D2O = SLD('D2O', rho=6.19)
```

### 3. Build the Layer Stack

Use the pipe `|` operator. 

**Order depends on the measurement geometry**:
- When the incident beam is coming from the substrate side (back reflection), the order is **surface → substrate** (left to right). Silicon is always last.
- When the beam is coming from the ambient (air or solvent) side, the order is reversed (substrate → surface). In that case, silicon is always first.

```python
# Bare metal in solvent
sample = D2O(0, 10) | Cu(500, 5) | Ti(35, 5) | Si

# Metal with oxide in solvent
sample = D2O(0, 10) | CuOx(30, 10) | Cu(500, 5) | Ti(35, 5) | Si

# Metal with coating and oxide
sample = H2O(0, 15) | ionomer(1000, 20) | CuOx(27, 3) | Cu(100, 5) | Ti(32, 1.3) | Si
```

The `SLD(thickness, roughness)` call creates a `Slab`. The first argument after the SLD name is thickness (Å), the second is interface roughness (Å).

### 4. Set Parameter Ranges

**CRITICAL**: Set ranges on `sample['LayerName']` (Slab), NOT on the SLD object.

```python
# Correct — set on sample slab
sample['Cu'].material.rho.range(5.0, 7.0)
sample['Cu'].thickness.range(400, 600)
sample['Cu'].interface.range(1.0, 12.0)

# WRONG — crashes with "'SLD' object has no attribute 'material'"
# Cu.material.rho.range(5.0, 7.0)
```

Use `.range(min, max)` for wide bounds or `.pm(delta)` for tight symmetric bounds.

Typical parameter ranges:

| Layer | Thickness | SLD range | Interface |
|-------|-----------|-----------|-----------|
| Solvent | 0 (fixed) | ±0.5–1.0 around nominal | 1–50 |
| Ionomer/coating | 10–2000 | -1 to 6.5 | 3–33 |
| Oxide (CuOx/PtOx) | 5–200 | -1 to 7 (CuOx), 5.7–6.2 (PtOx) | 3–25 |
| Cu | 50–750 | 5–7 | 1–25 |
| Pt | pm(0.2) | pm(0.07) | 5–10 |
| Ti | 10–60 | -4 to 2.1 | 1–22 |

### 5. Create the Experiment and FitProblem

```python
experiment = Experiment(sample=sample, probe=probe)
problem = FitProblem(experiment)
```

The `problem = FitProblem(...)` line **must** exist at module scope — bumps/refl1d requires it.

### 6. Complete Model Template

```python
""" Model description. Copy the user instructions template and fill in details. """
import os
import numpy as np
from refl1d.names import *


def create_fit_experiment(q, dq, data, errors):
    # Convert FWHM to sigma
    dq /= 2.355

    probe = QProbe(q, dq, data=(data, errors))
    probe.intensity = Parameter(value=1, name='intensity')
    probe.intensity.pm(0.05)

    # Define materials
    D2O = SLD('D2O', rho=6.19)
    Si = SLD('Si', rho=2.07)
    Ti = SLD('Ti', rho=-1.95)
    Cu = SLD('Cu', rho=6.4)

    # Build stack: surface | ... | substrate
    sample = D2O(0, 10) | Cu(500, 5) | Ti(35, 5) | Si

    M = Experiment(sample=sample, probe=probe)

    # Set parameter ranges
    sample['D2O'].material.rho.pm(0.5)
    sample['D2O'].interface.range(1, 25)

    sample['Cu'].thickness.range(400, 600)
    sample['Cu'].material.rho.range(5.0, 7.0)
    sample['Cu'].interface.range(1.0, 12.0)

    sample['Ti'].thickness.range(20.0, 50.0)
    sample['Ti'].material.rho.range(-4.0, 0.0)
    sample['Ti'].interface.range(1.0, 15.0)

    return M


# Data path — use relative path from repo root or absolute path to data directory
data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
data_file = os.path.join(data_dir, 'REFL_XXXXXX_combined_data_auto.txt')

_refl = np.loadtxt(data_file).T
experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])
problem = FitProblem(experiment)
```

## Running the Model

```bash
# Activate environment
source venv/bin/activate

# Preview (check layer stack and initial reflectivity)
bumps --preview model.py

# Fit with DREAM algorithm
bumps --fit=dream model.py --store=results/ --burn=1000 --steps=500
```

## References

- [Domain knowledge and refinement strategies](./references/domain-knowledge.md)
- [Multi-angle partial file probe creation](./references/multi-angle-probe.md)
