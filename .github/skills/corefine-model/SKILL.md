---
name: corefine-model
description: >
  Build a co-refinement refl1d model that fits multiple measurements together with
  shared constraints. Use when: co-refining multiple datasets, constraining parameters
  across experiments, fitting electrochemical sequences (OCV, applied potential),
  multi-state measurements, tying Cu/Ti/Pt substrate parameters, allowing surface
  layers to vary independently between states.
---

# Co-refinement Model

## When to Use

- Fitting 2+ measurements of the same sample under different conditions
- Constraining shared substrate parameters (Ti, Cu, Pt) across experiments
- Modeling an electrochemical sequence (OCV → applied potential → OCV recovery)
- Any case where physical parameters must be tied between datasets

## Procedure

### 1. Decide What to Share vs Free

In co-refinement, some parameters are physically identical across measurements
and must be constrained, while others change between states.

**Typically shared (tied) across all measurements:**
- Substrate properties: Ti thickness, Ti SLD, Ti interface
- Base metal intrinsic SLD: Cu or Pt SLD
- Base metal interface roughness (often)
- Solvent SLD (if same solvent throughout)

**Typically free (independent) per measurement:**
- Surface/coating thickness (e.g., oxide, SEI, ionomer)
- Surface/coating SLD and roughness
- Probe intensity
- Base metal thickness (can change due to dissolution/deposition)
- Solvent interface roughness

### 2. Use a Factory Function

Define a function that creates one experiment. Call it for each dataset:

```python
def create_fit_experiment(q, dq, data, errors):
    dq /= 2.355  # FWHM to sigma

    probe = QProbe(q, dq, data=(data, errors))
    probe.intensity = Parameter(value=1, name='intensity')
    probe.intensity.pm(0.05)

    D2O = SLD('D2O', rho=6.19)
    Si = SLD('Si', rho=2.07)
    Ti = SLD('Ti', rho=-2.0)
    Cu = SLD('Cu', rho=6.4)
    CuOx = SLD(name='CuOx', rho=5.0, irho=0.0)

    sample = D2O(0, 10) | CuOx(30, 10) | Cu(500, 5) | Ti(35, 5) | Si

    M = Experiment(sample=sample, probe=probe)

    # Set ranges for ALL parameters
    sample['D2O'].material.rho.pm(0.5)
    sample['D2O'].interface.range(1, 25)
    sample['CuOx'].thickness.range(5, 200)
    sample['CuOx'].material.rho.range(4.0, 7.0)
    sample['CuOx'].interface.range(3, 25)
    sample['Cu'].thickness.range(400, 600)
    sample['Cu'].material.rho.range(5.0, 7.0)
    sample['Cu'].interface.range(1.0, 12.0)
    sample['Ti'].thickness.range(20.0, 50.0)
    sample['Ti'].material.rho.range(-4.0, 0.0)
    sample['Ti'].interface.range(1.0, 15.0)

    return M
```

### 3. Load Data and Create Experiments

```python
data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')

# Load each measurement
_refl1 = np.loadtxt(os.path.join(data_dir, 'REFL_XXXXXX_combined_data_auto.txt')).T
experiment1 = create_fit_experiment(_refl1[0], _refl1[3], _refl1[1], _refl1[2])

_refl2 = np.loadtxt(os.path.join(data_dir, 'REFL_YYYYYY_combined_data_auto.txt')).T
experiment2 = create_fit_experiment(_refl2[0], _refl2[3], _refl2[1], _refl2[2])
```

### 4. Apply Constraints

Tie shared parameters by assigning one experiment's parameter to another:

```python
# Tie substrate parameters
experiment2.sample['Ti'].material.rho = experiment1.sample['Ti'].material.rho
experiment2.sample['Ti'].thickness = experiment1.sample['Ti'].thickness
experiment2.sample['Ti'].interface = experiment1.sample['Ti'].interface

# Tie base metal SLD and interface
experiment2.sample['Cu'].material.rho = experiment1.sample['Cu'].material.rho
experiment2.sample['Cu'].interface = experiment1.sample['Cu'].interface

# Tie solvent SLD (if same solvent)
experiment2.sample['D2O'].material.rho = experiment1.sample['D2O'].material.rho
```

**CRITICAL**: The assignment direction matters. Always assign TO the later experiment
FROM the reference experiment. The right-hand side becomes the master parameter.

### 5. Create FitProblem with List

```python
problem = FitProblem([experiment1, experiment2])
```

### 6. For Three or More Measurements

Chain constraints from the first experiment:

```python
for exp in [experiment2, experiment3]:
    exp.sample['Ti'].material.rho = experiment1.sample['Ti'].material.rho
    exp.sample['Ti'].thickness = experiment1.sample['Ti'].thickness
    exp.sample['Ti'].interface = experiment1.sample['Ti'].interface
    exp.sample['Cu'].material.rho = experiment1.sample['Cu'].material.rho
    exp.sample['Cu'].interface = experiment1.sample['Cu'].interface

problem = FitProblem([experiment1, experiment2, experiment3])
```

## Multi-Angle Co-refinement

When each measurement has multiple angle segments (partial files), create one
experiment per angle, share ALL structural parameters within an angle set,
and allow independent probe intensity per angle.

See the [multi-angle probe reference](../refl1d-model/references/multi-angle-probe.md)
for probe creation details.

### Single dataset, multiple angles

```python
thetas = [0.45, 1.2, 3.5]
probes, samples, experiments = [], [], []

for i, theta in enumerate(thetas):
    _probe = create_probe(data_files[i], theta)
    _sample = create_sample()

    probes.append(_probe)
    samples.append(_sample)
    experiments.append(Experiment(sample=_sample, probe=_probe))

    if i == 0:
        _probe.intensity = Parameter(1, name="Intensity 1")
        _probe.intensity.pm(0.15)
        _probe.sample_broadening.range(0, 0.05)
    else:
        _probe.intensity = Parameter(1, name="Intensity %d" % (i + 1))
        _probe.intensity.pm(0.15)
        _probe.sample_broadening = probes[0].sample_broadening
        _probe.theta_offset = probes[0].theta_offset

        # Share ALL structural parameters within same measurement
        for layer_name in ['Ti', 'Cu', 'CuOx', 'ionomer', 'H2O']:
            if layer_name in [s.name for s in _sample]:
                _sample[layer_name].material.rho = samples[0][layer_name].material.rho
                _sample[layer_name].thickness = samples[0][layer_name].thickness
                _sample[layer_name].interface = samples[0][layer_name].interface

problem = FitProblem(experiments)
```

### Multiple datasets × multiple angles (cross-dataset constraints)

**CRITICAL pitfall**: Within each dataset the within-angle constraint loop links
angles 1 and 2 to angle 0's parameter *objects*. When you later reassign a
parameter on angle 0's sample (e.g., `samples2[0]["Cu"].material.rho = ...`),
only that one reference changes — angles 1 and 2 still point to the **original**
parameter object. The cross-dataset constraint silently fails to propagate.

**Fix**: Apply cross-dataset constraints to **every** angle sample of the target
dataset, not just the reference sample.

```python
runs = [226642, 226652]
all_experiments = []
all_samples = {}  # run -> list of samples (one per angle)

for run in runs:
    probes, samples, experiments = [], [], []

    for i, theta in enumerate(thetas):
        _probe = create_probe(partial_files[run][i], theta)
        _sample = create_sample()

        probes.append(_probe)
        samples.append(_sample)
        experiments.append(Experiment(sample=_sample, probe=_probe))

        if i == 0:
            _probe.intensity = Parameter(1, name="Intensity %d-1" % run)
            _probe.intensity.pm(0.15)
            _probe.sample_broadening.range(0, 0.05)
        else:
            _probe.intensity = Parameter(1, name="Intensity %d-%d" % (run, i + 1))
            _probe.intensity.pm(0.15)
            _probe.sample_broadening = probes[0].sample_broadening
            _probe.theta_offset = probes[0].theta_offset

            # Within same dataset: share ALL structural params across angles
            for layer in ["D2O", "CuOx", "Cu", "Ti"]:
                _sample[layer].material.rho = samples[0][layer].material.rho
                _sample[layer].thickness = samples[0][layer].thickness
                _sample[layer].interface = samples[0][layer].interface

    all_samples[run] = samples
    all_experiments.extend(experiments)

# Cross-dataset constraints — apply to ALL angle samples
ref1 = all_samples[226642][0]  # master parameters from dataset 1, angle 0

for s2 in all_samples[226652]:  # every angle sample in dataset 2
    s2["D2O"].material.rho = ref1["D2O"].material.rho
    s2["D2O"].interface = ref1["D2O"].interface
    s2["Cu"].material.rho = ref1["Cu"].material.rho
    s2["Cu"].thickness = ref1["Cu"].thickness
    s2["Cu"].interface = ref1["Cu"].interface
    s2["Ti"].material.rho = ref1["Ti"].material.rho
    s2["Ti"].thickness = ref1["Ti"].thickness
    s2["Ti"].interface = ref1["Ti"].interface
    # CuOx NOT tied — independent per dataset

problem = FitProblem(all_experiments)
```

## References

- [Domain knowledge](../refl1d-model/references/domain-knowledge.md)
- [Multi-angle probe creation](../refl1d-model/references/multi-angle-probe.md)
