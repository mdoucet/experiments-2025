## Multi-Angle Probe Creation

When fitting partial files (one file per angle), use `make_probe` instead of `QProbe`.
This approach preserves angle-dependent resolution information.

### Single-Angle Probe

```python
from refl1d.probe import make_probe

def create_probe(data_file, theta):
    """Create a probe from a single partial data file.
    
    Args:
        data_file: Path to REFL_{run}_{seg}_{subrun}_partial.txt
        theta: Incident angle in degrees (TwoTheta / 2)
    """
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0.0 * q  # Or use moderator model (see below)

    probe = make_probe(
        T=theta, dT=dT, L=wl, dL=dL,
        data=(data, errors),
        radiation="neutron",
        resolution="uniform",
    )
    return probe
```

### Concatenated Multi-Angle Probe

Concatenate all angles into a single probe for one experiment:

```python
def get_probe_parts(data_file, theta):
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0.0 * q
    return data, errors, wl, dL, dT


def create_probe_from_list(data_files, thetas):
    theta_list, data_list, err_list = [], [], []
    wl_list, dL_list, dT_list = [], [], []

    for i, (f, th) in enumerate(zip(data_files, thetas)):
        _data, _error, _wl, _dL, _dT = get_probe_parts(f, th)
        data_list.append(_data)
        err_list.append(_error)
        wl_list.append(_wl)
        dL_list.append(_dL)
        dT_list.append(_dT)
        theta_list.append(np.ones(len(_wl)) * th)

    probe = make_probe(
        T=np.concatenate(theta_list),
        dT=np.concatenate(dT_list),
        L=np.concatenate(wl_list),
        dL=np.concatenate(dL_list),
        data=(np.concatenate(data_list), np.concatenate(err_list)),
        radiation="neutron",
        resolution="uniform",
    )
    return probe
```

### Wavelength Resolution Model (Don't use unless explicitely stated)

For higher fidelity, include the moderator emission time contribution:

```python
def delta_wl_over_wl(wl):
    """SNS moderator emission time → wavelength resolution."""
    dtof = 0.0148 * wl**3 - 0.5233 * wl**2 + 6.4797 * wl + 231.99
    dtof[wl > 2] = (
        392.31 * wl**6 - 3169.3 * wl**5 + 10445 * wl**4
        - 17872 * wl**3 + 16509 * wl**2 - 7448.4 * wl + 1280.5
    )
    dwl = 3.9560 * dtof / (1000000 * 15.282 * wl)
    return dwl

# Then in get_probe_parts:
dL = delta_wl_over_wl(wl) * q  # instead of 0.0 * q
```

### Probe Parameters for Multi-Angle Fits

When fitting multiple angles of the same sample, share structural parameters
but allow independent intensity per angle:

```python
probes[0].intensity = Parameter(1, name="Intensity 1")
probes[0].intensity.pm(0.15)
probes[0].sample_broadening.range(0, 0.05)

for i in range(1, len(probes)):
    probes[i].intensity = Parameter(1, name="Intensity %d" % (i + 1))
    probes[i].intensity.pm(0.15)
    probes[i].sample_broadening = probes[0].sample_broadening  # shared
    probes[i].theta_offset = probes[0].theta_offset             # shared
```
