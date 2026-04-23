"""Auto-generated analyzer model (ionomer-air-226622) — created by create-model."""

import os
import numpy as np
from bumps.fitters import fit
from refl1d.names import *

from refl1d.probe import make_probe

def create_sample(back_reflection=False):
    """Build a fresh sample stack for one state.

    Every probe in a given state is constructed with this stack,
    so all structural parameters (thickness / SLD / interface) are
    automatically tied across that state's segments via Python
    object identity. Each state gets its OWN call, so structural
    parameters are independent across states unless explicitly
    tied below with ``sample_B['X'].attr = sample_A['X'].attr``.

    ``back_reflection`` selects stack orientation so the default
    ``probe.back_reflectivity=False`` gives correct physics in both
    buried-interface and standard front-reflection geometries.
    """
    Air = SLD(name='Air', rho=0.0)
    Ionomer = SLD(name='Ionomer', rho=1.0)
    Copper = SLD(name='Copper', rho=6.55)
    Titanium = SLD(name='Titanium', rho=-1.95)
    Silicon = SLD(name='Silicon', rho=2.07)

    if back_reflection:
        sample = Air(0, 5.0) | Ionomer(500.0, 5.0) | Copper(150.0, 5.0) | Titanium(30.0, 5.0) | Silicon
        sample['Air'].interface.range(5.0, 15.0)
    else:
        sample = Silicon | Titanium(30.0, 5.0) | Copper(150.0, 5.0) | Ionomer(500.0, 5.0) | Air(0, 5.0)
        sample['Silicon'].interface.range(5.0, 15.0)

    # Parameter ranges
    sample['Air'].material.rho.range(-2.0, 2.0)
    sample['Ionomer'].thickness.range(450.0, 550.0)
    sample['Ionomer'].material.rho.range(-1.0, 3.0)
    sample['Ionomer'].interface.range(5.0, 15.0)
    sample['Copper'].thickness.range(135.0, 165.0)
    sample['Copper'].material.rho.range(4.55, 8.55)
    sample['Copper'].interface.range(5.0, 15.0)
    sample['Titanium'].thickness.range(27.0, 33.0)
    sample['Titanium'].material.rho.range(-4.95, 1.05)
    sample['Titanium'].interface.range(5.0, 15.0)

    return sample


def create_q_probe(data_file):
    """Angle-independent probe for a combined REF_L file."""
    q, data, errors, dq = np.loadtxt(data_file).T
    dq = dq / 2.355  # FWHM → 1-sigma
    probe = QProbe(q, dq, data=(data, errors))
    probe.intensity = Parameter(value=1.0, name="intensity")
    probe.intensity.range(0.0, 10.0)
    return probe

def create_angle_probe(data_file, theta):
    """Angle-based probe from one REF_L partial file."""
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0 * q  # wavelength resolution placeholder
    probe = make_probe(
        T=theta, dT=dT, L=wl, dL=dL,
        data=(data, errors),
        radiation="neutron",
        resolution="uniform",
    )
    probe.intensity = Parameter(value=1.0, name="intensity")
    probe.intensity.range(0.0, 10.0)
    return probe


# ── State: run_226622 (partials) ─────────────────────────
# All probes in this state share sample_run_226622, so every structural
# parameter (thickness, SLD, roughness) is tied across this state's
# segments by Python object identity — no explicit ties needed.
sample_run_226622 = create_sample(back_reflection=True)
theta_offset_run_226622 = Parameter(value=0.0, name="theta_offset_run_226622")
theta_offset_run_226622.range(-0.02, 0.02)
sample_broadening_run_226622 = Parameter(value=0.0, name="sample_broadening_run_226622")
sample_broadening_run_226622.range(0.0, 0.05)
probe_run_226622_1 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample4/Rawdata/REFL_226622_1_226622_partial.txt'), theta=0.370019)
probe_run_226622_1.theta_offset = theta_offset_run_226622
probe_run_226622_1.sample_broadening = sample_broadening_run_226622
experiment_run_226622_1 = Experiment(probe=probe_run_226622_1, sample=sample_run_226622)
probe_run_226622_2 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample4/Rawdata/REFL_226622_2_226623_partial.txt'), theta=1.20019)
probe_run_226622_2.theta_offset = theta_offset_run_226622
probe_run_226622_2.sample_broadening = sample_broadening_run_226622
experiment_run_226622_2 = Experiment(probe=probe_run_226622_2, sample=sample_run_226622)
probe_run_226622_3 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample4/Rawdata/REFL_226622_3_226624_partial.txt'), theta=3.500355)
probe_run_226622_3.theta_offset = theta_offset_run_226622
probe_run_226622_3.sample_broadening = sample_broadening_run_226622
experiment_run_226622_3 = Experiment(probe=probe_run_226622_3, sample=sample_run_226622)

problem = FitProblem([experiment_run_226622_1, experiment_run_226622_2, experiment_run_226622_3])
