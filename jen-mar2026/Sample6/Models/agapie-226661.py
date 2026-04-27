"""Auto-generated analyzer model (agapie-226661) — created by create-model."""

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
    Platinum = SLD(name='Platinum', rho=6.288)
    Titanium = SLD(name='Titanium', rho=-1.95)
    Silicon = SLD(name='Silicon', rho=2.07)

    if back_reflection:
        sample = Air(0, 5.0) | Ionomer(300.0, 5.0) | Platinum(150.0, 5.0) | Titanium(30.0, 5.0) | Silicon
        sample['Air'].interface.range(5.0, 30.0)
    else:
        sample = Silicon | Titanium(30.0, 5.0) | Platinum(150.0, 5.0) | Ionomer(300.0, 5.0) | Air(0, 5.0)
        sample['Silicon'].interface.range(5.0, 30.0)

    # Parameter ranges
    sample['Air'].material.rho.range(-2.0, 2.0)
    sample['Ionomer'].thickness.range(5.0, 500.0)
    sample['Ionomer'].material.rho.range(-1.0, 3.0)
    sample['Ionomer'].interface.range(5.0, 30.0)
    sample['Platinum'].thickness.range(5.0, 300.0)
    sample['Platinum'].material.rho.range(4.288, 8.288)
    sample['Platinum'].interface.range(5.0, 15.0)
    sample['Titanium'].thickness.range(5.0, 100.0)
    sample['Titanium'].material.rho.range(-4.95, 1.05)
    sample['Titanium'].interface.range(5.0, 30.0)

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


# ── State: run_226661 (partials) ─────────────────────────
# All probes in this state share sample_run_226661, so every structural
# parameter (thickness, SLD, roughness) is tied across this state's
# segments by Python object identity — no explicit ties needed.
sample_run_226661 = create_sample(back_reflection=True)
theta_offset_run_226661 = Parameter(value=0.0, name="theta_offset_run_226661")
theta_offset_run_226661.range(-0.02, 0.02)
sample_broadening_run_226661 = Parameter(value=0.0, name="sample_broadening_run_226661")
sample_broadening_run_226661.range(0.0, 0.05)
probe_run_226661_1 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample6/Rawdata/REFL_226661_1_226661_partial.txt'), theta=0.3701745)
probe_run_226661_1.theta_offset = theta_offset_run_226661
probe_run_226661_1.sample_broadening = sample_broadening_run_226661
experiment_run_226661_1 = Experiment(probe=probe_run_226661_1, sample=sample_run_226661)
probe_run_226661_2 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample6/Rawdata/REFL_226661_2_226662_partial.txt'), theta=1.20005)
probe_run_226661_2.theta_offset = theta_offset_run_226661
probe_run_226661_2.sample_broadening = sample_broadening_run_226661
experiment_run_226661_2 = Experiment(probe=probe_run_226661_2, sample=sample_run_226661)
probe_run_226661_3 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample6/Rawdata/REFL_226661_3_226663_partial.txt'), theta=3.500245)
probe_run_226661_3.theta_offset = theta_offset_run_226661
probe_run_226661_3.sample_broadening = sample_broadening_run_226661
experiment_run_226661_3 = Experiment(probe=probe_run_226661_3, sample=sample_run_226661)

problem = FitProblem([experiment_run_226661_1, experiment_run_226661_2, experiment_run_226661_3])
