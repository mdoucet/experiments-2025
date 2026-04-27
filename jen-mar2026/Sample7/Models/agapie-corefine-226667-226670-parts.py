"""Auto-generated analyzer model (agapie-corefine-226667-226670-parts) — created by create-model."""

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
    ambient = SLD(name='ambient', rho=0.0)
    ionomer = SLD(name='ionomer', rho=1.0)
    platinum = SLD(name='platinum', rho=6.288)
    titanium = SLD(name='titanium', rho=-1.95)
    silicon = SLD(name='silicon', rho=2.07)

    if back_reflection:
        sample = ambient(0, 5.0) | ionomer(300.0, 5.0) | platinum(150.0, 5.0) | titanium(30.0, 5.0) | silicon
        sample['ambient'].interface.range(5.0, 30.0)
    else:
        sample = silicon | titanium(30.0, 5.0) | platinum(150.0, 5.0) | ionomer(300.0, 5.0) | ambient(0, 5.0)
        sample['silicon'].interface.range(5.0, 15.0)

    # Parameter ranges
    sample['ambient'].material.rho.range(0.0, 6.5)
    sample['ionomer'].thickness.range(20.0, 1000.0)
    sample['ionomer'].material.rho.range(-1.0, 6.0)
    sample['ionomer'].interface.range(5.0, 30.0)
    sample['platinum'].thickness.range(10.0, 20.0)
    sample['platinum'].material.rho.range(4.288, 8.288)
    sample['platinum'].interface.range(5.0, 30.0)
    sample['titanium'].thickness.range(5.0, 50.0)
    sample['titanium'].material.rho.range(-4.95, 1.05)
    sample['titanium'].interface.range(5.0, 15.0)

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


# ── State: run_226667 (partials) ─────────────────────────
# All probes in this state share sample_run_226667, so every structural
# parameter (thickness, SLD, roughness) is tied across this state's
# segments by Python object identity — no explicit ties needed.
sample_run_226667 = create_sample(back_reflection=True)
theta_offset_run_226667 = Parameter(value=0.0, name="theta_offset_run_226667")
theta_offset_run_226667.range(-0.02, 0.02)
sample_broadening_run_226667 = Parameter(value=0.0, name="sample_broadening_run_226667")
sample_broadening_run_226667.range(0.0, 0.05)
probe_run_226667_1 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample7/Rawdata/REFL_226667_1_226667_partial.txt'), theta=0.370097)
probe_run_226667_1.theta_offset = theta_offset_run_226667
probe_run_226667_1.sample_broadening = sample_broadening_run_226667
experiment_run_226667_1 = Experiment(probe=probe_run_226667_1, sample=sample_run_226667)
probe_run_226667_2 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample7/Rawdata/REFL_226667_2_226668_partial.txt'), theta=1.200375)
probe_run_226667_2.theta_offset = theta_offset_run_226667
probe_run_226667_2.sample_broadening = sample_broadening_run_226667
experiment_run_226667_2 = Experiment(probe=probe_run_226667_2, sample=sample_run_226667)
probe_run_226667_3 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample7/Rawdata/REFL_226667_3_226669_partial.txt'), theta=3.500305)
probe_run_226667_3.theta_offset = theta_offset_run_226667
probe_run_226667_3.sample_broadening = sample_broadening_run_226667
experiment_run_226667_3 = Experiment(probe=probe_run_226667_3, sample=sample_run_226667)

# ── State: run_226670 (partials) ─────────────────────────
# All probes in this state share sample_run_226670, so every structural
# parameter (thickness, SLD, roughness) is tied across this state's
# segments by Python object identity — no explicit ties needed.
sample_run_226670 = create_sample(back_reflection=True)
theta_offset_run_226670 = Parameter(value=0.0, name="theta_offset_run_226670")
theta_offset_run_226670.range(-0.02, 0.02)
sample_broadening_run_226670 = Parameter(value=0.0, name="sample_broadening_run_226670")
sample_broadening_run_226670.range(0.0, 0.05)
probe_run_226670_1 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample7/Rawdata/REFL_226670_1_226670_partial.txt'), theta=0.3700725)
probe_run_226670_1.theta_offset = theta_offset_run_226670
probe_run_226670_1.sample_broadening = sample_broadening_run_226670
experiment_run_226670_1 = Experiment(probe=probe_run_226670_1, sample=sample_run_226670)
probe_run_226670_2 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample7/Rawdata/REFL_226670_2_226671_partial.txt'), theta=1.200395)
probe_run_226670_2.theta_offset = theta_offset_run_226670
probe_run_226670_2.sample_broadening = sample_broadening_run_226670
experiment_run_226670_2 = Experiment(probe=probe_run_226670_2, sample=sample_run_226670)
probe_run_226670_3 = create_angle_probe(os.path.join(os.path.expanduser('~'), 'git/experiments-2025/jen-mar2026/Sample7/Rawdata/REFL_226670_3_226672_partial.txt'), theta=3.500155)
probe_run_226670_3.theta_offset = theta_offset_run_226670
probe_run_226670_3.sample_broadening = sample_broadening_run_226670
experiment_run_226670_3 = Experiment(probe=probe_run_226670_3, sample=sample_run_226670)

# ── Shared structural parameters across states ──
sample_run_226670['ionomer'].thickness = sample_run_226667['ionomer'].thickness
sample_run_226670['ionomer'].material.rho = sample_run_226667['ionomer'].material.rho
sample_run_226670['ionomer'].interface = sample_run_226667['ionomer'].interface
sample_run_226670['platinum'].thickness = sample_run_226667['platinum'].thickness
sample_run_226670['platinum'].material.rho = sample_run_226667['platinum'].material.rho
sample_run_226670['platinum'].interface = sample_run_226667['platinum'].interface
sample_run_226670['titanium'].thickness = sample_run_226667['titanium'].thickness
sample_run_226670['titanium'].material.rho = sample_run_226667['titanium'].material.rho
sample_run_226670['titanium'].interface = sample_run_226667['titanium'].interface
sample_run_226670['silicon'].interface = sample_run_226667['silicon'].interface

problem = FitProblem([experiment_run_226667_1, experiment_run_226667_2, experiment_run_226667_3, experiment_run_226670_1, experiment_run_226670_2, experiment_run_226670_3])
