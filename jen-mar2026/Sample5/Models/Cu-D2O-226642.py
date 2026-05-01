import os
import numpy as np
from bumps.fitters import fit
from refl1d.names import *
from refl1d.probe import make_probe


def create_sample(back_reflection=False):
    D2O = SLD(name='D2O', rho=6.19)
    Si = SLD(name='Si', rho=2.07)
    Ti = SLD(name='Ti', rho=-1.35)
    Cu = SLD(name='Cu', rho=6.599)
    CuOx = SLD(name='CuOx', rho=4.889, irho=0.0)

    if back_reflection:
        sample = D2O(0, 9.27) | CuOx(37.36, 12.246) | Cu(380, 4.848) | Ti(36.17, 1.634) | Si
        sample['D2O'].interface.range(1.0, 50.0)
    else:
        sample = Si | Ti(36.17, 1.634) | Cu(380, 4.848) | CuOx(37.36, 12.246) | D2O(0, 9.27)
        sample['Si'].interface.range(1.0, 50.0)

    sample['D2O'].material.rho.range(6.15, 6.20)
    sample['Ti'].thickness.range(36.0, 36.2)
    sample['Ti'].material.rho.range(-4.0, -1.0)
    sample['Ti'].interface.range(1.0, 2.0)
    sample['CuOx'].thickness.range(30.0, 40.0)
    sample['CuOx'].material.rho.range(4.0, 6.0)
    sample['CuOx'].interface.range(9.0, 50.0)
    sample['Cu'].thickness.range(380, 550.0)
    sample['Cu'].material.rho.range(6.595, 6.60)
    sample['Cu'].interface.range(1.0, 12.0)

    return sample


def create_q_probe(data_file):
    q, data, errors, dq = np.loadtxt(data_file).T
    dq = dq / 2.355
    probe = QProbe(q, dq, data=(data, errors))
    probe.intensity = Parameter(value=1.0, name='intensity')
    probe.intensity.pm(0.05)
    return probe


def create_angle_probe(data_file, theta):
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0 * q
    probe = make_probe(
        T=theta,
        dT=dT,
        L=wl,
        dL=dL,
        data=(data, errors),
        radiation='neutron',
        resolution='uniform',
    )
    probe.intensity = Parameter(value=1.0, name='intensity')
    probe.intensity.pm(0.05)
    return probe


# ── State: run_226642 (partials) ─────────────────────────
sample_run_226642 = create_sample(back_reflection=True)
theta_offset_run_226642 = Parameter(value=0.0, name='theta_offset_run_226642')
theta_offset_run_226642.range(-0.02, 0.02)
sample_broadening_run_226642 = Parameter(value=0.0, name='sample_broadening_run_226642')
sample_broadening_run_226642.range(0.0, 0.05)

rawdata_dir = os.path.join(
    os.path.expanduser('~'),
    'OneDrive',
    'Desktop',
    'GitHub_2',
    'experiments-2025',
    'jen-mar2026',
    'Sample5',
    'Rawdata',
)

probe_run_226642_1 = create_angle_probe(
    os.path.join(rawdata_dir, 'REFL_226642_1_226642_partial.txt'),
    theta=0.3699945,
)
probe_run_226642_1.theta_offset = theta_offset_run_226642
probe_run_226642_1.sample_broadening = sample_broadening_run_226642
experiment_run_226642_1 = Experiment(probe=probe_run_226642_1, sample=sample_run_226642)

probe_run_226642_2 = create_angle_probe(
    os.path.join(rawdata_dir, 'REFL_226642_2_226643_partial.txt'),
    theta=1.20016,
)
probe_run_226642_2.theta_offset = theta_offset_run_226642
probe_run_226642_2.sample_broadening = sample_broadening_run_226642
experiment_run_226642_2 = Experiment(probe=probe_run_226642_2, sample=sample_run_226642)

probe_run_226642_3 = create_angle_probe(
    os.path.join(rawdata_dir, 'REFL_226642_3_226644_partial.txt'),
    theta=3.50005,
)
probe_run_226642_3.theta_offset = theta_offset_run_226642
probe_run_226642_3.sample_broadening = sample_broadening_run_226642
experiment_run_226642_3 = Experiment(probe=probe_run_226642_3, sample=sample_run_226642)

problem = FitProblem([
    experiment_run_226642_1,
    experiment_run_226642_2,
    experiment_run_226642_3,
])
