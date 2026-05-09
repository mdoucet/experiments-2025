import os
import numpy as np
from refl1d.names import *
from refl1d.probe import make_probe


def get_probe_parts(data_file, theta):
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0.0 * q
    return data, errors, wl, dL, dT


def create_probe_from_list(data_files, thetas):
    theta_list, data_list, err_list = [], [], []
    wl_list, dL_list, dT_list = [], [], []

    for f, th in zip(data_files, thetas):
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


def create_sample():
    D2O = SLD("D2O", rho=6.19)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-1.95)
    Cu = SLD("Cu", rho=6.55)
    CuOx = SLD(name="CuOx", rho=5.0, irho=0.0)

    # Stack: surface → substrate (beam enters from substrate side)
    sample = D2O(0, 10) | CuOx(20, 5) | Cu(400, 5) | Ti(35, 5) | Si

    sample["D2O"].material.rho.pm(0.5)
    sample["D2O"].interface.range(1, 25)

    sample["CuOx"].thickness.range(5, 200)
    sample["CuOx"].material.rho.range(3.0, 7.0)
    sample["CuOx"].interface.range(3, 25)

    sample["Cu"].thickness.range(300, 600)
    sample["Cu"].material.rho.range(5.0, 7.0)
    sample["Cu"].interface.range(1.0, 12.0)

    sample["Ti"].thickness.range(20.0, 50.0)
    sample["Ti"].material.rho.range(-4.0, 0.0)
    sample["Ti"].interface.range(1.0, 15.0)

    return sample


# --- Data paths ---
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Rawdata")
first_run = os.path.expanduser("~/data/226642_new_reduction.dat")

thetas = [0.37, 1.20, 3.50]

runs = [226642, 226652]
partial_files = {}
for run in runs:
    partial_files[run] = [
        os.path.join(data_dir, "REFL_%d_1_%d_partial.txt" % (run, run)),
        os.path.join(data_dir, "REFL_%d_2_%d_partial.txt" % (run, run + 1)),
        os.path.join(data_dir, "REFL_%d_3_%d_partial.txt" % (run, run + 2)),
    ]

# --- Build experiments: 3 angles × 2 datasets = 6 experiments ---
all_experiments = []

# Store all samples/probes per dataset for cross-dataset constraints
all_samples = {}
all_probes = {}

for run in runs:
    probes = []
    samples = []
    experiments = []

    for i, theta in enumerate(thetas):
        _probe = create_probe_from_list([partial_files[run][i]], [theta])
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

            # Within same dataset: share ALL structural parameters across angles
            for layer in ["D2O", "CuOx", "Cu", "Ti"]:
                _sample[layer].material.rho = samples[0][layer].material.rho
                _sample[layer].thickness = samples[0][layer].thickness
                _sample[layer].interface = samples[0][layer].interface

    all_samples[run] = samples
    all_probes[run] = probes
    all_experiments.extend(experiments)

# --- Cross-dataset constraints ---
# Tie everything EXCEPT CuOx between the two datasets.
# Must apply to ALL angle samples (not just angle 0) because each angle's
# sample holds its own reference to the parameter object.
ref1 = all_samples[226642][0]

for s2 in all_samples[226652]:
    # Tie D2O
    s2["D2O"].material.rho = ref1["D2O"].material.rho
    s2["D2O"].interface = ref1["D2O"].interface

    # Tie Cu
    s2["Cu"].material.rho = ref1["Cu"].material.rho
    #s2["Cu"].thickness = ref1["Cu"].thickness
    #s2["Cu"].interface = ref1["Cu"].interface

    # Tie Ti
    #s2["Ti"].material.rho = ref1["Ti"].material.rho
    s2["Ti"].thickness = ref1["Ti"].thickness
    s2["Ti"].interface = ref1["Ti"].interface

# CuOx is NOT tied — independent for each dataset

problem = FitProblem(all_experiments)
