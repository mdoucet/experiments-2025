from inspect import Parameter
import os
import numpy as np
from bumps.fitters import fit
from refl1d.names import *
from refl1d.probe import make_probe


def delta_wl_over_wl(wl):
    """
    The moderator emission time equation gives a time in microsecond.
    """
    dtof = 0.0148 * wl * wl * wl - 0.5233 * wl * wl + 6.4797 * wl + 231.99
    dtof[wl > 2] = (
        392.31 * wl**6
        - 3169.3 * wl**5
        + 10445 * wl**4
        - 17872 * wl * wl * wl
        + 16509 * wl * wl
        - 7448.4 * wl
        + 1280.5
    )
    dwl = 3.9560 * dtof / (1000000 * 15.282 * wl)
    return dwl


def get_probe_parts(data_file, theta):
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = delta_wl_over_wl(wl) * q
    return data, errors, wl, dL, dT


def create_probe_from_list(data_files, thetas, scale=1):
    theta_list = []
    data_list = []
    err_list = []
    wl_list = []
    dL_list = []
    dT_list = []

    for i in range(len(data_files)):
        _data, _error, _wl, _dL, _dT = get_probe_parts(data_files[i], thetas[i])
        if i == 2:
            _data *= scale
            # _error *= scale
            _error += abs(_data * (1 - scale) / 2.0)
        data_list.append(_data)
        err_list.append(_error)
        wl_list.append(_wl)
        dL_list.append(_dL)
        dT_list.append(_dT)
        theta_list.append(np.ones(len(_wl)) * thetas[i])

    theta = np.concatenate(theta_list)
    data = np.concatenate(data_list)
    errors = np.concatenate(err_list)
    wl = np.concatenate(wl_list)
    dL = np.concatenate(dL_list)
    dT = np.concatenate(dT_list)

    probe = make_probe(
        T=theta,
        dT=dT,
        L=wl,
        dL=dL,
        data=(data, errors),
        radiation="neutron",
        resolution="uniform",
    )
    return probe


def create_probe(data_file, theta):
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = delta_wl_over_wl(wl) * q

    # The following is how refl1d computes dQ
    # dQ = (4 * np.pi / wl) * np.sqrt((np.sin(np.pi/180*theta) * dL / wl) ** 2 + (np.cos(np.pi/180*theta) * dT * np.pi/180) ** 2)

    # dT and dL are FWHM
    probe = make_probe(
        T=theta,
        dT=dT,
        L=wl,
        dL=dL,
        data=(data, errors),
        radiation="neutron",
        resolution="uniform",
    )
    return probe


def create_sample(sld_h2o=5.5, has_oxide=False):
    H2O = SLD("H2O", rho=sld_h2o)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-4)
    Pt = SLD("Pt", rho=6.4)
    ionomer = SLD(name="ionomer", rho=1.5, irho=0.0)
    PtOx = SLD(name="PtOx", rho=5.0, irho=0.0)
    SiOx = SLD("SiOx", rho=2.90)

    if has_oxide:
        sample = (
            H2O(0, 25) | ionomer(1200, 20) | PtOx(50, 5) | Pt(150, 5) | Ti(40, 10) | Si
        )
    else:
        sample = H2O(0, 25) | ionomer(600, 20) | Pt(151, 12) | Ti(26.7, 4.38) | SiOx(18.5, 1) | Si

    # sample['H2O'].material.rho.range(4, 7)
    sample["H2O"].material.rho.pm(1)

    sample["Ti"].thickness.range(10.0, 60.0)
    sample["Ti"].material.rho.range(-4.0, 0)
    sample["Ti"].interface.range(1.0, 22.0)

    sample["Pt"].material.rho.range(5, 7)
    sample["Pt"].interface.range(1, 25)
    sample["Pt"].thickness.range(100, 250)

    sample["ionomer"].thickness.range(200.0, 2000.0)
    sample["ionomer"].material.rho.range(-1.0, 6.5)
    sample["ionomer"].interface.range(3.0, 33.0)

    if has_oxide:
        sample["PtOx"].thickness.range(5.0, 50.0)
        sample["PtOx"].material.rho.range(-1.0, 7.0)
        sample["PtOx"].interface.range(3.0, 25.0)
    else:
        sample["SiOx"].thickness.range(5.0, 40.0)
        sample["SiOx"].material.rho.range(2.0, 4.2)
        sample["SiOx"].interface.range(1.0, 15.0)

    return sample


# Auto-reduction directory
ar_dir = "/SNS/REF_L/IPTS-34347/shared/autoreduce/"
ar_dir = "/home/mat/git/analyzer/data/partial/"
# OCV 1
sld_h2o = 0
run = 223872 # air
run = 223875 # OCV
run = 223879 # -0.2V
run = 223883 # -0.4V
run = 223887 # OCV 2
run = 223890 # OCV 3 

run = 223901
run = 223912

data_files = [
    os.path.join(ar_dir, "REFL_%s_1_%s_partial.txt" % (run, run)),
    os.path.join(ar_dir, "REFL_%s_2_%s_partial.txt" % (run, run + 1)),
    os.path.join(ar_dir, "REFL_%s_3_%s_partial.txt" % (run, run + 2)),
]

thetas = [0.45, 1.2, 3.5]

has_oxide = False

probes = []
samples = []
experiments = []


for i in range(len(data_files)):
    _probe = create_probe_from_list([data_files[i]], [thetas[i]])  # , scale=1.30)
    _sample = create_sample(sld_h2o=sld_h2o, has_oxide=has_oxide)

    probes.append(_probe)
    samples.append(_sample)
    experiments.append(Experiment(sample=_sample, probe=_probe))

    if i == 0:
        _probe.intensity = Parameter(1, name="Intensity %d" % (i + 1))
        _probe.theta_offset.dev(0.05)
    else:
        _probe.intensity = probes[0].intensity
        _probe.sample_broadening = probes[0].sample_broadening
        _probe.theta_offset = probes[0].theta_offset

        _sample["Ti"].material.rho = samples[0]["Ti"].material.rho
        _sample["Ti"].thickness = samples[0]["Ti"].thickness
        _sample["Ti"].interface = samples[0]["Ti"].interface
        _sample["Pt"].material.rho = samples[0]["Pt"].material.rho
        _sample["Pt"].thickness = samples[0]["Pt"].thickness
        _sample["Pt"].interface = samples[0]["Pt"].interface

        if has_oxide:
            _sample["PtOx"].material.rho = samples[0]["PtOx"].material.rho
            _sample["PtOx"].thickness = samples[0]["PtOx"].thickness
            _sample["PtOx"].interface = samples[0]["PtOx"].interface
        else:
            _sample["SiOx"].material.rho = samples[0]["SiOx"].material.rho
            _sample["SiOx"].thickness = samples[0]["SiOx"].thickness
            _sample["SiOx"].interface = samples[0]["SiOx"].interface

        _sample["ionomer"].material.rho = samples[0]["ionomer"].material.rho
        _sample["ionomer"].thickness = samples[0]["ionomer"].thickness
        _sample["ionomer"].interface = samples[0]["ionomer"].interface

        _sample["H2O"].material.rho = samples[0]["H2O"].material.rho


problem = FitProblem(experiments)
