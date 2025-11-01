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


def create_sample(sld_thf=5.5, has_oxide=False):
    THF = SLD("THF", rho=sld_thf)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-2.8)
    Cu = SLD("Cu", rho=6.4)
    CuTop = SLD("CuTop", rho=6.4)
    material = SLD(name="material", rho=4.9, irho=0.0)

    sample = THF(0, 15) | material(25, 15) | CuTop(30, 5) | Cu(486, 9) | Ti(30, 2) | Si

    # sample['THF'].material.rho.range(4, 7)
    sample["THF"].material.rho.pm(1)
    sample["THF"].interface.range(3, 25)
    sample["Ti"].thickness.range(25.0, 60.0)
    sample["Ti"].material.rho.range(-4.0, 0)
    sample["Ti"].interface.range(1.0, 22.0)

    sample["Cu"].material.rho.range(5, 7)
    sample["Cu"].interface.range(1, 25)
    sample["Cu"].thickness.range(400, 600)
    sample["CuTop"].material.rho.range(5, 7)
    sample["CuTop"].interface.range(1, 25)
    sample["CuTop"].thickness.range(400, 600)

    sample["material"].thickness.range(10.0, 125.0)
    sample["material"].material.rho.range(1.0, 6.5)
    sample["material"].interface.range(3.0, 33.0)

    return sample


# Auto-reduction directory
ar_dir = "/SNS/REF_L/IPTS-34347/shared/autoreduce/"
ar_dir = "/home/mat/git/analyzer/data/partial/"
# OCV 1
sld_thf = 6.2
run = 223872 # air
run = 223875 # OCV
run = 223879 # -0.2V
run = 223883 # -0.4V
run = 223887 # OCV 2
run = 223890 # OCV 3 

run = 223901
run = 224011

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
    _sample = create_sample(sld_thf=sld_thf, has_oxide=has_oxide)

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
        _sample["Cu"].material.rho = samples[0]["Cu"].material.rho
        _sample["Cu"].thickness = samples[0]["Cu"].thickness
        _sample["Cu"].interface = samples[0]["Cu"].interface
        _sample["CuTop"].material.rho = samples[0]["CuTop"].material.rho
        _sample["CuTop"].thickness = samples[0]["CuTop"].thickness
        _sample["CuTop"].interface = samples[0]["CuTop"].interface

        if has_oxide:
            _sample["SiOx"].material.rho = samples[0]["SiOx"].material.rho
            _sample["SiOx"].thickness = samples[0]["SiOx"].thickness
            _sample["SiOx"].interface = samples[0]["SiOx"].interface

        _sample["material"].material.rho = samples[0]["material"].material.rho
        _sample["material"].thickness = samples[0]["material"].thickness
        _sample["material"].interface = samples[0]["material"].interface

        _sample["THF"].material.rho = samples[0]["THF"].material.rho

        _sample["THF"].interface = samples[0]["THF"].interface


problem = FitProblem(experiments)
