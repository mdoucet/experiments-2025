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


def create_sample(sld_h2o=5.5, has_oxide=True):
    H2O = SLD("H2O", rho=sld_h2o)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-4)
    Cu = SLD("Cu", rho=6.269)
    ionomer = SLD(name="ionomer", rho=2.7, irho=0.0)
    CuOx = SLD(name="CuOx", rho=5.0, irho=0.0)
    SiOx = SLD("SiOx", rho=2.90)
    SL1 = SLD("SL1", rho=0.129, irho=0)
    SL2 = SLD(name="SL2", rho=5, irho=0.0)

    if has_oxide:
        sample = H2O(0, 25) | SL2(20,3) | ionomer(950, 20) | SL1(30,3) | CuOx(50, 5) | Cu(151.9, 12) | Ti(26.7, 4.38) | SiOx(20.151, 1) | Si
    else:
        sample = H2O(0, 25) | SL2(157,49) | ionomer(950, 20) | SL1(19,9) | Cu(151.9, 12) | Ti(26.7, 4.38) | SiOx(20.151, 1) | Si

    # sample['H2O'].material.rho.range(4, 7)
    sample["H2O"].material.rho.range(6.3, 6.37)

    sample["Ti"].thickness.pm(0.15)
    sample["Ti"].material.rho.range(-4.0, 0)
    sample["Ti"].interface.range(1.0, 22.0)

    sample["Cu"].material.rho.pm(0.02)
    sample["Cu"].interface.range(1, 25)
    sample["Cu"].thickness.pm(0.5)

    sample["ionomer"].thickness.range(200.0, 2000.0)
    sample["ionomer"].material.rho.range(2.3, 3)
    sample["ionomer"].interface.range(3.0, 33.0)



    if has_oxide:
        sample["CuOx"].thickness.range(5.0, 50.0)
        sample["CuOx"].material.rho.range(5.7, 6.2)
        sample["CuOx"].interface.range(3.0, 25.0)

        sample["SiOx"].thickness.pm(0.5)
        sample["SiOx"].material.rho.range(2.0, 4.2)
        sample["SiOx"].interface.range(1.0, 15.0)

        sample["SL2"].thickness.range(0, 200.0)
        sample["SL2"].material.rho.range(0, 7)
        sample["SL2"].interface.range(1.0, 20.0)

        sample["SL1"].thickness.range(0,200)
        sample["SL1"].material.rho.range(0,7)
        sample["SL1"].interface.range(1,10)
    else:
        sample["SiOx"].thickness.pm(0.5)
        sample["SiOx"].material.rho.range(2.0, 4.2)
        sample["SiOx"].interface.range(1.0, 15.0)

        sample["SL2"].thickness.range(0, 200.0)
        sample["SL2"].material.rho.range(0, 7)
        sample["SL2"].interface.range(1.0, 20.0)

        sample["SL1"].thickness.range(0,200)
        sample["SL1"].material.rho.range(0,7)
        sample["SL1"].interface.range(1,10)

    return sample


# Auto-reduction directory
ar_dir = "/Users/jenni/OneDrive/Desktop/GitHub_2/experiments-2025/jen-mar2026/Sample4/Rawdata"

# OCV 1
sld_h2o = 0
run = 223915 # air
run = 223875 # OCV
run = 223879 # -0.2V
run = 223883 # -0.4V
run = 223887 # OCV 2
run = 223890 # OCV 3 
run = 223901
run = 223912
#jen removed %run from data files
data_files = [
    os.path.join(ar_dir, "REFL_226628_1_226628_partial.txt" ),
    os.path.join(ar_dir, "REFL_226628_2_226629_partial.txt"),
    os.path.join(ar_dir, "REFL_226628_3_226630_partial.txt"),
]

thetas = [0.45, 1.2, 3.5]

has_oxide = True

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
        _probe.intensity.range(0.6, 1)
        _probe.theta_offset.dev(0.05)
    else:
        # _probe.intensity = probes[0].intensity
        _probe.intensity = Parameter(1, name="Intensity %d" % (i + 1))
        _probe.intensity.range(0.6, 1)
        _probe.sample_broadening = probes[0].sample_broadening
        _probe.theta_offset = probes[0].theta_offset

        _sample["Ti"].material.rho = samples[0]["Ti"].material.rho
        _sample["Ti"].thickness = samples[0]["Ti"].thickness
        _sample["Ti"].interface = samples[0]["Ti"].interface
        _sample["Cu"].material.rho = samples[0]["Cu"].material.rho
        _sample["Cu"].thickness = samples[0]["Cu"].thickness
        _sample["Cu"].interface = samples[0]["Cu"].interface


        if has_oxide:
            _sample["CuOx"].material.rho = samples[0]["CuOx"].material.rho
            _sample["CuOx"].thickness = samples[0]["CuOx"].thickness
            _sample["CuOx"].interface = samples[0]["CuOx"].interface

            _sample["SiOx"].material.rho = samples[0]["SiOx"].material.rho
            _sample["SiOx"].thickness = samples[0]["SiOx"].thickness
            _sample["SiOx"].interface = samples[0]["SiOx"].interface
            
            _sample["SL2"].material.rho = samples[0]["SL2"].material.rho
            _sample["SL2"].thickness = samples[0]["SL2"].thickness
            _sample["SL2"].interface = samples[0]["SL2"].interface

            _sample["SL1"].material.rho = samples[0]["SL1"].material.rho
            _sample["SL1"].thickness = samples[0]["SL1"].thickness
            _sample["SL1"].interface = samples[0]["SL1"].interface
        else:
            _sample["SiOx"].material.rho = samples[0]["SiOx"].material.rho
            _sample["SiOx"].thickness = samples[0]["SiOx"].thickness
            _sample["SiOx"].interface = samples[0]["SiOx"].interface
            
            _sample["SL2"].material.rho = samples[0]["SL2"].material.rho
            _sample["SL2"].thickness = samples[0]["SL2"].thickness
            _sample["SL2"].interface = samples[0]["SL2"].interface

            _sample["SL1"].material.rho = samples[0]["SL1"].material.rho
            _sample["SL1"].thickness = samples[0]["SL1"].thickness
            _sample["SL1"].interface = samples[0]["SL1"].interface

        _sample["ionomer"].material.rho = samples[0]["ionomer"].material.rho
        _sample["ionomer"].thickness = samples[0]["ionomer"].thickness
        _sample["ionomer"].interface = samples[0]["ionomer"].interface



        _sample["H2O"].material.rho = samples[0]["H2O"].material.rho


problem = FitProblem(experiments)
