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


def create_sample(sld_D2O=6.3, has_oxide=True):
    D2O = SLD("D2O", rho=sld_D2O)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-4)
    Pt = SLD("Pt", rho=6.288)
    ionomer = SLD(name="ionomer", rho=1.5, irho=0.0)
    PtOx = SLD(name="PtOx", rho=5.0, irho=0.0)
    SiOx = SLD("SiOx", rho=2.8)
    #SL1 = SLD(name="SL1", rho=1)

    if has_oxide:
        sample = (
            D2O(0, 25) | ionomer(500, 20) | PtOx(50, 5) | Pt(145.318, 7.4) | Ti(26.797,11.4) | Si
        )
    else:
        sample = D2O(0, 25) | ionomer(500, 20) | Pt(145.318, 7.4) | Ti(26.797,11.4) | SiOx(18.1, 4.6) | Si

    #| SL1(50, 20) 
    # sample['D2O'].material.rho.range(4, 7)
    sample["D2O"].material.rho.range(5.5, 7.0)

    sample["Ti"].thickness.pm(0.5)
    sample["Ti"].material.rho.range(-4.0, 0)
    sample["Ti"].interface.pm(1.0)

    sample["Pt"].material.rho.pm(0.07)
    sample["Pt"].interface.range(5,10)
    sample["Pt"].thickness.pm(0.2)

    sample["ionomer"].thickness.range(10, 2000.0)
    sample["ionomer"].material.rho.range(-1.0, 6.5)
    sample["ionomer"].interface.range(3.0, 33.0)

    if has_oxide:
        sample["PtOx"].thickness.range(5.0, 50.0)
        sample["PtOx"].material.rho.range(5.7, 6.2)
        sample["PtOx"].interface.range(3.0, 25.0)
    else:
        sample["SiOx"].thickness.pm(0.5)
        sample["SiOx"].material.rho.range(2.0, 4.2)
        sample["SiOx"].interface.range(1.0, 10)

        #sample["SL1"].thickness.range(2.0, 2000.0)
        #sample["SL1"].material.rho.range(-1.0, 6.5)
        #sample["SL1"].interface.range(3.0, 33.0)

    return sample


# Auto-reduction directory
ar_dir = "/Users/jenni/OneDrive/Desktop/GitHub_2/experiments-2025/jen-mar2026/Sample6/Rawdata"

# OCV 1
sld_D2O = 0
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
    os.path.join(ar_dir, "REFL_226664_1_226664_partial.txt" ),
    os.path.join(ar_dir, "REFL_226664_2_226665_partial.txt"),
    os.path.join(ar_dir, "REFL_226664_3_226666_partial.txt"),
]

thetas = [0.45, 1.2, 3.5]

has_oxide = False

probes = []
samples = []
experiments = []


for i in range(len(data_files)):
    _probe = create_probe_from_list([data_files[i]], [thetas[i]])  # , scale=1.30)
    _sample = create_sample(sld_D2O=sld_D2O, has_oxide=has_oxide)

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

        #_sample["SL1"].material.rho = samples[0]["SL1"].material.rho
        #_sample["SL1"].thickness = samples[0]["SL1"].thickness
        #_sample["SL1"].interface = samples[0]["SL1"].interface

        _sample["D2O"].material.rho = samples[0]["D2O"].material.rho


problem = FitProblem(experiments)
