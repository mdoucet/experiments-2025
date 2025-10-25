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
            #_error *= scale
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


def create_sample(has_SEI=True, sld_thf=5.5):
    THF = SLD("THF", rho=sld_thf)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-3)
    Cu = SLD("Cu", rho=6.3)
    material = SLD(name="material", rho=5, irho=0.0)
    SEI = SLD(name="SEI", rho=5.0, irho=0.0)

    if has_SEI:
        sample = (
            THF(0, 25) | SEI(150, 25) | material(50, 25) | Cu(500, 5) | Ti(30, 10) | Si
        )
    else:
        sample = THF(0, 25) | material(60, 20) | Cu(500, 5) | Ti(40, 10) | Si

    # sample['THF'].material.rho.range(4, 7)
    sample["THF"].material.rho.pm(1)
    if has_SEI:
        sample["THF"].interface.range(3, 65)
    else:
        sample["THF"].interface.range(3, 25)
    sample["Ti"].thickness.range(25.0, 60.0)
    sample["Ti"].material.rho.range(-4.0, 0)
    sample["Ti"].interface.range(1.0, 22.0)

    sample["Cu"].material.rho.range(5, 7)
    sample["Cu"].interface.range(1, 25)
    sample["Cu"].thickness.range(400, 600)

    sample["material"].thickness.range(10.0, 125.0)
    sample["material"].material.rho.range(1.0, 6.5)
    sample["material"].interface.range(3.0, 33.0)

    if has_SEI:
        sample["SEI"].thickness.range(5.0, 350.0)
        sample["SEI"].material.rho.range(0.0, 7.0)
        sample["SEI"].interface.range(3.0, 25.0)

    return sample


# Auto-reduction directory
ar_dir = "/SNS/REF_L/IPTS-34347/shared/autoreduce/"
ar_dir = "/home/mat/git/analyzer/data/partial/"
# OCV 1
sld_thf = 6.2
run = 218386
#run = 218397
#run = 218346 # Expt 9
#run = 218295 # Expt 1
run = 223772

data_files = [
    os.path.join(ar_dir, "REFL_%s_1_%s_partial.txt" % (run, run)),
    os.path.join(ar_dir, "REFL_%s_2_%s_partial.txt" % (run, run + 1)),
    os.path.join(ar_dir, "REFL_%s_3_%s_partial.txt" % (run, run + 2)),
]

thetas = [0.45, 1.2, 3.5]

probe_OCV1 = create_probe_from_list(data_files, thetas) #, scale=1.30)
sample_OCV1 = create_sample(False, sld_thf=sld_thf)

# OCV 2
run = 218393 # Expt 11

#sld_thf = 1
#run = 218403 # Exp 11
#run = 218354 # Expt 9
#run = 218298 # Expt 1

data_files = [
    os.path.join(ar_dir, "REFL_%s_1_%s_partial.txt" % (run, run)),
    os.path.join(ar_dir, "REFL_%s_2_%s_partial.txt" % (run, run + 1)),
    os.path.join(ar_dir, "REFL_%s_3_%s_partial.txt" % (run, run + 2)),
]

probe_OCV2 = create_probe_from_list(data_files, thetas)
sample_OCV2 = create_sample(False, sld_thf=sld_thf)

probe_OCV1.intensity = Parameter(1, name='Intensity 1')
probe_OCV1.intensity.dev(0.1)

probe_OCV2.intensity = Parameter(1, name="Intensity 2")
probe_OCV2.intensity.dev(0.1)

#probe_OCV1.intensity = probe_OCV2.intensity

# probe_OCV2.sample_broadening.range(0, 0.03)
# probe_OCV1.sample_broadening = probe_OCV2.sample_broadening

sample_OCV1["Cu"].material.rho = sample_OCV2["Cu"].material.rho
sample_OCV1["THF"].material.rho = sample_OCV2["THF"].material.rho
sample_OCV1["Ti"].material.rho = sample_OCV2["Ti"].material.rho
sample_OCV1["Ti"].thickness = sample_OCV2["Ti"].thickness
sample_OCV1["Ti"].interface = sample_OCV2["Ti"].interface
sample_OCV1["Cu"].interface = sample_OCV2["Cu"].interface


exp_OCV1 = Experiment(sample=sample_OCV1, probe=probe_OCV1)
exp_OCV2 = Experiment(sample=sample_OCV2, probe=probe_OCV2)

#problem = FitProblem([exp_OCV1, exp_OCV2])
problem = FitProblem([exp_OCV1])
