import os
import numpy as np
from bumps.fitters import fit
from refl1d.names import *
from refl1d.probe import make_probe


def create_probe(data_file, theta):
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0 * q

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

def create_fit_experiment(q, dq, data, errors, probe=None):
    # Go from FWHM to sigma
    dq /= 2.355

    # The QProbe object represents the beam
    if not probe:
        probe = QProbe(q, dq, R=data, dR=errors)

    probe.intensity = Parameter(value=1, name="intensity")
    probe.intensity.pm(0.1)

    THF = SLD("THF", rho=6.0)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-2.69)
    SiOx = SLD("SiOx", rho=2.83)
    Cu = SLD("Cu", rho=6.344)
    material = SLD(name="material", rho=4.92, irho=0.0)
    SEI = SLD(name="SEI", rho=5.6, irho=0.0)

    sample = THF(0, 24.92) | material(53, 13.48) | Cu(500, 6.24) | Ti(35.38, 1.7) | Si

    M = Experiment(sample=sample, probe=probe)

    sample["THF"].material.rho.range(4.5, 6.4)
    sample["THF"].interface.range(1, 25)

    sample["Ti"].thickness.range(30.0, 60.0)
    sample["Ti"].material.rho.range(-3.0, -1)
    sample["Ti"].interface.range(1.0, 22.0)

    sample["material"].thickness.range(10.0, 200.0)
    sample["material"].material.rho.range(3.0, 12)
    sample["material"].interface.range(1.0, 33.0)

    sample["Cu"].thickness.range(1.0, 1000.0)
    sample["Cu"].material.rho.range(2.0, 12)
    sample["Cu"].interface.range(1.0, 12.0)

    return M


# Auto-reduction directory
ar_dir = "/SNS/REF_L/IPTS-34347/shared/autoreduce/"
tnr_dir = "/SNS/REF_L/IPTS-34347/shared/isaac/tNR/per-file/"
# Expt 10
data_file1 = os.path.join(ar_dir, "REFL_218386_combined_data_auto.txt")  # Expt 11 OCV
data_file2 = os.path.join(tnr_dir, "r218389_hold_initial_1.txt")  # Expt 11 OCV
#data_file1 = os.path.join(ar_dir, "REFL_218393_combined_data_auto.txt")  # Expt 11 OCV
#data_file2 = os.path.join(tnr_dir, "r218389_hold_gap_14_5.txt")  # Expt 11 OCV
#data_file2 = "/SNS/REF_L/IPTS-34347/shared/tNR/expt11-tNR-218389-120s/r218389_t006480.txt"
#data_file2 = "/SNS/REF_L/IPTS-34347/shared/tNR/expt11-tNR-218389-120s/r218389_t000000.txt"
# First measurement

_refl = np.loadtxt(data_file1).T

probe1 = QProbe(_refl[0], _refl[3], R=_refl[1], dR=_refl[2])
experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])#, probe=None)

# Second measurement
_refl = np.loadtxt(data_file2).T
probe2 = create_probe(data_file2, 0.6)
experiment2 = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2], probe=probe2)

# Constraints
if True:
    experiment2.sample["THF"].material.rho = experiment.sample["THF"].material.rho
    experiment2.sample["material"].material.rho = experiment.sample["material"].material.rho
    experiment2.sample["material"].interface = experiment.sample["material"].interface
    experiment2.sample["material"].thickness = experiment.sample["material"].thickness
    experiment2.sample["Cu"].material.rho = experiment.sample["Cu"].material.rho
    experiment2.sample["Cu"].interface = experiment.sample["Cu"].interface
    experiment2.sample["Cu"].thickness = experiment.sample["Cu"].thickness
    experiment2.sample["Ti"].material.rho = experiment.sample["Ti"].material.rho
    experiment2.sample["Ti"].thickness = experiment.sample["Ti"].thickness
    experiment2.sample["Ti"].interface = experiment.sample["Ti"].interface


problem = FitProblem([experiment, experiment2])
