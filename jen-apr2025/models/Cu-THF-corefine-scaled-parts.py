import os
import numpy as np
from bumps.fitters import fit
from refl1d.names import *


def create_fit_experiment(q, dq, data, errors):
    # Go from FWHM to sigma
    dq /= 2.355

    # The QProbe object represents the beam
    probe = QProbe(q, dq, data=(data, errors))
    probe.intensity = Parameter(value=1, name="intensity")
    probe.intensity.range(0.3, 1.5)

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
    sample["material"].material.rho.range(5.0, 12)
    sample["material"].interface.range(1.0, 33.0)

    sample["Cu"].thickness.range(1.0, 1000.0)
    sample["Cu"].material.rho.range(2.0, 12)
    sample["Cu"].interface.range(1.0, 12.0)

    return M


# Auto-reduction directory
ar_dir = "/SNS/REF_L/IPTS-34347/shared/autoreduce/"
ar_dir = "/home/mat/git/analyzer/data/partial/"

run = 223772 # First OCV
run = 223787

data_files = [
    os.path.join(ar_dir, "REFL_%s_1_%s_partial.txt" % (run, run)),
    os.path.join(ar_dir, "REFL_%s_2_%s_partial.txt" % (run, run + 1)),
    os.path.join(ar_dir, "REFL_%s_3_%s_partial.txt" % (run, run + 2)),
]


# Expt 10
data_file1 = data_files[0]
data_file2 = data_files[1]
data_file3 = data_files[2]

# First measurement
_refl = np.loadtxt(data_file1).T
experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])

# Second measurement
_refl = np.loadtxt(data_file2).T
experiment2 = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])

# Second measurement
_refl = np.loadtxt(data_file3).T
experiment3 = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])

# Constraints
experiment2.sample["THF"].material.rho = experiment.sample["THF"].material.rho
experiment2.sample["THF"].interface = experiment.sample["THF"].interface
experiment2.sample["Cu"].material.rho = experiment.sample["Cu"].material.rho
experiment2.sample["Cu"].interface = experiment.sample["Cu"].interface
experiment2.sample["Cu"].thickness = experiment.sample["Cu"].thickness
experiment2.sample["Ti"].material.rho = experiment.sample["Ti"].material.rho
experiment2.sample["Ti"].thickness = experiment.sample["Ti"].thickness
experiment2.sample["Ti"].interface = experiment.sample["Ti"].interface
experiment2.sample["material"].material.rho = experiment.sample["material"].material.rho
experiment2.sample["material"].thickness = experiment.sample["material"].thickness
experiment2.sample["material"].interface = experiment.sample["material"].interface


experiment3.sample["THF"].material.rho = experiment.sample["THF"].material.rho
experiment3.sample["THF"].interface = experiment.sample["THF"].interface
experiment3.sample["Cu"].material.rho = experiment.sample["Cu"].material.rho
experiment3.sample["Cu"].interface = experiment.sample["Cu"].interface
experiment3.sample["Cu"].thickness = experiment.sample["Cu"].thickness
experiment3.sample["Ti"].material.rho = experiment.sample["Ti"].material.rho
experiment3.sample["Ti"].thickness = experiment.sample["Ti"].thickness
experiment3.sample["Ti"].interface = experiment.sample["Ti"].interface
experiment3.sample["material"].material.rho = experiment.sample["material"].material.rho
experiment3.sample["material"].thickness = experiment.sample["material"].thickness
experiment3.sample["material"].interface = experiment.sample["material"].interface


problem = FitProblem([experiment, experiment2, experiment3])
