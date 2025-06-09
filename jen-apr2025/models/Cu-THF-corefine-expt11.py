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
    probe.intensity.pm(0.05)

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

# Expt 10
data_file1 = os.path.join(ar_dir, "REFL_218386_combined_data_auto.txt")  # Expt 11 OCV
data_file2 = os.path.join(ar_dir, "REFL_218393_combined_data_auto.txt")  # Expt 11 -0.8 V
data_file3 = os.path.join(ar_dir, "REFL_218397_combined_data_auto.txt")  # Expt 11 OCV

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
experiment2.sample["Cu"].material.rho = experiment.sample["Cu"].material.rho
experiment2.sample["Cu"].interface = experiment.sample["Cu"].interface
experiment2.sample["Ti"].material.rho = experiment.sample["Ti"].material.rho
experiment2.sample["Ti"].thickness = experiment.sample["Ti"].thickness
experiment2.sample["Ti"].interface = experiment.sample["Ti"].interface


experiment3.sample["THF"].material.rho = experiment.sample["THF"].material.rho
experiment3.sample["Cu"].material.rho = experiment.sample["Cu"].material.rho
experiment3.sample["Cu"].interface = experiment.sample["Cu"].interface
experiment3.sample["Ti"].material.rho = experiment.sample["Ti"].material.rho
experiment3.sample["Ti"].thickness = experiment.sample["Ti"].thickness
experiment3.sample["Ti"].interface = experiment.sample["Ti"].interface


problem = FitProblem([experiment, experiment2, experiment3])
