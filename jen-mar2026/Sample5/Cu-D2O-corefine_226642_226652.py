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

    D2O = SLD("D2O", rho=6.196)
    Si = SLD("Si", rho=2.07)
    Ti = SLD("Ti", rho=-2.555)
    SiOx = SLD("SiOx", rho=2.83)
    Cu = SLD("Cu", rho=6.429)
    CuOx = SLD(name="CuOx", rho=4.88, irho=0.0)
    SEI = SLD(name="SEI", rho=5.6, irho=0.0)

    sample = D2O(0, 9.77) | CuOx(31.77, 12.3) | Cu(491, 4.848) | Ti(36.688, 9.11) | Si

    M = Experiment(sample=sample, probe=probe)

    sample["D2O"].material.rho.pm(0.02)
    sample["D2O"].interface.range(1, 15)

    sample["Ti"].thickness.range(30.0, 60.0)
    sample["Ti"].material.rho.range(-3.0, -1)
    sample["Ti"].interface.range(1.0, 22.0)

    sample["CuOx"].thickness.range(0.0, 200.0)
    sample["CuOx"].material.rho.range(4.5, 12)
    sample["CuOx"].interface.range(1.0, 33.0)

    sample["Cu"].thickness.range(400, 550)
    sample["Cu"].material.rho.range(2.0, 12)
    sample["Cu"].interface.pm(0.2)

    return M


# Auto-reduction directory
ar_dir = "/Users/jenni/OneDrive/SLAC/ORNL/March_2025/Sample5"

# Expt 10
data_file1 = os.path.join(ar_dir, "REFL_226642_combined_data_auto.txt")
data_file2 = os.path.join(ar_dir, "REFL_226652_combined_data_auto.txt")
#data_file3 = os.path.join(ar_dir, "REFL_218397_combined_data_auto.txt")

# First measurement
_refl = np.loadtxt(data_file1).T
experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])

# Second measurement
_refl = np.loadtxt(data_file2).T
experiment2 = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])

# Third measurement
#_refl = np.loadtxt(data_file3).T
#experiment3 = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])

# Constraints
experiment2.sample["D2O"].material.rho = experiment.sample["D2O"].material.rho
experiment2.sample["Cu"].material.rho = experiment.sample["Cu"].material.rho
experiment2.sample["Cu"].interface = experiment.sample["Cu"].interface
experiment2.sample["Ti"].material.rho = experiment.sample["Ti"].material.rho
experiment2.sample["Ti"].thickness = experiment.sample["Ti"].thickness
experiment2.sample["Ti"].interface = experiment.sample["Ti"].interface

#experiment3.sample["Cu"].material.rho = experiment.sample["Cu"].material.rho
#experiment3.sample["Cu"].interface = experiment.sample["Cu"].interface
#experiment3.sample["Ti"].material.rho = experiment.sample["Ti"].material.rho
#experiment3.sample["Ti"].thickness = experiment.sample["Ti"].thickness
#experiment3.sample["Ti"].interface = experiment.sample["Ti"].interface

problem = FitProblem([experiment, experiment2])
