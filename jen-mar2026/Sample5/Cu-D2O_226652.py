import os
import numpy as np
from bumps.fitters import fit
from refl1d.names import *


def create_fit_experiment(q, dq, data, errors):
    # Go from FWHM to sigma
    dq /= 2.355

    # The QProbe object represents the beam
    probe = QProbe(q, dq, data=(data, errors))
    probe.intensity = Parameter(value=1, name='intensity')
    probe.intensity.pm(0.05)

    D2O = SLD('D2O', rho=6.19)
    Si = SLD('Si', rho=2.07)
    Ti = SLD('Ti', rho=-1.35)
    SiOx = SLD('SiOx', rho=2.83)
    Cu = SLD('Cu', rho=6.428)
    CuOx = SLD(name='CuOx', rho=4.889, irho=0.0)
    SEI = SLD(name='SEI', rho=5.6, irho=0.0)
    
    sample = ( D2O(0, 9.27) | CuOx(37.36, 12.246) | Cu(486.8, 4.848) | Ti(36.17, 1.634) | Si )

    M = Experiment(sample=sample, probe=probe)

    sample['D2O'].material.rho.range(6.18, 6.25)
    sample['D2O'].interface.range(1, 50)
    
    sample['Ti'].thickness.range(36.0, 36.2)
    sample['Ti'].material.rho.range(-4.0, -1)
    sample['Ti'].interface.range(1, 2)

    sample['CuOx'].thickness.range(30, 40)
    sample['CuOx'].material.rho.range(4, 6)
    sample['CuOx'].interface.range(9, 50)

    sample['Cu'].thickness.range(450, 550)
    sample['Cu'].material.rho.range(6.43, 6.44)
    sample['Cu'].interface.range(1.0, 12.0)
    
    return M

# Auto-reduction directory
ar_dir = '/Users/jenni/OneDrive/SLAC/ORNL/March_2025/Sample5/'

data_file = os.path.join(ar_dir, 'REFL_226652_combined_data_auto.txt')

_refl = np.loadtxt(data_file).T

experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])
problem = FitProblem(experiment)
