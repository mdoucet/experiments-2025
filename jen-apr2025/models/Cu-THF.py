import os
import numpy as np
from bumps.fitters import fit
from refl1d.names import *


def create_fit_experiment(q, dq, data, errors, has_SEI=True):
    # Go from FWHM to sigma
    dq /= 2.355

    # The QProbe object represents the beam
    probe = QProbe(q, dq, data=(data, errors))
    probe.intensity = Parameter(value=1, name='intensity')
    probe.intensity.pm(0.05)

    THF = SLD('THF', rho=6.0)
    Si = SLD('Si', rho=2.07)
    Ti = SLD('Ti', rho=-2.69)
    SiOx = SLD('SiOx', rho=2.83)
    Cu = SLD('Cu', rho=6.344)
    material = SLD(name='material', rho=4.92, irho=0.0)
    SEI = SLD(name='SEI', rho=5.6, irho=0.0)
    
    if has_SEI:
        sample = ( THF(0, 5) | SEI(30, 20) | material(42, 4) | Cu(500, 10) | Ti(44, 5) | Si )
    else:
        sample = ( THF(0, 24.92) | material(53, 13.48) | Cu(500, 6.24) | Ti(35.38, 1.7) | Si )

    M = Experiment(sample=sample, probe=probe)

    sample['THF'].material.rho.range(4.5, 6.4)
    sample['THF'].interface.range(1, 25)
    
    sample['Ti'].thickness.range(30.0, 60.0)
    sample['Ti'].material.rho.range(-3.0, -1)
    sample['Ti'].interface.range(1.0, 22.0)

    sample['material'].thickness.range(10.0, 200.0)
    sample['material'].material.rho.range(5.0, 12)
    sample['material'].interface.range(1.0, 33.0)

    sample['Cu'].thickness.range(1.0, 1000.0)
    sample['Cu'].material.rho.range(2.0, 12)
    sample['Cu'].interface.range(1.0, 12.0)
    
    if has_SEI:
        sample['SEI'].thickness.range(5.0, 150.0)
        sample['SEI'].material.rho.range(0.0, 12.0)
        sample['SEI'].interface.range(1.0, 25.0)

    #sample['SiOx'].thickness.range(5.0, 45.0)
    #sample['SiOx'].material.rho.range(0.0, 4.0)
    #sample['SiOx'].interface.range(1.0, 10.0)

    return M

# Auto-reduction directory
ar_dir = '/SNS/REF_L/IPTS-33612/shared/autoreduce/'
#ar_dir = '/Users/m2d/git/experiments-2024/val-sep24/data/'
ar_dir = '/SNS/REF_L/IPTS-34347/shared/autoreduce/'

data_file = os.path.join(ar_dir, 'REFL_218278_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_218281_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_218284_combined_data_auto.txt')
_refl = np.loadtxt(data_file).T

experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2], has_SEI=False)
problem = FitProblem(experiment)
