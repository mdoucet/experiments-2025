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

    THF = SLD('THF', rho=6.0)
    Si = SLD('Si', rho=2.07)
    Ti = SLD('Ti', rho=-2.69)
    SiOx = SLD('SiOx', rho=2.83)
    Pt = SLD('Pt', rho=6.344)
    material = SLD(name='material', rho=4.92, irho=0.0)
    SEI = SLD(name='SEI', rho=5.6, irho=0.0)
    
    sample = ( THF(0, 24.92) | material(53, 13.48) | Pt(500, 6.24) | Ti(35.38, 1.7) | Si )

    M = Experiment(sample=sample, probe=probe)

    sample['THF'].material.rho.range(4.5, 6.4)
    sample['THF'].interface.range(1, 25)
    
    sample['Ti'].thickness.range(30.0, 60.0)
    sample['Ti'].material.rho.range(-3.0, -1)
    sample['Ti'].interface.range(1.0, 22.0)

    sample['material'].thickness.range(10.0, 200.0)
    sample['material'].material.rho.range(5.0, 12)
    sample['material'].interface.range(1.0, 33.0)

    sample['Pt'].thickness.range(1.0, 1000.0)
    sample['Pt'].material.rho.range(2.0, 12)
    sample['Pt'].interface.range(1.0, 12.0)
    
    return M

# Auto-reduction directory
ar_dir = '/SNS/REF_L/IPTS-34347/shared/autoreduce/'

data_file = os.path.join(ar_dir, 'REFL_218301_combined_data_auto.txt')

_refl = np.loadtxt(data_file).T

experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])
problem = FitProblem(experiment)
