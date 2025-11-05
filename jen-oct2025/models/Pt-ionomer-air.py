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

    air = SLD('air', rho=0)
    H2O = SLD('H2O', rho=-0.2)
    Si = SLD('Si', rho=2.07)
    Pt = SLD('Pt', rho=6.2)
    Ti = SLD('Ti', rho=-2.9)
    ionomer = SLD('ionomer', rho=0.6)
    material = SLD(name='material', rho=5.46, irho=0.0)
    SEI = SLD(name='SEI', rho=5.6, irho=0.0)
    SiOx = SLD(name='SiOx', rho=3.2, irho=0.0)
    PtOx = SLD(name='PtOx', rho=5, irho=0)

    #sample = ( H2O(0, 50) | ionomer(1070, 4.6) | Pt(151, 6.1) | Ti(26, 9) | SiOx(20, 1) | Si )

    sample = ( Si(0, 1) | SiOx(20, 1) | Ti(26,9) | Pt(151, 6) | ionomer(600, 5) | air )
    M = Experiment(sample=sample, probe=probe)


    sample['Si'].interface.range(1, 25)

    sample['SiOx'].thickness.range(10.0, 60.0)
    sample['SiOx'].material.rho.range(2, 4.2)
    sample['SiOx'].interface.range(1.0, 22.0)


    #sample['PtOx'].thickness.range(10.0, 60.0)
    #sample['PtOx'].material.rho.range(-1, 8)
    #sample['PtOx'].interface.range(1.0, 22.0)

    sample['Pt'].thickness.range(100.0, 260.0)
    sample['Pt'].material.rho.range(5, 8)
    sample['Pt'].interface.range(1.0, 22.0)

    sample['ionomer'].thickness.range(800.0, 2000.0)
    sample['ionomer'].material.rho.range(2.0, 7)
    sample['ionomer'].interface.range(1.0, 12.0)

    sample['Ti'].thickness.range(10.0, 60.0)
    sample['Ti'].material.rho.range(-3.0, -1)
    sample['Ti'].interface.range(1.0, 22.0)

    return M

# Auto-reduction directory
ar_dir = '/SNS/REF_L/IPTS-35338/shared/autoreduce/'

data_file = os.path.join(ar_dir, 'REFL_223852_combined_data_auto.txt')
#data_file = os.path.join(ar_dir, 'REFL_223818_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223855_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223863_combined_data_auto.txt')

data_file = os.path.join(ar_dir, 'REFL_223866_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223869_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223872_combined_data_auto.txt')

data_file = os.path.join(ar_dir, 'REFL_223915_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223948_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223951_combined_data_auto.txt')
#data_file = os.path.join(ar_dir, 'REFL_223915_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223954_combined_data_auto.txt')
data_file = os.path.join(ar_dir, 'REFL_223957_combined_data_auto.txt')

#ar_dir = '/SNS/REF_L/IPTS-35338/shared/'
#data_file = os.path.join(ar_dir, 'REFL_223872_reduced_data.txt')
#data_file = os.path.join(ar_dir, 'REFL_223897_reduced_data.txt')



_refl = np.loadtxt(data_file).T

experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])
problem = FitProblem(experiment)
