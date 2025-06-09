import sys
import numpy as np
import os

sys.path.append(os.path.expanduser('~/git/TrON/'))
from refl1d.names import QProbe, Parameter, FitProblem
from tron.bayesian_analysis import model_utils


# Parse input arguments ########################################################
# First argument is the data file to use
reduced_file = sys.argv[1]

# Second argument is the starting model [experiment description]
expt_file = sys.argv[2]

# Third argument is the error information used for setting the prior
err_file = sys.argv[3]

# 0.1 was used so far (Jan 2023) with good results
prior_scale = 0.1

# Load data ####################################################################
q_min = 0.0
q_max = 0.4

try:
    Q, R, dR, dQ = np.loadtxt(reduced_file).T
except:
    Q, R, dR = np.loadtxt(reduced_file).T
    dQ = 0.028*Q

i_min = np.min([i for i in range(len(Q)) if Q[i]>q_min])
i_max = np.max([i for i in range(len(Q)) if Q[i]<q_max])+1

# SNS data is FWHM
dQ_std = dQ/2.35
probe = QProbe(Q[i_min:i_max], dQ_std[i_min:i_max], data=(R[i_min:i_max], dR[i_min:i_max]))

# Experiment ###################################################################
expt = model_utils.expt_from_json_file(expt_file, probe=probe,
                                       model_err_json_file=err_file,
                                       prior_scale=prior_scale, set_ranges=False)

sample = expt.sample
#sample['THF'].material.rho.range(4.5, 6.4)
sample['THF'].interface.range(5.0, 30.0)
sample['material'].material.rho.range(3.0, 12.0)
sample['material'].thickness.range(10.0, 200.0)
sample['material'].interface.range(5.0, 33.0)
#sample['Cu'].material.rho.range(2.0, 12.0)
sample['Cu'].thickness.range(1.0, 1000.0)
#sample['Cu'].interface.range(1.0, 12.0)
#sample['Ti'].material.rho.range(-3.0, -1.0)
#sample['Ti'].thickness.range(30.0, 60.0)
#sample['Ti'].interface.range(1.0, 22.0)


#probe.intensity.range(0.6, 1.1)
probe.intensity = Parameter(value=0.87, name="intensity")
probe.intensity.dev(0.05)

################################################################################
problem = FitProblem(expt)
