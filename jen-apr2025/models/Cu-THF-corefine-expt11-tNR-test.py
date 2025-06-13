import os
import json
import numpy as np
import copy

from bumps.fitters import fit
from bumps import serialize

import refl1d
from refl1d.names import QProbe, Parameter, SLD, Slab, Experiment, FitProblem
from refl1d.probe import make_probe

model_expt_json_file = "/SNS/users/m2d/git/experiments-2025/jen-apr2025/results/expt11-all-corefine/Cu-THF-corefine-expt11-1-expt.json"


def fix_all_parameters(expt, verbose=False):
    """
    Fix all the parameters within an Experiment object

    Parameters
    ----------
    expt : Experiment
        Experiment object to process
    verbose : bool
        If True, print out parameters that were not fixed
    """
    pars = expt.parameters()

    def _fix_parameters(item):
        if type(item) is Parameter:
            if verbose and not item.fixed:
                print("Found %s" % item)
            item.fixed = True
        elif type(item) is list:
            for p in item:
                _fix_parameters(p)
        elif type(item) is dict:
            for p, v in item.items():
                _fix_parameters(v)
        else:
            print("Found unknown parameter: %s" % item)

    _fix_parameters(pars)


def get_probe_parts(data_file, theta):
    q, data, errors, dq = np.loadtxt(data_file).T
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0 * q
    return data, errors, wl, dL, dT


def create_fit_experiment(q, dq, data, errors, theta=0.6):
    # Go from FWHM to sigma
    dq /= 2.355

    with open(model_expt_json_file, "rt") as input_file:
        serialized = input_file.read()
        serialized_dict = json.loads(serialized)
        expt = serialize.deserialize(serialized_dict, migration=True)

    # The QProbe object represents the beam
    wl = 4 * np.pi * np.sin(np.pi / 180 * theta) / q
    dT = dq / q * np.tan(np.pi / 180 * theta) * 180 / np.pi
    dL = 0 * q
    probe = make_probe(
        T=theta,
        dT=dT,
        L=wl,
        dL=dL,
        data=(data, errors),
        radiation="neutron",
        resolution="uniform",
    )

    expt = Experiment(probe=probe, sample=copy.copy(expt.sample))

    fix_all_parameters(expt, verbose=False)
    return expt


def limited_constraints(experiment_list):
    # The first experiment in the list is the reference
    experiment_list[0].probe.theta_offset.range(-0.05, 0.05)
    experiment_list[0].probe.intensity = Parameter(value=0.87, name="intensity")
    experiment_list[0].probe.intensity.pm(0.15)
    experiment_list[0].sample["Cu"].thickness.pm(100)

    for i, experiment2 in enumerate(experiment_list[1:], start=1):
        experiment2.probe.theta_offset = experiment_list[0].probe.theta_offset
        experiment2.probe.intensity = experiment_list[0].probe.intensity
        # Cu thickness
        experiment2.sample["Cu"].thickness = Parameter(
            name=f"Cu thickness {i}", value=experiment2.sample["Cu"].thickness.value
        )
        experiment2.sample["Cu"].thickness.pm(100)

        # THF interface
        experiment2.sample["THF"].interface = Parameter(
            name=f"THF interface {i}", value=experiment2.sample["THF"].interface.value
        )
        experiment2.sample["THF"].interface.pm(15)

        # material rho
        experiment2.sample["material"].material.rho = Parameter(
            name=f"material rho {i}",
            value=experiment2.sample["material"].material.rho.value,
        )
        experiment2.sample["material"].material.rho.pm(1.5)

        # material thickness
        experiment2.sample["material"].thickness = Parameter(
            name=f"material thickness {i}",
            value=experiment2.sample["material"].thickness.value,
        )
        experiment2.sample["material"].thickness.pm(50)

        # material interface
        experiment2.sample["material"].interface = Parameter(
            name=f"material interface {i}",
            value=experiment2.sample["material"].interface.value,
        )
        experiment2.sample["material"].interface.pm(15)


def linear_constraints(experiment_list):
    # The first experiment in the list is the reference
    # experiment_list[0].probe.theta_offset.range(-0.005, 0.005)
    offset = Parameter(name="theta offset", value=-0.013)
    experiment_list[0].probe.theta_offset = offset
    experiment_list[0].probe.intensity = Parameter(value=0.87, name="intensity")
    experiment_list[0].probe.intensity.pm(0.15)
    experiment_list[0].sample["Cu"].thickness.pm(100)

    for i, experiment2 in enumerate(experiment_list[1:], start=1):
        experiment2.probe.theta_offset = experiment_list[0].probe.theta_offset
        experiment2.probe.intensity = experiment_list[0].probe.intensity

        experiment2.sample["Cu"].thickness = Parameter(
            name=f"Cu thickness {i}",
            value=experiment2.sample["Cu"].thickness.value,
        )
        experiment2.sample["Cu"].thickness.pm(100)

    interfaces = ["THF", "material"]

    n_curves = len(experiment_list)
    for item in interfaces:
        if True:
            p_start = Parameter(
                name=f"{item} interface start",
                value=experiment_list[0].sample[item].interface.value,
            )
            p_end = Parameter(
                name=f"{item} interface end",
                value=experiment_list[-1].sample[item].interface.value,
            )
            # p_start.range(1, 25)
            p_end.range(1, 25)
            for i, experiment2 in enumerate(experiment_list):
                # Set the interface as a linear function of the index
                interface_value = p_start + (p_end - p_start) * i / (n_curves - 1)
                experiment2.sample[item].interface = interface_value
        else:
            for i, experiment2 in enumerate(experiment_list):
                experiment2.sample[item].interface = (
                    experiment2.sample["material"].thickness / 2.0
                )

    thicknesses = ["material"]  # , "Cu"]
    for item in thicknesses:
        p_start = Parameter(
            name=f"{item} thickness start",
            value=experiment_list[0].sample[item].thickness.value,
        )
        p_end = Parameter(
            name=f"{item} thickness end",
            value=experiment_list[-1].sample[item].thickness.value,
        )
        # p_start.pm(50)
        p_end.pm(50)

        # slope = (
        #    experiment_list[-1].sample[item].thickness.value
        #    - experiment_list[0].sample[item].thickness.value
        # ) / (n_curves - 1)
        # p_slope = Parameter(name=f"{item} thickness slope", value=slope)
        # p_slope.pm(50/ n_curves)

        # p_acc = Parameter(name=f"{item} thickness acc", value=0)
        # p_acc.pm(1)

        for i, experiment2 in enumerate(experiment_list):
            # Set the thickness as a linear function of the index
            thickness_value = p_start + (p_end - p_start) * i / (n_curves - 1)
            # thickness_value = p_start + p_slope * i + p_acc * i**2
            experiment2.sample[item].thickness = thickness_value

    rhos = ["material"]
    for item in rhos:
        p_start = Parameter(
            name=f"{item} rho start",
            value=experiment_list[0].sample[item].material.rho.value,
        )
        p_end = Parameter(
            name=f"{item} rho end",
            value=experiment_list[-1].sample[item].material.rho.value,
        )
        # p_start.pm(1.5)
        p_end.pm(1.5)
        for i, experiment2 in enumerate(experiment_list):
            # Set the rho as a linear function of the index
            rho_value = p_start + (p_end - p_start) * i / (n_curves - 1)
            experiment2.sample[item].material.rho = rho_value


# Auto-reduction directory

ar_dir = "/SNS/REF_L/IPTS-34347/shared/tNR/"
data_set = "expt11-tNR-218389-120s"

data_files = sorted(os.listdir(os.path.join(ar_dir, data_set)))
data_files = [
    os.path.join(ar_dir, data_set, f) for f in data_files if f.endswith(".txt")
]
print("Found %d data files in %s" % (len(data_files), os.path.join(ar_dir, data_set)))

data_files = data_files[::2]


expt_list = []

_refl = np.loadtxt(os.path.join(ar_dir, data_set, data_files[0])).T
experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])
expt_list.append(experiment)

_refl = np.loadtxt(os.path.join(ar_dir, "..", 
                                "autoreduce", "REFL_218386_combined_data_auto.txt")).T

experiment = create_fit_experiment(_refl[0], _refl[3], _refl[1], _refl[2])
expt_list.append(experiment)

problem = FitProblem(expt_list)
