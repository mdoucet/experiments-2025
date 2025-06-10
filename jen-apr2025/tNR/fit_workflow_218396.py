import os

from tron.bayesian_analysis import fitting_loop


# Data analysis directory
project_dir = os.path.expanduser("~/git/experiments-2025/jen-apr2025/")

# Upper-level data directory for the time-resolved data
storage_dir = "/SNS/REF_L/IPTS-34347/shared/tNR/"
storage_dir = os.path.expanduser("~/ORNL")

data_dir = os.path.join(storage_dir, "expt11-tNR-218396-120s")

# Directory where we store dynamic fit results
dyn_model_dir = os.path.join(storage_dir, "dyn-fitting")
if not os.path.exists(dyn_model_dir):
    os.makedirs(dyn_model_dir)

# Initial data set and model (starting point)
initial_err_file = os.path.join(
    project_dir, "results", "expt11-all-corefine", "Cu-THF-corefine-expt11-err.json"
)
initial_expt_file = os.path.join(
    project_dir, "results", "expt11-all-corefine", "Cu-THF-corefine-expt11-2-expt.json"
)

final_err_file = os.path.join(
    project_dir, "results", "expt11-all-corefine", "Cu-THF-corefine-expt11-err.json"
)
final_expt_file = os.path.join(
    project_dir, "results", "expt11-all-corefine", "Cu-THF-corefine-expt11-3-expt.json"
)

dynamic_run = 218396

# Create top-level directory for the dynamic fits
if not os.path.exists(os.path.join(dyn_model_dir, f"{dynamic_run}-dyn")):
    os.makedirs(os.path.join(dyn_model_dir, f"{dynamic_run}-dyn"))

store_basename = os.path.join(dyn_model_dir, f"{dynamic_run}-dyn/results-120s")

results_dir = os.path.join(dyn_model_dir, store_basename)

loop = fitting_loop.FittingLoop(
    data_dir,
    results_dir=results_dir,
    model_dir=project_dir,
    model_name="model-loop-218389",
    initial_err_file=initial_err_file,
    initial_expt_file=initial_expt_file,
    final_err_file=final_err_file,
    final_expt_file=final_expt_file,
)

print(loop)
loop.print_initial_final()

_file_list = sorted(os.listdir(data_dir))
_good_files = [_f for _f in _file_list if _f.startswith("r%d_t" % dynamic_run)]

loop.fit(_good_files, fit_forward=True)
