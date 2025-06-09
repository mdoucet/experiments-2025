import sys
import os

#sys.path.append(os.path.expanduser("~/git/LiquidsReflectometer/reduction"))


from lr_reduction import time_resolved

# NOTE: output directory must exist

reduced = time_resolved.reduce_slices(
    218389,
    "/SNS/REF_L/IPTS-34347/shared/autoreduce/template_down.xml",
    time_interval=120,
    output_dir="/SNS/REF_L/IPTS-34347/shared/tNR/expt11-tNR-218389-120s",
    scan_index=4,
    theta_offset=0,
    create_plot=True,
)
