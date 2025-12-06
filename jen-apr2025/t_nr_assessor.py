"""
Time-resolved neutron reflectometry (tNR) assessor module.

This module provides tools for analyzing time-series neutron reflectometry data.

https://numerics.mathdotnet.com/Distance
"""

import os
import re
import glob
from typing import Tuple, List

import numpy as np
import matplotlib.pyplot as plt


def load_tnr_data(filepath: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load time-resolved neutron reflectometry data from an ASCII file.

    Parameters
    ----------
    filepath : str
        Path to the data file.

    Returns
    -------
    tuple
        (Q, R, dR) arrays where Q is the momentum transfer, R is the reflectivity,
        and dR is the uncertainty on R.
    """
    data = np.loadtxt(filepath, usecols=(0, 1, 2))
    return data[:, 0], data[:, 1], data[:, 2]


def parse_timestamp_from_filename(filename: str) -> int:
    """
    Extract the timestamp from a tNR filename.

    Parameters
    ----------
    filename : str
        Filename in the format r<run_id>_t<timestamp>.txt

    Returns
    -------
    int
        The timestamp value extracted from the filename.
    """
    match = re.search(r"_t(\d+)\.txt$", filename)
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not parse timestamp from filename: {filename}")


def compute_distance_time_series(
    data_dir: str, pattern: str = "r*_t*.txt", output_plot: str = None, t_nr_run: str = ""
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the distance time series relative to the first measurement.

    This function loads all tNR data files matching the pattern, sorts them by
    timestamp, and calculates the distance between each curve and the
    reference curve (t=0).

    Parameters
    ----------
    data_dir : str
        Directory containing the tNR data files.
    pattern : str, optional
        Glob pattern to match data files. Default is "r*_t*.txt".
    output_plot : str, optional
        If provided, save the plot to this file path.

    Returns
    -------
    tuple
        (timestamps, distance_values) arrays where timestamps are the time values
        and distance_values are the corresponding distance values from the reference.
    """
    # Find all matching files
    file_pattern = os.path.join(data_dir, pattern)
    files = glob.glob(file_pattern)

    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {file_pattern}")

    # Parse timestamps and sort files
    file_data: List[Tuple[int, str]] = []
    for f in files:
        try:
            timestamp = parse_timestamp_from_filename(os.path.basename(f))
            file_data.append((timestamp, f))
        except ValueError:
            continue  # Skip files that don't match the expected pattern

    if not file_data:
        raise ValueError("No valid tNR files found with parseable timestamps")

    # Sort by timestamp
    file_data.sort(key=lambda x: x[0])

    timestamps = np.array([t for t, _ in file_data])

    # Load reference data (first file, t=0)
    _, ref_file = file_data[0]
    q_ref, r_ref, dr_ref = load_tnr_data(ref_file)

    # Calculate distance distance for each file
    distance_values = np.zeros(len(file_data))

    for i, (_, filepath) in enumerate(file_data):
        q, r, dr = load_tnr_data(filepath)

        # Verify Q arrays match
        if not np.allclose(q, q_ref):
            raise ValueError(f"Q array mismatch in file: {filepath}")

        distance_values[i] = np.mean((r_ref - r) ** 2 / (dr**2 + dr_ref**2))

    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(timestamps[1:-1], distance_values[1:-1], "o-", markersize=4, linewidth=1)
    plt.xlabel("Time (s)", fontsize=12)
    plt.ylabel("$\\chi^2$ from reference", fontsize=12)
    plt.title(f"Run {t_nr_run} -- $\\chi^2$ vs time relative to t=0", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_plot:
        plt.savefig(output_plot, dpi=150)
        print(f"Plot saved to: {output_plot}")

    plt.show()

    return timestamps, distance_values


if __name__ == "__main__":
    t_nr_run = "218389"
    #t_nr_run = "218396"
    data_directory = f"/home/mat/ORNL/expt11-tNR-{t_nr_run}-120s"

    timestamps, distance_values = compute_distance_time_series(
        data_directory, output_plot=f"tnr_{t_nr_run}_chi2_series.png",
        t_nr_run=t_nr_run,
    )

    print(f"\nAnalyzed {len(timestamps)} time points")
    print(f"Time range: {timestamps[0]} - {timestamps[-1]} s")
    print(f"Distance range: {distance_values.min():.6f} - {distance_values.max():.6f}")
