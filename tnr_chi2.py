"""
CLI for assessing temporal change in time-resolved neutron reflectometry
(tNR) data using several complementary analyses:

  1. Running chi-squared per interval, with per-type reference.
  2. Q-resolved standardized-residual heatmap.
  3. RQ^4 representation heatmap (log scale).
  4. PCA of the R(Q,t) matrix; first 3 PC time-loadings + Q-mode shapes.
  5. Symmetric KL divergence per Q (Gaussian per bin), per-type reference.

Ordering and timestamps come from the reduction JSON sidecar
``r<run>_eis_reduction.json`` produced alongside the reduced ASCII files.
Each interval entry's ``label`` maps to ``r<run>_<label>.txt``.

Pass ``--out-dir DIR`` to emit all analyses (PNGs + ASCII tables) using
``--label`` as the filename prefix. Otherwise pass individual ``--*-plot``
or ``--*-data`` options to selectively write outputs.
"""

import json
import os
import sys
from datetime import datetime
from typing import List, Tuple

import click
import matplotlib.pyplot as plt
import numpy as np


HOLD_COLOR = "#1f77b4"  # blue
EIS_COLOR = "#d62728"  # red
OTHER_COLOR = "#7f7f7f"  # grey


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_tnr_data(filepath: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load (Q, R, dR) from a tNR ASCII file."""
    data = np.loadtxt(filepath, usecols=(0, 1, 2))
    return data[:, 0], data[:, 1], data[:, 2]


def align_to_ref(
    q_ref: np.ndarray, q: np.ndarray, r: np.ndarray, dr: np.ndarray, rtol: float = 1e-4
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Align (q, r, dr) onto q_ref by nearest match within rtol relative tol.

    Returns (r_aligned, dr_aligned, mask) on the q_ref grid.
    """
    if len(q) == len(q_ref) and np.allclose(q, q_ref):
        return r, dr, np.ones(len(q_ref), dtype=bool)

    idx = np.searchsorted(q, q_ref)
    idx_left = np.clip(idx - 1, 0, len(q) - 1)
    idx_right = np.clip(idx, 0, len(q) - 1)
    diff_left = np.abs(q[idx_left] - q_ref)
    diff_right = np.abs(q[idx_right] - q_ref)
    use_left = diff_left < diff_right
    nearest = np.where(use_left, idx_left, idx_right)
    mindiff = np.minimum(diff_left, diff_right)
    mask = mindiff <= rtol * np.maximum(np.abs(q_ref), 1e-12)
    return r[nearest], dr[nearest], mask


def find_reduction_json(data_dir: str) -> str:
    candidates = [f for f in os.listdir(data_dir) if f.endswith("_eis_reduction.json")]
    if not candidates:
        raise FileNotFoundError(f"No *_eis_reduction.json file found in {data_dir}")
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple reduction JSON files found in {data_dir}: {candidates}. "
            "Specify one explicitly with --json."
        )
    return os.path.join(data_dir, candidates[0])


def load_run(data_dir: str, json_path: str, min_duration: float = 0.0):
    """Load all qualifying intervals onto a common Q grid.

    Returns
    -------
    times : (T,) float array — seconds since first interval's start
    types : list[str] of length T
    labels : list[str] of length T
    q : (Nq,) float array — reference Q grid (from first qualifying interval)
    R, dR : (T, Nq) arrays — aligned reflectivity and uncertainty
    valid : (T, Nq) bool array — True where a real measurement exists
    """
    with open(json_path) as fh:
        meta = json.load(fh)
    run_number = meta["run_number"]
    intervals = meta["intervals"]
    if not intervals:
        raise ValueError(f"No intervals listed in {json_path}")

    t0 = datetime.fromisoformat(intervals[0]["start"])

    times: List[float] = []
    types: List[str] = []
    labels: List[str] = []
    R_rows: List[np.ndarray] = []
    dR_rows: List[np.ndarray] = []
    valid_rows: List[np.ndarray] = []

    q_ref = None
    missing: List[str] = []
    n_short = 0

    for interval in intervals:
        label = interval["label"]
        itype = interval.get("interval_type", "unknown")
        start = datetime.fromisoformat(interval["start"])
        end = datetime.fromisoformat(interval["end"])
        duration = (end - start).total_seconds()
        if duration < min_duration:
            n_short += 1
            continue

        filename = f"r{run_number}_{label}.txt"
        filepath = os.path.join(data_dir, filename)
        if not os.path.isfile(filepath):
            missing.append(filename)
            continue

        q, r, dr = load_tnr_data(filepath)
        if q_ref is None:
            q_ref = q
            r_a, dr_a, mask = r.copy(), dr.copy(), np.ones(len(q_ref), dtype=bool)
        else:
            r_a, dr_a, mask = align_to_ref(q_ref, q, r, dr)
            if not mask.any():
                raise ValueError(f"No overlapping Q points with reference for: {filepath}")
            # Zero out invalid entries; valid mask records them.
            r_a = np.where(mask, r_a, 0.0)
            dr_a = np.where(mask, dr_a, 0.0)

        times.append((start - t0).total_seconds())
        types.append(itype)
        labels.append(label)
        R_rows.append(r_a)
        dR_rows.append(dr_a)
        valid_rows.append(mask & (dr_a > 0) if mask.any() else mask)

    if n_short:
        click.echo(f"Skipped {n_short} interval(s) with duration < {min_duration} s", err=True)
    if missing:
        click.echo(
            f"Warning: {len(missing)} interval file(s) listed in JSON not found "
            f"(first: {missing[0]})",
            err=True,
        )
    if not times:
        raise FileNotFoundError(
            f"None of the interval files referenced by {json_path} were found in {data_dir}"
        )

    return (
        np.array(times),
        types,
        labels,
        q_ref,
        np.vstack(R_rows),
        np.vstack(dR_rows),
        np.vstack(valid_rows),
    )


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def plot_by_type(ax, times, values, types, ylabel, title):
    """Scatter values vs times with per-type coloring + grey connecting line."""
    ax.plot(times, values, "-", color="lightgrey", linewidth=1, zorder=1)
    types_arr = np.array(types)
    for itype, color in [("hold", HOLD_COLOR), ("eis", EIS_COLOR)]:
        m = types_arr == itype
        if m.any():
            ax.plot(times[m], values[m], "o", color=color, markersize=5, label=itype, zorder=2)
    other = ~np.isin(types_arr, ["hold", "eis"])
    if other.any():
        ax.plot(times[other], values[other], "o", color=OTHER_COLOR, markersize=5,
                label="other", zorder=2)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()


# ---------------------------------------------------------------------------
# (1) chi-squared per interval, per-type reference
# ---------------------------------------------------------------------------

def analyze_chi2(times, types, R, dR, valid):
    """Per-interval chi-squared with paired-difference denominator.

    Returns chi2[t] = mean( (R - R_ref)^2 / (dR^2 + dR_ref^2) ) over valid Q.
    For two independent measurements of the same truth this has expectation 1.
    """
    types_arr = np.array(types)
    chi2 = np.full(len(times), np.nan)
    for itype in np.unique(types_arr):
        idx = np.where(types_arr == itype)[0]
        ref = idx[0]
        for i in idx:
            v = valid[i] & valid[ref]
            if not v.any():
                continue
            denom = dR[i, v] ** 2 + dR[ref, v] ** 2
            chi2[i] = float(np.mean((R[ref, v] - R[i, v]) ** 2 / denom))
    return chi2


def analyze_delta2(times, types, R, dR, valid):
    """Noise-corrected fractional change metric, comparable across count times.

    delta2[t] = mean( max(0, (R - R_ref)^2 - dR^2 - dR_ref^2) / R_ref^2 )

    The noise subtraction in the numerator removes the statistical floor; the
    R_ref^2 denominator makes the result dimensionless. delta is roughly the
    fractional RMS change of R relative to its reference value.
    """
    types_arr = np.array(types)
    delta2 = np.full(len(times), np.nan)
    for itype in np.unique(types_arr):
        idx = np.where(types_arr == itype)[0]
        ref = idx[0]
        for i in idx:
            v = valid[i] & valid[ref] & (R[ref] != 0)
            if not v.any():
                continue
            num = (R[ref, v] - R[i, v]) ** 2 - dR[i, v] ** 2 - dR[ref, v] ** 2
            num = np.maximum(num, 0.0)
            delta2[i] = float(np.mean(num / R[ref, v] ** 2))
    return delta2


def write_chi2_ascii(path, times, chi2, delta2, types, labels, data_dir, json_path):
    with open(path, "w") as fh:
        fh.write("# tNR chi-squared and noise-corrected delta^2 (per-type reference)\n")
        fh.write("# chi2  = mean((R-R_ref)^2 / (dR^2 + dR_ref^2))     (baseline 1)\n")
        fh.write("# delta2 = mean(max(0, (R-R_ref)^2 - dR^2 - dR_ref^2) / R_ref^2)\n")
        fh.write(f"# source_dir: {data_dir}\n# json: {json_path}\n# n_points: {len(times)}\n")
        fh.write("# time_s\tchi2\tdelta2\tinterval_type\tlabel\n")
        for t, c, d, ity, lab in zip(times, chi2, delta2, types, labels):
            fh.write(f"{t:.3f}\t{c:.8e}\t{d:.8e}\t{ity}\t{lab}\n")


def plot_split_by_type(times, values, types, ylabel, title, out_path,
                       hline=None):
    """Stacked subplot, one per interval_type present, sharing the time axis."""
    types_arr = np.array(types)
    present = [t for t in ["hold", "eis"] if (types_arr == t).any()]
    others = sorted(set(types_arr) - {"hold", "eis"})
    present += others
    if not present:
        present = ["all"]

    fig, axes = plt.subplots(len(present), 1, figsize=(10, 3 + 2.2 * len(present)),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    color_map = {"hold": HOLD_COLOR, "eis": EIS_COLOR}
    for ax, itype in zip(axes, present):
        m = types_arr == itype if itype != "all" else np.ones(len(times), dtype=bool)
        color = color_map.get(itype, OTHER_COLOR)
        ax.plot(times[m], values[m], "-", color="lightgrey", linewidth=1, zorder=1)
        ax.plot(times[m], values[m], "o", color=color, markersize=5, zorder=2,
                label=itype)
        if hline is not None:
            ax.axhline(hline, color="k", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)
    axes[0].set_title(title)
    axes[-1].set_xlabel("Time (s)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# (2) Standardized-residual heatmap
# ---------------------------------------------------------------------------

def plot_residual_heatmap(times, q, R, dR, valid, label, out_path):
    """z(t,Q) = (R(t,Q) - R_ref(Q)) / dR(t,Q), ref = first interval."""
    ref = 0
    z = np.zeros_like(R)
    common = valid[ref][None, :] & valid
    safe_dR = np.where(dR > 0, dR, np.nan)
    z = (R[ref][None, :] - R) / safe_dR
    z = np.where(common, z, np.nan)

    fig, ax = plt.subplots(figsize=(10, 6))
    # pcolormesh expects edges; use cell-centered approximation
    qmesh = np.concatenate([[q[0] - (q[1] - q[0]) / 2],
                            (q[:-1] + q[1:]) / 2,
                            [q[-1] + (q[-1] - q[-2]) / 2]])
    if len(times) > 1:
        dt = np.diff(times)
        tmesh = np.concatenate([[times[0] - dt[0] / 2],
                                times[:-1] + dt / 2,
                                [times[-1] + dt[-1] / 2]])
    else:
        tmesh = np.array([times[0] - 0.5, times[0] + 0.5])

    vmax = float(np.nanpercentile(np.abs(z), 98))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    pcm = ax.pcolormesh(tmesh, qmesh, z.T, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        shading="auto")
    ax.set_yscale("log")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Q (Å$^{-1}$)")
    title = f"Run {label} — standardized residual (R − R₀)/dR" if label else \
            "standardized residual (R − R₀)/dR"
    ax.set_title(title)
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("(R − R$_0$) / dR")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# (3) RQ^4 heatmap
# ---------------------------------------------------------------------------

def plot_rq4_heatmap(times, q, R, valid, label, out_path):
    """log10(R * Q^4) vs (time, Q)."""
    rq4 = R * (q[None, :] ** 4)
    with np.errstate(invalid="ignore", divide="ignore"):
        log_rq4 = np.where((rq4 > 0) & valid, np.log10(rq4), np.nan)

    fig, ax = plt.subplots(figsize=(10, 6))
    qmesh = np.concatenate([[q[0] - (q[1] - q[0]) / 2],
                            (q[:-1] + q[1:]) / 2,
                            [q[-1] + (q[-1] - q[-2]) / 2]])
    if len(times) > 1:
        dt = np.diff(times)
        tmesh = np.concatenate([[times[0] - dt[0] / 2],
                                times[:-1] + dt / 2,
                                [times[-1] + dt[-1] / 2]])
    else:
        tmesh = np.array([times[0] - 0.5, times[0] + 0.5])

    pcm = ax.pcolormesh(tmesh, qmesh, log_rq4.T, cmap="viridis", shading="auto")
    ax.set_yscale("log")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Q (Å$^{-1}$)")
    title = f"Run {label} — log$_{{10}}$(R·Q$^4$)" if label else "log$_{10}$(R·Q$^4$)"
    ax.set_title(title)
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("log$_{10}$(R·Q$^4$)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# (4) PCA
# ---------------------------------------------------------------------------

def analyze_pca(times, q, R, valid, n_components=3):
    """SVD-based PCA on log10(R), restricted to Q points common to all intervals.

    Operating on log10(R) prevents the low-Q (large R) region from dominating.

    Returns (scores [T, k], components [k, Nq_common], explained_variance_ratio,
    q_common, mean_curve [Nq_common]).
    """
    common = np.all(valid, axis=0) & np.all(R > 0, axis=0)
    if common.sum() < 3:
        raise ValueError("Too few Q points common to all intervals for PCA.")
    q_c = q[common]
    M = np.log10(R[:, common])
    mean = M.mean(axis=0)
    Mc = M - mean
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    k = min(n_components, len(S))
    scores = U[:, :k] * S[:k]
    components = Vt[:k]
    var = (S ** 2) / max(len(times) - 1, 1)
    evr = var[:k] / var.sum()
    return scores, components, evr, q_c, mean


def plot_pca(times, types, scores, components, evr, q_c, label, out_path):
    k = scores.shape[1]
    fig, axes = plt.subplots(k + 1, 1, figsize=(10, 3 + 2.2 * k), sharex=False)
    types_arr = np.array(types)

    # One subplot per component for clear time scores
    for j in range(k):
        ax = axes[j]
        ax.plot(times, scores[:, j], "-", color="lightgrey", linewidth=1, zorder=1)
        for itype, color in [("hold", HOLD_COLOR), ("eis", EIS_COLOR)]:
            m = types_arr == itype
            if m.any():
                ax.plot(times[m], scores[m, j], "o", color=color, markersize=4,
                        label=itype if j == 0 else None, zorder=2)
        ax.set_ylabel(f"PC{j+1} score\n({evr[j]*100:.1f}%)")
        ax.grid(True, alpha=0.3)
        if j == 0:
            ax.legend(fontsize=8, loc="best")
            ax.set_title(f"Run {label} — PCA time scores" if label else "PCA time scores")
        if j == k - 1:
            ax.set_xlabel("Time (s)")

    # Q-mode shapes
    ax = axes[-1]
    for j in range(k):
        ax.plot(q_c, components[j], label=f"PC{j+1} ({evr[j]*100:.1f}%)")
    ax.set_xscale("log")
    ax.set_xlabel("Q (Å$^{-1}$)")
    ax.set_ylabel("Component\namplitude (log$_{10}$R)")
    ax.set_title("PCA Q-mode shapes")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_pca_ascii(path, times, types, labels, scores, evr, data_dir, json_path):
    with open(path, "w") as fh:
        fh.write("# tNR PCA time scores\n")
        fh.write(f"# source_dir: {data_dir}\n# json: {json_path}\n")
        fh.write(f"# n_points: {len(times)}  n_components: {scores.shape[1]}\n")
        fh.write("# explained_variance_ratio: " + " ".join(f"{x:.6f}" for x in evr) + "\n")
        cols = "\t".join(f"PC{j+1}" for j in range(scores.shape[1]))
        fh.write(f"# time_s\t{cols}\tinterval_type\tlabel\n")
        for i, (t, ity, lab) in enumerate(zip(times, types, labels)):
            row = "\t".join(f"{scores[i, j]:.6e}" for j in range(scores.shape[1]))
            fh.write(f"{t:.3f}\t{row}\t{ity}\t{lab}\n")


# ---------------------------------------------------------------------------
# (5) Symmetric KL divergence (per-type ref, Gaussian per Q bin)
# ---------------------------------------------------------------------------

def analyze_kl(times, types, R, dR, valid):
    """Symmetric KL between Gaussian product distributions, per-type reference.

    For two Gaussians N(μ1, σ1), N(μ2, σ2):
        KL(1||2) = log(σ2/σ1) + (σ1^2 + (μ1-μ2)^2)/(2 σ2^2) - 1/2
    Symmetric KL summed over Q bins, divided by N to give per-bin value.
    """
    types_arr = np.array(types)
    kl = np.full(len(times), np.nan)
    for itype in np.unique(types_arr):
        idx = np.where(types_arr == itype)[0]
        ref = idx[0]
        mu_r = R[ref]
        sig_r = dR[ref]
        for i in idx:
            v = valid[i] & valid[ref] & (sig_r > 0) & (dR[i] > 0)
            if not v.any():
                continue
            mu_a, sig_a = R[i, v], dR[i, v]
            mu_b, sig_b = mu_r[v], sig_r[v]
            kl_ab = np.log(sig_b / sig_a) + (sig_a ** 2 + (mu_a - mu_b) ** 2) / (2 * sig_b ** 2) - 0.5
            kl_ba = np.log(sig_a / sig_b) + (sig_b ** 2 + (mu_b - mu_a) ** 2) / (2 * sig_a ** 2) - 0.5
            kl[i] = float(np.mean(0.5 * (kl_ab + kl_ba)))
    return kl


def write_kl_ascii(path, times, kl, types, labels, data_dir, json_path):
    with open(path, "w") as fh:
        fh.write("# tNR symmetric KL divergence per Q bin (per-type reference)\n")
        fh.write(f"# source_dir: {data_dir}\n# json: {json_path}\n# n_points: {len(times)}\n")
        fh.write("# time_s\tsym_kl\tinterval_type\tlabel\n")
        for t, k, ity, lab in zip(times, kl, types, labels):
            fh.write(f"{t:.3f}\t{k:.8e}\t{ity}\t{lab}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("data_dir", type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True))
@click.option("--json", "json_path", type=click.Path(exists=True, dir_okay=False, readable=True),
              default=None, help="Reduction JSON (defaults to *_eis_reduction.json in DATA_DIR).")
@click.option("--out-dir", type=click.Path(file_okay=False, writable=True), default=None,
              help="If set, run all analyses and write all outputs here using --label as prefix.")
@click.option("--label", default="", help="Label/prefix for filenames and plot titles.")
@click.option("--min-duration", type=float, default=5.0, show_default=True,
              help="Skip intervals shorter than this many seconds.")
@click.option("--no-show", is_flag=True, default=False, help="Do not display plots interactively.")
# Individual outputs (used when --out-dir is not given)
@click.option("--chi2-plot", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--chi2-data", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--delta2-plot", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--heatmap-plot", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--rq4-plot", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--pca-plot", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--pca-data", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--kl-plot", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--kl-data", type=click.Path(dir_okay=False, writable=True), default=None)
def main(data_dir, json_path, out_dir, label, min_duration, no_show,
         chi2_plot, chi2_data, delta2_plot, heatmap_plot, rq4_plot,
         pca_plot, pca_data, kl_plot, kl_data):
    """Assess temporal change in tNR data via multiple complementary analyses."""
    try:
        if json_path is None:
            json_path = find_reduction_json(data_dir)
        times, types, labels, q, R, dR, valid = load_run(
            data_dir, json_path, min_duration=min_duration
        )
    except (FileNotFoundError, ValueError, KeyError) as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    n_hold = sum(1 for t in types if t == "hold")
    n_eis = sum(1 for t in types if t == "eis")
    click.echo(f"Loaded {len(times)} intervals from {data_dir}")
    click.echo(f"  hold: {n_hold}  eis: {n_eis}  other: {len(times) - n_hold - n_eis}")
    click.echo(f"  time range: {times[0]:.1f} - {times[-1]:.1f} s")
    click.echo(f"  Q range:    {q.min():.4g} - {q.max():.4g}  (N_Q = {len(q)})")

    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        prefix = (label + "_") if label else ""
        chi2_plot = chi2_plot or os.path.join(out_dir, f"{prefix}chi2.png")
        chi2_data = chi2_data or os.path.join(out_dir, f"{prefix}chi2.txt")
        delta2_plot = delta2_plot or os.path.join(out_dir, f"{prefix}delta2.png")
        heatmap_plot = heatmap_plot or os.path.join(out_dir, f"{prefix}residual_heatmap.png")
        rq4_plot = rq4_plot or os.path.join(out_dir, f"{prefix}rq4_heatmap.png")
        pca_plot = pca_plot or os.path.join(out_dir, f"{prefix}pca.png")
        pca_data = pca_data or os.path.join(out_dir, f"{prefix}pca_scores.txt")
        kl_plot = kl_plot or os.path.join(out_dir, f"{prefix}kl.png")
        kl_data = kl_data or os.path.join(out_dir, f"{prefix}kl.txt")

    # --- chi^2 + delta^2 ---
    chi2 = analyze_chi2(times, types, R, dR, valid)
    delta2 = analyze_delta2(times, types, R, dR, valid)
    click.echo(f"  chi^2 range:   {np.nanmin(chi2):.4g} - {np.nanmax(chi2):.4g}")
    click.echo(f"  delta^2 range: {np.nanmin(delta2):.4g} - {np.nanmax(delta2):.4g}")
    if chi2_data:
        write_chi2_ascii(chi2_data, times, chi2, delta2, types, labels, data_dir, json_path)
        click.echo(f"  wrote {chi2_data}")
    if chi2_plot:
        plot_split_by_type(
            times, chi2, types,
            "$\\chi^2$",
            f"Run {label} — $\\chi^2$ vs time (per type)" if label
            else "$\\chi^2$ vs time (per type)",
            chi2_plot, hline=1.0,
        )
        click.echo(f"  wrote {chi2_plot}")
    if delta2_plot:
        plot_split_by_type(
            times, delta2, types,
            "$\\delta^2$ (fractional)",
            f"Run {label} — noise-corrected fractional change" if label
            else "noise-corrected fractional change",
            delta2_plot, hline=0.0,
        )
        click.echo(f"  wrote {delta2_plot}")

    # --- residual heatmap ---
    if heatmap_plot:
        plot_residual_heatmap(times, q, R, dR, valid, label, heatmap_plot)
        click.echo(f"  wrote {heatmap_plot}")

    # --- RQ^4 heatmap ---
    if rq4_plot:
        plot_rq4_heatmap(times, q, R, valid, label, rq4_plot)
        click.echo(f"  wrote {rq4_plot}")

    # --- PCA ---
    if pca_plot or pca_data:
        try:
            scores, components, evr, q_c, _ = analyze_pca(times, q, R, valid, n_components=3)
        except ValueError as exc:
            click.echo(f"PCA skipped: {exc}", err=True)
        else:
            click.echo(f"  PCA explained variance: " + ", ".join(f"{e*100:.1f}%" for e in evr))
            if pca_plot:
                plot_pca(times, types, scores, components, evr, q_c, label, pca_plot)
                click.echo(f"  wrote {pca_plot}")
            if pca_data:
                write_pca_ascii(pca_data, times, types, labels, scores, evr, data_dir, json_path)
                click.echo(f"  wrote {pca_data}")

    # --- KL ---
    if kl_plot or kl_data:
        kl = analyze_kl(times, types, R, dR, valid)
        click.echo(f"  KL range:    {np.nanmin(kl):.4g} - {np.nanmax(kl):.4g}")
        if kl_data:
            write_kl_ascii(kl_data, times, kl, types, labels, data_dir, json_path)
            click.echo(f"  wrote {kl_data}")
        if kl_plot:
            fig, ax = plt.subplots(figsize=(10, 6))
            plot_by_type(ax, times, kl, types,
                         "Symmetric KL per Q",
                         f"Run {label} — symmetric KL vs time" if label else "symmetric KL vs time")
            fig.tight_layout()
            fig.savefig(kl_plot, dpi=150)
            click.echo(f"  wrote {kl_plot}")
            if not no_show:
                plt.show()
            plt.close(fig)


if __name__ == "__main__":
    main()
