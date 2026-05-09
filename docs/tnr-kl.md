# Symmetric KL divergence of time-resolved neutron reflectometry data

This note explains the KL divergence outputs produced by
[tnr_chi2.py](../tnr_chi2.py) (`*_kl.png` and `*_kl.txt`) and how to read
them alongside the running $\chi^2$ output.

## Motivation

The running $\chi^2$ statistic answers *"how unlikely is this measurement
under the assumption that nothing has changed?"* It mixes the size of the
deviation with the inverse measurement variance, so a noisy short interval
inflates $\chi^2$ even when the underlying reflectivity has not moved.

The Kullback–Leibler (KL) divergence asks the more symmetric question
*"how distinguishable are these two measurements as probability
distributions?"* When both intervals are noisy, KL stays small; when the
mean reflectivity actually shifts beyond the noise level, KL grows. It
therefore complements $\chi^2$ as a change detector.

## What is computed

For each Q bin we treat the reflectivity as a Gaussian random variable

$$
R(Q) \sim \mathcal{N}\!\bigl(\mu(Q),\, \sigma(Q)^2\bigr)
$$

with $\mu = R$ (measured reflectivity) and $\sigma = dR$ (reported
uncertainty). Q bins are assumed independent, which lets us reason per-Q.

For two Gaussians $p = \mathcal{N}(\mu_1, \sigma_1^2)$ and
$q = \mathcal{N}(\mu_2, \sigma_2^2)$ the closed-form KL divergence is

$$
\mathrm{KL}(p \,\|\, q) \;=\; \log\!\frac{\sigma_2}{\sigma_1}
\;+\; \frac{\sigma_1^2 + (\mu_1 - \mu_2)^2}{2\,\sigma_2^2}
\;-\; \tfrac{1}{2}.
$$

KL is **asymmetric**, so we use the **symmetric** (Jensen) form

$$
\mathrm{KL}_{\text{sym}}(p, q) \;=\; \tfrac{1}{2}\bigl[\,\mathrm{KL}(p\,\|\,q) + \mathrm{KL}(q\,\|\,p)\,\bigr].
$$

For the reported scalar at time $t$ we average $\mathrm{KL}_{\text{sym}}$
across all valid Q bins (those with $dR > 0$ in both the current interval
and its reference):

$$
\overline{\mathrm{KL}}(t) \;=\; \frac{1}{N_Q^{\mathrm{valid}}}\sum_i \mathrm{KL}_{\text{sym}}\!\bigl(p_i(t),\, q_i(\text{ref})\bigr).
$$

The reference $q$ is **per-type**: the first qualifying interval of each
`interval_type` (e.g. the first `hold` is the reference for all subsequent
holds; the first `eis` for all subsequent eis intervals). This matches the
$\chi^2$ analysis and means each type starts at $\overline{\mathrm{KL}} = 0$.

## How to read the plot

`*_kl.png` plots $\overline{\mathrm{KL}}$ vs time on a single panel:

- **Grey line** connects all intervals in time order.
- **Blue circles** = `hold` intervals.
- **Red circles** = `eis` intervals.
- Reference points sit at exactly zero by construction.

Units are *nats per Q bin* (the per-bin average of a sum of KL terms
expressed in natural-log units). The absolute value is less informative than
the trend; what matters is the shape and how it compares to the
$\chi^2$ panel.

## Practical interpretation

- **KL ≈ $\chi^2 / 2$ (linear scaling)** → the mean shift dominates and
  the noise levels of the two intervals are comparable. Both metrics tell
  the same story; you can use either.
- **KL flat while $\chi^2$ rises** → the $\chi^2$ growth is being driven
  by the *noise term* (smaller $dR$ in newer intervals, or short-duration
  reference) rather than a real mean shift. This is the diagnostic the
  $\chi^2$ artifact discussion was about — KL stays honest because of the
  $\log(\sigma_2/\sigma_1)$ term.
- **KL rises while $\chi^2$ stays flat** → the *noise level* changed (e.g.
  a degraded reduction or a different counting time) without the mean
  shifting. Worth checking the reduction sidecars.
- **Step-wise KL changes between intervals** → real structural change
  rather than drift. Localize the step in time and inspect the
  residual heatmap and PCA at that interval.
- **Sustained linear KL drift** → continuous structural evolution; PCA's
  PC1 should show a matching monotonic time score.

## Caveats

- The Gaussian-per-Q-bin assumption is exact for the pointwise reflectivity
  values returned by the reduction and is a very good approximation in the
  high-count regime. It is **wrong** in the low-count tails (high-Q),
  where the count distribution is closer to Poisson; KL there can be
  optimistic.
- Q bins are treated as **independent**. Real reduced reflectivity is
  correlated across neighboring Q bins (resolution and binning), so the
  per-Q-bin sum overstates the total information content. We average rather
  than sum, which mitigates but does not remove this.
- Like the $\chi^2$ analysis, KL is sensitive to short noisy intervals.
  The `--min-duration` option (default 5 s) skips obviously truncated
  sub-intervals.
- KL is **always non-negative** and equals zero only at the reference. It
  cannot tell you the *direction* of the change. Use the residual heatmap
  or PCA scores for signed information.
- KL ignores Q-resolution effects. For comparing measurements with very
  different $dQ$ (e.g. different angles), align reductions onto a common
  resolution first.

## ASCII output format

`*_kl.txt` contains:

```
# tNR symmetric KL divergence per Q bin (per-type reference)
# source_dir: <DATA_DIR>
# json: <reduction JSON path>
# n_points: T
# time_s    sym_kl  interval_type   label
```

followed by one tab-separated row per interval. Time is seconds since the
first qualifying interval's start. `sym_kl` is the per-Q-bin mean of the
symmetric KL divergence (in nats); reference rows are exactly `0`.
Interval types are typically `hold` or `eis`; labels match the JSON
entries.
