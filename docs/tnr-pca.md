# PCA of time-resolved neutron reflectometry data

This note explains the PCA outputs produced by [tnr_chi2.py](../tnr_chi2.py)
(`*_pca.png` and `*_pca_scores.txt`) and how to read them.

## What is being decomposed

For a tNR run with `T` intervals measured on a common Q grid of `N_Q` points,
we form the data matrix

$$
M_{t,i} = \log_{10} R(t_i, Q_i)
$$

where rows are time intervals and columns are Q points. Working on
$\log_{10}R$ (rather than $R$ directly) is essential: $R$ spans many decades
between low and high Q, so a linear-$R$ PCA would be dominated by the
brightest low-$Q$ region and miss subtle structural changes elsewhere.

Only the Q points that are valid in **every** interval are used (after the
nearest-Q alignment performed in `align_to_ref`), so the matrix is
rectangular with no missing entries.

The matrix is mean-centered along time (i.e. for each Q we subtract the
average $\log_{10}R(Q)$ over all intervals) and then decomposed by SVD:

$$
M - \bar{M} = U \, S \, V^{\!T}
$$

We retain the first $k=3$ components and report:

- **Components** $V^{\!T}_j$ — vectors of length $N_Q$ describing how
  $\log_{10}R$ varies with $Q$ along PC$_j$. Plotted in the bottom panel.
- **Scores** $U_j S_j$ — vectors of length $T$ giving the amplitude of PC$_j$
  at each time. Plotted in the upper panels (one panel per PC).
- **Explained variance ratio (EVR)** — the fraction of total variance
  captured by each PC, $S_j^2 / \sum_l S_l^2$.

## How to read the plots

`*_pca.png` has $k+1$ panels stacked vertically:

1. **PC1 score vs time** — usually the dominant change axis. The legend marks
   the EVR percentage.
2. **PC2 score vs time** — second most important orthogonal mode.
3. **PC3 score vs time** — third mode.
4. **Q-mode shapes** — the components $V^{\!T}_j$ vs $Q$ (log-Q axis).

Hold and EIS intervals are colored differently in the time-score plots so you
can see whether the two interval types follow the same trajectory or split.

## Practical interpretation

- **A single dominant PC1** (e.g. EVR ≫ remaining PCs) with a smooth
  monotonic time score → the run is changing along **one effective
  reaction coordinate**. The shape of PC1 vs Q tells you *what* is changing:
  - amplitude concentrated at **low Q (≲ 0.03 Å⁻¹)** → thick or surface
    layers (e.g. ionomer swelling, electrolyte penetration).
  - amplitude at **high Q** → buried interfaces with short
    length scales (e.g. a Cu/Si oxide layer growing/dissolving).
  - oscillatory shape with a node spacing $\Delta Q$ → a thickness change of
    order $\pi/\Delta Q$ (Kiessig-fringe shift).
- **Two strong PCs of comparable EVR** → at least two independent processes
  are evolving (e.g. ionomer swelling at low Q and oxide growth at high Q).
  Plot PC2 score vs PC1 score to see whether they are coupled.
- **All EVRs small (each ~few %)** → no single coherent change pattern emerges
  above noise. The run is dominated by counting statistics, not structure.
- **PC1 jumps step-wise across EIS intervals** → the EIS measurement is
  introducing a systematic offset (different probe state) rather than tracking
  the same physics as the holds.

## Caveats

- PCA is **linear** and **unsupervised**: the PCs are mathematical axes of
  maximum variance, not physical layer parameters. They are diagnostic, not
  explanatory. Use them to decide whether (and where in Q) to focus a
  refl1d/bumps fit.
- The **mean** subtracted before SVD is not the reference curve from
  `analyze_chi2`; it is the per-Q time-average. PC scores are therefore signed
  deviations from that average.
- PCA does not propagate $dR$ uncertainties. Q points with very small
  reflectivity (large relative error) can still drive a component. If this is
  a concern, weight rows or columns inversely by their typical $dR/R$ before
  decomposition.
- The number of nonzero singular values is at most $\min(T, N_Q^{\text{common}})$.
  If you have more Q points than intervals, only $T-1$ PCs carry information
  after mean-centering.

## ASCII output format

`*_pca_scores.txt` contains, in order:

```
# tNR PCA time scores
# source_dir: <DATA_DIR>
# json: <reduction JSON path>
# n_points: T  n_components: k
# explained_variance_ratio: e1 e2 e3
# time_s    PC1     PC2     PC3     interval_type   label
```

followed by one tab-separated row per interval. Time is seconds since the
first qualifying interval's start. Interval types are typically `hold` or
`eis`; labels match the JSON entries.
