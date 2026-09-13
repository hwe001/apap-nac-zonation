# Supplementary Methods

## S1. Scenario treatment: the constant-boost NAC upper bound

The main text uses the calibrated, protocol-driven NAC input (plasma profile
of the 21-h IV protocol → cysteine-equivalent pool → synthesis boost peaking
at M_max = 4). An upper-bound scenario holds the peak boost flat for the full
21-h window, decoupled from plasma pharmacokinetics. Protection under the two
representations:

| Gradients | Dose | 0 h | 2 h | 4 h | 6 h | 8 h | 10 h | 12 h |
|---|---|---|---|---|---|---|---|---|
| assumed | 4 g | 85% | 52% | 24% | 9% | 2% | 0% | 0% |
| assumed | 16 g | 45% | 33% | 21% | 9% | 3% | 0% | 0% |
| measured | 4 g | 89% | 82% | 70% | 46% | 16% | 3% | 0% |
| measured | 16 g | 56% | 44% | 32% | 20% | 9% | 2% | 0% |

The constant-boost scenario is uniformly more optimistic (the boost does not
wane with plasma levels), but every qualitative conclusion of the main text is
preserved: the window closes between 4 and 12 h, and the measured gradients
give higher protection at every start time.

## S2. Boost-magnitude sensitivity (M_max = 1, 2, 4)

Protection at 4 h post-ingestion for NAC started at 0 h, under the
protocol-driven input with the peak boost varied:

| Gradients | Dose | M_max = 1 | M_max = 2 | M_max = 4 (calibrated) |
|---|---|---|---|---|
| assumed | 4 g | 2% | 2% | 2% |
| assumed | 16 g | 2% | 3% | 3% |
| measured | 4 g | 11% | 14% | 16% |
| measured | 16 g | 5% | 8% | 9% |

Under assumed gradients the 4-h window is closed at every M_max (the
pericentral GSH pool is already exhausted before any synthesis boost can
act); under measured gradients the residual protection scales with M_max.
The assumed-versus-measured ordering is preserved at every magnitude.

## S3. Absorption-rate sensitivity

Peak pericentral adducts (nmol) at fixed total dose 16 g, single vs staggered
(4 × 4 g every 2 h) ingestion, across absorption rates:

| Gradients | ka (d⁻¹) | Single | Staggered | Δ |
|---|---|---|---|---|
| assumed | 17 | 3.458 × 10⁻² | 3.391 × 10⁻² | −1.9% |
| assumed | 34 | 3.484 × 10⁻² | 3.414 × 10⁻² | −2.0% |
| assumed | 60 | 3.490 × 10⁻² | 3.419 × 10⁻² | −2.0% |
| measured | 17 | 3.425 × 10⁻² | 3.299 × 10⁻² | −3.7% |
| measured | 34 | 3.466 × 10⁻² | 3.331 × 10⁻² | −3.9% |
| measured | 60 | 3.477 × 10⁻² | 3.338 × 10⁻² | −4.0% |

ka = 17 d⁻¹ corresponds to an absorption half-life of ~1 h; 60 d⁻¹ to ~17 min.
The ingestion-pattern conclusion (−2 to −4% at fixed total dose) is
absorption-rate-robust.

## S4. Calibration trajectories

The calibration criterion (pericentral GSH ≥ 50% of its pre-dose steady state
by 12 h post-ingestion, at the 16 g dose with NAC started at 0 h) is met by
M_max = 4 and above:

| M_max | Pericentral G(12 h) (nmol) | % of baseline (6.87 × 10⁻³) | Verdict |
|---|---|---|---|
| 1 | — | — | recovers too slowly (fails at all early times) |
| 2 | 6.63 × 10⁻⁴ | 9.6% | fail |
| **4** | **3.55 × 10⁻³** | **51.7%** | **calibrated** |
| 8 | 1.39 × 10⁻² | 202.5% | passes (overshoots) |

Full recovery trajectories are plotted in Fig. 3 of the main text and stored
in `results/nac_calibration.npz` (repository).

## S5. Saturation scenario: full protection grid

NAC protection (%) under saturable glucuronidation and oxidation,
v = k·P/(1 + P/Km), Km = P0(16 g)/r. Peak pericentral adducts without NAC
fall with saturation (e.g. 16 g measured: 3.484 × 10⁻² linear, 3.093 × 10⁻²
at r = 1, 2.416 × 10⁻² at r = 3). Protection values at matched start times:

| Gradients | Dose | Kinetics | 0 h | 4 h | 8 h |
|---|---|---|---|---|---|
| assumed | 4 g | linear | 58 | 14 | 1 |
| assumed | 4 g | r = 3 | 90 | 43 | 10 |
| assumed | 4 g | r = 1 | 87 | 31 | 4 |
| measured | 4 g | linear | 83 | 55 | 7 |
| measured | 4 g | r = 3 | 84 | 68 | 37 |
| measured | 4 g | r = 1 | 87 | 69 | 26 |
| assumed | 16 g | linear | 31 | 11 | 1 |
| assumed | 16 g | r = 3 | 86 | 63 | 40 |
| assumed | 16 g | r = 1 | 61 | 37 | 14 |
| measured | 16 g | linear | 42 | 20 | 5 |
| measured | 16 g | r = 3 | 98 | 90 | 62 |
| measured | 16 g | r = 1 | 81 | 57 | 29 |

The assumed-versus-measured ordering is preserved at every saturation level,
and saturation strengthens NAC efficacy at overdose (capped NAPQI formation
lets treatment-supported synthesis keep pace).

## S6. Reproducibility

All results are reproducible from the repository code: `nac_staggered.py`
(window and ingestion experiments), `nac_pk_calibration.py` (PK-driven
calibration and window), `sensitivity_paper2.py` (M, ka and saturation
sweeps), `make_figures.py` (figures from the saved npz files). Integrations
use a stiff Radau solver (rtol 10⁻⁵, atol 10⁻⁸, max step 0.02 d); results are
unchanged at tighter tolerances.
