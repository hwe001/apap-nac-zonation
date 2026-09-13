# Paper 2 — NAC timing and staggered overdose under measured hepatic zonation

Working folder for the follow-up to the submitted JTB paper
("How much does assumed enzyme zonation matter for predicted acetaminophen
hepatotoxicity? A stability analysis against measured human profiles" —
companion repo: `github.com/hwe001/apap-zonation`, submission version v1.1.1).

## Motivation

The 2026 Ghosh et al. lobule model (Frontiers in Pharmacology) simulates only
single-bolus overdoses and does not model N-acetylcysteine (NAC) at all — it
appears only as a Discussion remark. This paper asks two clinically loaded
questions the existing zonation literature cannot answer:

1. **Does the predicted NAC efficacy window depend on the zonation gradients?**
   NAC is modeled as a temporary boost to hepatic GSH synthesis
   (bG → bG + M·B_G during a 21-h window, the standard IV protocol duration),
   started at variable times after ingestion. Metric: reduction in peak
   pericentral protein adducts vs no NAC.
2. **Does ingestion pattern (single vs staggered overdose) change the predicted
   injury?** The bolus is replaced by a first-order gut compartment
   (ka = 60 d⁻¹ ≈ 2.5 h⁻¹) with ingestion events at specified times
   (4 × 4 g q2h; 8 × 2 g q1h; total 16 g).

## First results (code/nac_staggered.py, 2026-09-13)

**NAC efficacy window** (M = 4 × basal synthesis; protection = reduction in
peak pericentral adducts):

| Gradients | Dose | 0 h | 8 h | 12 h |
|---|---|---|---|---|
| assumed | 4 g | 85% | 2% | 0% |
| assumed | 16 g | 45% | 3% | 0% |
| measured | 4 g | 89% | **16%** | 0% |
| measured | 16 g | 56% | **9%** | 0% |

- The model reproduces the clinical **8-hour rule**: protection collapses
  between 4 and 12 h post-ingestion at both doses and under both gradient
  schemes.
- **Measured gradients widen the late window**: at 8 h, NAC still achieves
  16% (4 g) / 9% (16 g) reduction under measured gradients vs 2% / 3% under
  assumed — because uniform GSH synthesis leaves pericentral GSH higher
  initially and depleting later, so late NAC still finds a treatable state.
  Early (0 h) protection is also higher under measured gradients.

**Staggered ingestion** (16 g total):

| Gradients | Pattern | peak C (nmol) | min pericentral GSH (nmol) |
|---|---|---|---|
| assumed | single | 3.490e-2 | 4.64e-6 |
| assumed | 4 × 4 g q2h | 3.419e-2 | 7.74e-6 |
| assumed | 8 × 2 g q1h | 3.415e-2 | 8.20e-6 |
| measured | single | 3.477e-2 | 4.95e-6 |
| measured | 4 × 4 g q2h | 3.338e-2 | 6.93e-6 |
| measured | 8 × 2 g q1h | 3.333e-2 | 7.24e-6 |

- At fixed total dose (16 g), staggered ingestion barely changes peak adducts
  (−2 to −4%) and slightly *raises* the GSH minimum — the model does **not**
  reproduce the clinical impression that staggered overdoses are more toxic
  through zonal metabolism itself. The clinically observed danger of staggered
  overdose is therefore attributable to **late presentation** (which the NAC
  window result captures directly), not to the ingestion pattern's metabolic
  footprint. This is an honest, useful negative result and a paper-2 talking
  point, not a bug.

## Caveats (carry into any manuscript)

- The NAC boost magnitude M is a scenario parameter (multiples of basal
  synthesis), not calibrated to mg/kg protocols; before submission, calibrate
  against published NAC pharmacokinetics (cysteine/GSH repletion data).
- ka = 60 d⁻¹ is an approximate human absorption rate; sensitivity-check it.
- Metric is peak pericentral adducts — a proxy, not necrosis (same limitation
  as paper 1).
- No inflammatory response, no transport, single hepatocyte type.

## Sensitivity and robustness (code/sensitivity_paper2.py, 2026-09-13)

**1. NAC boost magnitude M ∈ {1, 2, 4}** (`figures/sweep_M.png`): the
qualitative picture is M-robust — under *assumed* gradients the window is
closed by 4–6 h at every M (GSH is already gone before any boost can act);
under *measured* gradients the late window scales with M (4 h protection
11% → 14% → 16% at 4 g). The assumed-vs-measured ordering is preserved for
every M, so the uncalibrated boost size affects magnitude, not conclusion.

**2. Absorption rate ka ∈ {17, 34, 60} d⁻¹** (t½ ≈ 1 h, 30 min, 17 min): the
staggered-vs-single difference in peak adducts stays in a narrow band
(−1.9% to −4.0%) at every ka and both schemes — "staggered is metabolically
neutral at fixed total dose" is absorption-rate-robust.

**3. Saturation robustness** (kG, k450 → k·P/(1 + P/Km), Km = P0(16 g)/r):
saturation **strengthens** the clinical message rather than breaking it —

| 8-h protection (4 g / 16 g) | linear | r = 3 | r = 1 |
|---|---|---|---|
| assumed gradients | 2% / 3% | 10% / 40% | 4% / 14% |
| measured gradients | 16% / 9% | 37% / 62% | 26% / 29% |

- NAC becomes *more* effective under saturable kinetics at overdose (capped
  NAPQI formation lets GSH repletion keep up): 16 g early protection rises
  from 45–56% (linear) to 61–98% (r = 1–3).
- The **assumed-vs-measured ordering persists at every r** — the gradient
  assumption matters for treatment predictions under linear *and* saturable
  kinetics, and the gap at 16 g actually widens with saturation.
- Base adduct burden falls with saturation (3.49e-2 → 2.42e-2 at 16 g
  measured), consistent with capped oxidation.

## Roadmap (agreed 2026-09-13)

1. **This paper**: NAC window + staggered dosing (fast, distinctive — Ghosh
   models neither). Add NAC PK calibration, ka sensitivity, MM-kinetics
   robustness check.
2. **Transport-coupled stability sequel**: measured gradients + uncertainty in
   the Franiatte 2019 advection-diffusion framework (code:
   `E:\Google Drive\Student_Projects\internship\2018\Sylvain_Franiatte\`);
   also a real test of the glucuronidation anomaly.
3. **Michaelis–Menten saturation** × measured zonation (Ghosh's announced next
   step — our edge is combining it with measured gradients + stability test).
4. **Necrosis → clearance feedback + adduct release**: converts paper 1's
   calibrated k_clear into a genuine held-out serum-adduct validation
   (James 2009).
5. **Zonally-resolved virtual cohorts** (long game; needs donor-level enzyme
   distributions).

## Layout

```
code/
  intracellular_apap_model.py   single-hepatocyte kinetics (Reddyhoff 2015)
  zonated_apap_model.py         16-hepatocyte zonated model (paper 1 core)
  nac_staggered.py              NAC window + staggered-ingestion experiments
figures/                        nac_window.png, staggered_patterns.png
results/                        nac_window.npz
manuscript/                     (draft sections to come)
```
