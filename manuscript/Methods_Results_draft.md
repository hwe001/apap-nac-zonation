## Methods

### 2.1 Intracellular kinetics

The per-hepatocyte kinetics follow the five-species acetaminophen metabolism
model of Reddyhoff et al. (1), which tracks paracetamol *P*, the sulfation
co-substrate PAPS *S*, the reactive metabolite NAPQI *N*, glutathione *G*, and
covalent protein adducts *C*:

```
dP/dt = −kS·S·P − kG·P − k450·P + kN·N
dS/dt = −kS·S·P + bS − dS·S
dN/dt = k450·P − kN·N − kGSH·N·G − kPSH·N
dG/dt = −kGSH·N·G + bG − dG·G
dC/dt = kPSH·N − k_clear·C
```

The first five terms are those of Reddyhoff et al.; we add a first-order
elimination term to the adduct compartment (`k_clear = 0.42 d⁻¹`, from the
measured human serum adduct elimination rate; James et al. (12)) so that *C*
has the correct turnover. This rate is a **calibrated input**, not an
independent test, and it enters only the adduct turnover, not the treatment
analysis. Parameter values were taken from the original MATLAB source rather
than the published table (whose parameter names appear transposed relative to
their values) and reproduce the paper's stated dosing conversion exactly
(Table 1). The model is human-scaled: simulated doses of 4 g and 16 g,
1.6055 × 10¹¹ hepatocytes, hepatic penetration 0.8, and human-derived
glucuronidation and oxidation rates (1).

**Table 1 — Kinetic parameters (Reddyhoff et al. 2015, as in source code).**

| Parameter | Value | Meaning |
|---|---|---|
| kS | 2.26 × 10¹⁴ · 10⁻¹² cell·mol⁻¹·d⁻¹ | PAPS-dependent sulfation |
| kG | 2.99 d⁻¹ | glucuronidation |
| k450 | 0.315 d⁻¹ | CYP450 oxidation |
| kN | 0.0315 d⁻¹ | NAPQI → APAP reverse |
| bS | 2.65 × 10⁻¹⁴ / 10⁻¹² mol·cell⁻¹·d⁻¹ | PAPS synthesis |
| dS | 2.0 d⁻¹ | PAPS decay |
| kGSH | 1.6 × 10¹⁸ · 10⁻¹² cell·mol⁻¹·d⁻¹ | GSH conjugation of NAPQI |
| kPSH | 110 d⁻¹ | NAPQI–protein binding |
| bG | 1.374 × 10⁻¹⁴ / 10⁻¹² mol·cell⁻¹·d⁻¹ | GSH synthesis |
| dG | 2.0 d⁻¹ | GSH decay |
| k_clear | 0.42 d⁻¹ | adduct elimination (James 2009) |

### 2.2 Zonated extension

The intracellular model is run independently across *N* = 16 hepatocytes
arrayed along a periportal→pericentral sinusoid (hepatocyte 1 periportal, 16
pericentral; 80 state variables total). Each hepatocyte receives the same
initial APAP amount `P0 = dose · f / (MW · N_hep)`, with hepatic penetration
`f = 0.8` and molecular weight `MW = 151 g·mol⁻¹`. Enzyme rates are modulated
by a linear ramp along the sinusoid. For an enzyme with pericentral/periportal
fold *F*, the per-hepatocyte multiplier is

```
m(x) = a + (b − a)·x,   x = (i − 1)/15,   a = 2/(1+F), b = 2F/(1+F),
```

so that the *arithmetic* mean multiplier across the sinusoid is exactly 1 —
total enzyme is conserved and only its location varies. (A
geometric-mean-preserving ramp would not conserve total enzyme.)

### 2.3 Enzyme-gradient representations: assumed versus measured

Two gradient representations are compared throughout.

**Assumed.** The qualitative, Gebhardt-derived gradients carried by prior
zonation models (4, 5, 7): CYP450 5× pericentral; glucuronidation 9×
pericentral; sulfation, GSH conjugation and GSH synthesis 9× periportal.

**Measured.** Abundance folds re-derived from the single-cell spatial
proteomics of human liver (8). That study's published "zonation coefficient"
is the fixed-effect slope of a per-protein linear mixed model fitted to
robust-scaled intensities — a regression slope, *not* a portal-to-central
abundance ratio, and therefore not usable directly as a rate fold. The same
supplement reports each protein's mean intensity in 20 equidistant bins along
the porto-central axis (normalised to sum to 100%; bin 0 central, bin 19
portal). We take the fold as the ratio of the mean of bins 0–4 (central) to
the mean of bins 15–19 (portal); the choice of five bins per end rather than
the terminal bin reduces sensitivity to single-bin noise (single-bin ratios
are within 30% and change no conclusion; full per-protein inputs in
Supplementary Methods). Folds were cross-checked against mRNA
log2(portal/central) ratios from the spatial transcriptomics atlas (9): the
two datasets agree in direction but differ in magnitude, a genuine
mRNA–protein discordance that we carry into the saturation scenario rather
than average away. Because protein abundance is a proxy for catalytic
capacity (no lobule-resolved activity data exist), the measured representation
is a testable scenario, not a ground truth. Where several isoforms contribute
to one parameter, the fold is the isoform-weighted combination using
approximate relative contributions to human APAP oxidation (CYP2E1 ≈ 0.60,
CYP1A2 ≈ 0.25, CYP3A4 ≈ 0.15; Laine et al. (13)).

**Table 2 — Assumed versus measured zonation (pericentral/periportal fold).**

| Parameter | Enzyme | Assumed (prior models) | Measured (protein bin-ratio) |
|---|---|---|---|
| k450 | CYP450 | 5× pericentral | 3.4× (isoform-weighted; raw folds 1.9–6.9×) |
| kG | UGT | 9× pericentral | 1.4× (UGT2B7; 1.4–2.1× across isoforms) |
| kS | SULT | 9× periportal | unzonated (SULT1A1 0.99×) |
| kGSH | GST | 9× periportal | 1.2× pericentral (GSTA2; reversed direction) |
| bG | GSH synthesis | 9× periportal | unzonated (GCLC/GCLM/GSS ~1.0×) |

The two representations differ most consequentially in GSH synthesis: the
assumed gradient starves the pericentral zone of synthetic capacity (fold 0.2×
relative to the mean), whereas the measured enzymes are unzonated, leaving the
pericentral detoxification reserve intact. That difference is what the
treatment analysis below exploits.

### 2.4 NAC term: PK-calibrated GSH-synthesis boost

NAC is represented where it acts: as a cysteine donor that raises hepatic GSH
synthesis. The standard 21-h intravenous protocol drives the input —
150 mg/kg over 1 h, then 50 mg/kg over 4 h, then 100 mg/kg over 16 h. Plasma
NAC follows first-order elimination with the reported intravenous half-life of
~2 h (10); the modelled plasma peak is ≈ 127 mg/kg. Because GSH synthesis
support outlasts plasma NAC itself, the boost is driven by a
cysteine-equivalent pool fed by plasma NAC and eliminated with a half-life of
6 h (a stated assumption, the midpoint of the NAC half-life and hepatocyte
GSH turnover):

```
dA/dt  = I(t) − k_el·A,          A  = plasma NAC (mg/kg),  k_el = ln2 / 2 h
dCys/dt = k_in·A − k_out·Cys,    k_in = k_out = ln2 / 6 h
bG → bG + M_max · B_G · (Cys / Cys_ref)
```

where `I(t)` is the piecewise infusion rate and `Cys_ref` the pool peak for
the standard protocol. The single free magnitude, `M_max`, is **calibrated**,
not scanned arbitrarily. The criterion follows the classic mechanistic finding
of Lauterburg, Corcoran and Mitchell (2): NAC does *not* prevent the initial
GSH depletion caused by a toxic dose — it dramatically accelerates *recovery*
of detoxification capacity through de novo synthesis. We therefore require
that, at the simulated 16 g dose with NAC started at 0 h, pericentral GSH
recovers to ≥ 50% of its own pre-dose steady state by 12 h post-ingestion.
The smallest peak boost meeting the criterion is `M_max = 4×` basal synthesis
(recovery 51.7% of baseline; `M_max = 2×` recovers only 9.6% and fails; Fig.
3). All efficacy scans below use this calibrated, protocol-driven treatment;
a constant-boost scenario (the peak value held flat for 21 h) is reported in
the supplement as an upper bound.

### 2.5 Ingestion term: single and staggered dosing

The instantaneous bolus is generalised with a first-order gut compartment
(`ka = 60 d⁻¹`, absorption half-life ≈ 17 min; swept 17–60 d⁻¹, covering
typical immediate-release absorption): ingested amounts enter the gut and are
absorbed equally by all hepatocytes, with ingestion events applied at
specified times by piecewise integration. Regimens compared at fixed total
dose of 16 g: a single ingestion; 4 × 4 g every 2 h; 8 × 2 g every 1 h.

### 2.6 Saturation scenario and numerical analysis

The Reddyhoff kinetics are linear in substrate. To test whether the treatment
conclusions depend on that simplification — the same one flagged by the
largest current lobule model (7) — the glucuronidation and oxidation terms are
replaced by a saturable form `v = k·P / (1 + P/Km)` with `Km = P0(16 g)/r`,
where `r` interpolates between mild (r = 3) and strong (r = 1) saturation
relative to the 16 g dose; `r → ∞` recovers the linear model. The 82-variable
ODE system (80 cellular states plus the gut and plasma/pool compartments) was
integrated with a stiff Radau solver (`solve_ivp`, rtol 10⁻⁵, atol 10⁻⁸) over
3 days. Metrics: peak pericentral adducts (C_pk), NAC protection
(1 − C_pk(NAC)/C_pk(no NAC) at matched dose, gradient representation, and
start time), and the pericentral GSH minimum. One-at-a-time attribution and
multi-parameter sweeps are reported alongside.

---

## Results

### 3.1 Two gradient representations, one consequence for GSH

The measured human zonation disagrees with the literature-derived gradients in
direction and magnitude (Table 2, Fig. 1). The enzymes that generate and
detoxify NAPQI keep their assumed direction — CYP450 and glucuronidation
pericentral — but are far shallower than assumed (3.4× and 1.4× versus 5× and
9×). Sulfation and, critically, GSH synthesis are unzonated: the assumed
gradients starve the pericentral zone of synthetic capacity (fold 0.2×), the
measured enzymes leave it intact. GSH conjugation is weakly pericentral
rather than strongly periportal. Because pericentral GSH is the rate-limiting
defence in the zone where NAPQI is formed, this difference between
representations is exactly where a treatment term should act.

[[FIG fig1_gradients.png | Assumed versus measured zonation gradients for the five enzyme parameters, as fold-change relative to a uniform baseline (1×). Orange dashed: the Gebhardt-derived gradient carried by prior models; blue solid: the gradient re-derived from human spatial-proteomics bin profiles (isoform-weighted for CYP450). Sulfation and GSH synthesis are unzonated; GSH conjugation is reversed (weakly pericentral); CYP450 and glucuronidation retain direction but are shallower than assumed.]]

### 3.2 The model reproduces the clinical 8-hour rule

With the calibrated protocol-driven NAC term, protection — the reduction in
peak pericentral adducts versus no treatment — collapses as the treatment
start time is delayed (Fig. 2). At the simulated 4 g single dose, protection
falls from 58% (assumed gradients) and 83% (measured) at 0 h to 5% and 27% at
6 h, 1% and 7% at 8 h, and ~0% by 12 h. At the simulated 16 g single dose the
collapse is steeper still (31% and 42% at 0 h; 1% and 5% at 8 h). The model
therefore reproduces, from first principles, the clinical 8-hour rule: the
antidote must be started early because the pericentral GSH pool — once
exhausted by NAPQI — cannot be rescued quickly enough by synthesis alone,
however large the cysteine supply. The width of the useful window is set by
the race between NAPQI-driven GSH consumption and treatment-supported GSH
recovery, not by any parameter of the antidote's own pharmacokinetics.

[[FIG fig2_nac_window.png | Predicted NAC efficacy window: reduction in peak pericentral adducts versus no treatment, as a function of treatment start time after ingestion, at the simulated 4 g and 16 g single doses, under assumed (orange dashed) and measured (blue solid) gradient representations. Treatment is the calibrated, protocol-driven 21-h IV NAC input. Shading marks the late window beyond the clinical 8-hour rule.]]

### 3.3 Measured gradients widen the antidote window

At every start time and both doses, the measured representation predicts
*more* treatment benefit than the assumed one (Fig. 2). The gap is largest
where it matters most clinically: at 4 h, protection under measured gradients
is 55% versus 14% at 4 g, and 20% versus 11% at 16 g; at 8 h, 7% versus 1%
and 5% versus 1%. The mechanism follows directly from §3.1. Under the assumed
gradients the pericentral zone begins with a quarter of the mean GSH-synthesis
capacity and is drained to exhaustion within hours, after which no cysteine
supply can act on a pool that no longer exists. Under the measured,
unzonated-synthesis representation the pericentral reserve is larger from the
start and depletes later, so treatment started in the late window still finds
a pool worth refilling. Early protection is also higher (83% versus 58% at
0 h, 4 g) for the same reason. The zonation assumption a modeller makes
therefore propagates from injury prediction into treatment prediction — with
opposite sign: the representation that predicts *less* injury predicts *more*
treatability.

### 3.4 Staggered ingestion is metabolically neutral at fixed total dose

At a fixed total dose of 16 g, ingestion pattern barely changes the metabolic
burden. Peak pericentral adducts differ by −2.0% (assumed) and −4.0%
(measured) between a single ingestion and 4 × 4 g every 2 h, and by −2.2% and
−4.4% for 8 × 2 g every 1 h; the differences shrink or grow only marginally
across absorption rates of 17–60 d⁻¹ (Table 3). The pericentral GSH minimum
is in fact slightly *higher* under staggering (7.7 × 10⁻⁶ and 8.2 × 10⁻⁶ nmol
versus 4.6 × 10⁻⁶ for the single ingestion, assumed gradients) because each
smaller dose drains GSH less abruptly while synthesis continues between
events. The model does not reproduce the clinical impression that staggered
overdoses are more dangerous through pattern-specific metabolism. What it
does show is that the danger of staggered overdose is a *presentation-time*
phenomenon: staggered ingestions delay the patient's recognition that a toxic
total dose has been accumulated, pushing treatment start into the collapsed
part of the window in Fig. 2 — where, under measured gradients, a genuine
residual benefit remains that the assumed gradients would miss.

**Table 3 — Ingestion-pattern effect at fixed total dose (16 g; peak
pericentral adducts, nmol).**

| Gradients | ka (d⁻¹) | Single | 4 × 4 g q2h | 8 × 2 g q1h |
|---|---|---|---|---|
| assumed | 17 | 3.458 × 10⁻² | 3.391 × 10⁻² | — |
| assumed | 60 | 3.490 × 10⁻² | 3.419 × 10⁻² | 3.415 × 10⁻² |
| measured | 17 | 3.425 × 10⁻² | 3.299 × 10⁻² | — |
| measured | 60 | 3.477 × 10⁻² | 3.338 × 10⁻² | 3.333 × 10⁻² |

[[FIG fig4_staggered.png | Staggered versus single ingestion at fixed total dose (16 g): pericentral GSH (left column) and adducts (right column) over time, under assumed (top row) and measured (bottom row) gradient representations. Staggering slightly raises the GSH minimum and slightly lowers the peak adduct burden; the metabolic footprint of the pattern is small.]]

### 3.5 Robustness: boost magnitude, absorption, and saturable kinetics

The conclusions survive every scenario variation examined. First, the
calibrated boost magnitude is not fragile: under measured gradients the
late-window protection scales with `M_max` (protection at 4 h for 4 g: 11%,
14%, 16% for `M_max` = 1, 2, 4), while under assumed gradients the window
remains closed at 4–6 h for *every* `M_max` — when the pool is gone, no
synthesis boost can act (Fig. S2). Second, the ingestion-pattern result is
absorption-rate-robust (Table 3). Third, saturable kinetics *strengthen*
rather than break the treatment message (Table 4): with `Km = P0(16 g)/r`,
NAC protection at overdose increases steeply with saturation (protection at
0 h, 16 g: 31–42% linear, 61–81% at r = 1, 86–98% at r = 3), because capping
NAPQI formation lets treatment-supported synthesis keep pace with it — and
the assumed-versus-measured ordering persists at every r, widening at 16 g
(8-h protection: 62% versus 40% at r = 3). The linear model is thus the
conservative case for treatment efficacy.

**Table 4 — Saturation robustness (NAC protection, %, at matched start times;
calibrated treatment).**

| Gradients | Dose | Kinetics | 0 h | 4 h | 8 h |
|---|---|---|---|---|---|
| assumed | 4 g | linear | 58% | 14% | 1% |
| assumed | 4 g | r = 3 | 90% | 43% | 10% |
| assumed | 4 g | r = 1 | 87% | 31% | 4% |
| measured | 4 g | linear | 83% | 55% | 7% |
| measured | 4 g | r = 3 | 84% | 68% | 37% |
| measured | 4 g | r = 1 | 87% | 69% | 26% |
| assumed | 16 g | linear | 31% | 11% | 1% |
| assumed | 16 g | r = 3 | 86% | 63% | 40% |
| assumed | 16 g | r = 1 | 61% | 37% | 14% |
| measured | 16 g | linear | 42% | 20% | 5% |
| measured | 16 g | r = 3 | 98% | 90% | 62% |
| measured | 16 g | r = 1 | 81% | 57% | 29% |
