# Paper 2 — outline (independent work, different journal; revised 2026-09-13)

**Framing rule 1 — self-contained:** this manuscript presents the model, the
measured-zonation derivation, and the treatment question from first principles,
citing only primary literature (Gebhardt 1992; Reddyhoff 2015; Means & Ho 2019;
Franiatte 2019; Ghosh 2026; Weiss 2026; Yakubovsky 2026; Lauterburg 1983;
Prescott; Smilkstein 1988; James 2009). No result, method, or parameter is
inherited from the JTB submission by reference — everything the reader needs is
in this paper's Methods.

**Framing rule 2 — different journal:** the target is a pharmacometrics/
toxicology venue (see Target journal below), so the emphasis shifts from
"theoretical advance" (the JTB framing) toward **PK/PD rigour + clinical
translatability**: the calibrated NAC term, the regimen comparison, and the
8-hour-rule mechanism lead; the zonation-sensitivity result is the novel
modelling layer on top.

## Working title

**"When does N-acetylcysteine stop working? Treatment timing, ingestion
pattern, and hepatic zonation in a model of acetaminophen hepatotoxicity"**

Alternates:
- *The antidote window under measured zonation: modelling NAC timing and
  staggered acetaminophen overdose*
- *Treatment, not just injury: N-acetylcysteine timing in a zonated model of
  acetaminophen hepatotoxicity*

## Thesis (one sentence)

Zonation-resolved models of acetaminophen hepatotoxicity have been used to
predict injury but never treatment: embedding a PK-calibrated
N-acetylcysteine term and staggered ingestion into a 16-hepatocyte model whose
enzyme gradients are taken from measured human spatial proteomics shows that
the predicted antidote window reproduces the clinical 8-hour rule, is
substantially widened when literature-assumed gradients are replaced by
measured ones, and that the danger of staggered overdose arises from late
presentation rather than the ingestion pattern itself.

## Positioning (independent of any companion paper)

- The zonation-modelling lineage — Gebhardt-derived gradients in Means & Ho
  (2019), Franiatte et al. (2019), and the 5,114-cell lobule of Ghosh et al.
  (2026) — predicts injury; none of these models contains a treatment term.
  Ghosh et al. simulate only single-bolus overdoses and discuss NAC without
  modelling it.
- Whole-body QST platforms (DILIsym) model NAC but have no zonal resolution.
- No existing model asks whether the *antidote window* — the time after
  ingestion during which NAC still reduces the toxic-adduct burden — depends
  on how liver enzyme zonation is represented.
- The 2026 human spatial-omics datasets (Weiss et al., *Nat Metab*;
  Yakubovsky et al., *Nature*) make enzyme gradients measurable, so the
  zonation assumption can be treated as a testable variable rather than a
  fixed input.

## Contribution statement

1. Extend a zonation-resolved APAP metabolism model (16 hepatocytes, 5
   species, Reddyhoff kinetics) with a **PK-calibrated NAC term** and a
   **staggered-ingestion term** (first-order gut absorption with scheduled
   events).
2. Show the predicted NAC efficacy window **reproduces the clinical 8-hour
   rule** at both simulated doses and under both gradient representations.
3. Show the window is **substantially widened under measured human gradients**
   (protection at 8 h: 7% vs 1% at 4 g; 5% vs 1% at 16 g; protection at 4 h:
   55% vs 14%) — the assumption a modeller makes about enzyme zonation changes
   predicted *treatment* outcomes, not only injury outcomes.
4. Show that at fixed total dose (16 g), **staggered ingestion barely changes
   peak adduct burden** (−2 to −4%, robust across absorption rates) while
   raising the GSH minimum: the clinical danger of staggered overdose is
   attributable to late presentation (result 3), not to the ingestion
   pattern's metabolic footprint.
5. Robustness: the conclusions persist across saturable kinetics
   (Michaelis–Menten scenarios), and the NAC term is PK-calibrated (plasma
   t½ = 2 h; 21-h IV protocol; recovery criterion from Lauterburg 1983), not
   a free scenario parameter.

## Outline

### 1. Introduction
- APAP bioactivation, GSH detoxification, NAC as the standard antidote; the
  8-hour rule is empirical, from cohorts, not from mechanism.
- Zonation models predict injury; none models treatment (Ghosh 2026
  explicitly lacks a treatment term; DILIsym has NAC but no zonal resolution).
- The zonation-gradient problem: every existing model imposes
  literature-derived gradients; 2026 human spatial-omics make them measurable
  — and no one has asked how that choice propagates into *treatment*
  predictions.
- Aim statement (contributions 1–5).

### 2. Methods (fully self-contained)
- 2.1 **Intracellular kinetics**: the five-species Reddyhoff ODE system
  (P/S/N/G/C), parameters from the source code (transposed published table
  noted), adduct elimination k_clear from James 2009 (calibrated input,
  stated).
- 2.2 **Zonated extension**: 16 hepatocytes along a periportal→pericentral
  sinusoid; per-enzyme linear ramps with *arithmetic-mean-preserving*
  multipliers (total enzyme conserved; the geometric form does not conserve
  it).
- 2.3 **Measured zonation**: folds re-derived from the Weiss et al. 20-bin
  intensity profiles (mean of bins 0–4 vs 15–19), *not* from the published
  regression coefficients (which are slopes on robust-scaled data, not fold
  changes); cross-checked against Yakubovsky mRNA; isoform weighting for CYP
  (2E1/1A2/3A4 ≈ 0.60/0.25/0.15, Laine 2009); protein-abundance-as-proxy
  caveat. Gradient table: assumed (Gebhardt-derived, 5×/9×) vs measured
  (CYP 3.4×, UGT 1.4×, SULT/GSH-synthesis unzonated, GST 1.2× pericentral).
- 2.4 **NAC term (PK-calibrated)**: standard 21-h IV protocol
  (150/50/100 mg/kg over 1/4/16 h) → plasma NAC (t½ = 2 h, Prescott) →
  cysteine-equivalent pool (t½ = 6 h, stated assumption) → GSH-synthesis boost
  proportional to the pool, peak M_max calibrated to the Lauterburg recovery
  phenomenology (pericentral GSH recovers to ≥ 50% of baseline by 12 h at the
  16 g dose with NAC from 0 h → calibrated M_max = 4× basal synthesis).
- 2.5 **Ingestion term**: first-order gut compartment (ka ≈ 60 d⁻¹, swept
  17–60), ingestion events via piecewise integration; regimens: single 16 g;
  4 × 4 g q2h; 8 × 2 g q1h.
- 2.6 **Saturation scenario and numerics**: kG, k450 → k·P/(1 + P/Km),
  Km = P0(16 g)/r (r = 1–3); Radau solver; metrics (peak pericentral adducts;
  protection = 1 − Cpk(NAC)/Cpk(no NAC); GSH minimum).

### 3. Results
- 3.1 **The measured gradients are shallower than assumed** (brief, as context
  for treatment: gradient table + figure) — sets up why the two gradient
  representations differ as treatment scenarios.
- 3.2 **The predicted NAC window reproduces the 8-hour rule** (Fig. 1:
  protection vs start time; collapse between 4–12 h at both doses under both
  representations, PK-calibrated treatment).
- 3.3 **Measured gradients widen the window** (8 h: 7% vs 1% at 4 g; 5% vs 1%
  at 16 g; 4 h: 55% vs 14%): uniform GSH synthesis keeps the pericentral
  detoxification reserve higher for longer, which is what late NAC exploits.
- 3.4 **Staggered ingestion is metabolically neutral at fixed total dose**:
  peak adducts −2 to −4% (robust across ka); GSH minimum slightly higher. The
  clinical danger of staggered overdose is late presentation (link to 3.2),
  not the pattern itself.
- 3.5 **Robustness**: saturation (r = 1–3) *strengthens* NAC efficacy at
  overdose (16 g early protection up to 98%) while preserving the
  assumed-vs-measured ordering at every r; window shape robust to the
  cysteine-pool t½ and the recovery-criterion threshold.

### 4. Discussion
- Lead: the zonation assumption a modeller makes propagates from injury
  prediction into treatment prediction — unexamined in the literature (Ghosh
  2026 has no treatment term; DILIsym has NAC but no zonal resolution).
- Mechanism: uniform GSH synthesis (the measured picture) = a larger,
  longer-lived pericentral detoxification reserve; the antidote window is the
  time integral of that reserve.
- Clinical reading: the model supports early NAC above all, and quantifies
  residual benefit at 8 h; staggered-overdose danger = late presentation,
  arguing for aggressive treatment timelines rather than pattern-specific
  metabolism.
- Limitations: protein abundance as capacity proxy; pool t½ assumption;
  recovery-criterion threshold; adducts as proxy for injury; no transport, no
  inflammation, single cell type; the 8-h match is qualitative (no clinical
  outcome data fitted); rat/mouse-vs-human translation of the calibration
  anchor.
- Future: transport-coupled window; necrosis feedback interacting with NAC
  timing; population GSH distributions (virtual cohorts); other antidotes.

### 5. Conclusion
The antidote window is a modelled, measurable quantity: it inherits — with
opposite sign — the same zonation-assumption sensitivity that injury
predictions show, and measured human gradients give the more optimistic
treatment picture.

## Figures & tables

- Fig. 1 — assumed vs measured gradients for the five enzymes (context figure;
  to be generated standalone in this repo).
- Fig. 2 — NAC protection vs start time, 4 g and 16 g, assumed vs measured
  (the "8-hour rule" figure; PK-calibrated). *Data in README; figure to
  regenerate with the PK-driven values.*
- Fig. 3 — M_max calibration (GSH recovery trajectories vs the criterion).
- Fig. 4 — staggered vs single ingestion: pericentral G(t) and C(t)
  (2 × 2). *Exists: `figures/staggered_patterns.png`.*
- Table 1 — kinetic parameters.
- Table 2 — assumed vs measured gradient folds (with derivation provenance).
- Table 3 — PK-driven window (full grid; also constant-M comparison in
  supplement).
- Table 4 — ka sensitivity; Table 5 — saturation robustness.

## Target journal (different from paper 1)

Ranking for a PK-calibrated treatment-modelling paper with a clinical hook:

1. **CPT: Pharmacometrics & Systems Pharmacology** — the natural home: this is
   exactly their material (antidote PK/PD, regimen comparison, QSP with
   calibrated treatment terms; the DILIsym population publishes there). The
   zonation layer is the novel modelling contribution; the calibrated NAC term
   and 8-hour-rule mechanism are the translational contribution.
2. **Toxicological Sciences** — mechanism + toxicology audience; strong fit for
   the GSH-recovery mechanism and the adduct-burden endpoint; regular APAP
   coverage.
3. **Toxicology and Applied Pharmacology** — the historical APAP/NAC venue
   (Corcoran/Mitchell lineage); conservative but thematically perfect.
4. **Frontiers in Pharmacology** — open access, and where Ghosh 2026 sits;
   enables direct engagement with the lobule-modelling group.

Journal-specific items to check at submission: abstract structure/word limit,
whether Highlights are used (Elsevier only), graph-abstract requirements
(ToxSci), and whether the NAC PK calibration needs supplementary PK tables.

## Dependencies / risks of the independent framing

- None on the JTB submission: every method and number used here is derived
  inside this paper (fold derivation, gradient schemes, model). If paper 1 is
  accepted first, it can be cross-cited as a companion *without* dependence;
  if it is not, this paper stands alone.
- The two papers share the core model code; the repos are separate
  (`apap-zonation` for paper 1; this folder to become paper 2's repo).
