# Paper 2 — outline (drafted 2026-09-13, before the sensitivity sweeps finished)

## Working title

**"When does N-acetylcysteine stop working? Treatment timing, ingestion
pattern, and hepatic zonation in a model of acetaminophen hepatotoxicity"**

Alternates:
- *The antidote window under measured zonation: modelling NAC timing and
  staggered acetaminophen overdose*
- *Treatment, not just injury: N-acetylcysteine timing in a zonated model of
  acetaminophen hepatotoxicity*

## Thesis (one sentence)

Zonated acetaminophen models have been used to predict injury but never
treatment: embedding N-acetylcysteine — as a time-limited glutathione-synthesis
boost — and staggered ingestion into a zonation-resolved model shows that the
predicted antidote window reproduces the clinical 8-hour rule, is widened by
measured (rather than literature-assumed) enzyme gradients, and that the
danger of staggered overdose arises from late presentation rather than the
ingestion pattern itself.

## Positioning

- Ghosh et al. 2026 (the current frontier lobule model) simulates only
  single-bolus overdoses and has **no treatment term at all** — NAC appears
  only as a Discussion remark.
- Our paper 1 (submitted, JTB) supplies the measured gradients, the
  arithmetic-mean-conserving zonation machinery, and the dose-dependent
  stability view; paper 2 reuses that infrastructure and shifts the question
  from *injury prediction* to *treatment prediction*.
- Clinical hook: the "8-hour rule" for NAC is empirical; a mechanistic model
  that reproduces it — and predicts how it shifts with liver enzyme
  organization — is a genuinely new use of zonation modelling.

## Contribution statement

1. Extend a zonation-resolved APAP model with a **treatment term** (NAC as a
   time-limited GSH-synthesis boost) and a **staggered-ingestion term**
   (first-order gut absorption with scheduled events).
2. Show the predicted NAC efficacy window **reproduces the clinical 8-hour
   rule** at both doses and under both gradient schemes.
3. Show the window is **widened under measured human gradients** (16% vs 2%
   residual protection at 8 h for 4 g) — i.e., the assumption a modeller makes
   about enzyme zonation changes predicted treatment outcomes, not only
   injury outcomes.
4. Show that at fixed total dose (16 g), **staggered ingestion barely changes
   peak adduct burden** (−2 to −4%) while raising the GSH minimum: the clinical
   danger of staggered overdose is attributable to late presentation (result
   3), not to the ingestion pattern's metabolic footprint.
5. Robustness: the window and the assumed-vs-measured gap persist across NAC
   boost magnitudes (M = 1–4), absorption rates (ka = 17–60 d⁻¹), and a
   saturable-kinetics scenario (Michaelis–Menten with Km = P0(16 g)/r,
   r = 1–3).

## Outline

### 1. Introduction
- NAC works by repleting GSH; the 8-hour rule is empirical, derived from
  cohorts, not from mechanism.
- Zonation models predict injury; none models treatment (Ghosh 2026 explicitly
  lacks a treatment term).
- Paper 1's measured-gradient machinery provides the framework; the new
  question: does zonation matter for *treatment* the way it matters for
  *injury*?
- Aim statement (contributions 1–5).

### 2. Methods
- 2.1 Core model recap (brief — cite paper 1): 16 hepatocytes, 5 species,
  assumed vs measured gradients, arithmetic-mean conservation, k_clear.
- 2.2 **NAC term**: bG → bG + M·B_G for a 21-h window (standard IV protocol
  duration), start time scanned 0–12 h; M is a scenario parameter, swept
  1–4.
- 2.3 **Ingestion term**: first-order gut compartment (ka), ingestion events
  via piecewise integration; regimens: single 16 g; 4 × 4 g q2h; 8 × 2 g q1h.
- 2.4 **Saturation scenario**: kG, k450 → v = k·P/(1 + P/Km), Km = P0(16 g)/r.
- 2.5 Metrics: peak pericentral adducts; NAC protection = 1 − Cpk(NAC)/Cpk(no
  NAC); pericentral GSH minimum.

### 3. Results
- 3.1 The model reproduces the **8-hour rule** (Fig. 1: protection vs start
  time; collapse between 4–12 h at both doses, both schemes).
- 3.2 **Measured gradients widen the window** (late-window protection 16% vs
  2% at 4 g; 9% vs 3% at 16 g; early protection also higher). The gap is a
  treatment-relevant consequence of the gradient assumption.
- 3.3 **Staggered ingestion is metabolically neutral at fixed total dose**:
  peak adducts −2 to −4%; GSH minimum slightly higher. The clinical danger of
  staggered overdose is late presentation (link to 3.1/3.2), not the pattern
  itself.
- 3.4 **Robustness**: window shape and scheme gap persist across M = 1–4
  (Fig. 2), ka = 17–60 d⁻¹ (Table), and saturable kinetics r = 1–3 (Table);
  late-window protection is dose- and saturation-dependent in magnitude but
  the assumed-vs-measured ordering is preserved.

### 4. Discussion
- Lead: the zonation assumption a modeller makes propagates from injury
  prediction into treatment prediction — a consequence unexamined in the
  literature (Ghosh 2026 has no treatment term; DILIsym has NAC but no
  zonation).
- Mechanism: uniform GSH synthesis (the measured picture) keeps the pericentral
  detoxification reserve higher for longer, which is what late NAC exploits.
- Clinical reading: the model supports early NAC above all, and quantifies how
  much residual benefit remains at 8 h; staggered-overdose danger = late
  presentation, arguing for aggressive treatment timelines rather than
  pattern-specific metabolism.
- Limitations: M not calibrated to mg/kg NAC PK (scenario parameter, swept);
  ka approximate; adducts as proxy; no transport, no inflammation; single
  cell type; the 8-h "match" is qualitative (no clinical outcome data fitted).
- Future: transport-coupled window; necrosis feedback making NAC timing
  interact with cell death; population-level GSH distributions (virtual
  cohorts).

### 5. Conclusion
Zonation models can say something about *treatment*, not only about injury:
the antidote window is a modelled, measurable quantity, and it inherits —
in the opposite direction — the same gradient-assumption sensitivity that
paper 1 established for injury.

## Figures & tables

- Fig. 1 — NAC protection vs start time, 4 g and 16 g, assumed vs measured
  (the "8-hour rule" figure). *Exists: `figures/nac_window.png`.*
- Fig. 2 — M-sensitivity: protection vs start time for M = 1, 2, 4
  (4 panels: scheme × dose). *Exists: `figures/sweep_M.png`.*
- Fig. 3 — staggered vs single ingestion: pericentral G(t) and C(t)
  (2 × 2: scheme × species). *Exists: `figures/staggered_patterns.png`.*
- Table 1 — ka sensitivity (single vs staggered peak C at ka = 17/34/60).
- Table 2 — saturation robustness (protection at 0/4/8 h, linear vs r = 1, 3).

## Target journal

Journal of Theoretical Biology again (same framing family: theory producing a
biologically intelligible insight — here about treatment, not injury); fallback
CPT:PSP or Toxicology Letters for the clinical-adjacent reading.
