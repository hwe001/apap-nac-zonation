## Discussion

**The antidote window is a modelled quantity — and it inherits the zonation
assumption.** Embedding a PK-calibrated NAC term in a zonation-resolved model
produces three findings that no existing model in this field can generate,
because none of them represents treatment. First, the model reproduces the
clinical 8-hour rule from first principles: protection collapses between 4
and 12 h after ingestion at both simulated doses, because the width of the
useful window is set by the race between NAPQI-driven GSH consumption and
treatment-supported GSH recovery, not by any property of the antidote's own
pharmacokinetics. Second, the representation of hepatic enzyme zonation — to
date an injury-modelling decision — propagates into treatment prediction: the
measured human gradients, being shallower and unzonated for GSH synthesis,
maintain a larger pericentral detoxification reserve and therefore widen the
window (protection at 4 h: 55% versus 14% at 4 g; at 8 h, 7% versus 1%).
Third, at fixed total dose the ingestion pattern itself is metabolically
neutral (−2 to −4%), which relocates the clinical danger of staggered
overdose from metabolism to presentation time — the patient does not know to
seek treatment early, and arrives in the collapsed part of the window.

**Comparison with the state of the art.** The zonated-modelling lineage —
Gebhardt-derived gradients in Means & Ho (5) and Franiatte et al. (6), and
the 5,114-cell lobule coupled to whole-body pharmacokinetics of Ghosh et al.
(7) — has progressively refined the *injury* side: transport, cell states,
dose-dependent zonal sensitivity. None of it models treatment; Ghosh et al.
name NAC only as a clinical implication and identify biomarker-based
validation and Michaelis–Menten kinetics as their own future work, while
simulating only single-bolus overdoses. Conversely, whole-body quantitative
systems toxicology platforms such as DILIsym model NAC explicitly — with
dosing protocols, cysteine delivery and GSH repletion — but carry no zonal
resolution, and so cannot ask where in the lobule the antidote's benefit is
won or lost. The present model sits deliberately between the two: simple
enough to expose the mechanism (a GSH-replenishment race in the zone where
NAPQI is formed), calibrated enough to speak to the clinic (protocol-driven
PK, recovery-anchored magnitude), and first to make the antidote window an
output of a zonation-resolved model. Where our injury-side conclusions touch
those of Ghosh et al., they are consistent — pericentral NAPQI formation
drives the burden; GSH saturation governs overdose — but our treatment-side
results have no counterpart in either literature.

**Clinical reading — with its limits stated.** Three model-supported
statements are, we believe, clinically legible. (i) Early NAC dominates: the
model quantifies how steeply protection decays between 4 and 12 h and
identifies the exhaustion of the pericentral GSH pool as the mechanistic
deadline, supporting the existing emphasis on rapid treatment rather than
dose adjustment. (ii) Late treatment is not zero: under measured gradients a
residual 5–7% adduct reduction persists at 8 h, and saturable kinetics raise
late-window benefit further (up to 62% at 16 h-equivalent scenarios, Table
4) — a quantitative argument against therapeutic nihilism at the 8-hour
boundary, which the empirical literature itself treats as a soft threshold.
(iii) Staggered overdose calls for the same aggressive timeline as acute
overdose: the pattern is metabolically benign at fixed dose, so the observed
worse outcomes of staggered ingestion are consistent with delayed
presentation and treatment, not with a intrinsically more toxic metabolic
profile. We stress what these statements are not: they are not fitted to
clinical outcome data, and the model's adduct compartment is a proxy for
injury rather than necrosis itself; the 8-hour-rule agreement is a
qualitative, order-of-magnitude correspondence of time scales, not a
quantitative validation.

**Why the assumed and measured pictures diverge so sharply for treatment.**
The assumed Gebhardt-derived gradients, being rodent-derived and qualitative,
place GSH synthesis strongly periportal — precisely the wrong arrangement for
the treatment question, since they leave the pericentral zone (where NAPQI is
formed) with a quarter of the mean synthetic capacity and no reserve to
refill. The measured human enzymes are unzonated for GSH synthesis (GCLC,
GCLM, GSS folds ~1.0), consistent with the newer transcriptomic and
proteomic atlases (8, 9) and with classical microdissection data showing that
the classical "periportal GSH" pattern reflects substrate levels set by
consumption, not enzyme distribution. For injury prediction this difference
matters; for treatment prediction it matters more, because the treatment
acts *through* the synthesis pathway. A modeller choosing the assumed
representation is, in effect, simulating a liver whose antidote has less to
work with in the injured zone than a real human liver appears to have.

**Limitations.** Protein abundance is a proxy for catalytic capacity;
post-translational regulation, cofactor availability and isoform kinetics are
not represented, and the isoform weights for CYP450 are approximate
literature contributions (13). The cysteine-equivalent pool half-life (6 h)
is a stated assumption bridging the NAC plasma half-life (2 h) and
hepatocyte GSH turnover; the recovery criterion anchoring the calibrated
boost (50% of baseline by 12 h) is a modelling choice formalising the
Lauterburg et al. (2) phenomenology, not a fitted parameter. The model is a
16-cell sinusoid without transport, oxygen gradients or inflammatory
response; the adduct compartment is a proxy for injury; the adduct
elimination rate is a calibrated input (12); and the simulated 4 g and 16 g
single doses are scenario anchors rather than descriptions of any particular
clinical exposure pattern. The 8-hour-rule correspondence is qualitative —
matching the empirical time scale, not fitted to outcome curves. Donor-level
variability in the measured gradients is not represented; the proteomic folds
are cohort means.

**Outlook.** The treatment question generalises along the same axes as the
injury question. A transport-coupled version would ask whether portal-side
delivery shifts the window; a necrosis-feedback version would let NAC timing
interact with cell death rather than adduct accumulation; donor-level GSH and
CYP variability would turn the window into a population distribution,
identifying which patients lose treatability first. Each extension uses the
same calibrated treatment term built here, and each would move the field
from predicting injury toward predicting treatability — the quantity that
actually decides clinical outcomes.

---

## Code and data availability

The model code, the gradient-derivation inputs, the figures, and an
interactive browser viewer are publicly available at
[github.com/hwe001/apap-zonation](https://github.com/hwe001/apap-zonation)
(and will be mirrored to a dedicated repository for this manuscript). The
submission version will be archived with a DOI via Zenodo and cited here on
acceptance. The viewer displays precomputed simulations and is intended for
inspecting results, not as validation evidence.
