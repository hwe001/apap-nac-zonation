## Introduction

N-acetylcysteine (NAC) is the standard antidote in acetaminophen
(paracetamol) poisoning. Acetaminophen is cleared mainly by sulfation and
glucuronidation, but a minor oxidative pathway — mediated by CYP2E1, CYP1A2
and CYP3A4 — forms the reactive metabolite *N*-acetyl-*p*-benzoquinone imine
(NAPQI), which is normally detoxified by glutathione (GSH) conjugation (1).
When overdose exhausts hepatic GSH, NAPQI binds covalently to hepatocyte
proteins and centrilobular necrosis follows. NAC works upstream of that
endpoint: as a cysteine donor it fuels de novo GSH synthesis, restoring the
detoxification capacity of hepatocytes whose GSH pool NAPQI is draining (2).
Clinically, its benefit is strongly time-dependent: treatment started within
~8 h of ingestion is highly effective, whereas benefit declines sharply
thereafter — the "8-hour rule" (3). That rule, however, is an empirical
summary of clinical cohorts. No mechanistic model currently reproduces it,
predicts how it shifts, or explains what sets its width.

Zonation is the obvious candidate mechanism. The enzymes of acetaminophen
metabolism are zonated along the periportal–pericentral axis of the hepatic
lobule (4), and a sequence of zonation-resolved models — a finite-element
sinusoid (5), a transport-coupled extension (6), and most recently a
5,114-hepatocyte lobule coupled to whole-body pharmacokinetics (7) — has shown
that the spatial arrangement of enzyme activity shapes the predicted toxic
adduct burden. Two gaps are shared by this entire lineage. First, the enzyme
gradients are *imposed*: qualitative, literature-derived profiles ultimately
traceable to a rodent-derived review (4), with the largest current model
stating explicitly that no omics data were integrated (7). Second — and more
consequential for the present work — **none of these models contains a
treatment term**. NAC appears in the most recent of them only as a clinical
remark in the Discussion (7); whole-body quantitative systems toxicology
platforms represent NAC but carry no zonal resolution. A model of zonation
that cannot represent the antidote cannot ask the treatment question at all.

Two further developments make the question answerable now. First, human
spatial omics has made the enzyme gradients measurable: single-cell spatial
proteomics of human liver reports per-cell abundance along the porto-central
axis for the relevant enzymes (8), cross-checked by a multi-species spatial
transcriptomics atlas (9). Abundance-derived gradient folds differ from the
assumed ones — most GSH-related activity is unzonated rather than strongly
periportal — and protein abundance is a *proxy* for catalytic capacity, so
the representation is a testable scenario rather than a ground truth. Second,
the pharmacokinetics of the standard intravenous NAC protocol are reasonably
characterised (10, 11), so a treatment term can be calibrated rather than
invented.

Here we embed both pieces in a 16-hepatocyte model of the hepatic sinusoid
(80 state variables) built on the Reddyhoff intracellular kinetics (1): a
PK-calibrated NAC term (a time-limited, protocol-driven boost to GSH
synthesis) and a staggered-ingestion term (a first-order gut compartment
receiving scheduled doses). We ask three questions. Does the model reproduce
the clinical 8-hour rule? Does the choice between assumed and measured
enzyme gradients — an injury-modelling decision to date — propagate into
*treatment* predictions? And does the ingestion pattern itself (single bolus
versus staggered overdose at fixed total dose) alter the metabolic burden, or
is the clinical danger of staggered overdose fully explained by late
presentation? The answers — yes; yes, in the opposite direction to injury;
and the latter — carry a practical implication: the zonation assumption a
modeller makes changes when the antidote is predicted to stop working.
