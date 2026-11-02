# Precision study, round 2 — pre-registration

Committed **before** any structure was audited or any finding drawn.

## Question

Round 1 measured 38% (95% CI 27–51%) by model adjudication on `7b2dbbc`. Since
then four bugs were fixed and seven conventions adopted, all designed against the
round-1 sample, which therefore cannot measure them. **Did precision move on
structures the fixes never saw?**

## Version under test

The commit containing this file, with default thresholds. No threshold or rule
changes until round-2 precision is computed and published.

## Structures

All CIFs fetched by rsync from `cif/1/03/40`, `2/02/30`, `2/03/50`, `4/10/20`,
`4/31/00`, `7/02/40`, `7/04/10`, `8/10/30`, `9/00/30`, `9/01/20`. Any file whose
COD ID appears in an earlier sample (1,616 IDs across the 16, 500, 700 and 400
sets) is excluded; the overlap is reported.

## Sampling

Identical to round 1, so the two rounds are comparable by construction: the same
`scripts/draw_precision_sample.py`, structural findings only, one per (file,
check), 10 each for `NON_INTEGRAL_CELL_CONTENT`, `FORMULA_SITES_MISMATCH` and
`ISOLATED_ATOM`, then up to 4 per other firing check in descending firing rate
until 60. A check with fewer findings than its allocation contributes all of
them. New seed **20261101**.

Expected consequence, stated in advance: `NON_INTEGRAL_CELL_CONTENT` now fires
only when there is no formula to compare against, so it may contribute few or no
items, and the check mix will differ from round 1. The overall figure is
therefore compared alongside per-check figures, never alone.

## Adjudication

A fresh model adjudicator with none of the conversation that built cifaudit,
judging each item against the full CIF, recomputing what it needs, and barred
from running cifaudit. It receives the same rubric as round 1, including the two
expert clarifications from J. R. Ramirez that round 1 applied: an isolated O is
usually modelled solvent, and SQUEEZE'd solvent is recorded by authors and
normally included in the declared formula. Its verdicts are kept separate from
the human review page.

A human review page with the same 60 items is offered to the expert. If he marks
any, Cohen's kappa and raw agreement are reported over those items.

## Analysis

Per check and overall: precision = real / (real + false), Wilson 95% interval;
unsure reported and excluded. Round 1 and round 2 are shown side by side. With
2–10 items per check the intervals will be wide and overlapping intervals will
not be described as a change.
