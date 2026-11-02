# Precision study, round 4 — pre-registration

This was committed **before** any structure was fetched or audited, and before
any finding was drawn. The candidate directories were checked by rsync listing
only, which shows file names but no file contents.

## Question

Round 3 measured 26% (95% CI 16–38%) by model adjudication on `33f30cb`. Rounds
1 and 2 measured 38% and 36%. After round 3 was published, `6f1058f` fixed the
six defects it confirmed:

- the X–H reach
- `OCC_ZERO` at the tolerance boundary
- negative PART and full-weight members in the disorder sum
- the copy-merge tolerance for negative PART
- mismatched assembly letters

On the round-2 and round-3 samples, which were used to design the fixes, 37 of
81 false findings no longer fire and all 36 real ones still do. That is not a
precision figure. **Does precision move on structures the fixes never saw?**

## Version under test

The commit containing this file, with default thresholds. No threshold or rule
changes until round-4 precision is computed and published.

## Structures

All CIFs are fetched by rsync from eight COD directories, about 800 files:

- `cif/2/20/10`
- `4/13/10`, `4/34/20`
- `7/03/50`, `7/12/10`
- `8/10/60`
- `9/00/80`, `9/01/40`

Any file whose COD ID appears in an earlier set is excluded, and the overlap is
reported. That covers 3,216 IDs: the 16, 500, 700, 400, 800 and 800 sets.

## Sampling

Sampling is identical to rounds 1–3:

- the same `scripts/draw_precision_sample.py`
- structural findings only, one per (file, check)
- 10 each for `NON_INTEGRAL_CELL_CONTENT`, `FORMULA_SITES_MISMATCH` and
  `ISOLATED_ATOM`
- then up to 4 for each other firing check, in descending firing rate, until 60
- a check with fewer findings than its allocation contributes all of them

New seed: **20261103**.

## Adjudication

A fresh model adjudicator is used. It has none of the conversation that built
cifaudit, judges each item against the full CIF, recomputes what it needs, and is
barred from running or importing cifaudit.

Its brief is round 3's, with one change: the descriptions of the checks altered
by `6f1058f` are updated to what the code now does (X–H reach, `OCC_ZERO`, and
how disorder closure and alternatives are decided). The verdict rule is round 1's,
unchanged, with both expert clarifications from J. R. Ramirez:

- an isolated O is usually modelled solvent
- SQUEEZE'd solvent is recorded by authors and normally included in the declared
  formula

"Real" means the condition is present and worth a crystallographer's attention,
whatever the severity. `accurately_framed` is recorded for the separately
labelled lenient view.

## Analysis

Per check and overall, precision = real / (real + false), with a Wilson 95%
interval. Unsure verdicts are reported and excluded. Rounds 1–4 are shown side
by side. With 2–10 items per check the intervals will be wide, and overlapping
intervals will not be described as a change. The pooled per-check table is
updated and labelled exploratory, as before.

Every bug the adjudicator reports is recomputed independently before it is
published. None is fixed in the version measured here.

## Amendment 1 — 2026-11-01, at adjudicator launch, before any verdict

The fetch and draw went as registered. 800 files were fetched, with 0 overlap
against the 3,216 earlier IDs. Seed 20261103 drew 60 items. Auditing the 800
files took longer than in earlier rounds, but nothing failed.

Three checks were drawn whose descriptions round 3's brief did not contain:
`DENSITY_BELOW_ONE`, `NO_HYDROGENS_WITH_CARBON` and `SYMMETRY_POSSIBLY_MISSED`.
Each was described in one line taken from what the code does. The adjudicator
was also told spglib is available, so it can recompute the derived space group.
Nothing else in the brief differs from round 3's beyond the updates registered
above.
