# Precision study, round 3 — pre-registration

This was committed **before** any structure was fetched or audited, and before
any finding was drawn. The candidate directories were checked by rsync listing
only, which shows file names but no file contents.

## Question

Round 2 measured 36% (95% CI 25–48%) by model adjudication on `40fdaec`. After
that figure was published, `f38423e` fixed:

- three confirmed bugs: missing covalent radii, gemmi merging near-coincident
  copies, and X…H documentation
- their shared cause: disorder across a symmetry element
- two smaller defects

All of these were designed against the round-2 sample, so that sample cannot
measure them. **Did precision move on structures the fixes never saw?**

## Version under test

The commit containing this file, with default thresholds. No threshold or rule
changes until round-3 precision is computed and published.

## Structures

All CIFs are fetched by rsync from eight COD directories, about 800 files:

- `cif/2/01/40`, `2/10/20`
- `4/12/30`, `4/33/00`
- `7/01/30`, `7/11/00`
- `8/10/50`
- `9/00/70`

Any file whose COD ID appears in an earlier set is excluded, and the overlap is
reported. That covers 2,416 IDs: the 16, 500, 700, 400 and 800 sets.

## Sampling

Sampling is identical to rounds 1 and 2:

- the same `scripts/draw_precision_sample.py`
- structural findings only, one per (file, check)
- 10 each for `NON_INTEGRAL_CELL_CONTENT`, `FORMULA_SITES_MISMATCH` and
  `ISOLATED_ATOM`
- then up to 4 for each other firing check, in descending firing rate, until 60
- a check with fewer findings than its allocation contributes all of them

New seed: **20261102**.

The check mix will again differ from earlier rounds. The overall figure is
therefore compared alongside per-check figures, never alone.

## Adjudication

A fresh model adjudicator is used. It has none of the conversation that built
cifaudit, judges each item against the full CIF, recomputes what it needs, and is
barred from running or importing cifaudit.

It gets the round-1 rubric unchanged, including the two expert clarifications from
J. R. Ramirez:

- an isolated O is usually modelled solvent
- SQUEEZE'd solvent is recorded by authors and normally included in the declared
  formula

"Real" means the condition is present and worth a crystallographer's attention,
whatever the severity. As round-2 amendment 1 introduced, the adjudicator also
records `accurately_framed` per item, for a separate lenient view that is
labelled as such and never replaces the comparable figure.

## Analysis

Per check and overall, precision = real / (real + false), with a Wilson 95%
interval. Unsure verdicts are reported and excluded. Rounds 1, 2 and 3 are shown
side by side. With 2–10 items per check the intervals will be wide, and
overlapping intervals will not be described as a change.

Every bug the adjudicator reports is recomputed independently before it is
published. None is fixed in the version measured here.

## Amendment 1 — 2026-11-01, at adjudicator launch, before any verdict

The fetch and draw went as registered. 800 files were fetched, with 0 overlap
against the 2,416 earlier IDs. Seed 20261102 drew 60 items.

The verdict rule given to the adjudicator is round 1's, unchanged, and the
round-2 correction is written in from the start. The *check descriptions* differ
from round 2, in three ways:

1. **Updated to the version under test.** `ATOMS_TOO_CLOSE` treats an H and its
   only heavy neighbour as a bond. `SHORT_NONBONDED_CONTACT` covers X…H pairs,
   as the code always applied it. `OCC_ZERO` and `ISOLATED_H` are described,
   because they were drawn for the first time.
2. **A note on intended disorder handling.** One paragraph says how the tool is
   meant to treat disorder alternatives, including across a symmetry element.
   The adjudicator is told to judge against the actual structure, not that intent.
3. **A recomputation note.** gemmi's `get_all_unit_cell_sites()` merges copies
   within about 0.4–0.5 Å, so the adjudicator should apply the symmetry operators
   itself near symmetry elements. The round-2 adjudicator found this
   independently.

Items 2 and 3 tell the adjudicator more about the tool than round 2's brief did.
That could make it more likely to accept the tool's disorder handling. This is
recorded so the comparison with round 2 can be read with it in mind.
