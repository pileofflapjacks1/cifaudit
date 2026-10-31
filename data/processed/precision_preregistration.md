# Precision study — pre-registration

Committed **before** the sample of findings was drawn or looked at.

## Question

For each structural check, when it fires on a deposited structure, how often is
it right? Every earlier measurement was a firing rate; a firing rate is not a
precision.

## Version under test

Commit `7b2dbbc` with default thresholds. **No threshold or rule will change
until precision has been computed and published.**

## Structures

A fresh COD draw by rsync from `cif/1/02/20`, `2/00/10`, `4/00/50`, `4/20/30`,
`7/00/10`, `7/05/00`, `8/00/10`, `9/00/10` — directories disjoint from all three
earlier samples (16, 500, 700), with zero ID overlap checked against all of them.

## Sampling

1. Audit every file. Keep **structural** findings only. Reporting findings ("this
   field is absent") are checkable mechanically and not worth expert time.
2. Collapse to **one finding per (file, check)** — otherwise a single file with 28
   `H_MULTIPLE_CONNECTIONS` would dominate a check's estimate.
3. Allocation, target 60:
   - 10 each for the three checks with suspicious firing rates on the held-out
     700: `NON_INTEGRAL_CELL_CONTENT`, `FORMULA_SITES_MISMATCH`, `ISOLATED_ATOM`.
   - Then up to 4 each for every other structural check that fires, in descending
     order of firing rate, until the total reaches 60.
   - A check with fewer available findings than its allocation contributes all of
     them.
4. Selection within each check is uniform random with seed **20261031**,
   implemented in `scripts/draw_precision_sample.py`.

## Adjudication

One adjudicator (J. R. Ramirez, Ph.D., crystallographer). Each item is marked
**real**, **false**, or **unsure**, with an optional note. The adjudicator sees
the finding, its evidence, and the relevant CIF lines; blinding to the finding
is not possible, which is a stated limitation. A single adjudicator also means no
inter-rater agreement can be reported.

## Analysis

Per check: precision = real / (real + false), with a Wilson 95% interval.
Unsure counts are reported separately and excluded from the ratio. With 4–10
items per check the intervals will be wide; they will be published as they are.
