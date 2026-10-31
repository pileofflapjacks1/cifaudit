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

## Amendment 1 — 2026-10-31, after the draw, before any adjudication

A model adjudicator is added as an **independent second rater**. It is a fresh
agent with none of the conversation that built cifaudit, so that it carries no
investment in the checks. It judges each item against the full deposited CIF,
recomputing geometry or counts itself where useful, and is instructed **not** to
run cifaudit. Its verdicts are stored separately from the human review page, so
the human adjudicator is not shown them.

Reporting:
- Per-check precision is reported separately for each adjudicator, labelled as
  **human (expert)** or **model**. Model verdicts are never presented as expert
  adjudication.
- Where both have judged the same items, Cohen's kappa and raw agreement are
  reported over those items.
- If only the model has judged, the README says so in the headline.

The limitation this does not remove: a model judging a model-built tool's output
is weaker evidence than an expert's, whatever its stated confidence.

## Amendment 2 — 2026-10-31, rubric clarification during model adjudication

From the expert adjudicator (J. R. Ramirez), passed to the model adjudicator
while it was running and applied to every `ISOLATED_ATOM` item, including any
already judged:

> An unbound single O atom is commonly solvent — typically water whose H atoms
> were not located — or an O placed to account for electron density that may
> belong to a different molecule. Crystallographers routinely do this.

So an isolated O that fits modelled solvent (water-style label, partial
occupancy, a void, O···O/N contacts near 2.6–3.1 Å) is judged a normal modelling
choice and the flag `false`; an isolated O with no reasonable solvent
interpretation may still be `real`. Other elements are unaffected.

This is a clarification of how to judge, not a verdict on any item. The check
itself is **not** changed during the study; a revision is queued for after the
results are published.
