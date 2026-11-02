# Expert review of the precision findings — pre-registration

Committed before the expert has seen any item.

## Question

Every precision figure so far rests on one independent model adjudicator. **How
far does an expert crystallographer agree with it?** The answer bears on every
tier and every precision figure in this repository.

## Items

There are 20 items, drawn by `scripts/draw_expert_review.py` with seed
**20261104**. The draw is a simple random sample from all 240 findings sampled in
precision rounds 1–4. It is not stratified by the model's verdict, because
agreement statistics depend on how often each verdict occurs, and
stratification would distort them. The script was run once before this file
was committed, to check that it works. The seed was not changed afterwards.

## Procedure

The expert (J. R. Ramirez) marks each item on a private review page:

- **real**: the condition is present and worth a crystallographer's attention
- **false**: misread, or normal and expected for this structure
- **unsure**

This is the rubric the model adjudicator used. For each item the page shows:

- the finding as cifaudit reported it
- the CIF excerpt the adjudicator saw
- a link to the full entry at the Crystallography Open Database

**The model's verdict is not shown.** Neither is the tool's tier, nor any
precision figure, so neither can anchor the expert's judgement.

## Analysis

Over the items the expert marks:

- raw agreement and Cohen's kappa between expert and model, over real / false /
  unsure, using `scripts/score_precision.py --human`
- precision by the expert's marks, beside the model's on the same items
- every disagreement listed with both reasons

With 20 items, kappa will have a wide interval and will be reported as such.
If fewer than 20 are marked, the analysis covers those marked, and the count is
stated. Nothing in the code or the tiers changes on the basis of these marks
until the analysis is published.
