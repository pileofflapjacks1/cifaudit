# Second expert review — pre-registration

Committed before the second rater is recruited or has seen any item.

## Question

The first expert review ([`expert_review_results.md`](expert_review_results.md))
found no overlap between the findings the expert and the model judge called
worth acting on: κ = −0.18 on 20 items. That is uninterpretable without a
baseline. **Do two crystallographers agree with each other more than either
agrees with the model?**

- If the two experts agree closely and both diverge from the model, the model
  judge is miscalibrated for this task.
- If the experts also disagree with each other, whether a finding is worth
  acting on is not a stable judgement in this setting. Then no single judge,
  model or human, can set a precision figure for it.

## Items

The same 20 items as the first review (seed 20261104,
`scripts/draw_expert_review.py`), unchanged.

## Rater

One crystallographer who did not build the tool. The author (J. R. Ramirez) is
the first rater; the second is recruited by the author. Their experience is
recorded on the review page (years refining small-molecule structures, role) and
reported without their name unless they ask to be named.

## Procedure

The second rater marks each item on a separate private page:

- **The question is the first review's second-pass question, word for word:**
  should a crystallographer act on this, or look into it?
- **The choices are the same:** Real (worth acting on), False positive (misread,
  or normal and expected), Unsure. Each mark has an optional note and records
  whether an AI assistant helped.
- **The page shows the same material:** the finding as cifaudit reported it,
  the CIF excerpt, and a link to the COD entry. Olex2 or any crystallographic
  software may be used, but not cifaudit.
- **Nothing anchoring is shown:** neither the model's verdicts, nor the first
  rater's marks, nor any tier or precision figure.

**A known risk:** both sets of verdicts are public in this repository. The brief
asks the rater not to look at the repository, the README or the write-up until
they have finished. The rater is asked afterwards whether they saw any of it,
and the answer is reported.

## Analysis

Unchanged from the first review, with expert–expert agreement added. The second
pass of the first review is the first rater's comparable set.

| Comparison | Statistic |
|---|---|
| Rater 1 vs rater 2 | raw agreement, Cohen's κ, 4,000-resample bootstrap 95% interval |
| Rater 2 vs model | the same |
| Rater 1 vs model | as published: κ = −0.18 (−0.36 to 0.00) |

Also reported:

- each rater's precision on the 20 items
- the tier breakdown for rater 2
- a sensitivity row without AI-assisted marks
- every three-way disagreement, with all available reasons

**Decision rule, stated in advance.** "Experts agree with each other more than
with the model" is claimed only if the expert–expert κ interval lies entirely
above both expert–model κ intervals. The reverse claim, that experts agree no
better with each other than with the model, needs the expert–expert interval to
overlap at least one expert–model interval and its point estimate to fall below
0.40. Otherwise the result is reported as undetermined at n = 20. With 20 items
an undetermined result is likely, and is a finding in itself.

Nothing in the code, the tiers or the write-up's conclusions changes on the
basis of these marks until this analysis is published.
