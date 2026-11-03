# Expert review — results

Pre-registered in [`expert_review_preregistration.md`](expert_review_preregistration.md)
(`a936f0f`, amendment 1 `8967c77`). The 20 items were drawn at random (seed
20261104) from the 240 findings sampled in precision rounds 1–4.

J. R. Ramirez, a crystallographer, marked each item without seeing the model's
verdicts. The comparison is with the independent model adjudicator's published
verdicts on the same items.

## Result

On the registered question — *should a crystallographer act on this, or look
into it?* — **the expert and the model do not agree on which findings are real.**

| | items | agree | Cohen's κ (bootstrap 95%) | expert real | model real | both real |
|---|---|---|---|---|---|---|
| **pass 2, registered rubric** | 20 | 12 (60%) | **−0.18** (−0.36 to 0.00) | 2 | 6 | **0** |
| pass 2, without AI-assisted marks | 12 | 8 (67%) | −0.14 (−0.33 to 0.00) | 1 | 3 | 0 |
| pass 2, AI-assisted marks only | 8 | 4 (50%) | −0.23 (−0.56 to 0.00) | 1 | 3 | 0 |
| pass 1, "condition present" reading | 20 | 7 (35%) | +0.04 (0.00 to 0.17) | 19 | 6 | 6 |

All agreement in pass 2 is on false positives: both judges call 12 of the 20
items false. Each marks a different, small set as real, and the two sets do not
overlap. Kappa is below zero, and its interval does not reach above zero.

Precision on these 20 items:

- **by the expert: 2/20 · 10% (3–30%)**
- by the model: 6/20 · 30% (15–52%)

The expert is the stricter judge. If this expert is representative, the
model-adjudicated precision figures in this repository (pooled 35%) are
generous.

### The reliable tier does not survive

Seven of the 20 items come from checks that `tiers.py` now labels **reliable**:

- `FORMULA_SITES_MISMATCH` ×3
- `UEQ_RANGE_WIDE` ×2
- `ADP_EXTREME_ANISOTROPY`
- `NO_HYDROGENS_WITH_CARBON`

The model called 5 of the 7 real. **The expert called 0 of 7 real (0–35%).**

The tiers rest entirely on model adjudication, and on this sample an expert does
not share the judgements behind them. One expert and seven items cannot set new
tiers. But they are enough that "reliable" must not be presented as an expert
judgement, and the README now says so.

| tier (current) | items | expert real | model real |
|---|---|---|---|
| reliable | 7 | 0 | 5 |
| prompt | 9 | 1 | 1 |
| unmeasured | 4 | 1 | 0 |

### The rubric drives the result

17 of 20 marks changed between the expert's two passes. Pass 1 answered "is what
the tool says true?", and on that reading 19 of 20 were. Pass 2 answered the
registered question, "is it worth acting on?", and on that reading 2 of 20 were.

Most of cifaudit's findings are **accurate descriptions of things a
crystallographer would not act on**. Whether such a finding counts as correct
depends entirely on which question is asked. Precision figures for a validator
are meaningless without that question stated.

## The two items the expert marked real

| item | check | expert | model |
|---|---|---|---|
| r1-P06 (COD 7001035) | `NON_INTEGRAL_CELL_CONTENT` | real (no note) | false: C 181.5 came from the tool merging a disordered ethyl acetate atom with its inversion image 0.34 Å away |
| r2-Q40 (COD 4310039) | `ATOMS_TOO_CLOSE` | real (AI-assisted, no note) | false: a freely refined C–H of 0.719 Å is a short bond, not two atoms that cannot coexist |

On Q40 the expert's first-pass note gives the reason: a full-occupancy C–H of
0.72 Å is "far below any physically plausible value". The model accepted the
same number and objected only to calling it a clash.

## The six items only the model marked real

| item | check | model's reason (abridged) | expert pass-2 note, where given |
|---|---|---|---|
| r1-P19 | `FORMULA_SITES_MISMATCH` | one methyl H missing from the model (H 31 vs 32) | expected for a 1992 structure, and COD already records it |
| r1-P54 | `UEQ_RANGE_WIDE` | THF solvent with Uiso 0.27–0.40: poorly modelled disorder | a disordered carbon at Uiso 0.40 is ordinary |
| r2-Q23 | `ADP_EXTREME_ANISOTROPY` | ratio 4.67 in a modulated silicate, plus a non-positive-definite O | — |
| r3-R31 | `SHORT_XH` | O–H 0.695 Å, H pointing at Cu at 1.97 Å: misplaced | freely refined X-ray O–H runs short |
| r3-R42 | `UEQ_RANGE_WIDE` | Uiso near 1 Å² on a guest: physically meaningless | — |
| r4-T41 | `NO_HYDROGENS_WITH_CARBON` | an organic dabco-type component with no H in model or formula | — |

## AI assistance

The expert recorded that an AI assistant helped with 8 of the 20 second-pass
marks: P19, P54, Q40, R08, R19, R31, R32 and R47. Leaving those out does not
change the result: kappa is −0.14 and the real sets still do not overlap.
Pass-1 marks had no such flag. Several pass-1 notes contain LaTeX markup and
verbatim PLATON text, so pass 1 should be read as possibly assisted throughout.

## Expert points carried forward

- An occupancy of 0.001 is "for all practical purposes, zero" (R39, pass 1). That
  cuts against the round-3 change making `OCC_ZERO` fire only at exactly 0.
- A full-weight C–H of 0.72 Å is a real error (Q40). Since round 2,
  `ATOMS_TOO_CLOSE` skips an H and its only heavy neighbour and leaves the
  report to `SHORT_XH`, which still fires.

## Limits

- One expert, 20 items, a third of them from round-1 code.
- The rubric was restated between passes, after the expert had seen the
  pass-1 totals but no item-level model verdict.
- Eight second-pass marks were AI-assisted; pass 1 possibly throughout.
- Most second-pass marks have no note, so the expert's reasons are known for
  only a minority of items.

Data: [`expert_review_marks.json`](expert_review_marks.json) holds, per item,
both passes, the AI flag, the pass-2 note, the model's verdict and the check's
current tier. Pass-1 notes are not included.
