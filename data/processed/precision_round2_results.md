# Precision results, round 2 — model adjudication

Pre-registered in [`precision_round2_preregistration.md`](precision_round2_preregistration.md)
(`40fdaec`, amendment 1 `2ac59df`). Version measured: `40fdaec` (no change under
`src/` since `354e9fd`), default thresholds.

**This is model adjudication, not expert review.** A fresh agent with none of the
conversation that built cifaudit judged each finding against the full deposited
CIF. It recomputed geometry and counts itself with gemmi and numpy and never
imported or ran cifaudit. The human review has not been done, so no inter-rater
agreement can be reported.

There are 60 structural findings, drawn from 800 COD structures. None overlaps the
1,616 IDs used in earlier sets. Seed 20261101.

## Result

**Overall precision is 36% (95% CI 25–48%).** Round 1 measured 38% (27–51%) on
the code before the fixes. The intervals overlap almost entirely: **round 2 shows
no measurable change in overall precision.**

| check | round 1 (`7b2dbbc`) | round 2 (`40fdaec`) |
|---|---|---|
| `FORMULA_SITES_MISMATCH` | 6/10 · 60% (31–83%) | 7/10 · 70% (40–89%) |
| `ISOLATED_ATOM` | 1/10 · 10% (2–40%) | 3/10 · 30% (11–60%) |
| `NON_INTEGRAL_CELL_CONTENT` | 0/10 · 0% (0–28%) | not drawn (fired in no file, as predicted) |
| `ADP_EXTREME_ANISOTROPY` | 2/4 · 50% | 4/4 · 100% (51–100%) |
| `UEQ_RANGE_WIDE` | 4/4 · 100% | 3/4 · 75% |
| `SHORT_HH_CONTACT` | 2/4 · 50% | 1/4 · 25% |
| `ATOMS_TOO_CLOSE` | 1/4 · 25% | 0/4 · 0% (0–49%) |
| `SHORT_XH` | 2/4 · 50% | 0/4 · 0% (0–49%) |
| `H_MULTIPLE_CONNECTIONS` | 0/2 · 0% | 0/4 · 0% (0–49%) |
| `DISORDER_SUM_NOT_ONE` | — | 0/4 · 0% (0–49%) |
| `SCATTERING_TYPE_LABEL_MISMATCH` | — | 1/4 · 25% |
| `DENSITY_BELOW_ONE` | — | 2/4 · 50% |
| `SHORT_NONBONDED_CONTACT` | — | 0/3, 1 unsure |
| `NO_HYDROGENS_WITH_CARBON` | 3/4 · 75% | — |
| `DUPLICATE_LABELS` | 2/4 · 50% | — |
| **all** | **23/60 · 38% (27–51%)** | **21/59 · 36% (25–48%)**, 1 unsure |

No per-check difference is distinguishable at these sample sizes. The two rounds
also drew different check mixes. Round 2 has none of the fixed
`NON_INTEGRAL_CELL_CONTENT` items, which round 1 scored 0/10. It adds four checks
that round 1 never drew, and two of those (`DISORDER_SUM_NOT_ONE` and
`SHORT_NONBONDED_CONTACT`) score 0. The
overall figures therefore do not measure the same mixture.

Full confidence breakdown: false 29 high and 9 medium; real 9 high, 10 medium and
2 low; unsure 1 low.

### Lenient view (labelled; not comparable to round 1)

The adjudicator also recorded whether each finding is accurately worded and graded,
even where the condition is normal (amendment 1). Counting those as correct gives
**30/59 = 51% (38–63%)**. Nine items change:

- five ordinary lattice waters reported as an isolated O with the solvent message
  (Q11–Q15)
- three freely refined X-ray X–H shortenings (Q29, Q30, Q32)
- one nominal-formula mismatch (Q01)

Round 1 was not scored this way, so this view sits beside the comparable figure
and never replaces it.

## Why findings were false

1. **Disorder across a symmetry element** (8 items: Q33, Q35, Q37, Q38, Q41, Q57,
   Q58, Q60). This covers SHELX negative PART, and partial occupancy on a mirror,
   2-fold axis or inversion centre. An atom and its own symmetry image are mutually
   exclusive alternatives, but cifaudit treats them as two real atoms, giving
   clashes, short contacts and extra connections. This one cause produced most of
   the false geometric findings.
2. **Disorder labels that don't match the disorder** (Q25–Q28, Q31, Q34, Q36).
   Sometimes one disorder is split across several assembly letters, or several
   disorders share one letter. The real group sums are 1, but cifaudit sums by the
   labels as written.
3. **Normal features reported as problems.** These include hydrogen-bonded lattice
   water as an isolated O (Q11–Q15), low density in SQUEEZE'd frameworks (Q55,
   Q56), tunnel water with large Ueq (Q47), and an Fe-only magnetic phase (Q20).
4. **A bonded X–H treated as a clash, or a contact just under its threshold.**
   Q39 and Q40 are O–H 0.741 Å and C–H 0.719 Å, reported as "cannot both be real".
   Q29, Q30 and Q32 are freely refined X-ray X–H. In Q42 and Q44, H…H is 0.01–0.02 Å
   below 2.0 Å.
5. **Label conventions.** "Wat" was read as element "Wa" (Q49, Q51, Q52). In Q01 and
   Q08 the declared formula gives a nominal or analytical composition.

## Bugs reported, and checked

Before publishing, each was independently recomputed from the CIF outside the
adjudicator.

| item | defect | verification |
|---|---|---|
| Q16 (COD 4310078) | S1 reported isolated although it is 2.461 Å from four Ta1. | **Confirmed.** Ta has no entry in the covalent-radius table, which also lacks Hf and the lanthanides after Ce. The check records "connectivity involving Ta" as not tested, but then still reports S1 as isolated, instead of abstaining for an atom whose neighbours it cannot judge. |
| Q09 (COD 7041079) | H count 25 per formula; the CIF gives 28. | **Confirmed.** gemmi's unit-cell expansion merges H7B, H7C and H9A with their own mirror images 0.29–0.38 Å away, which are half-occupied PART −1 disorder. C1, 0.52 Å from its image, is kept. The exact merge tolerance is unverified. |
| Q57, Q58, Q60 | `SHORT_NONBONDED_CONTACT` reports C…H pairs. | **Confirmed as an inconsistency.** `audit.py` documents `nonbonded_min` as a heavy-atom limit, but `checks/contacts.py` applies it to every non-H…H pair. |
| systematic | Atoms disordered across a symmetry element are not treated as exclusive with their own images. | Consistent with cause 1 above. |
| Q39, Q40 | `ATOMS_TOO_CLOSE` does not exclude a bonded X–H pair. | Consistent with the messages in the findings. |
| Q49, Q51, Q52 | `SCATTERING_TYPE_LABEL_MISMATCH` treats a label prefix as an element without checking it is one. | Consistent with the messages. |

Per the pre-registration, none of these is fixed in the version measured here.

## Judgement calls

- Q02, Q03, Q06, Q07 and Q10 are formula mismatches caused by unlocated water, NH4
  or OH hydrogens. They are marked real at medium confidence: the model and the
  formula genuinely disagree, which is what PLAT041 alerts on. A stricter "normal
  and expected" reading would make them false and lower precision to 16/59.
- Q59 is unsure. It is an F…F contact of 2.14 Å between partly occupied atoms
  with no disorder labels, consistent with undeclared correlated disorder.
- Q43 is real at low confidence.

## Limits

- One model adjudicator and no human rater.
- 2–10 items per check.
- Different check mixes between rounds.
- The confirmed bugs bias round 2 downwards, relative to what the checks are meant
  to do.

Files:

- [`precision_round2_sample_items.json`](precision_round2_sample_items.json): the
  findings as drawn, without CIF text.
- [`precision_round2_model_verdicts.json`](precision_round2_model_verdicts.json)
