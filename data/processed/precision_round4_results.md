# Precision results, round 4 — model adjudication

Pre-registered in [`precision_round4_preregistration.md`](precision_round4_preregistration.md)
(`50f1789`, amendment 1 `ef4d316`). Version measured: `50f1789` (no change under
`src/` since `6f1058f`), default thresholds.

**This is model adjudication, not expert review.** A fresh agent with none of the
conversation that built cifaudit judged each finding against the full deposited
CIF. It recomputed geometry, counts and space groups itself with gemmi, numpy and
spglib, and never imported or ran cifaudit.

There are 60 structural findings, drawn from 800 COD structures. None overlaps the
3,216 IDs used in earlier sets. Seed 20261103.

## Result

**Overall precision is 40% (95% CI 29–53%).** The intervals of all four rounds
overlap, so by the pre-registered rule there is **still no measurable change**:

| round | version | precision | 95% CI |
|---|---|---|---|
| 1 | `7b2dbbc` | 23/60 · 38% | 27–51% |
| 2 | `40fdaec` | 21/59 · 36% | 25–48% |
| 3 | `33f30cb` | 15/58 · 26% | 16–38% |
| 4 | `50f1789` | 24/60 · 40% | 29–53% |

Per check:

| check | round 4 | rounds 1–3 |
|---|---|---|
| `FORMULA_SITES_MISMATCH` | 4/10 · 40% (17–69%) | 6/10, 7/10, 7/10 |
| `ISOLATED_ATOM` | 0/10 · 0% (0–28%) | 1/10, 3/10, 0/10 |
| `ADP_EXTREME_ANISOTROPY` | 4/4 | 2/4, 4/4, 4/4 |
| `UEQ_RANGE_WIDE` | 4/4 | 4/4, 3/4, 2/3 |
| `NO_HYDROGENS_WITH_CARBON` | 3/4 | 3/4, —, — |
| `SYMMETRY_POSSIBLY_MISSED` | 3/4 | first drawn |
| `SHORT_HH_CONTACT` | 2/4 | 2/4, 1/4, 0/4 |
| `SHORT_NONBONDED_CONTACT` | 2/4 | —, 0/3, 0/4 |
| `H_MULTIPLE_CONNECTIONS` | 2/4 | 0/2, 0/4, 0/4 |
| `DENSITY_BELOW_ONE` | 0/4 | —, 2/4, — |
| `DISORDER_SUM_NOT_ONE` | 0/4 | —, 0/4, 0/4 |
| `SHORT_XH` | 0/4 | 2/4, 0/4, 1/3 |

The check mix changes every round. This round drew three checks that score well
(`NO_HYDROGENS_WITH_CARBON`, `SYMMETRY_POSSIBLY_MISSED` and a full set of
`UEQ_RANGE_WIDE`), and none of the occupancy and H-reach checks that round 3 drew
and that scored 0/4. That alone can account for the movement from 26% to 40%.

Confidence breakdown: false 18 high, 17 medium, 1 low; real 8 high, 11 medium,
5 low; no unsure.

### Lenient view (labelled; not comparable to round 1)

Counting accurately worded findings as correct gives **46/60 = 77% (65–86%)**.
The 22 items that change are:

- nine isolated solvent, water or counter-ion atoms (T11, T13–T20)
- three SQUEEZE'd or porous low-density structures (T25–T27)
- five formula mismatches that are SQUEEZE or AMCSD conventions (T01, T03,
  T04, T06, T09)
- four freely refined short C–H (T45–T48)
- one symmetry finding (T60)

## Across all four rounds (exploratory, not pre-registered)

The table pools 237 judged findings from four code versions. It describes the
checks' behaviour across those versions, not one release.

| check | pooled | 95% CI |
|---|---|---|
| `ADP_EXTREME_ANISOTROPY` | 14/16 · 88% | 64–97% |
| `UEQ_RANGE_WIDE` | 13/15 · 87% | 62–96% |
| `NO_HYDROGENS_WITH_CARBON` | 6/8 · 75% | 41–93% |
| `SYMMETRY_POSSIBLY_MISSED` | 3/4 · 75% | 30–95% |
| `FORMULA_SITES_MISMATCH` | 24/40 · 60% | 45–74% |
| `SHORT_HH_CONTACT` | 5/16 · 31% | 14–56% |
| `DENSITY_BELOW_ONE` | 2/8 · 25% | 7–59% |
| `SHORT_XH` | 3/15 · 20% | 7–45% |
| `SHORT_NONBONDED_CONTACT` | 2/11 · 18% | 5–48% |
| `ATOMS_TOO_CLOSE` | 2/12 · 17% | 5–45% |
| `H_MULTIPLE_CONNECTIONS` | 2/14 · 14% | 4–40% |
| `ISOLATED_ATOM` | 4/40 · 10% | 4–23% |
| `DISORDER_SUM_NOT_ONE` | 0/12 · 0% | 0–24% |
| **all** | **83/237 · 35%** | **29–41%** |

The split from round 3 holds with more data:

- **Usually right when they fire:** anisotropy, the Ueq range, missing H with
  carbon, missed symmetry and, less reliably, the formula count.
- **Mostly wrong:** the distance, connectivity and disorder-closure checks.

`ISOLATED_ATOM` is wrong 36 times in 40, almost always because it reports
ordinary solvent or counter-ions, which it describes correctly and which the
round-1 rule scores false. `DISORDER_SUM_NOT_ONE` has never been right in 12
draws: real CIFs use assembly letters too loosely for the sum to mean anything.

## Why findings were false

1. **Ordinary counter-ions, solvent and water** (T11–T20).
2. **Expected effects of SQUEEZE or porosity:** formula mismatches (T01, T03,
   T10) and every density item (T25–T28).
3. **Disorder bookkeeping, not wrong occupancies.** T29 and T32 put
   alternatives under different assembly letters, so they are not summed
   together. T30 has two independent disorders (0.60/0.40 and 0.42/0.58) sharing
   one letter. T31 is disorder across a 2-fold axis.
4. **Symmetry-related alternatives treated as coexisting:** T33 (methanol across
   an inversion centre), T50 (Im-3m orientational disorder), T55 (correlated
   disorder across an inversion centre), and T51 (two 0.25-occupancy waters
   2.07 Å apart with no disorder group declared).
5. **Conventions:**
   - freely refined X-ray C–H of 0.80–0.85 Å, within about 3 esd (T45–T48)
   - AMCSD leaving H out of `_chemical_formula_sum` while declaring attached H
     (T04, T06, T09)

## Bugs reported, and checked

Before publishing, each was independently recomputed from the CIF outside the
adjudicator.

| item | defect | verification |
|---|---|---|
| T54 (COD 7035035) | A Cp-ring centroid placeholder, CT1, typed H with occupancy 0, was treated as a hydrogen with 5 heavy neighbours. | **Confirmed.** Its row is `CT1 H … 0.000 Uiso 0 1 d`. It has no `dum` calc flag, so it is not dropped as a dummy, and zero-occupancy sites take part in the geometry. |
| T05 (COD 4342091) | A formula mismatch was downgraded to info with "SQUEEZE explains it", although the model has *more* H than the formula (26 vs 24). | **Confirmed.** The SQUEEZE convention is applied whenever SQUEEZE is recorded, whichever way the counts differ. Removed solvent can only make the formula larger than the model. |
| T10 (COD 7035070) | SQUEEZE not recognised. | **Confirmed.** The `_platon_squeeze_*` names occur only inside the embedded `_shelx_fab_file` text (lines 483–15831), not as data names. |
| T33 (COD 4342081) | A methanol disordered across an inversion centre was not treated as alternatives. | **Confirmed, a known limit.** C15 (occupancy 0.5) sits 2.12 Å from its own image, so the rule does catch C15. The clash is H15B with H15C′, an H bonded to that image rather than the image itself, and the rule follows one bond only. The same case appeared in round 2 (COD 4102011). |
| T30 (COD 4131088) | Two independent disorders under one assembly letter are summed as one. | Consistent with the labels: 1.349 is what the labels give. Splitting assemblies by occupancy pattern would be a new heuristic, not a fix. |
| T50, T55 | Symmetry images in the same disorder group are not recognised as alternatives. | Not recomputed separately. Recorded as reported. |

Per the pre-registration, none of these is fixed in the version measured here.

## Notable real findings

These are the adjudicator's findings. They were not independently recomputed.

- **T49:** in (NH₄)₂SiF₆ the deposited N…F distance is 1.983 Å. Moving N from
  z = 0.25 to z = 0.75 gives 2.91 Å, so the N coordinate is almost certainly
  wrong.
- **T58:** tetradymite is declared R-3, but spglib finds R-3m, the accepted
  space group.

## Limits

- One model adjudicator and no human rater.
- 4–10 items per check per round.
- The check mix differs between rounds.
- The brief was extended for the three newly drawn checks (amendment 1).

Files:

- [`precision_round4_sample_items.json`](precision_round4_sample_items.json): the
  findings as drawn, without CIF text.
- [`precision_round4_model_verdicts.json`](precision_round4_model_verdicts.json)
