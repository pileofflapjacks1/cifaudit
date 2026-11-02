# Precision results, round 3 — model adjudication

Pre-registered in [`precision_round3_preregistration.md`](precision_round3_preregistration.md)
(`33f30cb`, amendment 1 `84756ae`). Version measured: `33f30cb` (no change under
`src/` since `f38423e`), default thresholds.

**This is model adjudication, not expert review.** A fresh agent with none of the
conversation that built cifaudit judged each finding against the full deposited
CIF. It recomputed geometry and counts itself with gemmi and numpy, applying the
symmetry operators itself, and never imported or ran cifaudit.

There are 60 structural findings, drawn from 800 COD structures. None overlaps the
2,416 IDs used in earlier sets. Seed 20261102.

## Result

**Overall precision is 26% (95% CI 16–38%).** Round 2 was 36% (25–48%) and
round 1 38% (27–51%). The intervals overlap, so by the pre-registered rule this
is **no measurable change**. It is certainly not the improvement the round-2
fixes were meant to bring. Under the same rule a fall can't be claimed either.

| check | round 1 (`7b2dbbc`) | round 2 (`40fdaec`) | round 3 (`33f30cb`) |
|---|---|---|---|
| `FORMULA_SITES_MISMATCH` | 6/10 | 7/10 | 7/10 · 70% (40–89%) |
| `ISOLATED_ATOM` | 1/10 | 3/10 | 0/10 · 0% (0–28%) |
| `ADP_EXTREME_ANISOTROPY` | 2/4 | 4/4 | 4/4 · 100% (51–100%) |
| `UEQ_RANGE_WIDE` | 4/4 | 3/4 | 2/3, 1 unsure |
| `SHORT_XH` | 2/4 | 0/4 | 1/3, 1 unsure |
| `ATOMS_TOO_CLOSE` | 1/4 | 0/4 | 1/4 |
| `SHORT_HH_CONTACT` | 2/4 | 1/4 | 0/4 |
| `H_MULTIPLE_CONNECTIONS` | 0/2 | 0/4 | 0/4 |
| `DISORDER_SUM_NOT_ONE` | — | 0/4 | 0/4 |
| `SHORT_NONBONDED_CONTACT` | — | 0/3, 1 unsure | 0/4 |
| `OCC_ZERO` | — | — | 0/4 |
| `ISOLATED_H` | — | — | 0/4 |
| other checks | 5/18 | 3/8 | — |
| **all** | **23/60 · 38% (27–51%)** | **21/59 · 36% (25–48%)** | **15/58 · 26% (16–38%)**, 2 unsure |

The check mix differs again. Round 3 drew `OCC_ZERO` and `ISOLATED_H` for the
first time, and both score 0/4. Its ten isolated atoms were all ordinary solvent
or counter-ions: water, Br⁻, Cl⁻ and NH₄⁺. Every one was described correctly,
and each counts as false under the round-1 rule.

The adjudicator's brief described more of the tool's intended disorder handling
than round 2's did (amendment 1). If that biased it, the bias would favour the
tool, and the figure is lower.

Confidence breakdown: false 28 high, 14 medium, 1 low; real 3 high, 6 medium,
6 low; unsure 2 low.

### Lenient view (labelled; not comparable to round 1)

Counting accurately worded findings as correct gives **28/58 = 48% (36–61%)**.
The 13 items that change are:

- ten isolated solvent or counter-ion atoms (R11–R20)
- unlocated water H (R03)
- a SQUEEZE'd formula (R05)
- a refined zero occupancy (R40)

## Across all three rounds (exploratory, not pre-registered)

Pooling the three rounds gives each check 7–30 items. The rounds measured
different code versions, so these figures describe the checks' behaviour across
versions, not one release.

| check | pooled | 95% CI |
|---|---|---|
| `ADP_EXTREME_ANISOTROPY` | 10/12 · 83% | 55–95% |
| `UEQ_RANGE_WIDE` | 9/11 · 82% | 52–95% |
| `FORMULA_SITES_MISMATCH` | 20/30 · 67% | 49–81% |
| `SHORT_XH` | 3/11 · 27% | 10–57% |
| `SHORT_HH_CONTACT` | 3/12 · 25% | 9–53% |
| `ATOMS_TOO_CLOSE` | 2/12 · 17% | 5–45% |
| `ISOLATED_ATOM` | 4/30 · 13% | 5–30% |
| `H_MULTIPLE_CONNECTIONS` | 0/10 · 0% | 0–28% |
| `DISORDER_SUM_NOT_ONE` | 0/8 · 0% | 0–32% |
| `SHORT_NONBONDED_CONTACT` | 0/7 · 0% | 0–35% |
| **all** | **59/177 · 33%** | **27–41%** |

The pattern is stable. Three checks are usually right when they fire: anisotropy,
the Ueq range, and the formula count. The geometric and connectivity checks are
mostly wrong. They depend on reconstructing disorder, which CIFs describe
inconsistently, and on fixed distance thresholds, which do not hold across the
periodic table. Fixing them one false positive at a time has not moved overall
precision in two rounds.

## Why findings were false

1. **H bonded beyond the flat 1.35 Å X–H cutoff.** P–H bonds of 1.35–1.38 Å and
   bridging metal hydrides (Lu, Pd, Sc; 1.69–2.13 Å) were read as unbonded H.
   That produced false `ISOLATED_H`, `SHORT_HH_CONTACT` and
   `SHORT_NONBONDED_CONTACT` findings (R34, R53–R55, R57–R60). R56 is a B–H…Ag
   agostic contact.
2. **Ordinary solvent and counter-ions** (R11–R20).
3. **Disorder labels the tool trusts.** The two groups of one disorder appear
   under different assembly letters, or one partner has no letter, so the tool
   treats them as coexisting: R28, R32, R46–R49, R51, R52.
4. **Negative PART in the disorder-sum check.** A methyl at PART −1 on a mirror
   sums to 0.5, because its other half is its own image (R25, R27).
5. **Expected or trivial conditions:**
   - an occupancy of 0.001 reported as 0 (R37–R39), and a refined zero occupancy
     in a site-mixing model (R40)
   - a SQUEEZE'd formula (R05) and unlocated water H (R03)
   - riding H…H contacts 0.005–0.08 Å under 2.0 Å (R33, R35, R36)
   - a disordered hexane with high U (R41)
   - an aromatic H judged against the sp³ C–H target (R30)

## Bugs reported, and checked

Before publishing, each was independently recomputed from the CIF outside the
adjudicator.

| item | defect | verification |
|---|---|---|
| R37–R39 (COD 9007098 and others) | Occupancy 0.001 reported as "occupancy 0.0". | **Confirmed, different cause.** Parsing is correct (0.001). `OCC_ZERO` fires at \|occ\| ≤ `occ_tol` = 0.001, so the boundary value counts as zero, and the message then states 0.0. |
| R34, R53 (COD 2014097, 4123099) | One X–H cutoff of 1.35 Å for every element. | **Confirmed.** P–H is 1.351 Å, just outside. The Lu–H hydrides are 1.971–2.143 Å, with no heavy atom within 1.6 Å of the H. |
| R25 (COD 7110006) | `DISORDER_SUM_NOT_ONE` ignores negative PART. | **Confirmed.** The only group in assembly A is PART −1 at 0.50. The other half is its mirror image. |
| R28 (COD 4123051) | A group's occupancy is taken from its largest member. | **Confirmed.** H74A/B are tagged group 1 with occupancy 1, so the sum is 1 + 0.25. The real C73/C11I pair is 0.75 + 0.25. |
| R08 (COD 7013013) | The merge tolerance for coincident copies is too tight. | **Confirmed.** C22 and C23 are 0.039 Å from their images and merge; H27A is 0.068 Å away and does not. All three are PART −1 at 0.5 on a mirror, so they should be treated alike. The 0.05 Å value set in `f38423e` is too tight; gemmi's ~0.4 Å was too loose. |
| systematic | Exclusivity depends on assembly letters being correct. | Consistent with cause 3 above. |

Per the pre-registration, none of these is fixed in the version measured here.

## Judgement calls

- R01 and R10 are counted as real: H missing on an organic group, or on a
  template, in a formula check. R03 (unlocated water H) is counted as false.
- R02 and R09 (NbO₂F, a shared O/F site with no occupancy column) are counted as
  real at low confidence. They could reasonably be called false.
- R29 and R44 are unsure.

## Limits

- One model adjudicator and no human rater.
- 3–10 items per check per round.
- The check mix differs between rounds.
- The brief differed from round 2's (amendment 1).

Files:

- [`precision_round3_sample_items.json`](precision_round3_sample_items.json): the
  findings as drawn, without CIF text.
- [`precision_round3_model_verdicts.json`](precision_round3_model_verdicts.json)
