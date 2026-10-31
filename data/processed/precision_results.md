# Precision results — model adjudication

Pre-registered in [`precision_preregistration.md`](precision_preregistration.md)
(`19ee1e2`, with amendments `b2f8352` and `36c2424`). Version measured:
`7b2dbbc`, default thresholds.

**This is model adjudication, not expert review.** A fresh agent with none of the
conversation that built cifaudit judged each finding against the full deposited
CIF, recomputing geometry and counts itself, without running cifaudit. The human
expert review page exists but has not been completed, so no inter-rater agreement
can be reported yet.

60 structural findings from 400 COD structures never used for tuning.

| check | real | false | unsure | precision | 95% CI |
|---|---|---|---|---|---|
| `NON_INTEGRAL_CELL_CONTENT` | 0 | 10 | 0 | 0% | 0–28% |
| `FORMULA_SITES_MISMATCH` | 6 | 4 | 0 | 60% | 31–83% |
| `ISOLATED_ATOM` | 1 | 9 | 0 | 10% | 2–40% |
| `SHORT_XH` | 2 | 2 | 0 | 50% | 15–85% |
| `ADP_EXTREME_ANISOTROPY` | 2 | 2 | 0 | 50% | 15–85% |
| `ATOMS_TOO_CLOSE` | 1 | 3 | 0 | 25% | 5–70% |
| `SHORT_HH_CONTACT` | 2 | 2 | 0 | 50% | 15–85% |
| `NO_HYDROGENS_WITH_CARBON` | 3 | 1 | 0 | 75% | 30–95% |
| `UEQ_RANGE_WIDE` | 4 | 0 | 0 | 100% | 51–100% |
| `DUPLICATE_LABELS` | 2 | 2 | 0 | 50% | 15–85% |
| `H_MULTIPLE_CONNECTIONS` | 0 | 2 | 0 | 0% | 0–66% |
| **all** | 23 | 37 | 0 | 38% | 27–51% |

**Overall precision is 38% (27–51%).** Fewer than half of the structural findings
on deposited structures are correct.

## Why the false ones were false

Spot-checked by recomputing independently; every claim checked held up.

**Three bugs in cifaudit, confirmed:**

1. **Unwrapped symmetry coordinates.** gemmi's unit-cell expansion leaves
   fractional coordinates outside [0, 1) (e.g. Zn at z = −0.75 and 1.25), and the
   distance search only tries lattice shifts of ±1 from the raw difference. Pairs
   whose copies sit more than ~1.5 cells apart lose their true nearest image. In
   COD 2... (P21) two Br atoms 2.42 Å from Zn came out "isolated". This affects
   every distance-based check and can hide real contacts as well as invent
   false ones.
2. **The anisotropy ratio uses raw CIF U_ij.** The earlier argument that no
   transformation is needed holds for the NPD *sign* test but not for eigenvalue
   *ratios*: in non-orthogonal cells sqrt(U3/U1) must be computed in Cartesian
   axes. P35 is 4.10 raw but 2.76 in Cartesian axes; P38 is 4.30 raw but 3.65.
   The 4.0 threshold was calibrated on the same wrong values and needs
   recalibrating. `ADP_NPD` is unaffected.
3. **Primes are stripped from labels.** Quote-stripping meant for CIF quoting
   turns `C1'` and `C1''` into `C1`, so distinct atoms read as duplicates.

**Checks that flag normal chemistry:**

- `NON_INTEGRAL_CELL_CONTENT` (0/10): fires on solid solutions whose CIF declares
  the non-integral formula, e.g. `Al3.96 Ca0.98 Na1.97 Si4.04` with Z = 4.
  Where content and formula disagree, `FORMULA_SITES_MISMATCH` already catches it.
- `ISOLATED_ATOM` (1/10): an isolated O is usually modelled solvent (expert
  rubric, amendment 2), and the 0.45 Å tolerance misses long secondary bonds
  such as Tl–N 2.74 Å, Pb–S 2.96 Å and axial Cu–Cl 2.90 Å.
- Old CIFs without an occupancy column make disorder alternatives count as full
  weight (`ATOMS_TOO_CLOSE`, `SHORT_HH_CONTACT`).
- `_atom_site_attached_hydrogens` is ignored, so declared implicit H read as
  missing.

**Dependable findings:** missing water or methyl H where the formula includes
them, organic CIFs with no H at all, aberrant Ueq (including a likely 0.437 typo
for 0.0437), and a genuine duplicated label.

## Limitations

One model adjudicator, not blinded to the finding. 2–10 items per check, so the
intervals are wide. Checks the allocation did not reach (`ADP_NPD`,
`DISORDER_SUM_NOT_ONE`, `ISOLATED_H` and rarer ones) are unmeasured. The fixes
above will be made after this commit and must be measured on a fresh sample;
re-scoring these 60 would be fitting.
