# cifaudit against real deposited structures

First contact with real data, 2026-10-30. Until this run, every CIF cifaudit had
ever seen was a fixture I wrote myself, including the "clean" reference. **That
turned out to be hiding five bugs.**

## Sample

16 pristine CIFs from the Crystallography Open Database, with COD provenance
headers intact, originally deposited via IUCr Journals. COD IDs: 1543089,
1543227, 1544230, 1544651, 1549607, 1550274, 2012000, 2013000, 2016430, 2016452,
2017775, 2100301, 2200000, 2221836, 2227862, 2233297.

**The sample is not random.** These were already on disk as the benchmark set for
an unrelated project (phase retrieval), so they are biased toward small-molecule
structures with deposited structure factors. Rates below are not population
estimates.

Reproduce with:

```bash
python3 scripts/run_cod_benchmark.py <dir-of-cifs> --json out.json
```

## Result

| | before | after |
|---|---|---|
| files with ≥1 finding | **16/16** | **2/16** |
| files completely clean | 0/16 | **14/16 (88%)** |
| total findings | 36 (19 warn, 17 info) | **3 (0 warn, 3 info)** |
| errors | 0 | 0 |
| crashes | 0 | 0 |
| unparseable | 0 | 0 |

Per-check firing rate after the fixes:

| check | files | rate |
|---|---|---|
| `SPACEGROUP_SETTING_DIFFERS` | 1/16 | 6% |
| `NON_INTEGRAL_CELL_CONTENT` | 1/16 | 6% |
| `PARTIAL_OCC_NO_DISORDER_GROUP` | 1/16 | 6% |

All three survivors look like true positives: COD 1543227 declares `P 1 21/n 1`
where spglib's standard setting for No. 14 is `P2_1/c`; COD 2017775 has three
partial-occupancy sites with no disorder-group labels and non-integral cell
contents.

## The five bugs real data found

**1. The missed-symmetry check fired on 16/16 — a 100% false-positive rate.**
spglib *derives* symmetry from the atoms it is given, and it was being given the
asymmetric unit, from which it correctly returns P1 for almost any real
structure. Expanding to the full unit cell first recovers the declared
space-group number on all 16. Same-number/different-symbol is now reported
separately as a setting difference (`PLAT128`), not as missed symmetry.

**2. Every short-contact finding was a hydrogen bond.** All 18 were H···O or
H···N between 1.3 and 2.2 Å — textbook hydrogen bonding. A D–H···A approach,
where the H is bonded to a donor, is now excluded; only an implausibly short one
(< 1.2 Å) is reported.

**3. 1,4 contacts were being reported as clashes.** Excluding bonded (1,2) and
geminal (1,3) pairs was not enough: a 1,4 C···H across a planar aromatic ring
runs about 2.18 Å. Exclusion is now by topological distance in the bond graph,
default 3 bonds.

**4. `HO1` was read as holmium.** COD 2013000 labels a hydroxyl hydrogen `HO1`.
Deriving the element from the label greedily matched two characters and produced
a scattering-type mismatch against its own `H`. A label is now consistent with
its type symbol when it *starts with* that symbol.

**5. Species numbering used `hash()` on element symbols.** Python salts string
hashing per process, so spglib's input — and therefore cifaudit's output — could
differ between runs of a tool whose entire claim is determinism. Now a sorted
index. There is a test that runs the audit under three `PYTHONHASHSEED` values
and requires byte-identical output.

A sixth change was not a bug: `MULTIPLE_FRAGMENTS` held for 14/16 files, because
molecular and ionic crystals have many fragments by construction. It was never a
defect indicator and is now a note rather than a finding.

Each of the five has a regression test in `tests/test_real_data_regressions.py`.

## What this does and does not show

**Sensitivity was not traded away.** All seven planted-defect fixtures still
fire after the changes — occupancy > 1, NPD ADPs, duplicate labels, cell-volume
mismatch, atoms below the hard floor, a hydrogen treatment declared with no H,
and a disorder assembly that does not close.

**But this is an in-sample result, and that is the main caveat.** The thresholds
and exclusion rules were changed *in response to* these 16 files, so "88% clean"
measures how well the tool was fitted to them, not how it behaves on structures
it has not seen. The defence is that each change is principled rather than a
threshold nudge — expand before deriving symmetry; a hydrogen bond is not a
clash; topological neighbours are not packing contacts; a label starting with
`H` is hydrogen — but **an out-of-sample run on a few hundred independently
drawn COD entries is the test that would actually count, and it has not been
done.**

**Speed is a known problem.** 6.1 s per file here, and the worst single file in
an earlier run took 16 s. Distance checks are O(n²) over 27 lattice images with
no spatial indexing, and abstain above 400 unit-cell sites rather than
truncating. A few hundred files is tolerable; a database sweep is not.

## Still untested against real data

Every `not_tested` area, plus: multi-datablock CIFs (all 16 had exactly one),
macromolecular CIFs, powder data, modulated or twinned structures, CIFs above
the 400-site cap, and anything that fails to parse — no file in this sample was
unparseable, so that path has only ever run against fixtures.
