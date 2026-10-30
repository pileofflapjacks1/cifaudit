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

---

# Out-of-sample run: 500 COD structures

2026-10-30, after the 16-file in-sample work above. **This is the run that
counts**, and it began by showing the in-sample result was indeed measuring fit.

## Sample

500 CIFs pulled via **rsync** — `rsync://www.crystallography.net/cif/` — which is
the channel COD documents for bulk use. COD's `robots.txt` is a blanket
`User-agent: * / Disallow: /`, so automated HTTP fetching was not used. Six
directories, deliberately disjoint from the 16 earlier files, spread across the
ID space for chemical diversity: `1/00/02`, `2/30/50` (empty), `4/32/10`,
`7/10/50`, `8/10/20`, `9/01/00`. Zero overlap with the earlier sample.

Not committed to this repository: no point redistributing COD's data.

## First result: the in-sample number did not hold

| | in-sample (16) | out-of-sample (500) |
|---|---|---|
| files with no findings | 88% | **8%** |
| errors | 0 | **1,055** |
| warnings | 0 | 3,887 |

## Bug six: the number parser rejected valid CIF

One bug accounted for 207 findings. `parse_number` required a digit after the
decimal point, so it rejected CIF's bare-dot form — `1.`, `0.`, `.5` — which is
valid and ubiquitous in older entries. Occupancies written `1.` read as
unparseable (99 files); `U13`/`U23` written `0.` read as *absent*, reporting
complete ADPs as incomplete (33 files). All 14 hand-written fixtures used
`0.2000` style, so nothing but real data could have caught it.

`OCC_UNPARSEABLE` 99 → 0 files. `ADP_UIJ_INCOMPLETE` 33 → 0.

## Then two conventions, both of which cifaudit had wrong

Researched against published practice rather than chosen (see the README section
"Conventions, and where they come from"):

- **Bonding** is an additive tolerance, `r_i + r_j + epsilon`, with epsilon
  0.35–0.45 Å — not the multiplicative 1.2× cifaudit used. Hydrogen is excluded
  from radius-based bonding and bonds to its single nearest heavy atom.
- **ADP anisotropy** is `SQRT(U3/U1)` on the principal axes (PLAT213), not the
  raw eigenvalue ratio. The old default of 10 raw equals 3.16 in axis units —
  only the 97.5th percentile of 12,899 real atoms. Now 4.0, the p99.
- **`UEQ_RANGE_WIDE`** belongs per element type, which is what PLAT220 is for:
  catching misassigned species.

## Where it stands

| | pre-fix | + parser fix | + conventions |
|---|---|---|---|
| no structural findings | — | — | **207/500 (41%)** |
| no findings of any kind | 42/500 (8%) | 42/500 (8%) | 47/500 (9%) |
| structural errors | 1,055 | 894 | 885 |
| seconds per file | 1.73 | 1.85 | **0.35** |

Per-check movement on the checks the conventions touched:

| check | pre | + parser | + conventions |
|---|---|---|---|
| `OCC_UNPARSEABLE` | 99 | **0** | 0 |
| `ADP_UIJ_INCOMPLETE` | 33 | **0** | 0 |
| `ADP_EXTREME_ANISOTROPY` | 114 | 121 | **61** |
| `UEQ_RANGE_WIDE` | 36 | 36 | **7** |
| `ISOLATED_H` | 14 | 14 | **10** |
| `SHORT_NONBONDED_CONTACT` | 24 | 24 | **20** |
| `SHORT_HH_CONTACT` | 30 | 30 | **27** |
| `ATOMS_TOO_CLOSE` | 137 | 137 | **137** |

The 5× speedup was incidental: the hydrogens check had its own O(H × heavy)
distance loop, and routing it through the shared bond graph removed the hot spot.

**41% is the honest headline, and it is not a false-positive rate.** It is the
fraction of deposited structures with no *structural* finding. Separating
reporting from structural is what made it legible at all: 74% of files omit
`_chemical_formula_moiety` and 40% omit a standard uncertainty on the cell
volume, and while those are real checkCIF tests, counting them next to NPD atoms
buried the signal.

## The largest remaining false-positive source: shared sites

`ATOMS_TOO_CLOSE` did not move: **137/500 (27%)**, and it is not a threshold.
Diagnosing the hydrogen findings printed `Sb1 0.00 A, W1 0.00 A` — one H
coincident with two *different* heavy atoms, because Sb and W **share one
crystallographic site** with partial occupancies. That is routine mixed-occupancy
inorganic chemistry, and cifaudit does not model shared sites at all, so every
one reads as atoms sitting on top of each other. The earlier deduplication only
caught identical labels at identical positions.

Until that is modelled, treat `ATOMS_TOO_CLOSE`, `H_MULTIPLE_CONNECTIONS` (13%)
and some of `FORMULA_SITES_MISMATCH` (19%) as unreliable on structures with
mixed-occupancy sites.

## Still in-sample

Fixes six, and both conventions, were applied in response to **these** 500 files.
The 41% is therefore an in-sample figure again. A third independent draw is the
test that would confirm it, and it has not been done — it now costs about three
minutes of compute, so there is no excuse other than that the conventions should
settle first.
