# cifaudit

**Offline consistency checks on a crystallographic CIF, as a library, a CLI and an
MCP server.**

Version 0.1.0 · Python ≥ 3.8 · MIT

> **cifaudit is not checkCIF.** It implements a small subset of the
> checkCIF/PLATON tests, asserts no alert level, and tells you in every response
> which tests it did not perform. Run
> [checkCIF](https://checkcif.iucr.org/) before you publish. **A clean cifaudit
> report does not mean the structure is correct.**

## Why this exists

[checkCIF](https://checkcif.iucr.org/) is the authority on CIF validation, and it
is a web form. There is no documented API and no documented batch access; the
report comes back as HTML or PDF, with "PDF by email" offered for slow checks. An
agent cannot call it, and using it means uploading a structure — possibly
unpublished — to a third party.

Meanwhile, of the 31 materials-and-chemistry MCP servers catalogued in
[awesome-materials-ai](https://github.com/dfieser/awesome-materials-ai), none
validates a crystal structure, and none reads the experimental quantities that
live in a deposited CIF: occupancies, anisotropic displacement parameters,
disorder groups. They treat a CIF as an idealised structure to compute on.

cifaudit treats it as an experimental report that might be wrong.

## What it checks

Nine checks, each mapped to the checkCIF/PLATON tests it corresponds to:

| Check | Looks for |
|---|---|
| `cell` | missing or non-positive cell parameters, angles outside (0, 180), a cell that is not realisable, `_cell_volume` disagreeing with the parameters |
| `labels` | duplicate or empty atom-site labels, `_atom_site_type_symbol` disagreeing with the label |
| `occupancy` | occupancy > 1, < 0, or = 0; partial occupancy with no disorder bookkeeping; disorder groups not summing to 1 |
| `adp` | **NPD atoms**, extreme U(i,j) anisotropy, incomplete U(i,j), `Uani` declared with no aniso loop entry, negative Uiso, wide Ueq spread |
| `symmetry` | missing or uninterpretable space group, missing symmetry operators, operator count disagreeing with the group, and **higher symmetry that the coordinates support but the CIF does not declare** (ADDSYM-equivalent, via spglib) |
| `contacts` | atoms below a hard floor, sites on their own lattice image, short non-bonded and H···H contacts — with **1,2 and 1,3 (geminal) pairs excluded**, as PLATON excludes them |
| `hydrogens` | no H in a carbon compound, a hydrogen treatment declared with no H present, isolated H, H bonded to more than one heavy atom, short X–H |
| `fragments` | isolated atoms (metals and O/N weighted more seriously), disconnected fragments in the cell |
| `formula` | `_chemical_formula_sum` against the symmetry-expanded atom-site list divided by Z, non-integral cell contents, density below 1.0 |

### The NPD criterion

The anisotropic Debye–Waller factor is
`T(h) = exp(-2π² Σ U_ij h_i h_j a*_i a*_j)`, so the quadratic form acts on
reciprocal indices scaled by the reciprocal axis lengths. That scaling is a
congruence by a positive diagonal matrix, which preserves positive-definiteness.
**Testing the eigenvalues of U_ij exactly as given in the CIF is therefore the
correct NPD test for any cell**, with no transformation to Cartesian axes.

## What it does *not* check

Every response carries a `not_tested` list. That is the point, not an apology.
Structure factors, refinement instruction files, solvent-accessible voids,
Hirshfeld rigid-bond tests, absolute structure, special-position site symmetry,
CCDC duplicate detection and network dimensionality are **not implemented at
all**. The IUCr list contains **526** PLATON tests; cifaudit references 60 of
them and implements a subset of those.

```
$ cifaudit mystructure.cif
...
NOT TESTED (8 areas):
  - occupancy < 1.0 for an atom on a special position
      requires site-symmetry analysis at each site, which is not implemented
  - Hirshfeld rigid-bond test; ADP restraint reporting; negative density at H positions
      requires refinement s.u. propagation and/or the structure-factor file, which
      cifaudit does not read
  ...
```

## Install

**Not on PyPI yet.** Install from the repository:

```bash
pip install "git+https://github.com/pileofflapjacks1/cifaudit"
pip install "cifaudit[symmetry] @ git+https://github.com/pileofflapjacks1/cifaudit"   # + spglib
pip install "cifaudit[mcp] @ git+https://github.com/pileofflapjacks1/cifaudit"        # + MCP server (needs Python >= 3.10)
```

Or from a clone:

```bash
git clone https://github.com/pileofflapjacks1/cifaudit
cd cifaudit
pip install -e ".[symmetry]"
```

Without `spglib` the missed-symmetry check **abstains** instead of silently
passing.

## Use

```bash
cifaudit structure.cif                        # human-readable report
cifaudit structure.cif --json                 # machine-readable
cifaudit multi.cif --list-datablocks          # do this first on an unknown file
cifaudit multi.cif --datablock second
cifaudit structure.cif --set hh_min=1.8       # override any threshold
cifaudit *.cif --fail-on warn                 # for CI
```

Exit codes: `0` clean, `1` findings at or above `--fail-on` (default `error`),
`2` the file could not be read or parsed.

```python
from cifaudit import audit
result = audit(open("structure.cif").read()).to_dict()
result["status"]        # "ok" | "not_enough_data"
result["findings"]      # each with .checkcif, .based_on, .sites, .evidence
result["not_tested"]    # always populated
```

### MCP

```bash
cifaudit-mcp            # stdio
```

| Tool | Purpose |
|---|---|
| `list_datablocks` | **Call first on an unfamiliar file.** A CIF may hold many structures; auditing the first silently is how you report on the wrong one. |
| `audit_cif` | Run the checks on one datablock. |
| `describe_alert` | Look up a PLAT code; returns PLATON's verbatim description. |
| `coverage` | What is and is not implemented, so a caller can decide whether checkCIF is required. |

## Design

**Abstention is a first-class result.** `status` is `ok` or `not_enough_data`, and
it describes whether the checks could run — never whether the structure is good.
Any check that cannot reach a verdict says so in `not_tested` with a reason,
rather than passing silently.

**No alert level is asserted.** checkCIF levels are A/B/C/G. The public PLATON
test list classifies each test as `Type_1`…`Type_5` and publishes no mapping
between the two, so cifaudit reports the PLAT code and PLATON's own `Type_N`
verbatim, alongside its own separate `error`/`warn`/`info` severity. Descriptions
are quoted from the
[IUCr list of PLATON validation tests](https://journals.iucr.org/services/cif/checking/platon_tests.html)
(read 2026-10-30).

**Every finding shows its work:** the CIF data names it used (`based_on`), the
sites involved, and the numbers (`evidence`). Thresholds are explicit, listed in
every response, and all overridable. They are cifaudit's defaults, not standards.

**Deterministic and offline.** No network access. No language model in the
validation path. A test asserts the package never references an HTTP client.

**A CIF is untrusted input.** It is a text file that may have come from anywhere,
and its free-text fields can contain prose addressed to a language model.
Therefore input is size-capped before parsing; free-text fields are never read,
returned or echoed; atom labels and space-group symbols are sanitised (control
characters stripped, length capped) before appearing in a report; and nothing in
a CIF is ever treated as an instruction. There is a prompt-injection test.

## Against real structures

Until 2026-10-30 every CIF cifaudit had seen was a fixture written by its author,
including the "clean" reference. The first run over **16 deposited COD
structures** produced findings on **16/16 files** and exposed five bugs.

| | before | after |
|---|---|---|
| files with ≥1 finding | **16/16** | **2/16** |
| files completely clean | 0/16 | **14/16 (88%)** |
| total findings | 36 (19 warn, 17 info) | **3 (0 warn, 3 info)** |
| errors / crashes / unparseable | 0 / 0 / 0 | 0 / 0 / 0 |

What real data found:

1. **The missed-symmetry check fired on 16/16.** spglib derives symmetry from
   the atoms it is given, and it was given the asymmetric unit — from which it
   returns P1 for almost any real structure. Expanding first recovers the
   declared space group on all 16.
2. **Every short-contact finding was a hydrogen bond** — all 18 were H···O or
   H···N between 1.3 and 2.2 Å.
3. **1,4 contacts were reported as clashes.** A 1,4 C···H across a planar
   aromatic ring is ~2.18 Å; exclusion now runs on topological distance.
4. **`HO1` was read as holmium** instead of a hydroxyl hydrogen.
5. **Species numbering used `hash()`**, which Python salts per process — so
   output could vary between runs of a tool whose claim is determinism.

All seven planted-defect fixtures still fire, so sensitivity was not traded for
the lower false-positive rate. Each bug has a regression test.

### Then 500, drawn independently

The in-sample caveat held. Pulled 500 more via rsync (COD's `robots.txt`
disallows automated HTTP; rsync is the channel they document for bulk use), from
six directories disjoint from the first 16: **8% clean, 1,055 errors.** 88% had
been measuring fit.

That run found a sixth bug and two wrong conventions:

- **The number parser rejected valid CIF.** `1.`, `0.` and `.5` are ubiquitous in
  older entries; requiring a digit after the point made occupancies unparseable
  (99 files) and read `0.` ADP components as absent (33 files). Both now 0.
- **Bonding used a multiplicative factor**; the field uses an additive tolerance.
- **ADP anisotropy used the raw eigenvalue ratio**; PLAT213 uses its square root.

Where it stands on those 500:

| | pre | + parser fix | + conventions |
|---|---|---|---|
| no structural findings | — | — | **207/500 (41%)** |
| structural errors | 1,055 | 894 | 885 |
| `ADP_EXTREME_ANISOTROPY` | 114 | 121 | **61** |
| `UEQ_RANGE_WIDE` | 36 | 36 | **7** |
| seconds per file | 1.73 | 1.85 | **0.35** |

**41% is the fraction of deposited structures with no *structural* finding** — not
a false-positive rate. And it is in-sample again, because those fixes were made
in response to these 500 files. A third independent draw is the test that would
confirm it; it costs about three minutes now and has not been done.

### Shared sites

Two elements can occupy one crystallographic position with partial occupancies —
routine inorganic chemistry, and modelling them as co-existing atoms made every
such site read as atoms 0.000 Å apart. They are now treated as mutually
exclusive **alternatives**: never bonded, excused from distance tests, but still
required to have occupancies that can coexist (`SHARED_SITE_OCC_GT_ONE`). Atoms
flagged `_atom_site_calc_flag = dum` are dropped from geometry and kept in the
composition.

| | before | after |
|---|---|---|
| no structural findings | 207/500 (41%) | **231/500 (46%)** |
| structural errors | 885 | **114** |
| `ATOMS_TOO_CLOSE` | 137 (27%) | **30 (6%)** |
| `H_MULTIPLE_CONNECTIONS` | 65 | **30** |

The first attempt at this made things *worse* (36% clean) in two ways worth
knowing: a disordered C–H at 0.930 Å with both occupancies 0.5 sums to exactly
1.0 and was classified as the hydrogen being an alternative to its own carbon, and
`_atom_site_symmetry_multiplicity` turns out to carry two incompatible meanings in
the wild. Both are described in the write-up, and both have regression tests.

Residual false positives against the previous baseline: `ISOLATED_H` 21 vs 10,
`SHORT_HH_CONTACT` 34 vs 27, `ISOLATED_ATOM` 38 vs 33 — not yet checked case by
case.

Full write-up in
[`data/processed/cod_benchmark.md`](data/processed/cod_benchmark.md).

```bash
python3 scripts/run_cod_benchmark.py <dir-of-cifs> --json out.json
```

## Conventions, and where they come from

Two judgements decide most of cifaudit's behaviour. Both follow published
practice rather than taste, and both were wrong in earlier versions.

### Bonding: additive tolerance, and hydrogen handled separately

```
d(i,j) <= r_i + r_j + epsilon          epsilon = 0.45 A  (bond_tolerance)
```

An **additive** tolerance is what the field uses. A set of bonding radii derived
from 519,021 COD records ([arXiv 2601.02017](https://arxiv.org/abs/2601.02017))
states epsilon is *"usually ranging between 0.35 Å and 0.45 Å"* and uses 0.45 Å;
Jmol-style/MaterialsCoord use 0.45 Å; the CCDC's recommended figure is about
0.4 Å; CoRE MOF 2019 uses Cordero radii plus 0.2 Å for metals under 1.95 Å and
0.4 Å above. cifaudit previously used a **multiplicative 1.2×**, which no
standard tool uses.

**Hydrogen is excluded from radius-based bonding**, because it does not work for
H — the standard epsilon would make H–O bonding *more* generous (1.42 Å) than the
old factor did (1.16 Å). The COD-derived radii set above has no hydrogen value at
all and discards every contact involving hydrogen. So an H is bonded to exactly
one heavy atom, its nearest within `xh_bond_max` (1.35 Å; metal hydrides reach
~1.8 Å and need this raised). Two or more candidates is reported as
information, following PLAT303: *"Hydrogen atoms are generally connected to only
one other atom. A hydrogen atom between two oxygen atoms is a special case.
Investigate whether this hydrogen atom is better described with a disorder model
with two partially occupied sites."* H···H is never a bond; dihydrogen complexes
are not modelled.

### ADP anisotropy: the square root of the eigenvalue ratio

PLATON's PLAT213, verbatim: *"The main axes values of the ADP(S) of the main
residue(s) are determined and ordered: U1 < U2 < U3. The value of **SQRT(U3/U1)**
main axis ADP ratio (Angstrom Units) is tested."* U is a mean-square displacement
in Å², so the square root is the ratio of ellipsoid **axis lengths**.

The threshold is set from the deposited corpus. Over **12,899 non-H anisotropic
atoms in 500 COD structures**:

| percentile | sqrt(U3/U1) |
|---|---|
| p50 | 1.66 |
| p95 | 2.81 |
| p99 | **4.01** |
| p99.9 | 9.56 |

`adp_axis_ratio_max` is **4.0** — the 99th percentile, so it flags the top ~1% of
atoms. The old default of 10 was on the *raw* ratio, equal to sqrt = 3.16, which
is only the 97.5th percentile and fired on 24% of files. Findings also report
PLATON's oblate/prolate distinction (*"Oblate criterium: U3 - U2 < U2 - U1.
Prolate otherwise."*).

`UEQ_RANGE_WIDE` is now compared **within one element type**, which is what
PLAT220 is for: *"Too high or too low Ueq's may be an indication for incorrectly
identified atomic species (i.e. O versus N)."* Comparing across all non-H atoms
asked a different question.

## Tests

```bash
pip install -e ".[test]" && pytest
```

97 tests. Every planted defect has a fixture; the clean reference must produce
zero findings. Several exist because they caught real bugs:

- `test_geminal_pairs_are_not_short_contacts` — the first implementation flagged a
  normal geminal H···O at 1.995 Å in formaldehyde as a short contact, because it
  excluded only bonded pairs and not 1,3 pairs.
- `test_no_alert_level_is_asserted` — guards against drifting back into claiming
  an A/B/C/G level that the source does not support.
- `tests/test_real_data_regressions.py` — eight tests, one per bug that a real
  deposited structure exposed (see above), including a determinism test that
  runs the audit under three `PYTHONHASHSEED` values and demands byte-identical
  output.

## Limitations

- The missed-symmetry comparison passes spglib the asymmetric unit and ignores
  disorder and partial occupancy. Treat a `SYMMETRY_DERIVED_DIFFERS` finding as a
  prompt to run PLATON's ADDSYM, not as a conclusion.
- Distance checks enumerate lattice images from −1 to +1 and are capped at 400
  unit-cell sites; larger structures abstain rather than truncate.
- Covalent radii are tabulated for a subset of elements
  ([Cordero et al., *Dalton Trans.* 2008, 2832](https://doi.org/10.1039/B801115J)).
  Elements outside that table cause the dependent checks to abstain.
- `cifaudit.server` is the one module not exercised by the test suite: the `mcp`
  SDK requires Python ≥ 3.10 and the machine this was written on has 3.8.

## Not affiliated

Not affiliated with or endorsed by the IUCr, the CCDC, or the authors of PLATON,
checkCIF, SHELX or Olex2. "checkCIF" and "PLATON" are referred to descriptively.
