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

## Tests

```bash
pip install -e ".[test]" && pytest
```

67 tests. Every planted defect has a fixture; the clean reference must produce
zero findings. Two tests exist because they caught real bugs:

- `test_geminal_pairs_are_not_short_contacts` — the first implementation flagged a
  normal geminal H···O at 1.995 Å in formaldehyde as a short contact, because it
  excluded only bonded pairs and not 1,3 pairs.
- `test_no_alert_level_is_asserted` — guards against drifting back into claiming
  an A/B/C/G level that the source does not support.

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
