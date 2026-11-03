# cifaudit

[![tests](https://github.com/pileofflapjacks1/cifaudit/actions/workflows/ci.yml/badge.svg)](https://github.com/pileofflapjacks1/cifaudit/actions/workflows/ci.yml)

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

### How far to trust each check

Every structural finding carries a **tier**: how often that check has been right
when its findings were sampled. Severity says how bad a finding would be if it
were true; the tier says how often it has been true. Four pre-registered
precision rounds supply the counts: 237 findings from 2,800 COD structures,
each judged against the full CIF by an independent model adjudicator. That is
model adjudication, not expert review.

| tier | rule | checks (right / judged) |
|---|---|---|
| **reliable** | judged ≥ 8 times, right in at least half | `ADP_EXTREME_ANISOTROPY` 14/16, `UEQ_RANGE_WIDE` 13/15, `NO_HYDROGENS_WITH_CARBON` 6/8, `FORMULA_SITES_MISMATCH` 24/40 |
| **prompt** | judged ≥ 8 times, right in fewer than half | `SHORT_HH_CONTACT` 5/16, `DENSITY_BELOW_ONE` 2/8, `SHORT_XH` 3/15, `SHORT_NONBONDED_CONTACT` 2/11, `ATOMS_TOO_CLOSE` 2/12, `H_MULTIPLE_CONNECTIONS` 2/14, `ISOLATED_ATOM` 4/40, `DISORDER_SUM_NOT_ONE` 0/12 |
| **unmeasured** | judged fewer than 8 times, never drawn, or redefined since | everything else, including `ADP_NPD`, `SYMMETRY_POSSIBLY_MISSED` (3/4) and `NON_INTEGRAL_CELL_CONTENT` (0/10, before it was redefined) |

**An expert does not share these tiers.** In a blind, pre-registered review of
20 of the sampled findings, a crystallographer marked **none** of the 7 from
"reliable" checks as worth acting on (the model marked 5). Overall the expert
found 2 of 20 worth acting on, against the model's 6, with no overlap between
the two sets (κ = −0.18). The tiers describe how an independent model judged
these checks. They are not an expert's judgement, and one expert on 20 items is
not enough to replace them. See
[`data/processed/expert_review_results.md`](data/processed/expert_review_results.md).

A prompt is something to look at, not a verdict. The text report lists prompts
in their own section, and they do not affect the exit code unless you pass
`--fail-on-prompts`. Each finding's `track_record` gives the counts behind its
tier, and `tests/test_tiers.py` recomputes every count from the published
verdict files. The MCP response adds a `tier_note`, so an assistant reading it
knows not to report a prompt as an error.

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

Exit codes: `0` clean, `1` a finding at or above `--fail-on` (default `error`)
in **any** file, `2` any file could not be read or parsed. Prompt-tier findings
(see [How far to trust each check](#how-far-to-trust-each-check)) are shown but
count towards exit 1 only with `--fail-on-prompts`. A batch fails if one
file fails; an unreadable file does not stop the rest being audited; and
`--fail-on never` suppresses exit 1 but not exit 2, since an unreadable file
means the tool could not do its job.

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

Verified end to end against `mcp` 2.3.0 (protocol 2025-11-25): all four tools
called over stdio by a real client, including abstention on bad input and the
prompt-injection fixture. Reproduce with `scripts/mcp_smoke.py`; both SDK
generations are supported, since 2.0 renamed `FastMCP` to `MCPServer`.

The implementations are importable directly, with no SDK and no MCP runtime:

```python
from cifaudit.tools import audit_cif, coverage, describe_alert, list_datablocks
list_datablocks(path="mystructure.cif")
```

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

### Held out: 700 structures the tool had never seen

Everything above was fitted to one 500-file sample over four rounds, so 46% was a
fitting statistic. Pre-registered a held-out draw against commit `1b705c9` with
all nineteen thresholds recorded and a commitment to change nothing in response:
700 CIFs from eight COD directories disjoint from both prior samples, zero ID
overlap.

| | in-sample (500) | **held-out (700)** | + bug fixes | + conventions | + round-2 fixes | + round-3 fixes |
|---|---|---|---|---|---|---|
| no structural findings | 46% | **47%** | 49% | 62% | 65% | **67%** |
| crashed / unparseable | 0 / 0 | **0 / 0** | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| seconds per file | 0.69 | 0.44 | 0.38 | 0.25 | 0.20 | 0.20 |

The last column re-runs the same 700 after the four bugs below were fixed
(`d72fecf`); nothing was tuned to them. Structural findings fell by about
two-thirds, mostly because each finding is now reported once rather than once
per symmetry copy. `DUPLICATE_LABELS` fell from 22 files to 2: twenty were
primed labels misread as duplicates.

The conventions column (`354e9fd`) re-runs the same 700 again. Most of the rise
to 62% is `NON_INTEGRAL_CELL_CONTENT` going from 143 files to 0: all 143 had a
formula to compare against, and that check was judged wrong on 10 of 10 sampled
findings. `FORMULA_SITES_MISMATCH` fell from 124 to 92 and `ISOLATED_ATOM` from
71 to 49. This is noise removed, not structures becoming more correct, and none
of these rates is a precision.

The round-2 fixes column re-runs the same 700 once more, after the defects found
by the second precision round were fixed (see below): 432 to 452 structurally
clean. All 16 `SCATTERING_TYPE_LABEL_MISMATCH` files were "Wat" water labels.
`ATOMS_TOO_CLOSE` fell from 40 to 16, nearly all of the difference being
disorder across a symmetry element, and `SYMMETRY_DERIVED_LOWER` fell from 14 to 2,
because gemmi's merging of near-coincident copies had broken the symmetry.
Findings that were new on this run were checked one by one. They are:

- formula counts that now match the declared formula
- partial waters that had been joined by a fake 1.1-1.7 A O-O "bond" and are now
  reported as solvent
- one heavily disordered structure with no occupancy column, whose overlapping
  copies are now visible

It held. Per-check rates move by a few points at most. I had predicted
materially worse and recorded that prediction before running it.

The likely reason: of nineteen parameters only two were set from the sample.
Everything else was a conceptual correction — symmetry must come from the
expanded cell, a hydrogen bond is not a clash, `1.` is a number, `HO1` is not
holmium, a dummy atom is not an atom, a shared site is one site. Principled
fixes generalise; tuned constants do not.

**47% is not a quality statistic about COD** and not a false-positive rate — it
is the fraction of structures with no structural finding. Rates still too high
for correct published work, and so still carrying false positives:
`NON_INTEGRAL_CELL_CONTENT` 20%, `FORMULA_SITES_MISMATCH` 18%, `ISOLATED_ATOM`
12%.

### Precision: when a check fires, is it right?

Firing rates say how often a check fires, not whether it is right. A
pre-registered sample of 60 structural findings from 400 fresh COD structures was
judged against the full CIFs by an independent model adjudicator. **This is model
adjudication, not expert review**; the expert pass has not been done.

**Overall precision: 38% (95% CI 27–51%).** Best: `UEQ_RANGE_WIDE` 4/4,
`NO_HYDROGENS_WITH_CARBON` 3/4, `FORMULA_SITES_MISMATCH` 6/10. Worst:
`NON_INTEGRAL_CELL_CONTENT` 0/10, `ISOLATED_ATOM` 1/10.

It found three bugs — unwrapped symmetry coordinates breaking the distance
search, the ADP anisotropy ratio computed without converting to Cartesian axes,
and primes stripped from atom labels — plus checks that flag normal chemistry.
Those three, and a fourth (findings repeated per symmetry copy), are fixed in
`d72fecf`. The 38% describes the code before the fixes.
Results, reasons and limitations:
[`data/processed/precision_results.md`](data/processed/precision_results.md).

**Round 2, after the fixes and conventions: 36% (95% CI 25–48%).** This used a
pre-registered fresh sample of 60 findings from 800 new COD structures, again
judged by model adjudication. That is no measurable change from round 1. The
fixes removed a failure (`NON_INTEGRAL_CELL_CONTENT` no longer fires), but the
sample now draws on geometric checks that fail on one shared cause: an atom
disordered across a symmetry element is treated as clashing with its own image.
Three more bugs were confirmed by independent recomputation. They were fixed
after the round-2 figure was published, together with the shared cause:

- **Radii.** The covalent-radius table now covers Cordero's full set. The
  isolated-atom test abstains when a neighbour has no radius.
- **Expansion.** cifaudit now does its own symmetry expansion. It merges only
  copies within 0.05 A, so atoms disordered across a symmetry element keep both
  copies. gemmi had merged any copy within about 0.4 A.
- **Disorder across a symmetry element.** A partially occupied atom close to its
  own image is an alternative to that image and to whatever is bonded to the
  image. A negative PART group is tested as a whole.
- **Smaller fixes.** A bonded X-H is no longer an "impossible" clash. "Wat1" is
  not read as an element. `nonbonded_min` is documented as it is applied
  (X...Y and X...H).

On the round-2 sample, 14 of the 38 false findings no longer fire and all 21
real ones still do. Those items were used to design the fixes, so this is not a
precision figure; that needs a third, fresh sample. A known gap remains: an H
two bonds from its own image (a CH2Cl2 about a 2-fold axis, COD 4102011) is
still reported.

Counting accurately worded info findings as correct gives 51% (38–63%). Round 1
was not scored that way, so that figure is not comparable.
[`data/processed/precision_round2_results.md`](data/processed/precision_round2_results.md).

**Round 3, after the round-2 fixes: 26% (95% CI 16–38%).** This used a
pre-registered fresh sample of 60 findings from 800 new COD structures, again
judged by model adjudication. The intervals overlap those of rounds 1 and 2, so
this is no measurable change, and not the improvement the fixes were meant to
bring.

Pooled across the three rounds (exploratory, 177 judged findings), the pattern
is stable:

- **Usually right when they fire:** `ADP_EXTREME_ANISOTROPY` (10/12),
  `UEQ_RANGE_WIDE` (9/11) and `FORMULA_SITES_MISMATCH` (20/30).
- **Mostly wrong:** the geometric and connectivity checks, at 0–27%. They depend
  on reconstructing disorder from labels that CIFs use inconsistently, and on
  fixed distance thresholds that do not hold across the periodic table. A flat
  1.35 Å X–H cutoff misses P–H and every bridging metal hydride.

Treat the first group as signal and the second as prompts to look.
[`data/processed/precision_round3_results.md`](data/processed/precision_round3_results.md).

After the round-3 figure was published, the six defects it confirmed were fixed:

- **X–H reach.** It scales with the heavy atom: r_X + r_H + 0.15 A, never below
  1.35 A. That covers P–H 1.53 A and Lu–H 2.33 A. Only the flat 1.35 A counts
  towards a "multiple connection".
- **Zero occupancy.** `OCC_ZERO` fires only at exactly zero.
- **Disorder closure.** Negative-PART groups are not summed. A group's occupancy
  comes from its partial members.
- **Copy merging.** A negative-PART atom merges with its own image within 0.10 A.
  Every other atom keeps the 0.05 A limit, because genuine copies 0.05–0.15 A
  apart exist.
- **Assembly letters.** Atoms in different positive disorder groups are treated
  as alternatives even when their assembly letters disagree, if their occupancies
  sum to at most 1 and they lie within 2.6 A.

On the 700 held-out files, structures with no structural finding went from 452
to 466. Every changed file was compared before and after: no finding was added,
and each removal was checked (refined occupancies of 0.0001–0.001, Hf hydrides,
an H–Cl bond, mislettered disorder). On the round-3 sample, 21 of 43 false
findings no longer fire and all 15 real ones still do. On the round-2 sample it
is 16 of 38, with all 21 real kept. Both samples were used to design the
fixes, so neither gives a precision figure.

**Round 4, after the round-3 fixes: 40% (95% CI 29–53%).** This used a
pre-registered fresh sample of 60 findings from 800 new COD structures, judged by
model adjudication. The four rounds give 38%, 36%, 26% and 40%, and every pair
of intervals overlaps: **no round has measurably changed overall precision.**
Because the check mix shifts from round to round, the overall figure moves with
which checks get drawn.

Pooled over all four rounds (exploratory, 237 judged findings):

| usually right when they fire | | mostly wrong | |
|---|---|---|---|
| `ADP_EXTREME_ANISOTROPY` | 14/16 | `SHORT_HH_CONTACT` | 5/16 |
| `UEQ_RANGE_WIDE` | 13/15 | `SHORT_XH` | 3/15 |
| `NO_HYDROGENS_WITH_CARBON` | 6/8 | `SHORT_NONBONDED_CONTACT` | 2/11 |
| `SYMMETRY_POSSIBLY_MISSED` | 3/4 | `ATOMS_TOO_CLOSE` | 2/12 |
| `FORMULA_SITES_MISMATCH` | 24/40 | `H_MULTIPLE_CONNECTIONS` | 2/14 |
| | | `ISOLATED_ATOM` | 4/40 |
| | | `DISORDER_SUM_NOT_ONE` | 0/12 |

Reliable or not, it is the same checks every time, as judged by the model. A
blind expert review of 20 of these findings disagreed with the model about
which ones were real, and marked none from the "reliable" checks as worth
acting on
([`expert_review_results.md`](data/processed/expert_review_results.md)). Rather than a fifth
fix-and-measure round, each structural finding now carries this track record as
a tier (see [How far to trust each check](#how-far-to-trust-each-check)). On
the 700 held-out files, 466 (67%) have no structural finding at all, and 528
(75%) have none outside the prompt tier.
[`data/processed/precision_round4_results.md`](data/processed/precision_round4_results.md).

After round 4 was published, the four defects it confirmed were fixed:

- **Zero occupancy.** A site with occupancy 0 is left out of the expansion, so a
  ring-centroid placeholder typed H is not treated as a hydrogen. `OCC_ZERO`
  still reports it.
- **SQUEEZE direction.** The SQUEEZE convention excuses a formula mismatch only
  when the formula has more atoms than the model. When the model has more,
  removed solvent cannot explain it, and the finding stays a warning that says
  so.
- **SQUEEZE in a .fab file.** An embedded `_shelx_fab_file` counts as evidence
  of solvent removal. The data name alone is used; its text is never read.
- **Two bonds from the image.** An atom bonded, through a partially occupied
  atom, to the image that clashes with its partner is now an alternative. This
  covers a methanol across an inversion centre (COD 4342081) and a CH2Cl2 about
  a 2-fold axis (COD 4102011).

On the 700 held-out files, a per-file comparison shows no finding added and six
removed. All six are disorder across a symmetry element, checked one by one.
The file-level figures are unchanged. Every real finding in the round-2, 3 and
4 samples still fires.

### Conventions adopted after the precision study

Each came from the expert adjudicator (J. R. Ramirez) or from a false positive
the model adjudicator explained, and each has a test in
`tests/test_conventions.py`:

- **An isolated O is information, described as probably modelled solvent** —
  water whose H were not located, or O placed to account for electron density.
- **"Isolated" means nothing within r_i + r_j + 0.90 Å**, a coordination-contact
  distance kept separate from the 0.45 Å covalent tolerance, so long secondary
  bonds such as Tl–N 2.74 Å and Pb–S 2.96 Å count.
- **SQUEEZE turns a formula mismatch into information**, detected from
  `_platon_squeeze_*` (or `_smtbx_masks_*`) data names, never from free text.
- **`_atom_site_attached_hydrogens` counts towards composition**, and implicit H
  is not reported as missing.
- **Non-integral cell contents are reported only when there is no formula to
  compare against.** A declared solid solution is not a defect.
- **No occupancy column** softens a close contact to a warning that says the
  atoms may be disorder alternatives.
- **Deuterium** is named D in messages, and an X–D distance is not judged when
  the CIF declares no radiation probe, since X-ray and neutron targets differ.

These were designed against the precision sample, so that sample can no longer
measure them. A fresh one is needed.

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

The ratio is computed in **Cartesian axes**, from
`U_cart = (A N) U_cif (A N)^T` with `A` the orthogonalisation matrix and
`N = diag(a*, b*, c*)`. The raw CIF tensor gives the right *sign* for the NPD
test but not the right *ratio*: in non-orthogonal cells it runs high (one atom
in the precision study read 4.10 raw and 2.76 in Cartesian axes). An earlier
version made exactly this mistake.

The threshold is set from the deposited corpus. Over **12,899 non-H anisotropic
atoms in 500 COD structures**:

| percentile | sqrt(U3/U1), Cartesian |
|---|---|
| p50 | 1.49 |
| p95 | 2.48 |
| p99 | **3.55** |
| p99.9 | 8.49 |

`adp_axis_ratio_max` is **3.55**, the 99th percentile, so it flags the top ~1% of
atoms. It was previously 4.0, the p99 of the raw ratios; before that, 10 on the
raw eigenvalue ratio. Findings also report
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

119 tests. Every planted defect has a fixture; the clean reference must produce
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
- The MCP server needs Python ≥ 3.10 (every published `mcp` release does), while
  the library and CLI run on 3.8+. It **has** now been run end to end against
  the real SDK — `mcp` 2.3.0, protocol 2025-11-25, all four tools called over
  stdio by a real client — via `scripts/mcp_smoke.py`. That is a script rather
  than a pytest test because the test environment is 3.8.

## Not affiliated

Not affiliated with or endorsed by the IUCr, the CCDC, or the authors of PLATON,
checkCIF, SHELX or Olex2. "checkCIF" and "PLATON" are referred to descriptively.
