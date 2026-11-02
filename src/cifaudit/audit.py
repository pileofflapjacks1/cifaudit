"""Orchestration: run every check over one datablock and aggregate."""

from typing import Any, Dict, List, Optional

from .checks import (adp, cell, contacts, formula, fragments, hydrogens, labels,
                     occupancy, symmetry)
from .load import CifError, Block, datablock_names, load
from .model import AuditResult, NOT_ENOUGH_DATA, NotTested, OK

CHECKS = (cell, labels, occupancy, adp, symmetry, contacts, hydrogens,
          fragments, formula)

# Every threshold is explicit and overridable. None is a standard; they are
# cifaudit's defaults and are stated in the result so a reader can disagree.
DEFAULTS: Dict[str, Any] = {
    "volume_rel_tol": 0.005,    # fractional disagreement allowed on _cell_volume
    "occ_tol": 1e-3,
    "disorder_sum_tol": 0.02,
    # PLATON PLAT213 tests SQRT(U3/U1) of the principal axes, not the raw
    # eigenvalue ratio: U is a mean-square displacement in A^2, so the square
    # root is the ratio of ellipsoid AXIS LENGTHS in A. Computed in Cartesian
    # axes (see checks/adp.py). Set as the 99th percentile over 12,899 non-H
    # anisotropic atoms in the 500-file calibration sample: p50 1.49, p95 2.48,
    # p99 3.55. An earlier value of 4.0 was the p99 of RAW ratios, which run
    # high in non-orthogonal cells.
    "adp_axis_ratio_max": 3.55,
    "ueq_ratio_max": 10.0,
    "ueq_min_atoms": 3,
    "symprec": 0.05,            # Angstrom, passed to spglib
    "hard_floor": 0.75,         # A; below this two atoms cannot both be real
    "nonbonded_min": 2.2,       # A; shortest acceptable non-bonded contact with at least one non-H atom (X...Y and X...H; H...H uses hh_min)
    "hh_min": 2.0,              # A; shortest acceptable H...H contact
    "hbond_min": 1.2,           # A; below this an H...acceptor approach is
                                # too short even for a hydrogen bond
    "bond_tolerance": 0.45,     # A; d <= r_cov_i + r_cov_j + tol => bonded
    "contact_tolerance": 0.90,  # A; an atom with a non-H neighbour within
                                # r_i + r_j + this is not "isolated" (secondary
                                # coordination such as Tl-N 2.74, Pb-S 2.96)
    "shared_site_max": 1.0,     # A; partial-occupancy atoms this close, whose
                                # occupancies cannot both be full, are treated
                                # as mutually exclusive alternatives
    "coincident_max": 0.10,     # A; closer than this counts as the SAME site
    "shared_site_bond_fraction": 0.5,
                                # alternatives must be closer than this fraction
                                # of (r_i + r_j): atoms sharing a site are never
                                # a bond length apart
    "xh_bond_max": 1.35,        # A; an H bonds to its nearest heavy atom within
                                # this. Metal hydrides reach ~1.8 A and will
                                # need this raised.
    "exclude_within_bonds": 3,  # skip contacts between atoms <=N bonds apart
    "xh_short_tol": 0.10,       # A below the X-ray X-H target before reporting
    "formula_tol": 0.05,
}


def audit(text: str, datablock: Optional[str] = None,
          options: Optional[Dict[str, Any]] = None) -> AuditResult:
    ctx = dict(DEFAULTS)
    if options:
        unknown = sorted(set(options) - set(DEFAULTS))
        if unknown:
            raise CifError("unknown option(s): %s" % ", ".join(unknown))
        ctx.update(options)

    blk: Block = load(text, datablock)
    res = AuditResult(status=OK, datablock=blk.name)

    names = datablock_names(text)
    if len(names) > 1 and datablock is None:
        res.notes.append(
            "This CIF has %d datablocks (%s). Only %r was audited, because none "
            "was requested. Pass a datablock name to audit another."
            % (len(names), ", ".join(names[:8]), blk.name))

    for mod in CHECKS:
        try:
            found, nt = mod.run(blk, ctx)
        except Exception as exc:       # a check must never take the audit down
            res.not_tested.append(NotTested(
                area=mod.NAME,
                reason="check raised %s: %s" % (type(exc).__name__, str(exc)[:160])))
            continue
        res.checked.append(mod.NAME)
        res.findings.extend(found)
        res.not_tested.extend(nt)
    for note in ctx.pop("_notes", []):
        res.notes.append(note)

    # Checks that walk the symmetry-expanded cell meet every copy of a site, so
    # one short N-D bond came back four times. Symmetry-equivalent copies give
    # byte-identical findings; keep the first of each. This also stops counts
    # and firing rates being inflated by a structure's multiplicity.
    seen = set()
    unique = []
    for f in res.findings:
        key = (f.code, f.message, tuple(f.sites))
        if key in seen:
            continue
        seen.add(key)
        unique.append(f)
    res.findings = unique

    if not blk.sites and not blk.cell.complete():
        res.status = NOT_ENOUGH_DATA
        res.notes.append("Neither a complete unit cell nor any atom sites were "
                         "found, so almost nothing could be tested.")

    res.notes.append("Thresholds used: %s" % ", ".join(
        "%s=%s" % (k, ctx[k]) for k in sorted(ctx)))
    return res
