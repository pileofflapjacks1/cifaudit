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
    "adp_ratio_max": 10.0,      # max/min eigenvalue of U(i,j)
    "ueq_ratio_max": 10.0,
    "ueq_min_atoms": 3,
    "symprec": 0.05,            # Angstrom, passed to spglib
    "hard_floor": 0.75,         # A; below this two atoms cannot both be real
    "nonbonded_min": 2.2,       # A; shortest acceptable non-bonded heavy contact
    "hh_min": 2.0,              # A; shortest acceptable H...H contact
    "bond_factor": 1.2,         # d <= bond_factor * (r_cov_i + r_cov_j) => bonded
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

    if not blk.sites and not blk.cell.complete():
        res.status = NOT_ENOUGH_DATA
        res.notes.append("Neither a complete unit cell nor any atom sites were "
                         "found, so almost nothing could be tested.")

    res.notes.append("Thresholds used: %s" % ", ".join(
        "%s=%s" % (k, ctx[k]) for k in sorted(ctx)))
    return res
