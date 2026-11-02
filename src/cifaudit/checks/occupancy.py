"""Site occupancies and disorder bookkeeping."""

from collections import defaultdict
from typing import List, Tuple

from ..model import ERROR, Finding, INFO, NotTested, WARN

NAME = "occupancy"


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    if not blk.sites:
        nt.append(NotTested(area="occupancy", reason="no atom sites in the CIF",
                            plat=["PLAT075", "PLAT074"]))
        return (f, nt)

    if "_atom_site_occupancy" not in blk.present:
        nt.append(NotTested(
            area="occupancy",
            reason="_atom_site_occupancy is absent; occupancies are conventionally "
                   "1.0 when omitted, but that is an assumption, so no occupancy "
                   "test was run",
            plat=["PLAT074", "PLAT075", "PLAT076", "PLAT077"]))
        return (f, nt)

    tol = ctx["occ_tol"]
    for s in blk.sites:
        if s.occupancy is None:
            f.append(Finding(code="OCC_UNPARSEABLE", severity=WARN,
                             message="Site %s has an occupancy that could not be "
                                     "parsed as a number." % s.label,
                             plat=["PLAT075"], based_on=["_atom_site_occupancy"],
                             sites=[s.label]))
            continue
        if s.occupancy > 1.0 + tol:
            f.append(Finding(code="OCC_GT_ONE", severity=ERROR,
                             message="Site %s has occupancy %.4f, greater than 1.0."
                                     % (s.label, s.occupancy),
                             plat=["PLAT075"], based_on=["_atom_site_occupancy"],
                             sites=[s.label], evidence={"occupancy": s.occupancy}))
        elif s.occupancy < 0.0:
            f.append(Finding(code="OCC_NEGATIVE", severity=ERROR,
                             message="Site %s has negative occupancy %.4f."
                                     % (s.label, s.occupancy),
                             plat=["PLAT074"], based_on=["_atom_site_occupancy"],
                             sites=[s.label], evidence={"occupancy": s.occupancy}))
        elif s.occupancy == 0.0:
            # Exactly zero only. A refined 0.001 (COD 9007098: FeM 0.00100,
            # formula Fe0.002) is a real, tiny occupancy; the old test,
            # |occ| <= occ_tol, caught it and then reported "occupancy 0.0".
            f.append(Finding(code="OCC_ZERO", severity=WARN,
                             message="Site %s has occupancy 0.0 and contributes "
                                     "nothing to the model (a dummy atom?)." % s.label,
                             plat=["PLAT074"], based_on=["_atom_site_occupancy"],
                             sites=[s.label]))

    # Partial occupancy without any disorder bookkeeping.
    partial = [s for s in blk.sites
               if s.occupancy is not None and tol < s.occupancy < 1.0 - tol]
    has_groups = any(s.disorder_group not in (None, "", ".", "?")
                     for s in blk.sites)
    has_asm = any(s.disorder_assembly not in (None, "", ".", "?")
                  for s in blk.sites)
    if partial and not (has_groups or has_asm):
        f.append(Finding(
            code="PARTIAL_OCC_NO_DISORDER_GROUP", severity=INFO,
            message="%d site(s) have partial occupancy but no "
                    "_atom_site_disorder_group or _atom_site_disorder_assembly is "
                    "given, so the disorder model cannot be checked for closure."
                    % len(partial),
            plat=["PLAT299", "PLAT300"],
            based_on=["_atom_site_occupancy", "_atom_site_disorder_group",
                      "_atom_site_disorder_assembly"],
            sites=[s.label for s in partial[:12]]))
        nt.append(NotTested(
            area="disorder group occupancy closure",
            reason="no disorder group/assembly labels to sum over",
            plat=["PLAT077"]))
    elif has_groups:
        by = defaultdict(list)
        for s in blk.sites:
            g = s.disorder_group
            if g in (None, "", ".", "?") or s.occupancy is None:
                continue
            by[(s.disorder_assembly or "-", g)].append(s)
        full = 1.0 - ctx["occ_tol"]
        sums = defaultdict(float)
        negative = set()
        for (asm, g), members in by.items():
            if str(g).strip().startswith("-"):
                # SHELX negative PART: the group is disordered across a
                # symmetry element and its other half is its own image, which
                # is not in the atom list. A methyl on a mirror at PART -1 and
                # occupancy 0.5 is complete (COD 7110006), so the assembly
                # cannot be summed from the listed groups.
                negative.add(asm)
                continue
            # A group's occupancy comes from its PARTIAL members. Full-weight
            # atoms carrying a group tag are bookkeeping: COD 4123051 tags two
            # occupancy-1 H with group 1, and taking the largest member gave
            # 1 + 0.25 instead of the real 0.75 + 0.25.
            partial = [m.occupancy for m in members if m.occupancy < full]
            sums[asm] += max(partial) if partial else max(
                m.occupancy for m in members)
        if negative:
            nt.append(NotTested(
                area="occupancy closure for assembl%s %s"
                     % ("y" if len(negative) == 1 else "ies",
                        ", ".join(sorted(negative))[:120]),
                reason="a negative-PART group's other half is its own symmetry "
                       "image, which is not in the atom list",
                plat=["PLAT077"]))
        for asm, total in sums.items():
            if asm in negative:
                continue
            if abs(total - 1.0) > ctx["disorder_sum_tol"]:
                f.append(Finding(
                    code="DISORDER_SUM_NOT_ONE", severity=WARN,
                    message="Disorder assembly %r: the group occupancies sum to "
                            "%.3f rather than 1.0." % (asm, total),
                    plat=["PLAT077"],
                    based_on=["_atom_site_occupancy", "_atom_site_disorder_group",
                              "_atom_site_disorder_assembly"],
                    evidence={"assembly": asm, "sum": round(total, 4)}))

    nt.append(NotTested(
        area="occupancy < 1.0 for an atom on a special position",
        reason="requires site-symmetry analysis at each site, which is not "
               "implemented",
        plat=["PLAT076"]))
    return (f, nt)
