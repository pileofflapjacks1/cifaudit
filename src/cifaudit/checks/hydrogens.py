"""Hydrogen atoms: presence, connectivity and X-H distances."""

from typing import List, Tuple

from ..elements import radius
from ..expand import ExpandError, expand
from ..geom import min_distance
from ..model import ERROR, Finding, INFO, NotTested, WARN

NAME = "hydrogens"

# Target X-H distances for X-ray data, from the PLATON test descriptions
# (PLAT350/352/354): "XRAY: 0.96" for C-H, 0.87 for N-H, 0.82 for O-H.
XH_TARGET_XRAY = {"C": 0.96, "N": 0.87, "O": 0.82}
_MISSING = (None, "", ".", "?")


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    h_sites = [s for s in blk.sites if s.is_hydrogen]
    has_c = any((s.element or "") == "C" for s in blk.sites)

    if not h_sites:
        if has_c:
            f.append(Finding(
                code="NO_HYDROGENS_WITH_CARBON", severity=WARN,
                message="No H atoms are present although carbon is. This is "
                        "unusual for a carbon-containing compound, but may be "
                        "correct.",
                plat=["PLAT040"],
                based_on=["_atom_site_type_symbol", "_atom_site_label"]))
        if blk.h_treatment not in _MISSING and \
                str(blk.h_treatment).lower() in ("geom", "mixed", "difmap", "calc"):
            f.append(Finding(
                code="H_TREATMENT_WITHOUT_H", severity=ERROR,
                message="_atom_sites_solution_hydrogens is %r but the CIF contains "
                        "no H atoms." % blk.h_treatment,
                plat=["PLAT078", "PLAT079"],
                based_on=["_atom_sites_solution_hydrogens"],
                evidence={"h_treatment": blk.h_treatment}))
        nt.append(NotTested(area="X-H distances and H connectivity",
                            reason="no H atoms in the model",
                            plat=["PLAT305", "PLAT303", "PLAT350"]))
        return (f, nt)

    if not blk.cell.complete():
        nt.append(NotTested(area="X-H distances and H connectivity",
                            reason="the unit cell is incomplete",
                            plat=["PLAT305", "PLAT303", "PLAT350"]))
        return (f, nt)
    try:
        sites, cell = expand(blk)
    except ExpandError as exc:
        nt.append(NotTested(area="X-H distances and H connectivity",
                            reason=str(exc),
                            plat=["PLAT305", "PLAT303", "PLAT350"]))
        return (f, nt)

    heavy = [s for s in sites
             if s.element and s.element not in ("H", "D")
             and all(v is not None for v in s.fract)]
    hs = [s for s in sites if s.element in ("H", "D")
          and all(v is not None for v in s.fract)]

    for h in hs:
        partners = []
        for x in heavy:
            rx = radius(x.element)
            if rx is None:
                continue
            d, _ = min_distance(cell, h.fract, x.fract)
            if d <= (rx + radius("H")) * ctx["bond_factor"]:
                partners.append((x, d))
        if not partners:
            f.append(Finding(
                code="ISOLATED_H", severity=WARN,
                message="H atom %s has no heavy-atom neighbour within bonding "
                        "distance." % h.label,
                plat=["PLAT305"],
                based_on=["_atom_site_fract_x", "_atom_site_fract_y",
                          "_atom_site_fract_z"],
                sites=[h.label]))
            continue
        if len(partners) > 1:
            f.append(Finding(
                code="H_MULTIPLE_CONNECTIONS", severity=WARN,
                message="H atom %s is within bonding distance of %d heavy atoms "
                        "(%s)." % (h.label, len(partners),
                                   ", ".join("%s %.2f A" % (p.label, d)
                                             for p, d in partners[:4])),
                plat=["PLAT303"],
                based_on=["_atom_site_fract_x"],
                sites=[h.label] + [p.label for p, _ in partners[:4]],
                evidence={"n_partners": len(partners)}))
        x, d = min(partners, key=lambda t: t[1])
        target = XH_TARGET_XRAY.get(x.element or "")
        if target is not None and d < target - ctx["xh_short_tol"]:
            f.append(Finding(
                code="SHORT_XH", severity=INFO,
                message="%s-%s is %.3f A, %.3f A shorter than the X-ray target of "
                        "%.2f A for %s-H." % (x.label, h.label, d, target - d,
                                              target, x.element),
                plat=["PLAT350" if x.element == "C" else
                      ("PLAT352" if x.element == "N" else "PLAT354")],
                based_on=["_atom_site_fract_x"], sites=[x.label, h.label],
                evidence={"distance": round(d, 4), "target_xray": target}))

    if blk.h_treatment in _MISSING:
        nt.append(NotTested(area="hydrogen treatment consistency",
                            reason="_atom_sites_solution_hydrogens is absent",
                            plat=["PLAT078", "PLAT079"]))
    nt.append(NotTested(area="unrefined donor-H reporting, refined C-H",
                        reason="requires the refinement instruction file",
                        plat=["PLAT007", "PLAT164"]))
    return (f, nt)
