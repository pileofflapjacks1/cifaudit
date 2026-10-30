"""Connectivity: isolated atoms and disconnected fragments."""

from typing import Dict, List, Set, Tuple

from ..elements import METALS, radius
from ..expand import ExpandError, expand
from ..geom import bond_graph
from ..model import Finding, INFO, NotTested, WARN

NAME = "fragments"
MAX_SITES = 400


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    if not blk.cell.complete():
        nt.append(NotTested(area="connectivity and isolated atoms",
                            reason="the unit cell is incomplete",
                            plat=["PLAT305", "PLAT004"]))
        return (f, nt)
    try:
        sites, cell = expand(blk)
    except ExpandError as exc:
        nt.append(NotTested(area="connectivity and isolated atoms",
                            reason=str(exc), plat=["PLAT305", "PLAT004"]))
        return (f, nt)

    sites = [s for s in sites if s.element and all(v is not None for v in s.fract)]
    if not sites:
        nt.append(NotTested(area="connectivity", reason="no usable atom sites",
                            plat=["PLAT305"]))
        return (f, nt)
    if len(sites) > MAX_SITES:
        nt.append(NotTested(area="connectivity and isolated atoms",
                            reason="%d unit-cell sites exceeds the %d-site cap"
                                   % (len(sites), MAX_SITES),
                            plat=["PLAT305", "PLAT004"]))
        return (f, nt)

    no_radius = sorted({s.element for s in sites if radius(s.element) is None})
    if no_radius:
        nt.append(NotTested(
            area="connectivity involving %s" % ", ".join(no_radius),
            reason="no covalent radius tabulated, so bonds to these atoms cannot "
                   "be decided",
            plat=["PLAT305"]))

    _dists, adj = bond_graph(sites, cell, ctx["bond_factor"], radius)

    seen: Set[int] = set()
    comps: List[List[int]] = []
    for i in range(len(sites)):
        if i in seen:
            continue
        stack, comp = [i], []
        seen.add(i)
        while stack:
            k = stack.pop()
            comp.append(k)
            for n in adj[k]:
                if n not in seen:
                    seen.add(n)
                    stack.append(n)
        comps.append(comp)

    singles = [c[0] for c in comps if len(c) == 1]
    for idx in singles:
        s = sites[idx]
        if radius(s.element) is None:
            continue
        if s.element in ("H", "D"):
            continue        # reported by the hydrogens check as ISOLATED_H
        if s.element in METALS:
            sev, hint = WARN, "for a metal this often means missing coordination"
        elif s.element in ("O", "N"):
            sev, hint = WARN, "for O or N this often means missing H"
        else:
            sev, hint = INFO, "it may be a counter-ion or free atom"
        f.append(Finding(
            code="ISOLATED_ATOM", severity=sev,
            message="%s (%s) has no bonded neighbour anywhere in the unit cell; "
                    "%s." % (s.label, s.element, hint),
            plat=["PLAT305"],
            based_on=["_atom_site_fract_x", "_atom_site_fract_y",
                      "_atom_site_fract_z"],
            sites=[s.label], evidence={"element": s.element}))

    if len(comps) > 1:
        sizes = sorted((len(c) for c in comps), reverse=True)
        f.append(Finding(
            code="MULTIPLE_FRAGMENTS", severity=INFO,
            message="The unit cell contains %d disconnected fragments (sizes %s). "
                    "Legitimate for molecular or ionic structures; suspicious for "
                    "a framework." % (len(comps), sizes[:10]),
            plat=["PLAT004"],
            based_on=["_atom_site_fract_x"],
            evidence={"n_fragments": len(comps), "sizes": sizes[:20]}))

    nt.append(NotTested(
        area="polymeric network dimensionality (0D/1D/2D/3D)",
        reason="components are found within one unit cell; distinguishing a "
               "1D chain from a 3D net requires propagation across cells, which "
               "is not implemented",
        plat=["PLAT004"]))
    return (f, nt)
