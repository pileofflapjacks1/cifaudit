"""Connectivity: isolated atoms and disconnected fragments."""

from typing import Dict, List, Set, Tuple

from ..elements import COVALENT, METALS, radius
from ..expand import ExpandError, expand
from ..geom import bond_graph
from ..model import Finding, INFO, NotTested, WARN

NAME = "fragments"
MAX_RADIUS = max(COVALENT.values())
UNDECIDED = object()
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

    dists, adj, _h_cands, excl = bond_graph(
        sites, cell, radius, ctx["bond_tolerance"], ctx["xh_bond_max"],
        ctx["shared_site_max"], ctx["occ_tol"],
        ctx["shared_site_bond_fraction"])

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

    def nearest_contact(i):
        """Closest non-H neighbour within secondary-contact distance, or None.

        "Isolated" means nothing nearby, not merely no covalent bond. Long
        secondary bonds - Tl-N 2.74 A, Pb-S 2.96 A, axial Cu-Cl 2.90 A in the
        precision study - fall outside the covalent tolerance but are plainly
        coordination. Mutually exclusive alternatives do not count.
        """
        ri = radius(sites[i].element)
        best = None
        for (a, b), (d, _img) in dists.items():
            if i not in (a, b) or (a, b) in excl:
                continue
            j = b if a == i else a
            if sites[j].element in ("H", "D"):
                continue
            rj = radius(sites[j].element)
            if ri is None:
                continue
            if rj is None:
                # A neighbour whose radius is unknown could be bonded to this
                # atom, so "isolated" cannot be decided: abstain. Any element
                # within reach of the LARGEST tabulated radius counts.
                if d <= ri + MAX_RADIUS + ctx["contact_tolerance"]:
                    return UNDECIDED
                continue
            if d <= ri + rj + ctx["contact_tolerance"] and (best is None or d < best[0]):
                best = (d, j)
        return best

    reported = set()
    undecided: List[str] = []
    for idx in singles:
        s = sites[idx]
        if radius(s.element) is None:
            continue
        if s.element in ("H", "D"):
            continue        # reported by the hydrogens check as ISOLATED_H
        if s.label in reported:
            continue
        near = nearest_contact(idx)
        if near is UNDECIDED:
            undecided.append(s.label)
            continue
        if near is not None:
            continue        # a secondary contact: coordinated, not isolated
        reported.add(s.label)
        if s.element == "O":
            # Expert convention (J. R. Ramirez): an unbound single O is routinely
            # modelled solvent - water whose H were not located - or an O placed
            # to account for electron density that may belong to another
            # molecule. Information, never a warning.
            sev = INFO
            msg = ("%s (O) has no bonded neighbour. An isolated O is usually "
                   "modelled solvent: water whose H were not located, or O placed "
                   "to account for electron density." % s.label)
        elif s.element in METALS:
            sev = WARN
            msg = ("%s (%s) has no neighbour within coordination distance; for a "
                   "metal this usually means missing coordination."
                   % (s.label, s.element))
        elif s.element == "N":
            sev = WARN
            msg = ("%s (N) has no neighbour within coordination distance; often a "
                   "solvent or guest molecule whose H were not modelled."
                   % s.label)
        else:
            sev = INFO
            msg = ("%s (%s) has no neighbour within coordination distance; it may "
                   "be a counter-ion or free atom." % (s.label, s.element))
        f.append(Finding(
            code="ISOLATED_ATOM", severity=sev, message=msg,
            plat=["PLAT305"],
            based_on=["_atom_site_fract_x", "_atom_site_fract_y",
                      "_atom_site_fract_z"],
            sites=[s.label],
            evidence={"element": s.element,
                      "contact_limit": "r_i + r_j + %.2f A" % ctx["contact_tolerance"]}))

    if undecided:
        labs = sorted(set(undecided))
        nt.append(NotTested(
            area="isolated-atom test for %s" % ", ".join(labs[:10]),
            reason="a neighbour within reach has no tabulated covalent radius, so "
                   "whether these atoms are bonded cannot be decided",
            plat=["PLAT305"]))

    # Fragment count is NOT reported as a finding. Measured on 16 deposited COD
    # structures it was >1 for 14 of them (88%): molecular and ionic crystals
    # have many fragments by construction, so a finding here is pure noise. The
    # count is still useful context, so it is attached to the isolated-atom
    # findings and returned as a note instead.
    if len(comps) > 1:
        sizes = sorted((len(c) for c in comps), reverse=True)
        ctx.setdefault("_notes", []).append(
            "connectivity: %d fragments in the unit cell (sizes %s). Normal for "
            "molecular or ionic structures; unexpected for a framework."
            % (len(comps), sizes[:10]))

    nt.append(NotTested(
        area="polymeric network dimensionality (0D/1D/2D/3D)",
        reason="components are found within one unit cell; distinguishing a "
               "1D chain from a 3D net requires propagation across cells, which "
               "is not implemented",
        plat=["PLAT004"]))
    return (f, nt)
