"""Hydrogen atoms: presence, connectivity and X-H distances."""

from typing import List, Tuple

from ..elements import radius
from ..expand import ExpandError, expand
from ..geom import bond_graph, min_distance
from ..model import ERROR, Finding, INFO, NotTested, WARN

NAME = "hydrogens"

# Target X-H distances for X-ray data, from the PLATON test descriptions
# (PLAT350/352/354): "XRAY: 0.96" for C-H, 0.87 for N-H, 0.82 for O-H.
XH_TARGET_XRAY = {"C": 0.96, "N": 0.87, "O": 0.82}
XH_TARGET_NEUTRON = {"C": 1.08, "N": 1.009, "O": 0.983}
_MISSING = (None, "", ".", "?")


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    h_sites = [s for s in blk.sites if s.is_hydrogen]
    has_c = any((s.element or "") == "C" for s in blk.sites)

    if not h_sites and blk.implicit_h > 0:
        nt.append(NotTested(
            area="X-H distances and H connectivity",
            reason="H atoms are declared only implicitly, through "
                   "_atom_site_attached_hydrogens, so they have no positions",
            plat=["PLAT305", "PLAT303", "PLAT350"]))
        return (f, nt)

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

    usable = [s for s in sites
              if s.element and all(x is not None for x in s.fract)]
    _dists, adj, h_cands, excl = bond_graph(
        usable, cell, radius, ctx["bond_tolerance"], ctx["xh_bond_max"],
        ctx["shared_site_max"], ctx["occ_tol"],
        ctx["shared_site_bond_fraction"])

    d_unknown_probe = set()
    for idx, s in enumerate(usable):
        if s.element not in ("H", "D"):
            continue
        cands = h_cands.get(idx, [])
        if not cands:
            f.append(Finding(
                code="ISOLATED_H", severity=WARN,
                message="H atom %s has no heavy-atom neighbour within bonding "
                        "reach (%.2f A for C, N and O; r_X + 0.46 A for heavier "
                        "atoms)." % (s.label, ctx["xh_bond_max"]),
                plat=["PLAT305"],
                based_on=["_atom_site_fract_x", "_atom_site_fract_y",
                          "_atom_site_fract_z"],
                sites=[s.label]))
            continue

        x_idx, d = min(cands, key=lambda t: t[1])
        x = usable[x_idx]

        # Two candidates that are alternatives for one another are one site,
        # not two connections. COD 1000241's H sat on the shared Sb/W position.
        def _alt(p, q):
            return (min(p, q), max(p, q)) in excl

        # Count connections within the flat X-H maximum, plus the bonded atom.
        # The longer reach for heavy atoms decides what an H with no light
        # neighbour is bonded to (a hydride); it must not turn every C-H near
        # a metal into a "multiple connection".
        counted = [(k, dd) for k, dd in cands
                   if dd <= ctx["xh_bond_max"] or k == x_idx]
        groups = []
        for k, dd in sorted(counted, key=lambda t: t[1]):
            for grp in groups:
                if all(_alt(k, m) for m, _ in grp):
                    grp.append((k, dd))
                    break
            else:
                groups.append([(k, dd)])

        if len(groups) > 1:
            # PLAT303: "Hydrogen atoms are generally connected to only one other
            # atom. A hydrogen atom between two oxygen atoms is a special case.
            # Investigate whether this hydrogen atom is better described with a
            # disorder model with two partially occupied sites."
            # So this is an invitation to look, not an error. The H is still
            # bonded to its nearest neighbour only.
            others = sorted(counted, key=lambda t: t[1])
            f.append(Finding(
                code="H_MULTIPLE_CONNECTIONS", severity=INFO,
                message="H atom %s has %d distinct heavy sites within %.2f A "
                        "(%s). It is bonded to the nearest. A hydrogen between "
                        "two electronegative atoms may be better modelled as "
                        "disorder over two partially occupied sites."
                        % (s.label, len(groups), ctx["xh_bond_max"],
                           ", ".join("%s %.2f A" % (usable[k].label, dd)
                                     for k, dd in others[:4])),
                plat=["PLAT303"],
                based_on=["_atom_site_fract_x"],
                sites=[s.label] + [usable[k].label for k, _ in others[:4]],
                evidence={"n_candidates": len(counted),
                          "n_distinct_sites": len(groups),
                          "bonded_to": x.label,
                          "distance": round(d, 4)}))

        # Neutron data locates nuclei, so X-H is genuinely longer. The PLATON
        # test descriptions give both: "XRAY: 0.96 NEUT 1.08" for C-H.
        probe = (blk.radiation_probe or "").lower()
        neutron = "neutron" in probe
        hname = "D" if s.element == "D" else "H"
        if hname == "D" and not probe:
            # A deuterated structure is usually a neutron study, but this CIF
            # does not say which probe was used, so neither target is safe.
            d_unknown_probe.add(s.label)
            continue
        table = XH_TARGET_NEUTRON if neutron else XH_TARGET_XRAY
        target = table.get(x.element or "")
        if target is not None and d < target - ctx["xh_short_tol"]:
            f.append(Finding(
                code="SHORT_XH", severity=INFO,
                message="%s-%s is %.3f A, %.3f A shorter than the %s target of "
                        "%.3f A for %s-%s." % (x.label, s.label, d, target - d,
                                               "neutron" if neutron else "X-ray",
                                               target, x.element, hname),
                plat=["PLAT350" if x.element == "C" else
                      ("PLAT352" if x.element == "N" else "PLAT354")],
                based_on=["_atom_site_fract_x"], sites=[x.label, s.label],
                evidence={"distance": round(d, 4), "target": target,
                          "probe": "neutron" if neutron else "x-ray"}))

    if d_unknown_probe:
        nt.append(NotTested(
            area="X-D distances for %s" % ", ".join(sorted(d_unknown_probe)[:6]),
            reason="deuterium positions, but no _diffrn_radiation_probe, so the "
                   "X-ray and neutron targets cannot be chosen between",
            plat=["PLAT350", "PLAT352", "PLAT354"]))

    if blk.h_treatment in _MISSING:
        nt.append(NotTested(area="hydrogen treatment consistency",
                            reason="_atom_sites_solution_hydrogens is absent",
                            plat=["PLAT078", "PLAT079"]))
    nt.append(NotTested(area="unrefined donor-H reporting, refined C-H",
                        reason="requires the refinement instruction file",
                        plat=["PLAT007", "PLAT164"]))
    return (f, nt)
