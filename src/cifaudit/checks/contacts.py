"""Interatomic distances: impossibly close atoms and short non-bonded contacts."""

from typing import List, Tuple

from ..elements import radius
from ..expand import ExpandError, expand
from ..geom import bond_graph, self_image_distance
from ..model import ERROR, Finding, NotTested, WARN

NAME = "contacts"
MAX_PAIR_SITES = 400


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    if not blk.cell.complete():
        nt.append(NotTested(area="all interatomic distances",
                            reason="the unit cell is incomplete",
                            plat=["PLAT310", "PLAT410", "PLAT411", "PLAT432"]))
        return (f, nt)
    try:
        sites, cell = expand(blk)
    except ExpandError as exc:
        nt.append(NotTested(area="all interatomic distances",
                            reason=str(exc),
                            plat=["PLAT310", "PLAT410", "PLAT411", "PLAT432"]))
        return (f, nt)

    sites = [s for s in sites if s.element and all(v is not None for v in s.fract)]
    if len(sites) > MAX_PAIR_SITES:
        nt.append(NotTested(
            area="interatomic distances",
            reason="%d unit-cell sites exceeds the %d-site pairwise cap"
                   % (len(sites), MAX_PAIR_SITES),
            plat=["PLAT310", "PLAT432"]))
        return (f, nt)

    hard = ctx["hard_floor"]
    dists, adj = bond_graph(sites, cell, ctx["bond_factor"], radius)

    # A site sitting on its own lattice image means the cell or the coordinates
    # are wrong, and it is not caught by pairwise comparison.
    for i, a in enumerate(sites):
        d_self = self_image_distance(cell, a.fract)
        if d_self is not None and d_self < hard:
            f.append(Finding(
                code="SITE_ON_OWN_IMAGE", severity=ERROR,
                message="Site %s lies %.3f A from its own lattice image, below "
                        "the %.2f A floor. The cell or the coordinates are wrong."
                        % (a.label, d_self, hard),
                plat=["PLAT310"],
                based_on=["_atom_site_fract_x", "_cell_length_a"],
                sites=[a.label], evidence={"distance": round(d_self, 4)}))

    reported = set()
    for (i, j), (d, img) in sorted(dists.items(), key=lambda kv: kv[1][0]):
        a, b = sites[i], sites[j]
        pair = tuple(sorted((a.label, b.label)))
        if d < hard:
            if ("hard", pair) in reported:
                continue
            reported.add(("hard", pair))
            f.append(Finding(
                code="ATOMS_TOO_CLOSE", severity=ERROR,
                message="%s and %s are %.3f A apart, below the %.2f A floor at "
                        "which no two atoms can both be real."
                        % (a.label, b.label, d, hard),
                plat=["PLAT310"],
                based_on=["_atom_site_fract_x", "_atom_site_fract_y",
                          "_atom_site_fract_z"],
                sites=[a.label, b.label],
                evidence={"distance": round(d, 4), "image": list(img),
                          "occupancies": [a.occ, b.occ]}))
            continue
        if j in adj[i]:
            continue                        # 1,2 - bonded
        if adj[i] & adj[j]:
            continue                        # 1,3 - geminal, share a bonded atom
        ra, rb = radius(a.element), radius(b.element)
        if ra is None or rb is None:
            continue
        both_h = a.element in ("H", "D") and b.element in ("H", "D")
        limit = ctx["hh_min"] if both_h else ctx["nonbonded_min"]
        if d >= limit:
            continue
        if ("short", pair) in reported:
            continue
        reported.add(("short", pair))
        f.append(Finding(
            code="SHORT_HH_CONTACT" if both_h else "SHORT_NONBONDED_CONTACT",
            severity=WARN,
            message="%s...%s is %.3f A, shorter than the %.2f A %s threshold. "
                    "They are neither bonded nor geminal (1,3)."
                    % (a.label, b.label, d, limit,
                       "H...H" if both_h else "non-bonded"),
            plat=["PLAT410", "PLAT411"] if both_h else ["PLAT430", "PLAT432"],
            based_on=["_atom_site_fract_x", "_atom_site_fract_y",
                      "_atom_site_fract_z"],
            sites=[a.label, b.label],
            evidence={"distance": round(d, 4), "image": list(img),
                      "both_full_weight": a.occ > 0.95 and b.occ > 0.95,
                      "note": "1,2 and 1,3 pairs are excluded from this test"}))

    nt.append(NotTested(
        area="hydrogen-bond geometry (D-H...A angles, donors without acceptors), "
             "solvent-accessible voids, CIF-reported-vs-computed distance "
             "consistency",
        reason="not implemented",
        plat=["PLAT420", "PLAT482", "PLAT601", "PLAT724"]))
    return (f, nt)
