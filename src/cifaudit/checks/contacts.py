"""Interatomic distances: impossibly close atoms and short non-bonded contacts."""

from typing import List, Tuple

from ..elements import radius
from ..expand import ExpandError, expand
from ..geom import bond_graph, self_image_distance, within_bonds
from ..model import ERROR, Finding, NotTested, WARN

NAME = "contacts"
MAX_PAIR_SITES = 400

# Elements that routinely act as hydrogen-bond donors and acceptors. An
# H...acceptor approach where the H is bonded to one of these is a hydrogen
# bond, which is a normal feature of a correct structure - not a clash.
# Measured on 16 deposited COD structures, treating these as short contacts
# produced 18 findings, ALL of them ordinary hydrogen bonds.
HBOND_ELEMENTS = frozenset(("N", "O", "F", "S", "Cl", "Br"))
# Carbon is a weak but entirely ordinary hydrogen-bond donor: C-H...O contacts
# around 2.2 A are common, and treating them as clashes produced a false
# positive on COD 2100301 (neutron, 15 K).
WEAK_DONORS = frozenset(("C", "N", "O", "F", "S", "Cl", "Br", "B"))


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
    near = within_bonds(adj, ctx["exclude_within_bonds"])

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

    # Pass A: find hydrogen bonds, so their secondary geometry can be excluded.
    # An N-H...O=C bond necessarily brings the H within ~2.2 A of the carboxyl
    # CARBON as well; reporting that as a separate clash is double counting
    # (observed on COD 2100301, N-H...C8 at 2.176 A).
    hbond_secondary = {}
    for (i, j), (d, _img) in dists.items():
        for h_idx, acc_idx in ((i, j), (j, i)):
            if sites[h_idx].element not in ("H", "D"):
                continue
            if sites[acc_idx].element not in HBOND_ELEMENTS:
                continue
            if d > ctx["nonbonded_min"] + 0.6:
                continue
            if not ({sites[k].element for k in adj[h_idx]} & WEAK_DONORS):
                continue
            hbond_secondary.setdefault(h_idx, set()).update(adj[acc_idx])

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
        if j in near.get(i, ()):
            continue        # within exclude_within_bonds of each other (1,2/1,3/1,4)
        if j in hbond_secondary.get(i, ()) or i in hbond_secondary.get(j, ()):
            continue        # geometry implied by a hydrogen bond, not a clash
        ra, rb = radius(a.element), radius(b.element)
        if ra is None or rb is None:
            continue
        both_h = a.element in ("H", "D") and b.element in ("H", "D")

        # Hydrogen bond? Exactly one H, bonded to a donor, approaching an acceptor.
        if not both_h:
            h_i, other = (i, b) if a.element in ("H", "D") else (
                (j, a) if b.element in ("H", "D") else (None, None))
            if h_i is not None and other.element in HBOND_ELEMENTS:
                donors = {sites[k].element for k in adj[h_i]}
                if donors & WEAK_DONORS:
                    if d < ctx["hbond_min"]:
                        f.append(Finding(
                            code="HBOND_IMPLAUSIBLY_SHORT", severity=WARN,
                            message="%s...%s is %.3f A. The geometry is a hydrogen "
                                    "bond, but %.2f A is shorter than any plausible "
                                    "one - check the H placement."
                                    % (a.label, b.label, d, d),
                            plat=["PLAT430"],
                            based_on=["_atom_site_fract_x"],
                            sites=[a.label, b.label],
                            evidence={"distance": round(d, 4),
                                      "hbond_min": ctx["hbond_min"]}))
                    continue        # a normal hydrogen bond: not a finding

        limit = ctx["hh_min"] if both_h else ctx["nonbonded_min"]
        if d >= limit:
            continue
        if ("short", pair) in reported:
            continue
        reported.add(("short", pair))
        f.append(Finding(
            code="SHORT_HH_CONTACT" if both_h else "SHORT_NONBONDED_CONTACT",
            severity=WARN,
            message="%s...%s is %.3f A, shorter than the %.2f A %s threshold, "
                    "and is not a bond, a hydrogen bond, or a pair within %d "
                    "bonds of each other."
                    % (a.label, b.label, d, limit,
                       "H...H" if both_h else "non-bonded",
                       ctx["exclude_within_bonds"]),
            plat=["PLAT410", "PLAT411"] if both_h else ["PLAT430", "PLAT432"],
            based_on=["_atom_site_fract_x", "_atom_site_fract_y",
                      "_atom_site_fract_z"],
            sites=[a.label, b.label],
            evidence={"distance": round(d, 4), "image": list(img),
                      "both_full_weight": a.occ > 0.95 and b.occ > 0.95,
                      "note": "pairs within %d bonds, and hydrogen bonds, "
                              "are excluded from this test"
                              % ctx["exclude_within_bonds"]}))

    nt.append(NotTested(
        area="hydrogen-bond geometry (D-H...A angles, donors without acceptors), "
             "solvent-accessible voids, CIF-reported-vs-computed distance "
             "consistency",
        reason="not implemented",
        plat=["PLAT420", "PLAT482", "PLAT601", "PLAT724"]))
    return (f, nt)
