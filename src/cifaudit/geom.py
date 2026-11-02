"""Geometry helpers. Periodic distances by explicit image search.

Images are enumerated over -1..+1 in each axis rather than using a minimum-image
convention, because minimum-image is wrong for skewed cells.
"""

import itertools
from typing import List, Optional, Tuple

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None

try:
    import numpy as _np
except ImportError:  # pragma: no cover
    _np = None

IMAGES = [t for t in itertools.product((-1, 0, 1), repeat=3)]
# Furthest an atom may sit from its own symmetry image and still count as part
# of a moiety disordered across that symmetry element (see exclusive_pairs).
SELF_IMAGE_MAX = 2.6
# Pairs further apart than this are never tested for exclusivity: no check
# looks further (the widest is the isolated-atom contact, at most
# 2 x 2.60 + 0.90 = 6.1 A for the largest tabulated radius).
EXCLUSIVITY_RANGE = 6.5
# Cordero covalent radius of H, and the margin added to r_X + r_H when that
# exceeds the flat X-H maximum (see xh_reach).
H_RADIUS = 0.31
XH_EXTRA = 0.15


def xh_reach(r_x, xh_max):
    """Longest X-H distance at which an H may be bonded to an atom of radius r_x.

    The flat `xh_max` (1.35 A) suits C, N and O, but not heavier atoms: P-H runs
    1.35-1.44 A and a bridging metal hydride 1.7-2.1 A. With a flat cutoff the
    third precision round read P-H at 1.351 A (COD 2014097) and Lu-H at 1.97 A
    (COD 4123099) as unbonded H, producing false isolated-H and short-contact
    findings. For heavier atoms the reach is r_X + r_H + 0.15 A - Lu 2.33, Pd
    1.85, P 1.53 - and never less than `xh_max`.
    """
    if r_x is None:
        return xh_max
    return max(xh_max, r_x + H_RADIUS + XH_EXTRA)


def min_distance(cell, f1, f2) -> Tuple[float, Tuple[int, int, int]]:
    """Shortest distance (Angstrom) between two fractional sites over images."""
    best = None
    best_img = (0, 0, 0)
    for img in IMAGES:
        p1 = cell.orthogonalize(gemmi.Fractional(f1[0], f1[1], f1[2]))
        p2 = cell.orthogonalize(gemmi.Fractional(f2[0] + img[0],
                                                 f2[1] + img[1],
                                                 f2[2] + img[2]))
        d = p1.dist(p2)
        if best is None or d < best:
            best = d
            best_img = img
    return (best, best_img)


def self_image_distance(cell, f) -> Optional[float]:
    """Shortest distance from a site to its own non-identity lattice images."""
    best = None
    for img in IMAGES:
        if img == (0, 0, 0):
            continue
        p1 = cell.orthogonalize(gemmi.Fractional(f[0], f[1], f[2]))
        p2 = cell.orthogonalize(gemmi.Fractional(f[0] + img[0], f[1] + img[1],
                                                 f[2] + img[2]))
        d = p1.dist(p2)
        if best is None or d < best:
            best = d
    return best


def _cartesian(cell, sites):
    """Fractional -> Cartesian for every site, once."""
    out = []
    for s in sites:
        p = cell.orthogonalize(gemmi.Fractional(*s.fract))
        out.append((p.x, p.y, p.z))
    return out


def _image_vectors(cell):
    """The 27 lattice translations -1..+1 as Cartesian vectors."""
    origin = cell.orthogonalize(gemmi.Fractional(0.0, 0.0, 0.0))
    vecs = []
    for img in IMAGES:
        p = cell.orthogonalize(gemmi.Fractional(float(img[0]), float(img[1]),
                                                float(img[2])))
        vecs.append((p.x - origin.x, p.y - origin.y, p.z - origin.z))
    return vecs


def _pairwise_min(cell, sites):
    """Shortest distance and image for every pair, vectorised when numpy exists.

    Identical results to the per-pair loop: the same 27 explicit images, in the
    same order, with ties resolved to the first image. Minimum-image convention
    is NOT used, because it is wrong for skewed cells.
    """
    n = len(sites)
    cart = _cartesian(cell, sites)
    vecs = _image_vectors(cell)
    dists = {}
    if _np is not None and n > 1:
        P = _np.asarray(cart, dtype=float)              # (n, 3)
        V = _np.asarray(vecs, dtype=float)              # (27, 3)
        for i in range(n - 1):
            # (n-i-1, 27, 3): every later atom, shifted by every image
            diff = (P[i + 1:][:, None, :] + V[None, :, :]) - P[i][None, None, :]
            d2 = (diff * diff).sum(axis=2)              # (n-i-1, 27)
            k = d2.argmin(axis=1)
            dmin = _np.sqrt(d2[_np.arange(d2.shape[0]), k])
            for off, (dd, kk) in enumerate(zip(dmin, k)):
                dists[(i, i + 1 + off)] = (float(dd), IMAGES[int(kk)])
        return dists
    for i in range(n):
        for j in range(i + 1, n):
            dists[(i, j)] = min_distance(cell, sites[i].fract, sites[j].fract)
    return dists


def exclusive_pairs(sites, dists, max_sep, occ_tol, radius_fn=None,
                    bond_fraction=0.5, bond_tolerance=0.45,
                    self_image_max=SELF_IMAGE_MAX, xh_max=1.35):
    """Pairs of sites that are mutually exclusive ALTERNATIVES, not neighbours.

    Two atoms can share one crystallographic position with partial occupancies -
    COD 1000241 has Sb1 and W1 both at (0,0,0) on Wyckoff 16c with occupancy 0.5
    each. In any given unit cell only one is present. Treating them as two
    co-existing atoms makes every such site read as atoms 0.000 A apart, which
    was 27% of 500 deposited structures.

    Two criteria, either sufficient:

    * DECLARED - different `_atom_site_disorder_group` within the same
      `_atom_site_disorder_assembly`. 274 of those 500 files declare groups.
    * GEOMETRIC - within `max_sep` of each other, AND closer than any bond
      between those two elements could be, AND occupancies that cannot both be
      fully present (sum <= 1). This covers substitution at one site (separation
      0) and split positions a few tenths of an Angstrom apart.

      The bond guard is essential. Without it, a disordered C-H pair at 0.930 A
      with both occupancies 0.5 sums to exactly 1.0 and was classified as the
      hydrogen being an "alternative" to the carbon it is bonded to - which
      orphaned 55 more files with ISOLATED_H than it fixed. Alternatives occupy
      the SAME position; they are never a bond length apart.

    * ACROSS A SYMMETRY ELEMENT - see _symmetry_disorder_pairs.

    Full-occupancy atoms on top of each other are NOT excluded: that is a real
    error and must still be reported.
    """
    excl = set()
    _missing = ("", ".", "?")
    for (i, j), (d, _img) in dists.items():
        a, b = sites[i], sites[j]
        ga = (getattr(a, "disorder_group", None) or "").strip()
        gb = (getattr(b, "disorder_group", None) or "").strip()
        if ga not in _missing and gb not in _missing and ga != gb:
            aa = (getattr(a, "disorder_assembly", None) or "").strip()
            ab = (getattr(b, "disorder_assembly", None) or "").strip()
            if aa == ab:
                excl.add((i, j))
                continue
            if (_positive_group(ga) and _positive_group(gb)
                    and a.occ + b.occ <= 1.0 + occ_tol
                    and d <= self_image_max):
                # Different groups under DIFFERENT assembly letters (or one
                # letter missing). Many CIFs file the two halves of one
                # disorder that way, so the letters cannot be trusted alone:
                # COD 7110032 has H25A (A, group 1) 0.57 A from C25' (no
                # letter, group 2), the alternative position of its own C25.
                # A first version kept pairs at a bond length bonded, and in
                # COD 4502089 an H of ring orientation A then bonded to a C of
                # orientation B at 0.94 A. Atoms of two INDEPENDENT disorders,
                # both partial, in different groups and bonded to each other
                # are rare; mislettered halves of one disorder are common.
                excl.add((i, j))
                continue
        if d > max_sep or (a.occ + b.occ) > 1.0 + occ_tol:
            continue
        if radius_fn is not None:
            ra, rb = radius_fn(a.element), radius_fn(b.element)
            if ra is not None and rb is not None and d > bond_fraction * (ra + rb):
                continue        # far enough apart to be a bond, not one site
        excl.add((i, j))
    if radius_fn is not None:
        excl |= _symmetry_disorder_pairs(sites, dists, excl, occ_tol, radius_fn,
                                         bond_tolerance, self_image_max, xh_max)
    return excl


def _positive_group(g):
    g = (g or "").strip()
    return bool(g) and not g.startswith("-") and g not in ("0", ".", "?")


def _negative_part(site):
    return (getattr(site, "disorder_group", None) or "").strip().startswith("-")


def _group_key(site):
    return ((getattr(site, "disorder_assembly", None) or "").strip(),
            (getattr(site, "disorder_group", None) or "").strip())


def _symmetry_disorder_pairs(sites, dists, already, occ_tol, radius_fn,
                             bond_tolerance, self_image_max, xh_max=1.35):
    """Pairs that are alternatives because a moiety is disordered ACROSS a
    symmetry element.

    A partially occupied group lying on a mirror, 2-fold axis or inversion
    centre overlaps its own symmetry image: in any one unit cell either the
    group or its image is present, never both. SHELX writes this as a negative
    PART ("do not connect to symmetry equivalents"); many CIFs only show it as
    occupancy 0.5 or lower. In the second precision round this was the single
    largest cause of false findings - an H 0.385 A from its own 2-fold image
    (COD 7041061), a C 0.75 A from its own mirror image (COD 4310047) - all
    reported as clashes, short contacts or extra connections.

    Sites a (copy A) and b (copy B) are alternatives when they lie in different
    symmetry copies (no operator in common), both are partially occupied with
    occupancies summing to at most one, and one of them, X, sits close to its OWN
    image X' in the other's copy - closer than X and X' could both be present -
    while the other IS X' or is bonded to X' (see alternative_via). "Close" is:

    * X in a negative PART: within `self_image_max` (SHELX declared it).
    * X is H: within 2.0 A, the H...H minimum.
    * otherwise: shorter than a bond could be (under 0.85 x 2r), or between a
      bond and `self_image_max` - too long to bond, too short to coexist.

    For a negative PART the test is made for the whole group: two atoms of one
    negative-PART group (same assembly and group) in different copies are
    alternatives when ANY atom of that group lies within `self_image_max` of its
    own image between those two copies. A ring disordered across a 2-fold axis
    has atoms near the axis close to their images, while H16A sits 0.385 A from
    the image of H18A, its counterpart in the other orientation, and far from
    its own (COD 7041061).

    A self-image at a bond length with no negative PART is left alone: a
    half-occupied ordered molecule can be bonded to its own image across a
    centre of symmetry. Atoms ON the symmetry element are in both copies, so
    they keep their bonds to either.
    """
    copy_index = {}
    for k, s in enumerate(sites):
        for o in getattr(s, "ops", ()) or ():
            copy_index[(s.label, o)] = k
    if not copy_index:
        return set()

    def qualifies(x, d_s):
        if _negative_part(x):
            return d_s <= self_image_max
        if (x.element or "").upper() in ("H", "D"):
            return d_s < 2.0
        r = radius_fn(x.element)
        if r is None:
            return False
        return d_s < 0.85 * 2 * r or 2 * r + bond_tolerance < d_s <= self_image_max

    groups = {}
    for s in sites:
        if _negative_part(s):
            groups.setdefault(_group_key(s), set()).add(s.label)
    overlap_cache = {}

    def group_overlaps(key, ka, kb):
        """Does any atom of negative group `key` sit near its own image
        between copies ka and kb?"""
        ck = (key, ka, kb)
        if ck not in overlap_cache:
            hit = False
            for lab in groups.get(key, ()):
                p, q = copy_index.get((lab, ka)), copy_index.get((lab, kb))
                if p is None or q is None or p == q:
                    continue
                if dists[(min(p, q), max(p, q))][0] <= self_image_max:
                    hit = True
                    break
            overlap_cache[ck] = hit
        return overlap_cache[ck]

    full = 1.0 - occ_tol
    ops = [frozenset(getattr(x, "ops", ()) or ()) for x in sites]
    neg = [_negative_part(x) for x in sites]
    gkey = [_group_key(x) for x in sites]
    self_cache = {}

    def pdist(p, q):
        return dists[(min(p, q), max(p, q))][0]

    def bond_reach(p, q):
        """Largest distance at which sites p and q could be bonded."""
        hp = (sites[p].element or "").upper() in ("H", "D")
        hq = (sites[q].element or "").upper() in ("H", "D")
        if hp and hq:
            return 0.0
        if hp or hq:
            return xh_reach(radius_fn(sites[q if hp else p].element), xh_max)
        rp, rq = radius_fn(sites[p].element), radius_fn(sites[q].element)
        if rp is None or rq is None:
            return 0.0
        return rp + rq + bond_tolerance

    def self_images(x_idx, other_ops):
        """Images of site x_idx, in the copies listed, that it cannot coexist
        with."""
        ck = (x_idx, other_ops)
        if ck not in self_cache:
            x = sites[x_idx]
            hits = []
            for m in other_ops:
                k = copy_index.get((x.label, m))
                if k is None or k == x_idx or k in hits:
                    continue
                if qualifies(x, pdist(x_idx, k)):
                    hits.append(k)
            self_cache[ck] = hits
        return self_cache[ck]

    def alternative_via(x_idx, y_idx, other_ops):
        """X's own image X' excludes X; does Y go with X' rather than with X?

        Y is an alternative to X when it IS X', or is bonded to X' and nearer
        to X' than to X. A water coordinating a half-occupied Ca is not an
        alternative to that Ca just because the Ca's image is too close
        (COD 9005067); an H bonded to the image of a disordered C is
        (COD 4310047, 1.305 A from C3, 0.96 A from C3').
        """
        for k in self_images(x_idx, other_ops):
            if k == y_idx:
                return True
            dy = pdist(y_idx, k)
            if dy <= bond_reach(y_idx, k) and dy < pdist(y_idx, x_idx):
                return True
            # Two bonds out: Y and X' share a partially occupied neighbour Z
            # that is not also bonded to X (an atom bonded to both alternatives
            # is shared, not part of either). In a methanol disordered across
            # an inversion centre (COD 4342081) H15B clashes with H15C', which
            # is bonded to C15' - the image of C15 - and the two methyl H are
            # geminal, so neither is X' or bonded to it. The same shape occurs
            # for a CH2Cl2 about a 2-fold axis (COD 4102011).
            for z in partial_nb.get(k, ()):
                if (z != y_idx and z != x_idx and z in partial_nb.get(y_idx, ())
                        and z not in partial_nb.get(x_idx, ())):
                    return True
        return False

    # Bonded neighbours among partially occupied sites, for alternative_via's
    # two-bond test. Built once; only partial sites can be part of a moiety
    # that alternates with its image.
    partial_nb = {}
    for (i, j), (d, _img) in dists.items():
        if d > 3.5 or sites[i].occ >= full or sites[j].occ >= full:
            continue
        if d <= bond_reach(i, j):
            partial_nb.setdefault(i, set()).add(j)
            partial_nb.setdefault(j, set()).add(i)

    out = set()
    for (i, j), (d, _img) in dists.items():
        if d > EXCLUSIVITY_RANGE or (i, j) in already:
            continue
        a, b = sites[i], sites[j]
        if a.occ >= full or b.occ >= full or a.occ + b.occ > 1.0 + occ_tol:
            continue
        oa, ob = ops[i], ops[j]
        if not oa or not ob or oa & ob:
            continue
        if (neg[i] and neg[j] and gkey[i] == gkey[j]
                and group_overlaps(gkey[i], min(oa), min(ob))):
            out.add((i, j))
            continue
        if alternative_via(i, j, ob) or alternative_via(j, i, oa):
            out.add((i, j))
    return out


def bond_graph(sites, cell, radius_fn, tolerance, xh_max,
               shared_site_max=0.0, occ_tol=1e-3, bond_fraction=0.5):
    """Pairwise shortest distances, the covalent bond graph, and H candidates.

    Heavy-atom bonding uses an ADDITIVE tolerance, which is the convention the
    field actually uses:

        d(i,j) <= r_i + r_j + epsilon

    epsilon is "usually ranging between 0.35 A and 0.45 A" (arXiv 2601.02017,
    bonding radii derived from 519,021 COD records); Jmol-style/MaterialsCoord
    use 0.45 A, CCDC's recommended figure is about 0.4 A, and CoRE MOF 2019 uses
    Cordero radii plus 0.2 A for metals under 1.95 A and 0.4 A above. cifaudit
    previously used a MULTIPLICATIVE 1.2x, which is not a convention anyone else
    follows.

    Hydrogen is handled separately and deliberately (the reach grows with the
    heavy atom's radius beyond C/N/O; see xh_reach). The COD-derived radii set
    above has no hydrogen value at all and discards every contact involving
    hydrogen, because radius-based bonding does not work for H: the standard
    epsilon would make H-O bonding MORE generous (1.42 A) than the old 1.2x rule
    (1.16 A), not less. So an H is bonded to exactly ONE heavy atom - its
    nearest within `xh_max` - which matches PLATON PLAT303: "Hydrogen atoms are
    generally connected to only one other atom."

    H...H is never treated as a bond. Dihydrogen complexes exist and are not
    modelled.

    Mutually exclusive alternatives (see exclusive_pairs) never bond to one
    another, and the set is returned so callers apply the same definition.

    Returns (dists, adj, h_candidates, exclusive). h_candidates maps an H index to
    every heavy atom within xh_max as (index, distance) - so PLAT303's
    "more than one connection" case can still be reported even though only the
    nearest is bonded.
    """
    dists = _pairwise_min(cell, sites)
    n = len(sites)
    adj = {i: set() for i in range(n)}
    radii = [radius_fn(s.element) for s in sites]
    is_h = [(s.element or "").upper() in ("H", "D") for s in sites]
    h_candidates = {}
    excl = exclusive_pairs(sites, dists, shared_site_max, occ_tol,
                           radius_fn, bond_fraction, tolerance,
                           xh_max=xh_max)

    for (i, j), (d, _img) in dists.items():
        if (i, j) in excl:
            continue        # mutually exclusive alternatives never bond
        if is_h[i] and is_h[j]:
            continue
        if is_h[i] or is_h[j]:
            h_idx, x_idx = (i, j) if is_h[i] else (j, i)
            if radii[x_idx] is not None and d <= xh_reach(radii[x_idx], xh_max):
                h_candidates.setdefault(h_idx, []).append((x_idx, d))
            continue
        ri, rj = radii[i], radii[j]
        if ri is not None and rj is not None and d <= ri + rj + tolerance:
            adj[i].add(j)
            adj[j].add(i)

    for h_idx, cands in h_candidates.items():
        x_idx, _d = min(cands, key=lambda t: t[1])
        adj[h_idx].add(x_idx)
        adj[x_idx].add(h_idx)

    return (dists, adj, h_candidates, excl)


def within_bonds(adj, max_depth):
    """For each atom, the set of atoms reachable within `max_depth` bonds.

    Contacts between atoms this close in the bond graph are a consequence of
    molecular geometry, not packing. PLATON separates intra- from
    inter-molecular contacts for the same reason. Excluding only 1,2 and 1,3
    was not enough: a 1,4 C...H contact across a planar aromatic ring runs
    about 2.1-2.2 A and was being reported as a clash (COD 2100301).
    """
    out = {}
    for start in adj:
        seen = {start}
        frontier = [start]
        for _ in range(max_depth):
            nxt = []
            for k in frontier:
                for n in adj[k]:
                    if n not in seen:
                        seen.add(n)
                        nxt.append(n)
            frontier = nxt
            if not frontier:
                break
        seen.discard(start)
        out[start] = seen
    return out
