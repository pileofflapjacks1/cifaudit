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


def bond_graph(sites, cell, radius_fn, tolerance, xh_max):
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

    Hydrogen is handled separately and deliberately. The COD-derived radii set
    above has no hydrogen value at all and discards every contact involving
    hydrogen, because radius-based bonding does not work for H: the standard
    epsilon would make H-O bonding MORE generous (1.42 A) than the old 1.2x rule
    (1.16 A), not less. So an H is bonded to exactly ONE heavy atom - its
    nearest within `xh_max` - which matches PLATON PLAT303: "Hydrogen atoms are
    generally connected to only one other atom."

    H...H is never treated as a bond. Dihydrogen complexes exist and are not
    modelled.

    Returns (dists, adj, h_candidates), where h_candidates maps an H index to
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

    for (i, j), (d, _img) in dists.items():
        if is_h[i] and is_h[j]:
            continue
        if is_h[i] or is_h[j]:
            h_idx, x_idx = (i, j) if is_h[i] else (j, i)
            if radii[x_idx] is not None and d <= xh_max:
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

    return (dists, adj, h_candidates)


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
