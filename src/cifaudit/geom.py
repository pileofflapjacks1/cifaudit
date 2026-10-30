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


def bond_graph(sites, cell, bond_factor, radius_fn):
    """Pairwise shortest distances plus the covalent bond graph.

    Returns (dists, adj) where dists maps (i, j) with i < j to (distance, image)
    and adj maps index -> set of bonded indices. Pairs involving an element with
    no tabulated radius are present in dists but never bonded.
    """
    dists = {}
    adj = {i: set() for i in range(len(sites))}
    for i in range(len(sites)):
        ri = radius_fn(sites[i].element)
        for j in range(i + 1, len(sites)):
            rj = radius_fn(sites[j].element)
            d, img = min_distance(cell, sites[i].fract, sites[j].fract)
            dists[(i, j)] = (d, img)
            if ri is not None and rj is not None and d <= (ri + rj) * bond_factor:
                adj[i].add(j)
                adj[j].add(i)
    return (dists, adj)


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
