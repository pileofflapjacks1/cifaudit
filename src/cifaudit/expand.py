"""Symmetry expansion to the full unit cell, via gemmi's small-structure path."""

from dataclasses import dataclass
from typing import List, Optional, Tuple

from .load import element_of, sanitise

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None

MAX_EXPANDED = 4000
# Two symmetry copies of one atom closer than this are ONE site (an atom on a
# special position, with coordinates rounded in the CIF). Rounding to 4-5
# decimals moves a copy by well under 0.01 A; 0.05 A leaves a margin. Raising it
# to 0.15 A for all atoms re-merged genuine copies 0.05-0.15 A apart in COD
# 9005049 and 7203095 and broke their formula counts again.
COINCIDENT_COPY_MAX = 0.05
# A SHELX negative-PART atom refined a few hundredths of an Angstrom off the
# symmetry element it is disordered across is still one site: COD 7013013 has
# PART -1 atoms 0.039 and 0.068 A from their mirror images, and its formula
# matches only when all of them are merged. A PART -33 water in COD 7203096
# has copies 0.15 A apart that ARE distinct (its formula needs both), so the
# limit sits between the two, at the 0.10 A used for coincident sites
# elsewhere.
COINCIDENT_COPY_MAX_NEGATIVE_PART = 0.10


@dataclass
class Expanded:
    label: str
    element: Optional[str]
    fract: Tuple[float, float, float]
    occ: float
    disorder_group: Optional[str] = None
    disorder_assembly: Optional[str] = None
    # Indices into the symmetry-operator list of every operator that maps the
    # asymmetric-unit atom onto this site: one for a general position, several
    # for a special position. Two sites with no operator in common belong to
    # different symmetry copies of the asymmetric unit.
    ops: Tuple[int, ...] = ()


class ExpandError(Exception):
    pass


def expand(blk, drop_dummy=True) -> Tuple[List[Expanded], Optional[object]]:
    """Return (all unit-cell sites, gemmi UnitCell).

    `drop_dummy` removes sites flagged `_atom_site_calc_flag` = `dum`. A dummy
    atom's POSITION is undetermined, so it must not take part in any distance or
    connectivity test - but it IS part of the composition, so the formula check
    passes drop_dummy=False. COD 1000241 declares 'H3 O7 Sb W' with Z=8 and
    carries its hydrogen as a dummy at (-1,-1,-1) with occupancy 0.75 on a
    32-fold site: 32 x 0.75 / 8 = 3.0 H, exactly the declared formula.

    Raises ExpandError when the CIF cannot be expanded, so callers abstain rather
    than silently auditing only the asymmetric unit.
    """
    if gemmi is None:
        raise ExpandError("gemmi is not installed")
    if blk._raw is None:
        raise ExpandError("no raw datablock retained")
    try:
        st = gemmi.make_small_structure_from_block(blk._raw)
    except Exception as exc:
        raise ExpandError("gemmi could not build a structure: %s"
                          % sanitise(str(exc), 160))
    if st.spacegroup is None and not blk.symops:
        raise ExpandError("no usable symmetry; only the asymmetric unit is known")
    ops = _operators(st, blk)
    if not ops:
        raise ExpandError("no usable symmetry operators")
    # Per-label disorder labels and dummy flags come from the CIF loop, since
    # gemmi's small-structure sites do not carry the assembly or the calc flag.
    meta = {}
    for s0 in blk.sites:
        meta[s0.label] = (s0.disorder_group, s0.disorder_assembly, s0.is_dummy)

    # The expansion is done here rather than with gemmi's
    # get_all_unit_cell_sites(), which merges any copy within a fixed tolerance
    # of about 0.4 A of another copy of the same atom. That is right for atoms
    # ON a special position but wrong for atoms disordered ACROSS one: in COD
    # 7041079 three half-occupied H (SHELX PART -1) lie 0.29-0.38 A from their
    # own mirror images, gemmi kept one copy of each, and the H count came out
    # 25 per formula instead of the declared 28. Only copies that coincide
    # within COINCIDENT_COPY_MAX are merged here (wider for negative PART).
    #
    # Every coordinate is wrapped into [0, 1). The distance search only tries
    # lattice shifts of -1, 0 and +1, so unwrapped copies more than ~1.5 cells
    # apart lost their true nearest image: in the precision study a Br 2.42 A
    # from Zn was reported as isolated.
    out = []
    for s in st.sites:
        lab = sanitise(s.label) or "?"
        dgroup, dasm, dummy = meta.get(lab, (None, None, False))
        if dummy and drop_dummy:
            continue
        element = element_of(s.type_symbol, s.label)
        occ = float(s.occ) if s.occ is not None else 1.0
        if occ == 0.0:
            # Zero occupancy contributes nothing to composition and is not an
            # atom for geometry. COD 7035035 carries a Cp-ring centroid CT1,
            # typed H with occupancy 0 and no 'dum' flag, which was read as a
            # hydrogen bonded to five ring carbons. OCC_ZERO still reports the
            # site from the atom list.
            continue
        merge_max = (COINCIDENT_COPY_MAX_NEGATIVE_PART
                     if (dgroup or "").strip().startswith("-")
                     else COINCIDENT_COPY_MAX)
        copies = []      # [fract, [op indices]]
        xyz = [s.fract.x, s.fract.y, s.fract.z]
        if any(v != v for v in xyz):
            continue        # no position (NaN): nothing to expand
        for k, op in enumerate(ops):
            f = [v % 1.0 for v in op.apply_to_xyz(xyz)]
            f = [0.0 if v >= 1.0 - 1e-12 else v for v in f]
            for c in copies:
                if _sep(st.cell, c[0], f) <= merge_max:
                    c[1].append(k)
                    break
            else:
                copies.append([f, [k]])
        for f, ks in copies:
            out.append(Expanded(
                label=lab, element=element, fract=(f[0], f[1], f[2]), occ=occ,
                disorder_group=dgroup, disorder_assembly=dasm, ops=tuple(ks)))
        if len(out) > MAX_EXPANDED:
            raise ExpandError("unit cell expands to more than %d sites"
                              % MAX_EXPANDED)
    return (out, st.cell)


def _operators(st, blk):
    """All symmetry operators, centring included, as gemmi Ops.

    The CIF's own operator list is used when it parses completely, since it is
    what the authors refined with; otherwise the space group's operators.
    """
    raw = list(getattr(st, "symops", None) or [])
    if raw:
        try:
            return [gemmi.parse_triplet(t) for t in raw]
        except Exception:
            pass
    if st.spacegroup is not None:
        return list(st.spacegroup.operations())
    return []


def _sep(cell, f1, f2):
    """Distance between two nearby fractional points, nearest lattice image."""
    d = [a - b for a, b in zip(f1, f2)]
    d = [v - round(v) for v in d]
    p = cell.orthogonalize(gemmi.Fractional(d[0], d[1], d[2]))
    return (p.x * p.x + p.y * p.y + p.z * p.z) ** 0.5
