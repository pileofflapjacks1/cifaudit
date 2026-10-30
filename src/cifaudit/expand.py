"""Symmetry expansion to the full unit cell, via gemmi's small-structure path."""

from dataclasses import dataclass
from typing import List, Optional, Tuple

from .load import element_of, sanitise

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None

MAX_EXPANDED = 4000


@dataclass
class Expanded:
    label: str
    element: Optional[str]
    fract: Tuple[float, float, float]
    occ: float
    disorder_group: Optional[str] = None
    disorder_assembly: Optional[str] = None


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
    try:
        sites = st.get_all_unit_cell_sites()
    except Exception as exc:
        raise ExpandError("symmetry expansion failed: %s" % sanitise(str(exc), 160))
    if len(sites) > MAX_EXPANDED:
        raise ExpandError("unit cell expands to %d sites, above the %d cap"
                          % (len(sites), MAX_EXPANDED))
    # Atoms on special positions can come back as coincident copies, which then
    # look like two atoms 0.000 A apart. Deduplicate on label plus rounded
    # fractional coordinates, wrapped into [0, 1).
    # Per-label disorder labels and dummy flags come from the CIF loop, since
    # gemmi's small-structure sites do not carry the assembly or the calc flag.
    meta = {}
    for s0 in blk.sites:
        meta[s0.label] = (s0.disorder_group, s0.disorder_assembly, s0.is_dummy)

    out = []
    seen = set()
    dropped_dummy = 0
    for s in sites:
        lab = sanitise(s.label) or "?"
        dgroup, dasm, dummy = meta.get(lab, (None, None, False))
        if dummy and drop_dummy:
            dropped_dummy += 1
            continue
        fx, fy, fz = s.fract.x, s.fract.y, s.fract.z
        key = (lab,
               round(fx % 1.0, 5), round(fy % 1.0, 5), round(fz % 1.0, 5))
        if key in seen:
            continue
        seen.add(key)
        out.append(Expanded(
            label=lab,
            element=element_of(s.type_symbol, s.label),
            fract=(fx, fy, fz),
            occ=float(s.occ) if s.occ is not None else 1.0,
            disorder_group=dgroup, disorder_assembly=dasm))
    return (out, st.cell)
