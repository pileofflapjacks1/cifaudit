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


class ExpandError(Exception):
    pass


def expand(blk) -> Tuple[List[Expanded], Optional[object]]:
    """Return (all unit-cell sites, gemmi UnitCell).

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
    out = []
    seen = set()
    for s in sites:
        fx, fy, fz = s.fract.x, s.fract.y, s.fract.z
        key = (sanitise(s.label) or "?",
               round(fx % 1.0, 5), round(fy % 1.0, 5), round(fz % 1.0, 5))
        if key in seen:
            continue
        seen.add(key)
        out.append(Expanded(
            label=sanitise(s.label) or "?",
            element=element_of(s.type_symbol, s.label),
            fract=(fx, fy, fz),
            occ=float(s.occ) if s.occ is not None else 1.0))
    return (out, st.cell)
