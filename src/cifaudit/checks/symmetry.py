"""Space-group symbol, symmetry operators, and higher symmetry missed."""

from typing import List, Tuple

from ..model import ERROR, Finding, INFO, NotTested, WARN

NAME = "symmetry"

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None

try:
    import spglib
    HAVE_SPGLIB = True
except ImportError:  # pragma: no cover
    spglib = None
    HAVE_SPGLIB = False

_MISSING = (None, "", ".", "?")


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    if blk.hm_symbol in _MISSING:
        f.append(Finding(
            code="SPACEGROUP_ABSENT", severity=ERROR,
            message="No space-group Hermann-Mauguin symbol is given.",
            plat=["PLAT122"],
            based_on=["_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt"]))
    if not blk.symops:
        f.append(Finding(
            code="SYMOPS_ABSENT", severity=ERROR,
            message="No symmetry operators are listed.",
            plat=["PLAT124"],
            based_on=["_symmetry_equiv_pos_as_xyz",
                      "_space_group_symop_operation_xyz"]))

    sg = None
    if blk.hm_symbol not in _MISSING and gemmi is not None:
        try:
            sg = gemmi.SpaceGroup(blk.hm_symbol)
        except Exception:
            sg = None
        if sg is None:
            f.append(Finding(
                code="SPACEGROUP_UNINTERPRETABLE", severity=ERROR,
                message="Space-group symbol %r could not be interpreted."
                        % blk.hm_symbol,
                plat=["PLAT121", "PLAT123"],
                based_on=["_symmetry_space_group_name_H-M"],
                evidence={"symbol": blk.hm_symbol}))

    if sg is not None and blk.symops:
        declared = len(blk.symops)
        expected = len(sg.operations())
        if declared != expected:
            f.append(Finding(
                code="SYMOP_COUNT_MISMATCH", severity=ERROR,
                message="Space group %s implies %d symmetry operators but %d "
                        "are listed." % (blk.hm_symbol, expected, declared),
                plat=["PLAT120"],
                based_on=["_symmetry_space_group_name_H-M",
                          "_symmetry_equiv_pos_as_xyz"],
                evidence={"declared": declared, "expected": expected,
                          "space_group": blk.hm_symbol}))

    if sg is not None and blk.hall_symbol not in _MISSING:
        try:
            if gemmi.SpaceGroup(blk.hall_symbol) is None:
                pass
        except Exception:
            pass

    # ADDSYM-equivalent: does a higher-symmetry group fit the coordinates?
    if not HAVE_SPGLIB:
        nt.append(NotTested(
            area="missed higher symmetry (ADDSYM-equivalent)",
            reason="spglib is not installed",
            plat=["PLAT104", "PLAT110", "PLAT111", "PLAT112", "PLAT113"]))
    elif not blk.cell.complete() or not blk.sites:
        nt.append(NotTested(
            area="missed higher symmetry (ADDSYM-equivalent)",
            reason="needs a complete cell and atom coordinates",
            plat=["PLAT110", "PLAT111", "PLAT112"]))
    else:
        usable = [s for s in blk.sites
                  if all(v is not None for v in s.fract) and s.element]
        if len(usable) < 1:
            nt.append(NotTested(
                area="missed higher symmetry (ADDSYM-equivalent)",
                reason="no atom sites with complete fractional coordinates",
                plat=["PLAT110", "PLAT111", "PLAT112"]))
        else:
            c = blk.cell
            lattice = _lattice(c)
            positions = [list(s.fract) for s in usable]
            species = [abs(hash(s.element)) % 1000 for s in usable]
            try:
                ds = spglib.get_symmetry_dataset(
                    (lattice, positions, species), symprec=ctx["symprec"])
            except Exception as exc:
                ds = None
                nt.append(NotTested(
                    area="missed higher symmetry (ADDSYM-equivalent)",
                    reason="spglib failed: %s" % str(exc)[:120],
                    plat=["PLAT110", "PLAT111", "PLAT112"]))
            if ds:
                derived_no = _ds(ds, "number")
                derived_sym = _ds(ds, "international")
                if sg is not None and derived_no and derived_no != sg.number:
                    sev = WARN if derived_no > sg.number else INFO
                    f.append(Finding(
                        code="SYMMETRY_DERIVED_DIFFERS", severity=sev,
                        message="The coordinates as deposited are consistent with "
                                "space group %s (No. %d) at symprec=%g, but the CIF "
                                "declares %s (No. %d). %s"
                                % (derived_sym, derived_no, ctx["symprec"],
                                   blk.hm_symbol, sg.number,
                                   "This is the direction that suggests missed "
                                   "symmetry." if derived_no > sg.number else
                                   "The declared group is the higher one."),
                        plat=["PLAT104", "PLAT110", "PLAT111", "PLAT112", "PLAT113"],
                        based_on=["_symmetry_space_group_name_H-M",
                                  "_atom_site_fract_x", "_atom_site_fract_y",
                                  "_atom_site_fract_z"],
                        evidence={"declared": blk.hm_symbol,
                                  "declared_number": sg.number,
                                  "derived": derived_sym,
                                  "derived_number": derived_no,
                                  "symprec": ctx["symprec"],
                                  "caveat": "spglib was given the asymmetric unit "
                                            "only; disorder and partial occupancy "
                                            "are ignored in this comparison"}))
                nt.append(NotTested(
                    area="non-standard space-group setting",
                    reason="not implemented", plat=["PLAT128"]))
    return (f, nt)


def _ds(ds, key):
    """spglib returns a dict in <=2.4 and a dataclass in newer versions."""
    if ds is None:
        return None
    if isinstance(ds, dict):
        return ds.get(key)
    return getattr(ds, key, None)


def _lattice(c):
    import math
    al, be, ga = (math.radians(c.alpha), math.radians(c.beta),
                  math.radians(c.gamma))
    cx = c.c * math.cos(be)
    cy = c.c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz2 = c.c * c.c - cx * cx - cy * cy
    cz = math.sqrt(cz2) if cz2 > 0 else 0.0
    return [[c.a, 0.0, 0.0],
            [c.b * math.cos(ga), c.b * math.sin(ga), 0.0],
            [cx, cy, cz]]
