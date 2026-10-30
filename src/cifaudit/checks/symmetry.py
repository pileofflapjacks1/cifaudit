"""Space-group symbol, symmetry operators, and higher symmetry missed."""

from typing import List, Tuple

from ..expand import ExpandError, expand
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
    #
    # spglib DERIVES symmetry from the atoms it is given, so it must be given the
    # whole unit cell. Passing the asymmetric unit makes it derive P1 every time:
    # measured on 16 deposited COD structures, that produced a finding on 16/16,
    # i.e. a 100% false-positive rate. Expanding first recovers the declared
    # space-group number on all of them.
    if not HAVE_SPGLIB:
        nt.append(NotTested(
            area="missed higher symmetry (ADDSYM-equivalent)",
            reason="spglib is not installed",
            plat=["PLAT104", "PLAT110", "PLAT111", "PLAT112", "PLAT113"]))
    elif not blk.cell.complete():
        nt.append(NotTested(
            area="missed higher symmetry (ADDSYM-equivalent)",
            reason="needs a complete cell",
            plat=["PLAT110", "PLAT111", "PLAT112"]))
    else:
        try:
            sites, _cell = expand(blk)
        except ExpandError as exc:
            sites = None
            nt.append(NotTested(
                area="missed higher symmetry (ADDSYM-equivalent)",
                reason="the unit cell could not be expanded, and the asymmetric "
                       "unit alone always derives P1: %s" % exc,
                plat=["PLAT110", "PLAT111", "PLAT112"]))
        if sites is not None:
            usable = [s for s in sites
                      if s.element and all(v is not None for v in s.fract)]
            if not usable:
                nt.append(NotTested(
                    area="missed higher symmetry (ADDSYM-equivalent)",
                    reason="no expanded sites with complete coordinates",
                    plat=["PLAT110", "PLAT111", "PLAT112"]))
            else:
                # Deterministic species numbering. str.hash is salted per process,
                # so hashing element symbols would make spglib's answer vary
                # between runs.
                order = sorted({s.element for s in usable})
                species = [order.index(s.element) + 1 for s in usable]
                positions = [list(s.fract) for s in usable]
                lattice = _lattice(blk.cell)
                try:
                    ds = spglib.get_symmetry_dataset(
                        (lattice, positions, species), symprec=ctx["symprec"])
                except Exception as exc:
                    ds = None
                    nt.append(NotTested(
                        area="missed higher symmetry (ADDSYM-equivalent)",
                        reason="spglib failed: %s" % str(exc)[:120],
                        plat=["PLAT110", "PLAT111", "PLAT112"]))
                derived_no = _ds(ds, "number")
                derived_sym = _ds(ds, "international")
                if ds is not None and derived_no and sg is not None:
                    if derived_no == sg.number:
                        # Same group. A different symbol is a setting difference.
                        if derived_sym and _norm(derived_sym) != _norm(blk.hm_symbol):
                            f.append(Finding(
                                code="SPACEGROUP_SETTING_DIFFERS", severity=INFO,
                                message="The declared and derived space groups are "
                                        "the same group (No. %d), but the symbols "
                                        "differ: CIF says %r, spglib's standard "
                                        "setting is %r. A setting difference, not "
                                        "a symmetry error."
                                        % (sg.number, blk.hm_symbol, derived_sym),
                                plat=["PLAT128"],
                                based_on=["_symmetry_space_group_name_H-M"],
                                evidence={"declared": blk.hm_symbol,
                                          "derived": derived_sym,
                                          "number": sg.number}))
                    elif derived_no > sg.number:
                        f.append(Finding(
                            code="SYMMETRY_POSSIBLY_MISSED", severity=WARN,
                            message="The expanded unit-cell contents are consistent "
                                    "with space group %s (No. %d) at symprec=%g, "
                                    "higher than the declared %s (No. %d). Run "
                                    "PLATON ADDSYM to confirm."
                                    % (derived_sym, derived_no, ctx["symprec"],
                                       blk.hm_symbol, sg.number),
                            plat=["PLAT104", "PLAT110", "PLAT111", "PLAT112",
                                  "PLAT113"],
                            based_on=["_symmetry_space_group_name_H-M",
                                      "_atom_site_fract_x", "_atom_site_fract_y",
                                      "_atom_site_fract_z"],
                            evidence={"declared": blk.hm_symbol,
                                      "declared_number": sg.number,
                                      "derived": derived_sym,
                                      "derived_number": derived_no,
                                      "symprec": ctx["symprec"],
                                      "caveat": "disorder and partial occupancy "
                                                "are ignored in this comparison"}))
                    else:
                        f.append(Finding(
                            code="SYMMETRY_DERIVED_LOWER", severity=INFO,
                            message="The expanded contents derive %s (No. %d) at "
                                    "symprec=%g, LOWER than the declared %s "
                                    "(No. %d). Usually disorder, partial occupancy "
                                    "or a loose symprec rather than an error."
                                    % (derived_sym, derived_no, ctx["symprec"],
                                       blk.hm_symbol, sg.number),
                            plat=["PLAT112"],
                            based_on=["_symmetry_space_group_name_H-M"],
                            evidence={"declared_number": sg.number,
                                      "derived_number": derived_no,
                                      "symprec": ctx["symprec"]}))
                nt.append(NotTested(area="non-standard space-group setting beyond "
                                         "the symbol comparison above",
                                    reason="not implemented", plat=["PLAT128"]))
    return (f, nt)


def _norm(sym):
    """Compare H-M symbols ignoring spaces, case and explicit '1' axes."""
    if not sym:
        return ""
    s = str(sym).replace(" ", "").replace("_", "").upper()
    return s.replace("1", "") if len(s.replace("1", "")) > 1 else s


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
