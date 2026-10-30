"""Atomic displacement parameters.

NPD criterion. The anisotropic Debye-Waller factor is
    T(h) = exp(-2*pi^2 * sum_ij U_ij h_i h_j a*_i a*_j),
so the quadratic form acts on reciprocal indices scaled by the reciprocal axis
lengths. That scaling is a congruence by a positive diagonal matrix, which
preserves positive-definiteness. Therefore testing the eigenvalues of the U_ij
matrix exactly as given in the CIF is the correct NPD test for ANY cell, with no
transformation to Cartesian axes required.
"""

import math
from collections import defaultdict
from typing import List, Tuple

from ..model import ERROR, Finding, INFO, NotTested, WARN

NAME = "adp"

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None

_UIJ = ("U11", "U22", "U33", "U12", "U13", "U23")


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    if not blk.sites:
        nt.append(NotTested(area="ADPs", reason="no atom sites in the CIF",
                            plat=["PLAT211", "PLAT220"]))
        return (f, nt)

    aniso = [s for s in blk.sites if s.uij]
    declared_aniso = [s for s in blk.sites
                      if (s.adp_type or "").lower() in ("uani", "bani")]

    for s in declared_aniso:
        if not s.uij:
            f.append(Finding(
                code="ADP_ANISO_DECLARED_NO_UIJ", severity=ERROR,
                message="Site %s is declared Uani but has no entry in the "
                        "_atom_site_aniso_* loop." % s.label,
                plat=["PLAT217", "PLAT806"],
                based_on=["_atom_site_adp_type", "_atom_site_aniso_label"],
                sites=[s.label]))

    for s in aniso:
        have = [k for k in _UIJ if k in s.uij]
        if len(have) < 6:
            f.append(Finding(
                code="ADP_UIJ_INCOMPLETE", severity=ERROR,
                message="Site %s has only %d of 6 U(i,j) components (%s)."
                        % (s.label, len(have), ", ".join(have)),
                plat=["PLAT217", "PLAT806"],
                based_on=["_atom_site_aniso_U_11", "_atom_site_aniso_U_22",
                          "_atom_site_aniso_U_33", "_atom_site_aniso_U_12",
                          "_atom_site_aniso_U_13", "_atom_site_aniso_U_23"],
                sites=[s.label], evidence={"present": have}))
            continue
        if gemmi is None:
            continue
        m = gemmi.SMat33d(s.uij["U11"], s.uij["U22"], s.uij["U33"],
                          s.uij["U12"], s.uij["U13"], s.uij["U23"])
        eig = sorted(m.calculate_eigenvalues())
        if eig[0] <= 0.0:
            f.append(Finding(
                code="ADP_NPD", severity=ERROR,
                message="Site %s is non-positive-definite (NPD): smallest "
                        "eigenvalue of U(i,j) is %.5f." % (s.label, eig[0]),
                plat=["PLAT211", "PLAT212"],
                based_on=["_atom_site_aniso_U_11", "_atom_site_aniso_U_22",
                          "_atom_site_aniso_U_33", "_atom_site_aniso_U_12",
                          "_atom_site_aniso_U_13", "_atom_site_aniso_U_23"],
                sites=[s.label],
                evidence={"eigenvalues": [round(e, 6) for e in eig]}))
        else:
            # PLAT213: "ordered: U1 < U2 < U3. The value of SQRT(U3/U1) main
            # axis ADP ratio (Angstrom Units) is tested". U is a mean-square
            # displacement, so the square root is the ratio of axis lengths.
            axis_ratio = math.sqrt(eig[2] / eig[0])
            if axis_ratio > ctx["adp_axis_ratio_max"]:
                # PLAT213: "Oblate criterium: U3 - U2 < U2 - U1. Prolate
                # otherwise."
                shape = ("oblate" if (eig[2] - eig[1]) < (eig[1] - eig[0])
                         else "prolate")
                f.append(Finding(
                    code="ADP_EXTREME_ANISOTROPY", severity=WARN,
                    message="Site %s has a principal-axis ratio "
                            "sqrt(U3/U1) = %.2f (%s ellipsoid). Above %.2f, "
                            "which is the 99th percentile over 12,899 non-H "
                            "anisotropic atoms in 500 deposited COD structures. "
                            "Large values may indicate unresolved disorder."
                            % (s.label, axis_ratio, shape,
                               ctx["adp_axis_ratio_max"]),
                    plat=["PLAT213", "PLAT215", "PLAT250"],
                    based_on=["_atom_site_aniso_U_11"], sites=[s.label],
                    evidence={"eigenvalues": [round(e, 6) for e in eig],
                              "axis_ratio_sqrt_U3_U1": round(axis_ratio, 3),
                              "raw_eigenvalue_ratio": round(eig[2] / eig[0], 3),
                              "shape": shape}))

    # Ueq spread, PER ELEMENT TYPE. PLAT220's purpose is detecting misassigned
    # atomic species, not overall spread: "Too high or too low Ueq's may be an
    # indication for incorrectly identified atomic species (i.e. O versus N)".
    # Comparing max/min across all non-H atoms together asked a different
    # question and fired on 7% of real files.
    by_element = defaultdict(list)
    for s in blk.sites:
        if s.is_hydrogen or s.uiso is None or s.uiso <= 0:
            continue
        if s.element:
            by_element[s.element].append((s.label, s.uiso))
    tested_any = False
    for element, vals in sorted(by_element.items()):
        if len(vals) < ctx["ueq_min_atoms"]:
            continue
        tested_any = True
        lo = min(vals, key=lambda t: t[1])
        hi = max(vals, key=lambda t: t[1])
        if hi[1] / lo[1] > ctx["ueq_ratio_max"]:
            f.append(Finding(
                code="UEQ_RANGE_WIDE", severity=WARN,
                message="Ueq across the %d %s atoms spans a factor of %.1f "
                        "(%s=%.4f to %s=%.4f). A wide range within one element "
                        "can mean a site has been assigned the wrong species."
                        % (len(vals), element, hi[1] / lo[1],
                           lo[0], lo[1], hi[0], hi[1]),
                plat=["PLAT220", "PLAT221"],
                based_on=["_atom_site_U_iso_or_equiv",
                          "_atom_site_type_symbol"],
                sites=[lo[0], hi[0]],
                evidence={"element": element, "n": len(vals),
                          "min": lo[1], "max": hi[1],
                          "ratio": round(hi[1] / lo[1], 3)}))
    if not tested_any:
        nt.append(NotTested(
            area="Ueq(max)/Ueq(min) range per element type",
            reason="no element has at least %d non-H atoms with a usable "
                   "_atom_site_U_iso_or_equiv" % ctx["ueq_min_atoms"],
            plat=["PLAT220", "PLAT221"]))

    for s in blk.sites:
        if s.uiso is not None and s.uiso < 0:
            f.append(Finding(
                code="UISO_NEGATIVE", severity=ERROR,
                message="Site %s has negative U_iso_or_equiv %.5f."
                        % (s.label, s.uiso),
                plat=["PLAT211"], based_on=["_atom_site_U_iso_or_equiv"],
                sites=[s.label], evidence={"u_iso": s.uiso}))

    nt.append(NotTested(
        area="Hirshfeld rigid-bond test; ADP restraint reporting; negative "
             "density at H positions",
        reason="requires refinement s.u. propagation and/or the structure-factor "
               "file, which cifaudit does not read",
        plat=["PLAT230", "PLAT003", "PLAT977"]))
    return (f, nt)
