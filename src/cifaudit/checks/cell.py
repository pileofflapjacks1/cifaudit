"""Unit-cell consistency."""

from typing import List, Tuple

from ..model import ERROR, Finding, INFO, NotTested, WARN

NAME = "cell"


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []
    c = blk.cell

    missing = [n for n in ("a", "b", "c", "alpha", "beta", "gamma")
               if getattr(c, n) is None]
    if missing:
        f.append(Finding(
            code="CELL_INCOMPLETE", severity=ERROR,
            message="Unit-cell parameters missing or unparseable: %s."
                    % ", ".join(missing),
            plat=["PLAT801"],
            based_on=["_cell_length_a", "_cell_length_b", "_cell_length_c",
                      "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"],
            evidence={"missing": missing}))
        nt.append(NotTested(area="cell volume and all geometry",
                            reason="cell is incomplete, so no distance or volume "
                                   "test can be performed",
                            plat=["PLAT150"]))
        return (f, nt)

    for n in ("a", "b", "c"):
        v = getattr(c, n)
        if v <= 0:
            f.append(Finding(code="CELL_AXIS_NONPOSITIVE", severity=ERROR,
                             message="Cell axis %s is %.4f; must be positive." % (n, v),
                             plat=["PLAT801"], based_on=["_cell_length_%s" % n],
                             evidence={"axis": n, "value": v}))
    for n in ("alpha", "beta", "gamma"):
        v = getattr(c, n)
        if not (0.0 < v < 180.0):
            f.append(Finding(code="CELL_ANGLE_OUT_OF_RANGE", severity=ERROR,
                             message="Cell angle %s is %.3f deg; must be in (0, 180)."
                                     % (n, v),
                             plat=["PLAT801"], based_on=["_cell_angle_%s" % n],
                             evidence={"angle": n, "value": v}))

    computed = c.computed_volume()
    if computed is None:
        f.append(Finding(code="CELL_NOT_REALISABLE", severity=ERROR,
                         message="The cell angles do not describe a realisable "
                                 "cell (the metric determinant is not positive).",
                         plat=["PLAT801"],
                         based_on=["_cell_angle_alpha", "_cell_angle_beta",
                                   "_cell_angle_gamma"],
                         evidence={"alpha": c.alpha, "beta": c.beta, "gamma": c.gamma}))
        return (f, nt)

    if c.volume is None:
        f.append(Finding(code="CELL_VOLUME_ABSENT", severity=INFO,
                         message="_cell_volume is not given; computed %.2f A^3."
                                 % computed,
                         plat=["PLAT150"], based_on=["_cell_volume"],
                         evidence={"computed": round(computed, 4)}))
    else:
        rel = abs(c.volume - computed) / computed
        if rel > ctx["volume_rel_tol"]:
            f.append(Finding(
                code="CELL_VOLUME_MISMATCH",
                severity=ERROR if rel > 0.01 else WARN,
                message="Reported _cell_volume %.2f A^3 disagrees with the volume "
                        "computed from the cell parameters, %.2f A^3 (%.2f%%)."
                        % (c.volume, computed, rel * 100.0),
                plat=["PLAT150"],
                based_on=["_cell_volume", "_cell_length_a", "_cell_length_b",
                          "_cell_length_c", "_cell_angle_alpha",
                          "_cell_angle_beta", "_cell_angle_gamma"],
                evidence={"reported": c.volume, "computed": round(computed, 4),
                          "relative_difference": round(rel, 6)}))
        if c.su.get("volume") is None:
            f.append(Finding(code="CELL_VOLUME_SU_ABSENT", severity=INFO,
                             message="_cell_volume has no standard uncertainty.",
                             plat=["PLAT151"], based_on=["_cell_volume"]))

    nt.append(NotTested(
        area="reduced/standard cell setting, cell-vs-diffraction temperature, "
             "axis rounding",
        reason="not implemented",
        plat=["PLAT180"]))
    return (f, nt)
