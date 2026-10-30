"""Atom-site labels and scattering types."""

from collections import Counter
from typing import List, Tuple

from ..elements import COVALENT
from ..load import element_of
from ..model import ERROR, Finding, NotTested, WARN

NAME = "labels"


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []
    if not blk.sites:
        nt.append(NotTested(area="atom labels", reason="no atom sites in the CIF",
                            plat=["PLAT070"]))
        return (f, nt)

    counts = Counter(s.label for s in blk.sites)
    dupes = sorted(l for l, n in counts.items() if n > 1 and l)
    if dupes:
        f.append(Finding(
            code="DUPLICATE_LABELS", severity=ERROR,
            message="Duplicate atom-site labels: %s. Any per-label association "
                    "(anisotropic U, bonds, H-bonds) is ambiguous."
                    % ", ".join(dupes[:10]),
            plat=["PLAT070"], based_on=["_atom_site_label"], sites=dupes[:10],
            evidence={"duplicates": {l: counts[l] for l in dupes[:10]}}))

    empty = [i for i, s in enumerate(blk.sites) if not s.label]
    if empty:
        f.append(Finding(code="EMPTY_LABEL", severity=ERROR,
                         message="%d atom site(s) have an empty label." % len(empty),
                         plat=["PLAT070"], based_on=["_atom_site_label"],
                         evidence={"row_indices": empty[:10]}))

    if "_atom_site_type_symbol" in blk.present:
        for s in blk.sites:
            if not s.type_symbol:
                continue
            from_type = element_of(s.type_symbol, None)
            # A label is consistent with the scattering type when it STARTS WITH
            # that element symbol. Deriving an element from the label greedily is
            # wrong: "HO1" is a hydroxyl H named for O1, not holmium, and that
            # misread produced a false positive on COD 2013000.
            lab = (s.label or "").upper()
            consistent = bool(from_type) and lab.startswith(from_type.upper())
            from_label = element_of(None, s.label)
            if from_type and not consistent:
                f.append(Finding(
                    code="SCATTERING_TYPE_LABEL_MISMATCH", severity=WARN,
                    message="Site %s has _atom_site_type_symbol %r, which does not "
                            "agree with the element implied by the label (%s)."
                            % (s.label, s.type_symbol, from_label or "?"),
                    plat=["PLAT017"],
                    based_on=["_atom_site_label", "_atom_site_type_symbol"],
                    sites=[s.label],
                    evidence={"from_type_symbol": from_type,
                              "from_label": from_label}))
    else:
        nt.append(NotTested(area="scattering type consistency",
                            reason="_atom_site_type_symbol is absent",
                            plat=["PLAT017"]))

    unknown = sorted({s.element for s in blk.sites
                      if s.element and s.element not in COVALENT})
    if unknown:
        nt.append(NotTested(
            area="any distance-based check involving %s" % ", ".join(unknown),
            reason="no covalent radius is tabulated for these elements, so "
                   "bonded/non-bonded cannot be decided for them",
            plat=["PLAT310", "PLAT432"]))
    return (f, nt)
