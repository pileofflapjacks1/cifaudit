"""Chemical formula, Z and density against the atom-site list."""

import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from ..expand import ExpandError, expand
from ..model import ERROR, Finding, INFO, NotTested, WARN

NAME = "formula"

_TOKEN = re.compile(r"([A-Z][a-z]?)\s*([0-9]*\.?[0-9]*)")
_MISSING = (None, "", ".", "?")


def parse_formula(text: Optional[str]) -> Dict[str, float]:
    """Parse a _chemical_formula_sum string such as 'C2 H2 O1'."""
    out: Dict[str, float] = defaultdict(float)
    if text in _MISSING:
        return {}
    for sym, num in _TOKEN.findall(str(text)):
        if not sym:
            continue
        out[sym] += float(num) if num not in ("", ".") else 1.0
    return dict(out)


def run(blk, ctx) -> Tuple[List[Finding], List[NotTested]]:
    f: List[Finding] = []
    nt: List[NotTested] = []

    if blk.formula_sum in _MISSING:
        f.append(Finding(code="FORMULA_SUM_ABSENT", severity=WARN,
                         message="_chemical_formula_sum is not given.",
                         plat=["PLAT047"], based_on=["_chemical_formula_sum"]))
    if blk.formula_moiety in _MISSING:
        f.append(Finding(code="FORMULA_MOIETY_ABSENT", severity=INFO,
                         message="_chemical_formula_moiety is not given.",
                         plat=["PLAT048"], based_on=["_chemical_formula_moiety"]))

    reported = parse_formula(blk.formula_sum)
    if not reported:
        nt.append(NotTested(area="formula vs atom sites",
                            reason="no parseable _chemical_formula_sum",
                            plat=["PLAT041"]))
    elif blk.z is None:
        nt.append(NotTested(
            area="formula vs atom sites",
            reason="_cell_formula_units_Z is absent, so the unit-cell content "
                   "cannot be divided into formula units",
            plat=["PLAT041", "PLAT045"]))
    else:
        # Unit-cell contents from symmetry expansion, with ONE exception:
        # dummy atoms. A dummy's coordinates are meaningless, so expanding them
        # gives the wrong multiplicity - COD 1000241 declares multiplicity 32
        # for a hydrogen placed at (-1,-1,-1), which lands on the shared Sb/W
        # site and expands to 16 copies. For those sites the declared
        # multiplicity is used instead.
        #
        # It is NOT used for ordinary sites: `_atom_site_symmetry_multiplicity`
        # is written with two incompatible meanings in the wild - Wyckoff
        # multiplicity in IUCr-style inorganic CIFs, site-symmetry order in
        # SHELX-style ones - and preferring it everywhere made this check fire
        # on 167 of 500 files instead of 94.
        cell_counts = None
        basis = None
        try:
            sites, _cell = expand(blk, drop_dummy=True)
        except ExpandError as exc:
            nt.append(NotTested(area="formula vs atom sites",
                                reason=str(exc), plat=["PLAT041"]))
            sites = None
        if sites is not None:
            cell_counts = defaultdict(float)
            for s in sites:
                if s.element:
                    cell_counts[s.element] += (s.occ if s.occ is not None else 1.0)
            basis = "symmetry expansion of the atom-site list"
            dummies = [s for s in blk.sites if s.is_dummy and s.element]
            if dummies:
                for s in dummies:
                    if s.multiplicity:
                        occ = s.occupancy if s.occupancy is not None else 1.0
                        cell_counts[s.element] += s.multiplicity * occ
                basis += (", plus declared multiplicity x occupancy for %d dummy "
                          "site(s) whose coordinates are undetermined"
                          % len(dummies))

        if cell_counts is not None:
            if blk.z and blk.z > 0:
                per_unit = {k: v / blk.z for k, v in cell_counts.items()}
                diffs = []
                for el in sorted(set(list(per_unit) + list(reported))):
                    got = per_unit.get(el, 0.0)
                    want = reported.get(el, 0.0)
                    if abs(got - want) > ctx["formula_tol"]:
                        diffs.append({"element": el,
                                      "from_sites": round(got, 3),
                                      "from_formula": round(want, 3)})
                if diffs:
                    f.append(Finding(
                        code="FORMULA_SITES_MISMATCH", severity=WARN,
                        message="The atom-site list, divided by Z=%g, does not "
                                "match _chemical_formula_sum for: %s. (Counts "
                                "from the %s.)"
                                % (blk.z, ", ".join(
                                    "%s %.2f vs %.2f" % (d["element"],
                                                         d["from_sites"],
                                                         d["from_formula"])
                                    for d in diffs[:8]), basis),
                        plat=["PLAT041", "PLAT043"],
                        based_on=["_chemical_formula_sum", "_cell_formula_units_Z",
                                  "_atom_site_occupancy", "_atom_site_type_symbol",
                                  "_atom_site_symmetry_multiplicity"],
                        evidence={"z": blk.z, "differences": diffs[:12],
                                  "basis": basis,
                                  "caveat": "rounding of occupancies in the CIF "
                                            "can produce small non-integral "
                                            "differences"}))
                non_integral = [
                    {"element": el, "count": round(v, 3)}
                    for el, v in cell_counts.items()
                    if abs(v - round(v)) > ctx["formula_tol"]]
                if non_integral:
                    f.append(Finding(
                        code="NON_INTEGRAL_CELL_CONTENT", severity=INFO,
                        message="Unit-cell contents are non-integral for: %s. "
                                "Common with disorder or SQUEEZE'd solvent."
                                % ", ".join("%s %.2f" % (d["element"], d["count"])
                                            for d in non_integral[:8]),
                        plat=["PLAT077"],
                        based_on=["_atom_site_occupancy"],
                        evidence={"counts": non_integral[:12]}))

    if blk.density_reported is not None and blk.density_reported < 1.0:
        f.append(Finding(
            code="DENSITY_BELOW_ONE", severity=WARN,
            message="Reported _exptl_crystal_density_diffrn is %.4f g/cm3. Values "
                    "below 1.0 usually mean an incomplete model or wrong symmetry."
                    % blk.density_reported,
            plat=["PLAT049"], based_on=["_exptl_crystal_density_diffrn"],
            evidence={"density": blk.density_reported}))

    nt.append(NotTested(
        area="moiety-formula string comparison, molecular weight, F(000), "
             "reported-vs-calculated density and Z",
        reason="not implemented; needs atomic weights and scattering tables",
        plat=["PLAT042", "PLAT043", "PLAT044", "PLAT045", "PLAT046"]))
    return (f, nt)
