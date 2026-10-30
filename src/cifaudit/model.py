"""Result types. Deliberately plain: dataclasses + dicts, no pydantic.

Severity is cifaudit's OWN judgement and is not a checkCIF alert level.
checkCIF/PLATON levels are A/B/C/G; the public PLATON test list gives each test a
``Type_1``..``Type_5`` instead, and publishes no mapping between the two. We
therefore report the PLAT code and PLATON's own Type verbatim and assert no level.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

# cifaudit severities. Not checkCIF levels.
ERROR = "error"    # internally inconsistent, or physically impossible as written
WARN = "warn"      # suspicious; a crystallographer should look
INFO = "info"      # worth noticing, commonly benign

SEVERITY_ORDER = {ERROR: 0, WARN: 1, INFO: 2}

OK = "ok"
NOT_ENOUGH_DATA = "not_enough_data"

# Two kinds of finding, which must not be mixed in one count.
#
# STRUCTURAL: something about the structure itself looks wrong.
# REPORTING:  the CIF does not say something, or says it inconsistently. These
#             are about completeness of the deposition, not about the crystal.
#
# Measured over 500 deposited COD structures, reporting findings dominated
# everything: _chemical_formula_moiety was absent in 74% of files and
# _cell_volume had no standard uncertainty in 40%. Counting those alongside NPD
# atoms and impossible contacts made the structural signal unreadable, so they
# are now reported separately. checkCIF makes the same distinction by alert
# level; cifaudit asserts no level, so it needs its own.
STRUCTURAL = "structural"
REPORTING = "reporting"

REPORTING_CODES = frozenset((
    "FORMULA_SUM_ABSENT", "FORMULA_MOIETY_ABSENT",
    "CELL_VOLUME_ABSENT", "CELL_VOLUME_SU_ABSENT",
    "PARTIAL_OCC_NO_DISORDER_GROUP",
    "SPACEGROUP_ABSENT", "SYMOPS_ABSENT", "SPACEGROUP_SETTING_DIFFERS",
    "ADP_ANISO_DECLARED_NO_UIJ", "ADP_UIJ_INCOMPLETE",
    "H_TREATMENT_WITHOUT_H", "OCC_UNPARSEABLE",
))


def kind_of(code):
    return REPORTING if code in REPORTING_CODES else STRUCTURAL


@dataclass
class Finding:
    code: str                      # cifaudit code, e.g. "OCC_GT_ONE"
    severity: str
    message: str
    plat: List[str] = field(default_factory=list)   # e.g. ["PLAT075"]
    based_on: List[str] = field(default_factory=list)  # CIF data names used
    sites: List[str] = field(default_factory=list)     # atom labels involved
    evidence: Dict[str, Any] = field(default_factory=dict)

    @property
    def kind(self) -> str:
        return kind_of(self.code)

    def to_dict(self) -> Dict[str, Any]:
        from .alerts import describe
        return {
            "code": self.code,
            "kind": self.kind,
            "severity": self.severity,
            "message": self.message,
            "checkcif": [describe(p) for p in self.plat],
            "based_on": self.based_on,
            "sites": self.sites,
            "evidence": self.evidence,
        }


@dataclass
class NotTested:
    area: str
    reason: str
    plat: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        from .alerts import describe
        return {
            "area": self.area,
            "reason": self.reason,
            "checkcif_not_covered": [describe(p) for p in self.plat],
        }


@dataclass
class AuditResult:
    status: str
    datablock: Optional[str]
    findings: List[Finding] = field(default_factory=list)
    not_tested: List[NotTested] = field(default_factory=list)
    checked: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    DISCLAIMER = (
        "cifaudit is not checkCIF. It implements a small subset of the "
        "checkCIF/PLATON tests offline and asserts no alert level. Run checkCIF "
        "(https://checkcif.iucr.org/) before publishing. A clean cifaudit report "
        "does not mean the structure is correct."
    )

    def counts(self) -> Dict[str, int]:
        """Severity totals for STRUCTURAL findings only.

        Reporting findings are counted separately in `reporting_counts`, so a
        file whose only issue is a missing optional data name does not read as
        a file with something wrong with its structure.
        """
        out = {ERROR: 0, WARN: 0, INFO: 0}
        for f in self.findings:
            if f.kind == STRUCTURAL:
                out[f.severity] = out.get(f.severity, 0) + 1
        return out

    def reporting_counts(self) -> Dict[str, int]:
        out = {ERROR: 0, WARN: 0, INFO: 0}
        for f in self.findings:
            if f.kind == REPORTING:
                out[f.severity] = out.get(f.severity, 0) + 1
        return out

    def to_dict(self) -> Dict[str, Any]:
        self.findings.sort(key=lambda f: (SEVERITY_ORDER.get(f.severity, 9), f.code))
        return {
            "status": self.status,
            "datablock": self.datablock,
            "counts": self.counts(),
            "reporting_counts": self.reporting_counts(),
            "findings": [f.to_dict() for f in self.findings],
            "checked": sorted(self.checked),
            "not_tested": [n.to_dict() for n in self.not_tested],
            "notes": self.notes,
            "disclaimer": self.DISCLAIMER,
        }
