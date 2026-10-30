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


@dataclass
class Finding:
    code: str                      # cifaudit code, e.g. "OCC_GT_ONE"
    severity: str
    message: str
    plat: List[str] = field(default_factory=list)   # e.g. ["PLAT075"]
    based_on: List[str] = field(default_factory=list)  # CIF data names used
    sites: List[str] = field(default_factory=list)     # atom labels involved
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        from .alerts import describe
        return {
            "code": self.code,
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
        out = {ERROR: 0, WARN: 0, INFO: 0}
        for f in self.findings:
            out[f.severity] = out.get(f.severity, 0) + 1
        return out

    def to_dict(self) -> Dict[str, Any]:
        self.findings.sort(key=lambda f: (SEVERITY_ORDER.get(f.severity, 9), f.code))
        return {
            "status": self.status,
            "datablock": self.datablock,
            "counts": self.counts(),
            "findings": [f.to_dict() for f in self.findings],
            "checked": sorted(self.checked),
            "not_tested": [n.to_dict() for n in self.not_tested],
            "notes": self.notes,
            "disclaimer": self.DISCLAIMER,
        }
