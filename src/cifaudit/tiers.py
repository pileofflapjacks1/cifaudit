"""How often each structural check has been right, and the tier that follows.

Severity says how bad a finding would be if it were true. The tier says how
often a finding from that check HAS been true when sampled. They are separate:
ATOMS_TOO_CLOSE is an error when real, but in four precision rounds it was real
2 times in 12.

Counts are pooled over the four pre-registered precision rounds
(data/processed/precision_*results.md): 237 findings drawn from 2,800 COD
structures never used for tuning, each judged against the full CIF by an
independent model adjudicator. That is model adjudication, not expert review.
The rounds measured four successive code versions, so these are the checks'
behaviour across versions, not one release. tests/test_tiers.py recomputes every
number from the published sample and verdict files.

Tiers:

* RELIABLE   - judged at least MIN_JUDGED times, real in at least half.
* PROMPT     - judged at least MIN_JUDGED times, real in fewer than half. A
               prompt to look, not a finding to act on.
* UNMEASURED - judged fewer than MIN_JUDGED times, never drawn, or changed in
               definition since it was measured. No claim either way.
"""

import math
from typing import Any, Dict, Optional

RELIABLE = "reliable"
PROMPT = "prompt"
UNMEASURED = "unmeasured"

MIN_JUDGED = 8
RELIABLE_MIN = 0.5

#: code -> (real, judged), unsure excluded, pooled over rounds 1-4.
PRECISION = {
    "ADP_EXTREME_ANISOTROPY": (14, 16),
    "UEQ_RANGE_WIDE": (13, 15),
    "NO_HYDROGENS_WITH_CARBON": (6, 8),
    "FORMULA_SITES_MISMATCH": (24, 40),
    "SHORT_HH_CONTACT": (5, 16),
    "SHORT_XH": (3, 15),
    "H_MULTIPLE_CONNECTIONS": (2, 14),
    "ATOMS_TOO_CLOSE": (2, 12),
    "DISORDER_SUM_NOT_ONE": (0, 12),
    "SHORT_NONBONDED_CONTACT": (2, 11),
    "ISOLATED_ATOM": (4, 40),
    "DENSITY_BELOW_ONE": (2, 8),
    "NON_INTEGRAL_CELL_CONTENT": (0, 10),
    "DUPLICATE_LABELS": (2, 4),
    "SCATTERING_TYPE_LABEL_MISMATCH": (1, 4),
    "OCC_ZERO": (0, 4),
    "ISOLATED_H": (0, 4),
    "SYMMETRY_POSSIBLY_MISSED": (3, 4),
}

#: Measured, but the check's definition has since changed enough that the old
#: figure no longer describes it. NON_INTEGRAL_CELL_CONTENT scored 0/10 when it
#: fired on every non-integral cell; since 354e9fd it fires only when there is no
#: formula to compare against, and it has not been drawn since.
REDEFINED = frozenset(("NON_INTEGRAL_CELL_CONTENT",))


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return (None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def tier_of(code: str) -> str:
    if code in REDEFINED or code not in PRECISION:
        return UNMEASURED
    real, judged = PRECISION[code]
    if judged < MIN_JUDGED:
        return UNMEASURED
    return RELIABLE if real / judged >= RELIABLE_MIN else PROMPT


def track_record(code: str) -> Optional[Dict[str, Any]]:
    """The sampled record behind a check's tier, or None if never sampled."""
    if code not in PRECISION:
        return None
    real, judged = PRECISION[code]
    lo, hi = wilson(real, judged)
    out = {"real": real, "judged": judged,
           "precision": round(real / judged, 3) if judged else None,
           "ci95": [round(lo, 3), round(hi, 3)] if judged else None,
           "basis": "model adjudication, precision rounds 1-4"}
    if code in REDEFINED:
        out["note"] = "measured before the check was redefined"
    return out
