"""cifaudit - offline, deterministic checks on a crystallographic CIF.

cifaudit is NOT checkCIF. It implements a small subset of the checkCIF/PLATON
tests, runs entirely offline, reads no free-text fields, puts no language model in
the validation path, and reports in every response which checkCIF tests it did not
perform. Run checkCIF (https://checkcif.iucr.org/) before publishing.
"""

from .audit import DEFAULTS, audit
from .load import CifError, datablock_names, load
from .model import AuditResult, ERROR, Finding, INFO, NotTested, WARN

__version__ = "0.1.0"
__all__ = ["audit", "load", "datablock_names", "DEFAULTS", "AuditResult",
           "Finding", "NotTested", "CifError", "ERROR", "WARN", "INFO",
           "__version__"]
