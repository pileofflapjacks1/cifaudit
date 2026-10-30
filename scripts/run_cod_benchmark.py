#!/usr/bin/env python3
"""Run cifaudit over a directory of real CIFs and report per-check firing rates.

The point is the FALSE-POSITIVE rate. Deposited, peer-reviewed structures should
mostly come back clean; a check that fires on most of them is measuring
chemistry, not defects.

    python3 scripts/run_cod_benchmark.py <dir-of-cifs> [--json out.json]

Writes a markdown table to stdout. Pass a directory of pristine deposited CIFs,
not fragments or tool output.
"""

import argparse
import collections
import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from cifaudit import audit, datablock_names          # noqa: E402
from cifaudit.load import CifError                   # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--json", dest="json_out")
    ap.add_argument("--glob", default="*.cif")
    ap.add_argument("--recursive", action="store_true",
                    help="walk subdirectories too")
    ap.add_argument("--limit", type=int, help="stop after N files")
    a = ap.parse_args(argv)

    if a.recursive:
        files = []
        for root, _dirs, names in os.walk(a.directory):
            files.extend(os.path.join(root, n) for n in names
                         if n.endswith(".cif"))
        files.sort()
    else:
        files = sorted(glob.glob(os.path.join(a.directory, a.glob)))
    if a.limit:
        files = files[:a.limit]
    if not files:
        print("no files matching %s in %s" % (a.glob, a.directory), file=sys.stderr)
        return 2

    tally = collections.Counter()
    sev = collections.Counter()
    per_file = []
    clean = crashed = unparseable = 0
    t0 = time.time()

    for path in files:
        name = os.path.basename(path)
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        rec = {"file": name}
        try:
            rec["datablocks"] = len(datablock_names(text))
            d = audit(text).to_dict()
        except CifError as exc:
            unparseable += 1
            rec["error"] = str(exc)
            per_file.append(rec)
            continue
        except Exception as exc:                      # a crash is a result
            crashed += 1
            rec["crash"] = "%s: %s" % (type(exc).__name__, exc)
            per_file.append(rec)
            continue
        found = sorted({f["code"] for f in d["findings"]})
        if not found:
            clean += 1
        for code in found:
            tally[code] += 1
        for f in d["findings"]:
            sev[f["severity"]] += 1
        rec.update(status=d["status"], counts=d["counts"], codes=found)
        per_file.append(rec)

    n = len(files)
    audited = n - crashed - unparseable
    elapsed = time.time() - t0

    print("# cifaudit over %d real CIFs in `%s`" % (n, a.directory))
    print()
    print("- audited: **%d**   unparseable: %d   crashed: **%d**" %
          (audited, unparseable, crashed))
    print("- files with zero findings: **%d/%d (%.0f%%)**" %
          (clean, audited, 100.0 * clean / audited if audited else 0))
    print("- findings by severity: %s" % (dict(sev) or "none"))
    print("- %.1fs total, %.2fs per file" % (elapsed, elapsed / n))
    print()
    print("| check | files | rate |")
    print("|---|---|---|")
    for code, k in tally.most_common():
        print("| `%s` | %d/%d | %.0f%% |" %
              (code, k, audited, 100.0 * k / audited if audited else 0))
    if not tally:
        print("| _none fired_ | 0 | 0% |")

    if a.json_out:
        with open(a.json_out, "w") as fh:
            json.dump({"directory": a.directory, "n": n, "audited": audited,
                       "clean": clean, "crashed": crashed,
                       "unparseable": unparseable, "seconds": elapsed,
                       "rates": dict(tally), "severity": dict(sev),
                       "files": per_file}, fh, indent=2)
        print("\nwrote %s" % a.json_out, file=sys.stderr)
    return 1 if crashed else 0


if __name__ == "__main__":
    sys.exit(main())
