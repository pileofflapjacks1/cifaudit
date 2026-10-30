"""Command line entry point: cifaudit <file.cif> [...]"""

import argparse
import json
import sys

from . import __version__
from .audit import DEFAULTS, audit
from .load import CifError, datablock_names
from .model import ERROR, WARN


def _text_report(d):
    out = []
    out.append("datablock: %s   status: %s" % (d["datablock"], d["status"]))
    c = d["counts"]
    out.append("findings:  %d error, %d warn, %d info"
               % (c["error"], c["warn"], c["info"]))
    out.append("")
    if not d["findings"]:
        out.append("  no findings from the checks that ran")
    for f in d["findings"]:
        cc = ", ".join(x["code"] for x in f["checkcif"])
        out.append("  [%-5s] %s" % (f["severity"], f["code"]))
        out.append("          %s" % f["message"])
        if cc:
            out.append("          checkCIF: %s" % cc)
        if f["sites"]:
            out.append("          sites: %s" % ", ".join(f["sites"][:8]))
    out.append("")
    out.append("NOT TESTED (%d areas):" % len(d["not_tested"]))
    for n in d["not_tested"]:
        out.append("  - %s" % n["area"])
        out.append("      %s" % n["reason"])
    out.append("")
    for note in d["notes"]:
        out.append("note: %s" % note)
    out.append("")
    out.append(d["disclaimer"])
    return "\n".join(out)


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="cifaudit",
        description="Offline consistency checks on a crystallographic CIF. "
                    "Not checkCIF; run checkCIF before publishing.")
    p.add_argument("files", nargs="+", help="CIF file(s)")
    p.add_argument("--datablock", help="datablock to audit (default: the first)")
    p.add_argument("--json", action="store_true", help="emit JSON")
    p.add_argument("--list-datablocks", action="store_true")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="override a threshold; repeatable. %s"
                        % ", ".join(sorted(DEFAULTS)))
    p.add_argument("--fail-on", choices=["error", "warn", "never"],
                   default="error",
                   help="exit non-zero when a finding at this level or above is "
                        "present (default: error)")
    p.add_argument("--version", action="version", version=__version__)
    a = p.parse_args(argv)

    options = {}
    for item in a.set:
        if "=" not in item:
            p.error("--set expects KEY=VALUE, got %r" % item)
        k, v = item.split("=", 1)
        if k not in DEFAULTS:
            p.error("unknown threshold %r; known: %s"
                    % (k, ", ".join(sorted(DEFAULTS))))
        try:
            options[k] = type(DEFAULTS[k])(v)
        except ValueError:
            p.error("value for %s must be %s" % (k, type(DEFAULTS[k]).__name__))

    worst = 2
    results = []
    for path in a.files:
        try:
            with open(path, "r") as fh:
                text = fh.read()
        except OSError as exc:
            print("%s: cannot read: %s" % (path, exc), file=sys.stderr)
            return 2
        try:
            if a.list_datablocks:
                names = datablock_names(text)
                if a.json:
                    results.append({"file": path, "datablocks": names})
                else:
                    print("%s: %s" % (path, ", ".join(names)))
                continue
            d = audit(text, a.datablock, options or None).to_dict()
        except CifError as exc:
            if a.json:
                results.append({"file": path, "error": str(exc)})
            else:
                print("%s: %s" % (path, exc), file=sys.stderr)
            worst = min(worst, 2)
            continue
        d["file"] = path
        results.append(d)
        if not a.json:
            if len(a.files) > 1:
                print("=" * 72)
                print(path)
                print("=" * 72)
            print(_text_report(d))
        c = d["counts"]
        if c[ERROR]:
            worst = min(worst, 1)
        elif c[WARN]:
            worst = min(worst, 1 if a.fail_on == "warn" else 0)
        else:
            worst = min(worst, 0)

    if a.json:
        print(json.dumps(results if len(results) != 1 else results[0], indent=2))
    if a.fail_on == "never":
        return 0
    return 0 if worst == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
