#!/usr/bin/env python3
"""Score the precision study, per check and per adjudicator.

    python scripts/score_precision.py <sample.json> --model <verdicts.json> [--human <verdicts.json>]

Each verdicts file is a JSON array of {"id", "verdict": real|false|unsure, ...}.
Precision = real / (real + false), Wilson 95% interval; unsure is reported and
excluded. Where both adjudicators judged an item, Cohen's kappa over the three
verdict classes and raw agreement are reported.
"""

import argparse
import json
import math
from collections import Counter, defaultdict


def wilson(k, n, z=1.96):
    if n == 0:
        return (None, None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, max(0.0, c - h), min(1.0, c + h))


def kappa(a, b):
    cats = sorted(set(a) | set(b))
    n = len(a)
    po = sum(1 for x, y in zip(a, b) if x == y) / n
    ca, cb = Counter(a), Counter(b)
    pe = sum((ca[c] / n) * (cb[c] / n) for c in cats)
    return (po, None if pe == 1 else (po - pe) / (1 - pe))


def table(label, items, verdicts):
    by = defaultdict(Counter)
    for it in items:
        v = verdicts.get(it["id"])
        if v:
            by[it["code"]][v] += 1
    rows, tot = [], Counter()
    for code in sorted(by, key=lambda c: -sum(by[c].values())):
        c = by[code]
        tot.update(c)
        p, lo, hi = wilson(c["real"], c["real"] + c["false"])
        rows.append((code, c["real"], c["false"], c["unsure"], p, lo, hi))
    p, lo, hi = wilson(tot["real"], tot["real"] + tot["false"])
    rows.append(("ALL", tot["real"], tot["false"], tot["unsure"], p, lo, hi))
    print("\n### %s\n" % label)
    print("| check | real | false | unsure | precision | 95% CI |")
    print("|---|---|---|---|---|---|")
    for code, r, f, u, p, lo, hi in rows:
        ps = "—" if p is None else "%.0f%%" % (100 * p)
        ci = "—" if p is None else "%.0f–%.0f%%" % (100 * lo, 100 * hi)
        name = "**all**" if code == "ALL" else "`%s`" % code
        print("| %s | %d | %d | %d | %s | %s |" % (name, r, f, u, ps, ci))


def load(path):
    if not path:
        return {}
    return {r["id"]: r["verdict"] for r in json.load(open(path))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample")
    ap.add_argument("--model")
    ap.add_argument("--human")
    a = ap.parse_args()
    items = json.load(open(a.sample))["items"]
    model, human = load(a.model), load(a.human)
    if model:
        table("Model adjudicator (%d of %d judged)" % (len(model), len(items)), items, model)
    if human:
        table("Human expert (%d of %d judged)" % (len(human), len(items)), items, human)
    both = [i["id"] for i in items if i["id"] in model and i["id"] in human]
    if both:
        po, k = kappa([model[i] for i in both], [human[i] for i in both])
        print("\nAgreement on %d items judged by both: raw %.0f%%, Cohen's kappa %s"
              % (len(both), 100 * po, "n/a" if k is None else "%.2f" % k))


if __name__ == "__main__":
    main()
