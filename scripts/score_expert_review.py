#!/usr/bin/env python3
"""Compare the expert's marks with the model adjudicator's verdicts.

    python scripts/score_expert_review.py <marks_dir_or_json>

<marks> is either a JSON array of {"id": "r3-R08", "verdict": ..., "note": ...}
or a directory of one JSON file per mark (as the review page's `marks`
collection exports, file name = id). Ids are "r<round>-<item>"; the model's
verdicts come from data/processed. Reports raw agreement and Cohen's kappa over
real / false / unsure, precision by each judge on the same items, and every
disagreement with both reasons. See data/processed/expert_review_preregistration.md.
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from score_precision import kappa, wilson  # noqa: E402

DATA = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed"
VERDICTS = {1: "precision_model_verdicts.json", 2: "precision_round2_model_verdicts.json",
            3: "precision_round3_model_verdicts.json",
            4: "precision_round4_model_verdicts.json"}
ITEMS = {1: "precision_sample_items.json", 2: "precision_round2_sample_items.json",
         3: "precision_round3_sample_items.json",
         4: "precision_round4_sample_items.json"}


def load_marks(src):
    p = pathlib.Path(src)
    if p.is_dir():
        out = []
        for f in sorted(p.rglob("*.json")):
            d = json.loads(f.read_text())
            d = d.get("data", d)
            d.setdefault("id", f.stem)
            out.append(d)
        return out
    return json.loads(p.read_text())


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    marks = [m for m in load_marks(argv[0])
             if m.get("verdict") in ("real", "false", "unsure")]
    model, code = {}, {}
    for r in VERDICTS:
        for v in json.loads((DATA / VERDICTS[r]).read_text()):
            model["r%d-%s" % (r, v["id"])] = v
        for it in json.loads((DATA / ITEMS[r]).read_text())["items"]:
            code["r%d-%s" % (r, it["id"])] = it["code"]
    pairs = [(m, model[m["id"]]) for m in marks if m["id"] in model]
    if not pairs:
        print("no marked items")
        return 1
    a = [m["verdict"] for m, _ in pairs]
    b = [v["verdict"] for _, v in pairs]
    po, k = kappa(a, b)
    print("### Expert vs model adjudicator (%d items marked)\n" % len(pairs))
    print("- raw agreement: %.0f%% (%d/%d)" % (100 * po, sum(x == y for x, y in zip(a, b)), len(a)))
    print("- Cohen's kappa: %s" % ("undefined" if k is None else "%.2f" % k))
    for name, vs in (("expert", a), ("model", b)):
        r_, f_ = vs.count("real"), vs.count("false")
        p, lo, hi = wilson(r_, r_ + f_)
        print("- precision by %s: %s" % (name, "—" if p is None else
              "%d/%d · %.0f%% (%.0f–%.0f%%)" % (r_, r_ + f_, 100 * p, 100 * lo, 100 * hi)))
    dis = [(m, v) for m, v in pairs if m["verdict"] != v["verdict"]]
    print("\n### Disagreements (%d)\n" % len(dis))
    print("| item | check | expert | model | expert note | model reason |")
    print("|---|---|---|---|---|---|")
    for m, v in dis:
        print("| %s | `%s` | %s | %s | %s | %s |" % (
            m["id"], code.get(m["id"], "?"), m["verdict"], v["verdict"],
            (m.get("note") or "").replace("|", "/")[:200],
            v.get("reason", "").replace("|", "/")[:200]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
