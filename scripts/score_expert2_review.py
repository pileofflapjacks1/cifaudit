#!/usr/bin/env python3
"""Score the second expert review: rater 1 vs rater 2 vs the model.

    python scripts/score_expert2_review.py <rater2 marks dir or json>

<marks> is a JSON array of {"id", "verdict", "note"?, "ai"?} or a directory of
one JSON file per mark (the second review page's `marks` collection, exported).
Rater 1 is the first review's second pass (data/processed/expert_review_marks.json).
Applies the decision rule in data/processed/expert2_review_preregistration.md.
"""

import json
import pathlib
import random
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from score_expert_review import load_marks  # noqa: E402
from score_precision import kappa, wilson  # noqa: E402

DATA = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed"
RESAMPLES, SEED = 4000, 1


def boot(a, b):
    rng = random.Random(SEED)
    ks = []
    for _ in range(RESAMPLES):
        s = [rng.randrange(len(a)) for _ in a]
        k = kappa([a[j] for j in s], [b[j] for j in s])[1]
        if k is not None:
            ks.append(k)
    ks.sort()
    return ks[int(0.025 * len(ks))], ks[int(0.975 * len(ks))]


def pair(label, a, b):
    po, k = kappa(a, b)
    lo, hi = boot(a, b)
    both = sum(1 for x, y in zip(a, b) if x == y == "real")
    print("| %s | %d | %d (%.0f%%) | %+.2f (%+.2f to %+.2f) | %d |"
          % (label, len(a), sum(x == y for x, y in zip(a, b)), 100 * po,
             k if k is not None else float("nan"), lo, hi, both))
    return k, lo, hi


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    r2 = {m["id"]: m for m in load_marks(argv[0])
          if m.get("verdict") in ("real", "false", "unsure") and m["id"] != "probe"}
    r1 = {m["id"]: m for m in json.loads((DATA / "expert_review_marks.json").read_text())}
    ids = sorted(i for i in r2 if i in r1)
    if len(ids) < len(r1):
        print("note: rater 2 marked %d of %d items\n" % (len(ids), len(r1)))
    for label, keep in (("all marks", lambda i: True),
                        ("without AI-assisted marks (either rater)",
                         lambda i: not r2[i].get("ai") and not r1[i]["ai_assisted_pass2"])):
        sel = [i for i in ids if keep(i)]
        if len(sel) < 2:
            continue
        a1 = [r1[i]["expert_pass2"] for i in sel]
        a2 = [r2[i]["verdict"] for i in sel]
        m = [r1[i]["model"] for i in sel]
        print("### %s\n" % label)
        print("| pair | n | agree | Cohen's kappa (bootstrap 95%) | both real |")
        print("|---|---|---|---|---|")
        ee = pair("rater 1 vs rater 2", a1, a2)
        e1 = pair("rater 1 vs model", a1, m)
        e2 = pair("rater 2 vs model", a2, m)
        for name, v in (("rater 1", a1), ("rater 2", a2), ("model", m)):
            r, f = v.count("real"), v.count("false")
            p, lo, hi = wilson(r, r + f)
            print("\n- precision by %s: %s" % (name, "—" if p is None else
                  "%d/%d = %.0f%% (%.0f-%.0f%%)" % (r, r + f, 100 * p, 100 * lo, 100 * hi)), end="")
        print("\n")
        if label == "all marks":
            if ee[1] > max(e1[2], e2[2]):
                verdict = "experts agree with each other more than with the model"
            elif ee[0] is not None and ee[0] < 0.40 and (ee[1] <= e1[2] or ee[1] <= e2[2]):
                verdict = "experts agree no better with each other than with the model"
            else:
                verdict = "undetermined at this n"
            print("**Pre-registered decision rule: %s.**\n" % verdict)
    print("### Items by verdict (rater 1 / rater 2 / model)\n")
    print("| item | rater 1 | rater 2 | model | rater 2 note |")
    print("|---|---|---|---|---|")
    for i in ids:
        print("| %s | %s | %s | %s | %s |" % (i, r1[i]["expert_pass2"], r2[i]["verdict"],
              r1[i]["model"], (r2[i].get("note") or "").replace("|", "/")[:160]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
