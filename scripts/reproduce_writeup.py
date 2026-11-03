#!/usr/bin/env python3
"""Recompute every figure in the write-up from the published files.

    python scripts/reproduce_writeup.py

Reads only data/processed. Prints each number the write-up states, under the
section that states it, so a reader can check the paper line by line. Exits 1
if any figure fails an internal consistency check.
"""

import json
import pathlib
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from score_precision import kappa, wilson  # noqa: E402

DATA = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed"
ROUNDS = {
    1: ("7b2dbbc", "precision_sample_items.json", "precision_model_verdicts.json"),
    2: ("40fdaec", "precision_round2_sample_items.json", "precision_round2_model_verdicts.json"),
    3: ("33f30cb", "precision_round3_sample_items.json", "precision_round3_model_verdicts.json"),
    4: ("50f1789", "precision_round4_sample_items.json", "precision_round4_model_verdicts.json"),
}
HOLDOUT = (
    ("held-out, as registered", "cod_holdout_700.json"),
    ("+ bug fixes", "cod_holdout_700_fixed.json"),
    ("+ conventions", "cod_holdout_700_conventions.json"),
    ("+ round-2 fixes", "cod_holdout_700_round2_fixes.json"),
    ("+ round-3 fixes", "cod_holdout_700_round3_fixes.json"),
    ("tiers, excluding prompts", "cod_holdout_700_tiers.json"),
)
BOOT_RESAMPLES = 4000
BOOT_SEED = 1
failures = []


def load(name):
    return json.loads((DATA / name).read_text())


def pct(k, n):
    p, lo, hi = wilson(k, n)
    return "%d/%d = %.0f%% (95%% CI %.0f-%.0f%%)" % (k, n, 100 * p, 100 * lo, 100 * hi)


def check(cond, what):
    if not cond:
        failures.append(what)
        print("  !! CHECK FAILED: %s" % what)


def boot_kappa(a, b):
    rng = random.Random(BOOT_SEED)
    ks = []
    for _ in range(BOOT_RESAMPLES):
        s = [rng.randrange(len(a)) for _ in a]
        k = kappa([a[j] for j in s], [b[j] for j in s])[1]
        if k is not None:
            ks.append(k)
    ks.sort()
    return ks[int(0.025 * len(ks))], ks[int(0.975 * len(ks))]


def main():
    print("# Figures in the write-up, recomputed from data/processed\n")

    # ---- Abstract and Results: four rounds -----------------------------------
    print("## Precision by round (model adjudicator)")
    pool = defaultdict(Counter)
    model = {}
    code_of = {}
    drawn = judged = real_all = 0
    for r, (version, items_f, verdicts_f) in ROUNDS.items():
        items = load(items_f)["items"]
        verdicts = {v["id"]: v for v in load(verdicts_f)}
        check({i["id"] for i in items} == set(verdicts), "round %d ids match" % r)
        c = Counter(verdicts[i["id"]]["verdict"] for i in items)
        print("- round %d (%s): %s, unsure %d, drawn %d, files audited %s"
              % (r, version, pct(c["real"], c["real"] + c["false"]), c["unsure"],
                 len(items), load(items_f).get("files_audited")))
        drawn += len(items)
        judged += c["real"] + c["false"]
        real_all += c["real"]
        for i in items:
            key = "r%d-%s" % (r, i["id"])
            model[key] = verdicts[i["id"]]
            code_of[key] = i["code"]
            pool[i["code"]][verdicts[i["id"]]["verdict"]] += 1
    print("- all rounds: %s; drawn %d, decided %d, unsure %d"
          % (pct(real_all, judged), drawn, judged, drawn - judged))
    check(drawn == 240 and judged == 237 and real_all == 83, "pooled 83/237 of 240")

    # ---- Results: pooled per check, and the tiers ----------------------------
    print("\n## Pooled precision per check, and the tier it implies")
    from cifaudit.tiers import PRECISION, tier_of  # noqa: E402
    rows = sorted(pool.items(), key=lambda kv: -(kv[1]["real"] / max(1, kv[1]["real"] + kv[1]["false"])))
    for code, c in rows:
        n = c["real"] + c["false"]
        print("- %-32s %s  tier: %s" % (code, pct(c["real"], n), tier_of(code)))
        check(PRECISION.get(code) == (c["real"], n), "tiers.PRECISION matches %s" % code)

    # ---- Results: held-out benchmark ------------------------------------------
    print("\n## Held-out 700: files with no structural finding")
    for label, f in HOLDOUT:
        d = load(f)
        key = ("structural_clean_excluding_prompts" if "excluding" in label
               else "structural_clean")
        print("- %-26s %d/%d = %.0f%%" % (label, d[key], d["audited"], 100 * d[key] / d["audited"]))

    # ---- Results: expert review -----------------------------------------------
    print("\n## Expert review (20 items, seed 20261104)")
    marks = load("expert_review_marks.json")
    check(len(marks) == 20, "20 expert items")
    for label, field, sub in (
            ("second pass, registered question", "expert_pass2", lambda m: True),
            ("second pass, without AI-assisted marks", "expert_pass2",
             lambda m: not m["ai_assisted_pass2"]),
            ("second pass, AI-assisted only", "expert_pass2", lambda m: m["ai_assisted_pass2"]),
            ("first pass, 'condition present'", "expert_pass1", lambda m: True)):
        ms = [m for m in marks if sub(m)]
        a = [m[field] for m in ms]
        b = [m["model"] for m in ms]
        check(all(m["model"] == model[m["id"]]["verdict"] for m in ms),
              "expert file's model verdicts match the round files")
        po, k = kappa(a, b)
        lo, hi = boot_kappa(a, b)
        both = sum(1 for x, y in zip(a, b) if x == y == "real")
        print("- %s: n=%d, agree %d (%.0f%%), kappa %+.2f (bootstrap %+.2f to %+.2f), "
              "expert real %d, model real %d, both real %d"
              % (label, len(ms), sum(x == y for x, y in zip(a, b)), 100 * po, k, lo, hi,
                 a.count("real"), b.count("real"), both))
    p2 = [m["expert_pass2"] for m in marks]
    print("- precision by expert (pass 2): %s" % pct(p2.count("real"), p2.count("real") + p2.count("false")))
    mv = [m["model"] for m in marks]
    print("- precision by model, same items: %s" % pct(mv.count("real"), mv.count("real") + mv.count("false")))
    conf = Counter((m["expert_pass2"], m["model"]) for m in marks)
    print("- pass 2 cross-table (expert, model): %s" % dict(sorted(conf.items())))
    changed = sum(m["expert_pass1"] != m["expert_pass2"] for m in marks)
    print("- marks changed between passes: %d of 20" % changed)
    print("- AI-assisted pass-2 marks: %d" % sum(m["ai_assisted_pass2"] for m in marks))
    for tier in ("reliable", "prompt", "unmeasured"):
        g = [m for m in marks if m["tier_now"] == tier]
        e = sum(m["expert_pass2"] == "real" for m in g)
        mo = sum(m["model"] == "real" for m in g)
        extra = (" (expert %s)" % pct(e, len(g))) if tier == "reliable" else ""
        print("- %-10s items %d: expert real %d, model real %d%s" % (tier, len(g), e, mo, extra))

    print("\n%s" % ("ALL CHECKS PASSED" if not failures else "%d CHECK(S) FAILED" % len(failures)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
    sys.exit(main())
