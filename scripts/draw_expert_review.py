#!/usr/bin/env python3
"""Draw the items for the expert review, reproducibly.

    python scripts/draw_expert_review.py [--n 20] [--seed 20261104]

A simple random sample, without stratification, from every finding drawn in
precision rounds 1-4 (data/processed/precision*_sample_items.json). Stratifying
by the model's verdict would distort agreement statistics, which depend on how
often each verdict occurs. Prints the selected round/id pairs as JSON.
"""

import argparse
import json
import pathlib
import random

DATA = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed"
ROUNDS = ("precision_sample_items.json", "precision_round2_sample_items.json",
          "precision_round3_sample_items.json",
          "precision_round4_sample_items.json")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20261104)
    a = ap.parse_args(argv)
    pool = []
    for r, name in enumerate(ROUNDS, 1):
        for it in json.loads((DATA / name).read_text())["items"]:
            pool.append({"round": r, "id": it["id"], "cod_id": it["cod_id"],
                         "code": it["code"]})
    pool.sort(key=lambda x: (x["round"], x["id"]))
    pick = random.Random(a.seed).sample(pool, a.n)
    pick.sort(key=lambda x: (x["round"], x["id"]))
    print(json.dumps({"seed": a.seed, "pool": len(pool), "items": pick}, indent=1))


if __name__ == "__main__":
    main()
