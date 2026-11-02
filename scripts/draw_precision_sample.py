#!/usr/bin/env python3
"""Draw the precision-study sample exactly as pre-registered.

See data/processed/precision_preregistration.md. Structural findings only, one
per (file, check), 10 each for the three suspicious checks, then up to 4 each for
every other firing check in descending firing rate until 60, seed 20261031.

    python scripts/draw_precision_sample.py <dir-of-cifs> <out.json> [--seed N]

Round 1 used the default seed 20261031; round 2 used --seed 20261101.
"""

import json
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from cifaudit import audit                      # noqa: E402
from cifaudit.load import CifError              # noqa: E402

SEED = 20261031
TARGET = 60
PRIORITY = ("NON_INTEGRAL_CELL_CONTENT", "FORMULA_SITES_MISMATCH", "ISOLATED_ATOM")
PRIORITY_N = 10
OTHER_N = 4

CONTEXT_TAGS = ("_cell_length_a", "_cell_length_b", "_cell_length_c",
                "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma",
                "_cell_volume", "_cell_formula_units_Z", "_chemical_formula_sum",
                "_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt",
                "_diffrn_radiation_probe", "_atom_sites_solution_hydrogens")


def context(text, sites, full_loop):
    """The CIF lines an expert needs: cell, formula, Z, symmetry, and the
    atom_site rows for the sites involved (or the whole loop, capped)."""
    lines = text.splitlines()
    head = [l for l in lines if l.split()[:1] and l.split()[0] in CONTEXT_TAGS]
    rows, header, in_loop = [], [], False
    for i, l in enumerate(lines):
        s = l.strip()
        if s == "loop_":
            nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
            in_loop = nxt.startswith("_atom_site_label")
            header = []
            continue
        if in_loop:
            if s.startswith("_atom_site_"):
                header.append(s)
                continue
            if not s or s.startswith("_") or s.startswith("loop_") or s.startswith("#"):
                in_loop = False
                continue
            tok = s.split()[0]
            if full_loop or tok in sites:
                rows.append(s)
    cap = 45
    more = max(0, len(rows) - cap)
    return {"header_items": head, "atom_site_columns": header,
            "atom_site_rows": rows[:cap],
            "atom_site_rows_omitted": more}


def main(src, out, seed=SEED):
    files = sorted(os.path.join(r, n) for r, _d, ns in os.walk(src)
                   for n in ns if n.endswith(".cif"))
    pool = defaultdict(list)          # code -> [item]
    fired = Counter()
    for path in files:
        try:
            text = open(path, encoding="utf-8", errors="replace").read()
            d = audit(text).to_dict()
        except (CifError, Exception):
            continue
        seen = set()
        for f in d["findings"]:
            if f["kind"] != "structural" or f["code"] in seen:
                continue
            seen.add(f["code"])
            fired[f["code"]] += 1
            pool[f["code"]].append({"path": path, "finding": f, "text": text})

    rng = random.Random(seed)
    chosen = []
    for code in PRIORITY:
        items = sorted(pool.get(code, []), key=lambda x: x["path"])
        chosen += [(code, x) for x in rng.sample(items, min(PRIORITY_N, len(items)))]
    for code, _n in sorted(fired.items(), key=lambda kv: (-kv[1], kv[0])):
        if code in PRIORITY or len(chosen) >= TARGET:
            continue
        items = sorted(pool[code], key=lambda x: x["path"])
        k = min(OTHER_N, len(items), TARGET - len(chosen))
        chosen += [(code, x) for x in rng.sample(items, k)]

    out_items = []
    for i, (code, x) in enumerate(chosen, 1):
        f = x["finding"]
        full = code in ("FORMULA_SITES_MISMATCH", "NON_INTEGRAL_CELL_CONTENT")
        out_items.append({
            "id": "P%02d" % i,
            "cod_id": os.path.basename(x["path"])[:-4],
            "code": code, "severity": f["severity"],
            "message": f["message"], "evidence": f["evidence"],
            "sites": f["sites"],
            "checkcif": [c["code"] for c in f["checkcif"]],
            "context": context(x["text"], set(f["sites"]), full),
        })
    json.dump({"seed": seed, "files_audited": len(files),
               "firing_files": dict(fired), "items": out_items},
              open(out, "w"), indent=1)
    print("files: %d | items drawn: %d" % (len(files), len(out_items)))
    print("per check:", dict(Counter(i["code"] for i in out_items)))
    print("available pool:", dict(fired))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=SEED)
    a = ap.parse_args()
    main(a.src, a.out, a.seed)
