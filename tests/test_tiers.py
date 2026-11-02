"""Tiers: each structural check's sampled track record, and what follows from it."""

import json
import pathlib
from collections import Counter, defaultdict

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"

ROUNDS = (
    ("precision_sample_items.json", "precision_model_verdicts.json"),
    ("precision_round2_sample_items.json", "precision_round2_model_verdicts.json"),
    ("precision_round3_sample_items.json", "precision_round3_model_verdicts.json"),
    ("precision_round4_sample_items.json", "precision_round4_model_verdicts.json"),
)


def test_track_records_match_the_published_verdicts():
    """Every number in tiers.PRECISION is recomputed from the four published
    sample and verdict files, so the tiers cannot drift from the evidence."""
    from cifaudit.tiers import PRECISION
    pool = defaultdict(Counter)
    for items_file, verdicts_file in ROUNDS:
        items = json.loads((DATA / items_file).read_text())["items"]
        verdicts = {v["id"]: v["verdict"]
                    for v in json.loads((DATA / verdicts_file).read_text())}
        assert {i["id"] for i in items} == set(verdicts)
        for i in items:
            pool[i["code"]][verdicts[i["id"]]] += 1
    recomputed = {code: (c["real"], c["real"] + c["false"])
                  for code, c in pool.items()}
    assert recomputed == PRECISION


def test_tier_rule():
    from cifaudit.tiers import PROMPT, RELIABLE, UNMEASURED, tier_of
    assert tier_of("ADP_EXTREME_ANISOTROPY") == RELIABLE     # 14/16
    assert tier_of("FORMULA_SITES_MISMATCH") == RELIABLE     # 24/40
    assert tier_of("NO_HYDROGENS_WITH_CARBON") == RELIABLE   # 6/8
    assert tier_of("ISOLATED_ATOM") == PROMPT                # 4/40
    assert tier_of("DISORDER_SUM_NOT_ONE") == PROMPT         # 0/12
    assert tier_of("DENSITY_BELOW_ONE") == PROMPT            # 2/8
    assert tier_of("SYMMETRY_POSSIBLY_MISSED") == UNMEASURED  # 3/4: too few
    assert tier_of("NON_INTEGRAL_CELL_CONTENT") == UNMEASURED  # redefined
    assert tier_of("ADP_NPD") == UNMEASURED                  # never drawn


def test_structural_findings_carry_a_tier_and_reporting_ones_do_not(run):
    r = run("h_across_centre")
    assert r["findings"]
    for f in r["findings"]:
        if f["kind"] == "structural":
            assert f["tier"] in ("reliable", "prompt", "unmeasured")
        else:
            assert f["tier"] is None and f["track_record"] is None
    assert sum(r["tier_counts"].values()) == sum(r["counts"].values())


def _clash(tmp_path, cif):
    """A file whose only structural error is ATOMS_TOO_CLOSE (prompt tier)."""
    text = cif("h_across_centre").replace("0.500 0.0500", "1.000 0.0500")
    path = tmp_path / "clash.cif"
    path.write_text(text)
    return str(path)


def test_a_prompt_tier_error_does_not_fail_the_run_by_default(tmp_path, cif,
                                                              capsys):
    from cifaudit.cli import main
    path = _clash(tmp_path, cif)
    assert main([path]) == 0
    out = capsys.readouterr().out
    assert "PROMPTS TO LOOK" in out and "ATOMS_TOO_CLOSE" in out


def test_fail_on_prompts_restores_it(tmp_path, cif):
    from cifaudit.cli import main
    assert main(["--fail-on-prompts", _clash(tmp_path, cif)]) == 1


def test_a_reliable_tier_warning_still_fails_under_fail_on_warn(tmp_path, cif):
    """FORMULA_SITES_MISMATCH is reliable (24/40)."""
    from cifaudit.cli import main
    path = tmp_path / "f.cif"
    path.write_text(cif("h_across_centre").replace("'H1'", "'H3'"))
    assert main(["--fail-on", "warn", str(path)]) == 1


def test_mcp_audit_explains_tiers_and_coverage_lists_them(cif):
    from cifaudit.tools import audit_cif, coverage
    r = audit_cif(cif_text=cif("h_across_centre"))
    assert "prompt" in r["tier_note"]
    t = coverage()["tiers"]
    assert t["ISOLATED_ATOM"] == {"tier": "prompt", "real": 4, "judged": 40}
