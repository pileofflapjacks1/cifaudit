"""The response contract: abstention, not_tested, alert mapping, no echo."""

import json

import pytest

from conftest import codes


def test_every_response_lists_what_was_not_tested(run):
    for name in ("clean", "npd", "no_cell", "triclinic"):
        r = run(name)
        assert r["not_tested"], name
        for n in r["not_tested"]:
            assert n["area"] and n["reason"]


def test_not_tested_entries_name_real_plat_codes(run):
    from cifaudit.alerts import PLAT
    for n in run("clean")["not_tested"]:
        for c in n["checkcif_not_covered"]:
            assert c["code"] in PLAT, c["code"]
            assert c["description"]


def test_findings_carry_a_real_checkcif_code(run):
    from cifaudit.alerts import PLAT
    for name in ("occ_gt_one", "npd", "bad_volume", "short_contact", "no_h"):
        for f in run(name)["findings"]:
            assert f["checkcif"], f["code"]
            for c in f["checkcif"]:
                assert c["code"] in PLAT
                assert c["description"] is not None


def test_no_alert_level_is_asserted(run):
    """checkCIF levels are A/B/C/G; the public PLATON list gives Type_N and
    publishes no mapping, so we must never emit an A/B/C/G level."""
    blob = json.dumps(run("npd"))
    assert '"level"' not in blob
    for f in run("npd")["findings"]:
        for c in f["checkcif"]:
            assert c["platon_type"] is None or c["platon_type"].startswith("Type_")


def test_disclaimer_present_and_defers_to_checkcif(run):
    d = run("clean")["disclaimer"]
    assert "not checkCIF" in d
    assert "checkcif.iucr.org" in d


def test_every_finding_cites_the_cif_items_it_used(run):
    for f in run("occ_gt_one")["findings"]:
        assert f["based_on"], f["code"]
        assert all(t.startswith("_") for t in f["based_on"])


def test_multiblock_is_announced_not_silently_picked(run, cif):
    from cifaudit import datablock_names
    r = run("multiblock")
    assert datablock_names(cif("multiblock")) == ["clean", "second"]
    assert any("2 datablocks" in n for n in r["notes"])


def test_named_datablock_is_honoured(run):
    assert run("multiblock", datablock="second")["datablock"] == "second"


def test_unknown_datablock_raises(cif):
    from cifaudit import audit
    from cifaudit.load import CifError
    with pytest.raises(CifError):
        audit(cif("multiblock"), datablock="nope")


def test_thresholds_are_disclosed(run):
    assert any("Thresholds used" in n for n in run("clean")["notes"])


def test_unknown_option_is_refused(cif):
    from cifaudit import audit
    from cifaudit.load import CifError
    with pytest.raises(CifError):
        audit(cif("clean"), options={"not_a_real_option": 1})


def test_thresholds_are_overridable(cif):
    """Loosening the volume tolerance must silence the volume finding.

    Note: raising nonbonded_min does NOT create contacts in the clean fixture,
    because every short pair there is 1,3 geminal and so excluded by topology
    rather than by threshold.
    """
    from cifaudit import audit
    default = audit(cif("bad_volume")).to_dict()
    assert "CELL_VOLUME_MISMATCH" in {f["code"] for f in default["findings"]}
    loose = audit(cif("bad_volume"), options={"volume_rel_tol": 1.0}).to_dict()
    assert "CELL_VOLUME_MISMATCH" not in {f["code"] for f in loose["findings"]}
