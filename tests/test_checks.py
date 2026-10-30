"""Each planted defect must be caught, and the clean reference must be silent."""

import pytest

from conftest import codes


def test_clean_has_no_findings(run):
    r = run("clean")
    assert r["status"] == "ok"
    assert r["findings"] == [], [f["code"] for f in r["findings"]]
    assert r["counts"] == {"error": 0, "warn": 0, "info": 0}


def test_clean_ran_every_check(run):
    r = run("clean")
    assert set(r["checked"]) == {
        "cell", "labels", "occupancy", "adp", "symmetry", "contacts",
        "hydrogens", "fragments", "formula"}


@pytest.mark.parametrize("fixture,code,severity", [
    ("occ_gt_one", "OCC_GT_ONE", "error"),
    ("occ_zero", "OCC_ZERO", "warn"),
    ("npd", "ADP_NPD", "error"),
    ("dup_labels", "DUPLICATE_LABELS", "error"),
    ("bad_volume", "CELL_VOLUME_MISMATCH", "error"),
    ("short_contact", "ATOMS_TOO_CLOSE", "error"),
    ("no_h", "NO_HYDROGENS_WITH_CARBON", "warn"),
    ("no_h", "H_TREATMENT_WITHOUT_H", "error"),
    ("disorder_bad", "DISORDER_SUM_NOT_ONE", "warn"),
])
def test_defect_is_caught(run, fixture, code, severity):
    r = run(fixture)
    hit = [f for f in r["findings"] if f["code"] == code]
    assert hit, "%s not found in %s; got %s" % (code, fixture, codes(r))
    assert hit[0]["severity"] == severity


def test_npd_eigenvalue_reported(run):
    r = run("npd")
    f = [x for x in r["findings"] if x["code"] == "ADP_NPD"][0]
    assert min(f["evidence"]["eigenvalues"]) < 0
    assert f["sites"] == ["C1"]


def test_disorder_summing_to_one_is_not_flagged(run):
    assert "DISORDER_SUM_NOT_ONE" not in codes(run("disorder_ok"))


def test_geminal_pairs_are_not_short_contacts(run):
    """H1...O1 is 1.995 A and H1...H2 is 1.861 A in the clean fixture.

    Both are below the thresholds but both are 1,3 (geminal) through C1, which
    PLATON also excludes. Treating them as contacts was a real bug.
    """
    assert "SHORT_NONBONDED_CONTACT" not in codes(run("clean"))
    assert "SHORT_HH_CONTACT" not in codes(run("clean"))


def test_triclinic_cell_produces_no_spurious_geometry_errors(run):
    r = run("triclinic")
    assert r["counts"]["error"] == 0, codes(r)


def test_no_cell_is_not_enough_data(run):
    r = run("no_cell")
    assert r["status"] == "not_enough_data"
    assert "CELL_INCOMPLETE" in codes(r)


def test_findings_are_sorted_most_severe_first(run):
    order = {"error": 0, "warn": 1, "info": 2}
    sev = [order[f["severity"]] for f in run("injection")["findings"]]
    assert sev == sorted(sev)
