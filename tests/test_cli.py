"""The CLI, including its exit codes."""

import json
import pathlib

import pytest

from cifaudit.cli import main

DATA = pathlib.Path(__file__).parent / "data"


def p(name):
    return str(DATA / (name + ".cif"))


def test_clean_exits_zero(capsys):
    assert main([p("clean")]) == 0
    assert "no findings from the checks that ran" in capsys.readouterr().out


def test_error_exits_one(capsys):
    assert main([p("npd")]) == 1
    out = capsys.readouterr().out
    assert "ADP_NPD" in out
    assert "PLAT211" in out


def test_fail_on_never_exits_zero(capsys):
    assert main([p("npd"), "--fail-on", "never"]) == 0


def test_fail_on_warn_catches_warnings(capsys):
    # occ_zero has warnings but no errors
    assert main([p("occ_zero")]) == 0
    capsys.readouterr()
    assert main([p("occ_zero"), "--fail-on", "warn"]) == 1


def test_json_output_is_valid(capsys):
    main([p("npd"), "--json"])
    d = json.loads(capsys.readouterr().out)
    assert d["datablock"] == "npd"
    assert d["file"].endswith("npd.cif")
    assert any(f["code"] == "ADP_NPD" for f in d["findings"])


def test_list_datablocks(capsys):
    main([p("multiblock"), "--list-datablocks"])
    assert "clean, second" in capsys.readouterr().out


def test_set_overrides_threshold(capsys):
    assert main([p("bad_volume")]) == 1
    capsys.readouterr()
    assert main([p("bad_volume"), "--set", "volume_rel_tol=1.0"]) == 0


def test_unknown_threshold_is_rejected():
    with pytest.raises(SystemExit):
        main([p("clean"), "--set", "nope=1"])


def test_missing_file_exits_two(capsys):
    assert main(["/nonexistent/path.cif"]) == 2


def test_report_always_shows_not_tested(capsys):
    main([p("clean")])
    assert "NOT TESTED" in capsys.readouterr().out


def test_report_always_shows_disclaimer(capsys):
    main([p("clean")])
    assert "not checkCIF" in capsys.readouterr().out


# ---- exit codes across MULTIPLE files ---------------------------------------
# No test ran more than one file, so four exit-code bugs went unnoticed until
# CI was being written. The worst: a batch with one failing and one clean file
# exited 0, so `cifaudit *.cif` in CI passed as long as any one file was clean.

def test_a_batch_fails_if_ANY_file_fails(capsys):
    assert main([p("npd"), p("clean")]) == 1
    capsys.readouterr()
    assert main([p("clean"), p("npd")]) == 1


def test_a_batch_of_clean_files_passes(capsys):
    assert main([p("clean"), p("clean")]) == 0


def test_list_datablocks_success_exits_zero(capsys):
    assert main([p("multiblock"), "--list-datablocks"]) == 0


def test_unparseable_cif_exits_two(capsys, tmp_path):
    bad = tmp_path / "broken.cif"
    bad.write_text("data_x\n_this is not valid cif ;;; [\n")
    assert main([str(bad)]) == 2


def test_a_missing_file_does_not_abandon_the_rest(capsys):
    assert main(["/nonexistent/nope.cif", p("clean")]) == 2
    out = capsys.readouterr().out
    assert "datablock: clean" in out, "clean.cif was never audited"


def test_bad_input_beats_fail_on_never(capsys, tmp_path):
    """'never' is about findings; it must not hide that a file was unreadable."""
    assert main(["/nonexistent/nope.cif", "--fail-on", "never"]) == 2


def test_json_batch_reports_every_file_including_errors(capsys):
    import json
    main(["/nonexistent/nope.cif", p("clean"), "--json"])
    rows = json.loads(capsys.readouterr().out)
    assert len(rows) == 2
    assert "error" in rows[0] and rows[1]["datablock"] == "clean"
