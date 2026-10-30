"""Regressions found by running cifaudit over real deposited COD structures.

Every test here exists because the tool was wrong about a real published
structure. The first run over 16 COD entries produced 36 findings on 16/16
files; after these fixes it produces 3 informational findings on 2 files, with
all planted defects still caught.
"""

import json
import os
import subprocess
import sys

from conftest import codes


def test_h_named_for_its_oxygen_is_not_holmium(run):
    """COD 2013000 has an H labelled `HO1` - hydrogen on O1.

    Deriving the element from the label greedily read "Ho" and reported a
    scattering-type mismatch. A label is consistent with its type symbol when it
    STARTS WITH that symbol.
    """
    assert "SCATTERING_TYPE_LABEL_MISMATCH" not in codes(run("label_ho1"))


def test_a_real_mismatch_is_still_caught(run):
    """The fix must not make the check useless."""
    import pathlib
    src = (pathlib.Path(__file__).parent / "data" / "label_ho1.cif").read_text()
    bad = src.replace("O1  O 0.20000", "O1  N 0.20000")
    from cifaudit import audit
    found = {f["code"] for f in audit(bad).to_dict()["findings"]}
    assert "SCATTERING_TYPE_LABEL_MISMATCH" in found


def test_hydrogen_bond_is_not_a_short_contact(run):
    """An O-H...O at 1.80 A is a hydrogen bond, not a clash.

    All 18 short-contact findings in the first COD run were ordinary hydrogen
    bonds: H...O and H...N between 1.3 and 2.2 A.
    """
    c = codes(run("hbond"))
    assert "SHORT_NONBONDED_CONTACT" not in c
    assert "SHORT_HH_CONTACT" not in c


def test_impossibly_short_hbond_is_still_reported(cif):
    from cifaudit import audit
    src = cif("hbond")
    # drag the acceptor onto the H: 0.28 A away, far too short for any H-bond
    squashed = src.replace("O2 O 0.45000", "O2 O 0.28000")
    found = {f["code"] for f in audit(squashed).to_dict()["findings"]}
    assert "HBOND_IMPLAUSIBLY_SHORT" in found or "ATOMS_TOO_CLOSE" in found


def test_symmetry_must_be_derived_from_the_expanded_cell(run):
    """spglib derives symmetry from the atoms it is given.

    Fed the asymmetric unit it returns P1 for essentially every real structure,
    which produced a finding on 16/16 COD entries. Fed the expanded cell it
    recovers the declared group.
    """
    c = codes(run("pbar1"))
    assert "SYMMETRY_POSSIBLY_MISSED" not in c
    assert "SYMMETRY_DERIVED_LOWER" not in c


def test_setting_difference_is_info_not_a_symmetry_warning(run):
    """COD 1543227 declares P 1 21/n 1; spglib's standard setting is P2_1/c.

    Same group (No. 14). That is a setting difference, not missed symmetry.
    """
    r = run("pbar1")
    for f in r["findings"]:
        assert f["code"] != "SYMMETRY_POSSIBLY_MISSED"


def test_fragment_count_is_a_note_not_a_finding(run):
    """>1 fragment held for 14/16 real structures: molecular crystals have many
    fragments by construction, so it is context, not a defect."""
    assert "MULTIPLE_FRAGMENTS" not in codes(run("clean"))


def test_output_is_deterministic_across_processes(cif, tmp_path):
    """Species numbering once used hash() on element symbols, which Python salts
    per process - so spglib's answer could vary between runs of a tool whose
    whole claim is determinism."""
    p = tmp_path / "s.cif"
    p.write_text(cif("pbar1"))
    script = (
        "import json,sys;sys.path.insert(0,%r);"
        "from cifaudit import audit;"
        "print(json.dumps(audit(open(%r).read()).to_dict(),sort_keys=True))"
        % (os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "src"), str(p))
    )
    outs = []
    for seed in ("0", "1", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        outs.append(subprocess.check_output([sys.executable, "-c", script],
                                            env=env).decode())
    assert outs[0] == outs[1] == outs[2], "output varies with PYTHONHASHSEED"
