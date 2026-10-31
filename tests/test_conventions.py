"""Crystallographic conventions, each adopted after the precision study.

Where a convention came from the expert adjudicator (J. R. Ramirez) the test
says so. Each fixture is synthetic; no COD text is included.
"""

from conftest import codes


def _find(run, fixture, code):
    return [f for f in run(fixture)["findings"] if f["code"] == code]


# ---- isolated atoms --------------------------------------------------------

def test_isolated_o_is_information_described_as_solvent(run):
    """Expert convention: an unbound single O is routinely modelled solvent
    (water with unlocated H, or O placed for electron density)."""
    hits = _find(run, "isolated_o", "ISOLATED_ATOM")
    assert hits and hits[0]["severity"] == "info"
    assert "solvent" in hits[0]["message"]


def test_a_secondary_contact_is_not_isolation(run):
    """Tl-N 2.74 A is outside the covalent tolerance but plainly coordination;
    the precision study had Tl-N 2.74, Pb-S 2.96 and axial Cu-Cl 2.90 A reported
    as isolated atoms. Expert choice: 'isolated' means nothing within
    r_i + r_j + 0.90 A."""
    assert "ISOLATED_ATOM" not in codes(run("secondary_contact"))


def test_contact_tolerance_is_the_agreed_value():
    from cifaudit.audit import DEFAULTS
    assert DEFAULTS["contact_tolerance"] == 0.90
    assert DEFAULTS["bond_tolerance"] == 0.45


# ---- composition -----------------------------------------------------------

def test_attached_hydrogens_count_towards_the_formula(run, cif):
    from cifaudit.load import load
    assert load(cif("attached_h")).sites[0].attached_h == 3
    c = codes(run("attached_h"))
    assert "FORMULA_SITES_MISMATCH" not in c
    assert "NO_HYDROGENS_WITH_CARBON" not in c


def test_squeeze_turns_a_formula_mismatch_into_information(run):
    """Expert convention: authors record SQUEEZE, and the declared formula then
    includes the solvent removed from the model."""
    sq = _find(run, "squeezed", "FORMULA_SITES_MISMATCH")
    assert sq and sq[0]["severity"] == "info"
    assert "SQUEEZE" in sq[0]["message"]
    assert sq[0]["evidence"]["solvent_removed"] == "SQUEEZE"


def test_the_same_mismatch_without_squeeze_is_still_a_warning(run):
    plain = _find(run, "unsqueezed", "FORMULA_SITES_MISMATCH")
    assert plain and plain[0]["severity"] == "warn"


def test_squeeze_is_detected_from_data_names_not_free_text(cif):
    from cifaudit.load import load
    assert load(cif("squeezed")).solvent_removed == "SQUEEZE"
    assert load(cif("unsqueezed")).solvent_removed is None


def test_a_declared_solid_solution_is_not_flagged(run):
    """0/10 sampled NON_INTEGRAL_CELL_CONTENT findings were right: every one was
    a solid solution whose CIF declares the non-integral formula."""
    c = codes(run("solid_solution"))
    assert "NON_INTEGRAL_CELL_CONTENT" not in c
    assert "FORMULA_SITES_MISMATCH" not in c


def test_non_integral_content_without_a_formula_is_still_reported(cif):
    from cifaudit import audit
    src = cif("solid_solution").replace(
        "_chemical_formula_sum            'Mg0.6 Fe0.4 O1'\n", "")
    c = {f["code"] for f in audit(src).to_dict()["findings"]}
    assert "NON_INTEGRAL_CELL_CONTENT" in c


# ---- old CIFs and deuterium ------------------------------------------------

def test_no_occupancy_column_softens_a_close_contact(run):
    hits = _find(run, "no_occupancy", "ATOMS_TOO_CLOSE")
    assert hits and hits[0]["severity"] == "warn"
    assert "declares no occupancies" in hits[0]["message"]
    assert hits[0]["evidence"]["occupancy_declared"] is False


def test_declared_occupancies_keep_a_close_contact_an_error(run):
    hits = _find(run, "short_contact", "ATOMS_TOO_CLOSE")
    assert hits and hits[0]["severity"] == "error"


def test_deuterium_without_a_probe_is_not_judged(run):
    r = run("deuterated")
    assert "SHORT_XH" not in codes(r)
    assert any("X-D distances" in n["area"] for n in r["not_tested"])


def test_deuterium_with_a_neutron_probe_says_c_d(cif):
    from cifaudit import audit
    src = cif("deuterated").replace(
        "_cell_formula_units_Z            1\n",
        "_cell_formula_units_Z            1\n_diffrn_radiation_probe neutron\n")
    hits = [f for f in audit(src).to_dict()["findings"] if f["code"] == "SHORT_XH"]
    assert hits, "0.80 A C-D is short against the 1.08 A neutron target"
    assert "for C-D" in hits[0]["message"] and "neutron" in hits[0]["message"]
