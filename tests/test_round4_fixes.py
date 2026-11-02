"""Defects found by the fourth precision round, each confirmed by independent
recomputation before it was fixed. Fixtures are synthetic; no COD text is
included.
"""


def _cif(rows, formula, sg="P 1", ops=("x, y, z",)):
    return """data_t
_chemical_formula_sum            '%s'
_cell_length_a                   10.0000
_cell_length_b                   10.0000
_cell_length_c                   10.0000
_cell_angle_alpha                90.00
_cell_angle_beta                 90.00
_cell_angle_gamma                90.00
_cell_formula_units_Z            1
_symmetry_space_group_name_H-M   '%s'
loop_
_symmetry_equiv_pos_as_xyz
%s
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
%s
""" % (formula, sg, "\n".join(" '%s'" % o for o in ops), "\n".join(rows))


def _found(text, code):
    from cifaudit import audit
    return [f for f in audit(text).to_dict()["findings"] if f["code"] == code]


# ---- zero-occupancy placeholders ------------------------------------------

def test_a_zero_occupancy_placeholder_is_not_an_atom():
    """COD 7035035: a Cp-ring centroid CT1, typed H with occupancy 0 and no
    'dum' flag, was read as an H bonded to five ring carbons."""
    rows = ["C1 C 0.50000 0.50000 0.50000 1.000 0.0300",
            "C2 C 0.62000 0.50000 0.50000 1.000 0.0300",
            "CT1 H 0.56000 0.50000 0.50000 0.000 0.0000"]
    text = _cif(rows, "C2")
    assert not _found(text, "H_MULTIPLE_CONNECTIONS")
    assert _found(text, "OCC_ZERO")
    from cifaudit.expand import expand
    from cifaudit.load import load
    assert "CT1" not in {s.label for s in expand(load(text))[0]}


# ---- SQUEEZE ---------------------------------------------------------------

def test_squeeze_does_not_excuse_a_model_larger_than_its_formula(cif):
    """COD 4342091: H 26 in the model, 24 in the formula, reported as
    'probably solvent removed'. Removed solvent only makes the formula larger."""
    text = cif("squeezed").replace("'C1 O2'", "'C1'")
    hits = [f for f in _found(text, "FORMULA_SITES_MISMATCH")]
    assert hits and hits[0]["severity"] == "warn"
    assert "cannot explain" in hits[0]["message"]


def test_squeeze_still_excuses_a_formula_larger_than_its_model(cif):
    hits = _found(cif("squeezed"), "FORMULA_SITES_MISMATCH")
    assert hits and hits[0]["severity"] == "info"


def test_a_fab_file_counts_as_solvent_removal(cif):
    """COD 7035070 records SQUEEZE only inside the embedded _shelx_fab_file
    text. The data name is the evidence; its content is never read."""
    text = cif("unsqueezed") + "\n_shelx_fab_file\n;\n   1    0    0    1.0  0.0\n;\n"
    hits = _found(text, "FORMULA_SITES_MISMATCH")
    assert hits and hits[0]["severity"] == "info"
    assert ".fab" in hits[0]["message"]
    from cifaudit.load import load
    assert load("data_x\n_cell_length_a 5\n").solvent_removed is None


# ---- two-bond symmetry disorder --------------------------------------------

_METHYL = ["C1 C 0.10500 0.00000 0.00000 %s 0.0400",
           "H1 H 0.08834 0.06388 0.06970 %s 0.0500",
           "H2 H 0.04207 0.02555 -0.06785 %s 0.0500"]


def _methyl(occ):
    return _cif([r % occ for r in _METHYL], "C H2", sg="P -1",
                ops=("x, y, z", "-x, -y, -z"))


def test_an_h_two_bonds_from_its_partners_image_is_an_alternative():
    """COD 4342081: in a methanol disordered across an inversion centre, H15B
    clashed with H15C', which is bonded to C15' - the image of C15 (0.5, 2.12 A
    from C15). H1 here is 1.58 A from H2'; neither H is close to its own image
    in a way the one-bond rule could use."""
    text = _methyl("0.500")
    assert not _found(text, "SHORT_HH_CONTACT")
    assert not _found(text, "ATOMS_TOO_CLOSE")


def test_the_same_geometry_at_full_occupancy_is_still_a_contact():
    assert _found(_methyl("1.000"), "SHORT_HH_CONTACT") or \
        _found(_methyl("1.000"), "ATOMS_TOO_CLOSE")
