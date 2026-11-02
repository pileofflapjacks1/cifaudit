"""Defects found by the third precision round, each confirmed by independent
recomputation before it was fixed. Fixtures are synthetic and built inline; no
COD text is included.
"""


def _cif(rows, sg="P 1", ops=("x, y, z",), formula="C1", groups=False):
    head = """data_t
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
""" % (formula, sg, "\n".join(" '%s'" % o for o in ops))
    if groups:
        head += "_atom_site_disorder_assembly\n_atom_site_disorder_group\n"
    return head + "\n".join(rows) + "\n"


def _run(text):
    from cifaudit import audit
    return audit(text).to_dict()


def _codes(text):
    return {f["code"] for f in _run(text)["findings"]}


# ---- X-H reach -------------------------------------------------------------

def test_a_long_p_h_bond_is_a_bond():
    """COD 2014097: P-H at 1.351 A was read as an unbonded H."""
    text = _cif(["P1 P 0.50000 0.50000 0.50000 1.000 0.0300",
                 "H1 H 0.63700 0.50000 0.50000 1.000 0.0400"], formula="H P")
    assert "ISOLATED_H" not in _codes(text)


def test_a_bridging_metal_hydride_is_a_bond():
    """COD 4123099: Lu-H hydrides at 1.97-2.14 A were read as unbonded H."""
    text = _cif(["Lu1 Lu 0.50000 0.50000 0.50000 1.000 0.0200",
                 "H1 H 0.70000 0.50000 0.50000 1.000 0.0400"], formula="H Lu")
    assert "ISOLATED_H" not in _codes(text)


def test_a_c_h_near_a_metal_is_not_a_multiple_connection():
    """The longer metal reach decides bonding for an H with no light neighbour;
    it must not make every C-H near a metal a 'multiple connection'."""
    text = _cif(["C1 C 0.50000 0.50000 0.50000 1.000 0.0300",
                 "H1 H 0.59600 0.50000 0.50000 1.000 0.0400",
                 "Ni1 Ni 0.59600 0.69000 0.50000 1.000 0.0200"],
                formula="C H Ni")
    assert "H_MULTIPLE_CONNECTIONS" not in _codes(text)


def test_xh_reach_values():
    from cifaudit.elements import radius
    from cifaudit.geom import xh_reach
    assert xh_reach(radius("C"), 1.35) == 1.35
    assert xh_reach(radius("O"), 1.35) == 1.35
    assert abs(xh_reach(radius("P"), 1.35) - 1.53) < 1e-9
    assert abs(xh_reach(radius("Lu"), 1.35) - 2.33) < 1e-9


# ---- zero occupancy --------------------------------------------------------

def test_a_tiny_refined_occupancy_is_not_zero():
    """COD 9007098: FeM 0.00100 was reported as 'occupancy 0.0'."""
    text = _cif(["Fe1 Fe 0.50000 0.50000 0.50000 0.001 0.0200"], formula="Fe0.001")
    assert "OCC_ZERO" not in _codes(text)


def test_an_exact_zero_occupancy_is_still_reported():
    text = _cif(["Fe1 Fe 0.50000 0.50000 0.50000 0.000 0.0200"], formula="Fe")
    assert "OCC_ZERO" in _codes(text)


# ---- disorder closure ------------------------------------------------------

def test_a_negative_part_group_is_not_summed():
    """COD 7110006: a methyl at PART -1 on a mirror, occupancy 0.5, was reported
    as summing to 0.5; its other half is its own image."""
    text = _cif(["C1 C 0.50000 0.50000 0.50000 1.000 0.0300 . .",
                 "H1 H 0.50000 0.60000 0.55000 0.500 0.0400 A -1",
                 "H2 H 0.50000 0.40000 0.55000 0.500 0.0400 A -1"],
                formula="C H", groups=True)
    r = _run(text)
    assert "DISORDER_SUM_NOT_ONE" not in {f["code"] for f in r["findings"]}
    assert any("occupancy closure" in n["area"] for n in r["not_tested"])


def test_full_weight_atoms_tagged_into_a_group_are_bookkeeping():
    """COD 4123051: two occupancy-1 H tagged group 1 gave a sum of 1.25 for a
    real 0.75 + 0.25 disorder."""
    text = _cif(["C1 C 0.20000 0.20000 0.20000 0.750 0.0300 B 1",
                 "H1 H 0.30000 0.20000 0.20000 1.000 0.0400 B 1",
                 "C2 C 0.60000 0.60000 0.60000 0.250 0.0300 B 2"],
                formula="C H", groups=True)
    assert "DISORDER_SUM_NOT_ONE" not in _codes(text)


def test_a_genuinely_wrong_sum_is_still_reported():
    text = _cif(["C1 C 0.20000 0.20000 0.20000 0.750 0.0300 B 1",
                 "C2 C 0.60000 0.60000 0.60000 0.750 0.0300 B 2"],
                formula="C", groups=True)
    assert "DISORDER_SUM_NOT_ONE" in _codes(text)


# ---- copy merging ----------------------------------------------------------

_P1BAR = dict(sg="P -1", ops=("x, y, z", "-x, -y, -z"))


def _n_sites(text):
    from cifaudit.expand import expand
    from cifaudit.load import load
    return len(expand(load(text))[0])


def test_a_negative_part_atom_just_off_its_element_is_one_site():
    """COD 7013013: PART -1 atoms 0.039 and 0.068 A from their mirror images
    must be treated alike; the formula matches only when they are merged."""
    text = _cif(["C1 C 0.00350 0.00000 0.00000 0.500 0.0300 . -1"],
                formula="C0.5", groups=True, **_P1BAR)
    assert _n_sites(text) == 1


def test_an_ordinary_atom_at_the_same_offset_keeps_both_copies():
    """Outside negative PART the 0.05 A limit stands: genuine copies 0.05-0.15 A
    apart exist (COD 9005049, 7203095)."""
    text = _cif(["C1 C 0.00350 0.00000 0.00000 0.500 0.0300 . ."],
                formula="C", groups=True, **_P1BAR)
    assert _n_sites(text) == 2


# ---- assembly labels -------------------------------------------------------

def test_halves_of_one_disorder_under_different_letters_are_alternatives():
    """COD 7110032: H25A (A, group 1) sat 0.57 A from C25' (no letter, group 2),
    the alternative position of its own C25."""
    text = _cif(["C1 C 0.50000 0.50000 0.50000 0.850 0.0300 A 1",
                 "H1 H 0.59800 0.50000 0.50000 0.850 0.0400 A 1",
                 "C2 C 0.59800 0.55700 0.50000 0.150 0.0300 . 2"],
                formula="C H0.85", groups=True)
    assert "ATOMS_TOO_CLOSE" not in _codes(text)


def test_oversubscribed_groups_under_different_letters_still_clash():
    text = _cif(["C1 C 0.50000 0.50000 0.50000 0.850 0.0300 A 1",
                 "H1 H 0.59800 0.50000 0.50000 0.850 0.0400 A 1",
                 "C2 C 0.59800 0.55700 0.50000 0.850 0.0300 . 2"],
                formula="C2 H", groups=True)
    assert "ATOMS_TOO_CLOSE" in _codes(text)
