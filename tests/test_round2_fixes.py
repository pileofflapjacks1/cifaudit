"""Defects found by the second precision round, each confirmed by independent
recomputation before it was fixed. Fixtures are synthetic; no COD text is
included.
"""

from conftest import codes


def _find(result, code):
    return [f for f in result["findings"] if f["code"] == code]


# ---- covalent radii --------------------------------------------------------

def test_an_atom_bonded_to_tantalum_is_not_isolated(run):
    """COD 4310078: a mu4-S 2.461 A from four Ta was reported isolated because
    Ta had no tabulated radius."""
    assert "ISOLATED_ATOM" not in codes(run("ta_neighbour"))


def test_radius_table_covers_cordero():
    from cifaudit.elements import radius
    for el, r in (("Ta", 1.70), ("Hf", 1.75), ("Gd", 1.96), ("Lu", 1.87),
                  ("Tc", 1.47), ("Th", 2.06), ("Cm", 1.69)):
        assert radius(el) == r


def test_isolation_abstains_next_to_an_element_with_no_radius(run):
    r = run("unknown_radius")
    assert "ISOLATED_ATOM" not in codes(r)
    assert any("isolated-atom test" in n["area"] and "S1" in n["area"]
               for n in r["not_tested"])


# ---- symmetry expansion ----------------------------------------------------

def test_a_site_disordered_across_a_centre_keeps_both_copies(cif):
    """COD 7041079: gemmi merged half-occupied H lying 0.29-0.38 A from their
    own mirror images, and the H count came out 25 instead of 28."""
    from cifaudit.expand import expand
    from cifaudit.load import load
    sites, _cell = expand(load(cif("h_across_centre")))
    assert len(sites) == 2
    assert all(s.occ == 0.5 for s in sites)


def test_disordered_copies_count_fully_towards_the_formula(run):
    assert "FORMULA_SITES_MISMATCH" not in codes(run("h_across_centre"))


def test_a_special_position_still_gives_one_site(cif):
    from cifaudit.expand import expand
    from cifaudit.load import load
    on_centre = cif("h_across_centre").replace("0.01500 0.00000", "0.00000 0.00000")
    sites, _cell = expand(load(on_centre))
    assert len(sites) == 1


# ---- disorder across a symmetry element ------------------------------------

def test_half_occupied_h_and_its_own_image_are_alternatives(run):
    """COD 7041061: an H 0.385 A from its own 2-fold image was an 'impossible'
    clash. Two half-occupied copies across a centre are never both present."""
    assert "ATOMS_TOO_CLOSE" not in codes(run("h_across_centre"))


def test_full_occupancy_on_its_own_image_is_still_an_error(cif):
    from cifaudit import audit
    full = cif("h_across_centre").replace("0.500 0.0500", "1.000 0.0500")
    found = {f["code"] for f in audit(full).to_dict()["findings"]}
    assert "ATOMS_TOO_CLOSE" in found


def _c1_pair_excluded(text):
    from cifaudit.audit import DEFAULTS as c
    from cifaudit.elements import radius
    from cifaudit.expand import expand
    from cifaudit.geom import bond_graph
    from cifaudit.load import load
    sites, cell = expand(load(text))
    _d, adj, _h, excl = bond_graph(
        sites, cell, radius, c["bond_tolerance"], c["xh_bond_max"],
        c["shared_site_max"], c["occ_tol"], c["shared_site_bond_fraction"])
    assert len(sites) == 2
    return (0, 1) in excl, 1 in adj[0]


def test_negative_part_is_not_bonded_to_its_own_image(cif):
    """SHELX PART -1 means 'do not connect to symmetry equivalents'. C1 is
    1.40 A from its inversion image - a bond length - but declared disordered
    across the centre."""
    excluded, bonded = _c1_pair_excluded(cif("negative_part"))
    assert excluded and not bonded


def test_without_negative_part_a_bond_length_image_stays_bonded(cif):
    """A half-occupied ORDERED molecule can be bonded to its own image across
    a centre; with no negative PART that bond is left alone."""
    text = cif("negative_part").replace("0.0400 . -1", "0.0400 . .")
    excluded, bonded = _c1_pair_excluded(text)
    assert bonded and not excluded


# ---- smaller fixes ---------------------------------------------------------

def test_a_short_bonded_xh_is_not_an_impossible_clash(run):
    """COD 4310039: a freely refined C-H of 0.719 A was reported as two atoms
    that cannot both be real. It is a short bond."""
    c = codes(run("short_bonded_xh"))
    assert "ATOMS_TOO_CLOSE" not in c
    assert "SHORT_XH" in c


def test_a_water_label_is_not_read_as_an_element(run):
    """'Wat1' with type O was reported as a label/type mismatch ('Wa')."""
    assert "SCATTERING_TYPE_LABEL_MISMATCH" not in codes(run("water_label"))


def test_label_element():
    from cifaudit.checks.labels import label_element
    assert label_element("Cl1") == "Cl"
    assert label_element("C1A") == "C"
    assert label_element("Fe1'") == "Fe"
    assert label_element("Wat1") is None
    assert label_element("OW1") is None
    assert label_element("N") == "N"


def test_a_real_label_type_mismatch_is_still_reported(cif):
    from cifaudit import audit
    bad = cif("water_label").replace("Wat1 O", "Fe1 O")
    found = {f["code"] for f in audit(bad).to_dict()["findings"]}
    assert "SCATTERING_TYPE_LABEL_MISMATCH" in found


def test_negative_part_counterparts_across_the_element_are_alternatives(run):
    """COD 7041061: H16A was 0.385 A from the image of H18A, its counterpart in
    the other orientation of a ring disordered across a 2-fold axis, while both
    were far from their OWN images. The group (C1 near the centre) overlaps its
    image, so every pair across the two copies is an alternative."""
    assert "ATOMS_TOO_CLOSE" not in codes(run("negative_part_group"))


def test_the_same_geometry_without_negative_part_is_still_a_clash(cif):
    from cifaudit import audit
    text = cif("negative_part_group").replace(". -1", ". .")
    found = {f["code"] for f in audit(text).to_dict()["findings"]}
    assert "ATOMS_TOO_CLOSE" in found
