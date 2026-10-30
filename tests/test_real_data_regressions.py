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


def test_bare_dot_cif_numbers_are_parsed(cif):
    """CIF writes "1.", "0." and ".5"; older COD entries use them throughout.

    Requiring a digit after the decimal point made parse_number reject them,
    which over 20 out-of-sample COD structures produced 171 spurious
    "unparseable occupancy" findings and 36 spurious "incomplete U(i,j)"
    findings (a "0." U13 read as absent).
    """
    from cifaudit.load import load
    b = load(cif("bare_dot_numbers"))
    assert b.cell.a == 10.0 and b.cell.volume == 1000.0
    assert len(b.sites) == 2
    for s in b.sites:
        assert s.occupancy == 1.0, s.label
        assert all(v is not None for v in s.fract), s.label
        assert s.uij is not None and len(s.uij) == 6, s.label


def test_bare_dot_file_produces_no_spurious_findings(run):
    c = codes(run("bare_dot_numbers"))
    assert "OCC_UNPARSEABLE" not in c
    assert "ADP_UIJ_INCOMPLETE" not in c
    assert "CELL_INCOMPLETE" not in c


def test_dot_and_question_mark_are_still_missing_values():
    """"." and "?" are CIF's inapplicable/unknown markers, not numbers."""
    from cifaudit.load import parse_number
    assert parse_number(".") == (None, None)
    assert parse_number("?") == (None, None)


def test_special_position_copies_are_deduplicated(cif):
    """gemmi's unit-cell expansion can emit coincident copies of an atom on a
    special position, which then read as two atoms 0.000 A apart (seen on COD
    1000200, Al1 at 0,0,0 with multiplicity 4)."""
    from cifaudit.expand import expand
    from cifaudit.load import load
    b = load(cif("bare_dot_numbers"))
    sites, _cell = expand(b)
    keys = [(s.label, round(s.fract[0] % 1, 5), round(s.fract[1] % 1, 5),
             round(s.fract[2] % 1, 5)) for s in sites]
    assert len(keys) == len(set(keys)), "coincident duplicates in expansion"


def test_bonding_uses_an_additive_tolerance_not_a_factor():
    """The field's convention is d <= r_i + r_j + epsilon, epsilon ~0.35-0.45 A.

    cifaudit previously used a multiplicative 1.2x, which no standard tool uses.
    arXiv 2601.02017 derives bonding radii from 519,021 COD records and states
    epsilon is "usually ranging between 0.35 A and 0.45 A"; Jmol-style and
    MaterialsCoord use 0.45 A; CoRE MOF 2019 uses Cordero + 0.2/0.4 A.
    """
    from cifaudit.audit import DEFAULTS
    assert "bond_factor" not in DEFAULTS
    assert 0.35 <= DEFAULTS["bond_tolerance"] <= 0.45


def test_hydrogen_bonds_to_exactly_one_heavy_atom(cif):
    """Radius-based bonding does not work for H.

    The COD-derived radii set has no hydrogen value and discards every contact
    involving hydrogen. PLATON PLAT303: "Hydrogen atoms are generally connected
    to only one other atom." So an H is bonded to its nearest heavy atom within
    xh_bond_max, and never to two.
    """
    from cifaudit.audit import DEFAULTS
    from cifaudit.elements import radius
    from cifaudit.expand import expand
    from cifaudit.geom import bond_graph
    from cifaudit.load import load
    b = load(cif("clean"))
    sites, cell = expand(b)
    sites = [s for s in sites if s.element and all(v is not None for v in s.fract)]
    _d, adj, _c, _x = bond_graph(sites, cell, radius,
                                 DEFAULTS["bond_tolerance"],
                                 DEFAULTS["xh_bond_max"],
                                 DEFAULTS["shared_site_max"],
                                 DEFAULTS["occ_tol"])
    for i, s in enumerate(sites):
        if (s.element or "").upper() in ("H", "D"):
            heavy = [k for k in adj[i]
                     if (sites[k].element or "").upper() not in ("H", "D")]
            assert len(heavy) <= 1, "%s bonded to %d heavy atoms" % (s.label,
                                                                     len(heavy))


def test_h_multiple_connections_is_informational_not_a_warning():
    """PLAT303 says to investigate a disorder model, not that it is an error."""
    from cifaudit.model import kind_of
    assert kind_of("H_MULTIPLE_CONNECTIONS") == "structural"
    import pathlib
    src = pathlib.Path("src/cifaudit/checks/hydrogens.py").read_text()
    i = src.index('code="H_MULTIPLE_CONNECTIONS"')
    assert "severity=INFO" in src[i:i + 120]


def test_adp_anisotropy_uses_sqrt_of_the_eigenvalue_ratio():
    """PLAT213 tests SQRT(U3/U1) of the principal axes, not the raw ratio.

    U is a mean-square displacement in A^2, so the square root is the ratio of
    ellipsoid axis lengths. The old default of 10 on the raw ratio equals
    sqrt = 3.16, only the 97.5th percentile over 12,899 non-H anisotropic atoms
    in 500 deposited COD structures, and it fired on 24% of files.
    """
    import math
    from cifaudit.audit import DEFAULTS
    assert "adp_ratio_max" not in DEFAULTS
    assert DEFAULTS["adp_axis_ratio_max"] == 4.0
    # 4.0 in axis units is 16.0 in raw eigenvalue units
    assert math.isclose(DEFAULTS["adp_axis_ratio_max"] ** 2, 16.0)


def test_adp_anisotropy_finding_reports_axis_ratio_and_shape(cif):
    """The evidence must carry sqrt(U3/U1) and PLATON's oblate/prolate call."""
    from cifaudit import audit
    src = cif("clean")
    # make C1 extremely prolate: U11 large, U22/U33 small
    stretched = src.replace("C1 0.0250 0.0250 0.0250 0.0000 0.0000 0.0000",
                            "C1 0.4000 0.0050 0.0050 0.0000 0.0000 0.0000")
    hits = [f for f in audit(stretched).to_dict()["findings"]
            if f["code"] == "ADP_EXTREME_ANISOTROPY"]
    assert hits, "expected an anisotropy finding"
    ev = hits[0]["evidence"]
    assert "axis_ratio_sqrt_U3_U1" in ev and "raw_eigenvalue_ratio" in ev
    assert ev["shape"] in ("oblate", "prolate")
    assert "PLAT213" in [c["code"] for c in hits[0]["checkcif"]]


def test_ueq_range_is_tested_per_element_type(cif):
    """PLAT220 exists to catch misassigned species ("O versus N"), so the Ueq
    range must be compared within one element, not across all non-H atoms."""
    import pathlib
    src = pathlib.Path("src/cifaudit/checks/adp.py").read_text()
    assert "by_element" in src
    i = src.index('code="UEQ_RANGE_WIDE"')
    assert "element" in src[i:i + 600]


def test_two_elements_sharing_one_site_are_not_a_clash(run):
    """COD 1000241 has Sb1 and W1 both at (0,0,0) on Wyckoff 16c, occupancy 0.5
    each - substitutional sharing. Only one is present in any unit cell.

    Modelling them as co-existing atoms made every such site read as atoms
    0.000 A apart: 137 of 500 deposited structures, 27%.
    """
    c = codes(run("shared_site"))
    assert "ATOMS_TOO_CLOSE" not in c
    assert "SITE_ON_OWN_IMAGE" not in c
    assert "SHARED_SITE_OCC_GT_ONE" not in c


def test_shared_site_with_impossible_occupancy_is_an_error(run):
    """Alternatives are excused from the clash test, but their occupancies must
    still be able to coexist. 0.5 + 0.9 = 1.4 cannot."""
    r = run("shared_site_bad")
    hits = [f for f in r["findings"] if f["code"] == "SHARED_SITE_OCC_GT_ONE"]
    assert hits, codes(r)
    assert hits[0]["severity"] == "error"
    assert hits[0]["evidence"]["sum"] > 1.0


def test_shared_site_atoms_do_not_bond_to_each_other(cif):
    from cifaudit.audit import DEFAULTS
    from cifaudit.elements import radius
    from cifaudit.expand import expand
    from cifaudit.geom import bond_graph
    from cifaudit.load import load
    b = load(cif("shared_site"))
    sites, cell = expand(b)
    _d, adj, _c, excl = bond_graph(sites, cell, radius,
                                   DEFAULTS["bond_tolerance"],
                                   DEFAULTS["xh_bond_max"],
                                   DEFAULTS["shared_site_max"],
                                   DEFAULTS["occ_tol"])
    idx = {s.label: i for i, s in enumerate(sites)}
    i, j = idx["Sb1"], idx["W1"]
    assert (min(i, j), max(i, j)) in excl
    assert j not in adj[i] and i not in adj[j]


def test_dummy_atoms_are_excluded_from_geometry(cif):
    """_atom_site_calc_flag 'dum' marks a placeholder whose position is
    undetermined. COD 1000241 carries its H at (-1,-1,-1), which wraps onto the
    shared Sb/W site and read as an H bonded to two heavy atoms 0.000 A away."""
    from cifaudit.expand import expand
    from cifaudit.load import load
    b = load(cif("dummy_atom"))
    assert [s.label for s in b.sites if s.is_dummy] == ["H1"]
    geom, _c = expand(b)
    assert "H1" not in {s.label for s in geom}


def test_dummy_atoms_still_count_towards_the_formula(run):
    """A dummy atom's position is undetermined but it IS part of the
    composition. COD 1000241 declares 'H3 O7 Sb W' with Z=8 and carries its
    hydrogen as a dummy with occupancy 0.75 on a 32-fold site: 32 x 0.75 / 8 = 3.
    """
    assert "FORMULA_SITES_MISMATCH" not in codes(run("dummy_atom"))


def test_composition_prefers_the_declared_multiplicity(cif):
    """The CIF states _atom_site_symmetry_multiplicity; expansion can disagree.

    For COD 1000241's dummy H the declared multiplicity is 32 while gemmi
    expands the placeholder coordinates to 16 copies, so geometric expansion
    gave 1.5 H per formula unit against the declared 3.
    """
    from cifaudit import audit
    src = cif("dummy_atom")
    d = audit(src).to_dict()
    assert "FORMULA_SITES_MISMATCH" not in {f["code"] for f in d["findings"]}
    # strip the multiplicities and the basis must fall back to expansion
    stripped = src.replace("_atom_site_symmetry_multiplicity\n", "")
    stripped = stripped.replace("O1 O 1 0.20000", "O1 O 0.20000")
    stripped = stripped.replace("H1 H 2 -1.00000", "H1 H -1.00000")
    audit(stripped).to_dict()      # must not raise


def test_a_disordered_bond_is_not_a_shared_site(run, cif):
    """COD 4321002 has H4 at 0.930 A from C4, both occupancy 0.5.

    Their occupancies sum to exactly 1.0, so the first version of the
    shared-site rule classified the hydrogen as an ALTERNATIVE to the carbon it
    is bonded to. That orphaned the H and, across the 500-structure sample,
    added ISOLATED_H to 55 more files than the rule fixed.

    Atoms sharing one site are never a bond length apart, so an alternative must
    also be closer than `shared_site_bond_fraction` x (r_i + r_j).
    """
    c = codes(run("disordered_ch"))
    assert "ISOLATED_H" not in c

    from cifaudit.audit import DEFAULTS
    from cifaudit.elements import radius
    from cifaudit.expand import expand
    from cifaudit.geom import bond_graph
    from cifaudit.load import load
    b = load(cif("disordered_ch"))
    sites, cell = expand(b)
    _d, adj, _h, excl = bond_graph(sites, cell, radius,
                                   DEFAULTS["bond_tolerance"],
                                   DEFAULTS["xh_bond_max"],
                                   DEFAULTS["shared_site_max"],
                                   DEFAULTS["occ_tol"],
                                   DEFAULTS["shared_site_bond_fraction"])
    idx = {s.label: i for i, s in enumerate(sites)}
    i, j = idx["C4"], idx["H4"]
    assert (min(i, j), max(i, j)) not in excl, "bonded pair marked as alternatives"
    assert j in adj[i], "the disordered C-H bond was dropped"


def test_composition_uses_expansion_except_for_dummies(cif):
    """`_atom_site_symmetry_multiplicity` carries two incompatible meanings in
    the wild - Wyckoff multiplicity in IUCr-style inorganic CIFs, site-symmetry
    order in SHELX-style ones. Preferring it for every site made
    FORMULA_SITES_MISMATCH fire on 167 of 500 files instead of 94, so it is now
    used only for dummy sites, whose coordinates cannot be expanded."""
    from cifaudit import audit
    d = audit(cif("dummy_atom")).to_dict()
    assert "FORMULA_SITES_MISMATCH" not in {f["code"] for f in d["findings"]}
    hits = [f for f in audit(cif("shared_site_bad")).to_dict()["findings"]
            if f["code"] == "FORMULA_SITES_MISMATCH"]
    if hits:
        assert "expansion" in hits[0]["evidence"]["basis"]
