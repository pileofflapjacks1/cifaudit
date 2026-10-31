"""checkCIF/PLATON test codes used by cifaudit.

Descriptions are quoted verbatim from the IUCr list of PLATON validation tests
(https://journals.iucr.org/services/cif/checking/platon_tests.html, read
2026-10-30), together with PLATON's own ``Type_N`` classification.

That page does NOT publish a mapping from Type_N to the checkCIF A/B/C/G alert
levels, so cifaudit does not claim one. ``type`` below is PLATON's Type_N as
printed; nothing here should be read as an alert level.
"""

from typing import Any, Dict

# code -> (PLATON Type_N, verbatim short description)
PLAT: Dict[str, Any] = {
    "PLAT003": (2, "Report the number of non-H atoms with Uiso or U(i,j) restraints"),
    "PLAT004": (5, "Report the dimensionality of a polymeric structure"),
    "PLAT007": (5, "Report on unrefined Donor-H-atoms"),
    "PLAT017": (1, "Check the consistency of scattering Type"),
    "PLAT040": (1, "Test for H-atoms [0,1]"),
    "PLAT041": (1, "Test SumFormula strings"),
    "PLAT042": (1, "Test MoietyFormula strings"),
    "PLAT043": (1, "Test for molecular weight values"),
    "PLAT044": (1, "Check reported against calculated density values"),
    "PLAT045": (1, "Check reported and calculated Z values"),
    "PLAT046": (1, "Check reported density with calculated density from Z * MW"),
    "PLAT047": (1, "Test sum formula given"),
    "PLAT048": (1, "Test moiety formula given"),
    "PLAT049": (1, "Check calculated density > 1.0"),
    "PLAT070": (1, "Test for duplicate atom labels"),
    "PLAT074": (1, "Test for occupancy = 0.0 (dummy atoms)"),
    "PLAT075": (1, "Test for occupancy > 1.0"),
    "PLAT076": (1, "Test for occupancy < 1.0 for atom on special position"),
    "PLAT077": (4, "Test for non-integral # of atoms in unit-cell"),
    "PLAT078": (1, "Test consistency 'geom'<->'no H-atoms' (_atom_sites_solution_hydrogens)"),
    "PLAT079": (1, "Test for inconsistency 'mixed' versus 'no H-atoms'"),
    "PLAT104": (1, "Test for additional translational symmetry [0, 1]"),
    "PLAT110": (2, "Test for additional translational symmetry [0, 1] (ADDSYM)"),
    "PLAT111": (2, "Test for additional centre of symmetry [0, 100] (ADDSYM)"),
    "PLAT112": (2, "Test for additional symmetry [0, 1] (ADDSYM)"),
    "PLAT113": (2, "Report New space-group suggested by ADDSYM"),
    "PLAT119": (1, "Report problem with symmetry operation syntax"),
    "PLAT120": (1, "Test for consistent _symmetry_space_group_name_H-M and Symm Opp"),
    "PLAT121": (1, "Test for valid _symmetry_space_group_name_H-M (+ aliases)"),
    "PLAT122": (1, "Test for ? _symmetry_space_group_name_H-M"),
    "PLAT123": (1, "Test for interpretable space-group symmetry"),
    "PLAT124": (1, "Test for _symmetry_equiv_pos_as_xyz present"),
    "PLAT128": (4, "Test for non-standard space_group settings"),
    "PLAT150": (1, "Check reported unit-cell Volume value against the calculated value"),
    "PLAT151": (1, "Check the reported s.u. on the unit-cell volume"),
    "PLAT164": (4, "Check for refined C-H H-atoms"),
    "PLAT180": (4, "Check rounding of cell axes and angles"),
    "PLAT211": (2, "Test for NPD ADP's in main residue(s)"),
    "PLAT212": (2, "Test for NPD ADP's in anion? & solvent ?"),
    "PLAT213": (2, "Test ratio adp max/min in main residue(s)"),
    "PLAT215": (3, "Test for unusual disordered atom ADP in main residue"),
    "PLAT217": (1, "Test for incomplete U(i,j) data"),
    "PLAT220": (2, "Test Ueq(max)/Ueq(Min) range for non-H atoms in non-solvent"),
    "PLAT221": (2, "Test Ueq(max)/Ueq(Min) range for non-H atoms in solvent"),
    "PLAT222": (3, "Test Uiso(max)/Uiso(Min) range for H atoms in non-solvent"),
    "PLAT230": (2, "Hirshfeld rigid-bond test [Acta Cryst. (1976), A32, 239-244]"),
    "PLAT250": (2, "Test for unusual anisotropic average U(i,j)"),
    "PLAT260": (2, "Test for high average Ueq for residue"),
    "PLAT299": (4, "Test for fixed partial occupancy sites"),
    "PLAT300": (4, "Test for fixed partial occupancy"),
    "PLAT303": (2, "Test for more than 1 connection to H-atoms"),
    "PLAT305": (2, "Test for isolated H-atoms"),
    "PLAT310": (2, "Test for 'too close' (symmetry related) full weight atoms"),
    "PLAT350": (3, "Test for short C - H (Angstrom difference) XRAY: 0.96 NEUT 1.08"),
    "PLAT352": (3, "Test for short N - H (Angstrom difference) XRAY: 0.87 NEUT 1.009"),
    "PLAT354": (3, "Test for short O - H (Angstrom difference) XRAY: 0.82 NEUT 0.983"),
    "PLAT360": (2, "Test for short C4 - C4 (Angstrom difference) XRAY: 1.54"),
    "PLAT410": (2, "Test for short non-bonding intra H..H contacts"),
    "PLAT411": (2, "Test for short non-bonding inter H..H contacts"),
    "PLAT420": (2, "Test for D-H bonds without acceptor"),
    "PLAT430": (2, "Test for short non-bonding inter D...A contacts"),
    "PLAT432": (2, "Test for short non-bonding inter X...Y contacts"),
    "PLAT482": (4, "Test for too small D-H...A angle"),
    "PLAT601": (None, "Structure Contains Solvent Accessible VOIDS"),
    "PLAT724": (2, "Test for consistency of contact distances and coordinates in the CIF"),
    "PLAT790": (4, "Test whether centre-of-gravity of a residue is within the unit-cell box"),
    "PLAT801": (4, "Test for missing, incomplete or out-of-order Cell data"),
    "PLAT806": (4, "Test for insufficient U(i,j) data in CIF loop"),
    "PLAT977": (2, "Test for negative density on H-atom positions"),
}

SOURCE = ("IUCr list of PLATON validation tests, "
          "https://journals.iucr.org/services/cif/checking/platon_tests.html "
          "(read 2026-10-30)")


def describe(code: str) -> Dict[str, Any]:
    """Return {code, platon_type, description} for a PLAT code.

    Unknown codes are returned with description None rather than invented.
    """
    entry = PLAT.get(code)
    if entry is None:
        return {"code": code, "platon_type": None, "description": None,
                "note": "not in cifaudit's PLAT table; look it up in the IUCr list"}
    ptype, desc = entry
    return {"code": code,
            "platon_type": ("Type_%d" % ptype) if ptype is not None else None,
            "description": desc}
