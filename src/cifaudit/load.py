"""CIF loading and a normalised view of one datablock.

A CIF is UNTRUSTED INPUT. It is a text file that may have come from anywhere, and
its free-text fields (_publ_*, _exptl_special_details, embedded .res files) can
contain arbitrary prose, including text addressed to a language model. Therefore:

* input is size-capped before parsing;
* free-text fields are never read, never returned and never echoed;
* atom labels and space-group symbols ARE returned, but sanitised: control
  characters stripped, length capped, and newlines removed, so a crafted label
  cannot inject structure into a report;
* nothing in a CIF is ever treated as an instruction.

No network access. No LLM in this path.
"""

import math
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

try:
    import gemmi
    HAVE_GEMMI = True
except ImportError:  # pragma: no cover - exercised by the import-guard test
    gemmi = None
    HAVE_GEMMI = False

MAX_BYTES = 8 * 1024 * 1024     # 8 MiB; deposited small-molecule CIFs are ~10-200 KiB
MAX_LABEL = 32
MAX_SITES = 100000

_CTRL = re.compile(r"[\x00-\x1f\x7f]")
_SU = re.compile(r"^\s*([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*(?:\(([0-9]+)\))?\s*$")
_ELEM = re.compile(r"^([A-Za-z]{1,2})")


class CifError(Exception):
    """Input could not be used. Carries a reason safe to show a caller."""


def sanitise(text: Optional[str], limit: int = MAX_LABEL) -> Optional[str]:
    """Make a CIF-supplied string safe to put in a report."""
    if text is None:
        return None
    out = _CTRL.sub("", str(text)).strip().strip("'\"")
    if len(out) > limit:
        out = out[:limit] + "..."
    return out


def parse_number(raw: Optional[str]) -> Tuple[Optional[float], Optional[float]]:
    """Parse a CIF numeric value, returning (value, standard_uncertainty).

    Handles ``0.0234(12)``. Returns (None, None) for ``.``, ``?``, '' or junk.
    The s.u. is scaled to the last decimal place of the value, as CIF defines it.
    """
    if raw is None:
        return (None, None)
    s = str(raw).strip()
    if s in (".", "?", ""):
        return (None, None)
    m = _SU.match(s)
    if not m:
        return (None, None)
    val = float(m.group(1))
    su = None
    if m.group(2) is not None:
        digits = m.group(1)
        if "." in digits:
            dp = len(digits.split(".", 1)[1])
        else:
            dp = 0
        su = int(m.group(2)) * (10.0 ** -dp)
    return (val, su)


def element_of(type_symbol: Optional[str], label: Optional[str]) -> Optional[str]:
    """Element from _atom_site_type_symbol, falling back to the label."""
    for src in (type_symbol, label):
        if not src:
            continue
        m = _ELEM.match(str(src).strip().strip("'\""))
        if m:
            sym = m.group(1)
            return sym[0].upper() + sym[1:].lower()
    return None


@dataclass
class Site:
    label: str
    type_symbol: Optional[str]
    element: Optional[str]
    fract: Tuple[Optional[float], Optional[float], Optional[float]]
    occupancy: Optional[float]
    occupancy_su: Optional[float]
    uiso: Optional[float]
    adp_type: Optional[str]
    multiplicity: Optional[float]
    disorder_group: Optional[str]
    disorder_assembly: Optional[str]
    uij: Optional[Dict[str, float]] = None   # U11..U23 when anisotropic

    @property
    def is_hydrogen(self) -> bool:
        return (self.element or "").upper() in ("H", "D")


@dataclass
class Cell:
    a: Optional[float] = None
    b: Optional[float] = None
    c: Optional[float] = None
    alpha: Optional[float] = None
    beta: Optional[float] = None
    gamma: Optional[float] = None
    volume: Optional[float] = None
    su: Dict[str, Optional[float]] = field(default_factory=dict)

    def complete(self) -> bool:
        return all(v is not None for v in
                   (self.a, self.b, self.c, self.alpha, self.beta, self.gamma))

    def computed_volume(self) -> Optional[float]:
        if not self.complete():
            return None
        ca = math.cos(math.radians(self.alpha))
        cb = math.cos(math.radians(self.beta))
        cg = math.cos(math.radians(self.gamma))
        t = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
        if t <= 0:
            return None
        return self.a * self.b * self.c * math.sqrt(t)


@dataclass
class Block:
    """One CIF datablock, normalised. Free-text fields are deliberately absent."""
    name: str
    cell: Cell
    sites: List[Site]
    hm_symbol: Optional[str]
    hall_symbol: Optional[str]
    symops: List[str]
    formula_sum: Optional[str]
    formula_moiety: Optional[str]
    formula_weight: Optional[float]
    z: Optional[float]
    density_reported: Optional[float]
    h_treatment: Optional[str]
    radiation_probe: Optional[str]
    present: set = field(default_factory=set)   # CIF data names actually present
    _raw: Any = None

    def has(self, *names: str) -> bool:
        return all(n in self.present for n in names)


_CELL_TAGS = (("a", "_cell_length_a"), ("b", "_cell_length_b"), ("c", "_cell_length_c"),
              ("alpha", "_cell_angle_alpha"), ("beta", "_cell_angle_beta"),
              ("gamma", "_cell_angle_gamma"), ("volume", "_cell_volume"))

# gemmi: the first tag is required, "?"-prefixed tags are optional. Index order
# below is relied on by the row unpacking in load().
_SITE_LOOP = [
    "_atom_site_label", "?_atom_site_type_symbol",
    "?_atom_site_fract_x", "?_atom_site_fract_y", "?_atom_site_fract_z",
    "?_atom_site_occupancy", "?_atom_site_U_iso_or_equiv",
    "?_atom_site_adp_type", "?_atom_site_site_symmetry_multiplicity",
    "?_atom_site_symmetry_multiplicity",
    "?_atom_site_disorder_group", "?_atom_site_disorder_assembly",
]

_ANISO_LOOP = ["_atom_site_aniso_label",
               "?_atom_site_aniso_U_11", "?_atom_site_aniso_U_22",
               "?_atom_site_aniso_U_33", "?_atom_site_aniso_U_12",
               "?_atom_site_aniso_U_13", "?_atom_site_aniso_U_23"]


def _val(block, *names):
    for n in names:
        v = block.find_value(n)
        if v is not None:
            return gemmi.cif.as_string(v)
    return None


def datablock_names(text: str) -> List[str]:
    """Names of every datablock, so the wrong one is never audited silently."""
    doc = _read(text)
    return [sanitise(b.name, 64) or "?" for b in doc]


def _read(text: str):
    if not HAVE_GEMMI:
        raise CifError("gemmi is required to parse CIF files; pip install gemmi")
    if text is None:
        raise CifError("no CIF content supplied")
    data = text.encode("utf-8", "replace") if isinstance(text, str) else bytes(text)
    if len(data) > MAX_BYTES:
        raise CifError("CIF exceeds the %d MiB cap (%d bytes); refusing to parse"
                       % (MAX_BYTES // (1024 * 1024), len(data)))
    try:
        doc = gemmi.cif.read_string(data.decode("utf-8", "replace"))
    except Exception as exc:
        raise CifError("CIF could not be parsed: %s" % sanitise(str(exc), 200))
    if len(doc) == 0:
        raise CifError("CIF contains no datablocks")
    return doc


def load(text: str, datablock: Optional[str] = None) -> Block:
    """Parse one datablock into a normalised Block."""
    doc = _read(text)
    blk = None
    if datablock:
        for b in doc:
            if b.name == datablock:
                blk = b
                break
        if blk is None:
            raise CifError("datablock %r not found; present: %s"
                           % (sanitise(datablock, 64), datablock_names(text)))
    else:
        blk = doc[0]

    present = set()
    for item in blk:
        if item.pair is not None:
            present.add(item.pair[0])
        elif item.loop is not None:
            for t in item.loop.tags:
                present.add(t)

    cell = Cell(su={})
    for attr, tag in _CELL_TAGS:
        v, su = parse_number(_val(blk, tag))
        setattr(cell, attr, v)
        cell.su[attr] = su

    # anisotropic Us, keyed by label
    aniso = {}
    atab = blk.find(_ANISO_LOOP)
    for row in atab:
        lab = sanitise(gemmi.cif.as_string(row[0])) if row.has(0) else None
        uij = {}
        for i, key in enumerate(("U11", "U22", "U33", "U12", "U13", "U23"), start=1):
            v, _ = parse_number(row[i] if row.has(i) else None)
            if v is not None:
                uij[key] = v
        if lab:
            aniso[lab] = uij

    sites: List[Site] = []
    tab = blk.find(_SITE_LOOP)
    for row in tab:
        if len(sites) >= MAX_SITES:
            break
        label = sanitise(gemmi.cif.as_string(row[0])) or ""
        tsym = sanitise(gemmi.cif.as_string(row[1])) if row.has(1) else None
        fx, _ = parse_number(row[2] if row.has(2) else None)
        fy, _ = parse_number(row[3] if row.has(3) else None)
        fz, _ = parse_number(row[4] if row.has(4) else None)
        occ, occ_su = parse_number(row[5] if row.has(5) else None)
        uiso, _ = parse_number(row[6] if row.has(6) else None)
        adp = sanitise(gemmi.cif.as_string(row[7])) if row.has(7) else None
        mult, _ = parse_number(row[8] if row.has(8) else None)
        if mult is None:
            mult, _ = parse_number(row[9] if row.has(9) else None)
        dgroup = sanitise(gemmi.cif.as_string(row[10])) if row.has(10) else None
        dasm = sanitise(gemmi.cif.as_string(row[11])) if row.has(11) else None
        sites.append(Site(
            label=label, type_symbol=tsym, element=element_of(tsym, label),
            fract=(fx, fy, fz), occupancy=occ, occupancy_su=occ_su, uiso=uiso,
            adp_type=adp, multiplicity=mult, disorder_group=dgroup,
            disorder_assembly=dasm, uij=aniso.get(label)))

    symops = []
    for tag in ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"):
        t = blk.find([tag])
        for row in t:
            s = sanitise(gemmi.cif.as_string(row[0]), 64)
            if s:
                symops.append(s)
        if symops:
            break

    fw, _ = parse_number(_val(blk, "_chemical_formula_weight"))
    zz, _ = parse_number(_val(blk, "_cell_formula_units_Z"))
    dens, _ = parse_number(_val(blk, "_exptl_crystal_density_diffrn"))

    return Block(
        name=sanitise(blk.name, 64) or "?",
        cell=cell, sites=sites,
        hm_symbol=sanitise(_val(blk, "_symmetry_space_group_name_H-M",
                                "_space_group_name_H-M_alt"), 32),
        hall_symbol=sanitise(_val(blk, "_symmetry_space_group_name_Hall",
                                  "_space_group_name_Hall"), 32),
        symops=symops,
        formula_sum=sanitise(_val(blk, "_chemical_formula_sum"), 200),
        formula_moiety=sanitise(_val(blk, "_chemical_formula_moiety"), 200),
        formula_weight=fw, z=zz, density_reported=dens,
        h_treatment=sanitise(_val(blk, "_atom_sites_solution_hydrogens"), 32),
        radiation_probe=sanitise(_val(blk, "_diffrn_radiation_probe"), 24),
        present=present, _raw=blk)
