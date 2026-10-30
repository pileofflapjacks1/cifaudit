"""Number, element and label parsing."""

import pytest

from cifaudit.load import CifError, element_of, parse_number, sanitise


@pytest.mark.parametrize("raw,val,su", [
    ("0.0234(12)", 0.0234, 0.0012),
    ("1.54", 1.54, None),
    ("0.5(2)", 0.5, 0.2),
    ("-0.003(11)", -0.003, 0.011),
    ("12(3)", 12.0, 3.0),
    ("1.2e-3", 0.0012, None),
])
def test_parse_number(raw, val, su):
    v, s = parse_number(raw)
    assert v == pytest.approx(val)
    if su is None:
        assert s is None
    else:
        assert s == pytest.approx(su)


@pytest.mark.parametrize("raw", [".", "?", "", "junk", None, "1.0(x)"])
def test_parse_number_missing(raw):
    assert parse_number(raw) == (None, None)


def test_su_scales_to_last_decimal():
    # 0.0234(12) means 0.0234 +/- 0.0012, not +/- 12.
    assert parse_number("0.0234(12)")[1] == pytest.approx(0.0012)
    assert parse_number("12.5(3)")[1] == pytest.approx(0.3)


@pytest.mark.parametrize("tsym,label,expect", [
    ("Zn2+", "Zn1", "Zn"), (None, "C12A", "C"), ("", "H3B", "H"),
    ("O", "O1W", "O"), (None, "FE1", "Fe"), (None, None, None),
])
def test_element_of(tsym, label, expect):
    assert element_of(tsym, label) == expect


def test_sanitise_strips_control_chars_and_caps_length():
    out = sanitise("bad\x00label\nwith\rnewlines")
    assert "\x00" not in out and "\n" not in out and "\r" not in out
    long = sanitise("A" * 500)
    assert len(long) <= 35


def test_oversized_input_refused():
    from cifaudit.load import MAX_BYTES, load
    with pytest.raises(CifError) as e:
        load("x" * (MAX_BYTES + 1))
    assert "cap" in str(e.value)


def test_empty_and_unparseable_input():
    from cifaudit.load import load
    with pytest.raises(CifError):
        load(None)
    with pytest.raises(CifError):
        load("this is not a cif at all")
