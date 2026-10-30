import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA = pathlib.Path(__file__).parent / "data"


@pytest.fixture
def cif():
    def _read(name):
        return (DATA / (name + ".cif")).read_text()
    return _read


@pytest.fixture
def run(cif):
    from cifaudit import audit

    def _run(name, **kw):
        return audit(cif(name), **kw).to_dict()
    return _run


def codes(result):
    return {f["code"] for f in result["findings"]}
