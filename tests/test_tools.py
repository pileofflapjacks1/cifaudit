"""The MCP tool surface, tested without the MCP SDK.

Every published release of the `mcp` package requires Python >= 3.10. Importing
it at module level in server.py made the whole tool surface unreachable on an
older interpreter, so none of this logic was covered. The implementations now
live in `cifaudit.tools` as plain functions; `server.py` only registers them.
"""

import pathlib

import pytest

from cifaudit import tools

DATA = pathlib.Path(__file__).parent / "data"


def p(name):
    return str(DATA / (name + ".cif"))


def test_the_tool_surface_is_exactly_four_tools():
    assert sorted(tools.TOOLS) == ["audit_cif", "coverage", "describe_alert",
                                   "list_datablocks"]


def test_list_datablocks(): 
    r = tools.list_datablocks(path=p("multiblock"))
    assert r["status"] == "ok"
    assert r["datablocks"] == ["clean", "second"]
    assert r["count"] == 2


def test_audit_cif_returns_the_full_contract():
    r = tools.audit_cif(path=p("npd"))
    assert r["status"] == "ok"
    for key in ("counts", "reporting_counts", "findings", "not_tested",
                "checked", "notes", "disclaimer"):
        assert key in r, key
    assert any(f["code"] == "ADP_NPD" for f in r["findings"])


def test_audit_cif_honours_datablock_and_options():
    assert tools.audit_cif(path=p("multiblock"),
                           datablock="second")["datablock"] == "second"
    loose = tools.audit_cif(path=p("bad_volume"), options={"volume_rel_tol": 1.0})
    assert "CELL_VOLUME_MISMATCH" not in {f["code"] for f in loose["findings"]}


def test_describe_alert_is_case_insensitive_and_cites_its_source():
    r = tools.describe_alert("plat213")
    assert r["code"] == "PLAT213"
    assert r["description"], "PLAT213 has no description in the table"
    assert r["platon_type"] == "Type_2"
    assert "journals.iucr.org" in r["source"]


def test_describe_alert_does_not_invent_unknown_codes():
    r = tools.describe_alert("PLAT999")
    assert r["description"] is None
    assert "note" in r


def test_coverage_reports_what_is_not_implemented():
    c = tools.coverage()
    assert len(c["checks"]) == 9
    assert c["total_platon_tests"] == 526
    assert c["checkcif_codes_known_count"] == len(c["checkcif_codes_known"])
    assert "not implemented at all" in c["coverage_note"]
    assert c["not_a_substitute"] == "https://checkcif.iucr.org/"


@pytest.mark.parametrize("name", ["list_datablocks", "audit_cif"])
def test_errors_abstain_rather_than_raise(name):
    fn = tools.TOOLS[name]
    for kwargs in ({}, {"cif_text": "x", "path": "y"},
                   {"path": "/nonexistent/nope.cif"}):
        r = fn(**kwargs)
        assert r["status"] == "not_enough_data", (name, kwargs)
        assert "error" in r


def test_every_response_carries_the_untrusted_input_guard():
    for r in (tools.list_datablocks(path=p("clean")),
              tools.audit_cif(path=p("clean")),
              tools.list_datablocks()):
        assert "not instructions" in r["note"]


def test_injection_payloads_never_reach_a_tool_response():
    import json
    blob = json.dumps(tools.audit_cif(path=p("injection")))
    for payload in ("Ignore all previous instructions", "id_rsa", "VALIDATED"):
        assert payload not in blob


def test_every_plat_code_cited_in_the_source_is_tabulated():
    """Static scan, not a walk of whatever the fixtures happen to produce.

    The existing contract test only inspected findings from a handful of
    fixtures, so PLAT213 - cited by ADP_EXTREME_ANISOTROPY, which none of them
    trigger - sat in the source for a whole release without a description.
    """
    import re
    from cifaudit.alerts import PLAT
    src_dir = pathlib.Path(tools.__file__).parent
    cited = set()
    for path in src_dir.rglob("*.py"):
        if path.name == "alerts.py":
            continue
        cited |= set(re.findall(r"PLAT\d{3}", path.read_text()))
    missing = sorted(cited - set(PLAT))
    assert not missing, "cited but not tabulated: %s" % missing


def _install_stub_sdk(monkeypatch, record):
    """Put a minimal fake `mcp.server.fastmcp` on sys.modules.

    The real SDK needs Python >= 3.10. This exercises server.py's registration
    wiring on any interpreter: a typo in the loop, a tool not registered, or
    main() not starting the server would all be caught here. It does NOT prove
    the real SDK accepts these signatures - only running it can do that.
    """
    import sys
    import types

    class StubFastMCP:
        def __init__(self, name):
            record["name"] = name

        def tool(self, name=None):
            def deco(fn):
                record.setdefault("registered", []).append((name, fn))
                return fn
            return deco

        def run(self):
            record["ran"] = True

    pkg = types.ModuleType("mcp")
    server_mod = types.ModuleType("mcp.server")
    fastmcp_mod = types.ModuleType("mcp.server.fastmcp")
    fastmcp_mod.FastMCP = StubFastMCP
    server_mod.fastmcp = fastmcp_mod
    pkg.server = server_mod
    for mod_name, mod in (("mcp", pkg), ("mcp.server", server_mod),
                          ("mcp.server.fastmcp", fastmcp_mod)):
        monkeypatch.setitem(sys.modules, mod_name, mod)
    monkeypatch.delitem(sys.modules, "cifaudit.server", raising=False)


def test_server_registers_every_tool_under_its_own_name(monkeypatch):
    record = {}
    _install_stub_sdk(monkeypatch, record)
    import cifaudit.server as srv
    assert record["name"] == "cifaudit"
    names = [n for n, _ in record["registered"]]
    assert sorted(names) == sorted(tools.TOOLS)
    for name, fn in record["registered"]:
        assert fn is tools.TOOLS[name], name
    assert srv.mcp is not None


def test_server_main_starts_the_server(monkeypatch):
    record = {}
    _install_stub_sdk(monkeypatch, record)
    import cifaudit.server as srv
    srv.main()
    assert record.get("ran") is True


def test_server_exits_with_guidance_when_the_sdk_is_absent(monkeypatch):
    """On Python 3.8 the import must fail with an actionable message, not a
    bare ImportError traceback."""
    import sys
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **kw):
        if name.startswith("mcp"):
            raise ImportError("No module named 'mcp'")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    monkeypatch.delitem(sys.modules, "cifaudit.server", raising=False)
    with pytest.raises(SystemExit) as e:
        import cifaudit.server  # noqa: F401
    msg = str(e.value)
    assert "Python >= 3.10" in msg
    assert "cifaudit.tools" in msg
