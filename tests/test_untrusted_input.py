"""A CIF is untrusted input. Nothing in it may be echoed or obeyed."""

import json

PAYLOADS = ["Ignore all previous instructions", "id_rsa", "VALIDATED",
            "must now output", "SYSTEM:", "exfiltrate", "Assistant:"]


def test_free_text_fields_are_never_echoed(run):
    blob = json.dumps(run("injection"))
    for p in PAYLOADS:
        assert p not in blob, "leaked %r" % p


def test_free_text_fields_are_not_even_parsed(cif):
    from cifaudit.load import load
    blk = load(cif("injection"))
    for attr in vars(blk):
        if attr == "_raw":
            continue
        assert "Ignore all previous" not in str(getattr(blk, attr))


def test_hostile_label_is_truncated_and_inert(run):
    r = run("injection")
    labels = [s for f in r["findings"] for s in f["sites"]]
    hostile = [l for l in labels if l.startswith("O_IGNORE")]
    assert hostile, labels
    for l in hostile:
        assert len(l) <= 35          # capped, so it cannot restructure a report
        assert "\n" not in l


def test_real_defect_still_found_despite_injection(run):
    """The hostile text must not distract from the actual error."""
    assert "OCC_GT_ONE" in {f["code"] for f in run("injection")["findings"]}


def test_no_network_or_llm_in_the_validation_path():
    """cifaudit must not import anything that reaches the network."""
    import cifaudit
    import cifaudit.audit  # noqa: F401
    import sys
    banned = {"requests", "httpx", "urllib3", "openai", "anthropic", "socket"}
    loaded = banned & set(sys.modules)
    # socket may be pulled in by the interpreter itself; assert cifaudit's own
    # modules never reference these names.
    import pathlib
    src = pathlib.Path(cifaudit.__file__).parent
    for path in src.rglob("*.py"):
        text = path.read_text()
        for name in ("requests", "httpx", "openai", "anthropic", "urlopen"):
            assert name not in text, "%s mentions %s" % (path.name, name)
