"""MCP server exposing cifaudit. A thin registration layer over `tools.py`.

Requires the `mcp` Python SDK, which needs Python >= 3.10 - every published
release does. The validation core, the CLI and all four tool implementations
work on 3.8+; only this file needs the newer interpreter.

    pip install "cifaudit[mcp]"
    cifaudit-mcp                    # stdio

Design notes:
* Tools are deterministic. No network access, no model in the validation path.
* `list_datablocks` should be called first on an unfamiliar file, so the wrong
  datablock is never audited silently.
* Every `audit_cif` response carries `not_tested` and a disclaimer. A clean
  report is not a statement that the structure is correct.
* CIF content is untrusted data, never instructions.

The logic lives in `cifaudit.tools` precisely so that it is testable without the
SDK present; this module contains no behaviour of its own.
"""

from . import tools

try:
    from mcp.server.fastmcp import FastMCP
except ImportError as exc:  # pragma: no cover - depends on the interpreter
    raise SystemExit(
        "The MCP server needs the 'mcp' package, which requires Python >= 3.10:\n"
        "    pip install \"cifaudit[mcp]\"\n"
        "The cifaudit library and CLI work without it on Python 3.8+, and the\n"
        "tool implementations are importable from cifaudit.tools.\n"
        "(%s)" % exc)

mcp = FastMCP("cifaudit")

for _name, _fn in tools.TOOLS.items():
    mcp.tool(name=_name)(_fn)


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
