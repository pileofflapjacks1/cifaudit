"""MCP server exposing cifaudit. A thin registration layer over `tools.py`.

Requires the `mcp` Python SDK, which needs Python >= 3.10 - every published
release does. The validation core, the CLI and all four tool implementations
work on 3.8+; only this file needs the newer interpreter.

    pip install "cifaudit[mcp]"
    cifaudit-mcp                    # stdio

Both SDK generations are supported. `mcp` 2.0 renamed `FastMCP` to `MCPServer`
and moved it to `mcp.server.mcpserver`; the v1 path is kept as a fallback
because plenty of environments still have it. Their constructor, `tool`
decorator and `run` method are shape-compatible for what this server needs.

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

SDK_API = None

try:                                            # mcp >= 2
    from mcp.server.mcpserver import MCPServer as _Server
    SDK_API = "mcp>=2 (mcp.server.mcpserver.MCPServer)"
except ImportError:
    try:                                        # mcp 1.x
        from mcp.server.fastmcp import FastMCP as _Server
        SDK_API = "mcp<2 (mcp.server.fastmcp.FastMCP)"
    except ImportError as exc:                  # pragma: no cover
        raise SystemExit(
            "The MCP server needs the 'mcp' package, which requires "
            "Python >= 3.10:\n"
            "    pip install \"cifaudit[mcp]\"\n"
            "The cifaudit library and CLI work without it on Python 3.8+, and\n"
            "the tool implementations are importable from cifaudit.tools.\n"
            "(%s)" % exc)

mcp = _Server(name="cifaudit")

for _name, _fn in tools.TOOLS.items():
    mcp.tool(name=_name)(_fn)


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
