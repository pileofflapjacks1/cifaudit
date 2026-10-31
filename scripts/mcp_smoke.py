"""End-to-end check: spawn cifaudit-mcp over stdio and call all four tools.

Needs the real SDK and therefore Python >= 3.10, so it is a script rather than a
pytest test. The unit tests cover the tool implementations and register them
against a STUB SDK; a stub verifies the code against the author's assumption of
the dependency, not the dependency. That is exactly how cifaudit shipped a
server written for `mcp` 1.x's FastMCP after 2.0 renamed it to MCPServer - the
stub obligingly provided `mcp.server.fastmcp.FastMCP` because that is what the
code imported. Only this script would have caught it.

    uv venv --python 3.12 /tmp/v312
    uv pip install --python /tmp/v312/bin/python --only-binary cryptography \
        -e ".[mcp,symmetry]"
    /tmp/v312/bin/python scripts/mcp_smoke.py \
        /tmp/v312/bin/cifaudit-mcp tests/data
"""
import asyncio, json, sys
from mcp import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

SERVER = sys.argv[1]
DATA = sys.argv[2]

def text_of(result):
    for c in result.content:
        if getattr(c, "type", None) == "text":
            return c.text
    return None

async def main():
    params = StdioServerParameters(command=SERVER, args=[])
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            init = await session.initialize()
            print("initialize OK -> server:", init.server_info.name,
                  "| protocol:", init.protocol_version)
            tools = await session.list_tools()
            print("list_tools  OK ->", sorted(t.name for t in tools.tools))
            print()

            r = await session.call_tool("list_datablocks",
                                        {"path": DATA + "/multiblock.cif"})
            d = json.loads(text_of(r))
            print("list_datablocks ->", d["status"], d["datablocks"])

            r = await session.call_tool("audit_cif", {"path": DATA + "/npd.cif"})
            d = json.loads(text_of(r))
            print("audit_cif       ->", d["status"], d["counts"],
                  "| findings:", [f["code"] for f in d["findings"]])
            print("                   not_tested areas:", len(d["not_tested"]))
            print("                   disclaimer:", d["disclaimer"][:46], "...")

            r = await session.call_tool("describe_alert", {"code": "PLAT213"})
            d = json.loads(text_of(r))
            print("describe_alert  ->", d["code"], d["platon_type"], "|",
                  d["description"][:52])

            r = await session.call_tool("coverage", {})
            d = json.loads(text_of(r))
            print("coverage        ->", len(d["checks"]), "checks,",
                  d["checkcif_codes_known_count"], "of",
                  d["total_platon_tests"], "PLATON tests")

            r = await session.call_tool("audit_cif", {})
            d = json.loads(text_of(r))
            print("abstention      ->", d["status"], "|", d["error"])

            r = await session.call_tool("audit_cif",
                                        {"path": DATA + "/injection.cif"})
            blob = text_of(r)
            leaked = [p for p in ("Ignore all previous instructions", "id_rsa",
                                  "VALIDATED") if p in blob]
            print("injection       -> leaked payloads:", leaked or "none")
            print()
            print("ALL FOUR TOOLS CALLED SUCCESSFULLY OVER STDIO")

asyncio.run(main())
