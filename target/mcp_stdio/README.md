# InterSpec as a stdio MCP server

This target builds `InterSpec-mcp-stdio`, an [MCP](https://modelcontextprotocol.io) server that an MCP
client (Claude Code, Claude Desktop, Cursor, ...) launches as a subprocess and talks to over stdin and
stdout.  It hosts a single InterSpec session of its own, with no window, no browser and no web server,
and offers InterSpec's analysis tools: loading spectra, fitting peaks, identifying sources, activity and
shielding fits, isotopics, detection limits, nuclear data, and so on.

It is meant to be installed next to an InterSpec desktop app, whose data directory and per-user
directory it uses.  It does not change anything there: its preferences database is in memory, so it
never contends with (or upgrades) the database of a running InterSpec.

What is here:
- `CMakeLists.txt` - the build: sets `BUILD_AS_MCP_STDIO_SERVER` (and the options it implies), adds the
  InterSpec library, and builds the executable.
- `mcp_stdio_main.cpp` - its `main()`.
- `PlutoPaintDevice.{h,cpp}` - draws charts to PNG/JPEG with plutovg.
- `QuickJsSandbox.{h,cpp}` - runs the assistant's `run_javascript` code with QuickJS.
- `embed_file.cmake` - compiles the font into the executable.
- `tests/stdio_test.py` - a test client (standard-library Python).
- `tests/test_quickjs_sandbox.cpp` - tests of the QuickJS sandbox (`ctest` in the build directory).

The build fetches plutovg (v1.3.3, MIT), quickjs-ng (v0.17.0, MIT), and the DejaVu fonts (2.37,
Bitstream Vera license), so none is added to the main InterSpec build; the font's license is installed to
`share/interspec/licenses/DejaVuFonts-LICENSE.txt`, and must accompany the executable.

The MCP message handling (`InterSpec/LlmMcpDispatcher.h`) is shared with the HTTP MCP endpoint of the
desktop builds and with the multi-session server (`target/mcp_http`).

## How it works
The session is hosted by Wt's test environment (`Wt::Test::WTestEnvironment`, from `libwttest`), which
creates a real InterSpec session without a browser; so the executable links `libwttest` instead of
Wt's HTTP connector, and builds its own copy of the InterSpec library rather than sharing the desktop
app's.  Tools run in the session on Wt's IO threads, as they do for the HTTP server.

The spectrum and time-chart images (`get_spectrum_image`, `get_time_chart_image`) are drawn by the
server: SVG with Wt's `WSvgImage`, and PNG/JPEG with [plutovg](https://github.com/sammycage/plutovg)
(`PlutoPaintDevice`), using the DejaVu Sans font compiled into the executable - so labels are in
English, whatever the user's language.  The chart drawing itself (`InterSpec/SpectrumChartImage.h`) is
shared with the macOS Quick Look extension.

InterSpec's built-in LLM assistant is offered (`assistant_submit_prompt`, `assistant_get_result`,
`assistant_cancel`, and `invoke_agent`) when `llm_config.xml` enables the LLM and sets
`<McpServer><AllowAssistantPrompt>` (or `<AllowAgentInvocation>`).  Its LLM requests are always made
from C++ (whatever `<HttpBackend>` says): NSURLSession on macOS, WinHTTP on Windows, and on Linux
libcurl - loaded at run time, so any installed `libcurl.so.4` (7.61 or newer) works, and the
executable does not depend on one.  libcurl uses `<HttpProxy>` if set, and otherwise the `https_proxy`,
`ALL_PROXY` and `NO_PROXY` environment variables; `<HttpCaBundle>` is added to the system's trust store.  Its `run_javascript`
tool runs in [QuickJS](https://github.com/quickjs-ng/quickjs) (`QuickJsSandbox`) rather than a browser
iframe, with the same results; each call gets a fresh engine limited to 64 MB, a 256 kB stack and its
time limit, with no modules, files, network or timers.

## Building
On Linux and Windows the executable is linked statically as far as the platform allows (only system
libraries are linked dynamically), so it can be copied next to an installed InterSpec.  On macOS it is
only built for local use.

Using the dependency prefix the other builds use (see `target/dep_build`, which also builds
`libwttest`):
```bash
cmake -S target/mcp_stdio -B target/mcp_stdio/build -G Ninja \
      -DCMAKE_PREFIX_PATH=/path/to/dep/prefix -DCMAKE_BUILD_TYPE=Release
cmake --build target/mcp_stdio/build
```
or with `-DInterSpec_FETCH_DEPENDENCIES=ON` to have CMake fetch and build them (which builds
`libwttest` for this target; not yet tested).

## Running
```
InterSpec-mcp-stdio [--datadir <dir>] [--userdatadir <dir>] [--llm-config-dir <dir>] [--log-file <file>]
```
- `--datadir`: InterSpec's `data` directory.  By default, that of an InterSpec installed next to the
  executable (`data/`, `resources/app/data/`, `../data/`, or `../share/interspec/data/`).
- `--userdatadir`: InterSpec's per-user directory, with things like `llm_config.xml`, an
  `llm_knowledge/` directory, user detector responses, and `serial_to_model.csv`.  By default, the one
  the desktop app uses, if it exists: `%APPDATA%\InterSpec` on Windows; `~/.local/share/InterSpec`,
  `~/.config/InterSpec`, or `~/.config/InterSpecAddOn` on Linux; `~/Library/Application Support/InterSpec`
  on macOS.
- `--llm-config-dir`: a directory with `llm_config.xml` (and `llm_tools_config.xml`, ...) to use first.
- `--log-file`: append the log here, rather than to stderr.

`<McpServer><Enabled>` in `llm_config.xml` does not need to be set: launching the program is the
opt-in.  The configuration is still read for the tool definitions.

Stdout carries only the MCP protocol; everything else InterSpec writes goes to stderr (or the log
file).  The server exits when its stdin is closed.

To use it from Claude Code:
```bash
claude mcp add interspec-stdio -- /path/to/InterSpec-mcp-stdio --datadir /path/to/InterSpec/data
```
or in an MCP client's JSON configuration:
```json
{ "mcpServers": { "interspec": { "command": "/path/to/InterSpec-mcp-stdio",
                                 "args": ["--datadir", "/path/to/InterSpec/data"] } } }
```

Spectra are loaded with `load_spectrum_file` and an absolute `path`.

## Testing
```bash
python3 -I target/mcp_stdio/tests/stdio_test.py --exe target/mcp_stdio/build/InterSpec-mcp-stdio \
    --datadir data --spectrum example_spectra/ba133_source_640s_20100317.n42 --energy 356
```
Add `--background example_spectra/background_20100317.n42 --time-spectrum example_spectra/passthrough.n42`
to also check the chart images (and `--save-images <dir>` to look at them), and
`--repeat 50 --alt-spectrum example_spectra/passthrough.n42` to load and fit repeatedly, reporting the
server's memory use.

`tests/stdio_torture.py` tries to break things: `sweep` calls every tool (with the call tables of
`target/mcp_http/tests/torture_test.py`) and checks none hangs - a tool that waits on a browser would;
`chaos` sends many concurrent, malformed, split, and oversized messages, stops reading stdout for a
while, and closes stdin or sends SIGTERM mid-calculation; `llm` runs heavy assistant turns, one
server process each, in parallel (needs `--llm-config-dir`, with `<AllowAssistantPrompt>`):
```bash
python3 -I target/mcp_stdio/tests/stdio_torture.py --exe target/mcp_stdio/build/InterSpec-mcp-stdio \
    --out /tmp/stdio_torture --phases sweep,chaos
```

When stdin closes, the server exits once the session is torn down, but waits at most 5 s for a
calculation still running (nothing needs saving: its database is in memory).
