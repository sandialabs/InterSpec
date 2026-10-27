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
The executable is linked statically as far as the platform allows (only system libraries and, on
macOS, system frameworks are linked dynamically).

Using the dependency prefix the other builds use (see `target/dep_build`, which also builds
`libwttest`):
```bash
cmake -S target/mcp_stdio -B target/mcp_stdio/build -G Ninja \
      -DCMAKE_PREFIX_PATH=/path/to/dep/prefix -DCMAKE_BUILD_TYPE=Release
cmake --build target/mcp_stdio/build
```
or with `-DInterSpec_FETCH_DEPENDENCIES=ON` to have CMake fetch and build them (which builds
`libwttest` for this target; not yet tested).  On macOS, add `-DCMAKE_OSX_ARCHITECTURES=arm64` (or
`"arm64;x86_64"` for a universal binary, if the prefix is universal).

To make a self-contained directory to give to others:
```bash
cmake --install target/mcp_stdio/build --prefix /path/to/InterSpec-mcp-stdio
```
which holds:
- `bin/InterSpec-mcp-stdio`, which finds the rest by itself, relative to its own location;
- `share/interspec/data/` - InterSpec's data (nuclear data, detector responses, the tool definitions,
  and the Claude skills in `data/skills/`); local `data/llm_knowledge/` skills are left out;
- `share/interspec/InterSpec_resources/app_text/` - InterSpec's text (without it, messages in tool
  results show as `??key??` placeholders);
- `share/interspec/{README.md, LICENSE.txt, NOTICE.html}` and `share/interspec/licenses/`, which must
  accompany the executable.

A macOS executable given to others must be signed with a Developer ID and notarized (as
`.github/workflows/build_macos.yml` does for the app), or Gatekeeper will refuse to run a downloaded copy.

On Windows this needs Visual Studio 2022 17.5 or newer (QuickJS uses C11 atomics), and the static-runtime
(`/MT`) prefix from `target/dep_build/dep_build_msvc2022.bat`.  The Visual Studio generator builds several
configurations, so pass `--config Release` to `cmake --build` (a Debug build needs Debug, `/MTd`,
dependencies).

## Running
```
InterSpec-mcp-stdio [--datadir <dir>] [--userdatadir <dir>] [--llm-config-dir <dir>] [--log-file <file>]
```
- `--datadir`: InterSpec's `data` directory.  By default, that of an InterSpec installed next to the
  executable (`data/`, `resources/app/data/`, `../data/`, or `../share/interspec/data/`).  To use the
  data of an installed macOS app instead, pass `/Applications/InterSpec.app/Contents/Resources/data`.
  `InterSpec_resources/` must be next to this directory.
- `--userdatadir`: InterSpec's per-user directory, with things like `llm_config.xml`, an
  `llm_knowledge/` directory, user detector responses, and `serial_to_model.csv`.  By default, the one
  the desktop app uses, if it exists: `%APPDATA%\InterSpec` on Windows; `~/.local/share/InterSpec`,
  `~/.config/InterSpec`, or `~/.config/InterSpecAddOn` on Linux; `~/Library/Application Support/InterSpec`
  or `~/Library/Application Support/sandia.InterSpec` on macOS.  The Mac App Store app keeps its in its
  container, which is not looked in (macOS asks the user before a program reads another app's
  container); pass
  `--userdatadir "$HOME/Library/Containers/gov.sandia.macOS.InterSpec/Data/Library/Application Support/sandia.InterSpec"`
  to use it.  None is needed: the preferences database is in memory either way.
- `--llm-config-dir`: a directory with `llm_config.xml` (and `llm_tools_config.xml`, ...) to use first.
- `--log-file`: append the log here, rather than to stderr.

`<McpServer><Enabled>` in `llm_config.xml` does not need to be set: launching the program is the
opt-in.  The analysis tools do not need a `llm_config.xml` at all (InterSpec does not ship one); the
assistant tools need one with an LLM provider configured, as above.

Stdout carries only the MCP protocol; everything else InterSpec writes goes to stderr (or the log
file).  The server exits when its stdin is closed.

To use it from Claude Code (`--datadir` is not needed for an installed directory, as above):
```bash
claude mcp add interspec -- /path/to/InterSpec-mcp-stdio/bin/InterSpec-mcp-stdio
```
or in an MCP client's JSON configuration (e.g., Claude Desktop's `claude_desktop_config.json`):
```json
{ "mcpServers": { "interspec": { "command": "/path/to/InterSpec-mcp-stdio/bin/InterSpec-mcp-stdio",
                                 "args": [] } } }
```

The server gives the client brief usage instructions (`data/mcp_stdio_instructions.md`) when it
connects.  For much better analyses, also give Claude InterSpec's skills - copy (or symlink) the
`interspec-*` directories of `share/interspec/data/skills/`, except `interspec-mcp-agent` (which is
for the multi-session server), into `~/.claude/skills/`:
```bash
for s in /path/to/InterSpec-mcp-stdio/share/interspec/data/skills/interspec-*; do
  [ "$(basename "$s")" = interspec-mcp-agent ] || cp -R "$s" ~/.claude/skills/
done
```

Spectra are loaded with `load_spectrum_file` and an absolute `path`.  Tool descriptions can differ
between the MCP servers: they are Inja templates (see `InterSpec/LlmPromptTemplate.h`), and
`mcp_server` is `"stdio"` here.

## Running LLM benchmarks
With `--run-benchmark`, the program is not an MCP server: it does not read stdin, and instead runs
`*_llm_benchmark.xml` question sets through InterSpec's LLM assistant - as the "Benchmarks / Quizzes"
menu of the assistant does in the GUI - prints a summary of each to stdout, and exits.
```
InterSpec-mcp-stdio [--datadir ...] [--userdatadir ...] [--llm-config-dir ...] [--log-file ...]
    --run-benchmark <xml> [--run-benchmark <xml> ...] [--skills all|none|<name>[,<name>...]]
    [--model <name>] [--judge-model <name>] [--output-dir <dir>]
InterSpec-mcp-stdio [--datadir ...] [--userdatadir ...] --list-benchmarks
```
- `--run-benchmark`: a benchmark XML; may be given more than once, to run several in order, in one session.
- `--skills`: the KnowledgeBase skills the assistant gets: `all` (the default), `none` (no KnowledgeBase
  agent at all, even with a RAG endpoint configured), or a comma-separated list of skill names (an
  unknown name is an error that lists the available ones).
- `--model`, `--judge-model`: the model to test, and the model that grades the answers, for this run
  only; both must be served by the active `<ApiProvider>` of `llm_config.xml` (they need not be listed
  in it).  By default the model tested is the active one, and the judge is the `<Model judge="true">`,
  or else the model being tested.  `llm_config.xml` is never written.
- `--output-dir`: where to write results; by default `<userdatadir>/llm_benchmarks`.  Required when
  there is no user data directory, so results never go into an installed data directory.
- `--list-benchmarks`: lists the benchmarks found (in `test_datasets/` and the KnowledgeBase skill
  directories, as the GUI menu does), one `<name><tab><path>` per line, and exits.

A run always starts fresh (an interrupted run's `_in_progress.json` is left alone; the GUI can resume
it), and shows no dialogs.  It needs an LLM provider configured in `llm_config.xml`, but not
`<McpServer><AllowAssistantPrompt>`.  For each benchmark the output directory gets
`<base>_<time>_<hash>_results.json`, and a `<base>_<time>_<hash>_logs/` directory with the conversation
of every attempt at every question (`<problem>_q<part>_try<N>.xml`, failed attempts included); the
results JSON lists each question's logs (`conversationLogs`), and the tool calls of its final attempt,
sub-agents' included (`toolCalls`).

Stdout carries only the summaries: for each benchmark, its name, model, judge model, skills, score,
correct and error counts, a line per question (score, error/skip marker, retries, tokens), and the path
of the results JSON.  Progress goes to stderr (or the `--log-file`).  Ctrl-C cancels the run and prints
what was done so far (the checkpoint is written after every question).

Exit code: 0 if every benchmark ran to the end (wrong answers and per-question errors are part of the
results), 1 for bad arguments or a setup problem, and 2 if a benchmark could not be run (or was
interrupted).

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
