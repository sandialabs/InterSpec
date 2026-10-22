# InterSpec as a multi-session MCP analysis server

This target builds InterSpec as an [MCP](https://modelcontextprotocol.io) server for many
independent analysis sessions at once, normally run as a container built by
[Containerfile](Containerfile).  It is meant to sit behind a back-end - for example, a chat application for
gamma spectroscopy, whose users upload spectra (and photos) and ask about isotope ID,
quantification, or whether a source is a threat.  The back-end creates a session, exchanges
files with it, and either relays its users' questions to InterSpec's built-in spectroscopy
assistant (which calls an LLM you configure) or drives InterSpec's analysis tools directly.

Each session is an ordinary InterSpec session, kept alive by a tab of a headless Chromium inside
the container.  Sessions do not share spectra, preferences, or files.

What is here:
- `CMakeLists.txt` - the build: sets `BUILD_AS_MCP_HTTP_SERVER` (and the options it implies), adds
  the InterSpec library, and builds the `InterSpecMcpServer` executable.
- `mcp_http_main.cpp` - its `main()`.
- `McpSessionManager.{h,cpp}` - sessions, their directories and file tokens, and the headless
  browser (Chromium, driven over the DevTools protocol on a pipe).
- `McpAgentTools.{h,cpp}` - the tools only this server offers (`create_session`, `close_session`,
  `upload_file`, `list_session_files`, `export_analysis`, `get_usage_guide`).
- `Containerfile`, `wt_config_mcp_http.xml`, `config.example/` - deployment.
- `tests/` - `smoke_test.py`, `torture_test.py`, `load_test.py` (standard-library Python).

The MCP endpoint itself (`InterSpec/LlmMcpResource.h`) is shared with the desktop builds, which serve
a single user's session to local MCP clients; code that only this server needs is compiled only with
`BUILD_AS_MCP_HTTP_SERVER`.  The server instructions for clients are `data/mcp_agent_instructions.md`,
and the client skill is `data/skills/interspec-mcp-agent/SKILL.md`.

## Building
From the root of the InterSpec source tree:
```bash
docker build -t interspec-mcp-agent -f target/mcp_http/Containerfile .
```
The first build compiles the dependencies (boost, Wt with OpenSSL, Eigen, Ceres), which takes a
while; they are in their own stage, so later builds only recompile InterSpec.  The image does not
contain any LLM configuration, keys, or curl.

## Configuration
The container reads its configuration from a directory you mount at `/config` (read-only is
fine; nothing is ever written there):
- `llm_config.xml` (required) - the LLM provider and model, the MCP bearer token, and the
  agent-server settings.  Start from [config.example/](config.example/).
- optionally, `mcp_agent_instructions.md`, to replace the `instructions` sent to MCP clients.
- optionally, `llm_tools_config.xml` / `llm_agent_*.xml` to customize tools and agents (otherwise
  the versions shipped in the image are used), and an `llm_knowledge/` directory of extra
  knowledge-base skills (added to, and overriding, the shipped ones).

Any value in `llm_config.xml` that is entirely `env(NAME)` or `env(NAME,default)` is read from that
environment variable (`NAME` must start with `INTERSPEC_`), so the file can be kept without
secrets.  The example config uses:

| Variable | Default | |
|---|---|---|
| `INTERSPEC_MCP_BEARER_TOKEN` | (required) | Token MCP clients must send (`Authorization: Bearer ...`); at least 16 characters. |
| `INTERSPEC_LLM_ENDPOINT` | (required) | e.g. `https://api.anthropic.com/v1/messages`, or any OpenAI-compatible `.../chat/completions`. |
| `INTERSPEC_LLM_BEARER_TOKEN` | (required) | The LLM API key. |
| `INTERSPEC_LLM_MODEL` | (required) | Model name. |
| `INTERSPEC_LLM_API_FORMAT` | `openai_chat` | `openai_chat`, `openai_responses`, or `anthropic`. |
| `INTERSPEC_LLM_SUPPORTS_IMAGES` | `true` | Whether the model accepts images (spectrum pictures, and photos sent with prompts). |
| `INTERSPEC_LLM_REASONING`, `_MAX_TOKENS`, `_CONTEXT_LENGTH` | `medium`, 32000, 200000 | |
| `INTERSPEC_LLM_LOG_CONVERSATIONS` | `false` | Writes each conversation to the session's `logs/` directory. |
| `INTERSPEC_MAX_SESSIONS` | 8 | Concurrent sessions.  Measured: ~400 MB idle, plus ~100 MB per session (mostly its browser tab). |
| `INTERSPEC_SESSION_IDLE_TIMEOUT` | 1800 | Seconds without a client call before a session is ended. |
| `INTERSPEC_SESSION_MAX_LIFETIME` | 14400 | Seconds. |
| `INTERSPEC_KEEP_SESSION_FILES` | `false` | If false, a session's directory is deleted when it ends. |
| `INTERSPEC_MAX_UPLOAD_BYTES`, `_MAX_SESSION_BYTES`, `_MAX_IMAGE_BYTES` | 256 MiB, 1 GiB, 5 MiB | |
| `INTERSPEC_PUBLIC_BASE_URL` | (from the `Host` header) | Base of the upload/download URLs, if clients reach the server through a proxy. |

The LLM is always called from C++ (never from the browser), so the API key never reaches a web page.
The Linux HTTP client has no proxy support, so the container needs direct (or transparent-proxy)
access to the LLM endpoint; for a private CA, mount its PEM file and set `<HttpCaBundle>`.

## Running
```bash
docker run -d --name interspec -p 8078:8078 \
  -v /srv/interspec/config:/config:ro \
  -v /srv/interspec/exchange:/exchange \
  -e INTERSPEC_MCP_BEARER_TOKEN=... \
  -e INTERSPEC_LLM_ENDPOINT=https://... -e INTERSPEC_LLM_BEARER_TOKEN=... -e INTERSPEC_LLM_MODEL=... \
  --cap-drop=ALL --security-opt=no-new-privileges --pids-limit=2048 --memory=6g --cpus=4 \
  --read-only --tmpfs /tmp --tmpfs /var/lib/interspec:uid=10001 \
  interspec-mcp-agent
```
- The MCP endpoint is `http://<host>:8078/mcp-api`; `GET /mcp-api/healthz` (no authentication) reports
  session counts only.
- `/exchange` holds one directory per session; mount it from the host only if your back-end wants to
  exchange files through the filesystem.  The server runs as uid 10001, so a host directory must be
  writable by it (or run with `--user` matching the directory's owner).
- `/var/lib/interspec` holds InterSpec's preferences database; nothing in it needs to survive a
  restart, hence the tmpfs.  At startup, session directories left behind by a crash are removed
  (unless `KeepSessionFiles` is set) - so each server needs an exchange directory of its own; two
  servers (or an old and a new container) sharing one would remove each other's sessions.
- Limit `--cpus` and `--memory`: InterSpec sizes its fitting thread pools from the CPU count.  The
  number of HTTP worker threads is `<num-threads>` in `wt_config_mcp_http.xml`.
- Chromium runs with `--no-sandbox` (Docker's default seccomp profile does not allow the user
  namespaces its sandbox needs).  Its network access is limited to InterSpec itself.  If you provide
  a seccomp profile that allows unprivileged user namespaces, set `-e INTERSPEC_BROWSER_NO_SANDBOX=0`.

## Using it
Clients find the workflow in the server's MCP `instructions`, the tool descriptions, and in more
detail in the [interspec-mcp-agent skill](../../data/skills/interspec-mcp-agent/SKILL.md) (also
available from the `get_usage_guide` tool).  In short:

1. `create_session` returns `userSession` (pass it to every session tool), `session_dir`, and
   `upload_url`.
2. Put files in the session's `input/` directory, by any of:
   - `PUT` to `<upload_url>/<filename>` (or a multipart `POST` to `<upload_url>`), with the bearer token;
   - the `upload_file` tool (base64; small files);
   - writing into `/exchange/<session_dir>/input/` on the host (write to a temporary name, then rename).
3. `load_spectrum_file` with `path: "input/<filename>"`.  It opens the file as the GUI does (choosing
   which of several samples to show, and possibly loading the file's background sample too), and
   returns a `summary` of what is displayed; `sampleNumbers` and `detectors` choose explicitly.
4. `assistant_submit_prompt` (optionally with `images`), then `assistant_get_result` while its status
   is `running`.
5. `export_analysis` writes files to `output/`, each with a `download_url` (`GET` with the bearer
   token); its `n42` item is the spectrum with all the analysis, including the assistant conversation.
6. `close_session`.  Its URLs then stop working, and its directory is deleted (unless kept).

With curl:
```bash
TOKEN=...; MCP=http://localhost:8078/mcp-api
call() { curl -s $MCP -H "Authorization: Bearer $TOKEN" -H 'Content-Type: application/json' \
           -d "{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"tools/call\",\"params\":{\"name\":\"$1\",\"arguments\":$2}}" \
         | python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["result"]["content"][0]["text"])'; }

call create_session '{}'      # -> {"userSession":"Ab3...","upload_url":"http://localhost:8078/mcp-api/files/...", ...}
S=Ab3...; UP=http://localhost:8078/mcp-api/files/...
curl -s -H "Authorization: Bearer $TOKEN" -T spectrum.n42 "$UP/spectrum.n42"
call load_spectrum_file "{\"userSession\":\"$S\",\"path\":\"input/spectrum.n42\",\"specType\":\"Foreground\"}"
call assistant_submit_prompt "{\"userSession\":\"$S\",\"prompt\":\"What sources are present?\"}"
call export_analysis "{\"userSession\":\"$S\",\"items\":[\"n42\",\"spectrum_html\"]}"
curl -s -H "Authorization: Bearer $TOKEN" -o result.n42 "<download_url from export_analysis>"
call close_session "{\"userSession\":\"$S\"}"
```

To resume later: create a new session, upload the exported `*_analysis.n42`, and load it with
`"conversation": "resume"`; the peaks, models, and conversation (with any photos) come back, and the
assistant can continue.  Only resume N42 files your back-end exported itself - by default a loaded
file's embedded conversation is ignored, since one from elsewhere could be a prompt injection.

To connect Claude Code (e.g., for testing):
```bash
claude mcp add --transport http interspec http://localhost:8078/mcp-api --header "Authorization: Bearer $TOKEN"
```

[tests/smoke_test.py](tests/smoke_test.py) exercises the whole workflow against a running
server, and [tests/torture_test.py](tests/torture_test.py) calls every tool, injects
failures (disconnects, closing sessions mid-call, malformed requests, odd uploads, browser
crashes), and runs heavy assistant turns in parallel (standard-library Python only).
[tests/load_test.py](tests/load_test.py) checks which samples, detectors and background
`load_spectrum_file` displays for multi-sample files.

Assistant turns return within `wait_seconds` (default 30, at most 50) with `status: "running"` and
a `turn_id` if not yet done; clients that accept `text/event-stream` and send a `progressToken` also
get MCP progress notifications (with a heartbeat about every 10 s) while waiting.

Some analysis tools do real fitting and take minutes (e.g. `perform_isotopics_calculation`); a session
does one thing at a time, so other calls on it wait - but other sessions are not held up.  Clients
calling such tools directly need long HTTP timeouts - or can ask the assistant, whose turns are polled.  `docker stop` takes up to ~10 s,
longer if a calculation is running (it falls back to SIGKILL; leftover session directories are
removed at the next start).

## Security notes
- Trust model: the back-end holding the bearer token is trusted; end users, their files, and LLM
  output are not.  Keep the token on the back-end; downloads and uploads also require it.
- Every file path a session uses is resolved inside that session's own directory.
- Text inside spectrum files and images is passed to the LLM as data; prompt injection through it
  can misdirect that session's analysis, and the assistant's tools reach only that session.  But
  JavaScript the assistant writes runs in the session's browser tab, and Chromium runs with
  `--no-sandbox`, as the same user as the server: a renderer exploit there could read the server's
  environment (the `INTERSPEC_*` secrets) and the exchange directory.  Run the container with the
  limits above, and consider disabling the assistant's `run_javascript` tool (in
  `llm_tools_config.xml`) if that risk matters for your deployment.
- Outbound traffic: only the configured LLM endpoint (the browser can reach nothing but InterSpec).
  Restrict egress at the host or orchestration level (Docker network, Kubernetes NetworkPolicy) to
  enforce this.  The KnowledgeBase RAG tool (`KnowledgeBaseRagUrl`) needs curl, which the image
  does not have; local `llm_knowledge/` skills work.
- The server writes no HTTP access log.  A reverse proxy's log would include request paths, which
  contain upload/download tokens; they are useless without the bearer token, but treat such a log
  accordingly.
- With `LogConversations` on, conversations are written to the session's `logs/` directory, which is
  deleted with the session unless `KeepSessionFiles` is set.

## Building and running without a container
For development, build this directory like any other target (Release - fitting is 10-100x slower
in debug builds):
```bash
cmake -S target/mcp_http -B build_mcp_http -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/deps
cmake --build build_mcp_http -j8
INTERSPEC_MCP_BEARER_TOKEN=... INTERSPEC_LLM_ENDPOINT=... INTERSPEC_LLM_BEARER_TOKEN=... INTERSPEC_LLM_MODEL=... \
INTERSPEC_EXCHANGE_DIR=/tmp/exchange \
./build_mcp_http/InterSpecMcpServer --llm-config-dir /path/to/config --docroot build_mcp_http/LibInterSpec \
            -c target/mcp_http/wt_config_mcp_http.xml --userdatadir /tmp/interspec_user --http-port 8078
```
It finds Chrome/Chromium in the usual places (or set `<BrowserPath>`); it listens on 127.0.0.1
unless given `--http-address`.  The build leaves out the spectra database, remote RID, batch tools,
the file-query tool, and maps, and does not force developer checks on.


