#!/usr/bin/env python3
"""Torture test of the InterSpec stdio MCP server (target/mcp_stdio).

Where stdio_test.py checks the main workflow, this tries to break things:

  sweep  - calls every MCP tool with plausible arguments, using the call tables of the HTTP
           server's torture test (target/mcp_http/tests/torture_test.py), against a Ba-133 spectrum
           with background, a uranium HPGe spectrum, and a passthrough file; plus calls that
           should fail.  Each call must return a result or a clean error within its deadline -
           a tool that waits on a browser would hang here - and the server must still respond
           afterwards.  Images are saved to --out.
  chaos  - many calls in flight at once, out-of-order replies, oversized, malformed, split and
           batched messages, invalid UTF-8, a client that stops reading stdout for a while,
           stdin closing during a long calculation, SIGTERM, and memory over many image calls.
  llm    - (needs --llm-config-dir with an LLM, and <AllowAssistantPrompt>) heavy assistant turns,
           each in its own server process, in parallel: nuclide ID, U enrichment, an activity fit,
           and a cancelled turn.

Standard-library Python only; run with `python3 -I`.

Example:
  python3 -I target/mcp_stdio/tests/stdio_torture.py --exe target/mcp_stdio/build/InterSpec-mcp-stdio \\
      --repo . --out /tmp/stdio_torture [--phases sweep,chaos,llm --llm-config-dir <dir>]
"""

import argparse
import importlib.util
import json
import os
import queue
import signal
import subprocess
import sys
import tempfile
import threading
import time


def load_http_torture(repo):
    """The HTTP server's torture test, for its call tables, report and helpers."""
    path = os.path.join(repo, 'target', 'mcp_http', 'tests', 'torture_test.py')
    sys.dont_write_bytecode = True   # no __pycache__ left in the source tree
    spec = importlib.util.spec_from_file_location('http_torture', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class StdioServer:
    """A server process, with replies routed to their requests by id, so many calls can be in flight
    from several threads.  Mirrors the HTTP torture client's call() and rpc()."""

    def __init__(self, exe, datadir, userdatadir, llm_config_dir=None, extra_env=None):
        cmd = [exe, '--datadir', datadir, '--userdatadir', userdatadir]
        if llm_config_dir:
            cmd += ['--llm-config-dir', llm_config_dir]
        env = dict(os.environ)
        env.update(extra_env or {})
        self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        self.lock = threading.Lock()
        self.write_lock = threading.Lock()
        self.next_id = 1
        self.waiting = {}         # id -> queue for its reply
        self.notifications = []
        self.unmatched = []       # replies with no request waiting (e.g. null-id errors)
        self.bad_lines = []       # anything on stdout that is not a JSON-RPC message
        self.stderr_lines = []
        self.reading = threading.Event()
        self.reading.set()
        threading.Thread(target=self._read_stdout, daemon=True).start()
        threading.Thread(target=self._read_stderr, daemon=True).start()
        self.rpc('initialize', {'protocolVersion': '2025-06-18', 'capabilities': {},
                                'clientInfo': {'name': 'stdio_torture', 'version': '1'}}, timeout=120)
        self.notify('notifications/initialized')

    def _read_stdout(self):
        for raw in self.proc.stdout:
            self.reading.wait()  # (cleared to simulate a client that stops reading)
            line = raw.decode('utf-8', errors='replace').rstrip('\n')
            try:
                msg = json.loads(line)
                if not isinstance(msg, dict) or msg.get('jsonrpc') != '2.0':
                    raise ValueError('not JSON-RPC')
            except ValueError:
                self.bad_lines.append(line[:300])
                continue
            if 'method' in msg:
                with self.lock:
                    self.notifications.append(msg)
                continue
            with self.lock:
                q = self.waiting.pop(msg.get('id'), None)
                if q is None:
                    self.unmatched.append(msg)
            if q is not None:
                q.put(msg)

    def _read_stderr(self):
        for raw in self.proc.stderr:
            self.stderr_lines.append(raw.decode('utf-8', errors='replace').rstrip('\n'))

    def write(self, data):
        with self.write_lock:
            self.proc.stdin.write(data)
            self.proc.stdin.flush()

    def send(self, msg):
        self.write((json.dumps(msg) + '\n').encode())

    def notify(self, method, params=None):
        msg = {'jsonrpc': '2.0', 'method': method}
        if params is not None:
            msg['params'] = params
        self.send(msg)

    def start(self, method, params=None):
        """Sends a request; returns (id, queue its reply will arrive on)."""
        q = queue.Queue()
        with self.lock:
            msg_id = self.next_id
            self.next_id += 1
            self.waiting[msg_id] = q
        msg = {'jsonrpc': '2.0', 'id': msg_id, 'method': method}
        if params is not None:
            msg['params'] = params
        self.send(msg)
        return msg_id, q

    def rpc(self, method, params=None, timeout=600):
        _, q = self.start(method, params)
        return q.get(timeout=timeout)

    def call(self, tool, args=None, timeout=600):
        """As the HTTP torture client: (result or None, [(mime, bytes)], error or None, seconds)."""
        import base64
        t0 = time.time()
        try:
            reply = self.rpc('tools/call', {'name': tool, 'arguments': args or {}}, timeout=timeout)
        except queue.Empty:
            return None, [], 'transport: no reply within %d s' % timeout, time.time() - t0
        except Exception as e:
            return None, [], 'transport: %r' % (e,), time.time() - t0
        elapsed = time.time() - t0
        if 'error' in reply:
            return None, [], reply['error'].get('message', str(reply['error'])), elapsed
        content = reply['result']['content']
        if reply['result'].get('isError'):
            return None, [], content[0]['text'], elapsed
        result = json.loads(content[0]['text'])
        images = [(c['mimeType'], base64.b64decode(c['data'])) for c in content if c.get('type') == 'image']
        return result, images, None, elapsed

    def rss_mb(self):
        try:
            out = subprocess.run(['ps', '-o', 'rss=', '-p', str(self.proc.pid)], capture_output=True, text=True).stdout
            return int(out.strip()) / 1024.0
        except (ValueError, OSError):
            return -1.0

    def close(self, timeout=120):
        """Closes stdin; returns (exit code or None if it had to be killed, seconds)."""
        t0 = time.time()
        try:
            self.proc.stdin.close()
        except OSError:
            pass
        try:
            code = self.proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            self.proc.kill()
            code = None
        return code, time.time() - t0


class ResponsiveServer(StdioServer):
    """Checks the server still answers after every swept call (the HTTP sweep only does this for
    agent sessions)."""
    report = None

    def call(self, tool, args=None, timeout=600):
        out = super().call(tool, args, timeout)
        if tool not in ('get_loaded_spectra', 'tools/list'):
            _, _, err, secs = super().call('get_loaded_spectra', {}, timeout=60)
            if err:
                self.report.add('sweep', tool + ' [then responsive?]', 'FAIL', secs, err)
        return out


# ------------------------------------------------------------------------------------------------

def sweep(T, report, args):
    with tempfile.TemporaryDirectory(prefix='stdio_torture_ud_') as ud:
        server = ResponsiveServer(args.exe, os.path.join(args.repo, 'data'), ud)
        server.report = report
        sweep_args = argparse.Namespace(repo=args.repo)
        T.sweep(server, report, sweep_args, args.out, False)
        code, secs = server.close()
        report.add('sweep', 'exit after sweep', 'ok' if code == 0 else 'FAIL', secs, 'exit code %s' % code)
        report.add('sweep', 'stdout carried only JSON-RPC', 'FAIL' if server.bad_lines else 'ok', 0,
                   json.dumps(server.bad_lines[:3]))


def chaos(T, report, args):
    phase = 'chaos'
    print('\n== Chaos', flush=True)
    data = os.path.join(args.repo, 'data')
    ex = os.path.join(args.repo, 'example_spectra')
    fg = os.path.join(ex, 'ba133_source_640s_20100317.n42')
    bg = os.path.join(ex, 'background_20100317.n42')
    upath = os.path.join(args.repo, 'target', 'testing', 'test_data', 'uranium_40pct_selfatten_HPGe_15cm.n42')

    def new_server(ud):
        s = StdioServer(args.exe, data, ud)
        s.call('load_spectrum_file', {'path': fg, 'specType': 'Foreground'})
        s.call('load_spectrum_file', {'path': bg, 'specType': 'Background'})
        return s

    with tempfile.TemporaryDirectory(prefix='stdio_torture_ud_') as ud:
        s = new_server(ud)

        # Many calls in flight at once: synchronous, async fits, images, nuclear data.
        t0 = time.time()
        starts = []
        for i in range(8):
            starts.append(('fit %d' % i, s.start('tools/call', {'name': 'add_analysis_peak',
                                                                'arguments': {'energy': [81, 276.4, 302.9, 356, 383.8, 160.6, 223.2, 53.2][i],
                                                                              'DoNotAddToAnalysisPeaks': True}})))
            starts.append(('info %d' % i, s.start('tools/call', {'name': 'get_spectrum_info', 'arguments': {}})))
            starts.append(('image %d' % i, s.start('tools/call', {'name': 'get_spectrum_image', 'arguments': {'max_size': 300}})))
            starts.append(('decay %d' % i, s.start('tools/call', {'name': 'decay_calculator',
                                                                  'arguments': {'nuclide': 'Co60', 'activity': '1 uCi', 'time_duration': '%d years' % (i + 1)}})))
        missing, errors = [], []
        for label, (msg_id, q) in starts:
            try:
                reply = q.get(timeout=300)
                if reply.get('id') != msg_id:
                    errors.append('%s: wrong id' % label)
                elif 'error' in reply:
                    errors.append('%s: %s' % (label, reply['error'].get('message')))
                elif reply.get('result', {}).get('isError'):
                    errors.append('%s: %s' % (label, reply['result']['content'][0]['text']))
            except queue.Empty:
                missing.append(label)
        report.add(phase, '32 calls in flight', 'FAIL' if (missing or errors) else 'ok', time.time() - t0,
                   'missing: %s; errors: %s' % (missing, errors[:3]) if (missing or errors) else 'all answered, ids matched')

        # Out of order: a slow calculation first, then a fast one, which must not wait for it.
        slow_id, slow_q = s.start('tools/call', {'name': 'add_analysis_peaks_for_source',
                                                 'arguments': {'source': None, 'options': ['FitNormPeaks']}})
        time.sleep(0.2)
        t0 = time.time()
        fast_id, fast_q = s.start('tools/call', {'name': 'decay_calculator',
                                                 'arguments': {'nuclide': 'Cs137', 'activity': '1 uCi', 'time_duration': '1 years'}})
        fast_reply = fast_q.get(timeout=300)
        fast_secs = time.time() - t0
        slow_done_first = not slow_q.empty()
        slow_reply = slow_q.get(timeout=600)
        report.add(phase, 'replies out of order', 'info', fast_secs,
                   'fast call answered in %.1f s; slow one %s first; both answered: %s'
                   % (fast_secs, 'finished' if slow_done_first else 'had not finished',
                      'error' not in fast_reply and 'result' in slow_reply or 'error' in slow_reply))

        # Malformed input; the server must answer what it can, and keep going.
        def expect_unmatched_error(label, payload, code):
            before = len(s.unmatched)
            t1 = time.time()
            s.write(payload)
            while len(s.unmatched) == before and time.time() - t1 < 30:
                time.sleep(0.05)
            got = s.unmatched[before:] if len(s.unmatched) > before else []
            ok = bool(got) and got[0].get('error', {}).get('code') == code
            report.add(phase, label, 'ok' if ok else 'FAIL', time.time() - t1, json.dumps(got[:1])[:200])

        expect_unmatched_error('not JSON', b'{not json\n', -32700)
        expect_unmatched_error('invalid UTF-8', b'{"jsonrpc":"2.0","id":"\xff\xfe","method":"ping"}\n', -32700)
        expect_unmatched_error('batch (array)', b'[{"jsonrpc":"2.0","id":1,"method":"ping"}]\n', -32600)
        expect_unmatched_error('no method', b'{"jsonrpc":"2.0","id":null}\n', -32600)

        # A 5 MB argument.
        t0 = time.time()
        r, _, err, _ = s.call('source_info', {'source': 'x' * (5 * 1024 * 1024)}, timeout=120)
        report.add(phase, '5 MB argument', 'ok' if (err or r is not None) else 'FAIL', time.time() - t0, (err or json.dumps(r))[:120])

        # One message split over writes, and several messages in one write.
        msg_id, q = (None, None)
        with s.lock:
            msg_id = s.next_id
            s.next_id += 1
            q = queue.Queue()
            s.waiting[msg_id] = q
        line = json.dumps({'jsonrpc': '2.0', 'id': msg_id, 'method': 'ping'}).encode()
        s.write(line[:10])
        time.sleep(0.5)
        s.write(line[10:] + b'\n')
        try:
            r = q.get(timeout=30)
            report.add(phase, 'message split over writes', 'ok' if r.get('result') == {} else 'FAIL', 0, json.dumps(r))
        except queue.Empty:
            report.add(phase, 'message split over writes', 'FAIL', 30, 'no reply')

        ids = []
        with s.lock:
            for _ in range(5):
                ids.append((s.next_id, queue.Queue()))
                s.waiting[s.next_id] = ids[-1][1]
                s.next_id += 1
        s.write(b''.join(json.dumps({'jsonrpc': '2.0', 'id': i, 'method': 'ping'}).encode() + b'\n' for i, _ in ids)
                + b'\n   \n\r\n')
        got = 0
        for _, q in ids:
            try:
                q.get(timeout=30)
                got += 1
            except queue.Empty:
                pass
        report.add(phase, '5 messages (and blank lines) in one write', 'ok' if got == 5 else 'FAIL', 0, '%d of 5 answered' % got)

        # notifications/cancelled for a call in flight: ignored, but harmless.
        msg_id, q = s.start('tools/call', {'name': 'add_analysis_peak', 'arguments': {'energy': 356, 'DoNotAddToAnalysisPeaks': True}})
        s.notify('notifications/cancelled', {'requestId': msg_id, 'reason': 'torture'})
        try:
            r = q.get(timeout=120)
            report.add(phase, 'notifications/cancelled', 'info', 0,
                       'ignored; the call still answered (%s)' % ('error' if 'error' in r else 'result'))
        except queue.Empty:
            report.add(phase, 'notifications/cancelled', 'FAIL', 120, 'the call never answered')

        # A client that stops reading stdout for a while, with large replies pending.
        s.reading.clear()
        pending = [s.start('tools/call', {'name': 'get_spectrum_image', 'arguments': {'format': 'png'}}) for _ in range(20)]
        time.sleep(10)
        alive = (s.proc.poll() is None)
        s.reading.set()
        t0 = time.time()
        answered = 0
        for _, q in pending:
            try:
                q.get(timeout=120)
                answered += 1
            except queue.Empty:
                pass
        r, _, err, _ = s.call('get_loaded_spectra', {}, timeout=60)
        report.add(phase, 'client stops reading stdout for 10 s', 'ok' if (alive and answered == 20 and not err) else 'FAIL',
                   time.time() - t0, '%d of 20 large replies after reading resumed; then %s' % (answered, err or 'responsive'))

        # Memory over many image calls.
        rss0 = s.rss_mb()
        for i in range(60):
            s.call('get_spectrum_image', {'energy_range': {'min_keV': 20 + i, 'max_keV': 500 + 10 * i}})
        rss1 = s.rss_mb()
        report.add(phase, 'memory over 60 PNG images', 'ok' if rss1 - rss0 < 100 else 'FAIL', 0,
                   'RSS %.0f -> %.0f MB' % (rss0, rss1))

        code, secs = s.close()
        report.add(phase, 'exit', 'ok' if code == 0 else 'FAIL', secs, 'exit code %s' % code)
        report.add(phase, 'stdout carried only JSON-RPC', 'FAIL' if s.bad_lines else 'ok', 0, json.dumps(s.bad_lines[:3]))

    # stdin closing while a long calculation runs: the server must exit (after it finishes).
    with tempfile.TemporaryDirectory(prefix='stdio_torture_ud_') as ud:
        s = StdioServer(args.exe, data, ud)
        s.call('load_spectrum_file', {'path': upath, 'specType': 'Foreground'})
        presets = s.call('list_isotopics_presets', {})[0] or {}
        names = [p.get('name') if isinstance(p, dict) else p for p in presets.get('presets', [])]
        preset = next((n for n in names if n and 'U' in n and 'Pu' not in n), names[0] if names else None)
        s.call('load_isotopics_preset', {'preset': preset})
        s.start('tools/call', {'name': 'perform_isotopics_calculation', 'arguments': {}})
        time.sleep(1.0)
        code, secs = s.close(timeout=600)
        report.add(phase, 'stdin closed during isotopics', 'ok' if (code == 0 and secs < 15) else 'FAIL', secs,
                   'exit code %s after %.1f s' % (code, secs))

    # SIGTERM mid-call.
    with tempfile.TemporaryDirectory(prefix='stdio_torture_ud_') as ud:
        s = new_server(ud)
        s.start('tools/call', {'name': 'add_analysis_peaks_for_source', 'arguments': {'source': 'Ba133'}})
        time.sleep(0.3)
        t0 = time.time()
        s.proc.send_signal(signal.SIGTERM)
        try:
            code = s.proc.wait(timeout=30)
        except subprocess.TimeoutExpired:
            s.proc.kill()
            code = None
        report.add(phase, 'SIGTERM mid-call', 'ok' if code is not None else 'FAIL', time.time() - t0, 'exit code %s' % code)


def llm(T, report, args):
    phase = 'llm'
    print('\n== Heavy assistant turns, one process each, in parallel', flush=True)
    if not args.llm_config_dir:
        report.add(phase, 'llm', 'FAIL', 0, 'needs --llm-config-dir')
        return
    data = os.path.join(args.repo, 'data')
    ex = os.path.join(args.repo, 'example_spectra')
    files = {
        'fg': (os.path.join(ex, 'ba133_source_640s_20100317.n42'), os.path.join(ex, 'background_20100317.n42')),
        'u': (os.path.join(args.repo, 'target', 'testing', 'test_data', 'uranium_40pct_selfatten_HPGe_15cm.n42'), None),
    }
    jobs = [
        ('nuclide ID', 'fg', 'Identify all radioactive sources in this spectrum, including any background (NORM). '
                             'Fit the peaks needed, and justify each identification briefly.'),
        ('U enrichment', 'u', 'This is an HPGe measurement of uranium. Determine the U-235 enrichment using the isotopics tools, '
                              'and give the result with its uncertainty.'),
        ('activity fit', 'fg', 'The source is 50 cm from the detector face, unshielded. Fit the Ba-133 peaks, then use the '
                               'activity fit to estimate the Ba-133 activity. Report the activity with its uncertainty.'),
    ]
    rss = {}

    def server_with(which, ud):
        s = StdioServer(args.exe, data, ud, llm_config_dir=args.llm_config_dir)
        fg_path, bg_path = files[which]
        s.call('load_spectrum_file', {'path': fg_path, 'specType': 'Foreground'})
        if bg_path:
            s.call('load_spectrum_file', {'path': bg_path, 'specType': 'Background'})
        return s

    def wait_turn(s, r, deadline):
        err = None
        while r.get('status') == 'running' and time.time() < deadline:
            r, _, err, _ = s.call('assistant_get_result', {'turn_id': r['turn_id'], 'wait_seconds': 45})
            if err:
                break
        return r, err

    def run(label, which, prompt):
        t0 = time.time()
        with tempfile.TemporaryDirectory(prefix='stdio_torture_ud_') as ud:
            s = server_with(which, ud)
            try:
                r, _, err, _ = s.call('assistant_submit_prompt', {'prompt': prompt, 'wait_seconds': 20})
                if err:
                    report.add(phase, label, 'FAIL', time.time() - t0, err)
                    return
                if r['status'] == 'running':
                    _, _, e1, _ = s.call('get_spectrum_info', {})
                    _, _, e2, _ = s.call('add_analysis_peak', {'energy': 356})
                    report.add(phase, label + ': busy rules', 'ok' if (not e1 and e2) else 'FAIL', 0,
                               'read-only: %s; state-changing: %s' % (e1 or 'ok', T.short(e2 or 'allowed!', 80)))
                r, err = wait_turn(s, r, time.time() + 1500)
                rss[label] = s.rss_mb()
                status = 'ok' if (not err and r.get('status') == 'done') else 'FAIL'
                report.add(phase, label, status, time.time() - t0,
                           err or ('%d tool calls; RSS %.0f MB; %s' % (len(r.get('tool_calls', [])), rss[label],
                                                                     T.short(r.get('response') or r.get('error'), 200))))
                with open(os.path.join(args.out, 'llm_%s.json' % label.replace(' ', '_')), 'w') as f:
                    json.dump(r, f, indent=1)
            finally:
                code, secs = s.close()
                if code != 0 or s.bad_lines:
                    report.add(phase, label + ': exit', 'FAIL', secs, 'exit %s; bad stdout %s' % (code, s.bad_lines[:2]))

    def run_cancel():
        label = 'cancel'
        t0 = time.time()
        with tempfile.TemporaryDirectory(prefix='stdio_torture_ud_') as ud:
            s = server_with('fg', ud)
            try:
                r, _, err, _ = s.call('assistant_submit_prompt', {'wait_seconds': 15, 'prompt':
                    'Identify every radioactive source in this spectrum. Fit all peaks you need, one source at a time.'})
                if err or r.get('status') != 'running':
                    report.add(phase, label, 'FAIL', time.time() - t0, err or 'turn finished before it could be cancelled')
                    return
                turn = r['turn_id']
                c, _, err, _ = s.call('assistant_cancel', {'turn_id': turn})
                report.add(phase, label + ': assistant_cancel', 'ok' if (not err and c.get('cancelled')) else 'FAIL', 0, err or json.dumps(c))
                r, _, err, secs = s.call('assistant_get_result', {'turn_id': turn, 'wait_seconds': 10})
                report.add(phase, label + ': result is cancelled', 'ok' if (not err and r.get('status') == 'cancelled') else 'FAIL',
                           secs, err or json.dumps(r)[:200])
                p1, _, _, _ = s.call('get_peaks', {})
                time.sleep(20)
                p2, _, _, _ = s.call('get_peaks', {})
                n1, n2 = len((p1 or {}).get('peaks', [])), len((p2 or {}).get('peaks', []))
                report.add(phase, label + ': turn stopped', 'ok' if n1 == n2 else 'FAIL', 20, 'peaks %d then %d, 20 s apart' % (n1, n2))
                r, _, err, _ = s.call('assistant_submit_prompt', {'wait_seconds': 45,
                    'prompt': 'In one short sentence, and without calling any tools: what is the half-life of Cs-137?'})
                if not err:
                    r, err = wait_turn(s, r, time.time() + 600)
                answer = (r or {}).get('response') or ''
                report.add(phase, label + ': next prompt answered', 'ok' if (not err and r.get('status') == 'done' and '30' in answer) else 'FAIL',
                           time.time() - t0, err or T.short(answer, 200))
            finally:
                s.close()

    threads = [threading.Thread(target=run, args=j) for j in jobs] + [threading.Thread(target=run_cancel)]
    [t.start() for t in threads]
    [t.join() for t in threads]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exe', required=True)
    parser.add_argument('--repo', default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
    parser.add_argument('--out', default='stdio_torture_out')
    parser.add_argument('--phases', default='sweep,chaos')
    parser.add_argument('--llm-config-dir', help='For the llm phase: a directory with an llm_config.xml.')
    args = parser.parse_args()
    args.exe = os.path.abspath(args.exe)
    args.repo = os.path.abspath(args.repo)

    os.makedirs(args.out, exist_ok=True)
    T = load_http_torture(args.repo)
    report = T.Report(args.out)
    t0 = time.time()
    phases = args.phases.split(',')
    if 'sweep' in phases:
        sweep(T, report, args)
    if 'chaos' in phases:
        chaos(T, report, args)
    if 'llm' in phases:
        llm(T, report, args)

    path = report.save()
    fails = report.failures()
    slow = [r for r in report.rows if r['status'] == 'SLOW']
    print('\n%d checks in %.0f s; %d failure(s), %d slow.  Report: %s' % (len(report.rows), time.time() - t0, len(fails), len(slow), path))
    for r in fails + slow:
        print('  %s %s/%s: %s' % (r['status'], r['phase'], r['name'], r['detail']))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
