#!/usr/bin/env python3
"""Tests the InterSpec stdio MCP server (target/mcp_stdio) as an MCP client would use it.

Uses only the Python standard library; run with `python3 -I`.  Checks: the initialize handshake
and protocol-version negotiation, ping, tools/list (no `userSession` argument, no tools that need a
browser), loading a spectrum, a synchronous and an async tool (with progress notifications),
error replies, that stdout carries nothing but JSON-RPC messages, and that the server exits
promptly when stdin is closed.  With --repeat N, also loads and fits N times, reporting the
server's memory use.

Example:
  python3 -I target/mcp_stdio/tests/stdio_test.py --exe target/mcp_stdio/build/InterSpec-mcp-stdio \\
      --datadir data --spectrum example_spectra/ba133_source_640s_20100317.n42 --energy 356
"""

import argparse
import json
import os
import queue
import subprocess
import sys
import threading
import time


class StdioClient:
    def __init__(self, cmd):
        self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
        self.next_id = 1
        self.messages = queue.Queue()
        self.bad_lines = []
        self.stderr_lines = []
        threading.Thread(target=self._read_stdout, daemon=True).start()
        threading.Thread(target=self._read_stderr, daemon=True).start()

    def _read_stdout(self):
        for raw in self.proc.stdout:
            line = raw.decode('utf-8', errors='replace').rstrip('\n')
            try:
                msg = json.loads(line)
                if not isinstance(msg, dict) or msg.get('jsonrpc') != '2.0':
                    raise ValueError('not a JSON-RPC 2.0 message')
                self.messages.put(msg)
            except ValueError:
                self.bad_lines.append(line)

    def _read_stderr(self):
        for raw in self.proc.stderr:
            self.stderr_lines.append(raw.decode('utf-8', errors='replace').rstrip('\n'))

    def send(self, msg):
        self.proc.stdin.write((json.dumps(msg) + '\n').encode())
        self.proc.stdin.flush()

    def send_raw(self, text):
        self.proc.stdin.write(text.encode())
        self.proc.stdin.flush()

    def request(self, method, params=None, timeout=300, notifications=None):
        """Sends a request, and returns its reply; notifications received meanwhile are appended
        to `notifications` (if given)."""
        msg_id = self.next_id
        self.next_id += 1
        msg = {'jsonrpc': '2.0', 'id': msg_id, 'method': method}
        if params is not None:
            msg['params'] = params
        self.send(msg)
        deadline = time.time() + timeout
        while True:
            remaining = deadline - time.time()
            if remaining <= 0:
                raise TimeoutError('no reply to %s within %d s' % (method, timeout))
            msg = self.messages.get(timeout=remaining)
            if msg.get('id') == msg_id and 'method' not in msg:
                return msg
            if notifications is not None:
                notifications.append(msg)

    def call(self, tool, arguments=None, progress_token=None, notifications=None, timeout=300):
        params = {'name': tool, 'arguments': arguments or {}}
        if progress_token is not None:
            params['_meta'] = {'progressToken': progress_token}
        return self.request('tools/call', params, timeout, notifications)

    def rss_mb(self):
        try:
            out = subprocess.run(['ps', '-o', 'rss=', '-p', str(self.proc.pid)],
                                 capture_output=True, text=True).stdout.strip()
            return int(out) / 1024.0
        except (ValueError, OSError):
            return -1.0


def png_size(data):
    """(width, height) from a PNG's IHDR chunk, or None if not a PNG."""
    if len(data) < 24 or data[:8] != b'\x89PNG\r\n\x1a\n' or data[12:16] != b'IHDR':
        return None
    return int.from_bytes(data[16:20], 'big'), int.from_bytes(data[20:24], 'big')


def image_of(reply):
    """(bytes, mimeType, description) of a tools/call reply's image block."""
    import base64
    if 'error' in reply:
        raise RuntimeError('tool error: %s' % reply['error'].get('message'))
    content = reply['result']['content']
    text = json.loads(content[0]['text'])
    images = [c for c in content if c.get('type') == 'image']
    if not images:
        raise RuntimeError('no image in the reply: %s' % json.dumps(text)[:200])
    return base64.b64decode(images[0]['data']), images[0]['mimeType'], text.get('description', '')


def tool_result(reply):
    """The tool's JSON result, from the text block of a tools/call reply."""
    if 'error' in reply:
        raise RuntimeError('tool error: %s' % reply['error'].get('message'))
    content = reply['result']['content']
    return json.loads(content[0]['text'])


class Checker:
    def __init__(self):
        self.failures = 0

    def check(self, name, ok, detail=''):
        print('%s: %s%s' % ('PASS' if ok else 'FAIL', name, (' - ' + detail) if detail else ''))
        sys.stdout.flush()
        if not ok:
            self.failures += 1
        return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exe', required=True)
    parser.add_argument('--datadir', required=True)
    parser.add_argument('--userdatadir', help='Defaults to an empty temporary directory.')
    parser.add_argument('--spectrum', required=True)
    parser.add_argument('--energy', type=float, required=True, help='Energy (keV) of a peak to fit.')
    parser.add_argument('--repeat', type=int, default=0,
                        help='Also load, fit and delete a peak this many times, reporting memory use.')
    parser.add_argument('--alt-spectrum', help='With --repeat, a second file to load between the fits.')
    parser.add_argument('--background', help='A background spectrum file, for the image checks.')
    parser.add_argument('--time-spectrum', help='A search-mode/passthrough file, for the time-chart image.')
    parser.add_argument('--save-images', help='Write the images rendered to this directory.')
    args = parser.parse_args()

    user_data_dir = args.userdatadir
    if not user_data_dir:
        import tempfile
        user_data_dir = tempfile.mkdtemp(prefix='interspec_stdio_test_')

    spectrum = os.path.abspath(args.spectrum)
    cmd = [args.exe, '--datadir', args.datadir, '--userdatadir', user_data_dir]
    client = StdioClient(cmd)
    c = Checker()

    try:
        # Handshake: a supported version is echoed; an unknown one gets our latest.
        reply = client.request('initialize', {'protocolVersion': '2025-06-18', 'capabilities': {},
                                              'clientInfo': {'name': 'stdio_test', 'version': '1'}}, timeout=120)
        c.check('initialize', reply.get('result', {}).get('protocolVersion') == '2025-06-18',
                json.dumps(reply.get('result', {}).get('serverInfo')))
        reply = client.request('initialize', {'protocolVersion': '1999-01-01'})
        c.check('initialize with an unknown version', reply['result']['protocolVersion'] == '2025-11-25')
        client.send({'jsonrpc': '2.0', 'method': 'notifications/initialized'})

        reply = client.request('ping')
        c.check('ping', reply.get('result') == {})

        reply = client.request('tools/list')
        tools = {t['name']: t for t in reply['result']['tools']}
        c.check('tools/list', len(tools) > 20, '%d tools' % len(tools))
        with_session = [n for n, t in tools.items() if 'userSession' in t['inputSchema'].get('properties', {})]
        c.check('no userSession argument', not with_session, ', '.join(with_session))
        # The default llm_config.xml does not let MCP clients use the assistant.
        assistant_tools = [n for n in ('assistant_submit_prompt', 'invoke_agent') if n in tools]
        c.check('no assistant tools with the default config', not assistant_tools, ', '.join(assistant_tools))
        c.check('image tools offered', 'get_spectrum_image' in tools and 'get_time_chart_image' in tools)

        reply = client.call('load_spectrum_file', {'path': spectrum, 'specType': 'Foreground'})
        result = tool_result(reply)
        c.check('load_spectrum_file', result.get('success', True) is not False and 'userSession' not in result,
                json.dumps(result)[:200])

        reply = client.call('get_spectrum_info', {})
        result = tool_result(reply)
        c.check('get_spectrum_info', 'foreground' in json.dumps(result).lower(), json.dumps(result)[:160])

        notifications = []
        reply = client.call('add_analysis_peak', {'energy': args.energy}, progress_token='fit-1',
                            notifications=notifications)
        result = tool_result(reply)
        c.check('add_analysis_peak (async)', result.get('status') == 'success'
                and abs(result.get('fitPeakEnergy', 0) - args.energy) < 5, json.dumps(result)[:200])
        c.check('no _visuals in results', '_visuals' not in result)

        # Chart images, drawn by the server.
        def save(name, data):
            if args.save_images:
                os.makedirs(args.save_images, exist_ok=True)
                with open(os.path.join(args.save_images, name), 'wb') as f:
                    f.write(data)

        if args.background:
            tool_result(client.call('load_spectrum_file', {'path': os.path.abspath(args.background),
                                                           'specType': 'Background'}))
        data, mime, desc = image_of(client.call('get_spectrum_image', {}))
        c.check('get_spectrum_image (png)', mime == 'image/png' and png_size(data) == (1200, 600),
                '%s %s; %s' % (mime, png_size(data), desc))
        save('spectrum.png', data)

        data, mime, desc = image_of(client.call('get_spectrum_image', {
            'energy_range': {'min_keV': args.energy - 60, 'max_keV': args.energy + 60},
            'y_axis_scale': 'linear', 'max_size': 800}))
        c.check('get_spectrum_image with energy_range, linear, max_size', png_size(data) == (800, 400),
                '%s; %s' % (png_size(data), desc))
        save('spectrum_zoom_linear.png', data)

        if args.background:
            data, mime, desc = image_of(client.call('get_spectrum_image', {'background_subtract': True}))
            c.check('get_spectrum_image background_subtract', png_size(data) is not None, desc)
            save('spectrum_bkg_subtracted.png', data)

        data, mime, desc = image_of(client.call('get_spectrum_image', {'format': 'jpeg', 'max_size': 640}))
        c.check('get_spectrum_image (jpeg)', mime == 'image/jpeg' and data[:3] == b'\xff\xd8\xff', desc)
        save('spectrum.jpg', data)

        data, mime, desc = image_of(client.call('get_spectrum_image', {'format': 'svg'}))
        c.check('get_spectrum_image (svg)', mime == 'image/svg+xml' and b'<svg' in data[:400], desc)
        save('spectrum.svg', data)

        if args.time_spectrum:
            tool_result(client.call('load_spectrum_file', {'path': os.path.abspath(args.time_spectrum),
                                                           'specType': 'Foreground'}))
            data, mime, desc = image_of(client.call('get_time_chart_image', {}))
            c.check('get_time_chart_image', png_size(data) == (1200, 600), '%s; %s' % (png_size(data), desc))
            save('time_chart.png', data)
            data, mime, desc = image_of(client.call('get_spectrum_image', {}))
            save('time_spectrum.png', data)
            tool_result(client.call('load_spectrum_file', {'path': spectrum, 'specType': 'Foreground'}))

        reply = client.call('no_such_tool', {})
        c.check('unknown tool is an error', reply.get('error', {}).get('code') == -32602)

        reply = client.request('no/such/method')
        c.check('unknown method is an error', reply.get('error', {}).get('code') == -32601)

        client.send_raw('{this is not json\n')
        msg = client.messages.get(timeout=30)
        c.check('malformed JSON is a parse error', msg.get('error', {}).get('code') == -32700)

        rss = []
        for i in range(args.repeat):
            if args.alt_spectrum:
                tool_result(client.call('load_spectrum_file', {'path': os.path.abspath(args.alt_spectrum),
                                                               'specType': 'Foreground'}))
            tool_result(client.call('load_spectrum_file', {'path': spectrum, 'specType': 'Foreground'}))
            client.call('edit_analysis_peak', {'energy': args.energy, 'action': 'DeletePeak'})
            tool_result(client.call('add_analysis_peak', {'energy': args.energy}))
            rss.append(client.rss_mb())
            if (i + 1) % 10 == 0:
                print('  after %d loads+fits: RSS %.0f MB' % (i + 1, rss[-1]))
        if rss:
            c.check('repeated loads and fits', True, 'RSS %.0f MB after the first, %.0f MB after %d'
                    % (rss[0], rss[-1], len(rss)))
    except Exception as e:
        c.check('no exceptions', False, repr(e))
    finally:
        # Closing stdin asks the server to exit.
        start = time.time()
        try:
            client.proc.stdin.close()
            code = client.proc.wait(timeout=60)
        except subprocess.TimeoutExpired:
            client.proc.kill()
            code = None
        c.check('exits when stdin closes', code == 0, 'exit code %s after %.1f s' % (code, time.time() - start))
        time.sleep(0.2)
        c.check('stdout carried only JSON-RPC', not client.bad_lines, repr(client.bad_lines[:3]))

    if c.failures:
        print('\n--- server stderr (last 30 lines) ---')
        print('\n'.join(client.stderr_lines[-30:]))
    print('\n%s' % ('ALL PASSED' if not c.failures else '%d FAILED' % c.failures))
    return 1 if c.failures else 0


if __name__ == '__main__':
    sys.exit(main())
