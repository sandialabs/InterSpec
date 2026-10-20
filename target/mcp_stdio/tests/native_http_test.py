#!/usr/bin/env python3
"""Tests InterSpec's native HTTP client (NativeHttp) against local servers.

Starts an HTTPS server with a self-signed certificate (made with the `openssl` command), an HTTP
CONNECT proxy, and a proxy that always asks for authentication (407); then runs the
test_native_http program (built from tests/test_native_http.cpp, for whichever backend the build
has) for each case, and checks the results.  Standard-library Python only; run with `python3 -I`.

Example:
  python3 -I target/mcp_stdio/tests/native_http_test.py --exe target/mcp_stdio/build/test_native_http
"""

import argparse
import http.server
import json
import os
import select
import socket
import socketserver
import ssl
import subprocess
import sys
import tempfile
import threading
import time

# A host name only the test proxy can reach (it maps it to 127.0.0.1).  Not localhost, which
#  NSURLSession never sends through a proxy; and as it can not be resolved, a request that does
#  not go through the proxy fails.
PROXIED_HOST = 'proxy-test.invalid'


class TestHandler(http.server.BaseHTTPRequestHandler):
    protocol_version = 'HTTP/1.1'

    def log_message(self, *args):
        pass

    def send_body(self, status, body, extra_headers=()):
        self.send_response(status)
        for name, value in extra_headers:
            self.send_header(name, value)
        self.send_header('Content-Type', 'text/plain')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == '/hello':
            self.send_body(200, b'hello world', [('X-Test', 'yes')])
        elif self.path == '/stream':
            # Chunked, with pauses, so it arrives in pieces.  (Typed, and over 512 bytes a piece, as
            #  NSURLSession holds back the start of an untyped body to sniff its type.)
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Transfer-Encoding', 'chunked')
            self.end_headers()
            for i in range(5):
                piece = ('part%d;' % i).encode() + b' ' * 1000
                self.wfile.write(b'%x\r\n%s\r\n' % (len(piece), piece))
                self.wfile.flush()
                time.sleep(0.4)
            self.wfile.write(b'0\r\n\r\n')
        elif self.path == '/slow':
            time.sleep(6)
            self.send_body(200, b'late')
        elif self.path == '/trickle':
            # A byte every half second, forever-ish: never idle, but never finished.
            self.send_response(200)
            self.send_header('Content-Length', '100000')
            self.end_headers()
            try:
                for _ in range(200):
                    self.wfile.write(b'x')
                    self.wfile.flush()
                    time.sleep(0.5)
            except OSError:
                pass
        elif self.path == '/429':
            self.send_body(429, b'{"error":"rate limited"}', [('Retry-After', '7')])
        elif self.path == '/big':
            self.send_body(200, b'y' * 50000)
        elif self.path == '/redirect':
            self.send_response(302)
            self.send_header('Location', '/hello')
            self.send_header('Content-Length', '0')
            self.end_headers()
        else:
            self.send_body(404, b'not found')

    def do_POST(self):
        length = int(self.headers.get('Content-Length', 0))
        data = self.rfile.read(length)
        self.send_body(200, b'echo:' + data + b';auth=' + (self.headers.get('Authorization') or '').encode())


class ThreadingServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True
    allow_reuse_address = True

    def handle_error(self, request, client_address):
        pass   # clients that hang up mid-response (the cancel tests, untrusted certificates) are expected


class ProxyHandler(socketserver.BaseRequestHandler):
    """A CONNECT proxy; with `require_auth`, it answers every request with 407."""
    require_auth = False
    connects = []

    def handle(self):
        data = b''
        while b'\r\n\r\n' not in data:
            more = self.request.recv(4096)
            if not more:
                return
            data += more
        request_line = data.split(b'\r\n')[0].decode()
        if self.require_auth:
            self.request.sendall(b'HTTP/1.1 407 Proxy Authentication Required\r\n'
                                 b'Proxy-Authenticate: Basic realm="test"\r\nContent-Length: 0\r\n\r\n')
            return
        method, target, _ = request_line.split(' ', 2)
        if method != 'CONNECT':
            self.request.sendall(b'HTTP/1.1 405 Method Not Allowed\r\nContent-Length: 0\r\n\r\n')
            return
        ProxyHandler.connects.append(target)
        host, port = target.rsplit(':', 1)
        if host.endswith('.invalid'):
            host = '127.0.0.1'   # the test's stand-in for a host only the proxy can reach
        upstream = socket.create_connection((host, int(port)))
        self.request.sendall(b'HTTP/1.1 200 Connection established\r\n\r\n')
        sockets = [self.request, upstream]
        try:
            while True:
                readable, _, _ = select.select(sockets, [], [], 30)
                if not readable:
                    break
                for s in readable:
                    chunk = s.recv(65536)
                    if not chunk:
                        return
                    (upstream if s is self.request else self.request).sendall(chunk)
        finally:
            upstream.close()


class AuthProxyHandler(ProxyHandler):
    require_auth = True


def start(server):
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server.server_address[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exe', required=True, help='The test_native_http program.')
    args = parser.parse_args()

    workdir = tempfile.mkdtemp(prefix='native_http_test_')
    cert, key = os.path.join(workdir, 'cert.pem'), os.path.join(workdir, 'key.pem')
    subprocess.run(['openssl', 'req', '-x509', '-newkey', 'rsa:2048', '-nodes', '-keyout', key, '-out', cert,
                    '-days', '30', '-subj', '/CN=localhost',
                    '-addext', 'subjectAltName=DNS:localhost,DNS:%s,IP:127.0.0.1' % PROXIED_HOST,
                    '-addext', 'extendedKeyUsage=serverAuth',
                    '-addext', 'basicConstraints=critical,CA:TRUE'],
                   check=True, capture_output=True)

    https = ThreadingServer(('127.0.0.1', 0), TestHandler)
    context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    context.load_cert_chain(cert, key)
    https.socket = context.wrap_socket(https.socket, server_side=True)
    https_port = start(https)
    http_port = start(ThreadingServer(('127.0.0.1', 0), TestHandler))
    proxy_port = start(socketserver.ThreadingTCPServer(('127.0.0.1', 0), ProxyHandler))
    auth_proxy_port = start(socketserver.ThreadingTCPServer(('127.0.0.1', 0), AuthProxyHandler))

    # A port with nothing listening.
    s = socket.socket()
    s.bind(('127.0.0.1', 0))
    closed_port = s.getsockname()[1]
    s.close()

    base = 'https://localhost:%d' % https_port
    failures = 0

    def run(name, check, *exe_args, env=None):
        nonlocal failures
        full_env = dict(os.environ)
        for var in ('https_proxy', 'HTTPS_PROXY', 'http_proxy', 'all_proxy', 'ALL_PROXY', 'no_proxy', 'NO_PROXY'):
            full_env.pop(var, None)
        full_env.update(env or {})
        proc = subprocess.run([args.exe] + list(exe_args), capture_output=True, text=True, timeout=120, env=full_env)
        try:
            result = json.loads(proc.stdout.strip().splitlines()[-1])
        except (ValueError, IndexError):
            result = {'stdout': proc.stdout, 'stderr': proc.stderr[-500:]}
        ok = False
        try:
            ok = bool(check(result))
        except Exception as e:
            result['checkException'] = repr(e)
        print('%s: %s%s' % ('PASS' if ok else 'FAIL', name, '' if ok else ' - ' + json.dumps(result)[:600]))
        sys.stdout.flush()
        if not ok:
            failures += 1
        return result

    first = run('HTTPS with the test CA', lambda r: r['error'] == 'Ok' and r['status'] == 200
                and r['body'] == 'hello world' and r['headers'].get('X-Test') == 'yes',
                '--url', base + '/hello', '--ca', cert)
    backend = first.get('backend', '?')
    print('backend:', backend)
    is_curl = backend.startswith('libcurl')

    run('HTTPS without the test CA is untrusted', lambda r: r['error'] == 'TlsCertificateUntrusted',
        '--url', base + '/hello')
    run('disableCertVerification', lambda r: r['error'] == 'Ok' and r['body'] == 'hello world',
        '--url', base + '/hello', '--insecure')
    run('POST body and headers', lambda r: r['error'] == 'Ok' and r['body'] == 'echo:{"q":1};auth=',
        '--url', base + '/echo', '--method', 'POST', '--body', '{"q":1}', '--ca', cert)
    run('plain http', lambda r: r['error'] == 'Ok' and r['body'] == 'hello world',
        '--url', 'http://127.0.0.1:%d/hello' % http_port)
    run('streamed chunked response', lambda r: r['error'] == 'Ok' and r['bodyBytes'] == 5 * 1006
        and r['body'].startswith('part0;') and r['chunks'] >= 3, '--url', base + '/stream', '--ca', cert)
    run('cancel mid-stream', lambda r: r['error'] == 'Cancelled' and r['elapsedMs'] < 1800,
        '--url', base + '/stream', '--ca', cert, '--cancel-after-bytes', '1')
    run('handler aborts mid-stream', lambda r: r['error'] == 'AbortedByHandler',
        '--url', base + '/stream', '--ca', cert, '--abort-after-bytes', '1')
    # (The Apple backend reports its timeout - which is an idle one - as Timeout.)
    run('idle timeout', lambda r: r['error'] in (('IdleTimeout',) if is_curl else ('IdleTimeout', 'Timeout'))
        and r['elapsedMs'] < 5500, '--url', base + '/slow', '--ca', cert, '--idle', '2')
    if is_curl:
        # (NSURLSession's resource timeout is not a whole-request limit, so this is libcurl-specific.)
        run('total timeout', lambda r: r['error'] == 'Timeout' and r['elapsedMs'] < 6000,
            '--url', base + '/trickle', '--ca', cert, '--total', '3', '--idle', '60')
    run('429 is a response, not an error', lambda r: r['error'] == 'Ok' and r['status'] == 429
        and r['headers'].get('Retry-After') == '7' and 'rate limited' in r['body'], '--url', base + '/429', '--ca', cert)
    run('response too large', lambda r: r['error'] == 'ResponseTooLarge',
        '--url', base + '/big', '--ca', cert, '--max-size', '1000')
    run('redirect not followed by default', lambda r: r['error'] == 'Ok' and r['status'] == 302,
        '--url', base + '/redirect', '--ca', cert)
    run('redirect followed when asked', lambda r: r['error'] == 'Ok' and r['status'] == 200 and r['body'] == 'hello world',
        '--url', base + '/redirect', '--ca', cert, '--follow')
    run('host not found', lambda r: r['error'] == 'HostNotFound', '--url', 'https://no-such-host.invalid/hello')
    run('connection refused', lambda r: r['error'] == 'ConnectFailed', '--url', 'https://127.0.0.1:%d/' % closed_port)

    proxied = 'https://%s:%d' % (PROXIED_HOST, https_port)
    ProxyHandler.connects.clear()
    run('through a CONNECT proxy (proxyUrl)', lambda r: r['error'] == 'Ok' and r['body'] == 'hello world'
        and ProxyHandler.connects, '--url', proxied + '/hello', '--ca', cert, '--proxy', 'http://127.0.0.1:%d' % proxy_port)
    if is_curl:
        ProxyHandler.connects.clear()
        run('through a CONNECT proxy (https_proxy)', lambda r: r['error'] == 'Ok' and r['body'] == 'hello world'
            and ProxyHandler.connects, '--url', proxied + '/hello', '--ca', cert,
            env={'https_proxy': 'http://127.0.0.1:%d' % proxy_port})
    run('proxy wants authentication (407)', lambda r: r['error'] == 'ProxyAuthRequired',
        '--url', proxied + '/hello', '--ca', cert, '--proxy', 'http://127.0.0.1:%d' % auth_proxy_port)
    run('proxy unreachable', lambda r: r['error'] in ('ProxyUnreachable', 'ConnectFailed'),
        '--url', proxied + '/hello', '--ca', cert, '--proxy', 'http://127.0.0.1:%d' % closed_port)

    print('\n%s' % ('ALL PASSED' if not failures else '%d FAILED' % failures))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
