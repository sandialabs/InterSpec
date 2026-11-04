#!/usr/bin/env python3
"""Smoke test of the InterSpec MCP analysis server (target/mcp_http).

Uses only the Python standard library.  Exercises: session creation, uploading by PUT and by
upload_file, loading spectra, images, exporting and downloading, resuming from an exported N42,
path confinement, authentication, and closing sessions (checking their files and URLs go away).

Example:
  INTERSPEC_MCP_BEARER_TOKEN=... ./smoke_test.py --url http://127.0.0.1:8078 \\
      --foreground example_spectra/ba133_source_640s_20100317.n42 \\
      --background example_spectra/background_20100317.n42 \\
      --exchange-dir /path/to/exchange [--prompt "What is in this spectrum?"]

--exchange-dir is optional: give it (the host side of the exchange directory) to also test
writing files directly into a session directory, and that directories are removed.
--prompt is optional: it needs a working LLM configured on the server.
"""

import argparse
import base64
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request


class McpClient:
    def __init__(self, url, token):
        self.url = url.rstrip('/')
        self.token = token
        self.next_id = 1

    def http(self, method, url, body=None, headers=None, auth=True):
        """Returns (status, headers, body bytes)."""
        hdrs = dict(headers or {})
        if auth and self.token:
            hdrs['Authorization'] = 'Bearer ' + self.token
        req = urllib.request.Request(url, data=body, method=method, headers=hdrs)
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return resp.status, dict(resp.headers), resp.read()
        except urllib.error.HTTPError as e:
            return e.code, dict(e.headers), e.read()

    def rpc(self, method, params=None, auth=True):
        msg = {'jsonrpc': '2.0', 'id': self.next_id, 'method': method}
        self.next_id += 1
        if params is not None:
            msg['params'] = params
        status, _, body = self.http('POST', self.url + '/mcp-api', json.dumps(msg).encode(),
                                    {'Content-Type': 'application/json',
                                     'Accept': 'application/json'}, auth)
        if status != 200:
            return {'http_status': status, 'body': body.decode(errors='replace')}
        return json.loads(body)

    def call(self, tool, **args):
        """Calls a tool; returns (result dict, error message)."""
        reply = self.rpc('tools/call', {'name': tool, 'arguments': args})
        if 'error' in reply:
            return None, reply['error'].get('message', str(reply['error']))
        if 'http_status' in reply:
            return None, 'HTTP %d: %s' % (reply['http_status'], reply['body'])
        text = reply['result']['content'][0]['text']
        if reply['result'].get('isError'):
            return None, text
        return json.loads(text), None


failures = []


def check(condition, what):
    print(('  ok   ' if condition else '  FAIL ') + what)
    if not condition:
        failures.append(what)
    return condition


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--url', default='http://127.0.0.1:8078')
    parser.add_argument('--token', default=os.environ.get('INTERSPEC_MCP_BEARER_TOKEN', ''))
    parser.add_argument('--foreground', required=True)
    parser.add_argument('--background')
    parser.add_argument('--exchange-dir')
    parser.add_argument('--prompt')
    args = parser.parse_args()

    client = McpClient(args.url, args.token)

    print('Server')
    status, _, body = client.http('GET', client.url + '/mcp-api/healthz', auth=False)
    check(status == 200 and json.loads(body)['status'] == 'ok', 'health check')
    init = client.rpc('initialize', {'protocolVersion': '2024-11-05', 'capabilities': {},
                                     'clientInfo': {'name': 'smoke', 'version': '1'}})
    check('instructions' in init.get('result', {}), 'initialize returns instructions')
    tools = {t['name'] for t in client.rpc('tools/list')['result']['tools']}
    check({'create_session', 'close_session', 'export_analysis', 'load_spectrum_file'} <= tools,
          'tools/list has the agent tools')
    check(McpClient(args.url, 'wrong-token').rpc('tools/list').get('http_status') == 401,
          'a wrong bearer token is refused')
    check(McpClient(args.url, '').rpc('tools/list').get('http_status') == 401,
          'no bearer token is refused')

    _, err = client.call('get_spectrum_info')
    check(err is not None and 'create_session' in err, 'session tools without userSession say to create_session')

    print('Session A')
    sess, err = client.call('create_session')
    if not check(err is None, 'create_session (%s)' % err):
        return 1
    handle = sess['userSession']
    upload_url = sess['upload_url']
    print('    session_dir:', sess['session_dir'])

    fg_name = os.path.basename(args.foreground)
    with open(args.foreground, 'rb') as f:
        status, _, body = client.http('PUT', upload_url + '/' + fg_name, f.read())
    check(status == 201, 'PUT upload of the foreground')
    fg_path = json.loads(body)['path'] if status == 201 else 'input/' + fg_name

    status, _, _ = McpClient(args.url, '').http('PUT', upload_url + '/x.n42', b'abc', auth=False)
    check(status == 401, 'upload without the bearer token is refused')

    result, err = client.call('load_spectrum_file', userSession=handle, path=fg_path, specType='Foreground')
    check(err is None and result.get('success'), 'load the foreground (%s)' % err)

    if args.background:
        with open(args.background, 'rb') as f:
            result, err = client.call('upload_file', userSession=handle, filename=os.path.basename(args.background),
                                      base64Data=base64.b64encode(f.read()).decode())
        check(err is None, 'upload_file of the background (%s)' % err)
        if err is None:
            result, err = client.call('load_spectrum_file', userSession=handle, path=result['path'],
                                      specType='Background')
            check(err is None, 'load the background (%s)' % err)

    if args.exchange_dir:
        direct = os.path.join(args.exchange_dir, sess['session_dir'], 'input', 'direct_copy.n42')
        shutil.copy(args.foreground, direct + '.tmp')
        os.rename(direct + '.tmp', direct)
        result, err = client.call('load_spectrum_file', userSession=handle, path='input/direct_copy.n42',
                                  specType='Secondary')
        check(err is None, 'load a file written directly into the session directory (%s)' % err)

    for bad_path in ['/etc/passwd', '../../etc/passwd', 'input/../../../etc/hosts']:
        _, err = client.call('load_spectrum_file', userSession=handle, path=bad_path, specType='Secondary')
        check(err is not None, 'path outside the session refused: %s' % bad_path)

    # Another session can not reach this one's files.
    other, err = client.call('create_session')
    if check(err is None, 'create another session'):
        for bad_path in ['../' + sess['session_dir'] + '/input/' + fg_name,
                         'input/../../' + sess['session_dir'] + '/input/' + fg_name]:
            _, err = client.call('load_spectrum_file', userSession=other['userSession'], path=bad_path,
                                 specType='Foreground')
            check(err is not None, "another session's file refused: %s" % bad_path)
        _, err = client.call('close_session', userSession=other['userSession'])
        check(err is None, 'close the other session')

    result, err = client.call('get_peaks', userSession=handle)
    check(err is None, 'get_peaks (%s)' % err)

    result, err = client.call('get_spectrum_image', userSession=handle)
    if check(err is None and 'image_file' in result, 'get_spectrum_image saves an image file (%s)' % err):
        status, hdrs, body = client.http('GET', result['image_file']['download_url'])
        check(status == 200 and body[:4] in (b'\x89PNG', b'<svg', b'<?xm'), 'download the spectrum image')

    if args.prompt:
        print('Assistant')
        result, err = client.call('assistant_submit_prompt', userSession=handle, prompt=args.prompt, wait_seconds=5)
        check(err is None and result.get('turn_id'), 'assistant_submit_prompt (%s)' % err)
        if err is None:
            turn_id = result['turn_id']
            _, busy_err = client.call('export_analysis', userSession=handle)
            if result['status'] == 'running':
                check(busy_err is not None, 'state-changing tools refused while the turn runs')
            deadline = time.time() + 900
            while result.get('status') == 'running' and time.time() < deadline:
                result, err = client.call('assistant_get_result', userSession=handle, turn_id=turn_id, wait_seconds=45)
                if err:
                    break
            check(err is None and result.get('status') == 'done', 'the turn finished: %s' % (err or result.get('status')))
            if result.get('response'):
                print('    response:', result['response'][:400].replace('\n', ' '), '...')

    print('Export')
    result, err = client.call('export_analysis', userSession=handle, items=['n42', 'peaks_csv', 'spectrum_html'])
    exported = {}
    if check(err is None and result.get('files'), 'export_analysis (%s)' % err):
        for f in result['files']:
            status, hdrs, body = client.http('GET', f['download_url'])
            check(status == 200 and len(body) == f['bytes'], 'download %s' % f['path'])
            exported[f['path']] = body

    result, err = client.call('list_session_files', userSession=handle)
    check(err is None and result['output'], 'list_session_files')

    analysis_n42 = next((body for path, body in exported.items() if path.endswith('_analysis.n42')), None)
    check(analysis_n42 is not None and b'<DHS:InterSpec' in analysis_n42 or (analysis_n42 and b'InterSpec' in analysis_n42),
          'the exported N42 has the InterSpec analysis')
    if args.prompt and analysis_n42:
        check(b'LlmConversationHistory' in analysis_n42, 'the exported N42 has the conversation')

    _, err = client.call('close_session', userSession=handle)
    check(err is None, 'close_session')
    _, err = client.call('get_spectrum_info', userSession=handle)
    check(err is not None and 'ended' in err, 'tools on a closed session say it ended')
    status, _, _ = client.http('PUT', upload_url + '/after.n42', b'abc')
    check(status == 404, 'upload URL is gone after close_session')
    if args.exchange_dir:
        time.sleep(3)
        check(not os.path.exists(os.path.join(args.exchange_dir, sess['session_dir'])),
              'session directory removed after close_session')

    if analysis_n42:
        print('Resume (session B)')
        sess_b, err = client.call('create_session')
        check(err is None, 'create a second session')
        if err is None:
            handle_b = sess_b['userSession']
            status, _, body = client.http('PUT', sess_b['upload_url'] + '/resume.n42', analysis_n42)
            check(status == 201, 'upload the exported N42')
            result, err = client.call('load_spectrum_file', userSession=handle_b, path='input/resume.n42',
                                      specType='Foreground', conversation='resume')
            check(err is None and 'restored' in result, 'load it with conversation="resume" (%s)' % err)
            if err is None:
                print('    restored:', json.dumps(result['restored']))
                if args.prompt:
                    check(result['restored'].get('conversation_turns', 0) > 0, 'the conversation was restored')
                    result, err = client.call('load_spectrum_file', userSession=handle_b, path='input/resume.n42',
                                              specType='Foreground', conversation='new')
                    check(err is None and result['restored'].get('conversation_turns') == 0,
                          'conversation="new" starts empty')
            _, err = client.call('close_session', userSession=handle_b)
            check(err is None, 'close the second session')

    print()
    print('%d failure(s)' % len(failures))
    for f in failures:
        print('  -', f)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
