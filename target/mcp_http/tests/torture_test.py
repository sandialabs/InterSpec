#!/usr/bin/env python3
"""Torture test of the InterSpec MCP analysis server (target/mcp_http).

Where smoke_test.py checks the main workflow, this tries to break things:

  sweep    - calls every MCP tool (all ~65 of them) with plausible arguments, against a Ba-133
             spectrum with background, a uranium HPGe spectrum (isotopics), and a passthrough
             (time history) file; plus calls that should fail.  Each call must return a result or a
             clean error within its deadline, and the session must still respond afterwards.
             Images (spectrum, zoomed, log/linear, background-subtracted, jpeg, time chart)
             are saved to --out for inspection.
  chaos    - failure injection: clients that disconnect mid-call, closing a session during a
             long call, killing the browser (and a single tab) while sessions are running,
             malformed JSON-RPC and uploads, odd file names, concurrent calls on one session,
             and rapid create/close cycles.
  llm      - (needs a working LLM on the server) several heavy assistant turns in parallel
             sessions: nuclide ID, U enrichment, and an activity fit.

Uses only the Python standard library (plus `pkill`/`pgrep` for the browser-kill tests, which are
skipped with --no-kill, e.g. for a remote server).

Example:
  INTERSPEC_MCP_BEARER_TOKEN=... ./torture_test.py --url http://127.0.0.1:8078 \\
      --repo /path/to/InterSpec --out /tmp/torture [--phases sweep,chaos,llm]
"""

import argparse
import base64
import json
import os
import socket
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request


def short(text, n=160):
    text = str(text).replace('\n', ' ')
    return text if len(text) <= n else text[:n] + '...'


class Client:
    def __init__(self, url, token):
        self.url = url.rstrip('/')
        self.token = token
        self.next_id = 1
        self.lock = threading.Lock()

    def http(self, method, url, body=None, headers=None, auth=True, timeout=600):
        hdrs = dict(headers or {})
        if auth and self.token:
            hdrs['Authorization'] = 'Bearer ' + self.token
        req = urllib.request.Request(url, data=body, method=method, headers=hdrs)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()

    def rpc(self, method, params=None, timeout=600):
        with self.lock:
            msg_id = self.next_id
            self.next_id += 1
        msg = {'jsonrpc': '2.0', 'id': msg_id, 'method': method}
        if params is not None:
            msg['params'] = params
        status, body = self.http('POST', self.url + '/mcp-api', json.dumps(msg).encode(),
                                 {'Content-Type': 'application/json', 'Accept': 'application/json'},
                                 timeout=timeout)
        if status != 200:
            return {'http_status': status, 'body': body.decode(errors='replace')}
        return json.loads(body)

    def call(self, tool, args=None, timeout=600):
        """Returns (result dict or None, images [(mime, bytes)], error or None, seconds)."""
        t0 = time.time()
        try:
            reply = self.rpc('tools/call', {'name': tool, 'arguments': args or {}}, timeout=timeout)
        except Exception as e:
            return None, [], 'transport: %r' % (e,), time.time() - t0
        elapsed = time.time() - t0
        if 'http_status' in reply:
            return None, [], 'HTTP %d %s' % (reply['http_status'], short(reply['body'])), elapsed
        if 'error' in reply:
            return None, [], reply['error'].get('message', str(reply['error'])), elapsed
        content = reply['result']['content']
        if reply['result'].get('isError'):
            return None, [], content[0]['text'], elapsed
        result = json.loads(content[0]['text'])
        images = [(c['mimeType'], base64.b64decode(c['data'])) for c in content if c.get('type') == 'image']
        return result, images, None, elapsed

    def healthy(self):
        status, body = self.http('GET', self.url + '/mcp-api/healthz', auth=False, timeout=10)
        return status == 200 and json.loads(body).get('status') == 'ok', (json.loads(body) if status == 200 else {})


class Report:
    def __init__(self, out_dir):
        self.rows = []
        self.out_dir = out_dir
        self.lock = threading.Lock()

    def add(self, phase, name, status, seconds, detail=''):
        with self.lock:
            self.rows.append({'phase': phase, 'name': name, 'status': status, 'seconds': round(seconds, 2),
                              'detail': short(detail, 300)})
        mark = {'ok': 'ok  ', 'expected-error': 'ok* ', 'FAIL': 'FAIL', 'SLOW': 'SLOW', 'info': 'info'}.get(status, status)
        print('  %s %-38s %6.1fs  %s' % (mark, name, seconds, short(detail, 110)), flush=True)

    def failures(self):
        return [r for r in self.rows if r['status'] in ('FAIL',)]

    def save(self):
        path = os.path.join(self.out_dir, 'torture_report.json')
        with open(path, 'w') as f:
            json.dump(self.rows, f, indent=1)
        return path


def upload(client, session, path, name=None):
    with open(path, 'rb') as f:
        status, body = client.http('PUT', session['upload_url'] + '/' + (name or os.path.basename(path)), f.read())
    if status != 201:
        raise RuntimeError('upload of %s failed: %d %s' % (path, status, body[:200]))
    return json.loads(body)['path']


# ------------------------------------------------------------------------------------------------
# Sweep
# ------------------------------------------------------------------------------------------------

def run_calls(client, report, phase, handle, calls, out_dir, context):
    """`calls`: list of (label, tool, args or callable(context)->args, expect), where expect is
    'ok', 'error', or 'any'.  A callable can read results of earlier calls from `context`."""
    for label, tool, args, expect in calls:
        if callable(args):
            try:
                args = args(context)
            except Exception as e:
                report.add(phase, label, 'FAIL', 0, 'could not build args: %r' % (e,))
                continue
        args = dict(args or {})
        if handle:
            args.setdefault('userSession', handle)
        result, images, err, secs = client.call(tool, args)
        context[label] = result
        for i, (mime, data) in enumerate(images):
            ext = {'image/png': '.png', 'image/jpeg': '.jpg', 'image/svg+xml': '.svg'}.get(mime, '.img')
            with open(os.path.join(out_dir, '%s_%s%s' % (phase, label, ext if i == 0 else '_%d%s' % (i, ext))), 'wb') as f:
                f.write(data)
        if err and err.startswith('transport'):
            report.add(phase, label, 'FAIL', secs, err)
        elif expect == 'ok' and err:
            report.add(phase, label, 'FAIL', secs, err)
        elif expect == 'error' and not err:
            report.add(phase, label, 'FAIL', secs, 'expected an error, got: ' + json.dumps(result))
        elif secs > 120:
            report.add(phase, label, 'SLOW', secs, err or json.dumps(result))
        else:
            report.add(phase, label, 'expected-error' if err else 'ok', secs,
                       err or (('%d image(s) ' % len(images)) if images else '') + json.dumps(result))

        # The session must still respond after every call.
        if handle:
            r, _, e, _ = client.call('get_loaded_spectra', {'userSession': handle}, timeout=60)
            if e and 'working on a turn' not in e:
                report.add(phase, label + ' [then responsive?]', 'FAIL', 0, e)


def checkpoint_name(context, label):
    res = context.get(label) or {}
    for key in ('checkpoint_name', 'checkpointName', 'name'):
        if key in res:
            return res[key]
    raise RuntimeError('no checkpoint name in ' + json.dumps(res))


def first_drf_name(context):
    res = context.get('available_detector_efficiency_functions') or {}
    text = json.dumps(res)
    for key in ('defaultAvailable', 'DefaultAvailable', 'available', 'detectors', 'drfs'):
        if key in res and res[key]:
            entry = res[key][0]
            return entry if isinstance(entry, str) else entry.get('name') or entry.get('identifier')
    # Fall back to a well-known shipped one
    return 'Detective-EX' if 'Detective-EX' in text else 'IdentiFINDER-NGH'


def first_u_preset(context):
    res = context.get('list_isotopics_presets') or {}
    names = []
    for p in res.get('presets', []):
        names.append(p.get('name') if isinstance(p, dict) else p)
    for n in names:
        if n and ('U' in n) and ('Pu' not in n):
            return n
    if names:
        return names[0]
    raise RuntimeError('no presets in ' + json.dumps(res))


def sweep(client, report, args, out_dir, mode_agent):
    phase = 'sweep'
    print('\n== Sweep of every tool', flush=True)
    ex = os.path.join(args.repo, 'example_spectra')
    td = os.path.join(args.repo, 'target', 'testing', 'test_data')
    files = {
        'fg': os.path.join(ex, 'ba133_source_640s_20100317.n42'),
        'bg': os.path.join(ex, 'background_20100317.n42'),
        'u': os.path.join(td, 'uranium_40pct_selfatten_HPGe_15cm.n42'),
        'tp': os.path.join(td, 'llm_time_history', 'Q345_Racetrack.n42'),
    }

    handle = None
    paths = dict(files)
    if mode_agent:
        sess, _, err, secs = client.call('create_session')
        if err:
            report.add(phase, 'create_session', 'FAIL', secs, err)
            return
        report.add(phase, 'create_session', 'ok', secs, sess['session_dir'])
        handle = sess['userSession']
        paths = {k: upload(client, sess, v) for k, v in files.items()}

    ctx = {}
    seen = set()

    ba_calls = [
        ('load fg', 'load_spectrum_file', {'path': paths['fg'], 'specType': 'Foreground'}, 'ok'),
        ('load bg', 'load_spectrum_file', {'path': paths['bg'], 'specType': 'Background'}, 'ok'),
        ('get_loaded_spectra', 'get_loaded_spectra', {}, 'ok'),
        ('get_spectrum_info', 'get_spectrum_info', {'fileInfo': True}, 'ok'),
        ('get_spectrum_info bg', 'get_spectrum_info', {'specType': 'Background'}, 'ok'),
        ('get_counts_in_energy_range', 'get_counts_in_energy_range', {'lowerEnergy': 345, 'upperEnergy': 365, 'maxChannels': 1}, 'ok'),
        ('get_counts_in_energy_range 32', 'get_counts_in_energy_range', {'lowerEnergy': 50, 'upperEnergy': 450, 'maxChannels': 32}, 'ok'),
        ('get_expected_fwhm', 'get_expected_fwhm', {'energy': 356}, 'ok'),
        ('get_automated_id_results', 'get_automated_id_results', {}, 'any'),
        ('get_peaks', 'get_peaks', {}, 'ok'),
        ('get_peaks unidentified', 'get_peaks', {'filter': 'unidentified', 'maxResults': 10}, 'ok'),
        ('get_peaks elevated', 'get_peaks', {'filter': 'elevated_above_background', 'lowerEnergy': 20, 'upperEnergy': 500}, 'ok'),
        ('get_peaks bg', 'get_peaks', {'specType': 'Background'}, 'ok'),
        ('add_analysis_peaks_for_source', 'add_analysis_peaks_for_source', {'source': 'Ba133'}, 'ok'),
        ('add_analysis_peaks_for_source NORM', 'add_analysis_peaks_for_source', {'source': None, 'options': ['FitNormPeaks']}, 'any'),
        ('add_analysis_peaks_for_source list', 'add_analysis_peaks_for_source', {'source': ['Ba133', 'Cs137'], 'options': ['RefitInterferingAnalysisPeaks']}, 'any'),
        ('add_analysis_peak', 'add_analysis_peak', {'energy': 1460.8, 'source': 'K40'}, 'any'),
        ('add_analysis_peak no-add', 'add_analysis_peak', {'energy': 609.3, 'DoNotAddToAnalysisPeaks': True}, 'any'),
        ('get_peaks analysis', 'get_peaks', {'filter': 'analysis'}, 'ok'),
        ('get_identified_sources', 'get_identified_sources', {}, 'ok'),
        ('create_peak_checkpoint', 'create_peak_checkpoint', {'name_hint': 'torture'}, 'ok'),
        ('edit_analysis_peak', 'edit_analysis_peak', {'energy': 356, 'refit': True, 'color': 'rgb(200,0,0)'}, 'ok'),
        ('edit_analysis_peak amp', 'edit_analysis_peak', {'energy': 356, 'amplitude': 100000, 'fitForAmplitude': False}, 'any'),
        ('set_peak_shape', 'set_peak_shape', {'energy': 356, 'continuumType': 'Quadratic'}, 'any'),
        ('edit_analysis_peak delete', 'edit_analysis_peak', {'energy': 383.8, 'action': 'DeletePeak'}, 'any'),
        ('restore_peaks_to_checkpoint', 'restore_peaks_to_checkpoint', lambda c: {'checkpoint_name': checkpoint_name(c, 'create_peak_checkpoint')}, 'ok'),
        ('compton_scatter_peak_check', 'compton_scatter_peak_check', {'energy': 356}, 'ok'),
        ('escape_peak_check', 'escape_peak_check', {'energy': 356}, 'ok'),
        ('sum_peak_check', 'sum_peak_check', {'energy': 437, 'distance': '25 cm'}, 'any'),
        ('beta_continuum_check', 'beta_continuum_check', {}, 'any'),
        ('search_sources_by_energy', 'search_sources_by_energy', {'energies': [{'energy': 356}, {'energy': 81}]}, 'ok'),
        ('search_sources_by_energy x', 'search_sources_by_energy', {'energies': [{'energy': 72.8}], 'source_category': 'Fluorescence x-rays'}, 'any'),
        ('sources_with_gammas_near_energy', 'sources_with_gammas_near_energy', {'energy': 356, 'window': 2}, 'ok'),
        ('source_info', 'source_info', {'source': 'Ba133'}, 'ok'),
        ('source_info reaction', 'source_info', {'source': 'H(n,g)'}, 'any'),
        ('nuclide_decay_chain', 'nuclide_decay_chain', {'nuclide': 'U238'}, 'ok'),
        ('decay_calculator', 'decay_calculator', {'nuclide': 'Cs137', 'activity': '1 uCi', 'time_duration': '10 years'}, 'ok'),
        ('source_photons', 'source_photons', {'source': 'Ba133', 'maxResults': 20}, 'ok'),
        ('source_photons U', 'source_photons', {'source': 'U238', 'age': '20 years', 'includeCascadeSumEnergies': True, 'maxResults': 10}, 'ok'),
        ('currie_mda_calc', 'currie_mda_calc', {'energy': 661.7, 'nuclide': 'Cs137', 'distance': '50 cm'}, 'any'),
        ('calculate_dose', 'calculate_dose', {'nuclide': 'Cs137', 'activity': '1 mCi', 'distance': '1 m'}, 'ok'),
        ('calculate_dose shielded', 'calculate_dose', {'nuclide': 'Co60', 'activity': '10 mCi', 'distance': '2 m', 'material': 'Fe', 'thickness': '2 cm'}, 'ok'),
        ('get_materials', 'get_materials', {}, 'ok'),
        ('get_material_info', 'get_material_info', {'material': 'Fe'}, 'ok'),
        ('available_detector_efficiency_functions', 'available_detector_efficiency_functions', {}, 'ok'),
        ('load_detector_efficiency_function', 'load_detector_efficiency_function', lambda c: {'identifier': first_drf_name(c)}, 'any'),
        ('load_detector_efficiency_function formula', 'load_detector_efficiency_function',
         {'formula': 'exp(-1.2 - 0.9*log(x/100))', 'detectorDiameter': '5 cm', 'geometryType': 'FarFieldIntrinsic', 'drfName': 'torture formula', 'energyUnits': 'keV'}, 'any'),
        ('detector_efficiency_function_info', 'detector_efficiency_function_info', {}, 'any'),
        ('photopeak_detection_efficiency', 'photopeak_detection_efficiency',
         {'Energies': [81, 356, 662], 'Distance': '25 cm', 'Shielding': [{'Material': 'Fe', 'Thickness': '1 cm'}, {'AD': 2.0, 'AN': 26}]}, 'any'),
        ('mark_peaks_for_activity_fit', 'mark_peaks_for_activity_fit', {'peak_energies': [81, 276.4, 302.9, 356, 383.8], 'use_for_fit': True}, 'any'),
        ('get_shielding_source_config', 'get_shielding_source_config', {}, 'any'),
        ('modify_shielding set_distance', 'modify_shielding_source_config', {'operation': 'set_distance', 'distance': '50 cm'}, 'any'),
        ('modify_shielding set_geometry', 'modify_shielding_source_config', {'operation': 'set_geometry', 'geometry': 'Spherical'}, 'any'),
        ('modify_shielding add_shielding', 'modify_shielding_source_config', {'operation': 'add_shielding', 'material': 'Fe', 'thickness': '0.2 cm'}, 'any'),
        ('modify_shielding set_fit_options', 'modify_shielding_source_config', {'operation': 'set_fit_options', 'attenuate_for_air': True}, 'any'),
        ('activity_fit', 'activity_fit', {'peak_energies': [81, 302.9, 356, 383.8]}, 'any'),
        ('activity_fit_one_off', 'activity_fit_one_off', {'distance': '50 cm', 'geometry': 'Spherical', 'peak_energies': [81, 356]}, 'any'),
        ('peak_based_relative_efficiency', 'peak_based_relative_efficiency', {'sources': ['Ba133']}, 'any'),
        ('save_energy_cal_checkpoint', 'save_energy_cal_checkpoint', {'name_hint': 'torture'}, 'ok'),
        ('fit_energy_calibration', 'fit_energy_calibration', {'specType': 'Foreground', 'peak_energies': [81, 302.9, 356, 383.8]}, 'any'),
        ('restore_energy_cal_checkpoint', 'restore_energy_cal_checkpoint', lambda c: {'checkpoint_name': checkpoint_name(c, 'save_energy_cal_checkpoint')}, 'ok'),
        ('get_spectrum_image', 'get_spectrum_image', {}, 'ok'),
        ('get_spectrum_image zoom linear', 'get_spectrum_image', {'energy_range': {'min_keV': 250, 'max_keV': 450}, 'y_axis_scale': 'linear'}, 'ok'),
        ('get_spectrum_image bgsub log', 'get_spectrum_image', {'background_subtract': True, 'y_axis_scale': 'log', 'energy_range': {'min_keV': 20, 'max_keV': 120}}, 'ok'),
        ('ERR svg image', 'get_spectrum_image', {'format': 'svg'}, 'error'),
        ('get_spectrum_image jpeg small', 'get_spectrum_image', {'format': 'jpeg', 'max_size': 400}, 'ok'),
        ('get_spectrum_image after zoom', 'get_spectrum_image', {}, 'ok'),
        ('unload secondary', 'load_spectrum_file', {'specType': 'Secondary'}, 'ok'),
    ]
    if mode_agent:
        ba_calls += [
            ('list_session_files', 'list_session_files', {}, 'ok'),
            ('export_analysis', 'export_analysis', {'items': ['n42', 'peaks_csv', 'spectrum_html'], 'basename': 'torture ba133'}, 'ok'),
            ('get_usage_guide', 'get_usage_guide', {}, 'ok'),
            ('upload_file', 'upload_file', {'filename': 'notes.txt', 'base64Data': base64.b64encode(b'hello').decode()}, 'ok'),
        ]

    u_calls = [
        ('load U', 'load_spectrum_file', {'path': paths['u'], 'specType': 'Foreground'}, 'ok'),
        ('unload bg', 'load_spectrum_file', {'specType': 'Background'}, 'ok'),
        ('list_isotopics_presets', 'list_isotopics_presets', {}, 'ok'),
        ('load_isotopics_preset', 'load_isotopics_preset', lambda c: {'preset': first_u_preset(c)}, 'ok'),
        ('get_isotopics_config', 'get_isotopics_config', {}, 'ok'),
        ('get_isotopics_config_schema', 'get_isotopics_config_schema', {}, 'ok'),
        ('modify_isotopics_options', 'modify_isotopics_options', {'note': 'torture test'}, 'ok'),
        ('modify_isotopics_curve_settings', 'modify_isotopics_curve_settings', {'rel_eff_eqn_order': 3}, 'any'),
        ('modify_isotopics_nuclides', 'modify_isotopics_nuclides', {'action': 'update', 'nuclides': [{'name': 'U235', 'age': '20 years'}]}, 'any'),
        ('modify_isotopics_constraints', 'modify_isotopics_constraints', {'action': 'clear_all', 'constraint_type': 'activity_ratio'}, 'any'),
        ('modify_isotopics_rois', 'modify_isotopics_rois', {'action': 'add', 'rois': [{'lower_energy': 1995, 'upper_energy': 2010}]}, 'any'),
        ('perform_isotopics_calculation', 'perform_isotopics_calculation', {}, 'any'),
        ('get_spectrum_image U', 'get_spectrum_image', {'energy_range': {'min_keV': 80, 'max_keV': 220}}, 'ok'),
        ('reset_isotopics_config', 'reset_isotopics_config', {}, 'ok'),
    ]

    tp_calls = [
        ('load passthrough', 'load_spectrum_file', {'path': paths['tp'], 'specType': 'Foreground'}, 'ok'),
        ('get_time_history_count_rates', 'get_time_history_count_rates', {'maxEntries': 40}, 'ok'),
        ('get_time_history_count_rates roi', 'get_time_history_count_rates', {'maxEntries': 40, 'energyMin': 600, 'energyMax': 700}, 'any'),
        ('get_displayed_time_ranges', 'get_displayed_time_ranges', {}, 'ok'),
        ('set_displayed_time_ranges', 'set_displayed_time_ranges', {'foreground': '0-20'}, 'any'),
        ('get_time_chart_image', 'get_time_chart_image', {}, 'ok'),
        ('get_spectrum_image passthrough', 'get_spectrum_image', {}, 'ok'),
    ]

    bad_calls = [
        ('ERR missing file', 'load_spectrum_file', {'path': 'input/no_such_file.n42', 'specType': 'Foreground'}, 'error'),
        ('ERR bad specType', 'load_spectrum_file', {'path': paths['fg'], 'specType': 'Primary'}, 'error'),
        ('ERR missing arg', 'get_expected_fwhm', {}, 'error'),
        ('ERR wrong type', 'get_expected_fwhm', {'energy': 'abc'}, 'error'),
        ('ERR bad nuclide', 'source_info', {'source': 'NotANuclide99'}, 'any'),
        ('ERR bad checkpoint', 'restore_peaks_to_checkpoint', {'checkpoint_name': 'no-such-checkpoint'}, 'error'),
        ('ERR bad format', 'get_spectrum_image', {'format': 'bmp'}, 'error'),
        ('ERR bad activity', 'decay_calculator', {'nuclide': 'Cs137', 'activity': 'lots'}, 'error'),
        ('ERR negative energy', 'add_analysis_peak', {'energy': -50}, 'error'),
        ('ERR unknown tool', 'no_such_tool', {}, 'error'),
        ('ERR bad preset', 'load_isotopics_preset', {'preset': '../../../etc/passwd'}, 'error'),
        ('ERR drf path', 'load_detector_efficiency_function', {'identifier': '/etc/passwd', 'source': 'FilePath'}, 'error'),
    ]
    if mode_agent:
        bad_calls += [
            ('ERR no userSession', 'get_peaks', {'userSession': ''}, 'error'),
            ('ERR bogus userSession', 'get_peaks', {'userSession': 'NotARealSessionHandle00000'}, 'error'),
            ('ERR traversal', 'load_spectrum_file', {'path': '../../../../etc/hosts', 'specType': 'Foreground'}, 'error'),
            ('ERR image traversal', 'assistant_submit_prompt', {'prompt': 'x', 'images': [{'path': '../../../../etc/hosts'}]}, 'any'),
            ('ERR export nothing', 'export_analysis', {'items': ['bogus']}, 'error'),
            ('ERR bad turn', 'assistant_get_result', {'turn_id': 'turn-nope'}, 'error'),
            ('ERR upload bad name', 'upload_file', {'filename': '../../evil.sh', 'base64Data': base64.b64encode(b'#!/bin/sh').decode()}, 'ok'),
        ]

    for group, calls in (('ba133', ba_calls), ('uranium', u_calls), ('passthrough', tp_calls), ('errors', bad_calls)):
        print('-- %s' % group, flush=True)
        run_calls(client, report, phase, handle, calls, out_dir, ctx)
        for c in calls:
            seen.add(c[1])

    # Which tools were never called?
    tools = {t['name'] for t in client.rpc('tools/list')['result']['tools']}
    untested = sorted(tools - seen - {'create_session', 'close_session', 'assistant_submit_prompt',
                                       'assistant_get_result', 'assistant_cancel', 'invoke_agent'})
    report.add(phase, 'coverage', 'info', 0, '%d tools; not swept here: %s' % (len(tools), ', '.join(untested) or 'none'))

    if mode_agent:
        # The 'evil.sh' upload must have landed, sanitized, inside input/
        files_now, _, err, _ = client.call('list_session_files', {'userSession': handle})
        names = [f['path'] for f in (files_now or {}).get('input', [])]
        ok = any(n == 'input/evil.sh' for n in names) and not any('..' in n for n in names)
        report.add(phase, 'sanitized upload name', 'ok' if ok else 'FAIL', 0, ', '.join(names))
        _, _, err, secs = client.call('close_session', {'userSession': handle})
        report.add(phase, 'close_session', 'FAIL' if err else 'ok', secs, err or '')


# ------------------------------------------------------------------------------------------------
# Chaos
# ------------------------------------------------------------------------------------------------

def raw_post_then_disconnect(url, token, body, wait_s):
    """Sends a request, waits a bit, and closes the socket without reading the response."""
    u = urllib.parse.urlparse(url)
    sock = socket.create_connection((u.hostname, u.port or 80), timeout=10)
    data = json.dumps(body).encode()
    req = ('POST %s HTTP/1.1\r\nHost: %s\r\nAuthorization: Bearer %s\r\nContent-Type: application/json\r\n'
           'Accept: application/json\r\nContent-Length: %d\r\nConnection: close\r\n\r\n'
           % (u.path or '/mcp-api', u.netloc, token, len(data))).encode() + data
    sock.sendall(req)
    time.sleep(wait_s)
    sock.close()


def new_session_with_spectrum(client, args, which='fg'):
    sess, _, err, _ = client.call('create_session')
    if err:
        raise RuntimeError('create_session: ' + err)
    path = {'fg': os.path.join(args.repo, 'example_spectra', 'ba133_source_640s_20100317.n42'),
            'u': os.path.join(args.repo, 'target', 'testing', 'test_data', 'uranium_40pct_selfatten_HPGe_15cm.n42')}[which]
    p = upload(client, sess, path)
    _, _, err, _ = client.call('load_spectrum_file', {'userSession': sess['userSession'], 'path': p, 'specType': 'Foreground'})
    if err:
        raise RuntimeError('load: ' + err)
    return sess


def chaos(client, report, args, out_dir):
    phase = 'chaos'
    print('\n== Failure injection', flush=True)
    url = client.url + '/mcp-api'

    def check_health(label):
        ok, h = client.healthy()
        report.add(phase, label, 'ok' if ok else 'FAIL', 0, json.dumps(h))
        return ok

    # 1. Malformed requests
    print('-- malformed requests', flush=True)
    cases = [
        ('invalid JSON', b'{"jsonrpc":"2.0","id":1,', 'application/json'),
        ('no method', b'{"jsonrpc":"2.0","id":1}', 'application/json'),
        ('params not object', b'{"jsonrpc":"2.0","id":1,"method":"tools/call","params":[1,2]}', 'application/json'),
        ('arguments not object', b'{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"get_materials","arguments":"x"}}', 'application/json'),
        ('no tool name', b'{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{}}', 'application/json'),
        ('unknown method', b'{"jsonrpc":"2.0","id":1,"method":"resources/list"}', 'application/json'),
        ('binary garbage', bytes(range(256)) * 4, 'application/octet-stream'),
        ('empty body', b'', 'application/json'),
    ]
    for label, body, ctype in cases:
        t0 = time.time()
        status, resp = client.http('POST', url, body, {'Content-Type': ctype}, timeout=30)
        ok = status in (200, 400, 405) and (status != 200 or b'"error"' in resp)
        report.add(phase, 'malformed: ' + label, 'ok' if ok else 'FAIL', time.time() - t0, '%d %s' % (status, short(resp.decode(errors='replace'), 120)))
    status, _ = client.http('POST', url, b'{}', {'Content-Type': 'application/json'}, auth=False)
    report.add(phase, 'no auth -> 401', 'ok' if status == 401 else 'FAIL', 0, str(status))
    check_health('healthy after malformed requests')

    # 2. Odd uploads
    print('-- uploads', flush=True)
    sess, _, err, _ = client.call('create_session')
    h = sess['userSession']
    up = sess['upload_url']
    odd_names = ['spaces in name.n42', 'ünïcödé-名前.spe', '..', '.hidden', 'a' * 300 + '.n42', '%2e%2e%2fescape.n42', 'nul%00byte.txt']
    for name in odd_names:
        status, body = client.http('PUT', up + '/' + urllib.parse.quote(name, safe='%'), b'data')
        ok = status in (201, 400, 404)
        stored = json.loads(body).get('path', '') if status == 201 else ''
        ok = ok and ('..' not in stored) and (stored == '' or stored.startswith('input/'))
        report.add(phase, 'upload name ' + short(name, 30), 'ok' if ok else 'FAIL', 0, '%d %s' % (status, stored or short(body, 80)))
    status, _ = client.http('PUT', up, b'data')
    report.add(phase, 'PUT without filename -> 400', 'ok' if status == 400 else 'FAIL', 0, str(status))
    boundary = 'XyZtorture'
    mp = (('--%s\r\nContent-Disposition: form-data; name="f1"; filename="a.txt"\r\n\r\nfirst\r\n'
           '--%s\r\nContent-Disposition: form-data; name="f2"; filename="../b.txt"\r\n\r\nsecond\r\n--%s--\r\n')
          % (boundary, boundary, boundary)).encode()
    status, body = client.http('POST', up, mp, {'Content-Type': 'multipart/form-data; boundary=' + boundary})
    report.add(phase, 'multipart POST of 2 files', 'ok' if status == 201 else 'FAIL', 0, '%d %s' % (status, short(body, 120)))
    status, _ = client.http('GET', up + '/a.txt')
    report.add(phase, 'GET on upload URL -> 405', 'ok' if status == 405 else 'FAIL', 0, str(status))
    r, _, err, _ = client.call('upload_file', {'userSession': h, 'filename': 'x.bin', 'base64Data': '!!!not base64!!!'})
    report.add(phase, 'upload_file invalid base64', 'ok', 0, err or json.dumps(r))
    big = os.urandom(300 * 1024 * 1024 // 100)  # 3 MB; limits are tested with the server's configured caps
    status, body = client.http('PUT', up + '/big.bin', big)
    report.add(phase, 'PUT 3 MB', 'ok' if status == 201 else 'FAIL', 0, str(status))
    r, _, err, _ = client.call('export_analysis', {'userSession': h})
    report.add(phase, 'export with no spectrum -> error', 'ok' if err else 'FAIL', 0, err or '')
    status, _ = client.http('GET', client.url + '/mcp-api/files/NotAToken/x')
    report.add(phase, 'unknown file token -> 404', 'ok' if status == 404 else 'FAIL', 0, str(status))
    client.call('close_session', {'userSession': h})

    # 3. Client disconnects mid-call
    print('-- disconnects', flush=True)
    sess = new_session_with_spectrum(client, args, 'u')
    h = sess['userSession']
    presets, _, _, _ = client.call('list_isotopics_presets', {'userSession': h})
    client.call('load_isotopics_preset', {'userSession': h, 'preset': first_u_preset({'list_isotopics_presets': presets})})
    # An extra ROI happens to make this problem solve in ~100 s rather than ~450 s (Release, M-series
    #  Mac); still long enough to hold the session through the tests below.
    client.call('modify_isotopics_rois', {'userSession': h, 'action': 'add', 'rois': [{'lower_energy': 1995, 'upper_energy': 2010}]})
    for tool, targs in (('perform_isotopics_calculation', {}), ('get_spectrum_image', {}),
                        ('add_analysis_peaks_for_source', {'source': 'U235'})):
        t0 = time.time()
        targs = dict(targs, userSession=h)
        raw_post_then_disconnect(url, client.token, {'jsonrpc': '2.0', 'id': 9, 'method': 'tools/call',
                                                     'params': {'name': tool, 'arguments': targs}}, 0.3)
        # The abandoned call keeps running (a synchronous one holds the session, possibly for many
        #  minutes); the session must answer once it is done.
        deadline = time.time() + 1200
        err = 'not tried'
        while time.time() < deadline:
            r, _, err, _ = client.call('get_loaded_spectra', {'userSession': h}, timeout=1200)
            if not err:
                break
            time.sleep(2)
        report.add(phase, 'disconnect during ' + tool, 'FAIL' if err else 'ok', time.time() - t0, err or 'session still responds')
    check_health('healthy after disconnects')

    # 4. A call made while a long synchronous call holds the session waits for it, then succeeds;
    #  and the session survives the long call.
    print('-- call during a long synchronous call', flush=True)
    box = {}

    def long_iso():
        box['iso'] = client.call('perform_isotopics_calculation', {'userSession': h}, timeout=1800)
    t = threading.Thread(target=long_iso)
    t0 = time.time()
    t.start()
    time.sleep(1)
    r, _, err, secs = client.call('get_loaded_spectra', {'userSession': h}, timeout=1800)
    t.join(1800)
    iso = box.get('iso', (None, [], 'no result', 0))
    report.add(phase, 'read call queued behind a long call', 'FAIL' if err else 'ok', secs,
               err or ('returned after %.0f s; the long call took %.0f s (%s)' % (secs, iso[3], iso[2] or 'ok')))
    r, _, err, _ = client.call('get_loaded_spectra', {'userSession': h}, timeout=120)
    report.add(phase, 'session alive after the long call', 'FAIL' if err else 'ok', time.time() - t0, err or '')

    # 5. close_session while a long call runs
    print('-- close during a call', flush=True)
    result_box = {}

    def long_call():
        result_box['r'] = client.call('perform_isotopics_calculation', {'userSession': h}, timeout=1800)
    t = threading.Thread(target=long_call)
    t.start()
    time.sleep(0.5)
    _, _, err, _ = client.call('close_session', {'userSession': h})
    t.join(1800)
    r, _, call_err, secs = result_box.get('r', (None, [], 'no result', 0))
    report.add(phase, 'close_session during a long call', 'FAIL' if t.is_alive() else 'ok', secs,
               'call returned: ' + (call_err or 'a result'))
    _, _, err, _ = client.call('get_peaks', {'userSession': h})
    report.add(phase, 'closed session refuses calls', 'ok' if err and 'ended' in err else 'FAIL', 0, err or '')
    check_health('healthy after close-during-call')

    # 6. Concurrent calls on one session
    print('-- concurrent calls on one session', flush=True)
    sess = new_session_with_spectrum(client, args, 'fg')
    h = sess['userSession']
    outcomes = []

    def hammer(i):
        tool, targs = [('get_spectrum_image', {}), ('get_peaks', {}), ('get_counts_in_energy_range', {'lowerEnergy': 300, 'upperEnergy': 400}),
                       ('add_analysis_peak', {'energy': 356, 'source': 'Ba133'}), ('get_spectrum_info', {})][i % 5]
        r, _, err, secs = client.call(tool, dict(targs, userSession=h), timeout=400)
        outcomes.append((tool, err, secs))
    ts = [threading.Thread(target=hammer, args=(i,)) for i in range(20)]
    t0 = time.time()
    [x.start() for x in ts]
    [x.join(450) for x in ts]
    hung = [x for x in ts if x.is_alive()]
    errs = [(tool, e) for tool, e, _ in outcomes if e]
    report.add(phase, '20 concurrent calls on one session', 'FAIL' if hung else 'ok', time.time() - t0,
               '%d hung, %d errors: %s' % (len(hung), len(errs), short(errs, 200)))
    client.call('close_session', {'userSession': h})

    # 7. Rapid create/close cycles
    print('-- create/close cycles', flush=True)
    t0 = time.time()
    failures = 0
    for i in range(15):
        s, _, err, _ = client.call('create_session')
        if err:
            failures += 1
            continue
        _, _, err, _ = client.call('close_session', {'userSession': s['userSession']})
        failures += bool(err)
    report.add(phase, '15 create/close cycles', 'FAIL' if failures else 'ok', time.time() - t0, '%d failures' % failures)
    time.sleep(5)
    ok, hinfo = client.healthy()
    report.add(phase, 'no sessions left', 'ok' if ok and hinfo.get('sessions') == 0 and hinfo.get('starting_sessions') == 0 else 'FAIL', 0, json.dumps(hinfo))

    if args.no_kill:
        return

    # 8. Kill one tab's renderer, then the whole browser, while sessions are running.
    print('-- browser crashes', flush=True)
    sess_a = new_session_with_spectrum(client, args, 'fg')
    sess_b = new_session_with_spectrum(client, args, 'fg')
    renderers = subprocess.run(['pgrep', '-f', 'interspec-browser.*--type=renderer'], capture_output=True, text=True).stdout.split()
    if renderers:
        subprocess.run(['kill', '-9', renderers[-1]])
        time.sleep(4)
        ended = []
        for s in (sess_a, sess_b):
            _, _, err, _ = client.call('get_loaded_spectra', {'userSession': s['userSession']})
            ended.append(bool(err))
        report.add(phase, 'one renderer killed', 'ok' if any(ended) else 'info', 0,
                   'sessions ended: %s (a killed renderer may host one tab, several, or a utility)' % ended)
    else:
        report.add(phase, 'one renderer killed', 'info', 0, 'no renderer process found (remote server?)')

    browser = subprocess.run(['pgrep', '-f', 'interspec-browser'], capture_output=True, text=True).stdout.split()
    main_pid = None
    for pid in browser:
        cmd = subprocess.run(['ps', '-o', 'command=', '-p', pid], capture_output=True, text=True).stdout
        if '--type=' not in cmd:
            main_pid = pid
    if main_pid:
        subprocess.run(['kill', '-9', main_pid])
        time.sleep(4)
        statuses = []
        for s in (sess_a, sess_b):
            _, _, err, _ = client.call('get_loaded_spectra', {'userSession': s['userSession']})
            statuses.append(err or 'still alive')
        report.add(phase, 'browser killed -> sessions end', 'ok' if all('ended' in x or 'no longer' in x for x in statuses) else 'FAIL', 0, short(statuses, 200))
        t0 = time.time()
        s, _, err, _ = client.call('create_session')
        report.add(phase, 'new session after browser crash', 'FAIL' if err else 'ok', time.time() - t0, err or 'browser restarted')
        if not err:
            client.call('close_session', {'userSession': s['userSession']})
    check_health('healthy after browser crashes')


# ------------------------------------------------------------------------------------------------
# Heavy LLM turns
# ------------------------------------------------------------------------------------------------

def llm(client, report, args, out_dir):
    phase = 'llm'
    print('\n== Heavy assistant turns, in parallel', flush=True)
    jobs = [
        ('nuclide ID', 'fg', 'Identify all radioactive sources in this spectrum, including any background (NORM). '
                             'Fit the peaks needed, and justify each identification briefly.'),
        ('U enrichment', 'u', 'This is an HPGe measurement of uranium. Determine the U-235 enrichment using the isotopics tools, '
                              'and give the result with its uncertainty.'),
        ('activity fit', 'fg', 'The source is 50 cm from the detector face, unshielded. Fit the Ba-133 peaks, then use the '
                               'activity fit to estimate the Ba-133 activity. Report the activity with its uncertainty.'),
    ]

    def run(label, which, prompt):
        t0 = time.time()
        try:
            sess = new_session_with_spectrum(client, args, which)
            h = sess['userSession']
            r, _, err, _ = client.call('assistant_submit_prompt', {'userSession': h, 'prompt': prompt, 'wait_seconds': 20})
            if err:
                report.add(phase, label, 'FAIL', time.time() - t0, err)
                return
            turn = r['turn_id']
            # While it runs: read-only calls work, state-changing ones are refused.
            if r['status'] == 'running':
                _, _, e1, _ = client.call('get_spectrum_info', {'userSession': h})
                _, _, e2, _ = client.call('add_analysis_peak', {'userSession': h, 'energy': 356})
                report.add(phase, label + ': busy rules', 'ok' if (not e1 and e2) else 'FAIL', 0,
                           'read-only: %s; state-changing: %s' % (e1 or 'ok', short(e2 or 'allowed!', 80)))
            deadline = time.time() + 1500
            while r.get('status') == 'running' and time.time() < deadline:
                r, _, err, _ = client.call('assistant_get_result', {'userSession': h, 'turn_id': turn, 'wait_seconds': 45})
                if err:
                    break
            status = 'ok' if (not err and r.get('status') == 'done') else 'FAIL'
            report.add(phase, label, status, time.time() - t0,
                       err or ('%d tool calls; %s' % (len(r.get('tool_calls', [])), short(r.get('response') or r.get('error'), 200))))
            with open(os.path.join(out_dir, 'llm_%s.json' % label.replace(' ', '_')), 'w') as f:
                json.dump(r, f, indent=1)
            client.call('close_session', {'userSession': h})
        except Exception as e:
            report.add(phase, label, 'FAIL', time.time() - t0, repr(e))

    def run_cancel():
        """A cancelled turn must stop (not keep calling tools), and the session take the next prompt."""
        label = 'cancel'
        t0 = time.time()
        try:
            sess = new_session_with_spectrum(client, args, 'fg')
            h = sess['userSession']
            r, _, err, _ = client.call('assistant_submit_prompt', {'userSession': h, 'wait_seconds': 15, 'prompt':
                'Identify every radioactive source in this spectrum. Fit all peaks you need, one source at a time.'})
            if err or r.get('status') != 'running':
                report.add(phase, label, 'FAIL', time.time() - t0, err or 'turn finished before it could be cancelled')
                return
            turn = r['turn_id']
            c, _, err, _ = client.call('assistant_cancel', {'userSession': h, 'turn_id': turn})
            report.add(phase, label + ': assistant_cancel', 'ok' if (not err and c.get('cancelled')) else 'FAIL',
                       0, err or json.dumps(c))
            r, _, err, secs = client.call('assistant_get_result', {'userSession': h, 'turn_id': turn, 'wait_seconds': 10})
            report.add(phase, label + ': result is cancelled', 'ok' if (not err and r.get('status') == 'cancelled') else 'FAIL',
                       secs, err or json.dumps(r)[:200])
            # The turn is over: nothing changes the peaks any more, and state-changing calls are allowed.
            p1, _, _, _ = client.call('get_peaks', {'userSession': h})
            time.sleep(20)
            p2, _, _, _ = client.call('get_peaks', {'userSession': h})
            n1, n2 = len((p1 or {}).get('peaks', [])), len((p2 or {}).get('peaks', []))
            report.add(phase, label + ': turn stopped', 'ok' if n1 == n2 else 'FAIL', 20,
                       'peaks %d then %d, 20 s apart' % (n1, n2))
            _, _, err, _ = client.call('add_analysis_peak', {'userSession': h, 'energy': 356})
            report.add(phase, label + ': session no longer busy', 'FAIL' if (err and 'working on a turn' in err) else 'ok',
                       0, err or '')
            r, _, err, _ = client.call('assistant_submit_prompt', {'userSession': h, 'wait_seconds': 45,
                'prompt': 'In one short sentence, and without calling any tools: what is the half-life of Cs-137?'})
            while not err and r.get('status') == 'running' and time.time() < t0 + 600:
                r, _, err, _ = client.call('assistant_get_result', {'userSession': h, 'turn_id': r['turn_id'], 'wait_seconds': 45})
            answer = (r or {}).get('response') or ''
            report.add(phase, label + ': next prompt answered', 'ok' if (not err and r.get('status') == 'done' and '30' in answer) else 'FAIL',
                       time.time() - t0, err or short(answer, 200))
            client.call('close_session', {'userSession': h})
        except Exception as e:
            report.add(phase, label, 'FAIL', time.time() - t0, repr(e))

    ts = [threading.Thread(target=run, args=j) for j in jobs] + [threading.Thread(target=run_cancel)]
    [t.start() for t in ts]
    [t.join() for t in ts]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--url', default='http://127.0.0.1:8078')
    parser.add_argument('--token', default=os.environ.get('INTERSPEC_MCP_BEARER_TOKEN', ''))
    parser.add_argument('--repo', default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
    parser.add_argument('--out', default='torture_out')
    parser.add_argument('--phases', default='sweep,chaos')
    parser.add_argument('--desktop', action='store_true',
                        help='Server is NOT an agent server (sweep only; uses the most recent session).')
    parser.add_argument('--no-kill', action='store_true', help='Skip the browser-kill tests.')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    client = Client(args.url, args.token)
    report = Report(args.out)
    phases = args.phases.split(',')
    t0 = time.time()

    if 'sweep' in phases:
        sweep(client, report, args, args.out, not args.desktop)
    if 'chaos' in phases and not args.desktop:
        chaos(client, report, args, args.out)
    if 'llm' in phases and not args.desktop:
        llm(client, report, args, args.out)

    path = report.save()
    fails = report.failures()
    slow = [r for r in report.rows if r['status'] == 'SLOW']
    print('\n%d checks in %.0f s; %d failure(s), %d slow.  Report: %s' % (len(report.rows), time.time() - t0, len(fails), len(slow), path))
    for r in fails + slow:
        print('  %s %s/%s: %s' % (r['status'], r['phase'], r['name'], r['detail']))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
