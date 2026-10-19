#!/usr/bin/env python3
"""Checks that `load_spectrum_file` loads files as the InterSpec GUI does, and reports what it did.

Writes small N42 files with a known layout (one background and ten foreground samples from three
gamma detectors and a neutron detector, as 60 s samples and as a search-mode portal file; a second
file from the same instrument; a file with two backgrounds; a file with a factory background and a
foreground), uploads them with real example files, and checks the samples, detectors and background
chosen for each load, plus the `summary`/`displayed`/`samplesNotDisplayed` the client gets back.

Example:
  INTERSPEC_MCP_BEARER_TOKEN=... ./load_test.py --url http://127.0.0.1:8078 --repo /path/to/InterSpec
"""

import argparse
import json
import math
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from torture_test import Client, upload  # noqa: E402

NCHAN = 1024
GAIN = 3.0


def spectrum(rng, live, cs137=0.0):
    """Counts for `live` seconds of a falling continuum, plus a 662 keV peak with `cs137` counts/s."""
    counts = []
    for ch in range(NCHAN):
        energy = (ch + 0.5) * GAIN
        mean = live * (40.0 * math.exp(-energy / 250.0) + 0.5)
        if cs137:
            sigma = 0.035 * 662.0 / 2.355
            mean += live * cs137 * GAIN * math.exp(-0.5 * ((energy - 662.0) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
        counts.append(max(0, int(round(rng.gauss(mean, math.sqrt(mean))))) if mean > 0 else 0)
    return counts


def write_n42(path, samples, serial='SN-1001', model='TestPortal', manufacturer='TestCo', seed=1):
    """`samples`: list of (class code, ISO start time, real time s, {det: cs137 rate}, neutron counts or None)."""
    rng = random.Random(seed)
    dets = sorted({d for s in samples for d in s[3]})
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<RadInstrumentData xmlns="http://physics.nist.gov/N42/2011/N42" n42DocUUID="%s">'
           % ('d0c0%04d-0000-4000-8000-%012d' % (seed, seed)),
           '<RadInstrumentInformation id="InstInfo1">',
           '<RadInstrumentManufacturerName>%s</RadInstrumentManufacturerName>' % manufacturer,
           '<RadInstrumentIdentifier>%s</RadInstrumentIdentifier>' % serial,
           '<RadInstrumentModelName>%s</RadInstrumentModelName>' % model,
           '<RadInstrumentClassCode>Portal Monitor</RadInstrumentClassCode>',
           '</RadInstrumentInformation>']
    for d in dets:
        out.append('<RadDetectorInformation id="%s"><RadDetectorCategoryCode>Gamma</RadDetectorCategoryCode>'
                   '<RadDetectorKindCode>NaI</RadDetectorKindCode></RadDetectorInformation>' % d)
    if any(s[4] is not None for s in samples):
        out.append('<RadDetectorInformation id="Nn1"><RadDetectorCategoryCode>Neutron</RadDetectorCategoryCode>'
                   '<RadDetectorKindCode>He3</RadDetectorKindCode></RadDetectorInformation>')
    out.append('<EnergyCalibration id="ECal"><CoefficientValues>0 %g 0</CoefficientValues></EnergyCalibration>' % GAIN)
    for i, (cls, start, real, det_rates, neutrons) in enumerate(samples, 1):
        out.append('<RadMeasurement id="Sample%d">' % i)
        out.append('<MeasurementClassCode>%s</MeasurementClassCode>' % cls)
        out.append('<StartDateTime>%s</StartDateTime>' % start)
        out.append('<RealTimeDuration>PT%gS</RealTimeDuration>' % real)
        for d in dets:
            if d not in det_rates:
                continue
            live = real * 0.98
            out.append('<Spectrum id="Sample%d_%s" radDetectorInformationReference="%s" energyCalibrationReference="ECal">'
                       % (i, d, d))
            out.append('<LiveTimeDuration>PT%gS</LiveTimeDuration>' % live)
            out.append('<ChannelData compressionCode="None">%s</ChannelData>'
                       % ' '.join(str(c) for c in spectrum(rng, live, det_rates[d])))
            out.append('</Spectrum>')
        if neutrons is not None:
            out.append('<GrossCounts id="Sample%d_Nn1" radDetectorInformationReference="Nn1">'
                       '<LiveTimeDuration>PT%gS</LiveTimeDuration><CountData>%d</CountData></GrossCounts>'
                       % (i, real, neutrons))
        out.append('</RadMeasurement>')
    out.append('</RadInstrumentData>')
    with open(path, 'w') as f:
        f.write('\n'.join(out))


def make_files(out_dir):
    dets = ['Aa1', 'Aa2', 'Ba1']
    paths = {}

    # One background, then ten 60 s foregrounds; sample 6 has a Cs137 source.
    multi = [('Background', '2024-01-01T11:30:00Z', 300, {d: 0.0 for d in dets}, 40)]
    for k in range(10):
        multi.append(('Foreground', '2024-01-01T12:%02d:00Z' % k, 60, {d: (30.0 if k == 4 else 0.0) for d in dets}, 10))
    paths['multi'] = os.path.join(out_dir, 'multi.n42')
    write_n42(paths['multi'], multi, seed=1)

    # The same, as a search-mode (passthrough) portal file: 4 s samples.
    portal = [multi[0]] + [('Foreground', '2024-01-01T12:00:%02dZ' % (k * 4), 4, r[3], 1) for k, r in enumerate(multi[1:])]
    paths['portal'] = os.path.join(out_dir, 'portal.n42')
    write_n42(paths['portal'], portal, seed=4)

    # Same instrument, but its background is further in time from its foreground than multi.n42's
    #  background - so loading it as the foreground keeps multi.n42's background.
    multi_b = [('Background', '2024-01-01T06:00:00Z', 300, {d: 0.0 for d in dets}, 40),
               ('Foreground', '2024-01-01T12:10:00Z', 20, {d: 5.0 for d in dets}, 3)]
    paths['multi_b'] = os.path.join(out_dir, 'multi_b.n42')
    write_n42(paths['multi_b'], multi_b, seed=2)

    # Two backgrounds: neither is loaded along with the foreground.
    two_backs = [('Foreground', '2024-01-01T12:00:00Z', 60, {d: 0.0 for d in dets}, 10),
                 ('Background', '2024-01-01T11:00:00Z', 300, {d: 0.0 for d in dets}, 40),
                 ('Background', '2024-01-01T13:00:00Z', 300, {d: 0.0 for d in dets}, 40),
                 ('Foreground', '2024-01-01T12:05:00Z', 60, {d: 0.0 for d in dets}, 10)]
    paths['two_backs'] = os.path.join(out_dir, 'two_backs.n42')
    write_n42(paths['two_backs'], two_backs, seed=5)

    # A different instrument's file: a factory background, then the measurement.
    two = [('Background', '2023-06-01T08:00:00Z', 600, {'Det1': 0.0}, None),
           ('Foreground', '2024-01-02T09:00:00Z', 300, {'Det1': 0.0}, None)]
    paths['bg_and_fg'] = os.path.join(out_dir, 'bg_and_fg.n42')
    write_n42(paths['bg_and_fg'], two, serial='OTHER-7', model='Handheld', seed=3)

    # Neutron data only; and a file whose second sample has only neutron data.
    paths['neutron_only'] = os.path.join(out_dir, 'neutron_only.n42')
    write_n42(paths['neutron_only'], [('Foreground', '2024-01-03T09:00:00Z', 60, {}, 25)], seed=6)
    paths['gamma_then_neutron'] = os.path.join(out_dir, 'gamma_then_neutron.n42')
    write_n42(paths['gamma_then_neutron'], [('Foreground', '2024-01-03T10:00:00Z', 60, {'Aa1': 0.0}, 10),
                                            ('Foreground', '2024-01-03T10:01:00Z', 60, {}, 10)], seed=7)
    return paths


class Checker:
    def __init__(self):
        self.failures = 0
        self.passes = 0

    def check(self, cond, label, detail=''):
        if cond:
            self.passes += 1
            print('  ok   ' + label)
        else:
            self.failures += 1
            print('  FAIL ' + label + ('  -- ' + detail if detail else ''))
        return cond


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--url', default='http://127.0.0.1:8078')
    ap.add_argument('--token', default=os.environ.get('INTERSPEC_MCP_BEARER_TOKEN', ''))
    ap.add_argument('--repo', default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
    ap.add_argument('--verbose', action='store_true', help='print every result')
    args = ap.parse_args()

    client = Client(args.url, args.token)
    ck = Checker()
    tmp = tempfile.mkdtemp(prefix='load_test_')
    files = make_files(tmp)
    ex = os.path.join(args.repo, 'example_spectra')
    td = os.path.join(args.repo, 'target', 'testing', 'test_data')
    files.update({'ba133': os.path.join(ex, 'ba133_source_640s_20100317.n42'),
                  'bg': os.path.join(ex, 'background_20100317.n42'),
                  'racetrack': os.path.join(td, 'llm_time_history', 'Q345_Racetrack.n42'),
                  'u': os.path.join(td, 'uranium_40pct_selfatten_HPGe_15cm.n42')})

    sessions = []

    def new_session():
        res, _, err, _ = client.call('create_session', {})
        if err:
            raise RuntimeError('create_session: ' + err)
        sessions.append(res['userSession'])
        return res, {k: upload(client, res, v) for k, v in files.items()}

    def load(handle, **kw):
        res, _, err, _ = client.call('load_spectrum_file', dict(kw, userSession=handle))
        if args.verbose or err:
            print('    load %s -> %s' % (json.dumps(kw), err or json.dumps(res, indent=1)))
        elif res:
            print('    summary: ' + res.get('summary', '(none)'))
            for n in res.get('notices', []):
                print('    notice:  ' + n)
        return res, err

    def disp(res, slot):
        return (res or {}).get('displayed', {}).get(slot, {})

    try:
        sess, p = new_session()
        h = sess['userSession']

        print('\nMulti-sample file as foreground (GUI choice)')
        r, err = load(h, path=p['multi'], specType='Foreground')
        if ck.check(err is None, 'loads', err or ''):
            f = disp(r, 'Foreground')
            ck.check(r['file']['numSamples'] == 11 and not r['file'].get('searchModeOrPassthrough'),
                     'file has 11 samples, not search-mode', json.dumps(r['file']))
            ck.check(r['file']['samplesBySourceType'].get('Foreground', {}).get('count') == 10
                     and r['file']['samplesBySourceType'].get('Background', {}).get('samples') == '1',
                     'file samples by type: 10 foreground, background is 1', json.dumps(r['file']['samplesBySourceType']))
            ck.check(r['file']['gammaDetectors'] == ['Aa1', 'Aa2', 'Ba1'] and r['file'].get('neutronDetectors') == ['Nn1'],
                     'three gamma detectors and a neutron detector', json.dumps(r['file']))
            ck.check(f.get('samples') == '2', 'GUI shows the FIRST foreground sample (2)', f.get('samples', ''))
            ck.check(f.get('gammaDetectors') == ['Aa1', 'Aa2', 'Ba1'] and f.get('neutronDetectors') == ['Nn1'],
                     'sums all detectors', json.dumps(f))
            b = disp(r, 'Background')
            ck.check(b.get('loaded') and b.get('samples') == '1' and b.get('changed') and b.get('fileJustLoaded'),
                     'also loaded its background sample 1 as background', json.dumps(b))
            ck.check(r.get('samplesNotDisplayed', {}).get('Foreground', {}).get('samples') == '3-11'
                     and 'Background' not in r.get('samplesNotDisplayed', {}),
                     'samples 3-11 reported as not displayed', json.dumps(r.get('samplesNotDisplayed')))
            ck.check('Displayed sample 2 (Foreground) as the foreground, summing all 3 gamma detectors' in r['summary']
                     and 'Also loaded sample 1 (Background) of this file as the background' in r['summary'],
                     'summary says what was shown, and the background', r['summary'])
            ck.check(any('first foreground sample is being shown' in n for n in r.get('notices', [])),
                     'GUI notice (now accurate) is passed on', json.dumps(r.get('notices')))
            ck.check(not any('Spectrum Files' in n or 'File Manager' in n for n in r.get('notices', [])),
                     'GUI navigation hints dropped from notices', json.dumps(r.get('notices')))

        print('\nSame file, sample 6 (the source), as foreground')
        r, err = load(h, path=p['multi'], specType='Foreground', sampleNumbers=[6])
        if ck.check(err is None, 'loads', err or ''):
            f = disp(r, 'Foreground')
            ck.check(r.get('reusedLoadedFile') == 'Foreground', 'reuses the displayed copy', json.dumps(r.get('reusedLoadedFile')))
            ck.check(f.get('samples') == '6' and f.get('changed'), 'shows sample 6', json.dumps(f))
            b = disp(r, 'Background')
            ck.check(b.get('samples') == '1' and not b.get('changed'), 'background untouched', json.dumps(b))
            ck.check(f.get('totalCounts', 0) > 0, 'has counts')

        print('\nOnly detector Aa1')
        r, err = load(h, path=p['multi'], specType='Foreground', sampleNumbers=[6], detectors=['Aa1'])
        if ck.check(err is None, 'loads', err or ''):
            f = disp(r, 'Foreground')
            ck.check(f.get('gammaDetectors') == ['Aa1'] and 'neutronDetectors' not in f, 'sums only Aa1', json.dumps(f))
            b = disp(r, 'Background')
            ck.check(b.get('gammaDetectors') == ['Aa1'] and b.get('samples') == '1', 'background follows the detector choice', json.dumps(b))
            ck.check('only gamma detector Aa1' in r['summary'] and 'The background now uses gamma detector Aa1' in r['summary']
                     and 'Also loaded' not in r['summary'], 'summary mentions the detectors', r['summary'])

        print('\nDetectors Aa1 + Ba1 + neutron')
        r, err = load(h, path=p['multi'], specType='Foreground', sampleNumbers=[6, 7], detectors=['Aa1', 'Ba1', 'Nn1'])
        if ck.check(err is None, 'loads', err or ''):
            f = disp(r, 'Foreground')
            ck.check(f.get('samples') == '6,7' and f.get('gammaDetectors') == ['Aa1', 'Ba1']
                     and f.get('neutronDetectors') == ['Nn1'], 'samples 6,7 with Aa1+Ba1+Nn1', json.dumps(f))

        print('\nNo detectors given: back to all of them')
        r, err = load(h, path=p['multi'], specType='Foreground', sampleNumbers=[6])
        ck.check(err is None and disp(r, 'Foreground').get('gammaDetectors') == ['Aa1', 'Aa2', 'Ba1'], 'all detectors',
                 err or json.dumps(disp(r, 'Foreground')))

        print('\nBad arguments')
        _, err = load(h, path=p['multi'], specType='Foreground', detectors=['Zz9'])
        ck.check(err is not None and 'Zz9' in err and 'Aa1' in err, 'unknown detector is an error listing detectors', err or '')
        _, err = load(h, path=p['multi'], specType='Foreground', detectors=['Nn1'])
        ck.check(err is not None and 'gamma' in err, 'neutron-only detector choice is an error', err or '')
        _, err = load(h, path=p['multi'], specType='Background', detectors=['Aa1'])
        ck.check(err is not None and 'Foreground' in err, 'detectors on a background is an error', err or '')
        _, err = load(h, path=p['multi'], specType='Foreground', sampleNumbers=[99])
        ck.check(err is not None and '1-11' in err, 'unknown sample is an error listing the samples', err or '')
        _, err = load(h, path=p['multi'], specType='Foreground', sampleNumbers='3')
        ck.check(err is not None, 'sampleNumbers as a string is an error', err or '')
        r, err = load(h, specType='Secondary')
        ck.check(err is None, 'still responsive', err or '')

        print('\nNo gamma data: an error, and nothing displayed changes')
        _, err = load(h, path=p['neutron_only'], specType='Foreground')
        ck.check(err is not None and 'no gamma spectra' in err, 'neutron-only file is an error', err or '')
        _, err = load(h, path=p['gamma_then_neutron'], specType='Foreground', sampleNumbers=[2])
        ck.check(err is not None and 'no gamma data' in err and 'gamma data are: 1,' in err,
                 'a neutron-only sample is an error, listing the gamma sample', err or '')
        r, err = load(h, specType='Secondary')
        ck.check(err is None and disp(r, 'Foreground').get('fileName') == 'multi.n42',
                 'the foreground is still multi.n42', err or json.dumps(disp(r, 'Foreground')))
        r, err = load(h, path=p['gamma_then_neutron'], specType='Foreground')
        ck.check(err is None and disp(r, 'Foreground').get('samples') == '1', 'its gamma sample 1 loads',
                 err or json.dumps(disp(r, 'Foreground')))
        r, err = load(h, path=p['multi'], specType='Foreground', sampleNumbers=[6])
        ck.check(err is None, 'multi.n42 back as the foreground', err or '')

        print('\nUnload the background, then load the same file as background (GUI choice)')
        r, err = load(h, specType='Background')
        ck.check(err is None and r.get('action') == 'unloaded' and not disp(r, 'Background').get('loaded')
                 and disp(r, 'Foreground').get('samples') == '6', 'unloads, reporting the display', err or json.dumps(r))
        r, err = load(h, path=p['multi'], specType='Background')
        if ck.check(err is None, 'loads', err or ''):
            ck.check(r.get('reusedLoadedFile') == 'Foreground', 'reuses the foreground copy', json.dumps(r.get('reusedLoadedFile')))
            ck.check(disp(r, 'Background').get('samples') == '1', 'same file: uses its Background sample',
                     json.dumps(disp(r, 'Background')))

        print('\nSecondary from the same file')
        r, err = load(h, path=p['multi'], specType='Secondary', sampleNumbers=[7])
        ck.check(err is None and disp(r, 'Secondary').get('samples') == '7' and r.get('reusedLoadedFile') == 'Foreground',
                 'secondary is sample 7 of the same copy', err or json.dumps(r))

        print('\nA background from another instrument (has a factory background, then a foreground)')
        r, err = load(h, path=p['bg_and_fg'], specType='Background')
        if ck.check(err is None, 'loads', err or ''):
            b = disp(r, 'Background')
            ck.check(b.get('samples') == '2' and 'Foreground' in b.get('samplesBySourceType', {}),
                     'different file: uses its Foreground sample as background', json.dumps(b))
            ck.check(any('foreground record' in n for n in r.get('notices', [])), 'GUI notice says so',
                     json.dumps(r.get('notices')))

        print('\nSecond file from the same instrument: its own background is further in time')
        r, err = load(h, path=p['multi'], specType='Background', sampleNumbers=[1])
        ck.check(err is None and disp(r, 'Background').get('samples') == '1', 'background back to multi.n42 sample 1', err or '')
        r, err = load(h, path=p['multi_b'], specType='Foreground')
        if ck.check(err is None, 'loads', err or ''):
            ck.check(disp(r, 'Foreground').get('samples') == '2', 'shows its foreground sample 2', json.dumps(disp(r, 'Foreground')))
            b = disp(r, 'Background')
            ck.check(b.get('fileName') == 'multi.n42' and not b.get('changed') and not b.get('fileJustLoaded'),
                     'keeps the closer-in-time background', json.dumps(b))
            ck.check('was not loaded, keeping the current background' in r['summary'], 'summary says why', r['summary'])
            ck.check(r.get('samplesNotDisplayed', {}).get('Background', {}).get('samples') == '1',
                     'its background sample 1 reported not displayed', json.dumps(r.get('samplesNotDisplayed')))
            ck.check(disp(r, 'Secondary').get('fileName') == 'multi.n42', 'secondary kept (same channel count)',
                     json.dumps(disp(r, 'Secondary')))

        print('\nTwo background samples: ambiguous, so neither is loaded with the foreground')
        r, err = load(h, specType='Background')
        r, err = load(h, path=p['two_backs'], specType='Foreground')
        if ck.check(err is None, 'loads', err or ''):
            ck.check(disp(r, 'Foreground').get('samples') == '1', 'shows the first foreground sample 1',
                     json.dumps(disp(r, 'Foreground')))
            ck.check(not disp(r, 'Background').get('loaded'), 'no background loaded', json.dumps(disp(r, 'Background')))
            ck.check(r.get('samplesNotDisplayed', {}).get('Background', {}).get('samples') == '2,3',
                     'backgrounds 2,3 reported not displayed', json.dumps(r.get('samplesNotDisplayed')))

        print('\nSearch-mode portal file (sums the foreground samples, loads its background)')
        r, err = load(h, path=p['portal'], specType='Foreground')
        if ck.check(err is None, 'loads', err or ''):
            ck.check(r['file'].get('searchModeOrPassthrough') is True, 'reported as search-mode', json.dumps(r['file']))
            ck.check(disp(r, 'Foreground').get('samples') == '2-11', 'sums samples 2-11', json.dumps(disp(r, 'Foreground')))
            ck.check(disp(r, 'Background').get('samples') == '1' and disp(r, 'Background').get('fileJustLoaded'),
                     'its background sample 1 loaded', json.dumps(disp(r, 'Background')))

        print('\nOther instruments')
        r, err = load(h, path=p['u'], specType='Foreground')
        if ck.check(err is None, 'uranium HPGe loads', err or ''):
            ck.check(disp(r, 'Foreground').get('numChannels') == 8192 and disp(r, 'Background').get('samples') == '2'
                     and disp(r, 'Background').get('fileJustLoaded'), 'with its own background sample',
                     json.dumps(r.get('displayed')))
            ck.check(not disp(r, 'Secondary').get('loaded') and disp(r, 'Secondary').get('changed'),
                     'the 1024-channel secondary was unloaded', json.dumps(disp(r, 'Secondary')))
        r, err = load(h, path=p['ba133'], specType='Foreground')
        if ck.check(err is None, 'ba133 loads', err or ''):
            ck.check(not disp(r, 'Background').get('loaded') and 'was unloaded, as its number of channels' in r['summary'],
                     'the 8192-channel background was unloaded', r['summary'])
            ck.check(r['file']['numSamples'] == 1 and not r.get('samplesNotDisplayed'), 'single sample file')
        r, err = load(h, path=p['bg'], specType='Background')
        if ck.check(err is None, 'background loads', err or ''):
            ck.check(disp(r, 'Background').get('fileName') == 'background_20100317.n42'
                     and not disp(r, 'Foreground').get('changed'), 'foreground unchanged, background set',
                     json.dumps(r.get('displayed')))
        r, err = load(h, path=p['racetrack'], specType='Foreground')
        if ck.check(err is None, 'racetrack (search-mode) loads', err or ''):
            ck.check(r['file'].get('searchModeOrPassthrough') is True and disp(r, 'Foreground').get('samples'),
                     'sums: ' + disp(r, 'Foreground').get('samples', ''), json.dumps(disp(r, 'Foreground')))

        print('\nBackground before foreground, with the same number of channels: the background is kept')
        sess3, p3 = new_session()
        h3 = sess3['userSession']
        r, err = load(h3, path=p3['bg'], specType='Background')
        ck.check(err is None and 'No foreground is loaded' in r.get('summary', ''), 'background alone', err or json.dumps(r))
        r, err = load(h3, path=p3['ba133'], specType='Foreground')
        ck.check(err is None and disp(r, 'Background').get('fileName') == 'background_20100317.n42'
                 and not disp(r, 'Background').get('changed'), 'background kept', err or json.dumps(r.get('displayed')))

        print('\nExported N42 remembers the displayed samples; sampleNumbers overrides that')
        load(h, path=p['multi'], specType='Foreground', sampleNumbers=[6])
        res, _, err, _ = client.call('export_analysis', {'userSession': h, 'items': ['n42'], 'basename': 'multi_s6'})
        if ck.check(err is None, 'exported', err or ''):
            sess2, _ = new_session()
            # Download it, and upload it to the second session
            dl = [x['download_url'] for x in res['files'] if x['path'].endswith('.n42')][0]
            status, body = client.http('GET', dl)
            local = os.path.join(tmp, 'multi_s6.n42')
            with open(local, 'wb') as f:
                f.write(body)
            p2 = upload(client, sess2, local)
            r, err = load(sess2['userSession'], path=p2, specType='Foreground')
            ck.check(err is None and disp(r, 'Foreground').get('samples') == '6'
                     and r['file']['samplesBySourceType'].get('Background', {}).get('samples') == '1',
                     'reload shows the saved sample 6, same numbering', err or json.dumps(r))
            r, err = load(sess2['userSession'], specType='Foreground')
            r, err = load(sess2['userSession'], path=p2, specType='Foreground', sampleNumbers=[3])
            ck.check(err is None and disp(r, 'Foreground').get('samples') == '3', 'sampleNumbers [3] wins over saved',
                     err or json.dumps(disp(r, 'Foreground')))
    finally:
        for s in sessions:
            client.call('close_session', {'userSession': s})

    print('\n%d passed, %d failed' % (ck.passes, ck.failures))
    sys.exit(1 if ck.failures else 0)


if __name__ == '__main__':
    main()
