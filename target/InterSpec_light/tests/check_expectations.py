#!/usr/bin/env python3
"""Checks the expectations a JSONL request script puts on its responses (extra request keys, which
the app ignores):
  "save": "<label>"         remember this response's peak list
  "expectPeaks": "<label>"  the peak list must match the saved one; "match" is "exact" (N42: all
                            fields, to float precision) or "csv" (SPE: the peak CSV's precision, and
                            without the fields a CSV does not hold)
  "expectPeakCount": N      the peak list has N peaks
  "expectSources": [...]    the peak list's source labels, in order ("" for no source)
  "expectColors": [...]     the peak list's colors, in order ("" for the default)
  "expectCandidates": [...] a peakInfoAt response suggests these sources, in order
  "expectSamples": [...]    the foreground shows these samples
  "expectMessage": "text"   a message contains this text
  "expectError": "text"     the request fails, with an error containing this text (the runners
                            do not count it as a failure)
Usage: check_expectations.py <script.jsonl> <responses> [first response index]"""
import sys
import json

def main():
  requests = [json.loads(l) for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
  responses = [json.loads(l) for l in open(sys.argv[2])]
  offset = int(sys.argv[3]) if len(sys.argv) > 3 else 1  # the reference-line library response
  if len(responses) != len(requests) + offset:
    sys.exit(f'{len(responses)} responses for {len(requests)} requests')

  saved, errors = {}, []
  for i, req in enumerate(requests):
    resp = responses[i + offset]
    where = f'request {i + 1} ({req["method"]})'
    peaks = resp.get('peakList')
    if 'save' in req:
      saved[req['save']] = peaks or []
    if 'error' in resp and 'expectError' not in req:
      errors.append(f'{where}: error {resp["error"]!r}')
    if 'expectError' in req and req['expectError'] not in resp.get('error', ''):
      errors.append(f'{where}: no error containing {req["expectError"]!r} (got {resp.get("error")!r})')
    if 'expectPeakCount' in req and len(peaks or []) != req['expectPeakCount']:
      errors.append(f'{where}: {len(peaks or [])} peaks, expected {req["expectPeakCount"]}')
    if 'expectSources' in req and [p['source'] for p in peaks or []] != req['expectSources']:
      errors.append(f'{where}: sources {[p["source"] for p in peaks or []]}, expected {req["expectSources"]}')
    if 'expectColors' in req and [p['color'] for p in peaks or []] != req['expectColors']:
      errors.append(f'{where}: colors {[p["color"] for p in peaks or []]}, expected {req["expectColors"]}')
    if 'expectCandidates' in req and ((resp.get('peak') or {}).get('candidates')) != req['expectCandidates']:
      errors.append(f'{where}: candidates {(resp.get("peak") or {}).get("candidates")}, expected {req["expectCandidates"]}')
    if 'expectSamples' in req:
      samples = ((resp.get('files') or {}).get('FOREGROUND') or {}).get('samples')
      if samples != req['expectSamples']:
        errors.append(f'{where}: samples {samples}, expected {req["expectSamples"]}')
    msgs = resp.get('messages', []) + ([resp['message']] if 'message' in resp else [])
    if 'expectMessage' in req and not any(req['expectMessage'] in m for m in msgs):
      errors.append(f'{where}: no message containing {req["expectMessage"]!r} in {msgs}')
    if 'expectPeaks' in req:
      errors += compare(saved[req['expectPeaks']], peaks or [], req.get('match', 'exact'), where)

  for e in errors:
    print('  ' + e)
  sys.exit(1 if errors else 0)


def compare(expected, actual, match, where):
  if len(expected) != len(actual):
    return [f'{where}: {len(actual)} peaks, expected {len(expected)}']
  csv = (match == 'csv')
  # (field, absolute tolerance, relative tolerance); CSV values are printed to ~0.01 keV / 0.1 counts
  numeric = [('mean', 0.006, 1e-6), ('fwhm', 0.006, 1e-6), ('area', 0.051, 1e-6), ('lower', 1e-3, 1e-5),
             ('upper', 1e-3, 1e-5)] if csv else \
            [('mean', 0, 2e-6), ('meanUnc', 1e-9, 2e-6), ('fwhm', 0, 2e-6), ('area', 1e-6, 2e-6),
             ('areaUnc', 1e-6, 2e-6), ('lower', 0, 2e-6), ('upper', 0, 2e-6), ('chi2dof', 1e-6, 2e-6)]
  same = ['source', 'continuumType', 'skewType', 'gaussian', 'color'] + ([] if csv else ['useForCal'])
  errors = []
  for n, (a, b) in enumerate(zip(expected, actual)):
    for key, atol, rtol in numeric:
      if abs(a[key] - b[key]) > max(atol, rtol * max(abs(a[key]), abs(b[key]))):
        errors.append(f'{where}: peak {n} ({a["mean"]:.2f} keV) {key} {b[key]} != {a[key]}')
    for key in same:
      if a[key] != b[key]:
        errors.append(f'{where}: peak {n} ({a["mean"]:.2f} keV) {key} {b[key]!r} != {a[key]!r}')
  return errors


if __name__ == '__main__':
  main()
