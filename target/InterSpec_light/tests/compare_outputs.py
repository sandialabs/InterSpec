#!/usr/bin/env python3
"""Compares native and WASM response logs: same structure, numbers equal to a relative 2e-5
(summing/rebinning is done in float32, whose rounding differs slightly between the builds).
File paths differ (the Node runner loads files from its in-memory filesystem), so "path" and
"name" strings are not compared."""
import sys
import json

def same(a, b, where):
  if isinstance(a, dict) and isinstance(b, dict):
    if set(a) != set(b):
      return f'{where}: keys differ {sorted(set(a) ^ set(b))}'
    for k in a:
      if k in ('path',):
        continue
      r = same(a[k], b[k], where + '.' + k)
      if r:
        return r
    return None
  if isinstance(a, list) and isinstance(b, list):
    if len(a) != len(b):
      return f'{where}: lengths {len(a)} vs {len(b)}'
    for i, (x, y) in enumerate(zip(a, b)):
      r = same(x, y, f'{where}[{i}]')
      if r:
        return r
    return None
  if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
    if abs(a - b) > 2e-5 * max(1.0, abs(a), abs(b)):
      return f'{where}: {a} vs {b}'
    return None
  return None if a == b else f'{where}: {str(a)[:80]!r} vs {str(b)[:80]!r}'

native = [json.loads(l) for l in open(sys.argv[1])]
wasm = [json.loads(l) for l in open(sys.argv[2])]
if len(native) != len(wasm):
  sys.exit(f'{len(native)} vs {len(wasm)} responses')
for i, (a, b) in enumerate(zip(native, wasm)):
  r = same(a, b, f'response[{i}]')
  if r:
    sys.exit(r)
