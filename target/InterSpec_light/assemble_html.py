#!/usr/bin/env python3
"""Assembles the self-contained InterSpec Light pages from the page template, the chart libraries
(read from the InterSpec repo), the app JS/CSS, the reference-line library, and the Emscripten
module (light_wasm.js + light_wasm.wasm).  Two variants are written:

  InterSpecLight.html               the WASM module and reference lines are gzip'ed then base64'ed
                                    (~2.2 MB); decompressed on load with the browser's
                                    DecompressionStream (Safari 16.4+, Chrome 80+, Firefox 113+).
  InterSpecLight_uncompressed.html  the WASM module is plain base64, and the reference lines plain
                                    JSON (~6 MB), for browsers or content scanners that object to
                                    the compressed data.

In both, the (minified) JavaScript is inline and readable; only data is encoded.

Usage: assemble_html.py --wasm-js <light_wasm.js> --ref-lines <ref_lines.json> --output-dir <dir>
       (light_wasm.wasm is expected beside light_wasm.js)
"""
import io
import os
import re
import sys
import gzip
import base64
import shutil
import argparse
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
WEB = os.path.join(HERE, 'web')

# Order matters: later scripts use earlier ones
APP_JS = [ 'core.js', 'charts.js', 'files.js', 'reflines.js', 'search.js', 'peaks.js', 'energycal.js', 'main.js' ]

def read(path):
  with open(path, 'r', encoding='utf-8') as f:
    return f.read()


def minify(js, name):
  """Minifies with the terser that ships with Emscripten (local names only; globals, and function
  arities - which SpectrumChartD3's WtEmit relies on - are kept).  Returns `js` unchanged if terser
  or node are not found (an activated emsdk sets EMSDK and EMSDK_NODE)."""
  emsdk = os.environ.get('EMSDK', '')
  terser = os.path.join(emsdk, 'upstream', 'emscripten', 'node_modules', 'terser', 'bin', 'terser')
  node = os.environ.get('EMSDK_NODE') or shutil.which('node')
  if not emsdk or not node or not os.path.exists(terser):
    print(f'  Warning: not minifying {name} (node or terser not found)')
    return js
  result = subprocess.run([node, terser, '--compress', '--mangle'], input=js.encode('utf-8'),
                          capture_output=True, check=True)
  return result.stdout.decode('utf-8')


def b64(data, compress):
  """base64, optionally of the gzip'ed data (max compression, fixed mtime for reproducible output)."""
  if isinstance(data, str):
    data = data.encode('utf-8')
  if compress:
    data = gzip.compress(data, compresslevel=9, mtime=0)
  return base64.b64encode(data).decode('ascii')


def png_data_uri(path, height):
  """The image as a PNG data URI, scaled to `height` pixels tall (2x its display size) if PIL is available."""
  try:
    from PIL import Image
    img = Image.open(path).convert('RGBA')
    width = max(1, round(img.width * height / img.height))
    buf = io.BytesIO()
    small = img.resize((width, height), Image.LANCZOS).quantize(colors=64, method=Image.Quantize.FASTOCTREE)
    small.save(buf, format='PNG', optimize=True)
    data = buf.getvalue()
  except ImportError:
    with open(path, 'rb') as f:
      data = f.read()
  return 'data:image/png;base64,' + base64.b64encode(data).decode('ascii')


def check_inline(name, content, closing):
  """An inlined <script>/<style> must not contain its closing tag, or an HTML comment opener."""
  bad = re.search(r'(?i)</' + closing + r'|<!--', content)
  if bad:
    sys.exit(f'ERROR: {name} contains "{bad.group(0)}", which would break the inlined <{closing}>')


def inline(html, marker, tag, name, content, attrs=''):
  check_inline(name, content, tag)
  token = f'<!-- {marker} -->'
  if html.count(token) != 1:
    sys.exit(f'ERROR: expected one "{token}" in the template')
  return html.replace(token, f'<{tag}{attrs}>\n{content}\n</{tag}>')


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--wasm-js', required=True)
  parser.add_argument('--ref-lines', required=True)
  parser.add_argument('--output-dir', required=True)
  args = parser.parse_args()

  d3dir = os.path.join(REPO, 'external_libs', 'SpecUtils', 'd3_resources')
  resdir = os.path.join(REPO, 'InterSpec_resources')
  wasm_path = os.path.splitext(args.wasm_js)[0] + '.wasm'

  html = read(os.path.join(WEB, 'index.html'))
  images = os.path.join(resdir, 'images')
  html = html.replace('%%INTERSPEC_LOGO%%', png_data_uri(os.path.join(images, 'InterSpec128.png'), 44))
  html = html.replace('%%SANDIA_LOGO%%', png_data_uri(os.path.join(images, 'SNL_Stacked_Black_Blue.png'), 40))
  if '%%' in html:
    sys.exit('ERROR: unreplaced %% token in the template')
  html = inline(html, 'CHART_CSS', 'style', 'SpectrumChartD3.css', read(os.path.join(d3dir, 'SpectrumChartD3.css')))
  html = inline(html, 'TIME_CSS', 'style', 'D3TimeChart.css', read(os.path.join(resdir, 'D3TimeChart.css')))
  html = inline(html, 'APP_CSS', 'style', 'style.css', read(os.path.join(WEB, 'style.css')))

  # The scripts, in the order they must run, inlined as ordinary <script> elements
  scripts = [ ('d3.v3.min.js', read(os.path.join(d3dir, 'd3.v3.min.js'))),
              ('SpectrumChartD3.js', minify(read(os.path.join(d3dir, 'SpectrumChartD3.js')), 'SpectrumChartD3.js')),
              ('D3TimeChart.js', minify(read(os.path.join(resdir, 'D3TimeChart.js')), 'D3TimeChart.js')),
              (os.path.basename(args.wasm_js), minify(read(args.wasm_js), os.path.basename(args.wasm_js))),
              ('app js', minify('\n\n'.join(read(os.path.join(WEB, name)) for name in APP_JS), 'app js')) ]
  script_html = ''
  for name, js in scripts:
    check_inline(name, js, 'script')
    script_html += f'<script>\n{js}\n</script>\n'
  if html.count('<!-- SCRIPTS -->') != 1 or html.count('<!-- ASSETS -->') != 1:
    sys.exit('ERROR: expected one "<!-- SCRIPTS -->" and one "<!-- ASSETS -->" in the template')
  html = html.replace('<!-- SCRIPTS -->', script_html)

  with open(wasm_path, 'rb') as f:
    wasm = f.read()
  ref_lines = read(args.ref_lines)
  # JSON can hold "</" only inside strings, where "<\/" is an equivalent escape
  ref_lines_escaped = ref_lines.replace('</', '<\\/')

  os.makedirs(args.output_dir, exist_ok=True)
  for filename, compress in [ ('InterSpecLight.html', True), ('InterSpecLight_uncompressed.html', False) ]:
    # Data blocks read by App.readAsset (web/main.js); they are not scripts, so never run
    encoding = 'gzip-base64' if compress else 'base64'
    wasm_block = f'<script type="application/octet-stream" id="light-wasm" data-encoding="{encoding}">{b64(wasm, compress)}</script>'
    if compress:
      ref_block = f'<script type="application/octet-stream" id="light-ref-lines" data-encoding="gzip-base64">{b64(ref_lines, True)}</script>'
    else:
      ref_block = f'<script type="application/json" id="light-ref-lines" data-encoding="text">{ref_lines_escaped}</script>'

    page = html.replace('<!-- ASSETS -->', wasm_block + '\n' + ref_block)
    out_path = os.path.join(args.output_dir, filename)
    with open(out_path, 'w', encoding='utf-8') as f:
      f.write(page)
    print(f'Wrote {out_path} ({len(page.encode("utf-8")):,} bytes)')


if __name__ == '__main__':
  main()
