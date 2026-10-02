#!/usr/bin/env python3
"""Seeds refgen/sources.txt with the nuclides, reactions, and x-ray elements that InterSpec's
data/ directory refers to.  The output is meant to be edited by hand afterwards; make_ref_lines
canonicalizes/validates names, and skips sources without photon lines.

Usage: collect_sources.py <InterSpec/data dir> <output sources.txt>
"""
import os
import re
import sys
import glob
import xml.etree.ElementTree as ET

# Common shielding / detector / sample elements, for fluorescence x-rays
DEFAULT_XRAY_ELEMENTS = [ 'Pb', 'W', 'Ta', 'Bi', 'U', 'Pu', 'Th', 'Np', 'Ra', 'Rn', 'Cd', 'Sn', 'In', 'Sb',
                          'Ba', 'Cs', 'I', 'Te', 'Fe', 'Cu', 'Zn', 'Ni', 'Ag', 'Au', 'Pt', 'Hg', 'Tl', 'Zr',
                          'Mo', 'Ge', 'Gd', 'Ga' ]

SKIPPED_CATEGORIES = { 'isbe-category-alphas', 'isbe-category-betas' }


def local(tag):
  return tag.split('}')[-1]


def main():
  if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)
  data_dir, out_path = sys.argv[1], sys.argv[2]

  nuclides, reactions, elements = [], [], list(DEFAULT_XRAY_ELEMENTS)

  def add(lst, name):
    name = name.strip()
    if name and name.lower() not in (n.lower() for n in lst):
      lst.append(name)

  # Reference-line definitions
  for fname in ('dynamic_ref_lines.xml', 'add_ref_line.xml'):
    root = ET.parse(os.path.join(data_dir, fname)).getroot()
    for el in root.iter():
      if local(el.tag) in ('IndividualSource', 'Nuc') and el.get('name'):
        add(nuclides, el.get('name'))

  # Nuclide search categories (not the alpha/beta-only ones)
  root = ET.parse(os.path.join(data_dir, 'NuclideSearchCatagories.xml')).getroot()
  for cat in root.iter('NucSearchCategory'):
    name = (cat.findtext('Name') or '').strip()
    if name in SKIPPED_CATEGORIES:
      continue
    for el in cat.iter():
      tag = local(el.tag)
      if tag == 'Nuclide' and el.text:
        add(nuclides, el.text)
      elif tag == 'Reaction' and el.text:
        add(reactions, el.text)
      elif tag == 'Element' and el.text:
        add(elements, el.text)

  # Nuclides with extra notes (FRMAC-derived list)
  root = ET.parse(os.path.join(data_dir, 'more_nuclide_info.xml')).getroot()
  for el in root.iter('Nuc'):
    if el.get('name'):
      add(nuclides, el.get('name'))

  # Characteristic gammas, and isotopics presets
  with open(os.path.join(data_dir, 'CharacteristicGammas.txt')) as f:
    for line in f:
      parts = line.split()
      if parts and not parts[0].startswith('#'):
        add(nuclides, parts[0])
  for fname in glob.glob(os.path.join(data_dir, 'rel_act', '*.xml')):
    for el in ET.parse(fname).getroot().iter('Nuclide'):
      if el.text:
        add(nuclides, el.text)

  # Reference spectra file names, e.g. "Ba133_Shielded.txt"
  for fname in glob.glob(os.path.join(data_dir, 'reference_spectra', '*', '*', '*.txt')):
    m = re.match(r'([A-Z][a-z]?-?\d{1,3}m?\d?)_', os.path.basename(fname))
    if m:
      add(nuclides, m.group(1))

  with open(out_path, 'w') as out:
    out.write('# Sources for make_ref_lines, seeded by collect_sources.py from InterSpec/data - edit freely.\n')
    out.write('# Format: <kind> <name>, where kind is nuclide, reaction, xray, or background.\n')
    out.write('background Background\n')
    for n in nuclides:
      out.write(f'nuclide {n}\n')
    for r in reactions:
      out.write(f'reaction {r}\n')
    for e in elements:
      out.write(f'xray {e}\n')

  print(f'Wrote {len(nuclides)} nuclides, {len(reactions)} reactions, {len(elements)} elements to {out_path}')


if __name__ == '__main__':
  main()
