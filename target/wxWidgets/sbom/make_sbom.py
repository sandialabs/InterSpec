#!/usr/bin/env python3
"""Generate a CycloneDX 1.6 SBOM for the InterSpec Windows (wxWidgets) build.

The list of components is curated in ``components.json`` (single source of truth), because
InterSpec has no package manifest and statically links everything into the final binary, so
neither a binary nor a source-tree scanner can produce an accurate, low-noise SBOM. See
README.md for the model and the deliberate exclusion list.

Two modes:
  * generate (default): resolve submodule commit SHAs from the checked-out tree, stamp the
    InterSpec version + git SHA, and write a CycloneDX 1.6 JSON SBOM to ``--output``.
  * --check (drift guard): verify components.json still matches the authoritative sources
    (dep_build_msvc2022.bat pinned versions, submodule paths, vendored dirs). Exits non-zero
    on any mismatch so a dependency bump forces an SBOM update.

Pure standard library; runs on the windows-2022 runner's Python 3 and locally on any OS.
"""

import argparse
import json
import os
import re
import subprocess
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COMPONENTS_JSON = os.path.join(SCRIPT_DIR, "components.json")
# Authoritative pinned versions for the prebuilt Windows dependencies.
DEP_BAT_REL = os.path.join("target", "patches", "dep_build_msvc2022.bat")

INTERSPEC_REPO = "https://github.com/sandialabs/InterSpec"


def load_components():
  with open(COMPONENTS_JSON, "r", encoding="utf-8") as f:
    return json.load(f)["components"]


def read_interspec_version(repo):
  """Read the project version from the top-level CMakeLists.txt (project(InterSpec VERSION x))."""
  cmake = os.path.join(repo, "CMakeLists.txt")
  try:
    with open(cmake, "r", encoding="utf-8") as f:
      text = f.read()
    m = re.search(r"project\s*\(\s*InterSpec\s+VERSION\s+([0-9][0-9.]*)", text)
    if m:
      return m.group(1)
  except OSError:
    pass
  return None


def git_sha(repo, path):
  """Return the HEAD commit SHA of the (sub)module checked out at repo/path, or None."""
  try:
    out = subprocess.run(
        ["git", "-C", os.path.join(repo, path), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=True)
    return out.stdout.strip()
  except (subprocess.CalledProcessError, OSError):
    return None


def external_references(comp):
  refs = []
  if comp.get("vcs"):
    refs.append({"type": "vcs", "url": comp["vcs"]})
  if comp.get("url"):
    refs.append({"type": "website", "url": comp["url"]})
  return refs


def to_cyclonedx_component(comp, repo):
  """Convert one curated entry into a CycloneDX component dict, resolving submodule SHAs."""
  src = comp["version_source"]
  version = comp.get("version")
  purl = comp.get("purl")

  if src.startswith("submodule:"):
    path = src.split(":", 1)[1]
    sha = git_sha(repo, path)
    if sha:
      version = sha
      base = comp.get("purl_base")
      if base:
        purl = "{0}@{1}".format(base, sha)

  out = {
      "type": comp.get("type", "library"),
      "name": comp["name"],
      "bom-ref": comp["name"] if not comp.get("group") else comp["group"] + "/" + comp["name"],
  }
  if comp.get("group"):
    out["group"] = comp["group"]
  if version:
    out["version"] = version
  if purl:
    out["purl"] = purl
  if comp.get("license"):
    out["licenses"] = [{"license": {"id": comp["license"]}}]
  refs = external_references(comp)
  if refs:
    out["externalReferences"] = refs
  if comp.get("notes"):
    out["description"] = comp["notes"]
  return out


def build_sbom(repo, git_sha_value, timestamp):
  components = load_components()
  version = read_interspec_version(repo) or (git_sha_value or "unknown")

  metadata_component = {
      "type": "application",
      "name": "InterSpec",
      "bom-ref": "InterSpec",
      "version": version,
      "externalReferences": [{"type": "vcs", "url": INTERSPEC_REPO}],
      "description": "Nuclear radiation spectral analysis application (Windows wxWidgets build).",
  }
  if git_sha_value:
    metadata_component["externalReferences"].append(
        {"type": "vcs", "url": "{0}/tree/{1}".format(INTERSPEC_REPO, git_sha_value),
         "comment": "Exact source commit for this build"})

  metadata = {
      "tools": [{"vendor": "InterSpec", "name": "make_sbom.py"}],
      "component": metadata_component,
  }
  if timestamp:
    metadata["timestamp"] = timestamp

  return {
      "bomFormat": "CycloneDX",
      "specVersion": "1.6",
      "version": 1,
      "metadata": metadata,
      "components": [to_cyclonedx_component(c, repo) for c in components],
  }


def check_drift(repo):
  """Validate components.json against authoritative sources. Returns list of problems."""
  problems = []
  components = load_components()

  # 1) Every bat: token must appear literally in dep_build_msvc2022.bat.
  bat_path = os.path.join(repo, DEP_BAT_REL)
  try:
    with open(bat_path, "r", encoding="utf-8", errors="replace") as f:
      bat_text = f.read()
  except OSError as e:
    return ["Cannot read {0}: {1}".format(DEP_BAT_REL, e)]

  for comp in components:
    src = comp["version_source"]
    kind, _, arg = src.partition(":")
    if kind == "bat":
      if arg not in bat_text:
        problems.append(
            "Component '{0}': token not found in {1} -> \"{2}\" "
            "(a dependency version likely changed; update components.json)."
            .format(comp["name"], DEP_BAT_REL, arg))
    elif kind in ("submodule", "vendored", "runtime"):
      full = os.path.join(repo, arg)
      if not os.path.exists(full):
        problems.append(
            "Component '{0}': {1} path does not exist -> {2}".format(comp["name"], kind, arg))
      if kind == "submodule" and git_sha(repo, arg) is None:
        problems.append(
            "Component '{0}': submodule not checked out (no git SHA) -> {1} "
            "(run with submodules initialized).".format(comp["name"], arg))
    else:
      problems.append(
          "Component '{0}': unknown version_source kind '{1}'".format(comp["name"], kind))

  return problems


def main(argv):
  ap = argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("--repo", default=os.getcwd(),
                  help="Path to the InterSpec repo checkout (default: cwd).")
  ap.add_argument("--check", action="store_true",
                  help="Drift-guard mode: validate components.json against sources and exit.")
  ap.add_argument("--git-sha", default=None, help="Build's git commit SHA (for metadata).")
  ap.add_argument("--timestamp", default=None,
                  help="ISO-8601 UTC timestamp for metadata.timestamp (pass one for reproducibility).")
  ap.add_argument("--output", default=None,
                  help="Path to write the SBOM JSON (default: stdout).")
  args = ap.parse_args(argv)

  repo = os.path.abspath(args.repo)

  if args.check:
    problems = check_drift(repo)
    if problems:
      sys.stderr.write("SBOM component list is out of sync with the build sources:\n")
      for p in problems:
        sys.stderr.write("  - {0}\n".format(p))
      sys.stderr.write(
          "\nUpdate target/wxWidgets/sbom/components.json to match, then re-run.\n")
      return 1
    sys.stdout.write("SBOM component list OK ({0} components).\n".format(len(load_components())))
    return 0

  sbom = build_sbom(repo, args.git_sha, args.timestamp)
  text = json.dumps(sbom, indent=2)
  if args.output:
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
      f.write(text + "\n")
    sys.stderr.write("Wrote SBOM with {0} components to {1}\n"
                     .format(len(sbom["components"]), args.output))
  else:
    sys.stdout.write(text + "\n")
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
