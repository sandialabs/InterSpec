#!/usr/bin/env python3
"""wt_resource_audit.py -- decide which Wt resources/ files InterSpec actually needs.

InterSpec bundles the whole Wt framework ``resources/`` directory (~5.9 MB) into
every packaged app, but uses only a slice of it.  This script is the *verification
step* that gates the configure-time prune (``cmake/PruneWtResources.cmake``): for
every file under Wt's ``resources/`` it works out, with evidence, whether the file
can be removed -- by tracing *how Wt would load it* and whether InterSpec ever
triggers that code path.

It is deliberately conservative: a file is only marked REMOVE when an app-wide
condition makes its loader unreachable (an unused CSS theme, IE-only stylesheets,
RTL-only images, or font-awesome which only ``Wt::WIcon`` pulls in).  Anything it
cannot positively justify removing is KEEP or REVIEW (both kept by the prune).

Outputs / modes
---------------
  (default)        print a per-file report and write the JSON manifest +
                   the CMake prune list.
  --check          recompute and compare against the committed manifest; exit
                   non-zero on drift (new/renamed/removed files, or a verdict
                   that flipped).  This is the staleness gate for Wt upgrades.
  --no-write       compute + report but do not write any file.

The loader facts (Wt file:line) are re-verified live against the Wt source on
every run, so a Wt upgrade that moves a loader is reflected automatically.
"""

import argparse
import json
import os
import re
import sys
import glob


# ---------------------------------------------------------------------------
# Locating inputs
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# target/tools/wt_resource_audit -> repo root
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(SCRIPT_DIR)))


def _newest_wt_src():
    """Best-effort: newest /…/install/wt-<ver>_src_code with resources/ + src/Wt."""
    cands = sorted(glob.glob(os.path.expanduser("~/install/wt-*_src_code")))
    for d in reversed(cands):
        if os.path.isdir(os.path.join(d, "resources")) and \
           os.path.isdir(os.path.join(d, "src", "Wt")):
            return d
    return None


def resolve_inputs(args):
    wt_res = args.wt_resources
    wt_src = args.wt_source
    # FetchContent layout: <wt-src>/resources sits beside <wt-src>/src/Wt.
    if not wt_src and wt_res:
        cand = os.path.dirname(os.path.abspath(wt_res))
        if os.path.isdir(os.path.join(cand, "src", "Wt")):
            wt_src = cand
    if not wt_src:
        wt_src = _newest_wt_src()  # e.g. an install-prefix build's sibling _src_code
    if not wt_res and wt_src:
        wt_res = os.path.join(wt_src, "resources")
    wt_wt = None
    if wt_src and os.path.isdir(os.path.join(wt_src, "src", "Wt")):
        wt_wt = os.path.join(wt_src, "src", "Wt")
    return wt_res, wt_wt


# ---------------------------------------------------------------------------
# Small grep helpers (read InterSpec / Wt source as text)
# ---------------------------------------------------------------------------

def _iter_source_files(root, exts):
    for base, _dirs, files in os.walk(root):
        # skip build / vcs noise
        if any(p in base for p in (os.sep + "build", os.sep + ".git", "_deps")):
            continue
        for f in files:
            if f.endswith(exts):
                yield os.path.join(base, f)


def grep(root, pattern, exts=(".cpp", ".h", ".hpp", ".C", ".cc")):
    """Return list of (relpath, lineno, line) matching the compiled regex."""
    rx = re.compile(pattern)
    hits = []
    if not root or not os.path.isdir(root):
        return hits
    for path in _iter_source_files(root, exts):
        try:
            with open(path, "r", errors="ignore") as fh:
                for n, line in enumerate(fh, 1):
                    if rx.search(line):
                        hits.append((os.path.relpath(path, root), n, line.rstrip()))
        except OSError:
            pass
    return hits


def first_hit(hits):
    return "{}:{}".format(hits[0][0], hits[0][1]) if hits else None


# ---------------------------------------------------------------------------
# InterSpec usage probes (the "does InterSpec trigger this Wt path" questions)
# ---------------------------------------------------------------------------

def probe_interspec(interspec_root):
    src_dirs = [os.path.join(interspec_root, d) for d in ("src", "InterSpec")]

    def g(pattern):
        out = []
        for d in src_dirs:
            out += grep(d, pattern)
        return out

    themes = sorted({m for (_f, _n, line) in g(r'setCssTheme\(')
                     for m in re.findall(r'setCssTheme\(\s*"([^"]+)"', line)})
    # bare Wt::WIcon (font-awesome glyph widget) -- exclude WIconPair
    wicon_hits = [h for h in g(r'\bWIcon\b') if "WIconPair" not in h[2]]
    fa_name_hits = g(r'(setIcon|addStyleClass|styleClass)\([^)]*\bfa[- ]')
    media_hits = g(r'\bWMediaPlayer\b')
    rtl_hits = g(r'RightToLeft|setLayoutDirection')

    return {
        "active_themes": themes or ["default"],
        "uses_wicon": bool(wicon_hits),
        "wicon_evidence": first_hit(wicon_hits),
        "uses_fa_names": bool(fa_name_hits),
        "uses_media_player": bool(media_hits),
        "media_evidence": first_hit(media_hits),
        "uses_rtl": bool(rtl_hits),
        "rtl_evidence": first_hit(rtl_hits),
        "supports_ie": False,  # macOS WKWebView / modern browsers only (see CLAUDE.md)
    }


# ---------------------------------------------------------------------------
# Wt loader probes (re-verified live so upgrades are reflected)
# ---------------------------------------------------------------------------

def probe_wt_loaders(wt_wt_dir):
    def g1(fname, pattern):
        p = os.path.join(wt_wt_dir, fname) if wt_wt_dir else None
        if not p or not os.path.isfile(p):
            return None
        rx = re.compile(pattern)
        with open(p, "r", errors="ignore") as fh:
            for n, line in enumerate(fh, 1):
                if rx.search(line):
                    return "{}:{}".format(fname, n)
        return None

    return {
        "font_awesome": g1("WIcon.C", r'font-awesome'),
        "ie_css": g1("WCssTheme.C", r'agentIsIElt'),
        "stripes": g1("WCssTheme.C", r'stripes/stripe|alternatingRowColors'),
        "jplayer": g1("WMediaPlayer.C", r'jPlayer/'),
    }


# ---------------------------------------------------------------------------
# CSS reachability (url(...) + @import) from the active-theme entry stylesheets
# ---------------------------------------------------------------------------

URL_RX = re.compile(r"""url\(\s*['"]?([^'")]+)['"]?\s*\)""")
IMPORT_RX = re.compile(r"""@import\s+(?:url\()?\s*['"]([^'"]+)['"]""")


def norm_relpath(p):
    return os.path.normpath(p).replace(os.sep, "/")


def css_reachable(wt_res, entry_css):
    """Return {resource_relpath: "referenced by <entry> via url()/@import"}."""
    reached = {}
    seen_css = set()
    work = [(c, c) for c in entry_css if os.path.isfile(os.path.join(wt_res, c))]
    while work:
        css_rel, origin = work.pop()
        if css_rel in seen_css:
            continue
        seen_css.add(css_rel)
        css_dir = os.path.dirname(css_rel)
        try:
            with open(os.path.join(wt_res, css_rel), "r", errors="ignore") as fh:
                text = fh.read()
        except OSError:
            continue
        for ref in URL_RX.findall(text) + IMPORT_RX.findall(text):
            ref = ref.split("?")[0].split("#")[0].strip()
            if not ref or ref.startswith(("data:", "http:", "https:", "//")):
                continue
            tgt = norm_relpath(os.path.join(css_dir, ref))
            if tgt.startswith("../"):
                continue  # escapes resources/ -- not ours
            full = os.path.join(wt_res, tgt)
            if tgt.endswith(".css") and os.path.isfile(full):
                work.append((tgt, origin))
            if os.path.isfile(full):
                reached.setdefault(tgt, "url()/@import from {}".format(origin))
    return reached


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

RTL_RX = re.compile(r'(^|[-_./])rtl([-_.]|$)', re.IGNORECASE)
IE_CSS = ("wt_ie.css", "wt_ie6.css")

# resources/<dir> bundles that exactly one Wt widget pulls in -- removable when
# that widget is unused (verified live against the current InterSpec source).
WIDGET_BUNDLES = {
    "font-awesome": {"widget": "Wt::WIcon", "loader": "font_awesome",
                     "used": lambda isp: isp["uses_wicon"] or isp["uses_fa_names"]},
    "jPlayer":      {"widget": "Wt::WMediaPlayer", "loader": "jplayer",
                     "used": lambda isp: isp["uses_media_player"]},
}


def _verdict_for(rel, parts, base, isp, wt, active, entry_set, reached):
    # ---- REMOVE rules (app-wide conditions override CSS reachability) ----
    if parts[0] == "themes" and len(parts) > 2 and parts[1] not in active:
        return "REMOVE", ("theme '{}' is never selected (InterSpec sets only {} "
                          "via setCssTheme)".format(parts[1], sorted(active)))
    if base in IE_CSS and not isp["supports_ie"]:
        return "REMOVE", ("IE<9-only stylesheet (Wt {} guards it with "
                          "agentIsIElt(9)/IE6); InterSpec runs WebKit/modern "
                          "only".format(wt["ie_css"]))
    if RTL_RX.search(rel) and not isp["uses_rtl"]:
        return "REMOVE", ("RTL-only asset; InterSpec never sets RightToLeft "
                          "layout, so the rule referencing it never matches")
    if parts[0] in WIDGET_BUNDLES:
        b = WIDGET_BUNDLES[parts[0]]
        if not b["used"](isp):
            return "REMOVE", ("only {} loads {} (Wt {}); current InterSpec source "
                              "uses no {}".format(b["widget"], parts[0],
                                                  wt.get(b["loader"]), b["widget"]))
        return "KEEP", "{} loads {}; used by InterSpec".format(b["widget"], parts[0])

    # ---- KEEP rules ----
    if rel in entry_set:
        return "KEEP", "Wt core/active-theme entry stylesheet"
    if rel in reached:
        return "KEEP", "referenced by active CSS ({})".format(reached[rel])
    if parts[0] == "themes" and parts[1] in active and \
            ("stripes" in parts or "no-stripes" in parts):
        return "KEEP", ("row-striping background set in C++ by "
                        "alternatingRowColors (Wt {}); height-parameterized -> "
                        "confirm exact set via runtime audit".format(wt["stripes"]))

    # ---- default: keep but flag for review ----
    return "REVIEW", ("not reached from active CSS and no rule matched; kept "
                      "pending runtime audit")


def classify(wt_res, files, isp, wt):
    active = set(isp["active_themes"])
    entry_css = (
        ["form.css", "html4_default.css",
         "transitions.css", "moz-transitions.css", "webkit-transitions.css"]
        + ["themes/{}/wt.css".format(t) for t in active]
    )
    reached = css_reachable(wt_res, entry_css)
    entry_set = set(c for c in entry_css if os.path.isfile(os.path.join(wt_res, c)))

    out = {}
    for rel in files:
        verdict, reason = _verdict_for(
            rel, rel.split("/"), os.path.basename(rel), isp, wt,
            active, entry_set, reached)
        out[rel] = {"verdict": verdict, "reason": reason}
    return out


# ---------------------------------------------------------------------------
# Prune-list compaction (whole-dir where every child is REMOVE)
# ---------------------------------------------------------------------------

def compact_remove(files_by_verdict, all_files):
    """Collapse REMOVE files to the shallowest fully-removable directories."""
    remove_set = set(p for p, v in files_by_verdict.items() if v["verdict"] == "REMOVE")

    # every ancestor directory that appears in the tree
    dirs = set()
    for f in all_files:
        d = os.path.dirname(f)
        while d:
            dirs.add(d)
            d = os.path.dirname(d)

    def members(d):
        return [f for f in all_files if f == d or f.startswith(d + "/")]

    full = set(d for d in dirs if members(d) and all(m in remove_set for m in members(d)))
    # keep only the maximal full dirs (parent not also full)
    maximal = sorted(d for d in full if os.path.dirname(d) not in full)

    prune = list(maximal)
    for f in sorted(remove_set):
        if not any(f == d or f.startswith(d + "/") for d in maximal):
            prune.append(f)
    return sorted(prune)


# ---------------------------------------------------------------------------
# Manifest build / IO
# ---------------------------------------------------------------------------

def wt_refs(basenames, wt_wt):
    """First Wt-source site mentioning each basename (single walk)."""
    out = {b: None for b in basenames}
    if not wt_wt:
        return out
    rxs = {b: re.compile(re.escape(b)) for b in basenames}
    for path in _iter_source_files(wt_wt, (".cpp", ".h", ".hpp", ".C", ".cc", ".css")):
        try:
            text = open(path, errors="ignore").read()
        except OSError:
            continue
        rel = os.path.relpath(path, wt_wt)
        for b, rx in rxs.items():
            if out[b] is None:
                m = rx.search(text)
                if m:
                    out[b] = "{}:{}".format(rel, text.count("\n", 0, m.start()) + 1)
    return out


def enrich_review(verdicts, wt_wt):
    """Attach the Wt loader site to REVIEW files so triage has a starting point."""
    review = [r for r, v in verdicts.items() if v["verdict"] == "REVIEW"]
    refs = wt_refs(sorted({os.path.basename(r) for r in review}), wt_wt)
    for r in review:
        site = refs.get(os.path.basename(r))
        if site:
            verdicts[r]["reason"] += "; Wt loader: {}".format(site)


def build_manifest(wt_res, wt_wt, interspec_root):
    files = sorted(
        norm_relpath(os.path.relpath(os.path.join(b, f), wt_res))
        for b, _d, fs in os.walk(wt_res) for f in fs
        if not f.startswith("."))
    isp = probe_interspec(interspec_root)
    wt = probe_wt_loaders(wt_wt)
    verdicts = classify(wt_res, files, isp, wt)
    enrich_review(verdicts, wt_wt)
    prune = compact_remove(verdicts, files)
    summary = {v: sum(1 for x in verdicts.values() if x["verdict"] == v)
               for v in ("KEEP", "REMOVE", "REVIEW")}
    return {
        "active_themes": isp["active_themes"],
        "interspec_conditions": {k: isp[k] for k in (
            "uses_wicon", "uses_fa_names", "uses_media_player",
            "uses_rtl", "supports_ie")},
        "wt_loaders": wt,
        "summary": {"total": len(files), **summary},
        "prune_list": prune,
        "files": verdicts,
    }


def write_outputs(manifest, manifest_path, cmake_path):
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    lines = [
        "# GENERATED by target/tools/wt_resource_audit/wt_resource_audit.py -- do not edit by hand.",
        "# Wt resources/ paths InterSpec does not use; removed from the bundled copy.",
        "# Re-generate after a Wt upgrade:  python3 target/tools/wt_resource_audit/wt_resource_audit.py",
        "set( WT_RESOURCES_PRUNE_LIST",
    ]
    lines += ['  "{}"'.format(p) for p in manifest["prune_list"]]
    lines += [")", ""]
    with open(cmake_path, "w") as fh:
        fh.write("\n".join(lines))


def report(manifest):
    s = manifest["summary"]
    print("Wt resource audit: {total} files | KEEP {KEEP}  REMOVE {REMOVE}  "
          "REVIEW {REVIEW}".format(**s))
    print("active theme(s): {}".format(manifest["active_themes"]))
    print("InterSpec conditions: {}".format(manifest["interspec_conditions"]))
    print("\nprune list ({} entries):".format(len(manifest["prune_list"])))
    for p in manifest["prune_list"]:
        print("  REMOVE  {}".format(p))
    review = sorted(p for p, v in manifest["files"].items()
                    if v["verdict"] == "REVIEW")
    if review:
        print("\nREVIEW (kept; confirm via runtime audit):")
        for p in review:
            print("  REVIEW  {}".format(p))


def check_drift(current, manifest_path):
    if not os.path.isfile(manifest_path):
        print("DRIFT: no committed manifest at {}".format(manifest_path))
        return 1
    with open(manifest_path) as fh:
        committed = json.load(fh)
    cur_f, com_f = current["files"], committed.get("files", {})
    drift = []
    for p in sorted(set(cur_f) - set(com_f)):
        drift.append("NEW/unclassified resource: {} ({})".format(
            p, cur_f[p]["verdict"]))
    for p in sorted(set(com_f) - set(cur_f)):
        drift.append("MISSING resource (was {}): {}".format(com_f[p]["verdict"], p))
    for p in sorted(set(cur_f) & set(com_f)):
        if cur_f[p]["verdict"] != com_f[p]["verdict"]:
            drift.append("VERDICT CHANGED {}: {} -> {}".format(
                p, com_f[p]["verdict"], cur_f[p]["verdict"]))
    if drift:
        print("Wt resource audit DRIFT detected ({}):".format(len(drift)))
        for d in drift:
            print("  " + d)
        print("\nRe-run without --check to regenerate the manifest, then re-triage.")
        return 1
    print("Wt resource audit: no drift vs committed manifest ({} files).".format(
        current["summary"]["total"]))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wt-resources", help="Path to Wt resources/ dir")
    ap.add_argument("--wt-source", help="Path to Wt source checkout (has src/Wt)")
    ap.add_argument("--interspec-root", default=REPO_ROOT)
    ap.add_argument("--manifest",
                    default=os.path.join(SCRIPT_DIR, "wt_resources_manifest.json"))
    ap.add_argument("--cmake-out",
                    default=os.path.join(REPO_ROOT, "cmake", "wt_resources_prune_list.cmake"))
    ap.add_argument("--check", action="store_true",
                    help="Compare to committed manifest; exit non-zero on drift")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args()

    wt_res, wt_wt = resolve_inputs(args)
    if not wt_res or not os.path.isdir(wt_res):
        sys.exit("ERROR: could not find Wt resources dir (use --wt-resources)")
    if not wt_wt:
        print("WARNING: no Wt src/Wt dir found; loader file:line evidence omitted",
              file=sys.stderr)
    print("Wt resources : {}".format(wt_res), file=sys.stderr)
    print("Wt source    : {}".format(wt_wt), file=sys.stderr)
    print("InterSpec    : {}".format(args.interspec_root), file=sys.stderr)

    manifest = build_manifest(wt_res, wt_wt, args.interspec_root)

    if args.check:
        sys.exit(check_drift(manifest, args.manifest))

    report(manifest)
    if not args.no_write:
        write_outputs(manifest, args.manifest, args.cmake_out)
        print("\nwrote {}\nwrote {}".format(args.manifest, args.cmake_out))


if __name__ == "__main__":
    main()
