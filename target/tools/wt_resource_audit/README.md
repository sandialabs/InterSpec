# Wt resource pruning & upgrade checks

InterSpec bundles Wt's `resources/` directory (~5.9 MB) into every packaged app
but only loads a small slice of it. `wt_resource_audit.py` decides — with
evidence — which files are unused, and the CMake build removes them from the
bundled copy at configure time (~5.4 MB / ~91 % smaller).

## Pieces

| File | Role |
|------|------|
| `target/tools/wt_resource_audit/wt_resource_audit.py` | Static analyzer: traces how Wt loads each resource and whether the current InterSpec source triggers that path. Emits the manifest + prune list; `--check` is the drift gate. |
| `target/tools/wt_resource_audit/wt_resources_manifest.json` | Generated, committed. Per-file `KEEP`/`REMOVE`/`REVIEW` verdict + evidence (Wt loader `file:line`, InterSpec condition). The baseline for drift detection. |
| `cmake/wt_resources_prune_list.cmake` | Generated, committed. `WT_RESOURCES_PRUNE_LIST` consumed by the prune. |
| `cmake/PruneWtResources.cmake` | `prune_wt_resources()` — copies Wt resources into the build tree minus the prune list. Hooked from the root `CMakeLists.txt`. |
| `target/testing/.../WtResourcePruneManifest` | CTest case running `--check`; fails on drift. |

The prune is controlled by the CMake option `PRUNE_WT_RESOURCES` (default `ON`).
Build with `-DPRUNE_WT_RESOURCES=OFF` to bundle the full Wt resources (e.g. to
experiment with the `bootstrap`/`polished` themes).

## What gets removed, and why it's safe

A file is only `REMOVE`d when an app-wide condition makes its Wt loader
unreachable — never a guess. Current removals (see the manifest for per-file
evidence): the `bootstrap` & `polished` themes (InterSpec sets only
`setCssTheme("default")`), `font-awesome` (only `Wt::WIcon` loads it; InterSpec
uses none), `jPlayer` (only `Wt::WMediaPlayer` loads it; the video feature was
removed), the IE-only stylesheets, and the RTL image variants (InterSpec is
never RightToLeft). `REVIEW` files are kept.

## Upgrading Wt (the part that can go stale)

The manifest is derived from the Wt source + InterSpec source, so it must be
re-validated whenever Wt changes. The CTest / `--check` gate fails loudly if you
forget. Procedure:

1. **Bump the Wt pin in both places:**
   - FetchContent: `GIT_TAG` (+ the `# ---- Wt <ver> ----` header) in
     `cmake/FetchInterSpecDeps.cmake`.
   - Prebuilt install scripts: `target/dep_build/dep_build_{macOS.sh,linux.sh,msvc2022.bat,msvc2022_x86.bat}`.
   Reconfigure so the new Wt is fetched/installed.

2. **See what changed in the resource set** (in a Wt checkout):
   ```
   git diff --stat <old-tag>..<new-tag> -- resources/
   git diff <old-tag>..<new-tag> -- resources/themes/default '*.css'   # url() changes
   ```

3. **Re-run the analyzer** (regenerates the manifest + prune list):
   ```
   python3 target/tools/wt_resource_audit/wt_resource_audit.py
   ```
   Review the diff of `target/tools/wt_resource_audit/wt_resources_manifest.json`. Triage any new files or
   `REVIEW` rows the resource diff surfaced.

4. **Confirm no drift, then re-run the runtime audit** (below) on a fresh build.

5. Commit the regenerated `target/tools/wt_resource_audit/wt_resources_manifest.json` and
   `cmake/wt_resources_prune_list.cmake`.

### Just checking (no regenerate)

```
python3 target/tools/wt_resource_audit/wt_resource_audit.py --check     # exit!=0 on drift
```
or run the test: `ctest -R WtResourcePruneManifest`. Inputs auto-detect; override
with `--wt-resources <dir>` / `--wt-source <wt-checkout>` if needed.

## Runtime audit (confirm nothing removed is actually requested)

Static analysis is the gate; confirm it empirically:

1. Build/run the local web server with `-DPRUNE_WT_RESOURCES=OFF` (full resources)
   and watch resource requests (server access log or browser devtools network).
2. Exercise the app comprehensively in **desktop, phone (`?isphone=1`) and tablet
   (`?istablet=1`)** modes — open every tool/dialog, tree/table views, etc.
3. No `resources/...` request may hit a `REMOVE` path. `KEEP`/`REVIEW` files never
   requested are candidates for further pruning (move to `REMOVE` only with
   evidence).
4. Rebuild with `-DPRUNE_WT_RESOURCES=ON` and smoke-test all three modes: zero
   404s for any pruned path.
