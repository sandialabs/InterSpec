#ifndef RelActAutoReport_h
#define RelActAutoReport_h
/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <string>
#include <ostream>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"
#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

// Forward declarations
namespace RelActCalcAuto
{
  struct RelActAutoSolution;
}

/** RelActAutoReport: Inja-driven report generation for `RelActCalcAuto::RelActAutoSolution`.

 Mirrors the design of `BatchInfoLog` for shielding/source and peak-fit reports.  Templates live
 at `InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/`.  The default templates
 (`std_rel_eff_summary.tmplt.html` and `std_rel_eff_summary.tmplt.txt`) reproduce the output of
 `RelActAutoSolution::print_html_report()` and `print_summary()` respectively.  Users may supply
 their own template via the `tmplt` (filename or absolute path) and `include_dir` (search root
 for `{% include %}`) parameters of `render_template`.
 */
namespace RelActAutoReport
{
  // ---- Template + Inja environment helpers (mirror BatchInfoLog) ----

  /** Returns the equivalent of `InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/`,
   with trailing path separator (as required by inja).  Returns "" if the static data directory
   has not been initialised.
   */
  std::string default_template_dir();

  /** Resolve a user-supplied template-include-dir string to an absolute path (with trailing
   separator).  Mirrors `BatchInfoLog::template_include_dir`.
   - "" or "none"   -> ""
   - "default"      -> `default_template_dir()`
   - any other path -> normalised with trailing path separator.
   */
  std::string template_include_dir( const std::string &user_path );

  /** Returns an inja environment configured with:
   - the (optional) user template include directory as its search root
   - `set_trim_blocks(true)`
   - `printFixed`, `printCompact` callbacks (delegated to `BatchInfoLog`)
   - `pct`, `safe_html`, `scientific`, `format`, `default` callbacks
   - the default templates pre-registered as named includes:
       "default-rel-act-auto-html-results"  -> std_rel_eff_summary.tmplt.html
       "default-rel-act-auto-txt-results"   -> std_rel_eff_summary.tmplt.txt
   */
  inja::Environment get_default_inja_env( const std::string &tmplt_dir );

  /** Loads and returns the JS+CSS assets needed for the embedded D3 spectrum chart.
   Returns key/value pairs, where the key is the JSON field name to use.  Same shape as
   `BatchInfoLog::load_spectrum_chart_js_and_css()`.
   */
  std::vector<std::pair<std::string,std::string>> load_spectrum_chart_js_and_css();

  /** Loads and returns the JS+CSS assets needed for the embedded RelEffPlot chart.
   Returns key/value pairs:
       "RelEffPlot_JS"   -> contents of InterSpec_resources/RelEffPlot.js
       "RelEffPlot_CSS"  -> contents of InterSpec_resources/RelEffPlot.css
   */
  std::vector<std::pair<std::string,std::string>> load_rel_eff_plot_js_and_css();

  /** Convert a `RelActAutoSolution` to the JSON payload consumed by the templates.
   Populates a comprehensive set of fields - both raw numerical values and pre-formatted "_str"
   convenience fields - so that custom templates can either format values themselves (via the
   `printFixed`/`printCompact` callbacks) or use the pre-formatted strings directly.

   The dump always contains the full set of keys, including the JS/CSS asset blobs (`assets.D3_JS`,
   `assets.SpectrumChart_JS`, `assets.SpectrumChart_CSS`, `assets.RelEffPlot_JS`,
   `assets.RelEffPlot_CSS`) and `spectrum_chart.set_js`. These are large strings (~hundreds of KB
   each); strip them with `jq 'del(.assets, .spectrum_chart.set_js)'` if you want a readable dump.
   */
  nlohmann::json solution_to_json( const RelActCalcAuto::RelActAutoSolution &solution );

  /** Render a template against a pre-built JSON payload.

   This is the lowest-level rendering entry point.  Typical usage:

   ```cpp
   const auto data = RelActAutoReport::solution_to_json( solution );
   auto env = RelActAutoReport::get_default_inja_env( include_dir );
   RelActAutoReport::render_template( out, env, data, "my.tmplt.html", include_dir );
   ```

   Pulling these three steps apart lets callers (a) build the JSON once and render multiple
   templates against it, and (b) extend `data` with extra keys before rendering.

   @param env an inja environment.  Should have been built by `get_default_inja_env(include_dir)`
     so that the bundled named templates and standard callbacks are registered, and so that any
     `{% include %}` directives in `tmplt` resolve relative to `include_dir`.
   @param data the JSON payload to render against; usually `solution_to_json(solution)`.
   @param tmplt selects the template to render. Accepted values:
     - `""`             -> built-in HTML report (same as `"html"`)
     - `"html"`         -> shorthand for the bundled
                           `InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/std_rel_eff_summary.tmplt.html`
     - `"txt"` / `"text"` -> shorthand for the bundled `std_rel_eff_summary.tmplt.txt`
     - `"json"`         -> `data.dump(2)` (skips Inja entirely)
     - any other string is treated as a path:
       1. if `include_dir` is non-empty, looked for under `<include_dir>/<tmplt>` first
       2. otherwise, accepted as an absolute path if it exists
       3. otherwise, looked for in the bundled
          `InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/` directory
   @param include_dir search root for resolving filenames in `tmplt`.  Should match the
     `include_dir` that was passed to `get_default_inja_env` to build `env` so that user-template
     `{% include %}` directives also resolve.  Special values `"default"` and `"none"` are
     accepted (see `template_include_dir`).
   */
  std::string render_template( inja::Environment &env,
                               const nlohmann::json &data,
                               const std::string &tmplt = "html",
                               const std::string &include_dir = "" );

  /** Render the template and write it to `out`.  See the string-returning overload for parameter docs. */
  void render_template( std::ostream &out,
                        inja::Environment &env,
                        const nlohmann::json &data,
                        const std::string &tmplt = "html",
                        const std::string &include_dir = "" );


  // ---- Inja callback helpers (exposed for testing / reuse) ----

  /** Inja callback: format `value*100` to `precision` significant figures (default 4). */
  std::string pct_callback( std::vector<const nlohmann::json *> &args );

  /** Inja callback: HTML-sanitize a string (escapes &, <, >, ', "). */
  std::string safe_html_callback( std::vector<const nlohmann::json *> &args );

  /** Format `value` in scientific notation with `precision` significant figures.  Used by the
   legacy `scientific` callback. */
  std::string format_scientific( double value, int precision = 6 );

  /** Format `value*100` as a percentage string with `precision` decimal places. */
  std::string format_percentage( double value, int precision = 2 );

  /** Format a duration in microseconds using `PhysicalUnits::printToBestTimeUnits`. */
  std::string format_duration( int64_t microseconds );

  /** Apply InterSpec's standard HTML-sanitisation to `val` in-place (escapes &, <, >, ', "). */
  void html_sanitize( std::string &val );

}//namespace RelActAutoReport

#endif // RelActAutoReport_h
