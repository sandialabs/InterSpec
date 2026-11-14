#ifndef DecayActivityCalc_h
#define DecayActivityCalc_h
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
#include <vector>
#include <optional>

#include <Wt/WString.h>

#include "SandiaDecay/SandiaDecay.h"


/** Non-GUI core of the "Nuclide Decay Info" tool (`DecayWindow` / `DecayActivityDiv`).

 Holds the tool's app-url ("interspec://decay/...") encoding and decoding, and the decay
 calculations the GUI needs, so they can be unit tested without a session.
 */
namespace DecayActivityCalc
{
  /** What the "Activity Chart" y-axis shows; the values are the indices of the GUI's combo box. */
  enum class ChartType : int { Activity = 0, Gamma = 1, Beta = 2, Alpha = 3 };

  /** The tool's tabs; the values are the tab indices. */
  enum class Tab : int { Chart = 0, Chain = 1, Calc = 2 };

  /** The `dispunits` URL value (and `ToolState::displayUnit`) for showing activities relative to
   the total initial activity.
   */
  extern const char * const sm_arbitraryDisplayUnit;

  /** A source nuclide the user has entered. */
  struct SourceNuclide
  {
    const SandiaDecay::Nuclide *nuclide = nullptr;

    /** Activity at time zero, in SandiaDecay units (becquerel). */
    double activity = 0.0;

    /** The activity as the user entered it (kept for display, and to track significant figures).
     Callers should keep it parsing to `activity`; see #normalizeActivityString.
     */
    std::string activityStr;

    /** Whether to display activities derived from this nuclide in curies (else becquerels). */
    bool useCurie = true;

    /** Age of the nuclide at time zero (i.e., how long its progeny have been building in). */
    double age = 0.0;
  };//struct SourceNuclide


  /** The state of the tool that is round-tripped through its app-url. */
  struct ToolState
  {
    Tab tab = Tab::Chart;
    bool gridLines = false;
    bool logY = false;
    ChartType chartType = ChartType::Activity;

    /** An exact name from `PhysicalUnits::sm_activityUnitNameValues` (e.g., "MBq" - note "mBq" is
     different), #sm_arbitraryDisplayUnit, or empty to use the default (MBq or mCi).
     */
    std::string displayUnit;

    /** Time span to show, as text; may be in half-lives of the first nuclide (e.g., "5 hl").
     Empty to use the default span.
     */
    std::string timeSpan;

    std::vector<SourceNuclide> nuclides;
  };//struct ToolState


  /** Parses the tool's app-url.

   Example URIs:
   - "interspec://decay/chain?nuclide=U238&activity=3uCi&initialage=20y&timespan=22y&actunits=ci"
   - "interspec://decay/chart?nuc=Ba133&act=3uCi&nuc=Cs137&act=2uci&actunits=ci&timespan=20y"
   - "interspec://decay/chain?nuc=232-th"

   @param path The path of the URI; one of "chart", "chain", or "calc" (case-insensitive), which
          selects the tab shown.  E.g., for "interspec://decay/chain?nuc=...", the path is "chain".
          Empty selects the chart.
   @param query The query portion of the URI, already url-decoded, in the form
          "key1=value1&key2=value2&..." (a leading '?' is allowed).  Order only matters when more
          than one nuclide is given: a nuclide's activity and age must follow it.
          Keys and values are case-insensitive, except `dispunits` values.
   @param defaultUseBq Whether to use becquerels for a nuclide when neither its `act` value, nor
          the `actunits` key, say which units to use.
   @param errors Messages for anything that could not be understood are appended here; the rest of
          the URL is still used.

   Possible keys:
   - "time", "timespan": the time span to show on the chart or calculator, using time units, or
     half-lives of the first nuclide.  Ex. "timespan=20y", "timespan=5hl", "timespan=200seconds".
   - "dispunits": the exact name of the activity unit to display the chart in (ex. "MBq", "uCi"),
     or "arbitrary".
   - "actunits": "becquerel", "bq", "curie", or "ci".  The units to use for nuclides whose
     activity does not specify them, and, if "dispunits" is not given, picks MBq or mCi for the
     chart.  In URLs without "dispunits" (as written before it existed, when "act" was always
     given in curies), it decides the units of every nuclide.
   - "charttype": "activity" (or "act"), "gamma", "beta", or "alpha"; what the chart y-axis shows.
   - "grid" or "gridlines": "1", "yes", or "true" shows grid lines on the chart.
   - "logy": "1", "yes", or "true" uses a log y-axis on the chart.
   - "iso", "nuc", "isotope", "nuclide": synonyms giving a nuclide; ex. "nuc=Am241", "iso=Co60m",
     "iso=Co-60", "iso=60Co".
   - "initialage", "age": the age of the immediately preceding nuclide at time zero.
   - "act", "activity": the activity of the immediately preceding nuclide, with units; ex.
     "act=5uci", "activity=10kbq".  Defaults to 1 mCi.

   Not currently supported: which nuclide the decay chain or calculator shows when there are
   several, the calculator dates, series the user has hidden on the chart, or the nuclide selected
   in the decay chain.
   */
  ToolState fromAppUrl( const std::string &path, std::string query, const bool defaultUseBq,
                        std::vector<Wt::WString> &errors );

  /** Encodes the state as an app-url, as accepted by #fromAppUrl.

   The result does not include the "interspec://decay/" prefix (so looks like
   "chart?grid=0&logy=0&..."), and is not url-encoded.
   */
  std::string toAppUrl( const ToolState &state );

  /** Returns whether the units written in an activity string are curies (true) or becquerels
   (false), or nothing if the string has neither.  Ex., "3 uCi" -> true, "10kBq" -> false.
   */
  std::optional<bool> activityStringIsCurie( const std::string &activity );

  /** Formats a time for a URL: English units (not the users locale), to six decimals, with
   trailing zeros dropped (ex., "40.42 h").
   */
  std::string printTimeForUrl( const double time );

  /** Makes `src.activityStr` parse to `src.activity`, by formatting it from the activity (in the
   source's units) if it is empty or does not match.
   */
  void normalizeActivityString( SourceNuclide &src );

  /** Returns the time at which the total activity of the mixture falls to `frac` of its initial
   value; used to pick the default time span to show.

   Steps forward in half the longest half-life (or a tenth of the range, if a search range is
   given), then searches the step that crosses `frac`.

   Returns 0.0 for a mixture with no activity, and -1.0 if no time was found.
   Throws if `frac` is not in (0,1], or `searchStart > searchEnd`.
   */
  double findTimeForActivityFrac( const SandiaDecay::NuclideMixture &mix, const double frac,
                                  const double searchStart = -1.0, const double searchEnd = -1.0 );

  /** Returns the number of particles of `type`, per second, a nuclide gives off at the activity
   (SandiaDecay units).

   Matches `SandiaDecay::NuclideMixture::gammas()`/`betas()`/`alphas()`, e.g., gamma records that
   are pure E0 transitions (and so emit no photon) are not counted.
   */
  double particleRate( const SandiaDecay::Nuclide &nuc, const double activity,
                       const SandiaDecay::ProductType type );
}//namespace DecayActivityCalc

#endif //DecayActivityCalc_h
