#ifndef RefLib_h
#define RefLib_h
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

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <optional>

#include "nlohmann/json.hpp"

#include "InterSpec/PeakDef.h"

namespace SpecUtils{ class Measurement; }

/** The precomputed reference-line library (see refgen/make_ref_lines.cpp for the JSON format),
 and the port of InterSpec's `PeakSearchGuiUtils::assign_nuc_from_ref_lines`.
 */
namespace RefLib
{
  struct Line
  {
    /** Energy the line is drawn at (for escape lines, the escape-peak energy). */
    double energy = 0.0;
    /** Intensity relative to the source's strongest photon. */
    double intensity = 0.0;
    bool is_xray = false;
    PeakDef::Source source;
  };//struct Line

  struct Source
  {
    std::string parent;
    std::vector<Line> lines;
  };//struct Source

  class Library
  {
  public:
    /** Parses the array written by refgen; throws on malformed input. */
    void load( const nlohmann::json &lib );

    /** Lookup ignoring case, spaces, and dashes; nullptr if not found. */
    const Source *find( const std::string &parent ) const;

    size_t size() const { return m_sources.size(); }

    /** All sources, keyed by parent name, lower-cased and without spaces or dashes. */
    const std::map<std::string,Source> &sources() const { return m_sources; }

  private:
    std::map<std::string,Source> m_sources;
  };//class Library


  /** A reference-line source currently shown in the GUI, with the color it is drawn in. */
  struct Shown
  {
    const Source *source = nullptr;
    std::string color;
  };//struct Shown

  /** Assigns the best-matching source from `shown` to `peak` (sets its source, and color).

   Port of `PeakSearchGuiUtils::assign_nuc_from_ref_lines`: lines are scored by
   (0.25*sigma + |mean - E|)/intensity (x-rays de-weighted), S.E./D.E. of high-energy lines are
   considered, and the expected energy must be within the peak's ROI.  If an existing peak has
   already claimed the best line, but the relative amplitudes suggest the two peaks should be
   swapped, returns that existing peak along with the source it should get instead.

   @param isHPGe Whether escape peaks should be considered for lines above 1255 keV.
   */
  std::unique_ptr<std::pair<std::shared_ptr<const PeakDef>,PeakDef::Source>>
    assign_source( PeakDef &peak,
                   const std::vector<std::shared_ptr<const PeakDef>> &previous_peaks,
                   const std::vector<Shown> &shown,
                   const bool isHPGe );

  /** Sources the shown reference lines suggest for `peak`, best first: of each nuclide, element, or
   reaction with a shown line within 4 sigma of the peak mean, the line that best matches the peak,
   scored as in `assign_source` (but without looking for S.E. or D.E. peaks of full-energy lines).
   Leaves out the nuclide, element, or reaction the peak already has.  Like InterSpec's
   `IsotopeId::peakCandidateSourceFromRefLines`, but with one line per nuclide.
   */
  std::vector<PeakDef::Source> suggest_sources( const PeakDef &peak, const std::vector<Shown> &shown );

  /** The source for `peak` that a user typed, read as InterSpec's
   `PeakModel::setNuclideXrayReaction` does, but using the library's lines.

   The text is a library source, nuclide, element (for its x-rays), or reaction, optionally with an
   energy, and "x-ray", "S.E.", "D.E.", or "annih."; e.g., "Cs137", "cs-137 661.66 keV",
   "Ba133 xray", "Th232 2614 S.E." (energies are of the photon, not of the escape peak), "Pb 75",
   or "H(n,g)".  A `source_label` reads back as that source.  The line is the one at the energy given
   (to its 0.01 keV), else of those within 4 sigma of the energy (or of the peak mean, plus 511 or
   1022 keV for S.E. or D.E.), the one with the smallest (0.1*window + distance)/intensity, else the
   nearest.

   @param from Set to the library source the line is from.
   @param note Set to a message for the user if the line is not near the energy, else cleared.
   @returns The source, or nullopt for empty text, or "none".  Throws std::runtime_error if the
            library has no such source.
   */
  std::optional<PeakDef::Source> source_from_text( const Library &lib, const PeakDef &peak,
                                                   const std::string &text, const Source *&from,
                                                   std::string &note );

  /** The label peaks show for a source, e.g., "Cs137 661.66 keV", "S.E. Th232 2614.53 keV",
   "Cs137 32.19 keV x-ray", or "Na22 511.00 keV annih."; empty for no source. */
  std::string source_label( const PeakDef::Source &src );

  /** Whether both are the same nuclide, element, or reaction (not necessarily the same line). */
  bool same_nuclide( const PeakDef::Source &a, const PeakDef::Source &b );
}//namespace RefLib

#endif //RefLib_h
