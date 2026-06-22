#ifndef PeakFitImproveData_h
#define PeakFitImproveData_h
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


#include <tuple>
#include <string>
#include <vector>
#include <istream>

#include "PeakFitImprove.h"

namespace SpecUtils
{
  class SpecFile;
  class Measurement;
}


namespace PeakFitImproveData
{
namespace G2k
{
struct G2kPeak
{
  char multiplet;
  int PeakNum;
  float Energy, PeakSig, FWHM, ContinuumCounts, NetPeakArea, NetAreaError, NetCountRate;
};

std::vector<G2kPeak> g2k_peaks_from_file( std::istream &strm );
}//namespace G2k


ExpectedPhotopeakInfo create_expected_photopeak( const InjectSourceInfo &info, const std::vector<PeakTruthInfo> &lines );


/** Filters expected photopeaks down to the set used for scoring, dropping peaks below a det-type
 low-energy threshold (30 keV for HPGe, 50 keV otherwise) where peaks are unreliable (poor
 resolution, fast-changing efficiency, x-ray overlap).  Exception: if the source's dominant
 (largest-area) expected peak is itself below the threshold, all peaks are kept (so a low-energy-
 dominant source is not discarded).  Mirrors the filter in FitPeaksForNuclideDev so the GA objective
 and the dev eval path score against the same peaks. */
std::vector<ExpectedPhotopeakInfo> filter_photopeaks_for_scoring(
    const std::vector<ExpectedPhotopeakInfo> &expected,
    PeakFitUtils::CoarseResolutionType det_type );


/** Fraction (0..1) of "definitely wanted" expected peak-AREA that was not detected, used for the GA
 miss-penalty term.  Denominator: summed area of def-wanted peaks (det-type gated:
 nsigma>def_want_nsigma and area>min_def_wanted_counts) in `scoring_peaks`.  Numerator: summed area of
 `missed_def_wanted` (i.e. CandidatePeakScore::def_expected_but_not_detected, after escape-peak
 correction).  Returns 0 when there is no def-wanted area. */
double missed_def_wanted_area_fraction(
    const std::vector<ExpectedPhotopeakInfo> &scoring_peaks,
    const std::vector<ExpectedPhotopeakInfo> &missed_def_wanted,
    PeakFitUtils::CoarseResolutionType det_type );

/**

 Expects:
 `base_name + ".pcf"`
 `base_name + "_gamma_lines.csv"`
 to exist.
 */
InjectSourceInfo parse_inject_source_files( const std::string &base_name,
                                           const std::vector<std::pair<float,float>> &deviation_pairs );


/**
 
 */
std::tuple<std::vector<DetectorInjectSet>,std::vector<DataSrcInfo>>
 load_inject_data_with_truth_info( const std::string &base_dir,
                                  const std::vector<std::string> &wanted_detectors,
                                  const std::vector<std::string> &live_times,
                                  const std::vector<std::string> &wanted_cities );


/** Write compacted data to output_dir, mirroring the original directory structure.

 The output PCF files contain only the 4 needed measurements (with deviation-corrected
 energy calibration baked in), and truth CSV files include pre-computed ExpectedPhotopeakInfo.
 */
void write_compact_data( const std::vector<DetectorInjectSet> &inject_sets,
                         const std::vector<DataSrcInfo> &input_srcs,
                         const std::string &output_dir );


/** Returns true if the base_dir appears to contain compacted data format. */
bool is_compact_data_directory( const std::string &base_dir );


/** Load compacted data from a compact directory.
 Produces the same result as load_inject_data_with_truth_info, but from the compacted format.
 */
std::tuple<std::vector<DetectorInjectSet>, std::vector<DataSrcInfo>>
  load_compact_data( const std::string &base_dir,
                     const std::vector<std::string> &wanted_detectors,
                     const std::vector<std::string> &live_times,
                     const std::vector<std::string> &wanted_cities );


/** Unified loader that auto-detects format (compact vs original) and delegates. */
std::tuple<std::vector<DetectorInjectSet>, std::vector<DataSrcInfo>>
  load_data( const std::string &base_dir,
             const std::vector<std::string> &wanted_detectors,
             const std::vector<std::string> &live_times,
             const std::vector<std::string> &wanted_cities );


#if( WRITE_ALL_SPEC_TO_HTML ) //Defined in PeakFitImprove.h
void write_html_summary( const std::vector<DataSrcInfo> &src_infos );
#endif
}//namespace PeakFitImproveData

#endif //PeakFitImproveData_h
