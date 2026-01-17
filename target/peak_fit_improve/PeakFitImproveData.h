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


#if( WRITE_ALL_SPEC_TO_HTML ) //Defined in PeakFitImprove.h
void write_html_summary( const std::vector<DataSrcInfo> &src_infos );
#endif
}//namespace PeakFitImproveData

#endif //PeakFitImproveData_h
