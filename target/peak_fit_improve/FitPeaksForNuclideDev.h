#ifndef FitPeaksForNuclideDev_h
#define FitPeaksForNuclideDev_h
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

#include <memory>
#include <string>
#include <functional>

#define OPENGA_EXTERN_LOCAL_VARS 1

#include "openGA.hpp"

#include <vector>

#include "InterSpec/PeakDef.h"

#include "PeakFitImprove.h"

namespace SpecUtils
{
  class Measurement;
}



namespace FitPeaksForNuclideDev
{
  void eval_peaks_for_nuclide( const std::vector<DataSrcInfo> &srcs_info );
}//namespace InitialFit_GA

#endif //FitPeaksForNuclideDev_h
