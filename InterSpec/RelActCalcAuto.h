#ifndef RelActCalcAuto_h
#define RelActCalcAuto_h
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

#include <vector>

#include "InterSpec/PeakDef.h" //for PeakContinuum::OffsetType


// Forward declerations
namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils

namespace SandiaDecay
{
  struct Nuclide;
}

/*
 We need to take in:
 
 foreground_spectrum;
 background_spectrum; //optional
 
 struct RoiRanges{ double lower_energy, upper_energy; continuum_type; bool force_full_range; }
 vector<RoiRanges>: roi_ranges;
 
 struct NucInputInfo{ const SandiaDecay::Nuclide *nuc; double age; bool fit_age; vector<double> gammas_to_release; }
 vector<NucInputInfo> nuclides;
 
 struct FreeFloatPeak{ double energy; bool release_fwhm; }
 vector<FreeFloatPeak> free_floating_peaks;
 
 options{ bool fit_energy_cal; RelEffEqnType; size_t RelEffEqnOrder; FwhmFcnForm/Order; }
 
 
 For deciding peak ROIs:
 - First use whole specified ROI, and fit for a solution
 - Eliminate statistically insignificant gamma lines; then use remaining gammas to decide peak ROI widths (and if should be further broken up?)
 
 To get rid of the degeneracy of rel-eff and rel-act:
 - Add an extra term to the Chi2 that forces the rel-act equation to be 1.0 at the lowest energy; maybe multiple difference of current rel-eff equation from 1.0 by the area of the largest peak in any of the ROIs
 
 Parameters to fit:
 - Energy offset; Energy gain adjust //will be fixed to begin with
 - Fwhm Parameters
 - RelEff coefs
 - Activities
 - Ages //will mostly be fixed
 - Free float peaks {amplitude, width} // if width is negative, then use Fwhm eqn
 
 */

namespace RelActCalcAuto
{

/** Struct to specify an energy range to consider for doing relative-efficiency/activity calc.
 */
struct RoiRange
{
  double lower_energy;
  double upper_energy;
  PeakContinuum::OffsetType continuum_type;
  bool force_full_range;
};//struct RoiRange


/** Struct to specify a nuclide to use for doing relative-efficiency/activity calc.
 */
struct NucInputInfo
{
  const SandiaDecay::Nuclide *nuclide;
  double age;
  bool fit_age;
  std::vector<double> gammas_to_release;
};//struct NucInputInfo

}//namespace RelActCalcAuto

#endif //RelActCalcAuto_h
