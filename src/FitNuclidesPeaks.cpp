/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <set>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

#include "InterSpec/FitNuclidesPeaks.h"


using namespace std;

class PeakDef;
class DetectorPeakResponse;
namespace SandiaDecay
{
  class Nuclide; 
}//namespace SandiaDecay


namespace FitNuclidesPeaks
{
 
/*
 
class FitNuclidesPeaksFcn
  : public ROOT::Minuit2::FCNBase
{
  //Detector calib offset add
  //Detector calib offset gain multiple
  //Detector resolution variable 1
  //Detector resolution variable 2
  //Detector resolution variable 3
  //Shielding Atomic Number
  //Shielding Areal Density
  //Detector Efficiency Response variable 1
  //Detector Efficiency Response variable 2
  //Detector Efficiency Response variable 3
  //Detector Efficiency Response variable 4
  //Nuclide 1 Strength
  //Nuclide 1 Age
  //Nuclide 2 Strength
  //Nuclide 2 Age
  //...
   
  //Estimate initial nuclide strength
   
public:
  FitNuclidesPeaksFcn( const std::shared_ptr<const Measurement> &data,
                       std::shared_ptr<const DetectorPeakResponse> detector,
                       const float lowerEnergy,
                       const float upperEnergy );
    
  FitNuclidesPeaksFcn &operator=( const FitNuclidesPeaksFcn &rhs );
    
  virtual double Up() const;
  virtual double operator()( const std::vector<double> &params ) const;
  
  
  //Creates the peaks from the given parameters, and if passed in the parameter
  //  errors, adds those to the peaks as well.
  std::vector< std::shared_ptr<PeakDef> > parametersToPeaks(
                                             const double *pars,
                                             const double *errors = 0 ) const;
  
  //If in initialAge <= 0.0, then the recomended age will be be used
  void addNuclide( const SandiaDecay::Nuclide *nuclide );
  
  ROOT::Minuit2::MnUserParameters startingFitParameters();
  
  std::shared_ptr<PeakDef> peaksFromPars( const double *pars, 
                                            const double *errors = 0 ) const;
  
protected:
  bool m_hasRetrievedStartPars;
  float m_lowerEnergy, m_upperEnergy;
  std::shared_ptr<const Measurement> m_data;
  std::shared_ptr<const DetectorPeakResponse> m_detector;
  
private:
  virtual double DoEval( const double *x ) const;
};//class FitNuclidesPeaksFcn

  
  */
  
  
  
  
  
}//namespace FitNuclidesPeaks

