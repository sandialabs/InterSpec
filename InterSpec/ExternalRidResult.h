#ifndef ExternalRidResult_h
#define ExternalRidResult_h
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
#include <variant>

#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/ReactionGamma.h"

/** Struct to represent a nuclide/source ID from an automated RID algorithm; practically either a detectors on-board RID, or GADRAS Full-Spectrum.
 
 TODO: make a number of these fields to use std::optional<>
 */
struct ExternalRidIsotope
{
  typedef std::variant<std::monostate,const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> SrcVariant;
  
  std::string name;
  std::string type;
  std::string confidenceStr;
  
  double countRate;
  double confidence;
  
  /** The source resolved to a nuclide, x-ray, or reaction, if possible. */
  SrcVariant source;
  
  const SandiaDecay::Nuclide *nuclide() const;
  
  const SandiaDecay::Element *element() const;
  
  const ReactionGamma::Reaction *reaction() const;
  
  bool is_null() const;
  
  std::string source_name() const;
  
  /** fills out `source`, and potentually `type` from `name` */
  void init();
  
  bool operator==(const ExternalRidIsotope& other) const;
};//struct ExternalRidIsotope

/** Struct to represent the results from an automated RID algorithm; practically either a detectors on-board RID, or GADRAS Full-Spectrum.
 
 TODO: make a number of these fields to use std::optional<>
 */
struct ExternalRidResults
{
  std::string algorithmName;
  
  int code = -1;
  int analysisError = 0;
  std::string errorMessage;
  std::vector<std::string> analysisWarnings;
  std::string drf;
  double chi2 = -1.0;
  //estimatedDose: in PhysicalUnits (i.e., 1.0E-6*PhysicalUnits::rem/PhysicalUnits::hour)
  double estimatedDose = -1.0;
  double stuffOfInterest = 1.0;
  double alarmBasisDuration = -1.0;
  
  std::string foregroundDescription;
  std::string foregroundTitle;
  std::string isotopeString;
  
  std::string analysisType;
  
  std::vector<ExternalRidIsotope> isotopes;
};//struct ExternalRidResults

#endif //ExternalRidResult_h
