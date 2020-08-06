#ifndef IsotopeId_h
#define IsotopeId_h
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

#include <deque>
#include <string>
#include <vector>

#include "InterSpec/ReferenceLineInfo.h"

//Forward Declarations
class PeakDef;
class PeakModel;
class DecayDataBaseServer;
class DetectorPeakResponse;
namespace SpecUtils
{
  class Measurement;
}

namespace SandiaDecay
{
  struct Nuclide;
  struct EnergyRatePair;
}


namespace IsotopeId
{
  void setDataDirectory( const std::string &dir );
  
  
struct NuclideStatWeightPair
{
  const SandiaDecay::Nuclide *nuclide;
  double weight;
  
  NuclideStatWeightPair()
    : nuclide( 0 ), weight( 0.0 ) {}
  NuclideStatWeightPair( const SandiaDecay::Nuclide * const n, const double w )
    : nuclide( n ), weight( w ){}
};//struct NuclideStatWeightPair
  
  
struct PeakToNuclideMatch
{
  double energy;
  std::vector<NuclideStatWeightPair> nuclideWeightPairs;
  
  PeakToNuclideMatch()
    : energy(0.0), nuclideWeightPairs(0) {}
};//struct PeakToNuclideMatch


//suggestNuclides(...): suggests nuclides for 'peak', based on a weighting
//  scheme that involves distance in energy, what other peaks are detected that
//  would be expected (or what expected peaks arent there) under a few shielding
//  scenarious, if it is a characteristic gamma of the nuclide, half-life, and
//  parent (or grandparent, etc.) half-life, and a few other criteria.
//Could still be improved more in the future, but not to bad right now.
void suggestNuclides( PeakToNuclideMatch &answer,
                      std::shared_ptr<const PeakDef> peak,
                      std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > all_peaks,
                      std::shared_ptr<const SpecUtils::Measurement> data,
                      std::shared_ptr<const DetectorPeakResponse> response );

  
//minDetectableCounts(): assumes no peak was detected, and that the entire
//  contribution to data is background.  We will return the number of counts,
//  from which we would be able to detect a peak 95% of the time, assuming
//  simple coutning statistics (for now).
//See eqns 5.52 through 5.56 in Practical Gamma-ray Spectrometry, by Gilmore.
double minDetectableCounts( double energy, double det_sigma, std::shared_ptr<const SpecUtils::Measurement> data );
double minDetectableCounts( std::shared_ptr<const PeakDef> peak, std::shared_ptr<const SpecUtils::Measurement> data );
  

//fractionDetectedWeight(...): used by suggestNuclides() to determine the
//  fraction of peaks detected, that would be expected for a given nuclide,
//  detector, observed data, and shielding configuration.
double fractionDetectedWeight( const std::vector<SandiaDecay::EnergyRatePair> &source_gammas,  //normalization doesnt matter
                           std::shared_ptr<const DetectorPeakResponse> response,
                           double shielding_an,
                           double shielding_ad,
                           std::shared_ptr<const SpecUtils::Measurement> data,
                           std::shared_ptr<const PeakDef> test_peak,
                           std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > all_peaks
                          );

  
//populateCandidateNuclides(...): populates candidate nuclides for 'peak' into
//  'candidates' in a form similar to "B133m 234.73 keV", and then posts
//  'doupdate' to WServer::post() using the session id 'sessionid', so this
//  way 'doupdate' will be executed in the main event loop.  As a concequency
//  of this, one of the bound arguments of 'doupdate' should be another
//  shared pointer for 'candidates'.
//Currently doesnt deal with single or double escape peaks.
void populateCandidateNuclides( std::shared_ptr<const SpecUtils::Measurement> data,
                               std::shared_ptr<const PeakDef> peak,
                               std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > hintpeaks,
                               std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > userpeaks,
                               const std::vector<ReferenceLineInfo> showingRefLines,
                               std::shared_ptr<const DetectorPeakResponse> detector,
                               const std::string sessionid,
                               std::shared_ptr< std::vector<std::string> > candidates,
                               boost::function<void(void)> doupdate );

  
//isotopesFromOtherPeaks(): looks at peaks other than 'peak' and places
//  nuclides, xrays, or reactions that have a gamma line near 'peak', into
//  'otherpeaknucs' in the format similar to "B133m 234.73 keV".
//Called from populateCandidateNuclides(...).
void isotopesFromOtherPeaks( std::vector<std::string> &otherpeaknucs,
                            std::shared_ptr<const PeakDef> peak,
                            std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > allpeaks );

  
//findCharacteristics(): Looks for gamma lines near 'peak' in
//  "data/PhotoPeak.lis", and inserts the relevant xray/nuclide/reaction into
//  'characteristicnucs' in the format similar to "B133m 234.73 keV".
//Called from populateCandidateNuclides(...).
void findCharacteristics( std::vector<std::string> &characteristicnucs,
                          std::shared_ptr<const PeakDef> peak );

  
//findCandidates(): calls suggestNuclides(...), and places its top 10 results
//  into 'suggestednucs' in the format similar to "B133m 234.73 keV".
//Called from populateCandidateNuclides(...).
void findCandidates( std::vector<std::string> &suggestednucs,
                     std::shared_ptr<const PeakDef> peak,
                     std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > allpeaks,
                     std::shared_ptr<const DetectorPeakResponse> detector,
                     std::shared_ptr<const SpecUtils::Measurement> data );

} //namespace IsotopeId

#endif //IsotopeId_h
