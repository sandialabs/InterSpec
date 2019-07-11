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

#include <set>
#include <map>
#include <deque>
#include <vector>
#include <sstream>

#include <Wt/WServer>
#include <Wt/WApplication>

#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/MassAttenuationTool.h"
#include "SpecUtils/SpectrumDataStructs.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergyModel.h"


#include "SandiaDecay/SandiaDecay.h"


using namespace Wt;
using namespace std;

namespace
{
  struct Sorter
  {
    const vector<double> &m_energies;
    IsotopeSearchByEnergyModel::Column m_column;
    Wt::SortOrder m_order;
    
    Sorter( const vector<double> &energies, int column, Wt::SortOrder order )
    : m_energies( energies ),
    m_column( IsotopeSearchByEnergyModel::Column(column) ),
    m_order( order )
    {}
    
    bool operator()( const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &lhsin,
                    const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &rhsin )
    {
      const bool asscending = (m_order==AscendingOrder);
      
      assert( lhsin.size()==rhsin.size() && lhsin.size()==m_energies.size() );
      
      if( m_energies.empty() )
        return asscending;
      
      const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &lhs = (asscending ? lhsin : rhsin);
      const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &rhs = (asscending ? rhsin : lhsin);
      
      switch( m_column )
      {
        case IsotopeSearchByEnergyModel::ParentIsotope:
        case IsotopeSearchByEnergyModel::SpecificIsotope:
          return (lhs[0].m_displayData[m_column] < rhs[0].m_displayData[m_column] );
          
        case IsotopeSearchByEnergyModel::Distance:
          return (lhs[0].m_distance < rhs[0].m_distance);
          
        case IsotopeSearchByEnergyModel::Energy:
        {
          return (std::stod(lhs[0].m_displayData[m_column].narrow())
                  < std::stod(rhs[0].m_displayData[m_column].narrow()));
        }
          
        case IsotopeSearchByEnergyModel::ProfileDistance:
          return (lhs[0].m_profileDistance < rhs[0].m_profileDistance);
          
        case IsotopeSearchByEnergyModel::BranchRatio:
        {
          double lhsval = 0.0, rhsval = 0.0;
          for( size_t i = 0; i < lhs.size(); ++i )
            lhsval += std::stod(lhs[i].m_displayData[m_column].narrow());
          for( size_t i = 0; i < rhs.size(); ++i )
            rhsval += std::stod(rhs[i].m_displayData[m_column].narrow());
          return lhsval < rhsval;
        }
          
        case IsotopeSearchByEnergyModel::ParentHalfLife:
          if( lhs[0].m_nuclide && rhs[0].m_nuclide )
            return (lhs[0].m_nuclide->halfLife < rhs[0].m_nuclide->halfLife);
          return (lhs[0].m_nuclide != NULL);
          
        case IsotopeSearchByEnergyModel::AssumedAge:
          return (lhs[0].m_age < rhs[0].m_age);
          
        case IsotopeSearchByEnergyModel::NumColumns:
          return false;
      }//switch( col )
      
      return false;
    }//bool operator()
  };//struct Sorter
  
  
  /** Conviencince typedef from (pointer to nuclide) to the set of energies of
     that nuclide
   */
  typedef map<const SandiaDecay::Nuclide *, set<double> > NuclideMatches;
  
  /** Returns nuclides filtered such that they have gammas/x-rays in the
      specified energy ranges, with at least the minimum branching ratio and
      half-lives specified.
   */
  NuclideMatches filter_nuclides( const double minbr,
                                  const double minHalfLife,
                                  const vector<double> &energies,
                                  const vector<double> &windows )
  {
    using namespace SandiaDecay;
    
    NuclideMatches filteredNuclides;
    
    //Get isotopes with gammas in all ranges...
    EnergyToNuclideServer::setLowerLimits( minHalfLife, minbr );
    
    auto nucnuc = EnergyToNuclideServer::energyToNuclide();
    if( !nucnuc )
      throw runtime_error( "Couldnt get EnergyToNuclideServer" );
    
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const float minenergy = static_cast<float>(energies[i] - windows[i]);
      const float maxenergy = static_cast<float>(energies[i] + windows[i]);
      
      const bool canBeAnnih = (510.99891f>=minenergy && 510.99891f<=maxenergy);
      
      EnergyToNuclideServer::EnergyNuclidePair enPair( minenergy, NULL );
      vector<EnergyToNuclideServer::EnergyNuclidePair>::const_iterator begin, end, pos;
      
      //nucnuc should be sorted by energy (EnergyNuclidePair::operator<)
      begin = lower_bound( nucnuc->begin(), nucnuc->end(), enPair );
      enPair.energy = maxenergy;
      end = upper_bound( begin, nucnuc->end(), enPair );
      for( pos = begin; pos != end; ++pos )
      {
        //nucnuc actually only contians the nuclides that they themselves give
        //  off the requested energies, meaning we have to go through and inspect
        //  all the nuclides that could decay to the nuclide in nucnuc to see
        //  if they are compatible with the search criteria.  For this we
        //  need to
        int stop = 0;
        double transbr = 1.0;
        for( const SandiaDecay::Transition *t : pos->nuclide->decaysToChildren )
        {
          for( const SandiaDecay::RadParticle &r : t->products )
          {
            if( r.type==SandiaDecay::GammaParticle && fabs(r.energy - pos->energy) < 0.001 )
            {
              transbr = t->branchRatio * r.intensity;
              stop = 1;
            }else if( canBeAnnih && r.type==SandiaDecay::PositronParticle )
            {
              //XXX - there is a slight issue with nuclides that have positrons and gammas near 511 keV, but whatever
              if( !stop )
                transbr = 0.0;
              transbr += 2.0 * t->branchRatio * r.intensity;
              stop = 2;
            }
          }
          
          if( stop==1 )
            break;
        }//for( const SandiaDecay::Transition *t : pos->nuclide->decaysToChildren )
        
        
        const vector<const Nuclide *> forebearers = pos->nuclide->forebearers();
        for( const Nuclide *nuc : forebearers )
        {
          if( nuc->halfLife < minHalfLife )
            continue;
          
          if( minbr<=0.0 || nuc->branchRatioToDecendant(pos->nuclide)*transbr>minbr )
            filteredNuclides[nuc].insert( energies[i] );
        }
      }
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    return filteredNuclides;
  }//filter_nuclides( )

  
  std::shared_ptr<const PeakDef> nearest_peak( const float energy,
                                              const std::vector<std::shared_ptr<const PeakDef>> &allpeaks )
  {
    std::shared_ptr<const PeakDef> nearest;
    double minDE = std::numeric_limits<double>::infinity();
    
    for( const auto &peak : allpeaks )
    {
      const double dE = fabs( peak->mean() - energy );
      if( (dE < minDE)
         && ((energy > peak->lowerX()) && (energy < peak->upperX())) )
      {
        minDE = dE;
        nearest = peak;
      }//if( dE < minDE )
    }//for( int row = 0; row < nrow; ++row )
    
    return nearest;
  }//nearest_peak(...)
  
  
  /** A first attempt to guess which nuclide most closely matches the searched
      energy ranges, using the "profile" of the spectrum, represented by the
      peaks fit (both aurtomated and user) - this is work in progress, and
      there is still very much room to come up with somethign better.
   */
  double profile_weight( std::shared_ptr<const DetectorPeakResponse> detector,
                       const std::shared_ptr<const Measurement> displayed_measurement,
                       const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                       const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                       std::vector<SandiaDecay::EnergyRatePair> srcgammas,
                       const vector<double> &energies,
                       const vector<double> &windows,
                       const IsotopeSearchByEnergyModel::IsotopeMatch &nucmatches,
                       double shielding_an,
                       double shielding_ad )
  {
    /* This function takes a kinda complex approach to seeing if the expected
       profile of a source nuclide matches the observed data - something much
       simpler could probably do just as well or maybe even better (but its been
       a while since I've gotten to think about physics, so we'll try the
       slightly more complicated approach first, especially since its not *that*
       complicated)
     
     [For shielding passed in] A rough overview of this functions is:
     -Add in effects of DRF and shielding for passed in gamma lines, then
      normalize sum of gamma line rates to 1.0.
     -Combine all peaks
     -Match source lines to peaks they may contribute towards
     -Create a sum of weighted gamma intensities for the source that are located
      within peaks.  The weight for each gamma line is determined comparing the
      expected yeild for each peak with the observed yeild.  The expected yeild
      is determined by normaling the gamma rate relative to the nearest search
      energy peak count rate (min detectable counts is used if no peak fit).
     -As a penalty: Sum the fraction of source lines that should have created
      a peak, but did not (but only if the expected peak area, determined
      similar to previous bullet, is larger than the minimum detectable area)
     
     See notes at the end of this function for additional improvements.
     
     See fractionDetectedWeight(...) function in IsotopeId.h for a similarish function
     
     */
    
    auto print_debug = [&]() -> bool {
      return nucmatches.m_displayData[IsotopeSearchByEnergyModel::Column::ParentIsotope].toUTF8() == "Fe59";
      //return nucmatches.m_nuclide->symbol && nucmatches.m_nuclide->symbol == "Fe59";
      //return false;
    };
    
    const string source_name = print_debug() ? nucmatches.m_displayData[IsotopeSearchByEnergyModel::Column::ParentIsotope].toUTF8() : "";
    
    //Combine user peaks with the automated search peaks
    const vector<shared_ptr<const PeakDef>> allpeaks
      = [&]() -> vector<shared_ptr<const PeakDef>> {
        vector<shared_ptr<const PeakDef>> answer;
        for( const auto &p : user_peaks )
        {
          if( p )
            answer.push_back( p );
        }
    
        //Add automated search peaks that arent already represented by the user peaks
        for( const auto &peak : automated_search_peaks ){
          if( !peak )
            continue;
          
          auto nearpeak = nearest_peak( peak->mean(), allpeaks );
          const double peak_sigma = peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth();
          if( !nearpeak || (fabs(nearpeak->mean() - peak->mean()) > peak_sigma) )
            answer.push_back( peak );
        }
        return answer;
      }();
    
    if( print_debug() )
    {
      cout << "Have peaks: ";
      for( auto i : allpeaks )
        cout << "{" << i->mean() << "keV,amp=" << (i->gausPeak() ? i->amplitude() : -1.0) << "},";
      cout << endl;
    }//if( print_debug() )
    
    
    if( srcgammas.empty() )
    {
      cerr << "IsotopeSearchByEnergyModel profile_weight: no source lines" << endl;
      return -1.0;
    }
    
    //Scale the yeilds for the detector response function and shielding specified
    for( auto &src : srcgammas )
    {
      const double det_sf = (!!detector ? detector->intrinsicEfficiency(src.energy) : 1.0f);
      const double xs = MassAttenuation::massAttenuationCoeficient( shielding_an, src.energy );
      const double shielding_sf = exp( -shielding_ad * xs );
      
      src.numPerSecond *= det_sf * shielding_sf;
    }//for( auto &src : srcgammas )
    
    
    const double min_energy = [&]() -> double {
        double answer = 80.0;  //80keV is arbitrary.  IF we are searching x-rays, perhaps we should allow down to 40 keV or so
        for( size_t i = 0; i < energies.size() && i < windows.size(); ++i )
          answer = std::min( answer, energies[i] - windows[i] );
        return answer;
      }();
    
    //Normalize expected gammas to sum to 1.0 above min_energy.
    const double srcgamma_sum = [&]() -> double {
      double sum = 0.0;
      for( const auto &src : srcgammas )
        if( src.energy >= min_energy )
          sum += src.numPerSecond;
      return sum;
    }();
    
    if( srcgamma_sum <= DBL_EPSILON )
      return -1.0;
    
    for( auto &src : srcgammas )
      src.numPerSecond /= srcgamma_sum;
    
    
    //Pair up search energie/windows with the fraction of gammas in that range and peak area
    vector<tuple<double,double,double,double>> search_energies;
    for( size_t i = 0; i < energies.size() && i < windows.size(); ++i )
    {
      const double energy = energies[i], w = windows[i];
      double frac_gammas_in_range = 0.0;
      for( auto &src : srcgammas )
      {
        if( (src.energy > (energy-w)) && (src.energy < (energy+w)) )
          frac_gammas_in_range += src.numPerSecond;
      }//for( auto &src : srcgammas )
      
      double peak_area = 0.0;
      for( const auto &p : allpeaks )
      {
        if( (p->mean() > (energy-w)) && (p->mean() < (energy+w)) )
          peak_area += p->gausPeak() ? p->amplitude() : p->areaFromData(displayed_measurement);
      }
      
      //If no peaks, use somethign like minimum detectable peak area
      if( peak_area < 0.1 && displayed_measurement )
        peak_area = 2.33*sqrt(displayed_measurement->gamma_integral(energy-w, energy+w));
      
      if( frac_gammas_in_range > 0.0 )
      {
        
        if( print_debug() )
          cout << "For " << source_name << " SearchEnergyRange "
               << energies[i] << "keV has frac_gammas_in_range="
               << frac_gammas_in_range << ", and peak_area=" << peak_area << endl;
        search_energies.emplace_back( energies[i], windows[i], frac_gammas_in_range, peak_area );
      }
    }//for( size_t i = 0; i < energies.size() && i < windows.size(); ++i )
    
    std::sort( begin(search_energies), end(search_energies),
              [](const tuple<double,double,double,double> &a,
                 const tuple<double,double,double,double> &b) -> bool
              { return get<0>(a) < get<0>(b); } );
    
    
    if( search_energies.size() == 0 )
    {
      if( print_debug() )
        cerr << "IsotopeSearchByEnergyModel: no search energies specified" << endl;
      return -1.0;
    }
    
    //Go through and see which lines can be accounted for with peaks in the spectrum.
    //  Lines in the search areas will always be considered to be in peaks.
    
    
    double in_peaks_frac = 0.0, not_in_peaks_frac = 0.0;
    vector<SandiaDecay::EnergyRatePair> out_of_peak_src_lines;
    std::map<std::shared_ptr<const PeakDef>,double> src_counts_in_peaks;
    
    for( auto &src : srcgammas )
    {
      if( src.energy < min_energy )
        continue;
      
      //Check if this srb energy is in a peak.
      //  Use Mean +- 0.75 FWHM as range (arbitrarily chosen)
      bool in_peak = false;
      auto peak = nearest_peak(src.energy, allpeaks );
      if( peak && (src.energy > (peak->mean() - 0.75*peak->fwhm()))
          && (src.energy < (peak->mean() + 0.75*peak->fwhm())) )
      {
        in_peak = true;
        //I think std::map requires its template type (e.g., the double) to be
        //  default constructable, which I think implies that the default
        //  constructor (e.g., value set to 0.0) will be called for the double
        //  in src_counts_in_peaks[peak] += src.numPerSecond;, but I'll be
        //  safe/faster anyway.
        auto pos = src_counts_in_peaks.find(peak);
        if( pos == end(src_counts_in_peaks) )
          src_counts_in_peaks.insert( make_pair(peak,src.numPerSecond) );
        else
          pos->second += src.numPerSecond;
        if( print_debug() )
          cout << "For " << source_name << " Peak " << peak->mean() << "keV "
               << " is contributted to by " << src.numPerSecond << endl;
        in_peaks_frac += src.numPerSecond;
      }//if( src.energy contributes to thie peak )
      
  
      bool in_search_range = false;
      for( const auto &range : search_energies )
      {
        const double energy = get<0>(range), window = get<1>(range);
        if( (src.energy > (energy - window)) && (src.energy < (energy + window)) )
        {
          in_search_range = true;
          break;
        }
      }//for( const auto &range : search_energies )
      
      if( !in_search_range && !in_peak )
      {
        out_of_peak_src_lines.push_back( src );
        not_in_peaks_frac += src.numPerSecond;
      }
    }//for( auto &src : srcgammas )
    
    
    double weighted_in_peak_frac = 0.0;
    
    for( const auto &peak_srcamp : src_counts_in_peaks )
    {
      const auto peak = peak_srcamp.first;
      const double peak_mean = peak->mean();
      const double peak_area = peak->amplitude();
      const double peak_range_rate = peak_srcamp.second;
      
      if( (peak_mean < min_energy) || (peak_area < DBL_EPSILON) || (peak_range_rate < DBL_EPSILON) )
        continue;  //there could be some edge effect here, but whatever for now.
      
      const tuple<double,double,double,double> &nearest_search
        = *std::min_element( begin(search_energies), end(search_energies),
          [=]( const tuple<double,double,double,double> &lhs, const tuple<double,double,double,double> &rhs ) -> bool {
            const double lhsdist = std::fabs( get<0>(lhs) - peak->mean() );
            const double rhsdist = std::fabs( get<0>(rhs) - peak->mean() );
            return lhsdist < rhsdist;
          } );
      
      const double nearest_search_energy = get<0>(nearest_search);
      const double nearest_search_window = get<1>(nearest_search);
      
      if( peak_mean > (nearest_search_energy-nearest_search_window)
         && peak_mean < (nearest_search_energy+nearest_search_window) )
      {
        weighted_in_peak_frac += peak_range_rate;
        continue;
      }
      
      const double search_range_rate = get<2>(nearest_search);
      if( search_range_rate <= DBL_EPSILON ) //this should never happen, or else we messed up before callign this function, but JIC
      {
//#if( PERFORM_DEVELOPER_CHECKS )
//        char buffer[256];
//        snprintf( buffer, sizeof(buffer), "Found a search_range_rate=%f for nuclide=%s - should not happen.",
//                  search_range_rate, source_name.c_str() );
//        log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
//#endif
//        continue;
        return -1.0;
      }//if( search_range_rate <= 0.0 )
      
      const double search_range_peak_area = get<3>(nearest_search);
      const double expected_peak_area = peak_range_rate * search_range_peak_area / search_range_rate;
      
      if( expected_peak_area > peak_area )
      {
        //There could be another nuclide adding to things, so only
        //  assymptotically punish up to a factor of 0.05 (arbitrarily chosen)
        weighted_in_peak_frac += peak_range_rate * std::max(0.05,1.0/sqrt(expected_peak_area/peak_area));
        
        if( print_debug() )
          cout << "For " << source_name << " Peak " << peak->mean() << "keV (smaller, expected="
               << expected_peak_area << " vs peak_area=" << peak_area << ") "
               << "contributing " << (peak_range_rate * std::max(0.05,1.0/sqrt(expected_peak_area/peak_area)))
               << endl;
      }else
      {
        //We arent seeing as many gammas here as expected, punish up to 0.2 (arbitrary)
        weighted_in_peak_frac += peak_range_rate * std::max( 0.05, expected_peak_area/peak_area );
        
        if( print_debug() )
          cout << "For " << source_name << " Peak " << peak->mean() << "keV (larger, expected="
               << expected_peak_area << " vs peak_area=" << peak_area << ") "
               << "contributing " << (peak_range_rate * std::max(0.05, expected_peak_area/peak_area))
               << endl;
      }
      if( print_debug() )
        cout << source_name << " at peak energy=" << peak->mean()
             << " gives expected_peak_area=" << expected_peak_area
             << ", and peak_area=" << peak_area << endl;
    }//for( const auto &peak_srcamp : src_counts_in_peaks )
    
    //Now punish for not having peaks you should probably have
    double weighted_out_of_peak_frac = 0.0;
    for( const auto &src : out_of_peak_src_lines )
    {
      if( src.energy < min_energy )
        continue;
      
      if( !displayed_measurement )
      {
        weighted_out_of_peak_frac += src.numPerSecond;
        continue;
      }
      
      const tuple<double,double,double,double> &nearest_search
      = *std::min_element( begin(search_energies), end(search_energies),
                          [&]( const tuple<double,double,double,double> &lhs, const tuple<double,double,double,double> &rhs ) -> bool {
                            const double lhsdist = std::fabs( get<0>(lhs) - src.energy );
                            const double rhsdist = std::fabs( get<0>(rhs) - src.energy );
                            return lhsdist < rhsdist;
                          } );
      
      //const double nearest_search_energy = get<0>(nearest_search);
      const double nearest_search_window = get<1>(nearest_search);
      const double search_range_rate = get<2>(nearest_search);
      const double search_range_peak_area = get<3>(nearest_search);
      
      const double peak_sigma = ((detector && detector->hasResolutionInfo())
                                  ? detector->peakResolutionSigma(src.energy) : nearest_search_window);
      
      const double background = gamma_integral( displayed_measurement, src.energy-3.0*peak_sigma, src.energy+3.0*peak_sigma );
      const double mindetcounts = 2.33 * sqrt(background);
      
      const double expected = src.numPerSecond * search_range_peak_area / search_range_rate;
      if( expected > mindetcounts )
      {
        if( print_debug() )
          cout << "For " << source_name << " out_of_peak_src_line " << src.energy << "keV "
               << "contributing " << src.numPerSecond << " (expected=" << expected << " vs mindetcounts=" << mindetcounts << ")" << endl;
        weighted_out_of_peak_frac += src.numPerSecond;  //could do a weighting similar to in-peak area
      }else
      {
        if( print_debug() )
          cout << "For " << source_name << " out_of_peak_src_line " << src.energy << "keV "
               << "is insignificant (expected=" << expected << " vs mindetcount=" << mindetcounts
                << "), not contributing" << endl;
      }
    }//for( const auto &src : out_of_peak_src_lines )
    
    
    if( print_debug() )
    {
      cout << source_name << " weighted_in_peak_frac=" << weighted_in_peak_frac << endl;
      cout << source_name << " weighted_out_of_peak_frac=" << weighted_out_of_peak_frac << endl;
    }
    
    
    //maybe make a penalty if the observed peaks in the search regions would
    //  show negative shielding or something.
    
    //Maybe reward this nuclide for being responsible for more peaks.  This
    //  would counteract nuclides that have just the one energy that you are
    //  searching for, even though there may be other nuclides that are
    //  responsible for many peaks.  It looks like something like +0.02 for
    //  every peak accounted for would be reasonably powerful for this.
    
    //Could reward for the search ranges having a higher BR of gammas in them
    
    //Right now src lines are included if they are within 0.75*FWHM of the peak
    //  mean.  We could weight lines according to how close they are to the peak
    //  mean; or even tighten up the 0.75 to like 0.5. This would be
    //  particularly useful for NaI detectors.  I think somewhere I had used
    //  a small constant plus fraction of FWHM, and that worked well; maybe look
    //  for that.
    
    //Should now check half-lives of parents and decendants and bias things
    //  towards preffering more likely to be seen isotopes (e.g., filter really
    //  short time periods out)
    
    //Maybe not in this function, but where value is returned, could preffer
    //  unshielded solutions.
    
    
    
    //Currently performs poorly on NaI for Np237
    
    return weighted_in_peak_frac - weighted_out_of_peak_frac;
  }//double profile_weight(...)

}//namespace


IsotopeSearchByEnergyModel::IsotopeMatch::IsotopeMatch()
: m_distance(0.0), m_age(0.0), m_branchRatio(0.0), m_profileDistance(-1.0),
m_nuclide(0), m_transition(0), m_particle(0),
m_sourceGammaType(PeakDef::NormalGamma),
m_element(0), m_xray(0), m_reaction(0)
{
}


IsotopeSearchByEnergyModel::IsotopeSearchByEnergyModel( Wt::WObject *parent )
: WAbstractItemModel( parent ),
  m_sortColumn( IsotopeSearchByEnergyModel::Column::ProfileDistance ),
  m_sortOrder( Wt::AscendingOrder )
{
  
}//IsotopeSearchByEnergyModel( constuctor )


IsotopeSearchByEnergyModel::~IsotopeSearchByEnergyModel()
{
  
}//IsotopeSearchByEnergyModel destructor


void IsotopeSearchByEnergyModel::nuclidesWithAllEnergies(
            const IsotopeSearchByEnergyModel::NucToEnergiesMap &filteredNuclides,
            const vector<double> &energies,
            const vector<double> &windows,
            const double minBR,
            const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
            const std::shared_ptr<const Measurement> displayed_measurement,
            const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
            const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
            std::vector< vector<IsotopeSearchByEnergyModel::IsotopeMatch> > &answer )
{
  if( energies.empty() )
    return;
  
  char buffer[32];
  
  //XXX - currently only taking 'most likely' combination of matching of
  //      energies - should implement getting all permutations!
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  for( const NucToEnergiesMap::value_type &nm : filteredNuclides )
  {
    //chack to see if this nuclide has gammas for each energy, not strictly
    //  necassarry, but probably computationally faster (do we care about this
    //  here though)
    bool hasAll = true;
    vector<size_t> energy_xray_is_for;
    vector<const SandiaDecay::EnergyIntensityPair *> matching_xrays;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double energy = energies[i];
      const double de = fabs( windows[i] );
      set<double>::const_iterator lb, up, iter;
      lb = nm.second.lower_bound( energy - de );
      up = nm.second.upper_bound( energy + de );
      bool found = false;
      for( iter = lb; iter != up; ++iter )
        found |= (fabs((*iter)-energy) <= de);
      
      
      if( !found && (energy < 115.0*PhysicalUnits::keV) )
      {//lets look through the x-rays for this
        double minxraydist = 999.9;
        SandiaDecay::EnergyIntensityPair closest_xray(0.0,0.0);
        const SandiaDecay::Element *el = db->element( nm.first->atomicNumber );
        const vector<SandiaDecay::EnergyIntensityPair> &xrays = el->xrays;
        for( size_t j = 0; j < xrays.size(); ++j )
        {
          const double dist = fabs(xrays[j].energy-energy);
          if( dist <= minxraydist )
          {
            minxraydist = dist;
            closest_xray = xrays[j];
          }//if( dist <= minxraydist )
          
          if( minxraydist <= de )
          {
            found = true;
            energy_xray_is_for.push_back( i );
            matching_xrays.push_back( &(xrays[j]) );
          }//if( minxraydist <= de )
        }//for( loop iver x-rays )
      }//if( !found && (energy < 115.0*PhysicalUnits::keV) )
      
      if( !found )
      {
        hasAll = false;
        break;
      }//if( !found )
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( !hasAll )
      continue;
    
    double dist = 0.0;
    vector<IsotopeMatch> nucmatches;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      vector<size_t>::const_iterator xraypos = find( energy_xray_is_for.begin(),
                                                    energy_xray_is_for.end(), i );
      
      if( xraypos != energy_xray_is_for.end() )
      {
        const size_t pos = xraypos - energy_xray_is_for.begin();
        const double energy = energies[i];
        const SandiaDecay::EnergyIntensityPair *xray = matching_xrays[pos];
        
        IsotopeMatch match;
        match.m_distance = fabs(energy - xray->energy);    //sum of distance over all energies searched
        dist += match.m_distance;
        
        match.m_age = 0.0;         //age assumed for listing things
        match.m_branchRatio = xray->intensity;
        
        //Only one of the following will be valid: m_nuclide, m_element, m_reaction
        match.m_nuclide = NULL;
        match.m_transition = NULL;
        match.m_particle = NULL;
        match.m_element = db->element( nm.first->atomicNumber );
        match.m_xray = xray;
        match.m_reaction = NULL;
        //        match.m_reactionEnergy;
        
        //        match.m_displayData[ParentHalfLife] = "";
        match.m_displayData[AssumedAge] = PhysicalUnits::printToBestTimeUnits( 0.0 );
        //        match.m_displayData[SpecificIsotope] = "";
        if( match.m_element )
          match.m_displayData[ParentIsotope] = match.m_element->symbol;
        
        if( xray )
        {
          snprintf( buffer, sizeof(buffer), "%.2f", xray->energy );
          match.m_displayData[Energy] = buffer;
        }//if( xray )
        
        snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
        match.m_displayData[Distance] = buffer;
        
        snprintf( buffer, sizeof(buffer), "%.2g", match.m_branchRatio );
        match.m_displayData[BranchRatio] = buffer;
        
        nucmatches.push_back( match );
      }else//if( xraypos != energy_xray_is_for.end() )
      {
        size_t trans_index = 0;
        const SandiaDecay::Transition *transition = NULL;
        IsotopeMatch match;
        
        PeakDef::SourceGammaType sourceGammaType = PeakDef::NormalGamma;
        
        PeakDef::findNearestPhotopeak( nm.first, energies[i],
                                      windows[i], transition, trans_index, sourceGammaType );
        if( !transition && (sourceGammaType!=PeakDef::AnnihilationGamma) )
          continue;
        
        match.m_nuclide = nm.first;
        match.m_transition = transition;
        if( transition )
          match.m_particle = &(transition->products[trans_index]);
        match.m_sourceGammaType = sourceGammaType;
        
        match.m_distance = 99999999.9;
        switch( sourceGammaType )
        {
          case PeakDef::NormalGamma:
          case PeakDef::XrayGamma:
            if( match.m_particle )
              match.m_distance = fabs(energies[i] - match.m_particle->energy);
            break;
            
          case PeakDef::AnnihilationGamma:
            match.m_distance = fabs(energies[i] - 510.99891*SandiaDecay::keV );
            break;
            
          case PeakDef::SingleEscapeGamma:
            if( match.m_particle )
              match.m_distance = fabs(energies[i] - (match.m_particle->energy - 510.99891));
            break;
            
          case PeakDef::DoubleEscapeGamma:
            if( match.m_particle )
              match.m_distance = fabs(energies[i] - (match.m_particle->energy - 2.0*510.99891) );
            break;
        }//switch( sourceGammaType )
        
        match.m_age = PeakDef::defaultDecayTime( nm.first );
        
        SandiaDecay::NuclideMixture mixture;
        mixture.addNuclide( SandiaDecay::NuclideActivityPair(nm.first,1.0) );
        const vector<SandiaDecay::EnergyRatePair> gammas
        = mixture.gammas( match.m_age,
                         SandiaDecay::NuclideMixture::OrderByAbundance, true );
        
        double nearestEnergy = 999999.9, nearestAbun = 0.0, maxAbund = -999.9;
        for( const SandiaDecay::EnergyRatePair &aep : gammas )
        {
          double d = 999999.9;
          
          if( match.m_sourceGammaType == PeakDef::AnnihilationGamma )
            d = fabs( aep.energy - 510.99891*SandiaDecay::keV );
          else if( match.m_particle )
            d = fabs( aep.energy - match.m_particle->energy );
          
          if( d < nearestEnergy )
          {
            nearestEnergy = d;
            nearestAbun = aep.numPerSecond;
          }//if( d < nearestEnergy )
          
          maxAbund = std::max( maxAbund, aep.numPerSecond );
        }//for( const SandiaDecay::AbundanceEnergyPair &aep : gammas )
        
        match.m_branchRatio = nearestAbun / maxAbund;
        
        if( match.m_branchRatio < minBR )
          continue;
        
        match.m_displayData[ParentIsotope] = match.m_nuclide->symbol;
        
        if( match.m_sourceGammaType == PeakDef::AnnihilationGamma )
          snprintf( buffer, sizeof(buffer), "510.99" );
        else if( match.m_particle )
          snprintf( buffer, sizeof(buffer), "%.2f", match.m_particle->energy );
        
        match.m_displayData[Energy] = buffer;
        
        snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
        match.m_displayData[Distance] = buffer;
        
        snprintf( buffer, sizeof(buffer), "%.2f", match.m_branchRatio );
        match.m_displayData[BranchRatio] = buffer;
        
        stringstream trnsitionstrm;
        if( match.m_transition && match.m_transition->parent && match.m_transition->child )
        {
          trnsitionstrm << match.m_transition->parent->symbol << "&rarr;"
          << match.m_transition->child->symbol;
        }else if( match.m_transition && match.m_transition->parent )
        {
          using namespace SandiaDecay;
          trnsitionstrm << match.m_transition->mode
          << " of " << match.m_transition->parent->symbol;
        }else if( !match.m_transition )
        {
          trnsitionstrm << "Annih. Gamma";
        }//if( match.m_transition->parent... ) / else
        
        if( match.m_transition && match.m_particle->type == SandiaDecay::XrayParticle )
          trnsitionstrm << " xray";
        
        match.m_displayData[SpecificIsotope] = trnsitionstrm.str();
        
        if( !i )
        {
          match.m_displayData[ParentHalfLife]
              = PhysicalUnits::printToBestTimeUnits(match.m_nuclide->halfLife);
          match.m_displayData[AssumedAge]
              = PhysicalUnits::printToBestTimeUnits(match.m_age);
        }//if( !i )
        
        dist += match.m_distance;
        
        nucmatches.push_back( match );
      }//if( xraypos != energy_xray_is_for.end() ) / else
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( nucmatches.size() != energies.size() )
      continue;
    
    if( !nucmatches[0].m_nuclide )
    {
      for( IsotopeMatch &match : nucmatches )
      {
        if( match.m_nuclide )
        {
          nucmatches[0].m_nuclide = match.m_nuclide;
          nucmatches[0].m_displayData[ParentIsotope] = match.m_nuclide->symbol;
          nucmatches[0].m_displayData[Energy]
                = nucmatches[0].m_displayData[Energy].narrow() + " (xray)";
          break;
        }//if( match.m_nuclide )
      }//for( IsotopeMatch &match : nmagicucmatches )
    }//if( !nucmatches[0].m_nuclide )
    
    nucmatches[0].m_distance = dist;
    
    const double gcm2 = PhysicalUnits::g / PhysicalUnits::cm2;
    const double atomic_nums[]   = { 1.0, 26.0, 74.0 };
    const double areal_density[] = { 0.0*gcm2, 10.0*gcm2, 25.0*gcm2 };
    static_assert( sizeof(atomic_nums) == 3*sizeof(atomic_nums[0]), "" );
    static_assert( sizeof(areal_density) == 3*sizeof(areal_density[0]), "" );
    
    double mw = -999.9;
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nucmatches[0].m_nuclide, 0.001*SandiaDecay::curie );
    vector<SandiaDecay::EnergyRatePair> srcgammas = mix.photons( nucmatches[0].m_age );
    
    for( size_t i = 0; i < 3; ++i )
    {
      const double weight = profile_weight( detector_response_function,
                                    displayed_measurement,
                                    user_peaks,
                                    automated_search_peaks, srcgammas,
                                    energies, windows, nucmatches[0],
                                    atomic_nums[i], areal_density[i] );
      mw = std::max(mw,weight);
    }
    nucmatches[0].m_profileDistance = mw;
    snprintf( buffer, sizeof(buffer), "%.2f", mw );
    nucmatches[0].m_displayData[ProfileDistance] = buffer;
  
    
    snprintf( buffer, sizeof(buffer), "%.2f", dist );
    nucmatches[0].m_displayData[Distance] = buffer;
    answer.push_back( nucmatches );
  }//for( const NuclideMatches::value_type &nm : filteredNuclides )
}//void nuclidesWithAllEnergies


void IsotopeSearchByEnergyModel::xraysWithAllEnergies(
                                                      const std::vector<double> &energies,
                                                      const std::vector<double> &windows,
                                                      const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
                                                      const std::shared_ptr<const Measurement> displayed_measurement,
                                                      const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                                      const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                                      SearchResults &answer )
{
  if( energies.empty() )
    return;
  
  char buffer[32];
  
  for( size_t i = 0; i < energies.size(); ++i )
    if( (energies[i]-windows[i]) > 120*PhysicalUnits::keV )
      return;
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const vector<const SandiaDecay::Element *> &elements = db->elements();
  
  for( const SandiaDecay::Element *el : elements )
  {
    vector<IsotopeMatch> nucmatches;
    vector<const SandiaDecay::EnergyIntensityPair *> xray_matches;
    const vector<SandiaDecay::EnergyIntensityPair> &xrays = el->xrays;
    if( xrays.size() < energies.size() )
      continue;
    
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double energy = energies[i];
      const double de = fabs( windows[i] );
      double minEnergy = 1000.0 * de;
      const SandiaDecay::EnergyIntensityPair *xray = NULL;
      
      for( const SandiaDecay::EnergyIntensityPair &e : xrays )
      {
        const double delta = fabs((e.energy)-energy);
        if( delta < minEnergy )
        {
          minEnergy = delta;
          xray = &e;
        }//if( delta < minEnergy )
      }//for( const SandiaDecay::EnergyIntensityPair &e : xrays )
      
      if( minEnergy < de )
        xray_matches.push_back( xray );
      else
        break;
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( xray_matches.size() != energies.size() )
      continue;
    
    double dist = 0.0;
    for( size_t i = 0; i < xray_matches.size(); ++i )
    {
      const double energy = energies[i];
      const SandiaDecay::EnergyIntensityPair *xray = xray_matches[i];
      
      IsotopeMatch match;
      match.m_distance = fabs(energy - xray->energy);    //sum of distance over all energies searched
      dist += match.m_distance;
      
      match.m_age = 0.0;         //age assumed for listing things
      match.m_branchRatio = xray->intensity;
      
      //Only one of the following will be valid: m_nuclide, m_element, m_reaction
      match.m_nuclide = NULL;
      match.m_transition = NULL;
      match.m_particle = NULL;
      match.m_element = el;
      match.m_xray = xray;
      match.m_reaction = NULL;
      match.m_displayData[AssumedAge] = PhysicalUnits::printToBestTimeUnits( 0.0 );
      match.m_displayData[ParentIsotope] = match.m_element->symbol;
      
      snprintf( buffer, sizeof(buffer), "%.2f", xray->energy );
      match.m_displayData[Energy] = buffer;
      
      snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
      match.m_displayData[Distance] = buffer;
      
      snprintf( buffer, sizeof(buffer), "%.2g", match.m_branchRatio );
      match.m_displayData[BranchRatio] = buffer;
      
      nucmatches.push_back( match );
    }//for( const SandiaDecay::EnergyIntensityPair *xray : xray_matches )
    
    nucmatches[0].m_distance = dist;
    
    {//begin code to get profile distance
      vector<SandiaDecay::EnergyRatePair> srcxrays;
      for( const auto &x : xrays )
        srcxrays.emplace_back( x.intensity, x.energy );
      
      nucmatches[0].m_profileDistance
         = profile_weight( detector_response_function, displayed_measurement,
                           user_peaks, automated_search_peaks, srcxrays,
                           energies, windows, nucmatches[0], 1.0, 0.0);
      
      snprintf( buffer, sizeof(buffer), "%.2f", nucmatches[0].m_profileDistance );
      nucmatches[0].m_displayData[ProfileDistance] = buffer;
    }//end code to get profile distance
    
    snprintf( buffer, sizeof(buffer), "%.2f", dist );
    nucmatches[0].m_displayData[Distance] = buffer;
    
    answer.push_back( nucmatches );
  }//for( const SandiaDecay::Element *el : elements )
}//void xraysWithAllEnergies(...)



void IsotopeSearchByEnergyModel::reactionsWithAllEnergies(
                                                          const std::vector<double> &energies,
                                                          const std::vector<double> &windows,
                                                          const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
                                                          const std::shared_ptr<const Measurement> displayed_measurement,
                                                          const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                                          const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                                          SearchResults &answer )
{
  if( energies.empty() )
    return;
  
  char buffer[32];
  
  const ReactionGamma *db = ReactionGammaServer::database();
  set<const ReactionGamma::Reaction *> reactions;
  
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const float energy = static_cast<float>(energies[i]);
    const float de = static_cast<float>( fabs( windows[i] ) );
    vector<const ReactionGamma::Reaction *> thesereaction;
    db->reactions( energy-de, energy+de, thesereaction );
    
    if( i == 0 )
    {
      reactions.insert( thesereaction.begin(), thesereaction.end() );
    }else
    {
      set<const ReactionGamma::Reaction *> surviving_reactions;
      
      for( const ReactionGamma::Reaction *rctn : thesereaction )
      {
        if( reactions.count( rctn ) )
          surviving_reactions.insert( rctn );
      }//for( const ReactionGamma::Reaction *rctn : thesereaction )
      
      reactions.swap( surviving_reactions );
      
      if( reactions.empty() )
        return;
    }//if( i == 0 ) / else
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  
  for( const ReactionGamma::Reaction *rctn : reactions )
  {
    double dist = 0.0;
    vector<IsotopeMatch> matches;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double energy = energies[i];
      double smallestDelta = 999999999.9;
      ReactionGamma::EnergyAbundance nearesteA;
      
      for( const ReactionGamma::EnergyAbundance &ea : rctn->gammas )
      {
        const double delta = fabs( ea.energy - energy );
        if( delta < smallestDelta  )
        {
          nearesteA = ea;
          smallestDelta = delta;
        }//if( delta < nearestEnergy  )
      }//for( const ReactionGamma::EnergyAbundance &ea : rctn->gammas )
      
      IsotopeMatch match;
      match.m_distance = energy - nearesteA.energy;
      dist += fabs(match.m_distance);
      match.m_age = 0.0;
      match.m_reaction = rctn;
      match.m_branchRatio = nearesteA.abundance;
      match.m_reactionEnergy = nearesteA;
      
      match.m_displayData[AssumedAge] = PhysicalUnits::printToBestTimeUnits( 0.0 );
      match.m_displayData[ParentIsotope] = rctn->name();
      
      snprintf( buffer, sizeof(buffer), "%.2f", nearesteA.energy );
      match.m_displayData[Energy] = buffer;
      
      //      snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
      //      match.m_displayData[Distance] = buffer;
      
      snprintf( buffer, sizeof(buffer), "%.2g", match.m_branchRatio );
      match.m_displayData[BranchRatio] = buffer;
      
      //      if( rctn->targetNuclide )
      //        match.m_displayData[SpecificIsotope] = rctn->targetNuclide->symbol;
      
      matches.push_back( match );
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    matches[0].m_distance = dist;
    
    {//begin code to get profile distance
      vector<SandiaDecay::EnergyRatePair> srcgammas;
      for( const ReactionGamma::EnergyAbundance &ea : rctn->gammas )
        srcgammas.emplace_back( ea.abundance, ea.energy );
      
      matches[0].m_profileDistance
        = profile_weight( detector_response_function, displayed_measurement,
                       user_peaks, automated_search_peaks, srcgammas,
                       energies, windows, matches[0], 1.0, 0.0);
      
      snprintf( buffer, sizeof(buffer), "%.2f", matches[0].m_profileDistance );
      matches[0].m_displayData[ProfileDistance] = buffer;
    }//end code to get profile distance
    
    snprintf( buffer, sizeof(buffer), "%.2g", matches[0].m_distance );
    matches[0].m_displayData[Distance] = buffer;
    
    answer.push_back( matches );
  }//for( const ReactionGamma::Reaction *rctn : reactions )
  
}//void reactionsWithAllEnergies(...)



void IsotopeSearchByEnergyModel::clearResults()
{
  beginRemoveRows( WModelIndex(), 0, rowCount()-1 );
  m_matches.clear();
  endRemoveRows();
}//void clearResults();


void IsotopeSearchByEnergyModel::updateSearchResults(
                                                     std::shared_ptr<IsotopeSearchByEnergyModel::SearchWorkingSpace> workingspace )
{
  const vector<double> &energies = workingspace->energies;
  const vector<double> &windows = workingspace->windows;
  
  vector< vector<IsotopeMatch> > &matches = workingspace->matches;
  
  m_sortColumn = workingspace->sortColumn;
  m_sortOrder = workingspace->sortOrder;
  
  if( m_matches.size() )
  {
    beginRemoveRows( WModelIndex(), 0, rowCount()-1 );
    m_matches.clear();
    endRemoveRows();
  }//if( m_matches.size() )
  
  if( matches.size() )
  {
    const int ninsert = static_cast<int>( matches.size() * energies.size() );
    beginInsertRows( WModelIndex(), 0, ninsert - 1 );
    m_windows = windows;
    m_energies = energies;
    m_matches.swap( matches );
    endInsertRows();
  }//if( matches.size() )
  
  workingspace->searchdoneCallback();
  
  wApp->triggerUpdate();
}//void updateSearchResults()



void IsotopeSearchByEnergyModel::setSearchEnergies(
                                                   std::shared_ptr<SearchWorkingSpace> workingspace,
                                                   const double minbr,
                                                   const double minHalfLife,
                                                   Wt::WFlags<IsotopeSearchByEnergyModel::RadSource> srcs,
                                                   const std::string appid,
                                                   boost::function< void(void) > updatefcn )
{
  if( !workingspace )
    throw runtime_error( "setSearchEnergies(...): invalid workingspace" );
  
  const vector<double> &energies = workingspace->energies;
  const vector<double> &windows = workingspace->windows;
  
  vector< vector<IsotopeMatch> > &matches = workingspace->matches;
  matches.clear();
  
  if( energies.size() != windows.size() )
    throw runtime_error( "setSearchEnergies(...): input error" );
  
  if( energies.empty() )
  {
    WServer::instance()->post(  appid, updatefcn );
    return;
  }//if( energies.empty() )
  
  
  using SandiaDecay::Element;
  using SandiaDecay::Nuclide;
  using SandiaDecay::EnergyIntensityPair;
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const vector<const SandiaDecay::Element *> &elements = db->elements();
  
  //Get x-rays with at least one of the energies
  map<const SandiaDecay::Element *, vector<EnergyIntensityPair> > filteredXrays;
  for( const Element *el : elements )
    for( const EnergyIntensityPair &xray : el->xrays )
      filteredXrays[el].push_back( xray );
  
  //Get reactions with at least one of the energies
  const ReactionGamma *rctnDb = ReactionGammaServer::database();
  vector<ReactionGamma::ReactionPhotopeak> reactions;
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const float minenergy = static_cast<float>(energies[i] - windows[i]);
    const float maxenergy = static_cast<float>(energies[i] + windows[i]);
    rctnDb->reactions( minenergy, maxenergy, reactions );
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  
  //Time to make all the pairings
  auto &user_peaks = workingspace->user_peaks;
  std::vector<std::shared_ptr<const PeakDef>> auto_peaks = workingspace->automated_search_peaks;
  
  if( auto_peaks.empty() && workingspace->foreground && workingspace->displayed_measurement )
  {
    //iOS/Android may not have auto-search peaks yet.  Also recently loaded
    //  spectra and it looks like sometimes when previous states were loaded
    const auto data = workingspace->displayed_measurement;
    auto userpeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>( begin(user_peaks), end(user_peaks) );
    const bool singleThreaded = false;
    auto_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks( data, userpeaksdeque, singleThreaded );
    
    
    const std::set<int> samplenums = workingspace->foreground_samplenums;
    auto autopeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>( begin(auto_peaks), end(auto_peaks) );
    workingspace->foreground->setAutomatedSearchPeaks( samplenums, autopeaksdeque );
  }//
  
  const auto &meas = workingspace->displayed_measurement;
  
  auto &drf = workingspace->detector_response_function;
  
  //Nuclides that match all energies
  if( srcs & kGamma )
  {
    const auto filteredNuclides = filter_nuclides( minbr, minHalfLife, energies, windows );
    nuclidesWithAllEnergies( filteredNuclides, energies, windows, minbr, drf, meas, user_peaks, auto_peaks, matches );
  }//if( srcs & kGamma )
  
  //Get elements with x-rays which match all energies
  if( srcs & kXRay )
    xraysWithAllEnergies( energies, windows, drf, meas, user_peaks, auto_peaks, matches );
  
  //Get elements with reactions which match all energies
  if( srcs & kReaction )
    reactionsWithAllEnergies( energies, windows, drf, meas, user_peaks, auto_peaks, matches );
  
  //Get elements with gamma+xrays which match all energies
  
  //Get elements with xrays+reactions which match all energies
  
  //Get elements with gamma+reactions which match all energies
  
  //sort the data
  sortData( matches, energies, workingspace->sortColumn, workingspace->sortOrder );
  
  WServer::instance()->post(  appid, updatefcn );
}//void setSearchEnergies( const vector<double> &energies, const double window )


int IsotopeSearchByEnergyModel::columnCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return NumColumns;
}//int columnCount( const Wt::WModelIndex &parent ) const


int IsotopeSearchByEnergyModel::rowCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return static_cast<int>( m_matches.size() * m_energies.size() );
}//int rowCount( const WModelIndex &parent ) const


WModelIndex IsotopeSearchByEnergyModel::parent( const WModelIndex & ) const
{
  return WModelIndex();
}//WModelIndex parent( const Wt::WModelIndex &index ) const;


const SandiaDecay::Nuclide *IsotopeSearchByEnergyModel::nuclide(
                                                                const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_nuclide;
}//const SandiaDecay::Nuclide *nuclide( const Wt::WModelIndex &index ) const


const SandiaDecay::Element *IsotopeSearchByEnergyModel::xrayElement(
                                                                    const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_element;
}//xrayElement( const Wt::WModelIndex &index ) const


const ReactionGamma::Reaction *IsotopeSearchByEnergyModel::reaction(
                                                                    const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_reaction;
}//reaction( const Wt::WModelIndex &index ) const


double IsotopeSearchByEnergyModel::assumedAge(
                                              const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return -1.0;
  
  return m_matches[matchNum].at(0).m_age;
}//double assumedAge( const Wt::WModelIndex &index ) const;


boost::any IsotopeSearchByEnergyModel::data( const WModelIndex &index,
                                            int role ) const
{
  const int row = index.row();
  const Column col = Column( index.column() );
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  const size_t energyNum = static_cast<size_t>( row % m_energies.size() );
  
  
  if( (role!=Wt::DisplayRole) || (matchNum>m_matches.size())
     || (col>=NumColumns) )
    return boost::any();
  
  const vector<IsotopeMatch> &match = m_matches[matchNum];
  const IsotopeMatch &iso = match[energyNum];
  
  
  switch( col )
  {
    case ParentIsotope:
    case Distance:
      if( energyNum )
        return boost::any();
      //fallthrough intentional (my first intentional use in like 5 years)
    case Energy: case BranchRatio: case ProfileDistance:
    case SpecificIsotope: case ParentHalfLife: case AssumedAge:
      return iso.m_displayData[col];
      break;
      
    case NumColumns:
      break;
  }//switch( col )
  
  return boost::any();
}//boost::any data(...)


boost::any IsotopeSearchByEnergyModel::headerData( int section,
                                                  Wt::Orientation orientation,
                                                  int role ) const
{
  if( orientation==Wt::Horizontal && role==Wt::LevelRole )
    return 0;
  
  if (role==Wt::DisplayRole)
  {
    switch( section )
    {
      case IsotopeSearchByEnergyModel::Column::ParentIsotope:
        return WString("Parent");
      case IsotopeSearchByEnergyModel::Column::Distance:
        return WString("Diff.");
      case IsotopeSearchByEnergyModel::Column::Energy:
        return WString("Energy (keV)");
      case IsotopeSearchByEnergyModel::Column::BranchRatio:
        return WString("B.R.");
      case IsotopeSearchByEnergyModel::Column::ProfileDistance:
        return WString("Profile");
      case IsotopeSearchByEnergyModel::Column::SpecificIsotope:
        return WString("Decay");
      case IsotopeSearchByEnergyModel::Column::ParentHalfLife:
        return WString("Parent H.L.");
      case IsotopeSearchByEnergyModel::Column::AssumedAge:
        return WString("Assumed Age");
      case IsotopeSearchByEnergyModel::Column::NumColumns:
        break;
    }//switch( col )
  }//DisplayRole
  else if (role==Wt::ToolTipRole)
  {
    switch( section )
    {
      case IsotopeSearchByEnergyModel::Column::ParentIsotope:
        return WString("Parent nuclide");
      case IsotopeSearchByEnergyModel::Column::Distance:
        return WString("Difference between selected nuclide's energy level and searched energy level");
      case IsotopeSearchByEnergyModel::Column::Energy:
        return WString("Nuclide energy");
      case IsotopeSearchByEnergyModel::Column::BranchRatio:
        return WString("Branching ratio of selected nuclide");
      case IsotopeSearchByEnergyModel::Column::ProfileDistance:
        return WString("A rough metric for how close the observed spectrum comes to having the expected peaks for the selected source");
      case IsotopeSearchByEnergyModel::Column::SpecificIsotope:
        return boost::any();
      case IsotopeSearchByEnergyModel::Column::ParentHalfLife:
        return WString("Parent half life");
      case IsotopeSearchByEnergyModel::Column::AssumedAge:
        return WString("Assumed age of nuclide");
      case IsotopeSearchByEnergyModel::Column::NumColumns:
        break;
    }//switch( col )
  }//ToolTipRole
  
  return boost::any();
}//boost::any headerData( int section, Orientation orientation, int role ) const


WModelIndex IsotopeSearchByEnergyModel::index( int row, int column,
                                              const WModelIndex &parent ) const
{
  if( (column>=NumColumns) || (row>=rowCount()) )
    return WModelIndex();
  
  void *ptr = (void *)&(m_matches[row]);  //ah, whatever
  return createIndex( row, column, ptr );
}//WModelIndex index( int row, int column, const WModelIndex &parent ) const



void IsotopeSearchByEnergyModel::sort( int column, Wt::SortOrder order )
{
  layoutAboutToBeChanged().emit();
  m_sortOrder = order;
  m_sortColumn = IsotopeSearchByEnergyModel::Column(column);
  sortData( m_matches, m_energies, column, order );
  layoutChanged().emit();
}//void sort(...)


WFlags<ItemFlag> IsotopeSearchByEnergyModel::flags( const WModelIndex &index ) const
{
  const Column col = Column( index.column() );
  const int level = (index.row() % m_energies.size());
  
  if( (col==ParentIsotope) && (level==0) )
    return WFlags<ItemFlag>( ItemIsSelectable | ItemIsXHTMLText );
  
  return WFlags<ItemFlag>( ItemIsXHTMLText );
}//WFlags<ItemFlag> flags( const Wt::WModelIndex &index ) const


void IsotopeSearchByEnergyModel::sortData( vector< vector<IsotopeMatch> > &data,
                                          const vector<double> &energies,
                                          int column, Wt::SortOrder order )
{
  std::stable_sort( data.begin(), data.end(), Sorter(energies, column, order) );
}//sortData(...)

