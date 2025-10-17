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
#include <deque>
#include <set>
#include <algorithm>
#include <limits>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AnalystChecks.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/RelActCalcAuto.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include <Wt/WApplication>

using namespace std;

namespace AnalystChecks
{
  std::string hello_world()
  {
    return "Hello World from AnalystChecks namespace!";
  }
  
  DetectedPeakStatus detected_peaks(const DetectedPeaksOptions& options, InterSpec* interspec)
  {
    if (!interspec) {
      throw std::runtime_error("No InterSpec session available");
    }
    
    // Get session ID from InterSpec (using Wt application)
    string user_session = "direct_call";
    if (Wt::WApplication* app = Wt::WApplication::instance()) {
      user_session = app->sessionId();
    }
    
    std::shared_ptr<SpecMeas> meas = interspec->measurment(options.specType);
    if (!meas) {
      throw std::runtime_error("No measurement loaded for " 
          + std::string(SpecUtils::descriptionText(options.specType)) 
          + " spectrum");
    }

    const set<int> sample_nums = interspec->displayedSamples(options.specType);
    if (sample_nums.empty()) {
      throw std::runtime_error("No samples displayed for " 
          + std::string(SpecUtils::descriptionText(options.specType)) 
          + " spectrum");
    }

    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(options.specType);
    if (!spectrum) {
      throw std::runtime_error("No spectrum displayed for " 
          + std::string(SpecUtils::descriptionText(options.specType)) 
          + " spectrum");
    }
    
    shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_peaks = meas->automatedSearchPeaks(sample_nums);
    shared_ptr<const deque<shared_ptr<const PeakDef>>> user_peaks = meas->peaks(sample_nums);
    
    if (!auto_peaks) {
      // Search for peaks
      const bool singleThreaded = false;
      const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
      const auto det = meas->detector();
      const vector<shared_ptr<const PeakDef>> found_auto_peaks
                = ExperimentalAutomatedPeakSearch::search_for_peaks(spectrum, det, user_peaks, singleThreaded, isHPGe);
        
      auto autopeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>(begin(found_auto_peaks), end(found_auto_peaks));
      auto_peaks = autopeaksdeque;
      meas->setAutomatedSearchPeaks(sample_nums, autopeaksdeque);
    }

    vector<shared_ptr<const PeakDef>> all_peaks;

    if (user_peaks)
      all_peaks.insert(all_peaks.end(), user_peaks->begin(), user_peaks->end());

    // We will add auto-search peaks only if there isnt already a user peak at essentially the same energy
    if (auto_peaks) {
      //Lambda to find the nearest peak, so far, to a given energy.
      auto nearest_peak = [&all_peaks](const float energy) -> std::shared_ptr<const PeakDef> {
        std::shared_ptr<const PeakDef> nearest;
        double minDE = std::numeric_limits<double>::infinity();

        for (const auto &peak : all_peaks) {
          const double dE = fabs(peak->mean() - energy);
          if ((dE < minDE)
              && ((energy > peak->lowerX()) && (energy < peak->upperX()))) {
            minDE = dE;
            nearest = peak;
          }
        }

        return nearest;
      };

      for (const shared_ptr<const PeakDef> &peak : *auto_peaks) {
        auto nearpeak = nearest_peak(peak->mean());
        const double peak_sigma = peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth();
        if (!nearpeak || (fabs(nearpeak->mean() - peak->mean()) > peak_sigma))
          all_peaks.push_back(peak);
      }
    }
    
    std::sort(begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMeanShrdPtr);
    
    DetectedPeakStatus result;
    result.userSession = user_session;
    result.peaks = std::move(all_peaks);
    
    return result;
  }
  
  
  FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec )
  {
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");
    
    // Get session ID from InterSpec (using Wt application)
    string user_session = "direct_call";
    if (Wt::WApplication* app = Wt::WApplication::instance()) {
      user_session = app->sessionId();
    }
    
    std::shared_ptr<SpecMeas> meas = interspec->measurment(options.specType);
    if (!meas) {
      throw std::runtime_error("No measurement loaded for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    const set<int> sample_nums = interspec->displayedSamples(options.specType);
    if (sample_nums.empty()) {
      throw std::runtime_error("No samples displayed for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(options.specType);
    if (!data) {
      throw std::runtime_error("No spectrum displayed for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    shared_ptr<const DetectorPeakResponse> det = meas->detector();
    
    shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks(sample_nums);
    assert( meas_peaks );
    if( !meas_peaks )
      throw runtime_error( "Unexpected error getting existing peak list" );
    
    
    const bool isHPGe = PeakFitUtils::is_likely_high_res( interspec );
    
    vector<shared_ptr<const PeakDef>> origPeaks;
    
    for( const shared_ptr<const PeakDef> &p : *meas_peaks )
    {
      //Avoid fitting a peak in the same area a data defined peak is.
      if( p->gausPeak() )
        origPeaks.push_back( p );
      else if( (options.energy >= p->lowerX()) && (options.energy <= p->upperX()) )
        throw runtime_error( "Can not fit a peak within the ROI of a data-defined peak." );
    }//if( pmodel->peaks() )
    
    double pixelPerKev = -1.0; //This triggers an "automed" peak fit, which has higher thresholds for keeping peak.
    pair<vector<shared_ptr<const PeakDef>>, vector<shared_ptr<const PeakDef>>> foundPeaks;
    foundPeaks = searchForPeakFromUser( options.energy, pixelPerKev, data, origPeaks, det, isHPGe );
    
    vector<shared_ptr<const PeakDef>> &peaks_to_add_in = foundPeaks.first;
    const vector<shared_ptr<const PeakDef>> &peaks_to_remove = foundPeaks.second;
    
    if( peaks_to_add_in.empty() )
      throw runtime_error( "Could not fit peak." );
    
    // The new peak may share a ROI with previously existing peaks, so we have to figure out
    //  what peak we just fit.
    //  TODO: For the moment, we'll just match up to nearest peak, but in a normal user double click for a peak, we match up peaks in ROI to previous peaks, and then take whats left-over, which is maybe a better way.
    shared_ptr<const PeakDef> fit_peak;
    for( shared_ptr<const PeakDef> peak : peaks_to_add_in )
    {
      if( peak->parentNuclide() || peak->xrayElement() || peak->reaction() )
        continue;  //We know new peak doesnt have a source assigned yet
      
      const double diff = fabs(peak->mean() - options.energy);
      if( !fit_peak || (diff < fabs(fit_peak->mean() - options.energy)) )
        fit_peak = peak;
    }//for( shared_ptr<const PeakDef> peak : peaks_to_add_in )
    
    
    assert( fit_peak );
    if( !fit_peak )
      throw runtime_error( "Could not identify newly fit peak." );
    
    // get fit peak
    
    if( options.source.has_value() && !options.source.value().empty() )
    {
      const string source = options.source.value();
      
      shared_ptr<PeakDef> new_fit_peak = make_shared<PeakDef>( *fit_peak );
      
      PeakModel::SetGammaSource res = PeakModel::setNuclideXrayReaction( *new_fit_peak, source, 4.0 );
      
      if( res == PeakModel::SetGammaSource::NoSourceChange )
        throw runtime_error( "Peak was fit, but source string '" + source + "' was not valid to set source." );
      
      bool replaced_peak = false;
      for( size_t index = 0; index < peaks_to_add_in.size(); ++index )
      {
        if( peaks_to_add_in[index] == fit_peak )
        {
          replaced_peak = true;
          fit_peak = new_fit_peak;
          peaks_to_add_in[index] = new_fit_peak;
          break;
        }
      }
      assert( replaced_peak );
      
      vector<Wt::WColor> src_colors;
      vector<shared_ptr<const PeakDef>> other_src_peaks;
      for( const shared_ptr<const PeakDef> &srcPeak : origPeaks )
      {
        if( (srcPeak->parentNuclide() || srcPeak->xrayElement() || srcPeak->reaction())
            && srcPeak->parentNuclide()==new_fit_peak->parentNuclide()
            && srcPeak->xrayElement()==new_fit_peak->xrayElement()
            && srcPeak->reaction()==new_fit_peak->reaction() )
        {
          other_src_peaks.push_back( srcPeak );
          
          if( !srcPeak->lineColor().isDefault()
             && (std::find( begin(src_colors), end(src_colors),srcPeak->lineColor()) == end(src_colors)) )
          {
            src_colors.push_back( srcPeak->lineColor() );
          }
        }
      }//for( const auto &p : *m_peaks )
          
      if( src_colors.size() == 1 )
      {
        new_fit_peak->setLineColor( src_colors[0] );
      }else if( other_src_peaks.empty() )
      {
        ReferencePhotopeakDisplay *ref_lines = interspec->referenceLinesWidget();
        if( ref_lines )
        {
          string src_name;
          if( new_fit_peak->parentNuclide() )
            src_name = new_fit_peak->parentNuclide()->symbol;
          else if( new_fit_peak->xrayElement() )
            src_name = new_fit_peak->xrayElement()->symbol;
          else if( new_fit_peak->reaction() )
            src_name = new_fit_peak->reaction()->name();
          
          Wt::WColor c = ref_lines->suggestColorForSource( src_name );
          new_fit_peak->setLineColor( c );
        }
      }
    }//if( options.source.has_value() )
    
    PeakModel *pmodel = interspec->peakModel();
    assert( pmodel );
    if( options.addToUsersPeaks && pmodel )
    {
      if( options.specType == SpecUtils::SpectrumType::Foreground )
      {
        pmodel->removePeaks( peaks_to_remove );
        pmodel->addPeaks( peaks_to_add_in );
      }else
      {
        // It would also be valid to handle the Foreground in this logic block, I think
        shared_ptr<deque<shared_ptr<const PeakDef>>> new_deque = make_shared<deque<shared_ptr<const PeakDef>>>(*meas_peaks);
        deque<shared_ptr<const PeakDef>> &peak_deque = *new_deque;
        
        for( shared_ptr<const PeakDef> old : peaks_to_remove )
        {
          auto pos = std::find( begin(peak_deque), end(peak_deque), old );
          assert( pos != end(peak_deque) );
          if( pos != end(peak_deque) )
            peak_deque.erase( pos );
        }
        
        for( shared_ptr<const PeakDef> new_peak : peaks_to_add_in )
          peak_deque.push_back( new_peak );
        std::sort( begin(peak_deque), end(peak_deque), &PeakDef::lessThanByMeanShrdPtr );
        
        pmodel->setPeaks( peak_deque, options.specType );
      }
    }//if( options.addToUsersPeaks )
    
    
    
    FitPeakStatus result;
    
    result.userSession = user_session;
    result.fitPeak = fit_peak;
    result.peaksInRoi = peaks_to_add_in;
    
    return result;
  }//FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec );
  
  
  GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec )
  {
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");
    
    // Get session ID from InterSpec (using Wt application)
    string user_session = "direct_call";
    if (Wt::WApplication* app = Wt::WApplication::instance()) {
      user_session = app->sessionId();
    }
    
    std::shared_ptr<SpecMeas> meas = interspec->measurment(options.specType);
    if (!meas) {
      throw std::runtime_error("No measurement loaded for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    const set<int> sample_nums = interspec->displayedSamples(options.specType);
    if (sample_nums.empty()) {
      throw std::runtime_error("No samples displayed for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks(sample_nums);
    assert( meas_peaks );
    if( !meas_peaks )
      throw runtime_error( "Unexpected error getting existing peak list" );


    GetUserPeakStatus status;
    status.userSession = user_session;

    if( meas_peaks )
    {
      // Filter peaks by energy range if specified
      const bool has_lower = options.lowerEnergy.has_value();
      const bool has_upper = options.upperEnergy.has_value();

      // Validate energy range if both bounds are specified
      if( has_lower && has_upper && (*options.upperEnergy < *options.lowerEnergy) )
      {
        throw std::runtime_error( "upperEnergy (" + std::to_string(*options.upperEnergy)
                                  + " keV) must be greater than or equal to lowerEnergy ("
                                  + std::to_string(*options.lowerEnergy) + " keV)" );
      }

      if( !has_lower && !has_upper )
      {
        // No filtering, return all peaks
        status.peaks.insert( end(status.peaks), begin(*meas_peaks), end(*meas_peaks) );
      }else
      {
        const double lower_energy = has_lower ? *options.lowerEnergy : -std::numeric_limits<double>::infinity();
        const double upper_energy = has_upper ? *options.upperEnergy : std::numeric_limits<double>::infinity();

        for( const shared_ptr<const PeakDef> &peak : *meas_peaks )
        {
          if( peak && peak->mean() >= lower_energy && peak->mean() <= upper_energy )
            status.peaks.push_back( peak );
        }
      }
    }

    return status;
  }
  
  vector<tuple<float,float,float>> cacl_estimated_gamma_importance( const vector<tuple<float,float>> &gamma_energies_and_yields )
  {
    float sum_yields_sqrt_energy = 0.0f;
    for( const tuple<float,float> &gamma : gamma_energies_and_yields )
    {
      sum_yields_sqrt_energy += std::get<1>(gamma) * sqrt(std::get<0>(gamma));
    }
    
    vector<tuple<float,float,float>> answer;
    answer.reserve(gamma_energies_and_yields.size());
    for( const tuple<float,float> &gamma : gamma_energies_and_yields )
    {
      float importance = std::get<1>(gamma) * sqrt(std::get<0>(gamma)) / sum_yields_sqrt_energy;
      answer.push_back({std::get<0>(gamma), std::get<1>(gamma), importance});
    }
    return answer;
  }
  
  std::vector<float> get_characteristic_gammas( const std::string &nuclide )
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    if( !db )
      throw std::runtime_error("No decay database available");

    const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide );
    if( !nuc )
    {
      // See if its a element for nuclide fluorescence
      const SandiaDecay::Element * const el = db->element( nuclide );
      if( !el )
      {
        if( SpecUtils::istarts_with( nuclide, "Neut" ) )
          return {478.0f, 847.0f, 2223.3f, 2235.3f};

        // See if its a reaction
        try
        {
          const ReactionGamma * const rctnDb = ReactionGammaServer::database();
          vector<ReactionGamma::ReactionPhotopeak> reactions;
          
          if( rctnDb )
            rctnDb->gammas( nuclide, reactions );

          if( reactions.empty() )
            throw std::runtime_error("No reactions found for '" + nuclide + "'");

          vector<tuple<float,float>> energies_and_yields;
          for( const auto &reaction : reactions )
            energies_and_yields.push_back({reaction.energy, reaction.abundance});

          vector<tuple<float,float,float>> importances = cacl_estimated_gamma_importance(energies_and_yields);
          std::sort( begin(importances), end(importances), [](const tuple<float,float,float> &a, const tuple<float,float,float> &b) {
              return std::get<2>(a) > std::get<2>(b);
          });

          const float max_importance = std::get<2>(importances[0]);
          vector<float> answer;
          for( const auto &importance : importances )
          {
            if( std::get<2>(importance) < max_importance * 0.25 )
              break;
            answer.push_back( std::get<0>(importance) );
          }
        

          if( answer.size() > 3 )
            answer.resize(3);

          return answer;
        }catch(...)
        {
          throw std::runtime_error("No nuclide, element, or reaction found for name '" + nuclide + "'");
        }//try / catch
      }//if( !el )

      // For uranium and lead we'll just hard-code the result
      if( el->atomicNumber == 92 ) //Uranium
        return {98.4f};
      if( el->atomicNumber == 82 ) //Lead
        return {75.0f};

      // We'll return the highest "importance" x-ray for the element
      vector<tuple<float,float>> gamma_energies_and_yields;
      for( const auto &xray : el->xrays )
        gamma_energies_and_yields.push_back({xray.energy, xray.intensity});

      vector<tuple<float,float,float>> importances = cacl_estimated_gamma_importance(gamma_energies_and_yields);
      std::sort( begin(importances), end(importances), [](const tuple<float,float,float> &a, const tuple<float,float,float> &b) {
        return std::get<2>(a) > std::get<2>(b);
      });

      if( importances.empty() )
        return {};

      return {std::get<0>(importances[0])};
    }//if( !nuc )
    
    // TODO: should we use `photopeak_lis()`?
    
    const vector<tuple<string, const SandiaDecay::Nuclide *, float>> &nuc_characteristics = IsotopeId::characteristicGammas();
    
    vector<float> answer;
    for( const auto &nuc_characteristic : nuc_characteristics )
    {
      if( get<0>(nuc_characteristic) == nuclide )
        answer.push_back( get<2>(nuc_characteristic) );
    }

    if( answer.empty() )
    {
      // We'll return up to a few of the highest "importance" gammas for the element.
      const double age = PeakDef::defaultDecayTime( nuc );
      SandiaDecay::NuclideMixture mixture;
      mixture.addNuclide( SandiaDecay::NuclideActivityPair(nuc,1.0E4) );
      const vector<SandiaDecay::EnergyRatePair> photons
                                        = mixture.photons(age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy);

      vector<tuple<float,float>> gamma_energies_and_yields;
      for( const auto &photon : photons )
        gamma_energies_and_yields.push_back({photon.energy, photon.numPerSecond});

      vector<tuple<float,float,float>> importances = cacl_estimated_gamma_importance(gamma_energies_and_yields);
      std::sort( begin(importances), end(importances), [](const tuple<float,float,float> &a, const tuple<float,float,float> &b) {
        return std::get<2>(a) > std::get<2>(b);
      });

      if( importances.empty() )
        return {};

      const float max_importance = std::get<2>(importances[0]);
      for( const auto &importance : importances )
      {
        if( std::get<2>(importance) < max_importance * 0.25 )
          break;
        answer.push_back( std::get<0>(importance) );
      }

      if( answer.size() > 3 )
        answer.resize(3);
    }//if( answer.empty() )

    return answer;
  }//std::vector<float> get_characteristic_gammas( const std::string &nuclide )
  
  
  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_nuclides_with_characteristics_in_energy_range( double lower_energy, double upper_energy, InterSpec *interspec  )
  {
    std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>> answer;
    
    if( lower_energy > upper_energy )
      std::swap( lower_energy, upper_energy );
    
    // First, we'll identify all the common nuclides, elements, and reactions that are 
    //  typically encountered.
    set<const SandiaDecay::Nuclide *> nuclides;
    set<const SandiaDecay::Element *> elements;
    set<const ReactionGamma::Reaction *> reactions;
    
    std::vector<IsotopeSearchByEnergy::NucSearchCategory> nuc_categories;
    if( interspec && interspec->nuclideSearch() )
      nuc_categories = interspec->nuclideSearch()->search_categories();
    else
      IsotopeSearchByEnergy::init_category_info( nuc_categories );

    for( const IsotopeSearchByEnergy::NucSearchCategory &category : nuc_categories )
    {
      nuclides.insert( begin(category.m_specific_nuclides), end(category.m_specific_nuclides) );
      elements.insert( begin(category.m_specific_elements), end(category.m_specific_elements) );
      reactions.insert( begin(category.m_specific_reactions), end(category.m_specific_reactions) );
    }

    const vector<tuple<string, const SandiaDecay::Nuclide *, float>> &nuc_characteristics = IsotopeId::characteristicGammas();
    for( const auto &nuc_characteristic : nuc_characteristics )
    {
      if( get<1>(nuc_characteristic) )
        nuclides.insert( get<1>(nuc_characteristic) );
    }

    // Now add all the nuclides that have characteristics in the energy range
    for( const auto &nuclide : nuclides )
    {
      const vector<float> gammas = get_characteristic_gammas( nuclide->symbol );
      for( const float gamma : gammas )
      {
        if( gamma >= lower_energy && gamma <= upper_energy )
          answer.push_back( nuclide );
      }
    }

    for( const auto &element : elements )
    {
      const vector<float> xrays = get_characteristic_gammas( element->symbol );
      for( const float xray : xrays )
      {
        if( xray >= lower_energy && xray <= upper_energy )
          answer.push_back( element );
      }
    }

    for( const auto &reaction : reactions )
    {
      const vector<float> gammas = get_characteristic_gammas( reaction->name() );
      for( const float gamma : gammas )
      {
        if( gamma >= lower_energy && gamma <= upper_energy )
          answer.push_back( reaction );
      }
    }

    return answer;
  }//get_nuclides_with_characteristics_in_energy_range(...)
  
  
  FitPeaksForNuclideStatus fit_peaks_for_nuclides( const FitPeaksForNuclideOptions &options, InterSpec *interspec )
  {
    // TODO: This function tries a few different RelActAuto configuration, and chooses the one with the best Chi2. It should actually be Chi2/DOF, and then maybe also take into account what peaks are chosen
    // TODO: The order of Ln(x) equation used should be based on how many peaks, and the energy range - should do something more inteligently
    // TODO: Should split range up - e.g. above and below 120 keV, if nuclide has peaks on both sides
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");
    
    // Get session ID from InterSpec (using Wt application)
    string user_session = "direct_call";
    if (Wt::WApplication* app = Wt::WApplication::instance()) {
      user_session = app->sessionId();
    }
    
    FitPeaksForNuclideStatus result;
    
    try {
      // Get the foreground spectrum
      std::shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
      if (!meas) {
        throw std::runtime_error("No foreground measurement available");
      }
      
      std::shared_ptr<const SpecUtils::Measurement> foreground = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
      if (!foreground) {
        throw std::runtime_error("No foreground spectrum available");
      }
      
      // Get background (optional)
      std::shared_ptr<const SpecUtils::Measurement> background = interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
      
      // Get detector response function
      std::shared_ptr<const DetectorPeakResponse> drf = meas->detector();
      
      // Get detected peaks to use as input
      DetectedPeaksOptions peak_options;
      peak_options.specType = SpecUtils::SpectrumType::Foreground;
      DetectedPeakStatus peak_status = detected_peaks(peak_options, interspec);
      
      std::vector<std::shared_ptr<const PeakDef>> all_peaks = peak_status.peaks;
      
      // Get nuclide names from options
      const std::vector<std::string> &nuclide_names = options.nuclides;
      
      if (nuclide_names.empty()) {
        throw std::runtime_error("No nuclides specified");
      }
      
      // Create base nuclide info that will be reused across all attempts
      std::vector<RelActCalcAuto::NucInputInfo> base_nuclides;
      for (const std::string &nuc_name : nuclide_names) {
        RelActCalcAuto::NucInputInfo nuc_info;
        nuc_info.source = RelActCalcAuto::source_from_string(nuc_name);
        
        if (RelActCalcAuto::is_null(nuc_info.source)) {
          throw std::runtime_error("Invalid nuclide name: " + nuc_name);
        }
        
        // Set reasonable default age for nuclides
        if (RelActCalcAuto::nuclide(nuc_info.source)) {
          const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(nuc_info.source);
          nuc_info.age = PeakDef::defaultDecayTime(nuc);
          nuc_info.fit_age = false; // Don't fit age by default
        }
        
        base_nuclides.push_back(nuc_info);
      }
      
      // Set up base ROI to cover the full spectrum energy range
      RelActCalcAuto::RoiRange base_roi;
      base_roi.lower_energy = std::max( 55.0f, foreground->gamma_energy_min() );
      base_roi.upper_energy = foreground->gamma_energy_max();
      base_roi.continuum_type = PeakContinuum::OffsetType::Linear;
      base_roi.force_full_range = false;
      base_roi.allow_expand_for_peak_width = true;
      
      // Create all trial options combinations to run in parallel
      std::vector<RelActCalcAuto::Options> trial_options;
      
      // Define skew types to try
      std::vector<PeakDef::SkewType> skew_types = {
        PeakDef::SkewType::NoSkew,
        PeakDef::SkewType::GaussExp
      };
      
      for (PeakDef::SkewType skew_type : skew_types) {
        // LnX equation form
        {
          RelActCalcAuto::Options solve_options;
          
          // Create relative efficiency curve with LnX
          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "LnX Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
          rel_eff_curve.rel_eff_eqn_order = 3;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;
          
          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois.push_back(base_roi);
          
          // Set other options
          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = RelActCalcAuto::FwhmForm::Gadras;
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;
          
          trial_options.push_back(std::move(solve_options));
        }
        
        // FramPhysicalModel with Aluminum shielding
        {
          RelActCalcAuto::Options solve_options;
          
          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "Physical Model Al Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;
          rel_eff_curve.phys_model_use_hoerl = true;
          
          // Create aluminum external shielding
          auto al_shield = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
          al_shield->atomic_number = 13.0; // Aluminum
          al_shield->fit_atomic_number = false;
          al_shield->areal_density = 0.1; // Starting value in g/cm2
          al_shield->fit_areal_density = true;
          al_shield->lower_fit_areal_density = 0.001;
          al_shield->upper_fit_areal_density = 10.0;
          
          rel_eff_curve.phys_model_external_atten.push_back(al_shield);
          
          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois.push_back(base_roi);
          
          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = RelActCalcAuto::FwhmForm::Gadras; //RelActCalcAuto::FwhmForm::Polynomial_3
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;
          
          trial_options.push_back(std::move(solve_options));
        }
        
        // FramPhysicalModel with Lead shielding
        {
          RelActCalcAuto::Options solve_options;
          
          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "Physical Model Pb Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;
          rel_eff_curve.phys_model_use_hoerl = true;
          
          // Create lead external shielding
          auto pb_shield = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
          pb_shield->atomic_number = 82.0; // Lead
          pb_shield->fit_atomic_number = false;
          pb_shield->areal_density = 0.1; // Starting value in g/cm2
          pb_shield->fit_areal_density = true;
          pb_shield->lower_fit_areal_density = 0.0;
          pb_shield->upper_fit_areal_density = 10.0;
          
          rel_eff_curve.phys_model_external_atten.push_back(pb_shield);
          
          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois.push_back(base_roi);
          
          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = RelActCalcAuto::FwhmForm::Gadras;
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;
          
          trial_options.push_back(std::move(solve_options));
        }
      }
      
      // Execute all trials in parallel using ThreadPool
      std::vector<RelActCalcAuto::RelActAutoSolution> solutions(trial_options.size());
      std::mutex solutions_mutex;
      
      SpecUtilsAsync::ThreadPool pool;
      for (size_t i = 0; i < trial_options.size(); ++i) {
        pool.post([&, i]() {
          try {
            RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
              trial_options[i], foreground, background, drf, all_peaks
            );
            
            std::lock_guard<std::mutex> lock(solutions_mutex);
            solutions[i] = std::move(solution);
          } catch (const std::exception &e) {
            std::lock_guard<std::mutex> lock(solutions_mutex);
            solutions[i].m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
            solutions[i].m_error_message = e.what();
          }
        });
      }
      
      pool.join();
      
      // Find the best solution from all parallel trials
      RelActCalcAuto::RelActAutoSolution best_solution;
      double best_chi2 = std::numeric_limits<double>::max();
      bool found_valid_solution = false;
      string last_error_message;
      
      for (size_t i = 0; i < solutions.size(); ++i) {
        const RelActCalcAuto::RelActAutoSolution &solution = solutions[i];
        
        if (solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success && 
            solution.m_chi2 < best_chi2) {
          best_chi2 = solution.m_chi2;
          best_solution = solution;
          found_valid_solution = true;
        } else if (solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success) {
          last_error_message = solution.m_error_message;
          cerr << "RelActCalcAuto::solve(trial " << i << "): failed: " << solution.m_error_message << endl;
        }
      }
      
      if (!found_valid_solution) {
        throw std::runtime_error("RelActCalcAuto::solve failed for all attempted configurations ('" + last_error_message + "')");
      }
      
      // Use the best solution
      RelActCalcAuto::RelActAutoSolution solution = std::move(best_solution);
      
      // Extract the fit peaks from the solution
      result.fitPeaks.clear();
      for (const PeakDef &peak : solution.m_fit_peaks_in_spectrums_cal ) {  //If we apply energy cal, we would choose solution.m_fit_peaks
        result.fitPeaks.push_back(std::make_shared<PeakDef>(peak));
      }
      
      // If doNotAddPeaksToUserSession is false, add the peaks to the user's session
      if (!options.doNotAddPeaksToUserSession) {
        PeakModel *pmodel = interspec->peakModel();
        if (pmodel) {
          std::vector<std::shared_ptr<const PeakDef>> peaks_to_add;
          for (const auto &peak : result.fitPeaks) {
            peaks_to_add.push_back(peak);
          }
          pmodel->addPeaks(peaks_to_add);
        }
      }
      
    } catch (const std::exception &e) {
      throw std::runtime_error("Error in fit_peaks_for_nuclides: " + string(e.what()));
    }
    
    return result;
  }//FitPeaksForNuclideStatus fit_peaks_for_nuclides( const FitPeaksForNuclideOptions &options, InterSpec *interspec )

  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_characteristics_near_energy( const double energy, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");
    
    shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !meas || !spectrum )
      throw std::runtime_error("No foreground loaded");
    
    
    float fwhm = -1.0;
    
    // Try the currently loaded peak detector response
    if( meas->detector() && meas->detector()->hasResolutionInfo() )
      fwhm = meas->detector()->peakResolutionFWHM( static_cast<float>(energy) );
    
    // Estimate the FWHM response...
    shared_ptr<deque<shared_ptr<const PeakDef>>> peak_deque;
    if( fwhm <= 0.0 )
    {
      try
      {
        DetectedPeaksOptions opts;
        opts.specType = SpecUtils::SpectrumType::Foreground;
        const DetectedPeakStatus peaks = detected_peaks( opts, interspec );
        peak_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(peaks.peaks), end(peaks.peaks) );
        
        DetectorPeakResponse drf;
        drf.setIntrinsicEfficiencyFormula( "1.0", 2.54*PhysicalUnits::cm, PhysicalUnits::keV,
                                               0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarField );
        
        drf.fitResolution( peak_deque, spectrum, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
        
        fwhm = drf.peakResolutionFWHM( static_cast<float>(energy) );
      }catch( std::exception & )
      {
      }
    }//if( fwhm <= 0.0 )
    
    // Find nearest peak, and if within 20% of the energy, use it.
    if( (fwhm <= 0.0) && peak_deque )
    {
      double smallest_energy_diff = std::numeric_limits<double>::max();
      double nearest_fwhm = 0.0;
      for( const auto &peak : *peak_deque )
      {
        double energy_diff = std::abs(peak->mean() - energy);
        if( energy_diff < smallest_energy_diff )
        {
          smallest_energy_diff = energy_diff;
          nearest_fwhm = peak->fwhm();
        }
      }

      if( smallest_energy_diff < 0.2*energy )
        fwhm = nearest_fwhm;
    }//if( (fwhm <= 0.0) && peak_deque )
    
    // Finally, a wag based on detection type.
    const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
    const vector<float> pars = isHPGe ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
    fwhm = DetectorPeakResponse::peakResolutionFWHM( energy, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, pars );

    if( fwhm <= 0.0 )
      throw std::runtime_error("Could not determine FWHM for energy " + std::to_string(energy) );
    
    return get_nuclides_with_characteristics_in_energy_range( energy - fwhm, energy + fwhm, interspec );
  }//get_characteristics_near_energy(...)
  
  
  float get_expected_fwhm( const double energy, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");
    
    shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !meas || !spectrum )
      throw std::runtime_error("No foreground loaded");
    
    float fwhm = -1.0;
    
    // Try the currently loaded peak detector response
    if( meas->detector() && meas->detector()->hasResolutionInfo() )
      fwhm = meas->detector()->peakResolutionFWHM( static_cast<float>(energy) );
    
    // Estimate the FWHM response using detected peaks...
    shared_ptr<deque<shared_ptr<const PeakDef>>> peak_deque;
    if( fwhm <= 0.0 )
    {
      try
      {
        DetectedPeaksOptions opts;
        opts.specType = SpecUtils::SpectrumType::Foreground;
        const DetectedPeakStatus peaks = detected_peaks( opts, interspec );
        peak_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(peaks.peaks), end(peaks.peaks) );
        
        DetectorPeakResponse drf;
        drf.setIntrinsicEfficiencyFormula( "1.0", 2.54*PhysicalUnits::cm, PhysicalUnits::keV,
                                               0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarField );
        
        drf.fitResolution( peak_deque, spectrum, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
        
        fwhm = drf.peakResolutionFWHM( static_cast<float>(energy) );
      }catch( std::exception & )
      {
      }
    }//if( fwhm <= 0.0 )
    
    // Find nearest peak, and if within 20% of the energy, use it.
    if( (fwhm <= 0.0) && peak_deque )
    {
      double smallest_energy_diff = std::numeric_limits<double>::max();
      double nearest_fwhm = 0.0;
      for( const auto &peak : *peak_deque )
      {
        double energy_diff = std::abs(peak->mean() - energy);
        if( energy_diff < smallest_energy_diff )
        {
          smallest_energy_diff = energy_diff;
          nearest_fwhm = peak->fwhm();
        }
      }

      if( smallest_energy_diff < 0.2*energy )
        fwhm = nearest_fwhm;
    }//if( (fwhm <= 0.0) && peak_deque )
    
    // Finally, a wag based on detection type.
    if( fwhm <= 0.0 )
    {
      const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
      const vector<float> pars = isHPGe ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
      fwhm = DetectorPeakResponse::peakResolutionFWHM( energy, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, pars );
    }

    if( fwhm <= 0.0 )
      throw std::runtime_error("Could not determine FWHM for energy " + std::to_string(energy) );
    
    return fwhm;
  }//get_expected_fwhm(...)
  
  
  SpectrumCountsInEnergyRange get_counts_in_energy_range( double lower_energy, double upper_energy, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");

    shared_ptr<const SpecUtils::Measurement> foreground = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !foreground )
      throw runtime_error( "get_counts_in_energy_range: no foreground loaded." );

    if( lower_energy > upper_energy )
      std::swap( lower_energy, upper_energy );

    SpectrumCountsInEnergyRange answer;
    answer.lower_energy = lower_energy;
    answer.upper_energy = upper_energy;

    answer.foreground_counts = foreground->gamma_integral( static_cast<float>(lower_energy), static_cast<float>(upper_energy) );
    if( foreground->live_time() > 0.0 )
      answer.foreground_cps = answer.foreground_counts / foreground->live_time();
    else 
      answer.foreground_cps = std::numeric_limits<double>::quiet_NaN();
    
    shared_ptr<const SpecUtils::Measurement> background = interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
    if( background )
    {
      SpectrumCountsInEnergyRange::CountsWithComparisonToForeground back;

      back.counts = background->gamma_integral( static_cast<float>(lower_energy), static_cast<float>(upper_energy) );
      if( background->live_time() > 0.0 )
        back.cps = back.counts / background->live_time();
      else
        back.cps = std::numeric_limits<double>::quiet_NaN();

      const double backSF = foreground->live_time() / background->live_time();
      const double nfore = answer.foreground_counts;
      const double nback = back.counts;
      const double scaleback = nback * backSF;
      const double backsigma = sqrt(nback);
      const double forsigma = sqrt(nfore);
      const double backscalesigma = backSF * backsigma;
      const double total_back_fore_sigma = sqrt(backscalesigma*backscalesigma + forsigma*forsigma);
      const double nsigma = (nfore - scaleback) / total_back_fore_sigma;

      back.num_sigma_rel_foreground = nsigma;

      answer.background_info = std::move( back );
    }//if( background )

    shared_ptr<const SpecUtils::Measurement> secondary = interspec->displayedHistogram(SpecUtils::SpectrumType::SecondForeground);
    if( secondary )
    {
      SpectrumCountsInEnergyRange::CountsWithComparisonToForeground sec;

      sec.counts = secondary->gamma_integral( static_cast<float>(lower_energy), static_cast<float>(upper_energy) );
      if( secondary->live_time() > 0.0 )
        sec.cps = sec.counts / secondary->live_time();
      else
        sec.cps = std::numeric_limits<double>::quiet_NaN();

      const double secSF = foreground->live_time() / secondary->live_time();
      const double nfore = answer.foreground_counts;
      const double nsec = sec.counts;
      const double scalesec = nsec * secSF;
      const double secsigma = sqrt(nsec);
      const double forsigma = sqrt(nfore);
      const double secscalesigma = secSF * secsigma;
      const double total_sec_fore_sigma = sqrt(secscalesigma*secscalesigma + forsigma*forsigma);
      const double nsigma = (nfore - scalesec) / total_sec_fore_sigma;

      sec.num_sigma_rel_foreground = nsigma;

      answer.secondary_info = std::move( sec );
    }//if( secondary )

    return answer;
  }//SpectrumCountsInEnergyRange get_counts_in_energy_range( const double lower_energy, const double upper_energy, InterSpec *interspec )


  const char* to_string( EditPeakAction action )
  {
    switch( action )
    {
      case EditPeakAction::SetEnergy:                    return "SetEnergy";
      case EditPeakAction::SetFwhm:                      return "SetFwhm";
      case EditPeakAction::SetAmplitude:                 return "SetAmplitude";
      case EditPeakAction::SetEnergyUncertainty:         return "SetEnergyUncertainty";
      case EditPeakAction::SetFwhmUncertainty:           return "SetFwhmUncertainty";
      case EditPeakAction::SetAmplitudeUncertainty:      return "SetAmplitudeUncertainty";
      case EditPeakAction::SetRoiLower:                  return "SetRoiLower";
      case EditPeakAction::SetRoiUpper:                  return "SetRoiUpper";
      case EditPeakAction::SetSkewType:                  return "SetSkewType";
      case EditPeakAction::SetContinuumType:             return "SetContinuumType";
      case EditPeakAction::SetSource:                    return "SetSource";
      case EditPeakAction::SetColor:                     return "SetColor";
      case EditPeakAction::SetUserLabel:                 return "SetUserLabel";
      case EditPeakAction::SetUseForEnergyCalibration:   return "SetUseForEnergyCalibration";
      case EditPeakAction::SetUseForShieldingSourceFit:  return "SetUseForShieldingSourceFit";
      case EditPeakAction::SetUseForManualRelEff:        return "SetUseForManualRelEff";
      case EditPeakAction::DeletePeak:                   return "DeletePeak";
      case EditPeakAction::SplitFromRoi:                 return "SplitFromRoi";
      case EditPeakAction::MergeWithLeft:                return "MergeWithLeft";
      case EditPeakAction::MergeWithRight:               return "MergeWithRight";
    }//switch( action )

    throw runtime_error( "Invalid EditPeakAction value" );
  }//const char* to_string( EditPeakAction action )


  EditPeakAction edit_peak_action_from_string( const std::string &str )
  {
    if( str == "SetEnergy" )                   return EditPeakAction::SetEnergy;
    if( str == "SetFwhm" )                     return EditPeakAction::SetFwhm;
    if( str == "SetAmplitude" )                return EditPeakAction::SetAmplitude;
    if( str == "SetEnergyUncertainty" )        return EditPeakAction::SetEnergyUncertainty;
    if( str == "SetFwhmUncertainty" )          return EditPeakAction::SetFwhmUncertainty;
    if( str == "SetAmplitudeUncertainty" )     return EditPeakAction::SetAmplitudeUncertainty;
    if( str == "SetRoiLower" )                 return EditPeakAction::SetRoiLower;
    if( str == "SetRoiUpper" )                 return EditPeakAction::SetRoiUpper;
    if( str == "SetSkewType" )                 return EditPeakAction::SetSkewType;
    if( str == "SetContinuumType" )            return EditPeakAction::SetContinuumType;
    if( str == "SetSource" )                   return EditPeakAction::SetSource;
    if( str == "SetColor" )                    return EditPeakAction::SetColor;
    if( str == "SetUserLabel" )                return EditPeakAction::SetUserLabel;
    if( str == "SetUseForEnergyCalibration" )  return EditPeakAction::SetUseForEnergyCalibration;
    if( str == "SetUseForShieldingSourceFit" ) return EditPeakAction::SetUseForShieldingSourceFit;
    if( str == "SetUseForManualRelEff" )       return EditPeakAction::SetUseForManualRelEff;
    if( str == "DeletePeak" )                  return EditPeakAction::DeletePeak;
    if( str == "SplitFromRoi" )                return EditPeakAction::SplitFromRoi;
    if( str == "MergeWithLeft" )               return EditPeakAction::MergeWithLeft;
    if( str == "MergeWithRight" )              return EditPeakAction::MergeWithRight;

    throw runtime_error( "Invalid edit action: " + str );
  }//EditPeakAction edit_peak_action_from_string( const std::string &str )


  EditAnalysisPeakStatus edit_analysis_peak( const EditAnalysisPeakOptions &options, InterSpec *interspec )
  {
    if( !interspec )
      throw runtime_error( "edit_analysis_peak: No InterSpec session available" );

    // Get session ID
    string user_session = "direct_call";
    if( Wt::WApplication *app = Wt::WApplication::instance() )
      user_session = app->sessionId();

    // Override with provided session if available
    if( options.userSession )
      user_session = *options.userSession;

    // Get the PeakModel for the specified spectrum type
    PeakModel *peakModel = interspec->peakModel();
    if( !peakModel )
      throw runtime_error( "edit_analysis_peak: No peak model available" );

    // Find nearest peak
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( options.energy );
    if( !peak )
    {
      return EditAnalysisPeakStatus{
        user_session,
        false,
        "No peak found near " + std::to_string(options.energy) + " keV",
        std::nullopt,
        {}
      };
    }//if( !peak )

    // Validate peak is within 1 FWHM of requested energy
    const double fwhm = get_expected_fwhm( options.energy, interspec );
    const double energy_diff = std::fabs( peak->mean() - options.energy );
    if( energy_diff > fwhm )
    {
      return EditAnalysisPeakStatus{
        user_session,
        false,
        "Nearest peak at " + std::to_string(peak->mean()) + " keV is more than 1 FWHM ("
          + std::to_string(fwhm) + " keV) away from requested " + std::to_string(options.energy) + " keV",
        std::nullopt,
        {}
      };
    }//if( energy_diff > fwhm )

    // Handle delete action separately since it doesn't need a modified peak
    if( options.editAction == EditPeakAction::DeletePeak )
    {
      peakModel->removePeak( peak );
      return EditAnalysisPeakStatus{
        user_session,
        true,
        "Deleted peak at " + std::to_string(peak->mean()) + " keV",
        std::nullopt,
        {}
      };
    }//if( DeletePeak )

    // Create a modified copy of the peak
    PeakDef modifiedPeak = *peak;

    // Perform the edit based on action
    try
    {
      switch( options.editAction )
      {
        case EditPeakAction::SetEnergy:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetEnergy requires doubleValue parameter" );
          modifiedPeak.setMean( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetFwhm:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetFwhm requires doubleValue parameter" );
          // Convert FWHM to sigma (sigma = FWHM / 2.35482)
          modifiedPeak.setSigma( (*options.doubleValue) / 2.35482 );
          break;
        }

        case EditPeakAction::SetAmplitude:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetAmplitude requires doubleValue parameter" );
          modifiedPeak.setAmplitude( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetEnergyUncertainty:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetEnergyUncertainty requires doubleValue parameter" );
          modifiedPeak.setMeanUncert( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetFwhmUncertainty:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetFwhmUncertainty requires doubleValue parameter" );
          // Convert FWHM uncertainty to sigma uncertainty
          modifiedPeak.setSigmaUncert( (*options.doubleValue) / 2.35482 );
          break;
        }

        case EditPeakAction::SetAmplitudeUncertainty:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetAmplitudeUncertainty requires doubleValue parameter" );
          modifiedPeak.setAmplitudeUncert( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetRoiLower:
        case EditPeakAction::SetRoiUpper:
        {
          if( !options.doubleValue )
            throw runtime_error( string(to_string(options.editAction)) + " requires doubleValue parameter" );

          shared_ptr<PeakContinuum> continuum = modifiedPeak.getContinuum();
          const double current_lower = continuum->lowerEnergy();
          const double current_upper = continuum->upperEnergy();

          const double new_lower = (options.editAction == EditPeakAction::SetRoiLower)
                                    ? *options.doubleValue : current_lower;
          const double new_upper = (options.editAction == EditPeakAction::SetRoiUpper)
                                    ? *options.doubleValue : current_upper;

          if( new_lower >= new_upper )
            throw runtime_error( "ROI lower energy must be less than upper energy" );

          continuum->setRange( new_lower, new_upper );
          break;
        }

        case EditPeakAction::SetSkewType:
        {
          if( !options.stringValue )
            throw runtime_error( "SetSkewType requires stringValue parameter" );

          // Parse skew type string
          PeakDef::SkewType skew_type;
          const string &skew_str = *options.stringValue;
          if( skew_str == "NoSkew" )                      skew_type = PeakDef::SkewType::NoSkew;
          else if( skew_str == "Bortel" )                 skew_type = PeakDef::SkewType::Bortel;
          else if( skew_str == "GaussExp" )               skew_type = PeakDef::SkewType::GaussExp;
          else if( skew_str == "CrystalBall" )            skew_type = PeakDef::SkewType::CrystalBall;
          else if( skew_str == "ExpGaussExp" )            skew_type = PeakDef::SkewType::ExpGaussExp;
          else if( skew_str == "DoubleSidedCrystalBall" ) skew_type = PeakDef::SkewType::DoubleSidedCrystalBall;
          else throw runtime_error( "Invalid skew type: " + skew_str );

          modifiedPeak.setSkewType( skew_type );
          break;
        }

        case EditPeakAction::SetContinuumType:
        {
          if( !options.stringValue )
            throw runtime_error( "SetContinuumType requires stringValue parameter" );

          // Parse continuum type string
          PeakContinuum::OffsetType continuum_type;
          const string &cont_str = *options.stringValue;
          if( cont_str == "None" )           continuum_type = PeakContinuum::OffsetType::NoOffset;
          else if( cont_str == "Constant" )  continuum_type = PeakContinuum::OffsetType::Constant;
          else if( cont_str == "Linear" )    continuum_type = PeakContinuum::OffsetType::Linear;
          else if( cont_str == "Quadratic" ) continuum_type = PeakContinuum::OffsetType::Quadratic;
          else if( cont_str == "Cubic" )     continuum_type = PeakContinuum::OffsetType::Cubic;
          else if( cont_str == "FlatStep" )  continuum_type = PeakContinuum::OffsetType::FlatStep;
          else if( cont_str == "LinearStep" )continuum_type = PeakContinuum::OffsetType::LinearStep;
          else if( cont_str == "BiLinearStep" )continuum_type = PeakContinuum::OffsetType::BiLinearStep;
          else if( cont_str == "External" )  continuum_type = PeakContinuum::OffsetType::External;
          else throw runtime_error( "Invalid continuum type: " + cont_str );

          shared_ptr<PeakContinuum> continuum = modifiedPeak.getContinuum();
          continuum->setType( continuum_type );
          break;
        }

        case EditPeakAction::SetSource:
        {
          if( !options.stringValue )
            throw runtime_error( "SetSource requires stringValue parameter" );

          // Use PeakModel::setNuclideXrayReaction to parse and set the source
          PeakModel::SetGammaSource result = PeakModel::setNuclideXrayReaction( modifiedPeak, *options.stringValue, 4.0 );
          if( result == PeakModel::SetGammaSource::NoSourceChange )
            throw runtime_error( "Failed to set source: " + *options.stringValue );
          break;
        }

        case EditPeakAction::SetColor:
        {
          if( !options.stringValue )
            throw runtime_error( "SetColor requires stringValue parameter" );

          // Parse CSS color string to Wt::WColor
          Wt::WColor color( *options.stringValue );
          modifiedPeak.setLineColor( color );
          break;
        }

        case EditPeakAction::SetUserLabel:
        {
          if( !options.stringValue )
            throw runtime_error( "SetUserLabel requires stringValue parameter" );

          modifiedPeak.setUserLabel( *options.stringValue );
          break;
        }

        case EditPeakAction::SetUseForEnergyCalibration:
        {
          if( !options.boolValue )
            throw runtime_error( "SetUseForEnergyCalibration requires boolValue parameter" );

          modifiedPeak.useForEnergyCalibration( *options.boolValue );
          break;
        }

        case EditPeakAction::SetUseForShieldingSourceFit:
        {
          if( !options.boolValue )
            throw runtime_error( "SetUseForShieldingSourceFit requires boolValue parameter" );

          modifiedPeak.useForShieldingSourceFit( *options.boolValue );
          break;
        }

        case EditPeakAction::SetUseForManualRelEff:
        {
          if( !options.boolValue )
            throw runtime_error( "SetUseForManualRelEff requires boolValue parameter" );

          modifiedPeak.useForManualRelEff( *options.boolValue );
          break;
        }

        case EditPeakAction::SplitFromRoi:
        {
          // Make the peak have its own unique continuum
          modifiedPeak.makeUniqueNewContinuum();
          break;
        }

        case EditPeakAction::MergeWithLeft:
        case EditPeakAction::MergeWithRight:
        {
          // Find the adjacent peak to merge with
          shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
          if( !all_peaks_ptr )
            throw runtime_error( "No peaks available in model" );

          const deque<shared_ptr<const PeakDef>> &all_peaks = *all_peaks_ptr;

          // Find current peak index
          auto peak_iter = find( all_peaks.begin(), all_peaks.end(), peak );
          if( peak_iter == all_peaks.end() )
            throw runtime_error( "Could not find peak in model" );

          shared_ptr<const PeakDef> adjacent_peak;
          if( options.editAction == EditPeakAction::MergeWithLeft )
          {
            if( peak_iter == all_peaks.begin() )
              throw runtime_error( "No peak to the left to merge with" );
            adjacent_peak = *(peak_iter - 1);
          }
          else // MergeWithRight
          {
            if( (peak_iter + 1) == all_peaks.end() )
              throw runtime_error( "No peak to the right to merge with" );
            adjacent_peak = *(peak_iter + 1);
          }

          // Share the continuum with the adjacent peak
          // Need to const_cast because setContinuum expects non-const shared_ptr
          modifiedPeak.setContinuum( const_pointer_cast<PeakContinuum>(adjacent_peak->continuum()) );
          break;
        }

        case EditPeakAction::DeletePeak:
          // Already handled above
          assert( false );
          break;
      }//switch( options.editAction )

      // Update the peak in the model
      peakModel->updatePeak( peak, modifiedPeak );

      // Get all peaks in the same ROI for the return value
      vector<shared_ptr<const PeakDef>> peaks_in_roi;
      shared_ptr<const PeakDef> updated_peak;
      shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
      if( all_peaks_ptr )
      {
        const deque<shared_ptr<const PeakDef>> &all_peaks = *all_peaks_ptr;
        shared_ptr<const PeakContinuum> continuum = modifiedPeak.continuum();

        for( const auto &p : all_peaks )
        {
          if( p->continuum() == continuum )
          {
            peaks_in_roi.push_back( p );
            // Find the updated peak by comparing energies
            if( !updated_peak && fabs(p->mean() - modifiedPeak.mean()) < 0.01 )
              updated_peak = p;
          }
        }
      }

      // If we couldn't find the updated peak in the model, use our local copy
      if( !updated_peak )
        updated_peak = make_shared<const PeakDef>(modifiedPeak);

      return EditAnalysisPeakStatus{
        user_session,
        true,
        "Successfully performed " + string(to_string(options.editAction)) + " on peak at " + std::to_string(modifiedPeak.mean()) + " keV",
        updated_peak,
        peaks_in_roi
      };
    }
    catch( const exception &e )
    {
      return EditAnalysisPeakStatus{
        user_session,
        false,
        string("Error performing ") + to_string(options.editAction) + ": " + e.what(),
        std::nullopt,
        {}
      };
    }//try / catch
  }//EditAnalysisPeakStatus edit_analysis_peak( const EditAnalysisPeakOptions &options, InterSpec *interspec )

} // namespace AnalystChecks
