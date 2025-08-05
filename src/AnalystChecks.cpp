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
    
    if( options.source.has_value() )
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
    
    
    if( options.addToUsersPeaks )
    {
      if( options.specType == SpecUtils::SpectrumType::Foreground )
      {
        PeakModel *pmodel = interspec->peakModel();
        assert( pmodel );
        if( pmodel )
        {
          pmodel->removePeaks( peaks_to_remove );
          pmodel->addPeaks( peaks_to_add_in );
        }//if( pmodel )
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
        interspec->setPeaks( options.specType, new_deque );
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
      status.peaks.insert( end(status.peaks), begin(*meas_peaks), end(*meas_peaks) );
    
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
  
  
  
} // namespace AnalystChecks
