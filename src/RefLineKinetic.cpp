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

#include <regex>

#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RefLineKinetic.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/ExternalRidResult.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/ReferenceLinePredef.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/IsotopeSearchByEnergy.h"

using namespace std;
using namespace Wt;


struct AlwaysSrcs
{
  const map<string,ReferenceLinePredef::NucMix> nuc_mixes;
  const map<string,ReferenceLinePredef::CustomSrcLines> custom_lines;
  const std::vector<ReferenceLinePredef::IndividualSource> individual_sources;
  
  AlwaysSrcs() = delete;
  AlwaysSrcs( map<string,ReferenceLinePredef::NucMix> &&mixes,
             map<string,ReferenceLinePredef::CustomSrcLines> &&lines,
             std::vector<ReferenceLinePredef::IndividualSource> &&srcs )
  : nuc_mixes( std::move(mixes) ),
  custom_lines( std::move(lines) ),
  individual_sources( std::move(srcs) )
  {
    
  }
};//struct AlwaysSrcs



RefLineKinetic::RefLineKinetic( D3SpectrumDisplayDiv *chart, InterSpec *interspec )
  : Wt::WObject( chart ),
  m_interspec( interspec ),
  m_chart( chart ),
  m_active( false ),
  m_has_inited( false ),
  m_init_error_msg{},
  m_always_srcs{},
  m_external_rid_results( nullptr ),
  m_renderFlags{}
{
  if( !m_interspec || !m_chart )
    throw std::runtime_error( "RefLineKinetic: null InterSpec parent or chart" );
  
  m_active = UserPreferences::preferenceValue<bool>( "KineticRefLine", m_interspec );
  m_interspec->preferences()->addCallbackWhenChanged( "KineticRefLine", this, &RefLineKinetic::setActive );
    
  m_interspec->hintPeaksSet().connect( boost::bind(&RefLineKinetic::autoSearchPeaksSet, this, boost::placeholders::_1) );
  m_interspec->displayedSpectrumChanged().connect( boost::bind(&RefLineKinetic::spectrumChanged, this,
    boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4
  ) );
  
  m_interspec->externalRidResultsRecieved().connect( boost::bind( &RefLineKinetic::autoRidResultsRecieved, this, boost::placeholders::_1 ) );
  
  m_chart->setKineticRefLineController( this );
  
  start_init_always_sources();
}//RefLineKinetic constructor


RefLineKinetic::~RefLineKinetic()
{
}


void RefLineKinetic::start_init_always_sources()
{
  try
  {
    // TODO: move this to being done in a background thread
    m_has_inited = true;
    
    string always_defs_file = SpecUtils::append_path(InterSpec::writableDataDirectory(), "kinetic_ref_lines.xml");
    if( !SpecUtils::is_file(always_defs_file) )
      always_defs_file = SpecUtils::append_path(InterSpec::staticDataDirectory(), "kinetic_ref_lines.xml");
    
    map<string,ReferenceLinePredef::NucMix> nuc_mixes;
    map<string,ReferenceLinePredef::CustomSrcLines> custom_lines;
    std::vector<ReferenceLinePredef::IndividualSource> indiv_sources;
    ReferenceLinePredef::load_ref_line_file( always_defs_file, nuc_mixes, custom_lines, &indiv_sources );
    
    m_always_srcs = make_unique<AlwaysSrcs>( std::move(nuc_mixes), std::move(custom_lines), std::move(indiv_sources) );
  }catch( std::exception &e )
  {
    cerr << "Failed to initialize RefLineKinetic: " << e.what() << endl;
    m_init_error_msg = e.what();
  }//try / catch
}//void start_init_always_sources()


void RefLineKinetic::setActive( bool active )
{
  if( m_active == active )
    return;
  
  m_active = active;
  
  m_renderFlags |= KineticRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderKineticRefLine();
}//void RefLineKinetic::setActive( bool active )


bool RefLineKinetic::isActive() const
{
  return m_active;
}


void RefLineKinetic::autoSearchPeaksSet( const SpecUtils::SpectrumType spectrum )
{
  m_renderFlags |= KineticRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderKineticRefLine();
}//void autoSearchPeaksSet(...)


void RefLineKinetic::spectrumChanged( const SpecUtils::SpectrumType spec_type,
                     const std::shared_ptr<SpecMeas> &measurement,
                     const std::set<int> &sample_numbers,
                     const std::vector<std::string> &detectors )
{
  if( spec_type == SpecUtils::SpectrumType::Foreground )
    m_external_rid_results.reset();
  
  m_renderFlags |= KineticRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderKineticRefLine();
}//void spectrumChanged(...)


void RefLineKinetic::autoRidResultsRecieved( const std::shared_ptr<const ExternalRidResults> &results )
{
  m_external_rid_results = results;
  
  m_renderFlags |= KineticRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderKineticRefLine();
}//void autoRidResultsRecieved( const std::shared_ptr<const ExternalRidResults> &results )


void RefLineKinetic::pushUpdates()
{
  if( m_renderFlags.testFlag(KineticRefLineRenderFlags::UpdateLines) )
    updateLines();
  
  m_renderFlags = 0;
}//void pushUpdates()


void RefLineKinetic::filterLines( ReferenceLineInfo &ref_lines,
                                 const RefLineKinetic::RefLineSrc src,
                                 const shared_ptr<const SpecUtils::Measurement> &meas )
{
  if( !meas )
    return;
  
  // TODO: implement
  // importance = yield(i)*sqrt(energy(i))/sum(yield*sqrt(energy))
}//void filterLines(...)


void RefLineKinetic::updateLines()
{
  if( !m_active )
  {
    m_chart->setKineticRefernceLines( {}, "" );
    return;
  }//if( !m_active )
  
  shared_ptr<const SpecUtils::Measurement> foreground = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  const set<int> &foreground_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
  const shared_ptr<const SpecMeas> foreground_meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  
  shared_ptr<const SpecUtils::Measurement> background = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
  const set<int> &background_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
  const shared_ptr<const SpecMeas> background_meas = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  
  RefLineInput base_input;
  base_input.m_lower_br_cutt_off = 0.0;
  base_input.m_promptLinesOnly = false;
  base_input.m_showGammas = true;
  base_input.m_showXrays = true;
  base_input.m_showAlphas = false;
  base_input.m_showBetas = false;
  base_input.m_showCascades = false;
  base_input.m_showEscapes = false;
  
  const shared_ptr<const DetectorPeakResponse> detector = foreground_meas ? foreground_meas->detector() : nullptr;
  if( detector && detector->isValid() )
  {
    base_input.m_detector_name = detector->name();
    base_input.m_det_intrinsic_eff = detector->intrinsicEfficiencyFcn();
  }
  
  vector<pair<double,ReferenceLineInfo>> ref_lines;
  
  if( m_always_srcs )
  {
    // Convert individual sources to ReferenceLineInfo objects
    for( const auto &src : m_always_srcs->individual_sources )
    {
      RefLineInput input = base_input;
      input.m_input_txt = src.m_name;
      
      if( src.m_age.has_value() )
        input.m_age = PhysicalUnits::printToBestTimeUnits( src.m_age.value() * PhysicalUnits::year );
      
      if( src.m_color.has_value() )
        input.m_color = src.m_color.value();
      
      if( src.is_background )
      {
        // TODO: transport through a soil sphere, using it as a trace source
      }
      
      if( src.shielding_material.has_value() )
      {
        const std::string &shielding_material = *src.shielding_material;
        
        if( src.shielding_thickness.has_value() )
        {
          input.m_shielding_name = *src.shielding_material;
          input.m_shielding_thickness = PhysicalUnits::printToBestLengthUnits( *src.shielding_thickness, 2 );
        }else
        {
          // Check if src.shielding_material is a string that contains AN=<number> and AD=<number> then is a generic shielding
          std::regex an_pattern( R"(an\s*=\s*([0-9]+(?:\.[0-9]+)?))", std::regex_constants::icase );
          std::regex ad_pattern( R"(ad\s*=\s*([0-9]+(?:\.[0-9]+)?))", std::regex_constants::icase );
          
          std::smatch an_match, ad_match;
          if( std::regex_search(shielding_material, an_match, an_pattern)
             && std::regex_search(shielding_material, ad_match, ad_pattern) )
          {
            input.m_shielding_an = an_match[1].str();
            input.m_shielding_ad = ad_match[1].str();
          }
        }//if( src.shielding_thickness.has_value() ) / else
        
        if( (!input.m_shielding_name.empty() && input.m_shielding_thickness.empty())
           || (!input.m_shielding_an.empty() && !input.m_shielding_ad.empty()) )
        {
          try
          {
            input.setShieldingAttFcn( m_interspec->materialDataBase() );
          }catch( std::exception &e )
          {
            cerr << "Error setting shielding for curve: " << e.what() << endl;
          }
        }//if( we have a defined shielding )
      }//if( src.shielding_material.has_value() )
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        filterLines( *ref_info, RefLineSrc::AlwaysShowing, foreground );
        ref_lines.emplace_back( src.m_weight, std::move(*ref_info) );
      }
      else
      {
        cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for individual source: " << src.m_name << endl;
      }
    }//for( const auto &src : m_always_srcs->individual_sources )
    
    // Convert nuclide mixtures to ReferenceLineInfo objects
    for( const auto &mix_pair : m_always_srcs->nuc_mixes )
    {
      const ReferenceLinePredef::NucMix &mix = mix_pair.second;
      
      RefLineInput input = base_input;
      input.m_input_txt = mix.m_name;
      input.m_age = mix.m_default_age_str;
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        filterLines( *ref_info, RefLineSrc::AlwaysShowing, foreground );
        ref_lines.emplace_back( mix.m_weight, std::move(*ref_info) );
      }
      else
      {
        cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for nuclide mixture: " << mix.m_name << endl;
      }
    }//for( const auto &mix_pair : m_always_srcs->nuc_mixes )
    
    // Convert custom source lines to ReferenceLineInfo objects
    for( const auto &custom_pair : m_always_srcs->custom_lines )
    {
      const ReferenceLinePredef::CustomSrcLines &custom = custom_pair.second;
      
      RefLineInput input = base_input;
      input.m_input_txt = custom.m_name;
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        filterLines( *ref_info, RefLineSrc::AlwaysShowing, foreground );
        ref_lines.emplace_back( custom.m_weight, std::move(*ref_info) );
      }
      else
      {
        cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for custom source: " << custom.m_name << endl;
      }
    }//for( const auto &custom_pair : m_always_srcs->custom_lines )
  }//if( m_always_srcs )
  
  // Now check which nuclides peaks have been assigned to, and add those lines, if they havent already been added
  PeakModel * const pmodel = m_interspec->peakModel();
  const shared_ptr<const deque<std::shared_ptr<const PeakDef>>> user_foreground_peaks = pmodel ? pmodel->peaks() : nullptr;
  if( user_foreground_peaks && user_foreground_peaks->size() )
  {
    // We will assign these nuclides a weight of 10
    for( const auto& peak : *user_foreground_peaks )
    {
      if( !peak )
        continue;
        
      std::string canonical_name;
      bool found_source = false;
      
      // Try to get the canonical name using PeakDef source methods
      if( const SandiaDecay::Nuclide *nuc = peak->parentNuclide() )
      {
        canonical_name = nuc->symbol;
        found_source = true;
      }
      else if( const SandiaDecay::Element *el = peak->xrayElement() )
      {
        canonical_name = el->symbol;
        found_source = true;
      }
      else if( const ReactionGamma::Reaction *rctn = peak->reaction() )
      {
        canonical_name = rctn->name();
        found_source = true;
      }
      
      if( !found_source )
        continue;
        
      // Check if this source is already in ref_lines using canonical name, and update weight if higher
      bool already_exists = false;
      double new_weight = 10.0;
      for( auto& existing : ref_lines )
      {
        if( existing.second.m_input.m_input_txt == canonical_name )
        {
          already_exists = true;
          if( new_weight > existing.first )
            existing.first = new_weight;
          break;
        }
      }
      
      if( already_exists )
        continue;
      
      // Create reference lines using canonical name
      RefLineInput input = base_input;
      input.m_input_txt = canonical_name;
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        filterLines( *ref_info, RefLineSrc::CharacteristicLine, foreground );
        ref_lines.emplace_back( new_weight, std::move(*ref_info) );
      }
      else
      {
        cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for user peak source: " << canonical_name << endl;
      }
    }
  }//if( user_foreground_peaks && user_foreground_peaks->size() )
  
  
  // Now check for external RID results, and add those lines that havent been added
  if( m_external_rid_results )
  {
    // We will assign these nuclides a weight of 5
    const auto& isotopes = m_external_rid_results->isotopes;
    
    for( const auto& isotope : isotopes )
    {
      if( isotope.is_null() )
        continue;
        
      std::string canonical_name;
      bool found_source = false;
      
      // Try to get the canonical name using the helper methods
      if( const SandiaDecay::Nuclide *nuc = isotope.nuclide() )
      {
        canonical_name = nuc->symbol;
        found_source = true;
      }
      else if( const SandiaDecay::Element *el = isotope.element() )
      {
        canonical_name = el->symbol;
        found_source = true;
      }
      else if( const ReactionGamma::Reaction *rctn = isotope.reaction() )
      {
        canonical_name = rctn->name();
        found_source = true;
      }
      
      if( !found_source )
        continue;
        
      // Check if this source is already in ref_lines using canonical name, and update weight if higher
      bool already_exists = false;
      double new_weight = 5.0;
      for( auto& existing : ref_lines )
      {
        if( existing.second.m_input.m_input_txt == canonical_name )
        {
          already_exists = true;
          if( new_weight > existing.first )
            existing.first = new_weight;
          break;
        }
      }
      
      if( already_exists )
        continue;
      
      // Create reference lines using canonical name
      RefLineInput input = base_input;
      input.m_input_txt = canonical_name;
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        filterLines( *ref_info, RefLineSrc::ExternalRid, foreground );
        ref_lines.emplace_back( new_weight, std::move(*ref_info) );
      }
      else
      {
        cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for external RID source: " << canonical_name << endl;
      }
    }
  }//if( m_external_rid_results )
  
  
  // Use the detectors on-board RID results
  if( foreground_meas && foreground_meas->detectors_analysis() )
  {
    // We will assign these nuclides a weight of 3
    const auto detector_analysis = foreground_meas->detectors_analysis();
    const auto& results = detector_analysis->results_;
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const ReactionGamma *reaction_db = nullptr;
    try
    {
      reaction_db = ReactionGammaServer::database();
    }catch(...)
    {
      // ReactionGamma database not available, skip reaction parsing
    }
    
    for( const auto& result : results )
    {
      if( result.nuclide_.empty() )
        continue;
        
      std::string canonical_name;
      bool found_source = false;
      
      // Try to get the nuclide from SandiaDecay database
      const SandiaDecay::Nuclide *nuc = db ? db->nuclide( result.nuclide_ ) : nullptr;
      if( nuc )
      {
        canonical_name = nuc->symbol;
        found_source = true;
      }
      else
      {
        // Try to get element from SandiaDecay database
        const SandiaDecay::Element *el = db ? db->element( result.nuclide_ ) : nullptr;
        if( el )
        {
          canonical_name = el->symbol;
          found_source = true;
        }
        else if( reaction_db )
        {
          // Try to get reaction from ReactionGamma database
          vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
          reaction_db->gammas( result.nuclide_, possible_rctns );
          
          if( !possible_rctns.empty() )
          {
            canonical_name = possible_rctns[0].reaction->name();
            found_source = true;
          }
        }
      }
      
      if( !found_source )
      {
        // TODO: check for algorithm specific names, like HEU, SNM, neutrons, etc
        continue;
      }
        
      // Check if this source is already in ref_lines using canonical name, and update weight if higher
      bool already_exists = false;
      double new_weight = 3.0;
      for( auto& existing : ref_lines )
      {
        if( existing.second.m_input.m_input_txt == canonical_name )
        {
          already_exists = true;
          if( new_weight > existing.first )
            existing.first = new_weight;
          break;
        }
      }
      
      if( already_exists )
        continue;
      
      // Create reference lines using canonical name
      RefLineInput input = base_input;
      input.m_input_txt = canonical_name;
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        filterLines( *ref_info, RefLineSrc::OnboardRid, foreground );
        ref_lines.emplace_back( new_weight, std::move(*ref_info) );
      }
      else
      {
        cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for on-board RID source: " << canonical_name << endl;
      }
    }
  }//if( foreground_meas && foreground_meas->detectors_analysis() )
  
  
  // If we want to be heavy-handed, we could use (some of the) suggestions from populateCandidateNuclides
  shared_ptr<const deque<shared_ptr<const PeakDef>>> foreground_autosearch_peaks = foreground_meas
            ? foreground_meas->automatedSearchPeaks(foreground_samples)
            : nullptr;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> background_autosearch_peaks = background_meas
            ? background_meas->automatedSearchPeaks(background_samples)
            : nullptr;
  const shared_ptr<const deque<std::shared_ptr<const PeakDef>>> user_background_peaks = 
            (background_meas && background_meas->sampleNumsWithPeaks().count(background_samples)) 
            ? background_meas->peaks(background_samples) 
            : nullptr;
  
  // We'll treat user peaks that dont have a source associated with them, as automated serach peaks.
  if( (!foreground_autosearch_peaks || foreground_autosearch_peaks->empty()) 
      && user_foreground_peaks && !user_foreground_peaks->empty() )
  {
    auto non_idd_peaks = make_shared<deque<shared_ptr<const PeakDef>>>();
    for( const auto &p : *user_foreground_peaks )
    {
      if( p && !p->parentNuclide() && !p->xrayElement() && !p->reaction() )
        non_idd_peaks->push_back( p );
    }

    if( !non_idd_peaks->empty() )
      foreground_autosearch_peaks = non_idd_peaks;
  }//

  // Process automated search peaks to find characteristic nuclides
  if( foreground_autosearch_peaks && foreground_autosearch_peaks->size() )
  {
    // Helper lambda to check if a peak has a source assigned (nuclide, xray, or reaction)
    auto peak_has_source = []( const std::shared_ptr<const PeakDef>& p ) -> bool {
      return p && ( p->parentNuclide() != nullptr || p->xrayElement() != nullptr || p->reaction() != nullptr );
    };
    
    // Helper lambda to check if two peaks overlap within 1 FWHM
    auto peaks_overlap = []( const std::shared_ptr<const PeakDef>& p1, const std::shared_ptr<const PeakDef>& p2 ) -> bool {
      if( !p1 || !p2 )
        return false;
      const double fwhm1 = p1->fwhm();
      const double fwhm2 = p2->fwhm();
      const double avg_fwhm = 0.5 * (fwhm1 + fwhm2);
      return fabs(p1->mean() - p2->mean()) <= avg_fwhm;
    };
    
    // Create a list of remaining automated search peaks after filtering
    vector<std::shared_ptr<const PeakDef>> remaining_autosearch_peaks;
    
    for( const auto& auto_peak : *foreground_autosearch_peaks )
    {
      if( !auto_peak )
        continue;
        
      bool should_skip = false;
      
      // Skip if this autosearch peak overlaps with a user foreground peak that has a source assigned
      if( user_foreground_peaks )
      {
        for( const auto& user_peak : *user_foreground_peaks )
        {
          if( peaks_overlap(auto_peak, user_peak) )
          {
            should_skip = true;
            break;
          }
        }
      }
      
      // Skip if this autosearch peak overlaps with any background peaks
      if( !should_skip && user_background_peaks )
      {
        for( const auto& bg_peak : *user_background_peaks )
        {
          if( peaks_overlap(auto_peak, bg_peak) )
          {
            should_skip = true;
            break;
          }
        }
      }
      
      if( !should_skip && background_autosearch_peaks )
      {
        for( const auto& bg_auto_peak : *background_autosearch_peaks )
        {
          if( peaks_overlap(auto_peak, bg_auto_peak) )
          {
            should_skip = true;
            break;
          }
        }
      }
      
      if( !should_skip )
        remaining_autosearch_peaks.push_back( auto_peak );
    }
    
    // Also process user foreground peaks that have no source assigned
    vector<std::shared_ptr<const PeakDef>> unsourced_user_peaks;
    if( user_foreground_peaks )
    {
      for( const auto& user_peak : *user_foreground_peaks )
      {
        if( user_peak && !peak_has_source(user_peak) )
          unsourced_user_peaks.push_back( user_peak );
      }
    }
    
    // Helper function to extract nuclide/element/reaction name from characteristic string
    auto extract_source_name = []( const std::string& characteristic ) -> std::string {
      // Format is typically: "Co60 1173.228 keV", "Fe xray 6.4 keV", "Fe(n,g) 7631.1 keV"
      size_t space_pos = characteristic.find(' ');
      if( space_pos != std::string::npos )
      {
        std::string source_name = characteristic.substr(0, space_pos);
        // Handle xray case: "Fe xray" -> "Fe"
        if( source_name.find("xray") != std::string::npos )
        {
          size_t xray_pos = characteristic.find(" xray");
          if( xray_pos != std::string::npos )
            return characteristic.substr(0, xray_pos);
        }
        return source_name;
      }
      return "";
    };
    
    // First, collect all source names from the top 5 characteristics of each peak
    std::set<std::string> candidate_source_names;
    
    // Process remaining automated search peaks
    if( !remaining_autosearch_peaks.empty() )
    {
      cout << "Finding characteristics for " << remaining_autosearch_peaks.size() << " remaining automated search peaks:" << endl;
      for( const auto& peak : remaining_autosearch_peaks )
      {
        vector<string> characteristicnucs;
        IsotopeId::findCharacteristics( characteristicnucs, peak );
        
        cout << "Peak at " << peak->mean() << " keV: ";
        if( characteristicnucs.empty() )
        {
          cout << "No characteristic nuclides found" << endl;
        }
        else
        {
          cout << "Found " << characteristicnucs.size() << " candidates: ";
          
          // Collect top 5 characteristics
          for( size_t i = 0; i < characteristicnucs.size() && i < 5; ++i )
          {
            if( i > 0 ) cout << ", ";
            cout << characteristicnucs[i];
            
            std::string source_name = extract_source_name( characteristicnucs[i] );
            if( !source_name.empty() )
              candidate_source_names.insert( source_name );
          }
          
          if( characteristicnucs.size() > 5 )
            cout << " (and " << (characteristicnucs.size() - 5) << " more)";
          cout << endl;
        }
      }
    }
    
    // Process unsourced user peaks
    if( !unsourced_user_peaks.empty() )
    {
      cout << "Finding characteristics for " << unsourced_user_peaks.size() << " user peaks without sources:" << endl;
      for( const auto& peak : unsourced_user_peaks )
      {
        vector<string> characteristicnucs;
        IsotopeId::findCharacteristics( characteristicnucs, peak );
        
        cout << "User peak at " << peak->mean() << " keV: ";
        if( characteristicnucs.empty() )
        {
          cout << "No characteristic nuclides found" << endl;
        }
        else
        {
          cout << "Found " << characteristicnucs.size() << " candidates: ";
          
          // Collect top 5 characteristics
          for( size_t i = 0; i < characteristicnucs.size() && i < 5; ++i )
          {
            if( i > 0 ) cout << ", ";
            cout << characteristicnucs[i];
            
            std::string source_name = extract_source_name( characteristicnucs[i] );
            if( !source_name.empty() )
              candidate_source_names.insert( source_name );
          }
          
          if( characteristicnucs.size() > 5 )
            cout << " (and " << (characteristicnucs.size() - 5) << " more)";
          cout << endl;
        }
      }
    }
    
    // Now filter the collected source names once
    std::set<std::string> characteristic_sources;
    
    if( !candidate_source_names.empty() )
    {
      IsotopeSearchByEnergy * const search = m_interspec->nuclideSearch();
      const std::vector<IsotopeSearchByEnergy::NucSearchCategory> categories = search ? search->search_categories() : std::vector<IsotopeSearchByEnergy::NucSearchCategory>{};
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      const ReactionGamma *reaction_db = nullptr;
      try
      {
        reaction_db = ReactionGammaServer::database();
      }catch(...){}
      
      for( const auto& source_name : candidate_source_names )
      {
        // Check if already in ref_lines
        bool already_in_ref_lines = false;
        for( const auto& existing : ref_lines )
        {
          if( existing.second.m_input.m_input_txt == source_name )
          {
            already_in_ref_lines = true;
            break;
          }
        }
        
        if( already_in_ref_lines || !search || categories.empty() || !db )
          continue;
        
        bool is_in_any_category = false;
        
        // Check if it's a nuclide in any category
        const SandiaDecay::Nuclide *nuc = db->nuclide( source_name );
        if( nuc )
        {
          // Filter out any nuclide that is a descendant of any nuclide already in ref_lines
          bool is_descendant = false;
          for( const auto& existing : ref_lines )
          {
            const SandiaDecay::Nuclide *existing_nuc = db->nuclide( existing.second.m_input.m_input_txt );
            if( existing_nuc && nuc->branchRatioFromForebear( existing_nuc ) > 0.0 )
            {
              is_descendant = true;
              break;
            }
          }
          
          if( !is_descendant )
          {
            is_in_any_category = (
              IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_medical_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_industrial_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_norm_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_snm_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_common_category_key, categories)
              // || IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_fission_category_key, categories)
            );
          }
        }
        else
        {
          // Check if it's an element
          const SandiaDecay::Element *el = db->element( source_name );
          if( el )
          {
            is_in_any_category = (
              IsotopeSearchByEnergy::is_in_category(el, IsotopeSearchByEnergy::sm_medical_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(el, IsotopeSearchByEnergy::sm_industrial_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(el, IsotopeSearchByEnergy::sm_norm_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(el, IsotopeSearchByEnergy::sm_snm_category_key, categories) ||
              IsotopeSearchByEnergy::is_in_category(el, IsotopeSearchByEnergy::sm_common_category_key, categories)
              // || IsotopeSearchByEnergy::is_in_category(el, IsotopeSearchByEnergy::sm_fission_category_key, categories )
            );
          }
          else if( reaction_db )
          {
            // Check if it's a reaction
            vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
            reaction_db->gammas( source_name, possible_rctns );
            
            if( !possible_rctns.empty() )
            {
              const ReactionGamma::Reaction *rctn = possible_rctns[0].reaction;
              
              // Check if we have enough peaks matching gamma lines for this reaction
              const auto& gammas = rctn->gammas;
              int required_matches = (gammas.size() > 1) ? 2 : 1;  // Single-gamma reactions need only 1 match
              int peak_matches = 0;
              
              // Collect all peaks (user foreground peaks + autosearch peaks)
              vector<std::shared_ptr<const PeakDef>> all_peaks;
              if( user_foreground_peaks )
              {
                for( const auto& p : *user_foreground_peaks )
                  if( p ) all_peaks.push_back( p );
              }
              if( foreground_autosearch_peaks )
              {
                for( const auto& p : *foreground_autosearch_peaks )
                  if( p ) all_peaks.push_back( p );
              }
              
              // Check how many peaks match gamma lines from this reaction
              for( const auto& gamma : gammas )
              {
                for( const auto& peak : all_peaks )
                {
                  const double peak_energy = peak->mean();
                  const double fwhm = peak->fwhm();
                  const double tolerance = fwhm > 0 ? fwhm : 5.0;  // Use FWHM or 5 keV default
                  
                  if( fabs(peak_energy - gamma.energy) <= tolerance )
                  {
                    peak_matches++;
                    break;  // Found a match for this gamma, move to next gamma
                  }
                }
                
                if( peak_matches >= required_matches )
                  break;  // We have enough matches
              }
              
              if( peak_matches >= required_matches )
              {
                is_in_any_category = (
                  IsotopeSearchByEnergy::is_in_category(rctn, IsotopeSearchByEnergy::sm_medical_category_key, categories) ||
                  IsotopeSearchByEnergy::is_in_category(rctn, IsotopeSearchByEnergy::sm_industrial_category_key, categories) ||
                  IsotopeSearchByEnergy::is_in_category(rctn, IsotopeSearchByEnergy::sm_norm_category_key, categories) ||
                  IsotopeSearchByEnergy::is_in_category(rctn, IsotopeSearchByEnergy::sm_snm_category_key, categories) ||
                  IsotopeSearchByEnergy::is_in_category(rctn, IsotopeSearchByEnergy::sm_common_category_key, categories)
                        // || IsotopeSearchByEnergy::is_in_category(rctn, IsotopeSearchByEnergy::sm_fission_category_key, categories )
                );
              }
            }
          }
        }
        
        if( is_in_any_category )
          characteristic_sources.insert( source_name );
      }
    }
    
    // Print the final set of characteristic sources
    if( !characteristic_sources.empty() )
    {
      cout << "Filtered characteristic sources for potential reference lines (" << characteristic_sources.size() << "): ";
      bool first = true;
      for( const auto& source : characteristic_sources )
      {
        if( !first ) cout << ", ";
        cout << source;
        first = false;
        
      
        //TODO: use following to score potential matches and filter further
        //double profile_weight( std::shared_ptr<const DetectorPeakResponse> detector,
        //                      const std::shared_ptr<const SpecUtils::Measurement> displayed_measurement,
        //                      const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
        //                      const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
        //                      std::vector<SandiaDecay::EnergyRatePair> srcgammas,
        //                      const vector<double> &energies,
        //                      const vector<double> &windows,
        //                      const IsotopeSearchByEnergyModel::IsotopeMatch &nucmatches,
        //                      double shielding_an,
        //                      double shielding_ad )
        
      }
      cout << endl;
    }
    else
    {
      cout << "No characteristic sources passed filtering criteria." << endl;
    }
  }
  


  // Get the FWHM JavaScript function from the foreground detector
  std::string js_fwhm_fcn;
  if( detector && detector->hasResolutionInfo() )
    try{ js_fwhm_fcn = detector->javaScriptFwhmFunction(); }catch(...){ assert(0); } //dont expect it to ever throw
  
  if( js_fwhm_fcn.empty() )//We'll put in a generic fcnt
  {
    const bool highres = PeakFitUtils::is_likely_high_res( m_interspec );
    const vector<float> coefs = highres ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
    js_fwhm_fcn = DetectorPeakResponse::javaScriptFwhmFunction( coefs, DetectorPeakResponse::kGadrasResolutionFcn );
  }//if( js_fwhm_fcn.empty() )
  
  
  // Note: since we are being called from within the render cycle most likely, this next function
  //       call will immediately send the appropriate JS to the client, as long as the chart has been rendered
  //       (i.e., it doesnt scheduleRender to lazily load to client, as long as its rendered).
  m_chart->setKineticRefernceLines( std::move(ref_lines), std::move(js_fwhm_fcn) );
}//void updateLines()
