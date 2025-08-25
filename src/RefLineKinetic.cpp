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
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/MoreNuclideInfo.h"

using namespace std;
using namespace Wt;


namespace
{
  // Helper function to check if two peaks overlap within 1 FWHM
  bool peaks_overlap( const shared_ptr<const PeakDef>& p1, const shared_ptr<const PeakDef>& p2 )
  {
    // See also/instead PeakDef::causilyConnected(...)
    if( !p1 || !p2 )
      return false;
    const double fwhm1 = p1->gausPeak() ? p1->fwhm() : 0.5*p1->roiWidth();
    const double fwhm2 = p2->gausPeak() ? p2->fwhm() : 0.5*p2->roiWidth();
    return fabs(p1->mean() - p2->mean()) <= (0.5*(fwhm1 + fwhm2));
  };
  
  // Function to combine users peaks with the auto-search peaks; if a auto-search peak overlaps with a user peak, it
  //  wont be added.
  vector<shared_ptr<const PeakDef>> combine_nonoverlapping_peaks(
                                            const shared_ptr<const deque<shared_ptr<const PeakDef>>> &user_peaks,
                                            const shared_ptr<const deque<shared_ptr<const PeakDef>>> &autosearch_peaks )
  {
    vector<shared_ptr<const PeakDef>> unique_peaks;
    
    if( user_peaks )
      unique_peaks.insert( end(unique_peaks), begin(*user_peaks), end(*user_peaks) );
    
    if( autosearch_peaks )
    {
      for( const auto& auto_peak : *autosearch_peaks )
      {
        bool overlaps_with_existing = false;
        for( size_t i = 0; !overlaps_with_existing && (i < unique_peaks.size()); ++i )
          overlaps_with_existing = peaks_overlap( auto_peak, unique_peaks[i] );
        
        if( !overlaps_with_existing )
          unique_peaks.push_back( auto_peak );
      }//for( const auto& auto_peak : *autosearch_peaks )
    }//if( autosearch_peaks )
    std::sort( begin(unique_peaks), end(unique_peaks), &PeakDef::lessThanByMeanShrdPtr );

    return unique_peaks;
  }
  
  template <typename T>
  bool is_in_a_category( const T *nuc, const vector<IsotopeSearchByEnergy::NucSearchCategory> &categories )
  {
    return (IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_medical_category_key, categories)
            || IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_industrial_category_key, categories)
            || IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_norm_category_key, categories)
            || IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_snm_category_key, categories)
            || IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_common_category_key, categories)
            // || IsotopeSearchByEnergy::is_in_category(nuc, IsotopeSearchByEnergy::sm_fission_category_key, categories)
            );
  }
}//namespace

struct AlwaysSrcs
{
  const map<string,ReferenceLinePredef::NucMix> nuc_mixes;
  const map<string,ReferenceLinePredef::CustomSrcLines> custom_lines;
  const vector<ReferenceLinePredef::IndividualSource> individual_sources;
  
  AlwaysSrcs() = delete;
  AlwaysSrcs( map<string,ReferenceLinePredef::NucMix> &&mixes,
             map<string,ReferenceLinePredef::CustomSrcLines> &&lines,
             vector<ReferenceLinePredef::IndividualSource> &&srcs )
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
  m_interspec->colorThemeChanged().connect( boost::bind( &RefLineKinetic::colorThemeChanged, this, boost::placeholders::_1 ) );
  
  // TODO: We also need to update lines after an energy calbration is done
  // TODO: We also need to make sure hint peaks are adjusted for energy changes
  
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
    vector<ReferenceLinePredef::IndividualSource> indiv_sources;
    ReferenceLinePredef::load_ref_line_file( always_defs_file, nuc_mixes, custom_lines, &indiv_sources );
    
    m_always_srcs = make_unique<AlwaysSrcs>( std::move(nuc_mixes), std::move(custom_lines), std::move(indiv_sources) );
  }catch( std::exception &e )
  {
    cerr << "Failed to initialize RefLineKinetic: " << e.what() << endl;
    m_init_error_msg = e.what();
  }//try / catch
}//void start_init_always_sources()


void RefLineKinetic::assignColorToInput( ReferenceLineInfo &lines ) const
{
  RefLineInput &input = lines.m_input;
  const string &source_name = input.m_input_txt;
  
  // If color is already set, don't override it
  if( !input.m_color.isDefault() )
    return;
  
  // Check ColorTheme::referenceLineColorForSources first
  const auto color_theme = m_interspec->getColorTheme();
  if( color_theme )
  {
    const auto specific_color_it = color_theme->referenceLineColorForSources.find( source_name );
    if( specific_color_it != color_theme->referenceLineColorForSources.end() && !specific_color_it->second.isDefault() )
    {
      input.m_color = specific_color_it->second;
      return;
    }
  }
  
  // Check if ReferencePhotopeakDisplay can suggest a color
  const auto refphoto = m_interspec->referenceLinesWidget();
  if( refphoto )
  {
    const Wt::WColor suggested_color = refphoto->suggestColorForSource( source_name );
    if( !suggested_color.isDefault() )
    {
      input.m_color = suggested_color;
      return;
    }
  }
  
  // Use the kinetic reference line default color from the ColorTheme
  if( color_theme && !color_theme->kineticRefLineDefaultColor.isDefault() )
  {
    input.m_color = color_theme->kineticRefLineDefaultColor;
  }
  else
  {
    // Fallback to the static default if ColorTheme is unavailable
    input.m_color = Wt::WColor( ColorTheme::sm_kinetic_ref_line_default_color );
  }
}//void RefLineKinetic::assignColorToInput(...)


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
                     const shared_ptr<SpecMeas> &measurement,
                     const set<int> &sample_numbers,
                     const vector<string> &detectors )
{
  if( spec_type == SpecUtils::SpectrumType::Foreground )
    m_external_rid_results.reset();
  
  m_renderFlags |= KineticRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderKineticRefLine();
}//void spectrumChanged(...)


void RefLineKinetic::autoRidResultsRecieved( const shared_ptr<const ExternalRidResults> &results )
{
  m_external_rid_results = results;
  
  m_renderFlags |= KineticRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderKineticRefLine();
}//void autoRidResultsRecieved( const shared_ptr<const ExternalRidResults> &results )


void RefLineKinetic::colorThemeChanged( const shared_ptr<const ColorTheme> &theme )
{
  // Color theme has changed, trigger an update to refresh colors
  m_renderFlags |= KineticRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderKineticRefLine();
}//void RefLineKinetic::colorThemeChanged(...)


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
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    return;

  const bool highres = PeakFitUtils::is_likely_high_res( m_interspec );
  
  shared_ptr<const SpecUtils::Measurement> foreground = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  const set<int> &foreground_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
  const shared_ptr<const SpecMeas> foreground_meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  
  shared_ptr<const SpecUtils::Measurement> background = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
  const set<int> &background_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
  const shared_ptr<const SpecMeas> background_meas = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  
  
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> foreground_autosearch_peaks = foreground_meas
                                  ? foreground_meas->automatedSearchPeaks(foreground_samples)
                                  : nullptr;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> background_autosearch_peaks = background_meas
                                  ? background_meas->automatedSearchPeaks(background_samples)
                                  : nullptr;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> user_background_peaks =
                                  (background_meas && background_meas->sampleNumsWithPeaks().count(background_samples))
                                  ? background_meas->peaks(background_samples)
                                  : nullptr;
  
  PeakModel * const pmodel = m_interspec->peakModel();
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> user_foreground_peaks = pmodel ? pmodel->peaks() : nullptr;

  const double foreground_lt = foreground ? foreground->live_time() : 0.0;
  const double background_lt = background ? background->live_time() : 0.0;
  
  // We'll create a combination of user and auto-search peaks, that dont overlap
  const vector<shared_ptr<const PeakDef>> unique_foreground_peaks
                                  = combine_nonoverlapping_peaks( user_foreground_peaks, foreground_autosearch_peaks );
  const vector<shared_ptr<const PeakDef>> unique_background_peaks
                                  = combine_nonoverlapping_peaks( user_background_peaks, background_autosearch_peaks );
  
  // Get peaks that are in the foreground, but not the background.
  shared_ptr<deque<shared_ptr<const PeakDef>>> nonbackground_peaks = make_shared<deque<shared_ptr<const PeakDef>>>();
  
  for( const shared_ptr<const PeakDef> &p: unique_foreground_peaks )
  {
    bool overlaps_background = false;
    // We will check if a peak overlaps the background peak, and if they do in energy, that the foreground peak is
    //  no bigger than 25% larger than background peak (arbitrarily chosen)
    for( size_t back_index = 0; !overlaps_background && (back_index < unique_background_peaks.size()); ++back_index )
    {
      if( peaks_overlap(p, unique_background_peaks[back_index]) )
      {
        if( (foreground_lt > 0.1) && (background_lt > 0.1) )
        {
          const double fore_rate = p->amplitude() / foreground_lt;
          const double back_rate = unique_background_peaks[back_index]->amplitude() / background_lt;
          overlaps_background = (fore_rate - back_rate) < 0.25*std::max(fore_rate, back_rate);
        }
      }
    }
    if( !overlaps_background )
      nonbackground_peaks->push_back( p );
  }//for( const shared_ptr<const PeakDef> &p: unique_foreground_peaks )
  
  shared_ptr<const DetectorPeakResponse> detector = foreground_meas ? foreground_meas->detector() : nullptr;
  if( !detector || !detector->isValid() )
  {
    // Load a generic detector efficiency function
    shared_ptr<DetectorPeakResponse> detPtr = make_shared<DetectorPeakResponse>();
    
    const string basename = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GenericGadrasDetectors" );
    const string csvfilename = SpecUtils::append_path( basename, highres ? "HPGe 40%/Efficiency.csv" : "NaI 3x3/Efficiency.csv" );
    const string datFilename = SpecUtils::append_path( basename, highres ? "HPGe 40%/Detector.dat" : "NaI 3x3/Detector.dat" );
    
#ifdef _WIN32
    ifstream csv( SpecUtils::convert_from_utf8_to_utf16(csvfilename).c_str(), ios_base::binary|ios_base::in );
    ifstream datFile( SpecUtils::convert_from_utf8_to_utf16(datFilename).c_str(), ios_base::binary|ios_base::in );
#else
    ifstream csv( csvfilename.c_str(), ios_base::binary|ios_base::in );
    ifstream datFile( datFilename.c_str(), ios_base::binary|ios_base::in );
#endif
    
    if( csv.good() && datFile.good() )
    {
      detPtr->fromGadrasDefinition( csv, datFile );
    }else
    {
      cerr << "findCandidates(...): error opening default detector file" << endl;
      detPtr->setIntrinsicEfficiencyFormula( "1.0", 3.0*PhysicalUnits::cm, PhysicalUnits::keV,
                                            0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarField );
    }
    
    detector = detPtr;
  }//if( !detector || !detector->isValid() )
  
  assert( detector );
  if( !detector->hasResolutionInfo() )
  {
    // TODO: fit this from unique_foreground_peaks, maybe using DetectorPeakResponse::fitResolution(...), or MakeDrfFit::performResolutionFit(...);
  }//if( !detector->hasResolutionInfo() )
  
  
  const ReactionGamma * const reaction_db = ([]()-> const ReactionGamma *{
    try{ return ReactionGammaServer::database(); }catch(...){}
    return static_cast<const ReactionGamma *>(nullptr); }
  )();
  
  
  RefLineInput base_input;
  base_input.m_lower_br_cutt_off = 0.0;
  base_input.m_promptLinesOnly = false;
  base_input.m_showGammas = true;
  base_input.m_showXrays = true;
  base_input.m_showAlphas = false;
  base_input.m_showBetas = false;
  base_input.m_showCascades = false;
  base_input.m_showEscapes = false;
  
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
        const string &shielding_material = *src.shielding_material;
        
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
      }else
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
      }else
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
  if( user_foreground_peaks && user_foreground_peaks->size() )
  {
    // We will assign these nuclides a weight of 10
    for( const auto& peak : *user_foreground_peaks )
    {
      if( !peak )
        continue;
        
      string canonical_name;
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
        
      string canonical_name;
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
        
      string canonical_name;
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
  
  
  // Add associated nuclides for existing reference lines
  shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> more_nuc_info = MoreNuclideInfo::MoreNucInfoDb::instance();
  if( more_nuc_info )
  {
    // Create a copy of current ref_lines to iterate over (to avoid modifying collection while iterating)
    vector<pair<double, ReferenceLineInfo>> current_ref_lines = ref_lines;
    
    for( const pair<double, ReferenceLineInfo> &ref_line_pair : current_ref_lines )
    {
      const double parent_weight = ref_line_pair.first;
      const string& ref_line_name = ref_line_pair.second.m_input.m_input_txt;
      
      const MoreNuclideInfo::NucInfo *nuc_info = more_nuc_info->info(ref_line_name);
      if( !nuc_info || nuc_info->m_associated.empty() )
        continue;
      
      for( const string &associated_name : nuc_info->m_associated )
      {
        // Check if this associated nuclide is already in ref_lines
        bool associated_already_exists = false;
        double associated_weight = 0.25 * parent_weight;
        
        for( auto& existing : ref_lines )
        {
          if( existing.second.m_input.m_input_txt == associated_name )
          {
            associated_already_exists = true;
            // Use max of current weight and 0.25 times parent weight
            existing.first = std::max(existing.first, associated_weight);
            break;
          }
        }//for( auto& existing : ref_lines )
        
        if( !associated_already_exists )
        {
          // Create reference lines for associated nuclide
          RefLineInput associated_input = base_input;
          associated_input.m_input_txt = associated_name;
          
          shared_ptr<ReferenceLineInfo> associated_ref_info = ReferenceLineInfo::generateRefLineInfo( associated_input );
          if( associated_ref_info && associated_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
          {
            filterLines( *associated_ref_info, RefLineSrc::CharacteristicLine, foreground );
            ref_lines.emplace_back( associated_weight, std::move(*associated_ref_info) );
          }
          else
          {
            cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for associated nuclide: " << associated_name << endl;
          }
        }//if( !associated_already_exists )
      }//for( const string &associated_name : nuc_info->m_associated )
    }//for( ref_line_pair in current_ref_lines )
  }//if( more_nuc_info )
  

  // For peaks above ~1460 keV, either user or auto-fit, add in escape peak lines
  // Add in escape peak lines, for peaks above a practical pair-production threshold
  const double pair_prod_thresh = highres ? 1255.0 : 2585; //The single_escape_sf and double_escape_sf give negative values below 1255; the value used for low-res is arbitrary
   //The efficiency of S.E. and D.E. peaks, relative to F.E. peak, for the 20% Generic GADRAS DRF
  //  included in InterSpec, is given pretty well by the following (energy in keV):
  const auto single_escape_sf = []( const double x ) -> double {
    return std::max( 0.0, (1.8768E-11 *x*x*x) - (9.1467E-08 *x*x) + (2.1565E-04 *x) - 0.16367 );
  };
  
  const auto double_escape_sf = []( const double x ) -> double {
    return std::max( 0.0, (1.8575E-11 *x*x*x) - (9.0329E-08 *x*x) + (2.1302E-04 *x) - 0.16176 );
  };
  
  vector<ReferenceLineInfo::RefLine> escape_lines;
  for( const shared_ptr<const PeakDef> &peak : unique_foreground_peaks )
  {
    if( !peak || peak->mean() < pair_prod_thresh )
      continue;
      
    const double peak_energy = peak->mean();
    const double peak_amplitude = peak->amplitude();
    
    // Single escape peak at energy - 511 keV
    const double se_energy = peak_energy - 510.9989;
    const double se_amplitude = peak_amplitude * single_escape_sf( peak_energy );
    const double se_peak_amp_threshold = 3.0; //TODO: base this more intelligently off of the spectrum itself
    if( se_amplitude >= se_peak_amp_threshold )
    {
      ReferenceLineInfo::RefLine se_line;
      se_line.m_energy = se_energy;
      se_line.m_normalized_intensity = se_amplitude;
      se_line.m_decaystr = "S.E.";
      se_line.m_parent_nuclide = peak->parentNuclide();
      se_line.m_transition = peak->nuclearTransition();
      se_line.m_reaction = peak->reaction();
      se_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::SingleEscape;
      escape_lines.push_back( se_line );
    }
    
    // Double escape peak at energy - 1022 keV
    const double de_energy = peak_energy - 2.0 * 510.9989;
    const double de_amplitude = peak_amplitude * double_escape_sf( peak_energy );
    const double de_peak_amp_threshold = 3.0; //TODO: base this more intelligently off of the spectrum itself
    if( de_amplitude >= de_peak_amp_threshold )
    {
      ReferenceLineInfo::RefLine de_line;
      de_line.m_energy = de_energy;
      de_line.m_normalized_intensity = de_amplitude;
      de_line.m_decaystr = "D.E.";
      de_line.m_parent_nuclide = peak->parentNuclide();
      de_line.m_transition = peak->nuclearTransition();
      de_line.m_reaction = peak->reaction();
      de_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape;
      escape_lines.push_back( de_line );
    }
  }//for( const auto& peak : unique_foreground_peaks )
  
  // Normalize escape lines so highest intensity is 1.0 and add to ref_lines
  if( !escape_lines.empty() )
  {
    double max_intensity = 0.0;
    for( const auto& line : escape_lines )
      max_intensity = std::max( max_intensity, line.m_normalized_intensity );
      
    if( max_intensity > 0.0 )
    {
      for( auto& line : escape_lines )
        line.m_normalized_intensity /= max_intensity;
        
      ReferenceLineInfo escape_info;
      escape_info.m_ref_lines = std::move( escape_lines );
      escape_info.m_input.m_input_txt = "Escape Peaks";
      //escape_info.m_input.m_color = ...;
      escape_info.m_validity = ReferenceLineInfo::InputValidity::Valid;
      
      ref_lines.emplace_back( 1.0, std::move(escape_info) );
    }
  }

  // Check if user or auto-fit peaks cooresponds to a cascade sum energy... maybe require the parent nuclide to not have a gamma at that energy, but how often does this happen?
  
  // We'll get the peaks that dont have a source assigned to them, and dont have a background peak, and then
  //  use those to find suggested characteristic - but we arent currently adding these lines;
  deque<shared_ptr<const PeakDef>> peaks_for_characteristics;
  for( const shared_ptr<const PeakDef> &peak : unique_foreground_peaks )
  {
    if( !peak || peak->parentNuclide() || peak->reaction() || peak->xrayElement() )
      continue;
    
    bool should_skip = false;
    for( size_t i = 0; !should_skip && (i < unique_background_peaks.size()); ++i )
      should_skip = peaks_overlap(peak, unique_background_peaks[i]);
    
    if( !should_skip )
      peaks_for_characteristics.push_back(peak);
  }//for( const shared_ptr<const PeakDef> &peak : unique_foreground_peaks )
  std::sort( begin(peaks_for_characteristics), end(peaks_for_characteristics), &PeakDef::lessThanByMeanShrdPtr );
    

  // Process automated search peaks to find characteristic nuclides
  // If we want to be heavy-handed, we could use (some of the) suggestions from populateCandidateNuclides
  if( !peaks_for_characteristics.empty() )
  {
    // This next section is experimental - it tries to match peaks to characteristic lines,
    // however, there is a concern that we might end up with just too many ref_lines and it
    // will be detrimental to the experience
    set<string> candidate_source_names;
    map<string,vector<shared_ptr<const PeakDef>>> candidate_to_peaks;
    
    cout << "Finding characteristics for " << peaks_for_characteristics.size() << " remaining automated search peaks:" << endl;
    for( const shared_ptr<const PeakDef> &peak : peaks_for_characteristics )
    {
      vector<string> characteristicnucs;
      //IsotopeId::findCharacteristics( characteristicnucs, peak );
      //IsotopeId::findCandidates( characteristicnucs, peak, nonbackground_peaks, detector, foreground );
      
      IsotopeId::PeakToNuclideMatch suggestedNucs;
      IsotopeId::suggestNuclides( suggestedNucs, peak, nonbackground_peaks, foreground, detector );
      
      const vector<IsotopeId::NuclideStatWeightPair> &sugestions = suggestedNucs.nuclideWeightPairs;
      for( size_t i = 0; i < 10 && i < sugestions.size(); ++i )
      {
        if( sugestions[i].nuclide )
        {
          if( sugestions[i].weight > 0.1*sugestions[0].weight ) //10% is arbitrary, but just looking through a few examples seems reasonable
          {
            characteristicnucs.push_back(sugestions[i].nuclide->symbol );
            candidate_to_peaks[sugestions[i].nuclide->symbol].push_back(peak);
          }
          //cout << "suggestNuclides(" << peak->mean() << ", " << i << ") gave " << sugestions[i].nuclide->symbol << ", with weight " << sugestions[i].weight << endl;
        }
      }
      
      cout << "Peak at " << peak->mean() << " keV: ";
      if( characteristicnucs.empty() )
      {
        cout << "No characteristic nuclides found" << endl;
      }else
      {
        cout << "Found " << characteristicnucs.size() << " candidates: ";
        
        // Collect top 5 characteristics
        for( size_t i = 0; i < characteristicnucs.size() && i < 5; ++i )
        {
          if( i > 0 ) cout << ", ";
          cout << characteristicnucs[i];
          
          candidate_source_names.insert( characteristicnucs[i] );
        }
        
        if( characteristicnucs.size() > 5 )
          cout << " (and " << (characteristicnucs.size() - 5) << " more)";
        cout << endl;
      }
    }//for( const shared_ptr<const PeakDef> &peak : peaks_for_characteristics )
  
    
    cout << "Starting candidate_source_names={";
    for( auto name : candidate_source_names )
      cout << name << ", ";
    cout << "}" << endl;
    
    // Now filter the collected source names once
    set<string> characteristic_sources;
    
    if( !candidate_source_names.empty() )
    {
      IsotopeSearchByEnergy * const search = m_interspec->nuclideSearch();
      const vector<IsotopeSearchByEnergy::NucSearchCategory> categories = search ? search->search_categories() : vector<IsotopeSearchByEnergy::NucSearchCategory>{};
      
      for( const string &source_name : candidate_source_names )
      {
        // Check if already in ref_lines
        bool already_in_ref_lines = false;
        for( size_t index = 0; !already_in_ref_lines && (index < ref_lines.size()); ++index )
          already_in_ref_lines = SpecUtils::iequals_ascii(ref_lines[index].second.m_input.m_input_txt, source_name);
        
        if( already_in_ref_lines || categories.empty() )
        {
          cout << "Skipping cadndiate source '" << source_name << "' since its already in ref_lines" << endl;
          continue;
        }
        
        bool is_in_any_category = false;
        
        // Check if it's a nuclide in any category
        const SandiaDecay::Nuclide *nuc = db->nuclide( source_name );
        if( nuc )
        {
          // Filter out any nuclide that is a descendant of any nuclide already in ref_lines
          bool is_descendant = false;
          for( size_t index = 0; !is_descendant && (index < ref_lines.size()); ++index )
          {
            const SandiaDecay::Nuclide *existing_nuc = db->nuclide( ref_lines[index].second.m_input.m_input_txt );
            is_descendant = (existing_nuc && (nuc->branchRatioFromForebear(existing_nuc) > 0.0));
          }
          
          if( !is_descendant )
          {
            is_in_any_category = is_in_a_category(nuc,categories);
            if( !is_in_any_category )
              cout << "Source '" << source_name << "' is not in NuclideSearchCatagories.xml" << endl;
          }else
          {
            cout << "Source '" << source_name << "' is a decendant of an existing ref line" << endl;
          }
        }else
        {
          //Currently, we wont actually get here, because the `candidate_source_names` are currently all nuclides.
          
          // Check if it's an element
          const SandiaDecay::Element *el = db->element( source_name );
          if( el )
          {
            is_in_any_category = is_in_a_category(el,categories);
          }else if( reaction_db )
          {
            // Check if it's a reaction
            vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
            reaction_db->gammas( source_name, possible_rctns );
            
            if( !possible_rctns.empty() )
            {
              const ReactionGamma::Reaction *rctn = possible_rctns[0].reaction;
              
              // Check if at least 25% of reaction abundance has matching peaks
              const vector<ReactionGamma::Reaction::EnergyYield> &gammas = rctn->gammas;
              const double required_percentage = 0.25;  // 25%
              
              // Calculate total abundance and matched abundance
              double total_abundance = 0.0, matched_abundance = 0.0;
              for( const ReactionGamma::Reaction::EnergyYield &gamma : gammas )
              {
                total_abundance += gamma.abundance;
                
                // Check if this gamma has a matching peak
                for( const shared_ptr<const PeakDef> &peak : unique_foreground_peaks )
                {
                  const double peak_energy = peak->mean();
                  const double fwhm = peak->fwhm();
                  const double tolerance = fwhm > 0 ? fwhm : 5.0;  // Use FWHM or 5 keV default
                  
                  if( fabs(peak_energy - gamma.energy) <= tolerance )
                  {
                    matched_abundance += gamma.abundance;
                    break;  // Found a match for this gamma, move to next gamma
                  }//if( matches )
                }//for( loop over unique_foreground_peaks )
              }//for( loop over rctn->gammas )
              
              const double abundance_percentage = (total_abundance > 0.0) ? (matched_abundance / total_abundance) : 0.0;
              is_in_any_category = (abundance_percentage >= required_percentage) && is_in_a_category(rctn,categories);
            }//if( !possible_rctns.empty() )
          }//if( el ) / else if( reaction_db )
        }//if( nuc ) / else
        
        if( is_in_any_category )
          characteristic_sources.insert( source_name );
        else
          cout << "Eliminating '" << source_name << "'" << endl;
      }//for( const string &source_name : candidate_source_names )
    }//if( !candidate_source_names.empty() )
    
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
        const vector<shared_ptr<const PeakDef>> user_peaks = user_foreground_peaks
             ? vector<shared_ptr<const PeakDef>>{begin(*user_foreground_peaks), end(*user_foreground_peaks)}
             : vector<shared_ptr<const PeakDef>>{};
        const vector<shared_ptr<const PeakDef>> automated_search_peaks = foreground_autosearch_peaks
             ? vector<shared_ptr<const PeakDef>>{begin(*foreground_autosearch_peaks), end(*foreground_autosearch_peaks)}
             : vector<shared_ptr<const PeakDef>>{};
        const shared_ptr<const deque<shared_ptr<const PeakDef>>> user_foreground_peaks = pmodel ? pmodel->peaks() : nullptr;
        
        const double atomic_nums[]   = { 1.0, 26.0, 74.0 };
        const double areal_density[] = { 0.0*PhysicalUnits::g_per_cm2, 10.0*PhysicalUnits::g_per_cm2, 25.0*PhysicalUnits::g_per_cm2 };
        
        double prof_weight = -999.9;
        
        vector<shared_ptr<const PeakDef>> peaks_for_candidate = candidate_to_peaks[source];
        
        // Remove any non-Gaussian peaks (we cant easily call `amplitude()` for them, so we'll ignore them for the moment.
        peaks_for_candidate.erase( std::remove_if( begin(peaks_for_candidate), end(peaks_for_candidate),
          [](const shared_ptr<const PeakDef> &p ){
            return !p || !p->gausPeak();
        }), end(peaks_for_candidate) );
        
        // Make sure all the peaks in `peaks_for_candidate` are unique (they should be)
        std::sort( begin(peaks_for_candidate), end(peaks_for_candidate),
          []( const shared_ptr<const PeakDef> &lhs, const shared_ptr<const PeakDef> &rhs ){
            return lhs.get() < rhs.get();
        } );
        peaks_for_candidate.erase( std::unique(begin(peaks_for_candidate), end(peaks_for_candidate)), end(peaks_for_candidate) );
        //Sort `peaks_for_candidate` by amplitude
        std::sort( begin(peaks_for_candidate), end(peaks_for_candidate),
                  []( const shared_ptr<const PeakDef> &lhs, const shared_ptr<const PeakDef> &rhs ){
          return lhs->amplitude() > rhs->amplitude();
        } );
        
        vector<double> energies, windows;
        const size_t max_peaks_per_source = 3;
        for( size_t i = 0; (i < max_peaks_per_source) && (i < peaks_for_candidate.size()); ++i )
        {
          energies.push_back( peaks_for_candidate[i]->mean() );
          windows.push_back( 2.5*peaks_for_candidate[i]->sigma() );
        }
        
        vector<SandiaDecay::EnergyRatePair> srcgammas;
        if( const SandiaDecay::Nuclide * const nuc = db->nuclide(source) )
        {
          SandiaDecay::NuclideMixture mix;
          mix.addNuclideByActivity( nuc, 0.001*PhysicalUnits::curie );
          srcgammas = mix.photons( PeakDef::defaultDecayTime(nuc) );
        }else if( const SandiaDecay::Element *el = db->element(source) )
        {
          for( const SandiaDecay::EnergyIntensityPair &p : el->xrays )
            srcgammas.emplace_back( p.intensity, p.energy );
        }else if( reaction_db )
        {
          try
          {
            vector<ReactionGamma::ReactionPhotopeak> rctn_lines;
            reaction_db->gammas(source, rctn_lines );
            for( const ReactionGamma::ReactionPhotopeak &line : rctn_lines )
              srcgammas.emplace_back( line.abundance, line.energy );
          }catch( std::exception &e )
          {
          }
        }//if(nuc) / else if( el ) / else if( rctn )
        
        // We dont need to normalize the intensities, I dont think
        //double max_rate = -9999.9;
        //for( const SandiaDecay::EnergyRatePair &erp : srcgammas )
        //  max_rate = std::max( max_rate, erp.numPerSecond );
        //if( max_rate <= 0.0 )
        //  srcgammas.clear();
        //for( SandiaDecay::EnergyRatePair &erp : srcgammas )
        //  erp.numPerSecond /= max_rate;
        
        for( size_t i = 0; !srcgammas.empty() && (i < 3); ++i )
        {
          const double weight = IsotopeId::profile_weight( detector, foreground, user_peaks, automated_search_peaks, srcgammas,
                                energies, windows, atomic_nums[i], areal_density[i], WString() );
          prof_weight = std::max( prof_weight, weight );
        }
        
        // A min allowed profile weight of 0.25 is arbitrary, but if anything, we could move it higher, to like 0.5
        const double min_allowed_profile_weight = 0.25;
        if( prof_weight < min_allowed_profile_weight )
        {
          cout << "Eliminating '" << source << "' from profile weight only " << prof_weight << endl;
          continue;
        }
        
        cout << "For source " << source << " the profile weight is " << prof_weight << endl;
        
        // Create reference lines for associated nuclide
        RefLineInput candidate_input = base_input;
        candidate_input.m_input_txt = source;
        const double candidate_weight = 0.1 + 4.9 * (2.0 * std::max(0.0, prof_weight - 0.5)); //
        
        shared_ptr<ReferenceLineInfo> candidate_ref_info = ReferenceLineInfo::generateRefLineInfo( candidate_input );
        if( candidate_ref_info && candidate_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
        {
          filterLines( *candidate_ref_info, RefLineSrc::CharacteristicLine, foreground );
          ref_lines.emplace_back( candidate_weight, std::move(*candidate_ref_info) );
        }else
        {
          cerr << "RefLineKinetic::updateLines(): Failed to generate valid reference line info for candidate nuclide: " << source << endl;
        }
      }
      
      cout << endl;
    }
    else
    {
      cout << "No characteristic sources passed filtering criteria." << endl;
    }
  }
  
  
  // If the deadtime is over 15% (low resolution) or 25% (HPGe), add in lines for random sum peaks of largest peaks
  const double max_live_time_frac = 1.0; //highres ? 0.75 : 0.85;
  const double rt = foreground ? foreground->real_time() : 0.0;
  const double lt = foreground ? foreground->live_time() : 0.0;
  if( (lt > 0.0) && (rt > 0.0) && (lt < max_live_time_frac*rt) )
  {
    const size_t max_summing_peaks = 8;
    const double rel_sum_prob_min = 0.01;
    
    // Get the highest amplitude peaks from unique_foreground_peaks, only considering Gaussian peaks
    vector<shared_ptr<const PeakDef>> summing_candidates;
    for( const shared_ptr<const PeakDef> &peak : unique_foreground_peaks )
    {
      if( peak && peak->gausPeak() && peak->amplitude() > 0.0 )
        summing_candidates.push_back( peak );
    }
    
    // Sort by amplitude (highest first)
    std::sort( begin(summing_candidates), end(summing_candidates),
              []( const shared_ptr<const PeakDef> &lhs, const shared_ptr<const PeakDef> &rhs ){
                return lhs->amplitude() > rhs->amplitude();
              });
    
    if( summing_candidates.size() > max_summing_peaks )
      summing_candidates.resize( max_summing_peaks );
    
    if( !summing_candidates.empty() )
    {
      // Calculate random-summing probabilities for all pairs
      vector<tuple<double,double,double,double>> sum_peaks; // energy, rel probability, energy_i, energy_j
      
      // Calculate maximum possible probability (highest peak with itself)
      const double max_sum_prob = summing_candidates[0]->amplitude() * summing_candidates[0]->amplitude();
      
      // Generate all unique pairs (including peak with itself)
      for( size_t i = 0; i < summing_candidates.size(); ++i )
      {
        for( size_t j = i; j < summing_candidates.size(); ++j )
        {
          const shared_ptr<const PeakDef> &peak_i = summing_candidates[i];
          const shared_ptr<const PeakDef> &peak_j = summing_candidates[j];
          const double sum_energy = peak_i->mean() + peak_j->mean();
          const double sum_prob = peak_i->amplitude() * peak_j->amplitude();
          const double rel_intensity = (sum_prob / max_sum_prob);
          if( rel_intensity >= rel_sum_prob_min )
            sum_peaks.emplace_back( sum_energy, rel_intensity, peak_i->mean(), peak_j->mean() );
        }//for( size_t j = i; j < summing_candidates.size(); ++j )
      }//for( size_t i = 0; i < summing_candidates.size(); ++i )
      
      // Create reference lines for significant sum peaks
      if( !sum_peaks.empty() )
      {
        ReferenceLineInfo sum_ref_info;
        sum_ref_info.m_input.m_input_txt = "Random Sum Peaks";
        //sum_ref_info.m_input.m_color = ...
        sum_ref_info.m_validity = ReferenceLineInfo::InputValidity::Valid;
        sum_ref_info.m_source_type = ReferenceLineInfo::SourceType::CustomEnergy;
        
        for( const auto &sum_peak : sum_peaks )
        {
          const double energy = std::get<0>( sum_peak );
          const double rel_intensity = std::get<1>( sum_peak );
          
          ReferenceLineInfo::RefLine ref_line;
          ref_line.m_energy = energy;
          ref_line.m_normalized_intensity = rel_intensity;
          ref_line.m_decay_intensity = rel_intensity;
          ref_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
          ref_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak;
          ref_line.m_attenuation_applies = false;
          char buffer[64] = { '\0' };
          snprintf( buffer, sizeof(buffer), "Random sum %.1f + %.1f keV", get<2>(sum_peak), get<3>(sum_peak) );
          ref_line.m_decaystr = buffer;
          
          sum_ref_info.m_ref_lines.push_back( ref_line );
        }//for( const auto &sum_peak : sum_peaks )
        
        const double randum_sum_lines_weight = 1.0; //totally arbitrary
        ref_lines.emplace_back( randum_sum_lines_weight, std::move(sum_ref_info) );
      }//if( !sum_peaks.empty() )
    }//if( summing_candidates.size() >= 2 )
  }//if( (lt > 0.0) && (rt > 0.0) && (lt < max_live_time_frac*rt) )


  // Get the FWHM JavaScript function from the foreground detector
  string js_fwhm_fcn;
  if( detector && detector->hasResolutionInfo() )
    try{ js_fwhm_fcn = detector->javaScriptFwhmFunction(); }catch(...){ assert(0); } //dont expect it to ever throw
  
  if( js_fwhm_fcn.empty() )//We'll put in a generic fcnt
  {
    const vector<float> coefs = highres ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
    js_fwhm_fcn = DetectorPeakResponse::javaScriptFwhmFunction( coefs, DetectorPeakResponse::kGadrasResolutionFcn );
  }//if( js_fwhm_fcn.empty() )
  
  
  for( pair<double,ReferenceLineInfo> &ref_line : ref_lines )
    assignColorToInput( ref_line.second );
  
  // Note: since we are being called from within the render cycle most likely, this next function
  //       call will immediately send the appropriate JS to the client, as long as the chart has been rendered
  //       (i.e., it doesnt scheduleRender to lazily load to client, as long as its rendered).
  m_chart->setKineticRefernceLines( std::move(ref_lines), std::move(js_fwhm_fcn) );
}//void updateLines()
