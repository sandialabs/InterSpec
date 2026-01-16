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
#include <chrono>
#include <Wt/WServer>
#include <Wt/WIOService>
#include <Wt/WApplication>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RefLineDynamic.h"
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
  
  template <typename T>
  Wt::WColor getCategoryColorForSource( const shared_ptr<const ColorTheme> &color_theme,
                                        const T *source,
                                        const vector<IsotopeSearchByEnergy::NucSearchCategory> &search_categories )
  {
    if( !color_theme || !source )
      return color_theme ? color_theme->dynamicRefLineOtherColor : Wt::WColor( ColorTheme::sm_dynamic_ref_line_other_color );
    
    const std::string &snm_key = IsotopeSearchByEnergy::sm_snm_category_key;
    const std::string &industrial_key = IsotopeSearchByEnergy::sm_industrial_category_key;
    const std::string &medical_key = IsotopeSearchByEnergy::sm_medical_category_key;
    const std::string &norm_key = IsotopeSearchByEnergy::sm_norm_category_key;
    const std::string &common_key = IsotopeSearchByEnergy::sm_common_category_key;
    
    if( IsotopeSearchByEnergy::is_in_category( source, snm_key, search_categories ) )
      return color_theme->dynamicRefLineSnmColor;
    
    if( IsotopeSearchByEnergy::is_in_category( source, industrial_key, search_categories ) )
      return color_theme->dynamicRefLineIndustrialColor;
    
    if( IsotopeSearchByEnergy::is_in_category( source, medical_key, search_categories ) )
      return color_theme->dynamicRefLineMedicalColor;
    
    if( IsotopeSearchByEnergy::is_in_category( source, norm_key, search_categories ) )
      return color_theme->dynamicRefLineNormColor;
    
    if( IsotopeSearchByEnergy::is_in_category( source, common_key, search_categories ) )
      return color_theme->dynamicRefLineCommonColor;
    
    // If no specific category found, use "Other" category
    return color_theme->dynamicRefLineOtherColor;
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



RefLineDynamic::RefLineDynamic( D3SpectrumDisplayDiv *chart, InterSpec *interspec )
  : Wt::WObject( chart ),
  m_interspec( interspec ),
  m_chart( chart ),
  m_active( false ),
  m_has_inited( false ),
  m_init_error_msg{},
  m_always_srcs{},
  m_external_rid_results( nullptr ),
  m_renderFlags{},
  m_current_calc_num( make_shared<atomic<size_t>>(0) ),
  m_current_ref_lines( nullptr )
{
  if( !m_interspec || !m_chart )
    throw std::runtime_error( "RefLineDynamic: null InterSpec parent or chart" );
  
  m_active = UserPreferences::preferenceValue<bool>( "DynamicRefLine", m_interspec );
  m_interspec->preferences()->addCallbackWhenChanged( "DynamicRefLine", this, &RefLineDynamic::setActive );
    
  m_interspec->hintPeaksSet().connect( boost::bind(&RefLineDynamic::autoSearchPeaksSet, this, boost::placeholders::_1) );
  m_interspec->displayedSpectrumChanged().connect( boost::bind(&RefLineDynamic::spectrumChanged, this,
    boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4
  ) );
  
  m_interspec->externalRidResultsRecieved().connect( boost::bind( &RefLineDynamic::autoRidResultsRecieved, this, boost::placeholders::_1 ) );
  m_interspec->colorThemeChanged().connect( boost::bind( &RefLineDynamic::colorThemeChanged, this, boost::placeholders::_1 ) );
  PeakModel * const pmodel = m_interspec->peakModel();
  if( pmodel )
  {
    pmodel->rowsInserted().connect( boost::bind( &RefLineDynamic::peaksAdded, this ) );
    pmodel->rowsRemoved().connect( boost::bind( &RefLineDynamic::peaksRemoved, this ) );
    pmodel->dataChanged().connect( boost::bind( &RefLineDynamic::peakModified, this ) );
  }

  // TODO: We also need to update lines after an energy calbration is done
  // TODO: We also need to make sure hint peaks are adjusted for energy changes
  
  m_chart->setDynamicRefLineController( this );
  
  start_init_always_sources();
}//RefLineDynamic constructor


RefLineDynamic::~RefLineDynamic()
{
}


void RefLineDynamic::start_init_always_sources()
{
  try
  {
    // TODO: move this to being done in a background thread
    m_has_inited = true;
    
    string always_defs_file;
    
    try
    {
      always_defs_file = SpecUtils::append_path(InterSpec::writableDataDirectory(), "dynamic_ref_lines.xml");
    }catch( std::exception &e )
    {
      //writableDataDirectory not set
    }
    
    if( always_defs_file.empty() || !SpecUtils::is_file(always_defs_file) )
      always_defs_file = SpecUtils::append_path(InterSpec::staticDataDirectory(), "dynamic_ref_lines.xml");
    
    map<string,ReferenceLinePredef::NucMix> nuc_mixes;
    map<string,ReferenceLinePredef::CustomSrcLines> custom_lines;
    vector<ReferenceLinePredef::IndividualSource> indiv_sources;
    ReferenceLinePredef::load_ref_line_file( always_defs_file, nuc_mixes, custom_lines, &indiv_sources );
    
    m_always_srcs = make_shared<AlwaysSrcs>( std::move(nuc_mixes), std::move(custom_lines), std::move(indiv_sources) );
  }catch( std::exception &e )
  {
    cerr << "Failed to initialize RefLineDynamic: " << e.what() << endl;
    m_init_error_msg = e.what();
  }//try / catch
}//void start_init_always_sources()


void RefLineDynamic::assignColorToInput( ReferenceLineInfo &lines ) const
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
  
  // Determine category-based color for dynamic reference lines
  if( color_theme )
  {
    Wt::WColor category_color;
    
    // Get the search categories to classify the source
    const auto nuclide_search = m_interspec->nuclideSearch();
    if( !nuclide_search )
    {
      input.m_color = Wt::WColor( ColorTheme::sm_dynamic_ref_line_other_color );
      return;
    }
    const auto &search_categories = nuclide_search->search_categories();
    
    // We need to check if any of the reference lines match categories
    // For now, we'll classify based on the primary source from the lines
    if( !lines.m_ref_lines.empty() )
    {
      const auto &first_line = lines.m_ref_lines.front();
      
      // Check if it's a nuclide source
      if( first_line.m_parent_nuclide )
        category_color = getCategoryColorForSource( color_theme, first_line.m_parent_nuclide, search_categories );
      else if( first_line.m_element )// Check if it's an element source (x-ray)
        category_color = getCategoryColorForSource( color_theme, first_line.m_element, search_categories );
      else if( first_line.m_reaction ) // Check if it's a reaction source
        category_color = getCategoryColorForSource( color_theme, first_line.m_reaction, search_categories );
    }//if( !lines.m_ref_lines.empty() )
    
    // If no specific category found, use "Other" category
    if( category_color.isDefault() )
      category_color = color_theme->dynamicRefLineOtherColor;
    
    // Set the color if valid, otherwise use fallback
    input.m_color = category_color.isDefault() ? WColor( ColorTheme::sm_dynamic_ref_line_other_color ) : category_color;
  }else
  {
    // Fallback to the static default if ColorTheme is unavailable
    input.m_color = Wt::WColor( ColorTheme::sm_dynamic_ref_line_other_color );
  }
}//void RefLineDynamic::assignColorToInput(...)


void RefLineDynamic::peaksAdded()
{
  if( !m_active )
    return;

  m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderDynamicRefLine();
  m_current_ref_lines.reset();
}


void RefLineDynamic::peaksRemoved()
{
  if( !m_active )
    return;

  m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderDynamicRefLine();
  m_current_ref_lines.reset();
}


void RefLineDynamic::peakModified()
{
  if( !m_active )
    return;

  m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderDynamicRefLine();
  m_current_ref_lines.reset();
}


void RefLineDynamic::setActive( bool active )
{
  if( m_active == active )
    return;
  
  m_active = active;

  m_current_ref_lines.reset();
  if( m_active )
  {
    m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
    m_chart->scheduleRenderDynamicRefLine();
  }else
  {
    m_chart->setDynamicRefernceLines( vector<pair<double,ReferenceLineInfo>>{}, "" );
  }
}//void RefLineDynamic::setActive( bool active )


bool RefLineDynamic::isActive() const
{
  return m_active;
}

shared_ptr<vector<pair<double,ReferenceLineInfo>>> RefLineDynamic::current_lines() const
{
  return m_current_ref_lines;
}//shared_ptr<vector<pair<double,ReferenceLineInfo>>> current_lines() const


void RefLineDynamic::autoSearchPeaksSet( const SpecUtils::SpectrumType spectrum )
{
  if( !m_active )
    return;
  
  m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderDynamicRefLine();
  m_current_ref_lines.reset();
}//void autoSearchPeaksSet(...)


void RefLineDynamic::spectrumChanged( const SpecUtils::SpectrumType spec_type,
                     const shared_ptr<SpecMeas> &measurement,
                     const set<int> &sample_numbers,
                     const vector<string> &detectors )
{
  if( spec_type == SpecUtils::SpectrumType::Foreground )
    m_external_rid_results.reset();
  
  m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderDynamicRefLine();
  m_current_ref_lines.reset();
}//void spectrumChanged(...)


void RefLineDynamic::autoRidResultsRecieved( const shared_ptr<const ExternalRidResults> &results )
{
  m_external_rid_results = results;
  
  m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderDynamicRefLine();
  m_current_ref_lines.reset();
}//void autoRidResultsRecieved( const shared_ptr<const ExternalRidResults> &results )


void RefLineDynamic::colorThemeChanged( const shared_ptr<const ColorTheme> &theme )
{
  // Color theme has changed, trigger an update to refresh colors
  m_renderFlags |= DynamicRefLineRenderFlags::UpdateLines;
  m_chart->scheduleRenderDynamicRefLine();
  m_current_ref_lines.reset();
}//void RefLineDynamic::colorThemeChanged(...)


void RefLineDynamic::startPushUpdates()
{
  if( m_renderFlags.testFlag(DynamicRefLineRenderFlags::UpdateLines) )
    startUpdateLines();
  
  m_renderFlags = 0;
}//void pushUpdates()


void RefLineDynamic::filterLines( ReferenceLineInfo &ref_lines,
                                 const RefLineDynamic::RefLineSrc src,
                                 const shared_ptr<const SpecUtils::Measurement> &meas,
                                 const std::shared_ptr<const DetectorPeakResponse> &detector )
{
  // This function tries to reduce the raw amount of information, and number of lines that
  // potentually get loaded to the client for the dynamic reference lines.
  //
  // This function currently filters the reflines by:
  // - Calculate a simple `importance = yield(i)*sqrt(energy(i))/sum(yield*sqrt(energy))`
  // - Order RefLines by this importance.
  // - derive an incredible rough, super terrible, "activity" estimate using the max value from the first 6
  //   most "important" lines.
  // - scale other lines by that "activity", and compute a (again, incredible rough, super terrible)
  //   significance estimate for its prediced counts relative to data at each energy.
  // - order RefLines by this "significance", and limit number of lines by this, by number of lines and/or
  //   actual significance value.
  //
  // This entire function is incredible not well tested (as of 20250826), and very, very, heuristic based, that
  // isnt based on much.
  //
  // Future things to look at/consider
  // - maybe a simple filtering by the simple "importance", defined above, would be good enough, or even better
  // - The "activity" estimate is amazingly crude - we could potentually put in more computation and
  // - It would be great to kinda consider the near-by gamma intensities - like if you have a giant gamma line, you
  //   probably wont see the tiny one right next to it - but also, if a source is really shielded, you would care
  //   about the high energy tiny lines, and not the low-energy big lines, so the selection should adapt to only
  //   consider the intensities some-what around each line
  // - Could take into account actual detected peaks
  // - Could derive a better "importance" value, using the reference spectra over all detectors.
   
  
  if( !meas )
    return;
  
  // The number of lines we are keeping is totally arbitrary!
  size_t max_lines = 125;
  size_t min_we_care_about = 50; //The number of lines where if we have less than this amount, we wont filter them
  size_t num_lines_before_sig_filter = 45;
  double min_coarse_significance = 0.25; //TODO: check out using this - see below
  
  //Number un-filtered lines: {Ra226: 426}, {U232: 238}, {Th232: 510}, {Eu154: 188}, {U235: 1031}, {Pu239: 1289} {Am241: 1356}, {Np237: 1221}
  
  switch( src )
  {
    case RefLineSrc::EscapeLines:
    case RefLineSrc::RandomSumLines:
      return;
      
    case RefLineSrc::AlwaysShowing:
      max_lines = 200;
      min_we_care_about = 25;
      min_coarse_significance = 0.25;
      break;
      
    case RefLineSrc::ExternalRid:
      max_lines = 250;
      min_we_care_about = 50;
      min_coarse_significance = 0.0;
      break;
      
    case RefLineSrc::OnboardRid:
      max_lines = 150;
      min_we_care_about = 50;
      min_coarse_significance = 0.01;
      break;
      
    case RefLineSrc::UserPeakLines:
      max_lines = 275;
      min_we_care_about = 50;
      min_coarse_significance = 0.0;
      break;
      
    case RefLineSrc::UserPeakAssociatedNucLines:
      max_lines = 50;
      min_we_care_about = 25;
      min_coarse_significance = 0.2;
      break;
      
    case RefLineSrc::CharacteristicLine:
      max_lines = 75;
      min_we_care_about = 15;
      min_coarse_significance = 0.25;
      break;
  }//switch( src )
  
  if( ref_lines.m_ref_lines.size() < min_we_care_about )
    return;
  
  const bool have_fwhm = (detector && detector->hasResolutionInfo());
  
  
  double sum_importances = 0.0;
  // tuple<importance, data_counts, activity, RefLine>
  vector<tuple<double,double,double,ReferenceLineInfo::RefLine>> ref_lines_importances( ref_lines.m_ref_lines.size() );
  for( size_t index = 0; index < ref_lines.m_ref_lines.size(); ++index )
  {
    const ReferenceLineInfo::RefLine &ref = ref_lines.m_ref_lines[index];
    double importance = 0.0;
    if( (ref.m_energy > 10.0) && (ref.m_normalized_intensity > 0.0) )
      importance = ref.m_normalized_intensity * sqrt(ref.m_energy);
    
    const double peak_sigma = have_fwhm ? detector->peakResolutionSigma(ref.m_energy) : 2.5;
    const double data_counts = meas->gamma_integral( ref.m_energy - 2.0*peak_sigma, ref.m_energy + 2.0*peak_sigma );
    const double sf = data_counts / ref.m_normalized_intensity;

    ref_lines_importances[index] = {importance, data_counts, sf, ref};
    sum_importances += importance;
  }//for( size_t index = 0; index < ref_lines.m_ref_lines.size(); ++i )
  
  // Sort by importance
  std::sort( begin(ref_lines_importances), end(ref_lines_importances), []( const auto &lhs, const auto &rhs ){
    return get<0>(lhs) > get<0>(rhs);
  } );
  
  double max_activity = -1.0;
  for( size_t i = 0 ; (i < 6) && (i < ref_lines_importances.size()); ++i )
    max_activity = std::max( max_activity, get<2>(ref_lines_importances[i]) );

  vector<pair<double,ReferenceLineInfo::RefLine>> sig_ref_lines;
  sig_ref_lines.reserve( ref_lines_importances.size() );
  for( size_t i = 0; i < ref_lines_importances.size(); ++i )
  {
    const double data = get<1>(ref_lines_importances[i]);
    const double pred = max_activity*get<3>(ref_lines_importances[i]).m_normalized_intensity;
    const double nsig = (pred / sqrt( data > 1.0 ? data : 1.0 ));
    sig_ref_lines.emplace_back( nsig, get<3>(ref_lines_importances[i]) );
  }
  
  // Sort by significance
  std::sort( begin(sig_ref_lines), end(sig_ref_lines), []( const auto &lhs, const auto &rhs ){
    return get<0>(lhs) > get<0>(rhs);
  } );
  
  vector<ReferenceLineInfo::RefLine> sig_enough_ref_lines;
  sig_enough_ref_lines.reserve( std::min(sig_ref_lines.size(),max_lines) );
  
  for( size_t i = 0; (i < max_lines) && (i < sig_ref_lines.size()); ++i )
  {
    if( (i >= min_we_care_about)
       && (i >= num_lines_before_sig_filter)
       && (get<0>(sig_ref_lines[i]) < min_coarse_significance) )
    {
      break;
    }
    
    sig_enough_ref_lines.push_back( get<1>(sig_ref_lines[i]) );
  }
  
  static mutex s_cout_mutex;
  
  std::lock_guard<std::mutex> lock( s_cout_mutex );
  
  if( sig_enough_ref_lines.size() < sig_ref_lines.size() )
    cout << "First eliminiated line from "<< ref_lines.m_input.m_input_txt
    << " at index=" << sig_enough_ref_lines.size()
    << " is E=" << sig_ref_lines[sig_enough_ref_lines.size()].second.m_energy
    << " with intensity=" << sig_ref_lines[sig_enough_ref_lines.size()].second.m_normalized_intensity
    << " and significance=" << get<0>(sig_ref_lines[sig_enough_ref_lines.size()])
    << endl;
  else
    cout << "Didnt eliminate any lines for " << ref_lines.m_input.m_input_txt << endl;
  
  std::sort( begin(sig_enough_ref_lines), end(sig_enough_ref_lines), []( const ReferenceLineInfo::RefLine &lhs, const ReferenceLineInfo::RefLine &rhs ){
    return lhs.m_energy < rhs.m_energy;
  } );
  
  cout << "Reducing " << ref_lines.m_input.m_input_txt << " from " << ref_lines.m_ref_lines.size() << " lines to " << sig_enough_ref_lines.size() << endl;
  ref_lines.m_ref_lines = std::move(sig_enough_ref_lines);
}//void filterLines(...)


void RefLineDynamic::finishUpdateLines( const std::shared_ptr<std::vector<std::pair<double,ReferenceLineInfo>>> &ref_lines_ptr,
                       const std::shared_ptr<std::string> &js_fwhm_fcn_ptr,
                                       const size_t calc_num )
{
  m_current_ref_lines = ref_lines_ptr;
  
  assert( ref_lines_ptr && js_fwhm_fcn_ptr );
  if( !ref_lines_ptr || !js_fwhm_fcn_ptr )
    return;
  
  if( calc_num != m_current_calc_num->load() )
  {
    cerr << "Got a stale RefLineUpdate - not setting." << endl;
    return;
  }
  
  vector<pair<double,ReferenceLineInfo>> ref_lines = *ref_lines_ptr;
  std::string &js_fwhm_fcn = *js_fwhm_fcn_ptr;
  
  for( pair<double,ReferenceLineInfo> &ref_line : ref_lines )
    assignColorToInput( ref_line.second );
  
  cout << "About to load " << ref_lines.size() << " dynamic ref lines to client" << endl;
  
  // Note: since we are being called from within the render cycle most likely, this next function
  //       call will immediately send the appropriate JS to the client, as long as the chart has been rendered
  //       (i.e., it doesnt scheduleRender to lazily load to client, as long as its rendered).
  m_chart->setDynamicRefernceLines( std::move(ref_lines), std::move(js_fwhm_fcn) );
  
  wApp->triggerUpdate();
}//void finishUpdateLines(...)


void RefLineDynamic::startUpdateLines()
{
  if( !m_active )
  {
    m_chart->setDynamicRefernceLines( {}, "" );
    return;
  }//if( !m_active )
  
  shared_ptr<const SpecUtils::Measurement> foreground_orig = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  if( !foreground_orig )
  {
    m_chart->setDynamicRefernceLines( {}, "" );
    return;
  }//if( !m_active )
  
  
  (*m_current_calc_num) += 1;
  shared_ptr<atomic<size_t>> calc_num = m_current_calc_num;
  const size_t starting_calc_num = calc_num->load();

  const bool highres = PeakFitUtils::is_likely_high_res( m_interspec );
  const MaterialDB * const materialDb = m_interspec->materialDataBase(); //MaterialDB itself is thread-safe
  
  const set<int> &foreground_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
  const shared_ptr<const SpecMeas> foreground_meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  
  // We'll make a copy of the foreground so we dont run into any mutlithreading issues with the original being changed while we are
  //  accessing it in the below worker thread
  shared_ptr<const SpecUtils::Measurement> foreground = foreground_orig ? make_shared<const SpecUtils::Measurement>(*foreground_orig)
                                                                        : nullptr;
  
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
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> user_foreground_peaks_gui = pmodel ? pmodel->peaks() : nullptr;
  
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> user_foreground_peaks = user_foreground_peaks_gui
          ? make_shared<deque<shared_ptr<const PeakDef>>>( *user_foreground_peaks_gui ) : nullptr;
  
  const double background_lt = background ? background->live_time() : 0.0;
  const double foreground_lt = foreground ? foreground->live_time() : 0.0;
  
  vector<SpecUtils::DetectorAnalysisResult> det_ana;
  if( foreground_meas && foreground_meas->detectors_analysis() )
    det_ana = foreground_meas->detectors_analysis()->results_;
  
  shared_ptr<const DetectorPeakResponse> detector_ptr = foreground_meas ? foreground_meas->detector() : nullptr;
  // Make a copy of the detector to not run into any multi-threading issues
  if( detector_ptr )
    detector_ptr = make_shared<DetectorPeakResponse>( *detector_ptr );
  
  // We'll create a combination of user and auto-search peaks, that dont overlap
  const vector<shared_ptr<const PeakDef>> unique_foreground_peaks
                                  = combine_nonoverlapping_peaks( user_foreground_peaks, foreground_autosearch_peaks );
  const vector<shared_ptr<const PeakDef>> unique_background_peaks
                                  = combine_nonoverlapping_peaks( user_background_peaks, background_autosearch_peaks );
  
  IsotopeSearchByEnergy * const search = m_interspec->nuclideSearch();
  const vector<IsotopeSearchByEnergy::NucSearchCategory> categories = search ? search->search_categories()
                                                                  : vector<IsotopeSearchByEnergy::NucSearchCategory>{};
  const vector<ExternalRidIsotope> ext_rid_isotopes = m_external_rid_results ? m_external_rid_results->isotopes
                                                                  : vector<ExternalRidIsotope>{};
  const shared_ptr<const AlwaysSrcs> always_srcs = m_always_srcs;
  
  shared_ptr<vector<pair<double,ReferenceLineInfo>>> ref_lines_answer = make_shared<vector<pair<double,ReferenceLineInfo>>>();
  shared_ptr<string> js_fwhm_fcn = make_shared<string>();
  
  const string sessionId = wApp->sessionId();
  boost::function<void()> updaterfcn = wApp->bind( boost::bind( &RefLineDynamic::finishUpdateLines,
                                              this, ref_lines_answer, js_fwhm_fcn, starting_calc_num
  ) );
  
  
  const auto do_work = [foreground, det_ana, unique_foreground_peaks, unique_background_peaks,
                        foreground_lt, background_lt, highres, always_srcs, materialDb, user_foreground_peaks,
                        ext_rid_isotopes, categories, detector_ptr,
                        calc_num, starting_calc_num, sessionId,
                        js_fwhm_fcn, ref_lines_answer,updaterfcn](){
    
//#define DEBUG_DYNAMIC_REF_LINES_TIMING
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    const auto lambda_start_time = std::chrono::high_resolution_clock::now();
    const auto get_elapsed_milliseconds = [lambda_start_time]() -> long long {
      return std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - lambda_start_time).count();
    };
#endif
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      return;
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Initial setup and database access completed" << endl;
#endif
    
    //auto scoped_lock = foreground_meas? make_unique<lock_guard<std::recursive_mutex>>( foreground_meas->mutex() ) : nullptr;
    
    vector<pair<RefLineSrc,pair<double,ReferenceLineInfo>>> working_answer;
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
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
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Filtered non-background peaks" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    shared_ptr<const DetectorPeakResponse> detector = detector_ptr;
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
        assert( 0 );
        cerr << "findCandidates(...): error opening default detector file" << endl;
        detPtr->setIntrinsicEfficiencyFormula( "1.0", 3.0*PhysicalUnits::cm, PhysicalUnits::keV,
                                              0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic);
      }
      
      detector = detPtr;
    }//if( !detector || !detector->isValid() )
    
    assert( detector );
    if( detector->isValid() && !detector->hasResolutionInfo() )
    {
      // TODO: fit this from unique_foreground_peaks, maybe using DetectorPeakResponse::fitResolution(...), or MakeDrfFit::performResolutionFit(...);
      
      auto new_drf = make_shared<DetectorPeakResponse>( *detector );
      try
      {
        //Internally, `DetectorPeakResponse::fitResolution` calls MakeDrfFit::performResolutionFit
        shared_ptr<deque<shared_ptr<const PeakDef>>> peaks
                   = make_shared<deque<shared_ptr<const PeakDef>>>( begin(unique_foreground_peaks), end(unique_foreground_peaks) );
        new_drf->fitResolution( peaks, foreground, DetectorPeakResponse::kSqrtPolynomial ); //DetectorPeakResponse::kSqrtEnergyPlusInverse
      }catch( std::exception & )
      {
        const vector<float> coefs = highres ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
        new_drf->setFwhmCoefficients( coefs, DetectorPeakResponse::kGadrasResolutionFcn );
      }
      
      detector = new_drf;
    }//if( !detector->hasResolutionInfo() )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Loaded detector efficiency function" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    const ReactionGamma * const reaction_db = ([]()-> const ReactionGamma *{
      try{ return ReactionGammaServer::database(); }catch(...){}
      return static_cast<const ReactionGamma *>(nullptr); }
    )();
    
    char buffer[64] = { '\0' };
    
    SpecUtilsAsync::ThreadPool pool;
    
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
    
    //cout << "Before always_srcs working_answer={";
    //for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    //  cout << "'" << ref_lines.second.second.m_input.m_input_txt << "', ";
    //cout << "}" << endl;
    
    if( always_srcs )
    {
      // Convert individual sources to ReferenceLineInfo objects
      for( const auto &src : always_srcs->individual_sources )
      {
        RefLineInput input = base_input;
        input.m_input_txt = src.m_name;
        
        if( src.m_age.has_value() )
        {
          input.m_age = PhysicalUnits::printToBestTimeUnits( src.m_age.value(), 6 );
        }else if( const SandiaDecay::Nuclide * const nuc = db->nuclide(src.m_name) )
        {
          input.m_age = PhysicalUnits::printToBestTimeUnits( PeakDef::defaultDecayTime(nuc), 6 );
        }
        
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
              input.setShieldingAttFcn( materialDb );
            }catch( std::exception &e )
            {
              cerr << "Error setting shielding for curve: " << e.what() << endl;
            }
          }//if( we have a defined shielding )
        }//if( src.shielding_material.has_value() )
        
        
        bool already_exists = false;
        for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        {
          auto &existing = working_answer[i];
          
          already_exists = (existing.second.second.m_input.m_input_txt == src.m_name);
          if( already_exists )
          {
            if( src.m_weight > existing.second.first )
              existing.second.first = src.m_weight;
            if( static_cast<int>(RefLineSrc::AlwaysShowing) < static_cast<int>(existing.first) )
              existing.first = RefLineSrc::AlwaysShowing;
          }
        }//for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        
        if( already_exists )
          continue;
        
        shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
        if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
        {
          assert( !std::count_if(begin(working_answer), end(working_answer), [&input](const auto &v ){
            return v.second.second.m_input.m_input_txt == input.m_input_txt;
          }) );
          
          working_answer.emplace_back( RefLineSrc::AlwaysShowing,
                                      pair<double,ReferenceLineInfo>{src.m_weight, std::move(*ref_info)} );
        }else
        {
          cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for individual source: " << src.m_name << endl;
        }
      }//for( const auto &src : always_srcs->individual_sources )
      
      // Convert nuclide mixtures to ReferenceLineInfo objects
      for( const auto &mix_pair : always_srcs->nuc_mixes )
      {
        const ReferenceLinePredef::NucMix &mix = mix_pair.second;
        
        RefLineInput input = base_input;
        input.m_input_txt = mix.m_name;
        input.m_age = mix.m_default_age_str;
        
        bool already_exists = false;
        for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        {
          auto &existing = working_answer[i];
          
          already_exists = (existing.second.second.m_input.m_input_txt == mix.m_name);
          if( already_exists )
          {
            if( mix.m_weight > existing.second.first )
              existing.second.first = mix.m_weight;
            if( static_cast<int>(RefLineSrc::AlwaysShowing) < static_cast<int>(existing.first) )
              existing.first = RefLineSrc::AlwaysShowing;
          }
        }//for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        
        if( already_exists )
          continue;
        
        shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
        if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
        {
          assert( !std::count_if(begin(working_answer), end(working_answer), [&input](const auto &v ){
            return v.second.second.m_input.m_input_txt == input.m_input_txt;
          }) );
          
          working_answer.emplace_back( RefLineSrc::AlwaysShowing,
                                      pair<double,ReferenceLineInfo>{mix.m_weight, std::move(*ref_info)} );
        }else
        {
          cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for nuclide mixture: " << mix.m_name << endl;
        }
      }//for( const auto &mix_pair : always_srcs->nuc_mixes )
      
      // Convert custom source lines to ReferenceLineInfo objects
      for( const auto &custom_pair : always_srcs->custom_lines )
      {
        const ReferenceLinePredef::CustomSrcLines &custom = custom_pair.second;
        
        RefLineInput input = base_input;
        input.m_input_txt = custom.m_name;
        
        bool already_exists = false;
        for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        {
          auto &existing = working_answer[i];
          
          already_exists = (existing.second.second.m_input.m_input_txt == custom.m_name);
          if( already_exists )
          {
            if( custom.m_weight > existing.second.first )
              existing.second.first = custom.m_weight;
            if( static_cast<int>(RefLineSrc::AlwaysShowing) < static_cast<int>(existing.first) )
              existing.first = RefLineSrc::AlwaysShowing;
          }
        }//for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        
        if( already_exists )
          continue;
        
        shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
        if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
        {
          assert( !std::count_if(begin(working_answer), end(working_answer), [&input](const auto &v ){
            return v.second.second.m_input.m_input_txt == input.m_input_txt;
          }) );
          
          working_answer.emplace_back( RefLineSrc::AlwaysShowing,
                                      pair<double,ReferenceLineInfo>{custom.m_weight, std::move(*ref_info)} );
        }else
        {
          cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for custom source: " << custom.m_name << endl;
        }
      }//for( const auto &custom_pair : always_srcs->custom_lines )
    }//if( always_srcs )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Processed 'always showing' sources" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    //cout << "Before user_foreground_peaks working_answer={";
    //for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    //  cout << "'" << ref_lines.second.second.m_input.m_input_txt << "', ";
    //cout << "}" << endl;
    
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
        for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        {
          auto &existing = working_answer[i];
          
          already_exists = (existing.second.second.m_input.m_input_txt == canonical_name);
          if( already_exists )
          {
            if( new_weight > existing.second.first )
              existing.second.first = new_weight;
            if( static_cast<int>(RefLineSrc::UserPeakLines) < static_cast<int>(existing.first) )
              existing.first = RefLineSrc::UserPeakLines;
          }
        }//for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
        
        if( already_exists )
          continue;
        
        // Create reference lines using canonical name
        RefLineInput input = base_input;
        input.m_input_txt = canonical_name;
        
        shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
        if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
        {
          assert( !std::count_if(begin(working_answer), end(working_answer), [&input](const auto &v ){
            return v.second.second.m_input.m_input_txt == input.m_input_txt;
          }) );
          
          working_answer.emplace_back( RefLineSrc::UserPeakLines,
                                      pair<double,ReferenceLineInfo>{new_weight, std::move(*ref_info)} );
        }else
        {
          cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for user peak source: " << canonical_name << endl;
        }
      }
    }//if( user_foreground_peaks && user_foreground_peaks->size() )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Processed user foreground peaks" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    // Now check for external RID results, and add those lines that havent been added
    // We will assign these nuclides a weight of 5
    for( const ExternalRidIsotope &isotope : ext_rid_isotopes )
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
      for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
      {
        auto &existing = working_answer[i];
        already_exists = (existing.second.second.m_input.m_input_txt == canonical_name);
        if( already_exists )
        {
          if( new_weight > existing.second.first )
            existing.second.first = new_weight;
          if( static_cast<int>(RefLineSrc::ExternalRid) < static_cast<int>(existing.first) )
            existing.first = RefLineSrc::ExternalRid;
        }//
      }//for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
      
      if( already_exists )
        continue;
      
      // Create reference lines using canonical name
      RefLineInput input = base_input;
      input.m_input_txt = canonical_name;
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        assert( !std::count_if(begin(working_answer), end(working_answer), [&input](const auto &v ){
          return v.second.second.m_input.m_input_txt == input.m_input_txt;
        }) );
        
        working_answer.emplace_back( RefLineSrc::ExternalRid,
                                    pair<double,ReferenceLineInfo>{new_weight, std::move(*ref_info)} );
      }else
      {
        cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for external RID source: " << canonical_name << endl;
      }
    }//for( const ExternalRidIsotope &isotope : ext_rid_isotopes )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Processed external RID isotopes" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    
    //cout << "Before det ana working_answer={";
    //for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    //  cout << "'" << ref_lines.second.second.m_input.m_input_txt << "', ";
    //cout << "}" << endl;
    
    
    // Use the detectors on-board RID results
    for( const SpecUtils::DetectorAnalysisResult &result : det_ana )
    {
      if( result.nuclide_.empty() )
        continue;
      
      string canonical_name;
      bool found_source = false;
      
      // Try to get the nuclide from SandiaDecay database
      if( const SandiaDecay::Nuclide *nuc = db->nuclide(result.nuclide_) )
      {
        canonical_name = nuc->symbol;
        found_source = true;
      }else
      {
        // Try to get element from SandiaDecay database
        if( const SandiaDecay::Element *el = db->element(result.nuclide_) )
        {
          canonical_name = el->symbol;
          found_source = true;
        }else if( reaction_db )
        {
          // Try to get reaction from ReactionGamma database
          try
          {
            vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
            reaction_db->gammas( result.nuclide_, possible_rctns );

            if( !possible_rctns.empty() )
            {
              canonical_name = possible_rctns[0].reaction->name();
              found_source = true;
            }
          }catch( std::exception & )
          {
            //we get here for example if nuclide did have paranthesis
          }
        }//
      }//if( result.nuclide_ is nuclide ) / else
      
      if( !found_source )
      {
        // TODO: check for algorithm specific names, like HEU, SNM, neutrons, etc
        continue;
      }
      
      // Check if this source is already in ref_lines using canonical name, and update weight if higher
      bool already_exists = false;
      double new_weight = 3.0;
      for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
      {
        auto &existing = working_answer[i];
        already_exists = (existing.second.second.m_input.m_input_txt == canonical_name);
        
        if( already_exists )
        {
          if( new_weight > existing.second.first )
            existing.second.first = new_weight;
          if( static_cast<int>(RefLineSrc::OnboardRid) < static_cast<int>(existing.first) )
            existing.first = RefLineSrc::OnboardRid;
        }
      }//for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
      
      if( already_exists )
        continue;
      
      // Create reference lines using canonical name
      RefLineInput input = base_input;
      input.m_input_txt = canonical_name;
      
      shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
      if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
      {
        assert( !std::count_if(begin(working_answer), end(working_answer), [&input](const auto &v ){
          return v.second.second.m_input.m_input_txt == input.m_input_txt;
        }) );
        
        working_answer.emplace_back( RefLineSrc::OnboardRid,
                                    pair<double,ReferenceLineInfo>{new_weight, std::move(*ref_info)} );
      }else
      {
        cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for on-board RID source: " << canonical_name << endl;
      }
    }//for( const SpecUtils::DetectorAnalysisResult &result : det_ana )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Processed on-board RID results" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    //cout << "Before more info working_answer={";
    //for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    //  cout << "'" << ref_lines.second.second.m_input.m_input_txt << "', ";
    //cout << "}" << endl;
    
    // Add associated nuclides for existing reference lines
    shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> more_nuc_info = MoreNuclideInfo::MoreNucInfoDb::instance();
    if( more_nuc_info )
    {
      // Create a copy of current ref_lines to iterate over (to avoid modifying collection while iterating)
      vector<pair<double,string>> current_ref_lines;
      current_ref_lines.resize( working_answer.size() );
      for( const auto &ref_line : working_answer )
        current_ref_lines.emplace_back( ref_line.second.first, ref_line.second.second.m_input.m_input_txt );
      
      for( const pair<double,string> &ref_line_pair : current_ref_lines )
      {
        const double parent_weight = ref_line_pair.first;
        const string &ref_line_name = ref_line_pair.second;
        
        const MoreNuclideInfo::NucInfo *nuc_info = more_nuc_info->info(ref_line_name);
        if( !nuc_info || nuc_info->m_associated.empty() )
          continue;
        
        for( string associated_name : nuc_info->m_associated )
        {
          // Check if this associated nuclide is already in ref_lines
          bool already_exists = false;
          double associated_weight = 0.25 * parent_weight;
          
          // We'll normalize the source name for nuclides, even though it would be fine if we didnt, so this
          //  way we catch duplicates a little earlier
          const SandiaDecay::Nuclide * const nuc = db->nuclide(associated_name);
          if( nuc )
            associated_name = nuc->symbol;
          
          for( size_t i = 0; !already_exists && (i < working_answer.size()); i += 1 )
          {
            auto &existing = working_answer[i];
            already_exists = (existing.second.second.m_input.m_input_txt == associated_name);
            if( already_exists )
            {
              // Use max of current weight and 0.25 times parent weight
              existing.second.first = std::max(existing.second.first, associated_weight);
              if( static_cast<int>(RefLineSrc::UserPeakAssociatedNucLines) < static_cast<int>(existing.first) )
                existing.first = RefLineSrc::UserPeakAssociatedNucLines;
            }
          }//for( auto& existing : ref_lines )
          
          if( already_exists )
            continue;
          
          // Create reference lines for associated nuclide
          RefLineInput associated_input = base_input;
          associated_input.m_input_txt = associated_name;
          
          shared_ptr<ReferenceLineInfo> associated_ref_info = ReferenceLineInfo::generateRefLineInfo( associated_input );
          if( !associated_ref_info || (associated_ref_info->m_validity != ReferenceLineInfo::InputValidity::Valid) )
          {
            cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for associated nuclide: "
            << associated_name << endl;
            continue;
          }
            
          already_exists = false;
          for( size_t i = 0; !already_exists && (i < working_answer.size()); ++i )
          {
            // Important: we need to compare the input text from the `associated_ref_info`, and not `associated_name`,
            //            because the text may have been changed.
            already_exists = (working_answer[i].second.second.m_input.m_input_txt == associated_ref_info->m_input.m_input_txt);
            if( already_exists )
            {
              working_answer[i].second.first = std::max(working_answer[i].second.first, associated_weight);
              if( static_cast<int>(RefLineSrc::UserPeakAssociatedNucLines) < static_cast<int>(working_answer[i].first) )
                working_answer[i].first = RefLineSrc::UserPeakAssociatedNucLines;
            }
          }//for( auto& existing : ref_lines )
            
          if( already_exists )
            continue;
            
          assert( !std::count_if(begin(working_answer), end(working_answer), [&associated_input](const auto &v ){
            return v.second.second.m_input.m_input_txt == associated_input.m_input_txt;
          }) );
              
          working_answer.emplace_back( RefLineSrc::UserPeakAssociatedNucLines,
                                        pair<double,ReferenceLineInfo>{associated_weight, std::move(*associated_ref_info)} );
        }//for( const string &associated_name : nuc_info->m_associated )
      }//for( ref_line_pair in current_ref_lines )
    }//if( more_nuc_info )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Added associated nuclides" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
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
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
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
        snprintf( buffer, sizeof(buffer), "S.E. of %.1f keV", peak_energy );
        se_line.m_decaystr = buffer;
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
        snprintf( buffer, sizeof(buffer), "D.E. of %.1f keV", peak_energy );
        de_line.m_decaystr = buffer;
        de_line.m_parent_nuclide = peak->parentNuclide();
        de_line.m_transition = peak->nuclearTransition();
        de_line.m_reaction = peak->reaction();
        de_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape;
        escape_lines.push_back( de_line );
      }
    }//for( const auto& peak : unique_foreground_peaks )
    
    //cout << "Before escape peaks working_answer={";
    //for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    //  cout << "'" << ref_lines.second.second.m_input.m_input_txt << "', ";
    //cout << "}" << endl;
    
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
        
        working_answer.emplace_back( RefLineSrc::EscapeLines,
                                    pair<double,ReferenceLineInfo>{1.0, std::move(escape_info)} );
      }
    }
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Generated escape peak lines" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
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
      std::mutex candidate_mutex;
      set<string> candidate_source_names;
      map<string,vector<shared_ptr<const PeakDef>>> candidate_to_peaks;
      
      //Before paralizing, in debug, took 7 seconds for a specific test case; after niave paralization was 1.4 seconds.
      
      cout << "Finding characteristics for " << peaks_for_characteristics.size() << " remaining automated search peaks:" << endl;
      for( const shared_ptr<const PeakDef> &peak : peaks_for_characteristics )
      {
        pool.post( [peak,&nonbackground_peaks,&foreground,&detector,&candidate_to_peaks,&candidate_source_names,&candidate_mutex](){
          vector<string> characteristicnucs;
          
          // TODO: `IsotopeId::suggestNuclides(...)` only suggests nuclide - should also call `IsotopeId::findCharacteristics(...)` to get reactions and x-rays.
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
                
                std::lock_guard<std::mutex> lock( candidate_mutex );
                candidate_to_peaks[sugestions[i].nuclide->symbol].push_back(peak);
              }
              //cout << "suggestNuclides(" << peak->mean() << ", " << i << ") gave " << sugestions[i].nuclide->symbol << ", with weight " << sugestions[i].weight << endl;
            }
          }
          
          cout << "Peak at " << peak->mean() << " keV: ";
          //cout << "Found " << characteristicnucs.size() << " candidates: ";
          // Collect top 5 characteristics
          for( size_t i = 0; i < characteristicnucs.size() && i < 5; ++i )
          {
            //cout << characteristicnucs[i] << ", ";
            std::lock_guard<std::mutex> lock( candidate_mutex );
            candidate_source_names.insert( characteristicnucs[i] );
          }
        } );
      }//for( const shared_ptr<const PeakDef> &peak : peaks_for_characteristics )
      
      pool.join();
      
      if( calc_num->load() != starting_calc_num )
      {
        cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
        return;
      }
      
      cout << "Starting candidate_source_names={";
      for( auto name : candidate_source_names )
        cout << name << ", ";
      cout << "}" << endl;
      
      // Now filter the collected source names once
      set<string> characteristic_sources;
      
      for( const string &source_name : candidate_source_names )
      {
        // Check if already in ref_lines
        bool already_in_ref_lines = false;
        for( size_t index = 0; !already_in_ref_lines && (index < working_answer.size()); ++index )
          already_in_ref_lines = SpecUtils::iequals_ascii(working_answer[index].second.second.m_input.m_input_txt, source_name);
        
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
          for( size_t index = 0; !is_descendant && (index < working_answer.size()); ++index )
          {
            const SandiaDecay::Nuclide *existing_nuc = db->nuclide( working_answer[index].second.second.m_input.m_input_txt );
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
      
      if( calc_num->load() != starting_calc_num )
      {
        cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
        return;
      }
      
      //cout << "Before characteristics working_answer={";
      //for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
      //  cout << "'" << ref_lines.second.second.m_input.m_input_txt << "', ";
      //cout << "}" << endl;
      
      // Print the final set of characteristic sources
      cout << "Filtered characteristic sources for potential reference lines (" << characteristic_sources.size() << "): ";
      for( const string &source : characteristic_sources )
      {
        cout << source << ", ";
        
        pool.post( [&source, &candidate_to_peaks, &candidate_mutex, db, reaction_db,
                     &detector, &foreground, &unique_foreground_peaks, &base_input, &working_answer](){
          const double atomic_nums[]   = { 1.0, 26.0, 74.0 };
          const double areal_density[] = { 0.0*PhysicalUnits::g_per_cm2, 10.0*PhysicalUnits::g_per_cm2, 25.0*PhysicalUnits::g_per_cm2 };
          
          double prof_weight = -999.9;
          
          vector<shared_ptr<const PeakDef>> peaks_for_candidate;
          
          {
            std::lock_guard<std::mutex> lock( candidate_mutex );
            peaks_for_candidate = candidate_to_peaks[source];
          }
          
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
            const double weight = IsotopeId::profile_weight( detector, foreground, unique_foreground_peaks, {}, srcgammas,
                                                            energies, windows, atomic_nums[i], areal_density[i], WString() );
            prof_weight = std::max( prof_weight, weight );
          }
          
          // A min allowed profile weight of 0.25 is arbitrary, but if anything, we could move it higher, to like 0.5
          const double min_allowed_profile_weight = 0.25;
          if( prof_weight < min_allowed_profile_weight )
          {
            //cout << "Eliminating '" << source << "' from profile weight only " << prof_weight << endl;
            return;
          }
          
          //cout << "For source " << source << " the profile weight is " << prof_weight << endl;
          
          // Create reference lines for associated nuclide
          RefLineInput candidate_input = base_input;
          candidate_input.m_input_txt = source;
          const double candidate_weight = 0.1 + 2.4 * (2.0 * std::max(0.0, prof_weight - 0.5)); //0.1 to 2.5 weight
          
          shared_ptr<ReferenceLineInfo> candidate_ref_info = ReferenceLineInfo::generateRefLineInfo( candidate_input );
          if( candidate_ref_info && (candidate_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid) )
          {
            std::lock_guard<std::mutex> lock( candidate_mutex );
            
            //Check if we already have this line
            bool have_source = false;
            for( size_t i = 0; !have_source && (i < working_answer.size()); ++i )
            {
              auto &lines = working_answer[i];
              have_source = (lines.second.second.m_input.m_input_txt == candidate_input.m_input_txt);
              if( have_source )
              {
                lines.second.first = std::max( lines.second.first, candidate_weight);
                if( static_cast<int>(RefLineSrc::CharacteristicLine) < static_cast<int>(lines.first) )
                  lines.first = RefLineSrc::CharacteristicLine;
              }
            }//for (size_t i = 0; !have_source && (i < working_answer.size()); ++i )
            
            if( !have_source )
            {
              assert( !std::count_if(begin(working_answer), end(working_answer), [&candidate_input](const auto &v ){
                return v.second.second.m_input.m_input_txt == candidate_input.m_input_txt;
              }) );
              
              working_answer.emplace_back( RefLineSrc::CharacteristicLine,
                                          pair<double,ReferenceLineInfo>{candidate_weight, std::move(*candidate_ref_info)} );
            }
          }else
          {
            cerr << "RefLineDynamic::updateLines(): Failed to generate valid reference line info for candidate nuclide: " << source << endl;
          }
        } );
      }//for( const auto& source : characteristic_sources )
      
      pool.join();
      
      cout << endl;
    }//if( !peaks_for_characteristics.empty() )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Found characteristic nuclides" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    // If the deadtime is over 15% (low resolution) or 25% (HPGe), add in lines for random sum peaks of largest peaks
    const double max_live_time_frac = highres ? 0.75 : 0.85;
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
        vector<tuple<double,double,shared_ptr<const PeakDef>,shared_ptr<const PeakDef>>> sum_peaks; // energy, rel probability, peak_i, peak_j
        
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
              sum_peaks.emplace_back( sum_energy, rel_intensity, peak_i, peak_j );
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
            const shared_ptr<const PeakDef> &peak_i = std::get<2>( sum_peak );
            const shared_ptr<const PeakDef> &peak_j = std::get<3>( sum_peak );
            const double peak_i_energy = peak_i->mean();
            const double peak_j_energy = peak_j->mean();
            
            ReferenceLineInfo::RefLine ref_line;
            ref_line.m_energy = energy;
            ref_line.m_normalized_intensity = rel_intensity;
            ref_line.m_decay_intensity = rel_intensity;
            ref_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
            ref_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak;
            ref_line.m_attenuation_applies = false;
            snprintf( buffer, sizeof(buffer), "Random sum %.1f + %.1f keV", peak_i_energy, peak_j_energy );
            ref_line.m_decaystr = buffer;
            
            if( peak_i->hasSourceGammaAssigned() || peak_j->hasSourceGammaAssigned() )
            {
              const string name_i = peak_i->sourceName();
              const string name_j = peak_j->sourceName();
              if( name_i == name_j )
              {
                ref_line.m_decaystr += " (" + name_i + ")";
              }else
              {
                ref_line.m_decaystr += " (" + (name_i.empty() ? string("?") : name_i)
                                       + " + " + (name_j.empty() ? string("?") : name_j) + ")";
              }//
            }//if( peak_i->hasSourceGammaAssigned() || peak_j->hasSourceGammaAssigned() )
            
            sum_ref_info.m_ref_lines.push_back( ref_line );
          }//for( const auto &sum_peak : sum_peaks )
          
          const double randum_sum_lines_weight = 1.0; //totally arbitrary
          working_answer.emplace_back( RefLineSrc::RandomSumLines,
                                      pair<double,ReferenceLineInfo>{randum_sum_lines_weight, std::move(sum_ref_info)} );
        }//if( !sum_peaks.empty() )
      }//if( summing_candidates.size() >= 2 )
    }//if( (lt > 0.0) && (rt > 0.0) && (lt < max_live_time_frac*rt) )
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Generated random sum peaks" << endl;
#endif
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread." << endl;
      return;
    }
    
    cout << "Starting to filter ref-lines" << endl;
    ref_lines_answer->reserve( working_answer.size() );
    
    //cout << "working_answer={";
    //for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    //  cout << "'" << ref_lines.second.second.m_input.m_input_txt << "', ";
    //cout << "}" << endl;
    
    set<string> added_lines;
    for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    {
      string ref_line_name = ref_lines.second.second.m_input.m_input_txt;
      
      if( !added_lines.count(ref_line_name) && !ref_lines.second.second.m_ref_lines.empty() )
      {
        added_lines.insert(ref_line_name);
        pool.post( [&ref_lines,&foreground,detector](){
          filterLines( ref_lines.second.second, ref_lines.first, foreground, detector );
        } );
      }else
      {
        ref_lines.second.second.reset();
        const size_t num_copies = std::count_if(begin(working_answer), end(working_answer), [&ref_line_name](const auto &v ){
          return v.second.second.m_input.m_input_txt == ref_line_name;
        });
        cerr << "Duplicate ref lines '" << ref_line_name << "' - have " << num_copies << " copies." << endl;
      }
    }//for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    pool.join();
    cout << "Done to filter ref-lines" << endl;
    
    
    added_lines.clear();
    for( pair<RefLineSrc,pair<double,ReferenceLineInfo>> &ref_lines : working_answer )
    {
      string ref_line_name = ref_lines.second.second.m_input.m_input_txt;
      if( !ref_line_name.empty() && !added_lines.count(ref_line_name) )
      {
        added_lines.insert( std::move(ref_line_name) );
        assert( !ref_lines.second.second.m_ref_lines.empty() );
        ref_lines_answer->emplace_back( ref_lines.second.first, std::move(ref_lines.second.second) );
      }else
      {
        assert( ref_line_name.empty() );
        assert( ref_lines.second.second.m_ref_lines.empty() );
      }
    }
    
    // Get the FWHM JavaScript function from the foreground detector
    if( detector && detector->hasResolutionInfo() )
      try{ *js_fwhm_fcn = detector->javaScriptFwhmFunction(); }catch(...){ assert(0); } //dont expect it to ever throw
    
    if( js_fwhm_fcn->empty() )//We'll put in a generic fcnt
    {
      const vector<float> coefs = highres ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
      *js_fwhm_fcn = DetectorPeakResponse::javaScriptFwhmFunction( coefs, DetectorPeakResponse::kGadrasResolutionFcn );
    }//if( js_fwhm_fcn.empty() )
    
    
    if( calc_num->load() != starting_calc_num )
    {
      cerr << "Calc of DynamicRefLines is stale, stopping this thread (at last possible minute)." << endl;
      return;
    }
    
    cout << "Will post update of DynamicRefLines" << endl;
    WServer::instance()->post(sessionId, updaterfcn );
    
#ifdef DEBUG_DYNAMIC_REF_LINES_TIMING
    cerr << "DEBUG: " << get_elapsed_milliseconds() << " ms - Final processing and JS function setup completed" << endl;
#endif
  };//do_work_lambda

  WServer::instance()->ioService().boost::asio::io_service::post( std::move(do_work) );
}//void startUpdateLines()
