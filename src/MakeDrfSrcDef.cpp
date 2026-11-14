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
#include <regex>

#include <boost/algorithm/string/regex.hpp>

#include <Wt/WDate.h>
#include <Wt/WLabel.h>
#include <Wt/WTable.h>
#include <Wt/WComboBox.h>
#include <Wt/WDateEdit.h>
#include <Wt/WDateTime.h>
#include <Wt/WLineEdit.h>
#include <Wt/WCheckBox.h>
#include <Wt/WMenuItem.h>
#include <Wt/WPopupMenu.h>
#include <Wt/WTableCell.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WDoubleSpinBox.h>
#include <Wt/WDateValidator.h>
#include <Wt/WDoubleValidator.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WSuggestionPopup.h>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/UserLibrary.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
//#include "InterSpec/IsotopeNameFilterModel.h"

using namespace std;
using namespace Wt;

namespace
{
  const int sm_distance_row        = 1;
  const int sm_distance_uncert_row = 2;
  const int sm_activity_row        = 3;
  const int sm_activity_uncert_row = 4;
  const int sm_assay_date_row      = 5;
  const int sm_spec_date_row       = 6;
  const int sm_age_at_assay_row    = 7;
  const int sm_decayed_info_row    = 8;
  const int sm_shield_material_row = 9;
  const int sm_options_row         = 10;
  const int sm_problems_row        = 11;

  /** The format the date fields are entered in. */
  const char * const sm_date_format = "yyyy-MM-dd";
  
  
  // Returns if the candidate source is already in the vector of sources
  bool is_already_in( const vector<SrcLibLineInfo> &srcs, const SrcLibLineInfo &candidate )
  {
    // Comparing SrcLibLineInfo::m_source_name would probably be enough, but oh well.
    auto pos = std::find_if( begin(srcs), end(srcs),
                            [&candidate]( const SrcLibLineInfo &lhs ) -> bool {
      return( candidate.m_activity == lhs.m_activity
         && candidate.m_nuclide == lhs.m_nuclide
         && candidate.m_activity_date == lhs.m_activity_date
             && candidate.m_source_name == lhs.m_source_name );
    });
    
    return (pos != end(srcs));
  }//bool is_already_in( const std::vector<SrcLibLineInfo> &srcs, const SrcLibLineInfo &candidate )

  // Same as above, just for vector of pointers
  bool is_already_in( const vector<shared_ptr<const SrcLibLineInfo>> &srcs, const SrcLibLineInfo &candidate )
  {
    // Comparing SrcLibLineInfo::m_source_name would probably be enough, but oh well.
    auto pos = std::find_if( begin(srcs), end(srcs),
                            [&candidate]( const shared_ptr<const SrcLibLineInfo> &lhs ) -> bool {
      return( candidate.m_activity == lhs->m_activity
         && candidate.m_nuclide == lhs->m_nuclide
         && candidate.m_activity_date == lhs->m_activity_date
             && candidate.m_source_name == lhs->m_source_name );
    });
    
    return (pos != end(srcs));
  }//bool is_already_in( const std::vector<SrcLibLineInfo> &srcs, const SrcLibLineInfo &candidate )
  
}//namespace


SrcLibLineInfo::SrcLibLineInfo()
  : m_activity( 0.0 ),
    m_nuclide( nullptr ),
    m_activity_date{},
    m_source_name{},
    m_comments{},
    m_line{},
    m_activity_uncert( -1.0 ),
    m_distance( -1.0 )
{
}


vector<SrcLibLineInfo> SrcLibLineInfo::sources_in_lib( std::istream &file )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  vector<SrcLibLineInfo> answer;
  
  if( !file || !db )
    return answer;
  
  string line;
  while( SpecUtils::safe_get_line( file, line) )
  {
    SpecUtils::trim( line );
    
    vector<string> fields;
    // We will allow for a space, a tab, or the U+2002 (Unicode En Space) character to separate
    //  fields.  In order to do this properly, we need to use `split_regex`, and then manually
    //  remove empty fields
    boost::algorithm::split_regex( fields, line, boost::regex( "\\s|\xe2\x80\x82" ) );
    fields.erase( std::remove_if(begin(fields), end(fields), 
                    [](const string &v){return v.empty();}), end(fields) );
    
    if( fields.size() < 3 )
      continue;
    
    const string::size_type underscore_pos = fields[0].find("_");
    if( underscore_pos == string::npos )
      continue;
    const string nuc_str = fields[0].substr(0,underscore_pos);
    
    SrcLibLineInfo src_info;
    src_info.m_nuclide = db->nuclide( nuc_str );
    if( !src_info.m_nuclide )
      continue;
    
    src_info.m_source_name = fields[0];
    
    try
    {
      src_info.m_activity = stod(fields[1]) * PhysicalUnits::bq;
    }catch(...)
    {
      continue;
    }
    
    if( (src_info.m_activity <= 0.0) || IsNan(src_info.m_activity) || IsInf(src_info.m_activity) )
      continue;
    
    const SpecUtils::time_point_t datestr = SpecUtils::time_from_string( fields[2] );
    src_info.m_activity_date = to_ptime(datestr);
    if( src_info.m_activity_date.is_special() )
      continue;
    
    for( size_t i = 3; i < fields.size(); ++i )
      src_info.m_comments += ((i==3) ? "" : " ") + fields[i];
    
    
    // By default PhysicalUnits::sm_distanceRegex has a "^" character at beginning, and "$"
    //  character at end - lets get rid of these
    string dist_regex = PhysicalUnits::sm_distanceRegex;
    SpecUtils::ireplace_all( dist_regex, "^", "" );
    SpecUtils::ireplace_all( dist_regex, "$", "" );
    
    std::smatch dist_mtch;
    std::regex dist_expr( string(".+([dist|distance]\\s*\\=\\s*(") 
                         + dist_regex
                         + ")).*?", std::regex::icase );
    if( std::regex_match( src_info.m_comments, dist_mtch, dist_expr ) )
    {
      try
      {
        src_info.m_distance = PhysicalUnits::stringToDistance( dist_mtch[2].str() );
      }catch( std::exception & )
      {
        src_info.m_distance = -1.0;
      }
    }//if( std::regex_match( remark, dist_mtch, dist_expr ) )
    
    std::smatch act_uncert_mtch;
    const string act_uncert_expr_str = string(".*((Act|Activity)(Uncert|Uncertainty)\\s*[\\=:]\\s*(")
                                        + PhysicalUnits::sm_positiveDecimalRegex + ")).*?";
    
    std::regex act_uncert_expr( act_uncert_expr_str, std::regex::icase );
    if( std::regex_match( src_info.m_comments, act_uncert_mtch, act_uncert_expr ) )
    {
      string strval = act_uncert_mtch[4].str();
      if( !SpecUtils::parse_double( strval.c_str(), strval.size(), src_info.m_activity_uncert ) )
        src_info.m_activity_uncert = -1.0;
    }//if( std::regex_match( remark, dist_mtch, dist_expr ) )
    
    
    src_info.m_line = std::move(line);
    
    answer.push_back( std::move(src_info) );
  }//while( SpecUtils::safe_get_line( file, line) )
  
  return answer;
}//sources_in_lib(...)


vector<SrcLibLineInfo> SrcLibLineInfo::sources_in_lib( const string &filename )
{
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream file( wfilename.c_str(), ios::in | ios::binary );
#else
  ifstream file( filename.c_str(), ios::in | ios::binary );
#endif
  
  return sources_in_lib( file );
}//vector<SrcLibLineInfo> SrcLibLineInfo::sources_in_lib( const string &filename )


bool SrcLibLineInfo::is_candidate_src_lib( std::istream &file )
{
  size_t linenum = 0;
  string line;
  while( SpecUtils::safe_get_line( file, line) )
  {
    ++linenum;
    if( linenum > 5 )
      return false;
    
    SpecUtils::trim( line );
    if( (line.size() < 5) || (line.front() == '#') )
      continue;
    
    vector<string> fields;
    SpecUtils::split( fields, line, " \t" );
    if( fields.size() < 3 )
      return false;
    
    const auto underscore_pos = fields[0].find( "_" );
    if( (underscore_pos == string::npos) || (underscore_pos < 2) )
      return false;
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *nuc = db->nuclide( fields[0].substr(0,underscore_pos) );
    if( !nuc )
      return false;
    
    try
    {
      double activity = stod(fields[1]);
      if( (activity <= 0) || IsNan(activity) || IsInf(activity) )
        return false;
    }catch(...)
    {
      return false;
    }
    
    const SpecUtils::time_point_t datestr = SpecUtils::time_from_string( fields[2] );
    if( SpecUtils::is_special( datestr ) )
      return false;
    
    return true;
  }//while( SpecUtils::safe_get_line( file, line) )
  
  return false;
}//bool is_candidate_src_lib( std::istream &strm );


std::vector<SrcLibLineInfo> SrcLibLineInfo::sources_in_all_libs()
{
  vector<SrcLibLineInfo> source_lib_srcs;
  
  // Libraries earlier in the list win a repeated source name.  First the user library's: the user's
  //  (`source_libs/*.lib`, and "Source.lib" in the user data folder), then packages'.
  vector<string> source_lib_files;
  set<string> skip_files;  //handled by the user library (incl. those the user disabled)
  string packages_dir;
#if( USE_USER_LIBRARY )
  try
  {
    if( UserLibrary::available() )
    {
      packages_dir = UserLibrary::packagesDirectory();
      for( const UserLibrary::Item &item : UserLibrary::items( "source_lib", true ) )
      {
        skip_files.insert( SpecUtils::lexically_normalize_path( item.found.path ) );
        if( !item.disabled && !item.needsNewer )
          source_lib_files.push_back( item.found.path );
      }
    }//if( UserLibrary::available() )
  }catch( std::exception &e )
  {
    cerr << "Failed to list check-source libraries: " << e.what() << endl;
  }
#endif

  // Then "Source.lib" files anywhere else in InterSpec's data, the working, or the user data folder;
  //  not in packages (the user library knows which of those are installed and enabled).
  vector<string> base_paths{ InterSpec::staticDataDirectory(), SpecUtils::get_working_path() };
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_LOCAL_SERVER )
  try{ base_paths.push_back( InterSpec::writableDataDirectory() ); } catch(...){}
#endif
  for( const string &path : base_paths )
  {
    for( const string &file : SpecUtils::recursive_ls(path, "Source.lib" ) )
    {
      const bool in_packages = !packages_dir.empty() && !UserLibrary::relativePath( packages_dir, file ).empty();
      if( !in_packages && !skip_files.count( SpecUtils::lexically_normalize_path( file ) ) )
      {
        skip_files.insert( SpecUtils::lexically_normalize_path( file ) );
        source_lib_files.push_back( file );
      }
    }//for( loop over Source.lib files )
  }//for( loop over base paths )
  
  for( const string &lib : source_lib_files )
  {
    const vector<SrcLibLineInfo> these_sources = SrcLibLineInfo::sources_in_lib( lib );
    for( const auto &src : these_sources )
    {
      // Only add source if we havent already added it
      const auto pos = std::find_if( begin(source_lib_srcs), end(source_lib_srcs),
                                 [&src]( const SrcLibLineInfo &rhs ) -> bool {
        return src.m_source_name == rhs.m_source_name;
      } );
      if( pos == end(source_lib_srcs) )
        source_lib_srcs.push_back( src );
    }//for( const auto &src : these_sources )
  }//for( const string &lib : source_lib_files )
  
  return source_lib_srcs;
}//sources_in_all_libs()


MakeDrfSrcDef::MakeDrfSrcDef( const SandiaDecay::Nuclide *nuc,
              const boost::posix_time::ptime &measDate )
: WContainerWidget(),
  m_table( nullptr ),
  m_effGeomType( static_cast<int>(DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic) ),
  m_problemsTxt( nullptr ),
  m_nuclide( nuc ),
  m_nuclideLabel( nullptr ),
  m_distanceLabel( nullptr ),
  m_distanceEdit( nullptr ),
  m_distanceUncertLabel( nullptr ),
  m_distanceUncertEdit( nullptr ),
  m_activityLabel( nullptr ),
  m_activityEdit( nullptr ),
  m_activityUncertainty( nullptr ),
  m_useAgeInfo( nullptr ),
  m_assayDate( nullptr ),
  m_drfMeasurementDate( nullptr ),
  m_sourceInfoAtMeasurement( nullptr ),
  m_sourceAgeAtAssay( nullptr ),
  m_useShielding( nullptr ),
  m_shieldingSelect( nullptr ),
  m_lib_src_btn( nullptr ),
  m_lib_src_menu( nullptr ),
  m_lib_srcs_for_nuc{},
  m_lib_srcs_from_file{},
  m_lib_srcs_added{},
  m_updated()
{
  wApp->useStyleSheet( "InterSpec_resources/MakeDrfSrcDef.css" );
  
  addStyleClass( "MakeDrfSrcDef" );
  
  create();
  
  setNuclide( m_nuclide );
  
  if( !measDate.is_special() )
    m_drfMeasurementDate->setDate( WDateTime::fromTime_t( boost::posix_time::to_time_t(measDate) ).date() );

  validate();
}//MakeDrfSrcDef constructor


MakeDrfSrcDef::~MakeDrfSrcDef()
{
}


Wt::Signal<> &MakeDrfSrcDef::updated()
{
  return m_updated;
}


void MakeDrfSrcDef::setSrcInfo( const SrcLibLineInfo &info )
{
  if( info.m_distance > 0.0 )
    m_distanceEdit->setText( PhysicalUnits::printToBestLengthUnits(info.m_distance, 5) );
  
  if( info.m_activity_uncert > 0.0 )
    m_activityUncertainty->setValue( 100 * info.m_activity_uncert / info.m_activity );
  
  setAssayInfo( info.m_activity, info.m_activity_date );
}//void setSrcInfo( const SrcLibLineInfo &info )


void MakeDrfSrcDef::updateSourceLibNuclides()
{
  m_lib_srcs_for_nuc.clear();
  if( m_nuclide )
  {
    //TODO: we are currently re-parsing all the Source.lib files, every time this function
    //      gets called, for every source... we should probably be a little friendlier
    if( m_lib_srcs_from_file.empty() )
      m_lib_srcs_from_file = SrcLibLineInfo::sources_in_all_libs();
    
    for( const shared_ptr<const SrcLibLineInfo> &info : m_lib_srcs_added )
    {
      if( (info->m_nuclide == m_nuclide) && !is_already_in(m_lib_srcs_for_nuc,*info) )
        m_lib_srcs_for_nuc.push_back( *info );
    }
    
    for( const SrcLibLineInfo &info : m_lib_srcs_from_file )
    {
      if( (info.m_nuclide == m_nuclide)  && !is_already_in(m_lib_srcs_for_nuc,info) )
        m_lib_srcs_for_nuc.push_back( info );
    }
  }//if( m_nuclide )
  
  if( m_lib_src_btn )
    m_lib_src_btn->setHidden( m_lib_srcs_for_nuc.empty() );
  
  WPopupMenu *menu = m_lib_src_menu;
  if( menu )
  {
    const vector<WMenuItem *> old_items = menu->items();
    for( WMenuItem *item : old_items )
      menu->removeItem( item );
    
    for( const SrcLibLineInfo &src : m_lib_srcs_for_nuc )
    {
      assert( src.m_nuclide == m_nuclide );
      
      WMenuItem *item = menu->addItem( src.m_source_name );
      item->triggered().connect( this, [this, src](){ setSrcInfo( src ); } );
    }//for( const SrcLibLineInfo &src : m_lib_srcs_for_nuc )
  }//if( menu )
}//void updateSourceLibNuclides()


void MakeDrfSrcDef::setNuclide( const SandiaDecay::Nuclide *nuc )
{
  m_nuclide = nuc;
  
  updateSourceLibNuclides();
  
  const bool notMuchEvolution = (!nuc || PeakDef::ageFitNotAllowed(nuc));
  
  m_table->rowAt(sm_age_at_assay_row)->setHidden( notMuchEvolution );
  
  if( notMuchEvolution || !nuc )
  {
    m_sourceAgeAtAssay->setText( "0s" );
  }else
  {
    double ageval = 5.0*nuc->halfLife;
    if( nuc->canObtainPromptEquilibrium() )
      ageval = log(10000.0)/log(2.0) * nuc->promptEquilibriumHalfLife();
    if( ageval > 20*PhysicalUnits::year )
      ageval = 20*PhysicalUnits::year;
    m_sourceAgeAtAssay->setText( PhysicalUnitsLocalized::printToBestTimeUnits(ageval) );
  }
  
  if( nuc )
  {
    const string label = "<span class=\"SrcTitleNuc\">" + nuc->symbol + "</span>,"
                         " <span class=\"SrcTitleHl\">T&frac12;="
                         + PhysicalUnitsLocalized::printToBestTimeUnits(nuc->halfLife) + "</span>";
    m_nuclideLabel->setText( WString::fromUTF8(label) );
    m_useAgeInfo->show();
  }else
  {
    m_nuclideLabel->setText( WString::tr("msd-non-specified-nuc") );
    m_useAgeInfo->setUnChecked();
    m_distanceEdit->setValueText( "25 cm" );
    
    const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    m_activityEdit->setValueText( useCi ? "1 uCi" : "37 kBq" );
    m_useAgeInfo->hide();
  }
  
  useAgeInfoUserToggled();
}//setNuclide(...)


void MakeDrfSrcDef::create()
{
  m_table = addNew<WTable>();
  m_table->addStyleClass( "SrcInputTable" );

  WTableCell *cell = m_table->elementAt(0,0);
  cell->setColumnSpan( 2 );
  cell->addStyleClass( "SrcNuclideTitle" );
  m_nuclideLabel = cell->addNew<WText>();
  
/*
  //Code to put a nuclide suggestion into a WLineEdit so user could select nuclide.
  string replacerJs, matcherJs;
  PhotopeakDelegate::EditWidget::replacerJs( replacerJs );
  PhotopeakDelegate::EditWidget::nuclideNameMatcherJs( matcherJs );
  
  WSuggestionPopup *nucSuggest = new WSuggestionPopup( matcherJs, replacerJs );
 #if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  nucSuggest->setJavaScriptMember("wtNoReparent", "true");
#endif
  
  nucSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::Unit::FontEm) );
  nucSuggest->forEdit( m_nuclideEdit,
                       PopupTrigger::Editing | PopupTrigger::DropDownIcon );
  
  IsotopeNameFilterModel *filterModel = new IsotopeNameFilterModel();
  
  filterModel->excludeNuclides( false );
  filterModel->excludeXrays( true );
  filterModel->excludeEscapes( true );
  filterModel->excludeReactions( true );
  
  filterModel->filter( "" );
  nucSuggest->setFilterLength( -1 );
  nucSuggest->setModel( filterModel );
  nucSuggest->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );
  nucSuggest->filterModel().connect( filterModel, &IsotopeNameFilterModel::filter );
*/
  
  cell = m_table->elementAt(sm_distance_row,0);
  m_distanceLabel = cell->addNew<WLabel>( WString::tr("msd-distance") );
  cell = m_table->elementAt(sm_distance_row,1);
  m_distanceEdit = cell->addNew<WLineEdit>();
  m_distanceEdit->setTextSize( 16 );
  
  m_distanceEdit->setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  m_distanceEdit->setAttributeValue( "autocorrect", "off" );
  m_distanceEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_distanceLabel->setBuddy( m_distanceEdit );
  auto distValidator = std::make_shared<WRegExpValidator>( PhysicalUnits::sm_distanceUnitOptionalRegex );
  distValidator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  distValidator->setInvalidNoMatchText( WString::tr("msd-tt-dist-format") );
  m_distanceEdit->setValidator( distValidator );
  m_distanceEdit->setText( "50 cm" );
  m_distanceEdit->changed().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_distanceEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserInput );

  // Optional distance uncertainty: fully correlated for all this source's peaks, it propagates
  //  into the efficiency uncertainty through the geometric factor's distance slope.
  cell = m_table->elementAt(sm_distance_uncert_row,0);
  m_distanceUncertLabel = cell->addNew<WLabel>( WString::tr("msd-distance-uncert") );
  cell = m_table->elementAt(sm_distance_uncert_row,1);
  m_distanceUncertEdit = cell->addNew<WLineEdit>();
  m_distanceUncertEdit->setTextSize( 16 );
  m_distanceUncertEdit->setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  m_distanceUncertEdit->setAttributeValue( "autocorrect", "off" );
  m_distanceUncertEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_distanceUncertLabel->setBuddy( m_distanceUncertEdit );
  m_distanceUncertEdit->setValidator( distValidator );
  m_distanceUncertEdit->setPlaceholderText( WString::tr("msd-distance-uncert-ph") );
  m_distanceUncertEdit->changed().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_distanceUncertEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserInput );
  
  
  cell = m_table->elementAt(sm_activity_row,0);
  m_activityLabel = cell->addNew<WLabel>( WString::tr("msd-activity") );

  cell = m_table->elementAt(sm_activity_row,1);
  m_activityEdit = cell->addNew<WLineEdit>();
  
  m_activityEdit->setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  m_activityEdit->setAttributeValue( "autocorrect", "off" );
  m_activityEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_activityEdit->setTextSize( 16 );
  m_activityLabel->setBuddy( m_activityEdit );
  
  auto val = std::make_shared<WRegExpValidator>( PhysicalUnits::sm_activityRegex );
  val->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  val->setInvalidNoMatchText( WString::tr("msd-tt-act-format") );
  m_activityEdit->setValidator( val );
  const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  m_activityEdit->setText( useCi ? "100 uCi" : "3.7 MBq" );
  m_activityEdit->changed().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_activityEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserInput );

  cell = m_table->elementAt(sm_activity_uncert_row,0);
  WLabel *label = cell->addNew<WLabel>( WString::tr("msd-act-uncert") );
  cell = m_table->elementAt(sm_activity_uncert_row,1);
  m_activityUncertainty = cell->addNew<WDoubleSpinBox>();
  m_activityUncertainty->setValue( 0.0 );
  m_activityUncertainty->setTextSize( 14 );
  m_activityUncertainty->setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  m_activityUncertainty->setAttributeValue( "autocorrect", "off" );
  m_activityUncertainty->setAttributeValue( "spellcheck", "off" );
#endif
  m_activityUncertainty->setSuffix( " %" );
  label->setBuddy( m_activityUncertainty );
  
  m_activityUncertainty->setDecimals( 1 );
  m_activityUncertainty->setRange( 0.0, 100.0 );
  m_activityUncertainty->setSingleStep( 1.0 );
  
  m_activityUncertainty->changed().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_activityUncertainty->enterPressed().connect( this, &MakeDrfSrcDef::handleUserInput );
  
  //WDoubleValidator *percentVal = new WDoubleValidator( this );
  //percentVal->setRange( 0.0, 100.0 );
  //m_activityUncertainty->setValidator( percentVal );
  
  cell = m_table->elementAt(sm_assay_date_row,0);
  label = cell->addNew<WLabel>( WString::tr("msd-assay-date") );
  cell = m_table->elementAt(sm_assay_date_row,1);
  m_assayDate = cell->addNew<WDateEdit>();
  WidgetUtils::setupDateEdit( m_assayDate, sm_date_format );
  label->setBuddy( m_assayDate );
  m_assayDate->changed().connect( this, [this](){ handleDateInput( m_assayDate ); } );
  
  cell = m_table->elementAt(sm_spec_date_row,0);
  label = cell->addNew<WLabel>( WString::tr("msd-spec-date") );
  cell = m_table->elementAt(sm_spec_date_row,1);
  m_drfMeasurementDate = cell->addNew<WDateEdit>();
  WidgetUtils::setupDateEdit( m_drfMeasurementDate, sm_date_format );
  label->setBuddy( m_drfMeasurementDate );
  m_drfMeasurementDate->changed().connect( this, [this](){ handleDateInput( m_drfMeasurementDate ); } );
  
  
  cell = m_table->elementAt(sm_age_at_assay_row,0);
  label = cell->addNew<WLabel>( WString::tr("msd-age-at-assay") );
  cell = m_table->elementAt(sm_age_at_assay_row,1);
  m_sourceAgeAtAssay = cell->addNew<WLineEdit>();
  
  m_sourceAgeAtAssay->setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  m_sourceAgeAtAssay->setAttributeValue( "autocorrect", "off" );
  m_sourceAgeAtAssay->setAttributeValue( "spellcheck", "off" );
#endif
  label->setBuddy( m_sourceAgeAtAssay );
  m_sourceAgeAtAssay->setPlaceholderText( WString::tr("msd-age-at-assay-ph") );
  m_sourceAgeAtAssay->changed().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_sourceAgeAtAssay->enterPressed().connect( this, &MakeDrfSrcDef::handleUserInput );

  val = std::make_shared<WRegExpValidator>( PhysicalUnitsLocalized::timeDurationRegex() );
  val->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  val->setInvalidNoMatchText( WString::tr("msd-tt-age-format") );
  m_sourceAgeAtAssay->setValidator( val );
  m_sourceAgeAtAssay->setText( "0s" );
  
  
  cell = m_table->elementAt(sm_decayed_info_row,0);
  label = cell->addNew<WLabel>( WString::tr("msd-aging-res") );
  cell = m_table->elementAt(sm_decayed_info_row,1);
  m_sourceInfoAtMeasurement = cell->addNew<WText>();
  m_sourceInfoAtMeasurement->addStyleClass( "SelectableTxt" );
  
  //Wt::WDateEdit *m_sourceCreationDate;
  
  cell = m_table->elementAt(sm_options_row,0);

  m_lib_src_btn = cell->addNew<WPushButton>();
  m_lib_src_btn->setIcon( "InterSpec_resources/images/db_small.png" );
  m_lib_src_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  m_lib_src_btn->setFloatSide( Wt::Side::Left );
  
  m_lib_src_menu = makePopupMenu( m_lib_src_btn );
  m_lib_src_menu->setMaximumSize( WLength::Auto, WLength(15, WLength::Unit::FontEm) );
  
  m_lib_src_menu->setAutoHide( true, 2500 );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  HelpSystem::attachToolTipOn( m_lib_src_btn, WString::tr("msd-tt-lib-src"), showToolTips );

  m_useAgeInfo = cell->addNew<WCheckBox>( WString::tr("msd-use-age-cb") );
  m_useAgeInfo->setFloatSide( Wt::Side::Right );
  m_useAgeInfo->setChecked( false );
  m_useAgeInfo->checked().connect( this, &MakeDrfSrcDef::useAgeInfoUserToggled );
  m_useAgeInfo->unChecked().connect( this, &MakeDrfSrcDef::useAgeInfoUserToggled );
  
  cell = m_table->elementAt(sm_options_row,1);
  m_useShielding = cell->addNew<WCheckBox>( WString::tr("msd-shielded-cb") );
  m_useShielding->setFloatSide( Wt::Side::Right );
  m_useShielding->setChecked( false );
  m_useShielding->checked().connect( this, &MakeDrfSrcDef::useShieldingInfoUserToggled );
  m_useShielding->unChecked().connect( this, &MakeDrfSrcDef::useShieldingInfoUserToggled );
  
  
  cell = m_table->elementAt(sm_shield_material_row,0);
  cell->setColumnSpan( 2 );
  
  m_shieldingSelect = cell->addNew<ShieldingSelect>();
  m_shieldingSelect->materialModified().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_shieldingSelect->materialChanged().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_shieldingSelect->addingIsotopeAsSource().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_shieldingSelect->removingIsotopeAsSource().connect( this, &MakeDrfSrcDef::handleUserInput );
  m_shieldingSelect->hide();

  // What is wrong with the inputs, if anything - see validate()
  cell = m_table->elementAt(sm_problems_row,0);
  cell->setColumnSpan( 2 );
  m_problemsTxt = cell->addNew<WText>( "", Wt::TextFormat::Plain );
  m_problemsTxt->addStyleClass( "FieldErrorMsg SrcProblems" );
  m_problemsTxt->setInline( false );
  m_table->rowAt(sm_problems_row)->hide();

  useAgeInfoUserToggled();
}//void create()


void MakeDrfSrcDef::useAgeInfoUserToggled()
{
  const bool useAge = m_useAgeInfo->isChecked();

  m_table->rowAt(sm_assay_date_row)->setHidden( !useAge );
  m_table->rowAt(sm_spec_date_row)->setHidden( !useAge );

  const bool notMuchEvolution = (!m_nuclide || PeakDef::ageFitNotAllowed(m_nuclide));
  m_table->rowAt(sm_age_at_assay_row)->setHidden( !useAge || notMuchEvolution );

  m_table->rowAt(sm_decayed_info_row)->setHidden( !useAge );

  m_assayDate->setDisabled( !useAge );
  m_drfMeasurementDate->setDisabled( !useAge );
  m_sourceInfoAtMeasurement->setDisabled( !useAge );
  m_sourceAgeAtAssay->setDisabled( !useAge );

  // While aging is used both dates are needed, so Wt marks an empty one just like bad text.
  m_assayDate->dateValidator()->setMandatory( useAge );
  m_drfMeasurementDate->dateValidator()->setMandatory( useAge );

  updateActivityLabel();
  validate();
  updateAgedText();

  m_updated.emit();
}//void useAgeInfoUserToggled()


void MakeDrfSrcDef::updateActivityLabel()
{
  const char *key = "msd-activity";
  switch( static_cast<DetectorPeakResponse::EffGeometryType>(m_effGeomType) )
  {
    case DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic:
    case DetectorPeakResponse::EffGeometryType::FarFieldAbsolute:
      key = (m_useAgeInfo && m_useAgeInfo->isChecked()) ? "msd-assay-act" : "msd-activity";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:   key = "msd-total-act";    break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:  key = "msd-act-per-cm2";  break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:   key = "msd-act-per-m2";   break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram: key = "msd-act-per-gram"; break;
  }//switch( effective geometry type )

  m_activityLabel->setText( WString::tr(key) );
}//void updateActivityLabel()


void MakeDrfSrcDef::updateAgedText()
{
  try
  {
    if( !m_nuclide )
      throw runtime_error( "" );
    
    double activity = activityAtSpectrumTime();
    double age = ageAtSpectrumTime();
    const string agestr = PhysicalUnitsLocalized::printToBestTimeUnits(age);
    
    const string enteredAct = m_activityEdit->text().toUTF8();
    
    const bool useCurie = (enteredAct.find_first_of( "cC" ) != string::npos);
    const string actstr = PhysicalUnits::printToBestActivityUnits(activity, 1, useCurie );
    
    string txt = actstr;
    if( !PeakDef::ageFitNotAllowed(m_nuclide) )
      txt += ", " + agestr;

    m_sourceInfoAtMeasurement->setText( WString::fromUTF8(txt) );
  }catch( std::exception & )
  {
    // The inputs it is computed from have a problem, which validate() shows below
    m_sourceInfoAtMeasurement->setText( WString::fromUTF8("\xe2\x80\x94") );  //em dash
  }
}//void updateAgedText()


void MakeDrfSrcDef::handleUserInput()
{
  validate();
  updateAgedText();
  m_updated.emit();
}//void handleUserInput()


void MakeDrfSrcDef::handleDateInput( Wt::WDateEdit *edit )
{
  // "2019-1-5" and the like become "2019-01-05"; text that isnt a date is left for validate() to flag
  WidgetUtils::normalizeDateEdit( edit );
  handleUserInput();
}//void handleDateInput( Wt::WDateEdit *edit )


std::vector<MakeDrfSrcDef::InputProblem> MakeDrfSrcDef::validate()
{
  vector<InputProblem> problems;
  auto add = [&problems]( WFormWidget *field, const WString &msg ){
    InputProblem p;
    p.field = field;
    p.message = msg;
    problems.push_back( std::move(p) );
  };//add

  auto trimmed = []( const WString &txt ) -> string {
    return SpecUtils::trim_copy( txt.toUTF8() );
  };

  // Distance, and its optional uncertainty - not shown for fixed geometries
  if( !m_distanceEdit->isHidden() )
  {
    double dist = -1.0;
    const string txt = trimmed( m_distanceEdit->text() );
    if( txt.empty() )
    {
      add( m_distanceEdit, WString::tr("msd-err-dist-missing") );
    }else
    {
      try
      {
        dist = PhysicalUnits::stringToDistance( txt );
        if( !(dist > 0.0) || IsInf(dist) )
          add( m_distanceEdit, WString::tr("msd-err-dist-nonpos") );
      }catch( std::exception & )
      {
        add( m_distanceEdit, WString::tr("msd-err-dist-invalid") );
      }
    }//if( txt.empty() ) / else

    const string uncert_txt = trimmed( m_distanceUncertEdit->text() );
    if( !uncert_txt.empty() )
    {
      try
      {
        const double uncert = PhysicalUnits::stringToDistance( uncert_txt );
        if( (uncert < 0.0) || IsNan(uncert) || IsInf(uncert) )
          add( m_distanceUncertEdit, WString::tr("msd-err-dist-uncert-invalid") );
        else if( (dist > 0.0) && (uncert >= dist) )
          add( m_distanceUncertEdit, WString::tr("msd-err-dist-uncert-range") );
      }catch( std::exception & )
      {
        add( m_distanceUncertEdit, WString::tr("msd-err-dist-uncert-invalid") );
      }
    }//if( !uncert_txt.empty() )
  }//if( !m_distanceEdit->isHidden() )

  // Activity: an activity of zero would make the efficiency infinite
  {
    const string txt = trimmed( m_activityEdit->text() );
    if( txt.empty() )
    {
      add( m_activityEdit, WString::tr("msd-err-act-missing") );
    }else
    {
      try
      {
        const double act = PhysicalUnits::stringToActivity( txt );
        if( !(act > 0.0) || IsInf(act) )
          add( m_activityEdit, WString::tr("msd-err-act-nonpos") );
      }catch( std::exception & )
      {
        add( m_activityEdit, WString::tr("msd-err-act-invalid") );
      }
    }//if( txt.empty() ) / else
  }

  try
  {
    fractionalActivityUncertainty();
  }catch( std::exception & )
  {
    add( m_activityUncertainty, WString::tr("msd-err-act-uncert") );
  }

  // The assay and spectrum dates, and age at assay, when the activity is decayed
  if( m_nuclide && m_useAgeInfo->isChecked() )
  {
    auto check_date = [&]( WDateEdit *edit, const char *missing_key, const char *label_key ) -> WDate {
      const string txt = trimmed( edit->text() );
      WDate date;
      switch( WidgetUtils::parseLenientDate( txt, sm_date_format, date ) )
      {
        case WidgetUtils::DateParse::Empty:
          add( edit, WString::tr(missing_key) );
          return WDate();

        case WidgetUtils::DateParse::Ambiguous:
          add( edit, WString::tr("msd-err-date-ambiguous").arg( WString::tr(label_key) ).arg( txt ) );
          return WDate();

        case WidgetUtils::DateParse::Invalid:
          add( edit, WString::tr("msd-err-date-invalid").arg( WString::tr(label_key) ) );
          return WDate();

        case WidgetUtils::DateParse::Valid:
          break;
      }//switch( parse the date )

      // Readable, but not (yet) re-written into the fields format - which is what the getters read
      if( !edit->date().isValid() )
      {
        add( edit, WString::tr("msd-err-date-invalid").arg( WString::tr(label_key) ) );
        return WDate();
      }

      return date;
    };//check_date

    const WDate assay = check_date( m_assayDate, "msd-err-assay-date-missing", "msd-assay-date" );
    const WDate meas = check_date( m_drfMeasurementDate, "msd-err-spec-date-missing", "msd-spec-date" );
    if( assay.isValid() && meas.isValid() && (assay > meas) )
      add( m_assayDate, WString::tr("msd-err-assay-after-spec")
                          .arg( assay.toString(sm_date_format) )
                          .arg( meas.toString(sm_date_format) ) );

    if( !m_table->rowAt(sm_age_at_assay_row)->isHidden() )
    {
      // Blank, or just zeros, means no age - see ageAtSpectrumTime()
      const string agestr = trimmed( m_sourceAgeAtAssay->text() );
      if( !agestr.empty() && (agestr.find_first_not_of("+-0.") != string::npos) )
      {
        try
        {
          const double age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, m_nuclide->halfLife );
          if( (age < 0.0) || IsNan(age) || IsInf(age) )
            add( m_sourceAgeAtAssay, WString::tr("msd-err-age-neg") );
        }catch( std::exception & )
        {
          add( m_sourceAgeAtAssay, WString::tr("msd-err-age-invalid") );
        }
      }//if( an age was entered )
    }//if( age at assay is shown )
  }//if( m_nuclide && m_useAgeInfo->isChecked() )

  // Shielding between source and detector
  if( m_useShielding->isChecked() && m_shieldingSelect )
  {
    if( m_shieldingSelect->isGenericMaterial() )
    {
      try
      {
        m_shieldingSelect->atomicNumber();
      }catch( std::exception & )
      {
        add( m_shieldingSelect->atomicNumberEdit(), WString::tr("msd-err-shield-an") );
      }

      try
      {
        m_shieldingSelect->arealDensity();
      }catch( std::exception & )
      {
        add( m_shieldingSelect->arealDensityEdit(), WString::tr("msd-err-shield-ad") );
      }
    }else if( !m_shieldingSelect->material() )
    {
      add( m_shieldingSelect->materialEdit(), WString::tr("msd-err-shield-material") );
    }else
    {
      try
      {
        m_shieldingSelect->thickness();
      }catch( std::exception & )
      {
        add( m_shieldingSelect->thicknessEdit(), WString::tr("msd-err-shield-thick") );
      }
    }//if( generic material ) / else
  }//if( shielding is used )

  // Mark the fields.  Text a fields own validator rejects is marked by Wt itself ("Wt-invalid", as
  //  the user types); a field whose text is well-formed, but whose value is not usable, gets
  //  "InputValueError".  The activity uncertainty spin-box is left to Wt: its validator does not
  //  give the right answer in the request that created it (see fractionalActivityUncertainty()).
  vector<WFormWidget *> fields{ m_distanceEdit, m_distanceUncertEdit, m_activityEdit,
                                m_assayDate, m_drfMeasurementDate, m_sourceAgeAtAssay };
  if( m_shieldingSelect )
  {
    fields.push_back( m_shieldingSelect->materialEdit() );
    fields.push_back( m_shieldingSelect->thicknessEdit() );
    fields.push_back( m_shieldingSelect->atomicNumberEdit() );
    fields.push_back( m_shieldingSelect->arealDensityEdit() );
  }//if( m_shieldingSelect )

  for( WFormWidget *field : fields )
  {
    if( !field )
      continue;

    const bool has_problem = std::any_of( begin(problems), end(problems),
                                [field]( const InputProblem &p ){ return p.field == field; } );
    const bool value_problem = has_problem && (field->validate() == ValidationState::Valid);
    WidgetUtils::setInputValueError( field, value_problem );
  }//for( WFormWidget *field : fields )

  string txt;
  for( const InputProblem &p : problems )
    txt += (txt.empty() ? "" : "\n") + p.message.toUTF8();

  m_problemsTxt->setText( WString::fromUTF8(txt) );
  m_table->rowAt(sm_problems_row)->setHidden( txt.empty() );
  toggleStyleClass( "SrcDefHasError", !problems.empty() );

  return problems;
}//std::vector<InputProblem> validate()


void MakeDrfSrcDef::useShieldingInfoUserToggled()
{
  m_shieldingSelect->setHidden( !m_useShielding->isChecked() );
  handleUserInput();
}//void useShieldingInfoUserToggled()


double MakeDrfSrcDef::enteredActivity() const
{
  string activitystr = m_activityEdit->text().toUTF8();
  
  SpecUtils::trim( activitystr );
  
  return PhysicalUnits::stringToActivity( activitystr );
}//double enteredActivity()


double MakeDrfSrcDef::distance() const
{
  assert( !m_distanceEdit->isHidden() );
  
  if( m_distanceEdit->isHidden() )
    throw runtime_error( "MakeDrfSrcDef: should not be calling distance, when a fixed geometry." );
  
  const string txt = m_distanceEdit->text().toUTF8();
  return PhysicalUnits::stringToDistance(txt);
}

const SandiaDecay::Nuclide *MakeDrfSrcDef::nuclide() const
{
  return m_nuclide;
}


double MakeDrfSrcDef::activityAtSpectrumTime() const
{
  const double userActivity = enteredActivity();
  if( !m_useAgeInfo->isChecked() || !m_nuclide )
    return userActivity;
  
  const WDate measDate = m_drfMeasurementDate->date();
  const WDate assayDate = m_assayDate->date();
  
  if( !measDate.isValid() )
    throw runtime_error( "Measurement date invalid" );
  
  if( !assayDate.isValid() )
    throw runtime_error( "Assay date invalid" );
  
  if( assayDate > measDate )
    throw runtime_error( "Assay date must be before measurement date" );
  
  const int numDays = assayDate.daysTo( measDate );
  
  SandiaDecay::NuclideMixture mix;
  mix.addNuclideByActivity( m_nuclide, userActivity );
  
  return mix.activity( 24.0*3600.0*numDays, m_nuclide );
}//double activityAtSpectrumTime() const


double MakeDrfSrcDef::fractionalActivityUncertainty() const
{
  switch( m_activityUncertainty->validate() )
  {
    case WValidator::State::Invalid:
    {
      //We actually get here if the value is asked for if we have the same
      //  main Wt threadlock as when m_activityUncertainty is created; So we
      //  will, check if the uncertainty string is about what we initially set
      //  it to ("0.0 %"), and if so, just return zero.
      const string value = m_activityUncertainty->valueText().toUTF8();
      if( SpecUtils::istarts_with( value, "0.0 ") )
        return 0.0;
      
      throw runtime_error( "Activity Uncertainty Invalid" );
    }
      
    case WValidator::State::InvalidEmpty:
      return 0.0;
      
    case WValidator::State::Valid:
      break;
  }//switch( m_activityUncertainty->validate() )
  
  return m_activityUncertainty->value() / 100.0;
}//double fractionalActivityUncertainty() const


double MakeDrfSrcDef::ageAtSpectrumTime() const
{
  if( !m_nuclide )
    return 0.0;

  if( !m_useAgeInfo->isChecked() || PeakDef::ageFitNotAllowed(m_nuclide) )
    return PeakDef::defaultDecayTime( m_nuclide, nullptr );
  
  string ageAtAssaystr = m_sourceAgeAtAssay->text().toUTF8();
  SpecUtils::trim( ageAtAssaystr );
  
  double ageAtAssay = 0.0;
  if( !ageAtAssaystr.empty() && (ageAtAssaystr.find_first_not_of("+-0.")!=string::npos) )
    ageAtAssay = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( ageAtAssaystr, m_nuclide->halfLife );
  
  const WDate measDate = m_drfMeasurementDate->date();
  const WDate assayDate = m_assayDate->date();
  
  if( !measDate.isValid() )
    throw runtime_error( "Measurement date invalid" );
  
  if( !assayDate.isValid() )
    throw runtime_error( "Assay date invalid" );
  
  if( assayDate > measDate )
    throw runtime_error( "Assay date must be before measurement date" );
  
  const double numDays = assayDate.daysTo( measDate );
  
  return ageAtAssay + 24.0*3600.0*numDays;
}//double ageAtSpectrumTime() const


ShieldingSelect *MakeDrfSrcDef::shielding()
{
  if( !m_useShielding->isChecked() )
    return nullptr;
  return m_shieldingSelect;
}//ShieldingSelect *shielding()


void MakeDrfSrcDef::setDistance( const double dist )
{
  m_distanceEdit->setText( PhysicalUnits::printToBestLengthUnits(dist) );
  validate();
}//void setDistance( const double dist );


double MakeDrfSrcDef::distanceUncertainty() const
{
  string txt = m_distanceUncertEdit->text().toUTF8();
  SpecUtils::trim( txt );
  if( txt.empty() )
    return 0.0;
  try
  {
    const double val = PhysicalUnits::stringToDistance( txt );
    return (val > 0.0) ? val : 0.0;
  }catch( std::exception & )
  {
    return 0.0;
  }
}//distanceUncertainty()


void MakeDrfSrcDef::setDistanceUncertainty( const double uncert )
{
  if( uncert > 0.0 )
    m_distanceUncertEdit->setText( PhysicalUnits::printToBestLengthUnits(uncert) );
  else
    m_distanceUncertEdit->setText( "" );
  validate();
}//setDistanceUncertainty(...)


std::string MakeDrfSrcDef::assayInfoString() const
{
  string answer;
  try
  {
    if( m_useAgeInfo->isChecked() )
    {
      const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      answer += "Assay " + m_assayDate->date().toString("yyyy-MM-dd").toUTF8()
                + ": " + PhysicalUnits::printToBestActivityUnits( enteredActivity(), 5, useCi );
      const string age = m_sourceAgeAtAssay->text().toUTF8();
      if( !age.empty() )
        answer += "; age at assay " + age;
      answer += "; spectrum " + m_drfMeasurementDate->date().toString("yyyy-MM-dd").toUTF8();
    }
  }catch( std::exception & )
  {
  }
  return answer;
}//assayInfoString()


void MakeDrfSrcDef::setActivity( const double act )
{
  const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  const int ndecimals = 4;
  m_activityEdit->setText( PhysicalUnits::printToBestActivityUnits(act, ndecimals, useCi) );
  validate();
  updateAgedText();
}//void MakeDrfSrcDef::setActivity( const double dist )


void MakeDrfSrcDef::setAssayInfo( const double activity,
                                  const boost::posix_time::ptime &assay_date )
{
  // We only want to update the UI at the end of this function; intermediate updates are harmless
  //  but we'll emit once more at the end to ensure final state is propagated.
  m_useAgeInfo->setChecked( !assay_date.is_special() );
  m_assayDate->setDate( WDateTime::fromTime_t( boost::posix_time::to_time_t(assay_date) ).date() );
  useAgeInfoUserToggled();

  if( activity > 0.0 )
  {
    const int ndecimals = 4;
    const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    m_activityEdit->setText( PhysicalUnits::printToBestActivityUnits(activity, ndecimals, useCi) );
  }

  validate();
  updateAgedText();

  m_updated.emit();
}//void setAssayInfo(..);


/*
void MakeDrfSrcDef::setAgeAtAssay( const double age )
{
  if( age >= 0.0 )
  {
    m_useAgeInfo->setChecked( true );
    useAgeInfoUserToggled();
    if( !m_assayDate->date().isValid() )
      m_assayDate->setDate( m_drfMeasurementDate->date() );
    m_sourceAgeAtAssay->setText( PhysicalUnitsLocalized::printToBestTimeUnits(age) );
  }else
  {
    setNuclide( m_nuclide );
  }
  
  updateAgedText();
}//void setAgeAtAssay( const double age );
*/

void MakeDrfSrcDef::setAgeAtMeas( const double age )
{
  if( age < 0.0 )
  {
    setNuclide( m_nuclide );
    updateAgedText();
  }
  
  m_useAgeInfo->setChecked( true );
  useAgeInfoUserToggled();
    
  WDate measDate = m_drfMeasurementDate->date();
  WDate assayDate = m_assayDate->date();
    
  if( !measDate.isValid() && !assayDate.isValid() )
  {
    m_drfMeasurementDate->setDate( WDate::currentDate() );
    m_assayDate->setDate( WDate::currentDate() );
    measDate = assayDate = WDate::currentDate();
  }else if( !assayDate.isValid() )
  {
    m_assayDate->setDate( measDate );
    assayDate = measDate;
  }else if( !measDate.isValid() )
  {
    m_drfMeasurementDate->setDate( assayDate );
    measDate = assayDate;
  }else if( assayDate > measDate )
  {
    m_assayDate->setDate( measDate );
    assayDate = measDate;
  }
    
  double assayToMeasTime = assayDate.daysTo(measDate) * PhysicalUnits::day;
  
  if( assayToMeasTime > age )
  {
    m_assayDate->setDate( measDate );
    assayDate = measDate;
    assayToMeasTime = 0.0;
  }
  
  m_sourceAgeAtAssay->setText( PhysicalUnitsLocalized::printToBestTimeUnits(age-assayToMeasTime) );

  validate();
  updateAgedText();
}//void setAgeAtMeas( const double age );


void MakeDrfSrcDef::setShielding( const float atomic_number, const float areal_density )
{
  m_useShielding->setChecked();
  m_shieldingSelect->setHidden( false );
  m_shieldingSelect->setAtomicNumberAndArealDensity( atomic_number, areal_density );
  validate();
}//void setShielding( const float atomic_number, const float areal_density )


void MakeDrfSrcDef::setIsEffGeometryType( const int drf_eff_geom_type )
{
  const auto type = static_cast<DetectorPeakResponse::EffGeometryType>( drf_eff_geom_type );
  
  assert( (type == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram) );
  
  m_effGeomType = drf_eff_geom_type;

  const bool hide_distance = ((type != DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic)
                              && (type != DetectorPeakResponse::EffGeometryType::FarFieldAbsolute));

  updateActivityLabel();
  m_distanceEdit->setHidden( hide_distance );
  m_distanceLabel->setHidden( hide_distance );
  m_distanceUncertEdit->setHidden( hide_distance );
  m_distanceUncertLabel->setHidden( hide_distance );

  validate();
}//void setIsEffGeometryType( const bool is_fixed )


std::string MakeDrfSrcDef::toGadrasLikeSourceString() const
{
  string answer;
  if( m_nuclide )
  {
    //m_nuclide->symbol == 'U235m2'
    const size_t numpos = m_nuclide->symbol.find_first_of("0123456789");
    assert( numpos != string::npos );
    string numbers = m_nuclide->symbol.substr(numpos);  //235m2
    string meta;
    const size_t metapos = numbers.find_first_not_of("0123456789");
    if( metapos != string::npos )
    {
      meta = numbers.substr(metapos);  //m2
      numbers = numbers.substr(0,metapos); //235
    }
    
    answer = numbers + m_nuclide->symbol.substr(0,numpos) + meta;
  }//if( m_nuclide )
  
  if( !answer.empty() )
    answer += ",";
  
  const double activity = activityAtSpectrumTime();
  const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  answer += PhysicalUnits::printToBestActivityUnits(activity,5,useCi);
  
  if( m_shieldingSelect && m_useShielding->isChecked() )
  {
    double an = 14, ad = 0.0;
    if( m_shieldingSelect->isGenericMaterial() )
    {
      if( !m_shieldingSelect->atomicNumberEdit()->text().empty()
         && !m_shieldingSelect->arealDensityEdit()->text().empty() )
      {
        an = m_shieldingSelect->atomicNumber();
        ad = m_shieldingSelect->arealDensity();
      }
    }else
    {
      std::shared_ptr<const Material> mat = m_shieldingSelect->material();
      if( mat )
      {
        an = mat->massWeightedAtomicNumber();
        ad = mat->density * m_shieldingSelect->thickness();
      }//if( mat )
    }//if( shield->isGenericMaterial() ) / else
    
    // GADRAS convention, and what MakeDrf reads back: {atomic number, areal density in g/cm2}
    const double ad_g_cm2 = ad / (PhysicalUnits::gram / PhysicalUnits::cm2);
    answer += "{" + std::to_string(an) + "," + std::to_string(ad_g_cm2) + "}";
  }//if( user is using shielding )
  
  const bool mayEvolve = (m_nuclide && !PeakDef::ageFitNotAllowed(m_nuclide));
  const double age = ageAtSpectrumTime();
  if( mayEvolve && age >= 0.0 )
    answer += " Age=" + PhysicalUnitsLocalized::printToBestTimeUnits(age);
  
  try
  {
    // A fixed geometry has no distance (and distance() asserts against being asked)
    const double dist = m_distanceEdit->isHidden() ? 0.0 : distance();
    if( dist > 0.0 )
      answer += " Distance=" + PhysicalUnits::printToBestLengthUnits( dist, 6 );
    const double dist_uncert = distanceUncertainty();
    if( (dist > 0.0) && (dist_uncert > 0.0) )
      answer += " DistanceUncert=" + PhysicalUnits::printToBestLengthUnits( dist_uncert, 6 );
  }catch( std::exception &e )
  {
    // Fixed geometry or something
  }
  
  return answer;
}//std::string toGadrasLikeSourceString() const


const vector<SrcLibLineInfo> &MakeDrfSrcDef::lib_srcs_for_nuc()
{
  return m_lib_srcs_for_nuc;
}


void MakeDrfSrcDef::addSourceLibrary( vector<shared_ptr<const SrcLibLineInfo>> srcs,
                      const bool auto_populate )
{
  if( auto_populate && m_nuclide )
  {
    for( const shared_ptr<const SrcLibLineInfo> &info : srcs )
    {
      if( info->m_nuclide == m_nuclide )
      {
        setSrcInfo( *info);
        break;
      }
    }
  }//if( auto_populate )
  
  
  // Put new sources in front
  for( const shared_ptr<const SrcLibLineInfo> &info : m_lib_srcs_added )
  {
    if( !is_already_in(srcs, *info) )
      srcs.push_back( info );
  }
  
  m_lib_srcs_added.swap( srcs );
  
  // Update DB menu
  updateSourceLibNuclides();
}//void addSourceLibrary(...)
