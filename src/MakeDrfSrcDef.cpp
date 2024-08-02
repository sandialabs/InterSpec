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

#include <boost/algorithm/string/regex.hpp>

#include <Wt/WDate>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WComboBox>
#include <Wt/WDateEdit>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WMenuItem>
#include <Wt/WPopupMenu>
#include <Wt/WTableCell>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WDoubleSpinBox>
#include <Wt/WDoubleValidator>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>
#include <Wt/WSuggestionPopup>

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
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
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
  const int sm_activity_row        = 2;
  const int sm_activity_uncert_row = 3;
  const int sm_assay_date_row      = 4;
  const int sm_spec_date_row       = 5;
  const int sm_age_at_assay_row    = 6;
  const int sm_decayed_info_row    = 7;
  const int sm_shield_material_row = 8;
  const int sm_options_row         = 9;
  
  
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
  
  vector<string> base_paths{ InterSpec::staticDataDirectory(), SpecUtils::get_working_path() };
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_LOCAL_SERVER )
  try{ base_paths.push_back( InterSpec::writableDataDirectory() ); } catch(...){}
#endif
  vector<string> source_lib_files;
  for( const string &path : base_paths )
  {
    const vector<string> files = SpecUtils::recursive_ls(path, "Source.lib" );
    source_lib_files.insert( end(source_lib_files), begin(files), end(files) );
  }
  
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
              const boost::posix_time::ptime &measDate,
              MaterialDB *materialDB,
              Wt::WSuggestionPopup *materialSuggest,
              Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_table( nullptr ),
  m_nuclide( nuc ),
  m_materialDB( materialDB ),
  m_materialSuggest( materialSuggest ),
  m_nuclideLabel( nullptr ),
  m_distanceLabel( nullptr ),
  m_distanceEdit( nullptr ),
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
  m_updated( this )
{
  wApp->useStyleSheet( "InterSpec_resources/MakeDrfSrcDef.css" );
  
  addStyleClass( "MakeDrfSrcDef" );
  
  create();
  
  setNuclide( m_nuclide );
  
  if( !measDate.is_special() )
  {
    m_drfMeasurementDate->setDate( WDateTime::fromPosixTime(measDate).date() );
    validateDateFields();
  }
}//MakeDrfSrcDef constructor


MakeDrfSrcDef::~MakeDrfSrcDef()
{
#if( WT_VERSION >= 0x3070000 )
  if( m_lib_src_menu )
    delete m_lib_src_menu;
  m_lib_src_menu = nullptr;
#endif
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
      item->triggered().connect( boost::bind( &MakeDrfSrcDef::setSrcInfo, this, src) );
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
    m_nuclideLabel->setText( "Non-specified Nuclide" );
    m_useAgeInfo->setUnChecked();
    m_distanceEdit->setValueText( "25 cm" );
    
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    m_activityEdit->setValueText( useCi ? "1 uCi" : "37 kBq" );
    m_useAgeInfo->hide();
  }
  
  useAgeInfoUserToggled();
}//setNuclide(...)


void MakeDrfSrcDef::create()
{
  m_table = new WTable( this );
  m_table->addStyleClass( "SrcInputTable" );
  
  WTableCell *cell = m_table->elementAt(0,0);
  cell->setColumnSpan( 2 );
  cell->addStyleClass( "SrcNuclideTitle" );
  m_nuclideLabel = new WText( cell );
  
/*
  //Code to put a nuclide suggestion into a WLineEdit so user could select nuclide.
  string replacerJs, matcherJs;
  PhotopeakDelegate::EditWidget::replacerJs( replacerJs );
  PhotopeakDelegate::EditWidget::nuclideNameMatcherJs( matcherJs );
  
  WSuggestionPopup *nucSuggest = new WSuggestionPopup( matcherJs, replacerJs );
 #if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  nucSuggest->setJavaScriptMember("wtNoReparent", "true");
#endif
  
  nucSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  nucSuggest->forEdit( m_nuclideEdit,
                       WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );
  
  IsotopeNameFilterModel *filterModel = new IsotopeNameFilterModel( this );
  
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
  m_distanceLabel = new WLabel( "Distance", cell );
  cell = m_table->elementAt(sm_distance_row,1);
  m_distanceEdit = new WLineEdit( cell );
  m_distanceEdit->setTextSize( 16 );
  
  m_distanceEdit->setAutoComplete( false );
  m_distanceEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_distanceEdit->setAttributeValue( "autocorrect", "off" );
  m_distanceEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_distanceLabel->setBuddy( m_distanceEdit );
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_distanceEdit->setValidator( distValidator );
  m_distanceEdit->setText( "50 cm" );
  m_distanceEdit->changed().connect( this, &MakeDrfSrcDef::handleUserChangedDistance );
  m_distanceEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedDistance );
  
  
  cell = m_table->elementAt(sm_activity_row,0);
  m_activityLabel = new WLabel( "Activity", cell );
  
  cell = m_table->elementAt(sm_activity_row,1);
  m_activityEdit = new WLineEdit( cell );
  
  m_activityEdit->setAutoComplete( false );
  m_activityEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_activityEdit->setAttributeValue( "autocorrect", "off" );
  m_activityEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_activityEdit->setTextSize( 16 );
  m_activityLabel->setBuddy( m_activityEdit );
  
  WRegExpValidator *val = new WRegExpValidator( PhysicalUnits::sm_activityRegex, this );
  val->setFlags( Wt::MatchCaseInsensitive );
  m_activityEdit->setValidator( val );
  const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  m_activityEdit->setText( useCi ? "100 uCi" : "3.7 MBq" );
  m_activityEdit->changed().connect( this, &MakeDrfSrcDef::handleUserChangedActivity );
  m_activityEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedActivity );

  cell = m_table->elementAt(sm_activity_uncert_row,0);
  WLabel *label = new WLabel( "Act. Uncert.&nbsp;", cell );  //The nbsp is to make this the longest label so when acti ity or shielding is shown, the width doesnt get changed
  cell = m_table->elementAt(sm_activity_uncert_row,1);
  m_activityUncertainty = new WDoubleSpinBox( cell );
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
  
  m_activityUncertainty->changed().connect( this, &MakeDrfSrcDef::handleUserChangedActivityUncertainty );
  m_activityUncertainty->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedActivityUncertainty );
  
  //WDoubleValidator *percentVal = new WDoubleValidator( this );
  //percentVal->setRange( 0.0, 100.0 );
  //m_activityUncertainty->setValidator( percentVal );
  
  cell = m_table->elementAt(sm_assay_date_row,0);
  label = new WLabel( "Assay Date", cell );
  cell = m_table->elementAt(sm_assay_date_row,1);
  m_assayDate = new WDateEdit( cell );
  //m_assayDate->setTextSize( 9 );
  label->setBuddy( m_assayDate );
  m_assayDate->changed().connect( this, &MakeDrfSrcDef::handleEnteredDatesUpdated );
  
  cell = m_table->elementAt(sm_spec_date_row,0);
  label = new WLabel( "Spec. Date", cell );
  cell = m_table->elementAt(sm_spec_date_row,1);
  m_drfMeasurementDate = new WDateEdit( cell );
  //m_drfMeasurementDate->setTextSize( 10 );
  //The right padding is 40px, could reduce down to 30.
  label->setBuddy( m_drfMeasurementDate );
  m_drfMeasurementDate->changed().connect( this, &MakeDrfSrcDef::handleEnteredDatesUpdated );
  
  
  cell = m_table->elementAt(sm_age_at_assay_row,0);
  label = new WLabel( "Age@Assay", cell );
  cell = m_table->elementAt(sm_age_at_assay_row,1);
  m_sourceAgeAtAssay = new WLineEdit( cell );
  
  m_sourceAgeAtAssay->setAutoComplete( false );
  m_sourceAgeAtAssay->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_sourceAgeAtAssay->setAttributeValue( "autocorrect", "off" );
  m_sourceAgeAtAssay->setAttributeValue( "spellcheck", "off" );
#endif
  label->setBuddy( m_sourceAgeAtAssay );
  m_sourceAgeAtAssay->changed().connect( this, &MakeDrfSrcDef::handleUserChangedAgeAtAssay );
  m_sourceAgeAtAssay->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedAgeAtAssay );
  
  val = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationRegex(), this );
  val->setFlags( Wt::MatchCaseInsensitive );
  m_sourceAgeAtAssay->setValidator( val );
  m_sourceAgeAtAssay->setText( "0s" );
  
  
  cell = m_table->elementAt(sm_decayed_info_row,0);
  label = new WLabel( "Aging Res.", cell );
  cell = m_table->elementAt(sm_decayed_info_row,1);
  m_sourceInfoAtMeasurement = new WText( cell );
  
  //Wt::WDateEdit *m_sourceCreationDate;
  
  cell = m_table->elementAt(sm_options_row,0);
  
  m_lib_src_btn = new WPushButton( cell );
  m_lib_src_btn->setIcon( "InterSpec_resources/images/db_small.png" );
  m_lib_src_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  m_lib_src_btn->setFloatSide( Wt::Left );
  
#if( WT_VERSION < 0x3070000 )
  m_lib_src_menu = new PopupDivMenu( m_lib_src_btn, PopupDivMenu::TransientMenu );
  m_lib_src_menu->setJavaScriptMember("wtNoReparent", "true");
#else
  // If we have the button own the popup menu, the menu will be placed in our current div holding
  //  all src info, that may have a significant amount of scroll in it, so we will then have to
  //  scroll this div, to see all the menu items, which is annoying.
  //  So instead we'll make the menu a global widget, and just pop it up at the clicked location.
  //  This is a little less than optimal, and make auto-hide of menu not quite as good, but maybe
  //  better than the alternative.
  m_lib_src_menu = new PopupDivMenu( nullptr, PopupDivMenu::TransientMenu );
  m_lib_src_menu->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  m_lib_src_btn->clicked().connect( boost::bind( 
          static_cast<void (WPopupMenu::*)(const WMouseEvent &)>(&WPopupMenu::popup),
          m_lib_src_menu, boost::placeholders::_1 ) );
#endif
  
  m_lib_src_menu->setAutoHide( true, 2500 );
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  const char *tooltip = "Sources defined in Source.lib file in your users data directory.<br/>"
  "When clicked, this button will display a menu with all sources for this nuclide - and when"
  " on of the items is selected, its information will be populated.";
  HelpSystem::attachToolTipOn( m_lib_src_btn, tooltip, showToolTips, HelpSystem::ToolTipPosition::Right );
  
  m_useAgeInfo = new WCheckBox( "Age?", cell );
  m_useAgeInfo->setFloatSide( Wt::Right );
  m_useAgeInfo->setChecked( false );
  m_useAgeInfo->checked().connect( this, &MakeDrfSrcDef::useAgeInfoUserToggled );
  m_useAgeInfo->unChecked().connect( this, &MakeDrfSrcDef::useAgeInfoUserToggled );
  
  cell = m_table->elementAt(sm_options_row,1);
  m_useShielding = new WCheckBox( "Shielded?", cell );
  m_useShielding->setFloatSide( Wt::Right );
  m_useShielding->setChecked( false );
  m_useShielding->checked().connect( this, &MakeDrfSrcDef::useShieldingInfoUserToggled );
  m_useShielding->unChecked().connect( this, &MakeDrfSrcDef::useShieldingInfoUserToggled );
  
  
  cell = m_table->elementAt(sm_shield_material_row,0);
  cell->setColumnSpan( 2 );
  
  m_shieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest, cell );
  m_shieldingSelect->materialModified().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->materialChanged().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->addingIsotopeAsSource().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->removingIsotopeAsSource().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->hide();
  
  useAgeInfoUserToggled();
}//void create()


void MakeDrfSrcDef::useAgeInfoUserToggled()
{
  const bool useAge = m_useAgeInfo->isChecked();
  m_activityEdit->label()->setText( useAge ? "Assay Act." : "Activity" );
  
  m_table->rowAt(sm_assay_date_row)->setHidden( !useAge );
  m_table->rowAt(sm_spec_date_row)->setHidden( !useAge );
  
  const bool notMuchEvolution = (!m_nuclide || PeakDef::ageFitNotAllowed(m_nuclide));
  m_table->rowAt(sm_age_at_assay_row)->setHidden( !useAge || notMuchEvolution );
  
  m_table->rowAt(sm_decayed_info_row)->setHidden( !useAge );
  
  m_assayDate->setDisabled( !useAge );
  m_drfMeasurementDate->setDisabled( !useAge );
  m_sourceInfoAtMeasurement->setDisabled( !useAge );
  m_sourceAgeAtAssay->setDisabled( !useAge );
  
  if( useAge )
  {
    // If the text fields are empty, the red error background wont show up, so we'll at least put
    //  a space there
    if( m_assayDate->text().empty() )
      m_assayDate->setText( " " );
    
    if( m_drfMeasurementDate->text().empty() )
      m_drfMeasurementDate->setText( " " );
    
    validateDateFields();
  }//if( useAge )
  
  m_updated.emit();
}//void useAgeInfoUserToggled()


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
    m_sourceInfoAtMeasurement->setText( "" );
  }
}//void updateAgedText()


void MakeDrfSrcDef::handleUserChangedDistance()
{
  try
  {
    distance();
    if( m_distanceEdit->hasStyleClass( "SrcInputError" ) )
      m_distanceEdit->removeStyleClass( "SrcInputError" );
  }catch( std::exception & )
  {
    if( !m_distanceEdit->hasStyleClass( "SrcInputError" ) )
      m_distanceEdit->addStyleClass( "SrcInputError" );
  }
  
  m_updated.emit();
}//void handleUserChangedDistance()


void MakeDrfSrcDef::handleUserChangedActivity()
{
  try
  {
    enteredActivity();
    if( m_activityEdit->hasStyleClass( "SrcInputError" ) )
      m_activityEdit->removeStyleClass( "SrcInputError" );
  }catch( std::exception & )
  {
    if( !m_activityEdit->hasStyleClass( "SrcInputError" ) )
      m_activityEdit->addStyleClass( "SrcInputError" );
  }
  
  if( m_useAgeInfo->isChecked() )
    updateAgedText();
  
  m_updated.emit();
}//void handleUserChangedActivity()


void MakeDrfSrcDef::handleUserChangedActivityUncertainty()
{
  if( WValidator::State::Valid == m_activityUncertainty->validate() )
  {
    if( m_activityUncertainty->hasStyleClass( "SrcInputError" ) )
      m_activityUncertainty->removeStyleClass( "SrcInputError" );
  }else
  {
    if( !m_activityUncertainty->hasStyleClass( "SrcInputError" ) )
      m_activityUncertainty->addStyleClass( "SrcInputError" );
  }
  
  m_updated.emit();
}//void handleUserChangedActivityUncertainty();


void MakeDrfSrcDef::handleUserChangedAgeAtAssay()
{
  string agestr = m_sourceAgeAtAssay->text().toUTF8();
  SpecUtils::trim( agestr );
  

  if( agestr.empty() || (agestr.find_first_not_of("+-0.")==string::npos) )
  {
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    m_sourceAgeAtAssay->setText( useCi ? "0 uCi" : "0 bq" );
  }else
  {
    try
    {
      const double hl = m_nuclide ? m_nuclide->halfLife : -1.0;
      PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, hl );
      if( m_activityEdit->hasStyleClass( "SrcInputError" ) )
        m_activityEdit->removeStyleClass( "SrcInputError" );
    }catch( std::exception &e )
    {
      if( !m_activityEdit->hasStyleClass( "SrcInputError" ) )
        m_activityEdit->addStyleClass( "SrcInputError" );
    }
  }//if( zero / else )
  
  updateAgedText();
  
  m_updated.emit();
}//void handleUserChangedAgeAtAssay()


void MakeDrfSrcDef::validateDateFields()
{
  // Only validate if we are actually using the date fields
  if( !m_useAgeInfo->isChecked() )
    return;
  
  string txt = m_drfMeasurementDate->text().toUTF8();
  if( (txt.length() > 1) && std::isspace(txt[0]) )
    m_drfMeasurementDate->setText( WString::fromUTF8( SpecUtils::trim_copy(txt) ) );
  
  txt = m_assayDate->text().toUTF8();
  if( (txt.length() > 1) && std::isspace(txt[0]) )
    m_assayDate->setText( WString::fromUTF8( SpecUtils::trim_copy(txt) ) );
  
  if( m_drfMeasurementDate->validate() == Wt::WValidator::Valid )
  {
    if( m_drfMeasurementDate->hasStyleClass( "SrcInputError" ) )
      m_drfMeasurementDate->removeStyleClass( "SrcInputError" );
  }else
  {
    if( !m_drfMeasurementDate->hasStyleClass( "SrcInputError" ) )
      m_drfMeasurementDate->addStyleClass( "SrcInputError" );
  }
  
  if( m_assayDate->validate() == Wt::WValidator::Valid )
  {
    if( m_assayDate->hasStyleClass( "SrcInputError" ) )
      m_assayDate->removeStyleClass( "SrcInputError" );
  }else
  {
    if( !m_assayDate->hasStyleClass( "SrcInputError" ) )
      m_assayDate->addStyleClass( "SrcInputError" );
  }
}//void validateDateFields();


void MakeDrfSrcDef::handleEnteredDatesUpdated()
{
  validateDateFields();
  
  updateAgedText();
  
  m_updated.emit();
}//void handleEnteredDatesUpdated()


void MakeDrfSrcDef::useShieldingInfoUserToggled()
{
  m_shieldingSelect->setHidden( !m_useShielding->isChecked() );
  m_updated.emit();
}//void useShieldingInfoUserToggled()


void MakeDrfSrcDef::handleUserChangedShielding()
{
  m_updated.emit();
}


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
}//void setDistance( const double dist );


void MakeDrfSrcDef::setActivity( const double act )
{
  const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  const int ndecimals = 4;
  m_activityEdit->setText( PhysicalUnits::printToBestActivityUnits(act, ndecimals, useCi) );
  updateAgedText();
}//void MakeDrfSrcDef::setActivity( const double dist )


void MakeDrfSrcDef::setAssayInfo( const double activity,
                                  const boost::posix_time::ptime &assay_date )
{
  // We only want to update the UI at the end of this function - so block emitting updates until
  //  then (and actually we'll possible get an error in calculation if we dont wait until
  //  everything is updated)
  const bool updateBlocked = m_updated.isBlocked();
  m_updated.setBlocked( true );
  
  m_useAgeInfo->setChecked( !assay_date.is_special() );
  m_assayDate->setDate( WDateTime::fromPosixTime(assay_date).date() );
  useAgeInfoUserToggled();
  
  if( activity > 0.0 )
  {
    const int ndecimals = 4;
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    m_activityEdit->setText( PhysicalUnits::printToBestActivityUnits(activity, ndecimals, useCi) );
  }
  
  validateDateFields();
  updateAgedText();
  
  m_updated.setBlocked( updateBlocked );
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

  validateDateFields();
  updateAgedText();
}//void setAgeAtMeas( const double age );


void MakeDrfSrcDef::setShielding( const float atomic_number, const float areal_density )
{
  m_useShielding->setChecked();
  m_shieldingSelect->setHidden( false );
  m_shieldingSelect->setAtomicNumberAndArealDensity( atomic_number, areal_density );
}//void setShielding( const float atomic_number, const float areal_density )


void MakeDrfSrcDef::setIsEffGeometryType( const int drf_eff_geom_type )
{
  const auto type = static_cast<DetectorPeakResponse::EffGeometryType>( drf_eff_geom_type );
  
  assert( (type == DetectorPeakResponse::EffGeometryType::FarField)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2)
         || (type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram) );
  
  bool hide_distance = false;
  const char *text = "Activity";
  switch( type )
  {
    case DetectorPeakResponse::EffGeometryType::FarField:
      hide_distance = false;
      text = "Activity";
      break;
      
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      hide_distance = true;
      text = "Total Act.";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      hide_distance = true;
      text = "Act./cm2";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      hide_distance = true;
      text = "Act./m2";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      hide_distance = true;
      text = "Act./gram";
      break;
  }//switch( type )
  
  m_activityLabel->setText( text );
  m_distanceEdit->setHidden( hide_distance );
  m_distanceLabel->setHidden( hide_distance );
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
  const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
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
    
    answer += "{" + std::to_string(an) + "," + std::to_string(ad) + "}";
  }//if( user is using shielding )
  
  const bool mayEvolve = (m_nuclide && !PeakDef::ageFitNotAllowed(m_nuclide));
  const double age = ageAtSpectrumTime();
  if( mayEvolve && age >= 0.0 )
    answer += " Age=" + PhysicalUnitsLocalized::printToBestTimeUnits(age);
  
  try
  {
    const double dist = distance();
    if( dist > 0.0 )
      answer += " Distance=" + PhysicalUnits::printToBestLengthUnits( dist, 6 );
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
