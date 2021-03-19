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

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string.hpp>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WBreak>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WTable>
#include <Wt/WSignal>
#include <Wt/WServer>
#include <Wt/WRandom>
#include <Wt/WDateTime>
#include <Wt/WTextArea>
#include <Wt/WTreeView>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WIOService>
#include <Wt/WTableCell>
#include <Wt/WTabWidget>
#include <Wt/WFileUpload>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WRadioButton>
#include <Wt/WButtonGroup>
#include <Wt/WEnvironment>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/Dbo/QueryModel>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectorEdit.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/DetectorPeakResponse.h"


// The regex in GCC 4.8.x does not have working regex...., so we will detect this via
//    https://stackoverflow.com/questions/12530406/is-gcc-4-8-or-earlier-buggy-about-regular-expressions#answer-41186162
#if defined(_MSC_VER) \
    || (__cplusplus >= 201103L &&                             \
        (!defined(__GLIBCXX__) || (__cplusplus >= 201402L) || \
            (defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) || \
            defined(_GLIBCXX_REGEX_STATE_LIMIT)           || \
                (defined(_GLIBCXX_RELEASE)                && \
                _GLIBCXX_RELEASE > 4))))
#define HAVE_WORKING_REGEX 1
#else
#define HAVE_WORKING_REGEX 0
#endif

#if( HAVE_WORKING_REGEX )
#include <regex>
namespace RegexNs = std;
#else
#include <boost/regex.hpp>
#warning "DetectorEdit using boost regex - support for this compiler will be dropped soon"
namespace RegexNs = boost;
#endif


using namespace std;
using namespace Wt;

using SpecUtils::Measurement;
using SpecUtils::DetectorType;

const char * const DetectorDisplay::sm_noDetectorTxt
  = "<font style=\"font-weight:100;color:#CFCFCF;\">&lt;click to select&gt;</font>";
const char * const DetectorDisplay::sm_noDetectorTxtMbl
= "<font style=\"font-weight:100;color:#CFCFCF;\">&lt;tap to select&gt;</font>";

namespace
{
  void right_select_item( WMenu *menu, WMenuItem *item )
  {
    menu->select( item );
    item->triggered().emit( item ); //
  }
  
  class UtcToLocalTimeDelegate : public Wt::WItemDelegate
  {
    int64_t m_now;
    int m_timeZoneOffset;
  public:
    UtcToLocalTimeDelegate( Wt::WObject *parent = 0 )
    : WItemDelegate( parent )
    {
      m_now = std::time(nullptr);
      if( wApp )
        m_timeZoneOffset = wApp->environment().timeZoneOffset();
      else
        m_timeZoneOffset = 0;
    }
    virtual ~UtcToLocalTimeDelegate(){}
    virtual Wt::WWidget *update( Wt::WWidget *widget,
                                const Wt::WModelIndex &index,
                                Wt::WFlags<Wt::ViewItemRenderFlag > flags )
    {
      if( flags & RenderEditing )
        throw runtime_error( "UtcToLocalTimeDelegate not for editing" );
      
      if( !(flags & RenderEditing) )
      {
        WText *text = dynamic_cast<WText *>( widget );
        
        if( !text )
          widget = text = new WText();
        
        int64_t val;
        try
        {
          val = boost::any_cast<int64_t>( index.data() );
          string valstr;
          if( val > 0 )
          {
            auto ptt = boost::posix_time::from_time_t( time_t(val) );
            ptt += boost::posix_time::seconds( 60*m_timeZoneOffset );
            valstr = SpecUtils::to_common_string(ptt, true);
          }
          
          text->setText( WString::fromUTF8( valstr ) );
        }catch( std::exception &e )
        {
          cerr << "UtcToLocalTimeDelegate caught: " << e.what() << endl;
          text->setText( "" );
          return widget;
        }//try / catch
      }//if( !(flags & RenderEditing) )
      
      widget->setStyleClass( asString( index.data(StyleClassRole) ) );
      if( flags & RenderSelected )
        widget->addStyleClass(  "Wt-selected" );

      return widget;
    }
  };//class UtcToLocalTimeDelegate

  
  //Allow the file to be comma or tab delimited, but if an individual field
  //  contains a comma or tab, then the field must be quoted by a double
  //  quote.  Note that if you just copy cells from Microsoft Excel, that
  //  contain a comma, and then past into a text editor, fields with a comma
  //  are not quoted.
  void split_escaped_csv( vector<string> &fields, const string &line )
  {
    fields.clear();
    
    typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
    boost::escaped_list_separator<char> separator("\\",",\t", "\"");
    Tokeniser t( line, separator );
    for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
      fields.push_back(*it);
  }//void split_escaped_csv(...)
  
}//namespace

class RelEffFile;

/** This class represents potaentially many CSV or TSV relative efficincy
    detector response files (that each may have multiple DRFs).
    Uses the "RelativeEffDRFPaths" user option to get and store the filenames
    of the CSV/TSV files (note that filenames are semicolon delimited).
    This widget allows users to add or remove files.
 */
class RelEffDetSelect : public Wt::WContainerWidget
{
  friend class RelEffFile;
  
protected:
  InterSpec *m_interspec;
  DetectorEdit *m_detectorEdit;
  WContainerWidget *m_files; //holds the RelEffFile objects..
  
public:
  RelEffDetSelect( InterSpec *interspec, DetectorEdit *detedit, WContainerWidget *parent = 0 );
  
  virtual ~RelEffDetSelect(){}
  
  void userSelectedRelEffDetSelect();
  
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS))
  void addFile();
  void removeFile( RelEffFile *fileWidget );
  void saveFilePathToUserPreferences();
#endif
  
  void trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  void detectorSelected( RelEffFile *file, std::shared_ptr<DetectorPeakResponse> det );
  
  void docreate();
  
  //Attempts to do a defered rednering, but it looks like thisisnt actually
  //  the case
  virtual void load();
};//class RelEffDetSelect


/** Represents a single CSV/TSV DRF file, which may have multiple DRFs.
 */
class RelEffFile : public Wt::WContainerWidget
{
  //If we are deploying on the web, we do not want to allow the user access to
  //  the filesystem(!), so we will only allow using detectors from
  //  "data/OUO_lanl_simplemass_detectors.tsv"
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  std::string m_file;
#else
  WLineEdit *m_fileEdit;
  WPushButton *m_setFileButton;
#endif
  RelEffDetSelect *m_relEffDetSelect;
  DetectorEdit *m_detectorEdit;
  
  Wt::WText *m_credits;
  Wt::WComboBox *m_detectorSelect;
  std::vector<std::shared_ptr<DetectorPeakResponse> > m_responses;
  
public:
  RelEffFile( std::string file,
             RelEffDetSelect *parentSelect,
             DetectorEdit *detectorEdit,
             WContainerWidget *parent );
  
  
  ~RelEffFile(){}
  
  const std::string filepath();
  
  void selectNone();
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  void filePathChanged();
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )
  
  void detectorSelectCallback();
  
  bool trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  //parseDetector() returns null on error
  static std::shared_ptr<DetectorPeakResponse> parseDetector( string &line );
  void initDetectors();
  void detectorSelected( const int index );
};//class RelEffFile


class GadrasDirectory : public Wt::WContainerWidget
{
  friend class GadrasDetSelect;
  
protected:
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  std::string m_directory;
#else
  WLineEdit *m_directoryEdit;
  WPushButton *m_setDirectoryButton;
#endif
  GadrasDetSelect *m_parent;
  DetectorEdit *m_detectorEdit;
  
  Wt::WComboBox *m_detectorSelect;
  std::vector<std::shared_ptr<DetectorPeakResponse> > m_responses;

  Wt::WText *m_msg;
  Wt::WPushButton *m_deleteBtn;
  
  GadrasDirectory( std::string dorectory, GadrasDetSelect *parentSelect,
                   DetectorEdit *detectorEdit, WContainerWidget *parent );
  
  
  ~GadrasDirectory(){}
  
  std::string directory();
  
  void selectNone();
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  void dirPathChanged();
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )
  
  void detectorSelectCallback();
  
  bool trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  //parseDetector() returns null on error
  static std::shared_ptr<DetectorPeakResponse> parseDetector( const string directory );
  static vector<string> recursive_list_gadras_drfs( const string &sourcedir );
  
  void initDetectors();
  void detectorSelected( const int index );
};//class GadrasDirectory


/** Handles selecting a GADRAS DRF
 */
class GadrasDetSelect : public Wt::WContainerWidget
{
protected:
  InterSpec *m_interspec;
  DetectorEdit *m_detectorEdit;
  
  WContainerWidget *m_directories; //holds the GadrasDirectory objects..
  
public:
  GadrasDetSelect( InterSpec *interspec, DetectorEdit *detedit, WContainerWidget *parent = 0 );
  
  virtual ~GadrasDetSelect(){}
  
  bool trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  std::shared_ptr<DetectorPeakResponse> selectedDetector();
  
  void addDirectory();
  
  void removeDirectory( GadrasDirectory *dirWidget );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  void saveFilePathToUserPreferences();
#endif
  
  void detectorSelected( GadrasDirectory *dir, std::shared_ptr<DetectorPeakResponse> det );
  
  void docreate();
  
  //Attempts to do a defered rednering, (not actually verified)
  virtual void load();
};//class GadrasDetSelect


////////////////////////////////////////////////////////////////////////////////
//////////////////   Begin RelEffFile implementation   /////////////////////////
////////////////////////////////////////////////////////////////////////////////
RelEffFile::RelEffFile( std::string file,
                       RelEffDetSelect *parentSelect,
                       DetectorEdit *detectorEdit,
                       WContainerWidget *parent )
: WContainerWidget( parent ),
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
m_file( file ),
#else
m_fileEdit( new WLineEdit() ),
m_setFileButton( new WPushButton("Set") ),
#endif
m_relEffDetSelect( parentSelect ),
m_detectorEdit( detectorEdit ),
m_credits( nullptr )
{
  addStyleClass( "RelEffFile" );
  
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
#else
  WContainerWidget *topdiv = new WContainerWidget( this );
  WPushButton *closeIcon = new WPushButton( topdiv );
  closeIcon->addStyleClass( "closeicon-wtdefault" );
  closeIcon->setToolTip( "Remove file from list of files InterSpec will look for detector response functions in." );
  //closeIcon->setAttributeValue( "style", "position: relative; top: 3px; right: 3px;" + closeIcon->attributeValue("style") );
  closeIcon->clicked().connect( boost::bind( &RelEffDetSelect::removeFile, parentSelect, this ) );
  
  topdiv->addWidget( m_fileEdit );
  m_fileEdit->setText( file );
  m_fileEdit->setTextSize( 48 );
  
  topdiv->addWidget( m_setFileButton );
  m_setFileButton->setMargin( 5, Wt::Left );
  m_setFileButton->disable();
  m_setFileButton->setToolTip( "Parse the specified file for detector response functions." );
  m_setFileButton->clicked().connect( this, &RelEffFile::filePathChanged );
  
  //m_fileEdit->changed().connect( m_setFileButton, &WPushButton::enable );
  m_fileEdit->textInput().connect( m_setFileButton, &WPushButton::enable );
  m_fileEdit->enterPressed().connect( this, &RelEffFile::filePathChanged );
#endif
  
  WContainerWidget *bottomDiv = new WContainerWidget( this );
  
  m_credits = new WText( "", Wt::XHTMLText, bottomDiv );
  m_credits->setInline( false );
  m_credits->addStyleClass( "RelEffFileCredits" );
  
  m_detectorSelect = new WComboBox( bottomDiv );
  m_detectorSelect->activated().connect( this, &RelEffFile::detectorSelectCallback );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  m_fileEdit->changed().connect( m_detectorSelect, &WPushButton::disable );
  //m_fileEdit->keyPressed().connect( m_detectorSelect, &WPushButton::disable );
#endif
  
  initDetectors();
}//RelEffFile


const std::string RelEffFile::filepath()
{
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  return m_file;
#else
  return m_fileEdit->text().toUTF8();
#endif
}

void RelEffFile::selectNone()
{
  if( m_detectorSelect && m_detectorSelect->count() )
    m_detectorSelect->setCurrentIndex( 0 );
}//void selectNone()

#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
void RelEffFile::filePathChanged()
{
  m_relEffDetSelect->saveFilePathToUserPreferences();
  
  initDetectors();
  
  m_setFileButton->disable();
  m_detectorSelect->enable();
}//void filePathChanged()
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )


void RelEffFile::detectorSelectCallback()
{
  std::shared_ptr<DetectorPeakResponse> det;
 
  const int index = m_detectorSelect->currentIndex();
  
  //Item at index 0 is always a non-detector.
  if( index > 0 && (index-1) < m_responses.size() )
    det = m_responses[index-1];
  
  if( m_relEffDetSelect )
    m_relEffDetSelect->detectorSelected( this, det );
  
  if( m_detectorEdit )
  {
    m_detectorEdit->setDetector( det );
    m_detectorEdit->emitChangedSignal();
  }//if( m_detectorEdit )
}//void detectorSelectCallback()


bool RelEffFile::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  if( !det )
  {
    m_detectorSelect->setCurrentIndex( 0 );
    return false;
  }//if( !det )
  
  for( size_t i = 0; i < m_responses.size(); ++i )
  {
    if( m_responses[i]->hashValue() == det->hashValue() )
    {
      m_detectorSelect->setCurrentIndex( i+1 );
      return true;
    }
  }//for( size_t i = 0; i < m_responses.size(); ++i )
  
  m_detectorSelect->setCurrentIndex( 0 );
  return false;
}//void trySelectDetector(...)


//parseDetector() returns null on error
std::shared_ptr<DetectorPeakResponse> RelEffFile::parseDetector( string &line )
{
  std::shared_ptr<DetectorPeakResponse> det;
  
  vector<string> fields;
  fields.reserve( 20 );
  
  split_escaped_csv( fields, line );
  
  if( fields.size() < 16 )
    return det;
  
  try
  {
    const string name = fields[0] + " (" + fields[1] + ")";
    
    vector<float> coefs;
    for( int i = 3; i < 11; ++i )
      coefs.push_back( static_cast<float>( std::stod( fields[i] ) ) );
    
    //Get rid of the zero coefficients
    for( size_t i = coefs.size()-1; i > 0; --i )
    {
      if( fabs(coefs[i]) > 1.0E-14 ) //1.0E-14 chosen arbitrarily
      {
        coefs.erase( coefs.begin()+i+1, coefs.end() );
        break;
      }
    }//for( size_t i = coefs.size()-1; i > 0; --i )
    
    const float dist = static_cast<float>( std::stod(fields[14])*PhysicalUnits::cm );
    const float diam = 2.0f*static_cast<float>( std::stod(fields[15])*PhysicalUnits::cm );
    const float eunits = float(PhysicalUnits::MeV);
    
    string description = fields[2] + " - from Relative Eff. File";
    det.reset( new DetectorPeakResponse( name, description ) );
    det->fromExpOfLogPowerSeriesAbsEff( coefs, dist, diam, eunits, 0.0f, 0.0f );
    det->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedRelativeEfficiencyDrf );
  }catch( std::exception &e )
  {
    det.reset();
    cerr << "RelEffDetSelect::parseDetector() caught: " << e.what() << endl;
  }//try /catch
  
  return det;
}//static std::shared_ptr<DetectorPeakResponse> parseDetector( string &line )


void RelEffFile::initDetectors()
{
  m_responses.clear();
  m_credits->setText( "" );
  while( m_detectorSelect->count() )
    m_detectorSelect->removeItem( 0 );
  
  string pathstr;
  vector<string> credits;
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  pathstr = m_file;
#else
  pathstr = m_fileEdit->text().toUTF8();
#endif
  
#ifdef _WIN32
  const std::wstring wpathstr = SpecUtils::convert_from_utf8_to_utf16(pathstr);
  std::ifstream input( wpathstr.c_str(), ios::in | ios::binary );
#else
  std::ifstream input( pathstr.c_str(), ios::in | ios::binary );
#endif
  
  const bool file_opened = input.is_open();
  
  if( file_opened )
  {
    vector<vector<float>> detcoefs;
    string line;
    while( SpecUtils::safe_get_line( input, line, 2048 ) )
    {
      SpecUtils::trim( line );
      
      if( SpecUtils::istarts_with( line, "#credit:") )
        credits.push_back( SpecUtils::trim_copy(line.substr(8)) );
      
      if( line.empty() || line[0]=='#' )
        continue;
      
      std::shared_ptr<DetectorPeakResponse> det = parseDetector( line );
      
      if( det )
      {
        detcoefs.push_back( det->efficiencyExpOfLogsCoeffs() );
        m_responses.push_back( det );
      }
    }//while( SpecUtils::safe_get_line( input, line ) )
    
    if( m_responses.empty() )
    {
      m_detectorSelect->addItem( "<no responses in file>" );
      m_detectorSelect->hide();
    }else
    {
      m_detectorSelect->addItem( "<select detector>" );
      m_detectorSelect->show();
    }
    
    if( m_responses.empty() )
      credits.push_back( "<span style=\"color:red;\">No valid DRFs in file.</span>" );
    
    for( const auto &det : m_responses )
      m_detectorSelect->addItem( det->name() );
  }else
  {
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
    credits.push_back( "<span style=\"color:red;\">No Rel. Eff. DRFs available.</span>" );
#else
    credits.push_back( "<span style=\"color:red;\">Could not open DRF file.</span>" );
#endif
    m_detectorSelect->addItem( "<no responses available>" );
    m_detectorSelect->hide();
  }//if( file_opened )
  
  if( file_opened && m_responses.empty() )
    credits.push_back( "<span style=\"color:red;\">File does not contain and DRFs</span>" );
  
  string creditHtml;
  for( const string &credit : credits )
    creditHtml += "<div>" + credit + "</div>";
  
  m_credits->setText( creditHtml );
  
  m_detectorSelect->setCurrentIndex( 0 );
}//initDetectors()

void RelEffFile::detectorSelected( const int index )
{
  std::shared_ptr<DetectorPeakResponse> det;
  if( index > 0 )
    det = m_responses.at( index - 1 );
  
  if( !m_detectorEdit )
    return;
  
  m_detectorEdit->setDetector( det );
  m_detectorEdit->emitChangedSignal();
}//void detectorSelected( const int index )
////////////////////////////////////////////////////////////////////////////////
///////////////////   End RelEffFile implementation   //////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////   Begin RelEffDetSelect implementation   //////////////////////
////////////////////////////////////////////////////////////////////////////////
RelEffDetSelect::RelEffDetSelect( InterSpec *interspec, DetectorEdit *detedit, WContainerWidget *parent )
    : WContainerWidget( parent ),
      m_interspec( interspec ),
      m_detectorEdit( detedit ),
      m_files( nullptr )
{
}


#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
void RelEffDetSelect::addFile()
{
  if( m_files )
    new RelEffFile( "", this, m_detectorEdit, m_files );
}//void removeFile( RelEffFile )


void RelEffDetSelect::removeFile( RelEffFile *fileWidget )
{
  if( !m_files || !fileWidget )
    return;
  
  for( WWidget *w : m_files->children() )
  {
    if( dynamic_cast<RelEffFile *>(w) == fileWidget )
    {
      delete w;
      break;
    }
  }//for( WWidget *w : m_files->children() )
  
  saveFilePathToUserPreferences();
}//void removeFile( RelEffFile )


void RelEffDetSelect::saveFilePathToUserPreferences()
{
  if( !m_files )
    return;
  
  auto children = m_files->children();
  
  vector<string> paths;
  for( auto w : children )
  {
    auto child = dynamic_cast<RelEffFile *>( w );
    const string filepath = child ? child->filepath() : string("");
    if( filepath.size() )
      paths.push_back( filepath );
      //concat_path += (concat_path.size() ? ";" : "") + filepath;
  }//string concat_path;
  
  //Make sure we only save unique paths
  paths.erase( std::unique( begin(paths), end(paths) ), end(paths) );

  string concat_path;
  for( const auto &filepath : paths )
    concat_path += (concat_path.size() ? ";" : "") + filepath;

  try
  {
    InterSpecUser::setPreferenceValue( m_interspec->m_user, "RelativeEffDRFPaths", concat_path, m_interspec );
  }catch( std::exception &e )
  {
    passMessage( "Error saving Rel. Eff. files to user preferences: " + string(e.what()), "DetectorEdit", WarningWidget::WarningMsgHigh );
  }
}//void saveFilePathToUserPreferences()
#endif //#if( !BUILD_FOR_WEB_DEPLOYMENT )


void RelEffDetSelect::detectorSelected( RelEffFile *file, std::shared_ptr<DetectorPeakResponse> det )
{
  if( !m_files || !file )
    return;
  
  for( WWidget *w : m_files->children() )
  {
    auto child = dynamic_cast<RelEffFile *>(w);
    if( child && (child != file) )
      child->selectNone();
  }//for( WWidget *w : m_files->children() )
}//void detectorSelected( RelEffFile *file, std::shared_ptr<DetectorPeakResponse> det );


void RelEffDetSelect::userSelectedRelEffDetSelect()
{
  
}//void userSelectedRelEffDetSelect()


void RelEffDetSelect::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  docreate();
  
  if( !m_files )
    return;
  
  auto children = m_files->children();
  
  bool found = false;
  for( auto w : children )
  {
    auto child = dynamic_cast<RelEffFile *>( w );
    if( !child )
      continue;
    
    if( found )
      child->selectNone();
    else
      found = child->trySelectDetector(det);
  }//for( size_t i = 0; i < children.size(); ++i )
}//void trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )


void RelEffDetSelect::docreate()
{
  if( m_files )
    return;
  
  m_files = new WContainerWidget( this );
  WContainerWidget *holder = new WContainerWidget( this );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  WPushButton *addIcon = new WPushButton();
  addIcon->setStyleClass( "AddRelEffFile Wt-icon" );
  addIcon->clicked().connect( this, &RelEffDetSelect::addFile );
  holder->addWidget( addIcon );
#endif
  
  holder->setWidth( WLength(100,WLength::Percentage) );
  holder->setToolTip( "Click the plus button to add an additonal relative efficiency detector response function file." );
  
  string pathstr;
  vector<string> paths;
#if( BUILD_FOR_WEB_DEPLOYMENT )
  pathstr = SpecUtils::append_path( InterSpec::staticDataDirectory(), "OUO_lanl_simplemass_detectors.tsv" );
#elif( !defined(IOS) )
  try
  {
    if( m_interspec )
      pathstr = InterSpecUser::preferenceValue<string>( "RelativeEffDRFPaths", m_interspec );
  }catch( std::exception & )
  {
    passMessage( "Error retrieving 'RelativeEffDRFPaths' preference.", "", WarningWidget::WarningMsgHigh );
  }
#endif

  
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__)) ) )
  try
  {
    const string userDir = InterSpec::writableDataDirectory();
    const vector<string> tsv_files = SpecUtils::recursive_ls( userDir, ".tsv");
    for( const auto &p : tsv_files )
      pathstr += (pathstr.empty() ? "" : ";") + p;
  }catch( std::exception & )
  {
    cerr << "Couldnt call into InterSpec::writableDataDirectory()" << endl;
  }
#endif
  
  SpecUtils::split( paths, pathstr, "\r\n;" );
  
  //Eliminate any duplicate entries in paths
  paths.erase( std::unique( begin(paths), end(paths) ), end(paths) );

  if( paths.empty() )
    new RelEffFile( "", this, m_detectorEdit, m_files );
  
  for( const string &path : paths )
    new RelEffFile( path, this, m_detectorEdit, m_files );
  
}//void docreate()


//Attempts to do a defered rednering, but it looks like thisisnt actually
//  the case
void RelEffDetSelect::load()
{
  if( !loaded() )
    docreate();
  WContainerWidget::load();
}//void load()
////////////////////////////////////////////////////////////////////////////////
/////////////////   End RelEffDetSelect implementation   ///////////////////////
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
////////////////   Begin GadrasDetSelect implementation   //////////////////////
////////////////////////////////////////////////////////////////////////////////


///************************ Begin Need to implement for GADRAS **************/////
GadrasDetSelect::GadrasDetSelect( InterSpec *interspec, DetectorEdit *detedit, WContainerWidget *parent )
  : WContainerWidget(parent),
    m_interspec( interspec ),
    m_detectorEdit( detedit ),
    m_directories( nullptr )
{
  //GADRAS-DRF URL: https://rsicc.ornl.gov/codes/psr/psr6/psr-610.html
}


void GadrasDetSelect::load()
{
  if( !loaded() )
    docreate();
  WContainerWidget::load();
}//void load()


void GadrasDetSelect::docreate()
{
  if( m_directories )
    return;
  
  m_directories = new WContainerWidget( this );
  
  //Need to point the GUI to the appropriate directory, and implement to an `ls` to find detctors with both Detcotr.dat and Efficy.csv.
  const string drfpaths = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", m_interspec );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  WPushButton *addIcon = new WPushButton();
  addIcon->setStyleClass( "AddRelEffFile Wt-icon" );
  addIcon->clicked().connect( this, &GadrasDetSelect::addDirectory );
  WContainerWidget *holder = new WContainerWidget( this );
  holder->addWidget( addIcon );
  holder->setWidth( WLength(100,WLength::Percentage) );
  holder->setToolTip( "Click the plus button to add an additonal directory to recursively look for detector response functions in." );
#endif
  
  string pathstr;
  vector<string> paths;
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  try
  {
    if( m_interspec )
      pathstr = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", m_interspec );
  }catch( std::exception & )
  {
    passMessage( "Error retrieving 'RelativeEffDRFPaths' preference.", "", WarningWidget::WarningMsgHigh );
  }
#endif

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__)) ) )
  try
  {
    //ToDo: do more testing and use only the Android implementation
#if( ANDROID )
    const string basestr = InterSpec::writableDataDirectory();
    const vector<string> subdirs = SpecUtils::ls_directories_in_directory( basestr );
    for( const string &subdirpath : subdirs )
    {
      auto subsubdirs = GadrasDirectory::recursive_list_gadras_drfs(subdirpath);
      if( subsubdirs.size() )
        pathstr += (pathstr.empty() ? "" : ";") + subdirpath;
    }//for( const string &subdir : subdirs )
#else
    using namespace boost::filesystem;
    auto itr = directory_iterator( InterSpec::writableDataDirectory() );
    
    for( ; itr != directory_iterator(); ++itr )
    {
      const boost::filesystem::path &p = itr->path();
      const string pstr = p.string<string>();
      if( SpecUtils::is_directory( pstr ) )
      {
        auto subdirs = GadrasDirectory::recursive_list_gadras_drfs(pstr);
        if( subdirs.size() )
          pathstr += (pathstr.empty() ? "" : ";") + pstr;
      }
    }//for( loop over
#endif //if ANDROID / else
  }catch( std::exception &e )
  {
    cerr << "Got exception looking for GADRAS DRFs in user document dir: " << e.what() << endl;
  }//try / catch
#endif
  
  SpecUtils::split( paths, pathstr, "\r\n;" );
  
  //Make sure we always at least have the default generic detectors available.
  bool hasGeneric = false;
  for( size_t i = 0; !hasGeneric && (i < paths.size()); ++i )
    hasGeneric = (paths[i].find("GenericGadrasDetectors") != string::npos);
  
  if( !hasGeneric )
  {
    const string drfpaths = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GenericGadrasDetectors" );
    paths.push_back( drfpaths );
  }
  
    
  for( const string &path : paths )
  {
    auto dir = new GadrasDirectory( path, this, m_detectorEdit, m_directories );
    
    if( path.find("GenericGadrasDetectors") != string::npos )
    {
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
      dir->m_directoryEdit->disable();
      dir->m_setDirectoryButton->hide();
      dir->m_deleteBtn->hide();
#endif
    }
    
#if( !defined(WIN32) )
    if( SpecUtils::istarts_with( path, "C:\\" ) )
      dir->hide();
#endif
  }

#if( IOS )
  const char *gadrasToolTip = "You can place detector response functions imported"
  " from GADRAS-DRF into <code>InterSpec</code>s directory within the <b>Files</b> app.";
#else
  const char *gadrasToolTip = "These are detector response functions imported"
  " from GADRAS-DRF.  Currently only the detector dimensions, "
  " absolute efficiencies, and FWHM resolutions are used.";
#endif
  
  WText *useInfo = new WText( gadrasToolTip, this );
  useInfo->setStyleClass("DetectorLabel");
  useInfo->setMargin( 5, Wt::Top );
  useInfo->setInline( false );
}//void docreate()



bool GadrasDetSelect::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  docreate();  //We need to initialize the widget in order for m_directories to be true
  
  if( !m_directories )
    return false;
  
  if( !det || det->name().empty() )
  {
    for( auto w : m_directories->children() )
    {
      GadrasDirectory *d = dynamic_cast<GadrasDirectory *>( w );
      if( !d )
        continue;  //shouldnt ever happen, but JIC
      d->m_detectorSelect->setCurrentIndex( -1 );
    }//for( auto w : m_directories->children() )
    
    return false;
  }//if( !det || det->name().empty() )
  
  for( auto w : m_directories->children() )
  {
    GadrasDirectory *d = dynamic_cast<GadrasDirectory *>( w );
    if( !d )
      continue;  //shouldnt ever happen, but JIC
    
    for( size_t i = 0; i < d->m_responses.size(); ++i )
    {
      auto response = d->m_responses[i];
      
      if( response == det || det->hashValue()==response->hashValue() )
      {
        d->m_detectorSelect->setCurrentIndex( i+1 );
        return true;
      }
    }//for( size_t i = 0; i < d->m_responses.size(); ++i )
    
    d->m_detectorSelect->setCurrentIndex( 0 );
  }//for( auto w : m_directories->children() )
  
  return false;
}//bool trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )


std::shared_ptr<DetectorPeakResponse> GadrasDetSelect::selectedDetector()
{
  if( !m_directories )
    return nullptr;
  
  for( auto w : m_directories->children() )
  {
    GadrasDirectory *d = dynamic_cast<GadrasDirectory *>( w );
    const int currentIndex = d ? d->m_detectorSelect->currentIndex() : -1;
    if( currentIndex >= 1 )
    {
      try
      {
        WAbstractItemModel *m = d->m_detectorSelect->model();
        const string p = boost::any_cast<std::string>( m->data( currentIndex, 0, Wt::UserRole ) );
        auto answer = DetectorEdit::initAGadrasDetectorFromDirectory( p, m_interspec );
        
        if( p.find("GenericGadrasDetectors") != string::npos )
          answer->setDrfSource( DetectorPeakResponse::DrfSource::DefaultGadrasDrf );
        else
          answer->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedGadrasDrf );
        
        return answer;
      }catch( std::exception &e )
      {
        passMessage( "Failed to parse a GADRAS detector: " + string(e.what()), "", WarningWidget::WarningMsgHigh );
      }//try / catch
    }//if( user has selected an index )
  }//for( auto w : m_directories->children() )
  
  return nullptr;
}//std::shared_ptr<DetectorPeakResponse> selectedDetector()


void GadrasDetSelect::addDirectory()
{
  if( m_directories )
    new GadrasDirectory( "", this, m_detectorEdit, m_directories );
}//void addDirectory()


void GadrasDetSelect::removeDirectory( GadrasDirectory *dirWidget )
{
  if( !m_directories || !dirWidget )
    return;
  
  for( WWidget *w : m_directories->children() )
  {
    if( dynamic_cast<GadrasDirectory *>(w) == dirWidget )
    {
      delete w;
      break;
    }
  }//for( WWidget *w : m_files->children() )
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  saveFilePathToUserPreferences();
#endif
}//void removeDirectory( GadrasDirectory *dirWidget )


#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
void GadrasDetSelect::saveFilePathToUserPreferences()
{
  if( !m_directories )
    return;
  
  string dirs;
  for( WWidget *w : m_directories->children() )
  {
    auto d = dynamic_cast<GadrasDirectory *>(w);
    if( !d )
      continue;
    const string p = d->m_directoryEdit->text().toUTF8();
    dirs += ((dirs.size() && !p.empty()) ? ";" : "") + p;
  }//for( WWidget *w : m_directories->children() )
  
  try
  {
    InterSpecUser::setPreferenceValue( m_interspec->m_user, "GadrasDRFPath", dirs, m_interspec );
  }catch( std::exception &e )
  {
    passMessage( "Error saving GADRAS path to user preferences: " + string(e.what()), "", WarningWidget::WarningMsgHigh );
  }
}//void saveFilePathToUserPreferences()
#endif


void GadrasDetSelect::detectorSelected( GadrasDirectory *dir, std::shared_ptr<DetectorPeakResponse> det )
{
  if( !m_directories || !dir )
    return;
  
  for( WWidget *w : m_directories->children() )
  {
    auto child = dynamic_cast<GadrasDirectory *>(w);
    if( child && (child != dir) )
      child->selectNone();
  }//for( WWidget *w : m_files->children() )
}//void detectorSelected( GadrasDirectory *dir, std::shared_ptr<DetectorPeakResponse> det )


GadrasDirectory::GadrasDirectory( std::string directory, GadrasDetSelect *parentSelect,
                  DetectorEdit *detectorEdit, WContainerWidget *parent )
: WContainerWidget( parent ),
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  m_directory( directory ),
#else
  m_directoryEdit( new WLineEdit( Wt::WString::fromUTF8(directory) ) ),
  m_setDirectoryButton( new WPushButton("Set") ),
#endif
  m_parent( parentSelect ),
  m_detectorEdit( detectorEdit ),
  m_detectorSelect( new WComboBox() ),
  m_msg( new WText() ),
  m_deleteBtn( new WPushButton() )
{
  setObjectName( "GadDir" + Wt::WRandom::generateId() );
  
  addStyleClass( "RelEffFile" );
  
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
#else
  WContainerWidget *topdiv = new WContainerWidget( this );
  topdiv->addWidget( m_deleteBtn );
  m_deleteBtn->addStyleClass( "closeicon-wtdefault" );
  m_deleteBtn->setToolTip( "Remove directory from list of directories InterSpec will look in for detector response functions in." );
  m_deleteBtn->clicked().connect( boost::bind( &GadrasDetSelect::removeDirectory, parentSelect, this ) );
  
  
//If we wanted to actually select directories, could do similar to file query widget... which webkitdirectory no longer seems to work, so see file query widget use of electrons dialog
//#if( BUILD_AS_ELECTRON_APP )
//  m_pathSelectedSignal.reset( new Wt::JSignal<std::string>( this, "BaseDirSelected", false ) );
//  const string uploadname = id() + "PathPicker";
//  const string uploadhtml = "<input id=\"" + uploadname + "\" type=\"file\" webkitdirectory=\"\" />";
  
//  WText *uploadtext = new WText( uploadhtml, XHTMLUnsafeText );
//  linelayout->addWidget( uploadtext, 0, 1 );
  
  //TODO: put in error handling!
//  wApp->doJavaScript( "document.getElementById('" + uploadname + "').onchange = function(event){"
//                     "var outputDir = document.getElementById('" + uploadname + "').files[0].path;"
//                     "Wt.emit( \"" + id() + "\", { name: 'BaseDirSelected' }, outputDir );"
//                     "};"
//                     );
  //m_pathSelectedSignal->connect( boost::bind( &SpecFileQueryWidget::newElectronPathSelected, this, _1 ) );
//#elif( BUILD_AS_OSX_APP )
//  SpecFileQuery::setIsSelectingDirectory( true );
//  setSearchDirectory( "" );
//  m_baseLocation = new WFileUpload();
//  m_baseLocation->changed().connect( this, &SpecFileQueryWidget::newMacOsPathSelected );
//  linelayout->addWidget( m_baseLocation, 0, 1 );
//#else
  
  //TODO: Check if 'directory' is a subdirectory of InterSpec::writableDataDirectory()
  //and if so replace the beggingin of directory with {DataDir}/... and make it
  //  so the edit cant be changed and set button is hidden, and then deal with
  //  in GadrasDirectory::directory()
/*
  bool isInDataDir = false;
  string dirCononical = directory;
  string appDataDir = InterSpec::writableDataDirectory();
  if( SpecUtils::make_canonical_path(appDataDir)
     && SpecUtils::make_canonical_path(dirCononical)
     && appDataDir.size()>2 && dirCononical.size()>2 )
  {
    while( !isInDataDir && (dirCononical.size()+1) >= appDataDir.size() )
    {
      if( dirCononical == appDataDir )
      {
        isInDataDir = true;
        dirCononical = directory;
        SpecUtils::make_canonical_path(dirCononical);
        dirCononical = fs_relative( appDataDir, dirCononical );
        dirCononical = SpecUtils::append_path( "{AppDataDir}", dirCononical );
      }else
      {
        dirCononical = SpecUtils::parent_path(dirCononical);
        SpecUtils::make_canonical_path( dirCononical ); //to make sure we get ending / or \ right
      }
    }
  }
*/
  
  
/*
 ToDo: custonize path picking for electron...
#if( BUILD_AS_ELECTRON_APP )
  m_pathSelectedSignal.reset( new Wt::JSignal<std::string>( this, "BaseDirSelected", false ) );
  
  const string uploadname = id() + "PathPicker";
  const string uploadhtml = "<input id=\"" + uploadname + "\" type=\"file\" webkitdirectory=\"\" />";
  
  WText *uploadtext = new WText( uploadhtml, XHTMLUnsafeText );
  linelayout->addWidget( uploadtext, 0, 1 );
  
  //TODO: put in error handling!
  wApp->doJavaScript( "document.getElementById('" + uploadname + "').onchange = function(event){"
                     "var outputDir = document.getElementById('" + uploadname + "').files[0].path;"
                     "Wt.emit( \"" + id() + "\", { name: 'BaseDirSelected' }, outputDir );"
                     "};"
                     );
  m_pathSelectedSignal->connect( boost::bind( &SpecFileQueryWidget::newElectronPathSelected, this, _1 ) );
#elif( BUILD_AS_OSX_APP )
 //For macOS dont saved picked directory to preferences in DB as sandboxing will mess this up.
  SpecFileQuery::setIsSelectingDirectory( true );
  setSearchDirectory( "" );
  m_baseLocation = new WFileUpload();
  m_baseLocation->changed().connect( this, &SpecFileQueryWidget::newMacOsPathSelected );
  linelayout->addWidget( m_baseLocation, 0, 1 );
#else
*/
  
  topdiv->addWidget( m_directoryEdit );
  m_directoryEdit->setText( directory );
  m_directoryEdit->setTextSize( 48 );
  
  topdiv->addWidget( m_setDirectoryButton );
  m_setDirectoryButton->setMargin( 5, Wt::Left );
  m_setDirectoryButton->disable();
  m_setDirectoryButton->setToolTip( "Set entered directory to recursively look in for detector response functions." );
  m_setDirectoryButton->clicked().connect( this, &GadrasDirectory::dirPathChanged );
  
  //m_fileEdit->changed().connect( m_setFileButton, &WPushButton::enable );
  m_directoryEdit->textInput().connect( m_setDirectoryButton, &WPushButton::enable );
  m_directoryEdit->enterPressed().connect( this, &GadrasDirectory::dirPathChanged );
#endif
  
  WContainerWidget *bottomDiv = new WContainerWidget( this );
  
  m_detectorSelect = new WComboBox( bottomDiv );
  m_detectorSelect->activated().connect( this, &GadrasDirectory::detectorSelectCallback );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  m_directoryEdit->changed().connect( m_detectorSelect, &WPushButton::disable );
#endif
  
  m_msg->setInline( false );
  bottomDiv->addWidget( m_msg );
  
  
  initDetectors();
}//GadrasDirectory constructor


std::string GadrasDirectory::directory()
{
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  return m_directory;
#else
  return m_directoryEdit->text().toUTF8();
#endif
}


void GadrasDirectory::selectNone()
{
  if( m_detectorSelect->count() )
    m_detectorSelect->setCurrentIndex( 0 );
}


#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
void GadrasDirectory::dirPathChanged()
{
  m_parent->saveFilePathToUserPreferences();
  
  initDetectors();
  
  m_setDirectoryButton->disable();
  m_detectorSelect->enable();
}//void dirPathChanged()
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )


void GadrasDirectory::detectorSelectCallback()
{
  std::shared_ptr<DetectorPeakResponse> det;
  
  const int index = m_detectorSelect->currentIndex();
  
  //Item at index 0 is always a non-detector.
  if( index > 0 && (index-1) < m_responses.size() )
    det = m_responses[index-1];
  
  if( m_parent )
    m_parent->detectorSelected( this, det );
  
  m_detectorEdit->setDetector( det );
  m_detectorEdit->emitChangedSignal();
}//void detectorSelectCallback()

  
bool GadrasDirectory::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  if( !det )
  {
    m_detectorSelect->setCurrentIndex( 0 );
    return false;
  }//if( !det )
  
  for( size_t i = 0; i < m_responses.size(); ++i )
  {
    if( m_responses[i]->hashValue() == det->hashValue() )
    {
      m_detectorSelect->setCurrentIndex( i+1 );
      return true;
    }
  }//for( size_t i = 0; i < m_responses.size(); ++i )
  
  m_detectorSelect->setCurrentIndex( 0 );
  return false;
}//bool trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )


//parseDetector() returns null on error
std::shared_ptr<DetectorPeakResponse> GadrasDirectory::parseDetector( string path )
{
  try
  {
    const string name = SpecUtils::filename( path );
    auto drf = make_shared<DetectorPeakResponse>( name, "" );
    drf->fromGadrasDirectory( path );
    
    return drf;
  }catch( std::exception &e )
  {
    cerr << "GadrasDirectory::parseDetector() caught: " << e.what() << endl;
  }
 
  return nullptr;
}//shared_ptr<DetectorPeakResponse> parseDetector( const string directory )


vector<string> GadrasDirectory::recursive_list_gadras_drfs( const string &sourcedir )
{
  vector<string> files;
  if( !SpecUtils::is_directory( sourcedir ) )
    return files;
  
  const string csv_file = SpecUtils::append_path( sourcedir, "Efficiency.csv");
  const string dat_file = SpecUtils::append_path( sourcedir, "Detector.dat");
  
  if( SpecUtils::is_file(csv_file) && SpecUtils::is_file(dat_file) )
    files.push_back( sourcedir );
  
  //ToDo: Maybe we should use the Android implemenatation always.
#if( ANDROID )
  const vector<string> subdirs = SpecUtils::ls_directories_in_directory( sourcedir );
  for( const string &subdirpath : subdirs )
  {
    auto subsubdirs = recursive_list_gadras_drfs( subdirpath );
    files.insert( end(files), begin(subsubdirs), end(subsubdirs) );
  }//for( const string &subdir : subdirs )
#else
using namespace boost::filesystem;
  directory_iterator end_itr; // default construction yields past-the-end
  
  directory_iterator itr;
  try
  {
    itr = directory_iterator( sourcedir );
  }catch( std::exception & )
  {
    //ex: boost::filesystem::filesystem_error: boost::filesystem::directory_iterator::construct: Permission denied: "..."
    return files;
  }
  
  for( ; itr != end_itr; ++itr )
  {
    const boost::filesystem::path &p = itr->path();
    const string pstr = p.string<string>();
    if( SpecUtils::is_directory( pstr ) )
    {
      auto subdirs = recursive_list_gadras_drfs(pstr);
      files.insert( end(files), begin(subdirs), end(subdirs) );
    }
  }//for( loop over
#endif //if ANDROID / else
  
  return files;
}//vector<string> recursive_list_gadras_drfs( const string &sourcedir )


void GadrasDirectory::initDetectors()
{
  m_detectorSelect->clear();
  m_responses.clear();
  
  m_msg->setText( "" );
  m_msg->hide();
  
  const std::string basedir = directory();
  
  if( !SpecUtils::is_directory( basedir ) )
  {
    m_msg->setText( "<span style=\"color:red;\">Not a valid directory.</span>" );
    m_msg->show();
    m_detectorSelect->addItem( "<invalid directory>" );
    m_detectorSelect->setCurrentIndex( 0 );
    m_detectorSelect->disable();
    return;
  }//if( !SpecUtils::is_directory( basedir ) )
  
  m_detectorSelect->disable();
  this->disable();
  
  
  const std::string objname = objectName();
  const std::string sessid = wApp->sessionId();
  
  auto updategui = [this,objname]( std::vector<std::shared_ptr<DetectorPeakResponse> > drfs ){
    
    //Lets make sure this widget hasnt been deleted, by looking for it in the DOM
    WWidget *w = wApp->findWidget(objname);
    if( !w )
    {
      cerr << "Couldnt find widget '" << objname << "' in the DOM." << endl;
      return;
    }
    
    //cout << "Found widget '" << objname << "' in the DOM!" << endl;
    
    for( auto drf : drfs )
    {
      if( drf )
        m_responses.push_back( drf );
    }
    
    if( m_responses.empty() )
    {
      m_detectorSelect->addItem( "<no responses in directory>" );
      m_detectorSelect->disable();
    }else
    {
      m_detectorSelect->addItem( "<select detector>" );
      m_detectorSelect->enable();
    }
    
    for( const auto drf : m_responses )
      m_detectorSelect->addItem( drf->name() );
    m_detectorSelect->setCurrentIndex( 0 );
    
    if( m_responses.empty() )
    {
      m_msg->setText( "<span style=\"color:red;\">No valid DRFs in directory, or its subdirectories.</span>" );
      m_msg->show();
    }else
    {
      m_detectorSelect->enable();
    }
    
    this->enable();
    
    
    wApp->triggerUpdate();
  };//updategui lambda
  
  auto searchpaths = [basedir, objname, sessid, updategui](){
    const vector<string> dirs = recursive_list_gadras_drfs( basedir );
    std::vector<std::shared_ptr<DetectorPeakResponse> > dets( dirs.size(), nullptr );
    
    SpecUtilsAsync::ThreadPool pool;
    for( size_t i = 0; i < dirs.size(); ++i )
      pool.post( [i,&dets,&dirs](){ dets[i] = GadrasDirectory::parseDetector( dirs[i] ); } );
    pool.join();
    
    //for( const auto &d : dirs )
    //{
    //  auto drf = parseDetector( d );
    //  if( drf )
    //    m_responses.push_back( drf );
    //}
    //end section that should be done off the main thread
    
    Wt::WServer *server = Wt::WServer::instance();
    if( server )
      server->post(sessid, std::bind( [updategui,dets](){ updategui(dets); } ) );
  };//searchpaths lamda
  
  
  Wt::WServer *server = Wt::WServer::instance();
  Wt::WIOService &io = server->ioService();
  io.post( searchpaths );
}//void initDetectors()


void GadrasDirectory::detectorSelected( const int index )
{
  std::shared_ptr<DetectorPeakResponse> det;
  if( index > 0 )
    det = m_responses.at( index - 1 );
  
  m_detectorEdit->setDetector( det );
  m_detectorEdit->emitChangedSignal();
}//void detectorSelected( const int index )





DetectorDisplay::DetectorDisplay( InterSpec *specViewer,
                                  SpectraFileModel *fileModel,
                                  WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_text( NULL ),
    m_interspec( specViewer ),
    m_fileModel( fileModel )
{
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  HelpSystem::attachToolTipOn( this, "The currently selected detector.  Click to select a"
                              " different detector, or modify this one.", showToolTips );
  
  addStyleClass( "DetectorDisplay" );  //In InterSpec.css since this widget is loaded almost always at initial load time anyway

  new WImage( "InterSpec_resources/images/detector_small_white.png", this );
  new WText( "Detector:&nbsp;", this );
  const char *txt = (m_interspec && m_interspec->isMobile()) ? sm_noDetectorTxtMbl : sm_noDetectorTxt;
  m_text = new WText( txt, XHTMLUnsafeText, this );

  std::shared_ptr<DetectorPeakResponse> detector;
  auto meas = specViewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( meas )
    detector = meas->detector();
  m_currentDetector = detector;

  if( detector && detector->isValid() )
    m_text->setText( detector->name() );

  specViewer->detectorChanged().connect( this, &DetectorDisplay::setDetector );
  specViewer->detectorModified().connect( this, &DetectorDisplay::setDetector );

  clicked().connect( this, &DetectorDisplay::editDetector );
}//DetectorDisplay constructor


DetectorDisplay::~DetectorDisplay()
{
} //~DetectorDisplay()


std::shared_ptr<DetectorPeakResponse> DetectorDisplay::detector()
{
  return m_currentDetector.lock();
}// std::shared_ptr<DetectorPeakResponse> DetectorDisplay::detector()


void DetectorDisplay::setDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  m_currentDetector = det;

  if( det && det->isValid() )
  {
    m_text->setText( det->name() );
  }else
  {
    const char *txt = (m_interspec && m_interspec->isMobile()) ? sm_noDetectorTxtMbl : sm_noDetectorTxt;
    m_text->setText( txt );
  }
}//void setDetector( std::shared_ptr<DetectorPeakResponse> det )

void DetectorDisplay::editDetector()
{
  new DetectorEditWindow( m_currentDetector.lock(), m_interspec, m_fileModel );
}//void editDetector()

std::shared_ptr<DetectorPeakResponse> DetectorEdit::detector()
{
  return m_detector;
} //std::shared_ptr<DetectorPeakResponse> DetectorEdit::detector()

DetectorEdit::DetectorEdit( std::shared_ptr<DetectorPeakResponse> currentDet,
                InterSpec *specViewer,
                SpectraFileModel *fileModel,
                AuxWindow* auxWindow )
  : WContainerWidget(),
    m_footer( nullptr ),
    m_interspec( specViewer ),
    m_fileModel( fileModel ),
    m_chart( 0 ),
    m_efficiencyModel( 0 ),
    m_detector( currentDet ),
    m_tabs( nullptr ),
    m_detectorDiameter( nullptr ),
    m_detectrDiameterDiv( nullptr ),
    m_efficiencyCsvUpload( nullptr ),
    m_detectrDotDatDiv( nullptr ),
    m_detectorDotDatUpload( nullptr ),
    m_acceptButton( nullptr ),
    m_cancelButton( nullptr ),
    m_noDrfButton( nullptr ),
    m_detectorManualFunctionName( nullptr ),
    m_detectorManualFunctionText( nullptr ),
    m_detectorManualDescription( nullptr ),
    m_eqnEnergyGroup( nullptr ),
    m_absOrIntrinsicGroup( nullptr ),
    m_detectorManualDiameterText( nullptr ),
    m_detectorManualDistText( nullptr ),
    m_detectorManualDistLabel( nullptr ),
    m_manualSetButton( nullptr ),
    m_gadrasDetSelect( nullptr ),
    m_relEffSelect( nullptr ),
    m_drfTypeMenu( nullptr ),
    m_drfTypeStack( nullptr ),
    m_deleteButton( nullptr ),
    m_DBtable( nullptr ),
    m_model( nullptr ),
    m_defaultForSerialNumber( nullptr ),
    m_defaultForDetectorModel( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/DetectorEdit.css" );
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  setOverflow(Wt::WContainerWidget::OverflowAuto);
  WGridLayout *mainLayout = new WGridLayout();
  
  setLayout( mainLayout );
  mainLayout->setContentsMargins(0, 0, 0, 0);
  mainLayout->setVerticalSpacing( 0 );
  mainLayout->setHorizontalSpacing( 0 );
  
  specViewer->detectorChanged().connect( this, &DetectorEdit::setDetector );
  specViewer->detectorModified().connect( this, &DetectorEdit::setDetector );

  if( currentDet )
  {
    m_originalDetectorCopy.reset( new DetectorPeakResponse() );
    (*m_originalDetectorCopy) = (*currentDet);
    m_originalDetector = currentDet;
    m_previousDetectorDef = *currentDet;
  }//if( currentDet )

  m_previousEmmittedDetector = m_detector;
  
  m_chart = new Wt::Chart::WCartesianChart();
  m_chart->setBackground(Wt::WColor(220, 220, 220));
  m_chart->setXSeriesColumn(0);
  m_chart->setLegendEnabled(false);
  m_chart->setType(Wt::Chart::ScatterPlot);

  //Before, m_efficiencyModel was actually a memory leak, because it didnt have
  //  a parent (m_chart in this case), so everytime you created a new one, the
  //  old one would be leaked since there was no longer a refernce to it anywhere
  m_efficiencyModel = new WStandardItemModel( m_chart );
  m_chart->setModel( m_efficiencyModel );
  
  m_chartEnergyLineColor = WColor("#B02B2C"); //Ruby on Rails red
  m_chartFwhmLineColor = WColor("#3F4C6B");  //Mozilla blue
  
  //We should check the color theme for colors
  auto theme = specViewer->getColorTheme();
  if( theme )
  {
    if( !theme->foregroundLine.isDefault() )
      m_chartEnergyLineColor = theme->foregroundLine;
    if( !theme->backgroundLine.isDefault() )
      m_chartFwhmLineColor = theme->backgroundLine;
    
    if( !theme->spectrumChartText.isDefault() )
    {
      WPen txtpen(theme->spectrumChartText);
      m_chart->setTextPen( txtpen );
      m_chart->axis(Chart::XAxis).setTextPen( txtpen );
      m_chart->axis(Chart::YAxis).setTextPen( txtpen );
      m_chart->axis(Chart::Y2Axis).setTextPen( txtpen );
    }
    
    if( theme->spectrumChartBackground.isDefault() )
      m_chart->setBackground( Wt::NoBrush );
    else
      m_chart->setBackground( WBrush(theme->spectrumChartBackground) );
    
    //From what I can tell, we cant change the legend text color easily, so
    //  we'll just cheat and back the legend background different enough from
    //  black so we can always read the text.  Lame, but whatever.
    m_chart->setLegendStyle( m_chart->legendFont(), m_chart->legendBorder(), WBrush(Wt::WColor(220, 220, 220, 120)) );
    
    if( (theme->spectrumChartMargins.isDefault() && !theme->spectrumChartBackground.isDefault()) )
    {
      //theme->spectrumChartBackground
    }else if( !theme->spectrumChartMargins.isDefault() )
    {
      //theme->spectrumChartMargins
    }
    
    if( !theme->spectrumAxisLines.isDefault() )
    {
      WPen defpen = m_chart->axis(Chart::XAxis).pen();
      defpen.setColor( theme->spectrumAxisLines );
      m_chart->axis(Chart::XAxis).setPen( defpen );
      m_chart->axis(Chart::Y1Axis).setPen( defpen );
      m_chart->axis(Chart::Y2Axis).setPen( defpen );
    }
  }//if( theme )
  
  /*
   * Provide ample space for the title, the X and Y axis and the legend.
   */
  m_chart->setPlotAreaPadding(5, Wt::Top);
  m_chart->setPlotAreaPadding(60, Wt::Bottom);
  m_chart->setPlotAreaPadding(55, Wt::Right | Wt::Left);

  m_chart->setMinimumSize(WLength(350), WLength(200));

  m_chart->axis(Wt::Chart::XAxis).setVisible(true);
  m_chart->axis(Wt::Chart::XAxis).setTitle("Energy (keV)");

  m_chart->axis(Wt::Chart::Y1Axis).setVisible(true);
  m_chart->axis(Wt::Chart::Y1Axis).setTitle("Efficiency");

#if( WT_VERSION >= 0x3030400 )
  m_chart->axis(Wt::Chart::Y1Axis).setTitleOrientation( Wt::Vertical );
  m_chart->axis(Wt::Chart::Y2Axis).setTitleOrientation( Wt::Vertical );
#endif

  m_chart->setLegendEnabled( true );

  m_chart->setLegendLocation(Wt::Chart::LegendInside, Wt::Top, Wt::AlignRight);
  //m_chart->setLegendLocation(Wt::Chart::LegendOutside, Wt::Bottom, Wt::AlignRight);
  
  //m_chart->setTitle("Detector Energy");
  mainLayout->addWidget( m_chart, 0, 0 );
  
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  distValidator->setInvalidBlankText( "0.0 cm" );
  
  WContainerWidget *lowerContent = new WContainerWidget();
  mainLayout->addWidget( lowerContent, 1, 0 );
  
  WGridLayout *lowerLayout = new WGridLayout();
  lowerContent->setLayout( lowerLayout );
  
  //Presize the lower content to accomidate the tallest DRF type ("Formula"),
  //  so everything wont change size when the user selects the different types.
  //  Although, as it stands now if you add a bunch of paths to GADRAS/Rel. Eff.
  //  you can make its content larger than the 190 px.
  //ToDo: improve the sizing of this layout to not have hardcoded sizes!
  lowerContent->resize( WLength::Auto, WLength(190,WLength::Pixel) );
  
  m_drfTypeStack = new Wt::WStackedWidget();
  m_drfTypeStack->addStyleClass( "UseInfoStack DetEditContent" );
  
  m_drfTypeMenu = new WMenu( m_drfTypeStack, Wt::Vertical );
  m_drfTypeMenu->addStyleClass( "VerticalMenu SideMenu DetEditMenu" );
  WContainerWidget *menuHolder = new WContainerWidget();
  menuHolder->addWidget( m_drfTypeMenu );
  lowerLayout->addWidget( menuHolder, 0, 0 );
  lowerLayout->addWidget( m_drfTypeStack, 0, 1 );
  
  mainLayout->setRowStretch( 0, 1 );
  lowerLayout->setColumnStretch( 1, 1 );
  
  //-------------------------------------
  //--- 1)  GADRAS
  //-------------------------------------
  
  m_gadrasDetSelect = new GadrasDetSelect( m_interspec, this );
  
  //-------------------------------------
  //--- 2)  Relative Efficiencies
  //-------------------------------------

  m_relEffSelect = new RelEffDetSelect( m_interspec, this );

  //-------------------------------------
  //--- 3)  Upload
  //-------------------------------------

  
  WContainerWidget *uploadDetTab = new WContainerWidget();
  
  const char *descrstr =
  "<div style=\"white-space:nowrap;\">Intrinsic Efficiency CSV File:</div>"
  "<div>The first column is energy in keV."
  " Second column is the % efficient (0 through 100) for gammas at that"
  " energy, incident on the detector, to be recorded in the photopeak.</div>"
  "<div>You will also need to upload a GADRAS Detector.dat file, or enter the"
  " detectors diameter, after uploading the efficiency file.</div>"
  "<div>You can also upload CSV exported from the <em>Make Detector Response</em> tool.</div>";
  
  WText *descrip = new WText( descrstr, uploadDetTab );
  descrip->setInline( false );
  descrip->setStyleClass("DetectorLabel");
  descrip->setInline( false );
  
  
  m_efficiencyCsvUpload = new WFileUpload( uploadDetTab );
  m_efficiencyCsvUpload->setInline( false );
  m_efficiencyCsvUpload->uploaded().connect( boost::bind( &DetectorEdit::fileUploadedCallback, this, UploadCallbackReason::EfficiencyCsvUploaded ) );
  m_efficiencyCsvUpload->fileTooLarge().connect( boost::bind(&SpecMeasManager::fileTooLarge, _1) );
  m_efficiencyCsvUpload->changed().connect( m_efficiencyCsvUpload, &WFileUpload::upload );
  m_efficiencyCsvUpload->setInline( false );
 
  m_detectrDiameterDiv = new WContainerWidget( uploadDetTab );
  m_detectrDiameterDiv->addStyleClass( "DetectorDiamDiv" );
  WLabel *label = new WLabel( "Enter detector diameter or upload Detector.dat file:", m_detectrDiameterDiv );
  label->setInline( false );

  new WLabel( "Detector Diameter:", m_detectrDiameterDiv );
  m_detectorDiameter = new WLineEdit( "0 cm", m_detectrDiameterDiv );
  
  m_detectorDiameter->setValidator( distValidator );
  m_detectorDiameter->setTextSize( 10 );
  m_detectorDiameter->changed().connect( boost::bind( &DetectorEdit::fileUploadedCallback, this, UploadCallbackReason::DetectorDiameterChanged ) );
  m_detectorDiameter->enterPressed().connect( boost::bind( &DetectorEdit::fileUploadedCallback, this, UploadCallbackReason::DetectorDiameterChanged ) );
  m_detectorDiameter->blurred().connect( boost::bind( &DetectorEdit::fileUploadedCallback, this, UploadCallbackReason::DetectorDiameterChanged ) );

  m_detectrDotDatDiv = new WContainerWidget( m_detectrDiameterDiv );
  m_detectrDotDatDiv->addStyleClass( "DetectorDotDatDiv" );

  label = new WLabel( "GADRAS Detector.Dat File:", m_detectrDotDatDiv );
  label->setInline( false );
  m_detectorDotDatUpload = new WFileUpload( m_detectrDotDatDiv );
  m_detectorDotDatUpload->setInline( false );
  m_detectorDotDatUpload->uploaded().connect( boost::bind( &DetectorEdit::fileUploadedCallback, this, UploadCallbackReason::DetectorDotDatUploaded ) );
  m_detectorDotDatUpload->fileTooLarge().connect( boost::bind(&SpecMeasManager::fileTooLarge, _1) );
  m_detectorDotDatUpload->changed().connect( m_detectorDotDatUpload, &WFileUpload::upload );
  m_detectrDiameterDiv->hide();
  m_detectrDiameterDiv->setHiddenKeepsGeometry( true );

    
  //-------------------------------------
  //--- 4)  Manual
  //-------------------------------------

  const char * const fcn = "exp(-343.63 + 269.10*log(x) -83.80*log(x)^2 "
                          "+ 13.00*log(x)^3 -1.01*log(x)^4 + 0.03*log(x)^5)";
  const char *const manualdetname = "User Defined Detector";
  const char *const diamtxt = "2.2 cm";
  const char *txt = "Manually define the efficiency of the detector,"
                    " given the energy x.  The efficiency values should stay"
                    " within 0.0~1.0 for the energy range of interest."
                    " Click the question mark icon for more info.";
    
  WContainerWidget *formulaDiv = new WContainerWidget();
  WText *selfDefineLabel = new WText( txt, formulaDiv );
  selfDefineLabel->setInline( false );
  selfDefineLabel->setStyleClass("DetectorLabel");
  
  WTable *formulaTable = new WTable( formulaDiv );
  formulaTable->addStyleClass( "FormulaDrfTbl" );
  WTableCell *cell = formulaTable->elementAt( 0, 0 );
  new WLabel("Detector name ", cell );
  cell = formulaTable->elementAt( 0, 1 );
  m_detectorManualFunctionName = new WLineEdit( cell );
  m_detectorManualFunctionName->setStyleClass("DetectorEditFunctionalFormText");
  m_detectorManualFunctionName->setText( manualdetname );
  m_detectorManualFunctionName->setEmptyText("Unique Detector Name");
  m_detectorManualFunctionName->setInline(false);
  m_detectorManualFunctionName->keyWentUp().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));

  cell = formulaTable->elementAt( 1, 0 );
  new WLabel("Efficiency f(x) = ", cell );
  cell = formulaTable->elementAt( 1, 1 );
  cell->setRowSpan( 2 );
  m_detectorManualFunctionText = new WTextArea( cell );
  //m_detectorManualFunctionText->setColumns(60);
  
  m_detectorManualFunctionText->setRows(3);
  m_detectorManualFunctionText->setStyleClass("DetectorEditFunctionalFormText");
  m_detectorManualFunctionText->setText(fcn);
  m_detectorManualFunctionText->setEmptyText(fcn);
  m_detectorManualFunctionText->setInline(false);
  m_detectorManualFunctionText->keyWentUp().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  
  cell = formulaTable->elementAt( 2, 0 );
  Wt::WContainerWidget *energyContainer = new Wt::WContainerWidget( cell );
  energyContainer->setMargin( 6, Wt::Top );
  m_eqnEnergyGroup = new WButtonGroup( energyContainer );
  WRadioButton *button = new WRadioButton( "keV", energyContainer );
  m_eqnEnergyGroup->addButton( button, 0 );
  button = new WRadioButton( "MeV", energyContainer );
  button->setMargin( 5, Wt::Left );
  m_eqnEnergyGroup->addButton( button, 1 );
  m_eqnEnergyGroup->checkedChanged().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  m_eqnEnergyGroup->setSelectedButtonIndex( 0 );
    
  HelpSystem::attachToolTipOn( energyContainer,"Energy unit for efficiency formula" ,
                              showToolTips, HelpSystem::ToolTipPosition::Top );

  cell = formulaTable->elementAt( 3, 0 );
  new WLabel("Description", cell );
  cell = formulaTable->elementAt( 3, 1 );
  m_detectorManualDescription = new WLineEdit( cell );
  m_detectorManualDescription->setWidth( WLength(100,WLength::Percentage) );
  m_detectorManualDescription->setEmptyText( "Description of this detector response" );
  m_detectorManualDescription->changed().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  
  cell = formulaTable->elementAt( 4, 0 );
  new WLabel("Detector diameter ", cell );
  cell = formulaTable->elementAt( 4, 1 );
  m_detectorManualDiameterText = new WLineEdit( cell );
  m_detectorManualDiameterText->setValidator( distValidator );
  m_detectorManualDiameterText->setText( diamtxt );
  m_detectorManualDiameterText->setEmptyText( diamtxt );
  m_detectorManualDiameterText->setValidator( distValidator );
  m_detectorManualDiameterText->blurred().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  m_detectorManualDiameterText->enterPressed().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  
  cell = formulaTable->elementAt( 5, 0 );
  cell->setColumnSpan( 2 );
  m_absOrIntrinsicGroup = new WButtonGroup( cell );
  button = new WRadioButton( "Intrinsic", cell );
  m_absOrIntrinsicGroup->addButton( button, 0 );
  button = new WRadioButton( "Absolute", cell );
  button->setMargin( 5, Wt::Left );
  m_absOrIntrinsicGroup->addButton( button, 1 );
  m_absOrIntrinsicGroup->checkedChanged().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  m_absOrIntrinsicGroup->setSelectedButtonIndex( 0 );

  m_detectorManualDistLabel = new WLabel( "Dist:", cell );
  m_detectorManualDistLabel->setMargin( 10, Wt::Left );
  
  m_detectorManualDistText = new Wt::WLineEdit( cell );
  m_detectorManualDistText->setValidator( distValidator );
  m_detectorManualDistText->setEmptyText("1 m");
  m_detectorManualDistText->setHiddenKeepsGeometry( true );
  m_detectorManualDistLabel->hide();
  m_detectorManualDistText->hide();
  m_detectorManualDistText->keyWentUp().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  
  m_manualSetButton = new WPushButton( "Set", cell );
  //m_manualSetButton->setInline(false);
  m_manualSetButton->setFloatSide( Wt::Right );
  //selfDefineLayout->setColumnStretch( 1, 1 );
  m_manualSetButton->clicked().connect(boost::bind(&DetectorEdit::setDefineDetector, this));
    
  //-------------------------------------
  //--- 5)  Recent
  //-------------------------------------
  WContainerWidget *recentDiv = new WContainerWidget( );
  WGridLayout* recentDivLayout = new WGridLayout( recentDiv );
  recentDivLayout->setContentsMargins( 0, 0, 0, 0 );

  Dbo::ptr<InterSpecUser> user = m_interspec->m_user;
  
  m_DBtable = new Wt::WTreeView();
  m_DBtable->setRootIsDecorated	(	false ); //makes the tree look like a table! :)
  
  m_DBtable->addStyleClass( "DbSpecFileSelectTable" );
  
  //We have to create a independant Dbo::Session for this class since the
  //  m_viewer->m_user.session() could be used in other threads, messing
  //  things up (Dbo::Session is not thread safe).
  m_sql.reset( new DataBaseUtils::DbSession( *m_interspec->sql() ) );

  m_model = new Dbo::QueryModel< Dbo::ptr<DetectorPeakResponse> >( m_DBtable );

  const auto userid = user.id();
  m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                           .where( "InterSpecUser_id = ? ").bind( userid )
                           .orderBy( "-1*m_lastUsedUtc" )
                    );
  
  /*  ToDo: - add column for when was last used, but give time in how long ago (remove diameter column)
            - Add bullets or dropdown to allow filtering by "All", "Uploaded", "Entered Formula", "Made From Data"
   */
  
  m_model->addColumn( "m_name" );
  m_model->addColumn( "m_description" );
  m_model->addColumn( "m_lastUsedUtc" );
  
  
  m_model->setHeaderData(  0, Horizontal, WString("Name"), DisplayRole );
  m_DBtable->setColumnWidth( 0, 200 );
  m_model->setHeaderData(  1, Horizontal, WString("Description"), DisplayRole );
  m_DBtable->setColumnWidth( 1, 150 );
  m_model->setHeaderData(  2, Horizontal, WString("Last Used"), DisplayRole );
  m_DBtable->setColumnWidth( 2, 150 );
 
  
  WItemDelegate *delegate = new UtcToLocalTimeDelegate( m_DBtable );
  m_DBtable->setItemDelegateForColumn( 2, delegate );

  m_DBtable->setModel( m_model );
  m_DBtable->setAlternatingRowColors( true );
  m_DBtable->setSelectionMode( SingleSelection );
  
  for( int col = 0; col < m_model->columnCount(); ++col )
  {
    m_DBtable->setColumnHidden( col, false );
    m_DBtable->setSortingEnabled( col, true );
  }//for( int col = 0; col < model->columnCount(); ++col )
  
  m_DBtable->selectionChanged().connect( this, &DetectorEdit::dbTableSelectionChanged );

  m_deleteButton = new WPushButton( "Delete" );
  m_deleteButton->setIcon( "InterSpec_resources/images/minus_min_white.png" );
  m_deleteButton->setDefault(true);
  m_deleteButton->disable();
  m_deleteButton->clicked().connect( this, &DetectorEdit::deleteDBTableSelected );

  WComboBox *filter = new WComboBox();
  filter->addItem( "All" );
  filter->addItem( "Uploaded" );
  filter->addItem( "Entered Formula" );
  filter->addItem( "Made From Data" );
  filter->addItem( "Included In N42 File" );
  filter->setCurrentIndex( 0 );
  
  filter->changed().connect( std::bind( [filter,this,userid](){
    m_DBtable->setSelectedIndexes( WModelIndexSet() );
    m_deleteButton->disable();
    
    if( filter->currentIndex() == 0 ) //"All"
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? ").bind( userid ), true );
    }else if( filter->currentIndex() == 1 ) //"Uploaded"
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND (m_efficiencySource = ? OR m_efficiencySource = ?)")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::UserImportedIntrisicEfficiencyDrf )
                        .bind( DetectorPeakResponse::DrfSource::UserImportedGadrasDrf ) , true );
    }else if( filter->currentIndex() == 2 ) //"Entered Formula"
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND m_efficiencySource = ?")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::UserSpecifiedFormulaDrf ) , true );
    }else if( filter->currentIndex() == 3 )
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND m_efficiencySource = ?")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::UserCreatedDrf ) , true );
    }else if( filter->currentIndex() == 4 )
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND m_efficiencySource = ?")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::FromSpectrumFileDrf ) , true );
    }else
    {
      return; //shouldnt happen
    }
    
    m_DBtable->sortByColumn( 2, Wt::SortOrder::DescendingOrder );
  } ) );
  
  recentDivLayout->addWidget( m_DBtable, 0, 0, 1, 3 );
  recentDivLayout->addWidget( m_deleteButton, 1, 0, AlignLeft);
  
  label = new WLabel( "Show:" );
  label->setBuddy( filter );
  recentDivLayout->addWidget( label, 1, 1, AlignRight | AlignMiddle );
  
  recentDivLayout->addWidget( filter, 1, 2 );
  recentDivLayout->setRowStretch( 0, 1 );
  //recentDiv->setOverflow(Wt::WContainerWidget::OverflowHidden);
  //recentDiv->setMaximumSize( WLength::Auto, 180 );
  recentDiv->setHeight( WLength( 100, WLength::Unit::Percentage ) );
  m_DBtable->sortByColumn( 2, Wt::SortOrder::DescendingOrder );
  
  //--------------------------------------------------------------------------------
  
  WMenuItem *item = m_drfTypeMenu->addItem( "GADRAS", m_gadrasDetSelect );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  item = m_drfTypeMenu->addItem( "Rel. Eff.", m_relEffSelect );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  item = m_drfTypeMenu->addItem( "Import", uploadDetTab );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  if( !m_interspec->isSupportFile() )
    item->setHidden(true);
  
  item = m_drfTypeMenu->addItem( "Formula", formulaDiv );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  item = m_drfTypeMenu->addItem( "Previous", recentDiv );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  m_drfTypeStack->setHeight( WLength(185.0) );
  
  
  auto meas = specViewer->measurment( SpecUtils::SpectrumType::Foreground );
  
  if( meas && !meas->instrument_id().empty() )
  {
    m_defaultForSerialNumber = new WCheckBox( "Use as default DRF for SN '" + meas->instrument_id() + "'" );
    m_defaultForSerialNumber->setMargin( 21, Wt::Left );
    mainLayout->addWidget( m_defaultForSerialNumber, mainLayout->rowCount(), 0, 1, mainLayout->columnCount() );
  }//if( have serial number )
  
  if( meas && (meas->detector_type()!=DetectorType::Unknown || !meas->instrument_model().empty()) )
  {
    string model;
    if( meas->detector_type() == DetectorType::Unknown )
      model = meas->instrument_model();
    else
      model = detectorTypeToString( meas->detector_type() );
    m_defaultForDetectorModel = new WCheckBox( "Use as default DRF for model '" + model + "'" );
    m_defaultForDetectorModel->setMargin( 21, Wt::Left );
    mainLayout->addWidget( m_defaultForDetectorModel, mainLayout->rowCount(), 0, 1, mainLayout->columnCount() );
  }//if( have
  
  
  if( auxWindow )
  {
    m_footer = auxWindow->footer();
  }else
  {
    m_footer = new WContainerWidget();
    mainLayout->addWidget( m_footer, mainLayout->rowCount(), 0, 1, 2 );
  }
  
  AuxWindow::addHelpInFooter( m_footer, "detector-edit-dialog" );
  
  
  m_noDrfButton = new WPushButton( "No Detector", m_footer );
  m_noDrfButton->clicked().connect( this, &DetectorEdit::finishWithNoDetector );
  if( specViewer && !specViewer->isPhone() )
    m_noDrfButton->addStyleClass( "NoDrfBtn" );
  
  if( auxWindow )
  {
    m_cancelButton = auxWindow->addCloseButtonToFooter( "Cancel" );
  }else
  {
    WPushButton *close = new WPushButton( "Close" );
    m_footer->insertWidget( m_footer->count(), close );
  }
  
  m_cancelButton->clicked().connect( this, &DetectorEdit::cancelAndFinish );

  
  if( specViewer && !specViewer->isPhone() )
  {
    m_acceptButton = new WPushButton( "Accept", m_footer );
    m_acceptButton->setFloatSide(Wt::Right);
    
    m_cancelButton->setIcon( "InterSpec_resources/images/reject.png" );
    HelpSystem::attachToolTipOn( m_cancelButton, "Remove all changes or selections made by this"
                                 " dialog, and close the dialog" , showToolTips );
    
    m_acceptButton->setIcon( "InterSpec_resources/images/accept.png" );
    HelpSystem::attachToolTipOn( m_acceptButton, "Accept all changes/selections made and close dialog" , showToolTips );
  }else
  {
    m_acceptButton = new WPushButton( "Use Detector", m_footer );
    m_acceptButton->addStyleClass( "CenterBtnInMblAuxWindowHeader" );
  }//if( isMobile() ) / else
  
  m_acceptButton->clicked().connect( this, &DetectorEdit::acceptAndFinish );
  if( !currentDet || !currentDet->isValid() )
    setAcceptButtonEnabled( false );
  
  m_drfTypeMenu->select( 0 );
  
  init();
}//DetectorEdit constructor


void DetectorEdit::setAcceptButtonEnabled( const bool enable )
{
  m_acceptButton->setEnabled( enable );
  m_noDrfButton->setHidden( !m_detector );
  if( m_defaultForSerialNumber )
    m_defaultForSerialNumber->setEnabled( enable );
  if( m_defaultForDetectorModel )
    m_defaultForDetectorModel->setEnabled( enable );
}


void DetectorEdit::deleteDBTableSelected()
{
  WModelIndexSet indices = m_DBtable->selectedIndexes();
  if( !indices.size() )
  {
    m_deleteButton->disable();
    return;
  }//if( !indices.size() )
  
  WModelIndex index = *indices.begin();
  DataBaseUtils::DbTransaction transaction( *m_sql );
  Dbo::ptr<DetectorPeakResponse> dbfile = m_model->stableResultRow( index.row() );
  dbfile.remove();
  transaction.commit();
  m_model->reload();
} //deleteDBTableSelected

// Grabs the detector if selected in DB table
void DetectorEdit::dbTableSelectionChanged()
{
  if( m_DBtable->selectedIndexes().size() )
    m_deleteButton->enable();
  
  bool failed = false;
  std::shared_ptr<DetectorPeakResponse> det;

  try
  {
    if( !m_model  )
      return;
    
    const WModelIndexSet indices = m_DBtable->selectedIndexes();
    
    if( indices.empty() )
    {
      m_deleteButton->disable();
      return;
    }//if( !indices.empty() )
    
    m_deleteButton->enable();
    const WModelIndex index = *indices.begin();
    Dbo::ptr<DetectorPeakResponse> dbfile = m_model->stableResultRow( index.row() );
    det = std::shared_ptr<DetectorPeakResponse>(new DetectorPeakResponse(*dbfile));
  }catch(...)
  {
    failed = true;
    passMessage( "Error getting from DetectorPeakResponse table", "DetectorEdit", WarningWidget::WarningMsgHigh );
  }//try / catch
  
  m_detector = det;
  setAcceptButtonEnabled( !failed );
  emitChangedSignal();
}//void dbTableSelectionChanged()

void DetectorEdit::init()
{
  //----initialize
  
  if( !m_detector )
  {
    //m_drfTypeMenu->select( 0 );
  }else
  {
    switch( m_detector->drfSource() )
    {
      case DetectorPeakResponse::DefaultGadrasDrf:
      case DetectorPeakResponse::UserAddedGadrasDrf:
      {
        const bool selected = m_gadrasDetSelect->trySelectDetector( m_detector );
        
        if( selected )
          m_drfTypeMenu->select( 0 );
        else
          m_drfTypeMenu->select( 2 );
        break;
      }
      
      case DetectorPeakResponse::UserAddedRelativeEfficiencyDrf:
      {
        m_drfTypeMenu->select( 1 );
        m_relEffSelect->trySelectDetector( m_detector );
        break;
      }
        
      case DetectorPeakResponse::UserSpecifiedFormulaDrf:
      {
        m_drfTypeMenu->select( 3 );
        m_detectorManualFunctionName->setText( m_detector->name() );
        m_detectorManualFunctionText->setText( m_detector->efficiencyFormula() );
        m_detectorManualDescription->setText( m_detector->description() );
        const float diam = m_detector->detectorDiameter();
        const string diamstr = PhysicalUnits::printToBestLengthUnits( diam );
        m_detectorManualDiameterText->setText( diamstr );
        const float units = m_detector->efficiencyEnergyUnits();
        int energygrp = 0;
        if( fabs(PhysicalUnits::MeV-units) < fabs(PhysicalUnits::keV-units) )
          energygrp = 1;
        m_eqnEnergyGroup->setSelectedButtonIndex( energygrp );
        m_absOrIntrinsicGroup->setSelectedButtonIndex( 0 );
        break;
      }//case kUserEfficiencyEquationSpecified:
        
      case DetectorPeakResponse::UserImportedIntrisicEfficiencyDrf:
      {
        m_drfTypeMenu->select( 2 );
        break;
      }
      
      case DetectorPeakResponse::UnknownDrfSource:
      case DetectorPeakResponse::UserImportedGadrasDrf:
      case DetectorPeakResponse::UserCreatedDrf:
      case DetectorPeakResponse::FromSpectrumFileDrf:
      {
        m_drfTypeMenu->select( 4 );
        try
        {
          const double start_time = SpecUtils::get_wall_time();
          for( int row = 0; row < m_model->rowCount(); ++row )
          {
            Wt::Dbo::ptr<DetectorPeakResponse> drf = m_model->resultRow(row);
            if( drf && (drf->hashValue() == m_detector->hashValue()) )
            {
              WModelIndexSet indexset;
              indexset.insert( m_model->index(row,0) );
              m_DBtable->setSelectedIndexes( indexset );
              break;
            }
            
            //On a quick test, it took about 25ms when the very first row
            //  matched the detector; even if we loop through all rows, it
            //  takes about the same amount of time
            const double now_time = SpecUtils::get_wall_time();
            if( (now_time-start_time) > 0.25 )
              break;
          }//for( loop over previous DRFs )
        }catch( std::exception &e )
        {
          cerr << "Caught exception trying to select previous DRF from database to select for the GUI: " << e.what() << endl;
        }//try / catch
        break;
      }
    };//switch( m_detector->drfSource() )
  }//if( m_detector == NULL ) / else
  
  selectButton( m_drfTypeStack, m_drfTypeMenu, false );
}// init()


DetectorEdit::~DetectorEdit()
{
} //~DetectorEdit()


void DetectorEdit::setDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  if( m_detector != det )
  {
    m_detector = det;
    init();
  }
}//void setDetector( std::shared_ptr<DetectorPeakResponse> det )


void DetectorEdit::verifyManualDefinition()
{ //if we change anything in the text fields, we should disable Accept
  setAcceptButtonEnabled( false );
  m_manualSetButton->enable();
  
  const bool hideAbs = (m_absOrIntrinsicGroup->checkedId() == 0);
  m_detectorManualDistLabel->setHidden( hideAbs );
  m_detectorManualDistText->setHidden( hideAbs );

/*
  try
  {
    string txt = m_detectorManualDiameterText->text().toUTF8();
    const float diam = (float)PhysicalUnits::stringToDistance( txt );
    if( m_absOrIntrinsicGroup->checkedId() != 0 )
      PhysicalUnits::stringToDistance( m_detectorManualDistText->text().toUTF8() );
    txt = m_detectorManualFunctionText->text().toUTF8();
    
    std::shared_ptr<DetectorPeakResponse> detec
                  = std::make_shared<DetectorPeakResponse>(
                                  new DetectorPeakResponse( "temp", "temp ") );
    
    detec->setIntrinsicEfficiencyFormula( txt, diam, 1.0f );
  }catch( std::exception & )
  {
    m_manualSetButton->disable();
  }//try / catch
*/
} //verifyManualDefinition


void DetectorEdit::setDefineDetector()
{
  std::shared_ptr<DetectorPeakResponse> detec = detector();
//  if( !!detec &&
//      detec->efficiencySource()==DetectorPeakResponse::kUserEfficiencyEquationSpecified )
//  {
//    m_absOrIntrinsicGroup->setSelectedButtonIndex( 0 );
//    m_detectorManualDistLabel->hide();
//    m_detectorManualDistText->hide();
//    detec->description();
//    const float diam = detec->detectorDiameter();
//    const string diamstr = PhysicalUnits::printToBestLengthUnits( diam );
//    m_detectorManualFunctionText->setText( detec->efficiencyFormula() );
//    m_detectorManualFunctionName->setText( detec->name() );
//    m_detectorManualDiameterText->setText( diamstr );
//    m_detectorManualDistText->setText( "0" );
//  
//    setAcceptButtonEnabled( true );
//    m_manualSetButton->disable();
//    return;
//  }//if(...)
  
  string fcn = m_detectorManualFunctionText->text().narrow();
  const string name = m_detectorManualFunctionName->text().toUTF8();
  string descr = m_detectorManualDescription->text().toUTF8();
  if( descr.empty() )
  {
    descr = "f(x)=" + fcn;
    if( descr.size() > 25 )
      descr = descr.substr( 0, 25 ) + "...";
    m_detectorManualDescription->setText( descr );
  }//if( descr.empty() )

  detec.reset( new DetectorPeakResponse( name, descr ) );

  try
  {
    const string diamtxt = m_detectorManualDiameterText->text().narrow();
    const float det_diam = (float)PhysicalUnits::stringToDistance( diamtxt );
    if( det_diam < 0.0001f )
      throw runtime_error( "Detector diameters must be positive" );
    
    if( m_absOrIntrinsicGroup->checkedId() != 0 )
    {
      const string disttxt = m_detectorManualDistText->text().narrow();
      const float dist = (float)PhysicalUnits::stringToDistance( disttxt );
      if( dist < 0.0001f )
        throw runtime_error( "Distance must be positive for absolute efficiencies" );
      
      
      const float gfactor = DetectorPeakResponse::fractionalSolidAngle(
                                                              det_diam, dist );
        
      fcn = std::to_string(1.0/gfactor) + "*(" + fcn + ")";
    }//if( m_absOrIntrinsicGroup->groupId() != 0 )
    
    const float energyUnits
              = static_cast<float>(m_eqnEnergyGroup->checkedId()==0
                                    ? PhysicalUnits::keV : PhysicalUnits::MeV);
    
    detec->setIntrinsicEfficiencyFormula( fcn, det_diam, energyUnits, 0.0f, 0.0f );
    detec->setDrfSource( DetectorPeakResponse::DrfSource::UserSpecifiedFormulaDrf );
    
    updateChart();
  }catch( std::exception &e )
  {
    passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    setAcceptButtonEnabled( false );
    m_manualSetButton->enable();
    updateChart();
    return;
  }//try / catch
  
  m_manualSetButton->disable();
  m_detector = detec;
  setAcceptButtonEnabled( true );

  emitChangedSignal();
} // setDefineDetector(const std::string fcn)



// This calls the action to initialize the charts to whatever
// is selected.  Each stack page should have something to update
void DetectorEdit::selectButton( WStackedWidget *stack,
                                 WMenu *menu,
                                 bool activateCallBack)
{
  //makes one of the stack visible
  stack->setCurrentIndex( menu->currentIndex() );  //I think this is now not needed
  
  if( activateCallBack )
  {
    
    //clear the chart when we change, regardless
    std::shared_ptr<DetectorPeakResponse> det;
    m_detector = det;
    
    switch( menu->currentIndex() )
    {
      case 0:
        gadrasDetectorSelectCallback();
        break;
    
      case 1:
        relEffDetectorSelectCallback();
        break;
        
      case 2:
        fileUploadedCallback( UploadCallbackReason::ImportTabChosen );
        break;
        
      case 3:
        setDefineDetector();
        break;
        
      case 4:
        dbTableSelectionChanged();
        break;
    } //switch
  }//if( activateCallBack )
  
  updateChart();
} //DetectorEdit::selectButton(Wt::WStackedWidget * stack, int index)


void DetectorEdit::updateChart()
{
  // clear series if any
  m_chart->removeSeries(1);
  m_chart->removeSeries(2);
  
  const bool hasEfficiency = !!m_detector && m_detector->isValid();
  
  if( hasEfficiency )
  {
    int nValidPoints = 0, nEffPoints = 0;
    try
    {
      const bool hasResloution = m_detector->hasResolutionInfo();
      const float minEnergy = 0.0f, maxEnergy = 3000.0f;
      const int numEnergyPoints = 1200; //whatever, you can pick this better using the chart width or whateever
    
      m_efficiencyModel->clear();
      m_efficiencyModel->insertRows( 0, numEnergyPoints );
      m_efficiencyModel->insertColumns( 0, hasResloution? 3 : 2 );
      float energy=0;
      float efficincy=0;
      for( int row = 0; row < numEnergyPoints; ++row )
      {
        energy = minEnergy + (float(row)/float(numEnergyPoints)) * (maxEnergy-minEnergy);
        efficincy = static_cast<float>( m_detector->intrinsicEfficiency( energy ) );
        m_efficiencyModel->setData(row, 0, energy );


        //Skip any points outside where we would expect.
        if( IsNan(efficincy) || IsInf(efficincy) || efficincy > 1.2f || efficincy < 0.0f )
        {
          m_efficiencyModel->setData(row, 1, boost::any() );
        }else
        {
          ++nValidPoints;
          m_efficiencyModel->setData(row, 1, efficincy );
        }
        
        if( hasResloution )
        {
          const float fwhm = m_detector->peakResolutionFWHM(energy);
          if( !IsNan(fwhm) && !IsInf(fwhm) && fwhm>=0.0f && fwhm<9999.9f )
          {
            ++nEffPoints;
            m_efficiencyModel->setData(row, 2, fwhm);
          }
        }
      }//for( int row = 0; row < numEnergyPoints; ++row )
     
      m_efficiencyModel->setHeaderData(0, Wt::WString("Energy"));
      m_efficiencyModel->setHeaderData(1, Wt::WString("Efficiency"));
      if( m_detector->hasResolutionInfo() )
        m_efficiencyModel->setHeaderData(2, Wt::WString("FWHM"));
    
      m_chart->setXSeriesColumn(0);  //Not having this line after creating a new model was why the x-axis was the wrong scale
      Wt::Chart::WDataSeries s1(1, Wt::Chart::LineSeries,Wt::Chart::Y1Axis);
      s1.setPen( WPen(m_chartEnergyLineColor) );
      m_chart->addSeries(s1);
    
      if( nValidPoints == 0 )
      {
        m_chart->axis(Wt::Chart::Y1Axis).setRange( 0.0, 1.0 );
      }else
      {
        //m_chart->axis(Wt::Chart::Y1Axis).setRoundLimits(<#WFlags<Wt::Chart::AxisValue> locations#>)
        m_chart->axis(Wt::Chart::Y1Axis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
      }
    
      if( nEffPoints )
      { //only if there is resolution FWHM
        Wt::Chart::WDataSeries s2(2, Wt::Chart::LineSeries,Wt::Chart::Y2Axis);
        s2.setPen( WPen(m_chartFwhmLineColor) );
        m_chart->addSeries(s2);
        m_chart->axis(Wt::Chart::Y2Axis).setTitle("FWHM");
        m_chart->axis(Wt::Chart::Y2Axis).setVisible(true);
      }else
      { //no ResolutionInfo
        m_chart->axis(Wt::Chart::Y2Axis).setVisible(false);
        m_chart->axis(Wt::Chart::Y2Axis).setTitle("");
      }//no ResolutionInfo
    }catch( std::exception &e )
    {
      cerr << "DetectorEdit::updateChart()\n\tCaught: " << e.what() << endl;
    }//try / catch
  }//if( hasEfficiency )
} //DetectorEdit::updateChart()


std::shared_ptr<DetectorPeakResponse> DetectorEdit::checkIfFileIsRelEff( const std::string filename )
{
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream csvfile( wfilename.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream csvfile( filename.c_str(), ios_base::binary | ios_base::in );
#endif
  
  if( !csvfile.is_open() )
    return nullptr;
  
  string line;
  int nlineschecked = 0;
  
  string drfname, drfdescrip;
  bool foundMeV = false, foundKeV = false;
  //ToDo: Need to implement getting lines safely where a quoted field may span
  //      several lines.
  while( SpecUtils::safe_get_line(csvfile, line, 2048) && (++nlineschecked < 100) )
  {
    foundKeV |= SpecUtils::icontains( line, "kev" );
    foundMeV |= SpecUtils::icontains( line, "mev" );
    
    vector<string> fields;
    split_escaped_csv( fields, line );
    
    if( fields.size() == 2 )
    {
      if( fields[0] == "# Name" )
        drfname = fields[1];
      if( fields[0] == "# Description" )
        drfdescrip = fields[1];
    }
    
    if( fields.size() < 16 )
      continue;
    
    if( !SpecUtils::iequals_ascii( fields[3], "c0")
       || !SpecUtils::iequals_ascii( fields[4], "c1")
       || !SpecUtils::iequals_ascii( fields[5], "c2")
       || !SpecUtils::iequals_ascii( fields[6], "c3")
       || !SpecUtils::icontains( fields[15], "radius") )
      continue;
    
    //Okay, next line should be
    if( !SpecUtils::safe_get_line(csvfile, line, 2048) )
      throw runtime_error( "Couldnt get next line" );
    
    split_escaped_csv( fields, line );
    
    try
    {
      vector<float> coefs( 8, 0.0f );
      for( int i = 0; i < 8; ++i )
      {
        string field = fields.at(3+i);
        coefs[i] = (field.empty() ? 0.0f : std::stof(field));
      }
      
      for( int i = 7; i >= 0; --i )
      {
        if( coefs[i] != 0.0f )
          break;
        else
          coefs.resize( coefs.size() - 1 );
      }
      
      if( coefs.size() < 1 )
        continue;
      
      const float dist = std::stof( fields.at(14) ) * PhysicalUnits::cm;
      const float radius = std::stof( fields.at(15) ) * PhysicalUnits::cm;
      
      const string name = (fields[0].empty() ? drfname : fields[0]);
      const float energUnits = ((foundKeV && !foundMeV) ? 1.0f : 1000.0f);
      float lowerEnergy = 0.0f, upperEnergy = 0.0f;
      auto det = std::make_shared<DetectorPeakResponse>( fields[0], drfdescrip );
      
      det->fromExpOfLogPowerSeriesAbsEff( coefs, dist, 2.0f*radius, energUnits,
                                          lowerEnergy, upperEnergy );
      
      //Look for the line that gives the appropriate energy range.
#define POS_DECIMAL_REGEX "\\+?\\s*((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?\\s*"
      const char * const rng_exprsn_txt = "Valid energy range:\\s*(" POS_DECIMAL_REGEX ")\\s*keV to\\s*(" POS_DECIMAL_REGEX ")\\s*keV.";
      
      RegexNs::regex range_expression( rng_exprsn_txt );
      
      while( SpecUtils::safe_get_line(csvfile, line, 2048) && (++nlineschecked < 100) )
      {
        if( SpecUtils::icontains( line, "Full width half maximum (FWHM) follows equation" ) )
        {
          const bool isSqrt = SpecUtils::icontains( line, "sqrt(" );  //Else contains "GadrasEqn"
          const auto form = isSqrt ? DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial : DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
          int nlinecheck = 0;
          while( SpecUtils::safe_get_line(csvfile, line, 2048)
                && !SpecUtils::icontains(line, "Values")
                && (++nlinecheck < 15) )
          {
          }
        
          vector<string> fwhm_fields;
          split_escaped_csv( fwhm_fields, line );
          
          if( fwhm_fields.size() > 1 && SpecUtils::icontains(fwhm_fields[0], "Values") )
          {
            try
            {
              vector<float> coefs;
              for( size_t i = 1; i < fwhm_fields.size(); ++i )
                coefs.push_back( stof(fwhm_fields[i]) );
              if( !coefs.empty() )
                det->setFwhmCoefficients( coefs,form );
            }catch(...)
            {
            }
            break;
          }
        }//if( start of FWHM section of CSV file )
        
        RegexNs::smatch range_matches;
        if( RegexNs::regex_search( line, range_matches, range_expression ) )
        {
          lowerEnergy = std::stof( range_matches[1] );
          upperEnergy = std::stof( range_matches[6] );
          det->setEnergyRange( lowerEnergy, upperEnergy );
          break;
        }
      }//while( getline )
      
      return det;
    }catch(...)
    {
      continue;
    }
  }//while( more lines )
  
  return nullptr;
}//checkIfFileIsRelEff(...)



void DetectorEdit::fileUploadedCallback( const UploadCallbackReason context )
{
  if( !m_efficiencyCsvUpload->empty() )
  {
    m_detectrDiameterDiv->show();
  }else
  {
    m_detectrDiameterDiv->hide();
  }

  switch( context )
  {
    case UploadCallbackReason::ImportTabChosen:
      break;
    case UploadCallbackReason::DetectorDiameterChanged:
      break;
    case UploadCallbackReason::DetectorDotDatUploaded:
      break;
      
    case UploadCallbackReason::EfficiencyCsvUploaded:
      if( !m_efficiencyCsvUpload->empty() )
      {
        const string filename = m_efficiencyCsvUpload->spoolFileName();
        auto det = DetectorEdit::checkIfFileIsRelEff( filename );
        
        if( det )
        {
          auto detDiamStr = PhysicalUnits::printToBestLengthUnits( det->detectorDiameter() );
          m_detectorDiameter->setText( detDiamStr );
          m_detectrDotDatDiv->hide();
          m_detectorDiameter->disable();
          m_detector = det;
          setAcceptButtonEnabled( true );
          emitChangedSignal();
          return;
        }else
        {
          m_detectrDotDatDiv->show();
          m_detectorDiameter->enable();
        }//if( det )
      }//if( we have a file to test )
      break;
  }//switch( context )

  const bool isGadrasDet = (!m_efficiencyCsvUpload->empty()
                            && !m_detectorDotDatUpload->empty()
                            && m_efficiencyCsvUpload->spoolFileName().size()
                            && m_detectorDotDatUpload->spoolFileName().size());

  float diameter = -1.0f;
  try
  {
    diameter = static_cast<float>( PhysicalUnits::stringToDistance( m_detectorDiameter->text().narrow() ) );
    if( diameter < float(0.001*PhysicalUnits::cm) )
      diameter = 0.0f;
  }catch(...)
  {
    m_detectorDiameter->setText( "" );
  }

  const bool isDiamDet = (!m_efficiencyCsvUpload->empty()
                          && diameter>0.0
                          && m_efficiencyCsvUpload->spoolFileName().size() );

  if( !isGadrasDet && !isDiamDet )
  {
    if( diameter > 0.0 )
      passMessage( "Need the detectors diameter", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( !isGadrasDet && !isDiamDet )

  const string csvFileName = m_efficiencyCsvUpload->spoolFileName();

#ifdef _WIN32
  const std::wstring wcsvFileName = SpecUtils::convert_from_utf8_to_utf16(csvFileName);
  ifstream csvfile( wcsvFileName.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream csvfile( csvFileName.c_str(), ios_base::binary|ios_base::in );
#endif
  
  if( !csvfile.is_open() )
    return;

  auto det = std::make_shared<DetectorPeakResponse>();
  try
  {
    if( isGadrasDet )
    {
      const string dotDatFileName = m_detectorDotDatUpload->spoolFileName();
#ifdef _WIN32
      const std::wstring wdotDatFileName = SpecUtils::convert_from_utf8_to_utf16(dotDatFileName);
      ifstream datfile( wdotDatFileName.c_str(), ios_base::binary|ios_base::in );
#else
      ifstream datfile( dotDatFileName.c_str(), ios_base::binary|ios_base::in );
#endif
      if( !datfile.is_open() )
        return;
      det->fromGadrasDefinition( csvfile, datfile );
      det->setDrfSource( DetectorPeakResponse::DrfSource::UserImportedGadrasDrf );
      if( context == UploadCallbackReason::DetectorDotDatUploaded )
        m_detectorDiameter->setText( PhysicalUnits::printToBestLengthUnits(det->detectorDiameter()) );
    }else
    {
      det->fromEnergyEfficiencyCsv( csvfile, diameter, float(PhysicalUnits::keV) );
      det->setDrfSource( DetectorPeakResponse::DrfSource::UserImportedIntrisicEfficiencyDrf );
    }//if( isGadrasDet ) / else
  }catch( std::exception &e )
  {
    setAcceptButtonEnabled( false );
    passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    return;
  }//try / catch
  
  string csvOrigName = m_efficiencyCsvUpload->clientFileName().toUTF8();
  if( SpecUtils::iends_with( csvOrigName, ".csv" ) )
    csvOrigName = csvOrigName.substr(0, csvOrigName.size()-4);
  if( det )
    det->setName( csvOrigName );
  
  m_detector = det;
  setAcceptButtonEnabled( true );
  
  emitChangedSignal();
}//void fileUploadedCallback();


void DetectorEdit::relEffDetectorSelectCallback()
{
  m_relEffSelect->userSelectedRelEffDetSelect();
}//void relEffDetectorSelectCallback()


void DetectorEdit::gadrasDetectorSelectCallback()
{
  std::shared_ptr<DetectorPeakResponse> det = m_gadrasDetSelect->selectedDetector();

  m_detector = det;
  setAcceptButtonEnabled( !!det );
  
  emitChangedSignal();
}//void gadrasDetectorSelectCallback()


void DetectorEdit::updateLastUsedTimeOrAddToDb( std::shared_ptr<DetectorPeakResponse> drf,
                                                long long db_user_id,
                                                std::shared_ptr<DataBaseUtils::DbSession> sql )
{
  if( !drf || !sql )
    return;
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    
    uint64_t hash = drf->hashValue();
    auto result = sql->session()->find<DetectorPeakResponse>()
                      .where("InterSpecUser_id = ? AND Hash = ?")
                      .bind( db_user_id )
                      .bind( reinterpret_cast<int64_t&>(hash) )
                      .resultList();
    int nupdated = 0;
    for( auto iter = result.begin(); iter != result.end(); ++iter )
    {
      (*iter).modify()->updateLastUsedTimeToNow();
      ++nupdated;
    }//
    
//    cout << "Have seen Detector '" << drf->name() << "' " << nupdated << " times." << endl;
    
    if( !nupdated )
    {
      //Create a separate DetectorPeakResponse because shared_ptr and
      //  dbo::ptr don't work well together
      DetectorPeakResponse *tempDetector = new DetectorPeakResponse( *drf );
      tempDetector->updateLastUsedTimeToNow();
      
      tempDetector->m_user = db_user_id; //adds the current user to the detectorpeakresponse insertion into DB
      sql->session()->add( tempDetector );
    }
    
    transaction.commit();
  }catch( std::exception &e )
  {
    passMessage( "Error writing DetectorPeakResponse to DB: " + std::string(e.what()), "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void updateLastUsedTimeOrAddToDb(...)



std::shared_ptr<DetectorPeakResponse> DetectorEdit::getUserPrefferedDetector(
                                std::shared_ptr<DataBaseUtils::DbSession> sql,
                                Wt::Dbo::ptr<InterSpecUser> user,
                                std::shared_ptr<const SpecUtils::SpecFile> meas )
{
  std::shared_ptr<DetectorPeakResponse> answer;
  if( !sql || !user || !meas )
    return answer;
  
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    
    const std::string &serial = meas->instrument_id();
    
    Wt::Dbo::ptr<UseDrfPref> pref;
    
    if( !serial.empty() )
    {
      auto results = sql->session()->find<UseDrfPref>()
                     .where( "InterSpecUser_id = ? AND MatchField = ? AND Criteria = ?" )
                     .bind( user.id() )
                     .bind( UseDrfPref::UseDrfType::UseDetectorSerialNumber )
                     .bind( serial )
                     .orderBy("id desc") //shouldnt have an effect because we should only get at most one result
                     .resultList();
      if( results.size() )
        pref = results.front();
    }
    
    if( !pref )
    {
      string model;
      if( meas->detector_type() == DetectorType::Unknown )
        model = meas->instrument_model();
      else
        model = detectorTypeToString( meas->detector_type() );
      if( !model.empty() )
      {
        auto results = sql->session()->find<UseDrfPref>()
                       .where( "InterSpecUser_id = ? AND MatchField = ? AND Criteria = ?" )
                       .bind( user.id() )
                       .bind( UseDrfPref::UseDrfType::UseDetectorModelName )
                       .bind( model )
                       .orderBy("id desc")  //shouldnt have an effect because we should only get at most one result
                       .resultList();
        
        if( results.size() )
          pref = results.front();
      }
    }//if( !results.size() )
    
    if( !pref )
    {
      //cout << "Could not find default DRF preference in database" << endl;
      transaction.commit();
      return answer;
    }
    
    if( pref->m_drfIndex < 0 )
    {
      //cout << "Default DRF preference in database was invalid" << endl;
      pref.remove();
      transaction.commit();
      return answer;
    }
    
    auto drflist = sql->session()->find<DetectorPeakResponse>()
    .where("InterSpecUser_id = ? AND id = ?")
    .bind( user.id() )
    .bind( pref->m_drfIndex )
    .resultList();
    
    for( auto iter = drflist.begin(); !answer && iter != drflist.end(); ++iter )
      answer = std::make_shared<DetectorPeakResponse>( **iter );
    
    transaction.commit();
  }catch( std::exception &e )
  {
    cerr << "Got exception getting user preffered detector from DB: " << e.what() << endl;
  }//try catch see if the user has a preference
  
  //if( answer )
  //  cout << "Got user preffered default det" << endl;
  //else
  //  cout << "Did not get user preffered default det" << endl;
  
  return answer;
};//getUserPrefferedDetector(...)


void DetectorEdit::setUserPrefferedDetector( std::shared_ptr<DetectorPeakResponse> drf,
                                     std::shared_ptr<DataBaseUtils::DbSession> sql,
                                     Wt::Dbo::ptr<InterSpecUser> user,
                                     UseDrfPref::UseDrfType prefType,
                                     std::shared_ptr<const SpecUtils::SpecFile> meas )
{
  //Check if there is a UseDrfPref for the measurement.  If so, we will delete
  //  it.  We will then add one to the database.
  
  if( !meas || !drf || !sql || !user )
    return;
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    
    string criteria;
    switch( prefType )
    {
      case UseDrfPref::UseDrfType::UseDetectorModelName:
        if( meas->detector_type() == DetectorType::Unknown )
          criteria = meas->instrument_model();
        else
          criteria = detectorTypeToString( meas->detector_type() );
        break;
        
      case UseDrfPref::UseDrfType::UseDetectorSerialNumber:
        criteria = meas->instrument_id();
        break;
    }//switch( prefType )
    
    //Remove preferences (there should only be a max of one) from the DB
    auto prevpref = user->drfPrefs().find().where( "MatchField = ? AND Criteria = ?" )
                         .bind( prefType )
                         .bind( criteria )
                         .resultList();
    for( auto iter = prevpref.begin(); iter != prevpref.end(); ++iter )
      iter->remove();
    
    
    //Note: adding the preference to the database, even if criteria is empty
    uint64_t uhash = drf->hashValue();
    int64_t ihash = reinterpret_cast<int64_t&>(uhash);
    auto drflist = sql->session()->find<DetectorPeakResponse>()
                       .where("InterSpecUser_id = ? AND Hash = ?")
                       .bind( user.id() )
                       .bind( ihash )
                       .resultList();
    if( drflist.size() < 1 )
      throw runtime_error( "Could not find detector with hash '" + std::to_string(uhash) + " in DB." );
    
    UseDrfPref *newprefraw = new UseDrfPref();
    Wt::Dbo::ptr<UseDrfPref> newpref( newprefraw );
      
    newprefraw->m_user = user;
    newprefraw->m_field = prefType;
    newprefraw->m_criteria = criteria;
    newprefraw->m_flags = 0x00;
    newprefraw->m_drfIndex = drflist.front().id();
    
    newpref = sql->session()->add( newpref );
    
    auto newid = newpref.id();
    auto drfindex = newpref->m_drfIndex;
    
    transaction.commit();
    
    cout << "Added prefered detector to DB for criteria='" << criteria
         << "', DRF index=" << drfindex << ", as UseDrfPref.id="
         << newid << endl;
  }catch( std::exception &e )
  {
    cerr << "setUserPrefferedDetector: Got exception setting user preffered detector to DB: " << e.what() << endl;
  }//try catch see if the user has a preference
  
}//void setUserPrefferedDetector(...)


void DetectorEdit::acceptAndFinish()
{
  updateLastUsedTimeOrAddToDb( m_detector, m_interspec->m_user.id(), m_sql );
  
  emitChangedSignal();
  emitModifiedSignal();
  
  auto meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  auto sql = m_interspec->sql();
  auto user = m_interspec->m_user;
  auto drf = m_detector;
  
  if( meas && drf && m_defaultForSerialNumber && m_defaultForSerialNumber->isChecked() )
  {
    UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorSerialNumber;
    WServer::instance()->ioService().post( std::bind( [=](){
      DetectorEdit::setUserPrefferedDetector( drf, sql, user, preftype, meas );
    } ) );
  }//if( m_defaultForSerialNumber and is checked )
  
  if( meas && drf && m_defaultForDetectorModel && m_defaultForDetectorModel->isChecked() )
  {
    UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorModelName;
    WServer::instance()->ioService().post( std::bind( [=](){
      DetectorEdit::setUserPrefferedDetector( drf, sql, user, preftype, meas );
    } ) );
  }//if( m_defaultForDetectorModel and is checked )
  
  
  done().emit();
}//void acceptAndFinish()


void DetectorEdit::cancelAndFinish()
{
  m_detector = m_originalDetector;
  if( m_detector && m_originalDetectorCopy )
    (*m_detector) = (*m_originalDetectorCopy);

  emitChangedSignal();
  emitModifiedSignal();
  done().emit();
}//void cancelAndFinish()


void DetectorEdit::finishWithNoDetector()
{
  m_detector.reset();
  emitChangedSignal();
  emitModifiedSignal();
  done().emit();
}//void finishWithNoDetector()


vector<pair<string,string>> DetectorEdit::avaliableGadrasDetectors() const
{
  vector<pair<string,string>> answer;

  //Look through paths specified in user prefernce 'GadrasDRFPath' and find
  //  directoires containing both a 'Detector.dat' and 'Efficiency.csv' file.
  
#if( BUILD_FOR_WEB_DEPLOYMENT )
  const string datadir = InterSpec::staticDataDirectory();
  const string drfpaths = SpecUtils::append_path( datadir, "GenericGadrasDetectors" )
                          + ";" + SpecUtils::append_path( datadir, "OUO_GadrasDetectors" );
#else
  const string drfpaths = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", m_interspec );
#endif
  
  vector<string> paths;
  SpecUtils::split( paths, drfpaths, "\r\n;" );
  
  for( std::string path : paths )
  {
    SpecUtils::trim( path );
    const vector<string> csvs = SpecUtils::recursive_ls( path, "Efficiency.csv" ); //ending string is case-insensitive
    
    for( const std::string &csv : csvs )
    {
      const string parent = SpecUtils::parent_path(csv);
      vector<string> files = SpecUtils::ls_files_in_directory(parent);
      bool found_dat = false;
      for( size_t i = 0; !found_dat && i < files.size(); ++i )
        found_dat = SpecUtils::iends_with( files[i], "Detector.dat" );
      if( found_dat )
      {
        try
        {
          //fs_relative probably shouldnt throw, but JIC.
          //fs_relative( const std::string &from_path, const std::string &to_path )
          string displayname = SpecUtils::fs_relative( path, parent );
          if( displayname.empty() || displayname == "." )
            displayname = SpecUtils::filename(parent);
          
          answer.push_back( make_pair(parent, displayname) );
        }catch(...)
        {
        }
      }//if( found_dat )
    }
  }//for( const std::string &path : paths )

  return answer;
}//vector<string> avaliableGadrasDetectors() const


std::shared_ptr<DetectorPeakResponse> DetectorEdit::initARelEffDetector( const SpecUtils::DetectorType type, InterSpec *interspec )
{
  using SpecUtils::DetectorType;
  
  string smname;
  switch( type )
  {
    case DetectorType::Exploranium:       smname = "GR135";             break;
    case DetectorType::IdentiFinder:      smname = "IdentiFINDER";      break;
    case DetectorType::IdentiFinderNG:    smname = "IdentiFINDER-NGH";  break;
    case DetectorType::IdentiFinderLaBr3: smname = "IdentiFINDER-LaBr"; break;
    case DetectorType::DetectiveUnknown:  smname = "Detective";         break;
    case DetectorType::DetectiveEx:       smname = "Detetive DX";       break;
    case DetectorType::DetectiveEx100:    smname = "Detective EX-100";  break;
    case DetectorType::DetectiveEx200:    smname = "Detective EX-200";  break;
    case DetectorType::DetectiveX:        smname = "Detective X";       break;
    case DetectorType::SAIC8:             smname = "";                  break;
    case DetectorType::Falcon5000:        smname = "Falcon 5000";       break;
    case DetectorType::Unknown:           smname = "";                  break;
    case DetectorType::MicroDetective:    smname = "Micro Detective";   break;
    
    default:
      break;
  }//switch( type )
  
  if( smname.empty() )
    throw std::runtime_error( "Could not find Rel. Eff. detector response functions" );
  
  string concat_filenames;
#if( BUILD_FOR_WEB_DEPLOYMENT )
  concat_filenames = SpecUtils::append_path( InterSpec::staticDataDirectory(), "OUO_lanl_simplemass_detectors.tsv" );
#else
  if( interspec )
    concat_filenames = InterSpecUser::preferenceValue<string>( "RelativeEffDRFPaths", interspec );
#endif
  
  vector<string> paths;
  SpecUtils::split( paths, concat_filenames, "\r\n;" );
  
  
  for( const string &filename : paths )
  {
#ifdef _WIN32
    const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
    std::ifstream input( wfilename.c_str() );
#else
    std::ifstream input( filename.c_str() );
#endif
    
    if( !input.is_open() )
      continue;
  
    string line;
    while( SpecUtils::safe_get_line( input, line, 2048 ) )
    {
      SpecUtils::trim( line );
      if( line.empty() || line[0]=='#' || line.size() < smname.size() )
        continue;
      
      if( SpecUtils::icontains( line, smname ) )
      {
        auto det = RelEffFile::parseDetector( line );
        
        if( det )
        {
          det->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedRelativeEfficiencyDrf );
          return det;
        }
      }//if( thisname == smname )
    }//while( SpecUtils::safe_get_line( input, line ) )
  }//for( const string &filename : paths )
  
  throw runtime_error( "Coudlnt find detector '" + smname + " in relative efficiency files" );
  
  return std::shared_ptr<DetectorPeakResponse>();
}//initARelEffDetector( int type )

std::shared_ptr<DetectorPeakResponse> DetectorEdit::initAGadrasDetector(
                                                            const SpecUtils::DetectorType type, InterSpec *interspec )

{
  using SpecUtils::DetectorType;
  
  string name;
  switch( type )
  {
    case DetectorType::Exploranium:       name = "GR135";              break;
    case DetectorType::IdentiFinder:      name = "identiFINDER-N";     break;
    case DetectorType::IdentiFinderNG:    name = "identiFINDER-NGH";   break;
    case DetectorType::IdentiFinderLaBr3: name = "identiFINDER-LaBr3"; break;
    case DetectorType::DetectiveUnknown:  name = "Detective";          break;
    case DetectorType::DetectiveEx:       name = "Detective-EX";       break;
    case DetectorType::DetectiveEx100:    name = "Detective-EX100";    break;
    case DetectorType::DetectiveEx200:    name = "Ortec IDM Portal";   break;
    case DetectorType::DetectiveX:        name = "Detective X";        break;
    case DetectorType::SAIC8:             name = "";                   break;
    case DetectorType::Falcon5000:        name = "Falcon 5000";        break;
    case DetectorType::Unknown:           name = "";                   break;
    case DetectorType::MicroDetective:    name = "MicroDetective";     break;
      
    //case DetectorType::RadHunterNaI:            name = ""; break;
    //case DetectorType::RadHunterLaBr3:          name = ""; break;
    case DetectorType::Rsi701:                    name = "NaI 2x4x16"; break;
    case DetectorType::Rsi705:                    name = "NaI 2x4x16"; break;
    case DetectorType::AvidRsi:                   name = "NaI 2x4x16"; break;
    case DetectorType::OrtecRadEagleNai:          name = "RadEagle"; break;
    //case DetectorType::OrtecRadEagleCeBr2Inch: name = ""; break;
    //case DetectorType::OrtecRadEagleCeBr3Inch: name = ""; break;
    //case DetectorType::OrtecRadEagleLaBr: name = ""; break;
    //case DetectorType::Sam940LaBr3: name = ""; break;
    case DetectorType::Sam940:                     name = "SAM-945"; break;
    case DetectorType::Sam945:                     name = "SAM-945"; break;
    case DetectorType::Srpm210:                    name = "SRPM-210"; break;
      
    case DetectorType::RIIDEyeNaI:                 name = "RIIDEyeX-GN1"; break;
    case DetectorType::Interceptor:                name = "Interceptor"; break;
      
    case DetectorType::MicroRaider:                name = "Raider"; break;
      
    case DetectorType::RadSeekerNaI:               name = "RadSeeker-NaI"; break;
    case DetectorType::RadSeekerLaBr:              name = "Radseeker-LaBr3"; break;
    
    case DetectorType::VerifinderNaI:              name = "Verifinder-NaI"; break;
    case DetectorType::VerifinderLaBr:             name = "Verifinder-LaBr3"; break;
      
    case DetectorType::IdentiFinderR500NaI:        name = "IdentiFINDER-R500-NaI";   break;
    case DetectorType::IdentiFinderR500LaBr:       name = "IdentiFINDER-R500-LaBr3"; break;
      
      
    case DetectorType::OrtecRadEagleCeBr2Inch:
    case DetectorType::OrtecRadEagleCeBr3Inch:
    case DetectorType::OrtecRadEagleLaBr:
    case DetectorType::RadHunterLaBr3:
    case DetectorType::RadHunterNaI:
    case DetectorType::Sam940LaBr3:
    case DetectorType::IdentiFinderTungsten:
    case DetectorType::IdentiFinderUnknown:
    case DetectorType::RIIDEyeLaBr:
      
      // \TODO: fill in these names
      break;
  }//switch( type )
  
  if( name.empty() )
    throw runtime_error( "There is no GADRAS detector response function for a "
                        + SpecUtils::detectorTypeToString(type) );

  std::shared_ptr<DetectorPeakResponse> det;
  
  try
  {
    det = initAGadrasDetector( name, interspec );
    if( det )
      return det;
  }catch(...)
  {
  }
  
  //We couldnt find an exact match - lets be a little looser
  string dettype, deteff;
  switch( type )
  {
    case DetectorType::Exploranium:
      break;
    case DetectorType::IdentiFinder:
      break;
    case DetectorType::IdentiFinderNG:
      break;
    case DetectorType::OrtecRadEagleNai:
      break;
      
    case DetectorType::IdentiFinderLaBr3:
      break;
      
    case DetectorType::DetectiveUnknown:
    case DetectorType::DetectiveEx:
    case DetectorType::MicroDetective:
      break;
    case DetectorType::DetectiveEx100:
      break;
    case DetectorType::Falcon5000:
      break;
  
    default:
      break;
  }//switch( type )
  
  const string secondname = dettype + " " + deteff;
  
  try
  {
    det = initAGadrasDetector( secondname, interspec );
    if( det )
      return det;
  }catch(...)
  {
  }
  
  throw runtime_error( "Could not find GADRAS detector names " + name + " or " + secondname  );
  return det;
}//initAGadrasDetector


std::shared_ptr<DetectorPeakResponse> DetectorEdit::initAGadrasDetector( const std::string &currentDetName, InterSpec *interspec )
{
  //Grab "GadrasDRFPath" and split it by semicolon and newlines, then go through
  //  and look for sub-folders with the given name that contain Detector.data
  //  and Efficiency.csv.
#if( BUILD_FOR_WEB_DEPLOYMENT )
  const string datadir = InterSpec::staticDataDirectory();
  const string drfpaths = SpecUtils::append_path( datadir, "GenericGadrasDetectors" )
                          + ";" + SpecUtils::append_path( datadir, "OUO_GadrasDetectors" );
#else
  const string drfpaths = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", interspec );
#endif
  
  
  vector<string> paths;
  SpecUtils::split( paths, drfpaths, "\r\n;" );
  for( string basepath : paths )
  {
    SpecUtils::trim( basepath );
    const string path = SpecUtils::append_path(basepath, currentDetName);
    if( !SpecUtils::is_directory(path) )
      continue;
    
    string thiscsv, thisdat;
    const vector<string> files = SpecUtils::ls_files_in_directory( path );
    for( size_t i = 0; (thiscsv.empty() || thisdat.empty()) && i < files.size(); ++i )
    {
      //We will be loose with upper/lower case
      if( SpecUtils::iends_with(files[i], "Efficiency.csv") )
        thiscsv = files[i];
      else if( SpecUtils::iends_with(files[i], "Detector.dat") )
        thisdat = files[i];
    }//for( size_t i = 0; i < files.size(); ++i )
    
    if( thiscsv.size() && thisdat.size() )
    {
      auto det = initAGadrasDetectorFromDirectory( path, interspec );
      if( det )
      {
        if( SpecUtils::icontains( basepath, "GenericGadrasDetectors") )
          det->setDrfSource( DetectorPeakResponse::DrfSource::DefaultGadrasDrf );
        else
          det->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedGadrasDrf );
      }
      
      return det;
    }
  }//for( const string &basepath : paths )
  
  throw runtime_error( "Could not find GADRAS detector names " + currentDetName );
}//std::shared_ptr<DetectorPeakResponse> initAGadrasDetector( const std::string &name );


std::shared_ptr<DetectorPeakResponse> DetectorEdit::initAGadrasDetectorFromDirectory( const std::string &path, InterSpec *interspec )
{
  auto det = std::make_shared<DetectorPeakResponse>();

  if( !SpecUtils::is_directory(path) )
    throw runtime_error( "'" + path + "' is not a directory" );
  
  //Lets be lienient on upper vs lower case.
  string thiscsv, thisdat;
  const vector<string> files = SpecUtils::ls_files_in_directory( path );
  for( size_t i = 0; (thiscsv.empty() || thisdat.empty()) && i < files.size(); ++i )
  {
    //We will be loose with upper/lower case
    if( SpecUtils::iends_with(files[i], "Efficiency.csv") )
      thiscsv = files[i];
    else if( SpecUtils::iends_with(files[i], "Detector.dat") )
      thisdat = files[i];
  }//for( size_t i = 0; i < files.size(); ++i )
  
  if( thiscsv.empty() )
    throw runtime_error( "Could not find a Efficiency.csv file in '" + path + "'" );
  
  if( thisdat.empty() )
    throw runtime_error( "Could not find a Detector.dat file in '" + path + "'" );
  
#ifdef _WIN32
  const std::wstring wthiscsv = SpecUtils::convert_from_utf8_to_utf16(thiscsv);
  ifstream effstrm( wthiscsv.c_str(), ios_base::binary|ios_base::in );
#else
  ifstream effstrm( thiscsv.c_str(), ios_base::binary|ios_base::in );
#endif
  
  if( !effstrm.is_open() )
    throw runtime_error( "Could not open '" + thiscsv + "'" );

#ifdef _WIN32
  const std::wstring wthisdat = SpecUtils::convert_from_utf8_to_utf16(thisdat);
  ifstream datstrm( wthisdat.c_str(), ios_base::binary|ios_base::in );
#else
  ifstream datstrm( thisdat.c_str(), ios_base::binary|ios_base::in );
#endif
  
  if( !datstrm.is_open() )
    throw runtime_error( "Could not open '" + thisdat + "'" );

  det->fromGadrasDefinition( effstrm, datstrm );
  det->setName( SpecUtils::filename( SpecUtils::parent_path(thiscsv) ) );

  return det;
}//std::shared_ptr<DetectorPeakResponse> initAGadrasDetectorFromDirectory()


void DetectorEdit::emitChangedSignal()
{
  //Make sure this is a necessary signal to emit
  if( m_previousEmmittedDetector == m_detector )
    return;

  //Now update m_previousDetectorDef and m_previousEmmittedDetector so we can
  //  appropriately filter out unecessary signal emits in the future
  if( m_detector )
    m_previousDetectorDef = *m_detector;
  else
    m_previousDetectorDef.reset();
  m_previousEmmittedDetector = m_detector;

  m_interspec->detectorChanged().emit( m_detector );
  setAcceptButtonEnabled( !!m_detector );
  
  updateChart();
}//void emitChangedSignal()


void DetectorEdit::emitModifiedSignal()
{
  //Make sure this is a necessary signal to emit
  if( !m_detector && !m_previousDetectorDef.isValid() )
    return;

  if( m_detector && ((*m_detector)==m_previousDetectorDef) )
    return;

  //Now update m_previousDetectorDef and m_previousEmmittedDetector so we can
  //  appropriately filter out unecessary signal emits in the future
  if( m_detector )
    m_previousDetectorDef = *m_detector;
  else
    m_previousDetectorDef.reset();
  m_previousEmmittedDetector = m_detector;

  m_interspec->detectorChanged().emit( m_detector );
}//void emitModifiedSignal()


Wt::Signal<> &DetectorEdit::done()
{
  return m_doneSignal;
}


DetectorEditWindow::DetectorEditWindow(
                                  std::shared_ptr<DetectorPeakResponse> det,
                                  InterSpec *specViewer,
                                  SpectraFileModel *fileModel )
  : AuxWindow("Detector Response Function Select",
              Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal) | AuxWindowProperties::TabletNotFullScreen),
    m_edit( NULL )
{
  disableCollapse();
  m_edit = new DetectorEdit( det, specViewer, fileModel, this );
  contents()->setPadding( 0 );
  stretcher()->addWidget( m_edit, 0, 0 );
  stretcher()->setContentsMargins( 6, 2, 6, 0 );
  m_edit->done().connect( boost::bind( &DetectorEditWindow::acceptAndDelete, this ) );
  finished().connect( boost::bind( &DetectorEditWindow::acceptAndDelete, this ) );
  rejectWhenEscapePressed();
  setMargin( 0 );
  if( specViewer->isTablet() )
    titleBar()->hide();
  resize( WLength(650,WLength::Pixel), WLength(610,WLength::Pixel));
  centerWindow();
  show();
  resizeToFitOnScreen();
  setResizable(true);
//  resizeWindow( 300, 400 );
}//DetectorEditWindow constructor


DetectorEditWindow::~DetectorEditWindow()
{
}


void DetectorEditWindow::acceptAndDelete( DetectorEditWindow *window )
{
  window->m_edit->emitChangedSignal();  //will only emit if necessary
  window->m_edit->emitModifiedSignal(); //will only emit if necessary
  AuxWindow::deleteAuxWindow( window );
}//void acceptAndDelete( DetectorEditWindow *window )
