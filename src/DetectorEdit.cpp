/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <boost/tokenizer.hpp>
#include <boost/assign/list_of.hpp>

#include <Wt/WText>
#include <Wt/WBreak>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WSignal>
#include <Wt/WTextArea>
#include <Wt/WTreeView>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WTabWidget>
#include <Wt/WFileUpload>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WRadioButton>
#include <Wt/WButtonGroup>
#include <Wt/WItemDelegate>
#include <Wt/WBorderLayout>
#include <Wt/WStackedWidget>
#include <Wt/Dbo/QueryModel>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>

#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectorEdit.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "SpecUtils/UtilityFunctions.h"
#include "SpecUtils/SpectrumDataStructs.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;

const char * const DetectorDisplay::sm_noDetectorTxt
  = "<font style=\"font-weight:100;color:#CFCFCF;\">&lt;click to select&gt;</font>";


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
  
  void addFile();
  
  void removeFile( RelEffFile *fileWidget );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  void saveFilePathToUserPreferences();
#endif
  
  void trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  void detectorSelected( RelEffFile *file, std::shared_ptr<DetectorPeakResponse> det );
  
  void docreate();
  
  //Attempts to do a defered rednering, but it looks like thisisnt actually
  //  the case
  void load();
};//class RelEffDetSelect


/** Represents a single CSV/TSV DRF file, which may have multiple DRFs.
 */
class RelEffFile : public Wt::WContainerWidget
{
  //If we are deploying on the web, we do not want to allow the user access to
  //  the filesystem(!), so we will only allow using detectors from
  //  "data/lanl_simplemass_detectors.tsv"
#if( BUILD_FOR_WEB_DEPLOYMENT )
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
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  void filePathChanged();
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )
  
  void detectorSelectCallback();
  
  
  bool trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  //parseDetector() returns null on error
  static std::shared_ptr<DetectorPeakResponse> parseDetector( string &line );
  void initDetectors();
  void detectorSelected( const int index );
};//class RelEffFile


class GadrasDetSelect : public Wt::WContainerWidget
{
protected:
  InterSpec *m_interspec;
  DetectorEdit *m_detectorEdit;
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  Wt::WText *m_message;
  Wt::WLineEdit *m_baseDir;
#endif
  Wt::WComboBox *m_detectorSelect;
  
public:
  GadrasDetSelect( InterSpec *interspec, DetectorEdit *detedit, WContainerWidget *parent = 0 );
  
  virtual ~GadrasDetSelect(){}
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  void baseDirChanged();
#endif
  
  void updateWidgetsFromInputPaths();
  
  bool trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  std::shared_ptr<DetectorPeakResponse> selectedDetector();
};//class RelEffDetSelect


////////////////////////////////////////////////////////////////////////////////
//////////////////   Begin RelEffFile implementation   /////////////////////////
////////////////////////////////////////////////////////////////////////////////
RelEffFile::RelEffFile( std::string file,
                       RelEffDetSelect *parentSelect,
                       DetectorEdit *detectorEdit,
                       WContainerWidget *parent )
: WContainerWidget( parent ),
#if( BUILD_FOR_WEB_DEPLOYMENT )
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
  
#if( BUILD_FOR_WEB_DEPLOYMENT )
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
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  m_fileEdit->changed().connect( m_detectorSelect, &WPushButton::disable );
  //m_fileEdit->keyPressed().connect( m_detectorSelect, &WPushButton::disable );
#endif
  
  initDetectors();
}//RelEffFile


const std::string RelEffFile::filepath()
{
#if( BUILD_FOR_WEB_DEPLOYMENT )
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

#if( !BUILD_FOR_WEB_DEPLOYMENT )
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
    if( m_responses[i]->name() == det->name() )
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
  
  //Allow the file to be comma or tab delimited, but if an individual field
  //  contains a comma or tab, then the field must be quoted by a double
  //  quote.  Note that if you just copy cells from Microsoft Excel, that
  //  contain a comma, and then past into a text editor, fields with a comma
  //  are not quoted.
  typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
  boost::escaped_list_separator<char> separator("\\",",\t", "\"");
  Tokeniser t( line, separator );
  for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
    fields.push_back(*it);
  
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
    
    det.reset( new DetectorPeakResponse() );
    det->fromExpOfLogPowerSeriesAbsEff( coefs, dist, diam, eunits );
    det->setEfficiencySource( DetectorPeakResponse::kRelativeEfficiencyDefintion );
    
    det->setName( name );
    det->setDescription( fields[2] + " - from Relative Eff. File" );
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
#if( BUILD_FOR_WEB_DEPLOYMENT )
  pathstr = m_file;
#else
  pathstr = m_fileEdit->text().toUTF8();
#endif
  
  std::ifstream input( pathstr.c_str(), ios::in | ios::binary );
  
  const bool file_opened = input.is_open();
  
  if( file_opened )
  {
    string line;
    while( UtilityFunctions::safe_get_line( input, line ) )
    {
      UtilityFunctions::trim( line );
      
      if( UtilityFunctions::istarts_with( line, "#credit:") )
        credits.push_back( UtilityFunctions::trim_copy(line.substr(8)) );
      
      if( line.empty() || line[0]=='#' )
        continue;
      
      std::shared_ptr<DetectorPeakResponse> det = parseDetector( line );
      
      if( det )
        m_responses.push_back( det );
    }//while( UtilityFunctions::safe_get_line( input, line ) )
    
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
    credits.push_back( "<span style=\"color:red;\">Could not open DRF file.</span>" );
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
      m_files( 0 )
{
}


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
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  saveFilePathToUserPreferences();
#endif
}//void removeFile( RelEffFile )


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


#if( !BUILD_FOR_WEB_DEPLOYMENT )
void RelEffDetSelect::saveFilePathToUserPreferences()
{
  auto children = m_files->children();
  
  string concat_path;
  for( auto w : children )
  {
    auto child = dynamic_cast<RelEffFile *>( w );
    const string filepath = child ? child->filepath() : string("");
    if( filepath.size() )
      concat_path += (concat_path.size() ? ";" : "") + filepath;
  }//string concat_path;
  
  try
  {
    InterSpecUser::setPreferenceValue( m_interspec->m_user, "RelativeEffDRFPaths", concat_path, m_interspec );
  }catch( std::exception &e )
  {
    passMessage( "Error saving Rel. Eff. files to user preferences: " + string(e.what()), "DetectorEdit", WarningWidget::WarningMsgHigh );
  }
}//void saveFilePathToUserPreferences()
#endif //#if( !BUILD_FOR_WEB_DEPLOYMENT )


void RelEffDetSelect::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
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
  
  WPushButton *addIcon = new WPushButton();
  addIcon->setStyleClass( "AddRelEffFile" );
  addIcon->clicked().connect( this, &RelEffDetSelect::addFile );
  WContainerWidget *holder = new WContainerWidget( this );
  holder->addWidget( addIcon );
  holder->setWidth( WLength(100,WLength::Percentage) );
  holder->setToolTip( "Click the plus button to add an additonal relative efficiency detector response function file." );
  
  string pathstr;
  vector<string> paths;
#if( BUILD_FOR_WEB_DEPLOYMENT )
  pathstr = UtilityFunctions::append_path( InterSpec::dataDirectory(), "lanl_simplemass_detectors.tsv" );
#else
  try
  {
    if( m_interspec )
      pathstr = InterSpecUser::preferenceValue<string>( "RelativeEffDRFPaths", m_interspec );
  }catch( std::exception &e )
  {
    passMessage( "Error retrieving 'RelativeEffDRFPaths' preference.", "", WarningWidget::WarningMsgHigh );
  }
#endif

  UtilityFunctions::split( paths, pathstr, "\r\n;" );
  
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
GadrasDetSelect::GadrasDetSelect( InterSpec *interspec, DetectorEdit *detedit, WContainerWidget *parent )
  : WContainerWidget(parent),
    m_interspec( interspec ),
    m_detectorEdit( detedit ),
#if( !BUILD_FOR_WEB_DEPLOYMENT )
    m_message( nullptr ),
    m_baseDir( nullptr ),
#endif
    m_detectorSelect( nullptr )
{
  //GADRAS-DRF URL: https://rsicc.ornl.gov/codes/psr/psr6/psr-610.html
  
  const char *gadrasToolTip = "These are detector response functions imported"
  " from GADRAS-DRF.  Currently only the detector dimensions, "
  " absolute efficiencies, and FWHM resolutions are"
  " used.";
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  //Need to point the GUI to the appropriate directory, and implement to an `ls` to find detctors with both Detcotr.dat and Efficy.csv.
  const string drfpaths = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", m_interspec );
  
  WContainerWidget *pathsDiv = new WContainerWidget( this );
  WLabel *label = new WLabel( "Paths to recursively look for DRFs.", pathsDiv );
  label->setInline( false );
  label = new WLabel( "Seperate multiple paths with semicolons:", pathsDiv );
  label->setInline( false );
  m_baseDir = new WLineEdit( pathsDiv );
  m_baseDir->setInline( false );
  m_baseDir->setWidth( WLength(95,WLength::Percentage) );
  m_baseDir->setText( drfpaths );
  m_baseDir->changed().connect( this, &GadrasDetSelect::baseDirChanged );
  m_baseDir->enterPressed().connect( this, &GadrasDetSelect::baseDirChanged );
  m_message = new WText( "", pathsDiv );
  m_message->setInline( false );
  
  pathsDiv->setMargin( 10, Wt::Top | Wt::Bottom );
#endif //#if( !BUILD_FOR_WEB_DEPLOYMENT )
  
  WContainerWidget *selectDiv = new WContainerWidget( this );
  label = new WLabel( "GADRAS detectors: ", selectDiv );
  
  m_detectorSelect = new WComboBox( selectDiv );
  label->setBuddy( m_detectorSelect );
  
  WText *useInfo = new WText( gadrasToolTip, selectDiv );
  useInfo->setStyleClass("DetectorLabel");
  useInfo->setMargin( 5, Wt::Top );
  useInfo->setInline( false );
  
  m_detectorSelect->activated().connect( m_detectorEdit, &DetectorEdit::gadrasDetectorSelectCallback );

  updateWidgetsFromInputPaths();
}//GadrasDetSelect(...)


#if( !BUILD_FOR_WEB_DEPLOYMENT )
void GadrasDetSelect::baseDirChanged()
{
  const string enteredPath = m_baseDir->text().toUTF8();
  
  try
  {
    InterSpecUser::setPreferenceValue( m_interspec->m_user, "GadrasDRFPath", enteredPath, m_interspec );
  }catch( std::exception &e )
  {
    passMessage( "Error saving GADRAS path to user preferences: " + string(e.what()), "", WarningWidget::WarningMsgHigh );
  }
  
  updateWidgetsFromInputPaths();
}//void baseDirChanged()
#endif


void GadrasDetSelect::updateWidgetsFromInputPaths()
{
  m_detectorSelect->clear();
  m_detectorSelect->addItem( "<None>" );
  
  WAbstractItemModel *m = m_detectorSelect->model();
  
  for( const auto &file_display : m_detectorEdit->avaliableGadrasDetectors() )
  {
    const int index = m_detectorSelect->count();
    m_detectorSelect->addItem( file_display.second );
    m->setData( m->index(index, 0), boost::any(file_display.first), Wt::UserRole );
  }
  
  m_detectorSelect->setCurrentIndex( -1 );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  const string drfpaths = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", m_interspec );
  string msg;
  vector<string> paths;
  UtilityFunctions::split( paths, drfpaths, "\r\n;" );
  for( string p : paths )
  {
    UtilityFunctions::trim( p );
    if( !UtilityFunctions::is_directory(p) )
      msg += "<div style=\"color:red;\">'" + p + "' is not a valid directory.</div>";
  }
  
  m_message->setText( msg );
#endif
}//void updateWidgetsFromInputPaths()

  

bool GadrasDetSelect::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  if( !det || det->name().empty() )
  {
    m_detectorSelect->setCurrentIndex(0);
    return false;
  }
  
  const int index = m_detectorSelect->findText( det->name() );
  if( index <= 0 )
  {
    m_detectorSelect->setCurrentIndex(0);
    return false;
  }
  
  m_detectorSelect->setCurrentIndex( index );
  
  return true;
}//void trySelectDetector( det )


std::shared_ptr<DetectorPeakResponse> GadrasDetSelect::selectedDetector()
{
  std::shared_ptr<DetectorPeakResponse> answer;
  const int currentIndex = m_detectorSelect->currentIndex();
  
  if( currentIndex < 1 )
    return answer;
  
  WAbstractItemModel *m = m_detectorSelect->model();
  
  try
  {
    string p = boost::any_cast<std::string>( m->data( currentIndex, 0, Wt::UserRole ) );
    answer = DetectorEdit::initAGadrasDetectorFromDirectory( p, m_interspec );
  }catch( std::exception &e )
  {
    passMessage( "Failed to parse a GADRAS detector: " + string(e.what()), "", WarningWidget::WarningMsgHigh );
  }//try / catch
  
  return answer;
}//std::shared_ptr<DetectorPeakResponse> selectedDetector()




DetectorDisplay::DetectorDisplay( InterSpec *specViewer,
                                  SpectraFileModel *fileModel,
                                  WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_text( NULL ),
    m_interspec( specViewer ),
    m_fileModel( fileModel )
{
  const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  HelpSystem::attachToolTipOn( this, "The currently selected detector.  Click to select a"
                              " different detector, or modify this one.", showToolTipInstantly );
  
  addStyleClass( "DetectorDisplay" );

  WText* label = new WText( "Detector: ", this );
  label->addStyleClass("CameraIcon");
  m_text = new WText( sm_noDetectorTxt, XHTMLUnsafeText, this );

  std::shared_ptr<DetectorPeakResponse> detector;
  if( m_interspec->measurment(kForeground) )
    detector = m_interspec->measurment(kForeground)->detector();
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
    m_text->setText( sm_noDetectorTxt );
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
  : WContainerWidget( auxWindow->contents() ),
    m_footer (auxWindow->footer()),
    m_interspec( specViewer ),
    m_fileModel( fileModel ),
    m_chart( 0 ),
    m_efficiencyModel( 0 ),
    m_detector( currentDet ),
    m_tabs( nullptr ),
    m_detectorDiameter( nullptr ),
    m_detectrDiameterDiv( nullptr ),
    m_efficiencyCsvUpload( nullptr ),
    m_detectorDotDatUpload( nullptr ),
    m_acceptButton( nullptr ),
    m_cancelButton( nullptr ),
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
    m_group( nullptr ),
    m_stack( nullptr ),
    m_deleteButton( nullptr ),
    m_DBtable( nullptr ),
    m_model( nullptr )
{
  const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  setOverflow(Wt::WContainerWidget::OverflowAuto);
  WBorderLayout * mainLayout = new WBorderLayout();
  setLayout( mainLayout );
  mainLayout->setSpacing( 0 );
  mainLayout->setContentsMargins(0, 0, 0, 0);
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
  
  /*
   * Provide ample space for the title, the X and Y axis and the legend.
   */
  m_chart->setPlotAreaPadding(40, Wt::Top);
  m_chart->setPlotAreaPadding(70, Wt::Bottom);
  m_chart->setPlotAreaPadding(55, Wt::Right | Wt::Left);

  m_chart->setMinimumSize(WLength(600), WLength(300));
  m_chart->setMargin(Wt::WLength::Auto, Wt::Left | Wt::Right);

  m_chart->axis(Wt::Chart::XAxis).setVisible(true);
  m_chart->axis(Wt::Chart::XAxis).setTitle("Energy (keV)");

  m_chart->axis(Wt::Chart::Y1Axis).setVisible(true);
  m_chart->axis(Wt::Chart::Y1Axis).setTitle("Efficiency");

#if( WT_VERSION >= 0x3030400 )
  m_chart->axis(Wt::Chart::Y1Axis).setTitleOrientation( Wt::Vertical );
  m_chart->axis(Wt::Chart::Y2Axis).setTitleOrientation( Wt::Vertical );
#endif

  m_chart->setLegendEnabled(true);

  m_chart->setLegendLocation(Wt::Chart::LegendOutside, Wt::Bottom, Wt::AlignRight);
  
  m_chart->setTitle("Detector Energy");
  mainLayout->addWidget(m_chart,Wt::WBorderLayout::Center);
  
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  distValidator->setInvalidBlankText( "0.0 cm" );

  
  Wt::WContainerWidget *southContainer = new Wt::WContainerWidget();
  
  mainLayout->addWidget( southContainer, Wt::WBorderLayout::South );
  
  WBorderLayout * southLayout =  new WBorderLayout();
  southContainer->setLayout(southLayout);
  
  Wt::WContainerWidget *radioContainer = new Wt::WContainerWidget();
  radioContainer->addStyleClass( "ManualDetectorRadio" );
  southLayout->addWidget( radioContainer, Wt::WBorderLayout::West );
  
  m_group = new Wt::WButtonGroup( radioContainer );
  Wt::WRadioButton * button = new Wt::WRadioButton( "GADRAS", radioContainer );
  button->setInline( false );
  button->setMargin( 10, Wt::Top );
  m_group->addButton(button);  
  button = new Wt::WRadioButton( "Rel. Eff.", radioContainer);
  button->setInline( false );
  button->setMargin( 10, Wt::Top );
  m_group->addButton(button);


  button = new Wt::WRadioButton( "Import", radioContainer);
  button->setInline( false );
  button->setMargin( 10, Wt::Top );
  m_group->addButton(button);
  
  if( !m_interspec->isSupportFile() )
    button->setHidden(true);

  button = new Wt::WRadioButton( "Formula", radioContainer);
  button->setInline( false );
  button->setMargin( 10, Wt::Top );
  m_group->addButton(button);
  button = new Wt::WRadioButton( "Recent", radioContainer);
  button->setInline( false );
  button->setMargin( 10, Wt::Top );
  m_group->addButton(button);
  
  m_stack = new Wt::WStackedWidget();
  southLayout->addWidget( m_stack, Wt::WBorderLayout::Center);
  m_group->checkedChanged().connect(boost::bind( &DetectorEdit::selectButton, this, m_stack, m_group,true ));
  
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

  
  WContainerWidget *uploadDetTab = new WContainerWidget(  );

  WLabel *label = new WLabel( "Intrinsic Efficiency File:", uploadDetTab );
  label->setInline( false );
  
  const char *descrstr =
  "<div>The first column of the efficiency CSV file should be the energy in keV"
  ", and the second column should be the percent efficient (decimal number 0"
  " through 100) for gammas at that energy that strike the detector face, to be"
  " recorded in the photopeak. Further columns will be ignored.</div>"
  "<div>You will also need to upload a GADRAS Detector.dat file, or enter the"
  " detectors diameter, after uploading the efficiency file.</div>";
  
  WText *descrip = new WText( descrstr, uploadDetTab );
  descrip->setInline( false );
  descrip->setStyleClass("DetectorLabel");
  descrip->setInline( false );
  
  
  m_efficiencyCsvUpload = new WFileUpload( uploadDetTab );
  m_efficiencyCsvUpload->setInline( false );
  m_efficiencyCsvUpload->uploaded().connect( boost::bind( &DetectorEdit::fileUploadedCallback, this ) );
  m_efficiencyCsvUpload->fileTooLarge().connect( boost::bind(&SpecMeasManager::fileTooLarge, _1) );
  m_efficiencyCsvUpload->changed().connect( m_efficiencyCsvUpload, &WFileUpload::upload );
  m_efficiencyCsvUpload->setInline( false );
 
  m_detectrDiameterDiv = new WContainerWidget( uploadDetTab );
  m_detectrDiameterDiv->addStyleClass( "DetectorDiamDiv" );
  label = new WLabel( "Enter detector diameter or upload Detector.dat file:", m_detectrDiameterDiv );
  label->setInline( false );

  new WLabel( "Detector Diameter:", m_detectrDiameterDiv );
  m_detectorDiameter = new WLineEdit( "0 cm", m_detectrDiameterDiv );
  
  m_detectorDiameter->setValidator( distValidator );
  m_detectorDiameter->setTextSize( 10 );
  m_detectorDiameter->changed().connect( this, &DetectorEdit::fileUploadedCallback );
  m_detectorDiameter->enterPressed().connect( this, &DetectorEdit::fileUploadedCallback );
  m_detectorDiameter->blurred().connect( this, &DetectorEdit::fileUploadedCallback );

  WContainerWidget *datDiv = new WContainerWidget( m_detectrDiameterDiv );
  datDiv->addStyleClass( "DetectorDotDatDiv" );

  label = new WLabel( "GADRAS Detector.Dat File:", datDiv );
  label->setInline( false );
  m_detectorDotDatUpload = new WFileUpload( datDiv );
  m_detectorDotDatUpload->setInline( false );
  m_detectorDotDatUpload->uploaded().connect( boost::bind( &DetectorEdit::fileUploadedCallback, this ) );
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
    
  WContainerWidget *selfDefine = new WContainerWidget();
  WGridLayout *selfDefineLayout = new WGridLayout();
  selfDefineLayout->setContentsMargins( 0, 0, 0, 0 );
  selfDefine->setLayout(selfDefineLayout);

  WText *selfDefineLabel = new WText( txt );
  selfDefineLayout->addWidget( selfDefineLabel, 0, 0, 1, 3 );
  selfDefineLabel->setStyleClass("DetectorLabel");
  
  selfDefineLayout->addWidget( new WLabel("Detector name "), 1, 0 );
  m_detectorManualFunctionName = new Wt::WLineEdit();
  selfDefineLayout->addWidget( m_detectorManualFunctionName, 1, 1, 1, 2 );
  m_detectorManualFunctionName->setStyleClass("DetectorEditFunctionalFormText");
  m_detectorManualFunctionName->setText( manualdetname );
  m_detectorManualFunctionName->setEmptyText("Unique Detector Name");
  m_detectorManualFunctionName->setInline(false);
  m_detectorManualFunctionName->keyWentUp().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));

  selfDefineLayout->addWidget( new WLabel("Efficiency f(x) = "), 2, 0 );
  m_detectorManualFunctionText = new Wt::WTextArea();
  selfDefineLayout->addWidget( m_detectorManualFunctionText, 2, 1, 1, 1 );
  m_detectorManualFunctionText->setColumns(60);
  m_detectorManualFunctionText->setRows(3);
  m_detectorManualFunctionText->setStyleClass("DetectorEditFunctionalFormText");
  m_detectorManualFunctionText->setText(fcn);
  m_detectorManualFunctionText->setEmptyText(fcn);
  m_detectorManualFunctionText->setInline(false);
  m_detectorManualFunctionText->keyWentUp().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
    
  Wt::WContainerWidget *energyContainer = new Wt::WContainerWidget();
  m_eqnEnergyGroup = new WButtonGroup( energyContainer );
  button = new WRadioButton( "keV", energyContainer );
  m_eqnEnergyGroup->addButton( button, 0 );
  button = new WRadioButton( "MeV", energyContainer );
  m_eqnEnergyGroup->addButton( button, 1 );
  m_eqnEnergyGroup->checkedChanged().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  m_eqnEnergyGroup->setSelectedButtonIndex( 0 );
  selfDefineLayout->addWidget( energyContainer, 2, 2 );
    
  HelpSystem::attachToolTipOn( energyContainer,"Energy unit for efficiency formula" , showToolTipInstantly, HelpSystem::Top );

  selfDefineLayout->addWidget( new WLabel("Description"), 3, 0 );
  m_detectorManualDescription = new WLineEdit();
  selfDefineLayout->addWidget( m_detectorManualDescription, 3, 1, 1, 2 );
  m_detectorManualDescription->setEmptyText( "Description of this detector response" );
  m_detectorManualDescription->changed().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  
  selfDefineLayout->addWidget( new WLabel("Detector diameter "), 4, 0 );
  m_detectorManualDiameterText = new Wt::WLineEdit( selfDefine );
  m_detectorManualDiameterText->setValidator( distValidator );
  
  selfDefineLayout->addWidget( m_detectorManualDiameterText, 4, 1, 1,2 );
  m_detectorManualDiameterText->setText(diamtxt);
  m_detectorManualDiameterText->setEmptyText(diamtxt);
  m_detectorManualDiameterText->setValidator( distValidator );
  m_detectorManualDiameterText->blurred().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  m_detectorManualDiameterText->enterPressed().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  
  Wt::WContainerWidget *formulaTypeContainer = new Wt::WContainerWidget();
  m_absOrIntrinsicGroup = new WButtonGroup( formulaTypeContainer );
  button = new WRadioButton( "Intrinsic", formulaTypeContainer );
  m_absOrIntrinsicGroup->addButton( button, 0 );
  button = new WRadioButton( "Absolute", formulaTypeContainer );
  button->setMargin( 5, Wt::Left );
  m_absOrIntrinsicGroup->addButton( button, 1 );
  m_absOrIntrinsicGroup->checkedChanged().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  m_absOrIntrinsicGroup->setSelectedButtonIndex( 0 );

  m_detectorManualDistLabel = new WLabel( "Dist:", formulaTypeContainer );
  m_detectorManualDistLabel->setMargin( 10, Wt::Left );
  
  m_detectorManualDistText = new Wt::WLineEdit( formulaTypeContainer );
  m_detectorManualDistText->setValidator( distValidator );
  m_detectorManualDistText->setEmptyText("1 m");
  m_detectorManualDistText->setHiddenKeepsGeometry( true );
  m_detectorManualDistLabel->hide();
  m_detectorManualDistText->hide();
  m_detectorManualDistText->keyWentUp().connect(boost::bind(&DetectorEdit::verifyManualDefinition, this));
  
  selfDefineLayout->addWidget( formulaTypeContainer, 5, 0, 1, 2 );
  

  m_manualSetButton = new Wt::WPushButton("Set");
  m_manualSetButton->setInline(false);
  selfDefineLayout->addWidget( m_manualSetButton, 5, 2 );
  selfDefineLayout->setColumnStretch( 1, 1 );
  m_manualSetButton->clicked().connect(boost::bind(&DetectorEdit::setDefineDetector, this));
    
  //-------------------------------------
  //--- 5)  Recent
  //-------------------------------------
  WContainerWidget *recentDiv = new WContainerWidget( );
  WGridLayout* recentDivLayout = new WGridLayout();
  recentDivLayout->setContentsMargins( 0, 0, 0, 0 );
  recentDiv->setLayout(recentDivLayout);

  Dbo::ptr<InterSpecUser> user = m_interspec->m_user;
  
  m_DBtable = new Wt::WTreeView();
  m_DBtable->setRootIsDecorated	(	false ); //makes the tree look like a table! :)
  
  m_DBtable->addStyleClass( "DbSpecFileSelectTable" );
  
  //We have to create a independant Dbo::Session for this class since the
  //  m_viewer->m_user.session() could be used in other threads, messing
  //  things up (Dbo::Session is not thread safe).
  m_sql.reset( new DataBaseUtils::DbSession( *m_interspec->sql() ) );

  m_model = new Dbo::QueryModel< Dbo::ptr<DetectorPeakResponse> >( m_DBtable );

  m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                           .where( "InterSpecUser_id = ? ").bind( user.id() )
                           .orderBy("-1*id") );
  
  m_model->addColumn( "m_name" );
  m_model->addColumn( "m_description" );
  m_model->addColumn( "m_detectorDiameter" );
  
  
  m_model->setHeaderData(  0, Horizontal, WString("Name"), DisplayRole );
  m_DBtable->setColumnWidth( 0, 200 );
  m_model->setHeaderData(  1, Horizontal, WString("Description"), DisplayRole );
  m_DBtable->setColumnWidth( 1, 150 );
  m_model->setHeaderData(  2, Horizontal, WString("Diameter"), DisplayRole );
  m_DBtable->setColumnWidth( 2, 150 );
 
  
  WItemDelegate *delegate = new WItemDelegate( m_DBtable );
  delegate->setTextFormat( "%.2f" );
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
  m_deleteButton->setStyleClass("CrossIcon");
  m_deleteButton->setDefault(true);
  m_deleteButton->disable();
  m_deleteButton->clicked().connect( this, &DetectorEdit::deleteDBTableSelected );

  recentDivLayout->addWidget(m_DBtable,0,0);
  recentDivLayout->setRowStretch(0, 1);
  recentDivLayout->addWidget(m_deleteButton,1,0,AlignLeft);
  recentDiv->setOverflow(Wt::WContainerWidget::OverflowHidden);
  recentDiv->setMaximumSize( WLength::Auto, 190 );
  
  //--------------------------------------------------------------------------------
  
  m_stack->addWidget( m_gadrasDetSelect );
  m_stack->addWidget( m_relEffSelect );
  m_stack->addWidget( uploadDetTab );
  m_stack->addWidget( selfDefine );
  m_stack->addWidget( recentDiv);
  m_stack->setMinimumSize(WLength::Auto, WLength(185.0));
  m_stack->setMaximumSize(WLength::Auto, WLength(185.0));
  
  if( !m_footer )
    m_footer = new WContainerWidget( southContainer );

 // m_footer->setStyleClass("modal-footer");
  
  AuxWindow::addHelpInFooter(m_footer, "detector-edit-dialog",auxWindow);
  m_cancelButton = auxWindow->addCloseButtonToFooter("Cancel");

  HelpSystem::attachToolTipOn( m_cancelButton,"Remove all changes or selections made by this dialog, and close the dialog" , showToolTipInstantly );
  
  m_cancelButton->clicked().connect( this, &DetectorEdit::cancelAndFinish );
  m_acceptButton = new WPushButton( "Accept", m_footer );
  m_acceptButton->setFloatSide(Wt::Right);
  m_acceptButton->addStyleClass( "AcceptIcon" );
  HelpSystem::attachToolTipOn( m_acceptButton,"Accept all changes/selections made and close dialog" , showToolTipInstantly );

  m_acceptButton->clicked().connect( this, &DetectorEdit::acceptAndFinish );
  if( !currentDet || !currentDet->isValid() )
    m_acceptButton->disable();
  
  m_group->setSelectedButtonIndex( 0 );
  
  init();
}//DetectorEdit constructor


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
  }//try
  catch(...)
  {
    failed = true;
    passMessage( "Error getting from DetectorPeakResponse table", "DetectorEdit", WarningWidget::WarningMsgHigh );
  }//try / catch
  
  if( failed )
    m_acceptButton->disable();
  else
    m_acceptButton->enable();
  
  m_detector = det;
  emitChangedSignal();
}//void dbTableSelectionChanged()

void DetectorEdit::init()
{
  //----initialize
  
  if( !m_detector )
  {
//    m_group->setSelectedButtonIndex(0);
  }else
  {
    switch( m_detector->efficiencySource() )
    {
      case DetectorPeakResponse::kGadrasEfficiencyDefintion:
      {
        const bool selected = m_gadrasDetSelect->trySelectDetector( m_detector );
        
        if( selected )
          m_group->setSelectedButtonIndex( 0 );
        else
          m_group->setSelectedButtonIndex( 2 );
        break;
      }//case kGadrasEfficiencyDefintion:
        
      case DetectorPeakResponse::kRelativeEfficiencyDefintion:
      {
        m_group->setSelectedButtonIndex( 1 );
        m_relEffSelect->trySelectDetector( m_detector );
        break;
      }//case kRelativeEfficiencyDefintion:
        
        
      case DetectorPeakResponse::kUserUploadedEfficiencyCsv:
      {
        m_group->setSelectedButtonIndex( 2 );
        break;
      }//case kUserUploadedEfficiencyCsv:
        
      case DetectorPeakResponse::kUserEfficiencyEquationSpecified:
      {
        m_group->setSelectedButtonIndex( 3 );
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
        
      case DetectorPeakResponse::kUnknownEfficiencySource:
        break;
    }//switch( m_detector->efficiencySource() )
  
  }//if( m_detector == NULL ) / else
  
  selectButton( m_stack, m_group, false );
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
  m_acceptButton->disable();
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
//    m_acceptButton->enable();
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
    
    detec->setIntrinsicEfficiencyFormula( fcn, det_diam, energyUnits );
    
    updateChart();
  }catch( std::exception &e )
  {
    passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    m_acceptButton->disable();
    m_manualSetButton->enable();
    updateChart();
    return;
  }//try / catch
  
  m_acceptButton->enable();
  m_manualSetButton->disable();
  m_detector = detec;

  emitChangedSignal();
} // setDefineDetector(const std::string fcn)



// This calls the action to initialize the charts to whatever
// is selected.  Each stack page should have something to update
void DetectorEdit::selectButton( WStackedWidget *stack,
                                 WButtonGroup *group,
                                 bool activateCallBack)
{
  
  //makes one of the stack visible
  stack->setCurrentIndex(group->checkedId());
  
  if( activateCallBack )
  {
    
    //clear the chart when we change, regardless
    std::shared_ptr<DetectorPeakResponse> det;
    m_detector = det;
    
    switch( group->checkedId() )
    {
      case 0:
        gadrasDetectorSelectCallback();
        break;
    
      case 1:
        relEffDetectorSelectCallback();
        break;
        
      case 2:
        fileUploadedCallback();
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

		//capping it at 1.2 (just outside of rendering range), so values don't crap out Wt
        if( efficincy>1.2 )
           efficincy = 1.2f; 
          
        m_efficiencyModel->setData(row, 1, efficincy );
        if( hasResloution )
          m_efficiencyModel->setData(row, 2, m_detector->peakResolutionFWHM(energy));
      }//for( int row = 0; row < numEnergyPoints; ++row )
     
      m_efficiencyModel->setHeaderData(0, Wt::WString("Energy"));
      m_efficiencyModel->setHeaderData(1, Wt::WString("Efficiency"));
      if( m_detector->hasResolutionInfo() )
        m_efficiencyModel->setHeaderData(2, Wt::WString("FWHM"));
    
      m_chart->setXSeriesColumn(0);  //Not having this line after creating a new model was why the x-axis was the wrong scale
      Wt::Chart::WDataSeries s1(1, Wt::Chart::LineSeries,Wt::Chart::Y1Axis);
      m_chart->addSeries(s1);
    
      //m_chart->axis(Wt::Chart::Y1Axis).setRange( 0.0, 1.0 );
      //m_chart->axis(Wt::Chart::Y1Axis).setRoundLimits(<#WFlags<Wt::Chart::AxisValue> locations#>)
      m_chart->axis(Wt::Chart::Y1Axis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
    
      if( hasResloution )
      { //only if there is resolution FWHM
        Wt::Chart::WDataSeries s2(2, Wt::Chart::LineSeries,Wt::Chart::Y2Axis);
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
      cerr << SRC_LOCATION << "\n\tCaught: " << e.what() << endl;
    }//try / catch
  }//if( hasEfficiency )
} //DetectorEdit::updateChart()


void DetectorEdit::fileUploadedCallback()
{
  if( !m_efficiencyCsvUpload->empty() )
    m_detectrDiameterDiv->show();
  else
    m_detectrDiameterDiv->hide();


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
  }catch(...){}

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

  ifstream csvfile( csvFileName.c_str(), ios_base::binary|ios_base::in );
  if( !csvfile.is_open() )
    return;

  std::shared_ptr<DetectorPeakResponse> det( new DetectorPeakResponse() );
  try
  {
    if( isGadrasDet )
    {
      const string dotDatFileName = m_detectorDotDatUpload->spoolFileName();
      ifstream datfile( dotDatFileName.c_str(), ios_base::binary|ios_base::in );
      if( !datfile.is_open() )
        return;
      det->fromGadrasDefinition( csvfile, datfile );
    }else
    {
      det->fromEnergyEfficiencyCsv( csvfile, diameter, float(PhysicalUnits::keV) );
    }//if( isGadrasDet ) / else
  }catch( std::exception &e )
  {
    m_acceptButton->disable();
    passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    return;
  }//try / catch
  
  string csvOrigName = m_efficiencyCsvUpload->clientFileName().toUTF8();
  if( UtilityFunctions::iends_with( csvOrigName, ".csv" ) )
    csvOrigName = csvOrigName.substr(0, csvOrigName.size()-4);
  if( det )
    det->setName( csvOrigName );

  m_acceptButton->enable();

  m_detector = det;
  emitChangedSignal();
}//void fileUploadedCallback();


void DetectorEdit::relEffDetectorSelectCallback()
{
  m_relEffSelect->userSelectedRelEffDetSelect();
}//void relEffDetectorSelectCallback()


void DetectorEdit::gadrasDetectorSelectCallback()
{
  std::shared_ptr<DetectorPeakResponse> det = m_gadrasDetSelect->selectedDetector();

  if( !det )
    m_acceptButton->disable();
  else
    m_acceptButton->enable();

  m_detector = det;
  emitChangedSignal();
}//void gadrasDetectorSelectCallback()


void DetectorEdit::acceptAndFinish()
{
  //Save detector to DB
  try
  {
    if( !! m_detector )
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      
      //Create a separate DetectorPeakResponse because shared_ptr and dbo::ptr don't work well together
      DetectorPeakResponse* tempDetector = new DetectorPeakResponse( *m_detector );
      tempDetector->m_user = m_interspec->m_user.id(); //adds the current user to the detectorpeakresponse insertion into DB
      m_sql->session()->add( tempDetector );
      transaction.commit();
    } //session
  }catch( std::exception &e )
  {
    passMessage( "Error writing DetectorPeakResponse to DB: " + std::string(e.what()), "DetectorEdit", WarningWidget::WarningMsgHigh );
  }//try / catch
  emitChangedSignal();
  emitModifiedSignal();
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


vector<pair<string,string>> DetectorEdit::avaliableGadrasDetectors() const
{
  vector<pair<string,string>> answer;

  //Look through paths specified in user prefernce 'GadrasDRFPath' and find
  //  directoires containing both a 'Detector.dat' and 'Efficiency.csv' file.
  
#if( BUILD_FOR_WEB_DEPLOYMENT )
  const string datadir = InterSpec::dataDirectory();
  const string drfpaths = UtilityFunctions::append_path( datadir, "GenericGadrasDetectors" )
                          + UtilityFunctions::append_path( datadir, "GadrasDetectors" );
#else
  const string drfpaths = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", m_interspec );
#endif
  
  vector<string> paths;
  UtilityFunctions::split( paths, drfpaths, "\r\n;" );
  
  for( std::string path : paths )
  {
    UtilityFunctions::trim( path );
    const vector<string> csvs = UtilityFunctions::recursive_ls( path, "Efficiency.csv" ); //ending string is case-insensitive
    
    for( const std::string &csv : csvs )
    {
      const string parent = UtilityFunctions::parent_path(csv);
      vector<string> files = UtilityFunctions::ls_files_in_directory(parent);
      bool found_dat = false;
      for( size_t i = 0; !found_dat && i < files.size(); ++i )
        found_dat = UtilityFunctions::iends_with( files[i], "Detector.dat" );
      if( found_dat )
      {
        try
        {
          //fs_relative probably shouldnt throw, but JIC.
          //fs_relative( const std::string &from_path, const std::string &to_path )
          string displayname = UtilityFunctions::fs_relative( path, parent );
          if( displayname.empty() || displayname == "." )
            displayname = UtilityFunctions::filename(parent);
          
          answer.push_back( make_pair(parent, displayname) );
        }catch(...)
        {
        }
      }//if( found_dat )
    }
  }//for( const std::string &path : paths )

  return answer;
}//vector<string> avaliableGadrasDetectors() const


std::shared_ptr<DetectorPeakResponse> DetectorEdit::initARelEffDetector( const int typeint, InterSpec *interspec )
{
  string smname;
  switch( typeint )
  {
    case kGR135Detector:             smname = "GR135";             break;
    case kIdentiFinderDetector:      smname = "IdentiFINDER";      break;
    case kIdentiFinderNGDetector:    smname = "IdentiFINDER-NGH";  break;
    case kIdentiFinderLaBr3Detector: smname = "IdentiFINDER-LaBr"; break;
    case kDetectiveDetector:         smname = "Detective";         break;
    case kDetectiveExDetector:       smname = "Detetive DX";       break;
    case kDetectiveEx100Detector:    smname = "Detective EX-100";  break;
    case kOrtecIDMPortalDetector:    smname = "Detective EX-200";  break;  //Havent actually seen anywhere
    case kSAIC8Detector:             smname = "";                  break;
    case kFalcon5000:                smname = "Falcon 5000";       break;
    case kUnknownDetector:           smname = "";                  break;
    case kMicroDetectiveDetector:    smname = "Micro Detective";   break;
  }//switch( type )
  
  if( smname.empty() )
    throw std::runtime_error( "Could not find Rel. Eff. detector response functions" );
  
  string concat_filenames;
#if( BUILD_FOR_WEB_DEPLOYMENT )
  concat_filenames = UtilityFunctions::append_path( InterSpec::dataDirectory(), "lanl_simplemass_detectors.tsv" );
#else
  if( interspec )
    concat_filenames = InterSpecUser::preferenceValue<string>( "RelativeEffDRFPaths", interspec );
#endif
  
  vector<string> paths;
  UtilityFunctions::split( paths, concat_filenames, "\r\n;" );
  
  for( const string &filename : paths )
  {
    std::ifstream input( filename.c_str() );
  
    if( !input.is_open() )
      continue;
  
    string line;
    while( UtilityFunctions::safe_get_line( input, line ) )
    {
      UtilityFunctions::trim( line );
      if( line.empty() || line[0]=='#' || line.size() < smname.size() )
        continue;
      
      if( UtilityFunctions::icontains( line, smname ) )
      {
        auto det = RelEffFile::parseDetector( line );
        if( det )
          return det;
      }//if( thisname == smname )
    }//while( UtilityFunctions::safe_get_line( input, line ) )
  }//for( const string &filename : paths )
  
  throw runtime_error( "Coudlnt find detector '" + smname + " in relative efficiency files" );
  
  return std::shared_ptr<DetectorPeakResponse>();
}//initARelEffDetector( int type )

std::shared_ptr<DetectorPeakResponse> DetectorEdit::initAGadrasDetector(
                                                            const int typeint, InterSpec *interspec )

{
  string name;
  switch( typeint )
  {
    case kGR135Detector:             name = "GR135";              break;
    case kIdentiFinderDetector:      name = "identiFINDER-N";     break;
    case kIdentiFinderNGDetector:    name = "identiFINDER-NGH";   break;
    case kIdentiFinderLaBr3Detector: name = "identiFINDER-LaBr3"; break;
    case kDetectiveDetector:         name = "Detective";          break;
    case kDetectiveExDetector:       name = "Detective-EX";       break;
    case kDetectiveEx100Detector:    name = "Detective-EX100";    break;
    case kOrtecIDMPortalDetector:    name = "Ortec IDM Portal";   break;
    case kSAIC8Detector:             name = "";                   break;
    case kFalcon5000:                name = "Falcon 5000";        break;
    case kUnknownDetector:           name = "";                   break;
    case kMicroDetectiveDetector:    name = "MicroDetective";     break;
      
    //case kRadHunterNaI:            name = ""; break;
    //case kRadHunterLaBr3:          name = ""; break;
    case kRsi701:                    name = "NaI 2x4x16"; break;
    case kRsi705:                    name = "NaI 2x4x16"; break;
    case kAvidRsi:                   name = "NaI 2x4x16"; break;
    case kOrtecRadEagleNai:          name = "RadEagle"; break;
    //case kOrtecRadEagleCeBr2Inch: name = ""; break;
    //case kOrtecRadEagleCeBr3Inch: name = ""; break;
    //case kOrtecRadEagleLaBr: name = ""; break;
    //case kSam940LaBr3: name = ""; break;
    case kSam940:                     name = "SAM-945"; break;
    case kSam945:                     name = "SAM-945"; break;
  }//switch( type )
  
  if( name.empty() )
    throw runtime_error( "There is no GADRAS detector response function for a "
                         + detectorTypeToString( DetectorType(typeint) ) );
  
  return initAGadrasDetector( name, interspec );
}//initAGadrasDetector


std::shared_ptr<DetectorPeakResponse> DetectorEdit::initAGadrasDetector( const std::string &currentDetName, InterSpec *interspec )
{
  //Grab "GadrasDRFPath" and split it by semicolon and newlines, then go through
  //  and look for sub-folders with the given name that contain Detector.data
  //  and Efficiency.csv.
#if( BUILD_FOR_WEB_DEPLOYMENT )
  const string datadir = InterSpec::dataDirectory();
  const string drfpaths = UtilityFunctions::append_path( datadir, "GenericGadrasDetectors" )
                          + UtilityFunctions::append_path( datadir, "GadrasDetectors" );
#else
  const string drfpaths = InterSpecUser::preferenceValue<string>( "GadrasDRFPath", interspec );
#endif
  
  
  vector<string> paths;
  UtilityFunctions::split( paths, drfpaths, "\r\n;" );
  for( string basepath : paths )
  {
    UtilityFunctions::trim( basepath );
    const string path = UtilityFunctions::append_path(basepath, currentDetName);
    if( !UtilityFunctions::is_directory(path) )
      continue;
    
    string thiscsv, thisdat;
    const vector<string> files = UtilityFunctions::ls_files_in_directory( path );
    for( size_t i = 0; (thiscsv.empty() || thisdat.empty()) && i < files.size(); ++i )
    {
      //We will be loose with upper/lower case
      if( UtilityFunctions::iends_with(files[i], "Efficiency.csv") )
        thiscsv = files[i];
      else if( UtilityFunctions::iends_with(files[i], "Detector.dat") )
        thisdat = files[i];
    }//for( size_t i = 0; i < files.size(); ++i )
    
    if( thiscsv.size() && thisdat.size() )
      return initAGadrasDetectorFromDirectory( path, interspec );
  }//for( const string &basepath : paths )
  
  throw runtime_error( "Could not find GADRAS detector names " + currentDetName );
}//std::shared_ptr<DetectorPeakResponse> initAGadrasDetector( const std::string &name );


std::shared_ptr<DetectorPeakResponse> DetectorEdit::initAGadrasDetectorFromDirectory( const std::string &path, InterSpec *interspec )
{
  auto det = std::make_shared<DetectorPeakResponse>();

  if( !UtilityFunctions::is_directory(path) )
    throw runtime_error( "'" + path + "' is not a directory" );
  
  //Lets be lienient on upper vs lower case.
  string thiscsv, thisdat;
  const vector<string> files = UtilityFunctions::ls_files_in_directory( path );
  for( size_t i = 0; (thiscsv.empty() || thisdat.empty()) && i < files.size(); ++i )
  {
    //We will be loose with upper/lower case
    if( UtilityFunctions::iends_with(files[i], "Efficiency.csv") )
      thiscsv = files[i];
    else if( UtilityFunctions::iends_with(files[i], "Detector.dat") )
      thisdat = files[i];
  }//for( size_t i = 0; i < files.size(); ++i )
  
  if( thiscsv.empty() )
    throw runtime_error( "Could not find a Efficiency.csv file in '" + path + "'" );
  
  if( thisdat.empty() )
    throw runtime_error( "Could not find a Detector.dat file in '" + path + "'" );
  
  ifstream effstrm( thiscsv.c_str(), ios_base::binary|ios_base::in );
  if( !effstrm.is_open() )
    throw runtime_error( "Could not open '" + thiscsv + "'" );

  ifstream datstrm( thisdat.c_str(), ios_base::binary|ios_base::in );
  if( !datstrm.is_open() )
    throw runtime_error( "Could not open '" + thisdat + "'" );

  det->fromGadrasDefinition( effstrm, datstrm );
  det->setName( UtilityFunctions::filename( UtilityFunctions::parent_path(thiscsv) ) );

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
  m_acceptButton->setDisabled( !m_detector );
  
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
  : AuxWindow("Detector Edit/Select", true),
    m_edit( NULL )
{
  disableCollapse();
  m_edit = new DetectorEdit( det, specViewer, fileModel, this );
  m_edit->done().connect( boost::bind( &DetectorEditWindow::acceptAndDelete, this ) );
  finished().connect( boost::bind( &DetectorEditWindow::acceptAndDelete, this ) );
  rejectWhenEscapePressed();
  setMargin(0);
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
