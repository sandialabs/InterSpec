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
#include <sstream>
#include <string>
#include <vector>

#include <boost/regex.hpp>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WSlider>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WTabWidget>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WDoubleValidator>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include "InterSpec/SpecMeas.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/CompactFileManager.h"
#if( !ANDROID && !IOS )
#include "InterSpec/FileDragUploadResource.h"
#endif

using namespace Wt;
using namespace std;


CompactFileManager::CompactFileManager( SpecMeasManager *fileManager,
                                        InterSpec *hostViewer,
                                        CompactFileManager::DisplayMode mode,
                                        WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_selects{ nullptr },
    m_displayedPreTexts{ nullptr },
    m_displayedPostTexts{ nullptr },
    m_sampleDivs{ nullptr },
    m_displaySampleNumEdits{ nullptr },
    m_nextSampleNumButtons{ nullptr },
    m_prevSampleNumButtons{ nullptr },
    m_scaleValueTxt{ nullptr },
    m_rescaleByLiveTime{ nullptr },
    m_previousSampleTxt{ WString() },
    m_titles{ nullptr },
    m_files( fileManager->model() ),
    m_interspec( hostViewer ),
    m_fullFileManager( fileManager ),
    m_displayMode( mode )
{
  wApp->useStyleSheet( "InterSpec_resources/CompactFileManager.css" );
  
  addStyleClass( "CompactFileManager" );

#if (USE_DB_TO_STORE_SPECTRA)
  std::shared_ptr<DataBaseUtils::DbSession> sql = hostViewer->sql();
#endif
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", hostViewer );
  
  WTabWidget *tabWidget = nullptr;
  switch( m_displayMode )
  {
    case LeftToRight:
      // We will use CSS Grid layout to position the Foreground/Background/Secondary portions
      addStyleClass( "LeftToRight" );
    break;
      
    case Tabbed:
      // We will use a tabbed widget to hold Foreground/Background/Secondary seperately
      addStyleClass( "Tabbed" );
      tabWidget = new WTabWidget( this );
    break;
  }//switch( m_displayMode )
  
  const char *val_regex = "(\\-?\\d+\\s*((\\-|to|through)\\s*\\d+)?,*\\s*)+";
  WRegExpValidator *validator = new WRegExpValidator( val_regex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  
  for( int j = 0; j < 3; ++j )
  {
    WContainerWidget *wrapper = new WContainerWidget();
    wrapper->addStyleClass( "DispType" );
    
    switch( m_displayMode )
    {
      case LeftToRight:
        addWidget( wrapper );
        break;
        
      case Tabbed:
      {
        const char *tablabel = 0;
        switch( j )
        {
          case 0: tablabel = "Foreground"; break;
          case 1: tablabel = "Background"; break;
          case 2: tablabel = "Second";     break;
        }//switch( i )
        
        tabWidget->addTab( wrapper, tablabel, WTabWidget::PreLoading );
        break;
      }//case Tabbed:
    }//switch( m_displayMode )
    
    SpecUtils::SpectrumType type;
    WLabel *label = NULL;

    switch( j )
    {
      case 0:
        type = SpecUtils::SpectrumType::Foreground;
        wrapper->addStyleClass( "Foreground" );
        label = new WLabel( "Foreground:", wrapper );
      break;
      
      case 1:
        type = SpecUtils::SpectrumType::Background;
        wrapper->addStyleClass( "Background" );
        label = new WLabel( "Background:", wrapper );
      break;
      
      case 2:
        type = SpecUtils::SpectrumType::SecondForeground;
        wrapper->addStyleClass( "Secondary" );
        label = new WLabel( "Second Foreground:", wrapper );
      break;
    }//switch( i )
    
    const int typeindex = static_cast<int>(type);
    
    label->setInline( false );
    
    m_selects[typeindex] = new WComboBox( wrapper );
    m_selects[typeindex]->setInline( false );
    m_selects[typeindex]->addStyleClass( "SpecFileSelect" );
    m_selects[typeindex]->setCurrentIndex( -1 );
    m_selects[typeindex]->activated().connect(
                  boost::bind( &CompactFileManager::handleFileChangeRequest, this, _1, type ) );
    
    m_sampleDivs[typeindex] = new WContainerWidget( wrapper );
    m_sampleDivs[typeindex]->setStyleClass( "SampleSelectRow" );

    m_prevSampleNumButtons[typeindex] = new WContainerWidget( m_sampleDivs[typeindex] );
    m_prevSampleNumButtons[typeindex]->setStyleClass( "PrevSampleNum" );
    
    m_prevSampleNumButtons[typeindex]->clicked().connect(
                 boost::bind( &CompactFileManager::handleUserIncrementSampleNum,
                              this, type, false) );
    
    m_displayedPreTexts[typeindex] = new WText( m_sampleDivs[typeindex] );
    WLineEdit *edit = new WLineEdit( m_sampleDivs[typeindex] );
    edit->addStyleClass( "SampleNumInput" );
    edit->setValidator( validator );
    edit->setAutoComplete( false );
    edit->setTextSize( 6 );
    
    const char *tooltip = "Enter the sample number you'de like displayed here"
                           ". You may enter a range of sample numbers similar"
                           " to '33-39' or '33 to 39'; or CSV sample numbers"
                           " like '34,39,84'";
    HelpSystem::attachToolTipOn( edit, tooltip, showToolTips, HelpSystem::ToolTipPosition::Bottom );
    
    m_displaySampleNumEdits[typeindex] = edit;
    

    m_displayedPostTexts[typeindex] = new WText( m_sampleDivs[typeindex] );
    m_nextSampleNumButtons[typeindex] = new WContainerWidget( m_sampleDivs[typeindex] );
    m_nextSampleNumButtons[typeindex]->addStyleClass( "NextNextSample" );
    
    
    m_nextSampleNumButtons[typeindex]->clicked().connect(
                 boost::bind( &CompactFileManager::handleUserIncrementSampleNum,
                              this, type, true) );
    
    m_displaySampleNumEdits[typeindex]->enterPressed().connect(
                 boost::bind( &CompactFileManager::handleUserChangeSampleNum,
                              this, type ) );
    
    m_displaySampleNumEdits[typeindex]->blurred().connect(
                 boost::bind( &CompactFileManager::handleSampleNumEditBlur,
                              this, type ) );
    
    m_sampleDivs[typeindex]->setHiddenKeepsGeometry( true );


    m_titles[typeindex] = new WText( "", wrapper );
    m_titles[typeindex]->addStyleClass( "SpecTitle" );
    m_titles[typeindex]->setInline( false );
    m_titles[typeindex]->hide();
    
    auto bottomrow = new WContainerWidget( wrapper );
    bottomrow->addStyleClass( "BottomRow" );
    
    if( type == SpecUtils::SpectrumType::Foreground )
    {
      m_scaleValueTxt[typeindex] = nullptr;
      m_rescaleByLiveTime[typeindex] = nullptr;
    }else
    {
      auto scalerow = new WContainerWidget( bottomrow );
      scalerow->addStyleClass( "ScaleFactorRow" );
      
      WLabel *label = new WLabel( "Scale Factor:&nbsp;", scalerow );
      m_scaleValueTxt[typeindex] = new NativeFloatSpinBox( scalerow );
      label->setBuddy( m_scaleValueTxt[typeindex] );
      //m_scaleValueTxt[typeindex]->setSingleStep(0.1);
      m_scaleValueTxt[typeindex]->setSpinnerHidden( true );
      m_scaleValueTxt[typeindex]->setRange( 0.0, 1000000.0 );
      m_scaleValueTxt[typeindex]->addStyleClass( "SpecNormTxt" );
      m_scaleValueTxt[typeindex]->setFormatString( "%.4G" );
      m_scaleValueTxt[typeindex]->valueChanged().connect( boost::bind( &CompactFileManager::handleUserEnterdScaleFactor, this, type) );
      m_scaleValueTxt[typeindex]->mouseWheel().connect( boost::bind( &CompactFileManager::handleUserEnterdScaleFactorWheel, this, type, _1) );
      
      HelpSystem::attachToolTipOn( m_scaleValueTxt[typeindex],
                                  "Factor to scale the spectrum by. Entering an"
                                  " empty string will cause spectrum to be live"
                                  " time normalized to foreground spectrum.",
                                  showToolTips );
      
      
      m_rescaleByLiveTime[typeindex] = new WPushButton( "Normalize", scalerow );
      m_rescaleByLiveTime[typeindex]->hide();
      m_rescaleByLiveTime[typeindex]->clicked().connect( boost::bind( &CompactFileManager::handleRenormalizeByLIveTime, this, type) );
      m_scaleValueTxt[typeindex]->disable();
    }//if( type == SpecUtils::SpectrumType::Foreground ) / else
  }//for( int j = 0; j < 3; ++j )
  
  
  
  //Lets add in a few more customizations based on the display type
  switch( m_displayMode )
  {
    case LeftToRight:
    {
      //Foreground, add in Manager and Library buttons for quick access
      WContainerWidget *buttons = new WContainerWidget( this );
      buttons->addStyleClass( "CompactManagerButtons" );
      
      auto helpBtn = new WContainerWidget( buttons );
      helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
      helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "compact-file-manager" ) );
      
      WPushButton *button = new WPushButton(  "Spectrum Manager...", buttons );
      //button->setIcon(Wt::WLink("InterSpec_resources/images/computer.png" ));
      button->clicked().connect( m_interspec->fileManager(), &SpecMeasManager::startSpectrumManager );
#if( USE_DB_TO_STORE_SPECTRA )
      WPushButton *button2 = new WPushButton( "Previous...", buttons );
      button2->clicked().connect( boost::bind( &SpecMeasManager::browseDatabaseSpectrumFiles,
                                              m_interspec->fileManager(),
                                              SpecUtils::SpectrumType::Foreground ) );
#endif
    break;
    }//case LeftToRight:
      
    case Tabbed:
    {
      //Add in a link to open files in the database, as
#if( USE_DB_TO_STORE_SPECTRA )
      WContainerWidget *buttons = new WContainerWidget( this );
      WPushButton *button = new WPushButton( "Prev. Saved Spectra", buttons );
        button->clicked().connect( boost::bind( &SpecMeasManager::browseDatabaseSpectrumFiles,
                                               fileManager,
                                               SpecUtils::SpectrumType::Foreground ) );
#endif
      break;
    }//case Tabbed:
  }//switch( m_displayMode )
    
  // Then actually pull in the available files. If there are none, say so.
  refreshContents();

  // Lastly, hook up the update handler so it'll keep it posted.
  m_files->rowsInserted().connect( this, &CompactFileManager::refreshContents );
  m_files->rowsRemoved().connect( this, &CompactFileManager::refreshContents );
  hostViewer->displayedSpectrumChanged().connect(
                          boost::bind( &CompactFileManager::handleDisplayChange,
                                       this, _1, _2, _3, _4 ) );
}//CompactFileManager constructor


CompactFileManager::~CompactFileManager()
{
  // no-op
}

void CompactFileManager::handleFileChangeRequest( int row, SpecUtils::SpectrumType type )
{
  const int typeindex = static_cast<int>( type );
  if( row < 0 || row >= m_files->rowCount() )
  {
    m_displayedPreTexts[typeindex]->setText( "" );
    m_displayedPostTexts[typeindex]->setText( "" );
    m_displaySampleNumEdits[typeindex]->setText( "" );
    std::shared_ptr<SpecMeas> blankptr;
    m_fullFileManager->displayFile( -1, blankptr, type, false, true, SpecMeasManager::VariantChecksToDo::None );
    if( m_titles[typeindex] )
    {
      m_titles[typeindex]->setText("");
      m_titles[typeindex]->hide();
    }//if( type == SpecUtils::SpectrumType::Foreground )
    return;
  }//if( we should unload current file )

  std::shared_ptr<const SpectraFileHeader> header = m_files->fileHeader(row);
  std::shared_ptr<SpecMeas> meas = header->parseFile();

  m_fullFileManager->displayFile( row, meas, type, false, true, SpecMeasManager::VariantChecksToDo::None );
}//void handleFileChangeRequest( SpecUtils::SpectrumType type );

/**
 Refactored static method so can be called outside of CompactFileManager
 
 cfm - optional, if provided, will call handleDisplayChange.  Set to NULL normally if called statically.  Otherwise, CompactFileManager calls will provide it.
 */
void CompactFileManager::changeToSampleNum( int sampleNum,
                                            SpecUtils::SpectrumType type,
                                            InterSpec *viewer,
                                            CompactFileManager* cfm )
{
  std::shared_ptr<SpecMeas> meas = viewer->measurment( type );

  if( meas )
  {
    const set<int> total_sample_nums = meas->sample_numbers();
    if( total_sample_nums.find(sampleNum) == total_sample_nums.end() )
    {
      if( cfm )
      {
        const auto dets = viewer->detectorsToDisplay(type);
        const auto &samples = viewer->displayedSamples(type);
        cfm->handleDisplayChange( type, meas, samples, dets );
      }
      
      passMessage( "Sample number requested doesnt exist.",
                   "", WarningWidget::WarningMsgHigh );
      return;
    }//if( total_sample_nums.find(sampleNum) == total_sample_nums.end() )

    
    set<int> sampleNumToLoad;
    sampleNumToLoad.insert( sampleNum );
    viewer->changeDisplayedSampleNums( sampleNumToLoad, type );
//    m_interspec->setSpectrum( meas, sampleNumToLoad, type, false );
  }else
  {
    passMessage( "Cant change sample numbers, there is no current measurment",
                 "", WarningWidget::WarningMsgHigh );
  }//if( meas ) / else
}//void changeToSampleNum(...)


void CompactFileManager::changeToSampleRange( int first, int last,
                                              SpecUtils::SpectrumType type )
{
  std::shared_ptr<SpecMeas> meas = m_interspec->measurment( type );

  if( meas )
  {
    const set<int> total_sample_nums = meas->sample_numbers();
    const set<int>::const_iterator firstPos = total_sample_nums.find(first);
    set<int>::const_iterator lastPos = total_sample_nums.find(last);

    if( firstPos == total_sample_nums.end()
        || lastPos == total_sample_nums.end() )
    {
      const auto dets = m_interspec->detectorsToDisplay(type);
      const auto &samples = m_interspec->displayedSamples(type);
      
      handleDisplayChange( type, meas, samples, dets );
      passMessage( "One of the sample number requested doesnt exist.",
                   "", WarningWidget::WarningMsgHigh );
      return;
    }//if( total_sample_nums.find(sampleNum) == total_sample_nums.end() )

    set<int> sampleNumToLoad;
    if( lastPos != total_sample_nums.end() )
      ++lastPos;
    sampleNumToLoad.insert( firstPos, lastPos );
    m_interspec->changeDisplayedSampleNums( sampleNumToLoad, type );
//    m_interspec->setSpectrum( meas, sampleNumToLoad, type, false );
  }else
  {
    passMessage( "Can't change sample numbers, there is no current measurment",
                 "",WarningWidget::WarningMsgHigh );
  }//if( meas ) / else
}//void changeToSampleRange(...)


void CompactFileManager::handleSampleNumEditBlur( SpecUtils::SpectrumType type )
{
  const int typeindex = static_cast<int>( type );
  if( m_previousSampleTxt[typeindex] == m_displaySampleNumEdits[typeindex]->text() )
    return;

  handleUserChangeSampleNum( type );
}//void handleSampleNumEditBlur( SpecUtils::SpectrumType spectrum_type )


void CompactFileManager::handleSpectrumScale( const double scale, SpecUtils::SpectrumType type )
{
  const int typeindex = static_cast<int>( type );
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      cerr << "CompactFileManager::handleSpectrumScale: cannot handle scaling of foreground - ignoring" << endl;
      return;
      
    case SpecUtils::SpectrumType::Background:
    case SpecUtils::SpectrumType::SecondForeground:
      m_interspec->setDisplayScaleFactor( scale, type );
      updateDisplayedScaleFactorNumbers( scale, type );
      if( m_rescaleByLiveTime[typeindex] )
        m_rescaleByLiveTime[typeindex]->show();
      break;
  }//switch( type )
}//void handleSpectrumScale( const double scale, SpecUtils::SpectrumType spectrum_type );


void CompactFileManager::handleUserChangeSampleNum( SpecUtils::SpectrumType type )
{
  const int typeindex = static_cast<int>( type );
  const WString text = m_displaySampleNumEdits[typeindex]->text();

  string fulltxt = text.toUTF8();

  //Lets replace numbers seperated by spaces, to be seperated by commas.
  //  We cant do a simple replace of spaces to commas, and using a regex would
  //  require using a lookahead or behind, and I dont think boost supports that
  //  always.  Note that below while loop is a little ineficient, but whatever
  boost::smatch mtch;
  boost::regex expr( ".*(\\d\\s+\\d).*" );
  while( boost::regex_match( fulltxt, mtch, expr ) )
    fulltxt[mtch[1].first - fulltxt.begin() + 1] = ',';
  
//Using a sregex_token_iterator doesnt work for single digit numbers, ex
//  '1 2 3 45 983 193'  will go to '1,2 3,45,983,193'
//  boost::regex expr( "\\d\\s+\\d" );
//  for( boost::sregex_token_iterator iter(fulltxt.begin(),fulltxt.end(),expr,0);
//      iter != boost::sregex_token_iterator(); ++iter )
//    fulltxt[iter->first - fulltxt.begin()+1] = ',';
  
  vector<string> sampleranges;
  SpecUtils::split( sampleranges, fulltxt, "," );

  std::shared_ptr<SpecMeas> meas = m_interspec->measurment( type );

  if( !meas )
  {
    passMessage( "Can't change sample numbers, there is no current measurment",
                 "", WarningWidget::WarningMsgHigh );
    m_displaySampleNumEdits[typeindex]->setText( "" );
    return;
  }//if( !meas )

  const set<int> samples = meas->sample_numbers();
  set<int> sampleNumToLoad;

  try
  {
    for( string textStr : sampleranges )
    {
      SpecUtils::trim( textStr );
      if( textStr.empty() )
        continue;

      boost::smatch matches;
      boost::regex range_expression( "(\\d+)\\s*(\\-|to|through)\\s*(\\d+)",
                                     boost::regex::perl|boost::regex::icase );
      if( boost::regex_match( textStr, matches, range_expression ) )
      {
        string firstStr = string( matches[1].first, matches[1].second );
        string lastStr = string( matches[3].first, matches[3].second );

        int first = std::stoi( firstStr );
        int last = std::stoi( lastStr );
        if( last < first )
          std::swap( last, first );

        const set<int>::const_iterator firstPos = samples.find(first);
        set<int>::const_iterator lastPos = samples.find(last);

        if( firstPos==samples.end() || lastPos==samples.end() )
        {
          stringstream msg;
          if( firstPos == samples.end() )
            msg << "Sample number " << firstStr << " doesnt exist.";
          else
            msg << "Sample number " << lastStr << " doesnt exist.";
          passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
          const auto dets = m_interspec->detectorsToDisplay(type);
          const set<int> &samples = m_interspec->displayedSamples(type);
          handleDisplayChange( type, meas, samples, dets );
          return;
        }//if( samples.find(sampleNum) == samples.end() )

        if( lastPos != samples.end() )
          ++lastPos;
        sampleNumToLoad.insert( firstPos, lastPos );
      }else
      {
        const int sample = std::stoi( textStr );
        
        if( !samples.count(sample) )
        {
          stringstream msg;
          msg << "Sample number " << sample << " does not exist in the measurement.";
          passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
          const auto dets = m_interspec->detectorsToDisplay(type);
          const set<int> &samples = m_interspec->displayedSamples(type);
          handleDisplayChange( type, meas, samples, dets );
          return;
        }//if( !samples.count(sample) )
          
//        changeToSampleNum( sample, type );
        sampleNumToLoad.insert( sample );
      }//if( is sample range ) / else
    }//for( string textStr : sampleranges )

//    m_interspec->setSpectrum( meas, sampleNumToLoad, type, false );
    m_interspec->changeDisplayedSampleNums( sampleNumToLoad, type );
    updateDisplayedScaleFactorNumbers( m_interspec->displayScaleFactor(type), type );
  }catch( exception &e )
  { 
    cerr << "CompactFileManager::handleUserChangeSampleNum( SpecUtils::SpectrumType type )" << "\n\t" << e.what() << endl;
    passMessage( "Error changing to requested sample number(s)",
                 "", WarningWidget::WarningMsgHigh );

    if( meas )
    {
      const auto dets = m_interspec->detectorsToDisplay(type);
      const auto &samples = m_interspec->displayedSamples(type);
      
      handleDisplayChange( type, meas, samples, dets );
    }
//    m_displaySampleNumEdits[type]->setText( "" );
  }//  try / catch
}//void handleUserChangeSampleNum(...)

/**
 Backward compatibility for previous CompactFileManager calls to call the refactored method
 */
void CompactFileManager::handleUserIncrementSampleNum( SpecUtils::SpectrumType type,
                                                      bool increment )
{
    handleUserIncrementSampleNum(type, increment, m_interspec, m_files, this);
}

/**
 Refactored code so can be called statically from outside of CompactFileManager
 */
void CompactFileManager::handleUserIncrementSampleNum( SpecUtils::SpectrumType type,
                                                       bool increment , InterSpec *hostViewer, SpectraFileModel *files, CompactFileManager* cfm)
{
  try
  {
    std::shared_ptr<SpecMeas> meas = hostViewer->measurment( type );

    if( !meas )
      throw std::runtime_error( "Unable to get measurment" );

    const set<int> total_sample_nums = meas->sample_numbers();
    const set<int> &currentSamples = hostViewer->displayedSamples( type );
    int currentSample = -1;

    if( currentSamples.empty() && total_sample_nums.size() )
      currentSample = *(total_sample_nums.begin());
    else if( currentSamples.size() > 1 )
      currentSample = *(currentSamples.rbegin());
    else if( currentSamples.size() )
      currentSample = *(currentSamples.begin());

    if( total_sample_nums.size() == 1 && currentSamples.size() )
      return;

    const WModelIndex index = files->index( meas );

    if( !index.isValid() )
      throw std::runtime_error( "Unable to get index" );

    set<int>::iterator pos = total_sample_nums.find( currentSample );
    if( pos == total_sample_nums.end() )
      throw std::runtime_error( "Unable to find current index" );

    if( increment )
    {
      ++pos;
    }else
    {
      if( pos == total_sample_nums.begin() )
        pos = total_sample_nums.end();
      --pos;
    }//if( increment ) / else

    if( pos == total_sample_nums.end() )
      pos = total_sample_nums.begin();

    changeToSampleNum( *pos, type , hostViewer, cfm);
  }catch( std::runtime_error e )
  {
    cerr << "CompactFileManager::handleUserIncrementSampleNum(...): caught "
         << e.what() << endl;
    passMessage( "Error advancing spectrum", "",WarningWidget::WarningMsgHigh );
  }catch(...)
  {
    passMessage( "Error advancing spectrum", "",WarningWidget::WarningMsgHigh );
  }//try / catch
}//void handleUserWantsNextSampleNum( SpecUtils::SpectrumType spectrum_type )


void CompactFileManager::handleDisplayChange( SpecUtils::SpectrumType spectrum_type,
                                  const std::shared_ptr<SpecMeas> meas,
                                  const set<int> &sample_numbers,
                                  const std::vector<std::string> &detectors )
{
  const int typeindex = static_cast<int>( spectrum_type );
  
  if( typeindex < 0 || typeindex >= 3 )
    throw runtime_error( "CompactFileManager::handleDisplayChange() - totally unexpected error" );

  if( m_titles[typeindex] )
    m_titles[typeindex]->hide();
  
  WComboBox *select = m_selects[typeindex];
  if( !meas )
  {
    select->setCurrentIndex( select->count()-1 );
    m_previousSampleTxt[typeindex] = "";
    m_displayedPreTexts[typeindex]->setText( "" );
    m_displayedPostTexts[typeindex]->setText( "" );
    m_displaySampleNumEdits[typeindex]->setText( "" );
    updateDisplayedScaleFactorNumbers( 1.0, spectrum_type );
    if( m_scaleValueTxt[typeindex] )
      m_scaleValueTxt[typeindex]->disable();
    if( m_rescaleByLiveTime[typeindex] )
      m_rescaleByLiveTime[typeindex]->hide();
    
    return;
  }//if( !meas )

  if( m_scaleValueTxt[typeindex] )
  {
    m_scaleValueTxt[typeindex]->enable();
  }
  
  if( m_rescaleByLiveTime[typeindex] )
    m_rescaleByLiveTime[typeindex]->hide();

  WModelIndex index = m_files->index( meas );
  if( !index.isValid() )
  {
    select->setCurrentIndex( select->count()-1 );
    m_previousSampleTxt[typeindex] = "";
    return;
  }//if( !index.isValid() )

  //if( meas )
  //  cout << "There are " << meas->sample_numbers().size() << " sample numbers, " << meas->measurements().size() << " meas" << endl;
  
  select->setCurrentIndex( index.row() );

  const set<int> total_sample_nums = meas->sample_numbers();

  if( sample_numbers.size() == 0 || total_sample_nums.size() == 0 )
  {
//    m_sampleDivs[typeindex]->hide();
  }else if( sample_numbers.size() == 1 )
  {
    const int lastNumber = *(total_sample_nums.rbegin());
    const int displayedNumber = *(sample_numbers.begin());

    m_displayedPostTexts[typeindex]->show();
    m_displaySampleNumEdits[typeindex]->show();
    
    const bool hideArrows = (total_sample_nums.size() < 2);
    m_nextSampleNumButtons[typeindex]->setHidden( hideArrows );
    m_prevSampleNumButtons[typeindex]->setHidden( hideArrows );
    
    
    WStringStream postMsg, editVal;
    postMsg << "/" << lastNumber; //(int)total_sample_nums.size();
    m_displayedPostTexts[typeindex]->setText( postMsg.str() );
    editVal << displayedNumber;
    m_displaySampleNumEdits[typeindex]->setText( editVal.str() );
    
    
    // Set the title, if there is one
    WString title;
    const vector<string> &dets = meas->detector_names();
    for( size_t i = 0; title.empty() && (i < dets.size()); ++i )
    {
      auto m = meas->measurement(displayedNumber, dets[i]);
      if( m && !m->title().empty() )
      {
        if( m->title().size() > 128 )  //128 chosen arbitrarily and not tested.
        {
          string t = m->title();
          SpecUtils::utf8_limit_str_size( t, 128 );
          title = WString::fromUTF8( m->title() );
        }else
        {
          title = WString::fromUTF8( m->title() );
        }
        
        
      }
    }//for( loop over detectors to find title )

    if( m_titles[typeindex] )
    {
      m_titles[typeindex]->setText( title );
      
      if( title.empty() )
      {
        if( !m_titles[typeindex]->text().empty() )
          m_titles[typeindex]->setText( "" );
      }else
      {
        m_titles[typeindex]->setText( "Title: " + title );
        m_titles[typeindex]->show();
      }
        
      if( title.narrow().length() > 70 )
        m_titles[typeindex]->setToolTip( title );
      else if( !m_titles[typeindex]->toolTip().empty() )
        m_titles[typeindex]->setToolTip( "" );
    }//if( m_titles[typeindex] )
  }else if( total_sample_nums.size() == sample_numbers.size() )
  {
    char buffer[64];
    const int totalNumber = *(total_sample_nums.rbegin());
    const int lastNumber = *(sample_numbers.rbegin());
    const int firstNumber = *(sample_numbers.begin());
    snprintf( buffer, sizeof(buffer), "/%i", totalNumber+1 );
    m_displayedPostTexts[typeindex]->setText( buffer );
    if( lastNumber == firstNumber )
      snprintf( buffer, sizeof(buffer), "%i", lastNumber );
    else if( sample_numbers.size() == 2 )
      snprintf( buffer, sizeof(buffer), "%i,%i", firstNumber, lastNumber );
    else
      snprintf( buffer, sizeof(buffer), "%i-%i", firstNumber, lastNumber );
    m_displaySampleNumEdits[typeindex]->setText( buffer );
    m_nextSampleNumButtons[typeindex]->hide();
    m_prevSampleNumButtons[typeindex]->hide();
  }else
  {
    const int lastNumber = *(total_sample_nums.rbegin());
      
    WStringStream postMsg;
    postMsg << "/" << lastNumber;
    m_displayedPostTexts[typeindex]->setText( postMsg.str() );
   
    const string displaySequence = SpecUtils::sequencesToBriefString( sample_numbers );
    const bool added = (displaySequence.find(',') != string::npos);
    m_displaySampleNumEdits[typeindex]->setText( displaySequence );
    
    m_nextSampleNumButtons[typeindex]->setHidden( !added );
    m_prevSampleNumButtons[typeindex]->setHidden( !added );
  }//else if( meas->passthrough() ){...}
  
  m_previousSampleTxt[typeindex] = m_displaySampleNumEdits[typeindex]->text();
  
  const double livetimeSF = m_interspec->displayScaleFactor(spectrum_type);
  
  updateDisplayedScaleFactorNumbers( livetimeSF, spectrum_type );
  
  if( spectrum_type == SpecUtils::SpectrumType::Foreground )
  {
    if( !!m_interspec->measurment(SpecUtils::SpectrumType::Background) )
      updateDisplayedScaleFactorNumbers( m_interspec->displayScaleFactor(SpecUtils::SpectrumType::Background), SpecUtils::SpectrumType::Background );
    if( !!m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground) )
      updateDisplayedScaleFactorNumbers( m_interspec->displayScaleFactor(SpecUtils::SpectrumType::SecondForeground), SpecUtils::SpectrumType::SecondForeground );
  }//if( spectrum_type == SpecUtils::SpectrumType::Foreground )
}//void handleDisplayChange(...)


void CompactFileManager::refreshContents()
{
  // Clear out the combo boxes
  for( int i = 0; i < 3; ++i )
    m_selects[i]->clear();

  if( m_files->rowCount() < 1 )
  {
    for( int i = 0; i < 3; ++i )
      m_selects[i]->addItem( "No uploaded/available spectra." );
  }else
  {
    for( int i = 0; i < m_files->rowCount(); ++i )
    {
      const WString name = m_files->fileHeader( i )->displayName();
      const WString time = m_files->fileHeader(i)->spectrumTime().toString("MMM d, yyyy");
      for( int i = 0; i < 3; ++i )
        m_selects[i]->addItem( name + " (" + time + ")" );
    }//for( loop over available files )

    m_selects[static_cast<int>(SpecUtils::SpectrumType::Foreground)]->addItem( "No foreground" );
    m_selects[static_cast<int>(SpecUtils::SpectrumType::Background)]->addItem( "No background" );
    m_selects[static_cast<int>(SpecUtils::SpectrumType::SecondForeground)]->addItem( "No second foreground" );

    // Lastly, select the proper file. If there is no active file in that
    // category, just select the default "No _______" option.
    for( int i = 0; i < 3; ++i )
    {
      const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
      std::shared_ptr<SpecMeas> meas = m_interspec->measurment( type );
      const auto &samples = m_interspec->displayedSamples(type);
      const auto detectors = m_interspec->detectorsToDisplay(type);
      handleDisplayChange( type, meas, samples, detectors );
    }//for( int i = 0; i < 3; ++i )
  }//if( m_files->rowCount() < 1 ) / else
}//void refreshContents()


void CompactFileManager::updateDisplayedScaleFactorNumbers( const double sf,
                                               const SpecUtils::SpectrumType type )
{
  const int typeindex = static_cast<int>(type);
  
  if( !m_scaleValueTxt[typeindex] )
    return;
  
  m_scaleValueTxt[typeindex]->setValue( static_cast<float>(sf) );
}//void updateDisplayedScaleFactorNumbers(...)


//void CompactFileManager::handleSliderChanged( const int slidervalue,
//                                              const SpecUtils::SpectrumType type )
//{
//  bool update = true;
//  const double multiple = (slidervalue-50.0) / 50.0;
//  
//  const double oldsf = m_interspec->displayScaleFactor(type);
//  
//  double sf = oldsf + (0.25 * multiple * oldsf);
//  
//  if( !m_interspec->measurment(type) )
//  {
//    sf = 1.0;
//    update = false;
//  }
//
//  updateDisplayedScaleFactorNumbers( sf, type );
//  
//  if( update )
//    m_interspec->setDisplayScaleFactor( sf, type );
//}//void handleSliderChanged(...);


void CompactFileManager::handleUserEnterdScaleFactor( const SpecUtils::SpectrumType type )
{
  const int typeindex = static_cast<int>(type);
  
  if( !m_scaleValueTxt[typeindex] )
    return;
  
  float sf;
  bool update = true;
  
  const string val = m_scaleValueTxt[typeindex]->text().toUTF8();
  
  if( val.empty() )
  {
    const float lt = m_interspec->liveTime(type);
    const float datalt = m_interspec->liveTime(SpecUtils::SpectrumType::Foreground);
    sf = ((lt>0.0f) ? (datalt/lt) : 1.0f);
//    char buffer[16];
//    snprintf(buffer, sizeof(buffer), "%.2f", sf );
  }else if( !(stringstream(val) >> sf) )
  {
    update = false;
    sf = m_interspec->displayScaleFactor(type);
  }//
  
  updateDisplayedScaleFactorNumbers( sf, type );
  
  
  if( update )
  {
    if( m_rescaleByLiveTime[typeindex] )
      m_rescaleByLiveTime[typeindex]->show();
    m_interspec->setDisplayScaleFactor( sf, type );
  }
}//void handleUserEnterdScaleFactor( const SpecUtils::SpectrumType type )



void CompactFileManager::handleUserEnterdScaleFactorWheel( const SpecUtils::SpectrumType type,  WMouseEvent e )
{
  int i = e.wheelDelta();
  const int typeindex = static_cast<int>(type);
  
  const float oldValue = m_scaleValueTxt[typeindex]->value();
  m_scaleValueTxt[typeindex]->setValue( oldValue + 0.05*oldValue*i);
  handleUserEnterdScaleFactor(type);
}//void handleUserEnterdScaleFactor( const SpecUtils::SpectrumType type )


void CompactFileManager::handleRenormalizeByLIveTime( const SpecUtils::SpectrumType type )
{
  const int typeindex = static_cast<int>(type);
  
  const float lt = m_interspec->liveTime(type);
  const float datalt = m_interspec->liveTime(SpecUtils::SpectrumType::Foreground);
  const double sf = ((lt>0.0f && datalt>0.0f) ? (datalt/lt) : 1.0f);
  updateDisplayedScaleFactorNumbers( sf, type );
  if( m_rescaleByLiveTime[typeindex] )
    m_rescaleByLiveTime[typeindex]->hide();
  m_interspec->setDisplayScaleFactor( sf, type );
}//void handleRenormalizeByLIveTime( const SpecUtils::SpectrumType type )

