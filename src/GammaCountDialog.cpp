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
#include <iostream>
#include <stdexcept>


#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WSelectionBox>
#include <Wt/WContainerWidget>


#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/GammaCountDialog.h"
#include "InterSpec/NativeFloatSpinBox.h"

using namespace std;
using namespace Wt;


GammaCountDialog::GammaCountDialog( InterSpec *specViewer )
: AuxWindow( "Energy Range Sum",
             (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen)
              | AuxWindowProperties::SetCloseable
              | AuxWindowProperties::DisableCollapse) ),
    m_specViewer( specViewer ),
    m_highlightRegionId( 0 ),
    m_lowerEnergy( NULL ),
    m_upperEnergy( NULL ),
    m_primaryGammaCount( NULL ),
    m_secondaryGammaCount( NULL ),
    m_backgroundGammaCount( NULL ),
    m_liveTimeScaleNote( NULL ),
    m_sigmaAboveBackground( NULL ),
    m_nsigmaHelp( NULL )
{
  init();
  
  handleEnergyRangeChange();
  
  m_specViewer->displayedSpectrumChanged().connect( this, &GammaCountDialog::handleEnergyRangeChange );
  m_specViewer->spectrumScaleFactorChanged().connect( this, &GammaCountDialog::handleEnergyRangeChange );
  
  rejectWhenEscapePressed();
  
  show();
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
  
  if( app && app->isPhone() )
  {
    if( app->viewer() )
    {
      titleBar()->hide();
      
      float safeAreas[4] = { 0.0f };
#if( IOS )
      InterSpecApp::DeviceOrientation orientation = InterSpecApp::DeviceOrientation::Unknown;
      app->getSafeAreaInsets( orientation, safeAreas[0], safeAreas[1], safeAreas[2], safeAreas[3] );
#endif
      repositionWindow( -32768, static_cast<int>(std::max(3.0f,0.5f*safeAreas[0])) );
      setMaximumSize( WLength::Auto, app->viewer()->renderedHeight() - std::max(0.5f*(safeAreas[0]+safeAreas[2]),6.0f) );
    }
  }else
  {
    resizeToFitOnScreen();
    centerWindowHeavyHanded();
  }//if( isPhone ) / else
}//GammaCountDialog constructor


GammaCountDialog::~GammaCountDialog()
{
  if( m_specViewer )
  {
    m_specViewer->removeHighlightedEnergyRange( m_highlightRegionId );
    m_highlightRegionId = 0;
  }
}//~GammaCountDialog()


void GammaCountDialog::init()
{
  wApp->useStyleSheet( "InterSpec_resources/GammaCountDialog.css" );
  wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_specViewer );
  
  if( !m_specViewer )
    throw runtime_error( "GammaCountDialog: you must pass in valid InterSpec pointer" );

  auto hist = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  const double xlow = hist ? std::min( hist->gamma_energy_min() - 10.0, 0.0 ) : 0.0;
  const double xhigh = hist ? 2.0 * hist->gamma_energy_max() : 15000.0;
  
  WContainerWidget *body = contents();
  body->addStyleClass( "GammaCountDialog" );

  WText *text = new WText( "Count the number of gammas in the specified energy range.", body );
  text->setInline( true );
  text->setStyleClass( "line-below" );
  
  WContainerWidget *inputs = new WContainerWidget( body );
  inputs->addStyleClass( "inputs" );

  WLabel *label = new WLabel( "Lower Energy", inputs );
  label->addStyleClass( "GridFirstRow GridFirstCol EnergyLabel" );
  m_lowerEnergy = new NativeFloatSpinBox( inputs );
  m_lowerEnergy->addStyleClass( "GridFirstRow GridSecondCol" );
  m_lowerEnergy->setSpinnerHidden();
  label->setBuddy( m_lowerEnergy );
  label = new WLabel( "keV", inputs );
  label->addStyleClass( "GridFirstRow GridThirdCol GridJustifyStart" );
  m_lowerEnergy->setRange( xlow, xhigh );
  
  label = new WLabel( "Upper Energy", inputs );
  label->addStyleClass( "GridSecondRow GridFirstCol EnergyLabel", inputs );
  m_upperEnergy = new NativeFloatSpinBox( inputs );
  label->setBuddy( m_upperEnergy );
  m_upperEnergy->addStyleClass( "GridSecondRow GridSecondCol" );
  m_upperEnergy->setSpinnerHidden();
  label = new WLabel( "keV", inputs );
  label->addStyleClass( "GridSecondRow GridThirdCol GridJustifyStart" );
  m_upperEnergy->setRange( xlow, xhigh );
  
  
  if( m_specViewer->isMobile() )
  {
    inputs->addStyleClass( "line-below" );
  }else
  {
    string key_sequence = "Shift-Alt-Drag";  // Linux / Windows
    const string &user_agent = wApp->environment().userAgent();
    if( SpecUtils::icontains( user_agent, "macos")
       || SpecUtils::icontains( user_agent, "macintel")
       || SpecUtils::icontains( user_agent, "macintosh")
       || SpecUtils::icontains( user_agent, "iphone")
       || SpecUtils::icontains( user_agent, "ipad") )
    {
      key_sequence = "Shift-Option-Drag"; //Apple
    }
    
    text = new WText("You can also <b>" + key_sequence + "</b> on the chart"
                              " to select the energy range", body);
    text->setInline( false );
    text->setStyleClass( "line-below shortcut-info" );
  }//if( !mobile )
  
    
//  if( m_specViewer->isMobile() )
//  {
//20150123: on android at least, calling setNativeControl() causes
//  a javascript exception (having to do with the validate)
//    m_lowerEnergy->setNativeControl(true); //mobile should not show spinner
//    m_upperEnergy->setNativeControl(true); //mobile should not show spinner
//  } //(isMobile())
  
  WContainerWidget *answers = new WContainerWidget( body );
  answers->addStyleClass( "answers" );
  
  label = new WLabel( "Foreground Counts:", answers );
  label->addStyleClass( "GridFirstRow GridFirstCol" );
  m_primaryGammaCount = new WText( answers );
  m_primaryGammaCount->addStyleClass( "GridFirstRow GridSecondCol" );
  
  label = new WLabel( "Secondary Counts:", answers );
  label->addStyleClass( "GridSecondRow GridFirstCol" );
  m_secondaryGammaCount = new WText( answers );
  m_secondaryGammaCount->addStyleClass( "GridSecondRow GridSecondCol" );
  
  label = new WLabel( "Background Counts:", answers );
  label->addStyleClass( "GridThirdRow GridFirstCol" );
  m_backgroundGammaCount = new WText( answers );
  m_backgroundGammaCount->addStyleClass( "GridThirdRow GridSecondCol" );
  
  
  m_sigmaAboveBackground = new WText( "", XHTMLText, answers );
  m_sigmaAboveBackground->addStyleClass( "line-above GridFourthRow GridFirstCol GridJustifyCenter GridSpanTwoCol nsigma" );
  
  m_nsigmaHelp = new WImage( answers );
  m_nsigmaHelp->setImageLink(Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
  m_nsigmaHelp->setStyleClass("Wt-icon GridFourthRow GridThirdCol GridJustifyEnd");
  m_nsigmaHelp->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
  m_nsigmaHelp->setHidden( true );
  
  const char *tooltip = "The number of sigma difference is calculated by dividing the"
                        " difference between the foreground and the (scaled) background,"
                        " divided by the uncertainty. The uncertainty is calculated by"
                        " summing the statistical uncertainties of the foreground and background in"
                        " quadrature, then taking the square root";
  HelpSystem::attachToolTipOn( m_nsigmaHelp, tooltip, true, HelpSystem::ToolTipPosition::Right );
  
  m_liveTimeScaleNote = new WText( "", XHTMLText, answers );
  m_liveTimeScaleNote->addStyleClass( "GridFifthRow GridFirstCol GridJustifyCenter GridSpanThreeCol line-above LiveTimeScaleNote" );  //""
  m_liveTimeScaleNote->hide();
  
  
  m_lowerEnergy->valueChanged().connect( this, &GammaCountDialog::handleEnergyRangeChange );
  m_upperEnergy->valueChanged().connect( this, &GammaCountDialog::handleEnergyRangeChange );
  m_specViewer->displayedSpectrumChanged().connect(
                boost::bind( &GammaCountDialog::handleSpectrumChange, this,
                            boost::placeholders::_1, boost::placeholders::_2,
                            boost::placeholders::_3, boost::placeholders::_4) );
  
 
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &GammaCountDialog::emitFinished );
  
  //Keep the keyboard form popping up
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(WApplication::instance());
  if( app && app->isMobile() )
    closeButton->setFocus();
}//GammaCountDialog::init()


void GammaCountDialog::emitFinished()
{
  finished().emit( WDialog::Rejected );
}

void GammaCountDialog::setEnergyRange( double lowEnergy, double highEnergy )
{
  if( highEnergy <= lowEnergy )
  {
    lowEnergy = 0.0;
    highEnergy = 0.0;
  }//if( highEnergy <= lowEnergy )

  auto hist = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( hist && (hist->num_gamma_channels() > 0) )
  {
    const double xmin = hist ? std::min( hist->gamma_energy_min(), 0.0f ) : 0.0f;
    const double xmax = hist ? hist->gamma_energy_max() : 15000.0;
    if( lowEnergy < xmin )
      lowEnergy = xmin;
    if( lowEnergy > xmax )
      lowEnergy = xmax;
    if( highEnergy < xmin )
      highEnergy = xmin;
    if( highEnergy > xmax )
      highEnergy = xmax;
  }//
  
  m_lowerEnergy->setValue( lowEnergy );
  m_upperEnergy->setValue( highEnergy );

  const WString lowerStrVal = m_lowerEnergy->text();
  const WString upperStrVal = m_upperEnergy->text();
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    const WString prevLowerEnergy = m_prevLowerEnergy;
    const WString prevUpperEnergy = m_prevUpperEnergy;
    
    auto setRange = [=]( const bool is_undo ){
      InterSpec *viewer = InterSpec::instance();
      GammaCountDialog *dialog = viewer ? viewer->showGammaCountDialog() : nullptr;
      if( dialog )
      {
        dialog->m_lowerEnergy->setText( is_undo ? prevLowerEnergy : lowerStrVal );
        dialog->m_upperEnergy->setText( is_undo ? prevUpperEnergy : upperStrVal );
        dialog->handleEnergyRangeChange();
      }
    };
    auto undo = [setRange](){ setRange(true); };
    auto redo = [setRange](){ setRange(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Set energy range to sum." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  
  m_prevLowerEnergy = lowerStrVal;
  m_prevUpperEnergy = upperStrVal;
  
  handleEnergyRangeChange();
}//void setEnergyRange( const double lowEnergy, const double highEnergy )


void GammaCountDialog::handleAppUrl( std::string query_str )
{
  const map<string,string> values = AppUtils::query_str_key_values( query_str );
  
  // Check version is appropriate
  const auto ver_iter = values.find( "V" );
  if( (ver_iter == end(values))
     || ((ver_iter->second != "1") && !SpecUtils::istarts_with(ver_iter->second, "1.")) )
    throw runtime_error( "Energy Range Sum: URI not compatible version." );
  
  const auto lower_iter = values.find( "LOW" );
  const auto upper_iter = values.find( "HIGH" );
  string lower_str = (lower_iter != end(values)) ? lower_iter->second : string();
  string upper_str = (upper_iter != end(values)) ? upper_iter->second : string();
  try{ std::stod(lower_str); }catch(...){ lower_str = ""; }
  try{ std::stod(upper_str); }catch(...){ upper_str = ""; }
  
  auto primary = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( lower_str.empty() && primary && (primary->num_gamma_channels() > 0) )
    lower_str = std::to_string( primary->gamma_channel_lower(0) );
  if( upper_str.empty() && primary && (primary->num_gamma_channels() > 0) )
    upper_str = std::to_string( primary->gamma_channel_upper(primary->num_gamma_channels() - 1) );
  
  m_lowerEnergy->setText( WString::fromUTF8(lower_str) );
  m_upperEnergy->setText( WString::fromUTF8(upper_str) );
    
  handleEnergyRangeChange();
}//void handleAppUrl( std::string query_str )


std::string GammaCountDialog::encodeStateToUrl() const
{
  string low, high;
  if( m_lowerEnergy->validate() == WValidator::State::Valid )
    low = m_lowerEnergy->text().toUTF8();
  if( m_upperEnergy->validate() == WValidator::State::Valid )
    high = m_upperEnergy->text().toUTF8();
  
  return "V=1&LOW=" + low + "&HIGH=" + high;
}//std::string encodeStateToUrl() const


void GammaCountDialog::handleEnergyRangeChange()
{
  try
  {
    float low = static_cast<float>( m_lowerEnergy->value() );
    float high = static_cast<float>( m_upperEnergy->value() );
    
    if( low > high )
    {
      std::swap( low, high );
      m_lowerEnergy->setValue( low );
      m_upperEnergy->setValue( high );
    }//if( low > high )
    
  
    std::shared_ptr<const SpecUtils::Measurement> hist = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if( hist && hist->num_gamma_channels() > 0 )
    {
      const size_t nchannel = hist->num_gamma_channels();
      const float xlow = hist->gamma_channel_lower( 0 );
      const float xhigh = hist->gamma_channel_upper( nchannel - 1 );
      if( low < (xlow-10.0f) )
        m_lowerEnergy->setValue( low = xlow );
      if( high > (xhigh+10.0f) )
        m_upperEnergy->setValue( high = xhigh );
    }//if( hist )

    
    m_specViewer->removeHighlightedEnergyRange( m_highlightRegionId );
    WColor color( 255, 255, 0, 155 );
    m_highlightRegionId = m_specViewer->addHighlightedEnergyRange( low, high, color );

    //Now need to do all the computations and such
    std::shared_ptr<const SpecUtils::Measurement> primary, secondary, background;
    primary    = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    secondary  = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::SecondForeground );
    background = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::Background );

    const double secondarySF = m_specViewer->displayScaleFactor( SpecUtils::SpectrumType::SecondForeground );
    const double backgroundSF = m_specViewer->displayScaleFactor( SpecUtils::SpectrumType::Background );

    setGammaCountText( m_primaryGammaCount,    primary, 1.0f, low, high );
    setGammaCountText( m_secondaryGammaCount,  secondary, secondarySF, low, high );
    setGammaCountText( m_backgroundGammaCount, background, backgroundSF, low, high );

    string secondarySfStr = "--", backgroundSfStr = "--";
    string livetimescale = "Note: ";
    char onstack[64];

    if( background )
    {
      if( backgroundSF > 0.1 )
        snprintf( onstack, sizeof(onstack), "%.2f", backgroundSF );
      else
        snprintf( onstack, sizeof(onstack), "%.2g", backgroundSF );
      backgroundSfStr = onstack;
      
      livetimescale += "Background is being scaled by " + backgroundSfStr ;
    }//if( background )

    if( secondary )
    {
      if( secondarySF > 0.1 )
        snprintf( onstack, sizeof(onstack), "%.2f", secondarySF );
      else
        snprintf( onstack, sizeof(onstack), "%.2g", secondarySF );
      
      secondarySfStr = onstack;
      
      if( background )
        livetimescale += " and secondary by ";
      else
        livetimescale += "Secondary background is being scaled by ";
      
      livetimescale += secondarySfStr;
    }//if( secondary )

      
    if( secondary || background )
    {
      livetimescale += ".";
      m_liveTimeScaleNote->setText( livetimescale );
      m_liveTimeScaleNote->show();
    }else
    {
      m_liveTimeScaleNote->setText("");
      m_liveTimeScaleNote->hide();
    }
    
    setUnceraintyText( primary, background, backgroundSF, low, high );
    
    double lowEnergy = 9999999.9, upEnergy = -9999999.9;
    if( primary )
    {
      const double loweste = primary->gamma_channel_lower( 0 );
      const double higheste = primary->gamma_channel_energies()->back();
      lowEnergy = std::min( lowEnergy, loweste );
      upEnergy = std::max( upEnergy, higheste );
    }else
    {
      m_sigmaAboveBackground->setText("");
      m_nsigmaHelp->hide();
    }

    if( secondary )
    {
      const double loweste = secondary->gamma_channel_lower( 0 );
      const double higheste = secondary->gamma_channel_energies()->back();
      lowEnergy = std::min( lowEnergy, loweste );
      upEnergy = std::max( upEnergy, higheste );
    }//if( secondary )

    if( background )
    {
      const double loweste = background->gamma_channel_lower( 0 );
      const double higheste = background->gamma_channel_energies()->back();
      lowEnergy = std::min( lowEnergy, loweste );
      upEnergy = std::max( upEnergy, higheste );
    }//if( background )

    upEnergy *= 2.0;
    m_lowerEnergy->setRange( lowEnergy-10.0, upEnergy );
    m_upperEnergy->setRange( lowEnergy-10.0, upEnergy );
  }catch( std::exception &e )
  {
    cerr << "\nGammaCountDialog::handleEnergyRangeChange()\n\tCaught:" << e.what() << endl << endl;
  }
}//void handleEnergyRangeChange()


void GammaCountDialog::setUnceraintyText( std::shared_ptr<const SpecUtils::Measurement> foreground,
                                          std::shared_ptr<const SpecUtils::Measurement> background,
                                          const double backSF,
                                          const float minEnergy,
                                          const float maxEnergy )
{
  if( !foreground || !background )
  {
    m_sigmaAboveBackground->setText("");
    m_nsigmaHelp->hide();
    return;
  }//if( !foreground || !background )
  
  const double nfore = gamma_integral( foreground, minEnergy, maxEnergy );
  const double nback = gamma_integral( background, minEnergy, maxEnergy );
  const double scaleback = nback * backSF;
  const double backsigma = sqrt(nback);
  const double forsigma = sqrt(nfore);
  const double backscalesigma = backSF * backsigma;
  const double total_back_fore_sigma = sqrt(backscalesigma*backscalesigma + forsigma*forsigma);
  const double nsigma = fabs(nfore - scaleback) / total_back_fore_sigma;
  const bool isneg = (scaleback > nfore);
  
  char onstack[32];
  snprintf( onstack, sizeof(onstack), "%.5g ", nsigma );
  string sigmatxt = onstack;
  
//  stringstream sigmatxtstrm;
//  if( nsigma > 100.0 || nsigma < 0.01 )
//    sigmatxtstrm << std::setprecision(4) << std::scientific << nsigma;
//  else
//    sigmatxtstrm << std::setprecision(3) << std::fixed << nsigma;
  
  WString txt = "Foreground is " + sigmatxt;
//  WString txt = "Foreground is " + sigmatxtstrm.str();
  
#ifndef WT_NO_STD_WSTRING
  txt += L"\x03C3";
#else
  txt += "&sigma;";
#endif
  txt += (isneg ? " below background." : " above background.");
  
  if( nback == 0 )
  {
    //Shouldnt normally happen
    m_sigmaAboveBackground->setText("");
    m_nsigmaHelp->hide();
  }else
  {
    m_sigmaAboveBackground->setText( txt );
    m_nsigmaHelp->show();
  }
}//void GammaCountDialog::setUnceraintyText(...)


void GammaCountDialog::setGammaCountText( Wt::WText *text, std::shared_ptr<const SpecUtils::Measurement> hist,
                                          const double scale_factor,
                                          const float minEnergy,
                                          const float maxEnergy )
{
  if( !hist || maxEnergy <= minEnergy )
  {
    text->setText( "--" );
    return;
  }//if( !hist )

  const double count = scale_factor * gamma_integral( hist, minEnergy, maxEnergy );
  
  char buffer[32];
  if( count > 1.0E5 )
    snprintf( buffer, sizeof(buffer), "%.5g", count );
  else
    snprintf( buffer, sizeof(buffer), "%.1f", count );
  
  text->setText( buffer );
}//void setGammaCountText( Wt::WText *text, std::shared_ptr<const SpecUtils::Measurement> hist, double minEnergy, double maxEnergy )






void GammaCountDialog::handleSpectrumChange( const SpecUtils::SpectrumType /*type*/,
                                             const std::shared_ptr<SpecMeas> &/*meas*/,
                                             const std::set<int> &/*displaySample*/,
                                             const std::vector<std::string> & /*detectors*/)
{
  handleEnergyRangeChange();
}//void handleSpectrumChange(...)



