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
#include <Wt/WApplication>
#include <Wt/WSelectionBox>
#include <Wt/WDoubleSpinBox>
#include <Wt/WPushButton>
#include <Wt/WContainerWidget>


#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/GammaCountDialog.h"
#include "InterSpec/SpectrumDisplayDiv.h"

using namespace std;
using namespace Wt;


GammaCountDialog::GammaCountDialog( InterSpec *specViewer )
: AuxWindow( "Energy Range Sum",
             (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneModal)
              | AuxWindowProperties::SetCloseable
              | AuxWindowProperties::DisableCollapse) ),
    m_specViewer( specViewer ),
    m_highlightRegionId( 0 ),
    m_lowerEnergy( NULL ),
    m_upperEnergy( NULL ),
    m_primaryGammaCount( NULL ),
    m_secondaryGammaCount( NULL ),
    m_backgroundGammaCount( NULL ),
    m_secondaryLiveTimeScale( NULL ),
    m_backgroundLiveTimeScale( NULL ),
    m_sigmaAboveBackground( NULL ),
    m_nsigmaHelp( NULL )
{
  init();
  
  handleEnergyRangeChange();
  
  m_specViewer->displayedSpectrumChanged().connect( this, &GammaCountDialog::handleEnergyRangeChange );
  
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
    centerWindow();
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
  
  if( !m_specViewer )
    throw runtime_error( "GammaCountDialog: you must pass in valid InterSpec pointer" );
 
  const bool showToolTipInstantly
         = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_specViewer );
  
  WGridLayout *layout = new WGridLayout();
  contents()->setLayout( layout );
  contents()->setOverflow( WContainerWidget::OverflowHidden );

  m_lowerEnergy = new WDoubleSpinBox();
  m_upperEnergy = new WDoubleSpinBox();
//  if( m_specViewer->isMobile() )
//  {
//20150123: on android at least, calling setNativeControl() causes
//  a javascript exception (having to do with the validate)
//    m_lowerEnergy->setNativeControl(true); //mobile should not show spinner
//    m_upperEnergy->setNativeControl(true); //mobile should not show spinner
//  } //(isMobile())
  
  m_primaryGammaCount       = new WText();
  m_secondaryGammaCount     = new WText();
  m_backgroundGammaCount    = new WText();
  m_secondaryLiveTimeScale  = new WText( "", XHTMLUnsafeText );
  m_backgroundLiveTimeScale = new WText( "", XHTMLUnsafeText );
  m_sigmaAboveBackground    = new WText( "", XHTMLUnsafeText );
  m_nsigmaHelp              = new WImage();
  
//  m_sigmaAboveBackground->setHiddenKeepsGeometry(true);
  m_lowerEnergy->setRange( 0.0, 10000.0 );
  m_upperEnergy->setRange( 0.0, 10000.0 );

  std::shared_ptr<const SpecUtils::Measurement> hist;
  if( m_specViewer )
    hist = m_specViewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( hist )
  {
    const double xlow = std::min( hist->gamma_energy_min() -10.0, 0.0 );
    const double xhigh = hist->gamma_energy_max();
    m_lowerEnergy->setRange( xlow, 2.0*xhigh );
    m_upperEnergy->setRange( xlow, 2.0*xhigh );
  }//if( hist )
  
  m_nsigmaHelp->setImageLink(Wt::WLink("InterSpec_resources/images/help.png") );
  m_nsigmaHelp->setStyleClass("helpIconTitleBar");
  m_nsigmaHelp->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
  m_nsigmaHelp->setHidden( true );
  
  const char *tooltip = "The number of sigma difference is calculated by dividing the"
                        " difference between the foreground and the (scaled) background,"
                        " divided by the uncertainty. The uncertainty is calculated by"
                        " summing the statistical uncertainties of the foreground and background in"
                        " quadrature, then taking the squareroot";
  HelpSystem::attachToolTipOn( m_nsigmaHelp, tooltip, true, HelpSystem::Right );

  WLabel *label = NULL;
  WText *text = NULL;
  text = new WText( "Count the number of gammas in the specified energy range." );
  text->setWordWrap(true);
  text->setWidth( WLength(10.5, WLength::FontEm) );
  int row=0;
  text->setStyleClass("line-below");
  layout->addWidget( text, row, 0, 1, 3/*, AlignCenter */);

  label = new WLabel( "Lower Energy" );
  layout->addWidget( label, ++row, 0, 1, 1, AlignLeft );
  layout->addWidget( m_lowerEnergy, row, 1, 1, 1, AlignLeft );
  label = new WLabel( "keV" );
  layout->addWidget( label, row, 2, 1, 1, AlignLeft );

  label = new WLabel( "Upper Energy" );
  layout->addWidget( label, ++row, 0, 1, 1, AlignLeft );
  layout->addWidget( m_upperEnergy, row, 1, 1, 1, AlignLeft );
  label = new WLabel( "keV" );
  layout->addWidget( label, row, 2, 1, 1, AlignLeft );

  if( m_specViewer && !m_specViewer->isMobile() )
  {
    WLabel* temp = new WLabel("<center><i><small>You can also <b>Shift-Alt-Drag</b> on the chart"
                              " to select the energy range</small></i></center>");
    temp->setWordWrap(true);
    layout->addWidget(temp,++row, 0 , 1 ,3,AlignCenter);
    temp->setStyleClass("line-below");
  }
    
  label = new WLabel( "Foreground Counts:" );
  layout->addWidget( label, ++row, 0, 1, 1, AlignLeft );
  layout->addWidget( m_primaryGammaCount, row, 1, 1, 2, AlignLeft );

  label = new WLabel( "Secondary Counts:" );
  layout->addWidget( label, ++row, 0, 1, 1, AlignLeft );
  layout->addWidget( m_secondaryGammaCount, row, 1, 1, 2, AlignLeft );

  label = new WLabel( "Background Counts:" );
  layout->addWidget( label, ++row, 0, 1, 1, AlignLeft );
  layout->addWidget( m_backgroundGammaCount, row, 1, 1, 2, AlignLeft );
    
  layout->addWidget( m_sigmaAboveBackground, ++row, 0, 1, 2, AlignCenter );
  layout->addWidget( m_nsigmaHelp, row, 2, 1, 1, AlignRight );
  
  layout->addWidget( m_backgroundLiveTimeScale, ++row, 0, 1, 3, AlignCenter );
  m_backgroundLiveTimeScale->setStyleClass("line-above");
    
  if( m_specViewer && !m_specViewer->isMobile() )
  {
    const char *txt = "Press <b>Enter</b> after typing in the energy";
//    text = new WText( txt, XHTMLText );
//    layout->addWidget( text, ++row, 0, 1, 3, AlignCenter );
      HelpSystem::attachToolTipOn(m_lowerEnergy, txt, showToolTipInstantly);
      HelpSystem::attachToolTipOn(m_upperEnergy, txt, showToolTipInstantly);
//    txt = "If you manually type an energy you will need to either press enter,"
//          " or click another field before changes take effect.";
//    HelpSystem::attachToolTipOn( contents(), txt, showToolTipInstantly );
  }//if( !ismobile )
  
  layout->setColumnStretch( 1, 1 );

  m_lowerEnergy->valueChanged().connect( this, &GammaCountDialog::handleEnergyRangeChange );
  m_upperEnergy->valueChanged().connect( this, &GammaCountDialog::handleEnergyRangeChange );
  m_specViewer->displayedSpectrumChanged().connect( boost::bind( &GammaCountDialog::handleSpectrumChange, this, _1, _2, _3, _4) );
  
 
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

void GammaCountDialog::setEnergyRange( double lowEnergy,
                                       double highEnergy )
{

  if( highEnergy <= lowEnergy )
  {
    lowEnergy = 0.0;
    highEnergy = 0.0;
  }//if( highEnergy <= lowEnergy )

  if( lowEnergy < m_lowerEnergy->minimum() )
    m_lowerEnergy->setMinimum( lowEnergy - 1.0 );
  if( highEnergy > m_upperEnergy->maximum()  )
    m_upperEnergy->setMaximum( highEnergy + 1.0 );

  m_lowerEnergy->setValue( lowEnergy );
  m_upperEnergy->setValue( highEnergy );

  handleEnergyRangeChange();
}//void setEnergyRange( const double lowEnergy, const double highEnergy )


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
    string livetimescale = "<i>Note: ";
    char onstack[64];
//      m_secondaryLiveTimeScale->hide();
//      m_backgroundLiveTimeScale->hide();

    if( background )
    {
        if( backgroundSF > 0.1 )
          snprintf( onstack, sizeof(onstack), "%.2f", backgroundSF );
        else
          snprintf( onstack, sizeof(onstack), "%.2g", backgroundSF );
        backgroundSfStr = onstack;

        livetimescale += "Background is being scaled by " + backgroundSfStr ;
    }//else backgroundSfStrm << "--";

    if( secondary )
    {
        if( secondarySF > 0.1 )
            snprintf( onstack, sizeof(onstack), "%.2f", secondarySF );
        else
            snprintf( onstack, sizeof(onstack), "%.2g", secondarySF );
          
        secondarySfStr = onstack;

        if (background)
            livetimescale += " and secondary by ";
        else
            livetimescale += "Secondary background is being scaled by ";
        
        livetimescale += secondarySfStr;
    }//else secondarySfStrm << "--";

      
      if (secondary||background) {
          livetimescale += ".</i>";
          m_backgroundLiveTimeScale->setText(livetimescale);
//          m_backgroundLiveTimeScale->show();
      } else {
          m_backgroundLiveTimeScale->setText("");
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



