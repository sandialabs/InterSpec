//  Decay.cpp, Created by Chan, Ethan on 12/18/13.

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

#include <math.h>
#include <iostream>

#include <Wt/Utils>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DecayWindow.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayActivityDiv.h"
#include "InterSpec/DecayDataBaseServer.h"

#if( USE_QR_CODES )
#include "InterSpec/QrCode.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#endif

using namespace Wt;
using namespace std;


DecayWindow::DecayWindow( InterSpec *viewer )
: AuxWindow( WString::tr("window-title-nuc-decay"),
             //(Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal) | AuxWindowProperties::DisableCollapse) ),
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::DisableCollapse) | AuxWindowProperties::EnableResize) ),
  m_activityDiv( 0 )
{
  if( viewer )
    viewer->useMessageResourceBundle( "DecayActivity" );
  
  m_activityDiv = new DecayActivityDiv( viewer );
  
  WGridLayout *layout = stretcher();
  layout->addWidget( m_activityDiv, 0, 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setRowStretch( 0, 1 );

  const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  const double actUnits = useBq ? PhysicalUnits::MBq : PhysicalUnits::microCi;
  const string actStr = useBq ? "1 MBq" : "1 uCi";
  m_activityDiv->addNuclide( 53, 135, 0, 1.0*actUnits, !useBq, 0.0, actStr);

  AuxWindow::addHelpInFooter( footer(), "decay-dialog" );
  Wt::WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );

  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( footer() );
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://decay/"
                                     + Wt::Utils::urlEncode(m_activityDiv->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("dw-qr-window-title"), WString::tr("dw-qr-window-text") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES

  rejectWhenEscapePressed();

  show();
  
  if( viewer && !viewer->isPhone() && (viewer->renderedWidth() > 100) && (viewer->renderedHeight() > 100) )
  {
    const int w = std::min( 800, viewer->renderedWidth() - 20 );
    const int h = std::min( 600, viewer->renderedHeight() );
    
    m_activityDiv->setMinimumSize( 450, std::min(450, h - 20) );
    resizeWindow( w, h );
    
    resizeToFitOnScreen();
    centerWindow();
  }else
  {
    // We get here if we are like opening a URL (e.g., "deep link") on application start
    //  On an iPhone, the AuxWindow is larger than the screen, and some of the contents are not
    //  visible.  The below tries to fix it, but it doesnt seem to work.
    //
    // As a note: iPhone 13 has about 812 x 375
    // TODO: check if we could use wApp->environment().screenWidth() // screenHeight()
    if( !viewer )
    {
      // good luck - we really shouldnt be here
    }else if( (viewer->renderedWidth() <= 100) || (viewer->renderedHeight() <= 100) )
    {
      if( !viewer->isPhone() )
        m_activityDiv->setMinimumSize( 450, 500 );
    }//
    
    // On phones, the below doesnt seem to work
    const bool wasPhone = m_isPhone;
    m_isPhone = false;
    if( wasPhone )
      resizeScaledWindow( 1.0, 1.0 ); //wouldnt have an effect if m_isPhone is true
    centerWindowHeavyHanded();        //wouldnt have an effect if m_isPhone is true
    m_isPhone = wasPhone;
  }//if( we know the screens width ) / else
}//Decay constructor


DecayWindow::~DecayWindow()
{
  //nothing to do here
}//Decay destructor


void DecayWindow::clearAllNuclides()
{
  if( m_activityDiv )
    m_activityDiv->clearAllNuclides();
}//void Decay::clearAllNuclides()


void DecayWindow::addNuclide( const int z, const int a, const int iso,
                        const double activity, const bool useCuries,
                        const double age, std::string activityStr, const double maxtime )
{
  if( m_activityDiv )
  {
    if (activityStr.empty())
    {
      activityStr = PhysicalUnits::printToBestActivityUnits(activity, 3, useCuries);
    }else
    {
      assert((activity < 1.0E-6)
        || (fabs(activity - PhysicalUnits::stringToActivity(activityStr)) < 0.001 * activity));
    }

    m_activityDiv->addNuclide( z, a, iso, activity, useCuries, age, activityStr);
    if( maxtime > 0.0 )
      m_activityDiv->setDecayChartTimeRange( maxtime );
  }//if( m_activityDiv )
}//void addNuclide(...)


void DecayWindow::handleAppUrl( std::string url )
{
  std::string host, path, query, frag;
  AppUtils::split_uri( url, host, path, query, frag );
  
  handleAppUrl( path.empty() ? host : path, query );
}//handleAppUrl( std::string url )


void DecayWindow::handleAppUrl( const std::string &path, const std::string &query_str )
{
  m_activityDiv->handleAppUrl( path, query_str );
}//void handleAppUrl( std::string path, std::string query_str )


std::string DecayWindow::encodeStateToUrl()
{
  return m_activityDiv->encodeStateToUrl();
}
