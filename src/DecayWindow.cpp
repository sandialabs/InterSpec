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

#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DecayWindow.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayActivityDiv.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace Wt;
using namespace std;


DecayWindow::DecayWindow( InterSpec *viewer )
: AuxWindow( "Nuclide Decay Information",
             //(Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal) | AuxWindowProperties::DisableCollapse) ),
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::DisableCollapse) | AuxWindowProperties::EnableResize) ),
  m_activityDiv( 0 )
{
  m_activityDiv = new DecayActivityDiv( viewer );
  
  WGridLayout *layout = stretcher();
  layout->addWidget( m_activityDiv, 0, 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setRowStretch( 0, 1 );

  const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  const double actUnits = useBq ? PhysicalUnits::MBq : PhysicalUnits::microCi;
  m_activityDiv->addNuclide( 53, 135, 0, 1.0*actUnits, !useBq, 0.0 );

  AuxWindow::addHelpInFooter( footer(), "decay-dialog" );
  Wt::WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );

  rejectWhenEscapePressed();

  if( !viewer || !viewer->isPhone() )
    m_activityDiv->setMinimumSize( 450, 500 );
  
  show();
  setMinimumSize( 460, 580 );
  resizeWindow( 800, 600 );
  resizeToFitOnScreen();
  centerWindow();
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
                        const double activity, const bool useCurries,
                        const double age, const double maxtime )
{
  if( m_activityDiv )
  {
    m_activityDiv->addNuclide( z, a, iso, activity, useCurries, age );
    if( maxtime > 0.0 )
      m_activityDiv->setDecayChartTimeRange( maxtime );
  }//if( m_activityDiv )
}//void addNuclide(...)


void DecayWindow::handleAppUrl( const std::string &path, const std::string &query_str )
{
  m_activityDiv->handleAppUrl( path, query_str );
}//void handleAppUrl( std::string path, std::string query_str )
