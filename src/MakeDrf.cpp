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

#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WSuggestionPopup>

#include "InterSpec/MakeDrf.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;
using namespace Wt;


MakeDrfWindow::MakeDrfWindow( InterSpec *viewer, MaterialDB *materialDB, Wt::WSuggestionPopup *materialSuggest )
  : AuxWindow( "Create DRF",
  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletModal)
  | AuxWindowProperties::SetCloseable
  | AuxWindowProperties::DisableCollapse) ),
  m_makeDrf( nullptr )
{
  m_makeDrf = new MakeDrf( viewer, materialDB, materialSuggest, contents() );
  
  AuxWindow::addHelpInFooter( footer(), "make-drf" );
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  
  //const int screenW = viewer->renderedWidth();
  //const int screenH = viewer->renderedHeight();
  //const int width = ((screenW < 600) ? screenW : 600);
  //const int height = ((screenH < 420) ? screenH : 420);
  //resizeWindow( width, height );
  
  resizeToFitOnScreen();
  centerWindow();
  rejectWhenEscapePressed( false );
}//GammaXsWindow(...) constrctor


MakeDrfWindow::~MakeDrfWindow()
{
}


MakeDrf::MakeDrf( InterSpec *viewer, MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggest,
                  Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_materialDB( materialDB ),
  m_materialSuggest( materialSuggest )
{
  assert( m_interspec );
  assert( m_materialDB );
  assert( m_materialSuggest );
  
  wApp->useStyleSheet( "InterSpec_resources/MakeDrf.css" );
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nuc = db->nuclide( "Cs137" );
  
  boost::posix_time::ptime measDate = UtilityFunctions::time_from_string( "2019-01-01T00:00:00" );
  
  
  new MakeDrfSrcDef( nuc, measDate, materialDB, materialSuggest, this );
  
}//MakeDrf( constructor )


MakeDrf::~MakeDrf()
{
  
}//~MakeDrf()
