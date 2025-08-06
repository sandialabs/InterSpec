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
#include <vector>
#include <sstream>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WContainerWidget>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/SimpleActivityCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

#if( USE_QR_CODES )
#include <Wt/Utils>
#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;

SimpleActivityCalcState::SimpleActivityCalcState()
  : peakEnergy( -1 ),
    nuclideAgeStr( "" ),
    distanceStr( "" ),
    geometryType( 0 ),
    shielding{}
{
}

bool SimpleActivityCalcState::operator==( const SimpleActivityCalcState &rhs ) const
{
  return (peakEnergy == rhs.peakEnergy)
         && (nuclideAgeStr == rhs.nuclideAgeStr)
         && (distanceStr == rhs.distanceStr)
         && (geometryType == rhs.geometryType)
         // TODO: need to implement equality check for `ShieldingSourceFitCalc::ShieldingInfo` && (shielding == rhs.shielding)
  ;
}

bool SimpleActivityCalcState::operator!=( const SimpleActivityCalcState &rhs ) const
{
  return !(*this == rhs);
}

std::string SimpleActivityCalcState::encodeToUrl() const
{
  return "dummy-url-implementation";
}

void SimpleActivityCalcState::decodeFromUrl( const std::string &uri )
{
}

void SimpleActivityCalcState::serialize( std::ostream &out ) const
{
}

void SimpleActivityCalcState::deserialize( std::istream &in )
{
}

SimpleActivityCalcWindow::SimpleActivityCalcWindow( MaterialDB *materialDB,
                                                  Wt::WSuggestionPopup *materialSuggestion,
                                                  InterSpec* viewer )
: AuxWindow( WString::tr("simple-activity-calc-title") ),
  m_tool( nullptr )
{
  setModal( false );
  setResizable( true );
  setClosable( true );

  m_tool = new SimpleActivityCalc( materialDB, materialSuggestion, viewer, contents() );

  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );

#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( footer() );
  qr_btn->setText( "QR Code" );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
#endif

  AuxWindow::addHelpInFooter( footer(), "simple-activity-calc" );

  centerWindow();
  resizeToFitOnScreen();
  show();
}

SimpleActivityCalcWindow::~SimpleActivityCalcWindow()
{
}

SimpleActivityCalc *SimpleActivityCalcWindow::tool()
{
  return m_tool;
}

SimpleActivityCalc::SimpleActivityCalc( MaterialDB *materialDB,
                                      Wt::WSuggestionPopup *materialSuggestion,
                                      InterSpec *specViewer,
                                      Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_renderFlags( 0x0 ),
  m_viewer( specViewer ),
  m_materialSuggest( materialSuggestion ),
  m_materialDB( materialDB ),
  m_peakSelect( nullptr ),
  m_nuclideInfo( nullptr ),
  m_ageEdit( nullptr ),
  m_distanceEdit( nullptr ),
  m_detectorDisplay( nullptr ),
  m_shieldingSelect( nullptr ),
  m_geometrySelect( nullptr ),
  m_resultText( nullptr ),
  m_errorText( nullptr ),
  m_advancedBtn( nullptr ),
  m_distanceRow( nullptr ),
  m_ageRow( nullptr ),
  m_geometryRow( nullptr ),
  m_stateUri( "" ),
  m_currentPeak( nullptr ),
  m_currentNuclide( nullptr ),
  m_currentAge( -1.0 )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
    app->useMessageResourceBundle( "SimpleActivityCalc" );

  addStyleClass( "SimpleActivityCalc" );
  
  init();
}

SimpleActivityCalc::~SimpleActivityCalc()
{
}

void SimpleActivityCalc::init()
{
  addWidget( new WText("Simple Activity Calculator - placeholder implementation") );
}

void SimpleActivityCalc::render( WFlags<RenderFlag> flags )
{
  WContainerWidget::render( flags );
}

void SimpleActivityCalc::setPeakFromEnergy( const double energy )
{
}

void SimpleActivityCalc::handleAppUrl( std::string uri )
{
}

std::string SimpleActivityCalc::encodeStateToUrl() const
{
  return "dummy-state-url";
}

SimpleActivityCalcState SimpleActivityCalc::currentState() const
{
  return SimpleActivityCalcState();
}

void SimpleActivityCalc::setState( const SimpleActivityCalcState &state )
{
}

void SimpleActivityCalc::handlePeakChanged()
{
}

void SimpleActivityCalc::handleDistanceChanged()
{
}

void SimpleActivityCalc::handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf )
{
}

void SimpleActivityCalc::handleGeometryChanged()
{
}

void SimpleActivityCalc::handleShieldingChanged()
{
}

void SimpleActivityCalc::handleSpectrumChanged()
{
}

void SimpleActivityCalc::updateResult()
{
}

void SimpleActivityCalc::handleOpenAdvancedTool()
{
}
