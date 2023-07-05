


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

#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WSpinBox>
#include <Wt/WCheckBox>
#include <Wt/WApplication>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/FeatureMarkerWidget.h"

using namespace Wt;
using namespace std;


FeatureMarkerWindow::FeatureMarkerWindow( InterSpec *viewer )
  : AuxWindow( "Feature Markers",
    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen)
     | AuxWindowProperties::SetCloseable
     | AuxWindowProperties::DisableCollapse) ),
    m_feature( nullptr )
{
  //ToDo: for mobile, replace the close icon with a larger button.
  
  addStyleClass( "FeatureMarkerWindow" );
  
  rejectWhenEscapePressed( false );
  
  m_feature = new FeatureMarkerWidget( viewer, contents() );
  //m_feature->setHeight( WLength(100,WLength::Percentage) );
  
  // This next call seems to help resize the window to show all the contents, otherwise "Sum Peak"
  //  will hang-off the bottom of the window.  Definitely a hack.
  wApp->doJavaScript( wApp->javaScriptClass() + ".TriggerResizeEvent();" );
  
  if( viewer->isMobile() )
    setModal( false );
  
  show();
  
  //const int screenW = viewer->renderedWidth();
  //const int screenH = viewer->renderedHeight();
  //const int width = ((screenW < 600) ? screenW : 600);
  //const int height = ((screenH < 420) ? screenH : 420);
  //resizeWindow( width, height );
  
  repositionWindow( 50, 25 );
}//GammaXsWindow(...) constrctor


FeatureMarkerWindow::~FeatureMarkerWindow()
{
}


void FeatureMarkerWindow::setFeatureMarkerChecked( const FeatureMarkerType option, const bool checked )
{
  if( m_feature )
    m_feature->setFeatureMarkerChecked( option, checked );
}


void FeatureMarkerWindow::setDisplayedComptonPeakAngle( const int angle )
{
  if( m_feature )
    m_feature->setDisplayedComptonPeakAngle( angle );
}


FeatureMarkerWidget::FeatureMarkerWidget( InterSpec* viewer, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_viewer( viewer ),
    m_escapePeaks( nullptr ),
    m_comptonPeak( nullptr ),
    m_comptonAngle( nullptr ),
    m_comptonEdge( nullptr ),
    m_sumPeaks( nullptr )
{
  init();
}

FeatureMarkerWidget::~FeatureMarkerWidget()
{
}


void FeatureMarkerWidget::setFeatureMarkerChecked( const FeatureMarkerType option, const bool checked )
{
  switch( option )
  {
    case FeatureMarkerType::EscapePeakMarker:
      m_escapePeaks->setChecked( checked );
      break;
      
    case FeatureMarkerType::ComptonPeakMarker:
      m_comptonPeak->setChecked( checked );
      break;
      
    case FeatureMarkerType::ComptonEdgeMarker:
      m_comptonEdge->setChecked( checked );
      break;
      
    case FeatureMarkerType::SumPeakMarker:
      m_sumPeaks->setChecked( checked );
      break;
      
    case FeatureMarkerType::NumFeatureMarkers:
      break;
  }//switch( option )
}//setFeatureMarkerChecked(...)


void FeatureMarkerWidget::setDisplayedComptonPeakAngle( const int angle )
{
  m_comptonAngle->setValue( angle );
}

void FeatureMarkerWidget::init()
{
  wApp->useStyleSheet( "InterSpec_resources/FeatureMarkerWidget.css" );
  addStyleClass( "FeatureMarkerWidget" );
  
  m_escapePeaks = new WCheckBox( "Escape Peaks", this );
  m_escapePeaks->setInline( false );
  m_escapePeaks->checked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                       FeatureMarkerType::EscapePeakMarker, true ) );
  m_escapePeaks->unChecked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                                FeatureMarkerType::EscapePeakMarker, false ) );
  
  m_comptonPeak = new WCheckBox( "Compton Peak", this );
  m_comptonPeak->setInline( false );
  m_comptonPeak->checked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                                FeatureMarkerType::ComptonPeakMarker, true ) );
  m_comptonPeak->unChecked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                                  FeatureMarkerType::ComptonPeakMarker, false ) );
  
  WContainerWidget *angleDiv = new WContainerWidget( this );
  angleDiv->addStyleClass( "AngleRow" );
  WLabel *label = new WLabel( "Angle", angleDiv );
  label->setMargin( WLength(1.8,WLength::FontEm), Wt::Left );
  label->setMargin( WLength(0.2,WLength::FontEm), Wt::Right );
  m_comptonAngle = new WSpinBox( angleDiv );
  m_comptonAngle->addStyleClass( "AngleInput" );

  //m_comptonAngle->setTextSize( 3 );
  //m_comptonAngle->setWidth( WLength( 4, WLength::FontEm ) );
  label->setBuddy( m_comptonAngle );
  m_comptonAngle->setRange( 0, 180 );
  m_comptonAngle->setValue( 180 );
  m_viewer->setComptonPeakAngle( 180 );
  m_comptonAngle->changed().connect( this, &FeatureMarkerWidget::handleComptonAngleChanged );
  m_comptonAngle->textInput().connect( this, &FeatureMarkerWidget::handleComptonAngleChanged );
  
  m_comptonEdge = new WCheckBox( "Compton Edge", this );
  m_comptonEdge->setInline( false );
  m_comptonEdge->checked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                                FeatureMarkerType::ComptonEdgeMarker, true ) );
  m_comptonEdge->unChecked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                                  FeatureMarkerType::ComptonEdgeMarker, false ) );
  
  
  m_sumPeaks = new WCheckBox( "Sum Peak", this );
  m_sumPeaks->setInline( false );
  m_sumPeaks->checked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                                FeatureMarkerType::SumPeakMarker, true ) );
  m_sumPeaks->unChecked().connect( boost::bind( &InterSpec::setFeatureMarkerOption, m_viewer,
                                                  FeatureMarkerType::SumPeakMarker, false ) );
}//init()


void FeatureMarkerWidget::handleComptonAngleChanged()
{
  const int angle = m_comptonAngle->value();
  m_viewer->setComptonPeakAngle( angle );
}//void handleComptonAngleChanged();

