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

#include <set>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WBorder>
#include <Wt/WString>
#include <Wt/WServer>
#include <Wt/WSpinBox>
#include <Wt/WCheckBox>
#include <Wt/WModelIndex>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/Chart/WGridData>
#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesian3DChart>
#include <Wt/Chart/WStandardColorMap>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SearchMode3DChart.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/SearchMode3DDataModel.h"

using namespace Wt;
using namespace std;

const int SearchMode3DChart::sm_maxTimeDivs = 512;
const int SearchMode3DChart::sm_maxEnergyDivs = 1024;
const int SearchMode3DChart::sm_minTimeOrEnergyDivs = 8;

SearchMode3DChart::SearchMode3DChart( InterSpec *viewer,
                                      WContainerWidget *parent )
 : WContainerWidget( parent ),
   m_viewer( viewer ),
   m_layout( nullptr ),
   m_width( 0 ),
   m_height( 0 ),
   m_loaded( false ),
   m_model( nullptr ),
   m_data( nullptr ),
   m_chart( nullptr ),
   m_inputMinEnergy( nullptr ),
   m_inputMaxEnergy( nullptr ),
   m_inputMinTime( nullptr ),
   m_inputMaxTime( nullptr ),
   m_timeDivisions( nullptr ),
   m_energyDivisions( nullptr ),
   m_loadingTxt( nullptr )
{
  init();
}//SearchMode3DChart constructor


SearchMode3DChart::~SearchMode3DChart()
{
}//SearchMode3DChart destructor


void SearchMode3DChart::init()
{
  wApp->useStyleSheet( "InterSpec_resources/SearchMode3DChart.css" );
  wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
  
  addStyleClass( "SearchMode3DChart" );
  
  setLayoutSizeAware( true );
  
  if( m_viewer )
    m_viewer->useMessageResourceBundle( "SearchMode3DChart" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>("ShowTooltips", InterSpec::instance());
  
  m_layout = new WGridLayout();
  m_layout->setContentsMargins( 9, 9, 9, 0 ); //left, top, right, bottom (default is 9 pixels on each side)
  setLayout( m_layout );
  
  m_loadingTxt = new WText( WString::tr("sm3dc-loading") );
  m_layout->addWidget( m_loadingTxt, 0, 0, AlignCenter | AlignMiddle );
  
  WContainerWidget *controlsDiv = new WContainerWidget();
  m_layout->addWidget( controlsDiv, 1, 0 );
  controlsDiv->addStyleClass( "SearchMode3DControls" );
  
  m_layout->setRowStretch( 0, 1 );
  m_layout->setColumnStretch( 0, 1 );
  
  /*
   FIXME: I cant seem to get the actual data display to properly become log, so disabling the checkbox for now
  //Creates Checkbox and connects with logScaleToggleSwitchOn and logScaleToggleSwitchOff
  WCheckBox *logScaleCheckBox = new WCheckBox( "Log Scale", controlsDiv );
  logScaleCheckBox->addStyleClass( "GridFirstRow GridFithCol" );
  logScaleCheckBox->setChecked( false );
  logScaleCheckBox->checked().connect( boost::bind( &SearchMode3DChart::setLogZ, this, true ) );
  logScaleCheckBox->unChecked().connect( boost::bind( &SearchMode3DChart::setLogZ, this, false ) );
  */
  
  //Creates a editable textbox that allows users to input the minimum energy they are interested in
  WLabel *label = new WLabel( WString::tr("sm3dc-min-energy"), controlsDiv );
  label->addStyleClass( "GridFirstRow GridFirstCol" );
  m_inputMinEnergy = new NativeFloatSpinBox( controlsDiv );
  m_inputMinEnergy->addStyleClass( "GridFirstRow GridSecondCol" );
  m_inputMinEnergy->setSpinnerHidden( true );
  label->setBuddy( m_inputMinEnergy );
  m_inputMinEnergy->valueChanged().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
    
  //Creates a editable textbox that allows users to input the maximum energy they are interested in
  label = new WLabel( WString::tr("sm3dc-max-energy"), controlsDiv );
  label->addStyleClass( "GridFirstRow GridThirdCol" );
  m_inputMaxEnergy = new NativeFloatSpinBox( controlsDiv );
  m_inputMaxEnergy->addStyleClass( "GridFirstRow GridFourthCol" );
  m_inputMaxEnergy->setSpinnerHidden( true );
  label->setBuddy( m_inputMaxEnergy );
  m_inputMaxEnergy->valueChanged().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
    
  //Creates a editable textbox that allows users to input the minimum energy they are interested in
  label = new WLabel( WString::tr("sm3dc-min-time"), controlsDiv );
  label->addStyleClass( "GridSecondRow GridFirstCol" );
  m_inputMinTime = new NativeFloatSpinBox( controlsDiv );
  m_inputMinTime->addStyleClass( "GridSecondRow GridSecondCol" );
  m_inputMinTime->setSpinnerHidden( true );
  label->setBuddy( m_inputMinTime );
  m_inputMinTime->valueChanged().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
    
  //Creates a editable textbox that allows users to input the maximum energy they are interested in
  label = new WLabel( WString::tr("sm3dc-max-time"), controlsDiv );
  label->addStyleClass( "GridSecondRow GridThirdCol" );
  m_inputMaxTime = new NativeFloatSpinBox( controlsDiv );
  m_inputMaxTime->addStyleClass( "GridSecondRow GridFourthCol" );
  m_inputMaxTime->setSpinnerHidden( true );
  label->setBuddy( m_inputMaxTime );
  m_inputMaxTime->valueChanged().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
  
  
  WText *spacer = new WText( "&nbsp;", controlsDiv );
  spacer->addStyleClass( "GridFirstRow GridFifthCol" );
  
  label = new WLabel( WString::tr("sm3dc-time-bins"), controlsDiv );
  label->addStyleClass( "GridFirstRow GridSixthCol" );
  m_timeDivisions = new WSpinBox( controlsDiv );
  m_timeDivisions->addStyleClass( "GridFirstRow GridSeventhCol" );
  m_timeDivisions->setValue( 60 );
  m_timeDivisions->setRange( sm_minTimeOrEnergyDivs, sm_maxTimeDivs );
  label->setBuddy( m_timeDivisions );
  m_timeDivisions->changed().connect( this, &SearchMode3DChart::handleNumTimeDivsChanged );
  m_timeDivisions->enterPressed().connect( this, &SearchMode3DChart::handleNumTimeDivsChanged );
  
  HelpSystem::attachToolTipOn( {label,m_timeDivisions}, WString::tr("sm3dc-tt-time-bins"),
                              showToolTips, HelpSystem::ToolTipPosition::Left );
  
  label = new WLabel( WString::tr("sm3dc-energy-bins"), controlsDiv );
  label->addStyleClass( "GridSecondRow GridSixthCol" );
  m_energyDivisions = new WSpinBox( controlsDiv );
  m_energyDivisions->addStyleClass( "GridSecondRow GridSeventhCol" );
  m_energyDivisions->setValue( 128 );
  m_energyDivisions->setRange( sm_minTimeOrEnergyDivs, sm_maxEnergyDivs );
  label->setBuddy( m_energyDivisions );
  m_energyDivisions->changed().connect( this, &SearchMode3DChart::handleNumEnergyDivsChanged );
  m_energyDivisions->enterPressed().connect( this, &SearchMode3DChart::handleNumEnergyDivsChanged );
  
  HelpSystem::attachToolTipOn( {label,m_energyDivisions}, WString::tr("sm3dc-tt-energy-bins"),
                              showToolTips, HelpSystem::ToolTipPosition::Left );
  
  m_viewer->displayedSpectrumChanged().connect(
            boost::bind( &SearchMode3DChart::newSpectralDataSet, this,
                         boost::placeholders::_1, boost::placeholders::_2,
                         boost::placeholders::_3, boost::placeholders::_4 )
  );
  
  // Disable any form widgets getting focus
  m_inputMinEnergy->setFocus(true);
  m_inputMinEnergy->setFocus(false);
  
  //If you do anything that involves the detector resolution, you might want
  //  to un-comment out the following line.
  //m_viewer->detectorChanged().connect( this, &SearchMode3DChart::updateDisplay );
}//void init()


void SearchMode3DChart::load()
{
  WContainerWidget::load();
  
  if( !m_loaded )
  {
    m_loaded = true;
    assert( !m_model );
    assert( !m_chart );
    WServer::instance()->post( wApp->sessionId(), wApp->bind( boost::bind( &SearchMode3DChart::initChart, this ) ) );
  }//if( !m_loaded )
}//void load()


void SearchMode3DChart::initChart()
{
  assert( !m_model );
  assert( !m_chart );
  
  if( m_loadingTxt )
  {
    delete m_loadingTxt;
    m_loadingTxt = nullptr;
  }
  
  setBinningLimits();
  
  m_model = new SearchMode3DDataModel( this );
  m_model->setMaxNumTimeSamples( m_timeDivisions->value() );
  m_model->setMaxNumEnergyChannels( m_energyDivisions->value() );
  
  m_chart = new Chart::WCartesian3DChart();
  m_layout->addWidget( m_chart, 0, 0 );
  m_chart->setType( Chart::ScatterPlot );
  
  m_chart->decorationStyle().setBorder( WBorder(WBorder::Solid, WBorder::Thin, Wt::black) );
  
  double w = 800.0, h = 600.0;
  if( m_viewer )
  {
    if( m_viewer->renderedWidth() > 100 )
      w = 0.9 * m_viewer->renderedWidth();
    if( m_viewer->renderedHeight() > 100 )
      h = 0.9 * m_viewer->renderedHeight();
  }//if( m_viewer )
  
  m_chart->resize( w, h );  //only an initial guess
  //m_chart->setTitle("3D Data View");
  
  //Creates X axis
  m_chart->axis(Chart::XAxis_3D).setTitle( WString::tr("sm3dc-time-axis-label") );
  m_chart->axis(Chart::XAxis_3D).setTitleOffset( 10 );
  m_chart->axis(Chart::XAxis_3D).setLabelFormat( "%1.1f" );
  m_chart->axis(Chart::XAxis_3D).setLabelBasePoint( 0 );
  m_chart->axis(Chart::XAxis_3D).setLabelAngle( 90 );
  
  //Creates Y axis
  m_chart->axis(Chart::YAxis_3D).setTitle( WString::tr("Energy (keV)") );
  m_chart->axis(Chart::YAxis_3D).setTitleOffset( 10 );
  m_chart->axis(Chart::YAxis_3D).setLabelFormat( "%.1f" );
  m_chart->axis(Chart::YAxis_3D).setLabelBasePoint( 0 );
  // m_chart->axis(Chart::YAxis_3D).setLabelInterval( 100 );
  m_chart->axis(Chart::YAxis_3D).setLabelAngle( 90 );
  
  //Creates Z axis
  m_chart->axis(Chart::ZAxis_3D).setTitle( WString::tr("Counts") );
  //m_chart->axis(Chart::ZAxis_3D).setTitle( "Counts per Channels" );
  m_chart->axis(Chart::ZAxis_3D).setTitleOffset( 20 );
  m_chart->axis(Chart::ZAxis_3D).setLabelFormat( "%.1f" );
  m_chart->axis(Chart::ZAxis_3D).setLabelBasePoint( 0 );
  // m_chart->axis(Wt::Chart::ZAxis_3D).setLabelInterval( 100 );
  
  //Create legend
  //m_chart->setLegendStyle(Wt::WFont(), Wt::WPen(), Wt::WBrush(Wt::WColor(Wt::lightGray)));
  //m_chart->setLegendEnabled(true);
  m_chart->setLegendEnabled( false );
  
  m_chart->setGridEnabled(Chart::XZ_Plane, Chart::ZAxis_3D, true);
  m_chart->setGridEnabled(Chart::YZ_Plane, Chart::ZAxis_3D, true);
  
  m_data = new Chart::WGridData( m_model );
  m_data->setTitle( WString::tr("sm3dc-counts-axis-label") );
  m_data->setType( Wt::Chart::SurfaceSeries3D );
  
  m_data->setSurfaceMeshEnabled( true );
  
  m_chart->addDataSeries( m_data );
  
  updateDisplay();
  updateRange();
  
  wApp->triggerUpdate();
}//void initChart();


void SearchMode3DChart::setLogZ( const bool log )
{
  if( !m_chart )
    return;
  
  m_chart->axis(Chart::ZAxis_3D).setScale( log ? Chart::LogScale : Chart::LinearScale );
  
  updateRange();
  
  //The client side data doesnt actually become log (but the axises do)... cant
  //  figure out how to force it.
  //m_model->dataChanged().emit( m_model->index(0, 0), m_model->index(m_model->rowCount()-1, m_model->columnCount()-1) );
  newSpectralDataSet( SpecUtils::SpectrumType::Foreground, nullptr, {}, {} );
  
  m_model->modelReset().emit();
}//void setLogZ( const bool log )


void SearchMode3DChart::setTimeLimits()
{
  const double minTime = m_model->minTime();
  const double maxTime = m_model->maxTime();
  
  m_inputMinTime->setRange( std::floor(minTime), std::ceil(maxTime) );
  m_inputMaxTime->setRange( std::floor(minTime), std::ceil(maxTime) );
  
  if( m_inputMinTime->validate() != WValidator::Valid )
    m_inputMinTime->setValue( minTime );
  
  if( m_inputMaxTime->validate() != WValidator::Valid )
    m_inputMaxTime->setValue( maxTime );
  
  if( fabs( m_inputMinTime->value() - m_inputMaxTime->value() ) < 1.0 )
  {
    m_inputMinTime->setValue( minTime );
    m_inputMaxTime->setValue( maxTime );
  }
  
  if( m_inputMinTime->value() > m_inputMaxTime->value() )
  {
    const double a = m_inputMinTime->value(), b = m_inputMaxTime->value();
    m_inputMinTime->setValue( b );
    m_inputMaxTime->setValue( a );
  }
}//void setTimeLimits()


void SearchMode3DChart::setEnergyLimits()
{
  const double minEnergy = m_model->minEnergy();
  const double maxEnergy = m_model->maxEnergy();
  
  m_inputMinEnergy->setRange( std::floor(minEnergy), std::ceil(maxEnergy) );
  m_inputMaxEnergy->setRange( std::floor(minEnergy), std::ceil(maxEnergy) );
  
  if( m_inputMinEnergy->validate() != WValidator::Valid )
    m_inputMinEnergy->setValue( minEnergy );
  
  if( m_inputMaxEnergy->validate() != WValidator::Valid )
    m_inputMaxEnergy->setValue( maxEnergy );
  
  if( fabs( m_inputMinEnergy->value() - m_inputMaxEnergy->value() ) < 1.0 )
  {
    m_inputMinEnergy->setValue( minEnergy );
    m_inputMaxEnergy->setValue( maxEnergy );
  }
  
  if( m_inputMinEnergy->value() > m_inputMaxEnergy->value() )
  {
    const double a = m_inputMinEnergy->value(), b = m_inputMaxEnergy->value();
    m_inputMinEnergy->setValue( b );
    m_inputMaxEnergy->setValue( a );
  }
}//void setEnergyLimits()


void SearchMode3DChart::setBinningLimits()
{
  auto meas = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( !meas )
    return;
  
  const int prevTimeDivs = std::min( std::max(m_timeDivisions->value(),sm_minTimeOrEnergyDivs), sm_maxTimeDivs );
  const int prevEnergyDivs = std::min( std::max(m_energyDivisions->value(),sm_minTimeOrEnergyDivs), sm_maxEnergyDivs );
  
  const int num_samples = static_cast<int>( meas->sample_numbers().size() );
  const int num_channels = static_cast<int>( meas->num_gamma_channels() );
  if( (num_samples > sm_minTimeOrEnergyDivs) && (num_channels > sm_minTimeOrEnergyDivs) )
  {
    m_timeDivisions->setRange( sm_minTimeOrEnergyDivs, std::min(num_samples,sm_maxTimeDivs) );
    m_energyDivisions->setRange( sm_minTimeOrEnergyDivs, std::min(num_channels,sm_maxEnergyDivs) );
    
    m_timeDivisions->setSingleStep( (m_timeDivisions->maximum() - m_timeDivisions->minimum())/5 );
    m_energyDivisions->setSingleStep( (m_energyDivisions->maximum() - m_energyDivisions->minimum())/5 );
    
    if( prevTimeDivs > m_timeDivisions->maximum() )
    {
      m_timeDivisions->setValue( m_timeDivisions->maximum() );
      m_model->setMaxNumTimeSamples( m_timeDivisions->value() );
    }
    
    if( prevEnergyDivs > m_energyDivisions->maximum() )
    {
      m_energyDivisions->setValue( m_energyDivisions->maximum() );
      m_model->setMaxNumEnergyChannels( m_energyDivisions->value() );
    }
  }//if( we have non-zero
}//void setBinningLimits();


void SearchMode3DChart::newSpectralDataSet( const SpecUtils::SpectrumType type,
                                           const std::shared_ptr<SpecMeas> &meas,
                                           const std::set<int> &sample_numbers,
                                           const std::vector<std::string> &detectors )
{
  if( !m_chart || (type != SpecUtils::SpectrumType::Foreground) )
    return;
  
  //Note that if the user changes the displayed sample numbers the updateDisplay()
  //  function will be called as well, so we could implement a check to see if
  //  things actually need to be re-rendered, but whatever for now.
    
  setBinningLimits();
  updateDisplay();
  updateRange();
}//void newSpectralDataSet()


void SearchMode3DChart::handleNumTimeDivsChanged()
{
  int timeDivs = std::min( std::max(m_timeDivisions->value(),sm_minTimeOrEnergyDivs), sm_maxTimeDivs );
  timeDivs = std::min( timeDivs, m_timeDivisions->maximum() );
  
  // If user entered a value outside of valid range, put it back in range
  if( timeDivs != m_timeDivisions->value() )
    m_timeDivisions->setValue( timeDivs );
  
  m_model->setMaxNumTimeSamples( timeDivs );
  updateDisplay();
  updateRange();
}//void handleNumTimeDivsChanged()


void SearchMode3DChart::handleNumEnergyDivsChanged()
{
  int energyDivs = std::min( std::max(m_energyDivisions->value(),sm_minTimeOrEnergyDivs), sm_maxEnergyDivs );
  energyDivs = std::min( energyDivs, m_energyDivisions->maximum() );
  
  // If user entered a value outside of valid range, put it back in range
  if( energyDivs != m_energyDivisions->value() )
    m_energyDivisions->setValue( energyDivs );
  
  m_model->setMaxNumEnergyChannels( energyDivs );
  updateDisplay();
  updateRange();
}//void handleNumEnergyDivsChanged()


void SearchMode3DChart::updateRange()
{
  if( !m_chart )
    return;
  
  setTimeLimits();
  setEnergyLimits();
  
  const double minenergy = m_inputMinEnergy->value();
  const double maxenergy = m_inputMaxEnergy->value();
  const double mintime = m_inputMinTime->value();
  const double maxtime = m_inputMaxTime->value();

  //Wt::Chart::WStandardColorMap *colormapUpdated =
  //                    new Wt::Chart::WStandardColorMap(minenergy,maxenergy, false);
  //m_data->setColorMap(colormapUpdated);
    
  // Update the X and Y axis to correspond to the user input
  m_chart->axis(Chart::YAxis_3D).setRange( minenergy, maxenergy );
  m_chart->axis(Chart::XAxis_3D).setRange( mintime, maxtime );
  
  std::pair<float,float> minmax_counts = m_model->minMaxCounts( mintime, maxtime, minenergy, maxenergy );
  double zmin = minmax_counts.first; // m_model->minCounts();
  double zmax = minmax_counts.second; // m_model->maxCounts();
  
  switch( m_chart->axis(Chart::ZAxis_3D).scale() )
  {
    case Chart::LinearScale:
      zmax = ((zmax < 1.0) ? 1.0 : 1.2*zmax);
      break;
      
    case Chart::LogScale:
      //Hmmm, it appears log scale on the Z-axis doesnt render correctly... maybe a bug in Wt?
      zmin = std::max( zmin, 0.1 );
      zmax = ((zmax < 1.0) ? 1.0 : 1.2*zmax);
      break;
      
    case Chart::CategoryScale:
    case Chart::DateScale:
    case Chart::DateTimeScale:
      //Shouldnt ever get here.
      break;
  }//switch( m_chart->axis(Wt::Chart::ZAxis_3D).scale() )
  
  m_chart->axis(Chart::ZAxis_3D).setRange( zmin, zmax );

  /*
   //Meh, not convinced its actually working
  const double x_range = maxtime - mintime;
  const double x_exponent = int(log(x_range));
  double x_magnitude = pow(10, x_exponent);
  if( x_magnitude <= 0.0 || x_magnitude > x_range )
    x_magnitude = 0.1*x_range;
  
  const double y_range = maxenergy - minenergy;
  const double y_exponent = int(log(y_range));
  double y_magnitude = pow(10, y_exponent);
  if( y_magnitude <= 0.0 || y_magnitude > y_range )
    y_magnitude = 0.1*y_range;
  
  m_chart->axis(Chart::XAxis_3D).setLabelInterval( x_magnitude );
  m_chart->axis(Chart::YAxis_3D).setLabelInterval( y_magnitude );
  //m_chart->axis(Chart::XAxis_3D).setLabelBasePoint( x_magnitude*ceil(mintime/x_magnitude) );
   */
}//updateRange()


void SearchMode3DChart::updateDisplay()
{
  if( !m_chart )
    return;
  
  //This function is called when first instantiating this SearchMode3DChart, as
  //  well as whenever the data the user is viewing is changed
  m_model->update( m_viewer );

  //You might want to customize the color map, perhaps based on sigma above
  //  background, with color generated by http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/
  Chart::WStandardColorMap *colormap =
    new Chart::WStandardColorMap( m_data->minimum(Wt::Chart::ZAxis_3D),
                                      m_data->maximum(Wt::Chart::ZAxis_3D),
                                      true );
  m_data->setColorMap( colormap );
  
  //TODO modify SearchMode3DDataModel to generate the colors.
}//void updateDisplay()


void SearchMode3DChart::layoutSizeChanged( int width, int height )
{
  m_width = width;
  m_height = height;
  
  //TODO: adjust number of time and energy samples based on the resolution.
  //m_model->setMaxNumTimeSamples()
  //m_model->setMaxNumEnergyChannels();
}//void layoutSizeChanged( int width, int height )

