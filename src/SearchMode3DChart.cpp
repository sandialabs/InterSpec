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
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WColor>
#include <Wt/WBorder>
#include <Wt/WString>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WModelIndex>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WDoubleSpinBox>
#include <Wt/Chart/WGridData>
#include <Wt/WContainerWidget>
#include <Wt/WDoubleValidator>
#include <Wt/WAbstractTableModel>
#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WCartesian3DChart>
#include <Wt/Chart/WStandardColorMap>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SearchMode3DChart.h"
#include "InterSpec/SearchMode3DDataModel.h"

using namespace Wt;
using namespace std;


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
   m_logScaleCheckBox( nullptr ),
   m_inputMinEnergy( nullptr ),
   m_inputMaxEnergy( nullptr ),
   m_inputMinTime( nullptr ),
   m_inputMaxTime( nullptr )
{
  init();
}//SearchMode3DChart constructor


SearchMode3DChart::~SearchMode3DChart()
{
}//SearchMode3DChart destructor


void SearchMode3DChart::init()
{
  wApp->useStyleSheet( "InterSpec_resources/SearchMode3DChart.css" );
  
  addStyleClass( "SearchMode3DChart" );
  
  setLayoutSizeAware( true );
  
  m_layout = new WGridLayout();
  setLayout( m_layout );
  
  m_layout->addWidget( new WText( "Loading..." ), 0, 0, 1, 5 );
  
  //Creates Checkbox and connects with logScaleToggleSwitchOn and logScaleToggleSwitchOff
  m_logScaleCheckBox = new WCheckBox( "Log Scale" );
  m_logScaleCheckBox->setChecked( false );
  m_logScaleCheckBox->checked().connect( boost::bind( &SearchMode3DChart::setLogZ, this, true ) );
  m_logScaleCheckBox->unChecked().connect( boost::bind( &SearchMode3DChart::setLogZ, this, false ) );
  m_layout->addWidget( m_logScaleCheckBox, 1, 4, AlignRight );
  
  //Creates a editable textbox that allows users to input the minimum energy they are interested in
  WLabel *label = new WLabel( "Min Energy" );
  m_layout->addWidget( label, 2, 0 );
  m_inputMinEnergy = new WDoubleSpinBox();
  m_inputMinEnergy->setWidth( WLength(10,WLength::FontEx) );
  label->setBuddy( m_inputMinEnergy );
  m_inputMinEnergy->changed().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
  //m_inputMinEnergy->enterPressed().connect( boost::bind(&SearchMode3DChart::updateDisplay, this) );
  //m_inputMinEnergy->blurred().connect( boost::bind(&SearchMode3DChart::updateDisplay, this) );
  //WDoubleValidator *val = new WDoubleValidator( container );
  //inputMinEnergy->setInputMask("999999");
  m_layout->addWidget( m_inputMinEnergy, 2, 1 );
    
  //Creates a editable textbox that allows users to input the maximum energy they are interested in
  //m_inputMaxEnergy = new Wt::WDoubleSpinBox( container );
  label = new WLabel( "Max Energy" );
  m_layout->addWidget( label, 2, 2 );
  m_inputMaxEnergy = new WDoubleSpinBox();
  m_inputMaxEnergy->setWidth( WLength(10,WLength::FontEx) );
  label->setBuddy( m_inputMaxEnergy );
  m_inputMaxEnergy->changed().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
  m_layout->addWidget( m_inputMaxEnergy, 2, 3 );
    
  //Creates a editable textbox that allows users to input the minimum energy they are interested in
  label = new WLabel( "Min Time" );
  m_layout->addWidget( label, 3, 0 );
  m_inputMinTime = new WDoubleSpinBox();
  m_inputMinTime->setWidth( WLength(10,WLength::FontEx) );
  label->setBuddy( m_inputMinTime );
  m_inputMinTime->changed().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
  m_layout->addWidget( m_inputMinTime, 3, 1 );
    
  //Creates a editable textbox that allows users to input the maximum energy they are interested in
  //m_inputMaxEnergy = new Wt::WDoubleSpinBox( container );
  label = new WLabel( "Max Time" );
  m_layout->addWidget( label, 3, 2 );
  m_inputMaxTime = new WDoubleSpinBox();
  m_inputMaxTime->setWidth( WLength(10,WLength::FontEx) );
  label->setBuddy( m_inputMaxTime );
  m_inputMaxTime->changed().connect( boost::bind(&SearchMode3DChart::updateRange, this) );
  m_layout->addWidget( m_inputMaxTime, 3, 3 );
  
  //If this chart window is opened, and the user changes the spectrum, lets
  //  update the chart.
  //TODO: if the user changes the displayed sample numbers, or the background
  //      or secondary spectrum, the updateDisplay() function will be called as
  //      well.
  m_viewer->displayedSpectrumChanged().connect( this, &SearchMode3DChart::newSpectralDataSet );
  
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
  
  m_model = new SearchMode3DDataModel( this );
  m_model->setMaxNumTimeSamples( 60 );
  m_model->setMaxNumEnergyChannels( 128 );
  
  m_chart = new Chart::WCartesian3DChart();
  m_layout->addWidget( m_chart, 0, 0, 1, 5 );
  m_layout->setRowStretch( 0, 1 );
  m_layout->setColumnStretch( 4, 1 );
  
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
  m_chart->axis(Chart::XAxis_3D).setTitle( "Time (seconds)" );
  m_chart->axis(Chart::XAxis_3D).setTitleOffset( 10 );
  m_chart->axis(Chart::XAxis_3D).setLabelFormat( "%1.1f" );
  m_chart->axis(Chart::XAxis_3D).setLabelBasePoint( 0 );
  m_chart->axis(Chart::XAxis_3D).setLabelAngle( 90 );
  
  //Creates Y axis
  m_chart->axis(Chart::YAxis_3D).setTitle( "Energy (keV)" );
  m_chart->axis(Chart::YAxis_3D).setTitleOffset( 10 );
  m_chart->axis(Chart::YAxis_3D).setLabelFormat( "%.1f" );
  m_chart->axis(Chart::YAxis_3D).setLabelBasePoint( 0 );
  // m_chart->axis(Chart::YAxis_3D).setLabelInterval( 100 );
  m_chart->axis(Chart::YAxis_3D).setLabelAngle( 90 );
  
  //Creates Z axis
  m_chart->axis(Chart::ZAxis_3D).setTitle( "Counts" );
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
  m_data->setTitle( "Counts per Channel" );
  m_data->setType( Wt::Chart::SurfaceSeries3D );
  
  m_data->setSurfaceMeshEnabled( true );
  
  m_chart->addDataSeries( m_data );
  
  updateDisplay();
  updateRange();
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
  newSpectralDataSet();
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


void SearchMode3DChart::newSpectralDataSet()
{
  if( !m_chart )
    return;
  
  //Note that if the user changes the displayed sample
  //  numbers, or the background/secondary spectrum, the updateDisplay()
  //  function will be called as well, so we could implement a check to see if
  //  things actually need to be re-rendered, but whatever for now.
    
  updateDisplay();
  updateRange();
}//void newSpectralDataSet()


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

