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

#include <cmath>
#include <string>
#include <vector>

#include <Wt/WPainter>
#include <Wt/Chart/WAxis>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/MakeDrfChart.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
//#include "SpecUtils/UtilityFunctions.h"

using namespace std;
using namespace Wt;

namespace
{
  const int sm_energy_col = 0;
  const int sm_data_eff_col = 1;
  const int sm_data_fwhm_col = 2;
  const int sm_equation_eff_col = 3;
  const int sm_equation_fwhm_col = 4;
  
  const int sm_num_model_cols = 5;
  /** The first `sm_num_eqn_energy_rows` rows in the model will all be points
   to chart the effificency and FWHM equations.  After that there will be
   MakeDrfChart::m_datapoints.size() more rows to represent the actual data
   points.
   */
  const int sm_num_eqn_energy_rows = 125; //ToDo: customize this based on chart size...
}//namespace


MakeDrfChart::MakeDrfChart( Wt::WContainerWidget *parent )
: Wt::Chart::WCartesianChart( parent ),
  m_det_diameter( 1.0*PhysicalUnits::cm ),
  m_det_lower_energy( 0.0 ),
  m_det_upper_energy( 3000.0 ),
  m_datapoints{},
  m_fwhmEnergyUnits( EqnEnergyUnits::keV ),
  m_fwhmCoefs{},
  m_efficiencyEnergyUnits( EqnEnergyUnits::keV ),
  m_efficiencyCoefs{},
  m_textPenColor( Wt::black ),
  m_xRangeChanged()
{
  //setAutoLayoutEnabled( true );
  setLayoutSizeAware( true );
  
  setPreferredMethod( WPaintedWidget::HtmlCanvas );
  
  setType( Chart::ScatterPlot );
  //m_chart->axis(Chart::XAxis).setLabelFormat( "%.2f" );
  axis(Chart::XAxis).setScale( Chart::LinearScale );
  axis(Chart::YAxis).setLabelFormat( "%.3g" );
  axis(Chart::Y2Axis).setLabelFormat( "%.3g" );
  
  axis(Chart::XAxis).setMinimum( m_det_lower_energy );
  axis(Chart::XAxis).setMaximum( m_det_upper_energy );
  
  WStandardItemModel *m = new WStandardItemModel( this );
  setModel( m );
  m->insertColumns( 0, sm_num_model_cols );
  m->insertRows( 0, sm_num_eqn_energy_rows );
  setXSeriesColumn( sm_energy_col );
  
  Chart::WDataSeries data_eff_series( sm_data_eff_col, Chart::SeriesType::PointSeries, Chart::YAxis );
  data_eff_series.setMarker( Chart::MarkerType::CircleMarker );
  data_eff_series.setMarkerSize( 6 );
  addSeries( data_eff_series );
  
  Chart::WDataSeries eqn_eff_series( sm_equation_eff_col, Chart::SeriesType::CurveSeries, Chart::YAxis );
  eqn_eff_series.setMarker( Wt::Chart::MarkerType::NoMarker );
  addSeries( eqn_eff_series );
  
  Chart::WDataSeries data_fwhm_series( sm_data_fwhm_col, Chart::SeriesType::PointSeries, Chart::Y2Axis );
  data_fwhm_series.setMarker( Chart::MarkerType::XCrossMarker );
  data_fwhm_series.setMarkerSize( 6 );
  addSeries( data_fwhm_series );
  
  Chart::WDataSeries eqn_fwhm_series( sm_equation_fwhm_col, Chart::SeriesType::CurveSeries, Chart::Y2Axis );
  eqn_fwhm_series.setMarker( Wt::Chart::MarkerType::NoMarker );
  addSeries( eqn_fwhm_series );
  
  //axis(Chart::XAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue/* | Chart::ZeroValue*/ );
  //axis(Chart::YAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  //axis(Chart::Y2Axis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  //axis(Chart::XAxis).setRoundLimits(...)
  
  setPlotAreaPadding(0, Wt::Top);
  setPlotAreaPadding(60, Wt::Bottom);
  setPlotAreaPadding(55, Wt::Right | Wt::Left);
  
  //axis(Chart::XAxis).setTitle( "Energy (keV)" );
  setPlotAreaPadding(20, Wt::Bottom);
  
  axis(Chart::YAxis).setVisible( true );
  axis(Chart::YAxis).setTitle( "Intrinsic Eff." );
  
  axis(Chart::Y2Axis).setVisible( true );
  axis(Chart::Y2Axis).setTitle( "FWHM" );

#if( WT_VERSION >= 0x3030400 )
  axis(Wt::Chart::Y1Axis).setTitleOrientation( Wt::Vertical );
  axis(Wt::Chart::Y2Axis).setTitleOrientation( Wt::Vertical );
#endif
  
  m->setHeaderData( sm_energy_col, Wt::Horizontal, boost::any(WString("Energy (keV)")), Wt::DisplayRole );
  m->setHeaderData( sm_data_eff_col, Wt::Horizontal, boost::any(WString("Data Intrinsic Eff.")), Wt::DisplayRole );
  m->setHeaderData( sm_data_fwhm_col, Wt::Horizontal, boost::any(WString("Data FWHM")), Wt::DisplayRole );
  m->setHeaderData( sm_equation_eff_col, Wt::Horizontal, boost::any(WString("Fit Intrinsic Eff.")), Wt::DisplayRole );
  m->setHeaderData( sm_equation_fwhm_col, Wt::Horizontal, boost::any(WString("Fit FWHM")), Wt::DisplayRole );
  
  setLegendEnabled( true );
  setLegendLocation( Wt::Chart::LegendLocation::LegendInside, Wt::Top, Wt::AlignmentFlag::AlignRight );
  setLegendColumns( 1, WLength(125,WLength::Pixel) );
  
  updateDataToModel();
  updateEqnEnergyToModel();
  updateEffEquationToModel();
  updateFwhmEquationToModel();
}//MakeDrfChart(...)


MakeDrfChart::~MakeDrfChart()
{
}


void MakeDrfChart::layoutSizeChanged( int width, int height )
{
  //ToDo: customize number of energy points based on chart size
}//void layoutSizeChanged( int width, int height )


void MakeDrfChart::setTextPenColor( const Wt::WColor &color )
{
  m_textPenColor = color;
}


void MakeDrfChart::updateYAxisRange()
{
  const double xmin = axis(Chart::XAxis).minimum();
  const double xmax = axis(Chart::XAxis).maximum();
  
  double miny = 9999.9f, maxy = -9999.9f, miny2 = 9999.9f, maxy2 = -9999.9f;
  
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
  const int nrows = m->rowCount();
  
  for( int row = 0; row < nrows; ++row )
  {
    const double x = Wt::asNumber( m->data(row, sm_energy_col) );
    if( x < xmin || x > xmax || isnan(x) )
      continue;
    
    double eff, fwhm;
    if( row < sm_num_eqn_energy_rows )
    {
      eff = Wt::asNumber( m->data(row, sm_equation_eff_col) );
      fwhm = Wt::asNumber( m->data(row, sm_equation_fwhm_col) );
    }else
    {
      eff = Wt::asNumber( m->data(row, sm_data_eff_col) );
      fwhm = Wt::asNumber( m->data(row, sm_data_fwhm_col) );
    }
    
    if( !isnan(eff) )
    {
      miny = std::min( miny, eff );
      maxy = std::max( maxy, eff );
    }
    
    if( !isnan(fwhm) )
    {
      miny2 = std::min( miny2, fwhm );
      maxy2 = std::max( maxy2, fwhm );
    }
  }//for( int row = 0; row < nrows; ++row )
  
  if( maxy < miny )
  {
    axis(Chart::YAxis).setRange(0.0, 1.0);
  }else
  {
    miny = (miny < 0.1) ? 0.0 : std::floor(9*miny)/10.0;
    maxy = std::ceil(11.0*maxy)/10.0;
    axis(Chart::YAxis).setRange(miny, maxy);
  }
  
  if( maxy2 < miny2 )
  {
    axis(Chart::Y2Axis).setRange(0.0, 1.0);
  }else
  {
    miny2 = (miny2 < 0.1) ? 0.0 : std::floor(9*miny2)/10.0;
    maxy2 = std::ceil(11.0*maxy2)/10.0;
    axis(Chart::Y2Axis).setRange(miny2, maxy2);
  }
}//void updateYAxisRange()


void MakeDrfChart::updateDataToModel()
{
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
  if( m->rowCount() > sm_num_eqn_energy_rows )
    m->removeRows( sm_num_eqn_energy_rows, (m->rowCount()-sm_num_eqn_energy_rows) );
  
  const int ndata = static_cast<int>( m_datapoints.size() );
  
  if( !ndata )
    return;
  
  m->insertRows( sm_num_eqn_energy_rows, ndata );
  
  for( int i = 0; i < ndata; ++i )
  {
    const DataPoint &data = m_datapoints[i];
    const int row = i + sm_num_eqn_energy_rows;
    const double energy = data.energy;
    
    m->setData( row, sm_energy_col, energy, Wt::DisplayRole );
    
    if( data.peak_area > 0.0f && data.source_count_rate > 0.0f
        && data.livetime > 0.0f && data.distance > 0.0f && m_det_diameter > 0.0f )
    {
      const double fracSolidAngle = DetectorPeakResponse::fractionalSolidAngle( m_det_diameter, data.distance );
      const double expected = data.source_count_rate * data.livetime * fracSolidAngle;
      const double eff = data.peak_area / expected;
      
      double fracUncert2 = 0.0;
      if( data.peak_area_uncertainty > 0.0f )
        fracUncert2 += std::pow( data.peak_area_uncertainty / data.peak_area, 2.0f );
      if( data.source_count_rate_uncertainty > 0.0f )
        fracUncert2 += std::pow( data.source_count_rate_uncertainty / data.source_count_rate, 2.0f );
      
      m->setData( row, sm_data_eff_col, boost::any( eff ), Wt::DisplayRole );
      if( fracUncert2 > 0.0 )
        m->setData( row, sm_data_eff_col, boost::any( std::sqrt(fracUncert2) ), Wt::UserRole );
      m->setData( row, sm_data_eff_col, boost::any(data.peak_color), Wt::MarkerPenColorRole );
      m->setData( row, sm_data_eff_col, boost::any(data.peak_color), Wt::MarkerBrushColorRole );
      if( data.source_information.size() )
        m->setData( row, sm_data_eff_col, boost::any( WString::fromUTF8(data.source_information)), Wt::ToolTipRole );
    }//if( we can calculate a effeiciency )
    
    if( data.peak_fwhm > 0.0f )
    {
      const double fwhm = data.peak_fwhm;
      m->setData( row, sm_data_fwhm_col, boost::any( fwhm ), Wt::DisplayRole );
      if( data.peak_fwhm_uncertainty > 0.0f )
        m->setData( row, sm_data_fwhm_col, boost::any( static_cast<double>(data.peak_fwhm_uncertainty) ), Wt::UserRole );
      m->setData( row, sm_data_fwhm_col, boost::any(data.peak_color), Wt::MarkerPenColorRole );
      m->setData( row, sm_data_fwhm_col, boost::any(data.peak_color), Wt::MarkerBrushColorRole );
    }//if( we have FWHM info )
  }//for( int i = 0; i < ndata; ++i )
  
  updateYAxisRange();
}//void updateDataToModel()


void MakeDrfChart::updateEqnEnergyToModel()
{
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
  assert( m );
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  
  if( m_det_upper_energy < m_det_lower_energy )
    std::swap( m_det_lower_energy, m_det_upper_energy );
  
  if( m_det_lower_energy == m_det_upper_energy )
  {
    //Shouldnt ever happen, but JIC
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
      m->setData( row, sm_energy_col, boost::any() );
    return;
  }//if( m_det_lower_energy == m_det_upper_energy )
  
  //ToDO: we could check the first/last energies and if equal, not bother to update
  
  for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
  {
    const double energy = m_det_lower_energy + ((m_det_upper_energy * row) / (sm_num_eqn_energy_rows - 1.0f));
    m->setData( row, sm_energy_col, boost::any(energy) );
  }//for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
}//void updateEqnEnergyToModel()


void MakeDrfChart::updateEffEquationToModel()
{
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
  assert( m );
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  
  if( m_efficiencyCoefs.empty() && m->data(0, sm_equation_eff_col).empty() )
    return;
  
  if( m_efficiencyCoefs.empty() )
  {
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
      m->setData( row, sm_equation_eff_col, boost::any() );
    return;
  }//if( no equation )
  
  const float units = ((m_efficiencyEnergyUnits==EqnEnergyUnits::keV) ? 1.0f : 0.001f);
  
  for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
  {
    const float energy = static_cast<float>( m_det_lower_energy + ((m_det_upper_energy * row) / (sm_num_eqn_energy_rows - 1.0)) );
    const double eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy*units, m_efficiencyCoefs );
    
    if( std::isnan(eff) || std::isinf(eff) )
      m->setData( row, sm_equation_eff_col, boost::any() );
    else
      m->setData( row, sm_equation_eff_col, boost::any(eff) );
  }//for( loop over eqn rows )
}//void updateEffEquationToModel()


void MakeDrfChart::updateFwhmEquationToModel()
{
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
  assert( m );
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  
  if( (m_fwhmCoefs.size()!=3) && m->data(0, sm_equation_eff_col).empty() )
    return;
  
  if( m_fwhmCoefs.size() != 3 )
  {
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
      m->setData( row, sm_equation_fwhm_col, boost::any() );
    return;
  }//if( no equation )
  
  const float units = ((m_fwhmEnergyUnits==EqnEnergyUnits::keV) ? 1.0f : 0.001f);
  const auto eqnType = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
  for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
  {
    const float energy = static_cast<float>( m_det_lower_energy + ((m_det_upper_energy * row) / (sm_num_eqn_energy_rows - 1.0)) );
    const double fwhm = DetectorPeakResponse::peakResolutionFWHM( units*energy, eqnType, m_fwhmCoefs );
    m->setData( row, sm_equation_fwhm_col, boost::any(fwhm) );
  }//for( loop over eqn rows )
}//void updateFwhmEquationToModel()


void MakeDrfChart::paint( Wt::WPainter &painter, const Wt::WRectF &rectangle ) const
{
  WCartesianChart::paint( painter, rectangle );
  
  painter.save();
  
  //Draw error bars
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
  const int nrows = m->rowCount();
  
  for( int row = sm_num_eqn_energy_rows; row < nrows; ++row )
  {
    const double energy = Wt::asNumber( m->data(row,sm_energy_col,Wt::DisplayRole) );
    const double eff = Wt::asNumber( m->data(row,sm_data_eff_col,Wt::DisplayRole) );
    const double uncert = Wt::asNumber( m->data(row,sm_data_eff_col,Wt::UserRole) );
    
    if( isnan(energy) || isnan(eff) || isnan(uncert) )
      continue;
  
    const WPointF upper_uncert = mapToDevice( energy, eff + uncert );
    const WPointF lower_uncert = mapToDevice( energy, eff - uncert );
    painter.setPen( WPen(Wt::black) );
    painter.drawLine( upper_uncert, lower_uncert );
  }
  
  //Paint some grey-ish rectangles over the energy we dont have peaks for
  if( !m_efficiencyCoefs.empty() && !m_datapoints.empty() )
  {
    float mindata = 9999999.9, maxdata = -9999999.9;
    for( const auto &d : m_datapoints )
    {
      mindata = std::min( mindata, d.energy );
      maxdata = std::max( maxdata, d.energy );
    }
  
    const WPointF left_ll = mapToDevice( m_det_lower_energy, axis(Chart::YAxis).minimum() );
    const WPointF left_ur = mapToDevice( mindata, axis(Chart::YAxis).maximum() );
    
    const WPointF right_ll = mapToDevice( maxdata, axis(Chart::YAxis).minimum() );
    const WPointF right_ur = mapToDevice( m_det_upper_energy, axis(Chart::YAxis).maximum() );
    
    painter.setBrush( WBrush( WColor(123,123,123,25) ) );
    painter.setPen( WPen(PenStyle::NoPen) );
    painter.drawRect(left_ll.x(), left_ll.y(), left_ur.x() - left_ll.x(), left_ur.y() - left_ll.y() );
    painter.drawRect(right_ll.x(), right_ll.y(), right_ur.x() - right_ll.x(), right_ur.y() - right_ll.y() );
    
    //ToDo: the above actually leaves some space on the outside edges - should fix this up
  }//if( we have equation and data )
  
  painter.restore();
}//Chi2Graphic::paint(


void MakeDrfChart::paintEvent( Wt::WPaintDevice *device )
{
  //calcAndSetAxisPadding( device->height().toPixels() );
  
  Wt::Chart::WCartesianChart::paintEvent( device );
}//void paintEvent( Wt::WPaintDevice *paintDevice )


void MakeDrfChart::setFwhmCoefficients( const std::vector<float> &coeffs, const EqnEnergyUnits units )
{
  m_fwhmCoefs = coeffs;
  m_fwhmEnergyUnits = units;
  updateFwhmEquationToModel();
}//void setFwhmCoefficients( const std::vector<double> &coeffs )


void MakeDrfChart::setEfficiencyCoefficients( const std::vector<float> &coeffs, const EqnEnergyUnits units )
{
  m_efficiencyCoefs = coeffs;
  m_efficiencyEnergyUnits = units;
  updateEffEquationToModel();
}//void setEfficiencyCoefficients( const std::vector<double> &coeffs )


void MakeDrfChart::setDataPoints( const std::vector<MakeDrfChart::DataPoint> &datapoints,
                                  const float det_diameter,
                                  const float lower_energy, const float upper_energy )
{
  m_datapoints = datapoints;
  m_det_diameter = det_diameter;
  
  const bool update_xrange = ((lower_energy != m_det_lower_energy) || (upper_energy != m_det_upper_energy));
  
  if( update_xrange )
  {
    m_det_lower_energy = lower_energy;
    m_det_upper_energy = upper_energy;
    
    updateEqnEnergyToModel();
    updateEffEquationToModel();
    updateFwhmEquationToModel();
    
    axis(Chart::XAxis).setRange( lower_energy, upper_energy );
    updateYAxisRange();
  }//if( energy range changed )
  
  updateDataToModel();
  
  if( update_xrange )
    m_xRangeChanged.emit(m_det_lower_energy,m_det_upper_energy);
}//void setDataPoints(...)


void MakeDrfChart::showFwhmPoints( const bool show )
{
  series(sm_data_fwhm_col).setHidden( !show );
  series(sm_equation_fwhm_col).setHidden( !show );
  axis(Chart::Y2Axis).setVisible( show );
  setPlotAreaPadding( (show ? 55 : 10), Wt::Right );
}//void showFwhmPoints( const bool show )


void MakeDrfChart::setXRange( double lower, double upper )
{
  if( upper < lower )
    std::swap( lower, upper );
  
  if( fabs(lower-upper) < 1.0 || isnan(lower) || isnan(upper) || isinf(lower) || isinf(upper) )
    return;
  
  axis(Chart::XAxis).setRange( lower, upper );
  updateYAxisRange();
}//void setXRange( double lower, upper )


Wt::Signal<double,double> &MakeDrfChart::xRangeChanged()
{
  return m_xRangeChanged;
}
