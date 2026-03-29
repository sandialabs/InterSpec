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
#include <memory>
#include <string>
#include <vector>

#include <Wt/WPainter.h>
#include <Wt/Chart/WAxis.h>
#include <Wt/cpp17/any.hpp>
#include <Wt/Chart/WCartesianChart.h>
#include <Wt/WStandardItemModel.h>
#include <Wt/Chart/WAbstractChartModel.h>
#include <Wt/Chart/WStandardChartProxyModel.h>

#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/MakeDrfChart.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;

namespace
{
  const int sm_energy_col = 0;
  const int sm_data_eff_col = 1;
  const int sm_data_fwhm_col = 2;
  const int sm_equation_eff_col = 3;
  const int sm_equation_fwhm_col = 4;
  const int sm_equation_eff_neg_uncert_col = 5; //Currently set data for this, but not displaying pending validation of values
  const int sm_equation_eff_pos_uncert_col = 6; //Currently set data for this, but not displaying pending validation of values
  
  const int sm_num_model_cols = 7;
  /** The first `sm_num_eqn_energy_rows` rows in the model will all be points
   to chart the effificency and FWHM equations.  After that there will be
   MakeDrfChart::m_datapoints.size() more rows to represent the actual data
   points.
   */
  const int sm_num_eqn_energy_rows = 125; //ToDo: customize this based on chart size...
}//namespace


DrfChartHolder::DrfChartHolder( std::unique_ptr<MakeDrfChart> chart )
  : WContainerWidget(),
    m_chart( chart.get() )
{
  setLayoutSizeAware( true );
  addWidget( std::move( chart ) );
}
  
DrfChartHolder::~DrfChartHolder()
{
};
  
void DrfChartHolder::layoutSizeChanged( int width, int height )
{
  m_chart->resize( width, height );
}




MakeDrfChart::MakeDrfChart()
: Wt::Chart::WCartesianChart(),
  m_det_diameter( 1.0*PhysicalUnits::cm ),
  m_det_setback( 0.0f ),
  m_det_lower_energy( 0.0 ),
  m_det_upper_energy( 3000.0 ),
  m_datapoints{},
  m_fwhmEnergyUnits( EqnEnergyUnits::keV ),
  m_fwhmCoefs{},
  m_fwhmCoefUncerts{},
  m_efficiencyEnergyUnits( EqnEnergyUnits::keV ),
  m_efficiencyCoefs{},
  m_efficiencyCoefUncerts{},
  m_fwhmEqnType( FwhmCoefType::Gadras ),
  m_xRangeChanged(),
  m_chartMarginBrush(),
  m_textPen( WColor(Wt::StandardColor::Black) )
{
  //setAutoLayoutEnabled( true );
  setLayoutSizeAware( true );
  
  setPreferredMethod( Wt::RenderMethod::HtmlCanvas );
  
  InterSpec::instance()->useMessageResourceBundle( "MakeDrf" );
    
  setType( Chart::ChartType::Scatter );
  //m_chart->axis(Chart::Axis::X).setLabelFormat( "%.2f" );
  axis(Chart::Axis::X).setScale( Chart::AxisScale::Linear );
  axis(Chart::Axis::Y).setLabelFormat( "%.3g" );
  axis(Chart::Axis::Y2).setLabelFormat( "%.3g" );

  axis(Chart::Axis::X).setMinimum( m_det_lower_energy );
  axis(Chart::Axis::X).setMaximum( m_det_upper_energy );
  
  auto itemModel = std::make_shared<WStandardItemModel>();
  auto proxyModel = std::make_shared<Chart::WStandardChartProxyModel>( itemModel );
  setModel( proxyModel );
  WStandardItemModel *m = itemModel.get();
    
  m->insertColumns( 0, sm_num_model_cols );
  m->insertRows( 0, sm_num_eqn_energy_rows );
  setXSeriesColumn( sm_energy_col );
  
  {
    auto data_eff_series = std::make_unique<Chart::WDataSeries>( sm_data_eff_col, Chart::SeriesType::Point, Chart::Axis::Y );
    data_eff_series->setMarker( Chart::MarkerType::Circle );
    data_eff_series->setMarkerSize( 6 );
    addSeries( std::move(data_eff_series) );
  }
  {
    auto eqn_eff_series = std::make_unique<Chart::WDataSeries>( sm_equation_eff_col, Chart::SeriesType::Curve, Chart::Axis::Y );
    eqn_eff_series->setMarker( Wt::Chart::MarkerType::None );
    addSeries( std::move(eqn_eff_series) );
  }
  
  //For the moment we wont show the errors since I havent verified they are actually reasonable and correct
  /*
  Chart::WDataSeries eqn_eff_neg_uncert_series( sm_equation_eff_neg_uncert_col, Chart::SeriesType::Curve, Chart::Axis::Y );
  eqn_eff_neg_uncert_series.setMarker( Wt::Chart::MarkerType::None );
  eqn_eff_neg_uncert_series.setPen( WPen(Wt::green) );
  addSeries( eqn_eff_neg_uncert_series );
  
  Chart::WDataSeries eqn_eff_pos_uncert_series( sm_equation_eff_pos_uncert_col, Chart::SeriesType::Curve, Chart::Axis::Y );
  eqn_eff_pos_uncert_series.setMarker( Wt::Chart::MarkerType::None );
  eqn_eff_pos_uncert_series.setPen( WPen(Wt::red) );
  addSeries( eqn_eff_pos_uncert_series );
  */
  
  
  {
    auto data_fwhm_series = std::make_unique<Chart::WDataSeries>( sm_data_fwhm_col, Chart::SeriesType::Point, Chart::Axis::Y2 );
    data_fwhm_series->setMarker( Chart::MarkerType::XCross );
    data_fwhm_series->setMarkerSize( 6 );
    addSeries( std::move(data_fwhm_series) );
  }
  {
    auto eqn_fwhm_series = std::make_unique<Chart::WDataSeries>( sm_equation_fwhm_col, Chart::SeriesType::Curve, Chart::Axis::Y2 );
    eqn_fwhm_series->setMarker( Wt::Chart::MarkerType::None );
    addSeries( std::move(eqn_fwhm_series) );
  }
  
  setPlotAreaPadding(0, Wt::Side::Top);
  setPlotAreaPadding(55, Wt::Side::Right | Wt::Side::Left);
  //axis(Chart::Axis::X).setTitle( "Energy (keV)" );
  setPlotAreaPadding(25, Wt::Side::Bottom);
  
  axis(Chart::Axis::Y).setVisible( true );
  axis(Chart::Axis::Y).setTitle( WString::tr("md-geom-intrinsic-eff") );

  axis(Chart::Axis::Y2).setVisible( true );
  axis(Chart::Axis::Y2).setTitle( WString::tr("FWHM") );

  axis(Wt::Chart::Axis::Y1).setTitleOrientation( Wt::Orientation::Vertical );
  axis(Wt::Chart::Axis::Y2).setTitleOrientation( Wt::Orientation::Vertical );

  setAxisPadding( 0 );
  axis(Chart::Axis::Y1).setMargin( 0 );
  axis(Chart::Axis::Y2).setMargin( 0 );
  
  m->setHeaderData( sm_energy_col, Wt::Orientation::Horizontal, Wt::cpp17::any( WString::tr("Energy (keV)") ), Wt::ItemDataRole::Display );
  m->setHeaderData( sm_data_eff_col, Wt::Orientation::Horizontal, Wt::cpp17::any(WString::tr("md-chart-data-intrinsic-eff-label")), Wt::ItemDataRole::Display );
  m->setHeaderData( sm_data_fwhm_col, Wt::Orientation::Horizontal, Wt::cpp17::any(WString::tr("md-chart-data-fwhm-label")), Wt::ItemDataRole::Display );
  m->setHeaderData( sm_equation_eff_col, Wt::Orientation::Horizontal, Wt::cpp17::any(WString::tr("md-chart-fit-intrinsic-eff-label")), Wt::ItemDataRole::Display );
  m->setHeaderData( sm_equation_eff_neg_uncert_col, Wt::Orientation::Horizontal, Wt::cpp17::any(WString::tr("md-chart-intrinsic-eff-plus-label")), Wt::ItemDataRole::Display );
  m->setHeaderData( sm_equation_eff_pos_uncert_col, Wt::Orientation::Horizontal, Wt::cpp17::any(WString::tr("md-chart-intrinsic-eff-minus-label")), Wt::ItemDataRole::Display );
  m->setHeaderData( sm_equation_fwhm_col, Wt::Orientation::Horizontal, Wt::cpp17::any(WString::tr("md-chart-fit-fwhm-label")), Wt::ItemDataRole::Display );
  
  setLegendEnabled( true );
  setLegendLocation( Wt::Chart::LegendLocation::Inside, Wt::Side::Top, Wt::AlignmentFlag::Right );
  setLegendColumns( 1, WLength(125,WLength::Unit::Pixel) );
  
  updateDataToModel();
  updateEqnEnergyToModel();
  updateEffEquationToModel();
  updateFwhmEquationToModel();
  
  auto viewer = InterSpec::instance();
  if( viewer )  //should always be true, but JIC
  {
    updateColorTheme( viewer->getColorTheme() );
    viewer->colorThemeChanged().connect( this, &MakeDrfChart::updateColorTheme );
  }
}//MakeDrfChart(...)


MakeDrfChart::~MakeDrfChart()
{
}


void MakeDrfChart::layoutSizeChanged( int width, int height )
{
  //ToDo: customize number of energy points based on chart size
}//void layoutSizeChanged( int width, int height )


void MakeDrfChart::updateYAxisRange()
{
  const double xmin = axis(Chart::Axis::X).minimum();
  const double xmax = axis(Chart::Axis::X).maximum();

  double miny = 9999.9f, maxy = -9999.9f, miny2 = 9999.9f, maxy2 = -9999.9f;

  auto proxyModel = std::dynamic_pointer_cast<Chart::WStandardChartProxyModel>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel().get() : nullptr;
  
  assert( m );
  if( !m )
    return;

  const int nrows = m->rowCount();
  
  for( int row = 0; row < nrows; ++row )
  {
    const double x = Wt::asNumber( m->data(row, sm_energy_col) );
    if( x < xmin || x > xmax || IsNan(x) || IsInf(x) )
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
    
    //Let the data efficiency be whatever, but clamp equation eff from 0 to 5.
    if( !IsNan(eff) && !IsInf(eff)
        && ((row >= sm_num_eqn_energy_rows) || (eff >= 0.0 && eff < 5.0)) )
    {
      miny = std::min( miny, eff );
      maxy = std::max( maxy, eff );
    }
    
    if( !IsNan(fwhm) && !IsInf(fwhm) )
    {
      miny2 = std::min( miny2, fwhm );
      maxy2 = std::max( maxy2, fwhm );
    }
  }//for( int row = 0; row < nrows; ++row )
  
  if( maxy <= miny )
  {
    axis(Chart::Axis::Y).setRange(0.0, 1.0);
  }else
  {
    miny = (miny < 0.1) ? 0.0 : std::floor(9*miny)/10.0;
    maxy = ((11.0*maxy) > 1.5) ? (std::ceil(11.0*maxy)/10.0) : 1.15*maxy;
    axis(Chart::Axis::Y).setRange(miny, maxy);
  }

  if( maxy2 < miny2 )
  {
    axis(Chart::Axis::Y2).setRange(0.0, 1.0);
  }else
  {
    miny2 = (miny2 < 0.1) ? 0.0 : std::floor(9*miny2)/10.0;
    maxy2 = std::ceil(11.0*maxy2)/10.0;
    axis(Chart::Axis::Y2).setRange(miny2, maxy2);
  }
}//void updateYAxisRange()


void MakeDrfChart::updateDataToModel()
{
  auto proxyModel = std::dynamic_pointer_cast<Chart::WStandardChartProxyModel>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel().get() : nullptr;
  
  assert( m );
  if( !m )
    return;
  
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
    
    m->setData( row, sm_energy_col, energy, Wt::ItemDataRole::Display );

    if( (data.peak_area > 0.0f) && (data.source_count_rate > 0.0f)
        && (data.livetime > 0.0f) && ((m_det_diameter > 0.0f) || (data.distance < 0.0)) )
    {
      const double fracSolidAngle
           = (data.distance < 0.0) ? 1.0
                                   : DetectorPeakResponse::fractionalSolidAngle( m_det_diameter, data.distance + m_det_setback );
      const double expected = data.source_count_rate * data.livetime * fracSolidAngle;
      const double eff = data.peak_area / expected;

      double fracUncert2 = 0.0;
      if( data.peak_area_uncertainty > 0.0f )
        fracUncert2 += std::pow( data.peak_area_uncertainty / data.peak_area, 2.0f );
      if( data.source_count_rate_uncertainty > 0.0f )
        fracUncert2 += std::pow( data.source_count_rate_uncertainty / data.source_count_rate, 2.0f );
      const double eff_uncert = eff * std::sqrt(fracUncert2);

      m->setData( row, sm_data_eff_col, Wt::cpp17::any(eff), Wt::ItemDataRole::Display );
      if( eff_uncert > 0.0 )
        m->setData( row, sm_data_eff_col, Wt::cpp17::any(eff_uncert), Wt::ItemDataRole::User );
      m->setData( row, sm_data_eff_col, Wt::cpp17::any(data.peak_color), Wt::ItemDataRole::MarkerPenColor );
      m->setData( row, sm_data_eff_col, Wt::cpp17::any(data.peak_color), Wt::ItemDataRole::MarkerBrushColor );
      if( data.source_information.size() )
        m->setData( row, sm_data_eff_col, Wt::cpp17::any( WString::fromUTF8(data.source_information)), Wt::ItemDataRole::ToolTip );
    }//if( we can calculate a efficiency )

    if( data.peak_fwhm > 0.0f )
    {
      const double fwhm = data.peak_fwhm;
      m->setData( row, sm_data_fwhm_col, Wt::cpp17::any( fwhm ), Wt::ItemDataRole::Display );
      if( data.peak_fwhm_uncertainty > 0.0f )
        m->setData( row, sm_data_fwhm_col, Wt::cpp17::any( static_cast<double>(data.peak_fwhm_uncertainty) ), Wt::ItemDataRole::User );
      m->setData( row, sm_data_fwhm_col, Wt::cpp17::any(data.peak_color), Wt::ItemDataRole::MarkerPenColor );
      m->setData( row, sm_data_fwhm_col, Wt::cpp17::any(data.peak_color), Wt::ItemDataRole::MarkerBrushColor );
    }//if( we have FWHM info )
  }//for( int i = 0; i < ndata; ++i )
  
  updateYAxisRange();
}//void updateDataToModel()


void MakeDrfChart::updateEqnEnergyToModel()
{
  auto proxyModel = std::dynamic_pointer_cast<Chart::WStandardChartProxyModel>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel().get() : nullptr;
  
  assert( m );
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  if( !m )
    return;
  
  if( m_det_upper_energy < m_det_lower_energy )
    std::swap( m_det_lower_energy, m_det_upper_energy );
  
  if( m_det_lower_energy == m_det_upper_energy )
  {
    //Shouldnt ever happen, but JIC
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
      m->setData( row, sm_energy_col, Wt::cpp17::any() );
    return;
  }//if( m_det_lower_energy == m_det_upper_energy )
  
  //ToDO: we could check the first/last energies and if equal, not bother to update
  
  for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
  {
    const double energy = m_det_lower_energy + ((m_det_upper_energy * row) / (sm_num_eqn_energy_rows - 1.0f));
    m->setData( row, sm_energy_col, Wt::cpp17::any(energy) );
  }//for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
}//void updateEqnEnergyToModel()


void MakeDrfChart::updateEffEquationToModel()
{
  auto proxyModel = std::dynamic_pointer_cast<Chart::WStandardChartProxyModel>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel().get() : nullptr;
  
  assert( m );
  if( !m )
    return;
  
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  
  if( m_efficiencyCoefs.empty() )
  {
    if( !Wt::cpp17::any_has_value(m->data(0, sm_equation_eff_col)) )
      return;
    
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
    {
      m->setData( row, sm_equation_eff_col, Wt::cpp17::any() );
      m->setData( row, sm_equation_eff_pos_uncert_col, Wt::cpp17::any() );
      m->setData( row, sm_equation_eff_neg_uncert_col, Wt::cpp17::any() );
    }
    return;
  }//if( no equation )
  
  const float units = ((m_efficiencyEnergyUnits==EqnEnergyUnits::keV) ? 1.0f : 0.001f);
  
  for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
  {
    const float energy = static_cast<float>( m_det_lower_energy + ((m_det_upper_energy * row) / (sm_num_eqn_energy_rows - 1.0)) );
    const double eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy*units, m_efficiencyCoefs );
    
    if( IsNan(eff) || IsInf(eff) )
    {
      m->setData( row, sm_equation_eff_col, Wt::cpp17::any() );
      m->setData( row, sm_equation_eff_pos_uncert_col, Wt::cpp17::any() );
      m->setData( row, sm_equation_eff_neg_uncert_col, Wt::cpp17::any() );
    }else
    {
      m->setData( row, sm_equation_eff_col, Wt::cpp17::any(eff) );
      
      Wt::cpp17::any lowerval, upperval;
      
      if( m_efficiencyCoefUncerts.size() == m_efficiencyCoefs.size() )
      {
        double neguncert = 0.0, posuncert = 0.0;
        
        for( size_t i = 0; i < m_efficiencyCoefUncerts.size(); ++i )
        {
          const float orig = m_efficiencyCoefs[i];
          m_efficiencyCoefs[i] = orig + 2.0f*m_efficiencyCoefUncerts[i];
          const double plus = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy*units, m_efficiencyCoefs ) - eff;
          m_efficiencyCoefs[i] = orig - 2.0f*m_efficiencyCoefUncerts[i];
          const double minus = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy*units, m_efficiencyCoefs ) - eff;
          m_efficiencyCoefs[i] = orig;
          if( plus > 0.0 )
            posuncert += plus*plus;
          else
            neguncert += plus*plus;
          if( minus > 0.0 )
            posuncert += minus*minus;
          else
            neguncert += minus*minus;
        }//for( size_t i = 0; i < m_efficiencyCoefUncerts.size(); ++i )
        
        posuncert = sqrt(posuncert);
        neguncert = sqrt(neguncert);
        
        if( IsNan(posuncert) || IsInf(posuncert) )
          upperval = Wt::cpp17::any( eff + posuncert );
        if( IsNan(neguncert) || IsInf(neguncert) )
          lowerval = Wt::cpp17::any( eff - neguncert );
      }//if( we have uncertainties )
      
      m->setData( row, sm_equation_eff_pos_uncert_col, upperval );
      m->setData( row, sm_equation_eff_neg_uncert_col, lowerval );
    }//if( valid eff value ) / else
  }//for( loop over eqn rows )
}//void updateEffEquationToModel()


void MakeDrfChart::updateFwhmEquationToModel()
{
  auto proxyModel = std::dynamic_pointer_cast<Chart::WStandardChartProxyModel>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel().get() : nullptr;

  assert( m );
  if( !m )
    return;
  
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  
  if( ((m_fwhmEqnType==FwhmCoefType::Gadras) && m_fwhmCoefs.size() != 3)
     || ((m_fwhmEqnType==FwhmCoefType::SqrtEnergyPlusInverse) && m_fwhmCoefs.size() != 3)
     || ((m_fwhmEqnType==FwhmCoefType::ConstantPlusSqrtEnergy) && m_fwhmCoefs.size() != 2)
     || (m_fwhmEqnType==FwhmCoefType::SqrtEqn && m_fwhmCoefs.size() < 1) )
  {
    if( !Wt::cpp17::any_has_value(m->data(0, sm_equation_fwhm_col)) )  //
      return;
    
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
      m->setData( row, sm_equation_fwhm_col, Wt::cpp17::any() );
    return;
  }//if( no equation )
  
  const float units = 1.0f; //For the moment always using keV  ((m_fwhmEnergyUnits==EqnEnergyUnits::keV) ? 1.0f : 0.001f);
  DetectorPeakResponse::ResolutionFnctForm eqnType;
  switch( m_fwhmEqnType )
  {
    case FwhmCoefType::Gadras:
      eqnType = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
      break;
      
    case FwhmCoefType::SqrtEnergyPlusInverse:
      eqnType = DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse;
      break;
    
    case FwhmCoefType::ConstantPlusSqrtEnergy:
      eqnType = DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy;
      break;
      
    case FwhmCoefType::SqrtEqn:
      eqnType = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      break;
  }//switch( m_fwhmEqnType )

  for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
  {
    const float energy = static_cast<float>( m_det_lower_energy + ((m_det_upper_energy * row) / (sm_num_eqn_energy_rows - 1.0)) );
    const double fwhm = DetectorPeakResponse::peakResolutionFWHM( units*energy, eqnType, m_fwhmCoefs );
    m->setData( row, sm_equation_fwhm_col, Wt::cpp17::any(fwhm) );
  }//for( loop over eqn rows )
}//void updateFwhmEquationToModel()


void MakeDrfChart::updateColorTheme( std::shared_ptr<const ColorTheme> theme )
{
  if( !theme )
    return; //Should reset colors to default, but whatever for now
  
  
  if( theme->spectrumAxisLines.isDefault() )
    m_textPen = WPen( WColor(Wt::StandardColor::Black) );
  else
    m_textPen = WPen( theme->spectrumAxisLines );
  setTextPen( m_textPen );
  
  if( theme->spectrumAxisLines.isDefault() )
  {
    axis(Chart::Axis::X).setPen( WPen(Wt::StandardColor::Black) );
    axis(Chart::Axis::Y).setPen( WPen(Wt::StandardColor::Black) );
    axis(Chart::Axis::Y2).setPen( WPen(Wt::StandardColor::Black) );
  }else
  {
    axis(Chart::Axis::X).setPen( WPen(theme->spectrumAxisLines) );
    axis(Chart::Axis::Y).setPen( WPen(theme->spectrumAxisLines) );
    axis(Chart::Axis::Y2).setPen( WPen(theme->spectrumAxisLines) );
  }

  if( theme->spectrumChartMargins.isDefault() )
    m_chartMarginBrush = WBrush();
  else
    m_chartMarginBrush = WBrush(theme->spectrumChartMargins);

  if( theme->spectrumChartBackground.isDefault() )
    setBackground( WBrush() );
  else
    setBackground( WBrush(theme->spectrumChartBackground) );
  
  auto setSeriesColor = []( Wt::Chart::WDataSeries &s, const WColor &color ){
    WPen pen = s.pen();
    WBrush brush = s.brush();
    pen.setColor( color );
    brush.setColor( color );
    s.setPen( pen );
    s.setBrush( brush );
  };//setSeriesColor lambda
  
  // Data series are drawn with a marker that is the same color as the peak, so we dont need to set
  //  those colors, but we do need to set the colors of the equation line.
  //Wt::Chart::WDataSeries &dataEffSeries = series(sm_data_eff_col);
  //Wt::Chart::WDataSeries &dataFwhmSeries = series(sm_data_fwhm_col);
  
  Wt::Chart::WDataSeries &eqnEffSeries = series(sm_equation_eff_col);
  WColor effcolor = theme->foregroundLine;
  if( effcolor.isDefault() )
    effcolor = WColor(Wt::StandardColor::Black);
  setSeriesColor( eqnEffSeries, effcolor );
  
  Wt::Chart::WDataSeries &eqnFwhmSeries = series(sm_equation_fwhm_col);
  WColor fwhmcolor = theme->backgroundLine;
  if( fwhmcolor.isDefault() )
    fwhmcolor = WColor(Wt::StandardColor::Gray);
  setSeriesColor( eqnFwhmSeries, fwhmcolor );
  
  //Efficiency series not yet implemented
  //WDataSeries &negEffUncertSeries = series(sm_equation_eff_neg_uncert_col);
  //WDataSeries &posEffUncertSeries = series(sm_equation_eff_pos_uncert_col);
  
  update(); //trigger re-render
}//void updateColorTheme();


void MakeDrfChart::paint( Wt::WPainter &painter, const Wt::WRectF &rectangle ) const
{
  painter.save();
  
  painter.setPen( m_textPen ); //TODO: Should make so legend text will be correct color - but doesnt
  
  WCartesianChart::paint( painter, rectangle );
  
  /// Start color theme support - this code is largely a copy of DecayActivityChart - so if you improve here, improve there
  /// TODO: should refactor into a common base class to support color theme
  const double height = painter.window().height();
  const double width = painter.window().width();
  
  auto plotArea = [&]() -> WRectF {
    int w, h;
    if( rectangle.isNull() || rectangle.isEmpty() )
    {
      w = static_cast<int>( width );
      h = static_cast<int>( height );
    }else
    {
      w = static_cast<int>( rectangle.width() );
      h = static_cast<int>( rectangle.height() );
    }
    
    const int padLeft = plotAreaPadding(Wt::Side::Left);
    const int padRight = plotAreaPadding(Wt::Side::Right);
    const int padTop = plotAreaPadding(Wt::Side::Top);
    const int padBottom = plotAreaPadding(Wt::Side::Bottom);

    WRectF area;
    if( orientation() == Wt::Orientation::Vertical )
      area = WRectF( padLeft, padTop, std::max(10, w - padLeft - padRight), std::max(10, h - padTop - padBottom) );
    else
      area = WRectF( padTop, padRight, std::max(10, w - padTop - padBottom), std::max(10, h - padRight - padLeft) );
    
    return area;
  };//plotArea
  
  auto paintMargins = [&](){
    auto area = plotArea();
    const double rx = area.x();
    const double ry = area.y();
    const double rw = area.width();
    const double rh = area.height();
    if( ry > 0 )  //Top strip
      painter.fillRect( hv(WRectF(rx,0,rw,ry)), m_chartMarginBrush );
    if( (ry+rh) < rectangle.height() ) //Bottom strip
      painter.fillRect( hv(WRectF(rx,ry+rh,rw,rectangle.height()-ry)), m_chartMarginBrush );
    if( rx > 0 )
      painter.fillRect( hv(WRectF(0,0,rx,rectangle.height())), m_chartMarginBrush );
    if( (rx+rw) < rectangle.width() )
      painter.fillRect( hv(WRectF(rx+rw,0,rectangle.width()-rw-rx,rectangle.height())), m_chartMarginBrush );
  };
  
  if( background().style() != Wt::BrushStyle::None && m_chartMarginBrush.style() != Wt::BrushStyle::None )
  {
    paintMargins();
    painter.fillRect( hv(plotArea()), background() );
  }else if( background().style() != Wt::BrushStyle::None )
  {
    painter.fillRect( hv(rectangle), background() );
  }else if( m_chartMarginBrush.style() != Wt::BrushStyle::None )
  {
    paintMargins();
  }
  /// end color theme support
  
  
  //Draw error bars
  auto proxyChartModel = std::dynamic_pointer_cast<Chart::WStandardChartProxyModel>( model() );
  assert( proxyChartModel );
  if( !proxyChartModel )
    return;

  std::shared_ptr<WAbstractItemModel> itemModelPtr = proxyChartModel->sourceModel();
  assert( itemModelPtr );
  if( !itemModelPtr )
    return;

  WAbstractItemModel *itemModel = itemModelPtr.get();
  const int nrows = proxyChartModel->rowCount();

  for( int row = sm_num_eqn_energy_rows; row < nrows; ++row )
  {
    const double energy = proxyChartModel->data(row,sm_energy_col);
    const double eff = proxyChartModel->data(row,sm_data_eff_col);
    const double uncert = Wt::asNumber( itemModel->data(row,sm_data_eff_col,Wt::ItemDataRole::User) );
      
    if( IsNan(energy) || IsNan(eff) || IsNan(uncert) )
      continue;
  
    const WPointF upper_uncert = mapToDevice( energy, eff + uncert );
    const WPointF lower_uncert = mapToDevice( energy, eff - uncert );
    painter.setPen( m_textPen ); // TODO: should we customize this?
    painter.drawLine( upper_uncert, lower_uncert );
  }
  
  //Paint some grey-ish rectangles over the energy we dont have peaks for
  if( !m_efficiencyCoefs.empty() && !m_datapoints.empty() )
  {
    float mindata = 9999999.9f, maxdata = -9999999.9f;
    for( const auto &d : m_datapoints )
    {
      mindata = std::min( mindata, d.energy );
      maxdata = std::max( maxdata, d.energy );
    }
  
    const WPointF left_ll = mapToDevice( m_det_lower_energy, axis(Chart::Axis::Y).minimum() );
    const WPointF left_ur = mapToDevice( mindata, axis(Chart::Axis::Y).maximum() );

    const WPointF right_ll = mapToDevice( maxdata, axis(Chart::Axis::Y).minimum() );
    const WPointF right_ur = mapToDevice( m_det_upper_energy, axis(Chart::Axis::Y).maximum() );

    painter.setBrush( WBrush( WColor(123,123,123,25) ) ); //TODO: get this from the color theme
    painter.setPen( WPen(Wt::PenStyle::None) );
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


void MakeDrfChart::setFwhmCoefficients( const std::vector<float> &coeffs,
                                        const std::vector<float> &uncerts,
                                        const FwhmCoefType eqnType,
                                        const EqnEnergyUnits units )
{
  m_fwhmCoefs = coeffs;
  m_fwhmCoefUncerts = uncerts;
  m_fwhmEqnType = eqnType;
  m_fwhmEnergyUnits = units;
  updateFwhmEquationToModel();
  updateYAxisRange();
}//void setFwhmCoefficients( const std::vector<double> &coeffs )


void MakeDrfChart::setEfficiencyCoefficients( const std::vector<float> &coeffs,
                                              const std::vector<float> &uncerts,
                                              const EqnEnergyUnits units )
{
  m_efficiencyCoefs = coeffs;
  m_efficiencyCoefUncerts = uncerts;
  m_efficiencyEnergyUnits = units;
  updateEffEquationToModel();
  updateYAxisRange();
}//void setEfficiencyCoefficients( const std::vector<double> &coeffs )


void MakeDrfChart::setDataPoints( const std::vector<MakeDrfChart::DataPoint> &datapoints,
                                  const float det_diameter,
                                  const float det_setback,
                                  const float lower_energy, const float upper_energy )
{
  m_datapoints = datapoints;
  m_det_diameter = det_diameter;
  m_det_setback = det_setback;
  
  const bool update_xrange = ((lower_energy != m_det_lower_energy) || (upper_energy != m_det_upper_energy));
  
  if( update_xrange )
  {
    m_det_lower_energy = lower_energy;
    m_det_upper_energy = upper_energy;
    
    updateEqnEnergyToModel();
    updateEffEquationToModel();
    updateFwhmEquationToModel();
    
    axis(Chart::Axis::X).setRange( lower_energy, upper_energy );
  }//if( energy range changed )

  updateDataToModel();

  if( det_diameter > 0.0 )
    axis(Chart::Axis::Y).setTitle( WString::tr("md-geom-intrinsic-eff") );
  else
    axis(Chart::Axis::Y).setTitle( WString::tr("Efficiency") );
  
  if( update_xrange )
    m_xRangeChanged.emit(m_det_lower_energy,m_det_upper_energy);
  
  updateYAxisRange();
}//void setDataPoints(...)


const std::vector<MakeDrfChart::DataPoint> &MakeDrfChart::currentDataPoints() const
{
  return m_datapoints;
}


float MakeDrfChart::currentDiameter() const
{
  return m_det_diameter;
}

void MakeDrfChart::showFwhmPoints( const bool show )
{
  series(sm_data_fwhm_col).setHidden( !show );
  series(sm_equation_fwhm_col).setHidden( !show );
  axis(Chart::Axis::Y2).setVisible( show );
  setPlotAreaPadding( (show ? 55 : 10), Wt::Side::Right );
}//void showFwhmPoints( const bool show )


void MakeDrfChart::showEfficiencyPoints( const bool show )
{
  series(sm_data_eff_col).setHidden( !show );
  series(sm_equation_eff_col).setHidden( !show );
  axis(Chart::Axis::Y1).setVisible( show );
  setPlotAreaPadding( (show ? 55 : 10), Wt::Side::Left );
}//void showEfficiencyPoints( const bool show )


void MakeDrfChart::setXRange( double lower, double upper )
{
  if( upper < lower )
    std::swap( lower, upper );
  
  if( fabs(lower-upper) < 1.0 || IsNan(lower) || IsNan(upper) || IsInf(lower) || IsInf(upper) )
    return;
  
  axis(Chart::Axis::X).setRange( lower, upper );
  updateYAxisRange();
}//void setXRange( double lower, upper )


Wt::Signal<double,double> &MakeDrfChart::xRangeChanged()
{
  return m_xRangeChanged;
}
