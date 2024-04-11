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
#if( WT_VERSION >= 0x3030600 )
#include <Wt/Chart/WAbstractChartModel>
#endif


#include <Wt/Chart/WCartesianChart>
#if( WT_VERSION >= 0x3030600 )
#include <Wt/Chart/WStandardChartProxyModel>
#endif

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


DrfChartHolder::DrfChartHolder( MakeDrfChart *chart, WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_chart( chart )
{
  setLayoutSizeAware( true );
  addWidget( m_chart );
}
  
DrfChartHolder::~DrfChartHolder()
{
};
  
void DrfChartHolder::layoutSizeChanged( int width, int height )
{
  m_chart->resize( width, height );
}




MakeDrfChart::MakeDrfChart( Wt::WContainerWidget *parent )
: Wt::Chart::WCartesianChart( parent ),
  m_det_diameter( 1.0*PhysicalUnits::cm ),
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
  m_textPen( WColor(GlobalColor::black) )
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
    
#if( WT_VERSION < 0x3030600 )
  setModel( m );
#else
  // TODO: make a class that derives from either Wt::WAbstractItemModel or
  //       Wt::Chart::WAbstractChartModel class, depending on Wt version,
  //       and implement all functionality needed in it
  auto *proxyModel = new Chart::WStandardChartProxyModel( m, this );
  setModel( proxyModel );
#endif
    
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
  
  //For the moment we wont show the errors since I havent verified they are actually reasonable and correct
  /*
  Chart::WDataSeries eqn_eff_neg_uncert_series( sm_equation_eff_neg_uncert_col, Chart::SeriesType::CurveSeries, Chart::YAxis );
  eqn_eff_neg_uncert_series.setMarker( Wt::Chart::MarkerType::NoMarker );
  eqn_eff_neg_uncert_series.setPen( WPen(Wt::green) );
  addSeries( eqn_eff_neg_uncert_series );
  
  Chart::WDataSeries eqn_eff_pos_uncert_series( sm_equation_eff_pos_uncert_col, Chart::SeriesType::CurveSeries, Chart::YAxis );
  eqn_eff_pos_uncert_series.setMarker( Wt::Chart::MarkerType::NoMarker );
  eqn_eff_pos_uncert_series.setPen( WPen(Wt::red) );
  addSeries( eqn_eff_pos_uncert_series );
  */
  
  
  Chart::WDataSeries data_fwhm_series( sm_data_fwhm_col, Chart::SeriesType::PointSeries, Chart::Y2Axis );
  data_fwhm_series.setMarker( Chart::MarkerType::XCrossMarker );
  data_fwhm_series.setMarkerSize( 6 );
  addSeries( data_fwhm_series );
  
  Chart::WDataSeries eqn_fwhm_series( sm_equation_fwhm_col, Chart::SeriesType::CurveSeries, Chart::Y2Axis );
  eqn_fwhm_series.setMarker( Wt::Chart::MarkerType::NoMarker );
  addSeries( eqn_fwhm_series );
  
  setPlotAreaPadding(0, Wt::Top);
  setPlotAreaPadding(55, Wt::Right | Wt::Left);
  //axis(Chart::XAxis).setTitle( "Energy (keV)" );
  setPlotAreaPadding(25, Wt::Bottom);
  
  axis(Chart::YAxis).setVisible( true );
  axis(Chart::YAxis).setTitle( "Intrinsic Eff." );
  
  axis(Chart::Y2Axis).setVisible( true );
  axis(Chart::Y2Axis).setTitle( "FWHM" );

#if( WT_VERSION >= 0x3030400 )
  axis(Wt::Chart::Y1Axis).setTitleOrientation( Wt::Vertical );
  axis(Wt::Chart::Y2Axis).setTitleOrientation( Wt::Vertical );
#endif
  
  setAxisPadding( 0 );
  axis(Chart::Y1Axis).setMargin( 0 );
  axis(Chart::Y2Axis).setMargin( 0 );
  
  m->setHeaderData( sm_energy_col, Wt::Horizontal, boost::any(WString("Energy (keV)")), Wt::DisplayRole );
  m->setHeaderData( sm_data_eff_col, Wt::Horizontal, boost::any(WString("Data Intrinsic Eff.")), Wt::DisplayRole );
  m->setHeaderData( sm_data_fwhm_col, Wt::Horizontal, boost::any(WString("Data FWHM")), Wt::DisplayRole );
  m->setHeaderData( sm_equation_eff_col, Wt::Horizontal, boost::any(WString("Fit Intrinsic Eff.")), Wt::DisplayRole );
  m->setHeaderData( sm_equation_eff_neg_uncert_col, Wt::Horizontal, boost::any(WString("Intrinsic Eff. +2&#963;")), Wt::DisplayRole );
  m->setHeaderData( sm_equation_eff_pos_uncert_col, Wt::Horizontal, boost::any(WString("Fit Intrinsic Eff. -2&#963;")), Wt::DisplayRole );
  m->setHeaderData( sm_equation_fwhm_col, Wt::Horizontal, boost::any(WString("Fit FWHM")), Wt::DisplayRole );
  
  setLegendEnabled( true );
  setLegendLocation( Wt::Chart::LegendLocation::LegendInside, Wt::Top, Wt::AlignmentFlag::AlignRight );
  setLegendColumns( 1, WLength(125,WLength::Pixel) );
  
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
  const double xmin = axis(Chart::XAxis).minimum();
  const double xmax = axis(Chart::XAxis).maximum();
  
  double miny = 9999.9f, maxy = -9999.9f, miny2 = 9999.9f, maxy2 = -9999.9f;
  
#if( WT_VERSION < 0x3030600 )
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
#else
  auto *proxyModel = dynamic_cast<Chart::WStandardChartProxyModel *>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel() : nullptr;
#endif
  
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
    axis(Chart::YAxis).setRange(0.0, 1.0);
  }else
  {
    miny = (miny < 0.1) ? 0.0 : std::floor(9*miny)/10.0;
    maxy = ((11.0*maxy) > 1.5) ? (std::ceil(11.0*maxy)/10.0) : 1.15*maxy;
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
#if( WT_VERSION < 0x3030600 )
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
#else
  auto *proxyModel = dynamic_cast<Chart::WStandardChartProxyModel *>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel() : nullptr;
#endif
  
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
    
    m->setData( row, sm_energy_col, energy, Wt::DisplayRole );
    
    if( (data.peak_area > 0.0f) && (data.source_count_rate > 0.0f)
        && (data.livetime > 0.0f) && ((m_det_diameter > 0.0f) || (data.distance < 0.0)) )
    {
      const double fracSolidAngle
           = (data.distance < 0.0) ? 1.0
                                   : DetectorPeakResponse::fractionalSolidAngle( m_det_diameter, data.distance );
      const double expected = data.source_count_rate * data.livetime * fracSolidAngle;
      const double eff = data.peak_area / expected;
      
      double fracUncert2 = 0.0;
      if( data.peak_area_uncertainty > 0.0f )
        fracUncert2 += std::pow( data.peak_area_uncertainty / data.peak_area, 2.0f );
      if( data.source_count_rate_uncertainty > 0.0f )
        fracUncert2 += std::pow( data.source_count_rate_uncertainty / data.source_count_rate, 2.0f );
      const double eff_uncert = eff * std::sqrt(fracUncert2);
      
      m->setData( row, sm_data_eff_col, boost::any(eff), Wt::DisplayRole );
      if( eff_uncert > 0.0 )
        m->setData( row, sm_data_eff_col, boost::any(eff_uncert), Wt::UserRole );
      m->setData( row, sm_data_eff_col, boost::any(data.peak_color), Wt::MarkerPenColorRole );
      m->setData( row, sm_data_eff_col, boost::any(data.peak_color), Wt::MarkerBrushColorRole );
      if( data.source_information.size() )
        m->setData( row, sm_data_eff_col, boost::any( WString::fromUTF8(data.source_information)), Wt::ToolTipRole );
    }//if( we can calculate a efficiency )
    
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
#if( WT_VERSION < 0x3030600 )
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
#else
  auto *proxyModel = dynamic_cast<Chart::WStandardChartProxyModel *>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel() : nullptr;
#endif
  
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
#if( WT_VERSION < 0x3030600 )
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
#else
  auto *proxyModel = dynamic_cast<Chart::WStandardChartProxyModel *>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel() : nullptr;
#endif
  
  assert( m );
  if( !m )
    return;
  
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  
  if( m_efficiencyCoefs.empty() )
  {
    if( m->data(0, sm_equation_eff_col).empty() )
      return;
    
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
    {
      m->setData( row, sm_equation_eff_col, boost::any() );
      m->setData( row, sm_equation_eff_pos_uncert_col, boost::any() );
      m->setData( row, sm_equation_eff_neg_uncert_col, boost::any() );
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
      m->setData( row, sm_equation_eff_col, boost::any() );
      m->setData( row, sm_equation_eff_pos_uncert_col, boost::any() );
      m->setData( row, sm_equation_eff_neg_uncert_col, boost::any() );
    }else
    {
      m->setData( row, sm_equation_eff_col, boost::any(eff) );
      
      boost::any lowerval, upperval;
      
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
          upperval = boost::any( eff + posuncert );
        if( IsNan(neguncert) || IsInf(neguncert) )
          lowerval = boost::any( eff - neguncert );
      }//if( we have uncertainties )
      
      m->setData( row, sm_equation_eff_pos_uncert_col, upperval );
      m->setData( row, sm_equation_eff_neg_uncert_col, lowerval );
    }//if( valid eff value ) / else
  }//for( loop over eqn rows )
}//void updateEffEquationToModel()


void MakeDrfChart::updateFwhmEquationToModel()
{
#if( WT_VERSION < 0x3030600 )
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
#else
  auto *proxyModel = dynamic_cast<Chart::WStandardChartProxyModel *>( model() );
  WAbstractItemModel *m = proxyModel ? proxyModel->sourceModel() : nullptr;
#endif

  assert( m );
  if( !m )
    return;
  
  assert( m->rowCount() >= sm_num_eqn_energy_rows );
  
  if( ((m_fwhmEqnType==FwhmCoefType::Gadras) && m_fwhmCoefs.size() != 3)
     || ((m_fwhmEqnType==FwhmCoefType::SqrtEnergyPlusInverse) && m_fwhmCoefs.size() != 3)
     || ((m_fwhmEqnType==FwhmCoefType::ConstantPlusSqrtEnergy) && m_fwhmCoefs.size() != 2)
     || (m_fwhmEqnType==FwhmCoefType::SqrtEqn && m_fwhmCoefs.size() < 1) )
  {
    if( m->data(0, sm_equation_fwhm_col).empty() )  //
      return;
    
    for( int row = 0; row < sm_num_eqn_energy_rows; ++row )
      m->setData( row, sm_equation_fwhm_col, boost::any() );
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
    m->setData( row, sm_equation_fwhm_col, boost::any(fwhm) );
  }//for( loop over eqn rows )
}//void updateFwhmEquationToModel()


void MakeDrfChart::updateColorTheme( std::shared_ptr<const ColorTheme> theme )
{
  if( !theme )
    return; //Should reset colors to default, but whatever for now
  
  
  if( theme->spectrumAxisLines.isDefault() )
    m_textPen = WPen( WColor(GlobalColor::black) );
  else
    m_textPen = WPen( theme->spectrumAxisLines );
  setTextPen( m_textPen );
  
  if( theme->spectrumAxisLines.isDefault() )
  {
    axis(Chart::XAxis).setPen( WPen(GlobalColor::black) );
    axis(Chart::YAxis).setPen( WPen(GlobalColor::black) );
    axis(Chart::Y2Axis).setPen( WPen(GlobalColor::black) );
  }else
  {
    axis(Chart::XAxis).setPen( WPen(theme->spectrumAxisLines) );
    axis(Chart::YAxis).setPen( WPen(theme->spectrumAxisLines) );
    axis(Chart::Y2Axis).setPen( WPen(theme->spectrumAxisLines) );
  }
  
  if( theme->spectrumChartMargins.isDefault() )
    m_chartMarginBrush = WBrush();
  else
    m_chartMarginBrush = WBrush(theme->spectrumChartMargins);
  
  if( theme->spectrumChartBackground.isDefault() )
    setBackground( Wt::NoBrush );
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
    effcolor = WColor(GlobalColor::black);
  setSeriesColor( eqnEffSeries, effcolor );
  
  Wt::Chart::WDataSeries &eqnFwhmSeries = series(sm_equation_fwhm_col);
  WColor fwhmcolor = theme->backgroundLine;
  if( fwhmcolor.isDefault() )
    fwhmcolor = WColor(GlobalColor::gray);
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
    
    const int padLeft = plotAreaPadding(Left);
    const int padRight = plotAreaPadding(Right);
    const int padTop = plotAreaPadding(Top);
    const int padBottom = plotAreaPadding(Bottom);
    
    WRectF area;
    if( orientation() == Vertical )
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
  
  if( background().style() != NoBrush && m_chartMarginBrush.style() != NoBrush )
  {
    paintMargins();
    painter.fillRect( hv(plotArea()), background() );
  }else if( background().style() != NoBrush )
  {
    painter.fillRect( hv(rectangle), background() );
  }else if( m_chartMarginBrush.style() != NoBrush )
  {
    paintMargins();
  }
  /// end color theme support
  
  
  //Draw error bars
#if( WT_VERSION < 0x3030600 )
  WStandardItemModel *m = static_cast<WStandardItemModel *>( model() );
#else
  auto *m = dynamic_cast<Chart::WStandardChartProxyModel *>( model() );
  assert( m );
  if( !m )
    return;
    
  WAbstractItemModel *itemModel = m->sourceModel();
  assert( itemModel );
  if( !itemModel )
    return;
#endif
  
  const int nrows = m->rowCount();
  
  for( int row = sm_num_eqn_energy_rows; row < nrows; ++row )
  {
#if( WT_VERSION < 0x3030600 )
    const double energy = Wt::asNumber( m->data(row,sm_energy_col,Wt::DisplayRole) );
    const double eff = Wt::asNumber( m->data(row,sm_data_eff_col,Wt::DisplayRole) );
    const double uncert = Wt::asNumber( m->data(row,sm_data_eff_col,Wt::UserRole) );
#else
    const double energy = m->data(row,sm_energy_col);
    const double eff = m->data(row,sm_data_eff_col);
    const double uncert = Wt::asNumber( itemModel->data(row,sm_data_eff_col,Wt::UserRole) );
#endif
      
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
  
    const WPointF left_ll = mapToDevice( m_det_lower_energy, axis(Chart::YAxis).minimum() );
    const WPointF left_ur = mapToDevice( mindata, axis(Chart::YAxis).maximum() );
    
    const WPointF right_ll = mapToDevice( maxdata, axis(Chart::YAxis).minimum() );
    const WPointF right_ur = mapToDevice( m_det_upper_energy, axis(Chart::YAxis).maximum() );
    
    painter.setBrush( WBrush( WColor(123,123,123,25) ) ); //TODO: get this from the color theme
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
  }//if( energy range changed )
  
  updateDataToModel();
  
  if( det_diameter > 0.0 )
    axis(Chart::YAxis).setTitle( "Intrinsic Eff." );
  else
    axis(Chart::YAxis).setTitle( "Efficiency" );
  
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
  axis(Chart::Y2Axis).setVisible( show );
  setPlotAreaPadding( (show ? 55 : 10), Wt::Right );
}//void showFwhmPoints( const bool show )


void MakeDrfChart::showEfficiencyPoints( const bool show )
{
  series(sm_data_eff_col).setHidden( !show );
  series(sm_equation_eff_col).setHidden( !show );
  axis(Chart::Y1Axis).setVisible( show );
  setPlotAreaPadding( (show ? 55 : 10), Wt::Left );
}//void showEfficiencyPoints( const bool show )


void MakeDrfChart::setXRange( double lower, double upper )
{
  if( upper < lower )
    std::swap( lower, upper );
  
  if( fabs(lower-upper) < 1.0 || IsNan(lower) || IsNan(upper) || IsInf(lower) || IsInf(upper) )
    return;
  
  axis(Chart::XAxis).setRange( lower, upper );
  updateYAxisRange();
}//void setXRange( double lower, upper )


Wt::Signal<double,double> &MakeDrfChart::xRangeChanged()
{
  return m_xRangeChanged;
}
