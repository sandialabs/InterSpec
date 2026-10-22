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

#include <map>
#include <set>
#include <cmath>
#include <mutex>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <Wt/WFont.h>
#include <Wt/WPainter.h>
#include <Wt/WSvgImage.h>
#include <Wt/WPaintDevice.h>
#include <Wt/Chart/WDataSeries.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/SpectrumChartImage.h"

using namespace std;
using namespace Wt;


namespace
{
  const float ns_large_render_font_size = 14;
  const float ns_compact_render_font_size = 10;

  /** Images narrower than this get compact axes (no tick labels on the y-axis). */
  const float ns_min_large_render_width = 480;

  std::mutex ns_raster_renderer_mutex;
  SpectrumChartImage::RasterRenderer ns_raster_renderer;


  /** A `SpectrumChart` that paints peaks given to it directly, so it works without a `PeakModel`
   (which needs an InterSpec session), and at any size.
   */
  class StandaloneSpectrumChart : public SpectrumChart
  {
  public:
    StandaloneSpectrumChart()
      : SpectrumChart(),
        m_directPeaks( nullptr )
    {
    }

    void setDirectPeaks( std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks )
    {
      m_directPeaks = peaks;
    }

    // The base class refuses to render below 50 px, which is too restrictive for thumbnails.
    bool isLargeEnough( Wt::WPainter & ) const override { return true; }

    void paintPeaks( Wt::WPainter &painter ) const override
    {
      if( !m_directPeaks || m_directPeaks->empty() )
        return;

      const double minx = axis( Chart::Axis::X ).minimum();
      const double maxx = axis( Chart::Axis::X ).maximum();

      std::vector<std::shared_ptr<const PeakDef>> nongauspeaks;

      // Gaussian peaks are grouped by the continuum (ROI) they share.
      std::map<const PeakContinuum *, std::vector<std::shared_ptr<const PeakDef>>> continuum_groups;
      for( const std::shared_ptr<const PeakDef> &peak : *m_directPeaks )
      {
        if( !peak || (peak->lowerX() > maxx) || (peak->upperX() < minx) )
          continue;

        if( !peak->gausPeak() )
          nongauspeaks.push_back( peak );
        else if( peak->continuum() )
          continuum_groups[peak->continuum().get()].push_back( peak );
      }

      for( auto &entry : continuum_groups )
      {
        std::vector<std::shared_ptr<const PeakDef>> &group = entry.second;
        std::sort( begin(group), end(group), &PeakDef::lessThanByMeanShrdPtr );
        paintGausPeaks( group, painter );
      }

      for( const std::shared_ptr<const PeakDef> &peak : nongauspeaks )
        paintNonGausPeak( *peak, painter );
    }//paintPeaks(...)

  private:
    std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> m_directPeaks;
  };//class StandaloneSpectrumChart


  /** Sets the chart's data series from the model, after its colors are set. */
  void set_series( SpectrumChart &chart, SpectrumDataModel &model, const double pen_width = -1.0 )
  {
    vector<Chart::WDataSeries> series = model.suggestDataSeries();
    std::vector<std::unique_ptr<Chart::WDataSeries>> series_ptrs;
    series_ptrs.reserve( series.size() );
    for( Chart::WDataSeries &s : series )
    {
      if( pen_width > 0.0 )
      {
        WPen pen = s.pen();
        pen.setWidth( WLength( pen_width ) );
        s.setPen( pen );
      }
      series_ptrs.push_back( std::make_unique<Chart::WDataSeries>( std::move(s) ) );
    }
    chart.setSeries( std::move(series_ptrs) );
  }//set_series(...)


  void apply_theme( SpectrumChart &chart, SpectrumDataModel &model, const ColorTheme *theme )
  {
    if( !theme )
      return;

    model.setForegroundSpectrumColor( theme->foregroundLine );
    model.setBackgroundSpectrumColor( theme->backgroundLine );
    model.setSecondarySpectrumColor( theme->secondaryLine );
    chart.setDefaultPeakColor( theme->defaultPeakLine );
    chart.setAxisLineColor( theme->spectrumAxisLines );
    chart.setChartMarginColor( theme->spectrumChartMargins );
    chart.setChartBackgroundColor( theme->spectrumChartBackground );
    chart.setTextColor( theme->spectrumChartText );
  }//apply_theme(...)


  /** Font sizes and padding for an image `width` pixels wide. */
  void set_axis_layout( SpectrumChart &chart, const float width, const std::string &x_title,
                        const std::string &compact_x_title, const bool has_title, const bool y2 )
  {
    const bool large = (width >= ns_min_large_render_width);

    WFont label_font = chart.axis( Chart::Axis::Y1 ).labelFont();
    label_font.setSize( WLength( large ? ns_large_render_font_size : ns_compact_render_font_size ) );

    chart.setPlotAreaPadding( 0, Wt::Side::Top );
    if( large )
    {
      chart.setPlotAreaPadding( 80, Wt::Side::Left );
      chart.setPlotAreaPadding( 42, Wt::Side::Bottom );
      chart.setPlotAreaPadding( y2 ? 20 : 10, Wt::Side::Right );
      if( has_title )
        chart.setPlotAreaPadding( ns_large_render_font_size + 4, Wt::Side::Top );
      chart.axis( Chart::Axis::X ).setTitle( x_title );
      chart.setCompactAxis( false );
    }else
    {
      chart.setPlotAreaPadding( 5, Wt::Side::Left );
      chart.setPlotAreaPadding( 20, Wt::Side::Bottom );
      chart.setPlotAreaPadding( 5, Wt::Side::Right );
      chart.axis( Chart::Axis::X ).setTitle( compact_x_title );
      chart.setCompactAxis( true );
    }

    // Room for labels of four or five digits.
    chart.setXAxisMinimumLabelSpacing( 5.5 * (large ? ns_large_render_font_size : ns_compact_render_font_size) );

    for( const Chart::Axis axis : { Chart::Axis::X, Chart::Axis::Y1, Chart::Axis::Y2 } )
    {
      chart.axis( axis ).setLabelFont( label_font );
      chart.axis( axis ).setTitleFont( label_font );
    }

    WFont title_font = chart.titleFont();
    title_font.setSize( label_font.sizeLength() );
    chart.setTitleFont( title_font );
  }//set_axis_layout(...)
}//namespace


namespace SpectrumChartImage
{

bool paintSpectrum( Wt::WPaintDevice &device, const SpectrumOptions &options )
{
  const std::shared_ptr<const SpecUtils::Measurement> &meas = options.foreground;
  if( !meas || (meas->num_gamma_channels() < 2) )
    return false;

  const float width = static_cast<float>( device.width().toPixels() );
  const bool large = (width >= ns_min_large_render_width);

  auto model = std::make_shared<SpectrumDataModel>();
  StandaloneSpectrumChart chart;
  chart.setModel( model );

  // Peaks are drawn on their fitted continuum, which is not background-subtracted (the D3 chart
  //  subtracts it; this chart can not), so they are left off a background-subtracted chart.
  const bool subtracting = (options.background && options.background_subtract);
  chart.setDirectPeaks( subtracting ? nullptr : options.peaks );

  model->setDataHistogram( meas );
  if( options.background )
  {
    model->setBackgroundHistogram( options.background );
    if( options.background_scale > 0.0 )
      model->setBackgroundDataScaleFactor( static_cast<float>( options.background_scale ) );
    model->setBackgroundSubtract( options.background_subtract );
  }

  if( options.secondary )
  {
    model->setSecondDataHistogram( options.secondary, false );
    if( options.secondary_scale > 0.0 )
      model->setSecondDataScaleFactor( static_cast<float>( options.secondary_scale ) );
  }

  apply_theme( chart, *model, options.theme.get() );
  set_series( chart, *model, options.thumbnail ? 0.5 : -1.0 );

  // A background-subtracted chart shows just the difference, as the GUI does.
  if( options.background && options.background_subtract && (model->backgroundColumn() >= 0) )
    chart.series( model->backgroundColumn() ).setHidden( true );

  const bool has_title = (!options.thumbnail && large && !options.title.empty());

  if( options.thumbnail )
  {
    // Axis lines, but no tick labels or titles.  (setCompactAxis(true) subtracts 17 px from the
    //  bottom padding, so the padding is set after it.)
    chart.setCompactAxis( true );
    chart.setPlotAreaPadding( 4, Wt::Side::Left );
    chart.setPlotAreaPadding( 4, Wt::Side::Bottom );
    chart.setPlotAreaPadding( 2, Wt::Side::Right );
    chart.setPlotAreaPadding( 2, Wt::Side::Top );
    chart.axis( Chart::Axis::X ).setLabelFormat( " " );
    chart.axis( Chart::Axis::Y1 ).setLabelFormat( " " );
  }else
  {
    set_axis_layout( chart, width, "Energy (keV)", "(keV)", has_title, false );
    if( has_title )
      chart.setTitle( WString::fromUTF8( options.title ) );
  }//if( thumbnail ) / else

  // Reference lines; the first is the "current" one, the others persisted.
  for( size_t i = options.reference_lines.size(); i > 1; --i )
  {
    if( options.reference_lines[i-1] )
    {
      chart.setReferncePhotoPeakLines( *options.reference_lines[i-1] );
      chart.persistCurrentReferncePhotoPeakLines();
    }
  }
  if( !options.reference_lines.empty() && options.reference_lines[0] )
    chart.setReferncePhotoPeakLines( *options.reference_lines[0] );

  const size_t nchannel = meas->num_gamma_channels();
  double lowx = meas->gamma_channel_lower( 0 );
  double upperx = meas->gamma_channel_upper( nchannel - 1 );
  if( options.energy_range && (options.energy_range->second > options.energy_range->first) )
  {
    lowx = options.energy_range->first;
    upperx = options.energy_range->second;
  }
  chart.setXAxisRange( lowx, upperx );


  // Combine channels so there are about two per pixel.
  const size_t displayed_nbin = meas->find_gamma_channel( static_cast<float>( upperx ) )
                                - meas->find_gamma_channel( static_cast<float>( lowx ) );
  const int plot_width = std::max( 1, static_cast<int>( width ) - chart.plotAreaPadding( Side::Left )
                                      - chart.plotAreaPadding( Side::Right ) );
  const float bins_per_pixel = float( displayed_nbin ) / float( plot_width );
  const int factor = std::max( static_cast<int>( std::ceil( bins_per_pixel / 2 ) ), 1 );
  model->setRebinFactor( factor );

  chart.setYAxisScale( options.log_y ? Wt::Chart::AxisScale::Log : Wt::Chart::AxisScale::Linear );
  chart.setAutoYAxisRange();

  if( !options.thumbnail && large )
    chart.axis( Chart::Axis::Y1 ).setTitle( (factor == 1) ? std::string( "Counts Per Channel" )
                                            : ("Counts Per " + std::to_string( factor ) + " Channels") );

  WPainter p( &device );
  chart.paint( p );

  // The logo, in the upper right, for thumbnails.
  if( options.thumbnail && !options.logo_path.empty() )
  {
    WPainter::Image logo( options.logo_path, options.logo_path );
    if( (logo.width() > 0) && (logo.height() > 0) )
    {
      const double img_w = device.width().toPixels();
      const double img_h = device.height().toPixels();
      const double logo_h = 0.35 * img_h;
      const double logo_w = logo_h * (double( logo.width() ) / double( logo.height() ));
      p.drawImage( WRectF( img_w - logo_w - 4, 4, logo_w, logo_h ), logo );
    }
  }

  p.end();

  return true;
}//paintSpectrum(...)


bool paintTimeSeries( Wt::WPaintDevice &device, const TimeSeriesOptions &options )
{
  const std::shared_ptr<const SpecMeas> &meas = options.meas;
  if( !meas )
    return false;

  const float width = static_cast<float>( device.width().toPixels() );

  const std::set<int> &sample_numbers = meas->sample_numbers();
  const std::vector<int> samples( begin(sample_numbers), end(sample_numbers) );
  const size_t num_samples = samples.size();
  if( num_samples < 2 )
    return false;

  const std::vector<std::string> &detectors = options.detectors.empty() ? meas->detector_names()
                                                                         : options.detectors;

  auto gamma_cps = std::make_shared<vector<float>>( num_samples, 0.0f );
  auto neutron_cps = std::make_shared<vector<float>>( num_samples, 0.0f );
  auto times = std::make_shared<vector<float>>( num_samples + 1, 0.0f );
  vector<bool> occupied( num_samples, false );

  bool has_neutrons = false;
  float time_sum = 0.0f;

  for( size_t i = 0; i < num_samples; ++i )
  {
    (*times)[i] = time_sum;

    int ndet = 0;
    bool is_occupied = false;
    for( const std::string &det : detectors )
    {
      const std::shared_ptr<const SpecUtils::Measurement> m = meas->measurement( samples[i], det );
      if( m && (m->real_time() > 0.0f) )
        ++ndet;
      if( m )
        is_occupied = (is_occupied || (m->occupied() == SpecUtils::OccupancyStatus::Occupied));
    }
    ndet = std::max( ndet, 1 );

    const std::shared_ptr<SpecUtils::Measurement> sum = meas->sum_measurements( { samples[i] }, detectors, nullptr );
    if( !sum )
      continue;

    float real_time = sum->real_time() / ndet;
    real_time = (real_time > 0.0f) ? real_time : 0.1f;

    has_neutrons = (has_neutrons || sum->contained_neutron());
    (*gamma_cps)[i] = static_cast<float>( sum->gamma_count_sum() / real_time );
    (*neutron_cps)[i] = static_cast<float>( sum->neutron_counts_sum() / real_time );
    occupied[i] = is_occupied;
    time_sum += real_time;
  }//for( loop over samples )
  (*times)[num_samples] = time_sum;

  // The "spectra" are count rates, with time as the x-axis.
  auto gamma_hist = std::make_shared<SpecUtils::Measurement>();
  auto neutron_hist = std::make_shared<SpecUtils::Measurement>();
  gamma_hist->set_title( "Gammas" );
  neutron_hist->set_title( "Neutrons" );
  gamma_hist->set_gamma_counts( gamma_cps, 0.0f, 0.0f );
  neutron_hist->set_gamma_counts( neutron_cps, 0.0f, 0.0f );

  auto time_cal = std::make_shared<SpecUtils::EnergyCalibration>();
  time_cal->set_lower_channel_energy( num_samples, *times );
  gamma_hist->set_energy_calibration( time_cal );
  neutron_hist->set_energy_calibration( time_cal );

  auto model = std::make_shared<SpectrumDataModel>();
  StandaloneSpectrumChart chart;
  chart.setModel( model );

  set_axis_layout( chart, width, "Time (s)", "(sec.)", false, has_neutrons );

  model->setDataHistogram( gamma_hist );
  if( has_neutrons )
  {
    model->setSecondDataHistogram( neutron_hist, true );
    chart.axis( Chart::Axis::Y2 ).setVisible( true );
  }

  // The highlighted time ranges: runs of consecutive samples in each set.
  const auto ranges_for = [&]( const std::function<bool(size_t)> &included ){
    vector<pair<double,double>> ranges;
    for( size_t i = 0; i < num_samples; ++i )
    {
      if( !included( i ) )
        continue;
      if( !ranges.empty() && (ranges.back().second == (*times)[i]) )
        ranges.back().second = (*times)[i+1];
      else
        ranges.emplace_back( (*times)[i], (*times)[i+1] );
    }
    return ranges;
  };//ranges_for

  if( options.highlight_occupied )
  {
    chart.setTimeHighLightRegions( ranges_for( [&]( size_t i ){ return occupied[i]; } ),
                                   SpectrumChart::ForegroundHighlight );
  }else
  {
    const std::pair<const std::set<int> *, SpectrumChart::HighlightRegionType> highlights[] = {
      { &options.foreground_samples, SpectrumChart::ForegroundHighlight },
      { &options.background_samples, SpectrumChart::BackgroundHighlight },
      { &options.secondary_samples, SpectrumChart::SecondForegroundHighlight }
    };
    for( const auto &[sample_set, type] : highlights )
    {
      if( !sample_set->empty() && (sample_set->size() < num_samples) )
        chart.setTimeHighLightRegions( ranges_for( [&]( size_t i ){ return sample_set->count( samples[i] ) > 0; } ),
                                       type );
    }
  }//if( highlight_occupied ) / else

  apply_theme( chart, *model, options.theme.get() );
  set_series( chart, *model );

  chart.setXAxisRange( times->front(), times->back() );
  chart.axis( Chart::Axis::Y1 ).setAutoLimits( Chart::AxisValue::Minimum | Chart::AxisValue::Maximum );
  if( has_neutrons )
    chart.axis( Chart::Axis::Y2 ).setAutoLimits( Chart::AxisValue::Minimum | Chart::AxisValue::Maximum );

  WPainter p( &device );
  chart.paint( p );
  p.end();

  return true;
}//paintTimeSeries(...)


void setRasterRenderer( RasterRenderer renderer )
{
  std::lock_guard<std::mutex> lock( ns_raster_renderer_mutex );
  ns_raster_renderer = std::move( renderer );
}


bool canRenderRaster()
{
  std::lock_guard<std::mutex> lock( ns_raster_renderer_mutex );
  return static_cast<bool>( ns_raster_renderer );
}


std::string renderImage( const int width, const int height, const std::string &format,
                         const std::function<void( Wt::WPaintDevice & )> &paint )
{
  if( (width <= 0) || (height <= 0) )
    throw std::runtime_error( "Invalid image size." );

  if( format == "svg" )
  {
    WSvgImage image( width, height );
    paint( image );
    std::ostringstream out;
    image.write( out );
    return out.str();
  }//if( format == "svg" )

  RasterRenderer renderer;
  {
    std::lock_guard<std::mutex> lock( ns_raster_renderer_mutex );
    renderer = ns_raster_renderer;
  }

  if( !renderer )
    throw std::runtime_error( "Rendering " + format + " images is not available in this build." );

  return renderer( width, height, format, paint );
}//renderImage(...)

}//namespace SpectrumChartImage
