/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.

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

// Refactored from SpecPreviewCommon.cpp (originally created 12/25/17), to use real
// InterSpec classes instead of the forked QL* versions, for macOS Quick Look thumbnail rendering.

#include "SpecPreviewCommon.h"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <iostream>

#include <hpdf.h>

#include <Wt/WFont>
#include <Wt/WPainter>
#include <Wt/WPdfImage>
#include <Wt/WPainterPath>
#include <Wt/Chart/WDataSeries>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/SpectrumDataModel.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

using namespace std;
using namespace Wt;
using SpecUtils::Measurement;


namespace
{
  const float ns_large_render_font_size = 14;
  const float ns_compact_render_font_size = 10;

  enum SpecRenderType
  {
    CompactSpecRender,
    LargeSpecRender
  };


  void HPDF_STDCALL hpdf_error_handler( HPDF_STATUS error_no, HPDF_STATUS detail_no, void * /*user_data*/ )
  {
    char buf[200];
    snprintf( buf, sizeof( buf ), "PDF image error: error_no=%04X, detail_no=%d",
              (unsigned int)error_no, (int)detail_no );
    throw std::runtime_error( buf );
  }//hpdf_error_handler(...)


  /** A minimal subclass of SpectrumChart that can paint peaks without needing PeakModel.

   The real SpectrumChart::paintPeaks() reads peaks from a PeakModel, which requires an
   InterSpec application instance. This subclass stores peaks directly and overrides paintPeaks().
   */
  class PreviewSpectrumChart : public SpectrumChart
  {
  public:
    PreviewSpectrumChart()
      : SpectrumChart()
      , m_directPeaks( nullptr )
    {
    }

    void setDirectPeaks( std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks )
    {
      m_directPeaks = peaks;
    }

    // Always allow rendering, even at small sizes (e.g., Finder thumbnails).
    // The base class refuses to render below 50px, which is too restrictive.
    bool isLargeEnough( Wt::WPainter & ) const override { return true; }

    // Expose chartArea() for manual axis line drawing in thumbnails.
    const Wt::WRectF &getChartArea() const { return chartArea(); }

    void paintPeaks( Wt::WPainter &painter ) const override
    {
      if( !m_directPeaks || m_directPeaks->empty() )
        return;

      // Get the data model to check axis ranges
      const double minx = axis( Chart::XAxis ).minimum();
      const double maxx = axis( Chart::XAxis ).maximum();

      std::vector<std::shared_ptr<const PeakDef>> gauspeaks, nongauspeaks;

      for( const auto &peak : *m_directPeaks )
      {
        if( !peak || peak->lowerX() > maxx || peak->upperX() < minx )
          continue;

        if( peak->gausPeak() )
          gauspeaks.push_back( peak );
        else
          nongauspeaks.push_back( peak );
      }

      // Group gaussian peaks by shared continuum (same ROI)
      std::map<const PeakContinuum *, std::vector<std::shared_ptr<const PeakDef>>> continuum_groups;
      for( const auto &peak : gauspeaks )
        continuum_groups[peak->continuum().get()].push_back( peak );

      for( auto &entry : continuum_groups )
      {
        std::vector<std::shared_ptr<const PeakDef>> &group = entry.second;
        std::sort( group.begin(), group.end(), &PeakDef::lessThanByMeanShrdPtr );
        paintGausPeaks( group, painter );
      }

      for( const auto &peak : nongauspeaks )
        paintNonGausPeak( *peak, painter );
    }//paintPeaks(...)

  private:
    std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> m_directPeaks;
  };//class PreviewSpectrumChart


  /** Renders a spectrum to a WPaintDevice (PDF or SVG image).

   @param thumbnail If true, renders with axis lines but no tick labels or titles, thinner
          spectrum line, and a watermark logo.  Used for Finder thumbnail generation.
   @param logo_path Path to the InterSpec logo PNG for watermark (only used when thumbnail=true).
          Pass empty string to skip the watermark.
   */
  template <class ImageType>
  bool renderSpectrum( ImageType &image,
                       std::shared_ptr<const SpecUtils::Measurement> foreground,
                       std::shared_ptr<const SpecUtils::Measurement> background,
                       std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks,
                       const bool thumbnail,
                       const std::string &logo_path = "" )
  {
    if( !foreground )
      return false;

    const std::string title = foreground->title();
    const float width = image.width().toPixels();

    std::shared_ptr<const SpecUtils::Measurement> meas = foreground;
    std::shared_ptr<const SpecUtils::Measurement> back = background;

    SpectrumDataModel dataModel;
    PreviewSpectrumChart chart;
    chart.setModel( &dataModel );
    chart.setDirectPeaks( peaks );

    dataModel.setDataHistogram( meas );
    if( back )
    {
      dataModel.setBackgroundHistogram( back );
      if( meas->live_time() > 0.0f && back->live_time() > 0.0f )
        dataModel.setBackgroundDataScaleFactor( meas->live_time() / back->live_time() );
    }

    vector<Chart::WDataSeries> series = dataModel.suggestDataSeries();

    if( thumbnail )
    {
      // Use a thinner line for thumbnails
      for( Chart::WDataSeries &s : series )
      {
        WPen pen = s.pen();
        pen.setWidth( WLength( 0.5 ) );
        s.setPen( pen );
      }
    }

    chart.setSeries( series );

    chart.setPlotAreaPadding( 0, Wt::Top );

    if( thumbnail )
    {
      // Thumbnail: axis lines visible, but no tick labels or titles
      chart.setPlotAreaPadding( 4, Wt::Left );
      chart.setPlotAreaPadding( 4, Wt::Bottom );
      chart.setPlotAreaPadding( 2, Wt::Right );
      chart.setPlotAreaPadding( 2, Wt::Top );
      chart.setCompactAxis( true );

      // Hide tick labels but keep axis lines
      chart.axis( Chart::XAxis ).setLabelFormat( " " );
      chart.axis( Chart::YAxis ).setLabelFormat( " " );
    }else
    {
      const SpecRenderType render_type = ((width < 641) ? CompactSpecRender : LargeSpecRender);

      WFont labelFont = chart.axis( Chart::YAxis ).labelFont();
      WFont titleFont = chart.titleFont();

      switch( render_type )
      {
        case LargeSpecRender:
          labelFont.setSize( WLength( ns_large_render_font_size ) );
          chart.setPlotAreaPadding( 62, Wt::Left );
          chart.setPlotAreaPadding( 42, Wt::Bottom );
          chart.setPlotAreaPadding( 10, Wt::Right );
          if( title.size() )
          {
            chart.setPlotAreaPadding( ns_large_render_font_size + 4, Wt::Top );
            titleFont.setSize( ns_large_render_font_size );
          }
          chart.axis( Chart::XAxis ).setTitle( "Energy (keV)" );
          chart.setCompactAxis( false );
          if( title.size() )
            chart.setTitle( title );
          break;

        case CompactSpecRender:
          labelFont.setSize( WLength( ns_compact_render_font_size ) );
          chart.setPlotAreaPadding( 5, Wt::Left );
          chart.setPlotAreaPadding( 20, Wt::Bottom );
          chart.setPlotAreaPadding( 0, Wt::Right );
          titleFont.setSize( ns_compact_render_font_size );
          chart.axis( Chart::XAxis ).setTitle( "(keV)" );
          chart.setCompactAxis( true );
          break;
      }//switch( render_type )

      chart.axis( Chart::XAxis ).setLabelFont( labelFont );
      chart.axis( Chart::YAxis ).setLabelFont( labelFont );
      chart.axis( Chart::XAxis ).setTitleFont( labelFont );
      chart.axis( Chart::YAxis ).setTitleFont( labelFont );
      chart.setTitleFont( titleFont );
    }//if( thumbnail ) / else

    const size_t nchannel = meas->num_gamma_channels();
    const float lowx = meas->gamma_channel_lower( 0 );
    const float upperx = meas->gamma_channel_upper( nchannel - 1 );
    chart.setXAxisRange( lowx, upperx );

    const size_t displayednbin = meas->find_gamma_channel( upperx ) - meas->find_gamma_channel( lowx );
    const int plotAreaWidth = static_cast<int>( image.width().toPixels() )
                              - chart.plotAreaPadding( Left )
                              - chart.plotAreaPadding( Right );
    const float bins_per_pixel = float( displayednbin ) / float( plotAreaWidth );
    const int factor = max( static_cast<int>( ceil( bins_per_pixel / 2 ) ), 1 );
    dataModel.setRebinFactor( factor );

    chart.setYAxisScale( Wt::Chart::LogScale );
    chart.setAutoYAxisRange();

    if( !thumbnail && (width >= 641) )
    {
      if( factor == 1 )
        chart.axis( Chart::YAxis ).setTitle( "Counts Per Channel" );
      else
        chart.axis( Chart::YAxis ).setTitle( "Counts Per " + std::to_string( factor ) + " Channels" );
    }

    WPainter p( &image );
    chart.paint( p );

    // Draw watermark logo in upper-right for thumbnails
    if( thumbnail && !logo_path.empty() )
    {
      WPainter::Image logo( logo_path, logo_path );
      if( logo.width() > 0 && logo.height() > 0 )
      {
        const double img_w = image.width().toPixels();
        const double img_h = image.height().toPixels();

        // Logo at ~35% of the image height, positioned in upper-right
        const double logo_h = 0.35 * img_h;
        const double logo_w = logo_h * (double( logo.width() ) / double( logo.height() ));
        const double logo_x = img_w - logo_w - 4;
        const double logo_y = 4;

        p.drawImage( WRectF( logo_x, logo_y, logo_w, logo_h ), logo );
      }
    }

    p.end();

    return true;
  }//renderSpectrum(...)


  /** Renders a time series chart to a WPaintDevice. */
  template <class ImageType>
  bool renderTimeSeries( ImageType &image, std::shared_ptr<SpecMeas> meas )
  {
    const float width = image.width().toPixels();
    const SpecRenderType render_type = ((width < 641) ? CompactSpecRender : LargeSpecRender);

    const std::set<int> &sample_numbers = meas->sample_numbers();
    const std::vector<int> sample_numbers_vec( sample_numbers.begin(), sample_numbers.end() );
    const size_t num_samples = sample_numbers_vec.size();
    const std::vector<int> &det_numbers = meas->detector_numbers();
    const vector<string> &det_names = meas->detector_names();

    if( num_samples < 2 )
      return false;

    auto gamma_counts = std::make_shared<vector<float>>( num_samples );
    auto neutron_counts = std::make_shared<vector<float>>( num_samples );
    auto time_values = std::make_shared<vector<float>>( num_samples + 1, 0.0f );
    vector<bool> is_occupied( num_samples, false );

    bool contained_neutrons = false;
    float time_sum = 0.0;

    for( size_t i = 0; i < sample_numbers_vec.size(); ++i )
    {
      bool is_occ = false;
      int ndet = 0;
      for( const int detnum : det_numbers )
      {
        const auto thismeas = meas->measurement( sample_numbers_vec[i], detnum );
        if( thismeas && thismeas->real_time() > 0.0 )
          ++ndet;
        if( thismeas )
          is_occ = (is_occ || (thismeas->occupied() == SpecUtils::OccupancyStatus::Occupied));
      }
      if( !ndet )
        ndet = 1;

      set<int> thisset;
      thisset.insert( sample_numbers_vec[i] );
      const auto sum = meas->sum_measurements( thisset, det_names, nullptr );
      if( !sum )
        continue;

      float real_time = sum->real_time() / ndet;
      real_time = (real_time > 0.0f ? real_time : 0.1f);

      contained_neutrons = (contained_neutrons || sum->contained_neutron());
      (*gamma_counts)[i] = sum->gamma_count_sum() / real_time;
      (*neutron_counts)[i] = sum->neutron_counts_sum() / real_time;
      is_occupied[i] = is_occ;
      (*time_values)[i] = time_sum;
      time_sum += real_time;
    }//for( size_t i = 0; ... )

    (*time_values)[num_samples] = time_sum;

    auto gammaH = std::make_shared<SpecUtils::Measurement>();
    auto neutronH = std::make_shared<SpecUtils::Measurement>();

    gammaH->set_title( "Gammas" );
    neutronH->set_title( "Neutrons" );

    gammaH->set_gamma_counts( gamma_counts, 0.0f, 0.0f );
    neutronH->set_gamma_counts( neutron_counts, 0.0f, 0.0f );

    auto timecal = make_shared<SpecUtils::EnergyCalibration>();
    timecal->set_lower_channel_energy( num_samples, *time_values );

    gammaH->set_energy_calibration( timecal );
    neutronH->set_energy_calibration( timecal );

    SpectrumDataModel dataModel;
    PreviewSpectrumChart chart;
    chart.setModel( &dataModel );
    chart.setPlotAreaPadding( 0, Wt::Top );

    WFont labelFont = chart.axis( Chart::YAxis ).labelFont();

    switch( render_type )
    {
      case LargeSpecRender:
        labelFont.setSize( WLength( ns_large_render_font_size ) );
        chart.setPlotAreaPadding( 62, Wt::Left );
        chart.setPlotAreaPadding( 42, Wt::Bottom );
        chart.setPlotAreaPadding( (contained_neutrons ? 20 : 5), Wt::Right );
        chart.axis( Chart::XAxis ).setTitle( "Time (s)" );
        chart.setCompactAxis( false );
        break;

      case CompactSpecRender:
        labelFont.setSize( WLength( ns_compact_render_font_size ) );
        chart.setPlotAreaPadding( 5, Wt::Left );
        chart.setPlotAreaPadding( 20, Wt::Bottom );
        chart.setPlotAreaPadding( 5, Wt::Right );
        chart.axis( Chart::XAxis ).setTitle( "(sec.)" );
        chart.setCompactAxis( true );
        break;
    }//switch( render_type )

    chart.axis( Chart::XAxis ).setLabelFont( labelFont );
    chart.axis( Chart::YAxis ).setLabelFont( labelFont );
    chart.axis( Chart::Y2Axis ).setLabelFont( labelFont );
    chart.axis( Chart::XAxis ).setTitleFont( labelFont );
    chart.axis( Chart::YAxis ).setTitleFont( labelFont );
    chart.axis( Chart::Y2Axis ).setTitleFont( labelFont );

    dataModel.setDataHistogram( gammaH );
    if( contained_neutrons )
    {
      dataModel.setSecondDataHistogram( neutronH, true );
      chart.axis( Chart::Y2Axis ).setVisible( true );
    }

    // Highlight occupied time regions
    vector<pair<double, double>> timeranges;
    {
      const size_t ntimes = std::min( time_values->size(), is_occupied.size() );
      bool inregion = (ntimes > 0) ? is_occupied[0] : false;
      double startt = (ntimes > 0) ? (*time_values)[0] : 0.0;
      for( size_t i = 1; i < ntimes; ++i )
      {
        if( inregion == is_occupied[i] )
          continue;
        if( !is_occupied[i] )
          timeranges.push_back( make_pair( startt, (*time_values)[i] ) );
        inregion = is_occupied[i];
        startt = (*time_values)[i];
      }
      if( inregion && time_values->size() && (!timeranges.empty() || is_occupied[0] != is_occupied.back()) )
        timeranges.push_back( make_pair( startt, time_values->back() ) );
    }

    chart.setTimeHighLightRegions( timeranges, SpectrumChart::ForegroundHighlight );

    const vector<Chart::WDataSeries> series = dataModel.suggestDataSeries();
    chart.setSeries( series );

    chart.setXAxisRange( time_values->at( 0 ), time_values->back() );
    chart.axis( Chart::YAxis ).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
    if( contained_neutrons )
      chart.axis( Chart::Y2Axis ).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );

    WPainter p( &image );
    chart.paint( p );
    p.end();

    return true;
  }//renderTimeSeries(...)
}//anonymous namespace


void render_spec_file_to_pdf( uint8_t **result, size_t *result_size,
                              const char * const filename,
                              const float width_pt, const float height_pt,
                              const enum SpecPreviewType type,
                              const char * const logo_path )
{
  (*result) = nullptr;
  (*result_size) = 0;

  auto spec = std::make_shared<SpecMeas>();

  if( !spec->load_file( filename, SpecUtils::ParserType::Auto, filename ) )
  {
    fprintf( stderr, "Failed to parse '%s' as a spectrum file.\n", filename );
    return;
  }

  if( !spec->num_measurements() )
  {
    fprintf( stderr, "No measurements in '%s'.\n", filename );
    return;
  }

  const std::set<int> &samplenums = spec->sample_numbers();

  std::set<int> foreground_sample_numbers, background_sample_numbers;
  std::shared_ptr<const SpecUtils::Measurement> foreground, background;

  // Determine foreground and background sample numbers
  if( spec->passthrough() )
  {
    for( const int sample : samplenums )
    {
      bool occupied = false;
      for( const int detnum : spec->detector_numbers() )
      {
        const auto m = spec->measurement( sample, detnum );
        occupied = (occupied
                    || (m && (m->occupied() == SpecUtils::OccupancyStatus::Occupied
                              || m->source_type() == SpecUtils::SourceType::Foreground)));
      }
      if( occupied )
        foreground_sample_numbers.insert( sample );
      else
        background_sample_numbers.insert( sample );
    }
  }else if( spec->num_measurements() == 1 || samplenums.size() == 1 )
  {
    foreground_sample_numbers = samplenums;
  }else
  {
    bool hasBackground = false, hasForeground = false;
    int childrow = -9999, backrow = -9999;

    for( const int sample : samplenums )
    {
      bool has_meas = false;
      SpecUtils::SourceType type = SpecUtils::SourceType::Unknown;
      for( const int detnum : spec->detector_numbers() )
      {
        const auto m = spec->measurement( sample, detnum );
        if( m && m->num_gamma_channels() )
        {
          has_meas = true;
          type = m->source_type();
        }
      }

      if( !has_meas )
        continue;

      switch( type )
      {
        case SpecUtils::SourceType::IntrinsicActivity:
        case SpecUtils::SourceType::Calibration:
          break;

        case SpecUtils::SourceType::Foreground:
          if( !hasForeground )
          {
            hasForeground = true;
            childrow = sample;
          }
          break;

        case SpecUtils::SourceType::Background:
          if( !hasBackground )
          {
            hasBackground = true;
            backrow = sample;
          }
          break;

        case SpecUtils::SourceType::Unknown:
          if( childrow < -999 )
            childrow = sample;
          break;
      }//switch( type )
    }//for( const int sample : samplenums )

    if( childrow < -999 )
      childrow = *samplenums.begin();

    foreground_sample_numbers.insert( childrow );
    if( hasBackground )
      background_sample_numbers.insert( backrow );
  }//if( passthrough ) / elif( single ) / else

  if( foreground_sample_numbers.empty() )
    swap( foreground_sample_numbers, background_sample_numbers );

  if( foreground_sample_numbers.empty() )
  {
    fprintf( stderr, "No foreground measurements in '%s'.\n", filename );
    return;
  }

  const vector<string> &detectors_to_use = spec->detector_names();

  if( foreground_sample_numbers.size() == 1 && detectors_to_use.size() == 1 )
  {
    foreground = spec->measurement( *foreground_sample_numbers.begin(), detectors_to_use.at( 0 ) );
  }else
  {
    foreground = spec->sum_measurements( foreground_sample_numbers, detectors_to_use, nullptr );
  }

  if( !background_sample_numbers.empty() )
    background = spec->sum_measurements( background_sample_numbers, detectors_to_use, nullptr );

  if( !foreground )
  {
    fprintf( stderr, "No measurements in '%s'.\n", filename );
    return;
  }

  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = spec->peaks( foreground_sample_numbers );

  // Render to PDF
  // Use 1:1 point-to-pixel ratio; PDFView auto-scales to fit the window
  const float width_px = width_pt;
  const float height_px = height_pt;

  HPDF_Doc pdfdoc = nullptr;

  try
  {
    pdfdoc = HPDF_New( hpdf_error_handler, nullptr );
    if( !pdfdoc )
      throw std::runtime_error( "Could not create libharu document." );

    HPDF_SetCompressionMode( pdfdoc, HPDF_COMP_ALL );
    HPDF_Page pdfpage = HPDF_AddPage( pdfdoc );
    HPDF_Page_SetWidth( pdfpage, width_px );
    HPDF_Page_SetHeight( pdfpage, height_px );
    HPDF_Page_GSave( pdfpage );
    HPDF_UseUTFEncodings( pdfdoc );

    const bool is_thumbnail = (type == SpectrumThumbnail);

    if( !is_thumbnail && spec->passthrough() )
    {
      // Preview of passthrough: spectrum on top 2/3, time series on bottom 1/3
      const float spec_height_px = 2.0f * height_px / 3.0f;

      auto spectrum_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, -spec_height_px / 2.0f, width_px, spec_height_px );
      if( !renderSpectrum( *spectrum_img, foreground, background, peaks, false ) )
        throw runtime_error( "Failed in renderSpectrum for passthrough" );

      auto time_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, 0, width_px, spec_height_px / 2.0f );
      if( !renderTimeSeries( *time_img, spec ) )
        throw runtime_error( "Failed in renderTimeSeries" );
    }else
    {
      // For thumbnails: always render just the spectrum filling the full area
      // For non-passthrough previews: same behavior
      const std::string logo = (logo_path ? logo_path : "");
      auto spectrum_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, 0, width_px, height_px );
      if( !renderSpectrum( *spectrum_img, foreground, background, peaks, is_thumbnail, logo ) )
        throw runtime_error( "Failed in renderSpectrum" );
    }

    // Stream PDF to buffer
    static_assert( sizeof( HPDF_BYTE ) == 1, "HPDF_BYTE must be exactly one byte" );
    vector<HPDF_BYTE> resultdata;
    HPDF_SaveToStream( pdfdoc );
    HPDF_ResetStream( pdfdoc );

    for( ;; )
    {
      HPDF_BYTE buf[4096];
      HPDF_UINT32 siz = 4096;
      HPDF_ReadFromStream( pdfdoc, buf, &siz );
      if( siz == 0 )
        break;
      resultdata.insert( resultdata.end(), buf, buf + siz );
    }

    if( resultdata.empty() )
      throw std::runtime_error( "Failed to stream PDF" );

    HPDF_Free( pdfdoc );
    pdfdoc = nullptr;

    (*result) = (uint8_t *)malloc( resultdata.size() );
    memcpy( (*result), &(resultdata[0]), resultdata.size() );
    (*result_size) = resultdata.size();
  }catch( std::exception &e )
  {
    fprintf( stderr, "Failed rendering PDF: %s\n", e.what() );
    if( pdfdoc )
      HPDF_Free( pdfdoc );
    return;
  }
}//render_spec_file_to_pdf(...)
