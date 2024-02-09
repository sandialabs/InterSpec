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
//
//  Created by Johnson, William C on 12/25/17.
//
#include "SpecPreviewCommon.h"

#include <stdio.h>

#include <string>
#include <sstream>

#include <Wt/Utils>
#include <Wt/WFont>
#include <Wt/WPainter>
#include <Wt/WSvgImage>
#if( RENDER_PREVIEWS_AS_PDF )
#include <hpdf.h>
#include <Wt/WPdfImage>
#endif

#include "QLSpecMeas.h"
#include "QLSpectrumChart.h"
#include "SpecUtils/SpecFile.h"
#include "QLSpectrumDataModel.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

using namespace std;
using namespace Wt;
using SpecUtils::Measurement;

enum SpecRenderType
{
  CompactSpecRender,
  LargeSpecRender
};
const float ns_large_render_font_size = 64;
const float ns_compact_render_font_size = 32;
const float ns_large_render_legend_font_size = 48;
const float ns_compact_render_legend_font_size = 24;



namespace
{
#if( RENDER_PREVIEWS_AS_PDF )
  void HPDF_STDCALL hpdf_error_handler( HPDF_STATUS error_no, HPDF_STATUS detail_no, void *user_data )
  {
    //Wt::WPdfImage *image = (Wt::WPdfImage *)user_data;
    //image->errorHandler(error_no, detail_no);
    printf("Pdf image error: error_no=%04X, detail_no=%d",
           (unsigned int) error_no, (int) detail_no);
    
    char buf[200];
    snprintf(buf, 200, "Pdf image error: error_no=%04X, detail_no=%d",
             (unsigned int) error_no, (int) detail_no);
    
    throw std::runtime_error( buf );
  }//hpdf_error_handler(...)
#endif
  
}//namespace

vector< pair<double,double> > timeRegionsToHighlight( const vector<float> &binning, const vector<bool> &status )
{
  vector< pair<double,double> > timeranges;
  
  const size_t ntimes = std::min( binning.size(), status.size() );
  
  if( !ntimes )
    return timeranges;
  
  bool inregion = status[0];
  double startt = binning[0];
  for( size_t i = 1; i < ntimes; ++i )
  {
    if( inregion == status[i] )
      continue;
    
    if( !status[i] )
      timeranges.push_back( make_pair(startt, binning[i]) );
    
    inregion = status[i];
    startt = binning[i];
  }//for( size_t i = 0; i < binning.size(); ++i )
  
  if( inregion && binning.size() && (!timeranges.empty() || status[0]!=status.back()) )
    timeranges.push_back( make_pair(startt, binning.back()) );
  
  return timeranges;
}//timeRegionsToHighlight(...)


//Returns true if successful, false otherwise.
template <class ImageType>
bool renderTimeSeries( ImageType &image, std::shared_ptr<QLSpecMeas> meas )
{
  const float width = image.width().toPixels();
  const float height = image.height().toPixels();
  
  const SpecRenderType render_type = ((width < 641) ? CompactSpecRender : LargeSpecRender);
  
  const std::set<int> &sample_numbers = meas->sample_numbers();
  const std::vector<int> sample_numbers_vec( sample_numbers.begin(), sample_numbers.end() );
  
  const size_t num_samples = sample_numbers_vec.size();
  const std::vector<int> &det_numbers = meas->detector_numbers();
  const vector<string> &det_names = meas->detector_names();
  
  if( num_samples < 2 )
  {
    printf( "Not enough time samples to render time series\n" );
    return false;
  }
  
  auto gamma_counts = std::make_shared< vector<float> >( num_samples );
  auto neutron_counts = std::make_shared< vector<float> >( num_samples );
  auto time_values = std::make_shared< vector<float> >( num_samples + 1, 0.0f );
  vector<bool> is_occupied( num_samples, false );
  
  bool contained_neutrons = false;
  float time_sum = 0.0;
  for( size_t i = 0; i < sample_numbers_vec.size(); ++i )
  {
    bool is_occ = false;
    int ndet = 0;
    for( int detnum : det_numbers )
    {
      auto thismeas = meas->measurement( sample_numbers_vec[i], detnum );
      if( thismeas && thismeas->real_time() > 0.0 )
        ++ndet;
      if( thismeas )
        is_occ = (is_occ || (thismeas->occupied() == SpecUtils::OccupancyStatus::Occupied));
    }
    if( !ndet )
      ndet = 1;
    
    set<int> thisset;
    thisset.insert( sample_numbers_vec[i] );
    auto sum = meas->sum_measurements( thisset, det_names, nullptr );
    if( !sum )
      continue;
    
    float real_time = sum->real_time() / ndet;
    real_time = (real_time>0.0f ? real_time : 0.1f);
    
    contained_neutrons = (contained_neutrons || sum->contained_neutron());
    (*gamma_counts)[i] = sum->gamma_count_sum() / real_time;
    (*neutron_counts)[i] = sum->neutron_counts_sum() / real_time;
    is_occupied[i] = is_occ;//(sum->occupied() == SpecUtils::Measurement::Occupied);
    (*time_values)[i] = time_sum;
    time_sum += real_time;
  }//for( const int sample : sample_numbers )
  
  (*time_values)[num_samples] = time_sum;
  
  auto gammaH   = std::make_shared<SpecUtils::Measurement>();
  auto neutronH = std::make_shared<SpecUtils::Measurement>();
    
  gammaH->set_title( "Gammas" );
  neutronH->set_title( "Neutrons" );
  
    
  gammaH->set_gamma_counts( gamma_counts, 0.0f, 0.0f );
  neutronH->set_gamma_counts( neutron_counts, 0.0f, 0.0f );
  
  auto timecal = make_shared<SpecUtils::EnergyCalibration>();
  timecal->set_lower_channel_energy( num_samples, *time_values );
  
  gammaH->set_energy_calibration( timecal );
  neutronH->set_energy_calibration( timecal );
  
  
  
  QLSpectrumDataModel dataModel;
  QLSpectrumChart chart;
  //chart.setPlotAreaPadding( 80, 10, 40, 42 );
  //chart.legendTextSizeCallback( const float legendwidth );
  
  chart.setModel( &dataModel );
  //chart.axis(Chart::YAxis).setTitle( "CPS" );
  //chart.setY2AxisTitle( "Neutron CPS" );
  //chart.enableLegend( true );
  chart.setPlotAreaPadding( 0, Wt::Top );
  
  WFont labelFont = chart.axis(Chart::YAxis).labelFont();
  
  switch( render_type )
  {
    case LargeSpecRender:
      labelFont.setSize( WLength(ns_large_render_font_size) );
      chart.setPlotAreaPadding( 165, Wt::Bottom );
      chart.setPlotAreaPadding( (contained_neutrons ? 75 : 5), Wt::Right );
      //chart.enableLegend( true );
      //chart.setPaintOnLegendFontSize( ns_large_render_legend_font_size );
      chart.axis(Chart::XAxis).setTitle( "Time (s)" );
      chart.setCompactAxis( false );
      break;
      
    case CompactSpecRender:
      labelFont.setSize( WLength(ns_compact_render_font_size) );
      chart.setPlotAreaPadding( 15, Wt::Left );
      chart.setPlotAreaPadding( 75, Wt::Bottom );
      chart.setPlotAreaPadding( (contained_neutrons ? 5 : 5), Wt::Right );
      chart.setPaintOnLegendFontSize( ns_compact_render_legend_font_size );
      chart.axis(Chart::XAxis).setTitle( "(sec.)" );
      chart.setCompactAxis( true );
      break;
  }//switch( type )
  
  chart.axis(Chart::XAxis).setLabelFont( labelFont );
  chart.axis(Chart::YAxis).setLabelFont( labelFont );
  chart.axis(Chart::Y2Axis).setLabelFont( labelFont );
  
  chart.axis(Chart::XAxis).setTitleFont( labelFont );
  chart.axis(Chart::YAxis).setTitleFont( labelFont );
  chart.axis(Chart::Y2Axis).setTitleFont( labelFont );
  

  dataModel.setDataHistogram( gammaH, -1, -1, -1 );
  if( contained_neutrons )
  {
    dataModel.setSecondDataHistogram( neutronH, -1.0, -1.0, -1.0, true );
    chart.axis(Chart::Y2Axis).setVisible( true );
  }
  
  auto regions = timeRegionsToHighlight( *time_values, is_occupied );
  chart.setTimeHighLightRegions( regions, QLSpectrumChart::ForegroundHighlight );
  
  const vector<Chart::WDataSeries> series = dataModel.suggestDataSeries();
  chart.setSeries( series );
  
  chart.setXAxisRange( time_values->at(0), time_values->back() );
  chart.axis(Chart::YAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  if( contained_neutrons )
    chart.axis(Chart::Y2Axis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  
#if(DYNAMICALLY_ADJUST_LEFT_CHART_PADDING)
  if( render_type == LargeSpecRender )
  {
    int pad_adj_status = chart.setLeftYAxisPadding( width, height );
    if( pad_adj_status < 0 )
      cerr << "Warning: failed to adjust time series left axis padding" << endl;
    
    if( contained_neutrons )
    {
      pad_adj_status = chart.setRightYAxisPadding( width, height );
      if( pad_adj_status < 0 )
        cerr << "Warning: failed to adjust time series right axis padding" << endl;
    }
  }
#endif
  
  WPainter p( &image );
  chart.paint( p );
  p.end();
  
  return true;
}//renderTimeSeries(....)



//Returns true if successful, false otherwise.
template <class ImageType>
bool renderSpectrum( ImageType &image,
                     std::shared_ptr<SpecUtils::Measurement> foreground,
                     std::shared_ptr<SpecUtils::Measurement> background,
                     std::shared_ptr< std::deque<std::shared_ptr<const QLPeakDef> > > peaks,
                     const bool legend )
{
  if( !foreground )
  {
    cerr << "renderSpectrum: invalid input" << endl;
    return false;
  }
  
  const std::string title = foreground->title();
  
  const float width = image.width().toPixels();
  const float height = image.height().toPixels();
  
  std::shared_ptr<SpecUtils::Measurement> meas, back;
  meas = std::make_shared<SpecUtils::Measurement>( *foreground );
  if( background )
    back = std::make_shared<SpecUtils::Measurement>( *background );
  
  if( legend )
  {
    //Although, I dont think this actually has an effect.
    if( !back )
    {
      meas->set_title( "" );
    }else
    {
      meas->set_title( "Foreground" );
      back->set_title( "Background" );
    }
  }//if( legend )
  
  //hack: calling qlmanage from command line always gives thumbnail
  //cout << "w=" << width << ", h=" << height << endl;
  const SpecRenderType render_type = ((width < 641) ? CompactSpecRender : LargeSpecRender);
  
  //printf( "128 pixels -> %f\n",  WLength(128,WLength::Point).toPixels() );
  //128 pixels -> 170.666667
  std::shared_ptr<QLSpecMeas> specmeas = std::make_shared<QLSpecMeas>();
  specmeas->add_measurement( meas, true );
  if( peaks )
    specmeas->setPeaks( *peaks, specmeas->sample_numbers() );
  
  QLSpectrumDataModel dataModel;
  QLSpectrumChart chart;
  chart.setModel( &dataModel );
  chart.setPeaks( peaks );
  
  const float liveTime = meas->live_time();
  const float realTime = meas->real_time();
  const double neutrons = meas->neutron_counts_sum();
  dataModel.setDataHistogram( meas, liveTime, realTime, neutrons );
  if( back )
  {
    float sf = 1.0;
    if( meas->live_time() > 0.0f && back->live_time() > 0.0f )
      sf = meas->live_time() / back->live_time();
    dataModel.setBackgroundHistogram( back, back->live_time(), back->real_time(), back->neutron_counts_sum() );
    dataModel.setBackgroundDataScaleFactor( sf );
  }//if( background )
  const vector<Chart::WDataSeries> series = dataModel.suggestDataSeries();
  chart.setSeries( series );
  
  chart.setPlotAreaPadding( 0, Wt::Top );
  
  WFont labelFont = chart.axis(Chart::YAxis).labelFont();
  WFont titleFont = chart.titleFont();
  switch( render_type )
  {
    case LargeSpecRender:
      labelFont.setSize( WLength(ns_large_render_font_size) );
      chart.setPlotAreaPadding( 250, Wt::Left );
      chart.setPlotAreaPadding( 165, Wt::Bottom );
      chart.setPlotAreaPadding( 35, Wt::Right );  //make room so the last energy label doesnt get truncated
      if( title.size() )
      {
        chart.setPlotAreaPadding( ns_large_render_font_size + 5, Wt::Top );
        titleFont.setSize( ns_large_render_font_size );
      }
      //chart.setPaintOnLegendFontSize( ns_large_render_font_size );
      chart.axis(Chart::XAxis).setTitle( "Energy (keV)" );
      chart.setCompactAxis( false );
      if( title.size() )
        chart.setTitle( title );
    break;
      
    case CompactSpecRender:
      labelFont.setSize( WLength(ns_compact_render_font_size) );
      chart.setPlotAreaPadding( 15, Wt::Left );
      chart.setPlotAreaPadding( 75, Wt::Bottom );
      chart.setPlotAreaPadding( 0, Wt::Right );
      titleFont.setSize( ns_compact_render_font_size );
      //chart.setPaintOnLegendFontSize( ns_compact_render_font_size );
      chart.axis(Chart::XAxis).setTitle( "(keV)" );
      chart.setCompactAxis( true );
    break;
  }//switch( render_type )
  
  chart.axis(Chart::XAxis).setLabelFont( labelFont );
  chart.axis(Chart::YAxis).setLabelFont( labelFont );
  
  chart.axis(Chart::XAxis).setTitleFont( labelFont );
  chart.axis(Chart::YAxis).setTitleFont( labelFont );
  
  chart.setTitleFont( titleFont );
  
  //if( legend )
  //  chart.enableLegend( true );
  
  const size_t nchannel = meas->num_gamma_channels();
  float lowx = meas->gamma_channel_lower( 0 );
  float upperx = meas->gamma_channel_upper( nchannel - 1 );
  
  chart.setXAxisRange( lowx, upperx );
  
  const size_t displayednbin = meas->find_gamma_channel( upperx )
  - meas->find_gamma_channel( lowx );
  
  const int plotAreaWidth = static_cast<int>( image.width().toPixels() )
  - chart.plotAreaPadding(Left)
  - chart.plotAreaPadding(Right);
  const float bins_per_pixel = float(displayednbin) / float(plotAreaWidth);
  const int factor = max( static_cast<int>(ceil(bins_per_pixel/2)), 1 );  //QuickLook specialization of allowing 2 bins per Px (e.g. retina screens or somethign maybe)
  
  //cout << "plotAreaWidth=" << plotAreaWidth << ", img->width().toPixels()=" << img->width().toPixels()
  //     << ", chart.plotAreaPadding(Left)=" << chart.plotAreaPadding(Left) << ", chart.plotAreaPadding(Right)="
  //     << chart.plotAreaPadding(Right) << endl;
  //cout << "bins_per_pixel=" << bins_per_pixel << ", factor=" << factor << endl;
  dataModel.setRebinFactor( factor );
  
  chart.setYAxisScale( Wt::Chart::LogScale );
  chart.setAutoYAxisRange();
  
  switch( render_type )
  {
    case CompactSpecRender:
      chart.setShowYLabels( false );
      break;
      
    case LargeSpecRender:
      if( factor == 1 )
        chart.axis(Chart::YAxis).setTitle( "Counts Per Channel" );
      else
        chart.axis(Chart::YAxis).setTitle( "Counts Per " + std::to_string(factor) + " Channels" );
  }//switch( render_type )
  
  //Normally padding adjust is called from QLSpectrumChart::paintEvent(), which
  //  doesnt get called in this application
#if(DYNAMICALLY_ADJUST_LEFT_CHART_PADDING)
  if( render_type == LargeSpecRender )
  {
    const int pad_adj_status = chart.setLeftYAxisPadding( width, height );
    if( pad_adj_status < 0 )
      cerr << "Warning: failed to adjust left axis padding" << endl;
  }
#endif
  
  WPainter p( &image );
  chart.paint( p );
  p.end();
  
  return true;
}//renderSpectrum(...)



void render_spec_file_to_preview( uint8_t **result, size_t *result_size,
                                  const char * const filename, float width_pt, float height_pt, const SpecPreviewType type )
{
  (*result) = nullptr;
  (*result_size) = 0;
  
  std::shared_ptr<QLSpecMeas> spec = std::make_shared<QLSpecMeas>();
  std::shared_ptr< std::deque<std::shared_ptr<const QLPeakDef> > > peaks;
  
  if( !spec->load_file( filename, SpecUtils::ParserType::Auto, filename ) )
  {
    printf( "Failed to parse '%s' as a spectrum file.\n", filename );
    return;
  }//

  if( !spec->num_measurements() )
  {
    printf( "No measurements in '%s'.\n", filename );
    return;
  }
  
  const std::set<int> &samplenums = spec->sample_numbers();
  
  std::set<int> foreground_sample_numbers, background_sample_numbers, intrinsic_sample_numbers;
  std::shared_ptr<SpecUtils::Measurement> foreground, background;

  if( spec->passthrough() )
  {
    for( int sample : samplenums )
    {
      bool occupied = false;
      for( int detnum : spec->detector_numbers() )
      {
        auto m = spec->measurement( sample, detnum );
        occupied = (occupied || (m && (m->occupied()== SpecUtils::OccupancyStatus::Occupied || m->source_type()==SpecUtils::SourceType::Foreground)));
      }
      if( occupied )
        foreground_sample_numbers.insert( sample );
      else
        background_sample_numbers.insert( sample );
    }//for( int sample : samplenums )
  }else
  {
    if( spec->num_measurements() == 1 || samplenums.size() == 1 )
    {
      foreground_sample_numbers = samplenums;
    }else
    {
      bool hasIntrinsic = false, hasBackground = false, hasForeground = false;
      int childrow = -9999, backrow = -9999, intrinsicrow = -9999;
      for( int sample : samplenums )
      {
        bool has_meas = false;
        SpecUtils::SourceType type = SpecUtils::SourceType::Unknown;
        
        for( int detnum : spec->detector_numbers() )
        {
          auto m = spec->measurement( sample, detnum );
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
            //lets try to not show intrinsic activity by default
            if( hasIntrinsic )
              continue;
            hasIntrinsic = true;
            intrinsicrow = sample;
            break;
            
          case SpecUtils::SourceType::Foreground:
            if( hasForeground )
              continue;
            hasForeground = true;
            childrow = sample;
            break;
            
          case SpecUtils::SourceType::Background:
            if( hasBackground )
              continue;
            hasBackground = true;
            backrow = sample;
            break;
            
          case SpecUtils::SourceType::Calibration:
            //do nothing
            break;
            
          case SpecUtils::SourceType::Unknown:
            if( childrow < -999 )
              childrow = sample;
            break;
        }//switch( header->m_samples[i].spectra_type )
      }//for( int sample : samplenums )
      
      if( (childrow < -999) )
        childrow = *samplenums.begin();
      
      foreground_sample_numbers.insert( childrow );
      if( hasBackground )
        background_sample_numbers.insert( backrow );
      if( hasIntrinsic )
        intrinsic_sample_numbers.insert( intrinsicrow );
    }//if( spec->num_measurements() == 1 || samplenums.size() == 1 ) / else
  }//if( spec->passthrough() ) / else
  
  
  if( foreground_sample_numbers.empty() )
    swap( foreground_sample_numbers, background_sample_numbers );
  
  if( foreground_sample_numbers.empty() )
  {
    printf( "No foreground measurements in '%s' spectrum file.\n", filename );
    return;
  }
  

  const vector<string> &detectors_to_use = spec->detector_names();
  
  if( foreground_sample_numbers.size()==1 && detectors_to_use.size()==1 )
  {
    auto m = spec->measurement( *foreground_sample_numbers.begin(), spec->detector_names().at(0) );
    if( m )
      foreground = std::make_shared<Measurement>( *m );
  }else
  {
    foreground = spec->sum_measurements( foreground_sample_numbers, detectors_to_use, nullptr );
  }
  
  
  if( background_sample_numbers.size() )
    background = spec->sum_measurements( background_sample_numbers, detectors_to_use, nullptr );
  
  if( !foreground )
  {
    printf( "No measurements in '%s' spectrum file.\n", filename );
    return;
  }
  
  peaks = spec->peaks( foreground_sample_numbers );
  
  if( foreground->title() == spec->filename() )
    foreground->set_title( "" );
  
  //My macbook 15" has a DPI of 220.
  // A point is 1/72 of an inch.
  const float resolution_multiple = 220.0f / 72.0f;  //5 works well...
  const float width_px = resolution_multiple * width_pt;
  const float height_px = resolution_multiple * height_pt;
  
  const bool show_legend = (type == SpectrumPreview && width_px >= 640);  //would work, but need to edit font...
  
#if( RENDER_PREVIEWS_AS_PDF )
  if( type == SpectrumThumbnail )
  {
#endif
    try
    {
      const float spec_height_px = (spec->passthrough() ? 2.0f*height_px/3.0f : height_px);
      
      std::shared_ptr<Wt::WSvgImage> spectrum_img, time_img;
      spectrum_img = std::make_shared<WSvgImage>( width_px, spec_height_px );
    
      if( !renderSpectrum( *spectrum_img, foreground, background, peaks, show_legend ) )
      {
        cerr << "Failed in renderSpectrum" << endl;
        return;
      }
      if( spec->passthrough() )
      {
        time_img = std::make_shared<WSvgImage>( width_px, spec_height_px/2.0f );
      
        if( !renderTimeSeries( *time_img, spec ) )
        {
          cerr << "Failed in renderTimeSeries" << endl;
          return;
        }
      }//if( spec->passthrough() )
    
      string fname = SpecUtils::filename(filename);
      stringstream resultstrm;
      fname = Wt::Utils::htmlEncode(fname);
    
      resultstrm << "<html>\n<head>\n\t"
      "<title>" << fname << "</title>\n"
      "</head>\n"
      "<style>"
      "</style>"
      "\n<body>\n"
      "<div style=\"text-align:center;\">\n";
      spectrum_img->write( resultstrm );
      resultstrm << "\n</div>\n";
    
      if( time_img )
      {
        resultstrm << "<div style=\"text-align:center;\">\n";
        time_img->write( resultstrm );
        resultstrm << "\n</div>\n";
      }
    
      resultstrm << "</body></html>" << endl;
    
      const string resultstr = resultstrm.str();
      (*result) = (uint8_t *)malloc(resultstr.size() + 1);
      memcpy( (*result), resultstr.c_str(), resultstr.size() + 1 );
      (*result_size) = resultstr.size();
      return;
    }catch( std::exception &e )
    {
      cerr << "Caught exception rendering SVG: " << e.what() << endl;
      return;
    }//try / catch
    
#if( RENDER_PREVIEWS_AS_PDF )
  }else
  {
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
      
      if( spec->passthrough() )
      {
        const float spec_height_px = 2.0f * height_px / 3.0f;
        
        auto spectrum_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, -spec_height_px/2.0f, width_px, spec_height_px );
        if( !renderSpectrum( *spectrum_img, foreground, background, peaks, show_legend ) )
          throw runtime_error( "Failed in renderSpectrum for passthrough" );
        
        auto time_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, 0, width_px, spec_height_px/2.0f );
        if( !renderTimeSeries( *time_img, spec ) )
          throw runtime_error( "Failed in renderTimeSeries" );
      }else
      {
        auto spectrum_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, 0, width_px, height_px );
        
        if( !renderSpectrum( *spectrum_img, foreground, background, peaks, show_legend ) )
          throw runtime_error( "Failed in renderSpectrum for passthrough" );
        
        set<int> remainingsamplenums = samplenums;
        for( int i : foreground_sample_numbers )
          remainingsamplenums.erase( i );
        for( int i : background_sample_numbers )
          remainingsamplenums.erase( i );
        //intrinsic_sample_numbers...
        
        if( remainingsamplenums.size() )
        {
          size_t npages = 1;
          size_t maxpages = remainingsamplenums.size() + 2;
          if( maxpages > 15 )
            maxpages = 10;
          
          for( const int samplenum : remainingsamplenums )
          {
            background.reset();
            foreground.reset();
            
            set<int> thissamplenums;
            thissamplenums.insert( samplenum );
            
            if( detectors_to_use.size()==1 )
            {
              auto m = spec->measurement( samplenum, spec->detector_names().at(0) );
              if( m )
                foreground = std::make_shared<Measurement>( *m );
            }else
            {
              foreground = spec->sum_measurements( thissamplenums, detectors_to_use, nullptr );
            }
            
            if( !foreground )
              continue;
            
            npages += 1;
            
            //Should set the title to indicate sample number and such
            
            peaks = spec->peaks( thissamplenums );
            
            pdfpage = HPDF_AddPage( pdfdoc );
            HPDF_Page_SetWidth( pdfpage, width_px );
            HPDF_Page_SetHeight( pdfpage, height_px );
            HPDF_Page_GSave( pdfpage );
            
            spectrum_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, 0, width_px, height_px );
            
            if( npages <= maxpages )
            {
              if( !renderSpectrum( *spectrum_img, foreground, background, peaks, show_legend ) )
                throw runtime_error( "Failed in renderSpectrum for passthrough" );
            }else
            {
              //Add in a page saying "And X more"
              const string message = "Plus " + to_string(remainingsamplenums.size() - maxpages) + " more spectra";
              
              QLSpectrumDataModel dataModel;
              QLSpectrumChart chart;
              chart.setModel( &dataModel );
              chart.setTextInMiddleOfChart( message.c_str() );
              chart.setShowYLabels( false );
              WPainter p( spectrum_img.get() );
              chart.paint( p );
              p.end();
              
              break;
              
              /*
               auto f = HPDF_GetFont( pdfdoc, "Helvetica", 0 );
               HPDF_REAL fsize = 48;
               pdfpage = HPDF_AddPage( pdfdoc );
               HPDF_Page_SetWidth( pdfpage, resolution_multiple * width_pt );
               HPDF_Page_SetHeight( pdfpage, resolution_multiple * height_pt );
               HPDF_Page_GSave( pdfpage );
               
               HPDF_Page_SetFontAndSize( pdfpage, f, fsize );
               
               // Undo the global inversion
               HPDF_Page_Concat(pdfpage, 1, 0, 0, -1, 0, resolution_multiple * height_pt/2);
               
               HPDF_Page_BeginText(pdfpage);
               
               // Need to fill text using pen color
               HPDF_Page_SetRGBFill( pdfpage, 0., 0., 0. ); //black
               
               HPDF_REAL left = 0.;
               HPDF_REAL top = 48;
               HPDF_REAL right = 10000;
               HPDF_REAL bottom = 0;
               
               //Currently text is printed upside down...
               HPDF_Page_TextRect( pdfpage, left, top, right, bottom, message.c_str(), HPDF_TALIGN_LEFT, 0);
               HPDF_Page_EndText(pdfpage);
               HPDF_Page_GRestore(pdfpage);
               
               break;
               */
            }
          }//for( const int samplenum : remainingsamplenums )
        }else
        {
          //Simple spectrum file with foreground/background
          auto spectrum_img = std::make_shared<WPdfImage>( pdfdoc, pdfpage, 0, 0, width_px, height_px );
          if( !renderSpectrum( *spectrum_img, foreground, background, peaks, show_legend ) )
            throw runtime_error( "Failed in renderSpectrum" );
        }
      }//if( spec->passthrough() ) / else
      
      // string fname = UtilityFunctions::filename(filename);
      
      static_assert( sizeof(HPDF_BYTE) == 1, "HPDF_BYTE must be exactly one byte" );
      vector<HPDF_BYTE> resultdata;
      HPDF_SaveToStream( pdfdoc );
      HPDF_ResetStream( pdfdoc );
      
      for( ; ; )
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
      
      //ba133_source_640s_20100317.n42 is 37090 bytes.
    }catch( std::exception &e )
    {
      cerr << "Failed rendering PDF: " << e.what() << endl;
      if( pdfdoc )
        HPDF_Free( pdfdoc );
      return;
    }
  }//if( Thumbnail ) / else
#endif
}//render_spec_file_to_html_preview(...)
