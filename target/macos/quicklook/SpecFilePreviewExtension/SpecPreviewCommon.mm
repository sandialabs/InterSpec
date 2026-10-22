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
// The charts themselves are drawn by SpectrumChartImage (shared with the stdio MCP server).

#include "SpecPreviewCommon.h"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <iostream>

#include "CGPaintDevice.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/SpectrumChartImage.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

using namespace std;
using namespace Wt;
using SpecUtils::Measurement;


namespace
{
  /** Shared helper: loads spectrum file and determines foreground/background measurements.

   @returns True on success, with foreground (and optionally background) populated.
   */
  bool load_spec_data( const char * const filename,
                       std::shared_ptr<SpecMeas> &spec,
                       std::shared_ptr<const SpecUtils::Measurement> &foreground,
                       std::shared_ptr<const SpecUtils::Measurement> &background,
                       std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &peaks )
  {
    spec = std::make_shared<SpecMeas>();

    if( !spec->load_file( filename, SpecUtils::ParserType::Auto, filename ) )
    {
      fprintf( stderr, "Failed to parse '%s' as a spectrum file.\n", filename );
      return false;
    }

    if( !spec->num_measurements() )
    {
      fprintf( stderr, "No measurements in '%s'.\n", filename );
      return false;
    }

    const std::set<int> &samplenums = spec->sample_numbers();

    std::set<int> foreground_sample_numbers, background_sample_numbers;

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
        SpecUtils::SourceType src_type = SpecUtils::SourceType::Unknown;
        for( const int detnum : spec->detector_numbers() )
        {
          const auto m = spec->measurement( sample, detnum );
          if( m && m->num_gamma_channels() )
          {
            has_meas = true;
            src_type = m->source_type();
          }
        }

        if( !has_meas )
          continue;

        switch( src_type )
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
        }//switch( src_type )
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
      return false;
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
      return false;
    }

    peaks = spec->peaks( foreground_sample_numbers );

    return true;
  }//load_spec_data(...)
}//anonymous namespace


CGImageRef render_spec_file_to_cgimage( const char * const filename,
                                        const float width, const float height,
                                        const enum SpecPreviewType type,
                                        const char * const logo_path )
{
  std::shared_ptr<SpecMeas> spec;
  std::shared_ptr<const SpecUtils::Measurement> foreground, background;
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks;

  if( !load_spec_data( filename, spec, foreground, background, peaks ) )
    return nullptr;

  try
  {
    CGPaintDevice device( width, height );

    // Fill with white background
    CGContextRef ctx = device.cgContext();
    CGContextSetRGBFillColor( ctx, 1.0, 1.0, 1.0, 1.0 );
    CGContextFillRect( ctx, CGRectMake( 0, 0, width, height ) );

    const bool is_thumbnail = (type == SpectrumThumbnail);

    if( !is_thumbnail && spec->passthrough() )
    {
      // Preview of passthrough: spectrum on top 2/3, time series on bottom 1/3
      const float spec_height = 2.0f * height / 3.0f;
      const float time_height = height - spec_height;

      SpectrumChartImage::SpectrumOptions spec_options;
      spec_options.foreground = foreground;
      spec_options.background = background;
      spec_options.peaks = peaks;
      spec_options.title = foreground->title();

      CGPaintDevice specDevice( device, 0, 0, width, spec_height );
      if( !SpectrumChartImage::paintSpectrum( specDevice, spec_options ) )
        throw std::runtime_error( "Failed in paintSpectrum for passthrough" );

      SpectrumChartImage::TimeSeriesOptions time_options;
      time_options.meas = spec;
      time_options.highlight_occupied = true;

      CGPaintDevice timeDevice( device, 0, spec_height, width, time_height );
      if( !SpectrumChartImage::paintTimeSeries( timeDevice, time_options ) )
        throw std::runtime_error( "Failed in paintTimeSeries" );
    }else
    {
      SpectrumChartImage::SpectrumOptions spec_options;
      spec_options.foreground = foreground;
      spec_options.background = background;
      spec_options.peaks = peaks;
      spec_options.title = foreground->title();
      spec_options.thumbnail = is_thumbnail;
      spec_options.logo_path = (logo_path ? logo_path : "");
      if( !SpectrumChartImage::paintSpectrum( device, spec_options ) )
        throw std::runtime_error( "Failed in paintSpectrum" );
    }

    return device.createCGImage();
  }catch( std::exception &e )
  {
    fprintf( stderr, "Failed rendering CGImage: %s\n", e.what() );
    return nullptr;
  }
}//render_spec_file_to_cgimage(...)
