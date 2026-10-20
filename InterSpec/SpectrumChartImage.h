#ifndef SpectrumChartImage_h
#define SpectrumChartImage_h
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
#include <deque>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <optional>
#include <functional>

// Forward declarations
class PeakDef;
class SpecMeas;
struct ColorTheme;
struct ReferenceLineInfo;

namespace SpecUtils
{
  class Measurement;
}

namespace Wt
{
  class WPaintDevice;
}


/** Draws spectrum and time-history charts with `Wt::WPainter`, onto any `Wt::WPaintDevice`, with no
 browser or `WApplication` needed: for the macOS Quick Look extension (CoreGraphics), and for chart
 images where no browser is attached (the stdio MCP server, plutovg).
 */
namespace SpectrumChartImage
{
  /** What to draw for a spectrum chart. */
  struct SpectrumOptions
  {
    std::shared_ptr<const SpecUtils::Measurement> foreground;
    std::shared_ptr<const SpecUtils::Measurement> background;
    std::shared_ptr<const SpecUtils::Measurement> secondary;

    /** Factors the background and secondary spectra are multiplied by when drawn (e.g.,
     `InterSpec::displayScaleFactor(...)`); if not positive, the ratio of live times to the foreground.
     */
    double background_scale = -1.0;
    double secondary_scale = -1.0;

    /** Peaks of the foreground. */
    std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks;

    /** Reference lines to draw; the first is drawn as the currently showing one. */
    std::vector<std::shared_ptr<const ReferenceLineInfo>> reference_lines;

    /** The energy range shown, in keV; the whole spectrum if not set. */
    std::optional<std::pair<double,double>> energy_range;

    bool log_y = true;

    /** Draw the foreground minus the scaled background. */
    bool background_subtract = false;

    /** Title drawn above the chart, if not empty (and the image is wide enough). */
    std::string title;

    /** For small thumbnails: axis lines without labels, a thin spectrum line, and (if
     `logo_path` is not empty) the logo in the upper right.
     */
    bool thumbnail = false;
    std::string logo_path;

    /** Colors; InterSpec's defaults if null. */
    std::shared_ptr<const ColorTheme> theme;
  };//struct SpectrumOptions

  /** Paints a spectrum chart filling `device`.  Returns false if there is no foreground. */
  bool paintSpectrum( Wt::WPaintDevice &device, const SpectrumOptions &options );


  /** What to draw for a time-history (gross count rate versus time) chart. */
  struct TimeSeriesOptions
  {
    std::shared_ptr<const SpecMeas> meas;

    /** Detectors to sum; all of them if empty. */
    std::vector<std::string> detectors;

    /** Sample numbers to highlight, as the foreground, background, and secondary spectra. */
    std::set<int> foreground_samples;
    std::set<int> background_samples;
    std::set<int> secondary_samples;

    /** If true, highlights the occupied samples as foreground, instead of `foreground_samples`. */
    bool highlight_occupied = false;

    std::shared_ptr<const ColorTheme> theme;
  };//struct TimeSeriesOptions

  /** Paints a time-history chart filling `device`.  Returns false if there are fewer than two
   samples.
   */
  bool paintTimeSeries( Wt::WPaintDevice &device, const TimeSeriesOptions &options );


  /** Renders `paint` (which paints onto the device it is given, of the given size) to an encoded
   raster image of `format` ("png" or "jpeg").  Set by programs that can (the stdio MCP server).
   */
  using RasterRenderer = std::function<std::string( int width, int height, const std::string &format,
                                                    const std::function<void( Wt::WPaintDevice & )> &paint )>;

  void setRasterRenderer( RasterRenderer renderer );

  /** If a `RasterRenderer` has been set. */
  bool canRenderRaster();

  /** Renders `paint` to an image of `format`: "svg" (with `Wt::WSvgImage`), or "png"/"jpeg" (with the
   `RasterRenderer`).  Returns the encoded image.
   @throws std::exception if the format can not be rendered.
   */
  std::string renderImage( int width, int height, const std::string &format,
                           const std::function<void( Wt::WPaintDevice & )> &paint );
}//namespace SpectrumChartImage

#endif //SpectrumChartImage_h
