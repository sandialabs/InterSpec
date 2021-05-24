#ifndef AssignPeaksToRefLine_h
#define AssignPeaksToRefLine_h
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

#include <deque>
#include <memory>
#include <vector>
#include <string>

#include <boost/function.hpp>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

//Forward declarations
class PeakDef;
class PeakModel;
class InterSpec;
struct ColorTheme;
namespace SpecUtils{ class Measurement; }
struct ReferenceLineInfo;
class ReferencePhotopeakDisplay;

namespace Wt{
  class WSvgImage;
};

/** Functions in this header/source are kinda go betweens of the GUI and the
 numerical code (although of course, seperation is never as clean as one would
 like).
 */
namespace PeakSearchGuiUtils
{
  
/** Renders the selcted energy range of the measurement to a SVG image.
 @param measurement What to plot.  Must be non-null pointer.
 @param peaks Peaks to include on the histogram.  May be empty.
 @param reflines Reference lines to show on the chart. May be empty
 @param lower_energy The lower energy to plot.  If equal to upper_energy than
        full energy range will be plotted.
 @param upper_energy The upper energy to plot.
 @param width_px Width of the SVG in pixels; used to decide re-binning and such.
 @param height_px Height of the SVG in pixels.
 @param theme The color theme of plot.  May be nullptr.
 @param compact Whether to plot the chart to maximize plot area.
 */
std::shared_ptr<Wt::WSvgImage> renderChartToSvg( std::shared_ptr<const SpecUtils::Measurement> measurement,
                                                std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef> > > peaks,
                                                const std::vector<std::shared_ptr<const ReferenceLineInfo>> &reflines,
                                              double lower_energy, double upper_energy,
                                              const int width_px, const int height_px,
                                              std::shared_ptr<const ColorTheme> theme,
                                                const bool compact );
  
/** Performs the automated search for peaks - setting the results to the GUI. */
void automated_search_for_peaks( InterSpec *interspec, const bool keep_old_peaks );


/** For all peaks currently fit for, that do not already have a
   nuclide/xray/reaction assigned, it attempts to assign them a source
   based on currently showing reference lines.
*/
void assign_peak_nuclides_from_reference_lines( InterSpec *viewer );
  
/** Attempts to assign a nuclide/xray/reaction from the currently displayed
   reference gamma lines.  This is called for instance when you are showing
   the reference lines for a nuclide and double click to create a peak.
   This function may change assignments for neighboring peaks in the PeakModel
   if it finds a better match given that there is a new peak in the region.
  */
void assign_nuclide_from_reference_lines( PeakDef &peak,
                                          PeakModel *peakModel,
                                          const std::shared_ptr<const SpecUtils::Measurement> &data,
                                          const ReferencePhotopeakDisplay *refLineDisp,
                                          const bool colorPeaksBasedOnReferenceLines,
                                          const bool showingEscapePeakFeature
                                         );
 
  
  
//searchForPeaksWorker(...): performs the actual search for peaks of the
//  entire spectrum.  Once search is done, results are placed into
//  'resultpeaks'  (which presumamble 'callback' also has a copy of) and
//  'callback' is posted to the WServer thread pool to be executed inside the
//  main event loop thread, and is done regardless of succesfulness of the
//  peak search.  'callback' is intended to be a bound function call to
//  setPeaksFromSearch(...) or setHintPeaks(...).
void search_for_peaks_worker( std::weak_ptr<const SpecUtils::Measurement> weak_data,
                             std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > existingPeaks,
                               const std::vector<ReferenceLineInfo> displayed,
                               const bool setColorFromRefLine,
                               std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks,
                               boost::function<void(void)> callback,
                               const std::string sessionID,
                               const bool singleThread);
  
/** Assigns peak nuclides/xrays/reactions from the reference photopeak lines by
   modifying the peaks passed in.  If a peak already has a nuclide set, it wont
   be changed unless a new peak is a better candidate for it, and there is
   another ref line that explains it.
 */
void assign_srcs_from_ref_lines( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                 std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > peaks,
                                const std::vector<ReferenceLineInfo> &displayed,
                                 const bool setColor,
                                const bool showingEscapePeakFeature );
  
/** Refits the peaks from a right-click refit request.
 Assumes you are in the Wt app primary thread.
 */
void refit_peaks_from_right_click( InterSpec * const interspec, const double rightClickEnergy );

/** Changes the continuum type and causes a refit of ROI.
 
 @param interspec The InterSpec instance to work with - it is assumed this function is being called from that apps primary thread.
 @param rightClickEnergy Energy used to identify a peak, that will in turn identify the ROI.
 @param continuum_type the #PeakContinuum::OffsetType
 */
void change_continuum_type_from_right_click( InterSpec * const interspec,
                                            const double rightClickEnergy,
                                            const int continuum_type );


/** Enum to tell #fit_template_peaks where the candidate peaks to fit are
   comming from.  This info will be propagated through to the GUI and influence
   what the user can do and see.
 */
enum class PeakTemplateFitSrc
{
  PreviousSpectrum,
  CsvFile
};//enum class PeakTemplateFitSrc
  

/** Fits the template peaks to the passed in data, and then will open a dialog
  to allow the user to select what peaks they want to keep.
 
 @param interspec Application instance this work is being performed for.
 @param data The gamma spectrum to fit the peaks to.  Must be valid spectrum.
 @param template_peaks Candidate peaks to fit.  The peak continuum and amplitude
        will be estimated based off of the data before fitting incase they are
        way off.  The continuum order, source nuclide, color, ROI extent, etc.,
        will all be directly preserved.  The FWHM will be fit for within some
        small range (so input peaks should be reasonably close), the amplitude
        and continuum coefficients will be freely fit for.
 @param original_peaks The original peaks in the spectrum.  These wont be
        altered.  These are the peaks in the spectrum before fitting for the
        template peaks.
 @param fitsrc The source of the template peaks.  Either from a previous
        spectrum or a peak CSV file.  Effects how results are displayed.
 @param sessionid The WApplication session ID to post to, to display results.
 
 As of 20191028 not particularly well tested.
 */
void fit_template_peaks( InterSpec *interspec,
                         std::shared_ptr<const SpecUtils::Measurement> data,
                         std::vector<PeakDef> template_peaks,
                         std::vector<PeakDef> original_peaks,
                         const PeakTemplateFitSrc fitsrc,
                         const std::string sessionid );
  
}//namespace PeakSearchGuiUtils
#endif //AssignPeaksToRefLine_h
