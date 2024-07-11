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
#include <tuple>
#include <memory>
#include <vector>
#include <string>

#include <boost/function.hpp>


//Forward declarations
class PeakDef;
class PeakModel;
class InterSpec;
struct ColorTheme;
class SpectrumChart;
struct ReferenceLineInfo;
class DetectorPeakResponse;
class ReferencePhotopeakDisplay;

namespace SpecUtils{
  class Measurement;
}

namespace Wt{
  class WSvgImage;
};

/** Functions in this header/source are kinda go between of the GUI and the
 numerical code (although of course, separation is never as clean as one would
 like).
 */
namespace PeakSearchGuiUtils
{
  
/** Renders the selected energy range of the measurement to a SVG image.
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

/** Same as `renderChartToSvg(...)`, but returns  a `SpectrumChart` that may be painted in a <canvas/> element,
 or whatever.
 (note: Wt 3.7.1 seems to have a issue in rendering to SVG, where sometimes there will be lines and artifacts, so this is a workaround)
 */
SpectrumChart *createFixedSpectrumDisplay( std::shared_ptr<const SpecUtils::Measurement> inmeas,
                            std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef> > > peaks,
                            const std::vector<std::shared_ptr<const ReferenceLineInfo>> &displayed,
                            double lowx, double upperx,
                            const int width, const int height,
                            std::shared_ptr<const ColorTheme> theme );
  
/** Function that is called when the user double-left-clicks on the spectrum.
 */
void fit_peak_from_double_click( InterSpec *interspec,
                                const double energy_clicked,
                                const double pixel_per_keV,
                                std::shared_ptr<const DetectorPeakResponse> det );
  
  
/** Performs the automated search for peaks - setting the results to the GUI. */
void automated_search_for_peaks( InterSpec *interspec, const bool keep_old_peaks );

/** Uses the currently displayed foreground to estimate the expected FWHM for the wanted energy.
 
 Value returned, is the first successful one of:
 - Interpolate between user-fit peaks
 - Use FWHM function of current detector efficiency function
 - Interpolate between auto-search peaks
 - sqrt(energy) interpolate from a single peak
 - Use a generic HPGe or NaI detector
 - return 0.0f - shouldnt happen.
 */
float estimate_FWHM_of_foreground( const float energy );

/** For all peaks currently fit for, this function attempts to assign them a source
 based on currently showing reference lines.
 
 @param only_peaks_with_no_src If true, only consider peaks that do not currently have an assigned nuclide/x-ray/reaction.
 @param only_current_ref_lines If true, then only use the current (e.g., top-most) reference lines, and not the "persisted" lines.
*/
void assign_peak_nuclides_from_reference_lines( InterSpec *viewer,
                                               const bool only_peaks_with_no_src,
                                               const bool only_current_ref_lines );
  
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
                             std::shared_ptr<const DetectorPeakResponse> drf,
                             std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > existingPeaks,
                               const std::vector<ReferenceLineInfo> displayed,
                               const bool setColorFromRefLine,
                               std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks,
                               boost::function<void(void)> callback,
                               const std::string sessionID,
                               const bool singleThread );
  
/** Assigns peak nuclides/xrays/reactions from the reference photopeak lines by
   modifying the peaks passed in.  
 
 @param only_peaks_with_no_src If false, If a peak already has a nuclide set, it wont be changed unless a new peak is a
        better candidate for it, and there is another ref line that explains it.  If true, peaks with any nuclide/x-ray-reaction assigned
        will be skipped.
 */
void assign_srcs_from_ref_lines( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                 std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > peaks,
                                const std::vector<ReferenceLineInfo> &displayed,
                                 const bool setColor,
                                const bool showingEscapePeakFeature,
                                const bool only_peaks_with_no_src );
  
/** Refits the peaks from a right-click refit request.
 Assumes you are in the Wt app primary thread.
 */
void refit_peaks_from_right_click( InterSpec * const interspec, const double rightClickEnergy );

/** Sets all peaks in the ROI pointed to by the passed in energy, to the FWHM specified by the DRF, then refits the peaks.
 
 Assumes you are in the Wt app primary thread.
 */
void refit_peaks_with_drf_fwhm( InterSpec * const interspec, const double rightClickEnergy );
  
/** Returns the energy of its assigned nuclides gamma line, or peak doesnt have assigned gamma line, returns the
 nearest showing Reference Photopeak line. Returns a negative value if neither can be found
 
 */
float reference_line_energy_near_peak( InterSpec * const interspec, 
                                      const PeakDef &peak );
  
/** Returns the reference line info, and index for the specific reference line
 
 @param only_nuclide If true, then only lines with a parent nuclide, who are a gamma or xray, will be considered
 
 Example usage:
 ```
 auto refline = reference_line_near_peak( interspec, peak, false );
 if( refline.first )
 {
   ReferenceLineInfo::RefLine &line = refline.first->m_ref_lines[refline.second];
   ...
 }else
 {
   // No line found
 }
 ```
 */
std::pair<std::unique_ptr<ReferenceLineInfo>,int> reference_line_near_peak( 
                                                                  InterSpec * const interspec,
                                                                  const PeakDef &peak,
                                                                  const bool only_nuclide );
  
/** Returns the reference line nuclide and energy near the given energy, that the user probably intends.
   
  Intended to be called with the energy of a users click - takes into account detectors resolution and nuclide yields.
   
  @returns The nuclide, its age, and energy of the reference line.  If not near a line, or no reference line, or a non-nuc reference
          line is showing, will return {nullptr, 0.0, 0.0f}.
   
  \sa PeakSearchGuiUtils::reference_line_energy_near_peak
*/
std::tuple<const SandiaDecay::Nuclide *, double /* age */, float /* energy */>
  nuclide_reference_line_near( InterSpec *viewer, const float energy );
  
  
/** Set the peak nearest `rightClickEnergy` to preferably its assigned gamma energy, or if none assigned, to
 nearest showing reference photopeak energy, and then refits peak.
 */
void refit_peak_with_photopeak_mean( InterSpec * const interspec, const double rightClickEnergy );
  
/** Changes the continuum type and causes a refit of ROI.
 
 @param interspec The InterSpec instance to work with - it is assumed this function is being called from that apps primary thread.
 @param rightClickEnergy Energy used to identify a peak, that will in turn identify the ROI.
 @param continuum_type the #PeakContinuum::OffsetType
 */
void change_continuum_type_from_right_click( InterSpec * const interspec,
                                            const double rightClickEnergy,
                                            const int continuum_type );

/** Changes the continuum type of a ROI, and does a refit.
 @param interspec The InterSpec instance to work with - it is assumed this function is being called from that apps primary thread.
 @param rightClickEnergy Energy used to identify a peak, that will in turn identify the ROI.
 @param skew_type the #PeakDef::SkewType
*/
void change_skew_type_from_right_click( InterSpec * const interspec,
                                       const double rightClickEnergy,
                                       const int skew_type );

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

void prepare_and_add_gadras_peaks(
                        std::shared_ptr<const SpecUtils::Measurement> data,
                        std::vector<PeakDef> gadras_peaks,
                        std::vector<PeakDef> original_peaks,
                        const std::string sessionid );



}//namespace PeakSearchGuiUtils
#endif //AssignPeaksToRefLine_h
