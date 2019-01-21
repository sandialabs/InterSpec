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
class Measurement;
class ReferencePhotopeakDisplay;


/** Functions in this header/source are kinda go betweens of the GUI and the
 numerical code (although of course, seperation is never as clean as one would
 like).
 */
namespace PeakSearchGuiUtils
{
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
                                          const std::shared_ptr<const Measurement> &data,
                                          const ReferencePhotopeakDisplay *refLineDisp,
                                          const bool colorPeaksBasedOnReferenceLines );
 
  
  
//searchForPeaksWorker(...): performs the actual search for peaks of the
//  entire spectrum.  Once search is done, results are placed into
//  'resultpeaks'  (which presumamble 'callback' also has a copy of) and
//  'callback' is posted to the WServer thread pool to be executed inside the
//  main event loop thread, and is done regardless of succesfulness of the
//  peak search.  'callback' is intended to be a bound function call to
//  setPeaksFromSearch(...) or setHintPeaks(...).
void search_for_peaks_worker( std::weak_ptr<const Measurement> weak_data,
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
void assign_srcs_from_ref_lines( const std::shared_ptr<const Measurement> &data,
                                 std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > peaks,
                                const std::vector<ReferenceLineInfo> &displayed,
                                 const bool setColor );
  
}//namespace PeakSearchGuiUtils
#endif //AssignPeaksToRefLine_h
