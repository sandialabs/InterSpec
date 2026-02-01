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

#ifndef FitPeaksForNuclidesGui_h
#define FitPeaksForNuclidesGui_h

#include "InterSpec_config.h"

class SimpleDialog;

/** Functions in this namespace handle the GUI workflow for fitting peaks
 for nuclides shown in the ReferencePhotopeakDisplay.

 These functions must be called from the Wt application thread.
 */
namespace FitPeaksForNuclidesGui
{
  /** Creates and shows the advanced options dialog for fitting peaks.

   Returns the created dialog so the caller can track it.
   The caller is responsible for managing the dialog lifecycle.

   Must be called from the Wt application thread.

   @return Pointer to the created SimpleDialog, or nullptr if creation failed.
   */
  SimpleDialog *showAdvancedDialog();

  /** Starts the peak fitting workflow for currently displayed nuclides.

   This function performs a 3-stage async operation:
   - Stage A (GUI thread): Peak detection using AnalystChecks
   - Stage B (Background thread): Calls FitPeaksForNuclides::fit_peaks_for_nuclides()
   - Stage C (GUI thread): Results dialog with preview, color assignment, accept/cancel

   @param from_advanced_dialog If true, indicates call is from the advanced dialog
          (currently unused, but reserved for future options).

   Must be called from the Wt application thread.
   */
  void startFitSources( const bool from_advanced_dialog );

}//namespace FitPeaksForNuclidesGui

#endif //FitPeaksForNuclidesGui_h
