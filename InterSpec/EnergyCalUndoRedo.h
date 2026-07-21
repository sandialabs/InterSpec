#ifndef EnergyCalUndoRedo_h
#define EnergyCalUndoRedo_h
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
#include <tuple>
#include <memory>
#include <vector>

class PeakDef;
class SpecMeas;

namespace SpecUtils
{
  class Measurement;
  struct EnergyCalibration;
}//namespace SpecUtils


namespace EnergyCalImp
{

/** For undo/redo, we store the mappings, for each measurement, from old to new energy
 calibration.  Using weak ptrs for Measurement, although I dont think it is strickly necassary.
 TODO: store energy calibration a little more compactly; could store just the coefficients, and
       not the channel lower-energies
 */
typedef std::vector<std::tuple<std::weak_ptr<const SpecUtils::Measurement>,
          std::shared_ptr<const SpecUtils::EnergyCalibration>,
          std::shared_ptr<const SpecUtils::EnergyCalibration>>>
        meas_old_new_cal_t;

/** For undo/redo, keep track of new and old peaks.
 TODO: translate peaks on-the-fly, to reduce memory use, but will need to track pre-post energy
       cal to enable this
 */
typedef std::vector<std::tuple<std::set<int>,
        std::deque<std::shared_ptr<const PeakDef>>,
        std::deque<std::shared_ptr<const PeakDef>>>>
      meas_old_new_peaks_t;


/** Aggregates the energy calibration (and associated peak) changes of a single user operation
 into one undo/redo step.

 We may make multiple calls to #EnergyCalTool::setEnergyCal and/or #EnergyCalTool::addDeviationPair,
 or other functions, for a single user-instigated change; we want to combine multiple of these
 calls into a single user undo/redo operation, so thread-local storage is used to aggregate all
 the changes (you may construct as many, possibly nested, sentries during the operation as is
 convenient), and the destructor of the last live sentry actually inserts the undo/redo step.

 Changes are recorded per affected SpecMeas (so restoration targets the right measurements), but
 the undo/redo step itself is associated with the foreground, as all InterSpec undo/redo is;
 changes to files not currently displayed as foreground/background/secondary do not get
 undone/redone.
 */
struct EnergyCalUndoRedoSentry
{
  EnergyCalUndoRedoSentry();
  ~EnergyCalUndoRedoSentry();

  meas_old_new_cal_t &cal_info( const std::shared_ptr<SpecMeas> &meas );
  meas_old_new_peaks_t &peak_info( const std::shared_ptr<SpecMeas> &meas );
  meas_old_new_peaks_t &hint_peak_info( const std::shared_ptr<SpecMeas> &meas );
};//struct EnergyCalUndoRedoSentry

}//namespace EnergyCalImp

#endif //EnergyCalUndoRedo_h
