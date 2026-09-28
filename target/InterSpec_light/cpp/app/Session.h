#ifndef LightSession_h
#define LightSession_h
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

#include <set>
#include <map>
#include <deque>
#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"

#include "RefLib.h"

class PeakDef;
struct PeakContinuum;
struct PeakFitDetPrefs;
namespace SpecUtils
{
  class SpecFile;
  class Measurement;
  class EnergyCalibration;
}

/** All state of the light app: loaded files, the foreground/background/secondary display slots,
 peaks, and reference lines.  Every public entry is `call(method, params)`, which returns a JSON
 object holding only the parts of the GUI state that changed (see `Part`), or throws.
 */
class Session
{
public:
  typedef std::deque<std::shared_ptr<const PeakDef>> PeakDeque;

  enum SpecType : int { Foreground = 0, Background = 1, Secondary = 2, NumSpecType = 3 };

  Session();
  ~Session();

  nlohmann::json call( const std::string &method, const nlohmann::json &params );

private:
  struct LoadedFile
  {
    int id = -1;
    std::string name;
    std::shared_ptr<SpecUtils::SpecFile> spec;
    /** Energy calibration of each Measurement as loaded, for "revert". */
    std::map<const SpecUtils::Measurement *,std::shared_ptr<const SpecUtils::EnergyCalibration>> orig_cals;
    /** Peaks, keyed by the foreground sample numbers they were fit on. */
    std::map<std::set<int>,PeakDeque> peaks;
    std::shared_ptr<PeakFitDetPrefs> prefs;
  };

  struct Slot
  {
    std::shared_ptr<LoadedFile> file;
    std::set<int> samples;
    /** Sum of the displayed samples and detectors; may be null (e.g., neutron-only selection). */
    std::shared_ptr<const SpecUtils::Measurement> summed;
  };

  /** Bit flags of the state parts to include in a response. */
  enum Part : unsigned
  {
    PartFiles       = 0x01,
    PartSpectra     = 0x02,
    PartTimeChart   = 0x04,
    PartPeaks       = 0x08,
    PartEnergyCal   = 0x10,
    PartDetectors   = 0x20,
    PartAll         = 0xFF
  };

  nlohmann::json state( const unsigned parts, const bool resetDomain = false ) const;

  // Files, samples, and detectors (Session.cpp)
  nlohmann::json loadFile( const nlohmann::json &p );
  nlohmann::json unload( const nlohmann::json &p );
  nlohmann::json setSamples( const nlohmann::json &p );
  nlohmann::json stepSample( const nlohmann::json &p );
  nlohmann::json timeDrag( const nlohmann::json &p );
  nlohmann::json setDetectors( const nlohmann::json &p );
  nlohmann::json exportFile( const nlohmann::json &p );

  void setSlot( const SpecType type, std::shared_ptr<LoadedFile> file, std::set<int> samples );
  void updateSummed( const SpecType type );
  void updateAllSummed();
  std::vector<std::string> displayedDetectors( const SpecType type ) const;
  std::set<int> defaultSamplesForType( const SpecType type ) const;

  nlohmann::json filesJson() const;
  nlohmann::json spectraJson() const;
  nlohmann::json timeChartJson() const;
  nlohmann::json timeHighlightsJson() const;
  nlohmann::json detectorsJson() const;

  // Peaks (SessionPeaks.cpp)
  nlohmann::json fitPeakAt( const nlohmann::json &p );
  nlohmann::json deletePeakAt( const nlohmann::json &p );
  nlohmann::json erasePeaks( const nlohmann::json &p );
  nlohmann::json roiDrag( const nlohmann::json &p );
  nlohmann::json fitRoiDrag( const nlohmann::json &p );
  nlohmann::json peakInfoAt( const nlohmann::json &p );
  nlohmann::json setContinuumType( const nlohmann::json &p );
  nlohmann::json setSkewType( const nlohmann::json &p );
  nlohmann::json refitRoi( const nlohmann::json &p );
  nlohmann::json setPeakProperty( const nlohmann::json &p );
  nlohmann::json clearPeaks( const nlohmann::json &p );
  nlohmann::json exportPeakCsv( const nlohmann::json &p );
  nlohmann::json searchPeaks( const nlohmann::json &p );
  nlohmann::json setRefLibrary( const nlohmann::json &p );
  nlohmann::json setShownRefLines( const nlohmann::json &p );

  PeakDeque &foregroundPeaks();
  const PeakDeque *foregroundPeaksConst() const;
  std::shared_ptr<const PeakDef> peakContaining( const double energy ) const;
  int effectiveDetType() const;
  bool isHighRes() const;
  void assignSourceToNewPeak( PeakDef &peak, const std::string &refLineParent, PeakDeque &peaks );
  std::string peaksToChartJson( const std::vector<std::shared_ptr<const PeakDef>> &peaks ) const;
  nlohmann::json peaksJson() const;
  nlohmann::json peakListJson() const;
  void resetDragCaches();

  // Energy calibration (SessionEnergyCal.cpp)
  nlohmann::json setEnergyCal( const nlohmann::json &p );
  nlohmann::json fitEnergyCal( const nlohmann::json &p );
  nlohmann::json revertEnergyCal( const nlohmann::json &p );
  nlohmann::json convertCalCoefs( const nlohmann::json &p );
  nlohmann::json energyCalJson() const;
  std::shared_ptr<const SpecUtils::EnergyCalibration> displayedEnergyCal() const;
  /** Applies a change of the displayed calibration to every gamma Measurement of the foreground
   file (all samples, displayed detectors), and moves all of that file's peaks to match. */
  void applyCalChange( const std::shared_ptr<const SpecUtils::EnergyCalibration> &disp_prev,
                       const std::shared_ptr<const SpecUtils::EnergyCalibration> &disp_new );

  int m_next_file_id;
  Slot m_slots[NumSpecType];
  /** Detectors shown for the foreground file (and files with the same detector names). */
  std::vector<std::string> m_shown_detectors;

  RefLib::Library m_ref_lib;
  std::vector<RefLib::Shown> m_shown_ref;

  // Reuse of fits while the user drags a ROI edge, or drags out a new ROI.
  std::shared_ptr<const PeakContinuum> m_drag_continuum;
  std::vector<std::shared_ptr<const PeakDef>> m_drag_last_peaks;
  std::chrono::steady_clock::time_point m_drag_time;
  std::vector<std::shared_ptr<const PeakDef>> m_create_roi_peaks;
};//class Session

#endif //LightSession_h
