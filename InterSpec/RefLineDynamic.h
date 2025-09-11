#ifndef RefLineDynamic_h
#define RefLineDynamic_h
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
#include <atomic>
#include <memory>
#include <string>
#include <vector>


#include <Wt/WObject>

namespace SpecUtils{ enum class SpectrumType : int; }

class SpecMeas;
class InterSpec;
struct AlwaysSrcs;
struct ColorTheme;
struct ReferenceLineInfo;
struct ExternalRidResults;
class D3SpectrumDisplayDiv;
class DetectorPeakResponse;

/**
TODO items:
 - Extend the energy lines to the full-height of the window - maybe either use a dotted line, or one with some alpha, or make the solid line, but a dot at the original height
 - For the random sum peaks, also show likely nuclides
 - When multiple nuclides are candidates, need to allow switching between them, or indicate things in some way to either cycle through them, and/or let the user
 - Each "catagory" of nuclides should get a different color, controllable by the color theme

 Current shortcommings:
 - Cascade sums not accounted for
 - The filtering the number of lines for each set could def use improvement
 - Lines not updated when there is an energy calibration change (the auto-serch peak should also be updated - and any being searched for need to be cancelled and restarted)
 - Should be able to handle RID output results like DU, Uranium, HEU, Neutrons, "Sr-85/Kr85", Annihilation, etc
 - Having a peak with a source associated with it should _really_ increase the weight (see HPGe Lu177 spectrum for example why).
 - Should pick the most important line, get the profile weight for the source, and adjust the source weight a bit using this
 */
class RefLineDynamic : public Wt::WObject
{
public:
  RefLineDynamic( D3SpectrumDisplayDiv *chart, InterSpec *viewer );
  virtual ~RefLineDynamic();
  
  bool successfully_initialized() const;
  
  void setActive( bool active );
  bool isActive() const;
  
  std::shared_ptr<std::vector<std::pair<double,ReferenceLineInfo>>> current_lines() const;
  
  /** Flags to indicate what actions need to be taken during pushUpdates. */
  enum DynamicRefLineRenderFlags
  {
    UpdateLines = 0x01
  };//enum DynamicRefLineRenderFlags
  
  /** Function to _start_ computing updated lines, in a background thread, and then when done will push to `m_chart`. */
  void startPushUpdates();
  
  enum class RefLineSrc : int
  {
    UserPeakLines,
    ExternalRid,
    OnboardRid,
    AlwaysShowing,
    UserPeakAssociatedNucLines,
    CharacteristicLine,
    EscapeLines,
    RandomSumLines,
  };//enum class RefLineSrc : int
  
  /** Attempts to filter the lines in `ref_lines` to just the most-likely relevant lines - we dont want like 500 lines the client has to deal with for a 30 second
   spectrum
   */
  static void filterLines( ReferenceLineInfo &ref_lines,
                          const RefLineSrc src,
                          const std::shared_ptr<const SpecUtils::Measurement> &meas,
                          const std::shared_ptr<const DetectorPeakResponse> &detector );
  
protected:
  void start_init_always_sources();
  
  void autoSearchPeaksSet( const SpecUtils::SpectrumType spectrum );
  void spectrumChanged( const SpecUtils::SpectrumType spec_type,
                       const std::shared_ptr<SpecMeas> &measurement,
                       const std::set<int> &sample_numbers,
                       const std::vector<std::string> &detectors );
  void autoRidResultsRecieved( const std::shared_ptr<const ExternalRidResults> &results );
  void colorThemeChanged( const std::shared_ptr<const ColorTheme> &theme );

  void startUpdateLines();
  void finishUpdateLines( const std::shared_ptr<std::vector<std::pair<double,ReferenceLineInfo>>> &ref_lines,
                         const std::shared_ptr<std::string> &js_fwhm_fcn,
                         const size_t calc_num );
  
  /** Helper method to assign color to ReferenceLineInfo, applying fallback logic. */
  void assignColorToInput( ReferenceLineInfo &input ) const;

  void peaksAdded();
  void peaksRemoved();
  void peakModified();

  InterSpec *m_interspec;
  D3SpectrumDisplayDiv *m_chart;
  
  bool m_active;
  
  /** Will be set to false until the "always" sources have been attempted to initialize. */
  bool m_has_inited;
  
  /** Will be empty if the "always" sources were succesfully initialized, and non-empty otherwise. */
  std::string m_init_error_msg;
  
  /** A pointer to information about source that should always be included in the dynamic reference lines.
   Note: this is a shared pointer to keep alive when doing work in a background thread (launched doing a lambda defined in
   `startUpdateLines()`), and the session is ended while the work is being done.
   */
  std::shared_ptr<const AlwaysSrcs> m_always_srcs;
  
  std::shared_ptr<const ExternalRidResults> m_external_rid_results;
  
  /** Flags to track what updates need to be made during render. */
  Wt::WFlags<DynamicRefLineRenderFlags> m_renderFlags;
  
  std::shared_ptr<std::atomic<size_t>> m_current_calc_num;
  
  std::shared_ptr<std::vector<std::pair<double,ReferenceLineInfo>>> m_current_ref_lines;
};//class RefLineDynamic

#endif // RefLineDynamic_h
