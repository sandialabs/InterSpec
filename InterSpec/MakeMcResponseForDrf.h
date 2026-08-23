#ifndef MakeMcResponseForDrf_h
#define MakeMcResponseForDrf_h
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

#include <atomic>
#include <memory>
#include <string>
#include <vector>

#include <Wt/WContainerWidget.h>

#include "InterSpec/AuxWindow.h"

class DrfChart;
class InterSpec;
class DetectorPeakResponse;
class DetectorGeometryInput;

namespace Wt
{
  class WText;
  class WGroupBox;
  class WComboBox;
  class WLineEdit;
  class WPushButton;
  class WProgressBar;
}//namespace Wt

namespace ceelo
{
  class DetectorResponse;
  struct GroundingPoint;
  struct GeometryDescriptor;
}//namespace ceelo

/** Characterizes a detectors response over its geometry (CeeLo): the user
 enters/edits the detector geometry and picks how the response is built -
 full Monte-Carlo characterization, a quick MC-backbone efficiency transfer,
 or an instant (no-MC) transfer anchored on the DRF's measured efficiency
 curve; the parameterized response can then be attached to the current DRF
 ("Use Response") - all off-axis/near-field/uncertainty-aware efficiency
 queries then dispatch to it.

 If the seed DRF carries raw measured efficiency points (from "Make Detector
 Response"), a MC-generated response is grounded to them (k(E) + covariance);
 otherwise, if it has an efficiency curve, the curve is sampled as a
 lower-quality grounding anchor.  The measured-curve transfer method instead
 uses those same points/curve directly as the transfer anchor (the anchor IS
 the grounding), and rebuilds automatically as the geometry is edited.
 */
class MakeMcResponseForDrf : public Wt::WContainerWidget
{
public:
  /** How the response is built - the entries of the "Method" combo. */
  enum class Method : int
  {
    /** Full MC characterization over the selected profile (existing path). */
    FullMc = 0,

    /** EFFTRAN-style transfer from MC at only the on-axis energy backbone
     (plus optionally a few cos-theta anchors) - roughly a tenth the MC of a
     full characterization. */
    QuickMc = 1,

    /** EFFTRAN-style transfer anchored on the DRF's measured efficiency
     points/curve - deterministic, no MC, rebuilt automatically on geometry
     edits. */
    CurveTransfer = 2
  };//enum class Method

  MakeMcResponseForDrf( InterSpec *viewer,
                        std::shared_ptr<const DetectorPeakResponse> seed_drf,
                        std::shared_ptr<const ceelo::GeometryDescriptor> geometry = nullptr );

  virtual ~MakeMcResponseForDrf() override;

  /** Pre-populates the geometry form from a physical descriptor (e.g. an
   imported ANGLE detector model), overriding the DRF-derived guess.  Safe to
   call after construction to re-seed an already-open tool. */
  void setGeometryFromDescriptor( const ceelo::GeometryDescriptor &geometry );

  /** Emitted when a generated response becomes available/unavailable
   (enables the windows "Use Response" button).
   */
  Wt::Signal<bool> &validationChanged();

  /** Whether a generated response is ready to accept. */
  bool hasResult() const;

  /** The most recently generated (and grounded) response; may be nullptr. */
  std::shared_ptr<const ceelo::DetectorResponse> generatedResponse() const;

  /** Emitted after every successful generation (before the user accepts) -
   e.g. so the Make Detector Response tool can attach the response to the
   DRF it assembles.  Deliberately non-const: the receiver may still need to
   ground the response to its own measured points.
   */
  Wt::Signal<std::shared_ptr<ceelo::DetectorResponse>> &responseGenerated();

  /** "Use Response": clones the seed DRF (or current foreground DRF),
   attaches the generated response, makes it the sessions detector, and
   emits #updatedDrf.
   */
  void acceptResponse();

  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &updatedDrf();

  /** Grounding anchors from a DRF: raw measured points when present
   (preferred), else the legacy efficiency curve sampled at a reference
   distance (lower quality: curve lack-of-fit leaks into k(E)).
   Returns whether the points are curve-derived through `curve_derived`.
   Public/static so batch tools and the cascade-truth tests can ground a
   generated response to an existing DRF the same way this widget does.
   */
  static std::vector<ceelo::GroundingPoint> groundingPointsForDrf(
                            const std::shared_ptr<const DetectorPeakResponse> &drf,
                            const ceelo::GeometryDescriptor &geom,
                            bool &curve_derived );

  /** The currently selected build method. */
  Method selectedMethod() const;

protected:
  void handleGeometryChanged();
  void handleMethodChanged();
  void handlePrecisionChanged();
  void updateEstimate();

  /** Refreshes the measured-curve anchor description row (source label +
   reference-distance edit) from the seed DRF and current geometry. */
  void updateAnchorInfo();

  /** Refreshes the response-preview chart (per-angle efficiency curves) from
   the currently generated response, or hides it when there is none. */
  void updateResponseChart();

  void startGeneration();
  void cancelGeneration();

  /** Called (on the session thread) with progress from the worker. */
  void updateProgress( const double frac, const std::string &stage,
                       const int generation_id );
  void handleGenerationFinished( std::shared_ptr<ceelo::DetectorResponse> result,
                                 const std::string &errmsg,
                                 const int generation_id );

  /** Per-node MC FEP precision (base target) from the GUI selection. */
  double selectedPrecision() const;

  /** User-entered measured-curve anchor reference distance, in cm, or a
   non-positive value when empty/unparseable (meaning: use the automatic
   distance the curve is pinned at). */
  double refDistanceOverrideCm() const;

  /** True when the "Balanced" precision option is selected, i.e. use the
      graded relax_mild precision map (ceelo::PrecisionProfile::RelaxMild). */
  bool selectedRelaxMild() const;

  InterSpec *m_interspec;
  std::shared_ptr<const DetectorPeakResponse> m_seedDrf;

  DetectorGeometryInput *m_geometry;

  Wt::WComboBox *m_method;      //Full MC | Quick MC (transfer) | From measured curve
  Wt::WComboBox *m_profile;     //Far-field | General | Contact
  Wt::WComboBox *m_precision;   //Fast (1%) | Normal (0.3%) | Balanced (relax_mild) | Thorough (0.1%) | Custom
  Wt::WLineEdit *m_customPrecision;
  Wt::WComboBox *m_anchorAngles;    //On-axis only | 3 cos-theta anchors
  Wt::WText *m_anchorInfo;          //measured-curve anchor source description
  Wt::WLineEdit *m_refDistance;     //measured-curve anchor reference distance
  Wt::WText *m_estimate;

  //Rows shown/hidden per selected method:
  Wt::WContainerWidget *m_profileRow;
  Wt::WContainerWidget *m_precRow;
  Wt::WContainerWidget *m_anchorAnglesRow;
  Wt::WContainerWidget *m_anchorRow;

  Wt::WPushButton *m_generate;
  Wt::WPushButton *m_cancelBtn;
  Wt::WProgressBar *m_progress;
  Wt::WText *m_status;

  //Response-preview chart (per-angle efficiency curves), shown once a response
  // has been generated.
  Wt::WGroupBox *m_chartBox;
  DrfChart *m_chart;
  Wt::WComboBox *m_chartMode;       //Absolute | Intrinsic
  Wt::WLineEdit *m_chartDistance;   //source distance for the absolute curves

  /** Identifies the current generation; results/progress from stale runs
   (user re-started or changed things) are discarded.  Only touched from the
   session thread.
   */
  int m_generationId;

  /** Cooperative cancel flag shared with the worker thread. */
  std::shared_ptr<std::atomic<bool>> m_cancelFlag;

  std::shared_ptr<const ceelo::DetectorResponse> m_result;

  Wt::Signal<bool> m_validationChanged;
  Wt::Signal<std::shared_ptr<ceelo::DetectorResponse>> m_responseGenerated;
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> m_updatedDrf;
};//class MakeMcResponseForDrf


/** AuxWindow holding a MakeMcResponseForDrf tool; create only via
 AuxWindow::make<MakeMcResponseForDrfWindow>( seed_drf ) (see
 InterSpec::showMcResponseWindow).
 */
class MakeMcResponseForDrfWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  MakeMcResponseForDrf *tool();

protected:
  // Constructor is protected; use AuxWindow::make<MakeMcResponseForDrfWindow>().
  explicit MakeMcResponseForDrfWindow(
                std::shared_ptr<const DetectorPeakResponse> seed_drf,
                std::shared_ptr<const ceelo::GeometryDescriptor> geometry = nullptr );

  MakeMcResponseForDrf *m_tool;
};//class MakeMcResponseForDrfWindow

#endif //MakeMcResponseForDrf_h
