#ifndef RelActAutoGuiEnergyRange_h
#define RelActAutoGuiEnergyRange_h
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

#include <string>
#include <utility>
#include <optional>

#include <Wt/WContainerWidget.h>

#include "InterSpec/PeakDef.h" //For PeakContinuum::OffsetType
#include "InterSpec/RelActCalcAuto.h"

//Forward declerations
namespace Wt
{
  class WText;
  class WLabel;
  class WComboBox;
  class WPushButton;
}


class NativeFloatSpinBox;


/** The GUI for one ROI (`RelActCalcAuto::RoiRange`) of the "Isotopics by nuclides" tool.

 Besides the energies and continuum type, the user picks how the energies are interpreted (a fixed
 range, the lines to include with the extent from the peak shape, or a range to split up by lines),
 and for the non-fixed types, may override how far the ROI extends past its lines.
 */
class RelActAutoGuiEnergyRange : public Wt::WContainerWidget
{
public:
  RelActAutoGuiEnergyRange();
  bool isEmpty() const;
  void handleRemoveSelf();
  void handleContinuumTypeChange();
  void handleRangeTypeChange();
  void handleEdgeChange();
  void enableSplitToIndividualRanges( const bool enable );
  void handleEnergyChange();
  void setEnergyRange( float lower, float upper );

  /** Returns true if the ROI's range type is `Fixed`. */
  bool forceFullRange() const;

  /** Sets the range type to `Fixed`, or if `force_full` is false, back to the (non-fixed) type the ROI
   was loaded with (`LineAnchored` if it has not been anything else).  A line-anchored ROI made fixed
   takes the range it was last fit over. */
  void setForceFullRange( const bool force_full );

  RelActCalcAuto::RoiRange::RangeLimitsType rangeLimitsType() const;
  void setRangeLimitsType( const RelActCalcAuto::RoiRange::RangeLimitsType type );

  void setContinuumType( const PeakContinuum::OffsetType type );
  void setHighlightRegionId( const size_t chart_id );
  size_t highlightRegionId() const;
  float lowerEnergy() const;
  float upperEnergy() const;
  void setFromRoiRange( const RelActCalcAuto::RoiRange &roi );
  RelActCalcAuto::RoiRange toRoiRange() const;

  /** Shows the energy range this ROI was actually fit over (which for the non-fixed range types
   comes from the peak shape).  Pass a `lower >= upper` (e.g., NaNs) to clear it. */
  void setFitRange( const double lower, const double upper );

  /** The range the ROI was fit over, if known, otherwise its lower and upper energy. */
  std::pair<double,double> displayRange() const;

  /** Sets the default edge values, shown as placeholder text for edges without an override. */
  void setEdgeDefaults( const RelActCalcAuto::RoiEdge &lower, const RelActCalcAuto::RoiEdge &upper );

  /** The value an edge field shows: for a coverage field (percent of the peak area covered, in
   (50,100)) the tail fraction, otherwise the sideband in FWHM (in [0,100)).  Empty or out-of-range
   fields give nullopt. */
  static std::optional<double> edgeFieldValue( NativeFloatSpinBox *field, const bool is_coverage );

  Wt::Signal<> &updated();
  Wt::Signal<> &remove();
  Wt::Signal<RelActAutoGuiEnergyRange *> &splitRangesRequested();

protected:
  /** Updates things after the range type changed; a line-anchored ROI changed to a range type
   takes the range it was last fit over as its energy range. */
  void rangeTypeChanged();

  /** Shows/hides/relabels things for the current range type. */
  void updateForRangeType();

  Wt::Signal<> m_updated;
  Wt::Signal<> m_remove_energy_range;
  Wt::Signal<RelActAutoGuiEnergyRange *> m_split_ranges_requested;

  Wt::WLabel *m_lower_label;
  Wt::WLabel *m_upper_label;
  NativeFloatSpinBox *m_lower_energy;
  NativeFloatSpinBox *m_upper_energy;

  /** Item 0 is "Auto"; item i+1 is `PeakContinuum::OffsetType(i)`. */
  Wt::WComboBox *m_continuum_type;

  /** The range limits type; see `range_type_index(...)` in the .cpp for the item order. */
  Wt::WComboBox *m_range_type;

  Wt::WPushButton *m_to_individual_rois;

  /** How far a non-fixed ROI extends past its lines: peak coverage (%) and continuum sideband (FWHM)
   below and above; empty means use the default. */
  Wt::WContainerWidget *m_edges;
  NativeFloatSpinBox *m_lower_coverage;
  NativeFloatSpinBox *m_lower_sideband;
  NativeFloatSpinBox *m_upper_coverage;
  NativeFloatSpinBox *m_upper_sideband;

  /** Shows the range the ROI was actually fit over. */
  Wt::WText *m_fit_range_txt;
  std::pair<double,double> m_fit_range;

  /** The continuum type an automatic continuum starts from. */
  PeakContinuum::OffsetType m_auto_continuum_seed;

  /** The non-fixed range type to go back to, if un-forcing the full range. */
  RelActCalcAuto::RoiRange::RangeLimitsType m_non_fixed_type;

  /** The range type the row is currently set up for, so a change of type knows what it was. */
  RelActCalcAuto::RoiRange::RangeLimitsType m_shown_type;

  /// Used to track the Highlight region this energy region corresponds to in D3SpectrumDisplayDiv
  size_t m_highlight_region_id;

  void emitSplitRangesRequested();
};//class RelActAutoGuiEnergyRange

#endif //RelActAutoGuiEnergyRange_h
