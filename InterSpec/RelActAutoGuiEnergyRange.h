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

#include <Wt/WContainerWidget>

#include "InterSpec/PeakDef.h" //For PeakContinuum::OffsetType
#include "InterSpec/RelActCalcAuto.h"

//Forward declerations
namespace Wt
{
  class WCheckBox;
  class WComboBox;
  class WPushButton;
}


class NativeFloatSpinBox;


class RelActAutoGuiEnergyRange : public Wt::WContainerWidget
{
public:
  RelActAutoGuiEnergyRange( Wt::WContainerWidget *parent = nullptr );
  bool isEmpty() const;
  void handleRemoveSelf();
  void handleContinuumTypeChange();
  void handleForceFullRangeChange();
  void enableSplitToIndividualRanges( const bool enable );
  void handleEnergyChange();
  void setEnergyRange( float lower, float upper );
  bool forceFullRange() const;
  void setForceFullRange( const bool force_full );
  void setContinuumType( const PeakContinuum::OffsetType type );
  void setHighlightRegionId( const size_t chart_id );
  size_t highlightRegionId() const;
  float lowerEnergy() const;
  float upperEnergy() const;
  void setFromRoiRange( const RelActCalcAuto::RoiRange &roi );
  RelActCalcAuto::RoiRange toRoiRange() const;
  Wt::Signal<> &updated();
  Wt::Signal<> &remove();
  Wt::Signal<RelActAutoGuiEnergyRange *> &splitRangesRequested();
  
protected:
  Wt::Signal<> m_updated;
  Wt::Signal<> m_remove_energy_range;
  Wt::Signal<RelActAutoGuiEnergyRange *> m_split_ranges_requested;
  
  NativeFloatSpinBox *m_lower_energy;
  NativeFloatSpinBox *m_upper_energy;
  Wt::WComboBox *m_continuum_type;
  Wt::WCheckBox *m_force_full_range;
  Wt::WPushButton *m_to_individual_rois;
  
  /// Used to track the Highlight region this energy region corresponds to in D3SpectrumDisplayDiv
  size_t m_highlight_region_id;
  
  void emitSplitRangesRequested();
};//class RelActAutoGuiEnergyRange

#endif //RelActAutoGuiEnergyRange_h
