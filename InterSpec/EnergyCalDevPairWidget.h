#ifndef EnergyCalDevPairWidget_h
#define EnergyCalDevPairWidget_h
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

#include <tuple>
#include <vector>
#include <utility>

#include <Wt/WSignal>
#include <Wt/WContainerWidget>

class NativeFloatSpinBox;

namespace Wt
{
  class WCheckBox;
}//namespace Wt


namespace EnergyCalImp
{

class DeviationPairDisplay;

/** A single deviation-pair row: an editable energy and offset (NativeFloatSpinBox), a delete
 button, and - optionally - a "Fit offset" checkbox (used by the multi-file calibration dialog).
 */
class DevPair : public Wt::WContainerWidget
{
public:
  DevPair( const bool show_fit_offset, Wt::WContainerWidget *parent = 0 );
  void setDevPair( const std::pair<float,float> &d );
  std::pair<float,float> devPair() const;
  bool fitOffset() const;         //false if the "Fit offset" checkbox isnt shown
  void setFitOffset( const bool fit );
  void visuallyIndicateChanged();

  NativeFloatSpinBox *m_energy, *m_offset;
  Wt::WCheckBox *m_fitOffset;     //nullptr when not shown
  Wt::WContainerWidget *m_delete;
  friend class DeviationPairDisplay;
};//class DevPair


/** Editable list of deviation pairs (add / edit / delete), as used in the Energy Calibration tab
 and the multi-file calibration dialog.
 */
class DeviationPairDisplay : public Wt::WContainerWidget
{
public:
  /** @param show_fit_offsets if true, each row gets a "Fit offset" checkbox. */
  DeviationPairDisplay( const bool show_fit_offsets = false, Wt::WContainerWidget *parent = 0 );

  void setDeviationPairs( std::vector< std::pair<float,float> > d );
  std::vector< std::pair<float,float> > deviationPairs() const;

  /** Like #deviationPairs, but each entry also carries the state of that rows "Fit offset"
   checkbox; sorted by energy.  (false where the checkbox isnt shown.) */
  std::vector< std::tuple<float,float,bool> > deviationPairsAndFit() const;

  void removeDevPair( DevPair *devpair );
  DevPair *newDevPair( const bool emitChangedNow );
  void setInvalidValues();
  void setValidValues();

  /** Makes the (inner) list of deviation-pair rows scroll vertically once it exceeds the given
   height (used by the multi-file dialog to keep the section compact). */
  void setPairsAreaMaxHeight( const int max_height_px );

  enum class UserFieldChanged : int
  {
    AddedDeviationPair,
    RemovedDeviationPair,
    EnergyChanged,
    OffsetChanged
  };//enum UserFieldChanged

  /** Signal that gets emited when dev pair is added, deleted, or changed.
   Int argument is of type UserFieldChanged enum.  (Toggling a "Fit offset" checkbox does NOT
   emit this - it is only read when a fit is performed.)
   */
  Wt::Signal<int> &changed();

protected:
  void sortDisplayOrder( const bool indicateVisually );
  void emitChanged( const UserFieldChanged whatChanged );

  Wt::Signal<int> m_changed;
  Wt::WContainerWidget *m_pairs;
  bool m_show_fit_offsets;
};//class DeviationPairDisplay

}//namespace EnergyCalImp

#endif //EnergyCalDevPairWidget_h
