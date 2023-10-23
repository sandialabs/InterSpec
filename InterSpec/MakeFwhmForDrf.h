#ifndef MakeFwhmForDrf_h
#define MakeFwhmForDrf_h
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

#include <memory>
#include <vector>

#include <Wt/WContainerWidget>

class PeakDef;
class AuxWindow;
class InterSpec;
class NativeFloatSpinBox;
class DetectorPeakResponse;

namespace Wt
{
  class WComboBox;
}

/* Adds FWHM functional information to the DRF passed in; resulting
 DRF can be set to the current one for the application
 */
class MakeFwhmForDrf : public Wt::WContainerWidget
{
public:
  MakeFwhmForDrf( InterSpec *viewer,
                 std::shared_ptr<const DetectorPeakResponse> drf,
                 Wt::WContainerWidget *parent = nullptr );
  
  virtual ~MakeFwhmForDrf();
  
  /** Create a AuxWindow with the MakeFwhmForDrf widget as the primary content.
      Returns pointer to the created AuxWindow, but is will already be shown,
      and have the signals to delete it when closed hooked up, so you probably
      wont need the window.
   */
  static AuxWindow *makeAddFwhmToDrfWindow();

  Wt::Signal<bool> &validationChanged();
  bool isValidFwhm() const;
  
  void setToDrf();
  
protected:
  InterSpec *m_interspec;
  std::shared_ptr<const DetectorPeakResponse> m_orig_drf;
  std::vector<std::shared_ptr<const PeakDef>> m_user_peaks;
  std::vector<std::shared_ptr<const PeakDef>> m_auto_fit_peaks;
  
  Wt::WComboBox *m_fwhmEqnType;
  Wt::WComboBox *m_sqrtEqnOrder;
  std::vector<NativeFloatSpinBox *> m_parEdits;
  
  Wt::Signal<bool> m_validationChanged;
};//class MakeFwhmForDrf


#endif //MakeFwhmForDrf_h

