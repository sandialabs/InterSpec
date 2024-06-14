#ifndef AddNewPeakDialog_h
#define AddNewPeakDialog_h
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

#include "InterSpec/AuxWindow.h"


//Forward declarations
class PeakModel;
class InterSpec;

namespace Wt
{
  class WText;
  class WCheckBox;
  class WComboBox;
  class WPushButton;
}//namespace Wt

namespace SpecUtils
{
  class Measurement;
}

class NativeFloatSpinBox;

/** Creates a window the user can add a peak to the spectrum, specifying the peak
 parameters, and seeing a preview of the peak as they go.
 
 A number of improvements should be made:
 - Support undo/redo
 - Switch to using d3-based spectrum plot
 - Add step-continuum types
 - Allow users to edit continuum parameters, and choose to fit each one individually or not
 */
class AddNewPeakDialog : public AuxWindow
{
public:
  static const float m_minfwhm;// = 0.05f;
  static const float m_maxfwhm;// = 450.0f;  //reasonable range of peak widths

protected:
  /** Flags for what to do in `render(flags)` funciton */
  enum RenderActions
  {
    UpdatePreview = 0x01,
  };//enum D3RenderActions
  
  Wt::WFlags<AddNewPeakDialog::RenderActions> m_renderFlags;
  
  
  InterSpec *m_viewer;
  PeakModel *m_peakModel;
  bool m_isPhone;
  std::shared_ptr<const SpecUtils::Measurement> m_meas;
  
  std::shared_ptr<PeakDef> m_candidatePeak;
  
  NativeFloatSpinBox *m_energySB;
  NativeFloatSpinBox *m_fwhmSB;
  NativeFloatSpinBox *m_areaSB;
  NativeFloatSpinBox *m_roiLowerSB;
  NativeFloatSpinBox *m_roiUpperSB;
  
  Wt::WComboBox *m_continuumType;
  
  Wt::WPushButton *m_fitBtn;
  Wt::WCheckBox *m_fitEnergy;
  Wt::WCheckBox *m_fitFWHM;
  Wt::WCheckBox *m_fitAmplitude;
  
  Wt::WText *m_chart;
  
public:
  
  AddNewPeakDialog( const float initialEnergy );
  
  /** Estimates the peak FWHM to initially use for an energy.
   
   To get the FWHM:
   1) see if there is a peak above and below the current one, if so interpolate
   2) see if there are auto-search peaks above and below the current one, if so interpolate
   3) see if DRF contains FWHM
   4) see if there is a peak above or below, and scale by sqrt
   5) Guess based on number of bins.
  */
  float estimateFWHM( const float energy );
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  void updateCandidatePeakPreview();

  void roiTypeChanged();
  
  void meanChanged();
  
  
  void fwhmChanged();
  
  
  void ampChanged();

  
  void roiRangeChanged();
  
  void doFit();
};//class AddNewPeakDialog : public AuxWindow

#endif //AddNewPeakDialog_h
