#ifndef GammaCountDialog_h
#define GammaCountDialog_h
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

#include <boost/function.hpp>

#include "InterSpec/AuxWindow.h"

class SpecMeas;
class InterSpec;
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SpectrumType : int; }


namespace Wt
{
  class WText;
  class WImage;
  class WDoubleSpinBox;
}//namespace Wt

class GammaCountDialog : public AuxWindow
{
public:
  GammaCountDialog( InterSpec *specViewer );
  virtual ~GammaCountDialog();

  //setEnergyRange(...): changes SpinBoxes, updates counts, and updates chart
  //  highlighted regions
  void setEnergyRange( double lowEnergy, double highEnergy );

protected:
  void init();
  void emitFinished();
  void handleEnergyRangeChange();
  void handleSpectrumChange( SpecUtils::SpectrumType type,
                             std::shared_ptr<SpecMeas> meas,
                             std::set<int> displaySample );

  static void setGammaCountText( Wt::WText *text, std::shared_ptr<const SpecUtils::Measurement> hist,
                                 const double scale_factor,
                                 const float minEnergy, const float maxEnergy );

  void setUnceraintyText( std::shared_ptr<const SpecUtils::Measurement> foreground,
                          std::shared_ptr<const SpecUtils::Measurement> background,
                          const double backSF,
                          const float minEnergy, const float maxEnergy );
  
protected:
  InterSpec *m_specViewer;

  size_t m_highlightRegionId;
  
  Wt::WDoubleSpinBox *m_lowerEnergy;
  Wt::WDoubleSpinBox *m_upperEnergy;

  Wt::WText *m_primaryGammaCount;
  Wt::WText *m_secondaryGammaCount;
  Wt::WText *m_backgroundGammaCount;
  Wt::WText *m_secondaryLiveTimeScale;
  Wt::WText *m_backgroundLiveTimeScale;
  Wt::WText *m_sigmaAboveBackground;
  Wt::WImage *m_nsigmaHelp;
};//class GammaCountDialog



#endif
