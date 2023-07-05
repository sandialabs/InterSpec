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
class NativeFloatSpinBox;
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SpectrumType : int; }


namespace Wt
{
  class WText;
  class WImage;
}//namespace Wt

class GammaCountDialog : public AuxWindow
{
public:
  GammaCountDialog( InterSpec *specViewer );
  virtual ~GammaCountDialog();

  //setEnergyRange(...): changes SpinBoxes, updates counts, and updates chart
  //  highlighted regions
  void setEnergyRange( double lowEnergy, double highEnergy );

  /** Handles receiving a "deep-link" url starting with "interspec://specsum/...".
   
   Example URIs:
   - "interspec://specsum?V=1&LOW=180&HIGH=190"
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://specsum?V=1&LOW=180&HIGH=190", then this string would be "V=1&LOW=180&HIGH=190".
          Assumes the string passed in has alaready been url-decoded.
          If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string does not include the
   "interspec://" protocol, or "specsum" path; so will look something like "V=1&LOW=180&HIGH=190",
   and it will not be url-encoded (although this shouldnt be needed).
   */
  std::string encodeStateToUrl() const;
  
protected:
  void init();
  void emitFinished();
  void handleEnergyRangeChange();
  void handleSpectrumChange( const SpecUtils::SpectrumType type,
                             const std::shared_ptr<SpecMeas> &meas,
                             const std::set<int> &displaySample,
                             const std::vector<std::string> &displayedDetectors );

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
  
  NativeFloatSpinBox *m_lowerEnergy;
  NativeFloatSpinBox *m_upperEnergy;

  Wt::WString m_prevLowerEnergy;
  Wt::WString m_prevUpperEnergy;
  
  Wt::WText *m_primaryGammaCount;
  Wt::WText *m_secondaryGammaCount;
  Wt::WText *m_backgroundGammaCount;
  Wt::WText *m_secondaryLiveTimeScale;
  Wt::WText *m_liveTimeScaleNote;
  Wt::WText *m_sigmaAboveBackground;
  Wt::WImage *m_nsigmaHelp;
};//class GammaCountDialog



#endif
