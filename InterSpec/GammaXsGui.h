#ifndef GammaXsGui_h
#define GammaXsGui_h
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

#include <array>
#include <vector>
#include <utility>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MassAttenuationTool.h" //for USE_SNL_GAMMA_ATTENUATION_VALUES
#include "InterSpec/DetectorPeakResponse.h"


namespace Wt
{
  class WText;
  class WLabel;
  class WLineEdit;
  class WGridLayout;
  class WDoubleValidator;
  class WSuggestionPopup;
}//namespace Wt

namespace SandiaDecay
{
  struct Element;
}

class MaterialDB;
class DetectorDisplay;

class GammaXsGui : public Wt::WContainerWidget
{
/*
  This class creates a Div that houses the gamma XS calculator widget
*/
public:
  //GammaXsGui(...) assumes `tool` is valid and remains accessible throughout
  //  the life of the GammaXsGui - also assumes no multithread issues with it
  GammaXsGui( MaterialDB *materialDB,
              Wt::WSuggestionPopup *materialSuggestion,
              InterSpec *specViewer,
              Wt::WContainerWidget *parent = 0 );
  virtual ~GammaXsGui();
  
  void updateDetectorCalc();
  void handleDetectorChange( std::shared_ptr<DetectorPeakResponse> det );
  void resetAnserFields();
  void handleMaterialChange();
  void calculateCrossSections();


  //parseMaterial(): throws on error
  std::vector<std::pair<const SandiaDecay::Element *, float> > parseMaterial();
  
  
  /** Handles receiving a "deep-link" url starting with "interspec://gammaxs?e=1001&m=Fe&d=1.2&t=0.1cm&r=1.2m".
   
   Example URIs:
   - "interspec://gammaxs?e=1001&m=Fe&d=1.2&t=0.1cm&r=1.2m"
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://gammaxs?e=1001&m=...", then this string would be "e=1001keV&m=...".
          This string is in standard URL format of "key1=value1&key2=value2&..." with ordering not mattering.
          Capitalization is not important.
          Assumes the string passed in has already been url-decoded.
          If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string is just the query portion of URL;
   so will look something like "e=1001&m=....", and it will not be url-encoded, and does not contain the leading '?'.
   */
  std::string encodeStateToUrl() const;
  
  
protected:
  
  std::array<Wt::WString,4> inputs() const;
  void fromInputs( const std::array<Wt::WString,4> &inputs );
  void checkAndAddUndoRedo();
  
  Wt::WLineEdit *m_energyEdit;
  Wt::WDoubleValidator *m_energyValidator;
  Wt::WLineEdit *m_materialEdit;
  Wt::WSuggestionPopup *m_materialSuggestion; //not owned by this object
  Wt::WText *m_effectiveZ;

  Wt::WText *m_totalAttenuation;
  Wt::WText *m_compton;
#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  Wt::WText *m_rayleigh;
#endif
  Wt::WText *m_photoElectric;
  Wt::WText *m_conversion;

  Wt::WLineEdit *m_density;
  Wt::WLineEdit *m_thickness;
  Wt::WText *m_transmissionFraction;  //labeled "Trans. Frac."
  float m_transmissionFractionVal;  //will be negative if not defined

  Wt::WGridLayout *m_layout;
  
  MaterialDB *m_materialDB;       //not owned by this object
  
  InterSpec* m_specViewer; //the reference holding the detector
  DetectorDisplay *m_detectorDisplay;
  Wt::WLabel *m_detectorDistanceLabel;
  Wt::WLineEdit *m_detectorDistance; //detector distance
  Wt::WLabel *m_efficiencyLabel;
  Wt::WText *m_efficiency; //detector efficiency
  Wt::WLabel *m_totalEfficiencyLabel;
  Wt::WText *m_totalEfficiency;
  Wt::WLabel *m_intrinsicEfficiencyLabel;
  Wt::WText *m_intrinsicEfficiency; //detector absolute efficiency
  Wt::WLabel *m_fractionalAngleLabel;
  Wt::WText *m_fractionalAngle; //detector fractional solid angle
  
  std::shared_ptr<const DetectorPeakResponse> m_detector;
  
  std::array<Wt::WString,4> m_prevInputs;
};//class GammaXsGui


class GammaXsWindow : public AuxWindow
{
/*
  This class creates a popup window to display the GammaXsGui.
  When the user closes the window all memmorry is cleaned up.
*/
public:
  GammaXsWindow( MaterialDB *materialDB,
                 Wt::WSuggestionPopup *materialSuggestion,
                 InterSpec* viewer );
  virtual ~GammaXsWindow();
  
  GammaXsGui *xstool();
protected:
  
  GammaXsGui *m_tool;
};//class GammaXsWindow


#endif  //ifndef GammaXsGui
