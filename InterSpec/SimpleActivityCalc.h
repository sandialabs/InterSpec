#ifndef SimpleActivityCalc_h
#define SimpleActivityCalc_h
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
#include <vector>
#include <memory>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/ShieldingSelect.h"

class DrfSelect;
class ShieldingSelect;
class SimpleActivityCalc;
class DetectorPeakResponse;
class MaterialDB;

namespace Wt
{
  class WText;
  class WLabel;
  class WComboBox;
  class WLineEdit;
  class WPushButton;
  class WGridLayout;
  class WSuggestionPopup;
}//namespace Wt

namespace SandiaDecay
{
  class Nuclide;
  struct Transition;
}

class PeakDef;

struct SimpleActivityCalcState
{
  double peakEnergy;
  std::string nuclideAgeStr;
  std::string distanceStr;
  int geometryType;
  ShieldingSourceFitCalc::ShieldingInfo shielding;
  
  SimpleActivityCalcState();
  
  bool operator==( const SimpleActivityCalcState &rhs ) const;
  bool operator!=( const SimpleActivityCalcState &rhs ) const;
  
  std::string encodeToUrl() const;
  void decodeFromUrl( const std::string &uri );
  
  void serialize( std::ostream &out ) const;
  void deserialize( std::istream &in );
};

class SimpleActivityCalcWindow : public AuxWindow
{
public:
  SimpleActivityCalcWindow( MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggestion,
                  InterSpec* viewer );
  
  virtual ~SimpleActivityCalcWindow();

  SimpleActivityCalc *tool();
  
protected:
  SimpleActivityCalc *m_tool;
};

class SimpleActivityCalc : public Wt::WContainerWidget
{
public:
  
  SimpleActivityCalc( MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggestion,
                  InterSpec *specViewer,
                  Wt::WContainerWidget *parent = 0 );
  
  virtual ~SimpleActivityCalc();
  
  void setPeakFromEnergy( const double energy );
  
  void handleAppUrl( std::string uri );
  std::string encodeStateToUrl() const;
  
  SimpleActivityCalcState currentState() const;
  void setState( const SimpleActivityCalcState &state );
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  void init();
  
  void handlePeakChanged();
  void handleDistanceChanged();
  void handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf );
  void handleGeometryChanged();
  void handleShieldingChanged();
  void handleSpectrumChanged();
  
  void updateResult();
  void handleOpenAdvancedTool();
  
  enum RenderActions
  {
    UpdateResult = 0x01,
    AddUndoRedoStep = 0x02,
    UpdateDisplayedSpectrum = 0x04
  };
  
  Wt::WFlags<RenderActions> m_renderFlags;
  
  InterSpec *m_viewer;
  Wt::WSuggestionPopup *m_materialSuggest;
  MaterialDB *m_materialDB;
  
  Wt::WComboBox *m_peakSelect;
  Wt::WText *m_nuclideInfo;
  Wt::WLineEdit *m_ageEdit;
  Wt::WLineEdit *m_distanceEdit;
  DrfSelect *m_detectorDisplay;
  ShieldingSelect *m_shieldingSelect;
  Wt::WComboBox *m_geometrySelect;
  Wt::WText *m_resultText;
  Wt::WText *m_errorText;
  Wt::WPushButton *m_advancedBtn;
  
  Wt::WContainerWidget *m_distanceRow;
  Wt::WContainerWidget *m_ageRow;
  Wt::WContainerWidget *m_geometryRow;
  
  std::string m_stateUri;
  
  std::shared_ptr<const PeakDef> m_currentPeak;
  const SandiaDecay::Nuclide *m_currentNuclide;
  double m_currentAge;
};

#endif //SimpleActivityCalc_h
