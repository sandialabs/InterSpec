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
#include <optional>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/GammaInteractionCalc.h"

class DrfSelect;
class ShieldingSelect;
class SimpleActivityCalc;
class DetectorPeakResponse;
class DetectorDisplay;
class MaterialDB;
class GadrasScatterTable;

namespace Wt
{
  class WText;
  class WLabel;
  class WComboBox;
  class WLineEdit;
  class WPushButton;
  class WGridLayout;
  class WSuggestionPopup;
  class WCheckBox;
}//namespace Wt

namespace rapidxml{ template<class Ch> class xml_node; }

namespace SandiaDecay
{
  class Nuclide;
  struct Transition;
}

class PeakDef;

namespace SpecUtils { class Measurement; }
namespace ShieldingSourceFitCalc { struct ShieldingInfo; }


enum class SimpleActivityGeometryType : int
{
  Point = 0,
  Plane = 1,
  SelfAttenuating = 2,
  TraceSrc = 3
};

struct SimpleActivityCalcState
{
  double peakEnergy = -1.0;
  std::string nuclideName;
  std::string nuclideAgeStr;
  std::string distanceStr;
  SimpleActivityGeometryType geometryType = SimpleActivityGeometryType::Point;
  std::optional<ShieldingSourceFitCalc::ShieldingInfo> shielding;
  bool backgroundSubtract = false;
  std::shared_ptr<const PeakDef> fit_background_peak;
  
  SimpleActivityCalcState();
  
  bool operator==( const SimpleActivityCalcState &rhs ) const;
  bool operator!=( const SimpleActivityCalcState &rhs ) const;
  
  std::string encodeToUrl() const;
  void decodeFromUrl( const std::string &uri, MaterialDB *materialDB = nullptr );
  
  void serialize( rapidxml::xml_node<char> * const parent_node ) const;
  void deSerialize( const ::rapidxml::xml_node<char> *src_node, MaterialDB *materialDB = nullptr );
  
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  static void equalEnough( const SimpleActivityCalcState &lhs, const SimpleActivityCalcState &rhs );
#endif
  
  static const int sm_xmlSerializationVersion = 0;
};

struct SimpleActivityCalcInput
{
  std::shared_ptr<const PeakDef> peak;
  std::shared_ptr<const PeakDef> background_peak;
  std::shared_ptr<const DetectorPeakResponse> detector;
  std::shared_ptr<const SpecUtils::Measurement> foreground;
  std::shared_ptr<const SpecUtils::Measurement> background;
  double distance = 0.0;
  SimpleActivityGeometryType geometryType = SimpleActivityGeometryType::Point;
  std::vector<ShieldingSourceFitCalc::ShieldingInfo> shielding;
  double age = 0.0;
};//struct SimpleActivityCalcInput


struct SimpleActivityCalcResult
{
  double activity;
  double activityUncertainty;
  double nuclideMass;
  bool isSelfAttenuating;
  double sourceDimensions;
  bool successful;
  std::string errorMessage;
  /** Dose of the source at the detector; only computed for simple point sources, for non-fixed-geometry detectors.
   Accounts for all gammas in source, as well as scatter continuum.
   */
  std::optional<double> source_dose;
  
  SimpleActivityCalcResult() 
    : activity(0.0), activityUncertainty(0.0), nuclideMass(0.0),
      isSelfAttenuating(false), sourceDimensions(0.0),
      successful(false), errorMessage("")
  {}
};//struct SimpleActivityCalcResult


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


/** The actual widget to display and calculate the activity from a single peak.
 
 TODO: implement creating the "advanced" ShieldingSourceDisplay display
 */
class SimpleActivityCalc : public Wt::WContainerWidget
{
public:
  
  SimpleActivityCalc( MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggestion,
                  InterSpec *specViewer,
                  Wt::WContainerWidget *parent = 0 );
  
  virtual ~SimpleActivityCalc();
  
  void setPeakFromEnergy( const double energy );
  
  /** Sets the state using the URI query (i.e., everything after the '?' character) description.
   Example input: "V=1&E=185.678635&N=U235&AGE=20.00 y&DIST=1 m&GEOM=Point"
   */
  void handleAppUrl( const std::string &query_str );
  std::string encodeStateToUrl() const;
  
  std::shared_ptr<SimpleActivityCalcState> currentState() const;
  void setState( const SimpleActivityCalcState &state );
  
  void addUndoRedoPoint();
  
  static SimpleActivityGeometryType geometryTypeFromString( const std::string& key );
  static std::string geometryTypeToStringKey( SimpleActivityGeometryType type );
  static std::string to_str( SimpleActivityGeometryType type );
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  void init();
  
  void handlePeakChanged();
  void handleDistanceChanged();
  void handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf );
  void handleGeometryChanged();
  void handleShieldingChanged();
  void handleSpectrumChanged();
  void handlePeaksChanged();
  void handleBackgroundSubtractChanged();
  void setFitBackgroundPeak( const std::shared_ptr<const PeakDef> peak );
  void updateResult();
  void handleAgeChanged();
  //void handleOpenAdvancedTool();
  
  void updatePeakList();
  void updateNuclideInfo();
  std::shared_ptr<const PeakDef> getCurrentPeak() const;
  
  int findBestReplacementPeak( std::shared_ptr<const PeakDef> targetPeak ) const;
  
  void updateBackgroundSubtractVisibility();
  std::shared_ptr<const PeakDef> findOverlappingBackgroundPeak() const;
  
  void handleAddUndoPoint();
  
  SimpleActivityCalcInput createCalcInput() const;
  static SimpleActivityCalcResult performCalculation( const SimpleActivityCalcInput& input,
                                                     const GadrasScatterTable * const scatter = nullptr );
  
  void updateGeometryOptions();
  
  enum RenderActions
  {
    UpdateResult = 0x01,
    AddUndoRedoStep = 0x02
  };
  
  void scheduleRender( RenderActions action );
  
  Wt::WFlags<RenderActions> m_renderFlags;
  bool m_haveRendered;
  
  InterSpec *m_viewer;
  Wt::WSuggestionPopup *m_materialSuggest;
  MaterialDB *m_materialDB;
  
  Wt::WComboBox *m_peakSelect;
  Wt::WText *m_nuclideInfo;
  Wt::WLineEdit *m_ageEdit;
  Wt::WLineEdit *m_distanceEdit;
  DetectorDisplay *m_detectorDisplay;
  ShieldingSelect *m_shieldingSelect;
  Wt::WComboBox *m_geometrySelect;
  Wt::WText *m_resultText;
  Wt::WText *m_errorText;
  Wt::WCheckBox *m_backgroundSubtractCheck;
  //Wt::WPushButton *m_advancedBtn;
  

  Wt::WContainerWidget *m_ageRow;
  Wt::WContainerWidget *m_distanceRow;
  Wt::WContainerWidget *m_geometryRow;
  Wt::WContainerWidget *m_backgroundSubtractRow;
  
  std::shared_ptr<const PeakDef> m_currentPeak;
  std::shared_ptr<const PeakDef> m_fit_background_peak;
  
  std::vector<double> m_peakEnergies;  // Stores energy values corresponding to m_peakSelect items
  
  std::shared_ptr<const SimpleActivityCalcState> m_previous_state;
};

#endif //SimpleActivityCalc_h
