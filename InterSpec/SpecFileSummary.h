#ifndef SpecFileSummary_h
#define SpecFileSummary_h
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
#include <memory>

#include "InterSpec/AuxWindow.h"

/*
 * This class gives a summary of file parameters.  Currently SpecUtils::SpecFile
 * isnt really geared towards keeping this information, so this class is a
 * little lacking...
 *XXX - this class should allow modifying of the spectrum informtion as well
*/

class SpecMeas;
class InterSpec;
class AnaResultDisplay;
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SpectrumType : int; }

namespace Wt
{
  class WText;
  class WComboBox;
  class WTextArea;
  class WGroupBox;
  class WLineEdit;
  class WPushButton;
  class WButtonGroup;
  class WIntValidator;
  class WInteractWidget;
  class WSelectionBox;
}//namespace Wt


class SpecFileSummary : public AuxWindow
{
public:
  enum AllowModifyStatus { kAllowModify, kDontAllowModify };

  enum EditableFields
  {
    kDisplayedLiveTime,
    kDisplayedRealTime,
    kTimeStamp,
    kPosition,
    kDescription,
    kSourceType,
    kMeasurmentRemarks,
    kFileRemarks,
    kFilename,
    kUuid,
    kLaneNumber,
    kMeasurement_location_name,
    kInspection,
    kInstrument_type,
    kManufacturer,
    kInstrument_model,
    kInstrument_id
  };//enum EditableFields

public:
  SpecFileSummary( const SpecUtils::SpectrumType spec_type, InterSpec *specViewer );
  virtual ~SpecFileSummary();

protected:
  void init();

  void showRiidAnalysis();
  void showMultimedia();
  
  void handleUserChangeSampleNum();
  void handleUserIncrementSampleNum( bool increment );

  void handleSpectrumTypeChanged();
  void handleAllowModifyStatusChange();

  void updateDisplayFromMemory();
  void updateMeasurmentFieldsFromMemory();

  void handleFieldUpdate( EditableFields field );

  void handleSpectrumChange( const SpecUtils::SpectrumType type,
                            const std::shared_ptr<SpecMeas> &meas,
                            const std::set<int> &displaySample,
                            const std::vector<std::string> &detectors );

  void reloadCurrentSpectrum();

  //currentMeasurment(): returns empty ptr on error
  std::shared_ptr<const SpecUtils::Measurement> currentMeasurment() const;

#if( USE_GOOGLE_MAP )
  void showGoogleMap();  
#elif( USE_LEAFLET_MAP )
  void showLeafletMap();
#endif
  
protected:
  InterSpec *m_specViewer;

  Wt::WButtonGroup *m_allowEditGroup;

  Wt::WContainerWidget *m_displaySampleDiv;
  Wt::WLineEdit        *m_displaySampleNumEdit;
  Wt::WText            *m_displayedPreText;
  Wt::WText            *m_displayedPostText;
  Wt::WInteractWidget  *m_nextSampleNumButton;
  Wt::WInteractWidget  *m_prevSampleNumButton;
  Wt::WIntValidator    *m_displaySampleNumValidator;

  //These following fields are specific to a single Measurment object
  Wt::WText *m_gammaCPS;
  Wt::WText *m_gammaSum;
  Wt::WText *m_neutronCPS;
  Wt::WText *m_neutronSum;
  Wt::WLineEdit *m_displayedLiveTime;
  Wt::WLineEdit *m_displayedRealTime;
  Wt::WLineEdit *m_timeStamp;

  Wt::WText *m_energyRange;
  Wt::WText *m_numberOfBins;

  Wt::WText *m_detector;
  Wt::WText *m_sampleNumber;
  Wt::WTextArea *m_measurmentRemarks;

  Wt::WLineEdit *m_longitude;
  Wt::WLineEdit *m_latitude;
  Wt::WLineEdit *m_gpsTimeStamp;
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  Wt::WPushButton *m_showMapButton;
#endif
  Wt::WLineEdit *m_title;
  Wt::WComboBox *m_source;

  Wt::WButtonGroup *m_spectraGroup;

  //Fields below here bellong to the whole file
  Wt::WTextArea *m_fileRemarks;
  Wt::WText     *m_sizeInMemmory;

  Wt::WLineEdit   *m_filename;
  Wt::WLineEdit   *m_uuid;
  Wt::WLineEdit   *m_laneNumber;
  Wt::WLineEdit   *m_measurement_location_name;
  Wt::WPushButton *m_ana_button;
  Wt::WPushButton *m_multimedia_button;
  Wt::WLineEdit   *m_inspection;
  Wt::WLineEdit   *m_instrument_type;
  Wt::WLineEdit   *m_manufacturer;
  Wt::WLineEdit   *m_instrument_model;
  Wt::WLineEdit   *m_instrument_id;

  Wt::WPushButton *m_reloadSpectrum;
};//class SpecFileSummary


#endif  //SpecFileSummary_h
