#ifndef FarmOptionsWindow_h
#define FarmOptionsWindow_h
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

#include <Wt/WSignal>
#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/FarmOptions.h"

#include <filesystem>
#include <regex>
#include <optional>
#include <cstdlib>
#include <string>

class InterSpec;

namespace Wt
{
  class WText;
  class WLabel;
  class WCheckBox;
  class WLineEdit;
  class WContainerWidget;
}//namespace Wt

/**
 * POD Class to store FRAM version number and path
 **/
struct FramInstall 
{
  std::filesystem::path path;
  int major;
  int minor;
};

/** Window for configuring FARM options.

 This dialog allows users to configure:
 - Master enable/disable for FARM analysis
 - GADRAS Full Spectrum Isotope ID analysis
 - RelActCalcAuto isotopics analysis
 - FRAM isotopics analysis
 - Statistical moments and a few other misc spectrum-derived quantities 
 - Output options for .farm.fertilized.n42 files (contains a remark with some JSON holding a number of values)

 Settings are stored in UserPreferences and emitted via the optionsChanged signal when saved.
 */
class FarmOptionsWindow : public AuxWindow
{
public:
  FarmOptionsWindow( InterSpec *viewer );
  virtual ~FarmOptionsWindow();

  /** Returns the current options configured in the UI. */
  Farm::FarmOptions currentOptions() const;

  /** Sets the UI to display the specified options. */
  void setOptions( const Farm::FarmOptions &opts );

  /** Signal emitted when user saves options. */
  Wt::Signal<Farm::FarmOptions> &optionsChanged() { return m_optionsChanged; }

private:
  /** Initializes the UI layout and widgets. */
  void init();

  /** Saves current options to UserPreferences and emits optionsChanged signal. */
  void saveToPreferences();

  void closeAndSave();
  
  /** Loads options from UserPreferences and updates the UI. */
  void loadFromPreferences();

  /** Updates enabled/disabled state of sub-options based on master enable. */
  void handleEnableChanged();

  /** Validates all path inputs and toggles "Wt-invalid" style accordingly.
   Empty paths are considered valid (fields are optional).
   */
  void validatePathInputs();

  /**
   * Slot to handle FRAM version toggle
   **/
  void onFramVersionToggled(Wt::WRadioButton* button);

  /**
   * Search the Program Files (x86) directory and return path to latest
   * release of requested major version.
   * e.g. User has 6.1, 7.1, and 7.2 installed
   * User selects v7 check box
   * Will return path to 7.2 executable
   **/
  std::optional<std::filesystem::path> findLatestFramExe(int requestedMajor);

  /**
   * Similar to findLatestFramExe() except it returns the directory where the
   * output is saved
   **/
  std::optional<std::filesystem::path> findLatestFramOutputDir(int requestedMajor);

  InterSpec *m_interspec;

  // Master enable
  Wt::WCheckBox *m_enableFarm;

  // GADRAS controls
  Wt::WCheckBox *m_enableGadras;
  Wt::WLineEdit *m_gadrasExePath;
  Wt::WCheckBox *m_synthesizeBackground;

  // RelActCalcAuto controls
  Wt::WCheckBox *m_enableRelAct;

  // FRAM controls
  Wt::WCheckBox *m_enableFram;
  Wt::WButtonGroup *m_framVersionGroup;
  Wt::WRadioButton *m_v6Fram;
  Wt::WRadioButton *m_v7Fram;
  Wt::WLineEdit *m_framExePath;
  Wt::WLineEdit *m_framOutputPath;

  // Output controls
  Wt::WCheckBox *m_writeFertilized;

  Wt::Signal<Farm::FarmOptions> m_optionsChanged;
};//class FarmOptionsWindow

#endif // FarmOptionsWindow_h
