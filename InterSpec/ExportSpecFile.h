#ifndef ExportSpecFile_h
#define ExportSpecFile_h
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
#include <string>
#include <vector>
#include <memory>
#include <cstdint>

#include <Wt/WContainerWidget>

#include "InterSpec/SimpleDialog.h"

class SpecMeas;
class InterSpec;


namespace Wt
{
  class WText;
  class WMenu;
  class WCheckBox;
  class WLineEdit;
}//namespace Wt

namespace SpecUtils
{
  enum class SpectrumType : int;
  enum class SaveSpectrumAsType : int;
}

namespace ExportSpecFileTool_imp
{
  class DownloadSpectrumResource;
}

/**
 TODO:
  - Add a checkbox for excluding serial number, or other sensitive info
  - Fix up displaying for undersized displays
  - Peaks aren't re-fit if background subtracted is selected - skeleton code is in there for this.
 */

class ExportSpecFileTool : public Wt::WContainerWidget
{
public:
  /** Constructor to create tool to select from any currently loaded spectrum files. */
  ExportSpecFileTool( InterSpec *viewer, Wt::WContainerWidget *parent );
  
  /** Constructor to create tool to export only the spectrum file passed in.
   
   @param spectrum The spectrum file to allow exporting.
   @param samples The samples to use for the export; if empty, and there are multiple, then user will be given option to select samples.
   @param detectors The detectors to use for the export; if empty, and there are multiple, then user will be given option to select samples.
   */
  ExportSpecFileTool( const std::shared_ptr<const SpecMeas> &spectrum,
                      const std::set<int> &samples,
                      const std::vector<std::string> &detectors,
                      InterSpec *viewer,
                      Wt::WContainerWidget *parent );
  
  Wt::Signal<bool> &done();
  
  void emitDone( const bool exported );
  
  /** Returns the currently selected file - as held in memory; i.e., the full file, not just selected samples or whatever. */
  std::shared_ptr<const SpecMeas> currentlySelectedFile() const;
  
  /** Returns sample numbers from `currentlySelectedFile()` that will be placed in the output. */
  std::set<int> currentlySelectedSamples() const;
  
  /** Returns detector names from `currentlySelectedFile()` that will be placed in the output. */
  std::vector<std::string> currentlySelectedDetectors() const;
  
  /** Returns the currently selected save-as file type.
   
   If QR-codes are enabled, returns SpecUtils::SaveSpectrumAsType::NumTypes for it.
   */
  SpecUtils::SaveSpectrumAsType currentSaveType() const;
  
  
  /** Returns 1 (ex., SPE, CHN, etc), 2 (QR-code), or `numeric_limits<uint16_t>::max()` (N42, PCF, etc),
   depending on the currently selected save-as type.
   
   @param spec The spectrum file to be saved; optional (but saves a call to `currentlySelectedFile`),
          and only used for QR-code format to determine one or two spectra
   */
  uint16_t maxRecordsInCurrentSaveType( std::shared_ptr<const SpecMeas> spec ) const;
  
  /** Returns if this option is selected in the GUI; the `generateFileToSave()` function will have removed
   most info, but if we are saving to N42, then we'll be a bit more careful and remove more at export time.
   */
  bool removeInterSpecInfo() const;
  
  /** Generates the final spectrums file, ready to be given to the user.
   
   E.g., if user selected, all detectors will be summed, or samples numbers, etc.
   
   Throws exception if there is issue generating the file.
   */
  std::shared_ptr<const SpecMeas> generateFileToSave();
  
  
  /** Handles receiving a "deep-link" url starting with "interspec://specexport/...".
   
   Example URIs:
   - "interspec://specexport?V=1&..."
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
   "interspec://specexport?V=1&FORMAT=N42-2012&Samples=1-2&...", then this string would be "V=1&FORMAT=N42-2012&Samples=1-2&...".
   Assumes the string passed in has alaready been url-decoded.
   If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string does not include the
   "interspec://" protocol, or "specexport" path; so will look something like "V=1&FORMAT=N42-2012&Samples=1,2&..",
   and it will not be url-encoded.
   */
  std::string encodeStateToUrl() const;
  
  
  /** Converts texts ranges like "1-5,8-9,11,12" into the numerical sample numbers. */
  static std::pair<std::set<int>,std::string> sampleNumbersFromTxtRange( std::string txt,
                                                 std::shared_ptr<const SpecMeas> meas,
                                                 const bool quiet );
  
protected:
  void init();
  void updateExportEnabled();
  void handleFileSelectionChanged();
  void updateInfoAboutSelectedFile();
  void refreshSampleAndDetectorOptions();
  void selectDisplayedFile( const SpecUtils::SpectrumType type );
  void handleAllSampleChanged();
  void handleCustomSampleChanged();
  void handleCustomSampleTxtChanged();
  void handleDisplaySampleChanged( const SpecUtils::SpectrumType type );
  void handleSamplesChanged();
  void handleFormatChange();
  void handleForePlusBackChanged();
  void handleFilterDetectorCbChanged();
  
  void handleSumToSingleRecordChanged();
  void handleSumTypeToSingleRecordChanged();
  void handleBackSubForeChanged();
  void handleSumDetPerSampleChanged();
  void handleIncludeInterSpecInfoChanged();
  
#if( USE_QR_CODES )
  void handleGenerateQrCode();
#endif
  
  void updateUndoRedo();
  void scheduleAddingUndoRedo();
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;
  
  InterSpec *m_interspec;
  
  /** Wether or not the dialog is for a specific file, or any loaded file. */
  const bool m_is_specific_file;
  
  std::shared_ptr<const SpecMeas> m_specific_spectrum;
  std::set<int> m_specific_samples;
  std::vector<std::string> m_specific_detectors;
  
  std::weak_ptr<const SpecMeas> m_current_file;
  
  Wt::Signal<bool> m_done;
  
  Wt::WComboBox *m_fileSelect;
  Wt::WCheckBox *m_forePlusBack;
  
  Wt::WText *m_fileInfo;
  
  Wt::WMenu *m_formatMenu;
  
  Wt::WContainerWidget *m_samplesHolder;
  Wt::WCheckBox *m_dispForeSamples;
  Wt::WCheckBox *m_dispBackSamples;
  Wt::WCheckBox *m_dispSecondSamples;
  Wt::WCheckBox *m_allSamples;
  Wt::WCheckBox *m_customSamples;
  Wt::WLineEdit *m_customSamplesEdit;
  Wt::WCheckBox *m_filterDetector;
  /** A <div /> that will only contain WCheckBox's, where the label of each CB will be the detectors name. */
  Wt::WContainerWidget *m_detectorFilterCbs;
  
  Wt::WContainerWidget *m_optionsHolder;
  Wt::WCheckBox *m_sumAllToSingleRecord;
  Wt::WCheckBox *m_sumForeToSingleRecord;
  Wt::WCheckBox *m_sumBackToSingleRecord;
  Wt::WCheckBox *m_sumSecoToSingleRecord;
  Wt::WCheckBox *m_backSubFore;
  Wt::WCheckBox *m_sumDetsPerSample;
  Wt::WCheckBox *m_excludeInterSpecInfo;
  Wt::WCheckBox *m_excludeGpsInfo;
  
  Wt::WText *m_sampleSelectNotAppTxt;
  Wt::WText *m_optionsNotAppTxt;
  
  Wt::WText *m_msg;
  
  Wt::WPushButton *m_export_btn;
#if( USE_QR_CODES )
  Wt::WPushButton *m_show_qr_btn;
#endif
  
  ExportSpecFileTool_imp::DownloadSpectrumResource *m_resource;
  
  /** For undo/redo purposes, we will keep track of last state as a URL. */
  std::shared_ptr<const std::string> m_last_state_uri;
};//class ExportSpecFileTool


class ExportSpecFileWindow : public SimpleDialog
{
public:
  ExportSpecFileWindow( InterSpec *viewer );
  
  ExportSpecFileWindow( const std::shared_ptr<const SpecMeas> &spectrum,
                       const std::set<int> &samples,
                       const std::vector<std::string> &detectors,
                       InterSpec *viewer );
  
  void setSpecificSpectrum( const std::shared_ptr<const SpecMeas> &spectrum,
                           const std::set<int> &samples,
                           const std::vector<std::string> &detectors,
                           InterSpec *viewer );
  
  void handleAppUrl( const std::string &query_str );
  std::string encodeStateToUrl() const;
  
  void scheduleDelete();
protected:
  
  ExportSpecFileTool *m_tool;
};//class MoreNuclideInfoWindow


#endif //ExportSpecFile_h
