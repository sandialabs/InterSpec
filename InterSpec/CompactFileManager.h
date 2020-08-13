#ifndef CompactFileManager_h
#define CompactFileManager_h
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

#include <Wt/WContainerWidget>

namespace Wt
{
  class WText;
  class WSlider;
  class WComboBox;
  class WLineEdit;
  class WDoubleSpinBox;
  class WPushButton;
  class WEvent;
} // namespace Wt


namespace SpecUtils{ enum class SpectrumType : int; }

class SpecMeas;
class InterSpec;
class SpectraFileModel;
class SpecMeasManager;


class CompactFileManager : public Wt::WContainerWidget
{
  /****************************************************************************\
  | CompactFileManager is a highly stripped down file manager for use in the
  | tool tab.
  | The main file manager is to big for tool tab space. It is intended to do only
  | simple functions: set fore/back ground and secondary spectrums, do uploads
  \****************************************************************************/

public:
  enum DisplayMode{ TopToBottom, LeftToRight, Tabbed };
  
  CompactFileManager( SpecMeasManager *fullFileManager,
                      InterSpec *hostViewer,
                      DisplayMode mode,
                      Wt::WContainerWidget *parent = 0 );
  ~CompactFileManager();

  void refreshContents();
  void handleDisplayChange( SpecUtils::SpectrumType spectrum_type,
                            const std::shared_ptr<SpecMeas> meas,
                            const std::set<int> &sample_numbers,
                            const std::vector<std::string> &detectors );

  void handleUserChangeSampleNum( SpecUtils::SpectrumType spectrum_type );
  void handleUserIncrementSampleNum( SpecUtils::SpectrumType spectrum_type, bool increment );
  static void handleUserIncrementSampleNum( SpecUtils::SpectrumType spectrum_type, bool increment, InterSpec *hostViewer, SpectraFileModel *files, CompactFileManager* cfm);
    
  static void changeToSampleNum( int sampleNum, SpecUtils::SpectrumType spectrum_type, InterSpec *hostViewer, CompactFileManager* cfm );

  void changeToSampleRange( int firstSample, int lastSample,
                            SpecUtils::SpectrumType spectrum_type );

  void handleFileChangeRequest( int row, SpecUtils::SpectrumType spectrum_type );

  //handleSampleNumEditBlur(): load-spuctrum if the desired sample nums changed
  void handleSampleNumEditBlur( SpecUtils::SpectrumType spectrum_type );

  /** Function called when the user scales a chart.  Updates displayed numbers,
   and shows "Normalize" button.
   */
  void handleSpectrumScale( const double scale, SpecUtils::SpectrumType spectrum_type );
  
protected:
  
  /** Updates the displayed scale factor numbers to user.
   Does not show the "Normalize" button, or notify the parent InterSpec instance.
   */
  void updateDisplayedScaleFactorNumbers( const double sf, const SpecUtils::SpectrumType type );
  
//  void handleSliderChanged( const int slidervalue, const SpecUtils::SpectrumType type );
  void handleUserEnterdScaleFactor( const SpecUtils::SpectrumType type );
  
  void handleUserEnterdScaleFactorWheel( const SpecUtils::SpectrumType type, Wt::WMouseEvent e);
  
  void handleRenormalizeByLIveTime( const SpecUtils::SpectrumType type );
  
private:
  Wt::WComboBox *m_selects[3];
  Wt::WText *m_displayedPreTexts[3];
  Wt::WText *m_displayedPostTexts[3];
  Wt::WContainerWidget *m_sampleDivs[3];
  Wt::WLineEdit *m_displaySampleNumEdits[3];
  Wt::WInteractWidget *m_nextSampleNumButtons[3];
  Wt::WInteractWidget *m_prevSampleNumButtons[3];
  Wt::WDoubleSpinBox *m_scaleValueTxt[3];  //could use a WInPlaceEdit
  Wt::WPushButton *m_rescaleByLiveTime[3];
  
  //We want to avoid un-necassarily re-loading the data on blur of the sample
  //  number edit, but we do wnat to re-load the data if the text shows
  //  something different than the spectrum, therefore we'll track the sample
  //  number text (see also handleSampleNumEditBlur(...) and
  //  handleDisplayChange(...)) and only reload the spectrum data if it has
  //  changed between the last spectrum data load and the blurr.
  Wt::WString m_previousSampleTxt[3];

  //When exactly one sample number is dispayed, when there is exactly one
  //  detector, if the measurment has a title, it will be displayed using
  //  m_foregroundTitle.  In the future it may be worth while to make this a
  //  WInPlaceEdit to allow the user an easy way to change it.
  Wt::WText *m_foregroundTitle;
  
  // A link to the file manager, as that's where everything's stored
  SpectraFileModel *m_files;
  // as well as the host viewer
  InterSpec *m_interspec;
  // and the file manager, which is used for rapid upload/application
  SpecMeasManager *m_fullFileManager;
  
  const DisplayMode m_displayMode;
}; //class CompactFileManager

#endif


