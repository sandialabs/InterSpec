#ifndef PeakInfoDisplay_h
#define PeakInfoDisplay_h
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

#include <Wt/WModelIndex>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/RowStretchTreeView.h"

class PeakModel;
class InterSpec;
class D3SpectrumDisplayDiv;

namespace Wt
{
  class WImage;
  class WLabel;
  class WWidget;
  class WMenuItem;
  class WTabWidget;
  class WGridLayout;
}//namespace Wt


//There is a `ColorDelegate` coded that starts to let the user select color in
//  the "Peak Manager", but it isnt fully working yet
#define ALLOW_PEAK_COLOR_DELEGATE 0

class PeakInfoDisplay : public Wt::WContainerWidget
{
public:
  PeakInfoDisplay( InterSpec *viewer,
                   D3SpectrumDisplayDiv *spectrumDisplayDiv,
                   PeakModel *peakModel,
                   Wt::WContainerWidget *parent = 0 );
  virtual ~PeakInfoDisplay();

  void enablePeakSearchButton( bool enable );
  
  void handleChartLeftClick( const double energy );
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  /** Currently just makes the buttons at the bottom be minimal. */
  void setNarrowPhoneLayout( const bool narrow );
#endif
  
protected:
  //init() must be called exactly once
  void init();

  void createNewPeak();
  void deleteSelectedPeak();
  void enablePeakDelete( Wt::WModelIndex index );
  void handleSelectionChanged();
  void disablePeakDelete();
  void confirmRemoveAllPeaks();
  void removeAllPeaks();
  void assignNuclidesFromRefLines();
  
  /** Copies the current peak CSV data to `$(this).data({'CsvFullData', 'CsvNoHeaderData', 'CsvCompactData'})`
   Called menu copy button is clicked, and copy options menu shown.
   */
  void copyCsvPeakDataToClient();
  /** Removes the data from the JS; called when copy menu is hidden. */
  void removeCsvPeakDatafromClient();
  
protected:
  PeakModel *m_model;
  //m_viewer is only necessary for calling guessIsotopesForPeaks(...), could
  //  change to instead have a method which sets a member variable
  //  DetectorPeakResponse when the InterSpec::detectorChanged() signal
  //  is emitted.
  InterSpec *m_viewer;
  D3SpectrumDisplayDiv *m_spectrumDisplayDiv;

  //variables which will be populated by createInfoTab();
  Wt::WGridLayout *m_infoLayout;
  RowStretchTreeView *m_infoView;

  Wt::WInteractWidget *m_deletePeak;
  
  Wt::WPushButton *m_searchForPeaks;
  
  Wt::WPushButton *m_clearPeaksButton;
  Wt::WPushButton *m_nucFromRefButton;
  Wt::WLabel *m_peakAddRemoveLabel;
};//class PeakInfoDisplay

#endif // #ifndef PeakInfoDisplay_h

