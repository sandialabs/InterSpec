#ifndef PeakInfoDisplay_h
#define PeakInfoDisplay_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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
#if ( USE_SPECTRUM_CHART_D3 )
class D3SpectrumDisplayDiv;
#else
class SpectrumDisplayDiv;
#endif
namespace Wt
{
  class WImage;
  class WWidget;
  class WMenuItem;
  class WTabWidget;
  class WBorderLayout;
}//namespace Wt


class PeakInfoDisplay : public Wt::WContainerWidget
{
public:
  PeakInfoDisplay( InterSpec *viewer,
#if ( USE_SPECTRUM_CHART_D3 )
                   D3SpectrumDisplayDiv *spectrumDisplayDiv,
#else
                   SpectrumDisplayDiv *spectrumDisplayDiv,
#endif
                   PeakModel *peakModel,
                   Wt::WContainerWidget *parent = 0 );
  virtual ~PeakInfoDisplay();

  void enablePeakSearchButton( bool enable );
  
protected:
  //init() must be called exactly once
  void init();

  void createNewPeak();
  void deleteSelectedPeak();
  void enablePeakDelete( Wt::WModelIndex index );
  void disablePeakDelete();
  void confirmRemoveAllPeaks();
  void assignNuclidesFromRefLines();
  
protected:
  PeakModel *m_model;
  //m_viewer is only necessary for calling guessIsotopesForPeaks(...), could
  //  change to instead have a method which sets a member variable
  //  DetectorPeakResponse when the InterSpec::detectorChanged() signal
  //  is emitted.
  InterSpec *m_viewer;
#if ( USE_SPECTRUM_CHART_D3 )
  D3SpectrumDisplayDiv *m_spectrumDisplayDiv;
#else
  SpectrumDisplayDiv *m_spectrumDisplayDiv;
#endif
  WContainerWidget *buttonDiv;

  //variables which will be populated by createInfoTab();
  Wt::WBorderLayout *m_infoLayout;
  RowStretchTreeView *m_infoView;

  Wt::WPushButton *m_deletePeak;
  
  Wt::WPushButton *m_searchForPeaks;
};//class PeakInfoDisplay

#endif // #ifndef PeakInfoDisplay_h

