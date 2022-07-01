#ifndef RelActManualGui_h
#define RelActManualGui_h
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

#include <map>
#include <string>
#include <memory>
#include <utility>

#include <Wt/WContainerWidget>

// Forward declarations
class PeakModel;
class AuxWindow;
class InterSpec;
class RelEffChart;
class NativeFloatSpinBox;
class RowStretchTreeView;

namespace Wt
{
  class WMenu;
  class WComboBox;
  class WResource;
  class WGridLayout;
}

namespace RelActCalcManual
{
  struct RelEffSolution;
}


class RelActManualGui : public Wt::WContainerWidget
{
public:
  RelActManualGui( InterSpec *viewer, Wt::WContainerWidget *parent = nullptr );
  
  void init();

  std::shared_ptr<const RelActCalcManual::RelEffSolution> currentSolution();
  
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  void calculateSolution();
  void updateGuiWithResults();
  
  void relEffEqnFormChanged();
  void relEffEqnOrderChanged();
  void nucDataSrcChanged();
  void addUncertChanged();
  
  void handlePeaksChanged();
  
  void displayedSpectrumChanged();
  
protected:
  enum RenderActions
  {
    UpdateCalc = 0x01
  };//enum D3RenderActions
  
  Wt::WFlags<RelActManualGui::RenderActions> m_renderFlags;
  
  
  std::shared_ptr<RelActCalcManual::RelEffSolution> m_currentSolution;
  
  InterSpec *m_interspec;
  
  Wt::WGridLayout *m_layout;
  Wt::WContainerWidget *m_optionsColumn;
  
  Wt::WComboBox *m_relEffEqnForm;
  Wt::WComboBox *m_relEffEqnOrder;
  
  Wt::WContainerWidget *m_nucDataSrcHolder;
  Wt::WComboBox *m_nucDataSrc;
  
  NativeFloatSpinBox *m_matchTolerance;
  
  // The additional uncertainty is a bit confusing, so to not overwhelm users, or have them worry
  //  about it to much, we'll provide a limited number of options.
  enum class AddUncert : int
  {
    Unweighted, StatOnly, OnePercent, FivePercent, TenPercent, TwentyFivePercent,
    FiftyPercent, SeventyFivePercent, OneHundredPercent, NumAddUncert
  };//enum class AddUncert
  
  Wt::WComboBox *m_addUncertainty;
  
#if( BUILD_AS_OSX_APP )
  Wt::WAnchor *m_downloadHtmlReport;
#else
  Wt::WPushButton *m_downloadHtmlReport;
#endif
  Wt::WResource *m_htmlResource;
  
  Wt::WContainerWidget *m_peakTableColumn;
  PeakModel *m_peakModel;
  RowStretchTreeView *m_peakTable;
  
  /// All entries in this next <div> will be of class ManRelEffNucDisp
  Wt::WContainerWidget *m_nuclidesDisp;
  
  /// Keep a cache of nuclide ages around incase the user removes a nuclide, but adds it in later.
  std::map<std::string,double> m_nucAge;
  
  Wt::WMenu *m_resultMenu;
  
  RelEffChart *m_chart;
  Wt::WContainerWidget *m_results;
};//class RelActManualGui


#endif //RelActManualGui_h
