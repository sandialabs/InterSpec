#ifndef FluxTool_h
#define FluxTool_h
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
#include <string>
#include <vector>
#include <memory>

#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

//  ToDo:
//    - Better and more-consistent printing to appropriate number of significant figures.
//    - Havent fully tested that copying to clipboard will work everywhere.
//    - Can maybe improve copying to clipboard using the clipboard API.
//    - Have some capability to automatically fit for a number of pre-defined
//      peaks for easier batch processing of like Co60 density measurements.

//I think copying to the clipboard is working well, but leaving this optional
//  for the moment until I do some more testing.
#define FLUX_USE_COPY_TO_CLIPBOARD 1

class InterSpec;
class FluxToolWidget;
class DetectorDisplay;
class RowStretchTreeView;

namespace Wt
{
  class WText;
  class WLineEdit;
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  class WPushButton;
#endif
  class WButtonGroup;
}//namespace Wt

namespace FluxToolImp
{
  class FluxModel;
  class FluxCsvResource;
}//namespace FluxToolImp


class FluxToolWindow : public AuxWindow
{
public:
  FluxToolWindow( InterSpec *viewer );
  
  virtual ~FluxToolWindow();
  
  /** See #FluxToolWidget::handleAppUrl */
  void handleAppUrl( const std::string &query_str );
  
  /** See #FluxToolWidget::encodeStateToUrl */
  std::string encodeStateToUrl() const;
  
protected:
  FluxToolWidget *m_fluxTool;
  
  friend class FluxToolWidget;
};//class FluxToolWindow


class FluxToolWidget : public Wt::WContainerWidget
{
public:
  enum class DisplayInfoLevel
  {
    /** Only energy, nuclide, and gammas into 4pi are shown.
     Gammas into 4pi uncertainty gets its own column in CSV, and is given as a percent uncertainty.
     */
    Simple,
    
    /** Nuclide, IntrinsicEff, GeometricEff, FluxOnDet columns are NOT shown.
     Uncertainties get own column in CSV, as actual value (e.g., not percent).
     */
    Normal,
    
    /** All columns are shown.
     Uncertainties are placed in their own column in CSV, as the actual value (e.g., not percent).
     */
    Extended
  };//enum class DisplayInfoLevel
  
  
public:
  FluxToolWidget( InterSpec *viewer,
                  Wt::WContainerWidget *parent = 0 );
  
  
  virtual ~FluxToolWidget();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  enum FluxColumns
  {
    FluxEnergyCol,
    FluxNuclideCol,
    FluxPeakCpsCol,
    FluxIntrinsicEffCol,
    FluxGeometricEffCol,
    FluxFluxOnDetCol,
    FluxFluxPerCm2PerSCol,
    FluxGammasInto4PiCol,
    FluxNumColumns
  };//enum FluxColumns
  
  Wt::Signal<> &tableUpdated();
  
  DisplayInfoLevel displayInfoLevel() const;
  
  
  /** Handles receiving a "deep-link" url starting with "interspec://flux?dist=1.2m&display=low".
   
   Example URIs:
   - "interspec://flux?VER=1&dist=1.2m&display=low"
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://flux?dist=1.2m&display=low", then this string would be "dist=1.2m&display=low".
          This string is in standard URL format of "key1=value1&key2=value2&..." with ordering not mattering.
          Capitalization is not important.
          Assumes the string passed in has already been url-decoded.
          If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string is just query portion of of URL,
   so will look something like "dist=1.2m&display=low", and it will not be url-encoded.
   */
  std::string encodeStateToUrl() const;
  
protected:
  void init();
  void distanceUpdated();
  void setTableNeedsUpdating();
  void refreshPeakTable();
  
  
  void setDisplayInfoLevel( const DisplayInfoLevel disptype );
  
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  void tableCopiedToCliboardCallback( const int copied );
#endif
  
  InterSpec *m_interspec;
  DetectorDisplay *m_detector;
  
  Wt::WText *m_msg;
  Wt::WLineEdit *m_distance;
  Wt::WString m_prevDistance; // For undo/redo
  RowStretchTreeView *m_table;
  
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  Wt::WPushButton *m_copyBtn;
  Wt::JSignal<int> m_infoCopied;
#endif
  
  /** We will only update dm_data and m_uncertainties from peak model right
      before rendering happens to avoid duplicate work if multiple peaks are
      being added.
   */
  bool m_needsTableRefresh;
  
  /** What columns to show. */
  DisplayInfoLevel m_displayInfoLevel;
  Wt::WButtonGroup *m_displayLevelButtons;
  
  Wt::Signal<> m_tableUpdated;
  
  std::array<Wt::WString,FluxColumns::FluxNumColumns> m_colnames;
  std::array<Wt::WString,FluxColumns::FluxNumColumns> m_colnamesCsv;
  
  std::vector<std::string> m_nucNames;
  std::vector<std::array<double,FluxColumns::FluxNumColumns>> m_data;
  std::vector<std::array<double,FluxColumns::FluxNumColumns>> m_uncertainties;
  
  friend class FluxToolWindow;
  friend class FluxToolImp::FluxModel;
  friend class FluxToolImp::FluxCsvResource;
};//class FluxToolWidget

#endif //FluxTool_h

