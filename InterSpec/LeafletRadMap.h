#ifndef LeafletRadMap_h
#define LeafletRadMap_h
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
#include <functional>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

// Forward declarations
class SpecMeas;
class InterSpec;
class SimpleDialog;
class LeafletRadMapWindow;
namespace SpecUtils{ enum class SpectrumType : int; }

/*
 TODO:
   - [ ] Sometimes the initial zoom for a bunch of points is way more zoomed out than it should be
 */

class LeafletRadMap : public Wt::WContainerWidget
{
public:
  /** Shows the map dialog for the passed in measurement.
   
   Depending on the "ShowMapDataWarning" user preference, may show a warning dialog first, before the LeafletRadMapWindow
   is shown.
   
   @param meas The SpecMeas to show GPS coordinates from
   @param sample_numbers The sample numbers to show markers for on the map.  If empty, all sample numbers will be used.
   @param detector_names The names of the detectors to include for ca=ounting CPS and such. If empty, all detectors will be used.
   @param on_create Optional callback to get the pointer to the created `LeafletRadMapWindow`.  Callback is necassary as the LeafletRadMap may not be immediately created if the "ShowMapDataWarning" is true.
   @returns Pointer to the warning dialog if user preferenence for "ShowMapDataWarning" is true.  If no warning dialog is shown, then will return nullptr.
   */
  static SimpleDialog *showForMeasurement( const std::shared_ptr<const SpecMeas> meas,
                                          const std::set<int> &sample_numbers,
                                          const std::vector<std::string> &detector_names,
                                  std::function<void(LeafletRadMapWindow *)> on_create = nullptr );
  
  LeafletRadMap( Wt::WContainerWidget *parent = nullptr );
  virtual ~LeafletRadMap();
  
  void displayMeasurementOnMap( const std::shared_ptr<const SpecMeas> &meas,
                                std::set<int> sample_numbers,
                                std::vector<std::string> detector_names );
  
  
  /** Override WWebWidget::doJavaScript() to wait until this widget has been rendered before
   executing javascript so we can be sure all the JS objects we need are created.
   */
  virtual void doJavaScript( const std::string &js );
  
  
  static std::string createGeoLocationJson( const std::shared_ptr<const SpecMeas> &meas,
                                           const std::set<int> &sample_to_include,
                                           const std::vector<std::string> &detector_names,
                                           const std::set<int> &foreground_samples,
                                           const std::set<int> &background_samples,
                                           const std::set<int> &secondary_samples );
  
protected:
  void defineJavaScript();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  void handleDisplayedSpectrumChanged();
  
  void handleLoadSamples( const std::string &samples, const std::string &meas_type );
  
  std::shared_ptr<const SpecMeas> m_meas;
  
  
  /** The javascript variable name used to refer to the LeafletRadMap JS object.
      Currently is `jsRef() + ".map"`.
   */
  const std::string m_jsmap;
  
  Wt::JSignal<std::string,std::string> m_displaySamples;
  Wt::Signal<SpecUtils::SpectrumType, std::shared_ptr<const SpecMeas>, std::set<int>> m_loadSelected;
  
  
  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
     Note that not all calls to the D3 chart before Wt's rendering need to go
     here as they will be options set to the D3 chart during first rendering.
   */
  std::vector<std::string> m_pendingJs;
};//class LeafletRadMap


/** This class creates a popup window to display the LeafletRadMap. */
class LeafletRadMapWindow : public AuxWindow
{

public:
  LeafletRadMapWindow();
  virtual ~LeafletRadMapWindow();
  
  LeafletRadMap *map();
  
protected:
  LeafletRadMap *m_map;
};//class GammaXsWindow

#endif //LeafletRadMap_h
