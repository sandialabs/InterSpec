#ifndef DetectionLimitSimple_h
#define DetectionLimitSimple_h
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

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

class DetectionLimitSimple;

namespace Wt
{
  class WMenu;
  class WText;
  class WComboBox;
  class WLineEdit;
  class WPushButton;
  class WGridLayout;
  class WButtonGroup;
  class WStackedWidget;
  class WSuggestionPopup;
}//namespace Wt


class DetectionLimitSimpleWindow : public AuxWindow
{
public:
  DetectionLimitSimpleWindow( MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggestion,
                  InterSpec* viewer );
  
  virtual ~DetectionLimitSimpleWindow();

  DetectionLimitSimple *tool();
  
protected:
  DetectionLimitSimple *m_tool;
};//class DetectionLimitSimple


/**
 */
class DetectionLimitSimple : public Wt::WContainerWidget
{
public:
  
  DetectionLimitSimple( MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggestion,
                  InterSpec *specViewer,
                  Wt::WContainerWidget *parent = 0 );
  
  
  virtual ~DetectionLimitSimple();
  
  /** Handles receiving a "deep-link" url starting with "interspec://simple-mda/...".
   
   Example URIs:
   - "interspec://simple-mda/convoluted?ver=1&nuc=u235&energy=185&dist=100cm&..."
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://simple-mda/curie?nuc=u238...", then this string would be "nuc=u238...",
          showing the Curie-style limit.
          This string is is in standard URL format of "key1=value1&key2=value2&..." with ordering not mattering.
          Capitalization is not important.
          Assumes the string passed in has already been url-decoded.
          If not a valid path or query_str, throws exception.
   */
  void handleAppUrl( std::string path, std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string does not include the
   "interspec://" protocol, or "dose" authority; so will look something like "curie?nuc=Cs137&energy=661&dist=100cm&...",
   The path part of the URI specifies tab the tool is on.
   and it will not be url-encoded.
   */
  std::string encodeStateToUrl() const;
protected:
  
  void init();
  
  InterSpec *m_viewer;
  Wt::WSuggestionPopup *m_materialSuggest;
  MaterialDB *m_materialDB;
  // ShieldingSelect *m_enterShieldingSelect;
  
  D3SpectrumDisplayDiv *m_spectrum;
  
  /** For tracking undo/redo, we will keep track of the widgets state, as a URI.
   \sa encodeStateToUrl
   \sa handleAppUrl
   */
  std::string m_stateUri;
};//class DoseCalcWidget

#endif //DetectionLimitSimple_h
