#ifndef InjaLogDialog_h
#define InjaLogDialog_h
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

#include <tuple>
#include <vector>
#include <string>
#include <utility>
#include <functional>

#include <Wt/WString>

#include <nlohmann/json.hpp>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "InterSpec/SimpleDialog.h"


// Forward declarations
namespace Wt
{
  class WText;
  class WComboBox;
  class WResource;
  class WPushButton;
  class WContainerWidget;
}// namespace Wt


/** A reusable dialog for displaying content rendered from inja templates.

 This dialog allows displaying HTML content generated from JSON data using inja templates.
 It supports multiple template options (selected via a combo box), displays the content in
 an iframe for isolation, and provides download functionality.

 Originally extracted from BatchGuiAnaWidget to be reusable across the application.
 */
class InjaLogDialog : public SimpleDialog
{
public:
  /** Enum to specify the type of log content. */
  enum class LogType
  {
    Html,  ///< HTML content that will be displayed directly in iframe
    Text   ///< Plain text content that will be HTML-escaped and wrapped for display
  };

  /** HTML wrapper tags for text log templates.
   Use these to wrap plain text content in HTML for display,
   so this way they can be stripped off when downloading the text version.
   */
  static const char * const sm_txt_log_pre_wrapper;
  static const char * const sm_txt_log_post_wrapper;

  /** Constructor for InjaLogDialog.

   @param title The dialog title
   @param data The JSON data to be passed to inja templates
   @param template_options A vector of template options. Each entry contains a tuple:
                          - get<0>: Display name for the template (shown in combo box)
                          - get<1>: Filename suffix to append to spectrum filename for downloads (e.g., "_act_fit.html" or "_std_fit_log.txt")
                          - get<2>: LogType enum indicating whether content is Html or Text
                          - get<3>: Function that takes inja::Environment and JSON data,
                                   and returns rendered content string
                            - For Html type: content must be a valid HTML page, as it will be displayed in an iframe
                            - For Text type: content will be HTML-escaped and wrapped in <pre> tags for display

   If template_options contains only one entry, no combo box is shown.
   If template_options is empty, an error is displayed.
   */
  InjaLogDialog( const Wt::WString &title,
                 const nlohmann::json &data,
                 std::vector<std::tuple<Wt::WString, std::string, LogType, std::function<std::string(inja::Environment&, const nlohmann::json&)>>> template_options );

  virtual ~InjaLogDialog();

  const std::string &current_content() const;
  const std::string current_suggested_name() const;

  /** Show or hide the toolbar containing the template selector and download button.
   @param show If true, toolbar is shown; if false, toolbar is hidden
   */
  void setToolbarVisible( const bool show );

protected:
  /** Update the displayed content using the currently selected template. */
  void updateDisplay();

  /** Handle template selection change from combo box. */
  void handleTemplateChange();

private:
  Wt::WResource *m_download_resource;


  // Member variables
  nlohmann::json m_data;
  std::vector<std::tuple<Wt::WString, std::string, LogType, std::function<std::string(inja::Environment&, const nlohmann::json&)>>> m_template_options;

  Wt::WContainerWidget *m_toolbar;
  Wt::WComboBox *m_template_selector;
  Wt::WPushButton *m_download_button;
  Wt::WText *m_iframe_holder;
  Wt::WResource *m_current_resource;
  std::string m_current_content;
};//class InjaLogDialog


#endif // InjaLogDialog_h
