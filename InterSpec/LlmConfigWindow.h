#ifndef LlmConfigWindow_h
#define LlmConfigWindow_h
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

#if( USE_LLM_INTERFACE )

#include <string>
#include <memory>
#include <vector>
#include <functional>

#include <Wt/WString>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/LlmConfig.h"

class InterSpec;

namespace Wt
{
  class WText;
  class WCheckBox;
  class WLineEdit;
  class WComboBox;
  class WPushButton;
  class WContainerWidget;
  // Wt::JSignal (a 6-arg template with defaults) is pulled in via AuxWindow.h.
}//namespace Wt

/** A popup window to graphically configure the LLM provider/model settings and the MCP server.

 The window edits a working copy of the full `LlmConfig` (seeded from the currently-loaded config, or
 the on-disk file, or a default template) and on "Accept" writes a documented `llm_config.xml` into
 `InterSpec::writableDataDirectory()`.  Only the GUI-exposed fields are edited; everything else in the
 working copy (e.g. DeepResearchUrl, CompactionSystemPrompt, DebugFile) is preserved so the config
 round-trips without data loss.

 After a successful save the supplied `onSaved` callback is invoked so the caller can re-read the
 config and refresh the LLM assistant.
 */
class LlmConfigWindow : public AuxWindow
{
public:
  /** @param viewer The InterSpec instance (used for i18n bundle + writable data dir).
      @param onApply Invoked when the user chooses to apply the saved config to the running session
                     (re-reads the config and refreshes the assistant).  May be empty.
      @param hasConversation Returns whether there is a non-empty conversation right now; used to
                     decide whether to warn before applying.  May be empty (treated as "no").
      @param importConfigPath If non-empty, the working copy is seeded from this on-disk
                     `llm_config.xml` (e.g. a file the user dragged onto / opened in the app) rather
                     than from the currently-loaded config.  In this "import" mode, Accept warns the
                     user before overwriting any existing config in the writable data directory.
   */
  LlmConfigWindow( InterSpec *viewer,
                   std::function<void()> onApply,
                   std::function<bool()> hasConversation,
                   const std::string &importConfigPath = std::string() );
  virtual ~LlmConfigWindow();

private:
  /** Transient per-provider UI state that is not part of the LlmConfig data model. */
  struct ProviderUiState
  {
    std::vector<std::string> discoveredModelIds; // model ids returned by the last "Fetch models"
    bool expanded = false;
    bool tokenRevealed = false;
    bool fetchInProgress = false;
    bool fetchOk = false;
    int fetchReqId = 0;          // id of the in-flight fetch request (matched in onFetchModelsResult)
    Wt::WString fetchStatusMsg;
  };//struct ProviderUiState

  void seedWorkingCopy();
  void buildUi();

  // Provider/model section (recreated wholesale on structural changes)
  void rebuildProviders();
  void buildProviderCard( size_t providerIndex );
  void buildModelEditor( Wt::WContainerWidget *parent, size_t providerIndex, size_t modelIndex );

  // Structural mutators (trigger a rebuild)
  void addProvider();
  void removeProvider( size_t providerIndex );
  void setActiveProvider( size_t providerIndex );
  void addModel( size_t providerIndex );
  void removeModel( size_t providerIndex, size_t modelIndex );
  void setActiveModel( size_t providerIndex, size_t modelIndex );

  // Lightweight refreshes (no rebuild -> no focus loss while typing)
  void updateSummary();
  void updateValidation();
  void refreshPreview();

  /** Resolved wire format for a provider: explicit apiFormat, else auto-detected from the endpoint. */
  LlmConfig::LlmApi::ApiFormat resolvedFormat( size_t providerIndex ) const;

  // "Fetch models" browser-fetch + JSignal round-trip
  void installFetchJs();
  void requestFetchModels( size_t providerIndex );
  void onFetchModelsResult( std::string payload );
  /** Derive the provider's models-listing URL (".../models") from its chat endpoint. */
  static std::string deriveModelsUrl( const std::string &chatEndpoint );

  // Preview + save
  void togglePreviewXml();
  bool saveConfig();
  void handleAccept();
  /** Save the file, optionally apply to the running session (m_onApply), then close the window. */
  void doAccept( bool applyToSession );

  static bool isPlaceholderToken( const std::string &token );

  InterSpec *m_interspec;
  std::function<void()> m_onApply;
  std::function<bool()> m_hasConversation;
  std::string m_writableDir; // empty => saving unavailable on this platform
  std::string m_importConfigPath; // non-empty => seed from this file and warn before overwriting

  LlmConfig m_working;                       // authoritative working state (llmApi + mcpServer used)
  std::vector<ProviderUiState> m_uiState;    // parallel to m_working.llmApi.providers

  // Top row
  Wt::WCheckBox *m_enableLlm;
  Wt::WText *m_activeSummary;
  Wt::WContainerWidget *m_providersArea;

  // MCP card
  Wt::WCheckBox *m_mcpEnable;
  Wt::WContainerWidget *m_mcpDetail;
#if( MCP_ENABLE_AUTH )
  Wt::WLineEdit *m_mcpToken;
  Wt::WPushButton *m_mcpTokenShow;
  Wt::WText *m_mcpTokenWarn;
#endif

  // Footer + XML preview
  Wt::WText *m_validationSummary;
  Wt::WPushButton *m_acceptBtn;
  Wt::WPushButton *m_previewToggle;
  Wt::WContainerWidget *m_previewContainer;
  Wt::WText *m_previewPre;
  bool m_showXml;

  // Fetch-models bridge: one JSON string payload {reqId, ok, status, body, error} from the browser.
  std::unique_ptr< Wt::JSignal<std::string> > m_fetchModelsResult;
  bool m_fetchJsInstalled;
  int m_nextFetchReqId;
};//class LlmConfigWindow

#endif // USE_LLM_INTERFACE
#endif // LlmConfigWindow_h
