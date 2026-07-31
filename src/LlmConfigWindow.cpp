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

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <functional>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WString>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WWebWidget>
#include <Wt/WJavaScript>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include <nlohmann/json.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmApiProtocol.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/LlmConfigWindow.h"

using namespace std;
using namespace Wt;

using ApiFormat = LlmConfig::LlmApi::ApiFormat;
using ReasoningEffort = LlmConfig::LlmApi::ReasoningEffort;
using InstructionVerbosity = LlmConfig::LlmApi::InstructionVerbosity;
using ApiProvider = LlmConfig::LlmApi::ApiProvider;
using ModelInfo = LlmConfig::LlmApi::ModelInfo;

namespace
{
  /** Extract the host portion of a URL for display (best-effort, no exceptions). */
  string host_from_url( const string &url )
  {
    string u = url;
    SpecUtils::trim( u );
    if( u.empty() )
      return "";
    const size_t scheme = u.find( "://" );
    if( scheme != string::npos )
      u = u.substr( scheme + 3 );
    const size_t slash = u.find( '/' );
    if( slash != string::npos )
      u = u.substr( 0, slash );
    return u;
  }//host_from_url

  WString format_label( ApiFormat fmt )
  {
    switch( fmt )
    {
      case ApiFormat::OpenAiChat:      return WString::tr("lcw-fmt-openai-chat");
      case ApiFormat::OpenAiResponses: return WString::tr("lcw-fmt-openai-responses");
      case ApiFormat::Anthropic:       return WString::tr("lcw-fmt-anthropic");
    }
    return WString::tr("lcw-fmt-openai-chat");
  }//format_label

  /** Reasoning variant -> combo index: 0=off, 1=on(bool), 2=low, 3=medium, 4=high. */
  int reasoning_to_index( const std::variant<bool,ReasoningEffort> &r )
  {
    if( std::holds_alternative<bool>(r) )
      return std::get<bool>(r) ? 1 : 0;
    switch( std::get<ReasoningEffort>(r) )
    {
      case ReasoningEffort::low:    return 2;
      case ReasoningEffort::medium: return 3;
      case ReasoningEffort::high:   return 4;
    }
    return 0;
  }//reasoning_to_index

  std::variant<bool,ReasoningEffort> index_to_reasoning( int idx )
  {
    switch( idx )
    {
      case 0: return false;
      case 1: return true;
      case 2: return ReasoningEffort::low;
      case 3: return ReasoningEffort::medium;
      case 4: return ReasoningEffort::high;
    }
    return false;
  }//index_to_reasoning

  /** InstructionVerbosity -> combo index: 0=terse, 1=normal, 2=verbose. */
  int verbosity_to_index( const InstructionVerbosity v )
  {
    switch( v )
    {
      case InstructionVerbosity::Terse:   return 0;
      case InstructionVerbosity::Normal:  return 1;
      case InstructionVerbosity::Verbose: return 2;
    }
    return 1;
  }//verbosity_to_index

  InstructionVerbosity index_to_verbosity( int idx )
  {
    switch( idx )
    {
      case 0: return InstructionVerbosity::Terse;
      case 1: return InstructionVerbosity::Normal;
      case 2: return InstructionVerbosity::Verbose;
    }
    return InstructionVerbosity::Normal;
  }//index_to_verbosity
}//namespace


LlmConfigWindow::LlmConfigWindow( InterSpec *viewer,
                                  std::function<void()> onApply,
                                  std::function<bool()> hasConversation,
                                  const std::string &importConfigPath )
  : AuxWindow( WString::tr("lcw-window-title"),
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                | AuxWindowProperties::SetCloseable
                | AuxWindowProperties::DisableCollapse) ),
    m_interspec( viewer ),
    m_onApply( std::move(onApply) ),
    m_hasConversation( std::move(hasConversation) ),
    m_importConfigPath( importConfigPath ),
    m_enableLlm( nullptr ),
    m_activeSummary( nullptr ),
    m_providersArea( nullptr ),
    m_mcpEnable( nullptr ),
    m_mcpDetail( nullptr ),
#if( MCP_ENABLE_AUTH )
    m_mcpToken( nullptr ),
    m_mcpTokenShow( nullptr ),
    m_mcpTokenWarn( nullptr ),
#endif
    m_validationSummary( nullptr ),
    m_acceptBtn( nullptr ),
    m_previewToggle( nullptr ),
    m_previewContainer( nullptr ),
    m_previewPre( nullptr ),
    m_showXml( false ),
    m_fetchJsInstalled( false ),
    m_nextFetchReqId( 1 )
{
  if( !m_interspec )
    throw runtime_error( "LlmConfigWindow: null InterSpec" );

  m_interspec->useMessageResourceBundle( "LlmConfigWindow" );
  wApp->useStyleSheet( "InterSpec_resources/LlmConfigWindow.css" );

  try
  {
    m_writableDir = InterSpec::writableDataDirectory();
  }catch( std::exception & )
  {
    m_writableDir.clear();
  }

  m_fetchModelsResult.reset( new JSignal<std::string>( this, "llmFetchModelsResult" ) );
  m_fetchModelsResult->connect( this, &LlmConfigWindow::onFetchModelsResult );

  seedWorkingCopy();
  buildUi();
  rebuildProviders();
  updateSummary();
  updateValidation();

  finished().connect( std::bind( &AuxWindow::deleteAuxWindow, this ) );

  const int w = m_interspec->renderedWidth();
  const int h = m_interspec->renderedHeight();
  if( (w > 100) && (h > 100) )
  {
    const int ww = std::min( 900, (95*w)/100 );
    const int wh = std::min( std::max( 480, (85*h)/100 ), 820 );
    resizeWindow( ww, wh );
  }else
  {
    resizeScaledWindow( 0.85, 0.9 );
  }

  rejectWhenEscapePressed();
  centerWindow();
  show();
}//LlmConfigWindow constructor


LlmConfigWindow::~LlmConfigWindow()
{
}


void LlmConfigWindow::seedWorkingCopy()
{
  // 0) Import mode: seed from the file the user opened/dropped, so the GUI shows that config.
  if( !m_importConfigPath.empty() )
  {
    try
    {
      std::pair<LlmConfig::LlmApi,LlmConfig::McpServer> parsed
                                    = LlmConfig::loadApiAndMcpConfigs( m_importConfigPath );
      m_working.llmApi = std::move( parsed.first );
      m_working.mcpServer = std::move( parsed.second );
      m_uiState.assign( m_working.llmApi.providers.size(), ProviderUiState() );
      return;
    }catch( std::exception & )
    {
      // The caller validates the file parses before opening this window, so this should not happen;
      // fall through to the normal seeding tiers rather than showing an empty window.
    }
  }//if( import mode )

  // 1) The in-memory (cached) config is the truth that the assistant uses; prefer it.
  try
  {
    shared_ptr<const LlmConfig> cached = InterSpecServer::llm_config();
    if( cached )
    {
      m_working.llmApi = cached->llmApi;
      m_working.mcpServer = cached->mcpServer;
      m_uiState.assign( m_working.llmApi.providers.size(), ProviderUiState() );
      return;
    }
  }catch( std::exception & )
  {
    // Fall through: no usable cached config (file missing/invalid).
  }

  // 2) Parse the on-disk file directly (writable dir first, then static).
  const char * const fname = "llm_config.xml";
  vector<string> candidates;
  if( !m_writableDir.empty() )
    candidates.push_back( SpecUtils::append_path( m_writableDir, fname ) );
  try
  {
    candidates.push_back( SpecUtils::append_path( InterSpec::staticDataDirectory(), fname ) );
  }catch( std::exception & ){ }

  for( const string &path : candidates )
  {
    if( !SpecUtils::is_file(path) )
      continue;
    try
    {
      std::pair<LlmConfig::LlmApi,LlmConfig::McpServer> parsed = LlmConfig::loadApiAndMcpConfigs( path );
      m_working.llmApi = std::move( parsed.first );
      m_working.mcpServer = std::move( parsed.second );
      m_uiState.assign( m_working.llmApi.providers.size(), ProviderUiState() );
      return;
    }catch( std::exception & )
    {
      // try next candidate
    }
  }//for( candidates )

  // 3) Nothing usable on disk: seed a minimal, enabled, single-provider template.
  m_working.llmApi = LlmConfig::LlmApi();
  m_working.llmApi.enabled = true;
  ApiProvider prov;
  prov.models.push_back( ModelInfo() );
  m_working.llmApi.providers.push_back( std::move(prov) );
  m_working.llmApi.activeProviderIndex = 0;
  m_working.mcpServer = LlmConfig::McpServer();
  m_uiState.assign( 1, ProviderUiState() );
}//seedWorkingCopy()


void LlmConfigWindow::buildUi()
{
  WContainerWidget *body = contents();
  body->addStyleClass( "LlmConfigBody" );
  body->setOverflow( WContainerWidget::OverflowAuto, Wt::Vertical );

  // --- Enable LLM row ---
  WContainerWidget *enableRow = new WContainerWidget( body );
  enableRow->addStyleClass( "LcwCard LcwEnableRow" );

  m_enableLlm = new WCheckBox( WString(), enableRow );
  m_enableLlm->addStyleClass( "LcwCheck" );
  m_enableLlm->setChecked( m_working.llmApi.enabled );
  m_enableLlm->changed().connect( std::bind( [this](){
    m_working.llmApi.enabled = m_enableLlm->isChecked();
    updateValidation();
    refreshPreview();
  } ) );

  WContainerWidget *enableText = new WContainerWidget( enableRow );
  enableText->addStyleClass( "LcwLabelStack" );
  WText *enTitle = new WText( WString::tr("lcw-enable-llm"), enableText );
  enTitle->addStyleClass( "LcwTitle" );
  WText *enSub = new WText( WString::tr("lcw-enable-llm-sub"), enableText );
  enSub->addStyleClass( "LcwSubtle" );

  m_activeSummary = new WText( enableRow );
  m_activeSummary->addStyleClass( "LcwActiveSummary" );

  // --- Providers header ---
  WContainerWidget *provHeader = new WContainerWidget( body );
  provHeader->addStyleClass( "LcwSectionHeader" );
  WText *provLbl = new WText( WString::tr("lcw-providers"), provHeader );
  provLbl->addStyleClass( "LcwSectionLabel" );
  WPushButton *addProv = new WPushButton( WString::tr("lcw-add-provider"), provHeader );
  addProv->addStyleClass( "LinkBtn" );  // app-standard link-style button (flat, no hover chrome)
  addProv->clicked().connect( this, &LlmConfigWindow::addProvider );

  // --- Providers list (recreated by rebuildProviders) ---
  m_providersArea = new WContainerWidget( body );
  m_providersArea->addStyleClass( "LcwProviders" );

  // --- MCP card ---
  WContainerWidget *mcpCard = new WContainerWidget( body );
  mcpCard->addStyleClass( "LcwCard LcwMcpCard" );

  WContainerWidget *mcpTop = new WContainerWidget( mcpCard );
  mcpTop->addStyleClass( "LcwEnableRow" );
  m_mcpEnable = new WCheckBox( WString(), mcpTop );
  m_mcpEnable->addStyleClass( "LcwCheck" );
  m_mcpEnable->setChecked( m_working.mcpServer.enabled );
  WContainerWidget *mcpText = new WContainerWidget( mcpTop );
  mcpText->addStyleClass( "LcwLabelStack" );
  WText *mcpTitle = new WText( WString::tr("lcw-mcp-server"), mcpText );
  mcpTitle->addStyleClass( "LcwTitle" );
  WText *mcpSub = new WText( WString::tr("lcw-mcp-sub"), mcpText );
  mcpSub->addStyleClass( "LcwSubtle" );

  m_mcpDetail = new WContainerWidget( mcpCard );
  m_mcpDetail->addStyleClass( "LcwMcpDetail" );

#if( MCP_ENABLE_AUTH )
  WLabel *mcpTokLbl = new WLabel( WString::tr("lcw-mcp-token"), m_mcpDetail );
  mcpTokLbl->addStyleClass( "LcwFieldLabel" );
  WContainerWidget *mcpTokRow = new WContainerWidget( m_mcpDetail );
  mcpTokRow->addStyleClass( "LcwTokenRow" );
  m_mcpToken = new WLineEdit( mcpTokRow );
  m_mcpToken->addStyleClass( "LcwInput LcwMono" );
  m_mcpToken->setEchoMode( WLineEdit::Password );
  m_mcpToken->setText( WString::fromUTF8( m_working.mcpServer.bearerToken ) );
  m_mcpToken->changed().connect( std::bind( [this](){
    m_working.mcpServer.bearerToken = m_mcpToken->text().toUTF8();
    updateValidation();
    refreshPreview();
  } ) );
  m_mcpTokenShow = new WPushButton( WString::tr("lcw-show"), mcpTokRow );
  m_mcpTokenShow->clicked().connect( std::bind( [this](){
    const bool reveal = (m_mcpToken->echoMode() == WLineEdit::Password);
    m_mcpToken->setEchoMode( reveal ? WLineEdit::Normal : WLineEdit::Password );
    m_mcpTokenShow->setText( reveal ? WString::tr("lcw-hide") : WString::tr("lcw-show") );
  } ) );
  m_mcpTokenWarn = new WText( WString::tr("lcw-mcp-token-invalid"), m_mcpDetail );
  m_mcpTokenWarn->addStyleClass( "LcwError" );
#endif

  m_mcpEnable->changed().connect( std::bind( [this](){
    m_working.mcpServer.enabled = m_mcpEnable->isChecked();
    m_mcpDetail->setHidden( !m_working.mcpServer.enabled );
    updateValidation();
    refreshPreview();
  } ) );
  m_mcpDetail->setHidden( !m_working.mcpServer.enabled );

  // --- XML preview (hidden until toggled) ---
  m_previewContainer = new WContainerWidget( body );
  m_previewContainer->addStyleClass( "LcwXmlPreview" );
  m_previewPre = new WText( m_previewContainer );
  m_previewPre->setTextFormat( Wt::PlainText );
  m_previewPre->addStyleClass( "LcwXmlPre" );
  m_previewContainer->setHidden( true );

  // --- Footer ---
  WContainerWidget *foot = footer();

  WContainerWidget *footLeft = new WContainerWidget( foot );
  footLeft->addStyleClass( "LcwFooterLeft" );
  m_validationSummary = new WText( footLeft );
  m_validationSummary->addStyleClass( "LcwValidation" );
  m_previewToggle = new WPushButton( WString::tr("lcw-preview-xml"), footLeft );
  m_previewToggle->addStyleClass( "LinkBtn" );
  m_previewToggle->clicked().connect( this, &LlmConfigWindow::togglePreviewXml );

  m_acceptBtn = new WPushButton( WString::tr("lcw-accept"), foot );
  m_acceptBtn->clicked().connect( this, &LlmConfigWindow::handleAccept );

  WPushButton *cancelBtn = new WPushButton( WString::tr("lcw-cancel"), foot );
  cancelBtn->clicked().connect( this, &AuxWindow::hide );

  if( m_writableDir.empty() )
  {
    m_acceptBtn->disable();
    m_acceptBtn->setToolTip( WString::tr("lcw-save-unavailable") );
  }
}//buildUi()


void LlmConfigWindow::rebuildProviders()
{
  if( !m_providersArea )
    return;

  m_providersArea->clear();

  // Keep the UI-state vector in lock-step with the providers.
  if( m_uiState.size() != m_working.llmApi.providers.size() )
    m_uiState.resize( m_working.llmApi.providers.size() );

  for( size_t pi = 0; pi < m_working.llmApi.providers.size(); ++pi )
    buildProviderCard( pi );
}//rebuildProviders()


void LlmConfigWindow::buildProviderCard( const size_t pi )
{
  const ApiProvider &prov = m_working.llmApi.providers[pi];
  ProviderUiState &ui = m_uiState[pi];
  const ApiFormat fmt = resolvedFormat( pi );
  const bool isActive = (pi == m_working.llmApi.activeProviderIndex);

  WContainerWidget *card = new WContainerWidget( m_providersArea );
  card->addStyleClass( "LcwProviderCard" );
  if( isActive )
    card->addStyleClass( "LcwActive" );

  // ---- Header (click to expand) ----
  WContainerWidget *header = new WContainerWidget( card );
  header->addStyleClass( "LcwProviderHeader" );
  header->clicked().connect( std::bind( [this,pi](){
    if( pi < m_uiState.size() )
      m_uiState[pi].expanded = !m_uiState[pi].expanded;
    rebuildProviders();
  } ) );

  WText *chev = new WText( ui.expanded ? "\xE2\x96\xBC" : "\xE2\x96\xB6", header ); // ▼ / ▶
  chev->addStyleClass( "LcwChevron" );

  WContainerWidget *titleStack = new WContainerWidget( header );
  titleStack->addStyleClass( "LcwProviderTitleStack" );
  const string host = host_from_url( prov.apiEndpoint );
  WText *hostLbl = new WText( host.empty() ? WString::tr("lcw-new-provider")
                                           : WString::fromUTF8(host), titleStack );
  hostLbl->addStyleClass( "LcwProviderHost" );

  const string activeModelName = (prov.activeModelIndex < prov.models.size())
                                   ? prov.models[prov.activeModelIndex].name : string();
  WString sub = WString::tr("lcw-provider-subtitle")
                  .arg( format_label(fmt) )
                  .arg( static_cast<int>(prov.models.size()) )
                  .arg( activeModelName.empty() ? WString::tr("lcw-active-none")
                                                : WString::fromUTF8(activeModelName) );
  WText *subLbl = new WText( sub, titleStack );
  subLbl->addStyleClass( "LcwProviderSub" );

  // Active checkbox (exclusive across providers).  Stop the header's expand-click from also firing.
  WContainerWidget *activeWrap = new WContainerWidget( header );
  activeWrap->addStyleClass( "LcwHeaderActive" );
  activeWrap->clicked().preventPropagation();
  WCheckBox *activeCb = new WCheckBox( WString::tr("lcw-active-label"), activeWrap );
  activeCb->setChecked( isActive );
  activeCb->changed().connect( std::bind( [this,pi](){ setActiveProvider(pi); } ) );

  WPushButton *removeBtn = new WPushButton( header );
  removeBtn->setStyleClass( "LcwRemoveIcon Wt-icon" );
  removeBtn->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  removeBtn->setToolTip( WString::tr("lcw-remove-provider") );
  removeBtn->clicked().preventPropagation();
  removeBtn->clicked().connect( std::bind( [this,pi](){ removeProvider(pi); } ) );

  if( !ui.expanded )
    return;

  // ---- Body ----
  WContainerWidget *bodyc = new WContainerWidget( card );
  bodyc->addStyleClass( "LcwProviderBody" );

  // Endpoint
  WLabel *epLbl = new WLabel( WString::tr("lcw-api-endpoint"), bodyc );
  epLbl->addStyleClass( "LcwFieldLabel" );
  WLineEdit *epEdit = new WLineEdit( bodyc );
  epEdit->addStyleClass( "LcwInput" );
  epEdit->setText( WString::fromUTF8(prov.apiEndpoint) );
  epEdit->setPlaceholderText( WString::tr("lcw-endpoint-ph") );
  epEdit->changed().connect( std::bind( [this,pi,epEdit](){
    if( pi >= m_working.llmApi.providers.size() )
      return;
    m_working.llmApi.providers[pi].apiEndpoint = epEdit->text().toUTF8();
    rebuildProviders();  // refresh format badge / host / anthropic warnings
    updateSummary();
    updateValidation();
    refreshPreview();
  } ) );

  // Wire format row
  WContainerWidget *fmtRow = new WContainerWidget( bodyc );
  fmtRow->addStyleClass( "LcwFormatRow" );
  WText *fmtLblTxt = new WText( WString::tr("lcw-wire-format"), fmtRow );
  fmtLblTxt->addStyleClass( "LcwFieldLabelInline" );
  WText *fmtBadge = new WText( format_label(fmt), fmtRow );
  fmtBadge->addStyleClass( "LcwBadge" );
  WText *fmtSrc = new WText( prov.apiFormat.has_value() ? WString::tr("lcw-fmt-source-manual")
                                                        : WString::tr("lcw-fmt-source-auto"), fmtRow );
  fmtSrc->addStyleClass( "LcwSubtle" );

  WContainerWidget *ovrWrap = new WContainerWidget( fmtRow );
  ovrWrap->addStyleClass( "LcwOverrideWrap" );
  new WText( WString::tr("lcw-override"), ovrWrap );
  WComboBox *ovr = new WComboBox( ovrWrap );
  ovr->addStyleClass( "LcwSelect" );
  ovr->addItem( WString::tr("lcw-fmt-auto") );             // 0 -> nullopt
  ovr->addItem( WString::tr("lcw-fmt-openai-chat") );      // 1 -> OpenAiChat
  ovr->addItem( WString::tr("lcw-fmt-openai-responses") ); // 2 -> OpenAiResponses
  ovr->addItem( WString::tr("lcw-fmt-anthropic") );        // 3 -> Anthropic
  int ovrIdx = 0;
  if( prov.apiFormat.has_value() )
  {
    switch( prov.apiFormat.value() )
    {
      case ApiFormat::OpenAiChat:      ovrIdx = 1; break;
      case ApiFormat::OpenAiResponses: ovrIdx = 2; break;
      case ApiFormat::Anthropic:       ovrIdx = 3; break;
    }
  }
  ovr->setCurrentIndex( ovrIdx );
  ovr->activated().connect( std::bind( [this,pi,ovr](){
    if( pi >= m_working.llmApi.providers.size() )
      return;
    ApiProvider &p = m_working.llmApi.providers[pi];
    switch( ovr->currentIndex() )
    {
      case 1:  p.apiFormat = ApiFormat::OpenAiChat;      break;
      case 2:  p.apiFormat = ApiFormat::OpenAiResponses; break;
      case 3:  p.apiFormat = ApiFormat::Anthropic;       break;
      default: p.apiFormat.reset();                      break;
    }
    rebuildProviders();
    updateSummary();
    updateValidation();
    refreshPreview();
  } ) );

  // Bearer token row
  WLabel *tokLbl = new WLabel( WString::tr("lcw-bearer-token"), bodyc );
  tokLbl->addStyleClass( "LcwFieldLabel" );
  WContainerWidget *tokRow = new WContainerWidget( bodyc );
  tokRow->addStyleClass( "LcwTokenRow" );
  WLineEdit *tokEdit = new WLineEdit( tokRow );
  tokEdit->addStyleClass( "LcwInput LcwMono" );
  tokEdit->setEchoMode( ui.tokenRevealed ? WLineEdit::Normal : WLineEdit::Password );
  tokEdit->setText( WString::fromUTF8(prov.bearerToken) );
  tokEdit->setPlaceholderText( WString::tr("lcw-token-ph") );
  tokEdit->changed().connect( std::bind( [this,pi,tokEdit](){
    if( pi >= m_working.llmApi.providers.size() )
      return;
    m_working.llmApi.providers[pi].bearerToken = tokEdit->text().toUTF8();
    updateValidation();
    refreshPreview();
  } ) );
  WPushButton *revealBtn = new WPushButton( ui.tokenRevealed ? WString::tr("lcw-hide")
                                                             : WString::tr("lcw-show"), tokRow );
  revealBtn->clicked().connect( std::bind( [this,pi,tokEdit,revealBtn](){
    if( pi >= m_uiState.size() )
      return;
    m_uiState[pi].tokenRevealed = !m_uiState[pi].tokenRevealed;
    tokEdit->setEchoMode( m_uiState[pi].tokenRevealed ? WLineEdit::Normal : WLineEdit::Password );
    revealBtn->setText( m_uiState[pi].tokenRevealed ? WString::tr("lcw-hide") : WString::tr("lcw-show") );
  } ) );
  WPushButton *fetchBtn = new WPushButton( ui.fetchInProgress ? WString::tr("lcw-fetching")
                                                              : WString::tr("lcw-fetch-models"), tokRow );
  if( ui.fetchInProgress )
    fetchBtn->disable();
  fetchBtn->clicked().connect( std::bind( [this,pi](){ requestFetchModels(pi); } ) );

  if( isPlaceholderToken( prov.bearerToken ) )
  {
    WText *tokWarn = new WText( WString::tr("lcw-token-warn"), bodyc );
    tokWarn->addStyleClass( "LcwWarn" );
  }
  if( !ui.fetchStatusMsg.empty() )
  {
    WText *fetchMsg = new WText( ui.fetchStatusMsg, bodyc );
    fetchMsg->addStyleClass( ui.fetchOk ? "LcwFetchOk" : "LcwFetchErr" );
  }

  // Models subsection
  WContainerWidget *modelsHeader = new WContainerWidget( bodyc );
  modelsHeader->addStyleClass( "LcwSectionHeader LcwModelsHeader" );
  WText *modelsLbl = new WText( WString::tr("lcw-models")
                                  .arg( static_cast<int>(prov.models.size()) ), modelsHeader );
  modelsLbl->addStyleClass( "LcwSectionLabel" );
  WPushButton *addModelBtn = new WPushButton( WString::tr("lcw-add-model"), modelsHeader );
  addModelBtn->addStyleClass( "LinkBtn" );
  addModelBtn->clicked().connect( std::bind( [this,pi](){ addModel(pi); } ) );

  if( prov.models.empty() )
  {
    WText *noModels = new WText( WString::tr("lcw-no-models"), bodyc );
    noModels->addStyleClass( "LcwNoModels" );
  }else
  {
    WContainerWidget *modelsList = new WContainerWidget( bodyc );
    modelsList->addStyleClass( "LcwModelsList" );
    for( size_t mi = 0; mi < prov.models.size(); ++mi )
      buildModelEditor( modelsList, pi, mi );
  }
}//buildProviderCard()


void LlmConfigWindow::buildModelEditor( WContainerWidget *parent, const size_t pi, const size_t mi )
{
  const ApiProvider &prov = m_working.llmApi.providers[pi];
  const ModelInfo &m = prov.models[mi];
  const ProviderUiState &ui = m_uiState[pi];
  const ApiFormat fmt = resolvedFormat( pi );
  const bool isActive = (mi == prov.activeModelIndex);

  WContainerWidget *card = new WContainerWidget( parent );
  card->addStyleClass( "LcwModelCard" );
  if( isActive )
    card->addStyleClass( "LcwActive" );

  // Row 1: active + name + remove
  WContainerWidget *row1 = new WContainerWidget( card );
  row1->addStyleClass( "LcwModelRow1" );

  WCheckBox *activeCb = new WCheckBox( WString::tr("lcw-active-label"), row1 );
  activeCb->setChecked( isActive );
  activeCb->changed().connect( std::bind( [this,pi,mi](){ setActiveModel(pi,mi); } ) );

  WContainerWidget *nameWrap = new WContainerWidget( row1 );
  nameWrap->addStyleClass( "LcwModelNameWrap" );

  // Quick-pick of fetched model ids (optional convenience), then the free-text name field.
  if( !ui.discoveredModelIds.empty() )
  {
    WComboBox *pick = new WComboBox( nameWrap );
    pick->addStyleClass( "LcwSelect" );
    pick->addItem( WString::tr("lcw-select-model") );
    for( const string &id : ui.discoveredModelIds )
      pick->addItem( WString::fromUTF8(id) );
    pick->setCurrentIndex( 0 );
    pick->activated().connect( std::bind( [this,pi,mi,pick](){
      const int idx = pick->currentIndex();
      if( idx <= 0 || pi >= m_working.llmApi.providers.size() )
        return;
      ApiProvider &p = m_working.llmApi.providers[pi];
      if( mi >= p.models.size() || (idx-1) >= static_cast<int>(m_uiState[pi].discoveredModelIds.size()) )
        return;
      p.models[mi].name = m_uiState[pi].discoveredModelIds[idx-1];
      rebuildProviders();
      updateSummary();
      updateValidation();
      refreshPreview();
    } ) );
  }

  WLineEdit *nameEdit = new WLineEdit( nameWrap );
  nameEdit->addStyleClass( "LcwInput LcwMono" );
  nameEdit->setText( WString::fromUTF8(m.name) );
  nameEdit->setPlaceholderText( WString::tr("lcw-model-name-ph") );
  nameEdit->changed().connect( std::bind( [this,pi,mi,nameEdit](){
    if( pi >= m_working.llmApi.providers.size() || mi >= m_working.llmApi.providers[pi].models.size() )
      return;
    m_working.llmApi.providers[pi].models[mi].name = nameEdit->text().toUTF8();
    rebuildProviders();  // refresh subtitle / summary
    updateSummary();
    updateValidation();
    refreshPreview();
  } ) );

  WPushButton *removeBtn = new WPushButton( row1 );
  removeBtn->setStyleClass( "LcwRemoveIcon Wt-icon" );
  removeBtn->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  removeBtn->setToolTip( WString::tr("lcw-remove-model") );
  removeBtn->clicked().connect( std::bind( [this,pi,mi](){ removeModel(pi,mi); } ) );

  // Row 2: max tokens / context limit / temperature
  WContainerWidget *grid = new WContainerWidget( card );
  grid->addStyleClass( "LcwNumGrid" );

  auto makeNumField = [this,grid]( const WString &label, const WString &ph, const string &val,
                                   std::function<void(WLineEdit *)> onChange ) -> WLineEdit *
  {
    WContainerWidget *cell = new WContainerWidget( grid );
    cell->addStyleClass( "LcwNumCell" );
    new WText( label, cell );
    WLineEdit *edit = new WLineEdit( cell );
    edit->addStyleClass( "LcwNumInput" );
    edit->setPlaceholderText( ph );
    edit->setText( WString::fromUTF8(val) );
    edit->changed().connect( std::bind( [edit,onChange](){ onChange(edit); } ) );
    return edit;
  };

  makeNumField( WString::tr("lcw-max-tokens"), WString::tr("lcw-api-default"),
    (m.maxTokens > 0 ? std::to_string(m.maxTokens) : string()),
    [this,pi,mi]( WLineEdit *e ){
      if( pi >= m_working.llmApi.providers.size() || mi >= m_working.llmApi.providers[pi].models.size() )
        return;
      ModelInfo &mm = m_working.llmApi.providers[pi].models[mi];
      const bool wasAnthroErr = (resolvedFormat(pi) == ApiFormat::Anthropic) && (mm.maxTokens <= 0);
      string s = e->text().toUTF8(); SpecUtils::trim(s);
      int v = 0;
      mm.maxTokens = (!s.empty() && SpecUtils::parse_int(s.c_str(), s.size(), v) && v > 0) ? v : 0;
      const bool nowAnthroErr = (resolvedFormat(pi) == ApiFormat::Anthropic) && (mm.maxTokens <= 0);
      updateValidation();
      refreshPreview();
      // The per-model "Anthropic requires Max tokens" inline warning is only (re)drawn on a rebuild;
      // rebuild only when that state actually flips so typing in other fields isn't disrupted.
      if( wasAnthroErr != nowAnthroErr )
        rebuildProviders();
    } );

  makeNumField( WString::tr("lcw-context-limit"), WString::tr("lcw-api-default"),
    (m.contextLengthLimit > 0 ? std::to_string(m.contextLengthLimit) : string()),
    [this,pi,mi]( WLineEdit *e ){
      if( pi >= m_working.llmApi.providers.size() || mi >= m_working.llmApi.providers[pi].models.size() )
        return;
      string s = e->text().toUTF8(); SpecUtils::trim(s);
      int v = 0;
      m_working.llmApi.providers[pi].models[mi].contextLengthLimit
        = (!s.empty() && SpecUtils::parse_int(s.c_str(), s.size(), v) && v > 0) ? v : 0;
      refreshPreview();
    } );

  makeNumField( WString::tr("lcw-temperature"), WString::tr("lcw-default"),
    (m.temperature.has_value() ? SpecUtils::printCompact(m.temperature.value(), 6) : string()),
    [this,pi,mi]( WLineEdit *e ){
      if( pi >= m_working.llmApi.providers.size() || mi >= m_working.llmApi.providers[pi].models.size() )
        return;
      string s = e->text().toUTF8(); SpecUtils::trim(s);
      double v = 0.0;
      if( !s.empty() && SpecUtils::parse_double(s.c_str(), s.size(), v) )
        m_working.llmApi.providers[pi].models[mi].temperature = std::min( 2.0, std::max(0.0, v) );
      else
        m_working.llmApi.providers[pi].models[mi].temperature.reset();
      refreshPreview();
    } );

  // Row 3: supports images / force tool / reasoning / verbosity
  WContainerWidget *row3 = new WContainerWidget( card );
  row3->addStyleClass( "LcwModelRow3" );

  WCheckBox *imgCb = new WCheckBox( WString::tr("lcw-supports-images"), row3 );
  imgCb->setChecked( m.supportsImages );
  imgCb->changed().connect( std::bind( [this,pi,mi,imgCb](){
    if( pi < m_working.llmApi.providers.size() && mi < m_working.llmApi.providers[pi].models.size() )
    {
      m_working.llmApi.providers[pi].models[mi].supportsImages = imgCb->isChecked();
      refreshPreview();
    }
  } ) );

  WCheckBox *toolCb = new WCheckBox( WString::tr("lcw-force-tool"), row3 );
  toolCb->setChecked( m.requireToolInStateMachine );
  toolCb->setToolTip( WString::tr("lcw-force-tool-tip") );
  toolCb->changed().connect( std::bind( [this,pi,mi,toolCb](){
    if( pi < m_working.llmApi.providers.size() && mi < m_working.llmApi.providers[pi].models.size() )
    {
      m_working.llmApi.providers[pi].models[mi].requireToolInStateMachine = toolCb->isChecked();
      refreshPreview();
    }
  } ) );

  WContainerWidget *reasonWrap = new WContainerWidget( row3 );
  reasonWrap->addStyleClass( "LcwReasonWrap" );
  new WText( WString::tr("lcw-reasoning"), reasonWrap );
  WComboBox *reason = new WComboBox( reasonWrap );
  reason->addStyleClass( "LcwSelect" );
  reason->addItem( WString::tr("lcw-reasoning-off") );
  reason->addItem( WString::tr("lcw-reasoning-on") );
  reason->addItem( WString::tr("lcw-reasoning-low") );
  reason->addItem( WString::tr("lcw-reasoning-medium") );
  reason->addItem( WString::tr("lcw-reasoning-high") );
  reason->setCurrentIndex( reasoning_to_index( m.reasoning ) );
  reason->activated().connect( std::bind( [this,pi,mi,reason](){
    if( pi < m_working.llmApi.providers.size() && mi < m_working.llmApi.providers[pi].models.size() )
    {
      m_working.llmApi.providers[pi].models[mi].reasoning = index_to_reasoning( reason->currentIndex() );
      refreshPreview();
    }
  } ) );

  WContainerWidget *verbWrap = new WContainerWidget( row3 );
  verbWrap->addStyleClass( "LcwReasonWrap" );
  new WText( WString::tr("lcw-verbosity"), verbWrap );
  WComboBox *verbosity = new WComboBox( verbWrap );
  verbosity->addStyleClass( "LcwSelect" );
  verbosity->addItem( WString::tr("lcw-verbosity-terse") );
  verbosity->addItem( WString::tr("lcw-verbosity-normal") );
  verbosity->addItem( WString::tr("lcw-verbosity-verbose") );
  verbosity->setCurrentIndex( verbosity_to_index( m.instructionVerbosity ) );
  verbosity->setToolTip( WString::tr("lcw-verbosity-tip") );
  verbosity->activated().connect( std::bind( [this,pi,mi,verbosity](){
    if( pi < m_working.llmApi.providers.size() && mi < m_working.llmApi.providers[pi].models.size() )
    {
      m_working.llmApi.providers[pi].models[mi].instructionVerbosity = index_to_verbosity( verbosity->currentIndex() );
      refreshPreview();
    }
  } ) );

  // Anthropic requires MaxTokens
  if( (fmt == ApiFormat::Anthropic) && (m.maxTokens <= 0) )
  {
    WText *needMax = new WText( WString::tr("lcw-anthropic-needs-max"), card );
    needMax->addStyleClass( "LcwError LcwModelError" );
  }
}//buildModelEditor()


LlmConfig::LlmApi::ApiFormat LlmConfigWindow::resolvedFormat( const size_t pi ) const
{
  const ApiProvider &p = m_working.llmApi.providers[pi];
  return p.apiFormat.value_or( LlmConfig::LlmApi::detectApiFormat( p.apiEndpoint ) );
}//resolvedFormat()


void LlmConfigWindow::addProvider()
{
  ApiProvider prov;
  prov.models.push_back( ModelInfo() );
  m_working.llmApi.providers.push_back( std::move(prov) );
  m_uiState.push_back( ProviderUiState() );
  m_uiState.back().expanded = true;
  rebuildProviders();
  updateSummary();
  updateValidation();
  refreshPreview();
}//addProvider()


void LlmConfigWindow::removeProvider( const size_t pi )
{
  if( pi >= m_working.llmApi.providers.size() )
    return;
  m_working.llmApi.providers.erase( m_working.llmApi.providers.begin() + pi );
  if( pi < m_uiState.size() )
    m_uiState.erase( m_uiState.begin() + pi );

  // Keep activeProviderIndex valid.
  if( m_working.llmApi.providers.empty() )
    m_working.llmApi.activeProviderIndex = 0;
  else if( m_working.llmApi.activeProviderIndex >= m_working.llmApi.providers.size() )
    m_working.llmApi.activeProviderIndex = m_working.llmApi.providers.size() - 1;
  else if( m_working.llmApi.activeProviderIndex > pi )
    m_working.llmApi.activeProviderIndex -= 1;

  rebuildProviders();
  updateSummary();
  updateValidation();
  refreshPreview();
}//removeProvider()


void LlmConfigWindow::setActiveProvider( const size_t pi )
{
  if( pi < m_working.llmApi.providers.size() )
    m_working.llmApi.activeProviderIndex = pi;
  rebuildProviders();   // re-checks exactly one box, restyles cards
  updateSummary();
  updateValidation();
  refreshPreview();
}//setActiveProvider()


void LlmConfigWindow::addModel( const size_t pi )
{
  if( pi >= m_working.llmApi.providers.size() )
    return;
  m_working.llmApi.providers[pi].models.push_back( ModelInfo() );
  rebuildProviders();
  updateValidation();
  refreshPreview();
}//addModel()


void LlmConfigWindow::removeModel( const size_t pi, const size_t mi )
{
  if( pi >= m_working.llmApi.providers.size() )
    return;
  ApiProvider &p = m_working.llmApi.providers[pi];
  if( mi >= p.models.size() )
    return;
  p.models.erase( p.models.begin() + mi );
  if( p.models.empty() )
    p.activeModelIndex = 0;
  else if( p.activeModelIndex >= p.models.size() )
    p.activeModelIndex = p.models.size() - 1;
  else if( p.activeModelIndex > mi )
    p.activeModelIndex -= 1;

  rebuildProviders();
  updateSummary();
  updateValidation();
  refreshPreview();
}//removeModel()


void LlmConfigWindow::setActiveModel( const size_t pi, const size_t mi )
{
  if( pi >= m_working.llmApi.providers.size() )
    return;
  ApiProvider &p = m_working.llmApi.providers[pi];
  if( mi < p.models.size() )
    p.activeModelIndex = mi;
  rebuildProviders();
  updateSummary();
  updateValidation();
  refreshPreview();
}//setActiveModel()


void LlmConfigWindow::updateSummary()
{
  if( !m_activeSummary )
    return;

  WString summary = WString::tr("lcw-active-none");
  if( m_working.llmApi.activeProviderIndex < m_working.llmApi.providers.size() )
  {
    const ApiProvider &p = m_working.llmApi.providers[m_working.llmApi.activeProviderIndex];
    const string host = host_from_url( p.apiEndpoint );
    string model;
    if( p.activeModelIndex < p.models.size() )
      model = p.models[p.activeModelIndex].name;
    if( !host.empty() || !model.empty() )
      summary = WString::tr("lcw-active-fmt")
                  .arg( host.empty() ? WString::tr("lcw-new-provider") : WString::fromUTF8(host) )
                  .arg( model.empty() ? WString::tr("lcw-active-none") : WString::fromUTF8(model) );
  }

  m_activeSummary->setText( WString::tr("lcw-active-prefix").arg( summary ) );
}//updateSummary()


void LlmConfigWindow::updateValidation()
{
  if( !m_validationSummary )
    return;

  // Hard errors make the file un-loadable; disable Accept.
  bool blocking = m_working.llmApi.providers.empty();
  for( const ApiProvider &p : m_working.llmApi.providers )
    blocking = blocking || p.models.empty();

  // Soft warnings (non-blocking; mirror the loader's warnings).
  vector<WString> warnings;
  if( m_working.llmApi.activeProviderIndex < m_working.llmApi.providers.size() )
  {
    const size_t pi = m_working.llmApi.activeProviderIndex;
    const ApiProvider &p = m_working.llmApi.providers[pi];
    if( p.apiEndpoint.empty() )
      warnings.push_back( WString::tr("lcw-warn-endpoint-empty") );
    if( isPlaceholderToken(p.bearerToken) )
      warnings.push_back( WString::tr("lcw-warn-token-placeholder") );
    if( p.activeModelIndex < p.models.size() )
    {
      const ModelInfo &am = p.models[p.activeModelIndex];
      if( am.name.empty() )
        warnings.push_back( WString::tr("lcw-warn-model-empty") );
      if( (resolvedFormat(pi) == ApiFormat::Anthropic) && (am.maxTokens <= 0) )
        warnings.push_back( WString::tr("lcw-warn-anthropic-max") );
    }
  }
#if( MCP_ENABLE_AUTH )
  if( m_working.mcpServer.enabled
     && (m_working.mcpServer.bearerToken == LlmConfig::McpServer::sm_invalid_bearer_token) )
  {
    warnings.push_back( WString::tr("lcw-warn-mcp-token") );
  }
  if( m_mcpTokenWarn )
    m_mcpTokenWarn->setHidden( m_working.mcpServer.bearerToken != LlmConfig::McpServer::sm_invalid_bearer_token );
#endif

  m_validationSummary->removeStyleClass( "LcwValidOk" );
  m_validationSummary->removeStyleClass( "LcwValidWarn" );
  if( warnings.empty() )
  {
    m_validationSummary->setText( WString::tr("lcw-config-valid") );
    m_validationSummary->addStyleClass( "LcwValidOk" );
    m_validationSummary->setToolTip( WString() );
  }else
  {
    m_validationSummary->setText( WString::tr("lcw-to-review").arg(static_cast<int>(warnings.size())) );
    m_validationSummary->addStyleClass( "LcwValidWarn" );
    WString tip;
    for( size_t i = 0; i < warnings.size(); ++i )
      tip += (i ? WString::fromUTF8("\n") : WString()) + warnings[i];
    m_validationSummary->setToolTip( tip );
  }

  if( m_acceptBtn )
  {
    const bool enable = !blocking && !m_writableDir.empty();
    m_acceptBtn->setEnabled( enable );
    if( blocking )
      m_acceptBtn->setToolTip( WString::tr("lcw-blocking-tip") );
    else if( m_writableDir.empty() )
      m_acceptBtn->setToolTip( WString::tr("lcw-save-unavailable") );
    else
      m_acceptBtn->setToolTip( WString() );
  }
}//updateValidation()


void LlmConfigWindow::refreshPreview()
{
  if( m_showXml && m_previewPre )
    m_previewPre->setText( WString::fromUTF8( LlmConfig::toXmlString( m_working ) ) );
}//refreshPreview()


void LlmConfigWindow::togglePreviewXml()
{
  m_showXml = !m_showXml;
  if( m_previewContainer )
    m_previewContainer->setHidden( !m_showXml );
  if( m_previewToggle )
    m_previewToggle->setText( m_showXml ? WString::tr("lcw-hide-xml") : WString::tr("lcw-preview-xml") );
  refreshPreview();
}//togglePreviewXml()


// ---------------- Fetch models ----------------

void LlmConfigWindow::installFetchJs()
{
  if( m_fetchJsInstalled )
    return;
  m_fetchJsInstalled = true;

  const string emit = m_fetchModelsResult->createCall( "payload" );

  // `deliver` is guarded so the timeout and the (aborted) fetch can't both report a result.
  const string js =
    "window.llmFetchModels = function(reqId, url, headersJson){"
      "var headers={};try{headers=JSON.parse(headersJson);}catch(e){}"
      "var done=false;"
      "var deliver=function(obj){if(done)return;done=true;var payload=JSON.stringify(obj);" + emit + ";};"
      "var ctrl=(typeof AbortController!=='undefined')?new AbortController():null;"
      "var to=setTimeout(function(){if(ctrl)ctrl.abort();"
        "deliver({reqId:reqId,ok:false,status:0,body:'',error:'request timed out'});},20000);"
      "try{"
        "fetch(url,{method:'GET',headers:headers,signal:ctrl?ctrl.signal:undefined})"
          ".then(function(r){var s=r.status;return r.text().then(function(b){clearTimeout(to);"
            "deliver({reqId:reqId,ok:(s>=200&&s<300),status:s,body:b,error:''});});})"
          ".catch(function(e){clearTimeout(to);deliver({reqId:reqId,ok:false,status:0,body:'',"
            "error:String((e&&e.message)?e.message:e)});});"
      "}catch(e){clearTimeout(to);deliver({reqId:reqId,ok:false,status:0,body:'',error:String(e)});}"
    "};";

  wApp->doJavaScript( js );
}//installFetchJs()


void LlmConfigWindow::requestFetchModels( const size_t pi )
{
  if( pi >= m_working.llmApi.providers.size() || pi >= m_uiState.size() )
    return;

  const ApiProvider &prov = m_working.llmApi.providers[pi];
  if( prov.apiEndpoint.empty() )
  {
    m_uiState[pi].fetchOk = false;
    m_uiState[pi].fetchStatusMsg = WString::tr("lcw-fetch-need-endpoint");
    rebuildProviders();
    return;
  }

  installFetchJs();

  const int reqId = m_nextFetchReqId++;
  m_uiState[pi].fetchReqId = reqId;
  m_uiState[pi].fetchInProgress = true;
  m_uiState[pi].fetchStatusMsg = WString();

  const string modelsUrl = deriveModelsUrl( prov.apiEndpoint );

  // Same auth/CORS headers the chat path uses for this wire format.
  nlohmann::json headersJson = nlohmann::json::object();
  try
  {
    std::unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( resolvedFormat(pi) );
    if( proto )
    {
      const vector<pair<string,string>> hdrs = proto->headers( modelsUrl, prov.bearerToken );
      for( const pair<string,string> &kv : hdrs )
        headersJson[kv.first] = kv.second;
    }
  }catch( std::exception & )
  {
    // proceed with no extra headers
  }

  rebuildProviders();  // reflect the "Fetching…" state

  const string call = "if(window.llmFetchModels){window.llmFetchModels("
    + std::to_string(reqId) + ","
    + WWebWidget::jsStringLiteral( modelsUrl ) + ","
    + WWebWidget::jsStringLiteral( headersJson.dump() ) + ");}";
  wApp->doJavaScript( call );
}//requestFetchModels()


void LlmConfigWindow::onFetchModelsResult( std::string payload )
{
  int reqId = 0;
  bool ok = false;
  int status = 0;
  string body, error;
  try
  {
    const nlohmann::json j = nlohmann::json::parse( payload );
    reqId  = j.value( "reqId", 0 );
    ok     = j.value( "ok", false );
    status = j.value( "status", 0 );
    body   = j.value( "body", string() );
    error  = j.value( "error", string() );
  }catch( std::exception & )
  {
    return;
  }

  // Find the provider awaiting this request (it may have been removed/reordered).
  size_t pi = m_uiState.size();
  for( size_t i = 0; i < m_uiState.size(); ++i )
  {
    if( m_uiState[i].fetchInProgress && (m_uiState[i].fetchReqId == reqId) )
    {
      pi = i;
      break;
    }
  }
  if( pi >= m_uiState.size() )
    return;

  m_uiState[pi].fetchInProgress = false;

  if( !ok )
  {
    m_uiState[pi].fetchOk = false;
    m_uiState[pi].fetchStatusMsg = (status >= 400)
        ? WString::tr("lcw-fetch-http-error").arg(status)
        : WString::tr("lcw-fetch-network-error");
    rebuildProviders();
    return;
  }

  // Collect model ids from a few common shapes: {data:[{id|name}]}, [..], {models:[{id|name}]}.
  vector<string> ids;
  auto add_id = [&ids]( const nlohmann::json &el ){
    string id;
    if( el.is_string() )
      id = el.get<string>();
    else if( el.is_object() )
    {
      if( el.contains("id") && el["id"].is_string() )       id = el["id"].get<string>();
      else if( el.contains("name") && el["name"].is_string() ) id = el["name"].get<string>();
    }
    SpecUtils::trim( id );
    if( !id.empty() && (std::find(ids.begin(), ids.end(), id) == ids.end()) )
      ids.push_back( id );
  };

  try
  {
    const nlohmann::json j = nlohmann::json::parse( body );
    if( j.is_object() && j.contains("data") && j["data"].is_array() )
    {
      for( const nlohmann::json &el : j["data"] ) add_id( el );
    }else if( j.is_array() )
    {
      for( const nlohmann::json &el : j ) add_id( el );
    }else if( j.is_object() && j.contains("models") && j["models"].is_array() )
    {
      for( const nlohmann::json &el : j["models"] ) add_id( el );
    }
  }catch( std::exception & )
  {
    m_uiState[pi].fetchOk = false;
    m_uiState[pi].fetchStatusMsg = WString::tr("lcw-fetch-parse-error");
    rebuildProviders();
    return;
  }

  std::sort( ids.begin(), ids.end() );
  m_uiState[pi].discoveredModelIds = ids;
  m_uiState[pi].fetchOk = !ids.empty();
  m_uiState[pi].fetchStatusMsg = ids.empty() ? WString::tr("lcw-fetch-none")
                                             : WString::tr("lcw-fetch-found").arg(static_cast<int>(ids.size()));
  rebuildProviders();
}//onFetchModelsResult()


std::string LlmConfigWindow::deriveModelsUrl( const std::string &chatEndpoint )
{
  string u = chatEndpoint;
  SpecUtils::trim( u );

  const size_t q = u.find_first_of( "?#" );
  if( q != string::npos )
    u = u.substr( 0, q );
  while( !u.empty() && (u.back() == '/') )
    u.pop_back();
  if( u.empty() )
    return u;

  // A chat/responses/messages endpoint has its models list as a sibling: replace the trailing verb.
  const char * const suffixes[] = { "/chat/completions", "/responses", "/messages" };
  for( const char * const suf : suffixes )
  {
    if( SpecUtils::iends_with( u, suf ) )
      return u.substr( 0, u.size() - strlen(suf) ) + "/models";
  }

  // No recognized suffix: treat the endpoint as a base and append "/models"
  // (e.g. ".../v1" -> ".../v1/models").
  return u + "/models";
}//deriveModelsUrl()


// ---------------- Save ----------------

bool LlmConfigWindow::saveConfig()
{
  if( m_writableDir.empty() )
    return false;
  const string path = SpecUtils::append_path( m_writableDir, "llm_config.xml" );
  return LlmConfig::saveToFile( m_working, path );
}//saveConfig()


void LlmConfigWindow::handleAccept()
{
  if( m_writableDir.empty() )
    return;  // Accept is disabled when there is nowhere to save.

  // Import mode: the user opened a config file to install it.  Warn before clobbering any existing
  //  config; there is nothing to preserve on the session side, so just save+apply on confirmation.
  if( !m_importConfigPath.empty() )
  {
    const string existing = SpecUtils::append_path( m_writableDir, "llm_config.xml" );
    if( SpecUtils::is_file(existing) )
    {
      SimpleDialog *dialog = new SimpleDialog( WString::tr("lcw-overwrite-title") );
      WText *txt = new WText( WString::tr("lcw-overwrite-msg"), dialog->contents() );
      txt->setTextFormat( Wt::XHTMLText );

      WPushButton *ok = dialog->addButton( WString::tr("lcw-overwrite-confirm") );
      dialog->addButton( WString::tr("lcw-cancel") );  // just dismisses; settings window stays open
      ok->clicked().connect( std::bind( [this](){ doAccept( true ); } ) );
      return;
    }//if( an existing config would be overwritten )

    doAccept( true );
    return;
  }//if( import mode )

  const bool hasConversation = m_hasConversation && m_hasConversation();
  if( !hasConversation )
  {
    // Nothing to protect - save and apply right away.
    doAccept( true );
    return;
  }

  // There is a conversation; figure out whether applying these settings would disrupt it by
  // comparing against the currently-active config.
  LlmInterface::ConfigChange change = LlmInterface::ConfigChange::FormatChanged;
  try
  {
    const std::shared_ptr<const LlmConfig> current = InterSpecServer::llm_config();
    if( current )
      change = LlmInterface::classifyConfigChange( *current, m_working );
  }catch( const std::exception & )
  {
    change = LlmInterface::ConfigChange::FormatChanged;
  }

  if( change == LlmInterface::ConfigChange::SameModel )
  {
    // Fully compatible - the conversation is preserved untouched, so no warning needed.
    doAccept( true );
    return;
  }

  // Warn the user, offering: apply now, save the file without touching the running session, or cancel.
  SimpleDialog *dialog = new SimpleDialog( WString::tr("lcw-apply-title") );
  const WString msg = (change == LlmInterface::ConfigChange::FormatChanged)
                        ? WString::tr("lcw-apply-reset-msg")
                        : WString::tr("lcw-apply-keep-msg");
  WText *txt = new WText( msg, dialog->contents() );
  txt->setTextFormat( Wt::XHTMLText );

  WPushButton *cont = dialog->addButton( WString::tr("lcw-apply-continue") );
  WPushButton *fileOnly = dialog->addButton( WString::tr("lcw-apply-file-only") );
  dialog->addButton( WString::tr("lcw-cancel") );  // just dismisses the dialog; settings window stays open

  cont->clicked().connect( std::bind( [this](){ doAccept( true ); } ) );
  fileOnly->clicked().connect( std::bind( [this](){ doAccept( false ); } ) );
}//handleAccept()


void LlmConfigWindow::doAccept( const bool applyToSession )
{
  if( !saveConfig() )
  {
    SimpleDialog *dialog = new SimpleDialog( WString::tr("lcw-save-failed-title") );
    new WText( WString::tr("lcw-save-failed"), dialog->contents() );
    dialog->addButton( WString::tr("lcw-ok") );
    return;  // keep the settings window open so the user can retry
  }

  if( applyToSession && m_onApply )
    m_onApply();

  // hide() emits finished() -> deleteAuxWindow(this); must be the last use of *this.
  // (Use hide() rather than accept() so we don't go through WDialog::done(), which would
  //  emit finished() twice for an AuxWindow.)
  hide();
}//doAccept()


bool LlmConfigWindow::isPlaceholderToken( const std::string &token )
{
  string t = token;
  SpecUtils::trim( t );
  if( t.empty() )
    return true;
  SpecUtils::to_lower_ascii( t );
  return (t == "...") || (t == "sk-...") || (t == "sk-ant-...") || (t == "invalid-bearer-token")
         || (t.find("put-your") != string::npos) || (t.find("put-") != string::npos);
}//isPlaceholderToken()

#endif // USE_LLM_INTERFACE
