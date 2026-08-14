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

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <Wt/WText.h>
#include <Wt/WTextArea.h>
#include <Wt/WPushButton.h>
#include <Wt/WServer.h>
#include <Wt/WApplication.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/LlmSubAgentFollowup.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace Wt;
using namespace std;


LlmSubAgentFollowup::LlmSubAgentFollowup( std::shared_ptr<LlmInteraction> interaction )
  : WObject(),
    m_originalInteraction( interaction ),
    m_llmInterface( nullptr ),
    m_requestPending( false )
{
  assert( interaction );
  if( !interaction )
    throw std::invalid_argument( "LlmSubAgentFollowup: null interaction" );

  InterSpec * const interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    throw std::runtime_error( "LlmSubAgentFollowup: no InterSpec instance" );

  shared_ptr<const LlmConfig> config;
  try
  {
    config = InterSpecServer::llm_config();
  }catch( std::exception &e )
  {
    throw std::runtime_error( string("LlmSubAgentFollowup: failed to get LLM config: ") + e.what() );
  }

  if( !config || !config->llmApi.enabled )
    throw std::runtime_error( "LlmSubAgentFollowup: LLM API not enabled" );

  // Create a private LlmInterface, seed it with the original interaction's history,
  // and block tool calls so LLM responds with text only.
  // We use shallowClone() so new turns appended by LlmInterface go into the clone's
  // independent responses vector and do not affect the original interaction.
  m_llmInterface = std::make_shared<LlmInterface>( interspec, config );
  m_llmInterface->setBlockToolCalls( true );

  // Write debug output to a timestamped file so the request JSON can be inspected
  //{
  //  const auto now = std::chrono::system_clock::now();
  //  const std::time_t t = std::chrono::system_clock::to_time_t( now );
  //  std::ostringstream oss;
  //  oss << "follow_up_debug_" << std::put_time( std::localtime( &t ), "%Y%m%d_%H%M%S" ) << ".txt";
  //  m_llmInterface->setDebugFile( oss.str() );
  //}

  auto history = std::make_shared<LlmConversationHistory>();
  history->getConversations().push_back( m_originalInteraction->shallowClone() );
  m_llmInterface->setHistory( history );

  m_llmInterface->conversationFinished().connect( this, &LlmSubAgentFollowup::handleConversationFinished );
}//LlmSubAgentFollowup constructor


LlmSubAgentFollowup::~LlmSubAgentFollowup()
{
  cout << "Deleting LlmSubAgentFollowup" << endl;
}


void LlmSubAgentFollowup::startDeleteSelf()
{
  // Deferred, because this is called from inside a button's clicked() emission: Wt4 frees the
  //  object through its owner's unique_ptr, and tearing it down mid-emit would leave the emitting
  //  signal walking freed connections.  `WServer::post` runs the removal on the session thread in
  //  a later pass of the event loop.  See SimpleDialog::startDeleteSelf().
  WApplication * const app = WApplication::instance();
  if( !app )
    return;

  const std::string sessionId = app->sessionId();
  WServer::instance()->post( sessionId, [this, sessionId](){
    WApplication * const inner_app = WApplication::instance();
    if( inner_app && (inner_app->sessionId() == sessionId) )
    {
      inner_app->removeChild( this );   //owner is wApp; see the creation site
      inner_app->triggerUpdate();
    }
  } );
}//startDeleteSelf()


void LlmSubAgentFollowup::showDialog()
{
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Ask Followup Question" );

  WTextArea *questionArea = dialog->contents()->addNew<WTextArea>();
  questionArea->setRows( 4 );
  questionArea->setPlaceholderText( "Enter your question about this conversation..." );
  questionArea->addStyleClass( "LlmFollowupInput" );
  questionArea->setWidth( WLength( 100, WLength::Unit::Percentage ) );

  WPushButton *sendBtn = dialog->addButton( "Send" );
  WPushButton *cancelBtn = dialog->addButton( "Cancel" );

  cancelBtn->clicked().connect( std::bind( [this, dialog]() {
    assert( !m_requestPending );
    if( !m_requestPending )
      startDeleteSelf();
  } ) );

  sendBtn->clicked().connect( std::bind( [this, dialog, questionArea]() {
    const string question = questionArea->text().toUTF8();
    if( question.empty() )
      return;

    // Set pending BEFORE accept() so that finished() handler sees it and won't delete this
    m_requestPending = true;
    m_llmInterface->sendUserMessage( question );
  } ) );

  dialog->show();
}//showDialog()


void LlmSubAgentFollowup::handleConversationFinished()
{
  m_requestPending = false;

  // Extract the response text from the last turn of the most recent conversation
  string responseText = "(No response received)";

  const shared_ptr<LlmConversationHistory> history = m_llmInterface->getHistory();
  if( history )
  {
    const vector<shared_ptr<LlmInteraction>> &convos = history->getConversations();
    // The most recently added conversation (after original) is at the back
    if( convos.size() >= 2 )
    {
      const shared_ptr<LlmInteraction> &lastConvo = convos.back();
      if( lastConvo )
      {
        // Walk backwards through responses to find the last FinalLlmResponse
        for( auto it = lastConvo->responses.rbegin(); it != lastConvo->responses.rend(); ++it )
        {
          const LlmInteractionFinalResponse *resp
            = dynamic_cast<const LlmInteractionFinalResponse *>( it->get() );
          if( resp )
          {
            responseText = resp->content();
            break;
          }
        }//for reverse iterate responses
      }//if( lastConvo )
    }//if( convos.size() >= 2 )
  }//if( history )

  showResponseDialog( responseText );
}//handleConversationFinished()


void LlmSubAgentFollowup::showResponseDialog( const std::string &responseText )
{
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Followup Response" );

  WTextArea *responseArea = dialog->contents()->addNew<WTextArea>();
  responseArea->setReadOnly( true );
  responseArea->setText( responseText );
  responseArea->setRows( 15 );
  responseArea->addStyleClass( "LlmFollowupResponse" );
  responseArea->setWidth( WLength( 100, WLength::Unit::Percentage ) );

  WText *label = dialog->contents()->addNew<WText>( "Ask another question:" );
  label->setPadding( WLength( 8 ), Wt::Side::Top );

  WTextArea *questionArea = dialog->contents()->addNew<WTextArea>();
  questionArea->setRows( 3 );
  questionArea->setPlaceholderText( "Enter another question, or close..." );
  questionArea->addStyleClass( "LlmFollowupInput" );
  questionArea->setWidth( WLength( 100, WLength::Unit::Percentage ) );

  WPushButton *sendBtn = dialog->addButton( "Send Followup" );
  WPushButton *closeBtn = dialog->addButton( "Close" );

  closeBtn->clicked().connect( std::bind( [this, dialog](){
    assert( !m_requestPending );
    if( !m_requestPending )
      startDeleteSelf();
  } ) );

  sendBtn->clicked().connect( std::bind( [this, dialog, questionArea](){
    const string question = questionArea->text().toUTF8();
    if( question.empty() )
      return;

    // Set pending BEFORE accept() so finished() handler won't delete this
    m_requestPending = true;
    m_llmInterface->sendUserMessage( question );
  } ) );

  dialog->show();
}//showResponseDialog()
