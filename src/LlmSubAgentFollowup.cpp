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

#include <Wt/WText>
#include <Wt/WTextArea>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

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
  : WObject( WApplication::instance() ),
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
  m_llmInterface = std::make_unique<LlmInterface>( interspec, config );
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


void LlmSubAgentFollowup::showDialog()
{
  SimpleDialog *dialog = new SimpleDialog( "Ask Followup Question" );

  WTextArea *questionArea = new WTextArea( dialog->contents() );
  questionArea->setRows( 4 );
  questionArea->setPlaceholderText( "Enter your question about this conversation..." );
  questionArea->addStyleClass( "LlmFollowupInput" );
  questionArea->setWidth( WLength( 100, WLength::Percentage ) );

  WPushButton *sendBtn = dialog->addButton( "Send" );
  WPushButton *cancelBtn = dialog->addButton( "Cancel" );

  cancelBtn->clicked().connect( std::bind( [this, dialog]() {
    assert( !m_requestPending );
    if( !m_requestPending )
      delete this;
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
  SimpleDialog *dialog = new SimpleDialog( "Followup Response" );

  WTextArea *responseArea = new WTextArea( dialog->contents() );
  responseArea->setReadOnly( true );
  responseArea->setText( responseText );
  responseArea->setRows( 15 );
  responseArea->addStyleClass( "LlmFollowupResponse" );
  responseArea->setWidth( WLength( 100, WLength::Percentage ) );

  WText *label = new WText( "Ask another question:", dialog->contents() );
  label->setPadding( WLength( 8 ), Wt::Top );

  WTextArea *questionArea = new WTextArea( dialog->contents() );
  questionArea->setRows( 3 );
  questionArea->setPlaceholderText( "Enter another question, or close..." );
  questionArea->addStyleClass( "LlmFollowupInput" );
  questionArea->setWidth( WLength( 100, WLength::Percentage ) );

  WPushButton *sendBtn = dialog->addButton( "Send Followup" );
  WPushButton *closeBtn = dialog->addButton( "Close" );

  closeBtn->clicked().connect( std::bind( [this, dialog](){
    assert( !m_requestPending );
    if( !m_requestPending )
      delete this;
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
