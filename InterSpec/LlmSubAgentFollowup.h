#ifndef LlmSubAgentFollowup_h
#define LlmSubAgentFollowup_h
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

#include <memory>
#include <string>

#include <Wt/WObject>

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

class LlmInterface;
struct LlmInteraction;


/** Lightweight utility for asking followup questions about a sub-agent conversation.

 Uses a private LlmInterface instance seeded with the original conversation history,
 so tool calls are blocked and no changes are made to the original conversation.

 Creates a new dialog for each question/response cycle; self-deletes when the user
 cancels a dialog with no request pending.
 */
class LlmSubAgentFollowup : public Wt::WObject
{
public:
  explicit LlmSubAgentFollowup( std::shared_ptr<LlmInteraction> interaction );

  virtual ~LlmSubAgentFollowup();

  /** Show the initial dialog prompting for a followup question. */
  void showDialog();

private:
  /** Called when the LlmInterface signals a finished conversation. */
  void handleConversationFinished();

  /** Show a dialog with the LLM response text and another question input. */
  void showResponseDialog( const std::string &responseText );

  std::shared_ptr<LlmInteraction> m_originalInteraction; // kept for shallowClone on construction
  std::unique_ptr<LlmInterface> m_llmInterface;          // private interface; tool calls blocked
  bool m_requestPending;
};//class LlmSubAgentFollowup


#endif // LlmSubAgentFollowup_h
