#ifndef LlmInteractionDisplay_h
#define LlmInteractionDisplay_h
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
#include <memory>
#include <string>
#include <functional>

#include <Wt/WPanel>
#include <Wt/WContainerWidget>

// Forward declarations
namespace Wt
{
  class WText;
  class WTextArea;
  class WPushButton;
}

class LlmInteraction;
class LlmInteractionTurn;
class LlmToolRequest;
class LlmToolResults;
class LlmInteractionError;
class LlmInteractionFinalResponse;


/** Base class for displaying a single LLM interaction turn (response, tool call, result, or error).
 *
 * Each turn is displayed as a WPanel with:
 * - A title bar with summary and menu icon
 * - A body with detailed content (created lazily on first expand)
 */
class LlmInteractionTurnDisplay : public Wt::WPanel
{
public:
  LlmInteractionTurnDisplay( std::shared_ptr<LlmInteractionTurn> turn,
                            Wt::WContainerWidget *parent = nullptr );

  virtual ~LlmInteractionTurnDisplay();

  /** Get the underlying interaction turn data */
  std::shared_ptr<LlmInteractionTurn> turn() const;

protected:
  /** Create the body content - called lazily on first expansion.
   *  Derived classes must implement this to create their specific content.
   */
  virtual void createBodyContent() = 0;

  /** Get the title text for this turn - derived classes should override */
  virtual Wt::WString getTitleText() const = 0;

  /** Create and show a dialog displaying JSON content */
  void showJsonDialog( const Wt::WString &title,
                      const std::string &jsonStr,
                      bool allowDownload = true );

  /** Create the menu icon and popup menu in the title bar.
   *  Derived classes can override to add custom menu items.
   */
  virtual void createMenuIcon();

  /** Handle panel expansion - creates body content if not yet created */
  void handleExpansion( bool expanded );

protected:
  std::shared_ptr<LlmInteractionTurn> m_turn;
  Wt::WPushButton *m_menuIcon;
  Wt::WContainerWidget *m_bodyContainer;
  bool m_bodyCreated;
};//class LlmInteractionTurnDisplay


/** Display for a final LLM response (assistant message with no tool calls) */
class LlmInteractionFinalResponseDisplay : public LlmInteractionTurnDisplay
{
public:
  LlmInteractionFinalResponseDisplay( std::shared_ptr<LlmInteractionFinalResponse> response,
                                      Wt::WContainerWidget *parent = nullptr );

  virtual ~LlmInteractionFinalResponseDisplay();

protected:
  virtual void createBodyContent() override;
  virtual Wt::WString getTitleText() const override;

private:
  std::shared_ptr<LlmInteractionFinalResponse> m_response;
};//class LlmInteractionFinalResponseDisplay


/** Display for LLM tool call requests */
class LlmToolRequestDisplay : public LlmInteractionTurnDisplay
{
public:
  LlmToolRequestDisplay( std::shared_ptr<LlmToolRequest> request,
                        Wt::WContainerWidget *parent = nullptr );

  virtual ~LlmToolRequestDisplay();

protected:
  virtual void createBodyContent() override;
  virtual Wt::WString getTitleText() const override;

private:
  std::shared_ptr<LlmToolRequest> m_request;
};//class LlmToolRequestDisplay


/** Display for tool call results */
class LlmToolResultsDisplay : public LlmInteractionTurnDisplay
{
public:
  LlmToolResultsDisplay( std::shared_ptr<LlmToolResults> results,
                        int nestingLevel = 0,
                        Wt::WContainerWidget *parent = nullptr );

  virtual ~LlmToolResultsDisplay();

protected:
  virtual void createBodyContent() override;
  virtual Wt::WString getTitleText() const override;
  virtual void createMenuIcon() override;

  /** Handle sub-agent conversation finishing */
  void handleSubAgentFinished( std::shared_ptr<LlmInteractionTurn> turn,
                               std::shared_ptr<LlmInteraction> subAgent );

private:
  std::shared_ptr<LlmToolResults> m_results;
  std::vector<class LlmInteractionDisplay*> m_subAgentDisplays;
  int m_nestingLevel;  // For passing to nested sub-agent displays

  // Map from sub-agent conversation to its result container div
  std::map<std::shared_ptr<LlmInteraction>, Wt::WContainerWidget*> m_subAgentResultDivs;
};//class LlmToolResultsDisplay


/** Display for interaction errors (timeouts, API errors, etc.) */
class LlmInteractionErrorDisplay : public LlmInteractionTurnDisplay
{
public:
  LlmInteractionErrorDisplay( std::shared_ptr<LlmInteractionError> error,
                             Wt::WContainerWidget *parent = nullptr );

  virtual ~LlmInteractionErrorDisplay();

  /** Set the retry callback - called when user clicks retry button */
  void setRetryCallback( std::function<void()> callback );

  /** Set the continue callback - called when user clicks "Continue Anyway" button */
  void setContinueCallback( std::function<void()> callback );

protected:
  virtual void createBodyContent() override;
  virtual Wt::WString getTitleText() const override;

  void handleRetry();
  void handleContinueAnyway();

private:
  std::shared_ptr<LlmInteractionError> m_error;
  Wt::WPushButton *m_retryBtn;
  Wt::WPushButton *m_continueBtn;
};//class LlmInteractionErrorDisplay


/** Display for a complete LLM interaction (user question + all responses/turns).
 *
 * This is the top-level widget for each conversation interaction, containing:
 * - User's question in the title
 * - Status indicator (pending, performing tools, complete)
 * - All interaction turns (responses, tool calls, results, errors) in the body
 */
class LlmInteractionDisplay : public Wt::WPanel
{
public:
  LlmInteractionDisplay( std::shared_ptr<LlmInteraction> interaction,
                        int nestingLevel = 0,
                        Wt::WContainerWidget *parent = nullptr );

  virtual ~LlmInteractionDisplay();

  /** Get the underlying interaction data */
  std::shared_ptr<LlmInteraction> interaction() const;

protected:
  /** Handle a new response being added to the interaction */
  void handleResponseAdded( std::shared_ptr<LlmInteractionTurn> turn );

  /** Handle a sub-agent conversation finishing */
  void handleSubAgentFinished( std::shared_ptr<LlmInteractionTurn> turn,
                              std::shared_ptr<LlmInteraction> subAgent );

  /** Handle the conversation finishing */
  void handleConversationFinished();

  /** Update the status indicator text and timer */
  void updateStatus();

  /** Update the summary section with start time, duration, and tokens */
  void updateSummary();

  /** Start JavaScript timer for status display */
  void startStatusTimer();

  /** Stop JavaScript timer and show final time */
  void stopStatusTimer();

  /** Create the title text with status */
  Wt::WString getTitleWithStatus() const;

  /** Create the menu icon and popup menu in the title bar */
  void createMenuIcon();

  /** Show a dialog displaying JSON content */
  void showJsonDialog( const Wt::WString &title,
                      const std::string &jsonStr,
                      bool allowDownload = true );

private:
  std::shared_ptr<LlmInteraction> m_interaction;
  Wt::WContainerWidget *m_turnContainer;
  Wt::WContainerWidget *m_summaryDiv;
  Wt::WText *m_statusText;
  Wt::WText *m_timerText;  // Separate widget for timer display
  Wt::WPushButton *m_menuIcon;

  int m_nestingLevel;  // For indentation of sub-agents
  bool m_isFinished;
  std::string m_timerId;  // Unique ID for JavaScript timer

  static int s_nextTimerId;
};//class LlmInteractionDisplay


#endif // LlmInteractionDisplay_h
