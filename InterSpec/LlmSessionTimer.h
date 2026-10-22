#ifndef LlmSessionTimer_h
#define LlmSessionTimer_h
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
#include <memory>
#include <functional>


/** A timer whose callback runs in the Wt session that created it (on a server thread, with the
 session lock held and `wApp` set), like `Wt::WTimer` - but driven by the server instead of by
 JavaScript in the browser, so it also works with no browser attached (the stdio MCP server), or a
 browser that has gone to sleep.

 Must be created, used, and destroyed in its session (i.e., with the session lock held); it may be
 destroyed from inside its own callback.  A stopped (or destroyed) timer's callback is never
 called, even if it was already due.  Pending waits are cancelled when stopped, so a long interval
 does not hold up the Wt server stopping.
 */
class LlmSessionTimer
{
public:
  /** Creates a stopped timer for the current session (`wApp` must be set). */
  explicit LlmSessionTimer( std::function<void()> callback );
  ~LlmSessionTimer();

  LlmSessionTimer( const LlmSessionTimer & ) = delete;
  LlmSessionTimer &operator=( const LlmSessionTimer & ) = delete;

  void setInterval( std::chrono::milliseconds interval );
  void setSingleShot( bool single_shot );

  /** Starts the timer, or restarts it (the interval then counts from now). */
  void start();

  void stop();

  bool isActive() const;

  /** Shared with the pending waits; see the .cpp. */
  struct State;

private:
  std::shared_ptr<State> m_state;
};//class LlmSessionTimer

#endif //LlmSessionTimer_h
