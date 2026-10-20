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
#include <cassert>

#include <boost/asio/steady_timer.hpp>

#include <Wt/WServer.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>

#include "InterSpec/LlmSessionTimer.h"

using namespace std;


/** The timer's state; owned by the `LlmSessionTimer`, and reached from its pending waits through weak
 pointers, so a wait that outlives the timer does nothing.

 Only used in the session (with its lock held), except `session_id`, which is constant.
 */
struct LlmSessionTimer::State
{
  std::string session_id;
  std::function<void()> callback;
  std::chrono::milliseconds interval{ 0 };
  bool single_shot = false;
  bool active = false;

  /** Incremented by every start() and stop(), so only the most recent wait fires. */
  uint64_t generation = 0;

  std::unique_ptr<boost::asio::steady_timer> timer;
};//struct LlmSessionTimer::State


namespace
{
  void arm( const std::shared_ptr<LlmSessionTimer::State> &state );

  /** Runs in the session, when a wait of `generation` has expired. */
  void fire( const std::weak_ptr<LlmSessionTimer::State> &weak_state, const uint64_t generation )
  {
    // Kept alive while the callback runs, in case it destroys the LlmSessionTimer.
    const std::shared_ptr<LlmSessionTimer::State> state = weak_state.lock();
    if( !state || !state->active || (state->generation != generation) )
      return;

    if( state->single_shot )
      state->active = false;
    else
      arm( state );

    const std::function<void()> callback = state->callback;
    if( callback )
      callback();
  }//fire(...)


  /** Starts a wait of `state->interval`; at its end, `fire` is scheduled into the session. */
  void arm( const std::shared_ptr<LlmSessionTimer::State> &state )
  {
    Wt::WServer * const server = Wt::WServer::instance();
    if( !server )
      return;

    if( !state->timer )
      state->timer = std::make_unique<boost::asio::steady_timer>( server->ioService() );

    const uint64_t generation = ++state->generation;
    const std::weak_ptr<LlmSessionTimer::State> weak_state = state;
    const std::string session_id = state->session_id;

    // (Setting the expiry cancels an earlier wait, whose handler then sees operation_aborted.)
    state->timer->expires_after( state->interval );
    state->timer->async_wait( [weak_state, generation, session_id]( const boost::system::error_code &ec ){
      if( ec || weak_state.expired() )
        return;

      Wt::WServer * const server = Wt::WServer::instance();
      if( server )
        server->schedule( std::chrono::milliseconds(1), session_id, [weak_state, generation](){
          fire( weak_state, generation );
        } );
    } );
  }//arm(...)
}//namespace


LlmSessionTimer::LlmSessionTimer( std::function<void()> callback )
  : m_state( std::make_shared<State>() )
{
  Wt::WApplication * const app = Wt::WApplication::instance();
  assert( app );
  m_state->session_id = app ? app->sessionId() : std::string();
  m_state->callback = std::move( callback );
}


LlmSessionTimer::~LlmSessionTimer()
{
  stop();
}


void LlmSessionTimer::setInterval( const std::chrono::milliseconds interval )
{
  m_state->interval = interval;
}


void LlmSessionTimer::setSingleShot( const bool single_shot )
{
  m_state->single_shot = single_shot;
}


void LlmSessionTimer::start()
{
  m_state->active = true;
  arm( m_state );
}


void LlmSessionTimer::stop()
{
  m_state->active = false;
  ++m_state->generation;
  if( m_state->timer )
    m_state->timer->cancel();
}


bool LlmSessionTimer::isActive() const
{
  return m_state->active;
}
