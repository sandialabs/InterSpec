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
#include <utility>
#include <algorithm>

#include <Wt/WServer>
#include <Wt/WApplication>

#include "InterSpec/NativeHttpClient.h"
#include "InterSpec/NativeHttpClientImpl.h"

using namespace std;

static_assert( USE_NATIVE_HTTP_CLIENT,
              "This file should only be compiled when USE_NATIVE_HTTP_CLIENT is enabled" );

namespace
{
/** The `std::error_category` for `NativeHttp::Error`.

 Messages here are the generic part of a failure; the specific part (the platform's own text)
 travels in `Completion::detail`.  These are worded for a user reading an error in the LLM
 panel, not for a developer reading a log.
 */
class NativeHttpCategory : public std::error_category
{
public:
  const char *name() const noexcept override
  {
    return "NativeHttp";
  }

  std::string message( int condition ) const override
  {
    switch( static_cast<NativeHttp::Error>(condition) )
    {
      case NativeHttp::Error::Ok:
        return "Success";
      case NativeHttp::Error::Cancelled:
        return "The request was cancelled";
      case NativeHttp::Error::AbortedByHandler:
        return "The request was stopped while the response was being read";
      case NativeHttp::Error::HostNotFound:
        return "The server name could not be resolved - check the endpoint URL, or whether DNS"
               " requires a VPN connection";
      case NativeHttp::Error::ConnectFailed:
        return "Could not connect to the server";
      case NativeHttp::Error::TlsHandshakeFailed:
        return "The secure connection could not be negotiated";
      case NativeHttp::Error::TlsCertificateUntrusted:
        return "The server's certificate was not trusted.  On a network that inspects TLS"
               " traffic, the inspecting proxy's root certificate needs to be installed in this"
               " machine's trust store, or given as a custom CA bundle in the LLM settings";
      case NativeHttp::Error::ProxyUnreachable:
        return "The configured proxy server could not be reached";
      case NativeHttp::Error::ProxyAuthRequired:
        return "The proxy server requires authentication that could not be provided"
               " automatically";
      case NativeHttp::Error::Timeout:
        return "The request took too long and was abandoned";
      case NativeHttp::Error::IdleTimeout:
        return "The server stopped sending data before the response was complete";
      case NativeHttp::Error::ResponseTooLarge:
        return "The response was larger than the configured limit";
      case NativeHttp::Error::OptionUnsupported:
        return "A requested network option is not supported by this platform's HTTP client";
      case NativeHttp::Error::BackendUnavailable:
        return "No native HTTP client is available in this build";
      case NativeHttp::Error::Unknown:
        break;
    }//switch( condition )

    return "An unrecognized network error occurred";
  }//message(...)
};//class NativeHttpCategory


/** Push queued changes to the browser.

 Callbacks arrive outside Wt's normal request/response cycle, so without this a handler that
 updates widgets would have no visible effect until the next user interaction.  Only worth doing
 for terminal events - doing it per chunk would flush once per network packet.
 */
void trigger_client_update()
{
  Wt::WApplication * const app = Wt::WApplication::instance();
  if( app )
    app->triggerUpdate();
}//trigger_client_update()

}//namespace


namespace NativeHttp
{

const std::error_category &error_category()
{
  static const NativeHttpCategory s_category;
  return s_category;
}//error_category()


std::error_code make_error_code( Error e )
{
  return std::error_code( static_cast<int>(e), error_category() );
}//make_error_code(...)


bool Detail::RequestState::shouldStop() const
{
  return detached.load() || stopRequested.load() || completed.load();
}//Detail::RequestState::shouldStop()


void Detail::RequestState::dispatch( std::function<void()> fcn )
{
  if( detached.load() )
    return;

  if( sessionId.empty() )
  {
    // No Wt session: run inline on whatever thread the backend called us from.  Used by unit
    //  tests and any non-GUI caller.
    //
    // The lock is released before invoking the handler.  Holding it across the call would
    //  deadlock the moment a handler called Call::cancel() on its own request, since that routes
    //  back into deliverComplete -> dispatch and the mutex is not recursive.  Serialization of
    //  backend threads is all that is needed here; the handler itself is the caller's business.
    {
      std::lock_guard<std::mutex> lock( handlerMutex );
      if( detached.load() )
        return;
    }
    fcn();
    return;
  }//if( sessionId.empty() )

  Wt::WServer * const server = Wt::WServer::instance();
  if( !server )
    return;

  // Capture the state by shared_ptr so it outlives the owning Call, and re-check `detached`
  //  inside the posted function rather than here - that is what makes the "no callbacks after
  //  ~Call() returns" guarantee hold for a callback that was already queued.
  //
  // Wt::WServer::post() with no delay goes through WIOService::schedule(0,...), which posts to
  //  an asio strand ("guarantees execution order"), and from there onto a per-session FIFO
  //  event queue - so these arrive in the order we post them.
  const std::shared_ptr<Detail::RequestState> self = shared_from_this();
  server->post( sessionId, [self, fcn](){
    if( self->detached.load() )
      return;
    std::lock_guard<std::mutex> lock( self->handlerMutex );
    if( !self->detached.load() )
      fcn();
  } );
}//Detail::RequestState::dispatch(...)


void Detail::RequestState::deliverHeaders( int status, HeaderList headers )
{
  httpStatus.store( status );

  if( !handler.onHeaders || shouldStop() )
    return;

  const std::function<void(int, const HeaderList &)> fcn = handler.onHeaders;
  dispatch( [fcn, status, headers](){ fcn( status, headers ); } );
}//Detail::RequestState::deliverHeaders(...)


bool Detail::RequestState::deliverChunk( std::string_view chunk )
{
  if( shouldStop() )
    return false;

  const std::size_t total = (bytesReceived += chunk.size());
  if( (request.maxResponseSize > 0) && (total > request.maxResponseSize) )
  {
    deliverComplete( Error::ResponseTooLarge,
                     "Response exceeded " + std::to_string(request.maxResponseSize) + " bytes" );
    return false;
  }//if( over the size cap )

  if( handler.onChunk )
  {
    // The view is only valid for this call, so the copy is not optional - the handler may not
    //  run until later, on the session thread.
    const std::function<bool(std::string_view)> fcn = handler.onChunk;
    const std::shared_ptr<Detail::RequestState> self = shared_from_this();
    std::string owned( chunk );
    dispatch( [self, fcn, owned = std::move(owned)](){
      if( !fcn( owned ) )
      {
        self->stoppedByHandler.store( true );
        self->stopRequested.store( true );
      }
    } );
  }//if( handler.onChunk )

  return !shouldStop();
}//Detail::RequestState::deliverChunk(...)


void Detail::RequestState::deliverComplete( Error e, std::string detail )
{
  // Exactly-once latch: several paths (backend error, our own timeout, cancellation) can race
  //  to finish a request, and all of them are allowed to just call this.
  if( completed.exchange( true ) )
    return;

  if( !handler.onComplete )
    return;

  Completion result;
  result.ec = make_error_code( e );
  result.status = httpStatus.load();
  result.detail = std::move( detail );
  result.backend = NativeHttp::backend();

  const std::function<void(const Completion &)> fcn = handler.onComplete;
  dispatch( [fcn, result](){
    fcn( result );
    trigger_client_update();
  } );
}//Detail::RequestState::deliverComplete(...)


Call::Call( std::shared_ptr<Detail::RequestState> state )
  : m_state( std::move(state) )
{
}//Call::Call(...)


Call::~Call()
{
  if( !m_state )
    return;

  // Silent teardown: mark detached first so any callback already queued on the session thread
  //  drops itself, then let the backend know it can stop.
  m_state->detached.store( true );
  m_state->stopRequested.store( true );
  Detail::cancelBackendRequest( m_state );
}//Call::~Call()


void Call::cancel()
{
  if( !m_state || m_state->completed.load() )
    return;

  m_state->stopRequested.store( true );
  Detail::cancelBackendRequest( m_state );
  m_state->deliverComplete( Error::Cancelled, "Cancelled by the application" );
}//Call::cancel()


std::unique_ptr<Call> start( const Request &req, StreamHandler handler,
                            const std::string &sessionId )
{
  const std::shared_ptr<Detail::RequestState> state = std::make_shared<Detail::RequestState>();
  state->request = req;
  state->handler = std::move( handler );
  state->sessionId = sessionId;

  std::unique_ptr<Call> call( new Call( state ) );

  // Everything that can go wrong before a single byte moves is still reported through
  //  onComplete, so the caller has one failure path rather than two.
  if( !available() )
  {
    state->deliverComplete( Error::BackendUnavailable,
                            std::string("No usable native HTTP backend: ") + backendName() );
    return call;
  }//if( !available() )

  if( req.url.empty() )
  {
    state->deliverComplete( Error::Unknown, "No endpoint URL was given" );
    return call;
  }//if( req.url.empty() )

  if( (req.url.compare(0, 7, "http://") != 0) && (req.url.compare(0, 8, "https://") != 0) )
  {
    state->deliverComplete( Error::Unknown, "Endpoint URL is not http(s): " + req.url );
    return call;
  }//if( not an http(s) URL )

  if( !req.proxyUrl.empty() && !supportsProxyOverride() )
  {
    state->deliverComplete( Error::OptionUnsupported,
                    std::string("An explicit proxy was configured, but the ") + backendName()
                    + " backend cannot honour one; clear it, or set the proxy system-wide." );
    return call;
  }//if( proxy asked for but unsupported )

  if( !req.caBundlePath.empty() && !supportsCaBundle() )
  {
    state->deliverComplete( Error::OptionUnsupported,
                    std::string("A custom CA bundle was configured, but the ") + backendName()
                    + " backend validates against the system certificate store; install the"
                      " certificate there instead." );
    return call;
  }//if( CA bundle asked for but unsupported )

  Detail::startBackendRequest( state );

  return call;
}//start(...)


// ---------------------------------------------------------------------------------------------
// Stub backend.
//
// Each platform's real backend replaces the definitions below from its own translation unit
// (NativeHttpClient_apple.mm, _winhttp.cpp, _wt.cpp), selected in CMakeLists.txt.  Until a given
// platform has one, this keeps USE_NATIVE_HTTP_CLIENT=ON building and linking, and makes the
// resulting behaviour explicit rather than mysterious.
// ---------------------------------------------------------------------------------------------
#if( !defined(NATIVE_HTTP_HAVE_BACKEND) )

bool available()
{
  return false;
}

Backend backend()
{
  return Backend::None;
}

const char *backendName()
{
  return "none (no backend compiled in)";
}

bool supportsProxyOverride()
{
  return false;
}

bool supportsCaBundle()
{
  return false;
}

namespace Detail
{
void startBackendRequest( const std::shared_ptr<Detail::RequestState> &state )
{
  state->deliverComplete( Error::BackendUnavailable,
                          "This build has no native HTTP backend for this platform" );
}

void cancelBackendRequest( const std::shared_ptr<Detail::RequestState> & )
{
}
}//namespace Detail

#endif //!defined(NATIVE_HTTP_HAVE_BACKEND)

}//namespace NativeHttp
