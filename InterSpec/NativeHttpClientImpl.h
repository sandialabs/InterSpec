#ifndef NATIVE_HTTP_CLIENT_IMPL_H
#define NATIVE_HTTP_CLIENT_IMPL_H
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

#include <mutex>
#include <atomic>
#include <memory>
#include <string>

#include "InterSpec/NativeHttpClient.h"

static_assert( USE_NATIVE_HTTP_CLIENT,
              "You should not include this header unless USE_NATIVE_HTTP_CLIENT is enabled" );

/** Shared state between `NativeHttp::Call` and whichever backend is compiled in.

 This header is deliberately NOT part of the public interface - only `NativeHttpClient.cpp` and
 the one backend translation unit include it.  Callers see `InterSpec/NativeHttpClient.h` only.

 The division of labour is: this file owns everything platform-independent (the handler, the
 one-shot completion latch, marshalling onto the Wt session thread), and the backend owns only
 the actual bytes-on-the-wire, calling `deliver*()` as things happen and checking
 `shouldStop()` to notice cancellation.
 */
namespace NativeHttp
{
namespace Detail
{
  /** The state a request needs, kept alive by `shared_ptr` so that a callback already in flight
   stays valid after the owning `Call` has gone.

   Derives from `enable_shared_from_this` because the marshalling helpers need to hand a
   reference to themselves to a lambda that will run later, on another thread. */
  struct RequestState : public std::enable_shared_from_this<RequestState>
  {
    Request request;
    StreamHandler handler;

    /** Wt session to marshal callbacks onto; empty means "run them on the backend thread". */
    std::string sessionId;

    /** Set by `~Call()`.  Every callback checks this *at the point it would run*, not before
     being posted, which is what makes "no callbacks after the destructor returns" true even
     for a callback already queued in the session's event queue. */
    std::atomic<bool> detached{ false };

    /** Set by `Call::cancel()`, or by `onChunk` returning false.  The backend polls this. */
    std::atomic<bool> stopRequested{ false };

    /** Why we are stopping, valid once `stopRequested` is set. */
    std::atomic<bool> stoppedByHandler{ false };

    /** Latch making `onComplete` fire exactly once, however many paths race to report it. */
    std::atomic<bool> completed{ false };

    /** Guards the `detached` re-check on the inline (empty-sessionId) path.

     Deliberately NOT held across the handler call: a handler that cancelled its own request would
     deadlock on it, since that routes back through deliverComplete() -> dispatch().  Handler
     invocations are therefore serialized only by the backend delivering them in order - which
     both real backends do (Wt emits on its strand, and the Apple delegate queue is forced
     serial).  A backend that delivered concurrently would need to serialize itself. */
    std::mutex handlerMutex;

    /** Bytes accumulated so far, for enforcing `Request::maxResponseSize` ourselves.  Atomic
     because a backend may deliver chunks from a thread pool rather than a single thread. */
    std::atomic<std::size_t> bytesReceived{ 0 };

    /** HTTP status from the response line, echoed into `Completion::status`; 0 until seen. */
    std::atomic<int> httpStatus{ 0 };

    /** Whatever the backend needs to find this transfer again in order to abort it (an
     NSURLSession task, a WinHTTP handle, a `Wt::Http::Client`).  Set by
     `startBackendRequest`, read by `cancelBackendRequest`, and guarded by `backendMutex`
     because those two can be called from different threads. */
    std::shared_ptr<void> backendData;
    std::mutex backendMutex;

    /** True once the backend should give up: cancelled, aborted by the handler, or finished. */
    bool shouldStop() const;

    /** Run `fcn` on the session thread (or inline when there is no session), dropping it if the
     `Call` has been detached or the session has died. */
    void dispatch( std::function<void()> fcn );

    void deliverHeaders( int status, HeaderList headers );

    /** Accounts the chunk against `maxResponseSize` and delivers it.  Returns false when the
     backend should stop - because the handler asked to, the cap was exceeded, or we were
     cancelled.  Note the handler's own false is observed one chunk late, since the handler runs
     asynchronously on the session thread; this is documented on `StreamHandler::onChunk`. */
    bool deliverChunk( std::string_view chunk );

    /** Deliver `onComplete`, at most once.  Later calls are ignored, so a backend may report a
     failure without first checking whether something else already did. */
    void deliverComplete( Error e, std::string detail );
  };//struct RequestState


  // The two entry points every backend provides.  Exactly one translation unit defines them,
  //  chosen by platform in CMakeLists.txt.  Until a given platform has a real backend,
  //  NativeHttpClient.cpp defines a stub that immediately reports Error::BackendUnavailable.

  /** Begin the transfer.  Must not block; must eventually cause exactly one
   `state->deliverComplete()`.  Takes a `shared_ptr` because the backend has to keep the state
   alive for as long as the transfer runs, independently of the owning `Call`. */
  void startBackendRequest( const std::shared_ptr<RequestState> &state );

  /** Ask the backend to abandon a transfer started by `startBackendRequest`.  Called with
   `state->stopRequested` already set, so a backend that only polls may leave this empty. */
  void cancelBackendRequest( const std::shared_ptr<RequestState> &state );
}//namespace Detail
}//namespace NativeHttp

#endif //NATIVE_HTTP_CLIENT_IMPL_H
