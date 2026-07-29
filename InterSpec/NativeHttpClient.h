#ifndef NATIVE_HTTP_CLIENT_H
#define NATIVE_HTTP_CLIENT_H
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
#include <string>
#include <vector>
#include <cstddef>
#include <utility>
#include <functional>
#include <string_view>
#include <system_error>
#include <type_traits>

static_assert( USE_NATIVE_HTTP_CLIENT,
              "You should not include this header unless USE_NATIVE_HTTP_CLIENT is enabled" );

/** A minimal HTTPS transport backed by the operating system's own HTTP stack.

 InterSpec normally issues LLM provider requests from the browser via `fetch()` (see
 `LlmInterface::makeApiCallWithId`).  That is deliberately the primary path: the browser
 already has the system trust store, the system proxy configuration, PAC/WPAD evaluation and
 often single-sign-on, none of which we want to reimplement.  It only works, though, for
 providers that return permissive CORS headers.  For a provider that does not, the request has
 to be made from C++ instead - and it still has to survive the same corporate networks.

 Hence one narrow interface with a per-platform backend behind it, each delegating to the
 system stack rather than to a bundled TLS/HTTP library:

 | Platform     | Backend            | Corporate coverage                        |
 |--------------|--------------------|-------------------------------------------|
 | Windows      | WinHTTP            | full, including Negotiate/NTLM SSO        |
 | macOS + iOS  | NSURLSession       | full except SSO (best effort)             |
 | Linux        | `Wt::Http::Client` | system trust only, no proxy support       |
 | Android      | `Wt::Http::Client` | system trust only, no proxy support       |

 Two caveats on the Windows row.  It was written without a Windows machine to develop on: it is
 verified to compile and link for x64 and x86 (mingw-w64 ships the same WinHTTP headers and
 import library), and its request/response, cancellation, size-cap and error-mapping paths have
 been exercised under CrossOver - but Wine reimplements WinHTTP, so that says nothing about
 Schannel, proxy discovery, or SSO.  And the Negotiate/NTLM handshake in particular is
 implemented per Microsoft's documented pattern but has never been run against a real proxy,
 because nothing here can host one.  Treat the whole row as needing a pass on real hardware.

 Nothing above this header knows which backend is in use, and no platform `#ifdef` reaches
 above it.  Not supporting proxies on Linux/Android is a deliberate, accepted limitation.
 */
namespace NativeHttp
{
  /** HTTP header field names paired with their values, in wire order. */
  using HeaderList = std::vector<std::pair<std::string,std::string>>;


  /** Which platform backend was compiled in; see `backendName()` for a display string. */
  enum class Backend
  {
    None,          ///< No backend available for this platform.
    NSUrlSession,  ///< macOS and iOS.
    WinHttp,       ///< Windows.
    WtHttpClient   ///< Linux and Android.
  };//enum class Backend


  /** Failure modes this transport reports.

   These are deliberately coarse: each maps to something a user or a help-desk ticket can act
   on, rather than mirroring any one platform's error space.  The precise platform message is
   carried separately in `Completion::detail`.
   */
  enum class Error
  {
    Ok = 0,                   ///< Success; `Completion::ec` is falsy.
    Cancelled,                ///< `Call::cancel()` was invoked, or the `Call` was destroyed.
    AbortedByHandler,         ///< `StreamHandler::onChunk` returned false.
    HostNotFound,             ///< DNS resolution failed.
    ConnectFailed,            ///< TCP connection refused, reset, or unreachable.
    TlsHandshakeFailed,       ///< TLS negotiation failed for a reason other than trust.
    TlsCertificateUntrusted,  ///< Server certificate did not validate; the TLS-inspection case.
    ProxyUnreachable,         ///< A proxy was configured or discovered but could not be reached.
    ProxyAuthRequired,        ///< Proxy returned 407 and no usable credentials were available.
    Timeout,                  ///< `Request::totalTimeout` elapsed.
    IdleTimeout,              ///< `Request::idleTimeout` elapsed with no bytes received.
    ResponseTooLarge,         ///< Accumulated body exceeded `Request::maxResponseSize`.
    OptionUnsupported,        ///< A requested escape hatch is not supported by this backend.
    BackendUnavailable,       ///< No backend compiled in, or it cannot do TLS (see below).
    Unknown                   ///< Anything not classified above; `detail` carries the specifics.
  };//enum class Error

  /** The `std::error_category` that `Error` values belong to. */
  InterSpec_API const std::error_category &error_category();

  /** Convert an `Error` to a `std::error_code`; found by ADL, so `std::error_code ec = Error::Timeout;`
   works (see the `std::is_error_code_enum` specialization at the bottom of this header). */
  InterSpec_API std::error_code make_error_code( Error e );


  /** Everything needed to issue one request.

   Field defaults are the values wanted for an LLM completion; only `url` and `body` normally
   need setting.  Every escape hatch defaults to "use whatever the platform would do", which is
   the behaviour that works for the overwhelming majority of users.
   */
  struct Request
  {
    std::string url;
    std::string method = "POST";

    /** Sent as given.  A `Content-Type` is *not* added for you.  Note these routinely carry the
     provider's bearer token, so they must never be logged verbatim.

     Do not set `Accept-Encoding`: this transport does not decompress response bodies, and not
     every backend does it for you. */
    HeaderList headers;

    std::string body;

    /** Redirects are not followed unless this is set, because the backends disagree by default
     (NSURLSession and WinHTTP follow, `Wt::Http::Client` does not) and an LLM API endpoint that
     redirects is nearly always a misconfiguration worth surfacing rather than chasing. */
    bool followRedirects = false;

    /** Reset on every byte received.

     An LLM can legitimately think for a minute or more before emitting its first byte, so a
     short overall timeout is wrong here.  An idle timeout is the right shape: it tolerates a
     slow start but still notices a connection that a proxy has silently reaped (commonly done
     somewhere between 60 s and 300 s).
     */
    std::chrono::seconds idleTimeout{ 120 };

    /** Wall-clock ceiling on the whole request; zero means no ceiling.

     Must stay below `LlmInterface::sm_watchdog_timeout_ms` (600 s), so that a transport failure
     is reported as a transport failure rather than being recovered by the conversation
     watchdog.  The default matches the 300 s ceiling the browser path already imposes.
     */
    std::chrono::seconds totalTimeout{ 300 };

    /** Cap on the accumulated response body; zero means unlimited.

     Enforced by this transport rather than being handed to the backend, because the backends do
     not agree on what a cap means.  Note in particular that `Wt::Http::Client` must NOT be given
     `setMaximumResponseSize(0)`: that stops it accumulating the body, and its "response complete"
     test is `body().size() >= contentLength_`, so a Content-Length response would never be seen
     as finished and every successful request would hang until the idle timeout.
     */
    std::size_t maxResponseSize = 64 * 1024 * 1024;

    // The remaining fields are user-facing escape hatches, for the corporate networks that no
    // HTTP library gets right every time.  Empty/false means "use the platform default".  A
    // backend that cannot honour one reports Error::OptionUnsupported rather than silently
    // ignoring it - see supportsProxyOverride() and supportsCaBundle().

    /** Explicit proxy, e.g. "http://proxy.example.com:8080"; empty means use the system proxy
     configuration (including WPAD/PAC discovery, where the backend supports it). */
    std::string proxyUrl;

    /** A PEM bundle file, or a directory of hash-named PEM files, to trust *in addition to*
     the system trust store; empty means use the system trust store alone. */
    std::string caBundlePath;

    /** Accept any server certificate.  This disables the protection TLS exists to provide and
     must only ever be reachable behind an explicit, clearly-worded user opt-in. */
    bool disableCertVerification = false;
  };//struct Request


  /** How a request ended; passed to `StreamHandler::onComplete` exactly once. */
  struct Completion
  {
    /** Falsy on success.  Compares equal to the corresponding `Error` enumerator.

     Explicitly initialized rather than left default-constructed: a default `std::error_code` is
     `(0, system_category)`, which is falsy but does *not* compare equal to `Error::Ok`. */
    std::error_code ec = make_error_code( Error::Ok );

    /** The HTTP status, or 0 if the request failed before a response line was read.

     Duplicated from `onHeaders` on purpose.  A proxy 407 or a gateway 502 arrives with a
     perfectly successful `ec` and an HTML body, and without this every caller would have to
     keep its own status variable plus a "did I ever get headers" flag just to tell that apart
     from a real response.
     */
    int status = 0;

    /** The platform's own description, already mapped to something readable, suitable for
     showing to a user or pasting into a bug report.

     This is the whole reason `onComplete` takes a struct rather than a bare `std::error_code`:
     a `std::error_code` is only an (int, category) pair and cannot carry the text.  For the
     WinHTTP backend in particular, remote diagnosis is the only diagnosis available.

     Guaranteed never to contain the bearer token, any other request header value, or the
     request URL - a URL can carry credentials in its userinfo, or a token in its query string.
     */
    std::string detail;

    /** Which backend actually ran the request. */
    Backend backend = Backend::None;
  };//struct Completion


  /** Callbacks for one request.  Any of them may be empty.

   Ordering is guaranteed: `onHeaders` (at most once) precedes every `onChunk`, and
   `onComplete` follows all of them and fires exactly once.  `onHeaders` is skipped entirely if
   the request fails before a response line is received.

   See `start()` for which thread these arrive on.
   */
  struct StreamHandler
  {
    /** The HTTP status and response headers.  A non-2xx status is *not* itself an error: the
     body is still delivered, because providers report errors in-band with a 200 as readily as
     with a 4xx, and the caller needs the body either way. */
    std::function<void(int status, const HeaderList &headers)> onHeaders;

    /** A slice of the response body.  The `string_view` is only valid for the duration of the
     call - copy anything you need to keep.  Chunks arrive on arbitrary byte boundaries and
     carry no relationship to any framing in the payload.

     Return false to abort the request; `onComplete` then reports `Error::AbortedByHandler`.

     Note that the abort is not synchronous when a session id was given: the handler runs on the
     session thread while the transfer continues on a backend thread, so the request is torn
     down at the next chunk boundary rather than the instant false is returned.  A few more
     chunks may therefore still be delivered.  Use `Call::cancel()` if you need to stop a
     transfer promptly.
     */
    std::function<bool(std::string_view chunk)> onChunk;

    std::function<void(const Completion &)> onComplete;
  };//struct StreamHandler


  /** A request in flight.

   Destroying a `Call` cancels the underlying request.  No callback will *begin* after the
   destructor returns; one already executing runs to completion.  Cancelling explicitly via
   `cancel()` instead *does* still deliver `onComplete`, with `Error::Cancelled`.  The split is
   deliberate: tearing down the owner (say, clearing the pending-request map) wants silence,
   whereas a user pressing Stop wants the conversation finalized.

   Destroy the `Call` on the session thread and the stronger "no callback runs at all after the
   destructor" property holds for free, because the session lock serializes the destructor
   against delivery.  That is the normal case.  With an empty session id, callbacks run on a
   backend thread and the caller must not destroy the `Call` concurrently with one.

   The destructor does not wait for the backend to wind down; it hands the platform object over
   to the backend's own thread to release.  That is what keeps it non-blocking, and it is what
   makes it safe for backends whose teardown is only valid on their own thread.

   The destructor does not block.  It flips an atomic flag shared with the backend, which keeps
   its own reference to the internal state and drops any callback that sees the flag set.  A
   destructor that waited on a backend thread could deadlock against that same thread trying to
   post into the Wt session it is being destroyed on.

   One `Call` handles exactly one request.  This is what lets the `Wt::Http::Client` backend
   work at all - a single `Wt::Http::Client` handles one request at a time - and it suits the
   caller, since `LlmInterface` can genuinely have several requests outstanding at once when
   sub-agents run in parallel.
   */
  namespace Detail
  {
    /** State shared by a `Call` and the backend running its transfer.  Defined in the private
     header `InterSpec/NativeHttpClientImpl.h`; opaque to callers. */
    struct RequestState;
  }//namespace Detail


  class InterSpec_API Call
  {
  public:
    ~Call();

    /** Cancel the request.  Safe to call from any thread, and safe to call more than once or
     after the request has already finished (both are no-ops).  Unless the request had already
     completed, `onComplete` is delivered with `Error::Cancelled`. */
    void cancel();

    Call( const Call & ) = delete;
    Call &operator=( const Call & ) = delete;

  private:
    explicit Call( std::shared_ptr<Detail::RequestState> state );

    /** Shared with the backend so in-flight callbacks stay safe after this `Call` is gone. */
    std::shared_ptr<Detail::RequestState> m_state;

    friend std::unique_ptr<Call> start( const Request &, StreamHandler, const std::string & );
  };//class Call


  /** Whether a backend is compiled in and usable on this platform.

   False does not only mean "no backend for this platform".  On Linux the backend is
   `Wt::Http::Client`, which can only do `https://` if Wt itself was built with SSL support -
   and the prefix build scripts default that off.  Callers should check this and fall back to
   the browser path rather than issuing a request that can only fail.
   */
  InterSpec_API bool available();

  /** Which backend is compiled in, regardless of whether it is presently usable. */
  InterSpec_API Backend backend();

  /** A stable diagnostic identifier for the compiled-in backend, e.g. "NSURLSession".

   Worth surfacing in the LLM configuration UI, since which backend is in use is the first thing
   wanted in any report of a network problem - but note this is deliberately NOT localized.  It
   is a fixed technical name, like a class name; a UI showing it should wrap it in a translated
   sentence rather than translating the name itself.
   */
  InterSpec_API const char *backendName();

  /** Whether this backend honours `Request::proxyUrl`.  False for `Wt::Http::Client`. */
  InterSpec_API bool supportsProxyOverride();

  /** Whether this backend honours `Request::caBundlePath`.  False for WinHTTP, which validates
   against the Windows certificate store - where a corporate root will already be installed. */
  InterSpec_API bool supportsCaBundle();


  /** Issue a request.

   @param req The request; copied, so the caller need not keep it alive.
   @param handler Callbacks; see `StreamHandler` for the ordering guarantees.
   @param sessionId A Wt session id, i.e. `wApp->sessionId()`.

   When `sessionId` is non-empty, every callback is delivered on that session's thread via
   `Wt::WServer::post()`, in order.  Callers therefore need no locking and no marshalling of
   their own, and may touch widgets directly from a callback.  This mirrors what
   `Wt::Http::Client` already does when constructed in session context, and doing it here for
   every backend is what keeps the three uniform.  If the session ends before the request
   finishes, the remaining callbacks are simply dropped.

   When `sessionId` is empty, callbacks arrive on a backend-owned thread instead.  That is for
   unit tests and any non-Wt context; it is not appropriate for UI code.

   Never returns null.  A request that cannot even be started - no backend, malformed URL, an
   escape hatch this backend does not support - reports that through `onComplete` like any
   other failure, so the caller has a single path to handle.

   @returns A handle that cancels the request when destroyed.
   */
  InterSpec_API std::unique_ptr<Call> start( const Request &req, StreamHandler handler,
                                           const std::string &sessionId );
}//namespace NativeHttp


namespace std
{
  /** Lets `NativeHttp::Error` convert implicitly to `std::error_code`. */
  template<>
  struct is_error_code_enum<NativeHttp::Error> : true_type {};
}//namespace std

#endif //NATIVE_HTTP_CLIENT_H
