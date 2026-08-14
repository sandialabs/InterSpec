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
#include <string>
#include <thread>
#include <memory>
#include <vector>
#include <sstream>
#include <iostream>

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
// Without this, windows.h defines min/max as macros and any later <algorithm>/<chrono> use breaks.
//  It survives a mingw build only because the standard headers above happen to come first there;
//  MSVC's transitive include set differs, which is exactly the class of break a cross-compile
//  check cannot catch.
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <winhttp.h>

#include "InterSpec/NativeHttpClient.h"
#include "InterSpec/NativeHttpClientImpl.h"

static_assert( USE_NATIVE_HTTP_CLIENT,
              "This file should only be compiled when USE_NATIVE_HTTP_CLIENT is enabled" );

/** The Windows backend, built on WinHTTP.

 WinHTTP ships with every Windows install, so this costs no third-party dependency.  It is also
 the strongest of the three backends on a managed network: it performs WPAD and PAC evaluation
 itself, validates through Schannel against the real Windows certificate store (so a
 group-policy-deployed corporate root just works), and can answer a proxy's Negotiate/NTLM
 challenge with the logged-on user's credentials without prompting.

 Structural choices worth knowing:

 - **Synchronous WinHTTP on a worker thread**, not `WINHTTP_FLAG_ASYNC`.  Results are marshalled
   onto the Wt session thread by `RequestState::dispatch()` regardless, so the async callback
   machinery would add a second state machine for no benefit.  One thread per request is fine
   here: requests are seconds-to-minutes of mostly-idle waiting, and there are only ever a
   handful outstanding.

 - **Cancellation is cooperative, not a handle close.**  It is tempting to close the request
   handle from another thread to break a blocked `WinHttpReadData`, and it usually appears to
   work - but Microsoft documents it as invalid: "An application should never call
   WinHttpCloseHandle on a synchronous request.  This can create a race condition", and separately
   "these HINTERNET handles cannot be closed while an API call using the handle is in progress".
   A mutex around the *handle variable* does not help, because the hazard is the *call* holding
   it.  So the worker polls `shouldStop()` at each chunk boundary and unwinds itself, and the
   handles are only ever touched by the thread that owns them.

   The cost is that a cancel arriving while the worker is blocked waiting for bytes is not acted
   on until the next chunk arrives, or the receive timeout expires.  That is not a user-visible
   hang: `Call::cancel()` delivers `onComplete` itself and `~Call()` is silent by design, so the
   conversation is finalized immediately either way - only the worker lingers, holding nothing the
   caller can see.  The public interface already promises no better than next-chunk-boundary
   teardown for the handler-abort case, for the same reason.

 - Everything in the API is UTF-16, so `widen()`/`narrow()` are used at every boundary.

 This backend was written without a Windows machine to test on, so it errs toward being explicit:
 every failure carries the `GetLastError()` code and a mapped message, because remote diagnosis
 is the only diagnosis available for it.

 It has been compile- and link-checked for both architectures without a Windows box, using
 mingw-w64 (which ships the same `winhttp.h` and a `libwinhttp.a` import library), e.g.:

 \code
   brew install mingw-w64
   x86_64-w64-mingw32-g++ -std=c++17 -c -Wall -Wextra -D_WIN32_WINNT=0x0601 \
     -I<build-dir> -I<repo-root> src/NativeHttpClient_winhttp.cpp -o /tmp/winhttp.o
 \endcode

 That catches wrong signatures, missing includes and unexported functions, and is worth re-running
 after any edit here.  It does NOT substitute for a run on real hardware: it compiles with GCC
 rather than MSVC, and nothing about Schannel, proxy discovery or SSO is exercised.
 */

namespace
{
using NativeHttp::Error;

/** UTF-8 -> UTF-16, for handing strings to WinHTTP. */
std::wstring widen( const std::string &str )
{
  if( str.empty() )
    return std::wstring();

  const int needed = MultiByteToWideChar( CP_UTF8, 0, str.c_str(),
                                          static_cast<int>(str.size()), nullptr, 0 );
  if( needed <= 0 )
    return std::wstring();

  std::wstring out( static_cast<size_t>(needed), L'\0' );
  MultiByteToWideChar( CP_UTF8, 0, str.c_str(), static_cast<int>(str.size()),
                       &out[0], needed );
  return out;
}//widen(...)


/** UTF-16 -> UTF-8, for handing strings back to the rest of InterSpec. */
std::string narrow( const wchar_t * const str, const size_t len )
{
  if( !str || !len )
    return std::string();

  const int needed = WideCharToMultiByte( CP_UTF8, 0, str, static_cast<int>(len),
                                          nullptr, 0, nullptr, nullptr );
  if( needed <= 0 )
    return std::string();

  std::string out( static_cast<size_t>(needed), '\0' );
  WideCharToMultiByte( CP_UTF8, 0, str, static_cast<int>(len), &out[0], needed,
                       nullptr, nullptr );
  return out;
}//narrow(...)


std::string narrow( const std::wstring &str )
{
  return narrow( str.c_str(), str.size() );
}


/** A readable description of a Win32/WinHTTP error code.

 WinHTTP's own codes live in winhttp.dll's message table, which FormatMessage will not find
 unless pointed at that module - hence the explicit `hWinhttp`.  The numeric code is always
 appended: it is stable across locales and is what someone actually searches for.
 */
std::string win_error_message( const DWORD code )
{
  LPWSTR buffer = nullptr;
  const HMODULE hWinhttp = GetModuleHandleW( L"winhttp.dll" );

  DWORD flags = FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM
                | FORMAT_MESSAGE_IGNORE_INSERTS;
  if( hWinhttp )
    flags |= FORMAT_MESSAGE_FROM_HMODULE;

  const DWORD len = FormatMessageW( flags, hWinhttp, code,
                                    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                                    reinterpret_cast<LPWSTR>(&buffer), 0, nullptr );

  std::string text;
  if( len && buffer )
  {
    // FormatMessage likes to append CR/LF; trim it so the message sits inline.
    size_t end = len;
    while( end && ((buffer[end-1] == L'\r') || (buffer[end-1] == L'\n')
                   || (buffer[end-1] == L' ')) )
      --end;
    text = narrow( buffer, end );
  }

  if( buffer )
    LocalFree( buffer );

  std::ostringstream ss;
  if( text.empty() )
    ss << "Windows error " << code;
  else
    ss << text << " [" << code << "]";
  return ss.str();
}//win_error_message(...)


/** Map a WinHTTP error onto our coarse set, so the user is told something actionable. */
Error map_win_error( const DWORD code )
{
  switch( code )
  {
    case ERROR_WINHTTP_NAME_NOT_RESOLVED:
      return Error::HostNotFound;

    case ERROR_WINHTTP_CANNOT_CONNECT:
    case ERROR_WINHTTP_CONNECTION_ERROR:
      return Error::ConnectFailed;

    case ERROR_WINHTTP_TIMEOUT:
      return Error::IdleTimeout;

    case ERROR_WINHTTP_OPERATION_CANCELLED:
      return Error::Cancelled;

    // The case this backend exists for on a managed network: a TLS-inspecting proxy re-signs
    //  certificates with a root that has to be in the Windows certificate store.
    case ERROR_WINHTTP_SECURE_FAILURE:
      return Error::TlsCertificateUntrusted;

    case ERROR_WINHTTP_SECURE_CERT_REV_FAILED:
    case ERROR_WINHTTP_SECURE_INVALID_CA:
    case ERROR_WINHTTP_SECURE_CERT_CN_INVALID:
    case ERROR_WINHTTP_SECURE_CERT_DATE_INVALID:
      return Error::TlsCertificateUntrusted;

    case ERROR_WINHTTP_CLIENT_AUTH_CERT_NEEDED:
    case ERROR_WINHTTP_SECURE_CHANNEL_ERROR:
      return Error::TlsHandshakeFailed;

    // Deliberately NOT ProxyAuthRequired: this covers a *server* 401 as well as a proxy 407, and
    //  telling someone with a stale API key to go ask IT about the proxy wastes their afternoon.
    case ERROR_WINHTTP_LOGIN_FAILURE:
      return Error::Unknown;

    case ERROR_WINHTTP_SHUTDOWN:
    case ERROR_WINHTTP_INTERNAL_ERROR:
      return Error::Unknown;

    default:
      break;
  }//switch( code )

  return Error::Unknown;
}//map_win_error(...)


/** The handles for one in-flight request.

 Touched only by the worker thread that created them - see the note on cancellation in the file
 header for why nothing else is allowed near them.  The worker is detached at creation rather than
 stored, so there is no `std::thread` here whose destructor could be reached from the worker
 itself.
 */
struct WinHttpRequestData
{
  HINTERNET m_session = nullptr;
  HINTERNET m_connect = nullptr;
  HINTERNET m_request = nullptr;

  ~WinHttpRequestData()
  {
    closeAll();
  }

  /** Only ever called on the worker thread. */
  void closeAll()
  {
    if( m_request ){ WinHttpCloseHandle( m_request ); m_request = nullptr; }
    if( m_connect ){ WinHttpCloseHandle( m_connect ); m_connect = nullptr; }
    if( m_session ){ WinHttpCloseHandle( m_session ); m_session = nullptr; }
  }
};//struct WinHttpRequestData


/** Split the raw CRLF-separated response header block into name/value pairs.

 WinHttpQueryHeaders(WINHTTP_QUERY_RAW_HEADERS_CRLF) hands back the status line followed by the
 headers, so the first line is skipped.
 */
NativeHttp::HeaderList parse_raw_headers( const std::wstring &raw )
{
  NativeHttp::HeaderList headers;

  size_t pos = raw.find( L"\r\n" );
  if( pos == std::wstring::npos )
    return headers;
  pos += 2;   // skip the status line

  while( pos < raw.size() )
  {
    size_t eol = raw.find( L"\r\n", pos );
    if( eol == std::wstring::npos )
      eol = raw.size();

    const std::wstring line = raw.substr( pos, eol - pos );
    pos = eol + 2;

    if( line.empty() )
      continue;

    const size_t colon = line.find( L':' );
    if( colon == std::wstring::npos )
      continue;

    std::wstring value = line.substr( colon + 1 );
    while( !value.empty() && ((value.front() == L' ') || (value.front() == L'\t')) )
      value.erase( value.begin() );

    headers.emplace_back( narrow( line.substr(0, colon) ), narrow( value ) );
  }//while( more header lines )

  return headers;
}//parse_raw_headers(...)


/** The whole transfer, run on a worker thread.

 Holds `data` and `state` by shared_ptr so both outlive the `Call` that started them.
 */
void run_request( std::shared_ptr<NativeHttp::Detail::RequestState> state,
                 std::shared_ptr<WinHttpRequestData> data )
{
  const NativeHttp::Request &req = state->request;

  // Anything that leaves this function must report exactly once; this makes that hard to forget.
  const auto fail = [&state,&data]( const Error err, const std::string &detail ){
    data->closeAll();
    state->deliverComplete( err, detail );
  };

  const auto failWin = [&fail]( const char * const what, const DWORD code ){
    fail( map_win_error(code), std::string(what) + ": " + win_error_message(code) );
  };

  URL_COMPONENTS parts;
  ZeroMemory( &parts, sizeof(parts) );
  parts.dwStructSize = sizeof(parts);
  parts.dwSchemeLength    = static_cast<DWORD>(-1);
  parts.dwHostNameLength  = static_cast<DWORD>(-1);
  parts.dwUrlPathLength   = static_cast<DWORD>(-1);
  parts.dwExtraInfoLength = static_cast<DWORD>(-1);

  const std::wstring wurl = widen( req.url );
  if( !WinHttpCrackUrl( wurl.c_str(), static_cast<DWORD>(wurl.size()), 0, &parts ) )
  {
    failWin( "Could not parse the endpoint URL", GetLastError() );
    return;
  }

  const std::wstring host = parts.dwHostNameLength
                              ? std::wstring( parts.lpszHostName, parts.dwHostNameLength )
                              : std::wstring();
  if( host.empty() )
  {
    fail( Error::Unknown, "The endpoint URL has no host" );
    return;
  }

  // Guarded: for "https://host" with no path WinHttpCrackUrl returns length 0 and may hand back a
  //  null pointer, and std::wstring(nullptr, 0) is undefined behaviour.
  std::wstring path = parts.dwUrlPathLength
                        ? std::wstring( parts.lpszUrlPath, parts.dwUrlPathLength )
                        : std::wstring();
  if( parts.dwExtraInfoLength )
    path.append( parts.lpszExtraInfo, parts.dwExtraInfoLength );
  if( path.empty() )
    path = L"/";

  const bool https = (parts.nScheme == INTERNET_SCHEME_HTTPS);

  // WINHTTP_ACCESS_TYPE_AUTOMATIC_PROXY (Win 8.1+) makes WinHTTP do WPAD and PAC evaluation
  //  itself, which is what a managed network needs and what we would otherwise have to
  //  reimplement via WinHttpGetIEProxyConfigForCurrentUser + WinHttpGetProxyForUrl.
  DWORD accessType = WINHTTP_ACCESS_TYPE_AUTOMATIC_PROXY;
  std::wstring proxyName;
  if( !req.proxyUrl.empty() )
  {
    // Normalize to a bare "server[:port]".
    //
    // Microsoft's grammar for WINHTTP_ACCESS_TYPE_NAMED_PROXY reads
    //  "([<scheme>=][<scheme>"://"]<server>[":"<port>])", i.e. the scheme parts are optional -
    //  but they are only optional in the sense that leaving them off always works.  Testing under
    //  Wine, "http://127.0.0.1:3129" failed immediately with ERROR_WINHTTP_NAME_NOT_RESOLVED and
    //  never reached the proxy at all, while "127.0.0.1:3129" tunnelled correctly (confirmed in
    //  the proxy's own access log).  Whether stock Windows is more forgiving is untested; the
    //  bare form is valid on both, so there is no reason to find out the hard way.
    //
    // This matters because "http://proxy.example.com:8080" is exactly the shape the settings
    //  tooltip tells the user to enter, and the shape anyone copying a proxy URL will paste.
    std::string proxy = req.proxyUrl;

    const size_t schemeEnd = proxy.find( "://" );
    if( schemeEnd != std::string::npos )
      proxy.erase( 0, schemeEnd + 3 );

    // Drop any path/query, and any "user@" - credentials cannot be passed this way.
    const size_t pathStart = proxy.find_first_of( "/?" );
    if( pathStart != std::string::npos )
      proxy.erase( pathStart );

    const size_t atPos = proxy.rfind( '@' );
    if( atPos != std::string::npos )
      proxy.erase( 0, atPos + 1 );

    if( proxy.empty() )
    {
      fail( Error::Unknown, "The proxy setting does not contain a server name" );
      return;
    }

    accessType = WINHTTP_ACCESS_TYPE_NAMED_PROXY;
    proxyName = widen( proxy );
  }

  HINTERNET session = WinHttpOpen( L"InterSpec", accessType,
                                   proxyName.empty() ? WINHTTP_NO_PROXY_NAME : proxyName.c_str(),
                                   WINHTTP_NO_PROXY_BYPASS, 0 );

  if( !session && (accessType == WINHTTP_ACCESS_TYPE_AUTOMATIC_PROXY) )
  {
    // Automatic mode needs Windows 8.1; fall back so this still works on older systems, where
    //  DEFAULT_PROXY picks up the machine-wide (netsh winhttp) configuration.
    accessType = WINHTTP_ACCESS_TYPE_DEFAULT_PROXY;
    session = WinHttpOpen( L"InterSpec", accessType, WINHTTP_NO_PROXY_NAME,
                           WINHTTP_NO_PROXY_BYPASS, 0 );
  }

  if( !session )
  {
    failWin( "Could not initialize WinHTTP", GetLastError() );
    return;
  }

  data->m_session = session;

  // WinHTTP's timeouts are per-phase.  Map the idle timeout onto the send/receive phases, which
  //  is the closest equivalent to "no progress for this long", and leave resolve/connect on
  //  something shorter so an unreachable host does not sit for two minutes.
  const int idleMs = (req.idleTimeout.count() > 0)
                       ? static_cast<int>(req.idleTimeout.count() * 1000) : 0;
  WinHttpSetTimeouts( session, 30000, 30000, idleMs ? idleMs : 120000, idleMs ? idleMs : 120000 );

  HINTERNET connect = WinHttpConnect( session, host.c_str(), parts.nPort, 0 );
  if( !connect )
  {
    failWin( "Could not connect", GetLastError() );
    return;
  }

  data->m_connect = connect;

  const DWORD requestFlags = (https ? WINHTTP_FLAG_SECURE : 0)
                             | WINHTTP_FLAG_ESCAPE_DISABLE_QUERY;

  HINTERNET request = WinHttpOpenRequest( connect, widen(req.method).c_str(), path.c_str(),
                                          nullptr, WINHTTP_NO_REFERER,
                                          WINHTTP_DEFAULT_ACCEPT_TYPES, requestFlags );
  if( !request )
  {
    failWin( "Could not create the request", GetLastError() );
    return;
  }

  data->m_request = request;

  // WinHttpSetTimeouts does NOT bound the wait for response headers; that is this separate
  //  option, and its default is 90 seconds.  Leaving it would cap time-to-first-byte at 90 s -
  //  which is exactly the case Request::idleTimeout exists to allow, since a non-streaming LLM
  //  emits no headers until the whole completion is ready.
  {
    DWORD headerWaitMs = (req.idleTimeout.count() > 0)
                           ? static_cast<DWORD>(req.idleTimeout.count() * 1000) : 300000;
    if( (req.totalTimeout.count() > 0)
       && (static_cast<DWORD>(req.totalTimeout.count() * 1000) > headerWaitMs) )
      headerWaitMs = static_cast<DWORD>(req.totalTimeout.count() * 1000);

    if( !WinHttpSetOption( request, WINHTTP_OPTION_RECEIVE_RESPONSE_TIMEOUT,
                           &headerWaitMs, sizeof(headerWaitMs) ) )
    {
      std::cerr << "NativeHttp: could not raise the response-header timeout ("
                << win_error_message( GetLastError() ) << "); slow completions may time out"
                << std::endl;
    }
  }

  // Redirects: WinHTTP follows them by default, so turn that off unless asked for - the other
  //  backends do not follow, and an LLM endpoint that redirects is usually a misconfiguration.
  //  The return value is checked because failing open here would replay the bearer token to
  //  whatever the redirect points at.
  if( !req.followRedirects )
  {
    DWORD policy = WINHTTP_DISABLE_REDIRECTS;
    if( !WinHttpSetOption( request, WINHTTP_OPTION_DISABLE_FEATURE, &policy, sizeof(policy) ) )
    {
      const DWORD code = GetLastError();
      fail( Error::OptionUnsupported,
            "Could not disable HTTP redirects, and following one would send the API key to the"
            " redirect target: " + win_error_message(code) );
      return;
    }
  }

  if( req.disableCertVerification )
  {
    // Explicitly opted into by the user, with a warning, for when nothing else works.
    DWORD securityFlags = SECURITY_FLAG_IGNORE_UNKNOWN_CA
                          | SECURITY_FLAG_IGNORE_CERT_DATE_INVALID
                          | SECURITY_FLAG_IGNORE_CERT_CN_INVALID
                          | SECURITY_FLAG_IGNORE_CERT_WRONG_USAGE;
    WinHttpSetOption( request, WINHTTP_OPTION_SECURITY_FLAGS,
                      &securityFlags, sizeof(securityFlags) );
  }

  // Single-sign-on against an intranet-zone proxy.  MEDIUM is the default and is deliberately
  //  kept: it offers the logged-on user's credentials to intranet servers only.  LOW would offer
  //  them to anything that challenges, including the origin server, which is not something to do
  //  by default just to save a round trip.
  DWORD logonPolicy = WINHTTP_AUTOLOGON_SECURITY_LEVEL_MEDIUM;
  WinHttpSetOption( request, WINHTTP_OPTION_AUTOLOGON_POLICY,
                    &logonPolicy, sizeof(logonPolicy) );

  std::wstring headerBlock;
  for( const std::pair<std::string,std::string> &kv : req.headers )
  {
    // Reject rather than sanitize: a CR or LF in a name or value would splice extra headers into
    //  the request, and these come from user-editable configuration.  Silently stripping would
    //  hide a broken config; failing names it.
    if( (kv.first.find_first_of("\r\n") != std::string::npos)
       || (kv.second.find_first_of("\r\n") != std::string::npos) )
    {
      fail( Error::Unknown, "Request header '" + kv.first + "' contains a line break" );
      return;
    }

    headerBlock += widen( kv.first ) + L": " + widen( kv.second ) + L"\r\n";
  }

  if( state->shouldStop() )
  {
    // Cancelled between start() and here; whoever cancelled owns the completion.
    data->closeAll();
    return;
  }

  // Send / receive, driving the authentication handshake.
  //
  // WinHTTP does not authenticate on its own.  The autologon policy set above only governs
  //  *whether default credentials may be used* once the application drives the loop; the loop
  //  itself - see a 401/407, ask which schemes are offered, set credentials, resend on the same
  //  handle - is the application's job, and without it a Negotiate/NTLM proxy simply fails on the
  //  first challenge.  This follows the pattern in Microsoft's "Authentication in WinHTTP".
  //
  // Bounded because a misconfigured proxy can otherwise challenge indefinitely.
  DWORD statusCode = 0;
  const int maxAuthAttempts = 5;

  for( int attempt = 0; ; ++attempt )
  {
    if( attempt >= maxAuthAttempts )
    {
      fail( Error::ProxyAuthRequired,
            "The network proxy kept asking for credentials that could not be supplied"
            " automatically.  Ask IT whether this machine is expected to authenticate to the"
            " proxy." );
      return;
    }

    if( state->shouldStop() )
    {
      data->closeAll();
      return;   // whoever set the flag owns the completion; see startBackendRequest
    }

    const BOOL sent = WinHttpSendRequest( request,
                                headerBlock.empty() ? WINHTTP_NO_ADDITIONAL_HEADERS
                                                    : headerBlock.c_str(),
                                headerBlock.empty() ? 0 : static_cast<DWORD>(-1),
                                req.body.empty() ? WINHTTP_NO_REQUEST_DATA
                                                 : const_cast<char *>(req.body.data()),
                                static_cast<DWORD>(req.body.size()),
                                static_cast<DWORD>(req.body.size()), 0 );
    if( !sent )
    {
      failWin( "Could not send the request", GetLastError() );
      return;
    }

    if( !WinHttpReceiveResponse( request, nullptr ) )
    {
      const DWORD code = GetLastError();
      // Documented on both SendRequest and ReceiveResponse: the handshake needs another round
      //  trip on this same handle.
      if( code == ERROR_WINHTTP_RESEND_REQUEST )
        continue;
      failWin( "No response from the server", code );
      return;
    }

    DWORD statusSize = sizeof(statusCode);
    if( !WinHttpQueryHeaders( request,
                              WINHTTP_QUERY_STATUS_CODE | WINHTTP_QUERY_FLAG_NUMBER,
                              WINHTTP_HEADER_NAME_BY_INDEX, &statusCode, &statusSize,
                              WINHTTP_NO_HEADER_INDEX ) )
    {
      failWin( "Could not read the response status", GetLastError() );
      return;
    }

    if( (statusCode != HTTP_STATUS_PROXY_AUTH_REQ) && (statusCode != HTTP_STATUS_DENIED) )
      break;   // a real response - including a 401 the caller should see, once we stop retrying

    // Which scheme to answer with, and whether we can answer at all.
    DWORD supported = 0, first = 0, target = 0;
    if( !WinHttpQueryAuthSchemes( request, &supported, &first, &target ) )
      break;   // nothing offered we can act on; let the caller see the 401/407 body

    // Preference order: Negotiate (Kerberos, then NTLM) is what a domain proxy uses, and is the
    //  only one that can succeed without prompting.  Basic/Digest would need credentials we do
    //  not have and must never invent, so they fall through to the caller as a plain 401/407.
    DWORD scheme = 0;
    if( supported & WINHTTP_AUTH_SCHEME_NEGOTIATE )
      scheme = WINHTTP_AUTH_SCHEME_NEGOTIATE;
    else if( supported & WINHTTP_AUTH_SCHEME_NTLM )
      scheme = WINHTTP_AUTH_SCHEME_NTLM;

    if( !scheme )
      break;

    // Null credentials = "use the logged-on user's", which is what the autologon policy gates.
    if( !WinHttpSetCredentials( request, target, scheme, nullptr, nullptr, nullptr ) )
    {
      failWin( "Could not use this machine's credentials for the proxy", GetLastError() );
      return;
    }
  }//for( send / authenticate )

  NativeHttp::HeaderList headers;
  {
    DWORD headerBytes = 0;
    WinHttpQueryHeaders( request, WINHTTP_QUERY_RAW_HEADERS_CRLF,
                         WINHTTP_HEADER_NAME_BY_INDEX, nullptr, &headerBytes,
                         WINHTTP_NO_HEADER_INDEX );
    if( (headerBytes > 0) && (GetLastError() == ERROR_INSUFFICIENT_BUFFER) )
    {
      std::wstring raw( headerBytes / sizeof(wchar_t), L'\0' );
      if( WinHttpQueryHeaders( request, WINHTTP_QUERY_RAW_HEADERS_CRLF,
                               WINHTTP_HEADER_NAME_BY_INDEX, &raw[0], &headerBytes,
                               WINHTTP_NO_HEADER_INDEX ) )
      {
        raw.resize( headerBytes / sizeof(wchar_t) );
        headers = parse_raw_headers( raw );
      }
    }
  }

  state->deliverHeaders( static_cast<int>(statusCode), std::move(headers) );

  // Body: WinHttpQueryDataAvailable / WinHttpReadData hand back de-chunked bytes, so what comes
  //  out here is the decoded entity body.
  std::vector<char> buffer;
  const std::chrono::steady_clock::time_point deadline
      = std::chrono::steady_clock::now() + req.totalTimeout;

  for( ; ; )
  {
    if( state->shouldStop() )
      break;

    // WinHTTP has no wall-clock ceiling of its own - its timeouts are all per-phase - so enforce
    //  Request::totalTimeout here, between reads.
    if( (req.totalTimeout.count() > 0) && (std::chrono::steady_clock::now() > deadline) )
    {
      fail( Error::Timeout, "The response did not finish within "
                            + std::to_string(req.totalTimeout.count()) + " seconds" );
      return;
    }

    DWORD available = 0;
    if( !WinHttpQueryDataAvailable( request, &available ) )
    {
      const DWORD code = GetLastError();
      if( state->stopRequested.load() )
        break;   // our own cancel closed the handle
      failWin( "Lost the connection while reading the response", code );
      return;
    }

    if( available == 0 )
      break;   // end of body

    if( buffer.size() < available )
      buffer.resize( available );

    DWORD read = 0;
    if( !WinHttpReadData( request, buffer.data(), available, &read ) )
    {
      const DWORD code = GetLastError();
      if( state->stopRequested.load() )
        break;
      failWin( "Could not read the response", code );
      return;
    }

    if( read == 0 )
      break;

    if( !state->deliverChunk( std::string_view( buffer.data(), read ) ) )
      break;   // cap exceeded, handler aborted, or cancelled - deliverChunk reports those
  }//for( read the whole body )

  data->closeAll();

  if( state->stopRequested.load() )
  {
    // Cancel and handler-abort both already have an owner for the completion; deliverComplete is
    //  latched anyway, so this is belt-and-braces rather than load-bearing.
    state->deliverComplete( state->stoppedByHandler.load() ? Error::AbortedByHandler
                                                           : Error::Cancelled,
                            "Cancelled by the application" );
    return;
  }

  state->deliverComplete( Error::Ok, std::string() );
}//run_request(...)

}//namespace


namespace NativeHttp
{

bool available()
{
  return true;
}

Backend backend()
{
  return Backend::WinHttp;
}

const char *backendName()
{
  return "WinHTTP";
}

bool supportsProxyOverride()
{
  return true;
}

bool supportsCaBundle()
{
  // Schannel validates against the Windows certificate store, and there is no supported way to
  //  hand WinHTTP an extra PEM bundle for server validation.  On Windows the answer to a
  //  corporate root is to install it in the store, which group policy normally already does.
  return false;
}


namespace Detail
{

void startBackendRequest( const std::shared_ptr<RequestState> &state )
{
  const std::shared_ptr<WinHttpRequestData> data = std::make_shared<WinHttpRequestData>();

  {
    std::lock_guard<std::mutex> lock( state->backendMutex );
    state->backendData = data;
  }

  // One thread per request, synchronous WinHTTP inside it.  startBackendRequest must not block,
  //  and the results are marshalled onto the session thread by dispatch() either way.
  //
  // Detached immediately rather than stored: keeping the std::thread inside `data` would close a
  //  cycle (state -> backendData -> data -> thread callable -> state), and when the owning Call
  //  goes first the worker ends up running ~WinHttpRequestData - and hence its own std::thread's
  //  destructor - from inside itself.
  try
  {
    std::thread( [state, data](){
      try
      {
        run_request( state, data );
      }catch( const std::exception &e )
      {
        data->closeAll();
        state->deliverComplete( Error::Unknown,
                                std::string("Unexpected failure in the HTTP worker: ") + e.what() );
      }catch( ... )
      {
        data->closeAll();
        state->deliverComplete( Error::Unknown, "Unexpected failure in the HTTP worker" );
      }

      // Break the cycle deterministically, and drop the request body - which carries the bearer
      //  token - as soon as the transfer is over rather than whenever the Call happens to die.
      {
        std::lock_guard<std::mutex> lock( state->backendMutex );
        state->backendData.reset();
      }
    } ).detach();
  }catch( const std::exception &e )
  {
    // Thread creation failed; nothing will ever report, so report here.
    state->deliverComplete( Error::Unknown,
                            std::string("Could not start a thread for the request: ") + e.what() );
  }
}//startBackendRequest(...)


void cancelBackendRequest( const std::shared_ptr<RequestState> & )
{
  // Nothing to do.  `stopRequested` is already set by the caller, and the worker polls
  //  `shouldStop()` at each chunk boundary and unwinds itself.
  //
  // Deliberately does NOT reach in and close the request handle: Microsoft documents that as a
  //  race on a synchronous request (see the note in this file's header), and a mutex cannot fix
  //  it because the hazard is the in-flight call, not the variable.  Nor does it join the worker
  //  - this runs on the session thread, and the completion has already been delivered by
  //  Call::cancel(), so there is nothing to wait for.
}//cancelBackendRequest(...)

}//namespace Detail

}//namespace NativeHttp
