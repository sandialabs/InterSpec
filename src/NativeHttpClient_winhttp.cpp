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

 - **Cancellation closes the handle.**  A synchronous `WinHttpReceiveResponse` /
   `WinHttpReadData` cannot be interrupted any other way; closing the request handle from another
   thread makes the in-flight call fail with `ERROR_WINHTTP_OPERATION_CANCELLED`, which the
   worker recognizes.  Everything touching the handle therefore goes through `m_mutex`.

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

    case ERROR_WINHTTP_LOGIN_FAILURE:
      return Error::ProxyAuthRequired;

    case ERROR_WINHTTP_UNRECOGNIZED_SCHEME:
    case ERROR_WINHTTP_INVALID_URL:
      return Error::Unknown;

    default:
      break;
  }//switch( code )

  return Error::Unknown;
}//map_win_error(...)


/** The handles for one in-flight request, plus the worker thread running it.

 `m_mutex` guards the handles specifically so that `cancelBackendRequest`, running on the session
 thread, can close the request handle out from under a blocking call on the worker thread without
 racing the worker's own close.
 */
struct WinHttpRequestData
{
  std::mutex m_mutex;
  HINTERNET m_session = nullptr;
  HINTERNET m_connect = nullptr;
  HINTERNET m_request = nullptr;
  std::thread m_worker;

  ~WinHttpRequestData()
  {
    // The worker owns the handles' lifetime; by the time it has been joined they are closed.
    if( m_worker.joinable() )
      m_worker.detach();
  }

  /** Close the request handle to interrupt a blocking WinHTTP call on the worker thread.  Safe
   to call from any thread, and safe if the worker has already closed it. */
  void cancel()
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    if( m_request )
    {
      WinHttpCloseHandle( m_request );
      m_request = nullptr;
    }
  }

  void closeAll()
  {
    std::lock_guard<std::mutex> lock( m_mutex );
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

  const std::wstring host( parts.lpszHostName, parts.dwHostNameLength );
  std::wstring path( parts.lpszUrlPath, parts.dwUrlPathLength );
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
    accessType = WINHTTP_ACCESS_TYPE_NAMED_PROXY;
    proxyName = widen( req.proxyUrl );
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

  {
    std::lock_guard<std::mutex> lock( data->m_mutex );
    data->m_session = session;
  }

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

  {
    std::lock_guard<std::mutex> lock( data->m_mutex );
    data->m_connect = connect;
  }

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

  {
    std::lock_guard<std::mutex> lock( data->m_mutex );
    data->m_request = request;
  }

  // Redirects: WinHTTP follows them by default, so turn that off unless asked for - the other
  //  backends do not follow, and an LLM endpoint that redirects is usually a misconfiguration.
  if( !req.followRedirects )
  {
    DWORD policy = WINHTTP_DISABLE_REDIRECTS;
    WinHttpSetOption( request, WINHTTP_OPTION_DISABLE_FEATURE, &policy, sizeof(policy) );
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
    headerBlock += widen( kv.first ) + L": " + widen( kv.second ) + L"\r\n";

  if( state->shouldStop() )
  {
    // Cancelled between start() and here; whoever cancelled owns the completion.
    data->closeAll();
    return;
  }

  BOOL ok = WinHttpSendRequest( request,
                                headerBlock.empty() ? WINHTTP_NO_ADDITIONAL_HEADERS
                                                    : headerBlock.c_str(),
                                headerBlock.empty() ? 0 : static_cast<DWORD>(-1),
                                req.body.empty() ? WINHTTP_NO_REQUEST_DATA
                                                 : const_cast<char *>(req.body.data()),
                                static_cast<DWORD>(req.body.size()),
                                static_cast<DWORD>(req.body.size()), 0 );
  if( !ok )
  {
    failWin( "Could not send the request", GetLastError() );
    return;
  }

  if( !WinHttpReceiveResponse( request, nullptr ) )
  {
    failWin( "No response from the server", GetLastError() );
    return;
  }

  DWORD statusCode = 0, statusSize = sizeof(statusCode);
  if( !WinHttpQueryHeaders( request,
                            WINHTTP_QUERY_STATUS_CODE | WINHTTP_QUERY_FLAG_NUMBER,
                            WINHTTP_HEADER_NAME_BY_INDEX, &statusCode, &statusSize,
                            WINHTTP_NO_HEADER_INDEX ) )
  {
    failWin( "Could not read the response status", GetLastError() );
    return;
  }

  // A 407 that got this far means WinHTTP could not satisfy the proxy challenge automatically -
  //  i.e. the SSO path above did not apply.  Say so specifically; "HTTP 407" on its own sends
  //  people looking in the wrong place.
  if( statusCode == HTTP_STATUS_PROXY_AUTH_REQ )
  {
    fail( Error::ProxyAuthRequired,
          "The network proxy requires credentials that Windows could not supply automatically."
          " Ask IT whether this machine is expected to authenticate to the proxy." );
    return;
  }

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
  for( ; ; )
  {
    if( state->shouldStop() )
      break;

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
  try
  {
    data->m_worker = std::thread( [state, data](){
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
    } );
  }catch( const std::exception &e )
  {
    // Thread creation failed; nothing will ever report, so report here.
    state->deliverComplete( Error::Unknown,
                            std::string("Could not start a thread for the request: ") + e.what() );
  }
}//startBackendRequest(...)


void cancelBackendRequest( const std::shared_ptr<RequestState> &state )
{
  std::shared_ptr<WinHttpRequestData> data;

  {
    std::lock_guard<std::mutex> lock( state->backendMutex );
    data = std::static_pointer_cast<WinHttpRequestData>( state->backendData );
  }

  if( !data )
    return;

  // Closing the request handle is the only way to interrupt a blocking WinHTTP call; the worker
  //  sees ERROR_WINHTTP_OPERATION_CANCELLED (or notices stopRequested) and unwinds.  Deliberately
  //  not joining the worker: this runs on the session thread, and blocking it on a socket that
  //  may take moments to unwind would stall the UI.  The worker holds its own shared_ptr to both
  //  the state and the handles, so letting it finish on its own is safe.
  data->cancel();
}//cancelBackendRequest(...)

}//namespace Detail

}//namespace NativeHttp
