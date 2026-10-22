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

/** Makes one request with the native HTTP client (`NativeHttp`, whichever backend this build has),
 and prints how it went as one line of JSON; driven by tests/native_http_test.py.

 Usage: test_native_http --url <url> [--method POST|GET] [--body <text>] [--ca <file|dir>]
          [--insecure] [--proxy <url>] [--idle <s>] [--total <s>] [--max-size <bytes>]
          [--cancel-after-bytes <n>] [--abort-after-bytes <n>] [--follow]
 */

#include <mutex>
#include <chrono>
#include <string>
#include <thread>
#include <cstdlib>
#include <iostream>
#include <condition_variable>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "InterSpec/NativeHttpClient.h"

using namespace std;
using json = nlohmann::json;

namespace
{
  const char *error_name( const NativeHttp::Error e )
  {
    switch( e )
    {
      case NativeHttp::Error::Ok:                      return "Ok";
      case NativeHttp::Error::Cancelled:               return "Cancelled";
      case NativeHttp::Error::AbortedByHandler:        return "AbortedByHandler";
      case NativeHttp::Error::HostNotFound:            return "HostNotFound";
      case NativeHttp::Error::ConnectFailed:           return "ConnectFailed";
      case NativeHttp::Error::TlsHandshakeFailed:      return "TlsHandshakeFailed";
      case NativeHttp::Error::TlsCertificateUntrusted: return "TlsCertificateUntrusted";
      case NativeHttp::Error::ProxyUnreachable:        return "ProxyUnreachable";
      case NativeHttp::Error::ProxyAuthRequired:       return "ProxyAuthRequired";
      case NativeHttp::Error::Timeout:                 return "Timeout";
      case NativeHttp::Error::IdleTimeout:             return "IdleTimeout";
      case NativeHttp::Error::ResponseTooLarge:        return "ResponseTooLarge";
      case NativeHttp::Error::OptionUnsupported:       return "OptionUnsupported";
      case NativeHttp::Error::BackendUnavailable:      return "BackendUnavailable";
      case NativeHttp::Error::Unknown:                 return "Unknown";
    }
    return "?";
  }//error_name(...)
}//namespace


int main( int argc, char **argv )
{
  NativeHttp::Request req;
  req.method = "GET";
  long long cancel_after = -1, abort_after = -1;

  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    const auto next = [&]() -> string {
      if( (i + 1) >= argc )
      {
        cerr << "Missing value for " << arg << endl;
        exit( 2 );
      }
      return argv[++i];
    };

    if( arg == "--url" )                     req.url = next();
    else if( arg == "--method" )             req.method = next();
    else if( arg == "--body" )               req.body = next();
    else if( arg == "--ca" )                 req.caBundlePath = next();
    else if( arg == "--insecure" )           req.disableCertVerification = true;
    else if( arg == "--proxy" )              req.proxyUrl = next();
    else if( arg == "--idle" )               req.idleTimeout = chrono::seconds( atoi( next().c_str() ) );
    else if( arg == "--total" )              req.totalTimeout = chrono::seconds( atoi( next().c_str() ) );
    else if( arg == "--max-size" )           req.maxResponseSize = static_cast<size_t>( atoll( next().c_str() ) );
    else if( arg == "--cancel-after-bytes" ) cancel_after = atoll( next().c_str() );
    else if( arg == "--abort-after-bytes" )  abort_after = atoll( next().c_str() );
    else if( arg == "--follow" )             req.followRedirects = true;
    else
    {
      cerr << "Unknown argument " << arg << endl;
      return 2;
    }
  }//for( loop over arguments )

  mutex result_mutex;
  condition_variable done_cv;
  bool done = false;
  json result = json::object();
  string body;
  int nchunks = 0;
  bool cancel_now = false;

  const chrono::steady_clock::time_point start = chrono::steady_clock::now();

  NativeHttp::StreamHandler handler;
  handler.onHeaders = [&]( int status, const NativeHttp::HeaderList &headers ){
    lock_guard<mutex> lock( result_mutex );
    result["headerStatus"] = status;
    json hdrs = json::object();
    for( const auto &h : headers )
      hdrs[h.first] = h.second;
    result["headers"] = hdrs;
  };
  handler.onChunk = [&]( string_view chunk ) -> bool {
    lock_guard<mutex> lock( result_mutex );
    body.append( chunk.data(), chunk.size() );
    ++nchunks;
    if( (cancel_after >= 0) && (static_cast<long long>( body.size() ) >= cancel_after) )
      cancel_now = true;
    return !((abort_after >= 0) && (static_cast<long long>( body.size() ) >= abort_after));
  };
  handler.onComplete = [&]( const NativeHttp::Completion &c ){
    lock_guard<mutex> lock( result_mutex );
    result["error"] = error_name( static_cast<NativeHttp::Error>( c.ec.value() ) );
    result["status"] = c.status;
    result["detail"] = c.detail;
    done = true;
    done_cv.notify_all();
  };

  // No Wt session, so the callbacks run on the backend's own thread.
  unique_ptr<NativeHttp::Call> call = NativeHttp::start( req, handler, "" );

  {
    unique_lock<mutex> lock( result_mutex );
    while( !done )
    {
      done_cv.wait_for( lock, chrono::milliseconds( 50 ) );
      if( cancel_now && !done )
      {
        cancel_now = false;
        lock.unlock();
        call->cancel();
        lock.lock();
      }
    }
  }

  // The backend thread may still be finishing up after onComplete.
  this_thread::sleep_for( chrono::milliseconds( 100 ) );

  lock_guard<mutex> lock( result_mutex );
  result["backend"] = NativeHttp::backendName();
  result["body"] = body.substr( 0, 2000 );
  result["bodyBytes"] = body.size();
  result["chunks"] = nchunks;
  result["elapsedMs"] = chrono::duration_cast<chrono::milliseconds>( chrono::steady_clock::now() - start ).count();
  cout << result.dump( -1, ' ', false, json::error_handler_t::replace ) << endl;
  return 0;
}//main(...)
