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
#include <chrono>
#include <string>
#include <thread>
#include <memory>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <iostream>

#include <dlfcn.h>

// Only the declarations and constants are used: libcurl is loaded at run time (see curl_api()).
#include <curl/curl.h>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/NativeHttpClient.h"
#include "InterSpec/NativeHttpClientImpl.h"

static_assert( USE_NATIVE_HTTP_CLIENT,
              "This file should only be compiled when USE_NATIVE_HTTP_CLIENT is enabled" );

/*
 The libcurl backend, for Linux (and, for testing, macOS).

 libcurl is not linked: it is loaded with dlopen() when first needed, trying `libcurl.so.4` and
 then `libcurl-gnutls.so.4` (the same ABI).  Linking it would tie the binary to the symbol
 versions of the build machine's libcurl (a manylinux build would require AlmaLinux's
 `CURL_OPENSSL_4`, which Debian and Ubuntu do not provide) - the same kind of coupling the
 dependency prefix avoids for OpenSSL.  If no libcurl is found, `available()` is false.

 Each request runs `curl_easy_perform()` on its own thread, which works with libcurl back to 7.61
 (RHEL 8; `curl_multi_poll`/`curl_multi_wakeup` would need 7.66/7.68), and matches the
 one-request-per-`Call` model.  Cancellation is polled by the progress callback, which libcurl
 calls at least once a second.

 Unlike `Wt::Http::Client`, libcurl honours proxies: an explicit `proxyUrl`, or else the usual
 `https_proxy`, `ALL_PROXY` and `NO_PROXY` environment variables, with any proxy authentication
 it supports.  Constants added after 7.61 are used by value, and only when the loaded libcurl is
 new enough.
 */

namespace
{
  /** The libcurl functions used, resolved from the loaded library. */
  struct CurlApi
  {
    bool ok = false;
    unsigned int version_num = 0;   // e.g., 0x075300 for 7.83.0
    std::string version;            // e.g., "7.61.1"
    std::string load_error;

    CURLcode (*global_init)( long ) = nullptr;
    CURL *(*easy_init)() = nullptr;
    CURLcode (*easy_setopt)( CURL *, CURLoption, ... ) = nullptr;
    CURLcode (*easy_perform)( CURL * ) = nullptr;
    CURLcode (*easy_getinfo)( CURL *, CURLINFO, ... ) = nullptr;
    void (*easy_cleanup)( CURL * ) = nullptr;
    const char *(*easy_strerror)( CURLcode ) = nullptr;
    struct curl_slist *(*slist_append)( struct curl_slist *, const char * ) = nullptr;
    void (*slist_free_all)( struct curl_slist * ) = nullptr;
    curl_version_info_data *(*version_info)( CURLversion ) = nullptr;
  };//struct CurlApi


  // Values from libcurl versions newer than 7.61, whose headers the build may not have.
  const CURLINFO ns_info_cainfo = static_cast<CURLINFO>( CURLINFO_STRING + 61 );  // 7.84
  const int ns_curle_proxy = 97;                                                  // 7.73
  const CURLoption ns_opt_suppress_connect_headers = static_cast<CURLoption>( CURLOPTTYPE_LONG + 265 );  // 7.54


  template<class T>
  bool resolve( void *lib, const char *name, T &fcn )
  {
    fcn = reinterpret_cast<T>( dlsym( lib, name ) );
    return (fcn != nullptr);
  }


  const CurlApi &curl_api()
  {
    static CurlApi api;
    static std::once_flag once;

    std::call_once( once, [](){
      const char * const candidates[] = {
#if defined(__APPLE__)
        "/usr/lib/libcurl.4.dylib", "libcurl.4.dylib"
#else
        "libcurl.so.4", "libcurl-gnutls.so.4"
#endif
      };

      void *lib = nullptr;
      for( const char * const name : candidates )
      {
        lib = dlopen( name, RTLD_NOW | RTLD_LOCAL );
        if( lib )
          break;
      }

      if( !lib )
      {
        const char * const err = dlerror();
        api.load_error = std::string( "libcurl could not be loaded" ) + (err ? (std::string(": ") + err) : std::string());
        return;
      }

      const bool resolved = resolve( lib, "curl_global_init", api.global_init )
                            && resolve( lib, "curl_easy_init", api.easy_init )
                            && resolve( lib, "curl_easy_setopt", api.easy_setopt )
                            && resolve( lib, "curl_easy_perform", api.easy_perform )
                            && resolve( lib, "curl_easy_getinfo", api.easy_getinfo )
                            && resolve( lib, "curl_easy_cleanup", api.easy_cleanup )
                            && resolve( lib, "curl_easy_strerror", api.easy_strerror )
                            && resolve( lib, "curl_slist_append", api.slist_append )
                            && resolve( lib, "curl_slist_free_all", api.slist_free_all )
                            && resolve( lib, "curl_version_info", api.version_info );
      if( !resolved )
      {
        api.load_error = "the loaded libcurl is missing functions";
        return;
      }

      const curl_version_info_data * const info = api.version_info( CURLVERSION_NOW );
      api.version_num = info ? info->version_num : 0;
      api.version = (info && info->version) ? info->version : "unknown";

      if( api.version_num < 0x073D00 )  // 7.61.0
      {
        api.load_error = "libcurl " + api.version + " is too old (7.61 or newer is needed)";
        return;
      }

      if( api.global_init( CURL_GLOBAL_DEFAULT ) != CURLE_OK )
      {
        api.load_error = "curl_global_init failed";
        return;
      }

      api.ok = true;
    } );

    return api;
  }//curl_api()


  /** Where distributions keep their CA bundle; used when libcurl can not say (before 7.84). */
  const char * const ns_ca_bundles[] = {
    "/etc/ssl/certs/ca-certificates.crt",   // Debian/Ubuntu
    "/etc/pki/tls/certs/ca-bundle.crt",     // RHEL/Fedora
    "/etc/ssl/ca-bundle.pem",               // openSUSE
    "/etc/ssl/cert.pem",                    // Alpine, macOS
  };


  /** One transfer: the request's state, and what the libcurl callbacks need. */
  struct CurlTransfer
  {
    std::shared_ptr<NativeHttp::Detail::RequestState> state;

    // The response being received; reset by each status line (redirects, 1xx responses).
    int status = 0;
    NativeHttp::HeaderList headers;
    bool headers_delivered = false;

    char error_buffer[CURL_ERROR_SIZE] = { '\0' };

    /** A CA bundle written for this request (the system's plus the user's), removed afterwards. */
    std::string temp_ca_file;
  };//struct CurlTransfer


  size_t header_callback( char *buffer, size_t size, size_t nitems, void *userdata )
  {
    CurlTransfer &transfer = *static_cast<CurlTransfer *>( userdata );
    const size_t nbytes = size * nitems;
    std::string line( buffer, nbytes );
    while( !line.empty() && ((line.back() == '\n') || (line.back() == '\r')) )
      line.pop_back();

    if( SpecUtils::starts_with( line, "HTTP/" ) )
    {
      // A new response: "HTTP/1.1 200 OK", "HTTP/2 429"
      transfer.headers.clear();
      transfer.status = 0;
      const size_t space = line.find( ' ' );
      if( space != std::string::npos )
        transfer.status = std::atoi( line.c_str() + space + 1 );
    }else if( line.empty() )
    {
      // The end of a response's headers.  Informational responses, and redirects that libcurl will
      //  follow, are not the final response.
      const bool informational = (transfer.status < 200);
      const bool followed_redirect = transfer.state->request.followRedirects
                                     && (transfer.status >= 300) && (transfer.status < 400);
      if( !informational && !followed_redirect && !transfer.headers_delivered )
      {
        transfer.headers_delivered = true;
        transfer.state->deliverHeaders( transfer.status, transfer.headers );
      }
    }else
    {
      const size_t colon = line.find( ':' );
      if( colon != std::string::npos )
      {
        std::string name = line.substr( 0, colon ), value = line.substr( colon + 1 );
        SpecUtils::trim( name );
        SpecUtils::trim( value );
        transfer.headers.emplace_back( std::move(name), std::move(value) );
      }
    }//if( status line ) / else if( end of headers ) / else

    return nbytes;
  }//header_callback(...)


  size_t write_callback( char *data, size_t size, size_t nmemb, void *userdata )
  {
    CurlTransfer &transfer = *static_cast<CurlTransfer *>( userdata );
    const size_t nbytes = size * nmemb;

    // Returning less than nbytes makes libcurl stop with CURLE_WRITE_ERROR.
    if( !transfer.state->deliverChunk( std::string_view( data, nbytes ) ) )
      return 0;
    return nbytes;
  }//write_callback(...)


  int progress_callback( void *clientp, curl_off_t, curl_off_t, curl_off_t, curl_off_t )
  {
    const CurlTransfer &transfer = *static_cast<const CurlTransfer *>( clientp );
    return transfer.state->shouldStop() ? 1 : 0;   // non-zero aborts the transfer
  }//progress_callback(...)


  /** True if a proxy would be used: an explicit one, or one from the environment. */
  bool uses_proxy( const NativeHttp::Request &req )
  {
    if( !req.proxyUrl.empty() )
      return true;
    for( const char * const name : { "https_proxy", "HTTPS_PROXY", "http_proxy", "all_proxy", "ALL_PROXY" } )
    {
      const char * const value = std::getenv( name );
      if( value && value[0] )
        return true;
    }
    return false;
  }//uses_proxy(...)


  /** `text` with the request URL and header values replaced, so it can not leak a token. */
  std::string scrub( std::string text, const NativeHttp::Request &req )
  {
    std::vector<std::string> secrets{ req.url };
    for( const std::pair<std::string,std::string> &header : req.headers )
      secrets.push_back( header.second );

    for( const std::string &secret : secrets )
    {
      if( secret.size() < 6 )
        continue;
      for( size_t pos = text.find( secret ); pos != std::string::npos; pos = text.find( secret, pos ) )
        text.replace( pos, secret.size(), "<redacted>" );
    }
    return text;
  }//scrub(...)


  /** Writes the system CA bundle plus the user's `ca_file` to a temporary file, as libcurl's
   CURLOPT_CAINFO replaces the system trust store rather than adding to it.  Returns its path.
   */
  std::string write_combined_ca_bundle( const CurlApi &api, CURL *easy, const std::string &ca_file )
  {
    std::string system_bundle;
    if( api.version_num >= 0x075400 )  // 7.84: libcurl can say what it would use
    {
      char *path = nullptr;
      if( (api.easy_getinfo( easy, ns_info_cainfo, &path ) == CURLE_OK) && path )
        system_bundle = path;
    }
    if( system_bundle.empty() || !SpecUtils::is_file( system_bundle ) )
    {
      system_bundle.clear();
      for( const char * const bundle : ns_ca_bundles )
      {
        if( SpecUtils::is_file( bundle ) )
        {
          system_bundle = bundle;
          break;
        }
      }
    }//if( libcurl did not give the bundle )

    const std::string temp_path = SpecUtils::temp_file_name( "interspec_ca_bundle", SpecUtils::temp_dir() );
    std::ofstream out( temp_path.c_str(), std::ios::binary );
    for( const std::string &path : { system_bundle, ca_file } )
    {
      if( path.empty() )
        continue;
      std::ifstream in( path.c_str(), std::ios::binary );
      out << in.rdbuf() << "\n";
    }
    if( !out )
      throw std::runtime_error( "Could not write the combined CA bundle" );
    return temp_path;
  }//write_combined_ca_bundle(...)


  NativeHttp::Error map_curl_error( const CurlApi &api, CURL *easy, const CURLcode res,
                                    const NativeHttp::Request &req, const double elapsed_s )
  {
    using NativeHttp::Error;

    long connect_code = 0;
    api.easy_getinfo( easy, CURLINFO_HTTP_CONNECTCODE, &connect_code );
    if( connect_code == 407 )
      return Error::ProxyAuthRequired;

    switch( static_cast<int>( res ) )
    {
      case CURLE_COULDNT_RESOLVE_HOST:
        return Error::HostNotFound;

      case CURLE_COULDNT_RESOLVE_PROXY:
        return Error::ProxyUnreachable;

      case CURLE_COULDNT_CONNECT:
        return uses_proxy( req ) ? Error::ProxyUnreachable : Error::ConnectFailed;

      case CURLE_RECV_ERROR:
      case CURLE_SEND_ERROR:
      case CURLE_GOT_NOTHING:
        return Error::ConnectFailed;

      case 51:  // CURLE_PEER_FAILED_VERIFICATION before 7.62
      case 60:  // CURLE_SSL_CACERT before 7.62, CURLE_PEER_FAILED_VERIFICATION since
      case 77:  // CURLE_SSL_CACERT_BADFILE
        return Error::TlsCertificateUntrusted;

      case CURLE_SSL_CONNECT_ERROR:
        return Error::TlsHandshakeFailed;

      case CURLE_OPERATION_TIMEDOUT:
      {
        // Both limits report this; the total one is the one whose time has run out.
        const double total_s = static_cast<double>( req.totalTimeout.count() );
        return ((total_s > 0.0) && (elapsed_s >= (total_s - 1.0))) ? Error::Timeout : Error::IdleTimeout;
      }

      default:
        if( static_cast<int>( res ) == ns_curle_proxy )
          return Error::ProxyUnreachable;
        return Error::Unknown;
    }//switch( res )
  }//map_curl_error(...)


  /** Sets up and runs the transfer, on its own thread; reports through `state`. */
  void run_transfer( const std::shared_ptr<CurlTransfer> &transfer )
  {
    using NativeHttp::Error;

    const CurlApi &api = curl_api();
    NativeHttp::Detail::RequestState &state = *transfer->state;
    const NativeHttp::Request &req = state.request;

    CURL * const easy = api.easy_init();
    if( !easy )
    {
      state.deliverComplete( Error::Unknown, "curl_easy_init failed" );
      return;
    }

    struct curl_slist *header_list = nullptr;
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    CURLcode res = CURLE_OK;

    try
    {
      api.easy_setopt( easy, CURLOPT_URL, req.url.c_str() );
      api.easy_setopt( easy, CURLOPT_NOSIGNAL, 1L );
      api.easy_setopt( easy, CURLOPT_ERRORBUFFER, transfer->error_buffer );

      if( req.method == "GET" )
      {
        api.easy_setopt( easy, CURLOPT_HTTPGET, 1L );
      }else
      {
        if( req.method != "POST" )
          api.easy_setopt( easy, CURLOPT_CUSTOMREQUEST, req.method.c_str() );
        api.easy_setopt( easy, CURLOPT_POST, 1L );
        api.easy_setopt( easy, CURLOPT_POSTFIELDSIZE_LARGE, static_cast<curl_off_t>( req.body.size() ) );
        api.easy_setopt( easy, CURLOPT_POSTFIELDS, req.body.data() );
      }

      // No "Expect: 100-continue" round trip; and no Accept-Encoding, as nothing is decompressed.
      header_list = api.slist_append( header_list, "Expect:" );
      for( const std::pair<std::string,std::string> &header : req.headers )
        header_list = api.slist_append( header_list, (header.first + ": " + header.second).c_str() );
      api.easy_setopt( easy, CURLOPT_HTTPHEADER, header_list );

      api.easy_setopt( easy, CURLOPT_HEADERFUNCTION, &header_callback );
      api.easy_setopt( easy, CURLOPT_HEADERDATA, transfer.get() );
      api.easy_setopt( easy, CURLOPT_WRITEFUNCTION, &write_callback );
      api.easy_setopt( easy, CURLOPT_WRITEDATA, transfer.get() );
      api.easy_setopt( easy, CURLOPT_NOPROGRESS, 0L );
      api.easy_setopt( easy, CURLOPT_XFERINFOFUNCTION, &progress_callback );
      api.easy_setopt( easy, CURLOPT_XFERINFODATA, transfer.get() );

      if( req.followRedirects )
        api.easy_setopt( easy, CURLOPT_FOLLOWLOCATION, 1L );

      // An idle timeout: less than 1 byte/s for this long.  (It also covers a slow first byte.)
      if( req.idleTimeout.count() > 0 )
      {
        api.easy_setopt( easy, CURLOPT_LOW_SPEED_LIMIT, 1L );
        api.easy_setopt( easy, CURLOPT_LOW_SPEED_TIME, static_cast<long>( req.idleTimeout.count() ) );
      }
      if( req.totalTimeout.count() > 0 )
        api.easy_setopt( easy, CURLOPT_TIMEOUT, static_cast<long>( req.totalTimeout.count() ) );

      // Without an explicit proxy, libcurl uses https_proxy / ALL_PROXY / NO_PROXY.
      if( !req.proxyUrl.empty() )
        api.easy_setopt( easy, CURLOPT_PROXY, req.proxyUrl.c_str() );
      api.easy_setopt( easy, CURLOPT_PROXYAUTH, static_cast<long>( CURLAUTH_ANY ) );
      api.easy_setopt( easy, ns_opt_suppress_connect_headers, 1L );  // the proxy's CONNECT reply is not ours

      if( req.disableCertVerification )
      {
        api.easy_setopt( easy, CURLOPT_SSL_VERIFYPEER, 0L );
        api.easy_setopt( easy, CURLOPT_SSL_VERIFYHOST, 0L );
      }else if( !req.caBundlePath.empty() )
      {
        // Added to the system trust store, not replacing it.
        if( SpecUtils::is_directory( req.caBundlePath ) )
        {
          api.easy_setopt( easy, CURLOPT_CAPATH, req.caBundlePath.c_str() );
        }else
        {
          transfer->temp_ca_file = write_combined_ca_bundle( api, easy, req.caBundlePath );
          api.easy_setopt( easy, CURLOPT_CAINFO, transfer->temp_ca_file.c_str() );
        }
      }//if( no verification ) / else if( an extra CA bundle )

      res = api.easy_perform( easy );
    }catch( std::exception &e )
    {
      res = CURLE_FAILED_INIT;
      std::strncpy( transfer->error_buffer, e.what(), sizeof(transfer->error_buffer) - 1 );
    }

    const double elapsed_s = std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count();

    if( state.stopRequested.load() )
    {
      state.deliverComplete( state.stoppedByHandler.load() ? Error::AbortedByHandler : Error::Cancelled,
                             "Cancelled by the application" );
    }else if( res == CURLE_OK )
    {
      // A response with no headers callback for the final response (should not happen), still
      //  needs its status reported.
      if( !transfer->headers_delivered )
      {
        long code = 0;
        api.easy_getinfo( easy, CURLINFO_RESPONSE_CODE, &code );
        transfer->headers_delivered = true;
        state.deliverHeaders( static_cast<int>( code ), transfer->headers );
      }
      state.deliverComplete( Error::Ok, std::string() );
    }else
    {
      // (A body over the size cap has already been reported, by deliverChunk().)
      std::string detail = transfer->error_buffer[0] ? std::string( transfer->error_buffer )
                                                     : std::string( api.easy_strerror( res ) );
      detail = scrub( detail, req ) + " (libcurl " + api.version + ", error " + std::to_string( static_cast<int>(res) ) + ")";
      state.deliverComplete( map_curl_error( api, easy, res, req, elapsed_s ), detail );
    }//if( cancelled ) / else if( ok ) / else

    api.slist_free_all( header_list );
    api.easy_cleanup( easy );
    if( !transfer->temp_ca_file.empty() )
      SpecUtils::remove_file( transfer->temp_ca_file );
  }//run_transfer(...)
}//namespace


namespace NativeHttp
{

bool available()
{
  return curl_api().ok;
}

Backend backend()
{
  return Backend::Curl;
}

const char *backendName()
{
  static const std::string name = curl_api().ok ? ("libcurl " + curl_api().version)
                                                : ("libcurl (not available: " + curl_api().load_error + ")");
  return name.c_str();
}

bool supportsProxyOverride()
{
  return true;
}

bool supportsCaBundle()
{
  return true;
}


namespace Detail
{

void startBackendRequest( const std::shared_ptr<RequestState> &state )
{
  if( !curl_api().ok )
  {
    state->deliverComplete( Error::BackendUnavailable, curl_api().load_error );
    return;
  }

  std::shared_ptr<CurlTransfer> transfer = std::make_shared<CurlTransfer>();
  transfer->state = state;

  // Detached rather than stored: the thread owns the transfer (and so the state) until it ends,
  //  and cancellation reaches it through `state->stopRequested`, which the progress callback polls.
  try
  {
    std::thread( [transfer](){ run_transfer( transfer ); } ).detach();
  }catch( std::exception &e )
  {
    state->deliverComplete( Error::Unknown, std::string( "Could not start the request thread: " ) + e.what() );
  }
}//startBackendRequest(...)


void cancelBackendRequest( const std::shared_ptr<RequestState> & )
{
  // Nothing to do: `stopRequested` is already set, and the transfer's progress callback - called at
  //  least once a second - aborts it.
}

}//namespace Detail

}//namespace NativeHttp
