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
#include <memory>
#include <vector>
#include <iostream>
#include <functional>


// WConfig.h has to come first: it is what defines WT_WITH_SSL, which everything below keys off.
#include <Wt/WConfig.h>
#include <Wt/WServer.h>
#include <Wt/WIOService.h>
#include <Wt/Http/Client.h>
#include <Wt/Http/Message.h>

#ifdef WT_WITH_SSL
#include <openssl/err.h>
#include <openssl/ssl.h>
#include <boost/asio/post.hpp>
#include <boost/asio/ssl/error.hpp>
#endif

#include "SpecUtils/Filesystem.h"

#include "InterSpec/NativeHttpClient.h"
#include "InterSpec/NativeHttpClientImpl.h"

static_assert( USE_NATIVE_HTTP_CLIENT,
              "This file should only be compiled when USE_NATIVE_HTTP_CLIENT is enabled" );

/** The Linux and Android backend, built on `Wt::Http::Client`.

 Wt is already a dependency, so this costs no new library.  It is also the weakest of the three
 backends: it has no proxy support at all, and its trust store is whatever OpenSSL is pointed at.
 That is an accepted limitation - see the table in NativeHttpClient.h.

 Several things about `Wt::Http::Client` need care, and each is commented where it is handled:

 - It can only do `https://` if Wt was compiled with SSL support, which the prefix build scripts
   leave off by default.  `available()` reports that rather than letting requests fail obscurely.
 - `setMaximumResponseSize(0)` looks like the way to stream, but it makes Wt stop accumulating the
   body, and the "is the response complete" test is `response_.body().size() >= contentLength_`.
   With a `Content-Length` response - which is what a non-streaming LLM POST returns - that test
   never becomes true and the request hangs until the idle timeout.  So the cap is left non-zero.
 - `abort()` and destruction are not thread-safe, and must happen on an io_service thread.
 - `done()` is emitted while the client holds its implementation mutex, which `abort()` and
   `~Client()` also take - so the `done()` handler must not take any mutex held around those calls.
 - `request()` can return false or throw without ever emitting `done()`.
 - Signals must be connected before the request is started, or they silently never fire.
 */

namespace
{
/** Candidate system trust stores, most common first.

 A statically linked OpenSSL has its default certificate directory compiled in, pointing at the
 prefix it was built in - a path that will not exist on a user's machine.  So the verify path has
 to be set explicitly; relying on OpenSSL's defaults (or on SSL_CERT_DIR / SSL_CERT_FILE) does not
 work for the way we build it.
 */
const char * const sm_ca_directories[] = {
  "/etc/ssl/certs",                       // Debian, Ubuntu, Arch, Alpine, and most others
  "/etc/pki/tls/certs",                   // RHEL, Fedora, CentOS
  "/system/etc/security/cacerts",         // Android (but see the note in available())
};

const char * const sm_ca_bundles[] = {
  "/etc/ssl/certs/ca-certificates.crt",   // Debian/Ubuntu
  "/etc/pki/tls/certs/ca-bundle.crt",     // RHEL/Fedora
  "/etc/ssl/ca-bundle.pem",               // openSUSE
  "/etc/ssl/cert.pem",                    // Alpine, and macOS when testing this backend there
};


/** Everything one request needs, owned by the state's `backendData` slot. */
struct WtRequestData
{
  /** Raw rather than a smart pointer because it must be deleted on the io_service thread; see
   `cancelBackendRequest`.  `Wt::Http::Client` is a `WObject`, hence a boost::signals2
   `trackable`, whose auto-disconnect on destruction is not thread-safe. */
  Wt::Http::Client *client = nullptr;
};//struct WtRequestData


/** Point the client at whichever system trust store this machine actually has. */
void apply_trust_settings( Wt::Http::Client &client, const NativeHttp::Request &req )
{
  if( req.disableCertVerification )
  {
    client.setSslCertificateVerificationEnabled( false );
    return;
  }

  client.setSslCertificateVerificationEnabled( true );

  // A user-supplied bundle wins for whichever slot it occupies, but is still additive: OpenSSL
  //  consults both the verify file and the verify path, so adding a corporate root does not stop
  //  ordinary public certificates from validating.
  bool haveFile = false, havePath = false;

  if( !req.caBundlePath.empty() )
  {
    if( SpecUtils::is_directory( req.caBundlePath ) )
    {
      client.setSslVerifyPath( req.caBundlePath );
      havePath = true;
    }else
    {
      client.setSslVerifyFile( req.caBundlePath );
      haveFile = true;
    }
  }//if( a custom CA bundle was configured )

  // Set a directory *and* a bundle where both exist, rather than stopping at the first hit.
  //  Merely existing is not enough for a hash directory: macOS, for instance, ships an empty
  //  /etc/ssl/certs alongside a fully populated /etc/ssl/cert.pem, and stopping at the directory
  //  would silently leave OpenSSL with nothing to verify against.
  if( !havePath )
  {
    for( const char * const dir : sm_ca_directories )
    {
      if( SpecUtils::is_directory(dir) && !SpecUtils::ls_files_in_directory(dir).empty() )
      {
        client.setSslVerifyPath( dir );
        havePath = true;
        break;
      }
    }
  }//if( no explicit verify directory yet )

  if( !haveFile )
  {
    for( const char * const bundle : sm_ca_bundles )
    {
      if( SpecUtils::is_file( bundle ) )
      {
        client.setSslVerifyFile( bundle );
        haveFile = true;
        break;
      }
    }
  }//if( no explicit verify file yet )

  if( !haveFile && !havePath )
  {
    // Wt falls back to OpenSSL's compiled-in defaults here - which, for the static OpenSSL the
    //  prefix scripts build, point at the build machine's prefix and will not exist.  Every
    //  request would then fail with an unhelpful verification error, so say why up front.
    std::cerr << "NativeHttp: found no system certificate store (looked in /etc/ssl/certs,"
                 " /etc/pki/tls/certs, /etc/ssl/cert.pem and friends).  HTTPS requests through"
                 " the Wt backend will fail certificate verification; set a CA bundle path in"
                 " the LLM settings." << std::endl;
  }
}//apply_trust_settings(...)


/** Map a boost::asio error onto our coarse set.

 Wt hands over a single error_code with no indication of which phase failed, so this is
 necessarily less precise than the other two backends.  Anything TLS-related lands in
 `TlsHandshakeFailed` unless OpenSSL's own category identifies a certificate problem.
 */
NativeHttp::Error map_asio_error( const boost::system::error_code &ec )
{
  if( !ec )
    return NativeHttp::Error::Ok;

  if( ec == boost::asio::error::host_not_found
     || ec == boost::asio::error::host_not_found_try_again )
    return NativeHttp::Error::HostNotFound;

  if( ec == boost::asio::error::connection_refused
     || ec == boost::asio::error::network_unreachable
     || ec == boost::asio::error::host_unreachable
     || ec == boost::asio::error::connection_reset
     || ec == boost::asio::error::broken_pipe
     || ec == boost::asio::error::eof )
    return NativeHttp::Error::ConnectFailed;

  if( ec == boost::asio::error::timed_out )
    return NativeHttp::Error::IdleTimeout;

  // Wt sets this on its own response-size overflow (Client.C:335), not just on a real
  //  message-size error - and it must not be treated as a transient failure worth retrying.
  if( ec == boost::asio::error::message_size )
    return NativeHttp::Error::ResponseTooLarge;

  // Deliberately NOT mapped to Error::Cancelled.  Wt 3.7.1's handleReadHeaders has an inverted
  //  aborted_ test (Client.C:458-462) relative to every other handler, so *any* failure while
  //  reading response headers - a proxy resetting the connection, or the idle timeout, which
  //  sets timed_out and then shuts the socket down - arrives here as operation_aborted.  A real
  //  cancellation is recognized by our own stopRequested flag before this function is consulted,
  //  so anything reaching here is a genuine failure and must be reported as one.
  if( ec == boost::asio::error::operation_aborted )
    return NativeHttp::Error::IdleTimeout;

#ifdef WT_WITH_SSL
  if( ec.category() == boost::asio::error::get_ssl_category() )
  {
    // OpenSSL reports verification failures through this reason code; everything else in the
    //  SSL category is a handshake or protocol problem.
    const int reason = ERR_GET_REASON( ec.value() );
    if( (reason == SSL_R_CERTIFICATE_VERIFY_FAILED)
       || (reason == SSL_R_SSLV3_ALERT_BAD_CERTIFICATE)
       || (reason == SSL_R_TLSV1_ALERT_UNKNOWN_CA) )
      return NativeHttp::Error::TlsCertificateUntrusted;

    return NativeHttp::Error::TlsHandshakeFailed;
  }
#endif

  return NativeHttp::Error::Unknown;
}//map_asio_error(...)


/** The io_service every request runs on.

 Normally the server's, so the client shares the thread pool Wt already runs.  Wt::Http::Client
 also uses it to decide whether to marshal signals into a session; we do our own marshalling in
 RequestState::dispatch(), so this is only about which thread the socket work happens on.
 */
Wt::WIOService *request_io_service()
{
  Wt::WServer * const server = Wt::WServer::instance();
  return server ? &server->ioService() : nullptr;
}//request_io_service()


/** Runs `fcn` on an io_service thread, posted to the io_service itself rather than through
 `WIOService::post()`: that goes through one strand shared by the whole server, so it waits behind
 any long-running posted work (in any session).  Returns false if there is no server.
 */
bool post_to_io( std::function<void()> fcn )
{
  Wt::WIOService * const ios = request_io_service();
  if( !ios )
    return false;
  boost::asio::post( *ios, std::move(fcn) );
  return true;
}//post_to_io(...)

}//namespace


namespace NativeHttp
{

bool available()
{
#ifdef WT_WITH_SSL
  // Wt::Http::Client needs a running WServer for its io_service.  Without one there is nothing to
  //  drive the sockets, so report unavailable rather than hanging.
  return (request_io_service() != nullptr);
#else
  // Wt was built with ENABLE_SSL=OFF, so Client::request() rejects any https:// URL outright.
  //  Saying so here means the caller falls back to the browser instead of issuing a request that
  //  can only fail - and the message names the build flag that fixes it.
  return false;
#endif
}

Backend backend()
{
  return Backend::WtHttpClient;
}

const char *backendName()
{
#ifdef WT_WITH_SSL
  return "Wt::Http::Client";
#else
  return "Wt::Http::Client (built without SSL - rebuild the prefix with WT_ENABLE_SSL=ON)";
#endif
}

bool supportsProxyOverride()
{
  return false;   // Wt::Http::Client has no proxy support of any kind
}

bool supportsCaBundle()
{
  return true;
}


namespace Detail
{

void startBackendRequest( const std::shared_ptr<RequestState> &state )
{
  Wt::WIOService * const ioService = request_io_service();
  if( !ioService )
  {
    state->deliverComplete( Error::BackendUnavailable,
                            "No Wt server is running, so there is no io_service to issue the"
                            " request on" );
    return;
  }

  const Request &req = state->request;

  if( req.method != "POST" && req.method != "GET" )
  {
    state->deliverComplete( Error::OptionUnsupported,
                            "The Wt::Http::Client backend only supports GET and POST, not "
                            + req.method );
    return;
  }

  const std::shared_ptr<WtRequestData> data = std::make_shared<WtRequestData>();
  data->client = new Wt::Http::Client( *ioService );

  Wt::Http::Client * const client = data->client;

  // setTimeout() restarts per async operation, so it is an idle timeout, not an overall one.
  //  There is no Wt equivalent of totalTimeout; the transport-level guard here is the idle one,
  //  and the caller's conversation watchdog backstops the rest.
  if( req.idleTimeout.count() > 0 )
    client->setTimeout( req.idleTimeout );   //Wt4 takes a std::chrono duration (Wt3 took int seconds)

  // Wt does follow redirects when asked, but only for GET (Client.C:1079-1082), so say so rather
  //  than silently not following a redirected POST.
  if( req.followRedirects )
  {
    if( req.method == "GET" )
    {
      client->setFollowRedirect( true );
    }else
    {
      state->deliverComplete( Error::OptionUnsupported,
                              "This backend can only follow redirects for GET requests" );
      delete client;
      return;
    }
  }//if( redirects were asked for )

  // NOT zero.  Zero disables Wt's accumulation, and its "response complete" test is
  //  `body().size() >= contentLength_` - which then never becomes true for a Content-Length
  //  response, so the request would hang until the idle timeout on every success.
  client->setMaximumResponseSize( req.maxResponseSize > 0 ? req.maxResponseSize
                                                          : (64u * 1024u * 1024u) );

  apply_trust_settings( *client, req );

  // Signals must be connected before request() is called: Wt only wires headersReceived and
  //  bodyDataReceived through to its internals `if (...isConnected())` at request time.
  client->headersReceived().connect( [state]( const Wt::Http::Message &response ){
      HeaderList headers;
      const std::vector<Wt::Http::Message::Header> &hdrs = response.headers();
      for( const Wt::Http::Message::Header &h : hdrs )
        headers.emplace_back( h.name(), h.value() );
      state->deliverHeaders( response.status(), std::move(headers) );
    } );

  client->bodyDataReceived().connect( [state]( const std::string &text ){
      state->deliverChunk( text );
    } );

  client->done().connect( [state, data]( const Wt::AsioWrapper::error_code &err,
                                         const Wt::Http::Message & ){
      // Wt documents that done() may still arrive with a successful response after abort(), so
      //  trust our own flag over the error code for cancellation.
      if( state->stopRequested.load() )
      {
        state->deliverComplete( state->stoppedByHandler.load() ? Error::AbortedByHandler
                                                               : Error::Cancelled,
                                "Cancelled by the application" );
      }else if( err )
      {
        state->deliverComplete( map_asio_error(err), err.message() );
      }else
      {
        state->deliverComplete( Error::Ok, std::string() );
      }

      // Clear and delete the client from a separate io_service task, not here: this runs inside
      //  Client::emitDone(), which holds the client's implementation mutex, and
      //  cancelBackendRequest() holds `backendMutex` while calling abort(), which takes that
      //  mutex - so taking `backendMutex` here would deadlock against a concurrent cancel.  The
      //  lock order is always `backendMutex`, then the client's mutex.  (~Client() takes the
      //  client's mutex too, so it waits for this emission to finish unwinding.)
      const bool posted = post_to_io( [state, data](){
        Wt::Http::Client *toDelete = nullptr;
        {
          std::lock_guard<std::mutex> lock( state->backendMutex );
          toDelete = data->client;
          data->client = nullptr;
        }
        delete toDelete;   // Once cleared under the mutex, no cancel can be using it
      } );

      // With no server left there is nothing to post to; leaking the client at shutdown is safer
      //  than deleting it inside its own signal emission.
      (void)posted;
    } );

  {
    std::lock_guard<std::mutex> lock( state->backendMutex );
    state->backendData = data;
  }

  bool started = false;
  std::string failure;

  try
  {
    if( req.method == "GET" )
    {
      std::vector<Wt::Http::Message::Header> hdrs;
      for( const std::pair<std::string,std::string> &kv : req.headers )
        hdrs.emplace_back( kv.first, kv.second );
      started = client->get( req.url, hdrs );
    }else
    {
      Wt::Http::Message message;
      for( const std::pair<std::string,std::string> &kv : req.headers )
        message.setHeader( kv.first, kv.second );
      message.addBodyText( req.body );
      started = client->post( req.url, message );
    }
  }catch( const std::exception &e )
  {
    // load_verify_file() uses the throwing overload, so a bad CA bundle path lands here rather
    //  than returning false.
    failure = e.what();
  }

  if( !started )
  {
    // request() returns false for a malformed URL, an unsupported scheme (including https when
    //  Wt has no SSL), or another request already in flight - and emits no done() signal at all,
    //  so the completion has to be synthesized here.
    Error err = Error::Unknown;
    std::string detail = failure;

    if( detail.empty() )
    {
#ifndef WT_WITH_SSL
      err = Error::BackendUnavailable;
      detail = "Wt was built without SSL support, so it cannot make https requests.  Rebuild the"
               " dependency prefix with WT_ENABLE_SSL=ON.";
#else
      detail = "Could not start the request; the endpoint URL may be malformed";
#endif
    }

    {
      std::lock_guard<std::mutex> lock( state->backendMutex );
      data->client = nullptr;
    }
    delete client;   // A cancel takes the mutex to use the client, so it is now done with it

    state->deliverComplete( err, detail );
  }//if( the request never started )
}//startBackendRequest(...)


void cancelBackendRequest( const std::shared_ptr<RequestState> &state )
{
  std::shared_ptr<WtRequestData> data;

  {
    std::lock_guard<std::mutex> lock( state->backendMutex );
    data = std::static_pointer_cast<WtRequestData>( state->backendData );
  }

  if( !data )
    return;

  // Abort from an io_service thread, not the caller's.  `data` is captured so the client cannot be
  //  deleted from under the post.
  const std::shared_ptr<RequestState> keepAlive = state;
  post_to_io( [data, keepAlive](){
    // Hold the mutex across abort(): the client is only deleted after being cleared under it (see
    //  the done() handler), so it cannot be freed while we use it.  abort() takes the client's own
    //  mutex, which may briefly wait for a done() emission in progress - that handler does not
    //  take `backendMutex`, so this cannot deadlock.
    std::lock_guard<std::mutex> lock( keepAlive->backendMutex );
    if( data->client )
      data->client->abort();
  } );
}//cancelBackendRequest(...)

}//namespace Detail

}//namespace NativeHttp
