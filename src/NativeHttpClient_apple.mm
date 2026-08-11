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

#import <Foundation/Foundation.h>
#import <Security/Security.h>
#import <CFNetwork/CFNetwork.h>

#include "InterSpec/NativeHttpClient.h"
#include "InterSpec/NativeHttpClientImpl.h"

static_assert( USE_NATIVE_HTTP_CLIENT,
              "This file should only be compiled when USE_NATIVE_HTTP_CLIENT is enabled" );

/** The macOS and iOS backend, built on NSURLSession.

 Identical for both platforms - there is deliberately no `#ifdef TARGET_OS_IPHONE` anywhere in
 this file.

 Why NSURLSession rather than a bundled TLS/HTTP library: it reads the system proxy settings and
 evaluates PAC scripts with no configuration from us, and it validates against
 Security.framework, which includes anchors installed by MDM and configuration profiles.  A
 bundled library sees none of that, which is exactly what breaks on a managed corporate Mac.

 Two implementation points matter:

 - We use the delegate API, not the `completionHandler:` convenience methods.  Those buffer the
   entire response before handing it over, so there would be no incremental delivery at all.

 - The delegate queue is an NSOperationQueue forced to be *serial*.  NSURLSession runs its own
   threads and does not need a `CFRunLoop` pumped by us, so it coexists with Wt's `io_service`
   without either interfering with the other - but a default-constructed NSOperationQueue is
   concurrent, and NSURLSession only guarantees callback ordering on a serial one.
 */

namespace
{
/** Owns the ObjC objects for one in-flight request.

 ARC permits `__strong` object references inside a C++ struct in Objective-C++, which is what
 lets this live in the `shared_ptr<void>` slot on `RequestState`.
 */
struct AppleRequestData
{
  __strong NSURLSession *session = nil;
  __strong NSURLSessionDataTask *task = nil;
};//struct AppleRequestData


std::string to_std_string( NSString * const str )
{
  if( !str )
    return std::string();
  const char * const utf8 = [str UTF8String];
  return utf8 ? std::string(utf8) : std::string();
}//to_std_string(...)


/** Map an NSError from NSURLSession onto our coarse error set.

 The point of the mapping is that the user is told something actionable: "the certificate was
 not trusted" leads somewhere, "error -1202" does not.
 */
NativeHttp::Error map_ns_error( NSError * const error )
{
  // Proxy failures arrive in CFNetwork's domain rather than NSURLErrorDomain, so they have to be
  //  matched separately - otherwise an unreachable or unauthenticated proxy is reported as an
  //  "unrecognized network error", which tells the one user who most needs a hint nothing at all.
  if( [[error domain] isEqualToString:(NSString *)kCFErrorDomainCFNetwork] )
  {
    switch( [error code] )
    {
      case kCFErrorHTTPProxyConnectionFailure:
      case kCFErrorHTTPSProxyConnectionFailure:
      case kCFStreamErrorHTTPSProxyFailureUnexpectedResponseToCONNECTMethod:
      case kCFErrorPACFileError:
        return NativeHttp::Error::ProxyUnreachable;

      case kCFErrorHTTPBadProxyCredentials:
      case kCFErrorPACFileAuth:
        return NativeHttp::Error::ProxyAuthRequired;

      default:
        break;
    }//switch( [error code] )
  }//if( a CFNetwork-domain error )

  if( ![[error domain] isEqualToString:NSURLErrorDomain] )
    return NativeHttp::Error::Unknown;

  switch( [error code] )
  {
    case NSURLErrorCancelled:
      return NativeHttp::Error::Cancelled;

    case NSURLErrorTimedOut:
      return NativeHttp::Error::Timeout;

    case NSURLErrorCannotFindHost:
    case NSURLErrorDNSLookupFailed:
      return NativeHttp::Error::HostNotFound;

    case NSURLErrorCannotConnectToHost:
    case NSURLErrorNetworkConnectionLost:
    case NSURLErrorNotConnectedToInternet:
    case NSURLErrorResourceUnavailable:
      return NativeHttp::Error::ConnectFailed;

    case NSURLErrorSecureConnectionFailed:
    case NSURLErrorClientCertificateRequired:
    case NSURLErrorClientCertificateRejected:
      return NativeHttp::Error::TlsHandshakeFailed;

    // The whole reason this backend exists on a managed network: a TLS-inspecting proxy
    //  re-signs certificates with a root that has to be in the system trust store.
    case NSURLErrorServerCertificateUntrusted:
    case NSURLErrorServerCertificateHasUnknownRoot:
    case NSURLErrorServerCertificateHasBadDate:
    case NSURLErrorServerCertificateNotYetValid:
      return NativeHttp::Error::TlsCertificateUntrusted;

    case NSURLErrorUserAuthenticationRequired:
      return NativeHttp::Error::ProxyAuthRequired;

    case NSURLErrorCannotLoadFromNetwork:
    case NSURLErrorBadServerResponse:
      return NativeHttp::Error::ConnectFailed;

    default:
      break;
  }//switch( [error code] )

  return NativeHttp::Error::Unknown;
}//map_ns_error(...)


/** Read a PEM bundle (or a directory of them) into SecCertificate objects.

 Used only for the custom-CA-bundle escape hatch.  Returns an empty array if nothing could be
 parsed, which the caller treats as a configuration error rather than silently continuing with
 the system anchors alone.
 */
NSArray *load_pem_certificates( const std::string &path )
{
  NSFileManager * const fm = [NSFileManager defaultManager];
  NSString * const nsPath = [NSString stringWithUTF8String:path.c_str()];

  BOOL isDir = NO;
  if( ![fm fileExistsAtPath:nsPath isDirectory:&isDir] )
    return @[];

  NSMutableArray * const files = [NSMutableArray array];
  if( isDir )
  {
    NSArray * const entries = [fm contentsOfDirectoryAtPath:nsPath error:nil];
    for( NSString *entry in entries )
      [files addObject:[nsPath stringByAppendingPathComponent:entry]];
  }else
  {
    [files addObject:nsPath];
  }

  NSMutableArray * const certs = [NSMutableArray array];

  for( NSString *file in files )
  {
    NSString * const contents = [NSString stringWithContentsOfFile:file
                                                         encoding:NSUTF8StringEncoding
                                                            error:nil];
    if( !contents )
      continue;

    // Pull out each -----BEGIN CERTIFICATE----- ... -----END CERTIFICATE----- block and turn
    //  its base64 payload into DER, which is what SecCertificateCreateWithData wants.
    NSScanner * const scanner = [NSScanner scannerWithString:contents];
    while( ![scanner isAtEnd] )
    {
      // NB: scanUpToString returns NO when the target is already AT the scan location - which is
      //  the normal case for a PEM file that begins with the BEGIN marker.  Treating that as "no
      //  more certificates" parses zero certificates out of a perfectly good bundle, so its
      //  result is deliberately ignored; only scanString failing means we are done.
      NSString *body = nil;
      [scanner scanUpToString:@"-----BEGIN CERTIFICATE-----" intoString:nil];
      if( ![scanner scanString:@"-----BEGIN CERTIFICATE-----" intoString:nil] )
        break;
      if( ![scanner scanUpToString:@"-----END CERTIFICATE-----" intoString:&body] )
        break;
      [scanner scanString:@"-----END CERTIFICATE-----" intoString:nil];

      if( !body )
        continue;

      NSData * const der = [[NSData alloc]
                            initWithBase64EncodedString:body
                                                options:NSDataBase64DecodingIgnoreUnknownCharacters];
      if( !der || ([der length] == 0) )
        continue;

      SecCertificateRef cert = SecCertificateCreateWithData( NULL, (__bridge CFDataRef)der );
      if( cert )
      {
        [certs addObject:(__bridge_transfer id)cert];
      }
    }//while( more certificates in this file )
  }//for( each candidate file )

  return certs;
}//load_pem_certificates(...)


/** Turn "http://user:pass@host:port" into the CFNetwork proxy dictionary NSURLSession wants. */
NSDictionary *make_proxy_dictionary( const std::string &proxyUrl )
{
  NSString * const nsUrl = [NSString stringWithUTF8String:proxyUrl.c_str()];
  NSURL * const url = [NSURL URLWithString:nsUrl];
  if( !url || ![url host] )
    return nil;

  NSNumber * const port = [url port] ? [url port] : @8080;

  // Route both http and https (i.e. CONNECT tunneling) through the same proxy; a corporate
  //  proxy that handles one invariably handles the other.
  return @{
    (NSString *)kCFNetworkProxiesHTTPEnable  : @YES,
    (NSString *)kCFNetworkProxiesHTTPProxy   : [url host],
    (NSString *)kCFNetworkProxiesHTTPPort    : port,
    @"HTTPSEnable"                           : @YES,
    @"HTTPSProxy"                            : [url host],
    @"HTTPSPort"                             : port,
  };
}//make_proxy_dictionary(...)

}//namespace


/** NSURLSession delegate for exactly one request. */
@interface InterSpecHttpDelegate : NSObject <NSURLSessionDataDelegate>
{
@public
  std::shared_ptr<NativeHttp::Detail::RequestState> _state;
}
@end


@implementation InterSpecHttpDelegate

- (void)URLSession:(NSURLSession *)session
          dataTask:(NSURLSessionDataTask *)dataTask
didReceiveResponse:(NSURLResponse *)response
 completionHandler:(void (^)(NSURLSessionResponseDisposition))completionHandler
{
  if( !_state )
  {
    completionHandler( NSURLSessionResponseCancel );
    return;
  }

  int status = 0;
  NativeHttp::HeaderList headers;

  if( [response isKindOfClass:[NSHTTPURLResponse class]] )
  {
    NSHTTPURLResponse * const http = (NSHTTPURLResponse *)response;
    status = static_cast<int>( [http statusCode] );

    NSDictionary * const fields = [http allHeaderFields];
    for( id key in fields )
    {
      if( [key isKindOfClass:[NSString class]] )
      {
        id value = [fields objectForKey:key];
        headers.emplace_back( to_std_string((NSString *)key),
                              to_std_string([value description]) );
      }
    }//for( each response header )
  }//if( an HTTP response )

  _state->deliverHeaders( status, std::move(headers) );

  completionHandler( _state->shouldStop() ? NSURLSessionResponseCancel
                                          : NSURLSessionResponseAllow );
}//didReceiveResponse


- (void)URLSession:(NSURLSession *)session
          dataTask:(NSURLSessionDataTask *)dataTask
    didReceiveData:(NSData *)data
{
  if( !_state )
    return;

  // NSData may be discontiguous, so walk its ranges rather than calling -bytes (which would
  //  flatten the whole thing into one allocation).
  __block bool keepGoing = true;
  [data enumerateByteRangesUsingBlock:^(const void *bytes, NSRange range, BOOL *stop){
    if( !keepGoing )
    {
      *stop = YES;
      return;
    }

    const std::string_view chunk( static_cast<const char *>(bytes), range.length );
    keepGoing = _state->deliverChunk( chunk );
    if( !keepGoing )
      *stop = YES;
  }];

  if( !keepGoing )
    [dataTask cancel];
}//didReceiveData


- (void)URLSession:(NSURLSession *)session
              task:(NSURLSessionTask *)task
didCompleteWithError:(NSError *)error
{
  const std::shared_ptr<NativeHttp::Detail::RequestState> state = _state;

  // There is a retain cycle here by construction: RequestState -> backendData -> NSURLSession
  //  -> (its strongly-held) delegate -> _state -> RequestState.  Two things break it, and both
  //  are deliberate: invalidating the session releases the delegate, and clearing backendData
  //  drops the C++ end of the cycle regardless of when NSURLSession gets around to that.
  //  Without the second, any future early-return added ahead of this point would permanently
  //  leak the session, the request body and the bearer token with it.
  [session finishTasksAndInvalidate];

  if( state )
  {
    std::lock_guard<std::mutex> lock( state->backendMutex );
    state->backendData.reset();
  }

  if( !state )
    return;

  if( !error )
  {
    state->deliverComplete( NativeHttp::Error::Ok, std::string() );
    return;
  }

  // A cancel we asked for is not a failure to report as one - distinguish the handler aborting
  //  from the user cancelling, both of which surface here as NSURLErrorCancelled.
  NativeHttp::Error mapped = map_ns_error( error );
  if( mapped == NativeHttp::Error::Cancelled )
  {
    if( state->stoppedByHandler.load() )
      mapped = NativeHttp::Error::AbortedByHandler;
  }

  // Carry the numeric code as well as the text: it is stable across locales, and it is what a
  //  search engine actually matches on when someone pastes an error into a bug report.
  std::string detail = to_std_string( [error localizedDescription] );
  detail += " [" + to_std_string([error domain]) + " "
            + std::to_string(static_cast<long long>([error code])) + "]";

  state->deliverComplete( mapped, detail );
}//didCompleteWithError


- (void)URLSession:(NSURLSession *)session
              task:(NSURLSessionTask *)task
willPerformHTTPRedirection:(NSHTTPURLResponse *)response
        newRequest:(NSURLRequest *)request
 completionHandler:(void (^)(NSURLRequest *))completionHandler
{
  // NSURLSession follows redirects by default; passing nil here declines, which is what we want
  //  unless the caller explicitly opted in (see Request::followRedirects).
  if( _state && _state->request.followRedirects )
    completionHandler( request );
  else
    completionHandler( nil );
}//willPerformHTTPRedirection


- (void)URLSession:(NSURLSession *)session
              task:(NSURLSessionTask *)task
didReceiveChallenge:(NSURLAuthenticationChallenge *)challenge
 completionHandler:(void (^)(NSURLSessionAuthChallengeDisposition, NSURLCredential *))completionHandler
{
  NSString * const method = [[challenge protectionSpace] authenticationMethod];

  if( ![method isEqualToString:NSURLAuthenticationMethodServerTrust] )
  {
    // Default handling is what lets CFNetwork attempt SPNEGO against the Heimdal credential cache
    //  on a domain-bound Mac, giving single-sign-on against an intranet proxy for free.  But it
    //  is only worth asking for where it can actually succeed silently.
    //
    // For Basic/Digest - a proxy wanting a username and password - default handling makes macOS
    //  put a system-modal "Proxy Authentication Required" dialog in front of the user.  That is
    //  the wrong thing to do for a background request the user did not initiate, and it cannot
    //  help anyway: InterSpec has no credentials to offer and deliberately does not collect any.
    //  Declining instead surfaces our own ProxyAuthRequired error, which says what to do.
    //
    // A non-zero previousFailureCount means the credentials we could offer have already been
    //  rejected once; retrying only produces the same prompt.
    const bool canSucceedSilently =
        ([method isEqualToString:NSURLAuthenticationMethodNegotiate]
         || [method isEqualToString:NSURLAuthenticationMethodNTLM]
         || [method isEqualToString:NSURLAuthenticationMethodClientCertificate]);

    if( canSucceedSilently && ([challenge previousFailureCount] == 0) )
      completionHandler( NSURLSessionAuthChallengePerformDefaultHandling, nil );
    else
      completionHandler( NSURLSessionAuthChallengeRejectProtectionSpace, nil );
    return;
  }

  SecTrustRef trust = [[challenge protectionSpace] serverTrust];
  if( !trust || !_state )
  {
    completionHandler( NSURLSessionAuthChallengePerformDefaultHandling, nil );
    return;
  }

  const NativeHttp::Request &req = _state->request;

  if( req.disableCertVerification )
  {
    // Explicitly opted into by the user, with a warning, for the case where nothing else works.
    completionHandler( NSURLSessionAuthChallengeUseCredential,
                       [NSURLCredential credentialForTrust:trust] );
    return;
  }

  if( req.caBundlePath.empty() )
  {
    completionHandler( NSURLSessionAuthChallengePerformDefaultHandling, nil );
    return;
  }

  NSArray * const anchors = load_pem_certificates( req.caBundlePath );
  if( [anchors count] == 0 )
  {
    completionHandler( NSURLSessionAuthChallengeCancelAuthenticationChallenge, nil );
    return;
  }

  SecTrustSetAnchorCertificates( trust, (__bridge CFArrayRef)anchors );
  // false => keep the system anchors as well, so adding a corporate root does not stop ordinary
  //  public certificates from validating.
  SecTrustSetAnchorCertificatesOnly( trust, false );

  CFErrorRef trustError = NULL;
  const bool trusted = SecTrustEvaluateWithError( trust, &trustError );
  if( trustError )
    CFRelease( trustError );

  if( trusted )
    completionHandler( NSURLSessionAuthChallengeUseCredential,
                       [NSURLCredential credentialForTrust:trust] );
  else
    completionHandler( NSURLSessionAuthChallengeCancelAuthenticationChallenge, nil );
}//didReceiveChallenge

@end


namespace NativeHttp
{

bool available()
{
  return true;
}

Backend backend()
{
  return Backend::NSUrlSession;
}

const char *backendName()
{
  return "NSURLSession";
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
  @autoreleasepool
  {
    const Request &req = state->request;

    NSString * const urlStr = [NSString stringWithUTF8String:req.url.c_str()];
    NSURL * const url = [NSURL URLWithString:urlStr];
    if( !url )
    {
      state->deliverComplete( Error::Unknown,
                              "Could not parse the endpoint URL; check the provider's"
                              " endpoint in the LLM settings" );
      return;
    }

    NSMutableURLRequest * const request = [NSMutableURLRequest requestWithURL:url];
    [request setHTTPMethod:[NSString stringWithUTF8String:req.method.c_str()]];

    for( const std::pair<std::string,std::string> &kv : req.headers )
    {
      [request setValue:[NSString stringWithUTF8String:kv.second.c_str()]
     forHTTPHeaderField:[NSString stringWithUTF8String:kv.first.c_str()]];
    }

    if( !req.body.empty() )
      [request setHTTPBody:[NSData dataWithBytes:req.body.data() length:req.body.size()]];

    // defaultSessionConfiguration is what picks up the system proxy settings and evaluates PAC.
    NSURLSessionConfiguration * const config =
        [NSURLSessionConfiguration defaultSessionConfiguration];

    // NSURLSession's two timeouts map exactly onto ours: timeoutIntervalForRequest restarts
    //  whenever new data arrives (an idle timeout), while timeoutIntervalForResource bounds the
    //  whole transfer.
    if( req.idleTimeout.count() > 0 )
      [config setTimeoutIntervalForRequest:static_cast<NSTimeInterval>(req.idleTimeout.count())];
    if( req.totalTimeout.count() > 0 )
      [config setTimeoutIntervalForResource:static_cast<NSTimeInterval>(req.totalTimeout.count())];

    // An LLM response is generated fresh every time; a cached one would be wrong.
    [config setRequestCachePolicy:NSURLRequestReloadIgnoringLocalCacheData];
    [config setURLCache:nil];

    if( !req.proxyUrl.empty() )
    {
      NSDictionary * const proxyDict = make_proxy_dictionary( req.proxyUrl );
      if( !proxyDict )
      {
        state->deliverComplete( Error::Unknown,
                                "Could not parse the proxy URL in the LLM network settings" );
        return;
      }
      [config setConnectionProxyDictionary:proxyDict];
    }//if( an explicit proxy was configured )

    InterSpecHttpDelegate * const delegate = [[InterSpecHttpDelegate alloc] init];
    delegate->_state = state;

    // The delegate queue MUST be serial.  A default-constructed NSOperationQueue is *concurrent*
    //  (maxConcurrentOperationCount defaults to NSOperationQueueDefaultMaxConcurrentOperationCount),
    //  and NSURLSession only guarantees callback ordering on a serial queue.  With a concurrent
    //  one, two didReceiveData: callbacks for the same task can run at once, so body chunks would
    //  be appended out of order - which corrupts exactly the large, multi-packet JSON responses
    //  this transport exists to carry.
    NSOperationQueue * const delegateQueue = [NSOperationQueue new];
    delegateQueue.maxConcurrentOperationCount = 1;

    NSURLSession * const session =
        [NSURLSession sessionWithConfiguration:config
                                      delegate:delegate
                                 delegateQueue:delegateQueue];

    NSURLSessionDataTask * const task = [session dataTaskWithRequest:request];

    {
      const std::shared_ptr<AppleRequestData> data = std::make_shared<AppleRequestData>();
      data->session = session;
      data->task = task;

      std::lock_guard<std::mutex> lock( state->backendMutex );
      state->backendData = data;
    }

    // Cancellation could have been requested between start() and here.
    if( state->shouldStop() )
    {
      [task cancel];
      [session finishTasksAndInvalidate];
      return;
    }

    [task resume];
  }//@autoreleasepool
}//startBackendRequest(...)


void cancelBackendRequest( const std::shared_ptr<RequestState> &state )
{
  std::shared_ptr<AppleRequestData> data;

  {
    std::lock_guard<std::mutex> lock( state->backendMutex );
    data = std::static_pointer_cast<AppleRequestData>( state->backendData );
  }

  if( !data )
    return;

  @autoreleasepool
  {
    // -cancel makes the task fail with NSURLErrorCancelled, which routes through
    //  didCompleteWithError like any other ending, so there is only one completion path.
    if( data->task )
      [data->task cancel];
  }
}//cancelBackendRequest(...)

}//namespace Detail

}//namespace NativeHttp
