/* Security regression tests for URL validation and desktop session tokens. */

#include "InterSpec_config.h"

#define BOOST_TEST_MODULE test_SecurityHardening_suite
#include <boost/test/included/unit_test.hpp>

#include <string>

#include "InterSpec/InterSpecServer.h"
#if( USE_REMOTE_RID )
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/RemoteRid.h"
#endif

using namespace std;

#if( USE_REMOTE_RID )
BOOST_AUTO_TEST_CASE( remote_rid_url_compatibility_and_rejection )
{
  const char *valid_urls[] = {
    "http://localhost:8080",
    "https://127.0.0.1:8443/api/v1",
    "https://rid.example.test/path?drf=NaI%203x3#result",
    "https://user:password@rid.example.test/api"
  };
  for( const char *url : valid_urls )
    BOOST_CHECK_MESSAGE( RemoteRid::isValidRestUrl(url), url );

  const char *invalid_urls[] = {
    "file:///tmp/service",
    "ftp://rid.example.test/api",
    "https://rid.example.test/'-alert(1)-'",
    "https://rid.example.test/\\payload",
    "https://rid.example.test/path\nnext",
    "https:///missing-host"
  };
  for( const char *url : invalid_urls )
    BOOST_CHECK_MESSAGE( !RemoteRid::isValidRestUrl(url), url );
}
#endif


BOOST_AUTO_TEST_CASE( desktop_tokens_must_be_preregistered )
{
  InterSpecServer::set_require_tokened_sessions( true );

  const string unknown = "security-test-unknown-token";
  BOOST_CHECK_EQUAL( InterSpecServer::session_status(unknown.c_str()), 0 );
  BOOST_CHECK_EQUAL( InterSpecServer::set_session_loaded(unknown.c_str()), 2 );
  BOOST_CHECK_EQUAL( InterSpecServer::session_status(unknown.c_str()), 0 );

  const string primary = "security-test-primary-token";
  BOOST_CHECK_EQUAL(
    InterSpecServer::add_allowed_session_token(
      primary.c_str(), InterSpecServer::SessionType::PrimaryAppInstance),
    0 );
  const auto primary_type = InterSpecServer::session_type( primary.c_str() );
  BOOST_REQUIRE( primary_type.first );
  BOOST_CHECK( primary_type.second == InterSpecServer::SessionType::PrimaryAppInstance );
  BOOST_CHECK_EQUAL( InterSpecServer::set_session_loaded(primary.c_str()), 0 );
  BOOST_CHECK_EQUAL( InterSpecServer::set_session_loaded(primary.c_str()), 1 );

  const string browser = "security-test-browser-token";
  BOOST_CHECK_EQUAL(
    InterSpecServer::add_allowed_session_token(
      browser.c_str(), InterSpecServer::SessionType::ExternalBrowserInstance),
    0 );
  const auto browser_type = InterSpecServer::session_type( browser.c_str() );
  BOOST_REQUIRE( browser_type.first );
  BOOST_CHECK( browser_type.second == InterSpecServer::SessionType::ExternalBrowserInstance );
  BOOST_CHECK_EQUAL( InterSpecServer::set_session_loaded(browser.c_str()), 0 );
}
