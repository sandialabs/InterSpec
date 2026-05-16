// libFuzzer driver for InterSpec::handleAppUrl - the top-level "interspec://"
// URL dispatcher.  A WTestEnvironment-backed headless InterSpec instance is
// created once in LLVMFuzzerInitialize, then each fuzz input is fed to
// handleAppUrl as a URL string.  This exercises every per-tool handleAppUrl
// (decay, simple-mda, dose, flux, gammaxs, drf, ...) through the same dispatch
// path that the EnterAppUrlWindow uses at runtime.
//
// State accumulation: each handleAppUrl call may open a window or change tools.
// After many iterations the app state grows.  Run with -rss_limit_mb=4096 (or
// similar) to bound memory; libFuzzer will report a leak if it crosses that.
//
// Process-exit handling: Wt 3.7.1's WTestEnvironment teardown is fragile - the
// WebController destroys its session and tries to lock a mutex that another
// static destructor has already destroyed, throwing std::system_error.  This
// would otherwise be reported as a fuzzer "crash" even though the iteration
// itself ran fine.  We install an atexit handler that calls _exit(0) to skip
// the broken teardown.  Real bugs during fuzz iterations are unaffected.
//
// See target/fuzzing/README.md for build & run instructions.

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>

#include <Wt/WApplication>
#include <Wt/Test/WTestEnvironment>

#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DecayDataBaseServer.h"


namespace
{
struct FuzzAppHolder
{
  // Lifetime mirrors target/testing/test_AnalystChecks.cpp's InterSpecTestFixture:
  // the InterSpecApp is intentionally raw-pointered and leaked at process exit
  // because deleting it after the WTestEnvironment is torn down has been
  // observed to cause issues.  For libFuzzer this is fine - one process, one
  // leak.
  std::unique_ptr<Wt::Test::WTestEnvironment> env;
  InterSpecApp *app = nullptr;
  std::unique_ptr<Wt::WApplication::UpdateLock> update_lock;
  InterSpec *interspec = nullptr;
};

std::unique_ptr<FuzzAppHolder> g_holder;


// Locate `data/sandia.decay.xml` by checking $INTERSPEC_DATA_DIR first, then
// walking upward from the cwd.  Mirrors the search the existing unit tests do.
std::string find_data_dir()
{
  if( const char *env = std::getenv("INTERSPEC_DATA_DIR") )
  {
    if( SpecUtils::is_file( SpecUtils::append_path(env, "sandia.decay.xml") ) )
      return env;
  }

  for( const char *d : {
        "data", "../data", "../../data", "../../../data", "../../../../data",
        "../../../../../data" } )
  {
    if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      return d;
  }
  return {};
}

}//anonymous namespace


extern "C" int LLVMFuzzerInitialize( int * /*argc*/, char *** /*argv*/ )
{
  // Skip Wt's WebController/session teardown at process exit - see file header.
  // _exit() bypasses static destructors entirely.  Registered before any of
  // Wt's static-init atexits would be effective, since atexit handlers run in
  // reverse order and LLVMFuzzerInitialize fires after static init.
  std::atexit( [](){ _exit(0); } );

  const std::string datadir = find_data_dir();
  if( datadir.empty() )
  {
    std::cerr << "fuzz_app_url: could not locate data/sandia.decay.xml.  "
                 "Set INTERSPEC_DATA_DIR or run from a directory near the "
                 "InterSpec checkout.\n";
    std::abort();
  }

  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path(datadir, "sandia.decay.xml") );
  if( !DecayDataBaseServer::database() )
  {
    std::cerr << "fuzz_app_url: failed to load decay database from " << datadir << "\n";
    std::abort();
  }

  InterSpec::setStaticDataDirectory( datadir );

  std::string wt_app_root = SpecUtils::lexically_normalize_path(
                              SpecUtils::append_path( datadir, "..") );

  g_holder = std::make_unique<FuzzAppHolder>();
  g_holder->env.reset( new Wt::Test::WTestEnvironment( "", "", Wt::Application ) );
  g_holder->env->setAppRoot( wt_app_root );

  g_holder->app = new InterSpecApp( *g_holder->env );
  g_holder->update_lock.reset( new Wt::WApplication::UpdateLock( g_holder->app ) );
  g_holder->interspec = g_holder->app->viewer();
  if( !g_holder->interspec )
  {
    std::cerr << "fuzz_app_url: InterSpecApp::viewer() returned null.\n";
    std::abort();
  }
  return 0;
}


extern "C" int LLVMFuzzerTestOneInput( const uint8_t *data, size_t size )
{
  if( size > 64 * 1024 )
    return 0;

  const std::string url( reinterpret_cast<const char *>(data), size );

  try
  {
    g_holder->interspec->handleAppUrl( url );
  }catch( const std::exception & )
  {
    // handleAppUrl throws runtime_error on bad input; that's expected.
  }

  // Process any deferred work the URL handler may have scheduled.
  Wt::WApplication::instance()->processEvents();
  return 0;
}
