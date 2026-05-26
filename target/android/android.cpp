//  Created by wcjohns on 20110322.
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
#include <fstream>
#include <iostream>

#include <boost/filesystem.hpp>

#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecServer.h"

#include <android/log.h>
#include <jni.h>


namespace
{
  // JVM + Java save-file callback registered via setFileSaveCallback().
  // Initialized once at app startup; subsequently read by Wt worker threads
  // that need to push a finished save back to the Android UI.
  std::mutex sm_jni_mutex;
  JavaVM* sm_jvm = nullptr;
  jobject sm_java_file_save_cb_obj = nullptr;     // GlobalRef
  jmethodID sm_java_file_save_cb_method_id = nullptr;


  // RAII wrapper around env->GetStringUTFChars / ReleaseStringUTFChars so that
  // every JNI string parameter is released regardless of exit path.
  class JStringUtf
  {
  public:
    JStringUtf( JNIEnv *env, jstring s )
      : m_env( env ), m_jstr( s ),
        m_chars( (env && s) ? env->GetStringUTFChars(s, nullptr) : nullptr )
    {}

    ~JStringUtf()
    {
      if( m_chars )
        m_env->ReleaseStringUTFChars( m_jstr, m_chars );
    }

    JStringUtf( const JStringUtf & ) = delete;
    JStringUtf &operator=( const JStringUtf & ) = delete;

    std::string str() const { return m_chars ? std::string(m_chars) : std::string{}; }

  private:
    JNIEnv * const m_env;
    const jstring m_jstr;
    const char * const m_chars;
  };//class JStringUtf


  // Attach the calling thread to the JVM and invoke the registered Java
  // "save file" callback with (tempPath, displayName).  Used by the native
  // file-save handler the lambda below registers with InterSpecApp.
  void invoke_java_save_callback( const std::string &tempFilePath,
                                  const std::string &displayName )
  {
    JavaVM *jvm = nullptr;
    jobject cb_obj = nullptr;
    jmethodID cb_mid = nullptr;
    {
      std::lock_guard<std::mutex> lock( sm_jni_mutex );
      jvm = sm_jvm;
      cb_obj = sm_java_file_save_cb_obj;
      cb_mid = sm_java_file_save_cb_method_id;
    }

    if( !jvm || !cb_obj || !cb_mid )
    {
      __android_log_print( ANDROID_LOG_ERROR, "save_file_callback",
                           "No Java save-file callback registered; dropping save of '%s'",
                           displayName.c_str() );
      return;
    }

    JNIEnv *env = nullptr;
    const jint status = jvm->AttachCurrentThread( &env, nullptr );
    if( status != JNI_OK || !env )
    {
      __android_log_print( ANDROID_LOG_ERROR, "save_file_callback",
                           "AttachCurrentThread failed (status=%d)", status );
      return;
    }

    const jstring jPath = env->NewStringUTF( tempFilePath.c_str() );
    const jstring jName = env->NewStringUTF( displayName.c_str() );
    env->CallVoidMethod( cb_obj, cb_mid, jPath, jName );

    if( env->ExceptionOccurred() )
    {
      env->ExceptionDescribe();
      env->ExceptionClear();
    }

    env->DeleteLocalRef( jPath );
    env->DeleteLocalRef( jName );

    jvm->DetachCurrentThread();
  }//invoke_java_save_callback
}//namespace


namespace AndroidUtils
{
  // Redirect cout/cerr to the Android logging system so it can be seen via
  // `adb logcat`.  Adapted from
  // http://stackoverflow.com/questions/8870174/is-stdcout-usable-in-android-ndk
  class androidbuf : public std::streambuf
  {
  public:
    enum Source { FromCout, FromCerr };
    enum { bufsize = 128 };

    androidbuf( Source src )
      : m_type( src ), m_origSrc( nullptr )
    {
      this->setp( buffer, buffer + bufsize - 1 );

      switch( m_type )
      {
        case FromCout:
          m_source = "cout";
          m_origSrc = std::cout.rdbuf( this );
          break;

        case FromCerr:
          m_source = "cerr";
          m_origSrc = std::cerr.rdbuf( this );
          break;
      }
    }

    ~androidbuf()
    {
      if( !m_origSrc )
        return;
      switch( m_type )
      {
        case FromCout: std::cout.rdbuf( m_origSrc ); break;
        case FromCerr: std::cerr.rdbuf( m_origSrc ); break;
      }
    }

  private:
    int overflow( int c ) override
    {
      if( c == traits_type::eof() )
      {
        *this->pptr() = traits_type::to_char_type( c );
        this->sbumpc();
      }
      return this->sync() ? traits_type::eof() : traits_type::not_eof( c );
    }

    int sync() override
    {
      int rc = 0;
      if( this->pbase() != this->pptr() )
      {
        rc = __android_log_write( ANDROID_LOG_INFO, m_source,
                                  std::string(this->pbase(), this->pptr()).c_str() );
        this->setp( buffer, buffer + bufsize - 1 );
      }
      return rc;
    }

    char buffer[bufsize];
    const char *m_source;
    const Source m_type;
    std::streambuf *m_origSrc;
  };//class androidbuf


  void print_cwd_info()
  {
    boost::filesystem::path cwd = boost::filesystem::current_path();
    __android_log_write( ANDROID_LOG_DEBUG, "print_cwd_info",
                         (std::string("Current path: '") + cwd.string<std::string>() + "' and contains:").c_str() );

    for( boost::filesystem::directory_iterator iter(cwd); iter != boost::filesystem::directory_iterator(); iter++ )
    {
      __android_log_write( ANDROID_LOG_DEBUG, "print_cwd_info",
                          (std::string("\t '") + iter->path().string<std::string>()
                          + (is_directory(iter->path()) ? "' dir " : "' file ")).c_str() );
    }
  }
}//namespace AndroidUtils


std::unique_ptr<AndroidUtils::androidbuf> g_stdbuf, g_errbuf;


extern "C"
{
JNIEXPORT void JNICALL Java_gov_sandia_InterSpec_InterSpec_initNative
  ( JNIEnv* env, jclass /*cls*/ )
{
  std::lock_guard<std::mutex> lock( sm_jni_mutex );
  env->GetJavaVM( &sm_jvm );
  __android_log_print( ANDROID_LOG_INFO, "initNative", "Cached JavaVM pointer" );
}


JNIEXPORT void JNICALL Java_gov_sandia_InterSpec_InterSpec_setFileSaveCallback
  ( JNIEnv *env, jclass /*cls*/, jobject impl_obj )
{
  std::lock_guard<std::mutex> lock( sm_jni_mutex );

  if( sm_java_file_save_cb_obj )
  {
    env->DeleteGlobalRef( sm_java_file_save_cb_obj );
    sm_java_file_save_cb_obj = nullptr;
    sm_java_file_save_cb_method_id = nullptr;
  }

  if( !impl_obj )
    return;

  jclass impl_cls = env->GetObjectClass( impl_obj );
  if( !impl_cls )
  {
    __android_log_write( ANDROID_LOG_ERROR, "setFileSaveCallback", "GetObjectClass returned null" );
    return;
  }

  // Java side: void callback(String tempFilePath, String displayName)
  jmethodID mid = env->GetMethodID( impl_cls, "callback",
                                    "(Ljava/lang/String;Ljava/lang/String;)V" );
  if( !mid )
  {
    __android_log_write( ANDROID_LOG_ERROR, "setFileSaveCallback",
                         "GetMethodID for callback(String,String) failed" );
    env->ExceptionClear();
    return;
  }

  sm_java_file_save_cb_obj = env->NewGlobalRef( impl_obj );
  sm_java_file_save_cb_method_id = mid;
}


JNIEXPORT jint JNICALL Java_gov_sandia_InterSpec_InterSpec_startServingInterSpec
  ( JNIEnv* env, jclass /*cls*/,
    jstring jprocess_name, jstring juserdatadir,
    jstring jbasedir, jstring jxml_config_path )
{
  const std::string process_name   = JStringUtf(env, jprocess_name  ).str();
  const std::string userdatadir    = JStringUtf(env, juserdatadir   ).str();
  const std::string basedir        = JStringUtf(env, jbasedir       ).str();
  const std::string xml_config_path= JStringUtf(env, jxml_config_path).str();

  if( !g_stdbuf )
    g_stdbuf.reset( new AndroidUtils::androidbuf( AndroidUtils::androidbuf::FromCout ) );
  if( !g_errbuf )
    g_errbuf.reset( new AndroidUtils::androidbuf( AndroidUtils::androidbuf::FromCerr ) );

  // 20220424: TODO: test whether this CWD change is still required.
  __android_log_print( ANDROID_LOG_INFO, "startServingInterSpec",
                       "Will change working directory to '%s'", basedir.c_str() );
  boost::filesystem::current_path( basedir );
  AndroidUtils::print_cwd_info();

  std::cout << std::showbase << std::hex << "Running with Wt version "
            << WT_VERSION << std::dec << ", from executable compiled on "
            << __DATE__ << std::endl;

#if( PERFORM_DEVELOPER_CHECKS )
  std::cout << "Developer tests are being performed" << std::endl;
#endif

  std::cout << std::endl;

  // Wt 4: WString is always UTF-8, setDefaultEncoding removed

  // Native file-save handler: Wt worker threads call this with the decoded
  // bytes of a download.  We write the bytes to a temp file, then JNI-callback
  // into Java with (tempPath, displayName) so the UI thread can present the
  // system Save-As dialog.
  InterSpecApp::setNativeFileSaveHandler(
    []( std::string data, std::string suggested_name ) {
      const std::string tempFileName
        = SpecUtils::temp_file_name( suggested_name, SpecUtils::temp_dir() );

      bool wrote_file = false;
      {
        std::ofstream tmp( tempFileName.c_str(),
                           std::ios_base::binary | std::ios_base::out );
        if( tmp.is_open() )
        {
          tmp.write( data.data(), data.size() );
          wrote_file = tmp.good();
        }
      }

      __android_log_print( ANDROID_LOG_INFO, "nativeFileSave",
                           "Wrote download to '%s' (success=%d, %zu bytes)",
                           tempFileName.c_str(), (int)wrote_file, data.size() );

      if( !wrote_file )
        return;

      invoke_java_save_callback( tempFileName, suggested_name );
    } );

  // Port 0 → let the OS choose an ephemeral port; the chosen port is returned.
  const int serverPort = InterSpecServer::start_server(
                              process_name.c_str(), userdatadir.c_str(),
                              basedir.c_str(), xml_config_path.c_str(),
                              0 /* server_port_num */ );

  if( serverPort <= 0 )
  {
    std::string errmsg;
    switch( serverPort )
    {
      case -1: errmsg = "Failed to create user data directory '" + userdatadir + "'"; break;
      case -2: errmsg = "InterSpecUserData.db: error in DataBaseUtils::setPreferenceDatabaseFile() or DbToFilesystemLink::setFileNumToFilePathDBNameBasePath()"; break;
      case -3: errmsg = "Error in set_detector_model_input_csv(), unexpected"; break;
      case -4: errmsg = "Error in InterSpec::setWritableDataDirectory('" + userdatadir + "')"; break;
      case -5: errmsg = "Error in DataBaseVersionUpgrade::checkAndUpgradeVersion()"; break;
      case -6: errmsg = "Error in ResourceUpdate::setupGlobalPrefsFromDb()"; break;
      case -7: errmsg = "Error in InterSpec::setStaticDataDirectory('" + basedir + "/data')"; break;
      case -8: errmsg = "Caught exception trying to start InterSpec server."; break;
      default: errmsg = "Unrecognized error code."; break;
    }

    __android_log_write( ANDROID_LOG_ERROR, "startServingInterSpec",
                         ("Failed to start server: '" + errmsg).c_str() );
  }

  return serverPort;
}


JNIEXPORT jint JNICALL Java_gov_sandia_InterSpec_InterSpec_openFile
  ( JNIEnv* env, jclass /*cls*/, jstring jsessionToken, jstring jfilepath )
{
  const std::string sessionToken = JStringUtf(env, jsessionToken).str();
  const std::string filepath     = JStringUtf(env, jfilepath    ).str();

  __android_log_write( ANDROID_LOG_INFO, "openFile",
                       ("Will open file '" + filepath + "' in session '" + sessionToken + "'").c_str() );

  return InterSpecServer::open_file_in_session(
            sessionToken.c_str(), ("[\"" + filepath + "\"]").c_str() );
}


JNIEXPORT jboolean JNICALL Java_gov_sandia_InterSpec_InterSpec_openAppUrl
  ( JNIEnv* env, jclass /*cls*/, jstring jsessionToken, jstring jurl )
{
  const std::string sessionToken = JStringUtf(env, jsessionToken).str();
  const std::string url          = JStringUtf(env, jurl         ).str();

  __android_log_write( ANDROID_LOG_INFO, "openAppUrl",
                       ("Will open url '" + url + "' in session '" + sessionToken + "'").c_str() );

  return InterSpecServer::pass_app_url_to_session( sessionToken.c_str(), url );
}


// Called by the Android WebView intercept (WebChromeClient.onCreateWindow or
// DownloadListener) when the user clicks a download link.  Fetches the URL
// in-process via Wt::Http::Client and pipes the bytes through the registered
// native file-save handler (the lambda set above in startServingInterSpec).
JNIEXPORT void JNICALL Java_gov_sandia_InterSpec_InterSpec_nativeDownloadUrl
  ( JNIEnv* env, jclass /*cls*/, jstring jurl )
{
  const std::string url = JStringUtf(env, jurl).str();
  if( url.empty() )
    return;

  __android_log_write( ANDROID_LOG_INFO, "nativeDownloadUrl",
                       ("Will download url='" + url + "'").c_str() );

  InterSpecServer::download_to_native_save( url );
}


JNIEXPORT jboolean JNICALL Java_gov_sandia_InterSpec_InterSpec_killServer
  ( JNIEnv* /*env*/, jclass /*cls*/ )
{
  InterSpecServer::killServer();
  __android_log_write( ANDROID_LOG_INFO, "killServer", "Killed server" );
  return true;
}


JNIEXPORT jboolean JNICALL Java_gov_sandia_InterSpec_InterSpec_setTempDir
  ( JNIEnv* env, jclass /*cls*/, jstring jtmpdir )
{
  const std::string tmpdir = JStringUtf(env, jtmpdir).str();

  setenv( "TMPDIR",    tmpdir.c_str(), 1 );
  setenv( "WT_TMP_DIR", tmpdir.c_str(), 1 );

  __android_log_write( ANDROID_LOG_INFO, "setTempDir",
                       ("Set tmp directory to " + tmpdir).c_str() );
  return true;
}


JNIEXPORT jboolean JNICALL Java_gov_sandia_InterSpec_InterSpec_setRequireSessionToken
  ( JNIEnv* /*env*/, jclass /*cls*/, jboolean jrequire )
{
  const bool require_token = jrequire;
  InterSpecServer::set_require_tokened_sessions( require_token );

  __android_log_write( ANDROID_LOG_INFO, "setRequireSessionToken",
                       ("Set require session token to " + std::to_string(require_token)).c_str() );
  return true;
}


JNIEXPORT jboolean JNICALL Java_gov_sandia_InterSpec_InterSpec_addPrimarySessionToken
  ( JNIEnv* env, jclass /*cls*/, jstring jtoken )
{
  const std::string token = JStringUtf(env, jtoken).str();

  const int status = InterSpecServer::add_allowed_session_token(
                          token.c_str(),
                          InterSpecServer::SessionType::PrimaryAppInstance );

  if( status != 0 )
    __android_log_write( ANDROID_LOG_ERROR, "addPrimarySessionToken",
                         ("Error adding primary session token='" + token + "'; code "
                          + std::to_string(status)).c_str() );
  else
    __android_log_write( ANDROID_LOG_INFO, "addPrimarySessionToken",
                         ("Added session token='" + token + "'").c_str() );

  return (status == 0);
}


JNIEXPORT jboolean JNICALL Java_gov_sandia_InterSpec_InterSpec_addExternalSessionToken
  ( JNIEnv* env, jclass /*cls*/, jstring jtoken )
{
  const std::string token = JStringUtf(env, jtoken).str();

  const int status = InterSpecServer::add_allowed_session_token(
                          token.c_str(),
                          InterSpecServer::SessionType::ExternalBrowserInstance );

  if( status != 0 )
    __android_log_write( ANDROID_LOG_ERROR, "addExternalSessionToken",
                         ("Error adding external session token='" + token + "'; code "
                          + std::to_string(status)).c_str() );
  else
    __android_log_write( ANDROID_LOG_INFO, "addExternalSessionToken",
                         ("Added session token='" + token + "'").c_str() );

  return (status == 0);
}


JNIEXPORT jint JNICALL Java_gov_sandia_InterSpec_InterSpec_removeSessionToken
  ( JNIEnv* env, jclass /*cls*/, jstring jtoken )
{
  const std::string token = JStringUtf(env, jtoken).str();

  const int removed = InterSpecServer::remove_allowed_session_token( token.c_str() );

  if( removed != 0 )
    __android_log_write( ANDROID_LOG_ERROR, "removeSessionToken",
                         ("Error removing session token='" + token + "'; code "
                          + std::to_string(removed)).c_str() );
  else
    __android_log_write( ANDROID_LOG_INFO, "removeSessionToken",
                         ("Removed session token='" + token + "'").c_str() );

  return removed;
}


JNIEXPORT jint JNICALL Java_gov_sandia_InterSpec_InterSpec_setInitialFileToLoad
  ( JNIEnv* env, jclass /*cls*/, jstring jtoken, jstring jfilepath )
{
  const std::string token    = JStringUtf(env, jtoken   ).str();
  const std::string filepath = JStringUtf(env, jfilepath).str();

  try
  {
    InterSpecServer::set_file_to_open_on_load( token.c_str(), filepath );
    __android_log_write( ANDROID_LOG_INFO, "setInitialFileToLoad",
                         ("Set initial file to open at filepath='" + filepath + "'").c_str() );
  }
  catch( std::exception &e )
  {
    std::cerr << "setInitialFileToLoad: " << e.what() << std::endl;
    __android_log_write( ANDROID_LOG_ERROR, "setInitialFileToLoad",
                         ("Error adding filepath='" + filepath + "'; msg: " + e.what()).c_str() );
    return -1;
  }

  return 0;
}

}  // extern "C"


// Wt 4's src/http/Android.C contains a `runMainAndCleanup` helper that calls
// `extern "C" int main(int, const char**)`.  That helper is only invoked
// through Wt's own JNI entry point (`Java_eu_webtoolkit_android_WtAndroid_startwt`)
// which InterSpec doesn't use -- we drive the embedded server via our own
// JNI bridge above -- but Android.C.o still gets pulled into libInterSpecAppLib.so
// because WServer::WServer references `preventRemoveOfSymbolsDuringLinking()`
// from the same TU.  So we have to supply a dummy `main` to satisfy the linker.
// It is unreachable at runtime; assert if anything ever calls it.
int main( int /*argc*/, const char ** /*argv*/ )
{
  assert( 0 );
  return -1;
}
