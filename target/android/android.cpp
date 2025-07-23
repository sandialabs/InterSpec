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

#include <set>
#include <mutex>
#include <string>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <iostream>
#include <condition_variable>

#include <Wt/WString>
#include <Wt/WResource>
#include <Wt/WApplication>
#include <Wt/WEnvironment>

#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecServer.h"
#include "SpecUtils/SerialToDetectorModel.h"
#include "InterSpec/DataBaseVersionUpgrade.h"


#include <boost/filesystem.hpp>

#include "InterSpec/InterSpecApp.h"

#include <android/log.h>
#include <jni.h>

// I cant get getting Java methods here in C++ to work; it would be great to have
//  \c android_save_file_in_temp directly call into Java to save the temporary file,
//  this way we could bypass having to implement \c WebViewClient::shouldOverrideUrlLoading
#define TRY_TO_RESOLVE_JAVA_FUNCTIONS_FROM_CPP 0

namespace
{
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  std::mutex sm_interop_mutex;

  JavaVM* sm_jvm = nullptr;
  jobject sm_activity = 0; // GlobalRef - doesnt appear to be valid or work, after we set it...

#if( !TRY_TO_RESOLVE_JAVA_FUNCTIONS_FROM_CPP )
  std::condition_variable sm_saved_condition;
  bool sm_haveSavedFile;
  std::string sm_lastSavedLocation;
  std::string sm_lastSavedDisplayName;

  // Since I cant get resolving java function working from the c++, we'll set callbacks (using the
  // #setFileSaveCallback function) that we can call from C++ into Java.
  // Right now this call is basically a no-op, but in the future we should make it trigger the
  // "Save As" dialog, instead of using the WebViewClient::shouldOverrideUrlLoading Java function
  // to intercept downloads.
  jobject sm_java_file_save_cb_obj;
  jmethodID sm_java_file_save_cb_method_id;
#endif
}

namespace AndroidUtils
{
//adapted from http://stackoverflow.com/questions/8870174/is-stdcout-usable-in-android-ndk
class androidbuf: public std::streambuf
{
  //A utility to redirect cout/cerr to the Android logging system so it can be
  // seen using 'adb logcat'
public:
    enum Source { FromCout, FromCerr  };
    enum { bufsize = 128 }; // ... or some other suitable buffer size
    androidbuf( Source src )
       : m_type( src ), m_origSrc( nullptr )
    {
      this->setp(buffer, buffer + bufsize - 1);
      
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
      }//switch( src )
    }//androidbuf( Source src )
  
   ~androidbuf()
  {
    //cout/cerr must be given back there original rdbufs or else there can be a
    //  problems with freeing resources
    if( !m_origSrc )
      return;
    switch( m_type )
    {
      case FromCout: std::cout.rdbuf( m_origSrc ); break;
      case FromCerr: std::cerr.rdbuf( m_origSrc ); break;
    }
  }//~androidbuf()
  
private:
    int overflow(int c) override
    {
        if (c == traits_type::eof()) {
            *this->pptr() = traits_type::to_char_type(c);
            this->sbumpc();
        }
        return this->sync() ? traits_type::eof() : traits_type::not_eof(c);
    }//int overflow(int c)

    int sync() override
    {
        int rc = 0;
        if (this->pbase() != this->pptr())
        {
          rc = __android_log_write( ANDROID_LOG_INFO, m_source,
                            std::string(this->pbase(), this->pptr()).c_str() );
          this->setp(buffer, buffer + bufsize - 1);
        }
        return rc;
    }//int sync()
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
  }//for( loop over files in current working directory  )
}//void print_cwd_info()

}//namespace AndroidUtils


std::unique_ptr<AndroidUtils::androidbuf> g_stdbuf, g_errbuf;



void android_save_file_in_temp( bool success, std::string filepath, std::string displayName )
{
  // Using hacked saving to temporary file in Android, instead of via network download of file.

  __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp",
                       "Temp Filepath '%s' and displayName '%s'", filepath.c_str(), displayName.c_str() );
#if( TRY_TO_RESOLVE_JAVA_FUNCTIONS_FROM_CPP )
  if( !sm_jvm )
  {
    __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","The JVM has not been set." );
    return;
  }

  JNIEnv* env = nullptr;
  jint status = sm_jvm->AttachCurrentThread( &env, nullptr );

  if( status != JNI_OK )
  {
    __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","The JVM has not been set." );
    return;
  }

  if( !env )
  {
    __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Failed to get env." );
    sm_jvm->DetachCurrentThread();
    return;
  }

  jstring jFilepath = env->NewStringUTF( filepath.c_str() );
  jstring jDisplayName = env->NewStringUTF( displayName.c_str() );
  jclass cls = env->GetObjectClass( sm_activity );
  if( cls ) {
    __android_log_print(ANDROID_LOG_INFO, "save_file_in_temp", "Got class using sm_activity.");
  }

  jclass second_cls = env->FindClass("gov/sandia/InterSpec/InterSpec");
if( second_cls )
  __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Got class gov.sandia.InterSpec.InterSpec." );
else
  __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Did not get class gov.sandia.InterSpec.InterSpec." );

  auto exc = env->ExceptionOccurred();
  if (exc)
  {
    env->ExceptionDescribe();
    env->ExceptionClear();
  }

  //jclass third_cls = env->FindClass("gov/sandia/InterSpec");
  //if( third_cls )
  //  __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Got class InterSpec." );
  //else
  //  __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Did not get class InterSpec." );

  if( cls )
  {
    __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Got class." );
    jmethodID methodID = env->GetMethodID( cls, "saveFileFromTempLocation", "(Ljava/lang/String;Ljava/lang/String;)V" );

    exc = env->ExceptionOccurred();
    if (exc)
    {
      env->ExceptionDescribe();
      env->ExceptionClear();
    }

    if( methodID )
    {
      __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Got methodID." );
    env->CallObjectMethod( sm_activity, methodID, jFilepath, jDisplayName );
    }else
    {
      // AAAARRRRG:  I always end up here - I have no idea why...
      __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Did not get methodID." );
    }
  }else
  {
    __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Failed to get class." );
  }

  exc = env->ExceptionOccurred();
  if (exc)
  {
    env->ExceptionDescribe();
    env->ExceptionClear();
  }


  sm_jvm->DetachCurrentThread();
#else
  {
    std::unique_lock<std::mutex> scoped_lock( sm_interop_mutex );
    sm_haveSavedFile = success;
    sm_lastSavedLocation = filepath;
    sm_lastSavedDisplayName = displayName;
    sm_saved_condition.notify_all();
  }

  // TODO:
  JNIEnv* env = nullptr;
  jint status = sm_jvm->AttachCurrentThread( &env, nullptr );

  auto exc = env->ExceptionOccurred();
  if (exc)
  {
    env->ExceptionDescribe();
    env->ExceptionClear();
  }

  if( status != JNI_OK )
  {
    __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","The JVM has not been set." );
    return;
  }

  if( !env )
  {
    __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp","Failed to get env." );
    sm_jvm->DetachCurrentThread();
    return;
  }

  __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp", "Will call!" );

  env->CallVoidMethod(sm_java_file_save_cb_obj, sm_java_file_save_cb_method_id );

  exc = env->ExceptionOccurred();
  if (exc)
  {
    env->ExceptionDescribe();
    env->ExceptionClear();
  }

  __android_log_print( ANDROID_LOG_INFO, "save_file_in_temp", "Have called!" );

  exc = env->ExceptionOccurred();
  if (exc)
  {
    env->ExceptionDescribe();
    env->ExceptionClear();
  }

  sm_jvm->DetachCurrentThread();
#endif
}//android_save_file_in_temp(...)


void android_download_workaround( Wt::WResource *resource, std::string description )
{
  // Using hacked saving to temporary file in Android, instead of via network download of file.

  bool wrote_file = false;
  std::string suggestedName;
  const std::string tempFileName = SpecUtils::temp_file_name( description, SpecUtils::temp_dir() );

  {//Begin write temp file
    std::ofstream tmp( tempFileName.c_str(), std::ios_base::binary | std::ios_base::out );
    if( resource && tmp.is_open() )
    {
      resource->write( tmp );
      suggestedName = resource->suggestedFileName().toUTF8();
      wrote_file = true;

      __android_log_print( ANDROID_LOG_INFO, "android_workaround",
                           "Wrote photopeak ref CSV to: '%s'", tempFileName.c_str() );
    }//if( tmp.is_open() )
  }//android_download_workaround(...)

  android_save_file_in_temp( wrote_file, tempFileName, suggestedName );
}//android_download_workaround(...)


extern "C" 
{
JNIEXPORT void JNICALL Java_gov_sandia_InterSpec_InterSpec_initNative( JNIEnv* env, jobject thiz )
{
  std::unique_lock<std::mutex> lock(sm_interop_mutex);

  __android_log_print( ANDROID_LOG_INFO, "initNative", "Will init pointers to JVM and Activity" );
  env->GetJavaVM( &sm_jvm );
  sm_activity = env->NewGlobalRef( thiz ); // This doesnt really appear to work...
  __android_log_print( ANDROID_LOG_INFO, "initNative", "done initing" );
}

JNIEXPORT void JNICALL Java_gov_sandia_InterSpec_InterSpec_setFileSaveCallback
        (JNIEnv *env, jclass cls, jobject impl_obj)
{
  std::unique_lock<std::mutex> lock(sm_interop_mutex);

  jclass impl_cls;
  jmethodID impl_cb_mid;

  impl_cls = env->GetObjectClass(impl_obj);
  if (!impl_cls)
  {
    __android_log_print( ANDROID_LOG_INFO, "nativeCallJavaMethod", "!impl_cls" );
    return;
  }

  impl_cb_mid = env->GetMethodID(impl_cls, "callback", "()V");
  if (!impl_cb_mid)
  {
    __android_log_print( ANDROID_LOG_INFO, "nativeCallJavaMethod", "!impl_cb_mid" );
    return;
  }

  sm_java_file_save_cb_obj = env->NewGlobalRef( impl_obj );
  sm_java_file_save_cb_method_id = impl_cb_mid;
}


JNIEXPORT
  jint
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_startServingInterSpec
  (JNIEnv* env, jobject thiz, jstring jprocess_name, jstring juserdatadir, jstring jbasedir, jstring jxml_config_path )
  {
    std::string process_name = std::string(env->GetStringUTFChars(jprocess_name, nullptr));
    std::string userdatadir = std::string(env->GetStringUTFChars(juserdatadir, nullptr));
    std::string basedir = std::string(env->GetStringUTFChars(jbasedir, nullptr));
    std::string xml_config_path = std::string(env->GetStringUTFChars(jxml_config_path, nullptr));

    if( !g_stdbuf )
      g_stdbuf.reset( new AndroidUtils::androidbuf( AndroidUtils::androidbuf::FromCout ) );
    if( !g_errbuf )
      g_errbuf.reset( new AndroidUtils::androidbuf( AndroidUtils::androidbuf::FromCerr ) );

    // 20220424: TODO: test if we still need to set the CWD to `basedir`
    __android_log_print( ANDROID_LOG_INFO, "startServingInterSpec",
                          "Will chang working directory to '%s'", basedir.c_str() );

    boost::filesystem::current_path( basedir );

    AndroidUtils::print_cwd_info();
  
    std::cout << std::showbase << std::hex << "Running with Wt version "
              << WT_VERSION << std::dec << ", from executable compiled on "
              << __DATE__ << std::endl;
  
#if( PERFORM_DEVELOPER_CHECKS )
    std::cout << "Developer tests are being performed" << std::endl;
#endif

    std::cout << std::endl;
  
    //Make it so WString defaults to assuming std::string or char * are UTF8
    //  encoded, rather than the system encoding.
    Wt::WString::setDefaultEncoding( Wt::UTF8 );


    int serverPort = InterSpecServer::start_server( process_name.c_str(), userdatadir.c_str(), 
                                             basedir.c_str(), xml_config_path.c_str() );

    if( serverPort <= 0 )
    {
      //ToDo: should probably add some mechanism so interspec_start_server() can poplulate the error string.
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
  }//startServingInterSpec



  JNIEXPORT
  jint
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_openFile(JNIEnv* env, jobject thiz, jstring jsessionToken, jstring jfilepath) 
  {
    const std::string sessionToken = std::string(env->GetStringUTFChars(jsessionToken, nullptr));
    const std::string filepath = std::string(env->GetStringUTFChars(jfilepath, nullptr));
    
    __android_log_write( ANDROID_LOG_INFO, "killServer", 
                         (std::string("Will open file '") + filepath + " in session '" + sessionToken + "'").c_str() );

    // We need to turn the file path into JSON.
    return InterSpecServer::open_file_in_session( sessionToken.c_str(), ("[\"" + filepath + "\"]").c_str() );
  }

  JNIEXPORT
  jboolean
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_openAppUrl(JNIEnv* env, jobject thiz, jstring jsessionToken, jstring jurl)
  {
    const std::string sessionToken = std::string(env->GetStringUTFChars(jsessionToken, nullptr));
    const std::string url = std::string(env->GetStringUTFChars(jurl, nullptr));

    __android_log_write( ANDROID_LOG_INFO, "openAppUrl",
                         (std::string("Will open url '") + url + " in session '" + sessionToken + "'").c_str() );

    // We need to turn the file path into JSON.
    return InterSpecServer::pass_app_url_to_session( sessionToken.c_str(), url );
  }




  JNIEXPORT
  jboolean
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_killServer(JNIEnv* env, jobject thiz) 
  {
    InterSpecServer::killServer();
    __android_log_write( ANDROID_LOG_INFO, "killServer", "Killed server"  );
    return true;
  }

  
  JNIEXPORT
  jboolean
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_setTempDir( JNIEnv* env, jobject thiz, jstring jtmpdir ) 
  {
    const std::string tmpdir = std::string(env->GetStringUTFChars(jtmpdir, 0));
    
    setenv("TMPDIR", tmpdir.c_str(), 1);
    setenv("WT_TMP_DIR", tmpdir.c_str(), 1);

    __android_log_write( ANDROID_LOG_INFO, "setTempDir",
                          ("Set tmp directory to " + tmpdir).c_str()  );

    return true;
  }


  JNIEXPORT
  jboolean
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_setRequireSessionToken(JNIEnv* env, jobject thiz, jboolean jrequire )
  {  
    const bool require_token = jrequire;
  
    InterSpecServer::set_require_tokened_sessions( require_token );
  
    __android_log_write( ANDROID_LOG_INFO, "setRequireSessionToken",
                          ("Set require session token to " + std::to_string(require_token)).c_str()  );

    return true;
  }


  JNIEXPORT
  jboolean
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_addPrimarySessionToken( JNIEnv* env, jobject thiz, jstring jtoken )
  {
    const std::string token = std::string(env->GetStringUTFChars(jtoken, nullptr));

    const int status = InterSpecServer::add_allowed_session_token( token.c_str(), InterSpecServer::SessionType::PrimaryAppInstance );

    if( status != 0 )
    {
      __android_log_write( ANDROID_LOG_ERROR, "addPrimarySessionToken",
                          (std::string("Error adding primary session token='") + token + "'; code " + std::to_string(status)).c_str()  );
    }else
    {
      __android_log_write( ANDROID_LOG_INFO, "addPrimarySessionToken",
                          (std::string("Added session token='") + token + "'").c_str()  );
    }

    return (status == 0);
  }


  JNIEXPORT
  jboolean
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_addExternalSessionToken( JNIEnv* env, jobject thiz, jstring jtoken )
  {
    const std::string token = std::string(env->GetStringUTFChars(jtoken, nullptr));

    const int status = InterSpecServer::add_allowed_session_token( token.c_str(), InterSpecServer::SessionType::ExternalBrowserInstance );

    if( status != 0 )
    {
      __android_log_write( ANDROID_LOG_ERROR, "addExternalSessionToken",
                          (std::string("Error adding external session token='") + token + "'; code " + std::to_string(status)).c_str()  );
    }else
    {
      __android_log_write( ANDROID_LOG_INFO, "addExternalSessionToken",
                          (std::string("Added session token='") + token + "'").c_str()  );
    }

    return (status == 0);
  }



  JNIEXPORT
  jint
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_removeSessionToken( JNIEnv* env, jobject thiz, jstring jtoken )
  {
    const std::string token = std::string(env->GetStringUTFChars(jtoken, nullptr));

    const int removed = InterSpecServer::remove_allowed_session_token( token.c_str() );

    if( removed != 0 )
    {
      __android_log_write( ANDROID_LOG_ERROR, "removeSessionToken",
                          (std::string("Error removing session token='") + token + "'; code " + std::to_string(removed)).c_str()  );
    }else
    {
      __android_log_write( ANDROID_LOG_INFO, "removeSessionToken",
                          (std::string("Removed session token='") + token + "'").c_str()  );
    }

    return removed;
  }

  
  JNIEXPORT
  jint
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_setInitialFileToLoad( JNIEnv* env, jobject thiz, jstring jtoken, jstring jfilepath )
  {
    const std::string token = std::string(env->GetStringUTFChars(jtoken, nullptr));
    const std::string filepath = std::string(env->GetStringUTFChars(jfilepath, nullptr));
    
    try
    {
      InterSpecServer::set_file_to_open_on_load( token.c_str(), filepath );
      __android_log_write( ANDROID_LOG_INFO, "interspec_set_initial_file_to_open",
                           (std::string("Set initial file to open at filepath='") + filepath + "'").c_str()  );
    }catch( std::exception &e )
    {
      std::cerr << "interspec_set_initial_file_to_open: " << e.what() << std::endl;
      __android_log_write( ANDROID_LOG_ERROR, "interspec_set_initial_file_to_open",
                          (std::string("Error adding filepath='") + filepath + "'; msg: " + std::string(e.what()) ).c_str()  );
      return -1;
    }

    return 0;
  }//setInitialFileToLoad


JNIEXPORT
jobjectArray
JNICALL
Java_gov_sandia_InterSpec_InterSpec_mostRecentSaveLocation( JNIEnv* env, jobject thiz )
{
  __android_log_write( ANDROID_LOG_DEBUG, "mostRecentSaveLocation", "In mostRecentSaveLocation");

  // Using hacked saving to temporary file in Android, instead of via network download of file.
  bool have_saved = false;
  std::string location, display;


  {
    std::unique_lock<std::mutex> lock(sm_interop_mutex);

    // Wait for up to 1.5 second to get the location; an arbirary amount of time, but should be good enough to write most files to disk
    if( sm_lastSavedLocation.empty() )
    {
      __android_log_write( ANDROID_LOG_DEBUG, "mostRecentSaveLocation", "sm_lastSavedLocation was empty.");

      // auto pred = [&](){return !sm_lastSavedLocation.empty();};
      sm_saved_condition.wait_for(lock, std::chrono::milliseconds(1500) );
    }

    have_saved = sm_haveSavedFile;
    location = sm_lastSavedLocation;
    display = sm_lastSavedDisplayName;

    __android_log_write( ANDROID_LOG_DEBUG, "mostRecentSaveLocation", ("Got sm_lastSavedLocation=" + location).c_str() );

    sm_haveSavedFile = false;
    sm_lastSavedLocation = sm_lastSavedDisplayName = "";
  }

  std::cout << "mostRecentSaveLocation: have_saved=" << have_saved << ", location='" << location << "', display='" << display << "'" << std::endl;
  __android_log_write( ANDROID_LOG_DEBUG, "mostRecentSaveLocation", ("have_saved=" + std::string(have_saved ? "true" : "false") + ", location=" + location + ", display=" + display).c_str() );

  jobjectArray ret;
  if( have_saved && !location.empty() )
  {
    ret = (jobjectArray)env->NewObjectArray(2,env->FindClass("java/lang/String"),env->NewStringUTF(""));
    env->SetObjectArrayElement( ret,0, env->NewStringUTF( location.c_str() ) );
    env->SetObjectArrayElement( ret, 1 , env->NewStringUTF( display.c_str() ) );
  }else
  {
    ret = (jobjectArray)env->NewObjectArray(0,env->FindClass("java/lang/String"),env->NewStringUTF(""));
  }

  return ret;
}


}  //extern "C"






int main( int argc, char **argv )
{ 
  assert( 0 );
}//int main( int argc, const char * argv[] )



