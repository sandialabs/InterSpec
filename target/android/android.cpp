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
#include <string>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <iostream>

#include <Wt/WString>
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
       : m_type( src ), m_origSrc( NULL )
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
    int overflow(int c)
    {
        if (c == traits_type::eof()) {
            *this->pptr() = traits_type::to_char_type(c);
            this->sbumpc();
        }
        return this->sync() ? traits_type::eof() : traits_type::not_eof(c);
    }//int overflow(int c)
    int sync()
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


extern "C" 
{

JNIEXPORT
  jint
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_startServingInterSpec
  (JNIEnv* env, jobject thiz, jstring jprocess_name, jstring juserdatadir, jstring jbasedir, jstring jxml_config_path )
  {
    std::string process_name = std::string(env->GetStringUTFChars(jprocess_name, 0));
    std::string userdatadir = std::string(env->GetStringUTFChars(juserdatadir, 0));
    std::string basedir = std::string(env->GetStringUTFChars(jbasedir, 0));
    std::string xml_config_path = std::string(env->GetStringUTFChars(jxml_config_path, 0));

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
    const std::string sessionToken = std::string(env->GetStringUTFChars(jsessionToken, 0));
    const std::string filepath = std::string(env->GetStringUTFChars(jfilepath, 0));
    
    __android_log_write( ANDROID_LOG_INFO, "killServer", 
                         (std::string("Will open file '") + filepath + " in session '" + sessionToken + "'").c_str() );

    // We need to turn the file path into JSON.
    return InterSpecServer::open_file_in_session( sessionToken.c_str(), ("[\"" + filepath + "\"]").c_str() );
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
    const std::string token = std::string(env->GetStringUTFChars(jtoken, 0));

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
    const std::string token = std::string(env->GetStringUTFChars(jtoken, 0));

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
  Java_gov_sandia_InterSpec_removeSessionToken( JNIEnv* env, jobject thiz, jstring jtoken )
  {
    const std::string token = std::string(env->GetStringUTFChars(jtoken, 0));

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
  Java_gov_sandia_InterSpec_setInitialFileToLoad( JNIEnv* env, jobject thiz, jstring jtoken, jstring jfilepath )
  {
    const std::string token = std::string(env->GetStringUTFChars(jtoken, 0));
    const std::string filepath = std::string(env->GetStringUTFChars(jfilepath, 0));
    
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
  jint
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_openfileininterppec
  (JNIEnv* env, jobject thiz, jstring path, jint type, jstring sessionid )
  {
    const std::string filepath = std::string(env->GetStringUTFChars(path, 0));
    const std::string id = std::string(env->GetStringUTFChars(sessionid, 0));
//   env->ReleaseStringUTFChars(...)
    __android_log_write( ANDROID_LOG_INFO, "openfileininterppec",
                         (std::string("Will try to open following file in InterSpec '") + filepath + "' in session " + id).c_str());
    
	 //
    InterSpecApp *app = InterSpecApp::instanceFromExtenalToken( id );
    if( !app )
    {
      __android_log_write( ANDROID_LOG_INFO, "openfileininterppec",
                          (std::string("Couldnt get app instance from id '") + id + "'").c_str() );
      
      std::set<InterSpecApp *> instances = InterSpecApp::runningInstances();
      if( instances.empty() )
        return 1;
      app = *instances.begin();
    }
	
	  Wt::WApplication::UpdateLock lock( app );
    const bool opened = app->userOpenFromFileSystem( filepath );
	  app->triggerUpdate();

	  if( !opened )
	  {
      __android_log_write( ANDROID_LOG_INFO, "openfileininterppec",
                           (std::string("InterSpec couldnt open file '") + filepath + "'").c_str() );
	    return 2;
    }
	
    __android_log_write( ANDROID_LOG_INFO, "openfileininterppec",
                         (std::string("InterSpec opened file '") + filepath + "'").c_str() );

    return 0;
  }//Java_gov_sandia_InterSpec_openfileininterppec
  
  
void Java_gov_sandia_InterSpec_InterSpec_settmpdir( JNIEnv* env, jobject thiz, jstring tmpPath )
{
  const char *path = env->GetStringUTFChars( tmpPath, 0 );
  setenv("TMPDIR", path, 1);
  __android_log_write( ANDROID_LOG_DEBUG, "settmpdir",
                       (std::string("Set TMPDIR='") + path + "'").c_str() );
}

}  //extern "C"






int main( int argc, char **argv )
{ 
  assert( 0 );
}//int main( int argc, const char * argv[] )



