#ifndef AndroidUtils_hpp
#define AndroidUtils_hpp
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>


#if(ALLOW_URL_TO_FILESYSTEM_MAP)
#include "InterSpec/DbToFilesystemLink.h"
#else 
#error "ALLOW_URL_TO_FILESYSTEM_MAP must be enabled for android" 
#endif

#include "InterSpec/InterSpecApp.h"


#ifdef ANDROID

#include <android/log.h>
#include <jni.h>

extern "C" 
{
  JNIEXPORT
  jint
  JNICALL
  Java_gov_sandia_InterSpec_InterSpec_addopenfiletodb
  (JNIEnv* env, jobject thiz, jstring path )
  {
     std::string filepath = std::string(env->GetStringUTFChars(path, 0));
//   env->ReleaseStringUTFChars(...)

    DbToFilesystemLink::FileIdToLocation requestinfo;
//    requestinfo.m_userName = "wcjohns";
    requestinfo.m_foregroundFilePath = filepath;
    
//    requestinfo.m_backgroundFilePath = "/path/to/background";
    const int entrynum = DbToFilesystemLink::addFileToOpenToDatabase( requestinfo );
    
    __android_log_write( ANDROID_LOG_DEBUG, "addopenfiletodb",
                         (std::string("Added filepath='") + filepath + "' to database as entrynum=" + std::to_string(entrynum)).c_str()  );
    
    return entrynum;
  }//Java_gov_sandia_InterSpec_InterSpec_addopenfiletodb
  
  
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
#endif  //#ifdef ANDROID

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
  using namespace std;
  using namespace boost;
  using namespace boost::filesystem;

  path cwd = current_path();
  __android_log_write( ANDROID_LOG_DEBUG, "print_cwd_info",
                       (std::string("Current path: '") + cwd.string<std::string>() + "' and contains:").c_str() );
  
  for( directory_iterator iter(cwd); iter != directory_iterator(); iter++ )
  {
    __android_log_write( ANDROID_LOG_DEBUG, "print_cwd_info",
                        (std::string("\t '") + iter->path().string<std::string>()
                        + (is_directory(iter->path()) ? "' dir " : "' file ")).c_str() );
  }//for( loop over files in current working directory  )
}//void print_cwd_info()

void set_anrdoid_cwd( int argc, char **argv )
{
  const std::string docroot_tag = "--docroot";
  for( int i = 0; i < argc; ++i )
  {
    if( (docroot_tag==argv[i]) && (i<(argc-1)) )
    {
      __android_log_print( ANDROID_LOG_INFO, "InterSpec",
                          "Changing working directory to %s", argv[i+1] );
      boost::filesystem::current_path( argv[i+1] );
      break;
    }//if( (docroot_tag==argv[i]) && (i<(argc-1)) )
  }//for( int i = 0; i < argc; ++i )

//  print_cwd_info();
}//void set_anrdoid_cwd( int argc, char **argv )

}//namespace AndroidUtils
#endif
