
/* cppsrc/main.cpp */
#include <napi.h>

#include <string>
#include <stdlib.h>

#include "target/electron/ElectronUtils.h"

//A good, simple, N-API tutorial is at:
//https://medium.com/@atulanand94/beginners-guide-to-writing-nodejs-addons-using-c-and-n-api-node-addon-api-9b3b718a9a7f



namespace functionexample
{
  int add(int a, int b)
  {
    return a + b;
  }

  std::string hello()
  {
    return "Hello World.js";
  }

  Napi::Number startServingInterSpec(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();

    if( info.Length() < 4 || !info[0].IsString() || !info[1].IsString() 
        || !info[2].IsString() || !info[3].IsString() ) 
    {
        Napi::TypeError::New(env, "startServingInterSpec: Expected 4 strings").ThrowAsJavaScriptException();
    } 

    std::string process_name = info[0].ToString();
    std::string userdatadir = info[1].ToString();
    std::string basedir = info[2].ToString();
    std::string xml_config_path = info[3].ToString();

    int serverPort = interspec_start_server( process_name.c_str(), userdatadir.c_str(), 
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
      Napi::TypeError::New(env, "startServingInterSpec: Received error code " + std::to_string(serverPort) + ": " + errmsg ).ThrowAsJavaScriptException();
    }

    return Napi::Number::New( env, serverPort );
  }

  Napi::Number openFile(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();

    if (info.Length() < 2 || !info[0].IsString() || !info[1].IsString()) {
        Napi::TypeError::New(env, "openFile: Expected two strings").ThrowAsJavaScriptException();
    } 

    std::string sessionToken = info[0].As<Napi::String>();
    std::string filepath = info[1].As<Napi::String>();
    
    int result = interspec_open_file( sessionToken.c_str(), filepath.c_str() );

    return Napi::Number::New( env, result );
  }

  Napi::Boolean killServer(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();

    interspec_kill_server();

    return Napi::Boolean::New( env, true );
  }

  
  Napi::Boolean setTempDir(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();

    if (info.Length() < 1 || !info[0].IsString() ) {
        Napi::TypeError::New(env, "setBaseDir: Expected one strings").ThrowAsJavaScriptException();
    } 

#ifdef WIN32
      std::wstring tmpdir = info[0].ToString()
      _wputenv_s( L"TMPDIR", tmpdir.c_str() );
      _wputenv_s( L"WT_TMP_DIR", tmpdir.c_str() );
#else
        std::string tmpdir = info[0].ToString();
        setenv("TMPDIR", tmpdir.c_str(), 1);
        setenv("WT_TMP_DIR", tmpdir.c_str(), 1);
#endif // WIN32


    return Napi::Boolean::New( env, true );
  }
  

  Napi::Boolean addSessionToken(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();

    if (info.Length() < 1 || !info[0].IsString() ) {
        Napi::TypeError::New(env, "addSessionToken: Expected one strings").ThrowAsJavaScriptException();
    } 

    std::string token = info[0].ToString();

    interspec_add_allowed_session_token( token.c_str() );

    return Napi::Boolean::New( env, true );
  }


  Napi::Number removeSessionToken(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();

    if (info.Length() < 1 || !info[0].IsString() ) {
        Napi::TypeError::New(env, "removeSessionToken: Expected one strings").ThrowAsJavaScriptException();
    } 

    std::string token = info[0].ToString();

    const int removed = interspec_remove_allowed_session_token( token.c_str() );

    return Napi::Number::New( env, removed );
  }


  Napi::Number AddWrapped(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();
    if (info.Length() < 2 || !info[0].IsNumber() || !info[1].IsNumber()) {
        Napi::TypeError::New(env, "Number expected").ThrowAsJavaScriptException();
    } 

    Napi::Number first = info[0].As<Napi::Number>();
    Napi::Number second = info[1].As<Napi::Number>();

    int returnValue = functionexample::add(first.Int32Value(), second.Int32Value());
    
    return Napi::Number::New(env, returnValue);
  }
}//namespace functionexample

Napi::Object InitAll(Napi::Env env, Napi::Object exports) {
  exports.Set( "startServingInterSpec", Napi::Function::New(env, functionexample::startServingInterSpec) );
  exports.Set( "killServer", Napi::Function::New(env, functionexample::killServer) );
  exports.Set( "openFile", Napi::Function::New(env, functionexample::openFile ));

  exports.Set( "setTempDir", Napi::Function::New(env, functionexample::setTempDir ));
  exports.Set( "addSessionToken", Napi::Function::New(env, functionexample::addSessionToken ));
  exports.Set( "removeSessionToken", Napi::Function::New(env, functionexample::removeSessionToken ));


  exports.Set("add", Napi::Function::New(env, functionexample::AddWrapped));
  return exports;
}


/**
* This code defines the entry-point for the Node addon, it tells Node where to go
* once the library has been loaded into active memory. The first argument must
* match the "target" in our *binding.gyp*. Using NODE_GYP_MODULE_NAME ensures
* that the argument will be correct, as long as the module is built with
* node-gyp (which is the usual way of building modules). The second argument
* points to the function to invoke. The function must not be namespaced.
*/
NODE_API_MODULE(testaddon, InitAll)
