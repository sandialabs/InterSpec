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

#include <napi.h>

#include <array>
#include <mutex>
#include <string>
#include <locale>
#include <codecvt>
#include <iostream>
#include <stdlib.h>
#include <condition_variable>

#include <boost/optional.hpp>

#include "target/electron/ElectronUtils.h"

//A good, simple, N-API tutorial is at:
//  https://blog.atulr.com/node-addon-guide/

namespace
{
void call_msg_send_in_node(Napi::Env env, Napi::Function callback, Napi::Reference<Napi::Value> *context, std::array<std::string,3> *data);
using SendMsg_TSFN = Napi::TypedThreadSafeFunction<Napi::Reference<Napi::Value>, std::array<std::string,3>, call_msg_send_in_node>;
boost::optional<SendMsg_TSFN> ns_message_to_node_js_callback;

enum class BrowseForDirStatus { NotStarted, Failed, Successful };

struct BrowseForDirData
{
  std::string token, window_title, window_msg;
  
  BrowseForDirStatus status;
  std::string result;
  
  std::mutex m;
  std::condition_variable cv;
};//struct BrowseForDirData

void call_dir_browse_in_node(Napi::Env env, Napi::Function callback, Napi::Reference<Napi::Value> *context, BrowseForDirData *data);
using BrowseDir_TSFN = Napi::TypedThreadSafeFunction<Napi::Reference<Napi::Value>, BrowseForDirData, call_dir_browse_in_node>;
boost::optional<BrowseDir_TSFN> ns_browse_for_directory_callback;


void call_msg_send_in_node(Napi::Env env, Napi::Function callback, Napi::Reference<Napi::Value> *context, std::array<std::string,3> *data)
{
  assert( data ); //we should always have valid data
  if( !data )
  {
    std::cerr << "call_msg_send_in_node: data is nullptr - unexpected." << std::endl;
    return;
  }
  
  std::unique_ptr<std::array<std::string,3>> data_holder( data );
  
  std::cout << "call_msg_send_in_node" << std::endl;
  
  // Check if TSFN has been aborted
  if( env == nullptr )
  {
    std::cerr << "call_msg_send_in_node: function has been aborted." << std::endl;
    return;
  }
    
  assert( callback != nullptr );
  
  if( callback == nullptr )
  {
    std::cerr << "call_msg_send_in_node: callback is null - unexpected." << std::endl;
    return;
  }
    
  std::array<std::string,3> &d = *data;
  callback.Call( context->Value(), { Napi::String::New(env,d[0]), Napi::String::New(env,d[1]),
                 Napi::String::New(env,d[2])} );
}//call_msg_send_in_node


void call_dir_browse_in_node(Napi::Env env, Napi::Function callback, Napi::Reference<Napi::Value> *context, BrowseForDirData *data)
{
  assert( data ); //we should always have valid data
  if( !data )
  {
    std::cerr << "call_dir_browse_in_node: data is nullptr - unexpected." << std::endl;
    return;
  }
  
  std::cout << "call_dir_browse_in_node" << std::endl;
  
  // Check if TSFN has been aborted
  if( env == nullptr )
  {
    std::cerr << "call_dir_browse_in_node: function has been aborted." << std::endl;
    return;
  }
  
  assert( callback != nullptr );
  
  if( callback == nullptr )
  {
    std::cerr << "call_dir_browse_in_node: callback is null - unexpected." << std::endl;
    return;
  }
  
  std::unique_lock<std::mutex> lk( data->m );
  
  data->status = BrowseForDirStatus::NotStarted;
  
  const std::string &token = data->token;
  const std::string &title = data->window_title;
  const std::string &msg = data->window_msg;
  std::string &result = data->result;
  
  Napi::Value rval = callback.Call( context->Value(), { Napi::String::New(env,token),
    Napi::String::New(env,title),
    Napi::String::New(env,msg)
  } );
  
  
  if( !rval.IsString() )
  {
    // This can happen if the user cancelled the operation.
    data->status = BrowseForDirStatus::Failed;
    std::cerr << "call_dir_browse_in_node: callback to JS didnt return a string." << std::endl;
  }else
  {
    data->status = BrowseForDirStatus::Successful;
    
    // For some reason directly using Napi calls to make UTF-8 string fails (empty string), so we
    //  will instead get the UTF-16 string, and convert to UTF-8
    //result = rval.As<Napi::String>().Utf8Value();
    
    const std::u16string val_as_utf16 = rval.As<Napi::String>().Utf16Value();
    std::wstring_convert<std::codecvt_utf8_utf16<char16_t>,char16_t> convert;
    result = convert.to_bytes(val_as_utf16);
    std::cout << "The call returned " << result.size() << " characters with value '" << result << "'" << std::endl;
  }
  
  lk.unlock();
  data->cv.notify_all();
}//call_dir_browse_in_node

}//namespace


namespace InterSpecAddOn
{
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
      std::u16string tmpdiru16 = info[0].ToString().Utf16Value();
      std::wstring tmpdir( std::begin(tmpdiru16), std::end(tmpdiru16) );
      _wputenv_s( L"TMPDIR", tmpdir.c_str() );
      _wputenv_s( L"WT_TMP_DIR", tmpdir.c_str() );
#else
        std::string tmpdir = info[0].ToString();
        setenv("TMPDIR", tmpdir.c_str(), 1);
        setenv("WT_TMP_DIR", tmpdir.c_str(), 1);
#endif // WIN32


    return Napi::Boolean::New( env, true );
  }

  Napi::Boolean setRequireSessionToken(const Napi::CallbackInfo& info)
  {
    Napi::Env env = info.Env();
  
    if (info.Length() < 1 || !info[0].IsBoolean() ) {
      Napi::TypeError::New(env, "setRequireSessionToken: Expected one bool").ThrowAsJavaScriptException();
    }
  
    const bool require = info[0].ToBoolean();
  
    interspec_set_require_session_token(require);
  
    return Napi::Boolean::New( env, true );
  }

  Napi::Boolean addPrimarySessionToken(const Napi::CallbackInfo& info)
  {
    Napi::Env env = info.Env();

    if (info.Length() < 1 || !info[0].IsString() ) {
        Napi::TypeError::New(env, "addPrimarySessionToken: Expected one strings").ThrowAsJavaScriptException();
    } 

    std::string token = info[0].ToString();

    interspec_add_allowed_primary_session_token( token.c_str() );

    return Napi::Boolean::New( env, true );
  }

  Napi::Boolean addExternalSessionToken(const Napi::CallbackInfo& info)
  {
    Napi::Env env = info.Env();
  
    if (info.Length() < 1 || !info[0].IsString() ) {
      Napi::TypeError::New(env, "addExternalSessionToken: Expected one strings").ThrowAsJavaScriptException();
    }
  
    std::string token = info[0].ToString();
  
    interspec_add_allowed_external_session_token( token.c_str() );
 
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


  Napi::Boolean usingElectronMenus(const Napi::CallbackInfo& info) 
  {
    Napi::Env env = info.Env();

    const bool isUsing = interspec_using_electron_menus();

    return Napi::Boolean::New( env, isUsing );
  }

  bool send_nodejs_message( const std::string &session_token, const std::string &msg_name, const std::string &msg_data )
  {
    std::cout << "In InterSpecAddOn::send_nodejs_message, for msg_name=" << msg_name << std::endl;
    
    if( !ns_message_to_node_js_callback.has_value() )
    {
      std::cerr << "InterSpecAddOn::send_nodejs_message: callback to JS not initialized" << std::endl;
      
      return false;
    }
    
    ns_message_to_node_js_callback->Acquire();
    
    //BlockingCall only blocks while putting into work queue - it does not block until the work
    //  is done, so we need to amke sure the next array will stick around for long enough, hence
    //  the call to new;  the \c call_msg_send_in_node function will take care of deleting it.
    std::array<std::string,3> *value = new std::array<std::string,3>{session_token,msg_name,msg_data};
    napi_status status = ns_message_to_node_js_callback->BlockingCall(value);
    
    ns_message_to_node_js_callback->Release();
    
    if( status == napi_ok )
    {
      std::cout << "InterSpecAddOn::send_nodejs_message: successfully called JS for msg_name=" << msg_name << std::endl;
      return true;
    }
    
    std::cerr << "InterSpecAddOn::send_nodejs_message: failed call to JS for msg_name=" << msg_name << std::endl;
    return false;
    
  }//send_nodejs_message(...)


  void setMessageToNodeJsCallback( const Napi::CallbackInfo &info )
  {
    Napi::Env env = info.Env();
    if( info.Length() != 1 || !info[0].IsFunction() )
    {
      Napi::TypeError::New(env, "setMessageToNodeJsCallback: expected function").ThrowAsJavaScriptException();
      return;
    }
    
    // See https://github.com/nodejs/node-addon-api/blob/main/doc/typed_threadsafe_function.md#example
    //  and also https://github.com/nodejs/node-addon-api/blob/main/doc/threadsafe.md
    //  for info about using this thread safe function callback in Napi.
    //  Basically, to call a node JS function, we have to do this from the Node main thread, but
    //  send_nodejs_message(...) is called from some random thread, so we have to use the
    //  TypedThreadSafeFunction mechanism (the ThreadSafeFunction mechanism didnt seem to work
    //  when I tried it - caused program aborts).
    Napi::Reference<Napi::Value> *context = new Napi::Reference<Napi::Value>( Persistent(info.This() ) );
    
    auto finalizer = []( Napi::Env, void *, Napi::Reference<Napi::Value> *ctx ) {
      // This finalizer is called at program exit.
      delete ctx;
    };//finalizer
    
    ns_message_to_node_js_callback = SendMsg_TSFN::New( env, info[0].As<Napi::Function>(),
                                                  "MessageToNodeJs",  0, 1, context, finalizer );

    std::cout << "setMessageToNodeJsCallback: Have set callback to JS" << std::endl;
  }//setMessageToNodeJsCallback


  void setBrowseForDirectoryCallback( const Napi::CallbackInfo &info )
  {
    Napi::Env env = info.Env();
    if( info.Length() != 1 || !info[0].IsFunction() )
    {
      Napi::TypeError::New(env, "setBrowseForDirectoryCallback: expected function").ThrowAsJavaScriptException();
      return;
    }

    Napi::Reference<Napi::Value> *context = new Napi::Reference<Napi::Value>( Persistent(info.This() ) );
  
    auto finalizer = []( Napi::Env, void *, Napi::Reference<Napi::Value> *ctx ) {
      delete ctx;
    };//finalizer
  
    ns_browse_for_directory_callback = BrowseDir_TSFN::New( env, info[0].As<Napi::Function>(),
                                             "BrowseForDirectory",  0, 1, context, finalizer );
  
    std::cout << "setBrowseForDirectoryCallback: Have set callback to JS" << std::endl;
  }//setBrowseForDirectoryCallback



  Napi::Boolean sendMessageToRenderer( const Napi::CallbackInfo &info )
  {
    Napi::Env env = info.Env();
    
    if( info.Length() < 2 || !(info[0].IsString() || info[0].IsNull() || info[0].IsUndefined()) || !info[1].IsString() ) {
      Napi::TypeError::New(env, "sendMessageToRenderer: Expected two or more strings").ThrowAsJavaScriptException();
    }
    
    std::string token;
    if( info[0].IsString() )
      token = info[0].ToString();
    else
      std::cerr << "sendMessageToRenderer: first argument is not string" << std::endl;
    
    std::string msg_name = info[1].ToString();
    std::string msg_data;
    if( info.Length() >= 3 && info[2].IsString() )
      msg_data = info[2].ToString();

    const bool success = ElectronUtils::handle_message_from_nodejs( token, msg_name, msg_data );
    
    return Napi::Boolean::New( env, success );
  }//void sendMessageToRenderer( const Napi::CallbackInfo &info )


bool browse_for_directory( const std::string &session_token,
                           const std::string &window_title,
                           const std::string &window_message,
                           std::function<void(std::string)> callback )
  {
    if( !ns_browse_for_directory_callback.has_value() )
    {
      std::cerr << "InterSpecAddOn::browse_for_directory: callback to JS not initialized" << std::endl;
      
      return false;
    }
    
    ns_browse_for_directory_callback->Acquire();
    
    BrowseForDirData data;
    data.token = session_token;
    data.window_title = window_title;
    data.window_msg = window_message;
    data.status = BrowseForDirStatus::NotStarted;
    
    std::unique_lock<std::mutex> lk(data.m);
    
    //BlockingCall only blocks while putting into work queue - it does not block until the work
    //  is done
    napi_status status = ns_browse_for_directory_callback->BlockingCall( &data );
    
    //std::cout << "About to start waiting on lock" << std::endl;
    
    data.cv.wait( lk, [&data]{ return data.status != BrowseForDirStatus::NotStarted; } );
    
    //std::cout << "Done waiting on lock" << std::endl;
    
    ns_browse_for_directory_callback->Release();
    
    if( status == napi_ok )
    {
      switch( data.status )
      {
        case BrowseForDirStatus::NotStarted:
          std::cerr << "InterSpecAddOn::browse_for_directory: call never started?  Not calling callback." << std::endl;
          return false;
          
        case BrowseForDirStatus::Failed:
          std::cerr << "InterSpecAddOn::browse_for_directory: call failed.  Not calling callback." << std::endl;
          return false;
          
        case BrowseForDirStatus::Successful:
          //std::cout << "InterSpecAddOn::browse_for_directory: successfully called=" << data.result << std::endl;
          callback( data.result );
          return true;
      }//switch( data.result )
    }
    
    std::cerr << "InterSpecAddOn::send_nodejs_message: failed call to JS (napi status=" << status << ")" << std::endl;
    return false;
  }
}//namespace InterSpecAddOn



Napi::Object InitAll(Napi::Env env, Napi::Object exports) {
  exports.Set( "startServingInterSpec", Napi::Function::New(env, InterSpecAddOn::startServingInterSpec) );
  exports.Set( "killServer", Napi::Function::New(env, InterSpecAddOn::killServer) );
  exports.Set( "openFile", Napi::Function::New(env, InterSpecAddOn::openFile ));

  exports.Set( "setTempDir", Napi::Function::New(env, InterSpecAddOn::setTempDir ));
  
  exports.Set( "setRequireSessionToken", Napi::Function::New(env, InterSpecAddOn::setRequireSessionToken ));
  exports.Set( "addPrimarySessionToken", Napi::Function::New(env, InterSpecAddOn::addPrimarySessionToken ));
  exports.Set( "addExternalSessionToken", Napi::Function::New(env, InterSpecAddOn::addExternalSessionToken ));
  
  
  exports.Set( "removeSessionToken", Napi::Function::New(env, InterSpecAddOn::removeSessionToken ));

  exports.Set( "usingElectronMenus", Napi::Function::New(env, InterSpecAddOn::usingElectronMenus ));
  
  exports.Set( "setMessageToNodeJsCallback", Napi::Function::New(env, InterSpecAddOn::setMessageToNodeJsCallback));
  exports.Set( "setBrowseForDirectoryCallback", Napi::Function::New(env, InterSpecAddOn::setBrowseForDirectoryCallback));
  
  exports.Set( "sendMessageToRenderer", Napi::Function::New(env, InterSpecAddOn::sendMessageToRenderer));
  
  return exports;
}


/**
* This code defines the entry-point for the Node addon, it tells Node where to go
* once the library has been loaded into active memory. The first argument must
* match the "target" in our *binding.gyp*. Using NODE_GYP_MODULE_NAME ensures
* that the argument will be correct, as long as the module is built with
* node-gyp (which is the usual way of building modules). The second argument
* points to the function to invoke. The function must not be namespaced.
 // Note that the actual name used doesnt seem to really matter, so I'll just use InsterSpecAddOn
*/
NODE_API_MODULE(InsterSpecAddOn, InitAll)
