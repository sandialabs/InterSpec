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
#include <iostream>
#include <stdexcept>

// Some includes to get terminal width (and UTF-8 cl arguments on Windows)
#if defined(__APPLE__) || defined(linux) || defined(unix) || defined(__unix) || defined(__unix__)
#include <sys/ioctl.h>
#include <unistd.h>
#elif defined(_WIN32)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN 1
#include <Windows.h>
#include <stdio.h>
#include <direct.h>
#include <shellapi.h>
#endif

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
#ifdef __APPLE__
#include <mach-o/dyld.h>  //for _NSGetExecutablePath
#elif( defined(_WIN32) )
#include <libloaderapi.h> //for GetModuleFileNameW
#include <Shlobj.h>       //for CSIDL_APPDATA and SHGetKnownFolderPath
#else
#include <limits.h>  //for PATH_MAX
#endif
#endif //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )

#include <boost/filesystem.hpp>  //for boost::filesystem::is_symlink and read_symlink


#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"

using namespace std;

namespace AppUtils
{
  void split_uri( string url, string &host, string &path, string &query, string &fragment )
  {
    // But also see rfc3986 for a regex that might work, or be the start of making one that
    //  will work, to parse the URL.
    
    host = path = query = fragment = "";
    
    //Get rid of (optional) leading "interspec://", so URL will look like: 'drf/specify?v=1&n=MyName&"My other Par"'
    const string interspec_scheme = "interspec://";
    const string spectrum_scheme = "RADDATA://";
    if( SpecUtils::istarts_with(url, interspec_scheme) )
      url = url.substr(interspec_scheme.size());
    else if( SpecUtils::istarts_with(url, spectrum_scheme) )
      url = url.substr(spectrum_scheme.size());
    
    
    const string::size_type q_pos = url.find( '?' );
    //string::size_type q_end_pos = q_pos + 1;
    
    const string::size_type hash_pos = url.find( '#' );
    const string::size_type path_end = std::min( hash_pos, q_pos );
    
    /*
    if( q_pos == string::npos )
    {
      // Since '?' isnt a QR alphanumeric code also allow the URL encoded value of it, "%3F"
     q_pos = url.find( "%3F" );
      if( q_pos == string::npos )
     q_pos = url.find( "%3f" );
      
      if( q_pos != string::npos )
        q_end_pos = q_pos + 3;
    }//
     */
      
    //if( q_pos == string::npos )
    //  throw runtime_error( "App URL did not contain the host/path component that specifies the intent"
    //                       " of the URL (e.g., what tool to use, or what info is contained in the URL)." );
    
    // the q_pos+1 should be safe, even if q_pos is last character in string.
    //if( url.find(q_end_pos, '?') != string::npos )
    //  throw runtime_error( "App URL contained more than one '?' character, which isnt allowed" );

    const string host_path = url.substr( 0, path_end );
    const string::size_type host_end = host_path.find( '/' );
    
    // If "host/path" is only one component (i.e., no slash), then populate host,
    //  and have path be empty
    host = (host_end == string::npos) ? host_path : host_path.substr(0,host_end);
    path = (host_end == string::npos) ? string("") : host_path.substr(host_end+1);
    
    if( (q_pos != string::npos) && (q_pos < hash_pos) )
      query = url.substr( q_pos+1, (hash_pos==string::npos) ? string::npos : (hash_pos - q_pos - 1) );
      
    if( hash_pos != string::npos )
      fragment = url.substr( hash_pos + 1 );
  }//void split_uri( string url, string &host, string &path, string &query )

  
  std::map<std::string,std::string> query_str_key_values( const std::string &query_str )
  {
    map<string,string> parts;
    vector<string> components;
    SpecUtils::split( components, query_str, "&" );
    
    for( string comp : components )
    {
      SpecUtils::trim( comp );
      if( comp.empty() )
        continue;
      
      auto pos = comp.find('=');
      if( pos == string::npos )
        pos = comp.find(':');
      
      if( pos == string::npos )
        continue;
      
      string key = comp.substr(0,pos);
      string value = comp.substr(pos+1);
      SpecUtils::trim(key);
      SpecUtils::to_upper_ascii(key);
      
      if( key.empty() )
        throw runtime_error( "fromAppUrl: query portion '" + comp + "' has empty name" );
      
      if( parts.count(key) )
        throw runtime_error( "fromAppUrl: query portion contains duplicate key '" + key + "'" );
      
      parts[key] = value;
    }//for( string comp : components )
    
    return parts;
  }//std::map<std::string,std::string> split_query_str( const std::string &query )
  
  /*
  vector<pair<string,string>> query_key_values( const string &query_str )
  {
    vector<pair<string,string>> parts;
    vector<string> components;
    SpecUtils::split( components, query_str, "&" );
    
    for( string comp : components )
    {
      SpecUtils::trim( comp );
      if( comp.empty() )
        continue;
      
      auto pos = comp.find('=');
      if( pos == string::npos )
        pos = comp.find(':');
      
      string key = comp.substr(0,pos);
      string value = (pos == string::npos) ? string("") : comp.substr(pos+1);
      SpecUtils::trim(key);
      SpecUtils::to_upper_ascii(key);
      
      parts.emplace_back( key, value );
    }//for( string comp : components )
    
    return parts;
  }//vector<pair<string,string>> query_key_values( const string &query );
   */
  
  
#if( USE_BATCH_TOOLS || BUILD_AS_LOCAL_SERVER )
#if defined(__APPLE__) || defined(unix) || defined(__unix) || defined(__unix__)
unsigned terminal_width()
{
  winsize ws = {};
  if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) <= -1)
    return 80;
  unsigned w = (ws.ws_col);
  return std::max( 40u, w );
}
#elif defined(_WIN32)
unsigned terminal_width()
{
  HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
  if( handle == INVALID_HANDLE_VALUE )
    return 80;
  
  CONSOLE_SCREEN_BUFFER_INFO info;
  if( !GetConsoleScreenBufferInfo(handle, &info) )
    return 80;
  
  return unsigned(info.srWindow.Right - info.srWindow.Left);
}
#else
static_assert( 0, "Not unix and not win32?  Unsupported getting terminal width" );
#endif
#endif //#if( USE_BATCH_TOOLS || BUILD_AS_LOCAL_SERVER )
  
  
uint32_t compile_date_as_int()
{
  //The below YEAR MONTH DAY macros are taken from
  //http://bytes.com/topic/c/answers/215378-convert-__date__-unsigned-int
  //  and I believe to be public domain code
#define YEAR ((((__DATE__ [7] - '0') * 10 + (__DATE__ [8] - '0')) * 10 \
+ (__DATE__ [9] - '0')) * 10 + (__DATE__ [10] - '0'))
#define MONTH (__DATE__ [2] == 'n' && __DATE__ [1] == 'a' ? 0 \
: __DATE__ [2] == 'b' ? 1 \
: __DATE__ [2] == 'r' ? (__DATE__ [0] == 'M' ? 2 : 3) \
: __DATE__ [2] == 'y' ? 4 \
: __DATE__ [2] == 'n' ? 5 \
: __DATE__ [2] == 'l' ? 6 \
: __DATE__ [2] == 'g' ? 7 \
: __DATE__ [2] == 'p' ? 8 \
: __DATE__ [2] == 't' ? 9 \
: __DATE__ [2] == 'v' ? 10 : 11)
#define DAY ((__DATE__ [4] == ' ' ? 0 : __DATE__ [4] - '0') * 10 + (__DATE__ [5] - '0'))
  
  return YEAR*10000 + (MONTH+1)*100 + DAY;
}//uint32_t compile_date_as_int()
  
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  
std::string current_exe_path()
{
#ifdef __APPLE__
    char path_buffer[PATH_MAX + 1] = { '\0' };
    uint32_t size = PATH_MAX + 1;
    
    if (_NSGetExecutablePath(path_buffer, &size) != 0)
      throw runtime_error( "_NSGetExecutablePath failed" );
    
    path_buffer[PATH_MAX] = '\0'; // JIC
    const string exe_path = path_buffer;
#elif( defined(_WIN32) )
    //static_assert( 0, "Need to test this EXE path stuff on Windows..." );
    wchar_t wbuffer[2*MAX_PATH];
    const DWORD len = GetModuleFileNameW( NULL, wbuffer, 2*MAX_PATH );
    if( len <= 0 )
      throw runtime_error( "Call to GetModuleFileName failed" );
    
    const string exe_path = SpecUtils::convert_from_utf16_to_utf8( wbuffer );
#else // if __APPLE__ / Win32 / else
    
    char path_buffer[PATH_MAX + 1] = { '\0' };
    const ssize_t ret = readlink("/proc/self/exe", path_buffer, PATH_MAX);
    
    if( (ret == -1) || (ret > PATH_MAX) )
      throw runtime_error( "Failed to read line" );
    
    assert( ret < PATH_MAX );
    path_buffer[ret] = '\0';
    path_buffer[PATH_MAX] = '\0'; // JIC
    const string exe_path = path_buffer;
#endif // else not __APPLE__
  
  return exe_path;
}//std::string current_exe_path()
  
  
bool locate_file( string &filename, const bool is_dir,
                 size_t max_levels_up, const bool include_path )
{
  auto check_exists = [is_dir]( const string &name ) -> bool {
    return is_dir ? SpecUtils::is_directory(name) : SpecUtils::is_file(name);
  };//auto check_exists
  
  if( SpecUtils::is_absolute_path(filename) )
    return check_exists(filename);
  
  // Check if path is there, relative to CWD
  if( check_exists(filename) )
    return true;
  
  if( SpecUtils::icontains(filename,"http://") || SpecUtils::icontains(filename,"https://") )
    return false;
  
  // We'll look relative to the executables path, but note that if we started from a symlink, I
  //  think it will resolve relative to actual executable
  try
  {
    const string exe_path = current_exe_path();
    
    string canonical_exe_path = exe_path;
    if( !SpecUtils::make_canonical_path(canonical_exe_path) )
      throw runtime_error( "Failed to make trial_parent_path canonical" );
    
    string trial_parent_path = canonical_exe_path;
    
    try
    {
      do
      {
        trial_parent_path = SpecUtils::parent_path( trial_parent_path );
        const string trialpath = SpecUtils::append_path( trial_parent_path, filename );
        
        if( check_exists(trialpath) )
        {
          filename = trialpath;
          return true;
        }
        
        if( boost::filesystem::is_symlink(trialpath) )
        {
          const string unsym_trialpath = boost::filesystem::read_symlink(trialpath).string<string>();
          if( check_exists(unsym_trialpath) )
          {
            filename = unsym_trialpath;
            return true;
          }
        }//if( is symlink )
        
        max_levels_up = (max_levels_up > 0) ? (max_levels_up - 1) : 0;
      }while( max_levels_up > 0 );
    }catch( std::exception & )
    {
      //cerr << "Caught exception searching exe's parent dirs: " << e.what() << endl;
    }//try / catch
  }catch( std::exception & )
  {
    //cerr << "Caught exception: " << e.what() << endl;
  }
  
  if( !include_path )
    return false;
  
  // Check relative to environments "PATH"
  const char *env_path = std::getenv("PATH");
  if( !env_path )
    return false;
  
  // grab the PATH system variable, and check if the input path is relative to any of those
  vector<string> paths;
#if( defined(_WIN32) )
  const char *delims = ";";
#else
  const char *delims = ":";
#endif
  SpecUtils::split( paths, env_path, ";" );
  
  for( string base : paths )
  {
    const string trialpath = SpecUtils::append_path( base, filename );
    
    if( check_exists(filename) )
    {
      filename = trialpath;
      return true;
    }
  }//for( string base : paths )
  
  // Out of luck, we failed if we're here.
  
  return false;
}//bool locate_file( ... )
#endif //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  
  
std::string file_contents( const std::string &filename )
{
  //Copied from SpecUtils::load_file_data( const char * const filename, std::vector<char> &data );
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  basic_ifstream<char> stream(wfilename.c_str(), ios::binary);
#else
  basic_ifstream<char> stream(filename.c_str(), ios::binary);
#endif

  if (!stream)
    throw runtime_error(string("cannot open file ") + filename);
  stream.unsetf(ios::skipws);

  // Determine stream size
  stream.seekg(0, ios::end);
  size_t size = static_cast<size_t>( stream.tellg() );
  stream.seekg(0);

  string data;
  data.resize( size );
  stream.read(&data.front(), static_cast<streamsize>(size));
  
  return data;
}//std::string file_contents( const std::string &filename )
  
  
#ifdef _WIN32
/** Get command line arguments encoded as UTF-8.
    This function just leaks the memory
 
 Note that environment variables are not in UTF-8, we could account for this
 similar to:
 wchar_t *wenvstrings = GetEnvironmentStringsW();
 ...
 */
void getUtf8Args( int &argc, char ** &argv )
{
  LPWSTR *argvw = CommandLineToArgvW( GetCommandLineW(), &argc );
  if( !argvw )
  {
    std::cout << "CommandLineToArgvW failed - good luck" << std::endl;
    return ;
  }
  
  argv = (char **)malloc(sizeof(char *)*argc);
    
  for( int i = 0; i < argc; ++i)
  {
    printf("Argument: %d: %ws\n", i, argvw[i]);
  
    const std::string asutf8 = SpecUtils::convert_from_utf16_to_utf8( argvw[i] );
    argv[i] = (char *)malloc( sizeof(char)*(asutf8.size()+1) );
    strcpy( argv[i], asutf8.c_str() );
  }//for( int i = 0; i < argc; ++i)

  // Free memory allocated for CommandLineToArgvW arguments.
  LocalFree(argvw);
}//void getUtf8Args()


void getUtf8Args( int &argc, wchar_t **argvw, char **&argv )
{
  argv = (char **)malloc( sizeof( char * ) * argc );

  for( int i = 0; i < argc; ++i )
  {
    //printf("Argument: %d: %ws\n", i, argvw[i]);
    const std::string asutf8 = SpecUtils::convert_from_utf16_to_utf8( argvw[i] );
    argv[i] = (char *)malloc( sizeof( char ) * (asutf8.size() + 1) );
    strcpy( argv[i], asutf8.c_str() );
  }//for( int i = 0; i < argc; ++i)
}//void processCustomArgs()


void cleanupUtf8Args( int &argc, char **&argv )
{
  for( int i = 0; i < argc; ++i )
    free( argv[i] );
  free( argv );
  argc = 0;
  argv = nullptr;
}//void cleanupUtf8Args( int &argc, char **&argv )


std::string user_data_dir()
{ 
  //  Should we use SHGetKnownFolderPath with FOLDERID_RoamingAppData instead?
  TCHAR ppszPath[MAX_PATH];
  HRESULT hr = ::SHGetFolderPath( nullptr, CSIDL_APPDATA, nullptr, SHGFP_TYPE_CURRENT, ppszPath );
  if( hr == E_FAIL ) // Folder doesnt exist, get default value?  
    hr = ::SHGetFolderPath( nullptr, CSIDL_APPDATA, nullptr, SHGFP_TYPE_DEFAULT, ppszPath );
 
  if( hr != S_OK )
    throw runtime_error( "Couldnt get APP data directory" );

  int len = MultiByteToWideChar( CP_ACP, 0, ppszPath, -1, NULL, 0 );
  std::wstring wstr( len, 0 );
  MultiByteToWideChar( CP_ACP, 0, ppszPath, -1, &wstr[0], len );


  int utf8Len = WideCharToMultiByte( CP_UTF8, 0, &wstr[0], -1, NULL, 0, NULL, NULL );
  std::string utf8Str( utf8Len, 0 );
  WideCharToMultiByte( CP_UTF8, 0, &wstr[0], -1, &utf8Str[0], utf8Len, NULL, NULL );
  
  return utf8Str;
}//std::string AppUtils::user_data_dir()
#endif
}//namespace AppUtils
