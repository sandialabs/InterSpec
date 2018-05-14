
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <Wt/WText>
#include <Wt/WEnvironment>
#include <Wt/WApplication>

#include "InterSpec/SrbActivityApp.h"
#include "InterSpec/SrbHeaderFooter.h"


using namespace Wt;
using namespace std;


string prependPathName( string dir, string name )
{
  //A translation of the TUnixSystem::PrependPathName(...) function
  if (name.empty() || name == ".")
  {
    if( !dir.empty() )
    {
      name = dir;
      if( dir[dir.length() - 1] != '/' )
        name += '/';
    } else name = "";
    return name;
  }

  if( dir.empty() ) dir = "/";
  else if( dir[dir.length() - 1] != '/' )
    name = '/' + name;
  name = dir + name;

  return name;
}//string prependPathName( string dir, string name )


string workingDirectory()
{
  char cwd[2024];
  if( ::getcwd(cwd, 2024) == 0 )
  {
    cerr << "workingDirectory(...): getcwd() failed" << endl;
    return "/";
  }
  return cwd;
}//string workingDirectory()

bool accessPathName( const string &path, int mode = R_OK )
{
  //XXX - !returns false if you can access the file!
  if (::access(path.c_str(), mode ) == 0)  //R_OK in unistd.h
    return false;
  return true;
}

const char *tempDirectory()
{
  //A translation of the TUnixSystem::TempDirectory(...) function
  //
  // Return a user configured or systemwide directory to create
  // temporary files in.

  const char *dir = ::getenv("TMPDIR");
  if (!dir || accessPathName(dir, W_OK))
    dir = "/tmp";

  return dir;
}//const char *tempDirectory() const


std::string getUserNameFromEnvironment()
{
  string remoteUser = wApp->environment().getCgiValue( "REMOTE_USER" );
  boost::algorithm::erase_all( remoteUser, " " );
  boost::algorithm::to_lower(remoteUser);

  return remoteUser;
}//std::string getUserNameFromEnvironment() const


FILE *tempFileName( string &base, const char *dir = 0 )
{
  //A translation of the TUnixSystem::TempFileName(...) function
  //
  // Create a secure temporary file by appending a unique
  // 6 letter string to base. The file will be created in
  // a standard (system) directory or in the directory
  // provided in dir. The full filename is returned in base
  // and a filepointer is returned for safely writing to the file
  // (this avoids certain security problems). Returns 0 in case
  // of error.

  string b = prependPathName(dir ? dir : tempDirectory(), base);
  base = b;
  base += "XXXXXX";


  char *arg = strdup(base.c_str());
  int fd = mkstemp(arg);
  base = arg;
  delete [] arg;

  if (fd == -1) {
    cerr << "tempFileName " << base << " error" << endl;
    return 0;
  } else {
    FILE *fp = fdopen(fd, "w+");
    if (fp == 0)
      cerr << "tempFileName: error converting filedescriptor (" << fd << ")" << endl;;
    return fp;
  }
}


int exectutePhpScript( std::string &results,
                                        const std::string &path_to_script,
                                        const std::string &script_name,
                                        const string &path_to_srb_common )
{

  results = "";

  string script = script_name;
  string script_location = path_to_script;
  if( path_to_script.empty() ) script_location = ".";

  script = prependPathName( script_location, script );

  if( script.empty() ||  script[0] != '/' )
  {
    const string work_dir = workingDirectory();
    script = prependPathName( work_dir, script );
  }//if( we need to make the path absolute )

  if( accessPathName( script ) ) //returns false if you CAN access the file
  {
    string errmsg = "Couln't access the file '" + script + "'";
    cerr << errmsg << endl; //so it will show up in apaches log
    return EXIT_FAILURE;
  }//if( we can touch the file ).

  const string host = wApp->environment().hostName();
  const string user = getUserNameFromEnvironment();
  string uri = host;
  const size_t ind = uri.find_first_of( "/" );
  if( ind != string::npos ) uri = uri.substr( ind );
  else                      uri = "";

  //TODO: It doesnt appear the below is actually working!
  const string file_contents = "<?php $PATH_TO_COMMON=\""
  + path_to_srb_common +
  "/\"; "
  "$HTTP_HOST=\""     + host + "\"; "  //Contents of the Host: header from the current request, if there is one.
  "$PHP_AUTH_USER=\"" + user + "\"; "  //When doing HTTP authentication this variable is set to the username provided by the user.
  "$REQUEST_URI=\""   + uri  + "\"; "  //The URI which was given in order to access this page; for instance, '/index.html'.
  "$_SERVER['" + host + "']=\"" + host + "\"; "
  "$_SERVER['" + user + "']=\"" + user + "\"; "
  "$_SERVER['" + uri + "']=\"" + uri + "\"; "
  " include_once( '"
  + script
  + "' ); \?>";

  string temp_file_name;
  FILE *temp_file = tempFileName( temp_file_name );

  if( temp_file == NULL )
  {
    string errmsg = "exectutePhpScript(...): couldnt open temp file";
    cerr << errmsg << endl;
    return EXIT_FAILURE;
  }//if( temp_file == NULL )

  fprintf( temp_file, "%s", file_contents.c_str() );

  fflush( temp_file );
  pclose( temp_file );

  string php_command = "php " + temp_file_name;
  FILE *phpheader = popen( php_command.c_str(), "r" );

  int close_status = EXIT_FAILURE;

  if( phpheader != NULL )
  {
    int character;
    while( (character = getc(phpheader)) != EOF ) results += (char)character;
    close_status = pclose( phpheader );

    if( close_status != EXIT_SUCCESS )
    {
      string errmsg = "The php--->c++/tex status gave a exit status of "
      + boost::lexical_cast<string>(close_status)
      + ", you should probably fix this '"
      + php_command + "'";
      cerr << errmsg << endl; //so it will show up in apaches log
    }//if( close_status != EXIT_SUCCESS )
  }else
  {
    string errmsg = "The exectutePhpScript failed to execute '" + php_command + "'";
    cerr << errmsg << endl;
  }//if( opened FILE output ) / else

  ::unlink( temp_file_name.c_str() );

  return close_status;

  return 0;
}//int exectutePhpScript( const std::string &command, std::string &result )



SrbHeader::SrbHeader( const std::string &appTitle, const std::string &baseurl,
                     WContainerWidget *parentParam )
: WContainerWidget(parentParam),
m_headerText( NULL ),
m_applicationTitle( NULL )
{
  setId( "SrbHeader" );
  //  setAttributeValue( "style", "" );
  //  setInline(true); //set to <span>
  
  string content = "";
  
  const string srb_common_dir = getSrbComonResourcePath();
  
  // popen is on POSIX systems. _popen is the Windows equivalent. FIREIRON041
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  // As the shell isn't globally available, it needs to be pointed to the active directory
  // for it to work properly.
  //const string php_command = "\"C:/PHP/php.exe\" " +
  boost::filesystem::current_path().parent_path().generic_string() + "/" +
  srb_common_dir + "/assets/includes/global-nav.php";
  //FILE *phpheader = _popen( php_command.c_str(), "r" );
#else
  const string php_command = "php " + srb_common_dir + "/assets/includes/global-nav.php";
  FILE *phpheader = popen( php_command.c_str(), "r" );
#endif // #if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  
  if( phpheader != NULL )
  {
    int character;
    while( (character = getc(phpheader)) != EOF ) content += (char)character;
    
    // pclose is on POSIX systems. _pclose is the Windows equivalent. FIREIRON041
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
    const int close_status = _pclose( phpheader );
#else
    const int close_status = pclose( phpheader );
#endif // #if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
    
    if( close_status != EXIT_SUCCESS )
    {
      string errmsg = "The php--->c++ srb header status gave a exit status of "
      + boost::lexical_cast<string>(close_status)
      + ", you should probably fix this - the php command was '"
      + php_command + "'";
      
      cerr << errmsg << endl; //so it will show up in apaches log
    }//if( close_status != EXIT_SUCCESS )
  } else {
    cerr << "\n\tFailed to execute '" << php_command << "'\n";
  }//if( opened FILE output ) / else
  
  if( content == "" ) //if for some reason the php command failed, lets default to something
  {
    string errmsg = "InterSpec is NOT using the php-header for some reason; this should be checked into";
    cerr << errmsg << endl;
    
    content =
    "\n  <ul>\n"\
    "    <li><a href=\"https://" + baseurl + "/srb/\" title=\"DHS SRB Webtools Homepage\" id=\"webtools-logo\"><div><span>DHS SRB Webtools Homepage</span></div></a></li> \n"\
    "    <li id=\"select-another-tool\">\n"\
    "      <div id=\"select-header\"><a id=\"select-header-link\">Select another tool</a></div>\n"\
    "      <ul id=\"select-list\" style=\"display:none;\">\n"\
    "        <li><a href=\"https://" + baseurl + "/srb/reports/cbp/\">Automated Reports</a></li>\n"\
    "        <li><a href=\"https://" + baseurl + "/Cambio/\">Cambio</a></li>\n"\
    "        <li><a href=\"https://" + baseurl + "/srb/dhs-server/\">DHS Data Management</a></li>\n"\
    "        <li><a href=\"https://" + baseurl + "/GADRAS\">GADRAS</a></li>\n"\
    "        <li><a href=\"https://" + baseurl + "/srb/viz/lsslogs/\">LSS logs: Google Earth Viz</a></li>\n"\
    "        <li><a href=\"https://" + baseurl + "/wiki/srb/siteinfo/\">POE Wiki</a></li>\n"\
    "      </ul>\n"\
    "    </li>\n"\
    "  </ul>\n\n";
  }//if( content == "" )
  
  m_headerText = new WText( content, Wt::XHTMLUnsafeText, this );
  m_headerText->setId( "global-nav" );
  m_headerText->setInline(false); //set to div
  
  if( appTitle != "" )
  {
    content =
    "\n  <span id=\"site-logo\"></span>\n"\
    "  <h1><a href=\"#\">" + appTitle + "</a></h1>\n\n";
    m_applicationTitle = new WText( content, Wt::XHTMLUnsafeText, this );
    m_applicationTitle->setId( "site-title" );
    m_applicationTitle->setInline(false); //set to div
  } // if( appTitle != "" )
} // SrbHeader constructor


SrbFooter::SrbFooter( WContainerWidget *parentParam )
: WContainerWidget(parentParam),
m_logoText( NULL ),
m_footerText( NULL ),
m_footerCloseText( NULL )
{
  //  setId( "SrbFooter" );
  setId( "footer" );
  setInline(false);
  
  // sponsors
  const string srb_common_url = getSrbComonResourceUrl();
  const string srb_common_dir = getSrbComonResourcePath();
  string php_dir = srb_common_dir + "/assets/includes/";
  
  string sponsors_html = "";
  string global_sponsors_html = "";
  exectutePhpScript( global_sponsors_html, php_dir, "global-sponsors.php", srb_common_url );
  sponsors_html += global_sponsors_html;
  WText *sponsors = new WText( sponsors_html, Wt::XHTMLUnsafeText, this );
  sponsors->setInline( false );
  sponsors->setId( "sponsors" );
  
  // footer
  string global_footer_html = "";
  const int rcode = exectutePhpScript( global_footer_html, php_dir, "global-footer.php", srb_common_url );
  if( rcode == EXIT_SUCCESS )
  {
    m_footerText = new WText( global_footer_html, Wt::XHTMLUnsafeText, this );
    m_footerText->setId( "footer-content" );
    m_footerText->setInline(false); //set to div
  }//if( rcode == EXIT_SUCCESS )
  
}//SrbFooter( constructor )



