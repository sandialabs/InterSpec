/**
 SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
 Copyright (C) 2016 William Johnson
 
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


#include "SpecUtils_config.h"

// Block out some warnings occurring in xutility.
// warning C4996: function call with parameters that may be unsafe -
#pragma warning(disable:4996)

#include <ctime>
#include <string>
#include <vector>
#include <locale>
#include <limits>
#include <time.h>
#include <ctype.h>
#include <stdio.h>
#include <fstream>
#include <errno.h>
#include <string.h>
#include <iostream>
#include <string.h>
#include <algorithm>


#include "rapidxml/rapidxml.hpp"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include <Lmcons.h>
#include <direct.h>
#include <io.h>
#elif __APPLE__
#include <sys/time.h>
#include <sys/sysctl.h>
#else
#include <sys/time.h>
#include <unistd.h>
#endif


//#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
#include <sys/stat.h>
//#endif

//#if( ANDROID )
//#include <boost/chrono/io/time_point_io.hpp>  //for internal_timegm
//#endif

#if( !SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
#include <boost/filesystem.hpp>
#endif


#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>

#include "SpecUtils/UtilityFunctions.h"

#if( ANDROID )
#include <jni.h>
#endif

using namespace std;

namespace
{
  //adapted from http://stackoverflow.com/questions/3152241/case-insensitive-stdstring-find
  template<typename charT>
  struct char_iequal
  {
    char_iequal( const std::locale &loc ) : m_loc(loc) {}
    bool operator()(charT ch1, charT ch2) {
      return std::toupper(ch1, m_loc) == std::toupper(ch2, m_loc);
    }
  private:
    const std::locale &m_loc;
  };
  
  //The bellow may not be instrinsic friendly; might be interesting to look and
  // see if we could get better machine code by making instrinic friendly (maybe
  // even have the float parsing functions inline)
  inline bool is_in( const char val, const char *delim )
  {
    while( *delim )
    {
      if( val == *delim )
        return true;
      ++delim;
    }
    return false;
  }
  
  inline const char *next_word( const char *input, const char *delim )
  {
    while( *input && is_in(*input, delim) )
      ++input;
    return input;
  }
  
  inline const char *end_of_word( const char *input, const char *delim )
  {
    while( *input && !is_in(*input, delim) )
      ++input;
    return input;
  }
}//namespace


namespace UtilityFunctions
{
  //Templated boost functions used multiple times tend to take up a ton of space
  //  in the final executable, so to keep it so that only one instaniation of
  //  each comonly used boost function is made, we'll force this with the below.
  //It also makes it a bit easier to eliminate boost from the dependancies later
  //  on if this is desired.
  //This results in a reduction of the generated code size for this file by
  // 71.1 kb on Win74 MinSizeRel
  bool istarts_with( const std::string &line, const char *label )
  {
    const size_t len1 = line.size();
    const size_t len2 = strlen(label);
    
    if( len1 < len2 )
      return false;
  
    const bool answer = ::rapidxml::internal::compare( line.c_str(), len2, label, len2, false );
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::istarts_with( line, label );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
                "Got %i when should have got %i for label '%s' and string '%s'",
                int(answer), int(correctAnswer), label, line.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }//istarts_with(...)
  
  bool istarts_with( const std::string &line, const std::string &label )
  {
    const size_t len1 = line.size();
    const size_t len2 = label.size();
    
    if( len1 < len2 )
      return false;
    
    const bool answer = ::rapidxml::internal::compare( line.c_str(), len2, label.c_str(), len2, false );
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::istarts_with( line, label );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), label.c_str(), line.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }//istarts_with(...)
  
  
  bool starts_with( const std::string &line, const char *label )
  {
    const size_t len1 = line.size();
    const size_t len2 = strlen(label);
    
    if( len1 < len2 )
      return false;
    
    const bool answer = ::rapidxml::internal::compare( line.c_str(), len2, label, len2, true );
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::istarts_with( line, label );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), label, line.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }//istarts_with(...)
  
  
  bool iends_with( const std::string &line, const std::string &label )
  {
    const size_t len1 = line.size();
    const size_t len2 = label.size();
    
    if( len1 < len2 )
      return false;
    
    const char * const lineend = line.c_str() + (len1 - len2);
    
    const bool answer = ::rapidxml::internal::compare( lineend, len2, label.c_str(), len2, false );
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::iends_with( line, label );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), label.c_str(), line.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }
  
  void erase_any_character( std::string &line, const char *chars_to_remove )
  {
    if( !chars_to_remove )
      return;
    
    auto should_remove = [chars_to_remove]( const char &c ) -> bool {
      for( const char *p = chars_to_remove; *p; ++p )
        if( c == *p )
          return true;
      return false;
    };//should_remove
    
    line.erase( std::remove_if(line.begin(), line.end(), should_remove), line.end() );
  }//void erase_any_character( std::string &line, const char *chars_to_remove );
  
  bool icontains( const char *line, const size_t length,
                  const char *label, const size_t labellen )
  {
    const char *start = line;
    const char *end = start + length;
    const char *it = std::search( start, end, label, label+labellen,
                                 char_iequal<char>(std::locale()) );
    const bool answer = (it != end);
    
#if(PERFORM_DEVELOPER_CHECKS)
    const string cppstr = string( start, end );
    const string cpplabel = string( label, label + labellen );
    const bool correctAnswer = boost::algorithm::icontains( cppstr, cpplabel );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), cpplabel.c_str(), cppstr.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }//icontains(...)
  
  bool icontains( const std::string &line, const char *label )
  {
    const size_t labellen = strlen(label);
    return icontains( line.c_str(), line.size(), label, labellen );
  }//icontains(...)
  
  
  bool icontains( const std::string &line, const std::string &label )
  {
    return icontains( line.c_str(), line.size(), label.c_str(), label.size() );
  }//icontains(...)
  
  
  bool contains( const std::string &line, const char *label )
  {
    const bool answer = (line.find(label) != string::npos);
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::contains( line, label );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), label, line.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }//contains(...)
  
  bool iequals( const char *str, const char *test )
  {
    const bool answer = ::rapidxml::internal::compare( str, strlen(str), test, strlen(test), false );
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::iequals( str, test );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), test, str );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }//bool iequals
  
  bool iequals( const std::string &str, const char *test )
  {
    const bool answer = ::rapidxml::internal::compare( str.c_str(), str.length(), test, strlen(test), false );
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::iequals( str, test );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), test, str.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }
  
  bool iequals( const std::string &str, const std::string &test )
  {
    const bool answer = ::rapidxml::internal::compare( str.c_str(), str.size(), test.c_str(), test.size(), false );
    
#if(PERFORM_DEVELOPER_CHECKS)
    const bool correctAnswer = boost::algorithm::iequals( str, test );
    
    if( answer != correctAnswer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Got %i when should have got %i for label '%s' and string '%s'",
               int(answer), int(correctAnswer), test.c_str(), str.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( answer != correctAnswer )
#endif
    
    return answer;
  }
  
  // trim from start
  static inline std::string &ltrim(std::string &s)
  {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    
    if( s.size() )
    {
      const size_t pos = s.find_first_not_of( '\0' );
      if( pos != 0 && pos != string::npos )
        s.erase( s.begin(), s.begin() + pos );
      else if( pos == string::npos )
        s.clear();
    }
    
    return s;
  }
  
  // trim from end
  static inline std::string &rtrim(std::string &s)
  {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    
    //remove null terminating characters.  Boost doesnt do this, but is
    //  necassary when reading fixed width binary data.
    const size_t pos = s.find_last_not_of( '\0' );
    if( pos != string::npos && (pos+1) < s.size() )
      s.erase( s.begin() + pos + 1, s.end() );
    else if( pos == string::npos )
      s.clear();  //string is all '\0' characters
    
    return s;
  }
  
  // trim from both ends
  void trim( std::string &s )
  {
#if(PERFORM_DEVELOPER_CHECKS)
    string copystr = s;
    boost::algorithm::trim( copystr );
#endif
    
    ltrim( rtrim(s) );
    
#if(PERFORM_DEVELOPER_CHECKS)
    
    size_t pos = copystr.find_first_not_of( '\0' );
    if( pos != 0 && pos != string::npos )
      copystr.erase( copystr.begin(), copystr.begin() + pos );
    else if( pos == string::npos )
      copystr.clear();
    pos = copystr.find_last_not_of( '\0' );
    if( pos != string::npos && (pos+1) < s.size() )
      copystr.erase( copystr.begin() + pos + 1, copystr.end() );
    
    if( copystr != s )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Trimmed strings not equal expect: '%s' (len %i), got: '%s' (len %i, from boost)",
               s.c_str(), int(s.size()), copystr.c_str(), int(copystr.size()) );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }
#endif
  }//trim(...)
  
  std::string trim_copy( std::string str )
  {
    trim( str );
    return str;
  }//trim_copy(...)
  
  void split( std::vector<std::string> &resutls,
             const std::string &input, const char *delims )
  {
    resutls.clear();
    
    //The bellow implementation is not well tested, but in principle might be faster
    //  than boost::algorithm::split (should be tested more and them used)
    //It would aslso reduce final binary size by about 9 kb
    size_t prev_delim_end = 0;
    size_t delim_start = input.find_first_of( delims, prev_delim_end );
    
    while( delim_start != std::string::npos )
    {
      if( (delim_start-prev_delim_end) > 0 )
        resutls.push_back( input.substr(prev_delim_end,(delim_start-prev_delim_end)) );
      
      prev_delim_end = input.find_first_not_of( delims, delim_start + 1 );
      if( prev_delim_end != std::string::npos )
        delim_start = input.find_first_of( delims, prev_delim_end + 1 );
      else
        delim_start = std::string::npos;
    }//while( this_pos < input.size() )
    
    if( prev_delim_end < input.size() )
      resutls.push_back( input.substr(prev_delim_end) );

#if(PERFORM_DEVELOPER_CHECKS)
    vector<string> coorectResults;
    boost::algorithm::split( coorectResults, input, boost::is_any_of(delims),
                            boost::token_compress_on );
    while( coorectResults.size() && coorectResults[0].empty() )
      coorectResults.erase( coorectResults.begin() );
    while( coorectResults.size() && coorectResults[coorectResults.size()-1].empty() )
      coorectResults.erase( coorectResults.end()-1 );
    
    if( resutls != coorectResults )
    {
      stringstream errormsg;
      
      errormsg << "Error splitting '";
      for( size_t i = 0; i < input.size(); ++i )
        errormsg << "\\x" << std::hex << int(input[i]);
      errormsg << "' by seperators '";
      for( size_t i = 0; i < strlen(delims); ++i )
        errormsg << "\\x" << std::hex << int(delims[i]);
      errormsg << "'\n\texpected [";
      for( size_t i = 0; i < coorectResults.size(); ++i )
      {
        if( i )
         errormsg << ", ";
        errormsg << "'";
        for( size_t j = 0; j < coorectResults[i].size(); ++j )
          errormsg << "\\x" << std::hex << int(coorectResults[i][j]);
        errormsg << "'";
      }
      
      errormsg << "]\n\tGot [";
      
      for( size_t i = 0; i < resutls.size(); ++i )
      {
        if( i )
          errormsg << ", ";
        errormsg << "'";
        for( size_t j = 0; j < resutls[i].size(); ++j )
          errormsg << "\\x" << std::hex << int(resutls[i][j]);
        errormsg << "'";
      }
      errormsg << "]";
      
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg.str().c_str() );
    }//if( resutls != coorectResults )
#endif
  }//void split(...)
  
  void to_lower( string &input )
  {
#if(PERFORM_DEVELOPER_CHECKS)
    string strcopy = input;
    boost::algorithm::to_lower( strcopy );
#endif
    
    for( size_t i = 0; i < input.size(); ++i )
      input[i] = static_cast<char>( tolower(input[i]) );
    
#if(PERFORM_DEVELOPER_CHECKS)
    if( strcopy != input )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Failed to lowercase string.  Expected: '%s', got: '%s'",
               strcopy.c_str(), input.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }
#endif
  }//to_lower(...)
  
  std::string to_lower_copy( std::string input )
  {
    to_lower( input );
    return input;
  }
  
  void to_upper( string &input )
  {
#if(PERFORM_DEVELOPER_CHECKS)
    string strcopy = input;
    boost::algorithm::to_upper( strcopy );
#endif
    
    for( size_t i = 0; i < input.size(); ++i )
      input[i] = static_cast<char>( toupper(input[i]) );

#if(PERFORM_DEVELOPER_CHECKS)
    if( strcopy != input )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Failed to uppercase string.  Expected: '%s', got: '%s'",
               strcopy.c_str(), input.c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }
#endif
  }//void to_upper( string &input )
  
  void ireplace_all( std::string &input, const char *pattern, const char *replacement )
  {
    //This function does not handle UTF8!
    if( input.empty() )
      return;

#if(PERFORM_DEVELOPER_CHECKS)
    string strcopy = input, original = input;
    size_t reslen = strcopy.size() + 1;
    while( reslen != strcopy.size() )
    {
      reslen = strcopy.size();
      boost::algorithm::ireplace_all( strcopy, pattern, replacement );
    }
#endif
    
    const size_t paternlen = strlen(pattern);
    if( !paternlen )
      return;
    
    //If we are replacing "XX" with "X" in the sequence "XXX" we want the
    //  result to be "X", however we have to protect against the case where
    //  the replacement string contains the search string
    const bool replace_contains_pattern = icontains( replacement, pattern );
    
    const size_t replacment_len = strlen(replacement);
    bool found = true;
    const char *start = input.c_str(), *end = input.c_str() + input.size();
    
    while( found )
    {
      const char * const it = std::search( start, end, pattern, pattern+paternlen,
                                    char_iequal<char>(std::locale()) );
      found = (it != end);
      if( found )
      {
        size_t delstart = it - input.c_str();
        input.erase( delstart, paternlen );
        input.insert( delstart, replacement );
        start = input.c_str() + delstart + (replace_contains_pattern ? replacment_len : size_t(0));
        end = input.c_str() + input.size();
      }//if( found )
    }//while( found )

#if(PERFORM_DEVELOPER_CHECKS)
    if( strcopy != input )
    {
      stringstream msg;
      msg << "Failed to replace '";
      for( size_t i = 0; i < strlen(pattern); ++i )
        msg << "\\x" << hex << int(pattern[i]);
      msg << "' with '";
      for( size_t i = 0; i < strlen(replacement); ++i )
        msg << "\\x" << hex << int(replacement[i]);
      msg << "' in '";
      for( size_t i = 0; i < original.size(); ++i )
        msg << "\\x" << hex << int(original[i]);
      msg << "'.\n\tExpected: '";
      for( size_t i = 0; i < strcopy.size(); ++i )
        msg << "\\x" << hex << int(strcopy[i]);
      msg << "' (length " << strcopy.size() << ")\n\tGot:      '";
      for( size_t i = 0; i < input.size(); ++i )
        msg << "\\x" << hex << int(input[i]);
      msg << "' (length " << input.size() << ").";
      
      log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
    }
#endif
  }//void ireplace_all(...)
  
/*
  std::string ireplace_all_copy( const std::string &input,
                                const char *pattern, const char *replacement )
  {
    string result;
    
    if( input.empty() )
      return result;

    const size_t paternlen = strlen(pattern);
    if( !paternlen )
      return result;
    
#if(PERFORM_DEVELOPER_CHECKS)
    string strcopy = input, original = input;
    size_t reslen = strcopy.size() + 1;
    while( reslen != strcopy.size() )
    {
      reslen = strcopy.size();
      boost::algorithm::ireplace_all( strcopy, pattern, replacement );
    }
#endif

    
    
    const size_t replacment_len = strlen(replacement);
    
    
    vector<const char *> good_begin, good_end;
    
    size_t newstrlen = 0;
    const char *start = input.c_str();
    const char * const end = input.c_str() + input.size();
    
    while( start < end )
    {
      const char * const it = std::search( start, end, pattern, pattern+paternlen,
                                          char_iequal<char>(std::locale()) );
      
      good_begin.push_back( start );
      good_end.push_back( it );
      
      newstrlen += (it - start);
      start = it;
      
      if( it != end )
      {
        newstrlen += replacment_len;
        start += paternlen;
      }else
      {
        good_begin.push_back( end );
        good_end.push_back( end );
      }
    }//while( found )
    
    if( good_begin.empty() )
      return input;
    
    result.reserve( newstrlen + 1 );  //+1, because I'm not sure if reserve takes into account '\0'
    
    for( size_t i = 0; i < good_begin.size(); ++i )
    {
      result.insert( result.end(), good_begin[i], good_end[i] );
      if( good_begin[i] != good_end[i] )
        result.insert( result.end(), replacement, replacement+replacment_len );
    }
    
#if(PERFORM_DEVELOPER_CHECKS)
    if( strcopy != result )
    {
      stringstream msg;
      msg << "Failed to replace '";
      for( size_t i = 0; i < strlen(pattern); ++i )
        msg << "\\x" << hex << int(pattern[i]);
      msg << "' with '";
      for( size_t i = 0; i < strlen(replacement); ++i )
        msg << "\\x" << hex << int(replacement[i]);
      msg << "' in '";
      for( size_t i = 0; i < original.size(); ++i )
        msg << "\\x" << hex << int(original[i]);
      msg << "'.\n\tExpected: '";
      for( size_t i = 0; i < strcopy.size(); ++i )
        msg << "\\x" << hex << int(strcopy[i]);
      msg << "' (length " << strcopy.size() << ")\n\tGot:      '";
      for( size_t i = 0; i < result.size(); ++i )
        msg << "\\x" << hex << int(result[i]);
      msg << "' (length " << result.size() << ").";
      
      log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
    }
#endif
    
    return result;
  }//ireplace_all_copy(...)
*/
  
  template<typename IterType>
  size_t utf8_iterate( IterType& it, const IterType& last )
  {
    if(it == last)
      return 0;
    
    unsigned char c;
    size_t res = 1;
    for(++it; last != it; ++it, ++res)
    {
      c = *it;
      
      //if highest value digit is not set, or if highest two digits are set
      //If the most significant bit isnt set, then its an ascii character.
      //  If the two most significant bits are set, then its the start of
      //  a UTF8 character; middle/end UTF8 bytes have most significant bit set,
      //  and second most significant bit not set
      //see: http://www.cprogramming.com/tutorial/unicode.html
      // 0x80 --> 10000000
      // 0xC0 --> 11000000
      if( !(c & 0x80) || ((c & 0xC0) == 0xC0))
        break;
    }
    
    return res;
  }//size_t utf8_str_iterator( IterType& it, const IterType& last )
  
  
  size_t utf8_iterate( const char * &it )
  {
    //Assumes null-terminated string
    if( !(*it) )
      return 0;
    
    unsigned char c;
    size_t res = 0;
    for( ; *it; ++it, ++res)
    {
      c = *it;
      if( !(c & 0x80) || ((c & 0xC0) == 0xC0))
        break;
    }
    
    return res;
  }//size_t utf8_str_iterator( IterType& it, const IterType& last )
  
  
  size_t utf8_str_len( const char * const str, size_t str_size_bytes )
  {
    size_t len = 0;
    
    if( !str_size_bytes )
    {
      for( const char *ptr = str; *ptr; utf8_iterate(ptr) )
        ++len;
    }else
    {
      for( const char *ptr = str, * const end = str + str_size_bytes;
           ptr != end; utf8_iterate(ptr, end) )
        ++len;
    }
    
    return len;
  }//size_t utf8_str_len( const char * const str, size_t str_size_bytes )
  
  
  void utf8_limit_str_size( std::string &str, const size_t max_bytes )
  {
    const size_t pos = utf8_str_size_limit( str.c_str(), str.size(), max_bytes );
    str = str.substr( 0, pos );
  }
  
  
  size_t utf8_str_size_limit( const char *str,
                              size_t len, const size_t max_bytes )
  {
    if( !len )
      len = strlen( str );
    
    if( !len || !max_bytes )
      return 0;
  
    if( len < max_bytes )
      return len;
    
    const char *iter = str + max_bytes - 1;
    for( ; iter != str; --iter )
    {
      if( ((*iter) & 0xC0) == 0xC0 )
        break;
    }
    
    return iter - str;
  }
  

  std::string to_common_string( const boost::posix_time::ptime &t )
  {
    if( t.is_special() )
      return "not-a-date-time";
    
    const int year = static_cast<int>( t.date().year() );
    const int day = static_cast<int>( t.date().day() );
    int hour = static_cast<int>( t.time_of_day().hours() );
    const int mins = static_cast<int>( t.time_of_day().minutes() );
    const int secs = static_cast<int>( t.time_of_day().seconds() );
    
    bool is_pm = (hour >= 12);
    if( is_pm )
      hour -= 12;
    if( hour == 0 )
      hour = 12;
    
    const char *month = "";
    switch( t.date().month() )
    {
      case 1: month = "Jan"; break;
      case 2: month = "Feb"; break;
      case 3: month = "Mar"; break;
      case 4: month = "Apr"; break;
      case 5: month = "May"; break;
      case 6: month = "Jun"; break;
      case 7: month = "Jul"; break;
      case 8: month = "Aug"; break;
      case 9: month = "Sep"; break;
      case 10: month = "Oct"; break;
      case 11: month = "Nov"; break;
      case 12: month = "Dec"; break;
    }
    
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%i-%s-%04i %02i:%02i:%02i %s",
              day, month, year, hour, mins, secs, (is_pm ? "PM" : "AM") );
    
    return buffer;
  }//to_common_string

  
  std::string print_to_iso_str( const boost::posix_time::ptime &t,
                               const bool extended )
  {
    //#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
    if( t.is_special() )
      return "not-a-date-time";
    //should try +-inf as well
    
    const int year = static_cast<int>( t.date().year() );
    const int month = static_cast<int>( t.date().month() );
    const int day = static_cast<int>( t.date().day() );
    const int hour = static_cast<int>( t.time_of_day().hours() );
    const int mins = static_cast<int>( t.time_of_day().minutes() );
    const int secs = static_cast<int>( t.time_of_day().seconds() );
    double frac = t.time_of_day().fractional_seconds()
    / double(boost::posix_time::time_duration::ticks_per_second());
    
    char buffer[256];
    if( extended ) //"2014-04-14T14:12:01.621543"
      snprintf( buffer, sizeof(buffer),
               "%i-%.2i-%.2iT%.2i:%.2i:%09.6f",
               year, month, day, hour, mins, (secs+frac) );
    else           //"20140414T141201.621543"
      snprintf( buffer, sizeof(buffer),
               "%i%.2i%.2iT%.2i%.2i%09.6f",
               year, month, day, hour, mins, (secs+frac) );

#if(PERFORM_DEVELOPER_CHECKS)
    string tmpval = buffer;
    if( tmpval.find(".") == string::npos )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Expected there to be a '.' in iso date time: '%s'", buffer );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }
#endif

    const char point = '.';
//    const char point = std::use_facet< std::numpunct<char> >(std::cout.getloc()).decimal_point();
    
    //Get rid of trailing zeros
    size_t result_len = strlen(buffer) - 1;
    while( result_len > 1 && buffer[result_len]=='0' )
      buffer[result_len--] = '\0';
    
    if( result_len > 1 && buffer[result_len]==point )
      buffer[result_len--] = '\0';
    
#if(PERFORM_DEVELOPER_CHECKS)
    string correctAnswer = extended
                            ? boost::posix_time::to_iso_extended_string( t )
                            : boost::posix_time::to_iso_string( t );
    
    if( correctAnswer.find(".") != string::npos )
    {
      result_len = correctAnswer.size() - 1;
      while( result_len > 1 && correctAnswer[result_len]=='0' )
        correctAnswer = correctAnswer.substr( 0, result_len-- );
    
      if( result_len > 1 && buffer[result_len]==point )
        correctAnswer = correctAnswer.substr( 0, result_len-- );
    }//if( correctAnswer.find(".") != string::npos )
    
    if( correctAnswer != buffer )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
               "Failed to format date correctly for %sextended iso format. Expected: '%s', got: '%s'",
               (extended ? "" : "non-"),
               correctAnswer.c_str(), buffer );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }
#endif
    
    
    return buffer;
  }//std::string to_extended_iso_string( const boost::posix_time::ptime &t )
  
  std::string to_extended_iso_string( const boost::posix_time::ptime &t )
  {
    return print_to_iso_str( t, true );
  }
  
  std::string to_iso_string( const boost::posix_time::ptime &t )
  {
    return print_to_iso_str( t, false );
  }//std::string to_iso_string( const boost::posix_time::ptime &t )
  

  boost::posix_time::ptime time_from_string( const char *time_string )
  {
#define CREATE_datetimes_TEST_FILE 0
    
#if( CREATE_datetimes_TEST_FILE )
    static std::mutex datetimes_file_mutex;
    static int ntimes_called = 0;
    string inputstr = time_string;
    UtilityFunctions::ireplace_all(inputstr, ",", "");
    if( inputstr.size() && inputstr[0]=='#' )
      inputstr = inputstr.substr(1);
    
//#ifndef _WIN32
    auto result = time_from_string_strptime( time_string, MiddleEndianFirst );
//#else
//    auto result = time_from_string_boost( time_string );
//#endif
    
    if( inputstr.empty() )
      return result;
    
    std::lock_guard<std::mutex> lock( datetimes_file_mutex );
    ++ntimes_called;
    if( (ntimes_called % 1000) == 0 )
      cerr << "Warning, time_from_string() is creating datetimes.txt test file" << endl;
    
    ofstream output( "datetimes.txt", ios::out | ios::app );
    
    if( output )
      output << inputstr << "," << UtilityFunctions::to_extended_iso_string(result) << "\r\n";
    else
      cerr << "Failed to open datetimes.txt for appending" << endl;
    
    return result;
#else  //CREATE_datetimes_TEST_FILE

//#ifndef _WIN32
    return time_from_string_strptime( time_string, MiddleEndianFirst );
//#else
//	return time_from_string_boost( time_string );
//#endif
#endif
  }
  

// #if( !(defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64)) )

#if( ANDROID )
  static_assert( 0, "Not sure about std::timegm existing for android c++11, you need to check, and then remove all traces of this hack" )
//inline time_t timegm(std::tm const *t){ return boost::chrono::detail::internal_timegm(t); }
#endif

#if( defined(WIN32) )
#define timegm _mkgmtime
#endif

bool strptime_wrapper( const char *s, const char *f, struct tm *t )
{
  //For the testTimeFromString unit test on my mac, the native take strptime
  //  takes 302188us to run, vs the c++11 version taking 4113835us.
  //  ~10 times slower, so preffer native strptime where available.
#ifdef _WIN32
  #define HAS_NATIVE_STRPTIME 0
  //static_assert( 0, "Have not tested c++11 implementation of strptime on Windows yet - you must do this (run testTimeFromString unit test at least)" );
  //After veriying all is okay, get rid of time_from_string_boost(...) and all traces of it
#else
  #define HAS_NATIVE_STRPTIME 1
#endif
  
#if( HAS_NATIVE_STRPTIME )
  return (strptime(s,f,t) != nullptr);
#else
  //*t = std::tm();  //should already be done
  //memset( t, 0, sizeof(*t) );  //Without this some tm fields dont make sense for some formats...
  
  std::istringstream input(s);
  input.imbue(std::locale(setlocale(LC_ALL, nullptr)));
  input >> std::get_time(t, f);
  if( input.fail() )
    return false;
  
  //cout << "Format '" << f << "' parsed '" << s << "' to " << std::put_time(t, "%c") << endl;
  /*
   cout
   << "seconds after the minute [0-60]: " << t.tm_sec << endl
   << "minutes after the hour [0-59]: " << t.tm_min << endl
   << "hours since midnight [0-23]: " << t.tm_hour << endl
   << "day of the month [1-31]: " << t.tm_mday << endl
   << "months since January [0-11]: " << t.tm_mon << endl
   << "years since 1900: " << t.tm_year << endl
   << "days since Sunday [0-6]: " << t.tm_wday << endl
   << "days since January 1 [0-365]: " << t.tm_yday << endl
   << "Daylight Savings Time flag: " << t.tm_isdst << endl
   << "offset from UTC in seconds: " << t.tm_gmtoff << endl
   //<< "timezone abbreviation: " << (t.tm_zone ? (const char *)t.tm_zone : "null")
   << endl;
   */
  
  return true;
#endif
}//char *strptime_wrapper(...)
  
  boost::posix_time::ptime time_from_string_strptime( std::string time_string,
                                                     const DateParseEndianType endian )
  {
    UtilityFunctions::to_upper( time_string );  //make sure strings like "2009-11-10t14:47:12" will work (some file parsers convert to lower case)
    UtilityFunctions::ireplace_all( time_string, "  ", " " );
    UtilityFunctions::ireplace_all( time_string, "_T", "T" );  //Smiths NaI HPRDS: 2009-11-10_T14:47:12Z
    UtilityFunctions::trim( time_string );
    
    //strptime(....) cant handle GMT offset (ex '2014-10-24T08:05:43-04:00')
    //  so we will do this manually.
    boost::posix_time::time_duration gmtoffset(0,0,0);
    const size_t offsetcolon = time_string.find_first_of( ':' );
    if( offsetcolon != string::npos )
    {
      const size_t signpos = time_string.find_first_of( "-+", offsetcolon+1 );
      if( signpos != string::npos )
      {
        //Note: the + and - symbols are in the bellow for dates like:
        //  '3-07-31T14:25:30-03:-59'
        const size_t endoffset = time_string.find_first_not_of( ":0123456789+-", signpos+1 );
        const string offset = endoffset==string::npos ? time_string.substr(signpos)
        : time_string.substr(signpos,endoffset-signpos-1);
        string normal = time_string.substr(0,signpos);
        if( endoffset != string::npos )
          normal += time_string.substr( endoffset );
        
        //normal will look like "2014-10-24T08:05:43"
        //offset will look like "-04:00"
        try
        {
          gmtoffset = boost::posix_time::duration_from_string(offset);
        }catch( std::exception &e )
        {
          cerr << "Failed to convert '" << offset << "' to time duration, from string '" << time_string << "'" << endl;
        }
        
        
        //Should also make sure gmtoffset is a reasonalbe value.
        
        //      cout << "offset for '" << time_string << "' is '" << offset << "' with normal time '" <<  normal<< "'" << endl;
        time_string = normal;
      }//if( signpos != string::npos )
    }//if( offsetcolon != string::npos )
    
    
    //strptime(....) cant handle fractions of second (because tm doesnt either)
    //  so we will manually convert fractions of seconds.
    boost::posix_time::time_duration fraction(0,0,0);
    //int fraction_nano_sec = 0;
    
    const size_t fraccolon = time_string.find_last_of( ':' );
    if( fraccolon != string::npos )
    {
      const size_t period = time_string.find( '.', fraccolon+1 );
      if( period != string::npos )
      {
        const size_t last = time_string.find_first_not_of( "0123456789", period+1 );
        string fracstr = ((last!=string::npos)
                          ? time_string.substr(period+1,last-period-1)
                          : time_string.substr(period+1));
        
        //Assume microsecond resolution at the best (note
        //  boost::posix_time::nanosecond isnt available on my OS X install)
        const size_t ndigits = 9;
        const size_t invfrac = 1E9;
        const size_t nticks = boost::posix_time::time_duration::ticks_per_second();
        
        if( fracstr.size() < ndigits )
          fracstr.resize( ndigits, '0' );
        else if( fracstr.size() > ndigits )
          fracstr.insert( ndigits, "." );
        
        int numres = 0;  //using int will get rounding wrong
        if( (stringstream(fracstr) >> numres) )
        {
          //fraction_nano_sec = numres;
          fraction = boost::posix_time::time_duration(0,0,0, numres*nticks/invfrac);
        }else
          cerr << "Failed to convert fraction '" << fracstr << "' to double" << endl;
        
        string normal = time_string.substr(0,period);
        if( last != string::npos )
          normal += time_string.substr(last);
        time_string = normal;
      }//if( period != string::npos )
    }//if( fraccolon != string::npos )
    
    //With years like 2070, 2096, and such, strptime seems to fail badly, so we
    //  will fix them up a bit.
    //Assumes the first time you'll get four numbers in a row, it will be the year
    bool add100Years = false;
    if( time_string.size() > 5 )
    {
      for( size_t i = 0; i < time_string.size()-4; ++i )
      {
        if( isdigit(time_string[i]) && isdigit(time_string[i+1])
           && isdigit(time_string[i+2]) && isdigit(time_string[i+3]) )
        {
          int value;
          if( stringstream(time_string.substr(i,4)) >> value )
          {
            //XXX - I havent yet determined what year this issue starts at
            if( value > 2030 && value < 2100 )
            {
              char buffer[8];
              snprintf( buffer, sizeof(buffer), "%i", value - 100 );
              time_string[i] = buffer[0];
              time_string[i+1] = buffer[1];
              time_string[i+2] = buffer[2];
              time_string[i+3] = buffer[3];
              add100Years = true;
            }
          }//if( stringstream(time_string.substr(i,4)) >> value )
          
          break;
        }//if( four numbers in a row )
      }//for( size_t i = 0; i < time_string.size(); ++i )
    }//if( time_string.size() > 5 )
    
    
    
    //middle: whether to try a middle endian (date start with month number)
    //  date decoding first, or alternetavely little endian (date start with day
    //  number).  Both endians will be tried, this just selects which one first.
    const bool middle = (endian == MiddleEndianFirst);
    const bool only = (endian == MiddleEndianOnly || endian == LittleEndianOnly);
    
    //Should probably just do a string replace of 'T' with a space to almost cut
    //  in half the number of string formats to try.
    //  Also, should go through list of formats listed in http://www.partow.net/programming/datetime/index.html
    //  and make sure they all parse.
    //
    const char * const formats[] =
    {
      "%d-%b-%y%n%r", //'15-MAY-14 08:30:44 PM'  (disambiguos: May 15 2014)
      "%Y-%m-%dT%H:%M:%SZ", //2010-01-15T23:21:15Z
      "%d-%b-%Y%n%r",       //1-Oct-2004 12:34:42 AM  //"%r" ~ "%I:%M:%S %p"
      (middle ? "%m/%d/%Y%n%r" : "%m/%d/%Y%n%r"),             //1/18/2008 2:54:44 PM
      (only ? "" : (middle ? "%m/%d/%Y%n%r" : "%m/%d/%Y%n%r")),
      (middle ? "%m/%d/%Y%n%H:%M:%S" : "%m/%d/%Y%n%H:%M:%S"), //08/05/2014 14:51:09
      (only ? "" : (middle ? "%m/%d/%Y%n%H:%M:%S" : "%m/%d/%Y%n%H:%M:%S")),
      (middle ? "%m-%d-%Y%n%H:%M:%S" : "%d-%m-%Y%n%H:%M:%S"), //14-10-2014 16:15:52
      (only ? "" : (middle ? "%d-%m-%Y%n%H:%M:%S" : "%m-%d-%Y%n%H:%M:%S" )),
      "%d-%b-%y%n%H:%M:%S", //16-MAR-06 13:31:02, or "12-SEP-12 11:23:30"
      "%d-%b-%Y%n%H:%M:%S", //31-Aug-2005 12:38:04,
      "%d %b %Y%n%H:%M:%S", //31 Aug 2005 12:38:04
      "%d-%b-%YT%H:%M:%S", //21-Feb-2015T07:35:32-05:00 //not working
      "%Y-%m-%dT%H:%M:%S",  //2014-03-27T08:58:22
      "%d.%m.%Y%n%H:%M:%S",
//      (middle ? "%m.%d.%Y%n%H:%M:%S" : "%d.%m.%Y%n%H:%M:%S"), //26.05.2010 02:53:49
      "%d.%m.%Y%n%H:%M:%S",
//      (only ? "" : (middle ? "%d.%m.%Y%n%H:%M:%S" : "%m.%d.%Y%n%H:%M:%S")),
      "%b. %d %Y%n%H:%M:%S",//May. 21 2013  07:06:42
      "%d.%m.%y%n%H:%M:%S",  //28.02.13 13:42:47
//      (middle ? "%m.%d.%y%n%H:%M:%S" : "%d.%m.%y%n%H:%M:%S"),  //28.02.13 13:42:47
      "%d.%m.%y%n%H:%M:%S",
//      (only ? "" : (middle ? "%d.%m.%y%n%H:%M:%S" : "%m.%d.%y%n%H:%M:%S")),
      "%d-%b-%YT%H:%M:%S%nZ",//9-Sep-2014T20:29:21 Z
      "%Y.%m.%d%n%H:%M:%S", //2012.07.28 16:48:02
      "%Y.%m.%dT%H:%M:%S", //2012.07.28T16:48:02
      "%Y-%m-%d%n%H:%M:%S", //2013-12-12 15:28:12
      "%Y-%m-%dT%H:%M:%S", //2013-12-12T15:28:12+00:00
      "%d.%b.%Y %H:%M:%S",//01.Nov.2010 21:43:35
      "%d.%b.%YT%H:%M:%S",//01.Nov.2010T21:43:35
      "%Y%m%dT%H:%M:%S",  //20100115T23:21:15
      "%Y%m%d %H:%M:%S", //20100115T23:21:15
      "%Y-%b-%d %H:%M:%S" //2017-Jul-07 09:16:37
    };
    
    
    const size_t nformats = sizeof(formats) / sizeof(formats[0]);
    
    const char *timestr = time_string.c_str();
    
    for( size_t i = 0; i < nformats; ++i )
    {
      struct tm t = std::tm();
      
      if( formats[i][0] && strptime_wrapper( timestr, formats[i], &t ) )
      {
        //if( add100Years )
          //t.tm_year += 100;
        //std::chrono::time_point tp = system_clock::from_time_t( std::mktime(&t) ) + std::chrono::nanoseconds(fraction_nano_sec);
        
        return boost::posix_time::from_time_t( timegm(&t) )
        + fraction
        + boost::gregorian::years( add100Years ? 100 : 0 )
        /*+ gmtoffset*/;  //ignore offset since we want time in local persons zone
      }
      //      return boost::posix_time::from_time_t( mktime(&tm) - timezone + 3600*daylight ) + fraction + gmtoffset;
    }//for( size_t i = 0; i < nformats; ++i )
    
    return boost::posix_time::ptime();
  }//boost::posix_time::ptime time_from_string_strptime( std::string time_string )
//#endif  //#ifndef _WIN32
  
#if(PERFORM_DEVELOPER_CHECKS)
  boost::posix_time::ptime time_from_string_boost( const char *time_string )
  {
    namespace pt = boost::posix_time;
    
    //  time_string = "2014-5-7T3:37:30";
    
    //  static std::mutex file_mutex;
    //  {
    //    std::lock_guard<std::mutex> lock( file_mutex );
    //    ofstream file( "datetimes.txt", ios::app );
    //    file << time_string << endl;
    //  }
    
    pt::ptime answer;
    
    if( !time_string )
      return answer;
    

    string edited_str;
    
    edited_str = time_string;
    UtilityFunctions::to_upper( edited_str );
    time_string = edited_str.c_str();
    
    
    //If the input has an AM or PM in it, then we have to manually manipulate the
    //  string so boost can canvert it
    const bool hasAm = icontains(time_string,"AM");
    const bool add_12_hours = icontains(time_string,"PM");
    if( hasAm || add_12_hours )
    {
      edited_str = time_string;
      if( hasAm )
        ireplace_all( edited_str, "AM", "" );
      else
        ireplace_all( edited_str, "PM", "" );
      
      trim( edited_str );
      
      //now make sure hour has two digits, or else boost will fail to convert
      //  due to %I flag not yet supported for input.
      //  eg. "9/29/2008 9:47:40" needs to become "9/29/2008 09:47:40"
      string::size_type pos = edited_str.find( ':' );
      if( pos > 2 && pos != string::npos
         && isdigit(edited_str[pos-1])
         && !isdigit(edited_str[pos-2]) )
      {
        edited_str = edited_str.substr(0,pos-1) + "0" + edited_str.substr(pos-1);
      }//if( pos > 2 && pos != string::npos )
      
      time_string = edited_str.c_str();
    }//if( hasAm || add_12_hours )
    
    
    //Check for: Smiths NaI HPRDS: 2009-11-10_T14:47:12Z
    const char *underscore_pos = strstr( time_string, "_T" );
    if( underscore_pos )
    {
      edited_str = time_string;
      UtilityFunctions::ireplace_all( edited_str, "_T", "T" );
      time_string = edited_str.c_str();
    }
    
    
    {
      //make dataes like '2012-11-23T9:28:36' into '2012-11-23T09:28:36'
      string tmpstr = time_string;
      size_t pos = 0;
      while( (pos = tmpstr.find_first_of( ":",pos+1)) != string::npos )
      {
        if( pos > 2 && !isdigit(tmpstr[pos-2]) )
        {
          edited_str = tmpstr = tmpstr.substr(0,pos-1) + "0" + tmpstr.substr(pos-1);
          time_string = edited_str.c_str();
        }
      }
      
      //make dates like '1/7/2008 08:05:44' into '01/07/2008 08:05:44'
      pos = 0;
      while( (pos = tmpstr.find_first_of( "/",pos+1)) != string::npos )
      {
        //Using alphanum instead isdigit because of dates like
        //  "15-Jun-2012 12:47:51"
        if( pos==1 || (pos > 1 && !isalnum(tmpstr[pos-2])) )
        {
          edited_str = tmpstr = tmpstr.substr(0,pos-1) + "0" + tmpstr.substr(pos-1);
          time_string = edited_str.c_str();
        }
      }
      
      //dates like 2014-5-7T04:40:13 to 2014-05-07T04:40:13
      size_t seppos = tmpstr.find_first_of( " T" );
      if( seppos != string::npos && seppos > 2 )
      {
        pos = 0;
        while( 1 )
        {
          pos = tmpstr.find_first_of( "-",pos+1);
          
          if( pos > seppos || pos == string::npos )
          {
            if( isdigit(tmpstr[seppos-1]) && tmpstr[seppos-2]=='-' )
            {
              edited_str = tmpstr = tmpstr.substr(0,seppos-1) + "0" + tmpstr.substr(seppos-1);
              time_string = edited_str.c_str();
              ++seppos;
            }
            
            break;
          }
          
          //Using alphanum instead isdigit because of dates like
          //  "15-Jun-2012 12:47:51"
          if( pos==1 || (pos > 1 && !isalnum(tmpstr[pos-2])) )
          {
            edited_str = tmpstr = tmpstr.substr(0,pos-1) + "0" + tmpstr.substr(pos-1);
            time_string = edited_str.c_str();
            ++seppos;
          }
        }
      }//if( seppos != string::npos )
      
    }
    
    
    //static variable initialization is not thread-safe in C++03 (but is in C++11
    //  I think).  So we have to actually take a mutex (because C++03 doesnt have
    //  atomics).  I'm sure I could do this in a more efficient manner (e.g.
    //  lockless, or sperately for C++11/C++03), but whatever for now since this
    //  function probably isnt a bottlneck anywhere.
    static std::mutex format_mutex;
    static std::unique_ptr< const std::vector<locale> > format_locals;
    
    const char *const formats[] =
    {
      "%Y-%b-%dT%H:%M:%S%FZ",
      "%Y-%m-%dT%H:%M:%S%F%Q",   //Default for Cambio I believe "2010-02-24T00:08:24+00:00"
      "%Y%m%dT%H%M%S%F%q",       //ISO Format: "20051015T131211-0700"
      "%Y-%m-%d %H:%M:%S%F%Q",   //Extended ISO format: "2005-10-15 13:12:11-07:00"
      "%d-%b-%Y %H:%M:%S%F",     //GADRAS format?
      "%Y-%b-%d %H:%M:%S%F %z",  //Default ptime output "2005-Oct-15 13:12:11 MST"
      "%Y-%b-%d %H:%M:%S%F %ZP", //Default ptime input "2005-Oct-15 13:12:11 MST-07"
      "%Y-%m-%dT%H:%M:%S%FZ",    //ORTEC IDM "2010-03-17T17:33:19.409Z"
      "%Y-%m-%d %H:%M:%S",
      "%Y/%m/%d %H:%M:%S",
      "%d.%m.%Y %H:%M:%S",
      "%d.%m.%y %H:%M:%S",       //used by identifinder (ascii spc files): "28.08.2012 16:12:26"
      "%m-%d-%y %H%M",           //used by DetectiveEx in its comments section (binary spc files) "09-06-12 0415"
      "%d-%b-%y %H:%M:%S",       //close to what is used by DetectiveEx (binary spc file) 06-SEP-12 04:17:27
      "%m/%d/%Y %H:%M:%S",       //used in some SPE files (like ORTEC IDM), ex '03/16/2010 12:01:04'
      "%m/%d/%Y %H:%M:%S%F",
      "%Y-%m-%d",
      "%Y-%b-%d %H:%M:%S", //PCF format 2017-Jul-07 09:43:07
      "%b. %d %Y  %H:%M:%S"     //Used in GR135 TXT files: "Oct. 09 2013  13:08:29"
    };//static const locale formats[] = {...}
    
    
    {//begin codeblock to check, and if necessary, initialize format_locals
      std::lock_guard<std::mutex> lock( format_mutex );
      if( !format_locals )
      {
        const size_t n_formats = sizeof(formats)/sizeof(formats[0]);
        std::vector<locale> *localsvec = new std::vector<locale>( n_formats );
        format_locals.reset( localsvec );
        
        for( size_t i = 0; i < n_formats; ++i )
        {
          pt::time_input_facet *facet = new pt::time_input_facet( formats[i] );
          (*localsvec)[i] = locale(locale::classic(), facet);
        }
      }//if( !format_locals )
    }//end codeblock to check, and if necessary, initialize format_locals
    
    const locale *dt_locals = &((*format_locals)[0]);
    const size_t ndt_locals = format_locals->size();
    
    for( size_t i = 0; i < ndt_locals; ++i )
    {
      stringstream ss( time_string );
      ss.imbue( dt_locals[i] );
      
      if( (ss >> answer) )
        break;
    }//for( loop over potential date formats )
    
    //  cout << "time_string-->" << time_string << endl;
    
    if( answer.is_special() && edited_str.size() )
      cerr << "Failed to convert '" << time_string << "' to datetime" << endl;
    else if( hasAm )
    {
      //12:30 noon is 12:30 pm, so if its 12:30 am, we need to subtract 12 hours
      if( answer.time_of_day().hours() == 12 )
        answer -= pt::hours(12);
    }else if( add_12_hours )
    {
      if( answer.time_of_day().hours() != 12 )
        answer += pt::hours(12);
    }
    
    return answer;
  }//boost::posix_time::ptime time_from_string( const char *str )
#endif
  
std::string temp_dir()
{
#if( ANDROID )
  {
    const char *val = std::getenv("TMPDIR");
    if( !val )
    {
      cerr << "Warning, unable to get \"TMPDIR\" environment variable;"
           << "returning: \"/data/local/tmp/\"" << endl;
      return "/data/local/tmp/";
    }
  
    return val;
  }
#endif
  
#if( !SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB && BOOST_VERSION >= 104700)
  try
  {
    return boost::filesystem::temp_directory_path().generic_string();
  }catch( std::exception &e )
  {
    cerr << "Warning, unable to get a temporary directory...: " << e.what()
         << "\nReturning '/tmp'" << endl;
  }//try/catch
#endif

#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  //Completely un-tested
  const DWORD len = GetTempPathW( 0, NULL );
  vector<TCHAR> buf( len );
  
  if( !len || !GetTempPath( len, &buf[0] ) )
  {
    const char *val = 0;
    (val = std::getenv("temp" )) ||
    (val = std::getenv("TEMP"));
    
    if( val )
      return val;
    
#if(PERFORM_DEVELOPER_CHECKS)
    log_developer_error( BOOST_CURRENT_FUNCTION, "Couldnt find temp path on Windows" );
#endif
    return "C:\\Temp";
  }
  
#if( defined(UNICODE) || defined(_UNICODE) )
  // TCHAR type is wchar_t
  const size_t utf8len = wcstombs( NULL, &buf[0], len );
  if( utf8len == 0 || static_cast<size_t>( int(-1) ) == utf8len )
  {
#if(PERFORM_DEVELOPER_CHECKS)
    log_developer_error( BOOST_CURRENT_FUNCTION, "Error converting temp path to UTF8 on Windows" );
#endif
    return "C:\\Temp";
  }
  vector<char> thepath( utf8len );
  wcstombs( &thepath[0], &buf[0], len );
  return string( thepath.begin(), thepath.end() );
#else
  // TCHAR type is char
  return string(buf.begin(), buf.end());
#endif
  
#else
  
  const char *val = NULL;
  (val = std::getenv("TMPDIR" )) ||
  (val = std::getenv("TMP"    )) ||
  (val = std::getenv("TEMP"   )) ||
  (val = std::getenv("TEMPDIR"));

  if( val && UtilityFunctions::is_directory(val) )
    return val;

  return "/tmp";
#endif
}//std::string temp_dir()

  
bool remove_file( const std::string &name )
{
#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
  return (0 == unlink(name.c_str()) );
#else
  try{ boost::filesystem::remove( name ); } catch(...){ return false; }
  return true;
#endif
}//bool remove_file( const std::string &name )

  
bool is_file( const std::string &name )
{
#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
//  struct stat fileInfo;
//  const int status = stat( name.c_str(), &fileInfo );
//  if( status != 0 )
//    return false;
//  return S_ISREG(fileinfo.st_mode);
  ifstream file( name.c_str() );
  return file.good();
#else
  bool isfile = false;
  try
  {
    isfile = (boost::filesystem::exists( name )
              && !boost::filesystem::is_directory( name ));
  }catch(...){}
  
  return isfile;
#endif
}//bool is_file( const std::string &name )
  

bool is_directory( const std::string &name )
{
#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
  struct stat statbuf;
  stat( name.c_str(), &statbuf);
  return S_ISDIR(statbuf.st_mode);
#else
  try{ return boost::filesystem::is_directory( name ); }catch(...){}
  return false;
#endif
}//bool is_directory( const std::string &name )

  
int create_directory( const std::string &name )
{
  if( is_directory(name) )
    return -1;
  
  int nError = 0;
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  nError = _mkdir(name.c_str()); // can be used on Windows
#else
  mode_t nMode = 0733; // UNIX style permissions
  nError = mkdir(name.c_str(),nMode); // can be used on non-Windows
#endif
  if (nError != 0) {
    // handle your error here
    return 0;
  }
  
  return 1;
}
  
bool can_rw_in_directory( const std::string &name )
{
  if( !is_directory(name) )
    return false;
    
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  const int can_access = _access( name.c_str(), 06 );
#else
  const int can_access = access( name.c_str(), R_OK | W_OK | X_OK );  
#endif
    
  return (can_access == 0);
}//bool can_rw_in_directory( const std::string &name )
  
  
std::string append_path( const std::string &base, const std::string &name )
{
#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  if( base.size() && (base[base.size()-1]=='\\'||base[base.size()-1]=='/') )
    return base + name;
  if( name.size() && (name[0]=='\\'||name[0]=='/') )
    return base + name;
  return base + '\\' + name;
#else
  if( base.size() && base[base.size()-1]=='/' )
    return base + name;
  if( name.size() && name[0]=='/' )
    return base + name;
  return base + '/' + name;
#endif
#else
  boost::filesystem::path p(base);
  p /= name;
#if( BOOST_VERSION < 106501 )
  return p.make_preferred().string<string>();
#else
  return p.make_preferred().lexically_normal().string<string>();
#endif
#endif
}//std::string append_path( const std::string &base, const std::string &name )
  

std::string filename( const std::string &path_and_name )
{
  
#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
  if( path_and_name.empty() )
    return path_and_name;
  
  string::size_type lastslash = path_and_name.find_last_of( "/\\" );
  if( lastslash == (path_and_name.size()-1) )
  {
    
  }
  
#error "UtilityFunctions::filename not implemented for SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB"
#else
  return boost::filesystem::path(path_and_name).filename().string<string>();
#endif
}//std::string filename( const std::string &path_and_name )
  
  
std::string parent_path( const std::string &path )
{
#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
#error "UtilityFunctions::parent_path not implemented for SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB"
#else
  return boost::filesystem::path(path).parent_path().string<string>();
#endif
}//std::string parent_path( const std::string &path )

  
std::string file_extension( const std::string &path )
{
  const string fn = filename( path );
  const size_t pos = fn.find_last_of( '.' );
  if( pos == string::npos )
    return "";
  return fn.substr(pos);
}
  
#if defined(WIN32) || defined(WIN64)
  // Copied from linux libc sys/stat.h:
  #define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
  #define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#endif
  
size_t file_size( const std::string &path )
{
  struct stat st;
  if( stat(path.c_str(), &st) < 0 )
    return 0;
  
  if( S_ISDIR(st.st_mode) )
    return 0;
  
  return st.st_size;
}
  
  
#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
std::string temp_file_name( std::string base, std::string directory )
{
  //see http://msdn.microsoft.com/en-us/library/aa363875%28VS.85%29.aspx
  // for windows
  //Or just grab from unique_path.cpp in boost.
  //Or, it looks like tmpnam should be definced in #include <cstdio> and in the std namespace
  
  if( directory.empty() || !is_directory(directory) )
    directory = UtilityFunctions::temp_dir();
    
  char buffer[L_tmpnam+1];
  tmpnam( buffer+1 );
  buffer[0] = '_';
  
  base += buffer;
  
  return append_path( directory, base );
}
#else
std::string temp_file_name( std::string bases, std::string temppaths )
{
  using boost::filesystem::path;
  using boost::filesystem::unique_path;

  boost::filesystem::path base = bases;
  boost::filesystem::path temppath = temppaths;
  
  base = base.filename();
  if( temppath.empty() || !boost::filesystem::is_directory(temppath) )
    temppath = UtilityFunctions::temp_dir();

  temppath /= (base.generic_string()  + "_%%%%-%%%%-%%%%-%%%%");

  return unique_path( temppath ).string<std::string>();
}//path temp_file_name( path base )

  /*
   Could replace recursive_ls_internal() with the following to help get rid of linking to boost
#ifdef _WIN32
   //https://stackoverflow.com/questions/2314542/listing-directory-contents-using-c-and-windows
   
  bool ListDirectoryContents(const char *sDir)
  {
    WIN32_FIND_DATA fdFile;
    HANDLE hFind = NULL;
    
    char sPath[2048];
    
    //Specify a file mask. *.* = We want everything!
    sprintf(sPath, "%s\\*.*", sDir);
    
    if((hFind = FindFirstFile(sPath, &fdFile)) == INVALID_HANDLE_VALUE)
    {
      printf("Path not found: [%s]\n", sDir);
      return false;
    }
    
    do
    {
      //Find first file will always return "."
      //    and ".." as the first two directories.
      if(strcmp(fdFile.cFileName, ".") != 0
         && strcmp(fdFile.cFileName, "..") != 0)
      {
        //Build up our file path using the passed in
        //  [sDir] and the file/foldername we just found:
        sprintf(sPath, "%s\\%s", sDir, fdFile.cFileName);
        
        //Is the entity a File or Folder?
        if(fdFile.dwFileAttributes &FILE_ATTRIBUTE_DIRECTORY)
        {
          printf("Directory: %s\n", sPath);
          ListDirectoryContents(sPath); //Recursion, I love it!
        }
        else{
          printf("File: %s\n", sPath);
        }
      }
    }
    while(FindNextFile(hFind, &fdFile)); //Find the next file.
    
    FindClose(hFind); //Always, Always, clean things up!
    
    return true;
  }
#else
  //From https://rosettacode.org/wiki/Walk_a_directory/Recursively#C
  enum {
    WALK_OK = 0,
    WALK_BADPATTERN,
    WALK_NAMETOOLONG,
    WALK_BADIO,
  };
  
#define WS_NONE    0
#define WS_RECURSIVE  (1 << 0)
#define WS_DEFAULT  WS_RECURSIVE
#define WS_FOLLOWLINK  (1 << 1)  // follow symlinks
#define WS_DOTFILES  (1 << 2)  // per unix convention, .file is hidden
#define WS_MATCHDIRS  (1 << 3)  // if pattern is used on dir names too
  
  int walk_recur(char *dname, regex_t *reg, int spec)
  {
    struct dirent *dent;
    DIR *dir;
    struct stat st;
    char fn[FILENAME_MAX];
    int res = WALK_OK;
    int len = strlen(dname);
    if (len >= FILENAME_MAX - 1)
      return WALK_NAMETOOLONG;
    
    strcpy(fn, dname);
    fn[len++] = '/';
    
    if (!(dir = opendir(dname))) {
      warn("can't open %s", dname);
      return WALK_BADIO;
    }
    
    errno = 0;
    while ((dent = readdir(dir))) {
      if (!(spec & WS_DOTFILES) && dent->d_name[0] == '.')
        continue;
      if (!strcmp(dent->d_name, ".") || !strcmp(dent->d_name, ".."))
        continue;
      
      strncpy(fn + len, dent->d_name, FILENAME_MAX - len);
      if (lstat(fn, &st) == -1) {
        warn("Can't stat %s", fn);
        res = WALK_BADIO;
        continue;
      }
      
      // don't follow symlink unless told so
      if (S_ISLNK(st.st_mode) && !(spec & WS_FOLLOWLINK))
        continue;
      
      // will be false for symlinked dirs
      if (S_ISDIR(st.st_mode)) {
        // recursively follow dirs
        if ((spec & WS_RECURSIVE))
          walk_recur(fn, reg, spec);
        
        if (!(spec & WS_MATCHDIRS)) continue;
      }
      
      // pattern match
      if (!regexec(reg, fn, 0, 0, 0)) puts(fn);
    }
    
    if (dir) closedir(dir);
    return res ? res : errno ? WALK_BADIO : WALK_OK;
  }
  
  int walk_dir(char *dname, char *pattern, int spec)
  {
    regex_t r;
    int res;
    if (regcomp(&r, pattern, REG_EXTENDED | REG_NOSUB))
      return WALK_BADPATTERN;
    res = walk_recur(dname, &r, spec);
    regfree(&r);
    
    return res;
  }
#endif
*/
  
vector<std::string> recursive_ls_internal( const std::string &sourcedir,
                                           file_match_function_t match_fcn,
                                           void *user_match_data,
                                           const size_t depth,
                                           const size_t numfiles )
{
  //See http://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
  //  (and scroll down) for how to implement not using boost.
  //OR see a straight C example https://rosettacode.org/wiki/Walk_a_directory/Recursively#C
  using namespace boost::filesystem;
  
  const size_t max_depth = 25;
  const size_t max_num_files = 100000;
  
  vector<string> files;

  /*
   //A shorter untested implementation, that might be better.
  for( recursive_directory_iterator iter(sourcedir), end; iter != end; ++iter )
  {
    const std::string name = iter->path().filename().string();
    const bool isdir = UtilityFunctions::is_directory( name );
    
    if( !isdir && (!match_fcn || match_fcn(name,user_match_data)) )
      files.push_back( name );
    
    if( files.size() >= max_num_files )
      break;
  }
  return files;
  */
  
  
  if( depth >= max_depth )
    return files;
  
  if ( !UtilityFunctions::is_directory( sourcedir ) )
    return files;
  
  directory_iterator end_itr; // default construction yields past-the-end
  
  directory_iterator itr;
  try
  {
    itr = directory_iterator( sourcedir );
  }catch( std::exception & )
  {
    //ex: boost::filesystem::filesystem_error: boost::filesystem::directory_iterator::construct: Permission denied: "..."
    return files;
  }
  
  for( ; (itr != end_itr) && ((files.size()+numfiles) < max_num_files); ++itr )
  {
    const boost::filesystem::path &p = itr->path();
    const string pstr = p.string<string>();
    
    const bool isdir = UtilityFunctions::is_directory( pstr );
    
    if( isdir )
    {
      const vector<string> r = recursive_ls_internal( pstr, match_fcn, user_match_data, depth + 1, files.size() );
      files.insert( files.end(), r.begin(), r.end() );
    }else if( UtilityFunctions::is_file( pstr ) ) //if ( itr->leaf() == patern ) // see below
    {
      if( !match_fcn || match_fcn(pstr,user_match_data) )
        files.push_back( pstr );
    }//
  }//for( loop over

  return files;
}//recursive_ls(...)
  
bool filter_ending( const std::string &path, void *user_match_data )
{
  const std::string *ending = (const std::string *)user_match_data;
  return iends_with(path, *ending);
}
  
vector<std::string> recursive_ls( const std::string &sourcedir,
                                  const std::string &ending  )
{
  if( ending.empty() )
    return recursive_ls_internal( sourcedir, (file_match_function_t)0, 0, 0, 0 );
  return recursive_ls_internal( sourcedir, &filter_ending, (void *)&ending, 0, 0 );
}
  
  
std::vector<std::string> recursive_ls( const std::string &sourcedir,
                                        file_match_function_t match_fcn,
                                        void *match_data )
{
  return recursive_ls_internal( sourcedir, match_fcn, match_data, 0, 0 );
}
  

vector<string> ls_files_in_directory( const std::string &sourcedir, const std::string &ending )
{
  if( ending.empty() )
    return ls_files_in_directory( sourcedir, (file_match_function_t)0, 0 );
  return ls_files_in_directory( sourcedir, &filter_ending, (void *)&ending );
}
  
std::vector<std::string> ls_files_in_directory( const std::string &sourcedir,
                                                 file_match_function_t match_fcn,
                                                 void *user_data )
{
  using namespace boost::filesystem;
  
  vector<string> files;
  if ( !UtilityFunctions::is_directory( sourcedir ) )
    return files;
  
  directory_iterator end_itr; // default construction yields past-the-end
  
  directory_iterator itr;
  try
  {
    itr = directory_iterator( sourcedir );
  }catch( std::exception & )
  {
    //ex: boost::filesystem::filesystem_error: boost::filesystem::directory_iterator::construct: Permission denied: "..."
    return files;
  }
  
  for( ; itr != end_itr; ++itr )
  {
    const boost::filesystem::path &p = itr->path();
    const string pstr = p.string<string>();
    const bool isfile = UtilityFunctions::is_file( pstr );
    
    if( isfile )
      if( !match_fcn || match_fcn(pstr,user_data) )
        files.push_back( pstr );
  }//for( loop over
  
  return files;
}//ls_files_in_directory(...)
  
  
#if( BOOST_VERSION < 106501 )
namespace
{
  namespace fs = boost::filesystem;
  
  //Get a relative path from 'from_path' to 'to_path'
  //  assert( make_relative( "/a/b/c/d", "/a/b/foo/bar" ) == "../../foo/bar" );
  // Return path when appended to from_path will resolve to same as to_path
  fs::path make_relative( fs::path from_path, fs::path to_path )
  {
    fs::path answer;
    
    //Make the paths absolute. Ex. turn 'some/path.txt' to '/current/working/directory/some/path.txt'
    from_path = fs::absolute( from_path );
    to_path = fs::absolute( to_path );
    
    fs::path::const_iterator from_iter = from_path.begin();
    fs::path::const_iterator to_iter = to_path.begin();
    
    //Loop through each each path until we have a component that doesnt match
    for( fs::path::const_iterator to_end = to_path.end(), from_end = from_path.end();
         from_iter != from_end && to_iter != to_end && *from_iter == *to_iter;
         ++from_iter, ++to_iter )
    {
    }
    
    //Add '..' to get from 'from_path' to our the base path we found
    for( ; from_iter != from_path.end(); ++from_iter )
    {
      if( (*from_iter) != "." )
        answer /= "..";
    }
    
    // Now navigate down the directory branch
    while( to_iter != to_path.end() )
    {
      answer /= *to_iter;
      ++to_iter;
    }
    
    return answer;
  }
}//namespace
#endif
  
std::string fs_relative( const std::string &from_path, const std::string &to_path )
{
#if( BOOST_VERSION < 106501 )
  return make_relative( from_path, to_path ).string<std::string>();
#else
  return boost::filesystem::relative( to_path, from_path ).string<std::string>();
#endif
}//std::string fs_relative( const std::string &target, const std::string &base )
  
#endif //if( !SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
  
//  Windows
#ifdef _WIN32
double get_wall_time()
{
  LARGE_INTEGER time,freq;
  if( !QueryPerformanceFrequency(&freq) )
    return -std::numeric_limits<double>::max();

  if( !QueryPerformanceCounter(&time) )
    return -std::numeric_limits<double>::max();
  
  return static_cast<double>(time.QuadPart) / freq.QuadPart;
}//double get_wall_time()
#else //  Posix/Linux
double get_wall_time()
{
  //\todo Test std::chrono implementation and then get rid of Windows specialization
  //return std::chrono::time_point_cast<std::chrono::microseconds>( std::chrono::system_clock::now() ).time_since_epoch().count() / 1.0E6;
  struct timeval time;
  if( gettimeofday(&time,NULL) )
    return -std::numeric_limits<double>::max();
  return static_cast<double>(time.tv_sec) + (0.000001 * time.tv_usec);
}
#endif

double get_cpu_time()
{
  return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
}
  
double landau_cdf(double x, double xi, double x0) {
  // implementation of landau distribution (from DISLAN)
  //The algorithm was taken from the Cernlib function dislan(G110)
  //Reference: K.S.Kolbig and B.Schorr, "A program package for the Landau
  //distribution", Computer Phys.Comm., 31(1984), 97-111
  //
  //Lifted from the root/math/mathcore/src/ProbFuncMathCore.cxx file
  //  by wcjohns 20120216

   static double p1[5] = {0.2514091491e+0,-0.6250580444e-1, 0.1458381230e-1, -0.2108817737e-2, 0.7411247290e-3};
   static double q1[5] = {1.0            ,-0.5571175625e-2, 0.6225310236e-1, -0.3137378427e-2, 0.1931496439e-2};

   static double p2[4] = {0.2868328584e+0, 0.3564363231e+0, 0.1523518695e+0, 0.2251304883e-1};
   static double q2[4] = {1.0            , 0.6191136137e+0, 0.1720721448e+0, 0.2278594771e-1};

   static double p3[4] = {0.2868329066e+0, 0.3003828436e+0, 0.9950951941e-1, 0.8733827185e-2};
   static double q3[4] = {1.0            , 0.4237190502e+0, 0.1095631512e+0, 0.8693851567e-2};

   static double p4[4] = {0.1000351630e+1, 0.4503592498e+1, 0.1085883880e+2, 0.7536052269e+1};
   static double q4[4] = {1.0            , 0.5539969678e+1, 0.1933581111e+2, 0.2721321508e+2};

   static double p5[4] = {0.1000006517e+1, 0.4909414111e+2, 0.8505544753e+2, 0.1532153455e+3};
   static double q5[4] = {1.0            , 0.5009928881e+2, 0.1399819104e+3, 0.4200002909e+3};

   static double p6[4] = {0.1000000983e+1, 0.1329868456e+3, 0.9162149244e+3, -0.9605054274e+3};
   static double q6[4] = {1.0            , 0.1339887843e+3, 0.1055990413e+4, 0.5532224619e+3};

   static double a1[4] = {0, -0.4583333333e+0, 0.6675347222e+0,-0.1641741416e+1};

   static double a2[4] = {0,  1.0            ,-0.4227843351e+0,-0.2043403138e+1};

   double v = (x - x0)/xi;
   double u;
   double lan;

   if (v < -5.5) {
      u = std::exp(v+1);
      lan = 0.3989422803*std::exp(-1./u)*std::sqrt(u)*(1+(a1[1]+(a1[2]+a1[3]*u)*u)*u);
   }
   else if (v < -1 ) {
      u = std::exp(-v-1);
      lan = (std::exp(-u)/std::sqrt(u))*(p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*v)*v)*v)*v)/
         (q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*v)*v)*v)*v);
   }
   else if (v < 1)
      lan = (p2[0]+(p2[1]+(p2[2]+p2[3]*v)*v)*v)/(q2[0]+(q2[1]+(q2[2]+q2[3]*v)*v)*v);
   else if (v < 4)
      lan = (p3[0]+(p3[1]+(p3[2]+p3[3]*v)*v)*v)/(q3[0]+(q3[1]+(q3[2]+q3[3]*v)*v)*v);
   else if (v < 12) {
      u = 1./v;
      lan = (p4[0]+(p4[1]+(p4[2]+p4[3]*u)*u)*u)/(q4[0]+(q4[1]+(q4[2]+q4[3]*u)*u)*u);
   }
   else if (v < 50) {
      u = 1./v;
      lan = (p5[0]+(p5[1]+(p5[2]+p5[3]*u)*u)*u)/(q5[0]+(q5[1]+(q5[2]+q5[3]*u)*u)*u);
   }
   else if (v < 300) {
      u = 1./v;
      lan = (p6[0]+(p6[1]+(p6[2]+p6[3]*u)*u)*u)/(q6[0]+(q6[1]+(q6[2]+q6[3]*u)*u)*u);
   }
   else {
      u = 1./(v-v*std::log(v)/(v+1));
      lan = 1-(a2[1]+(a2[2]+a2[3]*u)*u)*u;
   }
   return lan;

}//double landau_cdf(double x, double xi, double x0)

  
std::istream& safe_get_line(std::istream& is, std::string& t)
{
  return safe_get_line( is, t, 0 );
}
  
  
std::istream &safe_get_line( std::istream &is, std::string &t, const size_t maxlength )
{
  //from  http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
  //  adapted by wcjohns
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.
  std::istream::sentry se( is, true );
  std::streambuf *sb = is.rdbuf();

  for( ; !maxlength || (t.length() < maxlength); )
  {
    int c = sb->sbumpc(); //advances pointer to current location by one
    switch( c )
    {
      case '\r':
        c = sb->sgetc();  //does not advance pointer to current location
        if(c == '\n')
          sb->sbumpc();   //advances pointer to one current location by one
         return is;
      case '\n':
        return is;
      case EOF:
        is.setstate( ios::eofbit );
        return is;
      default:
        t += (char)c;
    }//switch( c )
  }//for(;;)

  return is;
}//safe_get_line(...)

  
bool split_to_ints( const char *input, const size_t length,
                     std::vector<int> &results )
{
  if( !input || !length )
    return true;
  
  namespace qi = boost::spirit::qi;
  
  const char *begin = input;
  const char *end = begin + length;
  
  const bool ok = qi::phrase_parse( begin, end, (*qi::int_) % qi::eol, qi::lit(",")|qi::space, results );

  //If we didnt consume the entire input, make sure only delimiters are left
  if( ok && begin != end )
  {
    for( const char *pos = begin; pos != end; ++pos )
    {
      if( !is_in( *pos, " \t\n\r," ) )
        return false;
    }
  }//if( ok && begin != end )
  
  return ok;
}//bool split_to_ints(...)
  
bool split_to_floats( const std::string &input, std::vector<float> &results )
{
  return split_to_floats( input.c_str(), input.length(), results );
}
  
bool split_to_floats( const char *input, const size_t length,
                      vector<float> &results )
{
  results.clear();
  
  if( !input || !length )
    return true;
  
  results.reserve( length/2 );
  
  namespace qi = boost::spirit::qi;
  
  const char *begin = input;
  const char *end = begin + length;
  
  //Note that adding the comma increases the parse time by about 15% over just
  //  qi::space. Could perform a check to see if a comma comes after the first
  //  number and then choose the parser expression to use...
  //Not sure difference between qi::space and qi::blank from testing various
  //  strings
  const bool ok = qi::phrase_parse( begin, end, (*qi::float_) % qi::eol, qi::lit(",")|qi::space, results );
  
#if(PERFORM_DEVELOPER_CHECKS)
  if( !ok )
  {
    if( *input && isdigit(*input) )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg), "Parsing failed: '%s'",
               string(begin,end).c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( *input && isdigit(*input) )
  }else
  {
    if( begin != end && results.size() )
    {
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg), "Trailing unpased string '%s'",
               string(begin,end).c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( begin != end )
  
  
    vector<float> checked_result;
    string datacopy( input, input+length );
    trim( datacopy );
    split_to_floats( (char *)datacopy.c_str(), checked_result, " \t\n\r,", false );
  
    if( checked_result.size() == results.size() )
    {
      for( size_t i = 0; i < results.size(); ++i )
      {
        const float a = results[i];
        const float b = checked_result[i];
        if( fabs(a-b) > 0.000001*std::max(fabs(a),fabs(b)) )
        {
          char errormsg[1024];
          snprintf( errormsg, sizeof(errormsg),
                   "Out of tolerance diffence for floats %.9g using "
                   "boost::spirit vs %.9g using alternative split_to_float on float %i",
                    a, b, int(i) );
          log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
        }
      }//for( size_t i = 0; i < results.size(); ++i )
    }else
    {
     char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
                "Parsed wrong number of floats %i using boost::spirit and %i "
                "using strtok for '%s'",
                 int(results.size()), int(checked_result.size()),
               string(input,end).c_str() );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    }//if( checked_result.size() == results.size() ) / else
  } //if( ok ) / else
#endif
  
  //If we didnt consume the entire input, make sure only delimiters are left
  if( ok && begin != end )
  {
    for( const char *pos = begin; pos != end; ++pos )
    {
      if( !is_in( *pos, " \t\n\r," ) )
        return false;
    }
  }//if( ok && begin != end )
  
  return ok;
}//bool split_to_floats(...)
  
  
bool parse_float( const char *input, const size_t length, float &result )
{
  namespace qi = boost::spirit::qi;
  
  const char *begin = input;
  const char *end = begin + length;
  
  result = 0.0f;
  
  bool ok = qi::phrase_parse( begin, end, qi::float_, qi::space, result );
  
//  if( ok && (begin != end) )
//    return false;
    
  return ok;
}
  
bool split_to_floats( const char *input, vector<float> &contents,
                      const char * const delims,
                      const bool cambio_zero_compress_fix )
{
#if(PERFORM_DEVELOPER_CHECKS)
  if( string(delims).find_first_of( ".0123456789+-.eE" ) != string::npos )
  {
    char errormsg[1024];
    snprintf( errormsg, sizeof(errormsg), "Invalid delimiter: '%s'", delims );
    log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
    return false;
  }
#endif
  
// PARSE_CONVERT_METHOD==1 is strtod
// PARSE_CONVERT_METHOD==2 is using boost::spirit::qi::parse
// PARSE_CONVERT_METHOD==3 is using strtok_r and atof
//  If option 3 is selected, then the signture of this function must be changed
//    from a const char *.
//  Also, methods 1 and 3 both implicitly add all the whitespaces to the
//    delimiters.
//  I think option 2 is the only 'correct' implementation, although the others
//  are close enough for parsing spectrum files. Options 1 and 3 became
//  depreciated 20151226
  
#define PARSE_CONVERT_METHOD 2
  
  const size_t input_size = contents.size();
  if( input_size )
  {
    contents.clear();
    contents.reserve( input_size / 2 );
  }//if( input_size )
    
  if( !input || !(*input) )
    return false;
  
#if( PARSE_CONVERT_METHOD == 1 )
  errno = 0;
  const char *pos = input;
  char *nextpos = input;
  pos = next_word( nextpos, delims );
  
  do
  {
    const char d = *pos;
    if( !isdigit(d) && d != '+' && d != '-' && d != '.' )
      break;
    
    float value = static_cast<float>( strtod( pos, &nextpos ) );
//    cerr << "pos=" << (void *)pos << ", nextpos=" << (void *)nextpos << " and " << *nextpos << endl;
    
    if( errno )
    {
#if(PERFORM_DEVELOPER_CHECKS)
      char errormsg[1024];
      char strpart[128] = { 0 };
      for( int c = 0; c < 127 && pos[c]; ++c )
        strpart[c] = pos[c];
      strpart[127] = '\0';
        
      snprintf( errormsg, sizeof(errormsg),
               "Couldnt convert string '%s' to a float using strtod(), error %i",
                strpart, errno );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
#endif
      return false;
    }//if( errno )
    
    if( cambio_zero_compress_fix && (value == 0.0) )
    {
      const char nextchar[2] = { *(pos+1), '\0' };
      if( !strstr(delims, nextchar) )  //If the next char is a delimeter, then we dont wantto apply the fix, otherwise apply the fix
        value = FLT_MIN;
    }//if( value == 0.0 )
    
    contents.push_back( value );
    
    pos = next_word( nextpos, delims );
  }while( pos && (*pos) && (pos != nextpos) );

#elif( PARSE_CONVERT_METHOD == 2 )
  
  const char *pos = input;
  
  while( *pos )
  {
    pos = next_word( pos, delims );
    if( *pos == '\0' )
      return true;
    
    const char * const start_pos = pos;
    const char *end = end_of_word( pos, delims );
    
    //Using a double here instead of a float causes about a 2.5% slow down.
    //  Using a float would be fine, but then you hit a limitation
    //  that the value before decimal point must be less than 2^32 (note, we are
    //  unlikely to ever see this in channel counts).
    double value;
    const bool ok = boost::spirit::qi::parse( pos, end, boost::spirit::qi::double_, value );
    
    
#if(PERFORM_DEVELOPER_CHECKS)
    if( !ok )
    {
      if( input && isdigit(*input) )
      {
        char errormsg[1024];
        snprintf( errormsg, sizeof(errormsg), "Parsing failed: '%s'",
                  string(start_pos, end).c_str() );
        log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
      }
    }else
    {
      if( pos != end )
      {
        char errormsg[1024];
        snprintf( errormsg, sizeof(errormsg), "Trailing unpased string '%s'",
                 string(pos,end).c_str() );
        log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
      }//if( begin != end )
      
      const float a = static_cast<float>( value );
      const float b = (float) atof( string(start_pos, end).c_str() );
      if( fabs(a-b) > 0.000001f*std::max(fabs(a),fabs(b)) )
      {
        char errormsg[1024];
        snprintf( errormsg, sizeof(errormsg),
                 "Out of tolerance diffence for floats %.9g using "
                 "boost::spirit vs %.9g using alternative split_to_float on float %i",
                 a, b, int(contents.size()) );
        log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
      }
    }//if( !ok )
#endif //PERFORM_DEVELOPER_CHECKS
    
    if( !ok )
      return false;
    
    //If the next char is a delimeter, then we dont wantto apply the fix, otherwise apply the fix
    if( cambio_zero_compress_fix && (value == 0.0)
        && !is_in( *(start_pos+1), delims ) )
    {
      value = FLT_MIN;
    }//if( value == 0.0 )
    
    contents.push_back( static_cast<float>(value) );
  }//while( *pos )
  
#elif( PARSE_CONVERT_METHOD == 3 )
  errno = 0;
  float value;
  char *pos_ptr = NULL;
  
  // Branches for Windows; strtok_r is on POSIX systems. strtok_s is the Windows equivalent.
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  char *value_str = strtok_s( input, delims, &pos_ptr );
#else
  char *value_str = strtok_r( input, delims, &pos_ptr );
#endif // #if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
    
  while( value_str != NULL )
  {
    //  XXX: using stringstream to convert the double makes split_to_floats(...)
    //       take about 40% of the time of parsing an ICD1 file (the example Anthony
    //       NM Passthrough) - using atof reduces this to 17.9% (and further down to
    //       14% if 'contents' is passed in with space already reserved).
    //       Running on a number of ICD1 files showed no errors using the
    //       stringstream implementation, so will continue using atof.
    //       Now approx 40% of time in this function is due to
    //       vector<float>::push_back(...), and 50% to atof(...); 11% to strtok_r
    //
    //    if( !(stringstream(value_str) >> value) )
    //      cerr << "Error converting '" << value_str << "' to float" << endl;
    //
    //Note that atof discards initial white spaces, meaning the functionality
    //  of this function is not correct if the delimiters dont include all the
    //  white space characters.
//    errno = (sscanf(value_str, "%f", &value) != 1);
//    value = fast_atof( value_str );
//    errno = !naive_atof( value, value_str );
    value = (float)atof( value_str );
//    value = (float)strtod( value_str, NULL );
/*
    {
      //find next delim or end
      const char *end = end_of_word( value_str, delims );
      errno = !boost::spirit::qi::parse( (const char *)value_str, end, boost::spirit::qi::float_, value );
      value_str = (char *)end;
    }
*/
    if( errno )
    {
#if(PERFORM_DEVELOPER_CHECKS)
      char errormsg[1024];
      snprintf( errormsg, sizeof(errormsg),
                "Couldnt convert string '%s' to a float using atof(), error %i",
                value_str, errno );
      log_developer_error( BOOST_CURRENT_FUNCTION, errormsg );
#endif
      return false;
    }//if( errno )

      
    //XXX
    //  In Cambio N42 files zeros written like "0" indicates zeroes compression,
    //  a zero written like "0.000", "-0.0", etc are all a sinle bin with
    //  content zero - so for this case well substiture in a really small number
    if( cambio_zero_compress_fix && (value == 0.0) )
    {
      if( strlen(value_str) > 1 )
        value = FLT_MIN; //std::numeric_limits<float>::min();
    }//if( value == 0.0 )
      
    contents.push_back( value );
      
    // Branches for Windows; strtok_r is on POSIX systems. strtok_s is the Windows equivalent.
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
    value_str = strtok_s( NULL, delims, &pos_ptr );
#else
    value_str = strtok_r( NULL, delims, &pos_ptr );
#endif // #if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  }//while( pch != NULL )
  
#else
#error Invalid parse method specified
#endif
  
  return true;
}//vector<float> split_to_floats( ... )

std::string sequencesToBriefString( const std::set<int> &sample_numbers )
{
  stringstream editVal;
  
  int added = 0;
  int firstInRange = *(sample_numbers.begin());
  int previous = firstInRange;
  
  for( set<int>::const_iterator iter = sample_numbers.begin();
      iter != sample_numbers.end(); ++iter )
  {
    const int thisval = *iter;
    
    if( (thisval > (previous+1)) )
    {
      editVal << string(added ? "," : "");
      if( previous == firstInRange )
        editVal << previous;
      else if( previous == (firstInRange+1) )
        editVal << firstInRange << "," << previous;
      else
        editVal << firstInRange << "-" << previous;
      
      ++added;
      firstInRange = thisval;
    }//if( thisval > (previous+1) )
      
    previous = thisval;
  }//for( loop over smaple_numbers )
      
  editVal << string(added ? "," : "");
  if( previous == firstInRange )
    editVal << previous;
  else if( previous == (firstInRange+1) )
    editVal << firstInRange << "," << previous;
  else
    editVal << firstInRange << "-" << previous;
  
  return editVal.str();
}
  
  
unsigned int levenshtein_distance( const string &source, const string &target )
{
  //This function largely derived from code found at:
  //  http://www.merriampark.com/ldcpp.htm  (by Anders Johnasen).
  //  There was no accompaning lincense information, and since I found many
  //  similar-but-seperate implementations on the internet, I take the code in
  //  this function to be licensed under a 'do-what-you-will' public domain
  //  license. --Will Johnson 20100824
  //This function is case insensitive.
  const size_t n = source.length();
  const size_t m = target.length();
  if( !n )
    return static_cast<unsigned int>(m);
  
  if( !m )
    return static_cast<unsigned int>(n);
  
  vector< vector<size_t> > matrix( n+1, vector<size_t>(m+1,0) );
  
  for( size_t i = 0; i <= n; i++)
    matrix[i][0]=i;
  
  for( size_t j = 0; j <= m; j++)
    matrix[0][j]=j;
  
  for( size_t i = 1; i <= n; i++)
  {
    const string::value_type s_i = std::toupper( source[i-1] );
    
    for( size_t j = 1; j <= m; j++)
    {
      const string::value_type t_j = std::toupper( target[j-1] );
      
      size_t cost = (s_i == t_j) ? 0 : 1;
      
      const size_t above = matrix[i-1][j];
      const size_t left = matrix[i][j-1];
      const size_t diag = matrix[i-1][j-1];
      size_t cell = min( above + 1, min(left + 1, diag + cost));
      
      // Step 6A: Cover transposition, in addition to deletion,
      // insertion and substitution. This step is taken from:
      // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's
      // Enhanced Dynamic Programming ASM Algorithm"
      // (http://www.acm.org/~hlb/publications/asm/asm.html)
      
      if( i>2 && j>2 )
      {
        size_t trans = matrix[i-2][j-2]+1;
        if( source[i-2] != t_j )
          trans++;
        if( s_i != target[j-2] )
          trans++;
        if( cell > trans )
          cell = trans;
      }//if()
      
      matrix[i][j]=cell;
    }//for( loop over 'target' letters j )
  }//for( loop over 'source' letters i )
    
  return static_cast<unsigned int>(matrix[n][m]);
}//unsigned int levenshtein_distance( const string &, const string &)

}//namespace UtilityFunctions
