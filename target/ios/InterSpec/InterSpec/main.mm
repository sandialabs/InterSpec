//
//  main.m
//  InterSpec
//
//  Created by Johnson, William C on 10/10/13.
//  Copyright (c) 2013 Sandia National Laboratories. All rights reserved.
//

#import <UIKit/UIKit.h>

#import "AppDelegate.h"


#include <string>
#include <iostream>

class iosbuf: public std::streambuf
{
  //A utility to redirect cout/cerr to the Android logging system so it can be
  // seen using 'adb logcat'
public:
  enum Source { FromCout, FromCerr  };
  enum { bufsize = 512 }; // ... or some other suitable buffer size
  iosbuf( Source src )
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
  }//iosbuf( Source src )
  
  ~iosbuf()
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
  }//~iosbuf()
  
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
      NSLog( @"%s", std::string(this->pbase(), this->pptr()).c_str() );
      this->setp(buffer, buffer + bufsize - 1);
    }
    return rc;
  }//int sync()
  char buffer[bufsize];
  const char *m_source;
  const Source m_type;
  std::streambuf *m_origSrc;
};//class iosbuf


int main(int argc, char *argv[])
{
  iosbuf stdbuf( iosbuf::FromCout );
  iosbuf errbuf( iosbuf::FromCerr );
      
  @autoreleasepool {
      return UIApplicationMain(argc, argv, nil, NSStringFromClass([AppDelegate class]));
  }
}
