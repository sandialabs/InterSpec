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
 
#import <UIKit/UIKit.h>
#import <Foundation/NSData.h>

#include <string>
#include <fstream>
#include <iostream>

#include "InterSpec/SpecMeas.h"
//#include "SpecUtils/UtilityFunctions.h"
#include "target/ios/InterSpec/FileHandling.h"
#include "target/ios/InterSpec/InterSpec/AppDelegate.h"
//#include "target/ios/InterSpec/InterSpec/ViewController.h"
/*
bool ends_with( const std::string &line, const std::string &label )
{
  if( line.length() < label.length() )
    return false;
  if( label.empty() || line.empty() )
    return false;

  std::string ending = line.substr( line.length() - label.length() );
  return (ending == label);
}
*/

void remove_file_extention( std::string &line )
{
  const size_t pos = line.find( '.' );
  if( pos != std::string::npos )
    line = line.substr( 0, pos );
}

void hand_spectrum_to_other_app( std::shared_ptr<SpecMeas> spec )
{
  NSLog( @"In exportForeground!" );

  if( !spec )
    return;

  AppDelegate *appDelegate = (AppDelegate *)[[UIApplication sharedApplication] delegate];

  //NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
  //NSString *outputfile = [NSString pathWithComponents: [NSArray arrayWithObjects:docDir, @"filetotransfer.n42", nil ] ];
  //NSString *docDir = [paths objectAtIndex:0];
  std::string filename = spec->filename();
  if( filename.empty() )
    filename = "spectrum";

  remove_file_extention( filename );
  filename += ".n42";

  NSString *filenameNSS = [NSString stringWithUTF8String: filename.c_str()];
  NSString *outputfile = [NSTemporaryDirectory() stringByAppendingPathComponent: filenameNSS];
  const char *cfilestr = [outputfile cStringUsingEncoding:NSUTF8StringEncoding];

  {
    std::ofstream outfile( cfilestr );
    spec->write_2012_N42( outfile );
  }

  dispatch_async(dispatch_get_main_queue(),^{
    [appDelegate sendSpectrumFileToOtherApp: outputfile];

    //Deleteing the file here will cause it to not actually get transfered.
    //  I think the file should be deleted by UIDocumentInteractionControllerDelegate
    //  but I cant actually get the appropriate methods called...
  });

}
