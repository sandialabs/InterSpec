#ifndef ElectronUtils_h
#define ElectronUtils_h
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

#include <string>

namespace ElectronUtils
{
  /** Must specify:
      - --userdatadir "/path/to/writable/app/data/dir"
      - --externalid "SomeUniqueString"
      - --ipc "<port number>"
   */
  int run_app( int argc, char *argv[] );
  
  /** The externalid passed into #run_app. */
  std::string external_id();
  
  /** Requests main.js to load a new clean session (i.e., don restore any state)
   This is a workaround to when the user requests a new session, the normal
   mechanism in c++ creates duplicate Electron menu items...
   
   @returns whether message was succesfully sent or not.
   */
  bool requestNewCleanSession();
  
  
  /** Sends a message through the websocket connection letting main.js know that
   the session has loaded.
   
   @param The 'externalid' URL argument so main.js knows which session has loaded.
   @returns whether message was succesfully sent or not.
   
   Note: main.js will wait till recieveing this notification before asking the
   session to open any files the OS requested.
   */
  bool notifyNodeJsOfNewSessionLoad( const std::string externalid );
}//namespace ElectronUtils


#endif  //#ifndef ElectronUtils_h
