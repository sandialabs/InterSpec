#ifndef RemoteRid_h
#define RemoteRid_h
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

#include <utility>
#include <functional>

#include <Wt/WContainerWidget>

class AuxWindow;
class InterSpec;

namespace SpecUtils { class SpecFile; }

namespace RestRidImp { class ExternalRidWidget; }

/** Interface to the Full-Spectrum nuclide ID algorithm, either to the REST web service, or to
 the command line interface of the executable.
 */
class RemoteRid : public Wt::WContainerWidget
{
public:
  /** Starts a remote RID dialog sequence.
   
   @param callback Function called once the Remote RID tool is created, optional.
   @returns Returns pointer to the "warning" dialog (about information being sent from your device),
            if it was created; will be nullptr if user has opted out of this.
   
   Checks if user has previously asked to never warn about remote RID capabilities, and if not
   will present a dialog to the user warning them about information leaving the application;
   if they cancel, will close dialog, and the callback will be called with nullptrs.
   If they accept, a AuxWindow with a RemoteRid widget will be created, and the callback will be
   called with relevant pointers.
   */
  static SimpleDialog *startRemoteRidDialog( InterSpec *viewer, std::function<void(AuxWindow *, RemoteRid *)> callback );
  
  /** Immediately creates a RemoteRid dialog, without asking the user, etc.
   
   The dialog will emit the `finished()` signal when done, and then you will need to delete it.
   */
  static std::pair<AuxWindow *, RemoteRid *> createDialog( InterSpec *viewer );
  
public:
  RemoteRid( InterSpec *viewer, Wt::WContainerWidget *parent );
  
  enum AnaFileOptions
  {
    OnlyDisplayedSearchSamples = 0x01
  };//enum AnaFileOptions
  
  static std::shared_ptr<SpecUtils::SpecFile> fileForAnalysis( InterSpec *interspec,
                                                      const Wt::WFlags<AnaFileOptions> flags = 0 );
  
  void alwaysCallRestAnaChecked();
  
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  void alwaysCallExeAnaChecked();
#endif
  
  RestRidImp::ExternalRidWidget *restRidTool();
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  RestRidImp::ExternalRidWidget *exeRidTool();
#endif
  
  static void startAutomatedOnLoadAnalysis( InterSpec *interspec,
                                            const Wt::WFlags<AnaFileOptions> flags = 0 );
  
  static void handleOpeningRemoteRidTool( InterSpec *interspec );
  static void disableAutoRemoteRid( InterSpec *interspec );
  static void handleShowingRemoteRidRefLines( InterSpec *interspec, std::string signal);
  
  /** Handles receiving a "deep-link" url starting with "interspec://remoterid?ver=1&url=...&always=1".
   
   Example URIs:
   - "interspec://remoterid?ver=1&url=https%3A%2F%2Ffull-spectrum.sandia.gov%2Fapi%2Fv1%2Finfo&always=1"
   - "interspec://remoterid?ver=1&path=C%3A%5Csome%5Cpath%5Cto%5Ca.exe&always=0"
   - "interspec://remoterid?ver=1&none=1"
   
   @param query_str The query portion of the URI.
   This string is is in standard URL format of "key1=value1&key2=value2&..." with ordering not mattering.
   Assumes the string passed in has alaready been url-decoded, but the url or path fields are URL encoded..
   If not a valid path or query_str, throws exception.
   */
  static void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.
   
   Returned string is just query portion of of URL, so will look something like:
    "ver=1&url=https%3A%2F%2Ffull-spectrum.sandia.gov%2Fapi%2Fv1%2Finfo&always=1",
   and should be url-encoded again before put in a URL or QR-code.
   */
  std::string encodeStateToUrl() const;
protected:
  InterSpec *m_interspec;
  
  RestRidImp::ExternalRidWidget *m_rest_rid;
  
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  Wt::WMenu *m_menu;
  RestRidImp::ExternalRidWidget *m_exe_rid;
#endif

};//class RemoteRid

#endif //RemoteRid_h
