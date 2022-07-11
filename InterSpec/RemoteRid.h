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
   
   Checks if user has previously asked to never warn about remote RID capabilities, and if not
   will present a dialog to the user warning them about information leaving the application;
   if they cancel, will close dialog, and the callback will be called with nullptrs.
   If they accept, a AuxWindow with a RemoteRid widget will be created, and the callback will be
   called with relevant pointers.
   */
  static void startRemoteRidDialog( InterSpec *viewer, std::function<void(AuxWindow *, RemoteRid *)> callback );
  
  /** Immediately creates a RemoteRid dialog, without asking the user, etc. */
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
  
  static void startAutomatedOnLoadAnalysis( InterSpec *interspec,
                                            const Wt::WFlags<AnaFileOptions> flags = 0 );
  
  static void handleOpeningRemoteRidTool( InterSpec *interspec );
  static void disableAutoRemoteRid( InterSpec *interspec );
  static void handleShowingRemoteRidRefLines( InterSpec *interspec, std::string signal);
  
protected:
  InterSpec *m_interspec;
  
  RestRidImp::ExternalRidWidget *m_rest_rid;
  
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  RestRidImp::ExternalRidWidget *m_exe_rid;
#endif

};//class RemoteRid

#endif //RemoteRid_h
