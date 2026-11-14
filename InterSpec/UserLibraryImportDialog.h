#ifndef UserLibraryImportDialog_h
#define UserLibraryImportDialog_h
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

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/UserLibraryPackage.h"

class InterSpec;
namespace Wt
{
  class WText;
  class WCheckBox;
  class WPushButton;
  class WRadioButton;
  class WContainerWidget;
}


/** The dialog shown when a user library package is dropped onto InterSpec (or chosen in the Packages
 window): it checks the package on a worker thread, shows what it holds and how it relates to
 what is installed, asks about name conflicts and sensitive items, and installs it.

 A multi-instance dialog (CLAUDE.md "Cat-B"): create it with #start, which wires `finished()` to
 AuxWindow::deleteAuxWindow.  All text from the package is shown as plain text.
 */
class UserLibraryImportDialog : public AuxWindow
{
  friend class AuxWindow;

public:
  /** Copies `archivePath` (so the caller may delete it right away), shows the dialog, and starts
   checking the package.  Must be called on the session thread.
   */
  static UserLibraryImportDialog *start( const std::string &archivePath, const std::string &displayName );

  virtual ~UserLibraryImportDialog();

protected:
  UserLibraryImportDialog( const std::string &displayName, const std::string &archiveCopy );

  /** Called on the session thread when #UserLibraryPackage::inspect finishes. */
  void inspectionDone( std::shared_ptr<UserLibraryPackage::Inspection> inspection );

  /** Called on the session thread when the install finishes; empty `error` on success. */
  void installDone( const std::string &error, const std::string &packageName );

  void render();
  void setAllConflicts( const bool useIncoming );
  void updateInstallButton();
  void startInstall();
  void showError( const Wt::WString &message );

  std::string m_displayName;
  std::string m_archiveCopy;  ///< deleted in the destructor
  std::shared_ptr<UserLibraryPackage::Inspection> m_inspection;
  bool m_installing;
  bool m_installed;

  Wt::WContainerWidget *m_body;
  Wt::WPushButton *m_install;
  Wt::WPushButton *m_cancel;

  /** (type id, name) -> the "use incoming" radio button of an open conflict. */
  std::map<std::pair<std::string,std::string>,Wt::WRadioButton *> m_useIncoming;
  std::map<std::pair<std::string,std::string>,Wt::WRadioButton *> m_keepExisting;
  std::vector<Wt::WCheckBox *> m_sensitiveConfirms;
};//class UserLibraryImportDialog

#endif //UserLibraryImportDialog_h
