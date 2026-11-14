#ifndef UserLibraryWindow_h
#define UserLibraryWindow_h
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
#include <cstdint>

#include <Wt/WString.h>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/UserLibrary.h"

class InterSpec;
namespace Wt
{
  class WMenu;
  class WText;
  class WAnchor;
  class WLineEdit;
  class WTextArea;
  class WComboBox;
  class WPushButton;
  class WStackedWidget;
  class WContainerWidget;
}


/** The Packages window (File -> Packages…): installed packages, every item by type (with
 which one is used for a name, and enabling/disabling them), and creating packages.

 A singleton tool window (CLAUDE.md "Cat-A"): InterSpec::showUserLibrary / InterSpec::closeUserLibrary.
 Listing reads files (and validates items), so it runs on a worker thread; results come back by widget
 id.  All text from packages is shown as plain text.
 */
class UserLibraryWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  virtual ~UserLibraryWindow();

  /** Re-lists everything; called when the library changes (InterSpec::userLibraryChanged). */
  void refresh();

  /** The translated name of a type (e.g. "Detector efficiency functions"). */
  static Wt::WString typeTitle( const std::string &typeId );

  /** "InterSpec", "My files", or "Package: <name>". */
  static Wt::WString sourceLabel( const UserLibrary::Source &source );

  /** The label of the place a type's items are used (e.g. "Detector Efficiency"); empty if there
   is no place to open.
   */
  static Wt::WString whereUsedLabel( const std::string &typeId );

  /** Opens the tool where a type's items are used.  Session thread. */
  static void showWhereUsed( const std::string &typeId );

  struct PackageInfo;
  struct ItemRow;
  struct Candidate;

protected:
  explicit UserLibraryWindow( InterSpec *viewer );

  // Packages page
  void refreshPackages();
  void packagesListed( std::shared_ptr<const std::vector<PackageInfo>> packages );
  void confirmUninstall( const std::string &packageId, const std::string &name );
  void startUpload();

  // By-type page
  void refreshItems();
  void itemsListed( const uint64_t request, std::shared_ptr<const std::vector<ItemRow>> rows );
  void renderItems();
  void confirmDelete( const size_t index );
  void libraryChangedHere();

  // Create page
  void refreshCandidates();
  void candidatesListed( const uint64_t request, std::shared_ptr<const std::vector<Candidate>> candidates );
  void updateStartFromOptions();
  void startFromChanged();
  void renderCandidates();
  void startCreate();
  void createDone( const std::string &error, std::shared_ptr<const std::string> data,
                   const UserLibrary::Manifest &manifest, const std::vector<std::string> &issues,
                   const std::vector<std::string> &warnings );

  InterSpec *m_interspec;
  Wt::WMenu *m_menu;
  Wt::WStackedWidget *m_stack;

  // Packages page
  Wt::WContainerWidget *m_packagesList;
  uint64_t m_packagesRequest;

  // By-type page
  Wt::WComboBox *m_typeSelect;
  Wt::WLineEdit *m_itemSearch;
  Wt::WPushButton *m_showMe;
  Wt::WContainerWidget *m_itemsList;
  std::shared_ptr<const std::vector<ItemRow>> m_rows;
  uint64_t m_itemsRequest;

  // Create page
  Wt::WComboBox *m_startFrom;
  Wt::WComboBox *m_minVersion;              ///< null before InterSpec 1.2 (packages need not say)
  std::vector<std::string> m_startFromIds;   ///< "" for a new package, else a created package's id
  std::map<std::string,Wt::WLineEdit *> m_fields;
  std::map<std::string,Wt::WTextArea *> m_areas;
  Wt::WLineEdit *m_candidateSearch;
  Wt::WText *m_numSelected;
  Wt::WContainerWidget *m_candidatesList;
  Wt::WContainerWidget *m_createStatus;
  Wt::WPushButton *m_createBtn;
  std::shared_ptr<const std::vector<Candidate>> m_candidates;
  std::map<std::string,bool> m_selected;          ///< by source key
  std::map<std::string,std::string> m_itemNotes;  ///< by source key
  std::map<std::string,bool> m_itemOverride;      ///< by source key
  std::map<std::string,std::string> m_changedSinceLast;  ///< source key -> previous version
  std::vector<std::string> m_movedSinceLast;      ///< package paths no longer found
  uint64_t m_candidatesRequest;
  bool m_creating;
};//class UserLibraryWindow

#endif //UserLibraryWindow_h
