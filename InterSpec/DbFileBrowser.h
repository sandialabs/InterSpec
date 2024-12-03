#ifndef DbFileBrowser_h
#define DbFileBrowser_h
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

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/SpectraFileModel.h"

namespace Wt
{
  class WPushButton;
  class WTree;
}//namespace Wt

namespace SpecUtils{ enum class SpectrumType : int; }

/**
 Widget to hold all the UI components to show the snapshot/spectra (used in both
 load dialog and welcome dialog).
 
 ToDo:
   - add signals when there is an error, or a selection is made
 
 */
class SnapshotBrowser : public Wt::WContainerWidget
{
public:
  SnapshotBrowser( SpecMeasManager *manager,
                   InterSpec *viewer,
                   std::shared_ptr<SpectraFileHeader> header,
                   Wt::WContainerWidget *buttonBar = nullptr,
                   Wt::WContainerWidget *parent = nullptr );
  
   Wt::Signal<> &finished();
  
  void loadSnapshotSelected();
  void loadSpectraSelected();

  int numSnaphots() const;
  
  /**  Returns query to find user states.
   
   @param user A pointer to the user the state must belong to
   @param session The DbSession to use
   @param header If a non-null (i.e., valid) header is passed in, then the query returns all of the
          snapshots that contain that specific spectrum. If a header is nullptr, then query returns
          all saved states for the user (just the upper level states, not their snapshots).
   @returns All states belonging to the user that contained the specified spectrum file, or if no
          spectrum file was specified, all user saved states.
   
   
   Note, I'm pretty sure you need an active transaction before calling this function - e.g.,
   \code
   DataBaseUtils::DbTransaction transaction( *m_session );
   \endcode
   */
  static Wt::Dbo::collection<Wt::Dbo::ptr<UserState> >
  get_user_states_collection( Wt::Dbo::ptr<InterSpecUser> user,
                              std::shared_ptr<DataBaseUtils::DbSession> &session,
                              std::shared_ptr<const SpectraFileHeader> header );
  
  /** Returns the number of saved states available. */
  static size_t num_saved_states( InterSpec *viewer,
                                 std::shared_ptr<DataBaseUtils::DbSession> session,
                                 std::shared_ptr<const SpectraFileHeader> header );
  
protected:
  void selectionChanged();
  void startDeleteSelected();
  void startEditSelected();
  
  void deleteSelected();
  
  std::map < Wt::WTreeNode *, Wt::Dbo::ptr < UserState > > m_UserStateLookup;
  std::map < Wt::WTreeNode *, Wt::Dbo::ptr < UserFileInDb > > m_UserFileInDbLookup;
  
  void addSpectraNodes(Wt::Dbo::collection< Wt::Dbo::ptr<UserState> >::const_iterator versionIterator, Wt::WTreeNode *versionNode);
  
  std::shared_ptr<SpecMeas> retrieveMeas( const int dbid );
  
  std::shared_ptr<DataBaseUtils::DbSession> m_session;
  SpecMeasManager  *m_manager;
  InterSpec        *m_viewer;
  Wt::WPushButton  *m_loadSnapshotButton;
  Wt::WPushButton  *m_loadSpectraButton;
  Wt::WPushButton  *m_renameButton;
  Wt::WButtonGroup *m_buttonGroup;
  Wt::WGroupBox    *m_buttonbox;
  Wt::WTree        *m_snapshotTable;
  Wt::WText        *m_descriptionLabel;
  Wt::WText        *m_timeLabel;
  Wt::WGridLayout  *m_relatedSpectraLayout;
  std::shared_ptr<SpectraFileHeader> m_header;
  
  Wt::Signal<> m_finished;
  AuxWindow *m_editWindow;
  
  int m_nrows;
};//class SnapshotBrowser


/** Widget that displays the passed-in spectra in the database, so the user can choose if
 they want to load any of them.
 
 @param modifiedFiles The files in the database to display
 @param type How to load the selected spectrum as (foreground, background, secondary)
 @param model The SpectraFileModel to use for performing operations
 @param manager The SpecMeasManager to use for performing operations
 @param header The optional SpectraFileHeader to remove when loading a database entry file.
 @param parent The Wt parent of this widget
 
 TODO: Make this a MVC widget (e.g., a WTableView) so we can potentially display all users spectra files in the database, without putting a million things into the DOM.  Not sure if we should then make it backed by a DB query, or a vector holding all informations.  This would let us have it be sortable, etc.  Would need to implement a custom delegate, at least for the load column.  After doing this we should make DbFileBrowser have the same imlementation as SpecMeasManager::showPreviousSpecFileUsesDialog(), so users can browse all their auto-saved spectra.  Then we should also maybe let users filter thier selection and such.
 */
class AutosavedSpectrumBrowser : public Wt::WContainerWidget
{
public:
  AutosavedSpectrumBrowser( const std::vector<Wt::Dbo::ptr<UserFileInDb>> &modifiedFiles,
                            SpecUtils::SpectrumType type,
                            SpectraFileModel *model,
                            SpecMeasManager *manager,
                            std::shared_ptr<SpectraFileHeader> header,
                            WContainerWidget *parent = nullptr );
  
  Wt::Signal<> &loadedASpectrum();
  
protected:
  void entryWasLoaded();
  
  Wt::Signal<> m_loadedASpectrum;
};//class AutosavedSpectrumBrowser


/**
 Main class to load snapshot/spectra.  Uses SnapshotBrowser to fill in the UI
 */
class DbFileBrowser : public AuxWindow
{
public:
  DbFileBrowser( SpecMeasManager *manager, InterSpec *viewer,
                 std::shared_ptr<SpectraFileHeader> header );
  
  int numSnapshots() const;
protected:
  SnapshotBrowser *m_factory;
};

#endif //DbFileBrowser_h
