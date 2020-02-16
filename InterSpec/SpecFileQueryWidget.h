#ifndef SpecFileQueryWidget_h
#define SpecFileQueryWidget_h
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
#include <string>
#include <vector>
#include <memory>
#include <atomic>
#include <iostream>

#include <boost/any.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#if( BUILD_AS_ELECTRON_APP )
#include <Wt/WEvent>
#endif
#include <Wt/WContainerWidget>

#include "InterSpec/SpecFileQuery.h"
#include "InterSpec/SpecFileQueryDbCache.h"

class InterSpec;

class PopupDivMenu;
class SpecMeasManager;
class ResultTableModel;
class RowStretchTreeView;

namespace SpecUtils
{
  class SpecFile;
}

namespace Wt
{
  class WText;
  class WAnchor;
  class WSpinBox;
  class WCheckBox;
  class WLineEdit;
  class WFileUpload;
  class WPushButton;
  class WApplication;
}

class SpecFileQueryDbCache;

/* Implimetation ideas
 -When user "hovers" over a row for abit, show a spectrum preview using the D3.js plotting (try to re-use logic from somewhere to what samples to plot).  Or similarish thing for when row is double-clicked.  (this all would be useful other places too...)
 -Help text! (could use this as an excuse to re-organize the help XML files).
*/

//Using recursive_directory_iterator saves the initial time to do the recursive_ls
// which can be substantial (for the large test directory on my mac like ~10
// seconds, and Windows like 30 or so)
#define USE_DIRECTORY_ITERATOR_METHOD 1


class SpecFileQueryWidget : public Wt::WContainerWidget
{
public:
  SpecFileQueryWidget( InterSpec *viewer, Wt::WContainerWidget *parent = 0 );
  virtual ~SpecFileQueryWidget();
  
  
  
  
  
protected:
  void init();
  
  /** Counts the number of candidate files to query, in the given (optionally
      recursive) source directory, filtering them by file max size (MB) and file
      extension preferences.
      Results are are posted to WServer using the specified WApplication
      sessionid, calling updateNumberFilesInGui to actually update the GUI.
      This function can run in any thread.
      \param widgetdeleted should be a copy of #m_widgetDeleted
   */
  static void updateNumberFiles( const std::string srcdir, const bool recursive,
                          const bool extfilter, const size_t maxsize,
                          SpecFileQueryWidget *querywidget,
                          const std::string sessionid,
                          std::shared_ptr< std::atomic<bool> > widgetdeleted,
                          std::shared_ptr<SpecFileQueryDbCache> database );
  
  /** Updates the GUI for the found number of candidate files, if the source
      directory, recursice and extension filter options are still all the same.
      This function must be run from the main WApplication thread for this
      instance.
      \param widgetdeleted should be a copy of #m_widgetDeleted to determine
      if this SpecFileQueryWidget had gotten deleted while the file search was
      happening.
   */
  static void updateNumberFilesInGui( const size_t nfiles,
                                      const bool completed,
                                      const std::string srcdir,
                                      const bool recursive,
                                      const bool extfilter,
                                      SpecFileQueryWidget *querywidget,
                               std::shared_ptr< std::atomic<bool> > widgetdeleted );
  
  void setResultFieldVisibility( const SpecFileQuery::FileDataField field, const bool visible );

#if( BUILD_AS_ELECTRON_APP )
  void newElectronPathSelected( std::string path );
#elif( BUILD_AS_OSX_APP )
  ///Called when the user selects a new path. \basePathChanged is called also
  /// when filter options change.
  void newMacOsPathSelected();
#endif
  
  void basePathChanged();
  void setResultsStale();
  void doCacheChanged();
  void doPersistCacheChanged();
  void queryChangedCallback( const std::string &queryJson );
  void searchRequestedCallback( const std::string &queryJson );
  
  void finishUpdate( std::shared_ptr< std::vector< std::vector<std::string> > > result,
                     const std::string description,
                     const double wallSeconds,
                     const bool wasCanceled,
                     std::shared_ptr< std::atomic<bool> > widgetDeleted );
  
  void cancelUpdate();
  void updateSearchStatus( const size_t nfilestotal, const size_t nfileschecked,
                           const std::string specialmsg,
                           std::shared_ptr< std::vector< std::vector<std::string> > > results,
                           std::shared_ptr< std::atomic<bool> > widgetDeleted );
  
  void doSearch( const std::string basedir, unsigned long options,
                 const size_t maxsize,
                 const SpecFileQuery::SpecLogicTest query,
                 const std::string sesssionID,
                 std::shared_ptr< SpecFileQueryDbCache > database,
                 std::shared_ptr< std::atomic<bool> > stopUpdate,
                 std::shared_ptr< std::atomic<bool> > widgetDeleted );
  
  void selectionChanged();
  void loadSelected();
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  void openSelectedFilesParentDir();
#endif
  
  std::string baseDirectory();
  
  /** Converts the JSON from the client side to a SpecFileQuery::SpecLogicTest.
   Throws std::exception if the query is not valid or other error encountered.
   */
  static SpecFileQuery::SpecLogicTest queryFromJson( const std::string &json,
                                                    const std::vector<EventXmlFilterInfo> &eventXml );
  
protected:
  
  std::string prepareEventXmlFilters();
  std::vector<EventXmlFilterInfo> m_eventXmlFilters;
  
  
  enum PrefilterOptions
  {
    FilterByFilename = 0x1,
    FilterDuplicates = 0x2,
    SearchRecursive  = 0x4
  };//enum PrefilterOptions
  
  std::shared_ptr< std::atomic<bool> > m_stopUpdate;
  std::shared_ptr< std::atomic<bool> > m_widgetDeleted;
  
  Wt::WApplication * const m_app;
  InterSpec * const m_viewer;
  
  Wt::WContainerWidget *m_conditions;
  Wt::WPushButton *m_update;
  Wt::WPushButton *m_cancelUpdate;
  Wt::WPushButton *m_loadSelectedFile;
  Wt::WPushButton *m_openSelectedDir;
  
  ResultTableModel *m_resultmodel;
  RowStretchTreeView *m_resultview;
  
#if( BUILD_AS_ELECTRON_APP )
  std::unique_ptr<Wt::JSignal<std::string>> m_pathSelectedSignal;
  std::string m_basePath;
#elif( BUILD_AS_OSX_APP )
  //TODO: currently using a WFileUpload to browse for a file - should use custom JS/html to get rid of some issues
  std::string m_basePath;
  Wt::WFileUpload *m_baseLocation;
#else
  Wt::WLineEdit *m_baseLocation;
#endif
  
  Wt::WCheckBox *m_recursive;
  Wt::WCheckBox *m_filterByExtension;
  Wt::WSpinBox *m_maxFileSize;
  Wt::WCheckBox *m_filterUnique;
  
  Wt::WCheckBox *m_cacheParseResults;
  Wt::WCheckBox *m_persistCacheResults;
  
  Wt::WPushButton *m_optionsBtn;
  PopupDivMenu *m_optionsMenu;
  
  Wt::WText *m_numberFiles;
  Wt::WText *m_numberResults;
  
#if( BUILD_AS_OSX_APP )
  Wt::WAnchor *m_csv;
#else
  Wt::WPushButton *m_csv;
#endif
  
  std::string m_queryJson;
  
  Wt::JSignal<std::string> m_queryChanged;
  Wt::JSignal<std::string> m_searchRequested;
  
  std::map<std::string,std::shared_ptr<SpecFileQueryDbCache>> m_path_caches;
};//class SpecFileQueryWidget


#endif
