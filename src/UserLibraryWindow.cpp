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

// The user library exists only where there is a user data folder (not in web deployments).
#if( USE_USER_LIBRARY )

#include <set>
#include <map>
#include <mutex>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <filesystem>

#include <Wt/WMenu.h>
#include <Wt/WText.h>
#include <Wt/Utils.h>
#include <Wt/WLabel.h>
#include <Wt/WAnchor.h>
#include <Wt/WCheckBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WMenuItem.h>
#include <Wt/WResource.h>
#include <Wt/WTextArea.h>
#include <Wt/WComboBox.h>
#include <Wt/WGridLayout.h>
#include <Wt/WFileUpload.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/Http/Response.h>
#include <Wt/WStackedWidget.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfLibrary.h"
#include "InterSpec/UserLibrary.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RefSpectraWidget.h"
#include "InterSpec/UserLibraryWindow.h"
#include "InterSpec/UserLibraryPackage.h"
#include "InterSpec/UserLibraryImportDialog.h"

#if( BUILD_AS_OSX_APP )
#include "target/macos/macOsUtils.h"
#endif

using namespace std;
using namespace Wt;


struct UserLibraryWindow::PackageInfo
{
  UserLibrary::Package pkg;
  map<string,size_t> numByType;
  size_t numNeedsNewer = 0;
  size_t numUnreadable = 0;
  size_t numModified = 0;
};//struct PackageInfo


struct UserLibraryWindow::ItemRow
{
  UserLibrary::Item item;
  UserLibrary::Validation validation;
  bool chosen = false;     ///< a recorded choice names this item
  bool hasChoice = false;  ///< a choice is recorded for one of its names
  bool hasRival = false;   ///< another usable item has one of its names
};//struct ItemRow


struct UserLibraryWindow::Candidate
{
  UserLibrary::Item item;
  string sourceKey;
  string packagePath;
  vector<string> includes;  ///< file names of support files (e.g. template includes)
  UserLibrary::Status status = UserLibrary::Status::Ok;

  /** For a folder the user added outside the user data folder (Detector Efficiency Locations, or a
   reference-spectra folder): that folder, and the folder package paths are taken relative to.
   */
  string locationPath, root;

  /** Saved in InterSpec's database: a detector efficiency (loaded, on a worker, when creating), or
   an item whose file content is made here (a colour theme).
   */
  bool fromDatabase = false;
  DrfLibrary::Entry dbEntry;
  shared_ptr<const string> generatedData;
  string generatedName;
};//struct Candidate


namespace
{
  /** Serves a created package file. */
  class PackageResource : public Wt::WResource
  {
  public:
    PackageResource( shared_ptr<const string> data, const string &filename )
      : WResource(), m_data( data )
    {
      suggestFileName( WString::fromUTF8( filename ), ContentDisposition::Attachment );
    }

    virtual ~PackageResource()
    {
      beingDeleted();
    }

    virtual void handleRequest( const Wt::Http::Request &, Wt::Http::Response &response ) override
    {
      // Not "application/zip": browsers and download managers may unzip it, or add ".zip".
      response.setMimeType( "application/octet-stream" );
      if( m_data )
        response.out().write( m_data->data(), static_cast<std::streamsize>( m_data->size() ) );
    }

  private:
    const shared_ptr<const string> m_data;
  };//class PackageResource


  WText *add_text( WContainerWidget *parent, const WString &txt, const char *style_class )
  {
    WText *t = parent->addNew<WText>( txt, TextFormat::Plain );
    t->setInline( false );
    if( style_class )
      t->addStyleClass( style_class );
    return t;
  }


  WString names_label( const UserLibrary::Item &item )
  {
    if( item.names.empty() )
      return WString::fromUTF8( SpecUtils::filename( item.found.path ) );

    string joined;
    for( size_t i = 0; (i < item.names.size()) && (i < 3); ++i )
      joined += (i ? ", " : "") + item.names[i];
    WString answer = WString::fromUTF8( joined );
    if( item.names.size() > 3 )
    {
      answer += WString::fromUTF8( " " );
      answer += WString::tr( "ul-item-more-names" ).arg( static_cast<int>( item.names.size() - 3 ) );
    }
    return answer;
  }//names_label(...)


  /** Whether "YYYY-MM-DD..." `date` is before today (UTC); false if not that format. */
  std::filesystem::path to_fs_path( const string &utf8 )
  {
#ifdef _WIN32
    return std::filesystem::path( SpecUtils::convert_from_utf8_to_utf16( utf8 ) );
#else
    return std::filesystem::path( utf8 );
#endif
  }


  /** Reference-spectra folders the user added in the Reference Spectra tool (outside the user data
   folder), as RefSpectraWidget finds them.  Reads a preference: session thread.
   */
  vector<string> added_reference_spectra_dirs( InterSpec *viewer )
  {
    vector<string> answer;
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
#if( BUILD_AS_OSX_APP )
    // Access is started once per process, and left on (a package being made may read the files);
    //  starting it on every refresh would leak a kernel resource each time.
    static std::mutex s_mutex;
    static std::map<std::string,std::string> s_started;  //bookmark key -> resolved path
    std::lock_guard<std::mutex> lock( s_mutex );
    for( const std::string &key : macOsUtils::getSecurityScopedBookmarkKeys( "SecBookmark_RefSpectra_" ) )
    {
      std::map<std::string,std::string>::const_iterator pos = s_started.find( key );
      if( pos == end(s_started) )
      {
        std::string path;
        if( macOsUtils::resolveAndStartAccessingSecurityScopedBookmark( key, path ) )
          pos = s_started.emplace( key, path ).first;
      }
      if( (pos != end(s_started)) && SpecUtils::is_directory( pos->second ) )
        answer.push_back( pos->second );
    }
#else
    try
    {
      const string dirs = UserPreferences::preferenceValue<string>( "RefSpectraDirs", viewer );
      vector<string> parts;
      SpecUtils::split( parts, dirs, ";" );
      for( const string &dir : parts )
      {
        if( SpecUtils::is_directory( dir ) )
          answer.push_back( dir );
      }
    }catch( std::exception & )
    {
    }
#endif
#endif
    return answer;
  }//added_reference_spectra_dirs(...)


  /** The InterSpec versions a package can say it needs: 1.1 (the first with the user library), each
   minor release of this major version since, and this exact version (the default).  Empty before 1.2:
   until then every InterSpec with a user library reads every package, so packages need not say.
   */
  vector<string> min_version_choices()
  {
    const string current = UserLibrary::interSpecVersion();
    bool numeric = false;
    if( (UserLibrary::compareVersions( current, "1.2", &numeric ) < 0) || !numeric )
      return {};

    vector<string> parts;
    SpecUtils::split( parts, current, "." );
    const unsigned long major = std::stoul( parts.at( 0 ) );
    const unsigned long minor = (parts.size() > 1) ? std::stoul( parts[1] ) : 0;

    vector<string> answer{ "1.1" };
    for( unsigned long m = ((major == 1) ? 2 : 0); m <= minor; ++m )
      answer.push_back( std::to_string( major ) + "." + std::to_string( m ) );
    if( UserLibrary::compareVersions( current, answer.back() ) != 0 )
      answer.push_back( current );
    return answer;
  }//min_version_choices()
}//namespace


Wt::WString UserLibraryWindow::typeTitle( const std::string &typeId )
{
  const UserLibrary::AssetType * const type = UserLibrary::assetType( typeId );
  return type ? WString::tr( type->titleKey ) : WString::fromUTF8( typeId );
}


Wt::WString UserLibraryWindow::sourceLabel( const UserLibrary::Source &source )
{
  switch( source.kind )
  {
    case UserLibrary::SourceKind::BuiltIn: return WString::tr( "ul-source-builtin" );
    case UserLibrary::SourceKind::User:    return WString::tr( "ul-source-user" );
    case UserLibrary::SourceKind::Package: return WString::tr( "ul-source-package" ).arg( WString::fromUTF8( source.label ) );
  }
  return WString();
}//sourceLabel(...)


Wt::WString UserLibraryWindow::whereUsedLabel( const std::string &typeId )
{
  if( typeId == "drf" )
    return WString::tr( "ul-show-drf" );
#if( USE_REL_ACT_TOOL )
  if( typeId == "rel_act_preset" )
    return WString::tr( "ul-show-rel-act" );
#endif
  if( typeId == "act_shield_template" )
    return WString::tr( "ul-show-act-shield" );
  if( typeId == "act_shield_setup" )
    return WString::tr( "ul-show-act-shield-setup" );
  if( typeId == "ref_spectra" )
    return WString::tr( "ul-show-ref-spectra" );
  if( typeId == "color_theme" )
    return WString::tr( "ul-show-color-theme" );
  if( typeId == "source_lib" )
    return WString::tr( "ul-show-make-drf" );
  if( typeId == "nuc_search_category" )
    return WString::tr( "ul-show-nuc-search" );
#if( USE_LLM_INTERFACE )
  if( typeId == "llm_skill" )
    return WString::tr( "ul-show-llm" );
#endif
  return WString();
}//whereUsedLabel(...)


void UserLibraryWindow::showWhereUsed( const std::string &typeId )
{
  InterSpec * const viewer = InterSpec::instance();
  if( !viewer )
    return;

  if( typeId == "drf" )
    viewer->showDrfSelectWindow();
#if( USE_REL_ACT_TOOL )
  else if( typeId == "rel_act_preset" )
    viewer->relActAutoWindow( true );
#endif
  else if( (typeId == "act_shield_template") || (typeId == "act_shield_setup") )
    viewer->shieldingSourceFit( true );
  else if( typeId == "ref_spectra" )
    RefSpectraDialog::createDialog( RefSpectraInitialBehaviour::Default, SpecUtils::SpectrumType::Foreground );
  else if( typeId == "color_theme" )
    viewer->showColorThemeWindow();
  else if( typeId == "source_lib" )
    viewer->showMakeDrfWindow();
  else if( typeId == "nuc_search_category" )
    viewer->showNuclideSearch();
#if( USE_LLM_INTERFACE )
  else if( typeId == "llm_skill" )
    viewer->createLlmTool();
#endif
}//showWhereUsed(...)


UserLibraryWindow::UserLibraryWindow( InterSpec *viewer )
  : AuxWindow( WString::tr("ul-window-title"),
               (AuxWindowProperties::TabletNotFullScreen | AuxWindowProperties::SetCloseable
                | AuxWindowProperties::EnableResize) ),
    m_interspec( viewer ),
    m_menu( nullptr ),
    m_stack( nullptr ),
    m_packagesList( nullptr ),
    m_packagesRequest( 0 ),
    m_typeSelect( nullptr ),
    m_itemSearch( nullptr ),
    m_showMe( nullptr ),
    m_itemsList( nullptr ),
    m_rows(),
    m_itemsRequest( 0 ),
    m_startFrom( nullptr ),
    m_minVersion( nullptr ),
    m_candidateSearch( nullptr ),
    m_numSelected( nullptr ),
    m_candidatesList( nullptr ),
    m_createStatus( nullptr ),
    m_createBtn( nullptr ),
    m_candidates(),
    m_candidatesRequest( 0 ),
    m_creating( false )
{
  assert( m_interspec );
  m_interspec->useMessageResourceBundle( "UserLibrary" );
  wApp->useStyleSheet( "InterSpec_resources/UserLibrary.css" );
  addStyleClass( "UserLibWindow" );

  const int w = m_interspec->renderedWidth(), h = m_interspec->renderedHeight();
  if( (w > 0) && (h > 0) )
    resizeWindow( std::min( std::max( w/2, 760 ), (9*w)/10 ), std::min( std::max( (6*h)/7, 480 ), (9*h)/10 ) );
  else
    resizeScaledWindow( 0.8, 0.85 );

  // Tab bar (a WMenu styled as tabs) above the pages; only the pages scroll.
  auto tab_bar = std::make_unique<WContainerWidget>();
  tab_bar->addStyleClass( "UserLibTabs" );
  auto stack = std::make_unique<WStackedWidget>();
  m_stack = stack.get();
  m_stack->addStyleClass( "UserLibStack" );
  m_menu = tab_bar->addNew<WMenu>( m_stack );
  m_menu->addStyleClass( "Wt-tabs" );  //InterSpec.css styles a `ul.Wt-tabs` as tabs

  WGridLayout *layout = stretcher();
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  layout->addWidget( std::move( tab_bar ), 0, 0 );
  layout->addWidget( std::move( stack ), 1, 0 );
  layout->setRowStretch( 1, 1 );

  const bool have_library = UserLibrary::available();

  // ---- Packages ----
  {
    auto page = std::make_unique<WContainerWidget>();
    page->addStyleClass( "UserLibPage" );
    page->setOverflow( Overflow::Auto, Orientation::Vertical );

    WContainerWidget *intro = page->addNew<WContainerWidget>();
    intro->addStyleClass( "UserLibIntro" );
    if( have_library )
    {
      intro->addNew<WText>( WString::tr( "ul-packages-intro" ) );
      WPushButton *install = intro->addNew<WPushButton>( WString::tr( "ul-install-btn" ) );
      install->addStyleClass( "LinkBtn" );
      install->clicked().connect( this, &UserLibraryWindow::startUpload );
    }else
    {
      intro->addNew<WText>( WString::tr( "ul-unavailable" ) );
    }

    m_packagesList = page->addNew<WContainerWidget>();
    m_packagesList->addStyleClass( "UserLibPackages" );
    m_menu->addItem( WString::tr( "ul-tab-packages" ), std::move( page ) );
  }

  // ---- By type ----
  {
    auto page = std::make_unique<WContainerWidget>();
    page->addStyleClass( "UserLibPage" );
    page->setOverflow( Overflow::Auto, Orientation::Vertical );

    WContainerWidget *controls = page->addNew<WContainerWidget>();
    controls->addStyleClass( "UserLibControls" );
    WLabel *type_label = controls->addNew<WLabel>( WString::tr( "ul-items-type" ) );
    m_typeSelect = controls->addNew<WComboBox>();
    type_label->setBuddy( m_typeSelect );
    for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
      m_typeSelect->addItem( WString::tr( type.titleKey ) );
    m_typeSelect->setCurrentIndex( 0 );
    m_typeSelect->changed().connect( this, &UserLibraryWindow::refreshItems );

    m_itemSearch = controls->addNew<WLineEdit>();
    m_itemSearch->setPlaceholderText( WString::tr( "ul-items-search" ) );
    m_itemSearch->addStyleClass( "UserLibSearch" );
    m_itemSearch->textInput().connect( this, &UserLibraryWindow::renderItems );

    m_showMe = controls->addNew<WPushButton>();
    m_showMe->addStyleClass( "LinkBtn" );
    m_showMe->clicked().connect( [this](){
      const int index = m_typeSelect->currentIndex();
      if( (index < 0) || (index >= static_cast<int>( UserLibrary::assetTypes().size() )) )
        return;
      // Out of the way of the tool it opens (closing deletes this window, later).
      const string type_id = UserLibrary::assetTypes()[index].id;
      hide();
      UserLibraryWindow::showWhereUsed( type_id );
    } );

    m_itemsList = page->addNew<WContainerWidget>();
    m_itemsList->addStyleClass( "UserLibItems" );
    m_menu->addItem( WString::tr( "ul-tab-items" ), std::move( page ) );
  }

  // ---- Create package ----
  if( have_library )
  {
    auto page = std::make_unique<WContainerWidget>();
    page->addStyleClass( "UserLibPage UserLibCreate" );
    page->setOverflow( Overflow::Auto, Orientation::Vertical );

    add_text( page.get(), WString::tr( "ul-create-intro" ), "UserLibIntro" );

    WContainerWidget *start_row = page->addNew<WContainerWidget>();
    start_row->addStyleClass( "UserLibControls" );
    WLabel *start_label = start_row->addNew<WLabel>( WString::tr( "ul-create-start" ) );
    m_startFrom = start_row->addNew<WComboBox>();
    start_label->setBuddy( m_startFrom );
    m_startFrom->changed().connect( this, &UserLibraryWindow::startFromChanged );

    add_text( page.get(), WString::tr( "ul-create-details" ), "UserLibSubhead" );
    WContainerWidget *form = page->addNew<WContainerWidget>();
    form->addStyleClass( "UserLibForm" );

    const pair<const char *,const char *> line_fields[] = {
      { "name", "ul-field-name" }, { "version", "ul-field-version" }, { "author", "ul-field-author" },
      { "organization", "ul-field-organization" }, { "contact", "ul-field-contact" },
      { "description", "ul-field-description" }, { "handling", "ul-field-handling" },
      { "review_by", "ul-field-review-by" }
    };
    for( const pair<const char *,const char *> &f : line_fields )
    {
      WLabel *label = form->addNew<WLabel>( WString::tr( f.second ) );
      WLineEdit *edit = form->addNew<WLineEdit>();
      label->setBuddy( edit );
      m_fields[f.first] = edit;
    }
    m_fields["review_by"]->setPlaceholderText( "YYYY-MM-DD" );

    // The InterSpec a package needs; only offered once there is more than one user library version.
    const vector<string> min_versions = min_version_choices();
    if( !min_versions.empty() )
    {
      WLabel *label = form->addNew<WLabel>( WString::tr( "ul-field-min-version" ) );
      m_minVersion = form->addNew<WComboBox>();
      label->setBuddy( m_minVersion );
      for( const string &v : min_versions )
        m_minVersion->addItem( WString::fromUTF8( v ) );
      m_minVersion->setCurrentIndex( m_minVersion->count() - 1 );
    }

    const pair<const char *,const char *> area_fields[] = {
      { "notes", "ul-field-notes" }, { "changes", "ul-field-changes" }
    };
    for( const pair<const char *,const char *> &f : area_fields )
    {
      WLabel *label = form->addNew<WLabel>( WString::tr( f.second ) );
      WTextArea *area = form->addNew<WTextArea>();
      area->setRows( 3 );
      label->setBuddy( area );
      m_areas[f.first] = area;
    }

    add_text( page.get(), WString::tr( "ul-create-items" ), "UserLibSubhead" );
    WContainerWidget *item_controls = page->addNew<WContainerWidget>();
    item_controls->addStyleClass( "UserLibControls" );
    m_candidateSearch = item_controls->addNew<WLineEdit>();
    m_candidateSearch->setPlaceholderText( WString::tr( "ul-items-search" ) );
    m_candidateSearch->addStyleClass( "UserLibSearch" );
    m_candidateSearch->textInput().connect( this, &UserLibraryWindow::renderCandidates );
    m_numSelected = item_controls->addNew<WText>( "", TextFormat::Plain );
    m_numSelected->addStyleClass( "UserLibPkgMeta" );

    m_candidatesList = page->addNew<WContainerWidget>();
    m_candidatesList->addStyleClass( "UserLibItems" );

    WContainerWidget *create_row = page->addNew<WContainerWidget>();
    create_row->addStyleClass( "UserLibCreateRow" );
    m_createBtn = create_row->addNew<WPushButton>( WString::tr( "ul-create-btn" ) );
    WidgetUtils::applyButtonRole( m_createBtn, WidgetUtils::ButtonRole::Affirm );
    m_createBtn->clicked().connect( this, &UserLibraryWindow::startCreate );

    m_createStatus = page->addNew<WContainerWidget>();
    m_createStatus->addStyleClass( "UserLibCreateStatus" );

    m_menu->addItem( WString::tr( "ul-tab-create" ), std::move( page ) );
  }//if( have_library )

  m_menu->itemSelected().connect( [this]( WMenuItem * ){
    if( m_menu->currentIndex() == 1 )
      renderItems();
  } );

  AuxWindow::addHelpInFooter( footer(), "user-library" );
  WPushButton *close = addCloseButtonToFooter();
  close->clicked().connect( this, &AuxWindow::hide );
  rejectWhenEscapePressed();

  m_interspec->userLibraryChanged().connect( this, &UserLibraryWindow::refresh );

  updateStartFromOptions();
  refresh();

  show();
  centerWindow();
}//UserLibraryWindow constructor


UserLibraryWindow::~UserLibraryWindow()
{
}


void UserLibraryWindow::refresh()
{
  refreshPackages();
  refreshItems();
  if( m_candidatesList )
    refreshCandidates();
}//refresh()


void UserLibraryWindow::libraryChangedHere()
{
  // Every open user library view, and tools that list items, re-read (this window included).
  m_interspec->userLibraryChanged().emit();
}


// -------------------------------------------------------------------------------------------------
//  Packages page
// -------------------------------------------------------------------------------------------------

void UserLibraryWindow::refreshPackages()
{
  if( !UserLibrary::available() )
    return;

  const uint64_t request = ++m_packagesRequest;
  if( m_packagesList->count() == 0 )
    add_text( m_packagesList, WString::tr( "ul-loading" ), "UserLibBusy" );

  const shared_ptr<vector<PackageInfo>> result = make_shared<vector<PackageInfo>>();
  const WidgetUtils::WidgetHandle self( this );
  AppUtils::run_on_ioservice( [result](){
      const vector<UserLibrary::Package> packages = UserLibrary::installedPackages();
      for( const UserLibrary::Package &pkg : packages )
      {
        PackageInfo info;
        info.pkg = pkg;
        for( const pair<const string,string> &f : pkg.files )
        {
          string sha;
          try
          {
            sha = UserLibrary::sha256OfFile( SpecUtils::append_path( pkg.path, f.first ) );
          }catch( std::exception & )
          {
          }
          info.numModified += (sha != f.second);
        }
        result->push_back( info );
      }//for( loop over packages )

      // Item counts and statuses, re-checked each time (a newer InterSpec may now read them).
      for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
      {
        for( const UserLibrary::Item &item : UserLibrary::items( type.id, true ) )
        {
          if( item.source.kind != UserLibrary::SourceKind::Package )
            continue;
          for( PackageInfo &info : *result )
          {
            if( info.pkg.id != item.source.packageId )
              continue;
            info.numByType[type.id] += 1;
            const UserLibrary::Validation v = UserLibrary::validateItem( type.id, item.found );
            info.numNeedsNewer += (item.needsNewer || (v.status == UserLibrary::Status::NeedsNewer));
            info.numUnreadable += (v.status == UserLibrary::Status::Unreadable);
          }
        }
      }//for( loop over types )
    },
    [self, result, request](){
      UserLibraryWindow * const window = self.resolve_as<UserLibraryWindow>();
      if( window && (window->m_packagesRequest == request) )
        window->packagesListed( result );
      if( wApp )
        wApp->triggerUpdate();
    } );
}//refreshPackages()


void UserLibraryWindow::packagesListed( std::shared_ptr<const std::vector<PackageInfo>> packages )
{
  m_packagesList->clear();
  if( !packages || packages->empty() )
  {
    add_text( m_packagesList, WString::tr( "ul-no-packages" ), "UserLibNone" );
    return;
  }

  for( const PackageInfo &info : *packages )
  {
    const UserLibrary::Package &pkg = info.pkg;
    const UserLibrary::Manifest &m = pkg.manifest;

    WContainerWidget *card = m_packagesList->addNew<WContainerWidget>();
    card->addStyleClass( "UserLibCard" );

    WContainerWidget *title_line = card->addNew<WContainerWidget>();
    title_line->addStyleClass( "UserLibPkgTitleLine" );
    title_line->addNew<WText>( WString::fromUTF8( m.name ), TextFormat::Plain )->addStyleClass( "UserLibPkgName" );
    title_line->addNew<WText>( WString::tr( "ul-pkg-version" ).arg( WString::fromUTF8( m.version ) ), TextFormat::Plain )
      ->addStyleClass( "UserLibPkgVersion" );

    if( !m.handling.empty() )
      add_text( card, WString::tr( "ul-pkg-handling" ).arg( WString::fromUTF8( m.handling ) ), "UserLibHandling" );

    if( !m.reviewBy.empty() )
    {
      if( m.reviewDue() )
        add_text( card, WString::tr( "ul-pkg-review-passed" ).arg( WString::fromUTF8( m.reviewBy ) ), "UserLibWarn" );
      else
        add_text( card, WString::tr( "ul-pkg-review-by" ).arg( WString::fromUTF8( m.reviewBy ) ), "UserLibPkgMeta" );
    }

    string by;
    for( const string &s : { m.author, m.organization } )
      if( !s.empty() )
        by += (by.empty() ? "" : " · ") + s;
    if( !by.empty() )
      add_text( card, WString::tr( "ul-pkg-by" ).arg( WString::fromUTF8( by ) ), "UserLibPkgMeta" );
    if( !m.contact.empty() )
      add_text( card, WString::tr( "ul-pkg-contact" ).arg( WString::fromUTF8( m.contact ) ), "UserLibPkgMeta" );

    string dates;
    if( !m.created.empty() )
      dates = WString::tr( "ul-pkg-created" ).arg( WString::fromUTF8( m.created.substr( 0, 10 ) ) ).toUTF8();
    const string installed = pkg.installs.empty() ? pkg.firstInstalled : pkg.installs.back().installed;
    if( !installed.empty() )
      dates += (dates.empty() ? "" : " · ")
               + WString::tr( "ul-pkg-installed" ).arg( WString::fromUTF8( installed.substr( 0, 10 ) ) ).toUTF8();
    if( !dates.empty() )
      add_text( card, WString::fromUTF8( dates ), "UserLibPkgMeta" );

    if( !m.description.empty() )
      add_text( card, WString::fromUTF8( m.description ), "UserLibPkgDesc" );

    // What it holds.
    string counts;
    for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
    {
      const map<string,size_t>::const_iterator pos = info.numByType.find( type.id );
      if( pos == end(info.numByType) )
        continue;
      counts += (counts.empty() ? "" : ", ")
                + WString::tr( "ul-pkg-count-type" ).arg( static_cast<int>( pos->second ) )
                    .arg( WString::tr( type.titleKey ) ).toUTF8();
    }
    if( !counts.empty() )
      add_text( card, WString::tr( "ul-pkg-counts" ).arg( WString::fromUTF8( counts ) ), "UserLibPkgMeta" );

    if( !pkg.recorded )
      add_text( card, WString::tr( "ul-pkg-not-recorded" ), "UserLibWarn" );
    if( info.numModified )
      add_text( card, WString::tr( "ul-pkg-modified" ).arg( static_cast<int>( info.numModified ) ), "UserLibWarn" );
    if( info.numNeedsNewer )
      add_text( card, WString::tr( "ul-pkg-needs-newer" ).arg( static_cast<int>( info.numNeedsNewer ) ), "UserLibWarn" );
    if( info.numUnreadable )
      add_text( card, WString::tr( "ul-pkg-unreadable" ).arg( static_cast<int>( info.numUnreadable ) ), "UserLibError" );

    if( !m.notes.empty() )
    {
      add_text( card, WString::tr( "ul-pkg-notes" ), "UserLibSubhead" );
      add_text( card, WString::fromUTF8( m.notes ), "UserLibPkgNotes" );
    }

    // Each version's change notes, newest first: those installed here, and earlier ones the manifest carries.
    vector<UserLibrary::ManifestPreviousVersion> history;
    set<string> versions_seen;
    for( auto it = pkg.installs.rbegin(); it != pkg.installs.rend(); ++it )
    {
      if( versions_seen.insert( it->manifest.version ).second )
        history.push_back( { it->manifest.version, it->manifest.created, it->manifest.changes } );
    }
    if( versions_seen.insert( m.version ).second )
      history.insert( begin(history), { m.version, m.created, m.changes } );
    for( const UserLibrary::ManifestPreviousVersion &v : m.previousVersions )
    {
      if( versions_seen.insert( v.version ).second )
        history.push_back( v );
    }
    bool any_changes = false;
    for( const UserLibrary::ManifestPreviousVersion &v : history )
      any_changes = any_changes || !v.changes.empty();
    if( any_changes )
    {
      WContainerWidget *changes = card->addNew<WContainerWidget>();
      changes->addStyleClass( "UserLibChanges" );
      add_text( changes, WString::tr( "ul-pkg-history" ), "UserLibSubhead" );
      for( const UserLibrary::ManifestPreviousVersion &v : history )
      {
        if( v.changes.empty() )
          continue;
        const string when = v.created.empty() ? string() : (" (" + v.created.substr( 0, 10 ) + ")");
        add_text( changes, WString::tr( "ul-pkg-history-entry" ).arg( WString::fromUTF8( v.version ) )
                             .arg( WString::fromUTF8( when ) ), "UserLibHistVersion" );
        add_text( changes, WString::fromUTF8( v.changes ), "UserLibHistText" );
      }
    }//if( any_changes )

    WContainerWidget *actions = card->addNew<WContainerWidget>();
    actions->addStyleClass( "UserLibActions" );
    WPushButton *recheck = actions->addNew<WPushButton>( WString::tr( "ul-recheck" ) );
    recheck->addStyleClass( "LinkBtn" );
    recheck->clicked().connect( this, &UserLibraryWindow::refreshPackages );
    WPushButton *uninstall = actions->addNew<WPushButton>( WString::tr( "ul-uninstall" ) );
    uninstall->addStyleClass( "LinkBtn UserLibDanger" );
    const string id = pkg.id, name = m.name;
    uninstall->clicked().connect( [this, id, name](){ confirmUninstall( id, name ); } );
  }//for( loop over packages )
}//packagesListed(...)


void UserLibraryWindow::confirmUninstall( const std::string &packageId, const std::string &name )
{
  SimpleDialog * const dialog = SimpleDialog::make( WString::tr( "ul-uninstall-title" ),
                         WString::tr( "ul-uninstall-msg" ).arg( WString::fromUTF8( Wt::Utils::htmlEncode( name ) ) ) );
  WPushButton *yes = dialog->addButton( WString::tr( "ul-uninstall-confirm" ), WidgetUtils::ButtonRole::Destructive );
  dialog->addButton( WString::tr( "Cancel" ), WidgetUtils::ButtonRole::Dismiss );

  yes->clicked().connect( [packageId, name](){
    try
    {
      UserLibraryPackage::uninstall( packageId );
      passMessage( WString::tr( "ul-uninstalled" ).arg( WString::fromUTF8( Wt::Utils::htmlEncode( name ) ) ),
                   WarningWidget::WarningMsgInfo );
    }catch( std::exception &e )
    {
      passMessage( WString::tr( "ul-uninstall-error" ).arg( WString::fromUTF8( Wt::Utils::htmlEncode( string( e.what() ) ) ) ),
                   WarningWidget::WarningMsgHigh );
    }

    // Re-lists this window too.
    InterSpec * const viewer = InterSpec::instance();
    if( viewer )
      viewer->userLibraryChanged().emit();
  } );
}//confirmUninstall(...)


void UserLibraryWindow::startUpload()
{
  AuxWindow *window = AuxWindow::make( WString::tr( "ul-install-btn" ),
                                       (AuxWindowProperties::IsModal | AuxWindowProperties::TabletNotFullScreen
                                        | AuxWindowProperties::DisableCollapse | AuxWindowProperties::SetCloseable) );
  window->finished().connect( [window](){ AuxWindow::deleteAuxWindow( window ); } );

  WFileUpload *upload = window->contents()->addNew<WFileUpload>();
  upload->setInline( false );
  upload->setFilters( string( UserLibraryPackage::FileExtension ) + ",.zip" );
  upload->changed().connect( upload, &WFileUpload::upload );
  upload->uploaded().connect( [window, upload](){
    const string spool = upload->spoolFileName();
    const string name = upload->clientFileName().toUTF8();
    UserLibraryImportDialog::start( spool, name );
    window->hide();
  } );
  upload->fileTooLarge().connect( [window, upload]( ::int64_t ){
    upload->hide();
    window->contents()->addNew<WText>( WString::tr( "ul-upload-too-large" ), TextFormat::Plain );
  } );

  WPushButton *cancel = window->addCloseButtonToFooter( WString::tr( "Cancel" ) );
  cancel->clicked().connect( window, &AuxWindow::hide );
  window->rejectWhenEscapePressed();
  window->show();
  window->centerWindow();
}//startUpload()


// -------------------------------------------------------------------------------------------------
//  By-type page
// -------------------------------------------------------------------------------------------------

void UserLibraryWindow::refreshItems()
{
  const int index = m_typeSelect->currentIndex();
  const vector<UserLibrary::AssetType> &types = UserLibrary::assetTypes();
  if( (index < 0) || (index >= static_cast<int>( types.size() )) )
    return;
  const string type_id = types[index].id;

  const WString show_label = whereUsedLabel( type_id );
  m_showMe->setText( show_label );
  m_showMe->setHidden( show_label.empty() );

  const uint64_t request = ++m_itemsRequest;
  m_itemsList->clear();
  add_text( m_itemsList, WString::tr( "ul-loading" ), "UserLibBusy" );

  const shared_ptr<vector<ItemRow>> result = make_shared<vector<ItemRow>>();
  const WidgetUtils::WidgetHandle self( this );
  AppUtils::run_on_ioservice( [result, type_id](){
      const UserLibrary::AssetType * const type = UserLibrary::assetType( type_id );
      const vector<UserLibrary::Item> all = UserLibrary::items( type_id, true );
      const map<string,UserLibrary::Choice> choices = UserLibrary::choices();
      for( const UserLibrary::Item &item : all )
      {
        ItemRow row;
        row.item = item;
        row.validation = UserLibrary::validateItem( type_id, item.found );
        for( const string &name : item.names )
        {
          for( const pair<const string,UserLibrary::Choice> &c : choices )
          {
            const string prefix = type_id + ":";
            if( !SpecUtils::starts_with( c.first, prefix.c_str() )
                || !UserLibrary::sameName( *type, c.first.substr( prefix.size() ), name ) )
              continue;
            row.hasChoice = true;
            row.chosen = row.chosen || (c.second.source == item.source.key());
          }

          for( const UserLibrary::Item &other : all )
          {
            if( (&other == &item) || other.disabled || other.needsNewer )
              continue;
            for( const string &n : other.names )
              row.hasRival = row.hasRival || UserLibrary::sameName( *type, n, name );
          }
        }//for( loop over names )
        result->push_back( row );
      }//for( loop over items )
    },
    [self, result, request](){
      UserLibraryWindow * const window = self.resolve_as<UserLibraryWindow>();
      if( window && (window->m_itemsRequest == request) )
        window->itemsListed( request, result );
      if( wApp )
        wApp->triggerUpdate();
    } );
}//refreshItems()


void UserLibraryWindow::itemsListed( const uint64_t, std::shared_ptr<const std::vector<ItemRow>> rows )
{
  m_rows = rows;
  renderItems();
}


void UserLibraryWindow::renderItems()
{
  if( !m_rows )
    return;

  m_itemsList->clear();
  const string query = SpecUtils::trim_copy( m_itemSearch->text().toUTF8() );

  // Kinds read when InterSpec starts: changes here apply after a restart.
  const int type_index = m_typeSelect->currentIndex();
  if( (type_index >= 0) && (type_index < static_cast<int>( UserLibrary::assetTypes().size() ))
      && (UserLibrary::assetTypes()[type_index].apply == UserLibrary::ApplyMode::Restart) )
    add_text( m_itemsList, WString::tr( "ul-items-restart" ), "UserLibWarn" );

  size_t num_shown = 0;
  for( size_t i = 0; i < m_rows->size(); ++i )
  {
    const ItemRow &row = (*m_rows)[i];
    const UserLibrary::Item &item = row.item;

    if( !query.empty() )
    {
      bool match = SpecUtils::icontains( item.relPath, query ) || SpecUtils::icontains( item.source.label, query );
      for( const string &n : item.names )
        match = match || SpecUtils::icontains( n, query );
      if( !match )
        continue;
    }
    ++num_shown;

    WContainerWidget *div = m_itemsList->addNew<WContainerWidget>();
    div->addStyleClass( "UserLibItem" );
    if( item.disabled )
      div->addStyleClass( "Disabled" );

    WContainerWidget *line = div->addNew<WContainerWidget>();
    line->addStyleClass( "UserLibItemLine" );
    line->addNew<WText>( names_label( item ), TextFormat::Plain )->addStyleClass( "UserLibItemName" );

    WString status;
    const char *status_class = "";
    if( item.disabled )
    {
      status = WString::tr( "ul-item-disabled" );
    }else if( item.needsNewer || (row.validation.status == UserLibrary::Status::NeedsNewer) )
    {
      status = WString::tr( "ul-item-needs-newer" );
      status_class = "Warn";
    }else if( row.validation.status == UserLibrary::Status::Unreadable )
    {
      status = WString::tr( "ul-item-unreadable" );
      status_class = "Error";
    }else if( item.inUse )
    {
      status = WString::tr( "ul-item-in-use" );
      if( row.chosen )
        status += WString::fromUTF8( " (" ) + WString::tr( "ul-item-your-choice" ) + WString::fromUTF8( ")" );
      status_class = "Ok";
    }else
    {
      status = WString::tr( "ul-item-not-used" );
      status_class = "Warn";
    }
    WText *badge = line->addNew<WText>( status, TextFormat::Plain );
    badge->addStyleClass( "UserLibBadge" );
    if( status_class[0] )
      badge->addStyleClass( status_class );

    add_text( div, sourceLabel( item.source ), "UserLibPkgMeta" );
    add_text( div, WString::fromUTF8( item.relPath ), "UserLibPath" );
    if( (row.validation.status != UserLibrary::Status::Ok) && !row.validation.message.empty() )
      add_text( div, WString::fromUTF8( row.validation.message ), "UserLibItemNotes" );

    WContainerWidget *actions = div->addNew<WContainerWidget>();
    actions->addStyleClass( "UserLibActions" );

    const bool usable = !item.disabled && !item.needsNewer && (row.validation.status == UserLibrary::Status::Ok);
    const UserLibrary::Item item_copy = item;

    if( usable && row.hasRival && !item.inUse )
    {
      WPushButton *use = actions->addNew<WPushButton>( WString::tr( "ul-item-use-this" ) );
      use->addStyleClass( "LinkBtn" );
      use->clicked().connect( [this, item_copy](){
        try
        {
          for( const string &name : item_copy.names )
            UserLibrary::setChoice( item_copy.typeId, name, item_copy );
        }catch( std::exception &e )
        {
          passMessage( WString::fromUTF8( Wt::Utils::htmlEncode( string( e.what() ) ) ), WarningWidget::WarningMsgHigh );
        }
        libraryChangedHere();
      } );
    }//if( could be the one used )

    if( row.hasChoice )
    {
      WPushButton *reset = actions->addNew<WPushButton>( WString::tr( "ul-item-reset-choice" ) );
      reset->addStyleClass( "LinkBtn" );
      reset->clicked().connect( [this, item_copy](){
        try
        {
          for( const string &name : item_copy.names )
            UserLibrary::clearChoice( item_copy.typeId, name );
        }catch( std::exception &e )
        {
          passMessage( WString::fromUTF8( Wt::Utils::htmlEncode( string( e.what() ) ) ), WarningWidget::WarningMsgHigh );
        }
        libraryChangedHere();
      } );
    }//if( row.hasChoice )

    if( item.source.kind != UserLibrary::SourceKind::BuiltIn )
    {
      WCheckBox *enabled = actions->addNew<WCheckBox>( WString::tr( "ul-item-enabled" ) );
      enabled->setChecked( !item.disabled );
      enabled->changed().connect( [this, item_copy, enabled](){
        try
        {
          UserLibrary::setDisabled( item_copy, !enabled->isChecked() );
        }catch( std::exception &e )
        {
          passMessage( WString::fromUTF8( Wt::Utils::htmlEncode( string( e.what() ) ) ), WarningWidget::WarningMsgHigh );
        }
        libraryChangedHere();
      } );
    }//if( not built-in )

    if( item.source.kind == UserLibrary::SourceKind::User )
    {
      WPushButton *del = actions->addNew<WPushButton>( WString::tr( "ul-item-delete" ) );
      del->addStyleClass( "LinkBtn UserLibDanger" );
      del->clicked().connect( [this, i](){ confirmDelete( i ); } );
    }
  }//for( loop over rows )

  if( !num_shown )
    add_text( m_itemsList, WString::tr( "ul-items-none" ), "UserLibNone" );
}//renderItems()


void UserLibraryWindow::confirmDelete( const size_t index )
{
  if( !m_rows || (index >= m_rows->size()) )
    return;

  const UserLibrary::Item item = (*m_rows)[index].item;
  if( item.source.kind != UserLibrary::SourceKind::User )
    return;

  SimpleDialog * const dialog = SimpleDialog::make( WString::tr( "ul-delete-title" ),
                  WString::tr( "ul-delete-msg" ).arg( WString::fromUTF8( Wt::Utils::htmlEncode( item.relPath ) ) ) );
  WPushButton *yes = dialog->addButton( WString::tr( "ul-item-delete-confirm" ), WidgetUtils::ButtonRole::Destructive );
  dialog->addButton( WString::tr( "Cancel" ), WidgetUtils::ButtonRole::Dismiss );

  yes->clicked().connect( [item](){
    try
    {
      UserLibrary::deleteUserItem( item );
    }catch( std::exception &e )
    {
      passMessage( WString::tr( "ul-delete-error" ).arg( WString::fromUTF8( Wt::Utils::htmlEncode( string( e.what() ) ) ) ),
                   WarningWidget::WarningMsgHigh );
    }
    if( item.typeId == "drf" )
      DrfLibrary::invalidatePath( item.found.path );

    InterSpec * const viewer = InterSpec::instance();
    if( viewer )
      viewer->userLibraryChanged().emit();
  } );
}//confirmDelete(...)


// -------------------------------------------------------------------------------------------------
//  Create page
// -------------------------------------------------------------------------------------------------

void UserLibraryWindow::updateStartFromOptions()
{
  if( !m_startFrom )
    return;

  const string current = ((m_startFrom->currentIndex() > 0)
                          && (m_startFrom->currentIndex() < static_cast<int>( m_startFromIds.size() )))
                         ? m_startFromIds[m_startFrom->currentIndex()] : string();

  m_startFrom->clear();
  m_startFromIds.clear();
  m_startFrom->addItem( WString::tr( "ul-create-new" ) );
  m_startFromIds.push_back( "" );

  int to_select = 0;
  for( const pair<const string,vector<UserLibrary::CreatedVersion>> &p : UserLibrary::createdPackages() )
  {
    if( p.second.empty() )
      continue;
    const UserLibrary::Manifest &last = p.second.back().manifest;
    m_startFrom->addItem( WString::tr( "ul-create-new-version" ).arg( WString::fromUTF8( last.name ) )
                            .arg( WString::fromUTF8( last.version ) ) );
    if( p.first == current )
      to_select = static_cast<int>( m_startFromIds.size() );
    m_startFromIds.push_back( p.first );
  }

  m_startFrom->setCurrentIndex( to_select );
  if( !to_select )
    startFromChanged();
}//updateStartFromOptions()


void UserLibraryWindow::startFromChanged()
{
  const int index = m_startFrom->currentIndex();
  const string id = ((index > 0) && (index < static_cast<int>( m_startFromIds.size() ))) ? m_startFromIds[index] : string();

  m_selected.clear();
  m_itemNotes.clear();
  m_itemOverride.clear();
  m_changedSinceLast.clear();
  m_movedSinceLast.clear();
  for( pair<const string,WLineEdit *> &f : m_fields )
    f.second->setText( "" );
  for( pair<const string,WTextArea *> &a : m_areas )
    a.second->setText( "" );

  if( id.empty() )
  {
    const UserLibrary::Creator creator = UserLibrary::lastCreator();
    m_fields["author"]->setText( WString::fromUTF8( creator.author ) );
    m_fields["organization"]->setText( WString::fromUTF8( creator.organization ) );
    m_fields["contact"]->setText( WString::fromUTF8( creator.contact ) );
    m_fields["version"]->setText( "1" );
    if( m_minVersion )
      m_minVersion->setCurrentIndex( m_minVersion->count() - 1 );
    renderCandidates();
    return;
  }

  // A new version of a package made before: its details, items, and per-item notes.
  const map<string,vector<UserLibrary::CreatedVersion>> created = UserLibrary::createdPackages();
  const map<string,vector<UserLibrary::CreatedVersion>>::const_iterator pos = created.find( id );
  if( (pos == end(created)) || pos->second.empty() )
    return;

  const UserLibrary::CreatedVersion &last = pos->second.back();
  const UserLibrary::Manifest &m = last.manifest;
  m_fields["name"]->setText( WString::fromUTF8( m.name ) );
  m_fields["version"]->setText( WString::fromUTF8( UserLibraryPackage::nextVersion( m.version ) ) );
  m_fields["author"]->setText( WString::fromUTF8( m.author ) );
  m_fields["organization"]->setText( WString::fromUTF8( m.organization ) );
  m_fields["contact"]->setText( WString::fromUTF8( m.contact ) );
  m_fields["description"]->setText( WString::fromUTF8( m.description ) );
  m_fields["handling"]->setText( WString::fromUTF8( m.handling ) );
  m_fields["review_by"]->setText( WString::fromUTF8( m.reviewBy ) );
  if( m_minVersion )
  {
    m_minVersion->setCurrentIndex( m_minVersion->count() - 1 );
    for( int i = 0; i < m_minVersion->count(); ++i )
    {
      if( m_minVersion->itemText( i ).toUTF8() == m.minInterSpecVersion )
        m_minVersion->setCurrentIndex( i );
    }
  }//if( m_minVersion )
  m_areas["notes"]->setText( WString::fromUTF8( m.notes ) );

  set<string> candidate_keys;
  if( m_candidates )
  {
    for( const Candidate &c : *m_candidates )
      candidate_keys.insert( c.sourceKey );
  }

  for( const UserLibrary::ManifestItem &mi : m.items )
  {
    const map<string,string>::const_iterator src = last.sources.find( mi.itemId );
    if( src == end(last.sources) )
      continue;
    const string &key = src->second;
    m_itemNotes[key] = mi.notes;
    m_itemOverride[key] = mi.intendedOverride;
    if( candidate_keys.count( key ) )
      m_selected[key] = true;
    else
      m_movedSinceLast.push_back( mi.path );
  }

  // Which of them changed since that version (hashing can be slow, so off the session thread).
  map<string,string> previous_sha;  // source key -> sha256
  for( const UserLibrary::ManifestItem &mi : m.items )
  {
    const map<string,string>::const_iterator src = last.sources.find( mi.itemId );
    if( src != end(last.sources) )
      previous_sha[src->second] = mi.sha256;
  }

  vector<UserLibrary::FoundItem> to_hash;
  vector<string> keys;
  if( m_candidates )
  {
    for( const Candidate &c : *m_candidates )
    {
      if( !c.fromDatabase && previous_sha.count( c.sourceKey ) )
      {
        to_hash.push_back( c.item.found );
        keys.push_back( c.sourceKey );
      }
    }
  }

  const string previous_version = m.version;
  const shared_ptr<map<string,string>> changed = make_shared<map<string,string>>();
  const WidgetUtils::WidgetHandle self( this );
  const uint64_t request = m_candidatesRequest;
  AppUtils::run_on_ioservice( [to_hash, keys, previous_sha, previous_version, changed](){
      for( size_t i = 0; i < to_hash.size(); ++i )
      {
        string sha;
        try
        {
          sha = UserLibrary::sha256OfItem( to_hash[i] );
        }catch( std::exception & )
        {
        }
        if( sha != previous_sha.at( keys[i] ) )
          (*changed)[keys[i]] = previous_version;
      }
    },
    [self, changed, request](){
      UserLibraryWindow * const window = self.resolve_as<UserLibraryWindow>();
      if( window && (window->m_candidatesRequest == request) )
      {
        window->m_changedSinceLast = *changed;
        window->renderCandidates();
      }
      if( wApp )
        wApp->triggerUpdate();
    } );

  renderCandidates();
}//startFromChanged()


void UserLibraryWindow::refreshCandidates()
{
  const uint64_t request = ++m_candidatesRequest;

  // Gathered on the session thread: the Detector Efficiency Library's locations (including folders
  //  the user added), the database detector efficiencies are saved in, and added reference-spectra
  //  folders.
  vector<DrfLibrary::Location> drf_locations;
  try
  {
    drf_locations = DrfSelect::libraryLocations( m_interspec );
  }catch( std::exception &e )
  {
    cerr << "UserLibraryWindow: could not get the detector locations: " << e.what() << endl;
  }
  const shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
  const long long user_id = m_interspec->user() ? m_interspec->user().id() : -1;
  const vector<string> ref_dirs = added_reference_spectra_dirs( m_interspec );

  // Colour themes saved in the database, as the JSON the Color Themes window exports.
  vector<Candidate> db_themes;
  try
  {
    if( sql && m_interspec->user() )
    {
      DataBaseUtils::DbTransaction transaction( *sql );
      for( const Wt::Dbo::ptr<ColorThemeInfo> &info : m_interspec->user()->colorThemes() )
      {
        ColorTheme theme;
        theme.setFromDataBase( *info );
        Candidate c;
        c.item.typeId = "color_theme";
        c.item.source.kind = UserLibrary::SourceKind::User;
        c.item.names.push_back( theme.theme_name.toUTF8() );
        c.item.relPath = theme.theme_name.toUTF8();
        c.fromDatabase = true;
        c.generatedData = make_shared<const string>( ColorTheme::toJson( theme ) );
        c.generatedName = UserLibrary::sanitizeFileName( theme.theme_name.toUTF8(), "theme" ) + ".json";
        c.sourceKey = "database-theme:" + std::to_string( info.id() );
        c.packagePath = "color_themes/" + c.generatedName;
        db_themes.push_back( c );
      }
      transaction.commit();
    }//if( sql && user )
  }catch( std::exception &e )
  {
    cerr << "UserLibraryWindow: could not read colour themes: " << e.what() << endl;
  }

  const shared_ptr<vector<Candidate>> result = make_shared<vector<Candidate>>();
  const WidgetUtils::WidgetHandle self( this );
  AppUtils::run_on_ioservice( [result, drf_locations, sql, user_id, ref_dirs, db_themes](){
      set<string> seen_paths;

      // The user's own files in the user data folder, and items of installed packages.
      for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
      {
        for( const UserLibrary::Item &item : UserLibrary::items( type.id, false ) )
        {
          if( item.source.kind == UserLibrary::SourceKind::BuiltIn )
            continue;

          Candidate c;
          c.item = item;
          c.sourceKey = UserLibraryPackage::sourceKey( item );
          try
          {
            c.packagePath = UserLibraryPackage::packagePathFor( item );
          }catch( std::exception & )
          {
            c.packagePath = item.relPath;
          }
          for( const string &s : item.found.support )
            c.includes.push_back( SpecUtils::filename( s ) );
          c.status = UserLibrary::validateItem( type.id, item.found ).status;
          seen_paths.insert( SpecUtils::lexically_normalize_path( item.found.path ) );
          result->push_back( c );
        }//for( loop over items )
      }//for( loop over types )

      // Detector efficiencies in folders added in the Detector Efficiency dialog, and saved in the database.
      map<string,const DrfLibrary::Location *> added_locations;
      for( const DrfLibrary::Location &loc : drf_locations )
      {
        const bool added = (loc.role == DrfLibrary::LocationRole::User)
                           || (loc.role == DrfLibrary::LocationRole::GadrasInstall)
                           || (loc.role == DrfLibrary::LocationRole::IsocsInstall)
                           || (loc.role == DrfLibrary::LocationRole::AngleInstall);
        if( added && !loc.hidden )
          added_locations[loc.id] = &loc;
      }

      DrfLibrary::ScanRequest scan_request;
      scan_request.locations = drf_locations;
      scan_request.sql = sql;
      scan_request.dbUserId = user_id;
      const DrfLibrary::Snapshot snapshot = DrfLibrary::scan( scan_request );

      map<string,Candidate> external;  //by normalized path
      for( const DrfLibrary::Entry &e : snapshot.entries )
      {
        if( !e.summary || e.hidden )
          continue;

        if( e.kind == DrfLibrary::EntryKind::DbRow )
        {
          // Copies of InterSpec's own detectors (recorded when used) are not the user's to share.
          if( e.origin == DrfLibrary::Origin::BuiltIn )
            continue;

          Candidate c;
          c.item.typeId = "drf";
          c.item.source.kind = UserLibrary::SourceKind::User;
          c.item.names.push_back( e.summary->name );
          c.item.relPath = e.summary->name;
          c.fromDatabase = true;
          c.dbEntry = e;
          c.sourceKey = "database:" + std::to_string( e.summary->hash );
          c.packagePath = "drfs/" + UserLibrary::sanitizeFileName( e.summary->name, "detector" ) + ".drf.xml";
          result->push_back( c );
          continue;
        }//if( saved in the database )

        const map<string,const DrfLibrary::Location *>::const_iterator loc = added_locations.find( e.locationId );
        if( (loc == end(added_locations)) || e.path.empty() )
          continue;

        const string norm_path = SpecUtils::lexically_normalize_path( e.path );
        if( seen_paths.count( norm_path ) )
          continue;  //e.g. a location inside the user data folder, already listed

        Candidate &c = external[norm_path];
        if( c.item.typeId.empty() )
        {
          c.item.typeId = "drf";
          c.item.source.kind = UserLibrary::SourceKind::User;
          c.item.found.path = e.path;
          c.item.found.isFolder = (e.kind == DrfLibrary::EntryKind::GadrasFolder);
          c.locationPath = loc->second->path;
          // Relative to the added folder's parent, so the folder's own name is kept in the package.
          c.root = SpecUtils::parent_path( loc->second->path );
          c.item.relPath = UserLibrary::relativePath( c.root, e.path );
          c.sourceKey = "location:" + e.path;
        }
        if( !e.companionPath.empty()
            && (std::find( begin(c.item.found.support), end(c.item.found.support), e.companionPath ) == end(c.item.found.support)) )
          c.item.found.support.push_back( e.companionPath );
        if( std::find( begin(c.item.names), end(c.item.names), e.summary->name ) == end(c.item.names) )
          c.item.names.push_back( e.summary->name );
      }//for( loop over entries )

      for( pair<const string,Candidate> &p : external )
      {
        Candidate &c = p.second;
        UserLibraryPackage::CreateItem ci;
        ci.item = c.item;
        ci.root = c.root;
        try
        {
          c.packagePath = UserLibraryPackage::packagePathFor( ci );
        }catch( std::exception & )
        {
          c.packagePath = c.item.relPath;
        }
        for( const string &s : c.item.found.support )
          c.includes.push_back( SpecUtils::filename( s ) );
        c.status = UserLibrary::validateItem( "drf", c.item.found ).status;
        result->push_back( c );
      }//for( loop over files in added folders )

      result->insert( end(*result), begin(db_themes), end(db_themes) );

      // Reference-spectra folders added in the Reference Spectra tool: each a set.
      for( const string &dir : ref_dirs )
      {
        if( seen_paths.count( SpecUtils::lexically_normalize_path( dir ) ) )
          continue;

        Candidate c;
        c.item.typeId = "ref_spectra";
        c.item.source.kind = UserLibrary::SourceKind::User;
        c.item.found.path = dir;
        c.item.found.isFolder = true;
        c.item.names.push_back( SpecUtils::filename( dir ) );
        c.item.relPath = SpecUtils::filename( dir );
        c.locationPath = dir;
        c.root = SpecUtils::parent_path( dir );
        c.sourceKey = "location:" + dir;
        c.packagePath = "reference_spectra/" + SpecUtils::filename( dir );
        c.status = UserLibrary::validateItem( "ref_spectra", c.item.found ).status;
        result->push_back( c );
      }//for( loop over added reference-spectra folders )
    },
    [self, result, request](){
      UserLibraryWindow * const window = self.resolve_as<UserLibraryWindow>();
      if( window && (window->m_candidatesRequest == request) )
        window->candidatesListed( request, result );
      if( wApp )
        wApp->triggerUpdate();
    } );
}//refreshCandidates()


void UserLibraryWindow::candidatesListed( const uint64_t, std::shared_ptr<const std::vector<Candidate>> candidates )
{
  const bool first_time = !m_candidates;
  m_candidates = candidates;

  // A "new version" selection made before the list arrived: match its items now.
  if( first_time && m_startFrom && (m_startFrom->currentIndex() > 0) )
    startFromChanged();
  else
    renderCandidates();
}//candidatesListed(...)


void UserLibraryWindow::renderCandidates()
{
  if( !m_candidatesList )
    return;

  m_candidatesList->clear();

  for( const string &path : m_movedSinceLast )
    add_text( m_candidatesList, WString::tr( "ul-create-moved" ).arg( WString::fromUTF8( path ) ), "UserLibWarn" );

  if( !m_candidates )
  {
    add_text( m_candidatesList, WString::tr( "ul-loading" ), "UserLibBusy" );
    return;
  }

  if( m_candidates->empty() )
  {
    add_text( m_candidatesList, WString::tr( "ul-create-no-items" ), "UserLibNone" );
    m_numSelected->setText( "" );
    return;
  }

  const string query = SpecUtils::trim_copy( m_candidateSearch->text().toUTF8() );
  size_t num_selected = 0;
  for( const Candidate &c : *m_candidates )
  {
    const map<string,bool>::const_iterator sel = m_selected.find( c.sourceKey );
    num_selected += ((sel != end(m_selected)) && sel->second);
  }
  m_numSelected->setText( WString::tr( "ul-create-num-selected" ).arg( static_cast<int>( num_selected ) ) );

  for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
  {
    WContainerWidget *group = nullptr;
    for( const Candidate &c : *m_candidates )
    {
      if( c.item.typeId != type.id )
        continue;

      const map<string,bool>::const_iterator sel = m_selected.find( c.sourceKey );
      const bool selected = (sel != end(m_selected)) && sel->second;

      if( !query.empty() && !selected )
      {
        bool match = SpecUtils::icontains( c.packagePath, query ) || SpecUtils::icontains( c.item.source.label, query );
        for( const string &n : c.item.names )
          match = match || SpecUtils::icontains( n, query );
        if( !match )
          continue;
      }

      if( !group )
      {
        group = m_candidatesList->addNew<WContainerWidget>();
        group->addStyleClass( "UserLibGroup" );
        add_text( group, WString::tr( type.titleKey ), "UserLibGroupTitle" );
      }

      WContainerWidget *div = group->addNew<WContainerWidget>();
      div->addStyleClass( "UserLibItem UserLibCandidate" );

      WCheckBox *cb = div->addNew<WCheckBox>( names_label( c.item ) );
      cb->addStyleClass( "UserLibItemName" );
      cb->setChecked( selected );
      const string key = c.sourceKey;
      cb->changed().connect( [this, key, cb](){
        m_selected[key] = cb->isChecked();
        renderCandidates();
      } );

      WString source_label = sourceLabel( c.item.source );
      if( c.fromDatabase )
        source_label = WString::tr( "ul-source-database" );
      else if( !c.locationPath.empty() )
        source_label = WString::tr( "ul-source-location" ).arg( WString::fromUTF8( c.locationPath ) );
      add_text( div, source_label, "UserLibPkgMeta" );
      add_text( div, WString::fromUTF8( c.packagePath ), "UserLibPath" );

      if( !c.includes.empty() )
      {
        string includes;
        for( const string &i : c.includes )
          includes += (includes.empty() ? "" : ", ") + i;
        add_text( div, WString::tr( "ul-create-includes" ).arg( WString::fromUTF8( includes ) ), "UserLibItemNotes" );
      }
      if( c.item.source.kind == UserLibrary::SourceKind::Package )
        add_text( div, WString::tr( "ul-create-from-package" ).arg( WString::fromUTF8( c.item.source.label ) ), "UserLibItemNotes" );
      if( c.status == UserLibrary::Status::Unreadable )
        add_text( div, WString::tr( "ul-item-unreadable" ), "UserLibError" );
      else if( c.status == UserLibrary::Status::NeedsNewer )
        add_text( div, WString::tr( "ul-item-needs-newer" ), "UserLibWarn" );
      if( m_changedSinceLast.count( key ) )
        add_text( div, WString::tr( "ul-create-changed" ), "UserLibWarn" );

      if( selected )
      {
        WContainerWidget *opts = div->addNew<WContainerWidget>();
        opts->addStyleClass( "UserLibItemOptions" );
        WLineEdit *notes = opts->addNew<WLineEdit>( WString::fromUTF8( m_itemNotes[key] ) );
        notes->setPlaceholderText( WString::tr( "ul-create-item-notes" ) );
        notes->changed().connect( [this, key, notes](){ m_itemNotes[key] = notes->text().toUTF8(); } );
        notes->textInput().connect( [this, key, notes](){ m_itemNotes[key] = notes->text().toUTF8(); } );
        WCheckBox *override_cb = opts->addNew<WCheckBox>( WString::tr( "ul-create-item-override" ) );
        override_cb->setInline( false );
        override_cb->setChecked( m_itemOverride[key] );
        override_cb->changed().connect( [this, key, override_cb](){ m_itemOverride[key] = override_cb->isChecked(); } );
      }//if( selected )
    }//for( loop over candidates )
  }//for( loop over types )
}//renderCandidates()


void UserLibraryWindow::startCreate()
{
  if( m_creating || !m_candidates )
    return;

  UserLibraryPackage::CreateRequest request;
  const int index = m_startFrom->currentIndex();
  if( (index > 0) && (index < static_cast<int>( m_startFromIds.size() )) )
    request.manifest.id = m_startFromIds[index];

  auto field = [this]( const char *key ) -> string {
    return SpecUtils::trim_copy( m_fields[key]->text().toUTF8() );
  };
  request.manifest.name = field( "name" );
  request.manifest.version = field( "version" );
  request.manifest.author = field( "author" );
  request.manifest.organization = field( "organization" );
  request.manifest.contact = field( "contact" );
  request.manifest.description = field( "description" );
  request.manifest.handling = field( "handling" );
  request.manifest.reviewBy = field( "review_by" );
  request.manifest.minInterSpecVersion = m_minVersion ? m_minVersion->currentText().toUTF8() : string();
  request.manifest.notes = m_areas["notes"]->text().toUTF8();
  request.manifest.changes = m_areas["changes"]->text().toUTF8();

  vector<pair<size_t,DrfLibrary::Entry>> to_load;  //database detector efficiencies, by request item
  for( const Candidate &c : *m_candidates )
  {
    const map<string,bool>::const_iterator sel = m_selected.find( c.sourceKey );
    if( (sel == end(m_selected)) || !sel->second )
      continue;
    UserLibraryPackage::CreateItem ci;
    ci.item = c.item;
    ci.root = c.root;
    ci.sourceKey = c.sourceKey;
    ci.notes = m_itemNotes[c.sourceKey];
    ci.intendedOverride = m_itemOverride[c.sourceKey];
    ci.generatedData = c.generatedData;
    ci.generatedName = c.generatedName;
    if( c.fromDatabase && !c.generatedData )
      to_load.emplace_back( request.items.size(), c.dbEntry );
    request.items.push_back( ci );
  }
  const shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();

  m_creating = true;
  m_createBtn->disable();
  m_createStatus->clear();
  add_text( m_createStatus, WString::tr( "ul-create-checking" ), "UserLibBusy" );

  struct Result
  {
    string error;
    shared_ptr<const string> data;
    UserLibrary::Manifest manifest;
    vector<string> issues, warnings;
  };
  const shared_ptr<Result> result = make_shared<Result>();
  const WidgetUtils::WidgetHandle self( this );
  AppUtils::run_on_ioservice( [request, result, to_load, sql]() mutable {
      try
      {
        for( const pair<size_t,DrfLibrary::Entry> &load : to_load )
          request.items[load.first].databaseDrf = DrfLibrary::loadDrf( load.second, sql );

        for( const UserLibraryPackage::PreflightIssue &issue : UserLibraryPackage::preflight( request ) )
          (issue.blocking ? result->issues : result->warnings).push_back( issue.message );
        if( !result->issues.empty() )
          return;

        std::ostringstream output;
        result->manifest = UserLibraryPackage::create( request, output );
        result->data = make_shared<const string>( output.str() );
      }catch( std::exception &e )
      {
        result->error = e.what();
      }
    },
    [self, result](){
      UserLibraryWindow * const window = self.resolve_as<UserLibraryWindow>();
      if( window )
        window->createDone( result->error, result->data, result->manifest, result->issues, result->warnings );
      if( wApp )
        wApp->triggerUpdate();
    } );
}//startCreate()


void UserLibraryWindow::createDone( const std::string &error, std::shared_ptr<const std::string> data,
                                    const UserLibrary::Manifest &manifest, const std::vector<std::string> &issues,
                                    const std::vector<std::string> &warnings )
{
  m_creating = false;
  m_createBtn->enable();
  m_createStatus->clear();

  if( !error.empty() )
  {
    add_text( m_createStatus, WString::tr( "ul-create-error" ).arg( WString::fromUTF8( error ) ), "UserLibError" );
    return;
  }

  if( !issues.empty() )
  {
    add_text( m_createStatus, WString::tr( "ul-create-issues" ), "UserLibError" );
    for( const string &msg : issues )
      add_text( m_createStatus, WString::fromUTF8( "• " + msg ), "UserLibError" );
  }
  if( !warnings.empty() )
  {
    add_text( m_createStatus, WString::tr( "ul-create-warnings" ), "UserLibWarn" );
    for( const string &msg : warnings )
      add_text( m_createStatus, WString::fromUTF8( "• " + msg ), "UserLibWarn" );
  }
  if( !issues.empty() || !data )
    return;

  add_text( m_createStatus, WString::tr( "ul-create-done" ).arg( WString::fromUTF8( manifest.name ) )
                              .arg( WString::fromUTF8( manifest.version ) ), "UserLibOk" );

  const string filename = UserLibrary::sanitizeFileName( manifest.name + " v" + manifest.version, "package" )
                          + UserLibraryPackage::FileExtension;
  WLink link( std::make_shared<PackageResource>( data, filename ) );
  link.setTarget( LinkTarget::NewWindow );
  WAnchor *download = m_createStatus->addNew<WAnchor>( link, WString::tr( "ul-create-download" ) );
  download->addStyleClass( "LinkBtn UserLibDownload" );

  // The package is now one this user made: the form becomes its next version (a later version
  //  number, the same id) - never another package with this one's id.
  const string id = manifest.id;
  updateStartFromOptions();
  for( size_t i = 0; i < m_startFromIds.size(); ++i )
  {
    if( m_startFromIds[i] == id )
      m_startFrom->setCurrentIndex( static_cast<int>( i ) );
  }
  startFromChanged();
}//createDone(...)

#endif //#if( USE_USER_LIBRARY )
