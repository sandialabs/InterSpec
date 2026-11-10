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
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

#include <Wt/WText.h>
#include <Wt/WCheckBox.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WButtonGroup.h>
#include <Wt/WRadioButton.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/UserLibrary.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/UserLibraryWindow.h"
#include "InterSpec/UserLibraryPackage.h"
#include "InterSpec/UserLibraryImportDialog.h"

using namespace std;
using namespace Wt;


namespace
{
  std::filesystem::path to_fs_path( const string &utf8 )
  {
#ifdef _WIN32
    return std::filesystem::path( SpecUtils::convert_from_utf8_to_utf16( utf8 ) );
#else
    return std::filesystem::path( utf8 );
#endif
  }


  WText *add_text( WContainerWidget *parent, const WString &txt, const char *style_class )
  {
    WText *t = parent->addNew<WText>( txt, TextFormat::Plain );
    t->setInline( false );
    if( style_class )
      t->addStyleClass( style_class );
    return t;
  }


  /** "A, B, C and 4 more", or the file name when there are no names. */
  WString names_label( const vector<string> &names, const string &rel_path )
  {
    if( names.empty() )
      return WString::fromUTF8( SpecUtils::filename( rel_path ) );

    string joined;
    for( size_t i = 0; (i < names.size()) && (i < 3); ++i )
      joined += (i ? ", " : "") + names[i];
    WString answer = WString::fromUTF8( joined );
    if( names.size() > 3 )
    {
      answer += WString::fromUTF8( " " );
      answer += WString::tr( "ul-item-more-names" ).arg( static_cast<int>( names.size() - 3 ) );
    }
    return answer;
  }//names_label(...)


  const char *state_key( const UserLibraryPackage::ItemState state )
  {
    using UserLibraryPackage::ItemState;
    switch( state )
    {
      case ItemState::New:         return "ul-import-state-new";
      case ItemState::AlreadyHave: return "ul-import-state-have";
      case ItemState::Unchanged:   return "ul-import-state-unchanged";
      case ItemState::Update:      return "ul-import-state-update";
      case ItemState::Conflict:    return "ul-import-state-conflict";
      case ItemState::NeedsNewer:  return "ul-import-state-newer";
      case ItemState::Unreadable:  return "ul-import-state-unreadable";
    }
    return "ul-import-state-new";
  }//state_key(...)


  /** The text of a sensitive item: an LLM skill's markdown files (SKILL.md first), or the file itself. */
  string sensitive_text( const UserLibrary::FoundItem &item )
  {
    string answer;
    vector<string> files = item.isFolder ? SpecUtils::recursive_ls( item.path, ".md" )
                                         : vector<string>{ item.path };
    std::sort( begin(files), end(files), []( const string &lhs, const string &rhs ){
      const bool lhs_main = (SpecUtils::filename( lhs ) == "SKILL.md");
      const bool rhs_main = (SpecUtils::filename( rhs ) == "SKILL.md");
      return (lhs_main != rhs_main) ? lhs_main : (lhs < rhs);
    } );

    for( const string &f : files )
    {
      std::ifstream input( to_fs_path( f ), ios::in | ios::binary );
      std::ostringstream strm;
      strm << input.rdbuf();
      if( item.isFolder )
        answer += "=== " + UserLibrary::relativePath( item.path, f ) + " ===\n";
      answer += strm.str() + "\n\n";
      if( answer.size() > 200*1024 )
      {
        answer.resize( 200*1024 );
        answer += "\n…";
        break;
      }
    }//for( loop over markdown files )

    return answer;
  }//sensitive_text(...)
}//namespace


UserLibraryImportDialog *UserLibraryImportDialog::start( const std::string &archivePath,
                                                         const std::string &displayName )
{
  // The upload's spool file goes away when the request finishes, so work on a copy.
  string copy = SpecUtils::temp_file_name( "interspec_pkg_" + UserLibrary::newUuid().substr( 0, 8 ),
                                           SpecUtils::temp_dir() );
  std::error_code ec;
  std::filesystem::copy_file( to_fs_path( archivePath ), to_fs_path( copy ),
                              std::filesystem::copy_options::overwrite_existing, ec );
  if( ec )
  {
    cerr << "UserLibraryImportDialog: could not copy '" << archivePath << "': " << ec.message() << endl;
    copy.clear();
  }

  UserLibraryImportDialog *dialog = AuxWindow::make<UserLibraryImportDialog>( displayName, copy );
  dialog->finished().connect( [dialog](){ AuxWindow::deleteAuxWindow( dialog ); } );

  if( copy.empty() )
  {
    dialog->showError( WString::tr( "ul-import-error" ).arg( WString::fromUTF8( displayName ) )
                         .arg( WString::fromUTF8( ec.message() ) ) );
    return dialog;
  }

  const WidgetUtils::WidgetHandle self( dialog );
  const shared_ptr<shared_ptr<UserLibraryPackage::Inspection>> result
                                    = make_shared<shared_ptr<UserLibraryPackage::Inspection>>();
  AppUtils::run_on_ioservice( [copy, displayName, result](){
      *result = UserLibraryPackage::inspect( copy, displayName );
    },
    [self, result](){
      UserLibraryImportDialog * const dlg = self.resolve_as<UserLibraryImportDialog>();
      if( dlg )
        dlg->inspectionDone( *result );
      else if( *result )
        UserLibraryPackage::discard( **result );  // the dialog was closed while checking
      if( wApp )
        wApp->triggerUpdate();
    } );

  return dialog;
}//start(...)


UserLibraryImportDialog::UserLibraryImportDialog( const std::string &displayName, const std::string &archiveCopy )
  : AuxWindow( WString::tr("ul-import-title"),
               (AuxWindowProperties::TabletNotFullScreen | AuxWindowProperties::DisableCollapse
                | AuxWindowProperties::SetCloseable) ),
    m_displayName( displayName ),
    m_archiveCopy( archiveCopy ),
    m_inspection(),
    m_installing( false ),
    m_installed( false ),
    m_body( nullptr ),
    m_install( nullptr ),
    m_cancel( nullptr )
{
  InterSpec * const viewer = InterSpec::instance();
  if( viewer )
    viewer->useMessageResourceBundle( "UserLibrary" );
  wApp->useStyleSheet( "InterSpec_resources/UserLibrary.css" );

  addStyleClass( "UserLibImportDialog" );
  contents()->setOverflow( Overflow::Auto, Orientation::Vertical );
  m_body = contents()->addNew<WContainerWidget>();
  m_body->addStyleClass( "UserLibImport" );

  add_text( m_body, WString::tr( "ul-import-checking" ).arg( WString::fromUTF8( displayName ) ), "UserLibBusy" );

  m_install = footer()->addNew<WPushButton>( WString::tr( "ul-install" ) );
  WidgetUtils::applyButtonRole( m_install, WidgetUtils::ButtonRole::Affirm );
  m_install->disable();
  m_install->clicked().connect( this, &UserLibraryImportDialog::startInstall );

  m_cancel = addCloseButtonToFooter( WString::tr( "Cancel" ), WidgetUtils::ButtonRole::Dismiss );
  m_cancel->clicked().connect( this, &AuxWindow::hide );

  rejectWhenEscapePressed();

  if( viewer && !viewer->isPhone() )
  {
    const int w = viewer->renderedWidth(), h = viewer->renderedHeight();
    if( (w > 0) && (h > 0) )
      resizeWindow( std::min( 720, (9*w)/10 ), std::min( 760, (9*h)/10 ) );
  }
  show();
  centerWindow();
}//UserLibraryImportDialog constructor


UserLibraryImportDialog::~UserLibraryImportDialog()
{
  if( m_inspection && !m_installed && !m_installing )
    UserLibraryPackage::discard( *m_inspection );

  if( !m_archiveCopy.empty() )
    SpecUtils::remove_file( m_archiveCopy );
}//~UserLibraryImportDialog()


void UserLibraryImportDialog::showError( const Wt::WString &message )
{
  m_body->clear();
  add_text( m_body, message, "UserLibError" );
  m_install->hide();
  m_cancel->setText( WString::tr( "Close" ) );
}//showError(...)


void UserLibraryImportDialog::inspectionDone( std::shared_ptr<UserLibraryPackage::Inspection> inspection )
{
  m_inspection = inspection;
  if( !inspection )
  {
    showError( WString::tr( "ul-import-error" ).arg( WString::fromUTF8( m_displayName ) )
                 .arg( WString::tr( "ul-import-unknown-error" ) ) );
    return;
  }

  if( !inspection->error.empty() )
  {
    showError( WString::tr( "ul-import-error" ).arg( WString::fromUTF8( m_displayName ) )
                 .arg( WString::fromUTF8( inspection->error ) ) );
    return;
  }

  render();
}//inspectionDone(...)


void UserLibraryImportDialog::render()
{
  using UserLibraryPackage::ItemState;
  using UserLibraryPackage::Relation;

  m_body->clear();
  m_useIncoming.clear();
  m_keepExisting.clear();
  m_sensitiveConfirms.clear();

  const UserLibraryPackage::Inspection &insp = *m_inspection;
  const UserLibrary::Manifest &m = insp.manifest;

  // --- The package ---
  WContainerWidget *head = m_body->addNew<WContainerWidget>();
  head->addStyleClass( "UserLibPkgHead" );
  WContainerWidget *title_line = head->addNew<WContainerWidget>();
  title_line->addStyleClass( "UserLibPkgTitleLine" );
  title_line->addNew<WText>( WString::fromUTF8( m.name ), TextFormat::Plain )->addStyleClass( "UserLibPkgName" );
  title_line->addNew<WText>( WString::tr( "ul-pkg-version" ).arg( WString::fromUTF8( m.version ) ), TextFormat::Plain )
    ->addStyleClass( "UserLibPkgVersion" );

  if( !m.handling.empty() )
    add_text( head, WString::tr( "ul-pkg-handling" ).arg( WString::fromUTF8( m.handling ) ), "UserLibHandling" );

  string by;
  for( const string &s : { m.author, m.organization } )
    if( !s.empty() )
      by += (by.empty() ? "" : " · ") + s;
  if( !by.empty() )
    add_text( head, WString::tr( "ul-pkg-by" ).arg( WString::fromUTF8( by ) ), "UserLibPkgMeta" );
  if( !m.contact.empty() )
    add_text( head, WString::tr( "ul-pkg-contact" ).arg( WString::fromUTF8( m.contact ) ), "UserLibPkgMeta" );
  if( !m.created.empty() )
    add_text( head, WString::tr( "ul-pkg-created" ).arg( WString::fromUTF8( m.created.substr( 0, 10 ) ) ), "UserLibPkgMeta" );
  if( !m.description.empty() )
    add_text( head, WString::fromUTF8( m.description ), "UserLibPkgDesc" );
  if( !m.notes.empty() )
    add_text( head, WString::fromUTF8( m.notes ), "UserLibPkgNotes" );
  if( !insp.hadManifest )
    add_text( head, WString::tr( "ul-import-no-manifest" ), "UserLibWarn" );

  // --- Compared with what is installed ---
  if( insp.installed )
  {
    const string old_version = insp.installed->manifest.version;
    switch( insp.relation )
    {
      case Relation::NewPackage:
        break;

      case Relation::Update:
        add_text( m_body, WString::tr( "ul-import-update" ).arg( WString::fromUTF8( old_version ) )
                            .arg( WString::fromUTF8( m.version ) ), "UserLibInfo" );
        break;

      case Relation::Reinstall:
        add_text( m_body, WString::tr( "ul-import-reinstall" ).arg( WString::fromUTF8( m.version ) ), "UserLibInfo" );
        break;

      case Relation::SameVersionDifferent:
        add_text( m_body, WString::tr( "ul-import-same-version" ).arg( WString::fromUTF8( m.version ) ), "UserLibWarn" );
        break;

      case Relation::Downgrade:
        add_text( m_body, WString::tr( "ul-import-downgrade" ).arg( WString::fromUTF8( m.version ) )
                            .arg( WString::fromUTF8( old_version ) ), "UserLibWarn" );
        break;
    }//switch( insp.relation )

    if( insp.relation != Relation::Reinstall )
      add_text( m_body, WString::tr( "ul-import-update-summary" ).arg( static_cast<int>( insp.added.size() ) )
                          .arg( static_cast<int>( insp.changed.size() ) ).arg( static_cast<int>( insp.removed.size() ) ),
                "UserLibPkgMeta" );

    if( !insp.locallyModified.empty() )
      add_text( m_body, WString::tr( "ul-import-modified" ).arg( static_cast<int>( insp.locallyModified.size() ) ), "UserLibWarn" );

    if( insp.authorChanged )
    {
      const UserLibrary::Manifest &old_m = insp.installed->manifest;
      const string was = old_m.author + (old_m.organization.empty() ? "" : (" · " + old_m.organization));
      const string now = m.author + (m.organization.empty() ? "" : (" · " + m.organization));
      add_text( m_body, WString::tr( "ul-import-author-changed" ).arg( WString::fromUTF8( was ) )
                          .arg( WString::fromUTF8( now ) ), "UserLibWarn" );
    }

    if( !insp.changesSince.empty() )
    {
      WContainerWidget *changes = m_body->addNew<WContainerWidget>();
      changes->addStyleClass( "UserLibChanges" );
      add_text( changes, WString::tr( "ul-import-changes-since" ), "UserLibSubhead" );
      for( const UserLibrary::ManifestPreviousVersion &v : insp.changesSince )
      {
        const string when = v.created.empty() ? string() : (" (" + v.created.substr( 0, 10 ) + ")");
        add_text( changes, WString::tr( "ul-pkg-history-entry" ).arg( WString::fromUTF8( v.version ) )
                             .arg( WString::fromUTF8( when ) ), "UserLibHistVersion" );
        if( !v.changes.empty() )
          add_text( changes, WString::fromUTF8( v.changes ), "UserLibHistText" );
      }
    }//if( !insp.changesSince.empty() )
  }//if( insp.installed )

  // --- Summary, and the conflict shortcuts ---
  size_t num_new = 0, num_have = 0, num_choice = 0;
  for( const UserLibraryPackage::IncomingItem &item : insp.items )
  {
    num_new += (item.state == ItemState::New) || (item.state == ItemState::Update);
    num_have += (item.state == ItemState::AlreadyHave) || (item.state == ItemState::Unchanged);
    num_choice += (item.state == ItemState::Conflict);
  }

  WContainerWidget *summary = m_body->addNew<WContainerWidget>();
  summary->addStyleClass( "UserLibSummary" );
  summary->addNew<WText>( WString::tr( "ul-import-summary" ).arg( static_cast<int>( num_new ) )
                            .arg( static_cast<int>( num_have ) ).arg( static_cast<int>( num_choice ) ),
                          TextFormat::Plain );
  if( num_choice )
  {
    WPushButton *keep_all = summary->addNew<WPushButton>( WString::tr( "ul-import-keep-all" ) );
    keep_all->addStyleClass( "LinkBtn" );
    keep_all->clicked().connect( [this](){ setAllConflicts( false ); } );
    WPushButton *use_all = summary->addNew<WPushButton>( WString::tr( "ul-import-use-all" ) );
    use_all->addStyleClass( "LinkBtn" );
    use_all->clicked().connect( [this](){ setAllConflicts( true ); } );
  }

  // --- The items, by type ---
  for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
  {
    vector<const UserLibraryPackage::IncomingItem *> of_type;
    for( const UserLibraryPackage::IncomingItem &item : insp.items )
      if( item.typeId == type.id )
        of_type.push_back( &item );
    if( of_type.empty() )
      continue;

    WContainerWidget *group = m_body->addNew<WContainerWidget>();
    group->addStyleClass( "UserLibGroup" );
    add_text( group, WString::tr( type.titleKey ), "UserLibGroupTitle" );

    for( const UserLibraryPackage::IncomingItem *item_ptr : of_type )
    {
      const UserLibraryPackage::IncomingItem &item = *item_ptr;
      WContainerWidget *row = group->addNew<WContainerWidget>();
      row->addStyleClass( "UserLibItem" );

      WContainerWidget *line = row->addNew<WContainerWidget>();
      line->addStyleClass( "UserLibItemLine" );
      line->addNew<WText>( names_label( item.validation.names, item.relPath ), TextFormat::Plain )
        ->addStyleClass( "UserLibItemName" );
      WText *badge = line->addNew<WText>( WString::tr( state_key( item.state ) ), TextFormat::Plain );
      badge->addStyleClass( "UserLibBadge" );
      if( (item.state == ItemState::Conflict) || (item.state == ItemState::NeedsNewer) )
        badge->addStyleClass( "Warn" );
      else if( item.state == ItemState::Unreadable )
        badge->addStyleClass( "Error" );

      add_text( row, WString::fromUTF8( item.relPath ), "UserLibPath" );
      if( !item.notes.empty() )
        add_text( row, WString::fromUTF8( item.notes ), "UserLibItemNotes" );
      if( item.hashMismatch )
        add_text( row, WString::tr( "ul-import-hash-mismatch" ), "UserLibWarn" );
      if( ((item.state == ItemState::Unreadable) || (item.state == ItemState::NeedsNewer))
          && !item.validation.message.empty() )
        add_text( row, WString::fromUTF8( item.validation.message ), "UserLibItemNotes" );

      // Name conflicts the user decides.
      set<string> asked;
      for( const UserLibraryPackage::NameConflict &c : item.conflicts )
      {
        if( c.identical || c.alreadyChosen || !asked.insert( c.name ).second )
          continue;

        WContainerWidget *conflict = row->addNew<WContainerWidget>();
        conflict->addStyleClass( "UserLibConflict" );

        WString source;
        switch( c.existing.source.kind )
        {
          case UserLibrary::SourceKind::BuiltIn: source = WString::tr( "ul-source-builtin" ); break;
          case UserLibrary::SourceKind::User:    source = WString::tr( "ul-source-user" );    break;
          case UserLibrary::SourceKind::Package:
            source = WString::tr( "ul-source-package" ).arg( WString::fromUTF8( c.existing.source.label ) );
            break;
        }

        add_text( conflict, c.overridesBuiltIn
                              ? WString::tr( "ul-import-conflict-builtin" ).arg( WString::fromUTF8( c.name ) )
                              : WString::tr( "ul-import-conflict-with" ).arg( WString::fromUTF8( c.name ) ).arg( source ),
                  nullptr );

        // The radio buttons keep their group alive (each holds a shared_ptr to it).
        const shared_ptr<WButtonGroup> button_group = std::make_shared<WButtonGroup>();
        WRadioButton *keep = conflict->addNew<WRadioButton>( WString::tr( "ul-import-keep-existing" ) );
        keep->setInline( false );
        WRadioButton *use = conflict->addNew<WRadioButton>( WString::tr( "ul-import-use-incoming" ) );
        use->setInline( false );
        button_group->addButton( keep );
        button_group->addButton( use );
        (c.useIncomingByDefault ? use : keep)->setChecked( true );

        const pair<string,string> key( item.typeId, c.name );
        m_useIncoming[key] = use;
        m_keepExisting[key] = keep;

        if( type.kind == UserLibrary::Kind::List )
          add_text( conflict, WString::tr( "ul-import-both-listed" ), "UserLibItemNotes" );
      }//for( loop over conflicts )

      // Items that steer the AI assistant: shown, and confirmed one by one.
      if( type.sensitive && (item.state != ItemState::Unreadable) )
      {
        WContainerWidget *sensitive = row->addNew<WContainerWidget>();
        sensitive->addStyleClass( "UserLibSensitive" );
        add_text( sensitive, WString::tr( "ul-import-sensitive" ), "UserLibWarn" );

        WPushButton *show_text = sensitive->addNew<WPushButton>( WString::tr( "ul-import-sensitive-show" ) );
        show_text->addStyleClass( "LinkBtn" );
        WText *text = add_text( sensitive, WString::fromUTF8( sensitive_text( item.found ) ), "UserLibSkillText" );
        text->hide();
        show_text->clicked().connect( [text, show_text](){ text->show(); show_text->hide(); } );

        WCheckBox *confirm = sensitive->addNew<WCheckBox>( WString::tr( "ul-import-sensitive-confirm" ) );
        confirm->setInline( false );
        confirm->changed().connect( this, &UserLibraryImportDialog::updateInstallButton );
        m_sensitiveConfirms.push_back( confirm );
      }//if( sensitive )
    }//for( loop over items of this type )
  }//for( loop over types )

  // --- What is not installed ---
  for( const bool refused : { true, false } )
  {
    vector<string> files;
    for( const UserLibraryPackage::SkippedFile &s : insp.skipped )
      if( s.refused == refused )
        files.push_back( s.relPath );
    if( files.empty() )
      continue;

    WContainerWidget *skipped = m_body->addNew<WContainerWidget>();
    skipped->addStyleClass( "UserLibSkipped" );
    add_text( skipped, WString::tr( "ul-import-not-installed" )
                         .arg( WString::tr( refused ? "ul-import-refused" : "ul-import-ignored" ) ),
              refused ? "UserLibWarn" : "UserLibSubhead" );
    for( size_t i = 0; (i < files.size()) && (i < 50); ++i )
      add_text( skipped, WString::fromUTF8( files[i] ), "UserLibPath" );
    if( files.size() > 50 )
      add_text( skipped, WString::tr( "ul-item-more-names" ).arg( static_cast<int>( files.size() - 50 ) ), "UserLibPath" );
  }//for( refused, then ignored )

  // --- Fingerprint ---
  WContainerWidget *sha = m_body->addNew<WContainerWidget>();
  sha->addStyleClass( "UserLibSha" );
  add_text( sha, WString::tr( "ul-import-sha" ), "UserLibPkgMeta" );
  add_text( sha, WString::fromUTF8( insp.archiveSha256 ), "UserLibPath UserLibMono" );

  m_install->show();
  updateInstallButton();
}//render()


void UserLibraryImportDialog::setAllConflicts( const bool useIncoming )
{
  for( pair<const pair<string,string>,WRadioButton *> &r : (useIncoming ? m_useIncoming : m_keepExisting) )
    r.second->setChecked( true );
  for( pair<const pair<string,string>,WRadioButton *> &r : (useIncoming ? m_keepExisting : m_useIncoming) )
    r.second->setChecked( false );
}//setAllConflicts(...)


void UserLibraryImportDialog::updateInstallButton()
{
  bool all_confirmed = true;
  for( const WCheckBox *cb : m_sensitiveConfirms )
    all_confirmed = all_confirmed && cb->isChecked();

  const bool can_install = m_inspection && m_inspection->error.empty() && !m_installing && !m_installed
                           && !m_inspection->items.empty() && all_confirmed;
  m_install->setEnabled( can_install );
}//updateInstallButton()


void UserLibraryImportDialog::startInstall()
{
  if( !m_inspection || m_installing || m_installed )
    return;

  UserLibraryPackage::Decisions decisions;
  for( const pair<const pair<string,string>,WRadioButton *> &r : m_useIncoming )
    decisions[r.first] = r.second->isChecked();

  m_installing = true;
  m_install->disable();
  m_cancel->disable();
  m_body->clear();
  add_text( m_body, WString::tr( "ul-loading" ), "UserLibBusy" );

  const shared_ptr<UserLibraryPackage::Inspection> inspection = m_inspection;
  const shared_ptr<string> error = make_shared<string>();
  const WidgetUtils::WidgetHandle self( this );
  AppUtils::run_on_ioservice( [inspection, decisions, error](){
      try
      {
        UserLibraryPackage::install( *inspection, decisions );
      }catch( std::exception &e )
      {
        *error = e.what();
        UserLibraryPackage::discard( *inspection );
      }
    },
    [self, error, inspection](){
      UserLibraryImportDialog * const dlg = self.resolve_as<UserLibraryImportDialog>();
      if( dlg )
        dlg->installDone( *error, inspection->manifest.name );
      InterSpec * const viewer = InterSpec::instance();
      if( viewer )
        viewer->userLibraryChanged().emit();
      if( wApp )
        wApp->triggerUpdate();
    } );
}//startInstall()


void UserLibraryImportDialog::installDone( const std::string &error, const std::string &packageName )
{
  m_installing = false;
  m_cancel->enable();

  if( !error.empty() )
  {
    showError( WString::tr( "ul-import-install-error" ).arg( WString::fromUTF8( error ) ) );
    return;
  }

  m_installed = true;
  m_install->hide();
  m_cancel->setText( WString::tr( "Close" ) );

  m_body->clear();
  add_text( m_body, WString::tr( "ul-import-installed" ).arg( WString::fromUTF8( packageName ) ), "UserLibOk" );

  // Where each kind of item it added can be used.
  set<string> types;
  for( const UserLibraryPackage::IncomingItem &item : m_inspection->items )
    types.insert( item.typeId );

  // Some kinds are read when InterSpec starts (reference lines, materials).
  for( const string &type_id : types )
  {
    const UserLibrary::AssetType * const type = UserLibrary::assetType( type_id );
    if( type && (type->apply == UserLibrary::ApplyMode::Restart) )
    {
      add_text( m_body, WString::tr( "ul-import-restart" ), "UserLibWarn" );
      break;
    }
  }

  if( !types.empty() )
  {
    add_text( m_body, WString::tr( "ul-import-where" ), "UserLibSubhead" );
    WContainerWidget *links = m_body->addNew<WContainerWidget>();
    links->addStyleClass( "UserLibShowMe" );
    for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
    {
      if( !types.count( type.id ) )
        continue;
      const WString label = UserLibraryWindow::whereUsedLabel( type.id );
      if( label.empty() )
        continue;
      WPushButton *btn = links->addNew<WPushButton>( label );
      btn->addStyleClass( "LinkBtn" );
      const string type_id = type.id;
      btn->clicked().connect( [this, type_id](){
        hide();  //out of the way of the tool it opens (closing deletes this dialog, later)
        UserLibraryWindow::showWhereUsed( type_id );
      } );
    }
  }//if( !types.empty() )

  WPushButton *open_lib = m_body->addNew<WPushButton>( WString::tr( "ul-open-library" ) );
  open_lib->addStyleClass( "LinkBtn" );
  open_lib->clicked().connect( [this](){
    InterSpec * const viewer = InterSpec::instance();
    if( viewer )
      viewer->showUserLibrary();
    hide();
  } );
}//installDone(...)

#endif //#if( USE_USER_LIBRARY )
