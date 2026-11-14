#ifndef UserLibrary_h
#define UserLibrary_h
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
#include <set>
#include <string>
#include <vector>
#include <cstdint>
#include <optional>
#include <functional>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json_fwd.hpp"

namespace inja{ class Environment; }


/** The user library: the user data InterSpec reads by type (detector efficiencies, isotopics presets,
 report templates, reference spectra, LLM skills, ...), wherever it comes from - InterSpec itself, an
 installed package, or the user's own files - and which item answers when a tool looks one up by name.

 - Packages live in `<user data>/packages/<folder>/`, laid out like the user data folder, and InterSpec
   never edits the files in them.  Installing an update replaces the folder; uninstalling deletes it.
 - `<user data>/user_library.json` (the "history") records what was installed and created, and the
   user's choices: which of several same-named items is used, and which items are disabled.  It is
   advisory: it is reconciled against what is on disk whenever it is read.
 - Precedence for same-named items, absent a recorded choice: the user's own files, then packages
   (earliest installed first), then InterSpec's own.

 Nothing here touches a widget, so every function may be called from any thread.  Without a user data
 folder (web deployments, or a process whose user data folder is the static data folder), there are
 no packages and no history, and only built-in and user items (if any) are found.
 */
namespace UserLibrary
{
  /** Whether this process has a user library: a user data folder, distinct from the static data folder. */
  bool available();

  /** `<user data>/packages`; throws if not #available. */
  std::string packagesDirectory();

  /** `<user data>/user_library.json`; throws if not #available. */
  std::string historyPath();

  /** The name of the manifest file at the root of a package. */
  extern const char * const ManifestFileName;


  // ---------------------------------------------------------------------------------------------
  //  Asset types
  // ---------------------------------------------------------------------------------------------

  /** How items of a type with the same name relate. */
  enum class Kind : int
  {
    List,         ///< the user picks from a list: same-named items are all listed, labelled by source
    MergeByName,  ///< only one item with a given name is used
    OneActive     ///< only one is active, by explicit activation (e.g. nuclear data; not used yet)
  };//enum class Kind

  /** When a change to items of a type takes effect. */
  enum class ApplyMode : int
  {
    Live,     ///< next time a tool lists or looks up items
    Reload,   ///< after a reload hook runs
    Restart   ///< after InterSpec restarts
  };//enum class ApplyMode

  enum class Status : int
  {
    Ok,
    NeedsNewer,  ///< written by a newer InterSpec; installed, but unused until InterSpec is upgraded
    Unreadable
  };//enum class Status

  /** The result of reading an item with the reader its tool uses. */
  struct Validation
  {
    Status status = Status::Unreadable;
    std::string message;              ///< why, for NeedsNewer or Unreadable (English, plain text)
    std::vector<std::string> names;   ///< the identity names the item defines (a DRF table: many)
  };//struct Validation

  /** An item found on disk: a file, or a folder that is one item (a GADRAS detector, an LLM skill). */
  struct FoundItem
  {
    std::string path;                  ///< absolute
    bool isFolder = false;
    std::vector<std::string> support;  ///< other files it needs (a .par grid's DETECTOR.txt); absolute
  };//struct FoundItem

  /** Where items of a type are found, and how they are read.  See #assetTypes. */
  struct AssetType
  {
    std::string id;         ///< also the `type` of a manifest item; e.g. "drf"
    std::string titleKey;   ///< i18n key (UserLibrary.xml), e.g. "ul-type-drf"
    std::string folder;     ///< folder in a package and in the user data folder, '/' separated

    /** Older folder names also read, in packages and in the user data folder. */
    std::vector<std::string> legacyFolders;

    /** Older single files also read from the user data folder (e.g. "add_ref_line.xml"), each one item. */
    std::vector<std::string> legacyUserFiles;

    /** InterSpec's own items: folders or files, relative to the static data folder (or the docroot,
     for paths starting with "InterSpec_resources/").
     */
    std::vector<std::string> builtInPaths;

    Kind kind = Kind::List;
    ApplyMode apply = ApplyMode::Live;
    bool sensitive = false;              ///< needs explicit confirmation to install, and shows its text
    bool caseInsensitiveNames = false;   ///< names come from file names

    /** Optional: two names are the same when this gives the same key - e.g. materials, which the
     material database matches ignoring case, spaces, '-' and '_'.  Overrides #caseInsensitiveNames.
     */
    std::function<std::string( const std::string &name )> nameKey;

    /** The items under a folder (or the item a single file is).  Never throws. */
    std::function<std::vector<FoundItem>( const std::string &root )> findItems;

    /** Reads an item with the real reader.  Never throws.  May be slow (parses the file). */
    std::function<Validation( const FoundItem &item )> validate;

    /** The identity names of an item; cheap where possible (file names), cached where not (DRFs). */
    std::function<std::vector<std::string>( const FoundItem &item )> names;

    /** Other files an item uses (e.g. a template's `{% include %}` targets), as absolute paths.
     Optional.
     */
    std::function<std::vector<std::string>( const FoundItem &item )> dependencies;
  };//struct AssetType

  /** Every type the user library knows, in display order. */
  const std::vector<AssetType> &assetTypes();

  /** The type with this id, or nullptr. */
  const AssetType *assetType( const std::string &typeId );

  /** The type whose folder (or legacy folder) a relative path (e.g. "drfs/a/Detector.dat") is in, or
   nullptr.  If found, `folder` is set to the matching folder (e.g. "drfs").
   */
  const AssetType *typeForRelativePath( const std::string &relPath, std::string &folder );


  // ---------------------------------------------------------------------------------------------
  //  Sources, items, and name lookups
  // ---------------------------------------------------------------------------------------------

  enum class SourceKind : int { BuiltIn, Package, User };

  struct Source
  {
    SourceKind kind = SourceKind::User;
    std::string packageId;   ///< for SourceKind::Package
    std::string label;       ///< the package's name (untrusted plain text); empty for others

    /** "builtin", "user", or "package:<id>": how choices refer to a source. */
    std::string key() const;

    bool operator==( const Source &rhs ) const;
  };//struct Source

  /** A place items of a type are looked for. */
  struct Root
  {
    Source source;

    /** Absolute folder, or file (InterSpec's own DRF table). */
    std::string path;

    /** The folder relative paths of items under this root are given from: the user data folder, the
     package folder, or the static data folder.
     */
    std::string base;

    /** This root is the user data folder itself; DRFs only (GADRAS folders and tables put there by
     hand - see DrfLibrary::LocationRole::UserData).
     */
    bool userDataTop = false;
  };//struct Root

  /** Where items of a type are looked for, highest default precedence first: the user's own folders,
   then each installed package (earliest installed first), then InterSpec's own.  Roots that do not
   exist are left out.
   */
  std::vector<Root> searchRoots( const std::string &typeId );

  /** One item, from one source. */
  struct Item
  {
    std::string typeId;
    Source source;
    FoundItem found;

    /** Relative to the source's base folder (`Root::base`), '/' separated; e.g. "drfs/Det A.xml".  */
    std::string relPath;

    std::string itemId;               ///< the package manifest's `item_id`, if any
    std::vector<std::string> names;
    bool disabled = false;

    /** A package item written for a newer InterSpec: listed only with `includeDisabled`, and never
     answers lookups (it starts working once InterSpec is upgraded).
     */
    bool needsNewer = false;

    /** False if every one of its names is answered by some other item. */
    bool inUse = true;
  };//struct Item

  /** Every item of a type, in precedence order (see #searchRoots).  Pickers pass
   `includeDisabled = false`, which also leaves out package items that need a newer InterSpec.
   Names are filled in, and `inUse` reflects the recorded choices.
   */
  std::vector<Item> items( const std::string &typeId, const bool includeDisabled = false );

  /** The item that answers a lookup by name: the recorded choice if that item is still there and not
   disabled, else the first by precedence.  Disabled items never answer.
   */
  std::optional<Item> resolve( const std::string &typeId, const std::string &name );

  /** Like #resolve, but chooses among `candidates` (as from #items), without listing again. */
  std::optional<Item> resolve( const std::string &typeId, const std::string &name,
                               const std::vector<Item> &candidates );

  /** The type's validation of an item, cached by the item's (and its support files') modification
   time and size - so it is re-checked when the files change, or in a newer InterSpec.
   */
  Validation validateItem( const std::string &typeId, const FoundItem &item );

  /** A name as the type compares names (see AssetType::nameKey, AssetType::caseInsensitiveNames). */
  std::string nameKey( const AssetType &type, const std::string &name );

  /** Whether two names are the same for a type (`nameKey` of both equal). */
  bool sameName( const AssetType &type, const std::string &lhs, const std::string &rhs );


  // ---------------------------------------------------------------------------------------------
  //  Package manifest
  // ---------------------------------------------------------------------------------------------

  struct ManifestItem
  {
    std::string itemId;     ///< stable across versions
    std::string path;       ///< package relative, '/' separated; a folder item has no trailing '/'
    std::string type;       ///< an AssetType::id
    std::string sha256;
    std::string notes;

    /** When the item was re-shared from someone else's package. */
    std::string derivedFromName, derivedFromId, derivedFromVersion;

    /** Meant to replace a built-in or same-named item: pre-selects it in the conflict dialog. */
    bool intendedOverride = false;
  };//struct ManifestItem

  /** An earlier version's change notes, carried forward so an update can show every version's. */
  struct ManifestPreviousVersion
  {
    std::string version, created, changes;
  };//struct ManifestPreviousVersion

  /** `interspec-package.json`.  Every text field is untrusted plain text. */
  struct Manifest
  {
    static const int sm_formatVersion = 1;

    int formatVersion = sm_formatVersion;
    std::string id, name, version, created, createdWith, minInterSpecVersion;
    std::string author, organization, contact;
    std::string description, notes, changes, handling, reviewBy;
    std::vector<ManifestPreviousVersion> previousVersions;
    std::vector<ManifestItem> items;

    /** Throws on invalid JSON or a missing/empty required field (`format_version`, `id`, `name`,
     `version`, and each item's `path` and `type`).  A `format_version` newer than this InterSpec
     reads is not an error here - see #needsNewerInterSpec.
     */
    static Manifest fromJson( const std::string &json );
    static Manifest fromJsonObject( const nlohmann::json &json );

    nlohmann::json toJsonObject() const;
    std::string toJson() const;

    /** The manifest item for a package-relative path, or nullptr. */
    const ManifestItem *itemForPath( const std::string &relPath ) const;

    /** Whether `format_version` or `min_interspec_version` is newer than this InterSpec. */
    bool needsNewerInterSpec() const;

    /** Whether `review_by` is a date (YYYY-MM-DD) before today (UTC). */
    bool reviewDue() const;
  };//struct Manifest


  // ---------------------------------------------------------------------------------------------
  //  Installed and created packages (from disk, and the history file)
  // ---------------------------------------------------------------------------------------------

  struct InstalledVersion
  {
    std::string installed;       ///< ISO 8601 UTC
    std::string archiveSha256;
    Manifest manifest;
  };//struct InstalledVersion

  struct Package
  {
    std::string id;
    std::string folder;          ///< folder name in packagesDirectory()
    std::string path;            ///< absolute
    Manifest manifest;           ///< as on disk now
    std::string firstInstalled;  ///< ISO 8601 UTC; empty if not recorded
    std::vector<InstalledVersion> installs;          ///< oldest first
    std::map<std::string,std::string> files;         ///< package-relative path -> sha256, as installed
    std::set<std::string> disabled;                  ///< item_id (or path, for items without one)
    bool recorded = true;        ///< false: on disk, but not in the history file ("added manually")
  };//struct Package

  /** The installed packages (folders in packagesDirectory() with a readable manifest), earliest
   installed first.  Cached; re-read when the folder or history file changes.
   */
  std::vector<Package> installedPackages();
  std::optional<Package> installedPackage( const std::string &packageId );

  /** Whether `path` is in packagesDirectory(), comparing real paths too (so a symbolic link in the user
   data path does not hide it).  Packaged templates must be rendered confined to their own folder.
   */
  bool isInPackages( const std::string &path );

  /** The installed package whose folder holds `path` (a file or folder in it), if any. */
  std::optional<Package> packageForPath( const std::string &path );

  struct CreatedVersion
  {
    std::string created;
    Manifest manifest;

    /** item_id -> where it came from: "user:<relPath>" or "package:<id>:<relPath>". */
    std::map<std::string,std::string> sources;
  };//struct CreatedVersion

  /** Packages this user created, by package id; versions oldest first. */
  std::map<std::string,std::vector<CreatedVersion>> createdPackages();

  struct LooseImport
  {
    std::string relPath, originalName, imported, sha256;
  };//struct LooseImport

  std::vector<LooseImport> looseImports();

  /** Records a single file (or folder) kept in the user data folder from a drop or import. */
  void recordLooseImport( const std::string &absPath, const std::string &originalName );

  struct Creator
  {
    std::string author, organization, contact;
  };//struct Creator

  /** Author/organization/contact of the package the user last created. */
  Creator lastCreator();


  // ---------------------------------------------------------------------------------------------
  //  Choices and disabled items
  // ---------------------------------------------------------------------------------------------

  struct Choice
  {
    std::string source;   ///< Source::key()
    std::string sha256;   ///< of the chosen item when chosen, to flag an update that changed it
    std::string chosen;   ///< ISO 8601 UTC
  };//struct Choice

  /** Recorded choices, by "<type id>:<name>". */
  std::map<std::string,Choice> choices();

  /** Makes `item` answer lookups by `name` (normally one of its names). */
  void setChoice( const std::string &typeId, const std::string &name, const Item &item );

  /** Back to default precedence for `name`. */
  void clearChoice( const std::string &typeId, const std::string &name );

  /** Hides an item from pickers and lookups (package items and the user's own; not built-in). */
  void setDisabled( const Item &item, const bool disabled );

  /** Absolute paths of the disabled items of a type; for tools that walk folders themselves. */
  std::set<std::string> disabledPaths( const std::string &typeId );

  /** Deletes one of the user's own items (a file, or a folder that is one item).  Refuses (throws) for
   anything not the user's, and for the user data folder itself, a type's folder, or a folder holding
   either or the packages.  A GADRAS detector folder loses only the files directly in it (and is
   removed if then empty), so detectors in folders below it are kept.
   */
  void deleteUserItem( const Item &item );


  // ---------------------------------------------------------------------------------------------
  //  History file
  // ---------------------------------------------------------------------------------------------

  /** Re-reads the history file, lets `fcn` change it, and writes it back atomically (a temporary
   file then a rename), then bumps #generation.  Throws if the file can not be written.  Serialized
   within this process; a rare race with another process is accepted.
   */
  void modifyHistory( const std::function<void( nlohmann::json &history )> &fcn );

  /** The history file's contents (an empty object if missing or unreadable). */
  nlohmann::json readHistory();

  /** Bumped by every install, uninstall, choice, or disable change in this process (and when a
   change by another process is noticed); a tool compares it to decide whether to re-list.
   */
  uint64_t generation();

  /** Forgets cached package listings and bumps #generation. */
  void notifyChanged();

  /** Tidies packagesDirectory() at startup: an `.old-*` folder left by an update that never finished is
   put back if its package's folder is missing; other `.incoming-*` and `.old-*` folders are removed.
   Folders less than an hour old are left alone (another InterSpec may be using them).
   */
  void cleanupStaleFolders();


  // ---------------------------------------------------------------------------------------------
  //  Utilities
  // ---------------------------------------------------------------------------------------------

  std::string sha256OfData( const void *data, const size_t size );

  /** Lower-case hex SHA-256 of a file's contents; throws if it can not be read. */
  std::string sha256OfFile( const std::string &path );

  /** A file item: its SHA-256.  A folder item: the SHA-256 of its sorted
   "<relative path>\t<file sha256>\n" lines.  Throws if a file can not be read.
   */
  std::string sha256OfItem( const FoundItem &item );

  /** Compares dotted versions ("1.10" > "1.9"); -1, 0, or 1.  Falls back to string comparison
   (setting `*numeric` false) when either is not dotted numbers.
   */
  int compareVersions( const std::string &lhs, const std::string &rhs, bool *numeric = nullptr );

  /** The version of this InterSpec, e.g. "1.0.14". */
  std::string interSpecVersion();

  std::string nowIsoUtc();
  std::string newUuid();

  /** The key two names collide under on a case-insensitive, normalization-insensitive file system:
   ASCII lower case, with combining marks (U+0300 to U+036F) removed.
   */
  std::string collisionKey( const std::string &utf8 );

  /** Whether a package-relative path is one InterSpec never exports or accepts (API keys, the user
   database, logs, ...).
   */
  bool isSecretPath( const std::string &relPath );

  /** Restricts the `{% include %}`/`{% extends %}` of the Inja template at `templatePath` to files in its
   own folder (and below): those are registered with `env` by their relative names, and looking up any
   other file is turned off - so a template from a package can not pull other files (e.g. API keys)
   into a report.  The built-in partials registered by name (e.g. "default-act-fit-html-results") still
   work.  Throws if a file can not be read.
   */
  void confineTemplateIncludes( inja::Environment &env, const std::string &templatePath );

  /** The `{% include %}`/`{% extends %}`/`{% import %}` targets of a template (or the files it includes)
   that are not files beside it - which a package would not carry, and a packaged template may not use.
   */
  std::vector<std::string> templateIncludesOutsideFolder( const std::string &templatePath );

  /** `path` relative to `base`, '/' separated; empty if not under it. */
  std::string relativePath( const std::string &base, const std::string &path );

  /** A name usable as a folder or file name on every platform, from untrusted text. */
  std::string sanitizeFileName( const std::string &name, const std::string &fallback );
}//namespace UserLibrary

#endif //UserLibrary_h
