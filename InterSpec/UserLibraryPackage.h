#ifndef UserLibraryPackage_h
#define UserLibraryPackage_h
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
#include <cstddef>
#include <iosfwd>
#include <optional>

#include "InterSpec/UserLibrary.h"

class DetectorPeakResponse;

/** Packages as files: checking and unpacking a package someone gives you, installing
 (or updating) it, uninstalling it, and creating one.

 A package is a zip archive (extension `.ispkg`, or a `.zip` with `interspec-package.json` at its root)
 whose layout mirrors the user data folder.  Nothing here touches a widget; #inspect is slow (it reads
 every item) and is meant to run off the session thread.
 */
namespace UserLibraryPackage
{
  /** ".ispkg" */
  extern const char * const FileExtension;

  /** What an archive may hold. */
  struct Limits
  {
    size_t maxEntries = 2000;
    size_t maxTotalBytes = size_t(512)*1024*1024;  ///< uncompressed, all files
    size_t maxFileBytes = size_t(128)*1024*1024;   ///< uncompressed, any one file
    size_t maxCompressionRatio = 200;
    size_t maxPathLength = 160;                    ///< bytes, relative to the package folder
    size_t maxComponentLength = 100;               ///< bytes, any one folder or file name
  };//struct Limits


  /** Whether a dropped file is a package: by extension, or a `.zip` with the manifest at its root
   (reads only the zip's directory).  Never throws.
   */
  bool looksLikePackage( const std::string &displayName, const std::string &path );

  /** Checks an archive entry's name: no absolute path or drive, no "." or ".." component, no names
   Windows forbids, within the length limits, valid UTF-8.  Returns it '/' separated (no trailing
   '/'); throws std::runtime_error with an (English) reason otherwise.
   */
  std::string checkEntryName( const std::string &name, const Limits &limits );


  enum class ItemState : int
  {
    New,          ///< nothing installed has it, or anything it is named
    AlreadyHave,  ///< identical content already installed (same hash)
    Unchanged,    ///< an update, and this item is the same as the installed version's
    Update,       ///< an update, and this item changed
    Conflict,     ///< something installed has one of its names, with different content
    NeedsNewer,   ///< written for a newer InterSpec: installed, but inactive until an upgrade
    Unreadable    ///< could not be read; installed with the package, but flagged
  };//enum class ItemState


  /** Something already installed with a name an incoming item also defines. */
  struct NameConflict
  {
    std::string name;
    UserLibrary::Item existing;
    bool identical = false;           ///< same content: nothing to decide
    bool overridesBuiltIn = false;    ///< the existing item is InterSpec's own
    bool alreadyChosen = false;       ///< a choice is recorded (e.g. from an earlier version): it stands
    bool useIncomingByDefault = false;
  };//struct NameConflict


  struct IncomingItem
  {
    UserLibrary::FoundItem found;     ///< in the staging folder
    std::string typeId;
    std::string relPath;              ///< in the package, '/' separated
    std::string itemId;
    std::string sha256;
    std::string notes;                ///< untrusted plain text
    bool intendedOverride = false;
    bool hashMismatch = false;        ///< the manifest's sha256 does not match the content
    UserLibrary::Validation validation;
    ItemState state = ItemState::New;
    std::vector<NameConflict> conflicts;
  };//struct IncomingItem


  /** A file in the archive that is not installed. */
  struct SkippedFile
  {
    std::string relPath;
    bool refused = false;   ///< never accepted (API keys, the user database, logs); else just not a known type
  };//struct SkippedFile


  enum class Relation : int
  {
    NewPackage,            ///< no package with this id is installed
    Update,                ///< a newer version of an installed package
    Reinstall,             ///< the same version, with the same content
    SameVersionDifferent,  ///< the same version number, but different content
    Downgrade              ///< an older version than is installed
  };//enum class Relation


  /** A package checked and unpacked, ready to install or discard. */
  struct Inspection
  {
    std::string displayName;          ///< the file's name, as dropped
    std::string archiveSha256;        ///< for the recipient to check with the sender
    UserLibrary::Manifest manifest;   ///< from the archive, or made up when it has none
    bool hadManifest = false;
    std::string stagingDir;           ///< `packages/.incoming-*`; removed by #install or #discard

    std::vector<IncomingItem> items;
    std::vector<SkippedFile> skipped;

    Relation relation = Relation::NewPackage;
    std::optional<UserLibrary::Package> installed;   ///< the installed package with the same id
    std::vector<std::string> added, changed, removed; ///< item paths, compared with the installed version
    std::vector<std::string> locallyModified;         ///< installed files changed on disk since installed
    bool authorChanged = false;                       ///< author or organization differ from installed
    std::vector<UserLibrary::ManifestPreviousVersion> changesSince;  ///< newest first

    /** Non-empty: the package can not be installed (English, plain text). */
    std::string error;
  };//struct Inspection


  /** Unpacks a package into a staging folder, checks every entry and item, and compares it with
   what is installed.  Never throws; a problem that stops the install is put in `error` (and the
   staging folder is then already removed).
   */
  std::shared_ptr<Inspection> inspect( const std::string &archivePath,
                                       const std::string &displayName,
                                       const Limits &limits = Limits() );


  /** For each conflict the user decided: (type id, name) -> use the incoming item.  Conflicts not
   listed take their default (`NameConflict::useIncomingByDefault`), unless a choice is recorded.
   */
  typedef std::map<std::pair<std::string,std::string>,bool> Decisions;

  /** Installs (or updates) an inspected package: moves the staging folder into place (an update
   swaps out the old folder, rolling back on failure), then records it and the conflict choices in
   the history.  Throws std::exception on failure, leaving any installed version as it was.
   */
  UserLibrary::Package install( Inspection &inspection, const Decisions &decisions );

  /** Removes the staging folder of an inspection that will not be installed. */
  void discard( Inspection &inspection );

  /** Deletes an installed package's folder, and its record and the choices that named it. */
  void uninstall( const std::string &packageId );


  // ---------------------------------------------------------------------------------------------
  //  Creating packages
  // ---------------------------------------------------------------------------------------------

  struct CreateItem
  {
    UserLibrary::Item item;          ///< from UserLibrary::items(); the user's own, or a package's
    std::string notes;
    bool intendedOverride = false;

    /** For a file or folder outside the user data folder and packages (a folder added in the
     Detector Efficiency dialog's Locations, or a reference-spectra folder added in that tool): the
     folder its package path is taken relative to - normally the added folder's parent, so the folder's
     own name is kept ("drfs/<folder>/...").  Empty otherwise.
     */
    std::string root;

    /** For a detector efficiency saved in InterSpec's database (no file): it is written into the
     package as "drfs/<name>.drf.xml", and `item.found` is unused.  Null otherwise.
     */
    std::shared_ptr<const DetectorPeakResponse> databaseDrf;

    /** For another item with no file (e.g. a colour theme saved in InterSpec's database): its content,
     written into the package as "<type folder>/<generatedName>".  Null otherwise.
     */
    std::shared_ptr<const std::string> generatedData;
    std::string generatedName;

    /** How the created-package record refers to where the item came from; empty to use
     #sourceKey of `item`.  (E.g. "database:<hash>" or "location:<path>".)
     */
    std::string sourceKey;
  };//struct CreateItem

  struct CreateRequest
  {
    /** The package's details; `id` empty for a new package (one is made), or the id of a package
     made before, for its new version.  `items`, `created`, `created_with`, and (unless set)
     `previous_versions` are filled in.
     */
    UserLibrary::Manifest manifest;
    std::vector<CreateItem> items;
  };//struct CreateRequest

  struct PreflightIssue
  {
    std::string message;   ///< English, plain text
    bool blocking = true;
  };//struct PreflightIssue

  /** Problems with a request: missing name or version, no items, built-in items, duplicate names or
   paths, unreadable items, secrets, paths too long.  Reads every item; may be slow.
   */
  std::vector<PreflightIssue> preflight( const CreateRequest &request, const Limits &limits = Limits() );

  /** Where an item goes in a package, e.g. "drfs/My Detectors/HPGe 40%". */
  std::string packagePathFor( const UserLibrary::Item &item );
  std::string packagePathFor( const CreateItem &item );

  /** The files an item is made of (absolute path -> path in the package): its file, or every file
   in its folder, plus support files (template includes, a .par's DETECTOR.txt).  Empty for a
   database detector efficiency (see CreateItem::databaseDrf), which is written from memory.
   */
  std::vector<std::pair<std::string,std::string>> packageFilesFor( const UserLibrary::Item &item );
  std::vector<std::pair<std::string,std::string>> packageFilesFor( const CreateItem &item );

  /** A detector efficiency as the XML InterSpec writes and reads (".drf.xml"). */
  std::string drfXml( const DetectorPeakResponse &drf );

  /** Writes a package to `output` and records it as created (not installed) in the history, along
   with the author/organization/contact for next time.  Throws if a blocking #preflight issue
   remains, or on a read/write error.  Returns the manifest written.
   */
  UserLibrary::Manifest create( const CreateRequest &request, std::ostream &output,
                                const Limits &limits = Limits() );

  /** A suggested next version: "3" -> "4", "1.2" -> "1.3"; else the version with ".1" appended. */
  std::string nextVersion( const std::string &version );

  /** "user:<relPath>" or "package:<id>:<relPath>": how a created package records where an item came
   from (to keep its item_id in the next version, and spot files that moved or changed).
   */
  std::string sourceKey( const UserLibrary::Item &item );
  std::string sourceKey( const CreateItem &item );
}//namespace UserLibraryPackage

#endif //UserLibraryPackage_h
