#ifndef DrfLibrary_h
#define DrfLibrary_h
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
#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include <functional>

#include "InterSpec/DetectorPeakResponse.h"

struct DrfAppearance;
namespace DataBaseUtils{ class DbSession; }


/** The data behind the Detector Efficiency dialog's Library: every detector the user can pick -
 shipped with InterSpec, in folders on disk, saved in the database, or embedded in the spectrum -
 as one list, with what the list needs to show for each.

 Nothing here touches a widget or a `Wt::WApplication`, so it is safe on any thread and testable on
 its own.  The flow is:
   1. The session gathers the inputs (`ScanRequest`: locations from the preferences, the database
      session and user, the detectors to suggest for the spectrum).
   2. `scan` runs on a worker thread.  It walks the locations and the user's database rows, and
      reuses each detector's `Summary` from the process-wide cache unless a file (by modification
      time and size) or a database row (by hash and size) changed - so reopening the dialog only
      re-reads what changed, and never shows a stale summary.
   3. `buildRows` turns the result plus the dialog's filter state into the rows to display.
   4. `loadDrf` gives the full DetectorPeakResponse for an entry, when the user previews or uses it.
 */
namespace DrfLibrary
{
  /** Where a detector comes from, as the Library's source chips group them. */
  enum class Origin : int
  {
    BuiltIn,      ///< ships with InterSpec (data/GenericGadrasDetectors, data/common_drfs.tsv)
    Gadras,       ///< a GADRAS-format folder (Detector.dat [+ Efficiency.csv])
    Isocs,        ///< an ISOCS .ecc, or a characterization .par + DETECTOR.txt
    Angle,        ///< an ANGLE .outx / .detx
    RelEffTable,  ///< a relative-efficiency table (.tsv/.csv, one detector per row)
    UserMade,     ///< saved inside InterSpec: formulas, Make DRF, one-off imports, customized copies
    ThisSpectrum  ///< embedded in the foreground spectrum file, and nowhere else in the Library
  };//enum class Origin

  /** The Origin a detector's DrfSource implies, for detectors not found in a location. */
  Origin originForSource( const DetectorPeakResponse::DrfSource source );


  /** What sort of source-to-detector placements a detector can be used for. */
  enum class Range : int
  {
    FarFieldOnly,          ///< a far-field efficiency scaled by solid angle
    NearFieldOk,           ///< has a CeeLo response, so near-field and off-axis are modeled
    FixedGeometry,         ///< only the setup it was made for
    NeedsCharacterization  ///< a crystal geometry with no efficiency yet (needs a Monte-Carlo run)
  };//enum class Range


  /** The detector crystal, as best it can be told. */
  struct Crystal
  {
    /** An entry of `CeeLoUtils::crystalMaterialNames()` (e.g. "HPGe", "NaI", "LaBr3"), or empty. */
    std::string material;

    /** Whether the dimensions below are the recorded crystal geometry. */
    bool fromGeometry = false;
    bool isBox = false;
    double diameter_cm = 0.0;  ///< cylinder diameter, or box width
    double height_cm = 0.0;    ///< box only
    double length_cm = 0.0;    ///< crystal length (depth)

    /** The DRF's stated face diameter (the flat-disk model's); zero for fixed geometry. */
    double faceDiameter_cm = 0.0;
  };//struct Crystal


  /** Display facts for one detector, computed once from its DetectorPeakResponse.

   Holds only numbers, flags and the DRF's own text, not translated strings, so the same summary can
   serve any session in any language.
   */
  struct Summary
  {
    std::string name;
    std::string description;
    DetectorPeakResponse::DrfSource source = DetectorPeakResponse::DrfSource::UnknownDrfSource;
    DetectorPeakResponse::EffGeometryType geometryType = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;

    Range range = Range::FarFieldOnly;
    Crystal crystal;

    /** IEEE 325 relative efficiency (see `DetectorPeakResponse::relativeEfficiencyIeee325`);
     negative when not meaningful.
     */
    double relEffPercent = -1.0;

    /** Whether the name already states a percent efficiency ("HPGe 40%", "Falcon (18%)"), so the
     list need not repeat one.
     */
    bool nameStatesRelEff = false;

    bool hasEfficiency = false;      ///< `DetectorPeakResponse::isValid()`
    DetectorPeakResponse::EfficiencyFnctForm efficiencyForm = DetectorPeakResponse::kNumEfficiencyFnctForms;
    size_t efficiencyParamCount = 0; ///< coefficients, or energy/efficiency pairs
    size_t measuredPointCount = 0;   ///< points the curve was fit to (Make DRF), if recorded
    bool importedGrid = false;       ///< efficiency is a characterization file's grid (.par)
    bool hasResponse = false;        ///< a CeeLo response is attached

    bool hasFwhm = false;
    float fwhmAt662 = -1.0f;         ///< keV; negative when no FWHM

    bool hasGeometry = false;        ///< a crystal geometry is in effect (or a response models one)
    bool hasUncertainty = false;     ///< reports an efficiency uncertainty
    bool uncertaintyIsAssumed = false; ///< ... but only an assumed one (see DrfModifyCalc::UncertSummary)
    bool hasTotalEfficiency = false; ///< `hasAnyTotalEfficiencyInfo()` - usable for summing corrections

    double lowerEnergy = 0.0, upperEnergy = 0.0;  ///< keV; both zero when unstated

    uint64_t hash = 0;               ///< `hashValue()`
    uint64_t parentHash = 0;         ///< `parentHashValue()`; non-zero for a customized copy

    /** For a detector from a GADRAS folder or a relative-efficiency table: its hash under the other
     source it has historically been saved with (`DefaultGadrasDrf` <-> `UserAddedGadrasDrf`,
     `DefaultRelativeEfficiencyDrf` <-> `UserAddedRelativeEfficiencyDrf`), so database rows and
     customized copies made under either still find it.  Empty for other sources.
     */
    std::vector<uint64_t> sourceAliasHashes;

    std::shared_ptr<const DrfAppearance> appearance;

    static std::shared_ptr<const Summary> fromDrf( const DetectorPeakResponse &drf );
  };//struct Summary


  /** Why a detector is suggested for the foreground spectrum; a bit set. */
  enum SuggestReason : uint32_t
  {
    SuggestNone = 0x0,
    SuggestForSerialNumber = 0x1,  ///< the user chose it for this instrument's serial number
    SuggestForModel = 0x2,         ///< the user chose it for this instrument model
    SuggestInSpectrumFile = 0x4,   ///< embedded in the spectrum file
    SuggestDefaultForType = 0x8    ///< InterSpec's default for this detector type
  };//enum SuggestReason


  /** A detector the controller wants suggested (or added, if the Library does not already have it). */
  struct Suggestion
  {
    std::shared_ptr<const DetectorPeakResponse> drf;
    uint32_t reasons = SuggestNone;
  };//struct Suggestion


  /** A place detectors are read from: a folder (searched recursively) or a single file. */
  /** What a location is, for how it is shown and managed. */
  enum class LocationRole : int
  {
    User,            ///< a folder or file the user added
    GenericGadras,   ///< data/GenericGadrasDetectors
    ExtraGadras,     ///< data/OUO_GadrasDetectors (web deployments)
    BuiltInTable,    ///< data/common_drfs.tsv
    KeptImports,     ///< the user-data `drfs` folder, where kept imports are copied
    UserData,        ///< the user-data folder (GADRAS folders and tables put there by hand)
    GadrasInstall,   ///< the standard GADRAS install folder
    IsocsInstall,    ///< the standard ISOCS (Genie 2000) data folder
    AngleInstall,    ///< ANGLE's standard folder
    Package          ///< the `drfs` folder of a user library package
  };//enum class LocationRole

  struct Location
  {
    /** Stable key; also how a user location is stored.  Usually the path. */
    std::string id;
    std::string path;
    LocationRole role = LocationRole::User;

    bool isFolder = true;
    bool isAuto = false;    ///< found automatically (built in, a standard install) vs added by the user
    bool builtIn = false;   ///< ships with InterSpec; its detectors get Origin::BuiltIn
    bool hidden = false;    ///< "Show in Library" unticked - still scanned, but its rows are hidden

    /** Which file formats to look for; a bit set of #LocationFormats. */
    uint32_t formats = 0x7;

    /** What to call it, when its role does not say (a package's name; untrusted plain text). */
    std::string label;

    /** Folders below `path` not to search (e.g. the user-data `packages` folder). */
    std::vector<std::string> excludePaths;

    /** Detectors found at these paths (a file, or a GADRAS folder) are hidden, like a hidden location's;
     normalized with SpecUtils::lexically_normalize_path.  (Items disabled in the user library.)
     */
    std::set<std::string> hiddenItems;
  };//struct Location

  enum LocationFormats : uint32_t
  {
    FormatGadras = 0x1,  ///< folders with a GADRAS Detector.dat [+ Efficiency.csv]
    FormatTable  = 0x2,  ///< relative-efficiency tables (.tsv, and .csv other than Efficiency.csv)
    FormatOther  = 0x4,  ///< .ecc, .par + DETECTOR.txt, .outx/.detx, and DRF .xml
    FormatAll    = 0x7
  };//enum LocationFormats


  /** How an entry's full DRF is (re)loaded. */
  enum class EntryKind : int
  {
    GadrasFolder,  ///< `path` is a folder holding Detector.dat [+ Efficiency.csv]
    TableRow,      ///< row `record` of the relative-efficiency table at `path`
    ImportFile,    ///< a file DrfImport reads (`companionPath` set for a .par + DETECTOR.txt pair)
    DbRow,         ///< database row `dbId`
    InMemory       ///< `inMemoryDrf` (e.g., the spectrum's own DRF)
  };//enum class EntryKind


  /** One detector in the Library. */
  struct Entry
  {
    /** Unique and stable for the session: "gadras:<path>", "row:<path>|<n>", "file:<path>",
     "db:<id>", "mem:<hash>".
     */
    std::string id;

    EntryKind kind = EntryKind::InMemory;
    Origin origin = Origin::UserMade;

    std::string locationId;     ///< the Location it was found in; empty if none
    std::string path;           ///< folder or file; empty for database rows and in-memory DRFs
    std::string companionPath;  ///< the DETECTOR.txt of a .par, if any
    size_t record = 0;          ///< table row

    long long dbId = -1;        ///< the user's database row for this detector, if there is one
    int64_t lastUsedUtc = 0;    ///< from the database; zero if never used

    std::shared_ptr<const Summary> summary;  ///< null only while it is being read, or on error
    bool pending = false;       ///< being read (see the `progress` callback of #scan)
    std::string error;          ///< why it could not be read

    bool hidden = false;        ///< its location is hidden from the Library, or it is disabled
    uint32_t suggestReasons = SuggestNone;

    /** Where it came from, when that is more than its Origin says: the package's name. */
    std::string sourceLabel;

    /** Every DRF hash this entry stands for: its own, plus those of database rows merged into it
     (which may differ in DrfSource only).  Customized copies find their parent through these.
     */
    std::vector<uint64_t> knownHashes;

    std::shared_ptr<const DetectorPeakResponse> inMemoryDrf;  ///< for EntryKind::InMemory
  };//struct Entry


  /** What a scan found at one location. */
  struct LocationStatus
  {
    std::string id;
    bool found = false;           ///< the folder or file exists
    size_t numDetectors = 0;      ///< detectors listed from it (including hidden)
    size_t numUnreadable = 0;     ///< detector-looking files that could not be read, or need input
  };//struct LocationStatus


  /** A detector-looking file (or GADRAS folder) from which no detector could be read. */
  struct UnreadableItem
  {
    std::string path;
    bool isFolder = false;
    std::string error;  ///< may be empty
  };//struct UnreadableItem

  struct Snapshot
  {
    std::vector<Entry> entries;
    std::vector<LocationStatus> locations;
    std::vector<UnreadableItem> unreadable;
  };//struct Snapshot


  struct ScanRequest
  {
    std::vector<Location> locations;

    /** The user's database rows are included when both are set. */
    std::shared_ptr<DataBaseUtils::DbSession> sql;
    long long dbUserId = -1;

    std::vector<Suggestion> suggestions;

    /** Set to abandon the scan (e.g., the dialog closed); `scan` then returns early. */
    std::shared_ptr<const std::atomic<bool>> cancel;

    /** If set, called (on the scanning thread) once the quick pass is done and some detectors
     still need reading: entries whose summary is not cached yet have `pending` set.  Lets the
     dialog show the list while the rest is read.
     */
    std::function<void(const Snapshot &)> progress;
  };//struct ScanRequest


  /** Walks the locations and database rows, reusing cached summaries where nothing changed.
   Blocking; run it off the session thread.  Never throws.
   */
  Snapshot scan( const ScanRequest &request );


  /** The standard locations for this build: InterSpec's own detectors, the user-data `drfs` folder
   (and any folder in the user-data directory holding GADRAS detectors), the `drfs` folder of each
   installed user library package, and - on Windows desktop builds - the default data folders of
   GADRAS, ISOCS and ANGLE.  Items disabled in the user library are in each location's `hiddenItems`.  Those are returned whether or not
   they exist; a scan sets `LocationStatus::found` false for a missing one, so the dialog can offer
   to locate it.

   Does not include user-added locations (those come from the preferences).
   */
  std::vector<Location> standardLocations();


  /** The user's changes to the Library's locations, as stored in a preference (see
   #parseLocationSettings): folders and files they added, standard locations they hid, and
   standard locations they pointed somewhere else ("Locate…" for one not found).
   */
  struct LocationSettings
  {
    std::vector<Location> userLocations;      ///< in the order added
    std::set<std::string> hiddenStandard;     ///< ids of hidden standard locations
    std::map<std::string,std::string> standardPaths;  ///< standard location id -> path chosen

    bool operator==( const LocationSettings &rhs ) const;
  };//struct LocationSettings

  /** Reads the preference string: ';'-separated entries, each
   - `path`, a folder or file the user added (`!path` if hidden from the Library), optionally as
     `gadras:path` for one holding GADRAS folders only (folders that came from the old "GadrasDRFPath"
     preference);
   - `@id` / `!@id` / `@id=path` / `!@id=path` for a standard location (by its `Location::id`): its
     path, if the user located it, and whether it is hidden.
   Whether a user location is a folder or a file is decided from the file system when scanned.
   Never throws; empty and malformed entries are skipped.
   */
  LocationSettings parseLocationSettings( const std::string &pref );

  /** The inverse of #parseLocationSettings.  Starts with a "#1" version marker (entries starting with
   '#' are ignored when read), so the result is never empty - telling a saved "nothing added" from
   never having been saved. */
  std::string serializeLocationSettings( const LocationSettings &settings );

  /** The locations to scan: `standard` with the user's hidden flags and paths applied, then the
   user's own locations (a user location duplicating a standard one, or another user location, is
   left out).  User locations' `isFolder` is set from the file system.
   */
  std::vector<Location> applyLocationSettings( std::vector<Location> standard,
                                               const LocationSettings &settings );


  /** The full DRF for an entry - from the process-wide cache when it was loaded recently.
   Blocking (reads files or the database); throws std::exception on failure.
   */
  std::shared_ptr<const DetectorPeakResponse> loadDrf( const Entry &entry,
                                                      std::shared_ptr<DataBaseUtils::DbSession> sql );


  /** Drops cached summaries and DRFs for a path (a file, or everything under a folder). */
  void invalidatePath( const std::string &path );

  /** Drops the cached summary for a database row. */
  void invalidateDbRow( const long long dbId );

  /** Drops everything cached; mainly for tests. */
  void clearCache();

  /** Cache statistics, for tests: how many detector summaries were computed (not reused). */
  size_t numSummariesComputed();


  /** Combines detectors found in locations, the user's database rows, and suggestions into the
   Library's entries (the last step of #scan; exposed for testing).

   - A database row whose stored hash is one of a located detector's `knownHashes` (its hash, or
     its `Summary::sourceAliasHashes`) is folded into it: the entry gets the row's `dbId`, last-used
     time and hash.  So is a row from a GADRAS folder or table (by its DrfSource; not a customized
     copy) that matches nothing by hash but has the name of a detector of the same kind found on
     disk - a usage record of the same detector from before its file, or InterSpec's modeling of it,
     changed.
   - Remaining rows become entries of their own: `Origin::UserMade` for customized copies and the
     sources InterSpec makes, else #originForSource.
   - Each suggestion marks the entry it matches the same way; one with no match is added as an
     `EntryKind::InMemory` entry (`Origin::ThisSpectrum` if it came from the spectrum file).
   */
  std::vector<Entry> mergeEntries( std::vector<Entry> located,
                                   std::vector<Entry> dbRows,
                                   const std::vector<Suggestion> &suggestions );


  /** The dialog's Library filter. */
  struct Filter
  {
    bool allOrigins = true;               ///< the "All" chip
    Origin origin = Origin::BuiltIn;      ///< used when !allOrigins
    bool withGeometry = false;            ///< the "With geometry" chip
    std::string query;                    ///< search text; case-insensitive, all words must match

    /** At most this many detectors (families) under "Recently used"; the rest are listed under
     "Everything else".
     */
    size_t maxRecent = 10;
  };//struct Filter


  enum class RowType : int { GroupHeader, Detector, CopiesToggle };
  enum class Group : int { Suggested, Recent, Everything, All };


  /** One row of the Library list. */
  struct Row
  {
    RowType type = RowType::Detector;
    Group group = Group::All;   ///< the group this row is in (and, for a header, which group)

    size_t entry = 0;           ///< Detector: its index in the entries; CopiesToggle: the parent's
    bool nested = false;        ///< Detector: a customized copy listed under its parent
    bool copyOfHidden = false;  ///< Detector: a copy listed on its own, as "Copy of ..."
    size_t numCopies = 0;       ///< CopiesToggle: how many copies it shows or hides
    bool expanded = false;      ///< CopiesToggle: whether the copies are shown
  };//struct Row


  /** Optional text to search, per entry, in addition to its name, description, crystal and path -
   e.g., the translated source label.
   */
  typedef std::function<std::string(const Entry &)> ExtraSearchText;

  /** The rows of the Library list: "Suggested for this spectrum", then "Recently used", then
   "Everything else" (alphabetical) - or a single unnamed group when the first two are empty.

   Customized copies (by `Summary::parentHash`) are listed under their parent, behind a
   CopiesToggle row; they are shown when `expanded[parent id]` is true, or - when not set - when
   the picked or in-use detector is one of them.  A copy whose parent is filtered out or missing is
   listed on its own, with `copyOfHidden` set.

   Hidden entries and entries without a summary that are not pending are left out.
   */
  std::vector<Row> buildRows( const std::vector<Entry> &entries,
                              const Filter &filter,
                              const std::map<std::string,bool> &expanded,
                              const std::string &pickedId,
                              const std::string &inUseId,
                              const ExtraSearchText &extraText = nullptr );

  /** Whether an entry passes a filter (origin chip, With geometry, and search text). */
  bool matchesFilter( const Entry &entry, const Filter &filter, const ExtraSearchText &extraText = nullptr );

  /** Index of the entry a customized copy was made from, or -1 (looks through `knownHashes`). */
  int parentIndex( const std::vector<Entry> &entries, const size_t child );
}//namespace DrfLibrary

#endif //DrfLibrary_h
