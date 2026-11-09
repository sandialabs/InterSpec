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
#include <list>
#include <mutex>
#include <regex>
#include <tuple>
#include <cassert>
#include <cctype>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <filesystem>

#if( defined(_WIN32) && !BUILD_FOR_WEB_DEPLOYMENT )
// As in AppUtils.cpp: keep Windows.h from defining min/max, or pulling in the old winsock.h.
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif
#include <Windows.h>
#include <Shlobj.h>  //for SHGetKnownFolderPath
#endif

#include <Wt/Dbo/Dbo.h>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "io/DetectorResponse.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfImport.h"
#include "InterSpec/DrfLibrary.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/DrfModifyCalc.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
  /** Folders deeper than this below a location are not searched. */
  const int ns_max_folder_depth = 12;

  /** Stop walking a location after this many files (a mis-chosen location, like a home folder). */
  const size_t ns_max_files_per_location = 20000;

  /** Full DRFs kept in memory for quick preview/use (a .par grid's response can be several MB). */
  const size_t ns_max_cached_drfs = 16;


  std::filesystem::path to_fs_path( const string &utf8 )
  {
#ifdef _WIN32
    return std::filesystem::path( SpecUtils::convert_from_utf8_to_utf16(utf8) );
#else
    return std::filesystem::path( utf8 );
#endif
  }//to_fs_path(...)


  /** A file's identity for cache validation: modification time and size. */
  struct FileStamp
  {
    string path;
    bool exists = false;
    int64_t mtime = 0;
    uint64_t size = 0;

    bool operator==( const FileStamp &rhs ) const
    {
      return (path == rhs.path) && (exists == rhs.exists) && (mtime == rhs.mtime) && (size == rhs.size);
    }
  };//struct FileStamp


  FileStamp stamp_file( const string &path )
  {
    FileStamp s;
    s.path = path;

    std::error_code ec;
    const std::filesystem::path p = to_fs_path( path );
    s.exists = std::filesystem::is_regular_file( p, ec ) && !ec;
    if( !s.exists )
      return s;

    s.size = static_cast<uint64_t>( std::filesystem::file_size( p, ec ) );
    if( ec )
      s.size = 0;

    ec.clear();
    const std::filesystem::file_time_type t = std::filesystem::last_write_time( p, ec );
    if( !ec )
      s.mtime = static_cast<int64_t>( t.time_since_epoch().count() );

    return s;
  }//stamp_file(...)


  string lower_ascii( string s )
  {
    SpecUtils::to_lower_ascii( s );
    return s;
  }


#if( defined(_WIN32) && !BUILD_FOR_WEB_DEPLOYMENT )
  /** The user's Documents folder (wherever it was redirected to, e.g. OneDrive); empty if unknown. */
  string documents_dir()
  {
    PWSTR path = nullptr;  // freed with CoTaskMemFree, whether or not the call succeeds
    const HRESULT hr = ::SHGetKnownFolderPath( FOLDERID_Documents, KF_FLAG_DEFAULT, nullptr, &path );
    const string answer = (SUCCEEDED(hr) && path) ? SpecUtils::convert_from_utf16_to_utf8( path ) : string();
    ::CoTaskMemFree( path );
    return answer;
  }//documents_dir()
#endif


  /** A unit of reading: one GADRAS folder, one table, or one import file (+ companion). */
  struct WorkItem
  {
    DrfLibrary::EntryKind kind = DrfLibrary::EntryKind::ImportFile;
    string key;            ///< cache key, and entry-id base
    string path, companion;
    vector<string> stampPaths;
    size_t locationIndex = 0;
  };//struct WorkItem


  /** What reading a WorkItem gave. */
  struct ItemResult
  {
    vector<FileStamp> stamps;
    vector<shared_ptr<const DrfLibrary::Summary>> summaries;
    size_t numUnreadable = 0;
    string error;
  };//struct ItemResult


  struct CachedDbRow
  {
    int64_t hash = 0;
    long long extraLength = 0;
    shared_ptr<const DrfLibrary::Summary> summary;
  };//struct CachedDbRow


  /** The process-wide cache.  Holds only immutable data, so it can be shared by every session and
   survives "Clear Session".  Reading happens outside the lock; two sessions scanning the same
   folder at once may both read a file, and the last to finish wins - harmless.
   */
  class Cache
  {
  public:
    static Cache &instance()
    {
      static Cache s_cache;
      return s_cache;
    }

    bool findItem( const string &key, const vector<FileStamp> &stamps, ItemResult &result )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      const auto pos = m_items.find( key );
      if( (pos == end(m_items)) || !(pos->second.stamps == stamps) )
        return false;
      result = pos->second;
      return true;
    }

    void storeItem( const string &key, const ItemResult &result )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      m_items[key] = result;
    }

    shared_ptr<const DrfLibrary::Summary> findRow( const long long id, const int64_t hash,
                                                   const long long extra_length )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      const auto pos = m_rows.find( id );
      if( (pos == end(m_rows)) || (pos->second.hash != hash) || (pos->second.extraLength != extra_length) )
        return nullptr;
      return pos->second.summary;
    }

    void storeRow( const long long id, const CachedDbRow &row )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      m_rows[id] = row;
    }

    shared_ptr<const DetectorPeakResponse> findDrf( const string &key, const uint64_t hash )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      for( auto iter = begin(m_drfs); iter != end(m_drfs); ++iter )
      {
        if( (iter->first == key) && iter->second && (iter->second->hashValue() == hash) )
        {
          m_drfs.splice( begin(m_drfs), m_drfs, iter );  //most recently used to the front
          return m_drfs.front().second;
        }
      }
      return nullptr;
    }

    void storeDrf( const string &key, shared_ptr<const DetectorPeakResponse> drf )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      m_drfs.remove_if( [&key]( const pair<string,shared_ptr<const DetectorPeakResponse>> &v ){
        return v.first == key;
      } );
      m_drfs.emplace_front( key, std::move(drf) );
      while( m_drfs.size() > ns_max_cached_drfs )
        m_drfs.pop_back();
    }

    void invalidatePath( const string &path )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      for( auto iter = begin(m_items); iter != end(m_items); )
      {
        bool under = false;
        for( const FileStamp &s : iter->second.stamps )
          under = under || SpecUtils::starts_with( s.path, path.c_str() );
        under = under || (iter->first.find(path) != string::npos);
        iter = under ? m_items.erase( iter ) : std::next( iter );
      }

      m_drfs.remove_if( [&path]( const pair<string,shared_ptr<const DetectorPeakResponse>> &v ){
        return v.first.find(path) != string::npos;
      } );
    }

    void invalidateRow( const long long id )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      m_rows.erase( id );
      const string key = "db:" + std::to_string(id);
      m_drfs.remove_if( [&key]( const pair<string,shared_ptr<const DetectorPeakResponse>> &v ){
        return v.first == key;
      } );
    }

    void clear()
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      m_items.clear();
      m_rows.clear();
      m_drfs.clear();
    }

    std::atomic<size_t> m_numComputed{ 0 };

  private:
    std::mutex m_mutex;
    map<string,ItemResult> m_items;
    map<long long,CachedDbRow> m_rows;
    list<pair<string,shared_ptr<const DetectorPeakResponse>>> m_drfs;
  };//class Cache


  shared_ptr<const DrfLibrary::Summary> summarize( const DetectorPeakResponse &drf )
  {
    Cache::instance().m_numComputed += 1;
    return DrfLibrary::Summary::fromDrf( drf );
  }


  bool is_table_extension( const string &lower_name )
  {
    if( lower_name == "efficiency.csv" )  // GADRAS; read with its Detector.dat
      return false;
    return SpecUtils::iends_with( lower_name, ".tsv" ) || SpecUtils::iends_with( lower_name, ".csv" );
  }


  bool is_other_extension( const string &lower_name )
  {
    return SpecUtils::iends_with( lower_name, ".ecc" )
           || SpecUtils::iends_with( lower_name, ".par" )
           || SpecUtils::iends_with( lower_name, ".outx" )
           || SpecUtils::iends_with( lower_name, ".detx" )
           || SpecUtils::iends_with( lower_name, ".xml" );
  }


  /** Whether a .xml file looks like an InterSpec DRF (cheap: looks at the start of the file). */
  bool looks_like_drf_xml( const string &path )
  {
#ifdef _WIN32
    ifstream input( SpecUtils::convert_from_utf8_to_utf16(path).c_str(), ios::in | ios::binary );
#else
    ifstream input( path.c_str(), ios::in | ios::binary );
#endif
    if( !input )
      return false;
    char buffer[4096] = { '\0' };
    input.read( buffer, sizeof(buffer) - 1 );
    const string header( buffer, static_cast<size_t>( input.gcount() ) );
    return header.find( "<DetectorPeakResponse" ) != string::npos;
  }//looks_like_drf_xml(...)


  /** A DETECTOR.txt (any case) in `dir` or in its `PARAMETERS` sub-folder; empty if neither. */
  string detector_txt_in( const string &dir )
  {
    const auto find_in = []( const string &folder ) -> string {
      for( const string &f : SpecUtils::ls_files_in_directory( folder, "detector.txt" ) )
      {
        if( lower_ascii( SpecUtils::filename( f ) ) == "detector.txt" )
          return f;
      }
      return string();
    };

    const string here = find_in( dir );
    if( !here.empty() )
      return here;

    for( const string &sub : SpecUtils::ls_directories_in_directory( dir ) )
    {
      if( lower_ascii( SpecUtils::filename( sub ) ) != "parameters" )
        continue;
      const string txt = find_in( sub );
      if( !txt.empty() )
        return txt;
    }

    return string();
  }//detector_txt_in(...)


  /** The DETECTOR.txt for a .par grid: in the grid's folder, or in that folder's or one of up to
   `max_up` parent folders' `PARAMETERS` sub-folder.  ISOCS keeps grids in `isocs\data\DCG\...` and
   the records for all of them in `isocs\data\PARAMETERS\DETECTOR.txt`.  `memo` holds per-folder
   answers for the rest of a walk.
   */
  string find_par_companion( const string &par_path, const int max_up, map<string,string> &memo )
  {
    string dir = SpecUtils::parent_path( par_path );
    for( int level = 0; (level <= max_up) && !dir.empty(); ++level )
    {
      map<string,string>::iterator pos = memo.find( dir );
      if( pos == end(memo) )
        pos = memo.emplace( dir, detector_txt_in( dir ) ).first;
      if( !pos->second.empty() )
        return pos->second;

      const string parent = SpecUtils::parent_path( dir );
      if( parent == dir )
        break;
      dir = parent;
    }//for( walk up from the grid )

    return string();
  }//find_par_companion(...)


  /** Walks one location, appending the work it holds.  Returns whether the location exists. */
  bool collect_work( const DrfLibrary::Location &loc, const size_t loc_index,
                     const std::atomic<bool> *cancel, vector<WorkItem> &work )
  {
    using DrfLibrary::EntryKind;

    map<string,string> detector_txts;  // folder -> its DETECTOR.txt (see find_par_companion)

    if( !loc.isFolder )
    {
      if( !SpecUtils::is_file( loc.path ) )
        return false;

      const string lname = lower_ascii( SpecUtils::filename( loc.path ) );
      WorkItem item;
      item.kind = is_table_extension( lname ) ? EntryKind::TableRow : EntryKind::ImportFile;
      item.key = (item.kind == EntryKind::TableRow ? "row:" : "file:") + loc.path;
      item.path = loc.path;
      item.stampPaths.push_back( loc.path );
      item.locationIndex = loc_index;

      if( SpecUtils::iends_with( lname, ".par" ) )
      {
        item.companion = find_par_companion( loc.path, 3, detector_txts );
        if( !item.companion.empty() )
          item.stampPaths.push_back( item.companion );
      }//if( a .par )

      work.push_back( item );
      return true;
    }//if( a single file )

    if( !SpecUtils::is_directory( loc.path ) )
      return false;

    size_t num_files_seen = 0;

    std::function<void(const string &, int)> walk;
    walk = [&]( const string &dir, const int depth ){
      if( (depth > ns_max_folder_depth) || (num_files_seen > ns_max_files_per_location) )
        return;
      if( cancel && cancel->load() )
        return;

      const vector<string> files = SpecUtils::ls_files_in_directory( dir );
      num_files_seen += files.size();

      // GADRAS: a folder with a Detector.dat is one detector (Efficiency.csv optional).
      string detector_dat, efficiency_csv;
      for( const string &f : files )
      {
        const string lname = lower_ascii( SpecUtils::filename( f ) );
        if( lname == "detector.dat" )
          detector_dat = f;
        else if( lname == "efficiency.csv" )
          efficiency_csv = f;
      }

      if( !detector_dat.empty() && (loc.formats & DrfLibrary::FormatGadras) )
      {
        WorkItem item;
        item.kind = EntryKind::GadrasFolder;
        item.key = "gadras:" + dir;
        item.path = dir;
        item.stampPaths.push_back( detector_dat );
        if( !efficiency_csv.empty() )
          item.stampPaths.push_back( efficiency_csv );
        item.locationIndex = loc_index;
        work.push_back( item );
      }//if( a GADRAS detector folder )

      for( const string &f : files )
      {
        const string lname = lower_ascii( SpecUtils::filename( f ) );
        if( !lname.empty() && (lname[0] == '.') )
          continue;

        if( is_table_extension( lname ) )
        {
          if( !(loc.formats & DrfLibrary::FormatTable) )
            continue;

          WorkItem item;
          item.kind = EntryKind::TableRow;
          item.key = "row:" + f;
          item.path = f;
          item.stampPaths.push_back( f );
          item.locationIndex = loc_index;
          work.push_back( item );
          continue;
        }//if( a table )

        if( !(loc.formats & DrfLibrary::FormatOther) || !is_other_extension( lname ) )
          continue;

        if( SpecUtils::iends_with( lname, ".xml" ) && !looks_like_drf_xml( f ) )
          continue;

        WorkItem item;
        item.kind = EntryKind::ImportFile;
        item.key = "file:" + f;
        item.path = f;
        item.stampPaths.push_back( f );
        if( SpecUtils::iends_with( lname, ".par" ) )
        {
          // Up to the location, and two above it - in case the user located `isocs\data\DCG`.
          item.companion = find_par_companion( f, depth + 2, detector_txts );
          if( !item.companion.empty() )
            item.stampPaths.push_back( item.companion );
        }
        item.locationIndex = loc_index;
        work.push_back( item );
      }//for( const string &f : files )

      for( const string &sub : SpecUtils::ls_directories_in_directory( dir ) )
      {
        const string lname = SpecUtils::filename( sub );
        if( !lname.empty() && (lname[0] == '.') )
          continue;
        walk( sub, depth + 1 );
      }
    };//walk

    walk( loc.path, 0 );

    return true;
  }//collect_work(...)


  /** Reads the DRFs of one work item.  Never throws. */
  vector<shared_ptr<DetectorPeakResponse>> read_item( const WorkItem &item,
                                                      const DrfLibrary::Location &loc,
                                                      size_t &num_unreadable, string &error )
  {
    using DrfLibrary::EntryKind;

    vector<shared_ptr<DetectorPeakResponse>> answer;

    try
    {
      switch( item.kind )
      {
        case EntryKind::GadrasFolder:
        {
          auto drf = make_shared<DetectorPeakResponse>( SpecUtils::filename( item.path ), "" );
          drf->fromGadrasDirectory( item.path, true /* allow a missing Efficiency.csv */ );

          // Matches what choosing a detector automatically for a spectrum has always recorded.
          if( drf->isValid() )
            drf->setDrfSource( loc.builtIn ? DetectorPeakResponse::DrfSource::DefaultGadrasDrf
                                           : DetectorPeakResponse::DrfSource::UserAddedGadrasDrf );
          answer.push_back( drf );
          break;
        }//case EntryKind::GadrasFolder:

        case EntryKind::TableRow:
        {
#ifdef _WIN32
          ifstream input( SpecUtils::convert_from_utf8_to_utf16(item.path).c_str(), ios::in | ios::binary );
#else
          ifstream input( item.path.c_str(), ios::in | ios::binary );
#endif
          if( !input )
            throw runtime_error( "Could not open file." );

          vector<string> credits;
          vector<shared_ptr<DetectorPeakResponse>> drfs;
          DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );

          for( const shared_ptr<DetectorPeakResponse> &drf : drfs )
          {
            if( !drf )
              continue;
            if( loc.builtIn )
              drf->setDrfSource( DetectorPeakResponse::DrfSource::DefaultRelativeEfficiencyDrf );
            answer.push_back( drf );
          }

          // Not a one-detector-per-row table: maybe another complete-DRF CSV (Make DRF, GammaQuant).
          if( answer.empty() )
          {
            const shared_ptr<const string> data = DrfImport::readFile( item.path );
            const shared_ptr<const DrfImport::ParsedFile> parsed
                               = DrfImport::parseFile( SpecUtils::filename(item.path), data, false );
            if( parsed && DrfImport::isCompleteDrfKind( parsed->kind ) )
            {
              for( const shared_ptr<const DetectorPeakResponse> &drf : parsed->drfs )
                if( drf && drf->isValid() )
                  answer.push_back( make_shared<DetectorPeakResponse>( *drf ) );
            }
          }//if( answer.empty() )

          if( answer.empty() )
            num_unreadable += 1;
          break;
        }//case EntryKind::TableRow:

        case EntryKind::ImportFile:
        {
          const shared_ptr<const string> data = DrfImport::readFile( item.path );
          const shared_ptr<const DrfImport::ParsedFile> parsed
                             = DrfImport::parseFile( SpecUtils::filename(item.path), data, false );
          if( !parsed )
            throw runtime_error( "Not a detector file." );

          if( DrfImport::isCompleteDrfKind( parsed->kind ) )
          {
            for( const shared_ptr<const DetectorPeakResponse> &drf : parsed->drfs )
              if( drf )
                answer.push_back( make_shared<DetectorPeakResponse>( *drf ) );
            break;
          }

          shared_ptr<const DrfImport::ParsedFile> companion;
          if( !item.companion.empty() )
            companion = DrfImport::parseFile( SpecUtils::filename(item.companion),
                                              DrfImport::readFile( item.companion ), false );

          const DrfImport::Source src = DrfImport::makeSource( parsed, companion, 0 );
          if( !src.recordNames.empty() )
          {
            // No DETECTOR.txt record names this grid; guessing one would show the wrong crystal.
            num_unreadable += 1;
            error = "No record in '" + SpecUtils::filename( item.companion ) + "' is for '"
                    + SpecUtils::filename( item.path ) + "'.";
            break;
          }

          if( src.status == DrfImport::Status::NeedsCharacterization )
          {
            if( src.base )
              answer.push_back( make_shared<DetectorPeakResponse>( *src.base ) );
            break;
          }

          if( src.status != DrfImport::Status::Ready )
          {
            // Needs a companion, a diameter, a distance, ...: the user has to import it.
            num_unreadable += 1;
            break;
          }

          DrfImport::Options options;
          options.interpretation = src.defaultInterpretation;
          const DrfImport::Result result = DrfImport::build( src, options );
          if( (result.status == DrfImport::Status::Ready) && result.drf )
            answer.push_back( result.drf );
          else
            num_unreadable += 1;
          break;
        }//case EntryKind::ImportFile:

        case EntryKind::DbRow:
        case EntryKind::InMemory:
          assert( 0 );
          break;
      }//switch( item.kind )
    }catch( std::exception &e )
    {
      num_unreadable += 1;
      error = e.what();
    }//try / catch

    return answer;
  }//read_item(...)


  /** Stable from the moment a file is listed: a file's first detector keeps the id it had while
   still being read (when how many it holds is not known yet). */
  string entry_id_for( const WorkItem &item, const size_t index, const size_t )
  {
    if( item.kind == DrfLibrary::EntryKind::TableRow )
      return item.key + "|" + std::to_string(index);
    return (index == 0) ? item.key : (item.key + "|" + std::to_string(index));
  }


  DrfLibrary::Entry make_located_entry( const WorkItem &item, const DrfLibrary::Location &loc,
                                        const size_t index, const size_t num,
                                        shared_ptr<const DrfLibrary::Summary> summary )
  {
    DrfLibrary::Entry e;
    e.id = entry_id_for( item, index, num );
    e.kind = item.kind;
    e.locationId = loc.id;
    e.path = item.path;
    e.companionPath = item.companion;
    e.record = index;
    e.hidden = loc.hidden;
    e.summary = summary;
    if( summary )
    {
      e.knownHashes.push_back( summary->hash );
      e.knownHashes.insert( end(e.knownHashes), begin(summary->sourceAliasHashes), end(summary->sourceAliasHashes) );
    }

    // Where a detector is found says what it is: the kind of file, rather than the source it is
    //  tagged with (a row of a table someone shared may be tagged "made by a user").  Only an
    //  InterSpec export (.xml/.csv) goes by its tag, since it can be any kind.
    const string lname = lower_ascii( item.path );
    DrfLibrary::Origin origin = summary ? DrfLibrary::originForSource( summary->source )
                                        : DrfLibrary::Origin::UserMade;
    if( item.kind == DrfLibrary::EntryKind::GadrasFolder )
      origin = DrfLibrary::Origin::Gadras;
    else if( item.kind == DrfLibrary::EntryKind::TableRow )
      origin = DrfLibrary::Origin::RelEffTable;
    else if( SpecUtils::iends_with( lname, ".ecc" ) || SpecUtils::iends_with( lname, ".par" ) )
      origin = DrfLibrary::Origin::Isocs;
    else if( SpecUtils::iends_with( lname, ".outx" ) || SpecUtils::iends_with( lname, ".detx" ) )
      origin = DrfLibrary::Origin::Angle;
    e.origin = loc.builtIn ? DrfLibrary::Origin::BuiltIn : origin;
    return e;
  }//make_located_entry(...)


  /** Origin of a database row not found in any location. */
  DrfLibrary::Origin db_row_origin( const DrfLibrary::Summary &summary )
  {
    using DrfSource = DetectorPeakResponse::DrfSource;

    if( summary.parentHash != 0 )
      return DrfLibrary::Origin::UserMade;

    switch( summary.source )
    {
      // A shipped detector no longer shipped, or a folder/table no longer listed: keep its chip.
      case DrfSource::DefaultGadrasDrf:
      case DrfSource::DefaultRelativeEfficiencyDrf:
        return DrfLibrary::Origin::BuiltIn;
      case DrfSource::UserAddedGadrasDrf:
        return DrfLibrary::Origin::Gadras;
      case DrfSource::UserAddedRelativeEfficiencyDrf:
        return DrfLibrary::Origin::RelEffTable;

      // Everything imported once, entered, or made in InterSpec is "Made by you".
      default:
        break;
    }//switch( summary.source )

    return DrfLibrary::Origin::UserMade;
  }//db_row_origin(...)


  /** The Library's sort order for names: case-insensitive, then exact. */
  bool name_less( const string &lhs, const string &rhs )
  {
    const string l = lower_ascii( lhs ), r = lower_ascii( rhs );
    return (l != r) ? (l < r) : (lhs < rhs);
  }
}//namespace


namespace DrfLibrary
{

Origin originForSource( const DetectorPeakResponse::DrfSource source )
{
  using DrfSource = DetectorPeakResponse::DrfSource;

  switch( source )
  {
    case DrfSource::DefaultGadrasDrf:
    case DrfSource::DefaultRelativeEfficiencyDrf:
      return Origin::BuiltIn;

    case DrfSource::UserAddedGadrasDrf:
    case DrfSource::UserImportedGadrasDrf:
    case DrfSource::GadrasDetectorDatOnly:
      return Origin::Gadras;

    case DrfSource::IsocsEcc:
    case DrfSource::CharacterizationParFile:
      return Origin::Isocs;

    case DrfSource::AngleOutx:
      return Origin::Angle;

    case DrfSource::UserAddedRelativeEfficiencyDrf:
      return Origin::RelEffTable;

    case DrfSource::FromSpectrumFileDrf:
      return Origin::ThisSpectrum;

    case DrfSource::UnknownDrfSource:
    case DrfSource::UserImportedIntrisicEfficiencyDrf:
    case DrfSource::UserSpecifiedFormulaDrf:
    case DrfSource::UserCreatedDrf:
    case DrfSource::UserImportedEfficiencyCsvDrf:
      return Origin::UserMade;
  }//switch( source )

  return Origin::UserMade;
}//Origin originForSource( const DetectorPeakResponse::DrfSource source )


std::shared_ptr<const Summary> Summary::fromDrf( const DetectorPeakResponse &drf )
{
  auto s = make_shared<Summary>();

  s->name = drf.name();
  s->description = drf.description();
  s->source = drf.drfSource();
  s->geometryType = drf.geometryType();
  s->hash = drf.hashValue();
  s->parentHash = drf.parentHashValue();

  {
    using DrfSource = DetectorPeakResponse::DrfSource;
    DrfSource alias = s->source;
    switch( s->source )
    {
      case DrfSource::DefaultGadrasDrf:             alias = DrfSource::UserAddedGadrasDrf;             break;
      case DrfSource::UserAddedGadrasDrf:           alias = DrfSource::DefaultGadrasDrf;               break;
      case DrfSource::DefaultRelativeEfficiencyDrf: alias = DrfSource::UserAddedRelativeEfficiencyDrf; break;
      case DrfSource::UserAddedRelativeEfficiencyDrf: alias = DrfSource::DefaultRelativeEfficiencyDrf; break;
      default: break;
    }
    if( alias != s->source )
      s->sourceAliasHashes.push_back( drf.hashValueWithSource( alias ) );
  }
  s->appearance = drf.appearance();
  s->lowerEnergy = drf.lowerEnergy();
  s->upperEnergy = drf.upperEnergy();

  s->hasEfficiency = drf.isValid();

  const shared_ptr<const ceelo::DetectorResponse> response = drf.ceeloResponse();
  s->hasResponse = !!response;
  s->hasGeometry = (drf.geometry() || response);

  if( !s->hasEfficiency )
    s->range = drf.storedGeometry() ? Range::NeedsCharacterization : Range::FarFieldOnly;
  else if( drf.isFixedGeometry() )
    s->range = Range::FixedGeometry;
  else if( response )
    s->range = Range::NearFieldOk;
  else
    s->range = Range::FarFieldOnly;

  // The crystal: the recorded geometry if there is one (stored, even if switched off, still says what
  //  the crystal is), else the text, else what the detector's own data implies.
  const ceelo::GeometryDescriptor *geom = nullptr;
  const shared_ptr<const ceelo::GeometryDescriptor> stored = drf.storedGeometry();
  if( stored )
    geom = stored.get();
  else if( response )
    geom = &response->descriptor;

  if( geom )
  {
    const int index = geom->crystal_material_index;
    if( (index >= 0) && (static_cast<size_t>(index) < geom->materials.size()) )
      s->crystal.material = CeeLoUtils::crystalNameFromText( geom->materials[index].name );

    const vector<double> &dims = geom->dimensions_cm;
    if( (geom->shape == ceelo::DetectorShape::Box) && (dims.size() >= 3) )
    {
      s->crystal.fromGeometry = true;
      s->crystal.isBox = true;
      s->crystal.diameter_cm = 2.0*dims[0];
      s->crystal.height_cm = 2.0*dims[1];
      s->crystal.length_cm = dims[2];
    }else if( (geom->shape != ceelo::DetectorShape::Box) && (dims.size() >= 2) )
    {
      s->crystal.fromGeometry = true;
      s->crystal.diameter_cm = 2.0*dims[0];
      s->crystal.length_cm = dims[1];
    }
  }//if( geom )

  if( s->crystal.material.empty() )
    s->crystal.material = CeeLoUtils::crystalNameFromText( drf.description() + " " + drf.name() );

  try
  {
    if( s->crystal.material.empty() )
      s->crystal.material = CeeLoUtils::crystalNameFromDrfEvidence( drf );
  }catch( std::exception & )
  {
  }

  if( !drf.isFixedGeometry() )
    s->crystal.faceDiameter_cm = drf.detectorDiameter() / PhysicalUnits::cm;

  // "HPGe 40%", "Canberra Falcon5000_LANL_025cm (18%)": the list need not repeat the number.
  static const std::regex percent_re( "[0-9]+(\\.[0-9]+)?\\s*%" );
  s->nameStatesRelEff = std::regex_search( drf.name(), percent_re );

  if( !s->hasEfficiency )
    return s;

  try
  {
    s->relEffPercent = drf.relativeEfficiencyIeee325();
  }catch( std::exception & )
  {
    s->relEffPercent = -1.0;
  }

  s->efficiencyForm = drf.efficiencyFcnType();
  s->importedGrid = drf.hasImportedGrid();

  const shared_ptr<const DetectorEfficiencyCurve> curve = drf.efficiencyCurve();
  switch( s->efficiencyForm )
  {
    case DetectorPeakResponse::kExpOfLogPowerSeries:
      s->efficiencyParamCount = drf.efficiencyExpOfLogsCoeffs().size();
      break;
    case DetectorPeakResponse::kEnergyEfficiencyPairs:
      s->efficiencyParamCount = curve ? curve->energyEfficiencies().size() : 0;
      break;
    case DetectorPeakResponse::kFunctialEfficienyForm:
    case DetectorPeakResponse::kNumEfficiencyFnctForms:
      break;
  }//switch( s->efficiencyForm )

  const shared_ptr<const MeasuredDrfPoints> points = drf.measuredPoints();
  s->measuredPointCount = (points && !points->empty()) ? points->points().size() : 0;

  try
  {
    s->hasFwhm = drf.hasResolutionInfo();
    if( s->hasFwhm )
      s->fwhmAt662 = drf.peakResolutionFWHM( 661.657f );
  }catch( std::exception & )
  {
    s->hasFwhm = false;
    s->fwhmAt662 = -1.0f;
  }

  try
  {
    const DrfModifyCalc::UncertSummary unc = DrfModifyCalc::uncertSummary( drf );
    s->hasUncertainty = (unc.valid && (unc.total > 0.0));
    s->uncertaintyIsAssumed = s->hasUncertainty && unc.dataIsAssumed;
  }catch( std::exception & )
  {
  }

  try
  {
    s->hasTotalEfficiency = drf.hasAnyTotalEfficiencyInfo();
  }catch( std::exception & )
  {
  }

  return s;
}//Summary::fromDrf(...)


std::vector<Location> standardLocations()
{
  vector<Location> locations;

  // Shipped with InterSpec
  try
  {
    const string data_dir = InterSpec::staticDataDirectory();

    Location generic;
    generic.path = SpecUtils::append_path( data_dir, "GenericGadrasDetectors" );
    generic.id = generic.path;
    generic.role = LocationRole::GenericGadras;
    generic.isAuto = generic.builtIn = true;
    generic.formats = FormatGadras;
    locations.push_back( generic );

    Location table;
    table.path = SpecUtils::append_path( data_dir, "common_drfs.tsv" );
    table.id = table.path;
    table.role = LocationRole::BuiltInTable;
    table.isFolder = false;
    table.isAuto = table.builtIn = true;
    locations.push_back( table );

#if( BUILD_FOR_WEB_DEPLOYMENT )
    // Web deployments may ship an extra set of GADRAS detectors.
    const string ouo = SpecUtils::append_path( data_dir, "OUO_GadrasDetectors" );
    if( SpecUtils::is_directory( ouo ) )
    {
      Location extra = generic;
      extra.path = extra.id = ouo;
      extra.role = LocationRole::ExtraGadras;
      locations.push_back( extra );
    }
#endif
  }catch( std::exception &e )
  {
    cerr << "DrfLibrary::standardLocations(): no static data directory: " << e.what() << endl;
  }//try / catch

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  // The user-data directory: `drfs/` holds files saved from the app (any format), and - as before -
  //  any GADRAS detector folders or relative-efficiency tables anywhere in it are listed too.
  try
  {
    const string user_dir = InterSpec::writableDataDirectory();

    Location drfs;
    drfs.path = SpecUtils::append_path( user_dir, "drfs" );
    drfs.id = drfs.path;
    drfs.role = LocationRole::KeptImports;
    drfs.isAuto = true;
    locations.push_back( drfs );

    Location user;
    user.path = user.id = user_dir;
    user.role = LocationRole::UserData;
    user.isAuto = true;
    user.formats = FormatGadras | FormatTable;
    locations.push_back( user );
  }catch( std::exception & )
  {
    // No writable directory set (e.g., some unit tests).
  }
#endif

  // Where the tools that make detector files keep them by default (they are Windows-only):
  //  - GADRAS: a folder per detector (Detector.dat + Efficiency.csv); the "GadrasDRFPath" default.
  //  - ISOCS (Genie 2000): .par grids in `DCG\` sub-folders, with every detector's record in
  //    `PARAMETERS\DETECTOR.txt` (see find_par_companion); .ecc files anywhere below are listed too.
  //  - ANGLE 5: results (.outx) and exported detectors (.detx) go in sub-folders of Documents\ANGLE
  //    (ANGLE 5 User's Manual, 931077 Rev D, 2022).  Its own detector list,
  //    %LOCALAPPDATA%\ANGLE\Detectors.datx, is not a format InterSpec reads.
#if( defined(_WIN32) && !BUILD_FOR_WEB_DEPLOYMENT )
  {
    vector<tuple<string,string,uint32_t,LocationRole>> candidates{
      { "GADRAS", "C:\\GADRAS\\Detector", FormatGadras, LocationRole::GadrasInstall },
      { "ISOCS", "C:\\GENIE2K\\isocs\\data", FormatOther, LocationRole::IsocsInstall }
    };

    const string documents = documents_dir();
    if( !documents.empty() )
      candidates.emplace_back( "ANGLE", SpecUtils::append_path( documents, "ANGLE" ), FormatOther,
                               LocationRole::AngleInstall );

    for( const tuple<string,string,uint32_t,LocationRole> &c : candidates )
    {
      Location loc;
      loc.id = std::get<0>( c );
      loc.path = std::get<1>( c );
      loc.formats = std::get<2>( c );
      loc.role = std::get<3>( c );
      loc.isAuto = true;
      locations.push_back( loc );
    }
  }
#endif

  return locations;
}//std::vector<Location> standardLocations()


bool LocationSettings::operator==( const LocationSettings &rhs ) const
{
  if( (hiddenStandard != rhs.hiddenStandard) || (standardPaths != rhs.standardPaths)
      || (userLocations.size() != rhs.userLocations.size()) )
    return false;
  for( size_t i = 0; i < userLocations.size(); ++i )
  {
    const Location &a = userLocations[i], &b = rhs.userLocations[i];
    if( (a.path != b.path) || (a.hidden != b.hidden) || (a.formats != b.formats) )
      return false;
  }
  return true;
}//LocationSettings::operator==


LocationSettings parseLocationSettings( const std::string &pref )
{
  LocationSettings settings;

  vector<string> entries;
  SpecUtils::split( entries, pref, ";\r\n" );
  for( string entry : entries )
  {
    SpecUtils::trim( entry );
    if( !entry.empty() && (entry[0] == '#') )  //the version marker #serializeLocationSettings writes
      continue;
    bool hidden = false;
    if( !entry.empty() && (entry[0] == '!') )
    {
      hidden = true;
      entry = SpecUtils::trim_copy( entry.substr( 1 ) );
    }
    if( entry.empty() )
      continue;

    if( entry[0] == '@' )
    {
      // A standard location: "@id" or "@id=path"
      const size_t eq = entry.find( '=' );
      const string id = SpecUtils::trim_copy( entry.substr( 1, (eq == string::npos) ? string::npos : (eq - 1) ) );
      if( id.empty() )
        continue;
      if( hidden )
        settings.hiddenStandard.insert( id );
      if( eq != string::npos )
      {
        const string path = SpecUtils::trim_copy( entry.substr( eq + 1 ) );
        if( !path.empty() )
          settings.standardPaths[id] = path;
      }
      continue;
    }//if( a standard location )

    Location loc;
    loc.formats = FormatAll;
    if( SpecUtils::istarts_with( entry, "gadras:" ) && (entry.size() > 7) )
    {
      loc.formats = FormatGadras;
      entry = SpecUtils::trim_copy( entry.substr( 7 ) );
    }
    loc.id = loc.path = entry;
    loc.hidden = hidden;

    bool duplicate = false;
    for( const Location &prev : settings.userLocations )
      duplicate = duplicate || (prev.path == loc.path);
    if( !duplicate && !loc.path.empty() )
      settings.userLocations.push_back( loc );
  }//for( string entry : entries )

  return settings;
}//parseLocationSettings(...)


std::string serializeLocationSettings( const LocationSettings &settings )
{
  vector<string> entries;
  for( const Location &loc : settings.userLocations )
  {
    if( loc.path.empty() )
      continue;
    entries.push_back( (loc.hidden ? "!" : "") + string( (loc.formats == FormatGadras) ? "gadras:" : "" ) + loc.path );
  }

  set<string> ids = settings.hiddenStandard;
  for( const pair<const string,string> &p : settings.standardPaths )
    ids.insert( p.first );
  for( const string &id : ids )
  {
    const auto path = settings.standardPaths.find( id );
    entries.push_back( (settings.hiddenStandard.count( id ) ? "!@" : "@") + id
                       + ((path != end(settings.standardPaths)) ? ("=" + path->second) : string()) );
  }

  // Always non-empty, so a saved "nothing added" differs from never having been saved.
  string answer = "#1";
  for( const string &e : entries )
    answer += ";" + e;
  return answer;
}//serializeLocationSettings(...)


std::vector<Location> applyLocationSettings( std::vector<Location> standard,
                                             const LocationSettings &settings )
{
  set<string> have;
  for( Location &loc : standard )
  {
    loc.hidden = settings.hiddenStandard.count( loc.id ) > 0;
    const auto path = settings.standardPaths.find( loc.id );
    if( path != end(settings.standardPaths) )
      loc.path = path->second;
    have.insert( SpecUtils::lexically_normalize_path( loc.path ) );
  }

  for( Location loc : settings.userLocations )
  {
    const string normalized = SpecUtils::lexically_normalize_path( loc.path );
    if( have.count( normalized ) )
      continue;
    have.insert( normalized );
    loc.role = LocationRole::User;
    loc.isAuto = loc.builtIn = false;
    loc.isFolder = !SpecUtils::is_file( loc.path );  //a missing one is shown as a folder
    standard.push_back( loc );
  }

  return standard;
}//applyLocationSettings(...)


Snapshot scan( const ScanRequest &request )
{
  Snapshot snapshot;
  Cache &cache = Cache::instance();
  const std::atomic<bool> *cancel = request.cancel.get();
  auto cancelled = [cancel]() -> bool { return cancel && cancel->load(); };

  // 1) Walk the locations.
  vector<WorkItem> work;
  for( size_t i = 0; i < request.locations.size(); ++i )
  {
    LocationStatus status;
    status.id = request.locations[i].id;
    status.found = collect_work( request.locations[i], i, cancel, work );
    snapshot.locations.push_back( status );
    if( cancelled() )
      return snapshot;
  }

  // The same file can be reached from two locations (one inside another): first one wins.
  {
    set<string> seen;
    vector<WorkItem> unique_work;
    for( WorkItem &item : work )
      if( seen.insert( item.key ).second )
        unique_work.push_back( std::move(item) );
    work.swap( unique_work );
  }

  // 2) Reuse what has not changed.
  vector<ItemResult> results( work.size() );
  vector<bool> have_result( work.size(), false );
  vector<size_t> to_read;
  for( size_t i = 0; i < work.size(); ++i )
  {
    vector<FileStamp> stamps;
    for( const string &p : work[i].stampPaths )
      stamps.push_back( stamp_file( p ) );

    if( cache.findItem( work[i].key, stamps, results[i] ) )
    {
      have_result[i] = true;
    }else
    {
      results[i].stamps = stamps;
      to_read.push_back( i );
    }
  }//for( size_t i = 0; i < work.size(); ++i )

  // The database rows: a column-only query, so no row's DrfExtra XML (which can hold a large CeeLo
  //  response) is parsed unless that row changed.
  struct RowInfo { long long id; int64_t hash, parent, last_used; long long extra_len; };
  vector<RowInfo> rows;
  vector<shared_ptr<const Summary>> row_summaries;
  vector<size_t> rows_to_read;

  if( request.sql && (request.dbUserId >= 0) )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *request.sql );

      typedef std::tuple<long long, long long, long long, long long, long long> RowTuple;
      const Wt::Dbo::collection<RowTuple> db_rows
          = request.sql->session()->query<RowTuple>( "SELECT id, Hash, ParentHash, m_lastUsedUtc, "
                                                     "COALESCE(LENGTH(m_drfExtra),0) FROM DetectorPeakResponse" )
                                        .where( "InterSpecUser_id = ?" )
                                        .bind( request.dbUserId )
                                        .resultList();
      for( const RowTuple &t : db_rows )
        rows.push_back( RowInfo{ std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t) } );

      transaction.commit();
    }catch( std::exception &e )
    {
      cerr << "DrfLibrary::scan(): error querying detector rows: " << e.what() << endl;
      rows.clear();
    }//try / catch

    row_summaries.resize( rows.size() );
    for( size_t i = 0; i < rows.size(); ++i )
    {
      row_summaries[i] = cache.findRow( rows[i].id, rows[i].hash, rows[i].extra_len );
      if( !row_summaries[i] )
        rows_to_read.push_back( i );
    }
  }//if( request.sql && (request.dbUserId >= 0) )

  if( cancelled() )
    return snapshot;

  // Builds entries from whatever is known so far.
  auto build_entries = [&]( const bool final_pass ) -> vector<Entry> {
    vector<Entry> located;
    for( size_t i = 0; i < work.size(); ++i )
    {
      const WorkItem &item = work[i];
      const Location &loc = request.locations[item.locationIndex];
      if( !have_result[i] )
      {
        if( !final_pass )
        {
          Entry e = make_located_entry( item, loc, 0, 1, nullptr );
          e.pending = true;
          located.push_back( e );
        }
        continue;
      }

      const vector<shared_ptr<const Summary>> &sums = results[i].summaries;
      for( size_t j = 0; j < sums.size(); ++j )
        located.push_back( make_located_entry( item, loc, j, sums.size(), sums[j] ) );
    }//for( size_t i = 0; i < work.size(); ++i )

    vector<Entry> db_entries;
    for( size_t i = 0; i < rows.size(); ++i )
    {
      if( !row_summaries[i] )
        continue;
      Entry e;
      e.id = "db:" + std::to_string( rows[i].id );
      e.kind = EntryKind::DbRow;
      e.dbId = rows[i].id;
      e.lastUsedUtc = rows[i].last_used;
      e.summary = row_summaries[i];
      e.knownHashes.push_back( row_summaries[i]->hash );
      db_entries.push_back( e );
    }

    return mergeEntries( std::move(located), std::move(db_entries), request.suggestions );
  };//build_entries

  if( request.progress && (!to_read.empty() || !rows_to_read.empty()) )
  {
    Snapshot partial = snapshot;
    partial.entries = build_entries( false );
    try
    {
      request.progress( partial );
    }catch( std::exception &e )
    {
      cerr << "DrfLibrary::scan(): progress callback threw: " << e.what() << endl;
    }
  }//if( report progress )

  // 3) Read what changed - files in parallel, database rows on this thread (one DB session).
  if( !to_read.empty() )
  {
    SpecUtilsAsync::ThreadPool pool;
    for( const size_t index : to_read )
    {
      pool.post( [index, &work, &results, &request, cancel](){
        if( cancel && cancel->load() )
          return;

        const WorkItem &item = work[index];
        ItemResult &result = results[index];
        const vector<shared_ptr<DetectorPeakResponse>> drfs
                 = read_item( item, request.locations[item.locationIndex], result.numUnreadable, result.error );
        for( const shared_ptr<DetectorPeakResponse> &drf : drfs )
        {
          try
          {
            result.summaries.push_back( summarize( *drf ) );
          }catch( std::exception &e )
          {
            result.numUnreadable += 1;
            result.error = e.what();
          }
        }//for( const auto &drf : drfs )
      } );
    }//for( const size_t index : to_read )
    pool.join();

    if( cancelled() )
      return snapshot;

    for( const size_t index : to_read )
    {
      have_result[index] = true;
      cache.storeItem( work[index].key, results[index] );
    }
  }//if( !to_read.empty() )

  if( !rows_to_read.empty() )
  {
    try
    {
      for( const size_t index : rows_to_read )
      {
        if( cancelled() )
          break;

        const RowInfo &info = rows[index];
        try
        {
          // A transaction per row, holding the (shared) DB session only to copy the row out - the
          //  summary (which can build a CeeLo response) is made without it.
          unique_ptr<DetectorPeakResponse> drf;
          {
            DataBaseUtils::DbTransaction transaction( *request.sql );
            const Wt::Dbo::ptr<DetectorPeakResponse> row
                      = request.sql->session()->find<DetectorPeakResponse>()
                                                .where( "id = ?" ).bind( info.id ).resultValue();
            if( row )
              drf = make_unique<DetectorPeakResponse>( *row );
            transaction.commit();
          }
          if( !drf )
            continue;

          row_summaries[index] = summarize( *drf );

          CachedDbRow cached;
          cached.hash = info.hash;
          cached.extraLength = info.extra_len;
          cached.summary = row_summaries[index];
          cache.storeRow( info.id, cached );
        }catch( std::exception &e )
        {
          cerr << "DrfLibrary::scan(): error reading detector row " << info.id << ": " << e.what() << endl;
        }
      }//for( const size_t index : rows_to_read )
    }catch( std::exception &e )
    {
      cerr << "DrfLibrary::scan(): error reading detector rows: " << e.what() << endl;
    }//try / catch
  }//if( !rows_to_read.empty() )

  if( cancelled() )
    return snapshot;

  // 4) Merge, and count per location.
  snapshot.entries = build_entries( true );

  for( size_t i = 0; i < work.size(); ++i )
  {
    LocationStatus &status = snapshot.locations[work[i].locationIndex];
    status.numDetectors += results[i].summaries.size();
    status.numUnreadable += results[i].numUnreadable;
  }

  return snapshot;
}//Snapshot scan( const ScanRequest &request )


std::vector<Entry> mergeEntries( std::vector<Entry> located,
                                 std::vector<Entry> dbRows,
                                 const std::vector<Suggestion> &suggestions )
{
  using DrfSource = DetectorPeakResponse::DrfSource;

  vector<Entry> entries = std::move( located );

  // Which kind of location-found detector a source says, for matching stale rows by name.
  auto family = []( const DrfSource source ) -> int {
    switch( source )
    {
      case DrfSource::DefaultGadrasDrf:
      case DrfSource::UserAddedGadrasDrf:
        return 1;
      case DrfSource::DefaultRelativeEfficiencyDrf:
      case DrfSource::UserAddedRelativeEfficiencyDrf:
        return 2;
      default:
        break;
    }
    return 0;
  };//family

  // Hash -> entry, and (family, name) -> entry; a visible entry is preferred when a detector is in
  //  two locations.
  map<uint64_t,size_t> by_hash;
  map<pair<int,string>,size_t> by_name;
  auto index_entry = [&]( const size_t i ){
    const Entry &e = entries[i];
    for( const uint64_t h : e.knownHashes )
    {
      const auto pos = by_hash.find( h );
      if( (pos == end(by_hash)) || (entries[pos->second].hidden && !e.hidden) )
        by_hash[h] = i;
    }

    const int fam = e.summary ? family( e.summary->source ) : 0;
    if( fam && (e.kind != EntryKind::DbRow) && (e.kind != EntryKind::InMemory) )
    {
      const pair<int,string> key( fam, lower_ascii( e.summary->name ) );
      const auto pos = by_name.find( key );
      if( (pos == end(by_name)) || (entries[pos->second].hidden && !e.hidden) )
        by_name[key] = i;
    }
  };//index_entry

  for( size_t i = 0; i < entries.size(); ++i )
    index_entry( i );

  // The entry a DRF (summary) stands for, or -1.
  auto find_match = [&]( const Summary &summary ) -> int {
    const auto by_h = by_hash.find( summary.hash );
    if( by_h != end(by_hash) )
      return static_cast<int>( by_h->second );

    const int fam = family( summary.source );
    if( fam && !summary.parentHash )
    {
      const auto by_n = by_name.find( pair<int,string>( fam, lower_ascii( summary.name ) ) );
      if( by_n != end(by_name) )
        return static_cast<int>( by_n->second );
    }

    return -1;
  };//find_match

  // Newest rows first, so an entry takes the newest matching row's id.
  std::stable_sort( begin(dbRows), end(dbRows), []( const Entry &lhs, const Entry &rhs ){
    return lhs.lastUsedUtc > rhs.lastUsedUtc;
  } );

  for( Entry &row : dbRows )
  {
    if( !row.summary )
      continue;

    const int match = find_match( *row.summary );
    if( match >= 0 )
    {
      Entry &e = entries[static_cast<size_t>(match)];
      if( e.dbId < 0 )
        e.dbId = row.dbId;
      e.lastUsedUtc = std::max( e.lastUsedUtc, row.lastUsedUtc );
      if( std::find( begin(e.knownHashes), end(e.knownHashes), row.summary->hash ) == end(e.knownHashes) )
      {
        e.knownHashes.push_back( row.summary->hash );
        by_hash.emplace( row.summary->hash, static_cast<size_t>(match) );
      }
      continue;
    }//if( already listed )

    row.origin = db_row_origin( *row.summary );
    entries.push_back( row );
    index_entry( entries.size() - 1 );
  }//for( Entry &row : dbRows )

  // A suggestion (the detector in use, a preferred one, the spectrum's own) marks the visible entry
  //  that IS that detector - by its own hash, not one of a database row folded into it (a different
  //  detector of the same name, or an older version) - else it is listed on its own, so the detector
  //  is always there to see and pick, even when its location is hidden.
  for( const Suggestion &suggestion : suggestions )
  {
    if( !suggestion.drf )
      continue;

    const shared_ptr<const Summary> summary = summarize( *suggestion.drf );
    int match = -1;
    for( size_t i = 0; (match < 0) && (i < entries.size()); ++i )
    {
      const Entry &e = entries[i];
      if( e.hidden || !e.summary )
        continue;
      const vector<uint64_t> &aliases = e.summary->sourceAliasHashes;
      if( (e.summary->hash == summary->hash)
          || (std::find( begin(aliases), end(aliases), summary->hash ) != end(aliases)) )
        match = static_cast<int>( i );
    }

    if( match >= 0 )
    {
      entries[static_cast<size_t>(match)].suggestReasons |= suggestion.reasons;
      continue;
    }

    Entry e;
    e.kind = EntryKind::InMemory;
    e.id = "mem:" + std::to_string( summary->hash );
    e.inMemoryDrf = suggestion.drf;
    e.summary = summary;
    e.knownHashes.push_back( summary->hash );
    e.suggestReasons = suggestion.reasons;
    e.origin = (suggestion.reasons & SuggestInSpectrumFile) ? Origin::ThisSpectrum : Origin::UserMade;
    entries.push_back( e );
    index_entry( entries.size() - 1 );
  }//for( const Suggestion &suggestion : suggestions )

  return entries;
}//mergeEntries(...)


std::shared_ptr<const DetectorPeakResponse> loadDrf( const Entry &entry,
                                                     std::shared_ptr<DataBaseUtils::DbSession> sql )
{
  if( entry.kind == EntryKind::InMemory )
  {
    if( !entry.inMemoryDrf )
      throw runtime_error( "No detector." );
    return entry.inMemoryDrf;
  }

  if( !entry.summary )
    throw runtime_error( entry.error.empty() ? string("Detector not read yet.") : entry.error );

  Cache &cache = Cache::instance();
  const string key = (entry.kind == EntryKind::DbRow) ? ("db:" + std::to_string(entry.dbId)) : entry.id;
  shared_ptr<const DetectorPeakResponse> drf = cache.findDrf( key, entry.summary->hash );
  if( drf )
    return drf;

  switch( entry.kind )
  {
    case EntryKind::DbRow:
    {
      if( !sql || (entry.dbId < 0) )
        throw runtime_error( "No database." );

      DataBaseUtils::DbTransaction transaction( *sql );
      const Wt::Dbo::ptr<DetectorPeakResponse> row
                 = sql->session()->find<DetectorPeakResponse>().where( "id = ?" ).bind( entry.dbId ).resultValue();
      if( !row )
        throw runtime_error( "The detector is no longer in the database." );
      drf = make_shared<DetectorPeakResponse>( *row );
      transaction.commit();
      break;
    }//case EntryKind::DbRow:

    case EntryKind::GadrasFolder:
    case EntryKind::TableRow:
    case EntryKind::ImportFile:
    {
      WorkItem item;
      item.kind = entry.kind;
      item.path = entry.path;
      item.companion = entry.companionPath;

      Location loc;
      loc.builtIn = (entry.origin == Origin::BuiltIn);

      size_t num_unreadable = 0;
      string error;
      const vector<shared_ptr<DetectorPeakResponse>> drfs = read_item( item, loc, num_unreadable, error );
      if( entry.record >= drfs.size() )
        throw runtime_error( error.empty() ? string("Could not read the detector file.") : error );
      drf = drfs[entry.record];
      break;
    }//case file kinds

    case EntryKind::InMemory:
      break;
  }//switch( entry.kind )

  if( !drf )
    throw runtime_error( "Could not load the detector." );

  cache.storeDrf( key, drf );
  return drf;
}//loadDrf(...)


void invalidatePath( const std::string &path )
{
  Cache::instance().invalidatePath( path );
}


void invalidateDbRow( const long long dbId )
{
  Cache::instance().invalidateRow( dbId );
}


void clearCache()
{
  Cache::instance().clear();
}


size_t numSummariesComputed()
{
  return Cache::instance().m_numComputed.load();
}


bool matchesFilter( const Entry &entry, const Filter &filter, const ExtraSearchText &extraText )
{
  if( !filter.allOrigins && (entry.origin != filter.origin) )
    return false;

  if( filter.withGeometry && !(entry.summary && entry.summary->hasGeometry) )
    return false;

  vector<string> words;
  SpecUtils::split( words, lower_ascii( filter.query ), " \t" );
  if( words.empty() )
    return true;

  // The file's and its folder's names - not the whole path, whose user and app folders would make
  //  "users" or "library" match everything kept or added.
  string haystack = SpecUtils::filename( entry.path ) + " "
                    + SpecUtils::filename( SpecUtils::parent_path( entry.path ) );
  if( entry.summary )
    haystack += " " + entry.summary->name + " " + entry.summary->description
                + " " + entry.summary->crystal.material;
  if( extraText )
    haystack += " " + extraText( entry );
  SpecUtils::to_lower_ascii( haystack );

  for( const string &w : words )
  {
    if( !w.empty() && (haystack.find( w ) == string::npos) )
      return false;
  }

  return true;
}//matchesFilter(...)


int parentIndex( const std::vector<Entry> &entries, const size_t child )
{
  if( (child >= entries.size()) || !entries[child].summary || !entries[child].summary->parentHash )
    return -1;

  const uint64_t parent_hash = entries[child].summary->parentHash;
  for( size_t i = 0; i < entries.size(); ++i )
  {
    if( i == child )
      continue;
    const vector<uint64_t> &hashes = entries[i].knownHashes;
    if( std::find( begin(hashes), end(hashes), parent_hash ) != end(hashes) )
      return static_cast<int>( i );
  }

  return -1;
}//int parentIndex(...)


std::vector<Row> buildRows( const std::vector<Entry> &entries,
                            const Filter &filter,
                            const std::map<std::string,bool> &expanded,
                            const std::string &pickedId,
                            const std::string &inUseId,
                            const ExtraSearchText &extraText )
{
  const size_t n = entries.size();

  vector<bool> visible( n, false ), matches( n, false );
  map<uint64_t,size_t> by_hash;
  for( size_t i = 0; i < n; ++i )
  {
    const Entry &e = entries[i];
    visible[i] = !e.hidden && (e.summary || e.pending);
    if( !visible[i] )
      continue;
    matches[i] = matchesFilter( e, filter, extraText );
    for( const uint64_t h : e.knownHashes )
      by_hash.emplace( h, i );
  }//for( size_t i = 0; i < n; ++i )

  // Each copy's top-level ancestor among the visible entries, so a copy of a copy nests under the
  //  original too.  -1 when it is not a copy, or its parent is not listed.
  vector<int> root( n, -1 );
  for( size_t i = 0; i < n; ++i )
  {
    if( !visible[i] || !entries[i].summary || !entries[i].summary->parentHash )
      continue;

    size_t current = i;
    set<size_t> seen{ i };
    while( entries[current].summary && entries[current].summary->parentHash )
    {
      const auto pos = by_hash.find( entries[current].summary->parentHash );
      if( (pos == end(by_hash)) || !visible[pos->second] || !seen.insert( pos->second ).second )
        break;
      current = pos->second;
    }

    if( current != i )
      root[i] = static_cast<int>( current );
  }//for( size_t i = 0; i < n; ++i )

  map<size_t,vector<size_t>> copies_of;
  for( size_t i = 0; i < n; ++i )
    if( root[i] >= 0 )
      copies_of[static_cast<size_t>(root[i])].push_back( i );

  // Top-level rows: matching entries that are not shown under a matching original.
  vector<size_t> tops;
  for( size_t i = 0; i < n; ++i )
  {
    if( !visible[i] || !matches[i] )
      continue;
    if( (root[i] >= 0) && matches[static_cast<size_t>(root[i])] )
      continue;
    tops.push_back( i );
  }

  auto family_last_used = [&]( const size_t i ) -> int64_t {
    int64_t used = entries[i].lastUsedUtc;
    if( root[i] < 0 )
    {
      const auto pos = copies_of.find( i );
      if( pos != end(copies_of) )
        for( const size_t k : pos->second )
          used = std::max( used, entries[k].lastUsedUtc );
    }
    return used;
  };//family_last_used

  auto by_name = [&]( const size_t lhs, const size_t rhs ) -> bool {
    const string &l = entries[lhs].summary ? entries[lhs].summary->name : entries[lhs].path;
    const string &r = entries[rhs].summary ? entries[rhs].summary->name : entries[rhs].path;
    return name_less( l, r );
  };

  vector<size_t> suggested, recent, rest;
  for( const size_t i : tops )
  {
    if( entries[i].suggestReasons != SuggestNone )
      suggested.push_back( i );
    else if( family_last_used( i ) > 0 )
      recent.push_back( i );
    else
      rest.push_back( i );
  }

  std::stable_sort( begin(suggested), end(suggested), by_name );
  std::stable_sort( begin(recent), end(recent), [&]( const size_t lhs, const size_t rhs ){
    return family_last_used( lhs ) > family_last_used( rhs );
  } );
  if( recent.size() > filter.maxRecent )
  {
    rest.insert( end(rest), begin(recent) + filter.maxRecent, end(recent) );
    recent.resize( filter.maxRecent );
  }
  std::sort( begin(rest), end(rest), by_name );

  vector<Row> rows;

  auto push_family = [&]( const size_t i, const Group group ){
    Row row;
    row.type = RowType::Detector;
    row.group = group;
    row.entry = i;
    row.copyOfHidden = (entries[i].summary && entries[i].summary->parentHash);
    rows.push_back( row );

    if( root[i] >= 0 )
      return;  // a copy listed on its own carries no copies of its own

    const auto pos = copies_of.find( i );
    if( (pos == end(copies_of)) || pos->second.empty() )
      return;

    vector<size_t> kids = pos->second;
    std::sort( begin(kids), end(kids), by_name );

    bool open = (entries[i].id == pickedId);
    for( const size_t k : kids )
      open = open || (entries[k].id == pickedId) || (entries[k].id == inUseId);
    const auto exp = expanded.find( entries[i].id );
    if( exp != end(expanded) )
      open = exp->second;

    Row toggle;
    toggle.type = RowType::CopiesToggle;
    toggle.group = group;
    toggle.entry = i;
    toggle.numCopies = kids.size();
    toggle.expanded = open;
    rows.push_back( toggle );

    if( !open )
      return;

    for( const size_t k : kids )
    {
      Row kid;
      kid.type = RowType::Detector;
      kid.group = group;
      kid.entry = k;
      kid.nested = true;
      rows.push_back( kid );
    }
  };//push_family

  auto push_group = [&]( const Group group, const vector<size_t> &members ){
    if( members.empty() )
      return;
    Row header;
    header.type = RowType::GroupHeader;
    header.group = group;
    rows.push_back( header );
    for( const size_t i : members )
      push_family( i, group );
  };//push_group

  push_group( Group::Suggested, suggested );
  push_group( Group::Recent, recent );
  push_group( (suggested.empty() && recent.empty()) ? Group::All : Group::Everything, rest );

  return rows;
}//buildRows(...)

}//namespace DrfLibrary
