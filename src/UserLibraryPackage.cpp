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

#include <set>
#include <map>
#include <mutex>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <filesystem>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/ZipArchive.h"
#include "InterSpec/UserLibrary.h"
#include "InterSpec/UserLibraryPackage.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using nlohmann::json;


namespace
{
  std::filesystem::path to_fs_path( const string &utf8 )
  {
#ifdef _WIN32
    return std::filesystem::path( SpecUtils::convert_from_utf8_to_utf16( utf8 ) );
#else
    return std::filesystem::path( utf8 );
#endif
  }//to_fs_path(...)


  string normalize( const string &path )
  {
    return SpecUtils::lexically_normalize_path( path );
  }


  bool is_dot_name( const string &path )
  {
    const string name = SpecUtils::filename( path );
    return name.empty() || (name[0] == '.');
  }


  /** Every file below a folder (skipping names starting with a dot), sorted. */
  vector<string> all_files( const string &root )
  {
    vector<string> answer;
    std::function<void(const string &, int)> walk;
    walk = [&]( const string &dir, const int depth ){
      if( depth > 16 )
        return;
      for( const string &f : SpecUtils::ls_files_in_directory( dir ) )
        if( !is_dot_name( f ) )
          answer.push_back( f );
      for( const string &sub : SpecUtils::ls_directories_in_directory( dir ) )
        if( !is_dot_name( sub ) )
          walk( sub, depth + 1 );
    };
    if( SpecUtils::is_directory( root ) )
      walk( root, 0 );
    std::sort( begin(answer), end(answer) );
    return answer;
  }//all_files(...)


  string read_whole_file( const string &path )
  {
    std::ifstream input( to_fs_path( path ), ios::in | ios::binary );
    if( !input )
      throw runtime_error( "Could not read '" + path + "'" );
    std::ostringstream strm;
    strm << input.rdbuf();
    if( !input && !input.eof() )
      throw runtime_error( "Could not read '" + path + "'" );
    return strm.str();
  }//read_whole_file(...)


  void remove_folder( const string &path )
  {
    if( path.empty() )
      return;
    std::error_code ec;
    std::filesystem::remove_all( to_fs_path( path ), ec );
    if( ec )
      cerr << "UserLibraryPackage: could not remove '" << path << "': " << ec.message() << endl;
  }//remove_folder(...)


  /** A sibling of `packages/<folder>` to move things through, e.g. ".incoming-1a2b3c4d". */
  /** Installs and uninstalls (from any session) run one at a time. */
  std::mutex ns_install_mutex;


  string hidden_folder_path( const string &prefix )
  {
    return SpecUtils::append_path( UserLibrary::packagesDirectory(), prefix + UserLibrary::newUuid().substr( 0, 8 ) );
  }


  void rename_folder( const string &from, const string &to )
  {
    std::error_code ec;
    std::filesystem::rename( to_fs_path( from ), to_fs_path( to ), ec );
    if( ec )
      throw runtime_error( "Could not move '" + from + "' to '" + to + "': " + ec.message() );
  }//rename_folder(...)


  /** The SHA-256 of whatever is at `path` (file or folder), or empty if it can not be read. */
  string hash_or_empty( const string &path )
  {
    try
    {
      UserLibrary::FoundItem item;
      item.path = path;
      item.isFolder = SpecUtils::is_directory( path );
      return UserLibrary::sha256OfItem( item );
    }catch( std::exception & )
    {
    }
    return string();
  }//hash_or_empty(...)


  /** The folder an item's package path is relative to: its type's folder (or legacy folder) in the
   user data folder or package; for DRFs placed elsewhere in the user data folder, that folder itself.
   */
  string item_root_dir( const UserLibrary::Item &item )
  {
    const UserLibrary::AssetType * const type = UserLibrary::assetType( item.typeId );
    if( !type )
      throw runtime_error( "Unknown item type '" + item.typeId + "'" );

    string base;
    switch( item.source.kind )
    {
      case UserLibrary::SourceKind::User:
#if( USE_USER_LIBRARY )
        base = InterSpec::writableDataDirectory();
#endif
        break;

      case UserLibrary::SourceKind::Package:
      {
        const optional<UserLibrary::Package> pkg = UserLibrary::installedPackage( item.source.packageId );
        if( !pkg )
          throw runtime_error( "The package of '" + item.relPath + "' is no longer installed." );
        base = pkg->path;
        break;
      }

      case UserLibrary::SourceKind::BuiltIn:
        throw runtime_error( "InterSpec's own items can not be put in a package." );
    }//switch( item.source.kind )

    string best;
    vector<string> folders = type->legacyFolders;
    folders.push_back( type->folder );
    for( const string &f : folders )
    {
      if( (f.size() > best.size()) && SpecUtils::istarts_with( item.relPath, f + "/" ) )
        best = f;
    }
    return best.empty() ? base : SpecUtils::append_path( base, best );
  }//item_root_dir(...)
}//namespace


namespace UserLibraryPackage
{
const char * const FileExtension = ".ispkg";


bool looksLikePackage( const std::string &displayName, const std::string &path )
{
  if( SpecUtils::iends_with( displayName, FileExtension ) )
    return true;

  if( !SpecUtils::iends_with( displayName, ".zip" ) )
    return false;

  try
  {
    std::ifstream input( to_fs_path( path ), ios::in | ios::binary );
    if( !input )
      return false;
    for( const shared_ptr<const ZipArchive::ZipFileHeader> &header : ZipArchive::list_zip_entries( input ) )
    {
      if( SpecUtils::iequals_ascii( header->filename, UserLibrary::ManifestFileName ) )
        return true;
    }
  }catch( std::exception & )
  {
  }
  return false;
}//looksLikePackage(...)


std::string checkEntryName( const std::string &name, const Limits &limits )
{
  if( name.empty() )
    throw runtime_error( "An entry has no name." );

  string path = name;
  SpecUtils::ireplace_all( path, "\\", "/" );

  if( (path[0] == '/') || ((path.size() >= 2) && std::isalpha( static_cast<unsigned char>(path[0]) ) && (path[1] == ':')) )
    throw runtime_error( "'" + name + "' is an absolute path." );

  while( !path.empty() && (path.back() == '/') )
    path.pop_back();

  vector<string> parts;
  SpecUtils::split_no_delim_compress( parts, path, "/" );
  if( parts.empty() )
    throw runtime_error( "An entry has no name." );

  for( const string &part : parts )
  {
    if( part.empty() )
      throw runtime_error( "'" + name + "' has an empty folder name." );
    if( (part == ".") || (part == "..") )
      throw runtime_error( "'" + name + "' refers outside its folder." );
    if( part.size() > limits.maxComponentLength )
      throw runtime_error( "'" + name + "' has a name that is too long." );
    if( (part.back() == '.') || (part.back() == ' ') )
      throw runtime_error( "'" + name + "' has a name ending in a dot or space, which Windows does not allow." );

    for( const char c : part )
    {
      const unsigned char uc = static_cast<unsigned char>( c );
      if( (uc < 0x20) || (uc == 0x7F) || (string( ":*?\"<>|" ).find( c ) != string::npos) )
        throw runtime_error( "'" + name + "' has a character file systems do not allow." );
    }

    string stem = part.substr( 0, part.find( '.' ) );
    SpecUtils::to_lower_ascii( stem );
    const bool reserved = (stem == "con") || (stem == "prn") || (stem == "aux") || (stem == "nul")
                          || ((stem.size() == 4) && (SpecUtils::starts_with( stem, "com" ) || SpecUtils::starts_with( stem, "lpt" ))
                              && (stem[3] >= '0') && (stem[3] <= '9'));
    if( reserved )
      throw runtime_error( "'" + name + "' uses a name Windows reserves for devices." );
  }//for( loop over path components )

  // Valid UTF-8 (zip names may be in a legacy code page; those are not accepted).
  for( size_t i = 0; i < path.size(); )
  {
    const unsigned char c = static_cast<unsigned char>( path[i] );
    size_t len = 0;
    if( c < 0x80 )
      len = 1;
    else if( (c & 0xE0) == 0xC0 )
      len = 2;
    else if( (c & 0xF0) == 0xE0 )
      len = 3;
    else if( (c & 0xF8) == 0xF0 )
      len = 4;
    if( !len || ((i + len) > path.size()) )
      throw runtime_error( "'" + name + "' is not a valid UTF-8 name." );
    for( size_t j = 1; j < len; ++j )
    {
      if( (static_cast<unsigned char>( path[i+j] ) & 0xC0) != 0x80 )
        throw runtime_error( "'" + name + "' is not a valid UTF-8 name." );
    }
    i += len;
  }//for( check UTF-8 )

  string answer;
  for( const string &part : parts )
    answer += (answer.empty() ? "" : "/") + part;

  if( answer.size() > limits.maxPathLength )
    throw runtime_error( "'" + name + "' is too long a path." );

  return answer;
}//checkEntryName(...)


std::shared_ptr<Inspection> inspect( const std::string &archivePath,
                                     const std::string &displayName,
                                     const Limits &limits )
{
  shared_ptr<Inspection> insp = make_shared<Inspection>();
  insp->displayName = displayName;

  try
  {
    if( !UserLibrary::available() )
      throw runtime_error( "This InterSpec has no user data folder to install packages into." );

    if( SpecUtils::file_size( archivePath ) > (limits.maxTotalBytes + 16*1024*1024) )
      throw runtime_error( "The package is larger than InterSpec accepts." );

    insp->archiveSha256 = UserLibrary::sha256OfFile( archivePath );

    std::ifstream input( to_fs_path( archivePath ), ios::in | ios::binary );
    if( !input )
      throw runtime_error( "Could not open the package file." );

    vector<shared_ptr<const ZipArchive::ZipFileHeader>> entries;
    try
    {
      entries = ZipArchive::list_zip_entries( input );
    }catch( std::exception & )
    {
      throw runtime_error( "The file is not a valid package (zip) archive." );
    }

    if( entries.size() > limits.maxEntries )
      throw runtime_error( "The package has more files than InterSpec accepts." );

    // Check every entry before writing anything: one bad entry rejects the whole package.
    vector<pair<shared_ptr<const ZipArchive::ZipFileHeader>,string>> accepted;
    set<string> collision_keys;
    size_t total_bytes = 0;
    for( const shared_ptr<const ZipArchive::ZipFileHeader> &header : entries )
    {
      if( header->is_symlink() )
        throw runtime_error( "The package contains a symbolic link ('" + header->filename + "')." );

      const string rel = checkEntryName( header->filename, limits );
      if( header->is_directory() )
        continue;

      if( !collision_keys.insert( UserLibrary::collisionKey( rel ) ).second )
        throw runtime_error( "The package has two files named '" + rel
                             + "' (or names differing only in capitalization or accents)." );

      if( header->uncompressed_size > limits.maxFileBytes )
        throw runtime_error( "'" + rel + "' is larger than InterSpec accepts." );
      total_bytes += header->uncompressed_size;
      if( total_bytes > limits.maxTotalBytes )
        throw runtime_error( "The package is larger, unpacked, than InterSpec accepts." );

      if( UserLibrary::isSecretPath( rel ) )
      {
        SkippedFile skipped;
        skipped.relPath = rel;
        skipped.refused = true;
        insp->skipped.push_back( skipped );
        continue;
      }

      // Hidden files (e.g. ".DS_Store", "x/.hidden/SKILL.md") and very deeply nested ones are not
      //  unpacked: nothing would list, check, or remove them (see all_files()), but some readers would.
      vector<string> components;
      SpecUtils::split( components, rel, "/" );
      bool hidden = false;
      for( const string &c : components )
        hidden = hidden || (!c.empty() && (c.front() == '.'));
      if( hidden || (components.size() > 17) )
      {
        SkippedFile skipped;
        skipped.relPath = rel;
        insp->skipped.push_back( skipped );
        continue;
      }

      accepted.emplace_back( header, rel );
    }//for( loop over entries )

    // Unpack into a hidden folder beside the installed packages (same volume, so installing is a rename).
    const string pkgs_dir = UserLibrary::packagesDirectory();
    {
      std::error_code ec;
      std::filesystem::create_directories( to_fs_path( pkgs_dir ), ec );
      if( !SpecUtils::is_directory( pkgs_dir ) )
        throw runtime_error( "Could not create the folder '" + pkgs_dir + "'." );
    }
    insp->stagingDir = hidden_folder_path( ".incoming-" );

    ZipArchive::ExtractionLimits extract_limits;
    extract_limits.max_uncompressed_bytes = limits.maxFileBytes;
    extract_limits.max_compression_ratio = limits.maxCompressionRatio;

    for( const pair<shared_ptr<const ZipArchive::ZipFileHeader>,string> &entry : accepted )
    {
      const string dest = SpecUtils::append_path( insp->stagingDir, entry.second );
      std::error_code ec;
      std::filesystem::create_directories( to_fs_path( SpecUtils::parent_path( dest ) ), ec );
      if( ec )
        throw runtime_error( "Could not unpack '" + entry.second + "': " + ec.message() );

      std::ofstream output( to_fs_path( dest ), ios::out | ios::binary | ios::trunc );
      if( !output )
        throw runtime_error( "Could not unpack '" + entry.second + "'." );
      if( entry.first->uncompressed_size > 0 )
      {
        try
        {
          ZipArchive::read_file_from_zip( input, entry.first, output, extract_limits );
        }catch( std::exception &e )
        {
          throw runtime_error( "Could not unpack '" + entry.second + "': " + string( e.what() ) );
        }
      }
      output.close();
      if( !output )
        throw runtime_error( "Could not write '" + entry.second + "'." );
    }//for( loop over accepted entries )

    // The manifest; a package made by hand may not have one.
    const string manifest_path = SpecUtils::append_path( insp->stagingDir, UserLibrary::ManifestFileName );
    if( SpecUtils::is_file( manifest_path ) )
    {
      if( SpecUtils::file_size( manifest_path ) > 4*1024*1024 )
        throw runtime_error( "The package's " + string( UserLibrary::ManifestFileName ) + " is too large." );
      try
      {
        insp->manifest = UserLibrary::Manifest::fromJson( read_whole_file( manifest_path ) );
      }catch( std::exception &e )
      {
        throw runtime_error( "The package's " + string( UserLibrary::ManifestFileName ) + " is not valid: " + e.what() );
      }
      insp->hadManifest = true;

      if( insp->manifest.formatVersion > UserLibrary::Manifest::sm_formatVersion )
        throw runtime_error( "This package was made by a newer InterSpec (package format "
                             + std::to_string( insp->manifest.formatVersion ) + "); upgrade InterSpec to install it." );
    }else
    {
      string name = SpecUtils::filename( displayName );
      const string ext = SpecUtils::file_extension( name );
      name = name.substr( 0, name.size() - ext.size() );
      insp->manifest.id = UserLibrary::newUuid();
      insp->manifest.name = name.empty() ? string( "Package" ) : name;
      insp->manifest.version = "1";
      insp->hadManifest = false;
    }//if( manifest ) / else

    // The items, by type folder.  Files that are no item of a known type are not installed.
    set<string> covered{ normalize( manifest_path ) };
    for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
    {
      vector<string> folders{ type.folder };
      folders.insert( end(folders), begin(type.legacyFolders), end(type.legacyFolders) );
      for( const string &folder : folders )
      {
        const string root = SpecUtils::append_path( insp->stagingDir, folder );
        if( !SpecUtils::is_directory( root ) || !type.findItems )
          continue;

        for( const UserLibrary::FoundItem &found : type.findItems( root ) )
        {
          IncomingItem item;
          item.found = found;
          item.typeId = type.id;
          item.relPath = UserLibrary::relativePath( insp->stagingDir, found.path );
          if( item.relPath.empty() || covered.count( normalize( found.path ) ) )
            continue;

          if( found.isFolder )
          {
            for( const string &f : all_files( found.path ) )
              covered.insert( normalize( f ) );
          }
          covered.insert( normalize( found.path ) );
          for( const string &s : found.support )
            covered.insert( normalize( s ) );

          insp->items.push_back( item );
        }//for( loop over found items )
      }//for( loop over the type's folders )
    }//for( loop over types )

    for( const string &f : all_files( insp->stagingDir ) )
    {
      if( covered.count( normalize( f ) ) )
        continue;
      SkippedFile skipped;
      skipped.relPath = UserLibrary::relativePath( insp->stagingDir, f );
      insp->skipped.push_back( skipped );
      SpecUtils::remove_file( f );
    }//for( loop over unpacked files )

    // Hashes, manifest details, and reading each item with its real reader.
    const bool too_new = insp->manifest.needsNewerInterSpec();
    for( IncomingItem &item : insp->items )
    {
      const UserLibrary::ManifestItem * const mi = insp->manifest.itemForPath( item.relPath );
      if( mi )
      {
        item.itemId = mi->itemId;
        item.notes = mi->notes;
        item.intendedOverride = mi->intendedOverride;
      }
      if( item.itemId.empty() && !insp->hadManifest )
        item.itemId = UserLibrary::newUuid();

      item.sha256 = UserLibrary::sha256OfItem( item.found );
      item.hashMismatch = mi && !mi->sha256.empty() && !SpecUtils::iequals_ascii( mi->sha256, item.sha256 );

      const UserLibrary::AssetType * const type = UserLibrary::assetType( item.typeId );
      item.validation = type->validate( item.found );
      if( too_new && (item.validation.status != UserLibrary::Status::Unreadable) )
      {
        item.validation.status = UserLibrary::Status::NeedsNewer;
        item.validation.message = "Needs InterSpec " + insp->manifest.minInterSpecVersion + " or newer.";
      }

      if( item.validation.status == UserLibrary::Status::NeedsNewer )
        item.state = ItemState::NeedsNewer;
      else if( item.validation.status == UserLibrary::Status::Unreadable )
        item.state = ItemState::Unreadable;
    }//for( loop over items )

    if( !insp->hadManifest )
    {
      for( const IncomingItem &item : insp->items )
      {
        UserLibrary::ManifestItem mi;
        mi.itemId = item.itemId;
        mi.path = item.relPath;
        mi.type = item.typeId;
        mi.sha256 = item.sha256;
        insp->manifest.items.push_back( mi );
      }
    }//if( !insp->hadManifest )

    // Compared with an installed version of the same package.
    insp->installed = UserLibrary::installedPackage( insp->manifest.id );
    if( insp->installed )
    {
      const UserLibrary::Package &installed = *insp->installed;
      const UserLibrary::Manifest &old_manifest = installed.manifest;

      // Each installed item's current content (by item_id, else path).
      map<string,string> old_hash_by_id, old_hash_by_path;
      set<string> old_matched;
      for( const UserLibrary::ManifestItem &mi : old_manifest.items )
      {
        const string hash = hash_or_empty( SpecUtils::append_path( installed.path, mi.path ) );
        if( !mi.itemId.empty() )
          old_hash_by_id[mi.itemId] = hash;
        old_hash_by_path[mi.path] = hash;
      }

      bool all_same = (insp->items.size() == old_manifest.items.size());
      for( IncomingItem &item : insp->items )
      {
        const UserLibrary::ManifestItem *old_item = nullptr;
        for( const UserLibrary::ManifestItem &mi : old_manifest.items )
        {
          if( (!item.itemId.empty() && (mi.itemId == item.itemId)) || (!old_item && (mi.path == item.relPath)) )
            old_item = &mi;
        }

        if( !old_item )
        {
          insp->added.push_back( item.relPath );
          all_same = false;
          continue;
        }

        old_matched.insert( old_item->itemId.empty() ? old_item->path : old_item->itemId );
        const string old_hash = old_hash_by_path[old_item->path];
        if( old_hash == item.sha256 )
        {
          if( item.state == ItemState::New )
            item.state = ItemState::Unchanged;
        }else
        {
          insp->changed.push_back( item.relPath );
          all_same = false;
          if( item.state == ItemState::New )
            item.state = ItemState::Update;
        }
      }//for( loop over incoming items )

      for( const UserLibrary::ManifestItem &mi : old_manifest.items )
      {
        if( !old_matched.count( mi.itemId.empty() ? mi.path : mi.itemId ) )
        {
          insp->removed.push_back( mi.path );
          all_same = false;
        }
      }

      const int cmp = UserLibrary::compareVersions( insp->manifest.version, old_manifest.version );
      if( cmp > 0 )
        insp->relation = Relation::Update;
      else if( cmp < 0 )
        insp->relation = Relation::Downgrade;
      else
        insp->relation = all_same ? Relation::Reinstall : Relation::SameVersionDifferent;

      for( const pair<const string,string> &f : installed.files )
      {
        if( hash_or_empty( SpecUtils::append_path( installed.path, f.first ) ) != f.second )
          insp->locallyModified.push_back( f.first );
      }

      insp->authorChanged = (SpecUtils::trim_copy( insp->manifest.author ) != SpecUtils::trim_copy( old_manifest.author ))
                  || (SpecUtils::trim_copy( insp->manifest.organization ) != SpecUtils::trim_copy( old_manifest.organization ));

      vector<UserLibrary::ManifestPreviousVersion> versions;
      versions.push_back( { insp->manifest.version, insp->manifest.created, insp->manifest.changes } );
      versions.insert( end(versions), begin(insp->manifest.previousVersions), end(insp->manifest.previousVersions) );
      for( const UserLibrary::ManifestPreviousVersion &v : versions )
      {
        if( UserLibrary::compareVersions( v.version, old_manifest.version ) > 0 )
          insp->changesSince.push_back( v );
      }
      std::stable_sort( begin(insp->changesSince), end(insp->changesSince),
        []( const UserLibrary::ManifestPreviousVersion &lhs, const UserLibrary::ManifestPreviousVersion &rhs ){
          return UserLibrary::compareVersions( lhs.version, rhs.version ) > 0;
      } );
    }//if( insp->installed )

    // Same-named items already installed (other than this package's own earlier version).
    const map<string,UserLibrary::Choice> choices = UserLibrary::choices();
    map<string,vector<UserLibrary::Item>> existing_by_type;
    for( IncomingItem &item : insp->items )
    {
      if( (item.validation.status != UserLibrary::Status::Ok) || item.validation.names.empty() )
        continue;

      const UserLibrary::AssetType * const type = UserLibrary::assetType( item.typeId );
      if( !existing_by_type.count( item.typeId ) )
        existing_by_type[item.typeId] = UserLibrary::items( item.typeId, true );
      const vector<UserLibrary::Item> &existing = existing_by_type[item.typeId];

      for( const string &name : item.validation.names )
      {
        for( const UserLibrary::Item &ex : existing )
        {
          if( ex.disabled || ex.needsNewer
              || ((ex.source.kind == UserLibrary::SourceKind::Package) && (ex.source.packageId == insp->manifest.id)) )
            continue;

          bool same_name = false;
          for( const string &n : ex.names )
            same_name = same_name || UserLibrary::sameName( *type, n, name );
          if( !same_name )
            continue;

          NameConflict conflict;
          conflict.name = name;
          conflict.existing = ex;
          conflict.overridesBuiltIn = (ex.source.kind == UserLibrary::SourceKind::BuiltIn);
          conflict.identical = (hash_or_empty( ex.found.path ) == item.sha256);
          conflict.useIncomingByDefault = item.intendedOverride;

          const string prefix = item.typeId + ":";
          for( const pair<const string,UserLibrary::Choice> &c : choices )
          {
            if( SpecUtils::starts_with( c.first, prefix.c_str() )
                && UserLibrary::sameName( *type, c.first.substr( prefix.size() ), name ) )
            {
              conflict.alreadyChosen = true;
              conflict.useIncomingByDefault = (c.second.source == ("package:" + insp->manifest.id));
            }
          }//for( loop over choices )

          item.conflicts.push_back( conflict );

          // Only the one in use now: items are in precedence order (the user's, packages', then
          //  InterSpec's), and keeping "the existing one" means keeping that one.
          break;
        }//for( loop over existing items )
      }//for( loop over the item's names )

      if( item.state == ItemState::New )
      {
        bool any_open = false, any_identical = false;
        for( const NameConflict &c : item.conflicts )
        {
          any_open = any_open || (!c.identical && !c.alreadyChosen);
          any_identical = any_identical || c.identical;
        }
        if( any_open )
          item.state = ItemState::Conflict;
        else if( any_identical )
          item.state = ItemState::AlreadyHave;
      }//if( item.state == ItemState::New )
    }//for( loop over items )
  }catch( std::exception &e )
  {
    insp->error = e.what();
    remove_folder( insp->stagingDir );
    insp->stagingDir.clear();
  }//try / catch

  return insp;
}//inspect(...)


UserLibrary::Package install( Inspection &insp, const Decisions &decisions )
{
  if( !insp.error.empty() )
    throw runtime_error( insp.error );
  if( insp.stagingDir.empty() || !SpecUtils::is_directory( insp.stagingDir ) )
    throw runtime_error( "The unpacked package is gone; drop the file again." );

  std::lock_guard<std::mutex> install_lock( ns_install_mutex );

  const string pkgs_dir = UserLibrary::packagesDirectory();
  const string id = insp.manifest.id;

  // What is installed now, not when the package was inspected: another install of this package
  //  (say, the same file opened twice) may have finished since, and this then updates it.
  insp.installed = UserLibrary::installedPackage( id );

  // A package made without a manifest gets the one made up for it, so its folder says what it is.
  if( !insp.hadManifest )
  {
    const string path = SpecUtils::append_path( insp.stagingDir, UserLibrary::ManifestFileName );
    std::ofstream output( to_fs_path( path ), ios::out | ios::binary | ios::trunc );
    output << insp.manifest.toJson();
    if( !output )
      throw runtime_error( "Could not write the package's manifest." );
  }

  // Each file's hash as installed, to notice later edits.
  json files = json::object();
  for( const string &f : all_files( insp.stagingDir ) )
  {
    const string rel = UserLibrary::relativePath( insp.stagingDir, f );
    if( rel != UserLibrary::ManifestFileName )
      files[rel] = UserLibrary::sha256OfFile( f );
  }

  // An update keeps its folder; a new package gets one named after it.
  string folder;
  if( insp.installed && SpecUtils::is_directory( insp.installed->path ) )
  {
    folder = insp.installed->folder;
  }else
  {
    const string base = UserLibrary::sanitizeFileName( insp.manifest.name, "Package" );
    folder = base;
    for( int i = 2; SpecUtils::is_directory( SpecUtils::append_path( pkgs_dir, folder ) )
                    || SpecUtils::is_file( SpecUtils::append_path( pkgs_dir, folder ) ); ++i )
      folder = base + " (" + std::to_string( i ) + ")";
  }

  const string final_path = SpecUtils::append_path( pkgs_dir, folder );
  string old_path;
  if( SpecUtils::is_directory( final_path ) )
  {
    old_path = hidden_folder_path( ".old-" );
    rename_folder( final_path, old_path );
  }

  try
  {
    rename_folder( insp.stagingDir, final_path );
  }catch( std::exception & )
  {
    if( !old_path.empty() )
    {
      try
      {
        rename_folder( old_path, final_path );
      }catch( std::exception &e )
      {
        cerr << "UserLibraryPackage::install(): could not put back '" << final_path << "': " << e.what() << endl;
      }
    }
    throw;
  }//try / catch

  insp.stagingDir.clear();
  remove_folder( old_path );

  // Hashes of the existing items the user kept, for the choices.
  map<string,string> kept_hashes;
  for( const IncomingItem &item : insp.items )
  {
    for( const NameConflict &c : item.conflicts )
    {
      if( !c.identical && !c.alreadyChosen )
        kept_hashes[c.existing.found.path] = hash_or_empty( c.existing.found.path );
    }
  }

  try
  {
    const string now = UserLibrary::nowIsoUtc();
    UserLibrary::modifyHistory( [&]( json &history ){
      json &pkgs = history["packages"];
      if( !pkgs.is_object() )
        pkgs = json::object();
      json &rec = pkgs[id];
      if( !rec.is_object() )
        rec = json::object();

      if( !rec.contains( "first_installed" ) || (rec.value( "folder", string() ) != folder) )
        rec["first_installed"] = now;
      rec["folder"] = folder;

      json &installs = rec["installs"];
      if( !installs.is_array() )
        installs = json::array();
      installs.push_back( { {"installed", now}, {"archive_sha256", insp.archiveSha256},
                            {"manifest", insp.manifest.toJsonObject()} } );
      while( installs.size() > 20 )  //the history file stays small; the package's own manifest is on disk
        installs.erase( installs.begin() );

      rec["files"] = files;
      if( !rec.contains( "disabled" ) || !rec["disabled"].is_array() )
        rec["disabled"] = json::array();

      json &choices = history["choices"];
      if( !choices.is_object() )
        choices = json::object();

      for( const IncomingItem &item : insp.items )
      {
        const UserLibrary::AssetType * const type = UserLibrary::assetType( item.typeId );
        for( const NameConflict &c : item.conflicts )
        {
          if( c.identical || c.alreadyChosen || !type )
            continue;

          const Decisions::const_iterator d = decisions.find( { item.typeId, c.name } );
          const bool use_incoming = (d != end(decisions)) ? d->second : c.useIncomingByDefault;

          vector<string> to_remove;
          const string prefix = item.typeId + ":";
          for( json::iterator it = choices.begin(); it != choices.end(); ++it )
          {
            if( SpecUtils::starts_with( it.key(), prefix.c_str() )
                && UserLibrary::sameName( *type, it.key().substr( prefix.size() ), c.name ) )
              to_remove.push_back( it.key() );
          }
          for( const string &k : to_remove )
            choices.erase( k );

          choices[prefix + c.name] = {
            {"source", use_incoming ? ("package:" + id) : c.existing.source.key()},
            {"sha256", use_incoming ? item.sha256 : kept_hashes[c.existing.found.path]},
            {"chosen", now}
          };
        }//for( loop over conflicts )
      }//for( loop over items )
    } );
  }catch( std::exception &e )
  {
    // The package is in place; it shows as "added manually" until recorded.
    cerr << "UserLibraryPackage::install(): installed '" << final_path << "', but could not record it: "
         << e.what() << endl;
  }

  UserLibrary::notifyChanged();

  const optional<UserLibrary::Package> pkg = UserLibrary::installedPackage( id );
  if( pkg )
    return *pkg;

  UserLibrary::Package answer;
  answer.id = id;
  answer.folder = folder;
  answer.path = final_path;
  answer.manifest = insp.manifest;
  answer.recorded = false;
  return answer;
}//install(...)


void discard( Inspection &insp )
{
  remove_folder( insp.stagingDir );
  insp.stagingDir.clear();
}//discard(...)


void uninstall( const std::string &packageId )
{
  std::lock_guard<std::mutex> install_lock( ns_install_mutex );

  const optional<UserLibrary::Package> pkg = UserLibrary::installedPackage( packageId );
  if( !pkg )
    throw runtime_error( "That package is not installed." );

  // Moved aside first, so a delete that fails half-way never leaves half a package in use.
  const string old_path = hidden_folder_path( ".old-" );
  rename_folder( pkg->path, old_path );
  remove_folder( old_path );

  UserLibrary::modifyHistory( [&]( json &history ){
    json &pkgs = history["packages"];
    if( pkgs.is_object() )
      pkgs.erase( packageId );

    json &choices = history["choices"];
    if( choices.is_object() )
    {
      vector<string> to_remove;
      for( json::iterator it = choices.begin(); it != choices.end(); ++it )
      {
        if( it->is_object() && (it->value( "source", string() ) == ("package:" + packageId)) )
          to_remove.push_back( it.key() );
      }
      for( const string &k : to_remove )
        choices.erase( k );
    }//if( choices.is_object() )
  } );
}//uninstall(...)


std::string sourceKey( const UserLibrary::Item &item )
{
  switch( item.source.kind )
  {
    case UserLibrary::SourceKind::User:    return "user:" + item.relPath;
    case UserLibrary::SourceKind::Package: return "package:" + item.source.packageId + ":" + item.relPath;
    case UserLibrary::SourceKind::BuiltIn: return "builtin:" + item.relPath;
  }
  return "user:" + item.relPath;
}//sourceKey(...)


std::string sourceKey( const CreateItem &item )
{
  return item.sourceKey.empty() ? sourceKey( item.item ) : item.sourceKey;
}


std::string packagePathFor( const UserLibrary::Item &item )
{
  CreateItem ci;
  ci.item = item;
  return packagePathFor( ci );
}


std::string packagePathFor( const CreateItem &ci )
{
  const UserLibrary::Item &item = ci.item;
  const UserLibrary::AssetType * const type = UserLibrary::assetType( item.typeId );
  if( !type )
    throw runtime_error( "Unknown item type '" + item.typeId + "'" );

  if( ci.databaseDrf )
    return type->folder + "/" + UserLibrary::sanitizeFileName( ci.databaseDrf->name(), "detector" ) + ".drf.xml";
  if( ci.generatedData )
    return type->folder + "/" + ci.generatedName;

  const string root = ci.root.empty() ? item_root_dir( item ) : ci.root;
  const string rel = UserLibrary::relativePath( root, item.found.path );
  if( rel.empty() )
    throw runtime_error( "Could not place '" + item.found.path + "' in the package." );
  return type->folder + "/" + rel;
}//packagePathFor(...)


std::vector<std::pair<std::string,std::string>> packageFilesFor( const UserLibrary::Item &item )
{
  CreateItem ci;
  ci.item = item;
  return packageFilesFor( ci );
}


std::string drfXml( const DetectorPeakResponse &drf )
{
  rapidxml::xml_document<char> doc;
  drf.toXml( &doc, &doc );
  string answer;
  rapidxml::print( std::back_inserter( answer ), doc, 0 );
  return answer;
}//drfXml(...)


std::vector<std::pair<std::string,std::string>> packageFilesFor( const CreateItem &ci )
{
  const UserLibrary::Item &item = ci.item;
  const UserLibrary::AssetType * const type = UserLibrary::assetType( item.typeId );
  if( !type )
    throw runtime_error( "Unknown item type '" + item.typeId + "'" );

  vector<pair<string,string>> answer;
  if( ci.databaseDrf || ci.generatedData )
    return answer;  //written from memory

  const string root = ci.root.empty() ? item_root_dir( item ) : ci.root;
  set<string> seen;
  auto add = [&]( const string &abs ){
    const string rel = UserLibrary::relativePath( root, abs );
    if( rel.empty() || !seen.insert( rel ).second )
      return;  // outside the item's folder (e.g. a DETECTOR.txt above it): not packaged
    answer.emplace_back( abs, type->folder + "/" + rel );
  };

  if( item.found.isFolder )
  {
    for( const string &f : all_files( item.found.path ) )
      add( f );
  }else
  {
    add( item.found.path );
  }

  for( const string &s : item.found.support )
    add( s );

  return answer;
}//packageFilesFor(...)


std::vector<PreflightIssue> preflight( const CreateRequest &request, const Limits &limits )
{
  vector<PreflightIssue> issues;
  auto add = [&issues]( const string &msg, const bool blocking ){
    PreflightIssue issue;
    issue.message = msg;
    issue.blocking = blocking;
    issues.push_back( issue );
  };

  if( SpecUtils::trim_copy( request.manifest.name ).empty() )
    add( "The package needs a name.", true );
  if( SpecUtils::trim_copy( request.manifest.version ).empty() )
    add( "The package needs a version.", true );
  if( request.items.empty() )
    add( "The package has no items.", true );

  // A new version of a package made before must have a later version number, or people who have
  //  the earlier one would see different content under the same version.
  if( !request.manifest.id.empty() )
  {
    const map<string,vector<UserLibrary::CreatedVersion>> created = UserLibrary::createdPackages();
    const map<string,vector<UserLibrary::CreatedVersion>>::const_iterator pos = created.find( request.manifest.id );
    if( (pos != end(created)) && !pos->second.empty() )
    {
      const string &last = pos->second.back().manifest.version;
      if( UserLibrary::compareVersions( request.manifest.version, last ) <= 0 )
        add( "Version " + last + " of this package was already made; give this one a later version (e.g. "
             + nextVersion( last ) + ").", true );
    }
  }//if( a new version of a package made before )

  map<pair<string,string>,string> names_seen;  // (type, name key) -> item path
  map<string,string> paths_seen;               // collision key -> source file
  size_t num_files = 0, total_bytes = 0;

  for( const CreateItem &ci : request.items )
  {
    const UserLibrary::Item &item = ci.item;
    const UserLibrary::AssetType * const type = UserLibrary::assetType( item.typeId );
    if( !type )
    {
      add( "'" + item.relPath + "' is of an unknown type.", true );
      continue;
    }

    if( item.source.kind == UserLibrary::SourceKind::BuiltIn )
    {
      add( "'" + item.relPath + "' comes with InterSpec, and can not be put in a package.", true );
      continue;
    }

    // A detector efficiency from the database: written as XML, so only its name and path to check.
    if( ci.databaseDrf )
    {
      const string name = ci.databaseDrf->name();
      const map<pair<string,string>,string>::const_iterator pos = names_seen.find( { type->id, name } );
      if( pos != end(names_seen) )
        add( "Two items are named '" + name + "' ('" + pos->second + "' and the one saved in InterSpec).", true );
      else
        names_seen[{ type->id, name }] = name;

      if( !ci.databaseDrf->isValid() )
        add( "The detector efficiency '" + name + "' saved in InterSpec is not complete.", true );

      try
      {
        const string path = checkEntryName( packagePathFor( ci ), limits );
        const string key = UserLibrary::collisionKey( path );
        if( paths_seen.count( key ) )
          add( "Two files would be '" + path + "' in the package.", true );
        paths_seen[key] = "database:" + name;
        num_files += 1;
      }catch( std::exception &e )
      {
        add( e.what(), true );
      }
      continue;
    }//if( ci.databaseDrf )

    // An item written from memory (e.g. a colour theme from the database): just its name and path.
    if( ci.generatedData )
    {
      for( const string &name : item.names )
      {
        const map<pair<string,string>,string>::const_iterator pos = names_seen.find( { type->id, name } );
        if( pos != end(names_seen) )
          add( "Two items are named '" + name + "'.", true );
        else
          names_seen[{ type->id, name }] = name;
      }
      try
      {
        const string path = checkEntryName( packagePathFor( ci ), limits );
        const string key = UserLibrary::collisionKey( path );
        if( paths_seen.count( key ) )
          add( "Two files would be '" + path + "' in the package.", true );
        paths_seen[key] = "generated:" + path;
        num_files += 1;
        total_bytes += ci.generatedData->size();
      }catch( std::exception &e )
      {
        add( e.what(), true );
      }
      continue;
    }//if( ci.generatedData )

    // A template's includes must be beside it, to come along in the package (and a packaged template
    //  may only include files in its own folder).
    if( SpecUtils::iends_with( type->id, "_template" ) )
    {
      for( const string &target : UserLibrary::templateIncludesOutsideFolder( item.found.path ) )
        add( "'" + item.relPath + "' includes '" + target + "', which is not in its folder, so would not be in the package.", true );
    }

    const UserLibrary::Validation v = UserLibrary::validateItem( item.typeId, item.found );
    if( v.status == UserLibrary::Status::Unreadable )
      add( "'" + item.relPath + "' can not be read" + (v.message.empty() ? string(".") : (": " + v.message)), true );
    else if( v.status == UserLibrary::Status::NeedsNewer )
      add( "'" + item.relPath + "' needs a newer InterSpec than this one.", false );

    for( const string &name : v.names )
    {
      const string key = UserLibrary::nameKey( *type, name );
      const map<pair<string,string>,string>::const_iterator pos = names_seen.find( { type->id, key } );
      if( pos != end(names_seen) )
        add( "Two items are named '" + name + "' ('" + pos->second + "' and '" + item.relPath + "').", true );
      else
        names_seen[{ type->id, key }] = item.relPath;
    }

    try
    {
      for( const pair<string,string> &f : packageFilesFor( ci ) )
      {
        if( UserLibrary::isSecretPath( f.second ) )
        {
          add( "'" + f.second + "' is a file InterSpec never shares.", true );
          continue;
        }

        try
        {
          checkEntryName( f.second, limits );
        }catch( std::exception &e )
        {
          add( e.what(), true );
        }

        const string key = UserLibrary::collisionKey( f.second );
        const map<string,string>::const_iterator pos = paths_seen.find( key );
        if( pos == end(paths_seen) )
        {
          paths_seen[key] = normalize( f.first );
          num_files += 1;
          total_bytes += SpecUtils::file_size( f.first );
        }else if( pos->second != normalize( f.first ) )
        {
          add( "Two files would be '" + f.second + "' in the package.", true );
        }
      }//for( loop over the item's files )
    }catch( std::exception &e )
    {
      add( e.what(), true );
    }
  }//for( loop over items )

  if( num_files > limits.maxEntries )
    add( "The package would have more files than InterSpec accepts.", true );
  if( total_bytes > limits.maxTotalBytes )
    add( "The package would be larger than InterSpec accepts.", true );

  return issues;
}//preflight(...)


UserLibrary::Manifest create( const CreateRequest &request, std::ostream &output, const Limits &limits )
{
  for( const PreflightIssue &issue : preflight( request, limits ) )
  {
    if( issue.blocking )
      throw runtime_error( issue.message );
  }

  UserLibrary::Manifest manifest = request.manifest;
  if( manifest.id.empty() )
    manifest.id = UserLibrary::newUuid();
  manifest.formatVersion = UserLibrary::Manifest::sm_formatVersion;
  manifest.created = UserLibrary::nowIsoUtc();

  string build;
#if( BUILD_AS_ELECTRON_APP )
  build = "desktop app";
#elif( BUILD_AS_OSX_APP )
  build = "macOS app";
#elif( IOS )
  build = "iOS app";
#elif( ANDROID )
  build = "Android app";
#elif( BUILD_AS_WX_WIDGETS_APP )
  build = "desktop app";
#else
  build = "local server";
#endif
  manifest.createdWith = "InterSpec " + UserLibrary::interSpecVersion() + " (" + build + ")";

  // Earlier versions of this package, as this user made them: their change notes, and item ids.
  map<string,string> previous_ids;  // source key -> item_id
  const map<string,vector<UserLibrary::CreatedVersion>> created = UserLibrary::createdPackages();
  const map<string,vector<UserLibrary::CreatedVersion>>::const_iterator prev = created.find( manifest.id );
  if( (prev != end(created)) && !prev->second.empty() )
  {
    const UserLibrary::CreatedVersion &last = prev->second.back();
    if( manifest.previousVersions.empty() )
    {
      manifest.previousVersions.push_back( { last.manifest.version, last.manifest.created, last.manifest.changes } );
      manifest.previousVersions.insert( end(manifest.previousVersions),
                                        begin(last.manifest.previousVersions), end(last.manifest.previousVersions) );
    }
    for( const pair<const string,string> &s : last.sources )
      previous_ids[s.second] = s.first;
  }//if( made before )

  manifest.items.clear();
  map<string,string> sources;  // item_id -> source key
  set<string> written;
  ZipArchive::ZipWriter zip( output );

  for( const CreateItem &ci : request.items )
  {
    const UserLibrary::Item &item = ci.item;

    UserLibrary::ManifestItem mi;
    const string key = sourceKey( ci );
    const map<string,string>::const_iterator id_pos = previous_ids.find( key );
    mi.itemId = (id_pos != end(previous_ids)) ? id_pos->second : UserLibrary::newUuid();
    mi.path = packagePathFor( ci );
    mi.type = item.typeId;
    mi.notes = ci.notes;

    // A detector efficiency from the database goes in as the XML InterSpec exports.
    if( ci.databaseDrf )
    {
      const string xml = drfXml( *ci.databaseDrf );
      mi.sha256 = UserLibrary::sha256OfData( xml.data(), xml.size() );
      mi.intendedOverride = ci.intendedOverride;
      if( written.insert( mi.path ).second )
        zip.add_file( mi.path, xml );
      sources[mi.itemId] = key;
      manifest.items.push_back( mi );
      continue;
    }//if( ci.databaseDrf )

    if( ci.generatedData )
    {
      mi.sha256 = UserLibrary::sha256OfData( ci.generatedData->data(), ci.generatedData->size() );
      mi.intendedOverride = ci.intendedOverride;
      if( written.insert( mi.path ).second )
        zip.add_file( mi.path, *ci.generatedData );
      sources[mi.itemId] = key;
      manifest.items.push_back( mi );
      continue;
    }//if( ci.generatedData )

    mi.sha256 = UserLibrary::sha256OfItem( item.found );
    mi.intendedOverride = ci.intendedOverride;

    // Re-sharing an item from someone else's package: say where it came from.
    if( item.source.kind == UserLibrary::SourceKind::Package )
    {
      const optional<UserLibrary::Package> pkg = UserLibrary::installedPackage( item.source.packageId );
      const UserLibrary::ManifestItem * const orig = pkg ? pkg->manifest.itemForPath( item.relPath ) : nullptr;
      if( orig && !orig->derivedFromId.empty() )
      {
        mi.derivedFromName = orig->derivedFromName;
        mi.derivedFromId = orig->derivedFromId;
        mi.derivedFromVersion = orig->derivedFromVersion;
      }else if( pkg )
      {
        mi.derivedFromName = pkg->manifest.name;
        mi.derivedFromId = pkg->manifest.id;
        mi.derivedFromVersion = pkg->manifest.version;
      }
    }//if( from a package )

    for( const pair<string,string> &f : packageFilesFor( ci ) )
    {
      if( !written.insert( f.second ).second )
        continue;  // e.g. a DETECTOR.txt two .par grids share
      zip.add_file( f.second, read_whole_file( f.first ) );
    }

    sources[mi.itemId] = key;
    manifest.items.push_back( mi );
  }//for( loop over items )

  zip.add_file( UserLibrary::ManifestFileName, manifest.toJson() );
  zip.finish();

  UserLibrary::modifyHistory( [&]( json &history ){
    json &c = history["created"];
    if( !c.is_object() )
      c = json::object();
    json &versions = c[manifest.id];
    if( !versions.is_array() )
      versions = json::array();
    versions.push_back( { {"created", manifest.created}, {"manifest", manifest.toJsonObject()}, {"sources", sources} } );
    while( versions.size() > 20 )  //the newest is what the next version starts from
      versions.erase( versions.begin() );

    history["last_creator"] = { {"author", manifest.author}, {"organization", manifest.organization},
                                {"contact", manifest.contact} };
  } );

  return manifest;
}//create(...)


std::string nextVersion( const std::string &version )
{
  const string trimmed = SpecUtils::trim_copy( version );
  vector<string> parts;
  SpecUtils::split( parts, trimmed, "." );
  bool numeric = !parts.empty();
  for( const string &p : parts )
    numeric = numeric && !p.empty() && (p.size() < 10)
              && std::all_of( begin(p), end(p), []( char c ){ return (c >= '0') && (c <= '9'); } );
  if( !numeric )
    return trimmed.empty() ? string( "1" ) : (trimmed + ".1");

  parts.back() = std::to_string( std::stoul( parts.back() ) + 1 );
  string answer;
  for( const string &p : parts )
    answer += (answer.empty() ? "" : ".") + p;
  return answer;
}//nextVersion(...)
}//namespace UserLibraryPackage
