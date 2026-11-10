#ifndef ZIP_ARCHIVE_READER_H
#define ZIP_ARCHIVE_READER_H
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

#include <map>
#include <cctype>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <stdint.h>

/** ZipArchive opens a ZIP file and allows you to extract files it contains.
   Its not incredibly well tested, and could definetly stand to use more error
   checking, but it does seem to work.
 
   To use this code, you must link to zlib as well as have zlib.h in your 
   include path.
*/

namespace ZipArchive
{
  struct ExtractionLimits
  {
    size_t max_uncompressed_bytes;
    size_t max_compression_ratio;
  };

  struct ZipFileHeader
  {
    uint16_t version;
    uint16_t flags;
    uint16_t compression_type;
    uint16_t stamp_date,stamp_time;
    uint32_t crc;
    uint32_t compressed_size;
    uint32_t uncompressed_size;
    std::string filename;
    uint32_t header_offset; // local header offset

    // Central-directory only (zero for a local header):
    uint16_t version_made_by = 0;     // high byte is the host system (3 = Unix)
    uint32_t external_attributes = 0; // for Unix hosts, the high 16 bits are the st_mode

    bool init( std::istream &istream, const bool global );

    /** Whether the entry is a symbolic link (only Unix-made archives record this). */
    bool is_symlink() const;

    /** Whether the entry is a directory (its name ends with '/'). */
    bool is_directory() const;
  };//struct ZipFileHeader
  
  
  //openZipFile(): throws exception with descriptive message upon error, or
  //  no files found.
  typedef std::map<std::string, std::shared_ptr<const ZipFileHeader> > \
          FilenameToZipHeaderMap;
  
  //open_zip_file(): throws std::exception with descriptive message upon error.
  //  Output is garunteed to have at least one entry.
  FilenameToZipHeaderMap open_zip_file( std::istream &instrm );

  /** Every entry of the central directory, in order - unlike #open_zip_file, entries with the same
   name are all kept, so a caller can reject such an archive.  Throws std::exception on error.
   */
  std::vector<std::shared_ptr<const ZipFileHeader>> list_zip_entries( std::istream &instrm );

  // Throws if any member exceeds the per-file limits, if declared aggregate
  // output exceeds max_aggregate_uncompressed_bytes, or if summing sizes would
  // overflow.
  void validate_archive_for_extraction(
                         const FilenameToZipHeaderMap &headers,
                         const ExtractionLimits &limits,
                         size_t max_aggregate_uncompressed_bytes );

  //read_file_from_zip():  throws std::exception upon error; in which case it
  //  is indeterminant if any data has been written to the output stream.
  //  Returns number of bytes written to output, which is garunteed to be
  //  greater than zero.
  size_t read_file_from_zip( std::istream &instrm,
                             std::shared_ptr<const ZipFileHeader> header,
                             std::ostream &output,
                             const ExtractionLimits &limits );

  /** Writes a zip archive (deflate, no Zip64: under 65535 entries and 4 GB).  Names are stored as
   UTF-8, '/' separated.

   Usage: `ZipWriter w( out ); w.add_file( "a/b.txt", data ); w.finish();`
   */
  class ZipWriter
  {
  public:
    explicit ZipWriter( std::ostream &output );

    /** Adds a file; throws on error (e.g. a write failure, or the archive getting too large). */
    void add_file( const std::string &name, const std::string &data, const bool compress = true );

    /** Writes the central directory; the archive is incomplete until called.  Throws on error. */
    void finish();

  private:
    struct Entry
    {
      std::string name;
      uint32_t crc = 0, compressed_size = 0, uncompressed_size = 0, offset = 0;
      uint16_t method = 0, dos_time = 0, dos_date = 0;
    };

    std::ostream &m_output;
    std::vector<Entry> m_entries;
    uint64_t m_position;
    bool m_finished;
  };//class ZipWriter
}//namespace ZipArchive
#endif
