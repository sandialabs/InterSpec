#ifndef ZIP_ARCHIVE_READER_H
#define ZIP_ARCHIVE_READER_H
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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
    
    bool init( std::istream &istream, const bool global );
  };//struct ZipFileHeader
  
  
  //openZipFile(): throws exception with descriptive message upon error, or
  //  no files found.
  typedef std::map<std::string, std::shared_ptr<const ZipFileHeader> > \
          FilenameToZipHeaderMap;
  
  //open_zip_file(): throws std::exception with descriptive message upon error.
  //  Output is garunteed to have at least one entry.
  FilenameToZipHeaderMap open_zip_file( std::istream &instrm );

  //read_file_from_zip():  throws std::exception upon error; in which case it
  //  is indeterminant if any data has been written to the output stream.
  //  Returns number of bytes written to output, which is garunteed to be
  //  greater than zero.
  size_t read_file_from_zip( std::istream &instrm,
                             std::shared_ptr<const ZipFileHeader> header,
                             std::ostream &output );
}//namespace ZipArchive
#endif
