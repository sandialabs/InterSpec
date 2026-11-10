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

extern "C"{
#include <zlib.h>
}

#include <ctime>
#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <limits>

#include "InterSpec/ZipArchive.h"

using namespace std;

//Code from this file was heavily influenced by
//  https://github.com/wdas/partio/blob/master/src/lib/io/ (Author Andrew Selle)
//  in terms of how to actually read from a zip archive.  The actual interface
//  here is totally differnent, as well as many implementation aspects, so I'm
//  assuming no formal attribution or licensing applies since it boils down to
//  zlib example code.


namespace ZipArchive
{

namespace
{
size_t validated_member_output_limit( const ZipFileHeader &header,
                                      const ExtractionLimits &limits )
{
  if( (limits.max_uncompressed_bytes == 0) || (limits.max_compression_ratio == 0) )
    throw runtime_error( "ZipArchive: invalid extraction limits" );

  const size_t uncompressed = header.uncompressed_size;
  const size_t compressed = header.compressed_size;
  if( uncompressed > limits.max_uncompressed_bytes )
    throw runtime_error( "ZipArchive: member exceeds uncompressed size limit" );

  if( uncompressed == 0 )
    return 0;

  if( compressed == 0 )
    throw runtime_error( "ZipArchive: non-empty member has zero compressed size" );

  const size_t quotient = uncompressed / compressed;
  const size_t remainder = uncompressed % compressed;
  if( (quotient > limits.max_compression_ratio)
     || ((quotient == limits.max_compression_ratio) && (remainder != 0)) )
    throw runtime_error( "ZipArchive: member exceeds compression ratio limit" );

  size_t ratio_limit = limits.max_uncompressed_bytes;
  if( compressed <= (std::numeric_limits<size_t>::max() / limits.max_compression_ratio) )
    ratio_limit = std::min( ratio_limit, compressed * limits.max_compression_ratio );

  return ratio_limit;
}
}//namespace


void validate_archive_for_extraction( const FilenameToZipHeaderMap &headers,
                                      const ExtractionLimits &limits,
                                      size_t max_aggregate_uncompressed_bytes )
{
  if( max_aggregate_uncompressed_bytes == 0 )
    throw runtime_error( "ZipArchive: invalid aggregate extraction limit" );

  size_t aggregate = 0;
  for( const auto &entry : headers )
  {
    if( !entry.second )
      throw runtime_error( "ZipArchive: archive contains a null member header" );

    validated_member_output_limit( *entry.second, limits );
    const size_t member_size = entry.second->uncompressed_size;
    if( member_size > (max_aggregate_uncompressed_bytes - aggregate) )
      throw runtime_error( "ZipArchive: archive exceeds aggregate extraction limit" );
    aggregate += member_size;
  }
}//validate_archive_for_extraction(...)

template<class T>
inline istream &binaryRead( istream &stream, T &x )
{
  return stream.read( (char *)(&x), sizeof(T) );
}


bool ZipFileHeader::init( istream& istream, const bool globalHeader )
{
  uint32_t sig;
  uint16_t version, flags;
  
  version_made_by = 0;
  external_attributes = 0;

  if( globalHeader )
  {
    binaryRead(istream,sig);
    if( sig != 0x02014b50 )
    {
      cerr<<"Did not find global header signature"<<endl;
      return false;
    }
    
    binaryRead( istream, version_made_by );
  }else
  {
    binaryRead(istream,sig);
    if( sig != 0x04034b50 )
    {
      cerr << "Did not find local header signature" << endl;
      return false;
    }
  }//if( globalHeader ) / else
  
  // Read rest of header
  binaryRead( istream, version );
  binaryRead( istream, flags );
  binaryRead( istream, compression_type );
  binaryRead( istream, stamp_date );
  binaryRead( istream, stamp_time );
  binaryRead( istream, crc );
  binaryRead( istream, compressed_size );
  binaryRead( istream, uncompressed_size );
  
  uint16_t filename_length,extra_length;
  binaryRead( istream, filename_length );
  binaryRead( istream, extra_length );
  
  uint16_t comment_length = 0;
  if( globalHeader )
  {
    binaryRead(istream,comment_length); // filecomment
    unsigned short disk_number_start,int_file_attrib;
    unsigned int ext_file_attrib;
    binaryRead( istream, disk_number_start ); // disk# start
    binaryRead( istream, int_file_attrib ); // internal file
    binaryRead( istream, ext_file_attrib ); // ext final
    binaryRead( istream, header_offset ); // rel offset
    external_attributes = ext_file_attrib;
  }//if( globalHeader )
  
  // A truncated or corrupt archive: the fields above may not have been read.
  if( !istream )
    return false;

  // size_t, so lengths of 0xFFFF don't wrap the buffer size to zero.
  const size_t bufflen = static_cast<size_t>( std::max( std::max( filename_length, extra_length ), comment_length ) ) + 1;
  vector<char> buff( bufflen, '\0' );
  
  istream.read( &buff[0], filename_length );
  if( !istream )
    return false;
  filename = string( &buff[0], &buff[0] + filename_length );
  const size_t nul_pos = filename.find( '\0' );
  if( nul_pos != string::npos )
    filename.resize( nul_pos );
  
  istream.read( &buff[0],extra_length);
  if( globalHeader )
    istream.read( &buff[0],comment_length);
  
  return !!istream;
}//bool init( istream& istream, const bool globalHeader )


size_t read_file_from_zip( std::istream &instrm,
                           std::shared_ptr<const ZipFileHeader> header,
                           std::ostream &output,
                           const ExtractionLimits &limits )
{
  if( !header )
    throw runtime_error( "ZipArchive: no zip file header passed in to read" );

  const size_t output_limit = validated_member_output_limit( *header, limits );
  
  instrm.seekg( header->header_offset );
  if( !instrm )
    throw runtime_error( "ZipArchive: failed seeking to member header" );
  
  size_t total_read = 0;
  size_t total_uncompressed = 0;
  const unsigned int buffer_size = 512;
  unsigned char in[buffer_size], out[buffer_size];
  const uint16_t DEFLATE = 8;
  const uint16_t UNCOMPRESSED = 0;
    
  z_stream strm = {};
  

  ZipFileHeader localheader;
  if( !localheader.init( instrm, false ) )
    throw runtime_error( "ZipArchive: error reading header from stream" );
  
  
  // initialize the inflate
  if( header->compression_type == DEFLATE ) //should maybye use localheader
  {
    const int result = inflateInit2( &strm, -MAX_WBITS );
    if( result != Z_OK )
      throw runtime_error( "ZipArchive: gzip inflateInit2 didnt return Z_OK" );

    struct InflateGuard
    {
      z_stream *stream;
      ~InflateGuard() { inflateEnd( stream ); }
    } guard{ &strm };
    
    bool end_of_file = false;
    while( !end_of_file )
    {
      strm.avail_out = buffer_size;
      strm.next_out = reinterpret_cast<Bytef *>(out);
      
      while( strm.avail_out != 0 )
      {
        if( strm.avail_in == 0 ) // buffer empty, read some more from file
        {
          if( header->compressed_size <= total_read )
            throw runtime_error( "ZipArchive: read size error" );
          
          const size_t nToRead = std::min<size_t>( buffer_size,
                                                   header->compressed_size - total_read );
             
          instrm.read( reinterpret_cast<char *>(in), static_cast<std::streamsize>(nToRead) );
          const std::streamsize count = instrm.gcount();
          if( count <= 0 )
            throw runtime_error( "ZipArchive: truncated compressed member" );

          strm.avail_in = static_cast<unsigned int>(count);
          total_read += strm.avail_in;
          strm.next_in = reinterpret_cast<Bytef *>(in);
        }//if( strm.avail_in == 0 )
        
        const int ret = inflate( &strm, Z_NO_FLUSH ); // decompress
        switch( ret )
        {
          case Z_STREAM_ERROR:
            throw runtime_error( "ZipArchive: libz error Z_STREAM_ERROR" );
          
          case Z_NEED_DICT:
          case Z_DATA_ERROR:
          case Z_MEM_ERROR:
            throw runtime_error( "ZipArchive: gzip error "
                                 + string(strm.msg ? strm.msg : "(no msg)") );
          case Z_BUF_ERROR:
            throw runtime_error( "ZipArchive: gzip made no decompression progress" );
          default:
            break;
        }//switch( ret )
        
        if( ret == Z_STREAM_END )
        {
          end_of_file = true;
          break;
        }
      }//while( strm.avail_out != 0 )
      
      const size_t unzip_count = buffer_size - strm.avail_out;
      if( unzip_count > (output_limit - total_uncompressed) )
        throw runtime_error( "ZipArchive: inflated member exceeds extraction limit" );

      if( unzip_count )
      {
        output.write( reinterpret_cast<const char *>(out),
                      static_cast<std::streamsize>(unzip_count) );
        if( !output )
          throw runtime_error( "ZipArchive: failed writing extracted member" );
        total_uncompressed += unzip_count;
      }
    }

    if( total_uncompressed != header->uncompressed_size )
      throw runtime_error( "ZipArchive: uncompressed size does not match archive header" );
    return total_uncompressed;
  }else if( (header->compression_type == UNCOMPRESSED)
           && (header->compressed_size == header->uncompressed_size) )
  {
    size_t num_written = 0;
    size_t num_bytes = header->compressed_size;
    while( num_bytes > 0 )
    {
      const size_t nread = std::min( buffer_size, static_cast<const unsigned int>(num_bytes) );
      instrm.read( (char *)in, nread );
      std::streamsize bytes_read = instrm.gcount();  // Get the number of bytes actually read
      if( bytes_read <= 0 )
        throw runtime_error( "ZipArchive: truncated uncompressed member" );

      if( static_cast<size_t>(bytes_read) > (output_limit - num_written) )
        throw runtime_error( "ZipArchive: stored member exceeds extraction limit" );

      output.write( (const char *)in, bytes_read );
      if( !output )
        throw runtime_error( "ZipArchive: failed writing extracted member" );
      num_bytes -= std::min( bytes_read, static_cast<std::streamsize>(num_bytes) );
      num_written += bytes_read;

      if( bytes_read < static_cast<std::streamsize>(nread) )
        throw runtime_error( "ZipArchive: truncated uncompressed member" );
    }//while( num_bytes > 0 )

    if( num_written != header->uncompressed_size )
      throw runtime_error( "ZipArchive: uncompressed size does not match archive header" );
    
    return num_written;
  }else
  {
    throw runtime_error( "ZipArchive: unrecognized compression" );
  }//if( header->compression_type == DEFLATE ) / else UNCOMPRESSED / else
  
}//size_t read_file_from_zip(...)



bool ZipFileHeader::is_symlink() const
{
  // Unix hosts keep st_mode in the high 16 bits; S_IFLNK is 0120000.
  const uint8_t host = static_cast<uint8_t>( version_made_by >> 8 );
  return (host == 3) && (((external_attributes >> 16) & 0170000) == 0120000);
}


bool ZipFileHeader::is_directory() const
{
  return !filename.empty() && ((filename.back() == '/') || (filename.back() == '\\'));
}


std::map<std::string, std::shared_ptr<const ZipFileHeader> >
                                          open_zip_file( std::istream &instrm )
{
  map<std::string, std::shared_ptr<const ZipFileHeader> > answer;
  for( const std::shared_ptr<const ZipFileHeader> &header : list_zip_entries( instrm ) )
    answer[header->filename] = header;
  return answer;
}//open_zip_file(...)


std::vector<std::shared_ptr<const ZipFileHeader>> list_zip_entries( std::istream &instrm )
{
  //Assumes the header is at the end of the file
  instrm.seekg( 0, ios_base::end );
  const ios::pos_type end_position = instrm.tellg();
  const size_t max_comment_size = 0xffff; // max size of header, 65535
  const size_t read_size_before_comment = 22;
  std::streamoff read_start = max_comment_size + read_size_before_comment;
  
  if( read_start > end_position)
    read_start = end_position;

  if( read_start < 0 )
    throw runtime_error( "ZipArchive: invalid file size" );

  instrm.seekg( end_position - read_start );
  
  //boost::shared_array<char> buf( new char[static_cast<size_t>(read_start)] );
  auto buf = std::unique_ptr<char []>{ new char[static_cast<size_t>(read_start)] };
  
  instrm.read( buf.get(),read_start);
  
  size_t headerstart = -1;
  for( size_t i = 0; static_cast<std::streamoff>(i) < read_start-3; i++)
  {
    if( buf[i]==0x50 && buf[i+1]==0x4b && buf[i+2]==0x05 && buf[i+3]==0x06 )
    {
      headerstart = i;
      break;
    }
  }
  buf.reset();
  
  if( headerstart == -1 )
    throw runtime_error( "ZipArchive: Couldnt find zip header" );
  
  const std::streamoff nbytesFromEnd = read_start - headerstart;
  instrm.seekg( end_position - nbytesFromEnd );
  
  uint32_t word;
  uint16_t this_disk_num, end_disk_num, num_files, num_files_this_disk;
  binaryRead( instrm, word ); // end of central
  binaryRead( instrm, this_disk_num ); // this disk number
  binaryRead( instrm, end_disk_num ); // this disk number
  
  if( this_disk_num != end_disk_num || this_disk_num != 0 )
    throw runtime_error( "ZipArchive: multi-disk zip files not supported" );
  
  binaryRead( instrm, num_files );
  binaryRead( instrm, num_files_this_disk );
  
  if( num_files != num_files_this_disk )
    throw runtime_error( "ZipArchive: multi-disk zip files not supported" );
  
  uint32_t header_size, header_offset;
  binaryRead( instrm, header_size ); // size of header
  binaryRead( instrm, header_offset ); // offset to header
  
  //lets read all the file headers
  std::vector<std::shared_ptr<const ZipFileHeader>> answer;
  instrm.seekg( header_offset );
  for( uint16_t i = 0; i < num_files; ++i )
  {
    std::shared_ptr<ZipFileHeader> header( new ZipFileHeader );
    // A bad entry leaves the stream out of step with the rest: the listing can't be trusted.
    if( !header->init( instrm, true ) )
      throw runtime_error( "ZipArchive: corrupt central directory (entry " + std::to_string( i + 1 )
                           + " of " + std::to_string( num_files ) + ")" );
    answer.push_back( header );
  }
  
  if( answer.empty() )
    throw runtime_error( "ZipArchive: no files found in archive" );
  
  return answer;
}//list_zip_entries(...)


namespace
{
  template<class T>
  void binaryWrite( std::ostream &stream, const T x )
  {
    stream.write( reinterpret_cast<const char *>( &x ), sizeof(T) );
  }

  /** MS-DOS date and time of now (local time), as zip files store them. */
  void dos_date_time( uint16_t &dos_time, uint16_t &dos_date )
  {
    const std::time_t t = std::time( nullptr );
    std::tm tm{};
#ifdef _WIN32
    localtime_s( &tm, &t );
#else
    localtime_r( &t, &tm );
#endif
    const int year = std::max( 1980, tm.tm_year + 1900 );
    dos_time = static_cast<uint16_t>( (tm.tm_hour << 11) | (tm.tm_min << 5) | (tm.tm_sec / 2) );
    dos_date = static_cast<uint16_t>( ((year - 1980) << 9) | ((tm.tm_mon + 1) << 5) | tm.tm_mday );
  }//dos_date_time(...)
}//namespace


ZipWriter::ZipWriter( std::ostream &output )
  : m_output( output ),
    m_entries(),
    m_position( 0 ),
    m_finished( false )
{
}


void ZipWriter::add_file( const std::string &name, const std::string &data, const bool compress )
{
  if( m_finished )
    throw runtime_error( "ZipWriter: archive already finished" );
  if( name.empty() || (name.size() > 0xFFFF) )
    throw runtime_error( "ZipWriter: invalid file name" );
  if( m_entries.size() >= 0xFFFF )
    throw runtime_error( "ZipWriter: too many files" );
  if( data.size() >= 0xFFFFFFFFull )
    throw runtime_error( "ZipWriter: file too large" );

  Entry entry;
  entry.name = name;
  entry.uncompressed_size = static_cast<uint32_t>( data.size() );
  entry.crc = static_cast<uint32_t>( crc32( 0L, reinterpret_cast<const Bytef *>( data.data() ),
                                            static_cast<uInt>( data.size() ) ) );
  dos_date_time( entry.dos_time, entry.dos_date );

  string compressed;
  entry.method = 0;
  if( compress && !data.empty() )
  {
    z_stream strm = {};
    if( deflateInit2( &strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -MAX_WBITS, 8, Z_DEFAULT_STRATEGY ) != Z_OK )
      throw runtime_error( "ZipWriter: deflateInit2 failed" );

    compressed.resize( deflateBound( &strm, static_cast<uLong>( data.size() ) ) );
    strm.next_in = reinterpret_cast<Bytef *>( const_cast<char *>( data.data() ) );
    strm.avail_in = static_cast<uInt>( data.size() );
    strm.next_out = reinterpret_cast<Bytef *>( &compressed[0] );
    strm.avail_out = static_cast<uInt>( compressed.size() );
    const int result = deflate( &strm, Z_FINISH );
    const uLong total_out = strm.total_out;
    deflateEnd( &strm );
    if( result != Z_STREAM_END )
      throw runtime_error( "ZipWriter: deflate failed" );
    compressed.resize( total_out );

    // Stored is better for data that does not compress (images, already-compressed files).  Data that
    //  compresses more than 100:1 is stored too, so readers that refuse high ratios as possible zip
    //  bombs (InterSpec refuses over 200:1) still accept it.
    if( (compressed.size() < data.size()) && ((compressed.size() * 100) >= data.size()) )
      entry.method = 8;
  }//if( compress )

  const string &payload = (entry.method == 8) ? compressed : data;
  entry.compressed_size = static_cast<uint32_t>( payload.size() );

  if( (m_position + 30 + name.size() + payload.size()) >= 0xFFFFFFFFull )
    throw runtime_error( "ZipWriter: archive too large" );
  entry.offset = static_cast<uint32_t>( m_position );

  // Local file header; general purpose bit 11: the name is UTF-8.
  binaryWrite<uint32_t>( m_output, 0x04034b50 );
  binaryWrite<uint16_t>( m_output, 20 );            // version needed to extract
  binaryWrite<uint16_t>( m_output, 0x0800 );        // flags
  binaryWrite<uint16_t>( m_output, entry.method );
  binaryWrite<uint16_t>( m_output, entry.dos_time );
  binaryWrite<uint16_t>( m_output, entry.dos_date );
  binaryWrite<uint32_t>( m_output, entry.crc );
  binaryWrite<uint32_t>( m_output, entry.compressed_size );
  binaryWrite<uint32_t>( m_output, entry.uncompressed_size );
  binaryWrite<uint16_t>( m_output, static_cast<uint16_t>( name.size() ) );
  binaryWrite<uint16_t>( m_output, 0 );             // extra field length
  m_output.write( name.data(), name.size() );
  m_output.write( payload.data(), payload.size() );
  if( !m_output )
    throw runtime_error( "ZipWriter: write failed" );

  m_position += 30 + name.size() + payload.size();
  m_entries.push_back( entry );
}//ZipWriter::add_file(...)


void ZipWriter::finish()
{
  if( m_finished )
    return;

  const uint64_t cd_start = m_position;
  for( const Entry &entry : m_entries )
  {
    binaryWrite<uint32_t>( m_output, 0x02014b50 );
    binaryWrite<uint16_t>( m_output, 0x031E );        // made by: Unix, zip spec 3.0
    binaryWrite<uint16_t>( m_output, 20 );            // version needed to extract
    binaryWrite<uint16_t>( m_output, 0x0800 );        // flags: UTF-8 name
    binaryWrite<uint16_t>( m_output, entry.method );
    binaryWrite<uint16_t>( m_output, entry.dos_time );
    binaryWrite<uint16_t>( m_output, entry.dos_date );
    binaryWrite<uint32_t>( m_output, entry.crc );
    binaryWrite<uint32_t>( m_output, entry.compressed_size );
    binaryWrite<uint32_t>( m_output, entry.uncompressed_size );
    binaryWrite<uint16_t>( m_output, static_cast<uint16_t>( entry.name.size() ) );
    binaryWrite<uint16_t>( m_output, 0 );             // extra field length
    binaryWrite<uint16_t>( m_output, 0 );             // comment length
    binaryWrite<uint16_t>( m_output, 0 );             // disk number start
    binaryWrite<uint16_t>( m_output, 0 );             // internal attributes
    binaryWrite<uint32_t>( m_output, 0100644u << 16 ); // external: a regular file, rw-r--r--
    binaryWrite<uint32_t>( m_output, entry.offset );
    m_output.write( entry.name.data(), entry.name.size() );
    m_position += 46 + entry.name.size();
  }//for( loop over entries )

  if( m_position >= 0xFFFFFFFFull )
    throw runtime_error( "ZipWriter: archive too large" );

  const uint64_t cd_size = m_position - cd_start;
  binaryWrite<uint32_t>( m_output, 0x06054b50 );
  binaryWrite<uint16_t>( m_output, 0 );               // this disk
  binaryWrite<uint16_t>( m_output, 0 );               // disk with the central directory
  binaryWrite<uint16_t>( m_output, static_cast<uint16_t>( m_entries.size() ) );
  binaryWrite<uint16_t>( m_output, static_cast<uint16_t>( m_entries.size() ) );
  binaryWrite<uint32_t>( m_output, static_cast<uint32_t>( cd_size ) );
  binaryWrite<uint32_t>( m_output, static_cast<uint32_t>( cd_start ) );
  binaryWrite<uint16_t>( m_output, 0 );               // comment length
  m_output.flush();
  if( !m_output )
    throw runtime_error( "ZipWriter: write failed" );

  m_finished = true;
}//ZipWriter::finish()
}//namespace ZipArchive
