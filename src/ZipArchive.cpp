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

#include <string>
#include <memory>
#include <iostream>
#include <algorithm>
#include <stdexcept>

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

template<class T>
inline istream &binaryRead( istream &stream, T &x )
{
  return stream.read( (char *)(&x), sizeof(T) );
}


bool ZipFileHeader::init( istream& istream, const bool globalHeader )
{
  uint32_t sig;
  uint16_t version, flags;
  
  if( globalHeader )
  {
    binaryRead(istream,sig);
    if( sig != 0x02014b50 )
    {
      cerr<<"Did not find global header signature"<<endl;
      return false;
    }
    
    binaryRead(istream,version);
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
  }//if( globalHeader )
  
  uint16_t bufflen = std::max( filename_length, extra_length );
  bufflen = std::max( bufflen, comment_length ) + 1;
      
  vector<char> buff( bufflen );
  
  istream.read( &buff[0], filename_length );
  buff[filename_length] = '\0';
  filename = string( &buff[0] );
  
  istream.read( &buff[0],extra_length);
  if( globalHeader )
    istream.read( &buff[0],comment_length);
  
  return true;
}//bool init( istream& istream, const bool globalHeader )


size_t read_file_from_zip( std::istream &instrm,
                           std::shared_ptr<const ZipFileHeader> header,
                           std::ostream &output )
{
  if( !header )
    throw runtime_error( "ZipArchive: no zip file header passed in to read" );
  
  instrm.seekg( header->header_offset );
  
  unsigned int total_read = 0;
  int total_uncompressed = 0;
  const unsigned int buffer_size = 512;
  unsigned char in[buffer_size], out[buffer_size];
  const uint16_t DEFLATE = 8;
  const uint16_t UNCOMPRESSED = 0;
    
  z_stream strm;
    
  strm.zalloc   = Z_NULL;
  strm.zfree    = Z_NULL;
  strm.opaque   = Z_NULL;
  strm.avail_in = 0;
  strm.next_in  = Z_NULL;
  

  ZipFileHeader localheader;
  if( !localheader.init( instrm, false ) )
    throw runtime_error( "ZipArchive: error reading header from stream" );
  
  
  // initialize the inflate
  if( header->compression_type == DEFLATE ) //should maybye use localheader
  {
    int result = inflateInit2( &strm, -MAX_WBITS );
    if( result != Z_OK )
      throw runtime_error( "ZipArchive: gzip inflateInit2 didnt return Z_OK" );
    
    bool end_of_file = false;
    while( !end_of_file )
    {
      strm.avail_out = buffer_size-4;
      strm.next_out = (Bytef*)(out+4);
      
      while( strm.avail_out != 0 )
      {
        if( strm.avail_in == 0 ) // buffer empty, read some more from file
        {
          if( header->compressed_size <= total_read )
            throw runtime_error( "ZipArchive: read size error" );
          
          const unsigned int nToRead = std::min( buffer_size,
                                        header->compressed_size - total_read );
             
          instrm.read( (char*)in, nToRead );
          strm.avail_in = static_cast<unsigned int>( instrm.gcount() );
          total_read += strm.avail_in;
          strm.next_in = (Bytef*)in;
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
          default:
            break;
        }//switch( ret )
        
        if( ret == Z_STREAM_END )
        {
          end_of_file = true;
          break;
        }
      }//while( strm.avail_out != 0 )
      
      const int unzip_count = buffer_size - strm.avail_out - 4;
      
      output.write( (const char *) (out + 4), unzip_count );
      
      total_uncompressed += unzip_count;
    }
    
    inflateEnd( &strm );
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
      output.write( (const char *)in, bytes_read );
      num_bytes -= std::min( bytes_read, static_cast<std::streamsize>(num_bytes) );
      num_written += bytes_read;

      if( bytes_read < static_cast<std::streamsize>(nread) )
        break;
    }//while( num_bytes > 0 )
    
    return num_written;
  }else
  {
    throw runtime_error( "ZipArchive: unrecognized compression" );
  }//if( header->compression_type == DEFLATE ) / else UNCOMPRESSED / else
  
}//size_t read_file_from_zip(...)



std::map<std::string, std::shared_ptr<const ZipFileHeader> >
                                          open_zip_file( std::istream &instrm )
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
  map<std::string, std::shared_ptr<const ZipFileHeader> > answer;
  instrm.seekg( header_offset );
  for( uint16_t i = 0; i < num_files; ++i )
  {
    std::shared_ptr<ZipFileHeader> header( new ZipFileHeader );
    const bool valid = header->init( instrm, true );
    if( valid )
      answer[header->filename] = header;
  }
  
  if( answer.empty() )
    throw runtime_error( "ZipArchive: no files found in archive" );
  
  return answer;
}
}//namespace ZipArchive

