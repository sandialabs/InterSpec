/* Security regression tests for bounded ZIP extraction. */

#define BOOST_TEST_MODULE test_ZipArchive_suite
#include <boost/test/included/unit_test.hpp>

#include <cstdint>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include <zlib.h>
}

#include "InterSpec/ZipArchive.h"

using namespace std;

namespace
{
template<typename T>
void append_le( string &bytes, T value )
{
  for( size_t i = 0; i < sizeof(T); ++i )
    bytes.push_back( static_cast<char>((value >> (8u * i)) & 0xffu) );
}

string local_member( const string &payload, uint16_t compression,
                     uint32_t compressed_size, uint32_t uncompressed_size )
{
  string bytes;
  append_le<uint32_t>( bytes, 0x04034b50u );
  append_le<uint16_t>( bytes, 20u );
  append_le<uint16_t>( bytes, 0u );
  append_le<uint16_t>( bytes, compression );
  append_le<uint16_t>( bytes, 0u );
  append_le<uint16_t>( bytes, 0u );
  append_le<uint32_t>( bytes, 0u );
  append_le<uint32_t>( bytes, compressed_size );
  append_le<uint32_t>( bytes, uncompressed_size );
  append_le<uint16_t>( bytes, 1u );
  append_le<uint16_t>( bytes, 0u );
  bytes.push_back( 'x' );
  bytes += payload;
  return bytes;
}

string raw_deflate( const string &input )
{
  z_stream stream = {};
  BOOST_REQUIRE_EQUAL( deflateInit2(&stream, Z_BEST_COMPRESSION, Z_DEFLATED,
                                    -MAX_WBITS, 8, Z_DEFAULT_STRATEGY), Z_OK );

  vector<unsigned char> output( compressBound(static_cast<uLong>(input.size())) );
  stream.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(input.data()));
  stream.avail_in = static_cast<uInt>(input.size());
  stream.next_out = output.data();
  stream.avail_out = static_cast<uInt>(output.size());
  BOOST_REQUIRE_EQUAL( deflate(&stream, Z_FINISH), Z_STREAM_END );
  const size_t produced = output.size() - stream.avail_out;
  deflateEnd( &stream );
  return string( reinterpret_cast<const char *>(output.data()), produced );
}

shared_ptr<ZipArchive::ZipFileHeader> header_for( uint16_t compression,
                                                  uint32_t compressed,
                                                  uint32_t uncompressed )
{
  auto header = make_shared<ZipArchive::ZipFileHeader>();
  header->compression_type = compression;
  header->compressed_size = compressed;
  header->uncompressed_size = uncompressed;
  header->header_offset = 0;
  header->filename = "x";
  return header;
}
}//namespace


BOOST_AUTO_TEST_CASE( stored_member_respects_limits )
{
  const string payload = "spectrum data";
  const string archive = local_member( payload, 0u, payload.size(), payload.size() );
  istringstream input( archive );
  ostringstream output;
  const auto header = header_for( 0u, payload.size(), payload.size() );

  BOOST_CHECK_EQUAL(
    ZipArchive::read_file_from_zip(input, header, output, {1024u, 100u}),
    payload.size() );
  BOOST_CHECK_EQUAL( output.str(), payload );
}


BOOST_AUTO_TEST_CASE( declared_member_and_ratio_limits_are_enforced )
{
  ZipArchive::FilenameToZipHeaderMap headers;
  headers["large"] = header_for( 8u, 100u, 1025u );
  BOOST_CHECK_THROW(
    ZipArchive::validate_archive_for_extraction(headers, {1024u, 100u}, 2048u),
    std::exception );

  headers.clear();
  headers["ratio"] = header_for( 8u, 1u, 101u );
  BOOST_CHECK_THROW(
    ZipArchive::validate_archive_for_extraction(headers, {1024u, 100u}, 2048u),
    std::exception );
}


BOOST_AUTO_TEST_CASE( aggregate_limit_is_enforced )
{
  ZipArchive::FilenameToZipHeaderMap headers;
  headers["one"] = header_for( 0u, 60u, 60u );
  headers["two"] = header_for( 0u, 60u, 60u );
  BOOST_CHECK_THROW(
    ZipArchive::validate_archive_for_extraction(headers, {100u, 100u}, 100u),
    std::exception );
}


BOOST_AUTO_TEST_CASE( actual_inflation_cannot_exceed_limit_or_declared_size )
{
  const string expanded( 4096u, 'A' );
  const string compressed = raw_deflate( expanded );
  const string archive = local_member( compressed, 8u, compressed.size(), 10u );
  istringstream input( archive );
  ostringstream output;
  const auto forged_header = header_for( 8u, compressed.size(), 10u );

  BOOST_CHECK_THROW(
    ZipArchive::read_file_from_zip(input, forged_header, output, {100u, 100u}),
    std::exception );
  BOOST_CHECK_LE( output.str().size(), 100u );
}


BOOST_AUTO_TEST_CASE( truncated_deflate_is_rejected )
{
  const string expanded( 1024u, 'B' );
  string compressed = raw_deflate( expanded );
  const uint32_t declared_compressed = static_cast<uint32_t>(compressed.size());
  compressed.pop_back();
  const string archive = local_member( compressed, 8u, declared_compressed, expanded.size() );
  istringstream input( archive );
  ostringstream output;
  const auto header = header_for( 8u, declared_compressed, expanded.size() );

  BOOST_CHECK_THROW(
    ZipArchive::read_file_from_zip(input, header, output, {2048u, 100u}),
    std::exception );
}
