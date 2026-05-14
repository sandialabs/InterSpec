// libFuzzer driver for DetectorPeakResponse::parseSingleCsvLineRelEffDrf.
// See target/testing/fuzzing/README.md for build & run instructions.

#include <cstdint>
#include <string>

#include "InterSpec/DetectorPeakResponse.h"

extern "C" int LLVMFuzzerTestOneInput( const uint8_t *data, size_t size )
{
  if( size > 64 * 1024 )
    return 0;

  std::string line( reinterpret_cast<const char *>(data), size );
  try
  {
    DetectorPeakResponse::parseSingleCsvLineRelEffDrf( line );
  }catch( const std::exception & )
  {
    // Expected on malformed inputs; only crashes / sanitizer reports fail.
  }
  return 0;
}
