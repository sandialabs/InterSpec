// libFuzzer driver for SpecMeas::load_from_N42 (the InterSpec extension of
// SpecUtils' N42 loader, which also restores peaks/DRF/state).
//
// See target/testing/fuzzing/README.md for build & run instructions.

#include <cstdint>
#include <string>
#include <sstream>

#include "InterSpec/SpecMeas.h"

extern "C" int LLVMFuzzerTestOneInput( const uint8_t *data, size_t size )
{
  if( size > 1 * 1024 * 1024 )
    return 0;

  std::string buf( reinterpret_cast<const char *>(data), size );
  std::istringstream is( buf );
  SpecMeas meas;
  try
  {
    meas.load_from_N42( is );
  }catch( const std::exception & )
  {
  }
  return 0;
}
