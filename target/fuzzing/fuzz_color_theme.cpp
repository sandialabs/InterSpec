// libFuzzer driver for ColorTheme::fromJson.
// See target/testing/fuzzing/README.md for build & run instructions.

#include <cstdint>
#include <string>

#include "InterSpec/ColorTheme.h"

extern "C" int LLVMFuzzerTestOneInput( const uint8_t *data, size_t size )
{
  if( size > 64 * 1024 )
    return 0;

  const std::string json( reinterpret_cast<const char *>(data), size );
  ColorTheme theme;
  try
  {
    ColorTheme::fromJson( json, theme );
  }catch( const std::exception & )
  {
    // Parse errors are expected; we're hunting for crashes / UB only.
  }
  return 0;
}
