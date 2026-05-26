// libFuzzer driver for RelActCalcManual::PeakCsvInput::peak_csv_to_peaks.
// See target/testing/fuzzing/README.md for build & run instructions.

#include <cstdint>
#include <string>
#include <sstream>

#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalcManual.h"

extern "C" int LLVMFuzzerTestOneInput( const uint8_t *data, size_t size )
{
  if( size > 256 * 1024 )
    return 0;

  std::string buf( reinterpret_cast<const char *>(data), size );
  std::istringstream is( buf );
  try
  {
    RelActCalcManual::PeakCsvInput::peak_csv_to_peaks( is );
  }catch( const std::exception & )
  {
  }
  return 0;
}
