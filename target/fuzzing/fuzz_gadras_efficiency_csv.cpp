// libFuzzer driver for DetectorPeakResponse::fromGadrasDefinition - mutates
// the Efficiency.csv side while the Detector.dat side is held to a known-good
// baseline.  Seed corpus: data/GenericGadrasDetectors/<detector>/Efficiency.csv.
//
// See target/testing/fuzzing/README.md for build & run instructions.

#include <cstdint>
#include <sstream>
#include <string>

#include "InterSpec/DetectorPeakResponse.h"

#include "gadras_baseline.h"

extern "C" int LLVMFuzzerTestOneInput( const uint8_t *data, size_t size )
{
  if( size > 256 * 1024 )
    return 0;

  std::string csv_buf( reinterpret_cast<const char *>(data), size );
  std::istringstream csv_stream( csv_buf );
  std::istringstream dat_stream( fuzz_gadras_baseline::k_detector_dat );

  DetectorPeakResponse drf;
  try
  {
    drf.fromGadrasDefinition( csv_stream, dat_stream );
  }catch( const std::exception & )
  {
  }
  return 0;
}
