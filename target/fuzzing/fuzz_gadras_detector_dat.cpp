// libFuzzer driver for DetectorPeakResponse::fromGadrasDefinition - mutates
// the Detector.dat side while the Efficiency.csv side is held to a known-good
// baseline.  Exercises both the old-style line-numbered parser and the XML
// parser branch (selected by the parser if the first line contains "xml").
// Seed corpus: data/GenericGadrasDetectors/<detector>/Detector.dat.
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

  std::string dat_buf( reinterpret_cast<const char *>(data), size );
  std::istringstream csv_stream( fuzz_gadras_baseline::k_efficiency_csv );
  std::istringstream dat_stream( dat_buf );

  DetectorPeakResponse drf;
  try
  {
    drf.fromGadrasDefinition( csv_stream, dat_stream );
  }catch( const std::exception & )
  {
  }
  return 0;
}
