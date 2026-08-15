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
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#include <set>
#include <memory>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE BatchSampleSelect_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/BatchSampleSelect.h"

using namespace std;
using namespace boost::unit_test;

namespace
{
/** Adds one Measurement, with usable gamma data, to `file`. */
void add_meas( SpecUtils::SpecFile &file,
               const int sample_number,
               const string &det_name,
               const SpecUtils::SourceType source_type )
{
  auto m = make_shared<SpecUtils::Measurement>();

  auto counts = make_shared<vector<float>>( 64, 1.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( counts->size(), {0.0f, 3.0f}, {} );

  m->set_gamma_counts( counts, 10.0f, 10.0f );
  m->set_energy_calibration( cal );
  m->set_sample_number( sample_number );
  m->set_detector_name( det_name );
  m->set_source_type( source_type );

  file.add_measurement( m, false );
}//void add_meas(...)


/** Convenience for building a file from a list of {sample number, source type}. */
shared_ptr<SpecUtils::SpecFile> make_file( const vector<pair<int,SpecUtils::SourceType>> &samples )
{
  auto file = make_shared<SpecUtils::SpecFile>();
  for( const pair<int,SpecUtils::SourceType> &s : samples )
    add_meas( *file, s.first, "Det1", s.second );
  file->cleanup_after_load();
  return file;
}//shared_ptr<SpecUtils::SpecFile> make_file(...)
}//namespace


BOOST_AUTO_TEST_CASE( SingleForegroundSample )
{
  const shared_ptr<SpecUtils::SpecFile> file
        = make_file( { {1, SpecUtils::SourceType::Foreground} } );

  BOOST_CHECK( BatchSampleSelect::classify_sample( *file, 1 )
               == BatchSampleSelect::SampleClass::Foreground );

  const set<int> answer = BatchSampleSelect::candidate_foreground_samples( *file );
  BOOST_CHECK_EQUAL( answer.size(), 1 );
  BOOST_CHECK( answer == set<int>{1} );

  BOOST_CHECK_EQUAL( BatchSampleSelect::single_foreground_sample( *file ), 1 );
}//BOOST_AUTO_TEST_CASE( SingleForegroundSample )


BOOST_AUTO_TEST_CASE( LoneBackgroundFallback )
{
  const shared_ptr<SpecUtils::SpecFile> file
        = make_file( { {1, SpecUtils::SourceType::Background} } );

  BatchSampleSelect::Options no_fallback;
  no_fallback.allow_lone_background_fallback = false;
  BOOST_CHECK( BatchSampleSelect::candidate_foreground_samples( *file, no_fallback ).empty() );

  BatchSampleSelect::Options with_fallback;
  with_fallback.allow_lone_background_fallback = true;
  BOOST_CHECK( BatchSampleSelect::candidate_foreground_samples( *file, with_fallback ) == set<int>{1} );

  // This is the long-standing behaviour of being able to peak-fit a background-marked file
  BOOST_CHECK_EQUAL( BatchSampleSelect::single_foreground_sample( *file ), 1 );
}//BOOST_AUTO_TEST_CASE( LoneBackgroundFallback )


BOOST_AUTO_TEST_CASE( ThreeForegroundsPlusBackgroundAndIntrinsic )
{
  // This is the shape that makes the "batch tool" link appear, and that `each` mode splits
  const shared_ptr<SpecUtils::SpecFile> file = make_file( {
    {1, SpecUtils::SourceType::Foreground},
    {2, SpecUtils::SourceType::Foreground},
    {3, SpecUtils::SourceType::Foreground},
    {4, SpecUtils::SourceType::Background},
    {5, SpecUtils::SourceType::IntrinsicActivity},
  } );

  const set<int> answer = BatchSampleSelect::candidate_foreground_samples( *file );
  BOOST_CHECK( answer == (set<int>{1,2,3}) );

  BOOST_CHECK_EQUAL( BatchSampleSelect::single_background_sample( *file ), 4 );

  // Ambiguous, so must throw
  BOOST_CHECK_THROW( BatchSampleSelect::single_foreground_sample( *file ), std::runtime_error );
}//BOOST_AUTO_TEST_CASE( ThreeForegroundsPlusBackgroundAndIntrinsic )


BOOST_AUTO_TEST_CASE( MixUnknownWithForeground )
{
  const shared_ptr<SpecUtils::SpecFile> file = make_file( {
    {1, SpecUtils::SourceType::Foreground},
    {2, SpecUtils::SourceType::Foreground},
    {3, SpecUtils::SourceType::Unknown},
    {4, SpecUtils::SourceType::Unknown},
  } );

  BatchSampleSelect::Options fore_only;
  fore_only.mix_unknown_with_foreground = false;
  BOOST_CHECK( BatchSampleSelect::candidate_foreground_samples( *file, fore_only ) == (set<int>{1,2}) );

  BatchSampleSelect::Options mixed;
  mixed.mix_unknown_with_foreground = true;
  BOOST_CHECK( BatchSampleSelect::candidate_foreground_samples( *file, mixed ) == (set<int>{1,2,3,4}) );
}//BOOST_AUTO_TEST_CASE( MixUnknownWithForeground )


BOOST_AUTO_TEST_CASE( OnlyUnknownSamples )
{
  const shared_ptr<SpecUtils::SpecFile> file = make_file( {
    {1, SpecUtils::SourceType::Unknown},
    {2, SpecUtils::SourceType::Unknown},
    {3, SpecUtils::SourceType::Unknown},
  } );

  BatchSampleSelect::Options fore_only;
  fore_only.mix_unknown_with_foreground = false;
  BOOST_CHECK( BatchSampleSelect::candidate_foreground_samples( *file, fore_only ) == (set<int>{1,2,3}) );

  BatchSampleSelect::Options mixed;
  mixed.mix_unknown_with_foreground = true;
  BOOST_CHECK( BatchSampleSelect::candidate_foreground_samples( *file, mixed ) == (set<int>{1,2,3}) );
}//BOOST_AUTO_TEST_CASE( OnlyUnknownSamples )


BOOST_AUTO_TEST_CASE( DetectorsDisagreeGiveUnknown )
{
  // The regression test for the classification change: when the Measurements making up one sample
  //  carry different source types, the answer must be `Unknown`, and must not depend on the order
  //  the detectors happen to be in.  This matches what `SpecUtils::SpecFile::sum_measurements()`
  //  assigns the summed Measurement.
  for( int order = 0; order < 2; ++order )
  {
    SpecUtils::SpecFile file;
    if( order == 0 )
    {
      add_meas( file, 1, "Det1", SpecUtils::SourceType::Background );
      add_meas( file, 1, "Det2", SpecUtils::SourceType::Foreground );
    }else
    {
      add_meas( file, 1, "Det1", SpecUtils::SourceType::Foreground );
      add_meas( file, 1, "Det2", SpecUtils::SourceType::Background );
    }
    file.cleanup_after_load();

    BOOST_CHECK( BatchSampleSelect::classify_sample( file, 1 )
                 == BatchSampleSelect::SampleClass::Unknown );

    // And it must agree with what SpecUtils itself would say
    const shared_ptr<const SpecUtils::Measurement> summed
          = file.sum_measurements( {1}, file.detector_names(), nullptr );
    BOOST_REQUIRE( !!summed );
    BOOST_CHECK( summed->source_type() == SpecUtils::SourceType::Unknown );
  }//for( int order = 0; order < 2; ++order )
}//BOOST_AUTO_TEST_CASE( DetectorsDisagreeGiveUnknown )


BOOST_AUTO_TEST_CASE( EmptyFile )
{
  SpecUtils::SpecFile file;
  file.cleanup_after_load();

  BOOST_CHECK( BatchSampleSelect::candidate_foreground_samples( file ).empty() );
  BOOST_CHECK_THROW( BatchSampleSelect::single_foreground_sample( file ), std::runtime_error );
}//BOOST_AUTO_TEST_CASE( EmptyFile )


BOOST_AUTO_TEST_CASE( MultiSampleHandlingStrings )
{
  BOOST_CHECK( BatchSampleSelect::multi_sample_handling_from_str( "auto" )
               == BatchSampleSelect::MultiSampleHandling::Auto );
  BOOST_CHECK( BatchSampleSelect::multi_sample_handling_from_str( "  EACH " )
               == BatchSampleSelect::MultiSampleHandling::EachSampleSeparately );
  BOOST_CHECK( BatchSampleSelect::multi_sample_handling_from_str( "Sum" )
               == BatchSampleSelect::MultiSampleHandling::SumAllSamples );
  BOOST_CHECK( BatchSampleSelect::multi_sample_handling_from_str( "" )
               == BatchSampleSelect::MultiSampleHandling::Auto );
  BOOST_CHECK_THROW( BatchSampleSelect::multi_sample_handling_from_str( "bogus" ), std::runtime_error );

  BOOST_CHECK_EQUAL( string( BatchSampleSelect::to_str( BatchSampleSelect::MultiSampleHandling::Auto ) ), "Auto" );
}//BOOST_AUTO_TEST_CASE( MultiSampleHandlingStrings )


BOOST_AUTO_TEST_CASE( ExpandInputFilesAutoIsPassThrough )
{
  // `Auto` must not parse anything, and must give exactly one work item per input file, with the
  //  output name and sample numbers left alone.
  const vector<string> files{ "/does/not/exist/a.n42", "/does/not/exist/b.n42" };

  const vector<BatchSampleSelect::InputWorkItem> items
        = BatchSampleSelect::expand_input_files( files, {}, BatchSampleSelect::MultiSampleHandling::Auto );

  BOOST_REQUIRE_EQUAL( items.size(), 2 );
  BOOST_CHECK_EQUAL( items[0].filename, files[0] );
  BOOST_CHECK_EQUAL( items[0].output_base_name, "a.n42" );
  BOOST_CHECK_EQUAL( items[0].label, "a.n42" );
  BOOST_CHECK( items[0].foreground_sample_numbers.empty() );
  BOOST_CHECK( !items[0].needs_private_copy );
  BOOST_CHECK( !items[0].source );
  BOOST_CHECK_EQUAL( items[1].output_base_name, "b.n42" );
}//BOOST_AUTO_TEST_CASE( ExpandInputFilesAutoIsPassThrough )
