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

#include "InterSpec_config.h"

#include <set>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchSampleSelect.h"

using namespace std;


namespace
{
  /** Inserts `suffix` immediately before the file extension of a leaf filename.
   ("foo.n42","_sample3") -> "foo_sample3.n42";  ("foo","_summed") -> "foo_summed".
   */
  string insert_filename_suffix( const string &leaf_name, const string &suffix )
  {
    const string ext = SpecUtils::file_extension( leaf_name );
    if( ext.empty() )
      return leaf_name + suffix;

    return leaf_name.substr( 0, leaf_name.size() - ext.size() ) + suffix + ext;
  }//insert_filename_suffix(...)
}//namespace


namespace BatchSampleSelect
{

SampleClass classify_sample( const SpecUtils::SpecFile &file, const int sample_number )
{
  const vector<shared_ptr<const SpecUtils::Measurement>> meass
                                                    = file.sample_measurements( sample_number );
  if( meass.empty() )
    return SampleClass::Empty;

  // Match `InterSpec::validForegroundSamples()`: derived data is only excluded when the file also
  //  holds non-derived data, so a file that is entirely derived data is still analyzable.
  const bool check_derived = (file.contains_derived_data() && file.contains_non_derived_data());

  // Replicate the rule `SpecUtils::SpecFile::sum_measurements()` uses to assign the summed
  //  Measurement its source type: collect the types of every Measurement, and only commit to one
  //  if they all agree.  Doing it this way, rather than actually calling `sum_measurements()`,
  //  avoids summing channel data - this function is called from GUI code that runs often.
  set<SpecUtils::SourceType> source_types;
  bool any_usable_gamma = false;

  for( const shared_ptr<const SpecUtils::Measurement> &m : meass )
  {
    if( !m )
      continue;

    if( check_derived && m->derived_data_properties() )
      return SampleClass::Derived;

    source_types.insert( m->source_type() );

    const shared_ptr<const vector<float>> counts = m->gamma_counts();
    if( counts && (counts->size() > 3) && m->energy_calibration() && m->energy_calibration()->valid() )
      any_usable_gamma = true;
  }//for( loop over Measurements of this sample )

  if( !any_usable_gamma )
    return SampleClass::Empty;

  if( source_types.size() != 1 )
    return SampleClass::Unknown;

  switch( *begin(source_types) )
  {
    case SpecUtils::SourceType::IntrinsicActivity:
    case SpecUtils::SourceType::Calibration:
      return SampleClass::Other;

    case SpecUtils::SourceType::Background:
      return SampleClass::Background;

    case SpecUtils::SourceType::Foreground:
      return SampleClass::Foreground;

    case SpecUtils::SourceType::Unknown:
      return SampleClass::Unknown;
  }//switch( *begin(source_types) )

  return SampleClass::Unknown;
}//SampleClass classify_sample( const SpecUtils::SpecFile &, const int )


SampleBuckets classify_samples( const SpecUtils::SpecFile &file )
{
  SampleBuckets answer;

  for( const int sample : file.sample_numbers() )
  {
    switch( classify_sample( file, sample ) )
    {
      case SampleClass::Foreground: answer.foreground.insert( sample ); break;
      case SampleClass::Unknown:    answer.unknown.insert( sample );    break;
      case SampleClass::Background: answer.background.insert( sample ); break;
      case SampleClass::Other:      answer.other.insert( sample );      break;
      case SampleClass::Derived:    break;
      case SampleClass::Empty:      break;
    }//switch( classify_sample( file, sample ) )
  }//for( const int sample : file.sample_numbers() )

  return answer;
}//SampleBuckets classify_samples( const SpecUtils::SpecFile & )


std::set<int> candidate_foreground_samples( const SpecUtils::SpecFile &file, const Options &options )
{
  const SampleBuckets buckets = classify_samples( file );

  if( !buckets.foreground.empty() )
  {
    if( !options.mix_unknown_with_foreground )
      return buckets.foreground;

    set<int> answer = buckets.foreground;
    answer.insert( begin(buckets.unknown), end(buckets.unknown) );
    return answer;
  }//if( !buckets.foreground.empty() )

  if( !buckets.unknown.empty() )
    return buckets.unknown;

  // A file whose sole record is marked background is still something we want to be able to
  //  analyze - e.g., fitting peaks in a background spectrum.
  if( options.allow_lone_background_fallback && (buckets.background.size() == 1) )
    return buckets.background;

  return set<int>{};
}//std::set<int> candidate_foreground_samples( const SpecUtils::SpecFile &, const Options & )


int single_foreground_sample( const SpecUtils::SpecFile &file )
{
  Options opts;
  opts.allow_lone_background_fallback = true;
  opts.mix_unknown_with_foreground = true;

  const set<int> samples = candidate_foreground_samples( file, opts );
  if( samples.size() != 1 )
    throw runtime_error( "Sample number to use could not be uniquely identified." );

  return *begin(samples);
}//int single_foreground_sample( const SpecUtils::SpecFile & )


int single_background_sample( const SpecUtils::SpecFile &file )
{
  const SampleBuckets buckets = classify_samples( file );

  if( buckets.background.size() != 1 )
    throw runtime_error( "Sample number to use could not be uniquely identified." );

  return *begin(buckets.background);
}//int single_background_sample( const SpecUtils::SpecFile & )


const char *to_str( const MultiSampleHandling handling )
{
  switch( handling )
  {
    case MultiSampleHandling::Auto:                 return "Auto";
    case MultiSampleHandling::EachSampleSeparately: return "EachSampleSeparately";
    case MultiSampleHandling::SumAllSamples:        return "SumAllSamples";
  }//switch( handling )

  return "Auto";
}//const char *to_str( const MultiSampleHandling )


MultiSampleHandling multi_sample_handling_from_str( std::string val )
{
  const std::string orig = val;
  SpecUtils::trim( val );
  SpecUtils::to_lower_ascii( val );

  if( val.empty() || (val == "auto") )
    return MultiSampleHandling::Auto;

  // The long forms are what `to_str()` emits into the report JSON, so accept them here too - a
  //  user copying a value out of a results JSON should be able to paste it on the command line.
  if( (val == "each") || (val == "each-sample") || (val == "separate")
     || (val == "eachsampleseparately") )
    return MultiSampleHandling::EachSampleSeparately;

  if( (val == "sum") || (val == "sum-all") || (val == "sumallsamples") )
    return MultiSampleHandling::SumAllSamples;

  throw runtime_error( "Invalid multi-sample-handling value ('" + orig
                       + "'); must be one of \"auto\", \"each\", or \"sum\"." );
}//MultiSampleHandling multi_sample_handling_from_str( std::string )


std::vector<InputWorkItem>
  expand_input_file( const std::string &filename,
                     const size_t input_index,
                     const std::shared_ptr<SpecMeas> &cached_file,
                     const MultiSampleHandling handling )
{
  vector<InputWorkItem> items;

  InputWorkItem base;
  base.input_index = input_index;
  base.filename = filename;
  base.output_base_name = SpecUtils::filename( filename );
  base.label = base.output_base_name;
  base.source = cached_file;

  // Fast path: byte-for-byte the historical behaviour, and nothing is parsed here, so `Auto`
  //  costs exactly what it always did.
  if( handling == MultiSampleHandling::Auto )
  {
    items.push_back( base );
    return items;
  }

  shared_ptr<SpecMeas> parsed = cached_file;
  if( !parsed )
  {
    parsed = make_shared<SpecMeas>();
    if( !parsed->load_file( filename, SpecUtils::ParserType::Auto, filename ) )
    {
      // Don't pre-empt the downstream error message - emit one pass-through item with a null
      //  source, so the fit function produces its canonical "Couldnt read in ..." warning.
      items.push_back( base );
      return items;
    }
  }//if( !parsed )

  // Passthrough / search-mode files are never split or summed.  Note we deliberately dont keep
  //  our parse of them - they are the largest files, and the fit function will load what it
  //  needs; holding it here would just inflate peak memory.
  if( parsed->passthrough() )
  {
    items.push_back( base );
    return items;
  }

  Options sel;
  sel.allow_lone_background_fallback = true;
  sel.mix_unknown_with_foreground = true;
  const set<int> candidates = candidate_foreground_samples( *parsed, sel );

  if( candidates.size() <= 1 )
  {
    // Nothing to split, so - as above - dont hold onto our parse.
    base.foreground_sample_numbers = candidates;  //possibly empty -> auto-detect downstream
    items.push_back( base );
    return items;
  }

  // From here on the file really is being split or summed, so hand our parse downstream rather
  //  than making the fit function read the file again.
  base.source = parsed;

  if( handling == MultiSampleHandling::SumAllSamples )
  {
    base.foreground_sample_numbers = candidates;
    base.output_base_name = insert_filename_suffix( base.output_base_name, "_summed" );
    base.label += " (summed samples)";
    items.push_back( base );
    return items;
  }

  assert( handling == MultiSampleHandling::EachSampleSeparately );

  size_t made = 0;
  for( const int sample : candidates )
  {
    InputWorkItem item = base;
    item.foreground_sample_numbers = set<int>{ sample };
    item.output_base_name = insert_filename_suffix( base.output_base_name,
                                                    "_sample" + std::to_string(sample) );
    item.label = base.label + " (sample " + std::to_string(sample) + ")";

    // All but the last item need their own copy, since the fit functions modify the file handed
    //  to them.  The last item can use `parsed` itself, as nothing else will refer to it.
    made += 1;
    item.needs_private_copy = (made < candidates.size());

    items.push_back( item );
  }//for( const int sample : candidates )

  return items;
}//std::vector<InputWorkItem> expand_input_file(...)


std::vector<InputWorkItem>
  expand_input_files( const std::vector<std::string> &files,
                      const std::vector<std::shared_ptr<SpecMeas>> &cached_files,
                      const MultiSampleHandling handling )
{
  vector<InputWorkItem> items;

  for( size_t i = 0; i < files.size(); ++i )
  {
    const shared_ptr<SpecMeas> cached = (i < cached_files.size()) ? cached_files[i] : nullptr;
    const vector<InputWorkItem> file_items = expand_input_file( files[i], i, cached, handling );
    items.insert( end(items), begin(file_items), end(file_items) );
  }

  return items;
}//std::vector<InputWorkItem> expand_input_files(...)

}//namespace BatchSampleSelect
