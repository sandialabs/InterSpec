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

#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"

#include "Session.h"
#if( LIGHT_PEAK_FILE_IO )
#include "PeakFileIo.h"
#endif

using namespace std;
using json = nlohmann::json;

namespace
{
  Session::SpecType type_from_json( const json &p )
  {
    return Session::typeFromName( p.value( "type", string("FOREGROUND") ) );
  }

  SpecUtils::SpectrumType to_specutils( const Session::SpecType type )
  {
    switch( type )
    {
      case Session::Background: return SpecUtils::SpectrumType::Background;
      case Session::Secondary:  return SpecUtils::SpectrumType::SecondForeground;
      default: break;
    }
    return SpecUtils::SpectrumType::Foreground;
  }

  /** Source type of a sample: that of its last measurement, like InterSpec's file manager. */
  SpecUtils::SourceType sample_source_type( const SpecUtils::SpecFile &spec, const int sample )
  {
    const auto meass = spec.sample_measurements( sample );
    return meass.empty() ? SpecUtils::SourceType::Unknown : meass.back()->source_type();
  }

  /** Groups sample numbers into contiguous runs, in the order of the file's sample numbers. */
  vector<pair<int,int>> sample_runs( const set<int> &all_samples, const set<int> &wanted )
  {
    vector<pair<int,int>> runs;
    bool in_run = false;
    int start = 0, prev = 0;
    for( const int s : all_samples )
    {
      const bool want = wanted.count( s );
      if( want && !in_run )
        start = s;
      else if( !want && in_run )
        runs.push_back( {start, prev} );
      in_run = want;
      prev = s;
    }
    if( in_run )
      runs.push_back( {start, prev} );
    return runs;
  }
}//namespace


const char *Session::typeName( const SpecType type )
{
  switch( type )
  {
    case Foreground: return "FOREGROUND";
    case Background: return "BACKGROUND";
    case Secondary:  return "SECONDARY";
    case NumSpecType: break;
  }
  return "";
}//typeName(...)


Session::SpecType Session::typeFromName( const string &name )
{
  for( int i = 0; i < NumSpecType; ++i )
  {
    if( SpecUtils::iequals_ascii( name, typeName( SpecType(i) ) ) )
      return SpecType(i);
  }
  throw runtime_error( "Invalid spectrum type '" + name + "'" );
}//typeFromName(...)


Session::Session()
  : m_next_file_id( 0 )
{
}


Session::~Session()
{
}


json Session::call( const string &method, const json &params )
{
  if( method == "getState" )
  {
    json answer = state( PartAll, true );
    answer["features"] = { {"peakFileIo", static_cast<bool>(LIGHT_PEAK_FILE_IO)} };
    return answer;
  }
  if( method == "loadFile" )          return loadFile( params );
  if( method == "unload" )            return unload( params );
  if( method == "setSamples" )        return setSamples( params );
  if( method == "stepSample" )        return stepSample( params );
  if( method == "timeDrag" )          return timeDrag( params );
  if( method == "setDetectors" )      return setDetectors( params );
  if( method == "exportFile" )        return exportFile( params );
  if( method == "setRefLibrary" )     return setRefLibrary( params );
  if( method == "setShownRefLines" )  return setShownRefLines( params );
  if( method == "fitPeakAt" )         return fitPeakAt( params );
  if( method == "deletePeakAt" )      return deletePeakAt( params );
  if( method == "erasePeaks" )        return erasePeaks( params );
  if( method == "roiDrag" )           return roiDrag( params );
  if( method == "fitRoiDrag" )        return fitRoiDrag( params );
  if( method == "peakInfoAt" )        return peakInfoAt( params );
  if( method == "setContinuumType" )  return setContinuumType( params );
  if( method == "setSkewType" )       return setSkewType( params );
  if( method == "refitRoi" )          return refitRoi( params );
  if( method == "setPeakProperty" )   return setPeakProperty( params );
  if( method == "clearPeaks" )        return clearPeaks( params );
  if( method == "exportPeakCsv" )     return exportPeakCsv( params );
  if( method == "searchPeaks" )       return searchPeaks( params );
  if( method == "setEnergyCal" )      return setEnergyCal( params );
  if( method == "fitEnergyCal" )      return fitEnergyCal( params );
  if( method == "revertEnergyCal" )   return revertEnergyCal( params );
  if( method == "convertCalCoefs" )   return convertCalCoefs( params );
  if( method == "exportCALp" )        return exportCALp( params );
  if( method == "importCALp" )        return importCALp( params );

  throw runtime_error( "Unknown method '" + method + "'" );
}//json call(...)


json Session::state( const unsigned parts, const bool resetDomain ) const
{
  json answer = json::object();

  if( parts & PartFiles )
    answer["files"] = filesJson();

  if( parts & PartDetectors )
    answer["detectors"] = detectorsJson();

  if( parts & PartSpectra )
  {
    answer["spectra"] = spectraJson();
    answer["resetDomain"] = resetDomain;
    answer["isHighRes"] = isHighRes();
  }

  if( parts & PartTimeChart )
  {
    answer["timeChart"] = timeChartJson();
    answer["timeHighlights"] = timeHighlightsJson();
  }

  if( parts & PartPeaks )
  {
    answer["peaks"] = peaksJson();
    answer["peakList"] = peakListJson();
  }

  if( parts & PartEnergyCal )
    answer["energyCal"] = energyCalJson();

  return answer;
}//json state(...)


json Session::loadFile( const json &p )
{
  const string path = p.at( "path" ).get<string>();
  const string name = p.value( "name", path );
  const SpecType type = type_from_json( p );

#if( LIGHT_PEAK_FILE_IO )
  // Files InterSpec wrote keep their sample numbers, which their saved peaks are keyed by
  const PeakFileIo::SavedPeaksLocation saved_loc = PeakFileIo::locate_saved_peaks( path );
  auto spec = make_shared<PeakFileIo::SampleKeepingSpecFile>( saved_loc.type == PeakFileIo::SavedPeaksLocation::Type::N42 );
#else
  auto spec = make_shared<SpecUtils::SpecFile>();
#endif
  if( !spec->load_file( path, SpecUtils::ParserType::Auto, name ) )
    throw runtime_error( "Could not parse '" + name + "' as a spectrum file." );

  if( spec->num_measurements() == 0 )
    throw runtime_error( "'" + name + "' contains no measurements." );

  vector<string> messages;

  // InterSpec asks the user about these; the light version just picks the common choice.
  if( spec->contains_derived_data() && spec->contains_non_derived_data() )
  {
    spec->keep_derived_data_variant( SpecUtils::SpecFile::DerivedVariantToKeep::NonDerived );
    messages.push_back( "File has derived and raw data; showing the raw data." );
  }

  const set<string> cal_variants = spec->energy_cal_variants();
  if( cal_variants.size() > 1 )
  {
    spec->keep_energy_cal_variants( { *begin(cal_variants) } );
    messages.push_back( "File has multiple energy binnings; using '" + *begin(cal_variants) + "'." );
  }

  auto file = make_shared<LoadedFile>();
  file->id = m_next_file_id++;
  file->name = name;
  file->spec = spec;
  for( const auto &m : spec->measurements() )
    file->orig_cals[m.get()] = m->energy_calibration();

  // For a file saved with peaks, the samples that were displayed (InterSpec shows them again)
  set<int> saved_display;
#if( LIGHT_PEAK_FILE_IO )
  try
  {
    PeakFileIo::FilePeaks saved = PeakFileIo::read_file_peaks( path, saved_loc, *spec, m_ref_lib );
    set<const PeakDef *> unique_peaks;
    for( const auto &sp : saved.peaks )
    {
      for( const auto &peak : sp.second )
        unique_peaks.insert( peak.get() );
    }
    file->peaks = std::move( saved.peaks );

    if( SpecUtils::iequals_ascii( saved.display_type, SpecUtils::descriptionText( to_specutils(type) ) ) )
      saved_display = saved.displayed_samples;

    const size_t npeaks = unique_peaks.size();
    if( npeaks )
      messages.push_back( "Loaded " + std::to_string(npeaks) + " peak" + ((npeaks == 1) ? "" : "s")
                          + " saved in the file." );
    messages.insert( end(messages), begin(saved.warnings), end(saved.warnings) );
  }catch( std::exception &e )
  {
    messages.push_back( string("Could not read the peaks saved in the file: ") + e.what() );
  }
#endif

  shared_ptr<const PeakFitDetPrefs> prefs = PeakFitDetPrefs::defaultForDetectorType(
                    static_cast<int>(spec->detector_type()), spec->manufacturer(), spec->instrument_model() );
  auto own_prefs = prefs ? make_shared<PeakFitDetPrefs>( *prefs ) : make_shared<PeakFitDetPrefs>();
  file->prefs = own_prefs;

  const set<int> &all_samples = spec->sample_numbers();
  const bool passthrough = spec->passthrough();
  set<int> selected, background;
  bool note_sample_choice = false;

  switch( type )
  {
    case Foreground:
    case Secondary:
    {
      if( !passthrough && (all_samples.size() > 1) )
      {
        // Port of SpecMeasManager::displayFile: prefer the last Foreground sample, then the first
        //  Unknown, then the first sample; note an unambiguous in-file background.
        int fore_sample = numeric_limits<int>::min(), back_sample = numeric_limits<int>::min();
        size_t num_fore = 0, num_back = 0, num_intrinsic = 0;
        bool have_choice = false;
        for( const int sample : all_samples )
        {
          switch( sample_source_type( *spec, sample ) )
          {
            case SpecUtils::SourceType::IntrinsicActivity:
              ++num_intrinsic;
              break;

            case SpecUtils::SourceType::Foreground:
              ++num_fore;
              fore_sample = sample;
              have_choice = true;
              break;

            case SpecUtils::SourceType::Background:
              ++num_back;
              back_sample = sample;
              break;

            case SpecUtils::SourceType::Calibration:
              break;

            case SpecUtils::SourceType::Unknown:
              if( !have_choice )
              {
                fore_sample = sample;
                have_choice = true;
              }
              break;
          }//switch( source type )

          if( num_fore && num_back )
            break;
        }//for( loop over samples )

        if( !have_choice )
          fore_sample = *begin(all_samples);
        selected.insert( fore_sample );

        if( num_intrinsic )
          messages.push_back( "File has intrinsic activity spectra; these are not shown by default." );
        note_sample_choice = true;

        const Slot &old_back = m_slots[Background];
        const bool new_back_better = !old_back.file
            || (old_back.file->spec->num_gamma_channels() != spec->num_gamma_channels())
            || (old_back.file->spec->instrument_id() != spec->instrument_id());

        if( (type == Foreground) && (num_back == 1) && (back_sample != fore_sample) && new_back_better )
          background.insert( back_sample );
      }else if( passthrough )
      {
        // All samples, except background and calibration; the file's background samples become
        //  the background.
        size_t ncalibration = 0;
        for( const int sample : all_samples )
        {
          const SpecUtils::SourceType st = sample_source_type( *spec, sample );
          const bool back = (st == SpecUtils::SourceType::Background);
          const bool calib = (st == SpecUtils::SourceType::Calibration);
          ncalibration += calib;
          if( back && (type != Secondary) )
            background.insert( sample );
          if( !back && !calib && (st != SpecUtils::SourceType::IntrinsicActivity) )
            selected.insert( sample );
        }//for( loop over samples )

        if( ncalibration )
          messages.push_back( "File contains " + std::to_string(ncalibration)
                              + " calibration spectra; these are not shown by default." );
      }else
      {
        selected = all_samples;
      }//if( !passthrough && multiple samples ) / else passthrough / else

      if( selected.empty() )
      {
        // e.g., passthrough data that is all marked background; take the first sample with counts
        for( const int sample : all_samples )
        {
          for( const auto &m : spec->sample_measurements( sample ) )
          {
            if( m->gamma_count_sum() > 0.0 )
              selected.insert( sample );
          }
          if( !selected.empty() )
            break;
        }
        background.clear();
      }//if( selected.empty() )
      break;
    }//case Foreground / Secondary

    case Background:
    {
      if( all_samples.size() > 1 )
      {
        // A different file than the foreground: prefer its Foreground record (a background taken as
        //  a normal measurement), then Unknown, then Background; else the first sample.
        const bool different_file = true; //a newly loaded file is always a different file
        int fore = numeric_limits<int>::min(), back = fore, unknown = fore;
        for( const int sample : all_samples )
        {
          switch( sample_source_type( *spec, sample ) )
          {
            case SpecUtils::SourceType::Foreground: if( fore == numeric_limits<int>::min() ) fore = sample; break;
            case SpecUtils::SourceType::Background: if( back == numeric_limits<int>::min() ) back = sample; break;
            case SpecUtils::SourceType::Unknown:    if( unknown == numeric_limits<int>::min() ) unknown = sample; break;
            default: break;
          }
        }//for( loop over samples )

        int chosen = *begin(all_samples);
        if( different_file && (fore != numeric_limits<int>::min()) )
          chosen = fore;
        else if( unknown != numeric_limits<int>::min() )
          chosen = unknown;
        else if( back != numeric_limits<int>::min() )
          chosen = back;
        selected.insert( chosen );
      }else
      {
        selected = all_samples;
      }
      break;
    }//case Background

    case NumSpecType:
      break;
  }//switch( type )

  if( !saved_display.empty() )
  {
    selected = saved_display;
    for( const int sample : selected )
      background.erase( sample );
  }//if( show what was displayed when the file was saved )

  if( note_sample_choice )
  {
    const string shown = (selected.size() == 1) ? ("sample " + std::to_string( *begin(selected) ))
                                                : (std::to_string( selected.size() ) + " samples");
    messages.push_back( "File has " + std::to_string(all_samples.size()) + " samples; showing " + shown
                        + (saved_display.empty() ? "." : ", as when it was saved.") );
  }//if( note_sample_choice )

  if( type == Foreground )
  {
    // A background/secondary showing samples of the previous foreground file goes with it
    const shared_ptr<LoadedFile> old_fore = m_slots[Foreground].file;
    for( const SpecType t : { Background, Secondary } )
    {
      if( old_fore && (m_slots[t].file == old_fore) )
        m_slots[t] = Slot();
    }

    // Default detectors: all, unless the file has multiple "virtual" detectors (sums of the real
    //  ones), in which case just those.
    const vector<string> &dets = spec->detector_names();
    vector<string> vds;
    for( const string &d : dets )
    {
      if( SpecUtils::istarts_with( d, "VD" ) )
        vds.push_back( d );
    }
    m_shown_detectors = (vds.size() > 1) ? vds : dets;
  }//if( type == Foreground )

  setSlot( type, file, selected );
  if( !background.empty() )
    setSlot( Background, file, background );

  // Classify the detector once, rather than on every fit
  if( (type == Foreground) && (own_prefs->m_det_type == PeakFitUtils::CoarseResolutionType::Unknown) )
    own_prefs->m_det_type = PeakFitUtils::coarse_det_type( m_slots[Foreground].summed, spec );

  resetDragCaches();

  json answer = state( PartAll, true );
  if( !messages.empty() )
    answer["messages"] = messages;
  return answer;
}//json loadFile( const json &p )


json Session::unload( const json &p )
{
  const SpecType type = type_from_json( p );
  if( type == Foreground )
  {
    for( Slot &s : m_slots )
      s = Slot();
    m_shown_detectors.clear();
  }else
  {
    m_slots[type] = Slot();
  }

  resetDragCaches();
  return state( PartAll, type == Foreground );
}//json unload( const json &p )


void Session::setSlot( const SpecType type, shared_ptr<LoadedFile> file, set<int> samples )
{
  Slot &slot = m_slots[type];
  if( !file )
  {
    slot = Slot();
    return;
  }

  // Drop invalid sample numbers
  const set<int> &valid = file->spec->sample_numbers();
  for( auto iter = begin(samples); iter != end(samples); )
    iter = valid.count( *iter ) ? std::next( iter ) : samples.erase( iter );
  if( samples.empty() )
    throw runtime_error( "No valid sample numbers to display." );

  slot.file = file;
  slot.samples = samples;
  updateAllSummed(); //Background scale factors depend on the foreground
}//void setSlot(...)


vector<string> Session::displayedDetectors( const SpecType type ) const
{
  const Slot &slot = m_slots[type];
  if( !slot.file )
    return {};

  const vector<string> &names = slot.file->spec->detector_names();
  const Slot &fore = m_slots[Foreground];
  const bool shares_selection = fore.file
          && ((slot.file == fore.file) || (fore.file->spec->detector_names() == names));
  if( !shares_selection )
    return names;

  vector<string> answer;
  for( const string &d : m_shown_detectors )
  {
    if( std::find( begin(names), end(names), d ) != end(names) )
      answer.push_back( d );
  }
  return answer;
}//vector<string> displayedDetectors( const SpecType type ) const


void Session::updateSummed( const SpecType type )
{
  Slot &slot = m_slots[type];
  slot.summed.reset();
  if( !slot.file || slot.samples.empty() )
    return;

  const vector<string> dets = displayedDetectors( type );
  if( dets.empty() )
    return;

  try
  {
    slot.summed = slot.file->spec->sum_measurements( slot.samples, dets, nullptr );
  }catch( std::exception &e )
  {
    cerr << "Failed to sum " << typeName(type) << ": " << e.what() << endl;
  }
}//void updateSummed( const SpecType type )


void Session::updateAllSummed()
{
  for( int i = 0; i < NumSpecType; ++i )
    updateSummed( SpecType(i) );
}


set<int> Session::defaultSamplesForType( const SpecType type ) const
{
  // Port of InterSpec::sampleNumbersForTypeFromForegroundFile
  set<int> answer;
  const Slot &fore = m_slots[Foreground];
  if( !fore.file )
    return answer;

  for( const auto &m : fore.file->spec->measurements() )
  {
    switch( m->source_type() )
    {
      case SpecUtils::SourceType::IntrinsicActivity:
      case SpecUtils::SourceType::Calibration:
        if( type == Secondary )
          answer.insert( m->sample_number() );
        break;
      case SpecUtils::SourceType::Background:
        if( type == Background )
          answer.insert( m->sample_number() );
        break;
      case SpecUtils::SourceType::Foreground:
      case SpecUtils::SourceType::Unknown:
        if( type == Foreground )
          answer.insert( m->sample_number() );
        break;
    }
  }//for( loop over measurements )

  return answer;
}//set<int> defaultSamplesForType(...)


json Session::setSamples( const json &p )
{
  const SpecType type = type_from_json( p );
  const set<int> samples = p.at( "samples" ).get<set<int>>();
  const Slot &slot = m_slots[type];
  shared_ptr<LoadedFile> file = slot.file ? slot.file : m_slots[Foreground].file;
  if( !file )
    throw runtime_error( "No file loaded." );

  setSlot( type, file, samples );
  resetDragCaches();
  return state( PartSpectra | PartTimeChart | PartPeaks | PartEnergyCal | PartFiles );
}//json setSamples( const json &p )


json Session::stepSample( const json &p )
{
  const SpecType type = type_from_json( p );
  const int delta = p.value( "delta", 1 );
  const Slot &slot = m_slots[type];
  if( !slot.file )
    throw runtime_error( "No file loaded." );

  const vector<int> all( begin(slot.file->spec->sample_numbers()), end(slot.file->spec->sample_numbers()) );
  if( all.empty() )
    throw runtime_error( "File has no samples." );

  const int current = slot.samples.empty() ? all.front() : *begin(slot.samples);
  const auto pos = std::lower_bound( begin(all), end(all), current );
  int index = static_cast<int>( pos - begin(all) );
  const int n = static_cast<int>( all.size() );
  index = ((index + delta) % n + n) % n;

  setSlot( type, slot.file, { all[index] } );
  resetDragCaches();
  return state( PartSpectra | PartTimeChart | PartPeaks | PartEnergyCal | PartFiles );
}//json stepSample( const json &p )


json Session::timeDrag( const json &p )
{
  // Port of InterSpec::timeChartDragged
  const Slot &fore = m_slots[Foreground];
  if( !fore.file )
    throw runtime_error( "No foreground loaded." );

  const int first = p.at( "first" ).get<int>();
  const int last = p.at( "last" ).get<int>();
  const int mods = p.value( "mods", 0 );

  enum class Action { Change, Add, Remove };
  const Action action = (mods & 0x1) ? Action::Add : ((mods & 0x2) ? Action::Remove : Action::Change);
  const SpecType type = (mods & 0x4) ? Background : ((mods & 0x8) ? Secondary : Foreground);

  const set<int> &all = fore.file->spec->sample_numbers();
  const auto start_iter = all.find( std::min(first,last) );
  auto end_iter = all.find( std::max(first,last) );
  if( (start_iter == end(all)) || (end_iter == end(all)) )
    throw runtime_error( "Invalid sample numbers from time chart." );
  ++end_iter;
  const set<int> interaction( start_iter, end_iter );

  const Slot &slot = m_slots[type];
  if( slot.file != fore.file )
  {
    if( action != Action::Remove )
      setSlot( type, fore.file, interaction );
  }else
  {
    set<int> samples = slot.samples;
    switch( action )
    {
      case Action::Change:
        samples = interaction;
        break;
      case Action::Add:
        samples.insert( begin(interaction), end(interaction) );
        break;
      case Action::Remove:
        for( const int s : interaction )
          samples.erase( s );
        break;
    }//switch( action )

    if( !samples.empty() )
      setSlot( type, fore.file, samples );
    else if( type == Foreground )
      setSlot( type, fore.file, defaultSamplesForType( Foreground ) );
    else
      setSlot( type, nullptr, {} );
  }//if( different file ) / else

  resetDragCaches();
  return state( PartSpectra | PartTimeChart | PartPeaks | PartEnergyCal | PartFiles );
}//json timeDrag( const json &p )


json Session::setDetectors( const json &p )
{
  m_shown_detectors = p.at( "shown" ).get<vector<string>>();
  updateAllSummed();
  resetDragCaches();
  return state( PartSpectra | PartTimeChart | PartPeaks | PartEnergyCal | PartDetectors );
}//json setDetectors( const json &p )


json Session::filesJson() const
{
  json slots = json::object();
  for( int i = 0; i < NumSpecType; ++i )
  {
    const Slot &slot = m_slots[i];
    if( !slot.file )
    {
      slots[typeName(SpecType(i))] = nullptr;
      continue;
    }

    const SpecUtils::SpecFile &spec = *slot.file->spec;
    json s;
    s["fileId"] = slot.file->id;
    s["name"] = slot.file->name;
    s["samples"] = slot.samples;
    s["numSamples"] = spec.sample_numbers().size();
    s["firstSample"] = spec.sample_numbers().empty() ? 0 : *begin(spec.sample_numbers());
    s["lastSample"] = spec.sample_numbers().empty() ? 0 : *spec.sample_numbers().rbegin();
    s["passthrough"] = spec.passthrough();
    s["sameFileAsForeground"] = (slot.file == m_slots[Foreground].file);
    if( slot.summed )
    {
      s["liveTime"] = slot.summed->live_time();
      s["realTime"] = slot.summed->real_time();
      s["numChannels"] = slot.summed->num_gamma_channels();
    }
    s["instrument"] = SpecUtils::trim_copy( spec.manufacturer() + " " + spec.instrument_model() );
    slots[typeName(SpecType(i))] = s;
  }//for( loop over slots )

  return slots;
}//json filesJson() const


json Session::detectorsJson() const
{
  json answer;
  answer["all"] = json::array();
  answer["shown"] = m_shown_detectors;

  const Slot &fore = m_slots[Foreground];
  if( !fore.file )
    return answer;

  const SpecUtils::SpecFile &spec = *fore.file->spec;
  const vector<string> &gamma = spec.gamma_detector_names();
  for( const string &name : spec.detector_names() )
  {
    const bool has_gamma = (std::find( begin(gamma), end(gamma), name ) != end(gamma));
    answer["all"].push_back( { {"name", name}, {"neutronOnly", !has_gamma} } );
  }

  return answer;
}//json detectorsJson() const


json Session::spectraJson() const
{
  json answer = json::array();
  const Slot &fore = m_slots[Foreground];
  const bool have_back = !!m_slots[Background].summed;

  for( int i = 0; i < NumSpecType; ++i )
  {
    const SpecType type = SpecType(i);
    const Slot &slot = m_slots[i];
    if( !slot.summed )
      continue;

    D3SpectrumExport::D3SpectrumOptions options;
    options.spectrum_type = to_specutils( type );
    options.title = slot.file->name;
    if( slot.file->spec->sample_numbers().size() > 1 )
    {
      if( slot.samples.size() == 1 )
        options.title += " (sample " + std::to_string( *begin(slot.samples) ) + ")";
      else
        options.title += " (" + std::to_string( slot.samples.size() ) + " samples)";
    }

    // Live-time normalize background and secondary to the foreground
    if( (type != Foreground) && fore.summed )
    {
      const double fore_lt = fore.summed->live_time();
      const double lt = slot.summed->live_time();
      if( (fore_lt > 0.0) && (lt > 0.0) )
        options.display_scale_factor = fore_lt / lt;
    }

    if( type == Foreground )
    {
      const PeakDeque *peaks = foregroundPeaksConst();
      if( peaks && !peaks->empty() )
        options.peaks_json = peaksToChartJson( vector<shared_ptr<const PeakDef>>( begin(*peaks), end(*peaks) ) );
    }

    const int specID = i;
    const int backID = ((type != Background) && have_back) ? static_cast<int>(Background) : -1;

    stringstream strm;
    D3SpectrumExport::write_spectrum_data_js( strm, *slot.summed, options, specID, backID );
    answer.push_back( json::parse( strm.str() ) );
  }//for( loop over spectrum types )

  return answer;
}//json spectraJson() const


json Session::timeChartJson() const
{
  // Trimmed port of D3TimeChart::setDataToClient - only for passthrough foregrounds.
  const Slot &fore = m_slots[Foreground];
  if( !fore.file || !fore.file->spec->passthrough() )
    return nullptr;

  const SpecUtils::SpecFile &spec = *fore.file->spec;
  const vector<string> dets = displayedDetectors( Foreground );
  if( dets.empty() )
    return nullptr;

  vector<double> real_times, gamma_counts, live_times, neutron_counts, neutron_live_times;
  vector<int> sample_numbers, source_types;
  vector<json> start_times;
  int64_t start_offset = 0;
  bool any_start_time = false, all_unknown = true, any_gamma = false, any_neutron = false;

  const auto epoch = std::chrono::system_clock::from_time_t( 0 );

  for( const int sample : spec.sample_numbers() )
  {
    double rt = 0.0, gammas = 0.0, lt = 0.0, neuts = 0.0, nlt = 0.0;
    bool have_data = false, have_gamma = false, have_neut = false;
    SpecUtils::time_point_t start{};
    SpecUtils::SourceType st = SpecUtils::SourceType::Unknown;

    for( const string &det : dets )
    {
      const auto m = spec.measurement( sample, det );
      if( !m )
        continue;

      have_data = true;
      rt = std::max( rt, static_cast<double>(m->real_time()) );
      if( SpecUtils::is_special(start) )
        start = m->start_time();
      if( m->source_type() != SpecUtils::SourceType::Unknown )
        st = m->source_type();

      if( m->num_gamma_channels() )
      {
        have_gamma = true;
        gammas += m->gamma_count_sum();
        lt += m->live_time();
      }

      if( m->contained_neutron() )
      {
        have_neut = true;
        neuts += m->neutron_counts_sum();
        nlt += (m->neutron_live_time() > 0.0f) ? m->neutron_live_time() : m->real_time();
      }
    }//for( loop over detectors )

    if( !have_data )
      continue;

    sample_numbers.push_back( sample );
    real_times.push_back( rt );
    source_types.push_back( static_cast<int>(st) );
    all_unknown = all_unknown && (st == SpecUtils::SourceType::Unknown);
    any_gamma |= have_gamma;
    any_neutron |= have_neut;
    gamma_counts.push_back( have_gamma ? gammas : std::numeric_limits<double>::quiet_NaN() );
    live_times.push_back( have_gamma ? lt : std::numeric_limits<double>::quiet_NaN() );
    neutron_counts.push_back( have_neut ? neuts : std::numeric_limits<double>::quiet_NaN() );
    neutron_live_times.push_back( have_neut ? nlt : std::numeric_limits<double>::quiet_NaN() );

    if( SpecUtils::is_special(start) )
    {
      start_times.push_back( nullptr );
    }else
    {
      any_start_time = true;
      const int64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>( start - epoch ).count();
      if( (start_offset <= 0) && (ms > 0) )
        start_offset = ms;
      start_times.push_back( ms - start_offset );
    }
  }//for( loop over samples )

  if( sample_numbers.empty() )
    return nullptr;

  const auto to_json_arr = []( const vector<double> &v ) -> json {
    json arr = json::array();
    for( const double d : v )
      arr.push_back( std::isfinite(d) ? json(d) : json(nullptr) );
    return arr;
  };

  json answer;
  answer["realTimes"] = real_times;
  if( any_start_time )
  {
    answer["startTimeOffset"] = start_offset;
    answer["startTimes"] = start_times;
  }
  answer["sampleNumbers"] = sample_numbers;
  if( !all_unknown )
    answer["sourceTypes"] = source_types;
  answer["isCountsRatio"] = false;

  if( any_gamma )
    answer["gammaCounts"] = json::array( { { {"detName", ""}, {"color", "#cfced2"},
                            {"counts", to_json_arr(gamma_counts)}, {"liveTimes", to_json_arr(live_times)} } } );
  if( any_neutron )
    answer["neutronCounts"] = json::array( { { {"detName", ""}, {"color", "#cfced2"},
                            {"counts", to_json_arr(neutron_counts)}, {"liveTimes", to_json_arr(neutron_live_times)} } } );

  // Occupied sample ranges (only if not every sample is occupied)
  set<int> occupied;
  for( const int sample : spec.sample_numbers() )
  {
    for( const auto &m : spec.sample_measurements( sample ) )
    {
      if( m->occupied() == SpecUtils::OccupancyStatus::Occupied )
        occupied.insert( sample );
    }
  }
  if( !occupied.empty() && (occupied != spec.sample_numbers()) )
  {
    json occ = json::array();
    for( const auto &run : sample_runs( spec.sample_numbers(), occupied ) )
      occ.push_back( { {"startSample", run.first}, {"endSample", run.second}, {"color", "rgb(128,128,128)"} } );
    answer["occupancies"] = occ;
  }

  return answer;
}//json timeChartJson() const


json Session::timeHighlightsJson() const
{
  json answer = json::array();
  const Slot &fore = m_slots[Foreground];
  if( !fore.file || !fore.file->spec->passthrough() )
    return answer;

  const char *colors[NumSpecType] = { "rgba(255,255,0,0.6)", "rgba(0,128,0,0.3)", "rgba(0,255,255,0.3)" };

  for( int i = 0; i < NumSpecType; ++i )
  {
    const Slot &slot = m_slots[i];
    if( (slot.file != fore.file) || slot.samples.empty() )
      continue;

    // Like InterSpec, dont highlight a type showing exactly its default samples
    if( slot.samples == defaultSamplesForType( SpecType(i) ) )
      continue;

    for( const auto &run : sample_runs( fore.file->spec->sample_numbers(), slot.samples ) )
      answer.push_back( { {"startSample", run.first}, {"endSample", run.second},
                          {"fillColor", colors[i]}, {"type", typeName(SpecType(i))} } );
  }//for( loop over slots )

  return answer;
}//json timeHighlightsJson() const


json Session::exportFile( const json &p )
{
  const Slot &fore = m_slots[Foreground];
  if( !fore.file )
    throw runtime_error( "No foreground loaded." );

  const string format = p.value( "format", string("N42-2012") );
  const string path = p.value( "path", string("/tmp/interspec_light_export") );

  struct Fmt { const char *name; SpecUtils::SaveSpectrumAsType type; bool whole_file; };
  const Fmt formats[] = {
    { "N42-2012", SpecUtils::SaveSpectrumAsType::N42_2012, true },
    { "N42-2006", SpecUtils::SaveSpectrumAsType::N42_2006, true },
    { "PCF",      SpecUtils::SaveSpectrumAsType::Pcf,      true },
    { "CSV",      SpecUtils::SaveSpectrumAsType::Csv,      false },
    { "TXT",      SpecUtils::SaveSpectrumAsType::Txt,      false },
    { "CHN",      SpecUtils::SaveSpectrumAsType::Chn,      false },
    { "SPE",      SpecUtils::SaveSpectrumAsType::SpeIaea,  false },
    { "CNF",      SpecUtils::SaveSpectrumAsType::Cnf,      false },
    { "TKA",      SpecUtils::SaveSpectrumAsType::Tka,      false },
    { "SPC",      SpecUtils::SaveSpectrumAsType::SpcBinaryInt, false }
  };

  const Fmt *fmt = nullptr;
  for( const Fmt &f : formats )
  {
    if( SpecUtils::iequals_ascii( format, f.name ) )
      fmt = &f;
  }
  if( !fmt )
    throw runtime_error( "Unknown export format '" + format + "'" );

  const SpecUtils::SpecFile &spec = *fore.file->spec;
  std::ofstream output( path.c_str(), ios::out | ios::binary | ios::trunc );
  if( !output )
    throw runtime_error( "Could not open export file." );

  bool written = false;
#if( LIGHT_PEAK_FILE_IO )
  // Like InterSpec: N42-2012 files hold the peaks of every set of samples, and which samples are
  //  displayed; SPE files hold the displayed peaks.
  const PeakDeque *shown_peaks = foregroundPeaksConst();
  if( fmt->type == SpecUtils::SaveSpectrumAsType::N42_2012 )
  {
    PeakFileIo::write_n42( output, spec, fore.file->peaks, fore.samples, displayedDetectors( Foreground ) );
    written = true;
  }else if( (fmt->type == SpecUtils::SaveSpectrumAsType::SpeIaea) && shown_peaks && !shown_peaks->empty() && fore.summed )
  {
    stringstream spe;
    spec.write( spe, fore.samples, displayedDetectors( Foreground ), fmt->type );
    output << PeakFileIo::add_spe_peaks( spe.str(), *shown_peaks, fore.summed );
    written = true;
  }
#endif

  if( !written && fmt->whole_file )
    spec.write( output, spec.sample_numbers(), spec.detector_names(), fmt->type );
  else if( !written )
    spec.write( output, fore.samples, displayedDetectors( Foreground ), fmt->type );

  if( !output )
    throw runtime_error( "Failed writing export file." );
  output.close();

  json answer;
  answer["path"] = path;
  answer["filename"] = foregroundBaseName() + "." + SpecUtils::suggestedNameEnding( fmt->type );
  return answer;
}//json exportFile( const json &p )


string Session::foregroundBaseName() const
{
  const Slot &fore = m_slots[Foreground];
  string base = fore.file ? fore.file->name : string();
  const size_t dot = base.find_last_of( '.' );
  if( (dot != string::npos) && (dot > 0) )
    base = base.substr( 0, dot );
  return base;
}//string foregroundBaseName() const
