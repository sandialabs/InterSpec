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
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/EnergyCal.h"

#include "Session.h"

using namespace std;
using json = nlohmann::json;

typedef shared_ptr<const SpecUtils::EnergyCalibration> CalPtr;

namespace
{
  const char *cal_type_str( const SpecUtils::EnergyCalType type )
  {
    switch( type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        return "Polynomial";
      case SpecUtils::EnergyCalType::FullRangeFraction:
        return "FullRangeFraction";
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        return "LowerChannelEdge";
      case SpecUtils::EnergyCalType::InvalidEquationType:
        break;
    }
    return "Invalid";
  }//cal_type_str(...)


  bool is_poly_or_frf( const CalPtr &cal )
  {
    if( !cal || !cal->valid() )
      return false;
    switch( cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
        return true;
      default:
        break;
    }
    return false;
  }//is_poly_or_frf(...)


  /** A calibration of the same type as `like`, but with different coefficients. */
  CalPtr make_cal( const SpecUtils::EnergyCalType type, const size_t nchannel,
                   const vector<float> &coefs, const vector<pair<float,float>> &dev_pairs )
  {
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    if( type == SpecUtils::EnergyCalType::FullRangeFraction )
      cal->set_full_range_fraction( nchannel, coefs, dev_pairs );
    else
      cal->set_polynomial( nchannel, coefs, dev_pairs );
    return cal;
  }//make_cal(...)


  /** "foreground", "background", or "secondary", for messages. */
  string type_label( const Session::SpecType type )
  {
    return SpecUtils::to_lower_ascii_copy( Session::typeName( type ) );
  }
}//namespace


CalPtr Session::displayedEnergyCal( const SpecType type ) const
{
  const Slot &slot = m_slots[type];
  if( !slot.file || slot.samples.empty() )
    return nullptr;

  const vector<string> dets = displayedDetectors( type );
  if( dets.empty() )
    return nullptr;

  try
  {
    return slot.file->spec->suggested_sum_energy_calibration( slot.samples, dets );
  }catch( std::exception & )
  {
  }
  return nullptr;
}//displayedEnergyCal(...)


json Session::energyCalJson() const
{
  const CalPtr cal = displayedEnergyCal();
  if( !cal || !cal->valid() )
    return nullptr;

  json answer;
  answer["type"] = cal_type_str( cal->type() );
  answer["coefs"] = cal->coefficients();
  answer["numChannels"] = cal->num_channels();
  answer["numDevPairs"] = cal->deviation_pairs().size();
  answer["editable"] = is_poly_or_frf( cal );
  answer["lowerEnergy"] = cal->lower_energy();
  answer["upperEnergy"] = cal->upper_energy();

  // Revert applies to every loaded file (a CALp file may have changed the background or secondary)
  bool changed = false;
  for( const Slot &slot : m_slots )
  {
    if( !slot.file || changed )
      continue;
    for( const auto &m : slot.file->spec->measurements() )
    {
      const auto pos = slot.file->orig_cals.find( m.get() );
      changed = changed || ((pos != end(slot.file->orig_cals)) && (pos->second != m->energy_calibration()));
    }
  }
  answer["changed"] = changed;

  return answer;
}//json energyCalJson() const


Session::CalChanges Session::calChangesFor( const SpecType type, const CalPtr &disp_prev, const CalPtr &disp_new ) const
{
  // Port of EnergyCalTool::applyCalChange, applied to the one file.
  if( !m_slots[type].file || !disp_prev || !disp_new || !disp_new->valid() )
    throw runtime_error( "Invalid energy calibration." );

  const SpecUtils::SpecFile &spec = *m_slots[type].file->spec;
  const vector<string> dets = displayedDetectors( type );

  const vector<float> &prev_coefs = disp_prev->coefficients();
  const vector<float> &new_coefs = disp_new->coefficients();
  bool offset_only = is_poly_or_frf( disp_prev ) && (disp_prev->type() == disp_new->type())
                     && (prev_coefs.size() == new_coefs.size())
                     && (disp_prev->deviation_pairs() == disp_new->deviation_pairs());
  for( size_t i = 1; offset_only && (i < prev_coefs.size()); ++i )
    offset_only = (prev_coefs[i] == new_coefs[i]);

  // Compute every new calibration (sharing a new one for each distinct old one) before changing
  //  anything, so an invalid result leaves the file untouched.
  map<CalPtr,CalPtr> old_to_new;
  CalChanges changes;

  for( const int sample : spec.sample_numbers() )
  {
    for( const string &det : dets )
    {
      const auto m = spec.measurement( sample, det );
      if( !m || (m->num_gamma_channels() <= 4) )
        continue;

      const CalPtr old_cal = m->energy_calibration();
      if( !old_cal || !old_cal->valid() )
        continue;

      CalPtr &new_cal = old_to_new[old_cal];
      if( !new_cal )
      {
        if( old_cal == disp_prev )
        {
          new_cal = disp_new;
        }else if( offset_only && is_poly_or_frf( old_cal ) )
        {
          vector<float> coefs = old_cal->coefficients();
          coefs[0] += (new_coefs[0] - prev_coefs[0]);
          new_cal = make_cal( old_cal->type(), old_cal->num_channels(), coefs, old_cal->deviation_pairs() );
        }else
        {
          new_cal = EnergyCal::propogate_energy_cal_change( disp_prev, disp_new, old_cal );
        }
      }//if( havent computed new cal yet )

      if( new_cal->num_channels() != m->num_gamma_channels() )
        throw runtime_error( "The new energy calibration is for " + std::to_string( new_cal->num_channels() )
                             + " channels, but detector '" + det + "' of the " + type_label( type ) + " has "
                             + std::to_string( m->num_gamma_channels() ) + "." );

      changes.push_back( { m, new_cal } );
    }//for( loop over detectors )
  }//for( loop over samples )

  return changes;
}//CalChanges calChangesFor(...)


void Session::setCals( const SpecType type, const CalChanges &changes )
{
  const Slot &slot = m_slots[type];
  SpecUtils::SpecFile &spec = *slot.file->spec;
  const vector<string> dets = displayedDetectors( type );

  // Where each displayed calibration goes; if detectors sharing one get different ones, the first wins
  map<CalPtr,CalPtr> old_to_new;
  for( const auto &mc : changes )
  {
    if( std::find( begin(dets), end(dets), mc.first->detector_name() ) != end(dets) )
      old_to_new.emplace( mc.first->energy_calibration(), mc.second );
  }

  // Move peaks of every sample-set of this file whose display calibration is being changed
  map<set<int>,PeakDeque> new_peaks;
  for( const auto &sp : slot.file->peaks )
  {
    if( sp.second.empty() )
      continue;

    const CalPtr old_cal = spec.suggested_sum_energy_calibration( sp.first, dets );
    const auto pos = old_to_new.find( old_cal );
    if( (pos == end(old_to_new)) || (pos->second == old_cal) )
      continue;

    PeakDeque moved = EnergyCal::translatePeaksForCalibrationChange( sp.second, old_cal, pos->second );
    std::sort( begin(moved), end(moved), &PeakDef::lessThanByMeanShrdPtr );
    new_peaks[sp.first] = std::move( moved );
  }//for( loop over peak sets )

  for( const auto &mc : changes )
    spec.set_energy_calibration( mc.second, mc.first );

  for( auto &sp : new_peaks )
    slot.file->peaks[sp.first] = std::move( sp.second );

  updateAllSummed();
  resetDragCaches();
}//void setCals(...)


json Session::setEnergyCal( const json &p )
{
  const CalPtr prev = displayedEnergyCal();
  if( !is_poly_or_frf( prev ) )
    throw runtime_error( "The displayed energy calibration can not be edited." );

  const string type = p.value( "type", string(cal_type_str( prev->type() )) );
  const vector<float> coefs = p.at( "coefs" ).get<vector<float>>();
  if( coefs.size() < 2 )
    throw runtime_error( "At least two coefficients are required." );

  const SpecUtils::EnergyCalType new_type = (type == "FullRangeFraction")
                                              ? SpecUtils::EnergyCalType::FullRangeFraction
                                              : SpecUtils::EnergyCalType::Polynomial;
  CalPtr newcal;
  try
  {
    newcal = make_cal( new_type, prev->num_channels(), coefs, prev->deviation_pairs() );
  }catch( std::exception &e )
  {
    throw runtime_error( string("Invalid energy calibration: ") + e.what() );
  }

  setCals( Foreground, calChangesFor( Foreground, prev, newcal ) );
  return state( PartSpectra | PartPeaks | PartEnergyCal );
}//json setEnergyCal( const json &p )


json Session::convertCalCoefs( const json &p )
{
  // Re-expresses the displayed calibration in the other (Polynomial <--> FRF) form.
  const CalPtr cal = displayedEnergyCal();
  if( !is_poly_or_frf( cal ) )
    throw runtime_error( "Can not convert this calibration." );

  const string to = p.at( "to" ).get<string>();
  const vector<float> coefs = p.at( "coefs" ).get<vector<float>>();
  const string from = p.value( "from", string(cal_type_str( cal->type() )) );

  vector<float> answer = coefs;
  if( (from == "Polynomial") && (to == "FullRangeFraction") )
    answer = SpecUtils::polynomial_coef_to_fullrangefraction( coefs, cal->num_channels() );
  else if( (from == "FullRangeFraction") && (to == "Polynomial") )
    answer = SpecUtils::fullrangefraction_coef_to_polynomial( coefs, cal->num_channels() );

  return { {"coefs", answer} };
}//json convertCalCoefs( const json &p )


json Session::fitEnergyCal( const json &p )
{
  // Port of EnergyCalTool::fitCoefficients
  const CalPtr orig_cal = displayedEnergyCal();
  if( !is_poly_or_frf( orig_cal ) )
    throw runtime_error( "Fitting is only supported for polynomial or full-range-fraction calibrations." );

  const vector<bool> fitfor_in = p.at( "fitFor" ).get<vector<bool>>();
  size_t num_fit = 0, max_order = 0;
  for( size_t i = 0; i < fitfor_in.size(); ++i )
  {
    if( fitfor_in[i] )
    {
      ++num_fit;
      max_order = i;
    }
  }
  if( !num_fit )
    throw runtime_error( "Select at least one coefficient to fit." );

  vector<EnergyCal::RecalPeakInfo> infos;
  for( const auto &peak : foregroundPeaks() )
  {
    if( !peak->useForEnergyCalibration() )
      continue;

    EnergyCal::RecalPeakInfo info;
    info.peakMean = peak->mean();
    info.peakMeanUncert = std::max( peak->meanUncert(), 0.25 );
    if( !std::isfinite( info.peakMeanUncert ) )
      info.peakMeanUncert = 0.5;
    info.photopeakEnergy = peak->gammaParticleEnergy();
    info.peakMeanBinNumber = orig_cal->channel_for_energy( peak->mean() );
    infos.push_back( info );
  }//for( loop over peaks )

  if( num_fit > infos.size() )
    throw runtime_error( "Fitting " + std::to_string(num_fit) + " coefficients needs at least that many"
                         " peaks with an assigned source (and 'use for cal' checked); have "
                         + std::to_string(infos.size()) + "." );

  const size_t nchannel = orig_cal->num_channels();
  const auto &devpairs = orig_cal->deviation_pairs();
  const size_t eqn_order = std::max( orig_cal->coefficients().size(), max_order + 1 );
  vector<bool> fitfor( eqn_order, false );
  for( size_t i = 0; i < fitfor_in.size() && (i < eqn_order); ++i )
    fitfor[i] = fitfor_in[i];

  vector<float> coefs = orig_cal->coefficients(), uncerts;
  coefs.resize( eqn_order, 0.0f );

  const bool is_frf = (orig_cal->type() == SpecUtils::EnergyCalType::FullRangeFraction);
  CalPtr answer;
  try
  {
    if( is_frf )
      EnergyCal::fit_energy_cal_frf( infos, fitfor, nchannel, devpairs, coefs, uncerts );
    else
      EnergyCal::fit_energy_cal_poly( infos, fitfor, nchannel, devpairs, coefs, uncerts );
    answer = make_cal( orig_cal->type(), nchannel, coefs, devpairs );
  }catch( std::exception & )
  {
    // Linear least squares failed; fall back to a Ceres fit
    EnergyCal::EnergyCalCeresFitSetup setup;
    setup.cal_type = orig_cal->type();
    setup.num_channels = nchannel;
    setup.fitfor = fitfor;
    setup.starting_coefs = orig_cal->coefficients();
    setup.starting_coefs.resize( eqn_order, 0.0f );
    setup.dev_pairs = devpairs;
    const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( infos, setup );
    answer = make_cal( orig_cal->type(), nchannel, result.coefs, devpairs );
  }//try / catch

  if( !answer || !answer->valid() )
    throw runtime_error( "Energy calibration fit failed." );

  double pre_dev = 0.0, post_dev = 0.0;
  for( const EnergyCal::RecalPeakInfo &info : infos )
  {
    pre_dev += fabs( info.peakMean - info.photopeakEnergy );
    post_dev += fabs( answer->energy_for_channel( info.peakMeanBinNumber ) - info.photopeakEnergy );
  }

  setCals( Foreground, calChangesFor( Foreground, orig_cal, answer ) );

  char msg[256];
  snprintf( msg, sizeof(msg), "Fit %i coefficient(s) using %i peaks; mean |peak - expected| went from %.3g to %.3g keV.",
            static_cast<int>(num_fit), static_cast<int>(infos.size()),
            pre_dev/infos.size(), post_dev/infos.size() );

  json result = state( PartSpectra | PartPeaks | PartEnergyCal );
  result["message"] = msg;
  return result;
}//json fitEnergyCal( const json &p )


json Session::revertEnergyCal( const json & )
{
  if( !m_slots[Foreground].file )
    throw runtime_error( "No foreground loaded." );

  // Every loaded file goes back to its calibrations as loaded (a CALp file may have changed the
  //  background or secondary too); a file shown in several slots is done once.
  set<const LoadedFile *> done;
  for( int i = 0; i < NumSpecType; ++i )
  {
    const Slot &slot = m_slots[i];
    if( !slot.file || !done.insert( slot.file.get() ).second )
      continue;

    CalChanges changes;
    for( const auto &m : slot.file->spec->measurements() )
    {
      const auto pos = slot.file->orig_cals.find( m.get() );
      if( (pos != end(slot.file->orig_cals)) && pos->second && (pos->second != m->energy_calibration()) )
        changes.push_back( { m, pos->second } );
    }

    if( !changes.empty() )
      setCals( SpecType(i), changes );
  }//for( loop over slots )

  return state( PartSpectra | PartPeaks | PartEnergyCal );
}//json revertEnergyCal( const json & )


json Session::exportCALp( const json &p )
{
  // Port of EnergyCalTool's CALpDownloadResource: the calibration of each gamma detector of the
  //  foreground file, preferring the displayed samples and detectors.  Unlike InterSpec, hidden
  //  detectors are written too, since a per-detector CALp needs every detector it is applied to.
  const Slot &fore = m_slots[Foreground];
  if( !fore.file )
    throw runtime_error( "No foreground loaded." );

  const SpecUtils::SpecFile &spec = *fore.file->spec;
  const string path = p.value( "path", string("/tmp/interspec_light.CALp") );
  std::ofstream output( path.c_str(), ios::out | ios::binary | ios::trunc );
  if( !output )
    throw runtime_error( "Could not open CALp file." );

  // Detector names are only written when needed to tell the calibrations apart
  const bool write_names = (spec.gamma_detector_names().size() > 1);
  set<string> written;
  const auto write_cal = [&]( const int sample, const string &det ){
    if( written.count( det ) )
      return;
    const shared_ptr<const SpecUtils::Measurement> m = spec.measurement( sample, det );
    const CalPtr cal = m ? m->energy_calibration() : nullptr;
    if( cal && cal->valid() && (cal->num_channels() >= 3)
        && SpecUtils::write_CALp_file( output, cal, write_names ? det : string() ) )
      written.insert( det );
  };

  const vector<string> dets = displayedDetectors( Foreground );
  for( const int sample : fore.samples )
  {
    for( const string &det : dets )
      write_cal( sample, det );
  }
  for( const int sample : spec.sample_numbers() )
  {
    for( const string &det : spec.gamma_detector_names() )
      write_cal( sample, det );
  }

  if( written.empty() )
    throw runtime_error( "There is no energy calibration to export." );
  if( !output )
    throw runtime_error( "Failed writing CALp file." );
  output.close();

  return { {"path", path}, {"filename", foregroundBaseName() + ".CALp"} };
}//json exportCALp( const json &p )


Session::CalChanges Session::calpChangesFor( const SpecType type, const string &calp, string &warning ) const
{
  // Port of SpecMeasManager::handleCALpFile and EnergyCalTool::applyCALpEnergyCal
  const CalPtr disp_cal = displayedEnergyCal( type );
  if( !disp_cal || !disp_cal->valid() )
    throw runtime_error( "The " + type_label( type ) + " has no energy calibration to change." );

  const Slot &slot = m_slots[type];
  const SpecUtils::SpecFile &spec = *slot.file->spec;
  const vector<string> &det_names = spec.detector_names();
  const size_t nchannel = disp_cal->num_channels();

  // A CALp file holds one calibration, or one for each detector (named)
  map<string,CalPtr> det_to_cal;
  istringstream input( calp );
  while( input.good() )
  {
    const std::streamoff entry_start = input.tellg();
    string name;
    CalPtr cal;
    try
    {
      cal = SpecUtils::energy_cal_from_CALp_file( input, nchannel, name );
    }catch( std::exception &e )
    {
      // Junk after the last calibration is ignored, but not a calibration that can not be read
      const size_t pos = static_cast<size_t>( std::max( entry_start, std::streamoff(0) ) );
      if( det_to_cal.empty() || SpecUtils::icontains( calp.substr( std::min( pos, calp.size() ) ), "CALp File" ) )
        throw runtime_error( "Not a valid CALp file (calibration " + std::to_string( det_to_cal.size() + 1 )
                             + "): " + e.what() );
      break;
    }

    if( !cal )
      break;

    // A named detector may have a different number of channels than the displayed spectrum
    if( is_poly_or_frf( cal ) && (std::find( begin(det_names), end(det_names), name ) != end(det_names)) )
    {
      for( const int sample : slot.samples )
      {
        const shared_ptr<const SpecUtils::Measurement> m = spec.measurement( sample, name );
        const size_t n = m ? m->num_gamma_channels() : size_t(0);
        if( n <= 3 )
          continue;

        try
        {
          if( n != nchannel )
            cal = make_cal( cal->type(), n, cal->coefficients(), cal->deviation_pairs() );
        }catch( std::exception & )
        {
          // The channel-count checks below report this
        }
        break;
      }//for( loop over samples )
    }//if( a named detector )

    det_to_cal[name] = cal;
  }//while( input.good() )

  if( det_to_cal.empty() )
    throw runtime_error( "Not a valid CALp file." );

  if( det_to_cal.size() == 1 )
  {
    // Like InterSpec, one calibration goes to all shown detectors: directly where they have the
    //  displayed calibration, and propagated as a change where they have another.
    const string &det = begin(det_to_cal)->first;
    const CalPtr &cal = begin(det_to_cal)->second;
    const bool det_in_file = (std::find( begin(det_names), end(det_names), det ) != end(det_names));
    const string show_one = det_in_file ? ("; show only detector '" + det + "' to apply it.")
                                        : string("; show one detector at a time to apply it.");
    if( cal->num_channels() != nchannel )
      throw runtime_error( "The CALp calibration for detector '" + det + "' is for "
                           + std::to_string( cal->num_channels() ) + " channels, but the displayed "
                           + type_label( type ) + " has " + std::to_string( nchannel ) + show_one );

    const CalChanges changes = calChangesFor( type, disp_cal, cal );
    set<string> dets;
    bool propagated = false;
    for( const auto &mc : changes )
    {
      dets.insert( mc.first->detector_name() );
      propagated = propagated || (mc.second != cal);
    }

    // A named calibration is only meant for its own detector
    if( propagated && !det.empty() && (dets.size() > 1) )
      throw runtime_error( "The CALp file only has a calibration for detector '" + det + "', but the shown"
                           " detectors of the " + type_label( type ) + " have different calibrations" + show_one );

    if( propagated && !cal->deviation_pairs().empty() )
      warning = "The shown detectors of the " + type_label( type ) + " have different calibrations, so only"
                " those with the displayed one got the CALp file's deviation pairs; to apply them to the"
                " others, show one detector at a time.";

    return changes;
  }//if( one calibration )

  // Per-detector calibrations: each displayed gamma detector must have one
  CalChanges changes;
  string missing;
  for( const string &det : displayedDetectors( type ) )
  {
    const auto pos = det_to_cal.find( det );
    for( const int sample : spec.sample_numbers() )
    {
      const shared_ptr<const SpecUtils::Measurement> m = spec.measurement( sample, det );
      if( !m || (m->num_gamma_channels() <= 4) )
        continue;

      if( (pos == end(det_to_cal)) && det.empty() )
        throw runtime_error( "The calibrations in the CALp file are for named detectors, but the detector of the "
                             + type_label( type ) + " has no name." );

      if( pos == end(det_to_cal) )
      {
        missing += (missing.empty() ? "'" : ", '") + det + "'";
        break;
      }

      if( pos->second->num_channels() != m->num_gamma_channels() )
        throw runtime_error( "The CALp calibration for detector '" + det + "' is for "
                             + std::to_string( pos->second->num_channels() ) + " channels, but the "
                             + type_label( type ) + " has " + std::to_string( m->num_gamma_channels() ) + "." );

      changes.push_back( { m, pos->second } );
    }//for( loop over samples )
  }//for( loop over displayed detectors )

  if( !missing.empty() )
    throw runtime_error( "The CALp file has no calibration for detector(s) " + missing + " of the "
                         + type_label( type ) + "." );

  return changes;
}//CalChanges calpChangesFor(...)


json Session::importCALp( const json &p )
{
  const string path = p.at( "path" ).get<string>();
  const string name = p.value( "name", string("the CALp file") );
  if( !m_slots[Foreground].file )
    throw runtime_error( "Load a foreground spectrum before applying a CALp file." );

  std::ifstream input( path.c_str(), ios::in | ios::binary );
  if( !input )
    throw runtime_error( "Could not open '" + name + "'." );
  const string calp( (std::istreambuf_iterator<char>( input )), std::istreambuf_iterator<char>() );

  set<SpecType> wanted;
  for( const json &t : p.value( "types", json::array( { "FOREGROUND" } ) ) )
    wanted.insert( typeFromName( t.get<string>() ) );

  // Check every file before changing any; a file shown in several slots is changed once.  Like
  //  InterSpec, they must all have the same number of channels, which a CALp file's coefficients assume.
  vector<pair<SpecType,CalChanges>> todo;
  string warnings;
  set<const LoadedFile *> files;
  for( int i = 0; i < NumSpecType; ++i )
  {
    const SpecType type = SpecType(i);
    const Slot &slot = m_slots[i];
    if( !wanted.count( type ) || !slot.file || !files.insert( slot.file.get() ).second )
      continue;

    const CalPtr disp_cal = displayedEnergyCal( type );
    const CalPtr first_cal = todo.empty() ? nullptr : displayedEnergyCal( todo.front().first );
    if( disp_cal && first_cal && (disp_cal->num_channels() != first_cal->num_channels()) )
      throw runtime_error( "The " + type_label( type ) + " has " + std::to_string( disp_cal->num_channels() )
                           + " channels, but the " + type_label( todo.front().first ) + " has "
                           + std::to_string( first_cal->num_channels() ) + "; a CALp file can only be applied"
                           " to spectra with the same number of channels." );

    string warning;
    todo.push_back( { type, calpChangesFor( type, calp, warning ) } );
    warnings += (warning.empty() ? "" : " ") + warning;
  }//for( loop over spectrum types )

  if( todo.empty() )
    throw runtime_error( "There is no spectrum to apply the CALp file to." );

  string applied;
  for( size_t i = 0; i < todo.size(); ++i )
  {
    setCals( todo[i].first, todo[i].second );
    applied += string( !i ? "" : ((i + 1) == todo.size()) ? " and " : ", " ) + type_label( todo[i].first );
  }

  json answer = state( PartSpectra | PartPeaks | PartEnergyCal );
  answer["message"] = "Applied the energy calibration from " + name + " to the " + applied + "." + warnings;
  return answer;
}//json importCALp( const json &p )
