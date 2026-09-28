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
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/SpecFile.h"
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
}//namespace


CalPtr Session::displayedEnergyCal() const
{
  const Slot &fore = m_slots[Foreground];
  if( !fore.file || fore.samples.empty() )
    return nullptr;

  const vector<string> dets = displayedDetectors( Foreground );
  if( dets.empty() )
    return nullptr;

  try
  {
    return fore.file->spec->suggested_sum_energy_calibration( fore.samples, dets );
  }catch( std::exception & )
  {
  }
  return nullptr;
}//displayedEnergyCal()


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

  bool changed = false;
  const Slot &fore = m_slots[Foreground];
  for( const auto &m : fore.file->spec->measurements() )
  {
    const auto pos = fore.file->orig_cals.find( m.get() );
    changed = changed || ((pos != end(fore.file->orig_cals)) && (pos->second != m->energy_calibration()));
  }
  answer["changed"] = changed;

  return answer;
}//json energyCalJson() const


void Session::applyCalChange( const CalPtr &disp_prev, const CalPtr &disp_new )
{
  // Port of EnergyCalTool::applyCalChange, applied to the foreground file only.
  if( !disp_prev || !disp_new || !disp_new->valid() )
    throw runtime_error( "Invalid energy calibration." );

  const Slot &fore = m_slots[Foreground];
  SpecUtils::SpecFile &spec = *fore.file->spec;
  const vector<string> dets = displayedDetectors( Foreground );

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
  vector<pair<shared_ptr<const SpecUtils::Measurement>,CalPtr>> to_set;

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

      to_set.push_back( { m, new_cal } );
    }//for( loop over detectors )
  }//for( loop over samples )

  // Move peaks of every sample-set of this file whose display calibration is being changed
  map<set<int>,PeakDeque> new_peaks;
  for( const auto &sp : fore.file->peaks )
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

  for( const auto &mc : to_set )
    spec.set_energy_calibration( mc.second, mc.first );

  for( auto &sp : new_peaks )
    fore.file->peaks[sp.first] = std::move( sp.second );

  updateAllSummed();
  resetDragCaches();
}//void applyCalChange(...)


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

  applyCalChange( prev, newcal );
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

  applyCalChange( orig_cal, answer );

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
  const Slot &fore = m_slots[Foreground];
  if( !fore.file )
    throw runtime_error( "No foreground loaded." );

  SpecUtils::SpecFile &spec = *fore.file->spec;
  const vector<string> all_dets = spec.detector_names();

  // Map each current calibration back to the original; peaks move with the display calibration.
  map<CalPtr,CalPtr> cur_to_orig;
  vector<pair<shared_ptr<const SpecUtils::Measurement>,CalPtr>> to_set;
  for( const auto &m : spec.measurements() )
  {
    const auto pos = fore.file->orig_cals.find( m.get() );
    if( (pos == end(fore.file->orig_cals)) || !pos->second || (pos->second == m->energy_calibration()) )
      continue;
    cur_to_orig[m->energy_calibration()] = pos->second;
    to_set.push_back( { m, pos->second } );
  }

  const vector<string> dets = displayedDetectors( Foreground );
  map<set<int>,PeakDeque> new_peaks;
  for( const auto &sp : fore.file->peaks )
  {
    if( sp.second.empty() )
      continue;
    const CalPtr cur = spec.suggested_sum_energy_calibration( sp.first, dets );
    const auto pos = cur_to_orig.find( cur );
    if( pos == end(cur_to_orig) )
      continue;
    PeakDeque moved = EnergyCal::translatePeaksForCalibrationChange( sp.second, cur, pos->second );
    std::sort( begin(moved), end(moved), &PeakDef::lessThanByMeanShrdPtr );
    new_peaks[sp.first] = std::move( moved );
  }

  for( const auto &mc : to_set )
    spec.set_energy_calibration( mc.second, mc.first );
  for( auto &sp : new_peaks )
    fore.file->peaks[sp.first] = std::move( sp.second );

  updateAllSummed();
  resetDragCaches();
  return state( PartSpectra | PartPeaks | PartEnergyCal );
}//json revertEnergyCal( const json & )
