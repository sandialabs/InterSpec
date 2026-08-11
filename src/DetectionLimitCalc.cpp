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

#include <cmath>
#include <mutex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cassert>
#include <limits>
#include <numeric>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>


#include "SpecUtils/SpecFile.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

#if( PERFORM_DEVELOPER_CHECKS )
#include <boost/tokenizer.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#endif //PERFORM_DEVELOPER_CHECKS

using namespace std;

namespace DetectionLimitCalc
{
 
#if( PERFORM_DEVELOPER_CHECKS )
  
float round_to_nearest_channel_edge( const float energy, const shared_ptr<const SpecUtils::Measurement> &m )
{
  if( !m )
    return energy;
    
  auto cal = m->energy_calibration();
    
  if( !cal || !cal->valid() )
    return energy;
    
  const double channel = std::max( 0.0, cal->channel_for_energy(energy) ); //std::max isnt necessary, but JIC
  const double whole = std::floor(channel);
  const double frac = channel - whole;
  
  if( (frac >= 0.5) && ((whole+1) < cal->num_channels()) )
    return static_cast<float>( cal->energy_for_channel(whole + 1.0) );
  
  return static_cast<float>( cal->energy_for_channel(whole) );
}//float round_to_nearest_channel_edge( float energy )
  
  
DetectionLimitCalc::CurrieMdaInput currie_input( const float energy,
                                                const shared_ptr<const SpecUtils::Measurement> &m,
                                                shared_ptr<const DetectorPeakResponse> &det,
                                               const double detection_probability )
{
  if( !m || (m->num_gamma_channels() < 16) || !det || !det->isValid() || !det->hasResolutionInfo() )
    throw runtime_error( "No measurement or no DRF." );
  
  const size_t nsidebin = 4;
  //const float num_fwhm = 2.5;
  const float nfwhm = 1.25; // recommended by ISO 11929:2010, could instead use 1.19
  const double confidence_level = detection_probability;
  
  const float fwhm = det->peakResolutionFWHM( energy );
  const float roi_lower_energy = round_to_nearest_channel_edge( energy - nfwhm*fwhm, m ) + 0.0001f;
  const float roi_upper_energy = round_to_nearest_channel_edge( energy + nfwhm*fwhm, m ) - 0.0001f;
  
  DetectionLimitCalc::CurrieMdaInput input;
  input.spectrum = m;
  input.gamma_energy = energy;
  input.roi_lower_energy = roi_lower_energy;
  input.roi_upper_energy = roi_upper_energy;
  input.num_lower_side_channels = nsidebin;
  input.num_upper_side_channels = nsidebin;
  input.detection_probability = confidence_level;
  input.additional_uncertainty = 0.0f;  // TODO: can we get the DRFs contribution to form this?
  
  return input;
}//CurrieMdaInput currie_input(...)
  
  
void batch_test()
{
  const string base_dir = "/Users/wcjohns/Downloads/MDA_calc_20230718/";
  const string spec_file = SpecUtils::append_path( base_dir, "Livermore_48_hour background 5-26-23.n42" );
  
  SpecMeas meas;
  if( !meas.load_N42_file( spec_file ) )
    throw runtime_error( "Couldn't open '" + spec_file + "'" );
  
  if( meas.num_measurements() != 1 )
    throw runtime_error( "Not exactly one measurement in the file" );
  
  vector<shared_ptr<const SpecUtils::Measurement>> meass = meas.measurements();
  assert( !meass.empty() );
  
  const shared_ptr<const SpecUtils::Measurement> spectrum = meass.empty() ? nullptr : meass[0];
  if( !spectrum || (spectrum->num_gamma_channels() < 128) )
    throw runtime_error( "No spectrum" );
  
  if( !spectrum->energy_calibration() || !spectrum->energy_calibration()->valid() )
    throw runtime_error( "No energy calibration" );
  
  shared_ptr<const DetectorPeakResponse> det = meas.detector();
  if( !det || !det->isValid() || !det->hasResolutionInfo() )
    throw runtime_error( "Invalid detector, or missing resolution information" );
  
  const double detection_probability = 0.95;
  const bool useCuries = false;
  const double distance = 2.54*PhysicalUnits::cm;
  const double shielding_thickness = 2.5*PhysicalUnits::mm;
  const double live_time = spectrum->live_time();
  
  const bool fixed_geom = det->isFixedGeometry();
  
  typedef boost::math::policies::policy<boost::math::policies::digits10<6> > my_pol_6;
  const boost::math::normal_distribution<float,my_pol_6> gaus_dist( 0.0f, 1.0f );
  // Will map 0.8414->1.00023, 0.95->1.64485, 0.975->1.95996, 0.995->2.57583, ...
  const float k = boost::math::quantile( gaus_dist, detection_probability );
  
  /* const string test_nucs[] = { "Co-60" }; */
  
  const string test_nucs[] = {
    "Ag-106", "Ag-106m", "Ag-108", "Ag-108m", "Ag-110", "Ag-110m", "Ag-111", "Al-28", "Al-29",
    "Ar-37", "Ar-41", "Ar-42", "As-74", "As-76", "As-77", "Au-196", "Au-197m", "Au-198", "B-12",
    "Be-10", "Bi-210m", "Bi-210", "Bi-211", "Br-83", "C-11", "C-14", "Ca-41", "Ca-45", "Ca-47",
    "Ca-49", "Cd-107", "Cd-109", "Cd-111", "Cd-113", "Cd-113", "Cd-115", "Cd-115", "Cd-117",
    "Cd-117", "Cl-36", "Cl-38", "Cl-38m", "Co-58", "Co-58m", "Co-60", "Co-60m", "Co-61", "Cr-49",
    "Cr-51", "Cr-55", "Cu-62", "Cu-64", "Cu-66", "Cu-67", "F-18", "F-20", "Fe-53", "Fe-55", "Fe-59",
    "Ge-75", "H-3", "He-6", "Hf-181", "Hf-182", "Hg-197", "Hg-197m", "Hg-199m", "Hg-203", "Hg-205",
    "I-129", "I-131", "K-40", "K-42", "K-43", "Li-8", "Lu-178", "Mg-27", "Mg-28", "Mn-54", "Mn-56",
    "Mo-101", "Mo-93", "Mo-93m", "Mo-99", "N-16", "N-17", "Na-22", "Na-24", "Na-24m", "Na-25",
    "Nb-92", "Nb-92m", "Nb-93m", "Nb-94", "Nb-94m", "Nb-95", "Nb-96", "Ne-23", "Ni-57", "Ni-59",
    "Ni-63", "Ni-65", "Ni-66", "O-15", "O-19", "P-32", "P-33", "P-34", "Pb-203", "Pb-204m",
    "Pb-205", "Pb-209", "Pd-103", "Pd-107", "Pd-107m", "Pd-109", "Pd-109m", "Pd-111m", "Pd-111t",
    "Po-210", "S-35", "S-37", "Sb-122", "Sb-122m", "Sb-124", "Sb-124m1", "Sb-124m2", "Sb-125",
    "Sc-46", "Sc-47", "Sc-48", "Sc-49", "Se-75", "Se-77m", "Se-79", "Se-79m", "Se-81", "Se-81m",
    "Se-83m", "Se-83", "Si-31", "Si-32", "Sn-113", "Sn-113m", "Sn-117m", "Sn-119m", "Sn-121",
    "Sn-121m", "Sn-123", "Sn-123m", "Sn-125m", "Sn-125", "Sn-126", "Sr-89", "Sr-90", "Ta-182",
    "Ta-182m", "Ta-183", "Tc-99", "Te-121", "Te-121m", "Te-123", "Te-123m", "Te-125m", "Te-127",
    "Te-127m", "Te-129m", "Te-129", "Te-131m", "Te-131", "Te-132", "Ti-45", "Ti-51", "V-52",
    "V-53", "W-181", "W-185", "W-185m", "W-187", "W-188", "Y-90", "Y-91", "Zn-65", "Zn-69",
    "Zn-69m", "Zn-71", "Zn-71m", "Zr-89", "Zr-93", "Zr-95", "Zr-95", "Zr-97"
  };
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  assert( matdb );

  const std::shared_ptr<const Material> stainless = matdb->material("stainless-steel NIST");
  assert( stainless );
  
  const auto shield_transmission = [stainless, shielding_thickness]( const double energy ) -> double {
    const double atten_coef = GammaInteractionCalc::transmition_coefficient_material( stainless.get(), energy, shielding_thickness );
    return exp( -1.0*atten_coef );
  };
  
  map<const SandiaDecay::Nuclide *,double> single_peak_sensitivity;
  
  for( const string &nuc_str : test_nucs )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( nuc_str );
    if( !nuc || nuc->isStable() )
    {
      cerr << "Warning: '" << nuc_str << "' is not a valid nuclide - skipping" << endl;
      continue;
    }
    
    // TODO: we could integrate over activation time-frame - the code is probably similar to
    //       nuclide decay during measurement, but we wont for the moment.
    double age = 5.0*nuc->halfLife;
    if( age > 10*SandiaDecay::year )
      age = 10*SandiaDecay::year;
    
    const double parent_act = 1.0E-3*SandiaDecay::curie; //Will get divided out, doesnt matter, as long as not too small.
    
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc, parent_act, age );
    
    const vector<SandiaDecay::EnergyRatePair> gammas
                  = mix.gammas( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
    
    vector<pair<double,double>> energy_activities;
    
    for( const SandiaDecay::EnergyRatePair &erp : gammas )
    {
      if( erp.energy < 10 )
        continue;
      
      try
      {
      const double det_eff = fixed_geom ? det->intrinsicEfficiency(erp.energy)
                                        : det->efficiency( erp.energy, distance );
      const double shield_trans = shield_transmission( erp.energy );
      const double gammas_per_bq_per_sec = erp.numPerSecond / parent_act;
      if( gammas_per_bq_per_sec < 1.0E-16 )
        continue;
      
      const CurrieMdaInput mda_input = currie_input( erp.energy, spectrum, det, detection_probability );
      
      const CurrieMdaResult result = DetectionLimitCalc::currie_mda_calc( mda_input );
      
      const double counts_per_bq = det_eff * shield_trans * gammas_per_bq_per_sec * live_time;
      
      // Since we dont know if we are right next to, or overlapping, a peak or something, we'll
      //  just require an excess of observed counts in the region
      const double decision_excess_counts = k * sqrt( result.peak_region_counts_sum );
      const double decision_act = decision_excess_counts / counts_per_bq;
      
      energy_activities.push_back( {erp.energy, decision_act} );
      
      /*
      const int label_width = 26;
      cout << "For " << nuc->symbol << ", at " << erp.energy << " keV:" << endl;
      
      cout << "\tdecision_act: " << PhysicalUnits::printToBestActivityUnits(decision_act, 4, useCuries) << endl;
      cout << endl;
      
      cout << std::left << std::setw(label_width) << "\tLower region channels:"
           << "[" << result.first_lower_continuum_channel << ", "
           << result.last_lower_continuum_channel << "]"
           << endl;
      cout << std::left << std::setw(label_width) << "\tLower region counts:"
           << SpecUtils::printCompact( result.lower_continuum_counts_sum, 5 )
           << endl;
      cout << std::left << std::setw(label_width) << "\tUpper region channels:"
           << "[" << result.first_upper_continuum_channel << ", "
           << result.last_upper_continuum_channel << "]"
           << endl;
      cout << std::left << std::setw(label_width) << "\tUpper region counts:"
           << SpecUtils::printCompact( result.upper_continuum_counts_sum, 5 ) << endl;
      cout << std::left << std::setw(label_width) << "\tPeak area channels:"
           << "[" << (result.last_lower_continuum_channel + 1) << ", "
           << (result.first_upper_continuum_channel - 1) << "]"
           << endl;
      cout << std::left << std::setw(label_width) << "\tPeak region counts:"
           << SpecUtils::printCompact( result.peak_region_counts_sum, 5 )
           << endl;
      cout << std::left << std::setw(label_width) << "\tPeak region null est.:"
           << SpecUtils::printCompact( result.estimated_peak_continuum_counts, 5 )
           << " +- " << SpecUtils::printCompact( result.estimated_peak_continuum_uncert, 5 )
           << endl;
      cout << std::left << std::setw(label_width) << "\tPeak critical limit:"
           << SpecUtils::printCompact( result.decision_threshold, 5 )
           << endl;
      cout << std::left << std::setw(label_width) << "\tPeak detection limit:"
           << SpecUtils::printCompact( result.detection_limit, 5 )
           << endl;
      cout << endl;
      const double intrinsic_eff = det->intrinsicEfficiency( erp.energy );
      const double geom_eff = det->fractionalSolidAngle( det->detectorDiameter(), distance + det->detectorSetback() );
      cout << std::left << std::setw(label_width) << "\tDetector Intrinsic Eff.:"
           << SpecUtils::printCompact( intrinsic_eff, 5 )
           << endl;
      cout << std::left << std::setw(label_width) << "\tSolid angle fraction:"
           << SpecUtils::printCompact( geom_eff, 5 )
           << endl;
      cout << std::left << std::setw(label_width) << "\tShielding transmission:"
           << SpecUtils::printCompact( shield_trans, 5 )
           << endl;
      cout << std::left << std::setw(label_width) << "\tNuclide branching ratio:"
           << SpecUtils::printCompact( gammas_per_bq_per_sec, 5 )
           << endl;
      cout << endl;
      */
        
      }catch( std::exception &e )
      {
        cerr << "Failed to calc limit for " << nuc_str << " at " << erp.energy << " keV: "
             << e.what() << endl;
      }// try / catch - for an energy
    }//for( const SandiaDecay::EnergyRatePair &erp : gammas )
    
    std::sort( begin(energy_activities), end(energy_activities),
      []( const pair<double,double> &lhs, const pair<double,double> &rhs ) -> bool {
        return lhs.second < rhs.second;
    } );
    
    // We could go through and use multi-peak MDA for the top ~3 or 5 peaks...
    if( !energy_activities.empty() )
    {
      cout << nuc_str << ", "
      << PhysicalUnits::printToBestActivityUnits(energy_activities.front().second, 4, useCuries)
      << ", " << energy_activities.front().first << " keV,"
      << PhysicalUnits::printToBestTimeUnits(nuc->halfLife,4)
      << endl;
      
      single_peak_sensitivity[nuc] = energy_activities.front().second;
    }//if( !energy_activities.empty() )
  }//for( const string &nuc_str : test_nucs )
  
  
  const double wanted_grams = 10;
  const double input_data_mass_grams = 35771;
  const string csv_dir = SpecUtils::append_path( base_dir, "csv_out_all_mats_fluence_318_fast_2" );
  const vector<string> csv_files = SpecUtils::recursive_ls( csv_dir, "35_years_exp.csv" );
  
  cout << "Nuclide, HalfLife (days), MinDetectableAct (bq), Expected Act (bq) per " << wanted_grams << " grams, Timespan detectable for (days), Material" << endl;
  
  for( const string csv_filename : csv_files )
  {
    ifstream file( csv_filename.c_str(), ios::in | ios::binary );
    assert( file.is_open() );
    
    string line;
    while( SpecUtils::safe_get_line( file, line ) )
    {
      if( SpecUtils::istarts_with(line, "CSV") )
      {
        assert( SpecUtils::icontains(line, "(in uCi)") );
      }
      
      if( SpecUtils::istarts_with(line, "Target") || SpecUtils::istarts_with(line, "CSV") )
        continue;
      
      vector<string> fields;
      
      // The reactions have a comma in them, but those fields are quoted, so we'll parse
      //  with just a little bit of care
      typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
      boost::escaped_list_separator<char> separator("\\",",", "\"");
      Tokeniser t( line, separator );
      for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
      {
          fields.push_back( *it );
      }
      
      if( fields.empty() )
        continue;
      
      if( fields.size() != 8 )
      {
        cout << "Line contains " << fields.size() << " fields: " << line << endl;
        continue;
      }
      
      string nuc_str = fields[1];
      SpecUtils::ireplace_all(nuc_str, "+", "" );
      SpecUtils::ireplace_all(nuc_str, "*", "" );
      if( SpecUtils::iends_with(nuc_str, "t") )
        nuc_str = nuc_str.substr(0, nuc_str.size() - 1);
      if( SpecUtils::iends_with(nuc_str, "s") )
        nuc_str = nuc_str.substr(0, nuc_str.size() - 1);
      
      const SandiaDecay::Nuclide * const nuc = db->nuclide( nuc_str );
      if( !nuc )
      {
        cerr << "Failed to get nuc from '" << nuc_str << "'" << endl;
        continue;
      }
      
      const string act_str_uci = fields[4];
      if( act_str_uci == "--" )
        continue;
      
      float csv_activity_uci;
      if( !SpecUtils::parse_float( act_str_uci.c_str(), act_str_uci.size(), csv_activity_uci) )
      {
        cerr << "Failed to parse activity '" << act_str_uci << "'" << endl;
        continue;
      }
      
      csv_activity_uci *= PhysicalUnits::microCi;
      
      const double wanted_mass_activity = wanted_grams * csv_activity_uci / input_data_mass_grams;
      
      if( !single_peak_sensitivity.count(nuc) )
      {
        // Happens if there isnt a gamma
        //cerr << "No activity limit available for " << nuc->symbol << endl;
        continue;
      }
      
      
      /** The decay constant that is defined as
          0.5 = exp( -decay_const*halfLife ), or put another way ln(0.5)/halfLife.
       */
      //min_det_act = wanted_mass_activity * exp(-nuc->decayConstant() * X )
      //min_det_act/wanted_mass_activity = exp(-nuc->decayConstant() * X )
      //ln( min_det_act/wanted_mass_activity ) = -nuc->decayConstant() * X
      //X = ln( min_det_act/wanted_mass_activity ) / -nuc->decayConstant()
      
      const double min_det_act = single_peak_sensitivity[nuc];
      if( wanted_mass_activity > min_det_act )
      {
        const double time_til_min_det = - std::log( min_det_act/wanted_mass_activity ) / nuc->decayConstant();
        
        cout << nuc->symbol << ","
             << (nuc->halfLife / PhysicalUnits::day) << ","
             << min_det_act << ","
             << wanted_mass_activity << ","
             << (time_til_min_det / PhysicalUnits::day) << ","
             << SpecUtils::filename(csv_filename)
             << endl;
      }
    }//while( SpecUtils::safe_get_line( file, line ) )
  }//for( const string csv_filename : csv_files )
  
}//void batch_test()
#endif

CurrieMdaInput::CurrieMdaInput()
  : spectrum(nullptr),
    gamma_energy(0.0f), roi_lower_energy(0.0f), roi_upper_energy(0.0f),
    num_lower_side_channels(0), num_upper_side_channels(0),
    detection_probability(0.0f), additional_uncertainty(0.0f)
{
}


CurrieMdaResult::CurrieMdaResult()
  : first_lower_continuum_channel(0), last_lower_continuum_channel(0), lower_continuum_counts_sum(0.0f),
    first_upper_continuum_channel(0), last_upper_continuum_channel(0), upper_continuum_counts_sum(0.0f),
    first_peak_region_channel(0), last_peak_region_channel(0), peak_region_counts_sum(0.0f),
    continuum_eqn{ 0.0f, 0.0f },
    estimated_peak_continuum_counts(0.0f), estimated_peak_continuum_uncert(0.0f),
    decision_threshold(0.0f), detection_limit(0.0f), source_counts(0.0f),
    lower_limit(0.0f), upper_limit(0.0f)
{
}


#if( PERFORM_DEVELOPER_CHECKS )

template<class T>
bool floats_equiv_enough( const T a, const T b )
{
  // This function checks the arguments are EITHER within 'abs_epsilon' of each other,
  //  OR 'rel_epsilon * max(a,b)' of each other.
  const T abs_epsilon = 1.0E-10f; //arbitrary, could pick like std::numeric_limits<float>::epsilon();
  const T rel_epsilon = 1.0E-5f;
  
  const auto diff = fabs(a - b);
  const auto maxval = std::max(fabs(a),fabs(b));
  
  return (diff < (rel_epsilon*maxval) || (diff < abs_epsilon));
};


void CurrieMdaResult::equal_enough( const CurrieMdaResult &test, const CurrieMdaResult &expected )
{
  //CurrieMdaInput::equal_enough( test.input, expected.input );
  
  vector<string> errors;
  char buffer[512];
  
  if( test.first_lower_continuum_channel != expected.first_lower_continuum_channel )
    errors.push_back( "Test first_lower_continuum_channel ("
                     + std::to_string(test.first_lower_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.first_lower_continuum_channel) );
  
  if( test.last_lower_continuum_channel != expected.last_lower_continuum_channel )
    errors.push_back( "Test last_lower_continuum_channel ("
                     + std::to_string(test.last_lower_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.last_lower_continuum_channel) );
  
  if( !floats_equiv_enough(test.lower_continuum_counts_sum, expected.lower_continuum_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
              "Test lower_continuum_counts_sum (%.6G) does not equal expected (%.6G)",
             test.lower_continuum_counts_sum, expected.lower_continuum_counts_sum );
    errors.push_back( buffer );
  }
  
  if( test.first_upper_continuum_channel != expected.first_upper_continuum_channel )
    errors.push_back( "Test first_upper_continuum_channel ("
                     + std::to_string(test.first_upper_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.first_upper_continuum_channel) );
  
  
  if( test.last_upper_continuum_channel != expected.last_upper_continuum_channel )
    errors.push_back( "Test last_upper_continuum_channel ("
                     + std::to_string(test.last_upper_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.last_upper_continuum_channel) );
  
  if( !floats_equiv_enough(test.upper_continuum_counts_sum, expected.upper_continuum_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test lower_continuum_counts_sum (%.6G) does not equal expected (%.6G)",
             test.upper_continuum_counts_sum, expected.upper_continuum_counts_sum );
    errors.push_back( buffer );
  }
  
  
  if( test.first_peak_region_channel != expected.first_peak_region_channel )
    errors.push_back( "Test first_peak_region_channel ("
                     + std::to_string(test.first_peak_region_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.first_peak_region_channel) );
  
  if( test.last_peak_region_channel != expected.last_peak_region_channel )
    errors.push_back( "Test last_peak_region_channel ("
                     + std::to_string(test.last_peak_region_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.last_peak_region_channel) );


  if( !floats_equiv_enough(test.peak_region_counts_sum, expected.peak_region_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test peak_region_counts_sum (%.6G) does not equal expected (%.6G)",
             test.peak_region_counts_sum, expected.peak_region_counts_sum );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.peak_region_counts_sum, expected.peak_region_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test peak_region_counts_sum (%.6G) does not equal expected (%.6G)",
             test.peak_region_counts_sum, expected.peak_region_counts_sum );
    errors.push_back( buffer );
  }
  
  /*
  for( int i = 0; i < 2; ++i )
  {
    if( !floats_equiv_enough(test.continuum_eqn[i], expected.continuum_eqn[i]) )
    {
      snprintf( buffer, sizeof(buffer),
               "Test continuum_eqn[%i] (%.6G) does not equal expected (%.6G)",
               i, test.continuum_eqn[i], expected.continuum_eqn[i] );
      errors.push_back( buffer );
    }
  }//for( int i = 0; i < 2; ++i )
  */
  
  if( !floats_equiv_enough(test.estimated_peak_continuum_counts, expected.estimated_peak_continuum_counts) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test estimated_peak_continuum_counts (%.6G) does not equal expected (%.6G)",
             test.estimated_peak_continuum_counts, expected.estimated_peak_continuum_counts );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.estimated_peak_continuum_uncert, expected.estimated_peak_continuum_uncert) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test estimated_peak_continuum_uncert (%.6G) does not equal expected (%.6G)",
             test.estimated_peak_continuum_uncert, expected.estimated_peak_continuum_uncert );
    errors.push_back( buffer );
  }
  
  

  if( !floats_equiv_enough(test.decision_threshold, expected.decision_threshold) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test decision_threshold (%.6G) does not equal expected (%.6G)",
             test.decision_threshold, expected.decision_threshold );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.detection_limit, expected.detection_limit) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test detection_limit (%.6G) does not equal expected (%.6G)",
             test.detection_limit, expected.detection_limit );
    errors.push_back( buffer );
  }
  

  if( !floats_equiv_enough(test.source_counts, expected.source_counts) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test source_counts (%.6G) does not equal expected (%.6G)",
             test.source_counts, expected.source_counts );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.lower_limit, expected.lower_limit) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test lower_limit (%.6G) does not equal expected (%.6G)",
             test.lower_limit, expected.lower_limit );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.upper_limit, expected.upper_limit) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test upper_limit (%.6G) does not equal expected (%.6G)",
             test.upper_limit, expected.upper_limit );
    errors.push_back( buffer );
  }
  
  if( errors.empty() )
    return;
  
  string err_msg = "CurrieMdaResult::equal_enough: test and expected values are not equal\n";
  for( size_t i = 0; i < errors.size(); ++i )
    err_msg += (i ? "\n\t" : "\t") + errors[i];
    
  throw runtime_error( err_msg );
}//CurrieMdaResult::equal_enough(...)
#endif //PERFORM_DEVELOPER_CHECKS



std::ostream &print_summary( std::ostream &strm, const CurrieMdaResult &result, const float w )
{
  /*
   {
   //If peak looks present, use 2.5 FWHM for peak width
   //If no peak is located, then use 1.2 FWHM for peak width
   const float fwhm = detector->peakResolutionFWHM(result.input.gamma_energy);
   const shared_ptr<const SpecUtils::EnergyCalibration> cal = result.input.spectrum->energy_calibration();
   const float peak_channel = cal->channel_for_energy(result.input.gamma_energy);
   const float lower_energy = result.input.roi_lower_energy;
   const float upper_energy = result.input.roi_upper_energy;
   
   strm << "Width is " << (upper_energy - lower_energy)/fwhm << " FWHMs" << endl;
   }
   */
  
  strm << "Activity less than ";
  if( w > 0.0 )
    strm << PhysicalUnits::printToBestActivityUnits(w*result.upper_limit) << " (" << result.upper_limit << " counts)" << endl;
  else
    strm << result.upper_limit << " counts" << endl;
  
  strm << "Activity greater than ";
  if( w > 0.0 )
    strm << PhysicalUnits::printToBestActivityUnits(w*result.lower_limit) << " (" << result.lower_limit << " counts)" << endl;
  else
    strm << result.lower_limit << " counts" << endl;
  
  strm << "Nominal activity estimate: ";
  if( w > 0.0 )
     strm << PhysicalUnits::printToBestActivityUnits(w*result.source_counts) <<" (" << result.source_counts << " counts)" << endl;
  else
    strm << result.source_counts << " counts" << endl;
  
  strm << "Continuum Starting Channel: " << result.first_lower_continuum_channel << endl;
  strm << "Continuum Ending Channel: " << result.last_upper_continuum_channel + 1 << endl;
  strm << "Peak Starting Channel: " << result.first_peak_region_channel << endl;
  strm << "Peak Ending Channel: " << result.last_peak_region_channel + 1 << endl;
  
  
  strm << "Critical Limit (decision threshold) in Peak Region: "<< result.decision_threshold << " counts";
  if( w > 0 )
    strm << " or " << PhysicalUnits::printToBestActivityUnits(w*result.decision_threshold);
  
  strm << "\nDetection limit: " << result.detection_limit << " counts";
  if( w > 0 )
    strm << ", or " << PhysicalUnits::printToBestActivityUnits(w*result.detection_limit) << endl;
  
  
  strm << "Gross Counts in Peak Foreground: " << result.source_counts << endl;
  //Gross Counts in Peak Background: 0
  strm << "Gross Counts in Continuum Foreground: " << result.estimated_peak_continuum_counts << endl;
  //Gross Counts in Continuum Background: 0
  strm << "Uncert in Peak Region: " << result.estimated_peak_continuum_uncert << endl;
  //Variance in Continuum Region: 14.17189
  //Detector Efficiency at Energy: 0.1304517
  //Range of Gammas Entering Detector: (0, 0.1656136)
  //Solid Angle: 0.0002390204
  //Range of Gammas from Source: (0, 692.8849)
  //Gamma Emission Rate per uCi of Source: 3.699356E+10
  //Transmission through Shielding: 1
  
  strm << "Continuum eqn (cnts/keV): " << result.continuum_eqn[0]
       << (result.continuum_eqn[1] > 0.0 ? " + " : " - ")
       << fabs(result.continuum_eqn[1]) << "*x"
       << ", with x=energy-" << result.input.gamma_energy << endl;
  
  return strm;
};//print_summary( CurrieMdaResult )

  
pair<size_t,size_t> round_roi_to_channels( shared_ptr<const SpecUtils::Measurement> spectrum,
                                  const float roi_lower_energy,
                                  const float roi_upper_energy )
{
  if( !spectrum )
    throw runtime_error( "mda_counts_calc: no spectrum" );
  
  shared_ptr<const SpecUtils::EnergyCalibration> cal = spectrum->energy_calibration();
  if( !cal || !cal->valid() )
    throw runtime_error( "mda_counts_calc: invalid energy calibration" );
  
  const double peak_region_lower_ch = std::max(0.0, cal->channel_for_energy( roi_lower_energy ) );
  const double peak_region_upper_ch = std::max(0.0, cal->channel_for_energy( roi_upper_energy ) );
  
  //if( (peak_region_lower_ch - num_lower_side_channels) < 0.0 )
  //  throw runtime_error( "mda_counts_calc: lower energy goes off spectrum" );
  
  //if( (peak_region_upper_ch + num_upper_side_channels) > cal->num_channels() )
  //  throw runtime_error( "mda_counts_calc: upper energy goes off spectrum" );
  
  size_t first_peak_region_channel = static_cast<size_t>( std::round(peak_region_lower_ch) );
  
  // If we pass in exactly the channel boundary, or really close to it, we want to round in the
  //  reasonable way, otherwise we need to makeup for the channel number defining the left side of
  //  each channel, so we will subtract off 0.5 from the channel we are supposed to go up through.
  size_t last_peak_region_channel;
  if( fabs(peak_region_upper_ch - std::floor(peak_region_upper_ch)) < 0.01 )
    last_peak_region_channel = static_cast<size_t>( std::floor(peak_region_upper_ch) - 1 );
  else
    last_peak_region_channel = static_cast<size_t>( std::round(peak_region_upper_ch - 0.5) );
  
  return make_pair( first_peak_region_channel, last_peak_region_channel );
}//round_roi_to_channels(...)
  
  
  
CurrieMdaResult currie_mda_calc( const CurrieMdaInput &input )
{
  using namespace SpecUtils;
  
  const shared_ptr<const Measurement> spec = input.spectrum;
  const shared_ptr<const EnergyCalibration> cal = spec ? spec->energy_calibration() : nullptr;
  const size_t nchannel = cal ? cal->num_channels() : size_t(0);
  
  if( !nchannel || !cal->valid() || !spec->gamma_counts() || !cal->channel_energies() )
    throw runtime_error( "mda_counts_calc: invalid spectrum passed in." );
  
  if( IsInf(input.roi_lower_energy) || IsNan(input.roi_lower_energy)
     ||  IsInf(input.roi_upper_energy) || IsNan(input.roi_upper_energy) )
    throw runtime_error( "mda_counts_calc: invalid ROI specified." );
  
  if( input.roi_lower_energy >= input.roi_upper_energy )
    throw runtime_error( "mda_counts_calc: upper ROI energy must be greater than lower energy." );
  
  if( (input.gamma_energy < input.roi_lower_energy)
      || (input.gamma_energy > input.roi_upper_energy) )
    throw runtime_error( "mda_counts_calc: gamma energy must be between lower and upper ROI." );
  
  if( ((input.num_lower_side_channels == 0) || (input.num_upper_side_channels == 0))
     && (input.num_lower_side_channels != input.num_upper_side_channels) )
    throw runtime_error( "mda_counts_calc: lower or upper side channels was zero, but not both." );
  
  if( input.num_lower_side_channels >= nchannel  )
    throw runtime_error( "mda_counts_calc: invalid num_lower_side_channels." );
  
  if( input.num_upper_side_channels >= nchannel )
    throw runtime_error( "mda_counts_calc: invalid num_upper_side_channels." );
  
  if( input.detection_probability <= 0.05 || input.detection_probability >= 1.0 )
    throw runtime_error( "mda_counts_calc: invalid detection_probability." );
  
  if( input.additional_uncertainty < 0.0f || input.additional_uncertainty >= 1.0f )
    throw runtime_error( "mda_counts_calc: invalid additional_uncertainty." );
  
  const vector<float> &gamma_counts = *spec->gamma_counts();
  const vector<float> &gamma_energies = *cal->channel_energies();
  
  assert( gamma_energies.size() == (gamma_counts.size() + 1) );
  
  const pair<size_t,size_t> channels = round_roi_to_channels( spec, input.roi_lower_energy, input.roi_upper_energy );
  
  CurrieMdaResult result;
  result.input = input;
  
  result.first_peak_region_channel = channels.first;
  result.last_peak_region_channel = channels.second;
  
  if( result.first_peak_region_channel < (input.num_lower_side_channels + 1) )
    throw std::runtime_error( "mda_counts_calc: lower peak region is outside spectrum energy range" );
  
  if( input.num_lower_side_channels == 0 )
  {
    result.first_lower_continuum_channel = 0;
    result.last_lower_continuum_channel  = 0;
    result.lower_continuum_counts_sum    = 0;
  }else
  {
    result.last_lower_continuum_channel = result.first_peak_region_channel - 1;
    result.first_lower_continuum_channel = result.last_lower_continuum_channel - input.num_lower_side_channels + 1;
    result.lower_continuum_counts_sum = spec->gamma_channels_sum(result.first_lower_continuum_channel, result.last_lower_continuum_channel);
  }
  
  if( input.num_upper_side_channels == 0 )
  {
    result.first_upper_continuum_channel = 0;
    result.last_upper_continuum_channel  = 0;
    result.upper_continuum_counts_sum    = 0;
  }else
  {
    result.first_upper_continuum_channel = result.last_peak_region_channel + 1;
    result.last_upper_continuum_channel = result.first_upper_continuum_channel + input.num_upper_side_channels - 1;
    
    if( result.last_upper_continuum_channel >= nchannel  )
      throw std::runtime_error( "mda_counts_calc: upper peak region is outside spectrum energy range" );
    
    result.upper_continuum_counts_sum = spec->gamma_channels_sum(result.first_upper_continuum_channel, result.last_upper_continuum_channel);
  }

  result.peak_region_counts_sum = spec->gamma_channels_sum(result.first_peak_region_channel, result.last_peak_region_channel);
  
  
  /*
   cout << "Lower region:\n\tChan\tEne\tCounts" << endl;
   for( size_t i = result.first_lower_continuum_channel; i <= result.last_lower_continuum_channel; ++i )
   cout << "\t" << i << "\t" << gamma_energies[i] << "\t" << gamma_counts[i] << endl;
   cout << "\tSum: " << result.lower_continuum_counts_sum << endl;
   
   cout << "\nPeak region:\n\tChan\tEne\tCounts" << endl;
   for( size_t i = result.first_peak_region_channel; i <= result.last_peak_region_channel; ++i )
   cout << "\t" << i << "\t" << gamma_energies[i] << "\t" << gamma_counts[i] << endl;
   cout << "\tSum: " << result.peak_region_counts_sum << endl;
   
   cout << "\nUpper region:\n\tChan\tEne\tCounts" << endl;
   for( size_t i = result.first_upper_continuum_channel; i <= result.last_upper_continuum_channel; ++i )
   cout << "\t" << i << "\t" << gamma_energies[i] << "\t" << gamma_counts[i] << endl;
   cout << "\tSum: " << result.upper_continuum_counts_sum << endl;
   */
  
  double peak_cont_sum_uncert = -999.9f, peak_cont_sum = -999.9f;
  
  if( input.num_upper_side_channels == 0 )
  {
    // Clamped for the same reason the fitted branch below is: a processed spectrum (background
    //  subtracted, or a corrected export) can hold negative channel counts, and `peak_cont_sigma`
    //  further down takes the square root of this.  Without it the limit comes out "nan".
    peak_cont_sum = (std::max)( 0.0, static_cast<double>(result.peak_region_counts_sum) );
    peak_cont_sum_uncert = sqrt( peak_cont_sum );
    
    const double peak_area_width = spec->gamma_channel_upper(result.last_peak_region_channel)
                                    - spec->gamma_channel_lower(result.first_peak_region_channel);
    result.continuum_eqn[1] = 0.0;
    result.continuum_eqn[0] = peak_cont_sum / peak_area_width;
  }else
  {
    const double lower_cont_counts = spec->gamma_channels_sum(result.first_lower_continuum_channel, result.last_lower_continuum_channel);
    const double upper_cont_counts = spec->gamma_channels_sum(result.first_upper_continuum_channel, result.last_upper_continuum_channel);
    const double lower_cont_width = spec->gamma_channel_upper(result.last_lower_continuum_channel)
                                    - spec->gamma_channel_lower(result.first_lower_continuum_channel);
    const double upper_cont_width = spec->gamma_channel_upper(result.last_upper_continuum_channel)
                                    - spec->gamma_channel_lower(result.first_upper_continuum_channel);
    
    // A width-weighted least-squares line through the side-band channels.
    //
    // The continuum is a density in counts/keV, so a channel of width `w` whose midpoint sits `x`
    //  from `input.gamma_energy` expects `w*(c0 + c1*x)` counts.  Weighting by 1/w - the Poisson
    //  variance of a slowly varying density - makes the normal equations pure moment matching,
    //  `sum(n) == sum(mu)` and `sum(n*x) == sum(mu*x)`.  Those are unbiased estimating equations for
    //  any binning and any count level, need no iteration, and never divide by an empty channel.
    //  Weighting by 1/n instead would bias the continuum low at low counts, which is the wrong
    //  direction for a detection limit; `decon_compute_peaks` documents the same trap for the
    //  sideband constraint of `DeconContinuumNorm::FixedByEdges`.
    //
    // This replaces averaging the two side-band *densities*, which equals the fitted line at the
    //  region centre only when the bands span equal energy.  The bands hold an equal number of
    //  *channels*, so any non-uniform binning biased the continuum by
    //  `slope*(upper_width - lower_width)/4*peak_width`: -256 counts, or -0.9 sigma of the
    //  continuum's own Poisson noise, on the 81 keV line of the bundled Ba-133 NaI spectrum, whose
    //  deviation pairs make the two bands 11.5 and 12.9 keV wide.  Projecting a dwell scales the
    //  bias by the projection factor while the noise only grows as its square root.

    // Walks the two side bands, handing each channel its width, its midpoint relative to
    //  `input.gamma_energy`, and its counts.
    const auto for_each_side_channel = [&]( auto &&fn ){
      const size_t bands[2][2] = {
        { result.first_lower_continuum_channel, result.last_lower_continuum_channel },
        { result.first_upper_continuum_channel, result.last_upper_continuum_channel }
      };

      for( const size_t (&band)[2] : bands )
      {
        for( size_t channel = band[0]; channel <= band[1]; ++channel )
        {
          const double lower = spec->gamma_channel_lower( channel );
          const double upper = spec->gamma_channel_upper( channel );
          fn( upper - lower, 0.5*(lower + upper) - input.gamma_energy,
             static_cast<double>( spec->gamma_channel_content(channel) ) );
        }
      }
    };//for_each_side_channel

    // Taken from the band totals rather than re-summed, so the fit and the reported per-band
    //  numbers cannot disagree.
    const double sum_w = lower_cont_width + upper_cont_width;
    const double sum_n = lower_cont_counts + upper_cont_counts;

    double sum_wx = 0.0;
    for_each_side_channel( [&sum_wx]( const double w, const double x, const double ){
      sum_wx += w*x;
    } );

    // The width-weighted centroid of the two bands.  Solving about it decouples the two
    //  coefficients, so there is no 2x2 system and no conditioning worry.  For uniform channels and
    //  an equal number of side channels each side it is exactly the peak region's centre, which is
    //  what collapses everything below to the mean-of-densities value this replaces.
    const double x_bar = (sum_w > 0.0) ? (sum_wx / sum_w) : 0.0;

    double sum_wuu = 0.0, sum_wuuu = 0.0, sum_nu = 0.0;
    for_each_side_channel( [&]( const double w, const double x, const double n ){
      const double u = x - x_bar;
      sum_wuu  += w*u*u;
      sum_wuuu += w*u*u*u;
      sum_nu   += n*u;
    } );

    // `sum_wuu` can only vanish if every side channel shares one midpoint, which cannot happen;
    //  the guard is for a corrupted calibration rather than a reachable case.
    const double cont_density_at_centroid = (sum_w > 0.0) ? (sum_n / sum_w) : 0.0;
    const double cont_slope = (sum_wuu > 0.0) ? (sum_nu / sum_wuu) : 0.0;

    const double peak_area_lower = spec->gamma_channel_lower(result.first_peak_region_channel);
    const double peak_area_upper = spec->gamma_channel_upper(result.last_peak_region_channel);
    const double peak_area_width = peak_area_upper - peak_area_lower;

    // How far the peak region's centre sits from the band centroid; zero for uniform binning.
    const double delta = 0.5*(peak_area_lower + peak_area_upper) - input.gamma_energy - x_bar;

    // Integrating a line over the peak region is its value at the region's centre times the width.
    const double peak_cont_integral = peak_area_width*(cont_density_at_centroid + cont_slope*delta);

    // Clamped because a fitted line can go negative where a mean of two densities could not, and
    //  `peak_cont_sigma` below takes the square root of this.
    peak_cont_sum = (std::max)( 0.0, peak_cont_integral );

    // The sandwich covariance `M^-1 (B^T V B) M^-1`, with V the per-channel Poisson variances taken
    //  at the fitted expectation.  In the centred basis it is three scalars.  Note the normal
    //  equations make `sum(mu) == sum(n)` and `sum(mu*u) == sum(n*u)` hold exactly, so a fitted and
    //  an observed variance estimate can only differ in the u^2 term.
    const double sum_mu_uu = cont_density_at_centroid*sum_wuu + cont_slope*sum_wuuu;
    const double var_density = (sum_w > 0.0) ? (sum_n / (sum_w*sum_w)) : 0.0;
    const double cov_density_slope = ((sum_w > 0.0) && (sum_wuu > 0.0))
                                       ? (sum_nu / (sum_w*sum_wuu)) : 0.0;
    const double var_slope = (sum_wuu > 0.0) ? (sum_mu_uu / (sum_wuu*sum_wuu)) : 0.0;
    const double var_peak_density = var_density + 2.0*delta*cov_density_slope
                                     + delta*delta*var_slope;

    peak_cont_sum_uncert = peak_area_width*sqrt( (std::max)( 0.0, var_peak_density ) );

    // `continuum_eqn` is written about `input.gamma_energy`, so shift the fit off the centroid.
    result.continuum_eqn[1] = cont_slope;
    result.continuum_eqn[0] = cont_density_at_centroid - cont_slope*x_bar;

#if( PERFORM_DEVELOPER_CHECKS )
    {// begin sanity check on continuum eqn
      // Now an algebraic identity rather than an approximation: the reported continuum counts *are*
      //  the reported continuum equation integrated over the peak region.  The previous check
      //  needed a 1E-3 fudge, and blamed "really bad numerical accuracy ... for NaI systems" - that
      //  was the unequal-band-width bias described above, not numerics, and it aborted Debug builds
      //  on the one non-linearly-calibrated spectrum in the test suite.
      const double peak_start_eq = peak_area_lower - input.gamma_energy;
      const double peak_end_eq = peak_area_upper - input.gamma_energy;

      const double peak_cont_eq_integral = result.continuum_eqn[0] * (peak_end_eq - peak_start_eq)
      + result.continuum_eqn[1] * 0.5 * (peak_end_eq*peak_end_eq - peak_start_eq*peak_start_eq);

      const double eq_diff = fabs( peak_cont_eq_integral - peak_cont_integral );
      assert( eq_diff <= 1.0E-8*(std::max)( 1.0, fabs(peak_cont_integral) ) );
    }// end sanity check on continuum eqn
#endif //PERFORM_DEVELOPER_CHECKS
  }//if( input.num_upper_side_channels == 0 ) / else
  
  assert( peak_cont_sum_uncert != -999.9f );
  assert( peak_cont_sum != -999.9f );
  result.estimated_peak_continuum_counts = static_cast<float>( peak_cont_sum );
  result.estimated_peak_continuum_uncert = static_cast<float>( peak_cont_sum_uncert );
  
  typedef boost::math::policies::policy<boost::math::policies::digits10<6> > my_pol_6;
  const boost::math::normal_distribution<double,my_pol_6> gaus_dist( 0.0, 1.0 );
  
  //  TODO: If/when we start having k_alpha != k_beta, then we probably need to be more careful
  //        around single vs double sided quantile.
  //   Will map 0.8414->1.00023, 0.95->1.64485, 0.975->1.95996, 0.995->2.57583, ...
  const double k = boost::math::quantile( gaus_dist, input.detection_probability );
  
  
  const double peak_cont_sigma = sqrt( peak_cont_sum_uncert*peak_cont_sum_uncert + peak_cont_sum );
  
  result.decision_threshold = k * peak_cont_sigma; //Note if using non-symmetric coverage, we would use k_alpha here
  
  // TODO: The calculation of result.detection_limit is using the simplified form requiring k_alpha == k_beta.
  //       If this is not the case it is an iterative solution (see eqn 129 on pg 47 of AQ-48)
  const double add_uncert = input.additional_uncertainty;
  if( k*k*add_uncert*add_uncert >= 1.0 )
  {
    result.detection_limit = -999; //TODO: indicate non-applicability better
  }else
  {
    // TODO: Using tbl 16 of AQ-48, I get a slightly high answer of 193.05 vs the expected answer of 191.906.
    //       If in the numerator I replace k*k with just k, then I get 191.983
    //       And if I plug those tables numbers into eqn 129 I get 0.474598 vs their 0.471705,
    //       so I am currently suspecting an error in the table.
    result.detection_limit = (2.0*result.decision_threshold + k*k) / (1.0 - k*k*add_uncert*add_uncert);
  }
  
  
  const float source_counts = result.peak_region_counts_sum - result.estimated_peak_continuum_counts;
  result.source_counts = source_counts;
  
  double region_sigma = peak_cont_sum_uncert*peak_cont_sum_uncert + result.peak_region_counts_sum;
  
  // TODO: I *think* this is right; e.g., use the nominal estimate of signal counts to estimate total uncertainty impact due to the "additional uncertainty" of the measurement, but I need to double check this.
  if( (source_counts > 0) && (add_uncert > 0) )
    region_sigma += source_counts*source_counts * add_uncert*add_uncert;
  region_sigma = sqrt( region_sigma );
  
  result.lower_limit = source_counts - k*region_sigma;
  result.upper_limit = source_counts + k*region_sigma;
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  assert( !IsInf(result.lower_limit) && !IsNan(result.lower_limit) );
  assert( !IsInf(result.upper_limit) && !IsNan(result.upper_limit) );
#endif
  
  /*
   cout << "lower_cont_counts=" << lower_cont_counts << endl;
   cout << "upper_cont_counts=" << upper_cont_counts << endl;
   cout << "peak_cont_density=" << peak_cont_density << endl;
   cout << "peak_cont_density_uncert=" << peak_cont_density_uncert << endl;
   cout << "peak_cont_sum=" << peak_cont_sum << endl;
   cout << "peak_cont_sum_uncert=" << peak_cont_sum_uncert << endl;
   cout << "peak_cont_sigma=" << peak_cont_sigma << endl;
   cout << "result.decision_threshold=" << result.decision_threshold << endl;
   cout << "result.detection_limit=" << result.detection_limit << endl;
   cout << "result.source_counts=" << result.source_counts << endl;
   cout << "result.lower_limit=" << result.lower_limit << endl;
   cout << "result.upper_limit=" << result.upper_limit << endl;
   */
  
  return result;
};//mda_counts_calc


const char *to_str( const PeakCurrieCheck::ResultType type )
{
  switch( type )
  {
    case PeakCurrieCheck::ResultType::NotDetected: return "NotDetected";
    case PeakCurrieCheck::ResultType::Detected:    return "Detected";
    case PeakCurrieCheck::ResultType::Deficit:     return "Deficit";
    case PeakCurrieCheck::ResultType::Error:       return "Error";
  }//switch( type )

  assert( 0 );
  return "Error";
}//const char *to_str( const PeakCurrieCheck::ResultType )


double peak_width_for_currie_check( const PeakDef &peak )
{
  if( peak.gausPeak() && (peak.fwhm() > 0.0) )
    return peak.fwhm();

  const double roi_width = peak.upperX() - peak.lowerX();

  return (roi_width > 0.0) ? (roi_width / 2.5) : 0.0;
}//double peak_width_for_currie_check( const PeakDef & )


double gaussian_fraction_in_roi( const double num_fwhm )
{
  if( num_fwhm <= 0.0 )
    return 0.0;

  // Half-width, in units of sigma, is 0.5*num_fwhm*2.35482; the fraction within +-x sigma of the
  //  mean is erf( x / sqrt(2) ).
  const double half_width_sigma = 0.5 * num_fwhm * 2.35482;

  return erf( half_width_sigma / sqrt(2.0) );
}//double gaussian_fraction_in_roi( const double )


std::string confidence_level_str( const double confidence_level )
{
  char buffer[64] = { '\0' };

  if( confidence_level < 0.999 )
    snprintf( buffer, sizeof(buffer), "%.4g%%", 100.0*confidence_level );
  else
    snprintf( buffer, sizeof(buffer), "1-%.2G", (1.0 - confidence_level) );

  return buffer;
}//std::string confidence_level_str( const double )


std::string confidence_level_pct_str( const double confidence_level )
{
  if( IsNan(confidence_level) || !(confidence_level > 0.0) || !(confidence_level < 1.0) )
    return "?";

  char buffer[64] = { '\0' };

  if( confidence_level < 0.999 )
  {
    snprintf( buffer, sizeof(buffer), "%.4g%%", 100.0*confidence_level );
  }else
  {
    // Show enough decimals for the complement to carry two significant figures.  The 4-sigma and
    //  5-sigma selections differ from 100% only in the 5th and 7th decimal place, so any fixed
    //  precision prints both as "100%" and makes them indistinguishable to the user.  Gives
    //  99.73% / 99.9937% / 99.999943% for 3/4/5 sigma.
    const double complement_pct = 100.0*(1.0 - confidence_level);
    const int decimals = (std::min)( 12,
                     (std::max)( 2, 1 - static_cast<int>( std::floor( std::log10(complement_pct) ) ) ) );
    snprintf( buffer, sizeof(buffer), "%.*f%%", decimals, 100.0*confidence_level );
  }

  return buffer;
}//std::string confidence_level_pct_str( const double )


std::pair<double,double> currie_check_energy_range( const CurrieMdaResult &result )
{
  const std::shared_ptr<const SpecUtils::Measurement> &spec = result.input.spectrum;
  assert( spec && spec->energy_calibration() );
  if( !spec || !spec->energy_calibration() )
    return { 0.0, 0.0 };

  const size_t first_channel = (result.input.num_lower_side_channels > 0)
                                  ? result.first_lower_continuum_channel
                                  : result.first_peak_region_channel;
  const size_t last_channel = (result.input.num_upper_side_channels > 0)
                                  ? result.last_upper_continuum_channel
                                  : result.last_peak_region_channel;

  return { spec->gamma_channel_lower(first_channel), spec->gamma_channel_upper(last_channel) };
}//std::pair<double,double> currie_check_energy_range( const CurrieMdaResult & )


PeakCurrieCheck currie_check_for_peak( const PeakDef &peak,
                                      const shared_ptr<const SpecUtils::Measurement> &spectrum,
                                      const PeakCurrieCheckOptions &options,
                                      const bool peak_was_fit )
{
  PeakCurrieCheck answer;

  // A side-channel count of zero puts `currie_mda_calc` into its "the spectrum is asserted to be
  //  background" mode, which is a different calculation than we are documenting here.
  const size_t num_side_channels = std::max( size_t(1), options.num_side_channels );

  // Record the options used up front, so they are still reported for peaks whose limit cant be
  //  computed.
  answer.result.input.detection_probability = options.confidence_level;
  answer.result.input.num_lower_side_channels = num_side_channels;
  answer.result.input.num_upper_side_channels = num_side_channels;
  answer.signal_fraction_in_roi = gaussian_fraction_in_roi( options.roi_num_fwhm );

  if( !spectrum || (spectrum->num_gamma_channels() < 4) )
  {
    answer.error_message = "no spectrum to evaluate the peak region in";
    answer.short_description = "Not computed";
    answer.result_summary = "Detection limit could not be computed: " + answer.error_message + ".";
    return answer;
  }//if( no usable spectrum )

  const double mean = peak.mean();
  const double fwhm = peak_width_for_currie_check( peak );

  CurrieMdaInput input;
  input.spectrum = spectrum;
  input.gamma_energy = static_cast<float>( mean );

  if( fwhm > 0.0 )
  {
    const double half_width = 0.5 * options.roi_num_fwhm * fwhm;
    input.roi_lower_energy = static_cast<float>( mean - half_width );
    input.roi_upper_energy = static_cast<float>( mean + half_width );
  }else
  {
    input.roi_lower_energy = static_cast<float>( peak.lowerX() );
    input.roi_upper_energy = static_cast<float>( peak.upperX() );
  }

  input.num_lower_side_channels = num_side_channels;
  input.num_upper_side_channels = num_side_channels;
  input.detection_probability = options.confidence_level;
  input.additional_uncertainty = 0.0f;

  try
  {
    answer.result = currie_mda_calc( input );
    answer.computed = true;
  }catch( std::exception &e )
  {
    // Keep the input around, so the options actually used are still reported for this peak
    answer.result.input = input;
    answer.error_message = e.what();
    answer.result_type = PeakCurrieCheck::ResultType::Error;
    answer.short_description = "Not computed";
    answer.result_summary = "Detection limit could not be computed: " + answer.error_message + ".";
    return answer;
  }//try / catch

  const CurrieMdaResult &res = answer.result;
  const string cl_str = confidence_level_str( options.confidence_level );

  // With (essentially) no counts anywhere in the region, the decision threshold and the upper
  //  limit both collapse to zero - the Gaussian statistics the calculation uses have broken
  //  down.  Left alone, that makes a single stray count read as a detection, and reports an
  //  upper limit of "less than 0 counts".  Ld does not degenerate (it tends to k^2, which is
  //  the right order for a zero-background Poisson limit), so fall back to quoting just that.
  answer.region_is_empty = (res.decision_threshold <= 0.0f);

  // Classify the result, using the same criteria as the GUI detection limit tools.
  if( answer.region_is_empty )
    answer.result_type = PeakCurrieCheck::ResultType::NotDetected;
  else if( res.source_counts > res.decision_threshold )
    answer.result_type = PeakCurrieCheck::ResultType::Detected;
  else if( res.upper_limit < 0.0f )
    answer.result_type = PeakCurrieCheck::ResultType::Deficit;
  else
    answer.result_type = PeakCurrieCheck::ResultType::NotDetected;

  switch( answer.result_type )
  {
    case PeakCurrieCheck::ResultType::NotDetected:
      if( answer.region_is_empty )
      {
        answer.short_description = "No counts in region";
        answer.result_summary = "Not detected; the peak region contains essentially no counts."
                                "  Minimum reliably detectable signal (Ld) is "
                                + SpecUtils::printCompact(res.detection_limit, 4) + " counts.";
        break;
      }//if( answer.region_is_empty )

      // A peak that was fit, yet sits below the decision threshold, passed the peak-fit
      //  significance tests but cannot be claimed as a detection - worth saying plainly.
      answer.short_description = peak_was_fit ? "Fit, but less than Lc" : "Less than Lc";
      answer.result_summary = (peak_was_fit
                    ? string("A peak was fit here, but the signal is below the decision threshold")
                    : string("Not detected.  Observed signal is below the decision threshold"))
                    + " (Lc = " + SpecUtils::printCompact(res.decision_threshold, 4) + " counts);"
                    " signal is less than " + SpecUtils::printCompact(res.upper_limit, 4)
                    + " counts at the " + cl_str + " confidence level."
                    "  Minimum reliably detectable signal (Ld) is "
                    + SpecUtils::printCompact(res.detection_limit, 4) + " counts.";
      break;

    case PeakCurrieCheck::ResultType::Detected:
      answer.short_description = "Greater than Lc";
      answer.result_summary = (peak_was_fit
                    ? string("Signal present: observed ")
                    : string("Signal present but peak not fit: observed "))
                    + SpecUtils::printCompact(res.source_counts, 4) + " counts is above the"
                    " decision threshold (Lc = "
                    + SpecUtils::printCompact(res.decision_threshold, 4) + " counts); signal"
                    " is between " + SpecUtils::printCompact(res.lower_limit, 4) + " and "
                    + SpecUtils::printCompact(res.upper_limit, 4) + " counts at the "
                    + cl_str + " confidence level.";
      break;

    case PeakCurrieCheck::ResultType::Deficit:
      answer.short_description = "Fewer counts than expected";
      answer.result_summary = "Fewer counts than expected from the neighboring continuum; observed"
                    " signal is consistent with less than 0 counts at the " + cl_str
                    + " confidence level.";
      break;

    case PeakCurrieCheck::ResultType::Error:
      assert( 0 );
      break;
  }//switch( answer.result_type )

  return answer;
}//PeakCurrieCheck currie_check_for_peak(...)


namespace
{
/** Smallest expected channel counts we will allow, as a fraction of the mean observed counts.

 A *relative* floor, unlike the modified-Neyman variance floor of 1.0 count that this calculation
 used to rely on: it keeps `log(observed/expected)` finite without imposing an absolute scale on
 spectra that have been scaled, background subtracted, or taken for very short dwells.
 */
const double sm_min_expected_fraction = 1.0E-8;

/** Poisson deviance of one channel, without input validation - for use in the optimizer's inner
 loops.  `expected` must already be known positive.
 */
inline double poisson_deviance_fast( const double observed, const double expected )
{
  // The `observed*log(observed/expected)` term is defined to be zero at observed == 0, which is
  //  the limit of x*log(x) as x -> 0.
  if( observed <= 0.0 )
    return 2.0*expected;

  return 2.0*( expected - observed + observed*std::log(observed/expected) );
}//poisson_deviance_fast(...)


/** The channel-integrated polynomial continuum basis.

 Element (i,k) is the integral of `(E - reference_energy)^k` over channel `i`, i.e.

 \code
   I(i,k) = ( (E1 - ref)^(k+1) - (E0 - ref)^(k+1) ) / (k+1)
 \endcode

 This is deliberately the same basis `PeakContinuum::offset_eqn_integral(...)` and
 `PeakFit::fit_amp_and_offset_imp(...)` use, so a continuum fit here and a continuum evaluated
 elsewhere cannot drift apart.  Column zero is the channel width, hence strictly positive, which
 is what lets us restore positivity by shifting the constant term.

 Each row is multiplied by its channel's `continuum_scale`, so that the expected counts stay
 `basis*coefficients + fixed_signal` for a projected measurement exactly as they do for the
 spectrum the continuum is defined on.  Every derivative in the optimizer is taken through
 `m_basis`, so nothing downstream needs to know about the scaling.
 */
Eigen::MatrixXd continuum_basis( const PoissonChannel * const channels,
                                const size_t nbin,
                                const size_t num_coefficients,
                                const double reference_energy )
{
  Eigen::MatrixXd basis( static_cast<Eigen::Index>(nbin),
                        static_cast<Eigen::Index>(num_coefficients) );

  for( size_t i = 0; i < nbin; ++i )
  {
    const double x0 = channels[i].lower_energy - reference_energy;
    const double x1 = channels[i].upper_energy - reference_energy;
    const double scale = channels[i].continuum_scale;

    double x0_power = x0, x1_power = x1;
    for( size_t k = 0; k < num_coefficients; ++k )
    {
      basis( static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(k) )
                              = scale * (x1_power - x0_power) / static_cast<double>(k + 1);
      x0_power *= x0;
      x1_power *= x1;
    }
  }//for( size_t i = 0; i < nbin; ++i )

  return basis;
}//continuum_basis(...)


/** The Poisson continuum objective, in an equilibrated (column-scaled) coefficient basis.

 Equilibration matters more than it looks: the raw basis columns span many orders of magnitude
 (column 0 is a channel width of a few keV, column 1 integrates an energy offset that can be
 hundreds of keV, column 2 its square), so the raw Newton system is badly conditioned purely
 because of the reference energy and ROI width.  Scaling each column to unit RMS removes that,
 and the scaling is undone before the coefficients are handed back.
 */
struct PoissonContinuumObjective
{
  Eigen::MatrixXd m_basis;          //!< equilibrated; (nbin, num_coefficients)
  Eigen::VectorXd m_scale;          //!< raw_coefficient(k) = scaled_coefficient(k) * m_scale(k)
  Eigen::VectorXd m_observed;
  Eigen::VectorXd m_fixed_signal;
  Eigen::VectorXd m_constraint_center;    //!< in scaled space; empty if unconstrained
  Eigen::MatrixXd m_constraint_precision; //!< in scaled space; empty if unconstrained
  double m_min_expected = 0.0;

  /** Log-barrier weight on the constraint `expected_i > m_min_expected`.

   The Poisson deviance is its own barrier for any channel that saw counts - `-2*n*log(E)` blows up
   as `E` falls - but a channel with zero counts contributes `2*E`, which simply decreases, so the
   minimizer pushes those channels down until the positivity constraint stops them.  That makes the
   constraint genuinely active for exactly the low-count, partly-empty regions this whole
   calculation is about, and a plain Newton method cannot handle an active constraint: once a step
   is blocked, there can still be a descent direction *along* the active face.

   So the constrained problem is solved as a short sequence of smooth barrier subproblems with
   decreasing #m_barrier_mu.  Every iterate stays strictly interior, the Newton system stays well
   posed, and the solution approaches the true constrained minimum from inside.  The duality gap is
   bounded by (number of channels)*mu, which is how the stopping point is chosen.
   */
  double m_barrier_mu = 0.0;

  size_t num_coefficients() const { return static_cast<size_t>( m_basis.cols() ); }

  Eigen::VectorXd expected( const Eigen::VectorXd &coefs ) const
  {
    return (m_basis * coefs) + m_fixed_signal;
  }

  /** The smallest expected channel counts; the fit is only valid while this exceeds
   #m_min_expected.
   */
  double min_expected( const Eigen::VectorXd &coefs ) const
  {
    return expected(coefs).minCoeff();
  }

  /** The Poisson deviance plus any Gaussian constraint penalty - the quantity actually reported.

   Excludes the barrier, which is only a device for solving the constrained problem.
   */
  double deviance( const Eigen::VectorXd &coefs ) const
  {
    const Eigen::VectorXd expect = expected( coefs );
    if( !expect.allFinite() || (expect.minCoeff() < m_min_expected) )
      return std::numeric_limits<double>::infinity();

    double answer = 0.0;
    for( Eigen::Index i = 0; i < expect.size(); ++i )
      answer += poisson_deviance_fast( m_observed(i), expect(i) );

    if( m_constraint_center.size() )
    {
      const Eigen::VectorXd offset = coefs - m_constraint_center;
      answer += offset.dot( m_constraint_precision * offset );
    }

    return std::isfinite(answer) ? answer : std::numeric_limits<double>::infinity();
  }//deviance(...)

  /** The quantity the optimizer minimizes: #deviance plus the log barrier.

   Infinite outside the strict interior, which is what keeps the line search feasible.
   */
  double value( const Eigen::VectorXd &coefs ) const
  {
    double answer = deviance( coefs );
    if( !std::isfinite(answer) || (m_barrier_mu <= 0.0) )
      return answer;

    const Eigen::VectorXd expect = expected( coefs );
    for( Eigen::Index i = 0; i < expect.size(); ++i )
    {
      const double slack = expect(i) - m_min_expected;
      if( !(slack > 0.0) )
        return std::numeric_limits<double>::infinity();
      answer -= m_barrier_mu * std::log( slack );
    }

    return std::isfinite(answer) ? answer : std::numeric_limits<double>::infinity();
  }//value(...)

  /** Gradient and Hessian of #value at \p coefs; \p coefs must give positive expected counts. */
  void derivatives( const Eigen::VectorXd &coefs,
                   Eigen::VectorXd &gradient,
                   Eigen::MatrixXd &hessian ) const
  {
    const Eigen::VectorXd expect = expected( coefs );

    // d(deviance)/dE = 2*(1 - n/E) and d2(deviance)/dE2 = 2*n/E^2, chained through the linear
    //  basis; the Hessian is positive semi-definite because every observed count is non-negative.
    Eigen::VectorXd d_expect( expect.size() );
    Eigen::VectorXd dd_expect( expect.size() );
    for( Eigen::Index i = 0; i < expect.size(); ++i )
    {
      const double ratio = m_observed(i) / expect(i);
      d_expect(i) = 2.0*(1.0 - ratio);
      dd_expect(i) = 2.0*ratio/expect(i);
    }

    gradient = m_basis.transpose() * d_expect;
    hessian = m_basis.transpose() * dd_expect.asDiagonal() * m_basis;

    if( m_constraint_center.size() )
    {
      const Eigen::VectorXd offset = coefs - m_constraint_center;
      gradient += 2.0 * (m_constraint_precision * offset);
      hessian += 2.0 * m_constraint_precision;
    }

    if( m_barrier_mu > 0.0 )
    {
      // d/dc_k [ -mu*log(E_i - f) ] = -mu*B_ik/(E_i - f), and the Hessian term
      //  +mu*B_ik*B_il/(E_i - f)^2 is positive semi-definite, so the subproblem stays convex.
      Eigen::VectorXd inverse_slack( expect.size() );
      Eigen::VectorXd inverse_slack_squared( expect.size() );
      for( Eigen::Index i = 0; i < expect.size(); ++i )
      {
        const double slack = (std::max)( expect(i) - m_min_expected,
                                        std::numeric_limits<double>::min() );
        inverse_slack(i) = m_barrier_mu / slack;
        inverse_slack_squared(i) = m_barrier_mu / (slack*slack);
      }

      gradient -= m_basis.transpose() * inverse_slack;
      hessian += m_basis.transpose() * inverse_slack_squared.asDiagonal() * m_basis;
    }

    assert( gradient.size() == m_basis.cols() );
    assert( (hessian.rows() == m_basis.cols()) && (hessian.cols() == m_basis.cols()) );
  }//derivatives(...)

  /** Expected (Fisher) information of the *scaled* coefficients, including the constraint.

   Uses `sum_i B(i,k)*B(i,l)/E_i` rather than the observed Hessian, so it stays positive definite
   even when most channels are empty.
   */
  Eigen::MatrixXd information( const Eigen::VectorXd &coefs ) const
  {
    const Eigen::VectorXd expect = expected( coefs );
    Eigen::VectorXd inverse_expected( expect.size() );
    for( Eigen::Index i = 0; i < expect.size(); ++i )
      inverse_expected(i) = 1.0 / expect(i);

    Eigen::MatrixXd answer = m_basis.transpose() * inverse_expected.asDiagonal() * m_basis;
    if( m_constraint_center.size() )
      answer += m_constraint_precision;

    return answer;
  }//information(...)

  /** Shifts the constant term up until every expected channel count clears \p target.

   Column zero of the basis is the channel width, so for an increasing energy calibration it is
   strictly positive and this always succeeds; the guard below covers the case where it is not.

   \p target should leave some slack above #m_min_expected: the barrier is undefined exactly on the
   constraint, and a point placed precisely on it has nowhere to step.
   */
  bool restore_positivity( Eigen::VectorXd &coefs, const double target ) const
  {
    if( !coefs.allFinite() )
      return false;

    const Eigen::VectorXd expect = expected( coefs );
    if( !expect.allFinite() )
      return false;

    double shift = 0.0;
    for( Eigen::Index i = 0; i < expect.size(); ++i )
    {
      const double basis_value = m_basis(i,0);
      if( basis_value <= 0.0 )
        continue;
      shift = (std::max)( shift, (target - expect(i)) / basis_value );
    }

    if( shift > 0.0 )
    {
      // Rounding in the division above can leave a channel a fraction of an ulp short of the
      //  target, which would then read as a failure; a hair of extra shift removes that.
      coefs(0) += shift*(1.0 + 1.0E-12) + std::numeric_limits<double>::min();
    }

    return expected(coefs).minCoeff() > m_min_expected;
  }//restore_positivity(...)

  /** An interior starting point: comfortably above the constraint, but small compared with the
   data, so the barrier path does not start somewhere silly.
   */
  double interior_target() const
  {
    const double mean_observed = m_observed.sum() / static_cast<double>( m_observed.size() );
    return (std::max)( 4.0*m_min_expected, 1.0E-3*(std::max)( 1.0, mean_observed ) );
  }
};//struct PoissonContinuumObjective


/** Damped-Newton minimization with adaptive Levenberg regularization and a backtracking line
 search that also enforces positivity.

 Returns true on convergence, and adds the iterations it used to \p num_iterations.
 */
bool minimize_poisson_newton( const PoissonContinuumObjective &objective,
                             Eigen::VectorXd &coefs,
                             double &objective_value,
                             size_t &num_iterations )
{
  const size_t max_iterations = 50;
  const size_t max_halvings = 40;
  const Eigen::Index npar = static_cast<Eigen::Index>( objective.num_coefficients() );

  objective_value = objective.value( coefs );
  if( !std::isfinite(objective_value) )
    return false;

  double lambda = 0.0;

  for( size_t iteration = 0; iteration < max_iterations; ++iteration )
  {
    ++num_iterations;

    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;
    objective.derivatives( coefs, gradient, hessian );

    if( !gradient.allFinite() || !hessian.allFinite() )
      return false;

    // Damp against a single curvature scale rather than each diagonal element: because the basis
    //  is equilibrated the coefficients are already comparable, and this keeps the step bounded
    //  even when the Hessian is entirely zero (which happens when every channel is empty).
    double hessian_scale = 0.0;
    for( Eigen::Index k = 0; k < npar; ++k )
      hessian_scale = (std::max)( hessian_scale, std::fabs(hessian(k,k)) );
    if( !(hessian_scale > 0.0) )
      hessian_scale = 1.0;

    Eigen::MatrixXd damped = hessian;
    for( Eigen::Index k = 0; k < npar; ++k )
      damped(k,k) += (lambda + 1.0E-12)*hessian_scale;

    const Eigen::LDLT<Eigen::MatrixXd> decomposition( damped );
    const Eigen::VectorXd step = decomposition.solve( -gradient );

    if( (decomposition.info() != Eigen::Success) || !step.allFinite() )
    {
      lambda = (lambda > 0.0) ? 10.0*lambda : 1.0E-4;
      if( lambda > 1.0E12 )
        return false;
      continue;
    }

    // Newton decrement; -g.step is (to second order) twice the improvement still available.
    const double decrement = -gradient.dot( step );
    if( (decrement >= 0.0) && (decrement < 1.0E-10*(std::max)(1.0, std::fabs(objective_value))) )
      return true;

    bool accepted = false;
    double scale = 1.0;
    for( size_t halving = 0; halving < max_halvings; ++halving )
    {
      const Eigen::VectorXd candidate = coefs + scale*step;
      const double candidate_value = objective.value( candidate );
      if( candidate_value < objective_value )
      {
        const double improvement = objective_value - candidate_value;
        coefs = candidate;
        objective_value = candidate_value;
        accepted = true;

        if( improvement < 1.0E-9*(std::max)(1.0, std::fabs(objective_value)) )
          return true;

        break;
      }//if( the candidate is an improvement )

      scale *= 0.5;
    }//for( size_t halving = 0; halving < max_halvings; ++halving )

    if( accepted )
    {
      lambda = (lambda <= 1.0E-4) ? 0.0 : 0.5*lambda;
    }else
    {
      // Deliberately NOT treated as convergence.  A blocked step means this *direction* leaves the
      //  feasible region, which says nothing about whether a descent direction remains along the
      //  active face - and with a log barrier in play it usually just means the step was too long.
      //  Shorten it by damping and try again; if damping runs out, report failure so a different
      //  starting estimate or the derivative-free fallback gets a turn.
      lambda = (lambda > 0.0) ? 10.0*lambda : 1.0E-4;
      if( lambda > 1.0E12 )
        return false;
    }
  }//for( size_t iteration = 0; iteration < max_iterations; ++iteration )

  // Ran out of iterations, but we may still be at a perfectly good minimum; let the caller decide
  //  by reporting non-convergence so a fallback is attempted.
  return false;
}//minimize_poisson_newton(...)


/** Bounded Nelder-Mead on the same objective, used only when every Newton attempt has failed.

 Slow and derivative free, but it does not care about curvature pathologies, so it rescues the
 handful of very-low-count configurations where the Newton system is hopeless.
 */
bool minimize_poisson_simplex( const PoissonContinuumObjective &objective,
                              Eigen::VectorXd &coefs,
                              double &objective_value,
                              size_t &num_iterations )
{
  const Eigen::Index npar = static_cast<Eigen::Index>( objective.num_coefficients() );
  const size_t max_iterations = 400;

  std::vector<Eigen::VectorXd> vertices;
  std::vector<double> values;

  vertices.push_back( coefs );
  for( Eigen::Index k = 0; k < npar; ++k )
  {
    Eigen::VectorXd vertex = coefs;
    const double magnitude = (std::max)( std::fabs(vertex(k)), 1.0E-3 );
    vertex(k) += 0.1*magnitude;
    if( !objective.restore_positivity( vertex, objective.interior_target() ) )
      return false;
    vertices.push_back( vertex );
  }

  for( const Eigen::VectorXd &vertex : vertices )
  {
    const double value = objective.value( vertex );
    if( !std::isfinite(value) )
      return false;
    values.push_back( value );
  }

  const auto reorder = [&vertices,&values]() {
    std::vector<size_t> order( values.size() );
    std::iota( begin(order), end(order), size_t(0) );
    std::sort( begin(order), end(order),
              [&values]( const size_t a, const size_t b ){ return values[a] < values[b]; } );

    std::vector<Eigen::VectorXd> sorted_vertices;
    std::vector<double> sorted_values;
    for( const size_t index : order )
    {
      sorted_vertices.push_back( vertices[index] );
      sorted_values.push_back( values[index] );
    }
    vertices.swap( sorted_vertices );
    values.swap( sorted_values );
  };

  for( size_t iteration = 0; iteration < max_iterations; ++iteration )
  {
    ++num_iterations;
    reorder();

    const double best = values.front();
    const double worst = values.back();
    if( (worst - best) < 1.0E-10*(std::max)(1.0, std::fabs(best)) )
      break;

    Eigen::VectorXd centroid = Eigen::VectorXd::Zero( npar );
    for( size_t i = 0; i + 1 < vertices.size(); ++i )
      centroid += vertices[i];
    centroid /= static_cast<double>( vertices.size() - 1 );

    const auto try_point = [&objective,&centroid,&vertices]( const double coefficient ) {
      return Eigen::VectorXd( centroid + coefficient*(centroid - vertices.back()) );
    };

    const Eigen::VectorXd reflected = try_point( 1.0 );
    const double reflected_value = objective.value( reflected );

    if( (reflected_value < values[values.size()-2]) && (reflected_value >= best) )
    {
      vertices.back() = reflected;
      values.back() = reflected_value;
      continue;
    }

    if( reflected_value < best )
    {
      const Eigen::VectorXd expanded = try_point( 2.0 );
      const double expanded_value = objective.value( expanded );
      if( expanded_value < reflected_value )
      {
        vertices.back() = expanded;
        values.back() = expanded_value;
      }else
      {
        vertices.back() = reflected;
        values.back() = reflected_value;
      }
      continue;
    }

    const Eigen::VectorXd contracted = try_point( -0.5 );
    const double contracted_value = objective.value( contracted );
    if( contracted_value < worst )
    {
      vertices.back() = contracted;
      values.back() = contracted_value;
      continue;
    }

    // Shrink toward the best vertex.
    for( size_t i = 1; i < vertices.size(); ++i )
    {
      vertices[i] = vertices.front() + 0.5*(vertices[i] - vertices.front());
      values[i] = objective.value( vertices[i] );
      if( !std::isfinite(values[i]) )
        return false;
    }
  }//for( size_t iteration = 0; iteration < max_iterations; ++iteration )

  reorder();
  if( !std::isfinite(values.front()) )
    return false;

  coefs = vertices.front();
  objective_value = values.front();

  return true;
}//minimize_poisson_simplex(...)


/** An input ROI resolved onto channel boundaries, after any overlapping ROIs have been combined.

 Field names deliberately mirror `DeconRoiInfo` so the peak-construction code does not care which
 it is looking at.
 */
struct MergedRoi
{
  size_t first_channel = 0;   //!< inclusive
  size_t last_channel = 0;    //!< inclusive
  float roi_start = 0.0f;     //!< lower edge of #first_channel
  float roi_end = 0.0f;       //!< upper edge of #last_channel
  PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;
  DeconContinuumNorm cont_norm_method = DeconContinuumNorm::Floating;
  size_t num_lower_side_channels = 0;
  size_t num_upper_side_channels = 0;
  std::vector<DeconRoiInfo::PeakInfo> peak_infos;
};//struct MergedRoi


/** Resolves every input ROI onto channel boundaries and combines any that share channels.

 Summing a per-channel likelihood over ROIs that overlap counts the same data twice and
 manufactures information that is not there - two identical ROIs would appear to improve a limit
 by sqrt(2).  Nothing upstream prevents it: the Detection Confidence Tool gives each gamma line its
 own ROI of mean +- 1.25 FWHM, so any doublet closer than 2.5 FWHM overlaps.

 Overlaps are merged rather than rejected, so that an input which used to produce a (wrong) answer
 still produces an answer - just the correct joint one.  Each merge is described in \p warnings.
 */
std::vector<MergedRoi> merge_overlapping_rois( const std::vector<DeconRoiInfo> &roi_info,
                                              const std::shared_ptr<const SpecUtils::Measurement> &measurement,
                                              std::vector<std::string> &warnings )
{
  const std::shared_ptr<const SpecUtils::EnergyCalibration> cal = measurement->energy_calibration();
  const size_t num_channels = measurement->num_gamma_channels();

  std::vector<MergedRoi> answer;

  for( const DeconRoiInfo &roi : roi_info )
  {
    if( !(roi.roi_end > roi.roi_start) )
      throw runtime_error( "decon_compute_peaks: ROI upper energy is not above its lower energy" );

    const std::pair<size_t,size_t> channels = round_roi_to_channels( measurement, roi.roi_start,
                                                                    roi.roi_end );
    if( (channels.second < channels.first) || (channels.second >= num_channels) )
      throw runtime_error( "decon_compute_peaks: ROI does not lie within the spectrum" );

    MergedRoi entry;
    entry.first_channel = channels.first;
    entry.last_channel = channels.second;
    entry.continuum_type = roi.continuum_type;
    entry.cont_norm_method = roi.cont_norm_method;
    entry.num_lower_side_channels = roi.num_lower_side_channels;
    entry.num_upper_side_channels = roi.num_upper_side_channels;
    entry.peak_infos = roi.peak_infos;

    answer.push_back( entry );
  }//for( const DeconRoiInfo &roi : roi_info )

  std::sort( begin(answer), end(answer),
            []( const MergedRoi &lhs, const MergedRoi &rhs ){
              return lhs.first_channel < rhs.first_channel;
            } );

  std::vector<MergedRoi> merged;
  for( const MergedRoi &entry : answer )
  {
    // Sharing an edge is not overlapping; only sharing a channel is.
    if( merged.empty() || (entry.first_channel > merged.back().last_channel) )
    {
      merged.push_back( entry );
      continue;
    }

    MergedRoi &target = merged.back();

    char buffer[256] = { '\0' };
    snprintf( buffer, sizeof(buffer),
             "Regions of interest %.1f-%.1f keV and %.1f-%.1f keV overlap, and were combined into"
             " one region so their shared channels are not counted twice.",
             cal->energy_for_channel(target.first_channel),
             cal->energy_for_channel(target.last_channel + 1),
             cal->energy_for_channel(entry.first_channel),
             cal->energy_for_channel(entry.last_channel + 1) );
    warnings.push_back( buffer );

    target.last_channel = (std::max)( target.last_channel, entry.last_channel );
    target.num_lower_side_channels = (std::max)( target.num_lower_side_channels,
                                                entry.num_lower_side_channels );
    target.num_upper_side_channels = (std::max)( target.num_upper_side_channels,
                                                entry.num_upper_side_channels );

    // Take the more flexible of the two.  Only compare the polynomial types by value - their enum
    //  values equal their parameter counts, so larger really is more flexible there - because the
    //  step and External types sit *above* them in the enum without that meaning holding.  Picking
    //  External here, for instance, would give `num_linear_fit_pars() == 0` and fail the fit.
    const bool target_is_polynomial = (target.continuum_type >= PeakContinuum::Constant)
                                      && (target.continuum_type <= PeakContinuum::Cubic);
    const bool entry_is_polynomial = (entry.continuum_type >= PeakContinuum::Constant)
                                     && (entry.continuum_type <= PeakContinuum::Cubic);
    if( target_is_polynomial && entry_is_polynomial )
      target.continuum_type = (std::max)( target.continuum_type, entry.continuum_type );
    else if( !target_is_polynomial && entry_is_polynomial )
      target.continuum_type = entry.continuum_type;

    if( target.cont_norm_method != entry.cont_norm_method )
    {
      target.cont_norm_method = DeconContinuumNorm::Floating;
      warnings.push_back( "The combined region of interest was given more than one continuum"
                         " normalization; the continuum was floated for it." );
    }

    target.peak_infos.insert( end(target.peak_infos),
                             begin(entry.peak_infos), end(entry.peak_infos) );
  }//for( const MergedRoi &entry : answer )

  for( MergedRoi &entry : merged )
  {
    entry.roi_start = cal->energy_for_channel( entry.first_channel );
    entry.roi_end = cal->energy_for_channel( entry.last_channel + 1 );
  }

  return merged;
}//merge_overlapping_rois(...)


/** Describes a contiguous run of a measurement's channels as independent Poisson observations.

 \p fixed_signal, when non-empty, must have one entry per channel.  \p continuum_scale is 1 for
 channels of the spectrum the continuum is defined on, and `T_s/T_b` for a projected measurement.
 */
std::vector<PoissonChannel>
measurement_channels( const std::shared_ptr<const SpecUtils::Measurement> &measurement,
                     const size_t first_channel,
                     const size_t last_channel,
                     const std::vector<double> &fixed_signal,
                     const double continuum_scale )
{
  assert( last_channel >= first_channel );
  const std::shared_ptr<const std::vector<float>> energies = measurement->channel_energies();
  const size_t nbin = 1 + last_channel - first_channel;
  assert( fixed_signal.empty() || (fixed_signal.size() == nbin) );

  std::vector<PoissonChannel> channels( nbin );
  for( size_t i = 0; i < nbin; ++i )
  {
    channels[i].lower_energy = static_cast<double>( energies->at(first_channel + i) );
    channels[i].upper_energy = static_cast<double>( energies->at(first_channel + i + 1) );
    channels[i].observed
        = (std::max)( 0.0, static_cast<double>(measurement->gamma_channel_content(first_channel + i)) );
    channels[i].fixed_signal = fixed_signal.empty() ? 0.0 : fixed_signal[i];
    channels[i].continuum_scale = continuum_scale;
  }//for( size_t i = 0; i < nbin; ++i )

  return channels;
}//measurement_channels(...)


/** The channels just below and just above an ROI, as the control observations that
 `DeconContinuumNorm::FixedByEdges` anchors its continuum with.

 The two blocks are disjoint - the ROI sits in the gap between them - which is exactly why the
 optimizer takes an explicit channel list rather than a contiguous range.  Requests that would run
 off either end of the spectrum are clipped rather than wrapping; `first_channel` is unsigned, so
 subtracting a side-channel count from an ROI near channel zero would otherwise underflow.
 */
std::vector<PoissonChannel>
sideband_channels( const std::shared_ptr<const SpecUtils::Measurement> &measurement,
                  const size_t first_channel,
                  const size_t last_channel,
                  const size_t num_lower,
                  const size_t num_upper,
                  const std::vector<std::pair<size_t,size_t>> &other_roi_channels )
{
  const size_t num_channels = measurement->num_gamma_channels();
  std::vector<PoissonChannel> channels;

  // A side channel that lies inside another region of interest is already contributing to the
  //  summed statistic as that region's Poisson data.  Letting it also enter here, through this
  //  region's constraint, counts the same counts twice - the very double counting
  //  `merge_overlapping_rois` exists to prevent, just arriving by a different route.  Reachable
  //  whenever two gamma lines sit closer together than the ROI widths plus their side channels.
  const auto in_another_roi = [&other_roi_channels]( const size_t channel ) -> bool {
    for( const std::pair<size_t,size_t> &roi : other_roi_channels )
    {
      if( (channel >= roi.first) && (channel <= roi.second) )
        return true;
    }
    return false;
  };

  const auto append_block = [&]( const size_t from, const size_t to ){
    for( size_t channel = from; channel <= to; ++channel )
    {
      if( in_another_roi( channel ) )
        continue;

      const std::vector<PoissonChannel> one
            = measurement_channels( measurement, channel, channel, {}, 1.0 );
      channels.insert( end(channels), begin(one), end(one) );
    }
  };

  if( (num_lower > 0) && (first_channel > 0) )
  {
    const size_t lower_end = first_channel - 1;
    const size_t lower_start = (first_channel > num_lower) ? (first_channel - num_lower) : 0;
    append_block( lower_start, lower_end );
  }

  if( (num_upper > 0) && ((last_channel + 1) < num_channels) )
  {
    const size_t upper_start = last_channel + 1;
    // Clamp the *count* rather than the sum, so an absurd `num_upper` from a corrupted saved state
    //  cannot wrap `last_channel + num_upper` around and produce an inverted range.
    const size_t available = num_channels - 1 - last_channel;
    const size_t upper_end = last_channel + (std::min)( num_upper, available );
    append_block( upper_start, upper_end );
  }

  return channels;
}//sideband_channels(...)


/** Sums the Poisson deviance over \p channels, given the continuum coefficients.

 Takes the same `PoissonChannel` list that `fit_continuum_poisson(...)` was handed, so the fit and
 the score cannot describe different data - which is the whole point of the Stage 1 rework.

 Deliberately evaluates the continuum straight from \p continuum_coefficients rather than through
 `PeakContinuum::offset_integral(...)`: that routine clamps a negative continuum to zero per
 channel, which would put a kink in the likelihood.  The optimizer keeps the expected counts
 positive anyway, so the two agree - but only the direct evaluation is guaranteed to.
 */
double decon_roi_statistic( const std::vector<PoissonChannel> &channels,
                           const double reference_energy,
                           const std::vector<double> &continuum_coefficients,
                           const double min_expected )
{
  const size_t nbin = channels.size();
  if( !nbin )
    return 0.0;

  const Eigen::MatrixXd basis = continuum_basis( channels.data(), nbin,
                                                continuum_coefficients.size(), reference_energy );

  double answer = 0.0;
  for( size_t i = 0; i < nbin; ++i )
  {
    double expected = channels[i].fixed_signal;
    for( size_t k = 0; k < continuum_coefficients.size(); ++k )
      expected += continuum_coefficients[k] * basis( static_cast<Eigen::Index>(i),
                                                    static_cast<Eigen::Index>(k) );

    // A fixed or externally supplied continuum can still dip non-positive; the floor keeps the
    //  statistic finite and is far below any physically meaningful count.
    expected = (std::max)( expected, min_expected );

    answer += poisson_deviance_fast( channels[i].observed, expected );
  }//for( size_t i = 0; i < nbin; ++i )

  return answer;
}//decon_roi_statistic(...)

}//anonymous namespace


double poisson_deviance( const double observed, const double expected )
{
  if( IsNan(observed) || IsInf(observed) || IsNan(expected) || IsInf(expected) )
    throw runtime_error( "poisson_deviance: non-finite input" );

  if( observed < 0.0 )
    throw runtime_error( "poisson_deviance: observed counts must be non-negative" );

  if( expected <= 0.0 )
    throw runtime_error( "poisson_deviance: expected counts must be strictly positive" );

  return poisson_deviance_fast( observed, expected );
}//poisson_deviance(...)


double min_expected_channel_counts( const double mean_observed_counts )
{
  return sm_min_expected_fraction * (std::max)( 1.0, mean_observed_counts );
}//min_expected_channel_counts(...)


PoissonContinuumFit fit_continuum_poisson( const float * const channel_edges,
                                          const double * const observed,
                                          const double * const fixed_signal,
                                          const size_t nbin,
                                          const size_t num_coefficients,
                                          const double reference_energy,
                                          const vector<double> &initial_coefficients,
                                          const vector<double> &constraint_center,
                                          const vector<double> &constraint_inverse_covariance )
{
  if( !channel_edges || !observed )
    throw runtime_error( "fit_continuum_poisson: null input array" );

  if( nbin < 1 )
    throw runtime_error( "fit_continuum_poisson: need at least one channel and one coefficient" );

  vector<PoissonChannel> channels( nbin );
  for( size_t i = 0; i < nbin; ++i )
  {
    channels[i].lower_energy = static_cast<double>( channel_edges[i] );
    channels[i].upper_energy = static_cast<double>( channel_edges[i+1] );
    channels[i].observed = observed[i];
    channels[i].fixed_signal = fixed_signal ? fixed_signal[i] : 0.0;
    channels[i].continuum_scale = 1.0;
  }

  return fit_continuum_poisson( channels.data(), nbin, num_coefficients, reference_energy,
                               initial_coefficients, constraint_center,
                               constraint_inverse_covariance );
}//fit_continuum_poisson(...)


PoissonContinuumFit fit_continuum_poisson( const PoissonChannel * const channels,
                                          const size_t nbin,
                                          const size_t num_coefficients,
                                          const double reference_energy,
                                          const vector<double> &initial_coefficients,
                                          const vector<double> &constraint_center,
                                          const vector<double> &constraint_inverse_covariance )
{
  if( !channels )
    throw runtime_error( "fit_continuum_poisson: null input array" );

  if( (nbin < 1) || (num_coefficients < 1) )
    throw runtime_error( "fit_continuum_poisson: need at least one channel and one coefficient" );

  if( num_coefficients > nbin )
    throw runtime_error( "fit_continuum_poisson: more coefficients than channels" );

  if( IsNan(reference_energy) || IsInf(reference_energy) )
    throw runtime_error( "fit_continuum_poisson: invalid reference energy" );

  if( !initial_coefficients.empty() && (initial_coefficients.size() != num_coefficients) )
    throw runtime_error( "fit_continuum_poisson: initial coefficient size mismatch" );

  const bool constrained = !constraint_center.empty();
  if( constrained
     && ((constraint_center.size() != num_coefficients)
         || (constraint_inverse_covariance.size() != num_coefficients*num_coefficients)) )
    throw runtime_error( "fit_continuum_poisson: constraint size mismatch" );

  if( !constrained && !constraint_inverse_covariance.empty() )
    throw runtime_error( "fit_continuum_poisson: constraint precision given without a center" );

  const Eigen::Index npar = static_cast<Eigen::Index>( num_coefficients );
  const Eigen::Index nchan = static_cast<Eigen::Index>( nbin );

  PoissonContinuumFit answer;
  answer.coefficients.assign( num_coefficients, 0.0 );

  for( size_t i = 0; i < nbin; ++i )
  {
    const PoissonChannel &channel = channels[i];

    // Only an *inverted* channel is rejected.  A zero-width channel is unhelpful but legal:
    //  `SpecUtils::EnergyCalibration` permits equal consecutive lower-channel energies (its check
    //  rejects only `edge[i] < edge[i-1]`), and some `LowerChannelEdge` spectra pad their tail
    //  that way.  Treating equality as invalid here would make spectra that previously produced a
    //  limit start throwing out of the whole scan.
    if( (channel.upper_energy < channel.lower_energy)
       || IsNan(channel.lower_energy) || IsInf(channel.lower_energy)
       || IsNan(channel.upper_energy) || IsInf(channel.upper_energy) )
      throw runtime_error( "fit_continuum_poisson: invalid channel energy range" );

    // A non-positive scale would flip the sign of basis column zero, which `restore_positivity`
    //  relies on being positive to shift the constant term into the feasible region.
    if( !(channel.continuum_scale > 0.0) || IsInf(channel.continuum_scale) )
      throw runtime_error( "fit_continuum_poisson: continuum scale must be strictly positive" );
  }//for( size_t i = 0; i < nbin; ++i )

  // Gather the data, and the equilibration scale for each basis column.
  const Eigen::MatrixXd raw_basis = continuum_basis( channels, nbin, num_coefficients,
                                                    reference_energy );
  if( !raw_basis.allFinite() )
    throw runtime_error( "fit_continuum_poisson: non-finite channel energies" );

  PoissonContinuumObjective objective;
  objective.m_observed.resize( nchan );
  objective.m_fixed_signal = Eigen::VectorXd::Zero( nchan );

  double observed_sum = 0.0, fixed_sum = 0.0;
  for( size_t i = 0; i < nbin; ++i )
  {
    const double counts = channels[i].observed;
    if( IsNan(counts) || IsInf(counts) || (counts < 0.0) )
      throw runtime_error( "fit_continuum_poisson: invalid observed counts" );

    objective.m_observed( static_cast<Eigen::Index>(i) ) = counts;
    observed_sum += counts;

    const double signal = channels[i].fixed_signal;
    if( IsNan(signal) || IsInf(signal) || (signal < 0.0) )
      throw runtime_error( "fit_continuum_poisson: invalid fixed signal counts" );
    objective.m_fixed_signal( static_cast<Eigen::Index>(i) ) = signal;
    fixed_sum += signal;
  }//for( size_t i = 0; i < nbin; ++i )

  const double mean_observed = observed_sum / static_cast<double>(nbin);
  objective.m_min_expected = min_expected_channel_counts( mean_observed );

  objective.m_scale.resize( npar );
  objective.m_basis.resize( nchan, npar );
  for( Eigen::Index k = 0; k < npar; ++k )
  {
    // Scale each column to unit RMS.  Expected counts must be unchanged, so
    //  B(:,k)*c_raw == (B(:,k)/rms)*c_scaled, hence c_raw == c_scaled/rms, which is what
    //  `m_scale` records.
    const double rms = std::sqrt( raw_basis.col(k).squaredNorm() / static_cast<double>(nbin) );
    const double inverse_rms = (rms > 0.0) ? (1.0/rms) : 1.0;
    objective.m_scale(k) = inverse_rms;
    objective.m_basis.col(k) = inverse_rms * raw_basis.col(k);
  }

  if( constrained )
  {
    // (c_raw - m)^T P (c_raw - m) with c_raw = S*c_scaled becomes
    //  (c_scaled - S^-1 m)^T (S P S) (c_scaled - S^-1 m).
    Eigen::VectorXd center( npar );
    for( Eigen::Index k = 0; k < npar; ++k )
      center(k) = constraint_center[static_cast<size_t>(k)] / objective.m_scale(k);

    Eigen::MatrixXd precision( npar, npar );
    for( Eigen::Index r = 0; r < npar; ++r )
    {
      for( Eigen::Index c = 0; c < npar; ++c )
      {
        const size_t index = static_cast<size_t>(r)*num_coefficients + static_cast<size_t>(c);
        precision(r,c) = constraint_inverse_covariance[index]
                         * objective.m_scale(r) * objective.m_scale(c);
      }
    }

    if( !center.allFinite() || !precision.allFinite() )
      throw runtime_error( "fit_continuum_poisson: non-finite constraint" );

    objective.m_constraint_center = center;
    objective.m_constraint_precision = precision;
  }//if( constrained )

  // Build the ladder of starting estimates, best first.
  vector<Eigen::VectorXd> seeds;

  if( !initial_coefficients.empty() )
  {
    Eigen::VectorXd seed( npar );
    for( Eigen::Index k = 0; k < npar; ++k )
      seed(k) = initial_coefficients[static_cast<size_t>(k)] / objective.m_scale(k);
    if( seed.allFinite() )
      seeds.push_back( seed );
  }

  {// Weighted least squares on `observed - fixed_signal`, with no clipping of that difference.
    Eigen::VectorXd weights( nchan );
    for( Eigen::Index i = 0; i < nchan; ++i )
      weights(i) = 1.0 / (std::max)( objective.m_observed(i), 1.0 );

    const Eigen::VectorXd residual = objective.m_observed - objective.m_fixed_signal;
    const Eigen::MatrixXd normal = objective.m_basis.transpose()
                                   * weights.asDiagonal() * objective.m_basis;
    const Eigen::VectorXd rhs = objective.m_basis.transpose()
                                * weights.asDiagonal() * residual;
    const Eigen::VectorXd seed
                       = normal.completeOrthogonalDecomposition().solve( rhs );
    if( seed.allFinite() )
      seeds.push_back( seed );
  }

  {// A flat continuum holding the observed excess; the most robust fallback there is.
    Eigen::VectorXd seed = Eigen::VectorXd::Zero( npar );
    const double basis_sum = objective.m_basis.col(0).sum();
    if( basis_sum > 0.0 )
    {
      seed(0) = (observed_sum - fixed_sum) / basis_sum;
      seeds.push_back( seed );

      // ...and a deliberately continuum-dominated variant, for when the trial signal is so large
      //  that the excess-preserving seed starts out at or below zero.
      Eigen::VectorXd low_seed = Eigen::VectorXd::Zero( npar );
      low_seed(0) = 0.5*observed_sum / basis_sum;
      seeds.push_back( low_seed );
    }
  }

  string last_error = "the continuum optimizer did not converge";

  // Every seed is run to completion and the best deviance kept, rather than accepting the first
  //  one that converges.  On a constrained problem different starting points can stop at genuinely
  //  different places, and the extra cost is a handful of Newton solves on a 2x2 or 3x3 system.
  bool have_best = false;
  Eigen::VectorXd best_coefs;
  double best_deviance = std::numeric_limits<double>::infinity();
  size_t best_restarts = 0;

  // When every channel saw at least one count the deviance is its own barrier - `-2*n*log(E)`
  //  diverges as any expected count falls to zero - so the minimum is strictly interior, the
  //  objective is strictly convex, and plain Newton finds the unique answer.  The positivity
  //  constraint can only become active if some channel is empty, and only then is the (several
  //  times more expensive) barrier path needed.  This is the common case for real spectra.
  const bool self_barriering = (objective.m_observed.minCoeff() > 0.0);

  for( size_t seed_index = 0; seed_index < seeds.size(); ++seed_index )
  {
    Eigen::VectorXd coefs = seeds[seed_index];
    if( !objective.restore_positivity( coefs, objective.interior_target() ) )
    {
      last_error = "could not find continuum coefficients giving positive expected counts";
      continue;
    }

    bool converged = false;
    size_t stage_restarts = 0;

    if( self_barriering )
    {
      double value = 0.0;
      converged = minimize_poisson_newton( objective, coefs, value, answer.num_iterations );
      if( !converged )
      {
        converged = minimize_poisson_simplex( objective, coefs, value, answer.num_iterations );
        stage_restarts = 1;
      }
    }else
    {
      // Sequence of barrier subproblems with decreasing mu.  The duality gap is bounded by
      //  nbin*mu, so stopping once that is far below the confidence threshold (2.7055) leaves the
      //  solution accurate where it matters.  Each subproblem starts from the previous solution,
      //  so the later, harder ones need only a few iterations.
      const double initial_deviance = objective.deviance( coefs );
      double mu = (std::max)( 1.0, std::isfinite(initial_deviance) ? fabs(initial_deviance) : 1.0 )
                  / static_cast<double>(nbin);
      const double min_mu = 1.0E-6 / static_cast<double>(nbin);

      while( true )
      {
        objective.m_barrier_mu = mu;

        double value = 0.0;
        bool stage_ok = minimize_poisson_newton( objective, coefs, value, answer.num_iterations );
        if( !stage_ok )
        {
          // The barrier subproblem is smooth and convex, so a Newton failure here means numerical
          //  trouble rather than an active constraint; let the derivative-free method try.
          stage_ok = minimize_poisson_simplex( objective, coefs, value, answer.num_iterations );
          stage_restarts = 1;
        }

        converged = converged || stage_ok;
        if( !stage_ok || (mu <= min_mu) )
          break;

        mu *= 0.1;
      }//while( barrier stages remain )

      objective.m_barrier_mu = 0.0;
    }//if( self_barriering ) / else

    const double deviance = objective.deviance( coefs );
    if( !converged || !std::isfinite(deviance) )
    {
      last_error = "the continuum optimizer did not converge";
      continue;
    }

    if( !have_best || (deviance < best_deviance) )
    {
      have_best = true;
      best_coefs = coefs;
      best_deviance = deviance;
      best_restarts = seed_index + stage_restarts;
    }

    // A strictly convex problem has one minimum, so there is nothing for the remaining starting
    //  estimates to find.  With an active constraint they can genuinely differ, so keep going.
    if( self_barriering )
      break;
  }//for( size_t seed_index = 0; seed_index < seeds.size(); ++seed_index )

  if( !have_best )
  {
    answer.converged = false;
    answer.error = last_error;
    answer.statistic = std::numeric_limits<double>::infinity();
    return answer;
  }

  answer.num_restarts = best_restarts;
  answer.converged = true;
  answer.statistic = best_deviance;

  for( Eigen::Index k = 0; k < npar; ++k )
    answer.coefficients[static_cast<size_t>(k)] = best_coefs(k) * objective.m_scale(k);

  // Covariance of the raw coefficients, from the expected information plus any constraint.
  const Eigen::MatrixXd information = objective.information( best_coefs );

  // With `c_raw = S*c_scaled`, `Cov_raw = S*Cov_scaled*S` and so `F_raw = S^-1*F_scaled*S^-1`.
  if( information.allFinite() )
  {
    answer.information.assign( num_coefficients*num_coefficients, 0.0 );
    for( Eigen::Index r = 0; r < npar; ++r )
    {
      for( Eigen::Index c = 0; c < npar; ++c )
      {
        const size_t index = static_cast<size_t>(r)*num_coefficients + static_cast<size_t>(c);
        answer.information[index] = information(r,c) / (objective.m_scale(r) * objective.m_scale(c));
      }
    }
  }//if( information.allFinite() )

  const Eigen::MatrixXd inverse = information.completeOrthogonalDecomposition().pseudoInverse();
  if( inverse.allFinite() )
  {
    answer.covariance.assign( num_coefficients*num_coefficients, 0.0 );
    for( Eigen::Index r = 0; r < npar; ++r )
    {
      for( Eigen::Index c = 0; c < npar; ++c )
      {
        const size_t index = static_cast<size_t>(r)*num_coefficients + static_cast<size_t>(c);
        answer.covariance[index] = inverse(r,c) * objective.m_scale(r) * objective.m_scale(c);
      }
    }
  }//if( inverse.allFinite() )

  return answer;
}//fit_continuum_poisson(...)


DeconRoiInfo::DeconRoiInfo()
: roi_start( 0.0f ),
  roi_end( 0.0f ),
  continuum_type( PeakContinuum::OffsetType::NoOffset ),
  cont_norm_method( DeconContinuumNorm::Floating ),
  num_lower_side_channels( 0 ),
  num_upper_side_channels( 0 ),
  peak_infos()
{
}


DeconRoiInfo::PeakInfo::PeakInfo()
 : energy( 0.0f ),
   fwhm( 0.0f ),
   counts_per_bq_into_4pi( 0.0 )
{
}



DeconComputeInput::DeconComputeInput()
 : distance( 0.0 ),
   activity( 0.0 ),
   include_air_attenuation( false ),
   shielding_thickness( 0.0 ),
   drf( nullptr ),
   measurement( nullptr ),
   roi_info()
{
}


DeconComputeResults decon_compute_peaks( const DeconComputeInput &input )
{
  DeconComputeResults result;
  result.input = input;
  result.chi2 = 0.0;
  result.statistic_name = "Cash";
  result.num_degree_of_freedom = 0;
  result.num_continuum_iterations = 0;
  result.num_continuum_restarts = 0;

  // Lets sanity check input
  if( (input.distance < 0.0) || IsNan(input.distance) || IsInf(input.distance) )
    throw runtime_error( "decon_compute_peaks: invalid input distance" );
  
  // Lets sanity check input
  if( (input.activity < 0.0) || IsNan(input.activity) || IsInf(input.activity) )
    throw runtime_error( "decon_compute_peaks: invalid input activity" );
  
  if( input.include_air_attenuation
      && ((input.shielding_thickness < 0.0)
          || IsNan(input.shielding_thickness)
          || IsInf(input.shielding_thickness)
          || (input.shielding_thickness >= input.distance)) )
    throw runtime_error( "decon_compute_peaks: invalid input shielding thickness" );
  
  if( !input.drf || !input.drf->isValid() || !input.drf->hasResolutionInfo() )
    throw runtime_error( "decon_compute_peaks: invalid DRF input" );
  
  
  if( !input.measurement || (input.measurement->num_gamma_channels() < 2)
     || !input.measurement->energy_calibration() //ptr should always be valid anyway
     || !input.measurement->energy_calibration()->valid() )
    throw runtime_error( "decon_compute_peaks: invalid spectrum input" );
  
  // Check if there is anything to do
  if( input.roi_info.empty() )
    return result;
    
  const bool fixed_geom = input.drf->isFixedGeometry();

  double total_chi2 = 0.0, total_DOF = 0.0;
  vector<PeakDef> inputPeaks, fittedPeaks;

  // A sample that has actually been measured replaces the expected-counts block; only meaningful under a
  //  background reference, and ignored (rather than quietly changing the answer) otherwise.  Its
  //  channels are read at the reference's channel indices, so the two must share a binning.
  const shared_ptr<const SpecUtils::Measurement> observed_sample
      = (input.measurement_model == DeconMeasurementModel::BackgroundReference)
          ? input.observed_sample : nullptr;

  if( observed_sample )
  {
    if( observed_sample->num_gamma_channels() != input.measurement->num_gamma_channels() )
      throw runtime_error( "decon_compute_peaks: observed sample and reference have different"
                          " numbers of channels" );

    const shared_ptr<const SpecUtils::EnergyCalibration> sample_cal
                                                    = observed_sample->energy_calibration();
    if( !sample_cal || !sample_cal->valid()
       || ((*sample_cal) != (*input.measurement->energy_calibration())) )
      throw runtime_error( "decon_compute_peaks: observed sample and reference have different"
                          " energy calibrations" );

    if( !(observed_sample->live_time() > 0.0f) )
      throw runtime_error( "decon_compute_peaks: observed sample has no live time" );
  }//if( observed_sample )

  // Depends only on the input, so it is settled here rather than inside the region loop: an empty
  //  region list, or a region with no peaks, would otherwise leave the field reading 1.0 under a
  //  background reference and quietly mislead anything that consumes it for something other than
  //  scaling peaks.  The loop's own `exposure_ratio` re-derives it with full validation.
  double overall_exposure_ratio = 1.0;
  if( (input.measurement_model == DeconMeasurementModel::BackgroundReference)
     && input.measurement
     && (input.measurement->live_time() > 0.0f) )
  {
    // A measured sample carries its own exposure, which then supersedes `sample_exposure` - the
    //  latter describes a measurement that has not been taken.
    const double sample_live = observed_sample ? static_cast<double>( observed_sample->live_time() )
                                               : input.sample_exposure;
    if( sample_live > 0.0 )
      overall_exposure_ratio = sample_live / static_cast<double>( input.measurement->live_time() );
  }

  // Resolve the ROIs onto channel boundaries once, and combine any that share channels, so that no
  //  channel contributes to the summed likelihood more than once.
  const vector<MergedRoi> merged_rois = merge_overlapping_rois( input.roi_info, input.measurement,
                                                               result.warnings );

  for( const MergedRoi &roi : merged_rois )
  {
    // A region with no gamma lines carries no information about the activity, and there would be
    //  no continuum object to score against.  Skipping it matches what the previous
    //  implementation effectively did (its chi2/DOF helper returned zero for an empty peak list).
    if( roi.peak_infos.empty() )
    {
      result.warnings.push_back( "A region of interest was given no gamma lines, and was skipped." );
      continue;
    }

    // Both the continuum fit and the score use exactly this channel range; they used to differ by
    //  one channel, which made the fitted continuum minimize a slightly different sum than the one
    //  reported.
    const size_t lower_channel = roi.first_channel;
    const size_t upper_channel = roi.last_channel;
    const float &roi_start = roi.roi_start;
    const float &roi_end = roi.roi_end;

    const DeconContinuumNorm &cont_norm_method = roi.cont_norm_method;
    PeakContinuum::OffsetType continuum_type = roi.continuum_type;
    switch( cont_norm_method )
    {
      case DeconContinuumNorm::Floating:
      case DeconContinuumNorm::FixedByFullRange:
        break;
        
      case DeconContinuumNorm::FixedByEdges:
        // Forced to Linear because the sideband constraint is only as good as the sidebands: a
        //  typical 4+4 side channels barely determine an offset and a slope, and asking them for a
        //  curvature term as well gives a constraint so weak that the "anchored" continuum is
        //  effectively floating again.  The constraint code itself is generic in the number of
        //  coefficients, so this is a defensible default rather than a limitation.
        assert( continuum_type == PeakContinuum::OffsetType::Linear );
        continuum_type = PeakContinuum::OffsetType::Linear;
        break;
    }//switch( cont_norm_method )
    
    const size_t &num_lower_side_channels = roi.num_lower_side_channels;
    const size_t &num_upper_side_channels = roi.num_upper_side_channels;

    const bool background_reference
                = (input.measurement_model == DeconMeasurementModel::BackgroundReference);

    // Exposure of the sample measurement, relative to the reference that was loaded.  A
    //  non-positive `sample_exposure` means "the same length again", giving a ratio of one.
    //  `PeakInfo::counts_per_bq_into_4pi` already carries the loaded spectrum's live time, so the
    //  trial signal has to be rescaled by this to describe the sample measurement.
    const double reference_exposure = static_cast<double>( input.measurement->live_time() );
    double exposure_ratio = 1.0;
    if( background_reference )
    {
      if( !(reference_exposure > 0.0) )
        throw runtime_error( "decon_compute_peaks: background reference spectrum has no live time" );

      // A sample that has been measured carries its own exposure; `sample_exposure` describes one
      //  that has not been taken, so the real thing supersedes it.
      const double sample_exposure = observed_sample
                                       ? static_cast<double>( observed_sample->live_time() )
                                       : input.sample_exposure;
      if( sample_exposure > 0.0 )
        exposure_ratio = sample_exposure / reference_exposure;

      if( !(exposure_ratio > 0.0) || IsInf(exposure_ratio) )
        throw runtime_error( "decon_compute_peaks: invalid sample exposure" );
    }//if( background_reference )

    overall_exposure_ratio = exposure_ratio;

    shared_ptr<PeakContinuum> peak_continuum;
    shared_ptr<const SpecUtils::Measurement> computed_global_cont;
    
    // Find the largest peak in the ROIs, energy to use as the continuum "reference energy"
    float reference_energy = 0.0f;
    for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
      reference_energy = std::max( reference_energy, peak_info.energy );
    
    assert( reference_energy != 0.0f );

    vector<PeakDef> roi_peaks;

    for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
    {
      const float &energy = peak_info.energy;
      const float fwhm = (peak_info.fwhm > 0.0f) ? peak_info.fwhm : input.drf->peakResolutionFWHM( energy );
      // FWHM = 2*sqrt(2*ln(2))*sigma; see `PeakDef::fwhm()`.
      const float sigma = fwhm / PhysicalUnits::fwhm_nsigma;
      const double det_eff = fixed_geom ? input.drf->intrinsicEfficiency(energy)
                                        : input.drf->efficiency( energy, input.distance );
      const double counts_4pi = peak_info.counts_per_bq_into_4pi;
      double air_atten = 1.0;
      
      if( input.include_air_attenuation && !fixed_geom )
      {
        const double air_len = input.distance - input.shielding_thickness;
        const double mu_air = GammaInteractionCalc::transmission_coefficient_air( energy, air_len );
        air_atten = exp( -1.0 * mu_air );
      }
      
      // `counts_4pi` carries the *loaded* spectrum's live time, so predicting a measurement of a
      //  different length means rescaling the signal by the exposure ratio; it is exactly 1 for
      //  the ordinary current-spectrum case.  Scaling the peak itself, rather than only the
      //  channel counts, keeps the returned peaks describing the measurement being reported on.
      const float amplitude
              = static_cast<float>( exposure_ratio * input.activity * counts_4pi * det_eff * air_atten );

      PeakDef peak( energy, sigma, amplitude );
      peak.setFitFor( PeakDef::CoefficientType::Mean, false );
      peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
      peak.setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
      
      if( peak_continuum )
      {
        peak.setContinuum( peak_continuum );
      }else
      {
        peak_continuum = peak.continuum();
        peak_continuum->setType( continuum_type );
        peak_continuum->setRange( roi_start, roi_end );
        peak_continuum->setType( continuum_type );

        // This next loop could be skipped
        for( size_t order = 0; order < peak_continuum->parameters().size(); ++order )
        {
          switch( cont_norm_method )
          {
            case DeconContinuumNorm::Floating:
            case DeconContinuumNorm::FixedByEdges:
              // FixedByEdges profiles its coefficients too - they are anchored by the sidebands'
              //  covariance, not held fixed - so they are fitted parameters in both cases.
              peak_continuum->setPolynomialCoefFitFor( order, true );
              break;

            case DeconContinuumNorm::FixedByFullRange:
              peak_continuum->setPolynomialCoefFitFor( order, false );
              break;
          }//switch( cont_norm_method )
        }//for( size_t order = 0; order < peak_continuum->parameters().size(); ++order )
        
        
        switch( continuum_type )
        {
          case PeakContinuum::Linear:
          case PeakContinuum::Quadratic:
            break;
          
          case PeakContinuum::NoOffset:
          case PeakContinuum::Constant:
          case PeakContinuum::Cubic:
          case PeakContinuum::FlatStep:
          case PeakContinuum::LinearStep:
          case PeakContinuum::BiLinearStep:
          case PeakContinuum::FlatStepCDF:
          case PeakContinuum::LinearStepCDF:
          case PeakContinuum::BiLinearStepCDF:
            assert( 0 );
            break;

          case PeakContinuum::External:
            if( !computed_global_cont )
              computed_global_cont = estimateContinuum( input.measurement );
            peak_continuum->setExternalContinuum( computed_global_cont );
            break;
        }//switch( assign continuum )
      }//if( peak_continuum ) / else


      roi_peaks.push_back( peak );
    }//for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )

    assert( upper_channel >= lower_channel );

    const shared_ptr<const vector<float>> &channel_energies = input.measurement->channel_energies();
    assert( channel_energies );
    const size_t nbin = (1 + upper_channel - lower_channel);

    // Counts each trial peak puts in each channel of the ROI; this is the fixed signal the
    //  continuum is profiled against, and it enters the score unchanged.
    vector<double> peak_channel_counts( nbin, 0.0 );
    for( const PeakDef &p : roi_peaks )
      p.gauss_integral( &(channel_energies->at(lower_channel)), peak_channel_counts.data(), nbin );

    // The channels the statistic is summed over.  Building this once and using it for *both* the
    //  continuum fit and the score is what guarantees the reported number is evaluated at the
    //  continuum that actually minimizes it.
    vector<PoissonChannel> channels;

    const size_t num_coefficients = PeakContinuum::num_linear_fit_pars( continuum_type );

    vector<double> initial_coefficients;
    if( peak_continuum && peak_continuum->parametersProbablySet()
       && (peak_continuum->parameters().size() >= num_coefficients)
       && (fabs(peak_continuum->referenceEnergy() - reference_energy) < 0.001) )
    {
      initial_coefficients.assign( begin(peak_continuum->parameters()),
                                  begin(peak_continuum->parameters()) + num_coefficients );
    }

    // `DeconContinuumNorm::FixedByEdges` uses the channels beside the region as background-control
    //  observations.  They enter the *same* Poisson likelihood as the region's own channels, with
    //  no trial signal in them - an exact joint likelihood rather than a Gaussian summary of a
    //  separate fit.
    //
    //  The earlier implementation fitted the sidebands alone and carried their covariance in as a
    //  penalty.  That is only a second-order approximation of this, and it had a failure mode
    //  precisely where the option is most used: the expected (Fisher) information weights each
    //  channel by 1/E, so an empty side channel driven onto the positivity floor contributed an
    //  enormous 1/E and over-constrained the continuum, biasing limits low below roughly three
    //  counts per channel.  Summing the deviance directly has no information matrix, no 1/E, and
    //  no quadratic approximation, so none of that arises.
    vector<PoissonChannel> sidebands;
    if( cont_norm_method == DeconContinuumNorm::FixedByEdges )
    {
      // Every other region's channels, so they can be kept out of this region's control sample: a
      //  side channel that is already another region's data would otherwise be counted twice.
      vector<pair<size_t,size_t>> other_roi_channels;
      for( const MergedRoi &other : merged_rois )
      {
        if( (other.first_channel != lower_channel) || (other.last_channel != upper_channel) )
          other_roi_channels.emplace_back( other.first_channel, other.last_channel );
      }

      sidebands = sideband_channels( input.measurement, lower_channel, upper_channel,
                                    num_lower_side_channels, num_upper_side_channels,
                                    other_roi_channels );

      if( sidebands.empty() )
      {
        result.warnings.push_back( "There were no usable channels beside the region of interest at "
                          + std::to_string( static_cast<int>(std::round(reference_energy)) )
                          + " keV, so its continuum was left free to float." );
      }
    }//if( FixedByEdges )

    if( background_reference )
    {
      const vector<PoissonChannel> roi_reference
            = measurement_channels( input.measurement, lower_channel, upper_channel, {}, 1.0 );

      // The sidebands are channels of the reference spectrum too, so they join the reference
      //  block.  (The sample would have its own sidebands; they are not modelled either way, which
      //  only makes the answer conservative.)
      vector<PoissonChannel> reference = roi_reference;
      reference.insert( end(reference), begin(sidebands), end(sidebands) );

      channels = reference;
      channels.reserve( reference.size() + nbin );

      if( observed_sample )
      {
        // The sample has actually been measured, so the second block is its real counts and the
        //  two blocks are a genuine joint likelihood.  No expected-counts step, and so no null fit.
        const vector<PoissonChannel> roi_sample
              = measurement_channels( observed_sample, lower_channel, upper_channel,
                                     peak_channel_counts, exposure_ratio );
        channels.insert( end(channels), begin(roi_sample), end(roi_sample) );
      }else
      {
        // The sample has not been measured, so it is taken at the counts it expects, with no measurement noise under
        //  the null.  Solving that null once, here, keeps the expected-counts data fixed across the whole
        //  activity scan - it describes the measurement being predicted, not the activity tested.
        const PoissonContinuumFit null_fit
            = fit_continuum_poisson( reference.data(), reference.size(), num_coefficients,
                                    reference_energy, initial_coefficients );

        if( !null_fit.converged )
          throw runtime_error( "decon_compute_peaks: background reference fit failed - "
                              + null_fit.error );

        result.num_continuum_iterations += null_fit.num_iterations;
        result.num_continuum_restarts += null_fit.num_restarts;

        const Eigen::MatrixXd basis = continuum_basis( roi_reference.data(), nbin, num_coefficients,
                                                      reference_energy );

        for( size_t i = 0; i < nbin; ++i )
        {
          double continuum_counts = 0.0;
          for( size_t k = 0; k < num_coefficients; ++k )
            continuum_counts += null_fit.coefficients[k] * basis( static_cast<Eigen::Index>(i),
                                                                 static_cast<Eigen::Index>(k) );

          PoissonChannel sample = roi_reference[i];
          sample.observed = (std::max)( 0.0, exposure_ratio*continuum_counts );
          sample.fixed_signal = peak_channel_counts[i];
          sample.continuum_scale = exposure_ratio;
          channels.push_back( sample );
        }//for( size_t i = 0; i < nbin; ++i )

        if( initial_coefficients.empty() )
          initial_coefficients = null_fit.coefficients;
      }//if( observed_sample ) / else
    }else
    {
      channels = measurement_channels( input.measurement, lower_channel, upper_channel,
                                      peak_channel_counts, 1.0 );
      channels.insert( end(channels), begin(sidebands), end(sidebands) );
    }//if( background_reference ) / else

    // The positivity floor is relative to the data actually being scored, so it is taken over the
    //  whole channel list rather than the region alone - that keeps the fit and the score using
    //  the same floor.
    double scored_observed_sum = 0.0;
    for( const PoissonChannel &channel : channels )
      scored_observed_sum += channel.observed;

    const double min_expected = min_expected_channel_counts(
                    scored_observed_sum / static_cast<double>(channels.size()) );

    // `channels` always carries the trial signal, because that is what the statistic is summed
    //  over.  FixedByFullRange is the one treatment whose *fit* sees a different model than its
    //  score: it solves the continuum with the signal left out, which is precisely what makes it
    //  an assertion that no signal is present.  Keeping the two lists separate rather than
    //  reusing one is what stops the score losing the signal term - without it the profile is
    //  flat in activity and never crosses the threshold at all.
    //  Only FixedByFullRange needs a second list, so only it pays for the copy; this runs once per
    //  trial activity, of which a scan evaluates a few hundred.
    vector<PoissonChannel> zero_signal_channels;
    if( cont_norm_method == DeconContinuumNorm::FixedByFullRange )
    {
      zero_signal_channels = channels;
      for( PoissonChannel &channel : zero_signal_channels )
        channel.fixed_signal = 0.0;
    }

    const vector<PoissonChannel> &fit_channels
        = (cont_norm_method == DeconContinuumNorm::FixedByFullRange) ? zero_signal_channels
                                                                     : channels;

    // Every continuum treatment solves its coefficients from the ROI channels, once per ROI.
    //  FixedByEdges differs only in carrying a Gaussian penalty from its sidebands, and
    //  FixedByFullRange only in leaving the trial signal out of the model while it solves.  (The
    //  FixedByEdges line used to be re-derived once per gamma line in the ROI, so the last peak's
    //  call silently won.)
    if( peak_continuum )
    {
      // Profile the continuum against the *same* Poisson likelihood the score uses, with the trial
      //  signal held fixed.  The previous weighted-least-squares solve minimized a different
      //  quantity (a modified-Neyman chi-square of `observed - signal`, clipped at zero), so the
      //  continuum it produced was not the minimizer of the statistic reported afterwards.
      const PoissonContinuumFit fit
          = fit_continuum_poisson( fit_channels.data(), fit_channels.size(), num_coefficients,
                                  reference_energy, initial_coefficients );

      if( !fit.converged )
        throw runtime_error( "decon_compute_peaks: continuum fit failed - " + fit.error );

      result.num_continuum_iterations += fit.num_iterations;
      result.num_continuum_restarts += fit.num_restarts;

      vector<double> continuum_coeffs = fit.coefficients;
      vector<double> continuum_coeffs_uncerts( num_coefficients, 0.0 );
      for( size_t k = 0; k < num_coefficients; ++k )
      {
        const size_t index = k*num_coefficients + k;
        if( index < fit.covariance.size() && (fit.covariance[index] > 0.0) )
          continuum_coeffs_uncerts[k] = sqrt( fit.covariance[index] );
      }

      // For FlatStepCDF/LinearStepCDF, `num_linear_fit_pars` excludes the step coefficient that
      //  `setParameters` expects.  BiLinearStepCDF has no step coefficient.
      if( PeakContinuum::is_peak_cdf_step_continuum( continuum_type )
         && (continuum_type != PeakContinuum::BiLinearStepCDF) )
      {
        continuum_coeffs.push_back( 0.0 );
        continuum_coeffs_uncerts.push_back( 0.0 );
      }
      peak_continuum->setType( continuum_type );
      peak_continuum->setParameters( reference_energy, continuum_coeffs, continuum_coeffs_uncerts );
    }//if( peak_continuum )

    for( PeakDef &p : roi_peaks )
    {
      // Set that we arent/didnt fit for Mean, Sigma, or Amplitude so we'll get DOF correct.
      p.setFitFor( PeakDef::CoefficientType::Mean, false );
      p.setFitFor( PeakDef::CoefficientType::Sigma, false );
      p.setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
    }

    const size_t num_continuum_pars = num_coefficients;
    vector<double> scoring_coefficients( begin(peak_continuum->parameters()),
                                        begin(peak_continuum->parameters()) + num_continuum_pars );

    // Scored over exactly the channel list the continuum was fit against - for a background
    //  reference that is the reference block and the projected sample together.
    const double chi2 = decon_roi_statistic( channels, reference_energy, scoring_coefficients,
                                            min_expected );


    // Degrees of freedom are goodness-of-fit metadata only; the confidence threshold is set by the
    //  one activity/distance parameter being profiled, not by this.  Counted as
    //  (channels carrying real counts) - (fitted parameters).
    //
    //  With the sidebands now entering the likelihood as ordinary observations this needs no
    //  special case: FixedByEdges simply has more channels.  The one exclusion is the background
    //  reference's projected-sample block, which is expected-counts pseudo-data built from the reference -
    //  counting it would inflate the DOF and pin chi2/dof near 0.5, because that block contributes
    //  ~zero deviance by construction.
    const size_t real_observations = background_reference ? (nbin + sidebands.size())
                                                          : channels.size();
    double dof = static_cast<double>( real_observations )
                 - static_cast<double>( num_continuum_pars );
    dof = (std::max)( 1.0, dof );

    total_chi2 += chi2;
    total_DOF += dof;

    const double chi2_dof = chi2 / ((dof > 0.0) ? dof : 1.0);
    for( PeakDef &p : roi_peaks )
      p.set_coefficient(chi2_dof, PeakDef::CoefficientType::Chi2DOF );

    // Hand the peaks back in the exposure of the spectrum they will be drawn over.
    //
    //  The Gaussian had to be projected to the future measurement to build `fixed_signal`, because
    //  that is the block whose continuum is also scaled by the ratio - inside the likelihood the
    //  two match.  The continuum written onto these peaks, though, is the one solved for the
    //  *reference* block, and `#measurement` is the reference spectrum.  Returning a projected
    //  Gaussian on it made the drawn peak stand too tall over both its own continuum and the data,
    //  by exactly the exposure ratio, and made `amplitude()` mean something different from every
    //  other field on the same peak.
    //
    //  So the amplitude is brought back to reference exposure here: what this activity would have
    //  looked like in the measurement actually in hand.  Consumers wanting counts for the predicted
    //  measurement multiply by `DeconComputeResults::exposure_ratio`.
    if( exposure_ratio != 1.0 )
    {
      for( PeakDef &p : roi_peaks )
        p.setAmplitude( p.amplitude() / exposure_ratio );
    }//if( exposure_ratio != 1.0 )

    fittedPeaks.insert( end(fittedPeaks), begin(roi_peaks), end(roi_peaks) );
  }//for( const DeconRoiInfo &roi : input.roi_info )

  result.chi2 = total_chi2;
  result.num_degree_of_freedom = static_cast<int>( std::round(total_DOF) );
  result.exposure_ratio = overall_exposure_ratio;
  result.fit_peaks = fittedPeaks;
  
  return result;
}//DeconComputeResults decon_compute_peaks( const DeconComputeInput &input )

  
DeconActivityOrDistanceLimitResult::DeconActivityOrDistanceLimitResult()
    : isDistanceLimit( false ),
    confidenceLevel( 0.0 ),
    minSearchValue( 0.0 ),
    maxSearchValue( 0.0 ),
    baseInput{},
    limitText(),
    quantityLimitStr(),
    bestCh2Text(),
    errorMessage(),
    overallBestChi2( 0.0 ),
    overallBestQuantity( 0.0 ),
    overallBestResults( nullptr ),
    foundUpperCl( false ),
    upperLimit( 0.0 ),
    upperLimitChi2( -1.0 ),
    upperLimitResults( nullptr ),
    foundLowerCl( false ),
    lowerLimit( 0.0 ),
    lowerLimitChi2( -1.0 ),
    lowerLimitResults( nullptr ),
    foundUpperDisplay( false ),
    upperDisplayRange( 0.0 ),
    foundLowerDisplay( false ),
    lowerDisplayRange( 0.0 ),
  chi2s{}
{
}

double decon_limit_delta( const double confidence_level, const DeconLimitType limit_type )
{
  if( !(confidence_level > 0.5) || !(confidence_level < 1.0) )
    throw runtime_error( "decon_limit_delta: confidence level must be between 0.5 and 1." );

  // Activity or distance is the one parameter of interest; `decon_compute_peaks(...)` re-fits the
  //  nuisance continuum at every trial value, so this is a one-parameter profile and a
  //  chi-square(1) threshold is the appropriate asymptotic approximation regardless of how many
  //  continuum coefficients there are.
  const boost::math::chi_squared chi_squared_dist( 1.0 );

  switch( limit_type )
  {
    case DeconLimitType::OneSidedUpperLimit:
      // A chi-square(1) variate is the square of a two-tailed standard normal, so its CDF has
      //  already folded both normal tails together.  A one-sided normal confidence level CL
      //  therefore corresponds to chi-square probability 2*CL - 1; using CL directly here would
      //  count the two-sidedness twice.
      return boost::math::quantile( chi_squared_dist, 2.0*confidence_level - 1.0 );

    case DeconLimitType::CentralInterval:
      // Both tails are wanted, which is exactly what the chi-square(1) CDF already describes.
      return boost::math::quantile( chi_squared_dist, confidence_level );
  }//switch( limit_type )

  throw runtime_error( "decon_limit_delta: unknown limit type." );
}//decon_limit_delta(...)


// The combo boxes are built in these orders, and the URL tokens decode through these, so a
//  reordering of any of the three enums has to be a deliberate act rather than a silent
//  re-selection of a different statistical treatment.
static_assert( static_cast<int>(DeconContinuumNorm::Floating) == 0, "DeconContinuumNorm reordered" );
static_assert( static_cast<int>(DeconContinuumNorm::FixedByEdges) == 1, "DeconContinuumNorm reordered" );
static_assert( static_cast<int>(DeconContinuumNorm::FixedByFullRange) == 2, "DeconContinuumNorm reordered" );
static_assert( static_cast<int>(DeconMeasurementModel::CurrentSpectrum) == 0, "DeconMeasurementModel reordered" );
static_assert( static_cast<int>(DeconMeasurementModel::BackgroundReference) == 1, "DeconMeasurementModel reordered" );
static_assert( static_cast<int>(DeconLimitType::OneSidedUpperLimit) == 0, "DeconLimitType reordered" );
static_assert( static_cast<int>(DeconLimitType::CentralInterval) == 1, "DeconLimitType reordered" );


int num_selectable_continuum_norms()
{
  return 2;  // Floating and FixedByEdges; FixedByFullRange is deprecated and not offered.
}


DeconContinuumNorm continuum_norm_from_index( const int index )
{
  switch( index )
  {
    case 0: return DeconContinuumNorm::Floating;
    case 1: return DeconContinuumNorm::FixedByEdges;
  }

  throw runtime_error( "continuum_norm_from_index: index " + std::to_string(index)
                       + " is not a selectable continuum treatment." );
}


int index_from_continuum_norm( const DeconContinuumNorm norm )
{
  switch( norm )
  {
    case DeconContinuumNorm::Floating:     return 0;
    case DeconContinuumNorm::FixedByEdges: return 1;

    case DeconContinuumNorm::FixedByFullRange:
      // Reaching here means a deprecated stored value was handed straight to a widget instead of
      //  being migrated.  Displaying it as some other treatment is the silent reinterpretation the
      //  migration exists to prevent.
      throw runtime_error( "index_from_continuum_norm: FixedByFullRange is deprecated and not"
                           " selectable; migrate it to BackgroundReference + Floating first." );
  }

  throw runtime_error( "index_from_continuum_norm: unknown continuum treatment." );
}


DeconMeasurementModel measurement_model_from_index( const int index )
{
  switch( index )
  {
    case 0: return DeconMeasurementModel::CurrentSpectrum;
    case 1: return DeconMeasurementModel::BackgroundReference;
  }

  throw runtime_error( "measurement_model_from_index: invalid index "+ std::to_string(index) );
}


int index_from_measurement_model( const DeconMeasurementModel model )
{
  switch( model )
  {
    case DeconMeasurementModel::CurrentSpectrum:     return 0;
    case DeconMeasurementModel::BackgroundReference: return 1;
  }

  throw runtime_error( "index_from_measurement_model: unknown measurement model." );
}


DeconLimitType limit_type_from_index( const int index )
{
  switch( index )
  {
    case 0: return DeconLimitType::OneSidedUpperLimit;
    case 1: return DeconLimitType::CentralInterval;
  }

  throw runtime_error( "limit_type_from_index: invalid index " + std::to_string(index) );
}


int index_from_limit_type( const DeconLimitType type )
{
  switch( type )
  {
    case DeconLimitType::OneSidedUpperLimit: return 0;
    case DeconLimitType::CentralInterval:    return 1;
  }

  throw runtime_error( "index_from_limit_type: unknown limit type." );
}


bool decode_continuum_norm_token( const string &token,
                                  DeconContinuumNorm &norm,
                                  DeconMeasurementModel &model,
                                  bool &migrated )
{
  // Prefix matching, because that is what the existing decoders did and abbreviated tokens are in
  //  the wild.  Note "FIX" and "FULL" are distinguished before any "F" prefix is tried.
  const auto starts_with = [&token]( const char * const prefix ) -> bool {
    return SpecUtils::istarts_with( token, prefix );
  };

  if( starts_with("UNK") || starts_with("FLOAT") )
  {
    norm = DeconContinuumNorm::Floating;
    model = DeconMeasurementModel::CurrentSpectrum;
    migrated = false;
    return true;
  }

  if( starts_with("EDGES") || starts_with("FIX") )
  {
    norm = DeconContinuumNorm::FixedByEdges;
    model = DeconMeasurementModel::CurrentSpectrum;
    migrated = false;
    return true;
  }

  if( starts_with("NOS") || starts_with("FULL") )
  {
    // The retired "assume no signal is present" option.  It was never a continuum treatment - it
    //  asserted the spectrum is signal-free and predicted a future measurement from it - so it maps
    //  onto the measurement model that says exactly that, with an ordinary floating continuum.
    norm = DeconContinuumNorm::Floating;
    model = DeconMeasurementModel::BackgroundReference;
    migrated = true;
    return true;
  }

  return false;
}//decode_continuum_norm_token(...)


std::string continuum_norm_token( const DeconContinuumNorm norm )
{
  switch( norm )
  {
    case DeconContinuumNorm::Floating:     return "FLOAT";
    case DeconContinuumNorm::FixedByEdges: return "EDGES";

    case DeconContinuumNorm::FixedByFullRange:
      throw runtime_error( "continuum_norm_token: FixedByFullRange is deprecated and must be"
                           " migrated rather than written back out." );
  }

  throw runtime_error( "continuum_norm_token: unknown continuum treatment." );
}//continuum_norm_token(...)


DeconLimitTextKind decon_limit_text_kind( const bool is_dist_limit,
                                          const bool found_lower_cl,
                                          const bool found_upper_cl,
                                          const bool is_predicted_sensitivity,
                                          const DeconLimitType limit_type )
{
  // Both crossings found.  Whether that is an interval depends entirely on which threshold was used
  //  to find them: two roots of a one-sided threshold are not a central interval.
  const DeconLimitTextKind both = (limit_type == DeconLimitType::CentralInterval)
                                    ? DeconLimitTextKind::CentralInterval
                                    : DeconLimitTextKind::TwoOneSidedBounds;

  if( is_dist_limit )
  {
    // A distance scan reports from its *lower* crossing - "you could be at least this far away and
    //  still see it".  An upper crossing on its own bounds nothing that can be stated.
    if( found_lower_cl )
      return found_upper_cl ? both : DeconLimitTextKind::DistanceLowerBound;
    return DeconLimitTextKind::None;
  }

  if( !found_upper_cl )
    return DeconLimitTextKind::None;

  // Checked BEFORE the two-crossing case: a predicted sensitivity describes a future measurement,
  //  and must never be worded as a bound on the loaded spectrum - see DeconMeasurementModel.
  //  `foundLowerCl` is forced false upstream so two crossings should be unreachable here, but that
  //  is enforced by an assert(), which is compiled out in Release.
  if( is_predicted_sensitivity )
    return DeconLimitTextKind::PredictedSensitivity;

  if( found_lower_cl )
    return both;

  // A central-interval scan that found only the upper crossing hit the physical boundary at zero.
  //  It is NOT a one-sided bound at this confidence: the root sits at `quantile(chi2(1), CL)`, so
  //  calling it one-sided at CL understates its coverage.  Report it as the interval [0, U].
  if( limit_type == DeconLimitType::CentralInterval )
    return DeconLimitTextKind::CentralIntervalUpperBound;

  return DeconLimitTextKind::OneSidedUpper;
}//decon_limit_text_kind(...)


namespace
{
/** Whether to print the profile-scan trace, and the lock that makes printing it safe.

 These writes used to be unconditional.  `decon_projected_limit` calls this scan from several
 threads at once, so they raced on whatever streambuf the caller had installed - a
 `std::ostringstream` in the tests, which is not thread-safe - and that reproduced as an
 intermittent SIGBUS in about one run in twenty-five.

 Off by default, for two reasons: a few hundred scans' worth of trace is noise during a Monte
 Carlo, and a library function should not write to the console uninvited.  The logs are still
 *built* either way - they cost a couple of `ostringstream`s and are exactly what you want when a
 scan misbehaves - so turning this on gives the same output as before.  The mutex keeps it safe if
 you do.
 */
bool log_decon_scan()
{
  // Function-local static: initialized exactly once, thread-safely, on first use.
  static const bool enabled = [](){
    const char * const value = getenv( "INTERSPEC_LOG_DECON_SCAN" );
    return value && (string(value) == "1");
  }();

  return enabled;
}//log_decon_scan()


std::mutex &decon_scan_log_mutex()
{
  static std::mutex mutex;
  return mutex;
}//decon_scan_log_mutex()

}//anonymous namespace


DeconActivityOrDistanceLimitResult get_activity_or_distance_limits( const double wantedCl,
                      const shared_ptr<const DetectionLimitCalc::DeconComputeInput> base_input,
                      const bool is_dist_limit,
                      const double min_search_quantity,
                      const double max_search_quantity,
                      const bool useCurie,
                      const DeconLimitType limit_type,
                      const DeconLimitDetail detail )
{
  const bool limit_only = (detail == DeconLimitDetail::LimitOnly);

  if( !base_input )
    throw runtime_error( "get_activity_or_distance_limits: invalid base input." );

  if( !(wantedCl > 0.5f) || !(wantedCl < 1.0f) )
    throw runtime_error( "get_activity_or_distance_limits: confidence level must be between 0.5 and 1." );

  DeconActivityOrDistanceLimitResult result;
  result.isDistanceLimit = is_dist_limit;

  result.confidenceLevel = wantedCl;
  result.limitType = limit_type;
  result.minSearchValue = min_search_quantity;
  result.maxSearchValue = max_search_quantity;
  result.baseInput = *base_input;

  const double cl_chi2_delta = decon_limit_delta( static_cast<double>(wantedCl), limit_type );
  result.limitDelta = cl_chi2_delta;

  // A background-reference scan answers "what could a future measurement detect?", not "how much
  //  signal could be in this spectrum?".  Recording it on the result means no display or report
  //  has to infer it, and it is what forbids a lower limit below.
  //
  // Unless the sample has actually been measured: `observed_sample` makes the two blocks a genuine
  //  joint likelihood over observed data, so the answer is an ordinary bound on that sample - which
  //  may legitimately have a lower limit - and must not be worded as a prediction.
  const bool predicted_sensitivity
        = (base_input->measurement_model == DeconMeasurementModel::BackgroundReference)
          && !base_input->observed_sample;
  result.is_predicted_sensitivity = predicted_sensitivity;
  if( (base_input->measurement_model == DeconMeasurementModel::BackgroundReference)
     && base_input->measurement )
  {
    // The sample's own times when it exists, otherwise the reference's dead-time fraction applied
    //  to the exposure being asked about.
    const shared_ptr<const SpecUtils::Measurement> times_from
        = base_input->observed_sample ? base_input->observed_sample : base_input->measurement;

    const double ref_real = (times_from->real_time() > 0.0f)
                              ? static_cast<double>(times_from->real_time())
                              : static_cast<double>(times_from->live_time());
    const double ref_live = (times_from->live_time() > 0.0f)
                              ? static_cast<double>(times_from->live_time()) : ref_real;

    if( base_input->observed_sample )
      result.sampleExposure = ref_live;
    else
      result.sampleExposure = (base_input->sample_exposure > 0.0)
                              ? base_input->sample_exposure : ref_live;

    // Reported in REAL time, which is what the user entered.  A detector reporting only one of the
    //  two times is taken to have zero dead time, so the conversion is the identity there.
    result.sampleRealTime = (ref_live > 0.0) ? (result.sampleExposure * ref_real / ref_live)
                                             : result.sampleExposure;
  }

  // Anything `decon_compute_peaks(...)` had to change about the input applies to every trial value,
  //  so pick it up once here rather than repeating it for each point of the scan.
  {
    DeconComputeInput probe = *base_input;
    probe.activity = 0.0;
    result.warnings = decon_compute_peaks( probe ).warnings;
  }

  // How far above the best-fit statistic the chart is drawn.  Must clear the threshold itself, or
  //  the 4- and 5-sigma choices draw a chart that does not include the level being reported.
  const double yrange = (std::max)( 15.0, cl_chi2_delta + 5.0 );

  const size_t nchi2 = 25;  //approx num chi2 to compute
  static_assert( nchi2 > 2, "We need at least two chi2" );
  
  vector<pair<double,double>> chi2s;
  double overallBestChi2 = 0.0, overallBestQuantity = 0.0;
  double upperLimit = 0.0, lowerLimit = 0.0;
  double upperLimitChi2 = -1.0, lowerLimitChi2 = -1.0;
  bool foundUpperCl = false, foundUpperDisplay = false, foundLowerCl = false, foundLowerDisplay = false;
  
  //
  double quantityRangeMin = 0.0, quantityRangeMax = 0.0;
  
  
  /// \TODO: currently all this stuff assumes a smooth continuously increasing Chi2 with increasing
  ///        activity, but this doesnt have to be the case, especially with quadratic continuums.
  
  std::atomic<size_t> num_iterations( 0 );
  
  //The boost::math::tools::bisect(...) function will make calls using the same value of activity,
  //  so we will cache those values to save some time.
  map<double,double> chi2cache;
  std::mutex chi2cache_mutex;
  
  if( !base_input )
    throw runtime_error( "missing input quantity." );
  
  auto compute_chi2 = [is_dist_limit,&base_input]( const double quantity, int *numDOF = nullptr ) -> double {
    assert( base_input );
    DetectionLimitCalc::DeconComputeInput input = *base_input;
    if( is_dist_limit )
      input.distance = quantity;
    else
      input.activity = quantity;
    const DetectionLimitCalc::DeconComputeResults results
                                                = DetectionLimitCalc::decon_compute_peaks( input );
    
    if( (results.num_degree_of_freedom == 0) && (results.chi2 == 0.0) )
      throw runtime_error( "No DOF" );
    
    if( numDOF )
      *numDOF = results.num_degree_of_freedom;
    
    return results.chi2;
  };//compute_chi2
  
  
  // This next lambda takes either distance or activity for its argument, depending which
  //  limit is being computed
  auto chi2ForQuantity = [&num_iterations,&chi2cache,&chi2cache_mutex,compute_chi2]( double const &quantity ) -> double {
    
    {//begin lock on achi2cache_mutex
      std::lock_guard<std::mutex> lock( chi2cache_mutex );
      const auto pos = chi2cache.find(quantity);
      if( pos != end(chi2cache) )
        return pos->second;
    }//end lock on chi2cache_mutex
    
    if( quantity < 0.0 )
      return std::numeric_limits<double>::max();
    
    const double chi2 = compute_chi2( quantity );
    
    ++num_iterations;
    
    {//begin lock on chi2cache_mutex
      std::lock_guard<std::mutex> lock( chi2cache_mutex );
      chi2cache.insert( std::pair<double,double>{quantity,chi2} );
    }//end lock on chi2cache_mutex
    
    return chi2;
  };//chi2ForQuantity
  
  
  // Locating the minimum with Brent over the whole range assumes the profile is unimodal.  It is
  //  not always: a quadratic continuum, or an unmodelled interference, can put a second basin in
  //  it, and Brent will happily return whichever one it happens to fall into.  So sample an
  //  ordered grid first, take the global sampled minimum, and only then refine locally.  The grid
  //  is spaced quadratically because the minimum is usually near the low end of the range, and it
  //  is reused for the displayed profile so that what is drawn and what was searched agree.
  double search_max = max_search_quantity;

  const size_t num_prescan = 33;
  vector<pair<double,double>> prescan;

  const auto evaluate_prescan = [&prescan,&chi2ForQuantity]( const double low, const double high ){
    prescan.clear();
    for( size_t i = 0; i < num_prescan; ++i )
    {
      const double t = static_cast<double>(i) / static_cast<double>(num_prescan - 1);
      const double quantity = low + (high - low)*t*t;
      prescan.emplace_back( quantity, chi2ForQuantity(quantity) );
    }
  };

  const auto best_prescan_index = [&prescan]() -> size_t {
    size_t best = 0;
    for( size_t i = 1; i < prescan.size(); ++i )
    {
      if( prescan[i].second < prescan[best].second )
        best = i;
    }
    return best;
  };

  evaluate_prescan( min_search_quantity, search_max );

  // If the whole range sits below the threshold there is no crossing to find in it.  Rather than
  //  depend entirely on the caller's range - which for activity is a Currie estimate multiplied by
  //  an arbitrary 50 - grow it geometrically up to a documented cap.
  if( !is_dist_limit )
  {
    const double range_cap = (std::max)( 1.0, max_search_quantity ) * 1.0E6;
    for( size_t expansion = 0; expansion < 20; ++expansion )
    {
      const size_t best = best_prescan_index();
      if( (prescan.back().second - prescan[best].second) >= cl_chi2_delta )
        break;
      if( search_max >= range_cap )
        break;
      search_max = (std::min)( range_cap, 4.0*(std::max)( search_max, 1.0E-12 ) );
      evaluate_prescan( min_search_quantity, search_max );
    }

    if( search_max > max_search_quantity )
    {
      // The caller's maximum is a heuristic (the GUI and batch both use a Currie estimate times
      //  50), so growing past it is intended - but the answer then comes from outside the range
      //  that was asked for, and that has to be visible rather than implied.
      result.maxSearchValue = search_max;
      result.warnings.push_back( "The requested search range did not contain the confidence"
                                " crossing, so it was extended beyond the value asked for.  A"
                                " limit far above the gross-counts estimate usually means the"
                                " gamma line carries little information in this spectrum." );
    }
  }//if( !is_dist_limit )

  {// Refine the global sampled minimum inside its own bracket.
    const size_t best = best_prescan_index();
    const double bracket_low = prescan[(best > 0) ? (best - 1) : 0].first;
    const double bracket_high = prescan[(std::min)(best + 1, prescan.size() - 1)].first;

    overallBestQuantity = prescan[best].first;
    overallBestChi2 = prescan[best].second;

    if( bracket_high > bracket_low )
    {
      //Float has 24 bits of mantissa; 12 bits gets us three significant figures.
      const int bits = 12;
      boost::uintmax_t refine_iter = 100;
      const pair<double, double> refined = boost::math::tools::brent_find_minima(
                                        chi2ForQuantity, bracket_low, bracket_high, bits, refine_iter );
      if( refined.second < overallBestChi2 )
      {
        overallBestQuantity = refined.first;
        overallBestChi2 = refined.second;
      }
    }//if( bracket_high > bracket_low )
  }

  // Note any additional threshold crossings inside the sampled range; the roots found below are
  //  the nearest ones outward from the minimum, which only describes the whole confidence set when
  //  the profile is well behaved.
  {
    size_t num_crossings = 0;
    for( size_t i = 1; i < prescan.size(); ++i )
    {
      const bool below_before = ((prescan[i-1].second - overallBestChi2) < cl_chi2_delta);
      const bool below_after = ((prescan[i].second - overallBestChi2) < cl_chi2_delta);
      num_crossings += (below_before != below_after);
    }

    if( num_crossings > 2 )
      result.warnings.push_back( "The profile crosses the confidence threshold more than twice, so"
                                " the reported limits bound only the interval containing the best"
                                " fit, not the whole confidence set.  This usually means the peak"
                                " model or continuum does not describe the data." );
  }

  const DetectorPeakResponse::EffGeometryType det_geom
  = base_input->drf ? base_input->drf->geometryType()
  : DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
  
  auto print_quantity = [is_dist_limit,det_geom,useCurie]( double quantity, int ndigits = 4 ) -> string {
    if( is_dist_limit )
      return PhysicalUnits::printToBestLengthUnits(quantity,ndigits);
    return PhysicalUnits::printToBestActivityUnits(quantity,ndigits,useCurie)
    + DetectorPeakResponse::det_eff_geom_type_postfix(det_geom);
  };//print_quantity
  
  if( log_decon_scan() )
  {
    std::lock_guard<std::mutex> lock( decon_scan_log_mutex() );
    cout << "Found min X2=" << overallBestChi2 << " with activity "
    << print_quantity(overallBestQuantity)
    << " and it took " << std::dec << num_iterations.load() << " iterations; searched from "
    << print_quantity(min_search_quantity)
    << " to " << print_quantity(search_max)
    << endl;
  }
  
  //boost::math::tools::bracket_and_solve_root(...)
  auto chi2ForRangeLimit = [&chi2ForQuantity, overallBestChi2, yrange]( double const &quantity ) -> double {
    return chi2ForQuantity(quantity) - overallBestChi2 - yrange;
  };
  
  auto chi2ForCL = [&chi2ForQuantity, overallBestChi2, cl_chi2_delta]( double const &quantity ) -> double {
    return chi2ForQuantity(quantity) - overallBestChi2 - cl_chi2_delta;
  };
  
  //Tolerance is called with two values of quantity (activity or distance, depending which limit
  //  is being found); one with the chi2 below root, and one above
  auto tolerance = [chi2ForCL](double quantity_1, double quantity_2) -> bool{
    const double chi2_1 = chi2ForCL(quantity_1);
    const double chi2_2 = chi2ForCL(quantity_2);
    
    // \TODO: make sure tolerance is being used correctly - when printing info out for every call I'm not sure it is being used right... (but answers seem reasonable, so...)
    //cout << "Tolerance called with quantity_1=" << PhysicalUnits::printToBestActivityUnits(quantity_1,false)  //PhysicalUnits::printToBestLengthUnits(quantity_1)
    //     << ", quantity_2=" << PhysicalUnits::printToBestActivityUnits(quantity_2,4,false)
    //     << " ---> chi2_1=" << chi2_1 << ", chi2_2=" << chi2_2 << endl;
    
    return fabs(chi2_1 - chi2_2) < 0.025;
  };//tolerance(...)
  
  //cout << "chi2ForCL(min_search_quantity)=" << chi2ForCL(min_search_quantity) << endl;
  
  SpecUtilsAsync::ThreadPool pool;
  
  //Before trying to find lower-bounding activity, make sure the best value isnt the lowest
  //  possible value (i.e., zero in this case), and that if we go to the lowest possible value,
  //  that the chi2 will increase at least by cl_chi2_delta
  // The two searches below run concurrently, so neither may write to `cout` directly: formatted
  //  output from two threads into one stream buffer is a data race, and when a caller has
  //  redirected `cout` (as the test harness does to silence these diagnostics) that race is on a
  //  plain `ostringstream` with no locking at all - which shows up as an intermittent heap
  //  corruption rather than as interleaved text.  Each task therefore fills its own buffer, and
  //  the main thread prints them once both have finished.
  stringstream lower_log, upper_log;

  pool.post( [&lowerLimit,&quantityRangeMin,&foundLowerCl,&lowerLimitChi2,&foundLowerDisplay,&num_iterations,&lower_log, //quantities we will modify
               min_search_quantity,overallBestQuantity,overallBestChi2,yrange, //values we can capture by value
               &tolerance,&chi2ForCL,&chi2ForQuantity,&print_quantity,&chi2ForRangeLimit //lambdas we will use
             ](){
    const double min_search_chi2 = chi2ForCL(min_search_quantity);
    if( (fabs(min_search_quantity - overallBestQuantity) > 0.001)
       && (min_search_chi2 > 0.0) )
    {
      pair<double,double> lower_val;
      
      boost::uintmax_t max_iter = 100;  //see note about needing to set before every use
      lower_val = boost::math::tools::bisect( chi2ForCL, min_search_quantity, overallBestQuantity, tolerance, max_iter );
      lowerLimit = 0.5*(lower_val.first + lower_val.second);
      foundLowerCl = true;
      lowerLimitChi2 = chi2ForQuantity(lowerLimit);
      lower_log << "lower_val CL activity="
      << print_quantity(lowerLimit)
      << " with Chi2(" << lowerLimit << ")=" << lowerLimitChi2
      << " (Best Chi2(" << overallBestQuantity << ")=" << overallBestChi2
      << "), num_iterations=" << std::dec << num_iterations.load() << " and search range from "
      << print_quantity(min_search_quantity)
      << " to "
      << print_quantity(overallBestQuantity)
      << endl;
      
      const double minActChi2 = chi2ForRangeLimit(min_search_quantity);
      if( minActChi2 < 0.0 )
      {
        quantityRangeMin = min_search_quantity;
        lower_log << "lower_val display activity being set to min_search_quantity ("
        << min_search_quantity << "): minActChi2=" << minActChi2
        << ", with Chi2(" << quantityRangeMin << ")=" << chi2ForQuantity(quantityRangeMin) << endl;
      }else
      {
        try
        {
          boost::uintmax_t max_iter = 100;
          lower_val = boost::math::tools::bisect( chi2ForRangeLimit, min_search_quantity, lowerLimit, tolerance, max_iter );
          quantityRangeMin = 0.5*(lower_val.first + lower_val.second);
          foundLowerDisplay = true;
          lower_log << "lower_val display quantity=" << print_quantity(quantityRangeMin)
          << " wih chi2=" << chi2ForQuantity(quantityRangeMin) << ", num_iterations=" << std::dec << num_iterations.load() << endl;
        }catch( std::exception &e )
        {
          const double delta_act = 0.1*(lowerLimit - quantityRangeMin);
          for( quantityRangeMin = lowerLimit; quantityRangeMin > 0; quantityRangeMin -= delta_act )
          {
            const double this_chi2 = chi2ForQuantity(quantityRangeMin);
            if( this_chi2 >= (overallBestChi2 + yrange) )
              break;
          }
          
          lower_log << "Couldnt find lower-limit of display range properly, so scanned down and found "
          << print_quantity(quantityRangeMin)
          << " where LowerLimit=" << print_quantity(lowerLimit)
          << " and ActRangeMin=" << print_quantity(quantityRangeMin)
          << " and BestActivity" << print_quantity(overallBestQuantity)
          << endl;
        }//try / catch
      }//
    }else
    {
      lowerLimit = 0.0;
      //quantityRangeMin = overallBestQuantity;
      quantityRangeMin = min_search_quantity;
      lower_log << "lower_val activity/distance already at min" << endl;
    }//if( fabs(min_search_quantity - overallBestQuantity) > 0.001*PhysicalUnits::bq ) / else
  } );//pool.post( ... find lower limit ...)
  
  pool.post( [&upperLimit,&quantityRangeMax,&foundUpperCl,&upperLimitChi2,&foundUpperDisplay,&num_iterations,&upper_log, //quantities we will modify
               search_max,overallBestQuantity,overallBestChi2,yrange,is_dist_limit,min_search_quantity, //values we can capture by value
               &tolerance,&chi2ForCL,&chi2ForQuantity,&print_quantity,&chi2ForRangeLimit //lambdas we will use
             ](){
    const double max_search_chi2 = chi2ForCL(search_max);
    if( (fabs(search_max - overallBestQuantity) > 0.001)
       && (max_search_chi2 > 0.0)  )
    {
      pair<double,double> upper_val;
      boost::uintmax_t max_iter = 100;
      upper_val = boost::math::tools::bisect( chi2ForCL, overallBestQuantity, search_max, tolerance, max_iter );
      upperLimit = 0.5*(upper_val.first + upper_val.second);
      foundUpperCl = true;
      upperLimitChi2 = chi2ForQuantity(upperLimit);
      
      upper_log << "upper_val CL activity=" << print_quantity(upperLimit)
      << " wih chi2=" << upperLimitChi2 << ", num_iterations=" << std::dec << num_iterations.load()
      << " and search range from " << print_quantity(overallBestQuantity)
      << " to "
      << print_quantity(search_max)
      << endl;
      
      const double maxSearchChi2 = chi2ForRangeLimit(search_max);
      if( maxSearchChi2 < 0.0 )
      {
        quantityRangeMax = search_max;
        upper_log << "upper_val display activity being set to search_max (" << search_max << "): maxSearchChi2=" << maxSearchChi2 << endl;
      }else
      {
        try
        {
          max_iter = 100;
          upper_val = boost::math::tools::bisect( chi2ForRangeLimit, upperLimit, search_max, tolerance, max_iter );
          quantityRangeMax = 0.5*(upper_val.first + upper_val.second);
          foundUpperDisplay = true;
          upper_log << "upper_val display quantity=" << print_quantity(quantityRangeMax)
          << " wih chi2=" << chi2ForQuantity(quantityRangeMax) << ", num_iterations="
          << std::dec << num_iterations.load() << endl;
        }catch( std::exception &e )
        {
          const double delta_act = std::max( 0.1*fabs(upperLimit - overallBestQuantity), 0.01*fabs(search_max - upperLimit) );
          for( quantityRangeMax = upperLimit; quantityRangeMax < search_max; quantityRangeMax += delta_act )
          {
            const double this_chi2 = chi2ForQuantity(quantityRangeMax);
            if( this_chi2 >= (overallBestChi2 + yrange) )
              break;
          }
          
          upper_log << "Couldn't find upper-limit of display range properly, so scanned up and found "
          << print_quantity(quantityRangeMax)
          << " where UpperLimit Chi2(" << upperLimit << ")="
          << print_quantity(upperLimit)
          << " and ActRangeMax Chi2(" << quantityRangeMax << ")="
          << print_quantity(quantityRangeMax)
          << " and BestActivity Chi2(" << overallBestQuantity << ")="
          << print_quantity(overallBestQuantity)
          << endl;
        }//try / catch
      }
    }else
    {
      upperLimit = overallBestQuantity;
      quantityRangeMax = search_max;
      
      if( is_dist_limit )
      {
        // We might be at a huge distance, so lets find the distance at which we would start to
        //  kinda see something, ever so slightly
        try
        {
          auto chi2ForMinDelta = [&chi2ForQuantity, overallBestChi2, yrange]( double const &quantity ) -> double {
            return chi2ForQuantity(quantity) - overallBestChi2 - 0.01;
          };
          
          boost::uintmax_t max_iter = 100;
          const auto effective_upper_val = boost::math::tools::bisect( chi2ForMinDelta, min_search_quantity, overallBestQuantity, tolerance, max_iter );
          upperLimit = 0.5*(effective_upper_val.first + effective_upper_val.second);
          quantityRangeMax = upperLimit;
        }catch( std::exception & )
        {
          // A flat distance profile has no useful display crossing.  Preserve the explicit
          // `foundUpperDisplay == false` result rather than aborting a Debug batch calculation.
        }
        //overallBestQuantity
      }//if( is_dist_limit )
      
      upper_log << "upper_val activity already at max" << endl;
    }
  } );//pool.post( ... find upper limit ...)

  // Both searches must finish before their diagnostics can safely reach a stream that anything
  //  else might also be writing to.  Emitted even when a task threw, since a scan that failed is
  //  exactly when the log is worth having.
  try
  {
    pool.join();
  }catch( ... )
  {
    if( log_decon_scan() )
    {
      std::lock_guard<std::mutex> lock( decon_scan_log_mutex() );
      cout << lower_log.str() << upper_log.str();
    }
    throw;
  }

  if( log_decon_scan() )
  {
    std::lock_guard<std::mutex> lock( decon_scan_log_mutex() );
    cout << lower_log.str() << upper_log.str();
    cout << "Found best chi2 and ranges with num_iterations=" << std::dec
         << num_iterations.load() << endl;
  }

  // The projected sample is taken at the counts it expects, with no measurement noise, so by construction it holds
  //  no excess over the null: there is nothing for a lower bound to exclude zero with, and any
  //  crossing found on that side is numerical noise rather than a detection.
  //
  //  This has to happen *before* the limit text and the display ranges are built, not after them.
  //  Those are driven by `foundLowerCl`, so suppressing the bound later would leave a result that
  //  says `foundLowerCl == false` while its own `limitText` reads "Between L and U at 95% CL" -
  //  precisely the two-sided-bound-on-the-current-spectrum reading this whole flag exists to
  //  prevent.  Reported rather than asserted, since a best fit landing a hair above zero can
  //  produce one and aborting a Debug run over a numerical artifact would be worse than saying so.
  if( predicted_sensitivity && foundLowerCl )
  {
    result.warnings.push_back( "A lower confidence bound was found for a predicted-sensitivity"
                              " calculation, where the null expectation contains no excess by"
                              " construction.  It has been discarded as a numerical artifact." );
  }

  if( predicted_sensitivity )
  {
    foundLowerCl = false;
    lowerLimit = 0.0;
    lowerLimitChi2 = -1.0;
    result.lowerLimitResults.reset();
  }
  
  if( quantityRangeMax < quantityRangeMin )
    std::swap( quantityRangeMin, quantityRangeMax );
  
  if( quantityRangeMax == quantityRangeMin )
  {
    quantityRangeMin = 0.9*quantityRangeMin;
    quantityRangeMax = 1.1*quantityRangeMin;
  }
  
  const double initialRangeDelta = fabs(quantityRangeMax - quantityRangeMin);
  if( is_dist_limit && !foundUpperCl )
  {
    // This is when there are nearly zero or negative counts so the Chi2 will just stay flat
    //  at larger and larger distances; in this case we have set quantityRangeMax to be approx
    //  where you start getting a little effect, so now we'll add in a little area after
    //  this point so you can see the Chi2 curve is flattened out
    quantityRangeMax += 0.33 * initialRangeDelta;
  }
  
  if( foundLowerDisplay )
    quantityRangeMin = std::max( min_search_quantity, quantityRangeMin - 0.1*initialRangeDelta );
  
  if( foundUpperDisplay )
    quantityRangeMax = std::min( search_max, quantityRangeMax + 0.1*initialRangeDelta );
  
  // If we didnt find an upper limit, then only display to a fwe multiples of lower limit,
  //  not entire range
  if( is_dist_limit && !foundUpperCl && !foundUpperDisplay && foundLowerCl )
    quantityRangeMax = std::min(quantityRangeMax, 3*lowerLimit ); //3 is arbirary
  
  // Points evenly covering the display range, plus every point of the ordered pre-scan that falls
  //  inside it.  Including the pre-scan means the drawn profile shows the same structure the root
  //  search actually saw - in particular a second basin, or an extra threshold crossing, cannot be
  //  present in the search but missing from the chart.
  if( !limit_only )
  {
    const double quantity_delta = fabs(quantityRangeMax - quantityRangeMin) / nchi2;
    vector<double> quantities;
    for( size_t i = 0; i <= nchi2; ++i )
      quantities.push_back( quantityRangeMin + quantity_delta*i );

    for( const pair<double,double> &point : prescan )
    {
      if( (point.first >= quantityRangeMin) && (point.first <= quantityRangeMax) )
        quantities.push_back( point.first );
    }

    std::sort( begin(quantities), end(quantities) );
    quantities.erase( std::unique( begin(quantities), end(quantities) ), end(quantities) );

    chi2s.resize( quantities.size() );
    for( size_t i = 0; i < quantities.size(); ++i )
    {
      const double quantity = quantities[i];
      pool.post( [i, quantity, &chi2s, &chi2ForQuantity](){
        chi2s[i].first = quantity;
        chi2s[i].second = chi2ForQuantity(quantity);
      } );
    }
    pool.join();
  }
  
  const double distance = is_dist_limit ? lowerLimit : base_input->distance;
  const double activity = is_dist_limit ? base_input->activity : upperLimit;
  const double other_quantity = is_dist_limit ? activity : distance;
  
  const auto localComputeForActivity = [base_input,limit_only]( const double activity,
                                              const double distance,
                                              double &chi2, int &numDOF )
      -> std::shared_ptr<const DetectionLimitCalc::DeconComputeResults> {
    chi2 = 0.0;
    numDOF = 0;

    // Under `DeconLimitDetail::LimitOnly` the returned peak sets are never read, and each of these
    //  is a full `decon_compute_peaks`.  Returning early is most of the saving.
    if( limit_only )
      return nullptr;
    std::vector<PeakDef> peaks;
    
    shared_ptr<DetectionLimitCalc::DeconComputeInput> input = make_shared<DetectionLimitCalc::DeconComputeInput>( *base_input );
    input->activity = activity;
    input->distance = distance;
    
    DetectionLimitCalc::DeconComputeResults results
                  = DetectionLimitCalc::decon_compute_peaks( *input );
    
    peaks = results.fit_peaks;
    chi2 = results.chi2;
    numDOF = results.num_degree_of_freedom;
    
    return make_shared<const DetectionLimitCalc::DeconComputeResults>(results);
  };//void localComputeForActivity(...)
  
  
  const string quantity_str = print_quantity(overallBestQuantity, 3);
  char buffer[256];

  if( !limit_only )
  {
    pool.post( [&result,is_dist_limit,&localComputeForActivity,other_quantity,overallBestQuantity](){
      double dummy_chi2;
      int dummy_numDOF;
      if( is_dist_limit )
        result.overallBestResults = localComputeForActivity( other_quantity, overallBestQuantity, dummy_chi2, dummy_numDOF );
      else
        result.overallBestResults = localComputeForActivity( overallBestQuantity, other_quantity, dummy_chi2, dummy_numDOF );
    } );
  }//if( !limit_only )
  
  // TODO: put all below computations into another thread
  string limit_str;
  bool inconsistent_recomputation = false;
  const auto chi2_recomputation_matches = []( const double first, const double second )
  {
    if( !std::isfinite(first) || !std::isfinite(second) )
      return false;
    const double difference = fabs( first - second );
    return (difference < 0.01)
           || (difference < 0.01*(std::max)(fabs(first), fabs(second)));
  };
  // Which sentence this scan gets to make.  Recorded on the result so the GUI can render the same
  //  statement localized rather than re-deriving the case from the flags.
  const DeconLimitTextKind text_kind = decon_limit_text_kind( is_dist_limit, foundLowerCl,
                                                             foundUpperCl, predicted_sensitivity,
                                                             limit_type );
  result.limitTextKind = text_kind;

  // `is_predicted_sensitivity` forces `foundLowerCl` false above, so a two-crossing result and a
  //  predicted sensitivity cannot both happen.
  assert( !predicted_sensitivity
          || ((text_kind != DeconLimitTextKind::CentralInterval)
              && (text_kind != DeconLimitTextKind::TwoOneSidedBounds)) );

  const string cl_str = confidence_level_pct_str( wantedCl );

  // Only the best-fit recomputation below is guaranteed to run, so the reported DOF is taken from
  //  it; these recomputations write into a scratch that is deliberately discarded.
  int scratch_dof = 0;

  switch( text_kind )
  {
    case DeconLimitTextKind::CentralInterval:
    case DeconLimitTextKind::TwoOneSidedBounds:
    {
      double lowerQuantityChi2 = -999.9, upperQuantityChi2 = -999.9;
      if( is_dist_limit )
      {
        result.lowerLimitResults = localComputeForActivity( other_quantity, lowerLimit, lowerQuantityChi2, scratch_dof );
        result.upperLimitResults = localComputeForActivity( other_quantity, upperLimit, upperQuantityChi2, scratch_dof );
      }else
      {
        result.lowerLimitResults = localComputeForActivity( lowerLimit, other_quantity, lowerQuantityChi2, scratch_dof );
        result.upperLimitResults = localComputeForActivity( upperLimit, other_quantity, upperQuantityChi2, scratch_dof );
      }

      inconsistent_recomputation = !limit_only
                                   && (!chi2_recomputation_matches( lowerQuantityChi2, lowerLimitChi2 )
                                       || !chi2_recomputation_matches( upperQuantityChi2, upperLimitChi2 ));

      limit_str = print_quantity( overallBestQuantity, 3 );
      const string lower_limit_str = print_quantity( lowerLimit, 2 );
      const string upper_limit_str = print_quantity( upperLimit, 2 );

      // Chi2 at upper and lower limits *should* be the same, but since I dont totally trust
      //  everything yet, we'll allow showing a discrepancy so we can see something is up
      if( fabs(lowerQuantityChi2 - upperQuantityChi2) < 0.05*(std::max)(lowerQuantityChi2, upperQuantityChi2) )
        snprintf( buffer, sizeof(buffer), "%.1f", 0.5*(lowerQuantityChi2 + upperQuantityChi2) );
      else
        snprintf( buffer, sizeof(buffer), "%.1f and %.1f", lowerQuantityChi2, upperQuantityChi2 );

      const string chi2_str = buffer;

      if( text_kind == DeconLimitTextKind::CentralInterval )
      {
        snprintf( buffer, sizeof(buffer), "Between %s and %s at %s central CL, &chi;<sup>2</sup>=%s",
                 lower_limit_str.c_str(), upper_limit_str.c_str(), cl_str.c_str(), chi2_str.c_str() );
      }else
      {
        // Two roots of a ONE-sided threshold, which is not a central interval: together they cover
        //  only about 2*CL-1 centrally.  Saying "between L and U at 95% CL" here - as this did
        //  before Increment C - overstates the coverage by about a factor of two in the tail.
        snprintf( buffer, sizeof(buffer),
                 "Two one-sided %s bounds: %s to %s (~%s central coverage), &chi;<sup>2</sup>=%s",
                 cl_str.c_str(), lower_limit_str.c_str(), upper_limit_str.c_str(),
                 confidence_level_pct_str( 2.0*wantedCl - 1.0 ).c_str(), chi2_str.c_str() );
      }

      break;
    }//case CentralInterval / TwoOneSidedBounds

    case DeconLimitTextKind::DistanceLowerBound:
    {
      double lowerQuantityChi2 = -999.9;
      result.lowerLimitResults = localComputeForActivity( other_quantity, lowerLimit, lowerQuantityChi2, scratch_dof );

      inconsistent_recomputation = !limit_only
                                   && !chi2_recomputation_matches( lowerQuantityChi2, lowerLimitChi2 );

      limit_str = print_quantity( lowerLimit, 3 );
      const string print_limit_str = print_quantity( lowerLimit, 2 );

      snprintf( buffer, sizeof(buffer),
               "Distance &ge;%s at %s CL (one-sided lower limit), &chi;<sup>2</sup>=%.1f",
               print_limit_str.c_str(), cl_str.c_str(), lowerQuantityChi2 );
      break;
    }//case DistanceLowerBound

    case DeconLimitTextKind::OneSidedUpper:
    case DeconLimitTextKind::PredictedSensitivity:
    case DeconLimitTextKind::CentralIntervalUpperBound:
    {
      double upperQuantityChi2 = -999.9;
      result.upperLimitResults = localComputeForActivity( upperLimit, other_quantity, upperQuantityChi2, scratch_dof );

      inconsistent_recomputation = !limit_only
                                   && !chi2_recomputation_matches( upperQuantityChi2, upperLimitChi2 );

      limit_str = print_quantity( upperLimit, 3 );
      const string print_limit_str = print_quantity( upperLimit, 2 );

      if( text_kind == DeconLimitTextKind::PredictedSensitivity )
      {
        // Not a bound on the loaded spectrum, and must never be worded as one.
        const string exposure_str = PhysicalUnits::printToBestTimeUnits( result.sampleRealTime, 3 );
        snprintf( buffer, sizeof(buffer),
                 "Predicted sensitivity: less than %s at %s CL for a %s measurement,"
                 " &chi;<sup>2</sup>=%.1f",
                 print_limit_str.c_str(), cl_str.c_str(), exposure_str.c_str(), upperQuantityChi2 );
      }else if( text_kind == DeconLimitTextKind::CentralIntervalUpperBound )
      {
        // Lower endpoint pinned at the physical boundary; see DeconLimitTextKind.
        snprintf( buffer, sizeof(buffer),
                 "Between 0 and %s at %s central CL, &chi;<sup>2</sup>=%.1f",
                 print_limit_str.c_str(), cl_str.c_str(), upperQuantityChi2 );
      }else
      {
        snprintf( buffer, sizeof(buffer),
                 "Less than %s at %s CL (one-sided upper limit), &chi;<sup>2</sup>=%.1f",
                 print_limit_str.c_str(), cl_str.c_str(), upperQuantityChi2 );
      }

      break;
    }//case OneSidedUpper / PredictedSensitivity

    case DeconLimitTextKind::None:
    {
      limit_str = print_quantity( overallBestQuantity, 3 );

      if( !foundUpperCl && !foundLowerCl )
      {
        snprintf( buffer, sizeof(buffer), "Error: failed upper or lower limits at %s", cl_str.c_str() );
      }else if( foundLowerCl )
      {
        // An activity scan that crossed going up but never came back down.
        result.errorMessage = "The activity scan found a lower crossing but did not bracket an upper limit.";
        snprintf( buffer, sizeof(buffer), "Error: Didn't find %s CL activity", cl_str.c_str() );
      }else if( is_dist_limit )
      {
        result.errorMessage = "The distance scan found an upper crossing but did not bracket a lower limit.";
        snprintf( buffer, sizeof(buffer), "Error: Didn't find %s CL det. distance", cl_str.c_str() );
      }else
      {
        result.errorMessage = "The profile scan did not bracket an upper confidence limit within the search range.";
        snprintf( buffer, sizeof(buffer), "Error: Didn't find %s CL upper limit", cl_str.c_str() );
      }

      break;
    }//case None
  }//switch( text_kind )

  if( inconsistent_recomputation )
  {
    foundLowerCl = false;
    foundUpperCl = false;
    result.lowerLimitResults.reset();
    result.upperLimitResults.reset();
    // The statement made above no longer holds, so retract the classification with it - otherwise a
    //  display driven by `limitTextKind` would still render "Less than X at 95% CL" from a result
    //  whose own `limitText` says the crossing was not reproducible.
    result.limitTextKind = DeconLimitTextKind::None;
    result.errorMessage = "The profile scan was not reproducible when its confidence crossing was recomputed.";
    limit_str = print_quantity( overallBestQuantity, 3 );
    snprintf( buffer, sizeof(buffer), "Error: profile scan confidence crossing was inconsistent" );
  }

  pool.join();

  result.limitText = buffer;
  result.quantityLimitStr = limit_str;
  
  result.overallBestChi2 = overallBestChi2;
  result.overallBestQuantity = overallBestQuantity;
  result.upperLimit = upperLimit;
  result.upperLimitChi2 = upperLimitChi2;
  result.lowerLimit = lowerLimit;
  result.lowerLimitChi2 = lowerLimitChi2;
  result.foundUpperCl = foundUpperCl;
  result.foundUpperDisplay = foundUpperDisplay;
  result.upperDisplayRange = quantityRangeMax;
  result.foundLowerCl = foundLowerCl;
  result.foundLowerDisplay = foundLowerDisplay;
  result.lowerDisplayRange = quantityRangeMin;

  if( !foundUpperCl && result.errorMessage.empty() )
    result.errorMessage = "The profile scan did not bracket an upper confidence limit within the search range.";


  // Take the DOF from the best-fit recomputation, which is the one that always runs.  Reading it
  //  from the per-branch recomputations left it at zero in the branches that do not recompute, and
  //  again whenever `inconsistent_recomputation` fired - so this line used to print ", 0 DOF".
  const int numDOF = result.overallBestResults ? result.overallBestResults->num_degree_of_freedom : 0;

  if( is_dist_limit )
  {
    if( foundUpperCl )
    {
      snprintf( buffer, sizeof(buffer), "Best &chi;<sup>2</sup> of %.1f at %s, %i DOF",
               overallBestChi2, quantity_str.c_str(), numDOF );
    }else
    {
      snprintf( buffer, sizeof(buffer), "&chi;<sup>2</sup> is %.1f at large distance, %i DOF",
               overallBestChi2, numDOF );
    }
  }else
  {
    snprintf( buffer, sizeof(buffer), "Best &chi;<sup>2</sup> of %.1f at %s, %i DOF",
             overallBestChi2, quantity_str.c_str(), numDOF );
  }//if( is_dist_limit ) / else
  
  
  result.bestCh2Text = buffer;
  result.chi2s = chi2s;

  return result;
};//get_activity_or_distance_limits(...).


namespace
{
/** Portable samplers for the projection Monte Carlo.

 `std::poisson_distribution` and `std::gamma_distribution` are implementation-defined, so the same
 seed gives different numbers on libstdc++, libc++ and MSVC.  A reported detection limit must not
 depend on which platform computed it, so the algorithms are fixed here.

 The same algorithms appear in `test_DetectionLimit.cpp`'s study helpers.  They are deliberately not
 shared: those are pinned to the seeds of studies already published in the review, and coupling them
 would mean any change here silently moved every one of those numbers.
 */
double projection_uniform( std::mt19937 &generator )
{
  // Straight from the engine's specified 32-bit output; `std::generate_canonical` is free to
  //  consume a different number of words on different standard libraries.
  return ( static_cast<double>( generator() ) + 0.5 ) / 4294967296.0;
}


double projection_normal( std::mt19937 &generator )
{
  const double u1 = projection_uniform( generator );
  const double u2 = projection_uniform( generator );
  return std::sqrt( -2.0*std::log(u1) ) * std::cos( 2.0*boost::math::constants::pi<double>()*u2 );
}


/** Gamma(shape,1) by Marsaglia-Tsang (2000), with their `shape < 1` boost. */
double projection_gamma( const double shape, std::mt19937 &generator )
{
  if( !(shape > 0.0) )
    return 0.0;

  if( shape < 1.0 )
  {
    // Both factors draw from `generator`, and the operands of `*` are unsequenced in C++17, so
    //  written as a single expression the compiler is free to choose which draws first - making the
    //  random stream, and so the reported limit, depend on the compiler.  That is precisely what
    //  these hand-rolled samplers exist to prevent, and this branch is the common one: a channel
    //  that recorded zero counts has shape `0 + alpha` < 1, which is most channels in the low-count
    //  regime a detection limit lives in.  Sequenced explicitly.
    const double boosted = projection_gamma( shape + 1.0, generator );
    const double uniform = projection_uniform( generator );
    return boosted * std::pow( uniform, 1.0/shape );
  }

  const double d = shape - 1.0/3.0;
  const double c = 1.0 / std::sqrt( 9.0*d );

  for( size_t iteration = 0; iteration < 1000; ++iteration )
  {
    double x = 0.0, v = 0.0;
    do
    {
      x = projection_normal( generator );
      v = 1.0 + c*x;
    }while( v <= 0.0 );

    v = v*v*v;
    const double u = projection_uniform( generator );
    const double x2 = x*x;

    if( u < (1.0 - 0.0331*x2*x2) )
      return d*v;
    if( std::log(u) < (0.5*x2 + d*(1.0 - v + std::log(v))) )
      return d*v;
  }

  return d;
}//projection_gamma(...)


/** Knuth's product method below 10, Hoermann's PTRS transformed rejection above it.  Both exact;
 the second because the first is O(mean) and channels can hold thousands of counts.
 */
double projection_poisson( const double mean, std::mt19937 &generator )
{
  if( !(mean > 0.0) )
    return 0.0;

  if( mean < 10.0 )
  {
    const double limit = std::exp( -mean );
    double product = 1.0;
    int count = 0;
    for( ; count < 10000; ++count )
    {
      product *= projection_uniform( generator );
      if( product <= limit )
        break;
    }
    return count;
  }//if( mean < 10.0 )

  const double b = 0.931 + 2.53*std::sqrt( mean );
  const double a = -0.059 + 0.02483*b;
  const double inverse_alpha = 1.1239 + 1.1328/( b - 3.4 );
  const double v_r = 0.9277 - 3.6224/( b - 2.0 );

  for( size_t iteration = 0; iteration < 10000; ++iteration )
  {
    const double u = projection_uniform( generator ) - 0.5;
    const double v = projection_uniform( generator );
    const double us = 0.5 - std::fabs( u );
    const double k = std::floor( (2.0*a/us + b)*u + mean + 0.43 );

    if( (us >= 0.07) && (v <= v_r) )
      return k;

    if( (k < 0.0) || ((us < 0.013) && (v > us)) )
      continue;

    if( std::log( v*inverse_alpha/(a/(us*us) + b) )
       <= (-mean + k*std::log(mean) - std::lgamma(k + 1.0)) )
    {
      return k;
    }
  }//for( bounded rejection loop )

  return mean;
}//projection_poisson(...)


/** A seed derived from the inputs, so the same spectrum and the same settings always give the same
 predicted limit.  A reported MDA that moves when nothing moved is worse than one that is slightly
 wrong, and a Monte Carlo median carries a few percent of sampling error however many trials it
 runs.
 */
uint32_t projection_seed( const shared_ptr<const SpecUtils::Measurement> &reference,
                          const double planned_real_time,
                          const size_t first_channel, const size_t last_channel )
{
  uint64_t hash = 14695981039346656037ull;
  const uint64_t prime = 1099511628211ull;

  const auto mix = [&hash,prime]( const double value ){
    // Bit pattern of the value, so the hash sees the whole mantissa rather than a rounded string.
    uint64_t bits = 0;
    static_assert( sizeof(bits) == sizeof(value), "double is not 64 bits" );
    memcpy( &bits, &value, sizeof(bits) );
    for( size_t byte = 0; byte < sizeof(bits); ++byte )
    {
      hash ^= ( (bits >> (8*byte)) & 0xFF );
      hash *= prime;
    }
  };

  mix( planned_real_time );
  mix( static_cast<double>( reference ? reference->live_time() : 0.0f ) );
  mix( static_cast<double>( reference ? reference->real_time() : 0.0f ) );
  mix( static_cast<double>( first_channel ) );
  mix( static_cast<double>( last_channel ) );

  const shared_ptr<const vector<float>> counts = reference ? reference->gamma_counts() : nullptr;
  if( counts )
  {
    for( size_t i = first_channel; (i <= last_channel) && (i < counts->size()); ++i )
      mix( static_cast<double>( (*counts)[i] ) );
  }

  return 20260810u + static_cast<uint32_t>( hash ^ (hash >> 32) );
}//projection_seed(...)


/** Median and the 16th/84th percentiles of \p values, which is sorted in place. */
void projection_quantiles( vector<double> &values, ProjectedLimit &answer )
{
  std::sort( begin(values), end(values) );

  const auto at = [&values]( const double quantile ) -> double {
    const size_t index = (std::min)( values.size() - 1,
                    static_cast<size_t>( quantile*static_cast<double>(values.size()) ) );
    return values[index];
  };

  answer.median = at( 0.50 );
  answer.lower = at( 0.16 );
  answer.upper = at( 0.84 );
}//projection_quantiles(...)

}//anonymous namespace


double projected_limit_prior_strength()
{
  // Jeffreys.  See `draw_projected_measurement` for what the two neighbouring choices cost.
  return 0.5;
}


shared_ptr<const SpecUtils::Measurement>
draw_projected_measurement( const shared_ptr<const SpecUtils::Measurement> &reference,
                            const float new_real_time,
                            const size_t first_channel,
                            const size_t last_channel,
                            std::mt19937 &generator )
{
  if( !reference )
    throw runtime_error( "draw_projected_measurement: null reference." );

  const float reference_real = (reference->real_time() > 0.0f) ? reference->real_time()
                                                               : reference->live_time();
  if( !(reference_real > 0.0f) )
    throw runtime_error( "draw_projected_measurement: reference has non-positive real time." );
  if( !(new_real_time > 0.0f) )
    throw runtime_error( "draw_projected_measurement: new real time must be > 0." );

  const shared_ptr<const vector<float>> observed = reference->gamma_counts();
  if( !observed || observed->empty() )
    throw runtime_error( "draw_projected_measurement: reference has no gamma counts." );

  if( (last_channel < first_channel) || (last_channel >= observed->size()) )
    throw runtime_error( "draw_projected_measurement: invalid channel window." );

  const shared_ptr<const SpecUtils::EnergyCalibration> cal = reference->energy_calibration();
  if( cal && cal->num_channels() && (observed->size() != cal->num_channels()) )
    throw runtime_error( "draw_projected_measurement: gamma_counts size does not match calibration." );

  const double ratio = static_cast<double>(new_real_time) / static_cast<double>(reference_real);
  const double alpha = projected_limit_prior_strength();

  auto counts = make_shared<vector<float>>( observed->size() );

  // Outside the window the expectation is written straight in: no caller looks there, and drawing
  //  a full 16k-channel spectrum costs ~200x what the limit itself does.
  for( size_t i = 0; i < observed->size(); ++i )
    counts->at(i) = static_cast<float>( ratio * (std::max)( 0.0, static_cast<double>((*observed)[i]) ) );

  for( size_t i = first_channel; i <= last_channel; ++i )
  {
    const double n = (std::max)( 0.0, static_cast<double>( (*observed)[i] ) );
    const double rate = projection_gamma( n + alpha, generator );
    counts->at(i) = static_cast<float>( projection_poisson( ratio*rate, generator ) );
  }

  auto scaled = make_shared<SpecUtils::Measurement>( *reference );

  const float reference_live = (reference->live_time() > 0.0f) ? reference->live_time()
                                                               : reference_real;
  scaled->set_gamma_counts( counts, static_cast<float>( reference_live*ratio ), new_real_time );

  return scaled;
}//draw_projected_measurement(...)


ProjectedLimit currie_projected_limit( const CurrieMdaInput &input,
                                       const double planned_real_time,
                                       const size_t num_trials )
{
  ProjectedLimit answer;

  if( !input.spectrum || !(planned_real_time > 0.0) || (num_trials < 8) )
    return answer;

  // The window the calculation reads: the peak region plus its side channels, and a channel of
  //  slack at each end so rounding cannot leave an undrawn channel inside it.
  size_t first_channel = 0, last_channel = 0;
  try
  {
    const pair<size_t,size_t> region = round_roi_to_channels( input.spectrum,
                                          input.roi_lower_energy, input.roi_upper_energy );
    const size_t below = input.num_lower_side_channels + 1;
    first_channel = (region.first > below) ? (region.first - below) : 0;
    last_channel = region.second + input.num_upper_side_channels + 1;

    const size_t nchannel = input.spectrum->num_gamma_channels();
    if( !nchannel )
      return answer;
    last_channel = (std::min)( last_channel, nchannel - 1 );
  }catch( std::exception & )
  {
    return answer;
  }

  std::mt19937 generator( projection_seed( input.spectrum, planned_real_time,
                                          first_channel, last_channel ) );

  vector<double> limits;
  limits.reserve( num_trials );

  for( size_t trial = 0; trial < num_trials; ++trial )
  {
    ++answer.num_attempted;

    try
    {
      CurrieMdaInput trial_input = input;
      trial_input.spectrum = draw_projected_measurement( input.spectrum,
                                static_cast<float>(planned_real_time),
                                first_channel, last_channel, generator );

      // `detection_limit` (the MDA), NOT `upper_limit`.  These are different quantities with very
      //  different spreads: the MDA depends only on the background level, while `upper_limit`
      //  carries the region's own observed fluctuation and swings several-fold between repeats.
      //  Both callers print this band beside the "minimum reliably detectable activity" line, so
      //  quantiling `upper_limit` attached a factor-of-three spread to a number that moves by a
      //  couple of percent.  The MDA is also the quantity a projection is *for*: planning against
      //  what a future dwell could detect.
      const CurrieMdaResult result = currie_mda_calc( trial_input );
      if( !IsNan(result.detection_limit) && !IsInf(result.detection_limit)
         && (result.detection_limit > 0.0) )
      {
        limits.push_back( result.detection_limit );
      }
    }catch( std::exception & )
    {
      // A realisation that cannot be evaluated is dropped and counted, not substituted for.
    }
  }//for( loop over trials )

  if( limits.size() < (num_trials/2) )
    return answer;

  answer.num_used = limits.size();
  projection_quantiles( limits, answer );
  answer.valid = true;

  return answer;
}//currie_projected_limit(...)


ProjectedLimit decon_projected_limit( const shared_ptr<const DeconComputeInput> &base_input,
                                      const double wantedCl,
                                      const double planned_real_time,
                                      const double max_search_quantity,
                                      const bool useCurie,
                                      const DeconLimitType limit_type,
                                      const ProjectedLimitScoring scoring,
                                      const size_t num_trials,
                                      const size_t num_threads,
                                      const shared_ptr<std::atomic_bool> &cancel )
{
  ProjectedLimit answer;

  if( !base_input || !base_input->measurement || !(planned_real_time > 0.0) || (num_trials < 8) )
    return answer;

  // Scoring each realisation against the reference only means anything if the reference is one.
  if( (scoring == ProjectedLimitScoring::JointWithReference)
     && (base_input->measurement_model != DeconMeasurementModel::BackgroundReference) )
  {
    return answer;
  }

  // Every channel any region reads, plus its side channels.
  const size_t nchannel = base_input->measurement->num_gamma_channels();
  if( !nchannel )
    return answer;

  size_t first_channel = nchannel, last_channel = 0;
  for( const DeconRoiInfo &roi : base_input->roi_info )
  {
    try
    {
      const pair<size_t,size_t> region = round_roi_to_channels( base_input->measurement,
                                                               roi.roi_start, roi.roi_end );
      const size_t below = roi.num_lower_side_channels + 1;
      first_channel = (std::min)( first_channel,
                                 (region.first > below) ? (region.first - below) : size_t(0) );
      last_channel = (std::max)( last_channel,
                        (std::min)( nchannel - 1,
                                   region.second + roi.num_upper_side_channels + 1 ) );
    }catch( std::exception & )
    {
      continue;
    }
  }//for( const DeconRoiInfo &roi : base_input->roi_info )

  if( first_channel > last_channel )
    return answer;

  // Drawn up front, single threaded, so the realisations do not depend on how many threads happen
  //  to be available - the same input must give the same limit on every machine.
  std::mt19937 generator( projection_seed( base_input->measurement, planned_real_time,
                                          first_channel, last_channel ) );

  vector<shared_ptr<const SpecUtils::Measurement>> draws;
  draws.reserve( num_trials );
  for( size_t trial = 0; trial < num_trials; ++trial )
  {
    // Each draw copies a whole Measurement, so on a 16k-channel spectrum this loop allocates and
    //  fills tens of megabytes before any scan starts; a superseded run should not sit through it.
    if( cancel && cancel->load() )
      return answer;

    try
    {
      draws.push_back( draw_projected_measurement( base_input->measurement,
                          static_cast<float>(planned_real_time),
                          first_channel, last_channel, generator ) );
    }catch( std::exception & )
    {
      break;
    }
  }

  if( draws.size() < (num_trials/2) )
    return answer;

  if( cancel && cancel->load() )
    return answer;

  answer.num_attempted = draws.size();

  // Each realisation is an independent profile scan, which is what makes this worth threading:
  //  one scan is milliseconds, and a few hundred of them would otherwise be a visible pause.
  vector<double> results( draws.size(), std::numeric_limits<double>::quiet_NaN() );
  const size_t threads = (std::max)( size_t(1), (std::min)( num_threads, draws.size() ) );

  {
    // A *private* pool, deliberately not `SpecUtilsAsync::ThreadPool`.  On every build except
    //  Apple/Android, `SpecUtils_USE_WT_THREADPOOL` is ON and that class submits into the same
    //  `WServer` io service that serves every request for every session - so a few hundred
    //  realisations would stall the whole application, which is exactly what running this off the
    //  GUI thread was meant to avoid.  Owning our own threads costs one io-service thread (the
    //  caller's) and no more.  The scans run with `DeconLimitDetail::LimitOnly`, so they do not
    //  nest a pool of their own either.
    //
    //  The scans' diagnostic trace is off unless INTERSPEC_LOG_DECON_SCAN is set, and locked when
    //  it is - see `log_decon_scan()` - so nothing here writes to a shared stream.
    boost::asio::thread_pool pool( threads );

    for( size_t thread = 0; thread < threads; ++thread )
    {
      boost::asio::post( pool, [&draws,&results,base_input,thread,threads,wantedCl,
                                max_search_quantity,useCurie,limit_type,scoring,cancel](){
        for( size_t i = thread; i < draws.size(); i += threads )
        {
          // Checked per realisation rather than only at the end: a superseded run would otherwise
          //  hold a thread for its full few hundred scans before anyone looked at the flag.
          if( cancel && cancel->load() )
            return;

          try
          {
            auto trial_input = make_shared<DeconComputeInput>( *base_input );

            if( scoring == ProjectedLimitScoring::JointWithReference )
            {
              // The reference stays the measurement and the realisation joins it as the observed
              //  sample, so this predicts exactly what `BackgroundReference` reports.
              //  `counts_per_bq_into_4pi` is quoted against the reference's live time, which is
              //  still the measurement, so it is left alone - `decon_compute_peaks` applies the
              //  exposure ratio itself.
              trial_input->observed_sample = draws[i];
            }else
            {
              // The realisation *is* the future measurement, so it is scored on its own: no
              //  reference block, no expected-counts step, no exposure ratio.
              trial_input->measurement = draws[i];
              trial_input->measurement_model = DeconMeasurementModel::CurrentSpectrum;
              trial_input->sample_exposure = 0.0;
              trial_input->observed_sample = nullptr;

              for( DeconRoiInfo &roi : trial_input->roi_info )
              {
                for( DeconRoiInfo::PeakInfo &peak : roi.peak_infos )
                {
                  // `counts_per_bq_into_4pi` carries the exposure of the spectrum it came from, so
                  //  it has to move with the realisation or the activity axis changes meaning.
                  if( (base_input->measurement->live_time() > 0.0f) && (draws[i]->live_time() > 0.0f) )
                    peak.counts_per_bq_into_4pi *= ( draws[i]->live_time()
                                                    / base_input->measurement->live_time() );
                }
              }
            }//if( JointWithReference ) / else

            // `LimitOnly`: this reads `upperLimit` and nothing else, so the chi2 display grid and
            //  the three result peak fits - well over half a scan's work - are skipped, and the
            //  scan does not spin up a nested pool.
            const DeconActivityOrDistanceLimitResult limit
                = get_activity_or_distance_limits( wantedCl, trial_input, false, 0.0,
                                                  max_search_quantity, useCurie, limit_type,
                                                  DeconLimitDetail::LimitOnly );
            if( limit.foundUpperCl && !IsNan(limit.upperLimit) && !IsInf(limit.upperLimit) )
              results[i] = limit.upperLimit;
          }catch( std::exception & )
          {
            // Dropped and counted below, not substituted for.
          }
        }//for( this thread's realisations )
      } );
    }//for( loop over threads )

    pool.join();
  }

  if( cancel && cancel->load() )
    return answer;

  vector<double> limits;
  limits.reserve( results.size() );
  for( const double value : results )
  {
    if( !IsNan(value) )
      limits.push_back( value );
  }

  if( limits.size() < (answer.num_attempted/2) )
    return answer;

  answer.num_used = limits.size();
  projection_quantiles( limits, answer );
  answer.valid = true;

  return answer;
}//decon_projected_limit(...)


shared_ptr<const SpecUtils::Measurement>
scale_spectrum_for_dwell( const shared_ptr<const SpecUtils::Measurement> &input,
                          const float new_real_time )
{
  if( !input )
    throw runtime_error( "scale_spectrum_for_dwell: null input." );
  // `!(x > 0)` so NaN is rejected too; `x <= 0` is false for NaN and would let it through.
  //  Some detectors report only one of the two times; treat that as zero dead time rather than
  //  refusing to project.
  const float input_real = (input->real_time() > 0.0f) ? input->real_time() : input->live_time();
  if( !(input_real > 0.0f) )
    throw runtime_error( "scale_spectrum_for_dwell: input has non-positive real time." );
  if( !(new_real_time > 0.0f) )
    throw runtime_error( "scale_spectrum_for_dwell: new real time must be > 0." );

  const shared_ptr<const vector<float>> orig_counts = input->gamma_counts();
  if( !orig_counts )
    throw runtime_error( "scale_spectrum_for_dwell: input has no gamma counts." );

  // Channel count must match the existing energy calibration so that
  // SpecUtils::Measurement::set_gamma_counts() preserves the calibration; otherwise
  // it silently drops it.
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = input->energy_calibration();
  if( cal && cal->num_channels() && (orig_counts->size() != cal->num_channels()) )
    throw runtime_error( "scale_spectrum_for_dwell: input gamma_counts size does not match calibration." );

  const double ratio = static_cast<double>(new_real_time) / static_cast<double>(input_real);

  shared_ptr<SpecUtils::Measurement> scaled = make_shared<SpecUtils::Measurement>( *input );

  shared_ptr<vector<float>> new_counts = make_shared<vector<float>>( orig_counts->size() );
  for( size_t i = 0; i < orig_counts->size(); ++i )
  {
    const double v = static_cast<double>( (*orig_counts)[i] ) * ratio;
    (*new_counts)[i] = static_cast<float>( v );
  }
  assert( new_counts->size() == orig_counts->size() );

  // A zero input live_time will produce a zero scaled live_time; downstream
  // counts-per-bq math (which multiplies by spec->live_time()) will then
  // silently zero out.  Caller is responsible for handling that case.
  const float input_live = (input->live_time() > 0.0f) ? input->live_time() : input_real;
  const float new_live_time = static_cast<float>( input_live * ratio );
  scaled->set_gamma_counts( new_counts, new_live_time, new_real_time );

  return scaled;
}//scale_spectrum_for_dwell(...)


PlannedMeasurement plan_measurement( const shared_ptr<const SpecUtils::Measurement> &reference,
                                     const double planned_real_time_s,
                                     const DeconMeasurementModel model )
{
  PlannedMeasurement answer;
  answer.currie = answer.decon = reference;

  // A non-positive time means the user did not ask for a projection; a reference that cannot be
  //  scaled (no real or live time) cannot support one either.
  // `!(x > 0)` rather than `x <= 0` so a NaN time is rejected rather than sailing through and
  //  poisoning every downstream count.
  if( !reference || !(planned_real_time_s > 0.0) )
    return answer;

  // Users think in real time and the result is reported in real time; the likelihood works in live
  //  time.  Some detectors report only one of the two - treat that as zero dead time.
  const double ref_real = (reference->real_time() > 0.0f) ? reference->real_time()
                                                          : reference->live_time();
  const double ref_live = (reference->live_time() > 0.0f) ? reference->live_time() : ref_real;
  if( !(ref_real > 0.0) || !(ref_live > 0.0) )
    return answer;

  answer.planned_real_time = planned_real_time_s;

  answer.currie = scale_spectrum_for_dwell( reference, static_cast<float>(planned_real_time_s) );

  if( model != DeconMeasurementModel::BackgroundReference )
  {
    // The limit describes one projected measurement, so both calculations see the same spectrum.
    answer.decon = answer.currie;
    return answer;
  }

  // The reference counts must stay real counts at their own exposure, or their counting statistics
  //  cannot propagate - that is the entire point of this measurement model.  The planned time is
  //  carried in `sample_exposure` instead, and `decon_compute_peaks` applies it.
  answer.decon = reference;

  // `sample_exposure` is compared against `Measurement::live_time()`, so it is a LIVE time, while
  //  `planned_real_time_s` is entered as a REAL time.  Preserving the dead-time fraction makes this
  //  exactly `answer.currie->live_time()`, which is what makes the Currie and deconvolution paths
  //  describe the same planned measurement rather than two dwells differing by that fraction.
  answer.sample_exposure = planned_real_time_s * (ref_live / ref_real);
  answer.exposure_ratio = answer.sample_exposure / ref_live;

  assert( fabs( answer.sample_exposure - static_cast<double>(answer.currie->live_time()) )
            <= 1.0E-3*(std::max)( 1.0, answer.sample_exposure ) );

  return answer;
}//plan_measurement(...)


bool projected_band_endpoints( const ProjectedLimit &limit,
                               string &lower_multiple,
                               string &upper_multiple )
{
  if( !limit.valid || !(limit.median > 0.0) )
    return false;

  char lower_buffer[32] = { '\0' }, upper_buffer[32] = { '\0' };
  snprintf( lower_buffer, sizeof(lower_buffer), "%.2f", limit.lower/limit.median );
  snprintf( upper_buffer, sizeof(upper_buffer), "%.2f", limit.upper/limit.median );

  // A band whose two ends print the same tells the user nothing they cannot see from the number
  //  itself, so it is left off rather than shown as "1.00 to 1.00".
  if( strcmp( lower_buffer, upper_buffer ) == 0 )
    return false;

  lower_multiple = lower_buffer;
  upper_multiple = upper_buffer;

  return true;
}//projected_band_endpoints(...)


}//namespace DetectionLimitCalc
