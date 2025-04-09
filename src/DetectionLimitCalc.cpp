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

#include <vector>
#include <iostream>
#include <stdexcept>

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"



#include "SpecUtils/SpecFile.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
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
  
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );
  
  const Material * const stainless = matdb.material("stainless-steel NIST");
  assert( stainless );
  
  const auto shield_transmission = [stainless, shielding_thickness]( const double energy ) -> double {
    const double atten_coef = GammaInteractionCalc::transmition_coefficient_material( stainless, energy, shielding_thickness );
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
      
      // result.detection_limit is: Estimate of the true number of signal counts, where we would
      //  reliably detect a signal above the decision threshold.
      const double det_limit_act = result.detection_limit / counts_per_bq;
      
      // result.detection_limit: This is the number of counts in the peak region, _above_ the
      //  predicted continuum number of counts, at which point we will consider signal to be present.
      const double decision_threshold_act = result.detection_limit / counts_per_bq;
      
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
      const double geom_eff = det->fractionalSolidAngle( det->detectorDiameter(), distance );
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
    peak_cont_sum = result.peak_region_counts_sum;
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
    
    const double lower_cont_density = lower_cont_counts / lower_cont_width;
    const double lower_cont_density_uncert = ((lower_cont_counts <= 0.0) ? 0.0 : (lower_cont_density / sqrt(lower_cont_counts)));
    
    const double upper_cont_density = upper_cont_counts / upper_cont_width;
    const double upper_cont_density_uncert = ((upper_cont_counts <= 0.0) ? 0.0 : (upper_cont_density / sqrt(upper_cont_counts)));
    
    const double peak_cont_density = 0.5*(lower_cont_density + upper_cont_density);
    const double peak_cont_density_uncert = 0.5*sqrt( upper_cont_density_uncert*upper_cont_density_uncert
                                                     + lower_cont_density_uncert*lower_cont_density_uncert );
    const double peak_cont_frac_uncert = ((peak_cont_density > 0.0) ? (peak_cont_density_uncert / peak_cont_density) : 1.0);
    
    const double peak_area_width = spec->gamma_channel_upper(result.last_peak_region_channel)
                                    - spec->gamma_channel_lower(result.first_peak_region_channel);
    
    peak_cont_sum = peak_cont_density * peak_area_width;
    peak_cont_sum_uncert = peak_cont_sum * peak_cont_frac_uncert;
    
    // The equation is centered around the input.gamma_energy with the density of counts at normal
    //  value at that point.  The Slope will be through the midpoints of each continuum.
    // TODO: should do a proper least-squares fit to the continuum; I think this will give us a slightly true-er answer
    
    const double lower_cont_mid_energy = spec->gamma_channel_lower(result.first_lower_continuum_channel) + 0.5*lower_cont_width;
    const double upper_cont_mid_energy = spec->gamma_channel_lower(result.first_upper_continuum_channel) + 0.5*upper_cont_width;
    
    result.continuum_eqn[1] = (upper_cont_density - lower_cont_density) / (upper_cont_mid_energy - lower_cont_mid_energy);
    result.continuum_eqn[0] = lower_cont_density - result.continuum_eqn[1]*(lower_cont_mid_energy - input.gamma_energy);
    
#if( PERFORM_DEVELOPER_CHECKS )
    {// begin sanity check on continuum eqn
      const double peak_start_eq = spec->gamma_channel_lower(result.first_peak_region_channel) - input.gamma_energy;
      const double peak_end_eq = spec->gamma_channel_upper(result.last_peak_region_channel) - input.gamma_energy;
      
      const double peak_cont_eq_integral = result.continuum_eqn[0] * (peak_end_eq - peak_start_eq)
      + result.continuum_eqn[1] * 0.5 * (peak_end_eq*peak_end_eq - peak_start_eq*peak_start_eq);
      const double upper_cont_eq = result.continuum_eqn[0] + (upper_cont_mid_energy - input.gamma_energy)*result.continuum_eqn[1];
      
      // Precision tests, for development - if we go down to a precision of 1E-4, instead of 1E-3,
      //  then these tests fail for NaI systems - I'm not sure if its something actually wrong, or
      //  just really bad numerical accuracy (although its hard to imagine going down to only 4
      //  or so, significant figures)
      const double eq_dens = fabs(peak_cont_eq_integral - peak_cont_sum);
      assert( (eq_dens < 0.1)
             || (eq_dens < 1.0E-3*std::max(peak_cont_eq_integral, peak_cont_sum)) );
      
      const double eq_diff = fabs(peak_cont_eq_integral - peak_cont_sum);
      assert( eq_diff < 0.1 || eq_diff < 1.0E-3*std::max(peak_cont_eq_integral, peak_cont_sum) );
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
  if( k*k*add_uncert >= 1.0 )
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
  result.num_degree_of_freedom = 0;
  
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
  
  // We should be good to go,
  vector<PeakDef> inputPeaks, fittedPeaks;
  
  for( const DeconRoiInfo &roi : input.roi_info )
  {
    const float  &roi_start = roi.roi_start; //This _should_ already be rounded to nearest bin edge; TODO: check that this is rounded
    const float  &roi_end = roi.roi_end; //This _should_ already be rounded to nearest bin edge; TODO: check that this is rounded
    
    const DeconContinuumNorm &cont_norm_method = roi.cont_norm_method;
    PeakContinuum::OffsetType continuum_type = roi.continuum_type;
    switch( cont_norm_method )
    {
      case DeconContinuumNorm::Floating:
      case DeconContinuumNorm::FixedByFullRange:
        break;
        
      case DeconContinuumNorm::FixedByEdges:
        assert( continuum_type == PeakContinuum::OffsetType::Linear );
        continuum_type = PeakContinuum::OffsetType::Linear;
        break;
    }//switch( cont_norm_method )
    
    const size_t &num_lower_side_channels = roi.num_lower_side_channels;
    const size_t &num_upper_side_channels = roi.num_upper_side_channels;
    
    
    shared_ptr<PeakContinuum> peak_continuum;
    std::shared_ptr<const SpecUtils::Measurement> computed_global_cont;
    
    // Find the largest peak in the ROIs, energy to use as the continuum "reference energy"
    float reference_energy = 0.0f;
    for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
      reference_energy = std::max( reference_energy, peak_info.energy );
    
    assert( reference_energy != 0.0f );
    
    for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
    {
      const float &energy = peak_info.energy;
      const float fwhm = (peak_info.fwhm > 0.0f) ? peak_info.fwhm : input.drf->peakResolutionFWHM( energy );
      const float sigma = fwhm / 2.634;
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
      
      const float amplitude = input.activity * counts_4pi * det_eff * air_atten;
    
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
        
        size_t nlowerside = num_lower_side_channels;
        size_t nupperside = num_upper_side_channels;
        if( roi.cont_norm_method != DeconContinuumNorm::FixedByEdges )
        {
          //if no value provided, use 4 channels
          if( !nlowerside )
            nlowerside = 4;
          if( !nupperside )
            nupperside = 4;
          
          // Clamp between 2 and 16 channels
          nupperside = ((nupperside < 2) ? size_t(2) : ((nupperside > 16) ? size_t(16) : nupperside)); //std::clamp(...), C++17
          nlowerside = ((nlowerside < 2) ? size_t(2) : ((nlowerside > 16) ? size_t(16) : nlowerside)); //std::clamp(...), C++17
        }else
        {
          assert( num_lower_side_channels > 0 );
          assert( num_upper_side_channels > 0 );
        }
        
        // First, we'll find a linear continuum as the starting point, and then go through
        //  and modify it how we need
        peak_continuum->calc_linear_continuum_eqn( input.measurement, reference_energy,
                                                  roi_start, roi_end, nlowerside, nupperside );
        
        peak_continuum->setType( continuum_type );
        
        
        if( cont_norm_method == DeconContinuumNorm::FixedByFullRange )
        {
          // We'll set a peaks amplitude for zero
          for( size_t order = 0; order < peak_continuum->parameters().size(); ++order )
            peak_continuum->setPolynomialCoefFitFor( order, true );
          
          const double mean = 0.5*(roi_start + roi_end);
          PeakDef worker_peak( 0.5*(roi_start + roi_end), sigma, 0.0 );
          worker_peak.setContinuum( peak_continuum );
          worker_peak.setFitFor( PeakDef::CoefficientType::Mean, false );
          worker_peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
          worker_peak.setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
          
          // Use peak sigma to determine if its a high or low resolution detector, but I dont think
          //  it matters much since the range is defined for the peak, and that is the only place
          //  `isHPGe` is used, I think.
          const double hpge_fwhm = PeakFitUtils::hpge_fwhm_fcn(mean);
          const double nai_fwhm = PeakFitUtils::nai_fwhm_fcn(mean);
          const bool isHPGe = (fabs(fwhm - hpge_fwhm) < fabs(fwhm - nai_fwhm));
          
          std::vector<PeakDef> fit_peak;
          fitPeaks( {worker_peak}, -1.0, -1.0, input.measurement, fit_peak, {}, false, isHPGe );
          
          assert( fit_peak.size() == 1 );
          if( fit_peak.size() == 1 )
          {
            peak_continuum = fit_peak[0].continuum();
            peak.setContinuum( peak_continuum );
          }else
          {
            string msg = "Error fitting DeconContinuumNorm::FixedByFullRange continuum - failed to"
                        " get a peak out - expected 1, got " + std::to_string(fit_peak.size());
            cerr << msg << endl;
#if( PERFORM_DEVELOPER_CHECKS )
            log_developer_error( __func__, msg.c_str() );
            throw runtime_error( msg );
#endif
          }//if( fit_peak.size() == 1 ) / else
        }//if( cont_norm_method == DeconContinuumNorm::FixedByFullRange )
        
        
        for( size_t order = 0; order < peak_continuum->parameters().size(); ++order )
        {
          switch( cont_norm_method )
          {
            case DeconContinuumNorm::Floating:
              peak_continuum->setPolynomialCoefFitFor( order, true );
              break;
              
            case DeconContinuumNorm::FixedByFullRange:
            case DeconContinuumNorm::FixedByEdges:
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
            assert( 0 );
            break;
            
          case PeakContinuum::External:
            if( !computed_global_cont )
              computed_global_cont = estimateContinuum( input.measurement );
            peak_continuum->setExternalContinuum( computed_global_cont );
            break;
        }//switch( assign continuum )
      }//if( peak_continuum ) / else
      
      inputPeaks.push_back( std::move(peak) );
    }//for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  if( inputPeaks.empty() )
    throw runtime_error( "decon_compute_peaks: No peaks given in ROI(s)" );

  // `isHPGe` is only used to define ROI, I think, so it doesnt matter much
  const bool isHPGe = ([=]() -> bool {
    vector<shared_ptr<const PeakDef>> peakv;
    for( const auto &p : inputPeaks )
      peakv.push_back( make_shared<const PeakDef>(p) );
    const auto m = PeakFitUtils::coarse_resolution_from_peaks(peakv);
    return (m == PeakFitUtils::CoarseResolutionType::High);
  })();
  
  
  ROOT::Minuit2::MnUserParameters inputFitPars;
  PeakFitChi2Fcn::addPeaksToFitter( inputFitPars, inputPeaks, input.measurement, PeakFitChi2Fcn::kFitForPeakParameters, isHPGe );
  
  const int npeaks = static_cast<int>( inputPeaks.size() );
  PeakFitChi2Fcn chi2Fcn( npeaks, input.measurement, nullptr );
  chi2Fcn.useReducedChi2( false );
  
  if( inputFitPars.VariableParameters() == 0 )
  {
    // If we choose to "Fix continuum" we can get here
    result.chi2 = chi2Fcn.chi2( inputFitPars.Params().data() );
    result.num_degree_of_freedom = 0;
    
    for( auto &peak : inputPeaks )
    {
      peak.setFitFor( PeakDef::CoefficientType::Mean, false );
      peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
    }
    
    result.fit_peaks = inputPeaks;
    
    return result;
  }//if( inputFitPars.VariableParameters() == 0 )
  
  assert( inputFitPars.VariableParameters() != 0 );
  
  ROOT::Minuit2::MnUserParameterState inputParamState( inputFitPars );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( chi2Fcn, inputParamState, strategy );
  
  unsigned int maxFcnCall = 5000;
  double tolerance = 2.5;
  tolerance = 0.5;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  const ROOT::Minuit2::MnUserParameters &fitParams = fitter.Parameters();
  //  minimum.IsValid()
  //      ROOT::Minuit2::MinimumState minState = minimum.State();
  //      ROOT::Minuit2::MinimumParameters minParams = minState.Parameters();
  
  //    cerr << endl << endl << "EDM=" << minimum.Edm() << endl;
  //    cerr << "MinValue=" <<  minimum.Fval() << endl << endl;
  
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
  
  if( !minimum.IsValid() )
  {
    //XXX - should we try to re-fit here? Or do something to handle the
    //      faliure in some reasonable way?
    cerr << endl << endl << "status is not valid"
    << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
    << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
    << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
    << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
    << "\n\tHasValidParameters: " << minimum.HasValidParameters()
    << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
    << endl;
    if( minimum.IsAboveMaxEdm() )
      cout << "\t\tEDM=" << minimum.Edm() << endl;
  }//if( !minimum.IsValid() )
  
  
  vector<double> fitpars = fitParams.Params();
  vector<double> fiterrors = fitParams.Errors();
  chi2Fcn.parametersToPeaks( fittedPeaks, &fitpars[0], &fiterrors[0] );
  
  double initialChi2 = chi2Fcn.chi2( &fitpars[0] );
  
  //Lets try to keep whether or not to fit parameters should be the same for
  //  the output peaks as the input peaks.
  //Note that this doesnt account for peaks swapping with each other in the fit
  assert( fittedPeaks.size() == inputPeaks.size() );
  
  //for( size_t i = 0; i < near_peaks.size(); ++i )
  //  fittedPeaks[i].inheritUserSelectedOptions( near_peaks[i], true );
  //for( size_t i = 0; i < fixedpeaks.size(); ++i )
  //  fittedPeaks[i+near_peaks.size()].inheritUserSelectedOptions( fixedpeaks[i], true );
  
  const double totalNDF = set_chi2_dof( input.measurement, fittedPeaks, 0, fittedPeaks.size() );
  
  result.chi2 = initialChi2;
  result.num_degree_of_freedom = static_cast<int>( std::round(totalNDF) );
  
  for( auto &peak : fittedPeaks )
  {
    peak.setFitFor( PeakDef::CoefficientType::Mean, false );
    peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
  }
  
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

DeconActivityOrDistanceLimitResult get_activity_or_distance_limits( const float wantedCl,
                      const shared_ptr<const DetectionLimitCalc::DeconComputeInput> base_input,
                      const bool is_dist_limit,
                      const double min_search_quantity,
                      const double max_search_quantity,
                      const bool useCurie )
{
  assert( base_input );
  if( !base_input )
    throw runtime_error( "get_activity_or_distance_limits: invalid base input." );
  
  DeconActivityOrDistanceLimitResult result;
  result.isDistanceLimit = is_dist_limit;
  
  result.confidenceLevel = wantedCl;
  result.minSearchValue = min_search_quantity;
  result.maxSearchValue = max_search_quantity;
  result.baseInput = *base_input;
  
  const double yrange = 15;
  
  // TODO: we are scanning activity or distance, which is a single degree of freedom - but does it matter that
  //       we are marginalizing over (i.e., fitting for) the nuisance parameters of the peaks and
  //       stuff?  I dont *think* so.
  const boost::math::chi_squared chi_squared_dist( 1.0 );
  
  // We want interval corresponding to 95%, where the quantile will give us CDF up to that
  //  point, so we actually want the quantile that covers 97.5% of cases.
  const float twoSidedCl = 0.5 + 0.5*wantedCl;
  
  const double cl_chi2_delta = boost::math::quantile( chi_squared_dist, twoSidedCl );
  
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
  
  
  /// `search_range` returns {best-chi2,best-quantity}
  auto search_range = [chi2ForQuantity]( double min_range, double max_range, boost::uintmax_t &max_iter ) -> pair<double, double> {
    // boost::brent_find_minima first evaluates the input at the range midpoint, then endpoints
    //  and so on - we could do a first pre-scan over the range to help make sure we dont miss
    //  a global minimum.
    
    
    //\TODO: if best activity is at min_search_quantity, it takes 50 iterations inside brent_find_minima
    //      to confirm; we could save this time by using just a little bit of intelligence...
    const int bits = 12; //Float has 24 bits of mantissa; should get us accurate to three significant figures
    
    return boost::math::tools::brent_find_minima( chi2ForQuantity, min_range, max_range, bits, max_iter );
  };//search_range lambda
  
  boost::uintmax_t max_iter = 100;  //this variable gets changed each use, so you need to reset it afterwards
  const pair<double, double> r = search_range( min_search_quantity, max_search_quantity, max_iter );
  
  overallBestChi2 = r.second;
  overallBestQuantity = r.first;
  
  const DetectorPeakResponse::EffGeometryType det_geom
  = base_input->drf ? base_input->drf->geometryType()
  : DetectorPeakResponse::EffGeometryType::FarField;
  
  auto print_quantity = [is_dist_limit,det_geom,useCurie]( double quantity, int ndigits = 4 ) -> string {
    if( is_dist_limit )
      return PhysicalUnits::printToBestLengthUnits(quantity,ndigits);
    return PhysicalUnits::printToBestActivityUnits(quantity,ndigits,useCurie)
    + DetectorPeakResponse::det_eff_geom_type_postfix(det_geom);
  };//print_quantity
  
  cout << "Found min X2=" << overallBestChi2 << " with activity "
  << print_quantity(overallBestQuantity)
  << " and it took " << std::dec << num_iterations.load() << " iterations; searched from "
  << print_quantity(min_search_quantity)
  << " to " << print_quantity(max_search_quantity)
  << endl;
  
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
  pool.post( [&lowerLimit,&quantityRangeMin,&foundLowerCl,&lowerLimitChi2,&foundLowerDisplay,&num_iterations, //quantities we will modify
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
      cout << "lower_val CL activity="
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
        cout << "lower_val display activity being set to min_search_quantity ("
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
          cout << "lower_val display quantity=" << print_quantity(quantityRangeMin)
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
          
          cout << "Couldnt find lower-limit of display range properly, so scanned down and found "
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
      cout << "lower_val activity/distance already at min" << endl;
    }//if( fabs(min_search_quantity - overallBestQuantity) > 0.001*PhysicalUnits::bq ) / else
  } );//pool.post( ... find lower limit ...)
  
  pool.post( [&upperLimit,&quantityRangeMax,&foundUpperCl,&upperLimitChi2,&foundUpperDisplay,&num_iterations, //quantities we will modify
               max_search_quantity,overallBestQuantity,overallBestChi2,yrange,is_dist_limit,min_search_quantity, //values we can capture by value
               &tolerance,&chi2ForCL,&chi2ForQuantity,&print_quantity,&chi2ForRangeLimit //lambdas we will use
             ](){
    const double max_search_chi2 = chi2ForCL(max_search_quantity);
    if( (fabs(max_search_quantity - overallBestQuantity) > 0.001)
       && (max_search_chi2 > 0.0)  )
    {
      pair<double,double> upper_val;
      boost::uintmax_t max_iter = 100;
      upper_val = boost::math::tools::bisect( chi2ForCL, overallBestQuantity, max_search_quantity, tolerance, max_iter );
      upperLimit = 0.5*(upper_val.first + upper_val.second);
      foundUpperCl = true;
      upperLimitChi2 = chi2ForQuantity(upperLimit);
      
      cout << "upper_val CL activity=" << print_quantity(upperLimit)
      << " wih chi2=" << upperLimitChi2 << ", num_iterations=" << std::dec << num_iterations.load()
      << " and search range from " << print_quantity(overallBestQuantity)
      << " to "
      << print_quantity(max_search_quantity)
      << endl;
      
      const double maxSearchChi2 = chi2ForRangeLimit(max_search_quantity);
      if( maxSearchChi2 < 0.0 )
      {
        quantityRangeMax = max_search_quantity;
        cout << "upper_val display activity being set to max_search_quantity (" << max_search_quantity << "): maxSearchChi2=" << maxSearchChi2 << endl;
      }else
      {
        try
        {
          max_iter = 100;
          upper_val = boost::math::tools::bisect( chi2ForRangeLimit, upperLimit, max_search_quantity, tolerance, max_iter );
          quantityRangeMax = 0.5*(upper_val.first + upper_val.second);
          foundUpperDisplay = true;
          cout << "upper_val display quantity=" << print_quantity(quantityRangeMax)
          << " wih chi2=" << chi2ForQuantity(quantityRangeMax) << ", num_iterations="
          << std::dec << num_iterations.load() << endl;
        }catch( std::exception &e )
        {
          const double delta_act = std::max( 0.1*fabs(upperLimit - overallBestQuantity), 0.01*fabs(max_search_quantity - upperLimit) );
          for( quantityRangeMax = upperLimit; quantityRangeMax < max_search_quantity; quantityRangeMax -= delta_act )
          {
            const double this_chi2 = chi2ForQuantity(quantityRangeMax);
            if( this_chi2 >= (overallBestChi2 + yrange) )
              break;
          }
          
          cout << "Couldn't find upper-limit of display range properly, so scanned up and found "
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
      quantityRangeMax = max_search_quantity;
      
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
        }catch( std::exception &e )
        {
          assert( 0 );
        }
        //overallBestQuantity
      }//if( is_dist_limit )
      
      cout << "upper_val activity already at max" << endl;
    }
  } );//pool.post( ... find upper limit ...)
  
  pool.join();
  
  cout << "Found best chi2 and ranges with num_iterations=" << std::dec << num_iterations.load() << endl;
  
  assert( quantityRangeMin <= quantityRangeMax );
  if( quantityRangeMax < quantityRangeMin )
    std::swap( quantityRangeMin, quantityRangeMax );
  
  if( quantityRangeMax == quantityRangeMin )
  {
    assert( !foundLowerCl && !foundUpperCl );
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
    quantityRangeMax = std::min( max_search_quantity, quantityRangeMax + 0.1*initialRangeDelta );
  
  // If we didnt find an upper limit, then only display to a fwe multiples of lower limit,
  //  not entire range
  if( is_dist_limit && !foundUpperCl && !foundUpperDisplay && foundLowerCl )
    quantityRangeMax = std::min(quantityRangeMax, 3*lowerLimit ); //3 is arbirary
  
  // TODO: be a little more intelligent in
  const double quantity_delta = fabs(quantityRangeMax - quantityRangeMin) / nchi2;
  chi2s.resize( nchi2 );
  
  for( size_t i = 0; i < nchi2; ++i )
  {
    const double quantity = quantityRangeMin + quantity_delta*i;
    pool.post( [i, quantity, &chi2s, &chi2ForQuantity](){
      chi2s[i].first = quantity;
      chi2s[i].second = chi2ForQuantity(quantity);
    } );
  }
  pool.join();
  
  const double distance = is_dist_limit ? lowerLimit : base_input->distance;
  const double activity = is_dist_limit ? base_input->activity : upperLimit;
  const double other_quantity = is_dist_limit ? activity : distance;
  
  const auto localComputeForActivity = [base_input]( const double activity, const double distance,
                                              double &chi2, int &numDOF )
      -> std::shared_ptr<const DetectionLimitCalc::DeconComputeResults> {
    chi2 = 0.0;
    numDOF = 0;
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
  
  
  int numDOF = 0;
  
  const string quantity_str = print_quantity(overallBestQuantity, 3);
  char buffer[128];
  
  pool.post( [&result,is_dist_limit,&localComputeForActivity,other_quantity,overallBestQuantity](){
    double dummy_chi2;
    int dummy_numDOF;
    if( is_dist_limit )
      result.overallBestResults = localComputeForActivity( other_quantity, overallBestQuantity, dummy_chi2, dummy_numDOF );
    else
      result.overallBestResults = localComputeForActivity( overallBestQuantity, other_quantity, dummy_chi2, dummy_numDOF );
  } );
  
  // TODO: put all below computations into another thread
  string limit_str;
  if( foundLowerCl && foundUpperCl )
  {
    double lowerQuantityChi2 = -999.9, upperQuantityChi2 = -999.9;
    std::vector<PeakDef> peaks;
    if( is_dist_limit )
    {
      result.lowerLimitResults = localComputeForActivity( other_quantity, lowerLimit, lowerQuantityChi2, numDOF );
      result.upperLimitResults = localComputeForActivity( other_quantity, upperLimit, upperQuantityChi2, numDOF );
    }else
    {
      result.lowerLimitResults = localComputeForActivity( lowerLimit, other_quantity, lowerQuantityChi2, numDOF );
      result.upperLimitResults = localComputeForActivity( upperLimit, other_quantity, upperQuantityChi2, numDOF );
    }
    
    assert( lowerQuantityChi2 == lowerLimitChi2 ); // TODO: check logic to make sure this is definitely true, then remove above computation
    assert( upperQuantityChi2 == upperLimitChi2 ); // TODO: check logic to make sure this is definitely true, then remove above computation
    
    limit_str = print_quantity( overallBestQuantity, 3 );
    const string lower_limit_str = print_quantity( lowerLimit, 2 );
    const string upper_limit_str = print_quantity( upperLimit, 2 );
    
    // Chi2 at upper and lower limits *should* be the same, but since I dont totally trust
    //  everything yet, we'll allow showing a discrepancy so we can see something is up
    if( fabs(lowerQuantityChi2 - upperQuantityChi2) < 0.05*std::max(lowerQuantityChi2, upperQuantityChi2) )
      snprintf( buffer, sizeof(buffer), "%.1f", 0.5*(lowerQuantityChi2 + upperQuantityChi2) );
    else
      snprintf( buffer, sizeof(buffer), "%.1f and %.1f", lowerQuantityChi2, upperQuantityChi2 );
    
    const string chi2_str = buffer;
    
    //snprintf( buffer, sizeof(buffer), "%.1f%% coverage in [%s, %s], &chi;<sup>2</sup>=%s",
    //         0.1*std::round(1000.0*wantedCl), lower_limit_str.c_str(), upper_limit_str.c_str(),
    //         chi2_str.c_str() );
    
    snprintf( buffer, sizeof(buffer), "Between %s and %s at %.1f%% CL, &chi;<sup>2</sup>=%s",
             lower_limit_str.c_str(), upper_limit_str.c_str(),
             0.1*std::round(1000.0*wantedCl),
             chi2_str.c_str() );
  }else if( !foundLowerCl && !foundUpperCl )
  {
    limit_str = print_quantity( overallBestQuantity, 3 );
    snprintf( buffer, sizeof(buffer), "Error: failed upper or lower limits at %.1f%%",
             0.1*std::round(1000.0*wantedCl) );
  }else if( foundLowerCl )
  {
    if( is_dist_limit )
    {
      double lowerQuantityChi2 = -999.9;
      result.lowerLimitResults = localComputeForActivity( other_quantity, lowerLimit, lowerQuantityChi2, numDOF );
      
      assert( lowerQuantityChi2 == lowerLimitChi2 ); // TODO: check logic to make sure this is definitely true, then remove above computation
      
      limit_str = print_quantity( lowerLimit, 3 );
      const string print_limit_str = print_quantity( lowerLimit, 2 );
      
      //More stat-nerd-esk language, maybe, if its even correct, but lets print something
      //  easier to interpret, for the commoners, like myself.
      //snprintf( buffer, sizeof(buffer), "%.1f%% coverage at %s with &chi;<sup>2</sup>=%.1f",
      //         0.1*std::round(1000.0*wantedCl), print_limit_str.c_str(), lowerQuantityChi2 );
      
      snprintf( buffer, sizeof(buffer), "Distance ≥%s at %.1f%% CL, &chi;<sup>2</sup>=%.1f",
               print_limit_str.c_str(), 0.1*std::round(1000.0*wantedCl), lowerQuantityChi2 );
    }else
    {
      assert( 0 );
      //double lowerQuantityChi2 = -999.9;
      //result.lowerLimitResults = localComputeForActivity( lowerLimit, other_quantity, peaks, lowerQuantityChi2, numDOF );
      snprintf( buffer, sizeof(buffer), "Error: Didn't find %.1f%% CL activity",
               0.1*std::round(1000.0*wantedCl));
    }
  }else
  {
    assert( foundUpperCl );
    
    if( is_dist_limit )
    {
      assert( 0 );
      snprintf( buffer, sizeof(buffer), "Error: Didn't find %.1f%% CL det. distance",
               0.1*std::round(1000.0*wantedCl) );
    }else
    {
      double upperQuantityChi2 = -999.9;
      result.upperLimitResults = localComputeForActivity( upperLimit, other_quantity, upperQuantityChi2, numDOF );
      
      // TODO: check logic to make sure this is definitely true, then remove above computation
      //  I think these quantities should be really close, but there may be a small amount of rounding
      cout << "upperQuantityChi2=" << upperQuantityChi2 << ", upperLimitChi2=" << upperLimitChi2 << endl;
      assert( (fabs(upperQuantityChi2 - upperLimitChi2) < 0.01)
             || (fabs(upperQuantityChi2 - upperLimitChi2) < 0.01*std::max(upperQuantityChi2,upperLimitChi2)) );
      
      limit_str = print_quantity(upperLimit,3);
      const string print_limit_str = print_quantity( upperLimit, 2 );
      //snprintf( buffer, sizeof(buffer), "%.1f%% coverage at %s with &chi;<sup>2</sup>=%.1f",
      //         0.1*std::round(1000.0*wantedCl), print_limit_str.c_str(), upperQuantityChi2 );
      snprintf( buffer, sizeof(buffer), "Less than %s at %.1f%% CL, &chi;<sup>2</sup>=%.1f",
               print_limit_str.c_str(), 0.1*std::round(1000.0*wantedCl), upperQuantityChi2 );
    }//if( is_dist_limit ) / else
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
  
  
}//namespace DetectionLimitCalc

