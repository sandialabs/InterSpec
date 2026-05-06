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

#include <map>
#include <set>
#include <regex>
#include <memory>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "ceres/jet.h"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/RelActAutoDev.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActAutoReport.h"
#include "InterSpec/BatchRelActAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/RelActCalcAuto_EnergyCal_imp.hpp"

using namespace std;

namespace RelActAutoDev
{

#if( PERFORM_DEVELOPER_CHECKS )
void check_auto_state_xml_serialization()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  assert( matdb );

  RelActCalcAuto::RelActAutoGuiState state;

  state.options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NonLinearFit;
  state.options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_3;
  state.options.spectrum_title = "Test";
  state.options.skew_type = PeakDef::SkewType::CrystalBall;
  state.options.additional_br_uncert = 0.011;

  RelActCalcAuto::RelEffCurveInput curve;
  curve.nucs_of_el_same_age = false;
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  curve.rel_eff_eqn_order = 0;

  // Fill out self-attenuation
  auto self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>();
  self_atten->atomic_number = 28.5;
  self_atten->material = nullptr;
  self_atten->areal_density = 1.0*PhysicalUnits::g_per_cm2;
  self_atten->fit_atomic_number = false;
  self_atten->lower_fit_atomic_number = 1.0;
  self_atten->upper_fit_atomic_number = 98.0;
  self_atten->fit_areal_density = true;
  self_atten->lower_fit_areal_density = 0.0;
  self_atten->upper_fit_areal_density = 500*PhysicalUnits::g_per_cm2;
  curve.phys_model_self_atten = self_atten;

  // Fill out external attenuation
  auto ext_atten = make_shared<RelActCalc::PhysicalModelShieldInput>();
  ext_atten->atomic_number = 0;
  const shared_ptr<const Material> pu = matdb->material("Pu (plutonium)");
  assert( pu );
  ext_atten->material = pu;
  ext_atten->areal_density = 1.0*PhysicalUnits::g_per_cm2;
  ext_atten->fit_atomic_number = true;
  ext_atten->lower_fit_atomic_number = 2.0;
  ext_atten->upper_fit_atomic_number = 97.0;
  ext_atten->fit_areal_density = false;
  ext_atten->lower_fit_areal_density = 1.0;
  ext_atten->upper_fit_areal_density = 10*PhysicalUnits::g_per_cm2;
  curve.phys_model_external_atten.push_back( ext_atten );
  curve.phys_model_use_hoerl = false;
  
  curve.pu242_correlation_method = RelActCalc::PuCorrMethod::ByPu239Only;

  RelActCalcAuto::NucInputInfo u235_input;
  u235_input.source = db->nuclide("U235");
  u235_input.age = 20.0*PhysicalUnits::year;
  u235_input.fit_age = false;
  u235_input.gammas_to_exclude = { 26.325 };
  u235_input.peak_color_css = "red";
  curve.nuclides.push_back( u235_input );

  RelActCalcAuto::NucInputInfo u238_input;
  u238_input.source = db->nuclide("U238");
  u238_input.age = 18.0*PhysicalUnits::year;
  u238_input.fit_age = true;
  u238_input.gammas_to_exclude = { 1001.1 };
  u238_input.peak_color_css = "blue";
  curve.nuclides.push_back( u238_input );

  
  state.options.rel_eff_curves.push_back( curve );

  RelActCalcAuto::RoiRange range;
  range.lower_energy = 500.0;
  range.upper_energy = 2000.0;
  range.continuum_type = PeakContinuum::OffsetType::Quadratic;
  range.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::CanBeBrokenUp;
  state.options.rois.push_back( range );

  // Add a second range
  range.lower_energy = 50.0;
  range.upper_energy = 200.0;
  range.continuum_type = PeakContinuum::OffsetType::Linear;
  range.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::CanExpandForFwhm;
  state.options.rois.push_back( range );

  RelActCalcAuto::FloatingPeak peak;
  peak.energy = 1000.0;
  peak.release_fwhm = true;
  peak.energy_origin = RelActCalcAuto::FloatingPeak::EnergyType::Known;
  state.options.floating_peaks.push_back( peak );

  state.background_subtract = true;
  
  state.show_ref_lines = true;
  state.lower_display_energy = 50.0;
  state.upper_display_energy = 2000.0;
  
  rapidxml::xml_document<char> doc;
  auto rel_act_node = state.serialize( &doc );
  
  RelActCalcAuto::RelActAutoGuiState state2;
  try
  {
    state2.deSerialize( rel_act_node );
  }catch( std::exception &e )
  {
    cerr << "Failed to deserialize state: " << e.what() << endl;
    return;
  }

  try
  {
    RelActCalcAuto::RelActAutoGuiState::equalEnough( state, state2 );
  }catch( std::exception &e )
  {
    cerr << "Failed to compare states: " << e.what() << endl;
    assert( 0 );
    exit( -1 );
    return;
  }

  cout << "States are equal" << endl;

// Now try multiple rel eff curves
  RelActCalcAuto::RelEffCurveInput curve2;
  curve2.nucs_of_el_same_age = true;
  curve2.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnXLnY;
  curve2.rel_eff_eqn_order = 2;

  RelActCalcAuto::NucInputInfo co60_input;
  co60_input.source = db->nuclide("Co60");
  co60_input.age = 18.0*PhysicalUnits::year;
  co60_input.fit_age = true;
  co60_input.gammas_to_exclude = { 1001.1 };
  co60_input.peak_color_css = "blue";
  curve2.nuclides.push_back( co60_input );

  state.options.rel_eff_curves.push_back( curve2 );

  auto rel_act_node2 = state.serialize( doc.first_node() );
  
  RelActCalcAuto::RelActAutoGuiState state3;  
  try
  {
    state3.deSerialize( rel_act_node2 );
  }catch( std::exception &e )
  {
    cerr << "Failed to deserialize state: " << e.what() << endl;
    return;
  }
  
  try
  {
    RelActCalcAuto::RelActAutoGuiState::equalEnough( state, state3 );
  }catch( std::exception &e )
  {
    cerr << "Failed to compare states: " << e.what() << endl;
    
    assert( 0 );
    exit( -1 );
    return;
  }
  
  cout << "States are equal" << endl;
}//void check_auto_state_xml_serialization()
#endif //PERFORM_DEVELOPER_CHECKS 
  
void check_physical_model_eff_function()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  
  const string data_dir = InterSpec::staticDataDirectory();
  assert( matdb );
  
  DetectorPeakResponse det;
  
  try
  {
    const string drf_dir = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GenericGadrasDetectors/HPGe 40%" );
    det.fromGadrasDirectory( drf_dir );
  }catch( std::exception &e )
  {
    cerr << "Failed to open directory." << endl;
    return;
  }
  
  const shared_ptr<const Material> pu = matdb->material("Pu (plutonium)");
  assert( pu );
  
  RelActCalc::PhysModelShield<double> self_atten;
  self_atten.material = pu;
  self_atten.areal_density = 1.98*PhysicalUnits::g_per_cm2;
  
  
  RelActCalc::PhysModelShield<double> ext_atten;
  ext_atten.material = matdb->material("stainless-steel NIST");
  ext_atten.areal_density = 0.8*PhysicalUnits::g_per_cm2;
  
  RelActCalc::PhysModelShield<double> generic_atten;
  generic_atten.atomic_number = 32.0; //generic AN - e.g., Germanium
  generic_atten.areal_density = 0.532*PhysicalUnits::g_per_cm2; //generic AD - e.g., 1 mm of Ge
  
  vector<RelActCalc::PhysModelShield<double>> external_attenuations;
  external_attenuations.push_back( ext_atten );
  external_attenuations.push_back( generic_atten );
  
  const double hoerl_b = 1.0; // Modified Hoerl b: pow(0.001*energy,b)
  const double hoerl_c = 1.0; // Modified Hoerl c: pow(c,1000/energy)
  
  
  cout << "Energy (keV), Counts" << endl;
  for( double energy = 50.0; energy < 3000; energy += 2 )
  {
    double eff = RelActCalc::eval_physical_model_eqn( energy, self_atten, external_attenuations,
                               &det, hoerl_b, hoerl_c );

    cout << energy << "," << eff << endl;
  }
}//void check_physical_model_eff_function()
  
void example_manual_phys_model()
{
  cout << "Running example_manual_phys_model" << endl;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  
  const string data_dir = InterSpec::staticDataDirectory();
  assert( matdb );
  
  const shared_ptr<const Material> uranium = matdb->material("U (uranium)");
  assert( uranium );
  

  const string input_file = "uranium_40%_HPGe_15cm_peaks_fit.n42";
  SpecMeas infile;
  const bool loaded = infile.load_file( input_file, SpecUtils::ParserType::Auto );
  if( !loaded )
  {
    cerr << "Failed to load '" << input_file << "', aborting." << endl;
    return;
  }
  std::shared_ptr<const DetectorPeakResponse> det = infile.detector();
  assert( det );
  
  shared_ptr<deque<shared_ptr<const PeakDef>>> orig_peaks = infile.peaks( {1} );
  assert( orig_peaks && orig_peaks->size() );

  assert( infile.detector_names().size() == 1 );
  const shared_ptr<const SpecUtils::Measurement> meas = infile.measurement( 1, infile.detector_names()[0] );
  assert( meas );
  
  const double base_rel_eff_uncert = 0.01;

  std::vector<RelActCalcManual::GenericPeakInfo> peaks;
  for( const auto &peak : *orig_peaks )
  {
    const SandiaDecay::Nuclide * const nuclide = peak->parentNuclide();
    if( !nuclide || (nuclide->atomicNumber != 92) )
      continue;
    
    RelActCalcManual::GenericPeakInfo info;
    info.m_energy = peak->gammaParticleEnergy();
    info.m_mean = peak->mean();
    info.m_fwhm = peak->fwhm();
    info.m_counts = peak->peakArea();
    info.m_counts_uncert = peak->peakAreaUncert();
    info.m_base_rel_eff_uncert = base_rel_eff_uncert;
      
    peaks.push_back( info );
  }//for( const auto &peak : *orig_peaks )

  vector<RelActCalcManual::PeakCsvInput::NucAndAge> isotopes;
  isotopes.emplace_back( "U235", 20.0*PhysicalUnits::year, false );
  isotopes.emplace_back( "U238", 20.0*PhysicalUnits::year, false );
  isotopes.emplace_back( "U234", 20.0*PhysicalUnits::year, false );

  
  const vector<pair<float,float>> energy_ranges{ { 50.0, 2000.0 } };
  const float energy_tolerance_sigma = 1.5;
  const std::vector<float> excluded_peak_energies;
  const float measurement_duration = meas->real_time();

  RelActCalcManual::PeakCsvInput::NucMatchResults matched_res 
              = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( peaks, RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay,
                                     energy_ranges, isotopes,
                                     energy_tolerance_sigma,
                                     excluded_peak_energies,
                                     measurement_duration );

  RelActCalcManual::RelEffInput rel_eff_solve_input;
  rel_eff_solve_input.peaks = matched_res.peaks_matched;
  rel_eff_solve_input.eqn_form = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  rel_eff_solve_input.eqn_order = 0;
  rel_eff_solve_input.phys_model_detector = det;
  rel_eff_solve_input.phys_model_use_hoerl = true;
  rel_eff_solve_input.use_ceres_to_fit_eqn = true;

  RelActCalc::PhysicalModelShieldInput self_atten_def;
  //self_atten_def.atomic_number = 80;
  self_atten_def.material = uranium;
  self_atten_def.areal_density = 1.25*PhysicalUnits::g_per_cm2;
  self_atten_def.fit_atomic_number = false;
  //self_atten_def.lower_fit_atomic_number = 1.0;
  //self_atten_def.upper_fit_atomic_number = 98.0;
  self_atten_def.fit_areal_density = true;
  self_atten_def.lower_fit_areal_density = 0.0;
  self_atten_def.upper_fit_areal_density = 500*PhysicalUnits::g_per_cm2;

  rel_eff_solve_input.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( self_atten_def );
/*
  RelActCalc::PhysicalModelShieldInput external_atten_def;
  external_atten_def.atomic_number = 26;
  external_atten_def.material = nullptr;
  external_atten_def.areal_density = 1.0*PhysicalUnits::g_per_cm2;
  external_atten_def.fit_atomic_number = false;
  external_atten_def.lower_fit_atomic_number = 1.0;
  external_atten_def.upper_fit_atomic_number = 98.0;
  external_atten_def.fit_areal_density = true;
  external_atten_def.lower_fit_areal_density = 0.0;
  external_atten_def.upper_fit_areal_density = 500*PhysicalUnits::g_per_cm2;
  //rel_eff_solve_input.phys_model_external_attens.push_back( make_shared<RelActCalc::PhysicalModelShieldInput>( external_atten_def ) );
*/

  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  assert( u235 && u238 );

  RelActCalcManual::ManualActRatioConstraint nuc_constraint;
  nuc_constraint.m_constrained_nuclide = "U235";
  nuc_constraint.m_controlling_nuclide = "U238";
  nuc_constraint.m_constrained_to_controlled_activity_ratio =  RelActCalc::mass_ratio_to_act_ratio(u235, u238, 0.2);
  rel_eff_solve_input.act_ratio_constraints.push_back( nuc_constraint );

  RelActCalcManual::RelEffSolution sol = RelActCalcManual::solve_relative_efficiency( rel_eff_solve_input );

  if( sol.m_status != RelActCalcManual::ManualSolutionStatus::Success )
  {
    cerr << "Failed to solve relative efficiency: " << sol.m_error_message << endl;
    return;
  }

  ofstream out_html( "manual_phys_model_result.html" );
  std::vector<std::shared_ptr<const PeakDef>> disp_peaks( orig_peaks->begin(), orig_peaks->end() );
  std::vector<std::shared_ptr<const PeakDef>> background_peaks;

  sol.print_html_report( out_html, "Manual Phys Model", meas, disp_peaks, nullptr, 0.0, background_peaks );


  cout << "Solution: " << endl;
  sol.print_summary( cout );
  cout << endl;
}//void example_manual_phys_model()

void run_u02_example()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const string data_dir = InterSpec::staticDataDirectory();
  assert( matdb );
  
  const string specfilename = "mixed_U02_sample.pcf";
  SpecUtils::SpecFile specfile;
  const bool loaded_spec = specfile.load_file(specfilename, SpecUtils::ParserType::Auto );
  if( !loaded_spec )
  {
    cerr << "Failed to load '" << specfilename << "', aborting." << endl;
    return;
  }

  //Sample 1 is background
  //Sample 38 is UO2_50%_50%
  assert( specfile.num_measurements() == 38 );
  shared_ptr<const SpecUtils::Measurement> background = specfile.measurement_at_index(0);
  assert( background );
  assert( background->title() == "Background" );
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement_at_index(37);
  assert( foreground );
  assert( SpecUtils::istarts_with( foreground->title(), "U02_50%_50% @ 25 cm H=100 cm" ) );
  
  
  const string setup_xml_path = "isotopics_by_nuclides_mixed_U02_sample-2.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );
  
  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );

  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );
  
  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node );
  
  //RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
  //constraint.nuclide = db->nuclide("U235");
  //constraint.lower_mass_fraction = 0.1;
  //constraint.upper_mass_fraction = 0.4;
  //state.options.rel_eff_curves[0].mass_fraction_constraints.push_back( constraint );



  const string detector_xml_path = "Detective-X_GADRAS_drf.xml";
  rapidxml::file<char> detector_input_file( detector_xml_path.c_str() );
  
  rapidxml::xml_document<char> detector_doc;
  detector_doc.parse<rapidxml::parse_trim_whitespace>( detector_input_file.data() );

  const rapidxml::xml_node<char> *detector_base_node = detector_doc.first_node( "DetectorPeakResponse" );
  assert( detector_base_node );
  
  
  auto det = make_shared<DetectorPeakResponse>();
  det->fromXml( detector_base_node );
  
  vector<shared_ptr<const PeakDef>> all_peaks{};
  
  const double start_cpu = SpecUtils::get_cpu_time();
  const double start_wall = SpecUtils::get_wall_time();
  
  const RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve( state.options,
                                              foreground, background, det, all_peaks,
                                              PeakFitUtils::CoarseResolutionType::High, nullptr );

  const double end_cpu = SpecUtils::get_cpu_time();
  const double end_wall = SpecUtils::get_wall_time();
  
  ofstream out_html( "U02_result.html" );
  solution.print_summary( cout );
  solution.print_html_report( out_html );
  
  pair<double,optional<double>> enrich = solution.mass_enrichment_fraction( db->nuclide("U235"), 0);
  
  cout << "Enrichment: " << enrich.first
  << " +- " << (enrich.second.has_value() ? std::to_string(enrich.second.value()) : string("UncertCalcError"))
  << endl;

  cout << "Took:\n"
  << "\tNum Function Calls: " << solution.m_num_function_eval_solution << endl
  << "\tFcn Evals total:    " << solution.m_num_function_eval_total << endl
  << "\tSeconds solving:    " << 1.0E-6*solution.m_num_microseconds_eval << endl
  << "\tSeconds in eval:    " << 1.0E-6*solution.m_num_microseconds_in_eval
    << " (" << (100.0*solution.m_num_microseconds_eval)/solution.m_num_microseconds_in_eval << "%)\n"
  << "\tMicrSec per eval:   " << ((1.0*solution.m_num_microseconds_in_eval)/solution.m_num_function_eval_solution) << "\n"
  << "\tWall Total (s):     " << (end_wall - start_wall) << endl
  << "\tCPU Total (s):      " << (end_cpu - start_cpu) << endl
  << endl;
  
  /* Compiled in Release mode on M1 MBP, 3 runs (not sure where randomness comes in...)
   
   20250127T11:30, git hash 94342ca74a9ca357bc8af709c46707279ccc7c49:
   Enrichment: 0.4355
   Took:
     Num Function Calls: 1837
     Fcn Evals total:    1872
     Seconds solving:    77.29
     Seconds in eval:    75.93 (101.8%)
     MicrSec per eval:   4.134e+04
     Wall Total (s):     77.29
     CPU Total (s):      15.52
   BUT: re-running resulted in different enrichments!
   
   After getting rid of the `pool.join()` that caused all the delay (about 25x faster, per evaluation):
   Enrichment: 0.4369
   Took:
     Num Function Calls: 2013
     Fcn Evals total:    2048
     Seconds solving:    4.429
     Seconds in eval:    3.463 (127.9%)e
     MicrSec per eval:   1720
     Wall Total (s):     4.43
     CPU Total (s):      13.9
   STILL gives different enrichments between runs
   
   20250128, git hash 8a6fe8550c86f0e422212a446d597dafb32dd30d:
   After changing to use Eigen to fit continuum, and caching the peak contributions in each channel
   (i.e., implementing and using `peaks_for_energy_range_imp(...)`)
   Enrichment: 0.4305
   Took:
     Num Function Calls: 2611
     Fcn Evals total:    2646
     Seconds solving:    0.591
     Seconds in eval:    0.3942 (149.9%)
     MicrSec per eval:   151
     Wall Total (s):     0.591
     CPU Total (s):      1.607
   Now gives same enrichment every time!
   Now takes 0.35% as long as originally!
   
   20250129
   After enabling auto-differentiation (e.g. `Jet<>`), and using stride of 4
   Enrichment: 0.4219
   Took:
     Num Function Calls: 187
     Fcn Evals total:    192
     Seconds solving:    0.1727
     Seconds in eval:    0.02955 (584.5%)
     MicrSec per eval:   158
     Wall Total (s):     0.1728
     CPU Total (s):      0.368
   Gives same enrichment every time, and maybe the Chi2 is better
   
   using stride = 8
   Enrichment: 0.4219
   Took:
     Num Function Calls: 119
     Fcn Evals total:    122
     Seconds solving:    0.1788
     Seconds in eval:    0.02016 (886.9%)
     MicrSec per eval:   169.4
     Wall Total (s):     0.1789
     CPU Total (s):      0.3668
   
   using stride = 16
   Enrichment: 0.4219
   Took:
     Num Function Calls: 85
     Fcn Evals total:    87
     Seconds solving:    0.1665
     Seconds in eval:    0.01722 (966.8%)
     MicrSec per eval:   202.6
     Wall Total (s):     0.1665
     CPU Total (s):      0.3361
   
   using stride = 32
   Enrichment: 0.4219
   Took:
     Num Function Calls: 51
     Fcn Evals total:    52
     Seconds solving:    0.1608
     Seconds in eval:    0.0118 (1363%)
     MicrSec per eval:   231.3
     Wall Total (s):     0.1608
     CPU Total (s):      0.3306
   */
}//void run_u02_example()
  






void czt_pu_example()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const string data_dir = InterSpec::staticDataDirectory();
  assert( matdb );
  
  const string specfilename = "czt_600keV__Pu84_60mm_1mmCd_ex.n42";
  
  SpecMeas specfile;
  const bool loaded_spec = specfile.load_file(specfilename, SpecUtils::ParserType::Auto );
  if( !loaded_spec )
  {
    cerr << "Failed to load '" << specfilename << "', aborting." << endl;
  }
  
  
  assert( specfile.num_measurements() == 1 );
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement_at_index( size_t(0) );
  assert( foreground );

  
  const string setup_xml_path = "isotopics_by_nuclides_CZT500_Pu84_@60mm_1mmCd_04-2.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );
  
  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );
  
  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );
  
  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node );
  
  assert( state.options.rel_eff_curves.size() == 1 );

  
  auto det = specfile.detector();
  
  
  //state.options.rel_eff_curves.resize(1);
  state.options.rel_eff_curves[0].phys_model_use_hoerl = false;
  
  
  // We can constrain the RelActivity.
  for( auto &nuc : state.options.rel_eff_curves[0].nuclides )
  {
    if( RelActCalcAuto::to_name(nuc.source) == "Pu239" )
    {
      nuc.starting_rel_act = 750*16228.0;
      nuc.min_rel_act = nuc.max_rel_act = nuc.starting_rel_act;
    }else if( RelActCalcAuto::to_name(nuc.source) == "Am241" )
    {
      nuc.starting_rel_act = 5*96521.1;
      nuc.min_rel_act = nuc.max_rel_act = nuc.starting_rel_act;
      
    }
  }


  //det->setFwhmCoefficients( {}, DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm );
  state.options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_4;
  state.options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;

  // Check serialization to/from XML
#if( PERFORM_DEVELOPER_CHECKS )
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *base_node = doc.allocate_node( rapidxml::node_element, "State" );
    doc.append_node( base_node );

    rapidxml::xml_node<char> *state_node = state.serialize( base_node );
    assert( state_node );

    RelActCalcAuto::RelActAutoGuiState state_cpy;
    state_cpy.deSerialize( state_node );

    RelActCalcAuto::RelActAutoGuiState::equalEnough( state, state_cpy );
  } 
#endif
  
  vector<shared_ptr<const PeakDef>> all_peaks{};
  const RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve( state.options,
                                                foreground, nullptr, det, all_peaks,
                                                PeakFitUtils::CoarseResolutionType::CZT, nullptr );
  ofstream out_html( "czt_Pu_rel_eff_result.html" );
  solution.print_summary( cout );
  solution.print_html_report( out_html );
  
  /*
  for( size_t i = 0; i < solution.m_rel_activities.size(); ++i )
  {
    const double enrichment = solution.mass_enrichment_fraction( u235, i );
    const double u235_counts = solution.nuclide_counts( u235, i );
    const double u238_counts = solution.nuclide_counts( u238, i );
    cout << "Enrichment " << i << std::left << ": " << setprecision(6) << setw(11) << enrichment
    << ", counts(u235)=" << setw(11) << u235_counts
    << ", counts(u238)=" << setw(11) << u238_counts << endl;
  }
  cout << "For sample: " << title << endl;
  */
}//void czt_pu_example()




void check_auto_nuclide_constraints_checks()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  const SandiaDecay::Nuclide * const u234 = db->nuclide("U234");
  
  RelActCalcAuto::RelEffCurveInput rel_eff_curve;

  RelActCalcAuto::NucInputInfo u235_input;
  u235_input.source = u235;
  u235_input.age = 20.0 * PhysicalUnits::year;
  u235_input.fit_age = false;
  u235_input.gammas_to_exclude = {};
  u235_input.peak_color_css = "rgb(0, 0, 255)";
  
  RelActCalcAuto::NucInputInfo u238_input;
  u238_input.source = u238;
  u238_input.age = 20.0 * PhysicalUnits::year;
  u238_input.fit_age = false;
  u238_input.gammas_to_exclude = {};
  u238_input.peak_color_css = "rgb(255, 69, 0)";

  RelActCalcAuto::NucInputInfo u234_input;
  u234_input.source = u234;
  u234_input.age = 20.0 * PhysicalUnits::year;
  u234_input.fit_age = false;
  u234_input.gammas_to_exclude = {};
  u234_input.peak_color_css = "rgb(34, 139, 34)";

  rel_eff_curve.nuclides.push_back( u235_input );
  rel_eff_curve.nuclides.push_back( u238_input );
  rel_eff_curve.nuclides.push_back( u234_input );

  rel_eff_curve.nucs_of_el_same_age = true;
  rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  rel_eff_curve.rel_eff_eqn_order = 0;
  
  // Check we dont throw an error for no constraints
  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;

    rel_eff_cpy.check_nuclide_constraints();
  }catch( std::exception &e )
  {
    cerr << "Failed constraint check when we shouldnt have: Error: " << e.what() << endl;
    assert( 0 );
  }

  // Check we dont throw an error for a valid constraint
  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.0072/(1.0 - 0.0072) ) );
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u238, u234, 0.00001 ) ); 

    rel_eff_cpy.check_nuclide_constraints();
  }catch( std::exception &e )
  {
    cerr << "Failed constraint check when we shouldnt have: Error: " << e.what() << endl;
    assert( 0 );
  }

  // Check we throw an error for a cycle
  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;

    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.0072/(1.0 - 0.0072) ) );
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u238, u234, 0.00001 ) ); 
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u234, u235, 0.00001 ) ); 
    
    rel_eff_cpy.check_nuclide_constraints();

    cerr << "Failed to detect cycle in nuclide constraints" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here because of the cycle
  }

  // Check we throw an error for a duplicate nuclide constraint
  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.02) );
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.00001 ) ); 

    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect duplicate nuclide constraints" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  // Check we throw an error for an invalid nuclide in a constraint
  try
  {
    const SandiaDecay::Nuclide * const u232 = db->nuclide("U232");

    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u232, u235, 0.00001 ) ); 
    
    rel_eff_cpy.check_nuclide_constraints();

    cerr << "Failed to detect invalid nuclide in constraint" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }

  // Check we throw an error for an invalid activity ratio
  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::ActRatioConstraint nuc_constraint;
    nuc_constraint.constrained_source = u235;
    nuc_constraint.controlling_source = u238;
    nuc_constraint.constrained_to_controlled_activity_ratio = -1.0;
    
    rel_eff_cpy.act_ratio_constraints.push_back( nuc_constraint );
    
    rel_eff_cpy.check_nuclide_constraints();

    cerr << "Failed to detect invalid activity ratio" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }

  // Check we throw an error for a short cycle
  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;

    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.0072/(1.0 - 0.0072) ) );
    rel_eff_cpy.act_ratio_constraints.push_back( RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u238, u235, 0.00001 ) ); 
    
    rel_eff_cpy.check_nuclide_constraints();

    cerr << "Failed to detect cycle in nuclide constraints" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here because of the cycle
  }


  // Now check mass fraction constraints
  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 0.0072;
    constraint.upper_mass_fraction = 0.0072;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    
    rel_eff_cpy.check_nuclide_constraints();
    
    // We are suppoest to get here
  }catch(const std::exception& e)
  {
    cerr << "Valid mass fraction constraint erroneously detected: Error: " << e.what() << endl;
    assert( 0 );
  }

  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;

    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = db->nuclide("U233");
    constraint.lower_mass_fraction = 0.0072;
    constraint.upper_mass_fraction = 0.0072;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    
    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect invalid nuclide in mass fraction constraint" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }


  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 0.0072;
    constraint.upper_mass_fraction = 1.00001;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    
    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect invalid nuclide mass fraction" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }

  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = -0.0001;
    constraint.upper_mass_fraction = 0.8;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    
    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect invalid nuclide mass fraction" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }

  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 1.0;
    constraint.upper_mass_fraction = 1.0;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    
    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect invalid nuclide mass fraction" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }


  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 0.0;
    constraint.upper_mass_fraction = 0.5;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    
    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect multiple mass fraction constraints for same nuclide" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }


  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 0.0;
    constraint.upper_mass_fraction = 0.5;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );

    constraint.nuclide = u238;
    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    
    constraint.nuclide = u234;
    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    

    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect all nuclides of element being mass fraction constrained" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }


  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 0.8;
    constraint.upper_mass_fraction = 0.9;

    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );

    constraint.nuclide = u234;
    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );
    

    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect sum of lower mass fraction values being larger than 1.0" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }


  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 0.1;
    constraint.upper_mass_fraction = 0.9;
    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );

    RelActCalcAuto::RelEffCurveInput::ActRatioConstraint nuc_constraint;
    nuc_constraint.constrained_source = u235;
    nuc_constraint.controlling_source = u238;
    nuc_constraint.constrained_to_controlled_activity_ratio = 0.1;
    rel_eff_cpy.act_ratio_constraints.push_back( nuc_constraint );

    rel_eff_cpy.check_nuclide_constraints();
    
    cerr << "Failed to detect nuclide being both activity and mass-fraction constrained." << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here
  }


  try
  {
    RelActCalcAuto::RelEffCurveInput rel_eff_cpy = rel_eff_curve;


    RelActCalcAuto::NucInputInfo co60_input;
    co60_input.source = db->nuclide("Co60");
    co60_input.age = 20.0 * PhysicalUnits::year;
    co60_input.fit_age = false;
    co60_input.gammas_to_exclude = {};
    rel_eff_cpy.nuclides.push_back( co60_input );


     RelActCalcAuto::NucInputInfo cs137_input;
    cs137_input.source = db->nuclide("Cs137");
    cs137_input.age = 20.0 * PhysicalUnits::year;
    cs137_input.fit_age = false;
    cs137_input.gammas_to_exclude = {};
    rel_eff_cpy.nuclides.push_back( cs137_input );

    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = u235;
    constraint.lower_mass_fraction = 0.1;
    constraint.upper_mass_fraction = 0.9;
    rel_eff_cpy.mass_fraction_constraints.push_back( constraint );

    RelActCalcAuto::RelEffCurveInput::ActRatioConstraint nuc_constraint;
    nuc_constraint.constrained_source = co60_input.source;
    nuc_constraint.controlling_source = u238;
    nuc_constraint.constrained_to_controlled_activity_ratio = 0.1;
    rel_eff_cpy.act_ratio_constraints.push_back( nuc_constraint );

    nuc_constraint.constrained_source = cs137_input.source;
    nuc_constraint.controlling_source = u235;
    nuc_constraint.constrained_to_controlled_activity_ratio = 0.1;
    rel_eff_cpy.act_ratio_constraints.push_back( nuc_constraint );

    rel_eff_cpy.check_nuclide_constraints();

    // We are suppoest to get here
  }catch(const std::exception& e)
  {
    cerr << "Failed to allow mass-constrained nuclide activity ratio control nuclide of other element." << endl;
    assert( 0 );
  }


  cout << "All act_ratio_constraints checks passed" << endl;
}//void check_auto_nuclide_constraints_checks()
  
void check_manual_nuclide_constraints_checks()
{
  
  RelActCalcManual::RelEffInput input;

  RelActCalcManual::GenericPeakInfo peak_info;
  peak_info.m_energy = 124.8;
  peak_info.m_counts = 1000.0;
  peak_info.m_counts_uncert = 10.0;

  RelActCalcManual::GenericLineInfo line_info;
  line_info.m_isotope = "U238";
  line_info.m_yield = 1.0;
  peak_info.m_source_gammas.push_back( line_info );
  input.peaks.push_back( peak_info );

  peak_info.m_energy = 185.0;
  peak_info.m_counts = 1000.0;
  peak_info.m_counts_uncert = 10.0;
  line_info.m_isotope = "U235";
  line_info.m_yield = 1.0;
  peak_info.m_source_gammas.clear();
  peak_info.m_source_gammas.push_back( line_info );
  input.peaks.push_back( peak_info );

  // Check we dont throw an error for no constraints
  try
  {
    RelActCalcManual::RelEffInput input_cpy = input;

    input_cpy.check_nuclide_constraints();
  }catch( std::exception &e )
  {
    cerr << "Failed constraint check when we shouldnt have: Error: " << e.what() << endl;
    assert( 0 );
  }

  // Try a valid constraint
  try
  {
    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.0072/(1.0 - 0.0072);

    RelActCalcManual::RelEffInput input_cpy = input;
    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    input_cpy.check_nuclide_constraints(); //shouldnt be any problems
  }catch( std::exception &e ) 
  {
    cerr << "Failed to create valid nuclide constraint: Error: " << e.what() << endl;
    assert( 0 );
  }

  // Make sure invalid constrained nuclide throws
  try
  {
    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "Co60";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.1;

    RelActCalcManual::RelEffInput input_cpy = input;
    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    input_cpy.check_nuclide_constraints(); //shouldnt be any problems

    cerr << "Failed to catch invalid constrained nuclide in constraint" << endl;
    assert( 0 );
  }catch( std::exception &e ) 
  {
    // We are suppoest to get here
  }

  // Make sure invalid controlling nuclide throws
  try
  {
    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "Co60";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.1;

    RelActCalcManual::RelEffInput input_cpy = input;
    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    input_cpy.check_nuclide_constraints(); //shouldnt be any problems

    cerr << "Failed to catch invalid constrilling nuclide in constraint" << endl;
    assert( 0 );
  }catch( std::exception &e ) 
  {
    // We are suppoest to get here
  }

  
  // Check we throw an error for a cycle A->B->C->A
  try
  {
    RelActCalcManual::RelEffInput input_cpy = input;

    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.1;

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );
    
    nuc_constraint.m_constrained_nuclide = "U238";
    nuc_constraint.m_controlling_nuclide = "U234";

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    nuc_constraint.m_constrained_nuclide = "U234";
    nuc_constraint.m_controlling_nuclide = "U235";

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    // Add a peak for U234
    peak_info.m_energy = 120.0;
    peak_info.m_counts = 100.0;
    peak_info.m_counts_uncert = 10.0;
    line_info.m_isotope = "U234";
    line_info.m_yield = 1.0;
    peak_info.m_source_gammas.clear();
    peak_info.m_source_gammas.push_back( line_info );
    input_cpy.peaks.push_back( peak_info );
    
    input_cpy.check_nuclide_constraints();

    cerr << "Failed to detect cycle in nuclide constraints" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here because of the cycle
  }

  // Check we throw an error for a cycle A->B->A
  try
  {
    RelActCalcManual::RelEffInput input_cpy = input;

    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.1;

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );
    
    nuc_constraint.m_constrained_nuclide = "U238";
    nuc_constraint.m_controlling_nuclide = "U235";

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    input_cpy.check_nuclide_constraints();

    cerr << "Failed to detect cycle in nuclide constraints" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here because of the cycle
  }

  // Check we throw an error for a duplicate nuclide constraint
  try
  {
    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.0072/(1.0 - 0.0072);

    RelActCalcManual::RelEffInput input_cpy = input;
    input_cpy.act_ratio_constraints.push_back( nuc_constraint );
    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    input_cpy.check_nuclide_constraints(); //shouldnt be any problems

    cerr << "Failed to detect duplicate constraints" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  // Check we throw an error for nuclide with two constraints
  try
  {
    RelActCalcManual::RelEffInput input_cpy = input;

    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.1;

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );
    
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U234";

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    // Add a peak for U234
    peak_info.m_energy = 120.0;
    peak_info.m_counts = 100.0;
    peak_info.m_counts_uncert = 10.0;
    line_info.m_isotope = "U234";
    line_info.m_yield = 1.0;
    peak_info.m_source_gammas.clear();
    peak_info.m_source_gammas.push_back( line_info );
    input_cpy.peaks.push_back( peak_info );
    
    input_cpy.check_nuclide_constraints();

    cerr << "Failed to detect nuclide with two constraints" << endl;
    assert( 0 );
  }catch(const std::exception& e)
  {
    // We are suppoest to get here because of the cycle
  }


  // Check we done throw an error for a nuclide controlling two other nuclides
  try
  {
    RelActCalcManual::RelEffInput input_cpy = input;

    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = 0.1;

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );
    
    nuc_constraint.m_constrained_nuclide = "U234";
    nuc_constraint.m_controlling_nuclide = "U238";

    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    // Add a peak for U234
    peak_info.m_energy = 120.0;
    peak_info.m_counts = 100.0;
    peak_info.m_counts_uncert = 10.0;
    line_info.m_isotope = "U234";
    line_info.m_yield = 1.0;
    peak_info.m_source_gammas.clear();
    peak_info.m_source_gammas.push_back( line_info );
    input_cpy.peaks.push_back( peak_info );
    
    input_cpy.check_nuclide_constraints();
  }catch(const std::exception& e)
  {
    cerr << "Failed to have two nuclides controlled by the same nuclide" << endl;
    assert( 0 );
  }

  // Check we throw an error for an invalid activity ratio
  try
  {
    RelActCalcManual::ManualActRatioConstraint nuc_constraint;
    nuc_constraint.m_constrained_nuclide = "U235";
    nuc_constraint.m_controlling_nuclide = "U238";
    nuc_constraint.m_constrained_to_controlled_activity_ratio = -1.0;

    RelActCalcManual::RelEffInput input_cpy = input;
    input_cpy.act_ratio_constraints.push_back( nuc_constraint );

    input_cpy.check_nuclide_constraints(); //shouldnt be any problems

    cerr << "Failed to detect invalid activity ratio" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  cout << "All manual act_ratio_constraints checks passed" << endl;
}//void check_manual_nuclide_constraints_checks()


void check_auto_hoerl_and_ext_shield_checks()
{
  RelActCalcAuto::Options options;
  options.same_hoerl_for_all_rel_eff_curves = true;
  options.same_external_shielding_for_all_rel_eff_curves = true;
  
  RelActCalcAuto::RelEffCurveInput rel_eff_curve;
  rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  rel_eff_curve.phys_model_use_hoerl = true;

  auto ext_shield = make_shared<RelActCalc::PhysicalModelShieldInput>( RelActCalc::PhysicalModelShieldInput() );  
  ext_shield->atomic_number = 26;
  ext_shield->material = nullptr;
  ext_shield->areal_density = 1.0 * PhysicalUnits::g_per_cm2;
  ext_shield->fit_atomic_number = false;
  ext_shield->fit_areal_density = true;
  ext_shield->lower_fit_areal_density = 0.0;
  ext_shield->upper_fit_areal_density = 500.0 * PhysicalUnits::g_per_cm2;
  rel_eff_curve.phys_model_external_atten.push_back( ext_shield );
 
  // Check we dont throw error if not sharing hoerl and external shielding
  try
  {
    RelActCalcAuto::Options options_cpy = options;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_2 = rel_eff_curve;
    options_cpy.same_hoerl_for_all_rel_eff_curves = false;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = false;
    options_cpy.check_same_hoerl_and_external_shielding_specifications();

    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();

    options_cpy.rel_eff_curves.push_back( rel_eff_curve_2 );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();
  }catch( std::exception &e )
  {
    cerr << "Failed to check same hoerl and external shielding specifications: " << e.what() << endl;
    assert( 0 );
  }
  
  // Check we dont throw an error for two rel_eff_curves with the same external shielding, 
  //  and sharing the same hoerl function and external shielding
  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = true;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = true;
    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();
  }catch( std::exception &e )
  {
    cerr << "Failed to check same hoerl and external shielding specifications: " << e.what() << endl;
    assert( 0 );
  }

  // Check we throw an error if we only have a single rel_eff_curve, but ask to share the same hoerl
  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = true;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = false;

    options_cpy.rel_eff_curves.clear();
    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();
    cerr << "Failed to throw an error for a single rel_eff_curve" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  // Check we throw an error if we only have a single rel_eff_curve, but ask to share the same external shielding
  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = false;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = true;

    options_cpy.rel_eff_curves.clear();
    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();
    cerr << "Failed to throw an error for a single rel_eff_curve" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  // Check we throw an error if we have two rel_eff_curves with different external shielding
  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = false;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = true;

    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    
    RelActCalcAuto::RelEffCurveInput rel_eff_curve_2;
    rel_eff_curve_2.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
    rel_eff_curve_2.phys_model_use_hoerl = true;
    
    auto ext_shield_2 = make_shared<RelActCalc::PhysicalModelShieldInput>( RelActCalc::PhysicalModelShieldInput() );  
    ext_shield_2->atomic_number = 52;
    ext_shield_2->material = nullptr;
    ext_shield_2->areal_density = 2.5 * PhysicalUnits::g_per_cm2;
    ext_shield_2->fit_atomic_number = false;
    ext_shield_2->fit_areal_density = true;
    ext_shield_2->lower_fit_areal_density = 0.0;
    ext_shield_2->upper_fit_areal_density = 500.0 * PhysicalUnits::g_per_cm2;
    rel_eff_curve_2.phys_model_external_atten.push_back( ext_shield_2 );
    
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_2 );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();
    cerr << "Failed to throw an error for different external shielding" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }
    
    
  // Check we throw an error if one of RelEff curves does not use a hoerl function
  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = true;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = false;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy = rel_eff_curve;
    rel_eff_curve_cpy.phys_model_use_hoerl = false;
    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy );

    options_cpy.check_same_hoerl_and_external_shielding_specifications();
    cerr << "Failed to throw an error for one rel_eff_curve not using a hoerl function" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  // Check we throw an error if we request to share external shieldings, but one of the 
  //  rel_eff_curves does not have an external shielding
  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = false;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = true;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy = rel_eff_curve;
    rel_eff_curve_cpy.phys_model_external_atten.clear();
    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();
    cerr << "Failed to throw an error for one rel_eff_curve not having an external shielding" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  // Check we dont throw an error if we request to share external shieldings, but one of the 
  //  rel_eff_curves is not a physical model
  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = false;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = true;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy = rel_eff_curve;
    rel_eff_curve_cpy.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
    rel_eff_curve_cpy.rel_eff_eqn_order = 2;
    rel_eff_curve_cpy.phys_model_use_hoerl = false;
    options_cpy.rel_eff_curves.push_back( rel_eff_curve );
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy );
    options_cpy.check_same_hoerl_and_external_shielding_specifications();

    cerr << "Failed to throw an error for one rel_eff_curve not being a physical model" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = true;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = true;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy = rel_eff_curve;
    rel_eff_curve_cpy.shielded_by_other_phys_model_curve_shieldings.insert( 1 );
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy );

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy2 = rel_eff_curve;
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy2 );

    options_cpy.check_same_hoerl_and_external_shielding_specifications();
    cerr << "Failed to throw an error for same_external_shielding_for_all_rel_eff_curves with one curve shielding another" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = true;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = false;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy = rel_eff_curve;
    rel_eff_curve_cpy.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield );
    rel_eff_curve_cpy.phys_model_external_atten.push_back( make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield ) );

    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy );

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy2 = rel_eff_curve;
    rel_eff_curve_cpy2.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield );
    rel_eff_curve_cpy2.shielded_by_other_phys_model_curve_shieldings.insert( 0 );
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy2 );

    options_cpy.check_same_hoerl_and_external_shielding_specifications();

    // We are suppoest to get here
  }catch( std::exception &e )
  {
    cerr << "Eroneously threw for shielded_by_other_phys_model_curve_shieldings" << endl;
    assert( 0 );
  }

  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = true;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = false;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy = rel_eff_curve;
    rel_eff_curve_cpy.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield );
    rel_eff_curve_cpy.phys_model_external_atten.push_back( make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield ) );
    rel_eff_curve_cpy.shielded_by_other_phys_model_curve_shieldings.insert( 1 );

    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy );

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy2 = rel_eff_curve;
    rel_eff_curve_cpy2.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield );
    rel_eff_curve_cpy2.shielded_by_other_phys_model_curve_shieldings.insert( 0 );
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy2 );

    options_cpy.check_same_hoerl_and_external_shielding_specifications();

    cerr << "Didnt throw for cyclical curves shielding eachother" << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  try
  {
    RelActCalcAuto::Options options_cpy = options;
    options_cpy.same_hoerl_for_all_rel_eff_curves = true;
    options_cpy.same_external_shielding_for_all_rel_eff_curves = false;

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy = rel_eff_curve;
    rel_eff_curve_cpy.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield );
    rel_eff_curve_cpy.phys_model_external_atten.push_back( make_shared<RelActCalc::PhysicalModelShieldInput>( *ext_shield ) );
    rel_eff_curve_cpy.shielded_by_other_phys_model_curve_shieldings.insert( 1 );

    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy );

    RelActCalcAuto::RelEffCurveInput rel_eff_curve_cpy2 = rel_eff_curve;
    options_cpy.rel_eff_curves.push_back( rel_eff_curve_cpy2 );

    options_cpy.check_same_hoerl_and_external_shielding_specifications();

    cerr << "Didnt throw for shielding curve not having any shieldings defined." << endl;
    assert( 0 );
  }catch( std::exception &e )
  {
    // We are suppoest to get here
  }

  cout << "All auto hoerl and ext shield checks passed" << endl;
}//void check_auto_hoerl_and_ext_shield_checks()


void utile_ana()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const string data_dir = InterSpec::staticDataDirectory();
  assert( matdb );
  
  const string specfilename = "HEU-025_DU-100_Detective.Chn";
  SpecUtils::SpecFile specfile;
  const bool loaded_spec = specfile.load_file(specfilename, SpecUtils::ParserType::Auto );
  if( !loaded_spec )
  {
    cerr << "Failed to load '" << specfilename << "', aborting." << endl;
  }
  
  const string backfilename = "bg_Detective.Chn";
  SpecUtils::SpecFile backfile;
  const bool loaded_back = backfile.load_file(backfilename, SpecUtils::ParserType::Auto );
  if( !loaded_back )
  {
    cerr << "Failed to load '" << backfilename << "', aborting." << endl;
  }
  
  
  assert( specfile.num_measurements() == 1 );
  assert( backfile.num_measurements() == 1 );
  
  const shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement_at_index(size_t(0));
  assert( foreground );
  
  const shared_ptr<const SpecUtils::Measurement> background = backfile.measurement_at_index(size_t(0));
  assert( background );
  
  
  const string setup_xml_path = "multienrich_u_detective_releff.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );
    
  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );
    
  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );
    
  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node );
    
  assert( state.options.rel_eff_curves.size() == 2 );
  
  const string detector_xml_path = "Detective-DX_LANL_25cm_fwhm.drf.xml";
  rapidxml::file<char> detector_input_file( detector_xml_path.c_str() );
    
  rapidxml::xml_document<char> detector_doc;
  detector_doc.parse<rapidxml::parse_trim_whitespace>( detector_input_file.data() );
    
  const rapidxml::xml_node<char> *detector_base_node = detector_doc.first_node( "DetectorPeakResponse" );
  assert( detector_base_node );
    
    
  auto det = make_shared<DetectorPeakResponse>();
  det->fromXml( detector_base_node );
    
  vector<shared_ptr<const PeakDef>> all_peaks{};
    
  const SandiaDecay::Nuclide * const u232 = db->nuclide("U232");
  const SandiaDecay::Nuclide * const u234 = db->nuclide("U234");
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  assert( u235 && u238 && u232 && u234 );
    
        
  const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint leu_constraint
    = RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.002/0.98 ); //DU - 0.2% enriched

  state.options.rel_eff_curves[0].act_ratio_constraints.push_back( leu_constraint );
    
    
  //const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint heu_constraint
  //  = RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.82467/0.046746 ); //a mas ratio of about 35 give ~90% enrichment
  //state.options.rel_eff_curves[1].act_ratio_constraints.push_back( heu_constraint );
    
  //state.options.rel_eff_curves.resize(1);
  const bool use_hoerl = false;
  state.options.rel_eff_curves[0].phys_model_use_hoerl = use_hoerl;
  state.options.rel_eff_curves[1].phys_model_use_hoerl = use_hoerl;
  state.options.same_hoerl_for_all_rel_eff_curves = use_hoerl;
  //state.options.same_external_shielding_for_all_rel_eff_curves = true;
    
    
  const RelActCalcAuto::RelActAutoSolution solution
                = RelActCalcAuto::solve( state.options, foreground, background, det, all_peaks,
                                         PeakFitUtils::CoarseResolutionType::High, nullptr );
  ofstream out_html( specfilename + "_releff_result.html" );

  solution.print_summary( cout );
  solution.print_html_report( out_html );

  for( size_t i = 0; i < solution.m_rel_activities.size(); ++i )
  {
    pair<double,optional<double>> enrich = solution.mass_enrichment_fraction( u235, i );
    const double enrichment = enrich.first;
    const double u235_counts = solution.nuclide_counts( u235, i );
    const double u238_counts = solution.nuclide_counts( u238, i );
    cout << "Enrichment " << i << std::left << ": "
    << setprecision(6) << setw(11) << enrichment
    << " +- " << setw(11) << (enrich.second.has_value() ? enrich.second.value() : -999.0)
      << ", counts(u235)=" << setw(11) << u235_counts
      << ", counts(u238)=" << setw(11) << u238_counts << endl;
  }
}//void utile_ana()


void leu_heu_ana()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const string data_dir = InterSpec::staticDataDirectory();
  assert( matdb );

  const string specfilename = "NBS900+U295_B95Dish_2h_COAX.CNF";  //LEU + HEU
  //const string specfilename = "NBS900_B95Dish_2h_COAX.CNF"; //HEU
  //const string specfilename = "U295_B95Dish_2h_COAX.CNF"; //LEU
  SpecUtils::SpecFile specfile;
  const bool loaded_spec = specfile.load_file(specfilename, SpecUtils::ParserType::Auto );
  if( !loaded_spec )
  {
    cerr << "Failed to load '" << specfilename << "', aborting." << endl;
  }

  assert( specfile.num_measurements() == 1 );

  const shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement_at_index(size_t(0));
  assert( foreground );

  const shared_ptr<const SpecUtils::Measurement> background;

  const string setup_xml_path = "isotopics_by_nuclides_jozsef_releff.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );

  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );

  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );

  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node );

  assert( state.options.rel_eff_curves.size() == 2 );

  const string detector_xml_path = "LAB14_100cm_ISOCS (39%).drf.xml";
  rapidxml::file<char> detector_input_file( detector_xml_path.c_str() );

  rapidxml::xml_document<char> detector_doc;
  detector_doc.parse<rapidxml::parse_trim_whitespace>( detector_input_file.data() );

  const rapidxml::xml_node<char> *detector_base_node = detector_doc.first_node( "DetectorPeakResponse" );
  assert( detector_base_node );


  auto det = make_shared<DetectorPeakResponse>();
  det->fromXml( detector_base_node );

  vector<shared_ptr<const PeakDef>> all_peaks{};

  const SandiaDecay::Nuclide * const u232 = db->nuclide("U232");
  const SandiaDecay::Nuclide * const u234 = db->nuclide("U234");
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  assert( u235 && u238 && u232 && u234 );


  //const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint leu_constraint
  //  = RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.02764 );

  //state.options.rel_eff_curves[0].act_ratio_constraints.push_back( leu_constraint );


  //const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint heu_constraint
  //  = RelActCalcAuto::RelEffCurveInput::ActRatioConstraint::from_mass_ratio( u235, u238, 0.82467/0.046746 ); //a mas ratio of about 35 give ~90% enrichment
  //state.options.rel_eff_curves[1].act_ratio_constraints.push_back( heu_constraint );

  //state.options.rel_eff_curves.resize(1);
  //const bool use_hoerl = false;
  //state.options.rel_eff_curves[0].phys_model_use_hoerl = use_hoerl;
  //state.options.rel_eff_curves[1].phys_model_use_hoerl = use_hoerl;
  //state.options.same_hoerl_for_all_rel_eff_curves = use_hoerl;
  //state.options.same_external_shielding_for_all_rel_eff_curves = true;


  const RelActCalcAuto::RelActAutoSolution solution
                = RelActCalcAuto::solve( state.options, foreground, background, det, all_peaks,
                                         PeakFitUtils::CoarseResolutionType::High, nullptr );
  ofstream out_html( specfilename + "_releff_result.html" );

  solution.print_summary( cout );
  solution.print_html_report( out_html );

  for( size_t i = 0; i < solution.m_rel_activities.size(); ++i )
  {
    pair<double,optional<double>> enrich = solution.mass_enrichment_fraction( u235, i );
    const double enrichment = enrich.first;
    const double u235_counts = solution.nuclide_counts( u235, i );
    const double u238_counts = solution.nuclide_counts( u238, i );
    cout << "Enrichment " << i << std::left << ": "
    << setprecision(6) << setw(11) << enrichment
    << " +- " << setw(11) << (enrich.second.has_value() ? enrich.second.value() : -999.0)
    << ", counts(u235)=" << setw(11) << u235_counts
    << ", counts(u238)=" << setw(11) << u238_counts << endl;
  }
}//void leu_heu_ana()


/** Reproducer for the RelActAuto free-Fe-shielding χ² bug.

 Reads the natural-U spectrum `unat.n42` (contains foreground, background, and the DRF used
 in the GUI) and the `HPGe U (120-1001 keV)` preset, and performs two fits:

   Run A: preset as-is (external Fe "fit" enabled).  Currently produces > 1.2 % enrichment
          and χ²/dof ≈ 411/241 despite Fe converging near zero - this is the symptom.
   Run B: same preset with Fe forced to 0 mm and not fit.  Produces χ²/dof ≈ 186/243.

 The free-parameter fit (Run A) cannot legitimately have worse χ² than Run B, since Run B's
 solution is in Run A's feasible set.  Run A is therefore terminating above the true minimum,
 and this function prints enough information to start bisecting why.
*/
void unat_hpge_fe_example()
{
  using namespace std;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );

  const string data_dir = InterSpec::staticDataDirectory();
  const string specfile_path = "/Users/wcjohns/coding/InterSpec_master/unat.n42";
  const string preset_xml_path = SpecUtils::append_path( data_dir, "rel_act/HPGe U (120-1001 keV).xml" );

  // -- Load spectrum --
  auto specmeas = make_shared<SpecMeas>();
  if( !specmeas->load_file( specfile_path, SpecUtils::ParserType::Auto ) )
  {
    cerr << "unat_hpge_fe_example: failed to load '" << specfile_path << "' - aborting." << endl;
    return;
  }

  shared_ptr<const SpecUtils::Measurement> foreground;
  shared_ptr<const SpecUtils::Measurement> background;
  for( const shared_ptr<const SpecUtils::Measurement> &m : specmeas->measurements() )
  {
    if( !m )
      continue;
    if( m->source_type() == SpecUtils::SourceType::Background )
      background = m;
    else if( !foreground )
      foreground = m;
  }
  assert( foreground );
  if( !foreground )
  {
    cerr << "unat_hpge_fe_example: no foreground measurement in '" << specfile_path << "' - aborting." << endl;
    return;
  }

  shared_ptr<const DetectorPeakResponse> det = specmeas->detector();
  assert( det );
  if( !det )
  {
    cerr << "unat_hpge_fe_example: no detector in '" << specfile_path << "' - aborting." << endl;
    return;
  }

  // -- Load preset XML --
  // The file object owns the buffer that rapidxml parses in-place, so it must outlive
  // every access to setup_doc.
  std::unique_ptr<rapidxml::file<char>> setup_input_file;
  rapidxml::xml_document<char> setup_doc;
  try
  {
    setup_input_file = std::make_unique<rapidxml::file<char>>( preset_xml_path.c_str() );
    setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file->data() );
  }catch( const std::exception &e )
  {
    cerr << "unat_hpge_fe_example: failed to parse '" << preset_xml_path << "': " << e.what() << endl;
    return;
  }

  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );
  if( !setup_base_node )
  {
    cerr << "unat_hpge_fe_example: preset XML missing RelActCalcAuto root node - aborting." << endl;
    return;
  }

  RelActCalcAuto::RelActAutoGuiState state_A;
  try
  {
    state_A.deSerialize( setup_base_node );
  }catch( const std::exception &e )
  {
    cerr << "unat_hpge_fe_example: deSerialize failed: " << e.what() << endl;
    return;
  }

  // Deep copy for Run B (the Options struct owns shared_ptr shields; we'll make a new one).
  RelActCalcAuto::RelActAutoGuiState state_B = state_A;

  // Helper to print the key numbers from a solution.
  const auto summarize = []( const char *label, const RelActCalcAuto::RelActAutoSolution &sol,
                             const SandiaDecay::Nuclide *u235 )
  {
    cout << "\n==== " << label << " ====" << endl;
    cout << "  Status:               " << static_cast<int>(sol.m_status) << endl;
    if( !sol.m_error_message.empty() )
      cout << "  Error message:        " << sol.m_error_message << endl;
    cout << "  chi2:                 " << sol.m_chi2 << endl;
    cout << "  dof:                  " << sol.m_dof << endl;
    if( sol.m_dof > 0 )
      cout << "  chi2/dof:             " << (sol.m_chi2 / static_cast<double>(sol.m_dof)) << endl;
    cout << "  num function evals:   " << sol.m_num_function_eval_solution << endl;

    for( size_t re_i = 0; re_i < sol.m_rel_activities.size(); ++re_i )
    {
      pair<double,std::optional<double>> enrich = sol.mass_enrichment_fraction( u235, re_i );
      cout << "  U-235 enrichment[" << re_i << "]: " << (100.0*enrich.first) << " %";
      if( enrich.second.has_value() )
        cout << "  (± " << (100.0*(*enrich.second)) << " %)";
      cout << endl;
    }

    for( size_t re_i = 0; re_i < sol.m_phys_model_results.size(); ++re_i )
    {
      const auto &phys = sol.m_phys_model_results[re_i];
      if( !phys.has_value() )
        continue;

      if( phys->self_atten.has_value() )
      {
        const auto &s = *phys->self_atten;
        const double ad_gpcm2 = s.areal_density / PhysicalUnits::g_per_cm2;
        cout << "  self-atten AD:        " << ad_gpcm2 << " g/cm^2";
        if( s.material )
          cout << "  (" << (s.areal_density / s.material->density) / PhysicalUnits::mm << " mm " << s.material->name << ")";
        cout << endl;
      }
      for( size_t ext_i = 0; ext_i < phys->ext_shields.size(); ++ext_i )
      {
        const auto &e = phys->ext_shields[ext_i];
        const double ad_gpcm2 = e.areal_density / PhysicalUnits::g_per_cm2;
        cout << "  ext-atten[" << ext_i << "] AD:      " << ad_gpcm2 << " g/cm^2";
        if( e.material )
          cout << "  (" << (e.areal_density / e.material->density) / PhysicalUnits::mm << " mm " << e.material->name << ")";
        cout << "  [fit=" << (e.areal_density_was_fit ? "yes" : "no") << "]" << endl;
      }
      if( phys->hoerl_b.has_value() || phys->hoerl_c.has_value() )
      {
        cout << "  Hoerl b=" << (phys->hoerl_b.has_value() ? *phys->hoerl_b : 0.0)
             << ", c=" << (phys->hoerl_c.has_value() ? *phys->hoerl_c : 1.0) << endl;
      }
    }
    cout.flush();
  };

  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  assert( u235 );

  // --- Run A: preset as-supplied ---
  cout << "\n================================================================" << endl;
  cout << "Run A: HPGe U (120-1001 keV) preset as-is (Fe external fit enabled)" << endl;
  cout << "================================================================" << endl;
  const RelActCalcAuto::RelActAutoSolution sol_A = RelActCalcAuto::solve(
      state_A.options, foreground, background, det, {},
      PeakFitUtils::CoarseResolutionType::High, nullptr );
  summarize( "Run A result", sol_A, u235 );

  pair<double,std::optional<double>> enrich_A = sol_A.mass_enrichment_fraction( u235, 0 );
  const double enrich_A_frac = enrich_A.first;
  const double chi2_A = sol_A.m_chi2;
  const size_t dof_A = sol_A.m_dof;

  // --- Run B: force external Fe to 0 and not-fit ---
  // Walk the external attens and, for each, rebuild as zero-thickness, not-fit shield.
  if( !state_B.options.rel_eff_curves.empty() )
  {
    auto &curve = state_B.options.rel_eff_curves[0];
    for( shared_ptr<const RelActCalc::PhysicalModelShieldInput> &ext : curve.phys_model_external_atten )
    {
      if( !ext )
        continue;
      auto modified = make_shared<RelActCalc::PhysicalModelShieldInput>( *ext );
      modified->areal_density = 0.0;
      modified->fit_areal_density = false;
      modified->lower_fit_areal_density = 0.0;
      modified->upper_fit_areal_density = 0.0;
      ext = modified;
    }
  }

  cout << "\n================================================================" << endl;
  cout << "Run B: same preset, external shielding forced to 0 and not-fit" << endl;
  cout << "================================================================" << endl;
  const RelActCalcAuto::RelActAutoSolution sol_B = RelActCalcAuto::solve(
      state_B.options, foreground, background, det, {},
      PeakFitUtils::CoarseResolutionType::High, nullptr );
  summarize( "Run B result", sol_B, u235 );

  const double enrich_B_frac = sol_B.mass_enrichment_fraction( u235, 0 ).first;
  const double chi2_B = sol_B.m_chi2;
  const size_t dof_B = sol_B.m_dof;

  cout << "\n================================================================" << endl;
  cout << "Summary" << endl;
  cout << "================================================================" << endl;
  cout << "Run A (Fe fit):    enrichment = " << (100.0*enrich_A_frac) << " %"
       << ",  chi2/dof = " << chi2_A << "/" << dof_A
       << " = " << (dof_A > 0 ? chi2_A/dof_A : 0.0) << endl;
  cout << "Run B (Fe fixed):  enrichment = " << (100.0*enrich_B_frac) << " %"
       << ",  chi2/dof = " << chi2_B << "/" << dof_B
       << " = " << (dof_B > 0 ? chi2_B/dof_B : 0.0) << endl;

  // Symptom asserts - these document the current (buggy) behavior so we know the reproducer is live.
  // User observed ~1.87% enrichment in the GUI for Run A; threshold at 1.2% gives headroom for minor variance.
  if( enrich_A_frac <= 0.012 )
    cerr << "WARNING: Run A enrichment <= 1.2% - symptom may not be reproducing." << endl;
  assert( enrich_A_frac > 0.012 );

  // Run B chi2 should be notably lower than Run A chi2 - the whole point of the reproducer.
  if( chi2_B >= 0.6 * chi2_A )
    cerr << "WARNING: Run B chi2 is not notably better than Run A - bug may not be reproducing." << endl;
  assert( chi2_B < 0.6 * chi2_A );

  cout << "\nReproducer verified: Run A has higher chi2 than Run B by factor of "
       << (chi2_A / std::max( 1.0, chi2_B )) << endl;
}//void unat_hpge_fe_example()



int dev_code()
{
  // IDB enrichment comparison
  {
    const std::string docroot = SpecUtils::append_path( InterSpec::staticDataDirectory(), ".." );
    const std::string idb_u_dir   = SpecUtils::append_path( docroot, "IDB_U" );
    const std::string idb_u_meta  = SpecUtils::append_path( idb_u_dir, "IDB_U_metadata.xml" );
    const std::string idb_pu_dir  = SpecUtils::append_path( docroot, "IDB_Pu" );
    const std::string idb_pu_meta = SpecUtils::append_path( idb_pu_dir, "IDB_Pu_metadata.xml" );
    run_idb_comparison( idb_u_dir, idb_u_meta, IdbMaterialType::Uranium,    "IDB_U_comparison.html" );
    run_idb_comparison( idb_pu_dir, idb_pu_meta, IdbMaterialType::Plutonium, "IDB_Pu_comparison.html" );
    return 0;
  }

  unat_hpge_fe_example();
  return 1;

  check_auto_nuclide_constraints_checks();
  check_manual_nuclide_constraints_checks();
  check_auto_hoerl_and_ext_shield_checks();
  return 1;
  
  //check_auto_state_xml_serialization();

  run_u02_example();
  return 1;
  
  //check_physical_model_eff_function();
  //return 1;
  
  //example_manual_phys_model();
  //return 1;

  //return RelActCalcAuto::run_test();
  
  //czt_pu_example();
  //return 1;
  
  //run_multi_enrich_u02_ex();
  //return 1;

  //run_multi_enrich_sensitivity_study_u02();
  //return 1;

  //utile_ana();
  leu_heu_ana();
  return 1;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
  {
    cerr << "Failed to open SandiaDecayDataBase" << endl;
    return EXIT_FAILURE;
  }
  
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  
  const string data_dir = InterSpec::staticDataDirectory();
  assert( matdb );
  
  RelActCalcAuto::Options options;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NonLinearFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_4;
  options.spectrum_title = "Dev Title";
  options.skew_type = PeakDef::SkewType::NoSkew;
  options.additional_br_uncert = 0.01;
  
  
  RelActCalcAuto::RelEffCurveInput rel_eff_curve;
  rel_eff_curve.nucs_of_el_same_age = true;
  rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  rel_eff_curve.rel_eff_eqn_order = 0;
  rel_eff_curve.pu242_correlation_method = RelActCalc::PuCorrMethod::ByPu239Only;
  
  RelActCalc::PhysicalModelShieldInput self_atten_def;
  self_atten_def.atomic_number = 0.0;
  const std::shared_ptr<const Material> pu = matdb->material("Pu (plutonium)");
  assert( pu );
  self_atten_def.material = pu;
  self_atten_def.areal_density = 19.8 * PhysicalUnits::g_per_cm2;
  //self_atten_def.fit_atomic_number = false;
  //self_atten_def.lower_fit_atomic_number = 1.0;
  //self_atten_def.upper_fit_atomic_number = 98.0;
  self_atten_def.fit_areal_density = true;
  self_atten_def.lower_fit_areal_density = 0.0;
  self_atten_def.upper_fit_areal_density = 500.0 * PhysicalUnits::g_per_cm2;
  
  rel_eff_curve.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( self_atten_def );
  
  
  RelActCalc::PhysicalModelShieldInput ext_shield;
  ext_shield.atomic_number = 26;
  ext_shield.material = nullptr;
  ext_shield.areal_density = 1.0 * PhysicalUnits::g_per_cm2;
  ext_shield.fit_atomic_number = false;
  //ext_shield.lower_fit_atomic_number = 1.0;
  //ext_shield.upper_fit_atomic_number = 98.0;
  ext_shield.fit_areal_density = true;
  ext_shield.lower_fit_areal_density = 0.0;
  ext_shield.upper_fit_areal_density = 500.0 * PhysicalUnits::g_per_cm2;
  rel_eff_curve.phys_model_external_atten.push_back( make_shared<RelActCalc::PhysicalModelShieldInput>( ext_shield ) );
  
  
  
  // Using example spectra 522 from https://nds.iaea.org/idb/ - randomly selected.
  const string specfilename = "IDB-all-spectra/named_spectrum_files/spec522_239Pu_84.8097.spe";
  SpecUtils::SpecFile specfile;
  const bool loaded_spec = specfile.load_file(specfilename, SpecUtils::ParserType::Auto );
  if( !loaded_spec )
  {
    cerr << "Failed to load '" << specfilename << "', aborting." << endl;
    return EXIT_FAILURE;
  }
  
  if( specfile.num_measurements() != 1 )
  {
    cerr << "Expected file to have one meas." << endl;
    return EXIT_FAILURE;
  }

  // Use "HPGe Pu 140 to 240 keV", adapted from FRAM
  vector<RelActCalcAuto::RoiRange> energy_ranges;

  auto make_roi = []( double lower, double upper, PeakContinuum::OffsetType cont_type ) {
    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = lower;
    roi.upper_energy = upper;
    roi.continuum_type = cont_type;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    return roi;
  };

  energy_ranges.push_back( make_roi( 127.9, 131.8, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 124.2, 127.8, PeakContinuum::LinearStep ) );
  energy_ranges.push_back( make_roi( 145.1, 150.9, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 150.9, 153.8, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 163.2, 166.2, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 158.9, 163.2, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 223.1, 230.4, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 200.9, 204.9, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 204.9, 211.0, PeakContinuum::LinearStep ) );
  energy_ranges.push_back( make_roi( 253.1, 257.1, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 309.5, 314.3, PeakContinuum::Linear ) );
  energy_ranges.push_back( make_roi( 338.9, 348.6, PeakContinuum::LinearStep ) );
  energy_ranges.push_back( make_roi( 330.2, 338.9, PeakContinuum::LinearStep ) );
  energy_ranges.push_back( make_roi( 363.2, 386.4, PeakContinuum::FlatStep ) );
  
  
  const vector<RelActCalcAuto::NucInputInfo> nuclides{ {
      db->nuclide("Pu238"),
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      std::nullopt, //fit_age_min
      std::nullopt, //fit_age_max
      std::nullopt, //min_rel_act
      std::nullopt, //max_rel_act
      std::nullopt, //starting_rel_act
      {}, //gammas to exclude
      "rgb(0, 0, 255)",
    }, {
      db->nuclide("Pu239"),
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      std::nullopt, //fit_age_min
      std::nullopt, //fit_age_max
      std::nullopt, //min_rel_act
      std::nullopt, //max_rel_act
      std::nullopt, //starting_rel_act
      {}, //gammas to exclude
      "rgb(255, 69, 0)",
    }, {
      db->nuclide("Pu240"),
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      std::nullopt, //fit_age_min
      std::nullopt, //fit_age_max
      std::nullopt, //min_rel_act
      std::nullopt, //max_rel_act
      std::nullopt, //starting_rel_act
      {}, //gammas to exclude
      "rgb(34, 139, 34)",
    }, {
      db->nuclide("Pu241"),
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      std::nullopt, //fit_age_min
      std::nullopt, //fit_age_max
      std::nullopt, //min_rel_act
      std::nullopt, //max_rel_act
      std::nullopt, //starting_rel_act
      {}, //gammas to exclude
      "rgb(204, 204, 0)",
    }
  };
  
  options.floating_peaks = vector<RelActCalcAuto::FloatingPeak>{};
  rel_eff_curve.nuclides = nuclides;
  options.rois = energy_ranges;
  
  options.rel_eff_curves.push_back( rel_eff_curve );

  rel_eff_curve.phys_model_external_atten.clear();
  self_atten_def.areal_density = 0.1 * PhysicalUnits::g_per_cm2;
  rel_eff_curve.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( self_atten_def );
  options.rel_eff_curves.push_back( rel_eff_curve );
  
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement_at_index( size_t(0) );
  assert( foreground );
  shared_ptr<const SpecUtils::Measurement> background = nullptr;
  shared_ptr<const DetectorPeakResponse> drf = nullptr; //We dont have this
  vector<shared_ptr<const PeakDef>> all_peaks{};
  
  const RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
    options, foreground, background, drf, all_peaks,
    PeakFitUtils::CoarseResolutionType::High, nullptr );
  
  ofstream out_html( "result.html" );
  solution.print_summary( cout );
  solution.print_html_report( out_html );

  return 1;
}//int dev_code()


namespace
{
  // Convert XML isotope name "234U" → SandiaDecay format "U234"
  std::string xml_to_sandia_nuc( const std::string &xml )
  {
    size_t i = 0;
    while( i < xml.size() && std::isdigit( static_cast<unsigned char>( xml[i] ) ) )
      ++i;
    return xml.substr( i ) + xml.substr( 0, i );
  }


  // Resolve a config file: try as-is, then under {staticDataDir}/rel_act/
  std::string resolve_rel_eff_config( const std::string &value )
  {
    if( value.empty() )
      throw std::runtime_error( "Empty RelActConfigFile value" );
    if( SpecUtils::is_file( value ) )
      return value;
    const std::string rel_act_dir = SpecUtils::append_path( InterSpec::staticDataDirectory(), "rel_act" );
    const std::string full = SpecUtils::append_path( rel_act_dir, value );
    if( SpecUtils::is_file( full ) )
      return full;
    throw std::runtime_error( "Cannot find config file '" + value + "'" );
  }


  // Load DRF from a tag value: try as file path, then as built-in name
  std::shared_ptr<DetectorPeakResponse> load_drf_from_tag( const std::string &value )
  {
    if( value.empty() )
      return nullptr;
    try
    {
      return BatchActivity::init_drf_from_name( value, "" );
    }
    catch( ... ) {}

    try
    {
      return BatchActivity::init_drf_from_name( "", value );
    }
    catch( ... ) {}

    std::cerr << "Warning: could not load DRF from tag value '" << value << "'\n";
    return nullptr;
  }


  std::string extract_head_content( const std::string &html )
  {
    const size_t h1 = html.find( "<head>" );
    const size_t h2 = html.find( "</head>" );
    if( h1 == std::string::npos || h2 == std::string::npos || h2 < h1 )
      return "";
    return html.substr( h1 + 6, h2 - h1 - 6 );
  }


  std::string extract_body_content( const std::string &html )
  {
    const size_t b1 = html.find( "<body>" );
    const size_t b2 = html.rfind( "</body>" );
    if( b1 == std::string::npos || b2 == std::string::npos || b2 < b1 )
      return html;
    return html.substr( b1 + 6, b2 - b1 - 6 );
  }


  // Replace the two hardcoded rel-eff chart IDs and wrap the init script in an IIFE
  // so variables don't leak into global scope when multiple reports are on one page.
  void fix_html_for_embedding( std::string &body, const int spec_id )
  {
    const std::string prefix = "s" + std::to_string( spec_id ) + "_";

    SpecUtils::ireplace_all( body, "\"chartsdiv\"",   ( "\"" + prefix + "chartsdiv\""  ).c_str() );
    SpecUtils::ireplace_all( body, "\"releffchart\"", ( "\"" + prefix + "releffchart\"" ).c_str() );

    // Wrap the RelEffPlot init <script> block in an IIFE to scope its variables
    const std::string marker = "new RelEffPlot(";
    const size_t marker_pos = body.find( marker );
    if( marker_pos == std::string::npos )
      return;

    const size_t script_open = body.rfind( "<script>", marker_pos );
    if( script_open == std::string::npos )
      return;

    const size_t content_start = script_open + 8;  // skip "<script>"
    const size_t script_close = body.find( "</script>", content_start );
    if( script_close == std::string::npos )
      return;

    const std::string old_content = body.substr( content_start, script_close - content_start );
    const std::string new_content = "\n(function(){\n" + old_content + "})();\n";
    body.replace( content_start, script_close - content_start, new_content );
  }


  std::string pct_diff_bg_color( const double pct_diff )
  {
    const double abs_diff = std::abs( pct_diff );
    if( abs_diff < 5.0 )  return "#c8e6c9"; // green
    if( abs_diff < 15.0 ) return "#fff9c4"; // yellow
    return "#ffccbc";                        // red
  }


  // Format a non-negative double; returns "&mdash;" when v < 0 (sentinel for "not available")
  std::string fmt( const double v, const int decimals = 4 )
  {
    if( v < 0.0 )
      return "&mdash;";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision( decimals ) << v;
    return oss.str();
  }

  // Format a signed double (e.g. percent difference — can legitimately be negative)
  std::string fmt_signed( const double v, const int decimals = 2 )
  {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision( decimals ) << v;
    return oss.str();
  }
}//anonymous namespace


void run_idb_comparison( const std::string &idb_dir,
                         const std::string &metadata_xml_path,
                         IdbMaterialType material_type,
                         const std::string &output_html_path,
                         size_t max_spectra )
{
  const bool is_uranium   = (material_type == IdbMaterialType::Uranium);
  const bool is_plutonium = (material_type == IdbMaterialType::Plutonium);
  const std::string material_name = is_uranium ? "Uranium" : (is_plutonium ? "Plutonium" : "MOX");
  const std::string primary_nuc   = is_uranium ? "U235" : "Pu240";

  const std::string default_config_name = is_uranium
      ? "HPGe U (120-1001 keV).xml"
      : "HPGe Pu (120-780 keV).xml";

  const std::string default_config_path = resolve_rel_eff_config( default_config_name );

  // Ordered list of nuclides for table columns
  const std::vector<std::string> nuclide_cols = is_uranium
      ? std::vector<std::string>{ "U234", "U235", "U236", "U238" }
      : std::vector<std::string>{ "Pu238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241" };

  // ---- Parse metadata XML ----
  std::vector<char> xml_data;
  try
  {
    SpecUtils::load_file_data( metadata_xml_path.c_str(), xml_data );
  }
  catch( std::exception &e )
  {
    throw std::runtime_error( "Failed to read metadata XML '" + metadata_xml_path + "': " + e.what() );
  }

  rapidxml::xml_document<char> doc;
  try
  {
    doc.parse<rapidxml::parse_trim_whitespace>( &xml_data[0] );
  }
  catch( std::exception &e )
  {
    throw std::runtime_error( "Failed to parse metadata XML '" + metadata_xml_path + "': " + e.what() );
  }

  const rapidxml::xml_node<char> *root = doc.first_node();
  if( !root )
    throw std::runtime_error( "Empty metadata XML" );

  // ---- Collect HPGe entries (up to max_spectra) ----
  struct SpectrumEntry
  {
    int id;
    std::string filename;
    std::string det_eff_tag;      // from <DetectorEfficiency>
    std::string config_file_tag;  // from <RelActConfigFile>
    std::map<std::string, double> known_mass_fracs;  // SandiaDecay format → wt%
  };

  std::vector<SpectrumEntry> entries;

  for( const rapidxml::xml_node<char> *sn = root->first_node( "spectrum" );
       sn; sn = sn->next_sibling( "spectrum" ) )
  {
    if( entries.size() >= max_spectra )
      break;

    const rapidxml::xml_attribute<char> *id_attr   = sn->first_attribute( "id" );
    const rapidxml::xml_attribute<char> *file_attr = sn->first_attribute( "file" );
    if( !id_attr || !file_attr )
      continue;

    const rapidxml::xml_node<char> *det_type_node = sn->first_node( "Detector_type" );
    if( !det_type_node || !SpecUtils::iequals_ascii( det_type_node->value(), "HPGe" ) )
      continue;

    const rapidxml::xml_node<char> *fracs_node =
        sn->first_node( "decay_corrected_mass_fractions_wt_pct" );
    if( !fracs_node )
      continue;

    SpectrumEntry entry;
    entry.id       = std::stoi( id_attr->value() );
    entry.filename = file_attr->value();

    for( const rapidxml::xml_node<char> *iso = fracs_node->first_node( "isotope" );
         iso; iso = iso->next_sibling( "isotope" ) )
    {
      const rapidxml::xml_attribute<char> *name_attr = iso->first_attribute( "name" );
      if( !name_attr || !iso->value() )
        continue;
      try
      {
        entry.known_mass_fracs[xml_to_sandia_nuc( name_attr->value() )] =
            std::stod( iso->value() );
      }
      catch( ... ) {}
    }

    const rapidxml::xml_node<char> *det_eff_node    = sn->first_node( "DetectorEfficiency" );
    const rapidxml::xml_node<char> *config_file_node = sn->first_node( "RelActConfigFile" );
    if( det_eff_node    && det_eff_node->value() )    entry.det_eff_tag      = det_eff_node->value();
    if( config_file_node && config_file_node->value() ) entry.config_file_tag = config_file_node->value();

    entries.push_back( std::move( entry ) );
  }//for each <spectrum>

  std::cout << "Found " << entries.size() << " HPGe spectra to process in "
            << metadata_xml_path << "\n";

  // ---- Process each spectrum ----
  struct SpectrumResult
  {
    SpectrumEntry entry;
    std::string   drf_name;
    std::string   rel_eff_xml_used;
    std::string   error_msg;
    std::shared_ptr<RelActCalcAuto::RelActAutoSolution> solution;
    std::map<std::string, std::pair<double, double>> fit_fracs;  // nuclide → (wt%, unc_wt%)
    double fit_age_s     = -1.0;
    double fit_age_unc_s = -1.0;
  };

  std::vector<SpectrumResult> results;
  results.reserve( entries.size() );

  for( const SpectrumEntry &entry : entries )
  {
    SpectrumResult result;
    result.entry = entry;

    std::cout << "  ID " << entry.id << ": " << entry.filename << " ... " << std::flush;

    try
    {
      // Resolve config file and load state
      const std::string config_path = entry.config_file_tag.empty()
          ? default_config_path
          : resolve_rel_eff_config( entry.config_file_tag );
      result.rel_eff_xml_used = config_path;

      const std::shared_ptr<RelActCalcAuto::RelActAutoGuiState> state =
          BatchRelActAuto::load_state_from_xml_file( config_path );

      RelActCalcAuto::Options options = state->options;
      options.foreground_filename = SpecUtils::filename( entry.filename );

      // Load DRF if tag present
      std::shared_ptr<DetectorPeakResponse> drf;
      if( !entry.det_eff_tag.empty() )
        drf = load_drf_from_tag( entry.det_eff_tag );
      if( drf )
        result.drf_name = drf->name();

      // Load the spectrum file
      const std::string spec_path = SpecUtils::append_path( idb_dir, entry.filename );
      SpecMeas spec_meas;
      if( !spec_meas.load_file( spec_path, SpecUtils::ParserType::Auto, spec_path ) )
        throw std::runtime_error( "Could not load spectrum file '" + spec_path + "'" );

      const std::shared_ptr<const SpecUtils::Measurement> foreground =
          spec_meas.sum_measurements( spec_meas.sample_numbers(),
                                      spec_meas.detector_names(), nullptr );
      if( !foreground || foreground->num_gamma_channels() < 16 )
        throw std::runtime_error( "Spectrum has too few gamma channels" );

      // Solve - this dev path filters to Detector_type=="HPGe" above (line ~2555),
      // so the coarse resolution is High.
      RelActCalcAuto::RelActAutoSolution sol = RelActCalcAuto::solve(
          options, foreground, nullptr, drf, {}, PeakFitUtils::CoarseResolutionType::High );

      if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        throw std::runtime_error( "Solve failed: " + sol.m_error_message );

      result.solution = std::make_shared<RelActCalcAuto::RelActAutoSolution>( std::move( sol ) );

      // Extract fit enrichment fractions
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      for( const auto &kv : entry.known_mass_fracs )
      {
        const SandiaDecay::Nuclide *nuc = db->nuclide( kv.first );
        if( !nuc )
          continue;
        try
        {
          const std::pair<double, std::optional<double>> frac_unc =
              result.solution->mass_enrichment_fraction( nuc, 0 );
          const double unc = frac_unc.second.has_value() ? frac_unc.second.value() : -1.0;
          result.fit_fracs[kv.first] = { frac_unc.first * 100.0, unc >= 0.0 ? unc * 100.0 : -1.0 };
        }
        catch( ... ) {}
      }

      // Find fit age (first nuclide with age_was_fit)
      if( !result.solution->m_rel_activities.empty() )
      {
        for( const RelActCalcAuto::NuclideRelAct &nra : result.solution->m_rel_activities[0] )
        {
          if( nra.age_was_fit && nra.age > 0.0 )
          {
            result.fit_age_s     = nra.age;
            result.fit_age_unc_s = nra.age_uncertainty;
            break;
          }
        }
      }

      std::cout << "OK\n";
    }
    catch( std::exception &e )
    {
      result.error_msg = e.what();
      std::cout << "FAILED: " << e.what() << "\n";
    }

    results.push_back( std::move( result ) );
  }//for each entry

  // ---- Generate HTML ----
  inja::Environment inja_env = RelActAutoReport::get_default_inja_env( "" );

  // Render per-spectrum HTML and extract head/body parts
  std::string shared_head;
  struct SectionHtml
  {
    int         id;
    std::string body;  // empty on render failure
  };
  std::vector<SectionHtml> sections;
  sections.reserve( results.size() );

  for( const SpectrumResult &result : results )
  {
    SectionHtml sec;
    sec.id = result.entry.id;

    if( result.solution )
    {
      try
      {
        const nlohmann::json data = RelActAutoReport::solution_to_json( *result.solution );
        std::string full_html = RelActAutoReport::render_template( inja_env, data, "html", "" );

        if( shared_head.empty() )
        {
          shared_head = extract_head_content( full_html );
          // Update the title for the combined page
          const size_t t1 = shared_head.find( "<title>" );
          const size_t t2 = shared_head.find( "</title>" );
          if( t1 != std::string::npos && t2 != std::string::npos && t2 > t1 )
          {
            shared_head = shared_head.substr( 0, t1 + 7 )
                        + "IDB " + material_name + " Enrichment Comparison"
                        + shared_head.substr( t2 );
          }
        }

        std::string body = extract_body_content( full_html );
        fix_html_for_embedding( body, result.entry.id );
        sec.body = std::move( body );
      }
      catch( std::exception &e )
      {
        sec.body = "<p><em>HTML rendering failed: " + std::string( e.what() ) + "</em></p>";
      }
    }

    sections.push_back( std::move( sec ) );
  }

  // Write combined HTML
  std::ofstream html_out( output_html_path );
  if( !html_out.is_open() )
    throw std::runtime_error( "Cannot open output file '" + output_html_path + "'" );

  html_out << "<!DOCTYPE html>\n<html>\n<head>\n";
  if( !shared_head.empty() )
  {
    html_out << shared_head << "\n";
  }
  else
  {
    html_out << "<meta charset=\"UTF-8\">\n"
             << "<title>IDB " << material_name << " Enrichment Comparison</title>\n";
  }
  // Additional CSS for the comparison tables
  html_out << "<style>\n"
           << "  .cmp-table { border-collapse: collapse; margin: 10px 0; font-size: 0.9em; }\n"
           << "  .cmp-table th, .cmp-table td { padding: 4px 8px; border: 1px solid #aaa; text-align: right; }\n"
           << "  .cmp-table th { background: #e0e0e0; text-align: center; }\n"
           << "  .cmp-table td.label { text-align: left; }\n"
           << "  .spec-section { margin: 20px 0; border: 1px solid #ccc; padding: 10px; }\n"
           << "  .spec-section h2 { margin-top: 0; font-size: 1.1em; }\n"
           << "  summary { font-size: 1.05em; font-weight: bold; cursor: pointer; }\n"
           << "</style>\n";
  html_out << "</head>\n<body>\n";

  html_out << "<h1>IDB " << material_name << " Enrichment Comparison</h1>\n";
  html_out << "<p>Processed <b>" << results.size() << "</b> HPGe spectra. "
           << "Default config: <code>" << default_config_name << "</code></p>\n";

  // ---- Summary table ----
  html_out << "<h2>Summary</h2>\n"
           << "<table class=\"cmp-table\">\n<thead>\n<tr>\n"
           << "  <th>ID</th><th>File</th>\n";
  for( const std::string &nuc : nuclide_cols )
    html_out << "  <th>Known " << nuc << " %</th>\n";
  html_out << "  <th>Fit " << primary_nuc << " %</th>\n"
           << "  <th>±1σ</th>\n"
           << "  <th>" << primary_nuc << " % diff</th>\n"
           << "  <th>Status</th>\n"
           << "</tr>\n</thead>\n<tbody>\n";

  for( const SpectrumResult &result : results )
  {
    const std::string link = "#spec_" + std::to_string( result.entry.id );

    if( !result.error_msg.empty() )
    {
      const int ncols = static_cast<int>( nuclide_cols.size() ) + 4;
      html_out << "<tr style=\"background:#e0e0e0\">\n"
               << "  <td><a href=\"" << link << "\">" << result.entry.id << "</a></td>\n"
               << "  <td>" << result.entry.filename << "</td>\n"
               << "  <td colspan=\"" << ncols << "\">" << result.error_msg << "</td>\n"
               << "</tr>\n";
      continue;
    }

    // Known fractions
    std::string known_primary = "";
    double known_primary_val = -1.0;
    {
      const auto it = result.entry.known_mass_fracs.find( primary_nuc );
      if( it != result.entry.known_mass_fracs.end() )
      {
        known_primary_val = it->second;
        known_primary = fmt( known_primary_val, 4 );
      }
    }

    // Fit primary fraction
    double fit_primary_val = -1.0;
    double fit_primary_unc  = -1.0;
    {
      const auto it = result.fit_fracs.find( primary_nuc );
      if( it != result.fit_fracs.end() )
      {
        fit_primary_val = it->second.first;
        fit_primary_unc = it->second.second;
      }
    }

    double pct_diff = 0.0;
    bool   have_diff = false;
    if( fit_primary_val >= 0.0 && known_primary_val > 0.0 )
    {
      pct_diff  = (fit_primary_val - known_primary_val) / known_primary_val * 100.0;
      have_diff = true;
    }

    const std::string bg = have_diff ? pct_diff_bg_color( pct_diff ) : "#e0e0e0";

    html_out << "<tr style=\"background:" << bg << "\">\n"
             << "  <td><a href=\"" << link << "\">" << result.entry.id << "</a></td>\n"
             << "  <td>" << result.entry.filename << "</td>\n";

    for( const std::string &nuc : nuclide_cols )
    {
      const auto it = result.entry.known_mass_fracs.find( nuc );
      if( it != result.entry.known_mass_fracs.end() )
        html_out << "  <td>" << fmt( it->second, 4 ) << "</td>\n";
      else
        html_out << "  <td>&mdash;</td>\n";
    }

    html_out << "  <td>" << fmt( fit_primary_val, 4 ) << "</td>\n"
             << "  <td>" << fmt( fit_primary_unc,  4 ) << "</td>\n"
             << "  <td>" << (have_diff ? fmt_signed( pct_diff, 2 ) : "&mdash;") << "</td>\n"
             << "  <td style=\"color:green\">OK</td>\n"
             << "</tr>\n";
  }//for each result

  html_out << "</tbody>\n</table>\n";

  // ---- Per-spectrum detail sections ----
  html_out << "<h2>Spectrum Details</h2>\n";

  for( size_t i = 0; i < results.size(); ++i )
  {
    const SpectrumResult &result = results[i];
    const SectionHtml    &sec    = sections[i];

    html_out << "<div id=\"spec_" << result.entry.id << "\" class=\"spec-section\">\n"
             << "<details>\n"
             << "<summary>ID " << result.entry.id << ": " << result.entry.filename;

    // Show primary nuclide summary in the details label
    {
      const auto known_it = result.entry.known_mass_fracs.find( primary_nuc );
      const auto fit_it   = result.fit_fracs.find( primary_nuc );
      if( known_it != result.entry.known_mass_fracs.end() )
        html_out << " &mdash; known " << primary_nuc << "=" << fmt( known_it->second, 4 ) << "%";
      if( fit_it != result.fit_fracs.end() && fit_it->second.first >= 0.0 )
        html_out << ", fit=" << fmt( fit_it->second.first, 4 ) << "%";
    }
    html_out << "</summary>\n";

    // Header info
    html_out << "<p>\n"
             << "  <b>DRF:</b> " << (result.drf_name.empty() ? "<em>none</em>" : result.drf_name)
             << " &nbsp; <b>RelEff config:</b> <code>"
             << SpecUtils::filename( result.rel_eff_xml_used ) << "</code>\n"
             << "</p>\n";

    if( !result.error_msg.empty() )
    {
      html_out << "<p style=\"color:red\"><b>Error:</b> " << result.error_msg << "</p>\n";
    }
    else
    {
      // Mass fraction comparison table
      html_out << "<table class=\"cmp-table\">\n"
               << "<thead><tr>"
               << "<th>Nuclide</th><th>Known wt%</th><th>Fit wt%</th>"
               << "<th>±1σ wt%</th><th>% Diff</th>"
               << "</tr></thead>\n<tbody>\n";

      for( const std::string &nuc : nuclide_cols )
      {
        const auto known_it = result.entry.known_mass_fracs.find( nuc );
        const auto fit_it   = result.fit_fracs.find( nuc );

        const double known_val = (known_it != result.entry.known_mass_fracs.end())
            ? known_it->second : -1.0;
        const double fit_val   = (fit_it != result.fit_fracs.end())
            ? fit_it->second.first : -1.0;
        const double fit_unc   = (fit_it != result.fit_fracs.end())
            ? fit_it->second.second : -1.0;

        double  row_pct_diff  = 0.0;
        bool    have_row_diff = false;
        if( fit_val >= 0.0 && known_val > 0.0 )
        {
          row_pct_diff  = (fit_val - known_val) / known_val * 100.0;
          have_row_diff = true;
        }

        const std::string row_bg = have_row_diff ? pct_diff_bg_color( row_pct_diff ) : "";
        const std::string style_attr = row_bg.empty() ? "" : " style=\"background:" + row_bg + "\"";

        html_out << "<tr" << style_attr << ">"
                 << "<td class=\"label\">" << nuc << "</td>"
                 << "<td>" << fmt( known_val, 4 ) << "</td>"
                 << "<td>" << fmt( fit_val,   4 ) << "</td>"
                 << "<td>" << fmt( fit_unc,   4 ) << "</td>"
                 << "<td>" << (have_row_diff ? fmt_signed( row_pct_diff, 2 ) : "&mdash;") << "</td>"
                 << "</tr>\n";
      }//for each nuclide column

      // Fit age row (Pu/MOX only, or any time age was fit)
      if( result.fit_age_s > 0.0 )
      {
        const std::string age_str = PhysicalUnits::printToBestTimeUnits( result.fit_age_s );
        const std::string unc_str = (result.fit_age_unc_s > 0.0)
            ? PhysicalUnits::printToBestTimeUnits( result.fit_age_unc_s )
            : "&mdash;";
        html_out << "<tr><td class=\"label\" colspan=\"2\">Fit age</td>"
                 << "<td colspan=\"2\">" << age_str << "</td>"
                 << "<td>±" << unc_str << "</td></tr>\n";
      }

      html_out << "</tbody></table>\n";

      // Embedded per-spectrum report (charts + tables from render_template)
      if( !sec.body.empty() )
        html_out << sec.body << "\n";
    }//if no error

    html_out << "</details>\n</div>\n\n";
  }//for each result section

  html_out << "</body>\n</html>\n";

  std::cout << "\nWrote comparison report to: " << output_html_path << "\n";
}//void run_idb_comparison(...)

}//namespace RelActAutoDev
