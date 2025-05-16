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

#include <regex>
#include <memory>
#include <vector>
#include <fstream>
#include <iostream>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActAutoDev.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace RelActAutoDev
{

#if( PERFORM_DEVELOPER_CHECKS )
void check_auto_state_xml_serialization()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

  MaterialDB matdb;
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );

  RelActCalcAuto::RelActAutoGuiState state;

  state.options.fit_energy_cal = true;
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
  const Material * const pu_db = matdb.material("Pu (plutonium)");
  assert( pu_db );
  const shared_ptr<const Material> pu = std::make_shared<Material>( *pu_db );
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
  u235_input.nuclide = db->nuclide("U235");
  u235_input.age = 20.0*PhysicalUnits::year;
  u235_input.fit_age = false;
  u235_input.gammas_to_exclude = { 26.325 };
  u235_input.peak_color_css = "red";
  curve.nuclides.push_back( u235_input );

  RelActCalcAuto::NucInputInfo u238_input;
  u238_input.nuclide = db->nuclide("U238");
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
  range.force_full_range = false;
  range.allow_expand_for_peak_width = false;
  state.options.rois.push_back( range );

  // Add a second range
  range.lower_energy = 50.0;
  range.upper_energy = 200.0;
  range.continuum_type = PeakContinuum::OffsetType::Linear;
  range.force_full_range = true;
  range.allow_expand_for_peak_width = true;
  state.options.rois.push_back( range );

  RelActCalcAuto::FloatingPeak peak;
  peak.energy = 1000.0;
  peak.release_fwhm = true;
  peak.apply_energy_cal_correction = true;
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
    state2.deSerialize( rel_act_node, &matdb );
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
  co60_input.nuclide = db->nuclide("Co60");
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
    state3.deSerialize( rel_act_node2, &matdb );
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
  
  MaterialDB matdb;
  
  const string data_dir = InterSpec::staticDataDirectory();
  const string materialfile = SpecUtils::append_path( data_dir, "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );
  
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
  
  const Material * const pu_db = matdb.material("Pu (plutonium)");
  assert( pu_db );
  const shared_ptr<const Material> pu = std::make_shared<Material>( *pu_db );
  
  RelActCalc::PhysModelShield<double> self_atten;
  self_atten.material = pu;
  self_atten.areal_density = 1.98*PhysicalUnits::g_per_cm2;
  
  
  RelActCalc::PhysModelShield<double> ext_atten;
  ext_atten.material = std::make_shared<Material>( *matdb.material("stainless-steel NIST") );
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

  MaterialDB matdb;
  
  const string data_dir = InterSpec::staticDataDirectory();
  const string materialfile = SpecUtils::append_path( data_dir, "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );
  
  const Material * const u_db = matdb.material("U (uranium)");
  assert( u_db );
  const shared_ptr<const Material> uranium = std::make_shared<Material>( *u_db );
  

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
  sol.print_html_report( out_html, "Manual Phys Model", meas, disp_peaks, nullptr, 0.0 );


  cout << "Solution: " << endl;
  sol.print_summary( cout );
  cout << endl;
}//void example_manual_phys_model()

void run_u02_example()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  MaterialDB matdb;
  const string data_dir = InterSpec::staticDataDirectory();
  const string materialfile = SpecUtils::append_path( data_dir, "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );
  
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
  shared_ptr<const SpecUtils::Measurement> background = specfile.measurement(size_t(0));
  assert( background );
  assert( background->title() == "Background" );
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement(size_t(37));
  assert( foreground );
  assert( SpecUtils::istarts_with( foreground->title(), "U02_50%_50% @ 25 cm H=100 cm" ) );
  
  
  const string setup_xml_path = "isotopics_by_nuclides_mixed_U02_sample-2.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );
  
  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );

  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );
  
  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node, &matdb );
  
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
                                              foreground, background, det, all_peaks, nullptr );
  
  const double end_cpu = SpecUtils::get_cpu_time();
  const double end_wall = SpecUtils::get_wall_time();
  
  ofstream out_html( "U02_result.html" );
  solution.print_summary( cout );
  solution.print_html_report( out_html );
  
  cout << "Enrichment: " << solution.mass_enrichment_fraction( db->nuclide("U235"), 0 ) << endl;
  
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
  MaterialDB matdb;
  const string data_dir = InterSpec::staticDataDirectory();
  const string materialfile = SpecUtils::append_path( data_dir, "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );
  
  const string specfilename = "czt_600keV__Pu84_60mm_1mmCd_ex.n42";
  
  SpecMeas specfile;
  const bool loaded_spec = specfile.load_file(specfilename, SpecUtils::ParserType::Auto );
  if( !loaded_spec )
  {
    cerr << "Failed to load '" << specfilename << "', aborting." << endl;
  }
  
  
  assert( specfile.num_measurements() == 1 );
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement( size_t(0) );
  assert( foreground );

  
  const string setup_xml_path = "isotopics_by_nuclides_CZT500_Pu84_@60mm_1mmCd_04-2.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );
  
  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );
  
  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );
  
  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node, &matdb );
  
  assert( state.options.rel_eff_curves.size() == 1 );

  
  auto det = specfile.detector();
  
  
  //state.options.rel_eff_curves.resize(1);
  state.options.rel_eff_curves[0].phys_model_use_hoerl = false;
  
  
  // We can constrain the RelActivity.
  for( auto &nuc : state.options.rel_eff_curves[0].nuclides )
  {
    if( nuc.nuclide->symbol == "Pu239" )
    {
      nuc.starting_rel_act = 750*16228.0;
      nuc.min_rel_act = nuc.max_rel_act = nuc.starting_rel_act;
    }else if( nuc.nuclide->symbol == "Am241" )
    {
      nuc.starting_rel_act = 5*96521.1;
      nuc.min_rel_act = nuc.max_rel_act = nuc.starting_rel_act;
      
    }
  }


  vector<shared_ptr<const PeakDef>> all_peaks{};
  const RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve( state.options,
                                                                            foreground, nullptr, det, all_peaks, nullptr );
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
  u235_input.nuclide = u235;
  u235_input.age = 20.0 * PhysicalUnits::year;
  u235_input.fit_age = false;
  u235_input.gammas_to_exclude = {};
  u235_input.peak_color_css = "rgb(0, 0, 255)";
  
  RelActCalcAuto::NucInputInfo u238_input;
  u238_input.nuclide = u238;
  u238_input.age = 20.0 * PhysicalUnits::year;
  u238_input.fit_age = false;
  u238_input.gammas_to_exclude = {};
  u238_input.peak_color_css = "rgb(255, 69, 0)";

  RelActCalcAuto::NucInputInfo u234_input;
  u234_input.nuclide = u234;
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
    co60_input.nuclide = db->nuclide("Co60");
    co60_input.age = 20.0 * PhysicalUnits::year;
    co60_input.fit_age = false;
    co60_input.gammas_to_exclude = {};
    rel_eff_cpy.nuclides.push_back( co60_input );


     RelActCalcAuto::NucInputInfo cs137_input;
    cs137_input.nuclide = db->nuclide("Cs137");
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
    nuc_constraint.constrained_source = co60_input.nuclide;
    nuc_constraint.controlling_source = u238;
    nuc_constraint.constrained_to_controlled_activity_ratio = 0.1;
    rel_eff_cpy.act_ratio_constraints.push_back( nuc_constraint );

    nuc_constraint.constrained_source = cs137_input.nuclide;
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

  cout << "All auto hoerl and ext shield checks passed" << endl;
}//void check_auto_hoerl_and_ext_shield_checks()


void utile_ana()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  MaterialDB matdb;
  const string data_dir = InterSpec::staticDataDirectory();
  const string materialfile = SpecUtils::append_path( data_dir, "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );
  
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
  
  const shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement(size_t(0));
  assert( foreground );
  
  const shared_ptr<const SpecUtils::Measurement> background = backfile.measurement(size_t(0));
  assert( background );
  
  
  const string setup_xml_path = "multienrich_u_detective_releff.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );
    
  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );
    
  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );
    
  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node, &matdb );
    
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
                = RelActCalcAuto::solve( state.options, foreground, background, det, all_peaks, nullptr );
  ofstream out_html( specfilename + "_releff_result.html" );
    
  solution.print_summary( cout );
  solution.print_html_report( out_html );
    
  for( size_t i = 0; i < solution.m_rel_activities.size(); ++i )
  {
    const double enrichment = solution.mass_enrichment_fraction( u235, i );
    const double u235_counts = solution.nuclide_counts( u235, i );
    const double u238_counts = solution.nuclide_counts( u238, i );
    cout << "Enrichment " << i << std::left << ": " << setprecision(6) << setw(11) << enrichment
      << ", counts(u235)=" << setw(11) << u235_counts
      << ", counts(u238)=" << setw(11) << u238_counts << endl;
  }
}//void utile_ana()


void leu_heu_ana()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  MaterialDB matdb;
  const string data_dir = InterSpec::staticDataDirectory();
  const string materialfile = SpecUtils::append_path( data_dir, "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );

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

  const shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement(size_t(0));
  assert( foreground );

  const shared_ptr<const SpecUtils::Measurement> background;

  const string setup_xml_path = "isotopics_by_nuclides_jozsef_releff.xml";
  rapidxml::file<char> setup_input_file( setup_xml_path.c_str() );

  rapidxml::xml_document<char> setup_doc;
  setup_doc.parse<rapidxml::parse_trim_whitespace>( setup_input_file.data() );

  const rapidxml::xml_node<char> *setup_base_node = setup_doc.first_node( "RelActCalcAuto" );
  assert( setup_base_node );

  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( setup_base_node, &matdb );

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
                = RelActCalcAuto::solve( state.options, foreground, background, det, all_peaks, nullptr );
  ofstream out_html( specfilename + "_releff_result.html" );

  solution.print_summary( cout );
  solution.print_html_report( out_html );

  for( size_t i = 0; i < solution.m_rel_activities.size(); ++i )
  {
    const double enrichment = solution.mass_enrichment_fraction( u235, i );
    const double u235_counts = solution.nuclide_counts( u235, i );
    const double u238_counts = solution.nuclide_counts( u238, i );
    cout << "Enrichment " << i << std::left << ": " << setprecision(6) << setw(11) << enrichment
      << ", counts(u235)=" << setw(11) << u235_counts
      << ", counts(u238)=" << setw(11) << u238_counts << endl;
  }
}//void leu_heu_ana()


int dev_code()
{
  check_auto_nuclide_constraints_checks();
  //check_manual_nuclide_constraints_checks();
  //check_auto_hoerl_and_ext_shield_checks();
  //return 1;
  
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
  
  MaterialDB matdb;
  
  const string data_dir = InterSpec::staticDataDirectory();
  const string materialfile = SpecUtils::append_path( data_dir, "MaterialDataBase.txt" );
  matdb.parseGadrasMaterialFile( materialfile, db, false );
  
  RelActCalcAuto::Options options;
  options.fit_energy_cal = true;
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
  const Material * const pu = matdb.material("Pu (plutonium)");
  assert( pu );
  self_atten_def.material = std::make_shared<const Material>( *pu );
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
  const vector<RelActCalcAuto::RoiRange> energy_ranges{ {
      127.9,                                 //Lower energy
      131.8,                                 //Upper energy
    PeakContinuum::Linear, //Continuum type
      true,                                  //force_full_range
      false                                  //allow_expand_for_peak_width
    },
    { 124.2, 127.8, PeakContinuum::LinearStep, true, false },
    { 145.1, 150.9, PeakContinuum::Linear, true, false },
    { 150.9, 153.8, PeakContinuum::Linear, true, false },
    { 163.2, 166.2, PeakContinuum::Linear, true, false },
    { 158.9, 163.2, PeakContinuum::Linear, true, false },
    { 223.1, 230.4, PeakContinuum::Linear, true, false },
    { 200.9, 204.9, PeakContinuum::Linear, true, false },
    { 204.9, 211.0, PeakContinuum::LinearStep, true, false },
    { 253.1, 257.1, PeakContinuum::Linear, true, false },
    { 309.5, 314.3, PeakContinuum::Linear, true, false },
    { 338.9, 348.6, PeakContinuum::LinearStep, true, false },
    { 330.2, 338.9, PeakContinuum::LinearStep, true, false },
    { 363.2, 386.4, PeakContinuum::FlatStep, true, false },
  };
  
  
  const vector<RelActCalcAuto::NucInputInfo> nuclides{ {
      db->nuclide("Pu238"),
      nullptr,
      nullptr,
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
      nullptr,
      nullptr,
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
      nullptr,
      nullptr,
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
      nullptr,
      nullptr,
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
  
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement( size_t(0) );
  assert( foreground );
  shared_ptr<const SpecUtils::Measurement> background = nullptr;
  shared_ptr<const DetectorPeakResponse> drf = nullptr; //We dont have this
  vector<shared_ptr<const PeakDef>> all_peaks{};
  
  const RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
    options, foreground, background, drf, all_peaks, nullptr );
  
  ofstream out_html( "result.html" );
  solution.print_summary( cout );
  solution.print_html_report( out_html );
  
  return 1;
}//int dev_code()
}//namespace RelActAutoDev
