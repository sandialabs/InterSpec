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

#include <memory>
#include <vector>
#include <fstream>
#include <iostream>

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
  
  const vector<shared_ptr<const Material>> external_attenuations{
    std::make_shared<Material>( *matdb.material("stainless-steel NIST") ),
    nullptr //generic material
  };
  
  
  const vector<double> paramaters{
    0.0,  //self-atten, atomic_number - must be zero since we are specifying a material
    1.98, //self-atten, areal density, in units of g/cm2 (ie, 1 mm Pu)
    0.0,  //stainless, AN - must be zero
    0.8,  //stainless, AD - e.g., 1 mm
    32.0, //generic AN - e.g., Germanium
    0.532, //generic AD - e.g., 1 mm of Ge
    1.0, // Modified Hoerl b: pow(0.001*energy,b)
    1.0  // Modified Hoerl c: pow(c,1000/energy)
  };
  
  function<double(double)> fcn = RelActCalc::physical_model_eff_function(
                            pu, external_attenuations, det, paramaters.data(), paramaters.size() );
  
  cout << "Energy (keV), Counts" << endl;
  for( double energy = 50.0; energy < 3000; energy += 2 )
    cout << energy << "," << fcn(energy) << endl;
  
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

  RelActCalc::PhysicalModelShieldInput self_atten_def;
  self_atten_def.atomic_number = 0.0;
  self_atten_def.material = uranium;
  self_atten_def.areal_density = 10.0*PhysicalUnits::g_per_cm2;
  self_atten_def.fit_atomic_number = false;
  self_atten_def.lower_fit_atomic_number = 1.0;
  self_atten_def.upper_fit_atomic_number = 98.0;
  self_atten_def.fit_areal_density = true;
  self_atten_def.lower_fit_areal_density = 0.0;
  self_atten_def.upper_fit_areal_density = 500*PhysicalUnits::g_per_cm2;

  rel_eff_solve_input.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( self_atten_def );

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
  rel_eff_solve_input.phys_model_external_attens.push_back( make_shared<RelActCalc::PhysicalModelShieldInput>( external_atten_def ) );


  RelActCalcManual::RelEffSolution sol = RelActCalcManual::solve_relative_efficiency( rel_eff_solve_input );

  if( sol.m_status != RelActCalcManual::ManualSolutionStatus::Success )
  {
    cerr << "Failed to solve relative efficiency" << endl;
    return;
  }

  ofstream out_html( "manual_phys_model_result.html" );
  std::vector<std::shared_ptr<const PeakDef>> disp_peaks( orig_peaks->begin(), orig_peaks->end() );
  sol.print_html_report( out_html, "Manual Phys Model", meas, disp_peaks, nullptr, 0.0 );


  cout << "Solution: " << endl;
  sol.print_summary( cout );
  cout << endl;
}//void example_manual_phys_model()

  
int dev_code()
{
  //check_physical_model_eff_function();
  //return 1;
  
  //example_manual_phys_model();
  //return 1;

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
  options.nucs_of_el_same_age = true;
  options.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  options.rel_eff_eqn_order = 0;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_4;
  options.spectrum_title = "Dev Title";
  options.pu242_correlation_method = RelActCalc::PuCorrMethod::ByPu239Only;
  options.skew_type = PeakDef::SkewType::NoSkew;
  
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
  
  options.phys_model_self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( self_atten_def );
  
  
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
  options.phys_model_external_atten.push_back( make_shared<RelActCalc::PhysicalModelShieldInput>( ext_shield ) );
  
  
  
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
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      {}, //gammas to exclude
      "rgb(0, 0, 255)",
    }, {
      db->nuclide("Pu239"),
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      {}, //gammas to exclude
      "rgb(255, 69, 0)",
    }, {
      db->nuclide("Pu240"),
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      {}, //gammas to exclude
      "rgb(34, 139, 34)",
    }, {
      db->nuclide("Pu241"),
      20.0*PhysicalUnits::year,  //Default age
      false, //fit age
      {}, //gammas to exclude
      "rgb(204, 204, 0)",
    }
  };
  
  
  const vector<RelActCalcAuto::FloatingPeak> extra_peaks{};
  
  
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurement( size_t(0) );
  assert( foreground );
  shared_ptr<const SpecUtils::Measurement> background = nullptr;
  shared_ptr<const DetectorPeakResponse> drf = nullptr; //We dont have this
  vector<shared_ptr<const PeakDef>> all_peaks{};
  
  const RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
    options, energy_ranges, nuclides, extra_peaks, foreground, background, drf, all_peaks, nullptr );
  
  ofstream out_html( "result.html" );
  solution.print_summary( cout );
  solution.print_html_report( out_html );
  
  return 1;
}//int dev_code()
}//namespace RelActAutoDev
