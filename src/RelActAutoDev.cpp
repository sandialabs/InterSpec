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

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActAutoDev.h"
#include "InterSpec/RelActCalcAuto.h"
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
    return EXIT_FAILURE;
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
  
  
int dev_code()
{
  //check_physical_model_eff_function();
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
  
  RelActCalcAuto::Options::PhysicalModelShieldInput self_atten_def;
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
  
  options.phys_model_self_atten = self_atten_def;
  
  
  RelActCalcAuto::Options::PhysicalModelShieldInput ext_shield;
  ext_shield.atomic_number = 26;
  ext_shield.material = nullptr;
  ext_shield.areal_density = 1.0 * PhysicalUnits::g_per_cm2;
  ext_shield.fit_atomic_number = false;
  //ext_shield.lower_fit_atomic_number = 1.0;
  //ext_shield.upper_fit_atomic_number = 98.0;
  ext_shield.fit_areal_density = true;
  ext_shield.lower_fit_areal_density = 0.0;
  ext_shield.upper_fit_areal_density = 500.0 * PhysicalUnits::g_per_cm2;
  options.phys_model_external_atten.push_back( ext_shield );
  
  
  
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
