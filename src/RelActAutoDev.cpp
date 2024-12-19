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

#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActAutoDev.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"

#include "SpecUtils/SpecFile.h"

using namespace std;

namespace RelActAutoDev
{
int dev_code()
{
  // Using example spectra from https://nds.iaea.org/idb/
  cout << "dev_code()" << endl;
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
  {
    cerr << "Failed to open SandiaDecayDataBase" << endl;
    return EXIT_FAILURE;
  }
  
  RelActCalcAuto::Options options;
  options.fit_energy_cal = true;
  options.nucs_of_el_same_age = true;
  options.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramEmpirical;
  options.rel_eff_eqn_order = 4;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_4;
  options.spectrum_title = "Dev Title";
  options.pu242_correlation_method = RelActCalc::PuCorrMethod::ByPu239Only;
  options.skew_type = PeakDef::SkewType::NoSkew;
  
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

  // "HPGe Pu 140 to 240 keV"
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
