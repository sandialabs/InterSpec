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
#include <random>
#include <string>
#include <vector>
#include <functional>

#define BOOST_TEST_MODULE RelActCalcAuto_FwhmForm_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"

using namespace std;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for the Rel. Act. tool is required for this test." );

namespace
{
  // Log-spaced energies (keV) over [lower, upper].
  vector<double> log_energies( const double lower, const double upper, const size_t num )
  {
    vector<double> energies;
    for( size_t i = 0; i < num; ++i )
      energies.push_back( lower * std::pow( upper/lower, i/(num - 1.0) ) );
    return energies;
  }
}//namespace


BOOST_AUTO_TEST_CASE( noise_plus_curved_power_basics )
{
  using RelActCalcAuto::FwhmForm;

  BOOST_CHECK_EQUAL( RelActCalcAuto::num_parameters( FwhmForm::NoisePlusCurvedPower ), size_t(6) );
  const string name = RelActCalcAuto::to_str( FwhmForm::NoisePlusCurvedPower );
  BOOST_CHECK( RelActCalcAuto::fwhm_form_from_str( name.c_str() ) == FwhmForm::NoisePlusCurvedPower );

  // At 661 keV the power-law term is exactly w1, whatever the slopes.
  const vector<double> pars{ 0.8, 1.2, 0.2, 0.7, 50.0, 3000.0 };
  BOOST_CHECK_CLOSE( RelActCalcAuto::eval_fwhm( 661.0, FwhmForm::NoisePlusCurvedPower, pars ),
                     std::sqrt( 0.8*0.8 + 1.2*1.2 ), 1.0E-9 );
}//BOOST_AUTO_TEST_CASE( noise_plus_curved_power_basics )


BOOST_AUTO_TEST_CASE( noise_plus_curved_power_is_positive_and_non_decreasing )
{
  using RelActCalcAuto::FwhmForm;

  // Any parameters inside the solve's bounds give a positive curve that never falls - inside the
  //  fitted range and on the end-slope continuation beyond it.
  std::mt19937 rng( 1234 );
  std::uniform_real_distribution<double> noise( 0.0, 20.0 ), width( 0.05, 100.0 ), slope( 0.0, 1.5 );
  const vector<double> energies = log_energies( 5.0, 10000.0, 400 );
  for( int trial = 0; trial < 500; ++trial )
  {
    const vector<double> pars{ noise(rng), width(rng), slope(rng), slope(rng), 40.0, 3000.0 };
    double previous = 0.0;
    for( const double energy : energies )
    {
      const double fwhm = RelActCalcAuto::eval_fwhm( energy, FwhmForm::NoisePlusCurvedPower, pars );
      BOOST_REQUIRE( std::isfinite(fwhm) && (fwhm > 0.0) );
      BOOST_REQUIRE_GE( fwhm, previous*(1.0 - 1.0E-12) );
      previous = fwhm;
    }
  }
}//BOOST_AUTO_TEST_CASE( noise_plus_curved_power_is_positive_and_non_decreasing )


BOOST_AUTO_TEST_CASE( noise_plus_curved_power_fit_recovers_its_own_curve )
{
  using RelActCalcAuto::FwhmForm;

  const vector<vector<double>> truths{
    { 0.9, 1.3, 0.35, 0.60, 40.0, 3000.0 },   // HPGe-like
    { 0.0, 48.0, 0.65, 0.55, 40.0, 3000.0 },  // NaI-like
    { 8.0, 13.0, 0.10, 0.45, 40.0, 3000.0 }   // CZT-like, noise dominated low
  };
  const vector<double> energies = log_energies( 40.0, 3000.0, 30 );
  for( const vector<double> &truth : truths )
  {
    vector<double> fwhms;
    for( const double energy : energies )
      fwhms.push_back( RelActCalcAuto::eval_fwhm( energy, FwhmForm::NoisePlusCurvedPower, truth ) );

    const vector<double> fit = RelActCalcAuto::fit_noise_plus_curved_power( energies, fwhms, 40.0, 3000.0 );
    BOOST_REQUIRE_EQUAL( fit.size(), size_t(6) );
    for( size_t i = 0; i < energies.size(); ++i )
      BOOST_CHECK_CLOSE( RelActCalcAuto::eval_fwhm( energies[i], FwhmForm::NoisePlusCurvedPower, fit ), fwhms[i], 0.1 );
  }
}//BOOST_AUTO_TEST_CASE( noise_plus_curved_power_fit_recovers_its_own_curve )


BOOST_AUTO_TEST_CASE( noise_plus_curved_power_follows_detector_class_curves )
{
  using RelActCalcAuto::FwhmForm;

  // The class resolution curves the fitter uses as its shape priors, 40 keV - 3 MeV, with the worst
  //  relative error allowed.  The HPGe curve is GADRAS's, whose low-energy term makes it dip from
  //  1.58 keV at 40 keV to 1.51 keV at 300 keV before rising - a dip no non-decreasing form follows.
  struct ClassCurve
  {
    string name;
    std::function<float(float)> fwhm;
    double max_rel_err;
  };
  const vector<ClassCurve> curves{
    { "HPGe", &PeakFitUtils::hpge_fwhm_fcn, 0.08 },
    { "NaI", &PeakFitUtils::nai_fwhm_fcn, 0.04 },
    { "LaBr3", &PeakFitUtils::labr_fwhm_fcn, 0.02 },
    { "CZT", &PeakFitUtils::czt_fwhm_fcn, 0.03 }
  };
  const vector<double> energies = log_energies( 40.0, 3000.0, 40 );
  for( const ClassCurve &curve : curves )
  {
    vector<double> fwhms;
    for( const double energy : energies )
      fwhms.push_back( curve.fwhm( static_cast<float>(energy) ) );

    const vector<double> fit = RelActCalcAuto::fit_noise_plus_curved_power( energies, fwhms, 40.0, 3000.0 );
    double max_rel_err = 0.0;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double model = RelActCalcAuto::eval_fwhm( energies[i], FwhmForm::NoisePlusCurvedPower, fit );
      max_rel_err = std::max( max_rel_err, std::fabs( model/fwhms[i] - 1.0 ) );
    }
    BOOST_TEST_MESSAGE( curve.name << ": worst relative error " << max_rel_err );
    BOOST_CHECK_MESSAGE( max_rel_err < curve.max_rel_err, curve.name << " class curve: worst relative error " << max_rel_err );
  }
}//BOOST_AUTO_TEST_CASE( noise_plus_curved_power_follows_detector_class_curves )
