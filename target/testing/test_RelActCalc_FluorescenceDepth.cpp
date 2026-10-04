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

/* `RelActCalc::fluorescence_self_atten_correction(...)` and the exponential integrals it uses, against references
 computed with scipy (special.expn/exp1, and adaptive quadrature of the slab production-escape integrals), and its
 Ceres-Jet derivative against finite differences.
 */

#include "InterSpec_config.h"

#include <cmath>
#include <vector>

#define BOOST_TEST_MODULE RelActCalc_FluorescenceDepth_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalc_imp.hpp"

using namespace std;


BOOST_AUTO_TEST_CASE( exponential_integrals )
{
  // x, E1, E2, E3 (scipy.special)
  const double ref[][4] = {
    { 1.0E-4, 8.63322470457, 0.999036682529, 0.499900050666 },
    { 0.3, 0.905676651676, 0.469115225179, 0.300041826564 },
    { 1.0, 0.219383934396, 0.148495506776, 0.109691967198 },
    { 1.7, 0.0746546444013, 0.0557706285706, 0.0439367277414 },
    { 5.0, 0.00114829559128, 0.000996469042709, 0.000877800892771 },
    { 20.0, 9.83552529065E-11, 9.40485643086E-11, 9.00911681335E-11 } };

  for( const auto &r : ref )
  {
    BOOST_CHECK_CLOSE( RelActCalc::expint_e1( r[0] ), r[1], 1.0E-7 );
    BOOST_CHECK_CLOSE( RelActCalc::expint_e2( r[0] ), r[2], 1.0E-7 );
    BOOST_CHECK_CLOSE( RelActCalc::expint_e3( r[0] ), r[3], 1.0E-7 );
  }

  BOOST_CHECK_EQUAL( RelActCalc::expint_e2( 0.0 ), 1.0 );
  BOOST_CHECK_EQUAL( RelActCalc::expint_e3( 0.0 ), 0.5 );
}//exponential_integrals


BOOST_AUTO_TEST_CASE( slab_correction_values )
{
  struct Case { double mu_f; vector<double> mu, w; double A, expected; };
  // mu per unit areal density; e.g. U3O8: 1.75 (U K-alpha1), 1.30 (U K-beta1), 1.33 (185.7 keV) cm2/g
  const vector<Case> cases = {
    { 1.75, {1.33}, {1.0}, 20.0, 0.8347969202 },
    { 1.30, {1.33}, {1.0}, 20.0, 0.8650293824 },
    { 1.75, {1.33}, {1.0}, 1.0, 0.9881364841 },
    { 1.20, {3.81, 0.09}, {1.0, 0.5}, 5.0, 0.9555307804 },
    { 1.92, {1.33, 2.22, 0.46}, {0.7, 0.2, 0.1}, 0.3, 0.9989791731 } };

  for( const Case &c : cases )
  {
    const double val = RelActCalc::fluorescence_self_atten_correction( c.mu_f, c.mu, c.w, c.A );
    BOOST_CHECK_CLOSE( val, c.expected, 0.01 );
  }

  // Thick U3O8 excited by 185.7 keV: K-beta1/K-alpha1 raised by ~3.6 %, as the micro-calorimeter U K pattern shows
  const double ratio = RelActCalc::fluorescence_self_atten_correction( 1.30, vector<double>{1.33}, vector<double>{1.0}, 20.0 )
                       / RelActCalc::fluorescence_self_atten_correction( 1.75, vector<double>{1.33}, vector<double>{1.0}, 20.0 );
  BOOST_CHECK_CLOSE( ratio, 0.8650293824/0.8347969202, 0.01 );

  // No layer, or nothing to excite: no correction
  BOOST_CHECK_EQUAL( RelActCalc::fluorescence_self_atten_correction( 1.75, vector<double>{1.33}, vector<double>{1.0}, 0.0 ), 1.0 );
  BOOST_CHECK_EQUAL( RelActCalc::fluorescence_self_atten_correction( 1.75, vector<double>{}, vector<double>{}, 5.0 ), 1.0 );
}//slab_correction_values


BOOST_AUTO_TEST_CASE( slab_correction_jet_derivative )
{
  typedef ceres::Jet<double,2> Jet;   // derivatives w.r.t. areal density and the first exciting weight
  const vector<double> mu = { 1.33, 0.46 };
  for( const double A : { 0.3, 2.0, 10.0 } )
  {
    const Jet areal_density( A, 0 );
    const vector<Jet> w = { Jet( 0.7, 1 ), Jet( 0.3 ) };
    const Jet val = RelActCalc::fluorescence_self_atten_correction( 1.30, mu, w, areal_density );

    const auto f = [&]( const double a, const double w0 ) {
      return RelActCalc::fluorescence_self_atten_correction( 1.30, mu, vector<double>{ w0, 0.3 }, a );
    };
    const double h = 1.0E-6;
    BOOST_CHECK_CLOSE( val.a, f( A, 0.7 ), 1.0E-9 );
    BOOST_CHECK_CLOSE( val.v[0], (f( A + h, 0.7 ) - f( A - h, 0.7 )) / (2*h), 1.0E-3 );
    BOOST_CHECK_CLOSE( val.v[1], (f( A, 0.7 + h ) - f( A, 0.7 - h )) / (2*h), 1.0E-3 );
  }
}//slab_correction_jet_derivative
