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
// Fix for Windows WinSock header ordering issue
// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE test_DetectorEfficiency_suite
#include <boost/test/included/unit_test.hpp>

#include <map>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
  bool close_enough( const double a, const double b, const double rel_tol = 1e-5 )
  {
    const double diff = fabs(a - b);
    const double max_val = (std::max)(fabs(a), fabs(b));
    return (diff <= rel_tol * max_val) || (diff <= 1e-10);
  }
}//namespace


BOOST_AUTO_TEST_CASE( test_point_uncerts_corr_limits )
{
  const vector<float> energies = { 60.0f, 120.0f, 400.0f, 1000.0f, 2000.0f };
  const vector<float> uncerts = { 0.10f, 0.07f, 0.05f, 0.04f, 0.06f };
  const size_t n = energies.size();

  // Very large correlation length: C[j][k] -> u_j * u_k (rank-1, fully correlated)
  {
    shared_ptr<DetectorEfficiencyUncert> uncert
                = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts, 1.0E6 );
    const vector<float> &cov = uncert->covarianceMatrix();
    BOOST_REQUIRE_EQUAL( cov.size(), n*n );

    for( size_t j = 0; j < n; ++j )
      for( size_t k = 0; k < n; ++k )
        BOOST_CHECK_MESSAGE( close_enough( cov[j*n + k], double(uncerts[j])*uncerts[k], 1e-3 ),
          "Fully-correlated C[" + to_string(j) + "][" + to_string(k) + "] should be u_j*u_k" );
  }

  // Non-positive correlation length: diagonal matrix
  {
    shared_ptr<DetectorEfficiencyUncert> uncert
                = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts, 0.0 );
    const vector<float> &cov = uncert->covarianceMatrix();
    BOOST_REQUIRE_EQUAL( cov.size(), n*n );

    for( size_t j = 0; j < n; ++j )
    {
      for( size_t k = 0; k < n; ++k )
      {
        if( j == k )
          BOOST_CHECK( close_enough( cov[j*n + k], double(uncerts[j])*uncerts[j], 1e-5 ) );
        else
          BOOST_CHECK_EQUAL( cov[j*n + k], 0.0f );
      }
    }
  }

  // Intermediate correlation length: off-diagonals match the Gaussian kernel
  {
    const double corr_len = 0.5;
    shared_ptr<DetectorEfficiencyUncert> uncert
                = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts, corr_len );
    const vector<float> &cov = uncert->covarianceMatrix();

    for( size_t j = 0; j < n; ++j )
    {
      for( size_t k = 0; k < n; ++k )
      {
        const double dlne = log(energies[j]) - log(energies[k]);
        const double rho = exp( -0.5 * pow(dlne/corr_len, 2.0) );
        const double expected = rho * uncerts[j] * uncerts[k];
        BOOST_CHECK_MESSAGE( close_enough( cov[j*n + k], expected, 1e-4 ),
          "Kernel covariance C[" + to_string(j) + "][" + to_string(k) + "]: got "
          + to_string(cov[j*n + k]) + ", expected " + to_string(expected) );
      }
    }
  }
}//test_point_uncerts_corr_limits


BOOST_AUTO_TEST_CASE( test_node_covariance_psd )
{
  // The covariance among arbitrary requested energies must remain positive
  //  semi-definite, including off-node energies and band contributions.
  std::mt19937 rng( 987654321u );
  std::uniform_real_distribution<double> energy_dist( 20.0, 3500.0 );
  std::uniform_real_distribution<double> uncert_dist( 0.005, 0.25 );
  std::uniform_real_distribution<double> corr_dist( 0.05, 3.0 );
  std::uniform_real_distribution<double> x_dist( -1.0, 1.0 );

  for( size_t trial = 0; trial < 25; ++trial )
  {
    const size_t nnode = 2 + (rng() % 12);
    vector<float> energies, uncerts;
    for( size_t i = 0; i < nnode; ++i )
    {
      energies.push_back( static_cast<float>(energy_dist(rng)) );
      uncerts.push_back( static_cast<float>(uncert_dist(rng)) );
    }

    shared_ptr<DetectorEfficiencyUncert> uncert
            = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts, corr_dist(rng) );

    // Add some piecewise bands too
    vector<EffUncertBand> bands;
    bands.push_back( EffUncertBand{ 0.0f, 100.0f, 0.05f } );
    bands.push_back( EffUncertBand{ 100.0f, 700.0f, 0.08f } );
    bands.push_back( EffUncertBand{ 800.0f, 3000.0f, 0.03f } );
    uncert->setBands( bands );

    // Random request energies, including beyond the node range
    const size_t nreq = 3 + (rng() % 10);
    vector<double> req_energies;
    for( size_t i = 0; i < nreq; ++i )
      req_energies.push_back( 5.0 + 4000.0 * (0.5 + 0.5*x_dist(rng)) );

    const vector<double> cov = uncert->efficiencyFracCovariance( req_energies );
    BOOST_REQUIRE_EQUAL( cov.size(), nreq*nreq );

    // Check x^T * M * x >= -1e-12 for random x
    for( size_t xtrial = 0; xtrial < 20; ++xtrial )
    {
      vector<double> x( nreq );
      for( size_t i = 0; i < nreq; ++i )
        x[i] = x_dist(rng);

      double xtmx = 0.0;
      for( size_t i = 0; i < nreq; ++i )
        for( size_t j = 0; j < nreq; ++j )
          xtmx += x[i] * cov[i*nreq + j] * x[j];

      BOOST_CHECK_MESSAGE( xtmx >= -1.0E-12,
        "Covariance not PSD: x^T*M*x = " + to_string(xtmx) );
    }//for( size_t xtrial = 0; xtrial < 20; ++xtrial )
  }//for( size_t trial = 0; trial < 25; ++trial )
}//test_node_covariance_psd


BOOST_AUTO_TEST_CASE( test_covariance_interpolation )
{
  const vector<float> energies = { 100.0f, 400.0f, 1600.0f };
  const vector<float> uncerts = { 0.10f, 0.05f, 0.08f };

  shared_ptr<DetectorEfficiencyUncert> uncert
              = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts, 0.7 );
  const vector<float> &node_cov = uncert->covarianceMatrix();
  const size_t n = energies.size();

  // Requesting the node energies returns the node covariance exactly
  {
    const vector<double> req = { 100.0, 400.0, 1600.0 };
    const vector<double> cov = uncert->nodeFracCovariance( req );
    for( size_t i = 0; i < n; ++i )
      for( size_t j = 0; j < n; ++j )
        BOOST_CHECK( close_enough( cov[i*n + j], node_cov[i*n + j], 1e-5 ) );
  }

  // The log-energy midpoint of two nodes blends them with weight 0.5
  {
    const double mid = sqrt( 100.0 * 400.0 );  //log-midpoint
    const vector<double> req = { mid };
    const vector<double> cov = uncert->nodeFracCovariance( req );
    BOOST_REQUIRE_EQUAL( cov.size(), 1u );

    // var = (0.5,0.5,0) * C * (0.5,0.5,0)^T
    const double expected = 0.25*node_cov[0] + 0.25*node_cov[1*n+1] + 0.5*node_cov[1];
    BOOST_CHECK_MESSAGE( close_enough( cov[0], expected, 1e-5 ),
      "Midpoint variance: got " + to_string(cov[0]) + ", expected " + to_string(expected) );
  }

  // Energies beyond the node range clamp to the end-node variance
  {
    const vector<double> req = { 10.0, 99999.0 };
    const vector<double> cov = uncert->nodeFracCovariance( req );
    BOOST_CHECK( close_enough( cov[0], node_cov[0], 1e-5 ) );          //below first node
    BOOST_CHECK( close_enough( cov[3], node_cov[(n-1)*n + (n-1)], 1e-5 ) ); //above last node
  }
}//test_covariance_interpolation


BOOST_AUTO_TEST_CASE( test_band_covariance_semantics )
{
  auto uncert = make_shared<DetectorEfficiencyUncert>();

  vector<EffUncertBand> bands;
  bands.push_back( EffUncertBand{ 0.0f, 50.0f, 0.05f } );
  bands.push_back( EffUncertBand{ 50.0f, 120.0f, 0.07f } );
  bands.push_back( EffUncertBand{ 120.0f, 3000.0f, 0.03f } );
  BOOST_REQUIRE_NO_THROW( uncert->setBands( bands ) );

  // Two energies in the same band: covariance u^2 everywhere in the block
  {
    const vector<double> req = { 60.0, 100.0 };
    const vector<double> cov = uncert->efficiencyFracCovariance( req );
    BOOST_CHECK( close_enough( cov[0], 0.07*0.07 ) );
    BOOST_CHECK( close_enough( cov[1], 0.07*0.07 ) );  //100% correlated
    BOOST_CHECK( close_enough( cov[2], 0.07*0.07 ) );
    BOOST_CHECK( close_enough( cov[3], 0.07*0.07 ) );
  }

  // Two energies in different bands: zero covariance between them
  {
    const vector<double> req = { 60.0, 661.0 };
    const vector<double> cov = uncert->efficiencyFracCovariance( req );
    BOOST_CHECK( close_enough( cov[0], 0.07*0.07 ) );
    BOOST_CHECK_EQUAL( cov[1], 0.0 );
    BOOST_CHECK_EQUAL( cov[2], 0.0 );
    BOOST_CHECK( close_enough( cov[3], 0.03*0.03 ) );
  }

  // An energy outside all bands has zero (band) uncertainty
  {
    const vector<double> req = { 5000.0 };
    const vector<double> cov = uncert->efficiencyFracCovariance( req );
    BOOST_CHECK_EQUAL( cov[0], 0.0 );
  }

  // Invalid bands throw
  {
    auto bad = make_shared<DetectorEfficiencyUncert>();
    vector<EffUncertBand> overlapping;
    overlapping.push_back( EffUncertBand{ 0.0f, 100.0f, 0.05f } );
    overlapping.push_back( EffUncertBand{ 50.0f, 200.0f, 0.05f } );
    BOOST_CHECK_THROW( bad->setBands( overlapping ), std::runtime_error );

    vector<EffUncertBand> inverted;
    inverted.push_back( EffUncertBand{ 100.0f, 50.0f, 0.05f } );
    BOOST_CHECK_THROW( bad->setBands( inverted ), std::runtime_error );

    vector<EffUncertBand> negative;
    negative.push_back( EffUncertBand{ 0.0f, 100.0f, -0.05f } );
    BOOST_CHECK_THROW( bad->setBands( negative ), std::runtime_error );
  }
}//test_band_covariance_semantics


BOOST_AUTO_TEST_CASE( test_uncert_xml_roundtrip )
{
  const vector<float> energies = { 59.5f, 122.0f, 661.7f, 1332.5f };
  const vector<float> uncerts = { 0.09f, 0.06f, 0.045f, 0.055f };

  shared_ptr<DetectorEfficiencyUncert> orig
              = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts, 0.5 );

  vector<EffUncertBand> bands;
  bands.push_back( EffUncertBand{ 50.0f, 122.0f, 0.08f } );
  bands.push_back( EffUncertBand{ 122.0f, 661.0f, 0.05f } );
  // Use a const_cast-free copy with bands set before sharing
  auto orig_mutable = make_shared<DetectorEfficiencyUncert>( *orig );
  orig_mutable->setBands( bands );
  orig_mutable->setCoefficientCovariance( { 1.0E-4f, -2.0E-5f, -2.0E-5f, 4.0E-5f } );

  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *parent = doc.allocate_node( rapidxml::node_element, "Parent" );
  doc.append_node( parent );
  orig_mutable->toXml( parent, &doc );

  const rapidxml::xml_node<char> *uncert_node = parent->first_node( "EfficiencyUncert" );
  BOOST_REQUIRE( uncert_node );

  DetectorEfficiencyUncert decoded;
  BOOST_REQUIRE_NO_THROW( decoded.fromXml( uncert_node ) );

  // Float parsing may differ in the last ULP, so use the tolerant comparison
  BOOST_CHECK_NO_THROW( DetectorEfficiencyUncert::equalEnough( decoded, *orig_mutable ) );
  BOOST_CHECK( decoded.hasBands() );
  BOOST_CHECK( decoded.hasNodeCovariance() );
  BOOST_CHECK_EQUAL( decoded.coefficientCovariance().size(), 4u );
  BOOST_CHECK( close_enough( decoded.correlationLength(), 0.5 ) );
}//test_uncert_xml_roundtrip


BOOST_AUTO_TEST_CASE( test_uncert_url_parts_roundtrip )
{
  const vector<float> energies = { 59.5f, 122.0f, 661.7f, 1332.5f };
  const vector<float> uncerts = { 0.09f, 0.06f, 0.045f, 0.055f };

  shared_ptr<DetectorEfficiencyUncert> tmp
              = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts, 0.5 );
  auto orig = make_shared<DetectorEfficiencyUncert>( *tmp );

  vector<EffUncertBand> bands;
  bands.push_back( EffUncertBand{ 50.0f, 122.0f, 0.08f } );
  orig->setBands( bands );

  map<string,string> parts;
  orig->toUrlParts( parts, "" );

  BOOST_CHECK( parts.count("EFUB") );
  BOOST_CHECK( parts.count("EFUE") );
  BOOST_CHECK( parts.count("EFUC") );
  BOOST_CHECK( parts.count("EFUL") );

  shared_ptr<DetectorEfficiencyUncert> decoded = DetectorEfficiencyUncert::fromUrlParts( parts, "" );
  BOOST_REQUIRE( decoded );

  // URL encoding uses limited significant figures, so compare with tolerance.
  BOOST_REQUIRE_EQUAL( decoded->bands().size(), orig->bands().size() );
  BOOST_CHECK( close_enough( decoded->bands()[0].fractionalUncert, 0.08, 1e-4 ) );

  const vector<float> &orig_cov = orig->covarianceMatrix();
  const vector<float> &dec_cov = decoded->covarianceMatrix();
  BOOST_REQUIRE_EQUAL( orig_cov.size(), dec_cov.size() );

  const size_t n = orig->covarianceEnergies().size();
  for( size_t i = 0; i < n; ++i )
  {
    for( size_t j = 0; j < n; ++j )
    {
      BOOST_CHECK_MESSAGE( close_enough( orig_cov[i*n + j], dec_cov[i*n + j], 1e-3 ),
        "URL round-trip covariance [" + to_string(i) + "][" + to_string(j) + "]" );

      // The reconstruction from the upper triangle must be exactly symmetric
      BOOST_CHECK_EQUAL( dec_cov[i*n + j], dec_cov[j*n + i] );
    }
  }

  // No relevant keys -> nullptr, not a throw
  const map<string,string> empty_parts{ {"NAME", "x"} };
  BOOST_CHECK( !DetectorEfficiencyUncert::fromUrlParts( empty_parts, "" ) );

  // Prefixed keys should be found with the matching prefix only
  map<string,string> prefixed;
  orig->toUrlParts( prefixed, "T" );
  BOOST_CHECK( prefixed.count("TEFUB") && prefixed.count("TEFUE") );
  BOOST_CHECK( !DetectorEfficiencyUncert::fromUrlParts( prefixed, "" ) );
  BOOST_CHECK( !!DetectorEfficiencyUncert::fromUrlParts( prefixed, "T" ) );
}//test_uncert_url_parts_roundtrip


BOOST_AUTO_TEST_CASE( test_efficiency_curve_basics )
{
  // Pairs-based curve evaluation matches DetectorPeakResponse::akimaInterpolate
  vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
  for( size_t i = 0; i < 12; ++i )
  {
    DetectorPeakResponse::EnergyEfficiencyPair p;
    p.energy = 50.0f + 250.0f*i;
    p.efficiency = 0.8f * exp( -0.0008f * p.energy ) + 0.05f;
    pairs.push_back( p );
  }

  DetectorEfficiencyCurve curve;
  BOOST_CHECK( !curve.isValid() );
  BOOST_CHECK_THROW( curve.efficiency( 661.0f ), std::runtime_error );

  curve.setFromPairs( pairs, static_cast<float>(PhysicalUnits::keV) );
  BOOST_CHECK( curve.isValid() );

  for( const float energy : { 75.0f, 661.0f, 1332.0f, 2614.0f } )
  {
    const float expected = DetectorPeakResponse::akimaInterpolate( energy, pairs );
    BOOST_CHECK( close_enough( curve.efficiency(energy), expected, 1e-6 ) );
  }

  // Exp-of-log-power-series curve matches the static evaluator
  DetectorEfficiencyCurve series_curve;
  const vector<float> coefs = { -2.0f, 0.6f, -0.12f };
  series_curve.setFromExpOfLogPowerSeries( coefs, {}, static_cast<float>(PhysicalUnits::keV) );
  BOOST_CHECK( series_curve.isValid() );
  for( const float energy : { 100.0f, 661.0f, 2000.0f } )
  {
    const float expected = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, coefs );
    BOOST_CHECK( close_enough( series_curve.efficiency(energy), expected, 1e-6 ) );
  }

  // Formula-based curve
  DetectorEfficiencyCurve formula_curve;
  formula_curve.setFromFormula( "0.5*exp(-0.001*x)", static_cast<float>(PhysicalUnits::keV) );
  BOOST_CHECK( formula_curve.isValid() );
  BOOST_CHECK( close_enough( formula_curve.efficiency(661.0f), 0.5*exp(-0.661), 1e-5 ) );

  BOOST_CHECK_THROW( formula_curve.setFromFormula( "not a (valid formula",
                            static_cast<float>(PhysicalUnits::keV) ), std::exception );
}//test_efficiency_curve_basics


BOOST_AUTO_TEST_CASE( test_efficiency_curve_xml_roundtrip )
{
  vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
  for( size_t i = 0; i < 8; ++i )
  {
    DetectorPeakResponse::EnergyEfficiencyPair p;
    p.energy = 60.0f + 333.0f*i;
    p.efficiency = 0.7f * exp( -0.0007f * p.energy );
    pairs.push_back( p );
  }

  auto orig = make_shared<DetectorEfficiencyCurve>();
  orig->setFromPairs( pairs, static_cast<float>(PhysicalUnits::keV) );

  const vector<float> uncert_energies = { 60.0f, 700.0f, 2300.0f };
  const vector<float> uncert_vals = { 0.1f, 0.05f, 0.07f };
  orig->setUncertainty( DetectorEfficiencyUncert::fromPointUncerts( uncert_energies, uncert_vals ) );

  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *parent = doc.allocate_node( rapidxml::node_element, "Parent" );
  doc.append_node( parent );
  orig->toXml( parent, &doc, "TotalEfficiency" );

  const rapidxml::xml_node<char> *node = parent->first_node( "TotalEfficiency" );
  BOOST_REQUIRE( node );

  DetectorEfficiencyCurve decoded;
  BOOST_REQUIRE_NO_THROW( decoded.fromXml( node ) );

  BOOST_CHECK_NO_THROW( DetectorEfficiencyCurve::equalEnough( decoded, *orig ) );
  BOOST_REQUIRE( decoded.uncertainty() );
  BOOST_CHECK_NO_THROW( DetectorEfficiencyUncert::equalEnough( *decoded.uncertainty(),
                                                          *orig->uncertainty() ) );

  // Formula round-trip
  auto formula_orig = make_shared<DetectorEfficiencyCurve>();
  formula_orig->setFromFormula( "0.5*exp(-0.001*x)", static_cast<float>(PhysicalUnits::keV) );

  rapidxml::xml_document<char> doc2;
  rapidxml::xml_node<char> *parent2 = doc2.allocate_node( rapidxml::node_element, "Parent" );
  doc2.append_node( parent2 );
  formula_orig->toXml( parent2, &doc2, "TotalEfficiency" );

  DetectorEfficiencyCurve formula_decoded;
  BOOST_REQUIRE_NO_THROW( formula_decoded.fromXml( parent2->first_node("TotalEfficiency") ) );
  BOOST_CHECK( formula_decoded == *formula_orig );
  BOOST_CHECK( close_enough( formula_decoded.efficiency(661.0f),
                             formula_orig->efficiency(661.0f), 1e-6 ) );
}//test_efficiency_curve_xml_roundtrip


BOOST_AUTO_TEST_CASE( test_efficiency_curve_url_roundtrip )
{
  // Exp-of-log series with uncertainty, prefixed as the total efficiency uses
  auto orig = make_shared<DetectorEfficiencyCurve>();
  const vector<float> coefs = { -2.0f, 0.6f, -0.12f };
  orig->setFromExpOfLogPowerSeries( coefs, {}, 1000.0f );  //MeV units

  const vector<float> uncert_energies = { 60.0f, 700.0f, 2300.0f };
  const vector<float> uncert_vals = { 0.1f, 0.05f, 0.07f };
  orig->setUncertainty( DetectorEfficiencyUncert::fromPointUncerts( uncert_energies, uncert_vals ) );

  map<string,string> parts;
  orig->toUrlParts( parts, "T" );

  BOOST_CHECK( parts.count("TEFT") && parts.count("TEFC") && parts.count("TEUNIT") );
  BOOST_CHECK( parts.count("TEFUE") && parts.count("TEFUC") );

  shared_ptr<DetectorEfficiencyCurve> decoded = DetectorEfficiencyCurve::fromUrlParts( parts, "T" );
  BOOST_REQUIRE( decoded );
  BOOST_CHECK( decoded->form() == DetectorPeakResponse::kExpOfLogPowerSeries );
  BOOST_CHECK( close_enough( decoded->energyUnits(), 1000.0, 1e-3 ) );
  BOOST_REQUIRE( decoded->uncertainty() );

  for( const float energy : { 100.0f, 661.0f, 2000.0f } )
    BOOST_CHECK( close_enough( decoded->efficiency(energy), orig->efficiency(energy), 1e-4 ) );

  // No relevant keys -> nullptr
  const map<string,string> unrelated{ {"EFFT", "E"}, {"EFFC", "1*2"} };
  BOOST_CHECK( !DetectorEfficiencyCurve::fromUrlParts( unrelated, "T" ) );
}//test_efficiency_curve_url_roundtrip


BOOST_AUTO_TEST_CASE( test_uncert_hash_and_equality )
{
  const vector<float> energies = { 60.0f, 661.0f, 2614.0f };
  const vector<float> uncerts = { 0.08f, 0.05f, 0.06f };

  shared_ptr<DetectorEfficiencyUncert> a = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts );
  shared_ptr<DetectorEfficiencyUncert> b = DetectorEfficiencyUncert::fromPointUncerts( energies, uncerts );

  BOOST_CHECK( *a == *b );

  size_t seed_a = 0, seed_b = 0;
  a->appendToHash( seed_a );
  b->appendToHash( seed_b );
  BOOST_CHECK_EQUAL( seed_a, seed_b );

  // Different content gives different hash
  auto c = make_shared<DetectorEfficiencyUncert>( *a );
  vector<EffUncertBand> bands;
  bands.push_back( EffUncertBand{ 50.0f, 200.0f, 0.05f } );
  c->setBands( bands );
  BOOST_CHECK( !(*a == *c) );

  size_t seed_c = 0;
  c->appendToHash( seed_c );
  BOOST_CHECK_NE( seed_a, seed_c );
}//test_uncert_hash_and_equality
