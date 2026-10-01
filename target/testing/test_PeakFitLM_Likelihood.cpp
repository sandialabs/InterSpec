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

/** Tests of PeakFitLM's Poisson-likelihood fits (by default of sparse ROIs; of every ROI with
 `PeakFitLMOptions::ForcePoissonLikelihood`), its marginal area uncertainties, and the deviance
 residual in InterSpec/PeakFitLMObjective_imp.hpp.  The likelihood fit is IRLS, or all-in-Ceres
 with `SPARSE_DATA_LIKELIHOOD_USE_CERES` (src/PeakFitLM.cpp); these tests hold for either.
 The statistical comparison over many configurations lives in the
 target/peak_fit_improve_ai/harness/peak_fit_objective_eval tool; these are the correctness checks.
 */

#include "InterSpec_config.h"

#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#define BOOST_TEST_MODULE PeakFitLM_Likelihood_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/PeakFitLMObjective_imp.hpp"

using namespace std;
using namespace boost::unit_test;

namespace
{

void set_data_dir()
{
  static bool s_set = false;
  if( s_set )
    return;
  s_set = true;

  string data_dir;
  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      data_dir = arg.substr( 10 );
  }

  if( data_dir.empty() || !SpecUtils::is_file( SpecUtils::append_path( data_dir, "sandia.decay.xml" ) ) )
  {
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
      {
        data_dir = d;
        break;
      }
    }
  }

  if( !data_dir.empty() )
    InterSpec::setStaticDataDirectory( data_dir );
}//set_data_dir()


/** Knuth's product method; portable (unlike std::poisson_distribution) and fine for the means used here. */
int poisson_deviate( const double mean, mt19937 &rng )
{
  if( !(mean > 0.0) )
    return 0;
  const double limit = std::exp( -mean );
  double product = 1.0;
  int count = 0;
  for( ; count < 100000; ++count )
  {
    product *= (static_cast<double>( rng() ) + 0.5) / 4294967296.0;
    if( product <= limit )
      break;
  }
  return count;
}


double deviance_term( const double n, const double m )
{
  return (n > 0.0) ? 2.0*(m - n + n*std::log(n/m)) : 2.0*m;
}


/** An HPGe-like spectrum around 661.657 keV: one Gaussian on a linear continuum. */
struct SyntheticPeak
{
  double area = 100.0;
  double continuum_per_channel = 1.0;
  double fwhm = 1.9;            // keV
  double fwhm_channels = 6.0;
  double mean = 661.657;

  /** Optional second peak this many FWHM above the first (0 = none), with area `area2`. */
  double second_offset_fwhm = 0.0;
  double area2 = 0.0;

  double sigma() const { return fwhm / 2.35482; }
  double mean2() const { return mean + second_offset_fwhm*fwhm; }
  double channel_width() const { return fwhm / fwhm_channels; }
  double roi_lower() const { return mean - 4.0*fwhm; }
  double roi_upper() const { return ((second_offset_fwhm > 0.0) ? mean2() : mean) + 3.5*fwhm; }

  /** Noise-free expected counts per channel, and the calibration they go with. */
  shared_ptr<SpecUtils::Measurement> expected() const
  {
    const double w = channel_width();
    const double lower = roi_lower() - 3.0*fwhm - 0.37*w;
    const size_t nchannel = static_cast<size_t>( std::ceil( (roi_upper() - roi_lower() + 6.0*fwhm)/w ) );
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( nchannel, { static_cast<float>(lower), static_cast<float>(w) }, {} );

    PeakDef peak( mean, sigma(), area );
    const shared_ptr<const vector<float>> energies = cal->channel_energies();
    vector<double> counts( nchannel, 0.0 );
    peak.gauss_integral( energies->data(), counts.data(), nchannel );
    if( second_offset_fwhm > 0.0 )
      PeakDef( mean2(), sigma(), area2 ).gauss_integral( energies->data(), counts.data(), nchannel );
    for( size_t i = 0; i < nchannel; ++i )
    {
      const double x = 0.5*((*energies)[i] + (*energies)[i+1]) - mean;
      counts[i] += continuum_per_channel * (1.0 - 0.3*x/(roi_upper() - roi_lower()));
    }

    auto meas = make_shared<SpecUtils::Measurement>();
    meas->set_gamma_counts( make_shared<const vector<float>>( counts.begin(), counts.end() ), 1.0f, 1.0f );
    meas->set_energy_calibration( cal );
    return meas;
  }//expected()

  /** A Poisson replica of `expected()`. */
  shared_ptr<const SpecUtils::Measurement> sample( const shared_ptr<const SpecUtils::Measurement> &expected,
                                                   mt19937 &rng ) const
  {
    vector<float> counts = *expected->gamma_counts();
    for( float &c : counts )
      c = static_cast<float>( poisson_deviate( c, rng ) );
    auto meas = make_shared<SpecUtils::Measurement>( *expected );
    meas->set_gamma_counts( make_shared<const vector<float>>( std::move(counts) ), 1.0f, 1.0f );
    return meas;
  }

  /** Starting peak for a fit: the truth, nudged by `mean_offset` sigma and `fwhm_factor`. */
  shared_ptr<const PeakDef> start_peak( const double mean_offset, const double fwhm_factor ) const
  {
    auto cont = make_shared<PeakContinuum>();
    cont->setType( PeakContinuum::OffsetType::Linear );
    cont->setRange( roi_lower(), roi_upper() );
    auto p = make_shared<PeakDef>( mean + mean_offset*sigma(), fwhm_factor*sigma(), area );
    p->setContinuum( cont );
    return p;
  }

  /** Both starting peaks of a doublet, sharing one continuum. */
  vector<shared_ptr<const PeakDef>> start_doublet( const double mean_offset, const double fwhm_factor ) const
  {
    auto first = make_shared<PeakDef>( *start_peak( mean_offset, fwhm_factor ) );
    auto second = make_shared<PeakDef>( mean2() - mean_offset*sigma(), fwhm_factor*sigma(), area2 );
    second->setContinuum( first->getContinuum() );
    return { first, second };
  }
};//struct SyntheticPeak


/** The Poisson deviance of a ROI's fit (peaks sharing one continuum), over the ROI's channels. */
double roi_deviance( const vector<shared_ptr<const PeakDef>> &peaks,
                     const shared_ptr<const SpecUtils::Measurement> &data )
{
  const shared_ptr<const PeakContinuum> cont = peaks[0]->continuum();
  const size_t ch0 = data->find_gamma_channel( static_cast<float>( cont->lowerEnergy() ) );
  const size_t ch1 = data->find_gamma_channel( static_cast<float>( cont->upperEnergy() ) );
  const size_t nchan = ch1 - ch0 + 1;
  const float * const x = data->channel_energies()->data() + ch0;

  vector<const PeakDef *> roi_peaks;
  for( const shared_ptr<const PeakDef> &p : peaks )
    roi_peaks.push_back( p.get() );

  vector<double> model( nchan, 0.0 );
  cont->offset_integral( x, model.data(), nchan, data, roi_peaks.data(), roi_peaks.size() );
  for( const shared_ptr<const PeakDef> &p : peaks )
    p->gauss_integral( x, model.data(), nchan );

  double deviance = 0.0;
  for( size_t i = 0; i < nchan; ++i )
  {
    const double n = std::max( 0.0, static_cast<double>( (*data->gamma_counts())[ch0 + i] ) );
    deviance += deviance_term( n, std::max( model[i], 1.0E-12 ) );
  }
  return deviance;
}//roi_deviance(...)


/** For a single-peak fit on a polynomial continuum: per parameter (area, mean, sigma, then each
 continuum coefficient), how far the Poisson deviance's minimum along that parameter lies from the
 fit value, in units of that parameter's own (conditional) uncertainty - `|g|/sqrt(2c)`, with g and c
 the deviance's first and second derivatives, by central differences.  Zero at the Poisson MLE.
 */
vector<double> distance_from_deviance_minimum( const shared_ptr<const PeakDef> &fit_peak,
                                               const shared_ptr<const SpecUtils::Measurement> &data )
{
  const size_t ncont = fit_peak->continuum()->parameters().size();
  vector<double> answer;
  for( size_t par = 0; par < 3 + ncont; ++par )
  {
    const auto shifted = [&]( const double frac ) -> double {
      auto p = make_shared<PeakDef>( *fit_peak );
      p->makeUniqueNewContinuum();
      double value = 0.0, uncert = 0.0;
      switch( par )
      {
        case 0: value = p->amplitude(); uncert = p->amplitudeUncert(); break;
        case 1: value = p->mean();      uncert = p->meanUncert();      break;
        case 2: value = p->sigma();     uncert = p->sigmaUncert();     break;
        default:
          value = p->continuum()->parameters()[par - 3];
          uncert = p->continuum()->uncertainties()[par - 3];
          break;
      }
      const double h = frac * 0.01 * uncert;
      switch( par )
      {
        case 0: p->setAmplitude( value + h ); break;
        case 1: p->setMean( value + h );      break;
        case 2: p->setSigma( value + h );     break;
        default: p->getContinuum()->setPolynomialCoef( par - 3, value + h ); break;
      }
      return roi_deviance( { p }, data );
    };

    const double plus = shifted( 1.0 ), minus = shifted( -1.0 ), center = shifted( 0.0 );
    const double g = (plus - minus) / 2.0;              // per step
    const double c = plus - 2.0*center + minus;         // per step^2
    answer.push_back( (c > 0.0) ? (fabs(g) / std::sqrt( 2.0*c )) : std::numeric_limits<double>::infinity() );
  }//for( par )

  return answer;
}//distance_from_deviance_minimum(...)

}//namespace


BOOST_AUTO_TEST_CASE( deviance_residual_values_and_derivatives )
{
  using Jet = ceres::Jet<double,1>;
  const double floor = 1.0E-8;

  for( const double n : { 0.0, 1.0, 3.0, 25.0, 1000.0 } )
  {
    vector<double> models{ 0.01, 0.5, 7.0 };
    if( n > 0.0 )
    {
      for( const double f : { 1.0, 1.0 + 1.0E-7, 1.0 - 1.0E-7, 1.0 + 5.0E-4, 1.0 - 5.0E-4, 1.0 + 2.0E-3, 0.5, 2.0 } )
        models.push_back( n*f );
    }

    for( const double m : models )
    {
      const double r = PeakFitLMObjective::poisson_deviance_residual( n, m, floor );
      const double dev = deviance_term( n, m );
      BOOST_CHECK_SMALL( r*r - dev, 1.0E-9*std::max( 1.0, dev ) + 1.0E-12 );
      if( fabs(n - m) > 1.0E-9*std::max( 1.0, n ) )
        BOOST_CHECK_MESSAGE( (r > 0.0) == (n > m), "sign of residual wrong for n=" << n << ", m=" << m );

      // Jet derivative vs a central difference (both sides on the same branch of the series switch).
      const Jet rj = PeakFitLMObjective::poisson_deviance_residual( n, Jet( m, 0 ), floor );
      BOOST_REQUIRE( std::isfinite( rj.v[0] ) );
      const double h = 1.0E-7 * std::max( m, 1.0E-3 );
      const double fd = (PeakFitLMObjective::poisson_deviance_residual( n, m + h, floor )
                         - PeakFitLMObjective::poisson_deviance_residual( n, m - h, floor )) / (2.0*h);
      BOOST_CHECK_MESSAGE( fabs( rj.v[0] - fd ) <= 1.0E-4*std::max( 1.0, fabs(fd) ),
                           "n=" << n << ", m=" << m << ": Jet derivative " << rj.v[0] << " vs FD " << fd );
    }//for( models )
  }//for( n )

  // A model at, or below, zero stays finite, with a finite derivative.
  for( const double m : { 0.0, -1.0, 1.0E-12 } )
  {
    const Jet rj = PeakFitLMObjective::poisson_deviance_residual( 5.0, Jet( m, 0 ), floor );
    BOOST_CHECK( std::isfinite( rj.a ) && std::isfinite( rj.v[0] ) );
  }
}//deviance_residual_values_and_derivatives


BOOST_AUTO_TEST_CASE( likelihood_fit_is_the_poisson_mle )
{
  // The likelihood fit's stationarity conditions are the Poisson likelihood equations, for every
  //  parameter: along each one, the deviance's minimum is (to the convergence tolerance) at the fit.
  //  The chi2 fit of sparse data is not, which shows the check can tell the difference.
  set_data_dir();
  const PeakFitUtils::CoarseResolutionType det = PeakFitUtils::CoarseResolutionType::High;
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> forced( PeakFitLM::ForcePoissonLikelihood );
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> chi2_only( PeakFitLM::NoSparseDataLikelihood );
  const char * const par_names[] = { "area", "mean", "sigma", "cont0", "cont1" };

  for( const pair<double,double> &area_cont : { pair<double,double>( 2000.0, 20.0 ), pair<double,double>( 60.0, 0.5 ) } )
  {
    SyntheticPeak spec;
    spec.area = area_cont.first;
    spec.continuum_per_channel = area_cont.second;
    mt19937 rng( 7 );
    const shared_ptr<const SpecUtils::Measurement> data = spec.sample( spec.expected(), rng );

    PeakFitLM::take_fit_objective_diagnostics();
    const vector<shared_ptr<const PeakDef>> mle = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.2, 1.05 ) }, data, det, forced );
    const PeakFitLM::FitObjectiveDiagnostics diag = PeakFitLM::take_fit_objective_diagnostics();
    const vector<shared_ptr<const PeakDef>> chi2 = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.2, 1.05 ) }, data, det, chi2_only );
    BOOST_REQUIRE( (mle.size() == 1) && (chi2.size() == 1) );
    BOOST_CHECK( !diag.fell_back_to_chi2 );
    BOOST_CHECK( (diag.irls_passes == 0) || diag.irls_converged );

    const vector<double> mle_dist = distance_from_deviance_minimum( mle[0], data );
    const vector<double> chi2_dist = distance_from_deviance_minimum( chi2[0], data );
    double chi2_max = 0.0;
    for( size_t i = 0; i < mle_dist.size(); ++i )
    {
      cout << "area " << spec.area << ", cont " << spec.continuum_per_channel << ": " << par_names[i]
           << " deviance minimum is " << mle_dist[i] << " sd from the likelihood fit (chi2 fit: "
           << chi2_dist[i] << ")" << endl;
      BOOST_CHECK_MESSAGE( mle_dist[i] < 0.05, "area " << spec.area << ": " << par_names[i]
                           << " is " << mle_dist[i] << " sd from the deviance minimum" );
      chi2_max = std::max( chi2_max, chi2_dist[i] );
    }
    if( spec.continuum_per_channel < 1.0 )
      BOOST_CHECK_MESSAGE( chi2_max > 0.2, "chi2 fit of sparse data unexpectedly at the deviance minimum" );
  }//for( statistics )
}//likelihood_fit_is_the_poisson_mle


BOOST_AUTO_TEST_CASE( likelihood_objectives_remove_low_count_bias )
{
  // 50 counts on half a count per channel: the modified-Neyman chi2 is well known to under-estimate
  //  the area here; the Poisson maximum-likelihood fits should not.
  set_data_dir();

  SyntheticPeak spec;
  spec.area = 50.0;
  spec.continuum_per_channel = 0.5;
  const shared_ptr<const SpecUtils::Measurement> expected = spec.expected();
  const PeakFitUtils::CoarseResolutionType det = PeakFitUtils::CoarseResolutionType::High;

  const vector<pair<string,Wt::WFlags<PeakFitLM::PeakFitLMOptions>>> objectives{
    { "chi2", PeakFitLM::NoSparseDataLikelihood },
    { "likelihood", PeakFitLM::ForcePoissonLikelihood },
    { "default", {} }   // sparse ROI: refit by likelihood
  };

  const size_t nreplica = 200;
  vector<vector<double>> areas( objectives.size() );
  mt19937 rng( 2024 );
  for( size_t r = 0; r < nreplica; ++r )
  {
    const shared_ptr<const SpecUtils::Measurement> data = spec.sample( expected, rng );
    const double mean_offset = ((r % 5) - 2.0)*0.1, fwhm_factor = 0.95 + 0.025*(r % 5);
    for( size_t i = 0; i < objectives.size(); ++i )
    {
      try
      {
        const vector<shared_ptr<const PeakDef>> fit
                     = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( mean_offset, fwhm_factor ) },
                                                       data, det, objectives[i].second );
        if( fit.size() == 1 )
          areas[i].push_back( fit[0]->amplitude() );
      }catch( std::exception & )
      {
      }
    }
  }//for( replicas )

  vector<double> z( objectives.size() );
  for( size_t i = 0; i < objectives.size(); ++i )
  {
    BOOST_REQUIRE( areas[i].size() > 0.95*nreplica );
    double sum = 0.0, sum2 = 0.0;
    for( const double a : areas[i] )
    {
      sum += a;
      sum2 += a*a;
    }
    const double n = static_cast<double>( areas[i].size() );
    const double mean = sum/n, sd = std::sqrt( std::max( sum2/n - mean*mean, 0.0 ) );
    z[i] = (mean - spec.area) / (sd/std::sqrt(n));
    cout << objectives[i].first << ": mean area " << mean << " (truth " << spec.area << "), bias "
         << z[i] << " standard errors" << endl;
  }

  BOOST_CHECK_MESSAGE( z[0] < -3.0, "chi2 should be biased low here; z=" << z[0] );
  BOOST_CHECK_MESSAGE( fabs( z[1] ) < 3.0, "likelihood fit biased; z=" << z[1] );
  BOOST_CHECK_MESSAGE( fabs( z[2] ) < 3.0, "default fit biased; z=" << z[2] );
}//likelihood_objectives_remove_low_count_bias


BOOST_AUTO_TEST_CASE( objective_option_handling )
{
  set_data_dir();

  SyntheticPeak spec;
  mt19937 rng( 5 );
  const shared_ptr<const SpecUtils::Measurement> data = spec.sample( spec.expected(), rng );
  const PeakFitUtils::CoarseResolutionType det = PeakFitUtils::CoarseResolutionType::High;

  // Forcing the likelihood and turning it off at once is a programming error.
  BOOST_CHECK_THROW( PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.0, 1.0 ) }, data, det,
                       Wt::WFlags<PeakFitLM::PeakFitLMOptions>( PeakFitLM::ForcePoissonLikelihood )
                         | PeakFitLM::NoSparseDataLikelihood ),
                     std::exception );

  // A negative channel (e.g., a background-subtracted spectrum) is not Poisson data: the fit is the
  //  chi2 one.
  vector<float> counts = *data->gamma_counts();
  counts[2] = -3.0f;
  auto subtracted = make_shared<SpecUtils::Measurement>( *data );
  subtracted->set_gamma_counts( make_shared<const vector<float>>( counts ), 1.0f, 1.0f );

  PeakFitLM::take_fit_objective_diagnostics();
  const vector<shared_ptr<const PeakDef>> chi2 = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.0, 1.0 ) }, subtracted, det, {} );
  const vector<shared_ptr<const PeakDef>> forced = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.0, 1.0 ) }, subtracted, det,
                                                                                  PeakFitLM::ForcePoissonLikelihood );
  BOOST_CHECK( PeakFitLM::take_fit_objective_diagnostics().fell_back_to_chi2 );
  BOOST_REQUIRE( (chi2.size() == 1) && (forced.size() == 1) );
  BOOST_CHECK_EQUAL( chi2[0]->amplitude(), forced[0]->amplitude() );
  BOOST_CHECK_EQUAL( chi2[0]->mean(), forced[0]->mean() );
}//objective_option_handling


BOOST_AUTO_TEST_CASE( marginal_area_uncertainty )
{
  // With the mean and FWHM fit, the area uncertainty must include their correlation with the area
  //  (`PeakFitLMOptions::ConditionalAreaUncertainties` gives the old, conditional, value): the
  //  reported uncertainty should match the actual spread of the fitted areas over Poisson replicas.
  set_data_dir();

  // A doublet 1.3 FWHM apart: the areas are strongly correlated with the means and width.
  SyntheticPeak spec;
  spec.area = 400.0;
  spec.area2 = 400.0;
  spec.second_offset_fwhm = 1.3;
  spec.continuum_per_channel = 2.0;
  const shared_ptr<const SpecUtils::Measurement> expected = spec.expected();
  const PeakFitUtils::CoarseResolutionType det = PeakFitUtils::CoarseResolutionType::High;
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> chi2( PeakFitLM::NoSparseDataLikelihood );
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> chi2_cond = chi2 | PeakFitLM::ConditionalAreaUncertainties;

  vector<double> areas, marg_uncerts, cond_uncerts;
  mt19937 rng( 31 );
  for( size_t r = 0; r < 300; ++r )
  {
    const shared_ptr<const SpecUtils::Measurement> data = spec.sample( expected, rng );
    const vector<shared_ptr<const PeakDef>> marg = PeakFitLM::fit_peaks_in_roi_LM( spec.start_doublet( 0.1, 1.05 ), data, det, chi2 );
    const vector<shared_ptr<const PeakDef>> cond = PeakFitLM::fit_peaks_in_roi_LM( spec.start_doublet( 0.1, 1.05 ), data, det, chi2_cond );
    BOOST_REQUIRE( (marg.size() == 2) && (cond.size() == 2) );
    BOOST_REQUIRE_EQUAL( marg[0]->amplitude(), cond[0]->amplitude() );   // only the uncertainty differs
    BOOST_CHECK( marg[0]->amplitudeUncert() >= cond[0]->amplitudeUncert() );
    areas.push_back( marg[0]->amplitude() );
    marg_uncerts.push_back( marg[0]->amplitudeUncert() );
    cond_uncerts.push_back( cond[0]->amplitudeUncert() );
  }

  double sum = 0.0, sum2 = 0.0;
  for( const double a : areas )
  {
    sum += a;
    sum2 += a*a;
  }
  const double n = static_cast<double>( areas.size() );
  const double sd = std::sqrt( (sum2 - sum*sum/n) / (n - 1.0) );
  std::sort( begin(marg_uncerts), end(marg_uncerts) );
  std::sort( begin(cond_uncerts), end(cond_uncerts) );
  const double marg_ratio = marg_uncerts[marg_uncerts.size()/2] / sd;
  const double cond_ratio = cond_uncerts[cond_uncerts.size()/2] / sd;
  cout << "reported/actual area spread: marginal " << marg_ratio << ", conditional " << cond_ratio << endl;

  // ~300 replicas pin the actual spread to about +-4%; the conditional value is well short of it.
  BOOST_CHECK_MESSAGE( (marg_ratio > 0.9) && (marg_ratio < 1.12), "marginal uncertainty / spread = " << marg_ratio );
  BOOST_CHECK_MESSAGE( cond_ratio < 0.85, "conditional uncertainty / spread = " << cond_ratio );

  // With the mean and FWHM held fixed there is nothing non-linear to propagate: identical.
  mt19937 rng2( 32 );
  const shared_ptr<const SpecUtils::Measurement> data = spec.sample( expected, rng2 );
  auto fixed_peak = make_shared<PeakDef>( *spec.start_peak( 0.0, 1.0 ) );
  fixed_peak->setFitFor( PeakDef::Mean, false );
  fixed_peak->setFitFor( PeakDef::Sigma, false );
  const vector<shared_ptr<const PeakDef>> a = PeakFitLM::fit_peaks_in_roi_LM( { fixed_peak }, data, det, {} );
  const vector<shared_ptr<const PeakDef>> b = PeakFitLM::fit_peaks_in_roi_LM( { fixed_peak }, data, det,
                                                                             PeakFitLM::ConditionalAreaUncertainties );
  BOOST_REQUIRE( (a.size() == 1) && (b.size() == 1) );
  BOOST_CHECK_EQUAL( a[0]->amplitudeUncert(), b[0]->amplitudeUncert() );
}//marginal_area_uncertainty


BOOST_AUTO_TEST_CASE( sparse_data_default )
{
  // By default a sparse ROI is refit by Poisson maximum likelihood and a well-populated one is left
  //  exactly as the chi2 fit; `NoSparseDataLikelihood` gives the plain chi2 fit.
  set_data_dir();
  const PeakFitUtils::CoarseResolutionType det = PeakFitUtils::CoarseResolutionType::High;
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> chi2_only( PeakFitLM::NoSparseDataLikelihood );
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> likelihood( PeakFitLM::ForcePoissonLikelihood );

  // Sparse: 60 counts on half a count per channel.
  {
    SyntheticPeak spec;
    spec.area = 60.0;
    spec.continuum_per_channel = 0.5;
    mt19937 rng( 11 );
    const shared_ptr<const SpecUtils::Measurement> data = spec.sample( spec.expected(), rng );

    PeakFitLM::take_fit_objective_diagnostics();
    const vector<shared_ptr<const PeakDef>> def = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.1, 1.05 ) }, data, det, {} );
    const PeakFitLM::FitObjectiveDiagnostics diag = PeakFitLM::take_fit_objective_diagnostics();
    const vector<shared_ptr<const PeakDef>> forced = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.1, 1.05 ) }, data, det, likelihood );
    const vector<shared_ptr<const PeakDef>> chi2 = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.1, 1.05 ) }, data, det, chi2_only );
    BOOST_REQUIRE( (def.size() == 1) && (forced.size() == 1) && (chi2.size() == 1) );
    BOOST_CHECK_EQUAL( diag.sparse_rois, size_t(1) );
    BOOST_CHECK_EQUAL( def[0]->amplitude(), forced[0]->amplitude() );
    BOOST_CHECK( def[0]->amplitude() != chi2[0]->amplitude() );
  }

  // Well populated: 20000 counts on 200 counts per channel - the chi2 fit, bit for bit.
  {
    SyntheticPeak spec;
    spec.area = 20000.0;
    spec.continuum_per_channel = 200.0;
    mt19937 rng( 12 );
    const shared_ptr<const SpecUtils::Measurement> data = spec.sample( spec.expected(), rng );

    PeakFitLM::take_fit_objective_diagnostics();
    const vector<shared_ptr<const PeakDef>> def = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.1, 1.05 ) }, data, det, {} );
    const PeakFitLM::FitObjectiveDiagnostics diag = PeakFitLM::take_fit_objective_diagnostics();
    const vector<shared_ptr<const PeakDef>> chi2 = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.1, 1.05 ) }, data, det, chi2_only );
    BOOST_REQUIRE( (def.size() == 1) && (chi2.size() == 1) );
    BOOST_CHECK_EQUAL( diag.sparse_rois, size_t(0) );
    BOOST_CHECK_EQUAL( def[0]->amplitude(), chi2[0]->amplitude() );
    BOOST_CHECK_EQUAL( def[0]->amplitudeUncert(), chi2[0]->amplitudeUncert() );
    BOOST_CHECK_EQUAL( def[0]->mean(), chi2[0]->mean() );
  }
}//sparse_data_default


BOOST_AUTO_TEST_CASE( refit_of_sparse_roi )
{
  // Refitting a sparse ROI that was fit by chi2 (as the peak search does) gives the likelihood fit:
  //  the refit is judged by the deviance it minimized, not by chi2 - which the chi2 fit it started
  //  from wins by construction, so the refit would be rejected.
  set_data_dir();
  const PeakFitUtils::CoarseResolutionType det = PeakFitUtils::CoarseResolutionType::High;

  SyntheticPeak spec;
  spec.area = 300.0;
  spec.continuum_per_channel = 0.3;
  mt19937 rng( 13 );
  const shared_ptr<const SpecUtils::Measurement> data = spec.sample( spec.expected(), rng );

  const vector<shared_ptr<const PeakDef>> chi2 = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.1, 1.05 ) }, data, det,
                                                                                PeakFitLM::NoSparseDataLikelihood );
  BOOST_REQUIRE( chi2.size() == 1 );

  PeakFitLM::take_fit_objective_diagnostics();
  const vector<shared_ptr<const PeakDef>> refit = PeakFitLM::refitPeaksThatShareROI_LM( data, nullptr, chi2, det );
  BOOST_CHECK_EQUAL( PeakFitLM::take_fit_objective_diagnostics().sparse_rois, size_t(1) );
  const vector<shared_ptr<const PeakDef>> likelihood = PeakFitLM::fit_peaks_in_roi_LM( chi2, data, det,
                                                                                      PeakFitLM::ForcePoissonLikelihood );
  BOOST_REQUIRE( refit.size() == 1 );
  BOOST_REQUIRE( likelihood.size() == 1 );
  BOOST_CHECK_EQUAL( refit[0]->amplitude(), likelihood[0]->amplitude() );
  BOOST_CHECK( refit[0]->amplitude() != chi2[0]->amplitude() );

  // Refitting the default (likelihood) fit lands back on the same minimum, by way of the chi2
  //  solution - within the convergence tolerance, which must not count as a worse fit.
  spec.area = 1000.0;
  spec.continuum_per_channel = 0.1;
  spec.fwhm_channels = 15.0;
  const shared_ptr<const SpecUtils::Measurement> expected = spec.expected();
  size_t nrefit = 0;
  for( size_t r = 0; r < 10; ++r )
  {
    const shared_ptr<const SpecUtils::Measurement> replica = spec.sample( expected, rng );
    PeakFitLM::take_fit_objective_diagnostics();
    const vector<shared_ptr<const PeakDef>> fit = PeakFitLM::fit_peaks_in_roi_LM( { spec.start_peak( 0.1, 1.05 ) }, replica, det );
    BOOST_REQUIRE( (fit.size() == 1) && (PeakFitLM::take_fit_objective_diagnostics().sparse_rois == 1) );
    nrefit += PeakFitLM::refitPeaksThatShareROI_LM( replica, nullptr, fit, det ).size();
  }
  BOOST_CHECK_EQUAL( nrefit, size_t(10) );
}//refit_of_sparse_roi


BOOST_AUTO_TEST_CASE( skew_of_fixed_amplitude_and_non_lls_peaks )
{
  // A skewed peak is in the model with its skew both when its amplitude is held fixed (it is then a
  //  fixed part of the linear solve) and when a pinned continuum coefficient takes the ROI off the
  //  linear solve - both used to be evaluated before the skew was set.  On a noise-free spectrum
  //  the fit then recovers the truth.
  set_data_dir();
  const PeakFitUtils::CoarseResolutionType det = PeakFitUtils::CoarseResolutionType::High;

  const double fwhm = 1.9, sigma = fwhm/2.35482, width = fwhm/6.0, cont_per_channel = 50.0;
  const double mean1 = 661.657, mean2 = mean1 + 1.5*fwhm, area1 = 5000.0, area2 = 3000.0;
  const double roi_lower = mean1 - 4.0*fwhm, roi_upper = mean2 + 3.5*fwhm;

  const auto skewed_peak = [&]( const double mean, const double area ) -> shared_ptr<PeakDef> {
    auto p = make_shared<PeakDef>( mean, sigma, area );
    p->setSkewType( PeakDef::SkewType::Bortel );
    p->set_coefficient( 1.0, PeakDef::SkewPar0 );
    p->setFitFor( PeakDef::SkewPar0, false );
    p->setFitFor( PeakDef::Sigma, false );
    return p;
  };

  // Noise-free spectrum: the two peaks on a flat continuum.
  const size_t nchannel = static_cast<size_t>( std::ceil( (roi_upper - roi_lower + 8.0*fwhm)/width ) );
  const float lower_edge = static_cast<float>( roi_lower - 4.0*fwhm - 0.37*width );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( nchannel, { lower_edge, static_cast<float>(width) }, {} );
  vector<double> expected( nchannel, cont_per_channel );
  skewed_peak( mean1, area1 )->gauss_integral( cal->channel_energies()->data(), expected.data(), nchannel );
  skewed_peak( mean2, area2 )->gauss_integral( cal->channel_energies()->data(), expected.data(), nchannel );
  auto data = make_shared<SpecUtils::Measurement>();
  data->set_gamma_counts( make_shared<const vector<float>>( expected.begin(), expected.end() ), 1.0f, 1.0f );
  data->set_energy_calibration( cal );

  const auto start_peaks = [&]( const bool fix_area2, const bool pin_continuum ) -> vector<shared_ptr<const PeakDef>> {
    auto cont = make_shared<PeakContinuum>();
    cont->setType( PeakContinuum::OffsetType::Linear );
    cont->setRange( roi_lower, roi_upper );
    if( pin_continuum )
    {
      cont->setParameters( roi_lower, { cont_per_channel/width, 0.0 }, {} );
      cont->setPolynomialCoefFitFor( 0, false );
      cont->setPolynomialCoefFitFor( 1, false );
    }
    shared_ptr<PeakDef> p1 = skewed_peak( mean1 + 0.1*sigma, 0.9*area1 );
    shared_ptr<PeakDef> p2 = skewed_peak( mean2, fix_area2 ? area2 : 1.1*area2 );
    p2->setFitFor( PeakDef::Mean, false );
    p2->setFitFor( PeakDef::GaussAmplitude, !fix_area2 );
    p1->setContinuum( cont );
    p2->setContinuum( cont );
    return { p1, p2 };
  };

  for( const bool pin_continuum : { false, true } )
  {
    const bool fix_area2 = !pin_continuum;
    const vector<shared_ptr<const PeakDef>> fit = PeakFitLM::fit_peaks_in_roi_LM( start_peaks( fix_area2, pin_continuum ),
                                                        data, det, PeakFitLM::NoSparseDataLikelihood );
    BOOST_REQUIRE( fit.size() == 2 );
    BOOST_CHECK_MESSAGE( fabs( fit[0]->amplitude()/area1 - 1.0 ) < 1.0E-3,
                         (pin_continuum ? "pinned continuum" : "fixed second area") << ": area "
                         << fit[0]->amplitude() << " vs truth " << area1 );
    BOOST_CHECK_MESSAGE( fabs( fit[0]->mean() - mean1 ) < 1.0E-3*sigma,
                         (pin_continuum ? "pinned continuum" : "fixed second area") << ": mean "
                         << fit[0]->mean() << " vs truth " << mean1 );
    if( pin_continuum )
      BOOST_CHECK_MESSAGE( fabs( fit[1]->amplitude()/area2 - 1.0 ) < 1.0E-3,
                           "pinned continuum: second area " << fit[1]->amplitude() << " vs truth " << area2 );
  }//for( pin_continuum )
}//skew_of_fixed_amplitude_and_non_lls_peaks
