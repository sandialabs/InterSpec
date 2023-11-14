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

#include <string>
#include <iostream>


//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PeakDists_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/DateTime.h"

#include "InterSpec/PeakDists.h"



using namespace std;
using namespace PeakDists;
using namespace boost::unit_test;


BOOST_AUTO_TEST_CASE( GaussianDist )
{
  // Check Gaussian has unit area
  {
    double mean = 100;
    double sigma = 5;
    double x0 = mean - 10*sigma;
    double x1 = mean + 10*sigma;
    
    double answer = gaussian_integral( mean, sigma, x0, x1 );
    // Note that `BOOST_CHECK_CLOSE(...)` checks for percentages, so 1E-10% is really fractionally 1E-12
    BOOST_CHECK_CLOSE( answer, 1.0, 1.0E-10 );
    
    
    x0 = mean - sigma;
    x1 = mean + sigma;
    answer = gaussian_integral( mean, sigma, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 0.682689492137, 1.0E-10 );
    
    x0 = mean - 2.0*sigma;
    x1 = mean + 2.0*sigma;
    answer = gaussian_integral( mean, sigma, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 0.954499736104, 1.0E-10 );
    
    x0 = mean - 3.0*sigma;
    x1 = mean + 3.0*sigma;
    answer = gaussian_integral( mean, sigma, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 0.997300203937, 1.0E-10 );
  }
  
  // Check fast Gaussian has unit area
  {
    double mean = 100;
    double sigma = 5;
    const double amplitude = 1.1231;
    double x0 = mean - 10*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 1024;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    gaussian_integral( mean, sigma, amplitude, &(energies[0]), &(counts[0]), num_channels );
    
    double answer_sum = std::accumulate( begin(counts), end(counts), 0.0 );
    BOOST_CHECK_CLOSE( answer_sum, amplitude, amplitude*1.0E-10 );
  }
  
}//BOOST_AUTO_TEST_CASE( GaussianDist )



BOOST_AUTO_TEST_CASE( BortelDist )
{
  // Check Bortel has unit area
  {
    double mean = 100;
    double sigma = 5;
    double skew = 0.5;
    double x0 = mean - 10*sigma;
    double x1 = mean + 10*sigma;
    
    // `bortel_integral(...)` calls `bortel_indefinite_integral(...)`
    double answer = bortel_integral( mean, sigma, skew, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 1.0, 1.0E-10 );
  }
  
  // Check fast Bortel has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew = 0.3187;
    const double amplitude = 1.2521;
    double x0 = mean - 30*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 1024;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    bortel_integral( mean, sigma, amplitude, skew, &(energies[0]), &(counts[0]), num_channels );
    
    double answer_sum = std::accumulate( begin(counts), end(counts), 0.0 );
    BOOST_CHECK_CLOSE( answer_sum, amplitude, amplitude*1.0E-10 );
  }
  
  
  // Check Bortel PDF has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew = 0.3187;
    const double amplitude = 1.1235;
    double x0 = mean - 30*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 2048;
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = bortel_pdf( mean, sigma, skew, x );
      
      answer_sum += dx * pdf_val;
    }
    
    BOOST_CHECK_CLOSE( answer_sum, 1.0, 1.0E-7 );
  }
  
  { //test bortel_coverage_limits
    double mean = 100;
    double sigma = 1;
    double skew = 0.5;
    double prob = 0.000000573303;
    
    pair<double,double> limits = bortel_coverage_limits( mean, sigma, skew, prob );
    double fraction = 1.0 - bortel_integral( mean, sigma, skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.317310507863;
    limits = bortel_coverage_limits( mean, sigma, skew, prob );
    fraction = 1.0 - bortel_integral( mean, sigma, skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.1 );
  }
  
}//BOOST_AUTO_TEST_CASE( BortelDist )


BOOST_AUTO_TEST_CASE( GaussExp )
{
  // Check GaussExp has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew = 1.9055;
    double x0 = mean - 20*sigma;
    double x1 = mean + 10*sigma;
    
    double answer = gauss_exp_integral( mean, sigma, skew, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 1.0, 1.0E-10 );
  }
  
  // Check fast GaussExp has unit area
  {
    
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew = 1.9055;
    const double amplitude = 1.2345;
    double x0 = mean - 30*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 2*2*1024;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    //const double start_cpu = SpecUtils::get_cpu_time();
    
    gauss_exp_integral( mean, sigma, amplitude, skew, &(energies[0]), &(counts[0]), num_channels );
    
    //const double end_cpu = SpecUtils::get_cpu_time();
    //cout << "GaussExp Took " << (end_cpu - start_cpu) << " seconds" << endl;
    
    double answer_sum = std::accumulate( begin(counts), end(counts), 0.0 );
    BOOST_CHECK_CLOSE( answer_sum, amplitude, amplitude*1.0E-10 );
  }
  
  
  // Check GaussExp PDF has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew = 1.9055;
    const double amplitude = 1.2345;
    double x0 = mean - 40*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 2048;
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = gauss_exp_pdf( mean, sigma, skew, x );
      
      answer_sum += dx * pdf_val;
    }
    
    BOOST_CHECK_CLOSE( answer_sum, 1.0, 1.0E-6 );
  }
  
  
  {
    double mean = 100;
    double sigma = 1;
    double skew = 0.5; //0.15, 3.25;
    double prob = 0.000000573303;
    
    pair<double,double> limits = gauss_exp_coverage_limits( mean, sigma, skew, prob );
    double fraction = 1.0 - gauss_exp_integral( mean, sigma, skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.317310507863;
    limits = gauss_exp_coverage_limits( mean, sigma, skew, prob );
    fraction = 1.0 - gauss_exp_integral( mean, sigma, skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
  }
  
  /*
  double gauss_exp_norm( const double sigma, const double skew );
  double gauss_exp_pdf(const double mean,
                       const double sigma,
                       const double skew,
                       const double x );
  
  double gauss_exp_tail_indefinite(const double mean,
                                  const double sigma,
                                  const double skew,
                                   const double x );
  
  double gauss_exp_indefinite(const double mean,
                             const double sigma,
                             const double skew,
                              const double x );
   */
}//BOOST_AUTO_TEST_CASE( GaussExp )


BOOST_AUTO_TEST_CASE( ExpGaussExp )
{
  // Check ExpGaussExp has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew_low = 1.9055;
    double skew_high = 1.2;
    double x0 = mean - 20*sigma;
    double x1 = mean + 20*sigma;
    
    double answer = exp_gauss_exp_integral( mean, sigma, skew_low, skew_high, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 1.0, 1.0E-8 );
  }
  
  // Check fast ExpGaussExp has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew_low = 1.9055;
    double skew_high = 1.2;
    const double amplitude = 1.23456;
    double x0 = mean - 30*sigma;
    double x1 = mean + 30*sigma;
    size_t num_channels = 1024;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    exp_gauss_exp_integral( mean, sigma, amplitude, skew_low, skew_high, &(energies[0]), &(counts[0]), num_channels );
    
    double answer_sum = std::accumulate( begin(counts), end(counts), 0.0 );
    BOOST_CHECK_CLOSE( answer_sum, amplitude, amplitude*1.0E-9 );
  }
  
  
  // Check ExpGaussExp PDF has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double skew_low = 1.9055;
    double skew_high = 1.2;
    const double amplitude = 1.23456;
    double x0 = mean - 40*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 2048;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = exp_gauss_exp_pdf( mean, sigma, skew_low, skew_high, x );
      
      answer_sum += dx * pdf_val;
    }
    
    BOOST_CHECK_CLOSE( answer_sum, 1.0, 1.0E-3 );
  }
  
  {
    double mean = 100;
    double sigma = 1.0;
    double left_skew = 1.05;
    double right_skew = 1.75;
    double prob = 0.000000573303;
    
    pair<double,double> limits = exp_gauss_exp_coverage_limits( mean, sigma, left_skew, right_skew, prob );
    double fraction = 1.0 - exp_gauss_exp_integral( mean, sigma, left_skew, right_skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.317310507863;
    limits = exp_gauss_exp_coverage_limits( mean, sigma, left_skew, right_skew, prob );
    fraction = 1.0 - exp_gauss_exp_integral( mean, sigma, left_skew, right_skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.1;
    limits = exp_gauss_exp_coverage_limits( mean, sigma, left_skew, right_skew, prob );
    fraction = 1.0 - exp_gauss_exp_integral( mean, sigma, left_skew, right_skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.2;
    limits = exp_gauss_exp_coverage_limits( mean, sigma, left_skew, right_skew, prob );
    fraction = 1.0 - exp_gauss_exp_integral( mean, sigma, left_skew, right_skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.4;
    limits = exp_gauss_exp_coverage_limits( mean, sigma, left_skew, right_skew, prob );
    fraction = 1.0 - exp_gauss_exp_integral( mean, sigma, left_skew, right_skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.001;
    limits = exp_gauss_exp_coverage_limits( mean, sigma, left_skew, right_skew, prob );
    fraction = 1.0 - exp_gauss_exp_integral( mean, sigma, left_skew, right_skew, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
  }
  /*
  double exp_gauss_exp_norm( const double sigma, const double skew_left, const double skew_right );

  
  double exp_gauss_exp_left_tail_indefinite( const double mean,
                                          const double sigma,
                                          const double skew_left,
                                          const double skew_right,
                                            const double x );

  double exp_gauss_exp_right_tail_indefinite( const double mean,
                                           const double sigma,
                                           const double skew_left,
                                           const double skew_right,
                                             const double x );

  double exp_gauss_exp_gauss_indefinite( const double mean,
                                      const double sigma,
                                      const double skew_left,
                                      const double skew_right,
                                        const double x );
  
  double exp_gauss_exp_indefinite( const double mean,
                                const double sigma,
                                const double skew_left,
                                const double skew_right,
                                  const double x );
   */
}//BOOST_AUTO_TEST_CASE( ExpGaussExp )


BOOST_AUTO_TEST_CASE( CrystalBall )
{
  // Check CrystalBall has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double alpha = 2.1699;
    double n = 2.5000;
    double x0 = mean - 100*sigma;
    double x1 = mean + 20*sigma;
    
    double answer = crystal_ball_integral( mean, sigma, alpha, n, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 1.0, 5.0E-3 );
  }
  
  // Check fast CrystalBall has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double alpha = 2.1699;
    double n = 2.5000;
    const double amplitude = 1.64;
    double x0 = mean - 100*sigma;
    double x1 = mean + 20*sigma;
    size_t num_channels = 1024; //65536;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    //const double start_cpu = SpecUtils::get_cpu_time();
    
    crystal_ball_integral( mean, sigma, amplitude, alpha, n, &(energies[0]), &(counts[0]), num_channels );
    
    //const double end_cpu = SpecUtils::get_cpu_time();
    //cout << "CrystalBall Took " << (end_cpu - start_cpu) << " seconds" << endl;
    
    double answer_sum = std::accumulate( begin(counts), end(counts), 0.0 );
    BOOST_CHECK_CLOSE( answer_sum/amplitude, 1.0, 5.0E-3 );
    
    for( size_t i = 0; i < num_channels; ++i )
    {
      if( counts[i] > amplitude*1.0E-8 )
      {
        const double val_check = amplitude*PeakDists::crystal_ball_integral( mean, sigma, alpha, n, energies[i], energies[i+1] );
        BOOST_CHECK_CLOSE( counts[i], val_check, 1.0E-6 );
      }
    }
  }
  
  
  // Check CrystalBall PDF has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double alpha = 2.1699;
    double n = 2.5000;
    double x0 = mean - 100*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 2*2*2048;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = crystal_ball_pdf( mean, sigma, alpha, n, x );
      
      answer_sum += dx * pdf_val;
    }
    
    BOOST_CHECK_CLOSE( answer_sum, 1.0, 5.0E-3 );
  }
  
  {
    double mean = 100;
    double sigma = 1.0;
    double alpha = 2.25;
    double n = 8;
    double prob = 0.000000573303;
    
    pair<double,double> limits = crystal_ball_coverage_limits( mean, sigma, alpha, n, prob );
    double fraction = 1.0 - crystal_ball_integral( mean, sigma, alpha, n, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.317310507863;
    limits = crystal_ball_coverage_limits( mean, sigma, alpha, n, prob );
    fraction = 1.0 - crystal_ball_integral( mean, sigma, alpha, n, limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
  }
  
  
  
  /*
  double crystal_ball_norm( const double sigma,
                           const double alpha,
                           const double n );

  double crystal_ball_tail_indefinite_t( const double sigma, const double alpha,
                                        const double n, const double t );
  */
  
}//BOOST_AUTO_TEST_CASE( CrystalBall )


BOOST_AUTO_TEST_CASE( DoubleSidedCrystalBall )
{
  // Check DoubleSidedCrystalBall has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double alpha_low = 2.1699;
    double n_low = 2.5000;
    double alpha_high = 1.8;
    double n_high = 3.5000;
    double x0 = mean - 250*sigma;
    double x1 = mean + 250*sigma;
    
    double answer = double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, x0, x1 );
    BOOST_CHECK_CLOSE( answer, 1.0, 5.0E-3 );
    
    
    // Make sure the indefinite integrals can have boundaries, not just right at their limits
    answer = 0;
    answer += double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, x0, mean-alpha_low*sigma - 3*sigma );
    answer += double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, mean-alpha_low*sigma - 3*sigma, mean-alpha_low*sigma );
    answer += double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, mean-alpha_low*sigma, mean-alpha_low*sigma + sigma );
    answer += double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, mean-alpha_low*sigma + sigma, mean+alpha_high*sigma );
    answer += double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, mean+alpha_high*sigma, mean+alpha_high*sigma + sigma );
    answer += double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, mean+alpha_high*sigma + sigma, x1 );
    BOOST_CHECK_CLOSE( answer, 1.0, 1.0E-2 );
    
    double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    double left_contrib = norm*DSCB_left_tail_indefinite_non_norm_t( alpha_low, n_low, (x0 - mean)/sigma );
    double right_contrib = -norm*DSCB_right_tail_indefinite_non_norm_t( alpha_high, n_high, (x1 - mean)/sigma );
    
    answer = double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, -1000, 3000.0 );
    BOOST_CHECK_CLOSE( answer, 1.0, 5.0E-3 );
  }
  
  // Check fast DoubleSidedCrystalBall has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double alpha_low = 2.1699;
    double n_low = 4.5000;
    double alpha_high = 1.8;
    double n_high = 4.1000;
    const double amplitude = 1.12345;
    double x0 = mean - 30*sigma;
    double x1 = mean + 30*sigma;
    size_t num_channels = 1024;
    vector<double> counts( num_channels, 0.0 );
    vector<float> energies( num_channels + 1, 0.0 );
    for( size_t i = 0; i < energies.size(); ++i )
      energies[i] = x0 + (i*(x1 - x0) / (num_channels + 1));
    
    //const double start_cpu = SpecUtils::get_cpu_time();
    
    double_sided_crystal_ball_integral( mean, sigma, amplitude, alpha_low, n_low, alpha_high, n_high, &(energies[0]), &(counts[0]), num_channels );
    
    //const double end_cpu = SpecUtils::get_cpu_time();
    //cout << "DoubleSidedCrystalBall Took " << (end_cpu - start_cpu) << " seconds" << endl;
    
    double answer_sum = std::accumulate( begin(counts), end(counts), 0.0 );
    
    const double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    answer_sum += norm*DSCB_left_tail_indefinite_non_norm_t( alpha_low, n_low, (x0 - mean)/sigma );
    answer_sum += -norm*DSCB_right_tail_indefinite_non_norm_t( alpha_high, n_high, (x1 - mean)/sigma );
    
    BOOST_CHECK_CLOSE( answer_sum, amplitude, 1.0E-3 );
    
    
    for( size_t i = 0; i < num_channels; ++i )
    {
      if( counts[i] > amplitude*1.0E-8 )
      {
        const double val_check = amplitude*PeakDists::double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, energies[i], energies[i+1] );
        BOOST_CHECK_CLOSE( counts[i], val_check, 1.0E-6 );
      }
    }//for( size_t i = 0; i < num_channels; ++i )
  }
  
  {
    double sigma = (1.9450 / 2.355);
    double alpha_low = 2.1699;
    double n_low = 2.5000;
    double alpha_high = 1.8;
    double n_high = 2.1000;
    
    const double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    double non_norm_sum = DSCB_left_tail_indefinite_non_norm_t(alpha_low, n_low, -alpha_low)
                          + DSCB_gauss_indefinite_non_norm_t(alpha_high)
                          - DSCB_gauss_indefinite_non_norm_t(-alpha_low)
                          - DSCB_right_tail_indefinite_non_norm_t(alpha_high, n_high, alpha_high);
    
    const double mult = non_norm_sum * norm;
    
    BOOST_CHECK_CLOSE( mult, 1.0, 1.0E-6 );
  }
  
  
  /*
   20231109 - the `DSCB_pdf_non_norm` function will fail the test with `alpha_low` and `alpha_high`
   both smaller, with smaller `n` values.  Area comes up to be like 0.82, instead of 1.
   If we  effectively only have either the left or right tail (but not both), then the areas
   are pretty close (which implies its probably in the normalization... but I dont see it).
   Totally unsure why.
   Since I cant figure out why, we just wont use the `DSCB_pdf_non_norm` anywhere...
   
  // Check near-Gaussian DoubleSidedCrystalBall PDF has unit area
  {
    double mean = 100;
    double sigma = 1;
    double alpha_low = 10;
    double n_low = 2.5;
    double alpha_high = 10;
    double n_high = 2.1000;
    double x0 = mean - 10*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 65536;
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = DSCB_pdf_non_norm( mean, sigma, alpha_low, n_low, alpha_high, n_high, x );
      
      answer_sum += dx * pdf_val;
    }
    
    answer_sum += DSCB_left_tail_indefinite_non_norm_t( alpha_low, n_low, (x0 - mean)/sigma );
    answer_sum += -DSCB_right_tail_indefinite_non_norm_t( alpha_high, n_high, (x1 - mean)/sigma );
    
    double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    const double area = norm*answer_sum;
    BOOST_CHECK_CLOSE( area, 1.0, 1.0E-3 );
  }
  
  // Check nearly-only-left-skew DoubleSidedCrystalBall PDF has unit area
  {
    double mean = 100;
    double sigma = 1;
    double alpha_low = 1.5;
    double n_low = 2.5;
    double alpha_high = 10;
    double n_high = 2.1000;
    double x0 = mean - 10*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 65536;
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = DSCB_pdf_non_norm( mean, sigma, alpha_low, n_low, alpha_high, n_high, x );
      
      answer_sum += dx * pdf_val;
    }
    
    answer_sum += DSCB_left_tail_indefinite_non_norm_t( alpha_low, n_low, (x0 - mean)/sigma );
    answer_sum += -DSCB_right_tail_indefinite_non_norm_t( alpha_high, n_high, (x1 - mean)/sigma );
    
    double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    const double area = norm*answer_sum;
    BOOST_CHECK_CLOSE( area, 1.0, 1.0E-6 );
  }
  
  // Check nearly-only-right-skew DoubleSidedCrystalBall PDF has unit area
  {
    double mean = 100;
    double sigma = 1;
    double alpha_low = 10;
    double n_low = 2.5;
    double alpha_high = 1.5;
    double n_high = 2.1000;
    double x0 = mean - 10*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 65536;
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = DSCB_pdf_non_norm( mean, sigma, alpha_low, n_low, alpha_high, n_high, x );
      
      answer_sum += dx * pdf_val;
    }
    
    double below_contrib = DSCB_left_tail_indefinite_non_norm_t( alpha_low, n_low, (x0 - mean)/sigma );
    double above_contrib = -DSCB_right_tail_indefinite_non_norm_t( alpha_high, n_high, (x1 - mean)/sigma );
    
    answer_sum += below_contrib;
    answer_sum += above_contrib;
    
    double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    const double area = norm*answer_sum;
    BOOST_CHECK_CLOSE( area, 1.0, 1.0E-4 );
  }
  
  
  // Check DoubleSidedCrystalBall PDF has unit area
  {
    double mean = 964.4019;
    double sigma = (1.9450 / 2.355);
    double alpha_low = 8.1699;
    double n_low = 2.5000;
    double alpha_high = 9.8;
    double n_high = 2.1000;
    double x0 = mean - 10*sigma;
    double x1 = mean + 10*sigma;
    size_t num_channels = 4096;
    
    double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    
    double answer_sum = 0.0;
    for( size_t i = 0; i < num_channels; ++i )
    {
      const double lower_energy = x0 + (i*(x1 - x0) / (num_channels + 1));
      const double upper_energy = x0 + ((i+1)*(x1 - x0) / (num_channels + 1));
      const double x = 0.5*(lower_energy + upper_energy);
      const double dx = upper_energy - lower_energy;
      
      double pdf_val = DSCB_pdf_non_norm( mean, sigma, alpha_low, n_low, alpha_high, n_high, x );
      const double area = dx * pdf_val;
      
      //double integrated = double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low, alpha_high, n_high, lower_energy, upper_energy );
      //double integrated_non_norm = integrated / norm;
      //const double diff_frac = fabs(integrated_non_norm - area) / std::max(integrated_non_norm,area);
      //if( diff_frac > 1.0E-4 )
      //{
      //  cout << "Starting at t=" << (x-mean)/sigma << ", diff=" << diff_frac << ", area=" << area << ", integrated_non_norm=" << integrated_non_norm << endl;
      //}
      answer_sum += area;
    }
    
    double t_0 = (x0 - mean)/sigma;
    double t_1 = (x1 - mean)/sigma;
    double below_contrib = DSCB_left_tail_indefinite_non_norm_t( alpha_low, n_low, t_0 );
    double above_contrib = -DSCB_right_tail_indefinite_non_norm_t( alpha_high, n_high, t_1 );
    
    answer_sum += below_contrib;
    answer_sum += above_contrib;
    
    
    const double area = norm*answer_sum;
    
    BOOST_CHECK_CLOSE( area, 1.0, 1.0E-3 );
  }
  */
  
  {
    // Check Gaussian component of DSCB is normalized correctly
    double gaus_limit_norm = DSCB_norm( 10, 5, 10, 5 );
    double gause_com = DSCB_gauss_indefinite_non_norm_t( 10 )
                        - DSCB_gauss_indefinite_non_norm_t( -10 );
    
    BOOST_CHECK_CLOSE( gause_com*gaus_limit_norm, 1.0, 1.0E-8 );
    
    
    // Check left and right tail integrals are implemented equivalently
    double left_indef = DSCB_left_tail_indefinite_non_norm_t( 1.5, 5, -1.6 );
    double right_indef = DSCB_right_tail_indefinite_non_norm_t( 1.5, 5, 1.6 );
    
    BOOST_CHECK_CLOSE( left_indef, -right_indef, 1.0E-8 );
    
    // Check that the normalization function is implemented correctly
    const double alpha_low = 1.3;
    const double n_low = 4;
    const double alpha_high = 2.3;
    const double n_high = 6.1;
    
    const double norm = DSCB_norm( alpha_low, n_low, alpha_high, n_high );
    double left_integral = DSCB_left_tail_indefinite_non_norm_t( alpha_low, n_low, -alpha_low );
    double right_integral = -DSCB_right_tail_indefinite_non_norm_t( alpha_high, n_high, alpha_high );
    double mid_integral = DSCB_gauss_indefinite_non_norm_t( alpha_high )
                          - DSCB_gauss_indefinite_non_norm_t( -alpha_low );
    
    BOOST_CHECK_CLOSE( left_integral + right_integral + mid_integral, 1.0/norm, 1.0E-8 );
  }
  
  
  {
    double mean = 100;
    double sigma = 1.5;
    double alpha_low = 1.7;
    double n_low = 8;
    double alpha_high = 2.5;
    double n_high = 20;
    
    double prob = 0.000000573303;
    pair<double,double> limits = double_sided_crystal_ball_coverage_limits( mean, sigma, alpha_low,
                                                                n_low, alpha_high, n_high, prob );
    double fraction = 1.0 - double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low,
                                                               alpha_high, n_high,
                                                               limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    
    prob = 1.0E-3;
    limits = double_sided_crystal_ball_coverage_limits( mean, sigma, alpha_low,
                                                                n_low, alpha_high, n_high, prob );
    fraction = 1.0 - double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low,
                                                        alpha_high, n_high,
                                                        limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 1.0E-2;
    limits = double_sided_crystal_ball_coverage_limits( mean, sigma, alpha_low,
                                                                n_low, alpha_high, n_high, prob );
    fraction = 1.0 - double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low,
                                                        alpha_high, n_high,
                                                        limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 1.0E-1;
    limits = double_sided_crystal_ball_coverage_limits( mean, sigma, alpha_low,
                                                                n_low, alpha_high, n_high, prob );
    fraction = 1.0 - double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low,
                                                        alpha_high, n_high,
                                                        limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.01 );
    
    prob = 0.317310507863;
    limits = double_sided_crystal_ball_coverage_limits( mean, sigma, alpha_low,
                                                                n_low, alpha_high, n_high, prob );
    fraction = 1.0 - double_sided_crystal_ball_integral( mean, sigma, alpha_low, n_low,
                                                        alpha_high, n_high,
                                                        limits.first, limits.second );
    BOOST_CHECK_CLOSE( fraction, prob, 0.1 );
  }
  
  
  /*

    double DSCB_norm( const double alpha_low,
                                          const double n_low,
                                          const double alpha_high,
                                          const double n_high );
    
    double DSCB_left_tail_indefinite_non_norm_t( const double alpha_low,
                                                            const double n_low,
                                                            const double t);
    
    double DSCB_right_tail_indefinite_non_norm_t( const double alpha_high,
                                                             const double n_high,
                                                             const double t);
    
    double DSCB_gauss_indefinite_non_norm_t( const double t);
   */
}//BOOST_AUTO_TEST_CASE( DoubleSidedCrystalBall )
