#ifndef BersteinPolynomial_hpp
#define BersteinPolynomial_hpp
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
#include <vector>
#include <cassert>
#include <stdexcept>
#include <algorithm>


#include <boost/math/special_functions/binomial.hpp>
#include "Eigen/Dense"

/** Berstein polynomials have the nice properties that when fitting them in a non-linear fitter, and you
 want to limit the range of values the function can fit to, you can just limit the range of each Berstein
 coefficient to that range.  They are also more numerically stable, I think.
 */
namespace BersteinPolynomial
{

// Precomputed binomial coefficients for degrees 0-6 (compile-time optimization)
// BINOMIAL_COEFFS[n][k] = C(n,k)
constexpr double BINOMIAL_COEFFS[7][7] = {
  {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},              // n=0: C(0,k)
  {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},              // n=1: C(1,k) 
  {1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0},              // n=2: C(2,k)
  {1.0, 3.0, 3.0, 1.0, 0.0, 0.0, 0.0},              // n=3: C(3,k)
  {1.0, 4.0, 6.0, 4.0, 1.0, 0.0, 0.0},              // n=4: C(4,k)
  {1.0, 5.0, 10.0, 10.0, 5.0, 1.0, 0.0},            // n=5: C(5,k)
  {1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0}            // n=6: C(6,k)
};

// Fast binomial coefficient lookup with fallback to boost for higher orders
template<typename T>
T binomial_coefficient( unsigned int n, unsigned int k )
{
  if( k > n ) return T(0);
  if( k == 0 || k == n ) return T(1);
  
  // Use precomputed table for common small degrees (much faster)
  if( n <= 6 )
  {
    return T( BINOMIAL_COEFFS[n][k] );
  }
  
  // Fallback to boost for higher degrees (rare case)
  return T( boost::math::binomial_coefficient<double>( n, k ) );
}

/** Evaluates a Bernstein polynomial of given degree with the provided coefficients.
 
 The Bernstein polynomial is defined as:
 P(x) = sum_{i=0}^n c_i * B_{i,n}(x)
 
 where B_{i,n}(x) = C(n,i) * x^i * (1-x)^{n-i} are the Bernstein basis functions
 and C(n,i) is the binomial coefficient.
 
 @param x The evaluation point, must be in the range [0, 1]
 @param coefficients Pointer to array of Bernstein coefficients
 @param num_coefficients Number of coefficients, determines polynomial degree
 @return The evaluated polynomial value
 
 Template parameter T can be double, float, or ceres::Jet<> types.
 */
template<typename T>
T evaluate( const T& x, const T * const coefficients, const size_t num_coefficients )
{
  if( num_coefficients == 0 )
    throw std::invalid_argument( "BersteinPolynomial::evaluate: coefficients cannot be empty" );
  
  const unsigned int n = static_cast<unsigned int>(num_coefficients - 1);  // degree
  T result = T(0.0);
  
  // Optimized evaluation using incremental power computation
  // This avoids O(n^2) power calculations while being cache-friendly
  const T one_minus_x = T(1.0) - x;
  
  // For small degrees (common case), use stack allocation to avoid heap overhead
  if( n <= 6 )  // Matches our precomputed binomial table
  {
    T x_powers[7], one_minus_x_powers[7];
    
    x_powers[0] = T(1.0);
    one_minus_x_powers[0] = T(1.0);
    
    for( unsigned int i = 1; i <= n; ++i )
    {
      x_powers[i] = x_powers[i-1] * x;
      one_minus_x_powers[i] = one_minus_x_powers[i-1] * one_minus_x;
    }
    
    for( unsigned int i = 0; i <= n; ++i )
    {
      const T binomial_coeff = binomial_coefficient<T>( n, i );
      const T basis_value = binomial_coeff * x_powers[i] * one_minus_x_powers[n-i];
      result += coefficients[i] * basis_value;
    }
  }
  else
  {
    // For higher degrees, use the original nested loop approach to avoid heap allocation
    // This is rare in practice (comment mentions expected orders up to 5-6)
    for( unsigned int i = 0; i <= n; ++i )
    {
      const T binomial_coeff = binomial_coefficient<T>( n, i );
      
      T x_pow_i = T(1);
      for( unsigned int j = 0; j < i; ++j )
        x_pow_i *= x;
      
      T one_minus_x_pow = T(1);
      for( unsigned int j = 0; j < (n - i); ++j )
        one_minus_x_pow *= one_minus_x;
      
      result += coefficients[i] * binomial_coeff * x_pow_i * one_minus_x_pow;
    }
  }
  
  // TODO: Consider scalar power tracking optimization (x_power *= x, one_minus_x_power /= one_minus_x)
  // but need to handle edge cases like x=0, x=1 more carefully to avoid division by zero
  
  return result;
}


/** Converts power series polynomial coefficients to Bernstein polynomial coefficients.
 
 Converts from power series representation: P(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n
 to Bernstein representation: P(x) = sum_{i=0}^n b_i * B_{i,n}(x)
 
 The input x range [x_min, x_max] is mapped to [0, 1] for the Bernstein representation.
 
 @param power_coeffs Pointer to array of power series coefficients [a_0, a_1, ..., a_n]
 @param num_coefficients Number of power series coefficients
 @param x_min Minimum x value for the power series domain
 @param x_max Maximum x value for the power series domain
 @return Vector of Bernstein coefficients
 */
template<typename T>
std::vector<T> power_series_to_bernstein( const T* power_coeffs,
                                          size_t num_coefficients,
                                          const T& x_min, 
                                          const T& x_max )
{
  if( num_coefficients == 0 )
    throw std::invalid_argument( "BersteinPolynomial::power_series_to_bernstein: coefficients cannot be empty" );
  
  if( x_max <= x_min )
    throw std::invalid_argument( "BersteinPolynomial::power_series_to_bernstein: x_max must be greater than x_min" );
  
  const unsigned int n = static_cast<unsigned int>(num_coefficients - 1);  // degree
  std::vector<T> bernstein_coeffs( n + 1, T(0) );
  
  // First, transform the polynomial coefficients from domain [x_min, x_max] to [0, 1]
  // If P(x) = sum a_j * x^j on [x_min, x_max], we want Q(t) = P(x_min + t*(x_max-x_min)) on [0, 1]
  const T x_range = x_max - x_min;
  std::vector<T> normalized_coeffs( n + 1, T(0) );
  
  // Q(t) = sum a_j * (x_min + t*x_range)^j
  // Expand using binomial theorem: (x_min + t*x_range)^j = sum_{k=0}^j C(j,k) * x_min^{j-k} * (t*x_range)^k
  for( unsigned int j = 0; j <= n; ++j )
  {
    for( unsigned int k = 0; k <= j; ++k )
    {
      const T binomial_j_k = binomial_coefficient<T>( j, k );
      
      // Compute x_min^{j-k}
      T x_min_pow = T(1);
      for( size_t p = 0; p < (j - k); ++p )
        x_min_pow *= x_min;
      
      // Compute x_range^k
      T x_range_pow = T(1);
      for( size_t p = 0; p < k; ++p )
        x_range_pow *= x_range;
      
      normalized_coeffs[k] += power_coeffs[j] * binomial_j_k * x_min_pow * x_range_pow;
    }
  }
  
  // Now convert from power basis to Bernstein basis on [0,1]  
  // Correct transformation matrix formula: b_i = sum_{j=0}^{i} a_j * C(i,j) / C(n,j)
  for( unsigned int i = 0; i <= n; ++i )
  {
    for( unsigned int j = 0; j <= i; ++j )  // Only j <= i terms contribute
    {
      const T binomial_i_j = binomial_coefficient<T>( i, j );
      const T binomial_n_j = binomial_coefficient<T>( n, j );
      
      bernstein_coeffs[i] += normalized_coeffs[j] * binomial_i_j / binomial_n_j;
    }
  }
  
  return bernstein_coeffs;
}


/** Converts Bernstein polynomial coefficients to power series polynomial coefficients.
 
 Converts from Bernstein representation: P(x) = sum_{i=0}^n b_i * B_{i,n}(x) defined on [0,1]
 to power series representation: P(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n defined on [x_min, x_max]
 
 @param bernstein_coeffs Pointer to array of Bernstein coefficients [b_0, b_1, ..., b_n]
 @param num_coefficients Number of Bernstein coefficients
 @param x_min Minimum x value for the output power series domain  
 @param x_max Maximum x value for the output power series domain
 @return Vector of power series coefficients
 */
template<typename T>
std::vector<T> bernstein_to_power_series( const T* bernstein_coeffs,
                                          size_t num_coefficients,
                                          const T& x_min,
                                          const T& x_max )
{
  if( num_coefficients == 0 )
    throw std::invalid_argument( "BersteinPolynomial::bernstein_to_power_series: coefficients cannot be empty" );
  
  if( x_max <= x_min )
    throw std::invalid_argument( "BersteinPolynomial::bernstein_to_power_series: x_max must be greater than x_min" );
  
  const unsigned int n = static_cast<unsigned int>(num_coefficients - 1);  // degree
  
  // First convert Bernstein coefficients to power series coefficients on [0,1]
  // Use the inverse of the transformation matrix M where M[i][j] = C(i,j)/C(n,j) for j<=i
  // The matrix is lower triangular, so use forward substitution
  std::vector<T> normalized_coeffs( n + 1, T(0) );
  
  // Forward substitution to solve the lower triangular system M * a = b
  for( unsigned int j = 0; j <= n; ++j )  // Work forwards from lowest degree
  {
    const T binomial_n_j = binomial_coefficient<T>( n, j );
    T sum = T(0);
    
    // Subtract off the contributions from lower-order terms already solved
    for( unsigned int k = 0; k < j; ++k )
    {
      const T binomial_j_k = binomial_coefficient<T>( j, k );
      const T binomial_n_k = binomial_coefficient<T>( n, k );
      sum += normalized_coeffs[k] * binomial_j_k / binomial_n_k;
    }
    
    // Solve for a_j: b_j = sum + a_j * (diagonal element)
    // The diagonal element M[j][j] = C(j,j)/C(n,j) = 1/C(n,j)
    const T diagonal = T(1) / binomial_n_j;
    normalized_coeffs[j] = (bernstein_coeffs[j] - sum) / diagonal;
  }
  
  // Now transform from domain [0,1] to [x_min, x_max]
  // If Q(t) = sum c_j * t^j on [0,1], we want P(x) = Q((x-x_min)/(x_max-x_min)) on [x_min, x_max]
  const T x_range = x_max - x_min;
  std::vector<T> power_coeffs( n + 1, T(0) );
  
  // P(x) = sum c_j * ((x - x_min) / x_range)^j
  // Expand using binomial theorem: ((x - x_min) / x_range)^j = sum_{k=0}^j C(j,k) * x^k * (-x_min)^{j-k} / x_range^j
  for( unsigned int j = 0; j <= n; ++j )
  {
    // Compute 1 / x_range^j
    T x_range_inv_pow = T(1);
    for( size_t p = 0; p < j; ++p )
      x_range_inv_pow /= x_range;
    
    for( unsigned int k = 0; k <= j; ++k )
    {
      const T binomial_j_k = binomial_coefficient<T>( j, k );
      
      // Compute (-x_min)^{j-k}
      T neg_x_min_pow = T(1);
      const T neg_x_min = -x_min;
      for( size_t p = 0; p < (j - k); ++p )
        neg_x_min_pow *= neg_x_min;
      
      power_coeffs[k] += normalized_coeffs[j] * binomial_j_k * neg_x_min_pow * x_range_inv_pow;
    }
  }
  
  return power_coeffs;
}


/** Evaluates a power series polynomial with domain transformation.
 
 Evaluates P(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n for x in [x_min, x_max],
 but maps the input x from [x_min, x_max] to [0, 1] and then evaluates the equivalent
 Bernstein polynomial.
 
 @param x The evaluation point in [x_min, x_max]
 @param power_coeffs Pointer to array of power series coefficients
 @param num_coefficients Number of power series coefficients
 @param x_min Minimum x value for the power series domain
 @param x_max Maximum x value for the power series domain  
 @return The evaluated polynomial value
 */
template<typename T>
T evaluate_power_series_via_bernstein( const T& x,
                                       const T* power_coeffs,
                                       size_t num_coefficients,
                                       const T& x_min,
                                       const T& x_max )
{
  if( num_coefficients == 0 )
    throw std::invalid_argument( "BersteinPolynomial::evaluate_power_series_via_bernstein: coefficients cannot be empty" );
  
  // Convert to Bernstein representation
  const std::vector<T> bernstein_coeffs = power_series_to_bernstein( power_coeffs, num_coefficients, x_min, x_max );
  
  // Transform x to [0, 1]
  const T x_normalized = (x - x_min) / (x_max - x_min);
  
  // Evaluate Bernstein polynomial
  return evaluate( x_normalized, bernstein_coeffs.data(), bernstein_coeffs.size() );
}


/** Performs linear least squares fitting to data using Bernstein polynomial basis.
 
 Fits a Bernstein polynomial of specified degree to the provided (x,y) data points
 using numerically stable SVD-based linear least squares, similar to RelActCalcAuto.
 Uses Eigen::JacobiSVD (or BDCSVD for older versions) for improved numerical stability
 compared to normal equations.
 
 The fitting minimizes: Σ[(y_i - P(x_i))² / sigma_i²] where P(x) is the Bernstein polynomial.
 
 @param x_values Vector of x data points
 @param y_values Vector of y data points  
 @param uncertainties Vector of uncertainties (standard deviations) for y values
 @param degree Degree of Bernstein polynomial to fit (determines number of coefficients)
 @param x_min Minimum x value for domain mapping to [0,1]
 @param x_max Maximum x value for domain mapping to [0,1] 
 @return Vector of fitted Bernstein coefficients
 @throws std::invalid_argument if input vectors have different sizes or invalid parameters
 @throws std::runtime_error if design matrix is rank deficient or solution fails
 
 Template parameter T can be double or float.
 */
template<typename T>
std::vector<T> fit_bernstein_lls( const std::vector<T>& x_values,
                                  const std::vector<T>& y_values,
                                  const std::vector<T>& uncertainties,
                                  unsigned int degree,
                                  const T& x_min,
                                  const T& x_max )
{
  // Input validation
  if( x_values.empty() || y_values.empty() || uncertainties.empty() )
    throw std::invalid_argument( "BersteinPolynomial::fit_bernstein_lls: input vectors cannot be empty" );
  
  if( x_values.size() != y_values.size() || x_values.size() != uncertainties.size() )
    throw std::invalid_argument( "BersteinPolynomial::fit_bernstein_lls: input vectors must have same size" );
  
  if( x_max <= x_min )
    throw std::invalid_argument( "BersteinPolynomial::fit_bernstein_lls: x_max must be greater than x_min" );
  
  if( x_values.size() < degree + 1 )
    throw std::invalid_argument( "BersteinPolynomial::fit_bernstein_lls: need at least degree+1 data points" );
  
  const size_t n_data = x_values.size();
  const size_t n_coeffs = degree + 1;
  const T x_range = x_max - x_min;
  
  // Build the design matrix A where A[i][j] = B_j,degree(x_i_normalized) / sigma_i
  std::vector<std::vector<T>> A( n_data, std::vector<T>( n_coeffs, T(0) ) );
  std::vector<T> b( n_data );
  
  for( size_t i = 0; i < n_data; ++i )
  {
    // Normalize x to [0,1]
    const T x_norm = (x_values[i] - x_min) / x_range;
    const T inv_sigma = T(1) / uncertainties[i];
    const T one_minus_x = T(1) - x_norm;
    
    // Fill row i of design matrix with weighted Bernstein basis functions
    // Use optimized power computation for common small degrees
    if( degree <= 6 )
    {
      T x_powers[7], one_minus_x_powers[7];
      
      x_powers[0] = T(1);
      one_minus_x_powers[0] = T(1);
      
      for( unsigned int k = 1; k <= degree; ++k )
      {
        x_powers[k] = x_powers[k-1] * x_norm;
        one_minus_x_powers[k] = one_minus_x_powers[k-1] * one_minus_x;
      }
      
      for( unsigned int j = 0; j <= degree; ++j )
      {
        const T binomial_coeff = binomial_coefficient<T>( degree, j );
        const T basis_value = binomial_coeff * x_powers[j] * one_minus_x_powers[degree-j];
        A[i][j] = basis_value * inv_sigma;
      }
    }
    else
    {
      // Fallback for higher degrees (rare case)
      for( unsigned int j = 0; j <= degree; ++j )
      {
        const T binomial_coeff = binomial_coefficient<T>( degree, j );
        
        T x_pow_j = T(1);
        for( unsigned int k = 0; k < j; ++k )
          x_pow_j *= x_norm;
        
        T one_minus_x_pow = T(1);
        for( unsigned int k = 0; k < (degree - j); ++k )
          one_minus_x_pow *= one_minus_x;
        
        A[i][j] = binomial_coeff * x_pow_j * one_minus_x_pow * inv_sigma;
      }
    }
    
    // Weighted y value
    b[i] = y_values[i] * inv_sigma;
  }
  
  // Build Eigen design matrix A and weighted y vector b
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigen_A( n_data, n_coeffs );
  Eigen::Matrix<T, Eigen::Dynamic, 1> eigen_b( n_data );
  
  for( size_t i = 0; i < n_data; ++i )
  {
    for( size_t j = 0; j < n_coeffs; ++j )
      eigen_A(i, j) = A[i][j];
    eigen_b(i) = b[i];
  }
  
  // Use SVD for numerically stable solution: A * c = b
  // This is more stable than normal equations (A^T A) * c = A^T b
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(
    eigen_A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#else
  Eigen::BDCSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(
    eigen_A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  // Check conditioning using singular values
  Eigen::Matrix<T, Eigen::Dynamic, 1> singular_values = svd.singularValues();
  const T threshold = T(1.0e-12) * singular_values.maxCoeff(); // Threshold for numerical singularity
  
  size_t rank = 0;
  for( int i = 0; i < singular_values.size(); ++i )
    rank += (singular_values(i) > threshold);
  
  if( rank < n_coeffs )
    throw std::runtime_error( "BersteinPolynomial::fit_bernstein_lls: design matrix is rank deficient (rank=" 
                             + std::to_string(rank) + ", expected=" + std::to_string(n_coeffs) + ")" );
  
  // Solve using SVD
  Eigen::Matrix<T, Eigen::Dynamic, 1> eigen_coeffs = svd.solve( eigen_b );
  
  // Convert back to std::vector
  std::vector<T> coeffs( n_coeffs );
  for( size_t i = 0; i < n_coeffs; ++i )
    coeffs[i] = eigen_coeffs(i);
  
  return coeffs;
}


// Convenience overloads for backward compatibility with std::vector

/** Convenience overload for evaluate() that accepts std::vector */
template<typename T>
T evaluate( const T& x, const std::vector<T>& coefficients )
{
  return evaluate( x, coefficients.data(), coefficients.size() );
}

/** Convenience overload for power_series_to_bernstein() that accepts std::vector */
template<typename T>
std::vector<T> power_series_to_bernstein( const std::vector<T>& power_coeffs, 
                                          const T& x_min, 
                                          const T& x_max )
{
  return power_series_to_bernstein( power_coeffs.data(), power_coeffs.size(), x_min, x_max );
}

/** Convenience overload for bernstein_to_power_series() that accepts std::vector */
template<typename T>
std::vector<T> bernstein_to_power_series( const std::vector<T>& bernstein_coeffs,
                                          const T& x_min,
                                          const T& x_max )
{
  return bernstein_to_power_series( bernstein_coeffs.data(), bernstein_coeffs.size(), x_min, x_max );
}

/** Convenience overload for evaluate_power_series_via_bernstein() that accepts std::vector */
template<typename T>
T evaluate_power_series_via_bernstein( const T& x,
                                       const std::vector<T>& power_coeffs,
                                       const T& x_min,
                                       const T& x_max )
{
  return evaluate_power_series_via_bernstein( x, power_coeffs.data(), power_coeffs.size(), x_min, x_max );
}

}  // namespace BersteinPolynomial

#endif  // BersteinPolynomial_hpp
