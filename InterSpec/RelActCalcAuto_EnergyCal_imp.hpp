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

#ifndef RelActCalcAuto_EnergyCal_imp_hpp
#define RelActCalcAuto_EnergyCal_imp_hpp

#include "InterSpec_config.h"

#include <cmath>
#include <vector>
#include <utility>
#include <cassert>
#include <algorithm>
#include <type_traits>

#include <Eigen/Dense>

#if( PERFORM_DEVELOPER_CHECKS )
#include "SpecUtils/SpecFile.h"
#endif

namespace ceres
{
  /* dummy namespace for when using this file for only doubles, and ceres.h hasnt been included */
}

namespace RelActCalcAutoImp
{

/** Node in cubic spline interpolation, templated for automatic differentiation.

 Represents a single interval in a cubic spline with polynomial coefficients.
 The template parameter T allows for use with automatic differentiation types
 (e.g., ceres::Jet) where derivatives need to be tracked.

 The polynomial form is: f(h) = a*h^3 + b*h^2 + c*h + y, where h = x - node.x
 */
template<typename T>
struct CubicSplineNodeT
{
  double x;  // Anchor x-value (energy) - always double since these are fixed
  T y;       // Value at anchor (can be Jet type for auto-diff)
  T a, b, c; // Polynomial coefficients: f(h) = a*h^3 + b*h^2 + c*h + y, where h = x - node.x
};


/** Create a natural cubic spline from deviation pair anchors and their offsets.

 The spline uses natural boundary conditions (second derivative = 0 at endpoints).
 This is compatible with SpecUtils deviation pair conventions.

 @param deviation_pairs Vector of (energy, offset) pairs, sorted by energy
 @returns Vector of spline nodes for interpolation
 */
template<typename T>
std::vector<CubicSplineNodeT<T>> create_cubic_spline(
  const std::vector<std::pair<double,T>> &deviation_pairs )
{
  const size_t n = deviation_pairs.size();

  if( n < 2 )
    return std::vector<CubicSplineNodeT<T>>{};

  // Verify sorted
#ifndef NDEBUG
  for( size_t i = 1; i < n; ++i )
    assert( deviation_pairs[i].first > deviation_pairs[i-1].first );
#endif

  // Extract anchor energies and offsets
  std::vector<double> anchor_energies( n );
  std::vector<T> offsets( n );
  for( size_t i = 0; i < n; ++i )
  {
    anchor_energies[i] = deviation_pairs[i].first;
    offsets[i] = deviation_pairs[i].second;
  }

  // CRITICAL: Transform x-coordinates to match SpecUtils convention
  // SpecUtils::create_cubic_spline_for_dev_pairs subtracts the offset from the energy
  // (see EnergyCalibration.cpp line 265: offsets[i].first -= offsets[i].second)
  // This is because deviation pairs represent: at true_energy, you observe it at (true_energy - offset)
  std::vector<double> transformed_energies( n );
  for( size_t i = 0; i < n; ++i )
  {
    // Extract the scalar value from T (could be double or ceres::Jet)
    double offset_val;
    if constexpr ( std::is_same_v<T, double> )
      offset_val = offsets[i];
    else
      offset_val = offsets[i].a;  // For ceres::Jet, .a is the scalar value

    transformed_energies[i] = anchor_energies[i] - offset_val;
  }

  // Build the tri-diagonal system to solve for second derivatives (b values)
  // Using natural spline boundary conditions (second derivative = 0 at ends)
  //
  // The matrix elements depend only on x-coordinates (doubles)
  // But the RHS and solution depend on y-coordinates (template type T)
  //
  // We use Eigen to build and invert the tri-diagonal matrix (all doubles),
  // then manually multiply the inverse by the RHS vector (type T) to get the solution.

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero( n, n );
  std::vector<T> rhs( n );

  // Interior points (use transformed_energies for spacing)
  for( size_t i = 1; i < n - 1; ++i )
  {
    const double h_prev = transformed_energies[i] - transformed_energies[i-1];
    const double h_next = transformed_energies[i+1] - transformed_energies[i];

    A( i, i - 1 ) = h_prev / 3.0;             // lower diagonal
    A( i, i )     = (h_prev + h_next) / 1.5;  // diagonal
    A( i, i + 1 ) = h_next / 3.0;             // upper diagonal

    rhs[i] = (offsets[i+1] - offsets[i]) / h_next - (offsets[i] - offsets[i-1]) / h_prev;
  }

  // Left boundary: natural spline (second derivative = 0)
  A( 0, 0 ) = 2.0;
  rhs[0] = T( 0.0 );

  // Right boundary: zero first derivative (matches SpecUtils deviation pair convention)
  // Condition: f'(x_n) = 0
  // This gives: 2*b_n + b_{n-1} = -6*(y_n - y_{n-1})/h_{n-1}^2
  const double h_last = transformed_energies[n-1] - transformed_energies[n-2];
  A( n - 1, n - 2 ) = h_last / 3.0;
  A( n - 1, n - 1 ) = 2.0 * h_last / 3.0;
  rhs[n - 1] = -(offsets[n-1] - offsets[n-2]) / h_last;

  // Compute the inverse of A using Eigen's partial-pivot LU decomposition
  const Eigen::MatrixXd A_inv = A.partialPivLu().inverse();

  // Multiply A_inv (double) by rhs (type T) to get b_vals (type T)
  // We do this manually to support ceres::Jet<> types
  std::vector<T> b_vals( n, T(0.0) );
  for( size_t i = 0; i < n; ++i )
  {
    for( size_t j = 0; j < n; ++j )
      b_vals[i] += A_inv( i, j ) * rhs[j];
  }

  // Calculate spline coefficients a, b, c for each interval
  // IMPORTANT: Use transformed_energies for the node x-coordinates to match SpecUtils
  std::vector<CubicSplineNodeT<T>> nodes( n );

  for( size_t i = 0; i < n - 1; ++i )
  {
    const double h = transformed_energies[i+1] - transformed_energies[i];

    nodes[i].x = transformed_energies[i];
    nodes[i].y = offsets[i];
    nodes[i].a = (b_vals[i+1] - b_vals[i]) / (3.0 * h);
    nodes[i].b = b_vals[i];
    nodes[i].c = (offsets[i+1] - offsets[i]) / h - (2.0 * b_vals[i] + b_vals[i+1]) * h / 3.0;
  }

  // Last node (for extrapolation - just store position, coefficients are zero)
  nodes[n-1].x = transformed_energies[n-1];
  nodes[n-1].y = offsets[n-1];
  nodes[n-1].a = T(0.0);
  nodes[n-1].b = T(0.0);
  nodes[n-1].c = T(0.0);

  return nodes;
}//create_cubic_spline(...)


/** Evaluate the templated cubic spline at a given energy.

 @param energy The energy to evaluate at (template type T)
 @param nodes The spline nodes from create_cubic_spline
 @returns The interpolated value (template type T)
 */
template<typename T>
T eval_cubic_spline( const T energy, const std::vector<CubicSplineNodeT<T>> &nodes )
{

  double lookup_energy;
  if constexpr ( std::is_same_v<T, double> )
    lookup_energy = energy;
  else
    lookup_energy = energy.a;  // Extract constant part from Jet


  if( nodes.empty() )
    return T(0.0);

  // Find the interval using binary search on x values
  const auto it = std::upper_bound( nodes.begin(), nodes.end(), lookup_energy,
    []( double e, const CubicSplineNodeT<T> &node ) { return e < node.x; } );

  // Clamp to first/last values for extrapolation
  if( it == nodes.begin() )
    return nodes.front().y;

  if( it == nodes.end() )
    return nodes.back().y;

  // Evaluate the cubic polynomial in the interval
  const CubicSplineNodeT<T> &node = *(it - 1);
  const T h = energy - node.x;

  // f(h) = ((a*h + b)*h + c)*h + y
  return ((node.a * h + node.b) * h + node.c) * h + node.y;
}//eval_cubic_spline(...)


/** Compute the inverse correction for deviation pairs (analogous to SpecUtils::correction_due_to_dev_pairs).

 Given a true_energy, this function solves for the offset that was applied:
   true_energy = observed_energy + deviation_correction(observed_energy)
   true_energy = (true_energy - answer) + deviation_correction(true_energy - answer)

 It uses an iterative Newton-Raphson style solver to find 'answer' such that:
   answer = deviation_correction(true_energy - answer)

 @param true_energy The true energy value
 @param nodes The spline nodes from create_cubic_spline
 @returns The offset correction that was applied
 */
template<typename T>
T correction_due_to_deviation_pairs( const T true_energy, const std::vector<CubicSplineNodeT<T>> &nodes )
{
  if( nodes.empty() )
    return T(0.0);

  // Initial guess: directly evaluate the spline at true_energy
  // This is what the inverse spline would give, but we'll refine it iteratively
  T answer = eval_cubic_spline( true_energy, nodes );

  // Verify by computing: check = deviation_correction(true_energy - answer)
  // We want: answer == check
  T check = eval_cubic_spline( true_energy - answer, nodes );
  T diff = answer - check;

  const double tolerance = 0.0001;  // 0.1 eV tolerance, same as SpecUtils

  if( abs(diff) < tolerance )
    return answer;

  // Iterate to refine the answer
  const int max_iterations = 15;
  int niters = 0;

  while( abs(diff) > tolerance && niters < max_iterations )
  {
    // Newton-Raphson update: answer -= diff
    answer = answer - diff;
    check = eval_cubic_spline( true_energy - answer, nodes );
    diff = answer - check;

    ++niters;
  }

#if( PERFORM_DEVELOPER_CHECKS )
  if( niters >= max_iterations )
  {
    double true_energy_val, diff_val;
    if constexpr ( std::is_same_v<T, double> )
    {
      diff_val = diff;
      true_energy_val = true_energy;
    }else
    {
      diff_val = diff.a;
      true_energy_val = true_energy.a;
    }

    char buffer[512];
    snprintf( buffer, sizeof(buffer),
      "correction_due_to_deviation_pairs: Failed to converge after %d iterations for true_energy=%.6f keV, final diff=%.6g keV",
      max_iterations, true_energy_val, diff_val );
    log_developer_error( __func__, buffer );
  }
#endif

  return answer;
}//correction_due_to_deviation_pairs(...)


/** Helper function to evaluate polynomial energy calibration at a given channel.

 @param channel The channel number (can be fractional)
 @param coeffs The polynomial coefficients
 @param deviation_pairs The deviation pairs (energy, offset)
 @returns Energy in keV at the given channel
 */
template<typename T>
T polynomial_energy( const T channel,
                               const std::vector<T> &coeffs,
                               const std::vector<CubicSplineNodeT<T>> &dev_pair_spline )
{
  T val = T(0.0);
  T channel_power = T(1.0);

  for( size_t i = 0; i < coeffs.size(); ++i )
  {
    val += coeffs[i] * channel_power;
    channel_power *= channel;
  }

  if( !dev_pair_spline.empty() )
    val += eval_cubic_spline( val, dev_pair_spline );

  return val;
}//polynomial_energy(...)


/** Helper function to evaluate full range fraction energy calibration at a given bin.

 For FRF calibration: E = C₀ + x*C₁ + x²*C₂ + x³*C₃ + C₄/(1+60*x) where x = bin/nchannel

 @param bin The bin number (can be fractional)
 @param coeffs The FRF coefficients (typically 2-5 coefficients)
 @param nchannel Total number of channels
 @param dev_pair_spline The deviation pair spline nodes
 @returns Energy in keV at the given bin
 */
template<typename T>
T fullrangefraction_energy( const T bin,
                                      const std::vector<T> &coeffs,
                                      const size_t nchannel,
                                      const std::vector<CubicSplineNodeT<T>> &dev_pair_spline )
{
  const T x = bin / static_cast<double>(nchannel);
  const size_t ncoeffs = std::min( coeffs.size(), size_t(4) );

  T val = T(0.0);
  T x_power = T(1.0);

  for( size_t c = 0; c < ncoeffs; ++c )
  {
    val += coeffs[c] * x_power;
    x_power *= x;
  }

  if( coeffs.size() > 4 )
    val += coeffs[4] / (T(1.0) + T(60.0)*x);

  if( !dev_pair_spline.empty() )
    val += eval_cubic_spline( val, dev_pair_spline );

  return val;
}//fullrangefraction_energy(...)


/** Find the channel corresponding to a given energy for polynomial calibration.

 Templated version of SpecUtils::find_polynomial_channel that supports automatic differentiation.
 Uses analytical solutions for low-order polynomials (order <= 2) and binary search for higher orders.

 TODO: when there are deviation pairs - the keeping derivatives around in `ceres::Jet<>` is kidna a hack; it seems to work, but it could likely be done better in some way
 
 @param energy The energy in keV to find the channel for
 @param coeffs The polynomial coefficients (E = C₀ + C₁*ch + C₂*ch² + ...)
 @param nchannel Number of channels in the spectrum
 @param deviation_pairs The deviation pair anchors and offsets
 @param accuracy The accuracy in keV for the binary search (default 0.001 keV)
 @returns The channel number (fractional) corresponding to the energy
 */
template<typename T>
T find_polynomial_channel( const T energy,
                           const std::vector<T> &coeffs,
                           const size_t nchannel,
                           const std::vector<std::pair<double,T>> &deviation_pairs,
                           const double accuracy = 0.001 )
{
  using namespace std;
  using namespace ceres;

  if( coeffs.empty() || nchannel < 2 )
    return T(0.0);

  // Count non-zero coefficients
  size_t ncoefs = 0;
  for( size_t i = 0; i < coeffs.size(); ++i )
  {
    double coef_val;
    if constexpr ( std::is_same_v<T, double> )
      coef_val = coeffs[i];
    else
      coef_val = coeffs[i].a;

    if( coef_val != 0.0 )
      ncoefs = i + 1;
  }

  if( ncoefs < 2 )
    return T(0.0);

  // Build deviation pair spline
  const std::vector<CubicSplineNodeT<T>> dev_pair_spline =
    deviation_pairs.empty() ? std::vector<CubicSplineNodeT<T>>{} : create_cubic_spline( deviation_pairs );

  // For low-order polynomials, use analytical solutions
  if( ncoefs < 4 && dev_pair_spline.empty() )
  {
    // Linear case
    if( ncoefs == 2 )
    {
      return (energy - coeffs[0]) / coeffs[1];
    }

    // Quadratic case: solve C₀ + C₁*ch + C₂*ch² = energy
    if( ncoefs == 3 )
    {
      const T a = coeffs[2];
      const T b = coeffs[1];
      const T c = coeffs[0] - energy;

      const T discriminant = b*b - T(4.0)*a*c;

      if( discriminant < 0.0 )
        return T(0.0);  // No real solution

      const T sqrt_disc = sqrt(discriminant);
      const T root1 = (-b + sqrt_disc) / (T(2.0) * a);
      const T root2 = (-b - sqrt_disc) / (T(2.0) * a);

      // Extract scalar values for comparison
      double root1_val, root2_val;
      if constexpr ( std::is_same_v<T, double> )
      {
        root1_val = root1;
        root2_val = root2;
      }
      else
      {
        root1_val = root1.a;
        root2_val = root2.a;
      }

      // Choose the root that's in valid range
      const bool root1_valid = (root1_val >= 0.0) && (root1_val < static_cast<double>(nchannel));
      const bool root2_valid = (root2_val >= 0.0) && (root2_val < static_cast<double>(nchannel));

      if( root1_valid && !root2_valid )
        return root1;
      if( root2_valid && !root1_valid )
        return root2;

      // Both valid - choose the one closer to linear solution
      if( root1_valid && root2_valid )
      {
        const T linear_sol = (energy - coeffs[0]) / coeffs[1];
        double linear_val;
        if constexpr ( std::is_same_v<T, double> )
          linear_val = linear_sol;
        else
          linear_val = linear_sol.a;

        const double dist1 = std::fabs( root1_val - linear_val );
        const double dist2 = std::fabs( root2_val - linear_val );

        return (dist1 < dist2) ? root1 : root2;
      }

      return T(0.0);  // Neither valid
    }
  }

  // For higher-order polynomials or with deviation pairs, use binary search

  // Find the energy range
  const T e_low = polynomial_energy( T(0.0), coeffs, dev_pair_spline );
  const T e_high = polynomial_energy( T(static_cast<double>(nchannel - 1)), coeffs, dev_pair_spline );


  // Check if energy is in range
  if( energy < (min)(e_low, e_high) || energy > (max)(e_low, e_high) )
  {
    // Extrapolate linearly
    if( energy < (min)(e_low, e_high) )
    {
      T answer = e_low;
      if constexpr ( !std::is_same_v<T, double> )
        answer.a = 0.0;
      else 
        answer = 0.0;
      return answer;
    }else
    {
      T answer = e_high;
      if constexpr ( !std::is_same_v<T, double> )
        answer.a = static_cast<double>(nchannel - 1);
      else 
        answer = static_cast<double>(nchannel - 1);
      return answer;
    }
  }

  // Binary search (and derivative-preserving implicit update for Jet types)
  if constexpr ( std::is_same_v<T, double> )
  {
    T low_ch(0.0);
    T high_ch( static_cast<double>(nchannel - 1) );

    const int max_iterations = 1000;
    int niter = 0;

    while( (high_ch - low_ch) > 1.0e-6 && niter < max_iterations )
    {
      const T mid_ch = 0.5 * (low_ch + high_ch);
      const T e_mid = polynomial_energy( T(mid_ch), coeffs, dev_pair_spline );

      if( abs(e_mid - energy) < accuracy )
        return T(mid_ch);

      if( (e_mid < energy && e_low < e_high) ||
          (e_mid > energy && e_low > e_high) )
      {
        low_ch = mid_ch;
      }
      else
      {
        high_ch = mid_ch;
      }

      ++niter;
    }

    return T(0.5 * (low_ch + high_ch));
  }else
  {
    // Scalar binary search for the channel value.
    std::vector<double> coeffs_scalar( coeffs.size(), 0.0 );
    for( size_t i = 0; i < coeffs.size(); ++i )
      coeffs_scalar[i] = coeffs[i].a;

    std::vector<std::pair<double,double>> dev_pairs_scalar;
    dev_pairs_scalar.reserve( deviation_pairs.size() );
    for( const auto &pair : deviation_pairs )
      dev_pairs_scalar.emplace_back( pair.first, pair.second.a );

    const std::vector<CubicSplineNodeT<double>> dev_pair_spline_scalar =
      dev_pairs_scalar.empty() ? std::vector<CubicSplineNodeT<double>>{} : create_cubic_spline( dev_pairs_scalar );

    auto energy_at = [&coeffs_scalar,&dev_pair_spline_scalar]( const double ch ) -> double {
      return polynomial_energy( ch, coeffs_scalar, dev_pair_spline_scalar );
    };

    const double e_low_val = energy_at( 0.0 );
    const double e_high_val = energy_at( static_cast<double>(nchannel - 1) );
    const double energy_val = energy.a;

    double low = 0.0;
    double high = static_cast<double>(nchannel - 1);
    const int max_iterations = 1000;
    int niter = 0;
    while( ((high - low) > 1.0e-6) && (niter < max_iterations) )
    {
      const double mid = 0.5 * (low + high);
      const double e_mid_val = energy_at( mid );
      if( std::fabs( e_mid_val - energy_val ) < accuracy )
      {
        low = mid;
        high = mid;
        break;
      }

      if( ((e_mid_val < energy_val) && (e_low_val < e_high_val))
          || ((e_mid_val > energy_val) && (e_low_val > e_high_val)) )
      {
        low = mid;
      }else
      {
        high = mid;
      }

      ++niter;
    }

    double ch_val = 0.5 * (low + high);

    const double eps = 1.0e-3;
    double df_dch = (energy_at( ch_val + eps ) - energy_at( ch_val - eps )) / (2.0 * eps);


#if( PERFORM_DEVELOPER_CHECKS )
    static int s_logged_small_df_dch = 0;
    if( (std::fabs( df_dch ) < 1.0e-14) && (s_logged_small_df_dch < 50) )
    {
      char buffer[256];
      snprintf( buffer, sizeof(buffer),
        "find_polynomial_channel: df_dch small (%.3e) at ch=%.6f energy=%.6f",
        df_dch, ch_val, energy.a );
      log_developer_error( __func__, buffer );
      s_logged_small_df_dch += 1;
    }
#endif

    const T f = polynomial_energy( T(ch_val), coeffs, dev_pair_spline ) - energy;
    T answer = T( ch_val );
    if( std::fabs( df_dch ) > 1.0e-14 )
    {
      for( size_t i = 0; i < static_cast<size_t>(T::DIMENSION); ++i )
        answer.v[i] = -f.v[i] / df_dch;
    }
    return answer;
  }
}//find_polynomial_channel(...)


/** Find the channel corresponding to a given energy for full range fraction calibration.

 Templated version of SpecUtils::find_fullrangefraction_channel that supports automatic differentiation.
 Uses analytical solutions for low-order equations and binary search for higher orders.

 TODO: when there are deviation pairs - the keeping derivatives around in `ceres::Jet<>` is kidna a hack; it seems to work, but it could likely be done better in some way
 
 @param energy The energy in keV to find the channel for
 @param coeffs The FRF coefficients
 @param nchannel Number of channels in the spectrum
 @param deviation_pairs The deviation pair anchors and offsets
 @param accuracy The accuracy in keV for the binary search (default 0.001 keV)
 @returns The channel number (fractional) corresponding to the energy
 */
template<typename T>
T find_fullrangefraction_channel( const T energy,
                                  const std::vector<T> &coeffs,
                                  const size_t nchannel,
                                  const std::vector<std::pair<double,T>> &deviation_pairs,
                                  const double accuracy = 0.001 )
{
  using namespace std;
  using namespace ceres;

  if( coeffs.empty() || nchannel < 2 )
    return T(0.0);

  // Count non-zero coefficients
  size_t ncoefs = 0;
  for( size_t i = 0; i < std::min(coeffs.size(), size_t(4)); ++i )
  {
    double coef_val;
    if constexpr ( std::is_same_v<T, double> )
      coef_val = coeffs[i];
    else
      coef_val = coeffs[i].a;

    if( coef_val != 0.0 )
      ncoefs = i + 1;
  }

  if( coeffs.size() > 4 )
    ncoefs = 5;  // Has the 1/(1+60x) term

  if( ncoefs < 2 )
    return T(0.0);

  // Build deviation pair spline
  const std::vector<CubicSplineNodeT<T>> dev_pair_spline =
    deviation_pairs.empty() ? std::vector<CubicSplineNodeT<T>>{} : create_cubic_spline( deviation_pairs );

  // For low-order FRF without deviation pairs, use analytical solutions
  if( ncoefs < 4 && dev_pair_spline.empty() )
  {
    // Linear case: E = C₀ + x*C₁, where x = bin/nchannel
    // Solve: bin = nchannel * (E - C₀) / C₁
    if( ncoefs == 2 )
    {
      return T(static_cast<double>(nchannel)) * (energy - coeffs[0]) / coeffs[1];
    }

    // Quadratic case: E = C₀ + x*C₁ + x²*C₂
    // Let y = bin/nchannel, solve: C₂*y² + C₁*y + (C₀ - E) = 0
    if( ncoefs == 3 )
    {
      const T a = coeffs[2];
      const T b = coeffs[1];
      const T c = coeffs[0] - energy;

      const T discriminant = b*b - T(4.0)*a*c;

      if( discriminant < 0.0 )
        return T(0.0);

      const T sqrt_disc = sqrt(discriminant);
      const T root1 = (-b + sqrt_disc) / (T(2.0) * a);
      const T root2 = (-b - sqrt_disc) / (T(2.0) * a);

      // Convert from x (fractional position) to bin number
      const T bin1 = root1 * T(static_cast<double>(nchannel));
      const T bin2 = root2 * T(static_cast<double>(nchannel));

      // Choose the root in valid range
      const bool bin1_valid = (bin1 >= 0.0) && (bin1 < static_cast<double>(nchannel));
      const bool bin2_valid = (bin2 >= 0.0) && (bin2 < static_cast<double>(nchannel));

      if( bin1_valid && !bin2_valid )
        return bin1;
      if( bin2_valid && !bin1_valid )
        return bin2;

      // Both valid - choose closer to linear solution
      if( bin1_valid && bin2_valid )
      {
        const T linear_sol = T(static_cast<double>(nchannel)) * (energy - coeffs[0]) / coeffs[1];
        const T dist1 = abs( bin1 - linear_sol );
        const T dist2 = abs( bin2 - linear_sol );

        return (dist1 < dist2) ? bin1 : bin2;
      }

      return T(0.0);
    }
  }

  // For higher-order or with deviation pairs, use binary search
  const T e_low = fullrangefraction_energy( T(0.0), coeffs, nchannel, dev_pair_spline );
  const T e_high = fullrangefraction_energy( T(static_cast<double>(nchannel - 1)), coeffs, nchannel, dev_pair_spline );

  // Check if energy is in range
  if( (energy < (min)(e_low, e_high)) || (energy > (max)(e_low, e_high)) )
  {
    if( energy < min(e_low, e_high) )
    {
      T answer = e_low;
      if constexpr ( !std::is_same_v<T, double> )
        answer.a = 0.0;
      else 
        answer = 0.0;
    return answer;
    }else
    {
      T answer = e_high;
      if constexpr ( !std::is_same_v<T, double> )
        answer.a = static_cast<double>(nchannel - 1);
      else 
        answer = static_cast<double>(nchannel - 1);

      return answer;
    }
  }

  // Binary search (and derivative-preserving implicit update for Jet types)
  if constexpr ( std::is_same_v<T, double> )
  {
    T low_bin( 0.0 );
    T high_bin( static_cast<double>(nchannel - 1) );

    const int max_iterations = 1000;
    int niter = 0;

    while( ((high_bin - low_bin) > 1.0e-6) && niter < max_iterations )
    {
      const T mid_bin = 0.5 * (low_bin + high_bin);
      const T e_mid = fullrangefraction_energy( mid_bin, coeffs, nchannel, dev_pair_spline );

      if( abs(e_mid - energy) < accuracy )
        return T(mid_bin);

      if( ((e_mid < energy) && (e_low < e_high))
          || ((e_mid > energy) && (e_low > e_high)) )
      {
        low_bin = mid_bin;
      }
      else
      {
        high_bin = mid_bin;
      }

      ++niter;
    }

    return T(0.5 * (low_bin + high_bin));
  }else
  {
    std::vector<double> coeffs_scalar( coeffs.size(), 0.0 );
    for( size_t i = 0; i < coeffs.size(); ++i )
      coeffs_scalar[i] = coeffs[i].a;

    std::vector<std::pair<double,double>> dev_pairs_scalar;
    dev_pairs_scalar.reserve( deviation_pairs.size() );
    for( const auto &pair : deviation_pairs )
      dev_pairs_scalar.emplace_back( pair.first, pair.second.a );

    const std::vector<CubicSplineNodeT<double>> dev_pair_spline_scalar =
      dev_pairs_scalar.empty() ? std::vector<CubicSplineNodeT<double>>{} : create_cubic_spline( dev_pairs_scalar );

    auto energy_at = [&coeffs_scalar,nchannel,&dev_pair_spline_scalar]( const double bin ) -> double {
      return fullrangefraction_energy( bin, coeffs_scalar, nchannel, dev_pair_spline_scalar );
    };

    const double e_low_val = energy_at( 0.0 );
    const double e_high_val = energy_at( static_cast<double>(nchannel - 1) );
    const double energy_val = energy.a;

    double low = 0.0;
    double high = static_cast<double>(nchannel - 1);
    const int max_iterations = 1000;
    int niter = 0;
    while( ((high - low) > 1.0e-6) && (niter < max_iterations) )
    {
      const double mid = 0.5 * (low + high);
      const double e_mid_val = energy_at( mid );
      if( std::fabs( e_mid_val - energy_val ) < accuracy )
      {
        low = mid;
        high = mid;
        break;
      }

      if( ((e_mid_val < energy_val) && (e_low_val < e_high_val))
          || ((e_mid_val > energy_val) && (e_low_val > e_high_val)) )
      {
        low = mid;
      }else
      {
        high = mid;
      }

      ++niter;
    }

    double bin_val = 0.5 * (low + high);

    const double eps = 1.0e-3;
    double df_dch = (energy_at( bin_val + eps ) - energy_at( bin_val - eps )) / (2.0 * eps);


#if( PERFORM_DEVELOPER_CHECKS )
    static int s_logged_small_df_dch = 0;
    if( (std::fabs( df_dch ) < 1.0e-14) && (s_logged_small_df_dch < 50) )
    {
      char buffer[256];
      snprintf( buffer, sizeof(buffer),
        "find_fullrangefraction_channel: df_dch small (%.3e) at bin=%.6f energy=%.6f",
        df_dch, bin_val, energy.a );
      log_developer_error( __func__, buffer );
      s_logged_small_df_dch += 1;
    }
#endif

    const T f = fullrangefraction_energy( T(bin_val), coeffs, nchannel, dev_pair_spline ) - energy;
    T answer = T( bin_val );
    if( std::fabs( df_dch ) > 1.0e-14 )
    {
      for( size_t i = 0; i < static_cast<size_t>(T::DIMENSION); ++i )
        answer.v[i] = -f.v[i] / df_dch;
    }
    return answer;
  }
}//find_fullrangefraction_channel(...)


/** Find the channel corresponding to a given energy for lower channel edge calibration.

 For lower channel edge calibration, the energies of each channel's lower edge are explicitly specified.
 This function performs binary search to find which channel contains the given energy.

 @param energy The energy in keV to find the channel for
 @param channel_energies The lower energies for each channel (should have nchannel or nchannel+1 entries)
 @param deviation_pairs The deviation pair anchors and offsets
 @param accuracy Not used for this calibration type (kept for API consistency)
 @returns The channel number (fractional) corresponding to the energy
 */
template<typename T>
T find_lowerchannel_channel( const T energy,
                             const std::vector<float> &channel_energies,
                             const std::vector<std::pair<double,T>> &deviation_pairs,
                             const double accuracy = 0.001 )
{
  if( channel_energies.size() < 2 )
    return T(0.0);

  const size_t nchannel = channel_energies.size() - 1;  // Last entry is upper edge of last channel

  // Build deviation pair spline
  const std::vector<CubicSplineNodeT<T>> dev_pair_spline =
    deviation_pairs.empty() ? std::vector<CubicSplineNodeT<T>>{} : create_cubic_spline( deviation_pairs );

  // Apply deviation pair correction to the energy
  T corrected_energy = energy;
  if( !dev_pair_spline.empty() )
    corrected_energy -= eval_cubic_spline( energy, dev_pair_spline );

  // Binary search to find the channel
  // Find largest i such that channel_energies[i] <= corrected_energy
  if( corrected_energy < static_cast<double>(channel_energies[0]) )
  {
    if constexpr ( !std::is_same_v<T, double> )
      corrected_energy.a = 0.0;
    else 
      corrected_energy = 0.0;
    return corrected_energy;
  }

  if( corrected_energy >= static_cast<double>(channel_energies[nchannel]) )
  {
    if constexpr ( !std::is_same_v<T, double> )
      corrected_energy.a = static_cast<double>(nchannel - 1);
    else 
      corrected_energy = static_cast<double>(nchannel - 1);
    return corrected_energy;
  }
  
  // Binary search
  size_t low = 0;
  size_t high = nchannel;

  while( high - low > 1 )
  {
    const size_t mid = (low + high) / 2;

    if( static_cast<double>(channel_energies[mid]) <= corrected_energy )
      low = mid;
    else
      high = mid;
  }

  // Interpolate within the channel
  const double e_low = channel_energies[low];
  const double e_high = channel_energies[low + 1];
  const T fraction = (corrected_energy - e_low) / (e_high - e_low);

  return T(static_cast<double>(low) + fraction);
}//find_lowerchannel_channel(...)


/** Find the channel corresponding to a given energy for lower channel edge calibration with offset and gain adjustments.

 This version applies offset_adj and gain_adj to the energy calibration before finding the channel.
 The energy adjustment formula is:
   adjusted_energy = orig_energy + offset_adj + (range_frac * gain_adj)
 where range_frac = (orig_energy - lower_energy) / range

 Instead of creating adjusted channel energies, we invert the adjustment and apply it to the input energy,
 preserving derivative information from offset_adj and gain_adj.

 @param energy The energy in keV to find the channel for
 @param channel_energies Lower channel edge energies (size = nchannel + 1)
 @param deviation_pairs The deviation pair anchors and offsets
 @param offset_adj Energy calibration offset adjustment in keV
 @param gain_adj Energy calibration gain adjustment in keV
 @param accuracy Not used for this calibration type (kept for API consistency)
 @returns The channel number (fractional) corresponding to the energy
 */
template<typename T>
T find_lowerchannel_channel( const T energy,
                             const std::vector<float> &channel_energies,
                             const std::vector<std::pair<double,T>> &deviation_pairs,
                             const T offset_adj,
                             const T gain_adj,
                             const double accuracy = 0.001 )
{
  using namespace std;
  using namespace ceres;

  if( channel_energies.size() < 2 )
    return T(0.0);

  const size_t nchannel = channel_energies.size() - 1;  // Last entry is upper edge of last channel

  // Get the energy range
  const double lower_energy = channel_energies[0];
  const double upper_energy = channel_energies[nchannel];
  const double range = upper_energy - lower_energy;

  if( range <= 0.0 )
    return T(0.0);

  // Build deviation pair spline
  const std::vector<CubicSplineNodeT<T>> dev_pair_spline =
    deviation_pairs.empty() ? std::vector<CubicSplineNodeT<T>>{} : create_cubic_spline( deviation_pairs );

  // Apply deviation pair correction to the energy first
  T corrected_energy = energy;
  if( !dev_pair_spline.empty() )
    corrected_energy -= eval_cubic_spline( energy, dev_pair_spline );

  // The adjustment transforms channel energies as:
  //   E_adj[i] = E_orig[i] + offset_adj + ((E_orig[i] - lower_energy) / range) * gain_adj
  // We need to invert this to find which original channel corresponds to our corrected_energy:
  //   corrected_energy = E_orig + offset_adj + ((E_orig - lower_energy) / range) * gain_adj
  //   corrected_energy = E_orig + offset_adj + (E_orig / range - lower_energy / range) * gain_adj
  //   corrected_energy = E_orig + offset_adj + E_orig * gain_adj / range - lower_energy * gain_adj / range
  //   corrected_energy = E_orig * (1 + gain_adj / range) + offset_adj - lower_energy * gain_adj / range
  //   E_orig = (corrected_energy - offset_adj + lower_energy * gain_adj / range) / (1 + gain_adj / range)

  const T unadjusted_energy = (corrected_energy - offset_adj + lower_energy * gain_adj / range) / (T(1.0) + gain_adj / range);

  // Now do binary search on the original channel_energies with the unadjusted energy
  double search_energy;
  if constexpr ( std::is_same_v<T, double> )
    search_energy = unadjusted_energy;
  else
    search_energy = unadjusted_energy.a;

  // Binary search to find the channel
  if( search_energy < static_cast<double>(channel_energies[0]) )
  {
    T answer = unadjusted_energy;
    if constexpr ( !std::is_same_v<T, double> )
      answer.a = 0.0;
    else
      answer = 0.0;
    return answer;
  }

  if( search_energy >= static_cast<double>(channel_energies[nchannel]) )
  {
    T answer = unadjusted_energy;
    if constexpr ( !std::is_same_v<T, double> )
      answer.a = static_cast<double>(nchannel - 1);
    else
      answer = static_cast<double>(nchannel - 1);
    return answer;
  }

  // Binary search
  size_t low = 0;
  size_t high = nchannel;

  while( (high - low) > 1 )
  {
    const size_t mid = (low + high) / 2;

    if( static_cast<double>(channel_energies[mid]) <= search_energy )
      low = mid;
    else
      high = mid;
  }

  // Interpolate within the channel
  const double e_low = channel_energies[low];
  const double e_high = channel_energies[low + 1];
  const T fraction = (unadjusted_energy - e_low) / (e_high - e_low);

  return T(static_cast<double>(low)) + fraction;
}//find_lowerchannel_channel(...) with offset_adj and gain_adj


}//namespace RelActCalcAutoImp

#endif //RelActCalcAuto_EnergyCal_imp_hpp
