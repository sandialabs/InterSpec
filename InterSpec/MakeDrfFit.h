#ifndef MakeDrfFit_h
#define MakeDrfFit_h
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

#include <deque>
#include <string>
#include <vector>
#include <memory>

#include "InterSpec/DetectorPeakResponse.h"

class PeakDef;

/** Fits of the detector FWHM and full-energy-peak efficiency functional forms.

 Both fits are Levenberg-Marquardt (Ceres, automatic differentiation through the templated
 evaluators in DetectorPeakResponse), seeded by closed-form linear least squares where the form
 allows it.  Results carry the full coefficient covariance.
 */
namespace MakeDrfFit
{
  /** Result of a FWHM functional-form fit. */
  struct FwhmFitResult
  {
    /** One entry per parameter of the functional form, including any held fixed (see
     #performResolutionFitEx for when parameters are held fixed). */
    std::vector<float> coefs;

    /** 1-sigma uncertainties; 0 for parameters held fixed. */
    std::vector<float> uncerts;

    /** Row-major coefs.size() x coefs.size() covariance; rows/columns of fixed parameters are
     zero.  Empty if it could not be computed. */
    std::vector<float> covRowMajor;

    /** Sum of squared normalized width residuals at the solution (not divided by dof). */
    double chi2 = 0.0;

    /** Number of (Gaussian) peaks minus number of free parameters; may be <= 0. */
    int dof = 0;

    /** Non-fatal issues (e.g. covariance unavailable); empty if none. */
    std::string warnings;
  };//struct FwhmFitResult


  /** Fits the FWHM functional form `fnctnlForm` to the widths of `peaks`.

   The residual for each Gaussian peak is (sigma_pred - sigma_meas)/sigma_uncert, with sigma_uncert
   the larger of the peaks fit sigma uncertainty and 1% of its sigma (5% of sigma when the peak has
   no uncertainty) - identical to #peak_width_chi2.

   With few peaks some parameters are held at their starting values, matching the long-standing
   behaviour callers rely on:
    - kGadrasResolutionFcn: 1 peak fits only B; 2 peaks fit B and C; 3+ fit all (A is searched
      separately over its negative and positive branches, since they are different functions).
    - kSqrtEnergyPlusInverse: 1 peak fits B (A=C=0); 2 peaks fit A,B (C=0); 3+ all.
    - kConstantPlusSqrtEnergy: 1 peak fits B (A=0); 2+ both.
    - kSqrtPolynomial: `sqrtEqnOrder` coefficients when the linear seed succeeds; otherwise 2
      coefficients (only the first free with < 3 peaks), or 3 with more than 3 peaks - i.e. it may
      return FEWER coefficients than `sqrtEqnOrder`.

   @param starting_coefs If sized for the form (3 for GADRAS) these seed the fit; else defaults.

   Throws std::runtime_error on invalid input or when the minimization fails.
   */
  FwhmFitResult performResolutionFitEx( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                                        const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                                        const int sqrtEqnOrder,
                                        const std::vector<float> &starting_coefs );

  /* If result passed into the function is of the proper size, then it will be
    used as the starting parameters for the fit, otherwise some default values
    will be chosen to start the fit with.  See #performResolutionFitEx.
   
    @param sqrtEqnOrder Only used if fnctnlForm==DetectorPeakResponse::kSqrtPolynomial
   
    @returns chi2 (not divided by dof).
   
    Throws exception on error with a kinda explanatory message.
  */
  double performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                             const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                             const int sqrtEqnOrder,
                             std::vector<float> &result,
                             std::vector<float> &uncerts );
  
  /** For use when fitting a FWHM equation from fit-peaks in a spectrum; use this function to filter peaks with anomolous widths from the set of peaks
   you are using to fit the FWHM.
   Generally you might call `performResolutionFit(...)`  with all of your peaks, then call this function to filter the outliers out, then call
   `performResolutionFit(...)` again.
   
   Does no work, and just retuns the input peaks if there are less than 5 input peaks.
   */
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>>
  removeOutlyingWidthPeaks( const std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &peaks,
                            const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                           const std::vector<float> &coefficients );
  
  /** Fits the FWHM as a sqrt( Sum_i {A_i *pow(x,i)} )
   
   @param peaks The peaks to use for the fit
   @param num_fit_coefficients The number of coefficients to fit _not_ the fit order (which would be num energy
          dependent coefficients).  Must be larger than zero
   @param include_inv_term If true, fits sqrt(A + B*x + C/x) instead of sqrt(A + B*x + C*x*x + D*x*x + ... )
          Right now, if true, assumes will fit for three coefficients.
   @param[out] coeffs The resulting fit coefficients.
   @param[out] coeffs The resulting fit coefficient uncertainties.
   @returns the chi2 of the fit
   
   Throws exception if fitting for more parameters than there are peaks, or if number of coefficients isnt at least 1, or the fit fails.
   */
  double fit_sqrt_poly_fwhm_lls( const std::deque< std::shared_ptr<const PeakDef> > &peaks,
                                       const int num_fit_coefficients,
                                       const bool include_inv_term,
                                       std::vector<float> &coeffs,
                                       std::vector<float> &coeff_uncerts );
  
  /** Fits `FWHM = A_0 + A_1*sqrt(x)
   
   Throws exception if fit fails.
   */
  double fit_constant_plus_sqrt_fwhm_lls( const std::deque< std::shared_ptr<const PeakDef> > &peaks,
                                  std::vector<float> &coeffs,
                                  std::vector<float> &coeff_uncerts );
  
  double peak_width_chi2( double predicted_sigma, const PeakDef &peak );
  
  
  /** One measured efficiency point for #performEfficiencyFit.

   The efficiency is whatever the fitted curve is meant to describe (far-field intrinsic for a
   normal DRF); the caller has already divided out the source geometry.  The three fractional
   uncertainties are combined into the data covariance by #effDataCovariance: the statistical part
   is independent between points, while the certificate and distance parts are 100% correlated
   among points sharing a #sourceKey (an error in a sources activity or position moves all of its
   peaks together).
   */
  struct EffFitPoint
  {
    /** Energy in the units the equation will be fit in (keV or MeV; detected by whether the
     largest energy exceeds 30). */
    float energy = 0.0f;

    float efficiency = 0.0f;

    /** Fractional 1-sigma uncertainties: independent (peak area), and correlated within
     `sourceKey` (source activity certificate, and source distance). */
    float fracStatUncert = 0.0f;
    float fracCertUncert = 0.0f;
    float fracDistUncert = 0.0f;

    std::string sourceKey;
  };//struct EffFitPoint


  /** Result of an efficiency functional-form fit. */
  struct EffFitResult
  {
    /** Coefficients of eff(x) = exp(c0 + c1*log(x) + c2*log(x)^2 + ...), in the energy units of
     the input points. */
    std::vector<float> coefs;

    /** 1-sigma coefficient uncertainties (sqrt of the covariance diagonal). */
    std::vector<float> uncerts;

    /** Row-major coefs.size() x coefs.size() coefficient covariance, already scaled by
     #birgeScale.  Empty if it could not be computed. */
    std::vector<float> covRowMajor;

    /** Sum of squared whitened residuals at the solution (not divided by dof). */
    double chi2 = 0.0;

    /** Number of points minus number of coefficients; may be 0. */
    int dof = 0;

    /** max(1, chi2/dof): when the scatter of the points exceeds their stated uncertainties, the
     covariance is inflated by this factor so the curve uncertainty reflects the actual scatter. */
    double birgeScale = 1.0;

    /** Non-fatal issues (e.g. covariance unavailable, non-convergence); empty if none. */
    std::string warnings;
  };//struct EffFitResult


  /** Row-major fractional data covariance of the points:
      C[i][j] = delta_ij*stat_i^2 + [same non-empty sourceKey]*(cert_i*cert_j + dist_i*dist_j)
   */
  std::vector<double> effDataCovariance( const std::vector<EffFitPoint> &data );


  /** Linear least squares fit of log(eff) to a polynomial in log(energy), each point weighted by
   its total fractional uncertainty (correlations ignored).  Used to seed #performEfficiencyFit.

   @param cov_row_major If non-null, receives the order x order coefficient covariance.
   @returns chi2 in linear (efficiency) space.

   Throws std::runtime_error on invalid input.
   */
  double fit_intrinsic_eff_least_linear_squares( const std::vector<EffFitPoint> &data,
                                                 const int order,
                                                 std::vector<float> &coeffs,
                                                 std::vector<float> &coeff_uncerts,
                                                 std::vector<double> *cov_row_major );


  /** Fits the data efficiencies to eff(x) = exp(A + B*log(x) + C*log(x)^2 + ...) with `fcnOrder`
   coefficients, using the full data covariance of #effDataCovariance (residuals are whitened by
   its Cholesky factor, so correlated source errors are handled correctly rather than being
   treated as independent scatter).

   The returned coefficients are in the same energy units as EffFitPoint::energy.

   Throws std::runtime_error on invalid input or when the minimization fails.
   */
  EffFitResult performEfficiencyFit( const std::vector<EffFitPoint> &data, const int fcnOrder );


  /** Legacy efficiency point: `efficiency_uncert` is absolute, treated as independent. */
  struct DetEffDataPoint{ float energy, efficiency, efficiency_uncert; };

  /** Legacy interface to #performEfficiencyFit: points are treated as uncorrelated (5% when no
   uncertainty is given), and only the diagonal uncertainties are returned.

   Detects if DetEffDataPoint::energy is in keV or MeV by testing if largest energy is greater
   than 30; if so, in keV, else MeV.  Returned coefficients will be in the same energy units as
   DetEffDataPoint::energy.

   @returns chi2/dof (or chi2 when there are exactly as many points as coefficients).

   Throws exception on error with a kinda explanatory message.
   */
  double performEfficiencyFit( const std::vector<DetEffDataPoint> data,
                             const int fcnOrder,
                             std::vector<float> &result,
                             std::vector<float> &uncerts );
  
}//namespace MakeDrfFit

#endif  //MakeDrfFit_h
