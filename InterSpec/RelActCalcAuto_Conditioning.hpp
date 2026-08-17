/* InterSpec: an application to analyze spectral gamma radiation data.
 *
 * Internal conditioning utilities for empirical RelActAuto efficiency curves.
 */
#ifndef RelActCalcAuto_Conditioning_hpp
#define RelActCalcAuto_Conditioning_hpp

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Eigen/Dense"

#include "InterSpec/RelActCalc.h"

namespace RelActCalcAutoImp
{

/** Value of one coefficient's legacy empirical basis function.

 The definitions intentionally mirror RelActCalc::eval_eqn_imp().  Keeping this
 helper separate lets the optimizer use a better-conditioned coordinate system
 while reports and persisted coefficients retain the legacy basis exactly.
 */
inline double empirical_legacy_basis_value( const RelActCalc::RelEffEqnForm form,
                                            const size_t coefficient,
                                            const double energy )
{
  if( !(energy > 0.0) || !std::isfinite(energy) )
    throw std::runtime_error( "Empirical efficiency basis requires a finite positive energy." );

  const double log_energy = std::log( energy );
  switch( form )
  {
    case RelActCalc::RelEffEqnForm::LnX:
    case RelActCalc::RelEffEqnForm::LnXLnY:
      if( coefficient == 0 ) return 1.0;
      if( coefficient == 1 ) return log_energy;
      return std::pow( log_energy, static_cast<double>(coefficient) );

    case RelActCalc::RelEffEqnForm::LnY:
      if( coefficient == 0 ) return 1.0;
      if( coefficient == 1 ) return energy;
      if( coefficient == 2 ) return 1.0 / energy;
      return 1.0 / std::pow( energy, static_cast<double>(coefficient - 1) );

    case RelActCalc::RelEffEqnForm::FramEmpirical:
      if( coefficient == 0 ) return 1.0;
      if( coefficient == 1 ) return 1.0 / (energy*energy);
      return std::pow( log_energy, static_cast<double>(coefficient - 1) );

    case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      throw std::runtime_error( "The physical relative-efficiency model has no empirical legacy basis." );
  }

  throw std::runtime_error( "Unknown empirical relative-efficiency form." );
}


inline bool empirical_form_is_exponential( const RelActCalc::RelEffEqnForm form )
{
  switch( form )
  {
    case RelActCalc::RelEffEqnForm::LnX: return false;
    case RelActCalc::RelEffEqnForm::LnY:
    case RelActCalc::RelEffEqnForm::LnXLnY:
    case RelActCalc::RelEffEqnForm::FramEmpirical: return true;
    case RelActCalc::RelEffEqnForm::FramPhysicalModel: break;
  }
  throw std::runtime_error( "The physical relative-efficiency model is not empirical." );
}


/** Fixed-pivot, QR-orthogonalized coordinates for one empirical efficiency curve.

 Slot zero is the (fixed) normalization gauge.  Slots 1..N-1 are coefficients
 of QR-orthonormal columns constructed from the centered legacy basis over the
 frozen ROI channel energies.  If B is that centered legacy design matrix and
 B=QR, optimizer coordinates are q=R*c_shape.  Consequently d(model)/dq is Q
 at the construction energies, while c_shape=R^-1*q maps exactly back to the
 legacy coefficient convention.

 The pivot constraint is exactly f(E_pivot)=1.  For LnX this is an affine
 constraint on the polynomial; for the exponential forms it is a zero log
 efficiency constraint.  The removed scale is carried by every activity of the
 curve, so predicted counts are unchanged.
 */
class EmpiricalBasisTransform
{
public:
  struct Seed
  {
    std::vector<double> orthogonal_coefficients;
    std::vector<double> normalized_legacy_coefficients;
    double activity_scale = 1.0;
  };

  EmpiricalBasisTransform() = default;

  EmpiricalBasisTransform( const RelActCalc::RelEffEqnForm form,
                           const size_t num_coefficients,
                           const double pivot_energy,
                           const std::vector<double> &frozen_energies )
    : m_form( form ),
      m_num_coefficients( num_coefficients ),
      m_pivot_energy( pivot_energy ),
      m_pivot_basis( num_coefficients, 0.0 ),
      m_shape_from_orthogonal(),
      m_orthogonal_from_shape()
  {
    if( form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      throw std::runtime_error( "Cannot construct an empirical basis transform for a physical model." );
    if( num_coefficients < 1 )
      throw std::runtime_error( "An empirical efficiency curve needs at least one coefficient." );
    if( !(pivot_energy > 0.0) || !std::isfinite(pivot_energy) )
      throw std::runtime_error( "Empirical efficiency pivot must be finite and positive." );

    for( size_t j = 0; j < num_coefficients; ++j )
      m_pivot_basis[j] = empirical_legacy_basis_value( form, j, pivot_energy );

    const size_t num_shape = num_coefficients - 1;
    if( num_shape == 0 )
      return;

    std::vector<double> energies;
    energies.reserve( frozen_energies.size() );
    for( const double energy : frozen_energies )
    {
      if( (energy > 0.0) && std::isfinite(energy) )
        energies.push_back( energy );
    }
    if( energies.size() < num_shape )
      throw std::runtime_error( "Not enough frozen ROI energies to orthogonalize the empirical basis." );

    Eigen::MatrixXd design( static_cast<Eigen::Index>(energies.size()),
                            static_cast<Eigen::Index>(num_shape) );
    for( size_t row = 0; row < energies.size(); ++row )
    {
      for( size_t col = 0; col < num_shape; ++col )
      {
        design( static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col) )
          = empirical_legacy_basis_value(form, col + 1, energies[row]) - m_pivot_basis[col + 1];
      }
    }

    // No column pivoting: the upper-triangular relationship preserves degree
    // nesting (q_last == 0 iff the highest legacy coefficient is zero), which
    // keeps model-order removal semantically exact.
    const Eigen::HouseholderQR<Eigen::MatrixXd> qr( design );
    Eigen::MatrixXd r = qr.matrixQR().topLeftCorner(
                          static_cast<Eigen::Index>(num_shape),
                          static_cast<Eigen::Index>(num_shape) )
                        .template triangularView<Eigen::Upper>();

    const double scale = (std::max)( 1.0, r.cwiseAbs().maxCoeff() );
    for( size_t diagonal = 0; diagonal < num_shape; ++diagonal )
    {
      const double value = std::fabs( r(static_cast<Eigen::Index>(diagonal),
                                        static_cast<Eigen::Index>(diagonal)) );
      if( !(value > 64.0*std::numeric_limits<double>::epsilon()*scale) )
        throw std::runtime_error( "Frozen ROI energies do not identify the requested empirical basis order." );
    }

    const Eigen::MatrixXd r_inverse
      = r.template triangularView<Eigen::Upper>().solve(
          Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(num_shape),
                                    static_cast<Eigen::Index>(num_shape)) );

    m_orthogonal_from_shape.assign( num_shape, std::vector<double>(num_shape, 0.0) );
    m_shape_from_orthogonal.assign( num_shape, std::vector<double>(num_shape, 0.0) );
    for( size_t row = 0; row < num_shape; ++row )
    {
      for( size_t col = 0; col < num_shape; ++col )
      {
        m_orthogonal_from_shape[row][col]
          = r(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
        m_shape_from_orthogonal[row][col]
          = r_inverse(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
      }
    }
  }

  RelActCalc::RelEffEqnForm form() const { return m_form; }
  size_t num_coefficients() const { return m_num_coefficients; }
  double pivot_energy() const { return m_pivot_energy; }

  /** Convert a legacy seed and return the common scale to absorb into activities. */
  Seed orthogonalize_legacy_seed( const std::vector<double> &legacy ) const
  {
    if( legacy.size() != m_num_coefficients )
      throw std::runtime_error( "Wrong number of empirical legacy seed coefficients." );

    Seed answer;
    answer.normalized_legacy_coefficients = legacy;
    answer.orthogonal_coefficients.assign( m_num_coefficients, 0.0 );

    double pivot_linear_value = 0.0;
    for( size_t j = 0; j < m_num_coefficients; ++j )
      pivot_linear_value += legacy[j] * m_pivot_basis[j];

    if( empirical_form_is_exponential(m_form) )
    {
      // Stay in the ordinary exponential regime.  A seed outside it is not a
      // meaningful physical initialization and the caller should use identity.
      if( !std::isfinite(pivot_linear_value) || (std::fabs(pivot_linear_value) >= 50.0) )
        throw std::runtime_error( "Empirical seed normalization is outside the finite exponential regime." );
      answer.activity_scale = std::exp( pivot_linear_value );
      answer.normalized_legacy_coefficients[0] -= pivot_linear_value;
    }else
    {
      if( !std::isfinite(pivot_linear_value) || !(pivot_linear_value > 0.0) )
        throw std::runtime_error( "LnX seed is not positive at the fixed normalization pivot." );
      answer.activity_scale = pivot_linear_value;
      for( double &coefficient : answer.normalized_legacy_coefficients )
        coefficient /= pivot_linear_value;
    }

    const size_t num_shape = m_num_coefficients - 1;
    for( size_t row = 0; row < num_shape; ++row )
    {
      double q = 0.0;
      for( size_t col = row; col < num_shape; ++col )
        q += m_orthogonal_from_shape[row][col]
             * answer.normalized_legacy_coefficients[col + 1];
      answer.orthogonal_coefficients[row + 1] = q;
    }
    return answer;
  }

  /** Convert optimizer coordinates to the exact legacy coefficient basis. */
  template<typename T>
  std::vector<T> legacy_coefficients( const T * const orthogonal,
                                      const size_t count ) const
  {
    if( !orthogonal || (count != m_num_coefficients) )
      throw std::runtime_error( "Wrong number of empirical orthogonal coefficients." );

    std::vector<T> legacy( m_num_coefficients, T(0.0) );
    const size_t num_shape = m_num_coefficients - 1;
    for( size_t row = 0; row < num_shape; ++row )
    {
      for( size_t col = 0; col < num_shape; ++col )
        legacy[row + 1] += T(m_shape_from_orthogonal[row][col]) * orthogonal[col + 1];
    }

    legacy[0] = empirical_form_is_exponential(m_form) ? T(0.0) : T(1.0);
    for( size_t row = 0; row < num_shape; ++row )
      legacy[0] -= T(m_pivot_basis[row + 1]) * legacy[row + 1];
    return legacy;
  }

  std::vector<double> legacy_coefficients( const std::vector<double> &orthogonal ) const
  {
    return legacy_coefficients( orthogonal.data(), orthogonal.size() );
  }

  /** Evaluate the anchored legacy linear predictor directly from optimizer coordinates.

   For LnX this is the efficiency itself.  For the exponential forms the caller applies the
   project's standard exponential/overflow continuation.  This allocation-free route is used in
   the hot residual loop; it is algebraically identical to legacy_coefficients() followed by the
   legacy basis evaluation.
   */
  template<typename T>
  T linear_predictor( const double energy, const T * const orthogonal,
                      const size_t count ) const
  {
    if( !orthogonal || (count != m_num_coefficients) )
      throw std::runtime_error( "Wrong number of empirical orthogonal coefficients." );

    const size_t num_shape = m_num_coefficients - 1;
    T answer = empirical_form_is_exponential(m_form) ? T(0.0) : T(1.0);
    for( size_t col = 0; col < num_shape; ++col )
    {
      double conditioned_basis = 0.0;
      for( size_t row = 0; row < num_shape; ++row )
      {
        const double centered_legacy
          = empirical_legacy_basis_value(m_form,row + 1,energy) - m_pivot_basis[row + 1];
        conditioned_basis += centered_legacy*m_shape_from_orthogonal[row][col];
      }
      answer += T(conditioned_basis)*orthogonal[col + 1];
    }
    return answer;
  }

  /** Jacobian d(legacy coefficients)/d(optimizer coefficients). */
  std::vector<std::vector<double>> legacy_jacobian() const
  {
    std::vector<std::vector<double>> jacobian(
      m_num_coefficients, std::vector<double>(m_num_coefficients, 0.0) );
    const size_t num_shape = m_num_coefficients - 1;
    for( size_t row = 0; row < num_shape; ++row )
    {
      for( size_t col = 0; col < num_shape; ++col )
      {
        jacobian[row + 1][col + 1] = m_shape_from_orthogonal[row][col];
        jacobian[0][col + 1] -= m_pivot_basis[row + 1]
                                  * m_shape_from_orthogonal[row][col];
      }
    }
    return jacobian;
  }

  /** Map an optimizer-coordinate covariance to legacy coefficient covariance. */
  std::vector<std::vector<double>> legacy_covariance(
                                      const std::vector<std::vector<double>> &orthogonal_covariance ) const
  {
    if( orthogonal_covariance.size() != m_num_coefficients )
      throw std::runtime_error( "Wrong empirical covariance size." );
    for( const std::vector<double> &row : orthogonal_covariance )
    {
      if( row.size() != m_num_coefficients )
        throw std::runtime_error( "Wrong empirical covariance row size." );
    }

    const std::vector<std::vector<double>> jacobian = legacy_jacobian();
    std::vector<std::vector<double>> answer(
      m_num_coefficients, std::vector<double>(m_num_coefficients, 0.0) );
    for( size_t row = 0; row < m_num_coefficients; ++row )
    {
      for( size_t col = 0; col < m_num_coefficients; ++col )
      {
        for( size_t i = 0; i < m_num_coefficients; ++i )
          for( size_t j = 0; j < m_num_coefficients; ++j )
            answer[row][col] += jacobian[row][i] * orthogonal_covariance[i][j]
                                  * jacobian[col][j];
      }
    }
    return answer;
  }

private:
  RelActCalc::RelEffEqnForm m_form = RelActCalc::RelEffEqnForm::LnX;
  size_t m_num_coefficients = 0;
  double m_pivot_energy = 0.0;
  std::vector<double> m_pivot_basis;
  std::vector<std::vector<double>> m_shape_from_orthogonal;
  std::vector<std::vector<double>> m_orthogonal_from_shape;
};

}//namespace RelActCalcAutoImp

#endif //RelActCalcAuto_Conditioning_hpp
