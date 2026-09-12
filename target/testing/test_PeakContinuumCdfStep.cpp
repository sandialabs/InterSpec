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
// Must be defined before Windows.h (or any header that includes it) is included; boost/test's
//  included/unit_test.hpp pulls <windows.h>, while PeakFit_imp.hpp pulls boost/asio.
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#include <cmath>
#include <deque>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>

#define BOOST_TEST_MODULE PeakContinuumCdfStep_suite
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include <Wt/WFlags.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/PeakFit_imp.hpp"

using namespace std;
using namespace boost::unit_test;

// Regression tests for the peak-CDF step continua (FlatStepCDF / LinearStepCDF / BiLinearStepCDF).
//
// The headline bug these guard against: `PeakContinuum::setType()` used to carry a *non-CDF* step
// magnitude into the *CDF* step-coefficient slot.  The two are numerically incompatible - the
// non-CDF step is a count density (counts/keV, O(1..100)) while the CDF step coefficient is in
// keV^-1 (O(1e-3..1e-6)), differing by the total ROI peak area, i.e. 1e3 to 1e6.  Carrying one
// into the other drove the amplitude least-squares solve to a ~zero peak amplitude (the peak
// vanished from the chart) and inflated the continuum to ~3x the true level; the chi2 surface out
// there is flat enough that Ceres could not climb back.
//
// `setType()` is now driven by a table of what each parameter slot *means*, so a value is carried
// across a type change only into a slot with the identical meaning.  The test below re-states that
// table independently - it is the specification, deliberately duplicated rather than shared with
// the implementation - and checks every ordered pair of continuum types against it.

namespace
{
/** Mirror of the `ContParMeaning` table in src/PeakDef.cpp; deliberately duplicated so this test
 states the expected transition behaviour independently of the implementation.
 */
enum class Meaning
{
  Poly0, Poly1, Poly2, Poly3,
  RightPoly0, RightPoly1,
  DataStep,   //!< FlatStep/LinearStep step magnitude; counts/keV
  CdfStep0,   //!< CDF step types' step coefficient; keV^-1
  CdfStep1    //!< BiLinearStepCDF's step-coefficient energy slope; keV^-2
};

const char *to_str( const Meaning m )
{
  switch( m )
  {
    case Meaning::Poly0:      return "Poly0";
    case Meaning::Poly1:      return "Poly1";
    case Meaning::Poly2:      return "Poly2";
    case Meaning::Poly3:      return "Poly3";
    case Meaning::RightPoly0: return "RightPoly0";
    case Meaning::RightPoly1: return "RightPoly1";
    case Meaning::DataStep:   return "DataStep";
    case Meaning::CdfStep0:   return "CdfStep0";
    case Meaning::CdfStep1:   return "CdfStep1";
  }
  return "?";
}


vector<Meaning> meanings_for( const PeakContinuum::OffsetType type )
{
  using M = Meaning;
  switch( type )
  {
    case PeakContinuum::NoOffset:        return {};
    case PeakContinuum::Constant:        return { M::Poly0 };
    case PeakContinuum::Linear:          return { M::Poly0, M::Poly1 };
    case PeakContinuum::Quadratic:       return { M::Poly0, M::Poly1, M::Poly2 };
    case PeakContinuum::Cubic:           return { M::Poly0, M::Poly1, M::Poly2, M::Poly3 };
    case PeakContinuum::FlatStep:        return { M::Poly0, M::DataStep };
    case PeakContinuum::LinearStep:      return { M::Poly0, M::Poly1, M::DataStep };
    case PeakContinuum::BiLinearStep:    return { M::Poly0, M::Poly1, M::RightPoly0, M::RightPoly1 };
    case PeakContinuum::FlatStepCDF:     return { M::Poly0, M::CdfStep0 };
    case PeakContinuum::LinearStepCDF:   return { M::Poly0, M::Poly1, M::CdfStep0 };
    // Reparameterised: a linear continuum plus a linearly varying step coefficient.  Algebraically
    //  the same family as blending two lines, but linear in the peak amplitudes.
    case PeakContinuum::BiLinearStepCDF: return { M::Poly0, M::Poly1, M::CdfStep0, M::CdfStep1 };
    case PeakContinuum::External:        return {};
  }
  return {};
}//meanings_for(...)


const vector<PeakContinuum::OffsetType> &all_offset_types()
{
  static const vector<PeakContinuum::OffsetType> s_types{
    PeakContinuum::NoOffset,     PeakContinuum::Constant,
    PeakContinuum::Linear,       PeakContinuum::Quadratic,
    PeakContinuum::Cubic,        PeakContinuum::FlatStep,
    PeakContinuum::LinearStep,   PeakContinuum::BiLinearStep,
    PeakContinuum::FlatStepCDF,  PeakContinuum::LinearStepCDF,
    PeakContinuum::BiLinearStepCDF, PeakContinuum::External
  };
  return s_types;
}


size_t index_of( const vector<Meaning> &meanings, const Meaning m )
{
  for( size_t i = 0; i < meanings.size(); ++i )
  {
    if( meanings[i] == m )
      return i;
  }
  return meanings.size();
}
}//namespace


// The meaning table must stay in lock-step with `PeakContinuum::num_parameters()`; if it does not,
//  every other assertion here is meaningless.
BOOST_AUTO_TEST_CASE( meaning_table_matches_num_parameters )
{
  for( const PeakContinuum::OffsetType type : all_offset_types() )
  {
    BOOST_CHECK_MESSAGE( meanings_for(type).size() == PeakContinuum::num_parameters(type),
      "Meaning table size " << meanings_for(type).size() << " != num_parameters() "
      << PeakContinuum::num_parameters(type) << " for type "
      << PeakContinuum::offset_type_str(type) );
  }
}


// The direct regression test for the headline bug, and for the three further transitions that
//  leaked a step magnitude into a polynomial slot (FlatStep*->Quadratic/Cubic,
//  LinearStep*->Cubic, FlatStepCDF->BiLinearStep*).
BOOST_AUTO_TEST_CASE( setType_preserves_only_same_meaning_parameters )
{
  for( const PeakContinuum::OffsetType from : all_offset_types() )
  {
    for( const PeakContinuum::OffsetType to : all_offset_types() )
    {
      const vector<Meaning> from_meanings = meanings_for( from );
      const vector<Meaning> to_meanings = meanings_for( to );

      PeakContinuum cont;
      cont.setType( from );
      cont.setRange( 100.0, 120.0 );

      // Distinct, easily recognisable values so a mis-copied slot is unambiguous.
      vector<double> values( from_meanings.size() ), uncerts( from_meanings.size() );
      for( size_t i = 0; i < from_meanings.size(); ++i )
      {
        values[i] = 100.0 + static_cast<double>(i);
        uncerts[i] = 1.0 + 0.125*static_cast<double>(i);
      }
      cont.setParameters( 110.0, values, uncerts );

      // Pin every other slot, so we can check a pin never rides into a slot that means something
      //  else - carrying a "do not fit" flag onto a different quantity is the subtle half of the
      //  bug this table exists to prevent.
      vector<bool> pinned( from_meanings.size(), false );
      for( size_t i = 1; i < from_meanings.size(); i += 2 )
      {
        pinned[i] = true;
        BOOST_REQUIRE( cont.setPolynomialCoefFitFor( i, false ) );
      }

      cont.setType( to );

      const string ctx = string("from ") + PeakContinuum::offset_type_str(from)
                         + " to " + PeakContinuum::offset_type_str(to);

      BOOST_REQUIRE_MESSAGE( cont.parameters().size() == PeakContinuum::num_parameters(to),
                             ctx << ": wrong parameter count " << cont.parameters().size() );
      BOOST_REQUIRE_MESSAGE( cont.uncertainties().size() == PeakContinuum::num_parameters(to),
                             ctx << ": wrong uncertainty count " << cont.uncertainties().size() );
      BOOST_REQUIRE_MESSAGE( cont.fitForParameter().size() == PeakContinuum::num_parameters(to),
                             ctx << ": wrong fitFor count " << cont.fitForParameter().size() );

      // Does the new type introduce a right-hand line the old type did not have?  If so the
      //  implementation deliberately seeds it from the left line rather than from zero.
      const bool seeds_right_line
        = (index_of( to_meanings, Meaning::RightPoly0 ) < to_meanings.size())
          && (index_of( from_meanings, Meaning::RightPoly0 ) >= from_meanings.size());

      for( size_t i = 0; i < to_meanings.size(); ++i )
      {
        const Meaning meaning = to_meanings[i];
        const size_t src = index_of( from_meanings, meaning );

        double expected = (src < from_meanings.size()) ? values[src] : 0.0;

        if( seeds_right_line && (meaning == Meaning::RightPoly0 || meaning == Meaning::RightPoly1) )
        {
          const Meaning left = (meaning == Meaning::RightPoly0) ? Meaning::Poly0 : Meaning::Poly1;
          const size_t left_src = index_of( from_meanings, left );
          expected = (left_src < from_meanings.size()) ? values[left_src] : 0.0;
        }

        BOOST_CHECK_MESSAGE( fabs(cont.parameters()[i] - expected) < 1.0E-9,
          ctx << ": slot " << i << " (" << to_str(meaning) << ") is "
              << cont.parameters()[i] << ", expected " << expected );

        // Uncertainties travel with their value, and a slot with no source starts at zero.
        double expected_uncert = (src < from_meanings.size()) ? uncerts[src] : 0.0;
        bool expected_pinned = (src < from_meanings.size()) ? pinned[src] : false;

        if( seeds_right_line && (meaning == Meaning::RightPoly0 || meaning == Meaning::RightPoly1) )
        {
          const Meaning left = (meaning == Meaning::RightPoly0) ? Meaning::Poly0 : Meaning::Poly1;
          const size_t left_src = index_of( from_meanings, left );
          expected_uncert = (left_src < from_meanings.size()) ? uncerts[left_src] : 0.0;
          expected_pinned = (left_src < from_meanings.size()) ? pinned[left_src] : false;
        }

        BOOST_CHECK_MESSAGE( fabs(cont.uncertainties()[i] - expected_uncert) < 1.0E-9,
          ctx << ": slot " << i << " (" << to_str(meaning) << ") uncertainty is "
              << cont.uncertainties()[i] << ", expected " << expected_uncert );

        BOOST_CHECK_MESSAGE( cont.fitForParameter()[i] == !expected_pinned,
          ctx << ": slot " << i << " (" << to_str(meaning) << ") fitFor is "
              << cont.fitForParameter()[i] << ", expected " << !expected_pinned
              << " - a pin must follow its meaning, never a slot index" );
      }//for( size_t i = 0; i < to_meanings.size(); ++i )
    }//for( to )
  }//for( from )
}


// The specific pairs the user hit.  A CDF step coefficient is ~1/(ROI peak area) times the
//  non-CDF one, so a step magnitude of order the continuum density must never survive the switch.
BOOST_AUTO_TEST_CASE( setType_zeroes_step_across_cdf_boundary )
{
  const vector<pair<PeakContinuum::OffsetType,PeakContinuum::OffsetType>> cross_family{
    { PeakContinuum::FlatStep,      PeakContinuum::FlatStepCDF },
    { PeakContinuum::FlatStep,      PeakContinuum::LinearStepCDF },
    { PeakContinuum::LinearStep,    PeakContinuum::FlatStepCDF },
    { PeakContinuum::LinearStep,    PeakContinuum::LinearStepCDF },
    { PeakContinuum::FlatStepCDF,   PeakContinuum::FlatStep },
    { PeakContinuum::FlatStepCDF,   PeakContinuum::LinearStep },
    { PeakContinuum::LinearStepCDF, PeakContinuum::FlatStep },
    { PeakContinuum::LinearStepCDF, PeakContinuum::LinearStep },
  };

  for( const auto &pr : cross_family )
  {
    const vector<Meaning> from_meanings = meanings_for( pr.first );
    const vector<Meaning> to_meanings = meanings_for( pr.second );

    PeakContinuum cont;
    cont.setType( pr.first );
    cont.setRange( 100.0, 120.0 );

    // A step magnitude typical of a real non-CDF fit; catastrophic if used as a CDF coefficient.
    vector<double> values( from_meanings.size(), 0.0 );
    for( size_t i = 0; i < from_meanings.size(); ++i )
      values[i] = (from_meanings[i] == Meaning::DataStep) ? 30.0 : 300.0;
    cont.setParameters( 110.0, values, {} );

    cont.setType( pr.second );

    const string ctx = string(PeakContinuum::offset_type_str(pr.first)) + " -> "
                       + PeakContinuum::offset_type_str(pr.second);

    const size_t step_idx = (index_of( to_meanings, Meaning::DataStep ) < to_meanings.size())
                            ? index_of( to_meanings, Meaning::DataStep )
                            : index_of( to_meanings, Meaning::CdfStep0 );
    BOOST_REQUIRE_MESSAGE( step_idx < to_meanings.size(), ctx << ": target has no step slot" );

    BOOST_CHECK_MESSAGE( cont.parameters()[step_idx] == 0.0,
      ctx << ": step slot survived the CDF boundary as " << cont.parameters()[step_idx]
          << " - it means a different physical quantity on each side" );
  }//for( const auto &pr : cross_family )
}


// The same-family switches must NOT zero the step - that is the whole point of keeping them
//  distinguishable from the cross-family ones.
BOOST_AUTO_TEST_CASE( setType_preserves_step_within_family )
{
  const vector<pair<PeakContinuum::OffsetType,PeakContinuum::OffsetType>> same_family{
    { PeakContinuum::FlatStep,      PeakContinuum::LinearStep },
    { PeakContinuum::LinearStep,    PeakContinuum::FlatStep },
    { PeakContinuum::FlatStepCDF,   PeakContinuum::LinearStepCDF },
    { PeakContinuum::LinearStepCDF, PeakContinuum::FlatStepCDF },
    { PeakContinuum::FlatStepCDF,   PeakContinuum::BiLinearStepCDF },
    { PeakContinuum::BiLinearStepCDF, PeakContinuum::FlatStepCDF },
  };

  for( const auto &pr : same_family )
  {
    const vector<Meaning> from_meanings = meanings_for( pr.first );
    const vector<Meaning> to_meanings = meanings_for( pr.second );
    const bool is_cdf = PeakContinuum::is_peak_cdf_step_continuum( pr.first );
    const Meaning step_meaning = is_cdf ? Meaning::CdfStep0 : Meaning::DataStep;
    const double step_val = is_cdf ? -1.5E-3 : -30.0;

    PeakContinuum cont;
    cont.setType( pr.first );
    cont.setRange( 100.0, 120.0 );

    const size_t src_step_idx = index_of( from_meanings, step_meaning );
    BOOST_REQUIRE_MESSAGE( src_step_idx < from_meanings.size(),
      string(PeakContinuum::offset_type_str(pr.first)) + " has no " + to_str(step_meaning) + " slot" );

    vector<double> values( from_meanings.size(), 300.0 );
    values[src_step_idx] = step_val;
    cont.setParameters( 110.0, values, {} );

    cont.setType( pr.second );

    const string ctx = string(PeakContinuum::offset_type_str(pr.first)) + " -> "
                       + PeakContinuum::offset_type_str(pr.second);
    const size_t step_idx = index_of( to_meanings, step_meaning );
    BOOST_REQUIRE_MESSAGE( step_idx < to_meanings.size(), ctx << ": target has no matching step slot" );

    BOOST_CHECK_MESSAGE( fabs(cont.parameters()[step_idx] - step_val) < 1.0E-12,
      ctx << ": step should have been preserved, but is " << cont.parameters()[step_idx]
          << " instead of " << step_val );
  }//for( const auto &pr : same_family )
}


namespace
{
/** A synthetic, noise-free HPGe-like ROI: one Gaussian on a constant continuum with a step. */
struct SyntheticRoi
{
  shared_ptr<SpecUtils::Measurement> data;
  vector<float> energies;   //!< nchannel+1 lower channel energies
  vector<float> counts;
  size_t nchannel = 0;
  double lower_energy = 0.0, upper_energy = 0.0, ref_energy = 0.0;

  double mean = 600.0, sigma = 1.0, amplitude = 20000.0;
  double continuum_density = 300.0;    //!< counts/keV at the ROI's left edge
  double step_density = -30.0;         //!< total counts/keV change across the peak
};


double gauss_cdf( const double x, const double mean, const double sigma )
{
  return 0.5*( 1.0 + erf( (x - mean)/(sigma*sqrt(2.0)) ) );
}


SyntheticRoi make_synthetic_roi( const double step_density = -30.0, const double amplitude = 20000.0 )
{
  SyntheticRoi roi;
  roi.step_density = step_density;
  roi.amplitude = amplitude;

  // A full, realistically sized spectrum - a ROI-only spectrum trips code that legitimately looks
  //  at channels outside the ROI.
  const double chan_width = 0.35;
  const double spec_start = 400.0;
  const size_t spec_nchannel = 1024;

  vector<float> spec_energies( spec_nchannel + 1 );
  for( size_t i = 0; i <= spec_nchannel; ++i )
    spec_energies[i] = static_cast<float>( spec_start + i*chan_width );

  // Truth model: constant continuum, a step proportional to the peak's ROI-anchored CDF, and the
  //  Gaussian itself.  The ROI spans mean +- 10 channels either side of +-3.5 sigma.
  roi.lower_energy = roi.mean - 3.5*roi.sigma;
  roi.upper_energy = roi.mean + 3.5*roi.sigma;

  const double f_lo = gauss_cdf( roi.lower_energy, roi.mean, roi.sigma );
  const double f_hi = gauss_cdf( roi.upper_energy, roi.mean, roi.sigma );

  vector<float> spec_counts( spec_nchannel, 0.0f );
  for( size_t i = 0; i < spec_nchannel; ++i )
  {
    const double lo = spec_energies[i], hi = spec_energies[i+1];
    const double center = 0.5*(lo + hi);
    const double cdfbar = (std::min)( (std::max)( gauss_cdf(center, roi.mean, roi.sigma), f_lo ), f_hi ) - f_lo;

    const double cont = (roi.continuum_density + (roi.step_density * cdfbar)) * (hi - lo);
    const double peak = roi.amplitude * (gauss_cdf(hi, roi.mean, roi.sigma) - gauss_cdf(lo, roi.mean, roi.sigma));
    spec_counts[i] = static_cast<float>( cont + peak );
  }

  shared_ptr<SpecUtils::EnergyCalibration> cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_lower_channel_energy( spec_nchannel, spec_energies );

  roi.data = make_shared<SpecUtils::Measurement>();
  roi.data->set_gamma_counts( make_shared<vector<float>>( spec_counts ), 1.0f, 1.0f );
  roi.data->set_energy_calibration( cal );

  // Snap the ROI to channel edges, the way the fitters do, and copy out its channel range so the
  //  direct fit_amp_and_offset_imp(...) calls below see exactly the channels PeakFitLM would.
  const size_t lower_channel = roi.data->find_gamma_channel( static_cast<float>(roi.lower_energy) );
  const size_t upper_channel = roi.data->find_gamma_channel( static_cast<float>(roi.upper_energy) );
  roi.nchannel = 1 + upper_channel - lower_channel;

  roi.energies.assign( spec_energies.begin() + lower_channel,
                       spec_energies.begin() + lower_channel + roi.nchannel + 1 );
  roi.counts.assign( spec_counts.begin() + lower_channel,
                     spec_counts.begin() + lower_channel + roi.nchannel );

  roi.lower_energy = roi.energies.front();
  roi.upper_energy = roi.energies.back();
  roi.ref_energy = roi.lower_energy;

  return roi;
}


/** Runs `fit_amp_and_offset_imp` over the ROI, then rebuilds the fitted peak+continuum and asks
 `PeakContinuum::offset_integral(...)` to evaluate it.  Returns the fitter's per-channel model in
 `fitter_counts` and the evaluator's in `evaluator_counts`.
 */
void fit_and_evaluate( const SyntheticRoi &roi,
                       const PeakContinuum::OffsetType type,
                       const vector<double> &step_coeffs,
                       const PeakDef::SkewType skew_type,
                       const double *skew_pars,
                       vector<double> &fitter_counts,
                       vector<double> &evaluator_counts,
                       vector<double> &fit_continuum_coeffs,
                       double &fit_amplitude,
                       const vector<PeakDef> &fixed_amp_peaks = {} )
{
  const vector<double> means{ roi.mean };
  const vector<double> sigmas{ roi.sigma };
  const vector<PeakDef> &no_fixed_peaks = fixed_amp_peaks;

  vector<double> amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts;
  fitter_counts.assign( roi.nchannel, 0.0 );

  PeakFit::fit_amp_and_offset_imp<PeakDef,double>(
      roi.energies.data(), roi.counts.data(), nullptr, roi.nchannel, type,
      step_coeffs.empty() ? nullptr : step_coeffs.data(), roi.ref_energy,
      means, sigmas, no_fixed_peaks, skew_type, skew_pars,
      amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts, fitter_counts.data() );

  BOOST_REQUIRE_EQUAL( amplitudes.size(), size_t(1) );
  fit_amplitude = amplitudes[0];

  // The solver returns only the polynomial terms; append the step coefficients we handed it.
  for( const double sc : step_coeffs )
    continuum_coeffs.push_back( sc );
  fit_continuum_coeffs = continuum_coeffs;

  shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
  cont->setType( type );
  cont->setRange( roi.lower_energy, roi.upper_energy );
  cont->setParameters( roi.ref_energy, continuum_coeffs, {} );

  PeakDef peak( roi.mean, roi.sigma, fit_amplitude );
  peak.setSkewType( skew_type );
  if( skew_pars )
  {
    for( size_t i = 0; i < PeakDef::num_skew_parameters(skew_type); ++i )
      peak.set_coefficient( skew_pars[i], static_cast<PeakDef::CoefficientType>(PeakDef::SkewPar0 + i) );
  }
  peak.setContinuum( cont );

  // The evaluator sums the step over EVERY peak in the ROI, fixed-amplitude ones included.
  vector<PeakDef> all_peaks{ peak };
  all_peaks.insert( all_peaks.end(), fixed_amp_peaks.begin(), fixed_amp_peaks.end() );

  vector<const PeakDef *> peak_ptrs;
  for( const PeakDef &p : all_peaks )
    peak_ptrs.push_back( &p );

  evaluator_counts.assign( roi.nchannel, 0.0 );
  cont->offset_integral( roi.energies.data(), evaluator_counts.data(), roi.nchannel,
                         roi.data, peak_ptrs.data(), peak_ptrs.size() );
  for( PeakDef &p : all_peaks )
    p.gauss_integral( roi.energies.data(), evaluator_counts.data(), roi.nchannel );
}
}//namespace


// The single most valuable test here: the model the fitter optimises must be the model InterSpec
//  draws and integrates.  Guards two independent ways they used to diverge:
//   - the fitter clamped only the polynomial at zero while the evaluator clamped polynomial+step,
//     so the fitted continuum could run deeply negative where the drawn one was floored at 0;
//   - the fitter's CDF was anchored at the ROI's lower edge while the evaluator's ran from
//     -infinity, offsetting the drawn continuum by step_coeff*SUM(amp_j*CDF_j(roi_lower)).
BOOST_AUTO_TEST_CASE( fitter_and_evaluator_agree_channel_by_channel )
{
  // `rel_tol`: the CDF step types agree to floating point by construction, because the evaluator
  //  averages the channel's edge CDFs exactly as unit_pdf_to_cdf's cumulative sum does.  The
  //  data-based step types are given a looser bound: their fitter and evaluator build `frac_data`
  //  by slightly different quadratures, a small pre-existing inconsistency this work did not
  //  change and does not address.
  struct Case { PeakContinuum::OffsetType type; vector<double> step_coeffs; const char *name;
                double rel_tol; };

  const SyntheticRoi roi = make_synthetic_roi();
  const double true_step_coeff = roi.step_density / roi.amplitude;

  const vector<Case> cases{
    { PeakContinuum::FlatStepCDF,     { true_step_coeff },         "FlatStepCDF",              1.0E-6 },
    { PeakContinuum::LinearStepCDF,   { true_step_coeff },         "LinearStepCDF",            1.0E-6 },
    { PeakContinuum::BiLinearStepCDF, { true_step_coeff, 1.0E-5 }, "BiLinearStepCDF",          1.0E-6 },
    // A large negative step drives the fitted continuum negative on the high side - the case that
    //  used to make the two disagree the most.
    { PeakContinuum::FlatStepCDF,     { -1.0E-1 },                 "FlatStepCDF (extreme step)",1.0E-6 },
    { PeakContinuum::Linear,          { },                         "Linear (must be unaffected)",1.0E-6 },
    { PeakContinuum::FlatStep,        { },                         "FlatStep (must be unaffected)",1.0E-4 },
  };

  for( const Case &c : cases )
  {
    vector<double> fitter_counts, evaluator_counts, cont_coeffs;
    double fit_amp = 0.0;
    fit_and_evaluate( roi, c.type, c.step_coeffs, PeakDef::SkewType::NoSkew, nullptr,
                      fitter_counts, evaluator_counts, cont_coeffs, fit_amp );

    BOOST_REQUIRE_EQUAL( fitter_counts.size(), evaluator_counts.size() );
    for( size_t i = 0; i < fitter_counts.size(); ++i )
    {
      const double tol = (std::max)( c.rel_tol * fabs(fitter_counts[i]), 1.0E-6 );
      BOOST_CHECK_MESSAGE( fabs(fitter_counts[i] - evaluator_counts[i]) < tol,
        c.name << ": channel " << i << " fitter=" << fitter_counts[i]
               << " evaluator=" << evaluator_counts[i]
               << " (diff " << (fitter_counts[i] - evaluator_counts[i]) << ")" );
    }
  }//for( const Case &c : cases )
}


// Same check with a low-energy skew tail, where the un-anchored CDF used to be worst: F0 can be
//  several percent for a Bortel tail rather than the ~2e-4 of a bare Gaussian.
BOOST_AUTO_TEST_CASE( fitter_and_evaluator_agree_with_skew )
{
  const SyntheticRoi roi = make_synthetic_roi();
  const double true_step_coeff = roi.step_density / roi.amplitude;
  const double skew_pars[1] = { 0.35 };   // Bortel skew: relative tail amplitude

  vector<double> fitter_counts, evaluator_counts, cont_coeffs;
  double fit_amp = 0.0;
  fit_and_evaluate( roi, PeakContinuum::FlatStepCDF, { true_step_coeff },
                    PeakDef::SkewType::Bortel, skew_pars,
                    fitter_counts, evaluator_counts, cont_coeffs, fit_amp );

  for( size_t i = 0; i < fitter_counts.size(); ++i )
  {
    const double tol = (std::max)( 1.0E-6 * fabs(fitter_counts[i]), 1.0E-6 );
    BOOST_CHECK_MESSAGE( fabs(fitter_counts[i] - evaluator_counts[i]) < tol,
      "Bortel: channel " << i << " fitter=" << fitter_counts[i]
                         << " evaluator=" << evaluator_counts[i] );
  }
}


// With the correct step coefficient the noise-free ROI must recover its injected peak area; with a
// step magnitude carried over from a non-CDF continuum - which is what setType() used to do - the
// amplitude least-squares solve drives the peak to ~0 and inflates the continuum.  That collapse is
// the mechanism behind both reported symptoms, so assert it directly rather than describing it.
BOOST_AUTO_TEST_CASE( amplitude_survives_correct_step_and_collapses_on_absurd_one )
{
  const SyntheticRoi roi = make_synthetic_roi();
  const double true_step_coeff = roi.step_density / roi.amplitude;

  vector<double> fitter_counts, evaluator_counts, cont_coeffs;
  double fit_amp = 0.0;

  fit_and_evaluate( roi, PeakContinuum::FlatStepCDF, { true_step_coeff },
                    PeakDef::SkewType::NoSkew, nullptr,
                    fitter_counts, evaluator_counts, cont_coeffs, fit_amp );

  BOOST_CHECK_MESSAGE( fabs(fit_amp - roi.amplitude) < 0.02*roi.amplitude,
    "Fitted amplitude " << fit_amp << " should be within 2% of " << roi.amplitude );

  BOOST_CHECK_MESSAGE( fabs(cont_coeffs[0] - roi.continuum_density) < 0.05*roi.continuum_density,
    "Fitted continuum density " << cont_coeffs[0] << " should be within 5% of "
    << roi.continuum_density );

  // Now the value setType() used to carry across: the non-CDF step magnitude (counts/keV) dropped
  //  into the CDF coefficient's keV^-1 slot, i.e. larger than correct by the ROI's peak area.
  const double carried_over = roi.step_density;
  double collapsed_amp = 0.0;
  vector<double> collapsed_cont;
  fit_and_evaluate( roi, PeakContinuum::FlatStepCDF, { carried_over },
                    PeakDef::SkewType::NoSkew, nullptr,
                    fitter_counts, evaluator_counts, collapsed_cont, collapsed_amp );

  BOOST_CHECK_MESSAGE( collapsed_amp < 0.01*roi.amplitude,
    "A step coefficient of " << carried_over << " (the non-CDF magnitude) should collapse the peak,"
    " but the amplitude came back as " << collapsed_amp );

  BOOST_CHECK_MESSAGE( collapsed_cont[0] > 2.0*roi.continuum_density,
    "...and should inflate the continuum well above " << roi.continuum_density
    << ", but it came back as " << collapsed_cont[0] );
}


namespace
{
/** Builds a single FlatStepCDF/LinearStepCDF/BiLinearStepCDF peak over the synthetic ROI, with
 the step coefficients left at zero - i.e. what `setType()` hands the fitter after a safe type
 change.
 */
shared_ptr<PeakDef> make_unfit_cdf_peak( const SyntheticRoi &roi,
                                         const PeakContinuum::OffsetType type )
{
  shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
  cont->setType( type );
  cont->setRange( roi.lower_energy, roi.upper_energy );

  vector<double> coeffs( PeakContinuum::num_parameters(type), 0.0 );
  coeffs[0] = roi.continuum_density;
  cont->setParameters( roi.ref_energy, coeffs, {} );

  shared_ptr<PeakDef> peak = make_shared<PeakDef>( roi.mean, roi.sigma, roi.amplitude );
  peak->setContinuum( cont );

  return peak;
}
}//namespace


// PeakFitLM used to seed the step coefficient at exactly whatever was stored - zero, after a type
//  change - with bounds of +-1e4.  chi2 is shallow enough in this parameter that starting at zero
//  usually finished at zero, quietly turning a FlatStepCDF into a plain Constant.  The warm start
//  pre-solves it by least-squares; the bounds are now the dimensionless +-2.
BOOST_AUTO_TEST_CASE( peakfitlm_recovers_step_from_zero_seed )
{
  const SyntheticRoi roi = make_synthetic_roi();
  const double true_step_coeff = roi.step_density / roi.amplitude;

  for( const PeakContinuum::OffsetType type : { PeakContinuum::FlatStepCDF,
                                                PeakContinuum::LinearStepCDF } )
  {
    const shared_ptr<PeakDef> peak = make_unfit_cdf_peak( roi, type );
    BOOST_REQUIRE( peak->continuum()->parameters().back() == 0.0 );

    vector<shared_ptr<const PeakDef>> results;
    PeakFitLM::fit_peaks_LM( results, { peak }, roi.data, 0.0, 0.0, true,
                             PeakFitUtils::CoarseResolutionType::High );

    const string ctx = PeakContinuum::offset_type_str( type );
    BOOST_REQUIRE_MESSAGE( results.size() == 1, ctx << ": fit returned " << results.size()
                                                    << " peaks instead of 1" );

    const shared_ptr<const PeakContinuum> cont = results[0]->continuum();
    const size_t num_poly = PeakContinuum::num_linear_fit_pars( type );
    const double fit_step = cont->parameters()[num_poly];

    BOOST_CHECK_MESSAGE( fit_step != 0.0, ctx << ": step coefficient stayed at exactly zero" );
    BOOST_CHECK_MESSAGE( fabs(fit_step - true_step_coeff) < 0.25*fabs(true_step_coeff),
      ctx << ": step coefficient " << fit_step << " should be within 25% of " << true_step_coeff );

    // The peak must survive, with something close to its injected area.
    BOOST_CHECK_MESSAGE( fabs(results[0]->amplitude() - roi.amplitude) < 0.05*roi.amplitude,
      ctx << ": amplitude " << results[0]->amplitude() << " should be within 5% of "
          << roi.amplitude );
  }//for( type )
}


// Holding a continuum coefficient fixed sends the ROI down PeakFitLM's non-LLS path, which had two
// defects: it threw from inside the Ceres cost function for CDF step continua (so the caller's
// catch cleared the result and the ROI's peaks silently disappeared), and - for EVERY continuum
// type, not just these - it never solved the peak amplitudes, leaving them at the internal 1.0
// placeholder the least-squares solve would otherwise have replaced.
//
// The amplitudes are now Ceres parameters whenever the LLS is not solving them, so a pinned
// coefficient must come back with both the pin honoured and a real peak area.
BOOST_AUTO_TEST_CASE( pinned_continuum_coefficient_still_fits_amplitude )
{
  const SyntheticRoi roi = make_synthetic_roi();

  // Non-CDF types are included deliberately: the amplitude defect was never specific to the CDF
  //  step continua, and Linear/FlatStep returned amplitude 1 just the same.
  for( const PeakContinuum::OffsetType type : { PeakContinuum::Linear,
                                                PeakContinuum::Quadratic,
                                                PeakContinuum::FlatStep,
                                                PeakContinuum::FlatStepCDF,
                                                PeakContinuum::LinearStepCDF,
                                                PeakContinuum::BiLinearStepCDF } )
  {
    const shared_ptr<PeakDef> peak = make_unfit_cdf_peak( roi, type );
    shared_ptr<PeakContinuum> cont = peak->continuum();

    // Pin the constant term at the value the synthetic ROI was built with.
    BOOST_REQUIRE( cont->setPolynomialCoefFitFor( 0, false ) );

    const string ctx = PeakContinuum::offset_type_str( type );

    vector<shared_ptr<const PeakDef>> results;
    BOOST_REQUIRE_NO_THROW( PeakFitLM::fit_peaks_LM( results, { peak }, roi.data, 0.0, 0.0, true,
                                                     PeakFitUtils::CoarseResolutionType::High ) );

    BOOST_REQUIRE_MESSAGE( results.size() == 1,
      ctx << " with a pinned continuum coefficient: got " << results.size() << " peaks, not 1" );

    // The pin must be honoured...
    BOOST_CHECK_MESSAGE( fabs(results[0]->continuum()->parameters()[0] - roi.continuum_density) < 1.0E-6,
      ctx << ": pinned coefficient moved from " << roi.continuum_density << " to "
          << results[0]->continuum()->parameters()[0] );

    // ...and the amplitude must be a real fit, not the 1.0 placeholder.
    BOOST_CHECK_MESSAGE( fabs(results[0]->amplitude() - roi.amplitude) < 0.05*roi.amplitude,
      ctx << ": amplitude " << results[0]->amplitude() << " should be within 5% of "
          << roi.amplitude );
  }//for( type )
}


// Pinning only a peak-CDF *step* coefficient does not need the non-LLS path at all: the step
// coefficients are Ceres parameters either way, so the polynomial terms stay least-squares solved
// and the pinned step simply becomes a constant parameter.
BOOST_AUTO_TEST_CASE( pinned_step_coefficient_keeps_lls_path )
{
  const SyntheticRoi roi = make_synthetic_roi();
  const double pinned_step = roi.step_density / roi.amplitude;

  for( const PeakContinuum::OffsetType type : { PeakContinuum::FlatStepCDF,
                                                PeakContinuum::LinearStepCDF } )
  {
    const size_t num_poly = PeakContinuum::num_linear_fit_pars( type );

    shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
    cont->setType( type );
    cont->setRange( roi.lower_energy, roi.upper_energy );
    vector<double> coeffs( PeakContinuum::num_parameters(type), 0.0 );
    coeffs[0] = roi.continuum_density;
    coeffs[num_poly] = pinned_step;
    cont->setParameters( roi.ref_energy, coeffs, {} );
    BOOST_REQUIRE( cont->setPolynomialCoefFitFor( num_poly, false ) );

    shared_ptr<PeakDef> peak = make_shared<PeakDef>( roi.mean, roi.sigma, roi.amplitude );
    peak->setContinuum( cont );

    const string ctx = PeakContinuum::offset_type_str( type );

    vector<shared_ptr<const PeakDef>> results;
    BOOST_REQUIRE_NO_THROW( PeakFitLM::fit_peaks_LM( results, { peak }, roi.data, 0.0, 0.0, true,
                                                     PeakFitUtils::CoarseResolutionType::High ) );
    BOOST_REQUIRE_MESSAGE( results.size() == 1, ctx << ": got " << results.size() << " peaks, not 1" );

    BOOST_CHECK_MESSAGE( fabs(results[0]->continuum()->parameters()[num_poly] - pinned_step)
                         < 1.0E-12,
      ctx << ": pinned step coefficient moved from " << pinned_step << " to "
          << results[0]->continuum()->parameters()[num_poly] );

    BOOST_CHECK_MESSAGE( fabs(results[0]->amplitude() - roi.amplitude) < 0.05*roi.amplitude,
      ctx << ": amplitude " << results[0]->amplitude() << " should be within 5% of "
          << roi.amplitude );
  }//for( type )
}


// Serialization version 3 changed what BiLinearStepCDF's four parameters mean, from two blended
// lines (left_const, left_linear, right_const, right_linear) to (const, linear, step0, step1).
// Files written before it must still describe the continuum they were saved with, so
// PeakContinuum::fromXml converts on read - which it can only do because the <Peak> nodes carrying
// the ROI's peak area sit alongside the <PeakContinuum> node in the same document.
BOOST_AUTO_TEST_CASE( legacy_bilinear_step_cdf_converts_on_read )
{
  const double total_amp = 20000.0;
  const double left0 = 300.0, left1 = -0.5, right0 = 270.0, right1 = -0.4;

  shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
  cont->setType( PeakContinuum::BiLinearStepCDF );
  cont->setRange( 590.0, 610.0 );
  cont->setParameters( 590.0, { left0, left1, right0, right1 }, {} );

  PeakDef peak( 600.0, 1.0, total_amp );
  peak.setContinuum( cont );

  // Serialize the way SpecMeas does: continuum and peak as siblings under a <Peaks> node.
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *peaks_node = doc.allocate_node( rapidxml::node_element, "Peaks" );
  doc.append_node( peaks_node );

  map<shared_ptr<PeakContinuum>,int> continuum_ids;
  peak.toXml( peaks_node, peaks_node, continuum_ids );

  rapidxml::xml_node<char> *cont_node = peaks_node->first_node( "PeakContinuum" );
  BOOST_REQUIRE( cont_node );
  BOOST_REQUIRE( peaks_node->first_node( "Peak" ) );

  // Pretend it was written by the pre-reparameterisation code.
  rapidxml::xml_attribute<char> *ver = cont_node->first_attribute( "version" );
  BOOST_REQUIRE( ver );
  ver->value( doc.allocate_string("2") );

  shared_ptr<PeakContinuum> read_cont = make_shared<PeakContinuum>();
  int cont_id = 0;
  read_cont->fromXml( cont_node, cont_id );

  const vector<double> &pars = read_cont->parameters();
  BOOST_REQUIRE_EQUAL( pars.size(), size_t(4) );

  // Tolerances are percentages; PeakContinuum::fromXml parses coefficients into a float.
  BOOST_CHECK_CLOSE( pars[0], left0, 1.0E-3 );
  BOOST_CHECK_CLOSE( pars[1], left1, 1.0E-3 );
  BOOST_CHECK_CLOSE( pars[2], (right0 - left0)/total_amp, 1.0E-3 );
  BOOST_CHECK_CLOSE( pars[3], (right1 - left1)/total_amp, 1.0E-3 );

  // A current-version node must be taken at face value, not converted a second time.
  ver->value( doc.allocate_string("3") );
  shared_ptr<PeakContinuum> modern = make_shared<PeakContinuum>();
  modern->fromXml( cont_node, cont_id );
  BOOST_CHECK_CLOSE( modern->parameters()[2], right0, 1.0E-3 );
  BOOST_CHECK_CLOSE( modern->parameters()[3], right1, 1.0E-3 );
}


// PeakFitChi2Fcn packs the shared-continuum index and the OffsetType into a single double that it
// hands Minuit2 as a parameter.  The OffsetType field was one decimal digit wide, so the two types
// numbered 10 and 11 - BiLinearStepCDF and External - decoded back as NoOffset, tripping the
// round-trip check in addPeaksToFitter(...) and aborting the fit for any ROI using them.
// This is the revival of the commented-out `PeakFitChi2Fcn::testOffsetConversions()`.
BOOST_AUTO_TEST_CASE( continuum_info_encoding_round_trips_every_type )
{
  for( const PeakContinuum::OffsetType type : all_offset_types() )
  {
    for( const int index : { -1, 0, 1, 2, 37, 500, 9998 } )
    {
      const string ctx = string(PeakContinuum::offset_type_str(type)) + " / index " + to_string(index);

      // Both orders of setting the two fields must work.
      double info = 0.0;
      PeakFitChi2Fcn::setSharedIndexToContinuumInfo( info, index );
      PeakFitChi2Fcn::setOffsetTypeToContinuumInfo( info, type );
      BOOST_CHECK_MESSAGE( PeakFitChi2Fcn::continuumInfoToSharedIndex(info) == index,
        ctx << ": index decoded as " << PeakFitChi2Fcn::continuumInfoToSharedIndex(info) );
      BOOST_CHECK_MESSAGE( PeakFitChi2Fcn::continuumInfoToOffsetType(info) == type,
        ctx << ": type decoded as "
            << PeakContinuum::offset_type_str( PeakFitChi2Fcn::continuumInfoToOffsetType(info) ) );

      info = 0.0;
      PeakFitChi2Fcn::setOffsetTypeToContinuumInfo( info, type );
      PeakFitChi2Fcn::setSharedIndexToContinuumInfo( info, index );
      BOOST_CHECK_MESSAGE( PeakFitChi2Fcn::continuumInfoToSharedIndex(info) == index,
        ctx << " (reversed order): index decoded as "
            << PeakFitChi2Fcn::continuumInfoToSharedIndex(info) );
      BOOST_CHECK_MESSAGE( PeakFitChi2Fcn::continuumInfoToOffsetType(info) == type,
        ctx << " (reversed order): type decoded as "
            << PeakContinuum::offset_type_str( PeakFitChi2Fcn::continuumInfoToOffsetType(info) ) );
    }//for( index )
  }//for( type )
}


// BiLinearStepCDF was reparameterised from two blended lines
//   (1-f)*(a0 + a1*E') + f*(b0 + b1*E'),   f = SUM_j(amp_j*CDFbar_j) / SUM_j(amp_j)
// to a linear continuum plus a linearly varying step coefficient
//   p0 + p1*E' + (s0 + s1*E')*SUM_j(amp_j*CDFbar_j)
// These span exactly the same set of shapes: substituting f = g/A and collecting terms gives
//   p0 = a0,  p1 = a1,  s0 = (b0-a0)/A,  s1 = (b1-a1)/A
// and the inverse b_k = p_k + s_k*A.  The equivalence survives channel integration because
// 0.5*(x1^2 - x0^2) == center*dx exactly, so folding E' into the step coefficient at the channel
// centre reproduces the exactly-integrated right-hand line.
//
// This test checks that identity numerically, to floating point, rather than trusting the algebra.
BOOST_AUTO_TEST_CASE( bilinear_step_cdf_reparameterization_is_exact )
{
  const SyntheticRoi roi = make_synthetic_roi();

  PeakDef peak( roi.mean, roi.sigma, roi.amplitude );
  const double total_amp = roi.amplitude;

  // Probe the evaluator with (p0,p1,s0,s1) = (0,0,1,0) to recover SUM_j(amp_j*CDFbar_j)*dx per
  //  channel, without having to re-implement the anchored CDF here.
  shared_ptr<PeakContinuum> probe = make_shared<PeakContinuum>();
  probe->setType( PeakContinuum::BiLinearStepCDF );
  probe->setRange( roi.lower_energy, roi.upper_energy );
  probe->setParameters( roi.ref_energy, { 0.0, 0.0, 1.0, 0.0 }, {} );
  peak.setContinuum( probe );

  const PeakDef *peak_ptr = &peak;
  vector<double> g_dx( roi.nchannel, 0.0 );
  probe->offset_integral( roi.energies.data(), g_dx.data(), roi.nchannel, roi.data, &peak_ptr, 1 );

  // An arbitrary pair of left/right lines, with a genuine slope difference.
  const double a0 = 300.0, a1 = -0.8, b0 = 268.0, b1 = 0.35;

  const double s0 = (b0 - a0)/total_amp;
  const double s1 = (b1 - a1)/total_amp;

  shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
  cont->setType( PeakContinuum::BiLinearStepCDF );
  cont->setRange( roi.lower_energy, roi.upper_energy );
  cont->setParameters( roi.ref_energy, { a0, a1, s0, s1 }, {} );
  peak.setContinuum( cont );

  vector<double> new_form( roi.nchannel, 0.0 );
  cont->offset_integral( roi.energies.data(), new_form.data(), roi.nchannel, roi.data, &peak_ptr, 1 );

  for( size_t i = 0; i < roi.nchannel; ++i )
  {
    const double x0_rel = roi.energies[i] - roi.ref_energy;
    const double x1_rel = roi.energies[i+1] - roi.ref_energy;
    const double dx = x1_rel - x0_rel;

    // The old form, evaluated the way the old code did: both lines integrated exactly across the
    //  channel, blended by the normalised CDF fraction.
    const double frac = (g_dx[i]/dx) / total_amp;
    const double left_int  = (a0 * dx) + (0.5 * a1 * ((x1_rel*x1_rel) - (x0_rel*x0_rel)));
    const double right_int = (b0 * dx) + (0.5 * b1 * ((x1_rel*x1_rel) - (x0_rel*x0_rel)));
    const double old_form = (std::max)( 0.0, ((1.0 - frac) * left_int) + (frac * right_int) );

    const double tol = (std::max)( 1.0E-11 * fabs(old_form), 1.0E-11 );
    BOOST_CHECK_MESSAGE( fabs(new_form[i] - old_form) < tol,
      "channel " << i << ": new form " << new_form[i] << " != old form " << old_form
                 << " (diff " << (new_form[i] - old_form) << ", frac=" << frac << ")" );
  }//for( size_t i = 0; i < roi.nchannel; ++i )

  // (The inverse map b_k = p_k + s_k*A is exact by construction of s_k above, so asserting it here
  //  would only be testing the test; what matters is the per-channel agreement checked above.)
}


namespace
{
/** The `Continuum_Type` value a pre-version-3 InterSpec would have written - i.e. no "(vN)" tag. */
const char * const s_legacy_bilinear_cdf_type = "BiLinearStepCDF";

/** Builds a peak CSV of the subset of columns `csv_to_candidate_fit_peaks` actually reads.

 `cont_coefs` is the whole `Continuum_Coefficients` cell: reference energy then the coefficients,
 space separated, exactly as `write_peak_csv` lays it out.
 */
struct CsvPeakRow
{
  double mean = 0.0, fwhm = 0.0, area = 0.0;
  double roi_lower = 0.0, roi_upper = 0.0;
  string cont_type;
  string cont_coefs;
  string skew_type = "None";
  string skew_coefs;
};


string make_peak_csv( const vector<CsvPeakRow> &rows )
{
  stringstream csv;
  csv << "Centroid,Net_Area,FWHM,ROI_Lower_Energy,ROI_Upper_Energy,Continuum_Type,Skew_Type,"
         "Continuum_Coefficients,Skew_Coefficients\r\n";

  csv.precision( 12 );
  for( const CsvPeakRow &r : rows )
  {
    csv << r.mean << ',' << r.area << ',' << r.fwhm << ',' << r.roi_lower << ',' << r.roi_upper
        << ',' << r.cont_type << ',' << r.skew_type << ',' << r.cont_coefs << ',' << r.skew_coefs
        << "\r\n";
  }

  return csv.str();
}


/** Formats a `Continuum_Coefficients` cell. */
string cont_coef_cell( const double ref_energy, const vector<double> &coefs )
{
  stringstream ss;
  ss.precision( 12 );
  ss << ref_energy;
  for( const double c : coefs )
    ss << ' ' << c;
  return ss.str();
}


/** Loads `csv` against `roi`'s spectrum, and returns the single continuum every peak shares. */
shared_ptr<const PeakContinuum> load_single_roi_csv( const SyntheticRoi &roi, const string &csv,
                                                     vector<PeakDef> &peaks )
{
  stringstream strm( csv );
  peaks = PeakModel::csv_to_candidate_fit_peaks( roi.data, strm );
  BOOST_REQUIRE( !peaks.empty() );

  shared_ptr<const PeakContinuum> cont = peaks.front().continuum();
  BOOST_REQUIRE( cont );
  for( const PeakDef &p : peaks )
    BOOST_REQUIRE_MESSAGE( p.continuum() == cont, "CSV peaks did not end up sharing one ROI" );

  return cont;
}


/** Evaluates `cont` over `roi`'s channels, with `peaks` as the ROI's peaks. */
vector<double> evaluate_over_roi( const SyntheticRoi &roi,
                                  const shared_ptr<const PeakContinuum> &cont,
                                  const vector<PeakDef> &peaks )
{
  vector<const PeakDef *> peak_ptrs;
  for( const PeakDef &p : peaks )
    peak_ptrs.push_back( &p );

  vector<double> counts( roi.nchannel, 0.0 );
  cont->offset_integral( roi.energies.data(), counts.data(), roi.nchannel, roi.data,
                         peak_ptrs.data(), peak_ptrs.size() );
  return counts;
}


/** The per-channel `SUM_j( amp_j * CDFbar_j ) * dx` the version-3 model multiplies its step by.

 Recovered by probing the evaluator with `(p0,p1,s0,s1) = (0,0,1,0)`, rather than re-implementing
 the anchored CDF here - the same trick `bilinear_step_cdf_reparameterization_is_exact` uses.
 */
vector<double> probe_anchored_cdf_sum( const SyntheticRoi &roi, vector<PeakDef> peaks )
{
  shared_ptr<PeakContinuum> probe = make_shared<PeakContinuum>();
  probe->setType( PeakContinuum::BiLinearStepCDF );
  probe->setRange( roi.lower_energy, roi.upper_energy );
  probe->setParameters( roi.ref_energy, { 0.0, 0.0, 1.0, 0.0 }, {} );
  for( PeakDef &p : peaks )
    p.setContinuum( probe );

  return evaluate_over_roi( roi, probe, peaks );
}


/** `SUM_j( amp_j * CDF_j(roi_lower) )` - the offset the un-anchored version-2 CDF carried. */
double roi_amp_cdf0( const SyntheticRoi &roi, const vector<PeakDef> &peaks )
{
  double amp_cdf0 = 0.0;
  for( const PeakDef &p : peaks )
  {
    const double *skew_pars = p.coefficients() + PeakDef::CoefficientType::SkewPar0;
    amp_cdf0 += p.amplitude() * PeakDists::peak_cdf( roi.lower_energy, p.mean(), p.sigma(),
                                                     p.skewType(), skew_pars );
  }
  return amp_cdf0;
}


/** The version-2 model's counts in channel `i`, built from the same CDF samples the version-3
 evaluator uses so that only the parameter map is under test (see the reparameterisation comment
 above `bilinear_step_cdf_reparameterization_is_exact`).

 Version 2 was `(1-f)*a(E) + f*b(E)` with `f = g_raw/A`, which rearranges to
 `a(E) + g_raw(E)*(s0 + s1*E')`; `g_raw = gbar + amp_cdf0` since version 2's CDF ran from -infinity.
 */
double legacy_channel_counts( const SyntheticRoi &roi, const size_t i,
                              const double gbar_dx, const double amp_cdf0, const double total_amp,
                              const double a0, const double a1, const double b0, const double b1 )
{
  const double x0_rel = roi.energies[i] - roi.ref_energy;
  const double x1_rel = roi.energies[i+1] - roi.ref_energy;
  const double dx = x1_rel - x0_rel;

  const double s0 = (b0 - a0)/total_amp;
  const double s1 = (b1 - a1)/total_amp;

  const double left_int = (a0 * dx) + (0.5 * a1 * ((x1_rel*x1_rel) - (x0_rel*x0_rel)));
  // `s1` multiplies E' at the channel centre; 0.5*(x1^2 - x0^2) == center*dx exactly.
  const double step_int = (s0 * dx + 0.5 * s1 * ((x1_rel*x1_rel) - (x0_rel*x0_rel)))
                          * ((gbar_dx/dx) + amp_cdf0);

  return (std::max)( 0.0, left_int + step_int );
}
}//namespace


// The peak CSV carries no file-level version field, so a `BiLinearStepCDF` cell with no "(vN)" tag
// is a pre-version-3 one and its four coefficients mean `(left_const, left_linear, right_const,
// right_linear)`.  Read under the current meaning instead, slot 2 - a count density of O(100) -
// would become a step coefficient whose natural scale is O(1e-3), i.e. the absurd step that
// `amplitude_survives_correct_step_and_collapses_on_absurd_one` shows collapsing the peak.
//
// Unlike the XML path, there is no continuum id to gather the ROI's peaks by; they are rows of the
// same file, so the conversion can only run once every row has been read.
BOOST_AUTO_TEST_CASE( legacy_bilinear_step_cdf_csv_converts_on_read )
{
  const SyntheticRoi roi = make_synthetic_roi();

  // Two peaks sharing the ROI, so the conversion's `SUM_j(amp_j)` is genuinely a sum over rows.
  const double amp0 = roi.amplitude, amp1 = 0.4*roi.amplitude;
  const double total_amp = amp0 + amp1;
  const double a0 = 300.0, a1 = -0.8, b0 = 268.0, b1 = 0.35;

  CsvPeakRow row;
  row.fwhm = 2.35482*roi.sigma;
  row.roi_lower = roi.lower_energy;
  row.roi_upper = roi.upper_energy;
  row.cont_type = s_legacy_bilinear_cdf_type;
  row.cont_coefs = cont_coef_cell( roi.ref_energy, { a0, a1, b0, b1 } );

  CsvPeakRow row0 = row, row1 = row;
  row0.mean = roi.mean;             row0.area = amp0;
  row1.mean = roi.mean + 1.5;       row1.area = amp1;

  vector<PeakDef> peaks;
  const shared_ptr<const PeakContinuum> cont
                  = load_single_roi_csv( roi, make_peak_csv( { row0, row1 } ), peaks );

  BOOST_REQUIRE_EQUAL( peaks.size(), size_t(2) );
  BOOST_REQUIRE( cont->type() == PeakContinuum::BiLinearStepCDF );

  const vector<double> &pars = cont->parameters();
  BOOST_REQUIRE_EQUAL( pars.size(), size_t(4) );

  // The step coefficients must be the differences over the ROI's total area - not the raw
  //  right-hand line, which is what an ungated read would have left here.
  BOOST_CHECK_CLOSE( pars[2], (b0 - a0)/total_amp, 1.0E-3 );
  BOOST_CHECK_CLOSE( pars[3], (b1 - a1)/total_amp, 1.0E-3 );
  BOOST_CHECK_MESSAGE( fabs(pars[2]) < 1.0,
    "step coefficient " << pars[2] << " is a count density, so the legacy coefficients were"
    " reinterpreted rather than converted" );

  // What actually matters: the continuum draws the same shape it did under version 2.
  const vector<double> loaded = evaluate_over_roi( roi, cont, peaks );
  const vector<double> gbar_dx = probe_anchored_cdf_sum( roi, peaks );
  const double amp_cdf0 = roi_amp_cdf0( roi, peaks );

  for( size_t i = 0; i < roi.nchannel; ++i )
  {
    const double expected = legacy_channel_counts( roi, i, gbar_dx[i], amp_cdf0, total_amp,
                                                   a0, a1, b0, b1 );
    // Coefficients go through a float in the CSV reader, hence 1e-6 rather than 1e-11.
    const double tol = (std::max)( 1.0E-6 * fabs(expected), 1.0E-6 );
    BOOST_CHECK_MESSAGE( fabs(loaded[i] - expected) < tol,
      "channel " << i << ": loaded " << loaded[i] << " != legacy shape " << expected
                 << " (diff " << (loaded[i] - expected) << ")" );
  }
}


// The version-2 CDF ran from -infinity, so it carried a constant `step_k * SUM_j(amp_j*CDF_j(roi
// lower))` that the ROI-anchored version-3 model does not.  For a bare Gaussian that offset is
// ~2e-4 and invisible; for a Bortel tail it is several percent, so this is the case that actually
// exercises the fold-in into the polynomial terms - in either the XML or the CSV path.
BOOST_AUTO_TEST_CASE( legacy_bilinear_step_cdf_csv_conversion_includes_cdf_anchor )
{
  const SyntheticRoi roi = make_synthetic_roi();

  const double total_amp = roi.amplitude;
  const double a0 = 300.0, a1 = -0.8, b0 = 268.0, b1 = 0.35;
  const double bortel_tau = 0.35;

  CsvPeakRow row;
  row.mean = roi.mean;
  row.area = total_amp;
  row.fwhm = 2.35482*roi.sigma;
  row.roi_lower = roi.lower_energy;
  row.roi_upper = roi.upper_energy;
  row.cont_type = s_legacy_bilinear_cdf_type;
  row.cont_coefs = cont_coef_cell( roi.ref_energy, { a0, a1, b0, b1 } );
  row.skew_type = PeakDef::to_string( PeakDef::SkewType::Bortel );
  row.skew_coefs = to_string( bortel_tau );

  vector<PeakDef> peaks;
  const shared_ptr<const PeakContinuum> cont
                  = load_single_roi_csv( roi, make_peak_csv( { row } ), peaks );

  BOOST_REQUIRE_EQUAL( peaks.size(), size_t(1) );
  BOOST_REQUIRE( peaks.front().skewType() == PeakDef::SkewType::Bortel );

  const double amp_cdf0 = roi_amp_cdf0( roi, peaks );
  const double s0 = (b0 - a0)/total_amp;

  // The offset has to be non-negligible, or this test is checking nothing.  Note that
  //  `PeakDists::peak_cdf` is not normalised the same way for every skew type - Bortel's runs over
  //  [-0.5,0.5] rather than [0,1] - but that cancels: version 2 blended with this same function, so
  //  the conversion reproduces whatever offset it carried.
  BOOST_REQUIRE_MESSAGE( fabs(amp_cdf0/total_amp) > 0.01,
    "Bortel CDF offset at the ROI edge is only " << (amp_cdf0/total_amp) );

  // The constant term must have absorbed the offset, rather than being carried over as-is.
  const vector<double> &pars = cont->parameters();
  BOOST_CHECK_CLOSE( pars[0], a0 + s0*amp_cdf0, 1.0E-3 );
  BOOST_CHECK_MESSAGE( fabs(pars[0] - a0) > 1.0E-6*fabs(a0),
    "constant term " << pars[0] << " is unchanged from the legacy left_const " << a0
                     << ", so the CDF-anchor offset was dropped" );

  const vector<double> loaded = evaluate_over_roi( roi, cont, peaks );
  const vector<double> gbar_dx = probe_anchored_cdf_sum( roi, peaks );

  double worst_without_foldin = 0.0;
  for( size_t i = 0; i < roi.nchannel; ++i )
  {
    const double expected = legacy_channel_counts( roi, i, gbar_dx[i], amp_cdf0, total_amp,
                                                   a0, a1, b0, b1 );
    const double tol = (std::max)( 1.0E-6 * fabs(expected), 1.0E-6 );
    BOOST_CHECK_MESSAGE( fabs(loaded[i] - expected) < tol,
      "Bortel: channel " << i << ": loaded " << loaded[i] << " != legacy shape " << expected
                         << " (diff " << (loaded[i] - expected) << ")" );

    // How wrong the shape would be had the fold-in been skipped.
    const double no_foldin = legacy_channel_counts( roi, i, gbar_dx[i], 0.0, total_amp,
                                                    a0, a1, b0, b1 );
    if( expected > 0.0 )
      worst_without_foldin = (std::max)( worst_without_foldin, fabs(no_foldin - expected)/expected );
  }

  BOOST_CHECK_MESSAGE( worst_without_foldin > 0.005,
    "dropping the CDF-anchor offset would only move the continuum by "
      << 100.0*worst_without_foldin << "%, so this test does not pin it" );
}


// The writer tags the current convention, so a round trip must come back unchanged - in particular
// not converted a second time, which would divide the step coefficients by the ROI area again.
BOOST_AUTO_TEST_CASE( bilinear_step_cdf_csv_round_trips )
{
  const SyntheticRoi roi = make_synthetic_roi();

  const double s0 = -30.0/roi.amplitude, s1 = 1.5e-5;
  const vector<double> coefs{ 300.0, -0.8, s0, s1 };

  shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
  cont->setType( PeakContinuum::BiLinearStepCDF );
  cont->setRange( roi.lower_energy, roi.upper_energy );
  cont->setParameters( roi.ref_energy, coefs, {} );

  shared_ptr<PeakDef> peak = make_shared<PeakDef>( roi.mean, roi.sigma, roi.amplitude );
  peak->setContinuum( cont );

  deque<shared_ptr<const PeakDef>> peaks{ peak };
  stringstream csv;
  PeakModel::write_peak_csv( csv, "synthetic", PeakModel::PeakCsvType::Full, peaks, roi.data );

  BOOST_REQUIRE_MESSAGE( csv.str().find("BiLinearStepCDF(v3)") != string::npos,
    "write_peak_csv did not tag the continuum convention:\n" << csv.str() );

  stringstream in( csv.str() );
  const vector<PeakDef> read_peaks = PeakModel::csv_to_candidate_fit_peaks( roi.data, in );
  BOOST_REQUIRE_EQUAL( read_peaks.size(), size_t(1) );

  const shared_ptr<const PeakContinuum> read_cont = read_peaks.front().continuum();
  BOOST_REQUIRE( read_cont->type() == PeakContinuum::BiLinearStepCDF );

  const vector<double> &pars = read_cont->parameters();
  BOOST_REQUIRE_EQUAL( pars.size(), size_t(4) );
  for( size_t i = 0; i < 4; ++i )
    BOOST_CHECK_CLOSE( pars[i], coefs[i], 1.0E-3 );
}


// A CSV claiming a convention this InterSpec does not know must not have its coefficients read
// under todays meaning; the row falls back to the continuum derived from the data.
BOOST_AUTO_TEST_CASE( future_continuum_convention_is_not_read_as_current )
{
  const SyntheticRoi roi = make_synthetic_roi();

  CsvPeakRow row;
  row.mean = roi.mean;
  row.area = roi.amplitude;
  row.fwhm = 2.35482*roi.sigma;
  row.roi_lower = roi.lower_energy;
  row.roi_upper = roi.upper_energy;
  row.cont_type = string(s_legacy_bilinear_cdf_type) + "(v99)";
  row.cont_coefs = cont_coef_cell( roi.ref_energy, { 300.0, -0.8, 268.0, 0.35 } );

  stringstream strm( make_peak_csv( { row } ) );
  const vector<PeakDef> peaks = PeakModel::csv_to_candidate_fit_peaks( roi.data, strm );
  BOOST_REQUIRE_EQUAL( peaks.size(), size_t(1) );

  const shared_ptr<const PeakContinuum> cont = peaks.front().continuum();
  BOOST_CHECK_MESSAGE( cont->type() != PeakContinuum::BiLinearStepCDF,
    "a future continuum convention was accepted as the current one" );
}


// A ROI holding one peak's amplitude fixed while fitting another.  The fixed peak contributes
// `amp*CDFbar` to the step exactly as the fitted one does - the evaluator and fit_continuum both
// sum over every peak in the ROI - so leaving it out of fit_amp_and_offset_imp's step basis makes
// the fitted continuum disagree with the drawn one.  That is the same class of divergence as the
// clamp and the CDF anchoring, and it only shows up when an ROI actually has a fixed-amp peak.
BOOST_AUTO_TEST_CASE( fitter_and_evaluator_agree_with_fixed_amplitude_peak )
{
  const SyntheticRoi roi = make_synthetic_roi();
  const double true_step_coeff = roi.step_density / roi.amplitude;

  // A second, weaker peak inside the ROI whose amplitude is not being fit.
  PeakDef fixed_peak( roi.mean + 2.0, roi.sigma, 0.25*roi.amplitude );
  const vector<PeakDef> fixed_amp_peaks{ fixed_peak };

  const vector<pair<PeakContinuum::OffsetType,vector<double>>> cases{
    { PeakContinuum::FlatStepCDF,     { true_step_coeff } },
    { PeakContinuum::LinearStepCDF,   { true_step_coeff } },
    { PeakContinuum::BiLinearStepCDF, { true_step_coeff, 1.0E-5 } },
  };

  for( const auto &c : cases )
  {
    vector<double> fitter_counts, evaluator_counts, cont_coeffs;
    double fit_amp = 0.0;
    fit_and_evaluate( roi, c.first, c.second, PeakDef::SkewType::NoSkew, nullptr,
                      fitter_counts, evaluator_counts, cont_coeffs, fit_amp, fixed_amp_peaks );

    const string ctx = PeakContinuum::offset_type_str( c.first );
    BOOST_REQUIRE_EQUAL( fitter_counts.size(), evaluator_counts.size() );
    for( size_t i = 0; i < fitter_counts.size(); ++i )
    {
      const double tol = (std::max)( 1.0E-6 * fabs(fitter_counts[i]), 1.0E-6 );
      BOOST_CHECK_MESSAGE( fabs(fitter_counts[i] - evaluator_counts[i]) < tol,
        ctx << " with fixed-amp peak: channel " << i << " fitter=" << fitter_counts[i]
            << " evaluator=" << evaluator_counts[i]
            << " (diff " << (fitter_counts[i] - evaluator_counts[i]) << ")" );
    }
  }//for( const auto &c : cases )
}


// The per-ROI Ceres parameter block gained an amplitude section (present only when the continuum is
// not being LLS-solved), which sits between the means and the per-ROI skew block.  The writer in
// setup_roi_parameters and the reader in process_one_roi must agree on where that skew block
// starts; when they did not, the skew seeds silently overwrote the amplitude seeds and the
// evaluator read zeros as skew - a failure that is invisible in Debug because every index still
// lands inside the ROI's own block.
//
// IndependentSkewValues is what puts a skew block inside each ROI, so this combination - pinned
// continuum coefficient + per-ROI skew + a skew type with parameters - is the one that pins the
// offsets against each other.
BOOST_AUTO_TEST_CASE( pinned_coefficient_with_independent_skew_keeps_blocks_aligned )
{
  const SyntheticRoi roi = make_synthetic_roi();

  for( const PeakContinuum::OffsetType type : { PeakContinuum::Linear,
                                                PeakContinuum::FlatStepCDF } )
  {
    const shared_ptr<PeakDef> peak = make_unfit_cdf_peak( roi, type );
    peak->setSkewType( PeakDef::SkewType::Bortel );
    peak->set_coefficient( 0.35, PeakDef::CoefficientType::SkewPar0 );
    peak->setFitFor( PeakDef::CoefficientType::SkewPar0, true );
    BOOST_REQUIRE( peak->continuum()->setPolynomialCoefFitFor( 0, false ) );

    const string ctx = PeakContinuum::offset_type_str( type );

    Wt::WFlags<PeakFitLM::PeakFitLMOptions> options;
    options |= PeakFitLM::PeakFitLMOptions::IndependentSkewValues;

    vector<shared_ptr<const PeakDef>> results;
    BOOST_REQUIRE_NO_THROW( results = PeakFitLM::fit_peaks_in_roi_LM( { peak }, roi.data,
                                          PeakFitUtils::CoarseResolutionType::High, options ) );

    BOOST_REQUIRE_MESSAGE( results.size() == 1,
      ctx << " (pinned coefficient + IndependentSkewValues): got " << results.size()
          << " peaks, not 1" );

    // The amplitude block must not have been clobbered by the skew block.
    BOOST_CHECK_MESSAGE( fabs(results[0]->amplitude() - roi.amplitude) < 0.10*roi.amplitude,
      ctx << ": amplitude " << results[0]->amplitude() << " should be within 10% of "
          << roi.amplitude );

    // ...and the skew parameter must be a real value, not the zero left behind by a misaligned read.
    const double skew0 = results[0]->coefficient( PeakDef::CoefficientType::SkewPar0 );
    BOOST_CHECK_MESSAGE( std::isfinite(skew0) && (skew0 > 0.0),
      ctx << ": Bortel skew parameter came back as " << skew0
          << ", which is what a misaligned skew block reads" );
  }//for( type )
}
