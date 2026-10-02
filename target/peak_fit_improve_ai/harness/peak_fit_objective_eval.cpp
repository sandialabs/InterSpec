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

/** peak_fit_objective_eval: compares the statistical behaviour of PeakFitLM's fit objectives - the
 modified-Neyman chi2, the default fit (chi2, with sparse ROIs refit by Poisson maximum likelihood),
 and the likelihood fit forced on every ROI - on Poisson-sampled spectra with known truth.  The
 likelihood fit is IRLS, or all-in-Ceres when PeakFitLM.cpp is built with
 `SPARSE_DATA_LIKELIHOOD_USE_CERES` set to 1.

 Every objective is run on the identical Poisson draw with the identical starting values, so the
 comparisons are paired.  For each (problem, peak, objective) it reports the area bias, the pull
 distribution (fit - reference)/reported_uncertainty, coverage, failure rate, the tail of bad fits,
 and CPU time; the same for the peak mean and FWHM.

 Two sources of problems:
  - `--mode=synthetic`: spectra generated from InterSpec's own peak + continuum model, so the truth
    is exact (no model mismatch).  A one-factor-at-a-time grid over detector class, FWHM in
    channels, peak area, continuum level, continuum type, skew and peak topology.
  - `--mode=inject`: the GADRAS-injected spectra of `peak_fit_accuracy_inject_compact`; PCF record 3
    is the noise-free source+background spectrum, so any number of Poisson replicas can be drawn.
    Bias is reported against both the GADRAS truth area and a "pseudo-truth" - the chi2 fit to the
    noise-free spectrum scaled up 1000x - which removes peak-shape model mismatch from the
    statistical comparison.

 Typical use (Release build):
   peak_fit_objective_eval --mode=synthetic --out=/tmp/objeval_syn --replicas=50
   peak_fit_objective_eval --mode=inject --out=/tmp/objeval_inj --replicas=10 --dets=Detective-X

 Run with --help for all options.
 */

#include "InterSpec_config.h"

#include <map>
#include <set>
#include <cmath>
#include <ctime>
#include <mutex>
#include <array>
#include <atomic>
#include <chrono>
#include <random>
#include <optional>
#include <string>
#include <thread>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include <Wt/WFlags.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"

#include "FitPeaksCorpusScore.h"

using namespace std;

namespace
{

// ------------------------------------------------------------------------------------------------
// Portable random numbers: copied from target/testing/test_DetectionLimit.cpp, because
//  std::poisson_distribution differs between standard libraries, and a run must be reproducible.
// ------------------------------------------------------------------------------------------------
double study_uniform( mt19937 &generator )
{
  return ( static_cast<double>( generator() ) + 0.5 ) / 4294967296.0;
}


int study_poisson( const double mean, mt19937 &generator )
{
  const double sm_knuth_mean_limit = 10.0;

  if( !(mean > 0.0) )
    return 0;

  if( mean < sm_knuth_mean_limit )
  {
    const double limit = std::exp( -mean );
    double product = 1.0;
    int count = 0;
    for( ; count < 10000; ++count )
    {
      product *= study_uniform( generator );
      if( product <= limit )
        break;
    }
    return count;
  }//if( mean < sm_knuth_mean_limit )

  // Hoermann (1993), "The transformed rejection method for generating Poisson random variables".
  const double b = 0.931 + 2.53 * std::sqrt( mean );
  const double a = -0.059 + 0.02483 * b;
  const double inverse_alpha = 1.1239 + 1.1328 / ( b - 3.4 );
  const double v_r = 0.9277 - 3.6224 / ( b - 2.0 );

  for( size_t iteration = 0; iteration < 10000; ++iteration )
  {
    const double u = study_uniform( generator ) - 0.5;
    const double v = study_uniform( generator );
    const double us = 0.5 - std::fabs( u );
    const double k = std::floor( ( 2.0 * a / us + b ) * u + mean + 0.43 );

    if( (us >= 0.07) && (v <= v_r) )
      return static_cast<int>( k );

    if( (k < 0.0) || ((us < 0.013) && (v > us)) )
      continue;

    const double log_mean = std::log( mean );
    const double lhs = std::log( v * inverse_alpha / ( a / ( us * us ) + b ) );
    const double rhs = -mean + k * log_mean - std::lgamma( k + 1.0 );
    if( lhs <= rhs )
      return static_cast<int>( k );
  }

  return static_cast<int>( std::floor( mean + 0.5 ) );  // unreachable in practice
}//study_poisson(...)


/** FNV-1a; unlike std::hash, stable across platforms and runs. */
uint64_t stable_hash( const string &s )
{
  uint64_t h = 14695981039346656037ULL;
  for( const char c : s )
  {
    h ^= static_cast<unsigned char>( c );
    h *= 1099511628211ULL;
  }
  return h;
}


uint32_t seed_for( const string &key, const size_t replica, const char *purpose, const uint64_t run_seed )
{
  const uint64_t h = stable_hash( key + "|" + std::to_string( replica ) + "|" + purpose
                                  + "|" + std::to_string( run_seed ) );
  return static_cast<uint32_t>( h ^ (h >> 32) );
}


double thread_cpu_seconds()
{
  timespec ts;
  clock_gettime( CLOCK_THREAD_CPUTIME_ID, &ts );
  return static_cast<double>( ts.tv_sec ) + 1.0E-9*static_cast<double>( ts.tv_nsec );
}


void parallel_for( const size_t count, const int threads, const std::function<void(size_t)> &fcn )
{
  const int nthreads = std::max( 1, std::min( threads, static_cast<int>(count) ) );
  if( nthreads == 1 )
  {
    for( size_t i = 0; i < count; ++i )
      fcn( i );
    return;
  }

  std::atomic<size_t> next{ 0 };
  vector<std::thread> workers;
  for( int t = 0; t < nthreads; ++t )
  {
    workers.emplace_back( [&](){
      while( true )
      {
        const size_t i = next.fetch_add( 1 );
        if( i >= count )
          break;
        fcn( i );
      }
    } );
  }
  for( std::thread &w : workers )
    w.join();
}//parallel_for(...)


double poisson_deviance_term( const double n, const double m )
{
  const double mm = std::max( m, 1.0E-12 );
  if( n <= 0.0 )
    return 2.0*mm;
  return 2.0*( mm - n + n*std::log( n/mm ) );
}


// ------------------------------------------------------------------------------------------------
// Objectives
// ------------------------------------------------------------------------------------------------
struct Objective
{
  string name;
  Wt::WFlags<PeakFitLM::PeakFitLMOptions> flags;
  /** Report the area uncertainty conditional on the fitted shapes (amplitude over
   `PeakFitLM::peak_detection_significance` with chi2 weights) instead of the fit's marginal one. */
  bool conditional_uncert = false;
};


vector<Objective> available_objectives()
{
  vector<Objective> answer;
  // "chi2" is the plain modified-Neyman fit; "default" is what a caller passing no options gets (chi2,
  //  with sparse ROIs refit by Poisson maximum likelihood).
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> neyman( PeakFitLM::PeakFitLMOptions::NoSparseDataLikelihood );
  answer.push_back( { "chi2", neyman } );
  answer.push_back( { "default", Wt::WFlags<PeakFitLM::PeakFitLMOptions>() } );
  // "likelihood": every ROI refit by Poisson maximum likelihood.
  answer.push_back( { "likelihood", Wt::WFlags<PeakFitLM::PeakFitLMOptions>( PeakFitLM::PeakFitLMOptions::ForcePoissonLikelihood ) } );

  // The chi2 fit with area uncertainties conditional on the means/widths/skew (the pre-2026-10
  //  reported uncertainty), for comparing uncertainty calibration.
  answer.push_back( { "chi2-cond", neyman, true } );
  return answer;
}


// ------------------------------------------------------------------------------------------------
// Problems
// ------------------------------------------------------------------------------------------------
struct TruthPeak
{
  double mean = 0.0, fwhm = 0.0, area = 0.0;
  double z = 0.0;          // S/sqrt(S+B), B within +-1 FWHM
  bool signal = true;      // inject: source line (else background line)

  // Per objective: that objective's fit to the noise-free spectrum scaled up by `sm_pseudo_scale`
  //  (areas divided back down); NaN if that fit failed.  The value the objective's estimate tends to
  //  with unlimited statistics - so the reference for its statistical bias and pulls - and, compared
  //  between objectives, a measure of how differently they respond to peak-shape model mismatch.
  vector<double> pseudo_area, pseudo_mean, pseudo_fwhm;

  double pseudo_area_for( const size_t obj ) const
  {
    return (obj < pseudo_area.size()) ? pseudo_area[obj] : std::numeric_limits<double>::quiet_NaN();
  }
  double pseudo_mean_for( const size_t obj ) const
  {
    return (obj < pseudo_mean.size()) ? pseudo_mean[obj] : std::numeric_limits<double>::quiet_NaN();
  }
  double pseudo_fwhm_for( const size_t obj ) const
  {
    return (obj < pseudo_fwhm.size()) ? pseudo_fwhm[obj] : std::numeric_limits<double>::quiet_NaN();
  }
};


struct Problem
{
  string id;
  vector<pair<string,string>> dims;        // report dimensions; same keys, in order, for every problem in a run
  PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::Unknown;
  double roi_lower = 0.0, roi_upper = 0.0;
  PeakContinuum::OffsetType fit_cont_type = PeakContinuum::OffsetType::Linear;
  PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;
  vector<double> skew_pars;                // truth (synthetic) / starting (inject) values
  bool fit_skew = true;
  vector<TruthPeak> truth;                 // sorted by mean
  shared_ptr<const SpecUtils::Measurement> expected;  // noise-free expected counts
  string spectrum_key;                     // problems sharing a spectrum share their Poisson draws
};


const double sm_pseudo_scale = 1000.0;


shared_ptr<SpecUtils::Measurement> measurement_with_counts( const shared_ptr<const SpecUtils::Measurement> &like,
                                                            vector<float> &&counts )
{
  auto meas = make_shared<SpecUtils::Measurement>( *like );
  const float lt = (like->live_time() > 0.0f) ? like->live_time() : 1.0f;
  const float rt = (like->real_time() > 0.0f) ? like->real_time() : lt;
  meas->set_gamma_counts( make_shared<const vector<float>>( std::move(counts) ), lt, rt );
  return meas;
}


shared_ptr<const SpecUtils::Measurement> poisson_sample( const shared_ptr<const SpecUtils::Measurement> &expected,
                                                         const double scale, mt19937 &rng )
{
  const vector<float> &mu = *expected->gamma_counts();
  vector<float> counts( mu.size(), 0.0f );
  for( size_t i = 0; i < mu.size(); ++i )
    counts[i] = static_cast<float>( study_poisson( scale*std::max( 0.0f, mu[i] ), rng ) );
  return measurement_with_counts( expected, std::move(counts) );
}


shared_ptr<const SpecUtils::Measurement> scaled_spectrum( const shared_ptr<const SpecUtils::Measurement> &expected,
                                                          const double scale )
{
  const vector<float> &mu = *expected->gamma_counts();
  vector<float> counts( mu.size(), 0.0f );
  for( size_t i = 0; i < mu.size(); ++i )
    counts[i] = static_cast<float>( scale*std::max( 0.0f, mu[i] ) );
  return measurement_with_counts( expected, std::move(counts) );
}


/** Starting peaks for one fit: all share one new continuum of the problem's fit type, over the ROI.

 With `perturb`, means move by up to +-0.25 sigma, the FWHM (one common factor) by +-10%, the
 amplitudes by +-30%, and fitted skew parameters by +-10% (clamped to their range).
 */
vector<shared_ptr<const PeakDef>> make_start_peaks( const Problem &prob, const double area_scale,
                                                    const bool perturb, mt19937 &rng )
{
  auto cont = make_shared<PeakContinuum>();
  cont->setType( prob.fit_cont_type );
  cont->setRange( prob.roi_lower, prob.roi_upper );

  const double fwhm_factor = perturb ? (0.9 + 0.2*study_uniform( rng )) : 1.0;

  vector<shared_ptr<const PeakDef>> peaks;
  for( const TruthPeak &t : prob.truth )
  {
    const double sigma = fwhm_factor * t.fwhm / 2.35482;
    const double mean = t.mean + (perturb ? (study_uniform( rng ) - 0.5)*0.5*(t.fwhm/2.35482) : 0.0);
    const double amp = area_scale * t.area * (perturb ? (0.7 + 0.6*study_uniform( rng )) : 1.0);

    auto p = make_shared<PeakDef>( mean, sigma, std::max( amp, 1.0 ) );
    p->setSkewType( prob.skew_type );
    const size_t nskew = PeakDef::num_skew_parameters( prob.skew_type );
    for( size_t i = 0; i < nskew; ++i )
    {
      const PeakDef::CoefficientType ct = static_cast<PeakDef::CoefficientType>( PeakDef::SkewPar0 + i );
      double val = (i < prob.skew_pars.size()) ? prob.skew_pars[i] : 0.0;
      if( perturb && prob.fit_skew )
      {
        double lower = 0.0, upper = 0.0, start = 0.0, step = 0.0;
        PeakDef::skew_parameter_range( prob.skew_type, ct, lower, upper, start, step );
        val *= (0.9 + 0.2*study_uniform( rng ));
        val = std::min( std::max( val, lower ), upper );
      }
      p->set_coefficient( val, ct );
      p->setFitFor( ct, prob.fit_skew );
    }
    p->setContinuum( cont );
    peaks.push_back( p );
  }//for( const TruthPeak &t : prob.truth )

  return peaks;
}//make_start_peaks(...)


/** Model counts of fitted peaks (continuum + peaks) in channels [ch0, ch1]. */
vector<double> model_counts( const vector<shared_ptr<const PeakDef>> &peaks,
                             const shared_ptr<const SpecUtils::Measurement> &data,
                             const size_t ch0, const size_t ch1 )
{
  const size_t nchan = ch1 - ch0 + 1;
  vector<double> model( nchan, 0.0 );
  if( peaks.empty() )
    return model;

  const shared_ptr<const vector<float>> &energies = data->channel_energies();
  const float *x = energies->data() + ch0;

  vector<const PeakDef *> roi_peaks;
  for( const shared_ptr<const PeakDef> &p : peaks )
    roi_peaks.push_back( p.get() );

  peaks.front()->continuum()->offset_integral( x, model.data(), nchan, data,
                                               roi_peaks.data(), roi_peaks.size() );
  for( const shared_ptr<const PeakDef> &p : peaks )
    p->gauss_integral( x, model.data(), nchan );

  return model;
}//model_counts(...)


// ------------------------------------------------------------------------------------------------
// Synthetic problems
// ------------------------------------------------------------------------------------------------
struct DetClass
{
  string name;
  PeakFitUtils::CoarseResolutionType type;
  double fwhm;   // keV at the 661.657 keV test energy
};


const DetClass sm_det_hpge{ "HPGe", PeakFitUtils::CoarseResolutionType::High, 1.9 };
const DetClass sm_det_labr{ "LaBr", PeakFitUtils::CoarseResolutionType::LaBr, 20.0 };
const DetClass sm_det_nai{ "NaI", PeakFitUtils::CoarseResolutionType::Low, 46.0 };


string fmt_num( const double v )
{
  ostringstream strm;
  strm << v;
  return strm.str();
}


/** `cont` is the fitted continuum name: NoOffset, Constant, Linear, Quadratic, FlatStep (generated
 as FlatStepCDF), FlatStepCDF or LinearStepCDF.  `topo` is single, doublet11 or doublet51 (second
 peak 1.5 FWHM above, with equal or one-fifth the area).
 */
Problem make_synthetic_problem( const string &suite, const DetClass &det, const double fwhm_ch,
                                const double area, const double cont_per_channel, const string &cont,
                                const PeakDef::SkewType skew, const bool fit_skew, const string &topo )
{
  const double e0 = 661.657;
  const double fwhm = det.fwhm;
  const double sigma = fwhm / 2.35482;
  const double chan_width = fwhm / fwhm_ch;

  Problem prob;
  prob.det_type = det.type;
  prob.skew_type = skew;
  prob.fit_skew = fit_skew;

  vector<pair<double,double>> mean_areas{ { e0, area } };
  if( topo == "doublet11" )
    mean_areas.push_back( { e0 + 1.5*fwhm, area } );
  else if( topo == "doublet51" )
    mean_areas.push_back( { e0 + 1.5*fwhm, 0.2*area } );
  else if( topo != "single" )
    throw runtime_error( "Unknown topology '" + topo + "'" );

  prob.roi_lower = mean_areas.front().first - 4.0*fwhm;
  prob.roi_upper = mean_areas.back().first + 3.5*fwhm;
  const double roi_width = prob.roi_upper - prob.roi_lower;

  // Skew truth values: moderate tails, well inside the fittable ranges.
  switch( skew )
  {
    case PeakDef::SkewType::NoSkew:
      break;
    case PeakDef::SkewType::GaussExp:
      prob.skew_pars = { 1.2 };
      break;
    case PeakDef::SkewType::DoubleSidedCrystalBall:
      prob.skew_pars = { 1.0, 3.0, 2.0, 6.0 };
      break;
    default:
      throw runtime_error( "make_synthetic_problem: unsupported skew type" );
  }//switch( skew )

  // The spectrum spans the ROI plus a margin; channel edges are deliberately not aligned to the ROI.
  const double spec_lower = prob.roi_lower - 3.0*fwhm - 0.37*chan_width;
  const size_t nchannel = static_cast<size_t>( std::ceil( (roi_width + 6.0*fwhm)/chan_width ) );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( nchannel, { static_cast<float>(spec_lower), static_cast<float>(chan_width) }, {} );

  // Truth continuum, as a density in counts/keV relative to e0.
  const double c0 = cont_per_channel / chan_width;
  double total_area = 0.0;
  for( const pair<double,double> &ma : mean_areas )
    total_area += ma.second;

  PeakContinuum::OffsetType gen_type = PeakContinuum::OffsetType::NoOffset;
  vector<double> gen_pars;
  if( cont == "NoOffset" || (cont_per_channel <= 0.0) )
  {
    gen_type = PeakContinuum::OffsetType::NoOffset;
  }else if( cont == "Constant" )
  {
    gen_type = PeakContinuum::OffsetType::Constant;
    gen_pars = { c0 };
  }else if( cont == "Linear" )
  {
    gen_type = PeakContinuum::OffsetType::Linear;
    gen_pars = { c0, -0.3*c0/roi_width };
  }else if( cont == "Quadratic" )
  {
    gen_type = PeakContinuum::OffsetType::Quadratic;
    gen_pars = { c0, -0.3*c0/roi_width, 0.6*c0/(roi_width*roi_width) };
  }else if( cont == "FlatStep" || cont == "FlatStepCDF" )
  {
    // density = p0 + s0*SUM(amp*CDFbar); the low-energy side is 25% higher than the high side.
    gen_type = PeakContinuum::OffsetType::FlatStepCDF;
    gen_pars = { 1.125*c0, -0.25*c0/total_area };
  }else if( cont == "LinearStepCDF" )
  {
    gen_type = PeakContinuum::OffsetType::LinearStepCDF;
    gen_pars = { 1.125*c0, -0.2*c0/roi_width, -0.25*c0/total_area };
  }else
  {
    throw runtime_error( "Unknown continuum '" + cont + "'" );
  }

  auto gen_cont = make_shared<PeakContinuum>();
  gen_cont->setType( gen_type );
  gen_cont->setRange( prob.roi_lower, prob.roi_upper );
  if( !gen_pars.empty() )
    gen_cont->setParameters( e0, gen_pars, {} );

  vector<shared_ptr<PeakDef>> gen_peaks;
  for( const pair<double,double> &ma : mean_areas )
  {
    auto p = make_shared<PeakDef>( ma.first, sigma, ma.second );
    p->setSkewType( skew );
    for( size_t i = 0; i < prob.skew_pars.size(); ++i )
      p->set_coefficient( prob.skew_pars[i], static_cast<PeakDef::CoefficientType>( PeakDef::SkewPar0 + i ) );
    p->setContinuum( gen_cont );
    gen_peaks.push_back( p );
  }

  const shared_ptr<const vector<float>> energies = cal->channel_energies();
  vector<double> expected( nchannel, 0.0 );
  vector<const PeakDef *> peak_ptrs;
  for( const shared_ptr<PeakDef> &p : gen_peaks )
    peak_ptrs.push_back( p.get() );
  if( gen_type != PeakContinuum::OffsetType::NoOffset )
    gen_cont->offset_integral( energies->data(), expected.data(), nchannel, nullptr,
                               peak_ptrs.data(), peak_ptrs.size() );
  for( const shared_ptr<PeakDef> &p : gen_peaks )
    p->gauss_integral( energies->data(), expected.data(), nchannel );

  auto meas = make_shared<SpecUtils::Measurement>();
  meas->set_gamma_counts( make_shared<const vector<float>>( expected.begin(), expected.end() ), 1.0f, 1.0f );
  meas->set_energy_calibration( cal );
  prob.expected = meas;

  // Truth peaks, with z = S/sqrt(S+B), B being the continuum within +-1 FWHM of the mean.
  for( const pair<double,double> &ma : mean_areas )
  {
    TruthPeak t;
    t.mean = ma.first;
    t.fwhm = fwhm;
    t.area = ma.second;
    const double cont_in_window = (gen_type == PeakContinuum::OffsetType::NoOffset)
                                  ? 0.0 : (c0 * 2.0 * fwhm);
    t.z = (0.76*t.area) / std::sqrt( std::max( 0.76*t.area + cont_in_window, 1.0E-9 ) );
    prob.truth.push_back( t );
  }

  // The fit continuum type.
  if( cont == "NoOffset" )
    prob.fit_cont_type = PeakContinuum::OffsetType::NoOffset;
  else if( cont == "Constant" )
    prob.fit_cont_type = PeakContinuum::OffsetType::Constant;
  else if( cont == "Linear" )
    prob.fit_cont_type = PeakContinuum::OffsetType::Linear;
  else if( cont == "Quadratic" )
    prob.fit_cont_type = PeakContinuum::OffsetType::Quadratic;
  else if( cont == "FlatStep" )
    prob.fit_cont_type = PeakContinuum::OffsetType::FlatStep;
  else if( cont == "FlatStepCDF" )
    prob.fit_cont_type = PeakContinuum::OffsetType::FlatStepCDF;
  else if( cont == "LinearStepCDF" )
    prob.fit_cont_type = PeakContinuum::OffsetType::LinearStepCDF;

  const string skew_str = PeakDef::to_string( skew );
  const string skew_fit_str = (skew == PeakDef::SkewType::NoSkew) ? "-" : (fit_skew ? "fit" : "fixed");
  prob.dims = {
    { "suite", suite }, { "det", det.name }, { "fwhm_ch", fmt_num(fwhm_ch) }, { "area", fmt_num(area) },
    { "cont_level", fmt_num(cont_per_channel) }, { "cont", cont }, { "skew", skew_str },
    { "skew_fit", skew_fit_str }, { "topo", topo }
  };

  prob.id = "syn";
  for( size_t i = 1; i < prob.dims.size(); ++i )   // the suite is not part of the identity
    prob.id += "_" + prob.dims[i].second;
  prob.spectrum_key = prob.id;

  return prob;
}//make_synthetic_problem(...)


/** One-factor-at-a-time synthetic design: a primary grid of (detector class/FWHM, area, continuum
 level) on a Linear continuum, then continuum-type, skew and topology variations on a reduced grid.
 */
vector<Problem> build_synthetic_problems( const set<string> &suites )
{
  vector<Problem> problems;
  set<string> seen;
  const auto add = [&]( Problem &&p ){
    if( seen.insert( p.id ).second )
      problems.push_back( std::move(p) );
  };

  const PeakDef::SkewType noskew = PeakDef::SkewType::NoSkew;

  if( suites.count( "primary" ) )
  {
    const vector<pair<DetClass,double>> det_fwhm{
      { sm_det_hpge, 2.5 }, { sm_det_hpge, 6.0 }, { sm_det_hpge, 15.0 },
      { sm_det_labr, 6.0 }, { sm_det_nai, 15.0 }
    };
    for( const pair<DetClass,double> &df : det_fwhm )
    {
      for( const double area : { 10.0, 30.0, 100.0, 300.0, 1.0E3, 1.0E4, 1.0E5 } )
      {
        for( const double level : { 0.0, 0.1, 1.0, 10.0, 100.0, 1000.0 } )
        {
          if( level <= 0.0 )
          {
            add( make_synthetic_problem( "primary", df.first, df.second, area, 0.0, "NoOffset", noskew, true, "single" ) );
            add( make_synthetic_problem( "primary", df.first, df.second, area, 0.0, "Constant", noskew, true, "single" ) );
          }else
          {
            add( make_synthetic_problem( "primary", df.first, df.second, area, level, "Linear", noskew, true, "single" ) );
          }
        }//for( level )
      }//for( area )
    }//for( det_fwhm )
  }//if( primary )

  const vector<pair<DetClass,double>> reduced_det{ { sm_det_hpge, 6.0 }, { sm_det_nai, 15.0 } };
  const vector<double> reduced_areas{ 30.0, 300.0, 1.0E4 };
  const vector<double> reduced_levels{ 0.1, 1.0, 10.0, 100.0 };

  if( suites.count( "cont" ) )
  {
    for( const pair<DetClass,double> &df : reduced_det )
      for( const double area : reduced_areas )
        for( const double level : reduced_levels )
          for( const char *cont : { "Constant", "Quadratic", "FlatStep", "LinearStepCDF" } )
            add( make_synthetic_problem( "cont", df.first, df.second, area, level, cont, noskew, true, "single" ) );
  }//if( cont )

  if( suites.count( "skew" ) )
  {
    for( const pair<DetClass,double> &df : reduced_det )
      for( const double area : reduced_areas )
        for( const double level : reduced_levels )
          for( const PeakDef::SkewType skew : { PeakDef::SkewType::GaussExp, PeakDef::SkewType::DoubleSidedCrystalBall } )
            for( const bool fit_skew : { true, false } )
              add( make_synthetic_problem( "skew", df.first, df.second, area, level, "Linear", skew, fit_skew, "single" ) );
  }//if( skew )

  if( suites.count( "topo" ) )
  {
    for( const pair<DetClass,double> &df : reduced_det )
      for( const double area : reduced_areas )
        for( const double level : reduced_levels )
          for( const char *topo : { "doublet11", "doublet51" } )
            add( make_synthetic_problem( "topo", df.first, df.second, area, level, "Linear", noskew, true, topo ) );
  }//if( topo )

  return problems;
}//build_synthetic_problems(...)


// ------------------------------------------------------------------------------------------------
// Inject problems
// ------------------------------------------------------------------------------------------------

// Copies of the mappings in target/peak_fit_improve (ClassifyDetType_GA / PeakFitImproveData), so
//  this tool does not need to link openGA.
PeakFitUtils::CoarseResolutionType true_det_type_for_name( const string &name )
{
  using PeakFitUtils::CoarseResolutionType;
  for( const char *s : { "Falcon", "Fulcrum", "LANL_X", "Detective", "HPGe_Planar" } )
    if( SpecUtils::icontains( name, s ) )
      return CoarseResolutionType::High;
  for( const char *s : { "1.5cm-2cm", "1cm-1cm", "CZT_H3D", "Kromek-GR1", "nanoRaider", "Interceptor", "Raider",
                         "InSpector 1000", "SAM-Eagle-LaBr", "LaBr3", "IdentiFINDER-LaBr", "Radseeker-LaBr" } )
    if( SpecUtils::icontains( name, s ) )
      return CoarseResolutionType::MedRes;
  return CoarseResolutionType::Low;
}


const char *det_type_str( const PeakFitUtils::CoarseResolutionType type )
{
  switch( type )
  {
    case PeakFitUtils::CoarseResolutionType::Low:         return "Low";
    case PeakFitUtils::CoarseResolutionType::LaBr:        return "LaBr";
    case PeakFitUtils::CoarseResolutionType::CZT:         return "CZT";
    case PeakFitUtils::CoarseResolutionType::MedRes:      return "MedRes";
    case PeakFitUtils::CoarseResolutionType::LowOrMedRes: return "LowOrMedRes";
    case PeakFitUtils::CoarseResolutionType::High:        return "High";
    case PeakFitUtils::CoarseResolutionType::Unknown:     return "Unknown";
  }
  return "Unknown";
}


PeakDef::SkewType skew_type_for_detector( const string &name )
{
  for( const char *s : { "1.5cm-2cm", "1cm-1cm", "CZT_H3D", "Kromek-GR1", "nanoRaider", "Interceptor", "Raider" } )
    if( SpecUtils::icontains( name, s ) )
      return PeakDef::SkewType::DoubleSidedCrystalBall;
  return PeakDef::SkewType::NoSkew;
}


struct InjectOptions
{
  string base_dir = "/Users/wcjohns/coding/InterSpec_peak_fit_improve/peak_fit_accuracy_inject_compact";
  vector<string> dets{ "Detective-X", "Falcon 5000", "HPGe_Planar_50%", "IdentiFINDER-R500-NaI",
                       "SAM-Eagle-NaI-3x3", "IdentiFINDER-LaBr3", "Kromek-GR1-CZT",
                       "CZT_H3D_M400_ORNL_25cm", "Radiacode-102", "RadEagle" };
  string location = "Livermore";
  vector<string> times{ "30", "300", "1800" };
  size_t max_sources = 50;
  size_t max_roi_peaks = 4;
  PeakContinuum::OffsetType cont = PeakContinuum::OffsetType::Linear;
  string source_filter;      // substring; empty = all
  std::optional<PeakDef::SkewType> skew;  // overrides the per-detector default

  // ROI extent below the lowest / above the highest peak, in FWHM; HPGe and other detectors.  Wide
  //  low-resolution ROIs take in structure (e.g., backscatter peaks) a polynomial cannot follow.
  double hpge_roi_lower_fwhm = 4.0, hpge_roi_upper_fwhm = 3.5;
  double lowres_roi_lower_fwhm = 2.5, lowres_roi_upper_fwhm = 2.5;
};


vector<Problem> build_inject_problems( const InjectOptions &opts )
{
  vector<Problem> problems;
  size_t num_skipped_groups = 0, num_files = 0;

  for( const string &det : opts.dets )
  {
    const PeakFitUtils::CoarseResolutionType det_type = true_det_type_for_name( det );
    const PeakDef::SkewType skew = opts.skew.has_value() ? *opts.skew : skew_type_for_detector( det );
    const bool is_hpge = (det_type == PeakFitUtils::CoarseResolutionType::High);
    const double min_energy = is_hpge ? 30.0 : 50.0;
    const double roi_lower_fwhm = is_hpge ? opts.hpge_roi_lower_fwhm : opts.lowres_roi_lower_fwhm;
    const double roi_upper_fwhm = is_hpge ? opts.hpge_roi_upper_fwhm : opts.lowres_roi_upper_fwhm;

    for( const string &time : opts.times )
    {
      const string dir = SpecUtils::append_path( SpecUtils::append_path( SpecUtils::append_path(
                                                    opts.base_dir, det ), opts.location ), time + "_seconds" );
      if( !SpecUtils::is_directory( dir ) )
      {
        cerr << "Warning: no directory '" << dir << "'" << endl;
        continue;
      }

      vector<string> truth_files = SpecUtils::ls_files_in_directory( dir, "_truth.csv" );
      std::sort( begin(truth_files), end(truth_files) );

      // Parse first, so the sub-sampling is over sources that have any truth at all.
      vector<pair<string,FitPeaksCorpus::InjectTruth>> candidates;
      for( const string &tf : truth_files )
      {
        if( !opts.source_filter.empty() && !SpecUtils::icontains( SpecUtils::filename(tf), opts.source_filter ) )
          continue;
        try
        {
          FitPeaksCorpus::InjectTruth truth = FitPeaksCorpus::parse_inject_truth_csv( tf );
          if( truth.merged_signal.empty() && truth.merged_background.empty() )
            continue;
          candidates.emplace_back( tf, std::move(truth) );
        }catch( std::exception &e )
        {
          cerr << "Warning: skipping '" << tf << "': " << e.what() << endl;
        }
      }//for( const string &tf : truth_files )

      vector<size_t> picks;
      if( candidates.size() <= opts.max_sources )
      {
        for( size_t i = 0; i < candidates.size(); ++i )
          picks.push_back( i );
      }else
      {
        for( size_t i = 0; i < opts.max_sources; ++i )
          picks.push_back( (i * candidates.size()) / opts.max_sources );
      }

      for( const size_t pick : picks )
      {
        const string &truth_path = candidates[pick].first;
        const FitPeaksCorpus::InjectTruth &truth = candidates[pick].second;
        const string src_name = SpecUtils::filename( truth_path ).substr( 0, SpecUtils::filename( truth_path ).size() - 10 );
        const string pcf_path = SpecUtils::append_path( dir, src_name + ".pcf" );

        SpecUtils::SpecFile spec;
        if( !spec.load_file( pcf_path, SpecUtils::ParserType::Pcf, "pcf" ) || (spec.num_measurements() != 4) )
        {
          cerr << "Warning: could not load 4-record PCF '" << pcf_path << "'" << endl;
          continue;
        }
        ++num_files;

        const shared_ptr<const SpecUtils::Measurement> expected = spec.measurements()[3];
        const float spec_lower = expected->gamma_energy_min();
        const float spec_upper = expected->gamma_energy_max();

        // Signal and background together, re-merged where they coincide.
        vector<FitPeaksCorpus::TruthPhotopeak> all_rows = truth.merged_signal;
        all_rows.insert( end(all_rows), begin(truth.merged_background), end(truth.merged_background) );
        std::sort( begin(all_rows), end(all_rows), []( const FitPeaksCorpus::TruthPhotopeak &a,
                                                      const FitPeaksCorpus::TruthPhotopeak &b ){
          return a.energy < b.energy;
        } );
        all_rows = FitPeaksCorpus::merge_unresolved_photopeaks( all_rows );

        vector<FitPeaksCorpus::TruthPhotopeak> rows;
        for( const FitPeaksCorpus::TruthPhotopeak &r : all_rows )
        {
          if( (r.energy >= min_energy) && (r.area > 0.0) && (r.fwhm > 0.0)
              && ((r.energy - 5.0*r.fwhm) > spec_lower) && ((r.energy + 5.0*r.fwhm) < spec_upper) )
            rows.push_back( r );
        }

        // Group into ROIs where the [E - roi_lower_fwhm*FWHM, E + roi_upper_fwhm*FWHM] windows overlap.
        vector<vector<FitPeaksCorpus::TruthPhotopeak>> groups;
        for( const FitPeaksCorpus::TruthPhotopeak &r : rows )
        {
          if( !groups.empty() )
          {
            const FitPeaksCorpus::TruthPhotopeak &prev = groups.back().back();
            if( (prev.energy + roi_upper_fwhm*prev.fwhm) > (r.energy - roi_lower_fwhm*r.fwhm) )
            {
              groups.back().push_back( r );
              continue;
            }
          }
          groups.push_back( { r } );
        }//for( rows )

        for( size_t gi = 0; gi < groups.size(); ++gi )
        {
          const vector<FitPeaksCorpus::TruthPhotopeak> &grp = groups[gi];
          if( grp.size() > opts.max_roi_peaks )
          {
            ++num_skipped_groups;
            continue;
          }

          Problem prob;
          prob.det_type = det_type;
          prob.skew_type = skew;
          prob.fit_skew = (skew != PeakDef::SkewType::NoSkew);
          prob.fit_cont_type = opts.cont;
          prob.expected = expected;
          prob.spectrum_key = det + "/" + opts.location + "/" + time + "/" + src_name;
          prob.roi_lower = grp.front().energy - roi_lower_fwhm*grp.front().fwhm;
          prob.roi_upper = grp.back().energy + roi_upper_fwhm*grp.back().fwhm;

          const size_t nskew = PeakDef::num_skew_parameters( skew );
          for( size_t i = 0; i < nskew; ++i )
          {
            double lower = 0.0, upper = 0.0, start = 0.0, step = 0.0;
            PeakDef::skew_parameter_range( skew, static_cast<PeakDef::CoefficientType>( PeakDef::SkewPar0 + i ),
                                           lower, upper, start, step );
            prob.skew_pars.push_back( start );
          }

          for( const FitPeaksCorpus::TruthPhotopeak &r : grp )
          {
            TruthPeak t;
            t.mean = r.energy;
            t.fwhm = r.fwhm;
            t.area = r.area;
            t.z = r.z_det();
            t.signal = r.signal;
            prob.truth.push_back( t );
          }

          char ebuf[32];
          snprintf( ebuf, sizeof(ebuf), "%.1f", grp.front().energy );
          prob.id = "inj_" + det + "_" + time + "_" + src_name + "_" + ebuf;
          SpecUtils::ireplace_all( prob.id, " ", "-" );
          prob.dims = {
            { "det", det }, { "det_type", det_type_str( det_type ) }, { "time", time },
            { "source", src_name }, { "roi_npeaks", std::to_string( grp.size() ) }
          };
          problems.push_back( std::move(prob) );
        }//for( groups )
      }//for( picks )
    }//for( times )
  }//for( dets )

  cerr << "Inject: " << num_files << " spectra, " << problems.size() << " ROIs ("
       << num_skipped_groups << " groups skipped for having more than " << opts.max_roi_peaks << " peaks)." << endl;

  return problems;
}//build_inject_problems(...)


// ------------------------------------------------------------------------------------------------
// Fitting and per-fit records
// ------------------------------------------------------------------------------------------------
struct FitRecord
{
  size_t problem = 0, replica = 0, objective = 0, peak = 0;
  string status;   // ok, throw, empty, lost
  double area = NAN, area_unc = NAN, mean = NAN, mean_unc = NAN, fwhm = NAN, fwhm_unc = NAN;
  double chi2dof_stamp = NAN, neyman_chi2 = NAN, deviance = NAN;
  size_t nchan = 0;
  double cpu_s = NAN;
  string message;

  // Candidate "is this ROI sparse" statistics, from this fit's model m over the ROI channels:
  //  min m, mean m, sqrt(sum 1/max(m,1)) over the ROI, the same over the peak region (+-1.5 FWHM of
  //  the fit peaks), and N/sqrt(T) over the peak region (channels / sqrt(model counts)).
  double sp_min = NAN, sp_mean = NAN, sp_roi = NAN, sp_peak = NAN, sp_nt = NAN;

  // From PeakFitLM::take_fit_objective_diagnostics()
  bool fell_back = false;
  size_t irls_passes = 0;
  bool irls_converged = false;
  size_t sparse_rois = 0;
};


struct RunOptions
{
  string mode = "synthetic";
  string out_dir;
  string datadir;
  size_t replicas = 20;
  int threads = 4;
  uint64_t seed = 1;
  string path = "core";        // core: fit_peaks_in_roi_LM; refit: peak-search chi2 fit, then refitPeaksThatShareROI_LM
  double scale = 1.0;          // inject: multiply expected counts
  vector<string> objectives;   // empty = all available
  set<string> suites{ "primary", "cont", "skew", "topo" };
  string only;                 // substring filter on problem ids
  long only_replica = -1;
  bool no_pseudo = false;
  string dump;                 // "SUBSTR:REPLICA": write channel data and each objective's model
  InjectOptions inject;
};


vector<shared_ptr<const PeakDef>> sorted_by_mean( vector<shared_ptr<const PeakDef>> peaks )
{
  std::sort( begin(peaks), end(peaks), []( const shared_ptr<const PeakDef> &a, const shared_ptr<const PeakDef> &b ){
    return a->mean() < b->mean();
  } );
  return peaks;
}


/** Each objective's fit of the noise-free spectrum, scaled by `sm_pseudo_scale`, from unperturbed
 starting values.
 */
void compute_pseudo_truth( Problem &prob, const double scale, const vector<Objective> &objectives )
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for( TruthPeak &t : prob.truth )
  {
    t.pseudo_area.assign( objectives.size(), nan );
    t.pseudo_mean.assign( objectives.size(), nan );
    t.pseudo_fwhm.assign( objectives.size(), nan );
  }

  const shared_ptr<const SpecUtils::Measurement> data = scaled_spectrum( prob.expected, sm_pseudo_scale*scale );
  for( size_t oi = 0; oi < objectives.size(); ++oi )
  {
    try
    {
      mt19937 rng( 0 );
      const vector<shared_ptr<const PeakDef>> start = make_start_peaks( prob, sm_pseudo_scale*scale, false, rng );
      const vector<shared_ptr<const PeakDef>> fit
           = sorted_by_mean( PeakFitLM::fit_peaks_in_roi_LM( start, data, prob.det_type, objectives[oi].flags ) );
      if( fit.size() != prob.truth.size() )
        continue;
      for( size_t i = 0; i < fit.size(); ++i )
      {
        prob.truth[i].pseudo_area[oi] = fit[i]->amplitude() / sm_pseudo_scale;
        prob.truth[i].pseudo_mean[oi] = fit[i]->mean();
        prob.truth[i].pseudo_fwhm[oi] = fit[i]->fwhm();
      }
    }catch( std::exception & )
    {
    }
  }//for( objectives )
}//compute_pseudo_truth(...)


string fmt( const double v );

vector<FitRecord> run_one( const vector<Problem> &problems, const size_t prob_index, const size_t replica,
                           const vector<Objective> &objectives, const RunOptions &opt )
{
  const Problem &prob = problems[prob_index];
  const double scale = (opt.mode == "inject") ? opt.scale : 1.0;

  mt19937 data_rng( seed_for( prob.spectrum_key, replica, "data", opt.seed ) );
  const shared_ptr<const SpecUtils::Measurement> data = poisson_sample( prob.expected, scale, data_rng );

  const size_t ch0 = data->find_gamma_channel( static_cast<float>(prob.roi_lower) );
  const size_t ch1 = data->find_gamma_channel( static_cast<float>(prob.roi_upper) );
  const vector<float> &counts = *data->gamma_counts();

  bool dump = false;
  if( !opt.dump.empty() )
  {
    const size_t colon = opt.dump.rfind( ':' );
    const string substr = opt.dump.substr( 0, colon );
    const long dump_replica = (colon == string::npos) ? 0 : std::stol( opt.dump.substr( colon + 1 ) );
    dump = (prob.id.find( substr ) != string::npos) && (static_cast<long>(replica) == dump_replica);
  }
  vector<vector<double>> dump_models( objectives.size() );

  vector<FitRecord> records;
  for( size_t oi = 0; oi < objectives.size(); ++oi )
  {
    // Identical starting values for every objective.
    mt19937 start_rng( seed_for( prob.id, replica, "start", opt.seed ) );
    const vector<shared_ptr<const PeakDef>> start = make_start_peaks( prob, scale, true, start_rng );

    vector<FitRecord> these( prob.truth.size() );
    for( size_t pi = 0; pi < these.size(); ++pi )
    {
      these[pi].problem = prob_index;
      these[pi].replica = replica;
      these[pi].objective = oi;
      these[pi].peak = pi;
      these[pi].nchan = ch1 - ch0 + 1;
    }

    vector<shared_ptr<const PeakDef>> fit;
    string error;
    double cpu = NAN;
    PeakFitLM::take_fit_objective_diagnostics();  //reset
    try
    {
      if( opt.path == "refit" )
      {
        // The peaks a refit typically starts from: the peak search's chi2 fit.
        const Wt::WFlags<PeakFitLM::PeakFitLMOptions> search_options( PeakFitLM::PeakFitLMOptions::NoSparseDataLikelihood );
        const vector<shared_ptr<const PeakDef>> initial = PeakFitLM::fit_peaks_in_roi_LM( start, data, prob.det_type,
                                                                                          search_options );
        PeakFitLM::take_fit_objective_diagnostics();  //only the refit's
        const double t0 = thread_cpu_seconds();
        fit = PeakFitLM::refitPeaksThatShareROI_LM( data, nullptr, initial, prob.det_type, objectives[oi].flags );
        cpu = thread_cpu_seconds() - t0;
      }else
      {
        const double t0 = thread_cpu_seconds();
        fit = PeakFitLM::fit_peaks_in_roi_LM( start, data, prob.det_type, objectives[oi].flags );
        cpu = thread_cpu_seconds() - t0;
      }
    }catch( std::exception &e )
    {
      error = e.what();
    }
    const PeakFitLM::FitObjectiveDiagnostics diag = PeakFitLM::take_fit_objective_diagnostics();

    fit = sorted_by_mean( fit );
    if( objectives[oi].conditional_uncert )
    {
      vector<shared_ptr<const PeakDef>> conditional;
      for( const shared_ptr<const PeakDef> &p : fit )
      {
        auto copy = make_shared<PeakDef>( *p );
        const double z = PeakFitLM::peak_detection_significance( *p, fit, data, /*chi2_weights=*/ true );
        if( (z != 0.0) && std::isfinite( z ) )
          copy->setAmplitudeUncert( std::fabs( p->amplitude() / z ) );
        conditional.push_back( copy );
      }
      fit = conditional;
    }

    double neyman = NAN, deviance = NAN;
    double sp_min = NAN, sp_mean = NAN, sp_roi = NAN, sp_peak = NAN, sp_nt = NAN;
    if( !fit.empty() )
    {
      try
      {
        const vector<double> model = model_counts( fit, data, ch0, ch1 );
        {
          double mn = std::numeric_limits<double>::infinity(), sum = 0.0, inv_roi = 0.0, inv_pk = 0.0, t_pk = 0.0;
          size_t n_pk = 0;
          for( size_t ch = ch0; ch <= ch1; ++ch )
          {
            const double m = model[ch - ch0];
            mn = std::min( mn, m );
            sum += m;
            inv_roi += 1.0 / std::max( m, 1.0 );
            const double e = 0.5*(data->gamma_channel_lower( ch ) + data->gamma_channel_upper( ch ));
            bool in_peak = false;
            for( const shared_ptr<const PeakDef> &p : fit )
              in_peak |= (fabs( e - p->mean() ) <= 1.5*p->fwhm());
            if( in_peak )
            {
              inv_pk += 1.0 / std::max( m, 1.0 );
              t_pk += std::max( m, 0.0 );
              ++n_pk;
            }
          }
          sp_min = mn;
          sp_mean = sum / static_cast<double>( ch1 - ch0 + 1 );
          sp_roi = std::sqrt( inv_roi );
          sp_peak = std::sqrt( inv_pk );
          sp_nt = (t_pk > 0.0) ? (static_cast<double>( n_pk ) / std::sqrt( t_pk )) : NAN;
        }
        if( dump )
          dump_models[oi] = model;
        neyman = deviance = 0.0;
        for( size_t ch = ch0; ch <= ch1; ++ch )
        {
          const double n = counts[ch], m = model[ch - ch0];
          neyman += (n - m)*(n - m) / std::max( n, 1.0 );
          deviance += poisson_deviance_term( n, m );
        }
      }catch( std::exception & )
      {
      }
    }//if( !fit.empty() )

    // Match fitted peaks to truth peaks: by order when the counts agree, else nearest mean within 1.5 sigma.
    vector<shared_ptr<const PeakDef>> matched( prob.truth.size() );
    if( fit.size() == prob.truth.size() )
    {
      matched = fit;
    }else
    {
      for( size_t pi = 0; pi < prob.truth.size(); ++pi )
      {
        double best = 1.5 * prob.truth[pi].fwhm / 2.35482;
        for( const shared_ptr<const PeakDef> &p : fit )
        {
          const double d = fabs( p->mean() - prob.truth[pi].mean );
          if( d < best )
          {
            best = d;
            matched[pi] = p;
          }
        }
      }
    }//if( same number of peaks ) / else

    for( size_t pi = 0; pi < these.size(); ++pi )
    {
      FitRecord &r = these[pi];
      r.cpu_s = cpu;
      r.sp_min = sp_min;
      r.sp_mean = sp_mean;
      r.sp_roi = sp_roi;
      r.sp_peak = sp_peak;
      r.sp_nt = sp_nt;
      r.fell_back = diag.fell_back_to_chi2;
      r.irls_passes = diag.irls_passes;
      r.irls_converged = diag.irls_converged;
      r.sparse_rois = diag.sparse_rois;
      r.neyman_chi2 = neyman;
      r.deviance = deviance;
      const shared_ptr<const PeakDef> &p = matched[pi];
      if( !error.empty() )
      {
        r.status = "throw";
        r.message = error;
      }else if( fit.empty() )
      {
        r.status = "empty";
      }else if( !p )
      {
        r.status = "lost";
      }else
      {
        r.status = "ok";
        r.area = p->amplitude();
        r.area_unc = p->amplitudeUncert();
        r.mean = p->mean();
        r.mean_unc = p->meanUncert();
        r.fwhm = p->fwhm();
        r.fwhm_unc = 2.35482 * p->sigmaUncert();
        r.chi2dof_stamp = p->chi2dof();
      }
    }//for( these )

    records.insert( end(records), begin(these), end(these) );
  }//for( objectives )

  if( dump )
  {
    string name = "dump_" + prob.id + "_r" + std::to_string( replica ) + ".tsv";
    SpecUtils::ireplace_all( name, "/", "_" );
    ofstream out( SpecUtils::append_path( opt.out_dir, name ) );
    const vector<float> &expected = *prob.expected->gamma_counts();
    out << "energy_lower\tenergy_upper\tdata\texpected";
    for( const Objective &o : objectives )
      out << "\tmodel_" << o.name;
    out << "\n";
    for( size_t ch = ch0; ch <= ch1; ++ch )
    {
      out << data->gamma_channel_lower( ch ) << "\t" << data->gamma_channel_upper( ch )
          << "\t" << counts[ch] << "\t" << fmt( scale*expected[ch] );
      for( const vector<double> &m : dump_models )
        out << "\t" << (m.empty() ? string("nan") : fmt( m[ch - ch0] ));
      out << "\n";
    }
  }//if( dump )

  return records;
}//run_one(...)


// ------------------------------------------------------------------------------------------------
// Statistics
// ------------------------------------------------------------------------------------------------
double median_of( vector<double> v )
{
  if( v.empty() )
    return NAN;
  std::sort( begin(v), end(v) );
  const size_t n = v.size();
  return (n % 2) ? v[n/2] : 0.5*(v[n/2 - 1] + v[n/2]);
}


double quantile_of( vector<double> v, const double q )
{
  if( v.empty() )
    return NAN;
  std::sort( begin(v), end(v) );
  const double pos = q * static_cast<double>( v.size() - 1 );
  const size_t lo = static_cast<size_t>( std::floor( pos ) );
  const size_t hi = std::min( lo + 1, v.size() - 1 );
  return v[lo] + (pos - static_cast<double>(lo))*(v[hi] - v[lo]);
}


struct Moments
{
  double n = 0.0, mean = NAN, sd = NAN;
};


Moments moments_of( const vector<double> &v )
{
  Moments m;
  m.n = static_cast<double>( v.size() );
  if( v.empty() )
    return m;
  double sum = 0.0;
  for( const double x : v )
    sum += x;
  m.mean = sum / m.n;
  if( v.size() > 1 )
  {
    double ss = 0.0;
    for( const double x : v )
      ss += (x - m.mean)*(x - m.mean);
    m.sd = std::sqrt( ss / (m.n - 1.0) );
  }
  return m;
}


string fmt( const double v )
{
  if( std::isnan( v ) )
    return "nan";
  ostringstream strm;
  strm << std::setprecision( 7 ) << v;
  return strm.str();
}


void write_outputs( const string &out_dir, const vector<Problem> &problems, const vector<Objective> &objectives,
                    const vector<FitRecord> &records, const RunOptions &opt )
{
  const double scale = (opt.mode == "inject") ? opt.scale : 1.0;

  // --- per_fit.tsv ---
  {
    ofstream out( SpecUtils::append_path( out_dir, "per_fit.tsv" ) );
    out << "problem\treplica\tobjective\tpeak\tstatus\ttruth_mean\ttruth_fwhm\ttruth_area\tpseudo_area"
           "\tarea\tarea_unc\tmean\tmean_unc\tfwhm\tfwhm_unc\tchi2dof_stamp\tneyman_chi2\tdeviance\tnchan\tcpu_s\tmessage"
           "\tfell_back\tirls_passes\tirls_converged\tsp_min\tsp_mean\tsp_roi\tsp_peak\tsp_nt\tsparse_rois\n";
    for( const FitRecord &r : records )
    {
      const Problem &prob = problems[r.problem];
      const TruthPeak &t = prob.truth[r.peak];
      string msg = r.message;
      SpecUtils::ireplace_all( msg, "\t", " " );
      SpecUtils::ireplace_all( msg, "\n", " " );
      out << prob.id << "\t" << r.replica << "\t" << objectives[r.objective].name << "\t" << r.peak
          << "\t" << r.status << "\t" << fmt(t.mean) << "\t" << fmt(t.fwhm) << "\t" << fmt(scale*t.area)
          << "\t" << fmt(scale*t.pseudo_area_for(r.objective))
          << "\t" << fmt(r.area) << "\t" << fmt(r.area_unc) << "\t" << fmt(r.mean) << "\t" << fmt(r.mean_unc)
          << "\t" << fmt(r.fwhm) << "\t" << fmt(r.fwhm_unc) << "\t" << fmt(r.chi2dof_stamp)
          << "\t" << fmt(r.neyman_chi2) << "\t" << fmt(r.deviance) << "\t" << r.nchan
          << "\t" << fmt(r.cpu_s) << "\t" << msg
          << "\t" << (r.fell_back ? 1 : 0) << "\t" << r.irls_passes << "\t" << (r.irls_converged ? 1 : 0)
          << "\t" << fmt(r.sp_min) << "\t" << fmt(r.sp_mean) << "\t" << fmt(r.sp_roi) << "\t" << fmt(r.sp_peak)
          << "\t" << fmt(r.sp_nt) << "\t" << r.sparse_rois << "\n";
    }
  }

  // Group records by (problem, peak, objective).
  map<array<size_t,3>,vector<const FitRecord *>> groups;
  for( const FitRecord &r : records )
    groups[{ r.problem, r.peak, r.objective }].push_back( &r );

  size_t chi2_index = objectives.size();
  for( size_t i = 0; i < objectives.size(); ++i )
    if( objectives[i].name == "chi2" )
      chi2_index = i;

  // --- summary.tsv ---
  ofstream summary( SpecUtils::append_path( out_dir, "summary.tsv" ) );
  summary << "problem";
  if( !problems.empty() )
    for( const pair<string,string> &d : problems.front().dims )
      summary << "\t" << d.first;
  summary << "\tobjective\tpeak\tsignal\ttruth_area\tpseudo_area\tpseudo_vs_chi2\ttruth_fwhm\tz\tn\tn_ok\tfail_rate"
             "\tarea_mean\tarea_sd\trel_bias_truth\trel_bias_truth_se\trel_bias_pseudo"
             "\tpull_mean\tpull_rms\tcov1\tcov2\tunc_ratio\ttail_rate"
             "\tmean_bias_sigma\tmean_pull_mean\tmean_pull_rms"
             "\tfwhm_rel_bias\tfwhm_pull_mean\tfwhm_pull_rms"
             "\tcpu_median\tcpu_p95\tcpu_mean\tirls_passes_mean\tirls_unconverged\tfell_back\n";

  // Paired comparison against the chi2 objective (index 0 when present).
  ofstream paired( SpecUtils::append_path( out_dir, "paired.tsv" ) );
  paired << "problem\treplica\tpeak\tobjective\tchi2_area\tobj_area\tchi2_sd\tdiff_in_sd\tchi2_status\tobj_status\n";

  for( const auto &g : groups )
  {
    const Problem &prob = problems[g.first[0]];
    const TruthPeak &t = prob.truth[g.first[1]];
    const Objective &obj = objectives[g.first[2]];
    const vector<const FitRecord *> &recs = g.second;

    const double truth_area = scale * t.area;
    const double pseudo_area = scale * t.pseudo_area_for( g.first[2] );
    // The statistical reference: this objective's pseudo-truth when we have it, else the exact truth.
    const double ref_area = std::isnan( pseudo_area ) ? truth_area : pseudo_area;
    const double ref_mean = std::isnan( t.pseudo_mean_for( g.first[2] ) ) ? t.mean : t.pseudo_mean_for( g.first[2] );
    const double ref_fwhm = std::isnan( t.pseudo_fwhm_for( g.first[2] ) ) ? t.fwhm : t.pseudo_fwhm_for( g.first[2] );
    const double truth_sigma = t.fwhm / 2.35482;

    vector<double> areas, pulls, mean_pulls, mean_offsets, fwhms, fwhm_pulls, cpus, unc, irls_passes;
    size_t n_ok = 0, n_cov1 = 0, n_cov2 = 0, n_unconverged = 0, n_fell_back = 0;
    for( const FitRecord *r : recs )
    {
      if( !std::isnan( r->cpu_s ) )
        cpus.push_back( r->cpu_s );
      irls_passes.push_back( static_cast<double>( r->irls_passes ) );
      n_unconverged += ((r->irls_passes > 0) && !r->irls_converged) ? 1 : 0;
      n_fell_back += r->fell_back ? 1 : 0;
      if( r->status != "ok" )
        continue;
      ++n_ok;
      areas.push_back( r->area );
      unc.push_back( r->area_unc );
      if( r->area_unc > 0.0 )
      {
        const double pull = (r->area - ref_area) / r->area_unc;
        pulls.push_back( pull );
        n_cov1 += (fabs(pull) < 1.0);
        n_cov2 += (fabs(pull) < 2.0);
      }
      mean_offsets.push_back( (r->mean - ref_mean) / truth_sigma );
      if( r->mean_unc > 0.0 )
        mean_pulls.push_back( (r->mean - ref_mean) / r->mean_unc );
      fwhms.push_back( r->fwhm );
      if( r->fwhm_unc > 0.0 )
        fwhm_pulls.push_back( (r->fwhm - ref_fwhm) / r->fwhm_unc );
    }//for( recs )

    const Moments am = moments_of( areas );
    const Moments pm = moments_of( pulls );
    const Moments mo = moments_of( mean_offsets );
    const Moments mpm = moments_of( mean_pulls );
    const Moments fm = moments_of( fwhms );
    const Moments fpm = moments_of( fwhm_pulls );

    const auto rms_of = []( const vector<double> &v ) -> double {
      if( v.empty() )
        return NAN;
      double ss = 0.0;
      for( const double x : v )
        ss += x*x;
      return std::sqrt( ss / static_cast<double>( v.size() ) );
    };

    // Tail: fits further than 5 robust-sigma from the median area.
    double tail_rate = NAN;
    if( areas.size() >= 5 )
    {
      const double med = median_of( areas );
      vector<double> absdev;
      for( const double a : areas )
        absdev.push_back( fabs( a - med ) );
      const double mad_sigma = std::max( 1.4826 * median_of( absdev ), 1.0E-9 );
      size_t ntail = 0;
      for( const double a : areas )
        ntail += (fabs( a - med ) > 5.0*mad_sigma);
      tail_rate = static_cast<double>( ntail ) / static_cast<double>( areas.size() );
    }

    const double n = static_cast<double>( recs.size() );
    summary << prob.id;
    for( const pair<string,string> &d : prob.dims )
      summary << "\t" << d.second;
    summary << "\t" << obj.name << "\t" << g.first[1] << "\t" << (t.signal ? 1 : 0)
            << "\t" << fmt(truth_area) << "\t" << fmt(pseudo_area)
            << "\t" << fmt( pseudo_area / (scale * t.pseudo_area_for( chi2_index )) - 1.0 )
            << "\t" << fmt(t.fwhm) << "\t" << fmt(t.z)
            << "\t" << recs.size() << "\t" << n_ok << "\t" << fmt( 1.0 - n_ok/n )
            << "\t" << fmt(am.mean) << "\t" << fmt(am.sd)
            << "\t" << fmt(am.mean/truth_area - 1.0) << "\t" << fmt(am.sd/std::sqrt(std::max(am.n,1.0))/truth_area)
            << "\t" << fmt(am.mean/pseudo_area - 1.0)
            << "\t" << fmt(pm.mean) << "\t" << fmt(rms_of(pulls))
            << "\t" << fmt( pulls.empty() ? NAN : (double(n_cov1)/pulls.size()) )
            << "\t" << fmt( pulls.empty() ? NAN : (double(n_cov2)/pulls.size()) )
            << "\t" << fmt( median_of(unc) / am.sd )
            << "\t" << fmt(tail_rate)
            << "\t" << fmt(mo.mean) << "\t" << fmt(mpm.mean) << "\t" << fmt(rms_of(mean_pulls))
            << "\t" << fmt(fm.mean/ref_fwhm - 1.0) << "\t" << fmt(fpm.mean) << "\t" << fmt(rms_of(fwhm_pulls))
            << "\t" << fmt(median_of(cpus)) << "\t" << fmt(quantile_of(cpus, 0.95)) << "\t" << fmt(moments_of(cpus).mean)
            << "\t" << fmt(moments_of(irls_passes).mean) << "\t" << fmt( n_unconverged/n ) << "\t" << fmt( n_fell_back/n )
            << "\n";

    // Paired rows: replicas where this objective differs from chi2 by more than 3 chi2-sd.
    if( (chi2_index < objectives.size()) && (g.first[2] != chi2_index) )
    {
      const auto chi2_iter = groups.find( { g.first[0], g.first[1], chi2_index } );
      if( chi2_iter != end(groups) )
      {
        vector<double> chi2_areas;
        map<size_t,const FitRecord *> chi2_by_replica;
        for( const FitRecord *r : chi2_iter->second )
        {
          chi2_by_replica[r->replica] = r;
          if( r->status == "ok" )
            chi2_areas.push_back( r->area );
        }
        const double chi2_sd = moments_of( chi2_areas ).sd;

        for( const FitRecord *r : recs )
        {
          const auto c = chi2_by_replica.find( r->replica );
          if( c == end(chi2_by_replica) )
            continue;
          const FitRecord *cr = c->second;
          const bool both_ok = (r->status == "ok") && (cr->status == "ok");
          const double diff = (both_ok && (chi2_sd > 0.0)) ? ((r->area - cr->area) / chi2_sd) : NAN;
          if( !both_ok || (fabs(diff) > 3.0) )
          {
            paired << prob.id << "\t" << r->replica << "\t" << r->peak << "\t" << obj.name
                   << "\t" << fmt(cr->area) << "\t" << fmt(r->area) << "\t" << fmt(chi2_sd) << "\t" << fmt(diff)
                   << "\t" << cr->status << "\t" << r->status << "\n";
          }
        }//for( recs )
      }//if( have chi2 group )
    }//if( compare to chi2 )
  }//for( groups )
}//write_outputs(...)


// ------------------------------------------------------------------------------------------------
// Command line
// ------------------------------------------------------------------------------------------------
void print_usage()
{
  cout << "peak_fit_objective_eval - compare PeakFitLM fit objectives on Poisson replicas with known truth\n"
          "\n"
          "  --mode=synthetic|inject       problem source (default synthetic)\n"
          "  --out=DIR                     output directory (required)\n"
          "  --replicas=N                  Poisson replicas per problem (default 20)\n"
          "  --threads=N                   worker threads (default 4)\n"
          "  --seed=N                      run seed (default 1)\n"
          "  --objectives=a,b,...          objectives to run (default all available):";
  for( const Objective &o : available_objectives() )
    cout << " " << o.name;
  cout << "\n"
          "  --path=core|refit             core: fit_peaks_in_roi_LM from perturbed truth; refit: a chi2 fit as\n"
          "                                the peak search makes, then refitPeaksThatShareROI_LM with the\n"
          "                                objective (default core)\n"
          "  --only=SUBSTR                 only problems whose id contains SUBSTR\n"
          "  --replica=K                   only replica K\n"
          "  --no-pseudo                   skip the noise-free pseudo-truth fits\n"
          "  --datadir=DIR                 InterSpec data directory (default: searched for)\n"
          " synthetic:\n"
          "  --suites=primary,cont,skew,topo  which parts of the grid (default all)\n"
          " inject:\n"
          "  --inject-dir=DIR              peak_fit_accuracy_inject_compact directory\n"
          "  --dets=A,B,...                detector directories (default: a 10-detector subset)\n"
          "  --location=NAME               (default Livermore)\n"
          "  --times=30,300,1800           live-time directories (default all three)\n"
          "  --max-sources=N               sources per directory, evenly sub-sampled (default 50)\n"
          "  --max-roi-peaks=N             skip ROIs with more truth peaks than this (default 4)\n"
          "  --source=SUBSTR               only sources containing SUBSTR\n"
          "  --cont=TYPE                   fit continuum type (default Linear)\n"
          "  --skew=TYPE                   skew type for every detector (default: DSCB for CZT, else NoSkew)\n"
          "  --hpge-roi-fwhm=LO,HI         HPGe ROI extent below/above the peaks, in FWHM (default 4,3.5)\n"
          "  --lowres-roi-fwhm=LO,HI       other detectors' ROI extent, in FWHM (default 2.5,2.5)\n"
          " output:\n"
          "  --dump=SUBSTR:K               for problems matching SUBSTR, replica K: write dump_<problem>.tsv\n"
          "                                with the channel data, expectation and each objective's model\n"
          "  --scale=K                     multiply the expected counts by K before sampling (default 1)\n"
          << endl;
}//print_usage()


vector<string> split_list( const string &s )
{
  vector<string> answer;
  SpecUtils::split( answer, s, "," );
  for( string &v : answer )
    SpecUtils::trim( v );
  answer.erase( std::remove( begin(answer), end(answer), string() ), end(answer) );
  return answer;
}


RunOptions parse_options( int argc, char **argv, bool &help )
{
  RunOptions opt;
  help = false;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( arg == "--help" || arg == "-h" )
    {
      help = true;
      continue;
    }

    const size_t eq = arg.find( '=' );
    if( (arg.substr( 0, 2 ) != "--") || (eq == string::npos) )
      throw runtime_error( "Unrecognized argument '" + arg + "'" );
    const string key = arg.substr( 2, eq - 2 );
    const string val = arg.substr( eq + 1 );

    if( key == "mode" ) opt.mode = val;
    else if( key == "out" ) opt.out_dir = val;
    else if( key == "datadir" ) opt.datadir = val;
    else if( key == "replicas" ) opt.replicas = static_cast<size_t>( std::stoul( val ) );
    else if( key == "threads" ) opt.threads = std::stoi( val );
    else if( key == "seed" ) opt.seed = std::stoull( val );
    else if( key == "path" ) opt.path = val;
    else if( key == "scale" ) opt.scale = std::stod( val );
    else if( key == "objectives" ) opt.objectives = split_list( val );
    else if( key == "suites" ) { const vector<string> s = split_list( val ); opt.suites = set<string>( begin(s), end(s) ); }
    else if( key == "only" ) opt.only = val;
    else if( key == "replica" ) opt.only_replica = std::stol( val );
    else if( key == "no-pseudo" ) opt.no_pseudo = (val != "0");
    else if( key == "inject-dir" ) opt.inject.base_dir = val;
    else if( key == "dets" ) opt.inject.dets = split_list( val );
    else if( key == "location" ) opt.inject.location = val;
    else if( key == "times" ) opt.inject.times = split_list( val );
    else if( key == "max-sources" ) opt.inject.max_sources = static_cast<size_t>( std::stoul( val ) );
    else if( key == "max-roi-peaks" ) opt.inject.max_roi_peaks = static_cast<size_t>( std::stoul( val ) );
    else if( key == "source" ) opt.inject.source_filter = val;
    else if( key == "skew" )
    {
      bool found = false;
      for( int t = 0; t < static_cast<int>( PeakDef::SkewType::NumSkewType ); ++t )
      {
        if( val == PeakDef::to_string( static_cast<PeakDef::SkewType>( t ) ) )
        {
          opt.inject.skew = static_cast<PeakDef::SkewType>( t );
          found = true;
        }
      }
      if( !found )
        throw runtime_error( "Unknown skew type '" + val + "'" );
    }
    else if( key == "dump" ) opt.dump = val;
    else if( key == "lowres-roi-fwhm" )
    {
      const vector<string> v = split_list( val );
      opt.inject.lowres_roi_lower_fwhm = std::stod( v.at(0) );
      opt.inject.lowres_roi_upper_fwhm = std::stod( v.at( v.size() > 1 ? 1 : 0 ) );
    }
    else if( key == "hpge-roi-fwhm" )
    {
      const vector<string> v = split_list( val );
      opt.inject.hpge_roi_lower_fwhm = std::stod( v.at(0) );
      opt.inject.hpge_roi_upper_fwhm = std::stod( v.at( v.size() > 1 ? 1 : 0 ) );
    }
    else if( key == "cont" )
    {
      opt.inject.cont = PeakContinuum::str_to_offset_type_str( val.c_str(), val.size() );
    }
    else throw runtime_error( "Unknown option '--" + key + "'" );
  }//for( args )

  if( opt.mode != "synthetic" && opt.mode != "inject" )
    throw runtime_error( "--mode must be synthetic or inject" );
  if( opt.path != "core" && opt.path != "refit" )
    throw runtime_error( "--path must be core or refit" );
  if( !help && opt.out_dir.empty() )
    throw runtime_error( "--out is required" );

  return opt;
}//parse_options(...)

}//namespace


int main( int argc, char **argv )
{
  RunOptions opt;
  bool help = false;
  try
  {
    opt = parse_options( argc, argv, help );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    print_usage();
    return 1;
  }

  if( help )
  {
    print_usage();
    return 0;
  }

  try
  {
    if( opt.datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data", "../../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
        {
          opt.datadir = d;
          break;
        }
      }
    }
    if( !opt.datadir.empty() )
      InterSpec::setStaticDataDirectory( opt.datadir );

    vector<Objective> objectives;
    const vector<Objective> available = available_objectives();
    if( opt.objectives.empty() )
    {
      objectives = available;
    }else
    {
      for( const string &name : opt.objectives )
      {
        const auto pos = std::find_if( begin(available), end(available), [&]( const Objective &o ){ return o.name == name; } );
        if( pos == end(available) )
          throw runtime_error( "Unknown objective '" + name + "'" );
        objectives.push_back( *pos );
      }
    }

    if( !SpecUtils::is_directory( opt.out_dir ) && !SpecUtils::create_directory( opt.out_dir ) )
      throw runtime_error( "Could not create '" + opt.out_dir + "'" );

    vector<Problem> problems = (opt.mode == "inject") ? build_inject_problems( opt.inject )
                                                      : build_synthetic_problems( opt.suites );
    if( !opt.only.empty() )
    {
      problems.erase( std::remove_if( begin(problems), end(problems), [&]( const Problem &p ){
        return p.id.find( opt.only ) == string::npos;
      } ), end(problems) );
    }
    if( problems.empty() )
      throw runtime_error( "No problems to run" );

    const auto wall_start = chrono::steady_clock::now();

    if( !opt.no_pseudo )
    {
      parallel_for( problems.size(), opt.threads, [&]( const size_t i ){
        compute_pseudo_truth( problems[i], (opt.mode == "inject") ? opt.scale : 1.0, objectives );
      } );
    }

    vector<pair<size_t,size_t>> work;
    for( size_t pi = 0; pi < problems.size(); ++pi )
      for( size_t r = 0; r < opt.replicas; ++r )
        if( (opt.only_replica < 0) || (static_cast<size_t>(opt.only_replica) == r) )
          work.emplace_back( pi, r );

    cerr << problems.size() << " problems, " << work.size() << " (problem, replica) units, "
         << objectives.size() << " objectives, " << opt.threads << " threads." << endl;

    vector<vector<FitRecord>> results( work.size() );
    std::atomic<size_t> ndone{ 0 };
    std::mutex progress_mutex;
    parallel_for( work.size(), opt.threads, [&]( const size_t i ){
      results[i] = run_one( problems, work[i].first, work[i].second, objectives, opt );
      const size_t done = ++ndone;
      if( (done % 500) == 0 )
      {
        std::lock_guard<std::mutex> lock( progress_mutex );
        cerr << "  " << done << " / " << work.size() << endl;
      }
    } );

    vector<FitRecord> records;
    for( vector<FitRecord> &r : results )
      records.insert( end(records), begin(r), end(r) );

    write_outputs( opt.out_dir, problems, objectives, records, opt );

    const double wall = chrono::duration<double>( chrono::steady_clock::now() - wall_start ).count();
    {
      ofstream meta( SpecUtils::append_path( opt.out_dir, "run_meta.txt" ) );
      meta << "command_line:";
      for( int i = 0; i < argc; ++i )
        meta << " " << argv[i];
      meta << "\nwall_seconds: " << wall << "\nproblems: " << problems.size()
           << "\nunits: " << work.size() << "\nobjectives:";
      for( const Objective &o : objectives )
        meta << " " << o.name;
      meta << "\n";
    }

    cerr << "Done in " << wall << " s; wrote " << opt.out_dir << endl;
  }catch( std::exception &e )
  {
    cerr << "Error: " << e.what() << endl;
    return 2;
  }

  return 0;
}//main(...)
