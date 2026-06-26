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
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#define BOOST_TEST_MODULE CrystalBallCztFit_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"

using namespace std;
using namespace boost::unit_test;

// Review item A26 / CZT validation.  CZT (e.g. H3D M400) photopeaks have a pronounced low-energy
//  tail from incomplete charge collection, so Crystal Ball (CB) and Double-Sided Crystal Ball (DSCB)
//  are the natural peak shapes - and are what InterSpec uses for CZT.  Fitting a real CZT Cs137
//  661.7 keV peak drives the rewritten (cancellation-free, no (n/alpha)^n) single-sided CB integral
//  AND the DSCB integral through the Ceres L-M (ceres::Jet) optimizer end-to-end - i.e. the
//  fit/derivative path that no other unit test exercises (CB/DSCB are never the default skew).

namespace
{
string g_data_dir;

void set_data_dir()
{
  static bool s_set = false;
  if( s_set )
    return;
  s_set = true;

  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      g_data_dir = arg.substr( 10 );
  }
  SpecUtils::ireplace_all( g_data_dir, "%20", " " );

  if( g_data_dir.empty()
     || !SpecUtils::is_file( SpecUtils::append_path(g_data_dir, "sandia.decay.xml") ) )
  {
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        g_data_dir = d;
        break;
      }
    }
  }

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(g_data_dir, "sandia.decay.xml") ),
                         "Could not locate the data directory (pass --datadir=...)" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
}//set_data_dir()


shared_ptr<const SpecUtils::Measurement> load_czt_cs137()
{
  set_data_dir();
  const string path = SpecUtils::append_path( g_data_dir,
        "reference_spectra/Common_Field_Nuclides/CZT_H3D_M400/Cs137_Unshielded.txt" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(path), "CZT Cs137 spectrum not at '" << path << "'" );

  SpecUtils::SpecFile specfile;
  BOOST_REQUIRE_MESSAGE( specfile.load_uri_file(path),
                         "Failed to decode RADDATA spectrum '" << path << "'" );

  const vector<shared_ptr<const SpecUtils::Measurement>> meas = specfile.measurements();
  BOOST_REQUIRE( !meas.empty() );
  const shared_ptr<const SpecUtils::Measurement> m = meas.front();
  BOOST_REQUIRE( m && (m->num_gamma_channels() > 16) );
  BOOST_REQUIRE( m->energy_calibration() && m->energy_calibration()->valid() );
  return m;
}//load_czt_cs137()


// Build a Cs137 661.7 keV peak with the given skew type over [roi_lower, roi_upper], fit it with the
//  Ceres L-M optimizer, and return the fitted peak (chi2dof populated).
shared_ptr<const PeakDef> fit_skew_peak( const shared_ptr<const SpecUtils::Measurement> &data,
                                         const PeakDef::SkewType skew_type,
                                         const double roi_lower, const double roi_upper,
                                         const PeakContinuum::OffsetType cont_type = PeakContinuum::Linear )
{
  const double mean0 = 661.657;
  const double sigma0 = 2.8;  // ~1% FWHM at 662 keV for the CZT H3D M400; the fit refines it.
  const double amp0 = (std::max)( 10.0, data->gamma_integral( static_cast<float>(roi_lower),
                                                            static_cast<float>(roi_upper) ) );

  const shared_ptr<PeakDef> peak = make_shared<PeakDef>( mean0, sigma0, amp0 );
  peak->setSkewType( skew_type );
  peak->setFitFor( PeakDef::Mean, true );
  peak->setFitFor( PeakDef::Sigma, true );
  peak->setFitFor( PeakDef::GaussAmplitude, true );

  // Seed each applicable skew coefficient at its starting value and mark it to be fit.
  for( int i = 0; i < 4; ++i )
  {
    const PeakDef::CoefficientType ct = static_cast<PeakDef::CoefficientType>( PeakDef::SkewPar0 + i );
    double lo = 0.0, up = 0.0, start = 0.0, step = 0.0;
    if( PeakDef::skew_parameter_range( skew_type, ct, lo, up, start, step ) )
    {
      peak->set_coefficient( start, ct );
      peak->setFitFor( ct, true );
    }
  }

  // A linear continuum over the narrow ROI reproduces the project owner's reported fit most closely
  //  (the DSCB/CB tails absorb the peak asymmetry); the type is a parameter so it can be varied.
  const shared_ptr<PeakContinuum> cont = peak->continuum();
  cont->setType( cont_type );
  cont->setRange( roi_lower, roi_upper );

  const bool isHPGe = PeakFitUtils::is_high_res( data );  // false for the low-channel CZT spectrum
  const vector<shared_ptr<const PeakDef>> input{ peak };
  const vector<shared_ptr<const PeakDef>> fit = PeakFitLM::fit_peaks_in_roi_LM( input, data, isHPGe );

  BOOST_REQUIRE_MESSAGE( fit.size() == 1, "Expected 1 fitted peak, got " << fit.size()
                         << " for skew " << PeakDef::to_string(skew_type) );
  return fit.front();
}//fit_skew_peak(...)


void check_skew_peak( const shared_ptr<const PeakDef> &peak, const PeakDef::SkewType skew_type,
                      const double roi_lower, const double roi_upper, const double chi2dof_ceiling,
                      const char * const label )
{
  BOOST_REQUIRE( peak );
  const double chi2dof = peak->chi2dof();
  const double mean = peak->mean();
  const double sigma = peak->sigma();
  const double amp = peak->amplitude();

  cout << "  " << label << ": chi2/dof=" << chi2dof << "  mean=" << mean << " keV  sigma=" << sigma
       << " keV  amp=" << amp << "  skew={";
  bool first = true;
  for( int i = 0; i < 4; ++i )
  {
    const PeakDef::CoefficientType ct = static_cast<PeakDef::CoefficientType>( PeakDef::SkewPar0 + i );
    double lo = 0.0, up = 0.0, start = 0.0, step = 0.0;
    if( PeakDef::skew_parameter_range( skew_type, ct, lo, up, start, step ) )
    {
      cout << (first ? "" : ", ") << peak->coefficient( ct );
      first = false;
    }
  }
  cout << "}" << endl;

  // The fit must converge to a sensible, finite, in-ROI peak with a good chi2/dof - a broken
  //  CB/DSCB integral or NaN gradient would blow chi2/dof up or move the peak out of the ROI.
  BOOST_CHECK_MESSAGE( std::isfinite(chi2dof) && (chi2dof > 0.0),
                       label << " chi2/dof not finite/positive: " << chi2dof );
  BOOST_CHECK_MESSAGE( chi2dof < chi2dof_ceiling,
                       label << " chi2/dof=" << chi2dof << " exceeds " << chi2dof_ceiling );
  BOOST_CHECK( std::isfinite(mean) && (mean > roi_lower) && (mean < roi_upper) );
  BOOST_CHECK_MESSAGE( std::fabs(mean - 661.657) < 5.0, label << " mean off: " << mean );
  BOOST_CHECK( std::isfinite(sigma) && (sigma > 0.0) );
  BOOST_CHECK( std::isfinite(amp) && (amp > 0.0) );

  // Every fitted skew coefficient must be finite and within its [lower, upper] bound.
  for( int i = 0; i < 4; ++i )
  {
    const PeakDef::CoefficientType ct = static_cast<PeakDef::CoefficientType>( PeakDef::SkewPar0 + i );
    double lo = 0.0, up = 0.0, start = 0.0, step = 0.0;
    if( !PeakDef::skew_parameter_range( skew_type, ct, lo, up, start, step ) )
      continue;
    const double v = peak->coefficient( ct );
    BOOST_CHECK_MESSAGE( std::isfinite(v) && (v >= (lo - 1.0E-6)) && (v <= (up + 1.0E-6)),
                         label << " skew par " << i << " = " << v << " outside [" << lo << "," << up << "]" );
  }
}//check_skew_peak(...)
}//namespace


BOOST_AUTO_TEST_CASE( CrystalBall_CZT_Cs137_fit )
{
  const shared_ptr<const SpecUtils::Measurement> data = load_czt_cs137();
  cout << "Loaded CZT Cs137 (H3D M400): " << data->num_gamma_channels() << " channels, "
       << data->gamma_energy_min() << "-" << data->gamma_energy_max() << " keV; is_high_res="
       << PeakFitUtils::is_high_res( data ) << endl;

  // Double-Sided Crystal Ball over [648, 672] keV.  Project owner reported chi2/dof ~ 1.5 in the
  //  GUI; this harness gets ~1.65 with a Linear continuum (minor ROI-edge/dof bookkeeping aside).
  const shared_ptr<const PeakDef> dscb
                      = fit_skew_peak( data, PeakDef::DoubleSidedCrystalBall, 648.0, 672.0 );
  check_skew_peak( dscb, PeakDef::DoubleSidedCrystalBall, 648.0, 672.0, 2.5, "DSCB [648,672]" );

  // Single-sided Crystal Ball over [648, 668.5] keV.  Project owner reported chi2/dof ~ 1.27; this
  //  harness gets ~1.4 with a Linear continuum.
  const shared_ptr<const PeakDef> cb
                      = fit_skew_peak( data, PeakDef::CrystalBall, 648.0, 668.5 );
  check_skew_peak( cb, PeakDef::CrystalBall, 648.0, 668.5, 2.0, "CB [648,668.5]" );
}//BOOST_AUTO_TEST_CASE( CrystalBall_CZT_Cs137_fit )
