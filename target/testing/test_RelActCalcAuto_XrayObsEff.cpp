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

/* With `Options::lorentzian_xrays`, the rel. eff. chart points ("observed efficiencies") of x-ray peaks must be
 found with the x-rays' own (Voigt) shapes.  Fluorescence x-rays have a free activity, so where they dominate a
 point, the freely-fit area must equal the solution's.  Before this was fixed, the free fit used the Gaussian-core
 shape of `Options::skew_type` (losing the Lorentzian wings), with skew parameters copied from the ROI's largest
 peak - here a Voigt x-ray, whose first skew parameter is its Lorentzian width.
 Needs the InterSpec "data" directory (`-- --datadir=...`).
 */

#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <iostream>

#define BOOST_TEST_MODULE RelActCalcAuto_XrayObsEff_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/XRayWidthServer.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;


namespace
{
  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = framework::master_test_suite().argc;
    char ** const argv = framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
    }

    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
        {
          datadir = d;
          break;
        }
      }
    }

    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path( datadir, "sandia.decay.xml" ) ),
                           "sandia.decay.xml not in '" << datadir << "'; pass -- --datadir=..." );
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
    BOOST_REQUIRE( DecayDataBaseServer::database() );
  }//set_data_dir()
}//namespace


BOOST_AUTO_TEST_CASE( fluorescence_xray_obs_eff_uses_voigt_shape )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Element * const uranium = db->element( "U" );
  BOOST_REQUIRE( uranium );

  // A planar-HPGe-like spectrum: 0.075 keV channels, 0.55 keV FWHM, with the U K-alpha fluorescence x-rays drawn
  //  as Voigt profiles with InterSpec's natural widths (U K-alpha HWHM ~50 eV), on a sloped continuum.
  const size_t nchannel = 4096;
  const float fwhm = 0.55f;
  const double sigma = fwhm / PhysicalUnits::fwhm_nsigma;
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( nchannel, { 0.0f, 0.075f }, {} );
  const shared_ptr<const vector<float>> energies = cal->channel_energies();
  BOOST_REQUIRE( energies && (energies->size() == nchannel + 1) );

  vector<double> model( nchannel, 0.0 );
  for( size_t i = 0; i < nchannel; ++i )
    model[i] = 2000.0 + 10.0*(150.0 - cal->energy_for_channel( i ));

  size_t num_xrays = 0;
  for( const SandiaDecay::EnergyIntensityPair &xray : uranium->xrays )
  {
    if( (xray.energy < 91.0) || (xray.energy > 100.0) )
      continue;

    const double hwhm = XRayWidths::get_xray_lorentzian_width( uranium, xray.energy, 0.5 );
    BOOST_REQUIRE_MESSAGE( hwhm > 0.0, "No natural width for U x-ray at " << xray.energy << " keV" );

    PeakDef peak( xray.energy, sigma, 1.0E7 * xray.intensity );
    peak.setSkewType( PeakDef::SkewType::VoigtPlusBortel );
    peak.set_coefficient( hwhm, PeakDef::CoefficientType::SkewPar0 );
    peak.set_coefficient( 0.0, PeakDef::CoefficientType::SkewPar1 );   // R = 0: a pure Voigt
    peak.set_coefficient( 1.0, PeakDef::CoefficientType::SkewPar2 );
    peak.gauss_integral( energies->data(), model.data(), nchannel );
    ++num_xrays;
  }//for( loop over U x-rays )
  BOOST_REQUIRE( num_xrays >= 2 );  // K-alpha1 and K-alpha2

  auto counts = make_shared<vector<float>>( nchannel, 0.0f );
  for( size_t i = 0; i < nchannel; ++i )
    (*counts)[i] = static_cast<float>( model[i] );

  auto spectrum = make_shared<SpecUtils::Measurement>();
  spectrum->set_gamma_counts( counts, 1000.0f, 1000.0f );
  spectrum->set_energy_calibration( cal );

  const PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::High;

  // Flat efficiency, with the spectrum's FWHM, so only the amplitudes and continuum are fit.
  auto drf = make_shared<DetectorPeakResponse>( "FlatPlanar", "Flat efficiency" );
  drf->fromExpOfLogPowerSeries( { -3.0f }, {}, 25.0*PhysicalUnits::cm, 1.0f*PhysicalUnits::cm,
                                static_cast<float>(PhysicalUnits::MeV), 0.0f, 3000.0f,
                                DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  drf->setFwhmCoefficients( { fwhm*fwhm, 0.0f }, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
  BOOST_REQUIRE( drf->hasResolutionInfo() );

  RelActCalcAuto::Options options;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::NotApplicable;
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency;
  options.skew_type = PeakDef::SkewType::GaussPlusBortel;
  options.lorentzian_xrays = true;
  options.additional_br_uncert = 0.0;
  options.auto_profile_weak_mass_fractions = false;
  options.auto_simplify_model = false;

  RelActCalcAuto::RelEffCurveInput curve;
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  curve.rel_eff_eqn_order = 0;
  RelActCalcAuto::NucInputInfo u_input;
  u_input.source = uranium;
  u_input.peak_color_css = "rgb(0,0,255)";
  curve.nuclides.push_back( u_input );
  options.rel_eff_curves.push_back( curve );

  RelActCalcAuto::RoiRange roi;
  roi.lower_energy = 91.0;
  roi.upper_energy = 101.0;
  roi.continuum_type = PeakContinuum::OffsetType::Linear;
  roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  options.rois.push_back( roi );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, spectrum, nullptr, drf, {}, det_type, nullptr, nullptr ) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                         "solve failed: " << sol.m_error_message );

  // The x-ray peaks of the solution are Voigts, as the data are.
  size_t num_voigt = 0;
  for( const PeakDef &peak : sol.m_fit_peaks )
    num_voigt += (peak.skewType() == PeakDef::SkewType::VoigtPlusBortel);
  BOOST_REQUIRE( num_voigt >= 2 );

  BOOST_REQUIRE_EQUAL( sol.m_obs_eff_for_each_curve.size(), 1 );
  size_t num_checked = 0;
  for( const RelActCalcAuto::RelActAutoSolution::ObsEff &obs : sol.m_obs_eff_for_each_curve[0] )
  {
    if( obs.fit_peaks.empty() || !obs.fit_peaks.front().xrayElement() )
      continue;

    BOOST_TEST_MESSAGE( "X-ray point at " << obs.energy << " keV: free-fit / solution area = "
                        << obs.observed_scale_factor );
    BOOST_CHECK_MESSAGE( fabs( obs.observed_scale_factor - 1.0 ) < 0.01,
                         "X-ray point at " << obs.energy << " keV has free-fit / solution area "
                         << obs.observed_scale_factor << "; expected 1" );
    ++num_checked;
  }//for( loop over observed efficiency points )

  BOOST_CHECK_MESSAGE( num_checked >= 2, "Expected the U K-alpha1 and K-alpha2 points; found " << num_checked );
}//fluorescence_xray_obs_eff_uses_voigt_shape
