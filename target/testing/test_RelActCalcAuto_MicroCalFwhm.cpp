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

/* `RelActCalcAuto::solve(...)` should fit the ~72 eV FWHM peaks of a finely binned (10 eV/channel)
 micro-calorimeter spectrum, instead of rejecting the FWHM found from the spectrum's peaks as narrower
 than any HPGe, and falling back to default HPGe widths.
 Needs the InterSpec "data" directory (`-- --datadir=...`).
 */

#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <iostream>

#define BOOST_TEST_MODULE RelActCalcAuto_MicroCalFwhm_suite
#include <boost/test/included/unit_test.hpp>

#include <boost/math/special_functions/erf.hpp>

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
#include "InterSpec/PeakFitDetPrefs.h"
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


BOOST_AUTO_TEST_CASE( micro_cal_fwhm_fit_from_spectrum_peaks )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  BOOST_REQUIRE( u235 );
  const double age = 20.0*PhysicalUnits::year;

  // 32768 channels of 10 eV (0 to 327.68 keV), with 72 eV FWHM peaks, like LANL's SOFIA spectra.
  const size_t nchannel = 32768;
  const double true_fwhm = 0.072;
  const double sigma = true_fwhm / PhysicalUnits::fwhm_nsigma;
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( nchannel, { 0.0f, 0.01f }, {} );

  // Noiseless spectrum: a gently falling continuum, plus a Gaussian for each U235 (and progeny)
  //  photon between 120 and 215 keV, with areas proportional to their yield (i.e., flat efficiency);
  //  the 143.77, 163.36, 185.72, and 205.32 keV lines dominate.
  SandiaDecay::NuclideMixture mix;
  mix.addNuclideByActivity( u235, 1.0 );
  const vector<SandiaDecay::EnergyRatePair> photons = mix.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );

  double rate_185 = 0.0, energy_185 = 0.0;
  for( const SandiaDecay::EnergyRatePair &p : photons )
  {
    if( (fabs( p.energy - 185.72 ) < 0.1) && (p.numPerSecond > rate_185) )
    {
      rate_185 = p.numPerSecond;
      energy_185 = p.energy;
    }
  }
  BOOST_REQUIRE( rate_185 > 0.0 );

  auto counts = make_shared<vector<float>>( nchannel, 0.0f );
  for( size_t i = 0; i < nchannel; ++i )
    (*counts)[i] = static_cast<float>( 30.0 + 0.1*(327.68 - cal->energy_for_channel( i )) );

  for( const SandiaDecay::EnergyRatePair &p : photons )
  {
    if( (p.energy < 120.0) || (p.energy > 215.0) )
      continue;

    const double area = 2.0E5 * p.numPerSecond / rate_185;
    const auto cdf = [&]( const double energy ) -> double {
      return 0.5*boost::math::erfc( -(energy - p.energy) / (sigma*std::sqrt(2.0)) );
    };

    const size_t first_channel = static_cast<size_t>( cal->channel_for_energy( p.energy - 10.0*sigma ) );
    const size_t last_channel = static_cast<size_t>( cal->channel_for_energy( p.energy + 10.0*sigma ) );
    for( size_t i = first_channel; (i <= last_channel) && (i < nchannel); ++i )
    {
      const double lower = cal->energy_for_channel( i ), upper = cal->energy_for_channel( i + 1 );
      (*counts)[i] += static_cast<float>( area*(cdf( upper ) - cdf( lower )) );
    }
  }//for( const SandiaDecay::EnergyRatePair &p : photons )

  auto spectrum = make_shared<SpecUtils::Measurement>();
  spectrum->set_gamma_counts( counts, 1000.0f, 1000.0f );
  spectrum->set_energy_calibration( cal );

  // The relaxed width limits are only for high-resolution detectors, which real micro-calorimeter
  //  spectra are classified as - but this sparse synthetic spectrum may not be, so we specify it.
  const PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::High;
  auto fit_prefs = make_shared<PeakFitDetPrefs>();
  fit_prefs->m_det_type = det_type;

  // A flat efficiency, with no FWHM information - so the FWHM must come from the spectrum's peaks.
  auto drf = make_shared<DetectorPeakResponse>( "FlatMicroCal", "Flat efficiency, no FWHM info" );
  drf->fromExpOfLogPowerSeries( { -3.0f }, {}, 25.0*PhysicalUnits::cm, 1.0f*PhysicalUnits::cm,
                                static_cast<float>(PhysicalUnits::MeV), 0.0f, 3000.0f,
                                DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  drf->setPeakFitDetPrefs( fit_prefs );
  BOOST_REQUIRE( !drf->hasResolutionInfo() );

  RelActCalcAuto::Options options;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_2;
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartingFromAllPeaksInSpectrum;
  options.skew_type = PeakDef::SkewType::NoSkew;
  options.additional_br_uncert = 0.0;
  options.auto_profile_weak_mass_fractions = false;
  options.auto_simplify_model = false;

  RelActCalcAuto::RelEffCurveInput curve;
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  curve.rel_eff_eqn_order = 0;
  RelActCalcAuto::NucInputInfo u235_input;
  u235_input.source = u235;
  u235_input.age = age;
  u235_input.peak_color_css = "rgb(0,0,255)";
  curve.nuclides.push_back( u235_input );
  options.rel_eff_curves.push_back( curve );

  for( const pair<double,double> &range : vector<pair<double,double>>{ {140.0, 148.0}, {160.0, 168.0}, {180.0, 210.0} } )
  {
    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = range.first;
    roi.upper_energy = range.second;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    options.rois.push_back( roi );
  }

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, spectrum, nullptr, drf, {}, det_type, nullptr, fit_prefs ) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                         "solve failed: " << sol.m_error_message );

  for( const string &warning : sol.m_warnings )
  {
    BOOST_TEST_MESSAGE( "Warning: " << warning );
    BOOST_CHECK_MESSAGE( warning.find( "Failed to estimate FWHM" ) == string::npos,
                         "Unexpected warning: " << warning );
  }

  const double fit_fwhm = RelActCalcAuto::eval_fwhm( energy_185, sol.m_fwhm_form, sol.m_fwhm_coefficients );
  BOOST_TEST_MESSAGE( "Fit FWHM at 185.7 keV: " << fit_fwhm << " keV" );
  BOOST_CHECK_MESSAGE( fabs( fit_fwhm - true_fwhm ) < 0.05*true_fwhm,
                       "Fit FWHM at 185.7 keV is " << fit_fwhm << " keV; expected " << true_fwhm );

  // The fit 185.7 keV peak should have the same width.
  const PeakDef *peak_185 = nullptr;
  for( const PeakDef &peak : sol.m_fit_peaks )
  {
    if( !peak_185 || (fabs( peak.mean() - energy_185 ) < fabs( peak_185->mean() - energy_185 )) )
      peak_185 = &peak;
  }
  BOOST_REQUIRE( peak_185 );
  BOOST_CHECK_SMALL( peak_185->mean() - energy_185, 0.002 );
  BOOST_CHECK_MESSAGE( fabs( peak_185->fwhm() - true_fwhm ) < 0.05*true_fwhm,
                       "Fit 185.7 keV peak FWHM is " << peak_185->fwhm() << " keV; expected " << true_fwhm );
}//micro_cal_fwhm_fit_from_spectrum_peaks
