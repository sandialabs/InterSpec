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

/* End-to-end `RelActCalcAuto::solve(...)` tests of ROIs sized from the peak shape
 (`RoiRange::RangeLimitsType::LineAnchored`), automatically chosen continua, and the re-try of
 severely misfit solutions, on synthetic Cs137 and Ba133 spectra; a "Split by lines" fit of a shipped NaI
 Eu152 reference spectrum; plus a check of the shipped presets' line anchors.  Needs the InterSpec "data" directory (`-- --datadir=...`).
 */

#include "InterSpec_config.h"

#include <cmath>
#include <limits>
#include <string>
#include <memory>
#include <vector>
#include <fstream>

#define BOOST_TEST_MODULE RelActCalcAuto_LineAnchoredFit_suite
#include <boost/test/included/unit_test.hpp>

#include <boost/math/special_functions/erf.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

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
#include "InterSpec/RelActCalcAuto_Roi.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;

namespace
{
  const double sm_cs137_energy = 661.657;

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


  /** A detector response with a flat efficiency, and a HPGe-like resolution. */
  shared_ptr<DetectorPeakResponse> make_drf()
  {
    auto drf = make_shared<DetectorPeakResponse>( "TestHPGe", "Flat efficiency, HPGe-like FWHM" );
    drf->fromExpOfLogPowerSeries( { -3.0f }, {}, 25.0*PhysicalUnits::cm, 5.0f*PhysicalUnits::cm,
                                  static_cast<float>(PhysicalUnits::MeV), 0.0f, 3000.0f,
                                  DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    // (Coefficients are for energy in MeV: FWHM = sqrt(0.81 + 1.7*E), ~1.4 keV at 662 keV.)
    drf->setFwhmCoefficients( { 0.81f, 1.7f }, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
    return drf;
  }


  /** A 8192 channel HPGe-like spectrum of a Cs137 peak on a flat continuum, with (if `step_counts` is
   non-zero) a step in the continuum proportional to the peak's cumulative area - the shape a
   `FlatStepCDF` continuum describes.  Counts are the expected values (no noise), so the result
   is deterministic. */
  shared_ptr<SpecUtils::Measurement> make_spectrum( const DetectorPeakResponse &drf,
                                                    const double step_counts )
  {
    const size_t nchannel = 8192;
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( nchannel, { 0.0f, 0.366f }, {} );

    const double fwhm = drf.peakResolutionFWHM( sm_cs137_energy );
    const double sigma = fwhm / 2.35482;
    const double peak_area = 2.0E5, base_counts = 60.0;

    const auto gauss_cdf = [&]( const double energy ) -> double {
      return 0.5*boost::math::erfc( -(energy - sm_cs137_energy) / (sigma*std::sqrt(2.0)) );
    };

    auto counts = make_shared<vector<float>>( nchannel, 0.0f );
    for( size_t i = 0; i < nchannel; ++i )
    {
      const double lower = cal->energy_for_channel( i ), upper = cal->energy_for_channel( i + 1 );
      const double mid = 0.5*(lower + upper);
      const double peak = peak_area*(gauss_cdf( upper ) - gauss_cdf( lower ));
      const double step = step_counts*(1.0 - gauss_cdf( mid ));
      (*counts)[i] = static_cast<float>( base_counts + step + peak );
    }

    auto meas = make_shared<SpecUtils::Measurement>();
    meas->set_gamma_counts( counts, 300.0f, 300.0f );
    meas->set_energy_calibration( cal );
    return meas;
  }//make_spectrum(...)


  RelActCalcAuto::Options make_options()
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

    RelActCalcAuto::Options options;
    options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
    options.fwhm_form = RelActCalcAuto::FwhmForm::NotApplicable;
    options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency;
    options.skew_type = PeakDef::SkewType::NoSkew;

    RelActCalcAuto::RelEffCurveInput curve;
    curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
    curve.rel_eff_eqn_order = 0;
    RelActCalcAuto::NucInputInfo cs137;
    cs137.source = db->nuclide( "Cs137" );
    cs137.age = 5.0*PhysicalUnits::year;
    cs137.peak_color_css = "rgb(0,0,255)";
    curve.nuclides.push_back( cs137 );
    options.rel_eff_curves.push_back( curve );

    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = roi.upper_energy = sm_cs137_energy;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.auto_continuum = true;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::LineAnchored;
    options.rois.push_back( roi );

    return options;
  }//make_options()


  RelActCalcAuto::RelActAutoSolution solve( const RelActCalcAuto::Options &options,
                                            const shared_ptr<const SpecUtils::Measurement> &spectrum,
                                            const shared_ptr<const DetectorPeakResponse> &drf )
  {
    return RelActCalcAuto::solve( options, spectrum, nullptr, drf, {},
                                  PeakFitUtils::CoarseResolutionType::High, nullptr );
  }
}//namespace


BOOST_AUTO_TEST_CASE( single_line_roi_extent_from_peak_shape )
{
  set_data_dir();
  const shared_ptr<DetectorPeakResponse> drf = make_drf();
  const shared_ptr<SpecUtils::Measurement> spectrum = make_spectrum( *drf, 0.0 );

  RelActCalcAuto::Options options = make_options();
  options.rois[0].auto_continuum = false;
  const RelActCalcAuto::RelActAutoSolution sol = solve( options, spectrum, drf );

  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                         "solve failed: " << sol.m_error_message );
  BOOST_REQUIRE_EQUAL( sol.m_final_roi_ranges.size(), 1 );
  BOOST_REQUIRE_EQUAL( sol.m_final_roi_info.size(), 1 );
  BOOST_CHECK( sol.m_input_rois == options.rois );

  // The fit ROI is fixed, and extends past the line by the default coverage (0.1% of the peak area)
  //  plus 1 FWHM - to within a channel.
  const RelActCalcAuto::RoiRange &roi = sol.m_final_roi_ranges[0];
  BOOST_CHECK( roi.range_limits_type == RelActCalcAuto::RoiRange::RangeLimitsType::Fixed );
  BOOST_CHECK_EQUAL( sol.m_final_roi_info[0].lower_anchor_energy, sm_cs137_energy );
  BOOST_CHECK( sol.m_final_roi_info[0].input_roi_indices == vector<size_t>{0} );

  const double fwhm = drf->peakResolutionFWHM( sm_cs137_energy );
  const double expected_half_width = 3.0902*fwhm/2.35482 + fwhm;
  BOOST_CHECK_SMALL( (sm_cs137_energy - roi.lower_energy) - expected_half_width, 0.4 );
  BOOST_CHECK_SMALL( (roi.upper_energy - sm_cs137_energy) - expected_half_width, 0.4 );

  // The options the fit used (which re-solves are built from) have the resolved ROI.
  BOOST_REQUIRE_EQUAL( sol.m_options.rois.size(), 1 );
  BOOST_CHECK( sol.m_options.rois[0] == roi );

  // Wider edges give a wider ROI.
  RelActCalcAuto::Options wide = options;
  wide.rois[0].lower_edge.sideband_fwhm = 3.0;
  wide.rois[0].upper_edge.tail_fraction = 1.0E-5;
  const RelActCalcAuto::RelActAutoSolution wide_sol = solve( wide, spectrum, drf );
  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status( wide_sol.m_status ) );
  BOOST_REQUIRE_EQUAL( wide_sol.m_final_roi_ranges.size(), 1 );
  BOOST_CHECK_SMALL( (roi.lower_energy - wide_sol.m_final_roi_ranges[0].lower_energy) - 2.0*fwhm, 0.4 );
  BOOST_CHECK_GT( wide_sol.m_final_roi_ranges[0].upper_energy, roi.upper_energy );
}//single_line_roi_extent_from_peak_shape


BOOST_AUTO_TEST_CASE( initial_rois_sized_from_spectrum_peaks )
{
  set_data_dir();

  // A Ba133 spectrum whose peaks have the resolution of `make_drf()`, fit using a DRF whose FWHM is
  //  twice that (as when a generic DRF is used for a better-resolution detector).  Before any fit,
  //  the ROIs should be sized from the peaks found in the spectrum, not the DRF - so the ROIs sized
  //  from the fit peak shape agree with them, and no re-fit is needed.
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const ba133 = db->nuclide( "Ba133" );
  BOOST_REQUIRE( ba133 );
  const double age = 1.0*PhysicalUnits::year;

  const shared_ptr<DetectorPeakResponse> true_drf = make_drf();
  auto wide_drf = make_shared<DetectorPeakResponse>( *true_drf );
  wide_drf->setFwhmCoefficients( { 4.0f*0.81f, 4.0f*1.7f }, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  // Noiseless spectrum: flat continuum, plus a Gaussian for each Ba133 gamma.
  const size_t nchannel = 16384;
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( nchannel, { 0.0f, 0.183f }, {} );

  SandiaDecay::NuclideMixture mix;
  mix.addNuclideByActivity( ba133, 1.0 );
  const vector<SandiaDecay::EnergyRatePair> gammas
                          = mix.gammas( age, SandiaDecay::NuclideMixture::OrderByEnergy, false );
  double max_rate = 0.0;
  for( const SandiaDecay::EnergyRatePair &g : gammas )
    max_rate = std::max( max_rate, g.numPerSecond );
  BOOST_REQUIRE( max_rate > 0.0 );

  auto counts = make_shared<vector<float>>( nchannel, 60.0f );
  for( const SandiaDecay::EnergyRatePair &g : gammas )
  {
    const double area = 2.0E5 * g.numPerSecond / max_rate;
    const double sigma = true_drf->peakResolutionFWHM( static_cast<float>(g.energy) ) / 2.35482;
    const auto cdf = [&]( const double energy ) -> double {
      return 0.5*boost::math::erfc( -(energy - g.energy) / (sigma*std::sqrt(2.0)) );
    };
    for( size_t i = 0; i < nchannel; ++i )
    {
      const double lower = cal->energy_for_channel( i ), upper = cal->energy_for_channel( i + 1 );
      if( (upper > (g.energy - 10.0*sigma)) && (lower < (g.energy + 10.0*sigma)) )
        (*counts)[i] += static_cast<float>( area*(cdf( upper ) - cdf( lower )) );
    }
  }//for( const SandiaDecay::EnergyRatePair &g : gammas )

  auto spectrum = make_shared<SpecUtils::Measurement>();
  spectrum->set_gamma_counts( counts, 300.0f, 300.0f );
  spectrum->set_energy_calibration( cal );

  RelActCalcAuto::Options options = make_options();
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_2;
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
  options.rel_eff_curves[0].nuclides[0].source = ba133;
  options.rel_eff_curves[0].nuclides[0].age = age;

  options.rois.clear();
  for( const pair<double,double> &lines : vector<pair<double,double>>{ {79.614, 80.997}, {276.399, 276.399},
                                                  {302.851, 302.851}, {356.013, 356.013}, {383.849, 383.849} } )
  {
    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = lines.first;
    roi.upper_energy = lines.second;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::LineAnchored;
    options.rois.push_back( roi );
  }

  const RelActCalcAuto::RelActAutoSolution sol = solve( options, spectrum, wide_drf );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                         "solve failed: " << sol.m_error_message );
  BOOST_REQUIRE_EQUAL( sol.m_final_roi_ranges.size(), options.rois.size() );
  BOOST_REQUIRE_EQUAL( sol.m_final_roi_info.size(), options.rois.size() );
  BOOST_CHECK_EQUAL( sol.m_num_roi_refits, size_t(0) );

  // And the ROIs are sized for the spectrum's peaks: the default 0.1% coverage, plus 1 FWHM.
  for( size_t i = 0; i < sol.m_final_roi_ranges.size(); ++i )
  {
    const RelActCalcAuto::RoiRange &roi = sol.m_final_roi_ranges[i];
    const RelActCalcAuto::RoiResolutionInfo &info = sol.m_final_roi_info[i];
    const double lower_fwhm = true_drf->peakResolutionFWHM( static_cast<float>(info.lower_anchor_energy) );
    const double upper_fwhm = true_drf->peakResolutionFWHM( static_cast<float>(info.upper_anchor_energy) );
    BOOST_CHECK_SMALL( (info.lower_anchor_energy - roi.lower_energy) - (3.0902/2.35482 + 1.0)*lower_fwhm, 0.3 );
    BOOST_CHECK_SMALL( (roi.upper_energy - info.upper_anchor_energy) - (3.0902/2.35482 + 1.0)*upper_fwhm, 0.3 );
  }
}//initial_rois_sized_from_spectrum_peaks


BOOST_AUTO_TEST_CASE( split_by_lines_on_low_resolution_many_line_source )
{
  // A shielded Eu152 NaI spectrum (one of the reference spectra InterSpec ships), fit with a single
  //  "Split by lines" range.  Windowing every Eu152 line before the first fit would merge them into one
  //  ROI spanning nearly the whole spectrum, whose continuum the fit could only follow by making the
  //  peaks many times too wide; windowing just the lines at peaks found in the spectrum keeps the ROIs,
  //  and so the fit, sensible.
  set_data_dir();
  const string dir = SpecUtils::append_path( InterSpec::staticDataDirectory(),
                                             "reference_spectra/Common_Field_Nuclides/IdentiFINDER-R500-NaI" );
  const auto load_spectrum = [&dir]( const string &filename ) -> shared_ptr<const SpecUtils::Measurement> {
    const string path = SpecUtils::append_path( dir, filename );
    SpecUtils::SpecFile spec;
    BOOST_REQUIRE_MESSAGE( spec.load_file( path, SpecUtils::ParserType::Auto, path ), "Could not load " << path );
    return spec.sum_measurements( spec.sample_numbers(), spec.detector_names(), nullptr );
  };
  const shared_ptr<const SpecUtils::Measurement> foreground = load_spectrum( "Eu152_Shielded.txt" );
  const shared_ptr<const SpecUtils::Measurement> background = load_spectrum( "background.txt" );
  BOOST_REQUIRE( foreground && background );

  shared_ptr<const DetectorPeakResponse> drf;
  {
    ifstream tsv( SpecUtils::append_path( InterSpec::staticDataDirectory(), "common_drfs.tsv" ).c_str(), ios::binary );
    vector<string> credits;
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    DetectorPeakResponse::parseMultipleRelEffDrfCsv( tsv, credits, drfs );
    for( const shared_ptr<DetectorPeakResponse> &d : drfs )
    {
      if( d && (d->name() == "IdentiFINDER-R500-NaI") )
        drf = d;
    }
  }
  BOOST_REQUIRE( drf && drf->hasResolutionInfo() );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  RelActCalcAuto::Options options = make_options();
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::LinearFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Berstein_4;
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
  options.rel_eff_curves[0].rel_eff_eqn_order = 3;
  options.rel_eff_curves[0].nuclides[0].source = db->nuclide( "Eu152" );
  options.rel_eff_curves[0].nuclides[0].age = 1.0*PhysicalUnits::year;
  options.rois[0].lower_energy = 50.0;
  options.rois[0].upper_energy = 3000.0;
  options.rois[0].range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::CanBeBrokenUp;

  const RelActCalcAuto::RelActAutoSolution sol = RelActCalcAuto::solve( options, foreground, background, drf, {},
                                                      PeakFitUtils::CoarseResolutionType::Low, nullptr );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                         "solve failed: " << sol.m_error_message );

  BOOST_CHECK_GE( sol.m_final_roi_ranges.size(), size_t(3) );
  for( const RelActCalcAuto::RoiRange &roi : sol.m_final_roi_ranges )
    BOOST_CHECK_LT( roi.upper_energy - roi.lower_energy, 1000.0 );

  BOOST_REQUIRE( sol.m_dof_data > 0 );
  BOOST_CHECK_MESSAGE( (sol.m_chi2_data / sol.m_dof_data) < 3.0,
                       "chi2/dof=" << (sol.m_chi2_data / sol.m_dof_data) );

  // The fit peak width at 344 keV agrees with the detector's.
  const double ref_energy = 344.279;
  const PeakDef *ref_peak = nullptr;
  for( const PeakDef &peak : sol.m_fit_peaks )
  {
    if( fabs( peak.mean() - ref_energy ) < 1.0 )
      ref_peak = &peak;
  }
  BOOST_REQUIRE( ref_peak );
  const double drf_fwhm = drf->peakResolutionFWHM( static_cast<float>(ref_energy) );
  BOOST_CHECK_MESSAGE( fabs( ref_peak->fwhm() - drf_fwhm ) < 0.25*drf_fwhm,
                       "FWHM at 344 keV is " << ref_peak->fwhm() << " keV; the DRF's is " << drf_fwhm );
}//split_by_lines_on_low_resolution_many_line_source


BOOST_AUTO_TEST_CASE( auto_continuum_picks_step_only_when_there_is_one )
{
  set_data_dir();
  const shared_ptr<DetectorPeakResponse> drf = make_drf();
  const RelActCalcAuto::Options options = make_options();

  // A continuum stepping down by 40 counts/channel across the peak: with enough continuum on either
  //  side, a step continuum is chosen.  (Within ~1 FWHM of the peak, as the default ROI extent gives,
  //  a straight line describes a step as well as a step does, and so is correctly kept.)
  RelActCalcAuto::Options wide_options = options;
  wide_options.rois[0].lower_edge.sideband_fwhm = 4.0;
  wide_options.rois[0].upper_edge.sideband_fwhm = 4.0;
  const RelActCalcAuto::RelActAutoSolution step_sol = solve( wide_options, make_spectrum( *drf, 40.0 ), drf );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( step_sol.m_status ),
                         "solve failed: " << step_sol.m_error_message );
  BOOST_REQUIRE_EQUAL( step_sol.m_final_roi_ranges.size(), 1 );
  BOOST_REQUIRE_EQUAL( step_sol.m_final_roi_info.size(), 1 );
  BOOST_CHECK( step_sol.m_final_roi_info[0].auto_continuum );
  string warnings;
  for( const string &w : step_sol.m_warnings )
    warnings += "\n    " + w;
  BOOST_CHECK_MESSAGE( PeakContinuum::is_step_continuum( step_sol.m_final_roi_ranges[0].continuum_type ),
                       "chose " << PeakContinuum::offset_type_str( step_sol.m_final_roi_ranges[0].continuum_type )
                       << "; warnings:" << warnings );

  // A flat continuum: the starting (linear) continuum is kept.
  const RelActCalcAuto::RelActAutoSolution flat_sol = solve( wide_options, make_spectrum( *drf, 0.0 ), drf );
  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status( flat_sol.m_status ) );
  BOOST_REQUIRE_EQUAL( flat_sol.m_final_roi_ranges.size(), 1 );
  BOOST_CHECK_MESSAGE( flat_sol.m_final_roi_ranges[0].continuum_type == PeakContinuum::OffsetType::Linear,
                       "chose " << PeakContinuum::offset_type_str( flat_sol.m_final_roi_ranges[0].continuum_type ) );
}//auto_continuum_picks_step_only_when_there_is_one


BOOST_AUTO_TEST_CASE( severe_misfit_gets_full_candidate_search )
{
  // A fit left far from describing the data (chi2/dof > 100) is re-tried from the starting points a
  //  robust solve uses; a fit that describes the data is not.
  set_data_dir();
  const shared_ptr<DetectorPeakResponse> drf = make_drf();

  RelActCalcAuto::Options options = make_options();
  options.rois[0].range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  options.rois[0].lower_energy = sm_cs137_energy - 10.0;
  options.rois[0].upper_energy = sm_cs137_energy + 10.0;
  options.rois[0].auto_continuum = false;

  const auto searched_for_misfit = []( const RelActCalcAuto::RelActAutoSolution &sol ) -> bool {
    for( const string &warning : sol.m_warnings )
    {
      if( SpecUtils::icontains( warning, "candidate search" ) && SpecUtils::icontains( warning, "severe misfit" ) )
        return true;
    }
    return false;
  };

  // A 20000 counts/channel continuum step, which a straight line can not describe.
  const RelActCalcAuto::RelActAutoSolution misfit_sol = solve( options, make_spectrum( *drf, 20000.0 ), drf );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( misfit_sol.m_status ),
                         "solve failed: " << misfit_sol.m_error_message );
  BOOST_REQUIRE( misfit_sol.m_dof_data > 0 );
  BOOST_CHECK_MESSAGE( (misfit_sol.m_chi2_data / misfit_sol.m_dof_data) > 100.0,
                       "chi2/dof=" << (misfit_sol.m_chi2_data / misfit_sol.m_dof_data) );
  BOOST_CHECK( searched_for_misfit( misfit_sol ) );

  // The same fit of a flat continuum.
  const RelActCalcAuto::RelActAutoSolution good_sol = solve( options, make_spectrum( *drf, 0.0 ), drf );
  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status( good_sol.m_status ) );
  BOOST_CHECK( !searched_for_misfit( good_sol ) );
}//severe_misfit_gets_full_candidate_search


BOOST_AUTO_TEST_CASE( preset_line_anchors_are_source_lines )
{
  // The energies a line-anchored ROI in a shipped preset is defined by must be photon lines of that
  //  preset's sources (in their decay chains, at their ages) - catching typos, and changes to the
  //  nuclear data.
  set_data_dir();
  const string preset_dir = SpecUtils::append_path( InterSpec::staticDataDirectory(), "rel_act" );
  const vector<string> presets = SpecUtils::recursive_ls( preset_dir, ".xml" );
  BOOST_REQUIRE( !presets.empty() );

  for( const string &preset : presets )
  {
    std::vector<char> data;
    SpecUtils::load_file_data( preset.c_str(), data );
    rapidxml::xml_document<char> doc;
    BOOST_REQUIRE_NO_THROW( doc.parse<rapidxml::parse_trim_whitespace>( &data[0] ) );
    RelActCalcAuto::RelActAutoGuiState state;
    BOOST_REQUIRE_NO_THROW( state.deSerialize( doc.first_node() ) );

    vector<double> line_energies;
    for( const RelActCalcAuto::RelEffCurveInput &curve : state.options.rel_eff_curves )
    {
      for( const RelActCalcAuto::NucInputInfo &nuc_input : curve.nuclides )
      {
        const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( nuc_input.source );
        if( !nuc )
          continue;
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nuc, 1.0*SandiaDecay::becquerel, nuc_input.age );
        for( const SandiaDecay::EnergyRatePair &line : mix.photons( 0.0 ) )
        {
          if( line.numPerSecond > 1.0E-7 )
            line_energies.push_back( line.energy );
        }
      }
    }//for( loop over rel. eff. curves )

    for( const RelActCalcAuto::RoiRange &roi : state.options.rois )
    {
      if( roi.range_limits_type != RelActCalcAuto::RoiRange::RangeLimitsType::LineAnchored )
        continue;

      for( const double anchor : { roi.lower_energy, roi.upper_energy } )
      {
        double nearest = std::numeric_limits<double>::max();
        for( const double energy : line_energies )
          nearest = std::min( nearest, fabs( energy - anchor ) );
        BOOST_CHECK_MESSAGE( nearest < 0.05, "In '" << SpecUtils::filename(preset) << "', the ROI anchored at "
                             << anchor << " keV is " << nearest << " keV from the nearest source line." );
      }
    }//for( loop over ROIs )
  }//for( const string &preset : presets )
}//preset_line_anchors_are_source_lines


BOOST_AUTO_TEST_CASE( options_with_new_roi_features_round_trip )
{
  set_data_dir();

  RelActCalcAuto::RelActAutoGuiState state;
  state.options = make_options();
  state.options.roi_settings.continuum_switch_min_improvement = 4.0;
  state.options.rois[0].upper_edge.sideband_fwhm = 1.5;

  rapidxml::xml_document<char> doc;
  BOOST_REQUIRE_NO_THROW( state.serialize( &doc ) );
  string xml;
  rapidxml::print( std::back_inserter(xml), doc, 0 );
  BOOST_CHECK( xml.find( "<Options version=\"7\"" ) != string::npos );
  BOOST_CHECK( xml.find( "RoiSettings" ) != string::npos );

  RelActCalcAuto::RelActAutoGuiState back;
  rapidxml::xml_document<char> doc2;
  doc2.parse<rapidxml::parse_trim_whitespace>( &xml[0] );
  BOOST_REQUIRE_NO_THROW( back.deSerialize( doc2.first_node() ) );
  BOOST_CHECK( back.options.roi_settings == state.options.roi_settings );
  BOOST_REQUIRE_EQUAL( back.options.rois.size(), 1 );
  BOOST_CHECK( back.options.rois[0] == state.options.rois[0] );

  // Default settings, and only fixed ROIs, keep the older (version 6 or less) format.
  RelActCalcAuto::RelActAutoGuiState old_style = state;
  old_style.options.roi_settings = RelActCalcAuto::Options::RoiSettings{};
  old_style.options.rois[0].range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  old_style.options.rois[0].upper_energy = sm_cs137_energy + 5.0;
  old_style.options.rois[0].auto_continuum = false;
  old_style.options.rois[0].upper_edge = RelActCalcAuto::RoiRange::EdgeOverride{};
  rapidxml::xml_document<char> doc3;
  old_style.serialize( &doc3 );
  string old_xml;
  rapidxml::print( std::back_inserter(old_xml), doc3, 0 );
  BOOST_CHECK( old_xml.find( "<Options version=\"7\"" ) == string::npos );
  BOOST_CHECK( old_xml.find( "RoiSettings" ) == string::npos );
}//options_with_new_roi_features_round_trip


BOOST_AUTO_TEST_CASE( default_roi_edges_depend_on_detector_type )
{
  using RelActCalcAuto::RoiEdge;
  using PeakFitUtils::CoarseResolutionType;

  // HPGe: 99.9% of the peak plus 1 FWHM of continuum, either side.
  BOOST_CHECK( RelActCalcAuto::default_roi_edge( CoarseResolutionType::High, true ) == (RoiEdge{ 1.0E-3, 1.0 }) );
  BOOST_CHECK( RelActCalcAuto::default_roi_edge( CoarseResolutionType::High, false ) == (RoiEdge{ 1.0E-3, 1.0 }) );

  // NaI and LaBr: half a FWHM of continuum.
  for( const CoarseResolutionType type : { CoarseResolutionType::Low, CoarseResolutionType::LaBr,
                                           CoarseResolutionType::LowOrMedRes } )
  {
    BOOST_CHECK( RelActCalcAuto::default_roi_edge( type, true ) == (RoiEdge{ 1.0E-3, 0.5 }) );
    BOOST_CHECK( RelActCalcAuto::default_roi_edge( type, false ) == (RoiEdge{ 1.0E-3, 0.5 }) );
  }

  // CZT: less of the peak below the line, where its long tail is.
  BOOST_CHECK( RelActCalcAuto::default_roi_edge( CoarseResolutionType::CZT, true ) == (RoiEdge{ 2.0E-2, 0.5 }) );
  BOOST_CHECK( RelActCalcAuto::default_roi_edge( CoarseResolutionType::CZT, false ) == (RoiEdge{ 1.0E-3, 0.5 }) );
}//default_roi_edges_depend_on_detector_type


