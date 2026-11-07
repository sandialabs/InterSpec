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
#include <deque>
#include <regex>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <optional>
#include <functional>

#define BOOST_TEST_MODULE test_SpectrumChartPeaks
#include <boost/test/included/unit_test.hpp>

#include <nlohmann/json.hpp>

#include <Wt/WPainter.h>
#include <Wt/WSvgImage.h>
#include <Wt/WPaintDevice.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/SpectrumChartImage.h"

using namespace std;

// Tests the peak drawing of `SpectrumChart` (a port of SpectrumChartD3.js's), and smoke-tests
//  rendering whole charts.  Pass `--outdir=<dir>` (after `--`) to write the rendered SVGs.

namespace
{
  const size_t sm_num_channels = 2048;
  const float sm_channel_width = 0.5f;  //keV

  /** The `--outdir=` command line argument, if given. */
  const std::string &output_dir()
  {
    static const std::string dir = [](){
      const int argc = boost::unit_test::framework::master_test_suite().argc;
      char ** const argv = boost::unit_test::framework::master_test_suite().argv;
      for( int i = 1; i < argc; ++i )
      {
        const string arg = argv[i];
        if( SpecUtils::istarts_with( arg, "--outdir=" ) )
          return arg.substr( 9 );
      }
      return string();
    }();

    return dir;
  }//output_dir()


  /** The counts of a Gaussian between `lo` and `hi`, computed independently of `PeakDef`. */
  double gaussian_area( const double mean, const double sigma, const double amp,
                        const double lo, const double hi )
  {
    const double scale = sigma * std::sqrt( 2.0 );
    return 0.5 * amp * (std::erfc( (lo - mean)/scale ) - std::erfc( (hi - mean)/scale ));
  }


  /** The counts of the linear continuum density `c0 + c1*(x - ref)` between `lo` and `hi`. */
  double linear_area( const double c0, const double c1, const double ref,
                      const double lo, const double hi )
  {
    return c0*(hi - lo) + 0.5*c1*((hi - ref)*(hi - ref) - (lo - ref)*(lo - ref));
  }


  /** A spectrum of a flat 100 counts per channel, plus a peak at 500 keV. */
  shared_ptr<SpecUtils::Measurement> make_spectrum()
  {
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( sm_num_channels, { 0.0f, sm_channel_width }, {} );

    auto counts = make_shared<vector<float>>( sm_num_channels, 100.0f );
    const PeakDef peak( 500.3, 1.2, 2.0E4 );
    for( size_t i = 0; i < sm_num_channels; ++i )
      (*counts)[i] += static_cast<float>( peak.gauss_integral( i*sm_channel_width, (i + 1)*sm_channel_width ) );

    auto meas = make_shared<SpecUtils::Measurement>();
    meas->set_gamma_counts( counts, 300.0f, 310.0f );
    meas->set_energy_calibration( cal );
    meas->set_title( "Test spectrum" );
    return meas;
  }//make_spectrum()


  /** The displayed channels of `meas`, combining `rebin` channels each, as `SpectrumDataModel` does. */
  void displayed_channels( const shared_ptr<const SpecUtils::Measurement> &meas, const size_t rebin,
                           vector<double> &edges, vector<double> &counts )
  {
    edges.clear();
    counts.clear();
    const size_t nrows = meas->num_gamma_channels() / rebin;
    for( size_t row = 0; row < nrows; ++row )
    {
      edges.push_back( meas->gamma_channel_lower( row*rebin ) );
      double sum = 0.0;
      for( size_t i = 0; i < rebin; ++i )
        sum += meas->gamma_channel_content( row*rebin + i );
      counts.push_back( sum );
    }
    edges.push_back( meas->gamma_channel_upper( nrows*rebin - 1 ) );
  }//displayed_channels(...)


  /** A ROI with a single Gaussian on a linear continuum. */
  vector<shared_ptr<const PeakDef>> make_linear_roi( const double lower, const double upper )
  {
    auto peak = make_shared<PeakDef>( 500.3, 1.2, 2.0E4 );
    shared_ptr<PeakContinuum> cont = peak->getContinuum();
    cont->setType( PeakContinuum::Linear );
    cont->setParameters( 495.0, vector<double>{ 200.0, -2.0 }, vector<double>{} );
    cont->setRange( lower, upper );
    return { peak };
  }//make_linear_roi(...)


  /** Calls `func(i, lo, hi)` for each sample exactly at a displayed-channel center, with the
   channels limits; returns how many there were.
   */
  size_t for_each_center_sample( const SpectrumChart::RoiCurves &curves, const vector<double> &edges,
                                 const std::function<void(size_t,double,double)> &func )
  {
    size_t nchecked = 0;
    for( size_t i = 1; (i + 1) < curves.energies.size(); ++i )  // (not the end samples, which are moved)
    {
      const double x = curves.energies[i];
      const auto pos = std::upper_bound( begin(edges), end(edges), x );
      if( (pos == begin(edges)) || (pos == end(edges)) )
        continue;
      const double lo = *(pos - 1), hi = *pos;
      if( fabs( x - 0.5*(lo + hi) ) > 1.0E-9*std::max( 1.0, fabs(x) ) )
        continue;
      func( i, lo, hi );
      ++nchecked;
    }
    return nchecked;
  }//for_each_center_sample(...)


  /** Some skew parameter values, for each skew type, within their valid ranges. */
  void set_reasonable_skew( PeakDef &peak, const PeakDef::SkewType type )
  {
    peak.setSkewType( type );
    for( size_t i = 0; i < PeakDef::num_skew_parameters( type ); ++i )
    {
      const auto coef = PeakDef::CoefficientType( PeakDef::SkewPar0 + i );
      double lower, upper, starting, step;
      if( PeakDef::skew_parameter_range( type, coef, lower, upper, starting, step ) )
        peak.set_coefficient( starting, coef );
    }
  }//set_reasonable_skew(...)


  std::string render_svg( const SpectrumChartImage::SpectrumOptions &options, const int width,
                          const int height, const std::string &name )
  {
    const std::string svg = SpectrumChartImage::renderImage( width, height, "svg",
                              [&options]( Wt::WPaintDevice &device ){
      SpectrumChartImage::paintSpectrum( device, options );
    } );

    if( !output_dir().empty() )
    {
      std::ofstream out( SpecUtils::append_path( output_dir(), name + ".svg" ), ios::binary );
      out << svg;
    }

    return svg;
  }//render_svg(...)
}//namespace


BOOST_AUTO_TEST_CASE( ValuesAtChannelCenters )
{
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  const vector<shared_ptr<const PeakDef>> roi = make_linear_roi( 490.0, 511.0 );
  const PeakDef &peak = *roi[0];

  // With the channels as they are, and combined in fours, the curves pass through the counts the
  //  displayed channel centered there holds.
  for( const size_t rebin : { size_t(1), size_t(4) } )
  {
    vector<double> edges, counts;
    displayed_channels( meas, rebin, edges, counts );

    SpectrumChart::RoiCurves curves;
    BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 480.0, 520.0, 25.0, curves ) );
    BOOST_REQUIRE( curves.peaks.size() == 1 );
    BOOST_REQUIRE( curves.continuum.size() == curves.energies.size() );

    const SpectrumChart::RoiCurves::PeakCurve &pc = curves.peaks[0];
    BOOST_CHECK_EQUAL( pc.first, 0 );  // A lone peak is drawn across its whole ROI
    BOOST_CHECK_EQUAL( pc.last, curves.energies.size() - 1 );

    const size_t nchecked = for_each_center_sample( curves, edges, [&]( size_t i, double lo, double hi ){
      const double cont = std::max( 0.0, linear_area( 200.0, -2.0, 495.0, lo, hi ) );
      const double pk = gaussian_area( peak.mean(), peak.sigma(), peak.amplitude(), lo, hi );
      BOOST_CHECK_SMALL( curves.continuum[i] - cont, 1.0E-6*std::max( 1.0, cont ) );
      BOOST_CHECK_SMALL( pc.counts[i] - pk, 1.0E-6*std::max( 1.0, pk ) );
    } );
    BOOST_CHECK_GT( nchecked, (rebin == 1) ? size_t(30) : size_t(6) );

    // Samples get closer than a channel, when channels are many pixels wide
    BOOST_CHECK_GT( curves.energies.size(), 2*nchecked );
  }//for( rebin of 1 and 4 )
}//BOOST_AUTO_TEST_CASE( ValuesAtChannelCenters )


BOOST_AUTO_TEST_CASE( RoiEdges )
{
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  vector<double> edges, counts;
  displayed_channels( meas, 1, edges, counts );

  // The curves end on the ROI edges, or the chart edges where the ROI is past them.
  const vector<shared_ptr<const PeakDef>> roi = make_linear_roi( 490.1, 510.7 );
  SpectrumChart::RoiCurves curves;
  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 480.0, 520.0, 25.0, curves ) );
  BOOST_CHECK_CLOSE( curves.energies.front(), 490.1, 1.0E-9 );
  BOOST_CHECK_CLOSE( curves.energies.back(), 510.7, 1.0E-9 );

  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 495.0, 505.0, 100.0, curves ) );
  BOOST_CHECK_CLOSE( curves.energies.front(), 495.0, 1.0E-9 );
  BOOST_CHECK_CLOSE( curves.energies.back(), 505.0, 1.0E-9 );

  for( size_t i = 1; i < curves.energies.size(); ++i )
    BOOST_CHECK_GT( curves.energies[i], curves.energies[i-1] );

  // Nothing to draw for a ROI out of view
  BOOST_CHECK( !SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 520.0, 600.0, 25.0, curves ) );
}//BOOST_AUTO_TEST_CASE( RoiEdges )


BOOST_AUTO_TEST_CASE( NarrowRoiInWideChannel )
{
  // Zoomed out, a ROI inside a single displayed channel is still drawn, about as tall as the bar.
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  vector<double> edges, counts;
  displayed_channels( meas, 32, edges, counts );  // 16 keV wide channels

  const vector<shared_ptr<const PeakDef>> roi = make_linear_roi( 499.0, 502.0 );
  SpectrumChart::RoiCurves curves;
  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 0.0, 1024.0, 1.0, curves ) );
  BOOST_REQUIRE( curves.energies.size() >= 2 );
  BOOST_CHECK_CLOSE( curves.energies.front(), 499.0, 1.0E-9 );
  BOOST_CHECK_CLOSE( curves.energies.back(), 502.0, 1.0E-9 );

  double maxpeak = 0.0;
  for( const double v : curves.peaks[0].counts )
    maxpeak = std::max( maxpeak, v );
  BOOST_CHECK_GT( maxpeak, 0.5*2.0E4 );
}//BOOST_AUTO_TEST_CASE( NarrowRoiInWideChannel )


BOOST_AUTO_TEST_CASE( MultiPeakExtents )
{
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  vector<double> edges, counts;
  displayed_channels( meas, 1, edges, counts );

  auto small = make_shared<PeakDef>( 495.0, 1.0, 100.0 );
  auto large = make_shared<PeakDef>( 505.0, 1.0, 1.0E6 );
  small->getContinuum()->setType( PeakContinuum::Linear );
  small->getContinuum()->setParameters( 485.0, vector<double>{ 50.0, 0.0 }, vector<double>{} );
  small->getContinuum()->setRange( 485.0, 520.0 );
  large->setContinuum( small->getContinuum() );

  const vector<shared_ptr<const PeakDef>> roi{ small, large };
  SpectrumChart::RoiCurves curves;
  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 480.0, 530.0, 20.0, curves ) );
  BOOST_REQUIRE( curves.peaks.size() == 2 );

  // In a multi-peak ROI, a Gaussian is drawn out to max(4, sqrt(2 ln(4000 amp))) sigma
  const auto check_extent = [&]( const SpectrumChart::RoiCurves::PeakCurve &pc, const PeakDef &peak ){
    const double nsigma = std::max( 4.0, std::sqrt( 2.0*std::log( 4.0E3*peak.amplitude() ) ) );
    BOOST_REQUIRE( pc.first <= pc.last );
    BOOST_CHECK_GT( pc.first, 0 );
    BOOST_CHECK_LT( pc.last, curves.energies.size() - 1 );
    BOOST_CHECK_LT( curves.energies[pc.first], peak.mean() - nsigma*peak.sigma() + 0.5 );
    BOOST_CHECK_GT( curves.energies[pc.first + 1], peak.mean() - nsigma*peak.sigma() - 1.0 );
    BOOST_CHECK_GT( curves.energies[pc.last], peak.mean() + nsigma*peak.sigma() - 0.5 );
    BOOST_CHECK_EQUAL( pc.counts[pc.first], 0.0 );  // starts and ends on what it sits on
    BOOST_CHECK_EQUAL( pc.counts[pc.last], 0.0 );
  };
  check_extent( curves.peaks[0], *small );
  check_extent( curves.peaks[1], *large );
  BOOST_CHECK_LT( curves.peaks[0].last, curves.peaks[1].last );
}//BOOST_AUTO_TEST_CASE( MultiPeakExtents )


BOOST_AUTO_TEST_CASE( BackgroundSubtract )
{
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  vector<double> edges, counts;
  displayed_channels( meas, 2, edges, counts );
  const vector<double> background( counts.size(), 30.0 );

  const vector<shared_ptr<const PeakDef>> roi = make_linear_roi( 490.0, 511.0 );
  SpectrumChart::RoiCurves curves;
  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, &background, 480.0, 520.0, 25.0, curves ) );

  const size_t nchecked = for_each_center_sample( curves, edges, [&]( size_t i, double lo, double hi ){
    const double cont = std::max( 0.0, linear_area( 200.0, -2.0, 495.0, lo, hi ) );
    BOOST_CHECK_SMALL( curves.continuum[i] - (cont - 30.0), 1.0E-6*std::max( 1.0, cont ) );
  } );
  BOOST_CHECK_GT( nchecked, size_t(10) );
}//BOOST_AUTO_TEST_CASE( BackgroundSubtract )


BOOST_AUTO_TEST_CASE( StepContinua )
{
  // Step continua are drawn from their channel-by-channel counts, so at displayed-channel centers
  //  they equal the sum of the channel counts the continuum defines.
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();

  for( const PeakContinuum::OffsetType type : { PeakContinuum::FlatStep, PeakContinuum::LinearStep,
                                                PeakContinuum::FlatStepCDF, PeakContinuum::LinearStepCDF,
                                                PeakContinuum::BiLinearStepCDF } )
  {
    auto peak = make_shared<PeakDef>( 500.3, 1.2, 2.0E4 );
    shared_ptr<PeakContinuum> cont = peak->getContinuum();
    cont->setType( type );
    vector<double> pars( cont->parameters().size(), 0.0 );
    pars[0] = 120.0;
    pars.back() = 15.0;
    cont->setParameters( 490.0, pars, vector<double>{} );
    cont->setRange( 490.0, 511.0 );
    const vector<shared_ptr<const PeakDef>> roi{ peak };

    vector<float> chan_energies;
    vector<double> chan_counts;
    BOOST_REQUIRE( PeakDef::continuum_channel_counts( roi, meas, chan_energies, chan_counts ) );
    BOOST_REQUIRE( chan_energies.size() == (chan_counts.size() + 1) );

    for( const size_t rebin : { size_t(1), size_t(4) } )
    {
      vector<double> edges, counts;
      displayed_channels( meas, rebin, edges, counts );

      SpectrumChart::RoiCurves curves;
      BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 480.0, 520.0, 25.0, curves ) );

      size_t nchecked = 0;
      for_each_center_sample( curves, edges, [&]( size_t i, double lo, double hi ){
        if( (lo < chan_energies.front()) || (hi > chan_energies.back()) )
          return;
        double expected = 0.0;
        for( size_t j = 0; j < chan_counts.size(); ++j )
        {
          if( (chan_energies[j] >= (lo - 1.0E-4)) && (chan_energies[j+1] <= (hi + 1.0E-4)) )
            expected += chan_counts[j];
        }
        BOOST_CHECK_SMALL( curves.continuum[i] - std::max( 0.0, expected ), 1.0E-4*std::max( 1.0, expected ) );
        ++nchecked;
      } );
      BOOST_CHECK_GT( nchecked, size_t(4) );
    }//for( rebin of 1 and 4 )
  }//for( loop over step continuum types )
}//BOOST_AUTO_TEST_CASE( StepContinua )


BOOST_AUTO_TEST_CASE( AllSkewTypes )
{
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  vector<double> edges, counts;
  displayed_channels( meas, 1, edges, counts );

  for( int i = 0; i < PeakDef::NumSkewType; ++i )
  {
    const auto type = PeakDef::SkewType( i );
    BOOST_TEST_CONTEXT( "Skew type " << PeakDef::to_string(type) )
    {
      auto peak = make_shared<PeakDef>( 500.3, 1.2, 2.0E4 );
      set_reasonable_skew( *peak, type );
      peak->getContinuum()->setType( PeakContinuum::Linear );
      peak->getContinuum()->setParameters( 495.0, vector<double>{ 200.0, -2.0 }, vector<double>{} );
      peak->getContinuum()->setRange( 480.0, 520.0 );

      auto other = make_shared<PeakDef>( 510.0, 1.2, 1.0E3 );
      set_reasonable_skew( *other, type );
      other->setContinuum( peak->getContinuum() );

      const vector<shared_ptr<const PeakDef>> roi{ peak, other };
      SpectrumChart::RoiCurves curves;
      BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 470.0, 530.0, 20.0, curves ) );

      double area = 0.0;
      for( const SpectrumChart::RoiCurves::PeakCurve &pc : curves.peaks )
      {
        BOOST_REQUIRE( pc.first <= pc.last );
        for( const double v : pc.counts )
        {
          BOOST_CHECK( std::isfinite( v ) && (v >= 0.0) );
          area = std::max( area, v );
        }
      }
      BOOST_CHECK_GT( area, 1000.0 );  // the 2E4 peak in 0.5 keV channels peaks well above this
    }//BOOST_TEST_CONTEXT
  }//for( loop over skew types )
}//BOOST_AUTO_TEST_CASE( AllSkewTypes )


BOOST_AUTO_TEST_CASE( DataDefinedPeak )
{
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  vector<double> edges, counts;
  displayed_channels( meas, 1, edges, counts );

  auto peak = make_shared<PeakDef>( 494.0, 507.0, 500.3, meas, nullptr );
  BOOST_REQUIRE( !peak->gausPeak() );

  const vector<shared_ptr<const PeakDef>> roi{ peak };
  SpectrumChart::RoiCurves curves;
  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 480.0, 520.0, 25.0, curves ) );
  BOOST_REQUIRE( curves.peaks.size() == 1 );

  double maxval = 0.0;
  for( const double v : curves.peaks[0].counts )
  {
    BOOST_CHECK( std::isfinite( v ) && (v >= 0.0) );
    maxval = std::max( maxval, v );
  }
  BOOST_CHECK_GT( maxval, 1000.0 );
}//BOOST_AUTO_TEST_CASE( DataDefinedPeak )


BOOST_AUTO_TEST_CASE( DataDefinedPeakBackgroundSubtracted )
{
  // A data-defined peak is the data above the continuum, so subtracting the background from both
  //  leaves it unchanged.
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();
  vector<double> edges, counts;
  displayed_channels( meas, 2, edges, counts );
  const vector<double> background( counts.size(), 30.0 );

  auto peak = make_shared<PeakDef>( 494.0, 507.0, 500.3, meas, nullptr );
  const vector<shared_ptr<const PeakDef>> roi{ peak };

  SpectrumChart::RoiCurves plain, subtracted;
  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 480.0, 520.0, 25.0, plain ) );
  BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, &background, 480.0, 520.0, 25.0, subtracted ) );
  BOOST_REQUIRE( plain.energies.size() == subtracted.energies.size() );

  double maxval = 0.0;
  for( size_t i = 0; i < plain.energies.size(); ++i )
  {
    BOOST_CHECK_SMALL( plain.peaks[0].counts[i] - subtracted.peaks[0].counts[i], 1.0E-9*std::max( 1.0, plain.peaks[0].counts[i] ) );
    BOOST_CHECK_SMALL( (plain.continuum[i] - 30.0) - subtracted.continuum[i], 1.0E-9*std::max( 1.0, plain.continuum[i] ) );
    maxval = std::max( maxval, subtracted.peaks[0].counts[i] );
  }
  BOOST_CHECK_GT( maxval, 1000.0 );
}//BOOST_AUTO_TEST_CASE( DataDefinedPeakBackgroundSubtracted )


BOOST_AUTO_TEST_CASE( ExternalContinuum )
{
  // An External continuum is drawn from its histograms channels, like the D3 chart does.
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();

  auto external = make_shared<SpecUtils::Measurement>( *meas );
  auto ext_counts = make_shared<vector<float>>( sm_num_channels );
  for( size_t i = 0; i < sm_num_channels; ++i )
    (*ext_counts)[i] = static_cast<float>( 20.0 + 0.05*i );
  external->set_gamma_counts( ext_counts, 300.0f, 310.0f );

  auto peak = make_shared<PeakDef>( 500.3, 1.2, 2.0E4 );
  shared_ptr<PeakContinuum> cont = peak->getContinuum();
  cont->setType( PeakContinuum::External );
  cont->setExternalContinuum( external );
  cont->setRange( 490.0, 511.0 );
  const vector<shared_ptr<const PeakDef>> roi{ peak };

  for( const size_t rebin : { size_t(1), size_t(2) } )
  {
    vector<double> edges, counts;
    displayed_channels( meas, rebin, edges, counts );

    SpectrumChart::RoiCurves curves;
    BOOST_REQUIRE( SpectrumChart::roiCurves( roi, meas, edges, counts, nullptr, 480.0, 520.0, 25.0, curves ) );

    const size_t nchecked = for_each_center_sample( curves, edges, [&]( size_t i, double lo, double hi ){
      const double expected = external->gamma_integral( static_cast<float>(lo), static_cast<float>(hi) );
      BOOST_CHECK_SMALL( curves.continuum[i] - expected, 1.0E-4*std::max( 1.0, expected ) );
    } );
    BOOST_CHECK_GT( nchecked, size_t(8) );
  }//for( rebin of 1 and 2 )
}//BOOST_AUTO_TEST_CASE( ExternalContinuum )


#if( SpecUtils_ENABLE_D3_CHART )
BOOST_AUTO_TEST_CASE( PeakJsonForD3Chart )
{
  // `PeakDef::gaus_peaks_to_json(...)` gives the D3 chart the channel counts of step and External
  //  continua (the ROI's channels, plus two each side), and the drawing range of skewed peaks;
  //  check these against independent computations.
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();

  // The channels the JSON should cover, and the counts in each
  const auto check_channels = []( const nlohmann::json &roi, const shared_ptr<const SpecUtils::Measurement> &hist,
                                   const double lower, const double upper,
                                   const std::function<double(size_t)> &expected_counts ){
    BOOST_REQUIRE( roi.contains( "continuumEnergies" ) && roi.contains( "continuumCounts" ) );
    const size_t nchannel = hist->num_gamma_channels();
    const size_t first = hist->find_gamma_channel( static_cast<float>(lower) ) - 2;
    const size_t last = std::min( hist->find_gamma_channel( static_cast<float>(upper) ) + 2, nchannel - 1 );

    const nlohmann::json &energies = roi["continuumEnergies"];
    const nlohmann::json &counts = roi["continuumCounts"];
    BOOST_REQUIRE_EQUAL( counts.size(), last - first + 1 );
    BOOST_REQUIRE_EQUAL( energies.size(), counts.size() + 1 );
    for( size_t i = 0; i < counts.size(); ++i )
    {
      BOOST_CHECK_SMALL( energies[i].get<double>() - hist->gamma_channel_lower( first + i ), 1.0E-3 );
      const double expected = expected_counts( first + i );
      BOOST_CHECK_SMALL( counts[i].get<double>() - expected, 1.0E-5*std::max( 1.0, std::fabs(expected) ) );
    }
    BOOST_CHECK_SMALL( energies.back().get<double>() - hist->gamma_channel_upper( last ), 1.0E-3 );
  };//check_channels

  // Step continua: counts from the continuum, channel by channel
  for( const PeakContinuum::OffsetType type : { PeakContinuum::FlatStep, PeakContinuum::LinearStepCDF } )
  {
    auto peak = make_shared<PeakDef>( 500.3, 1.2, 2.0E4 );
    shared_ptr<PeakContinuum> cont = peak->getContinuum();
    cont->setType( type );
    vector<double> pars( cont->parameters().size(), 0.0 );
    pars[0] = 120.0;
    pars.back() = 15.0;
    cont->setParameters( 490.0, pars, vector<double>{} );
    cont->setRange( 490.0, 511.0 );
    const vector<shared_ptr<const PeakDef>> roi_peaks{ peak };

    const nlohmann::json roi = nlohmann::json::parse(
                     PeakDef::gaus_peaks_to_json( roi_peaks, meas, Wt::WColor( 0, 51, 255 ), 255 ) );
    check_channels( roi, meas, 490.0, 511.0, [&]( size_t channel ){
      return cont->offset_integral( meas->gamma_channel_lower(channel), meas->gamma_channel_upper(channel),
                                    meas, roi_peaks );
    } );
    BOOST_CHECK( !roi["peaks"][0].contains( "visRange" ) );  // only for skewed peaks
  }//for( step continuum types )

  // External continuum: the external histograms channel contents
  {
    auto external = make_shared<SpecUtils::Measurement>( *meas );
    auto ext_counts = make_shared<vector<float>>( sm_num_channels );
    for( size_t i = 0; i < sm_num_channels; ++i )
      (*ext_counts)[i] = static_cast<float>( 20.0 + 0.05*i );
    external->set_gamma_counts( ext_counts, 300.0f, 310.0f );

    auto peak = make_shared<PeakDef>( 500.3, 1.2, 2.0E4 );
    set_reasonable_skew( *peak, PeakDef::GaussExp );
    peak->getContinuum()->setType( PeakContinuum::External );
    peak->getContinuum()->setExternalContinuum( external );
    peak->getContinuum()->setRange( 490.0, 511.0 );

    const nlohmann::json roi = nlohmann::json::parse(
                 PeakDef::gaus_peaks_to_json( { peak }, meas, Wt::WColor( 0, 51, 255 ), 255 ) );
    check_channels( roi, external, 490.0, 511.0, [&]( size_t channel ){
      return external->gamma_channel_content( channel );
    } );

    // The skewed peak's drawing range
    const optional<pair<double,double>> range = peak->visRange();
    BOOST_REQUIRE( range.has_value() );
    BOOST_REQUIRE( roi["peaks"][0].contains( "visRange" ) );
    BOOST_CHECK_SMALL( roi["peaks"][0]["visRange"][0].get<double>() - range->first, 1.0E-3 );
    BOOST_CHECK_SMALL( roi["peaks"][0]["visRange"][1].get<double>() - range->second, 1.0E-3 );
    BOOST_CHECK_LT( range->first, peak->mean() - 3.0*peak->sigma() );  // a low-side tail
    BOOST_CHECK_GT( range->second, peak->mean() + 3.0*peak->sigma() );
  }
}//BOOST_AUTO_TEST_CASE( PeakJsonForD3Chart )
#endif //SpecUtils_ENABLE_D3_CHART


BOOST_AUTO_TEST_CASE( RenderCharts )
{
  // Renders whole charts, with every kind of continuum, to make sure nothing throws, and no
  //  non-finite coordinates make it into the output.
  const shared_ptr<SpecUtils::Measurement> meas = make_spectrum();

  auto background = make_shared<SpecUtils::Measurement>( *meas );
  {
    auto bcounts = make_shared<vector<float>>( sm_num_channels, 90.0f );
    background->set_gamma_counts( bcounts, 300.0f, 310.0f );
    background->set_title( "Background" );
  }

  auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>();
  {
    // A lone peak on a linear continuum
    peaks->push_back( make_linear_roi( 490.0, 511.0 )[0] );

    // Several peaks sharing a stepped continuum, one with its own color
    double mean = 100.0;
    for( const PeakContinuum::OffsetType type : { PeakContinuum::FlatStep, PeakContinuum::LinearStepCDF,
                                                  PeakContinuum::Quadratic, PeakContinuum::NoOffset } )
    {
      auto first = make_shared<PeakDef>( mean, 1.5, 3.0E3 );
      auto second = make_shared<PeakDef>( mean + 6.0, 1.5, 1.0E3 );
      set_reasonable_skew( *second, PeakDef::GaussExp );
      second->setLineColor( Wt::WColor( 200, 30, 30 ) );
      shared_ptr<PeakContinuum> cont = first->getContinuum();
      cont->setType( type );
      if( type != PeakContinuum::NoOffset )
      {
        vector<double> pars( cont->parameters().size(), 0.0 );
        pars[0] = 50.0;
        pars.back() = (pars.size() > 1) ? 5.0 : pars[0];
        cont->setParameters( mean - 10.0, pars, vector<double>{} );
      }
      cont->setRange( mean - 10.0, mean + 16.0 );
      second->setContinuum( cont );
      peaks->push_back( first );
      peaks->push_back( second );
      mean += 100.0;
    }

    // A data-defined peak
    peaks->push_back( make_shared<PeakDef>( 694.0, 706.0, 700.0, meas, nullptr ) );
  }

  auto ref_lines = make_shared<ReferenceLineInfo>();
  {
    ref_lines->m_validity = ReferenceLineInfo::InputValidity::Valid;
    ref_lines->m_input.m_input_txt = "Test lines";
    ref_lines->m_input.m_color = Wt::WColor( 0, 128, 0 );
    const ReferenceLineInfo::RefLine::Particle particles[] = { ReferenceLineInfo::RefLine::Particle::Gamma,
      ReferenceLineInfo::RefLine::Particle::Xray, ReferenceLineInfo::RefLine::Particle::Beta,
      ReferenceLineInfo::RefLine::Particle::Alpha };
    for( int i = 0; i < 12; ++i )
    {
      ReferenceLineInfo::RefLine line;
      line.m_energy = 150.0 + 60.0*i;
      line.m_normalized_intensity = 1.0 / (1 + i);
      line.m_uncorrected_intensity = line.m_normalized_intensity * ((i % 3) ? 1.0 : 1.3);
      line.m_particle_type = particles[i % 4];
      line.m_major_line = (i < 3);
      ref_lines->m_ref_lines.push_back( line );
    }
  }

  const auto check_svg = []( const std::string &svg ){
    BOOST_CHECK_GT( svg.size(), size_t(1000) );
    // (Wt writes non-finite numbers as "NaN", "Infinity", or "-Infinity")
    BOOST_CHECK( !std::regex_search( svg, std::regex( "\\b(nan|inf|infinity)\\b", std::regex::icase ) ) );
  };

  SpectrumChartImage::SpectrumOptions options;
  options.foreground = meas;
  options.peaks = peaks;
  options.reference_lines.push_back( ref_lines );
  options.title = "Spectrum chart test";

  check_svg( render_svg( options, 1200, 600, "full_log" ) );

  options.log_y = false;
  options.energy_range = std::make_pair( 470.0, 530.0 );
  check_svg( render_svg( options, 1200, 600, "zoomed_linear" ) );

  options.log_y = true;
  options.energy_range = std::make_pair( 80.0, 520.0 );
  options.background = background;
  check_svg( render_svg( options, 1000, 500, "with_background" ) );

  options.background_subtract = true;
  options.log_y = false;
  check_svg( render_svg( options, 1000, 500, "background_subtracted" ) );

  options.background_subtract = false;
  options.log_y = true;
  options.energy_range.reset();
  check_svg( render_svg( options, 360, 220, "compact" ) );

  options.thumbnail = true;
  check_svg( render_svg( options, 160, 120, "thumbnail" ) );
}//BOOST_AUTO_TEST_CASE( RenderCharts )
