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

#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <boost/asio/io_service.hpp>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WServer.h>
#include <Wt/WIOService.h>
#include <Wt/WGridLayout.h>
#include <Wt/WComboBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WGroupBox.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WProgressBar.h>
#include <Wt/WContainerWidget.h>

// CeeLo (external_libs/CeeLo/src)
#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfChart.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/MakeMcResponseForDrf.h"
#include "InterSpec/DetectorGeometryInput.h"

using namespace Wt;
using namespace std;


vector<ceelo::GroundingPoint> MakeMcResponseForDrf::groundingPointsForDrf(
                            const shared_ptr<const DetectorPeakResponse> &drf,
                            const ceelo::GeometryDescriptor &geom,
                            bool &curve_derived )
{
    vector<ceelo::GroundingPoint> answer;
    curve_derived = false;

    if( !drf || !drf->isValid() )
      return answer;

    const shared_ptr<const MeasuredDrfPoints> raw = drf->measuredPoints();
    if( raw && !raw->empty() )
    {
      for( const MeasuredEffPoint &p : raw->points() )
      {
        if( (p.distance < 0.0f) || (p.efficiency <= 0.0f) )
          continue;  //fixed-geometry / invalid points cant anchor a geometry model

        ceelo::GroundingPoint gp;
        gp.energy_keV = p.energy;
        gp.measured_eff = p.efficiency;
        gp.frac_stat_sigma = p.fracStatUncert;
        gp.frac_cert_sigma = p.fracCertUncert;
        gp.source_key = p.sourceKey;
        gp.distance_cm = p.distance / PhysicalUnits::cm;
        gp.cos_theta = 1.0;  //Make Detector Response sources are on-axis
        answer.push_back( std::move(gp) );
      }//for( const MeasuredEffPoint &p : raw->points() )

      if( !answer.empty() )
        return answer;
    }//if( raw && !raw->empty() )

    if( drf->isFixedGeometry() )
      return answer;  //no geometry to map a fixed-geometry curve through

    // Curve fallback: sample the legacy intrinsic curve and turn it into
    //  absolute efficiencies.  When the curve was specified as absolute
    //  efficiency at a stated distance, anchor there - that is where the curve
    //  is actually pinned to data; otherwise use a comfortably far-field
    //  reference distance.  (intrinsicEfficiency() already backs out any air
    //  attenuation the absolute curve included, so the reconstructed absolute
    //  efficiencies below are in-vacuum - consistent with the ray-trace kernel
    //  the grounding fit compares against.)
    curve_derived = true;

    const double a_cm = geom.transverse_half_extent();
    double d_ref_cm = std::max( 50.0, 10.0 * a_cm );
    if( (drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute)
        && (drf->absoluteEfficiencyDistance() > 0.0) )
      d_ref_cm = drf->absoluteEfficiencyDistance() / PhysicalUnits::cm;
    const double diam = drf->detectorDiameter();
    const double dist = d_ref_cm * PhysicalUnits::cm;
    const double frac_solid_angle = DetectorPeakResponse::fractionalSolidAngle(
                                        diam, dist + drf->detectorSetback() );

    double e_lo = drf->lowerEnergy(), e_hi = drf->upperEnergy();
    if( (e_lo <= 0.0) || (e_hi <= e_lo) )
    {
      e_lo = 59.0;
      e_hi = 2614.0;
    }
    e_lo = std::max( e_lo, 20.0 );

    const int n_samples = 8;
    vector<double> energies;
    for( int i = 0; i < n_samples; ++i )
      energies.push_back( e_lo * std::pow( e_hi/e_lo, double(i)/(n_samples-1) ) );

    const vector<double> sigmas = drf->efficiencyFracCovariance( energies ).empty()
                          ? vector<double>()
                          : [&](){
                              vector<double> s;
                              const vector<double> cov = drf->efficiencyFracCovariance( energies );
                              for( size_t i = 0; i < energies.size(); ++i )
                                s.push_back( std::sqrt( std::max(0.0, cov[i*energies.size()+i]) ) );
                              return s;
                            }();

    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double intrinsic = drf->intrinsicEfficiency( static_cast<float>(energies[i]) );
      if( intrinsic <= 0.0 )
        continue;

      ceelo::GroundingPoint gp;
      gp.energy_keV = energies[i];
      gp.measured_eff = intrinsic * frac_solid_angle;
      gp.frac_stat_sigma = sigmas.empty() ? 0.05 : std::max( 0.01, sigmas[i] );
      gp.frac_cert_sigma = 0.0;
      gp.source_key = "legacy-curve";
      gp.distance_cm = d_ref_cm;
      gp.cos_theta = 1.0;
      answer.push_back( std::move(gp) );
    }//for( size_t i = 0; i < energies.size(); ++i )

    return answer;
}//groundingPointsForDrf(...)


MakeMcResponseForDrf::MakeMcResponseForDrf( InterSpec *viewer,
                            std::shared_ptr<const DetectorPeakResponse> seed_drf )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_seedDrf( seed_drf ),
    m_geometry( nullptr ),
    m_method( nullptr ),
    m_profile( nullptr ),
    m_precision( nullptr ),
    m_customPrecision( nullptr ),
    m_anchorAngles( nullptr ),
    m_anchorInfo( nullptr ),
    m_refDistance( nullptr ),
    m_estimate( nullptr ),
    m_profileRow( nullptr ),
    m_precRow( nullptr ),
    m_anchorAnglesRow( nullptr ),
    m_anchorRow( nullptr ),
    m_generate( nullptr ),
    m_cancelBtn( nullptr ),
    m_progress( nullptr ),
    m_status( nullptr ),
    m_chartBox( nullptr ),
    m_chart( nullptr ),
    m_chartMode( nullptr ),
    m_chartDistance( nullptr ),
    m_generationId( 0 ),
    m_cancelFlag( nullptr ),
    m_result( nullptr ),
    m_validationChanged(),
    m_responseGenerated(),
    m_updatedDrf()
{
  assert( m_interspec );
  wApp->useStyleSheet( "InterSpec_resources/MakeMcResponseForDrf.css" );
  m_interspec->useMessageResourceBundle( "MakeMcResponseForDrf" );

  addStyleClass( "MakeMcResponseForDrf" );

  //Geometry
  WGroupBox *geomBox = addNew<WGroupBox>( WString::tr("mmr-geometry-title") );
  geomBox->addStyleClass( "McGeomBox" );
  m_geometry = geomBox->addNew<DetectorGeometryInput>( m_interspec );
  m_geometry->seedFromDrf( m_seedDrf );
  m_geometry->changed().connect( this, &MakeMcResponseForDrf::handleGeometryChanged );

  //Characterization options
  WGroupBox *optsBox = addNew<WGroupBox>( WString::tr("mmr-options-title") );
  optsBox->addStyleClass( "McOptsBox" );

  WContainerWidget *methodRow = optsBox->addNew<WContainerWidget>();
  methodRow->addNew<WLabel>( WString::tr("mmr-method") );
  m_method = methodRow->addNew<WComboBox>();
  m_method->addItem( WString::tr("mmr-method-full-mc") );        //Method::FullMc
  m_method->addItem( WString::tr("mmr-method-quick-mc") );       //Method::QuickMc
  m_method->addItem( WString::tr("mmr-method-curve-transfer") ); //Method::CurveTransfer
  m_method->setCurrentIndex( 0 );
  m_method->activated().connect( this, &MakeMcResponseForDrf::handleMethodChanged );
  HelpSystem::attachToolTipOn( m_method, WString::tr("mmr-tt-method"), true );

  // When the DRF already carries a quick/transfer response, nudge toward the
  //  full characterization as the accuracy upgrade.
  const shared_ptr<const ceelo::DetectorResponse> seed_resp
                             = m_seedDrf ? m_seedDrf->ceeloResponse() : nullptr;
  if( seed_resp && seed_resp->model_transfer.has_value() )
  {
    WText *upgradeNote = optsBox->addNew<WText>( WString::tr("mmr-upgrade-note") );
    upgradeNote->addStyleClass( "McUpgradeNote" );
    upgradeNote->setInline( false );
  }

  m_profileRow = optsBox->addNew<WContainerWidget>();
  m_profileRow->addNew<WLabel>( WString::tr("mmr-profile") );
  m_profile = m_profileRow->addNew<WComboBox>();
  m_profile->addItem( WString::tr("mmr-profile-general") );
  m_profile->addItem( WString::tr("mmr-profile-farfield") );
  m_profile->addItem( WString::tr("mmr-profile-contact") );
  m_profile->setCurrentIndex( 0 );
  m_profile->activated().connect( this, &MakeMcResponseForDrf::updateEstimate );
  HelpSystem::attachToolTipOn( m_profile, WString::tr("mmr-tt-profile"), true );

  m_precRow = optsBox->addNew<WContainerWidget>();
  m_precRow->addNew<WLabel>( WString::tr("mmr-precision") );
  m_precision = m_precRow->addNew<WComboBox>();
  m_precision->addItem( WString::tr("mmr-precision-fast") );      //idx0: 1%
  m_precision->addItem( WString::tr("mmr-precision-normal") );    //idx1: 0.3%
  m_precision->addItem( WString::tr("mmr-precision-balanced") );  //idx2: relax_mild (0.3% base, relaxed high-E)
  m_precision->addItem( WString::tr("mmr-precision-thorough") );  //idx3: 0.1%
  m_precision->addItem( WString::tr("mmr-precision-custom") );    //idx4
  m_precision->setCurrentIndex( 1 );
  m_precision->activated().connect( this, &MakeMcResponseForDrf::handlePrecisionChanged );
  HelpSystem::attachToolTipOn( m_precision, WString::tr("mmr-tt-precision"), true );

  m_customPrecision = m_precRow->addNew<WLineEdit>();
  m_customPrecision->setTextSize( 5 );
  m_customPrecision->setPlaceholderText( "0.3%" );
  m_customPrecision->changed().connect( this, &MakeMcResponseForDrf::updateEstimate );
  m_customPrecision->hide();

  m_anchorAnglesRow = optsBox->addNew<WContainerWidget>();
  m_anchorAnglesRow->addNew<WLabel>( WString::tr("mmr-anchor-angles") );
  m_anchorAngles = m_anchorAnglesRow->addNew<WComboBox>();
  m_anchorAngles->addItem( WString::tr("mmr-anchor-angles-1") );  //on-axis only
  m_anchorAngles->addItem( WString::tr("mmr-anchor-angles-3") );  //3 cos-theta anchors
  m_anchorAngles->setCurrentIndex( 1 );  //3 anchors: off-axis FEP + real HPGe total
  m_anchorAngles->activated().connect( this, &MakeMcResponseForDrf::updateEstimate );
  HelpSystem::attachToolTipOn( m_anchorAngles, WString::tr("mmr-tt-anchor-angles"), true );
  m_anchorAnglesRow->hide();

  m_anchorRow = optsBox->addNew<WContainerWidget>();
  m_anchorInfo = m_anchorRow->addNew<WText>( "" );
  m_anchorInfo->addStyleClass( "McAnchorInfo" );
  m_anchorInfo->setInline( false );
  WContainerWidget *refDistRow = m_anchorRow->addNew<WContainerWidget>();
  refDistRow->addNew<WLabel>( WString::tr("mmr-ref-distance") );
  m_refDistance = refDistRow->addNew<WLineEdit>();
  m_refDistance->setTextSize( 10 );
  m_refDistance->changed().connect( this, &MakeMcResponseForDrf::handleGeometryChanged );
  HelpSystem::attachToolTipOn( m_refDistance, WString::tr("mmr-tt-ref-distance"), true );
  m_anchorRow->hide();

  m_estimate = optsBox->addNew<WText>( "" );
  m_estimate->addStyleClass( "McEstimate" );
  m_estimate->setInline( false );

  //Run row
  WContainerWidget *runRow = addNew<WContainerWidget>();
  runRow->addStyleClass( "McRunRow" );
  m_generate = runRow->addNew<WPushButton>( WString::tr("mmr-generate-btn") );
  m_generate->clicked().connect( this, &MakeMcResponseForDrf::startGeneration );
  m_cancelBtn = runRow->addNew<WPushButton>( WString::tr("Cancel") );
  m_cancelBtn->clicked().connect( this, &MakeMcResponseForDrf::cancelGeneration );
  m_cancelBtn->hide();
  m_progress = runRow->addNew<WProgressBar>();
  m_progress->setRange( 0.0, 1.0 );
  m_progress->hide();

  m_status = addNew<WText>( "" );
  m_status->addStyleClass( "McStatus" );
  m_status->setInline( false );

  // Response preview: per-angle efficiency curves for the generated response.
  //  Hidden until a response exists.
  m_chartBox = addNew<WGroupBox>( WString::tr("mmr-chart-title") );
  m_chartBox->addStyleClass( "McChartBox" );

  WContainerWidget *chartCtrlRow = m_chartBox->addNew<WContainerWidget>();
  chartCtrlRow->addStyleClass( "McChartCtrlRow" );
  chartCtrlRow->addNew<WLabel>( WString::tr("mmr-chart-mode") );
  m_chartMode = chartCtrlRow->addNew<WComboBox>();
  m_chartMode->addItem( WString::tr("mmr-chart-mode-absolute") );   //idx0
  m_chartMode->addItem( WString::tr("mmr-chart-mode-intrinsic") );  //idx1
  m_chartMode->setCurrentIndex( 0 );
  m_chartMode->activated().connect( this, &MakeMcResponseForDrf::updateResponseChart );

  WLabel *distLabel = chartCtrlRow->addNew<WLabel>( WString::tr("mmr-chart-distance") );
  m_chartDistance = chartCtrlRow->addNew<WLineEdit>( "25 cm" );
  m_chartDistance->setTextSize( 8 );
  distLabel->setBuddy( m_chartDistance );
  m_chartDistance->changed().connect( this, &MakeMcResponseForDrf::updateResponseChart );
  m_chartDistance->enterPressed().connect( this, &MakeMcResponseForDrf::updateResponseChart );
  HelpSystem::attachToolTipOn( m_chartDistance, WString::tr("mmr-tt-chart-distance"), true );

  m_chart = m_chartBox->addNew<DrfChart>();
  m_chart->addStyleClass( "McChart" );

  m_chartBox->hide();

  handleMethodChanged();
}//MakeMcResponseForDrf constructor


MakeMcResponseForDrf::~MakeMcResponseForDrf()
{
  // The worker only holds the cancel flag and plain strings - never `this` -
  //  so all we need to do is ask any in-flight generation to stop.
  if( m_cancelFlag )
    m_cancelFlag->store( true );
}//~MakeMcResponseForDrf()


Wt::Signal<bool> &MakeMcResponseForDrf::validationChanged()
{
  return m_validationChanged;
}


bool MakeMcResponseForDrf::hasResult() const
{
  return !!m_result;
}


std::shared_ptr<const ceelo::DetectorResponse> MakeMcResponseForDrf::generatedResponse() const
{
  return m_result;
}


Wt::Signal<std::shared_ptr<ceelo::DetectorResponse>> &MakeMcResponseForDrf::responseGenerated()
{
  return m_responseGenerated;
}


Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &MakeMcResponseForDrf::updatedDrf()
{
  return m_updatedDrf;
}


MakeMcResponseForDrf::Method MakeMcResponseForDrf::selectedMethod() const
{
  switch( m_method->currentIndex() )
  {
    case 1: return Method::QuickMc;
    case 2: return Method::CurveTransfer;
    default: return Method::FullMc;
  }
}//selectedMethod()


void MakeMcResponseForDrf::handleMethodChanged()
{
  const Method method = selectedMethod();

  m_profileRow->setHidden( method != Method::FullMc );
  m_precRow->setHidden( method == Method::CurveTransfer );
  m_anchorAnglesRow->setHidden( method != Method::QuickMc );
  m_anchorRow->setHidden( method != Method::CurveTransfer );

  // The measured-curve transfer is instant and rebuilds automatically - no
  //  explicit "Generate" step (this is what guarantees that entering geometry
  //  and accepting the dialog always yields an attached, distance-aware
  //  response).
  m_generate->setHidden( method == Method::CurveTransfer );

  // A result from a different method is not what the user is configuring;
  //  abandon any in-flight generation too (its finish handler is stale-guarded
  //  and balances the update lock itself), and take back the run-row UI.
  ++m_generationId;
  if( m_cancelFlag )
    m_cancelFlag->store( true );
  m_progress->hide();
  m_cancelBtn->hide();
  if( m_result )
  {
    m_result.reset();
    m_validationChanged.emit( false );
  }
  m_status->setText( "" );

  handleGeometryChanged();
}//handleMethodChanged()


void MakeMcResponseForDrf::handleGeometryChanged()
{
  //Any geometry change invalidates a previously generated response.
  if( m_result )
  {
    m_result.reset();
    m_validationChanged.emit( false );
    m_status->setText( WString::tr("mmr-status-stale") );
    updateResponseChart();
  }

  m_generate->setEnabled( m_geometry->isValid() );
  updateAnchorInfo();
  updateEstimate();

  // Auto-build the instant transfer whenever the inputs are usable.
  if( (selectedMethod() == Method::CurveTransfer) && m_geometry->isValid() )
    startGeneration();
}//handleGeometryChanged()


void MakeMcResponseForDrf::updateAnchorInfo()
{
  if( selectedMethod() != Method::CurveTransfer )
    return;

  if( !m_geometry->isValid() )
  {
    m_anchorInfo->setText( "" );
    return;
  }

  try
  {
    const ceelo::GeometryDescriptor gd = m_geometry->toDescriptor();
    const CeeLoUtils::TransferAnchor anchor
          = CeeLoUtils::transferAnchorForDrf( m_seedDrf, gd, refDistanceOverrideCm() );

    const string dist_str = PhysicalUnits::printToBestLengthUnits(
                                    anchor.ref_distance_cm * PhysicalUnits::cm );
    const int npoints = static_cast<int>( anchor.curve.energies_keV.size() );

    if( anchor.curve_derived )
    {
      // The reference distance is only meaningful input for a sampled curve.
      m_refDistance->setEnabled( true );
      if( m_refDistance->text().empty() )
        m_refDistance->setText( WString::fromUTF8(dist_str) );
      m_anchorInfo->setText( WString::tr("mmr-anchor-src-curve")
                              .arg( npoints ).arg( dist_str ) );
    }else
    {
      // Raw measured points pin the distance - display it, read-only.
      m_refDistance->setEnabled( false );
      m_refDistance->setText( WString::fromUTF8(dist_str) );
      m_anchorInfo->setText( WString::tr("mmr-anchor-src-points")
                              .arg( npoints ).arg( dist_str ) );
    }
  }catch( std::exception &e )
  {
    m_refDistance->setEnabled( false );
    m_anchorInfo->setText( WString::tr("mmr-anchor-unusable")
                            .arg( WString::fromUTF8(e.what()) ) );
  }
}//updateAnchorInfo()


double MakeMcResponseForDrf::refDistanceOverrideCm() const
{
  string txt = m_refDistance->text().toUTF8();
  SpecUtils::trim( txt );
  if( txt.empty() )
    return -1.0;

  try
  {
    const double dist = PhysicalUnits::stringToDistance( txt );
    return (dist > 0.0) ? (dist / PhysicalUnits::cm) : -1.0;
  }catch( std::exception & )
  {
    return -1.0;
  }
}//refDistanceOverrideCm()


void MakeMcResponseForDrf::updateResponseChart()
{
  if( !m_chart || !m_chartBox )
    return;

  if( !m_result )
  {
    m_chartBox->hide();
    m_chart->updateChart( nullptr );
    return;
  }

  // Wrap the generated response in a (throw-away) DetectorPeakResponse so the
  //  chart can sample it through the normal per-angle API.  Clone the seed/
  //  foreground DRF when available (keeps diameter / energy range), else a
  //  minimal far-field shell.
  shared_ptr<SpecMeas> foreground = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  shared_ptr<const DetectorPeakResponse> base = m_seedDrf
                                        ? m_seedDrf
                                        : (foreground ? foreground->detector() : nullptr);
  shared_ptr<DetectorPeakResponse> preview;
  if( base )
  {
    preview = make_shared<DetectorPeakResponse>( *base );
  }else
  {
    const double a_cm = m_result->transverse_half_extent();
    preview = make_shared<DetectorPeakResponse>( "preview", "" );
    preview->setIntrinsicEfficiencyFormula( "1.0", 2.0*a_cm*PhysicalUnits::cm,
                    PhysicalUnits::keV, 0.0f, 0.0f,
                    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  }
  preview->setCeeloResponse( m_result );

  // Chart source distance (absolute mode); default to 25 cm when unparseable.
  double dist = 25.0 * PhysicalUnits::cm;
  try
  {
    string txt = m_chartDistance->text().toUTF8();
    SpecUtils::trim( txt );
    if( !txt.empty() )
    {
      const double d = PhysicalUnits::stringToDistance( txt );
      if( d > 0.0 )
        dist = d;
    }
  }catch( std::exception & )
  {
  }

  m_chart->setShowResponseAngles( true );
  m_chart->setSourceDistance( dist );
  m_chart->setIntrinsicEfficiency( m_chartMode->currentIndex() == 1 );
  m_chart->updateChart( preview );
  m_chartBox->show();
}//updateResponseChart()


void MakeMcResponseForDrf::handlePrecisionChanged()
{
  m_customPrecision->setHidden( m_precision->currentIndex() != 4 );
  updateEstimate();
}//handlePrecisionChanged()


bool MakeMcResponseForDrf::selectedRelaxMild() const
{
  return (m_precision->currentIndex() == 2);
}//selectedRelaxMild()


double MakeMcResponseForDrf::selectedPrecision() const
{
  switch( m_precision->currentIndex() )
  {
    case 0: return 0.01;
    case 1: return 0.003;
    case 2: return 0.003;  //Balanced: relax_mild base precision (0.3%)
    case 3: return 0.001;
    default:
      break;
  }

  //Custom: accept "0.5%", "0.005", etc.
  string txt = m_customPrecision->text().toUTF8();
  SpecUtils::trim( txt );
  const bool percent = !txt.empty() && (txt.back() == '%');
  if( percent )
    txt.pop_back();

  try
  {
    double value = std::stod( txt );
    if( percent || (value >= 0.05) )  //"0.3" almost certainly means percent
      value /= 100.0;
    if( (value < 1.0E-4) || (value > 0.2) )
      throw runtime_error( "" );
    return value;
  }catch( std::exception & )
  {
    return 0.003;
  }
}//selectedPrecision()


void MakeMcResponseForDrf::updateEstimate()
{
  const Method method = selectedMethod();

  if( method == Method::CurveTransfer )
  {
    m_estimate->setText( WString::tr("mmr-estimate-instant") );
    return;
  }

  if( !m_geometry->isValid() )
  {
    m_estimate->setText( "" );
    return;
  }

  try
  {
    const ceelo::GeometryDescriptor gd = m_geometry->toDescriptor();
    ceelo::GenerationOptions opts;
    opts.node_fep_precision = selectedPrecision();
    opts.precision_profile = selectedRelaxMild()
                             ? ceelo::PrecisionProfile::RelaxMild
                             : ceelo::PrecisionProfile::Uniform;
    switch( m_profile->currentIndex() )
    {
      case 1: opts.profile = ceelo::ResponseProfile::FarField; break;
      case 2: opts.profile = ceelo::ResponseProfile::Contact; break;
      default: opts.profile = ceelo::ResponseProfile::General; break;
    }

    if( method == Method::QuickMc )
    {
      opts.transfer_mode = true;
      opts.n_anchor_angles = (m_anchorAngles->currentIndex() == 0) ? 1 : 3;
    }

    const int nodes = ceelo::ResponseGenerator::estimated_node_count( gd, opts );

    // Rough per-node wall time scales with 1/precision^2, capped at the
    //  per-node ceiling; the constants are ballpark (M1-class laptop). The
    //  relax_mild map cuts the expensive high-E node share ~2x (D0 memo) -
    //  but it only relaxes the angular/near scan stages, and transfer mode is
    //  (almost) all stage-1 backbone, so no relax credit is taken there.
    const double prec = opts.node_fep_precision;
    const double relax = (selectedRelaxMild() && (method == Method::FullMc)) ? 0.5 : 1.0;
    const double per_node_s = relax * std::min( 8.0, 0.05 * std::pow( 0.003/prec, 2.0 ) + 0.05 );
    const int total_min = std::max( 1, static_cast<int>( std::ceil( nodes * per_node_s / 60.0 ) ) );

    if( method == Method::QuickMc )
    {
      // Show the saving against a full General-profile characterization.
      ceelo::GenerationOptions full_opts = opts;
      full_opts.transfer_mode = false;
      full_opts.n_anchor_angles = 1;
      full_opts.profile = ceelo::ResponseProfile::General;
      const int full_nodes = ceelo::ResponseGenerator::estimated_node_count( gd, full_opts );
      m_estimate->setText( WString::tr("mmr-estimate-transfer")
                            .arg( nodes ).arg( total_min ).arg( full_nodes ) );
    }else
    {
      m_estimate->setText( WString::tr("mmr-estimate")
                            .arg( nodes ).arg( total_min ) );
    }
  }catch( std::exception & )
  {
    m_estimate->setText( "" );
  }
}//updateEstimate()


void MakeMcResponseForDrf::startGeneration()
{
  ceelo::GeometryDescriptor gd;
  try
  {
    gd = m_geometry->toDescriptor();
  }catch( std::exception &e )
  {
    m_status->setText( WString::fromUTF8( e.what() ) );
    return;
  }

  ++m_generationId;
  const int generation_id = m_generationId;

  m_result.reset();
  m_validationChanged.emit( false );

  const Method method = selectedMethod();

  if( method == Method::CurveTransfer )
  {
    // Instant, deterministic, no MC.  The anchor curves read the DRF, so they
    //  are built HERE on the session thread; the worker gets only value
    //  copies and pure-CeeLo inputs.
    CeeLoUtils::TransferAnchor anchor;
    try
    {
      anchor = CeeLoUtils::transferAnchorForDrf( m_seedDrf, gd, refDistanceOverrideCm() );
    }catch( std::exception &e )
    {
      m_status->setText( WString::tr("mmr-anchor-unusable")
                          .arg( WString::fromUTF8(e.what()) ) );
      return;
    }

    const ceelo::AnchorCurve tot_curve
                       = CeeLoUtils::totalTransferAnchorForDrf( m_seedDrf, anchor );
    const string det_name = m_seedDrf ? m_seedDrf->name() : string("user geometry");

    const string sessionId = wApp->sessionId();
    const string widgetId = id();

    auto worker = [gd,anchor,tot_curve,det_name,sessionId,widgetId,generation_id](){
      shared_ptr<ceelo::DetectorResponse> response;
      string errmsg;
      try
      {
        response = CeeLoUtils::makeTransferResponse( gd, anchor, tot_curve, det_name );
      }catch( std::exception &e )
      {
        errmsg = e.what();
      }

      WServer::instance()->post( sessionId, [widgetId,response,errmsg,generation_id](){
        auto *tool = dynamic_cast<MakeMcResponseForDrf *>( wApp->domRoot()->findById(widgetId) );
        if( tool )
          tool->handleGenerationFinished( response, errmsg, generation_id );
        wApp->triggerUpdate();
      } );
    };//worker

    m_status->setText( WString::tr("mmr-status-transfer-building") );

    wApp->enableUpdates( true );
    WServer::instance()->ioService().boost::asio::io_service::post( worker );
    return;
  }//if( method == Method::CurveTransfer )

  ceelo::GenerationOptions opts;
  opts.node_fep_precision = selectedPrecision();
  opts.precision_profile = selectedRelaxMild()
                           ? ceelo::PrecisionProfile::RelaxMild
                           : ceelo::PrecisionProfile::Uniform;
  switch( m_profile->currentIndex() )
  {
    case 1: opts.profile = ceelo::ResponseProfile::FarField; break;
    case 2: opts.profile = ceelo::ResponseProfile::Contact; break;
    default: opts.profile = ceelo::ResponseProfile::General; break;
  }

  if( method == Method::QuickMc )
  {
    // EFFTRAN-style transfer: MC only the on-axis energy backbone (plus a few
    //  forced cos-theta anchors when selected); the ray-traced kernel carries
    //  the distance/angle transfer.  generate() forces the FarField profile.
    opts.transfer_mode = true;
    opts.n_anchor_angles = (m_anchorAngles->currentIndex() == 0) ? 1 : 3;
  }

  opts.detector_name = m_seedDrf ? m_seedDrf->name() : string("user geometry");
  opts.base_seed = 1;  //deterministic; re-running the same setup reproduces

  m_cancelFlag = make_shared<std::atomic<bool>>( false );
  opts.cancel = m_cancelFlag;

  // Grounding anchors are captured NOW (value copies) - nothing from the
  //  widget tree crosses into the worker thread.
  bool curve_derived = false;
  const vector<ceelo::GroundingPoint> ground_pts
                      = groundingPointsForDrf( m_seedDrf, gd, curve_derived );

  const string sessionId = wApp->sessionId();
  const string widgetId = id();

  // Progress: throttled to whole-percent changes; findById(...) on the
  //  session thread is the only way widget access happens (never capture
  //  `this` or observing_ptr in the worker).
  auto last_pct = make_shared<std::atomic<int>>( -1 );
  opts.progress = [sessionId,widgetId,generation_id,last_pct]( double frac, const string &stage ){
    const int pct = static_cast<int>( 100.0 * frac );
    if( last_pct->exchange(pct) == pct )
      return;
    WServer::instance()->post( sessionId, [widgetId,frac,stage,generation_id](){
      auto *tool = dynamic_cast<MakeMcResponseForDrf *>( wApp->domRoot()->findById(widgetId) );
      if( tool )
        tool->updateProgress( frac, stage, generation_id );
      wApp->triggerUpdate();
    } );
  };

  auto worker = [gd,opts,ground_pts,curve_derived,sessionId,widgetId,generation_id](){
    shared_ptr<ceelo::DetectorResponse> response;
    string errmsg;
    try
    {
      response = ceelo::ResponseGenerator::generate( gd, opts );
      if( response && !ground_pts.empty() )
        ceelo::ResponseGenerator::ground_to_points( *response, ground_pts, curve_derived );
    }catch( ceelo::GenerationCancelled & )
    {
      errmsg = "cancelled";
    }catch( std::exception &e )
    {
      errmsg = e.what();
    }

    WServer::instance()->post( sessionId, [widgetId,response,errmsg,generation_id](){
      auto *tool = dynamic_cast<MakeMcResponseForDrf *>( wApp->domRoot()->findById(widgetId) );
      if( tool )
        tool->handleGenerationFinished( response, errmsg, generation_id );
      else
        cerr << "MakeMcResponseForDrf deleted while MC generation ran" << endl;
      wApp->triggerUpdate();
    } );
  };//worker

  m_generate->hide();
  m_cancelBtn->show();
  m_progress->setValue( 0.0 );
  m_progress->show();
  m_status->setText( WString::tr("mmr-status-running") );

  wApp->enableUpdates( true );

  WServer::instance()->ioService().boost::asio::io_service::post( worker );
}//startGeneration()


void MakeMcResponseForDrf::cancelGeneration()
{
  if( m_cancelFlag )
    m_cancelFlag->store( true );
  m_status->setText( WString::tr("mmr-status-cancelling") );
}//cancelGeneration()


void MakeMcResponseForDrf::updateProgress( const double frac, const std::string &stage,
                                           const int generation_id )
{
  if( generation_id != m_generationId )
    return;  //stale run

  m_progress->setValue( frac );
  m_status->setText( WString::fromUTF8(stage) );
}//updateProgress(...)


void MakeMcResponseForDrf::handleGenerationFinished(
                              std::shared_ptr<ceelo::DetectorResponse> result,
                              const std::string &errmsg,
                              const int generation_id )
{
  // Balance the enableUpdates(true) from startGeneration FIRST - every started
  //  generation posts exactly one finish, including runs made stale by a
  //  restart or method switch (Wt refcounts the update lock, so a stale
  //  finish must still decrement it or server push stays on forever).
  wApp->enableUpdates( false );

  if( generation_id != m_generationId )
    return;  //stale run - a newer run/state owns the UI

  m_generate->setHidden( selectedMethod() == Method::CurveTransfer );
  m_generate->setEnabled( m_geometry->isValid() );
  m_cancelBtn->hide();
  m_progress->hide();

  if( !result || !errmsg.empty() )
  {
    m_result.reset();
    m_validationChanged.emit( false );
    updateResponseChart();
    if( errmsg == "cancelled" )
      m_status->setText( WString::tr("mmr-status-cancelled") );
    else
      m_status->setText( WString::tr("mmr-status-error").arg( WString::fromUTF8(errmsg) ) );
    return;
  }//if( failed )

  m_result = result;
  m_validationChanged.emit( true );
  updateResponseChart();

  const bool grounded = !result->grounding.empty();
  if( result->model_transfer.has_value() )
  {
    // A transfer (quick-MC or measured-curve) response: state the validity
    //  floor - off-axis/near queries carry an honest, inflated uncertainty.
    const string min_dist = PhysicalUnits::printToBestLengthUnits(
                    result->provenance.min_distance_cm * PhysicalUnits::cm );
    m_status->setText( WString::tr("mmr-status-transfer-done")
                        .arg( min_dist ) );
  }else
  {
    m_status->setText( WString::tr( grounded ? "mmr-status-done-grounded"
                                             : "mmr-status-done" ) );
  }

  m_responseGenerated.emit( result );
}//handleGenerationFinished(...)


void MakeMcResponseForDrf::acceptResponse()
{
  if( !m_result )
    return;

  shared_ptr<SpecMeas> foreground = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  shared_ptr<DetectorPeakResponse> prev_det = foreground ? foreground->detector() : nullptr;

  const shared_ptr<const DetectorPeakResponse> base = m_seedDrf ? m_seedDrf : prev_det;

  shared_ptr<DetectorPeakResponse> new_det;
  if( base )
  {
    new_det = make_shared<DetectorPeakResponse>( *base );
    new_det->setParentHashValue( base->hashValue() );
  }else
  {
    // No existing DRF: make a minimal legacy shell around the MC response.
    const double a_cm = m_result->transverse_half_extent();
    new_det = make_shared<DetectorPeakResponse>( "MC characterized detector",
                                          "Monte-Carlo parameterized response" );
    new_det->setIntrinsicEfficiencyFormula( "1.0", 2.0*a_cm*PhysicalUnits::cm,
                    PhysicalUnits::keV, 0.0f, 0.0f,
                    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  }//if( base ) / else

  new_det->setCeeloResponse( m_result );

  m_interspec->detectorChanged().emit( new_det );

  UndoRedoManager * const undoManager = m_interspec->undoRedoManager();
  if( undoManager && undoManager->canAddUndoRedoNow() )
  {
    auto undo = [prev_det](){
      InterSpec::instance()->detectorChanged().emit( prev_det );
    };
    auto redo = [new_det](){
      InterSpec::instance()->detectorChanged().emit( new_det );
    };
    undoManager->addUndoRedoStep( undo, redo, "Attach MC response to DRF" );
  }//if( undoManager )

  m_updatedDrf.emit( new_det );
}//acceptResponse()


MakeMcResponseForDrf *MakeMcResponseForDrfWindow::tool()
{
  return m_tool;
}


MakeMcResponseForDrfWindow::MakeMcResponseForDrfWindow(
                          std::shared_ptr<const DetectorPeakResponse> seed_drf )
 : AuxWindow( WString::tr("window-title-mc-response"),
             (AuxWindowProperties::TabletNotFullScreen
              | AuxWindowProperties::SetCloseable
              | AuxWindowProperties::DisableCollapse
              | AuxWindowProperties::EnableResize
              | AuxWindowProperties::IsModal) ),
  m_tool( nullptr )
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( !viewer )
    return;

  viewer->useMessageResourceBundle( "MakeMcResponseForDrf" );

  const int ww = viewer->renderedWidth();
  const int wh = viewer->renderedHeight();
  if( ww > 100 && wh > 100 )
  {
    const int width = std::min( (8*ww)/9, 640 );
    const int height = std::min( 700, ((wh < 420) ? wh : (19*wh)/20 ) );
    resizeWindow( width, height );
    setMinimumSize( std::min(width,480), std::min(height,400) );
  }//if( ww > 100 && wh > 100 )

  if( !seed_drf )
  {
    shared_ptr<SpecMeas> foreground = viewer->measurment( SpecUtils::SpectrumType::Foreground );
    seed_drf = foreground ? foreground->detector() : nullptr;
  }

  {
    auto toolOwner = std::make_unique<MakeMcResponseForDrf>( viewer, seed_drf );
    m_tool = toolOwner.get();
    stretcher()->addWidget( std::move(toolOwner), 0, 0 );
  }
  stretcher()->setContentsMargins( 0, 0, 0, 0 );

  AuxWindow::addHelpInFooter( footer(), "make-mc-response" );

  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close") );
  closeButton->clicked().connect( this, &AuxWindow::hide );

  WPushButton *useBtn = footer()->addNew<WPushButton>( WString::tr("mmr-use-response-btn") );
  useBtn->clicked().connect( m_tool, &MakeMcResponseForDrf::acceptResponse );
  m_tool->validationChanged().connect( useBtn, [useBtn]( bool valid ){
    useBtn->setEnabled( valid );
  } );
  useBtn->disable();

  show();
  centerWindowHeavyHanded();
}//MakeMcResponseForDrfWindow constructor
