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
#include <atomic>
#include <chrono>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <functional>

#include <boost/asio/io_service.hpp>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WTable.h>
#include <Wt/WTimer.h>
#include <Wt/WServer.h>
#include <Wt/WCheckBox.h>
#include <Wt/WTableRow.h>
#include <Wt/WIOService.h>
#include <Wt/WTableCell.h>
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
#include "efficiency/EfficiencyCalculator.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfChart.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/MakeMcResponseForDrf.h"
#include "InterSpec/DetectorGeometryInput.h"

using namespace Wt;
using namespace std;


/** Measured Monte-Carlo throughput for one geometry, from which every planned node's cost is
 extrapolated: the events a node needs for a fractional precision p are `rel_var_per_event / p^2`
 (the probe's own (sigma/eps)^2 x N, which folds in CeeLo's variance-reduction biasing - the analog
 (1 - eps)/eps rule would not), and they run at `events_per_cpu_s` across `parallelism` threads.
 */
struct McTimeCalibration
{
  std::string geometry_key;        //ceelo::GeometryDescriptor::to_xml_string() it was measured for
  double rel_var_per_event = 0.0;
  double events_per_cpu_s = 0.0;
  double parallelism = 1.0;        //cpu_s / wall_s
  bool from_full_run = false;      //derived from a real generation's per-node costs, not the probe

  bool valid() const { return (rel_var_per_event > 0.0) && (events_per_cpu_s > 0.0); }
};//struct McTimeCalibration


/** Shared with the generation worker: atomics only, so no lock and nothing widget-side crosses. */
struct McProgressSnapshot
{
  std::atomic<int> nodes_done{ 0 };
  std::atomic<int> nodes_total{ 0 };
  std::atomic<int> stage{ 0 };   //ceelo::NodeProgress::stage of the last finished node
};//struct McProgressSnapshot


namespace
{
  /** The per-node target precision a run would use (an explicit map, the graded RelaxMild map, or
   the flat scalar) - built once, since the graded map is a std::function. */
  std::function<double(uint32_t,double)> node_precision_map( const ceelo::GenerationOptions &opts )
  {
    if( opts.node_precision )
      return opts.node_precision;
    if( opts.precision_profile == ceelo::PrecisionProfile::RelaxMild )
      return ceelo::relax_mild_precision_map( opts.node_fep_precision );
    const double base = opts.node_fep_precision;
    return [base]( uint32_t, double ) -> double { return base; };
  }//node_precision_map(...)


  /** Predicted wall-seconds for one node at fractional precision `prec`: from the measured
   throughput when there is one, else the historical ballpark (an M1-class laptop, capped like the
   per-node budget); either way the per-node event/CPU caps apply as in CeeLo's apply_node_budget.
   */
  double node_cost_s( const double prec, const ceelo::GenerationOptions &opts,
                      const McTimeCalibration *calib )
  {
    if( !calib || !calib->valid() )
    {
      // Measured, not guessed: a full 814-node run of a 3x3 NaI + 1 mm Al can on a 10-thread M1
      //  took 351 s, i.e. 0.431 s/node at the default 0.3% precision.  The 1/p^2 shape is right
      //  (node cost is dominated by the events needed to reach `prec`), the old 0.05 coefficient
      //  was not: it predicted 0.100 s/node, so the estimate read 1.4 min for a 5.9 min run and
      //  6.8 min for a 43 min one.  Better to be honest before the probe refines it.
      return std::min( 8.0, 0.22 * std::pow( 0.003/prec, 2.0 ) + 0.21 );
    }

    const double needed = calib->rel_var_per_event / (prec * prec);
    const double events = std::min( std::max( needed, double(opts.min_events_per_node) ),
                                    double(opts.max_events_per_node) );
    double cpu = events / calib->events_per_cpu_s;
    if( opts.max_cpu_seconds_per_node > 0.0 )
      cpu = std::min( cpu, opts.max_cpu_seconds_per_node );

    return cpu / std::max( 1.0, calib->parallelism ) + 0.03;  //+ per-node setup (calculator, quadrature)
  }//node_cost_s(...)


  /** Predicted cumulative wall-seconds after each planned node, in generation order (backbone
   energies ascending, then the angular nodes of each shape energy, then the near-field nodes of
   each); size total + 1, [0] = 0. */
  std::vector<double> prior_cumulative_cost( const ceelo::ResponseGenerator::NodePlan &plan,
                                             const ceelo::GenerationOptions &opts,
                                             const McTimeCalibration *calib )
  {
    const std::function<double(uint32_t,double)> prec = node_precision_map( opts );

    std::vector<double> cum;
    cum.reserve( static_cast<size_t>( std::max( 0, plan.total() ) ) + 1 );
    cum.push_back( 0.0 );

    auto push = [&]( const uint32_t stage, const double energy ){
      cum.push_back( cum.back() + node_cost_s( prec(stage, energy), opts, calib ) );
    };

    for( const double energy : plan.backbone_energies_keV )
      push( 1, energy );
    for( const double energy : plan.shape_energies_keV )
      for( int i = 0; i < plan.n_cos_theta * plan.n_phi; ++i )
        push( 2, energy );
    for( const double energy : plan.shape_energies_keV )
      for( int i = 0; i < plan.n_near_positions; ++i )
        push( 3, energy );

    return cum;
  }//prior_cumulative_cost(...)


  /** "0.4", "5.2", "12" minutes - pre-formatted, because WString::arg(double) is locale-formatted
   with no say over the decimals. */
  std::string minutes_str( const double seconds )
  {
    const double minutes = seconds / 60.0;
    char buf[32];
    snprintf( buf, sizeof(buf), (minutes < 10.0) ? "%.1f" : "%.0f", minutes );
    return buf;
  }//minutes_str(...)


  /** "under a minute" / "N minute(s)" / "H.h hour(s)" for a predicted duration. */
  WString duration_phrase( const double seconds )
  {
    if( seconds < 60.0 )
      return WString::tr("mmr-est-under-minute");

    const double minutes = seconds / 60.0;
    if( minutes >= 90.0 )
    {
      char buf[32];
      snprintf( buf, sizeof(buf), "%.1f", minutes / 60.0 );
      return WString::tr("mmr-est-hours").arg( string(buf) );
    }

    return WString::tr("mmr-est-minutes").arg( static_cast<int>( std::ceil(minutes) ) );
  }//duration_phrase(...)
}//namespace


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
    //  reference distance.  (farFieldIntrinsicEfficiency() already backs out any air
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

    // One covariance query: for a DRF carrying a response this traces a full aperture fan.
    const vector<double> cov = drf->efficiencyFracCovariance( energies );
    vector<double> sigmas;
    if( cov.size() == (energies.size()*energies.size()) )
    {
      for( size_t i = 0; i < energies.size(); ++i )
        sigmas.push_back( std::sqrt( std::max(0.0, cov[i*energies.size()+i]) ) );
    }

    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double intrinsic = drf->farFieldIntrinsicEfficiency( static_cast<float>(energies[i]) );
      if( intrinsic <= 0.0 )
        continue;

      ceelo::GroundingPoint gp;
      gp.energy_keV = energies[i];
      gp.measured_eff = intrinsic * frac_solid_angle;
      gp.frac_stat_sigma = sigmas.empty() ? CeeLoUtils::sm_default_anchor_frac_sigma
                                          : std::max( CeeLoUtils::sm_min_anchor_frac_sigma, sigmas[i] );
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
    m_seedProvider(),
    m_geometry( nullptr ),
    m_method( nullptr ),
    m_profile( nullptr ),
    m_precision( nullptr ),
    m_customPrecision( nullptr ),
    m_anchorAngles( nullptr ),
    m_groundCb( nullptr ),
    m_groundInfo( nullptr ),
    m_groundRow( nullptr ),
    m_groundInfoRow( nullptr ),
    m_anchorInfo( nullptr ),
    m_refDistance( nullptr ),
    m_estimate( nullptr ),
    m_profileRow( nullptr ),
    m_precRow( nullptr ),
    m_anchorAnglesRow( nullptr ),
    m_anchorInfoRow( nullptr ),
    m_anchorRow( nullptr ),
    m_generate( nullptr ),
    m_hideGenerateButton( false ),
    m_cancelBtn( nullptr ),
    m_progress( nullptr ),
    m_status( nullptr ),
    m_chartBox( nullptr ),
    m_chart( nullptr ),
    m_chartMode( nullptr ),
    m_chartDistance( nullptr ),
    m_chartRangeSet( false ),
    m_restoringState( false ),
    m_seedingFromDrf( false ),
    m_generationId( 0 ),
    m_cancelFlag( nullptr ),
    m_result( nullptr ),
    m_validationChanged(),
    m_geometryChanged(),
    m_hideChart( false ),
    m_responseGenerated(),
    m_updatedDrf(),
    m_calibration( nullptr ),
    m_calibrationId( 0 ),
    m_calibrating( false ),
    m_shown( false ),
    m_calibTimer( nullptr ),
    m_generating( false ),
    m_progressSnapshot( nullptr ),
    m_progressTimer( nullptr ),
    m_generationStart(),
    m_runGeometryKey(),
    m_nodesTotal( 0 ),
    m_priorCumulative(),
    m_transferAnchors( false )
{
  assert( m_interspec );
  wApp->useStyleSheet( "InterSpec_resources/MakeMcResponseForDrf.css" );
  m_interspec->useMessageResourceBundle( "MakeMcResponseForDrf" );

  addStyleClass( "MakeMcResponseForDrf" );

  //Geometry
  WGroupBox *geomBox = addNew<WGroupBox>( WString::tr("mmr-geometry-title") );
  geomBox->addStyleClass( "McGeomBox" );
  m_geometry = geomBox->addNew<DetectorGeometryInput>( m_interspec );
  m_geometry->seedFromDrf( m_seedDrf );  //uses the DRFs own geometry when it has one
  m_geometry->changed().connect( this, &MakeMcResponseForDrf::handleGeometryChanged );

  //Characterization options - a label/input table, so the inputs line up with each other and with
  //  the geometry form above (see `DgiTable` in DetectorGeometryInput).
  WGroupBox *optsBox = addNew<WGroupBox>( WString::tr("mmr-options-title") );
  optsBox->addStyleClass( "McOptsBox" );

  WTable *optsTable = optsBox->addNew<WTable>();
  optsTable->addStyleClass( "McOptsTable" );

  int opt_row = 0;
  auto add_row = [optsTable,&opt_row]( const WString &label ) -> WTableCell * {
    optsTable->elementAt( opt_row, 0 )->addNew<WLabel>( label );
    WTableCell *cell = optsTable->elementAt( opt_row, 1 );
    ++opt_row;
    return cell;
  };

  //A full-width row, for the notes/estimates that read as prose rather than a labeled field.
  auto add_wide_row = [optsTable,&opt_row]() -> WTableCell * {
    WTableCell *cell = optsTable->elementAt( opt_row, 0 );
    cell->setColumnSpan( 2 );
    ++opt_row;
    return cell;
  };

  {
    WTableCell *cell = add_row( WString::tr("mmr-method") );
    m_method = cell->addNew<WComboBox>();
  }
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
    WText *upgradeNote = add_wide_row()->addNew<WText>( WString::tr("mmr-upgrade-note") );
    upgradeNote->addStyleClass( "McUpgradeNote" );
    upgradeNote->setInline( false );
  }

  m_profile = add_row( WString::tr("mmr-profile") )->addNew<WComboBox>();
  m_profileRow = optsTable->rowAt( opt_row - 1 );
  m_profile->addItem( WString::tr("mmr-profile-general") );
  m_profile->addItem( WString::tr("mmr-profile-farfield") );
  m_profile->addItem( WString::tr("mmr-profile-contact") );
  m_profile->setCurrentIndex( 0 );
  m_profile->activated().connect( this, &MakeMcResponseForDrf::handleOptionChanged );
  HelpSystem::attachToolTipOn( m_profile, WString::tr("mmr-tt-profile"), true );

  WTableCell *precCell = add_row( WString::tr("mmr-precision") );
  m_precRow = optsTable->rowAt( opt_row - 1 );
  m_precision = precCell->addNew<WComboBox>();
  m_precision->addItem( WString::tr("mmr-precision-fast") );      //idx0: 1%
  m_precision->addItem( WString::tr("mmr-precision-normal") );    //idx1: 0.3%
  m_precision->addItem( WString::tr("mmr-precision-balanced") );  //idx2: relax_mild (0.3% base, relaxed high-E)
  m_precision->addItem( WString::tr("mmr-precision-thorough") );  //idx3: 0.1%
  m_precision->addItem( WString::tr("mmr-precision-custom") );    //idx4
  m_precision->setCurrentIndex( 1 );
  m_precision->activated().connect( this, &MakeMcResponseForDrf::handlePrecisionChanged );
  HelpSystem::attachToolTipOn( m_precision, WString::tr("mmr-tt-precision"), true );

  m_customPrecision = precCell->addNew<WLineEdit>();
  m_customPrecision->setTextSize( 5 );
  m_customPrecision->setPlaceholderText( "0.3%" );
  m_customPrecision->changed().connect( this, &MakeMcResponseForDrf::handleOptionChanged );
  m_customPrecision->hide();

  m_anchorAngles = add_row( WString::tr("mmr-anchor-angles") )->addNew<WComboBox>();
  m_anchorAnglesRow = optsTable->rowAt( opt_row - 1 );
  m_anchorAngles->addItem( WString::tr("mmr-anchor-angles-1") );  //on-axis only
  m_anchorAngles->addItem( WString::tr("mmr-anchor-angles-3") );  //3 cos-theta anchors
  m_anchorAngles->setCurrentIndex( 1 );  //3 anchors: off-axis FEP + real HPGe total
  m_anchorAngles->activated().connect( this, &MakeMcResponseForDrf::handleOptionChanged );
  HelpSystem::attachToolTipOn( m_anchorAngles, WString::tr("mmr-tt-anchor-angles"), true );
  m_anchorAnglesRow->hide();

  // Whether the MC response gets corrected to the detectors own measured efficiency, and what it
  //  would be corrected against - the difference between "this detector, measured" and "a detector
  //  of this shape, computed".
  {
    WTableCell *cell = add_row( WString::tr("mmr-grounding") );
    m_groundRow = optsTable->rowAt( opt_row - 1 );
    m_groundCb = cell->addNew<WCheckBox>( WString::tr("mmr-grounding-cb") );
    m_groundCb->setChecked( true );
    m_groundCb->changed().connect( this, &MakeMcResponseForDrf::handleOptionChanged );
    HelpSystem::attachToolTipOn( m_groundCb, WString::tr("mmr-tt-grounding"), true );
  }

  m_groundInfo = add_wide_row()->addNew<WText>( "" );
  m_groundInfoRow = optsTable->rowAt( opt_row - 1 );
  m_groundInfo->addStyleClass( "McAnchorInfo" );
  m_groundInfo->setInline( false );

  m_anchorInfo = add_wide_row()->addNew<WText>( "" );
  m_anchorInfoRow = optsTable->rowAt( opt_row - 1 );
  m_anchorInfo->addStyleClass( "McAnchorInfo" );
  m_anchorInfo->setInline( false );
  m_anchorInfoRow->hide();

  m_refDistance = add_row( WString::tr("mmr-ref-distance") )->addNew<WLineEdit>();
  m_anchorRow = optsTable->rowAt( opt_row - 1 );
  m_refDistance->setTextSize( 10 );
  m_refDistance->changed().connect( this, &MakeMcResponseForDrf::handleGeometryChanged );
  HelpSystem::attachToolTipOn( m_refDistance, WString::tr("mmr-tt-ref-distance"), true );
  m_anchorRow->hide();

  m_estimate = add_wide_row()->addNew<WText>( "" );
  m_estimate->addStyleClass( "McEstimate" );
  m_estimate->setInline( false );

  //Run row - inside the group box, since running is what these options configure.
  WContainerWidget *runRow = optsBox->addNew<WContainerWidget>();
  runRow->addStyleClass( "McRunRow" );
  m_status = runRow->addNew<WText>( "" );
  m_status->addStyleClass( "McStatus" );
  m_generate = runRow->addNew<WPushButton>( WString::tr("mmr-generate-btn") );
  //A lambda, not the slot directly: startGeneration() returns whether it actually started.
  m_generate->clicked().connect( this, [this](){ startGeneration(); } );
  m_cancelBtn = runRow->addNew<WPushButton>( WString::tr("Cancel") );
  m_cancelBtn->clicked().connect( this, &MakeMcResponseForDrf::cancelGeneration );
  m_cancelBtn->hide();
  m_progress = optsBox->addNew<WProgressBar>();
  m_progress->setRange( 0.0, 1.0 );
  m_progress->hide();

  // The calibration probe waits out a burst of edits; the progress timer only runs during a
  //  generation.  Plain WObjects owned here (no widget parent).
  m_calibTimer = make_unique<WTimer>();
  m_calibTimer->setSingleShot( true );
  m_calibTimer->setInterval( std::chrono::milliseconds(1500) );
  m_calibTimer->timeout().connect( this, &MakeMcResponseForDrf::startTimeCalibration );

  m_progressTimer = make_unique<WTimer>();
  m_progressTimer->setInterval( std::chrono::milliseconds(2000) );
  m_progressTimer->timeout().connect( this, &MakeMcResponseForDrf::refreshProgressText );

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
  m_chartMode->activated().connect( this, &MakeMcResponseForDrf::handleChartOptionChanged );

  WLabel *distLabel = chartCtrlRow->addNew<WLabel>( WString::tr("mmr-chart-distance") );
  m_chartDistance = chartCtrlRow->addNew<WLineEdit>( "25 cm" );
  m_chartDistance->setTextSize( 8 );
  distLabel->setBuddy( m_chartDistance );
  m_chartDistance->changed().connect( this, &MakeMcResponseForDrf::handleChartOptionChanged );
  m_chartDistance->enterPressed().connect( this, &MakeMcResponseForDrf::handleChartOptionChanged );
  HelpSystem::attachToolTipOn( m_chartDistance, WString::tr("mmr-tt-chart-distance"), true );

  m_chart = m_chartBox->addNew<DrfChart>();
  m_chart->addStyleClass( "McChart" );
  m_chart->setShowFwhm( false );  //this chart previews the geometrys efficiency, not the FWHM

  m_chartBox->hide();

  // Open on the support the detector actually has, rather than on a default the user then has to
  //  correct: the method that produced its response, or - for a detector that states its shape and
  //  has an efficiency but no response - the no-MC transfer, which is the only support that can be
  //  had without a run, and which builds itself the moment the geometry is usable.
  //
  // Deliberately NOT for a DRF whose geometry is only the cylinder guessed from its diameter: a
  //  transfer built on a guessed shape would be attached, silently, by the next "Use".
  if( seed_resp )
    m_method->setCurrentIndex( static_cast<int>(seed_resp->provenance.method) );
  else if( m_seedDrf && m_seedDrf->isValid() && m_seedDrf->geometry() )
    m_method->setCurrentIndex( static_cast<int>(Method::CurveTransfer) );

  m_seedingFromDrf = true;
  handleMethodChanged();   //row visibility for the selected method; clears m_result
  m_seedingFromDrf = false;

  if( seed_resp )
  {
    // The response the DRF already carries; showing it is the point of opening on its method.
    m_result = seed_resp;
    m_status->setText( WString::tr("mmr-status-existing") );
    m_validationChanged.emit( true );
    updateResponseChart();
  }else
  {
    handleGeometryChanged();  //lets the no-MC transfer build itself, now that seeding is done
  }
}//MakeMcResponseForDrf constructor


MakeMcResponseForDrf::~MakeMcResponseForDrf()
{
  // The worker only holds the cancel flag and plain strings - never `this` -
  //  so all we need to do is ask any in-flight generation to stop.
  if( m_cancelFlag )
    m_cancelFlag->store( true );
  if( m_calibTimer )
    m_calibTimer->stop();
  if( m_progressTimer )
    m_progressTimer->stop();
}//~MakeMcResponseForDrf()


void MakeMcResponseForDrf::setGeometryFromDescriptor( const ceelo::GeometryDescriptor &geometry,
                                                      const std::vector<std::string> &notes )
{
  m_geometry->setFromDescriptor( geometry, notes );
  handleGeometryChanged();   //refresh estimate / anchor / any generated result
}//setGeometryFromDescriptor(...)


void MakeMcResponseForDrf::setSeedDrf( std::shared_ptr<const DetectorPeakResponse> seed_drf )
{
  m_seedDrf = seed_drf;

  // Only refresh the rows that read the seed; deliberately NOT seedFromDrf()/method resets/
  //  auto-generate - the owner has already applied its edits and drives generation itself.
  updateAnchorInfo();
  updateGroundingInfo();
  updateResponseChart();
  updateEstimate();
}//setSeedDrf(...)


void MakeMcResponseForDrf::setSeedProvider( std::function<std::shared_ptr<const DetectorPeakResponse>()> provider )
{
  m_seedProvider = std::move( provider );
  refreshSeedFromProvider();
}//setSeedProvider(...)


void MakeMcResponseForDrf::refreshSeedFromProvider()
{
  if( !m_seedProvider )
    return;

  // Pulled before every generation, so an automatic rebuild cannot anchor on a seed that predates
  //  the owner's edits (which would produce a response that ignores them while looking current).
  const std::shared_ptr<const DetectorPeakResponse> seed = m_seedProvider();
  if( seed )
    setSeedDrf( seed );
}//refreshSeedFromProvider()


void MakeMcResponseForDrf::setGenerateButtonHidden( bool hidden )
{
  m_hideGenerateButton = hidden;
  if( m_generate && hidden )
    m_generate->hide();
}//setGenerateButtonHidden(...)


bool MakeMcResponseForDrf::generationReady() const
{
  return m_geometry->generationReady();
}//generationReady()


std::string MakeMcResponseForDrf::geometryProblem() const
{
  return m_geometry->problemDescription();
}//geometryProblem()


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


int MakeMcResponseForDrf::generationId() const
{
  return m_generationId;
}


bool MakeMcResponseForDrf::isGenerating() const
{
  return m_generating;
}


Wt::Signal<> &MakeMcResponseForDrf::userChangedNoRegen()
{
  return m_userChangedNoRegen;
}


Wt::Signal<std::shared_ptr<ceelo::DetectorResponse>> &MakeMcResponseForDrf::responseGenerated()
{
  return m_responseGenerated;
}


Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &MakeMcResponseForDrf::updatedDrf()
{
  return m_updatedDrf;
}


bool MakeMcResponseForDrf::State::operator==( const State &rhs ) const
{
  return (method == rhs.method)
         && (profile == rhs.profile)
         && (precision == rhs.precision)
         && (anchorAngles == rhs.anchorAngles)
         && (chartMode == rhs.chartMode)
         && (groundToMeasured == rhs.groundToMeasured)
         && (customPrecision == rhs.customPrecision)
         && (refDistance == rhs.refDistance)
         && (chartDistance == rhs.chartDistance)
         && (result == rhs.result)
         && (status == rhs.status)
         && (geometry == rhs.geometry);
}//State::operator==


MakeMcResponseForDrf::State MakeMcResponseForDrf::currentState() const
{
  State state;

  state.method = m_method->currentIndex();
  state.profile = m_profile->currentIndex();
  state.precision = m_precision->currentIndex();
  state.anchorAngles = m_anchorAngles->currentIndex();
  state.chartMode = m_chartMode->currentIndex();
  state.groundToMeasured = m_groundCb->isChecked();
  state.customPrecision = m_customPrecision->text().toUTF8();
  state.refDistance = m_refDistance->text().toUTF8();
  state.chartDistance = m_chartDistance->text().toUTF8();
  state.result = m_result;
  state.status = m_status->text().toUTF8();
  state.geometry = m_geometry->currentState();

  return state;
}//MakeMcResponseForDrf::currentState()


void MakeMcResponseForDrf::setState( const State &state )
{
  m_restoringState = true;

  // A generation that is still running would land on top of the state being restored; its finish
  //  handler is stale-guarded by the generation id.
  ++m_generationId;
  if( m_cancelFlag )
    m_cancelFlag->store( true );
  m_progress->hide();
  m_cancelBtn->hide();
  m_generating = false;
  if( m_progressTimer )
    m_progressTimer->stop();
  ++m_calibrationId;   //a probe in flight belongs to the state being replaced
  m_calibrating = false;
  if( m_calibTimer )
    m_calibTimer->stop();

  m_method->setCurrentIndex( state.method );
  m_profile->setCurrentIndex( state.profile );
  m_precision->setCurrentIndex( state.precision );
  m_anchorAngles->setCurrentIndex( state.anchorAngles );
  m_chartMode->setCurrentIndex( state.chartMode );
  m_groundCb->setChecked( state.groundToMeasured );
  m_customPrecision->setText( WString::fromUTF8(state.customPrecision) );
  m_refDistance->setText( WString::fromUTF8(state.refDistance) );
  m_chartDistance->setText( WString::fromUTF8(state.chartDistance) );

  m_geometry->setState( state.geometry );

  const bool had_result = !!m_result;
  m_result = state.result;

  // Re-sync the per-method row visibility and the estimate/anchor text without going through
  //  `handleMethodChanged`, which exists to *invalidate* a result on a user edit.
  const Method method = selectedMethod();
  m_profileRow->setHidden( method != Method::FullMc );
  m_precRow->setHidden( method == Method::CurveTransfer );
  m_anchorAnglesRow->setHidden( method != Method::QuickMc );
  m_anchorInfoRow->setHidden( method != Method::CurveTransfer );
  m_anchorRow->setHidden( method != Method::CurveTransfer );
  m_generate->setHidden( m_hideGenerateButton || (method == Method::CurveTransfer) );
  m_generate->setEnabled( m_geometry->generationReady() );
  m_customPrecision->setHidden( m_precision->currentIndex() != 4 );

  m_status->setText( WString::fromUTF8(state.status) );
  updateAnchorInfo();
  updateGroundingInfo();
  updateEstimate();
  updateResponseChart();

  if( had_result != !!m_result )
    m_validationChanged.emit( !!m_result );

  m_restoringState = false;
}//MakeMcResponseForDrf::setState(...)


Wt::Signal<> &MakeMcResponseForDrf::userChanged()
{
  return m_userChanged;
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
  m_anchorInfoRow->setHidden( method != Method::CurveTransfer );
  m_anchorRow->setHidden( method != Method::CurveTransfer );

  // The measured-curve transfer is instant and rebuilds automatically - no
  //  explicit "Generate" step (this is what guarantees that entering geometry
  //  and accepting the dialog always yields an attached, distance-aware
  //  response).
  m_generate->setHidden( m_hideGenerateButton || (method == Method::CurveTransfer) );

  // A result from a different method is not what the user is configuring;
  //  abandon any in-flight generation too (its finish handler is stale-guarded
  //  and balances the update lock itself), and take back the run-row UI.
  ++m_generationId;
  if( m_cancelFlag )
    m_cancelFlag->store( true );
  m_progress->hide();
  m_cancelBtn->hide();
  m_generating = false;
  if( m_progressTimer )
    m_progressTimer->stop();
  if( m_result )
  {
    m_result.reset();
    m_validationChanged.emit( false );
  }
  m_status->setText( "" );

  handleGeometryChanged();

  if( !m_restoringState )
    m_userChanged.emit();
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

  m_generate->setEnabled( m_geometry->generationReady() );
  updateAnchorInfo();
  updateGroundingInfo();
  updateEstimate();
  scheduleTimeCalibration();  //a changed geometry has a different Monte-Carlo cost

  // Auto-build the instant transfer whenever the inputs are usable - except while the constructor
  //  is still seeding, where a response the DRF already carries is about to be installed.
  if( !m_seedingFromDrf && (selectedMethod() == Method::CurveTransfer) && m_geometry->generationReady() )
    startGeneration();

  if( !m_restoringState )
    m_userChanged.emit();

  m_geometryChanged.emit();
}//handleGeometryChanged()


CeeLoUtils::TransferAnchor MakeMcResponseForDrf::transferAnchor( const ceelo::GeometryDescriptor &gd ) const
{
  // A DRF whose curve carries a covariance (a fitted equation's coefficient covariance, or per-point
  //  node covariance) is anchored on that curve with the covariance, so correlations between
  //  energies survive into the response; otherwise the raw-point / sampled-curve anchor.
  const shared_ptr<const DetectorEfficiencyCurve> curve = m_seedDrf ? m_seedDrf->efficiencyCurve() : nullptr;
  const bool have_curve_cov = (curve && curve->uncertainty() && !curve->uncertainty()->isEmpty());
  if( have_curve_cov )
    return CeeLoUtils::curveAnchorWithCovarianceForDrf( m_seedDrf, gd, refDistanceOverrideCm() );
  return CeeLoUtils::transferAnchorForDrf( m_seedDrf, gd, refDistanceOverrideCm() );
}//transferAnchor(...)


void MakeMcResponseForDrf::setMethod( const Method method )
{
  m_method->setCurrentIndex( static_cast<int>(method) );
  handleMethodChanged();
}//setMethod(...)


void MakeMcResponseForDrf::setChartHidden( const bool hidden )
{
  m_hideChart = hidden;
  if( m_chartBox && hidden )
    m_chartBox->hide();
  else if( m_chartBox && m_result )
    updateResponseChart();
}//setChartHidden(...)


Wt::Signal<> &MakeMcResponseForDrf::geometryChanged()
{
  return m_geometryChanged;
}


bool MakeMcResponseForDrf::geometryValid() const
{
  return m_geometry->isValid();
}


ceelo::GeometryDescriptor MakeMcResponseForDrf::geometryDescriptor() const
{
  return m_geometry->toDescriptor();
}


DetectorGeometryInput *MakeMcResponseForDrf::geometryInput()
{
  return m_geometry;
}


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
    const CeeLoUtils::TransferAnchor anchor = transferAnchor( gd );

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


bool MakeMcResponseForDrf::groundToMeasured() const
{
  return (m_groundCb && m_groundCb->isEnabled() && m_groundCb->isChecked());
}//bool groundToMeasured() const


void MakeMcResponseForDrf::updateGroundingInfo()
{
  if( !m_groundCb || !m_groundInfo || !m_groundRow || !m_groundInfoRow )
    return;

  // The measured-curve transfer IS anchored on the measured efficiency by construction, so a
  //  separate grounding choice would be meaningless there - it has its own anchor row.
  const bool applies = (selectedMethod() != Method::CurveTransfer);
  m_groundRow->setHidden( !applies );
  m_groundInfoRow->setHidden( !applies );
  if( !applies )
  {
    m_groundInfo->setText( "" );
    return;
  }//if( !applies )

  size_t npoints = 0;
  bool curve_derived = false;
  if( m_geometry->isValid() )
  {
    try
    {
      const ceelo::GeometryDescriptor gd = m_geometry->toDescriptor();
      npoints = groundingPointsForDrf( m_seedDrf, gd, curve_derived ).size();
    }catch( std::exception & )
    {
    }
  }//if( m_geometry->isValid() )

  m_groundCb->setEnabled( npoints > 0 );

  if( npoints == 0 )
    m_groundInfo->setText( WString::tr("mmr-ground-none") );
  else if( !groundToMeasured() )
    m_groundInfo->setText( WString::tr("mmr-ground-off") );
  else
    m_groundInfo->setText( WString::tr( curve_derived ? "mmr-ground-curve"
                                                      : "mmr-ground-points" )
                            .arg( static_cast<int>(npoints) ) );
}//void updateGroundingInfo()


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

  if( !m_result || m_hideChart )
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
    preview = make_shared<DetectorPeakResponse>( *base );

  // A seed with no efficiency curve of its own (a geometry-only import) is not "valid", and the
  //  chart refuses an invalid DRF for both the curve and the per-angle series - so give the preview
  //  the backbone curve sampled from the response, exactly as accepting it would.  Fall back to a
  //  unit formula shell (which is valid) if that fails, or when there is no seed at all.
  if( preview && !preview->isValid() )
  {
    try
    {
      CeeLoUtils::setLegacyEfficiencyFromResponse( *preview, m_result );
    }catch( std::exception & )
    {
      preview.reset();
    }
  }//if( preview && !preview->isValid() )

  if( !preview )
  {
    const double a_cm = m_result->transverse_half_extent();
    preview = make_shared<DetectorPeakResponse>( "preview", "" );
    preview->setIntrinsicEfficiencyFormula( "1.0", 2.0*a_cm*PhysicalUnits::cm,
                    PhysicalUnits::keV, 0.0f, 0.0f,
                    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  }//if( !preview )

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

  // Start on the DRFs valid range, but no further than 3 MeV - above that the curves are mostly
  //  empty chart.  Only the first showing picks a range (`updateChart` resets it from the DRFs
  //  JSON extent), so re-draws for a mode/distance change dont undo the users zooming.
  if( !m_chartRangeSet )
  {
    const double lower = std::max( 0.0, preview->lowerEnergy() );
    double upper = preview->upperEnergy();
    if( upper <= (lower + 100.0) )   //the DRF didnt declare a usable range
      upper = 3000.0;
    m_chart->setXAxisRange( lower, std::min( upper, 3000.0 ) );
    m_chartRangeSet = true;
  }//if( !m_chartRangeSet )

  m_chartBox->show();
}//updateResponseChart()


void MakeMcResponseForDrf::handleOptionChanged()
{
  updateEstimate();

  if( !m_restoringState )
    m_userChanged.emit();
}//handleOptionChanged()


void MakeMcResponseForDrf::handleChartOptionChanged()
{
  updateResponseChart();

  // `userChangedNoRegen`, not `userChanged`: the chart distance and the absolute/intrinsic mode only
  //  change how the Response preview is DRAWN.  Reporting them as edits made an owner mark the
  //  generated response stale, so looking at the result you just computed - at a different distance -
  //  lit up "Generate Response" and made "Use" offer to regenerate.  They are still part of the
  //  tool's state, so the owner still needs them for its undo step.
  if( !m_restoringState )
    m_userChangedNoRegen.emit();
}//handleChartOptionChanged()


void MakeMcResponseForDrf::handlePrecisionChanged()
{
  m_customPrecision->setHidden( m_precision->currentIndex() != 4 );
  updateEstimate();

  if( !m_restoringState )
    m_userChanged.emit();
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
    const ceelo::GenerationOptions opts = generationOptions();
    const ceelo::ResponseGenerator::NodePlan plan = ceelo::ResponseGenerator::plan_nodes( gd, opts );

    // The measured throughput applies only to the geometry it was measured on; until the probe for
    //  this one lands, the ballpark model stands in (and the text says which it is).
    const McTimeCalibration *calib = nullptr;
    if( m_calibration && (m_calibration->geometry_key == gd.to_xml_string()) )
      calib = m_calibration.get();

    const vector<double> cum = prior_cumulative_cost( plan, opts, calib );
    const double total_s = cum.empty() ? 0.0 : cum.back();

    // "40 on-axis energies + 14 angles x 9 energies + 72 near-field positions x 9 energies"
    const int nE = static_cast<int>( plan.shape_energies_keV.size() );
    WString parts = WString::tr("mmr-est-part-backbone").arg( plan.n_backbone() );
    if( (plan.n_cos_theta > 0) && (nE > 0) )
    {
      WString angular;
      if( plan.transfer_anchors )
        angular = WString::tr("mmr-est-part-anchors").arg( plan.n_cos_theta ).arg( nE );
      else if( plan.n_phi > 1 )
        angular = WString::tr("mmr-est-part-angular-box").arg( plan.n_cos_theta ).arg( plan.n_phi ).arg( nE );
      else
        angular = WString::tr("mmr-est-part-angular").arg( plan.n_cos_theta ).arg( nE );
      parts = WString::tr("mmr-est-part-join").arg( parts ).arg( angular );
    }//if( angular nodes )

    if( (plan.n_near_positions > 0) && (nE > 0) )
    {
      const WString near_part = WString::tr("mmr-est-part-near").arg( plan.n_near_positions ).arg( nE );
      parts = WString::tr("mmr-est-part-join").arg( parts ).arg( near_part  );
    }//if( near-field nodes )

    const char *durKey = "mmr-est-dur-model";
    if( m_calibrating )
      durKey = "mmr-est-dur-measuring";
    else if( calib )
      durKey = calib->from_full_run ? "mmr-est-dur-from-run" : "mmr-est-dur-calibrated";
    const WString duration = WString::tr(durKey).arg( duration_phrase(total_s) );

    if( method == Method::QuickMc )
    {
      // Show the saving against a full General-profile characterization.
      ceelo::GenerationOptions full_opts = opts;
      full_opts.transfer_mode = false;
      full_opts.n_anchor_angles = 1;
      full_opts.profile = ceelo::ResponseProfile::General;
      const int full_nodes = ceelo::ResponseGenerator::estimated_node_count( gd, full_opts );
      m_estimate->setText( WString::tr("mmr-estimate-transfer-v2")
                            .arg( plan.total() ).arg( parts ).arg( duration ).arg( full_nodes ) );
    }else
    {
      m_estimate->setText( WString::tr("mmr-estimate-v2")
                            .arg( plan.total() ).arg( parts ).arg( duration ) );
    }
  }catch( std::exception & )
  {
    m_estimate->setText( "" );
  }
}//updateEstimate()


bool MakeMcResponseForDrf::generationRunning() const
{
  return m_generating;
}//generationRunning()


bool MakeMcResponseForDrf::startGeneration()
{
  // One run at a time.  Without this a second click (the dialog's footer button re-enables itself
  //  while a run is in flight, since `m_result` is null then) started a second full-core Monte
  //  Carlo, and the line below would replace the cancel flag the first run is watching - leaving it
  //  burning every core to completion, uncancellable, even after the window is closed.
  if( m_generating )
    return false;

  // Whatever the owner's edits currently say - see #setSeedProvider.
  refreshSeedFromProvider();

  ceelo::GeometryDescriptor gd;
  try
  {
    gd = m_geometry->toDescriptor();
  }catch( std::exception &e )
  {
    m_status->setText( WString::fromUTF8( e.what() ) );
    return false;
  }

  // A valid-but-guessed geometry (length fabricated from the diameter of a legacy DRF) is not
  //  enough to characterize - require a real, user-confirmed shape.
  if( !m_geometry->generationReady() )
  {
    m_status->setText( WString::tr("mmr-status-geom-incomplete") );
    return false;
  }

  ++m_generationId;
  const int generation_id = m_generationId;

  // A calibration probe in flight (or waiting) is moot once the real thing runs.
  ++m_calibrationId;
  m_calibrating = false;
  if( m_calibTimer )
    m_calibTimer->stop();

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
      anchor = transferAnchor( gd );
    }catch( std::exception &e )
    {
      m_status->setText( WString::tr("mmr-anchor-unusable")
                          .arg( WString::fromUTF8(e.what()) ) );
      return false;
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
          tool->handleGenerationFinished( response, errmsg, generation_id, nullptr );
        else
          wApp->enableUpdates( false );  //balance startGeneration's enableUpdates(true) regardless
        wApp->triggerUpdate();
      } );
    };//worker

    m_status->setText( WString::tr("mmr-status-transfer-building") );

    m_generating = true;   //handleGenerationFinished clears it, as for the MC methods
    wApp->enableUpdates( true );
    WServer::instance()->ioService().boost::asio::io_service::post( worker );
    return true;
  }//if( method == Method::CurveTransfer )

  ceelo::GenerationOptions opts = generationOptions();
  opts.detector_name = m_seedDrf ? m_seedDrf->name() : string("user geometry");
  opts.base_seed = 1;  //deterministic; re-running the same setup reproduces

  // Belt and braces with the m_generating guard above: signal whatever the previous run was
  //  watching before letting go of it, so no worker can ever be left without a way to be stopped.
  if( m_cancelFlag )
    m_cancelFlag->store( true );
  m_cancelFlag = make_shared<std::atomic<bool>>( false );
  opts.cancel = m_cancelFlag;

  // Grounding anchors are captured NOW (value copies) - nothing from the
  //  widget tree crosses into the worker thread.
  bool curve_derived = false;
  const vector<ceelo::GroundingPoint> ground_pts = groundToMeasured()
                      ? groundingPointsForDrf( m_seedDrf, gd, curve_derived )
                      : vector<ceelo::GroundingPoint>{};

  const string sessionId = wApp->sessionId();
  const string widgetId = id();

  // What the run will do, and what each node is predicted to cost: the prior the ETA refines from
  //  the measured rate as nodes land (see refreshProgressText).
  m_runGeometryKey = gd.to_xml_string();
  const ceelo::ResponseGenerator::NodePlan plan = ceelo::ResponseGenerator::plan_nodes( gd, opts );
  const McTimeCalibration *calib
      = (m_calibration && (m_calibration->geometry_key == m_runGeometryKey)) ? m_calibration.get() : nullptr;
  m_nodesTotal = plan.total();
  m_transferAnchors = plan.transfer_anchors;
  m_priorCumulative = prior_cumulative_cost( plan, opts, calib );
  m_generationStart = std::chrono::steady_clock::now();
  m_progressSnapshot = make_shared<McProgressSnapshot>();
  m_progressSnapshot->nodes_total = m_nodesTotal;

  // Per-node progress: the worker writes the snapshot (atomics only - nothing widget-side crosses
  //  the thread boundary; findById(...) on the session thread is the only widget access), and the
  //  2 s timer paints it.  A post is made only when the stage changes or the last node lands, so
  //  a stage transition or the finish never waits on the timer.
  const shared_ptr<McProgressSnapshot> snapshot = m_progressSnapshot;
  auto last_stage = make_shared<std::atomic<int>>( 0 );
  opts.node_progress = [sessionId,widgetId,generation_id,snapshot,last_stage]( const ceelo::NodeProgress &p ){
    snapshot->nodes_done = p.nodes_done;
    snapshot->nodes_total = p.nodes_total;
    const int stage = static_cast<int>( p.stage );
    snapshot->stage = stage;

    const bool last_node = (p.nodes_done >= p.nodes_total);
    if( (last_stage->exchange(stage) == stage) && !last_node )
      return;

    WServer::instance()->post( sessionId, [widgetId,generation_id](){
      auto *tool = dynamic_cast<MakeMcResponseForDrf *>( wApp->domRoot()->findById(widgetId) );
      if( tool )
        tool->updateProgress( generation_id );
      wApp->triggerUpdate();
    } );
  };

  auto worker = [gd,opts,ground_pts,curve_derived,sessionId,widgetId,generation_id](){
    shared_ptr<ceelo::DetectorResponse> response;
    string errmsg;

    // Per-node costs of the run: the best possible calibration for this geometry's next estimate.
    auto stats = make_shared<ceelo::GenerationStats>();
    try
    {
      ceelo::GenerationOptions run_opts = opts;
      run_opts.stats_out = stats.get();
      response = ceelo::ResponseGenerator::generate( gd, run_opts );
      if( response && !ground_pts.empty() )
      {
        // Model efficiencies at each point's own geometry, with the point distances in InterSpec's
        //  face-referenced convention (rather than leaving ground_to_points to interpret them).
        vector<ceelo::GroundingPoint> pts = ground_pts;
        for( ceelo::GroundingPoint &p : pts )
        {
          if( p.model_eff > 0.0 )
            continue;
          const double theta = std::acos( std::min( std::max( p.cos_theta, -1.0 ), 1.0 ) );
          const double phi = p.phi_deg * 3.14159265358979323846 / 180.0;
          const Eigen::Vector3d pos = CeeLoUtils::sourcePositionFromFace( gd, theta, phi, p.distance_cm );
          p.model_eff = response->eps_fep_at( p.energy_keV, pos ).value;
        }
        ceelo::ResponseGenerator::ground_to_points( *response, pts, curve_derived );
      }
    }catch( ceelo::GenerationCancelled & )
    {
      errmsg = "cancelled";
    }catch( std::exception &e )
    {
      errmsg = e.what();
    }

    WServer::instance()->post( sessionId, [widgetId,response,errmsg,generation_id,stats](){
      auto *tool = dynamic_cast<MakeMcResponseForDrf *>( wApp->domRoot()->findById(widgetId) );
      if( tool )
      {
        tool->handleGenerationFinished( response, errmsg, generation_id, stats );
      }else
      {
        cerr << "MakeMcResponseForDrf deleted while MC generation ran" << endl;
        wApp->enableUpdates( false );  //balance startGeneration's enableUpdates(true) regardless
      }
      wApp->triggerUpdate();
    } );
  };//worker

  m_generate->hide();
  m_cancelBtn->show();
  m_progress->setValue( 0.0 );
  m_progress->show();
  m_status->setText( WString::tr("mmr-progress-starting").arg( m_nodesTotal ) );
  m_generating = true;
  m_progressTimer->start();

  wApp->enableUpdates( true );

  WServer::instance()->ioService().boost::asio::io_service::post( worker );

  return true;
}//startGeneration()


void MakeMcResponseForDrf::cancelGeneration()
{
  if( m_cancelFlag )
    m_cancelFlag->store( true );
  if( m_progressTimer )
    m_progressTimer->stop();   //so "Cancelling..." stays up
  m_status->setText( WString::tr("mmr-status-cancelling") );
}//cancelGeneration()


void MakeMcResponseForDrf::updateProgress( const int generation_id )
{
  if( generation_id != m_generationId )
    return;  //stale run

  refreshProgressText();
}//updateProgress(...)


void MakeMcResponseForDrf::refreshProgressText()
{
  if( !m_generating || !m_progressSnapshot )
    return;

  // A cancel is in flight: leave "Cancelling..." alone.
  if( m_cancelFlag && m_cancelFlag->load() )
    return;

  const int n = m_progressSnapshot->nodes_done.load();
  const int N = std::max( 1, m_progressSnapshot->nodes_total.load() );
  const int stage = m_progressSnapshot->stage.load();

  const double elapsed = std::chrono::duration<double>( std::chrono::steady_clock::now()
                                                        - m_generationStart ).count();

  // ETA: the prior's remaining cost, scaled by how the finished nodes ran against their
  //  prediction - blended toward 1 with a pseudo-count so the first few (cheap, noisy) nodes do
  //  not swing it, and clamped so a single slow node cannot either.
  double total = elapsed;
  if( m_priorCumulative.size() > static_cast<size_t>( std::max(n, 0) ) )
  {
    const double predicted_done = m_priorCumulative[static_cast<size_t>( std::max(n, 0) )];
    const double predicted_total = m_priorCumulative.back();
    const double rho = ((n >= 3) && (predicted_done > 0.0)) ? (elapsed / predicted_done) : 1.0;
    const double rho_b = std::min( 5.0, std::max( 0.2, (n*rho + 8.0) / (n + 8.0) ) );
    total = elapsed + rho_b * std::max( 0.0, predicted_total - predicted_done );
  }//if( have a prior )

  m_progress->setValue( std::min( 1.0, double(n) / N ) );

  if( n <= 0 )
  {
    m_status->setText( WString::tr("mmr-progress-starting").arg( N ) );
    return;
  }

  const char *stageKey = "mmr-stage-finishing";
  if( n < N )
  {
    switch( stage )
    {
      case 1:  stageKey = "mmr-stage-backbone"; break;
      case 2:  stageKey = m_transferAnchors ? "mmr-stage-anchors" : "mmr-stage-angular"; break;
      case 3:  stageKey = "mmr-stage-near"; break;
      default: break;
    }
  }//if( n < N )

  m_status->setText( WString::tr("mmr-progress").arg( n ).arg( N ).arg( WString::tr(stageKey) )
                       .arg( minutes_str(elapsed) ).arg( minutes_str(total) ) );
}//refreshProgressText()


void MakeMcResponseForDrf::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  WContainerWidget::render( flags );

  // Deliberately NOT starting the timing probe here.  Being rendered does not mean being looked at:
  //  the Geom & MC tab is `ContentLoading::Eager` (so findById can see this tool), and
  //  WStackedWidget merely `setHidden()`s the pages that are not current - they stay in the tree and
  //  do get a Full render.  Probing here spent up to 3 CPU-seconds x every core for a tab the user
  //  never opened.  `DrfModifyWidget::handleTabSelected` calls scheduleTimeCalibration() when this
  //  tab is actually selected, and the standalone window does so from its own constructor.
  if( flags.test( Wt::RenderFlag::Full ) )
    m_shown = true;
}//render(...)


void MakeMcResponseForDrf::setDisabled( bool disabled )
{
  // `WWidget::enable()/disable()` are stateless slots Wt may pre-learn and then replay client-side;
  //  this override changes server state (it can start a timer and a Monte-Carlo probe), so it has
  //  to opt out - see the setDisabled() note in Wt/WWidget.h.
  isNotStateless();

  WContainerWidget::setDisabled( disabled );
  if( !disabled )
    scheduleTimeCalibration();  //Flat Disk -> Geometry Modeled: the estimate now matters
}//setDisabled(...)


ceelo::GenerationOptions MakeMcResponseForDrf::generationOptions() const
{
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

  if( selectedMethod() == Method::QuickMc )
  {
    // EFFTRAN-style transfer: MC only the on-axis energy backbone (plus a few
    //  forced cos-theta anchors when selected); the ray-traced kernel carries
    //  the distance/angle transfer.  generate() forces the FarField profile.
    opts.transfer_mode = true;
    opts.n_anchor_angles = (m_anchorAngles->currentIndex() == 0) ? 1 : 3;
  }

  return opts;
}//generationOptions()


void MakeMcResponseForDrf::scheduleTimeCalibration()
{
  // m_calibrating: a probe already running is a full-core Monte Carlo on the shared server thread
  //  pool; queueing more of them behind a burst of edits would just take cores from the session.
  if( !m_calibTimer || !m_shown || m_generating || m_calibrating || !isEnabled()
      || (selectedMethod() == Method::CurveTransfer) || !m_geometry->isValid() )
  {
    return;
  }

  // Already measured for exactly this geometry: nothing to do but show it.
  try
  {
    if( m_calibration && (m_calibration->geometry_key == m_geometry->toDescriptor().to_xml_string()) )
    {
      updateEstimate();
      return;
    }
  }catch( std::exception & )
  {
    return;
  }

  m_calibTimer->stop();
  m_calibTimer->start();  //restarting is the debounce
}//scheduleTimeCalibration()


void MakeMcResponseForDrf::startTimeCalibration()
{
  if( m_generating || !isEnabled() || (selectedMethod() == Method::CurveTransfer) )
    return;

  ceelo::GeometryDescriptor gd;
  try
  {
    gd = m_geometry->toDescriptor();
  }catch( std::exception & )
  {
    return;
  }

  const string key = gd.to_xml_string();
  ++m_calibrationId;
  const int calibration_id = m_calibrationId;
  m_calibrating = true;
  updateEstimate();  //says "timing a short test run..."

  const string sessionId = wApp->sessionId();
  const string widgetId = id();

  // One short MC node of THIS geometry: on axis at the far-field backbone distance, 662 keV (a
  //  little above the log-midpoint of the range, since cost rises with energy), loose precision,
  //  and hard event/CPU/wall caps so it is a second or three at most.  Only value copies and the
  //  widget id cross into the worker.
  auto worker = [gd,key,sessionId,widgetId,calibration_id](){
    McTimeCalibration calib;
    calib.geometry_key = key;

    try
    {
      std::vector<std::unique_ptr<ceelo::Material>> owned;
      ceelo::EfficiencyCalculator calc;
      ceelo::ResponseGenerator::configure_calculator( calc, gd, owned );

      const double a = gd.transverse_half_extent();
      calc.set_point_source( Eigen::Vector3d( 0.0, 0.0, -std::max( 10.0*a, 10.0 ) ) );

      ceelo::SimulationConfig cfg;
      cfg.energy_keV = 661.7;
      cfg.termination.target_fep_rel_precision = 0.03;
      cfg.termination.min_events = 20000;
      cfg.termination.max_events = 300000;
      cfg.termination.max_cpu_seconds = 3.0;
      cfg.termination.max_wall_seconds = 6.0;
      cfg.seed = 7;
      const unsigned threads = std::max( 1u, std::thread::hardware_concurrency() );
      cfg.batch_size = std::max<uint64_t>( 2000, 20000 / threads );

      const ceelo::EfficiencyResult r = calc.compute( cfg );
      const double eps = r.full_energy_peak_efficiency;
      if( (eps > 0.0) && (r.num_events_simulated > 0)
          && (r.cpu_time_seconds > 0.0) && (r.wall_time_seconds > 0.0) )
      {
        const double rel = r.fep_uncertainty / eps;
        calib.rel_var_per_event = rel * rel * double(r.num_events_simulated);
        calib.events_per_cpu_s = double(r.num_events_simulated) / r.cpu_time_seconds;
        calib.parallelism = std::max( 1.0, r.cpu_time_seconds / r.wall_time_seconds );
      }
    }catch( std::exception & )
    {
      //calib stays invalid: the model estimate stands
    }

    WServer::instance()->post( sessionId, [widgetId,calib,calibration_id](){
      auto *tool = dynamic_cast<MakeMcResponseForDrf *>( wApp->domRoot()->findById(widgetId) );
      if( tool )
        tool->handleTimeCalibrationFinished( calib, calibration_id );
      else
        wApp->enableUpdates( false );  //balance the enableUpdates(true) below
      wApp->triggerUpdate();
    } );
  };//worker

  wApp->enableUpdates( true );
  WServer::instance()->ioService().boost::asio::io_service::post( worker );
}//startTimeCalibration()


void MakeMcResponseForDrf::handleTimeCalibrationFinished( const McTimeCalibration &calib,
                                                          const int calibration_id )
{
  wApp->enableUpdates( false );  //FIRST - every started probe posts exactly one finish

  if( calibration_id != m_calibrationId )
    return;  //superseded by an edit, a run, or a restore

  m_calibrating = false;
  if( calib.valid() )
    m_calibration = make_shared<const McTimeCalibration>( calib );

  updateEstimate();
}//handleTimeCalibrationFinished(...)


void MakeMcResponseForDrf::handleGenerationFinished(
                              std::shared_ptr<ceelo::DetectorResponse> result,
                              const std::string &errmsg,
                              const int generation_id,
                              std::shared_ptr<const ceelo::GenerationStats> stats )
{
  // Balance the enableUpdates(true) from startGeneration FIRST - every started
  //  generation posts exactly one finish, including runs made stale by a
  //  restart or method switch (Wt refcounts the update lock, so a stale
  //  finish must still decrement it or server push stays on forever).
  wApp->enableUpdates( false );

  if( generation_id != m_generationId )
    return;  //stale run - a newer run/state owns the UI

  m_generating = false;
  if( m_progressTimer )
    m_progressTimer->stop();

  // The run's own per-node costs are the best calibration there is for this geometry - even from a
  //  cancelled run, once a few nodes are in.
  if( stats && (stats->nodes.size() >= 3) && (stats->total_cpu_s > 0.0) && (stats->total_wall_s > 0.0) )
  {
    double sum_var = 0.0;
    int n_var = 0;
    for( const ceelo::NodeStat &ns : stats->nodes )
    {
      if( (ns.fep_rel_prec > 0.0) && (ns.events > 0) )
      {
        sum_var += ns.fep_rel_prec * ns.fep_rel_prec * double(ns.events);
        ++n_var;
      }
    }//for( each node )

    if( n_var >= 3 )
    {
      auto calib = make_shared<McTimeCalibration>();
      calib->geometry_key = m_runGeometryKey;
      calib->rel_var_per_event = sum_var / n_var;
      calib->events_per_cpu_s = double(stats->total_events) / stats->total_cpu_s;
      calib->parallelism = std::max( 1.0, stats->total_cpu_s / stats->total_wall_s );
      calib->from_full_run = true;

      // Only a usable measurement: an unusable one would both mislabel the estimate ("from the
      //  last run") and, being keyed on this geometry, stop it ever being measured again.
      if( calib->valid() )
        m_calibration = calib;
    }//if( n_var >= 3 )
  }//if( have run statistics )

  m_generate->setHidden( m_hideGenerateButton || (selectedMethod() == Method::CurveTransfer) );
  m_generate->setEnabled( m_geometry->generationReady() );
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

    updateEstimate();  //a cancelled run still measured a few nodes
    m_userChanged.emit();
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
                    CeeLoUtils::faceDistanceFromCrystalOrigin( result->descriptor,
                                                               result->provenance.min_distance_cm )
                    * PhysicalUnits::cm );
    m_status->setText( WString::tr("mmr-status-transfer-done")
                        .arg( min_dist ) );
  }else
  {
    m_status->setText( WString::tr( grounded ? "mmr-status-done-grounded"
                                             : "mmr-status-done" ) );
  }

  updateEstimate();       //now "from the last run"; done before the emits, see below

  // Nothing below may touch `this`.  A handler of either signal can apply the response and accept
  //  the dialog this tool lives in, which (Wt4) destroys the whole window tree synchronously -
  //  the hazard `AuxWindow::emitReject` documents.  Wt's signal ring is itself deletion-safe, so
  //  the emit call completes; it is the *frame* that must not carry on through freed members.
  //  Only the widget id is held across the emit - never an observing_ptr - and it is resolved on
  //  this same session thread.
  const WidgetUtils::WidgetHandle self( this );

  m_responseGenerated.emit( result );

  if( !self.resolve_as<MakeMcResponseForDrf>() )
    return;  //a handler accepted the dialog and took this tool with it

  m_userChanged.emit();   //a new response is state an undo/redo owner needs to capture
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
    new_det = make_shared<DetectorPeakResponse>( "MC characterized detector",
                                          "Monte-Carlo parameterized response" );
  }//if( base ) / else

  // Sample the response into an ordinary efficiency curve whenever the DRF has
  //  none of its own - a geometry-only characterization, or no seed DRF at all.
  //  This used to be a flat `setIntrinsicEfficiencyFormula("1.0")` placeholder,
  //  which left a DRF that answered every EffEval query correctly through the
  //  attached response but was not valid, serializable or exportable.
  if( !new_det->isValid() )
  {
    try
    {
      CeeLoUtils::setLegacyEfficiencyFromResponse( *new_det, m_result );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("mmr-err-no-backbone").arg(e.what()),
                   WarningWidget::WarningMsgHigh );
      return;
    }
  }//if( !new_det->isValid() )

  // Record the geometry the user described, as well as the response built from it.  The response
  //  carries its own descriptor, but detaching it later (Modify Detector Response -> Flat Disk) must
  //  not leave the detector with no statement of what it physically is.
  try
  {
    if( m_geometry && m_geometry->generationReady() )
      new_det->setGeometry( make_shared<const ceelo::GeometryDescriptor>( m_geometry->toDescriptor() ) );
  }catch( std::exception & )
  {
    //an incomplete form; the response's own descriptor is still there
  }

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

  // The run-time estimate is refined by a short test Monte Carlo.  The embedded copy starts it when
  //  its tab is selected; this window is, by definition, already being looked at.
  m_tool->scheduleTimeCalibration();

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
