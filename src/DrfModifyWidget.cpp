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

#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include <Wt/WText.h>
#include <Wt/WMenu.h>
#include <Wt/WTable.h>
#include <Wt/WLabel.h>
#include <Wt/WCheckBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WTextArea.h>
#include <Wt/WMenuItem.h>
#include <Wt/WGroupBox.h>
#include <Wt/WPushButton.h>
#include <Wt/WGridLayout.h>
#include <Wt/WStackedWidget.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/SpecFile.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/DrfModifyWidget.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/MakeMcResponseForDrf.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace Wt;
using namespace std;


DrfModifyWidget::DrfModifyWidget( InterSpec *viewer,
                                  std::shared_ptr<const DetectorPeakResponse> drf,
                                  std::shared_ptr<const ceelo::GeometryDescriptor> geometry )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_orig( drf ),
    m_geometry( geometry ),
    m_tabMenu( nullptr ),
    m_tabStack( nullptr ),
    m_name( nullptr ),
    m_description( nullptr ),
    m_mcTool( nullptr ),
    m_fwhmTool( nullptr ),
    m_fwhmTabItem( nullptr ),
    m_fwhmApply( nullptr ),
    m_fwhmTabVisited( false ),
    m_bandTable( nullptr ),
    m_addBand( nullptr ),
    m_removeBand( nullptr ),
    m_anchorTable( nullptr ),
    m_addAnchor( nullptr ),
    m_removeAnchor( nullptr ),
    m_anchorRefDistance( nullptr ),
    m_anchorDefaultUncert( nullptr ),
    m_updatedDrf()
{
  assert( m_interspec );
  if( m_interspec )
    m_interspec->useMessageResourceBundle( "DrfModifyWidget" );

  addStyleClass( "DrfModifyWidget" );

  WGridLayout *layout = setLayout( make_unique<WGridLayout>() );
  layout->setContentsMargins( 0, 0, 0, 0 );

  auto stackOwned = make_unique<WStackedWidget>();
  m_tabStack = stackOwned.get();
  m_tabStack->addStyleClass( "DrfModifyStack" );

  auto menuOwned = make_unique<WMenu>( m_tabStack );
  m_tabMenu = menuOwned.get();
  m_tabMenu->addStyleClass( "VerticalNavMenu HeavyNavMenu HorizontalMenu" );
  layout->addWidget( std::move(menuOwned), 0, 0 );
  layout->addWidget( std::move(stackOwned), 1, 0 );
  layout->setRowStretch( 1, 1 );

  // --- Tab: Name & Description ---------------------------------------------
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    panel->addStyleClass( "DrfModifyPanel" );

    WContainerWidget *nameRow = panel->addNew<WContainerWidget>();
    nameRow->addStyleClass( "DrfModifyRow" );
    nameRow->addNew<WLabel>( WString::tr("Name") );
    m_name = nameRow->addNew<WLineEdit>();
    m_name->setTextSize( 32 );
    m_name->setText( WString::fromUTF8( m_orig ? m_orig->name() : string() ) );

    WContainerWidget *descRow = panel->addNew<WContainerWidget>();
    descRow->addStyleClass( "DrfModifyRow" );
    descRow->addNew<WLabel>( WString::tr("Description") );
    m_description = descRow->addNew<WTextArea>();
    m_description->setColumns( 40 );
    m_description->setRows( 3 );
    m_description->setText( WString::fromUTF8( m_orig ? m_orig->description() : string() ) );

    // Measured-curve anchor editor: only when seeded with a physical geometry
    //  (e.g. an ANGLE import) AND the DRF carries raw measured points.  This is
    //  the source-of-truth editor for the reference curve the response is
    //  anchored to; edits are folded into the working DRF on apply().
    const shared_ptr<const MeasuredDrfPoints> measured
        = m_orig ? m_orig->measuredPoints() : nullptr;
    if( m_geometry && measured && !measured->empty() )
    {
      WGroupBox *box = panel->addNew<WGroupBox>( WString::tr("dmw-anchor-title") );
      box->addStyleClass( "DrfModifyAnchor" );
      box->addNew<WText>( WString::tr("dmw-anchor-help") )->setInline( false );

      // Reference distance the curve is anchored at (cm).
      const vector<MeasuredEffPoint> &pts = measured->points();
      double refDistCm = 0.0;
      for( const MeasuredEffPoint &p : pts )
      {
        if( p.distance > 0.0f )
        {
          refDistCm = p.distance / PhysicalUnits::cm;
          break;
        }
      }//for( find a positive reference distance )
      if( refDistCm <= 0.0 )
        refDistCm = m_orig->absoluteEfficiencyDistance() / PhysicalUnits::cm;

      WContainerWidget *distRow = box->addNew<WContainerWidget>();
      distRow->addStyleClass( "DrfModifyRow" );
      distRow->addNew<WLabel>( WString::tr("dmw-anchor-ref-dist") );
      m_anchorRefDistance = distRow->addNew<WLineEdit>();
      m_anchorRefDistance->setTextSize( 8 );
      if( refDistCm > 0.0 )
      {
        char buf[32];
        snprintf( buf, sizeof(buf), "%.4g", refDistCm );
        m_anchorRefDistance->setText( buf );
      }
      distRow->addNew<WLabel>( WString::tr("dmw-anchor-cm") );

      WContainerWidget *defRow = box->addNew<WContainerWidget>();
      defRow->addStyleClass( "DrfModifyRow" );
      defRow->addNew<WLabel>( WString::tr("dmw-anchor-default-uncert") );
      m_anchorDefaultUncert = defRow->addNew<WLineEdit>();
      m_anchorDefaultUncert->setTextSize( 6 );
      m_anchorDefaultUncert->setText( "5" );
      defRow->addNew<WLabel>( WString::tr("dmw-anchor-percent") );

      m_anchorTable = box->addNew<WTable>();
      m_anchorTable->addStyleClass( "DrfModifyAnchorTable" );
      m_anchorTable->elementAt(0,0)->addNew<WText>( WString::tr("dmw-anchor-energy") );
      m_anchorTable->elementAt(0,1)->addNew<WText>( WString::tr("dmw-anchor-eff") );
      m_anchorTable->elementAt(0,2)->addNew<WText>( WString::tr("dmw-anchor-uncert") );

      WContainerWidget *btns = box->addNew<WContainerWidget>();
      m_addAnchor = btns->addNew<WPushButton>( WString::tr("dmw-anchor-add") );
      m_addAnchor->addStyleClass( "LinkBtn" );
      m_addAnchor->clicked().connect( this, [this](){ addAnchorRow( 0.0f, 0.0f, 0.0f ); } );
      m_removeAnchor = btns->addNew<WPushButton>( WString::tr("dmw-anchor-remove") );
      m_removeAnchor->addStyleClass( "LinkBtn" );
      m_removeAnchor->clicked().connect( this, &DrfModifyWidget::removeAnchorRow );

      for( const MeasuredEffPoint &p : pts )
        addAnchorRow( p.energy, p.efficiency, p.fracStatUncert );
      m_removeAnchor->setEnabled( !m_anchors.empty() );
    }//if( seeded with geometry + measured points )

    m_tabMenu->addItem( WString::tr("dmw-tab-name"), std::move(panelOwned) );
  }

  // --- Tab: Geometry & MC characterization ---------------------------------
  {
    auto toolOwned = make_unique<MakeMcResponseForDrf>( m_interspec, m_orig, m_geometry );
    m_mcTool = toolOwned.get();

    // ContentLoading::Eager: this tool posts work to a worker thread and gets back to itself with
    //  `findById(...)`; a Lazy tab parks its contents in `WMenuItem::uContents_`, outside the
    //  widget tree, where `findById` can not see it - see the FWHM tab below.
    m_tabMenu->addItem( WString::tr("dmw-tab-geometry"), std::move(toolOwned),
                        ContentLoading::Eager );
  }

  // --- Tab: FWHM -----------------------------------------------------------
  {
    // MakeFwhmForDrf fits FWHM from the foreground's peaks; it (and its auto
    //  peak-search) require a loaded foreground spectrum.  Only build it when
    //  one is present - otherwise show a placeholder, so the rest of the
    //  Modify dialog (name, geometry/MC, uncertainty) still works with no
    //  spectrum loaded.
    const shared_ptr<const SpecUtils::Measurement> foreground
        = m_interspec ? m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground )
                      : nullptr;
    if( foreground )
    {
      auto panelOwned = make_unique<WContainerWidget>();
      WContainerWidget *panel = panelOwned.get();
      WGridLayout *fwhmLayout = panel->setLayout( make_unique<WGridLayout>() );
      fwhmLayout->setContentsMargins( 0, 0, 0, 0 );

      // MakeFwhmForDrf edits a (non-const) DRF; give it a private clone - we
      //  only read its fitted coefficients back out on apply().
      shared_ptr<DetectorPeakResponse> fwhm_seed
          = m_orig ? make_shared<DetectorPeakResponse>( *m_orig ) : nullptr;

      // `auto_fit_peaks == false`: the automated peak search is a multi-threaded fit of the whole
      //  spectrum, so dont pay for it every time this dialog is opened - #handleTabSelected kicks
      //  it off the first time the user actually looks at this tab.
      auto toolOwned = make_unique<MakeFwhmForDrf>( false, m_interspec, fwhm_seed );
      m_fwhmTool = toolOwned.get();
      fwhmLayout->addWidget( std::move(toolOwned), 0, 0 );
      fwhmLayout->setRowStretch( 0, 1 );

      // Looking at this tab must not silently cost the user the FWHM their detector already has,
      //  so replacing it is opt-in.  Starts un-checked even for a DRF with no FWHM at all, so that
      //  what "Use" does never depends on whether the user happened to visit this tab (the fit
      //  differs before and after the automated peak search): #handleTabSelected checks it, once,
      //  when the tab is first opened and the DRF has no FWHM to lose.  #apply honors this.
      m_fwhmApply = fwhmLayout->addWidget( make_unique<WCheckBox>( WString::tr("dmw-fwhm-apply-cb") ), 1, 0 );
      m_fwhmApply->addStyleClass( "DrfModifyFwhmApply" );
      HelpSystem::attachToolTipOn( m_fwhmApply, WString::tr("dmw-fwhm-apply-cb-tt"), true );

      // ContentLoading::Eager, so the tool is in the widget tree (and hence resolvable by
      //  `findById`) while its peak search runs; with the default Lazy policy the widget sits in
      //  `WMenuItem::uContents_` with no widget parent, the search completion silently no-ops, and
      //  the tab shows "Currently searching for peaks..." forever.
      m_fwhmTabItem = m_tabMenu->addItem( WString::tr("dmw-tab-fwhm"), std::move(panelOwned),
                                          ContentLoading::Eager );
    }else
    {
      auto placeholder = make_unique<WContainerWidget>();
      placeholder->addStyleClass( "DrfModifyPanel" );
      placeholder->addNew<WText>( WString::tr("dmw-fwhm-no-foreground") )->setInline( false );
      m_tabMenu->addItem( WString::tr("dmw-tab-fwhm"), std::move(placeholder) );
    }
  }

  // --- Tab: Baseline uncertainty (energy-range bands) ----------------------
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    panel->addStyleClass( "DrfModifyPanel DrfModifyUncert" );

    panel->addNew<WText>( WString::tr("dmw-uncert-help") )->setInline( false );

    m_bandTable = panel->addNew<WTable>();
    m_bandTable->addStyleClass( "DrfModifyBandTable" );
    m_bandTable->elementAt(0,0)->addNew<WText>( WString::tr("dmw-band-lower") );
    m_bandTable->elementAt(0,1)->addNew<WText>( WString::tr("dmw-band-upper") );
    m_bandTable->elementAt(0,2)->addNew<WText>( WString::tr("dmw-band-frac") );

    WContainerWidget *btns = panel->addNew<WContainerWidget>();
    m_addBand = btns->addNew<WPushButton>( WString::tr("dmw-add-band") );
    m_addBand->addStyleClass( "LinkBtn" );
    m_addBand->clicked().connect( this, [this](){ addBandRow( 0.0f, 0.0f, 0.05f ); } );
    m_removeBand = btns->addNew<WPushButton>( WString::tr("dmw-remove-band") );
    m_removeBand->addStyleClass( "LinkBtn" );
    m_removeBand->clicked().connect( this, &DrfModifyWidget::removeBandRow );

    // Seed rows from any existing piecewise-band uncertainty.
    const shared_ptr<const DetectorEfficiencyUncert> uncert
        = m_orig ? m_orig->efficiencyUncert() : nullptr;
    if( uncert )
    {
      for( const EffUncertBand &b : uncert->bands() )
        addBandRow( b.lowerEnergy, b.upperEnergy, b.fractionalUncert );
    }
    m_removeBand->setEnabled( !m_bands.empty() );

    m_tabMenu->addItem( WString::tr("dmw-tab-uncert"), std::move(panelOwned) );
  }

  m_tabMenu->itemSelected().connect( this, &DrfModifyWidget::handleTabSelected );

  m_tabMenu->select( 0 );

  WText *note = new WText( WString::tr("dmw-export-note") );
  note->addStyleClass( "DrfModifyNote" );
  note->setInline( false );
  layout->addWidget( std::unique_ptr<WText>(note), 2, 0 );
}//DrfModifyWidget constructor


DrfModifyWidget::~DrfModifyWidget()
{
}


void DrfModifyWidget::handleTabSelected( Wt::WMenuItem *item )
{
  // Only the first opening of the FWHM tab does anything; `WMenu::itemSelected()` also fires for
  //  the initial `select(0)`, and re-visiting the tab must not re-search or re-tick the checkbox.
  if( m_fwhmTabVisited || !item || (item != m_fwhmTabItem) || !m_fwhmTool )
    return;

  m_fwhmTabVisited = true;

  // Nothing to lose -> offer the fit by default, now that the user is looking at it.
  if( m_fwhmApply && (!m_orig || !m_orig->hasResolutionInfo()) )
    m_fwhmApply->setChecked( true );

  m_fwhmTool->startAutomatedPeakSearch();
}//void handleTabSelected( Wt::WMenuItem *item )


Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &DrfModifyWidget::updatedDrf()
{
  return m_updatedDrf;
}


void DrfModifyWidget::addBandRow( const float lowerEnergy, const float upperEnergy,
                                  const float fracUncert )
{
  const int row = m_bandTable->rowCount();  //row 0 is the header

  auto make_edit = [this,row]( const int col, const float value ) -> WLineEdit * {
    WLineEdit *edit = m_bandTable->elementAt(row,col)->addNew<WLineEdit>();
    edit->setTextSize( 8 );
    if( value != 0.0f )
    {
      char buf[32];
      snprintf( buf, sizeof(buf), "%.4g", value );
      edit->setText( buf );
    }
    return edit;
  };

  BandRow r;
  r.lower = make_edit( 0, lowerEnergy );
  r.upper = make_edit( 1, upperEnergy );
  r.frac  = make_edit( 2, fracUncert );
  m_bands.push_back( r );
  m_removeBand->setEnabled( !m_bands.empty() );
}//addBandRow(...)


void DrfModifyWidget::removeBandRow()
{
  if( m_bands.empty() )
    return;
  m_bandTable->removeRow( m_bandTable->rowCount() - 1 );
  m_bands.pop_back();
  m_removeBand->setEnabled( !m_bands.empty() );
}//removeBandRow()


void DrfModifyWidget::addAnchorRow( const float energy, const float efficiency,
                                    const float fracUncert )
{
  const int row = m_anchorTable->rowCount();  //row 0 is the header

  auto make_edit = [this,row]( const int col, const float value, const bool blankZero ) -> WLineEdit * {
    WLineEdit *edit = m_anchorTable->elementAt(row,col)->addNew<WLineEdit>();
    edit->setTextSize( 10 );
    if( !blankZero || (value != 0.0f) )
    {
      char buf[32];
      snprintf( buf, sizeof(buf), "%.6g", value );
      edit->setText( buf );
    }
    return edit;
  };

  AnchorRow r;
  r.energy = make_edit( 0, energy, true );
  r.eff    = make_edit( 1, efficiency, true );
  // Store the per-point uncertainty as a percentage; leave blank to inherit the
  //  default uncert % on apply.
  r.uncert = make_edit( 2, 100.0f*fracUncert, true );
  m_anchors.push_back( r );
  if( m_removeAnchor )
    m_removeAnchor->setEnabled( !m_anchors.empty() );
}//addAnchorRow(...)


void DrfModifyWidget::removeAnchorRow()
{
  if( m_anchors.empty() )
    return;
  m_anchorTable->removeRow( m_anchorTable->rowCount() - 1 );
  m_anchors.pop_back();
  m_removeAnchor->setEnabled( !m_anchors.empty() );
}//removeAnchorRow()


void DrfModifyWidget::applyAnchorEdits( DetectorPeakResponse &working )
{
  if( !m_anchorTable )
    return;  //anchor editor was not built

  // Default fractional uncertainty for blank per-point cells.
  float defaultFrac = 0.0f;
  {
    const string ds = m_anchorDefaultUncert->text().toUTF8();
    try{ defaultFrac = 0.01f * std::stof( ds ); }catch( std::exception & ){ defaultFrac = 0.0f; }
    if( defaultFrac < 0.0f )
      defaultFrac = 0.0f;
  }

  // Reference distance (cm) the curve is anchored at.
  double refDistCm = 0.0;
  {
    const string rs = m_anchorRefDistance->text().toUTF8();
    try{ refDistCm = std::stod( rs ); }catch( std::exception & ){ refDistCm = 0.0; }
  }
  if( refDistCm <= 0.0 )
    return;  //can't build absolute efficiency without a valid distance

  const double refDist = refDistCm * PhysicalUnits::cm;

  vector<DetectorPeakResponse::EnergyEffPoint> effpts;
  vector<MeasuredEffPoint> measpts;
  bool anyUncert = false;
  for( const AnchorRow &r : m_anchors )
  {
    const string es = r.energy->text().toUTF8();
    const string fs = r.eff->text().toUTF8();
    const string us = r.uncert->text().toUTF8();
    if( es.empty() && fs.empty() && us.empty() )
      continue;  //blank row

    float energy = 0.0f, eff = 0.0f;
    try
    {
      energy = std::stof( es );
      eff = std::stof( fs );
    }catch( std::exception & )
    {
      continue;  //skip malformed rows
    }
    if( (energy <= 0.0f) || (eff <= 0.0f) )
      continue;

    float frac = defaultFrac;
    if( !us.empty() )
    {
      try{ frac = 0.01f * std::stof( us ); }catch( std::exception & ){ frac = defaultFrac; }
    }
    if( frac < 0.0f )
      frac = 0.0f;

    DetectorPeakResponse::EnergyEffPoint e;
    e.energy = energy;
    e.efficiency = eff;
    if( frac > 0.0f )
    {
      e.efficiencyUncert = eff * frac;
      anyUncert = true;
    }
    effpts.push_back( e );

    MeasuredEffPoint m;
    m.energy = energy;
    m.efficiency = eff;
    m.fracStatUncert = frac;
    m.distance = static_cast<float>( refDist );
    measpts.push_back( m );
  }//for( const AnchorRow &r : m_anchors )

  if( effpts.size() < 2 )
    return;  //need at least two points for a curve; leave the DRF untouched

  (void)anyUncert;

  const float diameter = working.detectorDiameter();
  if( diameter <= 0.0f )
    return;  //FarFieldAbsolute needs a positive diameter; leave the DRF as-is

  try
  {
    working.setEfficiencyPoints( effpts, diameter, refDist,
                          DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

    auto meas = make_shared<MeasuredDrfPoints>();
    meas->setPoints( measpts );
    working.setMeasuredPoints( meas );
  }catch( std::exception & )
  {
    //Malformed edits (e.g. duplicate/degenerate energies): keep the seed curve.
  }
}//applyAnchorEdits(...)


void DrfModifyWidget::apply()
{
  // One working copy that every tab writes onto.
  shared_ptr<DetectorPeakResponse> working = m_orig
      ? make_shared<DetectorPeakResponse>( *m_orig )
      : make_shared<DetectorPeakResponse>();
  if( m_orig )
    working->setParentHashValue( m_orig->hashValue() );

  // Name / description.
  const string name = m_name->text().toUTF8();
  if( !name.empty() )
    working->setName( name );
  working->setDescription( m_description->text().toUTF8() );

  // Measured-curve anchor edits (energy/efficiency/uncert + reference distance):
  //  rebuild the far-field absolute efficiency and measured points before the
  //  MC/uncertainty steps read them.
  applyAnchorEdits( *working );

  // FWHM: pull the fitted coefficients (if any) without triggering the tool's
  //  own detector-changed emit.  Only when the user asked for it - see where the
  //  checkbox is created; un-checked leaves whatever FWHM the DRF came in with.
  const bool useFwhmFit = m_fwhmApply && m_fwhmApply->isChecked();
  const shared_ptr<MakeFwhmForDrf::ToolState> fwhm
      = (useFwhmFit && m_fwhmTool) ? m_fwhmTool->currentState() : nullptr;
  if( fwhm && !fwhm->m_parameters.empty() && m_fwhmTool->isValidFwhm() )
  {
    const DetectorPeakResponse::ResolutionFnctForm form
        = DetectorPeakResponse::ResolutionFnctForm( std::max( 0, fwhm->m_fwhm_index ) );

    // `setFwhmCoefficients` throws if the coefficients dont match the forms arity - dont let that
    //  escape into a signal handler and take the rest of the users edits with it.
    try
    {
      working->setFwhmCoefficients( fwhm->m_parameters, form );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("dmw-err-fwhm-not-applied").arg(e.what()),
                   WarningWidget::WarningMsgHigh );
    }
  }else if( useFwhmFit )
  {
    // The user asked for the fit, but there isnt a usable one (fit failed, or the peak search is
    //  still running) - say so rather than closing as if it had been applied.
    passMessage( WString::tr("dmw-err-no-fwhm-fit"), WarningWidget::WarningMsgHigh );
  }//if( have a valid FWHM fit ) / else if( the user wanted one )

  // Monte-Carlo response: attach + ground (to the DRF's measured points if it
  //  has them; otherwise the MC tool already grounded to a sampled curve).
  const shared_ptr<const ceelo::DetectorResponse> resp = m_mcTool->generatedResponse();
  if( resp )
    working->setCeeloResponse( resp );

  // Baseline uncertainty bands: keep any existing node covariance, replace the
  //  piecewise bands with the edited set.
  {
    vector<EffUncertBand> bands;
    for( const BandRow &r : m_bands )
    {
      float lo = 0.0f, hi = 0.0f, fr = 0.0f;
      const string ls = r.lower->text().toUTF8();
      const string us = r.upper->text().toUTF8();
      const string fs = r.frac->text().toUTF8();
      if( ls.empty() && us.empty() && fs.empty() )
        continue;  //blank row
      try
      {
        lo = std::stof( ls );
        hi = std::stof( us );
        fr = std::stof( fs );
      }catch( std::exception & )
      {
        continue;  //skip malformed rows rather than fail the whole apply
      }
      if( (hi > lo) && (fr >= 0.0f) )
        bands.push_back( EffUncertBand{ lo, hi, fr } );
    }//for( const BandRow &r : m_bands )

    std::stable_sort( begin(bands), end(bands),
      []( const EffUncertBand &a, const EffUncertBand &b ){ return a.lowerEnergy < b.lowerEnergy; } );

    const shared_ptr<const DetectorEfficiencyUncert> existing = working->efficiencyUncert();
    if( !bands.empty() || (existing && existing->hasBands()) )
    {
      auto uncert = existing ? make_shared<DetectorEfficiencyUncert>( *existing )
                             : make_shared<DetectorEfficiencyUncert>();
      uncert->setBands( bands );  //empty clears the bands, keeps node covariance
      working->setEfficiencyUncert( uncert->isEmpty() ? nullptr : uncert );
    }
  }//baseline uncertainty bands

  m_updatedDrf.emit( working );
}//apply()


// ---------------------------------------------------------------------------

DrfModifyWidget *DrfModifyWindow::tool()
{
  return m_tool;
}


DrfModifyWindow::DrfModifyWindow( InterSpec *viewer,
                                  std::shared_ptr<const DetectorPeakResponse> drf,
                                  std::shared_ptr<const ceelo::GeometryDescriptor> geometry )
 : AuxWindow( WString::tr("window-title-modify-drf"),
             (AuxWindowProperties::TabletNotFullScreen
              | AuxWindowProperties::SetCloseable
              | AuxWindowProperties::DisableCollapse
              | AuxWindowProperties::EnableResize
              | AuxWindowProperties::IsModal) ),
  m_tool( nullptr )
{
  assert( viewer );
  if( !viewer )
    return;

  const int ww = viewer->renderedWidth();
  const int wh = viewer->renderedHeight();
  if( ww > 100 && wh > 100 )
  {
    const int width = std::min( (8*ww)/9, 820 );
    const int height = std::min( 760, ((wh < 420) ? wh : (19*wh)/20) );
    resizeWindow( width, height );
    setMinimumSize( std::min(width,600), std::min(height,450) );
  }//if( ww > 100 && wh > 100 )

  {
    auto toolOwned = make_unique<DrfModifyWidget>( viewer, drf, geometry );
    m_tool = toolOwned.get();
    stretcher()->addWidget( std::move(toolOwned), 0, 0 );
  }
  stretcher()->setContentsMargins( 0, 0, 0, 0 );

  AuxWindow::addHelpInFooter( footer(), "modify-drf" );

  WPushButton *cancel = addCloseButtonToFooter( WString::tr("Cancel") );
  cancel->clicked().connect( this, &AuxWindow::hide );

  WPushButton *use = footer()->addNew<WPushButton>( WString::tr("dmw-use-btn") );
  use->clicked().connect( m_tool, &DrfModifyWidget::apply );

  show();
  resizeToFitOnScreen();
  centerWindow();
}//DrfModifyWindow constructor
