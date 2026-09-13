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
#include <numeric>
#include <functional>
#include <vector>
#include <algorithm>

#include <Wt/WText.h>
#include <Wt/WMenu.h>
#include <Wt/WTable.h>
#include <Wt/WLabel.h>
#include <Wt/WCheckBox.h>
#include <Wt/WComboBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WTextArea.h>
#include <Wt/WMenuItem.h>
#include <Wt/WGroupBox.h>
#include <Wt/WTableCell.h>
#include <Wt/WPushButton.h>
#include <Wt/WGridLayout.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include <Wt/WStackedWidget.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/SpecFile.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DrfModifyWidget.h"
#include "InterSpec/EccUncertOptions.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MakeMcResponseForDrf.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace Wt;
using namespace std;

namespace
{
/** A click landing on the body of a menu item, but outside its anchor, does not reach Wt's own
 wiring, so the tab does not change - which with the roomy `.SideMenu` padding is most of the
 button.  Connecting the item's `clicked()` to `WMenu::select(...)` covers those clicks.

 No explicit `triggered()`/`itemSelected()` emit: `WMenu::select()` emits when the selected index
 actually changes, and nothing when it does not, so doing it here would just double-fire.
 */
void make_item_selectable( WMenu *menu, WMenuItem *item )
{
  assert( menu && item );
  if( !menu || !item )
    return;

  item->clicked().connect( menu, [menu,item](){ menu->select( item ); } );
}//void make_item_selectable( WMenu *menu, WMenuItem *item )


/** The uncorrelated ("Stat. %") and correlated ("Corr. %") fractional uncertainties to seed a point
 row at `energy_kev` with.

 Prefers the split the covariance was actually built from, matching the node by energy rather than
 by index (only points that carried an uncertainty become nodes, so the node set can be a subset of
 the points).  Failing that - a covariance set directly, or one restored from a URL, which does not
 carry the split - the whole 1-sigma envelope is reported as correlated, which rebuilds the same
 covariance under a fully-correlated model.
 */
void uncert_for_energy( const shared_ptr<const DetectorEfficiencyUncert> &uncert,
                        const float energy_kev, float &statFrac, float &corrFrac )
{
  statFrac = corrFrac = 0.0f;

  if( !uncert || !uncert->hasNodeCovariance() || (energy_kev <= 0.0f) )
    return;

  if( uncert->hasComponentSplit() )
  {
    const vector<float> &nodes = uncert->covarianceEnergies();
    const vector<float> &corr = uncert->correlatedComponent();
    const vector<float> &uncorr = uncert->uncorrelatedComponent();

    for( size_t i = 0; i < nodes.size(); ++i )
    {
      // Node energies come from the same floats the rows do, so an exact match is the norm; the
      //  relative tolerance just absorbs a serialization round-trip.
      if( fabs(nodes[i] - energy_kev) <= (1.0E-5f * std::max(nodes[i], energy_kev)) )
      {
        corrFrac = (i < corr.size()) ? corr[i] : 0.0f;
        statFrac = (i < uncorr.size()) ? uncorr[i] : 0.0f;
        return;
      }
    }//for( size_t i = 0; i < nodes.size(); ++i )
  }//if( uncert->hasComponentSplit() )

  const vector<double> sig = uncert->fracUncertainties( { static_cast<double>(energy_kev) } );
  if( sig.size() == 1 )
    corrFrac = static_cast<float>( sig[0] );
}//uncert_for_energy(...)
}//namespace


DrfModifyWidget::DrfModifyWidget( InterSpec *viewer,
                                  std::shared_ptr<const DetectorPeakResponse> drf )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_orig( drf ),
    m_geometry( drf ? drf->geometry() : nullptr ),
    m_tabMenu( nullptr ),
    m_tabStack( nullptr ),
    m_name( nullptr ),
    m_description( nullptr ),
    m_mcTool( nullptr ),
    m_fwhmTool( nullptr ),
    m_fwhmTabItem( nullptr ),
    m_modeToggle( nullptr ),
    m_geometryModeled( false ),
    m_origHasPoints( false ),
    m_anchorHelp( nullptr ),
    m_pointsEditor( nullptr ),
    m_coefEditor( nullptr ),
    m_formulaEditor( nullptr ),
    m_anchorTable( nullptr ),
    m_anchorTableWrap( nullptr ),
    m_addAnchor( nullptr ),
    m_removeAnchor( nullptr ),
    m_anchorRefDistance( nullptr ),
    m_anchorDefaultUncert( nullptr ),
    m_anchors(),
    m_uncertOptions( nullptr ),
    m_anchorHasSourceCol( false ),
    m_anchorEnergyUnits( static_cast<float>(PhysicalUnits::keV) ),
    m_anchorIsAbsolute( false ),
    m_coefEdits(),
    m_coefParams( nullptr ),
    m_covTable( nullptr ),
    m_addCoef( nullptr ),
    m_removeCoef( nullptr ),
    m_coefCovMatrix(),
    m_formulaText( nullptr ),
    m_effEnergyUnits( nullptr ),
    m_effUnitsRow( nullptr ),
    m_seedState( nullptr ),
    m_generateBtn( nullptr ),
    m_changedSinceGenerate( false ),
    m_applyAfterGenerate( false ),
    m_suppressNextEditMark( false ),
    m_updatedDrf(),
    m_renderFlags(),
    m_currentState( nullptr ),
    m_restoringState( false )
{
  assert( m_interspec );
  if( m_interspec )
    m_interspec->useMessageResourceBundle( "DrfModifyWidget" );

  addStyleClass( "DrfModifyWidget" );

  // Own our stylesheet, rather than relying on the MakeMcResponseForDrf we happen to always create.
  wApp->useStyleSheet( "InterSpec_resources/MakeMcResponseForDrf.css" );

  // A side menu needs horizontal room; when there isnt any (phone, or just a narrow window), put
  //  the sections in a row of tabs across the top instead - same as DrfSelect does.
  int screen_width = m_interspec ? m_interspec->renderedWidth() : 0;
  if( (screen_width < 100) && m_interspec && m_interspec->isMobile() )
    screen_width = wApp->environment().screenWidth();  //not rendered yet - best guess
  const bool narrow_layout = ((screen_width > 100) && (screen_width < 600));

  WGridLayout *layout = setLayout( make_unique<WGridLayout>() );
  layout->setContentsMargins( 0, 0, 0, 0 );

  auto stackOwned = make_unique<WStackedWidget>();
  m_tabStack = stackOwned.get();
  m_tabStack->addStyleClass( "DrfModifyStack" );
  if( !narrow_layout )
    m_tabStack->addStyleClass( "UseInfoStack" );  //divider between the side menu and the content

  // Wt4's WStackedWidget ctor sets an inline `overflow:hidden` that beats the stylesheet, so the
  //  scrolling has to be set through the API - see the same note in DrfSelect.cpp.
  m_tabStack->setOverflow( Overflow::Auto, Wt::Orientation::Vertical );

  auto menuOwned = make_unique<WMenu>( m_tabStack );
  m_tabMenu = menuOwned.get();

  // The menu goes in its own container so the menu itself, and not the whole dialog, is what
  //  scrolls when the sections dont all fit.
  auto menuHolderOwned = make_unique<WContainerWidget>();
  WContainerWidget *menuHolder = menuHolderOwned.get();
  menuHolder->addWidget( std::move(menuOwned) );

  if( narrow_layout )
  {
    m_tabMenu->addStyleClass( "VerticalNavMenu HorizontalMenu HeavyNavMenu DrfModifyMenuHorizontal" );
    menuHolder->setOverflow( Overflow::Auto, Wt::Orientation::Horizontal );
    menuHolder->setOverflow( Overflow::Hidden, Wt::Orientation::Vertical );

    layout->addWidget( std::move(menuHolderOwned), 0, 0 );
    layout->addWidget( std::move(stackOwned), 1, 0 );
    layout->setRowStretch( 1, 1 );
  }else
  {
    m_tabMenu->addStyleClass( "VerticalNavMenu SideMenu HeavyNavMenu DrfModifyMenu" );
    menuHolder->setOverflow( Overflow::Auto, Wt::Orientation::Vertical );
    menuHolder->setOverflow( Overflow::Hidden, Wt::Orientation::Horizontal );

    layout->addWidget( std::move(menuHolderOwned), 0, 0 );
    layout->addWidget( std::move(stackOwned), 0, 1 );
    layout->setColumnStretch( 1, 1 );
    layout->setRowStretch( 0, 1 );
  }//if( narrow_layout ) / else

  // What kind of efficiency the DRF carries decides which Anchor editor is the right one, and
  //  (for far-field) whether it starts Flat Disk or Geometry Modeled:
  //
  //   - raw `measuredPoints` (an ANGLE / Make-Detector-Response import): absolute efficiencies at a
  //     reference distance, with per-point statistical + certificate uncertainties;
  //   - failing that, an energy/efficiency-pair efficiency curve (a GADRAS Efficiency.csv): the same
  //     measured table, just without any uncertainties, and intrinsic rather than absolute - so no
  //     reference distance applies.
  const shared_ptr<const MeasuredDrfPoints> measured
      = m_orig ? m_orig->measuredPoints() : nullptr;
  const bool have_measured = (measured && !measured->empty());

  const shared_ptr<const DetectorEfficiencyCurve> eff_curve
      = m_orig ? m_orig->efficiencyCurve() : nullptr;
  const bool have_pairs = (!have_measured && eff_curve
                  && (eff_curve->form() == DetectorPeakResponse::kEnergyEfficiencyPairs)
                  && (eff_curve->energyEfficiencies().size() >= 2));

  m_origHasPoints = (have_measured || have_pairs);
  m_anchorIsAbsolute = (m_orig && (m_orig->geometryType()
                                   == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute));

  // Two things that used to share m_anchorIsAbsolute.  A source column only makes sense for an
  //  absolute reference curve (its certificate uncertainty is blocked per source); the correlated
  //  column itself is shown for every curve, governed instead by m_uncertOptions.
  m_anchorHasSourceCol = m_anchorIsAbsolute;

  // The Energy column of a pairs curve is in the curve's own units - which a GADRAS CSV may set to
  //  MeV.  A formula curve's rows are covariance nodes instead, and those are keV by contract, so
  //  the column stays keV there no matter what units the formula itself is written in.
  //  Measured points are keV by contract (MeasuredEffPoint::energy) no matter what units the
  //  curve's equation uses - a MakeDrf detector can be an MeV equation fitted to keV points - so
  //  only a pairs curve, whose own numbers fill the column, sets this from the curve.
  const bool formula_form = (eff_curve && eff_curve->isValid()
                     && (eff_curve->form() == DetectorPeakResponse::kFunctialEfficienyForm));
  if( !formula_form && !have_measured && eff_curve && eff_curve->isValid()
      && (eff_curve->energyUnits() > 0.0f) )
    m_anchorEnergyUnits = eff_curve->energyUnits();

  // Fixed-geometry DRFs have no geometry to model, so no Geom & MC tab and no mode toggle.  A
  //  far-field DRF starts Geometry Modeled iff it already carries a Monte-Carlo response (or is a
  //  geometry-only import that has no efficiency yet, and so needs one).
  const bool fixed_geom = (m_orig && m_orig->isFixedGeometry());
  const bool has_geom_tab = !fixed_geom;
  m_geometryModeled = has_geom_tab
      && (!m_orig || m_orig->ceeloResponse() || !m_orig->isValid());

  // --- Tab: General (name / description) -----------------------------------
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    panel->addStyleClass( "DrfModifyPanel" );

    // Both label/field pairs go in a single grid, so the labels share one column and the fields
    //  line up on both edges - a row-per-pair cannot do that, since each row sizes independently.
    WContainerWidget *idGrid = panel->addNew<WContainerWidget>();
    idGrid->addStyleClass( "DrfModifyFieldGrid" );

    WLabel *nameLabel = idGrid->addNew<WLabel>( WString::tr("Name") );
    m_name = idGrid->addNew<WLineEdit>();
    m_name->setTextSize( 32 );
    m_name->setText( WString::fromUTF8( m_orig ? m_orig->name() : string() ) );
    m_name->changed().connect( this, &DrfModifyWidget::markEdited );
    nameLabel->setBuddy( m_name );

    WLabel *descLabel = idGrid->addNew<WLabel>( WString::tr("Description") );
    m_description = idGrid->addNew<WTextArea>();
    m_description->setColumns( 40 );
    m_description->setRows( 3 );
    m_description->setText( WString::fromUTF8( m_orig ? m_orig->description() : string() ) );
    m_description->changed().connect( this, &DrfModifyWidget::markEdited );
    descLabel->setBuddy( m_description );

    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-name"), std::move(panelOwned) );
    make_item_selectable( m_tabMenu, item );
  }

  // --- Tab: Geom & MC characterization (far-field only) --------------------
  if( has_geom_tab )
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    panel->addStyleClass( "DrfModifyGeomPanel" );

    // Flat Disk / Geometry Modeled toggle (checked == Geometry Modeled) above the tool.  In Flat
    //  Disk the tool is disabled wholesale (Wt propagates disabled to descendants).
    WContainerWidget *modeRow = panel->addNew<WContainerWidget>();
    modeRow->addStyleClass( "DrfModifyRow DrfModifyModeRow" );
    m_modeToggle = modeRow->addNew<SwitchCheckbox>( WString::tr("dmw-mode-flat"),
                                                    WString::tr("dmw-mode-geom") );
    m_modeToggle->setChecked( m_geometryModeled );
    HelpSystem::attachToolTipOn( modeRow, WString::tr("dmw-tt-mode"), true );
    m_modeToggle->checked().connect( this, &DrfModifyWidget::handleModeToggle );
    m_modeToggle->unChecked().connect( this, &DrfModifyWidget::handleModeToggle );

    auto toolOwned = make_unique<MakeMcResponseForDrf>( m_interspec, m_orig );
    m_mcTool = toolOwned.get();
    panel->addWidget( std::move(toolOwned) );

    // We drive generation from the footer "Generate Response" button (which first folds our edits
    //  into the seed), so hide the tool's own redundant generate button in Location Support.
    m_mcTool->setGenerateButtonHidden( true );

    // ContentLoading::Eager: this tool posts work to a worker thread and gets back to itself with
    //  `findById(...)`; a Lazy tab parks its contents in `WMenuItem::uContents_`, outside the
    //  widget tree, where `findById` can not see it - see the FWHM tab below.
    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-geometry"), std::move(panelOwned),
                                          ContentLoading::Eager );
    make_item_selectable( m_tabMenu, item );

    m_mcTool->setDisabled( !m_geometryModeled );  //Flat Disk greys the whole tool
  }//if( has_geom_tab )

  // --- Tab: FWHM -----------------------------------------------------------
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    WGridLayout *fwhmLayout = panel->setLayout( make_unique<WGridLayout>() );
    fwhmLayout->setContentsMargins( 0, 0, 0, 0 );

    // Fitting the FWHM needs a foreground to find peaks in, but viewing and hand-editing the one
    //  the detector already has does not - so the tool is always built, and only the "Fit FWHM"
    //  button (which the tool disables itself) needs a spectrum.
    const shared_ptr<const SpecUtils::Measurement> foreground
        = m_interspec ? m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground )
                      : nullptr;
    if( !foreground )
    {
      WText *note = fwhmLayout->addWidget( make_unique<WText>( WString::tr("dmw-fwhm-no-foreground") ), 0, 0 );
      note->addStyleClass( "DrfModifyNote" );
      note->setInline( false );
    }//if( no foreground )

    // MakeFwhmForDrf edits a (non-const) DRF; give it a private clone - we
    //  only read its coefficients back out on apply().
    shared_ptr<DetectorPeakResponse> fwhm_seed
        = m_orig ? make_shared<DetectorPeakResponse>( *m_orig ) : nullptr;

    // `ShowExisting`: open on whatever FWHM the DRF already has, so looking at this tab can never
    //  cost the user the one they had; the automated peak search (a multi-threaded fit of the whole
    //  spectrum) only runs when they ask for it with "Fit FWHM".
    // Always narrow: even at this dialog's widest, the tool shares its row with the equation
    //  controls, so the peak table never has room for the full-width column headers.
    auto toolOwned = make_unique<MakeFwhmForDrf>( MakeFwhmForDrf::InitialFit::ShowExisting,
                                                  m_interspec, fwhm_seed, true );
    m_fwhmTool = toolOwned.get();
    fwhmLayout->addWidget( std::move(toolOwned), 1, 0 );
    fwhmLayout->setRowStretch( 1, 1 );

    // ContentLoading::Eager, so the tool is in the widget tree (and hence resolvable by
    //  `findById`) while its peak search runs; with the default Lazy policy the widget sits in
    //  `WMenuItem::uContents_` with no widget parent, the search completion silently no-ops, and
    //  the tab shows "Currently searching for peaks..." forever.
    m_fwhmTabItem = m_tabMenu->addItem( WString::tr("dmw-tab-fwhm"), std::move(panelOwned),
                                        ContentLoading::Eager );
    make_item_selectable( m_tabMenu, m_fwhmTabItem );
  }

  // --- Tab: Anchor (efficiency representation + its uncertainty) ------------
  //  One of three editors, chosen by the efficiency form (see activeAnchorEditor):
  //   - the point table (energy / efficiency / stat % / corr % [/ cert % / source]) for a
  //     kEnergyEfficiencyPairs curve, and for Geometry Modeled, where the rows are the MC
  //     grounding anchors;
  //   - coefficient boxes plus their M*M covariance for a kExpOfLogPowerSeries curve;
  //   - the formula text for a kFunctialEfficienyForm curve, which also shows the point table with
  //     its Efficiency column hidden, so per-energy uncertainties can be given the same way.
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    panel->addStyleClass( "DrfModifyPanel DrfModifyAnchorPanel" );

    // Wording is set to match the visible editor in updateAnchorEditorVisibility().
    m_anchorHelp = panel->addNew<WText>();
    m_anchorHelp->setInline( false );

    const shared_ptr<const DetectorEfficiencyUncert> orig_uncert
        = m_orig ? m_orig->efficiencyUncert() : nullptr;

    // --- formula editor ---------------------------------------------------
    {
      auto formOwned = make_unique<WContainerWidget>();
      m_formulaEditor = formOwned.get();
      m_formulaEditor->addStyleClass( "DrfModifyFormula" );

      WContainerWidget *fcnRow = m_formulaEditor->addNew<WContainerWidget>();
      fcnRow->addStyleClass( "DrfModifyRow" );
      WLabel *fcnLabel = fcnRow->addNew<WLabel>( WString::tr("dmw-anchor-formula-label") );
      m_formulaText = fcnRow->addNew<WTextArea>();
      m_formulaText->setColumns( 40 );
      m_formulaText->setRows( 3 );
      m_formulaText->setPlaceholderText( WString::tr("dmw-anchor-formula-empty") );
      if( eff_curve && eff_curve->isValid()
          && (eff_curve->form() == DetectorPeakResponse::kFunctialEfficienyForm) )
        m_formulaText->setText( WString::fromUTF8( eff_curve->formula() ) );
      m_formulaText->changed().connect( this, [this](){ validateFormula(); markEdited(); } );
      fcnLabel->setBuddy( m_formulaText );
      HelpSystem::attachToolTipOn( fcnRow, WString::tr("dmw-tt-anchor-formula"), true );

      panel->addWidget( std::move(formOwned) );
    }//formula editor

    // --- coefficient editor -----------------------------------------------
    {
      auto coefOwned = make_unique<WContainerWidget>();
      m_coefEditor = coefOwned.get();
      m_coefEditor->addStyleClass( "DrfModifyCoefs" );

      WGroupBox *coefBox = m_coefEditor->addNew<WGroupBox>( WString::tr("dmw-anchor-coefs") );
      coefBox->addStyleClass( "DrfModifyCoefBox" );
      m_coefParams = coefBox->addNew<WContainerWidget>();
      m_coefParams->addStyleClass( "DrfModifyCoefParams" );

      WContainerWidget *coefBtns = coefBox->addNew<WContainerWidget>();
      m_addCoef = coefBtns->addNew<WPushButton>( WString::tr("dmw-anchor-add-term") );
      m_addCoef->addStyleClass( "LinkBtn" );
      m_addCoef->clicked().connect( this, [this](){ addCoefficient(); markEdited(); } );
      m_removeCoef = coefBtns->addNew<WPushButton>( WString::tr("dmw-anchor-remove-term") );
      m_removeCoef->addStyleClass( "LinkBtn" );
      m_removeCoef->clicked().connect( this, [this](){ removeCoefficient(); markEdited(); } );

      // The coefficient covariance; the same sigma/rho machinery the energy-node matrix used, which
      //  at the 3-6 terms an equation has is finally a size it is usable at.
      m_covTable = m_coefEditor->addNew<WTable>();
      m_covTable->addStyleClass( "DrfModifyBandTable DrfModifyCovTable" );

      if( eff_curve && eff_curve->isValid()
          && (eff_curve->form() == DetectorPeakResponse::kExpOfLogPowerSeries) )
      {
        for( const float coef : eff_curve->expOfLogPowerSeriesCoeffs() )
          addCoefficient( coef );
      }
      seedCoefCovFromUncert( orig_uncert );  //fills the shadow and does the initial rebuildCovTable()

      panel->addWidget( std::move(coefOwned) );
    }//coefficient editor

    // --- point-table editor -----------------------------------------------
    {
      auto ptsOwned = make_unique<WContainerWidget>();
      m_pointsEditor = ptsOwned.get();
      m_pointsEditor->addStyleClass( "DrfModifyAnchor" );

      // Reference distance the absolute curve is anchored at (cm); an intrinsic curve is
      //  per-gamma-striking-the-face, at no distance, so the row is hidden for it.
      double refDistCm = 0.0;
      if( have_measured )
      {
        for( const MeasuredEffPoint &p : measured->points() )
        {
          if( p.distance > 0.0f ){ refDistCm = p.distance / PhysicalUnits::cm; break; }
        }
      }//if( have_measured )
      if( (refDistCm <= 0.0) && m_orig )
        refDistCm = m_orig->absoluteEfficiencyDistance() / PhysicalUnits::cm;

      WContainerWidget *distRow = m_pointsEditor->addNew<WContainerWidget>();
      distRow->addStyleClass( "DrfModifyRow" );
      distRow->addNew<WLabel>( WString::tr("dmw-anchor-ref-dist") );
      m_anchorRefDistance = distRow->addNew<WLineEdit>();
      m_anchorRefDistance->setTextSize( 8 );
      m_anchorRefDistance->changed().connect( this, &DrfModifyWidget::markEdited );
      if( refDistCm > 0.0 )
      {
        char buf[32];
        snprintf( buf, sizeof(buf), "%.4g", refDistCm );
        m_anchorRefDistance->setText( buf );
      }
      distRow->addNew<WLabel>( WString::tr("dmw-anchor-cm") );
      distRow->setHidden( !m_anchorIsAbsolute );

      WContainerWidget *defRow = m_pointsEditor->addNew<WContainerWidget>();
      defRow->addStyleClass( "DrfModifyRow" );
      defRow->addNew<WLabel>( WString::tr("dmw-anchor-default-uncert") );
      m_anchorDefaultUncert = defRow->addNew<WLineEdit>();
      m_anchorDefaultUncert->setTextSize( 6 );
      m_anchorDefaultUncert->changed().connect( this, &DrfModifyWidget::markEdited );
      // Only a fallback for blank per-point stat cells - except for a formula curve with no rows at
      //  all, where it is the whole (flat) uncertainty.  Left blank when the source stated no
      //  uncertainties at all (a GADRAS Efficiency.csv), rather than inventing one.
      if( have_measured )
        m_anchorDefaultUncert->setText( "5" );
      defRow->addNew<WLabel>( WString::tr("dmw-anchor-percent") );

      m_anchorTableWrap = m_pointsEditor->addNew<WContainerWidget>();
      m_anchorTableWrap->addStyleClass( "DrfModifyAnchorTableWrap" );
      m_anchorTable = m_anchorTableWrap->addNew<WTable>();
      m_anchorTable->addStyleClass( "DrfModifyAnchorTable" );
      {
        int col = 0;
        m_anchorTable->elementAt(0,col++)->addNew<WText>( WString::tr("dmw-anchor-energy") );
        WTableCell * const effHeader = m_anchorTable->elementAt(0,col++);
        effHeader->addNew<WText>( WString::tr("dmw-anchor-eff") );
        effHeader->addStyleClass( "DrfEffCol" );
        m_anchorTable->elementAt(0,col++)->addNew<WText>( WString::tr("dmw-anchor-stat") );
        // The correlated component, shown for every curve.  With a source column it is the
        //  certificate uncertainty (100% correlated within a source); without one it is correlated
        //  across energy by the model below.
        m_anchorTable->elementAt(0,col++)->addNew<WText>(
                    WString::tr( m_anchorHasSourceCol ? "dmw-anchor-cert" : "dmw-anchor-corr" ) );
        if( m_anchorHasSourceCol )
          m_anchorTable->elementAt(0,col++)->addNew<WText>( WString::tr("dmw-anchor-source") );
      }

      WContainerWidget *btns = m_pointsEditor->addNew<WContainerWidget>();
      m_addAnchor = btns->addNew<WPushButton>( WString::tr("dmw-anchor-add") );
      m_addAnchor->addStyleClass( "LinkBtn" );
      m_addAnchor->clicked().connect( this, [this](){ addAnchorRow( 0.0f, 0.0f, 0.0f, 0.0f, string() ); markEdited(); } );
      m_removeAnchor = btns->addNew<WPushButton>( WString::tr("dmw-anchor-remove") );
      m_removeAnchor->addStyleClass( "LinkBtn" );
      m_removeAnchor->clicked().connect( this, [this](){ removeAnchorRow(); markEdited(); } );

      if( have_measured )
      {
        for( const MeasuredEffPoint &p : measured->points() )
          addAnchorRow( p.energy, p.efficiency, p.fracStatUncert, p.fracCertUncert, p.sourceKey );
      }else if( have_pairs )
      {
        // Seed the uncertainty columns from whatever the DRF actually carries: the split the
        //  covariance was built from when it has one, else the 1-sigma envelope as a purely
        //  correlated component (which reproduces the same covariance in fully-correlated mode).
        for( const DetectorPeakResponse::EnergyEfficiencyPair &p : eff_curve->energyEfficiencies() )
        {
          float stat = 0.0f, corr = 0.0f;
          uncert_for_energy( orig_uncert, p.energy * m_anchorEnergyUnits, stat, corr );
          addAnchorRow( p.energy, p.efficiency, stat, corr, string() );
        }
      }else if( orig_uncert && orig_uncert->hasNodeCovariance() )
      {
        // A formula curve carries no points of its own, so its rows come from the covariance nodes.
        for( const float energy_kev : orig_uncert->covarianceEnergies() )
        {
          float stat = 0.0f, corr = 0.0f;
          uncert_for_energy( orig_uncert, energy_kev, stat, corr );
          addAnchorRow( energy_kev / m_anchorEnergyUnits, 0.0f, stat, corr, string() );
        }
      }//if( have_measured ) / else if( have_pairs ) / else if( have a covariance )
      m_removeAnchor->setEnabled( !m_anchors.empty() );

      // How the correlated column is correlated across energy.  Shared with the .ecc import
      //  dialogs, so the mode -> correlation-length mapping and the example table live in one
      //  place, and the Modify dialog speaks the same language the import did.
      if( !m_anchorHasSourceCol )
      {
        vector<float> node_energies, corr_frac, uncorr_frac;
        for( const AnchorRow &r : m_anchors )
        {
          float energy = 0.0f, corr = 0.0f, stat = 0.0f;
          try{ energy = std::stof( r.energy->text().toUTF8() ); }catch( std::exception & ){ continue; }
          if( energy <= 0.0f )
            continue;
          try{ corr = 0.01f * std::stof( r.cert->text().toUTF8() ); }catch( std::exception & ){ corr = 0.0f; }
          try{ stat = 0.01f * std::stof( r.stat->text().toUTF8() ); }catch( std::exception & ){ stat = 0.0f; }
          node_energies.push_back( energy * m_anchorEnergyUnits );
          corr_frac.push_back( std::max( 0.0f, corr ) );
          uncorr_frac.push_back( std::max( 0.0f, stat ) );
        }//for( const AnchorRow &r : m_anchors )

        m_uncertOptions = m_pointsEditor->addNew<EccUncertOptions>( node_energies, corr_frac,
                                                                    uncorr_frac );
        m_uncertOptions->setImportToggleVisible( false );
        // correlationLength() <= 0 is ambiguous between "uncorrelated was chosen" and "the matrix
        //  was set directly", so fall back to fully correlated - the anchorsMatchSeed() guard is
        //  what keeps that guess from rewriting a covariance the user never touched.
        const double seed_len = (orig_uncert && (orig_uncert->correlationLength() > 0.0))
                       ? orig_uncert->correlationLength()
                       : DetectorEfficiencyUncert::sm_fullyCorrelatedLength;
        m_uncertOptions->setCorrelationLength( seed_len );
        m_uncertOptions->changed().connect( this, &DrfModifyWidget::markEdited );
      }//if( !m_anchorHasSourceCol )

      panel->addWidget( std::move(ptsOwned) );
    }//point-table editor

    // keV / MeV the coefficient equation or formula is written in.
    {
      m_effUnitsRow = panel->addNew<WContainerWidget>();
      m_effUnitsRow->addStyleClass( "DrfModifyRow DrfModifyEffUnitsRow" );
      WLabel *unitsLabel = m_effUnitsRow->addNew<WLabel>( WString::tr("dmw-anchor-energy-units") );
      m_effEnergyUnits = m_effUnitsRow->addNew<WComboBox>();
      m_effEnergyUnits->addItem( WString::tr("dmw-anchor-kev") );
      m_effEnergyUnits->addItem( WString::tr("dmw-anchor-mev") );
      const float curve_units = (eff_curve && eff_curve->isValid() && (eff_curve->energyUnits() > 0.0f))
                              ? eff_curve->energyUnits() : static_cast<float>(PhysicalUnits::keV);
      m_effEnergyUnits->setCurrentIndex( (curve_units > 10.0f) ? 1 : 0 );
      m_effEnergyUnits->activated().connect( this, &DrfModifyWidget::markEdited );
      unitsLabel->setBuddy( m_effEnergyUnits );
    }

    updateAnchorEditorVisibility();

    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-anchor"), std::move(panelOwned) );
    make_item_selectable( m_tabMenu, item );
  }

  m_tabMenu->itemSelected().connect( this, &DrfModifyWidget::handleTabSelected );

  // Both child tools hand their state changes here, so one undo step covers the whole dialog; an
  //  edit also marks any generated response stale (see markEdited).  m_mcTool is null for a
  //  fixed-geometry DRF (no Geom & MC tab).
  if( m_mcTool )
  {
    m_mcTool->userChanged().connect( this, &DrfModifyWidget::markEdited );
    m_mcTool->responseGenerated().connect( this, &DrfModifyWidget::handleResponseGenerated );
  }//if( m_mcTool )
  if( m_fwhmTool )
  {
    m_fwhmTool->setOwnerHandlesUndoRedo( true );
    m_fwhmTool->stateChanged().connect( this, &DrfModifyWidget::markEdited );
  }//if( m_fwhmTool )

  m_tabMenu->select( 0 );

  // Footer: the export note on the left, and (Geometry-Modeled only) a "Generate Response" button on
  //  the right that regenerates the MC response from the live edits.
  auto footerOwned = make_unique<WContainerWidget>();
  WContainerWidget *footerRow = footerOwned.get();
  footerRow->addStyleClass( "DrfModifyFooterRow" );

  WText *note = footerRow->addNew<WText>( WString::tr("dmw-export-note") );
  note->addStyleClass( "DrfModifyNote" );
  note->setInline( false );

  if( m_mcTool )
  {
    m_generateBtn = footerRow->addNew<WPushButton>( WString::tr("dmw-generate-btn") );
    m_generateBtn->addStyleClass( "DrfModifyGenerateBtn" );
    m_generateBtn->clicked().connect( this, &DrfModifyWidget::handleGenerateResponse );
  }//if( m_mcTool )

  if( narrow_layout )
    layout->addWidget( std::move(footerOwned), 2, 0 );        //below menu-row and stack-row
  else
    layout->addWidget( std::move(footerOwned), 1, 0, 1, 2 );  //below, spanning menu and stack columns

  updateGenerateButton();
  // What the Anchor tab opened with, so anchorsMatchSeed() can tell "the user changed nothing" from
  //  "the user re-typed the same numbers" - the apply paths leave the DRF alone in the first case.
  m_seedState = currentState();

}//DrfModifyWidget constructor


DrfModifyWidget::~DrfModifyWidget()
{
}


void DrfModifyWidget::handleTabSelected( Wt::WMenuItem *item )
{
  scheduleUndoRedoStep();

  // `startAutomatedPeakSearch` is idempotent, so re-visiting the tab costs nothing.
  if( item && (item == m_fwhmTabItem) && m_fwhmTool )
    m_fwhmTool->startAutomatedPeakSearch();
}//void handleTabSelected( Wt::WMenuItem *item )


Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &DrfModifyWidget::updatedDrf()
{
  return m_updatedDrf;
}


Wt::Signal<bool> &DrfModifyWidget::mcResponseAvailable()
{
  // Only reached when #needsMcResponse, which is false without an MC tool.
  assert( m_mcTool );
  return m_mcTool->validationChanged();
}


bool DrfModifyWidget::needsMcResponse() const
{
  // A fixed-geometry DRF has no MC tool, so a Monte-Carlo response can neither be generated nor is
  //  it what such a DRF would need - never hold "Use" hostage to one.
  return (m_mcTool && (!m_orig || !m_orig->isValid()));
}


std::shared_ptr<const DetectorPeakResponse> DrfModifyWidget::originalDrf() const
{
  return m_orig;
}


void DrfModifyWidget::addAnchorRow( const float energy, const float efficiency,
                                    const float fracStatUncert, const float fracCertUncert,
                                    const std::string &sourceKey )
{
  const int row = m_anchorTable->rowCount();  //row 0 is the header

  auto make_edit = [this,row]( const int col, const float value, const bool asPercent ) -> WLineEdit * {
    WTableCell * const cell = m_anchorTable->elementAt(row,col);
    // The Efficiency column is hidden by CSS for a formula curve; WTableColumn::setStyleClass does
    //  not reach the cells, so each one carries the class itself.
    if( col == 1 )
      cell->addStyleClass( "DrfEffCol" );
    WLineEdit *edit = cell->addNew<WLineEdit>();
    edit->setTextSize( 9 );
    const float shown = asPercent ? 100.0f*value : value;
    if( shown != 0.0f )
    {
      char buf[32];
      snprintf( buf, sizeof(buf), "%.6g", shown );
      edit->setText( buf );
    }
    edit->changed().connect( this, &DrfModifyWidget::markEdited );
    return edit;
  };

  // Uncertainties are shown/edited as percentages; a blank stat cell inherits the default % on
  //  apply.  Only the source column is conditional (see m_anchorHasSourceCol); the efficiency one
  //  is built always and hidden by CSS for a formula curve.
  int col = 0;
  AnchorRow r;
  r.energy = make_edit( col++, energy, false );
  r.eff    = make_edit( col++, efficiency, false );
  r.stat   = make_edit( col++, fracStatUncert, true );
  r.cert   = make_edit( col++, fracCertUncert, true );
  r.source = nullptr;
  if( m_anchorHasSourceCol )
  {
    r.source = m_anchorTable->elementAt(row,col++)->addNew<WLineEdit>();
    r.source->setTextSize( 10 );
    r.source->setText( WString::fromUTF8(sourceKey) );
    r.source->changed().connect( this, &DrfModifyWidget::markEdited );
  }//if( m_anchorHasSourceCol )
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


shared_ptr<DetectorEfficiencyUncert> DrfModifyWidget::buildUncertFromRows()
{
  if( !m_uncertOptions )
    return nullptr;

  // Default for blank per-point stat cells, and - with no rows at all - the whole uncertainty.
  float defaultFrac = 0.0f;
  {
    const string ds = m_anchorDefaultUncert->text().toUTF8();
    try{ defaultFrac = 0.01f * std::stof( ds ); }catch( std::exception & ){ defaultFrac = 0.0f; }
    if( defaultFrac < 0.0f )
      defaultFrac = 0.0f;
  }

  vector<float> energies, corr_fracs, stat_fracs;
  for( const AnchorRow &r : m_anchors )
  {
    const string es = r.energy->text().toUTF8();
    if( es.empty() )
      continue;

    float energy = 0.0f;
    try{ energy = std::stof( es ); }catch( std::exception & ){ continue; }
    if( energy <= 0.0f )
      continue;

    const string ss = r.stat->text().toUTF8();
    float statFrac = defaultFrac;
    if( !ss.empty() )
    {
      try{ statFrac = 0.01f * std::stof( ss ); }catch( std::exception & ){ statFrac = defaultFrac; }
    }

    float corrFrac = 0.0f;
    const string cs = r.cert->text().toUTF8();
    if( !cs.empty() )
    {
      try{ corrFrac = 0.01f * std::stof( cs ); }catch( std::exception & ){ corrFrac = 0.0f; }
    }

    // Covariance node energies are keV regardless of the units the curve (and hence the column)
    //  uses - see DetectorEfficiencyUncert.
    energies.push_back( energy * m_anchorEnergyUnits );
    corr_fracs.push_back( std::max( 0.0f, corrFrac ) );
    stat_fracs.push_back( std::max( 0.0f, statFrac ) );
  }//for( const AnchorRow &r : m_anchors )

  // No rows, but a default: one flat, fully-correlated node.  A single-node covariance extrapolates
  //  constantly, so this is the same fractional uncertainty at every energy - which is exactly what
  //  a flat "Default uncertainty" means.
  if( energies.empty() )
  {
    if( defaultFrac <= 0.0f )
      return nullptr;

    const float ref_energy = ((m_orig && (m_orig->lowerEnergy() > 0.0f)) ? m_orig->lowerEnergy()
                                                                         : 661.7f);
    try
    {
      return DetectorEfficiencyUncert::fromPointUncerts( { ref_energy }, { defaultFrac },
                             DetectorEfficiencyUncert::sm_fullyCorrelatedLength );
    }catch( std::exception & )
    {
      return nullptr;
    }
  }//if( energies.empty() )

  // Keep the widget's example table in step with the rows, then build directly: buildUncert()
  //  declines below two nodes, and one node is a legitimate (flat, fully-correlated) uncertainty.
  //  effectiveCorrLength() is still the single source of the mode -> length mapping.
  m_uncertOptions->setPoints( energies, corr_fracs, stat_fracs );

  try
  {
    return DetectorEfficiencyUncert::fromCorrelatedPlusDiagonal( energies, corr_fracs, stat_fracs,
                                              m_uncertOptions->effectiveCorrLength() );
  }catch( std::exception & )
  {
    return nullptr;
  }
}//buildUncertFromRows()


void DrfModifyWidget::applyAnchorEdits( DetectorPeakResponse &working )
{
  if( !m_anchorTable )
    return;  //the point editor was not built

  // Nothing was touched: leave the curve and uncertainty bit-identical.  Matters for a covariance
  //  the correlated+diagonal model cannot reproduce (a source-blocked MeasuredDrfPoints matrix, or
  //  a hand-authored one restored from a URL with no component split).
  if( anchorsMatchSeed() )
    return;

  // Default fractional statistical uncertainty for blank per-point stat cells.
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
  if( m_anchorIsAbsolute && (refDistCm <= 0.0) )
    return;  //can't build absolute efficiency without a valid distance

  const double refDist = refDistCm * PhysicalUnits::cm;

  vector<DetectorPeakResponse::EnergyEffPoint> effpts;
  vector<MeasuredEffPoint> measpts;
  for( const AnchorRow &r : m_anchors )
  {
    const string es = r.energy->text().toUTF8();
    const string fs = r.eff->text().toUTF8();
    if( es.empty() && fs.empty() )
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

    const string ss = r.stat->text().toUTF8();
    float statFrac = defaultFrac;
    if( !ss.empty() )
    {
      try{ statFrac = 0.01f * std::stof( ss ); }catch( std::exception & ){ statFrac = defaultFrac; }
    }
    if( statFrac < 0.0f )
      statFrac = 0.0f;

    // The correlated component is parsed for every curve; only its meaning differs (see
    //  m_anchorHasSourceCol).
    float certFrac = 0.0f;
    string srcKey;
    const string cs = r.cert->text().toUTF8();
    if( !cs.empty() )
    {
      try{ certFrac = 0.01f * std::stof( cs ); }catch( std::exception & ){ certFrac = 0.0f; }
    }
    if( certFrac < 0.0f )
      certFrac = 0.0f;
    if( r.source )
      srcKey = r.source->text().toUTF8();

    DetectorPeakResponse::EnergyEffPoint e;
    e.energy = energy;
    e.efficiency = eff;
    // The far-field curve carries the combined 1-sigma as a fallback per-point uncert; the richer
    //  per-source covariance is set from the measured points below (and overwrites it).
    const float combo = std::sqrt( statFrac*statFrac + certFrac*certFrac );
    if( combo > 0.0f )
      e.efficiencyUncert = eff * combo;
    effpts.push_back( e );

    MeasuredEffPoint m;
    m.energy = energy;
    m.efficiency = eff;
    m.fracStatUncert = statFrac;
    m.fracCertUncert = certFrac;
    m.sourceKey = srcKey;
    m.distance = m_anchorIsAbsolute ? static_cast<float>( refDist ) : -1.0f;
    measpts.push_back( m );
  }//for( const AnchorRow &r : m_anchors )

  if( effpts.size() < 2 )
    return;  //need at least two points for a curve; leave the DRF untouched

  try
  {
    // Points whose correlated column is blocked per source keep the far-field characterization
    //  path, which is what records `MeasuredDrfPoints` and needs the reference distance.
    if( m_anchorHasSourceCol )
    {
      const float diameter = working.detectorDiameter();
      if( diameter <= 0.0f )
        return;  //a far-field curve needs a positive diameter; leave the DRF as-is

      working.setEfficiencyPoints( effpts, diameter, refDist,
                                   DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

      auto meas = make_shared<MeasuredDrfPoints>();
      meas->setPoints( measpts );
      working.setMeasuredPoints( meas );

      // Drive the efficiency uncertainty from the edited stat/cert/source structure (C[i][j] =
      //  delta_ij*stat_i^2 + cert_i*cert_j*[same source]); a CeeLo response, when attached,
      //  overrides this at query time.
      const shared_ptr<DetectorEfficiencyUncert> uncert = meas->toEfficiencyUncert();
      if( uncert && !uncert->isEmpty() )
        working.setEfficiencyUncert( uncert );

      return;
    }//if( m_anchorHasSourceCol )

    // Everything else - an .ecc, a GADRAS Efficiency.csv, an ANGLE .outx - only changes its
    //  numbers.  setEfficiencyPoints would rewrite the geometry type, diameter and flags, which is
    //  what a fixed-geometry DRF has no spare copy of, and it rejects the zero diameter an .ecc has.
    vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
    pairs.reserve( effpts.size() );
    for( const DetectorPeakResponse::EnergyEffPoint &e : effpts )
    {
      DetectorPeakResponse::EnergyEfficiencyPair p;
      p.energy = e.energy;
      p.efficiency = e.efficiency;
      pairs.push_back( p );
    }

    auto curve = make_shared<DetectorEfficiencyCurve>();
    curve->setFromPairs( pairs, m_anchorEnergyUnits );
    curve->setUncertainty( buildUncertFromRows() );
    working.replaceEfficiencyCurve( curve );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("dmw-err-points-invalid").arg(e.what()), WarningWidget::WarningMsgHigh );
  }
}//applyAnchorEdits(...)


bool DrfModifyWidget::ToolState::operator==( const ToolState &rhs ) const
{
  if( (name != rhs.name) || (description != rhs.description) || (tabIndex != rhs.tabIndex)
     || (geometryModeled != rhs.geometryModeled)
     || (coefCovMatrix != rhs.coefCovMatrix)
     || (coefficients != rhs.coefficients)
     || (formula != rhs.formula)
     || (anchors != rhs.anchors)
     || (anchorRefDistance != rhs.anchorRefDistance)
     || (anchorDefaultUncert != rhs.anchorDefaultUncert)
     || (anchorCorrLength != rhs.anchorCorrLength)
     || (efficiencyEnergyUnits != rhs.efficiencyEnergyUnits)
     || (mc != rhs.mc) )
  {
    return false;
  }

  if( !fwhm || !rhs.fwhm )
    return (!fwhm == !rhs.fwhm);

  return ((*fwhm) == (*rhs.fwhm));
}//ToolState::operator==


std::shared_ptr<DrfModifyWidget::ToolState> DrfModifyWidget::currentState() const
{
  auto state = make_shared<ToolState>();

  state->name = m_name->text().toUTF8();
  state->description = m_description->text().toUTF8();
  state->tabIndex = m_tabMenu->currentIndex();
  state->geometryModeled = m_geometryModeled;

  // The covariance shadow is the authoritative numeric state (not the widgets).
  state->coefCovMatrix = m_coefCovMatrix;

  for( const NativeFloatSpinBox * const edit : m_coefEdits )
    state->coefficients.push_back( edit->text().toUTF8() );

  if( m_formulaText )
    state->formula = m_formulaText->text().toUTF8();

  for( const AnchorRow &r : m_anchors )
  {
    state->anchors.push_back( { r.energy->text().toUTF8(), r.eff->text().toUTF8(),
                                r.stat->text().toUTF8(),
                                r.cert ? r.cert->text().toUTF8() : string(),
                                r.source ? r.source->text().toUTF8() : string() } );
  }//for( const AnchorRow &r : m_anchors )

  if( m_anchorRefDistance )
    state->anchorRefDistance = m_anchorRefDistance->text().toUTF8();
  if( m_anchorDefaultUncert )
    state->anchorDefaultUncert = m_anchorDefaultUncert->text().toUTF8();
  if( m_uncertOptions )
    state->anchorCorrLength = m_uncertOptions->effectiveCorrLength();
  if( m_effEnergyUnits )
    state->efficiencyEnergyUnits = equationEnergyUnits();

  if( m_mcTool )
    state->mc = m_mcTool->currentState();
  if( m_fwhmTool )
    state->fwhm = m_fwhmTool->currentState();

  return state;
}//DrfModifyWidget::currentState()


void DrfModifyWidget::setState( const std::shared_ptr<const ToolState> &state )
{
  if( !state )
    return;

  m_restoringState = true;
  m_currentState = state;

  m_name->setText( WString::fromUTF8(state->name) );
  m_description->setText( WString::fromUTF8(state->description) );

  if( m_modeToggle )
  {
    m_geometryModeled = state->geometryModeled;
    m_modeToggle->setChecked( m_geometryModeled );
    if( m_mcTool )
      m_mcTool->setDisabled( !m_geometryModeled );
  }//if( m_modeToggle )

  // Coefficient boxes, then the covariance shadow and the table it drives.  Bound both loops by the
  //  widget count: removeCoefficient() keeps at least one term, so a target of zero would otherwise
  //  never be reached.
  while( (m_coefEdits.size() > state->coefficients.size()) && (m_coefEdits.size() > 1) )
    removeCoefficient();
  while( m_coefEdits.size() < state->coefficients.size() )
    addCoefficient();
  const size_t ncoef_set = std::min( m_coefEdits.size(), state->coefficients.size() );
  for( size_t i = 0; i < ncoef_set; ++i )
    m_coefEdits[i]->setText( WString::fromUTF8(state->coefficients[i]) );

  m_coefCovMatrix = state->coefCovMatrix;
  rebuildCovTable();

  if( m_formulaText )
    m_formulaText->setText( WString::fromUTF8(state->formula) );
  if( m_effEnergyUnits )
    m_effEnergyUnits->setCurrentIndex( (state->efficiencyEnergyUnits > 10.0f) ? 1 : 0 );

  if( m_anchorTable )
  {
    while( !m_anchors.empty() )
      removeAnchorRow();
    for( const std::array<string,5> &a : state->anchors )
    {
      addAnchorRow( 0.0f, 0.0f, 0.0f, 0.0f, string() );
      m_anchors.back().energy->setText( WString::fromUTF8(a[0]) );
      m_anchors.back().eff->setText( WString::fromUTF8(a[1]) );
      m_anchors.back().stat->setText( WString::fromUTF8(a[2]) );
      if( m_anchors.back().cert )
        m_anchors.back().cert->setText( WString::fromUTF8(a[3]) );
      if( m_anchors.back().source )
        m_anchors.back().source->setText( WString::fromUTF8(a[4]) );
    }//for( const std::array<string,5> &a : state->anchors )

    m_anchorRefDistance->setText( WString::fromUTF8(state->anchorRefDistance) );
    m_anchorDefaultUncert->setText( WString::fromUTF8(state->anchorDefaultUncert) );
  }//if( m_anchorTable )

  // setCorrelationLength does not emit changed(), so this stays out of the edit path.
  if( m_uncertOptions )
    m_uncertOptions->setCorrelationLength( state->anchorCorrLength );

  updateAnchorEditorVisibility();

  if( m_mcTool )
    m_mcTool->setState( state->mc );
  if( m_fwhmTool && state->fwhm )
    m_fwhmTool->setState( state->fwhm );

  if( (state->tabIndex >= 0) && (state->tabIndex < m_tabMenu->count()) )
    m_tabMenu->select( state->tabIndex );

  // Restoring is not a user edit, so nothing here should become the next undo step; a restored
  //  snapshot is not treated as newly edited (staleness is re-derived by the next real edit).
  m_renderFlags.clear( RenderActions::AddUndoRedoStep );
  m_changedSinceGenerate = false;
  updateGenerateButton();
  m_restoringState = false;
}//DrfModifyWidget::setState(...)


void DrfModifyWidget::scheduleUndoRedoStep()
{
  if( m_restoringState )
    return;

  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void scheduleUndoRedoStep()


void DrfModifyWidget::markEdited()
{
  if( m_restoringState )
    return;

  // The userChanged the MC tool emits right after a successful generation is not a user edit, so it
  //  must not re-flag the fresh response stale - consume the one-shot suppression instead.
  if( m_suppressNextEditMark )
    m_suppressNextEditMark = false;
  else
    m_changedSinceGenerate = true;

  updateGenerateButton();
  scheduleUndoRedoStep();
}//void markEdited()


void DrfModifyWidget::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool add_step = m_renderFlags.test( RenderActions::AddUndoRedoStep );
  const bool rebuild_cov = m_renderFlags.test( RenderActions::RebuildCovTable );
  m_renderFlags = Wt::WFlags<RenderActions>();

  // The covariance cell edits defer their table rebuild here, so the WLineEdit whose `changed()`
  //  fired is not deleted from inside its own event handler.
  if( rebuild_cov )
    rebuildCovTable();

  // Re-baseline on every render, not just flagged ones, so a change made without recording a step
  //  (a restore, a generation landing) doesnt leave a stale baseline for the next edit to diff.
  doAddUndoRedoStep( add_step );

  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void DrfModifyWidget::doAddUndoRedoStep( const bool add_step )
{
  const shared_ptr<const ToolState> prev = m_currentState;
  const shared_ptr<const ToolState> current = currentState();
  m_currentState = current;

  if( !add_step || !prev || !current || ((*prev) == (*current)) )
    return;

  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( !undoRedo || !undoRedo->canAddUndoRedoNow() )
    return;

  // Resolve the dialog when the step runs, rather than capturing it: this widget may well have been
  //  closed and re-created since (see the comments on `UndoRedoManager::addUndoRedoStep`).  A step
  //  whose dialog is gone simply does nothing.
  auto apply = []( const shared_ptr<const ToolState> &state ){
    InterSpec *viewer = InterSpec::instance();
    DrfModifyWidget *tool = viewer ? viewer->drfModifyWidget() : nullptr;
    if( tool )
      tool->setState( state );
  };

  undoRedo->addUndoRedoStep( [prev,apply](){ apply(prev); },
                             [current,apply](){ apply(current); },
                             "Modify-DRF tool change" );
}//void doAddUndoRedoStep( const bool add_step )


std::shared_ptr<DetectorPeakResponse> DrfModifyWidget::buildWorkingDrf( const bool includeMcResponse )
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

  // Anchor tab: apply whichever editor is showing.  Done before the MC step so a grounding sees
  //  the edited points.
  switch( activeAnchorEditor() )
  {
    case AnchorEditor::Points:       applyAnchorEdits( *working );      break;
    case AnchorEditor::Coefficients: applyCoefficientEdits( *working ); break;
    case AnchorEditor::Formula:      applyFormulaEdits( *working );     break;
  }//switch( activeAnchorEditor() )

  // FWHM: take whatever the FWHM tab is showing, without triggering the tool's own
  //  detector-changed emit.  The tab opens seeded from this DRF, so this is a no-op unless the user
  //  actually changed something there.
  if( m_fwhmTool )
  {
    const shared_ptr<MakeFwhmForDrf::ToolState> fwhm = m_fwhmTool->currentState();
    const DetectorPeakResponse::ResolutionFnctForm form
        = DetectorPeakResponse::ResolutionFnctForm( std::max( 0, fwhm->m_fwhm_index ) );

    // `setFwhmCoefficients` throws if the coefficients dont match the forms arity - dont let that
    //  escape into a signal handler and take the rest of the users edits with it.
    try
    {
      if( form == DetectorPeakResponse::kNumResolutionFnctForm )
      {
        working->setFwhmCoefficients( {}, form );   //the user chose "None": clear the FWHM
      }else if( !fwhm->m_parameters.empty() )
      {
        working->setFwhmCoefficients( fwhm->m_parameters, form );
      }else
      {
        // An equation form is selected, but it has no coefficients - the user changed the form and
        //  has not fit (or filled in) one yet.  Say so, rather than closing as if the equation
        //  showing had been applied; the DRF keeps whatever FWHM it came in with.
        passMessage( WString::tr("dmw-err-no-fwhm-fit"), WarningWidget::WarningMsgHigh );
      }
    }catch( std::exception &e )
    {
      passMessage( WString::tr("dmw-err-fwhm-not-applied").arg(e.what()),
                   WarningWidget::WarningMsgHigh );
    }
  }//if( m_fwhmTool )

  // Monte-Carlo response.  In Geometry-Modeled mode, attach the generated response (keeping any the
  //  DRF already carried when nothing new was generated); in Flat Disk, or as a regeneration seed,
  //  detach it - CRITICAL for a regeneration seed, since an attached CeeLo response overrides the
  //  manual points/covariance at query time (the very thing the regeneration is grounding on).
  const shared_ptr<const ceelo::DetectorResponse> resp
      = (includeMcResponse && m_geometryModeled && m_mcTool) ? m_mcTool->generatedResponse() : nullptr;
  if( includeMcResponse && m_geometryModeled && resp )
  {
    // A DRF built from geometry alone has no efficiency curve of its own, so sample the response
    //  into one - the Monte-Carlo "backbone" points.  Every EffEval query would already dispatch to
    //  the response, but without a curve the DRF is not valid, cannot be serialized, and cannot be
    //  exported; a Detector.dat or .detx import with no measured efficiency lands here.  Done BEFORE
    //  setCeeloResponse: setEfficiencyPoints recomputes the hash and resets flags, and the response
    //  should be attached to the finished curve.
    if( !working->isValid() )
    {
      try
      {
        CeeLoUtils::setLegacyEfficiencyFromResponse( *working, resp );
      }catch( std::exception &e )
      {
        passMessage( WString::tr("dmw-err-no-backbone").arg(e.what()),
                     WarningWidget::WarningMsgHigh );
      }
    }//if( !working->isValid() )

    working->setCeeloResponse( resp );
  }else if( includeMcResponse && !m_geometryModeled )
  {
    working->setCeeloResponse( nullptr );  //Flat Disk: detach any geometry-modeled response
  }else if( !includeMcResponse )
  {
    working->setCeeloResponse( nullptr );  //regeneration seed: manual points/covariance must drive
  }
  //else: Geometry Modeled with nothing newly generated - keep whatever response the DRF came with.

  return working;
}//buildWorkingDrf(...)


void DrfModifyWidget::apply()
{
  shared_ptr<DetectorPeakResponse> working = buildWorkingDrf( true );
  m_updatedDrf.emit( working );
}//apply()


void DrfModifyWidget::requestApply()
{
  // Flat Disk, but the original carried a geometry-modeled response: confirm the detach first.
  if( !m_geometryModeled && m_orig && m_orig->ceeloResponse() )
  {
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dmw-detach-title"),
                                                             WString::tr("dmw-detach-body") );
    WPushButton *ok = dialog->addButton( WString::tr("dmw-detach-accept") );
    dialog->addButton( WString::tr("Cancel") );
    ok->clicked().connect( this, &DrfModifyWidget::apply );
    return;
  }//if( detaching a geometry-modeled response )

  // Geometry Modeled with pending edits and a possible generation: offer to regenerate first, so the
  //  attached response reflects the edits rather than the state it was generated from.
  const bool canGen = (m_mcTool && m_mcTool->generationReady());
  if( m_geometryModeled && m_changedSinceGenerate && canGen )
  {
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dmw-regen-title"),
                                                             WString::tr("dmw-regen-body") );
    WPushButton *regen = dialog->addButton( WString::tr("dmw-regen-accept") );
    WPushButton *useAnyway = dialog->addButton( WString::tr("dmw-regen-use-anyway") );
    dialog->addButton( WString::tr("Cancel") );
    regen->clicked().connect( this, [this](){
      m_applyAfterGenerate = true;   //handleResponseGenerated applies once the response lands
      handleGenerateResponse();
    } );
    useAnyway->clicked().connect( this, &DrfModifyWidget::apply );
    return;
  }//if( pending edits in Geometry-Modeled mode )

  apply();
}//requestApply()


void DrfModifyWidget::handleModeToggle()
{
  m_geometryModeled = (m_modeToggle && m_modeToggle->isChecked());
  if( m_mcTool )
    m_mcTool->setDisabled( !m_geometryModeled );  //Flat Disk greys the whole tool
  updateAnchorEditorVisibility();
  markEdited();
}//handleModeToggle()


bool DrfModifyWidget::geometryModeled() const
{
  return m_geometryModeled;
}


DrfModifyWidget::AnchorEditor DrfModifyWidget::activeAnchorEditor() const
{
  // Points win over the curve's representation whenever the DRF has any.
  //  - Geometry Modeled: the rows are what the Monte-Carlo response is grounded to, whatever
  //    representation the seed curve happened to use.
  //  - m_origHasPoints: a MakeDrf detector is an exp-of-log CURVE that also carries the measured
  //    points it was fitted from, and it is those points - not the coefficients - that produce its
  //    (consumed) node covariance via MeasuredDrfPoints::toEfficiencyUncert.  Dispatching on the
  //    curve form alone would hide them behind the coefficient editor.
  if( m_geometryModeled || m_origHasPoints )
    return AnchorEditor::Points;

  const shared_ptr<const DetectorEfficiencyCurve> curve = m_orig ? m_orig->efficiencyCurve()
                                                                 : nullptr;
  if( curve && curve->isValid() )
  {
    switch( curve->form() )
    {
      case DetectorPeakResponse::kExpOfLogPowerSeries:  return AnchorEditor::Coefficients;
      case DetectorPeakResponse::kFunctialEfficienyForm: return AnchorEditor::Formula;
      case DetectorPeakResponse::kEnergyEfficiencyPairs:
      case DetectorPeakResponse::kNumEfficiencyFnctForms:
        break;
    }//switch( curve->form() )
  }//if( curve && curve->isValid() )

  return AnchorEditor::Points;
}//activeAnchorEditor()


void DrfModifyWidget::updateAnchorEditorVisibility()
{
  const AnchorEditor editor = activeAnchorEditor();

  // The formula editor also shows the point table, for per-energy uncertainties; its Efficiency
  //  column is hidden there, the efficiency coming from the formula.
  const bool show_points = (editor != AnchorEditor::Coefficients);
  const bool show_coefs = (editor == AnchorEditor::Coefficients);
  const bool show_formula = (editor == AnchorEditor::Formula);

  if( m_pointsEditor )
    m_pointsEditor->setHidden( !show_points );
  if( m_coefEditor )
    m_coefEditor->setHidden( !show_coefs );
  if( m_formulaEditor )
    m_formulaEditor->setHidden( !show_formula );

  // A formula supplies the efficiency, so its rows carry uncertainty only.
  if( m_anchorTable )
    m_anchorTable->toggleStyleClass( "DrfNoEffCol", show_formula );

  // Only the equation and formula are written in a choosable energy unit; a point table shows the
  //  energies in the column itself.
  if( m_effUnitsRow )
    m_effUnitsRow->setHidden( !(show_coefs || show_formula) );

  if( m_anchorHelp )
  {
    const char *help_id = "dmw-anchor-help-points-nosrc";
    if( show_coefs )
      help_id = "dmw-anchor-help-coefs";
    else if( show_formula )
      help_id = "dmw-anchor-help-formula";
    else if( m_anchorHasSourceCol )
      help_id = "dmw-anchor-help-points";  //wording follows the columns actually rendered

    m_anchorHelp->setText( WString::tr( help_id ) );
  }//if( m_anchorHelp )
}//updateAnchorEditorVisibility()


void DrfModifyWidget::seedCoefCovFromUncert( const std::shared_ptr<const DetectorEfficiencyUncert> &uncert )
{
  m_coefCovMatrix.clear();

  const size_t n = m_coefEdits.size();
  if( uncert && n )
  {
    const vector<float> &cov = uncert->coefficientCovariance();
    if( cov.size() == (n*n) )
      m_coefCovMatrix.assign( begin(cov), end(cov) );
  }//if( uncert && n )

  // No stored covariance: start from the legacy per-coefficient uncertainties when there are any,
  //  so the user has something to edit rather than a grid of zeros.  This is only a starting point -
  //  it assumes the coefficients are independent, which for a log-power-series fit they emphatically
  //  are not, so it over-estimates the band until a real fit covariance is stored.
  if( m_coefCovMatrix.empty() && n )
  {
    m_coefCovMatrix.assign( n*n, 0.0 );

    const shared_ptr<const DetectorEfficiencyCurve> curve = m_orig ? m_orig->efficiencyCurve()
                                                                   : nullptr;
    const vector<float> legacy = curve ? curve->expOfLogPowerSeriesUncerts() : vector<float>{};
    for( size_t i = 0; (i < n) && (legacy.size() == n); ++i )
      m_coefCovMatrix[i*n + i] = static_cast<double>(legacy[i]) * legacy[i];
  }//if( m_coefCovMatrix.empty() && n )

  rebuildCovTable();
}//seedCoefCovFromUncert(...)


void DrfModifyWidget::rebuildCovTable()
{
  if( !m_covTable )
    return;

  m_covTable->clear();  //wipes all cells and the WLineEdits/WTexts they hold

  const size_t n = m_coefEdits.size();
  if( m_coefCovMatrix.size() != (n*n) )
    m_coefCovMatrix.assign( n*n, 0.0 );

  if( !n )
    return;

  auto coef_name = []( const size_t i ) -> WString {
    return WString::fromUTF8( "A" + std::to_string(i) );
  };

  // Header row: a corner label, then one coefficient name per matrix column.
  m_covTable->elementAt(0,0)->addNew<WText>( WString::tr("dmw-cov-corner-coef") );
  for( size_t j = 0; j < n; ++j )
    m_covTable->elementAt(0, static_cast<int>(j)+1)->addNew<WText>( coef_name(j) );

  for( size_t i = 0; i < n; ++i )
  {
    const int trow = static_cast<int>(i) + 1;

    // Column 0: which coefficient this row is - a static label, since the row set follows the
    //  coefficient boxes above.
    m_covTable->elementAt(trow,0)->addNew<WText>( coef_name(i) );

    const double Cii = m_coefCovMatrix[i*n + i];
    const double sigma_i = (Cii > 0.0) ? std::sqrt(Cii) : 0.0;

    for( size_t j = 0; j < n; ++j )
    {
      const int tcol = static_cast<int>(j) + 1;
      const double Cjj = m_coefCovMatrix[j*n + j];
      const double sigma_j = (Cjj > 0.0) ? std::sqrt(Cjj) : 0.0;

      if( j == i )
      {
        // Diagonal: the 1-sigma uncertainty of the coefficient itself (absolute, not a percent -
        //  these coefficients are logs, so a percentage of one means nothing).
        char buf[32];
        snprintf( buf, sizeof(buf), "%.4g", sigma_i );
        WLineEdit *e = m_covTable->elementAt(trow,tcol)->addNew<WLineEdit>();
        e->setTextSize( 8 );
        e->addStyleClass( "DrfCovDiag" );
        e->setText( WString::fromUTF8(buf) );
        e->changed().connect( this, [this,i,e](){ covSigmaChanged( i, e->text().toUTF8() ); } );
      }else
      {
        const double Cij = m_coefCovMatrix[i*n + j];
        const double rho = ((sigma_i > 0.0) && (sigma_j > 0.0)) ? (Cij/(sigma_i*sigma_j)) : 0.0;
        char buf[32];
        snprintf( buf, sizeof(buf), "%.3g", rho );
        WLineEdit *e = m_covTable->elementAt(trow,tcol)->addNew<WLineEdit>();
        e->setTextSize( 6 );
        e->setText( WString::fromUTF8(buf) );

        if( j > i )
        {
          // Upper triangle: correlation (editable), only meaningful when both sigmas are non-zero.
          e->setEnabled( (sigma_i > 0.0) && (sigma_j > 0.0) );
          e->changed().connect( this, [this,i,j,e](){ covRhoChanged( i, j, e->text().toUTF8() ); } );
        }else
        {
          // Lower triangle: read-only mirror.
          e->addStyleClass( "DrfCovMirror" );
          e->setEnabled( false );
        }
      }//if( diagonal ) / else
    }//for( size_t j = 0; j < n; ++j )
  }//for( size_t i = 0; i < n; ++i )
}//rebuildCovTable()


void DrfModifyWidget::addCoefficient( const float value )
{
  const size_t n = m_coefEdits.size();

  WContainerWidget *div = m_coefParams->addNew<WContainerWidget>();
  div->addStyleClass( "ParDiv" );
  div->addNew<WLabel>( WString::fromUTF8( "A" + std::to_string(n) ) );
  NativeFloatSpinBox *edit = div->addNew<NativeFloatSpinBox>();
  edit->setSpinnerHidden( true );
  edit->setWidth( 90 );
  edit->setValue( value );
  edit->valueChanged().connect( this, &DrfModifyWidget::markEdited );
  m_coefEdits.push_back( edit );

  // Grow the covariance to match, leaving the existing block bit-exact.
  const size_t nn = m_coefEdits.size();
  vector<double> grown( nn*nn, 0.0 );
  for( size_t i = 0; i < n; ++i )
    for( size_t j = 0; j < n; ++j )
      grown[i*nn + j] = m_coefCovMatrix[i*n + j];
  m_coefCovMatrix = grown;

  if( m_removeCoef )
    m_removeCoef->setEnabled( nn > 1 );

  m_renderFlags |= RenderActions::RebuildCovTable;
  scheduleRender();
}//addCoefficient(...)


void DrfModifyWidget::removeCoefficient()
{
  const size_t n = m_coefEdits.size();
  if( n < 2 )
    return;  //an equation needs at least one term

  m_coefParams->removeWidget( m_coefEdits.back()->parent() );
  m_coefEdits.pop_back();

  const size_t nn = m_coefEdits.size();
  vector<double> shrunk( nn*nn, 0.0 );
  for( size_t i = 0; i < nn; ++i )
    for( size_t j = 0; j < nn; ++j )
      shrunk[i*nn + j] = m_coefCovMatrix[i*n + j];
  m_coefCovMatrix = shrunk;

  if( m_removeCoef )
    m_removeCoef->setEnabled( nn > 1 );

  m_renderFlags |= RenderActions::RebuildCovTable;
  scheduleRender();
}//removeCoefficient()


void DrfModifyWidget::covSigmaChanged( const std::size_t i, const std::string &text )
{
  const size_t n = m_coefEdits.size();
  if( (i >= n) || (m_coefCovMatrix.size() != (n*n)) )
    return;

  double s_new = 0.0;
  try{ s_new = std::stod( text ); }
  catch( std::exception & ){ m_renderFlags |= RenderActions::RebuildCovTable; scheduleRender(); return; }

  if( s_new < 0.0 )
    s_new = 0.0;

  const double Cii = m_coefCovMatrix[i*n + i];
  const double s_old = (Cii > 0.0) ? std::sqrt(Cii) : 0.0;

  // Scale row/col i so every correlation ρ_ik is held fixed, then set the variance exactly - this
  //  leaves C_jk (j,k != i) bit-identical.  When s_old == 0 the correlations were undefined, so the
  //  off-diagonals stay zero and only the variance is set.
  if( s_old > 0.0 )
  {
    const double f = s_new / s_old;
    for( size_t k = 0; k < n; ++k )
    {
      if( k == i )
        continue;
      m_coefCovMatrix[i*n + k] *= f;
      m_coefCovMatrix[k*n + i] *= f;
    }
  }//if( s_old > 0.0 )
  m_coefCovMatrix[i*n + i] = s_new * s_new;

  m_renderFlags |= RenderActions::RebuildCovTable;
  markEdited();
}//covSigmaChanged(...)


void DrfModifyWidget::covRhoChanged( const std::size_t i, const std::size_t j, const std::string &text )
{
  const size_t n = m_coefEdits.size();
  if( (i >= n) || (j >= n) || (i == j) || (m_coefCovMatrix.size() != (n*n)) )
    return;

  double rho = 0.0;
  try{ rho = std::stod( text ); }
  catch( std::exception & ){ m_renderFlags |= RenderActions::RebuildCovTable; scheduleRender(); return; }
  if( rho < -1.0 ) rho = -1.0;
  if( rho >  1.0 ) rho =  1.0;

  const double sigma_i = std::sqrt( std::max(0.0, m_coefCovMatrix[i*n + i]) );
  const double sigma_j = std::sqrt( std::max(0.0, m_coefCovMatrix[j*n + j]) );
  const double Cij = rho * sigma_i * sigma_j;
  m_coefCovMatrix[i*n + j] = Cij;
  m_coefCovMatrix[j*n + i] = Cij;

  m_renderFlags |= RenderActions::RebuildCovTable;
  markEdited();
}//covRhoChanged(...)


bool DrfModifyWidget::validateFormula()
{
  if( !m_formulaText )
    return false;

  const string fcn = m_formulaText->text().toUTF8();
  bool valid = !fcn.empty();

  if( valid )
  {
    // The same trial-parse the Detector Select "Formula" tab uses: build a throwaway DRF and see
    //  whether the formula compiles at all.
    try
    {
      const float diam = m_orig ? m_orig->detectorDiameter() : 0.0f;
      DetectorPeakResponse trial( "temp", "temp" );
      trial.setIntrinsicEfficiencyFormula( fcn, (diam > 0.0f) ? diam : 1.0f,
                          equationEnergyUnits(),
                          m_orig ? m_orig->lowerEnergy() : 0.0f,
                          m_orig ? m_orig->upperEnergy() : 0.0f,
                          DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    }catch( std::exception & )
    {
      valid = false;
    }
  }//if( valid )

  // An empty box is "not filled in yet", not "wrong" - dont flag it red.
  m_formulaText->toggleStyleClass( "Wt-invalid", !valid && !fcn.empty() );

  return valid;
}//validateFormula()


float DrfModifyWidget::equationEnergyUnits() const
{
  return (m_effEnergyUnits && (m_effEnergyUnits->currentIndex() == 1))
             ? static_cast<float>(PhysicalUnits::MeV)
             : static_cast<float>(PhysicalUnits::keV);
}//equationEnergyUnits()


void DrfModifyWidget::applyCoefficientEdits( DetectorPeakResponse &working )
{
  if( m_coefEdits.empty() )
    return;

  if( anchorsMatchSeed() )
    return;  //nothing touched; leave the curve and its uncertainty bit-identical

  vector<float> coefs;
  coefs.reserve( m_coefEdits.size() );
  for( NativeFloatSpinBox * const edit : m_coefEdits )
    coefs.push_back( edit->value() );

  const size_t n = coefs.size();
  vector<float> coefCov;
  if( m_coefCovMatrix.size() == (n*n) )
  {
    coefCov.reserve( n*n );
    for( const double val : m_coefCovMatrix )
      coefCov.push_back( static_cast<float>(val) );

    // All-zero is "no covariance given", not a covariance of zero.
    bool any = false;
    for( const float val : coefCov )
      any = (any || (val != 0.0f));
    if( !any )
      coefCov.clear();
  }//if( m_coefCovMatrix.size() == (n*n) )

  // The equation's diagonal doubles as the legacy per-coefficient uncertainty, which is what gets
  //  written to the app-URL and the older DB fields.
  vector<float> uncerts;
  if( !coefCov.empty() )
  {
    uncerts.resize( n );
    for( size_t i = 0; i < n; ++i )
      uncerts[i] = std::sqrt( std::max( 0.0f, coefCov[i*n + i] ) );
  }

  try
  {
    auto curve = make_shared<DetectorEfficiencyCurve>();
    curve->setFromExpOfLogPowerSeries( coefs, uncerts, equationEnergyUnits() );

    if( !coefCov.empty() )
    {
      const shared_ptr<const DetectorEfficiencyUncert> existing = working.efficiencyUncert();
      auto uncert = existing ? make_shared<DetectorEfficiencyUncert>( *existing )
                             : make_shared<DetectorEfficiencyUncert>();
      // For an exp-of-log curve this IS the uncertainty the fits see: DetectorEfficiencyCurve::
      //  fracCovariance propagates it as J*Sigma*J^T.  The node covariance is left untouched, so a
      //  detector that also carries one (a MakeDrf response, from its measured points) keeps it as
      //  provenance without the two being combined.
      uncert->setCoefficientCovariance( coefCov );
      curve->setUncertainty( uncert->isEmpty() ? nullptr : uncert );
    }else if( m_orig )
    {
      curve->setUncertainty( m_orig->efficiencyUncert() );
    }

    working.replaceEfficiencyCurve( curve );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("dmw-err-coefs-invalid").arg(e.what()), WarningWidget::WarningMsgHigh );
  }
}//applyCoefficientEdits(...)


void DrfModifyWidget::applyFormulaEdits( DetectorPeakResponse &working )
{
  if( !m_formulaText )
    return;

  if( anchorsMatchSeed() )
    return;  //nothing touched; leave the curve and its uncertainty bit-identical

  const string fcn = m_formulaText->text().toUTF8();
  if( fcn.empty() )
    return;

  try
  {
    auto curve = make_shared<DetectorEfficiencyCurve>();
    curve->setFromFormula( fcn, equationEnergyUnits() );
    // Without the correlation control there are no rows to rebuild from, so carry the existing
    //  uncertainty across rather than dropping it on a formula-only edit.
    curve->setUncertainty( m_uncertOptions ? buildUncertFromRows() : working.efficiencyUncert() );
    working.replaceEfficiencyCurve( curve );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("dmw-err-formula-invalid").arg(e.what()), WarningWidget::WarningMsgHigh );
  }
}//applyFormulaEdits(...)


bool DrfModifyWidget::anchorsMatchSeed() const
{
  if( !m_seedState )
    return false;

  const shared_ptr<ToolState> now = currentState();

  return (now->anchors == m_seedState->anchors)
         && (now->anchorRefDistance == m_seedState->anchorRefDistance)
         && (now->anchorDefaultUncert == m_seedState->anchorDefaultUncert)
         && (now->anchorCorrLength == m_seedState->anchorCorrLength)
         && (now->coefficients == m_seedState->coefficients)
         && (now->coefCovMatrix == m_seedState->coefCovMatrix)
         && (now->formula == m_seedState->formula)
         && (now->efficiencyEnergyUnits == m_seedState->efficiencyEnergyUnits);
}//anchorsMatchSeed()


void DrfModifyWidget::handleGenerateResponse()
{
  if( !m_mcTool || !m_geometryModeled )
    return;

  // Seed the MC tool with the live edits, CeeLo detached (buildWorkingDrf(false)), so the manual
  //  points/covariance ground the regeneration rather than being overridden by an attached response.
  const shared_ptr<DetectorPeakResponse> seed = buildWorkingDrf( false );
  m_mcTool->setSeedDrf( seed );

  // For a raw-measured-points DRF, the regeneration derives its uncertainty from those points, not
  //  from a hand-edited matrix - say so once, so an edited covariance is not silently discarded.
  if( m_origHasPoints && m_anchorIsAbsolute )
    passMessage( WString::tr("dmw-regen-from-points-note"), WarningWidget::WarningMsgInfo );

  m_mcTool->startGeneration();
  updateGenerateButton();
}//handleGenerateResponse()


void DrfModifyWidget::updateGenerateButton()
{
  if( !m_generateBtn )
    return;

  const bool haveResp = (m_mcTool && m_mcTool->generatedResponse());
  const bool canGen = (m_mcTool && m_mcTool->generationReady());

  // Shown only in Geometry-Modeled mode, once there is something to (re)generate from or a response
  //  to refresh.  Enabled whenever the geometry is complete enough to run and there is work to do:
  //  the first generation (no response yet - e.g. right after switching to Geometry Modeled with
  //  valid geometry) or a regeneration when an edit is pending.
  const bool show = m_geometryModeled && (canGen || haveResp);
  m_generateBtn->setHidden( !show );
  m_generateBtn->setEnabled( show && canGen && (!haveResp || m_changedSinceGenerate) );
}//updateGenerateButton()


void DrfModifyWidget::handleResponseGenerated( std::shared_ptr<ceelo::DetectorResponse> response )
{
  // A fresh response clears staleness.  The MC tool emits userChanged right after this, which routes
  //  to markEdited - suppress that one so it does not immediately re-flag the response stale.
  m_changedSinceGenerate = false;
  m_suppressNextEditMark = true;
  updateGenerateButton();

  if( m_applyAfterGenerate )
  {
    m_applyAfterGenerate = false;
    if( response )
      apply();   //the regenerate-then-use flow: the response now reflects the edits
  }
}//handleResponseGenerated(...)


// ---------------------------------------------------------------------------

DrfModifyWidget *DrfModifyWindow::tool()
{
  return m_tool;
}


DrfModifyWindow::DrfModifyWindow( InterSpec *viewer,
                                  std::shared_ptr<const DetectorPeakResponse> drf )
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
    auto toolOwned = make_unique<DrfModifyWidget>( viewer, drf );
    m_tool = toolOwned.get();
    stretcher()->addWidget( std::move(toolOwned), 0, 0 );
  }
  stretcher()->setContentsMargins( 0, 0, 0, 0 );

  AuxWindow::addHelpInFooter( footer(), "modify-drf" );

  WPushButton *cancel = addCloseButtonToFooter( WString::tr("Cancel") );
  cancel->clicked().connect( this, &AuxWindow::hide );

  WPushButton *use = footer()->addNew<WPushButton>( WString::tr("dmw-use-btn") );
  use->clicked().connect( m_tool, &DrfModifyWidget::requestApply );

  // A DRF with no efficiency curve of its own - a Detector.dat or .detx imported
  //  for its geometry alone - is not usable until the Monte Carlo has produced
  //  one.  Hold "Use" closed until it has, rather than letting the dialog be
  //  dismissed with an invalid detector.
  if( m_tool->needsMcResponse() )
  {
    use->disable();
    HelpSystem::attachToolTipOn( use, WString::tr("dmw-tt-use-needs-mc"), true );
    m_tool->mcResponseAvailable().connect( std::bind( [use]( const bool have ){
      use->setEnabled( have );
    }, std::placeholders::_1 ) );
  }//if( m_tool->needsMcResponse() )

  show();
  resizeToFitOnScreen();
  centerWindow();
}//DrfModifyWindow constructor
