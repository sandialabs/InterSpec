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

#include <set>
#include <array>
#include <cmath>
#include <cctype>
#include <memory>
#include <string>
#include <numeric>
#include <functional>
#include <vector>
#include <algorithm>

#include <Wt/Utils.h>
#include <Wt/WText.h>
#include <Wt/WMenu.h>
#include <Wt/WTable.h>
#include <Wt/WTableRow.h>
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
#include "SpecUtils/StringAlgo.h"

#include "io/DetectorResponse.h"

#include "InterSpec/DrfChart.h"
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
#include "InterSpec/DrfModifyCalc.h"
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


/** Parses a line edit as a number; `ok` is set false (and the value left alone) on a non-empty cell
 that is not a number, so a typo is reported rather than silently skipping the row.
 */
bool read_double( const WLineEdit * const edit, double &value, bool &was_blank )
{
  was_blank = true;
  if( !edit )
    return true;

  const string text = edit->text().toUTF8();
  if( text.empty() )
    return true;

  was_blank = false;

  // Whole-cell parse: `std::stod` alone would take "0.4x" as 0.4, which is a silent
  //  reinterpretation of something the user can see is wrong.
  size_t used = 0;
  try
  {
    value = std::stod( text, &used );
  }catch( std::exception & )
  {
    return false;
  }

  while( (used < text.size()) && std::isspace( static_cast<unsigned char>(text[used]) ) )
    ++used;

  if( used != text.size() )
    return false;

  return !(std::isnan(value) || std::isinf(value));
}//read_double(...)
}//namespace


DrfModifyWidget::DrfModifyWidget( InterSpec *viewer,
                                  std::shared_ptr<const DetectorPeakResponse> drf )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_orig( drf ),
    m_seedPoints( drf ? drf->measuredPoints() : nullptr ),
    m_editor( DrfModifyCalc::AnchorEditor::CurvePairs ),
    m_tabMenu( nullptr ),
    m_tabStack( nullptr ),
    m_name( nullptr ),
    m_description( nullptr ),
    m_generalChart( nullptr ),
    m_infoTable( nullptr ),
    m_infoValues{},
    m_generalTabItem( nullptr ),
    m_generalStale( true ),
    m_mcTool( nullptr ),
    m_fwhmTool( nullptr ),
    m_fwhmTabItem( nullptr ),
    m_geomTabItem( nullptr ),
    m_modeToggle( nullptr ),
    m_geometryModeled( false ),
    m_anchorHelp( nullptr ),
    m_responseNote( nullptr ),
    m_uncertSummary( nullptr ),
    m_pointsEditor( nullptr ),
    m_coefEditor( nullptr ),
    m_formulaEditor( nullptr ),
    m_anchorTable( nullptr ),
    m_anchorTableWrap( nullptr ),
    m_addAnchor( nullptr ),
    m_removeAnchor( nullptr ),
    m_anchorRefDistance( nullptr ),
    m_anchors(),
    m_uncertOptions( nullptr ),
    m_corrInertNote( nullptr ),
    m_anchorEnergyUnits( static_cast<float>(PhysicalUnits::keV) ),
    m_coefEdits(),
    m_coefParams( nullptr ),
    m_covTable( nullptr ),
    m_addCoef( nullptr ),
    m_removeCoef( nullptr ),
    m_coefSigmas(),
    m_coefRho(),
    m_covWarning( nullptr ),
    m_covPlaceholderNote( nullptr ),
    m_coefCovIsPlaceholder( false ),
    m_coefCovTouched( false ),
    m_formulaText( nullptr ),
    m_effEnergyUnits( nullptr ),
    m_effUnitsRow( nullptr ),
    m_seedState( nullptr ),
    m_generateBtn( nullptr ),
    m_generatedFromFingerprint( 0 ),
    m_pendingSeedFingerprint( 0 ),
    m_exportNote( nullptr ),
    m_generateHint( nullptr ),
    m_applyAfterGenerationId( -1 ),
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

  const shared_ptr<const DetectorEfficiencyCurve> eff_curve
      = m_orig ? m_orig->efficiencyCurve() : nullptr;

  // Which editor this DRF gets, and hence what its rows mean, comes from what the DRF carries - see
  //  DrfModifyCalc::editorForDrf.  The rows are then seeded from whatever that editor's apply writes
  //  back, which is the whole point: seeding a Create-DRF detector's rows from its ABSOLUTE measured
  //  points while the apply treated them as the INTRINSIC curve is what made a detector come out
  //  ~200x too insensitive.
  if( m_orig )
    m_editor = DrfModifyCalc::editorForDrf( *m_orig );

  const bool measured_rows = DrfModifyCalc::editorUsesMeasuredPoints( m_editor );

  // The Energy column of a pairs curve is in the curve's own units - which a GADRAS CSV may set to
  //  MeV.  Measured points and covariance nodes are keV by contract whatever units the curve's
  //  equation uses (a MakeDrf detector can be an MeV equation fitted to keV points), so only the
  //  pairs editor, whose own numbers fill the column, takes its units from the curve.
  if( (m_editor == DrfModifyCalc::AnchorEditor::CurvePairs) && eff_curve && eff_curve->isValid()
      && (eff_curve->energyUnits() > 0.0f) )
    m_anchorEnergyUnits = eff_curve->energyUnits();

  // Fixed-geometry DRFs have no geometry to model, so no Geom & MC tab and no mode toggle.  A
  //  far-field DRF starts Geometry Modeled iff it already carries a Monte-Carlo response, knows its
  //  physical shape (an ANGLE / Detector.dat import, or a geometry saved with it), or is a
  //  geometry-only import that has no efficiency yet, and so needs one.
  const bool fixed_geom = (m_orig && m_orig->isFixedGeometry());
  const bool has_geom_tab = !fixed_geom;
  m_geometryModeled = has_geom_tab
      && (!m_orig || m_orig->ceeloResponse() || m_orig->geometry() || !m_orig->isValid());

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
    m_name->changed().connect( this, &DrfModifyWidget::markEditedNoRegen );
    nameLabel->setBuddy( m_name );

    WLabel *descLabel = idGrid->addNew<WLabel>( WString::tr("Description") );
    m_description = idGrid->addNew<WTextArea>();
    m_description->setColumns( 40 );
    m_description->setRows( 3 );
    m_description->setText( WString::fromUTF8( m_orig ? m_orig->description() : string() ) );
    m_description->changed().connect( this, &DrfModifyWidget::markEditedNoRegen );
    descLabel->setBuddy( m_description );

    // What the detector carries and how it responds: a live chart of the efficiency + FWHM, and a
    //  summary of what it has (measured points, geometry, Monte-Carlo support, uncertainty, ...).
    //  Both follow every edit on the other tabs - see refreshGeneralTab().
    m_generalChart = panel->addNew<DrfChart>();
    m_generalChart->addStyleClass( "DrfModifyGeneralChart" );
    m_generalChart->setShowFwhm( true );

    buildInfoTable( panel );

    m_generalTabItem = m_tabMenu->addItem( WString::tr("dmw-tab-name"), std::move(panelOwned) );
    make_item_selectable( m_tabMenu, m_generalTabItem );
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

    // We drive generation from the footer "Generate Response" button, so hide the tool's own
    //  redundant generate button in Location Support.
    m_mcTool->setGenerateButtonHidden( true );

    // ContentLoading::Eager: this tool posts work to a worker thread and gets back to itself with
    //  `findById(...)`; a Lazy tab parks its contents in `WMenuItem::uContents_`, outside the
    //  widget tree, where `findById` can not see it - see the FWHM tab below.
    m_geomTabItem = m_tabMenu->addItem( WString::tr("dmw-tab-geometry"), std::move(panelOwned),
                                        ContentLoading::Eager );
    make_item_selectable( m_tabMenu, m_geomTabItem );

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
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    panel->addStyleClass( "DrfModifyPanel DrfModifyAnchorPanel" );

    // Wording is set to match the visible editor in updateAnchorEditorVisibility().
    m_anchorHelp = panel->addNew<WText>();
    m_anchorHelp->setInline( false );

    // While a geometry-modeled response is attached it - not the numbers below - answers every
    //  efficiency and uncertainty query.  Saying so is the difference between editing a curve and
    //  editing the detector.
    m_responseNote = panel->addNew<WText>( WString::tr("dmw-anchor-note-response") );
    m_responseNote->setInline( false );
    m_responseNote->addStyleClass( "DrfModifyNote" );

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

      m_covWarning = m_coefEditor->addNew<WText>();
      m_covWarning->setInline( false );
      m_covWarning->addStyleClass( "DrfModifyWarning" );
      m_covWarning->hide();

      // Says when the matrix on screen was manufactured from the stored per-coefficient sigmas
      //  rather than read from a covariance this DRF carries - the two look identical otherwise, and
      //  the manufactured one assumes an independence that a log-power-series fit never has.
      m_covPlaceholderNote = m_coefEditor->addNew<WText>( WString::tr("dmw-coef-cov-assumed-note") );
      m_covPlaceholderNote->setInline( false );
      m_covPlaceholderNote->addStyleClass( "DrfModifyNote" );
      m_covPlaceholderNote->hide();

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
      if( m_seedPoints )
      {
        for( const MeasuredEffPoint &p : m_seedPoints->points() )
        {
          if( p.distance > 0.0f ){ refDistCm = p.distance / PhysicalUnits::cm; break; }
        }
      }//if( m_seedPoints )
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
      // Only the absolute reference curve is anchored at one distance; a re-fit takes each point's
      //  own distance, and an intrinsic curve has none.
      distRow->setHidden( m_editor != DrfModifyCalc::AnchorEditor::AbsolutePoints );

      m_anchorTableWrap = m_pointsEditor->addNew<WContainerWidget>();
      m_anchorTableWrap->addStyleClass( "DrfModifyAnchorTableWrap" );
      m_anchorTable = m_anchorTableWrap->addNew<WTable>();
      m_anchorTable->addStyleClass( "DrfModifyAnchorTable" );
      {
        // Units in the header, because the column is not always keV (a GADRAS CSV curve may be MeV).
        const bool in_mev = (m_anchorEnergyUnits > 10.0f);
        const WString energy_units = WString::tr( in_mev ? "dmw-anchor-mev" : "dmw-anchor-kev" );

        // What the efficiency column means, which is what the apply must agree with: absolute at a
        //  distance for the measured-point editors and an absolute reference curve, intrinsic for a
        //  far-field curve of its own, and just "efficiency" for a fixed-geometry DRF.
        const char *eff_id = "dmw-anchor-eff";
        if( measured_rows
            || (m_orig && (m_orig->geometryType()
                           == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute)) )
          eff_id = "dmw-anchor-eff-abs";
        else if( m_orig && !m_orig->isFixedGeometry() )
          eff_id = "dmw-anchor-eff-intrinsic";

        int col = 0;
        m_anchorTable->elementAt(0,col++)->addNew<WText>(
                                        WString::tr("dmw-anchor-energy").arg(energy_units) );
        WTableCell * const effHeader = m_anchorTable->elementAt(0,col++);
        effHeader->addNew<WText>( WString::tr(eff_id) );
        effHeader->addStyleClass( "DrfEffCol" );
        // The two components, named the same way the control below them and the help text are:
        //  measured points split into statistical and per-source certificate, everything else into
        //  the part that is independent between energies and the part that is not.
        m_anchorTable->elementAt(0,col++)->addNew<WText>(
                    WString::tr( measured_rows ? "dmw-anchor-stat" : "dmw-anchor-uncorr" ) );
        m_anchorTable->elementAt(0,col++)->addNew<WText>(
                    WString::tr( measured_rows ? "dmw-anchor-cert" : "dmw-anchor-corr" ) );
        if( measured_rows )
        {
          m_anchorTable->elementAt(0,col++)->addNew<WText>( WString::tr("dmw-anchor-source") );
          // Per-point distance: characterization sources need not all sit at one distance
          m_anchorTable->elementAt(0,col++)->addNew<WText>( WString::tr("dmw-anchor-dist") );
        }
      }

      WContainerWidget *btns = m_pointsEditor->addNew<WContainerWidget>();
      m_addAnchor = btns->addNew<WPushButton>( WString::tr("dmw-anchor-add") );
      m_addAnchor->addStyleClass( "LinkBtn" );
      m_addAnchor->clicked().connect( this, [this](){ addAnchorRow( 0.0f, 0.0f, 0.0f, 0.0f, string() ); markEdited(); } );
      m_removeAnchor = btns->addNew<WPushButton>( WString::tr("dmw-anchor-remove") );
      m_removeAnchor->addStyleClass( "LinkBtn" );
      m_removeAnchor->clicked().connect( this, [this](){ removeAnchorRow(); markEdited(); } );

      if( measured_rows && m_seedPoints )
      {
        // The raw points, as they are: absolute efficiency at each point's own source distance.  The
        //  row remembers which point it came from, so editing its energy does not lose the point's
        //  peak area, live time, file name or distance uncertainty.
        const vector<MeasuredEffPoint> &points = m_seedPoints->points();
        for( size_t i = 0; i < points.size(); ++i )
        {
          const MeasuredEffPoint &p = points[i];
          addAnchorRow( p.energy, p.efficiency, p.fracStatUncert, p.fracCertUncert, p.sourceKey,
                        p.distance, static_cast<int>(i) );
        }
      }else if( (m_editor == DrfModifyCalc::AnchorEditor::CurvePairs) && eff_curve
                && eff_curve->isValid()
                && (eff_curve->form() == DetectorPeakResponse::kEnergyEfficiencyPairs) )
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
      }//if( measured rows ) / else if( pairs ) / else if( have a covariance )
      m_removeAnchor->setEnabled( !m_anchors.empty() );

      // How the correlated column is correlated across energy.  Shared with the .ecc import
      //  dialogs, so the mode -> correlation-length mapping and the example table live in one
      //  place, and the Modify dialog speaks the same language the import did.  Measured points do
      //  not get one: their correlated part is blocked per source, which a single length cannot say.
      if( !measured_rows )
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
        //  was set directly", so fall back to fully correlated - the "was this editor edited" guard
        //  is what keeps that guess from rewriting a covariance the user never touched.
        const double seed_len = (orig_uncert && (orig_uncert->correlationLength() > 0.0))
                       ? orig_uncert->correlationLength()
                       : DetectorEfficiencyUncert::sm_fullyCorrelatedLength;
        m_uncertOptions->setCorrelationLength( seed_len );
        m_uncertOptions->changed().connect( this, &DrfModifyWidget::markEdited );

        m_corrInertNote = m_pointsEditor->addNew<WText>( WString::tr("dmw-corr-inert-note") );
        m_corrInertNote->setInline( false );
        m_corrInertNote->addStyleClass( "DrfModifyNote" );
        m_corrInertNote->hide();
      }//if( !measured_rows )

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

    // What this detector reports, from the same call the activity fit uses - including how much of
    //  it is an envelope nobody measured.  See refreshUncertSummary().
    m_uncertSummary = panel->addNew<WText>();
    m_uncertSummary->setInline( false );
    m_uncertSummary->addStyleClass( "DrfModifyNote DrfModifyUncertSummary" );
    HelpSystem::attachToolTipOn( m_uncertSummary, WString::tr("dmw-tt-uncert-summary"), true );

    updateAnchorEditorVisibility();

    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-anchor"), std::move(panelOwned) );
    make_item_selectable( m_tabMenu, item );
  }

  m_tabMenu->itemSelected().connect( this, &DrfModifyWidget::handleTabSelected );

  // Both child tools hand their state changes here, so one undo step covers the whole dialog.
  //  m_mcTool is null for a fixed-geometry DRF (no Geom & MC tab).
  if( m_mcTool )
  {
    m_mcTool->userChanged().connect( this, &DrfModifyWidget::markEdited );
    m_mcTool->userChangedNoRegen().connect( this, &DrfModifyWidget::markEditedNoRegen );
    m_mcTool->responseGenerated().connect( this, &DrfModifyWidget::handleResponseGenerated );

    // EVERY generation - the footer button, and the automatic one a geometry change triggers for the
    //  instant curve-transfer method - anchors on the live edits, with the response detached so the
    //  manual points/covariance drive the grounding rather than being overridden by it.
    m_mcTool->setSeedProvider( [this]() -> shared_ptr<const DetectorPeakResponse> {
      vector<DrfModifyCalc::Problem> problems;
      const shared_ptr<DetectorPeakResponse> seed = buildWorkingDrf( false, problems );

      // An edit that could not be applied means this seed is NOT what the user is looking at, so
      //  whatever is generated from it must not come back looking current: a fingerprint of 0 never
      //  matches, so the response stays stale until the edit is fixed and it is rebuilt.
      m_pendingSeedFingerprint = (seed && !DrfModifyCalc::anyBlocking(problems))
                                 ? DrfModifyCalc::seedFingerprint( *seed ) : 0;
      return seed;
    } );
  }//if( m_mcTool )
  if( m_fwhmTool )
  {
    m_fwhmTool->setOwnerHandlesUndoRedo( true );
    m_fwhmTool->stateChanged().connect( this, &DrfModifyWidget::markEditedNoRegen );
  }//if( m_fwhmTool )

  m_tabMenu->select( 0 );

  // Footer: the export note on the left, and (Geometry-Modeled only) a "Generate Response" button on
  //  the right that regenerates the MC response from the live edits.
  auto footerOwned = make_unique<WContainerWidget>();
  WContainerWidget *footerRow = footerOwned.get();
  footerRow->addStyleClass( "DrfModifyFooterRow" );

  m_exportNote = footerRow->addNew<WText>( WString::tr("dmw-export-note") );
  m_exportNote->addStyleClass( "DrfModifyNote" );
  m_exportNote->setInline( false );

  if( m_mcTool )
  {
    // Why generation is blocked (an incomplete geometry, or one still guessed from the diameter);
    //  takes the export tip's place while it applies, so the reason sits right by the button.
    //  Plain text: the reason can quote what the user typed into a material field.
    m_generateHint = footerRow->addNew<WText>();
    m_generateHint->setTextFormat( TextFormat::Plain );
    m_generateHint->addStyleClass( "DrfModifyNote DrfModifyGenerateHint" );
    m_generateHint->setInline( false );
    m_generateHint->hide();

    m_generateBtn = footerRow->addNew<WPushButton>( WString::tr("dmw-generate-btn") );
    m_generateBtn->addStyleClass( "DrfModifyGenerateBtn" );
    m_generateBtn->clicked().connect( this, [this](){ handleGenerateResponse(); } );
  }//if( m_mcTool )

  if( narrow_layout )
    layout->addWidget( std::move(footerOwned), 2, 0 );        //below menu-row and stack-row
  else
    layout->addWidget( std::move(footerOwned), 1, 0, 1, 2 );  //below, spanning menu and stack columns

  // What the Anchor tab opened with, so the per-editor edited-checks can tell "the user changed
  //  nothing" from "the user re-typed the same numbers" - an untouched editor leaves the DRF alone.
  m_seedState = currentState();

  // The response the DRF arrived with was built from the DRF as it arrived.  Taking the fingerprint
  //  of a freshly built (un-edited) working copy, rather than of m_orig, means any difference from
  //  here on is a real edit and not a round-trip of the geometry form through the same code path.
  {
    vector<DrfModifyCalc::Problem> problems;
    const shared_ptr<DetectorPeakResponse> seed = buildWorkingDrf( false, problems );
    m_generatedFromFingerprint = seed ? DrfModifyCalc::seedFingerprint( *seed ) : 0;
  }

  updateGenerateButton();
  refreshUncertSummary();
  refreshGeneralTab();
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

  if( item && (item == m_generalTabItem) && m_generalStale )
    refreshGeneralTab();

  // The run-time estimate is timed on a short test Monte Carlo the first time the tab is looked at
  //  (and after the geometry changes); a tab never opened costs nothing.
  if( item && (item == m_geomTabItem) && m_mcTool )
    m_mcTool->scheduleTimeCalibration();
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


bool DrfModifyWidget::anchorHasSourceCols() const
{
  return DrfModifyCalc::editorUsesMeasuredPoints( m_editor );
}


void DrfModifyWidget::addAnchorRow( const float energy, const float efficiency,
                                    const float fracStatUncert, const float fracCertUncert,
                                    const std::string &sourceKey, const float distance,
                                    const int seedIndex )
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
  //  apply.  Only the source/distance columns are conditional (see anchorHasSourceCols); the
  //  efficiency one is built always and hidden by CSS for a formula curve.
  int col = 0;
  AnchorRow r;
  r.energy = make_edit( col++, energy, false );
  r.eff    = make_edit( col++, efficiency, false );
  r.stat   = make_edit( col++, fracStatUncert, true );
  r.cert   = make_edit( col++, fracCertUncert, true );
  r.source = nullptr;
  r.dist = nullptr;
  r.seedIndex = seedIndex;
  if( anchorHasSourceCols() )
  {
    r.source = m_anchorTable->elementAt(row,col++)->addNew<WLineEdit>();
    r.source->setTextSize( 10 );
    r.source->setText( WString::fromUTF8(sourceKey) );
    r.source->changed().connect( this, &DrfModifyWidget::markEdited );

    r.dist = m_anchorTable->elementAt(row,col++)->addNew<WLineEdit>();
    r.dist->setTextSize( 6 );
    if( distance > 0.0f )
    {
      char buf[32];
      snprintf( buf, sizeof(buf), "%.4g", distance / PhysicalUnits::cm );
      r.dist->setText( buf );
    }
    r.dist->changed().connect( this, &DrfModifyWidget::markEdited );
  }//if( anchorHasSourceCols() )
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


bool DrfModifyWidget::collectPointRows( std::vector<DrfModifyCalc::PointRow> &rows,
                                        std::vector<DrfModifyCalc::Problem> &problems ) const
{
  rows.clear();

  if( !m_anchorTable )
    return true;

  bool all_ok = true;
  for( size_t i = 0; i < m_anchors.size(); ++i )
  {
    const AnchorRow &r = m_anchors[i];
    const int row_number = static_cast<int>(i) + 1;

    DrfModifyCalc::PointRow row;
    row.rowNumber = row_number;
    row.seedIndex = r.seedIndex;

    bool blank = true;
    bool ok = true;

    ok = read_double( r.energy, row.energy, blank ) && ok;
    const bool energy_blank = blank;

    double efficiency = 0.0;
    ok = read_double( r.eff, efficiency, blank ) && ok;
    const bool eff_blank = blank;
    row.efficiency = efficiency;

    // A wholly blank row is the empty row an "Add point" left behind, not an error.
    if( energy_blank && eff_blank )
    {
      bool any_other = false;
      double scratch = 0.0;
      read_double( r.stat, scratch, blank );  any_other = any_other || !blank;
      read_double( r.cert, scratch, blank );  any_other = any_other || !blank;
      read_double( r.dist, scratch, blank );  any_other = any_other || !blank;
      if( !any_other && (!r.source || r.source->text().empty()) )
        continue;
    }//if( energy and efficiency are both blank )

    // A row that is blank in ONLY one of the two columns that make a point is not a point, and the
    //  calc layer would quietly drop it - which is the silent loss this function exists to prevent.
    //  (A formula curve's rows carry no efficiency, so only the energy is required there.)
    const bool need_efficiency = (m_editor != DrfModifyCalc::AnchorEditor::Formula);
    if( energy_blank || (need_efficiency && eff_blank) )
    {
      problems.push_back( DrfModifyCalc::Problem( "dmw-err-incomplete-row",
                                                  std::to_string(row_number) ) );
      all_ok = false;
      continue;
    }//if( only one of energy/efficiency was given )

    double stat_percent = 0.0;
    if( read_double( r.stat, stat_percent, blank ) )
    {
      if( !blank )
        row.fracStat = 0.01 * stat_percent;
    }else
    {
      ok = false;
    }

    double cert_percent = 0.0;
    if( read_double( r.cert, cert_percent, blank ) )
    {
      if( !blank )
        row.fracCert = 0.01 * cert_percent;
    }else
    {
      ok = false;
    }

    if( r.source )
      row.sourceKey = r.source->text().toUTF8();

    if( r.dist )
    {
      double dist_cm = 0.0;
      if( read_double( r.dist, dist_cm, blank ) )
      {
        if( !blank && (dist_cm > 0.0) )
          row.distance = dist_cm * PhysicalUnits::cm;
      }else
      {
        ok = false;
      }
    }//if( r.dist )

    if( !ok )
    {
      // Named, not skipped: silently dropping a malformed row is how an apply ends up writing fewer
      //  points than the user can see.
      problems.push_back( DrfModifyCalc::Problem( "dmw-err-bad-row",
                                                  std::to_string(row_number) ) );
      all_ok = false;
      continue;
    }

    // Values that parse but cannot describe a point.  The calc layer drops these, so if they are not
    //  caught here the apply quietly writes fewer points than are on screen.
    if( (row.energy <= 0.0) || (need_efficiency && (row.efficiency <= 0.0)) )
    {
      problems.push_back( DrfModifyCalc::Problem( "dmw-err-nonpositive-row",
                                                  std::to_string(row_number) ) );
      all_ok = false;
      continue;
    }

    if( (row.fracStat < 0.0) || (row.fracCert < 0.0) )
    {
      problems.push_back( DrfModifyCalc::Problem( "dmw-err-negative-uncert",
                                                  std::to_string(row_number) ) );
      all_ok = false;
      continue;
    }

    // A far-field point's absolute efficiency means nothing without the distance it was taken at,
    //  and for the re-fit editor there is no reference distance to fall back on.
    if( r.dist && (m_editor == DrfModifyCalc::AnchorEditor::RefitPoints)
        && (row.distance <= 0.0) && m_orig && !m_orig->isFixedGeometry() )
    {
      problems.push_back( DrfModifyCalc::Problem( "dmw-err-row-needs-distance",
                                                  std::to_string(row_number) ) );
      all_ok = false;
      continue;
    }

    rows.push_back( row );
  }//for( size_t i = 0; i < m_anchors.size(); ++i )

  return all_ok;
}//collectPointRows(...)


DrfModifyCalc::AnchorOptions DrfModifyWidget::anchorOptions() const
{
  DrfModifyCalc::AnchorOptions options;
  options.energyUnits = m_anchorEnergyUnits;

  // Only the absolute-reference editor shows (or means) a reference distance.  Handing one to the
  //  re-fit editor, whose distance column is per point, would let a blank cell silently take another
  //  point's distance - and would stop a distance-less row ever being refused.
  options.hasRowTable = (m_anchorTable != nullptr);

  if( m_anchorRefDistance && (m_editor == DrfModifyCalc::AnchorEditor::AbsolutePoints) )
  {
    bool blank = true;
    double dist_cm = 0.0;
    if( read_double( m_anchorRefDistance, dist_cm, blank ) && !blank && (dist_cm > 0.0) )
      options.refDistance = dist_cm * PhysicalUnits::cm;
  }

  if( m_uncertOptions )
    options.corrLength = m_uncertOptions->effectiveCorrLength();

  options.equationTerms = static_cast<int>( m_coefEdits.size() );

  return options;
}//anchorOptions()


bool DrfModifyWidget::applyAnchorTab( DetectorPeakResponse &working,
                                      std::vector<DrfModifyCalc::Problem> &problems )
{
  // Only the editor that is showing can apply, and only when it was actually edited.  Both matter:
  //  dispatching on anything that the Flat-Disk / Geometry-Modeled toggle can change would let a
  //  toggle send the apply down an editor the user never typed into, and rebuilding an untouched
  //  editor's numbers would replace a covariance the correlated+diagonal model cannot reproduce.
  switch( m_editor )
  {
    case DrfModifyCalc::AnchorEditor::Coefficients:
    {
      if( !coefficientsEdited() )
        return true;

      vector<float> coefs;
      coefs.reserve( m_coefEdits.size() );
      for( NativeFloatSpinBox * const edit : m_coefEdits )
        coefs.push_back( edit->value() );

      // What may be written back:
      //   - the user edited the table: their numbers, whatever they started from;
      //   - otherwise a real stored covariance is kept as-is, but ONLY while it still describes the
      //     equation.  Adding or removing a term resizes the shadow (which is why a value diff
      //     cannot tell "touched" from "resized"), and writing the resized matrix would claim the
      //     new coefficient is known exactly; leaving it to be dropped, with a note, is honest.
      //   - a manufactured assume-independent matrix the user never touched is not this DRF's
      //     covariance at all, so it is never written.
      const bool terms_changed = (m_seedState
                                  && (m_coefEdits.size() != m_seedState->coefficients.size()));
      const bool write_cov = m_coefCovTouched
                             || (!m_coefCovIsPlaceholder && !terms_changed);

      return DrfModifyCalc::applyCoefficients( working, coefs, m_coefSigmas, m_coefRho,
                                               equationEnergyUnits(), write_cov, problems );
    }//case Coefficients

    case DrfModifyCalc::AnchorEditor::Formula:
    {
      if( !formulaEdited() && !pointsEdited() )
        return true;

      vector<DrfModifyCalc::PointRow> rows;
      if( !collectPointRows( rows, problems ) )
        return false;

      return DrfModifyCalc::applyFormula( working, m_formulaText->text().toUTF8(),
                                          equationEnergyUnits(), rows, anchorOptions(), problems );
    }//case Formula

    case DrfModifyCalc::AnchorEditor::RefitPoints:
    case DrfModifyCalc::AnchorEditor::AbsolutePoints:
    case DrfModifyCalc::AnchorEditor::CurvePairs:
    {
      if( !pointsEdited() )
        return true;

      vector<DrfModifyCalc::PointRow> rows;
      if( !collectPointRows( rows, problems ) )
        return false;

      return DrfModifyCalc::applyPointRows( working, m_editor, rows, anchorOptions(),
                                            m_seedPoints, problems );
    }//case the three point editors
  }//switch( m_editor )

  return true;
}//applyAnchorTab(...)


void DrfModifyWidget::showProblems( const std::vector<DrfModifyCalc::Problem> &problems )
{
  for( const DrfModifyCalc::Problem &problem : problems )
  {
    WString message = WString::tr( problem.messageId );
    if( !problem.arg.empty() )
      message = message.arg( WString::fromUTF8(problem.arg) );

    passMessage( message, problem.blocking ? WarningWidget::WarningMsgHigh
                                           : WarningWidget::WarningMsgInfo );
  }//for( const Problem &problem : problems )
}//showProblems(...)


bool DrfModifyWidget::ToolState::operator==( const ToolState &rhs ) const
{
  if( (drfHash != rhs.drfHash)
     || (name != rhs.name) || (description != rhs.description) || (tabIndex != rhs.tabIndex)
     || (geometryModeled != rhs.geometryModeled)
     || (coefSigmas != rhs.coefSigmas)
     || (coefRho != rhs.coefRho)
     || (coefficients != rhs.coefficients)
     || (formula != rhs.formula)
     || (anchors != rhs.anchors)
     || (anchorRefDistance != rhs.anchorRefDistance)
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

  state->drfHash = m_orig ? m_orig->hashValue() : uint64_t(0);
  state->name = m_name->text().toUTF8();
  state->description = m_description->text().toUTF8();
  state->tabIndex = m_tabMenu->currentIndex();
  state->geometryModeled = m_geometryModeled;

  // The sigma/rho shadow is the authoritative numeric state (not the widgets).
  state->coefSigmas = m_coefSigmas;
  state->coefRho = m_coefRho;

  for( const NativeFloatSpinBox * const edit : m_coefEdits )
    state->coefficients.push_back( edit->text().toUTF8() );

  if( m_formulaText )
    state->formula = m_formulaText->text().toUTF8();

  for( const AnchorRow &r : m_anchors )
  {
    ToolState::RowState row;
    row.cells = { r.energy->text().toUTF8(), r.eff->text().toUTF8(),
                  r.stat->text().toUTF8(),
                  r.cert ? r.cert->text().toUTF8() : string(),
                  r.source ? r.source->text().toUTF8() : string(),
                  r.dist ? r.dist->text().toUTF8() : string() };
    row.seedIndex = r.seedIndex;
    state->anchors.push_back( row );
  }//for( const AnchorRow &r : m_anchors )

  if( m_anchorRefDistance )
    state->anchorRefDistance = m_anchorRefDistance->text().toUTF8();
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

  // An undo step is resolved against whatever Modify dialog is open when it runs, which may be a
  //  different detector than the one the step was recorded on (see ToolState::drfHash).
  const uint64_t this_drf = m_orig ? m_orig->hashValue() : uint64_t(0);
  if( state->drfHash != this_drf )
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

  // Coefficient boxes, then the sigma/rho shadow and the table it drives.  Bound both loops by the
  //  widget count: removeCoefficient() keeps at least one term, so a target of zero would otherwise
  //  never be reached.
  while( (m_coefEdits.size() > state->coefficients.size()) && (m_coefEdits.size() > 1) )
    removeCoefficient();
  while( m_coefEdits.size() < state->coefficients.size() )
    addCoefficient();
  const size_t ncoef_set = std::min( m_coefEdits.size(), state->coefficients.size() );
  for( size_t i = 0; i < ncoef_set; ++i )
    m_coefEdits[i]->setText( WString::fromUTF8(state->coefficients[i]) );

  m_coefSigmas = state->coefSigmas;
  m_coefRho = state->coefRho;
  rebuildCovTable();

  if( m_formulaText )
    m_formulaText->setText( WString::fromUTF8(state->formula) );
  if( m_effEnergyUnits )
    m_effEnergyUnits->setCurrentIndex( (state->efficiencyEnergyUnits > 10.0f) ? 1 : 0 );

  if( m_anchorTable )
  {
    while( !m_anchors.empty() )
      removeAnchorRow();
    for( const ToolState::RowState &a : state->anchors )
    {
      addAnchorRow( 0.0f, 0.0f, 0.0f, 0.0f, string(), -1.0f, a.seedIndex );
      m_anchors.back().energy->setText( WString::fromUTF8(a.cells[0]) );
      m_anchors.back().eff->setText( WString::fromUTF8(a.cells[1]) );
      m_anchors.back().stat->setText( WString::fromUTF8(a.cells[2]) );
      if( m_anchors.back().cert )
        m_anchors.back().cert->setText( WString::fromUTF8(a.cells[3]) );
      if( m_anchors.back().source )
        m_anchors.back().source->setText( WString::fromUTF8(a.cells[4]) );
      if( m_anchors.back().dist )
        m_anchors.back().dist->setText( WString::fromUTF8(a.cells[5]) );
    }//for( const ToolState::RowState &a : state->anchors )

    m_anchorRefDistance->setText( WString::fromUTF8(state->anchorRefDistance) );
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

  // Restoring is not a user edit, so nothing here should become the next undo step.  Nothing is done
  //  about response staleness either - it is re-derived from the restored content, so an undo/redo
  //  cannot strand it in the "fresh" state the way clearing a flag here used to.
  m_renderFlags.clear( RenderActions::AddUndoRedoStep );
  m_applyAfterGenerationId = -1;   //no run this snapshot describes is one we armed
  updateGenerateButton();
  refreshUncertSummary();
  markGeneralStale();
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

  updateGenerateButton();
  m_renderFlags |= RenderActions::RefreshSummary;
  markGeneralStale();
  scheduleUndoRedoStep();
}//void markEdited()


void DrfModifyWidget::markEditedNoRegen()
{
  if( m_restoringState )
    return;

  // For edits `DrfModifyCalc::seedFingerprint` does not cover - the name, the description, the FWHM
  //  - so they cannot make a response stale and there is no point re-deriving that (which rebuilds
  //  the working DRF).  The summary still has to be rebuilt (the name and the FWHM are both rows in
  //  it) and the edit is still an undo step.
  markGeneralStale();
  scheduleUndoRedoStep();
}//void markEditedNoRegen()


void DrfModifyWidget::markGeneralStale()
{
  m_generalStale = true;
  if( m_tabMenu && m_generalTabItem && (m_tabMenu->currentItem() == m_generalTabItem) )
    refreshGeneralTab();
}//void markGeneralStale()


void DrfModifyWidget::refreshGeneralTab()
{
  // Deliberately NOT called from render(): this replaces the info table's widgets, and doing that
  //  while Wt is rendering leaves it wiring event handlers to elements it has just replaced (which
  //  threw, and took the rest of the dialog's client-side setup with it).  The chart does not need
  //  deferring either - DrfChart re-sends whatever it holds when its client object is built.
  //
  // Nothing to do when nothing has changed - and it must NOT be done twice in one request, because
  //  the second `clear()` destroys widgets whose client-side wiring Wt has already queued.
  if( !m_generalStale )
    return;

  m_generalStale = false;
  if( !m_generalChart || !m_infoTable )
    return;

  // A preview must not nag about, e.g., an FWHM form the user has not filled in yet, so whatever
  //  could not be applied is collected and dropped - only `requestApply` shows any of it.
  vector<DrfModifyCalc::Problem> problems;
  const shared_ptr<DetectorPeakResponse> working = buildWorkingDrf( true, problems );
  m_generalChart->updateChart( working );   //draws nothing for a DRF with no efficiency yet
  fillInfoTable( working );
}//void refreshGeneralTab()


void DrfModifyWidget::buildInfoTable( Wt::WContainerWidget *parent )
{
  m_infoTable = parent->addNew<WTable>();
  m_infoTable->addStyleClass( "DrfModifyInfoTable" );

  const char *labelKeys[NumInfoRow] = {
    "dmw-info-eff", "dmw-info-geom", "dmw-info-support", "dmw-info-uncert",
    "dmw-info-fwhm", "dmw-info-total", "dmw-info-range", "dmw-info-diam"
  };

  for( int row = 0; row < NumInfoRow; ++row )
  {
    WText *label = m_infoTable->elementAt( row, 0 )->addNew<WText>( WString::tr( labelKeys[row] ) );
    label->addStyleClass( "DrfModifyInfoLabel" );

    //Plain: these values quote user text (a material name, a detector description).
    m_infoValues[row] = m_infoTable->elementAt( row, 1 )->addNew<WText>();
    m_infoValues[row]->setTextFormat( TextFormat::Plain );
  }//for( each row )
}//void buildInfoTable(...)


void DrfModifyWidget::fillInfoTable( const std::shared_ptr<const DetectorPeakResponse> &drf )
{
  if( !m_infoTable )
    return;

  // Only the text changes - see buildInfoTable.  A row with nothing to say is hidden.
  auto set_row = [this]( const InfoRow row, const WString &value ){
    if( !m_infoValues[row] )
      return;
    m_infoValues[row]->setText( value );
    if( m_infoTable->rowAt(row) )
      m_infoTable->rowAt(row)->setHidden( value.empty() );
  };

  auto length_str = []( const double cm ) -> WString {
    return WString::fromUTF8( PhysicalUnits::printToBestLengthUnits( cm * PhysicalUnits::cm, 3 ) );
  };

  const bool valid = (drf && drf->isValid());
  const shared_ptr<const ceelo::DetectorResponse> mc = drf ? drf->ceeloResponse() : nullptr;
  const shared_ptr<const MeasuredDrfPoints> points = drf ? drf->measuredPoints() : nullptr;
  const bool have_points = (points && !points->empty());

  // --- Efficiency: where the on-axis curve comes from -------------------------------------------
  {
    WString txt;
    if( !valid )
    {
      // The "generate one on the Geom & MC tab" wording only makes sense when that tab exists; a
      //  fixed-geometry DRF has no geometry to model and so no such tab.
      txt = WString::tr( m_geomTabItem ? "dmw-info-eff-none" : "dmw-info-eff-none-nogeom" );
    }else
    {
      switch( drf->efficiencyFcnType() )
      {
        case DetectorPeakResponse::kFunctialEfficienyForm:
          txt = WString::tr("dmw-info-eff-formula");
          break;

        case DetectorPeakResponse::kExpOfLogPowerSeries:
          txt = WString::tr("dmw-info-eff-exp-log")
                  .arg( static_cast<int>( drf->efficiencyExpOfLogsCoeffs().size() ) );
          break;

        case DetectorPeakResponse::kEnergyEfficiencyPairs:
          txt = WString::tr("dmw-info-eff-pairs")
                  .arg( static_cast<int>( drf->efficiencyCurve()->energyEfficiencies().size() ) );
          break;

        case DetectorPeakResponse::kNumEfficiencyFnctForms:
          txt = WString::tr( m_geomTabItem ? "dmw-info-eff-none" : "dmw-info-eff-none-nogeom" );
          break;
      }//switch( efficiency form )

      if( have_points )
        txt = WString::tr("dmw-info-eff-measured").arg( txt )
                .arg( static_cast<int>( points->points().size() ) );
    }//if( !valid ) / else

    set_row( InfoEfficiency, txt );
  }

  // --- Geometry: flat disk / fixed / physical shape ----------------------------------------------
  {
    WString txt;
    const shared_ptr<const ceelo::GeometryDescriptor> gd = drf ? drf->geometry() : nullptr;

    if( drf && drf->isFixedGeometry() )
    {
      const string &postfix = DetectorPeakResponse::det_eff_geom_type_postfix( drf->geometryType() );
      // `postfix` is an activity-UNIT suffix ("/cm2", "/m2", "/g"), so it has to be attached to a
      //  unit, not to the words "Fixed source geometry" - which read as "Fixed source geometry/cm2".
      txt = WString::tr("dmw-info-geom-fixed").arg( WString::fromUTF8(postfix) );
    }else if( gd )
    {
      string crystal;
      if( (gd->crystal_material_index >= 0)
          && (gd->crystal_material_index < static_cast<int>(gd->materials.size())) )
      {
        crystal = gd->materials[gd->crystal_material_index].name;
      }

      // Transverse extents are stored as halves, the length in full (CeeLo's convention).
      const vector<double> &dims = gd->dimensions_cm;
      if( (gd->shape == ceelo::DetectorShape::Box) && (dims.size() >= 3) )
        txt = WString::tr("dmw-info-geom-box").arg( WString::fromUTF8(crystal) )
                .arg( length_str(2.0*dims[0]) ).arg( length_str(2.0*dims[1]) ).arg( length_str(dims[2]) );
      else if( dims.size() >= 2 )
        txt = WString::tr("dmw-info-geom-cyl").arg( WString::fromUTF8(crystal) )
                .arg( length_str(2.0*dims[0]) ).arg( length_str(dims[1]) );
      else
        txt = WString::tr("dmw-info-geom-known");

      txt = WString::tr("dmw-info-geom-layers").arg( txt ).arg( static_cast<int>( gd->layers.size() ) );
    }else
    {
      const double diam_cm = drf ? (drf->detectorDiameter() / PhysicalUnits::cm) : 0.0;
      txt = WString::tr("dmw-info-geom-flat").arg( length_str(diam_cm) );
    }

    set_row( InfoGeometry, txt );
  }

  // --- Location support: how off-axis / near-field queries are answered ----------------------------
  {
    WString txt;
    if( mc )
    {
      const WString range = WString::tr("dmw-info-kev-range")
                              .arg( static_cast<int>( std::round(mc->provenance.valid_e_min_keV) ) )
                              .arg( static_cast<int>( std::round(mc->provenance.valid_e_max_keV) ) );
      switch( mc->provenance.method )
      {
        case ceelo::ProductionMethod::FullMc:
          txt = WString::tr("dmw-info-support-full")
                  .arg( WString::fromUTF8( ceelo::to_string(mc->provenance.profile) ) ).arg( range );
          break;

        case ceelo::ProductionMethod::QuickMcTransfer:
          txt = WString::tr("dmw-info-support-quick").arg( range );
          break;

        case ceelo::ProductionMethod::CurveTransfer:
          txt = WString::tr("dmw-info-support-curve").arg( range );
          break;
      }//switch( method )

      if( mc->model_transfer.has_value() )
      {
        // CeeLo quotes its floor from the crystal-face origin; the user is told a face distance.
        const string min_dist = PhysicalUnits::printToBestLengthUnits(
                CeeLoUtils::faceDistanceFromCrystalOrigin( mc->descriptor, mc->provenance.min_distance_cm )
                * PhysicalUnits::cm );
        txt = WString::tr("dmw-info-support-min-dist").arg( txt ).arg( WString::fromUTF8(min_dist) );
      }//if( a transfer response )

      if( !mc->grounding.empty() )
        txt = WString::tr( mc->grounding.curve_derived ? "dmw-info-grounded-curve"
                                                       : "dmw-info-grounded-points" ).arg( txt );
    }else
    {
      txt = WString::tr("dmw-info-support-none");
    }//if( mc ) / else

    set_row( InfoSupport, txt );
  }

  // --- Uncertainty: which of the (mutually exclusive) sources applies ----------------------------
  {
    WString txt;
    const shared_ptr<const DetectorEfficiencyUncert> uncert = drf ? drf->efficiencyUncert() : nullptr;

    if( mc )
    {
      txt = WString::tr( mc->grounding.empty() ? "dmw-info-uncert-mc"
                                               : "dmw-info-uncert-mc-grounded" );
    }else if( have_points )
    {
      std::set<string> sources;
      for( const MeasuredEffPoint &p : points->points() )
        sources.insert( p.sourceKey );
      txt = WString::tr("dmw-info-uncert-points")
              .arg( static_cast<int>( points->points().size() ) ).arg( static_cast<int>( sources.size() ) );
    }else if( uncert && uncert->hasNodeCovariance() )
    {
      txt = WString::tr("dmw-info-uncert-nodes").arg( static_cast<int>( uncert->covarianceEnergies().size() ) );
    }else if( uncert && !uncert->isEmpty() )
    {
      txt = WString::tr("dmw-info-uncert-coef");
    }else
    {
      txt = WString::tr("dmw-info-uncert-none");
    }

    set_row( InfoUncert, txt );
  }

  // --- FWHM -----------------------------------------------------------------------------------------
  {
    WString txt;
    if( drf && drf->hasResolutionInfo() )
    {
      const char *formKey = "dmw-info-fwhm-gadras";
      switch( drf->resolutionFcnType() )
      {
        case DetectorPeakResponse::kGadrasResolutionFcn:    formKey = "dmw-info-fwhm-gadras";     break;
        case DetectorPeakResponse::kSqrtPolynomial:         formKey = "dmw-info-fwhm-sqrt-poly";  break;
        case DetectorPeakResponse::kSqrtEnergyPlusInverse:  formKey = "dmw-info-fwhm-sqrt-inv";   break;
        case DetectorPeakResponse::kConstantPlusSqrtEnergy: formKey = "dmw-info-fwhm-const-sqrt"; break;
        case DetectorPeakResponse::kNumResolutionFnctForm:  break;
      }//switch( resolution form )

      txt = WString::tr("dmw-info-fwhm-form").arg( WString::tr(formKey) )
              .arg( static_cast<int>( drf->resolutionFcnCoefficients().size() ) );

      const float fwhm662 = drf->peakResolutionFWHM( 661.7f );
      if( fwhm662 > 0.0f )
        txt = WString::tr("dmw-info-fwhm-at-662").arg( txt )
                .arg( WString::fromUTF8( SpecUtils::printCompact( fwhm662, 3 ) ) );
    }else
    {
      txt = WString::tr("dmw-info-fwhm-none");
    }

    set_row( InfoFwhm, txt );
  }

  // --- Total efficiency (any energy deposited, not just the full peak) ---------------------------
  //  A generated Monte-Carlo response always carries a total-efficiency model - the same photon
  //  histories score both - so ask #hasAnyTotalEfficiencyInfo, not #hasTotalEfficiency (which only
  //  knows about an explicitly-set curve and answers "No" for every MC-backed detector).
  {
    WString txt;
    if( drf && drf->hasTotalEfficiency() )
      txt = WString::tr("dmw-info-toteff-curve");
    else if( drf && drf->hasAnyTotalEfficiencyInfo() )
      txt = WString::tr("dmw-info-toteff-mc");
    else
      txt = WString::tr("dmw-info-no");

    // The peak-to-total ratio is what decides whether cascade summing matters, so quote it at the
    //  usual reference energy when that energy is one this detector actually covers.
    const float ref_energy = 661.7f;
    const bool in_range = valid && (ref_energy >= drf->lowerEnergy())
                          && ((drf->upperEnergy() <= drf->lowerEnergy()) || (ref_energy <= drf->upperEnergy()));
    if( in_range && drf->hasAnyTotalEfficiencyInfo() )
    {
      try
      {
        const float fep = drf->farFieldIntrinsicEfficiency( ref_energy );
        const float tot = drf->totalIntrinsicEfficiencyAny( ref_energy );
        if( (fep > 0.0f) && (tot > fep) )
          txt = WString::tr("dmw-info-toteff-ratio").arg( txt )
                  .arg( WString::fromUTF8( SpecUtils::printCompact( tot/fep, 2 ) ) );
      }catch( std::exception & )
      {
        //An efficiency formula that will not evaluate here; the source text alone is enough.
      }
    }//if( worth quoting a ratio )

    set_row( InfoTotalEff, txt );
  }

  if( valid && (drf->upperEnergy() > (drf->lowerEnergy() + 1.0)) )
    set_row( InfoRange, WString::tr("dmw-info-kev-range")
                          .arg( static_cast<int>( std::round(drf->lowerEnergy()) ) )
                          .arg( static_cast<int>( std::round(drf->upperEnergy()) ) ) );
  else
    set_row( InfoRange, WString() );

  if( drf && (drf->detectorDiameter() > 0.0f) )
  {
    WString txt = length_str( drf->detectorDiameter() / PhysicalUnits::cm );
    if( drf->detectorSetback() > 0.0 )
      txt = WString::tr("dmw-info-diam-setback").arg( txt )
              .arg( length_str( drf->detectorSetback() / PhysicalUnits::cm ) );
    set_row( InfoDiameter, txt );
  }else
  {
    set_row( InfoDiameter, WString() );
  }
}//void fillInfoTable(...)


void DrfModifyWidget::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool add_step = m_renderFlags.test( RenderActions::AddUndoRedoStep );
  const bool rebuild_cov = m_renderFlags.test( RenderActions::RebuildCovTable );
  const bool refresh_summary = m_renderFlags.test( RenderActions::RefreshSummary );
  m_renderFlags = Wt::WFlags<RenderActions>();

  // The covariance cell edits defer their table rebuild here, so the WLineEdit whose `changed()`
  //  fired is not deleted from inside its own event handler.
  if( rebuild_cov )
    rebuildCovTable();

  if( refresh_summary )
    refreshUncertSummary();

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


std::shared_ptr<DetectorPeakResponse> DrfModifyWidget::buildWorkingDrf( const bool includeMcResponse,
                                            std::vector<DrfModifyCalc::Problem> &problems )
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

  // Anchor tab: apply whichever editor this DRF gets - which follows how it represents its
  //  efficiency, not the mode toggle, so there is no longer an editor the user could have typed
  //  into that the apply skips.  Done before the MC step so a grounding sees the edited points.
  applyAnchorTab( *working, problems );

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
        problems.push_back( DrfModifyCalc::Problem( "dmw-err-no-fwhm-fit", string(), false ) );
      }
    }catch( std::exception &e )
    {
      problems.push_back( DrfModifyCalc::Problem( "dmw-err-fwhm-not-applied", e.what(), false ) );
    }
  }//if( m_fwhmTool )

  // The geometry the user described.  Recorded whether or not a response is attached: the response
  //  carries its own descriptor, but a later switch to Flat Disk must not leave the detector with no
  //  statement of what it physically is (and hence unable to ever have a response again).
  if( m_mcTool && m_mcTool->geometryInput() && m_mcTool->generationReady() )
  {
    try
    {
      working->setGeometry( make_shared<const ceelo::GeometryDescriptor>(
                                                  m_mcTool->geometryInput()->toDescriptor() ) );
    }catch( std::exception & )
    {
      //an incomplete form; whatever the DRF already stated stands
    }
  }//if( the geometry form holds a real geometry )

  // Monte-Carlo response.  In Geometry-Modeled mode, attach the generated response (keeping any the
  //  DRF already carried when nothing new was generated); in Flat Disk, or as a regeneration seed,
  //  detach it - CRITICAL for a regeneration seed, since an attached CeeLo response overrides the
  //  manual points/covariance at query time (the very thing the regeneration is grounding on).
  const shared_ptr<const ceelo::DetectorResponse> resp
      = (includeMcResponse && m_geometryModeled && m_mcTool) ? m_mcTool->generatedResponse() : nullptr;

  // What shape this detector knows, taken BEFORE any detach below.  `geometry()` prefers an attached
  //  response's own descriptor, and the serializers write only one of the two - so for a DRF that
  //  has been through a file or the database the shape lives *only* in the response, and detaching
  //  it would erase the crystal outright.  Re-applied after the detach so Flat Disk keeps the
  //  geometry, as this function has always claimed to.
  //  Copied rather than aliased: `geometry()` hands back a pointer that shares ownership with the
  //  response, which would keep the whole (~100 KB) response alive behind a detached DRF.
  const shared_ptr<const ceelo::GeometryDescriptor> from_drf = working->geometry();
  const shared_ptr<const ceelo::GeometryDescriptor> known_geom
      = from_drf ? make_shared<const ceelo::GeometryDescriptor>( *from_drf ) : nullptr;
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
        problems.push_back( DrfModifyCalc::Problem( "dmw-err-no-backbone", e.what() ) );
      }
    }//if( !working->isValid() )

    working->setCeeloResponse( resp );
  }else if( includeMcResponse && !m_geometryModeled )
  {
    working->setCeeloResponse( nullptr );  //Flat Disk: detach any geometry-modeled response
    working->setGeometry( known_geom );    //but the detector is still the shape it always was
  }else if( !includeMcResponse )
  {
    working->setCeeloResponse( nullptr );  //regeneration seed: manual points/covariance must drive
    working->setGeometry( known_geom );    //the seed is re-characterized from this same shape
  }
  //else: Geometry Modeled with nothing newly generated - keep whatever response the DRF came with.
  //  Whether that response still describes these edits is #responseStale's job, not this one's.

  // Geometry Modeled without any response (none generated, none carried): keep what was typed into
  //  the form - the geometry is a fact about the detector, and is what a later generation, the
  //  General tab, and an export read.  Only a complete, user-confirmed geometry, never the
  //  length-equals-diameter guess.  (With a response attached, that carries its own geometry, and
  //  setGeometry is ignored anyway.)
  if( includeMcResponse && m_geometryModeled && m_mcTool && !working->ceeloResponse()
      && m_mcTool->generationReady() )
  {
    try
    {
      working->setGeometry( make_shared<const ceelo::GeometryDescriptor>( m_mcTool->geometryDescriptor() ) );
    }catch( std::exception & )
    {
      //generationReady() implies the form converts; nothing sensible to do if it somehow does not.
    }
  }//if( Geometry Modeled with no response )

  return working;
}//buildWorkingDrf(...)


void DrfModifyWidget::apply()
{
  vector<DrfModifyCalc::Problem> problems;
  shared_ptr<DetectorPeakResponse> working = buildWorkingDrf( true, problems );

  // requestApply() has already validated, so a blocking problem here would be a logic error - but
  // refuse rather than emit a DRF that does not hold what the user typed.
  if( DrfModifyCalc::anyBlocking( problems ) )
  {
    showProblems( problems );
    return;
  }

#if( PERFORM_DEVELOPER_CHECKS )
  // The editing invariants (see DrfModifyCalc), on the DRF actually being handed over - and only
  //  when this dialog is what made it inconsistent.  A DRF that arrived that way (a legacy stored
  //  one), or one whose equation the user deliberately hand-edited away from its points (which is a
  //  supported edit, reported as `dmw-note-points-stale`), is not a programming error.
  {
    bool points_now_stale = false;
    for( const DrfModifyCalc::Problem &problem : problems )
      points_now_stale = (points_now_stale || (problem.messageId == "dmw-note-points-stale"));

    std::string orig_why, why;
    const bool orig_ok = m_orig ? DrfModifyCalc::checkDrfSelfConsistent( *m_orig, orig_why ) : true;

    if( orig_ok && !points_now_stale && working->isValid()
        && !DrfModifyCalc::checkDrfSelfConsistent( *working, why ) )
    {
      log_developer_error( __func__,
        ("Modify-DRF produced an inconsistent DRF: " + why).c_str() );
      assert( 0 );
    }
  }
#endif

  showProblems( problems );   //the non-blocking notes: what else the apply did
  m_updatedDrf.emit( working );
}//apply()


void DrfModifyWidget::requestApply()
{
  // Nothing the user typed may be silently dropped: build once, up front, and refuse with the reason
  //  rather than closing the dialog looking successful.
  {
    vector<DrfModifyCalc::Problem> problems;
    buildWorkingDrf( true, problems );
    if( DrfModifyCalc::anyBlocking( problems ) )
    {
      showProblems( problems );
      return;
    }
  }

  // Flat Disk, but the original carried a geometry-modeled response: confirm the detach first.
  // Either the response the detector arrived with, or one generated in this session - a Monte Carlo
  //  the user just waited minutes for is exactly the thing not to discard without asking.
  const bool losing_response = !m_geometryModeled
                               && ((m_orig && m_orig->ceeloResponse())
                                   || (m_mcTool && m_mcTool->generatedResponse()));
  if( losing_response )
  {
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dmw-detach-title"),
                                                             WString::tr("dmw-detach-body") );
    WPushButton *ok = dialog->addButton( WString::tr("dmw-detach-accept") );
    dialog->addButton( WString::tr("Cancel") );
    ok->clicked().connect( this, &DrfModifyWidget::apply );
    return;
  }//if( detaching a geometry-modeled response )

  const bool canGen = (m_mcTool && m_mcTool->generationReady());

  // Geometry Modeled with no response at all - none generated in this session, none carried by the
  // original.  Using it silently would just give a flat-disk detector that happens to know its
  // shape, so offer to generate first, or say what blocks generating.
  const bool haveResp = (m_mcTool && m_mcTool->generatedResponse())
                        || (m_orig && m_orig->ceeloResponse());
  if( m_geometryModeled && m_mcTool && !haveResp )
  {
    if( canGen )
    {
      SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dmw-nogen-title"),
                                                               WString::tr("dmw-nogen-body") );
      WPushButton *gen = dialog->addButton( WString::tr("dmw-nogen-generate") );
      WPushButton *useAnyway = dialog->addButton( WString::tr("dmw-nogen-use-anyway") );
      dialog->addButton( WString::tr("Cancel") );
      gen->clicked().connect( this, [this](){
        //handleResponseGenerated applies once THIS run lands; a run that never starts, or that
        //  fails part way, must not leave a later unrelated generation armed.
        const bool started = handleGenerateResponse();
        m_applyAfterGenerationId = started ? m_mcTool->generationId() : -1;

        // Otherwise this dialog just closes and nothing happens: the reason lands on the Geom & MC
        //  tab's status line, which is not where the user is looking after pressing a footer button.
        if( !started )
          passMessage( WString::tr("dmw-nogen-not-started")
                         .arg( WString::fromUTF8( m_mcTool->geometryProblem() ) ),
                       WarningWidget::WarningMsgHigh );
      } );
      // No response is attached, so nothing is answering the queries in place of the edits - unlike
      //  the stale case below, using it as-is is a coherent (flat-disk) detector.
      useAnyway->clicked().connect( this, &DrfModifyWidget::apply );
    }else
    {
      // The reason can quote a typed material name, so encode it for the dialog's XHTML text.
      const string problem = Wt::Utils::htmlEncode( m_mcTool->geometryProblem() );
      SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dmw-geom-incomplete-title"),
                             WString::tr("dmw-geom-incomplete-body").arg( WString::fromUTF8(problem) ) );
      WPushButton *useAnyway = dialog->addButton( WString::tr("dmw-geom-incomplete-use-anyway") );
      dialog->addButton( WString::tr("Cancel") );
      useAnyway->clicked().connect( this, &DrfModifyWidget::apply );
    }//if( canGen ) / else

    return;
  }//if( Geometry Modeled with no response )

  // A response that IS attached but does not reflect the edits is not a cosmetic mismatch: while it
  //  is attached it answers every efficiency and uncertainty query, so the edit would simply be
  //  ignored.  There is therefore no "use anyway" here - either the response is rebuilt, or it goes.
  if( m_geometryModeled && responseStale() )
  {
    const bool instant = (m_mcTool
             && (m_mcTool->selectedMethod() == MakeMcResponseForDrf::Method::CurveTransfer));

    // The curve-transfer method is deterministic and instant, so there is nothing to ask about.
    if( canGen && instant )
    {
      // Keyed to the generation it arms: a run that never starts, or that fails part way, must not
      //  leave some later, unrelated generation applying and closing this dialog.
      const bool started = handleGenerateResponse();
      m_applyAfterGenerationId = started ? m_mcTool->generationId() : -1;
      if( !started )
        passMessage( WString::tr("dmw-err-regen-failed"), WarningWidget::WarningMsgHigh );
      return;
    }

    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dmw-regen-title"),
                   WString::tr( canGen ? "dmw-regen-required-body" : "dmw-regen-impossible-body" ) );
    if( canGen )
    {
      WPushButton *regen = dialog->addButton( WString::tr("dmw-regen-accept") );
      regen->clicked().connect( this, [this](){
        const bool started = handleGenerateResponse();
        m_applyAfterGenerationId = started ? m_mcTool->generationId() : -1;
        if( !started )
          passMessage( WString::tr("dmw-err-regen-failed"), WarningWidget::WarningMsgHigh );
      } );
    }
    WPushButton *detach = dialog->addButton( WString::tr("dmw-regen-detach") );
    detach->clicked().connect( this, &DrfModifyWidget::detachResponseAndApply );
    dialog->addButton( WString::tr("Cancel") );
    return;
  }//if( the attached response no longer describes the edits )

  apply();
}//requestApply()


void DrfModifyWidget::detachResponseAndApply()
{
  m_geometryModeled = false;
  if( m_modeToggle )
    m_modeToggle->setChecked( false );
  if( m_mcTool )
    m_mcTool->setDisabled( true );

  updateGenerateButton();
  apply();
}//detachResponseAndApply()


void DrfModifyWidget::handleModeToggle()
{
  m_geometryModeled = (m_modeToggle && m_modeToggle->isChecked());
  if( m_mcTool )
    m_mcTool->setDisabled( !m_geometryModeled );  //Flat Disk greys the whole tool

  // Deliberately does not touch which Anchor editor is showing: that follows how the efficiency is
  //  represented, and sending the apply down a different editor because of this toggle is what used
  //  to discard whichever editor's edits were not applied.
  //  Nor does it make a current response stale: the toggle decides WHETHER a response is attached,
  //  not what one built from these inputs would contain, and `responseStale()` fingerprints only
  //  the latter.  (markEdited refreshes the generate button, whose visibility follows the mode.)
  updateAnchorEditorVisibility();
  markEdited();
}//handleModeToggle()


bool DrfModifyWidget::geometryModeled() const
{
  return m_geometryModeled;
}


void DrfModifyWidget::updateAnchorEditorVisibility()
{
  const bool show_points = DrfModifyCalc::editorUsesPointTable( m_editor );
  const bool show_coefs = (m_editor == DrfModifyCalc::AnchorEditor::Coefficients);
  const bool show_formula = (m_editor == DrfModifyCalc::AnchorEditor::Formula);

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
    switch( m_editor )
    {
      case DrfModifyCalc::AnchorEditor::RefitPoints:    help_id = "dmw-anchor-help-refit";    break;
      case DrfModifyCalc::AnchorEditor::AbsolutePoints: help_id = "dmw-anchor-help-points";   break;
      case DrfModifyCalc::AnchorEditor::CurvePairs:     help_id = "dmw-anchor-help-points-nosrc"; break;
      case DrfModifyCalc::AnchorEditor::Coefficients:   help_id = "dmw-anchor-help-coefs";    break;
      case DrfModifyCalc::AnchorEditor::Formula:        help_id = "dmw-anchor-help-formula";  break;
    }//switch( m_editor )

    m_anchorHelp->setText( WString::tr( help_id ) );
  }//if( m_anchorHelp )

  // Whether a response is (going to be) in charge of every query.
  if( m_responseNote )
  {
    const bool have_resp = (m_geometryModeled
                            && ((m_mcTool && m_mcTool->generatedResponse())
                                || (m_orig && m_orig->ceeloResponse())));
    m_responseNote->setHidden( !have_resp );
  }
}//updateAnchorEditorVisibility()


void DrfModifyWidget::seedCoefCovFromUncert( const std::shared_ptr<const DetectorEfficiencyUncert> &uncert )
{
  m_coefSigmas.clear();
  m_coefRho.clear();
  m_coefCovIsPlaceholder = false;

  const size_t n = m_coefEdits.size();
  if( uncert && n )
  {
    const vector<float> &cov = uncert->coefficientCovariance();
    if( cov.size() == (n*n) )
      DrfModifyCalc::sigmaRhoFromCovariance( cov, m_coefSigmas, m_coefRho );
  }//if( uncert && n )

  // No stored covariance: start from the legacy per-coefficient uncertainties when there are any, so
  //  the user has something to edit rather than a grid of zeros.  It is only a starting point - it
  //  assumes the coefficients are independent, which for a log-power-series fit they emphatically are
  //  not - so it is labelled as such, and applyAnchorTab only writes it if the user edits it.
  if( m_coefSigmas.empty() && n && m_orig )
  {
    if( DrfModifyCalc::sigmaRhoFromLegacyUncerts( *m_orig, m_coefSigmas, m_coefRho ) )
      m_coefCovIsPlaceholder = true;
  }

  if( m_coefSigmas.size() != n )
  {
    m_coefSigmas.assign( n, 0.0 );
    m_coefRho.assign( n*n, 0.0 );
    for( size_t i = 0; i < n; ++i )
      m_coefRho[i*n + i] = 1.0;
  }

  if( m_covPlaceholderNote )
    m_covPlaceholderNote->setHidden( !m_coefCovIsPlaceholder );

  rebuildCovTable();
}//seedCoefCovFromUncert(...)


void DrfModifyWidget::rebuildCovTable()
{
  if( !m_covTable )
    return;

  m_covTable->clear();  //wipes all cells and the WLineEdits/WTexts they hold

  const size_t n = m_coefEdits.size();
  if( (m_coefSigmas.size() != n) || (m_coefRho.size() != (n*n)) )
  {
    m_coefSigmas.assign( n, 0.0 );
    m_coefRho.assign( n*n, 0.0 );
    for( size_t i = 0; i < n; ++i )
      m_coefRho[i*n + i] = 1.0;
  }

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

    for( size_t j = 0; j < n; ++j )
    {
      const int tcol = static_cast<int>(j) + 1;

      if( j == i )
      {
        // Diagonal: the 1-sigma uncertainty of the coefficient itself (absolute, not a percent -
        //  these coefficients are logs, so a percentage of one means nothing).
        char buf[32];
        snprintf( buf, sizeof(buf), "%.4g", m_coefSigmas[i] );
        WLineEdit *e = m_covTable->elementAt(trow,tcol)->addNew<WLineEdit>();
        e->setTextSize( 8 );
        e->addStyleClass( "DrfCovDiag" );
        e->setText( WString::fromUTF8(buf) );
        e->changed().connect( this, [this,i,e](){ covSigmaChanged( i, e->text().toUTF8() ); } );
      }else
      {
        char buf[32];
        snprintf( buf, sizeof(buf), "%.3g", m_coefRho[i*n + j] );
        WLineEdit *e = m_covTable->elementAt(trow,tcol)->addNew<WLineEdit>();
        e->setTextSize( 6 );
        e->setText( WString::fromUTF8(buf) );

        if( j > i )
        {
          // Upper triangle: correlation (editable).  Editable even when a sigma is zero - the
          //  correlation is kept in its own matrix, so it survives a sigma being zeroed and typed
          //  back, which a covariance-only shadow could not.
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

  // An impossible set of sigmas and correlations does not fail loudly downstream: the fit's Cholesky
  //  whitening quietly gives up and falls back to counting statistics only, so the user would get
  //  LESS uncertainty for entering something that cannot happen.  Say so here, and refuse on apply.
  if( m_covWarning )
  {
    const vector<double> cov = DrfModifyCalc::covarianceFromSigmaRho( m_coefSigmas, m_coefRho );
    std::string why;
    const bool usable = DetectorEfficiencyUncert::covarianceIsUsable( cov, &why );
    m_covWarning->setText( WString::tr("dmw-err-cov-not-psd") );
    m_covWarning->setHidden( usable );
    if( m_covTable )
      m_covTable->toggleStyleClass( "Wt-invalid", !usable );
  }//if( m_covWarning )
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

  // Grow the shadow to match, leaving the existing block bit-exact.
  const size_t nn = m_coefEdits.size();
  if( (m_coefSigmas.size() != n) || (m_coefRho.size() != (n*n)) )
  {
    m_coefSigmas.assign( n, 0.0 );
    m_coefRho.assign( n*n, 0.0 );
    for( size_t i = 0; i < n; ++i )
      m_coefRho[i*n + i] = 1.0;
  }

  m_coefSigmas.push_back( 0.0 );
  vector<double> grown( nn*nn, 0.0 );
  for( size_t i = 0; i < nn; ++i )
    grown[i*nn + i] = 1.0;
  for( size_t i = 0; i < n; ++i )
    for( size_t j = 0; j < n; ++j )
      grown[i*nn + j] = m_coefRho[i*n + j];
  m_coefRho = grown;

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
  if( m_coefSigmas.size() == n )
    m_coefSigmas.pop_back();
  if( m_coefRho.size() == (n*n) )
  {
    vector<double> shrunk( nn*nn, 0.0 );
    for( size_t i = 0; i < nn; ++i )
      for( size_t j = 0; j < nn; ++j )
        shrunk[i*nn + j] = m_coefRho[i*n + j];
    m_coefRho = shrunk;
  }

  if( m_removeCoef )
    m_removeCoef->setEnabled( nn > 1 );

  m_renderFlags |= RenderActions::RebuildCovTable;
  scheduleRender();
}//removeCoefficient()


void DrfModifyWidget::covSigmaChanged( const std::size_t i, const std::string &text )
{
  const size_t n = m_coefEdits.size();
  if( (i >= n) || (m_coefSigmas.size() != n) )
    return;

  m_coefCovTouched = true;

  double s_new = 0.0;
  try{ s_new = std::stod( text ); }
  catch( std::exception & ){ m_renderFlags |= RenderActions::RebuildCovTable; scheduleRender(); return; }

  // Only this sigma changes; the correlations live in their own matrix, so typing 0 here (and then
  //  typing the value back) cannot destroy them - which it did when the shadow was a covariance.
  m_coefSigmas[i] = std::max( 0.0, s_new );

  m_renderFlags |= RenderActions::RebuildCovTable;
  markEdited();
}//covSigmaChanged(...)


void DrfModifyWidget::covRhoChanged( const std::size_t i, const std::size_t j, const std::string &text )
{
  const size_t n = m_coefEdits.size();
  if( (i >= n) || (j >= n) || (i == j) || (m_coefRho.size() != (n*n)) )
    return;

  m_coefCovTouched = true;

  double rho = 0.0;
  try{ rho = std::stod( text ); }
  catch( std::exception & ){ m_renderFlags |= RenderActions::RebuildCovTable; scheduleRender(); return; }
  if( rho < -1.0 ) rho = -1.0;
  if( rho >  1.0 ) rho =  1.0;

  m_coefRho[i*n + j] = rho;
  m_coefRho[j*n + i] = rho;

  m_renderFlags |= RenderActions::RebuildCovTable;
  markEdited();
}//covRhoChanged(...)


bool DrfModifyWidget::pointsEdited() const
{
  if( !m_seedState )
    return false;

  const shared_ptr<ToolState> now = currentState();

  return (now->anchors != m_seedState->anchors)
         || (now->anchorRefDistance != m_seedState->anchorRefDistance)
         || (now->anchorCorrLength != m_seedState->anchorCorrLength);
}//pointsEdited()


bool DrfModifyWidget::coefficientsEdited() const
{
  if( !m_seedState )
    return false;

  const shared_ptr<ToolState> now = currentState();

  return (now->coefficients != m_seedState->coefficients)
         || (now->coefSigmas != m_seedState->coefSigmas)
         || (now->coefRho != m_seedState->coefRho)
         || (now->efficiencyEnergyUnits != m_seedState->efficiencyEnergyUnits);
}//coefficientsEdited()


bool DrfModifyWidget::coefCovarianceEdited() const
{
  if( !m_seedState )
    return false;

  const shared_ptr<ToolState> now = currentState();

  return (now->coefSigmas != m_seedState->coefSigmas) || (now->coefRho != m_seedState->coefRho);
}//coefCovarianceEdited()


bool DrfModifyWidget::formulaEdited() const
{
  if( !m_seedState )
    return false;

  const shared_ptr<ToolState> now = currentState();

  return (now->formula != m_seedState->formula)
         || (now->efficiencyEnergyUnits != m_seedState->efficiencyEnergyUnits);
}//formulaEdited()


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


bool DrfModifyWidget::responseStale()
{
  if( !m_geometryModeled || !m_mcTool )
    return false;

  // Nothing to be stale unless a response would actually be attached and answer the queries.
  const shared_ptr<const ceelo::DetectorResponse> resp = m_mcTool->generatedResponse()
                          ? m_mcTool->generatedResponse()
                          : (m_orig ? m_orig->ceeloResponse() : nullptr);
  if( !resp )
    return false;

  // A method change makes the attached response the wrong KIND of response, whatever it was built
  //  from.  (The two enumerations mirror each other; see MakeMcResponseForDrf's use of the same
  //  mapping when it opens on an existing response.)
  if( static_cast<int>(m_mcTool->selectedMethod()) != static_cast<int>(resp->provenance.method) )
    return true;

  vector<DrfModifyCalc::Problem> problems;
  const shared_ptr<DetectorPeakResponse> seed = buildWorkingDrf( false, problems );
  if( !seed )
    return false;

  return (DrfModifyCalc::seedFingerprint( *seed ) != m_generatedFromFingerprint);
}//responseStale()


bool DrfModifyWidget::handleGenerateResponse()
{
  if( !m_mcTool || !m_geometryModeled )
    return false;

  // The seed comes from the provider installed in the constructor (the live edits, with any response
  //  detached), so the manual points/covariance ground the regeneration rather than being overridden
  //  by the response we are replacing.
  //
  // For a raw-measured-points DRF, the regeneration derives its uncertainty from those points, not
  //  from a hand-edited matrix - say so once, so an edited covariance is not silently discarded.
  if( m_editor == DrfModifyCalc::AnchorEditor::RefitPoints )
    passMessage( WString::tr("dmw-regen-from-points-note"), WarningWidget::WarningMsgInfo );

  const bool started = m_mcTool->startGeneration();
  updateGenerateButton();

  return started;
}//handleGenerateResponse()


void DrfModifyWidget::updateGenerateButton()
{
  if( !m_generateBtn )
    return;

  const bool haveResp = (m_mcTool && m_mcTool->generatedResponse());
  const bool canGen = (m_mcTool && m_mcTool->generationReady());
  // `m_result` is cleared for the whole duration of a run, so without this the button re-enables
  //  itself the moment generation starts and a second click queues a second full-core Monte Carlo.
  const bool running = (m_mcTool && m_mcTool->isGenerating());

  // Shown whenever the mode is Geometry Modeled - a button that only appears once the geometry is
  //  complete leaves the user hunting for it.  Enabled when the geometry is complete enough to run
  //  and there is work to do: the first generation (no response yet), or a regeneration once the
  //  response no longer reflects the edits.  While the geometry blocks it, the reason is spelled
  //  out beside it, in the export tip's place.
  m_generateBtn->setHidden( !m_geometryModeled );
  m_generateBtn->setEnabled( m_geometryModeled && canGen && !running
                             && (!haveResp || responseStale()) );

  const string problem = (m_geometryModeled && !canGen) ? m_mcTool->geometryProblem() : string();
  if( m_generateHint )
  {
    m_generateHint->setText( WString::fromUTF8( problem ) );
    m_generateHint->setHidden( problem.empty() );
  }
  if( m_exportNote )
    m_exportNote->setHidden( !problem.empty() );
}//updateGenerateButton()


void DrfModifyWidget::handleResponseGenerated( std::shared_ptr<ceelo::DetectorResponse> response )
{
  // What the response that just landed was built from - the seed the provider handed the generation.
  //  Staleness is this comparison, rather than a flag some path could forget to clear.
  if( response )
    m_generatedFromFingerprint = m_pendingSeedFingerprint;

  updateGenerateButton();
  updateAnchorEditorVisibility();  //a response is now in charge of the queries
  refreshUncertSummary();
  markGeneralStale();

  // Only for the run "Generate & use" armed: a different run landing here is one the user started
  //  themselves (or the transfer method rebuilding itself), and must not apply and close the dialog.
  const bool armed = (m_applyAfterGenerationId >= 0) && m_mcTool
                     && (m_applyAfterGenerationId == m_mcTool->generationId());
  m_applyAfterGenerationId = -1;
  if( armed && response )
    apply();   //the regenerate-then-use flow: the response now reflects the edits
}//handleResponseGenerated(...)


void DrfModifyWidget::refreshUncertSummary()
{
  if( !m_uncertSummary )
    return;

  // Ask the DRF the apply would produce, through the same call the activity/shielding fit uses, so
  //  this line cannot drift from what the analysis propagates (and from the chart band).
  vector<DrfModifyCalc::Problem> problems;
  shared_ptr<DetectorPeakResponse> working;
  try
  {
    working = buildWorkingDrf( true, problems );
  }catch( std::exception & )
  {
    working = nullptr;
  }

  const DrfModifyCalc::UncertSummary summary = working
                      ? DrfModifyCalc::uncertSummary( *working )
                      : DrfModifyCalc::UncertSummary{};

  if( !summary.valid || (summary.total <= 0.0) )
  {
    m_uncertSummary->setText( WString::tr("dmw-uncert-summary-none") );
    m_uncertSummary->show();
    updateCorrelationNote();
    return;
  }

  char energy[32], total[32], data[32], model[32];
  snprintf( energy, sizeof(energy), "%.1f", summary.energy );
  snprintf( total, sizeof(total), "%.2g", 100.0*summary.total );
  snprintf( data, sizeof(data), "%.2g", 100.0*summary.data );
  snprintf( model, sizeof(model), "%.2g", 100.0*summary.model );

  // Three wordings, because "2.1% from this detector's data" is a lie for a detector that never
  //  stated an uncertainty - it is a default standing in for one (see UncertSummary::dataIsAssumed).
  const char *id = "dmw-uncert-summary";
  if( summary.dataIsAssumed )
    id = "dmw-uncert-summary-assumed";
  else if( summary.model > 0.0 )
    id = "dmw-uncert-summary-split";

  WString text = WString::tr( id ).arg( energy ).arg( total );
  if( summary.model > 0.0 || summary.dataIsAssumed )
    text = text.arg( data ).arg( model );

  m_uncertSummary->setText( text );
  m_uncertSummary->show();

  updateCorrelationNote();
}//refreshUncertSummary()


void DrfModifyWidget::updateCorrelationNote()
{
  if( !m_corrInertNote )
    return;

  // The correlation control scales the "Corr. %" column and nothing else, so with that column empty
  //  it has no effect on any number this detector reports - including the one on the line above,
  //  which a user who has just set "Fully correlated" will otherwise read as disagreeing with it.
  bool any_correlated = false;
  for( const AnchorRow &r : m_anchors )
  {
    if( !r.cert )
      continue;

    bool blank = true;
    double percent = 0.0;
    if( read_double( r.cert, percent, blank ) && !blank && (percent != 0.0) )
      any_correlated = true;
  }//for( const AnchorRow &r : m_anchors )

  const bool show = (m_uncertOptions && !m_uncertOptions->isHidden() && !m_anchors.empty()
                     && !any_correlated);
  m_corrInertNote->setHidden( !show );
}//updateCorrelationNote()


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
