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
#include <Wt/WLineEdit.h>
#include <Wt/WTextArea.h>
#include <Wt/WMenuItem.h>
#include <Wt/WGroupBox.h>
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
#include "InterSpec/DetectorEfficiency.h"
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
    m_uncertHelp( nullptr ),
    m_pointsEditor( nullptr ),
    m_covEditor( nullptr ),
    m_anchorTable( nullptr ),
    m_addAnchor( nullptr ),
    m_removeAnchor( nullptr ),
    m_anchorRefDistance( nullptr ),
    m_anchorDefaultUncert( nullptr ),
    m_anchorIsAbsolute( false ),
    m_covTable( nullptr ),
    m_addEnergy( nullptr ),
    m_removeEnergy( nullptr ),
    m_covEnergies(),
    m_covMatrix(),
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

  // What kind of efficiency the DRF carries decides which Uncertainty editor is the right one, and
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

  // --- Tab: Efficiency uncertainty -----------------------------------------
  //  Two editors, exactly one shown (see updateUncertEditorVisibility):
  //   - the measured-point editor (energy / efficiency / stat % / cert % / source), when the DRF
  //     carries measured points or the mode is Geometry Modeled - these feed a CeeLo grounding;
  //   - the σ/ρ node-covariance matrix, otherwise (an abstract Flat-Disk uncertainty).
  {
    auto panelOwned = make_unique<WContainerWidget>();
    WContainerWidget *panel = panelOwned.get();
    panel->addStyleClass( "DrfModifyPanel DrfModifyUncert" );

    // Wording is set to match the visible editor in updateUncertEditorVisibility().
    m_uncertHelp = panel->addNew<WText>();
    m_uncertHelp->setInline( false );

    // --- measured-point editor -------------------------------------------
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
      // Only a fallback for blank per-point stat cells; left blank when the source stated no
      //  uncertainties at all (a GADRAS Efficiency.csv), rather than inventing one.
      if( have_measured )
        m_anchorDefaultUncert->setText( "5" );
      defRow->addNew<WLabel>( WString::tr("dmw-anchor-percent") );

      m_anchorTable = m_pointsEditor->addNew<WTable>();
      m_anchorTable->addStyleClass( "DrfModifyAnchorTable" );
      m_anchorTable->elementAt(0,0)->addNew<WText>( WString::tr("dmw-anchor-energy") );
      m_anchorTable->elementAt(0,1)->addNew<WText>( WString::tr("dmw-anchor-eff") );
      m_anchorTable->elementAt(0,2)->addNew<WText>( WString::tr("dmw-anchor-stat") );
      if( m_anchorIsAbsolute )
      {
        // Certificate uncertainty is 100%-correlated within a source; the source column supplies
        //  the grouping key that makes that common-mode block meaningful.
        m_anchorTable->elementAt(0,3)->addNew<WText>( WString::tr("dmw-anchor-cert") );
        m_anchorTable->elementAt(0,4)->addNew<WText>( WString::tr("dmw-anchor-source") );
      }//if( m_anchorIsAbsolute )

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
        // Energies are stored in `PhysicalUnits::keV` units (== 1.0) - already keV, the units the
        //  rows use.  A pairs curve carries no uncertainties or source.
        for( const DetectorPeakResponse::EnergyEfficiencyPair &p : eff_curve->energyEfficiencies() )
          addAnchorRow( p.energy, p.efficiency, 0.0f, 0.0f, string() );
      }//if( have_measured ) / else if( have_pairs )
      m_removeAnchor->setEnabled( !m_anchors.empty() );

      panel->addWidget( std::move(ptsOwned) );
    }//measured-point editor

    // --- σ/ρ covariance-matrix editor ------------------------------------
    {
      auto covOwned = make_unique<WContainerWidget>();
      m_covEditor = covOwned.get();
      m_covEditor->addStyleClass( "DrfModifyCov" );

      m_covTable = m_covEditor->addNew<WTable>();
      m_covTable->addStyleClass( "DrfModifyBandTable DrfModifyCovTable" );

      WContainerWidget *btns = m_covEditor->addNew<WContainerWidget>();
      m_addEnergy = btns->addNew<WPushButton>( WString::tr("dmw-add-energy") );
      m_addEnergy->addStyleClass( "LinkBtn" );
      m_addEnergy->clicked().connect( this, [this](){ addEnergyRow(); markEdited(); } );
      m_removeEnergy = btns->addNew<WPushButton>( WString::tr("dmw-remove-energy") );
      m_removeEnergy->addStyleClass( "LinkBtn" );
      m_removeEnergy->clicked().connect( this, [this](){ removeEnergyRow(); markEdited(); } );

      const shared_ptr<const DetectorEfficiencyUncert> uncert
          = m_orig ? m_orig->efficiencyUncert() : nullptr;
      seedCovFromUncert( uncert );  //fills the shadow and does the initial rebuildCovTable()

      panel->addWidget( std::move(covOwned) );
    }//covariance-matrix editor

    updateUncertEditorVisibility();

    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-uncert"), std::move(panelOwned) );
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
    WLineEdit *edit = m_anchorTable->elementAt(row,col)->addNew<WLineEdit>();
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
  //  apply.  cert/source only exist for an absolute reference curve (see m_anchorIsAbsolute).
  AnchorRow r;
  r.energy = make_edit( 0, energy, false );
  r.eff    = make_edit( 1, efficiency, false );
  r.stat   = make_edit( 2, fracStatUncert, true );
  r.cert   = nullptr;
  r.source = nullptr;
  if( m_anchorIsAbsolute )
  {
    r.cert = make_edit( 3, fracCertUncert, true );
    r.source = m_anchorTable->elementAt(row,4)->addNew<WLineEdit>();
    r.source->setTextSize( 10 );
    r.source->setText( WString::fromUTF8(sourceKey) );
    r.source->changed().connect( this, &DrfModifyWidget::markEdited );
  }//if( m_anchorIsAbsolute )
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

    float certFrac = 0.0f;
    string srcKey;
    if( m_anchorIsAbsolute && r.cert && r.source )
    {
      const string cs = r.cert->text().toUTF8();
      if( !cs.empty() )
      {
        try{ certFrac = 0.01f * std::stof( cs ); }catch( std::exception & ){ certFrac = 0.0f; }
      }
      if( certFrac < 0.0f )
        certFrac = 0.0f;
      srcKey = r.source->text().toUTF8();
    }//if( m_anchorIsAbsolute )

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

  const float diameter = working.detectorDiameter();
  if( diameter <= 0.0f )
    return;  //a far-field curve needs a positive diameter; leave the DRF as-is

  try
  {
    // Write the points back as whatever they ARE - re-interpreting an intrinsic curve as absolute
    //  at a distance would rescale the detector by its solid angle.
    const DetectorPeakResponse::EffGeometryType geom_type = m_anchorIsAbsolute
                          ? DetectorPeakResponse::EffGeometryType::FarFieldAbsolute
                          : DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;

    working.setEfficiencyPoints( effpts, diameter, refDist, geom_type );

    // `MeasuredDrfPoints` are read elsewhere (grounding, transfer anchors) as absolute
    //  efficiencies at `distance`, so only absolute rows may be recorded as them.
    auto meas = make_shared<MeasuredDrfPoints>();
    meas->setPoints( measpts );
    if( m_anchorIsAbsolute )
      working.setMeasuredPoints( meas );

    // Drive the efficiency uncertainty from the edited stat/cert/source structure (C[i][j] =
    //  delta_ij*stat_i^2 + cert_i*cert_j*[same source]); a CeeLo response, when attached, overrides
    //  this at query time.  Keeps the fromPointUncerts covariance setEfficiencyPoints just built
    //  only when the richer form cannot be (fewer than two distinct energies).
    const shared_ptr<DetectorEfficiencyUncert> uncert = meas->toEfficiencyUncert();
    if( uncert && !uncert->isEmpty() )
      working.setEfficiencyUncert( uncert );
  }catch( std::exception & )
  {
    //Malformed edits (e.g. duplicate/degenerate energies): keep the seed curve.
  }
}//applyAnchorEdits(...)


bool DrfModifyWidget::ToolState::operator==( const ToolState &rhs ) const
{
  if( (name != rhs.name) || (description != rhs.description) || (tabIndex != rhs.tabIndex)
     || (geometryModeled != rhs.geometryModeled)
     || (covEnergies != rhs.covEnergies) || (covMatrix != rhs.covMatrix)
     || (anchors != rhs.anchors)
     || (anchorRefDistance != rhs.anchorRefDistance)
     || (anchorDefaultUncert != rhs.anchorDefaultUncert)
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
  state->covEnergies = m_covEnergies;
  state->covMatrix = m_covMatrix;

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

  // Covariance shadow, then the table it drives.
  m_covEnergies = state->covEnergies;
  m_covMatrix = state->covMatrix;
  rebuildCovTable();

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

  updateUncertEditorVisibility();

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

  // Uncertainty tab: apply whichever editor is showing.  The measured-point editor rebuilds the
  //  far-field efficiency + measured points + node covariance; the σ/ρ editor sets the node
  //  covariance directly.  Done before the MC step so a grounding sees the edited points.
  if( pointsEditorVisible() )
    applyAnchorEdits( *working );
  else
    applyCovarianceEdits( *working );

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
  updateUncertEditorVisibility();
  markEdited();
}//handleModeToggle()


bool DrfModifyWidget::geometryModeled() const
{
  return m_geometryModeled;
}


bool DrfModifyWidget::pointsEditorVisible() const
{
  return m_origHasPoints || m_geometryModeled;
}


void DrfModifyWidget::updateUncertEditorVisibility()
{
  const bool showPoints = pointsEditorVisible();
  if( m_pointsEditor )
    m_pointsEditor->setHidden( !showPoints );
  if( m_covEditor )
    m_covEditor->setHidden( showPoints );
  if( m_uncertHelp )
    m_uncertHelp->setText( WString::tr( showPoints ? "dmw-uncert-help-points"
                                                    : "dmw-uncert-help-cov" ) );
}//updateUncertEditorVisibility()


void DrfModifyWidget::seedCovFromUncert( const std::shared_ptr<const DetectorEfficiencyUncert> &uncert )
{
  m_covEnergies.clear();
  m_covMatrix.clear();

  if( uncert && uncert->hasNodeCovariance() )
  {
    const vector<float> &en = uncert->covarianceEnergies();
    const vector<float> &cov = uncert->covarianceMatrix();
    const size_t n = en.size();
    if( (n >= 1) && (cov.size() == n*n) )
    {
      m_covEnergies.assign( begin(en), end(en) );
      m_covMatrix.assign( begin(cov), end(cov) );
    }
  }//if( have node covariance )

  rebuildCovTable();
}//seedCovFromUncert(...)


void DrfModifyWidget::rebuildCovTable()
{
  if( !m_covTable )
    return;

  m_covTable->clear();  //wipes all cells and the WLineEdits/WTexts they hold

  const size_t n = m_covEnergies.size();

  auto fmt = []( char *buf, size_t buflen, const char *spec, double v ){
    snprintf( buf, buflen, spec, v );
  };

  // Header row: a corner label, then one echoed node energy per matrix column.
  m_covTable->elementAt(0,0)->addNew<WText>( WString::tr("dmw-cov-corner") );
  for( size_t j = 0; j < n; ++j )
  {
    char buf[32];
    fmt( buf, sizeof(buf), "%.4g", m_covEnergies[j] );
    m_covTable->elementAt(0, static_cast<int>(j)+1)->addNew<WText>( WString::fromUTF8(buf) );
  }

  for( size_t i = 0; i < n; ++i )
  {
    const int trow = static_cast<int>(i) + 1;

    // Column 0: editable node energy (keV).
    {
      char buf[32];
      fmt( buf, sizeof(buf), "%.4g", m_covEnergies[i] );
      WLineEdit *e = m_covTable->elementAt(trow,0)->addNew<WLineEdit>();
      e->setTextSize( 7 );
      e->setText( WString::fromUTF8(buf) );
      e->changed().connect( this, [this,i,e](){ covEnergyChanged( i, e->text().toUTF8() ); } );
    }

    const double Cii = m_covMatrix[i*n + i];
    const double sigma_i = (Cii > 0.0) ? std::sqrt(Cii) : 0.0;

    for( size_t j = 0; j < n; ++j )
    {
      const int tcol = static_cast<int>(j) + 1;
      const double Cjj = m_covMatrix[j*n + j];
      const double sigma_j = (Cjj > 0.0) ? std::sqrt(Cjj) : 0.0;

      if( j == i )
      {
        // Diagonal: 1-sigma fractional uncertainty as a percent (editable).
        char buf[32];
        fmt( buf, sizeof(buf), "%.4g", 100.0*sigma_i );
        WLineEdit *e = m_covTable->elementAt(trow,tcol)->addNew<WLineEdit>();
        e->setTextSize( 6 );
        e->addStyleClass( "DrfCovDiag" );
        e->setText( WString::fromUTF8(buf) );
        e->changed().connect( this, [this,i,e](){ covSigmaChanged( i, e->text().toUTF8() ); } );
      }else
      {
        const double Cij = m_covMatrix[i*n + j];
        const double rho = ((sigma_i > 0.0) && (sigma_j > 0.0)) ? (Cij/(sigma_i*sigma_j)) : 0.0;
        char buf[32];
        fmt( buf, sizeof(buf), "%.3g", rho );
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

  if( m_removeEnergy )
    m_removeEnergy->setEnabled( n > 0 );
  if( m_addEnergy )
    m_addEnergy->setEnabled( n < DetectorEfficiencyUncert::sm_maxCovarianceNodes );
}//rebuildCovTable()


void DrfModifyWidget::addEnergyRow()
{
  const size_t oldN = m_covEnergies.size();
  if( oldN >= DetectorEfficiencyUncert::sm_maxCovarianceNodes )
    return;
  const size_t n = oldN + 1;

  // Append past the current last (or a default when empty), so this adds rather than reorders - the
  //  shadow is display-ordered and only sorted at apply.
  const double newEnergy = (oldN > 0) ? (m_covEnergies.back() + 100.0) : 100.0;
  const double defSigma = 0.05;  //5% default, zero correlation to the existing nodes

  vector<double> newMat( n*n, 0.0 );
  for( size_t r = 0; r < oldN; ++r )
    for( size_t c = 0; c < oldN; ++c )
      newMat[r*n + c] = m_covMatrix[r*oldN + c];
  newMat[(n-1)*n + (n-1)] = defSigma*defSigma;

  m_covEnergies.push_back( newEnergy );
  m_covMatrix.swap( newMat );

  rebuildCovTable();
}//addEnergyRow()


void DrfModifyWidget::removeEnergyRow()
{
  const size_t oldN = m_covEnergies.size();
  if( oldN == 0 )
    return;
  const size_t n = oldN - 1;

  vector<double> newMat( n*n, 0.0 );
  for( size_t r = 0; r < n; ++r )
    for( size_t c = 0; c < n; ++c )
      newMat[r*n + c] = m_covMatrix[r*oldN + c];

  m_covEnergies.pop_back();
  m_covMatrix.swap( newMat );

  rebuildCovTable();
}//removeEnergyRow()


void DrfModifyWidget::covSigmaChanged( const std::size_t i, const std::string &text )
{
  const size_t n = m_covEnergies.size();
  if( i >= n )
    return;

  double pct = 0.0;
  try{ pct = std::stod( text ); }
  catch( std::exception & ){ m_renderFlags |= RenderActions::RebuildCovTable; scheduleRender(); return; }

  double s_new = 0.01 * pct;
  if( s_new < 0.0 )
    s_new = 0.0;

  const double Cii = m_covMatrix[i*n + i];
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
      m_covMatrix[i*n + k] *= f;
      m_covMatrix[k*n + i] *= f;
    }
  }//if( s_old > 0.0 )
  m_covMatrix[i*n + i] = s_new * s_new;

  m_renderFlags |= RenderActions::RebuildCovTable;
  markEdited();
}//covSigmaChanged(...)


void DrfModifyWidget::covRhoChanged( const std::size_t i, const std::size_t j, const std::string &text )
{
  const size_t n = m_covEnergies.size();
  if( (i >= n) || (j >= n) || (i == j) )
    return;

  double rho = 0.0;
  try{ rho = std::stod( text ); }
  catch( std::exception & ){ m_renderFlags |= RenderActions::RebuildCovTable; scheduleRender(); return; }
  if( rho < -1.0 ) rho = -1.0;
  if( rho >  1.0 ) rho =  1.0;

  const double sigma_i = std::sqrt( std::max(0.0, m_covMatrix[i*n + i]) );
  const double sigma_j = std::sqrt( std::max(0.0, m_covMatrix[j*n + j]) );
  const double Cij = rho * sigma_i * sigma_j;
  m_covMatrix[i*n + j] = Cij;
  m_covMatrix[j*n + i] = Cij;

  m_renderFlags |= RenderActions::RebuildCovTable;
  markEdited();
}//covRhoChanged(...)


void DrfModifyWidget::covEnergyChanged( const std::size_t i, const std::string &text )
{
  const size_t n = m_covEnergies.size();
  if( i >= n )
    return;

  double energy = 0.0;
  try{ energy = std::stod( text ); }
  catch( std::exception & ){ m_renderFlags |= RenderActions::RebuildCovTable; scheduleRender(); return; }
  if( energy <= 0.0 )
  {
    // Reject; the rebuild restores the previous value.  (Duplicate/non-increasing energies are
    //  caught at apply, since the shadow stays in display order until then.)
    m_renderFlags |= RenderActions::RebuildCovTable;
    scheduleRender();
    return;
  }//if( energy <= 0.0 )

  m_covEnergies[i] = energy;
  m_renderFlags |= RenderActions::RebuildCovTable;  //the header echo follows the edit
  markEdited();
}//covEnergyChanged(...)


void DrfModifyWidget::applyCovarianceEdits( DetectorPeakResponse &working )
{
  if( !m_covEditor )
    return;  //covariance editor was not built

  const size_t n = m_covEnergies.size();

  // Empty -> clear any node covariance (keeping other uncert content, e.g. coef covariance).
  if( n == 0 )
  {
    const shared_ptr<const DetectorEfficiencyUncert> existing = working.efficiencyUncert();
    if( existing && existing->hasNodeCovariance() )
    {
      auto uncert = make_shared<DetectorEfficiencyUncert>( *existing );
      uncert->setNodeCovariance( {}, {} );  //clears the node covariance
      working.setEfficiencyUncert( uncert->isEmpty() ? nullptr : uncert );
    }
    return;
  }//if( n == 0 )

  // Sort energies ascending, permuting the covariance to match; reject duplicate/non-increasing.
  vector<size_t> order( n );
  std::iota( begin(order), end(order), size_t(0) );
  std::sort( begin(order), end(order),
             [this]( size_t a, size_t b ){ return m_covEnergies[a] < m_covEnergies[b]; } );

  vector<float> energies( n );
  for( size_t k = 0; k < n; ++k )
    energies[k] = static_cast<float>( m_covEnergies[order[k]] );
  for( size_t k = 1; k < n; ++k )
  {
    if( !(energies[k] > energies[k-1]) )
    {
      passMessage( WString::tr("dmw-err-cov-energies"), WarningWidget::WarningMsgHigh );
      return;
    }
  }//for( check strictly increasing )

  vector<float> cov( n*n, 0.0f );
  for( size_t r = 0; r < n; ++r )
    for( size_t c = 0; c < n; ++c )
      cov[r*n + c] = static_cast<float>( m_covMatrix[order[r]*n + order[c]] );

  try
  {
    const shared_ptr<const DetectorEfficiencyUncert> existing = working.efficiencyUncert();
    auto uncert = existing ? make_shared<DetectorEfficiencyUncert>( *existing )
                           : make_shared<DetectorEfficiencyUncert>();
    uncert->setNodeCovariance( energies, cov );
    working.setEfficiencyUncert( uncert->isEmpty() ? nullptr : uncert );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("dmw-err-cov-invalid").arg(e.what()), WarningWidget::WarningMsgHigh );
  }
}//applyCovarianceEdits(...)


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
