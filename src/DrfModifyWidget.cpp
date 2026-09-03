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
#include "InterSpec/InterSpecApp.h"
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
    m_bandTable( nullptr ),
    m_addBand( nullptr ),
    m_removeBand( nullptr ),
    m_anchorTable( nullptr ),
    m_addAnchor( nullptr ),
    m_removeAnchor( nullptr ),
    m_anchorRefDistance( nullptr ),
    m_anchorDefaultUncert( nullptr ),
    m_anchorIsAbsolute( false ),
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

  // --- Tab: General (name / description / measured-curve anchors) ----------
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
    m_name->changed().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
    nameLabel->setBuddy( m_name );

    WLabel *descLabel = idGrid->addNew<WLabel>( WString::tr("Description") );
    m_description = idGrid->addNew<WTextArea>();
    m_description->setColumns( 40 );
    m_description->setRows( 3 );
    m_description->setText( WString::fromUTF8( m_orig ? m_orig->description() : string() ) );
    m_description->changed().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
    descLabel->setBuddy( m_description );

    // Measured-curve anchor editor: the source-of-truth editor for the reference curve the
    //  response is anchored to; edits are folded into the working DRF on apply().  Shown whenever
    //  the DRF states a geometry and its efficiency is a table of measured points:
    //
    //   - raw `measuredPoints` (an ANGLE import): absolute efficiencies at a reference distance;
    //   - failing that, an energy/efficiency-pair efficiency curve (a GADRAS Efficiency.csv):
    //     the same measured table, just without any uncertainties, and intrinsic rather than
    //     absolute - so no reference distance applies.
    const shared_ptr<const MeasuredDrfPoints> measured
        = m_orig ? m_orig->measuredPoints() : nullptr;
    const bool have_measured = (measured && !measured->empty());

    const shared_ptr<const DetectorEfficiencyCurve> eff_curve
        = m_orig ? m_orig->efficiencyCurve() : nullptr;
    const bool have_pairs = (!have_measured && eff_curve
                    && (eff_curve->form() == DetectorPeakResponse::kEnergyEfficiencyPairs)
                    && (eff_curve->energyEfficiencies().size() >= 2));

    m_anchorIsAbsolute = (m_orig && (m_orig->geometryType()
                                     == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute));

    if( m_geometry && (have_measured || have_pairs) )
    {
      WGroupBox *box = panel->addNew<WGroupBox>( WString::tr("dmw-anchor-title") );
      box->addStyleClass( "DrfModifyAnchor" );
      box->addNew<WText>( WString::tr("dmw-anchor-help") )->setInline( false );

      // Reference distance the curve is anchored at (cm) - only meaningful for absolute
      //  efficiencies; an intrinsic curve is per-gamma-striking-the-face, at no distance.
      double refDistCm = 0.0;
      if( have_measured )
      {
        for( const MeasuredEffPoint &p : measured->points() )
        {
          if( p.distance > 0.0f )
          {
            refDistCm = p.distance / PhysicalUnits::cm;
            break;
          }
        }//for( find a positive reference distance )
      }//if( have_measured )
      if( refDistCm <= 0.0 )
        refDistCm = m_orig->absoluteEfficiencyDistance() / PhysicalUnits::cm;

      WContainerWidget *distRow = box->addNew<WContainerWidget>();
      distRow->addStyleClass( "DrfModifyRow" );
      distRow->addNew<WLabel>( WString::tr("dmw-anchor-ref-dist") );
      m_anchorRefDistance = distRow->addNew<WLineEdit>();
      m_anchorRefDistance->setTextSize( 8 );
      m_anchorRefDistance->changed().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
      if( refDistCm > 0.0 )
      {
        char buf[32];
        snprintf( buf, sizeof(buf), "%.4g", refDistCm );
        m_anchorRefDistance->setText( buf );
      }
      distRow->addNew<WLabel>( WString::tr("dmw-anchor-cm") );
      distRow->setHidden( !m_anchorIsAbsolute );

      WContainerWidget *defRow = box->addNew<WContainerWidget>();
      defRow->addStyleClass( "DrfModifyRow" );
      defRow->addNew<WLabel>( WString::tr("dmw-anchor-default-uncert") );
      m_anchorDefaultUncert = defRow->addNew<WLineEdit>();
      m_anchorDefaultUncert->setTextSize( 6 );
      m_anchorDefaultUncert->changed().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
      // Only a fallback for blank per-point cells.  Left blank when the source stated no
      //  uncertainties at all (a GADRAS Efficiency.csv), rather than inventing one.
      if( have_measured )
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
      m_addAnchor->clicked().connect( this, [this](){ addAnchorRow( 0.0f, 0.0f, 0.0f ); scheduleUndoRedoStep(); } );
      m_removeAnchor = btns->addNew<WPushButton>( WString::tr("dmw-anchor-remove") );
      m_removeAnchor->addStyleClass( "LinkBtn" );
      m_removeAnchor->clicked().connect( this, [this](){ removeAnchorRow(); scheduleUndoRedoStep(); } );

      if( have_measured )
      {
        for( const MeasuredEffPoint &p : measured->points() )
          addAnchorRow( p.energy, p.efficiency, p.fracStatUncert );
      }else
      {
        // Energies are stored in `PhysicalUnits::keV` units, which is 1.0 - i.e. already keV,
        //  the same units the rows and MeasuredEffPoint use.
        for( const DetectorPeakResponse::EnergyEfficiencyPair &p : eff_curve->energyEfficiencies() )
          addAnchorRow( p.energy, p.efficiency, 0.0f );   //no uncertainty to show
      }//if( have_measured ) / else

      m_removeAnchor->setEnabled( !m_anchors.empty() );
    }//if( the DRF states a geometry, and its efficiency is a measured table )

    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-name"), std::move(panelOwned) );
    make_item_selectable( m_tabMenu, item );
  }

  // --- Tab: Geom & MC characterization -------------------------------------
  {
    auto toolOwned = make_unique<MakeMcResponseForDrf>( m_interspec, m_orig );
    m_mcTool = toolOwned.get();

    // ContentLoading::Eager: this tool posts work to a worker thread and gets back to itself with
    //  `findById(...)`; a Lazy tab parks its contents in `WMenuItem::uContents_`, outside the
    //  widget tree, where `findById` can not see it - see the FWHM tab below.
    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-geometry"), std::move(toolOwned),
                                          ContentLoading::Eager );
    make_item_selectable( m_tabMenu, item );
  }

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
    m_addBand->clicked().connect( this, [this](){ addBandRow( 0.0f, 0.0f, 0.05f ); scheduleUndoRedoStep(); } );
    m_removeBand = btns->addNew<WPushButton>( WString::tr("dmw-remove-band") );
    m_removeBand->addStyleClass( "LinkBtn" );
    m_removeBand->clicked().connect( this, [this](){ removeBandRow(); scheduleUndoRedoStep(); } );

    // Seed rows from any existing piecewise-band uncertainty.
    const shared_ptr<const DetectorEfficiencyUncert> uncert
        = m_orig ? m_orig->efficiencyUncert() : nullptr;
    if( uncert )
    {
      for( const EffUncertBand &b : uncert->bands() )
        addBandRow( b.lowerEnergy, b.upperEnergy, b.fractionalUncert );
    }
    m_removeBand->setEnabled( !m_bands.empty() );

    WMenuItem *item = m_tabMenu->addItem( WString::tr("dmw-tab-uncert"), std::move(panelOwned) );
    make_item_selectable( m_tabMenu, item );
  }

  m_tabMenu->itemSelected().connect( this, &DrfModifyWidget::handleTabSelected );

  // Both child tools hand their state changes here, so one undo step covers the whole dialog.
  m_mcTool->userChanged().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
  if( m_fwhmTool )
  {
    m_fwhmTool->setOwnerHandlesUndoRedo( true );
    m_fwhmTool->stateChanged().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
  }//if( m_fwhmTool )

  m_tabMenu->select( 0 );

  auto noteOwned = make_unique<WText>( WString::tr("dmw-export-note") );
  noteOwned->addStyleClass( "DrfModifyNote" );
  noteOwned->setInline( false );
  if( narrow_layout )
    layout->addWidget( std::move(noteOwned), 2, 0 );        //below menu-row and stack-row
  else
    layout->addWidget( std::move(noteOwned), 1, 0, 1, 2 );  //below, spanning menu and stack columns
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
  return m_mcTool->validationChanged();
}


bool DrfModifyWidget::needsMcResponse() const
{
  return (!m_orig || !m_orig->isValid());
}


std::shared_ptr<const DetectorPeakResponse> DrfModifyWidget::originalDrf() const
{
  return m_orig;
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
  for( WLineEdit *edit : { r.lower, r.upper, r.frac } )
    edit->changed().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
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
  for( WLineEdit *edit : { r.energy, r.eff, r.uncert } )
    edit->changed().connect( this, &DrfModifyWidget::scheduleUndoRedoStep );
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
  if( m_anchorIsAbsolute && (refDistCm <= 0.0) )
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
    // Write the points back as whatever they ARE - re-interpreting an intrinsic curve as absolute
    //  at a distance would rescale the detector by its solid angle.
    const DetectorPeakResponse::EffGeometryType geom_type = m_anchorIsAbsolute
                          ? DetectorPeakResponse::EffGeometryType::FarFieldAbsolute
                          : DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;

    working.setEfficiencyPoints( effpts, diameter, refDist, geom_type );

    // `MeasuredDrfPoints` are read elsewhere (grounding, transfer anchors) as absolute
    //  efficiencies at `distance`, so only absolute rows may be recorded as them.
    if( m_anchorIsAbsolute )
    {
      auto meas = make_shared<MeasuredDrfPoints>();
      meas->setPoints( measpts );
      working.setMeasuredPoints( meas );
    }
  }catch( std::exception & )
  {
    //Malformed edits (e.g. duplicate/degenerate energies): keep the seed curve.
  }
}//applyAnchorEdits(...)


bool DrfModifyWidget::ToolState::operator==( const ToolState &rhs ) const
{
  if( (name != rhs.name) || (description != rhs.description) || (tabIndex != rhs.tabIndex)
     || (bands != rhs.bands) || (anchors != rhs.anchors)
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

  for( const BandRow &r : m_bands )
    state->bands.push_back( { r.lower->text().toUTF8(), r.upper->text().toUTF8(),
                              r.frac->text().toUTF8() } );

  for( const AnchorRow &r : m_anchors )
    state->anchors.push_back( { r.energy->text().toUTF8(), r.eff->text().toUTF8(),
                                r.uncert->text().toUTF8() } );

  if( m_anchorRefDistance )
    state->anchorRefDistance = m_anchorRefDistance->text().toUTF8();
  if( m_anchorDefaultUncert )
    state->anchorDefaultUncert = m_anchorDefaultUncert->text().toUTF8();

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

  while( !m_bands.empty() )
    removeBandRow();
  for( const std::array<string,3> &b : state->bands )
  {
    addBandRow( 0.0f, 0.0f, 0.0f );
    m_bands.back().lower->setText( WString::fromUTF8(b[0]) );
    m_bands.back().upper->setText( WString::fromUTF8(b[1]) );
    m_bands.back().frac->setText( WString::fromUTF8(b[2]) );
  }//for( const std::array<string,3> &b : state->bands )

  if( m_anchorTable )
  {
    while( !m_anchors.empty() )
      removeAnchorRow();
    for( const std::array<string,3> &a : state->anchors )
    {
      addAnchorRow( 0.0f, 0.0f, 0.0f );
      m_anchors.back().energy->setText( WString::fromUTF8(a[0]) );
      m_anchors.back().eff->setText( WString::fromUTF8(a[1]) );
      m_anchors.back().uncert->setText( WString::fromUTF8(a[2]) );
    }//for( const std::array<string,3> &a : state->anchors )

    m_anchorRefDistance->setText( WString::fromUTF8(state->anchorRefDistance) );
    m_anchorDefaultUncert->setText( WString::fromUTF8(state->anchorDefaultUncert) );
  }//if( m_anchorTable )

  m_mcTool->setState( state->mc );
  if( m_fwhmTool && state->fwhm )
    m_fwhmTool->setState( state->fwhm );

  if( (state->tabIndex >= 0) && (state->tabIndex < m_tabMenu->count()) )
    m_tabMenu->select( state->tabIndex );

  // Restoring is not a user edit, so nothing here should become the next undo step.
  m_renderFlags.clear( RenderActions::AddUndoRedoStep );
  m_restoringState = false;
}//DrfModifyWidget::setState(...)


void DrfModifyWidget::scheduleUndoRedoStep()
{
  if( m_restoringState )
    return;

  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void scheduleUndoRedoStep()


void DrfModifyWidget::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool add_step = m_renderFlags.test( RenderActions::AddUndoRedoStep );
  m_renderFlags = Wt::WFlags<RenderActions>();

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

  // Monte-Carlo response: attach + ground (to the DRF's measured points if it
  //  has them; otherwise the MC tool already grounded to a sampled curve).
  const shared_ptr<const ceelo::DetectorResponse> resp = m_mcTool->generatedResponse();
  if( resp )
  {
    // A DRF built from geometry alone has no efficiency curve of its own, so
    //  sample the response into one - the Monte-Carlo "backbone" points.  Every
    //  EffEval query would already dispatch to the response, but without a curve
    //  the DRF is not valid, cannot be serialized, and cannot be exported; a
    //  Detector.dat or .detx import with no measured efficiency lands here.
    //  Done BEFORE setCeeloResponse: setEfficiencyPoints recomputes the hash and
    //  resets flags, and the response should be attached to the finished curve.
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
  }//if( resp )

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
  use->clicked().connect( m_tool, &DrfModifyWidget::apply );

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
