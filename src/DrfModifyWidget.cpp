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
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/DrfModifyWidget.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/MakeMcResponseForDrf.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace Wt;
using namespace std;


DrfModifyWidget::DrfModifyWidget( InterSpec *viewer,
                                  std::shared_ptr<const DetectorPeakResponse> drf )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_orig( drf ),
    m_tabMenu( nullptr ),
    m_tabStack( nullptr ),
    m_name( nullptr ),
    m_description( nullptr ),
    m_mcTool( nullptr ),
    m_fwhmTool( nullptr ),
    m_bandTable( nullptr ),
    m_addBand( nullptr ),
    m_removeBand( nullptr ),
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

    m_tabMenu->addItem( WString::tr("dmw-tab-name"), std::move(panelOwned) );
  }

  // --- Tab: Geometry & MC characterization ---------------------------------
  {
    auto toolOwned = make_unique<MakeMcResponseForDrf>( m_interspec, m_orig );
    m_mcTool = toolOwned.get();
    m_tabMenu->addItem( WString::tr("dmw-tab-geometry"), std::move(toolOwned) );
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
      // MakeFwhmForDrf edits a (non-const) DRF; give it a private clone - we
      //  only read its fitted coefficients back out on apply().
      shared_ptr<DetectorPeakResponse> fwhm_seed
          = m_orig ? make_shared<DetectorPeakResponse>( *m_orig ) : nullptr;
      auto toolOwned = make_unique<MakeFwhmForDrf>( true, m_interspec, fwhm_seed );
      m_fwhmTool = toolOwned.get();
      m_tabMenu->addItem( WString::tr("dmw-tab-fwhm"), std::move(toolOwned) );
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

  m_tabMenu->select( 0 );

  WText *note = new WText( WString::tr("dmw-export-note") );
  note->addStyleClass( "DrfModifyNote" );
  note->setInline( false );
  layout->addWidget( std::unique_ptr<WText>(note), 2, 0 );
}//DrfModifyWidget constructor


DrfModifyWidget::~DrfModifyWidget()
{
}


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

  // FWHM: pull the fitted coefficients (if any) without triggering the tool's
  //  own detector-changed emit.
  const shared_ptr<MakeFwhmForDrf::ToolState> fwhm
      = m_fwhmTool ? m_fwhmTool->currentState() : nullptr;
  if( fwhm && !fwhm->m_parameters.empty() && m_fwhmTool->isValidFwhm() )
  {
    const auto form = DetectorPeakResponse::ResolutionFnctForm(
                          std::max( 0, fwhm->m_fwhm_index ) );
    working->setFwhmCoefficients( fwhm->m_parameters, form );
  }//if( have a valid FWHM fit )

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

  show();
  resizeToFitOnScreen();
  centerWindow();
}//DrfModifyWindow constructor
