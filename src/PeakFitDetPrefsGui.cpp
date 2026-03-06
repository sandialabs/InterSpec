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
#include <cassert>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WApplication>

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/PeakFitDetPrefsGui.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;


PeakFitDetPrefsGui::PeakFitDetPrefsGui( InterSpec *viewer, const bool compactMode,
                                         WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_viewer( viewer ),
    m_compactMode( compactMode ),
    m_programmaticUpdate( false ),
    m_collapsedDiv( nullptr ),
    m_expandedDiv( nullptr ),
    m_detTypeCombo( nullptr ),
    m_fwhmMethodCombo( nullptr ),
    m_skewTypeCombo( nullptr ),
    m_skewParamsDiv( nullptr ),
    m_roiIndepCb( nullptr ),
    m_fitSkewLink( nullptr ),
    m_sourceLabel( nullptr ),
    m_pendingFwhmMethod( PeakFitDetPrefs::FwhmMethod::Normal )
{
  assert( m_viewer );

  for( int i = 0; i < 4; ++i )
  {
    m_lowerSkewSpin[i] = nullptr;
    m_upperSkewSpin[i] = nullptr;
  }

  addStyleClass( "PeakFitDetPrefsGui" );

  wApp->useStyleSheet( "InterSpec_resources/PeakFitDetPrefsGui.css" );

  

  init();
}//PeakFitDetPrefsGui constructor


PeakFitDetPrefsGui::~PeakFitDetPrefsGui()
{
}


void PeakFitDetPrefsGui::init()
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
    app->useMessageResourceBundle( "PeakFitDetPrefsGui" );
  
  if( m_compactMode )
  {
    // Collapsed strip: expand icon at top, vertical text, expand icon at bottom
    m_collapsedDiv = new WContainerWidget( this );
    m_collapsedDiv->addStyleClass( "PeakFitDetPrefsCollapsed" );
    m_collapsedDiv->clicked().connect( this, &PeakFitDetPrefsGui::toggleExpanded );

    WImage *topIcon = new WImage( WLink( "InterSpec_resources/images/expand_left.svg" ), m_collapsedDiv );
    topIcon->addStyleClass( "PfdpgExpandIcon PfdpgExpandIconTop" );

    WText *collapseLabel = new WText( WString::tr( "pfdpg-collapsed-label" ), m_collapsedDiv );
    collapseLabel->addStyleClass( "PeakFitDetPrefsCollapsedLabel" );

    WImage *bottomIcon = new WImage( WLink( "InterSpec_resources/images/expand_left.svg" ), m_collapsedDiv );
    bottomIcon->addStyleClass( "PfdpgExpandIcon PfdpgExpandIconBottom" );
  }//if( m_compactMode )

  // Expanded controls container
  m_expandedDiv = new WContainerWidget( this );
  m_expandedDiv->addStyleClass( "PeakFitDetPrefsExpanded" );

  if( m_compactMode )
  {
    m_expandedDiv->hide();

    // Clickable header to collapse - stays outside scrollable body
    WText *header = new WText( WString::tr( "pfdpg-collapsed-label" ), m_expandedDiv );
    header->addStyleClass( "PeakFitDetPrefsExpandedHeader" );
    header->clicked().connect( this, &PeakFitDetPrefsGui::toggleExpanded );
  }else
  {
    if( app )
      app->useMessageResourceBundle( "PeakEdit" );
  }//if( m_compactMode )

  // Scrollable body div - holds all controls so header stays fixed
  WContainerWidget *contentDiv = new WContainerWidget( m_expandedDiv );
  contentDiv->addStyleClass( "PeakFitDetPrefsExpandedBody" );

  // Detector type combo
  WLabel *detLabel = new WLabel( WString::tr( "pfdpg-det-type-label" ), contentDiv );
  detLabel->addStyleClass( "PfdpgLabel" );

  m_detTypeCombo = new WComboBox( contentDiv );
  m_detTypeCombo->addStyleClass( "PfdpgCombo" );
  m_detTypeCombo->addItem( WString::tr( "pfdpg-det-nai" ) );     // 0 => Low
  m_detTypeCombo->addItem( WString::tr( "pfdpg-det-labr" ) );    // 1 => LaBr
  m_detTypeCombo->addItem( WString::tr( "pfdpg-det-czt" ) );     // 2 => CZT
  m_detTypeCombo->addItem( WString::tr( "pfdpg-det-hpge" ) );    // 3 => High
  m_detTypeCombo->addItem( WString::tr( "pfdpg-det-unknown" ) ); // 4 => Unknown
  m_detTypeCombo->setCurrentIndex( 4 );
  m_detTypeCombo->activated().connect( std::bind( [this](){
    if( !m_programmaticUpdate )
      userChangedValue();
  }) );


  // FWHM method combo
  WLabel *fwhmLabel = new WLabel( WString::tr( "pfdpg-fwhm-method-label" ), contentDiv );
  fwhmLabel->addStyleClass( "PfdpgLabel" );

  m_fwhmMethodCombo = new WComboBox( contentDiv );
  m_fwhmMethodCombo->addStyleClass( "PfdpgCombo" );
  m_fwhmMethodCombo->addItem( WString::tr( "pfdpg-fwhm-normal" ) );      // 0 => Normal
  m_fwhmMethodCombo->addItem( WString::tr( "pfdpg-fwhm-det" ) );         // 1 => DetFwhm
  m_fwhmMethodCombo->addItem( WString::tr( "pfdpg-fwhm-det-refine" ) );  // 2 => DetPlusRefine
  m_fwhmMethodCombo->setCurrentIndex( 0 );
  m_fwhmMethodCombo->activated().connect( std::bind( [this](){
    if( m_programmaticUpdate )
      return;

    const int idx = m_fwhmMethodCombo->currentIndex();
    if( idx == 0 )
    {
      // Normal - just apply
      userChangedValue();
      return;
    }

    // DetFwhm or DetPlusRefine - check if DRF has FWHM info
    shared_ptr<const SpecMeas> meas = m_viewer
      ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground )
      : nullptr;
    shared_ptr<const DetectorPeakResponse> drf = meas ? meas->detector() : nullptr;

    if( drf && drf->hasResolutionInfo() )
    {
      userChangedValue();
      return;
    }

    // DRF lacks FWHM info - show dialog offering to fit FWHM or cancel
    SimpleDialog *dialog = new SimpleDialog(
      WString::tr( "pfdpg-no-fwhm-title" ),
      WString::tr( "pfdpg-no-fwhm-content" ) );

    WPushButton *fitBtn = dialog->addButton( WString::tr( "pfdpg-no-fwhm-fit" ) );
    WPushButton *cancelBtn = dialog->addButton( WString::tr( "pfdpg-no-fwhm-cancel" ) );

    const PeakFitDetPrefs::FwhmMethod pending = (idx == 1)
      ? PeakFitDetPrefs::FwhmMethod::DetFwhm
      : PeakFitDetPrefs::FwhmMethod::DetPlusRefine;

    cancelBtn->clicked().connect( std::bind( [this](){
      m_programmaticUpdate = true;
      m_fwhmMethodCombo->setCurrentIndex( 0 );
      m_programmaticUpdate = false;
    }) );

    fitBtn->clicked().connect( std::bind( [this, pending](){
      m_pendingFwhmMethod = pending;
      m_programmaticUpdate = true;
      m_fwhmMethodCombo->setCurrentIndex( 0 );
      m_programmaticUpdate = false;
      m_viewer->fwhmFromForegroundWindow( true );
    }) );
  }) );

  // When the DRF changes, check if we have a pending FWHM method to apply
  m_viewer->detectorModified().connect( std::bind( [this]( std::shared_ptr<DetectorPeakResponse> ){
    if( m_pendingFwhmMethod == PeakFitDetPrefs::FwhmMethod::Normal )
      return;

    shared_ptr<const SpecMeas> meas = m_viewer
      ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground )
      : nullptr;
    shared_ptr<const DetectorPeakResponse> drf = meas ? meas->detector() : nullptr;

    if( drf && drf->hasResolutionInfo() )
    {
      // FWHM is now available - apply the pending method
      m_programmaticUpdate = true;
      const int pendIdx = (m_pendingFwhmMethod == PeakFitDetPrefs::FwhmMethod::DetFwhm) ? 1 : 2;
      m_fwhmMethodCombo->setCurrentIndex( pendIdx );
      m_programmaticUpdate = false;
      m_pendingFwhmMethod = PeakFitDetPrefs::FwhmMethod::Normal;
      userChangedValue();
    }
    else
    {
      // Still no FWHM - clear pending
      m_pendingFwhmMethod = PeakFitDetPrefs::FwhmMethod::Normal;
    }
  }, std::placeholders::_1 ) );


  // Skew type combo
  WLabel *skewLabel = new WLabel( WString::tr( "pfdpg-skew-type-label" ), contentDiv );
  skewLabel->addStyleClass( "PfdpgLabel" );

  m_skewTypeCombo = new WComboBox( contentDiv );
  m_skewTypeCombo->addStyleClass( "PfdpgCombo" );
  for( int i = 0; i < static_cast<int>( PeakDef::NumSkewType ); ++i )
  {
    const PeakDef::SkewType st = static_cast<PeakDef::SkewType>( i );
    m_skewTypeCombo->addItem( WString::fromUTF8( PeakDef::to_label( st ) ) );
  }
  m_skewTypeCombo->setCurrentIndex( 0 );
  m_skewTypeCombo->activated().connect( std::bind( [this](){
    if( !m_programmaticUpdate )
    {
      updateSkewParamRows();
      userChangedValue();
    }
  }) );


  // ROI-independent skew checkbox
  WContainerWidget *cbRow = new WContainerWidget( contentDiv );
  cbRow->addStyleClass( "PfdpgCheckRow" );
  m_roiIndepCb = new WCheckBox( WString::tr( "pfdpg-roi-indep-skew" ), cbRow );
  m_roiIndepCb->setChecked( false );
  m_roiIndepCb->changed().connect( std::bind( [this](){
    if( !m_programmaticUpdate )
    {
      const bool hide = m_roiIndepCb->isChecked();
      m_skewParamsDiv->setHidden( hide );
      if( m_fitSkewLink )
        m_fitSkewLink->setHidden( hide );
      userChangedValue();
    }
  }) );

  // Skew parameters container (dynamic rows)
  m_skewParamsDiv = new WContainerWidget( contentDiv );
  m_skewParamsDiv->addStyleClass( "PfdpgSkewParams" );

  // "fit skew pars" link - visible when skew is selected and ROI-independent is unchecked
  m_fitSkewLink = new WText( WString::tr( "pfdpg-fit-skew-link" ), contentDiv );
  m_fitSkewLink->addStyleClass( "PfdpgFitSkewLink" );
  m_fitSkewLink->setHidden( true );
  m_fitSkewLink->clicked().connect( this, &PeakFitDetPrefsGui::showFitSkewDialog );

  // Source label
  m_sourceLabel = new WText( contentDiv );
  m_sourceLabel->addStyleClass( "PfdpgSourceLabel" );

  // Set initial state from spectrum (if any)
  updateFromSpecMeas();
}//void init()


void PeakFitDetPrefsGui::toggleExpanded()
{
  if( !m_compactMode || !m_collapsedDiv || !m_expandedDiv )
    return;

  const bool expanding = m_collapsedDiv->isVisible();
  m_collapsedDiv->setHidden( expanding );
  m_expandedDiv->setHidden( !expanding );
  
  if( expanding )
  {
    InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
    app->useMessageResourceBundle( "PeakEdit" ); //WMessageResourceBundle will check if its already been loaded, and if so, skip it
  }
}//void toggleExpanded()


namespace
{
  // Returns the PeakEdit.xml message IDs for label and tooltip of a given skew parameter.
  // The PeakEdit message resource bundle is already loaded by InterSpecApp.
  void skew_param_msg_ids( const PeakDef::SkewType skewType, const size_t paramIndex,
                           const char *&labelId, const char *&tooltipId )
  {
    labelId = nullptr;
    tooltipId = nullptr;

    switch( skewType )
    {
      case PeakDef::NoSkew:
      case PeakDef::NumSkewType:
        return;

      case PeakDef::Bortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-bortel-tau"; tooltipId = "pe-tt-skew-bortel-tau"; }
        return;

      case PeakDef::GaussExp:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-gaussexp-k"; tooltipId = "pe-tt-skew-gaussexp-k"; }
        return;

      case PeakDef::CrystalBall:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-crystalball-alpha"; tooltipId = "pe-tt-skew-crystalball-alpha"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-crystalball-n"; tooltipId = "pe-tt-skew-crystalball-n"; }
        return;

      case PeakDef::ExpGaussExp:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-expgaussexp-kl"; tooltipId = "pe-tt-skew-expgaussexp-kl"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-expgaussexp-kh"; tooltipId = "pe-tt-skew-expgaussexp-kh"; }
        return;

      case PeakDef::DoubleSidedCrystalBall:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-dscb-alphalow"; tooltipId = "pe-tt-skew-dscb-alphalow"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-dscb-nlow"; tooltipId = "pe-tt-skew-dscb-nlow"; }
        if( paramIndex == 2 ){ labelId = "pe-label-skew-dscb-alphahigh"; tooltipId = "pe-tt-skew-dscb-alphahigh"; }
        if( paramIndex == 3 ){ labelId = "pe-label-skew-dscb-nhigh"; tooltipId = "pe-tt-skew-dscb-nhigh"; }
        return;

      case PeakDef::VoigtPlusBortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-voigtplusbortel-gamma"; tooltipId = "pe-tt-skew-voigtplusbortel-gamma"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-voigtplusbortel-r"; tooltipId = "pe-tt-skew-voigtplusbortel-r"; }
        if( paramIndex == 2 ){ labelId = "pe-label-skew-voigtplusbortel-tau"; tooltipId = "pe-tt-skew-voigtplusbortel-tau"; }
        return;

      case PeakDef::GaussPlusBortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-gaussplusbortel-r"; tooltipId = "pe-tt-skew-gaussplusbortel-r"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-gaussplusbortel-tau"; tooltipId = "pe-tt-skew-gaussplusbortel-tau"; }
        return;

      case PeakDef::DoubleBortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-doublebortel-tau1"; tooltipId = "pe-tt-skew-doublebortel-tau1"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-doublebortel-deltatau2"; tooltipId = "pe-tt-skew-doublebortel-deltatau2"; }
        if( paramIndex == 2 ){ labelId = "pe-label-skew-doublebortel-eta"; tooltipId = "pe-tt-skew-doublebortel-eta"; }
        return;
    }//switch( skewType )
  }//void skew_param_msg_ids(...)


  /** Returns reasonable default skew parameter values for a given skew type and detector type.
   *  For scintillators (NaI/LaBr), values are chosen to produce near-Gaussian peaks.
   *  For HPGe, values reflect typical low-energy tailing that is stronger at lower energies.
   *  For CZT, values are intermediate.
   *  For Unknown, the starting values from PeakDef::skew_parameter_range are used.
   *
   *  The lower/upper arrays correspond to the low-energy/high-energy ends of the spectrum
   *  for energy-dependent parameters. For non-energy-dependent parameters, only lower[i] is used.
   */
  void default_skew_values( const PeakDef::SkewType skewType,
                             const PeakFitUtils::CoarseResolutionType detType,
                             std::optional<double> lower[4],
                             std::optional<double> upper[4] )
  {
    for( int i = 0; i < 4; ++i )
    {
      lower[i] = std::nullopt;
      upper[i] = std::nullopt;
    }

    const size_t nparams = PeakDef::num_skew_parameters( skewType );
    if( nparams == 0 )
      return;

    // For Unknown detector, use the starting values from skew_parameter_range
    if( detType == PeakFitUtils::CoarseResolutionType::Unknown )
    {
      for( size_t p = 0; p < nparams; ++p )
      {
        const PeakDef::CoefficientType coefType
          = static_cast<PeakDef::CoefficientType>(
            static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ) );
        double lo = 0, hi = 0, start = 0, step = 0;
        PeakDef::skew_parameter_range( skewType, coefType, lo, hi, start, step );
        lower[p] = start;
        if( PeakDef::is_energy_dependent( skewType, coefType ) )
          upper[p] = start;
      }
      return;
    }//if( Unknown )

    // Naming: L=Low(NaI), B=LaBr, C=CZT, H=HPGe
    // For energy-dep params, {low_E_val, high_E_val}; for non-energy-dep, just {val}
    // Key: larger tau/k/alpha = less skew; R/eta near 0 = more Gaussian

    switch( skewType )
    {
      case PeakDef::NoSkew:
      case PeakDef::NumSkewType:
        return;

      case PeakDef::Bortel: // tau: range [0,15], larger=less tail
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:  lower[0] = 8.0; upper[0] = 8.0; break;
          case PeakFitUtils::CoarseResolutionType::LaBr: lower[0] = 7.0; upper[0] = 7.0; break;
          case PeakFitUtils::CoarseResolutionType::CZT:  lower[0] = 3.0; upper[0] = 5.0; break;
          case PeakFitUtils::CoarseResolutionType::High: lower[0] = 0.8; upper[0] = 3.0; break;
          default: break;
        }
        return;

      case PeakDef::GaussExp: // k: range [0.15,3.25], larger=less tail
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:  lower[0] = 3.0;  upper[0] = 3.0;  break;
          case PeakFitUtils::CoarseResolutionType::LaBr: lower[0] = 3.0;  upper[0] = 3.0;  break;
          case PeakFitUtils::CoarseResolutionType::CZT:  lower[0] = 2.0;  upper[0] = 2.5;  break;
          case PeakFitUtils::CoarseResolutionType::High: lower[0] = 0.5;  upper[0] = 2.0;  break;
          default: break;
        }
        return;

      case PeakDef::CrystalBall: // alpha (energy-dep): [0.5,4], larger=less tail; n (not): [1.05,100]
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:  lower[0] = 3.5; upper[0] = 3.5; lower[1] = 5.0; break;
          case PeakFitUtils::CoarseResolutionType::LaBr: lower[0] = 3.5; upper[0] = 3.5; lower[1] = 5.0; break;
          case PeakFitUtils::CoarseResolutionType::CZT:  lower[0] = 2.5; upper[0] = 3.0; lower[1] = 3.0; break;
          case PeakFitUtils::CoarseResolutionType::High: lower[0] = 1.2; upper[0] = 2.5; lower[1] = 3.0; break;
          default: break;
        }
        return;

      case PeakDef::ExpGaussExp: // k_low, k_high: both energy-dep, [0.15,3.25]
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:
            lower[0] = 3.0; upper[0] = 3.0; lower[1] = 3.0; upper[1] = 3.0; break;
          case PeakFitUtils::CoarseResolutionType::LaBr:
            lower[0] = 3.0; upper[0] = 3.0; lower[1] = 3.0; upper[1] = 3.0; break;
          case PeakFitUtils::CoarseResolutionType::CZT:
            lower[0] = 2.0; upper[0] = 2.5; lower[1] = 3.0; upper[1] = 3.0; break;
          case PeakFitUtils::CoarseResolutionType::High:
            lower[0] = 0.5; upper[0] = 2.0; lower[1] = 3.0; upper[1] = 3.0; break;
          default: break;
        }
        return;

      case PeakDef::DoubleSidedCrystalBall: // a_lo(ed), n_lo, a_hi(ed), n_hi
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:
            lower[0] = 3.5; upper[0] = 3.5; lower[1] = 5.0;
            lower[2] = 3.5; upper[2] = 3.5; lower[3] = 5.0; break;
          case PeakFitUtils::CoarseResolutionType::LaBr:
            lower[0] = 3.5; upper[0] = 3.5; lower[1] = 5.0;
            lower[2] = 3.5; upper[2] = 3.5; lower[3] = 5.0; break;
          case PeakFitUtils::CoarseResolutionType::CZT:
            lower[0] = 2.5; upper[0] = 3.0; lower[1] = 3.0;
            lower[2] = 3.5; upper[2] = 3.5; lower[3] = 5.0; break;
          case PeakFitUtils::CoarseResolutionType::High:
            lower[0] = 1.2; upper[0] = 2.5; lower[1] = 3.0;
            lower[2] = 3.5; upper[2] = 3.5; lower[3] = 5.0; break;
          default: break;
        }
        return;

      case PeakDef::VoigtPlusBortel: // gamma(ed), R(ed), tau(ed) - all energy-dep
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:
            lower[0] = 0.01; upper[0] = 0.01;
            lower[1] = 0.01; upper[1] = 0.01;
            lower[2] = 5.0;  upper[2] = 5.0;  break;
          case PeakFitUtils::CoarseResolutionType::LaBr:
            lower[0] = 0.01; upper[0] = 0.01;
            lower[1] = 0.01; upper[1] = 0.01;
            lower[2] = 5.0;  upper[2] = 5.0;  break;
          case PeakFitUtils::CoarseResolutionType::CZT:
            lower[0] = 0.1;  upper[0] = 0.05;
            lower[1] = 0.1;  upper[1] = 0.05;
            lower[2] = 2.0;  upper[2] = 3.0;  break;
          case PeakFitUtils::CoarseResolutionType::High:
            lower[0] = 0.2;  upper[0] = 0.05;
            lower[1] = 0.3;  upper[1] = 0.1;
            lower[2] = 0.8;  upper[2] = 2.5;  break;
          default: break;
        }
        return;

      case PeakDef::GaussPlusBortel: // R(ed), tau(ed)
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:
            lower[0] = 0.01; upper[0] = 0.01;
            lower[1] = 5.0;  upper[1] = 5.0;  break;
          case PeakFitUtils::CoarseResolutionType::LaBr:
            lower[0] = 0.01; upper[0] = 0.01;
            lower[1] = 5.0;  upper[1] = 5.0;  break;
          case PeakFitUtils::CoarseResolutionType::CZT:
            lower[0] = 0.15; upper[0] = 0.05;
            lower[1] = 2.0;  upper[1] = 3.0;  break;
          case PeakFitUtils::CoarseResolutionType::High:
            lower[0] = 0.3;  upper[0] = 0.1;
            lower[1] = 0.8;  upper[1] = 2.5;  break;
          default: break;
        }
        return;

      case PeakDef::DoubleBortel: // tau1(ed), tau2_delta(ed), eta(not)
        switch( detType )
        {
          case PeakFitUtils::CoarseResolutionType::Low:
            lower[0] = 5.0; upper[0] = 5.0;
            lower[1] = 1.0; upper[1] = 1.0;
            lower[2] = 0.1; break;
          case PeakFitUtils::CoarseResolutionType::LaBr:
            lower[0] = 5.0; upper[0] = 5.0;
            lower[1] = 1.0; upper[1] = 1.0;
            lower[2] = 0.1; break;
          case PeakFitUtils::CoarseResolutionType::CZT:
            lower[0] = 2.0; upper[0] = 3.0;
            lower[1] = 1.0; upper[1] = 0.5;
            lower[2] = 0.3; break;
          case PeakFitUtils::CoarseResolutionType::High:
            lower[0] = 0.5; upper[0] = 2.0;
            lower[1] = 2.0; upper[1] = 1.0;
            lower[2] = 0.5; break;
          default: break;
        }
        return;
    }//switch( skewType )
  }//void default_skew_values(...)
}//anonymous namespace


void PeakFitDetPrefsGui::updateSkewParamRows()
{
  // Clear existing param widgets
  m_skewParamsDiv->clear();
  for( int i = 0; i < 4; ++i )
  {
    m_lowerSkewSpin[i] = nullptr;
    m_upperSkewSpin[i] = nullptr;
  }

  const int skewIdx = m_skewTypeCombo->currentIndex();
  if( skewIdx < 0 || skewIdx >= static_cast<int>( PeakDef::NumSkewType ) )
    return;

  const PeakDef::SkewType skewType = static_cast<PeakDef::SkewType>( skewIdx );
  const size_t nparams = PeakDef::num_skew_parameters( skewType );
  const bool hasSkew = (nparams > 0);

  // Hide checkbox, skew params, and fit link when NoSkew is selected
  if( m_roiIndepCb && m_roiIndepCb->parent() )
    m_roiIndepCb->parent()->setHidden( !hasSkew );
  m_skewParamsDiv->setHidden( !hasSkew || m_roiIndepCb->isChecked() );
  if( m_fitSkewLink )
    m_fitSkewLink->setHidden( !hasSkew || m_roiIndepCb->isChecked() );

  if( !hasSkew )
    return;

  // Check if any parameter is energy-dependent (needs two columns)
  bool hasEnergyDep = false;
  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType coefType
      = static_cast<PeakDef::CoefficientType>(
          static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ) );
    if( PeakDef::is_energy_dependent( skewType, coefType ) )
    {
      hasEnergyDep = true;
      break;
    }
  }//for( check energy dependence )

  // Use a single grid table for all params: 3-col if energy-dep, 2-col otherwise
  WContainerWidget *table = new WContainerWidget( m_skewParamsDiv );
  table->addStyleClass( hasEnergyDep ? "PfdpgParamTable" : "PfdpgParamTable PfdpgSingleCol" );

  // Single header row at top with "Low E" / "High E" column labels
  if( hasEnergyDep )
  {
    new WText( "", table ); // empty name-column cell
    WText *lowHeader = new WText( WString::tr( "pfdpg-lower-header" ), table );
    lowHeader->addStyleClass( "PfdpgColHeader" );
    WText *highHeader = new WText( WString::tr( "pfdpg-upper-header" ), table );
    highHeader->addStyleClass( "PfdpgColHeader" );
  }

  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType coefType
      = static_cast<PeakDef::CoefficientType>(
          static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ) );

    const bool energyDep = PeakDef::is_energy_dependent( skewType, coefType );

    double range_lower = 0, range_upper = 0, start_val = 0, step_size = 0;
    PeakDef::skew_parameter_range( skewType, coefType, range_lower, range_upper, start_val, step_size );

    const char *labelMsgId = nullptr;
    const char *tooltipMsgId = nullptr;
    skew_param_msg_ids( skewType, p, labelMsgId, tooltipMsgId );

    // Build tooltip text: base tooltip + energy-dependence note
    WString tooltipText;
    if( tooltipMsgId )
    {
      tooltipText = WString::tr( tooltipMsgId );
      tooltipText += energyDep ? WString::tr( "pfdpg-tt-energy-dep" )
                               : WString::tr( "pfdpg-tt-not-energy-dep" );
    }

    // Parameter label
    WText *paramLabel = new WText(
      labelMsgId ? WString::tr( labelMsgId ) : WString::tr( "pfdpg-param-label" ).arg( static_cast<int>( p ) ),
      table );
    paramLabel->addStyleClass( "PfdpgParamName" );

    if( !tooltipText.empty() )
    {
      HelpSystem::attachToolTipOn( paramLabel, tooltipText, true,
                                   HelpSystem::ToolTipPosition::Right );
    }

    if( energyDep )
    {
      // Energy-dependent: label | lowerSpin | upperSpin on one row
      NativeFloatSpinBox *lowerSpin = new NativeFloatSpinBox( table );
      lowerSpin->setRange( static_cast<float>( range_lower ), static_cast<float>( range_upper ) );
      //lowerSpin->setSingleStep( static_cast<float>( step_size ) );
      lowerSpin->setFormatString( "%.4G" );
      lowerSpin->setSpinnerHidden( true );
      lowerSpin->addStyleClass( "PfdpgSpin" );
      lowerSpin->setPlaceholderText( "fit" );
      m_lowerSkewSpin[p] = lowerSpin;
      lowerSpin->valueChanged().connect( std::bind( [this]( float ){
        if( !m_programmaticUpdate )
          userChangedValue();
      }, std::placeholders::_1 ) );

      NativeFloatSpinBox *upperSpin = new NativeFloatSpinBox( table );
      upperSpin->setRange( static_cast<float>( range_lower ), static_cast<float>( range_upper ) );
      //upperSpin->setSingleStep( static_cast<float>( step_size ) );
      upperSpin->setFormatString( "%.4G" );
      upperSpin->setSpinnerHidden( true );
      upperSpin->addStyleClass( "PfdpgSpin" );
      upperSpin->setPlaceholderText( "fit" );
      
      m_upperSkewSpin[p] = upperSpin;
      upperSpin->valueChanged().connect( std::bind( [this]( float ){
        if( !m_programmaticUpdate )
          userChangedValue();
      }, std::placeholders::_1 ) );

      if( !m_programmaticUpdate )
      {
        lowerSpin->setValue( static_cast<float>(start_val) );
        upperSpin->setValue( static_cast<float>(start_val) );
      }
      
      if( !tooltipText.empty() )
      {
        HelpSystem::attachToolTipOn( lowerSpin, tooltipText, true,
                                     HelpSystem::ToolTipPosition::Right );
        HelpSystem::attachToolTipOn( upperSpin, tooltipText, true,
                                     HelpSystem::ToolTipPosition::Right );
      }
    }
    else
    {
      // Non-energy-dependent: label | valSpin (spanning remaining columns if 3-col grid)
      NativeFloatSpinBox *valSpin = new NativeFloatSpinBox( table );
      valSpin->setRange( static_cast<float>( range_lower ), static_cast<float>( range_upper ) );
      //valSpin->setSingleStep( static_cast<float>( step_size ) );
      valSpin->setFormatString( "%.4G" );
      valSpin->setSpinnerHidden( true );
      valSpin->addStyleClass( hasEnergyDep ? "PfdpgSpin PfdpgSpinWide" : "PfdpgSpin" );
      valSpin->setPlaceholderText( "fit" );
      
      m_lowerSkewSpin[p] = valSpin;
      valSpin->valueChanged().connect( std::bind( [this]( float ){
        if( !m_programmaticUpdate )
          userChangedValue();
      }, std::placeholders::_1 ) );

      if( !m_programmaticUpdate )
        valSpin->setValue( static_cast<float>(start_val) );
      
      if( !tooltipText.empty() )
      {
        HelpSystem::attachToolTipOn( valSpin, tooltipText, true,
                                     HelpSystem::ToolTipPosition::Right );
      }
    }
  }//for( each skew param )

  // Fill in reasonable defaults when user changes skew type (not during programmatic update)
  if( !m_programmaticUpdate )
  {
    /*
    const int detIdx = m_detTypeCombo->currentIndex();
    PeakFitUtils::CoarseResolutionType detType = PeakFitUtils::CoarseResolutionType::Unknown;
    switch( detIdx )
    {
      case 0:  detType = PeakFitUtils::CoarseResolutionType::Low;  break;
      case 1:  detType = PeakFitUtils::CoarseResolutionType::LaBr; break;
      case 2:  detType = PeakFitUtils::CoarseResolutionType::CZT;  break;
      case 3:  detType = PeakFitUtils::CoarseResolutionType::High; break;
      default: break;
    }

    std::optional<double> defLower[4], defUpper[4];
    default_skew_values( skewType, detType, defLower, defUpper );

    for( size_t p = 0; p < nparams; ++p )
    {
      if( m_lowerSkewSpin[p] && defLower[p].has_value() )
        m_lowerSkewSpin[p]->setValue( static_cast<float>( defLower[p].value() ) );
      if( m_upperSkewSpin[p] && defUpper[p].has_value() )
        m_upperSkewSpin[p]->setValue( static_cast<float>( defUpper[p].value() ) );
    }
     */
  }//if( !m_programmaticUpdate )
}//void updateSkewParamRows()


void PeakFitDetPrefsGui::userChangedValue()
{
  if( m_programmaticUpdate )
    return;

  shared_ptr<SpecMeas> meas = m_viewer
    ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground )
    : nullptr;
  if( !meas )
    return;

  // Capture old prefs for undo
  shared_ptr<const PeakFitDetPrefs> oldPrefs = meas->peakFitDetPrefs();

  // Build new prefs from GUI
  shared_ptr<PeakFitDetPrefs> newPrefs = make_shared<PeakFitDetPrefs>();

  // Detector type: combo index maps to CoarseResolutionType
  const int detIdx = m_detTypeCombo->currentIndex();
  switch( detIdx )
  {
    case 0:  newPrefs->m_det_type = PeakFitUtils::CoarseResolutionType::Low;     break;
    case 1:  newPrefs->m_det_type = PeakFitUtils::CoarseResolutionType::LaBr;    break;
    case 2:  newPrefs->m_det_type = PeakFitUtils::CoarseResolutionType::CZT;     break;
    case 3:  newPrefs->m_det_type = PeakFitUtils::CoarseResolutionType::High;    break;
    default: newPrefs->m_det_type = PeakFitUtils::CoarseResolutionType::Unknown; break;
  }

  // Skew type
  const int skewIdx = m_skewTypeCombo->currentIndex();
  if( skewIdx >= 0 && skewIdx < static_cast<int>( PeakDef::NumSkewType ) )
    newPrefs->m_peak_skew_type = static_cast<PeakDef::SkewType>( skewIdx );
  else
    newPrefs->m_peak_skew_type = PeakDef::NoSkew;

  // Skew params: read from spin boxes if they exist and have valid text
  const size_t nparams = PeakDef::num_skew_parameters( newPrefs->m_peak_skew_type );
  for( size_t p = 0; p < nparams; ++p )
  {
    if( m_lowerSkewSpin[p] && !m_lowerSkewSpin[p]->text().empty()
       && m_lowerSkewSpin[p]->text().toUTF8() != "fit" )
    {
      newPrefs->m_lower_energy_skew[p] = static_cast<double>( m_lowerSkewSpin[p]->value() );
    }

    if( m_upperSkewSpin[p] && !m_upperSkewSpin[p]->text().empty()
       && m_upperSkewSpin[p]->text().toUTF8() != "fit" )
    {
      newPrefs->m_upper_energy_skew[p] = static_cast<double>( m_upperSkewSpin[p]->value() );
    }
  }

  newPrefs->m_roi_independent_skew = m_roiIndepCb->isChecked();

  // FWHM method
  const int fwhmIdx = m_fwhmMethodCombo->currentIndex();
  switch( fwhmIdx )
  {
    case 1:  newPrefs->m_fwhm_method = PeakFitDetPrefs::FwhmMethod::DetFwhm;       break;
    case 2:  newPrefs->m_fwhm_method = PeakFitDetPrefs::FwhmMethod::DetPlusRefine;  break;
    default: newPrefs->m_fwhm_method = PeakFitDetPrefs::FwhmMethod::Normal;         break;
  }

  newPrefs->m_source = PeakFitDetPrefs::LoadingSource::UserInputInGui;

  // Apply to SpecMeas
  meas->setPeakFitDetPrefs( newPrefs );

  // Update the source label
  m_sourceLabel->setText( WString::tr( "pfdpg-source-label" )
                         + WString::tr( "pfdpg-src-user" ) );

  // Register undo/redo
  UndoRedoManager *undoManager = m_viewer->undoRedoManager();
  if( undoManager && undoManager->canAddUndoRedoNow() )
  {
    weak_ptr<SpecMeas> weakMeas = meas;

    auto undo = [oldPrefs, weakMeas](){
      shared_ptr<SpecMeas> m = weakMeas.lock();
      if( m )
        m->setPeakFitDetPrefs( oldPrefs );
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        viewer->peakFitDetPrefsChanged().emit();
    };

    auto redo = [newPrefs, weakMeas](){
      shared_ptr<SpecMeas> m = weakMeas.lock();
      if( m )
        m->setPeakFitDetPrefs( newPrefs );
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        viewer->peakFitDetPrefsChanged().emit();
    };

    undoManager->addUndoRedoStep( undo, redo,
      WString::tr( "pfdpg-undo-det-type" ).toUTF8() );
  }//if( undoManager )
}//void userChangedValue()


void PeakFitDetPrefsGui::handlePrefsChanged()
{
  updateFromSpecMeas();
}//void handlePrefsChanged()


void PeakFitDetPrefsGui::updateFromSpecMeas()
{
  m_programmaticUpdate = true;

  shared_ptr<const SpecMeas> meas = m_viewer
    ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground )
    : nullptr;

  if( !meas )
  {
    setControlsEnabled( false );
    m_detTypeCombo->setCurrentIndex( 4 ); // Unknown
    m_fwhmMethodCombo->setCurrentIndex( 0 ); // Normal
    m_skewTypeCombo->setCurrentIndex( 0 ); // NoSkew
    updateSkewParamRows();
    m_sourceLabel->setText( "" );
    m_programmaticUpdate = false;
    return;
  }

  setControlsEnabled( true );

  shared_ptr<const PeakFitDetPrefs> prefs = meas->peakFitDetPrefs();

  if( !prefs )
  {
    m_detTypeCombo->setCurrentIndex( 4 );
    m_fwhmMethodCombo->setCurrentIndex( 0 );
    m_skewTypeCombo->setCurrentIndex( 0 );
    updateSkewParamRows();
    m_sourceLabel->setText( WString::tr( "pfdpg-source-label" )
                           + WString::tr( "pfdpg-src-default" ) );
    m_programmaticUpdate = false;
    return;
  }

  // Set detector type combo
  switch( prefs->m_det_type )
  {
    case PeakFitUtils::CoarseResolutionType::Low:     m_detTypeCombo->setCurrentIndex( 0 ); break;
    case PeakFitUtils::CoarseResolutionType::LaBr:    m_detTypeCombo->setCurrentIndex( 1 ); break;
    case PeakFitUtils::CoarseResolutionType::CZT:     m_detTypeCombo->setCurrentIndex( 2 ); break;
    case PeakFitUtils::CoarseResolutionType::High:    m_detTypeCombo->setCurrentIndex( 3 ); break;
    case PeakFitUtils::CoarseResolutionType::Unknown: m_detTypeCombo->setCurrentIndex( 4 ); break;
  }//switch( prefs->m_det_type )

  // Set FWHM method combo
  // If DetFwhm/DetPlusRefine but DRF lacks FWHM info, force Normal
  {
    PeakFitDetPrefs::FwhmMethod effective = prefs->m_fwhm_method;
    if( effective != PeakFitDetPrefs::FwhmMethod::Normal )
    {
      shared_ptr<const DetectorPeakResponse> drf = meas->detector();
      if( !drf || !drf->hasResolutionInfo() )
        effective = PeakFitDetPrefs::FwhmMethod::Normal;
    }

    switch( effective )
    {
      case PeakFitDetPrefs::FwhmMethod::Normal:        m_fwhmMethodCombo->setCurrentIndex( 0 ); break;
      case PeakFitDetPrefs::FwhmMethod::DetFwhm:       m_fwhmMethodCombo->setCurrentIndex( 1 ); break;
      case PeakFitDetPrefs::FwhmMethod::DetPlusRefine:  m_fwhmMethodCombo->setCurrentIndex( 2 ); break;
    }
  }

  // Set skew type combo
  const int skewIdx = static_cast<int>( prefs->m_peak_skew_type );
  if( skewIdx >= 0 && skewIdx < static_cast<int>( PeakDef::NumSkewType ) )
    m_skewTypeCombo->setCurrentIndex( skewIdx );
  else
    m_skewTypeCombo->setCurrentIndex( 0 );

  // Set ROI-independent skew checkbox
  m_roiIndepCb->setChecked( prefs->m_roi_independent_skew );

  // Rebuild skew param rows and update checkbox/params visibility
  updateSkewParamRows();

  // Fill in skew param values
  const size_t nparams = PeakDef::num_skew_parameters( prefs->m_peak_skew_type );
  for( size_t p = 0; p < nparams; ++p )
  {
    if( m_lowerSkewSpin[p] )
    {
      if( prefs->m_lower_energy_skew[p].has_value() )
        m_lowerSkewSpin[p]->setValue( static_cast<float>( prefs->m_lower_energy_skew[p].value() ) );
      else
        m_lowerSkewSpin[p]->setText( "" );
    }

    if( m_upperSkewSpin[p] )
    {
      if( prefs->m_upper_energy_skew[p].has_value() )
        m_upperSkewSpin[p]->setValue( static_cast<float>( prefs->m_upper_energy_skew[p].value() ) );
      else
        m_upperSkewSpin[p]->setText( "" );
    }
  }//for( each param )

  // Update source label
  WString srcTxt;
  switch( prefs->m_source )
  {
    case PeakFitDetPrefs::LoadingSource::Default:                 srcTxt = WString::tr( "pfdpg-src-default" );   break;
    case PeakFitDetPrefs::LoadingSource::UserInputInGui:          srcTxt = WString::tr( "pfdpg-src-user" );      break;
    case PeakFitDetPrefs::LoadingSource::FromDetectorPeakResponse: srcTxt = WString::tr( "pfdpg-src-drf" );      break;
    case PeakFitDetPrefs::LoadingSource::DefaultForDetectorType:  srcTxt = WString::tr( "pfdpg-src-det-type" );  break;
    case PeakFitDetPrefs::LoadingSource::FromSpectralData:        srcTxt = WString::tr( "pfdpg-src-spectral" );  break;
  }//switch
  m_sourceLabel->setText( WString::tr( "pfdpg-source-label" ) + srcTxt );

  m_programmaticUpdate = false;
}//void updateFromSpecMeas()


void PeakFitDetPrefsGui::setControlsEnabled( const bool enabled )
{
  m_detTypeCombo->setEnabled( enabled );
  m_fwhmMethodCombo->setEnabled( enabled );
  m_skewTypeCombo->setEnabled( enabled );
  m_roiIndepCb->setEnabled( enabled );

  for( int i = 0; i < 4; ++i )
  {
    if( m_lowerSkewSpin[i] )
      m_lowerSkewSpin[i]->setEnabled( enabled );
    if( m_upperSkewSpin[i] )
      m_upperSkewSpin[i]->setEnabled( enabled );
  }
}//void setControlsEnabled( bool )


shared_ptr<const PeakFitDetPrefs> PeakFitDetPrefsGui::currentPrefs() const
{
  shared_ptr<const SpecMeas> meas = m_viewer
    ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground )
    : nullptr;

  return meas ? meas->peakFitDetPrefs() : nullptr;
}//shared_ptr<const PeakFitDetPrefs> currentPrefs() const


void PeakFitDetPrefsGui::showFitSkewDialog()
{
  if( m_viewer )
    m_viewer->showFitSkewParamsWindow();
}//void showFitSkewDialog()
