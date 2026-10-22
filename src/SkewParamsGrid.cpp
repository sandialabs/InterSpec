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

#include <cassert>

#include <Wt/WText.h>
#include <Wt/WCheckBox.h>
#include <Wt/WApplication.h>

#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SkewParamsGrid.h"
#include "InterSpec/NativeFloatSpinBox.h"

using namespace std;
using namespace Wt;


namespace
{
  // The PeakEdit.xml message IDs for the label and tooltip of a skew parameter.
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

      case PeakDef::GadrasGeneric:
      case PeakDef::GadrasCZT:
      {
        static const char * const s_gadras_label_ids[6] = {
          "pe-label-skew-gadras-0", "pe-label-skew-gadras-1", "pe-label-skew-gadras-2",
          "pe-label-skew-gadras-3", "pe-label-skew-gadras-4", "pe-label-skew-gadras-5"
        };
        static const char * const s_gadras_tt_ids[6] = {
          "pe-tt-skew-gadras-0", "pe-tt-skew-gadras-1", "pe-tt-skew-gadras-2",
          "pe-tt-skew-gadras-3", "pe-tt-skew-gadras-4", "pe-tt-skew-gadras-5"
        };
        if( paramIndex < 6 ){ labelId = s_gadras_label_ids[paramIndex]; tooltipId = s_gadras_tt_ids[paramIndex]; }
        return;
      }
    }//switch( skewType )
  }//void skew_param_msg_ids(...)
}//namespace


SkewParamsGrid::SkewParamsGrid( const bool show_fit_column, const bool allow_blank )
  : WContainerWidget(),
    m_show_fit( show_fit_column ),
    m_allow_blank( allow_blank ),
    m_editable( true ),
    m_skew_type( PeakDef::SkewType::NoSkew ),
    m_result_header( nullptr )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
  {
    app->useMessageResourceBundle( "SkewParamsGrid" );
    app->useMessageResourceBundle( "PeakEdit" );
  }
  wApp->useStyleSheet( "InterSpec_resources/SkewParamsGrid.css" );

  addStyleClass( "SkewParamsGrid" );
  if( m_show_fit )
    addStyleClass( "SpgFit" );

  setSkewType( PeakDef::SkewType::NoSkew );
}//SkewParamsGrid constructor


void SkewParamsGrid::setSkewType( const PeakDef::SkewType type )
{
  m_skew_type = type;

  clear();
  for( size_t i = 0; i < 6; ++i )
  {
    m_lower[i] = m_upper[i] = nullptr;
    m_fit[i] = nullptr;
    m_result[i] = nullptr;
  }
  m_result_header = nullptr;
  removeStyleClass( "SpgV1" );
  removeStyleClass( "SpgV2" );
  removeStyleClass( "SpgRes" );

  const size_t nparams = PeakDef::num_skew_parameters( type );
  setHidden( nparams == 0 );
  if( nparams == 0 )
    return;

  bool has_energy_dep = false;
  for( size_t p = 0; p < nparams; ++p )
    has_energy_dep |= PeakDef::is_energy_dependent( type, PeakDef::CoefficientType( PeakDef::SkewPar0 + p ) );
  addStyleClass( has_energy_dep ? "SpgV2" : "SpgV1" );

  // Header row - only needed to label the energy columns, or the fit column
  if( has_energy_dep || m_show_fit )
  {
    addNew<WText>();
    if( has_energy_dep )
    {
      addNew<WText>( WString::tr( "spg-lower-header" ) )->addStyleClass( "SpgColHeader" );
      addNew<WText>( WString::tr( "spg-upper-header" ) )->addStyleClass( "SpgColHeader" );
    }else
    {
      addNew<WText>();
    }

    if( m_show_fit )
    {
      WText *fit_header = addNew<WText>( WString::tr( "spg-fit-header" ) );
      fit_header->addStyleClass( "SpgColHeader" );
      HelpSystem::attachToolTipOn( fit_header, WString::tr( "spg-tt-fit" ), true );

      m_result_header = addNew<WText>( WString::tr( "spg-result-header" ) );
      m_result_header->addStyleClass( "SpgColHeader" );
      m_result_header->setHidden( true );
    }
  }//if( need a header row )

  const auto add_spin = [this]( const double range_lower, const double range_upper, const bool wide ) -> NativeFloatSpinBox * {
    NativeFloatSpinBox *spin = addNew<NativeFloatSpinBox>();
    spin->setRange( static_cast<float>( range_lower ), static_cast<float>( range_upper ) );
    spin->setFormatString( "%.4G" );
    spin->setSpinnerHidden( true );
    spin->addStyleClass( wide ? "SpgSpin SpgSpinWide" : "SpgSpin" );
    if( m_allow_blank )
      spin->setPlaceholderText( WString::tr( "spg-blank-placeholder" ) );
    spin->setEnabled( m_editable );
    spin->valueChanged().connect( [this]( float ){ m_user_changed.emit(); } );
    return spin;
  };//add_spin lambda

  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType ct = PeakDef::CoefficientType( PeakDef::SkewPar0 + p );
    const bool energy_dep = PeakDef::is_energy_dependent( type, ct );

    double range_lower = 0.0, range_upper = 0.0, start_val = 0.0, step_size = 0.0;
    PeakDef::skew_parameter_range( type, ct, range_lower, range_upper, start_val, step_size );

    const char *label_id = nullptr, *tooltip_id = nullptr;
    skew_param_msg_ids( type, p, label_id, tooltip_id );

    WString tooltip;
    if( tooltip_id )
      tooltip = WString::tr( tooltip_id )
                + WString::tr( energy_dep ? "spg-tt-energy-dep" : "spg-tt-not-energy-dep" );

    WText *label = addNew<WText>( label_id ? WString::tr( label_id )
                                           : WString::tr( "spg-param-label" ).arg( static_cast<int>(p) ) );
    label->addStyleClass( "SpgParamName" );

    m_lower[p] = add_spin( range_lower, range_upper, has_energy_dep && !energy_dep );
    m_lower[p]->setValue( static_cast<float>( start_val ) );
    if( energy_dep )
    {
      m_upper[p] = add_spin( range_lower, range_upper, false );
      m_upper[p]->setValue( static_cast<float>( start_val ) );
    }

    if( !tooltip.empty() )
    {
      HelpSystem::attachToolTipOn( label, tooltip, true );
      HelpSystem::attachToolTipOn( m_lower[p], tooltip, true );
      if( m_upper[p] )
        HelpSystem::attachToolTipOn( m_upper[p], tooltip, true );
    }

    if( m_show_fit )
    {
      m_fit[p] = addNew<WCheckBox>();
      m_fit[p]->addStyleClass( "SpgFitCb" );
      m_fit[p]->setChecked( PeakDef::skew_parameter_fit_by_default( type, ct ) );
      m_fit[p]->setEnabled( m_editable );
      // Not `changed()`: on Wt 4.13 its value can be overwritten by a stale one from a pointer event
      //  batched with it, before the render that reads it.
      m_fit[p]->checked().connect( [this](){ m_user_changed.emit(); } );
      m_fit[p]->unChecked().connect( [this](){ m_user_changed.emit(); } );

      m_result[p] = addNew<WText>();
      m_result[p]->addStyleClass( "SpgResult" );
      m_result[p]->setHidden( true );
    }//if( m_show_fit )
  }//for( size_t p = 0; p < nparams; ++p )
}//void setSkewType( const PeakDef::SkewType type )


PeakDef::SkewType SkewParamsGrid::skewType() const
{
  return m_skew_type;
}


void SkewParamsGrid::setValue( const size_t index, const std::optional<double> &lower,
                               const std::optional<double> &upper )
{
  if( (index >= 6) || !m_lower[index] )
    return;

  assert( m_allow_blank || lower.has_value() );
  const auto set_spin = [this]( NativeFloatSpinBox *spin, const std::optional<double> &value ){
    if( value.has_value() )
      spin->setValue( static_cast<float>( value.value() ) );
    else if( m_allow_blank )
      spin->setText( "" );
  };

  set_spin( m_lower[index], lower );
  if( m_upper[index] )
    set_spin( m_upper[index], (upper.has_value() || m_allow_blank) ? upper : lower );
}//void setValue(...)


std::optional<double> SkewParamsGrid::lowerValue( const size_t index ) const
{
  // Without blanks allowed, `value()` puts back the last value for an emptied field
  if( (index >= 6) || !m_lower[index] || (m_allow_blank && m_lower[index]->text().empty()) )
    return std::nullopt;
  return static_cast<double>( m_lower[index]->value() );
}


std::optional<double> SkewParamsGrid::upperValue( const size_t index ) const
{
  if( (index >= 6) || !m_upper[index] || (m_allow_blank && m_upper[index]->text().empty()) )
    return std::nullopt;
  return static_cast<double>( m_upper[index]->value() );
}


bool SkewParamsGrid::isFit( const size_t index ) const
{
  return (index < 6) && m_fit[index] && m_fit[index]->isChecked();
}


void SkewParamsGrid::setFit( const size_t index, const bool fit )
{
  if( (index < 6) && m_fit[index] )
    m_fit[index]->setChecked( fit );
}


void SkewParamsGrid::setResults( const std::vector<Wt::WString> &results, const std::vector<Wt::WString> &tooltips )
{
  const bool show = !results.empty() && m_result_header;
  toggleStyleClass( "SpgRes", show );
  if( m_result_header )
    m_result_header->setHidden( !show );

  for( size_t i = 0; i < 6; ++i )
  {
    if( !m_result[i] )
      continue;
    m_result[i]->setText( (show && (i < results.size())) ? results[i] : WString() );
    m_result[i]->setToolTip( (show && (i < tooltips.size())) ? tooltips[i] : WString() );
    m_result[i]->setHidden( !show );
  }
}//void setResults(...)


void SkewParamsGrid::setEditable( const bool editable )
{
  m_editable = editable;
  for( size_t i = 0; i < 6; ++i )
  {
    if( m_lower[i] )
      m_lower[i]->setEnabled( editable );
    if( m_upper[i] )
      m_upper[i]->setEnabled( editable );
    if( m_fit[i] )
      m_fit[i]->setEnabled( editable );
  }
}//void setEditable( const bool editable )


Wt::Signal<> &SkewParamsGrid::userChanged()
{
  return m_user_changed;
}
