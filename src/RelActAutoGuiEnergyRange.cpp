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
#include <string>
#include <optional>
#include <algorithm>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WString.h>
#include <Wt/WComboBox.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RelActAutoGuiEnergyRange.h"

using namespace std;
using namespace Wt;



namespace
{
  using RangeType = RelActCalcAuto::RoiRange::RangeLimitsType;

  // The order of the items in the range type combo box.  The deprecated `CanExpandForFwhm` is not
  //  offered; a ROI with it is shown as (and becomes) `CanBeBrokenUp`.
  int range_type_index( const RangeType type )
  {
    switch( type )
    {
      case RangeType::Fixed:            return 0;
      case RangeType::LineAnchored:     return 1;
      case RangeType::CanBeBrokenUp:    return 2;
      case RangeType::CanExpandForFwhm: return 2;
    }
    assert( 0 );
    return 0;
  }//range_type_index(...)

  RangeType range_type_from_index( const int index )
  {
    switch( index )
    {
      case 1:  return RangeType::LineAnchored;
      case 2:  return RangeType::CanBeBrokenUp;
      default: return RangeType::Fixed;
    }
  }//range_type_from_index(...)



  void set_field( NativeFloatSpinBox *field, const std::optional<double> &value, const bool is_coverage )
  {
    if( !value.has_value() )
      field->setText( "" );
    else
      field->setValue( static_cast<float>( is_coverage ? (100.0*(1.0 - value.value())) : value.value() ) );
  }//set_field(...)
}//namespace


RelActAutoGuiEnergyRange::RelActAutoGuiEnergyRange()
  : WContainerWidget(),
  m_lower_label( nullptr ),
  m_upper_label( nullptr ),
  m_lower_energy( nullptr ),
  m_upper_energy( nullptr ),
  m_continuum_type( nullptr ),
  m_range_type( nullptr ),
  m_to_individual_rois( nullptr ),
  m_edges( nullptr ),
  m_lower_coverage( nullptr ),
  m_lower_sideband( nullptr ),
  m_upper_coverage( nullptr ),
  m_upper_sideband( nullptr ),
  m_fit_range_txt( nullptr ),
  m_fit_range( 0.0, 0.0 ),
  m_auto_continuum_seed( PeakContinuum::OffsetType::Linear ),
  m_non_fixed_type( RelActCalcAuto::RoiRange::RangeLimitsType::CanBeBrokenUp ),
  m_shown_type( RelActCalcAuto::RoiRange::RangeLimitsType::LineAnchored ),
  m_highlight_region_id( 0 )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
    app->useMessageResourceBundle( "RelActAutoGuiEnergyRange" );
  addStyleClass( "RelActAutoGuiEnergyRange" );
  
  wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
  
  InterSpec * const interspec = InterSpec::instance();
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", interspec );
  const bool isPhone = interspec && interspec->isPhone();

  m_lower_label = addNew<WLabel>( WString::tr("raager-lower-energy") );
  m_lower_label->addStyleClass( "GridFirstCol GridFirstRow" );

  m_lower_energy = addNew<NativeFloatSpinBox>();
  m_lower_energy->addStyleClass( "GridSecondCol GridFirstRow" );
  m_lower_energy->setSpinnerHidden( true );
  m_lower_label->setBuddy( m_lower_energy );

  WLabel *label = addNew<WLabel>( WString::tr("keV") );
  label->addStyleClass( "GridThirdCol GridFirstRow" );

  m_upper_label = addNew<WLabel>( WString::tr("raager-upper-energy") );
  m_upper_label->addStyleClass( "GridFirstCol GridSecondRow" );

  m_upper_energy = addNew<NativeFloatSpinBox>();
  m_upper_energy->addStyleClass( "GridSecondCol GridSecondRow" );
  m_upper_energy->setSpinnerHidden( true );
  m_upper_label->setBuddy( m_upper_energy );

  label = addNew<WLabel>( WString::tr("keV") );
  label->addStyleClass( "GridThirdCol GridSecondRow" );

  m_lower_energy->valueChanged().connect( this, &RelActAutoGuiEnergyRange::handleEnergyChange );
  m_upper_energy->valueChanged().connect( this, &RelActAutoGuiEnergyRange::handleEnergyChange );

  label = addNew<WLabel>( WString::tr(isPhone ? "raager-continuum-type-short" : "raager-continuum-type") );
  label->addStyleClass( "GridFourthCol GridFirstRow" );
  m_continuum_type = addNew<WComboBox>();
  m_continuum_type->addStyleClass( "GridFifthCol GridFirstRow" );
  label->setBuddy( m_continuum_type );

  // Item 0 chooses the continuum automatically; we wont allow "External" here
  m_continuum_type->addItem( WString::tr("raager-continuum-auto") );
  for( int i = 0; i < static_cast<int>(PeakContinuum::OffsetType::External); ++i )
  {
    const char *key = PeakContinuum::offset_type_label_tr( PeakContinuum::OffsetType(i) );
    m_continuum_type->addItem( WString::tr(key) );
  }//for( loop over PeakContinuum::OffsetType )

  m_continuum_type->setCurrentIndex( 0 ); //Auto, starting from `m_auto_continuum_seed`
  m_continuum_type->changed().connect( this, &RelActAutoGuiEnergyRange::handleContinuumTypeChange );
  HelpSystem::attachToolTipOn( m_continuum_type, WString::tr("raager-continuum-type-tt"), showToolTips );

  label = addNew<WLabel>( WString::tr("raager-range-type") );
  label->addStyleClass( "GridFourthCol GridSecondRow" );
  m_range_type = addNew<WComboBox>();
  m_range_type->addStyleClass( "GridFifthCol GridSecondRow" );
  label->setBuddy( m_range_type );
  m_range_type->addItem( WString::tr("raager-range-fixed") );
  m_range_type->addItem( WString::tr("raager-range-line-anchored") );
  m_range_type->addItem( WString::tr("raager-range-broken-up") );
  m_range_type->setCurrentIndex( range_type_index( RelActCalcAuto::RoiRange::RangeLimitsType::Fixed ) );
  m_range_type->changed().connect( this, &RelActAutoGuiEnergyRange::handleRangeTypeChange );
  HelpSystem::attachToolTipOn( m_range_type, WString::tr("raager-range-type-tt"), showToolTips );

  WPushButton *removeEnergyRange = addNew<WPushButton>();
  removeEnergyRange->setStyleClass( "DeleteEnergyRangeOrNuc GridSixthCol GridFirstRow Wt-icon" );
  removeEnergyRange->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  removeEnergyRange->clicked().connect( this, &RelActAutoGuiEnergyRange::handleRemoveSelf );

  m_to_individual_rois = addNew<WPushButton>();
  m_to_individual_rois->setStyleClass( "ToIndividualRois GridSixthCol GridSecondRow Wt-icon" );
  m_to_individual_rois->setIcon( "InterSpec_resources/images/expand_list.svg" );
  m_to_individual_rois->clicked().connect( this, &RelActAutoGuiEnergyRange::emitSplitRangesRequested );
  HelpSystem::attachToolTipOn( m_to_individual_rois, WString::tr("raager-split-to-individual-rois-tt"), showToolTips );

  // How far a non-fixed ROI extends past its lines; empty fields use the defaults.
  m_edges = addNew<WContainerWidget>();
  m_edges->addStyleClass( "RoiEdges GridFirstCol GridThirdRow GridSpanFiveCol" );

  const auto add_edge_fields = [this]( const char *side_key, NativeFloatSpinBox *&coverage,
                                       NativeFloatSpinBox *&sideband ){
    WContainerWidget *side = m_edges->addNew<WContainerWidget>();
    side->addStyleClass( "RoiEdge" );
    side->addNew<WLabel>( WString::tr(side_key) );

    // Empty fields (showing the default as placeholder text), rather than the spin box's initial 0,
    //  which would be taken as a zero sideband.
    coverage = side->addNew<NativeFloatSpinBox>();
    coverage->setSpinnerHidden( true );
    coverage->setText( "" );
    coverage->addStyleClass( "RoiEdgeCoverage" );
    coverage->valueChanged().connect( this, &RelActAutoGuiEnergyRange::handleEdgeChange );
    side->addNew<WLabel>( WString::tr("raager-edge-coverage-units") );

    sideband = side->addNew<NativeFloatSpinBox>();
    sideband->setSpinnerHidden( true );
    sideband->setText( "" );
    sideband->addStyleClass( "RoiEdgeSideband" );
    sideband->valueChanged().connect( this, &RelActAutoGuiEnergyRange::handleEdgeChange );
    side->addNew<WLabel>( WString::tr("raager-edge-sideband-units") );
  };//add_edge_fields

  add_edge_fields( "raager-edge-below", m_lower_coverage, m_lower_sideband );
  add_edge_fields( "raager-edge-above", m_upper_coverage, m_upper_sideband );
  HelpSystem::attachToolTipOn( m_edges, WString::tr("raager-edges-tt"), showToolTips );

  m_fit_range_txt = m_edges->addNew<WText>();
  m_fit_range_txt->addStyleClass( "RoiFitRange FainterTxt" );
  m_fit_range_txt->setHidden( true );

  setEdgeDefaults( RelActCalcAuto::RoiEdge{}, RelActCalcAuto::RoiEdge{} );
  updateForRangeType();
}//RelActAutoGuiEnergyRange constructor
  
  
bool RelActAutoGuiEnergyRange::isEmpty() const
{
  // A line-anchored ROI may be a single line (lower == upper)
  if( rangeLimitsType() == RelActCalcAuto::RoiRange::RangeLimitsType::LineAnchored )
    return (m_upper_energy->value() <= 0.0);

  return ((fabs(m_lower_energy->value() - m_upper_energy->value() ) < 1.0)
            || (m_upper_energy->value() <= 0.0));
}
  
  
void RelActAutoGuiEnergyRange::handleRemoveSelf()
{
  m_remove_energy_range.emit();
}//void handleRemoveSelf()
  
  
void RelActAutoGuiEnergyRange::handleContinuumTypeChange()
{
  m_updated.emit();
}
  

void RelActAutoGuiEnergyRange::handleRangeTypeChange()
{
  rangeTypeChanged();
  m_updated.emit();
}


void RelActAutoGuiEnergyRange::rangeTypeChanged()
{
  using RangeType = RelActCalcAuto::RoiRange::RangeLimitsType;
  const RangeType type = rangeLimitsType();

  // A line-anchored ROI's energies are its bounding lines, not a range (a single line would give a
  //  zero-width range, and a range to split by lines is clipped to its energies, which would cut the
  //  peaks of the end lines in half), so changing to another type starts from the range last fit.
  if( (m_shown_type == RangeType::LineAnchored) && (type != RangeType::LineAnchored) )
  {
    if( m_fit_range.first < m_fit_range.second )
    {
      m_lower_energy->setValue( static_cast<float>( m_fit_range.first ) );
      m_upper_energy->setValue( static_cast<float>( m_fit_range.second ) );
      setFitRange( 0.0, 0.0 );
    }else if( !(m_lower_energy->value() < m_upper_energy->value()) )
    {
      // Not fit yet; rather than an empty range, use the few keV around the line the chart shows.
      const float energy = m_lower_energy->value();
      m_lower_energy->setValue( std::max( 0.0f, energy - 2.0f ) );
      m_upper_energy->setValue( energy + 2.0f );
    }
  }

  if( !forceFullRange() )
    m_non_fixed_type = type;

  updateForRangeType();
}//void rangeTypeChanged()


void RelActAutoGuiEnergyRange::handleEdgeChange()
{
  // A value the fit can not use (e.g., covering 0% of the peak) would stay displayed while the default
  //  is what gets used; clear it, so the default placeholder shows instead.
  for( NativeFloatSpinBox * const field : { m_lower_coverage, m_upper_coverage } )
  {
    if( !field->text().empty() && !edgeFieldValue( field, true ).has_value() )
      field->setText( "" );
  }

  for( NativeFloatSpinBox * const field : { m_lower_sideband, m_upper_sideband } )
  {
    if( !field->text().empty() && !edgeFieldValue( field, false ).has_value() )
      field->setText( "" );
  }

  setFitRange( 0.0, 0.0 ); //No longer what this ROI would be fit over
  m_updated.emit();
}


void RelActAutoGuiEnergyRange::updateForRangeType()
{
  const RelActCalcAuto::RoiRange::RangeLimitsType type = rangeLimitsType();
  const bool is_fixed = (type == RelActCalcAuto::RoiRange::RangeLimitsType::Fixed);
  const bool line_anchored = (type == RelActCalcAuto::RoiRange::RangeLimitsType::LineAnchored);

  m_lower_label->setText( WString::tr( line_anchored ? "raager-lower-line" : "raager-lower-energy" ) );
  m_upper_label->setText( WString::tr( line_anchored ? "raager-upper-line" : "raager-upper-energy" ) );

  m_edges->setHidden( is_fixed );
  m_to_individual_rois->setHidden( m_to_individual_rois->isDisabled() || is_fixed );
  m_shown_type = type;
}//void updateForRangeType()
  
  
void RelActAutoGuiEnergyRange::enableSplitToIndividualRanges( const bool enable )
{
  m_to_individual_rois->setHidden( !enable || forceFullRange() );
  m_to_individual_rois->setEnabled( enable );
}
  
  
void RelActAutoGuiEnergyRange::handleEnergyChange()
{
  float lower = m_lower_energy->value();
  float upper = m_upper_energy->value();
  if( lower > upper )
  {
    m_lower_energy->setValue( upper );
    m_upper_energy->setValue( lower );
    
    std::swap( lower, upper );
  }//if( lower > upper )
  
  setFitRange( 0.0, 0.0 ); //No longer what this ROI would be fit over
  m_updated.emit();
}//void handleEnergyChange()
  
  
void RelActAutoGuiEnergyRange::setEnergyRange( float lower, float upper )
{
  if( lower > upper )
    std::swap( lower, upper );
  
  m_lower_energy->setValue( lower );
  m_upper_energy->setValue( upper );
  setFitRange( 0.0, 0.0 );
  
  m_updated.emit();
}//void setEnergyRange( float lower, float upper )

  
bool RelActAutoGuiEnergyRange::forceFullRange() const
{
  return (rangeLimitsType() == RelActCalcAuto::RoiRange::RangeLimitsType::Fixed);
}

  
void RelActAutoGuiEnergyRange::setForceFullRange( const bool force_full )
{
  if( force_full == forceFullRange() )
    return;

  setRangeLimitsType( force_full ? RelActCalcAuto::RoiRange::RangeLimitsType::Fixed : m_non_fixed_type );
  m_updated.emit();
}


RelActCalcAuto::RoiRange::RangeLimitsType RelActAutoGuiEnergyRange::rangeLimitsType() const
{
  return range_type_from_index( m_range_type->currentIndex() );
}


void RelActAutoGuiEnergyRange::setRangeLimitsType( const RelActCalcAuto::RoiRange::RangeLimitsType type )
{
  m_range_type->setCurrentIndex( range_type_index( type ) );
  rangeTypeChanged();
}//void setRangeLimitsType(...)

  
void RelActAutoGuiEnergyRange::setContinuumType( const PeakContinuum::OffsetType type )
{
  const int type_index = static_cast<int>( type );
  if( (type_index < 0)
     || (type_index >= static_cast<int>(PeakContinuum::OffsetType::External)) )
  {
    assert( 0 );
    return;
  }
  
  if( (1 + type_index) == m_continuum_type->currentIndex() )
    return;
  
  m_continuum_type->setCurrentIndex( 1 + type_index );
  m_updated.emit();
}//void setContinuumType( PeakContinuum::OffsetType type )

  
void RelActAutoGuiEnergyRange::setHighlightRegionId( const size_t chart_id )
{
  m_highlight_region_id = chart_id;
}

  
size_t RelActAutoGuiEnergyRange::highlightRegionId() const
{
  return m_highlight_region_id;
}


float RelActAutoGuiEnergyRange::lowerEnergy() const
{
  return m_lower_energy->value();
}


float RelActAutoGuiEnergyRange::upperEnergy() const
{
  return m_upper_energy->value();
}


void RelActAutoGuiEnergyRange::setFromRoiRange( const RelActCalcAuto::RoiRange &roi )
{
  m_lower_energy->setValue( roi.lower_energy );
  m_upper_energy->setValue( roi.upper_energy );

  // An "External" continuum can not be used by this tool; show it as linear.
  const PeakContinuum::OffsetType cont_type = (roi.continuum_type == PeakContinuum::OffsetType::External)
                                              ? PeakContinuum::OffsetType::Linear : roi.continuum_type;
  m_auto_continuum_seed = cont_type;
  m_continuum_type->setCurrentIndex( roi.auto_continuum ? 0 : (1 + static_cast<int>(cont_type)) );

  m_range_type->setCurrentIndex( range_type_index( roi.range_limits_type ) );
  if( !forceFullRange() )
    m_non_fixed_type = rangeLimitsType();

  set_field( m_lower_coverage, roi.lower_edge.tail_fraction, true );
  set_field( m_lower_sideband, roi.lower_edge.sideband_fwhm, false );
  set_field( m_upper_coverage, roi.upper_edge.tail_fraction, true );
  set_field( m_upper_sideband, roi.upper_edge.sideband_fwhm, false );

  setFitRange( 0.0, 0.0 );
  updateForRangeType();
  enableSplitToIndividualRanges( !forceFullRange() );
}//setFromRoiRange(...)


RelActCalcAuto::RoiRange RelActAutoGuiEnergyRange::toRoiRange() const
{
  RelActCalcAuto::RoiRange roi;

  roi.lower_energy = m_lower_energy->value();
  roi.upper_energy = m_upper_energy->value();

  const int cont_index = m_continuum_type->currentIndex();
  roi.auto_continuum = (cont_index <= 0);
  roi.continuum_type = roi.auto_continuum ? m_auto_continuum_seed : PeakContinuum::OffsetType( cont_index - 1 );

  roi.range_limits_type = rangeLimitsType();

  if( roi.range_limits_type != RelActCalcAuto::RoiRange::RangeLimitsType::Fixed )
  {
    roi.lower_edge.tail_fraction = edgeFieldValue( m_lower_coverage, true );
    roi.lower_edge.sideband_fwhm = edgeFieldValue( m_lower_sideband, false );
    roi.upper_edge.tail_fraction = edgeFieldValue( m_upper_coverage, true );
    roi.upper_edge.sideband_fwhm = edgeFieldValue( m_upper_sideband, false );
  }

  return roi;
}//RelActCalcAuto::RoiRange toRoiRange() const


std::optional<double> RelActAutoGuiEnergyRange::edgeFieldValue( NativeFloatSpinBox *field, const bool is_coverage )
{
  // Parse the text rather than use the float value, which would turn e.g. 99.9% into 9.99985E-4.
  const std::string text = field->text().toUTF8();
  double value = 0.0;
  if( text.empty() || !SpecUtils::parse_double( text.c_str(), text.size(), value ) )
    return std::nullopt;

  if( is_coverage )
  {
    if( !(value > 50.0) || !(value < 100.0) )
      return std::nullopt;
    return std::round( 1.0E10*(100.0 - value) ) / 1.0E12; //rounded, so 99.9% gives exactly 1E-3
  }

  if( !(value >= 0.0) || !(value < 100.0) )
    return std::nullopt;
  return value;
}//edgeFieldValue(...)


void RelActAutoGuiEnergyRange::setFitRange( const double lower, const double upper )
{
  const bool valid = (lower < upper);
  m_fit_range = valid ? std::make_pair( lower, upper ) : std::make_pair( 0.0, 0.0 );

  if( valid && !forceFullRange() )
    m_fit_range_txt->setText( WString::tr("raager-fit-range")
                              .arg( SpecUtils::printCompact(lower, 5) )
                              .arg( SpecUtils::printCompact(upper, 5) ) );
  m_fit_range_txt->setHidden( !valid || forceFullRange() );
}//void setFitRange(...)


std::pair<double,double> RelActAutoGuiEnergyRange::displayRange() const
{
  if( m_fit_range.first < m_fit_range.second )
    return m_fit_range;
  return { m_lower_energy->value(), m_upper_energy->value() };
}//displayRange()


void RelActAutoGuiEnergyRange::setEdgeDefaults( const RelActCalcAuto::RoiEdge &lower,
                                                const RelActCalcAuto::RoiEdge &upper )
{
  const auto coverage_str = []( const double tail_fraction ) -> WString {
    return WString::fromUTF8( SpecUtils::printCompact( 100.0*(1.0 - tail_fraction), 6 ) );
  };
  const auto sideband_str = []( const double sideband ) -> WString {
    return WString::fromUTF8( SpecUtils::printCompact( sideband, 4 ) );
  };

  m_lower_coverage->setPlaceholderText( coverage_str( lower.tail_fraction ) );
  m_lower_sideband->setPlaceholderText( sideband_str( lower.sideband_fwhm ) );
  m_upper_coverage->setPlaceholderText( coverage_str( upper.tail_fraction ) );
  m_upper_sideband->setPlaceholderText( sideband_str( upper.sideband_fwhm ) );
}//void setEdgeDefaults(...)


Wt::Signal<> &RelActAutoGuiEnergyRange::updated()
{
  return m_updated;
}


Wt::Signal<> &RelActAutoGuiEnergyRange::remove()
{
  return m_remove_energy_range;
}


Wt::Signal<RelActAutoGuiEnergyRange *> &RelActAutoGuiEnergyRange::splitRangesRequested()
{
  return m_split_ranges_requested;
}


void RelActAutoGuiEnergyRange::emitSplitRangesRequested()
{
  m_split_ranges_requested.emit(this);
}
