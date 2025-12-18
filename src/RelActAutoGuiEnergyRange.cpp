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

#include <string>

#include <Wt/WLabel>
#include <Wt/WString>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RelActAutoGuiEnergyRange.h"

using namespace std;
using namespace Wt;



RelActAutoGuiEnergyRange::RelActAutoGuiEnergyRange( WContainerWidget *parent )
  : WContainerWidget( parent ),
  m_updated( this ),
  m_remove_energy_range( this ),
  m_split_ranges_requested( this ),
  m_lower_energy( nullptr ),
  m_upper_energy( nullptr ),
  m_continuum_type( nullptr ),
  m_force_full_range( nullptr ),
  m_to_individual_rois( nullptr ),
  m_highlight_region_id( 0 )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
    app->useMessageResourceBundle( "RelActAutoGuiEnergyRange" );
  addStyleClass( "RelActAutoGuiEnergyRange" );
  
  wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  
  WLabel *label = new WLabel( WString::tr("raager-lower-energy"), this );
  label->addStyleClass( "GridFirstCol GridFirstRow" );
  
  m_lower_energy = new NativeFloatSpinBox( this );
  m_lower_energy->addStyleClass( "GridSecondCol GridFirstRow" );
  m_lower_energy->setSpinnerHidden( true );
  label->setBuddy( m_lower_energy );
  
  label = new WLabel( WString::tr("keV"), this );
  label->addStyleClass( "GridThirdCol GridFirstRow" );
  
  label = new WLabel( WString::tr("raager-upper-energy"), this );
  label->addStyleClass( "GridFirstCol GridSecondRow" );
  
  m_upper_energy = new NativeFloatSpinBox( this );
  m_upper_energy->addStyleClass( "GridSecondCol GridSecondRow" );
  m_upper_energy->setSpinnerHidden( true );
  label->setBuddy( m_upper_energy );
  
  label = new WLabel( WString::tr("keV"), this );
  label->addStyleClass( "GridThirdCol GridSecondRow" );
  
  m_lower_energy->valueChanged().connect( this, &RelActAutoGuiEnergyRange::handleEnergyChange );
  m_upper_energy->valueChanged().connect( this, &RelActAutoGuiEnergyRange::handleEnergyChange );
  
  label = new WLabel( WString::tr("raager-continuum-type"), this );
  label->addStyleClass( "GridFourthCol GridFirstRow" );
  m_continuum_type = new WComboBox( this );
  m_continuum_type->addStyleClass( "GridFifthCol GridFirstRow" );
  label->setBuddy( m_continuum_type );
  
  // We wont allow "External" here
  for( int i = 0; i < static_cast<int>(PeakContinuum::OffsetType::External); ++i )
  {
    const char *key = PeakContinuum::offset_type_label_tr( PeakContinuum::OffsetType(i) );
    m_continuum_type->addItem( WString::tr(key) );
  }//for( loop over PeakContinuum::OffsetType )
  
  m_continuum_type->setCurrentIndex( static_cast<int>(PeakContinuum::OffsetType::Linear) );
  m_continuum_type->changed().connect( this, &RelActAutoGuiEnergyRange::handleContinuumTypeChange );
  
  m_force_full_range = new WCheckBox( WString::tr("raager-force-full-range"), this );
  m_force_full_range->addStyleClass( "GridFourthCol GridSecondRow GridSpanTwoCol" );
  m_force_full_range->checked().connect( this, &RelActAutoGuiEnergyRange::handleForceFullRangeChange );
  m_force_full_range->unChecked().connect( this, &RelActAutoGuiEnergyRange::handleForceFullRangeChange );
  
  WPushButton *removeEnergyRange = new WPushButton( this );
  removeEnergyRange->setStyleClass( "DeleteEnergyRangeOrNuc GridSixthCol GridFirstRow Wt-icon" );
  removeEnergyRange->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  removeEnergyRange->clicked().connect( this, &RelActAutoGuiEnergyRange::handleRemoveSelf );
  
  m_to_individual_rois = new WPushButton( this );
  m_to_individual_rois->setStyleClass( "ToIndividualRois GridSixthCol GridSecondRow Wt-icon" );
  m_to_individual_rois->setIcon( "InterSpec_resources/images/expand_list.svg" );
  m_to_individual_rois->clicked().connect( this, &RelActAutoGuiEnergyRange::emitSplitRangesRequested );
  
  
  
  m_to_individual_rois->setHidden( m_force_full_range->isChecked() );
  
  HelpSystem::attachToolTipOn( m_to_individual_rois, WString::tr("raager-split-to-individual-rois-tt"), showToolTips );
}//RelActAutoGuiEnergyRange constructor
  
  
bool RelActAutoGuiEnergyRange::isEmpty() const
{
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
  

void RelActAutoGuiEnergyRange::handleForceFullRangeChange()
{
  m_to_individual_rois->setHidden( m_to_individual_rois->isDisabled() || m_force_full_range->isChecked() );
    
  m_updated.emit();
}
  
  
void RelActAutoGuiEnergyRange::enableSplitToIndividualRanges( const bool enable )
{
  m_to_individual_rois->setHidden( !enable || m_force_full_range->isChecked() );
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
  
  m_updated.emit();
}//void handleEnergyChange()
  
  
void RelActAutoGuiEnergyRange::setEnergyRange( float lower, float upper )
{
  if( lower > upper )
    std::swap( lower, upper );
  
  m_lower_energy->setValue( lower );
  m_upper_energy->setValue( upper );
  
  m_updated.emit();
}//void setEnergyRange( float lower, float upper )

  
bool RelActAutoGuiEnergyRange::forceFullRange() const
{
  return m_force_full_range->isChecked();
}

  
void RelActAutoGuiEnergyRange::setForceFullRange( const bool force_full )
{
  if( force_full == m_force_full_range->isChecked() )
    return;
  
  m_force_full_range->setChecked( force_full );
  m_updated.emit();
}

  
void RelActAutoGuiEnergyRange::setContinuumType( const PeakContinuum::OffsetType type )
{
  const int type_index = static_cast<int>( type );
  if( (type_index < 0)
     || (type_index >= static_cast<int>(PeakContinuum::OffsetType::External)) )
  {
    assert( 0 );
    return;
  }
  
  if( type_index == m_continuum_type->currentIndex() )
    return;
  
  m_continuum_type->setCurrentIndex( type_index );
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
  if( roi.continuum_type == PeakContinuum::OffsetType::External )
    throw runtime_error( "setFromRoiRange: External continuum type not supported" );

  if( roi.range_limits_type == RelActCalcAuto::RoiRange::RangeLimitsType::CanExpandForFwhm )
    throw runtime_error( "setFromRoiRange: CanExpandForFwhm range type not supported in GUI" );

  m_lower_energy->setValue( roi.lower_energy );
  m_upper_energy->setValue( roi.upper_energy );
  m_continuum_type->setCurrentIndex( static_cast<int>(roi.continuum_type) );

  const bool is_fixed = (roi.range_limits_type == RelActCalcAuto::RoiRange::RangeLimitsType::Fixed);
  m_force_full_range->setChecked( is_fixed );

  enableSplitToIndividualRanges( !is_fixed );
}//setFromRoiRange(...)


RelActCalcAuto::RoiRange RelActAutoGuiEnergyRange::toRoiRange() const
{
  RelActCalcAuto::RoiRange roi;

  roi.lower_energy = m_lower_energy->value();
  roi.upper_energy = m_upper_energy->value();
  roi.continuum_type = PeakContinuum::OffsetType( m_continuum_type->currentIndex() );

  const bool is_checked = m_force_full_range->isChecked();
  roi.range_limits_type = is_checked ? RelActCalcAuto::RoiRange::RangeLimitsType::Fixed
                                     : RelActCalcAuto::RoiRange::RangeLimitsType::CanBeBrokenUp;

  return roi;
}//RelActCalcAuto::RoiRange toRoiRange() const


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
