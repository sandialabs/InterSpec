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
#include <limits>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WContainerWidget>


#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RelActAutoGuiFreePeak.h"


using namespace std;
using namespace Wt;


RelActAutoGuiFreePeak::RelActAutoGuiFreePeak( WContainerWidget *parent )
  : WContainerWidget( parent ),
  m_energy( nullptr ),
  m_fwhm_constrained( nullptr ),
  m_apply_energy_cal( nullptr ),
  m_updated( this ),
  m_remove( this )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
    app->useMessageResourceBundle( "RelActAutoGuiFreePeak" );
  addStyleClass( "RelActAutoGuiFreePeak" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  
  WLabel *label = new WLabel( WString::tr("Energy"), this );
  label->addStyleClass( "GridFirstCol GridFirstRow" );
  
  m_energy = new NativeFloatSpinBox( this );
  label->setBuddy( m_energy );
  m_energy->setSpinnerHidden( true );
  m_energy->valueChanged().connect( this, &RelActAutoGuiFreePeak::handleEnergyChange );
  m_energy->addStyleClass( "GridSecondCol GridFirstRow" );
  
  // Things are a little cramped if we include the keV
  //label = new WLabel( "keV", this );
  //label->addStyleClass( "GridThirdCol GridFirstRow" );
  
  WPushButton *removeFreePeak = new WPushButton( this );
  removeFreePeak->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
  removeFreePeak->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  removeFreePeak->clicked().connect( this, &RelActAutoGuiFreePeak::handleRemoveSelf );
  removeFreePeak->addStyleClass( "GridThirdCol GridFirstRow" );
  
  m_fwhm_constrained = new WCheckBox( WString::tr("raagfp-constrain-fwhm"), this );
  m_fwhm_constrained->setChecked( true );
  m_fwhm_constrained->checked().connect( this, &RelActAutoGuiFreePeak::handleFwhmConstrainChanged );
  m_fwhm_constrained->unChecked().connect( this, &RelActAutoGuiFreePeak::handleFwhmConstrainChanged );
  m_fwhm_constrained->addStyleClass( "FreePeakConstrain GridFirstCol GridSecondRow GridSpanThreeCol" );
  
  HelpSystem::attachToolTipOn( m_fwhm_constrained, WString::tr("raagfp-constrain-fwhm-tt"), showToolTips );
  
  
  m_apply_energy_cal = new WCheckBox( WString::tr("raagfp-true-energy"), this );
  m_apply_energy_cal->setStyleClass( "CbNoLineBreak" );
  m_apply_energy_cal->setChecked( true );
  m_apply_energy_cal->checked().connect( this, &RelActAutoGuiFreePeak::handleApplyEnergyCalChanged );
  m_apply_energy_cal->unChecked().connect( this, &RelActAutoGuiFreePeak::handleApplyEnergyCalChanged );
  m_apply_energy_cal->addStyleClass( "FreePeakConstrain GridFirstCol GridThirdRow GridSpanThreeCol" );
  HelpSystem::attachToolTipOn( m_apply_energy_cal, WString::tr("raagfp-true-energy-tt"), showToolTips );
  
  
  m_invalid = new WText( WString::tr("raagfp-not-in-roi"), this );
  m_invalid->addStyleClass( "InvalidFreePeakEnergy GridFirstCol GridFourthRow GridSpanThreeCol" );
  m_invalid->hide();
}//RelActAutoGuiFreePeak constructor
  

float RelActAutoGuiFreePeak::energy() const
{
  return m_energy->value();
}


void RelActAutoGuiFreePeak::setEnergy( const float energy )
{
  m_energy->setValue( energy );
  m_updated.emit();
}


bool RelActAutoGuiFreePeak::fwhmConstrained() const
{
  return m_fwhm_constrained->isChecked();
}


void RelActAutoGuiFreePeak::setFwhmConstrained( const bool constrained )
{
  m_fwhm_constrained->setChecked( constrained );
  m_updated.emit();
}


bool RelActAutoGuiFreePeak::applyEnergyCal() const
{
  return (m_apply_energy_cal->isVisible() && m_apply_energy_cal->isChecked());
}


void RelActAutoGuiFreePeak::setApplyEnergyCal( const bool apply )
{
  m_apply_energy_cal->setChecked( apply );
}


void RelActAutoGuiFreePeak::setInvalidEnergy( const bool invalid )
{
  if( invalid != m_invalid->isVisible() )
    m_invalid->setHidden( !invalid );
}


void RelActAutoGuiFreePeak::handleRemoveSelf()
{
  m_remove.emit();
}


void RelActAutoGuiFreePeak::handleEnergyChange()
{
  m_updated.emit();
}//


void RelActAutoGuiFreePeak::handleFwhmConstrainChanged()
{
  m_updated.emit();
}


void RelActAutoGuiFreePeak::handleApplyEnergyCalChanged()
{
  m_updated.emit();
}


void RelActAutoGuiFreePeak::setApplyEnergyCalVisible( const bool visible )
{
  if( m_apply_energy_cal->isVisible() != visible )
    m_apply_energy_cal->setHidden( !visible );
}


Wt::Signal<> &RelActAutoGuiFreePeak::updated()
{
  return m_updated;
}


Wt::Signal<> &RelActAutoGuiFreePeak::remove()
{
  return m_remove;
}

