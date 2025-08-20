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

#include "InterSpec/InterSpec.h"
#include "InterSpec/RefLineKinetic.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/ExternalRidResult.h"

using namespace std;
using namespace Wt;

RefLineKinetic::RefLineKinetic( InterSpec *parent )
  : Wt::WObject( parent ),
  m_interspec( parent ),
  m_active( false )
{
  if( !m_interspec )
    throw std::runtime_error( "RefLineKinetic: null InterSpec parent" );
  
  m_active = UserPreferences::preferenceValue<bool>( "KineticRefLine", m_interspec );
  m_interspec->preferences()->addCallbackWhenChanged( "KineticRefLine", this, &RefLineKinetic::setActive );
    
  m_interspec->hintPeaksSet().connect( boost::bind(&RefLineKinetic::autoSearchPeaksSet, this, boost::placeholders::_1) );
  m_interspec->displayedSpectrumChanged().connect( boost::bind(&RefLineKinetic::spectrumChanged, this,
    boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4
  ) );
  
  m_interspec->externalRidResultsRecieved().connect( boost::bind( &RefLineKinetic::autoRidResultsRecieved, this, boost::placeholders::_1 ) );
  
  // Setup refreshing whenever `reflines->setExternalRidResults(...)` gets called - perhaps setup a more signal/slot mechanism for when results are avaialble
}//RefLineKinetic constructor


RefLineKinetic::~RefLineKinetic()
{
}


void RefLineKinetic::setActive( bool active )
{
  if( m_active == active )
    return;
  m_active = active;
}


bool RefLineKinetic::isActive() const
{
  return m_active;
}

void RefLineKinetic::handlePreferenceChange( bool active )
{
  m_active = active;
}


void RefLineKinetic::autoSearchPeaksSet( const SpecUtils::SpectrumType spectrum )
{
  
}//void autoSearchPeaksSet(...)


void RefLineKinetic::spectrumChanged( const SpecUtils::SpectrumType spec_type,
                     const std::shared_ptr<SpecMeas> &measurement,
                     const std::set<int> &sample_numbers,
                     const std::vector<std::string> &detectors )
{
  
}//void spectrumChanged(...)


void RefLineKinetic::autoRidResultsRecieved( const std::shared_ptr<const ExternalRidResults> &results )
{
  
}//void autoRidResultsRecieved( const std::shared_ptr<const ExternalRidResults> &results )
