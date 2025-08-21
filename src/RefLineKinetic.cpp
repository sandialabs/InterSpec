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

#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/RefLineKinetic.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/ExternalRidResult.h"
#include "InterSpec/ReferenceLinePredef.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

using namespace std;
using namespace Wt;


struct AlwaysSrcs
{
  const map<string,ReferenceLinePredef::NucMix> nuc_mixes;
  const map<string,ReferenceLinePredef::CustomSrcLines> custom_lines;
  const std::vector<ReferenceLinePredef::IndividualSource> individual_sources;
  
  AlwaysSrcs() = delete;
  AlwaysSrcs( map<string,ReferenceLinePredef::NucMix> &&mixes,
             map<string,ReferenceLinePredef::CustomSrcLines> &&lines,
             std::vector<ReferenceLinePredef::IndividualSource> &&srcs )
  : nuc_mixes( std::move(mixes) ),
  custom_lines( std::move(lines) ),
  individual_sources( std::move(srcs) )
  {
    
  }
};//struct AlwaysSrcs



RefLineKinetic::RefLineKinetic( D3SpectrumDisplayDiv *chart, InterSpec *parent )
  : Wt::WObject( parent ),
  m_interspec( parent ),
  m_chart( chart ),
  m_active( false ),
  m_has_inited( false ),
  m_init_error_msg{},
  m_always_srcs{}
{
  if( !m_interspec || m_chart )
    throw std::runtime_error( "RefLineKinetic: null InterSpec parent or chart" );
  
  m_active = UserPreferences::preferenceValue<bool>( "KineticRefLine", m_interspec );
  m_interspec->preferences()->addCallbackWhenChanged( "KineticRefLine", this, &RefLineKinetic::setActive );
    
  m_interspec->hintPeaksSet().connect( boost::bind(&RefLineKinetic::autoSearchPeaksSet, this, boost::placeholders::_1) );
  m_interspec->displayedSpectrumChanged().connect( boost::bind(&RefLineKinetic::spectrumChanged, this,
    boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4
  ) );
  
  m_interspec->externalRidResultsRecieved().connect( boost::bind( &RefLineKinetic::autoRidResultsRecieved, this, boost::placeholders::_1 ) );
  
  start_init_always_sources();
}//RefLineKinetic constructor


RefLineKinetic::~RefLineKinetic()
{
}


void RefLineKinetic::start_init_always_sources()
{
  try
  {
    // TODO: move this to being done in a background thread
    m_has_inited = true;
    
    string always_defs_file = SpecUtils::append_path(InterSpec::writableDataDirectory(), "kinetic_ref_lines.xml");
    if( SpecUtils::is_file(always_defs_file) )
      always_defs_file = SpecUtils::append_path(InterSpec::staticDataDirectory(), "kinetic_ref_lines.xml");
    
    map<string,ReferenceLinePredef::NucMix> nuc_mixes;
    map<string,ReferenceLinePredef::CustomSrcLines> custom_lines;
    std::vector<ReferenceLinePredef::IndividualSource> indiv_sources;
    ReferenceLinePredef::load_ref_line_file( always_defs_file, nuc_mixes, custom_lines, &indiv_sources );
    
    m_always_srcs = make_unique<AlwaysSrcs>( std::move(nuc_mixes), std::move(custom_lines), std::move(indiv_sources) );
  }catch( std::exception &e )
  {
    cerr << "Failed to initialize RefLineKinetic: " << e.what() << endl;
    m_init_error_msg = e.what();
  }//try / catch
}//void start_init_always_sources()


void RefLineKinetic::setActive( bool active )
{
  if( m_active == active )
    return;
  
  m_active = active;
  
  
}//void RefLineKinetic::setActive( bool active )


bool RefLineKinetic::isActive() const
{
  return m_active;
}


void RefLineKinetic::autoSearchPeaksSet( const SpecUtils::SpectrumType spectrum )
{
  updateLines();
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


void RefLineKinetic::updateLines()
{
  if( !m_active )
  {
    m_chart->setKineticRefernceLines( {}, "" );
    return;
  }//if( !m_active )
  
  if( m_always_srcs )
  {
    
  }//if( m_always_srcs )
}//void updateLines()
