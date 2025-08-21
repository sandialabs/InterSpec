#ifndef RefLineKinetic_h
#define RefLineKinetic_h
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

#include <set>
#include <memory>
#include <string>
#include <vector>


#include <Wt/WObject>

namespace SpecUtils{ enum class SpectrumType : int; }

class SpecMeas;
class InterSpec;
struct AlwaysSrcs;
struct ExternalRidResults;
class D3SpectrumDisplayDiv;

class RefLineKinetic : public Wt::WObject
{
public:
  RefLineKinetic( D3SpectrumDisplayDiv *chart, InterSpec *parent );
  virtual ~RefLineKinetic();
  
  bool successfully_initialized() const;
  
  void setActive( bool active );
  bool isActive() const;
  
protected:
  void start_init_always_sources();
  
  void autoSearchPeaksSet( const SpecUtils::SpectrumType spectrum );
  void spectrumChanged( const SpecUtils::SpectrumType spec_type,
                       const std::shared_ptr<SpecMeas> &measurement,
                       const std::set<int> &sample_numbers,
                       const std::vector<std::string> &detectors );
  void autoRidResultsRecieved( const std::shared_ptr<const ExternalRidResults> &results );

  void updateLines();
  
  InterSpec *m_interspec;
  D3SpectrumDisplayDiv *m_chart;
  
  bool m_active;
  
  /** Will be set to false until the "always" sources have been attempted to initialize. */
  bool m_has_inited;
  
  /** Will be empty if the "always" sources were succesfully initialized, and non-empty otherwise. */
  std::string m_init_error_msg;
  
  std::unique_ptr<const AlwaysSrcs> m_always_srcs;
};//class RefLineKinetic

#endif // RefLineKinetic_h
