#ifndef EnergyCalActions_h
#define EnergyCalActions_h
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

#include  <memory>
#include  <vector>

#include "InterSpec/AuxWindow.h"

//Forward Declarations
class EnergyCalTool;
struct MeasToApplyCoefChangeTo;
enum class MoreActionsIndex : int;



/**
 */
class EnergyCalAddActionsWindow : public AuxWindow
{
public:
  EnergyCalAddActionsWindow( const MoreActionsIndex actionType,
                             const std::vector<MeasToApplyCoefChangeTo> &measToChange,
                             EnergyCalTool *calibrator );
  
  virtual ~EnergyCalAddActionsWindow();
  
protected:
  const MoreActionsIndex m_actionType;
  std::shared_ptr<std::vector<MeasToApplyCoefChangeTo>> m_measToChange;
  EnergyCalTool *m_calibrator;
};//class EnergyCalPreserveWindow




#endif //EnergyCalActions_h
