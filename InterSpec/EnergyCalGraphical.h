#ifndef EnergyCalGraphical_h
#define EnergyCalGraphical_h
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

#include <time.h>


#include "InterSpec/AuxWindow.h"

//Forward includes
class EnergyCalTool;
class NativeFloatSpinBox;
namespace Wt
{
  class WCheckBox;
  class WButtonGroup;
}

class EnergyCalGraphicalConfirm : public AuxWindow
{
public:
  enum RecalTypes
  {
    kOffset, kLinear,
    /*kQuadratic,*/
    kDeviation, NumRecalTypes
  };//enum RecalTypes
  
  EnergyCalGraphicalConfirm( double lowe, double highe, EnergyCalTool *cal,
                            time_t lastCal, RecalTypes lastType, float lastEnergy );
  virtual ~EnergyCalGraphicalConfirm();
  void apply();
  void setEnergies( double xstart, double xfinish );
  
protected:
  EnergyCalTool *m_calibrator;
  Wt::WButtonGroup *m_typeButtons;
  Wt::WCheckBox *m_foregroundOnly;
  NativeFloatSpinBox *m_startE, *m_finalE;
  Wt::WCheckBox *m_preserveLastCal;
  const RecalTypes m_lastType;
  const float m_lastEnergy;
};//class EnergyCalGraphicalConfirm


#endif // EnergyCalGraphical_h
