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

#include "InterSpec/AuxWindow.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/EnergyCalMultiFile.h"
#include "InterSpec/EnergyCalAddActions.h"


using namespace std;
using namespace Wt;

EnergyCalAddActionsWindow::EnergyCalAddActionsWindow( const MoreActionsIndex actionType,
                             const std::vector<MeasToApplyCoefChangeTo> &measToChange,
                             EnergyCalTool *calibrator )
  : AuxWindow( "&nbsp", (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                      | AuxWindowProperties::DisableCollapse) ),
    m_actionType( actionType ),
    m_measToChange( make_shared<vector<MeasToApplyCoefChangeTo>>(measToChange)),
    m_calibrator( calibrator )
{
  switch( m_actionType )
  {
    case MoreActionsIndex::Linearize:
      break;
      
    case MoreActionsIndex::Truncate:
      break;
      
    case MoreActionsIndex::CombineChannels:
      break;
      
    case MoreActionsIndex::ConvertToFrf:
      break;
      
    case MoreActionsIndex::ConvertToPoly:
      break;
      
    case MoreActionsIndex::MultipleFilesCal:
      AuxWindow::setWindowTitle( "Multi-File Calibration" );
      new EnergyCalMultiFile( m_calibrator, this );
      break;
      
    case MoreActionsIndex::NumMoreActionsIndex:
      break;
  }//switch( m_actionType )
  
  rejectWhenEscapePressed();
  AuxWindow::show();
  setWidth( 400.0 );
  AuxWindow::resizeToFitOnScreen();
  AuxWindow::centerWindow();
}//EnergyCalAddActionsWindow constructor
  
EnergyCalAddActionsWindow::~EnergyCalAddActionsWindow()
{
}
 
