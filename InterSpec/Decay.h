#ifndef Decay_h
#define Decay_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

//
//  Decay.h
//  InterSpec
//
//  Created by Chan, Ethan on 12/18/13.
//
//

#include "InterSpec_config.h"

#include "InterSpec/AuxWindow.h"


class InterSpec;

namespace Wt
{
  class WText;
  class WLineEdit;
  class WDoubleSpinBox;
}//namespace Wt

class DecayActivityDiv;

class Decay : public AuxWindow
{
public:
  Decay( InterSpec *viewer );
  virtual ~Decay();
  
  void clearAllNuclides();
  void addNuclide( const int z, const int a, const int iso,
                  const double activity, const bool useCurries,
                  const double age, const double maxtime = -1.0 );
 protected:
  DecayActivityDiv *m_activityDiv;
};//class Decay


#endif //Decay
