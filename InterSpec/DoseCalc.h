#ifndef DoseCalc_h
#define DoseCalc_h
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
#include <vector>


class MaterialDB;
class DoseCalcWidget;
class InterSpec;
class ShieldingSelect;
class GadrasScatterTable;


namespace Wt
{
  class WText;
  class WComboBox;
  class WLineEdit;
  class WPushButton;
  class WGridLayout;
  class WButtonGroup;
  class WStackedWidget;
  class WSuggestionPopup;
}//namespace Wt

/* Todo as of 20151001
     -Add g, kg for self attenuating sources
     -Continuum groups for self attenuating sources
     -Continuum normalization
     -Neutron groups for a number of sources
     -Make shield amount solver
 */


namespace DoseCalc
{
  /** Computes the fluence-to-dose for gammas at a given energy, using 
   * ANSI/ANS-6.1.1-1991.
   *
   * \param energy Gamma energy in keV
   *
   * \returns Returns dose in units of PhysicalUnits.  That is, to convert
   *          into a human readable dose at a given distance, you would divide
   *          by PhysicalUnits::sievert/PhysicalUnits::hour and divide by
   *          4*pi*r*r (where r is in dimensions of PhysicalUnits to), to give
   *          sv/hr.
   */
  float gamma_dose( const float energy );
  
  
  /** Computes the fluence-to-dose for neutrons at a given energy, using 
   * ANSI/ANS-6.1.1-1991.
   *
   * \param energy Gamma energy in keV
   *
   * \returns Returns dose in units of PhysicalUnits.  That is, to convert
   *          into a human readable dose at a given distance, you would divide
   *          by PhysicalUnits::sievert/PhysicalUnits::hour and divide by
   *          4*pi*r*r (where r is in dimensions of PhysicalUnits to), to give
   *          sv/hr.
   */
  float neutron_dose( const float energy );
  
  
  /** \brief Compute the human dose in Seiverts based on radiation leakage and
   the distance from the source.
   
   The fluence-to-dose functions were taken from ANSI/ANS-6.1.1-1991 polynomial
   fits for the most conservative case: Anterior-Posterior orientation, Table 5 
   and Table 6. In this program the average fluence-to-dose factor for an energy 
   interval is taken to be the fluence-to-dose value of the average energy of 
   the interval.
   
   Adapted from GADRAS code iGetDose.f, routine IGETDOSE.
   
   \param energies Energies of incomming gammas in keV.
   \param intensity Intensities of cooresponding gammas. A weight of 1 equals
          1 source gama, before shielding.
   \param areal_density Areal density, in units of PhysicalUnits (e.g. would
          have to divide by g/cm2 to print out g/cm2) of shielding
   \param atomic_number Atomic numberof the shielding
   \param distance Distance from center of source to location of interest
   \param xstool The XS tool to use for attenuation calculations
   \param scatter The GadrasScatterTable object to use to create the continuum
   
   \returns Returns dose in units of PhysicalUnits.  That is, to convert
            into a human readable dose at a given distance, you would divide
            by PhysicalUnits::sievert/PhysicalUnits::hour to give sv/hr.
   */
  double gamma_dose_with_shielding( const std::vector<float> &energies,
                                    const std::vector<float> &intensity,
                                    const float areal_density,
                                    const float atomic_number,
                                    const float distance,
                                    const GadrasScatterTable &scatter );
}//namespace DoseCalc
#endif //DoseCalc_h
