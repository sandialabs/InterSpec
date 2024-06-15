#ifndef NuclideSourceEnter_h
#define NuclideSourceEnter_h
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

#include <map>
#include <string>
#include <vector>

#include <Wt/WString>
#include <Wt/WContainerWidget>


namespace Wt
{
  class WText;
  class WLineEdit;
}
  

namespace SandiaDecay 
{
  struct Nuclide;
}

/** Class to implement logic of reacting to user input change of nuclide or age.
 */
class NuclideSourceEnterController : public Wt::WObject
{
public:
  /** Takes pointers to the user display widgets - does not take ownership of them.
   
   Sets appropriate validators
   */
  NuclideSourceEnterController( Wt::WLineEdit *nuclideEdit,
                               Wt::WLineEdit *nuclideAgeEdit,
                               Wt::WText *halfLifeTxt,
                               Wt::WObject *parent );
  
  
  void handleNuclideUserInput();
  
  void handleAgeUserChange();
  
  double nuclideAge() const;
  Wt::WString nuclideAgeStr() const;
  const SandiaDecay::Nuclide *nuclide() const;
  
  void setNuclideText( const std::string &txt );
  
  void setNuclideAgeTxt( const std::string &txt );
  
  Wt::Signal<> &changed();
protected:
  Wt::WLineEdit *m_nuclideEdit;
  Wt::WLineEdit *m_nuclideAgeEdit;
  Wt::WText *m_halfLifeTxt;
  
  const SandiaDecay::Nuclide *m_currentNuc;
  std::map<const SandiaDecay::Nuclide *,std::string> m_prevAgeTxt;
  
  Wt::Signal<> m_changed;
};//class NuclideSourceEnterController



/** A class to handle user-input of nuclides, and their ages.
 * 
 * Currently (20240601) only used for DoseCalc tool, and simple activity limit.
 * Only allows user to enter nuclides, and is pretty basic, but has been separated 
 * out because it could be useful in more places - but probably needs some work
 * to make it a little more general and consistent with the Reference Photopeak tool
*/
class NuclideSourceEnter : public Wt::WContainerWidget
{
public:
  NuclideSourceEnter( const bool showHalfLife, const bool showToolTips, Wt::WContainerWidget *parent );
  virtual ~NuclideSourceEnter();
  
  void setNuclideText( const std::string &txt );
  
  void setNuclideAgeTxt( const std::string &txt );
  
  
  
  enum EquilibriumType
  {
    SecularEquilibrium, PromptEquilibrium, NoEquilibrium
  };
  
  
  static double getMaximumHalfLife( const SandiaDecay::Nuclide* nuclide );
  
  static void getMaximumHalfLifeInDescension( const SandiaDecay::Nuclide * const nuclide, double parentHalfLife, double &halfLife );
  
  static bool trySecularEquilibrium( const SandiaDecay::Nuclide * const nuclide, double &age );
  
  static bool tryPromptEquilibrium( const SandiaDecay::Nuclide * const nuclide, double &age );
  
  
  static EquilibriumType getEquilibrium( const SandiaDecay::Nuclide* nuclide, double &age, EquilibriumType eqCode );
  
  Wt::WString nuclideAgeStr() const;
  
  double nuclideAge() const;
  
  /** 
   * \returns empty results if no valid isotope, an invalid age, or negative or 
   *          zero activity. Other wise returns <energy,gamma/sec> pairs.
   */
  std::vector< std::pair<float,float> > photonEnergyAndIntensity( const double activity ) const;
  
  Wt::Signal<> &changed();
  
  const SandiaDecay::Nuclide *nuclide() const;
  
  
protected:
  Wt::WLineEdit *m_nuclideEdit;
  Wt::WLineEdit *m_nuclideAgeEdit;
  Wt::WText *m_halfLifeTxt;
  
  NuclideSourceEnterController *m_controller;
};//class NuclideSourceEnter


#endif //NuclideSourceEnter_h
