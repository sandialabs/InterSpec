#ifndef RelActAutoGuiNuclide_h
#define RelActAutoGuiNuclide_h
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
#include <variant>

#include <Wt/WColor>
#include <Wt/WContainerWidget>

//Forward declerations
namespace Wt
{
  class WLineEdit;
  class WLabel;
  class WCheckBox;
}

class ColorSelect;


class RelActAutoGuiNuclide : public Wt::WContainerWidget
{
public:
  RelActAutoGuiNuclide( Wt::WContainerWidget *parent = nullptr );
  
  
  void handleIsotopeChange();
  
  
  void setFitAgeVisible( bool visible, bool do_fit );
  
  
  void setAgeDisabled( bool disabled );
  
  
  void setAge( const std::string &age );
  
  void setNuclideEditFocus();
  
  void handleAgeChange();
  
  
  void handleFitAgeChange();
  
  
  void handleColorChange();
  
  /** Will return std::monostate if invalid text entered. */
  std::variant<std::monostate, const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> nuclide() const;
  
  Wt::WColor color() const;
  
  void setColor( const Wt::WColor &color );
  
  double age() const;
  
  RelActCalcAuto::NucInputInfo toNucInputInfo() const;
  
  void fromNucInputInfo( const RelActCalcAuto::NucInputInfo &info );
  
  void handleRemoveSelf();
  
  Wt::Signal<> &updated();
  
  Wt::Signal<> &remove();
  
  void addActRatioConstraint( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &constraint );
  void addMassFractionConstraint( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint );
  void setIsInCurves( const std::set<size_t> &curves_with_nuc, size_t num_rel_eff_curves );
  
protected:
  Wt::WLineEdit *m_nuclide_edit;
  Wt::WLabel *m_age_label;
  Wt::WLineEdit *m_age_edit;
  Wt::WCheckBox *m_fit_age;
  ColorSelect *m_color_select;
  
  Wt::Signal<> m_updated;
  Wt::Signal<> m_remove;
};//class RelActAutoGuiNuclide


#endif //RelActAutoGuiNuclide_h
