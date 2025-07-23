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

#include <Wt/WColor>
#include <Wt/WContainerWidget>

#include "InterSpec/RelActCalcAuto.h"

//Forward declerations
namespace Wt
{
  class WText;
  class WLabel;
  class WLineEdit;
  class WCheckBox;
}

class ColorSelect;
class RelActAutoGui;
class RelActAutoGuiNuclideConstraint;


class RelActAutoGuiNuclide : public Wt::WContainerWidget
{
public:
  RelActAutoGuiNuclide( RelActAutoGui *gui, Wt::WContainerWidget *parent = nullptr );

  /** Function called when the user changes the nuclide via the gui; not meant to be called if loading a serialized state. */
  void handleIsotopeChange();
  
  void setFitAge( const bool do_fit );
  bool fitAge() const;
  
  void setAge( const std::string &age );
  void setAge( const Wt::WString &age );


  /** Sets the min and max age range.  
    If fit age is checked, and the nominal age is outside this range, it will be changed as well. 

    Does not emit the `updated()` or `age_changed()` signals.
   */
  void setAgeRange( Wt::WString min_age, Wt::WString max_age );
  
  
  void setNuclideEditFocus();
  
  void handleAgeChange();
  void handleAgeRangeChange();
  void validateAndCorrectAgeRange();
  
  void handleFitAgeChange();
  
  void handleColorChange();
  
  /** Will return std::monostate if invalid text entered. */
  RelActCalcAuto::SrcVariant source() const;

  const SandiaDecay::Nuclide *nuclide() const;
  const SandiaDecay::Element *element() const;
  const ReactionGamma::Reaction *reaction() const;
  std::string source_name() const;

  /** Returns the Relative Effiiciency curve index of this source.

   Throws exception if could not be deterimind (e.g., during initial model loading, the widget hasnt been added to gui, etc.), and
   there is more than one rel. eff. curve
   */
  int relEffCurveIndex() const;

  bool hasActRatioConstraint() const;
  RelActCalcAuto::RelEffCurveInput::ActRatioConstraint actRatioConstraint() const;

  bool hasMassFractionConstraint() const;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint massFractionConstraint() const;


  Wt::WColor color() const;
  Wt::WColor getColorForSource( const RelActCalcAuto::SrcVariant &source ) const;

  void setColor( const Wt::WColor &color );
  
  double age() const;
  Wt::WString ageStr() const;
  
  std::pair<Wt::WString,Wt::WString> ageRangeStr() const;
  std::pair<std::optional<double>,std::optional<double>> ageRange() const;
  
  RelActCalcAuto::NucInputInfo toNucInputInfo() const;
  
  void fromNucInputInfo( const RelActCalcAuto::NucInputInfo &info );
  
  void handleRemoveSelf();
  
  Wt::Signal<> &updated();
  
  Wt::Signal<> &remove();

  Wt::Signal<RelActAutoGuiNuclide *,bool> &fit_age_changed();
  
  Wt::Signal<RelActAutoGuiNuclide *> &age_changed();

  void addRelActRangeConstraint( std::optional<double> min_rel_act, std::optional<double> max_rel_act );
  void addActRatioConstraint( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &constraint );
  void addMassFractionConstraint( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint );
  
  void handleAddConstraint();
  void handleRemoveConstraint();
  void handleConstraintChanged();

  void updateAllowedConstraints();

  void setSummaryText( const Wt::WString &text, const Wt::WString &tooltip );
protected:
  RelActAutoGui *m_gui;

  Wt::WLineEdit *m_nuclide_edit;
  Wt::WContainerWidget *m_age_container;
  Wt::WLineEdit *m_age_edit;  
  Wt::WCheckBox *m_fit_age;
  ColorSelect *m_color_select;

  Wt::WContainerWidget *m_lower_container;
  Wt::WPushButton *m_add_constraint_btn;
  RelActAutoGuiNuclideConstraint *m_constraint;
  Wt::WContainerWidget *m_age_range_container;
  Wt::WLineEdit *m_fit_age_min_edit;
  Wt::WLineEdit *m_fit_age_max_edit;

  Wt::WText *m_summary_text;

  Wt::Signal<> m_updated;
  Wt::Signal<> m_remove;
  Wt::Signal<RelActAutoGuiNuclide *,bool> m_fit_age_changed;
  Wt::Signal<RelActAutoGuiNuclide *> m_age_changed;

  RelActCalcAuto::SrcVariant m_src_info;
};//class RelActAutoGuiNuclide


#endif //RelActAutoGuiNuclide_h
