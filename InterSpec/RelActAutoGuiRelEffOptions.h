#ifndef RelActAutoGuiRelEffOptions_h
#define RelActAutoGuiRelEffOptions_h
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

#include <vector>
#include <memory>
#include <optional>

#include <Wt/WString>
#include <Wt/WContainerWidget>



#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcAuto.h"

namespace Wt
{
  class WMenu;
  class WText;
  class WCheckBox;
  class WComboBox;
  class WInPlaceEdit;
}//namespace Wt

class RelEffShieldWidget;
class RelActAutoGui;

class RelActAutoGuiRelEffOptions : public Wt::WContainerWidget
{
public:
  RelActAutoGuiRelEffOptions(RelActAutoGui *gui, Wt::WString name, Wt::WContainerWidget *parent = nullptr);
  ~RelActAutoGuiRelEffOptions();

  void showAndHideOptionsForEqnType();
  void initPhysModelShields();
  void setIsOnlyOneRelEffCurve(const bool is_only_rel_eff_curve);
  void setHasMultiplePhysicalModels(const bool has_multiple_phys_models);
  Wt::WString name() const;
  void setName(const Wt::WString &name);
  void setRelEffCurveInput(const RelActCalcAuto::RelEffCurveInput &rel_eff);
  void updatePuCorrelationOptions(const std::vector<RelActCalcAuto::NucInputInfo> &nuclides);
  void setEqnTxt( const Wt::WString &txt );

  Wt::Signal<RelActAutoGuiRelEffOptions *> &addRelEffCurve();
  Wt::Signal<RelActAutoGuiRelEffOptions *> &delRelEffCurve();
  Wt::Signal<RelActAutoGuiRelEffOptions *, Wt::WString> &nameChanged();
  Wt::Signal<RelActAutoGuiRelEffOptions *> &equationTypeChanged();
  Wt::Signal<RelActAutoGuiRelEffOptions *> &sameHoerlOnAllCurvesChanged();
  Wt::Signal<RelActAutoGuiRelEffOptions *> &sameExternalShieldingChanged();
  Wt::Signal<RelActAutoGuiRelEffOptions *> &optionsChanged();

  RelActCalc::RelEffEqnForm rel_eff_eqn_form() const;
  size_t rel_eff_eqn_order() const;
  std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> phys_model_self_atten() const;
  std::vector<std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>> phys_model_external_atten() const;
  bool phys_model_use_hoerl() const;
  void setPhysModelUseHoerl(const bool use_hoerl);
  bool physModelSameHoerlOnAllCurves() const;
  void setPhysModelSameHoerlOnAllCurves(const bool same_hoerl_all_curves);
  bool physModelSameExtShieldAllCurves() const;
  void setPhysModelSameExtShieldAllCurves(const bool same_ext_shield_all_curves);
  RelActCalc::PuCorrMethod pu242_correlation_method() const;

  void update_self_atten_shield_widget( const std::optional<RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo> &self_atten );
  void update_external_atten_shield_widget( const std::vector<RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo> &ext_shields);

  void update_self_atten_shield_widget( const RelEffShieldWidget *shield );
  void update_external_atten_shield_widget( const std::vector<const RelEffShieldWidget *> &ext_shields );
  
  const RelEffShieldWidget *selfAttenWidget() const;
  std::vector<const RelEffShieldWidget *> externalAttenWidgets() const;
  
protected:
  void update_shield_widget(const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo &shield,
                              RelEffShieldWidget *w);
  
  void handleRelEffEqnTypeChanged();
  void handleSameHoerlOnAllCurvesChanged();
  void handleSameExternalShieldingChanged();
  
protected:
  RelActAutoGui *const m_gui;

  Wt::WComboBox *m_rel_eff_eqn_form;
  Wt::WContainerWidget *m_eqn_order_div;
  Wt::WLabel *m_rel_eff_eqn_order_label;
  Wt::WComboBox *m_rel_eff_eqn_order;
  Wt::WInPlaceEdit *m_rel_eff_curve_name;
  /** This variable should always match the visibility state of `m_pu_corr_div`, but is necassary since we cant use
   the visibility state of the widget to track validity of this settings, because `isVisisble()` will return false when we
   are saving the state while closiing the RelActAuto tool.
   */
  bool m_pu_corr_enabled;
  Wt::WContainerWidget *m_pu_corr_div;
  Wt::WComboBox *m_pu_corr_method;
  Wt::WContainerWidget *m_phys_model_opts;
  Wt::WContainerWidget *m_phys_model_shields;
  RelEffShieldWidget *m_phys_model_self_atten;
  Wt::WContainerWidget *m_phys_ext_attens;
  Wt::WCheckBox *m_phys_model_use_hoerl;
  Wt::WCheckBox *m_phys_model_same_hoerl_on_all_curves;
  Wt::WCheckBox *m_phys_model_same_ext_shield_all_curves;
  Wt::WText *m_eqn_txt;
  Wt::WContainerWidget *m_add_del_rel_eff_div;
  Wt::WPushButton *m_add_rel_eff_btn;
  Wt::WPushButton *m_del_rel_eff_btn;
  Wt::Signal<RelActAutoGuiRelEffOptions *> m_add_rel_eff_curve_signal;
  Wt::Signal<RelActAutoGuiRelEffOptions *> m_del_rel_eff_curve_signal;
  Wt::Signal<RelActAutoGuiRelEffOptions *, Wt::WString> m_name_changed_signal;
  Wt::Signal<RelActAutoGuiRelEffOptions *> m_eqn_form_changed;
  Wt::Signal<RelActAutoGuiRelEffOptions *> m_same_hoerl_on_all_curves;
  Wt::Signal<RelActAutoGuiRelEffOptions *> m_same_ext_shield_on_all_curves;
  

  Wt::Signal<RelActAutoGuiRelEffOptions *> m_options_changed_signal;
  
  bool m_has_multiple_phys_models;

  void emitAddRelEffCurve();
  void emitDelRelEffCurve();
  void emitNameChanged();
  void emitOptionsChanged();
};//class RelActAutoGuiRelEffOptions

#endif // RelActAutoGuiRelEffOptions_h 
