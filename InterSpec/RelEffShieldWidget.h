#ifndef RelEffShieldWidget_H
#define RelEffShieldWidget_H

#include "InterSpec_config.h"

#include <Wt/WGroupBox>


#include <Wt/WContainerWidget>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WStackedWidget>
#include <Wt/WLabel>
#include <Wt/WSuggestionPopup>


//Forward declarations
namespace Wt
{ 
  class WLabel;
  class WCheckBox;
  class WLineEdit;
  class WStackedWidget;
  class WSuggestionPopup;
}//namespace Wt

struct Material;
class MaterialDB;
class SwitchCheckbox;
class NativeFloatSpinBox;

namespace RelActCalc
{
  struct PhysicalModelShieldInput;
}

struct RelEffShieldState 
{
  bool materialSelected;
  std::string material;
  std::string thickness;
  bool fitThickness;
  double atomicNumber;
  bool fitAtomicNumber;
  double arealDensity;  
  bool fitArealDensity;

  void toXml(rapidxml::xml_node<>* node) const;
  void fromXml(const rapidxml::xml_node<>* node);
};//struct RelEffShieldState

/** @brief A GUI widget roughly cooresponding to a RelActCalc::PhysicalModelShieldInput.
 * 
 * Currently, the widget is not fully complete, and does not allow for all the options 
 * in the PhysicalModelShieldInput, like limiting range of fit AN and AD.
 */
class RelEffShieldWidget : public Wt::WGroupBox
{
public:

  enum class ShieldType { SelfAtten, ExternalAtten };

  RelEffShieldWidget( ShieldType type, Wt::WContainerWidget *parent = nullptr );
  ~RelEffShieldWidget();

  bool isMaterialSelected() const;
  void setMaterialSelected( bool selected );
  const Material *material() const;
  /** The text showing in the material name input field - may not be a valid material name. */
  Wt::WString materialNameTxt() const;

  void setMaterial(const std::string &material);
  double thickness() const;
  Wt::WString thicknessTxt() const;
  void setThickness(double thickness);
  void setThickness( const Wt::WString &thickness );
  bool fitThickness() const;
  void setFitThickness(bool fit);
  double atomicNumber() const;
  void setAtomicNumber(double atomicNumber);
  bool fitAtomicNumber() const;
  void setFitAtomicNumber(bool fit);
  double arealDensity() const;
  void setArealDensity(double arealDensity);
  bool fitArealDensity() const;
  void setFitArealDensity(bool fit);

  bool nonEmpty() const;
  void resetState();

  /** For input to Rel Act shield fit.
   * 
   * If the widget is empty, returns nullptr.
   */
  std::shared_ptr<RelActCalc::PhysicalModelShieldInput> fitInput() const;

  std::unique_ptr<RelEffShieldState> state() const;
  void setState(const RelEffShieldState& state);

  Wt::Signal<void> &changed();

private:
  ShieldType m_type;
  Wt::Signal<void> m_changed;
  SwitchCheckbox *m_frameSwitch;
  Wt::WStackedWidget *m_stackedWidget;
  Wt::WContainerWidget *m_materialFrame;
  Wt::WLineEdit *m_materialEdit;
  Wt::WSuggestionPopup *m_materialSuggest;
  Wt::WLineEdit *m_thicknessEdit;
  Wt::WCheckBox *m_fitThickness;
  Wt::WContainerWidget *m_parametersFrame;
  Wt::WLabel *m_atomicNumberLabel;
  NativeFloatSpinBox *m_atomicNumber;
  Wt::WCheckBox *m_fitAtomicNumber;
  Wt::WLabel *m_arealDensityLabel;
  NativeFloatSpinBox *m_arealDensity;
  Wt::WCheckBox *m_fitArealDensity;

  void userUpdated();
  void materialUpdated();
  void materialTypeUpdated();
};//class RelEffShieldWidget

#endif // RelEffShieldWidget_H