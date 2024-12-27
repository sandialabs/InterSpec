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


class RelEffShieldWidget : public Wt::WGroupBox
{
public:
  RelEffShieldWidget( const Wt::WString &title, Wt::WContainerWidget *parent = nullptr );

  bool isMaterialSelected() const;
  std::string material() const;


//  std::shared_ptr<const Material> material();
//MaterialDB *m_materialDB;

  void setMaterial(const std::string &material);
  double thickness() const;
  void setThickness(double thickness);
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

private:
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

  void onFrameSwitchChanged();
};//class RelEffShieldWidget

#endif // RelEffShieldWidget_H