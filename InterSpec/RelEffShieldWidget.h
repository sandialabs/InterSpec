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

namespace rapidxml
{
  template<class Ch> class xml_node;
}//namespace rapidxml

namespace RelActCalc
{
  struct PhysicalModelShieldInput;
}

/** TODO: This class should probably be replaced by using only `RelActCalc::PhysicalModelShieldInput`, or maybe inherit from it, so this way the strings can be preserved. */
struct RelEffShieldState
{
  bool materialSelected;
  std::string material;
  std::string thickness;
  bool fitThickness;
  double atomicNumber;
  bool fitAtomicNumber;
  /** This arealDensity is not stored in units of `PhysicalUnits::g_per_cm2` - so you need to multiply by
   `PhysicalUnits::g_per_cm2` before using it anywhere, other than setting the GUI.
   */
  double arealDensity;
  bool fitArealDensity;
  /** Whether the optional weak prior pulling this shield's fit areal density toward 0 is enabled.
   Only meaningful for RelActCalcAuto physical-model fits (ignored by the Manual tool). */
  bool biasArealDensity = false;

  void toXml(rapidxml::xml_node<char>* node) const;
  void fromXml(const rapidxml::xml_node<char>* node);
  
  /** Will return nullptr if invalid or empty state, including the material name can not be parsed. */
  std::shared_ptr<RelActCalc::PhysicalModelShieldInput> fitInput() const;
  
  void setStateFromFitInput( const RelActCalc::PhysicalModelShieldInput &input );
};//struct RelEffShieldState

/** @brief A GUI widget roughly corresponding to a RelActCalc::PhysicalModelShieldInput.
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
  static std::shared_ptr<const Material> material( const std::string &mat_name );
  std::shared_ptr<const Material> material() const;
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

  /** Show/hide the optional "Bias AD" checkbox (hidden by default; the RelActCalcAuto GUI enables it). */
  void setArealDensityBiasVisible( bool visible );
  bool biasArealDensity() const;
  void setBiasArealDensity( bool bias );

  bool nonEmpty() const;
  /// Resets material/AN/AD and "Fit" checkboxes
  void resetState();
  /// Resets material/AN/AD, but not "Fit" checkboxes (incase you dont want to mess with users preferences)
  void resetMaterialEntryState();

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
  /** Optional "Bias AD" checkbox (own row at the bottom, floated right); hidden unless
   `setArealDensityBiasVisible(true)` is called (RelActCalcAuto only). */
  Wt::WCheckBox *m_biasArealDensity;

  void userUpdated();
  void materialUpdated();
  void materialTypeUpdated();
};//class RelEffShieldWidget

#endif // RelEffShieldWidget_H
