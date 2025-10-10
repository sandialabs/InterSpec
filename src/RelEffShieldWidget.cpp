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

#include "rapidxml/rapidxml.hpp"

#include <Wt/WCheckBox>
#include <Wt/WGroupBox>
#include <Wt/WLineEdit>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WRegExpValidator>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RelEffShieldWidget.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace Wt;
using namespace std;

RelEffShieldWidget::RelEffShieldWidget( ShieldType type, Wt::WContainerWidget *parent)
  : Wt::WGroupBox( Wt::WString::tr(type == ShieldType::SelfAtten ? "resw-self-atten-title" : "resw-ext-atten-title"), parent),
  m_type(type),
  m_changed( this )
{
  WApplication *app = WApplication::instance();
  InterSpec *interspec = InterSpec::instance();
  if( !interspec || !app )
    throw std::runtime_error( "RelEffShieldWidget: InterSpec instance not found" );
    
  if( interspec )
    interspec->useMessageResourceBundle( "RelEffShieldWidget" );

  app->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
  app->useStyleSheet( "InterSpec_resources/RelEffShieldWidget.css" );

  addStyleClass("RelEffShieldWidget");
    
  m_frameSwitch = new SwitchCheckbox( Wt::WString::tr("resw-material-frame-label"), 
                                    Wt::WString::tr("resw-generic-frame-label"), this );
  m_frameSwitch->checked().connect( this, &RelEffShieldWidget::materialTypeUpdated );
  m_frameSwitch->unChecked().connect( this, &RelEffShieldWidget::materialTypeUpdated );

  m_stackedWidget = new Wt::WStackedWidget( this );

  m_materialFrame = new Wt::WContainerWidget();
  m_stackedWidget->addWidget( m_materialFrame );
  m_materialFrame->addStyleClass("MaterialFrame");
  auto materialLabel = new Wt::WLabel( Wt::WString::tr("resw-material-frame-label"), m_materialFrame); //TODO: maybe make `resw-material-frame-label` its own string
  materialLabel->addStyleClass("GridFirstCol GridFirstRow");
    
  m_materialEdit = new Wt::WLineEdit( m_materialFrame );
  m_materialEdit->addStyleClass("GridSecondCol GridFirstRow");
  m_materialSuggest = interspec->shieldingSuggester();
  if( m_materialSuggest )
    m_materialSuggest->forEdit( m_materialEdit, WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );
  m_materialEdit->changed().connect( this, &RelEffShieldWidget::materialUpdated );  
  m_materialEdit->enterPressed().connect( this, &RelEffShieldWidget::materialUpdated );

  auto thicknessLabel = new Wt::WLabel( Wt::WString::tr("resw-thickness-label"), m_materialFrame);
  thicknessLabel->addStyleClass("GridFirstCol GridSecondRow");
    

  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_thicknessEdit = new Wt::WLineEdit( m_materialFrame );
  m_thicknessEdit->setValidator( distValidator );
  m_thicknessEdit->addStyleClass("GridSecondCol GridSecondRow");
  m_thicknessEdit->changed().connect( this, &RelEffShieldWidget::userUpdated );
  m_thicknessEdit->enterPressed().connect( this, &RelEffShieldWidget::userUpdated );

  m_fitThickness = new Wt::WCheckBox( WString::tr("Fit"), m_materialFrame);
  m_fitThickness->addStyleClass("GridThirdCol GridSecondRow GridJustifyEnd");
  m_fitThickness->unChecked().connect( this, &RelEffShieldWidget::userUpdated );
  m_fitThickness->checked().connect( this, &RelEffShieldWidget::userUpdated );
  m_fitThickness->setChecked( true );

  m_parametersFrame = new Wt::WContainerWidget();
  m_stackedWidget->addWidget( m_parametersFrame );

  m_parametersFrame->addStyleClass("GenericFrame");
  m_atomicNumberLabel = new Wt::WLabel( Wt::WString::tr("resw-atomic-number-label"), m_parametersFrame);
  m_atomicNumberLabel->addStyleClass("AN GridFirstCol GridFirstRow");
    
  m_atomicNumber = new NativeFloatSpinBox( m_parametersFrame );
  m_atomicNumber->setRange( 1.0f, 98.0f ); //dont set single-step, so it wont be invalid for a non-step value
  m_atomicNumber->setWidth( 50 );
  m_atomicNumber->setSpinnerHidden( true );
  m_atomicNumber->setFormatString( "%.3f" );
  m_atomicNumber->addStyleClass("AN GridSecondCol GridFirstRow");
  m_atomicNumber->setValue( 26.0f );
  m_atomicNumber->valueChanged().connect( this, &RelEffShieldWidget::userUpdated );
    
  m_fitAtomicNumber = new Wt::WCheckBox(WString::tr("Fit"), m_parametersFrame );
  m_fitAtomicNumber->addStyleClass("AN GridThirdCol GridFirstRow GridJustifyEnd");
  m_fitAtomicNumber->unChecked().connect( this, &RelEffShieldWidget::userUpdated );
  m_fitAtomicNumber->checked().connect( this, &RelEffShieldWidget::userUpdated );
  m_fitAtomicNumber->setChecked( false );

  m_arealDensityLabel = new Wt::WLabel(WString::tr("resw-areal-density-label"), m_parametersFrame);
  m_arealDensityLabel->addStyleClass("AD GridFirstCol GridSecondRow");
    
  m_arealDensity = new NativeFloatSpinBox( m_parametersFrame );
  m_arealDensity->setRange( 0.0f, 500.0f ); //dont set single-step, so it wont be invalid for a non-step value
  m_arealDensity->setWidth( 50 );
  m_arealDensity->setSpinnerHidden( true );
  m_arealDensity->setFormatString( "%.3f" );
  m_arealDensity->setPlaceholderText( Wt::WString::tr("resw-areal-density-placeholder") );
  m_arealDensity->addStyleClass("AD GridSecondCol GridSecondRow");
  m_arealDensity->setValue( 0.0 );
  m_arealDensity->valueChanged().connect( this, &RelEffShieldWidget::userUpdated );

  m_fitArealDensity = new Wt::WCheckBox(WString::tr("Fit"), m_parametersFrame );
  m_fitArealDensity->addStyleClass("AD GridThirdCol GridSecondRow GridJustifyEnd");
  m_fitArealDensity->setChecked( false );
  m_fitArealDensity->unChecked().connect( this, &RelEffShieldWidget::userUpdated );
  m_fitArealDensity->checked().connect( this, &RelEffShieldWidget::userUpdated );
  m_fitArealDensity->setChecked( true );

  m_stackedWidget->setCurrentIndex(0);

  userUpdated();
}//RelEffShieldWidget constructor


RelEffShieldWidget::~RelEffShieldWidget()
{
  InterSpec *interspec = InterSpec::instance();
  if( !interspec )
    return;
    
  const vector<Wt::WObject *> &kids = interspec->Wt::WObject::children();
  // kids.size() is usually just 2 or 3, so this isnt that heavy of an operation.
  auto pos = std::find( begin(kids), end(kids), static_cast<WObject *>(m_materialSuggest) );
  if( pos != end(kids) )
    m_materialSuggest->removeEdit( m_materialEdit );
}//~RelEffShieldWidget()


bool RelEffShieldWidget::isMaterialSelected() const
{
  return !m_frameSwitch->isChecked();
}


void RelEffShieldWidget::setMaterialSelected( bool selected )
{
  m_stackedWidget->setCurrentIndex( selected ? 0 : 1 );
  m_frameSwitch->setChecked( !selected );
}


const Material *RelEffShieldWidget::material( const std::string &text, MaterialDB *matDB )
{
  assert( matDB );
  if( !matDB )
    return nullptr;
  
  if( text.empty() )
    return nullptr;
  
  try
  {
    const Material *answer = matDB->material( text );
    return answer;
  }catch(...)
  {
    //material wasnt in the database
  }
  
  //See if 'text' is a chemical formula, if so add it to possible suggestions
  try
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const Material *mat = matDB->parseChemicalFormula( text, db );
    
    // Update the material suggestions widget, if we are in an InterSpec
    InterSpec *interspec = InterSpec::instance();
    Wt::WSuggestionPopup *suggester = interspec ? interspec->shieldingSuggester() : nullptr;
    if( suggester )
      suggester->addSuggestion( mat->name, mat->name );
    return mat;
  }catch(...)
  {
    // No luck
  }

  return nullptr;
}//const Material *material( const std::string &text )


const Material *RelEffShieldWidget::material() const
{
  InterSpec * const viewer = InterSpec::instance();
  MaterialDB * const matDB = viewer ? viewer->materialDataBase() : nullptr;
  
  const string text = m_materialEdit->text().toUTF8();
  return material( text, matDB );
}//const Material *material() const


Wt::WString RelEffShieldWidget::materialNameTxt() const
{
  return m_materialEdit->text();
}


void RelEffShieldWidget::setMaterial(const std::string &material)
{
  m_materialEdit->setText(material);
}


double RelEffShieldWidget::thickness() const
{
  const string txt = m_thicknessEdit->text().toUTF8();
  try
  {
    return PhysicalUnits::stringToDistance( txt );
  }catch(...)
  {
  }

  return 0.0;
}


Wt::WString RelEffShieldWidget::thicknessTxt() const
{
  return m_thicknessEdit->text();
}


void RelEffShieldWidget::setThickness( double thickness )
{
  const string txt = PhysicalUnits::printToBestLengthUnits(thickness, 5);
  m_thicknessEdit->setText( txt );
}


void RelEffShieldWidget::setThickness( const Wt::WString &thickness )
{
  m_thicknessEdit->setText( thickness );
}


bool RelEffShieldWidget::fitThickness() const
{
  return m_fitThickness->isChecked();
}


void RelEffShieldWidget::setFitThickness(bool fit)
{
  m_fitThickness->setChecked( fit );
  m_thicknessEdit->setDisabled( fit );
}


double RelEffShieldWidget::atomicNumber() const
{
  return m_atomicNumber->value();
}


void RelEffShieldWidget::setAtomicNumber(double atomicNumber)
{
  m_atomicNumber->setValue(atomicNumber);
}


bool RelEffShieldWidget::fitAtomicNumber() const
{
  return m_fitAtomicNumber->isChecked();
}


void RelEffShieldWidget::setFitAtomicNumber(bool fit)
{
  m_fitAtomicNumber->setChecked( fit );
  m_atomicNumber->setDisabled( fit );
}


double RelEffShieldWidget::arealDensity() const
{
  return m_arealDensity->value();
}


void RelEffShieldWidget::setArealDensity(double arealDensity)
{
  m_arealDensity->setValue(arealDensity);
}


bool RelEffShieldWidget::fitArealDensity() const
{
  return m_fitArealDensity->isChecked();
}


void RelEffShieldWidget::setFitArealDensity(bool fit)
{
  m_fitArealDensity->setChecked( fit );
  m_arealDensity->setDisabled( fit );
}


bool RelEffShieldWidget::nonEmpty() const
{
  if( isMaterialSelected() )
  {
    const Material *mat = material();
    if( !mat )
      return false;
    
    if( fitThickness() )
      return true;

    //For the moment, even if AD is zero, we'll say its non-empty, as long as there is some text
    return !m_thicknessEdit->valueText().empty();
    //const double thick = thickness();
    //return (thick > 0.0);
  }//if( isMaterialSelected() )

  const double an = m_atomicNumber->value();
  return ( (((an >= 1.0) && (an <= 98.0)) || fitAtomicNumber())
          && (!m_arealDensity->valueText().empty() || fitArealDensity()) );
  
  //return (fitAtomicNumber() || ((m_atomicNumber->value() >= 1.0) && (m_atomicNumber->value() <= 98.0)))
  //      && (fitArealDensity() || (m_arealDensity->value() > 0.0));
}


void RelEffShieldWidget::resetState()
{
  setMaterial("");
  setThickness(0.0);
  setFitThickness(true);
  setAtomicNumber(26.0);
  setFitAtomicNumber(false);
  setArealDensity(0.0);
  setFitArealDensity(true);
}


void RelEffShieldWidget::userUpdated()
{
  m_arealDensity->setDisabled( m_fitArealDensity->isChecked() );
  m_atomicNumber->setDisabled( m_fitAtomicNumber->isChecked() );
  m_thicknessEdit->setDisabled( m_fitThickness->isChecked() );
 
  m_changed.emit();
}


void RelEffShieldWidget::materialUpdated()
{ 
  const Material * const mat = material();
  m_materialEdit->setText( mat ? mat->name : string("") );
  
  userUpdated();
}


void RelEffShieldWidget::materialTypeUpdated()
{
  m_stackedWidget->setCurrentIndex( m_frameSwitch->isChecked() ? 1 : 0 );
  
  userUpdated();
}//materialTypeUpdated()


std::unique_ptr<RelEffShieldState> RelEffShieldWidget::state() const
{
  auto s = std::make_unique<RelEffShieldState>();
  s->materialSelected = isMaterialSelected();
  s->material = materialNameTxt().toUTF8();
  s->thickness = m_thicknessEdit->text().toUTF8();
  s->fitThickness = fitThickness();
  s->atomicNumber = atomicNumber();
  s->fitAtomicNumber = fitAtomicNumber();
  s->arealDensity = arealDensity();
  s->fitArealDensity = fitArealDensity();

  return s;
}


void RelEffShieldWidget::setState(const RelEffShieldState& s)
{
  m_frameSwitch->setChecked( !s.materialSelected );
  setMaterial(s.material);
  setThickness(s.thickness);
  setFitThickness(s.fitThickness);
  setAtomicNumber(s.atomicNumber);
  setFitAtomicNumber(s.fitAtomicNumber);
  setArealDensity(s.arealDensity);
  setFitArealDensity(s.fitArealDensity);
  
  //materialTypeUpdated(); //Not calling, because we dont want to emit that user changed things
  m_stackedWidget->setCurrentIndex( m_frameSwitch->isChecked() ? 1 : 0 );
}


Wt::Signal<void> &RelEffShieldWidget::changed()
{
  return m_changed;
}


void RelEffShieldState::toXml( rapidxml::xml_node<> *node ) const
{
  using namespace rapidxml;

  xml_document<> *doc = node->document();
  assert( doc );
  xml_node<>* root = doc->allocate_node(node_element, "RelEffShield");
    
  auto allocateNode = [&](const char* name, const char* value) {
    char* nodeName = doc->allocate_string(name);
    char* nodeValue = doc->allocate_string(value);
    return doc->allocate_node(node_element, nodeName, nodeValue);
  };

  root->append_node(allocateNode("MaterialDefined", materialSelected ? "true" : "false"));
  root->append_node(allocateNode("material", material.c_str()));
  root->append_node(allocateNode("thickness", thickness.c_str()));
  root->append_node(allocateNode("fitThickness", fitThickness ? "true" : "false"));
  root->append_node(allocateNode("atomicNumber", std::to_string(atomicNumber).c_str()));
  root->append_node(allocateNode("fitAtomicNumber", fitAtomicNumber ? "true" : "false"));
  root->append_node(allocateNode("arealDensity", std::to_string(arealDensity).c_str()));
  root->append_node(allocateNode("fitArealDensity", fitArealDensity ? "true" : "false"));

  node->append_node(root);
}//void RelEffShieldState::toXml(rapidxml::xml_node<>* node)


void RelEffShieldState::fromXml(const rapidxml::xml_node<>* node)
{
  const rapidxml::xml_node<>* val = node->first_node("MaterialDefined");
  if( !val )
    throw std::runtime_error( "RelEffShieldState: missing required node 'MaterialDefined'" );
  materialSelected = XML_VALUE_ICOMPARE(val, "true");
  
  val = XML_FIRST_NODE(node, "material");
  material = SpecUtils::xml_value_str(val);
  
  val = XML_FIRST_NODE(node, "thickness");
  thickness = SpecUtils::xml_value_str(val);
  
  val = XML_FIRST_NODE(node, "fitThickness");
  fitThickness = val ? XML_VALUE_ICOMPARE(val, "true") : false;
  
  val = XML_FIRST_NODE(node, "atomicNumber");
  const string an_str = SpecUtils::xml_value_str(val);
  if( !SpecUtils::parse_double(an_str.c_str(), an_str.size(), atomicNumber) )
    atomicNumber = 26.0;
  
  val = XML_FIRST_NODE(node, "fitAtomicNumber");
  fitAtomicNumber = val ? XML_VALUE_ICOMPARE(val, "true") : false;
  
  val = XML_FIRST_NODE(node, "arealDensity");
  const string ad_str = SpecUtils::xml_value_str(val);  
  if( !SpecUtils::parse_double(ad_str.c_str(), ad_str.size(), arealDensity) )
    arealDensity = 0.0;
  
  val = XML_FIRST_NODE(node, "fitArealDensity");
  fitArealDensity = val ? XML_VALUE_ICOMPARE(val, "true") : false;
}//void RelEffShieldState::fromXml(const rapidxml::xml_node<>* node)


std::shared_ptr<RelActCalc::PhysicalModelShieldInput> RelEffShieldState::fitInput( MaterialDB *materialDB ) const
{
  auto self_atten = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
  if( materialSelected )
  {
    const Material * const mat = RelEffShieldWidget::material( material, materialDB );
    assert( mat );
    if( !mat )
      return nullptr;
    
    self_atten->material = make_shared<Material>( *mat );
    const double distance = PhysicalUnits::stringToDistance(thickness);
    self_atten->areal_density = distance * mat->density;
    
    cout << "RelEffShieldState::fitInput: dist=" << distance/PhysicalUnits::cm << " cm, AD="
    << self_atten->areal_density / PhysicalUnits::g_per_cm2 << " g/cm2, for '" << mat->name << "'"
    << " wit hdensity " << self_atten->material->density * PhysicalUnits::cm3 / PhysicalUnits::g << " g/cm3."
    << endl;
    
    // The rest of this scope is more to remind me in the future
    self_atten->fit_atomic_number = false;
    self_atten->fit_areal_density = fitThickness;
    
    if( fitThickness )
    {
      //RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2 * PhysicalUnits::g_per_cm2;
      self_atten->lower_fit_areal_density = 0.0;
      self_atten->upper_fit_areal_density = 0.0;
    }
    
    self_atten->atomic_number = 0.0;
  }else
  {
    self_atten->atomic_number = atomicNumber;
    self_atten->areal_density = arealDensity * PhysicalUnits::g_per_cm2;
    
    self_atten->fit_atomic_number = fitAtomicNumber;
    self_atten->fit_areal_density = fitArealDensity;
    
    // The rest of this scope is more to remind me in the future
    self_atten->lower_fit_atomic_number = 1.0;
    self_atten->upper_fit_atomic_number = 98.0;
    
    if( fitArealDensity )
    {
      self_atten->lower_fit_areal_density = 0.0;
      self_atten->upper_fit_areal_density = 0.0;
    }
  }//if( materialSelected ) / else
  
  return self_atten;
}//fitInput()



void RelEffShieldState::setStateFromFitInput( const RelActCalc::PhysicalModelShieldInput &input )
{
  if( input.material )
  {
    materialSelected = true;
    material = input.material->name;
    const double dist = input.areal_density / input.material->density;
    thickness = PhysicalUnits::printToBestLengthUnits( dist, 6 );
  }else
  {
    materialSelected = false;
    material.clear();
    thickness = "";
  }
  
  atomicNumber = input.atomic_number;
  arealDensity = input.areal_density / PhysicalUnits::g_per_cm2;
  fitArealDensity = input.fit_areal_density;
  fitAtomicNumber = input.fit_atomic_number;
  fitThickness = input.fit_areal_density;
}//void setStateFromFitInput( const RelActCalc::PhysicalModelShieldInput &input )


std::shared_ptr<RelActCalc::PhysicalModelShieldInput> RelEffShieldWidget::fitInput() const
{
  if( !nonEmpty() )
    return nullptr;
  
  auto s = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
  
  if( isMaterialSelected() )
  {
    const Material * const mat = material();
    assert( mat );
    if( !mat )
      return nullptr;

    s->material = make_shared<Material>(*mat);
    s->atomic_number = 0.0;
    s->fit_atomic_number = false;
    s->fit_areal_density = fitThickness();

    if( s->material && !s->fit_areal_density )
    {
      const double thick = thickness();
      if( thick < 0.0 )
        s->areal_density = 0.0;
      else
        s->areal_density = thick * s->material->density;
    }else
    {
      s->areal_density = 1.25*PhysicalUnits::g_per_cm2; //Arbitrary default value, just to make sure it's not zero, btu also not to go wild incase last fit value was big
    }
  }else
  {
    const double an = atomicNumber();
    assert( an >= 0.9 && an <= 98.1 );
    if( an < 0.9 || an > 98.1 )
      return nullptr;

    s->atomic_number = an;
    s->fit_atomic_number = fitAtomicNumber();

    s->fit_areal_density = fitArealDensity();
    if( s->fit_areal_density )
      s->areal_density = 1.25*PhysicalUnits::g_per_cm2;
    else
      s->areal_density = arealDensity()*PhysicalUnits::g_per_cm2;
  }//if( isMaterialSelected() ) / else

  return s;
}//shared_ptr<RelActCalc::PhysicalModelShieldInput> fitInput() const
