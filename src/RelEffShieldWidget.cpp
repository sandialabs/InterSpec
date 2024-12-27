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
  auto materialLabel = new Wt::WLabel( Wt::WString::tr("resw-material-label"), m_materialFrame);
  materialLabel->addStyleClass("GridFirstCol GridFirstRow");
    
  m_materialEdit = new Wt::WLineEdit( m_materialFrame );
  m_materialEdit->addStyleClass("GridSecondCol GridFirstRow");
  m_materialSuggest = interspec->shieldingSuggester();
  if( m_materialSuggest )
    m_materialSuggest->forEdit( m_materialEdit );
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
  m_atomicNumber->setRange( 1.0f, 98.0f );
  m_atomicNumber->setSingleStep( 1.0f );
  m_atomicNumber->setWidth( 50 );
  m_atomicNumber->setSpinnerHidden( true );
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
  m_arealDensity->setRange( 0.0f, 500.0f );
  m_arealDensity->setSingleStep( 1.0f );
  m_arealDensity->setWidth( 50 );
  m_arealDensity->setSpinnerHidden( true );
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

bool RelEffShieldWidget::isMaterialSelected() const
{
  return m_stackedWidget->currentIndex() == 0;
}

const Material *RelEffShieldWidget::material() const
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return nullptr;

  MaterialDB *matDB = interspec->materialDataBase();
  assert( matDB );
  if( !matDB )
    return nullptr;
  
  const string text = m_materialEdit->text().toUTF8();
  
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
    if( m_materialSuggest )
      m_materialSuggest->addSuggestion( mat->name, mat->name );
    return mat;
  }catch(...)
  {
    // No luck
  }

  return nullptr;
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
  return std::stod(m_thicknessEdit->text().toUTF8());
}

void RelEffShieldWidget::setThickness(double thickness)
{
  m_thicknessEdit->setText(std::to_string(thickness));
}

bool RelEffShieldWidget::fitThickness() const
{
  return m_fitThickness->isChecked();
}

void RelEffShieldWidget::setFitThickness(bool fit)
{
  m_fitThickness->setChecked(fit);
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
  m_fitAtomicNumber->setChecked(fit);
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
  m_fitArealDensity->setChecked(fit);
}

void RelEffShieldWidget::userUpdated()
{
  m_arealDensity->setDisabled( !m_fitArealDensity->isChecked() );
  m_atomicNumber->setDisabled( !m_fitAtomicNumber->isChecked() );
  m_thicknessEdit->setDisabled( !m_fitThickness->isChecked() );
 
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
if (m_frameSwitch->isChecked())
    m_stackedWidget->setCurrentIndex(1);
  else
    m_stackedWidget->setCurrentIndex(0);

  userUpdated();
}//materialTypeUpdated()


std::unique_ptr<RelEffShieldState> RelEffShieldWidget::state() const
{
  auto s = std::make_unique<RelEffShieldState>();
  s->materialSelected = isMaterialSelected();
  s->material = materialNameTxt().toUTF8();
  s->thickness = thickness();
  s->fitThickness = fitThickness();
  s->atomicNumber = atomicNumber();
  s->fitAtomicNumber = fitAtomicNumber();
  s->arealDensity = arealDensity();
  s->fitArealDensity = fitArealDensity();

  return s;
}

void RelEffShieldWidget::setState(const RelEffShieldState& s)
{
  m_frameSwitch->setChecked(s.materialSelected);
  setMaterial(s.material);
  setThickness(s.thickness);
  setFitThickness(s.fitThickness);
  setAtomicNumber(s.atomicNumber);
  setFitAtomicNumber(s.fitAtomicNumber);
  setArealDensity(s.arealDensity);
  setFitArealDensity(s.fitArealDensity);
  
  materialTypeUpdated();
}

Wt::Signal<void> &RelEffShieldWidget::changed()
{
  return m_changed;
}

void RelEffShieldState::toXml( rapidxml::xml_node<> *node ) const
{
  using namespace rapidxml;

  xml_document<> *doc = node->document();
  xml_node<>* root = doc->allocate_node(node_element, "RelEffShieldState");
    
  auto allocateNode = [&](const char* name, const char* value) {
    char* nodeName = doc->allocate_string(name);
    char* nodeValue = doc->allocate_string(value);
    return doc->allocate_node(node_element, nodeName, nodeValue);
  };

  root->append_node(allocateNode("materialSelected", materialSelected ? "true" : "false"));
  root->append_node(allocateNode("material", material.c_str()));
  root->append_node(allocateNode("thickness", std::to_string(thickness).c_str()));
  root->append_node(allocateNode("fitThickness", fitThickness ? "true" : "false"));
  root->append_node(allocateNode("atomicNumber", std::to_string(atomicNumber).c_str()));
  root->append_node(allocateNode("fitAtomicNumber", fitAtomicNumber ? "true" : "false"));
  root->append_node(allocateNode("arealDensity", std::to_string(arealDensity).c_str()));
  root->append_node(allocateNode("fitArealDensity", fitArealDensity ? "true" : "false"));

  node->append_node(root);
}

void RelEffShieldState::fromXml(const rapidxml::xml_node<>* node)
{
  const rapidxml::xml_node<>* val = node->first_node("materialSelected");
  if( !val )
    throw std::runtime_error( "Missing node: materialSelected" );
  materialSelected = XML_VALUE_ICOMPARE(val, "true");
  
  val = node->first_node("material");
  if( !val )
    throw std::runtime_error( "Missing node: material" );
  material = SpecUtils::xml_value_str(val);
  
  val = node->first_node("thickness");
  if( !val )
    throw std::runtime_error( "Missing node: thickness" );
  thickness = std::stod(SpecUtils::xml_value_str(val));
  
  val = node->first_node("fitThickness");
  if( !val )
    throw std::runtime_error( "Missing node: fitThickness" );
  fitThickness = XML_VALUE_ICOMPARE(val, "true");
  
  val = node->first_node("atomicNumber");
  if( !val )
    throw std::runtime_error( "Missing node: atomicNumber" );
  atomicNumber = std::stod(SpecUtils::xml_value_str(val));
  
  val = node->first_node("fitAtomicNumber");
  if( !val )
    throw std::runtime_error( "Missing node: fitAtomicNumber" );
  fitAtomicNumber = XML_VALUE_ICOMPARE(val, "true");
  
  val = node->first_node("arealDensity");
  if( !val )
    throw std::runtime_error( "Missing node: arealDensity" );
  arealDensity = std::stod(SpecUtils::xml_value_str(val));
  
  val = node->first_node("fitArealDensity");
  if( !val )
    throw std::runtime_error( "Missing node: fitArealDensity" );
  val = node->first_node("fitArealDensity");
  fitArealDensity = XML_VALUE_ICOMPARE(val, "true");
}

std::shared_ptr<RelActCalc::PhysicalModelShieldInput> RelEffShieldWidget::fitInput() const
{
  auto s = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
  
  if( isMaterialSelected() )
  {
    const Material * const mat = material();
    s->material = mat ? make_shared<Material>(*mat) : nullptr;
    s->fit_areal_density = fitThickness();

    if( s->material && !s->fit_areal_density )
      s->areal_density = thickness() * s->material->density;
    else
      s->areal_density = 1.25*PhysicalUnits::g_per_cm2; //Arbitrary default value, just to make sure it's not zero, btu also not to go wild incase last fit value was big
  }else
  {
    s->fit_areal_density = fitArealDensity();
    if( s->fit_areal_density )
      s->areal_density = 1.25*PhysicalUnits::g_per_cm2;
    else
      s->areal_density = arealDensity()*PhysicalUnits::g_per_cm2;
    
    s->atomic_number = atomicNumber();
    s->fit_atomic_number = fitAtomicNumber();
  }//if( isMaterialSelected() ) / else

  return s;
}//shared_ptr<RelActCalc::PhysicalModelShieldInput> fitInput() const
