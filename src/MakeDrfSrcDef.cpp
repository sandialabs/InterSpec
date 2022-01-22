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

#include <Wt/WDate>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WComboBox>
#include <Wt/WDateEdit>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WTableCell>
#include <Wt/WApplication>
#include <Wt/WDoubleSpinBox>
#include <Wt/WDoubleValidator>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>
#include <Wt/WSuggestionPopup>

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/IsotopeSelectionAids.h"
//#include "InterSpec/IsotopeNameFilterModel.h"

using namespace std;
using namespace Wt;

namespace
{
  const int sm_distance_row        = 1;
  const int sm_activity_row        = 2;
  const int sm_activity_uncert_row = 3;
  const int sm_assay_date_row      = 4;
  const int sm_spec_date_row       = 5;
  const int sm_age_at_assay_row    = 6;
  const int sm_decayed_info_row    = 7;
  const int sm_shield_material_row = 8;
  const int sm_options_row         = 9;
}//namespace


MakeDrfSrcDef::MakeDrfSrcDef( const SandiaDecay::Nuclide *nuc,
              const boost::posix_time::ptime &measDate,
              MaterialDB *materialDB,
              Wt::WSuggestionPopup *materialSuggest,
              Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_table( nullptr ),
  m_nuclide( nuc ),
  m_materialDB( materialDB ),
  m_materialSuggest( materialSuggest ),
  m_nuclideLabel( nullptr ),
  m_distanceEdit( nullptr ),
  m_activityEdit( nullptr ),
  m_activityUncertainty( nullptr ),
  m_useAgeInfo( nullptr ),
  m_assayDate( nullptr ),
  m_drfMeasurementDate( nullptr ),
  m_sourceInfoAtMeasurement( nullptr ),
  m_sourceAgeAtAssay( nullptr ),
  m_useShielding( nullptr ),
  m_shieldingSelect( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/MakeDrfSrcDef.css" );
  
  addStyleClass( "MakeDrfSrcDef" );
  
  create();
  
  setNuclide( m_nuclide );
  
  if( !measDate.is_special() )
  {
    m_drfMeasurementDate->setDate( WDateTime::fromPosixTime(measDate).date() );
    validateDateFields();
  }
}//MakeDrfSrcDef constructor


MakeDrfSrcDef::~MakeDrfSrcDef()
{
}


Wt::Signal<> &MakeDrfSrcDef::updated()
{
  return m_updated;
}

void MakeDrfSrcDef::setNuclide( const SandiaDecay::Nuclide *nuc )
{
  m_nuclide = nuc;
  const bool notMuchEvolution = (!nuc || PeakDef::ageFitNotAllowed(nuc));
  
  m_table->rowAt(sm_age_at_assay_row)->setHidden( notMuchEvolution );
  
  if( notMuchEvolution || !nuc )
  {
    m_sourceAgeAtAssay->setText( "0s" );
  }else
  {
    double ageval = 5.0*nuc->halfLife;
    if( nuc->canObtainPromptEquilibrium() )
      ageval = log(10000.0)/log(2.0) * nuc->promptEquilibriumHalfLife();
    if( ageval > 20*PhysicalUnits::year )
      ageval = 20*PhysicalUnits::year;
    m_sourceAgeAtAssay->setText( PhysicalUnits::printToBestTimeUnits(ageval) );
  }
  
  if( nuc )
  {
    const string label = "<span class=\"SrcTitleNuc\">" + nuc->symbol + "</span>,"
                         " <span class=\"SrcTitleHl\">&lambda;<sub>&frac12;</sub>="
                         + PhysicalUnits::printToBestTimeUnits(nuc->halfLife) + "</span>";
    m_nuclideLabel->setText( WString::fromUTF8(label) );
    m_useAgeInfo->show();
  }else
  {
    m_nuclideLabel->setText( "Non-specified Nuclide" );
    m_useAgeInfo->setUnChecked();
    m_distanceEdit->setValueText( "25 cm" );
    m_activityEdit->setValueText( "1 uCi" );
    m_useAgeInfo->hide();
  }
  
  useAgeInfoUserToggled();
}//setNuclide(...)


void MakeDrfSrcDef::create()
{
  m_table = new WTable( this );
  m_table->addStyleClass( "SrcInputTable" );
  
  WTableCell *cell = m_table->elementAt(0,0);
  cell->setColumnSpan( 2 );
  cell->addStyleClass( "SrcNuclideTitle" );
  m_nuclideLabel = new WText( cell );
  
/*
  //Code to put a nuclide suggestion into a WLineEdit so user could select nuclide.
  string replacerJs, matcherJs;
  PhotopeakDelegate::EditWidget::replacerJs( replacerJs );
  PhotopeakDelegate::EditWidget::nuclideNameMatcherJs( matcherJs );
  
  WSuggestionPopup *nucSuggest = new WSuggestionPopup( matcherJs, replacerJs );
  nucSuggest->setJavaScriptMember("wtNoReparent", "true");
  
  nucSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  nucSuggest->forEdit( m_nuclideEdit,
                       WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );
  
  IsotopeNameFilterModel *filterModel = new IsotopeNameFilterModel( this );
  
  filterModel->excludeNuclides( false );
  filterModel->excludeXrays( true );
  filterModel->excludeEscapes( true );
  filterModel->excludeReactions( true );
  
  filterModel->filter( "" );
  nucSuggest->setFilterLength( -1 );
  nucSuggest->setModel( filterModel );
  nucSuggest->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );
  nucSuggest->filterModel().connect( filterModel, &IsotopeNameFilterModel::filter );
*/
  
  cell = m_table->elementAt(sm_distance_row,0);
  WLabel *label = new WLabel( "Distance", cell );
  cell = m_table->elementAt(sm_distance_row,1);
  m_distanceEdit = new WLineEdit( cell );
  m_distanceEdit->setTextSize( 16 );
  m_distanceEdit->setAutoComplete( false );
  label->setBuddy( m_distanceEdit );
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_distanceEdit->setValidator( distValidator );
  m_distanceEdit->setText( "50 cm" );
  m_distanceEdit->changed().connect( this, &MakeDrfSrcDef::handleUserChangedDistance );
  m_distanceEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedDistance );
  
  
  cell = m_table->elementAt(sm_activity_row,0);
  label = new WLabel( "Activity", cell );
  
  cell = m_table->elementAt(sm_activity_row,1);
  m_activityEdit = new WLineEdit( cell );
  m_activityEdit->setAutoComplete( false );
  m_activityEdit->setTextSize( 16 );
  label->setBuddy( m_activityEdit );
  
  WRegExpValidator *val = new WRegExpValidator( PhysicalUnits::sm_activityRegex, this );
  val->setFlags( Wt::MatchCaseInsensitive );
  m_activityEdit->setValidator( val );
  m_activityEdit->setText( "100 uCi" );
  m_activityEdit->changed().connect( this, &MakeDrfSrcDef::handleUserChangedActivity );
  m_activityEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedActivity );

  cell = m_table->elementAt(sm_activity_uncert_row,0);
  label = new WLabel( "Act. Uncert.&nbsp;", cell );  //The nbsp is to make this the longest label so when acti ity or shielding is shown, the width doesnt get changed
  cell = m_table->elementAt(sm_activity_uncert_row,1);
  m_activityUncertainty = new WDoubleSpinBox( cell );
  m_activityUncertainty->setValue( 0.0 );
  m_activityUncertainty->setTextSize( 14 );
  m_activityUncertainty->setAutoComplete( false );
  m_activityUncertainty->setSuffix( " %" );
  label->setBuddy( m_activityUncertainty );
  
  m_activityUncertainty->setDecimals( 1 );
  m_activityUncertainty->setRange( 0.0, 100.0 );
  m_activityUncertainty->setSingleStep( 1.0 );
  
  m_activityUncertainty->changed().connect( this, &MakeDrfSrcDef::handleUserChangedActivityUncertainty );
  m_activityUncertainty->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedActivityUncertainty );
  
  //WDoubleValidator *percentVal = new WDoubleValidator( this );
  //percentVal->setRange( 0.0, 100.0 );
  //m_activityUncertainty->setValidator( percentVal );
  
  cell = m_table->elementAt(sm_assay_date_row,0);
  label = new WLabel( "Assay Date", cell );
  cell = m_table->elementAt(sm_assay_date_row,1);
  m_assayDate = new WDateEdit( cell );
  //m_assayDate->setTextSize( 9 );
  label->setBuddy( m_assayDate );
  m_assayDate->changed().connect( this, &MakeDrfSrcDef::handleEnteredDatesUpdated );
  
  cell = m_table->elementAt(sm_spec_date_row,0);
  label = new WLabel( "Spec. Date", cell );
  cell = m_table->elementAt(sm_spec_date_row,1);
  m_drfMeasurementDate = new WDateEdit( cell );
  //m_drfMeasurementDate->setTextSize( 10 );
  //The right padding is 40px, could reduce down to 30.
  label->setBuddy( m_drfMeasurementDate );
  m_drfMeasurementDate->changed().connect( this, &MakeDrfSrcDef::handleEnteredDatesUpdated );
  
  
  cell = m_table->elementAt(sm_age_at_assay_row,0);
  label = new WLabel( "Age@Assay", cell );
  cell = m_table->elementAt(sm_age_at_assay_row,1);
  m_sourceAgeAtAssay = new WLineEdit( cell );
  m_sourceAgeAtAssay->setAutoComplete( false );
  label->setBuddy( m_sourceAgeAtAssay );
  m_sourceAgeAtAssay->changed().connect( this, &MakeDrfSrcDef::handleUserChangedAgeAtAssay );
  m_sourceAgeAtAssay->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedAgeAtAssay );
  
  val = new WRegExpValidator( PhysicalUnits::sm_timeDurationRegex, this );
  val->setFlags( Wt::MatchCaseInsensitive );
  m_sourceAgeAtAssay->setValidator( val );
  m_sourceAgeAtAssay->setText( "0s" );
  
  
  cell = m_table->elementAt(sm_decayed_info_row,0);
  label = new WLabel( "Aging Res.", cell );
  cell = m_table->elementAt(sm_decayed_info_row,1);
  m_sourceInfoAtMeasurement = new WText( cell );
  
  //Wt::WDateEdit *m_sourceCreationDate;
  
  cell = m_table->elementAt(sm_options_row,0);
  m_useAgeInfo = new WCheckBox( "Age?", cell );
  m_useAgeInfo->setFloatSide( Wt::Right );
  m_useAgeInfo->setChecked( false );
  m_useAgeInfo->checked().connect( this, &MakeDrfSrcDef::useAgeInfoUserToggled );
  m_useAgeInfo->unChecked().connect( this, &MakeDrfSrcDef::useAgeInfoUserToggled );
  
  cell = m_table->elementAt(sm_options_row,1);
  m_useShielding = new WCheckBox( "Shielded?", cell );
  m_useShielding->setFloatSide( Wt::Right );
  m_useShielding->setChecked( false );
  m_useShielding->checked().connect( this, &MakeDrfSrcDef::useShieldingInfoUserToggled );
  m_useShielding->unChecked().connect( this, &MakeDrfSrcDef::useShieldingInfoUserToggled );
  
  
  cell = m_table->elementAt(sm_shield_material_row,0);
  cell->setColumnSpan( 2 );
  
  m_shieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest, cell );
  m_shieldingSelect->materialModified().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->materialChanged().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->addingIsotopeAsSource().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->removingIsotopeAsSource().connect( this, &MakeDrfSrcDef::handleUserChangedShielding );
  m_shieldingSelect->hide();
  
  useAgeInfoUserToggled();
}//void create()


void MakeDrfSrcDef::useAgeInfoUserToggled()
{
  const bool useAge = m_useAgeInfo->isChecked();
  m_activityEdit->label()->setText( useAge ? "Assay Act." : "Activity" );
  
  m_table->rowAt(sm_assay_date_row)->setHidden( !useAge );
  m_table->rowAt(sm_spec_date_row)->setHidden( !useAge );
  
  const bool notMuchEvolution = (!m_nuclide || PeakDef::ageFitNotAllowed(m_nuclide));
  m_table->rowAt(sm_age_at_assay_row)->setHidden( !useAge || notMuchEvolution );
  
  m_table->rowAt(sm_decayed_info_row)->setHidden( !useAge );
  
  m_assayDate->setDisabled( !useAge );
  m_drfMeasurementDate->setDisabled( !useAge );
  m_sourceInfoAtMeasurement->setDisabled( !useAge );
  m_sourceAgeAtAssay->setDisabled( !useAge );
  
  if( useAge )
  {
    // If the text fields are empty, the red error background wont show up, so we'll at least put
    //  a space there
    if( m_assayDate->text().empty() )
      m_assayDate->setText( " " );
    
    if( m_drfMeasurementDate->text().empty() )
      m_drfMeasurementDate->setText( " " );
    
    validateDateFields();
  }//if( useAge )
  
  m_updated.emit();
}//void useAgeInfoUserToggled()


void MakeDrfSrcDef::updateAgedText()
{
  try
  {
    if( !m_nuclide )
      throw runtime_error( "" );
    
    double activity = activityAtSpectrumTime();
    double age = ageAtSpectrumTime();
    const string agestr = PhysicalUnits::printToBestTimeUnits(age);
    
    const string enteredAct = m_activityEdit->text().toUTF8();
    
    const bool useCurrie = (enteredAct.find_first_of( "cC" ) != string::npos);
    const string actstr = PhysicalUnits::printToBestActivityUnits(activity, 1, useCurrie );
    
    string txt = actstr;
    if( !PeakDef::ageFitNotAllowed(m_nuclide) )
      txt += ", " + agestr;

    m_sourceInfoAtMeasurement->setText( WString::fromUTF8(txt) );
  }catch( std::exception & )
  {
    m_sourceInfoAtMeasurement->setText( "" );
  }
}//void updateAgedText()


void MakeDrfSrcDef::handleUserChangedDistance()
{
  try
  {
    distance();
    if( m_distanceEdit->hasStyleClass( "SrcInputError" ) )
      m_distanceEdit->removeStyleClass( "SrcInputError" );
  }catch( std::exception & )
  {
    if( !m_distanceEdit->hasStyleClass( "SrcInputError" ) )
      m_distanceEdit->addStyleClass( "SrcInputError" );
  }
  
  m_updated.emit();
}//void handleUserChangedDistance()


void MakeDrfSrcDef::handleUserChangedActivity()
{
  try
  {
    enteredActivity();
    if( m_activityEdit->hasStyleClass( "SrcInputError" ) )
      m_activityEdit->removeStyleClass( "SrcInputError" );
  }catch( std::exception & )
  {
    if( !m_activityEdit->hasStyleClass( "SrcInputError" ) )
      m_activityEdit->addStyleClass( "SrcInputError" );
  }
  
  if( m_useAgeInfo->isChecked() )
    updateAgedText();
  
  m_updated.emit();
}//void handleUserChangedActivity()


void MakeDrfSrcDef::handleUserChangedActivityUncertainty()
{
  if( WValidator::State::Valid == m_activityUncertainty->validate() )
  {
    if( m_activityUncertainty->hasStyleClass( "SrcInputError" ) )
      m_activityUncertainty->removeStyleClass( "SrcInputError" );
  }else
  {
    if( !m_activityUncertainty->hasStyleClass( "SrcInputError" ) )
      m_activityUncertainty->addStyleClass( "SrcInputError" );
  }
  
  m_updated.emit();
}//void handleUserChangedActivityUncertainty();


void MakeDrfSrcDef::handleUserChangedAgeAtAssay()
{
  string agestr = m_sourceAgeAtAssay->text().toUTF8();
  SpecUtils::trim( agestr );
  
  double age = 0.0;
  if( agestr.empty() || (agestr.find_first_not_of("+-0.")==string::npos) )
  {
    m_sourceAgeAtAssay->setText( "0 uCi" );
  }else
  {
    try
    {
      const double hl = m_nuclide ? m_nuclide->halfLife : -1.0;
      age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( agestr, hl );
      if( m_activityEdit->hasStyleClass( "SrcInputError" ) )
        m_activityEdit->removeStyleClass( "SrcInputError" );
    }catch( std::exception &e )
    {
      if( !m_activityEdit->hasStyleClass( "SrcInputError" ) )
        m_activityEdit->addStyleClass( "SrcInputError" );
    }
  }//if( zero / else )
  
  updateAgedText();
  
  m_updated.emit();
}//void handleUserChangedAgeAtAssay()


void MakeDrfSrcDef::validateDateFields()
{
  // Only validate if we are actually using the date fields
  if( !m_useAgeInfo->isChecked() )
    return;
  
  string txt = m_drfMeasurementDate->text().toUTF8();
  if( (txt.length() > 1) && std::isspace(txt[0]) )
    m_drfMeasurementDate->setText( WString::fromUTF8( SpecUtils::trim_copy(txt) ) );
  
  txt = m_assayDate->text().toUTF8();
  if( (txt.length() > 1) && std::isspace(txt[0]) )
    m_assayDate->setText( WString::fromUTF8( SpecUtils::trim_copy(txt) ) );
  
  if( m_drfMeasurementDate->validate() == Wt::WValidator::Valid )
  {
    if( m_drfMeasurementDate->hasStyleClass( "SrcInputError" ) )
      m_drfMeasurementDate->removeStyleClass( "SrcInputError" );
  }else
  {
    if( !m_drfMeasurementDate->hasStyleClass( "SrcInputError" ) )
      m_drfMeasurementDate->addStyleClass( "SrcInputError" );
  }
  
  if( m_assayDate->validate() == Wt::WValidator::Valid )
  {
    if( m_assayDate->hasStyleClass( "SrcInputError" ) )
      m_assayDate->removeStyleClass( "SrcInputError" );
  }else
  {
    if( !m_assayDate->hasStyleClass( "SrcInputError" ) )
      m_assayDate->addStyleClass( "SrcInputError" );
  }
}//void validateDateFields();


void MakeDrfSrcDef::handleEnteredDatesUpdated()
{
  validateDateFields();
  
  updateAgedText();
  
  m_updated.emit();
}//void handleEnteredDatesUpdated()


void MakeDrfSrcDef::useShieldingInfoUserToggled()
{
  m_shieldingSelect->setHidden( !m_useShielding->isChecked() );
  m_updated.emit();
}//void useShieldingInfoUserToggled()


void MakeDrfSrcDef::handleUserChangedShielding()
{
  m_updated.emit();
}


double MakeDrfSrcDef::enteredActivity() const
{
  string activitystr = m_activityEdit->text().toUTF8();
  
  SpecUtils::trim( activitystr );
  
  return PhysicalUnits::stringToActivity( activitystr );
}//double enteredActivity()


double MakeDrfSrcDef::distance() const
{
  const string txt = m_distanceEdit->text().toUTF8();
  return PhysicalUnits::stringToDistance(txt);
}

const SandiaDecay::Nuclide *MakeDrfSrcDef::nuclide() const
{
  return m_nuclide;
}


double MakeDrfSrcDef::activityAtSpectrumTime() const
{
  const double userActivity = enteredActivity();
  if( !m_useAgeInfo->isChecked() || !m_nuclide )
    return userActivity;
  
  const WDate measDate = m_drfMeasurementDate->date();
  const WDate assayDate = m_assayDate->date();
  
  if( !measDate.isValid() )
    throw runtime_error( "Measurement date invalid" );
  
  if( !assayDate.isValid() )
    throw runtime_error( "Assay date invalid" );
  
  if( assayDate > measDate )
    throw runtime_error( "Assay date must be before measurement date" );
  
  const int numDays = assayDate.daysTo( measDate );
  
  SandiaDecay::NuclideMixture mix;
  mix.addNuclideByActivity( m_nuclide, userActivity );
  
  return mix.activity( 24.0*3600.0*numDays, m_nuclide );
}//double activityAtSpectrumTime() const


double MakeDrfSrcDef::fractionalActivityUncertainty() const
{
  switch( m_activityUncertainty->validate() )
  {
    case WValidator::State::Invalid:
    {
      //We actually get here if the value is asked for if we have the same
      //  main Wt threadlock as when m_activityUncertainty is created; So we
      //  will, check if the uncertainty string is about what we initially set
      //  it to ("0.0 %"), and if so, just return zero.
      const string value = m_activityUncertainty->valueText().toUTF8();
      if( SpecUtils::istarts_with( value, "0.0 ") )
        return 0.0;
      
      throw runtime_error( "Activity Uncertainty Invalid" );
    }
      
    case WValidator::State::InvalidEmpty:
      return 0.0;
      
    case WValidator::State::Valid:
      break;
  }//switch( m_activityUncertainty->validate() )
  
  return m_activityUncertainty->value() / 100.0;
}//double fractionalActivityUncertainty() const


double MakeDrfSrcDef::ageAtSpectrumTime() const
{
  if( !m_nuclide )
    return 0.0;

  if( !m_useAgeInfo->isChecked() || PeakDef::ageFitNotAllowed(m_nuclide) )
    return PeakDef::defaultDecayTime( m_nuclide, nullptr );
  
  string ageAtAssaystr = m_sourceAgeAtAssay->text().toUTF8();
  SpecUtils::trim( ageAtAssaystr );
  
  double ageAtAssay = 0.0;
  if( !ageAtAssaystr.empty() && (ageAtAssaystr.find_first_not_of("+-0.")!=string::npos) )
    ageAtAssay = PhysicalUnits::stringToTimeDurationPossibleHalfLife( ageAtAssaystr, m_nuclide->halfLife );
  
  const WDate measDate = m_drfMeasurementDate->date();
  const WDate assayDate = m_assayDate->date();
  
  if( !measDate.isValid() )
    throw runtime_error( "Measurement date invalid" );
  
  if( !assayDate.isValid() )
    throw runtime_error( "Assay date invalid" );
  
  if( assayDate > measDate )
    throw runtime_error( "Assay date must be before measurement date" );
  
  const double numDays = assayDate.daysTo( measDate );
  
  return ageAtAssay + 24.0*3600.0*numDays;
}//double ageAtSpectrumTime() const


ShieldingSelect *MakeDrfSrcDef::shielding()
{
  if( !m_useShielding->isChecked() )
    return nullptr;
  return m_shieldingSelect;
}//ShieldingSelect *shielding()


void MakeDrfSrcDef::setDistance( const double dist )
{
  m_distanceEdit->setText( PhysicalUnits::printToBestLengthUnits(dist) );
}//void setDistance( const double dist );


void MakeDrfSrcDef::setActivity( const double act )
{
  const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  const int ndecimals = 4;
  m_activityEdit->setText( PhysicalUnits::printToBestActivityUnits(act, ndecimals, useCi) );
  updateAgedText();
}//void MakeDrfSrcDef::setActivity( const double dist )


void MakeDrfSrcDef::setAssayInfo( const double activity,
                                  const boost::posix_time::ptime &assay_date )
{
  m_useAgeInfo->setChecked( !assay_date.is_special() );
  useAgeInfoUserToggled();
  m_assayDate->setDate( WDateTime::fromPosixTime(assay_date).date() );
  
  if( activity > 0.0 )
  {
    const int ndecimals = 4;
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    m_activityEdit->setText( PhysicalUnits::printToBestActivityUnits(activity, ndecimals, useCi) );
  }
  
  validateDateFields();
  updateAgedText();
}//void setAssayInfo(..);


/*
void MakeDrfSrcDef::setAgeAtAssay( const double age )
{
  if( age >= 0.0 )
  {
    m_useAgeInfo->setChecked( true );
    useAgeInfoUserToggled();
    if( !m_assayDate->date().isValid() )
      m_assayDate->setDate( m_drfMeasurementDate->date() );
    m_sourceAgeAtAssay->setText( PhysicalUnits::printToBestTimeUnits(age) );
  }else
  {
    setNuclide( m_nuclide );
  }
  
  updateAgedText();
}//void setAgeAtAssay( const double age );
*/

void MakeDrfSrcDef::setAgeAtMeas( const double age )
{
  if( age < 0.0 )
  {
    setNuclide( m_nuclide );
    updateAgedText();
  }
  
  m_useAgeInfo->setChecked( true );
  useAgeInfoUserToggled();
    
  WDate measDate = m_drfMeasurementDate->date();
  WDate assayDate = m_assayDate->date();
    
  if( !measDate.isValid() && !assayDate.isValid() )
  {
    m_drfMeasurementDate->setDate( WDate::currentDate() );
    m_assayDate->setDate( WDate::currentDate() );
    measDate = assayDate = WDate::currentDate();
  }else if( !assayDate.isValid() )
  {
    m_assayDate->setDate( measDate );
    assayDate = measDate;
  }else if( !measDate.isValid() )
  {
    m_drfMeasurementDate->setDate( assayDate );
    measDate = assayDate;
  }else if( assayDate > measDate )
  {
    m_assayDate->setDate( measDate );
    assayDate = measDate;
  }
    
  double assayToMeasTime = assayDate.daysTo(measDate) * PhysicalUnits::day;
  
  if( assayToMeasTime > age )
  {
    m_assayDate->setDate( measDate );
    assayDate = measDate;
    assayToMeasTime = 0.0;
  }
  
  m_sourceAgeAtAssay->setText( PhysicalUnits::printToBestTimeUnits(age-assayToMeasTime) );

  validateDateFields();
  updateAgedText();
}//void setAgeAtMeas( const double age );


void MakeDrfSrcDef::setShielding( const float atomic_number, const float areal_density )
{
  m_useShielding->setChecked();
  m_shieldingSelect->setHidden( false );
  m_shieldingSelect->setAtomicNumberAndArealDensity( atomic_number, areal_density );
}//void setShielding( const float atomic_number, const float areal_density )


std::string MakeDrfSrcDef::toGadrasLikeSourceString() const
{
  string answer;
  if( m_nuclide )
  {
    //m_nuclide->symbol == 'U235m2'
    const size_t numpos = m_nuclide->symbol.find_first_of("0123456789");
    assert( numpos != string::npos );
    string numbers = m_nuclide->symbol.substr(numpos);  //235m2
    string meta;
    const size_t metapos = numbers.find_first_not_of("0123456789");
    if( metapos != string::npos )
    {
      meta = numbers.substr(metapos);  //m2
      numbers = numbers.substr(0,metapos); //235
    }
    
    answer = numbers + m_nuclide->symbol.substr(0,numpos) + meta;
  }//if( m_nuclide )
  
  if( !answer.empty() )
    answer += ",";
  
  const double activity = activityAtSpectrumTime();
  const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  answer += PhysicalUnits::printToBestActivityUnits(activity,5,useCi);
  
  if( m_shieldingSelect && m_useShielding->isChecked() )
  {
    double an = 14, ad = 0.0;
    if( m_shieldingSelect->isGenericMaterial() )
    {
      if( !m_shieldingSelect->atomicNumberEdit()->text().empty()
         && !m_shieldingSelect->arealDensityEdit()->text().empty() )
      {
        an = m_shieldingSelect->atomicNumber();
        ad = m_shieldingSelect->arealDensity();
      }
    }else
    {
      std::shared_ptr<Material> mat = m_shieldingSelect->material();
      if( mat )
      {
        an = mat->massWeightedAtomicNumber();
        ad = mat->density * m_shieldingSelect->thickness();
      }//if( mat )
    }//if( shield->isGenericMaterial() ) / else
    
    answer += "{" + std::to_string(an) + "," + std::to_string(ad) + "}";
  }//if( user is using shielding )
  
  const bool mayEvolve = (m_nuclide && !PeakDef::ageFitNotAllowed(m_nuclide));
  const double age = ageAtSpectrumTime();
  if( mayEvolve && age >= 0.0 )
    answer += " Age=" + PhysicalUnits::printToBestTimeUnits(age);
  
  return answer;
}//std::string toGadrasLikeSourceString() const
