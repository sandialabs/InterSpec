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

#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/IsotopeNameFilterModel.h"

using namespace std;
using namespace Wt;

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
  m_nuclideEdit( nullptr ),
  m_distanceEdit( nullptr ),
  m_activityEdit( nullptr ),
  m_activityUnits( nullptr ),
  m_activityUncertainty( nullptr ),
  m_useAgeInfo( nullptr ),
  m_activityDate( nullptr ),
  m_drfMeasurementDate( nullptr ),
  m_sourceActivityAtMeasurement( nullptr ),
  m_sourceAgeAtMeasurement( nullptr ),
//  m_sourceCreationDate( nullptr ),
  m_useShielding( nullptr ),
  m_shieldingSelect( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/MakeDrfSrcDef.css" );
  
  create();
  
  if( m_nuclide )
    m_nuclideEdit->setText( WString::fromUTF8(m_nuclide->symbol) );
  
  if( !measDate.is_special() )
  {
    const string timestr = UtilityFunctions::to_extended_iso_string( measDate );
    m_drfMeasurementDate->setDate( WDate::fromString( WString::fromUTF8(timestr), "%Y-%m-%d%n%H:%M:%S") );  //Not tested yet
  }
}//MakeDrfSrcDef constructor


MakeDrfSrcDef::~MakeDrfSrcDef()
{
}


void MakeDrfSrcDef::create()
{
  m_table = new WTable( this );
  
  WTableCell *cell = m_table->elementAt(0,0);
  WLabel *label = new WLabel( "Nuclide", cell );
  cell = m_table->elementAt(0,1);
  cell->setColumnSpan( 2 );
  m_nuclideEdit = new WLineEdit( cell );
  m_nuclideEdit->setAutoComplete( false );
  label->setBuddy( m_nuclideEdit );
  
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
  nucSuggest->filterModel().connect( filterModel, &IsotopeNameFilterModel::filter );
  
  m_nuclideEdit->changed().connect( this, &MakeDrfSrcDef::handleNuclideChanged );
  
  string tooltip = "ex. <b>U235</b>, <b>235 Uranium</b>, "
                   "<b>U-235m</b> (meta stable state), <b>Cs137</b>, etc.";
  HelpSystem::attachToolTipOn( m_nuclideEdit, tooltip, false );
  
  cell = m_table->elementAt(1,0);
  label = new WLabel( "Distance", cell );
  cell = m_table->elementAt(1,1);
  cell->setColumnSpan( 2 );
  m_distanceEdit = new WLineEdit( cell );
  m_distanceEdit->setTextSize( 5 );
  m_distanceEdit->setAutoComplete( false );
  label->setBuddy( m_distanceEdit );
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_distanceEdit->setValidator( distValidator );
  m_distanceEdit->setText( "50 cm" );
  
  cell = m_table->elementAt(2,1);
  cell->setColumnSpan( 2 );
  m_useAgeInfo = new WCheckBox( "Enter Aging Info", cell );
  m_useAgeInfo->setChecked( false );
  m_useAgeInfo->changed().connect( this, &MakeDrfSrcDef::useAgeInfoUserToggled );
  
  
  cell = m_table->elementAt(3,0);
  label = new WLabel( "Activity", cell );
  
  cell = m_table->elementAt(3,1);
  m_activityEdit = new WLineEdit( cell );
  m_activityEdit->setAutoComplete( false );
  label->setBuddy( m_activityEdit );
  
  WRegExpValidator *val = new WRegExpValidator( PhysicalUnits::sm_activityUnitOptionalRegex, this );
  val->setFlags( Wt::MatchCaseInsensitive );
  m_activityEdit->setValidator( val );
  m_activityEdit->setText( "100" );
  m_activityEdit->changed().connect( this, &MakeDrfSrcDef::handleUserChangedActivity );
  m_activityEdit->enterPressed().connect( this, &MakeDrfSrcDef::handleUserChangedActivity );
  
  cell = m_table->elementAt(3,2);
  m_activityUnits = new WComboBox( cell );
  
  for( const auto &u : PhysicalUnits::sm_activityUnitHtmlNameValues )
  {
    string val = u.first;
    const unsigned char utf8mu[] = { 0xCE, 0xBC, 0 };
    UtilityFunctions::ireplace_all( val, "&mu;", (const char *)utf8mu /*"\u03BC"*/ );
    m_activityUnits->addItem( WString::fromUTF8(val) );
  }
  
  m_activityUnits->setCurrentIndex( 6 ); //uci

  cell = m_table->elementAt(4,0);
  label = new WLabel( "Act. Uncert.", cell );
  cell = m_table->elementAt(4,1);
  cell->setColumnSpan( 2 );
  m_activityUncertainty = new WDoubleSpinBox( cell );
  m_activityUncertainty->setText( "0" );
  m_activityUncertainty->setTextSize( 6 );
  m_activityUncertainty->setAutoComplete( false );
  m_activityUncertainty->setSuffix( " %" );
  label->setBuddy( m_activityUncertainty );
  
  m_activityUncertainty->setDecimals( 1 );
  m_activityUncertainty->setRange( 0.0, 100.0 );
  m_activityUncertainty->setSingleStep( 1.0 );
  
  //WDoubleValidator *percentVal = new WDoubleValidator( this );
  //percentVal->setRange( 0.0, 100.0 );
  //m_activityUncertainty->setValidator( percentVal );

  
  cell = m_table->elementAt(5,0);
  label = new WLabel( "Assay Date", cell );
  cell = m_table->elementAt(5,1);
  cell->setColumnSpan( 2 );
  m_activityDate = new WDateEdit( cell );
  label->setBuddy( m_activityDate );
  
  cell = m_table->elementAt(6,0);
  label = new WLabel( "Spectrum Date", cell );
  cell = m_table->elementAt(6,1);
  cell->setColumnSpan( 2 );
  m_drfMeasurementDate = new WDateEdit( cell );
  label->setBuddy( m_drfMeasurementDate );
  
  
  cell = m_table->elementAt(7,0);
  label = new WLabel( "At Meas. Time", cell );
  cell = m_table->elementAt(7,1);
  m_sourceActivityAtMeasurement = new WText( cell );
  cell = m_table->elementAt(7,2);
  m_sourceAgeAtMeasurement = new WText( cell );
  
  //Wt::WDateEdit *m_sourceCreationDate;
  
  cell = m_table->elementAt(8,1);
  cell->setColumnSpan( 2 );
  m_useShielding = new WCheckBox( "Shielded?", cell );
  m_useShielding->setChecked( false );
  m_useShielding->changed().connect( this, &MakeDrfSrcDef::useShiledingInfoUserToggled );
  
  
  cell = m_table->elementAt(9,0);
  cell->setColumnSpan( 3 );
  m_shieldingSelect = new ShieldingSelect( m_materialDB, nullptr, m_materialSuggest, false, cell );
  m_shieldingSelect->hide();


  useAgeInfoUserToggled();
}//void create()


void MakeDrfSrcDef::handleNuclideChanged()
{
  if( m_useAgeInfo->isChecked() )
  {
    //update m_sourceActivityAtMeasurement
    //update m_sourceAgeAtMeasurement
  }
}//void handleNuclideChanged()


void MakeDrfSrcDef::useAgeInfoUserToggled()
{
  const bool useAge = m_useAgeInfo->isChecked();
  m_activityEdit->label()->setText( useAge ? "Assay Act." : "Activity" );
  
  m_table->rowAt(5)->setHidden( !useAge );
  m_table->rowAt(6)->setHidden( !useAge );
  m_table->rowAt(7)->setHidden( !useAge );
  
  m_activityDate->setDisabled( !useAge );
  m_drfMeasurementDate->setDisabled( !useAge );
  m_sourceActivityAtMeasurement->setDisabled( !useAge );
  m_sourceAgeAtMeasurement->setDisabled( !useAge );
  //m_sourceCreationDate->setDisabled( !useAge );
}//void useAgeInfoUserToggled()


void MakeDrfSrcDef::handleUserChangedActivity()
{
  enteredActivity();
  
  if( m_useAgeInfo->isChecked() )
  {
    //update m_sourceActivityAtMeasurement
  }
}//void handleUserChangedActivity()


void MakeDrfSrcDef::useShiledingInfoUserToggled()
{
  m_shieldingSelect->setHidden( !m_useShielding->isChecked() );
}//void useShiledingInfoUserToggled()


double MakeDrfSrcDef::enteredActivity()
{
  double activity = 0.0;
  string activitystr = m_activityEdit->text().toUTF8();
  
  UtilityFunctions::trim( activitystr );
  
  const size_t pos = activitystr.find_first_not_of( "0123456789Ee.-+" );  //a regex would be better
  
  if( pos == string::npos )
  {
    if( !(stringstream(activitystr) >> activity) )
      throw runtime_error( "Invalid number" );
    const int unitsInd = m_activityUnits->currentIndex();
    activity *= PhysicalUnits::sm_activityUnitHtmlNameValues.at(unitsInd).second;
  }else
  {
    const bool hasb = UtilityFunctions::icontains( activitystr, "b" );
    const bool hasc = UtilityFunctions::icontains( activitystr, "c" );
    
    if( hasb && hasc )
      throw runtime_error( "Invalid activity string, couldnt determine units" );
    
    activity = PhysicalUnits::stringToActivity( activitystr );
    
    const PhysicalUnits::UnitNameValuePair &bestunit
    = PhysicalUnits::bestActivityUnitHtml( activity, hasc );
    
    PhysicalUnits::UnitNameValuePairV::const_iterator pos =
    std::find( PhysicalUnits::sm_activityUnitHtmlNameValues.begin(), PhysicalUnits::sm_activityUnitHtmlNameValues.end(), bestunit );
    assert( pos != PhysicalUnits::sm_activityUnitHtmlNameValues.end() );
    m_activityUnits->setCurrentIndex( pos - PhysicalUnits::sm_activityUnitHtmlNameValues.begin() );
    
    const double txtdblval = activity / bestunit.second;
    char txtval[32];
    snprintf( txtval, sizeof(txtval), "%f", floor(1000*txtdblval + 0.5)/1000 );  //fix this to 3 sig figs, not 3 decimals
    
    //Get rid of unwanted decimals (they should be all zeros from rounding above)
    //(this should probably be a function call somewhere)
    string durstr = txtval;
    const size_t decpos = durstr.find_last_of( '.' );
    if( decpos != string::npos )
    {
      const string predec = durstr.substr( 0, decpos );
      string postdec = durstr.substr( decpos + 1 );
      size_t lastpos = postdec.find_last_not_of( "0" );
      if( lastpos == string::npos )
        postdec = "";
      else
        postdec = postdec.substr( 0, lastpos + 1 );
      
      durstr = predec;
      if( postdec.size() )
        durstr += "." + postdec;
    }//if( decpos != string::npos )
    
    m_activityEdit->setText( durstr );
  }//if( user entered only numbers ) / else user entered units to
  
  return activity;
}//double enteredActivity()


