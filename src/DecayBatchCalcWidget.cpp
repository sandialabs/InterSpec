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

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <Wt/Utils.h>
#include <Wt/WLink.h>
#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WAnchor.h>
#include <Wt/WSpinBox.h>
#include <Wt/WComboBox.h>
#include <Wt/WCheckBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WResource.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/Http/Response.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WDoubleValidator.h>
#include <Wt/WSuggestionPopup.h>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/GroupBox.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DecayBatchCalc.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DecayBatchCalcWidget.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/FileDragUploadResource.h"

#if( USE_QR_CODES )
#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;


namespace
{
  // Safety cap: number of table cells beyond which we skip building the HTML preview (export still
  //  works).  Sized to comfortably clear the motivating grouped file - 35 locations x 103 nuclides x
  //  6 columns = 21630 cells - while still catching a runaway.
  const size_t sm_max_preview_cells = 40000;

  // See CopyFluxDataTextToClipboard in FluxTool.cpp / CopyUrlToClipboard in QrCode.cpp.
  WT_DECLARE_WT_MEMBER
  (CopyDecayBatchToClipboard, Wt::JavaScriptFunction, "CopyDecayBatchToClipboard",
   function( sender, event, id )
  {
    var el = document.getElementById(id);
    var text = el && el._isData ? el._isData.TableData : null;
    if( !text )
      return false;

    var didcopy = 0;
    try
    {
      function listener(e){
        e.clipboardData.setData("text/plain", text);
        didcopy = 1;
        e.preventDefault();
      }
      document.addEventListener("copy", listener);
      document.execCommand("copy");
      document.removeEventListener("copy", listener);
      if( didcopy )
        return didcopy;
    }catch(error){
      console.warn( 'Failed to copy text to clipboard' );
    }

    if( window.clipboardData && window.clipboardData.setData ){
      return window.clipboardData.setData("Text", text);
    }else if( document.queryCommandSupported && document.queryCommandSupported("copy") ){
      var temparea = document.createElement("textarea");
      temparea.textContent = text;
      temparea.style.position = "fixed";
      document.body.appendChild(temparea);
      temparea.select();
      try{
        var copysuccess = document.execCommand("copy");
        return copysuccess ? 2 : 0;
      }catch( ex ){
        return 0;
      }finally{
        document.body.removeChild( temparea );
      }
    }
    return 0;
  }
  );


  /** WResource serving the current result table as a downloadable CSV. */
  class DecayBatchCalcCsvResource : public Wt::WResource
  {
    DecayBatchCalcWidget *m_widget;
    Wt::WApplication *m_app;

  public:
    DecayBatchCalcCsvResource( DecayBatchCalcWidget *widget )
    : WResource(), m_widget( widget ), m_app( WApplication::instance() )
    {
      assert( m_app );
      // Hold the session lock across handleRequest so we don't race ~WResource; see FluxCsvResource.
      setTakesUpdateLock( true );
    }

    virtual ~DecayBatchCalcCsvResource()
    {
      beingDeleted();
    }

    virtual void handleRequest( const Wt::Http::Request &, Wt::Http::Response &response )
    {
      suggestFileName( "batch_decay.csv", ContentDisposition::Attachment );
      response.setMimeType( "text/csv" );
      if( m_widget )
        response.out() << m_widget->currentResultCsv();
    }
  };//class DecayBatchCalcCsvResource


  /** Formats an activity for a row's edit box, in the given unit token (e.g. "uCi"), keeping enough
   digits that reading the text back gives the same value.

   `PhysicalUnits::printToBestActivityUnits` is not usable here: its digit count is *decimal places*,
   not significant figures, so a measured 9.53E-15 uCi comes back as "0.0000 fCi" - and since the row
   text is what `gatherInputs` re-parses, that silently turns a real measurement into zero (152 of the
   1400 rows of the motivating file).  Values here span ~30 decades, so use %.<n>G, which keeps
   significant figures wherever the exponent lands.
   */
  std::string activity_text_in_unit( const double activity, const std::string &unit_token )
  {
    // 9 significant figures: more than any input file supplies (the motivating one gives 7), and well
    //  inside a double's ~15-17, so the text round-trips to the same value the file held.
    double scale = 1.0;
    try
    {
      if( !unit_token.empty() )
        scale = PhysicalUnits::stringToActivity( "1" + unit_token );
    }catch( std::exception & )
    {
      scale = 0.0;
    }

    if( scale <= 0.0 )
      scale = PhysicalUnits::becquerel;   // unparseable unit; fall back to the SandiaDecay unit

    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.9G", activity / scale );

    std::string answer = buffer;
    if( !unit_token.empty() )
      answer += " " + unit_token;

    return answer;
  }//activity_text_in_unit(...)


  /** Parses an activity string that may carry units ("5 uCi") or be a bare scalar ("5", interpreted
   as becquerel).  Returns false if it does not resolve.  On success `value` is in SandiaDecay units
   (becquerel), `has_unit` says whether an explicit unit was given, and `unit_token` is that unit
   (e.g. "uCi") or empty for a bare scalar.
   Allowing either form lets the user work in relative amounts (all scalars) or absolute activities
   (all unit-bearing); DecayBatchCalcWidget keeps the whole set to a single style. */
  bool parse_activity_text( const std::string &input, double &value, bool &has_unit,
                            std::string &unit_token )
  {
    value = 0.0;
    has_unit = false;
    unit_token.clear();

    std::string txt = input;
    SpecUtils::trim( txt );
    if( txt.empty() )
      return false;

    // Prefer an activity value with explicit units (e.g. "5 uCi", "3.2mCi").
    try
    {
      value = PhysicalUnits::stringToActivity( txt );
      has_unit = true;

      // The unit is whatever follows the leading numeric value.  Use std::stod to find where the
      //  number (including any scientific-notation exponent) ends, so a leading unit letter (e.g. the
      //  'E' of a hypothetical exa-prefix) is never mistaken for part of the number.
      size_t num_end = 0;
      try{ std::stod( txt, &num_end ); }catch( std::exception & ){ num_end = 0; }
      unit_token = txt.substr( num_end );
      SpecUtils::trim( unit_token );

      return true;
    }catch( std::exception & )
    {
    }

    // Otherwise accept a bare scalar, interpreted as becquerel.
    try
    {
      size_t end_pos = 0;
      const double val = std::stod( txt, &end_pos );
      if( txt.find_first_not_of( " \t", end_pos ) != std::string::npos )
        return false;
      value = val * PhysicalUnits::becquerel;
      return true;
    }catch( std::exception & )
    {
    }

    return false;
  }//parse_activity_text(...)
}//anonymous namespace


DecayBatchCalcNuclide::DecayBatchCalcNuclide( WSuggestionPopup *nuclideSuggest )
  : WContainerWidget(),
    m_nuclideSuggest( nuclideSuggest ),
    m_nuclideEdit( nullptr ),
    m_ageContainer( nullptr ),
    m_ageEdit( nullptr ),
    m_activityEdit( nullptr ),
    m_removeBtn( nullptr ),
    m_unitLabel(),
    m_prevAgeText(),
    m_prevActivityText( "1" ),
    m_changed(),
    m_remove()
{
  addStyleClass( "DecayBatchCalcNuclide" );

  m_nuclideEdit = addNew<WLineEdit>();
  m_nuclideEdit->setAutoComplete( false );
  m_nuclideEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_nuclideEdit->setAttributeValue( "autocorrect", "off" );
  m_nuclideEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_nuclideEdit->setWidth( WLength(80.0, WLength::Unit::Pixel) );
  m_nuclideEdit->setPlaceholderText( WString::tr("dbc-nuc-placeholder") );
  if( nuclideSuggest )
    nuclideSuggest->forEdit( m_nuclideEdit, PopupTrigger::Editing );
  m_nuclideEdit->changed().connect( this, &DecayBatchCalcNuclide::handleNuclideChange );
  m_nuclideEdit->enterPressed().connect( this, &DecayBatchCalcNuclide::handleNuclideChange );

  // Age input (label + edit); only meaningful for ageable nuclides, so the whole group is hidden for
  //  non-ageable ones (matching the "Isotopics by nuclides" row behavior).
  m_ageContainer = addNew<WContainerWidget>();
  m_ageContainer->addStyleClass( "DbcAgeContainer" );
  WLabel *ageLabel = m_ageContainer->addNew<WLabel>( WString::tr("dbc-age-label") );
  m_ageEdit = m_ageContainer->addNew<WLineEdit>();
  ageLabel->setBuddy( m_ageEdit );
  m_ageEdit->setWidth( WLength(70.0, WLength::Unit::Pixel) );
  auto ageValidator = make_shared<WRegExpValidator>( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex() );
  ageValidator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  m_ageEdit->setValidator( ageValidator );
  m_ageEdit->setAutoComplete( false );
  m_ageEdit->setAttributeValue( "ondragstart", "return false" );
  m_ageEdit->changed().connect( this, &DecayBatchCalcNuclide::handleAgeChange );
  m_ageEdit->enterPressed().connect( this, &DecayBatchCalcNuclide::handleAgeChange );

  // Activity input (label + edit); validated manually so a bare scalar (no units) is accepted.
  WContainerWidget *actContainer = addNew<WContainerWidget>();
  actContainer->addStyleClass( "DbcActContainer" );
  WLabel *actLabel = actContainer->addNew<WLabel>( WString::tr("dbc-act-label") );
  m_activityEdit = actContainer->addNew<WLineEdit>( "1" );
  actLabel->setBuddy( m_activityEdit );
  m_activityEdit->setWidth( WLength(80.0, WLength::Unit::Pixel) );
  m_activityEdit->setAutoComplete( false );
  m_activityEdit->setAttributeValue( "ondragstart", "return false" );
  m_activityEdit->changed().connect( this, &DecayBatchCalcNuclide::handleActivityChange );
  m_activityEdit->enterPressed().connect( this, &DecayBatchCalcNuclide::handleActivityChange );

  m_removeBtn = addNew<WPushButton>();
  m_removeBtn->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
  m_removeBtn->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  m_removeBtn->clicked().connect( [this](){ m_remove.emit(this); } );

  updateAgeEnabledState();
}//DecayBatchCalcNuclide constructor


DecayBatchCalcNuclide::~DecayBatchCalcNuclide()
{
  // The suggestion popup is shared across all rows and keeps a raw pointer to each registered edit;
  //  it does not drop the edit when the edit is destroyed.  Detach ours here (while it is still
  //  alive, before the base WContainerWidget destructor deletes our children) so a later suggestion
  //  activation does not dereference this deleted row's edit.  m_nuclideSuggest is an observing_ptr,
  //  so this is skipped if the popup was already torn down.
  if( m_nuclideSuggest && m_nuclideEdit )
    m_nuclideSuggest->removeEdit( m_nuclideEdit );
}//DecayBatchCalcNuclide destructor


const SandiaDecay::Nuclide *DecayBatchCalcNuclide::nuclide() const
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    return nullptr;
  return db->nuclide( m_nuclideEdit->text().toUTF8() );
}//nuclide()


double DecayBatchCalcNuclide::age() const
{
  const SandiaDecay::Nuclide * const nuc = nuclide();
  if( !nuc || !m_ageEdit->isEnabled() )
    return 0.0;

  const string age_str = m_ageEdit->text().toUTF8();
  if( age_str.empty() )
    return 0.0;

  try
  {
    const double age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( age_str, nuc->halfLife );
    return (age >= 0.0) ? age : 0.0;
  }catch( std::exception & )
  {
    return 0.0;
  }
}//age()


double DecayBatchCalcNuclide::activity() const
{
  double value = 0.0;
  bool has_unit = false;
  string unit_token;
  if( parse_activity_text( m_activityEdit->text().toUTF8(), value, has_unit, unit_token ) )
    return value;
  return 0.0;
}//activity()


bool DecayBatchCalcNuclide::isValid() const
{
  const SandiaDecay::Nuclide * const nuc = nuclide();
  if( !nuc || nuc->isStable() )
    return false;

  // A zero activity is a legitimate input - a measurement of "none detected" - so only unparseable or
  //  negative text is invalid.  (Grouped input relies on this: dropping a location's zero rows would
  //  leave its mixture missing nuclides the other locations have, and a location measured as all
  //  zeroes would disappear from the output entirely.)
  string txt = m_activityEdit->text().toUTF8();
  SpecUtils::trim( txt );
  if( txt.empty() )
    return false;

  double value = 0.0;
  bool has_unit = false;
  string unit_token;
  if( !parse_activity_text( txt, value, has_unit, unit_token ) )
    return false;

  return (value >= 0.0);
}//isValid()


bool DecayBatchCalcNuclide::activityHasUnit() const
{
  double value = 0.0;
  bool has_unit = false;
  string unit_token;
  if( parse_activity_text( m_activityEdit->text().toUTF8(), value, has_unit, unit_token ) )
    return has_unit;
  return false;
}//activityHasUnit()


string DecayBatchCalcNuclide::activityUnitStr() const
{
  double value = 0.0;
  bool has_unit = false;
  string unit_token;
  if( parse_activity_text( m_activityEdit->text().toUTF8(), value, has_unit, unit_token ) )
    return unit_token;
  return string();
}//activityUnitStr()


void DecayBatchCalcNuclide::setNuclide( const SandiaDecay::Nuclide *nuc, double activity,
                                        bool useCurie, double age, const string &activityStr )
{
  m_unitLabel.clear();  // fresh population; any prior areal/volumetric label no longer applies
  m_nuclideEdit->setText( nuc ? WString::fromUTF8(nuc->symbol) : WString() );

  const string actStr = activityStr.empty()
    ? PhysicalUnits::printToBestActivityUnits( activity, 4, useCurie )
    : activityStr;
  m_activityEdit->setText( WString::fromUTF8(actStr) );
  m_prevActivityText = actStr;

  updateAgeEnabledState();
  if( m_ageEdit->isEnabled() && (age > 0.0) )
    m_ageEdit->setText( PhysicalUnits::printToBestTimeUnits( age, 4 ) );
  else
    m_ageEdit->setText( "" );
  m_prevAgeText = m_ageEdit->text().toUTF8();
}//setNuclide(...)


string DecayBatchCalcNuclide::nuclideText() const
{
  return m_nuclideEdit->text().toUTF8();
}


string DecayBatchCalcNuclide::activityText() const
{
  return m_activityEdit->text().toUTF8();
}


string DecayBatchCalcNuclide::ageText() const
{
  return m_ageEdit->text().toUTF8();
}


void DecayBatchCalcNuclide::setActivityText( const string &txt )
{
  m_activityEdit->setText( WString::fromUTF8(txt) );
  m_prevActivityText = txt;
}//setActivityText(...)


void DecayBatchCalcNuclide::setUnitLabel( const string &label )
{
  m_unitLabel = label;
}


string DecayBatchCalcNuclide::unitLabel() const
{
  return m_unitLabel;
}


void DecayBatchCalcNuclide::focusNuclideEdit()
{
  m_nuclideEdit->setFocus( true );
}


Wt::Signal<> &DecayBatchCalcNuclide::changed()
{
  return m_changed;
}


Wt::Signal<DecayBatchCalcNuclide *> &DecayBatchCalcNuclide::removeRequested()
{
  return m_remove;
}


void DecayBatchCalcNuclide::handleNuclideChange()
{
  updateAgeEnabledState();
  m_changed.emit();
}


void DecayBatchCalcNuclide::handleAgeChange()
{
  const string txt = m_ageEdit->text().toUTF8();

  bool ok = txt.empty();  // an empty age means "as-decayed / zero age"
  if( !ok )
  {
    try
    {
      const SandiaDecay::Nuclide * const nuc = nuclide();
      PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( txt, nuc ? nuc->halfLife : 0.0 );
      ok = true;
    }catch( std::exception & )
    {
    }
  }//if( !txt.empty() )

  if( ok )
  {
    m_prevAgeText = txt;
    m_changed.emit();
  }else
  {
    m_ageEdit->setText( WString::fromUTF8(m_prevAgeText) );  // revert: the rejected entry disappears
  }
}//handleAgeChange()


void DecayBatchCalcNuclide::handleActivityChange()
{
  double value = 0.0;
  bool has_unit = false;
  string unit_token;
  const string txt = m_activityEdit->text().toUTF8();

  // Zero is accepted: it is a valid "none detected" measurement (see isValid()).
  if( parse_activity_text( txt, value, has_unit, unit_token ) && (value >= 0.0) )
  {
    m_prevActivityText = txt;
    m_changed.emit();
  }else
  {
    m_activityEdit->setText( WString::fromUTF8(m_prevActivityText) );  // revert
  }
}//handleActivityChange()


void DecayBatchCalcNuclide::updateAgeEnabledState()
{
  const SandiaDecay::Nuclide * const nuc = nuclide();
  const bool ageable = (nuc && !PeakDef::ageFitNotAllowed(nuc));
  // Keep the age label + input present (so rows stay aligned); just disable them when not applicable.
  m_ageContainer->setDisabled( !ageable );
  m_ageEdit->setEnabled( ageable );
  if( !ageable )
  {
    m_ageEdit->setText( "" );
    m_prevAgeText.clear();
  }
}//updateAgeEnabledState()




DecayBatchCalcLocation::DecayBatchCalcLocation( const string &location )
  : GroupBox(),
    m_location(),
    m_productSuffix(),
    m_extraColumns(),
    m_activityUnit(),
    m_rows( nullptr )
{
  addStyleClass( "DbcLocation" );

  // GroupBox has no way to un-set its legend, so an untitled group is styled away instead (see
  //  #DecayBatchCalcLocation).  setTitle() creates the legend as child 0; keep the rows below it in
  //  their own container so row indices are independent of it.
  setTitle( WString() );
  m_rows = addNew<WContainerWidget>();
  m_rows->addStyleClass( "DbcLocationRows" );

  setLocation( location );
}//DecayBatchCalcLocation constructor


const string &DecayBatchCalcLocation::location() const
{
  return m_location;
}


void DecayBatchCalcLocation::setLocation( const string &location )
{
  // Trimmed so a whitespace-only name (a stray space in a spreadsheet export) counts as untitled
  //  rather than drawing a border with an empty legend and locking "Mix inputs" on.  This also makes
  //  the CSV path agree with handleAppUrl(), which trims its values before they get here.
  m_location = location;
  SpecUtils::trim( m_location );
  setTitle( WString::fromUTF8(m_location) );
  toggleStyleClass( "DbcUntitled", m_location.empty() );
}//setLocation(...)


void DecayBatchCalcLocation::setPassthrough( const string &product_suffix,
                                             const vector<pair<string,string>> &extra_columns,
                                             const string &activity_unit )
{
  m_productSuffix = product_suffix;
  m_extraColumns = extra_columns;
  m_activityUnit = activity_unit;
}//setPassthrough(...)


const string &DecayBatchCalcLocation::productSuffix() const
{
  return m_productSuffix;
}


const vector<pair<string,string>> &DecayBatchCalcLocation::extraColumns() const
{
  return m_extraColumns;
}


const string &DecayBatchCalcLocation::activityUnit() const
{
  return m_activityUnit;
}


int DecayBatchCalcLocation::numRows() const
{
  return m_rows->count();
}


DecayBatchCalcNuclide *DecayBatchCalcLocation::row( const int index ) const
{
  if( (index < 0) || (index >= m_rows->count()) )
    return nullptr;
  return dynamic_cast<DecayBatchCalcNuclide *>( m_rows->widget(index) );
}//row(...)


DecayBatchCalcNuclide *DecayBatchCalcLocation::addRow( WSuggestionPopup *nuclideSuggest )
{
  return m_rows->addNew<DecayBatchCalcNuclide>( nuclideSuggest );
}


int DecayBatchCalcLocation::removeRow( DecayBatchCalcNuclide *row )
{
  // Rows ask to be removed from their own remove-button signal, so the destruction has to wait for
  //  that emit to unwind (see WidgetUtils::removeWidgetLater).
  WidgetUtils::removeWidgetLater( m_rows, row );
  return m_rows->count();
}//removeRow(...)




DecayBatchCalcWidget::DecayBatchCalcWidget( InterSpec *viewer )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_nuclideRows( nullptr ),
    m_addRowBtn( nullptr ),
    m_nuclideSuggest( nullptr ),
    m_timeEdit( nullptr ),
    m_stepsSpin( nullptr ),
    m_minActEdit( nullptr ),
    m_actUnitsCombo( nullptr ),
    m_mixInput( nullptr ),
    m_showProgeny( nullptr ),
    m_incActivity( nullptr ),
    m_incXrays( nullptr ),
    m_incGammas( nullptr ),
    m_incAlphas( nullptr ),
    m_incBetas( nullptr ),
    m_dropArea( nullptr ),
    m_uploadResource( nullptr ),
    m_resultContainer( nullptr ),
    m_csvDownload( nullptr ),
    m_copyBtn( nullptr ),
    m_infoCopied( this, "dbcCopied", true ),
    m_lastResult(),
    m_suppressUpdate( false ),
    m_previousStateUri(),
    m_trackUndo( false )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
  {
    app->useMessageResourceBundle( "DecayBatchCalc" );
    app->useStyleSheet( "InterSpec_resources/DecayBatchCalc.css" );
    // Provides BatchInputDropUploadSetup(...) / setupOnDragEnterDom(...) used by the CSV drop area.
    app->require( "InterSpec_resources/BatchGuiWidget.js" );
    // Register Wt.WT.CopyDecayBatchToClipboard (declared at the top of this file) for the Copy button.
    LOAD_JAVASCRIPT( app, "DecayBatchCalcWidget.cpp", "DecayBatchCalcWidget", wtjsCopyDecayBatchToClipboard );
  }

  addStyleClass( "DecayBatchCalcWidget" );

  const bool showToolTips = m_interspec
    ? UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec ) : false;

  // --- Nuclide input area -------------------------------------------------------------------
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  auto isoSuggestModel = make_shared<IsotopeNameFilterModel>();
  isoSuggestModel->excludeXrays( true );
  isoSuggestModel->excludeEscapes( true );
  isoSuggestModel->excludeReactions( true );

  m_nuclideSuggest = addNew<WSuggestionPopup>( matcherJs, replacerJs );
  m_nuclideSuggest->addStyleClass( "nuclide-suggest" );
  m_nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::Unit::FontEm) );
  IsotopeNameFilterModel::setQuickTypeFixHackjs( m_nuclideSuggest );
  isoSuggestModel->filter( "" );
  m_nuclideSuggest->setFilterLength( -1 );
  m_nuclideSuggest->setModel( isoSuggestModel );
  m_nuclideSuggest->filterModel().connect( isoSuggestModel.get(), &IsotopeNameFilterModel::filter );

  WContainerWidget *nucArea = addNew<WContainerWidget>();
  nucArea->addStyleClass( "DbcNuclideArea" );

  m_nuclideRows = nucArea->addNew<WContainerWidget>();
  m_nuclideRows->addStyleClass( "DbcNuclideRows" );

  // Footer bar with the "+" button and (logically grouped with it) the CSV drop/upload area.  The
  //  bar's top border draws the horizontal line above it, like the "Isotopics by nuclides" tool.
  WContainerWidget *addBar = nucArea->addNew<WContainerWidget>();
  addBar->addStyleClass( "DbcAddBar" );
  m_addRowBtn = addBar->addNew<WPushButton>();
  m_addRowBtn->setStyleClass( "AddEnergyRangeOrNuc Wt-icon" );
  m_addRowBtn->setIcon( "InterSpec_resources/images/plus_min_black.svg" );
  m_addRowBtn->clicked().connect( [this](){ addEmptyNuclideRow()->focusNuclideEdit(); } );
  HelpSystem::attachToolTipOn( m_addRowBtn, WString::tr("dbc-add-nuc-tt"), showToolTips );

  // Disable the main app's spectrum-file drag/drop while this tool is open.
  doJavaScript( "window._IS=window._IS||{};window._IS.BlockFileDrops=true;" );

  m_uploadResource = make_unique<FileDragUploadResource>();
  m_uploadResource->fileDrop().connect( this, &DecayBatchCalcWidget::handleFileDrop );

  m_dropArea = addBar->addNew<WContainerWidget>();
  m_dropArea->addStyleClass( "DbcDropArea" );
  m_dropArea->addNew<WText>( WString::tr("dbc-drop-hint") );
  m_dropArea->doJavaScript( "BatchInputDropUploadSetup(" + m_dropArea->jsRef()
                            + ", '" + m_uploadResource->url() + "');" );
  doJavaScript( "setupOnDragEnterDom(['" + m_dropArea->id() + "']);" );

  // --- Options area -------------------------------------------------------------------------
  WContainerWidget *optsArea = addNew<WContainerWidget>();
  optsArea->addStyleClass( "DbcOptionsArea" );

  WContainerWidget *timeRow = optsArea->addNew<WContainerWidget>();
  timeRow->addStyleClass( "DbcOptRow" );
  WLabel *label = timeRow->addNew<WLabel>( WString::tr("dbc-decay-time") );
  m_timeEdit = timeRow->addNew<WLineEdit>( "1 y" );
  m_timeEdit->setWidth( WLength(45.0, WLength::Unit::Pixel) );
  label->setBuddy( m_timeEdit );
  // A negative time is allowed: it looks *backwards*, solving for what the activities were then.
  auto timeValidator = make_shared<WRegExpValidator>( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalPosOrNegRegex() );
  timeValidator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  m_timeEdit->setValidator( timeValidator );
  m_timeEdit->setAutoComplete( false );
  m_timeEdit->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  m_timeEdit->enterPressed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  HelpSystem::attachToolTipOn( m_timeEdit, WString::tr("dbc-decay-time-tt"), showToolTips );

  label = timeRow->addNew<WLabel>( WString::tr("dbc-num-steps") );
  m_stepsSpin = timeRow->addNew<WSpinBox>();
  m_stepsSpin->setRange( 1, 1000 );
  m_stepsSpin->setValue( 1 );
  m_stepsSpin->setWidth( WLength(30.0, WLength::Unit::Pixel) );
  label->setBuddy( m_stepsSpin );
  m_stepsSpin->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  m_stepsSpin->valueChanged().connect( [this](int){ scheduleResultUpdate(); } );
  HelpSystem::attachToolTipOn( m_stepsSpin, WString::tr("dbc-num-steps-tt"), showToolTips );

  label = timeRow->addNew<WLabel>( WString::tr("dbc-activity-units") );
  m_actUnitsCombo = timeRow->addNew<WComboBox>();
  m_actUnitsCombo->addItem( WString::tr("dbc-curie") );
  m_actUnitsCombo->addItem( WString::tr("dbc-becquerel") );
  const bool useBq = m_interspec
    ? UserPreferences::preferenceValue<bool>( "DisplayBecquerel", m_interspec ) : false;
  m_actUnitsCombo->setCurrentIndex( useBq ? 1 : 0 );
  label->setBuddy( m_actUnitsCombo );
  m_actUnitsCombo->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );

  WContainerWidget *cutRow = optsArea->addNew<WContainerWidget>();
  cutRow->addStyleClass( "DbcOptRow" );
  label = cutRow->addNew<WLabel>( WString::tr("dbc-min-activity") );
  m_minActEdit = cutRow->addNew<WLineEdit>();
  m_minActEdit->setWidth( WLength(70.0, WLength::Unit::Pixel) );
  m_minActEdit->setPlaceholderText( WString::tr("dbc-min-activity-placeholder") );
  label->setBuddy( m_minActEdit );
  // Unit-less on purpose: it is compared against the number shown in the result table (see
  //  `BatchDecayOptions::min_activity`), whatever unit that happens to be in.
  auto cutValidator = make_shared<WDoubleValidator>();
  cutValidator->setBottom( 0.0 );
  m_minActEdit->setValidator( cutValidator );
  m_minActEdit->setAutoComplete( false );
  m_minActEdit->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  m_minActEdit->enterPressed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  HelpSystem::attachToolTipOn( m_minActEdit, WString::tr("dbc-min-activity-tt"), showToolTips );

  WContainerWidget *cbRow = optsArea->addNew<WContainerWidget>();
  cbRow->addStyleClass( "DbcOptRow" );
  m_mixInput = cbRow->addNew<WCheckBox>( WString::tr("dbc-mix-input") );
  m_mixInput->addStyleClass( "CbNoLineBreak" );
  m_mixInput->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  HelpSystem::attachToolTipOn( m_mixInput, WString::tr("dbc-mix-input-tt"), showToolTips );

  m_showProgeny = cbRow->addNew<WCheckBox>( WString::tr("dbc-show-progeny") );
  m_showProgeny->setChecked( true );
  m_showProgeny->addStyleClass( "CbNoLineBreak" );
  m_showProgeny->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  HelpSystem::attachToolTipOn( m_showProgeny, WString::tr("dbc-show-progeny-tt"), showToolTips );

  WContainerWidget *incRow = optsArea->addNew<WContainerWidget>();
  incRow->addStyleClass( "DbcOptRow" );
  incRow->addNew<WLabel>( WString::tr("dbc-include") );
  m_incActivity = incRow->addNew<WCheckBox>( WString::tr("dbc-inc-activity") );
  m_incActivity->addStyleClass( "CbNoLineBreak" );
  m_incActivity->setChecked( true );
  m_incXrays = incRow->addNew<WCheckBox>( WString::tr("dbc-inc-xrays") );
  m_incXrays->addStyleClass( "CbNoLineBreak" );
  m_incGammas = incRow->addNew<WCheckBox>( WString::tr("dbc-inc-gammas") );
  m_incGammas->addStyleClass( "CbNoLineBreak" );
  m_incAlphas = incRow->addNew<WCheckBox>( WString::tr("dbc-inc-alphas") );
  m_incAlphas->addStyleClass( "CbNoLineBreak" );
  m_incBetas = incRow->addNew<WCheckBox>( WString::tr("dbc-inc-betas") );
  m_incBetas->addStyleClass( "CbNoLineBreak" );
  for( WCheckBox *cb : { m_incActivity, m_incXrays, m_incGammas, m_incAlphas, m_incBetas } )
    cb->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );

  // --- Results area -------------------------------------------------------------------------
  m_resultContainer = addNew<WContainerWidget>();
  m_resultContainer->addStyleClass( "DbcResultArea" );
  m_resultContainer->setOverflow( Wt::Overflow::Auto );

  // --- Export row (children of this widget, NOT the window footer) --------------------------
  WContainerWidget *exportRow = addNew<WContainerWidget>();
  exportRow->addStyleClass( "DbcExportRow" );

  auto csvResource = make_shared<DecayBatchCalcCsvResource>( this );
  WLink csvLink( csvResource );
  csvLink.setTarget( LinkTarget::NewWindow );
  m_csvDownload = exportRow->addNew<WAnchor>( csvLink, WString::tr("dbc-export-csv") );
  m_csvDownload->setStyleClass( "LinkBtn DownloadLink" );
  m_csvDownload->setText( WString::tr("dbc-export-csv") );

  m_copyBtn = exportRow->addNew<WPushButton>( WString::tr("dbc-copy-clipboard") );
  m_copyBtn->setStyleClass( "LinkBtn" );
  m_copyBtn->clicked().connect( "function(s,e){"
    "var success = Wt.WT.CopyDecayBatchToClipboard(s,e,'" + m_copyBtn->id() + "');"
    "Wt.emit( '" + id() + "', {name:'dbcCopied', eventObject:e}, success );"
    "}" );
  m_infoCopied.connect( this, [this]( const int success ){ tableCopiedToClipboardCallback( success ); } );

  // Start with a single, empty input row, then set the initial result/export state.
  addEmptyNuclideRow();
  scheduleResultUpdate();

  // Establish the undo baseline as the just-built state, and begin tracking changes.  Nuclides seeded
  //  by the caller right after construction fold into the "open tool" undo step (same event loop).
  m_previousStateUri = encodeStateToUrl();
  m_trackUndo = true;
}//DecayBatchCalcWidget constructor


DecayBatchCalcWidget::~DecayBatchCalcWidget()
{
  // Re-enable the main app's spectrum-file drag/drop.
  WApplication *app = wApp;
  if( app )
  {
    app->doJavaScript( "if(window._IS)window._IS.BlockFileDrops=null;" );
    if( m_dropArea )
      app->doJavaScript( "removeOnDragEnterDom(['" + m_dropArea->id() + "']);" );
  }
}//~DecayBatchCalcWidget()


vector<DecayBatchCalcLocation *> DecayBatchCalcWidget::locationGroups() const
{
  vector<DecayBatchCalcLocation *> groups;
  for( int i = 0; i < m_nuclideRows->count(); ++i )
  {
    DecayBatchCalcLocation *group = dynamic_cast<DecayBatchCalcLocation *>( m_nuclideRows->widget(i) );
    if( group )
      groups.push_back( group );
  }
  return groups;
}//locationGroups()


vector<DecayBatchCalcNuclide *> DecayBatchCalcWidget::allRows() const
{
  vector<DecayBatchCalcNuclide *> rows;
  for( DecayBatchCalcLocation * const group : locationGroups() )
  {
    for( int i = 0; i < group->numRows(); ++i )
    {
      DecayBatchCalcNuclide *row = group->row( i );
      if( row )
        rows.push_back( row );
    }
  }
  return rows;
}//allRows()


DecayBatchCalcLocation *DecayBatchCalcWidget::groupOf( DecayBatchCalcNuclide *row ) const
{
  for( DecayBatchCalcLocation * const group : locationGroups() )
  {
    for( int i = 0; i < group->numRows(); ++i )
    {
      if( group->row(i) == row )
        return group;
    }
  }
  return nullptr;
}//groupOf(...)


DecayBatchCalcLocation *DecayBatchCalcWidget::addLocationGroup( const string &location )
{
  return m_nuclideRows->addNew<DecayBatchCalcLocation>( location );
}


void DecayBatchCalcWidget::connectRow( DecayBatchCalcNuclide *row )
{
  row->changed().connect( this, &DecayBatchCalcWidget::scheduleResultUpdate );
  row->removeRequested().connect( this, &DecayBatchCalcWidget::handleRemoveRow );
}//connectRow(...)


bool DecayBatchCalcWidget::hasLocationGroups() const
{
  for( DecayBatchCalcLocation * const group : locationGroups() )
  {
    if( !group->location().empty() )
      return true;
  }
  return false;
}//hasLocationGroups()


void DecayBatchCalcWidget::updateMixInputState()
{
  const bool grouped = hasLocationGroups();

  // Grouped input co-decays each location as its own mixture, so mixing is implied within a location
  //  and impossible across them.  Show that as a checked-but-disabled box rather than quietly
  //  ignoring whatever the user had set.
  if( grouped )
    m_mixInput->setChecked( true );
  m_mixInput->setEnabled( !grouped );

  // Only the text changes; whether the tooltip shows at all was decided by the constructor's
  //  HelpSystem::attachToolTipOn (which respects the user's "ShowTooltips" preference).
  m_mixInput->setToolTip( WString::tr( grouped ? "dbc-mix-input-grouped-tt" : "dbc-mix-input-tt" ),
                         Wt::TextFormat::XHTML );
}//updateMixInputState()


DecayBatchCalcLocation *DecayBatchCalcWidget::resetToSingleGroup()
{
  m_nuclideRows->clear();
  return addLocationGroup( string() );
}//resetToSingleGroup()


DecayBatchCalcNuclide *DecayBatchCalcWidget::addRowToGroup( DecayBatchCalcLocation *group )
{
  assert( group );
  DecayBatchCalcNuclide * const row = group->addRow( m_nuclideSuggest );
  row->setActivityText( defaultNewActivityText() );
  connectRow( row );
  return row;
}//addRowToGroup(...)


DecayBatchCalcNuclide *DecayBatchCalcWidget::addEmptyNuclideRow()
{
  // New rows join the last group, which for the ordinary (ungrouped) case is the only one.
  const vector<DecayBatchCalcLocation *> groups = locationGroups();
  return addRowToGroup( groups.empty() ? addLocationGroup( string() ) : groups.back() );
}//addEmptyNuclideRow()


string DecayBatchCalcWidget::defaultNewActivityText() const
{
  // If any existing row carries an activity unit, the new row joins that units style with "1 <unit>".
  for( DecayBatchCalcNuclide * const row : allRows() )
  {
    if( row->activityHasUnit() )
    {
      const string unit = row->activityUnitStr();
      if( !unit.empty() )
        return "1 " + unit;
    }
  }//for( each row )

  return "1";
}//defaultNewActivityText()


void DecayBatchCalcWidget::harmonizeActivityUnits()
{
  const vector<DecayBatchCalcNuclide *> rows = allRows();

  // Split valid rows into unit-bearing and bare-scalar, tracking their positions.
  vector<int> unitIdx, scalarIdx;
  vector<string> unitTok;
  for( size_t i = 0; i < rows.size(); ++i )
  {
    DecayBatchCalcNuclide * const row = rows[i];

    // A zero activity is a valid input (see isValid()), and its unit still labels the output - the
    //  core takes a location's unit from its first row - so it must be harmonized like any other.
    if( !row->nuclide() || !row->isValid() )
      continue;

    if( row->activityHasUnit() )
    {
      unitIdx.push_back( static_cast<int>(i) );
      unitTok.push_back( row->activityUnitStr() );
    }else
    {
      scalarIdx.push_back( static_cast<int>(i) );
    }
  }//for( each row )

  // Already consistent when everything is unit-bearing, or everything is a bare scalar.
  if( unitIdx.empty() || scalarIdx.empty() )
    return;

  // Give each bare scalar the unit of its nearest unit-bearing sibling (by row position), reflecting
  //  it in the GUI so the interpretation is explicit.
  for( const int s : scalarIdx )
  {
    int best = -1;
    for( size_t k = 0; k < unitIdx.size(); ++k )
    {
      if( (best < 0) || (std::abs(unitIdx[k] - s) < std::abs(unitIdx[best] - s)) )
        best = static_cast<int>( k );
    }

    if( best >= 0 )
    {
      DecayBatchCalcNuclide * const row = rows[s];
      string txt = row->activityText();
      SpecUtils::trim( txt );
      row->setActivityText( txt + " " + unitTok[best] );
    }
  }//for( each scalar row )
}//harmonizeActivityUnits()


void DecayBatchCalcWidget::handleRemoveRow( DecayBatchCalcNuclide *row )
{
  DecayBatchCalcLocation * const group = row ? groupOf( row ) : nullptr;
  if( !group )
    return;

  const vector<DecayBatchCalcLocation *> groups = locationGroups();

  // Always keep at least one row present: with nothing else left, blank the row instead of removing
  //  it (and drop the group's title, so the tool is back to its untitled default look).
  if( (groups.size() == 1) && (group->numRows() <= 1) )
  {
    row->setNuclide( nullptr, 0.0, false, 0.0, string() );
    group->setLocation( string() );
    updateMixInputState();
    scheduleResultUpdate();
    return;
  }

  // Emptying a group removes it - unless it is the last one, handled above.
  if( group->removeRow( row ) == 0 )
  {
    WidgetUtils::removeWidgetLater( m_nuclideRows, group );

    // That may have left a single group behind, which should look like the untitled default rather
    //  than keep a location title (`removeWidgetLater` has already detached it, so it is not counted).
    const vector<DecayBatchCalcLocation *> remaining = locationGroups();
    if( remaining.size() == 1 )
      remaining.front()->setLocation( string() );
  }

  updateMixInputState();
  scheduleResultUpdate();
}//handleRemoveRow(...)


void DecayBatchCalcWidget::addNuclide( const int z, const int a, const int iso,
                                       const double activity, const bool useCurie,
                                       const double age, const string &activityStr )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nuc = db ? db->nuclide( z, a, iso ) : nullptr;
  if( !nuc )
    return;

  const vector<DecayBatchCalcNuclide *> rows = allRows();

  // Re-use a lone empty row if present, otherwise append to the last group.
  DecayBatchCalcNuclide *row = nullptr;
  if( (rows.size() == 1) && !rows.front()->nuclide() )
  {
    row = rows.front();
  }else
  {
    const vector<DecayBatchCalcLocation *> groups = locationGroups();
    DecayBatchCalcLocation *group = groups.empty() ? addLocationGroup( string() ) : groups.back();
    row = group->addRow( m_nuclideSuggest );
    connectRow( row );
  }

  row->setNuclide( nuc, activity, useCurie, age, activityStr );
  scheduleResultUpdate();
}//addNuclide(...)


void DecayBatchCalcWidget::clearNuclides()
{
  resetToSingleGroup();
  addEmptyNuclideRow();

  // The cut is a plain number compared against whatever units the results are in, so keeping it
  //  across a change of inputs would silently filter a table it was never meant for.
  m_minActEdit->setText( WString() );

  updateMixInputState();
  scheduleResultUpdate();
}//clearNuclides()


void DecayBatchCalcWidget::scheduleResultUpdate()
{
  if( m_suppressUpdate )
    return;
  updateResult();
}//scheduleResultUpdate()


vector<DecayBatchCalc::BatchNuclide> DecayBatchCalcWidget::gatherInputs() const
{
  vector<DecayBatchCalc::BatchNuclide> inputs;
  for( DecayBatchCalcLocation * const group : locationGroups() )
  {
    for( int i = 0; i < group->numRows(); ++i )
    {
      DecayBatchCalcNuclide * const row = group->row( i );
      if( !row || !row->isValid() )
        continue;

      DecayBatchCalc::BatchNuclide bn;
      bn.nuclide = row->nuclide();
      bn.nuclide_str = bn.nuclide ? bn.nuclide->symbol : string();
      bn.age = row->age();
      bn.activity = row->activity();
      bn.unit_label = row->unitLabel();
      // A non-empty location puts the core into grouped/long output; the untitled default group
      //  leaves it empty, so the ordinary case is byte-for-byte what it always was.  The passthrough
      //  cells belong to the location, not the row (the core reads one set per group).
      bn.location = group->location();
      bn.product_suffix = group->productSuffix();
      bn.extra_columns = group->extraColumns();
      bn.activity_unit = group->activityUnit();
      inputs.push_back( bn );
    }//for( each row of the group )
  }//for( each group )

  return inputs;
}//gatherInputs()


bool DecayBatchCalcWidget::minActivityCut( double &cut ) const
{
  cut = 0.0;

  string cut_str = m_minActEdit->text().toUTF8();
  SpecUtils::trim( cut_str );
  if( cut_str.empty() )
    return true;   // blank means "no cut", which is a valid state

  // stod() alone would stop at the first junk character and silently turn "1e-3 Bq" into a cut of
  //  0.001, so require the whole (trimmed) text to be consumed.
  try
  {
    size_t used = 0;
    const double value = std::stod( cut_str, &used );
    if( (used != cut_str.size()) || IsNan(value) || IsInf(value) || (value < 0.0) )
      return false;
    cut = value;
  }catch( std::exception & )
  {
    return false;
  }

  return true;
}//minActivityCut(...)


DecayBatchCalc::BatchDecayOptions DecayBatchCalcWidget::gatherOptions() const
{
  DecayBatchCalc::BatchDecayOptions opts;

  const string time_str = m_timeEdit->text().toUTF8();
  try
  {
    opts.time_span = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( time_str, 0.0 );
  }catch( std::exception & )
  {
    opts.time_span = 0.0;
  }

  // Grouped output labels use the time exactly as typed (e.g. "48h"), so they match the source tool's.
  opts.time_span_str = time_str;
  SpecUtils::trim( opts.time_span_str );
  SpecUtils::erase_any_character( opts.time_span_str, " \t" );

  double cut = 0.0;
  opts.min_activity = minActivityCut( cut ) ? cut : 0.0;

  opts.num_steps = static_cast<size_t>( std::max( 1, m_stepsSpin->value() ) );
  opts.use_curie = (m_actUnitsCombo->currentIndex() == 0);
  opts.mix_input = m_mixInput->isChecked();
  opts.show_progeny = m_showProgeny->isChecked();
  opts.include_activity = m_incActivity->isChecked();
  opts.include_xrays = m_incXrays->isChecked();
  opts.include_gammas = m_incGammas->isChecked();
  opts.include_alphas = m_incAlphas->isChecked();
  opts.include_betas = m_incBetas->isChecked();

  return opts;
}//gatherOptions()


void DecayBatchCalcWidget::updateResult()
{
  m_resultContainer->clear();
  m_lastResult = DecayBatchCalc::BatchDecayResult();

  // Keep all rows to a single activity-units style before reading them (so cross-nuclide
  //  normalization is consistent).  setActivityText does not emit changed(), so this won't recurse.
  harmonizeActivityUnits();

  // Record an undo/redo step for any state change (nuclide/age/activity/option), captured after
  //  harmonization so the stored state matches what is shown.  Runs regardless of the validity
  //  early-returns below, so clearing a nuclide (making the state "empty") is still undoable.
  handleAddUndoPoint();

  const vector<DecayBatchCalc::BatchNuclide> inputs = gatherInputs();
  const DecayBatchCalc::BatchDecayOptions opts = gatherOptions();

  // A row whose nuclide text is non-empty but doesn't resolve to a decayable nuclide is an error; a
  //  blank row is just ignored.  Don't let a partial export slip out while such a row is present.
  bool hasInvalidRow = false;
  for( DecayBatchCalcNuclide * const row : allRows() )
  {
    string nucTxt = row->nuclideText();
    SpecUtils::trim( nucTxt );
    if( !nucTxt.empty() && !row->isValid() )
    {
      hasInvalidRow = true;
      break;
    }
  }//for( each row )

  const bool disableExport = inputs.empty() || hasInvalidRow;
  m_csvDownload->setHidden( disableExport );
  m_copyBtn->setHidden( disableExport );

  if( inputs.empty() )
  {
    m_resultContainer->addNew<WText>( WString::tr("dbc-no-inputs") );
    return;
  }

  // Zero (or unparseable) is rejected; a negative time is meaningful - it looks backwards.
  if( opts.time_span == 0.0 )
  {
    m_resultContainer->addNew<WText>( WString::tr("dbc-invalid-time") );
    m_csvDownload->hide();
    m_copyBtn->hide();
    return;
  }

  // Say so rather than computing with no cut: results that quietly ignore what the user typed look
  //  just like results that honoured it.
  double cut = 0.0;
  if( !minActivityCut( cut ) )
  {
    m_resultContainer->addNew<WText>( WString::tr("dbc-invalid-min-activity") );
    m_csvDownload->hide();
    m_copyBtn->hide();
    return;
  }

  try
  {
    m_lastResult = DecayBatchCalc::decay( inputs, opts );
  }catch( std::exception &e )
  {
    m_resultContainer->addNew<WText>( WString::tr("dbc-calc-error").arg( e.what() ) );
    m_csvDownload->hide();
    m_copyBtn->hide();
    return;
  }

  // Preview still shows the resolvable rows, but export stays disabled until the bad row is fixed.
  if( hasInvalidRow )
  {
    WText *warn = m_resultContainer->addNew<WText>( WString::tr("dbc-invalid-nuclide") );
    warn->addStyleClass( "DbcWarning" );
    warn->setInline( false );
  }

  if( !m_lastResult.warnings.empty() )
  {
    WText *warn = m_resultContainer->addNew<WText>( WString::fromUTF8(m_lastResult.warnings) );
    warn->addStyleClass( "DbcWarning" );
    warn->setInline( false );
  }

  // Safety cap on preview size; CSV/clipboard still work.
  if( m_lastResult.num_data_cells > sm_max_preview_cells )
  {
    m_resultContainer->addNew<WText>( WString::tr("dbc-too-large") );
    updateCopyToClipboardText();
    return;
  }

  // Build the HTML table.
  stringstream html;
  html << "<table class=\"DbcResultTable\">";
  html << "<thead><tr>";
  for( const string &h : m_lastResult.column_headers )
    html << "<th>" << Wt::Utils::htmlEncode(h) << "</th>";
  html << "</tr></thead><tbody>";
  for( const vector<string> &row : m_lastResult.rows )
  {
    html << "<tr>";
    for( size_t i = 0; i < row.size(); ++i )
      html << (i ? "<td>" : "<th>") << Wt::Utils::htmlEncode(row[i]) << (i ? "</td>" : "</th>");
    html << "</tr>";
  }
  html << "</tbody></table>";

  m_resultContainer->addNew<WText>( WString::fromUTF8(html.str()), Wt::TextFormat::XHTML );

  updateCopyToClipboardText();
}//updateResult()


void DecayBatchCalcWidget::handleAddUndoPoint()
{
  // Skip during the initial build (see #m_trackUndo) and while suppressing bulk updates.
  if( !m_trackUndo || m_suppressUpdate )
    return;

  const string current = encodeStateToUrl();
  const string prev = m_previousStateUri;

  // Nothing meaningful changed (e.g. an empty row was added, which doesn't affect the encoded state).
  if( current == prev )
    return;

  // Keep the baseline in sync with what is now shown, so the next change diffs correctly - including
  //  when we are mid-restore (undo/redo) and won't record a step below.
  m_previousStateUri = current;

  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( !undoRedo || undoRedo->isInUndoOrRedo() || !undoRedo->canAddUndoRedoNow() )
    return;

  // Restore a captured state by re-looking-up the tool (it may have been closed/reopened since) and
  //  replaying the URI.  createDecayBatchCalcWindow() returns the existing window when open, or a
  //  fresh one otherwise; it adds no nested step while an undo/redo is executing.
  auto restore = []( const string &uri ){
    InterSpec *viewer = InterSpec::instance();
    DecayBatchCalcWindow *win = viewer ? viewer->createDecayBatchCalcWindow() : nullptr;
    if( win )
      win->handleAppUrl( uri );
  };//restore

  undoRedo->addUndoRedoStep( [restore,prev](){ restore( prev ); },
                             [restore,current](){ restore( current ); },
                             "Batch decay state change" );
}//handleAddUndoPoint()


std::string DecayBatchCalcWidget::currentResultCsv() const
{
  return DecayBatchCalc::result_to_csv( m_lastResult );
}


void DecayBatchCalcWidget::updateCopyToClipboardText()
{
  const string csv = currentResultCsv();
  // Stash on the element's own `_isData` expando, which CopyDecayBatchToClipboard reads.
  m_copyBtn->doJavaScript( "var el=document.getElementById('" + m_copyBtn->id() + "');"
                           " el._isData = el._isData || {};"
                           " el._isData.TableData = " + Wt::WWebWidget::jsStringLiteral(csv, '\'') + ";" );
}//updateCopyToClipboardText()


void DecayBatchCalcWidget::tableCopiedToClipboardCallback( const int copied )
{
  if( copied )
  {
    passMessage( WString::tr("dbc-copied-ok"), WarningWidget::WarningMsgInfo );
  }else
  {
    passMessage( WString::tr("dbc-copied-fail"), WarningWidget::WarningMsgHigh );
  }
}//tableCopiedToClipboardCallback(...)


void DecayBatchCalcWidget::handleFileDrop( const string &/*display_name*/, const string &spool_name )
{
  try
  {
#ifdef _WIN32
    const std::wstring wname = SpecUtils::convert_from_utf8_to_utf16( spool_name );
    std::ifstream input( wname.c_str(), ios::in | ios::binary );
#else
    std::ifstream input( spool_name.c_str(), ios::in | ios::binary );
#endif
    if( !input )
      throw runtime_error( WString::tr("dbc-file-read-error").toUTF8() );

    stringstream buffer;
    buffer << input.rdbuf();
    loadCsvContents( buffer.str() );
  }catch( std::exception &e )
  {
    SimpleDialog *dialog = SimpleDialog::make( WString::tr("dbc-csv-error-title"), WString::fromUTF8(e.what()) );
    dialog->addButton( WString::tr("Close"), WidgetUtils::ButtonRole::Dismiss );
  }

  if( m_uploadResource )
    m_uploadResource->clearSpooledFiles();
}//handleFileDrop(...)


void DecayBatchCalcWidget::loadCsvContents( const string &contents )
{
  // Throws with a human-readable message on malformed input.
  const vector<DecayBatchCalc::BatchNuclide> parsed = DecayBatchCalc::parse_csv( contents );

  const bool useCurie = (m_actUnitsCombo->currentIndex() == 0);

  // Coalesce the per-row additions into a single recompute.
  m_suppressUpdate = true;
  m_nuclideRows->clear();

  // A cut from the previous inputs is meaningless against a new file's values (see clearNuclides()).
  m_minActEdit->setText( WString() );

  // One group per location, in first-seen order.  The source file is typically nuclide-major (all
  //  locations for one nuclide, then the next nuclide), so the rows of a location are gathered rather
  //  than assumed contiguous.  A file with no locations gives a single untitled group, which looks
  //  exactly like the hand-entry case.
  vector<DecayBatchCalcLocation *> groups;
  for( const DecayBatchCalc::BatchNuclide &bn : parsed )
  {
    if( !bn.nuclide )
      continue;

    DecayBatchCalcLocation *group = nullptr;
    for( DecayBatchCalcLocation * const g : groups )
    {
      if( g->location() == bn.location )
      {
        group = g;
        break;
      }
    }

    if( !group )
    {
      group = addLocationGroup( bn.location );
      group->setPassthrough( bn.product_suffix, bn.extra_columns, bn.activity_unit );
      groups.push_back( group );
    }

    DecayBatchCalcNuclide *row = group->addRow( m_nuclideSuggest );
    connectRow( row );

    // Show the value in the file's own unit ("uCi"), to enough digits to round-trip: the row text is
    //  what gatherInputs() re-parses, so a lossy rendering here would change the measurement (see
    //  #activity_text_in_unit).  A file giving no unit keeps the user's Ci/Bq preference.
    const string actStr = bn.activity_unit.empty()
      ? activity_text_in_unit( bn.activity, useCurie ? "Ci" : "Bq" )
      : activity_text_in_unit( bn.activity, bn.activity_unit );
    row->setNuclide( bn.nuclide, bn.activity, useCurie, bn.age, actStr );
    row->setUnitLabel( bn.unit_label );  // after setNuclide, which clears it
  }//for( each parsed nuclide )

  if( groups.empty() )
  {
    resetToSingleGroup();
    addEmptyNuclideRow();
  }

  updateMixInputState();
  m_suppressUpdate = false;
  updateResult();
}//loadCsvContents(...)


std::string DecayBatchCalcWidget::encodeStateToUrl() const
{
  string query = "time=" + m_timeEdit->text().toUTF8();
  query += "&steps=" + std::to_string( m_stepsSpin->value() );
  query += "&actunits=" + string( (m_actUnitsCombo->currentIndex() == 0) ? "ci" : "bq" );
  query += "&mix=" + string( m_mixInput->isChecked() ? "1" : "0" );
  query += "&progeny=" + string( m_showProgeny->isChecked() ? "1" : "0" );

  // Always emitted (even when blank), so restoring a state fully defines the field.  Url-encoded so a
  //  '+' in an exponent ("1E+20") survives - it would otherwise come back from urlDecode as a space -
  //  and so nothing typed into the field can break the query split.
  string cut = m_minActEdit->text().toUTF8();
  SpecUtils::trim( cut );
  query += "&mincut=" + Wt::Utils::urlEncode( cut );

  string inc;
  auto add_inc = [&inc]( const char *name ){ inc += (inc.empty() ? "" : ",") + string(name); };
  if( m_incActivity->isChecked() ) add_inc( "act" );
  if( m_incXrays->isChecked() )    add_inc( "xray" );
  if( m_incGammas->isChecked() )   add_inc( "gamma" );
  if( m_incAlphas->isChecked() )   add_inc( "alpha" );
  if( m_incBetas->isChecked() )    add_inc( "beta" );
  query += "&inc=" + inc;

  const vector<DecayBatchCalcLocation *> groups = locationGroups();

  // Encode the total row count (including blank/unresolved rows) so that adding or removing a row is
  // itself a state change for undo/redo, even before the row is filled out.  handleAppUrl pads with
  // blank rows to match.
  int num_rows = 0;
  for( DecayBatchCalcLocation * const group : groups )
    num_rows += group->numRows();
  query += "&rows=" + std::to_string( num_rows );

  for( DecayBatchCalcLocation * const group : groups )
  {
    // The ordinary hand-entry case is a single untitled group with nothing to say about itself, and
    //  encodes exactly as it did before groups existed - so URIs stay short, and old ones still read.
    const bool needs_header = (groups.size() > 1)
                              || !group->location().empty()
                              || !group->productSuffix().empty()
                              || !group->extraColumns().empty()
                              || !group->activityUnit().empty();

    if( needs_header )
    {
      // Group-level values are url-encoded individually: they come from a CSV, so they can hold
      //  spaces, '&', '=' or ',' (e.g. "Deposition at 28 hrs").  handleAppUrl decodes each one.
      query += "&loc=" + Wt::Utils::urlEncode( group->location() );
      query += "&grows=" + std::to_string( group->numRows() );
      if( !group->productSuffix().empty() )
        query += "&psuf=" + Wt::Utils::urlEncode( group->productSuffix() );
      if( !group->activityUnit().empty() )
        query += "&aunit=" + Wt::Utils::urlEncode( group->activityUnit() );
      for( const pair<string,string> &xc : group->extraColumns() )
      {
        query += "&xcol=" + Wt::Utils::urlEncode( xc.first );
        query += "&xval=" + Wt::Utils::urlEncode( xc.second );
      }
    }//if( needs_header )

    for( int i = 0; i < group->numRows(); ++i )
    {
      DecayBatchCalcNuclide * const row = group->row( i );
      const SandiaDecay::Nuclide * const nuc = row ? row->nuclide() : nullptr;
      if( !nuc )
        continue;
      query += "&nuc=" + nuc->symbol;
      // Encoded for the same reason the group values are: the activity text holds a space before its
      //  unit, and `label` comes straight from a CSV's Unit column.
      query += "&act=" + Wt::Utils::urlEncode( row->activityText() );
      const double age = row->age();
      if( age > 0.0 )
        query += "&initialage=" + PhysicalUnits::printToBestTimeUnits( age, 6 );
      // Carry any opaque areal/volumetric suffix (e.g. "/m2") so state round-trips losslessly.
      const string label = row->unitLabel();
      if( !label.empty() )
        query += "&label=" + Wt::Utils::urlEncode( label );
    }//for( each row of the group )
  }//for( each group )

  return "calc?" + query;
}//encodeStateToUrl()


void DecayBatchCalcWidget::handleAppUrl( std::string /*path*/, std::string query_str )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

  if( !query_str.empty() && (query_str[0] == '?') )
    query_str = query_str.substr(1);

  auto isNucKey = []( const string &key ) -> bool {
    return (key == "nuc") || (key == "nuclide") || (key == "iso") || (key == "isotope");
  };

  vector<string> args;
  SpecUtils::split( args, query_str, "&" );

  vector<pair<string,string>> fields;
  for( string arg : args )
  {
    SpecUtils::trim( arg );
    if( arg.empty() )
      continue;
    const string::size_type pos = arg.find('=');
    string key = (pos == string::npos) ? arg : arg.substr(0,pos);
    string value = (pos == string::npos) ? string() : arg.substr(pos+1);
    SpecUtils::trim( key );
    SpecUtils::trim( value );
    SpecUtils::to_lower_ascii( key );
    fields.push_back( { key, value } );
  }//for( each arg )

  bool useCurie = (m_actUnitsCombo->currentIndex() == 0);
  int desiredRows = 0;  // 0 == unspecified (older URLs); otherwise pad blank rows to this count.

  // First pass: options.
  for( const pair<string,string> &kv : fields )
  {
    const string &key = kv.first;
    string value = kv.second;
    if( key == "rows" )
    {
      try{ desiredRows = std::max(0, std::stoi(value)); }catch(...){}
    }else if( key == "mincut" )
    {
      m_minActEdit->setText( WString::fromUTF8( Wt::Utils::urlDecode(value) ) );
    }else if( key == "time" || key == "timespan" )
    {
      m_timeEdit->setText( value );
    }else if( key == "steps" )
    {
      try{ m_stepsSpin->setValue( std::max(1, std::stoi(value)) ); }catch(...){}
    }else if( key == "actunits" )
    {
      SpecUtils::to_lower_ascii( value );
      useCurie = !((value == "bq") || (value == "becquerel"));
      m_actUnitsCombo->setCurrentIndex( useCurie ? 0 : 1 );
    }else if( key == "mix" )
    {
      m_mixInput->setChecked( (value == "1") || (value == "true") || (value == "yes") );
    }else if( key == "progeny" )
    {
      m_showProgeny->setChecked( (value == "1") || (value == "true") || (value == "yes") );
    }else if( key == "inc" )
    {
      vector<string> parts;
      SpecUtils::split( parts, value, "," );
      auto has = [&parts]( const string &n ){
        for( const string &p : parts ){ if( SpecUtils::iequals_ascii(p,n) ) return true; }
        return false;
      };
      m_incActivity->setChecked( has("act") );
      m_incXrays->setChecked( has("xray") );
      m_incGammas->setChecked( has("gamma") );
      m_incAlphas->setChecked( has("alpha") );
      m_incBetas->setChecked( has("beta") );
    }
  }//for( first pass )

  // Second pass: location groups and their nuclides (act/age/label apply to the preceding nuc).  A
  //  "loc" key starts a group and is followed by that group's own passthrough cells; rows join the
  //  group most recently started, or an implicit untitled one - which is both the ordinary hand-entry
  //  case and what URLs written before groups existed look like.  Coalesce into one recompute.
  m_suppressUpdate = true;
  m_nuclideRows->clear();

  DecayBatchCalcLocation *group = nullptr;
  vector<int> groupRowCounts;  // Per group, in creation order; 0 == unspecified (see "grows").
  bool any_group_counts = false;  // Whether the URL gave per-group counts at all.

  // Puts a group in place for a row that arrived without the URL naming one.
  auto current_group = [this,&group,&groupRowCounts]() -> DecayBatchCalcLocation * {
    if( !group )
    {
      group = addLocationGroup( string() );
      groupRowCounts.push_back( 0 );
    }
    return group;
  };

  // Accumulated passthrough cells of the current group; reset by each "loc".
  string productSuffix, activityUnit;
  vector<pair<string,string>> extraColumns;

  for( size_t i = 0; i < fields.size(); ++i )
  {
    const string &key = fields[i].first;
    const string &value = fields[i].second;

    // Group-level values were url-encoded individually (they come from a CSV, so may hold spaces,
    //  '&', '=' or ','), so each needs decoding here.
    if( key == "loc" )
    {
      group = addLocationGroup( Wt::Utils::urlDecode(value) );
      groupRowCounts.push_back( 0 );
      productSuffix.clear();
      activityUnit.clear();
      extraColumns.clear();
      continue;
    }//if( key == "loc" )

    if( key == "grows" )
    {
      current_group();
      try
      {
        groupRowCounts.back() = std::max( 0, std::stoi(value) );
        any_group_counts = true;
      }catch( ... ){}
      continue;
    }//if( key == "grows" )

    if( (key == "psuf") || (key == "aunit") || (key == "xcol") || (key == "xval") )
    {
      if( key == "psuf" )
        productSuffix = Wt::Utils::urlDecode( value );
      else if( key == "aunit" )
        activityUnit = Wt::Utils::urlDecode( value );
      else if( key == "xcol" )
        extraColumns.push_back( { Wt::Utils::urlDecode(value), string() } );
      else if( !extraColumns.empty() )  //"xval" - always written right after its "xcol"
        extraColumns.back().second = Wt::Utils::urlDecode( value );

      current_group()->setPassthrough( productSuffix, extraColumns, activityUnit );
      continue;
    }//if( a group passthrough key )

    if( !isNucKey(key) )
      continue;

    const SandiaDecay::Nuclide * const nuc = db ? db->nuclide( value ) : nullptr;
    if( !nuc || nuc->isStable() )
      continue;

    string act_str, age_str, label_str;
    for( size_t j = i + 1; j < fields.size(); ++j )
    {
      if( isNucKey(fields[j].first) || (fields[j].first == "loc") )
        break;
      // `act` and `label` are url-encoded by encodeStateToUrl (an activity has a space before its
      //  unit, and a label comes from a CSV), so they decode like the group values.
      if( (fields[j].first == "act") || (fields[j].first == "activity") )
        act_str = Wt::Utils::urlDecode( fields[j].second );
      else if( (fields[j].first == "age") || (fields[j].first == "initialage") )
        age_str = fields[j].second;
      else if( fields[j].first == "label" )
        label_str = Wt::Utils::urlDecode( fields[j].second );
    }//for( trailing act/age/label )

    double act = 0.0, age = 0.0;
    try{ act = PhysicalUnits::stringToActivity( act_str ); }catch(...){}
    try{ age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( age_str, nuc->halfLife ); }catch(...){}

    DecayBatchCalcNuclide * const row = current_group()->addRow( m_nuclideSuggest );
    connectRow( row );
    row->setNuclide( nuc, act, useCurie, age, act_str );
    row->setUnitLabel( label_str );  // after setNuclide, which clears it
  }//for( second pass )

  // Pad blank rows so each group's row count - and the total - match the encoded state (see
  // encodeStateToUrl); keeps undo/redo of "add a blank row" solid.
  const vector<DecayBatchCalcLocation *> groups = locationGroups();
  for( size_t i = 0; (i < groups.size()) && (i < groupRowCounts.size()); ++i )
  {
    while( groups[i]->numRows() < groupRowCounts[i] )
      addRowToGroup( groups[i] );
  }

  // Older URLs (and the plain hand-entry case) carry only the overall "rows" count, so honour it when
  //  no per-group counts were given.  With them, they are authoritative: topping up to a "rows" total
  //  that disagrees would append the difference to the last group, silently moving rows between
  //  locations if the URI was truncated or hand-edited.
  const int minRows = any_group_counts ? 1 : std::max( 1, desiredRows );
  while( static_cast<int>(allRows().size()) < minRows )
    addEmptyNuclideRow();

  updateMixInputState();
  m_suppressUpdate = false;
  updateResult();
}//handleAppUrl(...)




DecayBatchCalcWindow::DecayBatchCalcWindow( InterSpec *viewer )
  : AuxWindow( WString::tr("dbc-window-title"),
              (AuxWindowProperties::DisableCollapse
               | AuxWindowProperties::EnableResize
               | AuxWindowProperties::SetCloseable) ),
    m_calc( nullptr )
{
  if( viewer )
    viewer->useMessageResourceBundle( "DecayBatchCalc" );

  WGridLayout *layout = stretcher();
  m_calc = layout->addWidget( make_unique<DecayBatchCalcWidget>( viewer ), 0, 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setRowStretch( 0, 1 );

  AuxWindow::addHelpInFooter( footer(), "decay-batch-dialog" );
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );

#if( USE_QR_CODES )
  WPushButton *qr_btn = footer()->addNew<WPushButton>();
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( this, [this](){
    try
    {
      const string url = "interspec://decaybatch/" + Wt::Utils::urlEncode( m_calc->encodeStateToUrl() );
      QrCode::displayTxtAsQrCode( url, WString::tr("dbc-window-title"), WString::tr("dbc-qr-window-text") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  } );
#endif //USE_QR_CODES

  rejectWhenEscapePressed();
  show();

  if( viewer && !viewer->isPhone() && (viewer->renderedWidth() > 100) && (viewer->renderedHeight() > 100) )
  {
    const int w = std::min( 550, viewer->renderedWidth() - 20 );
    const int h = std::min( 600, viewer->renderedHeight() );
    // A modest minimum height: below the natural size the inputs and results share the space (see
    //  DecayBatchCalc.css), and a larger floor would just have the dialog clip the tool's bottom.
    m_calc->setMinimumSize( std::min(480, w), std::min(320, h - 20) );
    resizeWindow( w, h );
    resizeToFitOnScreen();
    centerWindow();
  }else
  {
    // We get here when opened before the layout is known (e.g. restoring app state / a deep link on
    //  application start), where renderedWidth()/Height() aren't available yet.  Without an explicit
    //  minimum size the window collapses to its footer; give it a sensible floor (mirrors DecayWindow).
    //  The screen size is unknown here, so keep the floor low enough to fit a short display.
    if( viewer && !viewer->isPhone()
       && ((viewer->renderedWidth() <= 100) || (viewer->renderedHeight() <= 100)) )
    {
      m_calc->setMinimumSize( 480, 420 );
    }

    const bool wasPhone = m_isPhone;
    m_isPhone = false;
    if( wasPhone )
      resizeScaledWindow( 1.0, 1.0 ); //wouldnt have an effect if m_isPhone is true
    centerWindowHeavyHanded();        //wouldnt have an effect if m_isPhone is true
    m_isPhone = wasPhone;
  }
}//DecayBatchCalcWindow constructor


DecayBatchCalcWindow::~DecayBatchCalcWindow()
{
}


void DecayBatchCalcWindow::addNuclide( const int z, const int a, const int iso,
                                       const double activity, const bool useCurie,
                                       const double age, const string &activityStr )
{
  if( m_calc )
    m_calc->addNuclide( z, a, iso, activity, useCurie, age, activityStr );
}


void DecayBatchCalcWindow::clearNuclides()
{
  if( m_calc )
    m_calc->clearNuclides();
}


void DecayBatchCalcWindow::handleAppUrl( const string &url )
{
  string host, path, query, frag;
  AppUtils::split_uri( url, host, path, query, frag );
  handleAppUrl( path.empty() ? host : path, query );
}


void DecayBatchCalcWindow::handleAppUrl( const string &path, const string &query_str )
{
  if( m_calc )
    m_calc->handleAppUrl( path, query_str );
}


std::string DecayBatchCalcWindow::encodeStateToUrl() const
{
  return m_calc ? m_calc->encodeStateToUrl() : string();
}


void DecayBatchCalcWindow::loadCsvContents( const std::string &contents )
{
  if( m_calc )
    m_calc->loadCsvContents( contents );
}
