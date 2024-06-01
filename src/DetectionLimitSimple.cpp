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
#include <vector>
#include <sstream>


#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WMenuItem>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WStackedWidget>
#include <Wt/WDoubleValidator>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/GammaXsGui.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectionLimitSimple.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"

#if( USE_QR_CODES )
#include <Wt/Utils>

#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;


namespace
{
  
}//namespace



DetectionLimitSimpleWindow::DetectionLimitSimpleWindow( MaterialDB *materialDB,
                                Wt::WSuggestionPopup *materialSuggestion,
                                InterSpec *viewer )
: AuxWindow( WString::tr("window-title-simple-mda"),
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
             | AuxWindowProperties::SetCloseable
             | AuxWindowProperties::DisableCollapse) )
{
  rejectWhenEscapePressed( true );
  
  m_tool = new DetectionLimitSimple( materialDB, materialSuggestion, viewer, contents() );
  m_tool->setHeight( WLength(100,WLength::Percentage) );
  
  AuxWindow::addHelpInFooter( footer(), "simple-mda-dialog" );
  
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( footer() );
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://simple-mda/" + Wt::Utils::urlEncode(m_tool->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("dlsw-qr-tool-state-title"),
                                 WString::tr("dlsw-qr-tool-state-txt") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES

  
  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close") );
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  show();
  
  // If we are loading this widget, as we  are creating the InterSpec session,
  //  the screen width and height wont be available, so we'll just assume its
  //  big enough, which it should be.
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  int width = 625, height = 435;
  if( (screenW > 100) && (screenW < width) )
    width = screenW;
  if( (screenH > 100) && (screenH < height) )
    height = screenH;
  
  resizeWindow( width, height );
  
  // But I think this next call should fix things up, even if we do have a tiny screen
  resizeToFitOnScreen();
  
  centerWindowHeavyHanded();
}//DetectionLimitSimpleWindow(...) constructor


DetectionLimitSimpleWindow::~DetectionLimitSimpleWindow()
{
}


DetectionLimitSimple *DetectionLimitSimpleWindow::tool()
{
  return m_tool;
}





DetectionLimitSimple::DetectionLimitSimple( MaterialDB *materialDB,
                                 Wt::WSuggestionPopup *materialSuggestion,
                                 InterSpec *specViewer,
                                 Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
   m_viewer( specViewer ),
   m_materialSuggest( materialSuggestion ),
   m_materialDB( materialDB ),
   m_spectrum( nullptr ),
   m_stateUri()
{
  init();
}//DoseCalcWidget constructor


void DetectionLimitSimple::init()
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  wApp->useStyleSheet( "InterSpec_resources/DetectionLimitSimple.css" );
  m_viewer->useMessageResourceBundle( "DetectionLimitSimple" );
      
  addStyleClass( "DetectionLimitSimple" );
 
  m_spectrum = new D3SpectrumDisplayDiv( this );
  m_spectrum->addStyleClass( "SimpleMdaChart" );
  m_spectrum->setXAxisTitle( "" );
  m_spectrum->setYAxisTitle( "", "" );
  m_spectrum->setYAxisLog( false );
  m_spectrum->applyColorTheme( m_viewer->getColorTheme() );
  m_viewer->colorThemeChanged().connect( boost::bind( &D3SpectrumDisplayDiv::applyColorTheme, m_spectrum, boost::placeholders::_1 ) );
  m_spectrum->disableLegend();
  m_spectrum->setShowPeakLabel( SpectrumChart::PeakLabels::kShowPeakUserLabel, true );
  
  //shared_ptr<const SpecUtils::Measurement> hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  //m_spectrum->setData( hist, true );
  //m_spectrum->setXAxisRange( lower_lower_energy - 0.5*dx, upper_upper_energy + 0.5*dx );
  
  WContainerWidget *generalInput = new WContainerWidget( this );
  
  
}//void DetectionLimitSimple::init()


DetectionLimitSimple::~DetectionLimitSimple()
{
  //nothing to do here
}//~DoseCalcWidget()


void DetectionLimitSimple::handleAppUrl( std::string path, std::string query_str )
{
  //blah blah blah handle all this
  /*
#if( PERFORM_DEVELOPER_CHECKS )
  const string expected_uri = path + "?" + query_str;
#endif
  
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  Quantity calcQuantity;
  if( SpecUtils::iequals_ascii(path, "dose") )
    calcQuantity = Quantity::Dose;
  else if( SpecUtils::iequals_ascii(path, "act") )
    calcQuantity = Quantity::Activity;
  else if( SpecUtils::iequals_ascii(path, "dist") )
    calcQuantity = Quantity::Distance;
  else if( SpecUtils::iequals_ascii(path, "shield") )
    calcQuantity = Quantity::Shielding;
  else if( SpecUtils::iequals_ascii(path, "intro") )
  {
    handleQuantityClick( Quantity::NumQuantity );
    return;
  }else
    throw runtime_error( "Dose Calc tool: invalid URI path." );
  
    
  string inshielduri, outshielduri;
  const size_t shieldpos = query_str.find( "&INSHIELD=&" );
  
  if( shieldpos != string::npos )
  {
    inshielduri = query_str.substr(shieldpos + 11);
    query_str = query_str.substr(0, shieldpos);
  }
  
  size_t outshieldpos = inshielduri.find("&OUTSHIELD=&");
  if( outshieldpos != string::npos )
  {
    outshielduri = inshielduri.substr(outshieldpos + 12);
    inshielduri = inshielduri.substr(0, outshieldpos);
  }else if( (outshieldpos = query_str.find("&OUTSHIELD=&")) != string::npos )
  {
    outshielduri = query_str.substr(outshieldpos + 12);
    query_str = query_str.substr(0, outshieldpos);
  }
  
  if( inshielduri.empty() )
    m_enterShieldingSelect->setToNoShielding();
  else
    m_enterShieldingSelect->handleAppUrl( inshielduri );
  
  if( outshielduri.empty() )
    m_answerShieldingSelect->setToNoShielding();
  else
    m_answerShieldingSelect->handleAppUrl( outshielduri );
  
  SpecUtils::ireplace_all( query_str, "%23", "#" );
  SpecUtils::ireplace_all( query_str, "%26", "&" );
  SpecUtils::ireplace_all( query_str, "curries", "curies" ); //fix up me being a bad speller
  
  const map<string,string> parts = AppUtils::query_str_key_values( query_str );
  const auto ver_iter = parts.find( "VER" );
  if( ver_iter == end(parts) )
    Wt::log("warn") << "No 'VER' field in Dose Calc tool URI.";
  else if( ver_iter->second != "1" && !SpecUtils::starts_with(ver_iter->second, "1.") )
    throw runtime_error( "Can not read Dose Calc tool URI version '" + ver_iter->second + "'" );
  
  auto findUnitIndex = [&parts]( const string &key, Wt::WComboBox *combo ) -> int {
    const auto act_unit_iter = parts.find(key);
    if( act_unit_iter == end(parts) )
      return -1;
    
    for( int i = 0; i < combo->count(); ++i )
    {
      if( SpecUtils::iequals_ascii( combo->itemText(i).toUTF8(), act_unit_iter->second ) )
        return i;
    }
    assert( 0 );
    return -1;
  };//findUnitIndex
  
  const int act_in_unit_index = findUnitIndex( "ACTINUNIT", m_activityEnterUnits );
  const int act_out_unit_index = findUnitIndex( "ACTOUTUNIT", m_activityAnswerUnits );
  const int dose_in_unit_index = findUnitIndex( "DOSEINUNIT", m_doseEnterUnits );
  const int dose_out_unit_index = findUnitIndex( "DOSEOUTUNIT", m_doseAnswerUnits );
  
  const auto act_iter = parts.find("ACT");
  if( act_iter != end(parts) && !act_iter->second.empty() && (act_in_unit_index < 0) )
    throw runtime_error( "Dose Calc tool URI does not contain activity unit info." );
  
  const auto dose_iter = parts.find("DOSE");
  if( dose_iter != end(parts) && !dose_iter->second.empty() && (dose_in_unit_index < 0) )
    throw runtime_error( "Dose Calc tool URI does not contain dose unit info." );
  
  m_menu->select( static_cast<int>(calcQuantity) );
  handleQuantityClick( calcQuantity );
  
  if( act_in_unit_index >= 0 )
    m_activityEnterUnits->setCurrentIndex( act_in_unit_index );
  if( act_out_unit_index >= 0 )
    m_activityAnswerUnits->setCurrentIndex( act_out_unit_index );
  if( dose_in_unit_index >= 0 )
    m_doseEnterUnits->setCurrentIndex( dose_in_unit_index );
  if( dose_out_unit_index >= 0 )
    m_doseAnswerUnits->setCurrentIndex( dose_out_unit_index );
  
  const auto nuc_iter = parts.find("NUC");
  if( nuc_iter != end(parts) )
    m_gammaSource->setNuclideText( nuc_iter->second );
  
  const auto age_iter = parts.find("AGE");
  if( age_iter != end(parts) )
    m_gammaSource->setNuclideAgeTxt( age_iter->second );
  
  if( act_iter != end(parts) )
    m_activityEnter->setText( WString::fromUTF8(act_iter->second) );
  else
    m_activityEnter->setText( "" );
  
  if( dose_iter != end(parts) )
    m_doseEnter->setText( WString::fromUTF8(dose_iter->second) );
  else
    m_doseEnter->setText( "" );
  
  const auto dist_iter = parts.find("DIST");
  if( dist_iter != end(parts) )
    m_distanceEnter->setText( WString::fromUTF8(dist_iter->second) );
  else
    m_distanceEnter->setText( "" );
  
  updateResult();
  
#if( PERFORM_DEVELOPER_CHECKS )
  if( m_stateUri != expected_uri )
  {
    Wt::log("warn") << "DoseCalcWidget::handleAppUrl: input URI doesnt match current URI.\n\t input: '"
                    << expected_uri.c_str() << "'\n\tresult: '" << m_stateUri.c_str() << "'";
  }
#endif
   */
}//void handleAppUrl( std::string uri )


std::string DetectionLimitSimple::encodeStateToUrl() const
{
  /*
  // "interspec://dose/act?nuc=u238&dose=1.1ur/h&dist=100cm&..."
  
  string answer;
  
  switch( m_currentCalcQuantity )
  {
    case Dose:        answer += "dose";   break;
    case Activity:    answer += "act";    break;
    case Distance:    answer += "dist";   break;
    case Shielding:   answer += "shield"; break;
    case NumQuantity: answer += "intro"; break;
  }//switch( m_currentCalcQuantity )
  
  answer += "?VER=1";

  if( m_currentCalcQuantity == NumQuantity )
    return answer;
  
  // We could limit what info we put in the URL, based on current m_currentCalcQuantity,
  //  but we might as well put the full state into the URI.
  
  auto addTxtField = [&answer]( const string &key, Wt::WLineEdit *edit ){
    string txt = edit->text().toUTF8();
    SpecUtils::ireplace_all( txt, "#", "%23" );
    SpecUtils::ireplace_all( txt, "&", "%26" );
    answer += "&" + key + "=" + txt;
  };

  const SandiaDecay::Nuclide *nuc = m_gammaSource->nuclide();
  if( nuc )
    answer += "&NUC=" + nuc->symbol;
  
  if( nuc && !PeakDef::ageFitNotAllowed(nuc) )
    answer += "&AGE=" + m_gammaSource->nuclideAgeStr().toUTF8();
  
  addTxtField( "ACT", m_activityEnter );
  answer += "&ACTINUNIT=" + m_activityEnterUnits->currentText().toUTF8();
  answer += "&ACTOUTUNIT=" + m_activityAnswerUnits->currentText().toUTF8();
  
  addTxtField( "DOSE", m_doseEnter );
  answer += "&DOSEINUNIT=" + m_doseEnterUnits->currentText().toUTF8();
  answer += "&DOSEOUTUNIT=" + m_doseAnswerUnits->currentText().toUTF8();
  
  addTxtField( "DIST", m_distanceEnter );
  
  // We'll mark shielding URL starting with the below, and everything after this is the shielding.
  //  Having an order-independent method would be better, but for the moment...
  switch( m_currentCalcQuantity )
  {
    case Dose:
    case Activity:
    case Distance:
      if( !m_enterShieldingSelect->isNoShielding() )
        answer += "&INSHIELD=&" + m_enterShieldingSelect->encodeStateToUrl();
      break;
      
    case Shielding:
      if( !m_answerShieldingSelect->isNoShielding() )
        answer += "&OUTSHIELD=&" + m_answerShieldingSelect->encodeStateToUrl();
      break;
      
    case NumQuantity:
      break;
  }//switch( m_currentCalcQuantity )
  
  
  
  return answer;
   */
  //blah blah blah handle all this
  return "";
}//std::string encodeStateToUrl() const;


