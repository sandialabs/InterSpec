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
#include <sstream>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WServer>
#include <Wt/WComboBox>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WPopupMenu>
#include <Wt/WPushButton>
#include <Wt/WSplitButton>
#include <Wt/WApplication>
#include <Wt/WSuggestionPopup>
#include <Wt/WDoubleValidator>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakEdit.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeNameFilterModel.h"


using namespace Wt;
using namespace std;
using SpecUtils::Measurement;

namespace
{
  
  PeakEdit *get_session_peak_editor()
  {
    InterSpec *viewer = InterSpec::instance();
    PeakEditWindow *editWindow = viewer ? viewer->peakEdit() : nullptr;
    PeakEdit *edit = editWindow ? editWindow->peakEditor() : nullptr;
    assert( edit );
    return edit;
  }//PeakEdit *get_session_peak_editor()
  
}//namespace


PeakEditWindow::PeakEditWindow( const double energy,
                                PeakModel *peakmodel,
                                InterSpec *viewer )
  : AuxWindow( WString::tr("window-title-peak-editor"), AuxWindowProperties::PhoneNotFullScreen | AuxWindowProperties::DisableCollapse  )
{
  wApp->useStyleSheet( "InterSpec_resources/PeakEdit.css" );
  
  addStyleClass( "PeakEditWindow" );
  
  AuxWindow::addHelpInFooter( footer(), "peak-editor" );
  
  m_edit = new PeakEdit( energy, peakmodel, viewer, this );
  rejectWhenEscapePressed();
  
  if( viewer->isPhone() )
  {
    setMaximumSize( WLength::Auto, viewer->renderedHeight()-25 );
    contents()->setOverflow( WContainerWidget::OverflowAuto, Wt::Vertical );
  }else
  {
    resizeToFitOnScreen();
  }
  
  setClosable( true );
  setResizable( false );
  centerWindow();
  AuxWindow::show();
}//PeakEditWindow

PeakEdit *PeakEditWindow::peakEditor()
{
  return m_edit;
}

Wt::Signal<> &PeakEditWindow::editingDone()
{
  return m_edit->done();
}


PeakDef::CoefficientType PeakEdit::row_to_peak_coef_type( const PeakEdit::PeakPars type )
{
  switch( type )
  {
    case PeakPars::Mean:             return PeakDef::CoefficientType::Mean;
    case PeakPars::Sigma:            return PeakDef::CoefficientType::Sigma;
    case PeakPars::GaussAmplitude:   return PeakDef::CoefficientType::GaussAmplitude;
    case PeakPars::SkewPar0:  return PeakDef::CoefficientType::SkewPar0;
    case PeakPars::SkewPar1:       return PeakDef::CoefficientType::SkewPar1;
    case PeakPars::SkewPar2:      return PeakDef::CoefficientType::SkewPar2;
    case PeakPars::SkewPar3:      return PeakDef::CoefficientType::SkewPar3;
    case PeakPars::Chi2DOF:          return PeakDef::CoefficientType::Chi2DOF;
    case PeakPars::SigmaDrfPredicted:
    case PeakPars::RangeStartEnergy:
    case PeakPars::RangeEndEnergy:
    case PeakPars::OffsetPolynomial0:
    case PeakPars::OffsetPolynomial1:
    case PeakPars::OffsetPolynomial2:
    case PeakPars::OffsetPolynomial3:
    case PeakPars::OffsetPolynomial4:
    case PeakPars::SetContinuumToLinear:
    case PeakPars::PeakColor:
    case PeakPars::NumPeakPars:
      break;
  }//switch( type )
  
  throw logic_error( "Invalid conversion from PeakEdit::PeakPars to PeakDef::CoefficientType" );
  return PeakDef::CoefficientType::NumCoefficientTypes;
}//PeakDef::CoefficientType PeakEdit::row_to_peak_coef_type( const PeakEdit::PeakPars type )


PeakEdit::PeakEdit( const double energy,
                    PeakModel *peakmodel,
                    InterSpec *viewer,
                    AuxWindow *aux )
  : WContainerWidget( aux->contents() ),
    m_energy( energy ),
    m_peakModel( peakmodel ),
    m_viewer( viewer ),
    m_peakIndex(),
    m_originalPeak(),
    m_originalPeaks{},
    m_currentPeak(),
    m_blockInfoRefresh( false ),
    m_valueTable( NULL ),
    m_nuclide( NULL ),
    m_suggestions( NULL ),
    m_filterModel( NULL ),
    m_photoPeakEnergy( NULL ),
    m_userLabel( NULL ),
    m_color( nullptr ),
    m_applyColorForAllNuc( nullptr ),
    m_peakType( NULL ),
    m_continuumType( NULL ),
    m_skewType( NULL ),
    m_resetLinearContinuum( NULL ),
    m_apply( NULL ),
    m_accept( NULL ),
    m_cancel( NULL ),
    m_refit( NULL ),
    m_otherPeaksDiv( NULL ),
    m_otherPeakTxt( NULL ),
    m_prevPeakInRoi( NULL ),
    m_nextPeakInRoi( NULL ),
    m_drfFwhm( nullptr ),
    m_footer( aux->footer() ),
    m_aux( aux )
{
  init();
  changePeak( energy );
}//PeakEdit( constructor )


Wt::WString PeakEdit::rowLabel( const PeakPars t )
{
  switch( t )
  {
    case PeakEdit::Mean:              return WString::tr("pe-label-centroid");
    case PeakEdit::Sigma:             return WString::tr("FWHM");
    case PeakEdit::SigmaDrfPredicted: return WString::tr("pe-label-drf-fwhm");
    case PeakEdit::GaussAmplitude:    return WString::tr("pe-label-amp");
    case PeakEdit::SkewPar0:          return WString::tr("pe-label-skew0");
    case PeakEdit::SkewPar1:          return WString::tr("pe-label-skew1");
    case PeakEdit::SkewPar2:          return WString::tr("pe-label-skew2");
    case PeakEdit::SkewPar3:          return WString::tr("pe-label-skew3");
    case PeakEdit::OffsetPolynomial0: return WString::tr("pe-label-cont-p0");
    case PeakEdit::OffsetPolynomial1: return WString::tr("pe-label-cont-p1");
    case PeakEdit::OffsetPolynomial2: return WString::tr("pe-label-cont-p2");
    case PeakEdit::OffsetPolynomial3: return WString::tr("pe-label-cont-p3");
    case PeakEdit::OffsetPolynomial4: return WString::tr("pe-label-cont-p4");
    case PeakEdit::RangeStartEnergy:  return WString::tr("pe-label-roi-start");
    case PeakEdit::RangeEndEnergy:    return WString::tr("pe-label-roi-end");
    case PeakEdit::Chi2DOF:           return WString::tr("pe-label-chi2-dof");
    case PeakEdit::PeakColor:         return WString::tr("pe-label-peak-color");
    case PeakPars::SetContinuumToLinear:
    case PeakEdit::NumPeakPars:
      return "";
      break;
  }//case( t )

  return "";
}//const char *PeakEdit::rowLabel( const PeakPars t )


void PeakEdit::init()
{
  wApp->useStyleSheet( "InterSpec_resources/PeakEdit.css" );
  addStyleClass( "PeakEdit" );
  
  if( m_viewer )
    m_viewer->useMessageResourceBundle( "PeakEdit" );
  
  m_valueTable = new WTable( this );
  m_valueTable->setHeaderCount( 1, Horizontal );
  m_valueTable->addStyleClass( "PeakEditTable" );
  
  WLabel *label = new WLabel( WString::tr("pe-parameter"), m_valueTable->elementAt(0,0) );
  label = new WLabel( WString::tr("pe-value"), m_valueTable->elementAt(0,1) );
  label = new WLabel( WString::tr("pe-uncertainty"), m_valueTable->elementAt(0,2) );
  label = new WLabel( WString::tr("Fit"), m_valueTable->elementAt(0,3) );
  
  for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
  {
   if( (t == PeakEdit::PeakColor)
      || (t == PeakEdit::PeakPars::SigmaDrfPredicted)
      || (t == PeakEdit::PeakPars::SetContinuumToLinear) )
   {
     m_values[t] = m_uncertainties[t] = nullptr;
     m_fitFors[t] = nullptr;
     m_valueTable->rowAt(t+1)->setHidden(true);
     
     continue;
   }//if( t == PeakEdit::PeakColor )
    
    label = new WLabel( rowLabel(t), m_valueTable->elementAt(t+1,0) );
    m_values[t] = new WLineEdit( m_valueTable->elementAt(t+1,1) );
    m_uncertainties[t] = new WLineEdit( m_valueTable->elementAt(t+1,2) );
    
    m_values[t]->setAttributeValue( "ondragstart", "return false" );
    label->setBuddy( m_values[t] );
    m_uncertainties[t]->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_values[t]->setAttributeValue( "autocorrect", "off" );
    m_values[t]->setAttributeValue( "spellcheck", "off" );
    m_uncertainties[t]->setAttributeValue( "autocorrect", "off" );
    m_uncertainties[t]->setAttributeValue( "spellcheck", "off" );
#endif
    
    
    switch( t )
    {
      case PeakEdit::Mean: case PeakEdit::Sigma: case PeakEdit::GaussAmplitude:
      case PeakEdit::SkewPar0: case PeakEdit::SkewPar1:
      case PeakEdit::SkewPar2: case PeakEdit::SkewPar3:
      case PeakEdit::OffsetPolynomial0:
      case PeakEdit::OffsetPolynomial1: case PeakEdit::OffsetPolynomial2:
      case PeakEdit::OffsetPolynomial3: case PeakEdit::OffsetPolynomial4:
        m_fitFors[t] = new WCheckBox( m_valueTable->elementAt(t+1,3) );
      break;
      
      case PeakEdit::RangeStartEnergy: case PeakEdit::RangeEndEnergy:
      case PeakEdit::Chi2DOF: case PeakEdit::NumPeakPars:
        m_fitFors[t] = NULL;
      break;
        
      case PeakEdit::PeakPars::SigmaDrfPredicted:
      case PeakEdit::PeakPars::PeakColor:
      case PeakPars::SetContinuumToLinear:
        assert( 0 );
        break;
    }//case( t )
    
    
    
/*
  //Setting the ToolTip on the input form dowsnt seem to work
    txt = "";
    switch( t )
    {
      case PeakEdit::Mean:
        txt = "This uncertainty is used during fitting for calibration"
              " parameters, and during automatic isotope identification."
              " Altering this value will effect those algorithms.";
      break;
      case PeakEdit::Sigma:
        txt = "This uncertainty is used during fitting for the detector"
              " resolution response.";
      break;
      case PeakEdit::GaussAmplitude:
        txt = "The amplitude uncertainty is used when fitting for the source"
              " activity and shielding thicknesses.";
      break;
      case PeakEdit::SkewPar0:   case PeakEdit::SkewPar1:
      case PeakEdit::SkewPar2: case PeakEdit::SkewPar3:
      case PeakEdit::OffsetPolynomial0: case PeakEdit::OffsetPolynomial1:
      case PeakEdit::OffsetPolynomial2: case PeakEdit::OffsetPolynomial3:
      case PeakEdit::OffsetPolynomial4:
 case PeakEdit::PeakColor: case PeakEdit::PeakPars::SetContinuumToLinear:
      case PeakEdit::RangeStartEnergy:  case PeakEdit::RangeEndEnergy:
      case PeakPars::SigmaDrfPredicted:
      case PeakEdit::Chi2DOF: case PeakEdit::NumPeakPars:
        break;
    }//case( t )

    if( !txt.empty() )
      m_uncertainties[t]->setToolTip( txt );
*/
    
    
    switch( t )
    {
      case PeakEdit::Mean:
        m_values[t]->changed().connect( this, &PeakEdit::validateMeanOrRoiChange );
        m_values[t]->enterPressed().connect( this, &PeakEdit::validateMeanOrRoiChange );
      break;
        
      case PeakEdit::Sigma: case PeakEdit::GaussAmplitude:
      case PeakEdit::SkewPar0: case PeakEdit::SkewPar1:
      case PeakEdit::SkewPar2: case PeakEdit::SkewPar3:
      case PeakEdit::PeakPars::SigmaDrfPredicted:
      case PeakPars::SetContinuumToLinear:
      case PeakEdit::NumPeakPars:
        break;
      
      case PeakEdit::RangeStartEnergy: case PeakEdit::RangeEndEnergy:
        m_values[t]->changed().connect( this, &PeakEdit::validateMeanOrRoiChange );
        m_values[t]->enterPressed().connect( this, &PeakEdit::validateMeanOrRoiChange );
        //note, purposeful fall-through
        
      case PeakEdit::OffsetPolynomial0: case PeakEdit::OffsetPolynomial1:
      case PeakEdit::OffsetPolynomial2: case PeakEdit::OffsetPolynomial3:
      case PeakEdit::OffsetPolynomial4:
      case PeakEdit::Chi2DOF:
      case PeakEdit::PeakColor:
        m_uncertainties[t]->disable();
      break;
    }//case( t )
    
    WDoubleValidator *validator = new WDoubleValidator( m_values[t] );
    m_values[t]->setValidator( validator );
    m_values[t]->addStyleClass( "numberValidator"); //used to detect mobile keyboard
    validator = new WDoubleValidator( m_uncertainties[t] );
    m_uncertainties[t]->setValidator( validator );
    m_uncertainties[t]->addStyleClass( "numberValidator"); //used to detect mobile keyboard
//    m_values[t]->changed().connect( boost::bind( &PeakEdit::checkIfDirty, this, t, false ) );
//    m_uncertainties[t]->changed().connect( boost::bind( &PeakEdit::checkIfDirty, this, t, true ) );
    m_values[t]->textInput()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, false ) );
    m_values[t]->enterPressed()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, false ) );
    m_uncertainties[t]->textInput()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, true ) );
    m_uncertainties[t]->enterPressed()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, true ) );
    
    if( m_fitFors[t] )
    {
      m_fitFors[t]->checked().connect( boost::bind( &PeakEdit::fitTypeChanged, this, t ) );
      m_fitFors[t]->unChecked().connect( boost::bind( &PeakEdit::fitTypeChanged, this, t ) );
    }//if( m_fitFors[t] )
  }//for(...)

  WTableRow *row = nullptr;
  
  
  row = m_valueTable->rowAt( PeakEdit::PeakPars::SigmaDrfPredicted + 1 );
  row->elementAt(0)->setColumnSpan(4);
  m_drfFwhm = new WText( "", row->elementAt(0) );
  m_drfFwhm->addStyleClass( "PeakEditDrfFwhm" );
  row->setHidden( true );
  
  // When the detector changed/modified signal is emitted, the DRF for the foreground
  //  may not be updated yet, so we will delay until after current event loop
  auto drfUpdateWorker = wApp->bind( boost::bind(&PeakEdit::updateDrfFwhmTxt, this) );
  auto drfUpdater = [drfUpdateWorker](){
    WServer::instance()->schedule(100, wApp->sessionId(), [drfUpdateWorker](){
      drfUpdateWorker();
      wApp->triggerUpdate();
    });
  };
  m_viewer->detectorChanged().connect( std::bind(drfUpdater) );
  m_viewer->detectorModified().connect( std::bind(drfUpdater) );
  
  
  row = m_valueTable->rowAt( PeakEdit::PeakPars::SetContinuumToLinear + 1 );
  row->elementAt(0)->setColumnSpan(4);
  m_resetLinearContinuum = new WPushButton( WString::tr("pe-btn-cont-to-linear"), row->elementAt(0) );
  m_resetLinearContinuum->setToolTip( WString::tr("pe-tt-btn-cont-to-linear") );
  m_resetLinearContinuum->addStyleClass( "LightButton PeakEditEstLinCont" );
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer );
  HelpSystem::attachToolTipOn( m_resetLinearContinuum, WString::tr("pe-tt-btn-cont-to-linear"),
                              showToolTips, HelpSystem::ToolTipPosition::Right );
  row->setHidden( true );
  m_resetLinearContinuum->clicked().connect( boost::bind( &PeakEdit::estimateLinearContinuumFromData, this ) );
  
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+1 );
  label = new WLabel( WString::tr("Nuclide"), row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_nuclide = new WLineEdit( row->elementAt(1) );
  
  m_nuclide->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_nuclide->setAttributeValue( "autocorrect", "off" );
  m_nuclide->setAttributeValue( "spellcheck", "off" );
#endif
  
  string replacerJs, matcherJs;
  PhotopeakDelegate::EditWidget::replacerJs( replacerJs );
  PhotopeakDelegate::EditWidget::nuclideNameMatcherJs( matcherJs );
    
  m_suggestions = new WSuggestionPopup( matcherJs, replacerJs );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  m_suggestions->setJavaScriptMember("wtNoReparent", "true");
#endif
  
  m_suggestions->addStyleClass( "nuclide-suggest" );
  m_suggestions->forEdit( m_nuclide, WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );

  
//  std::shared_ptr<const PeakDef> dummypeak;
//  m_filterModel = new PeakIsotopeNameFilterModel( dummypeak, this );
  m_filterModel = new IsotopeNameFilterModel( this );
  m_filterModel->filter( "" );
  m_suggestions->setFilterLength( -1 );
  m_suggestions->setModel( m_filterModel );
//  m_suggestions->filterModel().connect( m_filterModel, &PeakIsotopeNameFilterModel::filter );
  m_suggestions->filterModel().connect( m_filterModel, &IsotopeNameFilterModel::filter );
  
  IsotopeNameFilterModel::setQuickTypeFixHackjs( m_suggestions );
  IsotopeNameFilterModel::setEnterKeyMatchFixJs( m_suggestions, m_nuclide );
  
  m_nuclide->enterPressed().connect( this, &PeakEdit::isotopeChanged );
  m_nuclide->blurred().connect( this, &PeakEdit::isotopeChanged );
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+2 );
  label = new WLabel( WString::tr("Photopeak"), row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_photoPeakEnergy = new WComboBox( row->elementAt(1) );
  m_photoPeakEnergy->setWidth( WLength(13.5,WLength::FontEm) );
  m_photoPeakEnergy->activated().connect( this, &PeakEdit::transitionChanged );
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+3 );
  label = new WLabel( WString::tr("pe-label"), row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_userLabel = new WLineEdit( row->elementAt(1) );
  
  m_userLabel->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_userLabel->setAttributeValue( "autocorrect", "off" );
  m_userLabel->setAttributeValue( "spellcheck", "off" );
#endif
  
  m_userLabel->setWidth( WLength(100,WLength::Percentage) );
  m_userLabel->blurred()
              .connect( boost::bind( &PeakEdit::checkIfUserLabelDirty, this ) );
  m_userLabel->enterPressed()
              .connect( boost::bind( &PeakEdit::checkIfUserLabelDirty, this ) );
  m_userLabel->textInput().connect( boost::bind( &PeakEdit::checkIfUserLabelDirty, this ) );
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+4 );
  label = new WLabel( WString::tr("pe-label-peak-color"), row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(1);
  m_color = new ColorSelect( ColorSelect::AllowNoColor, row->elementAt(1) );
  m_color->cssColorChanged().connect( boost::bind( &PeakEdit::checkIfColorDirty, this ) );
  row->elementAt(2)->setColumnSpan(2);
  m_applyColorForAllNuc = new WCheckBox( "dummy", row->elementAt(2) ); //"dummy" is needed or else later custom labels wont render; Wt bug?
  m_applyColorForAllNuc->setHidden( true );
  m_applyColorForAllNuc->setUnChecked();
  m_applyColorForAllNuc->checked().connect( boost::bind( &PeakEdit::checkIfColorDirty, this ) );
  m_applyColorForAllNuc->unChecked().connect( boost::bind( &PeakEdit::checkIfColorDirty, this ) );
  
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+5 );
  label = new WLabel( WString::tr("pe-peak-type"), row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  
  m_peakType = new WComboBox( row->elementAt(1) );
  m_peakType->addItem( WString::tr("pe-peak-model-gaussian") );    //PeakDef::GaussianDefined
  m_peakType->addItem( WString::tr("pe-peak-model-data") ); //PeakDef::DataDefined
  m_peakType->activated().connect( this, &PeakEdit::peakTypeChanged );
//

  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+6 );
  label = new WLabel( WString::tr("Continuum"), row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_continuumType = new WComboBox( row->elementAt(1) );
  m_continuumType->activated().connect( this, &PeakEdit::contnuumTypeChanged );
  
  for( PeakContinuum::OffsetType t = PeakContinuum::OffsetType(0);
       t <= PeakContinuum::External; t = PeakContinuum::OffsetType(t+1) )
  {
    m_continuumType->addItem( WString::tr(PeakContinuum::offset_type_label_tr(t)) );
  }//for( loop over PeakContinuum::OffsetType )

  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+7 );
  label = new WLabel( WString::tr("pe-skew-type"), row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_skewType = new WComboBox( row->elementAt(1) );
  m_skewType->activated().connect( this, &PeakEdit::skewTypeChanged );
  
  static_assert( PeakDef::DoubleSidedCrystalBall > PeakDef::ExpGaussExp, "SkewType ordering changed DSCB should be largest" );
  static_assert( PeakDef::DoubleSidedCrystalBall > PeakDef::CrystalBall, "SkewType ordering changed DSCB should be largest" );
  static_assert( PeakDef::DoubleSidedCrystalBall > PeakDef::GaussExp, "SkewType ordering changed DSCB should be largest" );
  static_assert( PeakDef::DoubleSidedCrystalBall > PeakDef::Bortel, "SkewType ordering changed DSCB should be largest" );
  
  for( PeakDef::SkewType t = PeakDef::SkewType(0);
      t <= PeakDef::DoubleSidedCrystalBall; t = PeakDef::SkewType(t+1) )
  {
    m_skewType->addItem( PeakDef::to_label(t) );
  }//for( loop over PeakDef::OffsetType )
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+8 );
  row->elementAt(0)->setColumnSpan(4);
  m_otherPeaksDiv = new WContainerWidget( row->elementAt(0) );
  m_otherPeaksDiv->addStyleClass( "PeakEditOtherPeaks" );
  m_otherPeakTxt = new WText( m_otherPeaksDiv );
  
  m_prevPeakInRoi = new WPushButton( m_otherPeaksDiv );
  m_prevPeakInRoi->setStyleClass( "PeakEditPreviousPeak" );
  m_prevPeakInRoi->clicked().connect( boost::bind( &PeakEdit::changeToNextPeakInRoi, this, true ) );
  
  m_nextPeakInRoi = new WPushButton( m_otherPeaksDiv );
  m_nextPeakInRoi->setStyleClass( "PeakEditNextPeak" );
  m_nextPeakInRoi->clicked().connect( boost::bind( &PeakEdit::changeToNextPeakInRoi, this, false ) );
  
  
  
  
  m_cancel = m_aux->addCloseButtonToFooter( WString::tr("Cancel"), false, m_footer );
  m_refit  = new WSplitButton( WString::tr("pe-btn-refit"), m_footer );
  m_apply  = new WPushButton( WString::tr("Apply"),  m_footer );
  m_accept = new WPushButton( WString::tr("Accept"), m_footer );
  //if( m_viewer && !m_viewer->isMobile() )
  //  m_accept->setIcon( "InterSpec_resources/images/accept.png" );
  
  //Add class to give padding on left side (or modify current style class)
  
  WPushButton *deleteButton = new WPushButton( WString::tr("Delete"), m_footer );
//  deleteButton->setFloatSide( Wt::Right );
  
  m_cancel->clicked().connect( this, &PeakEdit::cancel );
  m_refit->actionButton()->clicked().connect( boost::bind( &PeakEdit::refit, this, RefitOption::Normal ) );
  m_apply->clicked().connect(  this, &PeakEdit::apply  );
  m_accept->clicked().connect( this, &PeakEdit::accept );
  deleteButton->clicked().connect( this, &PeakEdit::deletePeak );
  
  m_peakModel->dataChanged().connect( this, &PeakEdit::refreshPeakInfo );
  m_peakModel->layoutChanged().connect( this, &PeakEdit::refreshPeakInfo );
  m_peakModel->rowsRemoved().connect( this, &PeakEdit::peakModelRowsRemoved );
  m_peakModel->rowsInserted().connect( this, &PeakEdit::peakModelRowsAdded );
  
  
  // Add ROI refit options to `m_refit`.
  if( m_viewer->isMobile() )
  {
    WPopupMenu *menu = new WPopupMenu();
    
    WMenuItem *item = menu->addItem( Wt::WString::tr("pe-btn-refit-roi-standard") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::Normal ) );
    item = menu->addItem( Wt::WString::tr("pe-btn-refit-roi-independent") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::RoiIndependentFwhms ) );
    item = menu->addItem( Wt::WString::tr("pe-btn-refit-roi-fine") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::RoiSmallRefinement ) );
    item = menu->addItem( Wt::WString::tr("pe-btn-refit-roi-independent-fine") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::RoiSmallRefinementIndependentFwhm ) );
    
    m_refit->setMenu( menu );
  }else
  {
    PopupDivMenu *menu = new PopupDivMenu( nullptr, PopupDivMenu::MenuType::TransientMenu);
    
    PopupDivMenuItem *item = menu->addMenuItem( Wt::WString::tr("pe-btn-refit-roi-standard") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::Normal ) );
    item = menu->addMenuItem( Wt::WString::tr("pe-btn-refit-roi-independent") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::RoiIndependentFwhms ) );
    item = menu->addMenuItem( Wt::WString::tr("pe-btn-refit-roi-fine") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::RoiSmallRefinement ) );
    item = menu->addMenuItem( Wt::WString::tr("pe-btn-refit-roi-independent-fine") );
    item->triggered().connect( boost::bind( &PeakEdit::refit, this, RefitOption::RoiSmallRefinementIndependentFwhm ) );
    
    m_refit->setMenu( menu );
  }//if( mobile ) / else
  
  if( m_refit->dropDownButton() )
    m_refit->dropDownButton()->hide();
}//void init()


PeakEdit::~PeakEdit()
{
  if( m_suggestions )
    delete m_suggestions;
}//~PeakEdit()


void PeakEdit::validateMeanOrRoiChange()
{
  try
  {
    if( m_values[Mean]->validate() != WValidator::Valid )
    {
      passMessage( "Invalid mean value - must be numeric.", WarningWidget::WarningMsgHigh );
      throw runtime_error( "" );
    }
    
    if( m_values[RangeStartEnergy]->validate() != WValidator::Valid )
    {
      passMessage( "Invalid ROI start value - must be numeric.", WarningWidget::WarningMsgHigh );
      throw runtime_error( "" );
    }
    
    if( m_values[Mean]->validate() != WValidator::Valid )
    {
      passMessage( "Invalid ROI end value - must be numeric.", WarningWidget::WarningMsgHigh );
      throw runtime_error( "" );
    }
    
    const double newmean = std::stod( m_values[Mean]->text().toUTF8() );
    const double newroilow = std::stod( m_values[RangeStartEnergy]->text().toUTF8() );
    const double newroiup = std::stod( m_values[RangeEndEnergy]->text().toUTF8() );
    
    if( newmean <= newroilow || newmean >= newroiup )
    {
      passMessage( "Peak mean must be between ROI start and end.", WarningWidget::WarningMsgHigh );
      throw runtime_error( "" );
    }
    
    if( newroilow >= newroiup )
    {
      //Probably cant get here, but JIC
      passMessage( "ROI start must be less than ROI end.", WarningWidget::WarningMsgHigh );
      throw runtime_error( "" );
    }
    
  }catch(...)
  {
    refreshPeakInfo();
  }//try / catch
}//void validateMeanOrRoiChange()


void PeakEdit::changeToNextPeakInRoi( const bool back )
{
  size_t peakn = 0;
  vector<std::shared_ptr<const PeakDef> > peaksInRoi;
  std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > allpeaks
                                                         = m_peakModel->peaks();
  
  if( !allpeaks )
    return;

  std::shared_ptr<const PeakContinuum> continuum = m_currentPeak.continuum();
  
  for( std::deque<PeakModel::PeakShrdPtr>::const_iterator i = allpeaks->begin();
        i != allpeaks->end(); ++i )
  {
    if( (*i)->continuum() == continuum )
    {
      if( fabs((*i)->mean() - m_currentPeak.mean()) < 0.1 )
        peakn = peaksInRoi.size();
      peaksInRoi.push_back( *i );
    }
  }//for( loop over peaks )
  
  const size_t npeak = peaksInRoi.size();
  
  if( peakn == 0 && back )
    return;
  
  if( peakn == (npeak-1) && !back )
    return;
  
  //If we're here, we will change peaks; lets check if we have changes to accept
  //  and if so, accept them.
  if( m_apply->isEnabled() )
    apply();
    
  peakn += (back ? -1 : +1);
  
  const double newmean = peaksInRoi[peakn]->mean();
  changePeak( newmean );
}//void changeToNextPeakInRoi( const bool back )


bool PeakEdit::isEditingValidPeak() const
{
  return m_peakIndex.isValid();
}//bool PeakEdit::isEditingValidPeak() const


void PeakEdit::changePeak( const double energy )
{
  WModelIndex nearestIndex;
  PeakModel::PeakShrdPtr nearPeak;
  double minDE = std::numeric_limits<double>::infinity();
  
  const int nrow = m_peakModel->rowCount();
  for( int row = 0; row < nrow; ++row )
  {
    WModelIndex index = m_peakModel->index( row, 0 );
    const PeakModel::PeakShrdPtr &peak = m_peakModel->peak( index );
    const double dE = fabs( peak->mean() - energy );
    if( dE < minDE )
    {
      minDE = dE;
      nearestIndex = index;
      nearPeak = peak;
    }//if( dE < minDE )
  }//for( int row = 0; row < nrow; ++row )
  
  double lowerx(0.0), upperx(0.0);
  if( nearPeak )
  {
    lowerx = nearPeak->lowerX();
    upperx = nearPeak->upperX();
  }
  
  if( nearPeak && (energy<lowerx || energy>upperx) )
    nearPeak.reset();

//  if( !nearPeak )
//  {
//    m_peakIndex = WModelIndex();
//    m_currentPeak.reset();
//    m_originalPeak.reset();
//    m_originalPeaks.clear();
//  }
  if( nearPeak )
  {
    m_energy = energy;
    m_peakIndex = nearestIndex;
    m_currentPeak = m_originalPeak = *nearPeak;
    
    m_originalPeaks.clear();
    const auto orig_peaks = m_peakModel->peaks();
    if( orig_peaks )
    {
      for( const auto &p : *orig_peaks )
        m_originalPeaks.push_back( p );
    }
    
    for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
    {
      m_valIsDirty[t] = m_uncertIsDirty[t] = false;
    }
    
    refreshPeakInfo();
  }//if( !nearPeak ) / else
}//void changePeak( double energy )


double PeakEdit::currentPeakEnergy() const
{
  return m_energy;
}//double currentPeakEnergy() const;


void PeakEdit::isotopeChanged()
{
  m_photoPeakEnergy->clear();
  string txt = m_nuclide->text().narrow();
  
  
  PeakDef::SourceGammaType gammaType;
  PeakDef::gammaTypeFromUserInput( txt, gammaType );
  
  SpecUtils::trim( txt );
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( txt );
  const SandiaDecay::Element *el = db->element( txt );

  if( nuc )
  {
    const double mean = m_currentPeak.mean();
    const double sigma = m_currentPeak.gausPeak() ? m_currentPeak.sigma() : 0.25*m_currentPeak.roiWidth();
    setNuclideFields( nuc, mean, sigma, NULL, -1, gammaType );
  }else if( el )
  {
    setNuclideFields( el, m_currentPeak.mean() );
  }else
  {
    const ReactionGamma *rctndb = ReactionGammaServer::database();
    vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
    try
    {
      rctndb->gammas( txt, possible_rctns );
    }catch(...){}
    //XXX - just taking first reaction, however there may be multiple
    const ReactionGamma::Reaction *rctn = NULL;
    if( possible_rctns.size() )
      rctn = possible_rctns[0].reaction;

    if( rctn )
    {
      setNuclideFields( rctn, m_currentPeak.mean() );
    }else
    {
      if( txt.empty() || SpecUtils::iequals_ascii(txt.c_str(), "none") )
      {
        m_photoPeakEnergy->clear();
        m_nuclide->setText("");
      }else
      {
        string nuclabel = "";
        if( m_currentPeak.parentNuclide() )
          nuclabel = m_currentPeak.parentNuclide()->symbol;
        m_nuclide->setText( nuclabel );
      }
      return;
    }//if( rctn ) / else
  }//if( !nuc )
  
  
  m_apply->enable();
  m_accept->enable();
}//void isotopeChanged()


void PeakEdit::transitionChanged()
{
  m_apply->enable();
  m_accept->enable();
}//void transitionChanged()


Wt::Signal<> &PeakEdit::done()
{
  return m_doneSignal;
}

void PeakEdit::peakModelRowsRemoved( WModelIndex, int firstRowToBeRemoved, int lastRowToBeRemoved )
{
  if( m_blockInfoRefresh )
    return;
  
  if( !m_peakIndex.isValid() )
    return;
  
  // We might be here if the user right-clicked on the peak and asked for a re-fit, which we would
  //  then, in principle want to pick the new peak up to edit.
  const int currentRow = m_peakIndex.row();
  const int nremoved = std::max( 1 + lastRowToBeRemoved - firstRowToBeRemoved, 1 );
  if( (currentRow >= firstRowToBeRemoved) && (currentRow <= lastRowToBeRemoved) )
  {
    m_peakIndex = WModelIndex();
    refreshPeakInfo();
  }else if( (currentRow > lastRowToBeRemoved) && (currentRow > 0) )
  {
    m_peakIndex = m_peakModel->index( currentRow - nremoved, m_peakIndex.column() );
  }
}//void peakModelRowsRemoved()


void PeakEdit::peakModelRowsAdded( Wt::WModelIndex , int /* firstRowInserted */ , int /* lastRowInserted */ )
{
  if( m_blockInfoRefresh )
    return;
  
  // We may be here if the user refit either the current peak, or any other peak by right-clicking
  //  on them (which causes all peaks to be removed, and then added back in).
  //  There are probably a number of other ways to get here.
  
  //const int nRowsInserted = std::max( 1 + (lastRowInserted - firstRowInserted), 1 );
  
  
  float mean = m_currentPeak.mean();
  if( mean < 1.0 )
    mean = m_energy;
  
  if( mean >= 1.0 )
  {
    PeakModel::PeakShrdPtr nearest = m_peakModel->nearestPeak( mean );
    
    double fwhm = 1.0;
    if( nearest )
      fwhm = nearest->gausPeak() ? nearest->fwhm() : 0.25*nearest->roiWidth();
    
    // Allowing new peak to be within 0.5*FWHM of original peak is an arbitrary choice, but should
    //  about uniquely identify the originally intended peak region
    if( nearest && (fabs(nearest->mean() - mean) < 0.5*fwhm) )
    {
      m_peakIndex = m_peakModel->indexOfPeak( nearest );
    }else
    {
      m_peakIndex = WModelIndex();
    }
  }else
  {
    m_peakIndex = WModelIndex();
  }
  
  refreshPeakInfo();
}//void PeakEdit::peakModelRowsAdded()


void PeakEdit::refreshPeakInfo()
{
  if( m_blockInfoRefresh )
    return;
  
  PeakModel::PeakShrdPtr updated_peak;
  
  try
  {
    if( m_peakIndex.isValid() )
      updated_peak = m_peakModel->peak( m_peakIndex );
  }catch( std::exception &e )
  {
    // We can get here if m_peakIndex.row() is out of range.
  }//try / catch
  
  
  //Check if updated_peak is the same mean as m_currentPeak, and if not, if check if we can
  //  find m_currentPeak in m_peakModel.  Some edge cases we are ignoring are: multiple peaks
  //  with same, or close together means, or when its the mean value that changes, and theres
  //  another nearby peak
  const double prev_mean = m_currentPeak.mean() > 1.0 ? m_currentPeak.mean() : m_energy;
  if( prev_mean >= 1.0 )
  {
    PeakModel::PeakShrdPtr nearest = m_peakModel->nearestPeak( prev_mean );
    if( nearest && (nearest != updated_peak) && (fabs(nearest->mean() - prev_mean) < 0.001) )
    {
      updated_peak = nearest;
    }
  }//if( prev_mean >= 1.0 )
  
  if( updated_peak )
    m_peakIndex = m_peakModel->indexOfPeak( updated_peak );
  
  
  if( !updated_peak )
  {
    m_peakIndex = WModelIndex();
    
    for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
    {
      m_valTxts[t] = m_uncertTxts[t] = "";
      if( m_values[t] )
        m_values[t]->setText( "" );
      if( m_uncertainties[t] )
        m_uncertainties[t]->setText( "" );
      if( m_fitFors[t] )
        m_fitFors[t]->setUnChecked();
    }//for( PeakPars t )
    
    m_nuclide->setText( "" );
    m_photoPeakEnergy->clear();
    m_userLabel->setText( "" );
    
    m_cancel->enable();
    m_refit->disable();
    m_apply->disable();
    m_accept->disable();
    
    m_otherPeaksDiv->hide();
    m_refit->dropDownButton()->hide();
    m_otherPeakTxt->setText( "" );
    
    m_drfFwhm->setText( "" );
    
    WTableRow *drfFwhmRow = m_valueTable->rowAt( PeakEdit::PeakPars::SigmaDrfPredicted );
    if( drfFwhmRow )
      drfFwhmRow->setHidden( true );
    
    return;
  }//if( !nrow )
  
  m_currentPeak = *updated_peak;
  
  std::shared_ptr<const PeakContinuum> continuum = m_currentPeak.continuum();
  
  for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
  {
    m_valIsDirty[t] = false;
    m_uncertIsDirty[t] = false;
    
    WTableRow *row = m_valueTable->rowAt( t + 1 );
    
    const int rownum = static_cast<int>(t);
    double val = 0.0, uncert = 0.0;
    
    switch( t )
    {
      case PeakPars::Mean: case PeakPars::Sigma: case PeakPars::GaussAmplitude:
      case PeakPars::SkewPar0: case PeakPars::SkewPar1:
      case PeakEdit::SkewPar2: case PeakEdit::SkewPar3:
      case PeakPars::Chi2DOF:
      {
        const PeakDef::CoefficientType ct = row_to_peak_coef_type( t );
        val = m_currentPeak.coefficient( ct );
        uncert = m_currentPeak.uncertainty( ct );
        break;
      }
        
      case PeakPars::SigmaDrfPredicted:
      case PeakPars::RangeStartEnergy: case PeakPars::RangeEndEnergy:
      case PeakPars::OffsetPolynomial0: case PeakPars::OffsetPolynomial1:
      case PeakPars::OffsetPolynomial2: case PeakPars::OffsetPolynomial3:
      case PeakPars::OffsetPolynomial4: case PeakPars::PeakColor:
      case PeakPars::SetContinuumToLinear:
      case PeakPars::NumPeakPars:
        break;
    }//switch( t )
    
    
    switch( t )
    {
      case PeakEdit::Mean:
      break;
      
      case PeakEdit::Sigma:
        val *= 2.3548201;
        uncert *= 2.3548201; //note purposeful fall through
        
      case PeakEdit::GaussAmplitude:
        switch( m_currentPeak.m_type )
        {
          case PeakDef::GaussianDefined:
            row->show();
            row->elementAt(1)->enable();
            row->elementAt(2)->enable();
            if( t == PeakEdit::GaussAmplitude )
            {
              val = m_currentPeak.peakArea();
              uncert = m_currentPeak.peakAreaUncert();
            }//if( t == PeakDef::GaussAmplitude )
          break;
          
          case PeakDef::DataDefined:
          {
            std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
            row->elementAt(1)->disable();
            row->elementAt(2)->disable();
            
            if( data )
            {
              assert( m_currentPeak.continuum()->energyRangeDefined() );
              if( m_currentPeak.continuum()->energyRangeDefined() )
              {
                const size_t lower_channel = data->find_gamma_channel( m_currentPeak.lowerX() );
                const size_t upper_channel = data->find_gamma_channel( m_currentPeak.upperX() - 0.00001 );
                const double dataArea = data->gamma_channels_sum( lower_channel, upper_channel );
                
                val = data->gamma_channels_sum( lower_channel, upper_channel );
                
                // \TODO: The next line I think would be equivalent (or maybe more correct) than
                //        previous line, but I think edge-cases need checking
                //val = data->gamma_integral( m_currentPeak.lowerX(), m_currentPeak.upperX() );
              }else
              {
                assert( 0 );
                val = -1;
              }
            }else
            {
              val = 0.0;
            }//if( data ) / else
            break;
          }//case PeakDef::DataDefined:
        }//switch( m_currentPeak.m_type )
      break;
        
      case PeakEdit::SkewPar0:
      case PeakEdit::SkewPar1:
      case PeakEdit::SkewPar2:
      case PeakEdit::SkewPar3:
      {
        const int i = t - PeakEdit::PeakPars::SkewPar0;
        const PeakDef::CoefficientType ct = row_to_peak_coef_type( t );
        auto type = m_currentPeak.skewType();
        double lower, upper, starting_val, step_size;
        if( PeakDef::skew_parameter_range( type, ct, lower, upper, starting_val, step_size ) )
        {
          row->setHidden( false );
          
          val = m_currentPeak.coefficient( ct );
          uncert = m_currentPeak.uncertainty( ct );
          
          if( (val < lower) || (val > upper) )
          {
            val = starting_val;
            uncert = 0.0;
          }
          
          m_values[t]->setHidden( false );
          auto validator = dynamic_cast<WDoubleValidator *>( m_values[t]->validator() );
          assert( validator );
          if( validator )
            validator->setRange( lower, upper );
        }else
        {
          row->setHidden( true );
        }
        
        break;
      }
        
      case PeakEdit::OffsetPolynomial0: case PeakEdit::OffsetPolynomial1:
      case PeakEdit::OffsetPolynomial2: case PeakEdit::OffsetPolynomial3:
      {
        const int coefnum = (t - PeakEdit::OffsetPolynomial0);
        const PeakContinuum::OffsetType type = continuum->type();
        
        bool hide = true;
        switch( type )
        {
          case PeakContinuum::NoOffset:   case PeakContinuum::External:
            hide = true;
            break;
            
          case PeakContinuum::Constant:   case PeakContinuum::Linear:
          case PeakContinuum::Quadratic:  case PeakContinuum::Cubic:
            hide = (coefnum >= type);
            break;
            
          case PeakContinuum::FlatStep:   case PeakContinuum::LinearStep:
          case PeakContinuum::BiLinearStep:
            hide = (coefnum >= (2 + (type - PeakContinuum::FlatStep)));
            break;
        }//switch( type )
        
        row->setHidden( hide );
        
        if( hide )
        {
          val = uncert = 0.0;
        }else
        {
          val = continuum->parameters().at( coefnum );
          uncert = continuum->uncertainties().at( coefnum );
        }//if( hide ) / else
        
        break;
      };//case (a polynomial coeffieienct )
        
      case PeakEdit::OffsetPolynomial4:
        row->setHidden( true );
        val = uncert = 0.0;
      break;
        
      case PeakEdit::PeakColor:
      {
        const WColor &c = m_currentPeak.lineColor();
        m_color->setColor( c );
        
        const bool couldPropogate = checkNuclideForDiffColors();
        m_applyColorForAllNuc->setChecked( false );
        m_applyColorForAllNuc->setHidden( !couldPropogate );
        break;
      }//case PeakEdit::PeakColor:
        
        
      case PeakEdit::RangeStartEnergy:
      case PeakEdit::RangeEndEnergy:
        uncert = 0.0;
        assert( continuum->energyRangeDefined() );
        
        if( continuum->energyRangeDefined() )
        {
          val = (t==RangeEndEnergy) ? continuum->upperEnergy()
                                    : continuum->lowerEnergy();
        }else
        {
          // We shouldnt ever get here, I dont think
          const shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);

          if( data )
          {
            const bool isHPGe = PeakFitUtils::is_likely_high_res( m_viewer );
            
            size_t lower_channel, upper_channel;
            estimatePeakFitRange( m_currentPeak, data, isHPGe, lower_channel, upper_channel );
            
            if( t == PeakEdit::RangeEndEnergy )
              val = data->gamma_channel_upper( upper_channel );
            else
              val = data->gamma_channel_lower( lower_channel );
          }//if( data )
        }//if( !m_currentPeak.xRangeDefined() )
      break;
      
      case PeakEdit::Chi2DOF:
        row->elementAt(1)->disable();
        row->elementAt(2)->disable();
        uncert = 0.0;
        row->setHidden( !m_currentPeak.chi2Defined() );
      break;
        
      case PeakPars::SigmaDrfPredicted:
        updateDrfFwhmTxt();
        break;
        
      case PeakPars::SetContinuumToLinear:
        row->setHidden( m_currentPeak.gausPeak() );
        break;
        
      case PeakEdit::NumPeakPars:
      break;
    }//case( t )
    
    if( m_fitFors[t] )
    {
      switch( t )
      {
        case PeakPars::Mean: case PeakPars::Sigma: case PeakPars::GaussAmplitude:
        case PeakPars::SkewPar0: case PeakPars::SkewPar1:
        case PeakPars::SkewPar2: case PeakPars::SkewPar3:
        case PeakPars::Chi2DOF:
        {
          const PeakDef::CoefficientType ct = row_to_peak_coef_type( t );
          m_fitFors[t]->setChecked( m_currentPeak.fitFor(ct) );
          break;
        }
          
        case PeakPars::RangeStartEnergy: case PeakPars::RangeEndEnergy:
        case PeakPars::SigmaDrfPredicted:
        case PeakPars::SetContinuumToLinear:
        case PeakPars::NumPeakPars:
          break;
          
        case PeakPars::OffsetPolynomial0: case PeakPars::OffsetPolynomial1:
        case PeakPars::OffsetPolynomial2: case PeakPars::OffsetPolynomial3:
        case PeakPars::OffsetPolynomial4: case PeakPars::PeakColor:
        {
          const int coefnum = rownum - static_cast<int>(PeakPars::OffsetPolynomial0);
          const vector<bool> fitfor = continuum->fitForParameter();
          if( coefnum >= 0 && static_cast<size_t>(coefnum) < fitfor.size() )
            m_fitFors[t]->setChecked( fitfor[coefnum] );
          break;
        }
      }//switch( t )
    }//if( m_fitFors[t] )
    
    char valTxt[32], uncertTxt[32];
    
    snprintf( valTxt, sizeof(valTxt), "%.4f", val );
    snprintf( uncertTxt, sizeof(uncertTxt), "%.4f", uncert );
    
    switch( t )
    {
      case PeakPars::Mean: case PeakPars::Sigma: case PeakPars::GaussAmplitude:
      case PeakPars::SkewPar0: case PeakPars::SkewPar1:
      case PeakPars::SkewPar2: case PeakPars::SkewPar3:
      case PeakPars::OffsetPolynomial0: case PeakPars::OffsetPolynomial1:
      case PeakPars::OffsetPolynomial2: case PeakPars::OffsetPolynomial3:
      case PeakPars::OffsetPolynomial4: case PeakPars::PeakColor:
        break;
    
      case PeakPars::Chi2DOF:
      case PeakPars::RangeStartEnergy: case PeakPars::RangeEndEnergy:
      case PeakPars::SigmaDrfPredicted:
      case PeakPars::SetContinuumToLinear:
      case PeakPars::NumPeakPars:
        uncertTxt[0] = '\0';
        break;
    }//switch( t )
    
    
    if( row->isHidden() )
      uncertTxt[0] = valTxt[0] = '\0';
    
    m_valTxts[t] = valTxt;
    m_uncertTxts[t] = uncertTxt;
    if( m_values[t] )
      m_values[t]->setText( valTxt );
    
    if( m_uncertainties[t] )
      m_uncertainties[t]->setText( uncertTxt );
  }//for( PeakPars t )
  
  m_userLabel->setText( WString::fromUTF8( m_currentPeak.userLabel() ) );
  
  const SandiaDecay::Nuclide *nuclide = m_currentPeak.parentNuclide();
  const SandiaDecay::Element *el = m_currentPeak.xrayElement();
  const ReactionGamma::Reaction *rctn = m_currentPeak.reaction();
  
  if( nuclide )
    setNuclideFields( nuclide,
                      m_currentPeak.mean(),
                      (m_currentPeak.gausPeak() ? m_currentPeak.sigma() : m_currentPeak.roiWidth()),
                      m_currentPeak.m_transition,
                      m_currentPeak.m_radparticleIndex,
                      m_currentPeak.sourceGammaType() );
  
  else if( el )
    setNuclideFields( el, m_currentPeak.xrayEnergy() );
  else if( rctn )
    setNuclideFields( rctn, m_currentPeak.reactionEnergy() );
  else
    setNuclideFields( nuclide, m_currentPeak.mean(),
                     (m_currentPeak.gausPeak() ? m_currentPeak.sigma() : m_currentPeak.roiWidth()),
                     m_currentPeak.m_transition, m_currentPeak.m_radparticleIndex,
                     m_currentPeak.sourceGammaType() );
  
  m_peakType->setCurrentIndex( m_currentPeak.type() );
  m_continuumType->setCurrentIndex( continuum->type() );
  m_skewType->setCurrentIndex( m_currentPeak.skewType() );
  
  m_valueTable->rowAt( PeakEdit::NumPeakPars+7 )->setHidden( m_currentPeak.type()==PeakDef::DataDefined );
  
  
  size_t thispeak = 0;
  vector<std::shared_ptr<const PeakDef> > peaksInRoi;
  std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > allpeaks
                                                        = m_peakModel->peaks();
  if( !!allpeaks )
  {
    for( std::deque< PeakModel::PeakShrdPtr >::const_iterator i = allpeaks->begin();
        i != allpeaks->end(); ++i )
    {
      if( (*i)->continuum() == continuum )
      {
        if( fabs((*i)->mean() - m_currentPeak.mean()) < 0.1 )
          thispeak = peaksInRoi.size();
        peaksInRoi.push_back( *i );
      }
    }
  }//if( !!allpeaks )
  
  if( peaksInRoi.size() < 2 )
  {
    m_otherPeaksDiv->hide();
    m_otherPeakTxt->setText( "" );
    m_refit->dropDownButton()->hide();
  }else
  {
    const size_t npeak = peaksInRoi.size();
    m_prevPeakInRoi->setHidden( thispeak == 0 );
    m_nextPeakInRoi->setHidden( thispeak == (npeak-1) );
    
    m_otherPeaksDiv->show();
    m_refit->dropDownButton()->show();
    WString txt = WString::tr("pe-multipeak-in-roi")
        .arg( static_cast<int>(thispeak+1) )
        .arg( static_cast<int>(npeak) );
    m_otherPeakTxt->setText( txt );
  }//if( peaksInRoi.size() < 2 ) / else
  
  updateSkewParameterLabels( m_currentPeak.skewType() );
  
  m_apply->disable();
//  m_accept->disable();
  m_accept->enable();
  m_cancel->enable();
  m_refit->setDisabled( m_currentPeak.type() != PeakDef::GaussianDefined );
}//void refreshPeakInfo()


void PeakEdit::setNuclideFields( const ReactionGamma::Reaction *rctn,
                                 const float mean )
{
  
  m_nuclide->setText( "" );
  m_photoPeakEnergy->clear();
  
  if( !rctn )
    return;
  
  const size_t ngammas = rctn->gammas.size();
  size_t nearest = 0;
  double nearestDE = 999999.9;
  
  for( size_t i = 0; i < ngammas; ++i )
  {
    const float thisDE = fabs( rctn->gammas[i].energy - mean );
    if( thisDE < nearestDE )
    {
      nearest = i;
      nearestDE = thisDE;
    }//if( thisDE < nearestE )
    
    char text[128];
    snprintf( text, sizeof(text), "%.4f keV I=%.1e",
              rctn->gammas[i].energy, rctn->gammas[i].abundance );
    m_photoPeakEnergy->addItem( text );
  }//for( size_t i = 0; i < nxray; ++i )
  
  if( nearestDE > 999999.0 )
  {
    m_photoPeakEnergy->clear();
    return;
  }//if( nearestDE > 999999.0 )
  
  m_nuclide->setText( rctn->name() );
  m_photoPeakEnergy->setCurrentIndex( static_cast<int>(nearest) );  
}//void setNuclideFields( ReactionGamma::Reaction *rctn, double mean )


void PeakEdit::setNuclideFields( const SandiaDecay::Element *el,
                                 const float mean )
{
  m_nuclide->setText( "" );
  m_photoPeakEnergy->clear();
  
  if( !el )
    return;
  
  const size_t nxray = el->xrays.size();
  size_t nearest = 0;
  float nearestDE = 999999.9f;
  
  for( size_t i = 0; i < nxray; ++i )
  {
    const float thisDE = fabs( el->xrays[i].energy - mean );
    if( thisDE < nearestDE )
    {
      nearest = i;
      nearestDE = thisDE;
    }//if( thisDE < nearestE )
    
    char text[128];
    snprintf( text, sizeof(text), "%.4f keV I=%.1e",
              el->xrays[i].energy, el->xrays[i].intensity );
    m_photoPeakEnergy->addItem( text );
  }//for( size_t i = 0; i < nxray; ++i )
  
  if( nearestDE > 999999.0 )
  {
    m_photoPeakEnergy->clear();
    return;
  }//if( nearestDE > 999999.0 )
  
  m_nuclide->setText( el->symbol + " xray" );
  m_photoPeakEnergy->setCurrentIndex( static_cast<int>(nearest) );
}//void PeakEdit::setNuclideFields( const SandiaDecay::Element *el, double mean )


void PeakEdit::setNuclideFields( const SandiaDecay::Nuclide *nuclide,
                                 const double mean, const double sigma,
                                 const SandiaDecay::Transition *transition,
                                 int particle_index,
                                 PeakDef::SourceGammaType sourceGammaType )
{
  m_nuclide->setText( "" );
  m_photoPeakEnergy->clear();
  
  if( !nuclide )
    return;

  if( (!transition || particle_index<0) && (sourceGammaType!=PeakDef::AnnihilationGamma) )
  {
    double extraEnergy = 0.0;
    switch( sourceGammaType )
    {
      case PeakDef::NormalGamma: break;
      case PeakDef::AnnihilationGamma: break;
      case PeakDef::SingleEscapeGamma: extraEnergy = 510.99891; break;
      case PeakDef::DoubleEscapeGamma: extraEnergy = 2.0*510.99891; break;
      case PeakDef::XrayGamma: break;
    }//switch( sourceGammaType )
    
    const bool xrayOnly = (sourceGammaType == PeakDef::XrayGamma);
    
    size_t transition_index = 0;
    PeakDef::findNearestPhotopeak( nuclide, mean + extraEnergy, 4.0*sigma, xrayOnly, -1.0,
                                    transition, transition_index, sourceGammaType );
    particle_index = static_cast<int>( transition_index );
  }//if( !transition || particle_index<0 )
  
  double bestEnergy = mean;
  
  switch( sourceGammaType )
  {
    case PeakDef::NormalGamma:
    case PeakDef::SingleEscapeGamma:
    case PeakDef::DoubleEscapeGamma:
    case PeakDef::XrayGamma:
      if( transition && particle_index < static_cast<int>(transition->products.size()) && particle_index>=0 )
        bestEnergy = transition->products[particle_index].energy;
    break;
    case PeakDef::AnnihilationGamma:
      bestEnergy = 510.99891*SandiaDecay::keV;
    break;
  }//switch( sourceGammaType )
  
  double minDE = DBL_MAX;
  int currentIndex = -1;
  const vector< IsotopeSelectionAids::NearGammaInfo > energies
             = IsotopeSelectionAids::equilibriumGammasByNearestEnergy( nuclide,
                                                        mean, 1.0E-10, true, true /*(sourceGammaType==PeakDef::XrayGamma)*/ );
  
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const IsotopeSelectionAids::NearGammaInfo &info = energies[i];
    
    double energy = info.gamma_energy;
    const double intensity = info.relative_intensity;
    
    const double dE = fabs( bestEnergy - energy );
    energy = floor(10000.0*energy + 0.5)/10000.0;
    
    if( dE < minDE )
    {
      minDE = dE;
      currentIndex = static_cast<int>( i );
    }
    
//    m_viewer->measurment(SpecUtils::SpectrumType::Foreground)->detector()
    
    char text[128];
    
    switch( info.gamma_type )
    {
      case PeakDef::NormalGamma: case PeakDef::AnnihilationGamma:
        snprintf( text, sizeof(text), "%.4f keV I=%.1e", energy, intensity );
      break;
      
      case PeakDef::SingleEscapeGamma:
        snprintf( text, sizeof(text), "S.E. %.4f keV", energy );
      break;
      
      case PeakDef::DoubleEscapeGamma:
        snprintf( text, sizeof(text), "D.E. %.4f keV", energy );
      break;
        
      case PeakDef::XrayGamma:
        snprintf( text, sizeof(text), "%.4f keV xray I=%.1e", energy, intensity );
      break;
    }//switch( info.gamma_type )
    
    m_photoPeakEnergy->addItem( text );
  }//for each( const double energy, energies )
  
  //PeakDef::findNearestPhotopeak() doesnt consider escape peaks, if the
  //  absolute closest peaj is a escape peak, and the chosen peak isnt too on,
  //  choose the escape peak
  if( energies.size()
      && (energies[0].gamma_type==PeakDef::SingleEscapeGamma
          || energies[0].gamma_type==PeakDef::DoubleEscapeGamma)
     && ((minDE/sigma) > 2.0) )
  {
    currentIndex = 0;
  }//
  
  m_nuclide->setText( nuclide->symbol );
  m_photoPeakEnergy->setCurrentIndex( currentIndex );
}//void setNuclideFields(...)


void PeakEdit::checkIfDirty( PeakPars t, bool uncert )
{
  if( uncert )
  {
    m_uncertIsDirty[t] = (m_uncertainties[t]->text() != m_uncertTxts[t]);
    if( m_uncertIsDirty[t] )
    {
      m_apply->enable();
      m_accept->enable();
    }//if( m_uncertIsDirty[type] )
  }else
  {
    m_valIsDirty[t] = (m_values[t] && (m_values[t]->text() != m_valTxts[t]));
    if( m_valIsDirty[t] )
    {
      m_apply->enable();
      m_accept->enable();
    }//if( m_valIsDirty[type] )
  }//if( uncert )
}//void checkIfDirty( PeakPars type, bool uncert )


bool PeakEdit::checkNuclideForDiffColors()
{
  const Wt::WColor c = m_color->color();
  
  //Only offer to apply this color to all peaks of a given source if there
  //  is more than one color for that source.  The logic to decide this
  //  isnt quite right yet.
  int nsamesrc = 0;
  auto allpeaks = m_peakModel->peaks();
  if( !allpeaks )
    return false;
  
  if( auto nuc = m_currentPeak.parentNuclide() )
  {
    m_applyColorForAllNuc->setText( WString::tr("pe-apply-all-of-type").arg(nuc->symbol) );
    for( const auto &p : *allpeaks )
      nsamesrc += (p->parentNuclide()==nuc && p->lineColor()!=c);
  }else if( auto el = m_currentPeak.xrayElement() )
  {
    m_applyColorForAllNuc->setText( WString::tr("pe-apply-all-of-type").arg(el->symbol) );
    for( const auto & p : *allpeaks )
      nsamesrc += (p->xrayElement()==el && p->lineColor()!=c);
  }else if( auto rctn = m_currentPeak.reaction() )
  {
    m_applyColorForAllNuc->setText( WString::tr("pe-apply-all-of-type").arg(rctn->name()) );
    for( const auto & p : *allpeaks )
      nsamesrc += (p->reaction()==rctn && p->lineColor()!=c);
  }else
  {
    m_applyColorForAllNuc->setText( WString::tr("pe-apply-all-of-type").arg("no src") );
    for( const auto & p : *allpeaks )
      nsamesrc += (!p->parentNuclide() && !p->xrayElement() && !p->reaction()
                   && p->lineColor()!=c);
  }
  
  if( m_originalPeak.lineColor() != c )
    nsamesrc -= 1;
  
  return nsamesrc>=1;
}//bool checkNuclideForDiffColors()


void PeakEdit::checkIfColorDirty()
{
  const Wt::WColor &origColor = m_originalPeak.lineColor();
  const Wt::WColor selectedColor = m_color->color();
  
  m_valIsDirty[PeakColor] = (m_applyColorForAllNuc->isChecked()
                             || (origColor!=selectedColor));
  
  const bool couldPropogate = checkNuclideForDiffColors();
  m_applyColorForAllNuc->setHidden( !couldPropogate );
  
  if( m_valIsDirty[PeakColor] )
  {
    m_apply->enable();
    m_accept->enable();
  }
}//void checkIfColorDirty()


void PeakEdit::checkIfUserLabelDirty()
{
  const string val = m_userLabel->text().toUTF8();
  const bool dirty = (val != m_originalPeak.userLabel());
  
  if( dirty )
  {
    m_apply->enable();
    m_accept->enable();
  }//if( dirty )
}//void checkIfUserLabelDirty()


void PeakEdit::fitTypeChanged( PeakPars t )
{
  if( !m_peakIndex.isValid() )
    return;
  
  if( !m_fitFors[t] )
    throw runtime_error( "PeakEdit::fitTypeChanged(...): logic error" );
  
  const bool isChecked = m_fitFors[t]->isChecked();
  std::shared_ptr<PeakContinuum> continuum = m_currentPeak.continuum();
  
  switch( t )
  {
    case Mean: case Sigma: case GaussAmplitude:
    case SkewPar0: case SkewPar1:
    case PeakPars::SkewPar2: case PeakPars::SkewPar3:
    case Chi2DOF:
    {
      const PeakDef::CoefficientType coef = row_to_peak_coef_type( t );
      m_peakModel->setPeakFitFor( m_peakIndex, coef, isChecked );
      //  m_originalPeak.setFitFor( coef, isChecked );
      //m_currentPeak.setFitFor( coef, isChecked );
      break;
    }//case for one of the PeakDef fitCoefficients
      
    case OffsetPolynomial0: case OffsetPolynomial1:
    case OffsetPolynomial2: case OffsetPolynomial3: case OffsetPolynomial4:
    {
      const int num = t - OffsetPolynomial0;
      m_peakModel->setContinuumPolynomialFitFor( m_peakIndex, num, isChecked );
      //continuum->setPolynomialCoefFitFor( num, isChecked );
      break;
    }
    
    case SigmaDrfPredicted:
    case PeakPars::SetContinuumToLinear:
    case RangeStartEnergy: case RangeEndEnergy:
    case PeakColor: case NumPeakPars:
    break;
  };//switch( t )
  
  m_currentPeak = *m_peakModel->peak( m_peakIndex );
  
//  m_apply->enable();
  m_accept->enable();
}//void fitTypeChanged( PeakPars t )


void PeakEdit::peakTypeChanged()
{
  PeakDef::DefintionType type
                        = PeakDef::DefintionType( m_peakType->currentIndex() );
  std::shared_ptr<const PeakContinuum> continuum = m_currentPeak.continuum();
  
  switch( type )
  {
    case PeakDef::GaussianDefined:
      m_refit->enable();
      for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
      {
        WTableRow *row = m_valueTable->rowAt(t+1);
        switch ( t )
        {
          case PeakEdit::Mean:              row->setHidden( false ); break;
          case PeakEdit::Sigma:             row->setHidden( false );  break;
          case PeakEdit::GaussAmplitude:    row->setHidden( false );  break;
          
          case PeakEdit::SkewPar0:
            row->setHidden( (m_currentPeak.skewType()==PeakDef::NoSkew) );
            break;
            
          case PeakEdit::SkewPar1:
            row->setHidden( (m_currentPeak.skewType()!=PeakDef::CrystalBall)
                           && (m_currentPeak.skewType()!=PeakDef::DoubleSidedCrystalBall)
                           && (m_currentPeak.skewType()!=PeakDef::ExpGaussExp) );
            break;
            
          case PeakEdit::SkewPar2:
          case PeakEdit::SkewPar3:
            row->setHidden( (m_currentPeak.skewType()!=PeakDef::DoubleSidedCrystalBall) );
            break;
            
          case PeakEdit::OffsetPolynomial0: case PeakEdit::OffsetPolynomial1:
          case PeakEdit::OffsetPolynomial2: case PeakEdit::OffsetPolynomial3:
          {
            const PeakContinuum::OffsetType type = continuum->type();
            const bool hide = (type==PeakContinuum::External)
                              || ((t-PeakEdit::OffsetPolynomial0)>=type);
            row->setHidden( hide );
            break;
          }//case( polynomial coefficient )
            
          case PeakEdit::OffsetPolynomial4:    row->setHidden( true ); break;

          case PeakEdit::RangeStartEnergy:     row->setHidden( false ); break;
          case PeakEdit::RangeEndEnergy:       row->setHidden( false ); break;
          case PeakEdit::Chi2DOF:
            row->setHidden( m_currentPeak.chi2Defined() );
          break;
            
          case PeakEdit::PeakColor:            row->setHidden( false ); break;
          case PeakPars::SigmaDrfPredicted:    row->setHidden( false ); break;
          case PeakPars::SetContinuumToLinear: row->setHidden( true );  break;
          case PeakEdit::NumPeakPars:                                   break;
        }//switch ( t )
      }//for( PeakPars t )
    break;
      
    case PeakDef::DataDefined:
    {
      m_refit->disable();
      
      for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
      {
        WTableRow *row = m_valueTable->rowAt(t+1);
        switch ( t )
        {
          case PeakEdit::Mean:              row->setHidden( false ); break;
          case PeakEdit::Sigma:             row->setHidden( true );  break;
          case PeakEdit::GaussAmplitude:    row->setHidden( true );  break;
          case PeakEdit::SkewPar0:
          case PeakEdit::SkewPar1:
          case PeakPars::SkewPar2:
          case PeakPars::SkewPar3:
            row->setHidden( true );
            break;
            
          case PeakEdit::OffsetPolynomial0:
          case PeakEdit::OffsetPolynomial1:
          case PeakEdit::OffsetPolynomial2:
          case PeakEdit::OffsetPolynomial3:
          case PeakEdit::OffsetPolynomial4:
//            row->setHidden( true );
          break;
            
          case PeakEdit::RangeStartEnergy:     row->setHidden( false ); break;
          case PeakEdit::RangeEndEnergy:       row->setHidden( false ); break;
          case PeakEdit::Chi2DOF:              row->setHidden( m_currentPeak.chi2Defined() ); break;
          case PeakEdit::PeakColor:            row->setHidden( true );  break;
          case PeakPars::SigmaDrfPredicted:    row->setHidden( true );  break;
          case PeakPars::SetContinuumToLinear: row->setHidden( false ); break;
          case PeakEdit::NumPeakPars: break;
        }//switch ( t )
      }//for( PeakPars t )    break;
    }//case PeakDef::DataDefined
  };//switch( PeakDef::DefintionType type )

  contnuumTypeChanged();
  skewTypeChanged();
  
  m_apply->enable();
  m_accept->enable();
}//void peakTypeChanged()


void PeakEdit::contnuumTypeChanged()
{
  PeakContinuum::OffsetType type = PeakContinuum::OffsetType( m_continuumType->currentIndex() );
  
  switch( type )
  {
    case PeakContinuum::NoOffset:
    case PeakContinuum::External:
//      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial0)->setHidden( true );
//      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial1)->setHidden( true );
//      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial2)->setHidden( true );
//      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial3)->setHidden( true );
//      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial4)->setHidden( true );
    break;
    
    case PeakContinuum::Constant:
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial0)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial1)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial2)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial3)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial4)->setHidden( true );
    break;
    
    case PeakContinuum::Linear:
    case PeakContinuum::FlatStep:
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial0)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial1)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial2)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial3)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial4)->setHidden( true );
    break;
  
    case PeakContinuum::Quadratic:
    case PeakContinuum::LinearStep:
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial0)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial1)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial2)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial3)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial4)->setHidden( true );
    break;
    
    case PeakContinuum::Cubic:
    case PeakContinuum::BiLinearStep:
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial0)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial1)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial2)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial3)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::OffsetPolynomial4)->setHidden( true );
    break;
  }//switch( type )
  
  m_apply->enable();
  m_accept->enable();
}//void contnuumTypeChanged()


void PeakEdit::updateSkewParameterLabels( const PeakDef::SkewType skewType )
{
  
  WLabel *label[4] = { nullptr, nullptr, nullptr, nullptr };
  for( int i = 0; i < 4; ++i )
  {
    WTableCell *cell = m_valueTable->elementAt(1+PeakEdit::SkewPar0+i, 0 );
    assert( cell && (cell->children().size() == 1) );
    if( !cell || cell->children().empty() )
      continue;
    
    label[i] = dynamic_cast<WLabel *>( cell->children()[0] );
    assert( label[i] );
  }//for( int i = 0; i < 4; ++i )
  
  
  switch( skewType )
  {
    case PeakDef::NoSkew:
      return;
      
    case PeakDef::SkewType::Bortel:
      if( label[0] )
        label[0]->setText( "Skew &alpha;" );
      break;
      
    case PeakDef::SkewType::CrystalBall:
      if( label[0] )
        label[0]->setText( "Skew &alpha;" );
      if( label[1] )
        label[1]->setText( "Skew n" );
      break;
      
    case PeakDef::SkewType::DoubleSidedCrystalBall:
      if( label[0] )
        label[0]->setText( "Skew &alpha;<sub>low</sub>" );
      if( label[1] )
        label[1]->setText( "Skew n<sub>low</sub>" );
      if( label[2] )
        label[2]->setText( "Skew &alpha;<sub>high</sub>" );
      if( label[3] )
        label[3]->setText( "Skew n<sub>high</sub>" );
      break;
      
    case PeakDef::SkewType::GaussExp:
      if( label[0] )
        label[0]->setText( "Skew k" );
      break;
      
    case PeakDef::SkewType::ExpGaussExp:
      if( label[0] )
        label[0]->setText( "Skew k<sub>L</sub>" );
      if( label[1] )
        label[1]->setText( "Skew k<sub>H</sub>" );
      break;
  }//switch( skewType )
}//void updateSkewParameterLabels();


void PeakEdit::skewTypeChanged()
{
  PeakDef::SkewType type = PeakDef::SkewType( m_skewType->currentIndex() );
  std::shared_ptr<const PeakContinuum> continuum = m_currentPeak.continuum();
  
  if( m_currentPeak.type() == PeakDef::DefintionType::DataDefined )
  {
    type = PeakDef::NoSkew;
    m_skewType->setCurrentIndex(PeakDef::NoSkew);
  }
  
  switch( type )
  {
    case PeakDef::NoSkew:
      m_valueTable->rowAt(1+PeakEdit::SkewPar0)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::SkewPar1)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::SkewPar2)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::SkewPar3)->setHidden( true );
    break;
      
    case PeakDef::SkewType::Bortel:
    case PeakDef::SkewType::CrystalBall:
    case PeakDef::SkewType::DoubleSidedCrystalBall:
    case PeakDef::SkewType::GaussExp:
    case PeakDef::SkewType::ExpGaussExp:
    {
      const int num_skew_pars = static_cast<int>(PeakDef::CoefficientType::SkewPar3)
                                - static_cast<int>(PeakDef::CoefficientType::SkewPar0);
      
      for( int i = 0; i <= num_skew_pars; ++i )
      {
        const auto index = static_cast<PeakEdit::PeakPars>( PeakEdit::PeakPars::SkewPar0 + i );
        const PeakDef::CoefficientType ct = row_to_peak_coef_type( index );
        
        double lower, upper, starting_val, step_size;
        if( PeakDef::skew_parameter_range( type, ct, lower, upper, starting_val, step_size ) )
        {
          m_values[index]->setHidden( false );
          m_valueTable->rowAt(1+index)->setHidden( false );
          m_values[index]->setText( SpecUtils::printCompact(starting_val, 4) );
          m_fitFors[index]->setChecked( true );
          auto validator = dynamic_cast<WDoubleValidator *>( m_values[index]->validator() );
          assert( validator );
          if( validator )
            validator->setRange( lower, upper );
        }else
        {
          m_values[index]->setHidden( true );
          m_uncertainties[index]->setText( "" );
          m_valueTable->rowAt(1+index)->setHidden( true );
        }
      }//for( int i = 0; i <= num_skew_pars; ++i )
      
      break;
    }//case Bortel, CrystalBall, DoubleSidedCrystalBall, GaussExp, ExpGaussExp
  }//switch( type )
  
  updateSkewParameterLabels( type );
  setSkewInputValueRanges( type );
  
  m_apply->enable();
  m_accept->enable();
}//void skewTypeChanged() const


bool PeakEdit::isDirty() const
{
  if( m_peakType->currentIndex() != m_currentPeak.type() )
    return true;
  
  std::shared_ptr<const PeakContinuum> continuum = m_currentPeak.continuum();
  if( m_continuumType->currentIndex() != continuum->type() )
    return true;
  
  if( m_skewType->currentIndex() != m_currentPeak.skewType() )
    return true;
  
  for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
    if( m_valIsDirty[t] || m_uncertIsDirty[t] )
      return true;

  if( nuclideInfoIsDirty() )
    return true;

  return (m_userLabel->text().toUTF8() != m_originalPeak.userLabel());
}//bool isDirty() const


void PeakEdit::updateDrfFwhmTxt()
{
  const double mean = m_currentPeak.mean();
  WTableRow *row = m_valueTable->rowAt( SigmaDrfPredicted + 1 );
  shared_ptr<const SpecMeas> meas = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const DetectorPeakResponse> drf = meas ? meas->detector() : nullptr;
  if( !drf || !drf->hasResolutionInfo() || (mean < 10.0) )
  {
    m_drfFwhm->setText( "" );
    row->setHidden( true );
  }else
  {
    const double drf_fwhm = drf->peakResolutionFWHM( mean );
    char text[32] = { '\0' };
    snprintf( text, sizeof(text), "%.2f", drf_fwhm );
    m_drfFwhm->setText( WString::tr("pe-drf-predict-fwhm-txt").arg(text) );
    row->setHidden( false );
  }//case PeakPars::SigmaDrfPredicted:
}//void updateDrfFwhmTxt();


void PeakEdit::setSkewInputValueRanges( const PeakDef::SkewType type )
{
  const int num_skew_pars = static_cast<int>(PeakDef::CoefficientType::SkewPar3)
                            - static_cast<int>(PeakDef::CoefficientType::SkewPar0);
  
  for( int i = 0; i <= num_skew_pars; ++i )
  {
    const auto index = static_cast<PeakEdit::PeakPars>( PeakEdit::PeakPars::SkewPar0 + i );
    const PeakDef::CoefficientType ct = row_to_peak_coef_type( index );
    
    double lower, upper, starting_val, step_size;
    if( PeakDef::skew_parameter_range( type, ct, lower, upper, starting_val, step_size ) )
    {
      m_values[index]->setHidden( false );
      m_valueTable->rowAt(index+1)->setHidden( false );
      //m_values[index]->setText( SpecUtils::printCompact(starting_val, 4) );
      
      auto validator = dynamic_cast<WDoubleValidator *>( m_values[index]->validator() );
      assert( validator );
      if( validator )
        validator->setRange( lower, upper );
      
      double val = 0.0;
      if( (m_values[index]->validate() != WValidator::Valid)
         || !(stringstream(m_values[index]->text().toUTF8()) >> val)
         || (val < lower) || (val > upper) )
      {
        m_values[index]->setValueText( SpecUtils::printCompact(starting_val, 4) );
      }
    }else
    {
      m_valueTable->rowAt(index+1)->setHidden( true );
      m_uncertainties[index]->setText( "" );
    }
  }//for( int i = 0; i <= num_skew_pars; ++i )
}//void setSkewInputValueRanges( const PeakDef::SkewType skewType )


void PeakEdit::estimateLinearContinuumFromData()
{
  const double ref_energy = m_currentPeak.mean();
  
  vector<PeakModel::PeakShrdPtr> prev_peaks;
  shared_ptr<const deque< PeakModel::PeakShrdPtr > > prev_peaks_ptr = m_peakModel->peaks();
  if( prev_peaks_ptr )
    prev_peaks.insert( end(prev_peaks), begin(*prev_peaks_ptr), end(*prev_peaks_ptr) );
  
  shared_ptr<const SpecUtils::Measurement> meas = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  if( !meas )
  {
    //Shouldnt really get here, I wouldnt think.
    passMessage( WString::tr("pe-est-lin-cont-no-data"), WarningWidget::WarningMsgHigh );
    return;
  }
  
  
  try
  {
    apply();
  }catch( std::exception &e )
  {
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
    return;
  }
  
  
  
  try
  {
    const double x0 = m_currentPeak.lowerX();
    const double x1 = m_currentPeak.upperX();
    
    const size_t start_channel      = meas->find_gamma_channel( x0 );
    const size_t end_channel        = meas->find_gamma_channel( x1 );
    
    const size_t num_side_bins = 3;
    double coefficients[2] = { 0.0, 0.0 };
    
    PeakContinuum::eqn_from_offsets( start_channel, end_channel, ref_energy,
                                    meas, num_side_bins, num_side_bins,
                                    coefficients[1], coefficients[0] );
    
    m_currentPeak.continuum()->setType( PeakContinuum::Linear );
    m_currentPeak.continuum()->setParameters( ref_energy, coefficients, nullptr );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("pe-failed-est-lin-cont").arg(e.what()), WarningWidget::WarningMsgHigh );
  }//try /catch
  
  refreshPeakInfo();
  
  UndoRedoManager *undo_manager = m_viewer->undoRedoManager();
  if( !undo_manager || !undo_manager->canAddUndoRedoNow() )
    return;
  
  vector<PeakModel::PeakShrdPtr> after_peaks;
  shared_ptr<const deque< PeakModel::PeakShrdPtr > > after_peaks_ptr = m_peakModel->peaks();
  if( after_peaks_ptr )
    after_peaks.insert( end(after_peaks), begin(*after_peaks_ptr), end(*after_peaks_ptr) );
  
  auto undo = [prev_peaks,ref_energy](){
    PeakEdit *edit = get_session_peak_editor();
    if( !edit )
      return;
    edit->m_peakModel->setPeaks( prev_peaks );
    edit->changePeak( ref_energy );
  };//undo
  
  auto redo = [after_peaks,ref_energy](){
    PeakEdit *edit = get_session_peak_editor();
    if( !edit )
      return;
    edit->m_energy = ref_energy;
    edit->m_peakModel->setPeaks( after_peaks );
    edit->changePeak( ref_energy );
  };
  
  undo_manager->addUndoRedoStep( undo, redo, "Est. linear continuum." );
}//void estimateLinearContinuumFromData()


bool PeakEdit::nuclideInfoIsDirty() const
{
  string currentTrans = m_photoPeakEnergy->currentText().narrow();
  PeakDef::SourceGammaType sourceType;
  PeakDef::gammaTypeFromUserInput( currentTrans, sourceType );
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( m_nuclide->text().narrow() );
  if( nuc != m_currentPeak.parentNuclide() )
    return true;
  
  if( nuc )
  {
    const SandiaDecay::Transition *transition = m_currentPeak.m_transition;
    const int rad_part_index = m_currentPeak.m_radparticleIndex;
  
    if( currentTrans.empty() && transition )
      return true;
    if( !transition && !currentTrans.empty() )
      return true;
    if( transition && m_currentPeak.sourceGammaType() != sourceType )
      return true;
  
    if( transition && (rad_part_index>=0) )
    {
      double energy;
      stringstream dblstrm( currentTrans );
    
      if( dblstrm >> energy )
      {
        //XXX - energies are good to 4 decimal places, we're assuming this is
        //      good eneough to uniqly identify gammas from tranitions, this
        //      could easily not be the case
        const double expectedEnergy = transition->products[rad_part_index].energy;
        if( fabs(energy-expectedEnergy) > 0.0006 )
          return true;
      }else
      {
        cerr << "Sholdnt be here dude, couldnt turn '" << currentTrans
             << "' into an energy." << endl;
        return true;
      }
    }//if( currentTrans.size() )
    
    return false;
  }//if( nuc )
  
  //Check to se if an element
  string nuctxt = m_nuclide->text().narrow();
  SpecUtils::ireplace_all( nuctxt, "xray", "" );
  SpecUtils::ireplace_all( nuctxt, "x-ray", "" );
  SpecUtils::trim( nuctxt );
  const SandiaDecay::Element *el = db->element( nuctxt );
  if( el )
  {
    if( el != m_currentPeak.xrayElement() )
      return true;
    
    double energy;
    if( !(stringstream(currentTrans) >> energy) )
    {
      cerr << "Sholdnt be here dude, couldnt turn '" << currentTrans
           << "' into an energy." << endl;
      return true;
    }//if( !(dblstrm >> energy) )
    
    return (fabs(energy-m_currentPeak.xrayEnergy()) > 0.0006f);
  }//if( el )

  //Check to se if a reaction
  const ReactionGamma *rctndb = ReactionGammaServer::database();
  vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
  try
  {
    rctndb->gammas( nuctxt, possible_rctns );
  }catch(...){}
  
  //XXX - just taking first reaction, however there may be multiple
  const ReactionGamma::Reaction *rctn = NULL;
  if( possible_rctns.size() )
    rctn = possible_rctns[0].reaction;

  if( rctn != m_currentPeak.reaction() )
    return true;
  
  if( rctn )
  {
    double energy;
    if( !(stringstream(currentTrans) >> energy) )
    {
      cerr << "Sholdnt be here dude, couldnt turn '" << currentTrans
            << "' into an energy." << endl;
      return true;
    }//if( !(dblstrm >> energy) )
    return (fabs(energy-m_currentPeak.reactionEnergy()) > 0.0006f);
  }//if( rctn )
  
  return false;
}//bool nuclideInfoIsDirty() const


void PeakEdit::refit( const PeakEdit::RefitOption type )
{
  try
  {
    apply();
  }catch( std::exception &e )
  {
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
    return;
  }
  
  vector<PeakDef> inputPeak, peaks_outside_roi, outputPeak;
  vector< std::shared_ptr<const PeakDef> > inpkptrs;
  
  const std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  const std::shared_ptr<const Measurement> continuum;
  std::shared_ptr<const SpecMeas> meas = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  
  
  if( !data || !meas )
    return;
  
  const PeakModel::PeakShrdPtr &thispeak = m_peakModel->peak(m_peakIndex);
//  inputPeak.push_back( *thispeak );
  
  const int npeak = static_cast<int>(m_peakModel->npeaks());
  for( int peak = 0; peak < npeak; ++peak )
  {
    WModelIndex index = m_peakModel->index( peak, 0 );
    const PeakModel::PeakShrdPtr &p = m_peakModel->peak(index);
    
    if( p->continuum() == thispeak->continuum() )
    {
      inputPeak.push_back( *p );
      inpkptrs.push_back( p );
    }else
    {
      peaks_outside_roi.push_back( *p );
    }
  }//for( int peak = 0; peak < npeak; ++peak )

  if( inputPeak.size() > 1 )  //JIC
    std::sort( inputPeak.begin(), inputPeak.end(), &PeakDef::lessThanByMean );
  
  //  const double lowE = peak->mean() - 0.1;
  //  const double upE = peak->mean() + 0.1;
  const double lowE = inputPeak.front().mean() - 0.1;
  const double upE = inputPeak.back().mean() + 0.1;
  const double ncausalitysigma = 0.0;
  const double stat_threshold  = 0.0;
  const double hypothesis_threshold = 0.0;
  
  
  if( inputPeak.size() > 1 )
  {
    const std::shared_ptr<DetectorPeakResponse> &detector
                                = m_viewer->measurment(SpecUtils::SpectrumType::Foreground)->detector();
    
    Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
    switch( type )
    {
      case RefitOption::Normal:
        // Use default fitting options
        break;
        
      case RefitOption::RoiIndependentFwhms:
        fit_options |= PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent;
        break;
        
      case RefitOption::RoiSmallRefinement:
        fit_options |= PeakFitLM::PeakFitLMOptions::SmallRefinementOnly;
        break;
        
      case RefitOption::RoiSmallRefinementIndependentFwhm:
        fit_options |= PeakFitLM::PeakFitLMOptions::SmallRefinementOnly;
        fit_options |= PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent;
        break;
    }//switch( type )
    
    
    const PeakShrdVec outp = refitPeaksThatShareROI( data, detector, inpkptrs, fit_options );
    for( size_t i = 0; i < outp.size(); ++i )
      outputPeak.push_back( *outp[i] );
  }else
  {
    assert( type == RefitOption::Normal ); //
    Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;  //No options - full refit...
    const bool isHPGe = PeakFitUtils::is_likely_high_res( m_viewer );
    outputPeak = fitPeaksInRange( lowE, upE, ncausalitysigma, stat_threshold,
                                  hypothesis_threshold, inputPeak, data, fit_options, isHPGe );
  }
  
  
  if( outputPeak.size() != inputPeak.size() )
  {
    passMessage( WString::tr("pe-failed-refit-insig"), WarningWidget::WarningMsgHigh );
    return;
  }//if( outputPeak.size() != inputPeak.size() )
  
  
  const double origEnergy = m_energy;
  const double newEnergy = (inputPeak.size() > 1) ? m_currentPeak.mean() : outputPeak[0].mean();
  vector<PeakModel::PeakShrdPtr> prev_peaks;
  shared_ptr<const deque< PeakModel::PeakShrdPtr > > prev_peaks_ptr = m_peakModel->peaks();
  if( prev_peaks_ptr )
    prev_peaks.insert( end(prev_peaks), begin(*prev_peaks_ptr), end(*prev_peaks_ptr) );

  vector<PeakDef> all_peaks = peaks_outside_roi;
  all_peaks.insert( end(all_peaks), begin(outputPeak), end(outputPeak) );
  std::sort( begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMean );

  m_energy = newEnergy;
  m_peakModel->setPeaks( all_peaks );
  changePeak( m_energy );
    
  auto undo = [prev_peaks,origEnergy](){
    PeakEdit *edit = get_session_peak_editor();
    if( !edit )
      return;
    edit->m_peakModel->setPeaks( prev_peaks );
    edit->changePeak( origEnergy );
  };//undo
    
  auto redo = [newEnergy,all_peaks](){
    PeakEdit *edit = get_session_peak_editor();
    if( !edit )
      return;

    edit->m_energy = newEnergy;
    edit->m_peakModel->setPeaks( all_peaks );
    edit->changePeak( newEnergy );
  };
    
  UndoRedoManager *undo_manager = m_viewer->undoRedoManager();
  if( undo_manager )
    undo_manager->addUndoRedoStep( undo, redo, "Peak refit." );
}//void refit()


void PeakEdit::setAmplitudeForDataDefinedPeak()
{
  std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  std::shared_ptr<PeakContinuum> continuum = m_currentPeak.continuum();
  if( !data || !continuum )
    return;
  
  const float roilower = static_cast<float>( continuum->lowerEnergy() );
  const float roiupper = static_cast<float>( continuum->upperEnergy() );
  const size_t lower_channel = data->find_gamma_channel( roilower );
  const size_t upper_channel = data->find_gamma_channel( roiupper );
  
  double peaksum = 0.0, datasum = 0.0;
  for( size_t channel = lower_channel; channel <= upper_channel; ++channel )
  {
    const float lowere = max( data->gamma_channel_lower(channel), roilower );
    const float uppere = min( data->gamma_channel_upper(channel), roiupper );
    const double dataarea = gamma_integral( data, lowere, uppere );
    const double contarea = continuum->offset_integral( lowere, uppere, data );
    datasum += max( 0.0, dataarea);
    peaksum += max( 0.0, (dataarea-contarea));
  }//for( int bin = lowerbin; bin <= upperbin; ++bin )
  
  m_currentPeak.set_coefficient( peaksum, PeakDef::GaussAmplitude );
  m_currentPeak.set_uncertainty( sqrt(datasum), PeakDef::GaussAmplitude );
}//void setAmplitudeForDataDefinedPeak()


void PeakEdit::apply()
{
  if( !isDirty() )
    return;
  
  // I dont think we need the `UndoRedoManager::PeakModelChange` since we will manage peaks
  // UndoRedoManager::PeakModelChange peak_undo_creator;
  
  PeakDef revertPeak; //used to restore peak if catches an exception.
  shared_ptr<const PeakDef> orig_peak_ptr;
    
  try
  {
    orig_peak_ptr = m_peakModel->peak( m_peakIndex );
    revertPeak = *orig_peak_ptr;
    
    // Grab a few things for undo/redo
    const double starting_energy = m_energy;
    const PeakDef starting_current_peak = m_currentPeak;
    vector< shared_ptr<const PeakDef> > starting_peaks;
    const auto start_peaks_ptr = m_peakModel->peaks();
    if( start_peaks_ptr )
      starting_peaks.insert( end(starting_peaks), begin(*start_peaks_ptr), end(*start_peaks_ptr) );
    
    
    // The PeakContinuum may be shared among multiple peaks, so if we change the continuum,
    //  we actually will need to make a copy of it, and change all affected peaks to use new one.
    //  If we dont make a copy of the continuum, then the copy of the original peaks that we
    //  use for undo/redo, will have their continuums changed, so the undo/redo wont actually
    //  appear to take place.
    //  The other thing we could have done, is make a deep-copy of the starting peaks, and use them
    //  for undo/redo, but doing it this way is possibly a little less error prone incase we
    //  have shared the PeakContinuum with some other peak, not in the PeakModel.
    //  Sharing the continuum between peaks to define a ROI is a bad design, but at this point
    //  we just got to go with it, until we can afford to make a large change of the ROI owning
    //  the peaks.
    shared_ptr<const PeakContinuum> orig_continuum_ptr = m_currentPeak.continuum();
    const PeakContinuum orig_continuum = *orig_continuum_ptr;
    shared_ptr<PeakContinuum> continuum = make_shared<PeakContinuum>( orig_continuum );
    
    const PeakContinuum::OffsetType offset = PeakContinuum::OffsetType(m_continuumType->currentIndex());
    if( offset != continuum->type() )
    {
      continuum->setType( offset );
      
      switch( offset )
      {
        case PeakContinuum::NoOffset:
          break;
          
        case PeakContinuum::External:
        {
          if( !!continuum->externalContinuum() )
            break;
          
          if( !!m_peakModel->peaks() )
          {
            for( const PeakModel::PeakShrdPtr &p : *m_peakModel->peaks() )
            {
              std::shared_ptr<const PeakContinuum> thiscont = p->continuum();
              if( !!thiscont->externalContinuum() )
              {
                continuum->setExternalContinuum( thiscont->externalContinuum() );
                break;
              }
            }//for( const PeakModel::PeakShrdPtr &p : *peaks )
          }//if( !!peaks && !continuum->externalContinuum() )
          
          if( !continuum->externalContinuum() )
          {
            std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
            std::shared_ptr<Measurement> background = estimateContinuum( data );
            continuum->setExternalContinuum( background );
          }//if( !continuum->externalContinuum() )
          
          break;
        }//case PeakContinuum::External:
          
        case PeakContinuum::Constant:
        case PeakContinuum::Linear:
        case PeakContinuum::Quadratic:
        case PeakContinuum::Cubic:
        case PeakContinuum::FlatStep:
        case PeakContinuum::LinearStep:
        case PeakContinuum::BiLinearStep:
          break;
      }//switch( offset )
    }//if( offset != m_currentPeak.offsetType() )
    
    
    
    const PeakDef::SkewType skewType = PeakDef::SkewType( m_skewType->currentIndex() );
    if( skewType != m_currentPeak.skewType() )
    {
      bool valid_skew = false;
      switch( skewType )
      {
        case PeakDef::SkewType::NoSkew:
        case PeakDef::SkewType::Bortel:
        case PeakDef::SkewType::GaussExp:
        case PeakDef::SkewType::CrystalBall:
        case PeakDef::SkewType::ExpGaussExp:
        case PeakDef::SkewType::DoubleSidedCrystalBall:
          valid_skew = true;
          break;
      }//switch( skewType )
      assert( valid_skew );
      if( !valid_skew )
        throw std::logic_error( "invalid skew type" );
      
      m_currentPeak.setSkewType( skewType );
      setSkewInputValueRanges( skewType );
      updateSkewParameterLabels( skewType );
      
      for( PeakPars t : {SkewPar0, SkewPar1, SkewPar2, SkewPar3} )
      {
        const PeakDef::CoefficientType ct = row_to_peak_coef_type( t );
      
        double lower, upper, starting_val, step_size;
        if( PeakDef::skew_parameter_range( skewType, ct, lower, upper, starting_val, step_size ) )
        {
          double val = 0.0, uncert = 0.0;
          if( !(stringstream(m_values[t]->text().toUTF8()) >> val)
             || (val < lower) || (val > upper) )
          {
            val = starting_val;
          }
          
          if( !(stringstream(m_uncertainties[t]->text().toUTF8()) >> uncert) )
            m_uncertainties[t]->setText( "0" );
          
          m_currentPeak.set_coefficient( val, ct );
          m_currentPeak.set_uncertainty( uncert, ct );
        }//if( use this parameter )
      }//for( loop over {SkewPar0, SkewPar1, SkewPar2, SkewPar3} )
      
    }//if( skewType != m_currentPeak.skewType() )
    
    
    for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
    {
      if( !m_values[t] && (t != PeakEdit::PeakColor) )
        continue;
      
      if( m_valIsDirty[t] )
      {
        double val = 0.0;
        
        if( t != PeakEdit::PeakColor )
        {
          if( m_values[t]->validate() != WValidator::Valid )
            throw runtime_error( "Value for '"
                                + rowLabel(t).toUTF8() + "' is not a valid number" );
          
          const string valtxt = m_values[t]->text().toUTF8();
          if( !(stringstream(valtxt) >> val) && !valtxt.empty() )
            throw runtime_error( "Error converting " + valtxt + " to float" );
        }//if( t != PeakEdit::PeakColor )
        
        switch( t )
        {
          case PeakEdit::GaussAmplitude:
            if( val < 0.0 )
              val = -val;
            m_currentPeak.setPeakArea( val );
            break;
            
          case PeakEdit::Sigma:
            val /= 2.3548201; //note: purposeful fall-through
            if( val < 0.0 )
              val = -val;
          
          case PeakEdit::Mean:
            m_energy = m_currentPeak.mean();
            //note: fall-through intentional
            
          case PeakEdit::SkewPar0:
          case PeakEdit::SkewPar1:
          case PeakEdit::SkewPar2:
          case PeakEdit::SkewPar3:
            m_currentPeak.set_coefficient( val, row_to_peak_coef_type(t) );
            break;
            
          case PeakEdit::OffsetPolynomial0:
          case PeakEdit::OffsetPolynomial1:
          case PeakEdit::OffsetPolynomial2:
          case PeakEdit::OffsetPolynomial3:
          case PeakEdit::OffsetPolynomial4:
            continuum->setPolynomialCoef( t-OffsetPolynomial0, val );
            if( m_currentPeak.type()==PeakDef::DataDefined )
              setAmplitudeForDataDefinedPeak();
            break;
            
          case PeakEdit::RangeStartEnergy:
          case PeakEdit::RangeEndEnergy:
          {
            double start(0.0), end(0.0);
            string starttxt = m_values[PeakEdit::RangeStartEnergy]->text().narrow();
            string endtxt = m_values[PeakEdit::RangeEndEnergy]->text().narrow();
            if( !(stringstream(starttxt) >> start) && !starttxt.empty() )
              throw runtime_error( "Error converting " + starttxt + " to float" );
            if( !(stringstream(endtxt) >> end) && !endtxt.empty() )
              throw runtime_error( "Error converting " + endtxt + " to float" );
            
            continuum->setRange( start, end );
            
            if( m_currentPeak.type() == PeakDef::DataDefined )
              setAmplitudeForDataDefinedPeak();
            
            break;
          }//case RangeStartEnergy/RangeEndEnergy
            
          case PeakEdit::PeakColor:
          {
            const WColor newColor = m_color->color();
            m_currentPeak.setLineColor( newColor );
            
            if( m_applyColorForAllNuc->isVisible() && m_applyColorForAllNuc->isChecked() )
            {
              m_applyColorForAllNuc->setUnChecked();
              
              std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > peaks = m_peakModel->peaks();
              if( peaks )
              {
                for( const auto &p : *peaks )
                {
                  if( p->parentNuclide()==m_currentPeak.parentNuclide()
                     && p->xrayElement()==m_currentPeak.xrayElement()
                     && p->reaction()==m_currentPeak.reaction()
                     && p->lineColor()!=newColor )
                  {
                    Wt::WModelIndex peakindex = m_peakModel->indexOfPeak( p );
                    if( peakindex.isValid() )
                    {
                      Wt::WModelIndex index = m_peakModel->index( peakindex.row(), PeakModel::kPeakLineColor );
                      m_peakModel->setData( index, WString(newColor.isDefault() ? "none" : newColor.cssText(false)) );
                    }
                  }//if( this is the same nuclide/element/reaction )
                }//for( const auto &p : *peaks )
              }//if( peaks )
            }//if( m_applyColorForAllNuc->isChecked() )
            
            break;
          }//case PeakEdit::PeakColor:
            
          case PeakEdit::Chi2DOF:
          case PeakPars::SigmaDrfPredicted:
          case PeakPars::SetContinuumToLinear:
          case PeakEdit::NumPeakPars:
            break;
        }//case( t )
      }//if( m_valIsDirty[t] )
      
      if( m_uncertIsDirty[t] )
      {
        double val = 0.0;
        
        if( m_uncertainties[t]->validate() != WValidator::Valid )
          throw runtime_error( "Uncertainty for '"
                              + rowLabel(t).toUTF8() + "' is not a valid number" );
        
        string valtxt = m_uncertainties[t]->text().narrow();
        if( !(stringstream(valtxt) >> val) && !valtxt.empty() )
          throw runtime_error( "Error converting " + valtxt + " to float" );
        
        if( val < 0.0 )
        {
          val = -val;
          //This will be propagated to the GUI by refreshPeakInfo();
          //        char uncertTxt[16];
          //        snprintf( uncertTxt, sizeof(uncertTxt), "%.4f", val );
          //        m_uncertainties[t]->setText( uncertTxt );
          
          
          
        }//if( val < 0.0 )
        
        switch( t )
        {
          case PeakEdit::GaussAmplitude:
            m_currentPeak.setPeakAreaUncert( val );
            break;
            
          case PeakEdit::Sigma:
            val /= 2.3548201; //note: purposfull fall-through
          case PeakEdit::Mean:
          case PeakEdit::SkewPar0:
          case PeakEdit::SkewPar1:
          case PeakEdit::SkewPar2:
          case PeakEdit::SkewPar3:
            m_currentPeak.set_uncertainty( val, row_to_peak_coef_type(t) );
            break;
            
          case PeakEdit::OffsetPolynomial0:
          case PeakEdit::OffsetPolynomial1:
          case PeakEdit::OffsetPolynomial2:
          case PeakEdit::OffsetPolynomial3:
          case PeakEdit::OffsetPolynomial4:
            continuum->setPolynomialUncert( t-OffsetPolynomial0, val );
            break;
            
          case PeakEdit::RangeStartEnergy:  case PeakEdit::RangeEndEnergy:
          case PeakEdit::Chi2DOF:           case PeakEdit::PeakColor:
          case PeakEdit::SigmaDrfPredicted: case PeakEdit::NumPeakPars:
          case PeakPars::SetContinuumToLinear:
            break;
        }//case( t )
      }//if( m_uncertIsDirty[t] )
    }//for( PeakPars t )
    
    //We have to set peak type after setting continuum parameters incase user
    //  changes a continuum type/polynomial-coefficient AND peak type
    if( m_peakType->currentIndex() != m_currentPeak.type() )
    {
      m_currentPeak.m_type = PeakDef::DefintionType(m_peakType->currentIndex());
      
      switch( m_peakType->currentIndex() )
      {
        case PeakDef::GaussianDefined:
          break;
          
        case PeakDef::DataDefined:
        {
          if( !continuum->energyRangeDefined() )
          {
            std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
            continuum->setRange( m_currentPeak.lowerX(), m_currentPeak.upperX() );
          }//if( !m_currentPeak.xRangeDefined() )
          
          const double middle = 0.5*(m_currentPeak.lowerX() + m_currentPeak.upperX() );
          m_currentPeak.set_coefficient( middle, PeakDef::Mean );
          m_currentPeak.set_uncertainty( 0.0, PeakDef::Mean );
          
          setAmplitudeForDataDefinedPeak();
          
          break;
        }//case PeakDef::DataDefined:
      }//switch( m_peakType->currentIndex() )
    }//if( m_peakType->currentIndex() != m_currentPeak.type() )
    
    
    m_currentPeak.setUserLabel( m_userLabel->text().toUTF8() );
    
    if( nuclideInfoIsDirty() )
    {
      string nuctxt = m_nuclide->text().narrow();
      
      PeakDef::SourceGammaType srcType;
      PeakDef::gammaTypeFromUserInput( nuctxt, srcType );
      
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      const SandiaDecay::Nuclide *nuc = db->nuclide( nuctxt );
      
      double energy;
      string currentTrans = m_photoPeakEnergy->currentText().narrow();
      PeakDef::SourceGammaType photopeakSourceType;
      PeakDef::gammaTypeFromUserInput( currentTrans, photopeakSourceType );
      
      switch( photopeakSourceType )
      {
        case PeakDef::NormalGamma:
        case PeakDef::AnnihilationGamma:
        case PeakDef::XrayGamma:
          break;
          
        case PeakDef::SingleEscapeGamma: case PeakDef::DoubleEscapeGamma:
          srcType = photopeakSourceType;
          break;
      }//switch( photopeakSourceType )
      
      stringstream dblstrm( currentTrans );
      
      if( !(dblstrm >> energy) )
      {
        //Shouldnt ever make it here
        cerr << "Misread '" << currentTrans << "' as " << energy << endl;
        energy = -1.0;
      }
      
      if( nuc )
      {
        if( energy < 0.0 )
          nuc = NULL;
        
        size_t transition_index = 0;
        const SandiaDecay::Transition *transition = NULL;
        
        const bool xrayOnly = (srcType == PeakDef::XrayGamma);
        
        PeakDef::SourceGammaType nearestGammaType;
        PeakDef::findNearestPhotopeak( nuc, energy, 0.0, xrayOnly, -1.0,
                                      transition, transition_index, nearestGammaType );
        
        switch( srcType )
        {
          case PeakDef::NormalGamma:
          case PeakDef::AnnihilationGamma:
          case PeakDef::XrayGamma:
            break;
            
          case PeakDef::SingleEscapeGamma:
          case PeakDef::DoubleEscapeGamma:
            nearestGammaType = srcType;
            break;
        }//switch( srcType )
        
        m_currentPeak.setNuclearTransition( nuc, transition,
                                           static_cast<int>(transition_index),
                                           nearestGammaType );
      }else if( !nuc )
      {
        m_currentPeak.setNuclearTransition( NULL, NULL, -1, PeakDef::NormalGamma );
        //try for an xray
        if( nuctxt.find_first_of( "0123456789" ) == string::npos )
        {
          SpecUtils::ireplace_all( nuctxt, "xray", "" );
          SpecUtils::ireplace_all( nuctxt, "x-ray", "" );
          SpecUtils::trim( nuctxt );
          const SandiaDecay::Element *el = db->element( nuctxt );
          const SandiaDecay::EnergyIntensityPair *nearpair
          = PeakDef::findNearestXray( el, energy );
          
          if( nearpair )
          {
            m_currentPeak.setXray( el, nearpair->energy );
          }else
          {
            //Try for a reaction
            const ReactionGamma *rctndb = ReactionGammaServer::database();
            vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
            try
            {
              rctndb->gammas( nuctxt, possible_rctns );
            }catch(...){}
            
            //XXX - just taking first reaction, however there may be multiple
            const ReactionGamma::Reaction *rctn = NULL;
            if( possible_rctns.size() )
              rctn = possible_rctns[0].reaction;
            
            if( rctn )
            {
              double nearestE = -999.9;
              for( const ReactionGamma::Reaction::EnergyYield &eip : rctn->gammas )
                if( fabs(eip.energy-energy) < fabs(nearestE-energy) )
                  nearestE = eip.energy;
              m_currentPeak.setReaction( rctn, static_cast<float>(nearestE), srcType );
            }else
              m_currentPeak.setReaction( NULL, 0.0f, PeakDef::NormalGamma );
          }//if( near ) / else (try for a reaction )
        }//if( nuctxt.find_first_of( "0123456789" ) == string::npos )
      }//if( nuc ) / else
    }//if( nuclideInfoIsDirty() )
    
    m_blockInfoRefresh = true;
    
    if( !(orig_continuum == *continuum) )
    {
      bool found_orig_peak = false;
      vector<PeakDef> new_peaks;
      for( const auto p : starting_peaks )
      {
        if( p->continuum() == orig_continuum_ptr )
        {
          PeakDef new_peak( *p );
          
          if( p == orig_peak_ptr )
          {
            found_orig_peak = true;
            new_peak = m_currentPeak;
          }
          
          new_peak.setContinuum( continuum );
          
          new_peaks.push_back( new_peak );
          m_peakModel->removePeak( p );
        }
      }//for( const auto p : starting_peaks )
      
      assert( found_orig_peak );
      assert( new_peaks.size() >= 1 );
        
      m_peakModel->addPeaks( new_peaks );
    }else
    {
      m_peakModel->removePeak( m_peakIndex );
      m_peakIndex = m_viewer->addPeak( m_currentPeak, false, SpecUtils::SpectrumType::Foreground );
      m_currentPeak = *m_peakModel->peak( m_peakIndex );
    }//if( orig_continuum != *m_currentPeak.continuum() )
    
    m_energy = m_currentPeak.mean();
    m_blockInfoRefresh = false;
    
    refreshPeakInfo();
    
    const double ending_energy = m_energy;
    vector<PeakModel::PeakShrdPtr> ending_peaks;
    const auto ending_peaks_ptr = m_peakModel->peaks();
    if( ending_peaks_ptr )
      ending_peaks.insert( end(ending_peaks), begin(*ending_peaks_ptr), end(*ending_peaks_ptr) );
    
    m_originalPeaks = ending_peaks;
    
  
    auto redo = [this, ending_peaks, ending_energy](){
      PeakEdit *editor = get_session_peak_editor();
      if( !editor )
        return;
      
      editor->m_peakModel->setPeaks( ending_peaks );
      editor->changePeak( ending_energy );
    };//redo
    
    auto undo = [starting_peaks, starting_energy](){
      PeakEdit *editor = get_session_peak_editor();
      if( !editor )
        return;
      
      editor->m_peakModel->setPeaks( starting_peaks );
      editor->changePeak( starting_energy );
    };//undo
    
    UndoRedoManager *undo_manager = m_viewer->undoRedoManager();
    if( undo_manager && undo_manager->canAddUndoRedoNow() )
      undo_manager->addUndoRedoStep( undo, redo, "Peak edit." );
  }catch ( std::exception &e)
  {
    if( m_peakIndex.isValid() )
      m_peakModel->removePeak( m_peakIndex );
    
    if( orig_peak_ptr )
    {
      m_peakIndex = m_viewer->addPeak( revertPeak, false , SpecUtils::SpectrumType::Foreground);
      m_currentPeak = *m_peakModel->peak( m_peakIndex );
      m_energy = m_currentPeak.mean();
    }//if( orig_peak_ptr )
    
    const string what = e.what();
    throw runtime_error( "Error applying changes to peak.<br/>Fix error(s) before retrying: " + what );
  }//try / catch
}//void apply()


void PeakEdit::accept()
{
  try
  {
    apply();
  }catch( std::exception &e )
  {
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
    return;
  }
  
  m_originalPeak = m_currentPeak;
  m_doneSignal.emit();
}//void accept()


void PeakEdit::cancel()
{
  if( m_peakIndex.isValid() )
  {
    shared_ptr<const deque<PeakModel::PeakShrdPtr>> current_peaks = m_peakModel->peaks();
    
    bool any_changes = false;
    if( !current_peaks )
    {
      any_changes = true;
    }else
    {
      const size_t num_orig_peaks = current_peaks ? current_peaks->size() : size_t(0);
      
      any_changes = (num_orig_peaks != m_originalPeaks.size());
      assert( any_changes || (num_orig_peaks == m_originalPeaks.size()) );
      for( size_t i = 0; !any_changes && i < std::min(num_orig_peaks, m_originalPeaks.size()); ++i )
      {
        auto orig = m_originalPeaks[i];
        auto current = (*current_peaks)[i];
        any_changes = !((*current) == (*orig));
      }
    }//if( !current_peaks ) / else
    
    if( any_changes )
    {
      m_peakModel->setPeaks( m_originalPeaks );
      PeakModel::PeakShrdPtr nearest = m_peakModel->nearestPeak( m_originalPeak.mean() );
      assert( nearest );
      
      if( nearest )
      {
        m_peakIndex = m_peakModel->indexOfPeak( nearest );
        m_currentPeak = *nearest;
      }else
      {
        m_peakIndex = WModelIndex();
        m_originalPeak.reset();
        m_currentPeak.reset();
      }
    }//if( any_changes )
  }//if( m_peakIndex.isValid() )
  
  m_doneSignal.emit();
}//void cancel()


void PeakEdit::deletePeak()
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  if( m_peakIndex.isValid() )
    m_peakModel->removePeak( m_peakIndex );
  m_peakIndex = WModelIndex();
  m_originalPeak.reset();
  m_originalPeaks.clear();
  m_currentPeak.reset();
  refreshPeakInfo();
  
  m_doneSignal.emit();
}//void deletePeak()



