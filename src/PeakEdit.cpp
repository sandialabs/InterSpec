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
#include <Wt/WComboBox>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WSuggestionPopup>
#include <Wt/WDoubleValidator>
#include <Wt/WContainerWidget>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakEdit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeNameFilterModel.h"


using namespace Wt;
using namespace std;
using SpecUtils::Measurement;

static_assert( PeakDef::Mean == 0,
               "PeakDef::CoefficientType::Mean typedef value has unexpected value" );
static_assert( (PeakDef::Chi2DOF+1) == int(PeakDef::NumCoefficientTypes),
               "PeakDef::CoefficientType::Chi2DOF or NumCoefficientTypes typedef values have unexpected values" );

static_assert( int(PeakEdit::Mean)            == int(PeakDef::Mean),
               "PeakEdit::Mean and PeakDef::Mean are not equal as expected" );
static_assert( int(PeakEdit::Sigma)           == int(PeakDef::Sigma),
               "PeakEdit::Sigma and PeakDef::Sigma are not equal as expected" );
static_assert( int(PeakEdit::GaussAmplitude)  == int(PeakDef::GaussAmplitude),
               "PeakEdit::GaussAmplitude and PeakDef::GaussAmplitude are not equal as expected" );
static_assert( int(PeakEdit::LandauAmplitude) == int(PeakDef::LandauAmplitude),
              "PeakEdit::LandauAmplitude and PeakDef::LandauAmplitude are not equal as expected" );
static_assert( int(PeakEdit::LandauMode)      == int(PeakDef::LandauMode),
               "PeakEdit::LandauMode and PeakDef::LandauMode are not equal as expected" );
static_assert( int(PeakEdit::LandauSigma)     == int(PeakDef::LandauSigma),
               "PeakEdit::LandauSigma and PeakDef::LandauSigma are not equal as expected" );
static_assert( int(PeakEdit::Chi2DOF)         == int(PeakDef::Chi2DOF),
               "PeakEdit::Chi2DOF and PeakDef::Chi2DOF are not equal as expected" );


PeakEditWindow::PeakEditWindow( const double energy,
                                PeakModel *peakmodel,
                                InterSpec *viewer )
  : AuxWindow( "Peak Editor", WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen) | AuxWindowProperties::DisableCollapse  )
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
    m_apply( NULL ),
    m_accept( NULL ),
    m_cancel( NULL ),
    m_refit( NULL ),
    m_otherPeaksDiv( NULL ),
    m_otherPeakTxt( NULL ),
    m_prevPeakInRoi( NULL ),
    m_nextPeakInRoi( NULL ),
    m_footer( aux->footer() ),
    m_aux( aux )
{
  init();
  changePeak( energy );
}//PeakEdit( constructor )


const char *PeakEdit::rowLabel( const PeakPars t )
{
  switch( t )
  {
    case PeakEdit::Mean:              return "Centroid";
    case PeakEdit::Sigma:             return "FWHM";
    case PeakEdit::GaussAmplitude:    return "Amplitude";
    case PeakEdit::LandauAmplitude:   return "Skew Amp.";
    case PeakEdit::LandauMode:        return "Skew Mode";
    case PeakEdit::LandauSigma:       return "Skew Width";
    case PeakEdit::OffsetPolynomial0: return "Cont. P0";
    case PeakEdit::OffsetPolynomial1: return "Cont. P1";
    case PeakEdit::OffsetPolynomial2: return "Cont. P2";
    case PeakEdit::OffsetPolynomial3: return "Cont. P3";
    case PeakEdit::OffsetPolynomial4: return "Cont. P4";
    case PeakEdit::RangeStartEnergy:  return "ROI Start (keV)";
    case PeakEdit::RangeEndEnergy:    return "ROI End (keV)";
    case PeakEdit::Chi2DOF:           return "&chi;2/DOF";
    case PeakEdit::PeakColor:         return "Peak Color";
    case PeakEdit::NumPeakPars:       return "";
      break;
  }//case( t )

  return "";
}//const char *PeakEdit::rowLabel( const PeakPars t )


void PeakEdit::init()
{
  addStyleClass( "PeakEdit" );
  
  WDoubleValidator *validator = new WDoubleValidator( this );
  
  m_valueTable = new WTable( this );
  m_valueTable->setHeaderCount( 1, Horizontal );
  m_valueTable->addStyleClass( "PeakEditTable" );
  
  WLabel *label = new WLabel( "Parameter", m_valueTable->elementAt(0,0) );
  label = new WLabel( "Value", m_valueTable->elementAt(0,1) );
  label = new WLabel( "Uncertainty", m_valueTable->elementAt(0,2) );
  label = new WLabel( "Fit", m_valueTable->elementAt(0,3) );
  
  for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
  {
    label = new WLabel( rowLabel(t), m_valueTable->elementAt(t+1,0) );
    m_values[t] = new WLineEdit( m_valueTable->elementAt(t+1,1) );
    m_uncertainties[t] = new WLineEdit( m_valueTable->elementAt(t+1,2) );
#if( BUILD_AS_OSX_APP || IOS )
    m_values[t]->setAttributeValue( "autocorrect", "off" );
    m_values[t]->setAttributeValue( "spellcheck", "off" );
    m_uncertainties[t]->setAttributeValue( "autocorrect", "off" );
    m_uncertainties[t]->setAttributeValue( "spellcheck", "off" );
#endif
    
    
    switch( t )
    {
      case PeakEdit::Mean: case PeakEdit::Sigma: case PeakEdit::GaussAmplitude:
      case PeakEdit::LandauAmplitude: case PeakEdit::LandauMode:
      case PeakEdit::LandauSigma:
      case PeakEdit::OffsetPolynomial0:
      case PeakEdit::OffsetPolynomial1: case PeakEdit::OffsetPolynomial2:
      case PeakEdit::OffsetPolynomial3: case PeakEdit::OffsetPolynomial4:
        m_fitFors[t] = new WCheckBox( m_valueTable->elementAt(t+1,3) );
      break;
      
      case PeakEdit::RangeStartEnergy: case PeakEdit::RangeEndEnergy:
      case PeakEdit::Chi2DOF: case PeakEdit::PeakColor:
      case PeakEdit::NumPeakPars:
        m_fitFors[t] = NULL;
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
      case PeakEdit::LandauAmplitude:   case PeakEdit::LandauMode:
      case PeakEdit::LandauSigma:
      case PeakEdit::OffsetPolynomial0: case PeakEdit::OffsetPolynomial1:
      case PeakEdit::OffsetPolynomial2: case PeakEdit::OffsetPolynomial3:
      case PeakEdit::OffsetPolynomial4:
 case PeakEdit::PeakColor:
      case PeakEdit::RangeStartEnergy:  case PeakEdit::RangeEndEnergy:
      case PeakEdit::Chi2DOF: case PeakEdit::NumPeakPars:
        break;
    }//case( t )

    if( !txt.empty() )
      m_uncertainties[t]->setToolTip( txt );
*/
    
    
    switch( t )
    {
      case PeakEdit::Mean:
//        m_values[t]->changed().connect( this, &PeakEdit::validateMeanOrRoiChange ); //doesnt seem to reliably work for some reason
        m_values[t]->blurred().connect( this, &PeakEdit::validateMeanOrRoiChange );
        m_values[t]->enterPressed().connect( this, &PeakEdit::validateMeanOrRoiChange );
      break;
        
      case PeakEdit::Sigma: case PeakEdit::GaussAmplitude:
      case PeakEdit::LandauAmplitude: case PeakEdit::LandauMode:
      case PeakEdit::LandauSigma:
        break;
      
      case PeakEdit::RangeStartEnergy: case PeakEdit::RangeEndEnergy:
//        m_values[t]->changed().connect( this, &PeakEdit::validateMeanOrRoiChange );
        m_values[t]->blurred().connect( this, &PeakEdit::validateMeanOrRoiChange );
        m_values[t]->enterPressed().connect( this, &PeakEdit::validateMeanOrRoiChange );
        //note, purposful fallthrough
        
      case PeakEdit::OffsetPolynomial0: case PeakEdit::OffsetPolynomial1:
      case PeakEdit::OffsetPolynomial2: case PeakEdit::OffsetPolynomial3:
      case PeakEdit::OffsetPolynomial4:
      case PeakEdit::Chi2DOF:
      case PeakEdit::PeakColor:
      case PeakEdit::NumPeakPars:
        m_uncertainties[t]->disable();
      break;
    }//case( t )
    
    
    m_values[t]->setValidator( validator );
    m_values[t]->addStyleClass( "numberValidator"); //used to detect mobile keyboard
    m_uncertainties[t]->setValidator( validator );
    m_uncertainties[t]->addStyleClass( "numberValidator"); //used to detect mobile keyboard
//    m_values[t]->changed().connect( boost::bind( &PeakEdit::checkIfDirty, this, t, false ) );
//    m_uncertainties[t]->changed().connect( boost::bind( &PeakEdit::checkIfDirty, this, t, true ) );    
    m_values[t]->blurred()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, false ) );
    m_values[t]->enterPressed()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, false ) );
    m_uncertainties[t]->blurred()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, true ) );
    m_uncertainties[t]->enterPressed()
            .connect( boost::bind( &PeakEdit::checkIfDirty, this, t, true ) );
    
    if( m_fitFors[t] )
    {
      m_fitFors[t]->checked().connect( boost::bind( &PeakEdit::fitTypeChanged, this, t ) );
      m_fitFors[t]->unChecked().connect( boost::bind( &PeakEdit::fitTypeChanged, this, t ) );
    }//if( m_fitFors[t] )
  }//for(...)

  WTableRow *row = NULL;
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+1 );
  label = new WLabel( "Nuclide", row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_nuclide = new WLineEdit( row->elementAt(1) );
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
  
  m_suggestions->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  m_suggestions->forEdit( m_nuclide,
                  WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );

  
//  std::shared_ptr<const PeakDef> dummypeak;
//  m_filterModel = new PeakIsotopeNameFilterModel( dummypeak, this );
  m_filterModel = new IsotopeNameFilterModel( this );
  m_filterModel->filter( "" );
  m_suggestions->setFilterLength( -1 );
  m_suggestions->setModel( m_filterModel );
  m_suggestions->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );
//  m_suggestions->filterModel().connect( m_filterModel, &PeakIsotopeNameFilterModel::filter );
  m_suggestions->filterModel().connect( m_filterModel, &IsotopeNameFilterModel::filter );
  
  m_nuclide->enterPressed().connect( this, &PeakEdit::isotopeChanged );
  m_nuclide->blurred().connect( this, &PeakEdit::isotopeChanged );
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+2 );
  label = new WLabel( "Photopeak", row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_photoPeakEnergy = new WComboBox( row->elementAt(1) );
  m_photoPeakEnergy->setWidth( WLength(13.5,WLength::FontEm) );
  m_photoPeakEnergy->activated().connect( this, &PeakEdit::transitionChanged );
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+3 );
  label = new WLabel( "Label", row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_userLabel = new WLineEdit( row->elementAt(1) );
#if( BUILD_AS_OSX_APP || IOS )
  m_userLabel->setAttributeValue( "autocorrect", "off" );
  m_userLabel->setAttributeValue( "spellcheck", "off" );
#endif
  
  m_userLabel->setWidth( WLength(100,WLength::Percentage) );
  m_userLabel->blurred()
              .connect( boost::bind( &PeakEdit::checkIfUserLabelDirty, this ) );
  m_userLabel->enterPressed()
              .connect( boost::bind( &PeakEdit::checkIfUserLabelDirty, this ) );
  
  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+4 );
  label = new WLabel( "Peak Color", row->elementAt(0) );
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
  label = new WLabel( "Peak Type", row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  
  m_peakType = new WComboBox( row->elementAt(1) );
  m_peakType->addItem( "Gaussian" );    //PeakDef::GaussianDefined
  m_peakType->addItem( "Data Region" ); //PeakDef::DataDefined
  m_peakType->activated().connect( this, &PeakEdit::peakTypeChanged );
//

  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+6 );
  label = new WLabel( "Continuum", row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_continuumType = new WComboBox( row->elementAt(1) );
  m_continuumType->activated().connect( this, &PeakEdit::contnuumTypeChanged );
  
  for( PeakContinuum::OffsetType t = PeakContinuum::OffsetType(0);
       t <= PeakContinuum::External; t = PeakContinuum::OffsetType(t+1) )
  {
    m_continuumType->addItem( PeakContinuum::offset_type_label(t) );
  }//for( loop over PeakContinuum::OffsetType )

  row = m_valueTable->rowAt( PeakEdit::NumPeakPars+7 );
  label = new WLabel( "Skew Type", row->elementAt(0) );
  row->elementAt(1)->setColumnSpan(2);
  m_skewType = new WComboBox( row->elementAt(1) );
  m_skewType->activated().connect( this, &PeakEdit::skewTypeChanged );
  
  for( PeakDef::SkewType t = PeakDef::SkewType(0);
      t <= PeakDef::LandauSkew; t = PeakDef::SkewType(t+1) )
  {
    switch ( t )
    {
      case PeakDef::NoSkew: m_skewType->addItem( "None" );            break;
      case PeakDef::LandauSkew: m_skewType->addItem( "Landau Skew" ); break;
    }//switch ( t )
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
  
  
  
  
//  WContainerWidget *buttonDiv = new WContainerWidget( m_footer );
  m_cancel = m_aux->addCloseButtonToFooter( "Cancel", false, m_footer );//new WPushButton( "Cancel", m_footer );
  
  m_refit  = new WPushButton( "Refit",  m_footer );
  m_apply  = new WPushButton( "Apply",  m_footer );

  m_accept = new WPushButton( "Accept", m_footer );
  if( m_viewer && !m_viewer->isMobile() )
    m_accept->setIcon( "InterSpec_resources/images/accept.png" );
  
  //Add class to give padding on left side (or modify current style class)
  
  WPushButton *deleteButton = new WPushButton( "Delete", m_footer );
//  deleteButton->setFloatSide( Wt::Right );
  
  m_cancel->clicked().connect( this, &PeakEdit::cancel );
  m_refit->clicked().connect(  this, &PeakEdit::refit  );
  m_apply->clicked().connect(  this, &PeakEdit::apply  );
  m_accept->clicked().connect( this, &PeakEdit::accept );
  deleteButton->clicked().connect( this, &PeakEdit::deletePeak );
  
  m_peakModel->dataChanged().connect( this, &PeakEdit::refreshPeakInfo );
  m_peakModel->layoutChanged().connect( this, &PeakEdit::refreshPeakInfo );
  m_peakModel->rowsRemoved().connect( this, &PeakEdit::peakModelRowsRemoved );
  m_peakModel->rowsInserted().connect( this, &PeakEdit::peakModelRowsAdded );
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
    findROIEnergyLimits( lowerx, upperx, *nearPeak,
                         m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground) );
  
  if( nearPeak && (energy<lowerx || energy>upperx) )
    nearPeak.reset();

//  if( !nearPeak )
//  {
//    m_peakIndex = WModelIndex();
//    m_currentPeak.reset();
//    m_originalPeak.reset();
//  }
  if( nearPeak )
  {
    m_energy = energy;
    m_peakIndex = nearestIndex;
    m_currentPeak = m_originalPeak = *nearPeak;
    for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
    {
      m_valIsDirty[t] = m_uncertIsDirty[t] = false;
    }
    
    refreshPeakInfo();
  }//if( !nearPeak ) / else
}//void changePeak( double energy )


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
    
    // Allowing new peak to be within 0.5*FWHM of original peak is an arbitrary choice, but should
    //  about uniquely identify the originally intended peak region
    if( nearest && (fabs(nearest->mean() - mean) < 0.5*nearest->fwhm()) )
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
      m_values[t]->setText( "" );
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
    m_otherPeakTxt->setText( "" );
    
    return;
  }//if( !nrow )
  
  m_currentPeak = *updated_peak;
  
  std::shared_ptr<const PeakContinuum> continuum = m_currentPeak.continuum();
  
  for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
  {
    WTableRow *row = m_valueTable->rowAt( t + 1 );
    
    const int rownum = static_cast<int>(t);
    double val = 0.0, uncert = 0.0;
    if( rownum < PeakDef::NumCoefficientTypes )
    {
      const PeakDef::CoefficientType ct = PeakDef::CoefficientType( rownum );
      val = m_currentPeak.coefficient( ct );
      uncert = m_currentPeak.uncertainty( ct );
    }
    
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
              size_t lower_channel, upper_channel;
              estimatePeakFitRange( m_currentPeak, data, lower_channel, upper_channel );
              val = data->gamma_channels_sum( lower_channel, upper_channel );
              
              // \TODO: The next line I think would be equivalent (or maybe more correct) than
              //        previous line, but I think edge-cases need checking
              //val = data->gamma_integral( m_currentPeak.lowerX(), m_currentPeak.upperX() );
            }else
            {
              val = 0.0;
            }//if( data ) / else
            break;
          }//case PeakDef::DataDefined:
        }//switch( m_currentPeak.m_type )
      break;
        
      case PeakEdit::LandauAmplitude:
      case PeakEdit::LandauMode:
      case PeakEdit::LandauSigma:
        row->setHidden( (m_currentPeak.skewType()==PeakDef::NoSkew) );
      break;
        
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
        
        if( continuum->energyRangeDefined() )
        {
          val = (t==RangeEndEnergy) ? continuum->upperEnergy()
                                    : continuum->lowerEnergy();
        }else
        {
          const std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);

          if( data )
          {
            size_t lower_channel, upper_channel;
            estimatePeakFitRange( m_currentPeak, data, lower_channel, upper_channel );
            if( t == PeakEdit::RangeEndEnergy )
              val += data->gamma_channel_upper( upper_channel );
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
        
      case PeakEdit::NumPeakPars:
      break;
    }//case( t )
    
    if( m_fitFors[t] )
    {
      if( rownum < PeakDef::NumCoefficientTypes )
      {
        const PeakDef::CoefficientType ct = PeakDef::CoefficientType( rownum );
        m_fitFors[t]->setChecked( m_currentPeak.fitFor(ct) );
      }else
      {
        const int coefnum = rownum - 2 - PeakDef::NumCoefficientTypes;
        const vector<bool> fitfor = continuum->fitForParameter();
        if( coefnum >= 0 && static_cast<size_t>(coefnum) < fitfor.size() )
          m_fitFors[t]->setChecked( fitfor[coefnum] );
      }//if( rownum < PeakDef::NumCoefficientTypes ) / else
    }//if( m_fitFors[t] )
    
    char valTxt[32], uncertTxt[32];
    snprintf( valTxt, sizeof(valTxt), "%.4f", val );
    snprintf( uncertTxt, sizeof(uncertTxt), "%.4f", uncert );
    
    if( row->isHidden() )
      uncertTxt[0] = valTxt[0] = '\0';
    
    m_valTxts[t] = valTxt;
    m_uncertTxts[t] = uncertTxt;
    m_values[t]->setText( valTxt );
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
  }else
  {
    char text[128];
    const size_t npeak = peaksInRoi.size();
    snprintf( text, sizeof(text), "Peak %i of %i in ROI", int(thispeak+1), int(npeak) );
    m_prevPeakInRoi->setHidden( thispeak == 0 );
    m_nextPeakInRoi->setHidden( thispeak == (npeak-1) );
    
    m_otherPeaksDiv->show();
    m_otherPeakTxt->setText( text );
  }//if( peaksInRoi.size() < 2 ) / else
  
  
  m_apply->disable();
//  m_accept->disable();
  m_cancel->enable();
  m_refit->setEnabled( m_currentPeak.type() == PeakDef::GaussianDefined );
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
    PeakDef::findNearestPhotopeak( nuclide, mean + extraEnergy, 4.0*sigma, xrayOnly,
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
    m_uncertIsDirty[t] = (m_uncertainties[t]->text()!=m_uncertTxts[t]);
    if( m_uncertIsDirty[t] )
    {
      m_apply->enable();
      m_accept->enable();
    }//if( m_uncertIsDirty[type] )
  }else
  {
    m_valIsDirty[t] = (m_values[t]->text()!=m_valTxts[t]);
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
    m_applyColorForAllNuc->setText( "All " + nuc->symbol + " peaks" );
    for( const auto &p : *allpeaks )
      nsamesrc += (p->parentNuclide()==nuc && p->lineColor()!=c);
  }else if( auto el = m_currentPeak.xrayElement() )
  {
    m_applyColorForAllNuc->setText( "All " + el->symbol + " peaks" );
    for( const auto & p : *allpeaks )
      nsamesrc += (p->xrayElement()==el && p->lineColor()!=c);
  }else if( auto rctn = m_currentPeak.reaction() )
  {
    m_applyColorForAllNuc->setText( "All " + rctn->name() + " peaks" );
    for( const auto & p : *allpeaks )
      nsamesrc += (p->reaction()==rctn && p->lineColor()!=c);
  }else
  {
    m_applyColorForAllNuc->setText( "All no src peaks" );
    for( const auto & p : *allpeaks )
      nsamesrc += (!p->parentNuclide() && !p->xrayElement() && !p->reaction()
                   && p->lineColor()!=c);
  }
  
  if( m_originalPeak.lineColor()!=c )
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
    case Mean: case Sigma: case GaussAmplitude: case LandauAmplitude:
    case LandauMode: case LandauSigma: case Chi2DOF:
    {
      PeakDef::CoefficientType coef = PeakDef::CoefficientType(static_cast<int>(t));
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
    
    case RangeStartEnergy: case RangeEndEnergy:
    case PeakColor: case NumPeakPars:
    break;
  };//switch( t )
  
  m_currentPeak = *m_peakModel->peak( m_peakIndex );
  
//  m_apply->enable();
  m_accept->enable();
}//void fitTypeChanged( PeakDef::CoefficientType t )


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
          
          case PeakEdit::LandauAmplitude:
          case PeakEdit::LandauMode:
          case PeakEdit::LandauSigma:
            row->setHidden( (m_currentPeak.skewType()==PeakDef::LandauSkew) );
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
            
          case PeakEdit::OffsetPolynomial4:
            row->setHidden( true );
          break;

          case PeakEdit::RangeStartEnergy:  row->setHidden( false ); break;
          case PeakEdit::RangeEndEnergy:    row->setHidden( false ); break;
          case PeakEdit::Chi2DOF:
            row->setHidden( m_currentPeak.chi2Defined() );
          break;
            
          case PeakEdit::PeakColor:         row->setHidden( false ); break;
          case PeakEdit::NumPeakPars: break;
        }//switch ( t )
      }//for( PeakDef::CoefficientType t )
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
          case PeakEdit::LandauAmplitude:   row->setHidden( true );  break;
          case PeakEdit::LandauMode:        row->setHidden( true );  break;
          case PeakEdit::LandauSigma:       row->setHidden( true );  break;
          case PeakEdit::OffsetPolynomial0:
          case PeakEdit::OffsetPolynomial1:
          case PeakEdit::OffsetPolynomial2:
          case PeakEdit::OffsetPolynomial3:
          case PeakEdit::OffsetPolynomial4:
//            row->setHidden( true );
          break;
          case PeakEdit::RangeStartEnergy:  row->setHidden( false ); break;
          case PeakEdit::RangeEndEnergy:    row->setHidden( false ); break;
          case PeakEdit::Chi2DOF:
            row->setHidden( m_currentPeak.chi2Defined() );
          break;
          case PeakEdit::PeakColor:         row->setHidden( false ); break;
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
      m_valueTable->rowAt(1+PeakEdit::LandauAmplitude)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::LandauMode)->setHidden( true );
      m_valueTable->rowAt(1+PeakEdit::LandauSigma)->setHidden( true );
    break;
      
    case PeakDef::LandauSkew:
      if( !m_currentPeak.gausPeak() )
      {
        passMessage( "Only Gaussian peaks can have a skew.", WarningWidget::WarningMsgHigh );
        m_skewType->setCurrentIndex( PeakDef::NoSkew );
        return;
      }//if( !m_currentPeak.gausPeak() )
      
      m_valueTable->rowAt(1+PeakEdit::LandauAmplitude)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::LandauMode)->setHidden( false );
      m_valueTable->rowAt(1+PeakEdit::LandauSigma)->setHidden( false );
      
      
      if( m_currentPeak.coefficient(PeakDef::LandauMode) <= 0.0 )
      {
        const double val = (2.5+0.3*0.22278)*m_currentPeak.sigma();
        
        char valTxt[32];
        snprintf( valTxt, sizeof(valTxt), "%.6f", val );
        m_currentPeak.set_coefficient( val, PeakDef::LandauMode );
        m_values[PeakDef::LandauMode]->setText( valTxt );
        m_uncertainties[PeakDef::LandauMode]->setText( "" );
      }//if( LandauMode not specficied )
      
      if( m_currentPeak.coefficient(PeakDef::LandauSigma) <= 0.0 )
      {
        const double val = 3.0*m_currentPeak.sigma();
        char valTxt[32];
        snprintf( valTxt, sizeof(valTxt), "%.6f", val );
        m_currentPeak.set_coefficient( val, PeakDef::LandauSigma );
        m_values[PeakDef::LandauSigma]->setText( valTxt );
        m_uncertainties[PeakDef::LandauSigma]->setText( "" );
      }//if( LandauMode not specficied )

      
      if( m_currentPeak.coefficient(PeakDef::LandauAmplitude) <= 0.0 )
      {
        //Integrate area between continuum and peak to guess multiple
        double amp = 0.003;
        std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
//        std::shared_ptr<const Measurement> continuum = m_viewer->continuum();
        
        if( data )
        {
          double lowe, highe;
          findROIEnergyLimits( lowe, highe, m_currentPeak, data );
          const double upe =  m_currentPeak.mean() - 3.5*m_currentPeak.sigma();
          if( upe > lowe )
          {
            double contArea = -1.0;
            const double dataArea = gamma_integral( data, lowe, upe );
            
            switch( continuum->type() )
            {
              case PeakContinuum::NoOffset:
                break;
                
              case PeakContinuum::External:
//                if( continuum )
//                  contArea = integral( continuum, lowe, upe );
                contArea = continuum->offset_integral( lowe, upe, nullptr );
                break;
                
              case PeakContinuum::Constant: case PeakContinuum::Linear:
              case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
//                contArea = m_currentPeak.offset_integral( lowe, upe );
                contArea = continuum->offset_integral( lowe, upe, nullptr );
                break;
                
              case PeakContinuum::FlatStep:
              case PeakContinuum::LinearStep:
              case PeakContinuum::BiLinearStep:
              {
                std::shared_ptr<const Measurement> data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
                contArea = continuum->offset_integral( lowe, upe, data );
              }
            }//switch( m_currentPeak.offsetType() )
            
            if( dataArea > contArea )
              amp = (dataArea-contArea) / m_currentPeak.peakArea();
            
            PeakDef tempPeak = m_currentPeak;
            tempPeak.setSkewType( type );
            tempPeak.set_coefficient( 1.0, PeakDef::LandauAmplitude );
            tempPeak.set_coefficient( 1.0, PeakDef::GaussAmplitude );
            const double frac = tempPeak.skew_integral( lowe, upe );
            
            if( (frac > 0.0) && !IsInf(frac) && !IsNan(frac) )
              amp /= frac;
            
            m_currentPeak.set_coefficient( amp, PeakDef::LandauAmplitude );
          }//if( upp > lowe )
        }//if( data )
        
        
        const string val = std::to_string( amp );
        m_values[PeakDef::LandauAmplitude]->setText( val );
        m_uncertainties[PeakDef::LandauAmplitude]->setText( "" );
      }//if( LandauAmplitude not specficied )
    break;
  }//switch( type )
  
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


void PeakEdit::refit()
{
  try
  {
    apply();
  }catch( std::exception &e )
  {
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
    return;
  }

  
  vector<PeakDef> inputPeak, fixedPeaks, outputPeak;
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
      fixedPeaks.push_back( *p );
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
  
  
  //if peak
  if( (inputPeak.size()>1) && thispeak->continuum()->isPolynomial() && !!data )
  {
    const std::shared_ptr<DetectorPeakResponse> &detector
                                = m_viewer->measurment(SpecUtils::SpectrumType::Foreground)->detector();
    const PeakShrdVec outp = refitPeaksThatShareROI( data, detector, inpkptrs, 0.25 );
    for( size_t i = 0; i < outp.size(); ++i )
      outputPeak.push_back( *outp[i] );
  }else
  {
    const bool isRefit = true;
    outputPeak = fitPeaksInRange( lowE, upE, ncausalitysigma, stat_threshold,
                                  hypothesis_threshold, inputPeak, data,
                                  fixedPeaks, isRefit );
  }
  
  
  if( outputPeak.size() != inputPeak.size() )
  {
    passMessage( "Failed to refit peak (became insignificant), so not doing "
                 "anything", WarningWidget::WarningMsgHigh );
    return;
  }//if( outputPeak.size() != inputPeak.size() )
  
  if( inputPeak.size() > 1 )
  {
    m_energy = m_currentPeak.mean();
    fixedPeaks.insert( fixedPeaks.end(), outputPeak.begin(), outputPeak.end() );
    std::sort( fixedPeaks.begin(), fixedPeaks.end(), &PeakDef::lessThanByMean );
    m_peakModel->setPeaks( fixedPeaks );
    
    changePeak( m_energy );
  }else
  {
    for( PeakDef::CoefficientType t = PeakDef::CoefficientType(0);
        t < PeakDef::NumCoefficientTypes;
        t = PeakDef::CoefficientType(t+1) )
    {
      outputPeak[0].setFitFor(t, m_currentPeak.fitFor(t));
    }//for( loop over PeakDef::CoefficientType )
    
    m_currentPeak = outputPeak[0];
    m_energy = m_currentPeak.mean();
    m_blockInfoRefresh = true;
    m_peakModel->removePeak( m_peakIndex );
    m_peakIndex   = m_viewer->addPeak( m_currentPeak, false );
    m_currentPeak = *m_peakModel->peak( m_peakIndex );
    m_blockInfoRefresh = false;
    
    refreshPeakInfo();
  }//if( inputPeak.size() > 1 )
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
  
  PeakDef revertPeak = *m_peakModel->peak( m_peakIndex ); //used to restore peak if catches an exception.
    
  try
  {
    std::shared_ptr<PeakContinuum> continuum = m_currentPeak.continuum();
    
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
      cerr << "PeakEdit::apply(): handle skew type change better" << endl;
      
      switch( skewType )
      {
        case PeakDef::NoSkew:
          m_currentPeak.setSkewType( skewType );
          break;
          
        case PeakDef::LandauSkew:
          m_currentPeak.setSkewType( skewType );
          break;
      }//switch( skewType )
    }//if( skewType != m_currentPeak.skewType() )
    
    
    for( PeakPars t = PeakPars(0); t < NumPeakPars; t = PeakPars(t+1) )
    {
      if( m_valIsDirty[t] )
      {
        if( m_values[t]->validate() != WValidator::Valid )
          throw runtime_error( "Value for '"
                              + string(rowLabel(t)) + "' is not a valid number" );
        
        double val = 0.0;
        string valtxt = m_values[t]->text().narrow();
        
        //      if( valtxt..size() && (sscanf( valtxt.c_str(), "%lf", &val ) != 1) )
        //        throw runtime_error( "Error converting " + valtxt + " to float" );
        if( !(stringstream(valtxt) >> val) && !valtxt.empty() )
          throw runtime_error( "Error converting " + valtxt + " to float" );
        
        switch( t )
        {
          case PeakEdit::GaussAmplitude:
            if( val < 0.0 )
              val = -val;
            m_currentPeak.setPeakArea( val );
            break;
            
          case PeakEdit::Sigma:
            val /= 2.3548201; //note: purposfull fall-through
            if( val < 0.0 )
              val = -val;
          
          case PeakEdit::Mean:
            m_energy = m_currentPeak.mean();
            //note: fall-through intentional
            
          case PeakEdit::LandauAmplitude:
          case PeakEdit::LandauMode:
          case PeakEdit::LandauSigma:
            m_currentPeak.set_coefficient( val, PeakDef::CoefficientType(static_cast<int>(t)) );
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
          case PeakEdit::NumPeakPars:
            break;
        }//case( t )
      }//if( m_valIsDirty[t] )
      
      if( m_uncertIsDirty[t] )
      {
        double val = 0.0;
        
        if( m_uncertainties[t]->validate() != WValidator::Valid )
          throw runtime_error( "Uncertainty for '"
                              + string(rowLabel(t)) + "' is not a valid number" );
        
        string valtxt = m_uncertainties[t]->text().narrow();
        if( !(stringstream(valtxt) >> val) && !valtxt.empty() )
          throw runtime_error( "Error converting " + valtxt + " to float" );
        
        if( val < 0.0 )
        {
          val = -val;
          //THis will be propogated to the GUI by refreshPeakInfo();
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
          case PeakEdit::LandauAmplitude:
          case PeakEdit::LandauMode:
          case PeakEdit::LandauSigma:
            m_currentPeak.set_uncertainty( val, PeakDef::CoefficientType(static_cast<int>(t)) );
            break;
            
          case PeakEdit::OffsetPolynomial0:
          case PeakEdit::OffsetPolynomial1:
          case PeakEdit::OffsetPolynomial2:
          case PeakEdit::OffsetPolynomial3:
          case PeakEdit::OffsetPolynomial4:
            continuum->setPolynomialUncert( t-OffsetPolynomial0, val );
            break;
            
          case PeakEdit::RangeStartEnergy: case PeakEdit::RangeEndEnergy:
          case PeakEdit::Chi2DOF:          case PeakEdit::PeakColor:
          case PeakEdit::NumPeakPars:
            break;
        }//case( t )
      }//if( m_uncertIsDirty[t] )
    }//for( PeakPars t )
    
    //We have to set peak type after setingin continuum parameters incase user
    //  changes a continuum type/polynomial-coefficent AND peak type
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
            double lowe, highe;
            findROIEnergyLimits( lowe, highe, m_currentPeak, data );
            continuum->setRange( lowe, highe );
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
        cerr << "Missread '" << currentTrans << "' as " << energy << endl;
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
        PeakDef::findNearestPhotopeak( nuc, energy, 0.0, xrayOnly,
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
              for( const ReactionGamma::EnergyAbundance &eip : rctn->gammas )
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
    m_peakModel->removePeak( m_peakIndex );
    m_peakIndex = m_viewer->addPeak( m_currentPeak, false );
    m_currentPeak = *m_peakModel->peak( m_peakIndex );
    m_blockInfoRefresh = false;
    
    refreshPeakInfo();
  }catch ( std::exception &e)
  {
    if( m_peakIndex.isValid() )
      m_peakModel->removePeak( m_peakIndex );
    m_peakIndex = m_viewer->addPeak( revertPeak, false );
    m_currentPeak = *m_peakModel->peak( m_peakIndex );
    m_energy = m_currentPeak.mean();
    
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
    m_peakModel->removePeak( m_peakIndex );
    m_peakIndex = m_viewer->addPeak( m_originalPeak, false );
    m_currentPeak = *m_peakModel->peak( m_peakIndex );
  }//if( m_peakIndex.isValid() )
  
  m_doneSignal.emit();
}//void cancel()


void PeakEdit::deletePeak()
{
  if( m_peakIndex.isValid() )
    m_peakModel->removePeak( m_peakIndex );
  m_peakIndex = WModelIndex();
  m_originalPeak.reset();
  m_currentPeak.reset();
  refreshPeakInfo();
  
  m_doneSignal.emit();
}//void deletePeak()



