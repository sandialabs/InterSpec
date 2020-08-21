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

#include <set>
#include <deque>
#include <limits>
#include <vector>
#include <sstream>
#include <iostream>
#include <iostream>

#include <Wt/WText>
#include <Wt/WTime>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WTextArea>
#include <Wt/WGroupBox>
#include <Wt/WTabWidget>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WItemDelegate>
#include <Wt/WDoubleSpinBox>
#include <Wt/WDoubleValidator>
#include <Wt/WCssDecorationStyle>


// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/Recalibrator.h"
#include "InterSpec/InterSpecApp.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpectraFileModel.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/IsotopeSelectionAids.h"


using namespace Wt;
using namespace std;

using SpecUtils::Measurement;
using SpecUtils::SpectrumType;

const int ForeIndex = static_cast<int>( SpecUtils::SpectrumType::Foreground );
const int BackIndex = static_cast<int>( SpecUtils::SpectrumType::Background );
const int SecondIndex = static_cast<int>( SpecUtils::SpectrumType::SecondForeground );

static_assert( ForeIndex   == 0, "SpecUtils::SpectrumType has changed and now code is invalid" );
static_assert( SecondIndex == 1, "SpecUtils::SpectrumType has changed and now code is invalid" );
static_assert( BackIndex   == 2, "SpecUtils::SpectrumType has changed and now code is invalid" );


namespace
{
  
  
  template<class T> struct index_compare_assend
  {
    index_compare_assend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
    bool operator()(const size_t a, const size_t b) const
    {
      return arr[a] < arr[b];
    }
    const T arr;
  };//struct index_compare

  
  
  
  
  
}//namespace

          

          

Recalibrator::Recalibrator( InterSpec *viewer, PeakModel *peakModel )
  : WContainerWidget( 0 ),
    m_interspec( viewer ),
    m_peakModel( peakModel ),
    m_isotopeView( 0 ),
    m_fitCoefButton( 0 ),
    m_multiFileButton( nullptr ),
    m_graphicalRecal( 0 ),
    m_applyToLabel( nullptr ),
    m_applyTo{ nullptr },
    m_convertToPolynomialLabel( nullptr ),
    m_convertToPolynomial( nullptr ),
    m_polyConverDiv( nullptr ),
    m_convertToPolyDialog( nullptr ),
    m_revert( 0 ),
    m_acceptText( 0 ),
    m_coeffEquationType( SpecUtils::EnergyCalType::InvalidEquationType ),
    m_devPairs( nullptr ),
    m_layout( nullptr )
{
  for( auto &v : m_uncertainties )
    v = -1.0;

  const Recalibrator::LayoutStyle style = m_interspec->isPhone() ? kTall : kWide;
  
  initWidgets( style );
  
  m_peakModel->rowsInserted().connect( this, &Recalibrator::checkIfCanFitCoefs );
  m_peakModel->rowsRemoved().connect( this, &Recalibrator::checkIfCanFitCoefs );
  
  m_peakModel->rowsInserted().connect( this, &Recalibrator::checkIfCanMultiFileCal );
  m_peakModel->rowsRemoved().connect( this, &Recalibrator::checkIfCanMultiFileCal );
  
  m_peakModel->dataChanged().connect( this, &Recalibrator::checkIfCanFitCoefs );
  
  
  // Wt::Signal<SpecUtils::SpectrumType,std::shared_ptr<SpecMeas>,std::set<int>>& m_interspec->displayedSpectrumChanged();
  
  
}//Recalibrator::Recalibrator(...)


Recalibrator::~Recalibrator()
{
  // no-op
} // Recalibrator::~Recalibrator()



void Recalibrator::setWasGraphicalRecal( int type, float energy )
{
  time( &m_lastGraphicalRecal );
  m_lastGraphicalRecalType = type;
  m_lastGraphicalRecalEnergy = energy;
}//void setWasGraphicalRecal( int type, double energy )

/*
 Passes in the parent when in Tall mode.  Necessary for help feature to know
 who the parent is
*/
void Recalibrator::initWidgets( Recalibrator::LayoutStyle style, AuxWindow *parent)
{
  m_currentLayout = style;
  
  const char *styleClasses[2] = {"RecalibratorTall", "RecalibratorWide"};
  const int styleIndToAdd = (style==Recalibrator::LayoutStyle::kTall ?  0 : 1);
  const int styleIndToRemove = (style==Recalibrator::LayoutStyle::kTall ?  1 : 0);
  if( !hasStyleClass(styleClasses[styleIndToAdd]) )
    addStyleClass( styleClasses[styleIndToAdd] );
  if( hasStyleClass(styleClasses[styleIndToRemove]) )
    removeStyleClass( styleClasses[styleIndToRemove] );
          
  const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
  if( !m_layout )
  {
    m_layout = new WGridLayout();
    setLayout( m_layout );
  }
          
  m_layout->clear();
  if( style == kTall )
    m_layout->setContentsMargins(9,9,9,9);
  else
    m_layout->setContentsMargins(5,3,4,1);
  
  
  // Set up the numeric boxes.
  WLabel *termlabels[sm_numCoefs] = { nullptr };
  for( size_t i = 0; i < sm_numCoefs; ++i )
  {
    m_coefficientDisplay[i] = new NativeFloatSpinBox();
    
    //m_coefficientDisplay[i]->setDecimals( 6 );

    m_coefficientDisplay[i]->setValue( 0.0 );
    //m_coefficientDisplay[i]->setRange( -1 * ( 2 << 20 ), 2 << 20 );

    // The spinboxes are disabled until a file is added~
    m_coefficientDisplay[i]->setDisabled( true );

    m_coefficientDisplay[i]->setMinimumSize( 95.0, WLength(1.0,WLength::FontEm) );
    
    // The buttons should all be enabled if the tick boxes are changed
    m_coefficientDisplay[i]->valueChanged().connect( boost::bind( &Recalibrator::engageRecalibration, this, ApplyRecal ) );
          
    termlabels[i] = new WLabel();
    if( i == 0 )
      termlabels[i]->setText( "Offset Term " );
    else if( i == 1 )
      termlabels[i]->setText( "Linear Term " );
    else if( i == 2 )
      termlabels[i]->setText( "Quad. Term " );
    else if( i == 3 )
      termlabels[i]->setText( "Cubic Term " );
    else if( i == 4 )
      termlabels[i]->setText( "Quart. Term " );
          
    m_termSuffix[i] = new WLabel();
          
    m_fitFor[i] = new WCheckBox( "fit" );
    m_fitFor[i]->setChecked( (i<2) );
    
    m_fitFor[i]->checked().connect( this, &Recalibrator::checkIfCanFitCoefs );
    m_fitFor[i]->unChecked().connect( this, &Recalibrator::checkIfCanFitCoefs );
    
    HelpSystem::attachToolTipOn( m_fitFor[i], "Fit for the value of this coefficient when using "
                                  "peaks associated with isotopes, to determine "
                                  "calibration." , showToolTipInstantly );
  }//for( size_t i = 0; i < sm_numCoefs; ++i )


  m_fitCoefButton = new WPushButton( "Fit Coeffs" );
  m_fitCoefButton->setIcon( "InterSpec_resources/images/ruler_small.png" );
  HelpSystem::attachToolTipOn( m_fitCoefButton,"Will use the expected energy of photopeaks "
                              "associated with the observed peak, in order to "
                              "fit calibration coefficients", showToolTipInstantly );
  //XXX - see notes above about postioning of m_isotopeView

  m_fitCoefButton->clicked().connect( boost::bind( &Recalibrator::recalibrateByPeaks, this ) );
  
  m_isotopeView = new RowStretchTreeView();
  m_isotopeView->setRootIsDecorated( false ); //makes the tree look like a table! :)
  
  m_isotopeView->addStyleClass( "RecalIsotopeView" );

  m_isotopeView->setModel( m_peakModel );
  const int numModelCol = m_peakModel->columnCount();
  for( int col = 0; col < numModelCol; ++col )
    m_isotopeView->setColumnHidden( col, true );

  m_isotopeView->setSortingEnabled( true );
  m_isotopeView->setAlternatingRowColors( true );
  m_isotopeView->setSelectable( true );
  m_isotopeView->setSelectionMode( SingleSelection );
  m_isotopeView->setEditTriggers( WAbstractItemView::SingleClicked
                                  | WAbstractItemView::DoubleClicked );

  m_isotopeView->setColumnHidden( PeakModel::kUseForCalibration, false );
  m_isotopeView->setColumnHidden( PeakModel::kMean, false );
  m_isotopeView->setColumnHidden( PeakModel::kIsotope, false );
  m_isotopeView->setColumnHidden( PeakModel::kPhotoPeakEnergy, false );
  m_isotopeView->setColumnHidden( PeakModel::kDifference, false );
  
  if( style == kTall )
  {
    m_isotopeView->setColumnWidth( PeakModel::kUseForCalibration, WLength(3.7, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kMean, WLength(3.3, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kIsotope, WLength(3.5, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(5.25, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kDifference, WLength(3.5, WLength::FontEm) );
  } //kTall
  else
  {
    m_isotopeView->setColumnWidth( PeakModel::kUseForCalibration, WLength(3.7, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kMean, WLength(4.5, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kIsotope, WLength(4.5, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(7.25, WLength::FontEm) );
    m_isotopeView->setColumnWidth( PeakModel::kDifference, WLength(6, WLength::FontEm) );
  } //kWide
  

  WItemDelegate *dblDelagate = new WItemDelegate( m_isotopeView );
  dblDelagate->setTextFormat( "%.2f" );
  m_isotopeView->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );

  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_isotopeView );
  m_isotopeView->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );

  PhotopeakDelegate *photopeakDelegate = new PhotopeakDelegate( PhotopeakDelegate::GammaEnergyDelegate, true, m_isotopeView );
  m_isotopeView->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );


  // Add the check boxes for what to apply to; initially all disabled
  m_applyToLabel = new WText( "Apply to:" );
  
  for( int i = 0; i < 3; ++i )
  {
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
    const char *name = 0;
    switch( type )
    {
      case SpectrumType::Foreground:       name = "Foreground"; break;
      case SpectrumType::Background:       name = "Background"; break;
      case SpectrumType::SecondForeground: name = "2nd Spec"; break;
    }
    
    m_applyTo[i] = new WCheckBox( name );
    m_applyTo[i]->setChecked( true );
    m_applyTo[i]->setDisabled(true);
    m_applyTo[i]->checked().connect( boost::bind(&Recalibrator::specTypeCheckedCallback, this, type, true ) );
    m_applyTo[i]->unChecked().connect( boost::bind(&Recalibrator::specTypeCheckedCallback, this, type, false ) );
  }//for( loop over SpecUtils::SpectrumType )
  
  
  
  
 WContainerWidget *footer = nullptr;
  if( style == kTall )
    footer = (parent ? parent->footer() : new WContainerWidget());
          
  if( footer )
    AuxWindow::addHelpInFooter( footer, "energy-calibration-dialog" );
  
  m_multiFileButton = new WPushButton( "Mult. Files...", footer );
  
  if( style == kTall )
  {
    WPushButton *closeBtn = nullptr;
    if( parent )
      closeBtn = parent->addCloseButtonToFooter( "Close", false, footer );
    else
      closeBtn = new WPushButton( "Close", footer );
    closeBtn->clicked().connect( m_interspec, &InterSpec::handEnergyCalWindowClose );
    closeBtn->addStyleClass( "RecalCancelbtn" );
  }//if( style == kTall )
  
  m_revert = new WPushButton( "Revert", footer );
  if( style == kTall )
    m_revert->addStyleClass( "RecalCancelbtn" );
  m_revert->setIcon( "InterSpec_resources/images/arrow_undo.png" );

  // And disable all the buttons
  m_revert->disable();

  //Add a Polynomial/FullWidthFraction/LowerEdge select here
  //Or add a button for when LowerEdge that converts to Polynomial
  
  //Hook the buttons up
  m_revert->clicked().connect( boost::bind( &Recalibrator::engageRecalibration, this, RevertRecal ) );
  
  int row = 0; //keeps track of row for adding widgets to layout
  for( size_t i = 0; i < sm_numCoefs; ++i )
  {
    m_layout->addWidget( termlabels[i],           row, 0, Wt::AlignMiddle );
    m_layout->addWidget( m_coefficientDisplay[i], row, 1, Wt::AlignMiddle );
    m_layout->addWidget( m_termSuffix[i],         row, 2, Wt::AlignMiddle );
    m_layout->addWidget( m_fitFor[i],             row, 3, Wt::AlignCenter | Wt::AlignMiddle );
    termlabels[i]->setBuddy( m_coefficientDisplay[i] );
    m_termSuffix[i]->setBuddy( m_coefficientDisplay[i] );
    row++;
  }
  
  m_devPairs = new DeviationPairDisplay( nullptr );
  //m_devPairs->changed().connect( boost::bind( &Recalibrator::engageRecalibration, this, ApplyRecal ) );
  m_devPairs->changed().connect( boost::bind( &Recalibrator::userChangedDeviationPairs, this ) );
  
  
  if( style == kTall )
  {
    m_devPairs->addStyleClass( "DevPairDisplayTall" );
          
    m_layout->addWidget( m_devPairs, row, 0, 1, 4 );
    m_devPairs->setMinimumSize(WLength::Auto, WLength(150));
    m_devPairs->setMaximumSize(WLength::Auto, WLength(150));
    row++;

    m_layout->addWidget( m_isotopeView, row, 0, 1, 4 );
    m_isotopeView->setMinimumSize(WLength::Auto, WLength(150));
    m_layout->setRowStretch( row, 1 );
    m_layout->setColumnStretch( 1, 1 );
    row++;

    m_layout->addWidget( m_fitCoefButton, row, 3, 1, 1 );
    m_layout->addWidget( m_multiFileButton, row, 2, 1, 1 );
    row++;
  }//kTall
  
  m_applyTo[ForeIndex]->setStyleClass( "RecalApplyToCbL" );
  m_applyTo[SecondIndex]->setStyleClass( "RecalApplyToCbR" );
  m_applyTo[BackIndex]->setStyleClass( "RecalApplyToCbL" );

  m_layout->addWidget( m_applyToLabel,         row, 0 );
  m_layout->addWidget( m_applyTo[ForeIndex],   row, 1 );
  m_layout->addWidget( m_applyTo[SecondIndex], row, 2 );
  row++;
          
  m_layout->addWidget( m_applyTo[BackIndex],   row, 1 );
  row++;
  
  if( style == kTall )
  {
    m_polyConverDiv = new WContainerWidget();
    m_layout->addWidget( m_polyConverDiv, row, 0, 1, 3 );
    row++;
          
    m_convertToPolynomialLabel = new WText( "Calibration is specified by channel energy" );
    m_convertToPolynomialLabel->addStyleClass( "RecalConvertToPolyTxt" );
                        
    m_convertToPolynomial = new WPushButton( "Convert To Polynomial..." );
    m_convertToPolynomial->clicked().connect( this, &Recalibrator::startConvertToPolynomial );
          
    //m_layout->addWidget( m_convertToPolynomialLabel, row, 0, 1, 3 );
    //row++;
    //m_layout->addWidget( m_convertToPolynomial, row, 0, 1, 3, Wt::AlignCenter | Wt::AlignTop );
    //row++;
          
    m_convertToPolynomialLabel->setInline( false );
    m_polyConverDiv->addWidget( m_convertToPolynomialLabel );
    m_polyConverDiv->addWidget( m_convertToPolynomial );
          
    m_convertToPolynomialLabel->hide();
    m_convertToPolynomial->hide();
    m_polyConverDiv->hide();
  }//if( style == kTall )
  
  m_multiFileButton->setIcon( "InterSpec_resources/images/page_white_stack.png" );
  HelpSystem::attachToolTipOn( m_multiFileButton, "Tool to use peaks from multiple files to fit for calibration",
                              showToolTipInstantly, HelpSystem::Top );
  m_multiFileButton->clicked().connect( this, &Recalibrator::createMultifileFitter );
  
  
  if( style == kWide )
  {
    m_devPairs->addStyleClass( "DevPairDisplayWide" );
    
    //for( int i = 0; i < (row-1); ++i )
    //  m_layout->setRowStretch( i, 1 );
    
    //m_layout->addWidget( m_revert,        row, 0, 1,   2, Wt::AlignLeft | Wt::AlignMiddle );
    //row++;
            
    //m_layout->addWidget( new WContainerWidget(), row, 0 );
    m_layout->addWidget( m_fitCoefButton,  row, 2, 1, 2, Wt::AlignCenter | Wt::AlignBottom );
    
    auto helpBtn = new WContainerWidget();
    helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
    helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "energy-calibration-dialog" ) );
    m_layout->addWidget( helpBtn,  row, 0, Wt::AlignLeft | Wt::AlignBottom );
    row++;
    
    m_layout->addWidget( m_isotopeView,     0,   4, row, 1 );
    m_layout->addWidget( m_devPairs,        0,   5, row, 1 );
          
    m_layout->addWidget( m_revert,          0, 6 );
    m_layout->addWidget( m_multiFileButton, 1, 6 );
    
    m_polyConverDiv = new WContainerWidget();
    m_convertToPolynomialLabel = new WText( "Cal. is by channel energy", m_polyConverDiv );
    m_convertToPolynomialLabel->setInline( false );
    m_convertToPolynomialLabel->decorationStyle().font().setSize( WLength(60,WLength::Percentage) );
          
    m_convertToPolynomial = new WPushButton( "To Polynomial...", m_polyConverDiv );
    m_convertToPolynomial->clicked().connect( this, &Recalibrator::startConvertToPolynomial );
    
    //m_layout->addWidget( m_convertToPolynomialLabel, 2, 6, AlignBottom | AlignCenter );
    //m_layout->addWidget( m_convertToPolynomial,      3, 6, AlignTop );
    m_layout->addWidget( m_polyConverDiv, 2, 6, 2, 1 );
    
    m_convertToPolynomialLabel->hide();
    m_convertToPolynomial->hide();
    m_polyConverDiv->hide();
          
    m_layout->setColumnStretch( 4, 1 );
    m_layout->setRowStretch( row-1, 1 );
  }//if( style == kWide )
  
  if( style == kTall )
  {
    removeStyleClass( "RecalibratorWide" );
    addStyleClass( "RecalibratorTall" );
  }else
  {
    removeStyleClass( "RecalibratorTall" );
    addStyleClass( "RecalibratorWide" );
  }//if( style == kTall ) / else
}//void Recalibrator::initWidgets(...)


Recalibrator::LayoutStyle Recalibrator::currentLayoutStyle() const
{
  return m_currentLayout;
}


void Recalibrator::setWideLayout()
{
  initWidgets( kWide );
  refreshRecalibrator();
}//void Recalibrator::setWideLayout()


void Recalibrator::setTallLayout( AuxWindow *parent )
{
  initWidgets( kTall, parent );
  refreshRecalibrator();
}//void Recalibrator::setTallLayout()


void Recalibrator::deleteGraphicalRecalConfirmWindow()
{
  if( m_graphicalRecal )
  {
    delete m_graphicalRecal;
    m_graphicalRecal = NULL;
  }//if( m_graphicalRecal )
  
  const bool showToolTipInstantly
      = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  if( showToolTipInstantly )
  {
#if( USE_SPECTRUM_CHART_D3 )
    passMessage( "If you recalibrate again by ALT+CTRL+DRAG on"
                " another portion of the spectrum, you will be given the"
                " option of preserving the effects of this calibration"
                " as well.", "", WarningWidget::WarningMsgInfo );
#else
    passMessage( "If you recalibrate again by right-click dragging on"
                " another portion of the spectrum, you will be given the"
                " option of preserving the effects of this calibration"
                " as well.", "", WarningWidget::WarningMsgInfo );
#endif
  }
  
}//deleteGraphicalRecalConfirmWindow()



void Recalibrator::specTypeCheckedCallback( const SpecUtils::SpectrumType type,
                                            const bool checked )
{
  //This function makes sure that the "Apply To" checkboxes are consistent for
  //  the foreground/background/second if they share the same file (e.g., if
  //  the background and foreground are from the same file, then the
  //  recalibration will always be applied to both if either are selected)
  
  const int typeindex = static_cast<int>(type);
  std::shared_ptr<SpecMeas> meas = m_interspec->measurment(type);
  
  for( int i = 0; i < 3; ++i )
  {
    const SpecUtils::SpectrumType t = SpecUtils::SpectrumType(i);
    if( (t != type) && (meas == m_interspec->measurment(t)) )
      m_applyTo[i]->setChecked( m_applyTo[typeindex]->isChecked() );
  }//
}//void specTypeCheckedCallback(...)


void Recalibrator::shiftPeaksForEnergyCalibration( PeakModel *peakmodel,
                                                   const std::vector<float> &new_pars,
                                                   const std::vector< std::pair<float,float> > &new_devpairs,
                                                   SpecUtils::EnergyCalType new_eqn_type,
                                                   std::shared_ptr<SpecMeas> meas,
                                                   const SpecUtils::SpectrumType spectype,
                                                   std::vector<float> old_pars,
                                                   const std::vector< std::pair<float,float> > &old_devpairs,
                                                   SpecUtils::EnergyCalType old_eqn_type )
{
  
  if( old_eqn_type==SpecUtils::EnergyCalType::LowerChannelEdge
     || old_eqn_type==SpecUtils::EnergyCalType::InvalidEquationType
     || new_eqn_type==SpecUtils::EnergyCalType::LowerChannelEdge
     || new_eqn_type==SpecUtils::EnergyCalType::InvalidEquationType )
  {
    //ToDo: We can only currently handle Polynomial or Full Range Fraction
    //      calibrations.  This could be fixed.
    return;
  }
    
  if( !meas )
    throw runtime_error( "shiftPeaksForEnergyCalibration: invalid input" );
  
  if( spectype != SpectrumType::Foreground )
  {
    assert( 0 );
    //meas->shiftPeaksForRecalibration( old_pars, old_devpairs, old_eqn_type,
    //                                  new_pars, new_devpairs, new_eqn_type );
    return;
  }
  
  if( !peakmodel )
    throw runtime_error( "shiftPeaksForEnergyCalibration: no peak model passed in" );
  
  const set<int> samples = meas->displayedSampleNumbers();
  std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > > peaks;
  peaks = meas->peaks( samples );
  if( peaks != peakmodel->peaks() )
    throw logic_error( "Recalibrator::shiftPeaksForEnergyCalibration: logic error in retireving peaks" );
  
  peakmodel->setPeakFromSpecMeas( std::shared_ptr<SpecMeas>(), set<int>() );
  
  assert( 0 );
  //meas->shiftPeaksForRecalibration( old_pars, old_devpairs, old_eqn_type,
  //                                 new_pars, new_devpairs, new_eqn_type );
  
  
  peakmodel->setPeakFromSpecMeas( meas, samples );
}//shiftPeaksForEnergyCalibration(...)




void Recalibrator::handleGraphicalRecalRequest( double xstart, double xfinish )
{
  try
  {
    if( !m_interspec->displayedHistogram(SpectrumType::Foreground) )
      return;
  
    if( m_graphicalRecal )
      m_graphicalRecal->setEnergies( xstart, xfinish );
    else
      m_graphicalRecal = new GraphicalRecalConfirm( xstart, xfinish, this,
                                                 m_lastGraphicalRecal,
                    GraphicalRecalConfirm::RecalTypes(m_lastGraphicalRecalType),
                                                 m_lastGraphicalRecalEnergy );
  }catch( std::runtime_error & )
  {
    passMessage( "Internal error doing graphical recal; sorry :(", "",
                 WarningWidget::WarningMsgHigh );
  }
}//void handleGraphicalRecalRequest( double xstart, double xfinish )


void Recalibrator::createMultifileFitter()
{
  new MultiFileCalibFit( this );
}//void createMultifileFitter()


void Recalibrator::startConvertToPolynomial()
{
  if( m_convertToPolyDialog )
  {
    passMessage( "Convert dialog is already showing!", "", WarningWidget::WarningMsgHigh );
    m_convertToPolyDialog->show();
    return;
  }
  
  auto foreground = m_interspec->measurment(SpectrumType::Foreground);
  auto displ_foreground = m_interspec->displayedHistogram(SpectrumType::Foreground);
  
  if( !foreground || !displ_foreground )
  {
    passMessage( "No foreground to convert to polynomial binning.", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  if( displ_foreground->energy_calibration_model() != SpecUtils::EnergyCalType::LowerChannelEdge )
  {
    passMessage( "Foreground calibration is not specified by lower channel - wont convert to polynomial.", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  const vector<string> used_dets = m_interspec->detectorsToDisplay(SpectrumType::Foreground);
  const set<int> samples = m_interspec->displayedSamples(SpectrumType::Foreground);
  auto energy_cal = foreground->suggested_sum_energy_calibration( samples, used_dets );
  
  if( !energy_cal )
  {
    passMessage( "Issue getting channel energies of foreground.", "", WarningWidget::WarningMsgHigh );
    return;
  }

  
  m_convertToPolyDialog = new AuxWindow( "Confirm",
                                        WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneModal)
                                                | DisableCollapse | SetCloseable | IsAlwaysModal );
  m_convertToPolyDialog->rejectWhenEscapePressed();
  
  m_convertToPolyDialog->setWidth( 300 );
  WContainerWidget *body = m_convertToPolyDialog->contents();
  
  new WText( "Converting from energy calibration specified by each channels"
             " lower energy to polynomial will lose any non-linearity in the"
             " energy calibration.<br />"
             "This same calibration will be applied to all spectra in the "
             " foreground, background, and secondary files.<br />"
             "<br />"
             "Are you sure you would like to proceed?",
             TextFormat::XHTMLText, body );
  
  Wt::WPushButton *cancel = m_convertToPolyDialog->addCloseButtonToFooter( "No" );
  cancel->clicked().connect( boost::bind( &Recalibrator::finishConvertToPolynomial, this, false ) );
  
  Wt::WPushButton *proceed = m_convertToPolyDialog->addCloseButtonToFooter( "Yes" );
  proceed->clicked().connect( boost::bind( &Recalibrator::finishConvertToPolynomial, this, true ) );
  
  m_convertToPolyDialog->finished().connect( boost::bind( &Recalibrator::finishConvertToPolynomial, this, false ) );
  
  m_convertToPolyDialog->centerWindowHeavyHanded();
}//void startConvertToPolynomial()


void Recalibrator::finishConvertToPolynomial( const bool followThrough )
{
  if( !m_convertToPolyDialog )
    return;
  
  AuxWindow::deleteAuxWindow( m_convertToPolyDialog );
  m_convertToPolyDialog = nullptr;
  
  if( !followThrough )
    return;
  
  std::shared_ptr<SpecMeas> foreground = m_interspec->measurment(SpectrumType::Foreground);
  std::shared_ptr<SpecMeas> background = m_interspec->measurment(SpectrumType::Background);
  std::shared_ptr<SpecMeas> secondary = m_interspec->measurment(SpectrumType::SecondForeground);
  
  auto displ_foreground = m_interspec->displayedHistogram(SpectrumType::Foreground);
  
  if( !foreground || !displ_foreground
      || displ_foreground->num_gamma_channels() < 16
      || !displ_foreground->channel_energies()
      || displ_foreground->channel_energies()->empty() )
  {
    passMessage( "Recalibrator::finishConvertToPolynomial: invalid foreground.", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  //ToDo: If any of of the SpecMeas objects have Measurements with a differing
  //      number of channels, we will get an exception.  Should handle this
  //      correctly...
  map<pair<size_t,float>,shared_ptr<SpecUtils::EnergyCalibration>> newcals;
  
  auto set_cal = [&newcals]( std::shared_ptr<SpecMeas> &meas ){
    if( !meas )
      return;
    for( auto m : meas->measurements() )
    {
      const size_t nchannel = m->num_gamma_channels();
      const float upper_energy = m->gamma_energy_max();
      if( nchannel < 2 || upper_energy < 1.0f )
        continue;
      
      shared_ptr<SpecUtils::EnergyCalibration> &cal = newcals[{nchannel,upper_energy}];
      if( !cal )
      {
        cal = make_shared<SpecUtils::EnergyCalibration>();
        cal->set_polynomial( nchannel, {0.0f,upper_energy/nchannel}, {} );
      }
      
      if( cal != m->energy_calibration() )
        meas->set_energy_calibration( cal, m );
    }//for( auto m : foreground->measurements() )
  };//set_cal lamda
  
  set_cal( foreground );
  set_cal( background );
  set_cal( secondary );
  
  m_interspec->refreshDisplayedCharts();
  
  refreshRecalibrator();
}//void finishConvertToPolynomial()


void Recalibrator::checkIfCanFitCoefs()
{
  //An incredibly minor optimization: we will only call enable/disable if its different than the current state
  //  (from a superficial glance at Wt source code I think it effectively does this check anyway to avoid re-rendering
  //   but I guess I'll just leave this check in for the moment).
  const bool wasEnabled = m_fitCoefButton->isEnabled();
  
  std::shared_ptr<const SpecUtils::Measurement> disp_foreground
             = m_interspec->displayedHistogram(SpectrumType::Foreground);
  
  if( !disp_foreground || !m_applyTo[ForeIndex]->isChecked() )
  {
    if( wasEnabled )
      m_fitCoefButton->disable();
    return;
  }
  
  const SpecUtils::EnergyCalType cal_type = disp_foreground->energy_calibration_model();
  
  switch( cal_type )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    break;
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      if( wasEnabled )
        m_fitCoefButton->disable();
      return;
    break;
  }//switch( cal_type )
  
  int num_coeff_fit = 0;
  for( size_t i = 0; i < sm_numCoefs; ++i )
    num_coeff_fit += m_fitFor[i]->isChecked();
  
  if( !num_coeff_fit )
  {
    if( wasEnabled )
      m_fitCoefButton->disable();
    return;
  }
  
  int npeaks_use = 0;
  const size_t npeaks = m_peakModel->npeaks();
  for( size_t peakn = 0; peakn < npeaks; ++peakn )
  {
    const PeakDef &peak = m_peakModel->peak(peakn);
    npeaks_use += peak.useForCalibration();
  }//for( loop over peakn )
  
  if( !npeaks_use )
  {
    if( wasEnabled )
      m_fitCoefButton->disable();
    return;
  }
  
  if( !wasEnabled )
    m_fitCoefButton->enable();
}//void checkIfCanFitCoefs()


void Recalibrator::checkIfCanMultiFileCal()
{
  //See notes in #Recalibrator::checkIfCanFitCoefs about this minor optimization
  const bool wasEnabled = m_multiFileButton->isEnabled();
  
  SpectraFileModel *fileModel = m_interspec->fileManager()->model();
  if( !fileModel )
  {
    if( wasEnabled )
      m_multiFileButton->disable();
    return;
  }//if( !fileModel )
  
  const int nfile = fileModel->rowCount();
  if( nfile > 1 )
  {
    if( !wasEnabled )
      m_multiFileButton->enable();
    return;
  }
  
  auto meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  if( !meas )
  {
    if( wasEnabled )
      m_multiFileButton->disable();
    return;
  }//if( !fileModel )
  
  const set<set<int> > peakSamples = meas->sampleNumsWithPeaks();
  if( peakSamples.size() < 2 )
  {
    if( wasEnabled )
      m_multiFileButton->disable();
    return;
  }
  
  if( !wasEnabled )
    m_multiFileButton->enable();
}//void checkIfCanMultiFileCal()

          
void Recalibrator::userChangedDeviationPairs()
{
  try  //begin check validity of equation {note: not a rigorous check}
  {
    //Lets be optimistic and and assume all will be well
    m_devPairs->setValidValues();
          
    vector<pair<float,float>> devpairs = m_devPairs->deviationPairs();
    if( devpairs.empty() )
    {
      engageRecalibration(ApplyRecal);
      return;
    }
          
    std::shared_ptr<const SpecUtils::Measurement> disp_foreground
               = m_interspec->displayedHistogram( SpectrumType::Foreground );
    if( !disp_foreground )
    {
      m_devPairs->setDeviationPairs( vector<pair<float,float>>{} );
      throw std::runtime_error( "You must have a foreground loaded to set deviation pairs." );
    }
        
    auto &channel_counts = disp_foreground->gamma_counts();
    if( !channel_counts || channel_counts->empty() )
    {
      m_devPairs->setDeviationPairs( disp_foreground->deviation_pairs() );
      throw std::runtime_error( "Foreground does not have a defined spectrum." );
    }
    
    vector<float> frfcoef;
    const vector<float> eqn = disp_foreground->calibration_coeffs();
    const size_t nchannel = channel_counts->size();
    const SpecUtils::EnergyCalType cal_type = disp_foreground->energy_calibration_model();
          
    switch( cal_type )
    {
      case SpecUtils::EnergyCalType::InvalidEquationType:
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        m_devPairs->setDeviationPairs( disp_foreground->deviation_pairs() );
        throw std::runtime_error( "Deviation pairs can only be defined for polynomial or full range fraction calibrations." );
        break;
          
      case SpecUtils::EnergyCalType::FullRangeFraction:
          frfcoef = eqn;
        break;
          
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        frfcoef = SpecUtils::polynomial_coef_to_fullrangefraction( eqn, nchannel );
        break;
    }//switch( calibration type )
          
    
    const bool valid = checkFullRangeFractionCoeffsValid( frfcoef, devpairs, nchannel );
    
    if( !valid )
    {
      m_devPairs->setInvalidValues();
      m_revert->setDisabled( false );
      stringstream msg;
      msg << "The coeffiecients {" << eqn[0] << ", " << eqn[1] << ", " << eqn[2]
          << "} are invalid because they will cause higher numbered bins to"
          << " have lower energy values.";

      throw std::runtime_error( msg.str() );
    }//if( (nearend >= end) || (begin >= nearbegin) )
          
    engageRecalibration(ApplyRecal);
  }catch( std::exception & e )//end codeblock
  {
    passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    
    return;
  }//try / catch to check validity of equation
}//void userChangedDeviationPairs()


void Recalibrator::engageRecalibration( RecalAction action )
{
  cout << "engageRecalibration" << endl;
  std::shared_ptr<SpecMeas> foreground, back, second;
  foreground = m_interspec->measurment(SpectrumType::Foreground);
  back       = m_interspec->measurment(SpectrumType::Background);
  second     = m_interspec->measurment(SpectrumType::SecondForeground);

  std::shared_ptr<const SpecUtils::Measurement> disp_foreground, disp_back, disp_second;
  disp_foreground = m_interspec->displayedHistogram( SpectrumType::Foreground );
  disp_back = m_interspec->displayedHistogram( SpectrumType::Background );
  disp_second = m_interspec->displayedHistogram( SpectrumType::SecondForeground );
  
  if( !foreground || !foreground->num_gamma_channels()
     || !disp_foreground || !disp_foreground->num_gamma_channels() )
    return;
 
#if( PERFORM_DEVELOPER_CHECKS )
  const double origForSum    = ((!!foreground) ? foreground->deep_gamma_count_sum() : 0.0);
  const double origBackSum   = ((!!back) ? back->deep_gamma_count_sum() : 0.0);
  const double origSecondSum = ((!!second) ? second->deep_gamma_count_sum() : 0.0);
#endif

  
  const std::vector<std::pair<float,float>> devpairs = m_devPairs->deviationPairs();
  const std::vector<std::pair<float,float>> &olddevpairs = disp_foreground->deviation_pairs();
  
  const vector< float > originalPars = disp_foreground->calibration_coeffs();
  const SpecUtils::EnergyCalType old_calib_type = disp_foreground->energy_calibration_model();
  
  if( action == ApplyRecal )
  {
    m_revert->setDisabled( false );
    
    for( auto w : m_applyTo )
      w->setDisabled( true );
  }else
  {
    m_revert->setDisabled( true );
  }

  // If they just want to RevertRecal, pass the original data back in, and then just refresh.
  if( action == RevertRecal )
  {
    passMessage( "Calibration reverted.", "", 0 );

    if( !!second
        && m_applyTo[SecondIndex]->isChecked()
        && second!=foreground && second!=back )
    {
      const CalibrationInformation &info = m_originalCal[SecondIndex];
      
      try
      {
        for( const auto &iter : info.m_cal_to_meas )
        {
          const shared_ptr<const SpecUtils::EnergyCalibration> &cal = iter.first;
          for( const shared_ptr<const SpecUtils::Measurement> &m : iter.second )
            second->set_energy_calibration( cal, m );
        }
        
        //blurg - when we shift the peaks, we have to know from which calibration to which calibration,
        //  which means we need to consistently define whcih energy calbiration defines the display,
        //  or at least keep track of it for when we change things.
        assert( 0 );
        //second->shiftPeaksForRecalibration( disp_second->calibration_coeffs(),
        //disp_second->deviation_pairs(),
        //disp_second->energy_calibration_model(),
        //info.coefficients,
        //info.deviationpairs,
        //info.type );
      }catch( std::exception &e )
      {
        
      }
      
      // Fix this for new calibration! //second->recalibrate_by_eqn( info.coefficients, info.deviationpairs, info.type );
    }
    
    if( !!back
        && m_applyTo[BackIndex]->isChecked() && back!=foreground )
    {
      const CalibrationInformation &info = m_originalCal[BackIndex];
      assert( 0 );
      //back->shiftPeaksForRecalibration( disp_back->calibration_coeffs(),
      //                                  disp_back->deviation_pairs(),
      //                                  disp_back->energy_calibration_model(),
      //                                  info.coefficients,
      //                                  info.deviationpairs,
      //                                  info.type );
      // Fix this for new calibration! //back->recalibrate_by_eqn( info.coefficients, info.deviationpairs, info.type );
    }

    //Have to treat the foreground measurment a bit special since the PeakModel
    //  is connected ot its peaks values
    if( !!foreground && m_applyTo[ForeIndex]->isChecked() )
    {
      const CalibrationInformation &info = m_originalCal[ForeIndex];
      // Fix this for new calibration! //foreground->recalibrate_by_eqn( info.coefficients, info.deviationpairs, info.type );
      shiftPeaksForEnergyCalibration( m_peakModel,
                                      info.coefficients, info.deviationpairs, info.type,
                                     foreground, SpectrumType::Foreground,
                                     originalPars, olddevpairs, old_calib_type );
    }

    // Redraw everything
    m_interspec->refreshDisplayedCharts();
    refreshRecalibrator();

    return;
  }//if( action == RevertRecal )

  if( foreground )
  {
    const auto &alldets = foreground->detector_names();
    vector<string> detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
    
    bool displayAll = true;
    for( const auto &det : alldets )
    {
      auto pos = std::find( begin(detectors), end(detectors), det );
      displayAll = displayAll && (pos != end(detectors));
    }

    if( !displayAll )
      passMessage( "Calibration only applied to displayed detectors", "", 1 );
  }//if( foreground )
  
  // Extract the values, noting the exponents.
  vector< float > eqn( sm_numCoefs );
  for( size_t i = 0; i < sm_numCoefs; ++i )
    eqn[i] = m_coefficientDisplay[i]->value();

  try  //begin check validity of equation {note: not a rigourous check}
  {
    const std::shared_ptr<const std::vector<float>> &binning = disp_foreground->channel_energies();
    const int nbin = static_cast<int>( binning->size() );
    vector<float> poly = eqn;
    
    vector<float> frfcoef = eqn;
    if( m_coeffEquationType == SpecUtils::EnergyCalType::Polynomial
       || m_coeffEquationType == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
      frfcoef = SpecUtils::polynomial_coef_to_fullrangefraction( frfcoef, nbin );
    
    bool valid = checkFullRangeFractionCoeffsValid( frfcoef, devpairs, nbin );
    
    if( !valid )
    {
      static_assert( sm_numCoefs >= 3, "" );
      stringstream msg;
      msg << "The coeffiecients {" << eqn[0] << ", " << eqn[1] << ", " << eqn[2]
          << "} are invalid because they will cause higher numbered bins to"
          << " have lower energy values.";

      throw std::runtime_error( msg.str() );
    }//if( (nearend >= end) || (begin >= nearbegin) )
  }catch( std::exception & e )//end codeblock
  {
    passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    m_revert->setDisabled( false );
    
    return;
  }//try / catch to check validity of equation
  
  vector<string> displayed_detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  
  switch( m_coeffEquationType )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    {
      //Recalibrate the foreground
      if( !!foreground && m_applyTo[ForeIndex]->isChecked() )
      {
        const vector<float> origeqn = disp_foreground->calibration_coeffs();
        const vector< pair<float,float> > origdevpars = disp_foreground->deviation_pairs();
        const SpecUtils::EnergyCalType origtype = disp_foreground->energy_calibration_model();
        
        // Fix this for new calibration! //foreground->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType, displayed_detectors, false );
        
        if( action == ApplyRecal )
          shiftPeaksForEnergyCalibration( m_peakModel,
                                          eqn, devpairs, m_coeffEquationType,
                                          foreground, SpectrumType::Foreground,
                                          origeqn, origdevpars, origtype );
      }//if( !!foreground && m_applyTo[ForeIndex]->isChecked() )


      // Recalibrate for the secondary data measurement
      if( !!second && m_applyTo[SecondIndex]->isChecked() )
      {
        const vector<string> &foredet = foreground->detector_names();
        const vector<string> &backdet = second->detector_names();
        const set<string> forenameset( foredet.begin(), foredet.end() );
        const set<string> backnameset( backdet.begin(), backdet.end() );
        
        assert( 0 );
        //if( second != foreground || !m_applyTo[ForeIndex]->isChecked() )
        //  second->shiftPeaksForRecalibration( disp_second->calibration_coeffs(),
        //                                      disp_second->deviation_pairs(),
        //                                      disp_second->energy_calibration_model(),
        //                                      eqn, devpairs, m_coeffEquationType );

        if( forenameset == backnameset )
        {
          // Fix this for new calibration! //second->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType, displayed_detectors, false );
        }else
        {
          // Fix this for new calibration! //second->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType );
        }
      }//if( !!second && m_applyTo[SecondIndex]->isChecked() )
      
      
      // Recalibrate for the background
      if( !!back && m_applyTo[BackIndex]->isChecked() )
      {
        const vector<string> &foredet = foreground->detector_names();
        const vector<string> &backdet = back->detector_names();
        const set<string> forenameset( foredet.begin(), foredet.end() );
        const set<string> backnameset( backdet.begin(), backdet.end() );

        assert( 0 );
        //if( (back != foreground || !m_applyTo[ForeIndex]->isChecked())
        //    && (back != second || !m_applyTo[SecondIndex]->isChecked()) )
        //  back->shiftPeaksForRecalibration( disp_back->calibration_coeffs(),
        //                                    disp_back->deviation_pairs(),
        //                                    disp_back->energy_calibration_model(),
        //                                    eqn, devpairs, m_coeffEquationType );
        if( forenameset == backnameset )
        {
          // Fix this for new calibration! //back->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType, displayed_detectors, false  );
        }else
        {
          // Fix this for new calibration! //back->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType );
        }//if( forenameset == backnameset )
        
      }//if( !!back && m_applyTo[BackIndex]->isChecked() )
      
      
      const bool keepCurrentEnergyRange = true;
      m_interspec->refreshDisplayedCharts();

      break;
    }//case Polynomial, FullRangeFraction:


    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
    {
      passMessage( "Recalibration can only be done for polynomial or full range"
                   " fraction calibrated spectrums", "engageRecalibration",
                   WarningWidget::WarningMsgHigh );
      break;
    }//case LowerChannelEdge and InvalidEquationType
  }//switch( m_coeffEquationType )


  // Recalibrate for the continuum?
//  m_interspec->displayTimeSeriesData( false ); // This can be removed, right?
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  const double newForSum    = ((!!foreground) ? foreground->deep_gamma_count_sum() : 0.0);
  const double newBackSum   = ((!!back) ? back->deep_gamma_count_sum() : 0.0);
  const double newSecondSum = ((!!second) ? second->deep_gamma_count_sum() : 0.0);
  
  char buffer[512];
  if( fabs(origForSum-newForSum) > 1.0E-6*std::max(fabs(origForSum),fabs(newForSum)) )
  {
    snprintf( buffer, sizeof(buffer),
              "Gamma sum of foreground changed after calibration, "
              "pre=%1.8e, post=%1.8e", origForSum, newForSum );
    log_developer_error( __func__, buffer );
  }
  
  if( fabs(origBackSum-newBackSum) > 1.0E-6*std::max(fabs(origBackSum),fabs(newBackSum)) )
  {
    snprintf( buffer, sizeof(buffer),
             "Gamma sum of background changed after calibration, "
             "pre=%1.8e, post=%1.8e", origBackSum, newBackSum );
    log_developer_error( __func__, buffer );
  }
  
  if( fabs(origSecondSum-newSecondSum) > 1.0E-6*std::max(fabs(origSecondSum),fabs(newSecondSum)) )
  {
    snprintf( buffer, sizeof(buffer),
             "Gamma sum of second foreground changed after calibration, "
             "pre=%1.8e, post=%1.8e", origSecondSum, newSecondSum );
    log_developer_error( __func__, buffer );
  }
#endif
}//void engageRecalibration( RECAL_ACTION action )


void Recalibrator::recalibrateByPeaks()
{
  try
  {
    const size_t npeaks = m_peakModel->npeaks();

    int num_coeff_fit = 0;
    for( size_t i = 0; i < sm_numCoefs; ++i )
      num_coeff_fit += m_fitFor[i]->isChecked();
              
    if( num_coeff_fit < 1 )
    {
      const char *msg = "You must select at least one coefficient to fit for";
      passMessage( msg, "recalibrateByPeaks", WarningWidget::WarningMsgHigh );
      return;
    }//if( num_coeff_fit < 1 )

    std::shared_ptr<const SpecMeas> meas;
    std::shared_ptr<const Measurement> displ_meas;
    
    
    for( int i = 0; (!meas || !displ_meas) && (i < 3); ++i )
    {
      if( m_applyTo[i]->isChecked() )
      {
        const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
        meas = m_interspec->measurment(type);
        displ_meas = m_interspec->displayedHistogram(type);
      }//if( m_applyTo[i]->isChecked() )
    }//for( int i = 0; (!meas || !displ_meas) && (i < 3); ++i )
    
    if( !meas || !displ_meas )
      throw std::runtime_error( "ErrorMsgNo spectrum available for calibration" );

    const std::shared_ptr<const std::vector<float>> &binning = displ_meas->channel_energies();
    
    if( !binning || binning->size() < 16 )
      throw std::runtime_error( "ErrorMsgNot enough binns for rebinning" );

    const size_t nbin = binning->size();
    const SpecUtils::EnergyCalType calibration_type = displ_meas->energy_calibration_model();

    vector<float> calib_coefs = displ_meas->calibration_coeffs();

    switch( calibration_type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      break;

      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
      {
        passMessage( "Changing calibration to be full width fraction, "
                     "some information will be lost.",
                     "", WarningWidget::WarningMsgHigh );

        const float lower_energy = binning->at(0);
        const float upper_energy = 2*binning->at(nbin-1) - binning->at(nbin-2);
        const float channel_width = (upper_energy-lower_energy) / nbin;

        vector<float> poly_eqn( 2, 0.0 );
        poly_eqn[0] = lower_energy - 0.5*channel_width;
        poly_eqn[1] = channel_width;

        std::shared_ptr<SpecMeas> foreground = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
        std::shared_ptr<const Measurement> displ_foreground
                                = m_interspec->displayedHistogram(SpectrumType::Foreground);
        if( foreground && displ_foreground )
        {
          std::shared_ptr<SpecMeas> second = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
          std::shared_ptr<SpecMeas> back = m_interspec->measurment(SpecUtils::SpectrumType::Background);
          const std::vector<std::pair<float,float>> &devpairs = displ_foreground->deviation_pairs();
        
          // Fix this for new calibration! //if( m_applyTo[ForeIndex]->isChecked() )
          // Fix this for new calibration! //  foreground->rebin_by_eqn( poly_eqn, devpairs, SpecUtils::EnergyCalType::Polynomial );

          // Fix this for new calibration! //if( !!second && m_applyTo[SecondIndex]->isChecked() )
          // Fix this for new calibration! //  second->rebin_by_eqn( poly_eqn, devpairs, SpecUtils::EnergyCalType::Polynomial );

          // Fix this for new calibration! //if( !!back && m_applyTo[BackIndex]->isChecked() )
          // Fix this for new calibration! //  back->rebin_by_eqn( poly_eqn, devpairs, SpecUtils::EnergyCalType::Polynomial );

          m_coeffEquationType = SpecUtils::EnergyCalType::FullRangeFraction;
          calib_coefs = SpecUtils::polynomial_coef_to_fullrangefraction( poly_eqn, nbin );
          
          updateCoefLabels();
        }//if( !!foreground )
        
        break;
      }//case LowerChannelEdge and InvalidEquationType
    }//switch( calibration_type )

    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contain the error of the fit mean for that peak
    vector<EnergyCal::RecalPeakInfo> peakInfos;

    for( size_t peakn = 0; peakn < npeaks; ++peakn )
    {
      const PeakDef &peak = m_peakModel->peak( peakn );

      if( peak.useForCalibration() )
      {
        const double wantedEnergy = peak.gammaParticleEnergy();

        EnergyCal::RecalPeakInfo peakInfo;
        peakInfo.peakMean = peak.mean();
        peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
        if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
          peakInfo.peakMeanUncert = 0.5;
        
        peakInfo.photopeakEnergy = wantedEnergy;
        if( m_coeffEquationType == SpecUtils::EnergyCalType::FullRangeFraction )
        {
          peakInfo.peakMeanBinNumber = SpecUtils::find_fullrangefraction_channel( peak.mean(),
                                           calib_coefs, nbin,
                                           displ_meas->deviation_pairs(),
                                           0.001f );
        }else
        {
          const vector<float> fwfcoef = SpecUtils::polynomial_coef_to_fullrangefraction( calib_coefs, nbin );
          peakInfo.peakMeanBinNumber = SpecUtils::find_fullrangefraction_channel( peak.mean(),
                                                                  fwfcoef, nbin,
                                                                  displ_meas->deviation_pairs(),
                                                                  0.001f );
        }

        if( IsNan(peakInfo.peakMeanBinNumber)
            || IsInf(peakInfo.peakMeanBinNumber) )
          throw runtime_error( "Invalid result fromm "
                               "find_fullrangefraction_channel(...)" );

        peakInfos.push_back( peakInfo );
      }//if( peak.useForCalibration() )
    }//for( int col = 0; col < numModelCol; ++col )

    if( num_coeff_fit > static_cast<int>( peakInfos.size() ) )
    {
      const char *msg = "You must use at least as many peaks associated with "
                        "nuclides, as the number of coefficients you want to "
                        "fit for.";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }//if( order > static_cast<int>( peakInfos.size() ) )

    
    bool fit_coefs = false;
    vector<double> parValues, parErrors;
    
    try
    {
      //We can currently only so this if there are no non-linear devication pairs...
      vector<float> lls_fit_coefs = calib_coefs, lls_fit_coefs_uncert;
      
      //Convert eqation type to polynomial incase there are any fixed paramters.
      if( m_coeffEquationType == SpecUtils::EnergyCalType::FullRangeFraction )
        lls_fit_coefs = SpecUtils::fullrangefraction_coef_to_polynomial(calib_coefs, binning->size() );
      lls_fit_coefs.resize( sm_numCoefs, 0.0f );
      
      vector<bool> fitfor( lls_fit_coefs.size(), false );
      
      for( size_t i = 0; i < calib_coefs.size() && i < sm_numCoefs; ++i )
        fitfor[i] = m_fitFor[i]->isChecked();
      
      const double chi2 = EnergyCal::fit_energy_cal_poly( peakInfos, fitfor,
                                              displ_meas->deviation_pairs(),
                                              lls_fit_coefs, lls_fit_coefs_uncert );
      
      stringstream msg;
      msg << "\nfit_energy_cal_poly gave chi2=" << chi2 << " with coefs={";
      for( size_t i = 0; i < lls_fit_coefs.size(); ++i )
        msg << lls_fit_coefs[i] << "+-" << lls_fit_coefs_uncert[i] << ", ";
      msg << "}\n";
      cout << msg.str() << endl;
      
      
      switch( m_coeffEquationType )
      {
        case SpecUtils::EnergyCalType::Polynomial:
          break;
          
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          m_coeffEquationType = SpecUtils::EnergyCalType::Polynomial;
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          lls_fit_coefs = SpecUtils::polynomial_coef_to_fullrangefraction(lls_fit_coefs, binning->size());
          lls_fit_coefs_uncert = SpecUtils::polynomial_coef_to_fullrangefraction(lls_fit_coefs_uncert, binning->size());
          break;
          
        case SpecUtils::EnergyCalType::InvalidEquationType:
        case SpecUtils::EnergyCalType::LowerChannelEdge:
          assert( 0 );  //shouldnt ever get here
          break;
      }//switch( m_coeffEquationType )
      
      
      parValues.resize( lls_fit_coefs.size() );
      parErrors.resize( lls_fit_coefs.size() );
      for( size_t i = 0; i < lls_fit_coefs.size(); ++i )
      {
        parValues[i] = lls_fit_coefs[i];
        parErrors[i] = lls_fit_coefs_uncert[i];
      }
      
      fit_coefs = true;
    }catch( std::exception &e )
    {
      cerr << "fit_energy_cal_poly threw: " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      char buffer[512] = { '\0' };
      snprintf( buffer, sizeof(buffer)-1, "fit_energy_cal_poly threw: %s", e.what() );
      log_developer_error( __func__, buffer );
#endif
      fit_coefs = false;
    }//try / catch fit for coefficents using least linear squares
    
    
    if( !fit_coefs )
    {
      if( calib_coefs.size() < sm_numCoefs )
        calib_coefs.resize( sm_numCoefs, 0.0 );
     
      vector<bool> fitfor( calib_coefs.size(), false );
      for( size_t i = 0; i < sm_numCoefs; ++i )
        fitfor[i] = m_fitFor[i]->isChecked();
      
      const size_t nchannel = displ_meas->num_gamma_channels();
      const auto &devpairs = displ_meas->deviation_pairs();
      
      std::string warning_msg;
      std::vector<float> coefs, coefs_uncert;
      
      EnergyCal::fit_energy_cal_iterative( peakInfos, nchannel, m_coeffEquationType, fitfor,
                                          calib_coefs, devpairs, coefs, coefs_uncert, warning_msg );
      
      if( warning_msg.size() )
        passMessage( warning_msg, "", WarningWidget::WarningMsgHigh );
      
      parValues.resize( coefs.size() );
      parErrors.resize( coefs.size() );
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        parValues[i] = coefs[i];
        parErrors[i] = coefs_uncert[i];
      }
    }//if( !fit_coefs )
    
    
    for( int i = 0; i < parValues.size(); ++i )
    {
      m_coefficientDisplay[i]->setValue( parValues[i] );
      m_uncertainties[i] = parErrors[i];
    }//for( size_t i = 0; i < parValues.size(); ++i )
    
    engageRecalibration( ApplyRecal );
  }catch( std::exception &e )
  {
    string exceptionmsg = e.what();
    string msg = "Failed calibration by fitting peak means.";
    
    if( SpecUtils::starts_with( exceptionmsg, "ErrorMsg" ) )
      msg = exceptionmsg.substr(8);
    
    cerr << "Recalibrator::recalibrateByPeaks():\n\tCaught: " << exceptionmsg << endl;
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void recalibrateByPeaks()


void Recalibrator::updateCoefLabels()
{
  for( size_t i = 0; i < sm_numCoefs; ++i )
  {
    string txt;
    switch( m_coeffEquationType )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        if( i == 0 )
          txt = "keV";
        else if( i == 1 )
          txt = "<sup>keV</sup>/<sub>chnl</sub>";
        else
          txt = "<sup>keV</sup>/<sub>chnl<sup>" + std::to_string(i) + "</sup></sub>";
      break;
          
      case SpecUtils::EnergyCalType::FullRangeFraction:
        if( i == 0 )
          txt = "keV";
        else if( i == 1 )
          txt = "<sup>keV</sup>/<sub>FWF</sub>";
        else
          txt = "<sup>keV</sup>/<sub>FWF<sup>" + std::to_string(i) + "</sup></sub>";
      break;
          
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      break;
    }//switch( m_coeffEquationType )
      
    m_termSuffix[i]->setText( txt );
  }//for( size_t i = 0; i < sm_numCoefs; ++i )
};//void updateCoefLabels()

          
void Recalibrator::refreshRecalibrator()
{
  //Needed to check, because Recalibrator always exists whether tool tabs are
  //  shown or not.  refreshRecalibrator() will be called again when it is
  //  visible
  //wcjohns 20140802: I'm not sure this is a great practice... should at least
  //  overide setHidden() to call refreshRecalibrator() to be sure.
  if( isHidden() )
    return;
  
  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = GraphicalRecalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  std::shared_ptr<SpecMeas> foreground, background, second;
  foreground = m_interspec->measurment(SpectrumType::Foreground);
  background = m_interspec->measurment(SpectrumType::Background);
  second     = m_interspec->measurment(SpectrumType::SecondForeground);
  
  std::shared_ptr<const Measurement> displ_foreground, displ_background, displ_second;
  displ_foreground = m_interspec->displayedHistogram(SpectrumType::Foreground);
  displ_background = m_interspec->displayedHistogram(SpectrumType::Background);
  displ_second     = m_interspec->displayedHistogram(SpectrumType::SecondForeground);
  
  const vector<string> detectors = m_interspec->detectorsToDisplay(SpectrumType::Foreground);
  const set<int> samples = m_interspec->displayedSamples(SpectrumType::Foreground);
  const auto energy_cal = foreground ? foreground->suggested_sum_energy_calibration(samples, detectors) : nullptr;
  
  if( m_convertToPolyDialog )
  {
    AuxWindow::deleteAuxWindow( m_convertToPolyDialog );
    m_convertToPolyDialog = nullptr;
  }
  
  if( energy_cal && displ_foreground )
  {
    assert( energy_cal == displ_foreground->energy_calibration() );
  }
  
  // Limiting the data binning size protects against rebinning low channel count
  // detectors that it doesnt make sense to recalibrate, because you wouldnt be
  // able to see peaks anyway
  if( !foreground || !displ_foreground || !energy_cal
      || energy_cal->num_channels() < 16 || (!!foreground
        && (displ_foreground->energy_calibration_model()==SpecUtils::EnergyCalType::LowerChannelEdge
            || displ_foreground->energy_calibration_model()==SpecUtils::EnergyCalType::InvalidEquationType)) )
  {
    m_fitCoefButton->disable();
    m_multiFileButton->disable();
    
    m_applyToLabel->hide();
    for( size_t i = 0; i < sm_numCoefs; ++i )
    {
      m_coefficientDisplay[i]->setValue( 0.0 );
      m_fitFor[i]->hide();
      m_applyTo[i]->hide();
      m_coefficientDisplay[i]->hide();
      m_termSuffix[i]->hide();
    }
    
    if( energy_cal && foreground
        && energy_cal->num_channels() >= 16
        && energy_cal->type()==SpecUtils::EnergyCalType::LowerChannelEdge )
    {
      m_convertToPolynomialLabel->show();
      m_convertToPolynomial->show();
      m_polyConverDiv->show();
    }else
    {
      m_convertToPolynomialLabel->hide();
      m_convertToPolynomial->hide();
      m_polyConverDiv->hide();
    }
    
    return;
  }//if( we cant recalibrate this detector anyway )
  
  if( m_coefficientDisplay[0]->isHidden() )
  {
    m_convertToPolynomialLabel->hide();
    m_convertToPolynomial->hide();
    m_polyConverDiv->hide();
    m_fitCoefButton->enable();
    m_applyToLabel->show();
    for( size_t i = 0; i < sm_numCoefs; ++i )
    {
      m_fitFor[i]->show();
      m_applyTo[i]->show();
      m_coefficientDisplay[i]->show();
      m_termSuffix[i]->show();
    }
  }//if( m_coefficientDisplay[0]->isHidden() )
  
  m_applyTo[ForeIndex]->enable();
  m_applyTo[BackIndex]->enable();
  m_applyTo[SecondIndex]->enable();

  if( energy_cal )
  {
    assert( energy_cal->channel_energies() );
    assert( energy_cal->channel_energies()->size() > 2 );
  }//if( energy_cal )
  
  const size_t nbin = energy_cal->num_channels();
  const vector<float> &binning = *energy_cal->channel_energies();

  // Store the original information, then create a copy and operate on the copy. Note objects pointed to
  // in other places by m_dataMeasurement etc. may be used by other locations in the code.
  // The backup is made by copy instead of via pointer shuffling.

//  m_originalPeakDefs.clear();
//  const std::deque< PeakShrdPtr > &peaks = m_peakModel->peaks();
//  for( const PeakShrdPtr &p : peaks )
//    m_originalPeakDefs.push_back( *p );

  for( int i = 0; i < 3; ++i )
  {
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
    
    std::shared_ptr<SpecMeas> meas = m_interspec->measurment(type);
    std::shared_ptr<const Measurement> displ_meas
                                       = m_interspec->displayedHistogram(type);
    if( !displ_meas )
    {
      m_originalCal[i].reset();
    }else
    {
      m_originalCal[i].type = displ_meas->energy_calibration_model();
      m_originalCal[i].coefficients = displ_meas->calibration_coeffs();
      m_originalCal[i].deviationpairs = displ_meas->deviation_pairs();
      
      m_originalCal[i].sample_numbers = m_interspec->displayedSamples( type );
      // Fix this for new calibration! //m_originalCal[i].detectors_numbers = m_interspec->displayedDetectorNumbers();
    }
  }//for( loop over SpecUtils::SpectrumType );
  
  
  vector<float> equationCoefficients = displ_foreground->calibration_coeffs();
  m_coeffEquationType = m_originalCal[ForeIndex].type;
  if( m_coeffEquationType == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
    m_coeffEquationType = SpecUtils::EnergyCalType::Polynomial;
  
  switch( m_coeffEquationType )
  {
    // Note: fall-through intentional
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    break;

    case SpecUtils::EnergyCalType::LowerChannelEdge:
    //Note: we will never make it here!
    case SpecUtils::EnergyCalType::InvalidEquationType:
      equationCoefficients.clear();
      equationCoefficients.resize( sm_numCoefs, 0.0f );
      if( binning.size() != 0 )
      {
        equationCoefficients[0] = binning[0];
        equationCoefficients[1] = ( binning[ binning.size() - 1 ] - binning[0] ) / nbin;
      }

    break;
  }//switch( m_coeffEquationType )

  updateCoefLabels();
  
  // Set up the little tick/spin/whatever boxes
  for( size_t i = 0; i < sm_numCoefs; ++i )
  {
    m_coefficientDisplay[i]->setDisabled( false );
    const float val = (i < equationCoefficients.size()) ? equationCoefficients[i] : 0.0f;
    m_coefficientDisplay[i]->setValue( val );
          
    float stepsize = 1.0f;
    switch( m_coeffEquationType )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        stepsize = 1.0f / std::pow(nbin,i);
      break;
          
      case SpecUtils::EnergyCalType::FullRangeFraction:
        stepsize = 1.0;
      break;
          
      case SpecUtils::EnergyCalType::InvalidEquationType:
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        stepsize = 0.0f;
      break;
    }//switch( m_coeffEquationType )
      
    m_coefficientDisplay[i]->setSingleStep( stepsize );
  }//for( int i = 0; i < sm_numCoefs; ++i )
  
  m_devPairs->setDeviationPairs( m_originalCal[ForeIndex].deviationpairs );
  
  // Disable whichever check boxes correspond to null
  m_applyTo[ForeIndex]->setDisabled( !foreground );
  m_applyTo[ForeIndex]->setChecked( !!foreground);
  
  m_applyTo[BackIndex]->setDisabled( !background );
  m_applyTo[BackIndex]->setChecked( !!background);
  
  m_applyTo[SecondIndex]->setDisabled( !second );
  m_applyTo[SecondIndex]->setChecked( !!second);
  
  checkIfCanFitCoefs();
  checkIfCanMultiFileCal();
}//void refreshRecalibrator()




DevPair::DevPair( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    //m_energy( new WDoubleSpinBox() ),
    //m_offset( new WDoubleSpinBox() ),
    m_energy( new NativeFloatSpinBox() ),
    m_offset( new NativeFloatSpinBox() ),
    m_delete( new WContainerWidget() )
{
  WGridLayout* layout = new WGridLayout();
  layout->setContentsMargins(0, 0, 0, 0);
  layout->setVerticalSpacing( 0 );
  
  setLayout(layout);
  setStyleClass( "DevPair" );
  
  layout->addWidget(m_energy,0,0);
  layout->addWidget(m_offset,0,1);
  layout->addWidget(m_delete,0,2);
  layout->setColumnStretch(0, 1);
  layout->setColumnStretch(1, 1);
  layout->setColumnStretch(2, 0);
          
  m_energy->setStyleClass( "DevEnergy" );
  m_offset->setStyleClass( "DevOffset" );
  m_energy->setPlaceholderText( "Energy (keV)" );
  m_offset->setPlaceholderText( "Offset (keV)" );
  m_energy->setValueText( "" );
  m_offset->setValueText( "" );
          
  m_delete->addStyleClass( "Wt-icon DeleteDevPair" );
}//DevPair constructor


void DevPair::setDevPair( const std::pair<float,float> &d )
{
  m_energy->setValue( d.first );
  m_offset->setValue( d.second );
  
  auto printval = []( float val ) -> std::string {
    char buffer[64];
    const float fraction = val - std::floor(val);
    if( fraction == 0.0 )
      snprintf( buffer, sizeof(buffer), "%.0f", val );
    else if( fraction == 0.1 )
      snprintf( buffer, sizeof(buffer), "%.1f", val );
    else
      snprintf( buffer, sizeof(buffer), "%.2f", val );
    return buffer;
  };
 
  m_energy->setText( printval(d.first) );
  m_offset->setText( printval(d.second) );
}


std::pair<float,float> DevPair::devPair() const
{
  const float energy = m_energy->value();
  const float offset = m_offset->value();
  return std::make_pair( energy, offset );
}


void DevPair::visuallyIndicateChanged()
{
  doJavaScript( "$('#" + id() + "').fadeIn(100).fadeOut(100).fadeIn(100).fadeOut(100).fadeIn(100);" );
}//void visuallyIndicateChanged();


DeviationPairDisplay::DeviationPairDisplay( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_pairs( NULL )
{
  addStyleClass( "DevPairDisplay" );
  WLabel *title = new WLabel( "Deviation Pairs", this );
  title->setStyleClass( "Wt-itemview Wt-header Wt-label DevPairTitle" );
  title->setInline( false );

  m_pairs = new WContainerWidget(this);
  m_pairs->setStyleClass( "DevPairsContainer" );
          
  auto footer = new WContainerWidget(this);
  footer->setStyleClass( "DevPairsFooter" );
  
  auto addBtn = new WContainerWidget( footer );
  addBtn->addStyleClass( "Wt-icon AddDevPair" );
  addBtn->clicked().connect( boost::bind(&DeviationPairDisplay::newDevPair, this, true) );
  addBtn->setToolTip( "Add another deviation pair" );
}//DeviationPairDisplay constructor


void DeviationPairDisplay::setDeviationPairs( vector< pair<float,float> > d )
{
  m_pairs->clear();
  std::sort( d.begin(), d.end() );
  
  for( size_t i = 0; i < d.size(); ++i )
  {
    DevPair *dev = newDevPair( false );
    dev->setDevPair( d[i] );
  }//for( size_t i = 0; i < d.size(); ++i )
  
  sortDisplayOrder(false);
}//setDeviationPairs(...)


void DeviationPairDisplay::sortDisplayOrder( const bool indicateVisually )
{
  vector<size_t> sort_indices;
  vector<DevPair *> displays;
  
  vector<float> offsets;
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    DevPair *p = dynamic_cast<DevPair *>( t );
    if( p )
    {
      const size_t index = sort_indices.size();
      offsets.push_back( p->devPair().first );
      displays.push_back( p );
      sort_indices.push_back( index );
    }
  }//for( WWidget *t : childs )
  
  std::stable_sort( sort_indices.begin(), sort_indices.end(),
                    index_compare_assend<vector<float>&>(offsets) );

  bool order_changed = false;
  for( size_t i = 0; i < sort_indices.size(); ++i )
    order_changed |= (sort_indices[i] != i );
  
  if( !order_changed )
    return;
          
  for( size_t i = 0; i < displays.size(); ++i )
    m_pairs->removeWidget( displays[i] );
    
  for( size_t i = 0; i < displays.size(); ++i )
  {
    DevPair *p = displays[ sort_indices[i] ];
    m_pairs->addWidget( p );
    if( indicateVisually && (i != sort_indices[i]) )
      p->visuallyIndicateChanged();
  }
}//void sortDisplayOrder()


void DeviationPairDisplay::emitChanged()
{
  m_changed.emit();
}


vector< pair<float,float> > DeviationPairDisplay::deviationPairs() const
{
  vector< pair<float,float> > answer;
  const vector<WWidget *> childs = m_pairs->children();
  for( const WWidget *t : childs )
  {
    const DevPair *p = dynamic_cast<const DevPair *>( t );
    if( p )
      answer.push_back( p->devPair() );
  }//for( WWidget *t : childs )
  
  std::sort( answer.begin(), answer.end() );
  
  return answer;
}//deviationPairs()


void DeviationPairDisplay::removeDevPair( DevPair *devpair )
{
  if( !devpair )
    return;
  
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    if( devpair == dynamic_cast<DevPair *>( t ) ) //dynamic_cast prob not necessary
    {
      delete devpair;
      sortDisplayOrder(false);
      emitChanged();
      return;
    }
  }//for( WWidget *t : childs )
}//removeDevPair(...)


void DeviationPairDisplay::setInvalidValues()
{
  addStyleClass( "InvalidDevPairs" );
}


void DeviationPairDisplay::setValidValues()
{
  removeStyleClass( "InvalidDevPairs" );
}


DevPair *DeviationPairDisplay::newDevPair( const bool emitChangedNow )
{
  DevPair *dev = new DevPair( m_pairs );
  dev->m_delete->clicked().connect( boost::bind( &DeviationPairDisplay::removeDevPair, this, dev ) );
  //dev->m_energy->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  //dev->m_offset->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_energy->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_offset->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_energy->blurred().connect( boost::bind(&DeviationPairDisplay::sortDisplayOrder, this, true) );
  
  if( emitChangedNow )
    emitChanged();
  return dev;
}//newDevPair()

Wt::Signal<> &DeviationPairDisplay::changed()
{
  return m_changed;
}



void Recalibrator::CalibrationInformation::reset()
{
  type = SpecUtils::EnergyCalType::InvalidEquationType;
  coefficients.clear();
  deviationpairs.clear();
  sample_numbers.clear();
  detectors_numbers.clear();
}


Recalibrator::CalibrationInformation::CalibrationInformation()
{
  reset();
}


Recalibrator::GraphicalRecalConfirm::GraphicalRecalConfirm( double lowe,
                                                            double highe,
                                                            Recalibrator *cal,
                                                            time_t lastCal,
                                    GraphicalRecalConfirm::RecalTypes lastType,
                                                            float lastEnergy
                                                            )
: AuxWindow( "Confirm Recalibration",
             (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
               | AuxWindowProperties::SetCloseable | AuxWindowProperties::DisableCollapse) ),
  m_calibrator( cal ), m_typeButtons( NULL), m_foregroundOnly( NULL ),
  m_startE( NULL ), m_finalE( NULL ),
  m_preserveLastCal( NULL ),
  m_lastType( lastType ),
  m_lastEnergy( lastEnergy )
{
  setResizable( false );
  rejectWhenEscapePressed();
  centerWindow();
  show();

  finished().connect( cal, &Recalibrator::deleteGraphicalRecalConfirmWindow );
  
  if( !cal->m_interspec->displayedHistogram(SpectrumType::Foreground) )
  {
    throw runtime_error( "You must have a dispayed spectrum to do a recalubration" );
    return;
  }//if( !spectrum->m_interspec->displayedHistogram(SpectrumType::Foreground) )
  
  WTable *table = new WTable( contents() );
  //    table->setHeaderCount( 1, Wt::Vertical );
  new WLabel( "Original Energy", table->elementAt(0,0) );
  m_startE = new WDoubleSpinBox( table->elementAt(0,1) );
  new WLabel( "Modified Energy", table->elementAt(1,0) );
  m_finalE = new WDoubleSpinBox( table->elementAt(1,1) );
  m_startE->setMaximum( 20000 );
  m_finalE->setMaximum( 20000 );
  setEnergies( lowe, highe );
  
  
  m_typeButtons = new WButtonGroup( contents() );
  WGroupBox *buttonBox = new WGroupBox( "Parameter to adjust", contents() );
  
  for( RecalTypes t = RecalTypes(0); t < NumRecalTypes; t = RecalTypes(t+1) )
  {
    const char *txt = "";
    switch( t )
    {
      case kOffset:       txt = "Offset";             break;
      case kLinear:       txt = "Linear";             break;
//      case kQuadratic:    txt = "Quadratic";          break;
      case kDeviation:    txt = "Add Deviation Pair"; break;
      case NumRecalTypes: txt = "";                   break;
    }//switch( t )
    
    WRadioButton *button = new WRadioButton( txt, buttonBox );
    button->setInline( false );
    m_typeButtons->addButton( button, t );
  }//for( loop over recal types )
  m_typeButtons->setSelectedButtonIndex( kLinear );
  
  
  InterSpec *viewer = m_calibrator->m_interspec;
  std::shared_ptr<const Measurement> secondaryH  = viewer->displayedHistogram(SpecUtils::SpectrumType::SecondForeground);
  std::shared_ptr<const Measurement> backgroundH = viewer->displayedHistogram(SpecUtils::SpectrumType::Background);
  
  if( secondaryH || backgroundH )
  {
    m_foregroundOnly = new WCheckBox( "Foreground Only", contents() );
    m_foregroundOnly->setInline( false );
    m_foregroundOnly->setAttributeValue( "style",
                "margin-left:1.25em;margin-top:0.25em;margin-bottom=0.25em;" );
  }//if( secondaryH || backgroundH )
  
  time_t currenttime;
  time( &currenttime );
  const double secondsLastCal = difftime( currenttime, lastCal );
  if( secondsLastCal < 120.0
     && (m_lastType==kOffset || m_lastType==kLinear)
     && (m_lastEnergy > 0.0) )
  {
    char msg[128];
    snprintf(msg, sizeof(msg), "Preserve %.1f keV Cal.", m_lastEnergy );
    m_preserveLastCal = new WCheckBox( msg, contents() );
    m_preserveLastCal->setInline( false );
    m_preserveLastCal->setStyleClass( "PreserveLastCalCb" );
    
    const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", cal->m_interspec );

    HelpSystem::attachToolTipOn( m_preserveLastCal,"This is only possible if a offset or"
                                " linear term adjustment was previously"
                                " made within the last 2 minutes.", showToolTipInstantly );
    m_preserveLastCal->setChecked();
    buttonBox->disable();
    m_preserveLastCal->checked().connect( buttonBox, &WWidget::disable );
    m_preserveLastCal->unChecked().connect( buttonBox, &WWidget::enable );
  }//if( preserve last cal possibly )
  
  const shared_ptr<const SpecMeas> meas = cal->m_interspec->measurment(SpectrumType::Foreground);
  const vector<string> detectors = cal->m_interspec->detectorsToDisplay(SpectrumType::Foreground);
  
  bool displayAll = true;
  if( meas )
  {
    for( const string &name : meas->detector_names() )
    {
      const auto pos = std::find( begin(detectors), end(detectors), name );
      displayAll = displayAll && (pos != end(detectors));
    }
  }//if( meas )
  
  if( !displayAll )
  {
    WText *t = new WText( "(Calibration only applied to<br>displayed detectors)",
                         XHTMLUnsafeText, contents() );
    t->setInline( false );
    t->setAttributeValue( "style", "color: #737373; width: auto; text-align: center;" );
  }//if( !displayAll )
  
  
  
  
  
  AuxWindow::addHelpInFooter( footer(), "graphical-recal-dialog" );
  
  
  WPushButton *button = new WPushButton( "Cancel", footer() );
  button->clicked().connect( this, &AuxWindow::hide );
  
  button = new WPushButton( "Accept", footer()  );
  button->setIcon( "InterSpec_resources/images/accept.png" );
  button->clicked().connect( this, &GraphicalRecalConfirm::apply );
  
  //Could place dialog near where the mouse is
  //  stringstream moveJs;
  //  moveJs << "var el=" << this->jsRef() << ";el.style.top='" << y << "px';"
  //         << "el.style.left='" << x << "px';";
  //  doJavaScript( moveJs.str();
}//GraphicalRecalConfirm


Recalibrator::GraphicalRecalConfirm::~GraphicalRecalConfirm()
{
}


bool Recalibrator::checkFullRangeFractionCoeffsValid( const std::vector<float> &eqn,
                                              const std::vector< std::pair<float,float> > &devpairs,
                                              size_t nbin )
{
  return SpecUtils::calibration_is_valid( SpecUtils::EnergyCalType::FullRangeFraction,
                                              eqn, devpairs, nbin );

//  const float nearend   = SpecUtils::fullrangefraction_energy( nbin-2, eqn, nbin, devpairs );
//  const float end       = SpecUtils::fullrangefraction_energy( nbin-1, eqn, nbin, devpairs );
//  const float begin     = SpecUtils::fullrangefraction_energy( 0,      eqn, nbin, devpairs );
//  const float nearbegin = SpecUtils::fullrangefraction_energy( 1,      eqn, nbin, devpairs );
//  return !( (nearend >= end) || (begin >= nearbegin) );
}//checkFullRangeFractionCoeffsValid



void Recalibrator::GraphicalRecalConfirm::apply()
{
  const float startE = static_cast<float>(m_startE->value());
  const float finalE = static_cast<float>(m_finalE->value());
  const RecalTypes type = RecalTypes( m_typeButtons->selectedButtonIndex() );
  
  if( startE == finalE )
  {
    finished().emit(WDialog::Accepted);
    return;
  }//if( startE == finalE )
  
  
  //Need to implement recalibrating logic, and also logic to check how long
  //  ago user last calibrated this spectra...
  InterSpec *viewer = m_calibrator->m_interspec;
  std::shared_ptr<SpecMeas> foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  std::shared_ptr<const Measurement> displ_foreground
                                      = viewer->displayedHistogram(SpectrumType::Foreground);
  if( !foreground || !displ_foreground )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  
  float shift = finalE - startE;
  vector<float> eqn = displ_foreground->calibration_coeffs();
  SpecUtils::EnergyCalType eqnType = displ_foreground->energy_calibration_model();
  std::vector<std::pair<float,float>> deviationPairs = displ_foreground->deviation_pairs();
  const std::vector<std::pair<float,float>> origdev = deviationPairs;
  const size_t nbin = foreground->num_gamma_channels();
  
  switch( eqnType )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      eqnType = SpecUtils::EnergyCalType::FullRangeFraction;
      eqn = SpecUtils::polynomial_coef_to_fullrangefraction( eqn, nbin );
    break;
    case SpecUtils::EnergyCalType::FullRangeFraction:
    break;
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      throw runtime_error( "You cannot recalibrate spectrums with an unknown "
                           "calibration or one defined by lower bin energees" );
  }//switch( eqnType )
  
  const vector<float> orig_eqn = eqn;
  
  //m_preserveLastCal will only be valid if m_lastType is kOffset or kLinear,
  // and m_lastEnergy > 0.0
  const bool preserveLast = (m_preserveLastCal && m_preserveLastCal->isChecked()
                             && (m_lastEnergy!=startE));
  
  if( preserveLast )
  {
    const float lowbin = SpecUtils::find_fullrangefraction_channel( m_lastEnergy, orig_eqn,
                                                  nbin, deviationPairs, 0.001f );
    const float upbin = SpecUtils::find_fullrangefraction_channel( startE, orig_eqn,
                                                  nbin, deviationPairs, 0.001f );
    
    //From a quick empiracal test (so by no means definitive), the below
    //  behaves well for the case deviation pairs are defined.
    const float x1 = upbin / nbin;
    const float x2 = lowbin / nbin;
    const float E1 = eqn[0] + eqn[1]*x1;  // startE, if no deviation pairs
    const float E2 = eqn[0] + eqn[1]*x2;  // m_lastEnergy
    eqn[1] = (E1 - E2 + shift) / (x1-x2);
    eqn[0] = E2 - x2*eqn[1];
  }else
  {
    switch( type )
    {
      case kOffset:
        eqn[0] += shift;
      break;
      
      case kLinear:
      {
        //From a quick empiracal test (so by no means definitive), the below
        //  behaves well for the case deviation pairs are defined.
        const float binnum = SpecUtils::find_fullrangefraction_channel( startE, eqn, nbin,
                                                        deviationPairs, 0.001f );
        const float binfrac = binnum / nbin;
        eqn[1] += (shift / binfrac);
        break;
      }//case kLinear:
            
//      case kQuadratic:
//      break;
      
      case kDeviation:
      {
        if( deviationPairs.empty() )
        {
          deviationPairs.push_back( std::pair<float,float>(0.0f,0.0f) );
          deviationPairs.push_back( std::pair<float,float>(startE,shift) );
        }else
        {
          //We are going to to find the bin number startE cooresponds to, then
          //  find the energy this would have cooresponded to if there were no
          //  deviation pairs, then we will add a deviation
          const float bin = SpecUtils::find_fullrangefraction_channel( startE, orig_eqn, nbin,
                                                        deviationPairs, 0.001f );
          
          float originalE = 0.0;
          const float x = bin/static_cast<float>(nbin);
          for( size_t c = 0; c < eqn.size(); ++c )
            originalE += eqn[c] * pow(x,static_cast<float>(c) );
          shift = finalE - originalE;
          deviationPairs.push_back( std::pair<float,float>(originalE,shift) );
        }//if( deviationPairs.empty() )
      
        std::sort( deviationPairs.begin(), deviationPairs.end() );
        break;
      }//case kDeviation:
      
      case NumRecalTypes:
      break;
    }//switch( type )
  }//if( preserveLast ) / else
  
//  if( eqnType == SpecUtils::EnergyCalType::Polynomial || eqnType == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
//    eqn = fullrangefraction_coef_to_polynomial( eqn, nbin );
  
  bool valid = checkFullRangeFractionCoeffsValid( eqn, deviationPairs, nbin );
  
  if( !valid )
  {
    stringstream msg;
    msg << "The coeffiecients {" << eqn[0] << ", " << eqn[1] << ", " << eqn[2]
        << "} are invalid because they will cause higher numbered bins to"
        << " have lower energy values. Calibration not applied.";
    
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    finished().emit(WDialog::Accepted);
    return;
  }//if( (nearend >= end) || (begin >= nearbegin) )
  
  const vector<string> detectors = m_calibrator->m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  // Fix this for new calibration! //foreground->recalibrate_by_eqn( eqn, deviationPairs, eqnType, detectors, false );

  std::shared_ptr<SpecMeas> background = viewer->measurment(SpecUtils::SpectrumType::Background);
  std::shared_ptr<SpecMeas> secondary = viewer->measurment(SpecUtils::SpectrumType::SecondForeground);
  std::shared_ptr<const Measurement> displ_background
                                      = viewer->displayedHistogram(SpectrumType::Background);
  std::shared_ptr<const Measurement> displ_secondary
                                = viewer->displayedHistogram(SpectrumType::SecondForeground);
  
  if( !m_foregroundOnly || !m_foregroundOnly->isChecked() )
  {    
    if( background && displ_background && (background != foreground) )
    {
      const vector<string> &foredet = foreground->detector_names();
      const vector<string> &backdet = background->detector_names();
      const set<string> forenameset(foredet.begin(),foredet.end());
      const set<string> backnameset(backdet.begin(),backdet.end());
      
      assert( 0 );
      //background->shiftPeaksForRecalibration( displ_background->calibration_coeffs(),
      //                                        displ_background->deviation_pairs(),
      //                                        displ_background->energy_calibration_model(),
      //                                        eqn, deviationPairs, eqnType );

      
      // Fix this for new calibration! //if( forenameset == backnameset )
      // Fix this for new calibration! //  background->recalibrate_by_eqn( eqn, deviationPairs, eqnType, detectors, false );
      // Fix this for new calibration! //else
      // Fix this for new calibration! //  background->recalibrate_by_eqn( eqn, deviationPairs, eqnType );
    }//if( background )
    
    if( secondary && displ_secondary && (secondary != foreground) )
    {
      const vector<string> &foredet = foreground->detector_names();
      const vector<string> &secodet = secondary->detector_names();
      const set<string> forenameset(foredet.begin(),foredet.end());
      const set<string> seconameset(secodet.begin(),secodet.end());
      
      assert( 0 );
      //secondary->shiftPeaksForRecalibration( displ_secondary->calibration_coeffs(),
      //                                       displ_secondary->deviation_pairs(),
      //                                       displ_secondary->energy_calibration_model(),
      //                                       eqn, deviationPairs, eqnType );
      
      // Fix this for new calibration! //if( forenameset == seconameset )
      // Fix this for new calibration! //  secondary->recalibrate_by_eqn( eqn, deviationPairs, eqnType, detectors, false );
      // Fix this for new calibration! //else
      // Fix this for new calibration! //  secondary->recalibrate_by_eqn( eqn, deviationPairs, eqnType );
    }//if( secondary )
  }
//  else
//  {
//    if( !!background )
//      background->rebin_by_eqn( eqn, deviationPairs, eqnType );
//    if( !!secondary )
//      secondary->rebin_by_eqn( eqn, deviationPairs, eqnType );
//  }//if( we should recalibrate background/secondary/continuum )

  
  Recalibrator::shiftPeaksForEnergyCalibration( m_calibrator->m_peakModel,
                         eqn, deviationPairs, eqnType,
                         foreground, SpectrumType::Foreground, orig_eqn, origdev, eqnType );
  
  // Redraw everything
  viewer->refreshDisplayedCharts();
  m_calibrator->refreshRecalibrator();
  m_calibrator->setWasGraphicalRecal( static_cast<int>(type), finalE );

  finished().emit(WDialog::Accepted);
}//void apply()

void Recalibrator::GraphicalRecalConfirm::setEnergies( double xstart,
                                                      double xfinish )
{
  m_startE->setValue( floor(100.0*xstart+0.5)/100.0 );
  m_finalE->setValue( floor(100.0*xfinish+0.5)/100.0 );
}//void setEnergies(...)




Recalibrator::MultiFileModel::MultiFileModel( SpectraFileModel *fileModel, Wt::WObject *parent )
: WAbstractItemModel( parent ),
  m_fileModel( fileModel )
{
  refreshData();
  m_fileModel->rowsInserted().connect( boost::bind(&Recalibrator::MultiFileModel::refreshData, this) );
  m_fileModel->rowsRemoved().connect( boost::bind(&Recalibrator::MultiFileModel::refreshData, this) );
  
  //Should add in a listener here to the PeakModel to see if peaks are
  //  added/removed; this might neccesitate tracking 
}//MultiFileModel constructor


Recalibrator::MultiFileModel::~MultiFileModel()
{
}//~MultiFileModel()


WModelIndex Recalibrator::MultiFileModel::index( int row, int column,
                                              const WModelIndex &parent ) const
{
  if( parent.isValid() )
  {
    if( column < 0 || column >= 3 )
      return WModelIndex();
    
    if( parent.row() < 0 || parent.row() >= static_cast<int>(m_peaks.size()) )
      return WModelIndex();
    
    return createIndex( row, column, ::uint64_t(parent.row()) );
  }//if( parent.isValid() )
  
  if( row < 0 || row >= static_cast<int>(m_peaks.size()) )
    return WModelIndex();
  
  if( column != 0 )
    return WModelIndex();
  return createIndex( row, column, numeric_limits<unsigned int>::max() );
}//index(...)


WModelIndex Recalibrator::MultiFileModel::parent( const WModelIndex &index ) const
{
  if( !index.isValid() )
    return WModelIndex();
  
  const ::uint64_t filenum = index.internalId();
  if( filenum == numeric_limits<unsigned int>::max() )
    return WModelIndex();
  
  if( filenum < m_peaks.size() )
    return createIndex( filenum, 0, numeric_limits<unsigned int>::max() );
  
  return WModelIndex();
}//parent(...)


int Recalibrator::MultiFileModel::rowCount( const WModelIndex &parent ) const
{
  if( !parent.isValid() )
    return static_cast<int>( m_peaks.size() );
  
  if( parent.parent().isValid() )
    return 0;
  
  const int filenum = parent.row();
  if( filenum >= 0 && filenum < static_cast<int>(m_peaks.size()) )
    return static_cast<int>( m_peaks[filenum].size() );
  
  return 0;
}//rowCount(...)


int Recalibrator::MultiFileModel::columnCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 3;
  return 1;
}//columnCount(...)

  
boost::any Recalibrator::MultiFileModel::data( const Wt::WModelIndex &index,
                                               int role ) const
{
  if( (role != Wt::DisplayRole) && (role != Wt::CheckStateRole) )
    return boost::any();
  
  if( !index.isValid() )
    return boost::any();
  
  if( !index.parent().isValid() && index.column() != 0 )
    return boost::any();
  
  const int row = index.row();
  const int col = index.column();
  
  if( index.parent().isValid() )
  {
    const int filenum = index.parent().row();
    if( filenum < 0 || filenum >= static_cast<int>(m_peaks.size()) )
      return boost::any();
    if( row < 0 || row >= static_cast<int>(m_peaks[filenum].size()) )
      return boost::any();
    
    const bool checked = m_peaks[filenum][row].first;
    std::shared_ptr<const PeakDef> peak = m_peaks[filenum][index.row()].second;
    
    if( !peak )
      return boost::any();
    
    if( !peak->parentNuclide() && !peak->xrayElement() && !peak->reaction() )
      return boost::any();
    
    if( role == Wt::CheckStateRole )
    {
      if( col == 0 )
        return boost::any( checked );
      return boost::any();
    }
    
    if( col == 0 )
    {
      if( peak->parentNuclide() )
        return boost::any( WString(peak->parentNuclide()->symbol) );
      if( peak->xrayElement() )
        return boost::any( WString(peak->xrayElement()->symbol) );
      if( peak->reaction() )
        return boost::any( WString(peak->reaction()->name()) );
      return boost::any( WString("--") );
    }//if( col == 0 )
    
    char msg[128];
    msg[0] = '\0';
    
    if( col == 1 )
      snprintf(msg, sizeof(msg), "%.2f keV Detected", peak->mean() );
    
    if( col == 2 )
    {
      try
      {
        const float wantedEnergy = peak->gammaParticleEnergy();
        snprintf(msg, sizeof(msg), "%.2f keV Expected", wantedEnergy );
      }catch(...)
      {
        snprintf(msg, sizeof(msg), "%s", "--" );
      }
    }//if( col == 2 )
    
    if( !strlen(msg) )
      return boost::any();
    return boost::any( WString(msg) );
  }//if( index.parent().isValid() )
  
  if( role != Wt::DisplayRole )
    return boost::any();
  
  if( row < 0 || row >= static_cast<int>(m_peaks.size()) )
    return boost::any();
  
  std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( row );
  if( header )
    return boost::any( header->displayName() );
  return boost::any();
}//data(...)


bool Recalibrator::MultiFileModel::setData( const Wt::WModelIndex &index,
                                            const boost::any &value, int role )
{
  if( role != Wt::CheckStateRole
      || !index.isValid() || !index.parent().isValid() )
    return false;
  
  const int filenum = index.parent().row();
  const int row = index.row();
  const int col = index.column();
  
  if( col != 0 )
    return false;
  if( filenum >= static_cast<int>(m_peaks.size()) )
    return false;
  if( row >= static_cast<int>(m_peaks[filenum].size()) )
    return false;
  
  try
  {
    m_peaks[filenum][row].first = boost::any_cast<bool>( value );
  }catch(...)
  {
    cerr << "Got bad boost::any_cast<bool>( value )" << endl;
    return false;
  }
  
  return true;
}//setData(...)



WFlags<ItemFlag> Recalibrator::MultiFileModel::flags( const WModelIndex &index ) const
{
  if( index.isValid() && index.parent().isValid() && index.column()==0 )
    return WFlags<ItemFlag>(ItemIsUserCheckable);
  return WFlags<ItemFlag>(ItemIsXHTMLText);
}//flags(...)


boost::any Recalibrator::MultiFileModel::headerData( int section,
                                                    Wt::Orientation orientation,
                                                     int role ) const
{
  if( orientation != Wt::Horizontal )
    return boost::any();
  if( role != DisplayRole )
    return boost::any();
  
  if( section == 0 )
    boost::any( WString("Nuclide") );
  if( section == 1 )
    boost::any( WString("Obs. Energy") );
  if( section == 2 )
    boost::any( WString("Photopeak Energy") );
  return boost::any();
}//headerData(...);


void Recalibrator::MultiFileModel::refreshData()
{
  if(!m_peaks.empty())
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks.size())-1 );
    m_peaks.clear();
    endRemoveRows();
  }//if( m_peaks.size() )
  
  vector<std::vector< std::pair<bool,std::shared_ptr<const PeakDef> > > > newdata;
  
  const int nfile = m_fileModel->rowCount();
  for( int filenum = 0; filenum < nfile; ++filenum )
  {
    vector< pair<bool,std::shared_ptr<const PeakDef> > > peaks;
    std::shared_ptr<SpectraFileHeader> header
                                           = m_fileModel->fileHeader( filenum );
    if( !header )
    {
      newdata.push_back( peaks );
      continue;
    }//if( !header )
    
    std::shared_ptr<SpecMeas> spec = header->parseFile();
    
    if( !spec )
    {
      newdata.push_back( peaks );
      continue;
    }//if( !header )
    
    typedef set<int> IntSet;
    typedef std::shared_ptr<const PeakDef> PeakPtr;
    const set< set<int> > peaksamplenums = spec->sampleNumsWithPeaks();
    for( const IntSet &samplnums : peaksamplenums )
    {
      std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > measpeaks;
      measpeaks = spec->peaks( samplnums );
      if( !measpeaks )
        continue;
      
      for( const PeakPtr &peak : *measpeaks )
      {
        const bool use = peak->useForCalibration();
        pair<bool,std::shared_ptr<const PeakDef> > peakpair( use, peak );
        peaks.push_back( peakpair );
      }//for( const PeakPtr &peak : peaks )
      
      newdata.push_back( peaks );
    }//for( const IntSet &samplnums : peaksamplenums )
  }//for( int filenum = 0; filenum < nfile; ++filenum )

  if(!newdata.empty())
  {
    beginInsertRows( WModelIndex(), 0, static_cast<int>(newdata.size()) -1 );
    m_peaks.swap( newdata );
    endInsertRows();
  }//if( newdata.size() )
}//refreshData()


Recalibrator::MultiFileCalibFit::MultiFileCalibFit( Recalibrator *cal )
: AuxWindow( "Multi-File Calibration" ),
  m_calibrator( cal ),
  m_model( 0 ),
  m_use( 0 ),
  m_cancel( 0 ),
  m_fit( 0 ),
  m_fitSumary( 0 )
{
  InterSpec *viewer = m_calibrator->m_interspec;
  SpectraFileModel *fileModel = viewer->fileManager()->model();
  
  m_model = new MultiFileModel( fileModel, this );

  RowStretchTreeView *tree = new RowStretchTreeView();
  tree->setSortingEnabled(false);
  tree->setModel( m_model );
  tree->setColumn1Fixed( false );
  tree->setHeaderHeight( 2 );
  tree->setColumnWidth( 0, 200 );
  tree->setColumnWidth( 1, 150 );
  tree->setColumnWidth( 2, 150 );
  
  WGroupBox *fitFor = new WGroupBox( "Coefficents to fit for" );
  WGridLayout *fitForLayout = new WGridLayout( fitFor );
  
  for( size_t i = 0; i < sm_numCoefs; ++i )
  {
    m_calVal[i] = m_calibrator->m_coefficientDisplay[i]->value();
    m_calUncert[i] = -1.0;
    
    WLabel *label = 0;
    switch( i )
    {
      case 0: label = new WLabel( "Offset" ); break;
      case 1: label = new WLabel( "Linear" ); break;
      case 2: label = new WLabel( "Quadratic" ); break;
      case 3: label = new WLabel( "Quartic" ); break;
      default: label = new WLabel( std::to_string(i+1) + "th" ); break;
    }//switch( i )
    
    m_coefvals[i] = new WLineEdit();
    m_fitFor[i] = new WCheckBox( "Fit" );
    m_coefvals[i]->disable();
    
    fitForLayout->addWidget( label,         i, 0 );
    fitForLayout->addWidget( m_coefvals[i], i, 1 );
    fitForLayout->addWidget( m_fitFor[i],   i, 2 );
  }//for( int i = 0; i < sm_numCoefs; ++i )
  
  fitForLayout->setColumnStretch( 1, 1 );
  m_fitFor[0]->setChecked( true );
  m_fitFor[1]->setChecked( true );
  
 
  m_fitSumary = new WTextArea();
  m_fitSumary->setHeight( 75 );
  m_fitSumary->setMaximumSize( WLength::Auto, 75 );
  
  WContainerWidget *instructions = new WContainerWidget();
  WText *line = new WText( "Select peaks to use from each file then click &quot;Fit&quot;.", instructions );
  line->setInline( false );
  line = new WText( "If satisfied, click &quot;Use&quot; to set calibration for involved files.", instructions );
  line->setInline( false );
  line = new WText( "Calibration will be applied to all files with at least one selected peak.", instructions );
  line->setInline( false );

  
  WGridLayout *layout = stretcher();
  layout->setContentsMargins( 0, 0, 0, 0 );
  
  layout->addWidget( instructions, 0, 0 );
  layout->addWidget( tree, 1, 0 );
  layout->setRowStretch( 1, 1 );
  layout->addWidget( fitFor, 2, 0 );
  layout->addWidget( m_fitSumary, 3, 0 );
  
  AuxWindow::addHelpInFooter( footer(), "multi-file-calibration-dialog" );
  
  m_cancel = new WPushButton( "Cancel", footer() );
  m_fit    = new WPushButton( "Fit", footer() );
  m_use    = new WPushButton( "Use", footer());
  
  m_use->disable();
  m_cancel->clicked().connect( boost::bind( &Recalibrator::MultiFileCalibFit::handleFinish, this, WDialog::Rejected ) );
  m_use->clicked().connect( boost::bind( &Recalibrator::MultiFileCalibFit::handleFinish, this, WDialog::Accepted ) );
  m_fit->clicked().connect( this, &Recalibrator::MultiFileCalibFit::doFit );
  
  finished().connect( this, &Recalibrator::MultiFileCalibFit::handleFinish );
  
  m_fitSumary->disable();
  m_fitSumary->hide();
  
  const int w = 600 < viewer->renderedHeight() ? 600 : viewer->renderedHeight();
  const int h = static_cast<int>(0.8*viewer->renderedHeight());
  resizeWindow( w, h );
  
  updateCoefDisplay();
  rejectWhenEscapePressed();
  
  centerWindow();
}//MultiFileCalibFit constructor


Recalibrator::MultiFileCalibFit::~MultiFileCalibFit()
{
  
}//~MultiFileCalibFit()


void Recalibrator::MultiFileCalibFit::doFit()
{
  //XXX - The logic of this function is very similar to
  //      Recalibrator::recalibrateByPeaks() - so a refactoriztion should be
  //      done
  
  vector< std::shared_ptr<const PeakDef> > peakstouse;
  for( size_t i = 0; i < m_model->m_peaks.size(); ++i )
  {
     const vector< pair<bool,std::shared_ptr<const PeakDef> > > &peaks
                                                          = m_model->m_peaks[i];
     for( size_t j = 0; j < peaks.size(); ++j )
     {
       if( peaks[j].first && peaks[j].second )
         peakstouse.push_back( peaks[j].second );
     }//for( size_t j = 0; j < m_peaks.size(); ++j )
  }//for( size_t i = 0; i < m_peaks.size(); ++i )
  
  try
  {
    const size_t npeaks = peakstouse.size();
    const int num_coeff_fit = m_fitFor[0]->isChecked()
                              + m_fitFor[1]->isChecked()
                              + m_fitFor[2]->isChecked();
    
    if( num_coeff_fit < 1 )
    {
      const char *msg = "You must select at least one coefficient to fit for";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }//if( num_coeff_fit < 1 )
    
    if( num_coeff_fit > static_cast<int>(npeaks) )
    {
      const char *msg = "You must select at least as many peaks as coeficents to fit for";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }//if( num_coeff_fit < 1 )
    
    
    std::shared_ptr<const SpecMeas> meas = m_calibrator->m_interspec->measurment(SpectrumType::Foreground);
    
    if( !meas )
    {
      const char *msg = "You need to be displaying a foreground spectrum to do a calibration fit";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }
    
    
    std::shared_ptr<const Measurement> eqnmeas;
    const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
    
    for( size_t i = 0; !eqnmeas && i < meass.size(); ++i )
      if( meass[i]->num_gamma_channels() )
        eqnmeas = meass[i];
    if( !eqnmeas )
    {
      const char *msg = "You need to be displaying a foreground spectrum to do a calibration fit (unexpected error)";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }//if( !eqnmeas )
    
    const std::shared_ptr<const std::vector<float>> &binning = eqnmeas->channel_energies();
    const size_t nchannel = eqnmeas->num_gamma_channels();
    
    if( !binning || nchannel < 16 )
    {
      const char *msg = "The spectrum isnt high enough resolution to fit for a calibration";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }
    
    const SpecUtils::EnergyCalType calibration_type = eqnmeas->energy_calibration_model();
    
    vector<float> calib_coefs = eqnmeas->calibration_coeffs();
    
    switch( calibration_type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        break;
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
      {
        passMessage( "Invalid starting calibration type: unknown or lower bin "
                     "edge energy not allowed", "", WarningWidget::WarningMsgHigh );
        return;
        break;
      }//case LowerChannelEdge or InvalidEquationType
    }//switch( calibration_type )
    
    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contian the error of the fit mean for that peak
    vector<EnergyCal::RecalPeakInfo> peakInfos;
    
    for( size_t peakn = 0; peakn < npeaks; ++peakn )
    {
      const PeakDef &peak = *peakstouse[peakn];
      
      const double wantedEnergy = peak.gammaParticleEnergy();
      
      EnergyCal::RecalPeakInfo peakInfo;
      peakInfo.peakMean = peak.mean();
      peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
      if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
        peakInfo.peakMeanUncert = 0.5;
        
      peakInfo.photopeakEnergy = wantedEnergy;
      if( m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::FullRangeFraction )
      {
        peakInfo.peakMeanBinNumber = SpecUtils::find_fullrangefraction_channel( peak.mean(),
                                                                  calib_coefs, nchannel,
                                                                  eqnmeas->deviation_pairs(),
                                                                  0.001f );
      }else
      {
        const vector<float> fwfcoef = SpecUtils::polynomial_coef_to_fullrangefraction( calib_coefs, nchannel );
        peakInfo.peakMeanBinNumber = SpecUtils::find_fullrangefraction_channel( peak.mean(),
                                                                  fwfcoef, nchannel,
                                                                  eqnmeas->deviation_pairs(),
                                                                  0.001f );
      }//if( FullRangeFraction ) / else
        
      if( IsNan(peakInfo.peakMeanBinNumber)
          || IsInf(peakInfo.peakMeanBinNumber) )
        throw runtime_error( "Invalid result fromm "
                             "find_fullrangefraction_channel(...)" );
      peakInfos.push_back( peakInfo );
    }//for( int col = 0; col < numModelCol; ++col )
    
    
    bool fit_coefs = false;
    vector<double> parValues, parErrors;
    
    try
    { 
      vector<float> lls_fit_coefs = calib_coefs, lls_fit_coefs_uncert;
      
      //Convert eqation type to polynomial incase there are any fixed paramters.
      if( m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::FullRangeFraction )
        lls_fit_coefs = SpecUtils::fullrangefraction_coef_to_polynomial(calib_coefs, nchannel );
      lls_fit_coefs.resize( sm_numCoefs, 0.0f );
      
      vector<bool> fitfor( lls_fit_coefs.size(), false );
      
      for( size_t i = 0; i < calib_coefs.size() && i < sm_numCoefs; ++i )
        fitfor[i] = m_fitFor[i]->isChecked();
      
      const double chi2 = EnergyCal::fit_energy_cal_poly( peakInfos, fitfor,
                                              eqnmeas->deviation_pairs(),
                                              lls_fit_coefs, lls_fit_coefs_uncert );
      
      stringstream msg;
      msg << "\nfit_energy_cal_poly gave chi2=" << chi2 << " with coefs={";
      for( size_t i = 0; i < lls_fit_coefs.size(); ++i )
        msg << lls_fit_coefs[i] << "+-" << lls_fit_coefs_uncert[i] << ", ";
      msg << "}\n";
      cout << msg.str() << endl;
      
      parValues.resize( lls_fit_coefs.size() );
      parErrors.resize( lls_fit_coefs.size() );
      for( size_t i = 0; i < lls_fit_coefs.size(); ++i )
      {
        parValues[i] = lls_fit_coefs[i];
        parErrors[i] = lls_fit_coefs_uncert[i];
      }
      
      switch( m_calibrator->m_coeffEquationType )
      {
        case SpecUtils::EnergyCalType::Polynomial:
          break;
          
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          m_calibrator->m_coeffEquationType = SpecUtils::EnergyCalType::Polynomial;
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          lls_fit_coefs = SpecUtils::polynomial_coef_to_fullrangefraction(lls_fit_coefs, nchannel);
          lls_fit_coefs_uncert = SpecUtils::polynomial_coef_to_fullrangefraction(lls_fit_coefs_uncert, nchannel);
          break;
          
        case SpecUtils::EnergyCalType::InvalidEquationType:
        case SpecUtils::EnergyCalType::LowerChannelEdge:
          assert( 0 );  //shouldnt ever get here
          break;
      }//switch( m_coeffEquationType )
      
      fit_coefs = true;
    }catch( std::exception &e )
    {
      cerr << "fit_energy_cal_poly threw: " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      char buffer[512] = { '\0' };
      snprintf( buffer, sizeof(buffer)-1, "fit_energy_cal_poly threw: %s", e.what() );
      log_developer_error( __func__, buffer );
#endif
      fit_coefs = false;
    }//try / catch fit for coefficents using least linear squares
    
    
    if( !fit_coefs )
    {
      //This Minuit based fitting methodolgy is depreciated I think; the LLS code
      //  should work better, and seems to be releiabel, but leaving this code
      //  in for a while as a backup
      const auto eqntype = m_calibrator->m_coeffEquationType;
      const auto &devpairs = eqnmeas->deviation_pairs();
      if( calib_coefs.size() < sm_numCoefs )
        calib_coefs.resize( sm_numCoefs, 0.0 );
      
      vector<bool> fitfor( calib_coefs.size(), false );
      for( size_t i = 0; i < sm_numCoefs; ++i )
        fitfor[i] = m_fitFor[i]->isChecked();
        
      
      string warning_msg;
      vector<float> coefs, coefs_uncert;
      EnergyCal::fit_energy_cal_iterative( peakInfos, nchannel, eqntype, fitfor, calib_coefs,
                                 devpairs, coefs, coefs_uncert, warning_msg );
      
      if( warning_msg.size() )
        passMessage( warning_msg, "", WarningWidget::WarningMsgHigh );
      
      assert( coefs.size() == coefs_uncert.size() );
      parValues.resize( coefs.size(), 0.0 );
      parErrors.resize( coefs_uncert.size(), 0.0 );
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        parValues[i] = coefs[i];
        parErrors[i] = coefs_uncert[i];
      }
    }//if( !fit_coefs )
    
    for( size_t i = 0; i < parValues.size(); ++i )
      if( IsInf(parValues[i]) || IsNan(parValues[i]) )
        throw runtime_error( "Invalid calibration parameter from fit :(" );
    
    assert( parValues.size() >= sm_numCoefs );
    
    m_eqnType = calibration_type;
    for( int i = 0; i < sm_numCoefs; ++i )
    {
      m_calVal[i] = parValues[i];
      m_calUncert[i] = parErrors[i];
    }//for( size_t i = 0; i < parValues.size(); ++i )
    
    //Try to loop over peaks to give chi2 values and such
    vector<float> float_coef;
    for( const double d : parValues )
      float_coef.push_back( static_cast<float>(d) );
    
    if( m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::Polynomial
       || m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
      float_coef = SpecUtils::polynomial_coef_to_fullrangefraction( float_coef, nchannel );
    
    stringstream msg;
    for( const EnergyCal::RecalPeakInfo &info : peakInfos )
    {
      const double predictedMean = SpecUtils::fullrangefraction_energy( info.peakMeanBinNumber, float_coef, nchannel, eqnmeas->deviation_pairs() );
      double uncert = ((info.peakMeanUncert<=0.0) ? 1.0 : info.peakMeanUncert );
      double chi2 = pow(predictedMean - info.photopeakEnergy, 2.0 ) / (uncert*uncert);
      
      msg << "-Peak originally at " << info.peakMean << " +- "
          << info.peakMeanUncert << " keV for photopeak at "
          << info.photopeakEnergy << " keV ended up at " << predictedMean
          << " keV and contributed " << chi2 << " towards the chi2.\n";
    }//for( const &RecalPeakInfo info : peakInfo )
    
    m_fitSumary->setText( msg.str() );
    m_fitSumary->show();
    updateCoefDisplay();
    m_use->enable();
  }catch( std::exception &e )
  {
    m_use->disable();
    string exceptionmsg = e.what();
    string msg = "Failed calibration by fitting peak means.";
    
    if( SpecUtils::starts_with( exceptionmsg, "ErrorMsg" ) )
      msg = exceptionmsg.substr(8);
    
    cerr << "Recalibrator::MultiFileCalibFit::doFit(): \n\tCaught: " << exceptionmsg << endl;
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void doFit()


void Recalibrator::MultiFileCalibFit::updateCoefDisplay()
{
  for( size_t i = 0; i < sm_numCoefs; ++i )
  {
    char msg[32];
    snprintf( msg, sizeof(msg), "%.4g", m_calVal[i] );
    
    WString val(msg);
    if( m_calUncert[i] > 0.0 )
    {
#ifndef WT_NO_STD_WSTRING
      val += L" \x00B1 ";  //plusminus
#else
      val += " +- ";
#endif
      snprintf( msg, sizeof(msg), "%.4g", m_calUncert[i] );
      val += msg;
    }//if( uncertainty is available )
    
    m_coefvals[i]->setText( val );
  }//for( int i = 0; i < sm_numCoefs; ++i )
}//void updateCoefDisplay()


void Recalibrator::MultiFileCalibFit::handleFinish( WDialog::DialogCode result )
{
  switch( result )
  {
    case WDialog::Rejected:
      cerr << "\nRejected Recalibrator::MultiFileCalibFit" << endl;
    break;
      
    case WDialog::Accepted:
    {
      InterSpec *viewer = m_calibrator->m_interspec;
      
      // @TODO: the equation could totally be invalid
      vector<float> eqn( m_calVal, m_calVal + sm_numCoefs );
      while( !eqn.empty() && eqn.back()==0.0 )
        eqn.resize( eqn.size()-1 );
      
      static_assert( sm_numCoefs >= 3, "" );
      cerr << "\n\nm_calVal={" << m_calVal[0] << ", " << m_calVal[1] << ", " << m_calVal[2] << "}" << endl;
      
      vector<float> oldcalibpars;
      std::vector<std::pair<float,float>> devpairs;
      SpecUtils::EnergyCalType oldEqnType;
      
      std::shared_ptr<SpecMeas> fore = viewer->measurment(SpectrumType::Foreground);
      std::shared_ptr<SpecMeas> back = viewer->measurment(SpectrumType::Background);
      std::shared_ptr<SpecMeas> second = viewer->measurment(SpectrumType::SecondForeground);
      std::shared_ptr<const Measurement> displ_foreground
                                     = viewer->displayedHistogram(SpectrumType::Foreground);

      
      std::shared_ptr<const Measurement> eqnmeass[3];
      for( int i = 0; i < 3; ++i )
      {
        const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
        std::shared_ptr<SpecMeas> meas = viewer->measurment(type);
        if( meas )
        {
          const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
          for( size_t j = 0; !eqnmeass[i] && j < meass.size(); ++j )
            if( meass[j]->num_gamma_channels() )
              eqnmeass[i] = meass[j];
        }
      }//for( int i = 0; i < 3; ++i )
      
      if( fore && eqnmeass[ForeIndex] )
      {
        oldcalibpars = eqnmeass[ForeIndex]->calibration_coeffs();
        oldEqnType = eqnmeass[ForeIndex]->energy_calibration_model();
        devpairs = eqnmeass[ForeIndex]->deviation_pairs();
      }else if( second && eqnmeass[SecondIndex] )
      {
        oldcalibpars = eqnmeass[SecondIndex]->calibration_coeffs();
        oldEqnType = eqnmeass[SecondIndex]->energy_calibration_model();
        devpairs = eqnmeass[SecondIndex]->deviation_pairs();
      }else if( back && eqnmeass[BackIndex] )
      {
        oldcalibpars = eqnmeass[BackIndex]->calibration_coeffs();
        oldEqnType = eqnmeass[BackIndex]->energy_calibration_model();
        devpairs = eqnmeass[BackIndex]->deviation_pairs();
      }else
      {
        break;
      }
      
      vector<string> displayed_detectors = viewer->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
      
      
      if( !m_calibrator->m_applyTo[ForeIndex]->isChecked() )
        fore.reset();
      if( !m_calibrator->m_applyTo[BackIndex]->isChecked() )
        back.reset();
      if( !m_calibrator->m_applyTo[SecondIndex]->isChecked() )
        second.reset();
      
      //XXX - applying the calibration parameters does not shift peaks!
      for( size_t i = 0; i < m_model->m_peaks.size(); ++i )
      {
        const vector< pair<bool,std::shared_ptr<const PeakDef> > > &peaks
                                                          = m_model->m_peaks[i];
        bool used = false;
        for( size_t j = 0; j < peaks.size(); ++j )
          used |= peaks[j].first;
        
        if( !used )
          continue;
        
        SpectraFileModel *fileModel = viewer->fileManager()->model();
        std::shared_ptr<SpectraFileHeader> header
                                = fileModel->fileHeader( static_cast<int>(i) );
        if( !header )
          continue;
        
        std::shared_ptr<SpecMeas> meas = header->parseFile();
        if( !meas )
          continue;
        
        // Fix this for new calibration! //meas->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
        cerr << "\n\nRecalled " << header->displayName().toUTF8()
             << " to m_calVal={" << m_calVal[0] << ", " << m_calVal[1]
             << ", " << m_calVal[2] << "}" << endl;
        
        if( meas == fore )
        {
          fore.reset();
          Recalibrator::shiftPeaksForEnergyCalibration(
                                            m_calibrator->m_peakModel,
                                            eqn, devpairs, m_eqnType,
                                            meas, SpectrumType::Foreground,
                                            oldcalibpars, devpairs, oldEqnType );
        }else
        {
          assert( 0 );
          //meas->shiftPeaksForRecalibration( oldcalibpars, devpairs, oldEqnType,
          //                                 eqn, devpairs, m_eqnType );
        }//if( meas == fore ) / else
        
        if( meas == back )
          back.reset();
        if( meas == second )
          second.reset();
      }//for( size_t i = 0; i < m_peaks.size(); ++i )
      
      // Fix this for new calibration! //if( fore )
      // Fix this for new calibration! //  fore->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      // Fix this for new calibration! //if( back )
      // Fix this for new calibration! //  back->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      // Fix this for new calibration! //if( second )
      // Fix this for new calibration! //  second->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      
      viewer->refreshDisplayedCharts();
      m_calibrator->refreshRecalibrator();
      
      cerr << "\nAccepted Recalibrator::MultiFileCalibFit" << endl;
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  delete this;
}//void handleFinish(...)




