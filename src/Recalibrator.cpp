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
#include <iostream>
#include <sstream>
#include <vector>

//Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnUserParameterState.h"

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
#include <Wt/WCssDecorationStyle>

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/Recalibrator.h"
#include "InterSpec/WarningWidget.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/InterSpecApp.h"
#if ( USE_SPECTRUM_CHART_D3 )
#include "InterSpec/D3SpectrumDisplayDiv.h"
#else
#include "InterSpec/SpectrumDisplayDiv.h"
#endif
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/HelpSystem.h"

#include "SandiaDecay/SandiaDecay.h"

using namespace Wt;
using namespace std;


Recalibrator::Recalibrator(
#if ( USE_SPECTRUM_CHART_D3 )
                            D3SpectrumDisplayDiv *spectrumDisplayDiv,
#else
                            SpectrumDisplayDiv *spectrumDisplayDiv,
#endif
                            InterSpec *m_hostViewer,
                            PeakModel *peakModel)
  : WContainerWidget( 0 ),
    m_spectrumDisplayDiv( spectrumDisplayDiv ),
    m_hostViewer( m_hostViewer ),
    m_peakModel( peakModel ),
    m_isotopeView( 0 ),
    m_fitCoefButton( 0 ),
    m_graphicalRecal( 0 ),
    m_applyToLabel( nullptr ),
    m_applyTo{ nullptr },
    m_convertToPolynomialLabel( nullptr ),
    m_convertToPolynomial( nullptr ),
    m_convertToPolyDialog( nullptr ),
    m_revert( 0 ),
    m_acceptText( 0 ),
    m_acceptButtonDiv( 0 ),
    m_coeffEquationType( Measurement::InvalidEquationType ),
    m_devPairs( nullptr ),
    m_layout( nullptr )
{
  m_uncertainties[0] = m_uncertainties[1] = m_uncertainties[2] = -1.0;
  
  //A temporary layout is given to initialize the widgets
  WGridLayout * tempLayout = new WGridLayout();
  setRecalibratorLayout(tempLayout);

  Recalibrator::LayoutStyle style = kWide;
  try
  {
    if( !m_hostViewer->m_user->preferenceValue<bool>( "StartDocked" ) )
      style = kTall;
  } catch( std::exception & )
  {
  }
  
  
  initWidgets( style );
} // Recalibrator::Recalibrator(...)


Recalibrator::~Recalibrator()
{
  // no-op
} // Recalibrator::~Recalibrator()


void Recalibrator::setRecalibratorLayout(Wt::WGridLayout *layout)
{
  m_layout = layout;
} //  Recalibrator::setRecalibratorLayout(Wt::WGridLayout *layout)

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
void Recalibrator::initWidgets( Recalibrator::LayoutStyle style, AuxWindow* parent)
{
  m_currentLayout = style;
  
  const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_hostViewer );
  
  m_layout->clear();
  if( style == kTall )
  {
    m_layout->setContentsMargins(9,9,9,9);
  }
  else
  {
    m_layout->setContentsMargins(5,3,4,1);
  }
  
  
  // Set up the numeric boxes.
  for( int i = 0; i < 3; ++i )
  {
    m_coefficientDisplay[i] = new WDoubleSpinBox();
    
    // This needs to be 6 because Wt has some bad code--see WebUtils.C, 180-207, namely
    //  static const int exp[] = { 1, 10, 100, 1000, 10000, 100000, 1000000 };
    //  long long i = static_cast<long long>(d * exp[digits] + (d > 0 ? 0.49 : -0.49));
    // Note: this might be fixed at some point, possibly the next Wt release?
    m_coefficientDisplay[i]->setDecimals( 6 );
  }

  for( int i = 0; i < 3; ++i )
  {
    m_coefficientDisplay[i]->setValue( 0.0 );
    m_coefficientDisplay[i]->setRange( -1 * ( 2 << 20 ), 2 << 20 );

    // The spinboxes are disabled until a file is added~
    m_coefficientDisplay[i]->setDisabled( true );

    m_coefficientDisplay[i]->setMinimumSize( 95.0, WLength(1.0,WLength::FontEm) );
    
    // The buttons should all be enabled if the tick boxes are changed
    m_coefficientDisplay[i]->changed().connect( boost::bind( &Recalibrator::engageRecalibration, this, ApplyRecal ) );
  }//for( int i = 0; i < 3; ++i )

  // Set up the full layout
  WLabel *offsetLabel = new WLabel( "Offset Term " );
  m_exponentLabels[0] = new WLabel( "x10<sup>0</sup> keV" );

  WLabel *linearLabel = new WLabel( "Linear Term " );
  m_exponentLabels[1] = new WLabel( "x10<sup>0</sup> <sup>keV</sup>/<sub>chnl</sub>" );

  WLabel *quadraticLabel = new WLabel( "Quad. Term " );
  m_exponentLabels[2] = new WLabel( "x10<sup>0</sup> <sup>keV</sup>/<sub>chnl<sup>2</sup></sub>" );


  for( int i = 0; i < 3; ++i )
  {
    m_fitFor[i] = new WCheckBox( "fit" );
    m_fitFor[i]->setChecked( (i<2) );
    HelpSystem::attachToolTipOn( m_fitFor[i], "Fit for the value of this coefficient when using "
                            "peaks associated with isotopes, to determine "
                            "calibration." , showToolTipInstantly );
  }//for( int i = 0; i < 3; ++i )

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
    const SpectrumType type = SpectrumType(i);
    const char *name = 0;
    switch( type )
    {
      case kForeground:       name = "Foreground"; break;
      case kBackground:       name = "Background"; break;
      case kSecondForeground: name = "2nd Spec"; break;
    }
    
    m_applyTo[type] = new WCheckBox( name );
    m_applyTo[type]->setChecked( true );
    m_applyTo[type]->setDisabled(true);
    m_applyTo[type]->checked().connect( boost::bind(&Recalibrator::specTypeCheckedCallback, this, type, true ) );
    m_applyTo[type]->unChecked().connect( boost::bind(&Recalibrator::specTypeCheckedCallback, this, type, false ) );
  }//for( loop over SpectrumType )
  
  
  m_convertToPolynomialLabel = new WText( "Calibration is specified by channel energy" );
  m_convertToPolynomialLabel->setInline( false );
  m_convertToPolynomialLabel->addStyleClass( "RecalConvertToPolyTxt" );
  
  m_convertToPolynomial = new WPushButton( "Convert To Polynomial..." );
  m_convertToPolynomial->clicked().connect( this, &Recalibrator::startConvertToPolynomial );
  
  // Lastly, add in the buttons
  m_acceptButtonDiv = new WContainerWidget();
  
  if (parent)
    AuxWindow::addHelpInFooter( m_acceptButtonDiv, "energy-calibration-dialog" );
  
  WPushButton *multFiles = NULL; //need to reposition depending tall or wide mode
  
  if( style == kWide )
  {
    multFiles = new WPushButton("Mult. Files...", m_acceptButtonDiv );
  }else
  {
    multFiles = new WPushButton("Mult. Files...");
  } //kTall
  
  if (style==kTall && parent)
  {
      //add close button if it is in a wide state (AuxWindow)
      m_closeButton = parent->addCloseButtonToFooter("Close", false, m_acceptButtonDiv);
  } //kTall
  else
  {
      m_closeButton = new WPushButton( "Close", m_acceptButtonDiv);
  } //kWide
    
  m_closeButton->clicked().connect( m_hostViewer, &InterSpec::handRecalibratorWindowClose );
  m_closeButton->addStyleClass( "RecalCancelbtn" );
  
  m_revert = new WPushButton( "Revert", m_acceptButtonDiv );
  m_revert->setStyleClass( "RecalCancelbtn" );
  m_revert->setIcon( "InterSpec_resources/images/arrow_undo.png" );

  // And disable all the buttons
  m_revert->disable();

  //Add a Polynomial/FullWidthFraction/LowerEdge select here
  //Or add a button for when LowerEdge that converts to Polynomial
  
  //Hook the buttons up
  m_revert->clicked().connect( boost::bind( &Recalibrator::engageRecalibration, this, RevertRecal ) );
  
  int row = 0; //keeps track of row for adding widgets to layout
  m_layout->addWidget( offsetLabel,             row, 0 );
  m_layout->addWidget( m_coefficientDisplay[0], row, 1 );
  m_layout->addWidget( m_exponentLabels[0],     row, 2 );
  m_layout->addWidget( m_fitFor[0],             row, 3, AlignCenter );
  
  row++; //1
  m_layout->addWidget( linearLabel,             row, 0 );
  m_layout->addWidget( m_coefficientDisplay[1], row, 1 );
  m_layout->addWidget( m_exponentLabels[1],     row, 2 );
  m_layout->addWidget( m_fitFor[1],             row, 3, AlignCenter );

  row++; //2
  m_layout->addWidget( quadraticLabel,          row, 0 );
  m_layout->addWidget( m_coefficientDisplay[2], row, 1 );
  m_layout->addWidget( m_exponentLabels[2],     row, 2 );
  m_layout->addWidget( m_fitFor[2],             row, 3, AlignCenter );
  
  
  offsetLabel->setBuddy( m_coefficientDisplay[0] );
  linearLabel->setBuddy( m_coefficientDisplay[1] );
  quadraticLabel->setBuddy( m_coefficientDisplay[2] );
  
  
  m_devPairs = new DeviationPairDisplay( nullptr );
  m_devPairs->changed().connect( boost::bind( &Recalibrator::engageRecalibration, this, ApplyRecal ) );
  
  
  WPushButton *addButton = new WPushButton( "Add Pair" );
  addButton->setIcon( "InterSpec_resources/images/plus_min_white.png" );
  addButton->clicked().connect( boost::bind(&DeviationPairDisplay::newDevPair, m_devPairs, true) );
  
  if( style == kTall )
  {
    row++; //3
    m_layout->addWidget( m_devPairs, row, 0, 1, 4 );
    m_devPairs->setMinimumSize(WLength::Auto, WLength(150));
    m_devPairs->setMaximumSize(WLength::Auto, WLength(150));

    row++; //4
    m_layout->addWidget(addButton, row, 3, 1, 1 );
    m_devPairs->setStyleClass( "DeviationPairDisplayTall" );

    row++; //5
    m_layout->addWidget( m_isotopeView, row, 0, 1, 4 );
    m_isotopeView->setMinimumSize(WLength::Auto, WLength(150));
    m_layout->setRowStretch( row, 1 );
    m_layout->setColumnStretch( 1, 1 );

    row++; //6
    m_fitCoefButton->setStyleClass( "RecalFitBtn" );
    m_layout->addWidget( m_fitCoefButton, row, 3, 1, 1 );
    m_layout->addWidget( multFiles, row, 2, 1, 1 );
  }//kTall
  
  m_applyTo[kForeground]->setStyleClass( "RecalApplyToCbL" );
  m_applyTo[kSecondForeground]->setStyleClass( "RecalApplyToCbR" );
  m_applyTo[kBackground]->setStyleClass( "RecalApplyToCbL" );

  row++; //tall=7, wide=3
  m_layout->addWidget( m_applyToLabel,               row, 0 );
  m_layout->addWidget( m_applyTo[kForeground],       row, 1 );
  m_layout->addWidget( m_applyTo[kSecondForeground], row, 2 );

  row++; //tall=8, wide=4
  m_layout->addWidget( m_applyTo[kBackground], row, 1 );
  
  //For some reason if we put m_convertToPolynomialLabel and m_convertToPolynomial
  //  into their own rows and show/hide them, there are some wierd artificats
  //  introduced into the layout, so for now we'll just have things a bit.
  //row++;
  //m_layout->addWidget( m_convertToPolynomialLabel, row, 0, 1, 3 );
  //row++;
  //m_layout->addWidget( m_convertToPolynomial, row, 0, 1, 3, Wt::AlignCenter | Wt::AlignTop );
  
  WContainerWidget *polyDiv = new WContainerWidget();
  polyDiv->addWidget( m_convertToPolynomialLabel );
  polyDiv->addWidget( m_convertToPolynomial );
  m_convertToPolynomialLabel->hide();
  m_convertToPolynomial->hide();
  
  row++;
  m_layout->addWidget( polyDiv, row, 0, 1, 3 );
  
  
  multFiles->setIcon( "InterSpec_resources/images/page_white_stack.png" );
  if( style == kWide )
    multFiles->addStyleClass( "MultCalBtnWide" );
  
  HelpSystem::attachToolTipOn( multFiles,"Tool to use peaks from multiple files to fit for calibration", showToolTipInstantly, HelpSystem::Top );
  multFiles->clicked().connect( this, &Recalibrator::createMultifileFitter );
  
  if( style == kWide )
  {
    m_devPairs->setStyleClass( "DeviationPairDisplayWide" );
    
    for( int i = 0; i <= row; ++i )
      m_layout->setRowStretch( row, 1 );
    
    row++;
    m_layout->addWidget( m_acceptButtonDiv, row, 0, 1, 4 );
    
    m_layout->addWidget( m_isotopeView, 0, 4, row, 1 );
    m_layout->addWidget( m_devPairs, 0, 5, row, 1 );
    
    WContainerWidget* container = new WContainerWidget();
    m_fitCoefButton->addStyleClass( "RecalFitBtn" );
    container->addWidget(m_fitCoefButton);
    m_layout->addWidget( container, row, 4 , 1 , 1 );
    
    container = new WContainerWidget();
    container->addWidget(addButton);
    addButton->addStyleClass( "RecalFitBtn" );
    m_layout->addWidget( container, row, 5, 1, 1);

    m_layout->setColumnStretch( 4, 1 );
    
    //remove the floating
    multFiles->removeStyleClass( "FloatLeft" );
    m_closeButton->setHidden(true);
    m_revert->removeStyleClass( "RecalCancelbtn" );
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


void Recalibrator::setTallLayout(AuxWindow* parent)
{
  initWidgets( kTall, parent);
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
      = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_hostViewer );
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



void Recalibrator::specTypeCheckedCallback( const SpectrumType type,
                                            const bool checked )
{
  std::shared_ptr<SpecMeas> meas = m_hostViewer->measurment(type);
  for( int i = 0; i < 3; ++i )
  {
    const SpectrumType t = SpectrumType(i);
    
    if( t == i )
      continue;
    std::shared_ptr<SpecMeas> thismeas = m_hostViewer->measurment(t);
    if( thismeas == meas )
      m_applyTo[t]->setChecked( m_applyTo[type]->isChecked() );
  }//
}//void specTypeCheckedCallback(...)



WContainerWidget* Recalibrator::getFooter()
{
  return m_acceptButtonDiv;
}

void Recalibrator::shiftPeaksForEnergyCalibration( PeakModel *peakmodel,
                                                   const std::vector<float> &new_pars,
                                                   const std::vector< std::pair<float,float> > &new_devpairs,
                                                   Measurement::EquationType new_eqn_type,
                                                   std::shared_ptr<SpecMeas> meas,
                                                   const SpectrumType spectype,
                                                   std::vector<float> old_pars,
                                                   const std::vector< std::pair<float,float> > &old_devpairs,
                                                   Measurement::EquationType old_eqn_type )
{
  
  if( old_eqn_type==Measurement::EquationType::LowerChannelEdge
     || old_eqn_type==Measurement::EquationType::InvalidEquationType
     || new_eqn_type==Measurement::EquationType::LowerChannelEdge
     || new_eqn_type==Measurement::EquationType::InvalidEquationType )
  {
    //ToDo: We can only currently handle Polynomial or Full Range Fraction
    //      calibrations.  This could be fixed.
    return;
  }
    
  if( !meas )
    throw runtime_error( "shiftPeaksForEnergyCalibration: invalid input" );
  
  if( spectype != kForeground )
  {
    meas->shiftPeaksForRecalibration( old_pars, old_devpairs, old_eqn_type,
                                      new_pars, new_devpairs, new_eqn_type );
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
  
  meas->shiftPeaksForRecalibration( old_pars, old_devpairs, old_eqn_type,
                                   new_pars, new_devpairs, new_eqn_type );
  
  
  peakmodel->setPeakFromSpecMeas( meas, samples );
}//shiftPeaksForEnergyCalibration(...)




void Recalibrator::handleGraphicalRecalRequest( double xstart, double xfinish )
{
  try
  {
    if( !m_hostViewer->displayedHistogram(kForeground) )
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
  
  auto foreground = m_hostViewer->measurment(kForeground);
  auto displ_foreground = m_hostViewer->displayedHistogram(kForeground);
  
  if( !foreground || !displ_foreground )
  {
    passMessage( "No foreground to convert to polynomial binning.", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  if( displ_foreground->energy_calibration_model() != Measurement::LowerChannelEdge )
  {
    passMessage( "Foreground calibration is not specified by lower channel - wont convert to polynomial.", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  const vector<bool> useGamma = m_hostViewer->detectors_to_display();
  const set<int> samples = m_hostViewer->displayedSamples(kForeground);
  ShrdConstFVecPtr dataBinning = m_hostViewer->getBinning( samples, useGamma, foreground );
  
  if( !dataBinning )
  {
    passMessage( "Issue getting channel energies of foreground.", "", WarningWidget::WarningMsgHigh );
    return;
  }

  
  m_convertToPolyDialog = new AuxWindow( "Confirm", WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneModal) | DisableCollapse | SetCloseable | IsAlwaysModal );
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
  
  std::shared_ptr<SpecMeas> foreground = m_hostViewer->measurment(kForeground);
  std::shared_ptr<const Measurement> displ_foreground = m_hostViewer->displayedHistogram(kForeground);
  
  if( !foreground || !displ_foreground
      || displ_foreground->num_gamma_channels() < 16
      || !displ_foreground->channel_energies()
      || displ_foreground->channel_energies()->empty() )
  {
    passMessage( "Recalibrator::finishConvertToPolynomial: invalid foreground.", "", WarningWidget::WarningMsgHigh );
    return;
  }
  

  vector<float> coefs( 2, 0.0f );
  DeviationPairVec devpairs;
  coefs[1] = displ_foreground->channel_energies()->back() / displ_foreground->num_gamma_channels();
  
  //ToDo: If any of of the SpecMeas objects have Measurements with a differing
  //      number of channels, we will get an exception.  Should handle this
  //      correctly...
  if( foreground )
    foreground->recalibrate_by_eqn( coefs, devpairs, Measurement::EquationType::Polynomial );
  
  std::shared_ptr<SpecMeas> background = m_hostViewer->measurment(kBackground);
  if( background )
    background->recalibrate_by_eqn( coefs, devpairs, Measurement::EquationType::Polynomial );
  
  std::shared_ptr<SpecMeas> secondary = m_hostViewer->measurment(kSecondForeground);
  if( secondary )
    secondary->recalibrate_by_eqn( coefs, devpairs, Measurement::EquationType::Polynomial );
  
  if( foreground )
    m_hostViewer->displayForegroundData( true );
  if( secondary )
    m_hostViewer->displaySecondForegroundData();
  if( background )
    m_hostViewer->displayBackgroundData();
  
  refreshRecalibrator();
}//void finishConvertToPolynomial()


void Recalibrator::engageRecalibration( RecalAction action )
{
  std::shared_ptr<SpecMeas> foreground, back, second;
  foreground = m_hostViewer->measurment(kForeground);
  back       = m_hostViewer->measurment(kBackground);
  second     = m_hostViewer->measurment(kSecondForeground);

  std::shared_ptr<const Measurement> disp_foreground, disp_back, disp_second;
  disp_foreground = m_hostViewer->displayedHistogram( kForeground );
  disp_back = m_hostViewer->displayedHistogram( kBackground );
  disp_second = m_hostViewer->displayedHistogram( kSecondForeground );
  
  if( !foreground || !foreground->num_gamma_channels()
     || !disp_foreground || !disp_foreground->num_gamma_channels() )
    return;
 
#if( PERFORM_DEVELOPER_CHECKS )
  const double origForSum    = ((!!foreground) ? foreground->deep_gamma_count_sum() : 0.0);
  const double origBackSum   = ((!!back) ? back->deep_gamma_count_sum() : 0.0);
  const double origSecondSum = ((!!second) ? second->deep_gamma_count_sum() : 0.0);
#endif

  
  const DeviationPairVec devpairs = m_devPairs->deviationPairs();
  const DeviationPairVec &olddevpairs = disp_foreground->deviation_pairs();
  
  const vector< float > originalPars = disp_foreground->calibration_coeffs();
  const Measurement::EquationType old_calib_type = disp_foreground->energy_calibration_model();
  
  if( action == ApplyRecal )
  {
    m_revert->setDisabled( false );
    
    m_applyTo[kForeground]->setDisabled( true );
    m_applyTo[kBackground]->setDisabled( true );
    m_applyTo[kSecondForeground]->setDisabled( true );
  }else
  {
    m_revert->setDisabled( true );
  }

  // If they just want to RevertRecal, pass the original data back in, and then just refresh.
  if( action == RevertRecal )
  {
    passMessage( "Calibration reverted.", "", 0 );

    if( !!second
        && m_applyTo[kSecondForeground]->isChecked()
        && second!=foreground && second!=back )
    {
      const CalibrationInformation &info = m_originalCal[kSecondForeground];
      second->shiftPeaksForRecalibration( disp_second->calibration_coeffs(),
                                          disp_second->deviation_pairs(),
                                          disp_second->energy_calibration_model(),
                                          info.coefficients,
                                          info.deviationpairs,
                                          info.type );
      second->recalibrate_by_eqn( info.coefficients, info.deviationpairs, info.type );
    }
    
    if( !!back
        && m_applyTo[kBackground]->isChecked() && back!=foreground )
    {
      const CalibrationInformation &info = m_originalCal[kBackground];
      back->shiftPeaksForRecalibration( disp_back->calibration_coeffs(),
                                        disp_back->deviation_pairs(),
                                        disp_back->energy_calibration_model(),
                                        info.coefficients,
                                        info.deviationpairs,
                                        info.type );
      back->recalibrate_by_eqn( info.coefficients, info.deviationpairs, info.type );
    }

    //Have to treat the foreground measurment a bit special since the PeakModel
    //  is connected ot its peaks values
    if( !!foreground && m_applyTo[kForeground]->isChecked() )
    {
      const CalibrationInformation &info = m_originalCal[kForeground];
      foreground->recalibrate_by_eqn( info.coefficients, info.deviationpairs, info.type );
      shiftPeaksForEnergyCalibration( m_peakModel,
                                      info.coefficients, info.deviationpairs, info.type,
                                     foreground, kForeground,
                                     originalPars, olddevpairs, old_calib_type );
    }

    // Redraw everything
    if( foreground )
      m_hostViewer->displayForegroundData( true );
    if( second )
      m_hostViewer->displaySecondForegroundData();
    if( back )
      m_hostViewer->displayBackgroundData();

    refreshRecalibrator();

    return;
  }//if( action == RevertRecal )

  
  bool displayAll = true;
  for( bool b : m_hostViewer->detectors_to_display() )
    displayAll = (displayAll && b);
  if( !displayAll )
    passMessage( "Calibration only applied to displayed detectors", "", 1 );

  // Extract the values, noting the exponents.
  vector< float > eqn( 3 );
  eqn[0] = static_cast< float >( pow( 10.0, m_coeffExps[0] ) * m_coefficientDisplay[0]->value() );
  eqn[1] = static_cast< float >( pow( 10.0, m_coeffExps[1] ) * m_coefficientDisplay[1]->value() );
  eqn[2] = static_cast< float >( pow( 10.0, m_coeffExps[2] ) * m_coefficientDisplay[2]->value() );


  try  //begin check validity of equation {note: not a rigourous check}
  {
    const ShrdConstFVecPtr &binning = disp_foreground->channel_energies();
    const int nbin = static_cast<int>( binning->size() );
    vector<float> poly = eqn;
    
    vector<float> frfcoef = eqn;
    if( m_coeffEquationType == Measurement::Polynomial
       || m_coeffEquationType == Measurement::UnspecifiedUsingDefaultPolynomial )
      frfcoef = polynomial_coef_to_fullrangefraction( frfcoef, nbin );
    
    bool valid = checkFullRangeFractionCoeffsValid( frfcoef, devpairs, nbin );
    
    if( !valid )
    {
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
  
  vector<string> displayed_detectors = m_hostViewer->displayed_detector_names();
  
  switch( m_coeffEquationType )
  {
    case Measurement::Polynomial:
    case Measurement::UnspecifiedUsingDefaultPolynomial:
    case Measurement::FullRangeFraction:
    {
      //Recalibrate the foreground
      if( !!foreground && m_applyTo[kForeground]->isChecked() )
      {
        const vector<float> origeqn = disp_foreground->calibration_coeffs();
        const vector< pair<float,float> > origdevpars = disp_foreground->deviation_pairs();
        const Measurement::EquationType origtype = disp_foreground->energy_calibration_model();
        
        foreground->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType,
                                      displayed_detectors, false );
        
        if( action == ApplyRecal )
          shiftPeaksForEnergyCalibration( m_peakModel,
                                          eqn, devpairs, m_coeffEquationType,
                                          foreground, kForeground,
                                          origeqn, origdevpars, origtype );
      }//if( !!foreground && m_applyTo[kForeground]->isChecked() )


      // Recalibrate for the secondary data measurement
      if( !!second && m_applyTo[kSecondForeground]->isChecked() )
      {
        const vector<string> &foredet = foreground->detector_names();
        const vector<string> &backdet = second->detector_names();
        const set<string> forenameset( foredet.begin(), foredet.end() );
        const set<string> backnameset( backdet.begin(), backdet.end() );
        
        if( second != foreground || !m_applyTo[kForeground]->isChecked() )
          second->shiftPeaksForRecalibration( disp_second->calibration_coeffs(),
                                              disp_second->deviation_pairs(),
                                              disp_second->energy_calibration_model(),
                                              eqn, devpairs, m_coeffEquationType );

        if( forenameset == backnameset )
        {
          second->recalibrate_by_eqn( eqn, devpairs,
                                    m_coeffEquationType,
                                    displayed_detectors, false );
        }else
        {
          second->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType );
        }
      }//if( !!second && m_applyTo[kSecondForeground]->isChecked() )
      
      
      // Recalibrate for the background
      if( !!back && m_applyTo[kBackground]->isChecked() )
      {
        const vector<string> &foredet = foreground->detector_names();
        const vector<string> &backdet = back->detector_names();
        const set<string> forenameset( foredet.begin(), foredet.end() );
        const set<string> backnameset( backdet.begin(), backdet.end() );

        if( (back != foreground || !m_applyTo[kForeground]->isChecked())
            && (back != second || !m_applyTo[kSecondForeground]->isChecked()) )
          back->shiftPeaksForRecalibration( disp_back->calibration_coeffs(),
                                            disp_back->deviation_pairs(),
                                            disp_back->energy_calibration_model(),
                                            eqn, devpairs, m_coeffEquationType );
        if( forenameset == backnameset )
        {
          back->recalibrate_by_eqn( eqn, devpairs,
                                  m_coeffEquationType,
                                  displayed_detectors, false  );
        }else
        {
          back->recalibrate_by_eqn( eqn, devpairs, m_coeffEquationType );
        }//if( forenameset == backnameset )
        
      }//if( !!back && m_applyTo[kBackground]->isChecked() )
      
      
      const bool keepCurrentEnergyRange = true;
      m_hostViewer->displayForegroundData( keepCurrentEnergyRange );
      m_hostViewer->displaySecondForegroundData();
      m_hostViewer->displayBackgroundData();
 

      break;
    }//case Polynomial, FullRangeFraction:


    case Measurement::LowerChannelEdge:
    case Measurement::InvalidEquationType:
    {
      passMessage( "Recalibration can only be done for polynomial or full range"
                   " fraction calibrated spectrums", "engageRecalibration",
                   WarningWidget::WarningMsgHigh );
      break;
    }//case LowerChannelEdge and InvalidEquationType
  }//switch( m_coeffEquationType )


  // Recalibrate for the continuum?
//  m_hostViewer->displayTimeSeriesData( false ); // This can be removed, right?
  
  
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
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }
  
  if( fabs(origBackSum-newBackSum) > 1.0E-6*std::max(fabs(origBackSum),fabs(newBackSum)) )
  {
    snprintf( buffer, sizeof(buffer),
             "Gamma sum of background changed after calibration, "
             "pre=%1.8e, post=%1.8e", origBackSum, newBackSum );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }
  
  if( fabs(origSecondSum-newSecondSum) > 1.0E-6*std::max(fabs(origSecondSum),fabs(newSecondSum)) )
  {
    snprintf( buffer, sizeof(buffer),
             "Gamma sum of second foreground changed after calibration, "
             "pre=%1.8e, post=%1.8e", origSecondSum, newSecondSum );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }
#endif
}//void engageRecalibration( RECAL_ACTION action )


void Recalibrator::recalibrateByPeaks()
{
  try
  {
    const size_t npeaks = m_peakModel->npeaks();


    const int num_coeff_fit = m_fitFor[0]->isChecked()
                              + m_fitFor[1]->isChecked()
                              + m_fitFor[2]->isChecked();

    if( num_coeff_fit < 1 )
    {
      const char *msg = "You must select at least one coefficient to fit for";
      passMessage( msg, "recalibrateByPeaks", WarningWidget::WarningMsgHigh );
      return;
    }//if( num_coeff_fit < 1 )

    std::shared_ptr<const Measurement> data = m_spectrumDisplayDiv->histUsedForXAxis();
    std::shared_ptr<const SpecMeas> meas;
    std::shared_ptr<const Measurement> displ_meas;
    
    
    for( int i = 0; (!meas || !displ_meas) && (i < 3); ++i )
    {
      if( m_applyTo[i]->isChecked() )
      {
        const SpectrumType type = SpectrumType(i);
        meas = m_hostViewer->measurment(type);
        displ_meas = m_hostViewer->displayedHistogram(type);
      }//if( m_applyTo[i]->isChecked() )
    }//for( int i = 0; (!meas || !displ_meas) && (i < 3); ++i )
    
    if( !meas || !displ_meas || !data )
      throw std::runtime_error( "ErrorMsgNo spectrum available for calibration" );

    const ShrdConstFVecPtr &binning = displ_meas->channel_energies();
    
    if( !binning || binning->size() < 16 )
      throw std::runtime_error( "ErrorMsgNot enough binns for rebinning" );

    const size_t nbin = binning->size();
    const Measurement::EquationType calibration_type = displ_meas->energy_calibration_model();

    vector<float> calib_coefs = displ_meas->calibration_coeffs();

    switch( calibration_type )
    {
      case Measurement::Polynomial:
      case Measurement::UnspecifiedUsingDefaultPolynomial:
//        calib_coefs = polynomial_coef_to_fullrangefraction( calib_coefs, nbin );
      break;

      case Measurement::FullRangeFraction:
//        calib_coefs = fullrangefraction_coef_to_polynomial( calib_coefs );
      break;

      case Measurement::LowerChannelEdge:
      case Measurement::InvalidEquationType:
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

        std::shared_ptr<SpecMeas> &foreground = m_hostViewer->m_dataMeasurement;
        std::shared_ptr<const Measurement> displ_foreground
                                = m_hostViewer->displayedHistogram(kForeground);
        if( !!foreground && !!displ_foreground )
        {
          std::shared_ptr<SpecMeas> &second = m_hostViewer->m_secondDataMeasurement;
          std::shared_ptr<SpecMeas> &back = m_hostViewer->m_backgroundMeasurement;
          const DeviationPairVec &devpairs = displ_foreground->deviation_pairs();
        
          if( m_applyTo[kForeground]->isChecked() )
            foreground->rebin_by_eqn( poly_eqn, devpairs, Measurement::Polynomial );

          if( !!second && m_applyTo[kSecondForeground]->isChecked() )
            second->rebin_by_eqn( poly_eqn, devpairs, Measurement::Polynomial );

          if( !!back && m_applyTo[kBackground]->isChecked() )
            back->rebin_by_eqn( poly_eqn, devpairs, Measurement::Polynomial );

          m_coeffEquationType = Measurement::FullRangeFraction;
          calib_coefs = polynomial_coef_to_fullrangefraction( poly_eqn, nbin );
        }//if( !!foreground )
        
        break;
      }//case LowerChannelEdge and InvalidEquationType
    }//switch( calibration_type )

    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contian the error of the fit mean for that peak
    vector<RecalPeakInfo> peakInfos;

    for( size_t peakn = 0; peakn < npeaks; ++peakn )
    {
      const PeakDef &peak = m_peakModel->peak( peakn );

      if( peak.useForCalibration() )
      {
        const double wantedEnergy = peak.gammaParticleEnergy();

        RecalPeakInfo peakInfo;
        peakInfo.peakMean = peak.mean();
        peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
        if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
          peakInfo.peakMeanUncert = 0.5;
        
        peakInfo.photopeakEnergy = wantedEnergy;
        if( m_coeffEquationType == Measurement::FullRangeFraction )
        {
          peakInfo.peakMeanBinNumber = find_bin_fullrangefraction( peak.mean(),
                                           calib_coefs, nbin,
                                           displ_meas->deviation_pairs(),
                                           0.001f );
        }else
        {
          const vector<float> fwfcoef = polynomial_coef_to_fullrangefraction( calib_coefs, nbin );
          peakInfo.peakMeanBinNumber = find_bin_fullrangefraction( peak.mean(),
                                                                  fwfcoef, nbin,
                                                                  displ_meas->deviation_pairs(),
                                                                  0.001f );
        }

        if( IsNan(peakInfo.peakMeanBinNumber)
            || IsInf(peakInfo.peakMeanBinNumber) )
          throw runtime_error( "Invalid result fromm "
                               "find_bin_fullrangefraction(...)" );

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

    PolyCalibCoefMinFcn chi2Fcn( peakInfos,
                                 binning->size(),
                                 m_coeffEquationType,
                                 displ_meas->deviation_pairs() );
    ROOT::Minuit2::MnUserParameters inputPrams;

    if( calib_coefs.size() < 3 )
      calib_coefs.resize( 3, 0.0 );

//    for( int i = 0; i < calib_coefs.size(); ++i )
//      cerr << "initial calib_coefs[" << i << "]=" << calib_coefs[i] << endl;

    string name = "0";
    double delta = 1.0/pow( static_cast<double>(nbin), 0.0 );
    if( m_coeffEquationType == Measurement::FullRangeFraction )  //set delats to 0.5 keV for middle bin
      delta = 0.5;

    if( m_fitFor[0]->isChecked() )
      inputPrams.Add( name, calib_coefs[0], delta );
    else
      inputPrams.Add( name, calib_coefs[0] );

    name = "1";
    delta = 1.0/pow( static_cast<double>(data->GetNbinsX()), 1.0 );
    if( m_coeffEquationType == Measurement::FullRangeFraction )
      delta = 1.0;

    if( m_fitFor[1]->isChecked() )
    {
      inputPrams.Add( name, calib_coefs[1], delta );
      inputPrams.SetLowerLimit( name, 0.0 );
    }else
    {
      inputPrams.Add( name, calib_coefs[1] );
    }

    name = "2";
    delta = 1.0/pow( static_cast<double>(data->GetNbinsX()), 2.0 );
    if( m_coeffEquationType == Measurement::FullRangeFraction )
      delta = 0.1;

    if( m_fitFor[2]->isChecked() )
      inputPrams.Add( name, calib_coefs[2], delta );
    else
      inputPrams.Add( name, calib_coefs[2] );

    ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
    ROOT::Minuit2::MnStrategy strategy( 1 ); //0 low, 1 medium, >=2 high

    const int npars = inputPrams.VariableParameters();
    const unsigned int maxFcnCall = 200 + 100*npars + 5*npars*npars;
    const double tolerance = 0.1; //0.05 * peakInfos.size();

    ROOT::Minuit2::CombinedMinimizer fitter;
    ROOT::Minuit2::FunctionMinimum minimum
                           = fitter.Minimize( chi2Fcn, inputParamState,
                                              strategy, maxFcnCall, tolerance );
    
    //Not sure why Minuit2 doesnt like converging on the minumum verry well, but
    //  rather than showing the user an error message, we'll give it anither try
    if( minimum.IsAboveMaxEdm() )
    {
      ROOT::Minuit2::MnMigrad fitter( chi2Fcn, inputParamState, strategy );
      minimum = fitter( maxFcnCall, tolerance );
    }//if( minimum.IsAboveMaxEdm() )
    
    if( !minimum.IsValid() )
    {
      if( !minimum.HasValidCovariance() || !minimum.HasValidParameters() )
      {
        string msg = "Fit for calibration parameters failed.";
        if( m_fitFor[2] )
          msg += " you might try not fitting for quadratic term.";
        passMessage( msg, "", WarningWidget::WarningMsgHigh );
        return;
      }
      
      cerr << minimum << endl;
      stringstream msg;
      msg << "Warning: calibration coefficient fit results "
      "may be invalid, please check, and if necessary "
      "revert this calibration instead.";
      
      if( minimum.IsAboveMaxEdm() )
      {
        msg << " The estimated distance to optimal calibration parameters is "
               "too large.";
        
        cerr << "EDM=" << minimum.Edm() << ", tolerance=" << tolerance << ", npeaks=" << peakInfos.size() << endl;
      }
      if( minimum.HasReachedCallLimit() )
        msg << " To many calls to chi2 routine were made.";
      
      
      passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    }//if( !minimum.IsValid() )

    const ROOT::Minuit2::MnUserParameters &fitPrams = minimum.UserState().Parameters();
    const vector<double> parValues = fitPrams.Params();
    const vector<double> parErrors = fitPrams.Errors();

//    for( int i = 0; i < parValues.size(); ++i )
//      cerr << "parValues[" << i << "]=" << parValues[1] << endl;

    for( size_t i = 0; i < parValues.size(); ++i )
      if( IsInf(parValues[i]) || IsNan(parValues[i]) )
        throw runtime_error( "Invalid calibration parameter from fit :(" );
        
    assert( parValues.size() == 3 );
    
    for( int i = 0; i < 3; ++i )
    {
      m_coeffExps[i] = 0.0;
      m_coeffExps[i] = static_cast<int>( floor( log10( abs( parValues[i] ) ) ) );
      m_coefficientDisplay[i]->setValue( parValues[i] * pow( 10.0, -m_coeffExps[i] ) );
      m_uncertainties[i] = parErrors[i];
    }//for( size_t i = 0; i < parValues.size(); ++i )

    engageRecalibration( ApplyRecal );
  }catch( std::exception &e )
  {
    string exceptionmsg = e.what();
    string msg = "Failed calibration by fitting peak means.";
    
    if( UtilityFunctions::starts_with( exceptionmsg, "ErrorMsg" ) )
      msg = exceptionmsg.substr(8);
    
    cerr << SRC_LOCATION << "\n\tCaught: " << exceptionmsg << endl;
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//try / catch

}//void recalibrateByPeaks()


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
  foreground = m_hostViewer->measurment(kForeground);
  background = m_hostViewer->measurment(kBackground);
  second     = m_hostViewer->measurment(kSecondForeground);
  
  std::shared_ptr<const Measurement> displ_foreground, displ_background, displ_second;
  displ_foreground = m_hostViewer->displayedHistogram(kForeground);
  displ_background = m_hostViewer->displayedHistogram(kBackground);
  displ_second     = m_hostViewer->displayedHistogram(kSecondForeground);
  
  const vector< bool > useGamma = m_hostViewer->detectors_to_display();
  const set<int> samples = m_hostViewer->displayedSamples(kForeground);
  ShrdConstFVecPtr dataBinning = m_hostViewer->getBinning( samples, useGamma, foreground );

  if( m_convertToPolyDialog )
  {
    AuxWindow::deleteAuxWindow( m_convertToPolyDialog );
    m_convertToPolyDialog = nullptr;
  }
  
  // Limiting the data binning size protects against rebinning low channel count
  // detectors that it doesnt make sense to recalibrate, because you wouldnt be
  // able to see peaks anyway
  if( !foreground || !displ_foreground || !dataBinning
      || dataBinning->size() < 16 || (!!foreground
        && (displ_foreground->energy_calibration_model()==Measurement::LowerChannelEdge
            || displ_foreground->energy_calibration_model()==Measurement::InvalidEquationType)) )
  {
    m_fitCoefButton->disable();
    m_applyToLabel->hide();
    for( int i = 0; i < 3; ++i )
    {
      m_coefficientDisplay[i]->setValue( 0.0 );
      m_fitFor[i]->hide();
      m_applyTo[i]->hide();
      m_coefficientDisplay[i]->hide();
      m_exponentLabels[i]->hide();
    }
    
    if( dataBinning && foreground
        && dataBinning->size() >= 16
        && displ_foreground->energy_calibration_model()==Measurement::LowerChannelEdge )
    {
      m_convertToPolynomialLabel->show();
      m_convertToPolynomial->show();
    }else
    {
      m_convertToPolynomialLabel->hide();
      m_convertToPolynomial->hide();
    }
    
    return;
  }//if( we cant recalibrate this detector anyway )
  
  if( m_coefficientDisplay[0]->isHidden() )
  {
    m_convertToPolynomialLabel->hide();
    m_convertToPolynomial->hide();
    m_fitCoefButton->enable();
    m_applyToLabel->show();
    for( int i = 0; i < 3; ++i )
    {
      m_fitFor[i]->show();
      m_applyTo[i]->show();
      m_coefficientDisplay[i]->show();
      m_exponentLabels[i]->show();
    }
  }//if( m_coefficientDisplay[0]->isHidden() )
  
  m_applyTo[kForeground]->enable();
  m_applyTo[kBackground]->enable();
  m_applyTo[kSecondForeground]->enable();

  const size_t nbin = dataBinning->size();
  const vector< float > &binning = *( dataBinning );

  // Store the original information, then create a copy and operate on the copy. Note objects pointed to
  // in other places by m_dataMeasurement etc. may be used by other locations in the code.
  // The backup is made by copy instead of via pointer shuffling.

//  m_originalPeakDefs.clear();
//  const std::deque< PeakShrdPtr > &peaks = m_peakModel->peaks();
//  for( const PeakShrdPtr &p : peaks )
//    m_originalPeakDefs.push_back( *p );

  for( int i = 0; i < 3; ++i )
  {
    const SpectrumType type = SpectrumType(i);
    
    std::shared_ptr<SpecMeas> meas = m_hostViewer->measurment(type);
    std::shared_ptr<const Measurement> displ_meas
                                       = m_hostViewer->displayedHistogram(type);
    if( !displ_meas )
    {
      m_originalCal[i].reset();
    }else
    {
      m_originalCal[i].type = displ_meas->energy_calibration_model();
      m_originalCal[i].coefficients = displ_meas->calibration_coeffs();
      m_originalCal[i].deviationpairs = displ_meas->deviation_pairs();
      
      m_originalCal[i].sample_numbers = m_hostViewer->displayedSamples( type );
      m_originalCal[i].detectors_numbers
                                    = m_hostViewer->displayedDetectorNumbers();
    }
  }//for( loop over SpectrumType );
  
  
  vector< float > equationCoefficients = displ_foreground->calibration_coeffs();
  m_coeffEquationType = m_originalCal[kForeground].type;
  if( m_coeffEquationType == Measurement::UnspecifiedUsingDefaultPolynomial )
    m_coeffEquationType = Measurement::Polynomial;
  
  switch( m_coeffEquationType )
  {
    // Note: fall-through intentional
  case Measurement::Polynomial:
  case Measurement::UnspecifiedUsingDefaultPolynomial:
  case Measurement::FullRangeFraction:
    break;

  case Measurement::LowerChannelEdge:
    //Note: we will never make it here!
  case Measurement::InvalidEquationType:
    equationCoefficients.resize( 3 ); // Third degree polynomial *wonk*
    if( binning.size() != 0 )
    {
      equationCoefficients[0] = binning[0];
      equationCoefficients[1] = ( binning[ binning.size() - 1 ] - binning[0] ) / nbin;
    }
    else
    {
      equationCoefficients[1] = 0.0;
    }
    equationCoefficients[2] = 0.0;

    break;
  } // switch( m_coeffEquationType )

  // Find how much the things should tick by as a base
  double tickLevel[3];
  if( m_coeffEquationType == Measurement::Polynomial
      || m_coeffEquationType == Measurement::UnspecifiedUsingDefaultPolynomial )
  {
    tickLevel[0] = 1.0;
    tickLevel[1] = 1.0 / nbin;
    tickLevel[2] = 1.0 / (nbin * nbin);
  }else
  {
    //FullRangeFraction
    tickLevel[0] = 1.0;
    tickLevel[1] = 2.0;
    tickLevel[2] = 4.0;
  }


  // Update the powers. If they're tiny (error-zone), treat them as zero.
  if ( equationCoefficients.size() < 1 || abs( equationCoefficients[0] ) < pow( 10.0, -20.0 ) ) {
    m_coeffExps[0] = 0; // log10( 1.0 ), log10( tickLevel[0] )
  } else {
    m_coeffExps[0] = static_cast<int>( floor( log10( abs( equationCoefficients[0] ) ) ) );
  }

  if ( equationCoefficients.size() < 2 || abs( equationCoefficients[1] ) < pow( 10.0, -20.0 ) ) {
    m_coeffExps[1] = static_cast<int>( floor( log10( tickLevel[1] ) ) );
  } else {
    m_coeffExps[1] = static_cast<int>( floor( log10( abs( equationCoefficients[1] ) ) ) );
  }

  if ( equationCoefficients.size() < 3 || abs( equationCoefficients[2] ) < pow( 10.0, -20.0 ) ) {
    m_coeffExps[2] = static_cast<int>( floor( log10( tickLevel[2] ) ) );
  } else {
    m_coeffExps[2] = static_cast<int>( floor( log10( abs( equationCoefficients[2] ) ) ) );
  }
  // And display them.
  stringstream offsetExp, linearExp, quadraticExp;
  offsetExp    << "x10<sup>" << m_coeffExps[0] << "</sup>";
  linearExp    << "x10<sup>" << m_coeffExps[1] << "</sup>";
  quadraticExp << "x10<sup>" << m_coeffExps[2] << "</sup>";

  switch( m_coeffEquationType )
  {
    case Measurement::Polynomial:
    case Measurement::UnspecifiedUsingDefaultPolynomial:
    case Measurement::InvalidEquationType:
      offsetExp    << " keV";
      linearExp    << " <sup>keV</sup>/<sub>chnl</sub>";
      quadraticExp << " <sup>keV</sup>/<sub>chnl<sup>2</sup></sub>";
    break;

    case Measurement::FullRangeFraction:
      offsetExp    << " keV";
      linearExp    << " <sup>keV</sup>/<sub>FWF</sub>";
      quadraticExp << " <sup>keV</sup>/<sub>FWF<sup>2</sup></sub>";
    break;

    case Measurement::LowerChannelEdge:
    break;
      offsetExp    << " keV";
      linearExp    << " keV";
      quadraticExp << " keV";
    break;
  }//switch( m_coeffEquationType )


  m_exponentLabels[0]->setText(    offsetExp.str() );
  m_exponentLabels[1]->setText(    linearExp.str() );
  m_exponentLabels[2]->setText( quadraticExp.str() );

  // Set up the little tick/spin/whatever boxes
  //m_coefficientDisplay[0]->setSingleStep( 1.0 );
  //m_linear->setSingleStep( 1.0 / nbin );
  //m_quadratic->setSingleStep( 1.0 / ( nbin * nbin ) );
  for( int i = 0; i < 3; ++i )
  {
    m_coefficientDisplay[i]->setSingleStep( tickLevel[i] * pow( 10.0, -m_coeffExps[i] ) );
    m_coefficientDisplay[i]->setDisabled( false );
  }//for( int i = 0; i < 3; ++i )
  
  if( equationCoefficients.size() > 0 )
    m_coefficientDisplay[0]->setValue( equationCoefficients[0] * pow( 10.0, -m_coeffExps[0] ) );
  else
    m_coefficientDisplay[0]->setValue( 0.0 );
  
  if( equationCoefficients.size() > 1 )
    m_coefficientDisplay[1]->setValue( equationCoefficients[1] * pow( 10.0, -m_coeffExps[1] ) );
  else
    m_coefficientDisplay[1]->setValue( 0 );

  if( equationCoefficients.size() > 2 )
    m_coefficientDisplay[2]->setValue( equationCoefficients[2] * pow( 10.0, -m_coeffExps[2] ) );
  else
    m_coefficientDisplay[2]->setValue( 0 );
  
  m_devPairs->setDeviationPairs( m_originalCal[kForeground].deviationpairs );
  
  // Disable whichever check boxes correspond to null
  m_applyTo[kForeground]->setDisabled( !foreground );
  m_applyTo[kForeground]->setChecked( !!foreground);
  
  m_applyTo[kBackground]->setDisabled( !background );
  m_applyTo[kBackground]->setChecked( !!background);
  
  m_applyTo[kSecondForeground]->setDisabled( !second );
  m_applyTo[kSecondForeground]->setChecked( !!second);
}// void refreshRecalibrator( InterSpec *m_hostViewer )




DevPair::DevPair( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_energy( new WDoubleSpinBox() ),
    m_offset( new WDoubleSpinBox() ),
    m_delete( new WContainerWidget() )
{
  WGridLayout* layout = new WGridLayout();
  layout->setContentsMargins(0, 0, 0, 0);
  
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
  m_energy->setMinimum( -10000.0 );
  m_energy->setMaximum( 1000000.0 );
  m_offset->setMinimum( -10000.0 );
  m_offset->setMaximum( 1000000.0 );
  m_delete->addStyleClass( "Wt-icon DeleteDevPair" );
}//DevPair constructor


void DevPair::setDevPair( const std::pair<float,float> &d )
{
  m_energy->setValue( d.first );
  m_offset->setValue( d.second );
}


std::pair<float,float> DevPair::devPair() const
{
  return pair<float,float>( float(m_energy->value()), float(m_offset->value()) );
}


DeviationPairDisplay::DeviationPairDisplay( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_pairs( NULL )
{
  setStyleClass( "DeviationPairDisplay" );
  WLabel *title = new WLabel( "Deviation Pairs", this );
  title->setStyleClass( "DevPairTitle" );
  title->setInline( false );
  WContainerWidget *headerDiv = new WContainerWidget( this );
  headerDiv->setStyleClass( "DevPairHeader" );
  WLabel *header = new WLabel( "Energy (keV)", headerDiv );
  header->setStyleClass(  "DevPairEnergyHeader" );
  header = new WLabel( "Offset (keV)", headerDiv );
  header->setStyleClass(  "DevPairOffsetHeader" );
  
  m_pairs = new WContainerWidget( this );
  m_pairs->setStyleClass( "DevPairsContainer" ); 
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
}//setDeviationPairs(...)


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
      emitChanged();
      return;
    }
  }//for( WWidget *t : childs )
}//removeDevPair(...)

DevPair *DeviationPairDisplay::newDevPair( const bool emitChangedNow )
{
  DevPair *dev = new DevPair( m_pairs );
  dev->m_delete->clicked().connect( boost::bind( &DeviationPairDisplay::removeDevPair, this, dev ) );
  dev->m_energy->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_offset->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  if( emitChangedNow )
    emitChanged();
  return dev;
}//newDevPair()

Wt::Signal<> &DeviationPairDisplay::changed()
{
  return m_changed;
}



Recalibrator::PolyCalibCoefMinFcn::PolyCalibCoefMinFcn(
                     const vector<RecalPeakInfo> &peakInfo,
                     const size_t nbin,
                     Measurement::EquationType eqnType,
                     const std::vector< std::pair<float,float> > &devpair )
  : ROOT::Minuit2::FCNBase(),
    m_nbin( nbin ),
    m_eqnType( eqnType ),
    m_peakInfo( peakInfo ),
    m_devpair( devpair )
{
  switch( m_eqnType )
  {
    case Measurement::Polynomial:
    case Measurement::FullRangeFraction:
    case Measurement::UnspecifiedUsingDefaultPolynomial:
      break;
    case Measurement::LowerChannelEdge:
    case Measurement::InvalidEquationType:
      throw runtime_error( "PolyCalibCoefMinFcn can only wotrk with Full "
                          "Range Fraction and Polynomial binnings" );
      break;
  }
}//PolyCalibCoefMinFcn( constructor )


Recalibrator::PolyCalibCoefMinFcn::~PolyCalibCoefMinFcn()
{
  // no-op
}


double Recalibrator::PolyCalibCoefMinFcn::Up() const
{
  return 1.0;
}


double Recalibrator::PolyCalibCoefMinFcn::operator()( const vector<double> &coef ) const
{
  double chi2 = 0.0;

  vector<float> float_coef;
  for( const double d : coef )
  {
    if( IsInf(d) || IsNan(d) )
    {
      cerr << "Recalibrator::PolyCalibCoefMinFcn::operator(): "
           << "invalid input paramater" << endl;
      return 99999999.0;
    }
    float_coef.push_back( static_cast<float>(d) );
  }

//  if( m_eqnType == Measurement::FullRangeFraction )
//    float_coef = fullrangefraction_coef_to_polynomial( float_coef, m_nbin );
  if( m_eqnType == Measurement::Polynomial
      || m_eqnType == Measurement::UnspecifiedUsingDefaultPolynomial )
    float_coef = polynomial_coef_to_fullrangefraction( float_coef, m_nbin );

  
  for( const float d : float_coef )
  {
    if( IsInf(d) || IsNan(d) )
    {
      cerr << "Recalibrator::PolyCalibCoefMinFcn::operator(): invalid"
              " conversion to Polynomial from full width fraction" << endl;
       return 99999999.0;
    }
  }//for( float d : float_coef )
  
//  if( !MeasurementInfo::calibration_is_valid( m_eqnType, float_coef, m_devpair, m_nbin ) )
//  {
//    cerr << "Recalibrator::PolyCalibCoefMinFcn::operator(): "
//         << "invalid input paramater" << endl;
//    return 99999999.0;
//  }
  
  const float nearend = fullrangefraction_energy( m_nbin-2, float_coef, m_nbin, m_devpair );
  const float end = fullrangefraction_energy( m_nbin-1, float_coef, m_nbin, m_devpair );
  const float begin = fullrangefraction_energy( 0, float_coef, m_nbin, m_devpair );
  const float nearbegin = fullrangefraction_energy( 1, float_coef, m_nbin, m_devpair );
  
//  const float nearend = bin_number_to_energy_polynomial( m_nbin-2, float_coef, m_nbin );
//  const float end = bin_number_to_energy_polynomial( m_nbin-1, float_coef, m_nbin );
//  const float begin = bin_number_to_energy_polynomial( 0, float_coef, m_nbin );
//  const float nearbegin = bin_number_to_energy_polynomial( 1, float_coef, m_nbin );
  
  if( (nearend >= end) || (begin >= nearbegin) )
  {
    cerr << "Got nearend=" << nearend << ", end= " << end
         << ", begin=" << begin << ", nearbegin" << nearbegin << endl;
    return 99999999.0;
  }//if( almostLastEnergy > lastEnergy )
  
  
  for( const RecalPeakInfo &info : m_peakInfo )
  {
//    const double predictedMean = bin_number_to_energy_polynomial( info.peakMeanBinNumber, float_coef, m_nbin );
    const double predictedMean = fullrangefraction_energy( info.peakMeanBinNumber, float_coef, m_nbin, m_devpair );
    double uncert = ((info.peakMeanUncert<=0.0) ? 1.0 : info.peakMeanUncert );
    chi2 += pow(predictedMean - info.photopeakEnergy, 2.0 ) / (uncert*uncert);
  }//for( const &RecalPeakInfo info : peakInfo )
  
  if( IsInf(chi2) || IsNan(chi2) )
  {
    cerr << "Recalibrator::PolyCalibCoefMinFcn::operator(): invalid result chi2" << endl;
    return 1000.0;
  }
  
  return chi2;
}//operator()


void Recalibrator::CalibrationInformation::reset()
{
  type = Measurement::InvalidEquationType;
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
  
  if( !cal->m_hostViewer->displayedHistogram(kForeground) )
  {
    throw runtime_error( "You must have a dispayed spectrum to do a recalubration" );
    return;
  }//if( !spectrum->m_hostViewer->displayedHistogram(kForeground) )
  
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
  
  
#if ( USE_SPECTRUM_CHART_D3 )
  const D3SpectrumDisplayDiv *spectrum = m_calibrator->m_spectrumDisplayDiv;
#else
  const SpectrumDisplayDiv *spectrum = m_calibrator->m_spectrumDisplayDiv;
#endif
  std::shared_ptr<const Measurement> secondaryH  = spectrum->secondData();
  std::shared_ptr<const Measurement> backgroundH = spectrum->background();
  
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
    
    const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", cal->m_hostViewer );

    HelpSystem::attachToolTipOn( m_preserveLastCal,"This is only possible if a offset or"
                                " linear term adjustment was previously"
                                " made within the last 2 minutes.", showToolTipInstantly );
    m_preserveLastCal->setChecked();
    buttonBox->disable();
    m_preserveLastCal->checked().connect( buttonBox, &WWidget::disable );
    m_preserveLastCal->unChecked().connect( buttonBox, &WWidget::enable );
  }//if( preserve last cal possibly )
  
  bool displayAll = true;
  for( bool b : cal->m_hostViewer->detectors_to_display() )
    displayAll = (displayAll && b);
  
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
  return MeasurementInfo::calibration_is_valid( Measurement::FullRangeFraction,
                                              eqn, devpairs, nbin );

//  const float nearend   = fullrangefraction_energy( nbin-2, eqn, nbin, devpairs );
//  const float end       = fullrangefraction_energy( nbin-1, eqn, nbin, devpairs );
//  const float begin     = fullrangefraction_energy( 0,      eqn, nbin, devpairs );
//  const float nearbegin = fullrangefraction_energy( 1,      eqn, nbin, devpairs );
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
  InterSpec *viewer = m_calibrator->m_hostViewer;
  std::shared_ptr<SpecMeas> foreground = viewer->m_dataMeasurement;
  std::shared_ptr<const Measurement> displ_foreground
                                      = viewer->displayedHistogram(kForeground);
  if( !foreground || !displ_foreground )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  
  float shift = finalE - startE;
  vector<float> eqn = displ_foreground->calibration_coeffs();
  Measurement::EquationType eqnType = displ_foreground->energy_calibration_model();
  DeviationPairVec deviationPairs = displ_foreground->deviation_pairs();
  const DeviationPairVec origdev = deviationPairs;
  const size_t nbin = foreground->num_gamma_channels();
  
  switch( eqnType )
  {
    case Measurement::Polynomial:
    case Measurement::UnspecifiedUsingDefaultPolynomial:
      eqnType = Measurement::FullRangeFraction;
      eqn = polynomial_coef_to_fullrangefraction( eqn, nbin );
    break;
    case Measurement::FullRangeFraction:
    break;
    case Measurement::LowerChannelEdge:
    case Measurement::InvalidEquationType:
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
    const float lowbin = find_bin_fullrangefraction( m_lastEnergy, orig_eqn,
                                                  nbin, deviationPairs, 0.001f );
    const float upbin = find_bin_fullrangefraction( startE, orig_eqn,
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
        const float binnum = find_bin_fullrangefraction( startE, eqn, nbin,
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
          const float bin = find_bin_fullrangefraction( startE, orig_eqn, nbin,
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
  
//  if( eqnType == Measurement::Polynomial || eqnType == Measurement::UnspecifiedUsingDefaultPolynomial )
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
  
  const vector<string> detectors
                       = m_calibrator->m_hostViewer->displayed_detector_names();
  foreground->recalibrate_by_eqn( eqn, deviationPairs, eqnType, detectors, false );

  std::shared_ptr<SpecMeas> background = viewer->m_backgroundMeasurement;
  std::shared_ptr<SpecMeas> secondary = viewer->m_secondDataMeasurement;
  std::shared_ptr<const Measurement> displ_background
                                      = viewer->displayedHistogram(kBackground);
  std::shared_ptr<const Measurement> displ_secondary
                                = viewer->displayedHistogram(kSecondForeground);
  
  if( !m_foregroundOnly || !m_foregroundOnly->isChecked() )
  {    
    if( background && displ_background && (background != foreground) )
    {
      const vector<string> &foredet = foreground->detector_names();
      const vector<string> &backdet = background->detector_names();
      const set<string> forenameset(foredet.begin(),foredet.end());
      const set<string> backnameset(backdet.begin(),backdet.end());
      
      background->shiftPeaksForRecalibration( displ_background->calibration_coeffs(),
                                              displ_background->deviation_pairs(),
                                              displ_background->energy_calibration_model(),
                                              eqn, deviationPairs, eqnType );

      
      if( forenameset == backnameset )
        background->recalibrate_by_eqn( eqn, deviationPairs, eqnType, detectors, false );
      else
        background->recalibrate_by_eqn( eqn, deviationPairs, eqnType );
    }//if( background )
    
    if( secondary && displ_secondary && (secondary != foreground) )
    {
      const vector<string> &foredet = foreground->detector_names();
      const vector<string> &secodet = secondary->detector_names();
      const set<string> forenameset(foredet.begin(),foredet.end());
      const set<string> seconameset(secodet.begin(),secodet.end());
      
      secondary->shiftPeaksForRecalibration( displ_secondary->calibration_coeffs(),
                                             displ_secondary->deviation_pairs(),
                                             displ_secondary->energy_calibration_model(),
                                             eqn, deviationPairs, eqnType );
      
      if( forenameset == seconameset )
        secondary->recalibrate_by_eqn( eqn, deviationPairs, eqnType, detectors, false );
      else
        secondary->recalibrate_by_eqn( eqn, deviationPairs, eqnType );
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
                         foreground, kForeground, orig_eqn, origdev, eqnType );
  
  // Redraw everything
  if( viewer->m_dataMeasurement )
    viewer->displayForegroundData( true );
  if( viewer->m_secondDataMeasurement )
    viewer->displaySecondForegroundData();
  if( viewer->m_backgroundMeasurement )
    viewer->displayBackgroundData();
  
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



PreserveCalibWindow::PreserveCalibWindow(
                             std::shared_ptr<SpecMeas> newmeas,
                             const SpectrumType newtype,
                             std::shared_ptr<SpecMeas> oldmeas,
                             const SpectrumType oldtype,
                             Recalibrator *calibrator )
: AuxWindow( "Keep Calibration?",
             Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneModal) ),
  m_calibrator( calibrator ),
  m_newmeas( newmeas ),
  m_newtype( newtype ),
  m_oldtype( oldtype )
{
  if( !newmeas || !oldmeas || !calibrator )
    throw runtime_error( "PreserveCalibWindow: invalid input" );
  
  std::shared_ptr<const Measurement> eqnoldmeas, eqnnewmeas;
  const vector< MeasurementConstShrdPtr > oldmeass = oldmeas->measurements();
  const vector< MeasurementConstShrdPtr > newmeass = newmeas->measurements();
  
  for( size_t i = 0; !eqnoldmeas && i < oldmeass.size(); ++i )
    if( oldmeass[i]->num_gamma_channels() )
      eqnoldmeas = oldmeass[i];
  for( size_t i = 0; !eqnnewmeas && i < newmeass.size(); ++i )
    if( newmeass[i]->num_gamma_channels() )
      eqnnewmeas = newmeass[i];
  
  if( !eqnoldmeas || !eqnnewmeas )
    throw runtime_error( "PreserveCalibWindow: invalid input 2" );

  
  m_type = eqnoldmeas->energy_calibration_model();
  m_coeffs = eqnoldmeas->calibration_coeffs();
  m_devPairs = eqnoldmeas->deviation_pairs();
  
  switch( m_type )
  {
    case Measurement::Polynomial:
    case Measurement::FullRangeFraction:
    case Measurement::UnspecifiedUsingDefaultPolynomial:
      break;
      
    case Measurement::LowerChannelEdge:
    case Measurement::InvalidEquationType:
      throw runtime_error( "PreserveCalibWindow: invalid type" );
    break;
  }//switch( m_devPairs )
  
  string msg = "<div>It looks like you have loaded a ";
  
  msg += descriptionText(newtype);
  msg += " spectrum";
  
  if( m_newtype == m_oldtype )
  {
    msg += " that is from the same detector as the previous spectrum";
  }else
  {
    msg += " that is from the same detector as the ";
    msg += descriptionText(m_oldtype);
    msg += " spectrum";
  }

  msg += ", but with a different calibration.</div>"
         "<table class=\"RecalibCoefTable\"><tr><th></th><th>Previous (";
  
  switch( eqnnewmeas->energy_calibration_model() )
  {
    case Measurement::Polynomial:
    case Measurement::UnspecifiedUsingDefaultPolynomial:
      msg += "Polynomial";
      break;
    case Measurement::FullRangeFraction: msg += "FRF"; break;
    case Measurement::LowerChannelEdge:
    case Measurement::InvalidEquationType:
      break;
  }//switch( newmeas->energy_calibration_model() )
  
  msg += ")</th><th>New (";
  switch( eqnoldmeas->energy_calibration_model() )
  {
    case Measurement::Polynomial:
    case Measurement::UnspecifiedUsingDefaultPolynomial:
      msg += "Polynomial";
      break;
    case Measurement::FullRangeFraction: msg += "FRF"; break;
    case Measurement::LowerChannelEdge:
    case Measurement::InvalidEquationType:
      break;
  }//switch( oldmeas->energy_calibration_model() )
  msg += ")</th></tr>";
  
  
  
  const vector<float> &newcoefs = eqnnewmeas->calibration_coeffs();
  const size_t ncoefs = std::max( newcoefs.size(), m_coeffs.size() );
  
  for( size_t i = 0; i < ncoefs; ++i )
  {
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "<tr><td>%i</td>", static_cast<int>(i) );
    msg += buffer;
    
    if( i < m_coeffs.size() )
    {
      snprintf( buffer, sizeof(buffer), "<td>%.4g</td>", m_coeffs[i] );
      msg += buffer;
    }else
      msg += "<td></td>";
    
    if( i < newcoefs.size() )
    {
      snprintf( buffer, sizeof(buffer), "<td>%.4g</td>", newcoefs[i] );
      msg += buffer;
    }else
      msg += "<td></td>";
    msg += "</tr>";
  }//for( size_t i = 0; i < ncoefs; ++i )
  
  msg += "</table>";
  
  msg += "<div style=\"margin-top: 10px;\"><b>"
         "Would you like to use the previous calibration with the new file?"
         "</b></div>";
  
  new WText( msg, XHTMLUnsafeText, contents() );
  
  
  WContainerWidget *bottom = footer();
  //bottom->setStyleClass("modal-footer");
  WPushButton *cancel = new WPushButton( "No", bottom );
  cancel->clicked().connect( this, &PreserveCalibWindow::emitReject );
  cancel->setWidth( WLength(47.5,WLength::Percentage) );
  cancel->setFloatSide(Wt::Left);
  
  WPushButton *use = new WPushButton( "Yes", bottom );
  use->clicked().connect( this, &PreserveCalibWindow::doRecalibration );
  use->setWidth( WLength(47.5,WLength::Percentage) );
  use->setFloatSide( Wt::Right );
  
  rejectWhenEscapePressed();
  AuxWindow::disableCollapse();
  AuxWindow::show();
  AuxWindow::resize( 400.0, WLength::Auto );
  AuxWindow::centerWindow();
}//PreserveCalibWindow


PreserveCalibWindow::~PreserveCalibWindow()
{
}


bool PreserveCalibWindow::candidate( std::shared_ptr<SpecMeas> newmeas,
                                     std::shared_ptr<SpecMeas> oldmeas )
{
  if( !newmeas || !oldmeas || newmeas==oldmeas )
    return false;
  
  
  std::shared_ptr<const Measurement> eqnoldmeas, eqnnewmeas;
  const vector< MeasurementConstShrdPtr > oldmeass = oldmeas->measurements();
  const vector< MeasurementConstShrdPtr > newmeass = newmeas->measurements();
  
  for( size_t i = 0; !eqnoldmeas && i < oldmeass.size(); ++i )
    if( oldmeass[i]->num_gamma_channels() )
      eqnoldmeas = oldmeass[i];
  for( size_t i = 0; !eqnnewmeas && i < newmeass.size(); ++i )
    if( newmeass[i]->num_gamma_channels() )
      eqnnewmeas = newmeass[i];
  
  if( !eqnoldmeas || !eqnnewmeas )
    return false;
  
  if( newmeas->instrument_id() != oldmeas->instrument_id() )
    return false;
  
  if( newmeas->num_gamma_channels() != oldmeas->num_gamma_channels() )
    return false;
  
  
  
  const vector<float> &oldcoefs = eqnoldmeas->calibration_coeffs();
  const vector< pair<float,float> > &olddevpairs = eqnoldmeas->deviation_pairs();
  const Measurement::EquationType oldtype = eqnoldmeas->energy_calibration_model();
  
  const vector<float> &newcoefs = eqnnewmeas->calibration_coeffs();
  const vector< pair<float,float> > &newdevpairs = eqnnewmeas->deviation_pairs();
  const Measurement::EquationType newtype = eqnnewmeas->energy_calibration_model();

  if( oldcoefs.size() != newcoefs.size() )
    return true;
  if( olddevpairs.size() != newdevpairs.size() )
    return true;
  if( oldtype != newtype )
    return true;
  
  for( size_t i = 0; i < oldcoefs.size(); ++i )
  {
    const float a = oldcoefs[i];
    const float b = newcoefs[i];
    const float diff = fabs(a - b);
    
    //FLT_EPSILON the minimum positive number such that 1.0 + FLT_EPSILON != 1.0.
    //FLT_EPSILON==1.19209290E-07F
    //FLT_MIN==1.17549435E-38F
    if( (diff > (1.0E-5*std::max(fabs(a),fabs(b))))
        && diff > (std::pow(1.0E-3f,static_cast<float>(i+1))) )
      return true;
  }//for( size_t i = 0; i < oldcoefs.size(); ++i )
  
  for( size_t i = 0; i < olddevpairs.size(); ++i )
  {
    float a = olddevpairs[i].first;
    float b = newdevpairs[i].first;
    float diff = fabs(a - b);
    
    if( (diff > (1.0E-5*std::max(fabs(a),fabs(b))))
       && diff > (std::pow(1.0E-3f,static_cast<float>(i+1))) )
      return true;
    
    a = olddevpairs[i].second;
    b = newdevpairs[i].second;
    diff = fabs(a - b);
    
    if( (diff > (1.0E-5*std::max(fabs(a),fabs(b))))
       && diff > (std::pow(1.0E-3f,static_cast<float>(i+1))) )
      return true;
  }//for( size_t i = 0; i < olddevpairs.size(); ++i )
  
  return false;
}//candidate(...)

void PreserveCalibWindow::doRecalibration()
{
  InterSpec *viewer = m_calibrator->m_hostViewer;
  std::shared_ptr<SpecMeas> meas = viewer->measurment(m_newtype);
  std::shared_ptr<const Measurement> displ_foreground
                                      = viewer->displayedHistogram(kForeground);

  
  if( m_newmeas != meas )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  
  std::shared_ptr<const Measurement> eqnmeas;
  const vector< MeasurementConstShrdPtr > meass = meas->measurements();
  
  for( size_t i = 0; !eqnmeas && i < meass.size(); ++i )
    if( meass[i]->num_gamma_channels() )
      eqnmeas = meass[i];
  if( !eqnmeas )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  
  const vector<float> oldcoefs = eqnmeas->calibration_coeffs();
  const vector< pair<float,float> > olddevpairs = eqnmeas->deviation_pairs();
  const Measurement::EquationType oldtype = eqnmeas->energy_calibration_model();
  
  const vector<string> detectors = viewer->displayed_detector_names();
  meas->recalibrate_by_eqn( m_coeffs, m_devPairs, m_type, detectors, false );

  Recalibrator::shiftPeaksForEnergyCalibration( m_calibrator->m_peakModel,
                                                m_coeffs, m_devPairs, m_type,
                                                meas, m_newtype, oldcoefs, olddevpairs, oldtype );
  
  viewer->refreshDisplayedCharts();
  m_calibrator->refreshRecalibrator();
  
  finished().emit(WDialog::Accepted);
}//void doRecalibration()



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
  InterSpec *viewer = m_calibrator->m_hostViewer;
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
  
  for( int i = 0; i < 3; ++i )
  {
    m_calVal[i] = pow( 10.0, m_calibrator->m_coeffExps[i] )
                       * m_calibrator->m_coefficientDisplay[i]->value();
    m_calUncert[i] = -1.0;
    
    WLabel *label = 0;
    switch( i )
    {
      case 0: label = new WLabel( "Offset" ); break;
      case 1: label = new WLabel( "Linear" ); break;
      case 2: label = new WLabel( "Quadratic" ); break;
    }//switch( i )
    
    m_coefvals[i] = new WLineEdit();
    m_fitFor[i] = new WCheckBox( "Fit" );
    m_coefvals[i]->disable();
  
    fitForLayout->addWidget( label,         i, 0 );
    fitForLayout->addWidget( m_coefvals[i], i, 1 );
    fitForLayout->addWidget( m_fitFor[i],   i, 2 );
  }//for( int i = 0; i < 3; ++i )
  
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
    
    
    std::shared_ptr<const Measurement> data = m_calibrator->m_spectrumDisplayDiv->histUsedForXAxis();
    std::shared_ptr<const SpecMeas> meas = m_calibrator->m_hostViewer->measurment(kForeground);
    
    if( !meas || !data )
    {
      const char *msg = "You need to be displaying a foreground spectrum to do a calibration fit";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }
    
    
    std::shared_ptr<const Measurement> eqnmeas;
    const vector< MeasurementConstShrdPtr > meass = meas->measurements();
    
    for( size_t i = 0; !eqnmeas && i < meass.size(); ++i )
      if( meass[i]->num_gamma_channels() )
        eqnmeas = meass[i];
    if( !eqnmeas )
    {
      const char *msg = "You need to be displaying a foreground spectrum to do a calibration fit (unexpected error)";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }//if( !eqnmeas )
    
    const ShrdConstFVecPtr &binning = eqnmeas->channel_energies();
    
    if( !binning || binning->size() < 16 )
    {
      const char *msg = "The spectrum isnt high enough resolution to fit for a calibration";
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
      return;
    }
    
    const size_t nbin = binning->size();
    const Measurement::EquationType calibration_type = eqnmeas->energy_calibration_model();
    
    vector<float> calib_coefs = eqnmeas->calibration_coeffs();
    
    switch( calibration_type )
    {
      case Measurement::Polynomial:
      case Measurement::FullRangeFraction:
      case Measurement::UnspecifiedUsingDefaultPolynomial:
        break;
        
      case Measurement::LowerChannelEdge:
      case Measurement::InvalidEquationType:
      {
        passMessage( "Invalid starting calibration type: unknown or lower bin "
                     "edge energy not allowed", "", WarningWidget::WarningMsgHigh );
        return;
        break;
      }//case LowerChannelEdge or InvalidEquationType
    }//switch( calibration_type )
    
    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contian the error of the fit mean for that peak
    vector<RecalPeakInfo> peakInfos;
    
    for( size_t peakn = 0; peakn < npeaks; ++peakn )
    {
      const PeakDef &peak = *peakstouse[peakn];
      
      const double wantedEnergy = peak.gammaParticleEnergy();
      
      RecalPeakInfo peakInfo;
      peakInfo.peakMean = peak.mean();
      peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
      if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
        peakInfo.peakMeanUncert = 0.5;
        
      peakInfo.photopeakEnergy = wantedEnergy;
      if( m_calibrator->m_coeffEquationType == Measurement::FullRangeFraction )
      {
        peakInfo.peakMeanBinNumber = find_bin_fullrangefraction( peak.mean(),
                                                                  calib_coefs, nbin,
                                                                  eqnmeas->deviation_pairs(),
                                                                  0.001f );
      }else
      {
        const vector<float> fwfcoef = polynomial_coef_to_fullrangefraction( calib_coefs, nbin );
        peakInfo.peakMeanBinNumber = find_bin_fullrangefraction( peak.mean(),
                                                                  fwfcoef, nbin,
                                                                  eqnmeas->deviation_pairs(),
                                                                  0.001f );
      }//if( FullRangeFraction ) / else
        
      if( IsNan(peakInfo.peakMeanBinNumber)
          || IsInf(peakInfo.peakMeanBinNumber) )
        throw runtime_error( "Invalid result fromm "
                             "find_bin_fullrangefraction(...)" );
      peakInfos.push_back( peakInfo );
    }//for( int col = 0; col < numModelCol; ++col )
    
  
    PolyCalibCoefMinFcn chi2Fcn( peakInfos,
                                binning->size(),
                                m_calibrator->m_coeffEquationType,
                                eqnmeas->deviation_pairs() );
    ROOT::Minuit2::MnUserParameters inputPrams;
    
    if( calib_coefs.size() < 3 )
      calib_coefs.resize( 3, 0.0 );
    
    //    for( int i = 0; i < calib_coefs.size(); ++i )
    //      cerr << "initial calib_coefs[" << i << "]=" << calib_coefs[i] << endl;
    
    string name = "0";
    double delta = 1.0/pow( static_cast<double>(nbin), 0.0 );
    if( m_calibrator->m_coeffEquationType == Measurement::FullRangeFraction )
      delta = 0.5;
    
    if( m_fitFor[0]->isChecked() )
      inputPrams.Add( name, calib_coefs[0], delta );
    else
      inputPrams.Add( name, calib_coefs[0] );
    
    name = "1";
    delta = 1.0/pow( static_cast<double>(data->GetNbinsX()), 1.0 );
    if( m_calibrator->m_coeffEquationType == Measurement::FullRangeFraction )
      delta = 1.0;
    
    if( m_fitFor[1]->isChecked() )
    {
      inputPrams.Add( name, calib_coefs[1], delta );
      inputPrams.SetLowerLimit( name, 0.0 );
    }else
    {
      inputPrams.Add( name, calib_coefs[1] );
    }
    
    name = "2";
    delta = 1.0/pow( static_cast<double>(data->GetNbinsX()), 2.0 );
    if( m_calibrator->m_coeffEquationType == Measurement::FullRangeFraction )
      delta = 0.1;
    
    if( m_fitFor[2]->isChecked() )
      inputPrams.Add( name, calib_coefs[2], delta );
    else
      inputPrams.Add( name, calib_coefs[2] );
    
    ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
    ROOT::Minuit2::MnStrategy strategy( 1 ); //0 low, 1 medium, >=2 high
    
    
    const int npars = inputPrams.VariableParameters();
    const unsigned int maxFcnCall = 200 + 100*npars + 5*npars*npars;
    const double tolerance = 1.0;
    
//    ROOT::Minuit2::MnMigrad fitter( chi2Fcn, inputParamState, strategy );
//    ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
    
    ROOT::Minuit2::CombinedMinimizer fitter;
    ROOT::Minuit2::FunctionMinimum minimum
                          = fitter.Minimize( chi2Fcn, inputParamState,
                                            strategy, maxFcnCall, tolerance );
    
    
    //Not sure why Minuit2 doesnt like converging on the minumum verry well, but
    //  rather than showing the user an error message, we'll give it anither try
    if( minimum.IsAboveMaxEdm() )
    {
      ROOT::Minuit2::MnMigrad fitter( chi2Fcn, inputParamState, strategy );
      minimum = fitter( maxFcnCall, tolerance );
    }
    
    if( !minimum.IsValid() )
    {
      if( !minimum.HasValidCovariance() || !minimum.HasValidParameters() )
      {
        string msg = "Fit for calibration parameters failed.";
        if( m_fitFor[2] )
          msg += " you might try not fitting for quadratic term.";
        passMessage( msg, "", WarningWidget::WarningMsgHigh );
        return;
      }//if( !minimum.HasValidCovariance() || !minimum.HasValidParameters() )
      
      cerr << minimum << endl;
      stringstream msg;
      msg << "Warning: calibration coefficient fit results "
      "may be invalid, please check, and if necessary "
      "cancel this calibration.";
      
      if( minimum.IsAboveMaxEdm() )
        msg << " The estimated distance to optimal calibration parameters is "
        "too large.";
      if( minimum.HasReachedCallLimit() )
        msg << " To many calls to chi2 routine were made.";
      
      passMessage( msg.str(), "",WarningWidget::WarningMsgHigh );
    }//if( !minimum.IsValid() )
    
    const ROOT::Minuit2::MnUserParameters &fitPrams = minimum.UserState().Parameters();
    const vector<double> parValues = fitPrams.Params();
    const vector<double> parErrors = fitPrams.Errors();
    
    //    for( int i = 0; i < parValues.size(); ++i )
    //      cerr << "parValues[" << i << "]=" << parValues[1] << endl;
    
    for( size_t i = 0; i < parValues.size(); ++i )
      if( IsInf(parValues[i]) || IsNan(parValues[i]) )
        throw runtime_error( "Invalid calibration parameter from fit :(" );
    
    assert( parValues.size() == 3 );
    
    m_eqnType = calibration_type;
    for( int i = 0; i < 3; ++i )
    {
      m_calVal[i] = parValues[i];
      m_calUncert[i] = parErrors[i];
    }//for( size_t i = 0; i < parValues.size(); ++i )
    
    //Try to loop over peaks to give chi2 values and such
    vector<float> float_coef;
    for( const double d : parValues )
      float_coef.push_back( static_cast<float>(d) );
    
    if( m_calibrator->m_coeffEquationType == Measurement::Polynomial
       || m_calibrator->m_coeffEquationType == Measurement::UnspecifiedUsingDefaultPolynomial )
      float_coef = polynomial_coef_to_fullrangefraction( float_coef, nbin );
    
    stringstream msg;
    for( const RecalPeakInfo &info : peakInfos )
    {
      const double predictedMean = fullrangefraction_energy( info.peakMeanBinNumber, float_coef, nbin, eqnmeas->deviation_pairs() );
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
    
    if( UtilityFunctions::starts_with( exceptionmsg, "ErrorMsg" ) )
      msg = exceptionmsg.substr(8);
    
    cerr << SRC_LOCATION << "\n\tCaught: " << exceptionmsg << endl;
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void doFit()


void Recalibrator::MultiFileCalibFit::updateCoefDisplay()
{
  for( int i = 0; i < 3; ++i )
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
  }//for( int i = 0; i < 3; ++i )
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
      InterSpec *viewer = m_calibrator->m_hostViewer;
      
      //XXX - the equation could totally be invalid
      vector<float> eqn;
      eqn.push_back( m_calVal[0] );
      eqn.push_back( m_calVal[1] );
      if( m_calVal[2] != 0.0 )
        eqn.push_back( m_calVal[2] );
      
      
      cerr << "\n\nm_calVal={" << m_calVal[0] << ", " << m_calVal[1] << ", " << m_calVal[2] << "}" << endl;
      
      vector<float> oldcalibpars;
      DeviationPairVec devpairs;
      Measurement::EquationType oldEqnType;
      
      std::shared_ptr<SpecMeas> fore = viewer->measurment(kForeground);
      std::shared_ptr<SpecMeas> back = viewer->measurment(kBackground);
      std::shared_ptr<SpecMeas> second = viewer->measurment(kSecondForeground);
      std::shared_ptr<const Measurement> displ_foreground
                                     = viewer->displayedHistogram(kForeground);

      
      std::shared_ptr<const Measurement> eqnmeass[3];
      for( int i = 0; i < 3; ++i )
      {
        const SpectrumType type = SpectrumType(i);
        std::shared_ptr<SpecMeas> meas = viewer->measurment(type);
        if( meas )
        {
          const vector< MeasurementConstShrdPtr > meass = meas->measurements();
          for( size_t j = 0; !eqnmeass[i] && j < meass.size(); ++j )
            if( meass[j]->num_gamma_channels() )
              eqnmeass[i] = meass[j];
        }
      }//for( int i = 0; i < 3; ++i )
      
      if( fore && eqnmeass[kForeground] )
      {
        oldcalibpars = eqnmeass[kForeground]->calibration_coeffs();
        oldEqnType = eqnmeass[kForeground]->energy_calibration_model();
        devpairs = eqnmeass[kForeground]->deviation_pairs();
      }else if( second && eqnmeass[kSecondForeground] )
      {
        oldcalibpars = eqnmeass[kSecondForeground]->calibration_coeffs();
        oldEqnType = eqnmeass[kSecondForeground]->energy_calibration_model();
        devpairs = eqnmeass[kSecondForeground]->deviation_pairs();
      }else if( back && eqnmeass[kBackground] )
      {
        oldcalibpars = eqnmeass[kBackground]->calibration_coeffs();
        oldEqnType = eqnmeass[kBackground]->energy_calibration_model();
        devpairs = eqnmeass[kBackground]->deviation_pairs();
      }else
      {
        break;
      }
      
      vector<string> displayed_detectors = viewer->displayed_detector_names();
      
      
      if( !m_calibrator->m_applyTo[kForeground]->isChecked() )
        fore.reset();
      if( !m_calibrator->m_applyTo[kBackground]->isChecked() )
        back.reset();
      if( !m_calibrator->m_applyTo[kSecondForeground]->isChecked() )
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
        
        meas->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
        cerr << "\n\nRecalled " << header->displayName().toUTF8()
             << " to m_calVal={" << m_calVal[0] << ", " << m_calVal[1]
             << ", " << m_calVal[2] << "}" << endl;
        
        if( meas == fore )
        {
          fore.reset();
          Recalibrator::shiftPeaksForEnergyCalibration(
                                            m_calibrator->m_peakModel,
                                            eqn, devpairs, m_eqnType,
                                            meas, kForeground,
                                            oldcalibpars, devpairs, oldEqnType );
        }else
        {
          meas->shiftPeaksForRecalibration( oldcalibpars, devpairs, oldEqnType,
                                           eqn, devpairs, m_eqnType );
        }//if( meas == fore ) / else
        
        if( meas == back )
          back.reset();
        if( meas == second )
          second.reset();
      }//for( size_t i = 0; i < m_peaks.size(); ++i )
      
      if( fore )
        fore->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      if( back )
        back->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      if( second )
        second->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      
      viewer->displayForegroundData( true );
      viewer->displaySecondForegroundData();
      viewer->displayBackgroundData();
      
      m_calibrator->refreshRecalibrator();
      
      cerr << "\nAccepted Recalibrator::MultiFileCalibFit" << endl;
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  delete this;
}//void handleFinish(...)




