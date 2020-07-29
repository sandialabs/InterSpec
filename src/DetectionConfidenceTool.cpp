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

#include <iostream>

#include <boost/math/distributions/chi_squared.hpp>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WComboBox>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WTableCell>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>



//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MinosError.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FumiliMinimizer.h"
#include "Minuit2/ScanMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMinimize.h"



#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/DetectorEdit.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/DetectionConfidenceTool.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

using namespace std;
using namespace Wt;

#define MU_CHARACTER "\xCE\xBC"

namespace
{

class MdaPeakRow : public WContainerWidget
{
public:
  const float m_energy;
  const double m_gammas_per_bq;
  const float m_fwhm;
  WCheckBox *m_use;
  NativeFloatSpinBox *m_roi_start;
  NativeFloatSpinBox *m_roi_end;
  WComboBox *m_continuum;
  WText *m_poisonLimit;
  double m_simpleMda;
  shared_ptr<const SpecUtils::Measurement> m_meas;
  
  Wt::Signal<> m_changed;
  
  
  void setSimplePoisonTxt()
  {
    if( !m_meas || (m_meas->num_gamma_channels() < 16) )
      return;
    
    const float roi_lower_energy = m_roi_start->value();
    const float roi_upper_energy = m_roi_end->value();
    
    const size_t nsidebin = 1;
    const size_t nchannels = m_meas->num_gamma_channels();
    
    // \TODO: check that we are getting the right channels, and not off by one or something.
    size_t lower_bin = m_meas->find_gamma_channel(roi_lower_energy);
    lower_bin = ((lower_bin > (2*nsidebin+1)) ? (lower_bin - nsidebin) : (nsidebin + 1));
    
    size_t upper_bin = m_meas->find_gamma_channel(roi_lower_energy);
    upper_bin = ((upper_bin < (nchannels - 2*nsidebin - 1)) ? (upper_bin + nsidebin) : (nchannels - nsidebin - 1));
    
    lower_bin = std::min( lower_bin, nchannels - nsidebin - 1 );
    upper_bin = std::max( upper_bin, nsidebin + 1 );
    
    double coefs[2];
    PeakContinuum::eqn_from_offsets( lower_bin, upper_bin, m_energy, m_meas, nsidebin, coefs[1], coefs[0] );
    
    double continuum = PeakContinuum::offset_eqn_integral( coefs,
                                       PeakContinuum::OffsetType::Linear,
                                       roi_lower_energy, roi_upper_energy, m_energy );
    const double data_area = m_meas->gamma_integral(roi_lower_energy, roi_upper_energy);
    
    
    cout << "Predicted continuum area=" << continuum << " vs data " << data_area << endl;
    
    continuum = data_area;
    
    m_simpleMda = (2.71 + 4.65*sqrt(continuum))/m_gammas_per_bq;
    if( IsInf(m_simpleMda) || IsInf(m_simpleMda) )
      m_simpleMda = -1.0;
    const string mdastr = PhysicalUnits::printToBestActivityUnits(m_simpleMda);
    cout << "mda(" << m_energy << ") = " << mdastr << endl;
    
    m_poisonLimit->setText( "Simple MDA=" + mdastr );
  }//void setSimplePoisonTxt()
  
  void emitChanged()
  {
    setSimplePoisonTxt();
    m_changed.emit();
  }
  
public:
  MdaPeakRow( const float energy, const double gammasPerBq ,
              const float fwhm,
              float roi_start, float roi_end,
              shared_ptr<const SpecUtils::Measurement> meas,
             Wt::WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
  m_energy( energy ),
  m_gammas_per_bq( gammasPerBq ),
  m_fwhm( fwhm ),
  m_simpleMda( -1.0 ),
  m_meas( meas ),
  m_changed( this )
  {
    addStyleClass( "MdaPeakRow" );
    
    m_use = new WCheckBox( "", this );
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "&nbsp;%.2f keV, FWHM=%.2f, Eff=%.2g/bq&nbsp;",
              m_energy, m_fwhm, m_gammas_per_bq );
    auto label = new WLabel( buffer, this );
    label->addStyleClass( "FixedInfo" );
    
    label = new WLabel( "&nbsp;ROI Start:", this );
    m_roi_start = new NativeFloatSpinBox( this );
    m_roi_start->setValue( roi_start );
    //m_roi_start->setTextSize(7); //doesnt seem to work
    m_roi_start->setWidth(75);
    
    label = new WLabel( "keV, &nbsp;ROI End:", this );
    m_roi_end = new NativeFloatSpinBox( this );
    m_roi_end->setValue( roi_end );
    //m_roi_end->setTextSize(7); //doesnt seem to work
    m_roi_end->setWidth(75);
    
    label = new WLabel( "keV,&nbsp;Continuum", this );
    m_continuum = new WComboBox( this );
    m_continuum->addItem( "Linear" );
    m_continuum->addItem( "Quadratic" );
    m_continuum->setCurrentIndex( 0 );
    
    m_roi_start->valueChanged().connect( this, &MdaPeakRow::emitChanged );
    m_roi_end->valueChanged().connect( this, &MdaPeakRow::emitChanged );
    m_continuum->activated().connect( this, &MdaPeakRow::emitChanged );
    m_continuum->changed().connect( this, &MdaPeakRow::emitChanged );
    
    m_poisonLimit = new WText( "&nbsp;", this );
    m_poisonLimit->addStyleClass( "Poisson");
    setSimplePoisonTxt();
  }
};//class MdaPeakRow

}//namespace



DetectionConfidenceWindow::DetectionConfidenceWindow( InterSpec *viewer,
                                                     MaterialDB *materialDB,
                                                     WSuggestionPopup *materialSuggest )
: AuxWindow( "Detection Confidence Tool",
  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletModal)
   | AuxWindowProperties::SetCloseable
   | AuxWindowProperties::DisableCollapse) ),
  m_tool( nullptr )
{
  assert( viewer );
  rejectWhenEscapePressed( true );
  
  m_tool = new DetectionConfidenceTool( viewer, materialDB, materialSuggest );
  //stretcher()->addWidget( m_tool, 0, 0 );
  contents()->addWidget( m_tool );
  
  AuxWindow::addHelpInFooter( footer(), "detection-confidence-tool" );
  
  //WContainerWidget *buttonDiv = footer();
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  const int width = ((screenW < 906) ? screenW : 0.75*screenW );
  const int height = ((screenH < 525) ? screenH : 0.95*screenH );
  resizeWindow( width, height );
  
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//DetectionConfidenceWindow(...) constrctor


DetectionConfidenceWindow::~DetectionConfidenceWindow()
{
}




DetectionConfidenceTool::DetectionConfidenceTool( InterSpec *viewer,
                                                  MaterialDB *materialDB,
                                                  Wt::WSuggestionPopup *materialSuggest,
                                                  WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_interspec( viewer ),
    m_chart( nullptr ),
    m_peakModel( nullptr ),
    m_nuclideEdit( nullptr ),
    m_nuclideSuggest( nullptr ),
    m_detectorDisplay( nullptr ),
    m_distanceEdit( nullptr ),
    m_materialDB( materialDB ),
    m_materialSuggest( materialSuggest ),
    m_shieldingSelect( nullptr )
{
  /** \TODO:
   - Put the Chi2 chart to the right of the spectrum, when it should exist
   - Give the user a choice about using continuum fixed at null hypothesis
   - Allow combining ROI with neghboring peaks
   - Have the energy rows fold down to show more information, similar to Steves tool, for each energy
   - Allow users to select CL, not just 95%
   */
  wApp->useStyleSheet( "InterSpec_resources/DetectionConfidenceTool.css" );
  
  addStyleClass( "DetectionConfidenceTool" );
  
  //new WLabel( "DetectionConfidenceTool", this );
  const WLength labelWidth(3.5,WLength::FontEm), fieldWidth(4,WLength::FontEm);
  const WLength optionWidth(5.25,WLength::FontEm), buttonWidth(5.25,WLength::FontEm);

  
  WGridLayout *grid = new WGridLayout();
  setLayout( grid );
  
  m_chart = new D3SpectrumDisplayDiv();
  
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  const int width = ((screenW < 906) ? screenW : 0.75*screenW ); // maybe cap width off at 1800 px?
  const int height = ((screenH < 525) ? screenH : 0.8*screenH );
  //m_chart->resize( WLength(width-20,WLength::Unit::Pixel), WLength(std::min(0.3*height,300.0),WLength::Unit::Pixel) );
  //m_chart->resize( WLength::Auto, WLength(300,WLength::Unit::Pixel) );
  cout << "height-->" << 1.0*height << endl;
  m_chart->setMinimumSize( WLength::Auto, WLength( std::max(200.0,0.4*height),WLength::Unit::Pixel) );
  
  
  m_chart->setCompactAxis( true );
  m_chart->setXAxisTitle( "Energy (keV)" );
  m_chart->setYAxisTitle( "Counts/Channel" );
  m_chart->enableLegend( false );
  m_chart->showHistogramIntegralsInLegend( true );
  
  auto theme = m_interspec->getColorTheme();
  if( theme )
  {
    m_chart->setForegroundSpectrumColor( theme->foregroundLine );
    m_chart->setBackgroundSpectrumColor( theme->backgroundLine );
    m_chart->setSecondarySpectrumColor( theme->secondaryLine );
    m_chart->setDefaultPeakColor( theme->defaultPeakLine );
  
    m_chart->setAxisLineColor( theme->spectrumAxisLines );
    m_chart->setChartMarginColor( theme->spectrumChartMargins );
    m_chart->setChartBackgroundColor( theme->spectrumChartBackground );
    m_chart->setTextColor( theme->spectrumChartText );
  }//if( theme )
  
  
  m_peakModel = new PeakModel( this );
  m_chart->setPeakModel( m_peakModel );
  
  grid->addWidget( m_chart, 0, 0 );
  
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  
  
  WTable *inputTable = new WTable();
  inputTable->addStyleClass( "UserInputTable" );
  grid->addWidget( inputTable, 1, 0, AlignLeft );
  
  WTableCell *cell = inputTable->elementAt(0, 0);
  WLabel *label = new WLabel( "Distance:", cell );
  label->setMinimumSize( labelWidth, WLength::Auto );
  //grid->addWidget( label, 1, 0 );
  
  cell = inputTable->elementAt(0, 1);
  m_distanceEdit = new WLineEdit( "100 cm", cell );
  //m_distanceEdit->setTextSize( 5 );
  m_distanceEdit->setValidator( distValidator );
  m_distanceEdit->setMinimumSize( fieldWidth, WLength::Auto );
  //grid->addWidget( m_distanceEdit, 1, 1, AlignLeft );
  label->setBuddy( m_distanceEdit );
  m_distanceEdit->changed().connect( this, &DetectionConfidenceTool::handleInputChange );
  m_distanceEdit->enterPressed().connect( this, &DetectionConfidenceTool::handleInputChange );
  
  cell = inputTable->elementAt(1, 0);
  label = new WLabel( "Nuclide:", cell );
  label->setMinimumSize( labelWidth, WLength::Auto );
  //grid->addWidget( label, 2, 0 );
  
  cell = inputTable->elementAt(1, 1);
  m_nuclideEdit = new WLineEdit( "", cell );
  m_nuclideEdit->setMargin( 1 );
  m_nuclideEdit->setMinimumSize( fieldWidth, WLength::Auto );
  m_nuclideEdit->setAutoComplete( false );
  m_nuclideEdit->changed().connect( boost::bind( &DetectionConfidenceTool::handleInputChange, this ) );
  label->setBuddy( m_nuclideEdit );
  //grid->addWidget( m_nuclideEdit, 2, 1, AlignLeft );
  

  
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
  isoSuggestModel->excludeReactions( true );
  isoSuggestModel->excludeEscapes( true );
  isoSuggestModel->excludeXrays( true );
  m_nuclideSuggest = new WSuggestionPopup( matcherJs, replacerJs, this );
  m_nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
  m_nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );

  IsotopeNameFilterModel::setQuickTypeFixHackjs( m_nuclideSuggest );
  
  isoSuggestModel->filter( "" );
  m_nuclideSuggest->setFilterLength( -1 );
  m_nuclideSuggest->setModel( isoSuggestModel );
  m_nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
  m_nuclideSuggest->forEdit( m_nuclideEdit, WSuggestionPopup::Editing );  // | WSuggestionPopup::DropDownIcon

  ReferencePhotopeakDisplay *reflines = viewer->referenceLinesWidget();
  if( reflines )
  {
    const ReferenceLineInfo &current = reflines->currentlyShowingNuclide();
    if( current.nuclide )
      m_nuclideEdit->setText( current.nuclide->symbol );
  }
  
  SpectraFileModel *specFileModel = viewer->fileManager()->model();
  m_detectorDisplay = new DetectorDisplay( viewer, specFileModel );
  viewer->detectorChanged().connect( boost::bind( &DetectionConfidenceTool::handleInputChange, this ) );
  viewer->detectorModified().connect( boost::bind( &DetectionConfidenceTool::handleInputChange, this ) );
  
  cell = inputTable->elementAt(0, 3);
  cell->setRowSpan( 2 );
  cell->addWidget( m_detectorDisplay );
  
  //grid->addWidget( m_detectorDisplay, 3, 0, 1, 2, Wt::AlignLeft );
  
  cell = inputTable->elementAt(0, 4);
  label = new WLabel( "Peaks disp. act.:", cell );
  cell = inputTable->elementAt(0, 5);
  m_displayActivity = new WLineEdit( cell );
  WRegExpValidator *actvalidator = new WRegExpValidator( PhysicalUnits::sm_activityRegex, m_displayActivity );
  actvalidator->setFlags(Wt::MatchCaseInsensitive);
  m_displayActivity->setValidator( actvalidator );
  m_displayActivity->setTextSize( 10 );
  m_displayActivity->setText( "0 uCi" );
  m_displayActivity->enterPressed().connect( this, &DetectionConfidenceTool::updateShownPeaks );
  m_displayActivity->changed().connect( this, &DetectionConfidenceTool::updateShownPeaks );
  
  
  
  m_shieldingSelect = new ShieldingSelect( m_materialDB, NULL, m_materialSuggest, false );
  m_shieldingSelect->materialEdit()->setEmptyText( "<shielding material>" );
  m_shieldingSelect->materialChanged().connect( this, &DetectionConfidenceTool::handleInputChange );
  m_shieldingSelect->materialModified().connect( this, &DetectionConfidenceTool::handleInputChange );
  m_shieldingSelect->setMinimumSize( WLength(250), WLength::Auto );
  
  cell = inputTable->elementAt(0,2);
  cell->addWidget( m_shieldingSelect );
  cell->setRowSpan( 2 );
  //grid->addWidget( m_shieldingSelect, 4, 0, 1, 2, Wt::AlignLeft );
  
  const auto primaryMeas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  auto spec = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( primaryMeas && spec )
  {
    //Make a copy of things so we dont mess the real stuff up - eventually we may want to
    //  Fully copy the SpecMeas 
    auto ourspec = make_shared<SpecUtils::Measurement>( *spec );
    ourspec->set_detector_name( "Aa1" );
    m_our_meas = make_shared<SpecMeas>();
    m_our_meas->setDetector( primaryMeas->detector() );
    m_our_meas->add_measurement( ourspec, true );
    m_peakModel->setPeakFromSpecMeas( m_our_meas, {ourspec->sample_number()} );
    m_chart->setData( ourspec, false );
  }//if( spec )
  
  

  
  m_results = new WContainerWidget();
  m_results->addStyleClass( "MdaResults" );
  m_results->hide();
  grid->addWidget( m_results, 2, 0 );
  
  ///////////////////////////////////////////////////////////////////////
  // \TODO: put shart stuff, especially theme, into its own function(s)
  m_chi2Chart = new Wt::Chart::WCartesianChart( m_results );
  m_chi2Chart->setBackground(Wt::WColor(220, 220, 220));
  m_chi2Chart->setAutoLayoutEnabled();
  m_chi2Chart->setType(Wt::Chart::ScatterPlot);
  m_chi2Model = new WStandardItemModel( m_chart );
  m_chi2Chart->setModel( m_chi2Model );
  m_chi2Chart->setType( Chart::ScatterPlot );
  m_chi2Chart->setXSeriesColumn(0);
        
  if( theme )
  {
    if( !theme->spectrumChartText.isDefault() )
    {
      WPen txtpen(theme->spectrumChartText);
      m_chi2Chart->setTextPen( txtpen );
      m_chi2Chart->axis(Chart::XAxis).setTextPen( txtpen );
      m_chi2Chart->axis(Chart::YAxis).setTextPen( txtpen );
      m_chi2Chart->axis(Chart::Y2Axis).setTextPen( txtpen );
    }
        
    if( theme->spectrumChartBackground.isDefault() )
      m_chi2Chart->setBackground( Wt::NoBrush );
    else
      m_chi2Chart->setBackground( WBrush(theme->spectrumChartBackground) );
        
        
    if( (theme->spectrumChartMargins.isDefault() && !theme->spectrumChartBackground.isDefault()) )
    {
      //theme->spectrumChartBackground
    }else if( !theme->spectrumChartMargins.isDefault() )
    {
      //theme->spectrumChartMargins
    }
        
    if( !theme->spectrumAxisLines.isDefault() )
    {
      WPen defpen = m_chi2Chart->axis(Chart::XAxis).pen();
      defpen.setColor( theme->spectrumAxisLines );
      m_chi2Chart->axis(Chart::XAxis).setPen( defpen );
      m_chi2Chart->axis(Chart::Y1Axis).setPen( defpen );
      m_chi2Chart->axis(Chart::Y2Axis).setPen( defpen );
    }
  }//if( theme )
      
  m_chi2Chart->setPlotAreaPadding(5, Wt::Top);
  m_chi2Chart->setPlotAreaPadding(55, Wt::Bottom);
  m_chi2Chart->setPlotAreaPadding(55, Wt::Left);
  m_chi2Chart->setPlotAreaPadding(10, Wt::Right);
  m_chi2Chart->axis(Wt::Chart::XAxis).setTitle("Activity (" MU_CHARACTER "Ci)");
  
#define CHI_CHARACTER "\xCF\x87"
#define SQ_CHARACTER "\xC2\xB2"
  m_chi2Chart->axis(Wt::Chart::Y1Axis).setTitle( WString::fromUTF8( CHI_CHARACTER SQ_CHARACTER ) );
      
#if( WT_VERSION >= 0x3030400 )
  m_chi2Chart->axis(Wt::Chart::Y1Axis).setTitleOrientation( Wt::Vertical );
#endif
      
  m_chi2Chart->setLegendEnabled( false );
  m_chi2Chart->setHeight( WLength( std::max(250.0,0.35*height),WLength::Unit::Pixel) );
  m_chi2Chart->setWidth( 0.5*width );
  //////////////////////////////////////////////////////////////////////
  m_bestChi2Act = new WText( "&nbsp;", m_results );
  m_bestChi2Act->setInline( false );
  
  m_upperLimit = new WText( "&nbsp;", m_results );
  m_upperLimit->setInline( false );
  //
  
  
  m_errorMsg = new WText("&nbsp;");
  m_errorMsg->addStyleClass( "MdaErrMsg" );
  grid->addWidget( m_errorMsg, 3, 0 /*, AlignCenter | AlignMiddle */ );
  m_errorMsg->hide();
  
  m_peaks = new WContainerWidget();
  grid->addWidget( m_peaks, 4, 0 );
  m_peaks->addStyleClass( "MdaPeaks" );
  
  //grid->setRowResizable( 0, screenH/3 );
  //grid->setColumnStretch( 1, 1 );
  //grid->setRowStretch(grid->rowCount() -1, 1 );
  //grid->setRowStretch( 0, 1 );
  
  handleInputChange();
}//DetectionConfidenceTool constructor
  
  
DetectionConfidenceTool::~DetectionConfidenceTool()
{
}


void DetectionConfidenceTool::handleInputChange()
{
  vector<tuple<float,bool,float,float>> oldvalues;
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    oldvalues.push_back( {rw->m_energy, rw->m_use->isChecked(), rw->m_roi_start->value(), rw->m_roi_end->value()} );
  }
  
  
  m_peaks->clear();
  m_chi2Model->clear();
  m_results->hide();
  m_errorMsg->setText( "&nbsp;" );
  m_errorMsg->hide();
  
  auto spec = m_chart->data();
  if( !m_our_meas || !spec )
  {
    m_errorMsg->setText( "No Data Loaded" );
    m_errorMsg->show();
    return;
  }
  
  std::shared_ptr<DetectorPeakResponse> drf = m_our_meas->detector();
  if( !drf )
  {
    m_errorMsg->setText( "No DRF Loaded" );
    m_errorMsg->show();
    return;
  }
  
  if( !drf->hasResolutionInfo() || !drf->isValid() )
  {
    m_errorMsg->setText( "DRF does not have FWHM info" );
    m_errorMsg->show();
    return;
  }
  
  const string isotxt = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotxt );

  if( !nuc )
  {
    m_errorMsg->setText( "No valid nuclide" );
    m_errorMsg->show();
    return;
  }
  
  const string disttxt = m_distanceEdit->text().toUTF8();
  double distance = 0.0;
  try
  {
    distance = PhysicalUnits::stringToDistance(disttxt);
    
    if( distance <= 0.0 )
      throw runtime_error( "Distance cant be zero or negative" );
  }catch( std::exception & )
  {
    m_errorMsg->setText( "Distance not valid" );
    m_errorMsg->show();
    return;
  }//try / catch
  
  
  setRefLinesAndGetLineInfo();
  
  string agestr;
  double age = PeakDef::defaultDecayTime( nuc, &agestr );
  
  float shielding_an = 0.0f, shielding_ad = 0.0f, shielding_thickness = 0.0f;
  std::shared_ptr<Material> shielding_material;
  const bool generic_shielding = m_shieldingSelect->isGenericMaterial();
  if( generic_shielding )
  {
    shielding_an = m_shieldingSelect->atomicNumber();
    shielding_ad = m_shieldingSelect->arealDensity();
  }else
  {
    shielding_material = m_shieldingSelect->material();
    shielding_thickness = m_shieldingSelect->thickness();
  }//if( generic shielding ) / else
  
  
  
  std::vector<std::tuple<double,double>> lines;
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( nuc, 1.0E-3 * SandiaDecay::curie );

  const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( age );
  const double parent_activity = mixture.activity( age, nuc );
  
  vector<SandiaDecay::EnergyRatePair> gammas = mixture.gammas( age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
  
  
  try
  {
    boost::function<double(float)> att_coef_fcn;
    
    if( generic_shielding )
    {
      att_coef_fcn = boost::bind( &GammaInteractionCalc::transmition_coefficient_generic,
                                 shielding_an, shielding_ad, _1 );
    }else if( shielding_material && shielding_thickness > 0.0 )
    {
      att_coef_fcn = boost::bind( &GammaInteractionCalc::transmition_coefficient_material,
                                 shielding_material.get(), _1, shielding_thickness );
    }
    
    
    for( const auto &erp : gammas )
    {
      const double energy = erp.energy;
      double br = erp.numPerSecond / parent_activity;
      br *= drf->efficiency( static_cast<float>(energy), static_cast<float>(distance) );
      if( !att_coef_fcn.empty() )
        br *= exp( -1.0 * att_coef_fcn( energy ) );
      
      lines.push_back( {energy, br} );
    }//for( const auto &erp : gammas )
  }catch(...)
  {
    m_errorMsg->setText( "Error Calculating effects of shieldign" );
    m_errorMsg->show();
    return;
  }
  
  
  for( const auto &line : lines )
  {
    const float energy = get<0>(line);
    const double countsPerBq = get<1>(line)*spec->live_time();
    const float fwhm = drf->peakResolutionFWHM(energy);
    float roi_start = energy - 1.5*fwhm;
    float roi_end = energy + 1.5*fwhm;
    bool use = false;
    
   for( const auto &oldval : oldvalues )
   {
     const float old_energy = get<0>(oldval);
     if( old_energy != energy )
       continue;
     
     use = get<1>(oldval);
     roi_start = get<2>(oldval);
     roi_end = get<3>(oldval);
   }//for( const auto &oldval : oldvalues )
    
    
    auto row = new MdaPeakRow( energy, countsPerBq , fwhm, roi_start, roi_end, spec, m_peaks );
    row->m_use->setChecked(use);
    row->m_changed.connect( this, &DetectionConfidenceTool::doCalc );
  }//for( const auto &line : lines )
  
  doCalc();
}//void handleInputChange()


void DetectionConfidenceTool::doCalc()
{
  //m_simpleMda
  double minActivity = std::numeric_limits<double>::infinity(), maxActivity = 0.0;
  
  int nused = 0;
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    if( !rw->m_use->isChecked() )
      continue;
    
    ++nused;
    if( (rw->m_simpleMda > 0.0) && !IsInf(rw->m_simpleMda) && !IsNan(rw->m_simpleMda)  )
    {
      minActivity = std::min( minActivity, rw->m_simpleMda );
      maxActivity = std::max( maxActivity, rw->m_simpleMda );
    }
  }//for( auto w : m_peaks->children() )
  
  if( !nused )
    return;
  
  if( IsInf(minActivity) || maxActivity==0.0 )
  {
    minActivity = 0.0;
    maxActivity = 0.1*PhysicalUnits::curie;
  }else
  {
    minActivity *= 0.25;
    maxActivity *= 4.0;
  }
  
  size_t num_iterations = 0;
  int numDOF = 0;
  double overallBestChi2 = std::numeric_limits<double>::infinity(), overallBestActivity = 0.0;
  vector<pair<double,double>> chi2s;
  const size_t maxNumActivities = 50;
  while( chi2s.size() < (maxNumActivities/2) )
  {
    cout << "Searching range " << PhysicalUnits::printToBestActivityUnits(minActivity) << " to "
         << PhysicalUnits::printToBestActivityUnits(maxActivity) << endl;
    ++num_iterations;
    if( num_iterations > 100 )
    {
      cerr << "Too many iterations trying to find chart range" << endl;
      break;
    }
    size_t best_index = maxNumActivities + 1;
    vector<pair<double,double>> these_chi2s;
    const double step = (maxActivity - minActivity) / maxNumActivities;
    double bestChi2 = std::numeric_limits<double>::infinity(), bestChi2Activity = 0.0;
    for( double activity = minActivity; activity <= maxActivity; activity += step )
    {
      double chi2;
      std::vector<shared_ptr<PeakDef>> peaks;
      computeForAcivity( activity, peaks, chi2, numDOF );
      
      if( chi2 < bestChi2 )
      {
        bestChi2 = chi2;
        bestChi2Activity = activity;
        best_index = these_chi2s.size();
      }
      
      if( chi2 < overallBestChi2 )
      {
        overallBestChi2 = chi2;
        overallBestActivity = activity;
      }
      
      these_chi2s.push_back( {activity,chi2} );
    
      if( chi2==0.0 && numDOF == 0 )
        return;
    }//for( double activity = minActivity; activity <= maxActivity; activity += step )
    
    
    if( best_index > maxNumActivities )
    {
      cerr << "Totally unexpected ( best_index < 0 )" << endl;
      return;
    }
    
    size_t lower_index = best_index, upper_index = best_index;
    while( (lower_index > 0) && (these_chi2s[lower_index].second < (bestChi2 + 20)) )
      --lower_index;
    
    while( (upper_index < these_chi2s.size()) && (these_chi2s[upper_index].second < (bestChi2 + 20)) )
      ++upper_index;
    
    chi2s.clear();
    chi2s.insert( end(chi2s), begin(these_chi2s)+lower_index, begin(these_chi2s)+upper_index );
    if( chi2s.size() > 1 )
    {
      minActivity = chi2s.front().first;
      maxActivity = chi2s.back().first;
    }
  }//while( chi2s.size() && chi2s.size() < (maxNumActivities/2) )
  
  
  cout << "Fit Chi2s:" << endl;
  for( const auto &p : chi2s )
  {
    cout << PhysicalUnits::printToBestActivityUnits(p.first) << " --> chi2=" << p.second << endl;
  }
  cout << "And there were " << std::dec << numDOF << " DOFs" << endl;
  
  
  if( chi2s.empty() )
  {
    m_results->hide();
    return;
  }
  
  size_t bestIndex = 0;
  for( size_t i = 1; i < chi2s.size(); ++i )
  {
    if( chi2s[i].second < chi2s[bestIndex].second )
      bestIndex = i;
  }
  
  const double delta = 4.0; //95.45% CL, or something
  double upperActivity = chi2s[bestIndex].first, upperActivtyChi2 = chi2s[bestIndex].second;
  for( size_t i = bestIndex + 1; i < chi2s.size(); ++i )
  {
    upperActivity = chi2s[i].first;
    upperActivtyChi2 = chi2s[i].second;
    if( chi2s[i].second > (delta + chi2s[bestIndex].second) )
      break;
  }
  
  
  char buffer[128];
  snprintf( buffer, sizeof(buffer), "Best &chi;<sup>2</sup> of %.1f at activity %s, %i DOF",
            overallBestChi2,
            PhysicalUnits::printToBestActivityUnits(overallBestActivity).c_str(),
            numDOF );
  m_bestChi2Act->setText( buffer );
  
  string upperLimitActStr = PhysicalUnits::printToBestActivityUnits(upperActivity,3,true);
  snprintf( buffer, sizeof(buffer), "95%% coverage at %s with &chi;<sup>2</sup> of %.1f",
            upperLimitActStr.c_str(), upperActivtyChi2 );
  
  m_upperLimit->setText( buffer );
  
  
  const bool useCurries = true;
  const bool avrgAct = 0.5*(chi2s.front().first + chi2s.back().first);
  const auto &units = PhysicalUnits::bestActivityUnitHtml( avrgAct, useCurries);
  
  m_chi2Model->insertColumn( 0 );
  m_chi2Model->insertColumn( 1 );
  m_chi2Model->insertRows(0, static_cast<int>(chi2s.size()) );
  
  for( size_t i = 0; i < chi2s.size(); ++i )
  {
    m_chi2Model->setData( static_cast<int>(i), 0, (chi2s[i].first / units.second) );
    m_chi2Model->setData( static_cast<int>(i), 1, chi2s[i].second );
  }
  
  m_chi2Chart->axis(Wt::Chart::XAxis).setTitle("Activity (" + units.first + ")");
  m_chi2Model->setHeaderData(0, WString("Activity") );
  
  m_chi2Chart->setXSeriesColumn(0);
  Wt::Chart::WDataSeries series(1, Wt::Chart::LineSeries,Wt::Chart::Y1Axis);
  //series.setPen( WPen(m_chartFwhmLineColor) );
  m_chi2Chart->addSeries(series);
  
  //m_chi2Chart->axis(Wt::Chart::Y1Axis).setRoundLimits( Chart::AxisValue::MinimumValue | Chart::AxisValue::MaximumValue );
  m_chi2Chart->axis(Wt::Chart::Y1Axis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  m_chi2Chart->axis(Wt::Chart::XAxis).setRange( (chi2s.front().first/units.second), (chi2s.back().first/units.second) );
  //m_chi2Chart->axis(Wt::Chart::XAxis).setAutoLimits( Chart::MaximumValue );
  //m_chi2Chart->axis(Wt::Chart::XAxis).setRoundLimits( Chart::AxisValue::MinimumValue | Chart::AxisValue::MaximumValue );

  
  // \TODO: draw line at upperActivity.
  
  m_results->show();
  
  m_displayActivity->setText( WString::fromUTF8(upperLimitActStr) );
  updateShownPeaks();
}//void doCalc()


void DetectionConfidenceTool::updateShownPeaks()
{
  m_peakModel->removeAllPeaks();
  
  string actStr = m_displayActivity->text().toUTF8();
  double activity = 0;
  try
  {
    activity = PhysicalUnits::stringToActivity(actStr);
  }catch(...)
  {
    return;
  }//try / catch
  
  int numDOF;
  double chi2;
  std::vector<shared_ptr<PeakDef>> peaks;
  computeForAcivity( activity, peaks, chi2, numDOF );
  for( auto &p : peaks )
    m_peakModel->addNewPeak( *p );
  
  cout << "Done in updateShownPeaks()" << endl;
}//void DetectionConfidenceTool::updateShownPeaks()


void DetectionConfidenceTool::computeForAcivity( const double activity,
                                                 std::vector<shared_ptr<PeakDef>> &peaks,
                                                 double &chi2, int &numDOF )
{
  peaks.clear();
  chi2 = 0.0;
  numDOF = 0;
  
  auto spec = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !spec )
  {
    cerr << "No displayed histogram!" << endl;
    return;
  }
  
  vector<PeakDef> inputPeaks, fitPeaks;
  
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    if( !rw->m_use->isChecked() )
      continue;
    
    const float mean = rw->m_energy;
    const float fwhm = rw->m_fwhm;
    const float sigma = fwhm / 2.634;
    const float amplitude = activity * rw->m_gammas_per_bq;
    const float roi_start = rw->m_roi_start->value();
    const float roi_end = rw->m_roi_end->value();
    
    PeakDef peak( mean, sigma, amplitude );
    peak.setFitFor( PeakDef::CoefficientType::Mean, false );
    peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
    peak.setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
    shared_ptr<PeakContinuum> cont = peak.continuum();
    
    switch( rw->m_continuum->currentIndex() )
    {
      case 0:
        cont->setType( PeakContinuum::OffsetType::Linear );
        break;
      case 1:
        cont->setType( PeakContinuum::OffsetType::Quadratic );
        break;
    }//switch( assign continuum )
    
    cont->setRange( roi_start, roi_end );
    cont->calc_linear_continuum_eqn( spec, roi_start, roi_end, 2 );
    cont->setPolynomialCoefFitFor( 0, true );
    cont->setPolynomialCoefFitFor( 1, true );
    
    inputPeaks.push_back( std::move(peak) );
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  if( inputPeaks.empty() )
  {
    cerr << "No peaks to do calc for!" << endl;
    return;
  }
  
  ROOT::Minuit2::MnUserParameters inputFitPars;
  PeakFitChi2Fcn::addPeaksToFitter( inputFitPars, inputPeaks, spec, PeakFitChi2Fcn::kFitForPeakParameters );
      
  const int npeaks = static_cast<int>( inputPeaks.size() );
  PeakFitChi2Fcn chi2Fcn( npeaks, spec, nullptr );
  chi2Fcn.useReducedChi2( false );
      
  assert( inputFitPars.VariableParameters() != 0 );
          
  ROOT::Minuit2::MnUserParameterState inputParamState( inputFitPars );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( chi2Fcn, inputParamState, strategy );
      
  unsigned int maxFcnCall = 5000;
  double tolerance = 2.5;
  tolerance = 0.5;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  const ROOT::Minuit2::MnUserParameters &fitParams = fitter.Parameters();
  //  minimum.IsValid()
  //      ROOT::Minuit2::MinimumState minState = minimum.State();
  //      ROOT::Minuit2::MinimumParameters minParams = minState.Parameters();
    
  //    cerr << endl << endl << "EDM=" << minimum.Edm() << endl;
  //    cerr << "MinValue=" <<  minimum.Fval() << endl << endl;
      
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
        
  if( !minimum.IsValid() )
  {
      //XXX - should we try to re-fit here? Or do something to handle the
      //      faliure in some reasonable way?
      cerr << endl << endl << "status is not valid"
            << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
            << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
            << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
            << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
            << "\n\tHasValidParameters: " << minimum.HasValidParameters()
            << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
            << endl;
      if( minimum.IsAboveMaxEdm() )
        cout << "\t\tEDM=" << minimum.Edm() << endl;
  }//if( !minimum.IsValid() )
      
  
  vector<double> fitpars = fitParams.Params();
  vector<double> fiterrors = fitParams.Errors();
  chi2Fcn.parametersToPeaks( fitPeaks, &fitpars[0], &fiterrors[0] );
      
  double initialChi2 = chi2Fcn.chi2( &fitpars[0] );
  
  //Lets try to keep whether or not to fit parameters should be the same for
  //  the output peaks as the input peaks.
  //Note that this doesnt account for peaks swapping with eachother in the fit
  assert( fitPeaks.size() == inputPeaks.size() );
  
  //for( size_t i = 0; i < near_peaks.size(); ++i )
  //  fitpeaks[i].inheritUserSelectedOptions( near_peaks[i], true );
  //for( size_t i = 0; i < fixedpeaks.size(); ++i )
  //  fitpeaks[i+near_peaks.size()].inheritUserSelectedOptions( fixedpeaks[i], true );
        
  const double totalNDF = set_chi2_dof( spec, fitPeaks, 0, fitPeaks.size() );
  
  chi2 = initialChi2;
  numDOF = static_cast<int>( std::round(totalNDF) );
  for( const auto &peak : fitPeaks )
    peaks.push_back( make_shared<PeakDef>( peak ) );
}//void DetectionConfidenceTool::computeForAcivity(...)




void DetectionConfidenceTool::setRefLinesAndGetLineInfo()
{
  // \TODO: this function is essentually the same as #ReferencePhotopeakDisplay::updateDisplayChange
  //        and given its un-cleaness, it may be worth refactoring.
  auto spec = m_chart->data();
  if( !m_our_meas || !spec )
    return;
    
  std::shared_ptr<DetectorPeakResponse> drf = m_our_meas->detector();
  if( !drf || !drf->hasResolutionInfo() || !drf->isValid() )
    return;
  
  
  const string isotxt = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotxt );

  if( !nuc )
    return;
    
  const string disttxt = m_distanceEdit->text().toUTF8();
  double distance = 0.0;
  try
  {
    distance = PhysicalUnits::stringToDistance(disttxt);
    
    if( distance <= 0.0 )
      throw runtime_error( "Distance cant be zero or negative" );
  }catch( std::exception & )
  {
    return;
  }//try / catch
  
  
  const double brCutoff = 0.0;
  double shielding_an = 0.0, shielding_ad = 0.0, shielding_thickness = 0.0;
  std::shared_ptr<Material> shielding_material;
  const bool generic_shielding = m_shieldingSelect->isGenericMaterial();
  if( generic_shielding )
  {
    shielding_an = m_shieldingSelect->atomicNumber();
    shielding_ad = m_shieldingSelect->arealDensity();
  }else
  {
    shielding_material = m_shieldingSelect->material();
    shielding_thickness = m_shieldingSelect->thickness();
  }//if( generic shielding ) / else
  
  string agestr;
  double age = PeakDef::defaultDecayTime( nuc, &agestr );
  
  auto theme = m_interspec->getColorTheme();
  
  ReferenceLineInfo reflines;
  reflines.nuclide = nuc;
  if( theme && theme->referenceLineColor.size() )
    reflines.lineColor = theme->referenceLineColor[0];
  
  reflines.showGammas = true;
  reflines.showXrays = false;
  reflines.showAlphas = false;
  reflines.showBetas = false;
  reflines.showLines = true;
  reflines.promptLinesOnly = false;
  reflines.isBackground = false;
  reflines.isReaction = false;
  reflines.displayLines = true;
  reflines.age = age;
  reflines.lowerBrCuttoff = 0.0;
  
  reflines.labelTxt = isotxt;
  
  if( generic_shielding )
  {
    const double an = m_shieldingSelect->atomicNumber();
    const double ad = m_shieldingSelect->arealDensity();
    reflines.shieldingName = "AN=" + std::to_string(an) + ", AD=" + std::to_string(ad) + " g.cm2";
  }else if( shielding_material )
  {
    reflines.shieldingName = shielding_material->name;
  }
  reflines.shieldingThickness = shielding_thickness;
  reflines.detectorName = drf->name();
  
  
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( nuc, 1.0E-3 * SandiaDecay::curie );
  
  const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( age );
  const double parent_activity = nuc ? mixture.activity( age, nuc ) : 0.0;
  
  vector<double> energies, branchratios;
  vector<SandiaDecay::ProductType> particle_type;
  vector<const SandiaDecay::Transition *> transistions;
  vector<SandiaDecay::ProductType> types;
  
  
  
  types.push_back( SandiaDecay::GammaParticle );
  types.push_back( SandiaDecay::PositronParticle );
  types.push_back( SandiaDecay::XrayParticle );
  
  
  const float positron_energy = static_cast<float>( 510.9989 * PhysicalUnits::keV );
  int positiron_decayMode = 0;
  float positron_branchRatio = 0.0f;
  const SandiaDecay::Nuclide *positronNuc = 0;
  
  
  std::set<const SandiaDecay::Nuclide *> positronparents;
  std::set<const SandiaDecay::Transition *> positrontrans;
  
  for( SandiaDecay::ProductType type : types )
  {
    for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
    {
      const SandiaDecay::Nuclide *nuclide = activities[nucIndex].nuclide;
      const double activity = activities[nucIndex].activity;
      
      const size_t n_decaysToChildren = nuclide->decaysToChildren.size();
      
      for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
      {
        const SandiaDecay::Transition *transition
        = nuclide->decaysToChildren[decayIndex];
        const size_t n_products = transition->products.size();
        
        for( size_t productNum = 0; productNum < n_products; ++productNum )
        {
          const SandiaDecay::RadParticle &particle
          = transition->products[productNum];
          if( type == SandiaDecay::PositronParticle && particle.type == SandiaDecay::PositronParticle )
          {
            const double br = activity * particle.intensity
            * transition->branchRatio / parent_activity;
            positron_branchRatio += 2.0*br;
            positiron_decayMode   = transition->mode;
            positronparents.insert( transition->parent );
            positrontrans.insert( transition );
          }else if( (particle.type == type) && (particle.type == SandiaDecay::XrayParticle) )
          {
            size_t index = 0;
            for( ; index < energies.size(); ++index )
            {
              if( fabs(energies[index] - particle.energy) < 1.0E-6 )
                break;
            }
            
            const double br = activity * particle.intensity
            * transition->branchRatio / parent_activity;
            
            if( index < energies.size() )
            {
              branchratios[index] += br;
            }else
            {
              transistions.push_back( NULL );
              energies.push_back( particle.energy );
              branchratios.push_back( br );
              particle_type.push_back( SandiaDecay::XrayParticle );
            }
          }else if( particle.type == type )
          {
            transistions.push_back( transition );
            energies.push_back( particle.energy );
            const double br = activity * particle.intensity
            * transition->branchRatio / parent_activity;
            branchratios.push_back( br );
            particle_type.push_back( type );
          }//if( particle.type == type )
        }//for( size_t productNum = 0; productNum < n_products; ++productNum )
      }//for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
    }//for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
  }//for( SandiaDecay::ProductType type : types )
  
  if( positron_branchRatio > 0.0 )
  {
    if( positronparents.size() == 1 )
      positronNuc = *positronparents.begin();
    
    if( positrontrans.size() == 1 )
      transistions.push_back( *positrontrans.begin() );
    else
      transistions.push_back( NULL );
    energies.push_back( positron_energy );
    branchratios.push_back( positron_branchRatio );
    particle_type.push_back( SandiaDecay::GammaParticle );
  }//if( positronrow.branchRatio > 0.0 )
  
  
  
  //fold in detector response
  std::shared_ptr<DetectorPeakResponse> det = m_detectorDisplay->detector();
  
  //Wider peaks mean not as large value of 'y' for the peaks
  if( det && det->isValid() )
  {
    for( size_t i = 0; i < branchratios.size(); ++i )
      if( (particle_type[i] == SandiaDecay::GammaParticle) || (particle_type[i] == SandiaDecay::XrayParticle) )
        branchratios[i] *= det->efficiency( energies[i], distance );
  }//if( detector )
  
  
  //fold in shielding here....
  try
  {
    std::shared_ptr<const Material> material;
    boost::function<double(float)> att_coef_fcn;
    
    if( m_shieldingSelect->isGenericMaterial() )
    {
      const float atomic_number = static_cast<float>(m_shieldingSelect->atomicNumber());
      const float areal_density = static_cast<float>(m_shieldingSelect->arealDensity());
      att_coef_fcn = boost::bind( &GammaInteractionCalc::transmition_coefficient_generic,
                                 atomic_number, areal_density, _1 );
    }else
    {
      material = m_shieldingSelect->material();
      if( !!material )
      {
        const float thick = static_cast<float>(m_shieldingSelect->thickness());
        att_coef_fcn
        = boost::bind( &GammaInteractionCalc::transmition_coefficient_material,
                      material.get(), _1, thick );
      }//if( !!material )
    }//if( isGenericMaterial ) / else
    
    if( !att_coef_fcn.empty() )
    {
      for( size_t i = 0; i < branchratios.size(); ++i )
        if( (particle_type[i] == SandiaDecay::GammaParticle) || (particle_type[i] == SandiaDecay::XrayParticle) )
          branchratios[i] *= exp( -1.0 * att_coef_fcn( energies[i] ) );
    }//if( att_coef_fcn )
  }catch( MassAttenuation::ErrorLoadingDataException & )
  {
    throw runtime_error( "Failed to open gamma XS data file" );
  }catch( std::exception &e )
  {
    cerr << "ReferencePhotopeakDisplay::updateDisplayChange(): caught error " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
    char msg[512];
    snprintf( msg, sizeof(msg), "Error caclulating attenuation: %s", e.what() );
    log_developer_error( BOOST_CURRENT_FUNCTION, msg );
#endif
  }
  
  
  //Peak height is: area*(1/(sigma*sqrt(2*pi)))*exp( -(x-mean)^2 / (2*sigma^2) ),
  //  therefore peak height is proportianal to area/sigma, lets correct for this
  /*
  if( det && det->isValid() && det->hasResolutionInfo() )
  {
    const vector<double> origbr = branchratios;
    
    try
    {
      for( size_t i = 0; i < branchratios.size(); ++i )
      {
        const double sigma = det->peakResolutionSigma( energies[i] );
        if( sigma <= 0.0 )
          throw exception();
        branchratios[i] /= sigma;
      }//for( size_t i = 0; i < branchratios.size(); ++i )
    }catch(...)
    {
      branchratios = origbr;
      cerr << "Encountered a negative or zero peak width, not taking detector "
      << "resolution into account, sorry :(" << endl;
    }//try / catch
  }//if( detector->hasResolutionInfo() )
   */
  
  //Some decays may not produce gammas, but do produce xrays (not verified) so
  //  we want to normalize gammas and xrays relative to the largest branching
  //  ratio gamma or xray.
  double max_gamma_xray = 0.0;
  
  map<SandiaDecay::ProductType,double> maxbrs;
  
  for( size_t i = 0; i < branchratios.size(); ++i )
  {
    const SandiaDecay::ProductType type = particle_type[i];
    
    if( !maxbrs.count(type) )
      maxbrs[type] = 0.0;
    //    if( (transistions[i]
    //         || inforows[i].decayMode==DecayParticleModel::RowData::ReactionToGammaMode) )
    maxbrs[type] = max(maxbrs[type],branchratios[i]);
    
    if( type == SandiaDecay::GammaParticle || type == SandiaDecay::XrayParticle )
      max_gamma_xray = std::max( max_gamma_xray, branchratios[i] );
  }//for( size_t i = 0; i < branchratios.size(); ++i )
  
  for( size_t i = 0; i < branchratios.size(); ++i )
  {
    const double energy = energies[i];
    
    //If this is an xray caused by a decay, lets normalize its amplitude relative
    //  to the gamma amplitudes.  If we are displaying just the xrays of an
    //  element, than we will normalize them to go between zero and one.
    const bool is_gamma = (particle_type[i] == SandiaDecay::GammaParticle);
    const bool is_xray = (particle_type[i] == SandiaDecay::XrayParticle);
    const bool is_decay_xray_gamma = (nuc && (is_gamma || is_xray));
    const double br = branchratios[i] / (is_decay_xray_gamma ? max_gamma_xray : maxbrs[particle_type[i]]);
    
    const SandiaDecay::Transition *transition = transistions[i];
    
    if( transition && (br <= brCutoff || IsInf(br) || IsNan(br)) )
      continue;
    
    string particlestr, decaystr, elstr;
    if( !is_xray )
    {
      const SandiaDecay::ProductType parttype
      = SandiaDecay::ProductType( particle_type[i] );
      particlestr = SandiaDecay::to_str( parttype );
      
      if( transition )
      {
        if( transition->parent )
          decaystr = transition->parent->symbol;
        if( transition->child )
          decaystr += " to " + transition->child->symbol;
        decaystr += string(" via ") + SandiaDecay::to_str(transition->mode);
      }//if( transition ) / else ...
    }else
    {
      const SandiaDecay::Element *element = db->element( nuc->atomicNumber );
      
      particlestr = "xray";
      decaystr = "xray";
      if( element )
        elstr = element->name;
    }//if( xray ) / else
        
    reflines.energies.push_back(     energy );
    reflines.intensities.push_back(  br );
    reflines.particlestrs.push_back( particlestr );
    reflines.decaystrs.push_back(    decaystr );
    reflines.elementstrs.push_back(  elstr );
  }//for( size_t i = 0; i < branchratios.size(); ++i )
  
  typedef std::map<SandiaDecay::ProductType,double>::const_iterator MaxBrIter;
  
  for( MaxBrIter iter = maxbrs.begin(); iter != maxbrs.end(); ++iter )
  {
    const char *typestr = SandiaDecay::to_str( SandiaDecay::ProductType(iter->first) );
    reflines.particle_sf[typestr] = iter->second;
    
    if( iter->first == SandiaDecay::GammaParticle || iter->first == SandiaDecay::XrayParticle )
      reflines.particle_sf[typestr] = max_gamma_xray;
  }
  
  
  m_chart->setReferncePhotoPeakLines( reflines );
}//std::vector<std::tuple<float,double,bool>> setRefLinesAndGetLineInfo();






void DetectionConfidenceTool::do_development()
{
  auto specfile = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !specfile )
  {
    new WLabel( "No Foreground File", this );
    return;
  }
  
  auto spec = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !spec )
  {
    new WLabel( "No foreground spectrum", this );
    return;
  }
  
  std::shared_ptr<DetectorPeakResponse> drf = specfile->detector();
  if( !drf )
  {
    new WLabel( "No DRF loaded", this );
    return;
  }
  
  if( !drf->hasResolutionInfo() || !drf->isValid() )
  {
    new WLabel( "DRF not valid or doesnt have FWHM info.", this );
    return;
  }
  
  
  vector<PeakDef> inputPeaks, fitPeaks;
  
  //vector<float> energies = { 80.9971, 276.4, 302.853, 356.017, 383.848 };
  vector<float> energies = { 356.017 };
  //vector<tuple<float,float>> peak_info = { {356.017f,0.0f} };
  
  //We provide <energy,gammas_per_uci> and optionally activity, and get out chi2, DOF, and activity if not provided
  //First, fit for preffered activity.
  //Then increase the activity until the Chi2 increases by quantile(chi_squared(n_pars),0.95)
  
  
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const float mean = energies[i];
    const float sigma = drf->peakResolutionSigma(mean);
    const float amplitude = 0.0f;
    PeakDef peak( mean, sigma, amplitude );
    peak.setFitFor( PeakDef::CoefficientType::Mean, false );
    peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
    peak.setFitFor( PeakDef::CoefficientType::GaussAmplitude, true );
    shared_ptr<PeakContinuum> cont = peak.continuum();
    cont->setType( PeakContinuum::OffsetType::Linear );
    
    double lowerEnengy, upperEnergy;
    findROIEnergyLimits( lowerEnengy, upperEnergy, peak, spec );
    cont->setRange( lowerEnengy, upperEnergy );
    cont->calc_linear_continuum_eqn( spec, lowerEnengy, upperEnergy, 2 );
    cont->setPolynomialCoefFitFor( 0, true );
    cont->setPolynomialCoefFitFor( 1, true );
    
    inputPeaks.push_back( std::move(peak) );
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  ROOT::Minuit2::MnUserParameters inputFitPars;
  PeakFitChi2Fcn::addPeaksToFitter( inputFitPars, inputPeaks, spec, PeakFitChi2Fcn::kFitForPeakParameters );
  
      
  const int npeaks = static_cast<int>( inputPeaks.size() );
  PeakFitChi2Fcn chi2Fcn( npeaks, spec, nullptr );
  chi2Fcn.useReducedChi2( false );
      
  assert( inputFitPars.VariableParameters() != 0 );
          
  ROOT::Minuit2::MnUserParameterState inputParamState( inputFitPars );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( chi2Fcn, inputParamState, strategy );
      
  unsigned int maxFcnCall = 5000;
  double tolerance = 2.5;
  tolerance = 0.5;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  const ROOT::Minuit2::MnUserParameters &fitParams = fitter.Parameters();
  //  minimum.IsValid()
  //      ROOT::Minuit2::MinimumState minState = minimum.State();
  //      ROOT::Minuit2::MinimumParameters minParams = minState.Parameters();
    
  //    cerr << endl << endl << "EDM=" << minimum.Edm() << endl;
  //    cerr << "MinValue=" <<  minimum.Fval() << endl << endl;
      
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
        
  if( !minimum.IsValid() )
  {
      //XXX - should we try to re-fit here? Or do something to handle the
      //      faliure in some reasonable way?
      cerr << endl << endl << "status is not valid"
            << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
            << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
            << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
            << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
            << "\n\tHasValidParameters: " << minimum.HasValidParameters()
            << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
            << endl;
      if( minimum.IsAboveMaxEdm() )
        cout << "\t\tEDM=" << minimum.Edm() << endl;
  }//if( !minimum.IsValid() )
      
  
  vector<double> fitpars = fitParams.Params();
  vector<double> fiterrors = fitParams.Errors();
  chi2Fcn.parametersToPeaks( fitPeaks, &fitpars[0], &fiterrors[0] );
      
  double initialChi2 = chi2Fcn.chi2( &fitpars[0] );
  
  //Lets try to keep whether or not to fit parameters should be the same for
  //  the output peaks as the input peaks.
  //Note that this doesnt account for peaks swapping with eachother in the fit
  assert( fitPeaks.size() == inputPeaks.size() );
  
  //for( size_t i = 0; i < near_peaks.size(); ++i )
  //  fitpeaks[i].inheritUserSelectedOptions( near_peaks[i], true );
  //for( size_t i = 0; i < fixedpeaks.size(); ++i )
  //  fitpeaks[i+near_peaks.size()].inheritUserSelectedOptions( fixedpeaks[i], true );
        
  set_chi2_dof( spec, fitPeaks, 0, fitPeaks.size() );
      
  
  auto label = new WLabel( "Initial Chi2=" + std::to_string(initialChi2), this );
  cout << "Intiial Chi2=" << initialChi2 << endl;
  label->setInline( false );
  
  using boost::math::chi_squared;
  using boost::math::quantile;
  using boost::math::complement;
  using boost::math::cdf;
  
  
  std::sort( fitPeaks.begin(), fitPeaks.end(), &PeakDef::lessThanByMean );
  for( auto &peak : fitPeaks )
  {
    const double amp = peak.amplitude();
    string msg = "Peak at " + std::to_string(peak.mean()) + " keV"
                 + " fit amplitude: " + std::to_string(amp)
                 + " and Chi2/dof=" + std::to_string(peak.chi2dof());
    auto label = new WLabel( msg, this );
    label->setInline( false );
    cout << msg << endl;
    
    std::vector<PeakDef *> peakptrs;
    for( auto &peak : fitPeaks )
      peakptrs.push_back( &peak );
      
    double chi2, dof;
    get_chi2_and_dof_for_roi( chi2, dof, spec, peakptrs );
    
    chi_squared dist( dof );
    const double prob = boost::math::cdf(dist,chi2); //Probability we would have seen a chi2 this large.
    double p_value = 1.0 - prob; //Probability we would have observed this good of a chi2, or better
    cout << "There are " << dof << " DOF" << endl;
    cout << "This chi2=" << chi2 << " and prob of observing a Chi2 at least this extreme " << p_value << endl;
    cout << "95% of the time, we would see a Chi2 value of " << quantile(dist, 0.95) << " or less" << endl;
    
    for( double alpha = 0.0; alpha < 1.0; alpha += 0.05 )
    {
      cout << "alpha=" << alpha << " --> quantile=" << quantile(dist, alpha) << endl;
    }
    
    //We are estimating 1 parameter (amplitude), so the amplitude of the peak that we are 95% confident
    //  it would be lower than is the amplitude that causes a chi2 change of 3.84
    //const double n_parameters_being_fit = 1;
    //chi_squared deltadist( n_parameters_being_fit );
    //cout << "quantile(0.95)=" << quantile(deltadist, 0.95) << " (should be 3.84)" << endl;
    
    //We wa
    
  }//for( const auto &peak : fitPeaks )
  
  
  
  
  
  
  
  /*
  vector<PeakDef> *all_peaks = &fitpeaks;
      std::unique_ptr< vector<PeakDef> > all_peaks_ptr;
      if( fixedpeaks.size() )
      {
        all_peaks = new vector<PeakDef>( fitpeaks );
        all_peaks_ptr.reset( all_peaks );
        all_peaks->insert( all_peaks->end(), fixedpeaks.begin(), fixedpeaks.end() );
        std::sort( all_peaks->begin(), all_peaks->end(), &PeakDef::lessThanByMean );
      }//if( fixedpeaks.size() )
      
      for( size_t i = 1; i <= fitpeaks.size(); ++i ) //Note weird convntion of index
      {
        const PeakDef *peak = &(fitpeaks[i-1]);
        const bool significant = chi2_significance_test(
                                                            *peak, stat_threshold, hypothesis_threshold,
                                                            *all_peaks, data );
        if( !significant )
        {
  #if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          DebugLog(cerr) << "\tPeak at mean=" << peak->mean()
                         << "is being discarded for not being significant"
                         << "\n";
  #endif
          fitpeaks.erase( fitpeaks.begin() + --i );
        }//if( !significant )
      }//for( size_t i = 1; i < fitpeaks.size(); ++i )
          
      bool removed_peak = false;
      for( size_t i = 1; i < fitpeaks.size(); ++i ) //Note weird convntion of index
      {
        PeakDef *this_peak = &(fitpeaks[i-1]);
        PeakDef *next_peak = &(fitpeaks[i+1-1]);
            
        const double min_sigma = min( this_peak->sigma(), next_peak->sigma() );
        const double mean_diff = next_peak->mean() - this_peak->mean();
            
        //In order to remove a gaussian, the peaks must both be within a sigma
        //  of eachother.  Note that this proccess doesnt care about the widths
        //  of the gaussians because we are assuming that the width of the gaussian
        //  should only be dependant on energy, so should only have one width of
        //  gaussian for a given energy
        if( (mean_diff/min_sigma) < 1.0 ) //XXX 1.0 chosen arbitrarily, and not checked
        {
  #if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          DebugLog(cerr) << "Removing duplicate peak at x=" << this_peak->mean() << " sigma="
              << this_peak->sigma() << " in favor of mean=" << next_peak->mean()
              << " sigma=" << next_peak->sigma() << "\n";
  #endif
              
          removed_peak = true;
              
          //Delete the peak with the worst chi2
          if( this_peak->chi2dof() > next_peak->chi2dof() )
            fitpeaks.erase( fitpeaks.begin() + i - 1 );
          else
            fitpeaks.erase( fitpeaks.begin() + i );
              
          i = i - 1; //incase we have multiple close peaks in a row
        }//if( (mean_diff/min_sigma) < 1.0 ) / else
      }//for( size_t i = 1; i < fitpeaks.size(); ++i )
          
      if( removed_peak )
        fitPeaks( fitpeaks, stat_threshold, hypothesis_threshold,
                  data, fitpeaks, fixedpeaks, false );
          
      if( datadefined_peaks.size() )
      {
        fitpeaks.insert( fitpeaks.end(),
                        datadefined_peaks.begin(), datadefined_peaks.end() );
        std::sort( fitpeaks.begin(), fitpeaks.end(), &PeakDef::lessThanByMean );
      }
          
      return;
    }catch( std::exception &e )
    {
      cerr << "fitPeaks(...)\n\tSerious programming logic error: caught"
           << " exception where I really shouldnt have.  what()=" << e.what()
           << endl;
    }catch(...)
    {
      cerr << "fitPeaks(...)\n\tSerious programming logic error: caught"
           << " unknown exception where I really shouldnt have." << endl;
    }//try/catch()
          
    //We will only reach here if there was no exception, so since never expect
    //  this to actually happen, just assign the results to be same as the input
    fitpeaks = all_near_peaks;
  }//vector<PeakDef> fitPeaks(...);
   */
}//void do_development()
