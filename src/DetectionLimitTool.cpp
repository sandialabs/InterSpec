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

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WComboBox>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WTableCell>
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
#include "InterSpec/DetectionLimitTool.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

using namespace std;
using namespace Wt;

#define MU_CHARACTER "\xCE\xBC"

namespace
{

void local_eqn_from_offsets( size_t lowchannel,
                                     size_t highchannel,
                                     const double peakMean,
                                     const std::shared_ptr<const SpecUtils::Measurement> &data,
                                     const size_t num_side_channel,
                                     double &m, double &b )
{
  double x1  = data->gamma_channel_lower( lowchannel ) - peakMean;
  const double dx1 = data->gamma_channel_width( lowchannel );
  
  //XXX - hack! below 'c' will be inf if the following check would fail; I
  //      really should work out how to fix this properly with math, but for
  //      right now I'll fudge it, which will toss the answer off a bit, but
  //      whatever.
  if( fabs(2.0*x1*dx1 + dx1*dx1) < DBL_EPSILON )
  {
    static int ntimes = 0;
    if( ntimes++ < 4 )
      cerr << "Should fix math in local_eqn_from_offsets" << endl;
    x1 = data->gamma_channel_lower( lowchannel+1 ) - peakMean;
  }
  
  const double x2  = data->gamma_channel_lower( highchannel ) - peakMean;
  
  const double dx2 = data->gamma_channel_width( highchannel );
  
  if( (2*num_side_channel + 1) > data->num_gamma_channels() )
    throw std::runtime_error( "Too many side channels specified" );
  
  if( lowchannel <= num_side_channel )
    lowchannel = num_side_channel + 1;
  highchannel = std::min( highchannel, data->num_gamma_channels() - 1 - num_side_channel );
  
  if( lowchannel >= highchannel )
    throw runtime_error( "Lower channel is above upper channel" );
  
  const double nbinInv = 1.0 / num_side_channel;
  const double lower_region_sum = data->gamma_channels_sum( lowchannel - num_side_channel, lowchannel - 1 );
  const double upper_region_sum = data->gamma_channels_sum( highchannel + 1, highchannel + num_side_channel);
  const double y1 = nbinInv * lower_region_sum;
  const double y2 = nbinInv * upper_region_sum;
  
  cout << "lowchannel=" << lowchannel << ", highchannel=" << highchannel << ", num_side_channel=" << num_side_channel << endl;
  cout << "Pre lower ROI: ";
  for( auto i = lowchannel - num_side_channel; i <= (lowchannel - 1); ++i )
    cout << i << ", ";
  cout << "  --> " << lower_region_sum << endl;
  cout << "Post upper ROI: ";
  for( auto i = highchannel + 1; i <= (highchannel + num_side_channel); ++i )
    cout << i << ", ";
  cout << "  --> " << upper_region_sum << endl;
  
  
  const double c = (2.0*x2*dx2 + dx2*dx2)/(2.0*x1*dx1 + dx1*dx1);
  b = (y2-y1*c)/(dx2-dx1*c);
  m = 2.0*(y1-b*dx1)/(2.0*x1*dx1+dx1*dx1);
  
  if( IsNan(m) || IsInf(m) || IsNan(b) || IsInf(b) )
  {
    cerr << "local_eqn_from_offsets(...): Invalid results" << endl;
    m = b = 0.0;
  }
}//eqn_from_offsets(...)


class MdaPeakRow : public WContainerWidget
{
public:
  const float m_energy;
  const double m_distance;
  const double m_gammas_per_bq;
  const float m_fwhm;
  std::shared_ptr<DetectorPeakResponse> m_drf;
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
    
    const size_t nsidebin = 5;
    const size_t nchannels = m_meas->num_gamma_channels();
    
    // \TODO: check that we are getting the right channels, and not off by one or something.
    size_t lower_bin = m_meas->find_gamma_channel(roi_lower_energy);
    //lower_bin = ((lower_bin > (2*nsidebin+1)) ? (lower_bin - nsidebin) : (nsidebin + 1));
    
    size_t upper_bin = m_meas->find_gamma_channel(roi_upper_energy);
    //upper_bin = ((upper_bin < (nchannels - 2*nsidebin - 1)) ? (upper_bin + nsidebin) : (nchannels - nsidebin - 1));
    
    lower_bin = std::min( lower_bin, nchannels - nsidebin - 1 );
    upper_bin = std::max( upper_bin, nsidebin + 1 );
    
    double coefs[2];
    local_eqn_from_offsets( lower_bin, upper_bin, m_energy, m_meas, nsidebin, coefs[1], coefs[0] );
    
    //double continuum = PeakContinuum::offset_eqn_integral( coefs,
    //                                   PeakContinuum::OffsetType::Linear,
    //                                   roi_lower_energy, roi_upper_energy, m_energy );
    //const double data_area = m_meas->gamma_integral(roi_lower_energy, roi_upper_energy);
    
    double continuum = PeakContinuum::offset_eqn_integral( coefs,
                                           PeakContinuum::OffsetType::Linear,
                                           m_meas->gamma_channel_lower(lower_bin),
                                           m_meas->gamma_channel_upper(upper_bin),
                                           m_energy );
    const double data_area = m_meas->gamma_channels_sum(lower_bin, upper_bin);
    
    /*
    auto energy_cal = m_meas->energy_calibration();
    const double meanchannel = energy_cal->channel_for_energy(m_energy);
    const double gad_peak_region = 10.39;
    const double gad_lower_channel = meanchannel - 0.5*gad_peak_region;
    const double gad_upper_channel = meanchannel + 0.5*gad_peak_region;
    const double gad_lower_energy = energy_cal->energy_for_channel(gad_lower_channel);
    const double gad_upper_energy = energy_cal->energy_for_channel(gad_upper_channel);
    
    cout << "(gad_upper_energy - gad_lower_energy)/fwhm(" << m_energy << " keV)=" << (gad_upper_energy - gad_lower_energy)/m_drf->peakResolutionFWHM(m_energy) << endl;
    
    cout << "For " << m_energy << " keV, ROI is " << (1 + upper_bin - lower_bin)
    << " channels centered around channel " << meanchannel
    << ". Gadras region would be from " << gad_lower_energy << " to " << gad_upper_energy << " keV"
    << " (channels " << gad_lower_channel << " to " << gad_upper_channel << ")"
    << endl;
    
    cout << "Predicted continuum area=" << continuum << " vs data " << data_area << endl;
    */
    //continuum = data_area;
    
    m_simpleMda = (2.71 + 4.65*sqrt(continuum))/m_gammas_per_bq;
    if( IsInf(m_simpleMda) || IsInf(m_simpleMda) )
      m_simpleMda = -1.0;
    const string mdastr = PhysicalUnits::printToBestActivityUnits(m_simpleMda);
    cout << "mda(" << m_energy << ") = " << mdastr << endl;
    
    m_poisonLimit->setText( "Simple MDA=" + mdastr );
    
    cout << endl << endl;
    cout << "Gross Counts in Peak Foreground: " << data_area << endl;
    //cout << "Gross Counts in Peak Background: 0
    cout << "Gross Counts in Continuum Foreground: " << continuum << endl;
    //Gross Counts in Continuum Background: 0
    //Variance in Peak Region: 9.062885
    //Variance in Continuum Region: 8.400412
    cout << "Detector Efficiency at Energy: " << m_drf->intrinsicEfficiency(m_energy) << endl;
    //Range of Gammas Entering Detector: (0, 0.4153258)
    cout << "Solid Angle: " << m_drf->fractionalSolidAngle(m_drf->detectorDiameter(), m_distance) << endl;
    //Range of Gammas from Source: (0, 1607.631)
    cout << "Gamma Emission Rate per Ci of Source: " << (PhysicalUnits::curie * m_gammas_per_bq / m_meas->live_time())/m_drf->efficiency(m_energy, m_distance) << endl;
    //Transmission through Shielding: 1
    cout << "Number of counts limit: " << m_simpleMda*m_gammas_per_bq << endl;
    
    // Need to propagate continuum uncertainty to get variance in peak region.
    // Make a single-peak limit section that gives some of the information above, and maybe draws the limiting peak.
    // Also need to make limit used settable.
    // Need to add in if it passes the critical limit what it is.
    
    /*
    boost::uintmax_t max_iter = 1000;
    const double mda = m_simpleMda;
    const double givenActivity_ci = 10; //10 Ci
    const double origdist = 1.0;
    auto distMinForActivity = [mda,givenActivity_ci,origdist]( double const &x ) -> double {
      return fabs( ((1.0/(x*x)) * exp(-0.0094276*x)) - (mda/PhysicalUnits::curie)/givenActivity_ci/origdist );
    };
      
    const pair<double, double> r = boost::math::tools::brent_find_minima( distMinForActivity, 1.0, 100000.0, 10, max_iter );
    
    cout << "For energy(" << m_energy << "), detection " << givenActivity_ci << " Ci, dist=" << r.first
         << "m with error " << r.second << endl;
     */
  }//void setSimplePoisonTxt()
  
  void emitChanged()
  {
    setSimplePoisonTxt();
    m_changed.emit();
  }
  
public:
  MdaPeakRow( const float energy, const double gammasPerBq,
             const double distance,
              float roi_start, float roi_end,
              shared_ptr<const SpecUtils::Measurement> meas,
              std::shared_ptr<DetectorPeakResponse> drf,
             Wt::WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
  m_energy( energy ),
  m_distance( distance ),
  m_gammas_per_bq( gammasPerBq ),
  m_fwhm( drf->peakResolutionFWHM(energy) ),
  m_drf( drf ),
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
    
    m_use->checked().connect( this, &MdaPeakRow::emitChanged );
    m_use->unChecked().connect( this, &MdaPeakRow::emitChanged );
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



DetectionLimitWindow::DetectionLimitWindow( InterSpec *viewer,
                                                     MaterialDB *materialDB,
                                                     WSuggestionPopup *materialSuggest )
: AuxWindow( "Detection Confidence Tool",
  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
   | AuxWindowProperties::SetCloseable
   | AuxWindowProperties::DisableCollapse) ),
  m_tool( nullptr )
{
  assert( viewer );
  rejectWhenEscapePressed( true );
  
  m_tool = new DetectionLimitTool( viewer, materialDB, materialSuggest );
  WContainerWidget *content = contents();
  content->addStyleClass( "DetectionLimitWindowContent" );
  content->addWidget( m_tool );
  
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
}//DetectionLimitWindow(...) constrctor


DetectionLimitWindow::~DetectionLimitWindow()
{
}




DetectionLimitTool::DetectionLimitTool( InterSpec *viewer,
                                                  MaterialDB *materialDB,
                                                  Wt::WSuggestionPopup *materialSuggest,
                                                  WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_interspec( viewer ),
    m_needsUpdate( true ),
    m_chart( nullptr ),
    m_peakModel( nullptr ),
    m_nuclideEdit( nullptr ),
    m_currentNuclide( nullptr ),
    m_nuclideSuggest( nullptr ),
    m_detectorDisplay( nullptr ),
    m_distanceEdit( nullptr ),
    m_materialDB( materialDB ),
    m_materialSuggest( materialSuggest ),
    m_shieldingSelect( nullptr ),
    m_displayActivity( nullptr ),
    m_results( nullptr ),
    m_chi2Chart( nullptr ),
    m_bestChi2Act( nullptr ),
    m_upperLimit( nullptr ),
    m_errorMsg( nullptr )
{
  /** TODO:
   - [ ] Make test cases that will quickly iterate through, to test things
   - [ ] Allow user to choose activity limit, or distance limit
   - [ ] Make a "by eye" equivalent CL
   - [ ] All adding in a scale factor, so if spectrum is of 30 minutes, allow making limit for 30 second spectra
   - [x] Make the Chi2 plot a D3 based plot
   - [x] Put the Chi2 chart to the right of the spectrum, when it should exist
   - [ ] Give the user a choice about using continuum fixed at null hypothesis
   - [ ] Allow combining ROI with neghboring peaks
   - [x] Make it so when user change ROI on chart, it updates the text input
   - [ ] Allow minor gamma-lines overlapping with primary gamma lines to contribute to peak area
   - [ ] Allow user to pick Currie limit ranges, and improve clarity of this stuff, like maybe have each peak be a WPanel or something
   - [ ] For each row show plot of current peak in that row
   - [ ] Have the energy rows fold down to show more information, similar to Steves tool, for each energy
   - [ ] Allow users to select CL, not just 95%
   - [ ] Add in allowing to age nuclide (didnt I generalize inputting a nuclide somewhere?  Hopefully just re-use that)
   - [ ] Allow users to double click on the spectrum to add a peak to the limit, or similarly for erasing a peak
   - [ ] If user clicks on a result row, have chart zoom to that general region
   - [x] Default fill in reference lines/shielding/age as user has in Reference PhotoPeak tool
   - [ ] Have an explicit option for including attenuation in air
   
  
   Add in options for "limit live time", Distance/Activity, Confidence Limit to use, nuclide age
   
   
   To calculate the Currie
   For each Peak, fit the most likely area, and and its limits.
   If lower confidence is above zero, assume there is a peak
   */
  WApplication * const app = wApp;
  assert( app );
  
  app->useStyleSheet( "InterSpec_resources/DetectionLimitTool.css" );
  app->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  app->require( "InterSpec_resources/DetectionLimitTool.js" );
  
  addStyleClass( "DetectionLimitTool" );
  
  
  //new WLabel( "DetectionLimitTool", this );
  const WLength labelWidth(3.5,WLength::FontEm), fieldWidth(4,WLength::FontEm);
  const WLength optionWidth(5.25,WLength::FontEm), buttonWidth(5.25,WLength::FontEm);
  
  // Container to hold and m_chart and m_results (which itself holds m_chi2Chart,  m_bestChi2Act,
  //  and m_upperLimit).
  WContainerWidget *upperCharts = new WContainerWidget( this );
  upperCharts->addStyleClass( "MdaCharts" );
  

  m_chart = new D3SpectrumDisplayDiv( upperCharts );
  
  /*
  //const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  //const int width = ((screenW < 906) ? screenW : 0.75*screenW ); // maybe cap width off at 1800 px?
  const int height = ((screenH < 525) ? screenH : 0.8*screenH );
  //m_chart->resize( WLength(width-20,WLength::Unit::Pixel), WLength(std::min(0.3*height,300.0),WLength::Unit::Pixel) );
  //m_chart->resize( WLength::Auto, WLength(300,WLength::Unit::Pixel) );
  cout << "height-->" << 1.0*height << endl;
  m_chart->setMinimumSize( WLength::Auto, WLength( std::max(200.0,0.4*height),WLength::Unit::Pixel) );
  */
  
  /** \TODO: Right now the user can drag the ROI around, and if the peak amplitude becaomes
             insignificant or under an area of 5, then the old peak is used, just with ROI range
             altered.  So what we'll do is just catch when the peak ROI drag is final, and then
             adjust the ROI, and re-due our fit - we should probably do something more intuitive
             so the user knows whats going on.
   */
  m_chart->roiDragUpdate().connect( boost::bind( &DetectionLimitTool::roiDraggedCallback, this, _1, _2, _3, _4, _5, _6 ) );
  
  m_chart->setCompactAxis( true );
  m_chart->setXAxisTitle( "Energy (keV)" );
  m_chart->setYAxisTitle( "Counts/Channel" );
  m_chart->disableLegend();
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
  
  
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  
  
  WTable *inputTable = new WTable( this );
  inputTable->addStyleClass( "UserInputTable" );
  
  
  WTableCell *cell = inputTable->elementAt(0, 0);
  WLabel *label = new WLabel( "Distance:", cell );
  label->setMinimumSize( labelWidth, WLength::Auto );
  
  cell = inputTable->elementAt(0, 1);
  m_distanceEdit = new WLineEdit( "100 cm", cell );
  //m_distanceEdit->setTextSize( 5 );
  m_distanceEdit->setValidator( distValidator );
  m_distanceEdit->setMinimumSize( fieldWidth, WLength::Auto );
  label->setBuddy( m_distanceEdit );
  m_distanceEdit->changed().connect( this, &DetectionLimitTool::handleInputChange );
  m_distanceEdit->enterPressed().connect( this, &DetectionLimitTool::handleInputChange );
  
  cell = inputTable->elementAt(1, 0);
  label = new WLabel( "Nuclide:", cell );
  label->setMinimumSize( labelWidth, WLength::Auto );
  
  cell = inputTable->elementAt(1, 1);
  m_nuclideEdit = new WLineEdit( "", cell );
  m_nuclideEdit->setMargin( 1 );
  m_nuclideEdit->setMinimumSize( fieldWidth, WLength::Auto );
  m_nuclideEdit->setAutoComplete( false );
  m_nuclideEdit->changed().connect( boost::bind( &DetectionLimitTool::handleNuclideChange, this ) );
  label->setBuddy( m_nuclideEdit );

  
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
  
  // \TODO: allow entering age, and see ReferencePhotopeakDisplay::handleIsotopeChange( const bool useCurrentAge )
  
  
  isoSuggestModel->filter( "" );
  m_nuclideSuggest->setFilterLength( -1 );
  m_nuclideSuggest->setModel( isoSuggestModel );
  m_nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
  m_nuclideSuggest->forEdit( m_nuclideEdit, WSuggestionPopup::Editing );  // | WSuggestionPopup::DropDownIcon
  
  SpectraFileModel *specFileModel = viewer->fileManager()->model();
  m_detectorDisplay = new DetectorDisplay( viewer, specFileModel );
  viewer->detectorChanged().connect( boost::bind( &DetectionLimitTool::handleInputChange, this ) );
  viewer->detectorModified().connect( boost::bind( &DetectionLimitTool::handleInputChange, this ) );
  
  cell = inputTable->elementAt(0, 3);
  cell->setRowSpan( 2 );
  cell->addWidget( m_detectorDisplay );

  
  cell = inputTable->elementAt(0, 4);
  label = new WLabel( "Peaks disp. act.:", cell );
  cell = inputTable->elementAt(0, 5);
  m_displayActivity = new WLineEdit( cell );
  WRegExpValidator *actvalidator = new WRegExpValidator( PhysicalUnits::sm_activityRegex, m_displayActivity );
  actvalidator->setFlags(Wt::MatchCaseInsensitive);
  m_displayActivity->setValidator( actvalidator );
  m_displayActivity->setTextSize( 10 );
  m_displayActivity->setText( "0 uCi" );
  m_displayActivity->enterPressed().connect( this, &DetectionLimitTool::updateShownPeaks );
  m_displayActivity->changed().connect( this, &DetectionLimitTool::updateShownPeaks );
  
  
  
  m_shieldingSelect = new ShieldingSelect( m_materialDB, NULL, m_materialSuggest, false );
  m_shieldingSelect->materialEdit()->setEmptyText( "<shielding material>" );
  m_shieldingSelect->materialChanged().connect( this, &DetectionLimitTool::handleInputChange );
  m_shieldingSelect->materialModified().connect( this, &DetectionLimitTool::handleInputChange );
  m_shieldingSelect->setMinimumSize( WLength(250), WLength::Auto );
  
  cell = inputTable->elementAt(0,2);
  cell->addWidget( m_shieldingSelect );
  cell->setRowSpan( 2 );
  
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

  m_results = new WContainerWidget( upperCharts );
  m_results->addStyleClass( "MdaResults" );
  m_results->hide();
  
  m_chi2Chart = new WContainerWidget( m_results );
  m_chi2Chart->addStyleClass( "MdaChi2Chart" );
  
  //////////////////////////////////////////////////////////////////////
  m_bestChi2Act = new WText( "&nbsp;", m_results );
  m_bestChi2Act->setInline( false );
  m_bestChi2Act->addStyleClass( "MdaResultTxt" );
    
  m_upperLimit = new WText( "&nbsp;", m_results );
  m_upperLimit->setInline( false );
  m_upperLimit->addStyleClass( "MdaResultTxt" );
  
  
  m_errorMsg = new WText("&nbsp;", this );
  m_errorMsg->addStyleClass( "MdaErrMsg" );
  m_errorMsg->hide();
  
  m_peaks = new WContainerWidget( this );
  m_peaks->addStyleClass( "MdaPeaks" );
  
  
  ReferencePhotopeakDisplay *reflines = viewer->referenceLinesWidget();
  if( reflines )
  {
    const ReferenceLineInfo &current = reflines->currentlyShowingNuclide();
    if( current.nuclide )
      m_nuclideEdit->setText( current.nuclide->symbol );
    
    // \TODO: age
    
    const ShieldingSelect *shielding = reflines->shieldingSelect();
    if( shielding )
    {
      const bool isGeneric = shielding->isGenericMaterial();
      if( isGeneric && (shielding->arealDensity() > 0.0) )
      {
        const double an = shielding->atomicNumber();
        const double ad = shielding->arealDensity();
        m_shieldingSelect->setAtomicNumberAndArealDensity( an, ad );
      }
      
      if( !isGeneric && shielding->materialEdit() && shielding->thicknessEdit() )
      {
        const string material = shielding->materialEdit()->text().toUTF8();
        const string thick = shielding->thicknessEdit()->text().toUTF8();
        
        if( !material.empty() && !thick.empty() )
          m_shieldingSelect->setMaterialNameAndThickness( material, thick );
      }//if( generic ) / else
    }//if( shielding )
  }//if( reflines )
  
  initChi2Chart();
  handleNuclideChange();  //JIC we picked up the one from RefLines
  handleInputChange();
}//DetectionLimitTool constructor
  

  
DetectionLimitTool::~DetectionLimitTool()
{
}


void DetectionLimitTool::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  if( renderFull || m_needsUpdate )
    doCalc();
  
  m_needsUpdate = false;
  
  WContainerWidget::render( flags );
}//render( flags )


void DetectionLimitTool::initChi2Chart()
{
  m_chi2Chart->setJavaScriptMember( "chart", "new MdaChi2Chart(" + m_chi2Chart->jsRef() + ", {});");
  const string jsgraph = m_chi2Chart->jsRef() + ".chart";
  
  m_chi2Chart->setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + m_chi2Chart->id() + "') )"
          + jsgraph + ".redraw();"
        "}"
      "});"
  );
  
  m_chi2Chart->callJavaScriptMember( "resizeObserver.observe", m_chi2Chart->jsRef() );
}



void DetectionLimitTool::roiDraggedCallback( double new_roi_lower_energy,
                                                  double new_roi_upper_energy,
                                                  double new_roi_lower_px,
                                                  double new_roi_upper_px,
                                                  double original_lower_energy,
                                                  bool is_final_range )
{
  if( !is_final_range )
    return;
  
  if( !m_peakModel )  //Shouldnt ever happen
    return;
  
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> origpeaks = m_peakModel->peaks();
  if( !origpeaks || origpeaks->empty() )
    return;  //shouldnt ever happen
  
  
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    
    const float this_roi_start = rw->m_roi_start->value();
    if( fabs(this_roi_start - original_lower_energy) < 0.1 )
    {
      rw->m_roi_start->setValue( new_roi_lower_energy );
      rw->m_roi_end->setValue( new_roi_upper_energy );
      rw->emitChanged();
      
      cout << "Updated ROI of peak at " << rw->m_energy << " keV to go from "
           << new_roi_lower_energy << " to " << new_roi_upper_energy << " keV" << endl;
      return;
    }
  }//for( auto w : m_peaks->children() )
  
  #if( PERFORM_DEVELOPER_CHECKS )
    char msg[512];
    snprintf( msg, sizeof(msg), "Failed to find a ROI that started at %f keV", new_roi_lower_energy );
    log_developer_error( BOOST_CURRENT_FUNCTION, msg );
  #endif

  cerr << "Failed to find a ROI that started at " << new_roi_lower_energy <<  " keV" << endl;
}//void roiDraggedCallback(...)

void DetectionLimitTool::handleNuclideChange()
{
  const string isotxt = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotxt );

  if( nuc == m_currentNuclide )
    return;
  
  m_currentNuclide = nuc;
  handleInputChange();
}//void handleNuclideChange()


void DetectionLimitTool::handleInputChange()
{
  vector<tuple<float,bool,float,float>> oldvalues;
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    oldvalues.push_back( {rw->m_energy, rw->m_use->isChecked(), rw->m_roi_start->value(), rw->m_roi_end->value()} );
  }
  
  
  m_peaks->clear();
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
  

  if( !m_currentNuclide )
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
  double age = PeakDef::defaultDecayTime( m_currentNuclide, &agestr );
  
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
  mixture.addNuclideByActivity( m_currentNuclide, 1.0E-3 * SandiaDecay::curie );

  const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( age );
  const double parent_activity = mixture.activity( age, m_currentNuclide );
  
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
    float roi_start = energy - 1.125*fwhm;
    float roi_end = energy + 1.125*fwhm;
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
    
    
    auto row = new MdaPeakRow( energy, countsPerBq, distance, roi_start, roi_end, spec, drf, m_peaks );
    row->m_use->setChecked(use);
    row->m_changed.connect( this, &DetectionLimitTool::scheduleCalcUpdate );
  }//for( const auto &line : lines )
  
  scheduleCalcUpdate();
}//void handleInputChange()


void DetectionLimitTool::scheduleCalcUpdate()
{
  m_needsUpdate = true;
  scheduleRender();
}//void scheduleCalcUpdate()


void DetectionLimitTool::doCalc()
{
  //m_simpleMda
  double minSearchActivity = std::numeric_limits<double>::infinity(), maxSearchActivity = 0.0;
  
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
      minSearchActivity = std::min( minSearchActivity, rw->m_simpleMda );
      maxSearchActivity = std::max( maxSearchActivity, rw->m_simpleMda );
    }
  }//for( auto w : m_peaks->children() )
  
  if( !nused )
  {
    m_results->hide();
    if( m_errorMsg->isHidden() || m_errorMsg->text().empty() )
    {
      m_errorMsg->setText( "No peaks are selected" );
      m_errorMsg->show();
    }
    
    return;
  }//if( !nused )
  
  if( IsInf(minSearchActivity) || maxSearchActivity==0.0 )
  {
    minSearchActivity = 0.0;
    maxSearchActivity = PhysicalUnits::curie;
  }else
  {
    //minSearchActivity *= 0.25;
    minSearchActivity = 0.0;
    maxSearchActivity *= 10.0;
  }
  
  m_errorMsg->hide();
  
  const double yrange = 15;
  const double cl_chi2_delta = 4.0; //95.45% CL, or something
  
  const size_t nchi2 = 25;  //approx num chi2 to compute
  vector<pair<double,double>> chi2s;
  double overallBestChi2 = 0.0, overallBestActivity = 0.0, upperLimit = 0.0, lowerLimit = 0.0, activityRangeMin = 0.0, activityRangeMax = 0.0;
  bool foundUpperCl = false, foundUpperDisplay = false, foundLowerCl = false, foundLowerDisplay = false;
  
  /// \TODO: currently all this stuff assumes a smooth continuosly increasing Chi2 with increasing
  ///        activity, but this doesnt have to be the case, especially with quadratic continuums.
  
  try
  {
    size_t num_iterations = 0;
    
    //The boost::math::tools::bisect(...) function will make calls using the same value of activity,
    //  so we will cache those values to save some time.
    map<double,double> chi2cache;
    auto chi2ForAct = [this,&num_iterations,&chi2cache]( double const &activity ) -> double {
      
      const auto pos = chi2cache.find(activity);
      if( pos != end(chi2cache) )
        return pos->second;
      
      int numDOF = 0;
      double chi2;
      vector<PeakDef> peaks;
      computeForAcivity( activity, peaks, chi2, numDOF );
      
      ++num_iterations;
      
      if( (numDOF == 0) && (chi2 == 0.0) )
        throw runtime_error( "No DOF" );
      
      chi2cache.insert( std::pair<double,double>{activity,chi2} );
      
      return chi2;
    };//chi2ForAct
    
    
    //\TODO: if best activity is at minSearchActivity, it takes 50 iterations inside brent_find_minima
    //      to confirm; we could save this time by using just a little bit of intelligence...
    const int bits = 8; //Float has 24 bits of mantisa, so 12 would be its max precision; 8 should get us accurate to almost three significant figures
    boost::uintmax_t max_iter = 100;  //this variable gets changed each use, so you need to reset it afterwards
    const pair<double, double> r = boost::math::tools::brent_find_minima( chi2ForAct, minSearchActivity, maxSearchActivity, bits, max_iter );
  
    overallBestChi2 = r.second;
    overallBestActivity = r.first;
    
    cout << "Found min X2=" << overallBestChi2 << " with activity "
         << PhysicalUnits::printToBestActivityUnits(overallBestActivity)
         << " and it took " << std::dec << num_iterations << " iterations" << endl;
    
    //boost::math::tools::bracket_and_solve_root(...)
    auto chi2ForRangeLimit = [&chi2ForAct, overallBestChi2, yrange]( double const &activity ) -> double {
      return chi2ForAct(activity) - overallBestChi2 - yrange;
    };
    
    auto chi2ForCL = [&chi2ForAct, overallBestChi2, cl_chi2_delta]( double const &activity ) -> double {
      return chi2ForAct(activity) - overallBestChi2 - cl_chi2_delta;
    };
    
    //Tolerance is called with two values of activity; one with the chi2 bellow root, and one above
    auto tolerance = [chi2ForCL](double act1, double act2) -> bool{
      const double chi2_1 = chi2ForCL(act1);
      const double chi2_2 = chi2ForCL(act2);
      
      // \TODO: make sure tolerance is being used correctly - when pringting info out for every call I'm not sure it is being used right... (but answers seem reasonable, so...)
      cout << "Tolerance called with act1=" << PhysicalUnits::printToBestActivityUnits(act1)
           << ", act2=" << PhysicalUnits::printToBestActivityUnits(act2)
           << " ---> chi2_1=" << chi2_1 << ", chi2_2=" << chi2_2 << endl;
      
      return fabs(chi2_1 - chi2_2) < 0.025;
    };//tolerance(...)
    
    //cout << "chi2ForCL(minSearchActivity)=" << chi2ForCL(minSearchActivity) << endl;
    
    //Before trying to find lower-bounding activity, make sure the best value isnt the lowest
    //  possible value (i.e., zero in this case), and that if we go to the lowest possible value,
    //  that the chi2 will increase at least by cl_chi2_delta
    if( (fabs(minSearchActivity - overallBestActivity) > PhysicalUnits::nCi)
       && (chi2ForCL(minSearchActivity) > 0.0) )
    {
      pair<double,double> lower_val;

      max_iter = 100;  //see note about needing to set before every use
      lower_val = boost::math::tools::bisect( chi2ForCL, minSearchActivity, overallBestActivity, tolerance, max_iter );
      lowerLimit = 0.5*(lower_val.first + lower_val.second);
      foundLowerCl = true;
      cout << "lower_val CL activity=" << PhysicalUnits::printToBestActivityUnits(lower_val.first)
           << " wih chi2=" << chi2ForAct(lowerLimit) << ", num_iterations=" << std::dec << num_iterations
           << " and search range from " << PhysicalUnits::printToBestActivityUnits(minSearchActivity) << " to "
           << PhysicalUnits::printToBestActivityUnits(overallBestActivity)
           << endl;
      
      const double lowerLimitChi2 = chi2ForRangeLimit(lowerLimit);
      if( lowerLimitChi2 < 0.0 )
      {
        activityRangeMin = minSearchActivity;
        cout << "lower_val display activity being set to minSearchActivity (" << minSearchActivity << "): lowerLimitChi2=" << lowerLimitChi2 << endl;
      }else
      {
        max_iter = 100;
        lower_val = boost::math::tools::bisect( chi2ForRangeLimit, minSearchActivity, lowerLimit, tolerance, max_iter );
        activityRangeMin = 0.5*(lower_val.first + lower_val.second);
        foundLowerDisplay = true;
        cout << "lower_val display activity=" << PhysicalUnits::printToBestActivityUnits(activityRangeMin)
             << " wih chi2=" << chi2ForAct(activityRangeMin) << ", num_iterations=" << std::dec << num_iterations << endl;
      }
    }else
    {
      lowerLimit = 0.0;
      activityRangeMin = overallBestActivity;
      cout << "lower_val activity already at min" << endl;
    }//if( fabs(minSearchActivity - overallBestActivity) > PhysicalUnits::nCi )
    
    if( (fabs(maxSearchActivity - overallBestActivity) > PhysicalUnits::nCi)
         && (chi2ForCL(maxSearchActivity) > 0.0)  )
    {
      pair<double,double> upper_val;
      max_iter = 100;
      upper_val = boost::math::tools::bisect( chi2ForCL, overallBestActivity, maxSearchActivity, tolerance, max_iter );
      upperLimit = 0.5*(upper_val.first + upper_val.second);
      foundUpperCl = true;
      cout << "upper_val CL activity=" << PhysicalUnits::printToBestActivityUnits(upperLimit)
           << " wih chi2=" << chi2ForAct(upperLimit) << ", num_iterations=" << std::dec << num_iterations
           << " and search range from " << PhysicalUnits::printToBestActivityUnits(overallBestActivity) << " to "
           << PhysicalUnits::printToBestActivityUnits(maxSearchActivity)
           << endl;
      
      const double maxSearchChi2 = chi2ForRangeLimit(maxSearchActivity);
      if( maxSearchChi2 < 0.0 )
      {
        activityRangeMax = maxSearchActivity;
        cout << "upper_val display activity being set to maxSearchActivity (" << maxSearchActivity << "): maxSearchChi2=" << maxSearchChi2 << endl;
      }else
      {
        max_iter = 100;
        upper_val = boost::math::tools::bisect( chi2ForRangeLimit, upperLimit, maxSearchActivity, tolerance, max_iter );
        activityRangeMax = 0.5*(upper_val.first + upper_val.second);
        foundUpperDisplay = true;
        cout << "upper_val display activity=" << PhysicalUnits::printToBestActivityUnits(activityRangeMax)
             << " wih chi2=" << chi2ForAct(activityRangeMax) << ", num_iterations=" << std::dec << num_iterations << endl;
      }
    }else
    {
      upperLimit = overallBestActivity;
      activityRangeMax = overallBestActivity;
      cout << "upper_val activity already at max" << endl;
    }
    
    cout << "Found best chi2 and ranges with num_iterations=" << std::dec << num_iterations << endl;
    
    const double act_delta = fabs(activityRangeMax - activityRangeMin) / nchi2;
    for( size_t i = 0; i < nchi2; ++i )
    {
      const double act = activityRangeMin + act_delta*i;
      chi2s.push_back( {act,chi2ForAct(act)} );
    }
  }catch( std::exception &e )
  {
    m_bestChi2Act->setText( "" );
    m_upperLimit->setText( "" );
    m_results->hide();
    m_errorMsg->setText( "Error calculating Chi2: " + string(e.what()) );
    m_errorMsg->show();
    return;
  }//try / catch
  

  int numDOF = 0;
  double upperActivtyChi2 = -999.9; //upperActivtyChi2=overallBestChi2 + cl_chi2_delta
  std::vector<PeakDef> peaks;
  computeForAcivity( upperLimit, peaks, upperActivtyChi2, numDOF );
  
  
  char buffer[128];
  snprintf( buffer, sizeof(buffer), "Best &chi;<sup>2</sup> of %.1f at activity %s, %i DOF",
            overallBestChi2,
            PhysicalUnits::printToBestActivityUnits(overallBestActivity).c_str(),
            numDOF );
  m_bestChi2Act->setText( buffer );
  
  string upperLimitActStr = PhysicalUnits::printToBestActivityUnits(upperLimit,3,true);
  snprintf( buffer, sizeof(buffer), "95%% coverage at %s with &chi;<sup>2</sup> of %.1f",
            upperLimitActStr.c_str(), upperActivtyChi2 );
  
  if( foundUpperCl )
    m_upperLimit->setText( buffer );
  else
    m_upperLimit->setText( "Error: Didnt find upper 95%% limit" );
  
  if( foundLowerCl && (lowerLimit > 0.0) )
  {
    //display lower limit to user...
  }


  
  const bool useCurries = true;
  const double avrgAct = 0.5*(chi2s.front().first + chi2s.back().first);
  const auto &units = PhysicalUnits::bestActivityUnitHtml( avrgAct, useCurries );
  
  string datajson = "{\n\t\"data\": [";

  double minchi2 = std::numeric_limits<double>::infinity();
  double maxchi2 = -std::numeric_limits<double>::infinity();
  
  for( size_t i = 0; i < chi2s.size(); ++i )
  {
    minchi2 = std::min( minchi2, chi2s[i].second );
    maxchi2 = std::max( maxchi2, chi2s[i].second );
    datajson += string(i ? "," : "")
                + "{\"x\": " + std::to_string( (chi2s[i].first / units.second) )
                + ", \"y\": " + std::to_string( chi2s[i].second ) + "}";
  }
  datajson += "]";
  
  
  double chi2range = maxchi2 - minchi2;
  if( IsInf(minchi2) || IsInf(maxchi2) ) //JIC
  {
    minchi2 = 0;
    maxchi2 = 100;
    chi2range = 0;
  }
  
  //datajson += ",\n\t\"xunits\": \"" + units.first + "\"";
  
  datajson += ",\n\t\"xtitle\": \"";
  if( units.first == "&mu;Ci" )  /// \TODO: fix PhysicalUnits::bestActivityUnitHtml(..) to just us u character
    datajson += "Activity (" MU_CHARACTER "Ci)\"";
  else
    datajson += "Activity (" + units.first + ")\"";

  datajson += ",\n\t\"ytitle\": \"\xCF\x87\xC2\xB2\""; //chi^2
  
  double xstart = (chi2s.front().first/units.second);
  double xend = (chi2s.front().first/units.second);
  //Blah blah blah, coarse labels to be round numbers
  //const double initialrange = xend - xstart;
  
  //datajson += ",\n\t\"xrange\": [" + std::to_string(xstart) + ", " + std::to_string(xend) + "]";
  //datajson += ",\n\t\"yrange\": [" + std::to_string(minchi2 - 0.1*chi2range) + ", " + std::to_string(maxchi2 + 0.1*chi2range) + "]";
    
  // TODO: draw line at upperLimit.
  
  
  std::shared_ptr<const ColorTheme> theme = m_interspec ? m_interspec->getColorTheme() : nullptr;
  if( theme )
  {
    auto csstxt = []( const Wt::WColor &c, const char * defcolor ) -> string {
      return "\"" + (c.isDefault() ? string(defcolor) : c.cssText()) + "\"";
    };
     
    datajson += ",\n\t\"lineColor\": " + csstxt(theme->foregroundLine, "black");
    //datajson += ",\n\t\"axisColor\": " + csstxt(theme->spectrumAxisLines, "black");  //should be taken care of by CSS from D3SpectrumChart
    datajson += ",\n\t\"chartBackgroundColor\": " + csstxt(theme->spectrumChartBackground, "rgba(0,0,0,0)");
    //datajson += ",\n\t\"textColor\": " + csstxt(theme->spectrumChartText, "black"); //should be taken care of by CSS from D3SpectrumChart
  }//if( theme )
  
  datajson += "\n}";
  
  m_results->show();
  
  const string jsgraph = m_chi2Chart->jsRef() + ".chart";
  m_chi2Chart->doJavaScript( jsgraph + ".setData(" + datajson + ")" );
  
  m_displayActivity->setText( WString::fromUTF8(upperLimitActStr) );
  updateShownPeaks();
}//void doCalc()


void DetectionLimitTool::updateShownPeaks()
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
  std::vector<PeakDef> peaks;
  computeForAcivity( activity, peaks, chi2, numDOF );
  m_peakModel->setPeaks( peaks );
    
  cout << "Done in updateShownPeaks()" << endl;
}//void DetectionLimitTool::updateShownPeaks()


void DetectionLimitTool::computeForAcivity( const double activity,
                                                 std::vector<PeakDef> &peaks,
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
    
    cont->setRange( roi_start, roi_end );
    cont->calc_linear_continuum_eqn( spec, roi_start, roi_end, 2 );  //sets to linear continuum type
    
    const int continuumSelected = rw->m_continuum->currentIndex();
    switch( continuumSelected )
    {
      case 0:
        cont->setType( PeakContinuum::OffsetType::Linear );
        for( size_t order = 0; order < 2; ++order )
          cont->setPolynomialCoefFitFor( order, true );
        break;
        
      case 1:
        cont->setType( PeakContinuum::OffsetType::Quadratic );
        for( size_t order = 0; order < 3; ++order )
          cont->setPolynomialCoefFitFor( order, true );
        break;
    }//switch( assign continuum )
    
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
  
  for( auto &peak : fitPeaks )
  {
    peak.setFitFor( PeakDef::CoefficientType::Mean, false );
    peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
  }
  
  peaks.swap( fitPeaks );
}//void DetectionLimitTool::computeForAcivity(...)


void DetectionLimitTool::setRefLinesAndGetLineInfo()
{
  // \TODO: this function is essentually the same as #ReferencePhotopeakDisplay::updateDisplayChange
  //        and given its un-cleaness, it may be worth refactoring.
  
  if( !m_currentNuclide )
    return;
  
  auto spec = m_chart->data();
  if( !m_our_meas || !spec )
    return;
    
  std::shared_ptr<DetectorPeakResponse> drf = m_our_meas->detector();
  if( !drf || !drf->hasResolutionInfo() || !drf->isValid() )
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
  double age = PeakDef::defaultDecayTime( m_currentNuclide, &agestr );
  
  auto theme = m_interspec->getColorTheme();
  
  ReferenceLineInfo reflines;
  reflines.nuclide = m_currentNuclide;
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
  
  reflines.labelTxt = m_currentNuclide->symbol;
  
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
  mixture.addNuclideByActivity( m_currentNuclide, 1.0E-3 * SandiaDecay::curie );
  
  const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( age );
  const double parent_activity = m_currentNuclide ? mixture.activity( age, m_currentNuclide ) : 0.0;
  
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
    const bool is_decay_xray_gamma = (m_currentNuclide && (is_gamma || is_xray));
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
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      const SandiaDecay::Element *element = db->element( m_currentNuclide->atomicNumber );
      
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






void DetectionLimitTool::do_development()
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
