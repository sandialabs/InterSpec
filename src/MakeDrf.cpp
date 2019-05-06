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

#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WSuggestionPopup>

#include "InterSpec/MakeDrf.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MakeDrfChart.h"
#include "InterSpec/InterSpecApp.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;
using namespace Wt;


MakeDrfWindow::MakeDrfWindow( InterSpec *viewer, MaterialDB *materialDB, Wt::WSuggestionPopup *materialSuggest )
  : AuxWindow( "Create DRF",
  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletModal)
  | AuxWindowProperties::SetCloseable
  | AuxWindowProperties::DisableCollapse) ),
  m_makeDrf( nullptr )
{
  m_makeDrf = new MakeDrf( viewer, materialDB, materialSuggest );
  
  stretcher()->addWidget( m_makeDrf, 0, 0 );
  stretcher()->setContentsMargins( 0, 0, 0, 0 );
  
  AuxWindow::addHelpInFooter( footer(), "make-drf" );
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  
  //const int screenW = viewer->renderedWidth();
  //const int screenH = viewer->renderedHeight();
  //const int width = ((screenW < 600) ? screenW : (3*screenW)/4);
  //const int height = ((screenH < 420) ? screenH : (5*screenH)/6 );
  //resizeWindow( width, height );
  
  resizeToFitOnScreen();
  centerWindow();
  rejectWhenEscapePressed( false );
}//MakeDrfWindow(...) constrctor


MakeDrfWindow::~MakeDrfWindow()
{
}


MakeDrf::MakeDrf( InterSpec *viewer, MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggest,
                  Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_materialDB( materialDB ),
  m_materialSuggest( materialSuggest ),
  m_chart( nullptr )
{
  assert( m_interspec );
  assert( m_materialDB );
  assert( m_materialSuggest );
  
  wApp->useStyleSheet( "InterSpec_resources/MakeDrf.css" );
  
  addStyleClass( "MakeDrf" );
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( layout );
  
  WContainerWidget *optionsDiv = new WContainerWidget();
  optionsDiv->addStyleClass( "MakeDrfOptions" );
  layout->addWidget( optionsDiv, 0, 0 );
  
  m_chart = new MakeDrfChart();
  layout->addWidget( m_chart, 0, 1 );
  
  //m_chart->setMinimumSize( 250, 250 );
  const int wh = viewer->renderedHeight();
  const int ww = viewer->renderedWidth();
  const int chartHeight = std::min( wh/3, 350 );
  const int chartWidth = std::min( static_cast<int>(0.75*ww - 150), 800 );
  m_chart->resize( chartWidth, chartHeight );
  
  layout->addWidget( new WContainerWidget(), 1, 0, 1, 2 );
  
  //layout->setRowStretch( 0, 1 );
  //layout->setRowStretch( 1, 3 );
  //layout->setColumnStretch( 1, 1 );
  

  //Make a chart to show relative efficiency
  //  make energy extents as large as any of the input files peaks are being used from
  //  gray out above and below first/last points.
  
  {
    vector<float> effCoefs{ -1.96451f, -1.1967f, -0.1079f, 0.0955f, 0.0381f, 0.0065f };
    m_chart->setEfficiencyCoefficients( effCoefs, MakeDrfChart::EqnEnergyUnits::MeV );
    
    const float diameter = 1.780*2*PhysicalUnits::cm;
    const float distance = 25*PhysicalUnits::cm;
    vector<MakeDrfChart::DataPoint> datas;
    
    MakeDrfChart::DataPoint point;
    point.distance = distance;
    point.livetime = 297.0f * PhysicalUnits::second;
    
    point.energy = 367.94;
    point.peak_area = 26037;
    point.source_count_rate = 171982;
    point.source_information = "Tl200";//0.390728 abs eff, 0.0126 solid angle, Tl200
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    point.energy = 439.51;
    point.peak_area = 27463;
    point.source_count_rate = 217047;
    point.source_information = "Tl202"; //0.335635 abs eff,
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    point.energy = 579.3;
    point.peak_area = 4332.19;
    point.source_count_rate = 67956.9;
    point.source_information = "Tl200";  //0.257736, intrinsic eff
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    point.energy = 661.657;
    point.peak_area = 4200;
    point.source_count_rate = 49708;
    point.source_information = "Cs137"; //int. eff. 0.224391
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    point.energy = 828.27;
    point.peak_area = 3506;
    point.source_count_rate = 79106.9;
    point.source_information = "Tl200"; //0.174915 int eff
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    point.energy = 1205.75;
    point.peak_area = 15589;
    point.source_count_rate = 342305;
    point.source_information = "Tl200"; //Int. Eff0.111747
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    point.energy = 1363.2;
    point.peak_area = 1628;
    point.source_count_rate = 66023;
    point.source_information = "Tl200"; //Int Eff 0.0960932
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    point.energy = 1514.9;
    point.peak_area = 10771;
    point.source_count_rate = 49665.3;
    point.source_information = "Tl200"; //Int Eff 0.844061
    point.peak_area_uncertainty = sqrt( point.peak_area );
    point.source_count_rate_uncertainty = 0.25f * point.source_count_rate;
    datas.push_back( point );
    
    m_chart->setDataPoints( datas, diameter, 0, 3000.0 );
  }
  
  shared_ptr<SpecMeas> meas = viewer->measurment(kForeground);
  if( meas )
  {
    const set<set<int>> sampsWithPeaks = meas->sampleNumsWithPeaks();
    for( const set<int> &peakSamps : sampsWithPeaks )
    {
      std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > > peaks = meas->peaks( peakSamps );
      
    }//for( const set<int> &peakSamps : sampsWithPeaks )
    
    
    shared_ptr<DetectorPeakResponse> det = meas->detector();
    if( det && (det->efficiencyFcnType()==DetectorPeakResponse::EfficiencyFnctForm::kExpOfLogPowerSeries) )
    {
      if( det->resolutionFcnType() == DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn )
      {
        const vector<float> &resCoef = det->resolutionFcnCoefficients();
        m_chart->setFwhmCoefficients( resCoef, MakeDrfChart::EqnEnergyUnits::keV );
      }
    }//if( meas )
  }//if( meas )
  
  
  //Make an input for detector dimensions
  
  //Make an input for order of efficiency fit
  
  //Make an input for order of FWHM fit
  
  //Make Text display of DRF equations for copy/paste
  
  //Make display for each file that is loaded with any peaks fit
  //  Make display for any combination of sample numbers that have a peak fit for each file
  //    Make option to use/not-use each peak
  //      It would be nice to show a 40px by 40px peak preview for each peak
  //    Make source MakeDrfSrcDef for each nuclide used
  //Make way to save DRF...
  
  /*
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( "Cs137" );
  
  boost::posix_time::ptime measDate = UtilityFunctions::time_from_string( "2019-01-01T00:00:00" );
  
  
  MakeDrfSrcDef *srcinput = new MakeDrfSrcDef( nuc, measDate, materialDB, materialSuggest, this );
  
  new WText( "<br />", this );
  
  nuc = db->nuclide( "Am241" );
  srcinput = new MakeDrfSrcDef( nuc, measDate, materialDB, materialSuggest, this );
  */
}//MakeDrf( constructor )


MakeDrf::~MakeDrf()
{
  
}//~MakeDrf()
