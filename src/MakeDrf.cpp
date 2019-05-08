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
#include <fstream>

#include <Wt/WText>
#include <Wt/WPanel>
#include <Wt/WCheckBox>
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
#include "InterSpec/SpecMeasManager.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;
using namespace Wt;

namespace
{
  bool source_info_from_lib_file( string srcname, const string &filename,
                                  double &activity, boost::posix_time::ptime &activityDate, string &comments )
  {
    ifstream file( filename.c_str(), ios::in | ios::binary );
    if( !file )
      return false;
   
    UtilityFunctions::trim( srcname );
    
    string line;
    while( UtilityFunctions::safe_get_line( file, line) )
    {
      vector<string> fields;
      UtilityFunctions::split( fields, line, " \t" );
      if( fields.size() < 3 )
        continue;
      if( !UtilityFunctions::iequals(fields[0], srcname) )
        continue;
      try
      {
        activity = stod(fields[1]) * PhysicalUnits::bq;
      }catch(...)
      {
        continue;
      }
      
      activityDate = UtilityFunctions::time_from_string( fields[2] .c_str());
      comments = "";
      for( size_t i = 3; i < fields.size(); ++i )
        comments += ((i==3) ? "" : " ") + fields[i];
      
      return true;
    }//while( UtilityFunctions::safe_get_line( file, line) )
    
    return false;
  }//source_info_from_lib_file(...)
  
  class DrfPeak : public WContainerWidget
  {
  public:
    const std::shared_ptr<const PeakDef> m_peak;
    const double m_livetime;
    Wt::WCheckBox *m_useCb;
    
  
    DrfPeak( std::shared_ptr<const PeakDef> peak, double livetime, WContainerWidget *parent = nullptr )
    : WContainerWidget( parent ),
      m_peak( peak ),
      m_livetime( livetime ),
      m_useCb( nullptr )
    {
      addStyleClass( "DrfPeak" );
      m_useCb = new WCheckBox( "Use", this );
      
      char buffer[512];
      if( peak && peak->parentNuclide() && peak->nuclearTransition() )
      {
        m_useCb->setChecked();
        snprintf( buffer, sizeof(buffer), "%s: %.2f keV peak with %.1f cps for %.2f keV gamma.",
                  peak->parentNuclide()->symbol.c_str(),
                  peak->mean(),
                  (peak->amplitude() / livetime),
                  peak->gammaParticleEnergy() );
      }else
      {
        m_useCb->setUnChecked();
        m_useCb->disable();
        disable();
        snprintf( buffer, sizeof(buffer), "%.2f keV peak - no source associated", peak->mean() );
      }
      
      WText *txt = new WText( WString::fromUTF8(buffer), this );
    }//DrfPeak constructor;
  };//class DrfPeak
  
  
  class DrfSpecFileSample : public WPanel
  {
    std::shared_ptr<const SpecMeas> m_meas;
    set<int> m_samples;
    MaterialDB *m_materialDB;
    Wt::WSuggestionPopup *m_materialSuggest;
    WContainerWidget *m_peaks;
    WContainerWidget *m_sources;
    Wt::Signal<> m_srcInfoUpdated;
    
  public:
    DrfSpecFileSample( std::shared_ptr<const SpecMeas> meas, set<int> samples,
                      MaterialDB *materialDB,
                      Wt::WSuggestionPopup *materialSuggest,
                      WContainerWidget *parent = nullptr )
    : WPanel( parent ),
      m_meas( meas ),
      m_samples( samples ),
      m_materialDB( materialDB ),
      m_materialSuggest( materialSuggest ),
      m_peaks( nullptr ),
      m_sources( nullptr )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      
      addStyleClass( "DrfSpecFileSample" );
      
      WContainerWidget *body = new WContainerWidget();
      setCentralWidget( body );
      
      m_peaks = new WContainerWidget( body );
      m_sources = new WContainerWidget( body );
      
      std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaksptr = meas->peaks( samples );
      
      //Make a widget representing sample number combination with peaks in file
      auto peaks = make_shared< deque<shared_ptr<const PeakDef>> >( *peaksptr );
      assert( peaks && !peaks->empty() );
      
      string title;
      double livetime = 0.0;
      //ToDo: The peaks may not belong to all the detectors (the user may have
      //      chosen to not show some detectors)
      const vector<string> &detnames = meas->detector_names();
      for( const int sample : samples )
      {
        for( const string &detname : detnames )
        {
          shared_ptr<const Measurement> m = meas->measurement( sample, detname );
          if( m && !m->title().empty() )
            title = m->title();
          if( m && m->num_gamma_channels() > 5 )
            livetime += m->live_time();
        }//
      }//for( const int sample : samples )
      
      string filename = UtilityFunctions::filename( meas->filename() );
      if( meas->sample_numbers().size()==1 && !title.empty() )
      {
        title = filename + ": " + title;
      }else if( samples.size() == 1 && !title.empty() )
      {
        title = filename + ": Sample " + std::to_string(*samples.begin()) + ", &quot;" + title + "&quot;";
      }else if( samples.size() == 1 && title.empty() )
      {
        title = filename + ": Sample " + std::to_string(*samples.begin());
      }else if( title.empty() )
      {
        title = filename + ": Samples {";
        for( auto iter = begin(samples); iter != end(samples); ++iter )
          title += (iter==begin(samples) ? "" : ",") + std::to_string(*iter);
        title += "}";
      }else  //we have a title and multiple sample numbers
      {
        title = filename + ": &quot;" + title + "&quot; Samples {";
        for( auto iter = begin(samples); iter != end(samples); ++iter )
          title += (iter==begin(samples) ? "" : ",") + std::to_string(*iter);
        title += "}";
      }//if / else: choose how to structure the title
      
      setTitle( WString::fromUTF8(title) );
      
      for( auto peak : *peaks )
      {
        DrfPeak *p = new DrfPeak( peak, livetime, m_peaks );
        p->m_useCb->changed().connect( this, &DrfSpecFileSample::refreshSourcesVisible );
      }
      
      refreshSourcesVisible();
      
      //See if the file conforms to GADRAS conventions
      //if( samples.size() == 1 )
      {
        double distance = -1.0, activity = -1.0;
        boost::posix_time::ptime activityDate;
        const SandiaDecay::Nuclide *nuc = nullptr;
        
        string srcinfo;
        const vector<string> &dets = meas->detector_names();
        for( const string &det : dets )
        {
          for( const int sample : samples )
          {
            auto m = meas->measurement(sample, det);
            if( !m )
              continue;
            
            string spectitle = m->title();
            size_t pos = spectitle.find( "@" );
            if( pos != string::npos )
            {
              //Look for distance and source info in title. Examples:
              //  "U-232 @ 100 cm, H=100 cm"
              //  "Background, H=100 cm"
              string dist = spectitle.substr( pos+1 );
              try
              {
                UtilityFunctions::trim( dist );
                pos = dist.find( ',' ); //ToDo: make finding the end of the distnace more robust - like with a regex
                if( pos != string::npos )
                  dist = dist.substr(0, pos);
                distance = PhysicalUnits::stringToDistance( dist );
              }catch(...){}
            }//if( pos != string::npos )
            
            //make work for strings like "133BA_blahblahblah"
            const size_t first_space = spectitle.find_first_of( " \t_" );
            //ToDo: make this more robost for finding end of the nuclide, like requiring both numbers and letters.
            
            nuc = db->nuclide( spectitle.substr(0,first_space) );
            
            
            for( string remark : m->remarks() )
            {
              if( !UtilityFunctions::istarts_with(remark, "Source:") )
                continue;
              
              remark = remark.substr(7);
              UtilityFunctions::trim( remark );
              
              if( !nuc )
                nuc = db->nuclide( remark );  //Probably shouldnt give to many false positives?
              
              const size_t commapos = remark.find(',');
              const size_t underscorpos = remark.find('_'); //For a source database name like "133BA_<serial number>"
              
              if( (commapos < underscorpos)
                 && (nuc == db->nuclide(remark.substr(0,commapos))) )
              {
                //We have something like "133Ba,10uCi"
                try
                {
                  activity = PhysicalUnits::stringToActivity( remark.substr(commapos+1) );
                }catch(...)
                {}
              }//
              
              if( (activity < 0.0) && (underscorpos != string::npos)
                 && (nuc == db->nuclide(remark.substr(0,underscorpos)) ))
              {
                //- Check for "*Source.lib" file in application directory,
                //  and if found, scan through for the source identified in previous line, and put in date
                vector<string> base_paths{ InterSpec::staticDataDirectory(), UtilityFunctions::get_working_path() };
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__)) ) )
                try{ base_paths.push_back( InterSpec::writableDataDirectory() ); } catch(...){}
#endif
                vector<string> source_lib_files;
                for( const string &path : base_paths )
                {
                  const vector<string> files = UtilityFunctions::recursive_ls(path, "Source.lib" );
                  source_lib_files.insert( end(source_lib_files), begin(files), end(files) );
                }
                
                string comments;
                for( const string &lib : source_lib_files )
                {
                  if( source_info_from_lib_file(remark, lib, activity, activityDate, comments) )
                    break;
                }//for( const string &lib : source_lib_files )
              }//if( activty < 0.0 )
              
            }//for( string remark : m->remarks() )
          }//for( const int sample : samples )
          
        }//for( const string &det : dets )
        
        if( nuc )
        {
          for( auto w : m_sources->children() )
          {
            auto src = dynamic_cast<MakeDrfSrcDef *>(w);
            if( src && src->nuclide()==nuc )
            {
              if( distance > 0.0 )
                src->setDistance( distance );
              if( activity > 0.0 && activityDate.is_special() )
                src->setActivity( activity );
              else if( activity > 0.0 )
                src->setAssayInfo( activity, activityDate );
            }//if( we found the source )
          }//for( auto w : m_sources->children() )
        }//if( nuc )
        
        //Look through comments of each spectra for a comment similar to "Source: 232U_NIST0623220{92,0.07}", and if curly brackets add shielding to source
        //  - If source is like "Source: 57CO,16.79uC" then intpret appropriately.
        //
        //If measurement has a valid date, put into source widget.
        //Implement background subtraction (for instrinsic Cs137, or U232)
        
        
        if( samples.size() == 1 )
        {
          //check if background, and if so, deslect all peaks.
        }
        
        
      }//if( samples.size() == 1 )
      
      
      
      
    }//DrfSpecFileSample constructor
    
    
    void refreshSourcesVisible()
    {
      //ToDo: order nuclides by some predictable reasonable way, instead of sorted by pointer
      set<const SandiaDecay::Nuclide *> selectedNucs;
      map<const SandiaDecay::Nuclide *,MakeDrfSrcDef *> nucToWidget;
      
      for( auto w : m_peaks->children() )
      {
        auto p = dynamic_cast<DrfPeak *>( w );
        if( p && p->m_useCb->isChecked() )
          selectedNucs.insert( p->m_peak->parentNuclide() );
      }
      
      for( auto w : m_sources->children() )
      {
        auto p = dynamic_cast<MakeDrfSrcDef *>( w );
        if( p && p->nuclide() )
          nucToWidget[p->nuclide()] = p;
      }
      
      for( auto n : selectedNucs )
      {
        if( !nucToWidget.count(n) )
        {
          boost::posix_time::ptime measDate;
          const vector<string> &detNames = m_meas->detector_names();
          for( size_t i = 0; measDate.is_special() && i < detNames.size(); ++i )
          {
            for( const int sample : m_samples )
            {
              auto m = m_meas->measurement( sample, detNames[i] );
              if( m && !m->start_time().is_special() )
              {
                measDate = m->start_time();
                break;
              }
            }//for( const int sample : m_samples )
          }//for( size_t i = 0; i < detNames.size(); ++i )
          
          MakeDrfSrcDef *src = new MakeDrfSrcDef( n, measDate, m_materialDB, m_materialSuggest, m_sources );
          src->updated().connect( std::bind([this](){ m_srcInfoUpdated.emit(); }) );
          
          nucToWidget[n] = src;
        }//if( we dont have a source widget for a nuclide )
      }//for( auto n : selectedNucs )
      
      for( auto &p : nucToWidget )
        p.second->setHidden( !selectedNucs.count(p.first) );
      
      m_srcInfoUpdated.emit();
    }//void refreshSourcesVisible()
    
    
    Wt::Signal<> &srcInfoUpdated(){ return m_srcInfoUpdated; };
    
  };//class DrfSpecFileSample
  
  
  class DrfSpecFile : public WPanel
  {
    std::shared_ptr<const SpecMeas> m_meas;
    Wt::Signal<> m_updated;
    
  public:
    DrfSpecFile( std::shared_ptr<const SpecMeas> meas,
                MaterialDB *materialDB,
                Wt::WSuggestionPopup *materialSuggest,
                WContainerWidget *parent = nullptr )
      : WPanel( parent ),
        m_meas( meas )
    {
      addStyleClass( "DrfSpecFile" );
      
      const set<set<int>> sampsWithPeaks = meas->sampleNumsWithPeaks();

      WContainerWidget *sampleWidgets = new WContainerWidget();
      setCentralWidget( sampleWidgets );
      
      for( const set<int> &peakSamps : sampsWithPeaks )
      {
        std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaksptr = meas->peaks( peakSamps );
        if( !peaksptr || peaksptr->empty() )
          continue;
        
        auto sample = new DrfSpecFileSample( meas, peakSamps, materialDB, materialSuggest, sampleWidgets );
        sample->srcInfoUpdated().connect( std::bind( [this](){ m_updated.emit(); }) );
      }//for( const set<int> &peakSamps : sampsWithPeaks )
    }//DrfSpecFile(...)
    
    Wt::Signal<> &updated() { return m_updated; }
  };//class DrfSpecFile
  
  
  
}//namespace to implement DrfSpecFile and DrfSpecFileSample


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
  
  WContainerWidget *files = new WContainerWidget();
  files->setMaximumSize( WLength::Auto, std::max((wh - chartHeight - 100), 250) );
  files->setOverflow( Overflow::OverflowAuto, Wt::Vertical );
  
  layout->addWidget( files, 1, 0, 1, 2 );
  
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
  
  SpecMeasManager *manager = viewer->fileManager();
  SpectraFileModel *measmodel = manager->model();
  
  const int nLoadedFiles = measmodel->rowCount();
  for( int filenum = 0; filenum < nLoadedFiles; ++filenum )
  {
    std::shared_ptr<SpectraFileHeader> header = measmodel->fileHeader( filenum );
    if( !header )
      continue;
    
    std::shared_ptr<SpecMeas> meas = header->parseFile();
    if( !meas )
      continue; //JIC
    
    const set<set<int>> sampsWithPeaks = meas->sampleNumsWithPeaks();
    if( sampsWithPeaks.empty() )
      continue;
    
    //Make a widget representing each file
    DrfSpecFile *fileWidget = new DrfSpecFile( meas, materialDB, materialSuggest, files );
    fileWidget->updated().connect( this, &MakeDrf::handleSourcesUpdates );
  }//for( loop over opened files )
  
  
  /*
  shared_ptr<SpecMeas> meas = viewer->measurment(kForeground);
  if( meas )
  {
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
  */
  
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


void MakeDrf::handleSourcesUpdates()
{
  
}//void MakeDrf::handleSourcesUpdates()
