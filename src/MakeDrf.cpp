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
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WSvgImage>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WPopupWidget>
#include <Wt/WDoubleSpinBox>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>

#include "InterSpec/MakeDrf.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MakeDrfChart.h"
#include "InterSpec/InterSpecApp.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SpecMeasManager.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ShieldingSourceDisplay.h"


using namespace std;
using namespace Wt;

namespace
{
  bool source_info_from_lib_file( string srcname, const string &filename,
                                  double &activity, boost::posix_time::ptime &activityDate, string &comments )
  {
    cout << "Will try to find '" << srcname << "' from " << filename << endl;
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
    WText *m_descTxt;
    WText *m_backSubTxt;
    WDoubleSpinBox *m_userBr;
    WPushButton *m_previewBtn;
    WPopupWidget *m_previewPopup;
    std::shared_ptr<const Measurement> m_summed_meas;
    
    DrfPeak( std::shared_ptr<const PeakDef> peak, double livetime,
             std::shared_ptr<const Measurement> summed_meas,
             WContainerWidget *parent = nullptr )
    : WContainerWidget( parent ),
      m_peak( peak ),
      m_livetime( livetime ),
      m_useCb( nullptr ),
      m_descTxt( nullptr ),
      m_backSubTxt( nullptr ),
      m_userBr( nullptr ),
      m_previewBtn( nullptr ),
      m_previewPopup( nullptr ),
      m_summed_meas( summed_meas )
    {
      addStyleClass( "DrfPeak" );
      m_useCb = new WCheckBox( "Use", this );
      
      char buffer[512];
      if( peak && peak->parentNuclide() && peak->nuclearTransition() )
      {
        const bool use = (peak->PeakDef::sourceGammaType()==PeakDef::SourceGammaType::NormalGamma);
        m_useCb->setChecked( use );
        snprintf( buffer, sizeof(buffer), "%s: %.2f keV peak with %.1f cps for %.2f keV gamma.",
                  peak->parentNuclide()->symbol.c_str(),
                  peak->mean(),
                  (peak->amplitude() / livetime),
                  peak->gammaParticleEnergy() );
      }else
      {
        snprintf( buffer, sizeof(buffer), "%.2f keV peak - no nuc. associated", peak->mean() );
        m_userBr = new WDoubleSpinBox();
        //m_userBr->addStyleClass( "DrfPeakNoNucBr" );
        m_userBr->setRange( 0.0, 1.0 );
        m_userBr->setSingleStep( 0.01 );
        m_userBr->setValue( 1.0 );
        m_userBr->setMargin( 5, Wt::Left );
        m_userBr->setTextSize( 8 );
        
        m_useCb->setUnChecked();
        m_userBr->hide();
        m_useCb->changed().connect( std::bind( [this](){
          m_userBr->setHidden( !m_useCb->isChecked() );
        } ) );
      }
      
      m_descTxt = new WText( WString::fromUTF8(buffer), this );
      m_descTxt->addStyleClass( "DrfPeakInfoTxt" );
      m_backSubTxt = new WText( "", this );
      m_backSubTxt->addStyleClass( "DrfPeakBackSubTxt" );
      m_backSubTxt->hide();
      
      if( m_userBr )
      {
        WLabel *l = new WLabel( "BR=", this );
        l->addStyleClass( "DrfPeakNoNucBrLbl" );
        l->setBuddy( m_userBr );
        addWidget( m_userBr );
      }
      
      m_previewBtn = new WPushButton( "", this );
      m_previewBtn->setStyleClass( "DrfPeakShowChart" );
      m_previewBtn->setIcon( "InterSpec_resources/images/peak_small.png" );
      //m_previewBtn = new WImage( "InterSpec_resources/images/peak_small.png", this );
      
      m_previewBtn->clicked().connect( this, &DrfPeak::togglePeakPreview );
    }//DrfPeak constructor;
    
    
    void setBackgroundBeingSubtractedInfo( const bool beingsubtracted, const float background_cps )
    {
      const bool wasSubtractedFrom = !m_backSubTxt->isHidden();
      if( wasSubtractedFrom == beingsubtracted )
        return;
      
      if( beingsubtracted )
      {
        char buffer[64];
        snprintf( buffer, sizeof(buffer)-1, "will sub. %.3f cps for bckgrnd", background_cps );
        m_backSubTxt->setText( WString::fromUTF8(buffer) );
        m_backSubTxt->show();
      }else
      {
        m_backSubTxt->setText( "" );
        m_backSubTxt->hide();
      }
    }
    
    void togglePeakPreview()
    {
      //WAnimation animation(WAnimation::AnimationEffect::Pop, WAnimation::Linear, 250 );
      if( m_previewPopup )
      {
        if( m_previewPopup->isHidden() )
        {
          m_previewBtn->addStyleClass( "active" );
          m_previewPopup->setHidden( false );
        }else
        {
          m_previewBtn->removeStyleClass( "active" );
          m_previewPopup->setHidden( true );
        }
        
        return;
      }//if( m_previewPopup )
      
      try
      {
        if( !m_peak )
          throw runtime_error( "No peak for preview" );
        
        auto peakdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>( 1, m_peak );
        vector<shared_ptr<const ReferenceLineInfo>> reflines;  //ToDo: actually fill these out
        const double lower_energy = m_peak->mean() - 0.75*m_peak->roiWidth();
        const double upper_energy = m_peak->mean() + 0.75*m_peak->roiWidth();
        
        const int chart_width_px = 225;
        const int chart_height_px = 125;
        std::shared_ptr<const ColorTheme> theme; //ToDo: actually get this from the InterSpec class...
        
        shared_ptr<WSvgImage> svg
        = PeakSearchGuiUtils::renderChartToSvg( m_summed_meas, peakdeque,
                                               reflines, lower_energy, upper_energy,
                                               chart_width_px, chart_height_px, theme, true );
        if( !svg )
          throw runtime_error( "No preview svg chart" );
        
        stringstream strm;
        svg->write( strm );
        
        WText *chart = new WText( strm.str(), Wt::XHTMLUnsafeText /*, this*/ );
        chart->addStyleClass( "DrfPeakChart" );
        m_previewPopup = new WPopupWidget( chart, this );
      }catch( std::exception & )
      {
        WText *msg = new WText( "Unable to generate preview" );
        m_previewPopup = new WPopupWidget( msg, this );
      }//try / catch make a chart
      
      m_previewPopup->setAnchorWidget( m_previewBtn, Wt::Orientation::Vertical );
      m_previewPopup->setJavaScriptMember("wtNoReparent", "true");
      m_previewPopup->show();
      m_previewBtn->addStyleClass( "active" );
    }//void togglePeakPreview()
    
  };//class DrfPeak
  
  
  class DrfSpecFileSample : public WPanel
  {
    std::shared_ptr<const SpecMeas> m_meas;
    set<int> m_samples;
    MaterialDB *m_materialDB;
    WSuggestionPopup *m_materialSuggest;
    WContainerWidget *m_peaks;
    WContainerWidget *m_sources;
    WCheckBox *m_background;
    WText *m_backgroundTxt;
    WCheckBox *m_allNoneSome;
    Wt::Signal<> m_srcInfoUpdated;
    
    void updateTitle()
    {
      if( !m_meas )
        return;
      
      string title;
      
      const vector<string> &detnames = m_meas->detector_names();
      for( const int sample : m_samples )
      {
        for( const string &detname : detnames )
        {
          shared_ptr<const Measurement> m = m_meas->measurement( sample, detname );
          if( m && !m->title().empty() )
            title = m->title();
        }//
      }//for( const int sample : samples )
      
      
      string filename = m_meas->filename();
      filename = UtilityFunctions::filename( filename );
      if( m_meas->sample_numbers().size()==1 && !title.empty() )
      {
        title = filename + (filename.empty() ? "" : ": ") + title;
      }else if( m_samples.size() == 1 && !title.empty() )
      {
        title = filename + (filename.empty() ? "" : ": ") + "Sample " + std::to_string(*m_samples.begin()) + ", &quot;" + title + "&quot;";
      }else if( m_samples.size() == 1 && title.empty() )
      {
        title = filename + (filename.empty() ? "" : ": ") + " Sample " + std::to_string(*m_samples.begin());
      }else if( title.empty() )
      {
        title = filename + (filename.empty() ? "" : ": ") + " Samples {";
        for( auto iter = begin(m_samples); iter != end(m_samples); ++iter )
          title += (iter==begin(m_samples) ? "" : ",") + std::to_string(*iter);
        title += "}";
      }else  //we have a title and multiple sample numbers
      {
        title = filename + ": &quot;" + title + "&quot; Samples {";
        for( auto iter = begin(m_samples); iter != end(m_samples); ++iter )
          title += (iter==begin(m_samples) ? "" : ",") + std::to_string(*iter);
        title += "}";
      }//if / else: choose how to structure the title
      
      //Count the number of peaks being used
      int npeaks = 0, npeaksused = 0;
      
      for( auto w : m_peaks->children() )
      {
        auto p = dynamic_cast<DrfPeak *>( w );
        npeaks += (p ? 1 : 0);
        npeaksused += ((p && p->m_useCb->isChecked()) ? 1 : 0);
      }
      
      title += " (using " + std::to_string(npeaksused) + " of " + std::to_string(npeaks) + " peaks)";
      
      if( m_allNoneSome )
      {
        if( npeaksused == npeaks )
          m_allNoneSome->setCheckState( CheckState::Checked );
        else if( npeaksused )
          m_allNoneSome->setCheckState( CheckState::PartiallyChecked );
        else
          m_allNoneSome->setCheckState( CheckState::Unchecked );
      }//if( m_allNoneSome )
      
      setTitle( WString::fromUTF8(title) );
    }//void updateTitle()
    
    
    void isBackgroundToggled()
    {
      m_backgroundTxt->setHidden( !m_background->isChecked() );
      m_sources->setHidden( m_background->isChecked() );
      m_srcInfoUpdated.emit();
    }//void isBackgroundToggled()
    
    
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
      m_sources( nullptr ),
      m_background( nullptr ),
      m_backgroundTxt( nullptr ),
      m_allNoneSome( nullptr )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      
      addStyleClass( "DrfSpecFileSample" );
      
      WContainerWidget *body = new WContainerWidget();
      setCentralWidget( body );
      
      m_peaks = new WContainerWidget( body );
      
      WContainerWidget *bckgrndDiv = new WContainerWidget( body );
      bckgrndDiv->addStyleClass( "DrfSpecFileSampleBackground" );
      m_background = new WCheckBox( "Is Background?", bckgrndDiv );
      m_background->changed().connect( this, &DrfSpecFileSample::isBackgroundToggled );
      m_background->setChecked( false );
      m_backgroundTxt = new WText( "These peaks will be subtracted from other samples", bckgrndDiv );
      m_backgroundTxt->addStyleClass( "DrfSpecFileSampleBackgroundTxt" );
      m_backgroundTxt->hide();
      
      
      m_sources = new WContainerWidget( body );
      
      setCollapsible( true );
      
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
      
      std::shared_ptr<const Measurement> summed_meas;
      if( samples.size()==1 && detnames.size()==1 )
        summed_meas = meas->measurement( *begin(samples), detnames[0] );
      if( !summed_meas )
        summed_meas = meas->sum_measurements( samples, meas->detector_numbers() );
      
      //ToDo: inprinciple summed_meas->live_time() should equal livetime; should check
      
      int npeaks = 0;
      for( auto peak : *peaks )
      {
        DrfPeak *p = new DrfPeak( peak, livetime, summed_meas, m_peaks );
        ++npeaks;
        p->m_useCb->changed().connect( this, &DrfSpecFileSample::refreshSourcesVisible );
        if( p->m_userBr )
        {
          p->m_userBr->changed().connect( this, &DrfSpecFileSample::refreshSourcesVisible );
          p->m_userBr->enterPressed().connect( this, &DrfSpecFileSample::refreshSourcesVisible );
        }
      }//for( auto peak : *peaks )
      
      
      if( npeaks > 1 )
      {
        m_allNoneSome = new WCheckBox( "All Peaks", titleBarWidget() );
        m_allNoneSome->clicked().preventPropagation();
        m_allNoneSome->addStyleClass( "DrfSpecFileAllNoneSomeCb" );
        m_allNoneSome->changed().connect( this, &DrfSpecFileSample::handleUserToggleAllNoneSum );
      }
      
      
      refreshSourcesVisible();
      
      //See if the file conforms to GADRAS conventions
      //if( samples.size() == 1 )
      {
        string shielding;
        bool is_background = false;
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
            
            if( samples.size() == 1 )
            {
              is_background |= (m->source_type()==Measurement::SourceType::Background);
              is_background |= UtilityFunctions::icontains( spectitle, "back" );
              is_background |= UtilityFunctions::icontains( spectitle, "bgr" );
            }//if( samples.size() == 1 )
            
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
              
              const size_t openCurlyPos = remark.find( '{' );
              if( shielding.empty() && (openCurlyPos != string::npos) )
              {
                const size_t closeCurlyPos = remark.find( '}', openCurlyPos );
                if( (closeCurlyPos != string::npos) && ((openCurlyPos+4)<closeCurlyPos) )
                {
                  shielding = remark.substr(openCurlyPos+1, closeCurlyPos-openCurlyPos-1);
                  remark = remark.substr(0,openCurlyPos);
                  UtilityFunctions::trim(remark);
                }
              }//if( source had shielding defined )
              
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
              
              if( !shielding.empty() )
              {
                try
                {
                  vector<float> an_ad;
                  UtilityFunctions::split_to_floats( shielding, an_ad );
                  if( an_ad.size() >= 2 && an_ad[0] >= 1.0f && an_ad[0] <= 100.0f
                     && an_ad[1] >= 0.0f && an_ad[1] <= 500.0f )
                  {
                    src->setShielding( an_ad[0], an_ad[1]*PhysicalUnits::gram/PhysicalUnits::cm2 );
                  }
                }catch( std::exception &e )
                {
                  cerr << "DrfSpecFileSample: Caught exception setting shielding from '" << shielding << "': " << e.what() << endl;
                }
              }//if( !shielding.empty() )
            }//if( we found the source )
          }//for( auto w : m_sources->children() )
        }//if( nuc )
        
        if( is_background )
        {
          m_background->setChecked( true );
          m_backgroundTxt->show();
          m_sources->hide();
        }
        
        //Look through comments of each spectra for a comment similar to "Source: 232U_NIST0623220{92,0.07}", and if curly brackets add shielding to source
        //  - If source is like "Source: 57CO,16.79uC" then intpret appropriately.
        //
        //If measurement has a valid date, put into source widget.
        //Implement background subtraction (for instrinsic Cs137, or U232)
        
        
        if( samples.size() == 1 )
        {
          //check if background, and if so, deselect all peaks.
        }
        
      }//if( samples.size() == 1 )
    }//DrfSpecFileSample constructor
    
    void handleUserToggleAllNoneSum()
    {
      assert( m_allNoneSome );
      
      bool useAll = true;
      switch( m_allNoneSome->checkState() )
      {
        case Wt::Checked:   useAll = true;  break;
        case Wt::Unchecked: useAll = false; break;
        case Wt::PartiallyChecked:
          return;
      }//switch( m_allNoneSome->checkState() )
      
      
      for( auto w : m_peaks->children() )
      {
        auto p = dynamic_cast<DrfPeak *>( w );
        if( p )
        {
          p->m_useCb->setChecked( useAll );
          if( p->m_userBr )
            p->m_userBr->setHidden( !useAll );
        }
      }//for( auto w : m_peaks->children() )
      
      refreshSourcesVisible();
    }//void handleUserToggleAllNoneSum()
    
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
        if( p )
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
      
      updateTitle();
      
      m_srcInfoUpdated.emit();
    }//void refreshSourcesVisible()
    
    bool isBackground()
    {
      return m_background->isChecked();
    }
    
    vector<DrfPeak *> peaks()
    {
      vector<DrfPeak *> answer;
      for( auto w : m_peaks->children() )
      {
        auto p = dynamic_cast<DrfPeak *>( w );
        if( p )
          answer.push_back( p );
      }
      
      return answer;
    }//vector<DrfPeak *> peaks()
    
    vector<MakeDrfSrcDef *> sources()
    {
      vector<MakeDrfSrcDef *> answer;
      for( auto w : m_sources->children() )
      {
        auto p = dynamic_cast<MakeDrfSrcDef *>( w );
        if( p )
          answer.push_back( p );
      }
      
      return answer;
    }//vector<MakeDrfSrcDef *> sources()
    
    
    /** Returns only peaks selected for use with a source nuclide. */
    vector<pair<DrfPeak *,MakeDrfSrcDef *>> selected_peak_to_sources()
    {
      vector<pair<DrfPeak *,MakeDrfSrcDef *>> answer;
      
      const auto peakWidgets = peaks();
      const auto sourceWidgets = sources();
      
      for( auto peakw : peakWidgets )
      {
        if( !peakw->m_useCb->isChecked() )
          continue;
        
        const auto nuc = peakw->m_peak->parentNuclide();
        
        MakeDrfSrcDef *src = nullptr;
        for( auto sourcew : sourceWidgets )
        {
          if( sourcew->nuclide() == nuc )
          {
            src = sourcew;
            break;
          }
        }//for( auto sourcew : sourceWidgets )
        
        if( src )
          answer.push_back( make_pair(peakw, src) );
      }//for( auto peakw : peakWidgets )
      
      return answer;
    }//vector<pair<DrfPeak *,MakeDrfSrcDef *>> peak_to_sources()
    
    Wt::Signal<> &srcInfoUpdated(){ return m_srcInfoUpdated; };
    
    
    std::shared_ptr<const SpecMeas> measurement() { return m_meas; }
    const set<int> &samples(){ return m_samples; }
    
  };//class DrfSpecFileSample
  
  
  class DrfSpecFile : public WPanel
  {
    std::shared_ptr<const SpecMeas> m_meas;
    WContainerWidget *m_sampleWidgets;
    Wt::Signal<> m_updated;
    
  public:
    DrfSpecFile( std::shared_ptr<const SpecMeas> meas,
                MaterialDB *materialDB,
                Wt::WSuggestionPopup *materialSuggest,
                WContainerWidget *parent = nullptr )
      : WPanel( parent ),
        m_meas( meas ),
        m_sampleWidgets( nullptr )
    {
      addStyleClass( "DrfSpecFile" );
      
      const set<set<int>> sampsWithPeaks = meas->sampleNumsWithPeaks();

      m_sampleWidgets = new WContainerWidget();
      setCentralWidget( m_sampleWidgets );
      
      for( const set<int> &peakSamps : sampsWithPeaks )
      {
        std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaksptr = meas->peaks( peakSamps );
        if( !peaksptr || peaksptr->empty() )
          continue;
        
        auto sample = new DrfSpecFileSample( meas, peakSamps, materialDB, materialSuggest, m_sampleWidgets );
        sample->srcInfoUpdated().connect( std::bind( [this](){ m_updated.emit(); }) );
      }//for( const set<int> &peakSamps : sampsWithPeaks )
    }//DrfSpecFile(...)
    
    Wt::Signal<> &updated() { return m_updated; }
    
    vector<DrfSpecFileSample *> fileSamples()
    {
      vector<DrfSpecFileSample *> answer;
      for( auto w : m_sampleWidgets->children() )
      {
        auto p = dynamic_cast<DrfSpecFileSample *>(w);
        if( p )
          answer.push_back( p );
      }
      return answer;
    }//fileSamples()
    
    std::shared_ptr<const SpecMeas> measurement(){ return m_meas; }
  };//class DrfSpecFile
  
  
  
}//namespace to implement DrfSpecFile and DrfSpecFileSample


MakeDrfWindow::MakeDrfWindow( InterSpec *viewer, MaterialDB *materialDB, Wt::WSuggestionPopup *materialSuggest )
  : AuxWindow( "Create Detector Response Function",
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
  m_chart( nullptr ),
  m_files( nullptr ),
  m_detDiameter( nullptr ),
  m_showFwhmPoints( nullptr ),
  m_fwhmEqnType( nullptr ),
  m_sqrtEqnOrder( nullptr ),
  m_effEqnOrder( nullptr ),
  m_effEqnUnits( nullptr ),
  m_chartLowerE( nullptr ),
  m_chartUpperE( nullptr ),
  m_errorMsg( nullptr ),
  m_intrinsicEffAnswer( nullptr ),
  m_fwhmFitId( 0 ),
  m_effEqnFitId( 0 )
{
  assert( m_interspec );
  assert( m_materialDB );
  assert( m_materialSuggest );
  
  wApp->useStyleSheet( "InterSpec_resources/MakeDrf.css" );
  
  addStyleClass( "MakeDrf" );
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( layout );
  
  WContainerWidget *fitOptionsDiv = new WContainerWidget();
  fitOptionsDiv->addStyleClass( "MakeDrfOptions" );
  layout->addWidget( fitOptionsDiv, 0, 0, 2, 1 );
  
  m_detDiameter = new WLineEdit( fitOptionsDiv );
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_detDiameter->setValidator( distValidator );
  m_detDiameter->setText( "2.54 cm" );
  m_detDiameter->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  m_detDiameter->enterPressed().connect( this, &MakeDrf::handleSourcesUpdates );
  
  
  m_fwhmEqnType = new WComboBox( fitOptionsDiv );
  m_fwhmEqnType->setInline( false );
  m_fwhmEqnType->addItem( "Gadras Eqation" );
  m_fwhmEqnType->addItem( "Sqrt Power Series" );
  m_fwhmEqnType->setCurrentIndex( 0 );
  m_fwhmEqnType->changed().connect( this, &MakeDrf::handleFwhmTypeChanged );
  
  m_sqrtEqnOrder = new WComboBox( fitOptionsDiv );
  m_sqrtEqnOrder->setInline( false );
  m_sqrtEqnOrder->hide();
  m_sqrtEqnOrder->changed().connect( this, &MakeDrf::handleSqrtEqnOrderChange );
  
  
  m_effEqnOrder = new WComboBox( fitOptionsDiv );
  m_effEqnOrder->setInline( false );
  m_effEqnOrder->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  
  m_effEqnUnits = new WComboBox( fitOptionsDiv );
  m_effEqnUnits->setInline( false );
  m_effEqnUnits->addItem( "Int. Eff in keV" );
  m_effEqnUnits->addItem( "Int. Eff in MeV" );
  m_effEqnUnits->setCurrentIndex( 1 );
  m_effEqnUnits->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  
  m_chart = new MakeDrfChart();
  layout->addWidget( m_chart, 0, 1 );
  
  //m_chart->setMinimumSize( 250, 250 );
  const int wh = viewer->renderedHeight();
  const int ww = viewer->renderedWidth();
  const int chartHeight = std::min( wh/3, 350 );
  const int chartWidth = std::min( static_cast<int>(0.75*ww - 150), 800 );
  m_chart->resize( chartWidth, chartHeight );
  
  WContainerWidget *chartOptionsDiv = new WContainerWidget();
  layout->addWidget( chartOptionsDiv, 1, 1 );
  
  m_showFwhmPoints = new WCheckBox( "Show FWHM points", chartOptionsDiv );
  m_showFwhmPoints->setChecked( true );
  m_showFwhmPoints->changed().connect( this, &MakeDrf::handleShowFwhmPointsToggled );
  
  WLabel *label = new WLabel( "Display Energy, Lower:", chartOptionsDiv );
  label->setMargin(15,Wt::Left);
  m_chartLowerE = new WDoubleSpinBox( chartOptionsDiv );
  label->setBuddy( m_chartLowerE );
  
  label = new WLabel( "Upper:", chartOptionsDiv );
  label->setMargin(5,Wt::Left);
  m_chartUpperE = new WDoubleSpinBox( chartOptionsDiv );
  m_chartUpperE->setRange(0, 10000);
  m_chartUpperE->setValue( 3000.0 );
  label->setBuddy( m_chartUpperE );
  
  auto updateChartRange = [this](){
    if( m_chartLowerE->validate() != WValidator::State::Valid
        || m_chartUpperE->validate() != WValidator::State::Valid )
      return;
    m_chart->setXRange(m_chartLowerE->value(), m_chartUpperE->value());
  };
  
  m_chartLowerE->changed().connect( std::bind(updateChartRange) );
  m_chartUpperE->changed().connect( std::bind(updateChartRange) );
  m_chart->xRangeChanged().connect( this, &MakeDrf::chartEnergyRangeChangedCallback );
  
  
  m_files = new WContainerWidget();
  m_files->setMaximumSize( WLength::Auto, std::max((wh - chartHeight - 150), 250) );
  m_files->setOverflow( Overflow::OverflowAuto, Wt::Vertical );
  
  m_errorMsg = new WText( "", Wt::XHTMLText );
  m_errorMsg->addStyleClass( "MakeDrfErrTxt" );
  layout->addWidget( m_errorMsg, 1, 0, 1, 2 );
  m_errorMsg->hide();
  
  
  m_intrinsicEffAnswer = new WText( "", Wt::XHTMLText );
  m_intrinsicEffAnswer->addStyleClass( "MakeDrfEqnAnswer" );
  layout->addWidget( m_intrinsicEffAnswer, 2, 0, 1, 2 );
  m_intrinsicEffAnswer->hide();
  
  layout->addWidget( m_files, 3, 0, 1, 2 );
  
  //layout->setRowStretch( 0, 1 );
  //layout->setRowStretch( 1, 3 );
  //layout->setColumnStretch( 1, 1 );
  
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
    DrfSpecFile *fileWidget = new DrfSpecFile( meas, materialDB, materialSuggest, m_files );
    fileWidget->updated().connect( this, &MakeDrf::handleSourcesUpdates );
  }//for( loop over opened files )
  
  //If we directly call handleSourcesUpdates() now, getting the activity
  //  uncertainty will throw an exception because they wont validate.... whatever.
  WServer::instance()->post( wApp->sessionId(), wApp->bind( boost::bind(&MakeDrf::handleSourcesUpdates,this) ) );
}//MakeDrf( constructor )


MakeDrf::~MakeDrf()
{

}//~MakeDrf()


void MakeDrf::handleSourcesUpdates()
{
  size_t numchan = 0;
  vector< std::shared_ptr<const PeakDef> > peaks;
  vector<MakeDrfFit::DetEffDataPoint> effpoints;
  
  bool detDiamInvalid = false;
  double diameter = 2.54*PhysicalUnits::cm;
  try
  {
    diameter = PhysicalUnits::stringToDistance( m_detDiameter->text().toUTF8() );
  }catch(...)
  {
    detDiamInvalid = true;
  }
  
  //Go through and and grab background peaks
  vector<DrfPeak *> backgroundpeaks;
  for( auto w : m_files->children() )
  {
    auto fileWidget = dynamic_cast<DrfSpecFile *>( w );
    if( !fileWidget )
      continue;
    
    numchan = std::max( numchan, fileWidget->measurement()->num_gamma_channels() );
    
    for( DrfSpecFileSample *sample : fileWidget->fileSamples() )
    {
      if( !sample->isBackground() )
        continue;
      
      for( auto peak : sample->peaks() )
      {
        peak->setBackgroundBeingSubtractedInfo( false, 0.0 );
        if( peak->m_useCb->isChecked() )
          backgroundpeaks.push_back( peak );
      }
    }//for( DrfSpecFileSample *sample : sampleWidgets )
  }//for( auto w : m_files->children()  )
  
  //Go through and make MakeDrfChart::DataPoint points, and also make sure to
  //  indicate which peaks are being background subtracted.
  bool not_all_sources_used = false;
  float minenergy = 999999.0f, maxenergy = -9999999.0f;
  vector<MakeDrfChart::DataPoint> datapoints;
  for( auto w : m_files->children()  )
  {
    auto fileWidget = dynamic_cast<DrfSpecFile *>( w );
    if( !fileWidget )
      continue;
    
    for( DrfSpecFileSample *sample : fileWidget->fileSamples() )
    {
      if( sample->isBackground() )
        continue;
      
      auto meas = sample->measurement();
      for( auto det : meas->detector_names() )
      {
        for( auto samplenum : sample->samples() )
        {
          auto m = meas->measurement( samplenum, det );
          if( m && m->num_gamma_channels() > 7 )
          {
            minenergy = std::min( m->gamma_energy_min(), minenergy );
            maxenergy = std::max( m->gamma_energy_max(), maxenergy );
          }
        }
      }//for( auto det : meas->detector_names() )
      
      
      //Check if we will background subtract
      for( auto pp : sample->selected_peak_to_sources() )
      {
        const std::shared_ptr<const PeakDef> peak = pp.first->m_peak;
        double back_peak_area = 0.0, back_peak_area_uncert = 0.0, back_peak_lt = 0.0;
        for( auto b : backgroundpeaks )
        {
          const std::shared_ptr<const PeakDef> backpeak = b->m_peak;
          if( fabs(backpeak->mean() - peak->mean()) < 1.5*peak->sigma()
             || (backpeak->parentNuclide() && peak->parentNuclide() && (fabs(backpeak->gammaParticleEnergy() - peak->gammaParticleEnergy()) < 1.0)) )
          {
            back_peak_area += backpeak->peakArea();
            //ToDo: Check (when I'm not so tired) to make sure this is the right way to handle uncert.
            back_peak_area_uncert = sqrt( back_peak_area_uncert*back_peak_area_uncert + backpeak->peakAreaUncert()*backpeak->peakAreaUncert() );
            back_peak_lt += b->m_livetime;
          }//
        }//for( auto b : backgroundpeaks )
        
        const bool subBack = (back_peak_area > DBL_EPSILON && back_peak_lt > DBL_EPSILON);
        pp.first->setBackgroundBeingSubtractedInfo( subBack, (subBack ? (back_peak_area/back_peak_lt) : 0.0) );
        
        MakeDrfChart::DataPoint point;
        point.energy = peak->parentNuclide() ? peak->gammaParticleEnergy() : peak->mean();
        point.livetime = pp.first->m_livetime;
        point.peak_area = peak->peakArea();
        point.peak_area_uncertainty = peak->peakAreaUncert();
        
        if( subBack )
        {
          const double frac_back_uncert = back_peak_area_uncert / back_peak_area;
          const double scaled_back = point.livetime * back_peak_area / back_peak_lt;
          
          cout << "Setting backsub for " << point.energy << " keV; point.peak_area=" << point.peak_area << ", scaled_back=" << scaled_back << endl;
          
          point.peak_area -= scaled_back;
          //ToDo: check this uncertainty is actually correct.
          point.peak_area_uncertainty = sqrt( point.peak_area_uncertainty*point.peak_area_uncertainty
                                             + scaled_back*frac_back_uncert*scaled_back*frac_back_uncert );
        }
        
        if( peak->type() == PeakDef::GaussianDefined )
        {
          point.peak_fwhm = peak->fwhm();
          point.peak_fwhm_uncertainty = 2.35482*peak->uncertainty(PeakDef::CoefficientType::Sigma);
        }else
        {
          point.peak_fwhm = point.peak_fwhm_uncertainty = 0.0f;
        }
        
        const auto nuc = pp.second->nuclide();
        
        try
        {
          point.distance = pp.second->distance();
          const double activity = pp.second->activityAtSpectrumTime();
          const double age = pp.second->ageAtSpectrumTime();
          
          double transmittion_factor = 1.0;
          ShieldingSelect *shield = pp.second->shielding();
          if( shield )
          {
            //ToDo: Attenuation calculation not checked!
            double an = 14, ad = 0.0;
            if( shield->isGenericMaterial() )
            {
              if( !shield->atomicNumberEdit()->text().empty()
                 && !shield->arealDensityEdit()->text().empty() )
              {
                an = shield->atomicNumber();
                ad = shield->arealDensity();
              }
            }else
            {
              std::shared_ptr<Material> mat = shield->material();
              if( mat )
              {
                an = mat->massWeightedAtomicNumber();
                ad = mat->density * shield->thickness();
              }//if( mat )
            }//if( shield->isGenericMaterial() ) / else
            
            const double mu = MassAttenuation::massAttenuationCoeficient( an, point.energy );
            transmittion_factor = exp( -mu * ad );
          }//if( shield )
          
          point.source_count_rate = 0.0;
          const double width = 1.25*(peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth());
          
          if( nuc )
          {
            //ToDo: only create one mixture per nuclide (instead of one per peak),
            //  and can probably also access relevant gammas more efficiently.
            SandiaDecay::NuclideMixture mix;
            mix.addAgedNuclideByActivity( nuc, activity, age );
            const auto rates = mix.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
            for( const auto &r : rates )
            {
              if( fabs(r.energy - peak->gammaParticleEnergy()) < width )
                point.source_count_rate += r.numPerSecond;
            }
          
            point.source_count_rate *= transmittion_factor;
          }else
          {
            if( !pp.first->m_userBr || (pp.first->m_userBr->validate()!=WValidator::State::Valid) )
              throw runtime_error( "Logic error - no nuclide associated with " + std::to_string(peak->mean()) + " keV peak, and no user BR." );
            point.source_count_rate += transmittion_factor * pp.first->m_userBr->value() * activity;
          }
          
          //If we get here on initial GUI load, we seem to get error validating the uncertainty (hence a workaround is used to delay first reading)
          point.source_count_rate_uncertainty = point.source_count_rate * pp.second->fractionalActivityUncertainty();
        }catch( std::exception &e )
        {
          cerr << "handleSourcesUpdates: got exception: " << e.what() << endl;
          not_all_sources_used = true;
          continue;
        }
        
        point.peak_color = peak->lineColor();
        point.background_peak_area = back_peak_area;
        point.background_peak_live_time = back_peak_lt;
        
        const double particleEnergy = peak->parentNuclide() ? peak->gammaParticleEnergy() : peak->mean();
        
        char buffer[256] = { '\0' };
        if( nuc )
        {
          snprintf( buffer, sizeof(buffer)-1, "%s %.1f keV, Peak %.1f counts",
                   nuc->symbol.c_str(), particleEnergy, point.peak_area );
        }else
        {
          snprintf( buffer, sizeof(buffer)-1, "No assoc. nuclide, %.1f keV, Peak %.1f counts",
                    particleEnergy, point.peak_area );
        }
        //ToDo: Add more/better source info.
        point.source_information = buffer;
        
        
        {
          const double fracSolidAngle = DetectorPeakResponse::fractionalSolidAngle( diameter, pp.second->distance() );
          const double expected = point.source_count_rate * point.livetime * fracSolidAngle;
          const double eff = point.peak_area / expected;
          
          double fracUncert2 = 0.0;
          if( point.peak_area_uncertainty > 0.0f )
            fracUncert2 += std::pow( point.peak_area_uncertainty / point.peak_area, 2.0f );
          if( point.source_count_rate_uncertainty > 0.0f )
            fracUncert2 += std::pow( point.source_count_rate_uncertainty / point.source_count_rate, 2.0f );
          
          MakeDrfFit::DetEffDataPoint effpoint;
          effpoint.energy = point.energy;
          effpoint.efficiency = eff;
          effpoint.efficiency_uncert = sqrt(fracUncert2);
          
          effpoints.push_back( effpoint );
        }
        
        
        peaks.push_back( peak );
        datapoints.push_back( point );
      }
    }//for( DrfSpecFileSample *sample : sampleWidgets )
  }//for( auto w : m_files->children()  )
  
  
  if( maxenergy < minenergy )
  {
    minenergy = 0.0;
    maxenergy = 3000.0;
  }
  
  m_chartLowerE->setRange( std::min(0.0,minenergy-10.0), maxenergy );
  m_chartUpperE->setRange( minenergy, maxenergy+10 );
  
  string msg;
  if( not_all_sources_used )
    msg += "Some sources had errors.";
  if( detDiamInvalid )
    msg += string(msg.empty() ? "" : "  ") + "Detector diameter is invalid - assumin 2.54 cm.";
  
  m_errorMsg->setText( msg );
  m_errorMsg->setHidden( msg.empty() );
  
  
  const int numPeaks = static_cast<int>( peaks.size() );
  const int currentOrder = m_effEqnOrder->currentIndex();
  const int nCurrentOrders = m_effEqnOrder->count();
  
  if( nCurrentOrders == 8 && numPeaks >= 8 )
  {
    //Dont need to do anything
  }else if( nCurrentOrders != numPeaks )
  {
    m_effEqnOrder->clear();
    for( int order = 0; order < numPeaks && order < 8; ++order )
      m_effEqnOrder->addItem( std::to_string(order+1) + " Fit Params"  );
    
    if( currentOrder > 0 )
    {
      if( (currentOrder < (numPeaks/2)) && (currentOrder < 5) )
        m_effEqnOrder->setCurrentIndex( std::min(numPeaks/2,5) );
      else if( currentOrder >= numPeaks )
        m_effEqnOrder->setCurrentIndex( numPeaks-1 );
      else
        m_effEqnOrder->setCurrentIndex( currentOrder );
    }else if( numPeaks )
    {
      const int index = std::min(numPeaks-1,6); //Default to 6.
      m_effEqnOrder->setCurrentIndex( index );
    }else
      m_effEqnOrder->setCurrentIndex( -1 );
  }//
  
  const int currentSqrtOrder = m_effEqnOrder->currentIndex();
  const int nCurrentSqrtOrders = m_sqrtEqnOrder->count();
  if( nCurrentSqrtOrders==6 && numPeaks>=6 )
  {
    //Dont need to do anything
  }else if( nCurrentSqrtOrders != numPeaks )
  {
    m_sqrtEqnOrder->clear();
    for( int order = 0; order < numPeaks && order < 6; ++order )
      m_sqrtEqnOrder->addItem( "Order " + std::to_string(order) + " FWHM"  );
    
    if( currentSqrtOrder > 0 )
    {
      if( (currentSqrtOrder < (numPeaks/2)) && (currentSqrtOrder < 4) )
        m_sqrtEqnOrder->setCurrentIndex( std::min(numPeaks/2,4) );
      else if( currentSqrtOrder >= numPeaks )
        m_sqrtEqnOrder->setCurrentIndex( numPeaks-1 );
      else
        m_sqrtEqnOrder->setCurrentIndex( nCurrentSqrtOrders );
    }else if( numPeaks )
    {
      const int index = std::min(numPeaks-1,4); //Default to 5 fit params.
      m_sqrtEqnOrder->setCurrentIndex( index );
    }else
      m_sqrtEqnOrder->setCurrentIndex( -1 );
  }//
  
  m_chart->setDataPoints( datapoints, diameter, minenergy, maxenergy );
  
  fitEffEqn( effpoints );
  fitFwhmEqn( peaks, numchan );
}//void handleSourcesUpdates()


void MakeDrf::handleFwhmTypeChanged()
{
  switch( m_fwhmEqnType->currentIndex() )
  {
    case 0:
      m_sqrtEqnOrder->hide();
      break;
      
    case 1:
      m_sqrtEqnOrder->show();
      break;
  }//switch( m_fwhmEqnType->currentIndex() )
  
  //We could update *just* the FWHM equation like:
  //size_t numchan = 0;
  //vector< std::shared_ptr<const PeakDef> > peaks;
  // ... get peaks and numchan ...
  //fitFwhmEqn( peaks, numchan );
  //But instead we'll be lazy and update everything
  
  handleSourcesUpdates();
}//void handleFwhmTypeChanged()


void MakeDrf::handleSqrtEqnOrderChange()
{
  //We could update *just* the FWHM equation like:
  //size_t numchan = 0;
  //vector< std::shared_ptr<const PeakDef> > peaks;
  // ... get peaks and numchan ...
  //fitFwhmEqn( peaks, numchan );
  //But instead we'll be lazy and update everything
  
  handleSourcesUpdates();
}//void handleSqrtEqnOrderChange()


void MakeDrf::handleShowFwhmPointsToggled()
{
  m_chart->showFwhmPoints( m_showFwhmPoints->isChecked() );
}//void handleShowFwhmPointsToggled()


void MakeDrf::chartEnergyRangeChangedCallback( double lower, double upper )
{
  m_chartLowerE->setValue( lower );
  m_chartUpperE->setValue( upper );
}


void MakeDrf::fitFwhmEqn( std::vector< std::shared_ptr<const PeakDef> > peaks,
                          const size_t num_gamma_channels )
{
  ++m_fwhmFitId;
  
  const int fitid = m_fwhmFitId;
  const string sessionId = wApp->sessionId();
  
  m_fwhmCoefs.clear();
  m_fwhmCoefUncerts.clear();
  m_chart->setFwhmCoefficients( vector<float>{}, vector<float>{},
                                MakeDrfChart::FwhmCoefType::Gadras,
                                MakeDrfChart::EqnEnergyUnits::keV );

  int sqrtEqnOrder = -1;
  DetectorPeakResponse::ResolutionFnctForm fnctnlForm;
  
  switch( m_fwhmEqnType->currentIndex() )
  {
    case 0:
      fnctnlForm = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
      break;
      
    case 1:
      fnctnlForm = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      sqrtEqnOrder = m_sqrtEqnOrder->currentIndex() + 1;
      break;
      
    default:
      return;
  }//switch( m_fwhmEqnType->currentIndex() )
  
  
  //ToDo: I'm not entirely sure the next line protects against updateFwhmEqn()
  //  not being called if this widget is deleted before fit is done.
  auto updater = boost::bind( &MakeDrf::updateFwhmEqn, this, _1, _2, static_cast<int>(fnctnlForm), fitid );
  
  auto worker = [sessionId,fnctnlForm,peaks,num_gamma_channels,sqrtEqnOrder,updater]() {
    try
    {
      auto peakdequ = std::make_shared<std::deque< std::shared_ptr<const PeakDef> > >( peaks.begin(), peaks.end() );
    
      //Takes between 5 and 500ms for a HPGe detector
      const double start_time = UtilityFunctions::get_wall_time();
      vector<float> fwhm_coefs, fwhm_coefs_uncert;
      MakeDrfFit::performResolutionFit( peakdequ, num_gamma_channels, fnctnlForm, sqrtEqnOrder, fwhm_coefs, fwhm_coefs_uncert );
    
      const double end_time = UtilityFunctions::get_wall_time();
    
      assert( fwhm_coefs.size() == fwhm_coefs_uncert.size() );
      cout << "Fit FWHM: {";
      for( size_t i = 0; i < fwhm_coefs.size(); ++i )
        cout << fwhm_coefs[i] << "+-" << fwhm_coefs_uncert[i] << ", ";
      cout << "}; took " << (end_time-start_time) << " seconds" << endl;
      
      WServer::instance()->post( sessionId, std::bind( [updater,fwhm_coefs,fwhm_coefs_uncert](){
        updater(fwhm_coefs,fwhm_coefs_uncert);
      } ) );
    }catch( std::exception &e )
    {
      cout << "Failed to fit FWHM coefs: " << e.what() << endl;
    }//try / catch fit FWHM
  };
  
  WServer::instance()->ioService().post( worker );
}//void fitFwhmEqn( std::vector< std::shared_ptr<const PeakDef> > peaks )


void MakeDrf::updateFwhmEqn( std::vector<float> coefs,
                             std::vector<float> uncerts,
                             const int functionalForm,
                             const int fitid )
{
  if( fitid != m_fwhmFitId )
    return;
  
  MakeDrfChart::FwhmCoefType eqnType = MakeDrfChart::FwhmCoefType::Gadras;
  const auto fcntform = DetectorPeakResponse::ResolutionFnctForm( functionalForm );
  switch( fcntform )
  {
    case DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn:
      eqnType = MakeDrfChart::FwhmCoefType::Gadras;
      m_fwhmEqnType->setCurrentIndex( 0 );
      break;
    
    case DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial:
      eqnType = MakeDrfChart::FwhmCoefType::SqrtEqn;
      m_fwhmEqnType->setCurrentIndex( 1 );
      break;
      
    case DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm:
    default:
      m_fwhmEqnType->setCurrentIndex( -1 );
      return;
  }//switch( fcntform )
  
  m_chart->setFwhmCoefficients( coefs, uncerts, eqnType, MakeDrfChart::EqnEnergyUnits::keV );
  
  wApp->triggerUpdate();
}//void updateFwhmEqn(...)


void MakeDrf::fitEffEqn( std::vector<MakeDrfFit::DetEffDataPoint> data )
{
  ++m_effEqnFitId;
  
  const int fitid = m_effEqnFitId;
  const string sessionId = wApp->sessionId();
  
  m_effEqnCoefs.clear();
  m_effEqnCoefUncerts.clear();
  m_chart->setEfficiencyCoefficients( m_effEqnCoefs, m_effEqnCoefUncerts, MakeDrfChart::EqnEnergyUnits::keV );
  m_intrinsicEffAnswer->setText( "" );
  
  const int eqnOrderIndex = m_effEqnOrder->currentIndex();
  if( data.empty() )
    return;
  
  if( eqnOrderIndex < 0 )
    return;
  
  const int nfitpars = std::min( eqnOrderIndex+1, 8 );
  
  
  const bool inMeV = (m_effEqnUnits->currentIndex() == 1);
  if( inMeV )
  {
    for( auto &a : data )
      a.energy /= 1000.0f;
  }//if( inMeV )
  
  //ToDo: I'm not entirely sure the next line protects against updateEffEqn()
  //  not being called if this widget is deleted before fit is done.
  auto updater = boost::bind( &MakeDrf::updateEffEqn, this, _1, _2, fitid, _3 );
  
  auto worker = [sessionId,data,nfitpars,updater]() {
    try
    {
      //Takes between 5 and 500ms for a HPGe detector
      const double start_time = UtilityFunctions::get_wall_time();
      vector<float> result, uncerts;
      MakeDrfFit::performEfficiencyFit( data, nfitpars, result, uncerts );
      
      const double end_time = UtilityFunctions::get_wall_time();
      
      assert( result.size() == uncerts.size() );
      cout << "Fit Eff: {";
      for( size_t i = 0; i < result.size(); ++i )
        cout << result[i] << "+-" << uncerts[i] << ", ";
      cout << "}; took " << (end_time-start_time) << " seconds" << endl;
      
      WServer::instance()->post( sessionId, std::bind( [updater,result,uncerts](){
        updater( result, uncerts, string("") );
      } ) );
    }catch( std::exception &e )
    {
      const string errmsg = e.what();
      cout << "Failed to fit intrinsic eff coefs: " << errmsg << endl;
      WServer::instance()->post( sessionId, std::bind( [updater,errmsg](){
        updater( vector<float>(), vector<float>(), errmsg );
      } ) );
    }//try / catch fit FWHM
  };
  
  WServer::instance()->ioService().post( worker );
}//void fitEffEqn( std::vector<MakeDrfFit::DetEffDataPoint> data )


void MakeDrf::updateEffEqn( std::vector<float> coefs, std::vector<float> uncerts, const int fitid, const string errmsg )
{
  const bool isMeV = (m_effEqnUnits->currentIndex()==1);
  const auto units = (isMeV ? MakeDrfChart::EqnEnergyUnits::MeV : MakeDrfChart::EqnEnergyUnits::keV);

  m_effEqnCoefs = coefs;
  m_effEqnCoefUncerts = uncerts;
  m_chart->setEfficiencyCoefficients( coefs, uncerts, units );
  
  if( !errmsg.empty() )
  {
    m_intrinsicEffAnswer->setText( "" );
    m_intrinsicEffAnswer->setHidden( true );
    string errormsg = m_errorMsg->text().toUTF8();
    const bool hadErrorMsg = !errormsg.empty();
    if( hadErrorMsg )
      errormsg = "<div>" + errormsg + "</div><div>";
    errormsg += errormsg;
    if( hadErrorMsg )
      errormsg = "</div>";
    m_errorMsg->setText( WString::fromUTF8(errormsg) );
    m_errorMsg->setHidden( false );
  }
  
  if( coefs.size() > 1 )
  {
    string eqn = "Eff<sub>int.</sub>(x) = exp( ";
    for( size_t i = 0; i < coefs.size(); ++i )
    {
      const float val = fabs(coefs[i]);
      char buffer[64] = {'\0'};
      
      //We will print to 5 significant digits, AFTER the decimal place since
      //  this is a sum of terms.
      const string frmtstr = (val > 1.0) ? ("%." + std::to_string( static_cast<int>(std::ceil(std::log10(val)))+5 ) + "g") : string("%.5g");
      snprintf( buffer, sizeof(buffer), frmtstr.c_str(), val );
      
      eqn += (coefs[i]>=0.0) ? "+" : "-";
      eqn += i ? " " : "";
      eqn += buffer;
      eqn += i ? "*log(x)" : "";
      eqn += (i > 1) ? ("<sup>" + std::to_string(i) + "</sup> ") : string(" ");
    }
    eqn += ")";
    
    //eqn += " (x in ";
    //eqn += (isMeV ? "MeV)" : "keV)");
    
    m_intrinsicEffAnswer->setText( eqn );
    m_intrinsicEffAnswer->setHidden( false );
  }else
  {
    m_intrinsicEffAnswer->setText( "" );
    m_intrinsicEffAnswer->setHidden( true );
  }//if( coefs.size() > 1 ) / else
  
  wApp->triggerUpdate();
}//void updateEffEqn(...)
