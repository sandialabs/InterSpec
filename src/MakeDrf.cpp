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
#include <regex>
#include <deque>
#include <fstream>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WPanel>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WAnchor>
#include <Wt/WServer>
#include <Wt/WResource>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WSvgImage>
#include <Wt/WGroupBox>
#include <Wt/WTemplate>
#include <Wt/WTextArea>
#include <Wt/WTableCell>
#include <Wt/WTabWidget>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WPopupWidget>
#include <Wt/Http/Response>
#include <Wt/WDoubleSpinBox>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>

#include "InterSpec/MakeDrf.h"
#include "InterSpec/SpecMeas.h"
#include "SpecUtils/DateTime.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfSelect.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/MakeDrfChart.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"


// The regex in GCC 4.8.x does not have working regex
#if( defined(__GLIBCXX__) && (__cplusplus < 201402L) )
static_assert( defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) \
              || defined(_GLIBCXX_REGEX_STATE_LIMIT) \
              || (defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE > 4), "GCC 4.8 is not supported due to buggy regex implementation" );
#endif

using namespace std;
using namespace Wt;

namespace
{
  const float ns_NaI3x3IntrinsicEff = 0.47096f; //linear interpolation based on Efficiency.csv for generic 3x3. So could be improved...

  bool source_info_from_lib_file( string srcname, const string &filename,
                                  double &activity, boost::posix_time::ptime &activityDate, string &comments )
  {
    cout << "Will try to find '" << srcname << "' from " << filename << endl;
#ifdef _WIN32
    const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
    ifstream file( wfilename.c_str(), ios::in | ios::binary );
#else
    ifstream file( filename.c_str(), ios::in | ios::binary );
#endif
    if( !file )
      return false;
   
    SpecUtils::trim( srcname );
    
    string line;
    while( SpecUtils::safe_get_line( file, line) )
    {
      vector<string> fields;
      SpecUtils::split( fields, line, " \t" );
      if( fields.size() < 3 )
        continue;
      if( !SpecUtils::iequals_ascii(fields[0], srcname) )
        continue;
      try
      {
        activity = stod(fields[1]) * PhysicalUnits::bq;
      }catch(...)
      {
        continue;
      }
      
      activityDate = SpecUtils::time_from_string( fields[2] .c_str());
      comments = "";
      for( size_t i = 3; i < fields.size(); ++i )
        comments += ((i==3) ? "" : " ") + fields[i];
      
      return true;
    }//while( SpecUtils::safe_get_line( file, line) )
    
    return false;
  }//source_info_from_lib_file(...)
  
  
  //ToDo: this function is also implemented in SpecFileQueryWidget.cpp - should
  //      put in the same place
  void output_csv_field( std::ostream &out, std::string s)
  {
    //adapted from http://darrendev.blogspot.com/2009/11/escaping-csv-in-c.html 20160623
    //There are two escaping rules for each field in a comma-separated value row:
    //  1. Change each double quote to two double quotes.
    //  2. Surround with double quotes if the field contains a comma or double quote.
    
    //ToDo: see if this agrees with boost::io::quoted(...)
    
    if( s.find('"') != std::string::npos ) //Escape double-quotes
    {
      std::string::size_type pos = 0;
      while( 1 )
      {
        pos = s.find( '"', pos );
        if( pos == std::string::npos )
          break;
        s.replace(pos,1,"\"\"");
        pos += 2; //Need to skip over those two quotes, to avoid an infinite loop!
      }
      
      out << '"' << s << '"';
    }else if( s.find(',') != std::string::npos
             || s.find('\n') != std::string::npos
             || s.find('\r') != std::string::npos )  //Need to surround with "..."
    {
      out << '"' << s << '"';
    }else
      out << s;
  }//output_csv_field
  
  
  
  /** Class to download either a N42-2012, or PCF file. */
  class CalFileDownloadResource : public Wt::WResource
  {
    const bool m_pcf;
    MakeDrf * const m_makedrf;
    string m_filename;
    Wt::WApplication *m_app;
    
  public:
    CalFileDownloadResource( const bool pcf, MakeDrf *parent )
     : WResource( parent ), m_pcf(pcf), m_makedrf( parent ), m_filename( "" ),
       m_app( WApplication::instance() )
    {
      assert( m_app );
    }
    
    virtual ~CalFileDownloadResource()
    {
      beingDeleted();
    }
    
    void setSuggestFileName( string filename )
    {
      m_filename = filename;
    }
    
    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in CalFileDownloadResource.";
        
        response.out() << "Error grabbing application lock to form CalFileDownloadResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        
        return;
      }//if( !lock )
      
      if( !m_makedrf )
        return;
      
      shared_ptr<SpecMeas> calfile = m_makedrf->assembleCalFile();
      if( !calfile )
        return;
      
      string filename = m_filename.empty() ? calfile->filename() : m_filename;
      if( filename.empty() )
        filename = "MakeDrfData";
      
      //Remove bad filename characters
      const string notallowed = "\\/:?\"<>|*";
      for( auto it = begin(filename) ; it < end(filename) ; ++it )
      {
        if( notallowed.find(*it) != string::npos )
          *it = ' ';
      }
      
      const char * const extension = (m_pcf ? ".pcf" : ".n42");
      if( !SpecUtils::iequals_ascii( SpecUtils::file_extension(filename), extension ) )
        filename += extension;
      
      suggestFileName( filename, WResource::Attachment );
      response.setMimeType( "application/octet-stream" );
      
      if( m_pcf)
        calfile->write_pcf( response.out() );
      else
        calfile->write_2012_N42( response.out() );
    }
  };//class CalFileDownloadResource
  
  
  float effEqnUncert( const float energy, vector<float> coefs, const vector<float> &uncerts )
  {
    assert( coefs.size() == uncerts.size() );
    
    const float eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, coefs );
    
    double neguncert = 0.0, posuncert = 0.0;
    
    for( size_t i = 0; i < uncerts.size(); ++i )
    {
      const float orig = coefs[i];
      coefs[i] = orig + uncerts[i];
      const double plus = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, coefs ) - eff;
      coefs[i] = orig - uncerts[i];
      const double minus = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, coefs ) - eff;
      coefs[i] = orig;
      if( plus > 0.0 )
        posuncert += plus*plus;
      else
        neguncert += plus*plus;
      if( minus > 0.0 )
        posuncert += minus*minus;
      else
        neguncert += minus*minus;
    }//for( size_t i = 0; i < uncerts.size(); ++i )
    
    posuncert = sqrt(posuncert);
    neguncert = sqrt(neguncert);
    
    if( IsNan(posuncert) || IsInf(posuncert)
       || IsNan(neguncert) || IsInf(neguncert) )
      return -999.9f;
    
    return 0.5f*(posuncert + neguncert);
  }//float effEqnUncert(...)
  
  
  
  /** Class to downlaod CSV that contains info fit for. */
  class DrfSummaryBase : public Wt::WResource
  {
  protected:
    string m_filename;
    string m_description;
    MakeDrf * const m_makedrf;
    Wt::WApplication *m_app;
    
  public:
    DrfSummaryBase( MakeDrf *parent )
    : WResource( parent ), m_filename(""), m_makedrf( parent ), m_app( WApplication::instance() )
    {
      assert( m_app );
    }
    
    virtual ~DrfSummaryBase()
    {
      beingDeleted();
    }
    
    void setSuggestFileName( string filename )
    {
      m_filename = filename;
      
      if( filename.empty() )
        filename = "drf";
      filename += "_";
      filename += fileEnding();
      suggestFileName( filename, WResource::Attachment );
    }
    
    void setDescription( const string &desc )
    {
      m_description = desc;
    }
    
    virtual const char *fileEnding() const = 0; //"drfinfo.csv"
    
    virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response ) = 0;
  };//class DrfSummaryBase


  class CsvDrfDownloadResource : public DrfSummaryBase
  {
  public:
    CsvDrfDownloadResource( MakeDrf *parent )
    : DrfSummaryBase( parent )
    {
      suggestFileName( fileEnding(), WResource::Attachment );
    }
    
    virtual const char *fileEnding() const { return "drfinfo.csv"; };
    
    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in CsvDrfDownloadResource.";
        
        response.out() << "Error grabbing application lock to form CsvDrfDownloadResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        
        return;
      }//if( !lock )
      
      response.setMimeType( "text/csv" );
      if( m_makedrf )
        m_makedrf->writeCsvSummary( response.out(), m_filename, m_description );
    }//void handleRequest(...)
  };//class CsvDrfDownloadResource
  
  

  /** Class to download CSV that contains info fit for. */
  class RefSheetDownloadResource : public DrfSummaryBase
  {
  
  public:
    RefSheetDownloadResource( MakeDrf *parent )
     : DrfSummaryBase( parent )
    {
      suggestFileName( fileEnding(), WResource::Attachment );
    }
  
    virtual const char *fileEnding() const { return "ref_sheet.html"; };
  
    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
    
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in RefSheetDownloadResource.";
      
        response.out() << "Error grabbing application lock to form RefSheetDownloadResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
      
        return;
      }//if( !lock )
    
      response.setMimeType( "text/csv" );
      if( m_makedrf )
        m_makedrf->writeRefSheet( response.out(), m_filename, m_description );
    }//void handleRequest(...)
  };//class RefSheetDownloadResource
  

  class DrfPeak : public WContainerWidget
  {
  public:
    const std::shared_ptr<const PeakDef> m_peak;
    const double m_livetime;
    Wt::WCheckBox *m_useForEffCb;
    WText *m_descTxt;
    WText *m_backSubTxt;
    WDoubleSpinBox *m_userBr;
    WPushButton *m_previewBtn;
    WPopupWidget *m_previewPopup;
    std::shared_ptr<const SpecUtils::Measurement> m_summed_meas;
    Wt::Signal<DrfPeak *> m_peakPreviewShow;
    bool m_isBackground;

    
    DrfPeak( std::shared_ptr<const PeakDef> peak, double livetime,
             std::shared_ptr<const SpecUtils::Measurement> summed_meas,
             WContainerWidget *parent = nullptr )
    : WContainerWidget( parent ),
      m_peak( peak ),
      m_livetime( livetime ),
      m_useForEffCb( nullptr ),
      m_descTxt( nullptr ),
      m_backSubTxt( nullptr ),
      m_userBr( nullptr ),
      m_previewBtn( nullptr ),
      m_previewPopup( nullptr ),
      m_summed_meas( summed_meas ),
      m_isBackground( false )
    {
      addStyleClass( "DrfPeak" );
      m_useForEffCb = new WCheckBox( "Use", this );
      m_useForEffCb->addStyleClass( "DrfPeakUseCb" );
      
      char buffer[512];
      if( peak && peak->parentNuclide() && peak->nuclearTransition() )
      {
        const bool use = peak->useForDrfIntrinsicEffFit();
        
        const char *gammatype = "";
        switch( peak->sourceGammaType() )
        {
          case PeakDef::SourceGammaType::NormalGamma:
            break;
          case PeakDef::SourceGammaType::AnnihilationGamma:
            gammatype = " (annih.)";
            break;
          case PeakDef::SourceGammaType::DoubleEscapeGamma:
            gammatype = " (D.E.)";
            break;
          case PeakDef::SourceGammaType::SingleEscapeGamma:
            gammatype = " (S.E.)";
            break;
          case PeakDef::SourceGammaType::XrayGamma:
            gammatype = " (x-ray)";
            break;
        }//switch( peak->sourceGammaType() )
        
        m_useForEffCb->setChecked( use );
        snprintf( buffer, sizeof(buffer), "%s: %.2f keV peak with %.1f cps for %.2f keV gamma%s.",
                  peak->parentNuclide()->symbol.c_str(),
                  peak->mean(),
                  (peak->amplitude() / livetime),
                  peak->gammaParticleEnergy(),
                  gammatype );
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
        
        m_useForEffCb->setUnChecked();
        m_userBr->hide();
        m_useForEffCb->checked().connect( std::bind( [this](){
          m_userBr->setHidden( !m_useForEffCb->isChecked() || m_isBackground );
        } ) );
        m_useForEffCb->unChecked().connect( std::bind( [this](){
          m_userBr->setHidden( !m_useForEffCb->isChecked() || m_isBackground );
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
    
    bool useForEffFit() const
    {
      return m_useForEffCb->isChecked();
    }
    
    void setUse( const bool use )
    {
      m_useForEffCb->setChecked( use );
      if( m_userBr )
        m_userBr->setHidden( !use || m_isBackground );
   
#ifndef _WIN32
#warning "Use peak model (after modding it) to set whether to use intrinsic eff or not"
#endif
      //m_peakModel->setUseForDrfIntrinsicEffFit( use );
      //m_peakModel->setData( )
    }
    
    void setIsBackground( const bool back )
    {
      m_isBackground = back;
      
      //if( m_userBr && back )
      //  m_userBr->setValue( 1.0 );
      
      if( m_userBr )
        m_userBr->setHidden( !m_useForEffCb->isChecked() || m_isBackground );
    }//void setIsBackground( const bool back )
    
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
    
    void hidePeakPreview()
    {
      if( m_previewPopup && !m_previewPopup->isHidden() )
      {
        m_previewBtn->removeStyleClass( "active" );
        m_previewPopup->setHidden( true );
      }
    }//void hidePeakPreview()
    
    void togglePeakPreview()
    {
      //WAnimation animation(WAnimation::AnimationEffect::Pop, WAnimation::Linear, 250 );
      if( m_previewPopup )
      {
        if( m_previewPopup->isHidden() )
        {
          m_previewBtn->addStyleClass( "active" );
          m_previewPopup->setHidden( false );
          m_peakPreviewShow.emit( this );
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
        
        //ToDo: add a "X" in the corner so users can close the showing window.
      }catch( std::exception & )
      {
        WText *msg = new WText( "Unable to generate preview" );
        m_previewPopup = new WPopupWidget( msg, this );
      }//try / catch make a chart
      
      m_previewPopup->setAnchorWidget( m_previewBtn, Wt::Orientation::Vertical );
      m_previewPopup->setJavaScriptMember("wtNoReparent", "true");
      m_previewPopup->show();
      m_previewBtn->addStyleClass( "active" );
      m_peakPreviewShow.emit( this );
    }//void togglePeakPreview()
    
  };//class DrfPeak
  
  
  class DrfSpecFileSample : public WPanel
  {
    std::shared_ptr<SpecMeas> m_meas;
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
          shared_ptr<const SpecUtils::Measurement> m = m_meas->measurement( sample, detname );
          if( m && !m->title().empty() )
            title = m->title();
        }//
      }//for( const int sample : samples )
      
      
      string filename = m_meas->filename();
      filename = SpecUtils::filename( filename );
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
        npeaksused += ((p && p->useForEffFit()) ? 1 : 0);
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
      
      for( auto w : m_peaks->children() )
      {
        if( auto p = dynamic_cast<DrfPeak *>( w ) )
        {
          if( p->m_userBr )
            p->setIsBackground( m_background->isChecked() );
        }//
      }//for( auto w : m_peaks->children() )
      
      m_srcInfoUpdated.emit();
    }//void isBackgroundToggled()
    
    
  public:
    DrfSpecFileSample( std::shared_ptr<SpecMeas> meas, set<int> samples,
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
      m_background->checked().connect( this, &DrfSpecFileSample::isBackgroundToggled );
      m_background->unChecked().connect( this, &DrfSpecFileSample::isBackgroundToggled );
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
          shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, detname );
          if( m && !m->title().empty() )
            title = m->title();
          if( m && m->num_gamma_channels() > 5 )
            livetime += m->live_time();
        }//
      }//for( const int sample : samples )
      
      std::shared_ptr<const SpecUtils::Measurement> summed_meas;
      if( samples.size()==1 && detnames.size()==1 )
        summed_meas = meas->measurement( *begin(samples), detnames[0] );
      if( !summed_meas )
        summed_meas = meas->sum_measurements( samples, meas->detector_names(), nullptr );
      
      //ToDo: inprinciple summed_meas->live_time() should equal livetime; should check
      
      int npeaks = 0;
      for( auto peak : *peaks )
      {
        DrfPeak *p = new DrfPeak( peak, livetime, summed_meas, m_peaks );
        ++npeaks;
        p->m_useForEffCb->changed().connect( this, &DrfSpecFileSample::refreshSourcesVisible );
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
        m_allNoneSome->checked().connect( this, &DrfSpecFileSample::handleUserToggleAllNoneSum );
        m_allNoneSome->unChecked().connect( this, &DrfSpecFileSample::handleUserToggleAllNoneSum );
      }
      
      
      refreshSourcesVisible();
      
      //See if the file conforms to GADRAS conventions
      //if( samples.size() == 1 )
      {
        string shielding;
        bool is_background = false;
        double distance = -1.0, activity = -1.0, age_at_meas = -1.0;
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
              is_background |= (m->source_type()==SpecUtils::SourceType::Background);
              is_background |= SpecUtils::icontains( spectitle, "back" );
              is_background |= SpecUtils::icontains( spectitle, "bgr" );
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
                SpecUtils::trim( dist );
                pos = dist.find_first_of( ",@" ); //ToDo: make finding the end of the distance more robust - like with a regex
                if( pos != string::npos )
                  dist = dist.substr(0, pos);
                distance = PhysicalUnits::stringToDistance( dist );
              }catch(...){}
            }//if( pos != string::npos )
            
            //make work for strings like "133BA_something_something"
            const size_t first_space = spectitle.find_first_of( " \t_" );
            //ToDo: make this more robost for finding end of the nuclide, like requiring both numbers and letters.
            
            nuc = db->nuclide( spectitle.substr(0,first_space) );
            
            vector<string> remarks = m->remarks();
            if( !spectitle.empty() )
              remarks.push_back( "Source:" + spectitle );
            
            
            for( string remark : remarks )
            {
              if( !SpecUtils::istarts_with(remark, "Source:") )
                continue;
              
              remark = remark.substr(7);
              
              SpecUtils::trim( remark );
              
              //Look for a string like "Age= 13y 5d 3s something something something"
              std::smatch mtch;
              std::regex expr( string(".+(Age\\s*\\=\\s*(") + PhysicalUnits::sm_timeDurationRegex + ")).*?", std::regex::icase );
                
              if( std::regex_match( remark, mtch, expr ) )
              {
                try
                {
                  age_at_meas = PhysicalUnits::stringToTimeDurationPossibleHalfLife( mtch[2].str(), (nuc ? nuc->halfLife : -1.0) );
                }catch( std::exception & )
                {
                  cerr << "Failed to convert '" << mtch[2].str() << "' to an age." << endl;
                }
               
                const string str_to_remove = mtch[1].str();
                const size_t removepos = remark.find( str_to_remove );
                if( removepos != string::npos )
                {
                  remark = remark.substr(0,removepos) + " " + remark.substr( removepos + str_to_remove.length() );
                  cout << "After removing matching substr, string = '" << remark << "'" << endl;
                }else
                {
                  //This shouldnt ever happen right?
                  cerr << "Failed to find regex group '" << str_to_remove << "' in '" << remark << "' - ignoring" << endl;
                }
              }//if( found "age = ..." )
              
              const size_t openCurlyPos = remark.find( '{' );
              if( shielding.empty() && (openCurlyPos != string::npos) )
              {
                const size_t closeCurlyPos = remark.find( '}', openCurlyPos );
                if( (closeCurlyPos != string::npos) && ((openCurlyPos+4)<closeCurlyPos) )
                {
                  shielding = remark.substr(openCurlyPos+1, closeCurlyPos-openCurlyPos-1);
                  remark = remark.substr(0,openCurlyPos);
                  SpecUtils::trim(remark);
                }
              }//if( source had shielding defined )
              
              /*
               //Untested regex quivalent (or maybe a bit stricter) of the above.
               std::smatch shieldmtch;
               std::regex expr( ".+({\\s*[+]?[0-9]*\\.?[0-9]+\\s*,\\s*[+]?[0-9]*\\.?[0-9]+\\s*}).*?" );  //could be optimized, and extended to arbitrary number of floats
               if( std::regex_match( remark, shieldmtch, expr ) )
               {
                 shielding = shieldmtch[1].str();
                 const size_t removepos = remark.find( shielding );
                 if( removepos != string::npos )
                   remark = remark.substr(0,removepos) + " " + remark.substr( removepos + shielding.length() );
               }
               */
              
              if( !nuc )
                nuc = db->nuclide( remark );  //Doesnt look like this gives a nuclide if there is anything invalid trailing the nuclide (ex, 'Ba133,10uCi')
              
              const size_t commapos = remark.find(',');
              const size_t underscorpos = remark.find('_'); //For a source database name like "133BA_<serial number>"
              
              if( !nuc && commapos!=string::npos )
                nuc = db->nuclide( remark.substr(0,commapos) );
              
              if( !nuc && underscorpos!=string::npos )
                nuc = db->nuclide( remark.substr(0,underscorpos) );
              
              if( !nuc )
              {
                //ToDo: come up with a Regex that will find any reasonable nuclide "Cs137, Cs-137, 137Cs, Cesium 137, U 238m 238mU, etc."
                const size_t first_space = spectitle.find_first_of( " \t_" );
                if( first_space != string::npos )
                  nuc = db->nuclide( remark.substr(0,first_space) );
              }
              
              if( (commapos < underscorpos)
                 && (nuc == db->nuclide(remark.substr(0,commapos))) )
              {
                //We have something like "133Ba,10uCi"
                //ToDo: There may be multiple nuclides specified, for example:
                //      "232U,10uC{26,10}+137Cs,3.2mC{13,5}" - should handle
                try
                {
                  string activitystr = remark.substr(commapos+1);
                  SpecUtils::trim( activitystr );
                  
                  //search for a positive decimal number foloowed by some letters; take fir occirance
                  std::smatch mtch;
                  std::regex expr( "(\\+?\\s*((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?)\\s*([a-zA-Z \\-]+)" );
                  if( std::regex_search( activitystr, mtch, expr ) )
                    activitystr = SpecUtils::trim_copy( mtch[0] );
                  
                  activity = PhysicalUnits::stringToActivity( activitystr );
                }catch( std::exception & )
                {
                }
              }//
              
              if( (activity < 0.0) && (underscorpos != string::npos)
                 && (nuc == db->nuclide(remark.substr(0,underscorpos)) ))
              {
                //- Check for "*Source.lib" file in application directory,
                //  and if found, scan through for the source identified in previous line, and put in date
                vector<string> base_paths{ InterSpec::staticDataDirectory(), SpecUtils::get_working_path() };
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__)) ) )
                try{ base_paths.push_back( InterSpec::writableDataDirectory() ); } catch(...){}
#endif
                vector<string> source_lib_files;
                for( const string &path : base_paths )
                {
                  const vector<string> files = SpecUtils::recursive_ls(path, "Source.lib" );
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
              
              if( age_at_meas >= 0.0 )
                src->setAgeAtMeas( age_at_meas );
              
              if( !shielding.empty() )
              {
                try
                {
                  vector<float> an_ad;
                  SpecUtils::split_to_floats( shielding, an_ad );
                  
                  if( an_ad.size() >= 2
                     && an_ad[0] >= 1.0f && an_ad[0] <= 100.0f
                     && an_ad[1] >= 0.0f && an_ad[1] <= GammaInteractionCalc::sm_max_areal_density_g_cm2 )
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
        
          for( auto w : m_peaks->children() )
          {
            if( auto p = dynamic_cast<DrfPeak *>( w ) )
              p->setIsBackground( is_background );
          }//for( auto w : m_peaks->children() )
        }//if( is_background )
      }//if( samples.size() == 1 )
    }//DrfSpecFileSample constructor
    
    void handleUserToggleAllNoneSum()
    {
      assert( m_allNoneSome );
      
      bool useAll = true;
      switch( m_allNoneSome->checkState() )
      {
        case Wt::Checked:          useAll = true;  break;
        case Wt::Unchecked:        useAll = false; break;
        case Wt::PartiallyChecked: useAll = false; break;
      }//switch( m_allNoneSome->checkState() )
      
      for( auto w : m_peaks->children() )
      {
        auto p = dynamic_cast<DrfPeak *>( w );
        if( p )
          p->setUse( useAll );
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
        if( p && p->useForEffFit() )
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
        if( !peakw->useForEffFit() )
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
    std::shared_ptr<SpecMeas> m_meas;
    WContainerWidget *m_sampleWidgets;
    Wt::Signal<> m_updated;
    
  public:
    DrfSpecFile( std::shared_ptr<SpecMeas> meas,
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
      
      int nsamples = 0;
      for( const set<int> &peakSamps : sampsWithPeaks )
      {
        std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaksptr = meas->peaks( peakSamps );
        if( !peaksptr || peaksptr->empty() )
          continue;
        
        ++nsamples;
        auto sample = new DrfSpecFileSample( meas, peakSamps, materialDB, materialSuggest, m_sampleWidgets );
        sample->srcInfoUpdated().connect( std::bind( [this](){ m_updated.emit(); }) );
      }//for( const set<int> &peakSamps : sampsWithPeaks )
      
      if( nsamples < 2 )
      {
        addStyleClass( "SingleSample" );
      }else
      {
        auto t = WString::fromUTF8( meas->filename() );
        if( !t.empty() )
          setTitle( t );
      }
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
  
  /** MakeDrfChart need resize() explicitly called when its size changes, and
     its not happening when placed directly into the WGridLayout, so this class
     just does that explicitly
   */
  class DrfChartHolder : public WContainerWidget
  {
    MakeDrfChart *m_chart;
    
  public:
    DrfChartHolder( MakeDrfChart *chart, WContainerWidget *parent )
    : WContainerWidget( parent ),
      m_chart( chart )
    {
      setLayoutSizeAware( true );
      addWidget( m_chart );
    }
    
    virtual ~DrfChartHolder(){};
    
    virtual void layoutSizeChanged(int width, int height)
    {
      m_chart->resize( width, height );
    }
  };//class DrfChartHolder
  
}//namespace to implement DrfSpecFile and DrfSpecFileSample


MakeDrf::MakeDrf( InterSpec *viewer, MaterialDB *materialDB,
                  Wt::WSuggestionPopup *materialSuggest,
                  Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_materialDB( materialDB ),
  m_materialSuggest( materialSuggest ),
  m_intrinsicEfficiencyIsValid( this ),
  m_finished( this ),
  m_chart( nullptr ),
  m_files( nullptr ),
  m_detDiameter( nullptr ),
  m_showFwhmPoints( nullptr ),
  m_fwhmOptionGroup( nullptr ),
  m_fwhmEqnType( nullptr ),
  m_sqrtEqnOrder( nullptr ),
  m_effEqnOrder( nullptr ),
  m_effEqnUnits( nullptr ),
  m_effOptionGroup( nullptr ),
  m_airAttenuate( nullptr ),
  m_chartLowerE( nullptr ),
  m_chartUpperE( nullptr ),
  m_errorMsg( nullptr ),
  m_intrinsicEffAnswer( nullptr ),
  m_fwhmFitId( 0 ),
  m_effEqnFitId( 0 ),
  m_fwhmEqnChi2( -999.9 ),
  m_effEqnChi2( -999.9 ),
  m_effLowerEnergy( 0.0f ),
  m_effUpperEnergy( 0.0f )
{
  assert( m_interspec );
  assert( m_materialDB );
  assert( m_materialSuggest );
  
  wApp->useStyleSheet( "InterSpec_resources/MakeDrf.css" );
  
  addStyleClass( "MakeDrf" );
  
  WGridLayout *upperLayout = new WGridLayout();
  upperLayout->setContentsMargins( 0, 0, 0, 0 );
  upperLayout->setVerticalSpacing( 0 );
  upperLayout->setHorizontalSpacing( 0 );
  
  WContainerWidget *fitOptionsDiv = new WContainerWidget();
  fitOptionsDiv->addStyleClass( "MakeDrfOptions" );
  upperLayout->addWidget( fitOptionsDiv, 0, 0, 2, 1 );
  
  WGroupBox *detDiamGroup = new WGroupBox( "Det. Diameter", fitOptionsDiv );
  m_detDiameter = new WLineEdit( detDiamGroup );
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_detDiameter->setValidator( distValidator );
  m_detDiameter->setText( "2.54 cm" );
  m_detDiameter->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  m_detDiameter->enterPressed().connect( this, &MakeDrf::handleSourcesUpdates );
  
  
  m_effOptionGroup = new WGroupBox( "Intrinsic Eff.", fitOptionsDiv );
  m_effEqnOrder = new WComboBox( m_effOptionGroup );
  m_effEqnOrder->setInline( false );
  m_effEqnOrder->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  
  m_effEqnUnits = new WComboBox( m_effOptionGroup );
  m_effEqnUnits->setInline( false );
  m_effEqnUnits->addItem( "Eqn. in keV" );
  m_effEqnUnits->addItem( "Eqn. in MeV" );
  m_effEqnUnits->setCurrentIndex( 1 );
  m_effEqnUnits->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  
  
  m_fwhmOptionGroup = new WGroupBox( "FWHM Eqn.", fitOptionsDiv );
  m_fwhmEqnType = new WComboBox( m_fwhmOptionGroup );
  m_fwhmEqnType->setInline( false );
  m_fwhmEqnType->addItem( "Gadras Eqation" );
  m_fwhmEqnType->addItem( "Sqrt Power Series" );
  m_fwhmEqnType->setCurrentIndex( 0 );
  m_fwhmEqnType->changed().connect( this, &MakeDrf::handleFwhmTypeChanged );
  
  m_sqrtEqnOrder = new WComboBox( m_fwhmOptionGroup );
  m_sqrtEqnOrder->setInline( false );
  m_sqrtEqnOrder->hide();
  m_sqrtEqnOrder->changed().connect( this, &MakeDrf::handleSqrtEqnOrderChange );
  
  m_effOptionGroup->hide();
  m_fwhmOptionGroup->hide();
  
  WGroupBox *genOpts = new WGroupBox( fitOptionsDiv );
  m_airAttenuate = new WCheckBox( "Atten. for air", genOpts );
  m_airAttenuate->setChecked( true );
  m_airAttenuate->setInline( false );
  m_airAttenuate->addStyleClass( "AirAttenCb" );
  m_airAttenuate->checked().connect( this, &MakeDrf::handleSourcesUpdates );
  m_airAttenuate->unChecked().connect( this, &MakeDrf::handleSourcesUpdates );
  
  m_chart = new MakeDrfChart();
  DrfChartHolder *chartholder = new DrfChartHolder( m_chart, nullptr );
  upperLayout->addWidget( chartholder, 0, 1 );
  
  WContainerWidget *chartOptionsDiv = new WContainerWidget();
  chartOptionsDiv->addStyleClass( "MakeDrfChartOptions" );
  upperLayout->addWidget( chartOptionsDiv, 1, 1 );
  m_showFwhmPoints = new WCheckBox( "Show FWHM points", chartOptionsDiv );
  m_showFwhmPoints->setChecked( true );
  m_showFwhmPoints->checked().connect( this, &MakeDrf::handleShowFwhmPointsToggled );
  m_showFwhmPoints->unChecked().connect( this, &MakeDrf::handleShowFwhmPointsToggled );
  
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
  
  WContainerWidget *fileHolder = new WContainerWidget();
  fileHolder->setOverflow( Overflow::OverflowAuto, Wt::Vertical );
  fileHolder->addStyleClass( "MakeDrfFiles" );
  
  m_files = new WContainerWidget( fileHolder );
  
  m_errorMsg = new WText( "", Wt::XHTMLText );
  m_errorMsg->addStyleClass( "MakeDrfErrTxt" );
  upperLayout->addWidget( m_errorMsg, 2, 0, 1, 2 );
  m_errorMsg->hide();
  
  m_intrinsicEffAnswer = new WText( "", Wt::XHTMLText );
  m_intrinsicEffAnswer->addStyleClass( "MakeDrfEqnAnswer" );
  upperLayout->addWidget( m_intrinsicEffAnswer, 3, 0, 1, 2 );
  
  upperLayout->setRowStretch( 0, 1 );
  upperLayout->setColumnStretch( 1, 1 );
  
  WGridLayout *layout = new WGridLayout();
  setLayout( layout );
  layout->setContentsMargins( 0, 0, 0, 0 );
  
  if( viewer->isPhone() )
  {
    //Phone layout untested!
    WTabWidget *tab = new WTabWidget();
    layout->addWidget( tab, 0, 0 );
    WContainerWidget *chartTab = new WContainerWidget();
    chartTab->setLayout( upperLayout );
    tab->addTab( chartTab, "Chart/Options" );
    tab->addTab( fileHolder, "Peaks To Use" );
  }else
  {
    layout->addLayout( upperLayout, 0, 0 );
    layout->addWidget( fileHolder, 1, 0 );
    layout->setRowResizable( 0, true, 250 );
  }//if( is phone ) / else
  
  SpecMeasManager *manager = viewer->fileManager();
  SpectraFileModel *measmodel = manager->model();
  
  vector<DrfSpecFile *> added;
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
    added.push_back( fileWidget );
    

    for( DrfSpecFileSample *fs : fileWidget->fileSamples() )
    {
      for( DrfPeak *p : fs->peaks() )
        p->m_peakPreviewShow.connect( boost::bind(&MakeDrf::peakPreviewShown, this, _1)  );
    }
    
  }//for( loop over opened files )
  
  if( added.size() == 1 )
  {
    added[0]->addStyleClass( "SingleFile" );
    if( !added[0]->title().empty() )
      added[0]->setTitle( "" );
  }
  
  
  //If we directly call handleSourcesUpdates() now, getting the activity
  //  uncertainty will throw an exception because they wont validate.... whatever.
  WServer::instance()->post( wApp->sessionId(), wApp->bind( boost::bind(&MakeDrf::handleSourcesUpdates,this) ) );
}//MakeDrf( constructor )


MakeDrf::~MakeDrf()
{

}//~MakeDrf()


AuxWindow *MakeDrf::makeDrfWindow( InterSpec *viewer, MaterialDB *materialDB, Wt::WSuggestionPopup *materialSuggest )
{
  AuxWindow *window = new AuxWindow( "Create Detector Response Function",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                     | AuxWindowProperties::SetCloseable
                                     | AuxWindowProperties::DisableCollapse
                                     | AuxWindowProperties::EnableResize) );
  
  const int ww = viewer->renderedWidth();
  const int wh = viewer->renderedHeight();
  if( ww > 100 && wh > 100 )
  {
    const int width = std::min( 3*ww/4, 900 );
    const int height = ((wh < 420) ? wh : (19*wh)/20 );
    
    window->resizeWindow( width, height );
    window->setMinimumSize( std::min(width,640), std::min(height,480) );
  }//if( ww > 100 && wh > 100 )
    
  MakeDrf *makeDrfWidget = new MakeDrf( viewer, materialDB, materialSuggest );
    
  window->stretcher()->addWidget( makeDrfWidget, 0, 0 );
  window->stretcher()->setContentsMargins( 0, 0, 0, 0 );
    
  AuxWindow::addHelpInFooter( window->footer(), "make-drf" );
    
  WPushButton *closeButton = window->addCloseButtonToFooter( "Close" );
  closeButton->clicked().connect( window, &AuxWindow::hide );
    
  WPushButton *saveAs = new WPushButton( "Store/Export...", window->footer() );
  saveAs->clicked().connect( makeDrfWidget, &MakeDrf::startSaveAs );
  makeDrfWidget->intrinsicEfficiencyIsValid().connect( boost::bind( &WPushButton::setEnabled, saveAs, _1 ) );
  saveAs->disable();
    
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
  window->show();
    
  window->resizeToFitOnScreen();
  window->centerWindow();
  window->rejectWhenEscapePressed( false );
  
  makeDrfWidget->m_finished.connect( window, &AuxWindow::hide );
  
  return window;
}//AuxWindow *makeDrfWindow(...)


void MakeDrf::startSaveAs()
{
  Wt::WFlags<AuxWindowProperties> windowprop = Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                                               | AuxWindowProperties::PhoneNotFullScreen
                                               | AuxWindowProperties::SetCloseable
                                               | AuxWindowProperties::DisableCollapse;
  if( m_effEqnCoefs.empty() )
  {
    AuxWindow *w = new AuxWindow( "Error", windowprop );
    WText *t = new WText( "Sorry,&nbsp;DRF&nbsp;not&nbsp;valid.", Wt::XHTMLText, w->contents() );
    t->setInline( false );
    WPushButton *b = new WPushButton( "Close", w->footer() );
    b->clicked().connect( w, &AuxWindow::hide );
    w->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, w) );
    w->rejectWhenEscapePressed();
    w->show();
    w->centerWindow();
    return;
  }//if( m_effEqnCoefs.empty() )
  
  
  //ToDo: - add download link to CSV/TSV absolute eff file (e.g., for S.M.)
  //      - consider creating a log file for user to see all the details; maybe
  //        put this on a popup of the main widget
  //        Could get most info from const m_chart->currentDataPoints();
  //      - add in a interface to select previously characterized detectors, and
  //        also upload them
  
  AuxWindow *w = new AuxWindow( "Store/Export DRF", windowprop );

  WTable *table = new WTable( w->contents() );
  table->addStyleClass( "MakeDrfSaveAsTable" );
  
  WTableCell *cell = table->elementAt(0, 0);
  WLabel *label = new WLabel( "Name", cell );
  
  cell = table->elementAt(0, 1);
  WLineEdit *name = new WLineEdit( cell );
  label->setBuddy( name );
  name->setTextSize( 32 );
  name->setMaxLength( 255 );
  name->setAutoComplete( false );
  //const char *valid_file_name_regex = "(^(?!\\.)(?!com[0-9]$)(?!con$)(?!lpt[0-9]$)(?!nul$)(?!prn$)[^\\|*\\?\\\\:<>/$\"]*[^\\.\\|*\\?\\\\:<>/$\"]+$";
  const char *valid_file_name_regex = R"MyRegexDelim(^[^\\\/\:\*\?\"\'\,\;\<\>\|]+$)MyRegexDelim";
  WRegExpValidator *validator = new WRegExpValidator( valid_file_name_regex, name );
  name->setValidator( validator );
  
  cell = table->elementAt(0, 2);
  WImage *help = new WImage(Wt::WLink("InterSpec_resources/images/help_mobile.svg"), cell);
  help->addStyleClass( "MakeDrfSaveHelp" );
  const char *tooltip = "Name of the detector response function to save to <em>InterSpecs</em>"
                        " internal database so you can use the response function later."
                        "  Name must be a valid filename (e.g., no \\\\, /, :, etc. characters).";
  HelpSystem::attachToolTipOn( help, tooltip, true, HelpSystem::ToolTipPosition::Left );
  
  cell = table->elementAt(1, 0);
  label = new WLabel( "Description", cell );
  cell = table->elementAt(1, 1);
  cell->setRowSpan( 2 );
  WTextArea *description = new WTextArea( cell );
  label->setBuddy( description );
  description->setWidth( WLength(97,WLength::Percentage) );
  
  cell = table->elementAt(1, 2);
  help = new WImage(Wt::WLink("InterSpec_resources/images/help_mobile.svg"), cell);
  help->addStyleClass( "MakeDrfSaveHelp" );
  tooltip = "Optional description of the DRF to help remind yourself of details when you use the DRF in the future.";
  HelpSystem::attachToolTipOn( help, tooltip, true, HelpSystem::ToolTipPosition::Left );
  
  std::shared_ptr<const SpecMeas> representative_meas;
  
  for( auto w : m_files->children() )
  {
    auto f = dynamic_cast<DrfSpecFile *>( w );
    if( !f )
      continue;
    size_t npeaks = 0;
    for( auto sample : f->fileSamples() )
    {
      for( auto peak : sample->peaks() )
        npeaks += peak->useForEffFit();
    }
    
    if( npeaks )
    {
      representative_meas = f->measurement();
      break;
    }
  }//for( auto w : m_files->children() )

  WCheckBox *def_for_serial_cb = nullptr;
  WCheckBox *def_for_model_cb = nullptr;
  
  int currentRow = table->rowCount();
  if( representative_meas )
  {
    string serial_number = representative_meas->instrument_id();
    
    if( !serial_number.empty() )
    {
      cell = table->elementAt(currentRow, 0);
      cell->setColumnSpan( 2 );
    
      def_for_serial_cb = new WCheckBox( "Make default DRF for serial number '" + serial_number + "'", cell );
      cell = table->elementAt(currentRow, 2);
      help = new WImage(Wt::WLink("InterSpec_resources/images/help_mobile.svg"), cell);
      help->addStyleClass( "MakeDrfSaveHelp" );
      tooltip = "Make it so when spectra from a detector with a matching serial"
                " number is loaded into InterSpec, this detector response function"
                " will automatically be loaded and used.";
      HelpSystem::attachToolTipOn( help, tooltip, true, HelpSystem::ToolTipPosition::Left );
    }//if( serial_number.size() )
    
    string model;
    if( representative_meas->detector_type() == SpecUtils::DetectorType::Unknown )
      model = representative_meas->instrument_model();
    else
      model = detectorTypeToString( representative_meas->detector_type() );
    
    if( !model.empty() )
    {
      currentRow = table->rowCount();
      cell = table->elementAt(currentRow, 0);
      cell->setColumnSpan( 2 );
      def_for_model_cb = new WCheckBox( "Make default DRF for model '" + model + "'", cell );
      cell = table->elementAt(currentRow, 2);
      help = new WImage(Wt::WLink("InterSpec_resources/images/help_mobile.svg"), cell);
      help->addStyleClass( "MakeDrfSaveHelp" );
      tooltip = "Make it so when spectra from this model detector is loaded"
                " into InterSpec, this detector response function will"
                " automatically be loaded and used.";
      HelpSystem::attachToolTipOn( help, tooltip, true, HelpSystem::ToolTipPosition::Left );
    }//if( !model.empty() )
  }//if( representative_meas )
  
  
  currentRow = table->rowCount();
  cell = table->elementAt(currentRow, 0);
  cell->setColumnSpan( 2 );
  CalFileDownloadResource *n42Resource = new CalFileDownloadResource( false, this );
  WAnchor *n42anchor = new WAnchor( n42Resource, "Export data as N42-2012 file.", cell );
  n42anchor->setTarget( AnchorTarget::TargetNewWindow );
  
  cell = table->elementAt(currentRow, 2);
  help = new WImage(Wt::WLink("InterSpec_resources/images/help_mobile.svg"), cell);
  help->addStyleClass( "MakeDrfSaveHelp" );
  tooltip = "Export a N42 file containing all spectra and peaks used to create this DRF."
            "  Source and shielding information is also usually stored into the file as well to,"
            " in principal, except for detector diameter and equation orders,"
            " provide all the input information used to create the DRF.";
  HelpSystem::attachToolTipOn( help, tooltip, true, HelpSystem::ToolTipPosition::Left );
  
  currentRow = table->rowCount();
  cell = table->elementAt(currentRow, 0);
  cell->setColumnSpan( 2 );
  CsvDrfDownloadResource *csvResource = new CsvDrfDownloadResource( this );
  WAnchor *csvanchor = new WAnchor( csvResource, "Export DRF as CSV.", cell );
  csvanchor->setTarget( AnchorTarget::TargetNewWindow );
  

  cell = table->elementAt(currentRow, 2);
  help = new WImage(Wt::WLink("InterSpec_resources/images/help_mobile.svg"), cell);
  help->addStyleClass( "MakeDrfSaveHelp" );
  tooltip = "Exports the DRF into a CSV file that contains all of the information of the DRF."
  " Especially useful for using with other tools.";
  HelpSystem::attachToolTipOn( help, tooltip, true, HelpSystem::ToolTipPosition::Left );
  
  
  currentRow = table->rowCount();
  cell = table->elementAt(currentRow, 0);
  cell->setColumnSpan( 2 );
  RefSheetDownloadResource *refSheetResource = new RefSheetDownloadResource( this );
  WAnchor *refSheetAnchor = new WAnchor( refSheetResource, "Export quick ref card", cell );
  refSheetAnchor->setTarget( AnchorTarget::TargetNewWindow );
  
  cell = table->elementAt(currentRow, 2);
  help = new WImage(Wt::WLink("InterSpec_resources/images/help_mobile.svg"), cell);
  help->addStyleClass( "MakeDrfSaveHelp" );
  tooltip = "Exports a 3x5 style card that has a brief summary of detector performance.";
  HelpSystem::attachToolTipOn( help, tooltip, true, HelpSystem::ToolTipPosition::Left );
  
  
  auto updateName = [name,csvResource,n42Resource,refSheetResource](){
    if( name->validate() == Wt::WValidator::Valid )
    {
      n42Resource->setSuggestFileName( name->text().toUTF8() );
      csvResource->setSuggestFileName( name->text().toUTF8() );
      refSheetResource->setSuggestFileName( name->text().toUTF8() );
    }else
    {
      n42Resource->setSuggestFileName( "" );
      csvResource->setSuggestFileName( "" );
      refSheetResource->setSuggestFileName( "" );
    }
  };
  
  name->changed().connect( std::bind(updateName) );
  name->enterPressed().connect( std::bind(updateName) );
  name->blurred().connect( std::bind(updateName) );
  
  
  auto updateDesc = [description,csvResource,refSheetResource](){
    csvResource->setDescription( description->text().toUTF8() );
    refSheetResource->setDescription( description->text().toUTF8() );
  };
  
  description->changed().connect( std::bind(updateDesc) );
  description->enterPressed().connect( std::bind(updateDesc) );
  description->blurred().connect( std::bind(updateDesc) );

  
  
  WPushButton *b = w->addCloseButtonToFooter( "Cancel" );
  b->clicked().connect( w, &AuxWindow::hide );
  
  // TODO: Need to display messsage to user about how it is saved, how to access it again, and figure out why menus go away for people (maybe only on Windows?) - also should consider not closing AuxWindow - also fix assay date being invalid error (seems to cause uncaught exception).
  auto doSave = [this, w, validator, name, description,
                 def_for_serial_cb, def_for_model_cb,
                 representative_meas ](){
    auto state = validator->validate(name->text()).state();
    if( name->text().empty() || state!=WValidator::Valid )
    {
      passMessage( "Detector name not valid.", "", 3 );
      return;
    }
    
    if( m_effEqnCoefs.empty() )
    {
      passMessage( "Intrinsic Efficiency Equation is not valid.", "", 3 );
      return;
    }
    
    string drfname = name->text().toUTF8();
    string drfdescrip = description->text().toUTF8();
    SpecUtils::trim( drfname );
    SpecUtils::trim( drfdescrip );
    
    auto drf = make_shared<DetectorPeakResponse>( drfname, drfdescrip );
    
    float diameter;
    try
    {
      const double diam = PhysicalUnits::stringToDistance( m_detDiameter->text().toUTF8() );
      diameter = static_cast<float>( diam );
    }catch(...)
    {
      passMessage( "Detector diameter entered is not a valid distance.", "", 3 );
      return;
    }
    
    const bool inMeV = isEffEqnInMeV();
    const float eqnEnergyUnits = inMeV ? 1000.0f : 1.0f;
    
    float lowerEnergy = 0.0f, upperEnergy = 0.0f;
    const std::vector<MakeDrfChart::DataPoint> &data = m_chart->currentDataPoints();
    if( data.size() >= 2 )
    {
      lowerEnergy = data.front().energy;
      upperEnergy = data.back().energy;
    }
    
    drf->fromExpOfLogPowerSeriesAbsEff( m_effEqnCoefs, m_effEqnCoefUncerts,
                                        0.0f, diameter, eqnEnergyUnits, lowerEnergy, upperEnergy );
    drf->setDrfSource( DetectorPeakResponse::DrfSource::UserCreatedDrf );
    
    if( !m_fwhmCoefs.empty() )
    {
      switch( m_fwhmEqnType->currentIndex() )
      {
        case 0:
          drf->setFwhmCoefficients( m_fwhmCoefs, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn );
          break;
          
        case 1:
          drf->setFwhmCoefficients( m_fwhmCoefs, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
          break;
      }//switch( m_fwhmEqnType->currentIndex() )
    }//if( !m_fwhmCoefs.empty() )
  
    //Need something here to indicate this is a created DRF.
    
    try
    {
      auto sql = std::make_shared<DataBaseUtils::DbSession>( *m_interspec->sql() );
      
      //negladgeable chance we have another detector with same hash, so we'll
      // skip checking for it.
      
      DataBaseUtils::DbTransaction transaction( *sql );
      //Create a separate DetectorPeakResponse because shared_ptr and dbo::ptr don't work well together
      DetectorPeakResponse *tempDetector = new DetectorPeakResponse( *drf );
      tempDetector->m_user = m_interspec->m_user.id();
      auto newDbDet = sql->session()->add( tempDetector );
      
      transaction.commit();
    }catch( std::exception &e )
    {
      passMessage( "Error saving DetectorPeakResponse to data base: " + std::string(e.what()), "", WarningWidget::WarningMsgHigh );
      return;
    }//try / catch
    
    m_interspec->detectorChanged().emit( drf );
    
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    Wt::Dbo::ptr<InterSpecUser> user = m_interspec->m_user;
    
    if( def_for_serial_cb && def_for_serial_cb->isChecked() && representative_meas )
    {
      UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorSerialNumber;
      WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
        DrfSelect::setUserPrefferedDetector( drf, sql, user, preftype, representative_meas );
      } ) );
    }//if( def_for_serial_cb and is checked )
    
    if( def_for_model_cb && def_for_model_cb->isChecked() && representative_meas )
    {
      UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorModelName;
      WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
        DrfSelect::setUserPrefferedDetector( drf, sql, user, preftype, representative_meas );
      } ) );
    }//if( def_for_serial_cb and is checked )
    
    w->hide();
    
    //m_finished.emit();
    passMessage( "Saved '" + drfname + "' to the internal database.", "", 0 );
  };//auto doSave
  
  WPushButton *save = w->addCloseButtonToFooter( "Save" );
  save->clicked().connect( std::bind(doSave) );
  save->disable();
  
  auto updateSaveBtn = [name,save,validator](){
    bool enable = false;
    switch( validator->validate(name->text()).state() )
    {
      case WValidator::Invalid:      enable = false; break;
      case WValidator::InvalidEmpty: enable = false; break;
      case WValidator::Valid:        enable = !name->text().empty(); break;
    }
    save->setEnabled( enable );
  };
  
  name->validated().connect( std::bind(updateSaveBtn) );  //doesnt seem to ever get called...
  //name->changed().connect( std::bind(updateSaveBtn) );
  //name->enterPressed().connect( std::bind(updateSaveBtn) );
  name->keyPressed().connect( std::bind(updateSaveBtn) ); //Uhg, I guess we'll just use this until I figure out whats wrong with above so save button will update in real time
  
  
  w->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, w) );
  w->rejectWhenEscapePressed();
  w->show();
  w->centerWindow();
}//void startSaveAs();


Wt::Signal<bool> &MakeDrf::intrinsicEfficiencyIsValid()
{
  return m_intrinsicEfficiencyIsValid;
}


bool MakeDrf::isIntrinsicEfficiencyValid()
{
  return !m_effEqnCoefs.empty();
}//bool MakeDrf::isDrfValid()


void MakeDrf::handleSourcesUpdates()
{
  size_t numchan = 0;
  vector< std::shared_ptr<const PeakDef> > peaks;
  vector<MakeDrfFit::DetEffDataPoint> effpoints;
  
  bool highres = false;
  
  bool detDiamInvalid = false;
  double diameter = 2.54*PhysicalUnits::cm;
  try
  {
    diameter = PhysicalUnits::stringToDistance( m_detDiameter->text().toUTF8() );
  }catch(...)
  {
    detDiamInvalid = true;
  }
  
  const bool doAirAtten = m_airAttenuate->isChecked();
  
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
        if( peak->useForEffFit() )
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
            highres = (highres || PeakFitUtils::is_high_res(m));
            
            minenergy = std::min( m->gamma_energy_min(), minenergy );
            maxenergy = std::max( m->gamma_energy_max(), maxenergy );
          }
        }
      }//for( auto det : meas->detector_names() )
      

      const vector<pair<DrfPeak *,MakeDrfSrcDef *>> peaks_to_sources = sample->selected_peak_to_sources();
      map<const SandiaDecay::Nuclide *,vector<SandiaDecay::EnergyRatePair>> mixtures;
      
      // ageAtSpectrumTime() or activityAtSpectrumTime() can through an exception if the entered
      //  values are invalid; we will skip the entire sample if we encounter this.
      try
      {
        for( auto pp : peaks_to_sources )
        {
          const SandiaDecay::Nuclide * const nuc = pp.second->nuclide();
          if( !pp.first || !nuc || mixtures.count(nuc) )
            continue;
          
          const double activity = pp.second->activityAtSpectrumTime();
          const double age = pp.second->ageAtSpectrumTime();
          
          SandiaDecay::NuclideMixture mix;
          mix.addAgedNuclideByActivity( nuc, activity, age );
          
          mixtures[nuc] = mix.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
        }//for( auto pp : peaks_to_sources )
      }catch( std::exception & )
      {
        not_all_sources_used = true;
        continue;
      }//try / catch
      
      for( auto pp : peaks_to_sources )
      {
        MakeDrfSrcDef * const srcDef = pp.second;
        const SandiaDecay::Nuclide * const nuc = srcDef->nuclide();
        
        DrfPeak * const drfPeak = pp.first;
        const std::shared_ptr<const PeakDef> peak = drfPeak->m_peak;
        
        //Check if we will background subtract
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
        drfPeak->setBackgroundBeingSubtractedInfo( subBack, (subBack ? (back_peak_area/back_peak_lt) : 0.0) );
        
        MakeDrfChart::DataPoint point;
        point.energy = peak->parentNuclide() ? peak->gammaParticleEnergy() : peak->mean();
        point.livetime = drfPeak->m_livetime;
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
        
        
        
        auto trans_frac = [doAirAtten]( const float energy, MakeDrfSrcDef *src ) -> double {
          
          double airTransFrac = 1.0;
          if( doAirAtten )
          {
            const double distance = src->distance();
            const double mu = GammaInteractionCalc::transmission_length_coefficient_air( energy );
            airTransFrac = exp( -mu * distance );
          }//if( doAirAtten )
          
          ShieldingSelect *shield = src ? src->shielding() : nullptr;
          if( !shield )
            return airTransFrac;
          
          //ToDo: Attenuation calculation not checked!
          double an = 14, ad = 0.0;
          if( shield->isGenericMaterial() )
          {
            an = shield->atomicNumber();
            ad = shield->arealDensity();
            
            if( an < 1.0 || ad < 0.0 )
            {
              an = 14;
              ad = 0.0;
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
          
          const double mu = MassAttenuation::massAttenuationCoeficient( an, energy );
          return airTransFrac * exp( -mu * ad );
        };//trans_frac lambda
        
        
        auto src_rate = [&mixtures,trans_frac,diameter]( const float energy, const float width,
                                                MakeDrfSrcDef * const src, DrfPeak * const peak,
                                                float &source_count_rate,
                                                float &incidentGamma, float &fracUncert ) {
          source_count_rate = incidentGamma = fracUncert = 0.0;
          
          const double distance = src->distance();
          const double transmittion_factor = trans_frac( energy, src );
          const SandiaDecay::Nuclide * const nuc = src->nuclide();
          
          if( nuc )
          {
            const vector<SandiaDecay::EnergyRatePair> &rates = mixtures[nuc];
            for( const auto &r : rates )
            {
              if( fabs(r.energy - energy) < width )
                source_count_rate += (r.numPerSecond * transmittion_factor);
            }
          }else
          {
            if( !peak->m_userBr || (peak->m_userBr->validate()!=WValidator::State::Valid) )
              throw runtime_error( "Logic error - no nuclide associated with " + std::to_string(peak->m_peak->mean()) + " keV peak, and no user BR." );
            
            const double activity = src->activityAtSpectrumTime();
            const double br = peak->m_userBr->value();
            source_count_rate += transmittion_factor * br * activity;
          }
          
          //If we get here on initial GUI load, we seem to get error validating
          //  the uncertainty (hence a workaround is used to delay first reading)
          
          if( source_count_rate > 0.0f )
          {
            const double livetime = peak->m_livetime;
            const double fracSolidAngle = DetectorPeakResponse::fractionalSolidAngle( diameter, distance );
            incidentGamma = source_count_rate * livetime * fracSolidAngle;
            fracUncert = src->fractionalActivityUncertainty();
          }//if( source_count_rate > 0.0f )
        };//src_rate labmda
        
        MakeDrfFit::DetEffDataPoint effpoint;
        const double width = 1.25*(peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth());
        
        try
        {
          point.distance = srcDef->distance();
          
          float source_count_rate, expected, fracUncert;
          src_rate( point.energy, width, srcDef, drfPeak, source_count_rate, expected, fracUncert );
          
          
          //Now loop through all the other nuclides and get their expected contribtion
          for( const auto &nuc_to_rates : mixtures )
          {
            const SandiaDecay::Nuclide * const other_nuc = nuc_to_rates.first;
            if( other_nuc == nuc || !other_nuc )
              continue;
            
            MakeDrfSrcDef *otherSrcDef = nullptr;
            for( size_t i = 0; !otherSrcDef && i < peaks_to_sources.size(); ++i )
            {
              if( peaks_to_sources[i].second->nuclide() == other_nuc )
                otherSrcDef = peaks_to_sources[i].second;
            }
            
            if( !otherSrcDef )
              continue;
            
            assert( otherSrcDef != srcDef );
            
            float other_source_count_rate, other_expected, other_fracUncert;
            src_rate( point.energy, width, otherSrcDef, drfPeak, other_source_count_rate, other_expected, other_fracUncert );
            
            if( other_source_count_rate > 0.0 )
            {
              fracUncert = sqrt( pow(expected*fracUncert,2.0f) + pow(other_expected*other_fracUncert,2.0f) );
              expected += other_expected;
              source_count_rate += other_source_count_rate;
              
              //ToDo: To MakeDrfChart::DataPoint add a other source count rate and uncertainty field.
              
              fracUncert /= expected;
            }
          }//for( const auto &nuc_to_rates : mixtures )
          
          
          point.source_count_rate = source_count_rate;
          point.source_count_rate_uncertainty = source_count_rate*fracUncert;
          
          const double eff = point.peak_area / expected;
          
          double fracUncert2 = 0.0;
          if( point.peak_area_uncertainty > 0.0f )
            fracUncert2 += std::pow( point.peak_area_uncertainty / point.peak_area, 2.0f );
          if( point.source_count_rate_uncertainty > 0.0f )
            fracUncert2 += std::pow( point.source_count_rate_uncertainty / point.source_count_rate, 2.0f );
          
          effpoint.energy = point.energy;
          effpoint.efficiency = eff;
          effpoint.efficiency_uncert = sqrt(fracUncert2);
        }catch( std::exception &e )
        {
          cerr << "handleSourcesUpdates: got exception: " << e.what() << endl;
          not_all_sources_used = true;
          continue;
        }
        
        point.peak_color = peak->lineColor();
        point.background_peak_area = back_peak_area;
        point.background_peak_live_time = back_peak_lt;
        
        //ToDo: Add more/better source info.
        char buffer[256] = { '\0' };
        if( nuc )
        {
          snprintf( buffer, sizeof(buffer)-1, "%s %.1f keV, Peak %.1f counts",
                   nuc->symbol.c_str(), point.energy, point.peak_area );
        }else
        {
          snprintf( buffer, sizeof(buffer)-1, "No assoc. nuclide, %.1f keV, Peak %.1f counts",
                    point.energy, point.peak_area );
        }
        
        point.source_information = buffer;
        
        
        effpoints.push_back( effpoint );
        peaks.push_back( peak );
        datapoints.push_back( point );
      }//for( auto pp : peaks_to_sources )
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
  
  const int currentSqrtOrder = m_sqrtEqnOrder->currentIndex();
  const int nCurrentSqrtOrders = m_sqrtEqnOrder->count();
  if( nCurrentSqrtOrders==6 && numPeaks>=6 )
  {
    //Dont need to do anything
  }else if( nCurrentSqrtOrders != numPeaks )
  {
    m_sqrtEqnOrder->clear();
    for( int order = 0; order < numPeaks && order < 6; ++order )
      m_sqrtEqnOrder->addItem( std::to_string(order+1) + " Fit Params"  );
    
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
  
  m_effOptionGroup->setHidden( datapoints.empty() );
  m_fwhmOptionGroup->setHidden( datapoints.empty() );
  
  
  fitEffEqn( effpoints );
  fitFwhmEqn( peaks, highres );
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


void MakeDrf::peakPreviewShown( DrfPeak *peak )
{
  //Hide all other previews showing.
  for( auto w : m_files->children() )
  {
    auto fileWidget = dynamic_cast<DrfSpecFile *>( w );
    if( !fileWidget )
      continue;
    
    for( DrfSpecFileSample *fs : fileWidget->fileSamples() )
    {
      for( DrfPeak *p : fs->peaks() )
      {
        if( p != peak )
          p->hidePeakPreview();
      }//for( DrfPeak *p : fs->peaks() )
    }//for( DrfSpecFileSample *fs : fileWidget->fileSamples() )
  }//for( auto w : m_files->children() )
}//void peakPreviewShown( DrfPeak *peak )


void MakeDrf::fitFwhmEqn( std::vector< std::shared_ptr<const PeakDef> > peaks,
                          const bool isHighResolution )
{
  ++m_fwhmFitId;
  
  const int fitid = m_fwhmFitId;
  const string sessionId = wApp->sessionId();
  
  m_fwhmEqnChi2 = -999.9;
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
  auto updater = boost::bind( &MakeDrf::updateFwhmEqn, this, _1, _2, _3, static_cast<int>(fnctnlForm), fitid );
  
  const string thisid = id();
  
  auto worker = [sessionId,fnctnlForm,peaks,isHighResolution,sqrtEqnOrder,updater,thisid]() {
    try
    {
      auto peakdequ = std::make_shared<std::deque< std::shared_ptr<const PeakDef> > >( peaks.begin(), peaks.end() );
    
      //Takes between 5 and 500ms for a HPGe detector
      const double start_time = SpecUtils::get_wall_time();
      vector<float> fwhm_coefs, fwhm_coefs_uncert;
      const double chi2 = MakeDrfFit::performResolutionFit( peakdequ, fnctnlForm, isHighResolution, sqrtEqnOrder, fwhm_coefs, fwhm_coefs_uncert );
    
      const double end_time = SpecUtils::get_wall_time();
    
      assert( fwhm_coefs.size() == fwhm_coefs_uncert.size() );
      cout << "Fit FWHM: {";
      for( size_t i = 0; i < fwhm_coefs.size(); ++i )
        cout << fwhm_coefs[i] << "+-" << fwhm_coefs_uncert[i] << ", ";
      cout << "}; took " << (end_time-start_time) << " seconds" << endl;
      
      WServer::instance()->post( sessionId, std::bind( [updater,fwhm_coefs,fwhm_coefs_uncert,thisid,chi2](){
        if( wApp->domRoot() && dynamic_cast<MakeDrf *>(wApp->domRoot()->findById(thisid)) )
          updater(fwhm_coefs,fwhm_coefs_uncert,chi2);
        else
          cerr << "MakeDrf widget was deleted while calculating FWHM coefs" << endl;
      } ) );
    }catch( std::exception &e )
    {
      cout << "Failed to fit FWHM coefs: " << e.what() << endl;
    }//try / catch fit FWHM
  };
  
  WServer::instance()->ioService().boost::asio::io_service::post( worker );
}//void fitFwhmEqn( std::vector< std::shared_ptr<const PeakDef> > peaks )


void MakeDrf::updateFwhmEqn( std::vector<float> coefs,
                             std::vector<float> uncerts,
                             const double chi2,
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
  
  m_fwhmEqnChi2 = chi2;
  m_fwhmCoefs = coefs;
  m_fwhmCoefUncerts = uncerts;
  m_chart->setFwhmCoefficients( coefs, uncerts, eqnType, MakeDrfChart::EqnEnergyUnits::keV );
  
  wApp->triggerUpdate();
}//void updateFwhmEqn(...)


void MakeDrf::fitEffEqn( std::vector<MakeDrfFit::DetEffDataPoint> data )
{
  ++m_effEqnFitId;
  
  const int fitid = m_effEqnFitId;
  const string sessionId = wApp->sessionId();
  
  m_effEqnChi2 = -999.9;
  m_effEqnCoefs.clear();
  m_effEqnCoefUncerts.clear();
  m_intrinsicEfficiencyIsValid.emit( !m_effEqnCoefs.empty() );
  m_chart->setEfficiencyCoefficients( m_effEqnCoefs, m_effEqnCoefUncerts, MakeDrfChart::EqnEnergyUnits::keV );
  m_intrinsicEffAnswer->setText( "" );
  
  const int eqnOrderIndex = m_effEqnOrder->currentIndex();
  if( data.empty() )
    return;
  
  if( eqnOrderIndex < 0 )
    return;
  
  const int nfitpars = std::min( eqnOrderIndex+1, 8 );
  
  
  const bool inMeV = isEffEqnInMeV();
  if( inMeV )
  {
    for( auto &a : data )
      a.energy /= 1000.0f;
  }//if( inMeV )
  
  //ToDo: I'm not entirely sure the next line protects against updateEffEqn()
  //  not being called if this widget is deleted before fit is done. (it doesnt!)
  auto updater = boost::bind( &MakeDrf::updateEffEqn, this, _1, _2, _3, _4, _5, fitid, _6 );
  const string thisid = id();
  
  
  auto worker = [sessionId,thisid,data,nfitpars,updater,inMeV]() {
    try
    {
      //Takes between 5 and 500ms for a HPGe detector
      const double start_time = SpecUtils::get_wall_time();
      vector<float> result, uncerts;
      const double chi2 = MakeDrfFit::performEfficiencyFit( data, nfitpars, result, uncerts );
      
      const double end_time = SpecUtils::get_wall_time();
      
      assert( result.size() == uncerts.size() );
      cout << "Fit Eff: {";
      for( size_t i = 0; i < result.size(); ++i )
        cout << result[i] << "+-" << uncerts[i] << ", ";
      cout << "}; took " << (end_time-start_time) << " seconds" << endl;
      
      float lowestEnergy = 99999.0f, highestEnergy = -99999.0f;
      for( const MakeDrfFit::DetEffDataPoint &p : data )
      {
        lowestEnergy = std::min( lowestEnergy, (inMeV ? 1000.0f : 1.0f) * p.energy );
        highestEnergy = std::max( highestEnergy, (inMeV ? 1000.0f : 1.0f) * p.energy );
      }
      
      WServer::instance()->post( sessionId, std::bind( [updater,thisid,result,uncerts,chi2,lowestEnergy,highestEnergy](){
        //Make sure *this is still in the widget tree (incase user closed window while computation was being done)
        if( wApp->domRoot() && dynamic_cast<MakeDrf *>(wApp->domRoot()->findById(thisid) ) )
          updater( result, uncerts, chi2, lowestEnergy, highestEnergy, string("") );
        else
          cerr << "MakeDrf widget was deleted while efficiency was being calculated" << endl;
      } ) );
    }catch( std::exception &e )
    {
      const string errmsg = e.what();
      cout << "Failed to fit intrinsic eff coefs: " << errmsg << endl;
      WServer::instance()->post( sessionId, std::bind( [updater,errmsg,thisid](){
        if( wApp->domRoot() && dynamic_cast<MakeDrf *>(wApp->domRoot()->findById(thisid) ) )
          updater( vector<float>(), vector<float>(), -999.9, 0.0f, 0.0f, errmsg );
        else
          cerr << "MakeDrf widget was deleted while efficiency was being calculated" << endl;
      } ) );
    }//try / catch fit FWHM
  };
  
  WServer::instance()->ioService().boost::asio::io_service::post( worker );
}//void fitEffEqn( std::vector<MakeDrfFit::DetEffDataPoint> data )


void MakeDrf::updateEffEqn( std::vector<float> coefs, std::vector<float> uncerts,
                            const double chi2,
                            const float lowestEnergy, const float highestEnergy,
                            const int fitid, const string errmsg )
{
  const bool isMeV = isEffEqnInMeV();
  const auto units = (isMeV ? MakeDrfChart::EqnEnergyUnits::MeV : MakeDrfChart::EqnEnergyUnits::keV);

  m_effEqnChi2 = chi2;
  m_effLowerEnergy = lowestEnergy;
  m_effUpperEnergy = highestEnergy;
  m_effEqnCoefs = coefs;
  m_effEqnCoefUncerts = uncerts;
  m_intrinsicEfficiencyIsValid.emit( !m_effEqnCoefs.empty() );
  m_chart->setEfficiencyCoefficients( coefs, uncerts, units );
  
  if( !errmsg.empty() )
  {
    m_intrinsicEffAnswer->setText( "" );
    //m_intrinsicEffAnswer->setHidden( true );
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
  
  if( !coefs.empty() )
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
    //m_intrinsicEffAnswer->setHidden( false );
  }else
  {
    m_intrinsicEffAnswer->setText( "" );
    //m_intrinsicEffAnswer->setHidden( true );
  }//if( !coefs.empty() ) / else
  
  wApp->triggerUpdate();
}//void updateEffEqn(...)


std::shared_ptr<SpecMeas> MakeDrf::assembleCalFile()
{
  auto answer = std::make_shared<SpecMeas>();
  
  try
  {
    //Go through and and grab background peaks
    map<std::shared_ptr<SpecUtils::Measurement>,vector<shared_ptr<PeakDef>>> meas_to_peaks;
    
    for( auto w : m_files->children() )
    {
      auto fileWidget = dynamic_cast<DrfSpecFile *>( w );
      if( !fileWidget )
        continue;
      
      for( DrfSpecFileSample *sample : fileWidget->fileSamples() )
      {
        auto meas = sample->measurement();
        set<int> samplenums = sample->samples();
        if( !meas )
          continue;  //shouldnt happen
        
        vector<shared_ptr<PeakDef>> newpeaks;
        
        bool beingUsed = false;
        for( auto peak : sample->peaks() )
        {
          shared_ptr<PeakDef> newpeak = make_shared<PeakDef>( *peak->m_peak );
#ifndef _WIN32
#warning "newpeak->setUseForDrfIntrinsicEffFit not needed anymore"
#endif
          newpeak->setUseForDrfIntrinsicEffFit( peak->useForEffFit() );
          //blah blah blah I dont think this will be needed after we add UseForDrfIntrinsicEffFit as a peak model, assuming m_peak is same as the peak model has...
          newpeaks.push_back( newpeak );
          beingUsed = (beingUsed || peak->useForEffFit());
        }
        
        if( !beingUsed )
          continue;
        
        auto detnames = meas->detector_names();
        std::shared_ptr<SpecUtils::Measurement> m;
        if( detnames.size() == 1 && samplenums.size()==1 )
        {
          auto origm = meas->measurement( *samplenums.begin(), detnames[0] );
          if( origm )
            m = make_shared<SpecUtils::Measurement>( *origm );
        }
        
        if( !m )
          m = meas->sum_measurements( samplenums, detnames, nullptr );
        if( !m )
          continue;  //shouldnt happen
        
        answer->add_measurement( m, false );
        
        const auto srctype = sample->isBackground() ? SpecUtils::SourceType::Background : SpecUtils::SourceType::Foreground;
        answer->set_source_type( srctype, m );
        
        //ToDo: look through remarks for something like "Source:" and remove it.
        //ToDo: Add a remark starting with "Source:" and put source information in (or properly deal with storing source info).
        
        meas_to_peaks[m] = newpeaks;
        
        //These next values will just end up with whatever the last measurement
        //  used had; usually all files used will have the same info.
        answer->set_manufacturer( meas->manufacturer() );
        answer->set_detector_type( meas->detector_type() );
        answer->set_instrument_id( meas->instrument_id() );
        answer->set_instrument_type( meas->instrument_type() );
        answer->set_instrument_model( meas->instrument_model() );
        
        double distance = -1.0;
        vector<string> newremarks;
        for( MakeDrfSrcDef *source : sample->sources() )
        {
          try
          {
            newremarks.push_back( "Source: " + source->toGadrasLikeSourceString() );
          }catch( std::exception &e )
          {
            //ToDo: Do we need to handle the error better here?
            cerr << "Caught exception calling MakeDrfSrcDef::toGadrasLikeSourceString(): " << e.what() << endl;
          }
          try
          {
            distance = source->distance();
          }catch(...){}
        }//
        
        for( const string &remark : m->remarks() )
        {
          if( !SpecUtils::icontains(remark, "Source:") )
            newremarks.push_back( remark );
        }//for( MakeDrfSrcDef *source : sample->sources() )
        
        answer->set_remarks( newremarks, m );
        
        //If we have a distance, put it in the title
        if( distance > 0.0 )
        {
          string title = m->title();
          std::smatch mtch;
          
          string regexstr = PhysicalUnits::sm_distanceRegex;
          if( regexstr.size()>1 && regexstr[0]=='^')
            regexstr = regexstr.substr(1);
          if( regexstr.size()>1 && regexstr[regexstr.size()-1]=='$')
            regexstr = regexstr.substr(0,regexstr.size()-1);
          regexstr = string(".+(@\\s*") + regexstr + ").*?";
          
          std::regex expr( regexstr, std::regex::icase );
          
          if( std::regex_match( title, mtch, expr ) )
          {
            size_t distpos = title.find( mtch[1].str() );
            title = SpecUtils::trim_copy( title.substr(0,distpos) )
                    + " @ " + PhysicalUnits::printToBestLengthUnits(distance)
                    + " " + SpecUtils::trim_copy( title.substr(distpos+mtch[1].str().size()) );
            //the last space we added is sometimes a waste, and results in something like 'U232 @ 33.33 cm , H=100 cm'
          }else if( title.find('@') == string::npos )
          {
            //maybe the regex isnt complete; we'll just skip distance if there is a '@' symbol
            //  anywhere in the title
            
            title = SpecUtils::trim_copy(title)
                    + " @ " + PhysicalUnits::printToBestLengthUnits(distance);
          }
          
          SpecUtils::trim( title );
          answer->set_title( title, m );
        }//if( distance > 0.0 )
        
      }//for( DrfSpecFileSample *sample : sampleWidgets )
    }//for( auto w : m_files->children()  )
    
    answer->cleanup_after_load( 0x0 /*SpecUtils::SpecFile::CleanupAfterLoadFlags::DontChangeOrReorderSamples*/ );
    
    auto newmeasvec = answer->measurements();
    
    for( const auto &measpeak : meas_to_peaks )
    {
      set<int> samples;
      samples.insert( measpeak.first->sample_number() );
      deque<shared_ptr<const PeakDef>> newpeakdeque( begin(measpeak.second), end(measpeak.second) );
      answer->setPeaks( newpeakdeque, samples );
      cout << "Set " << newpeakdeque.size() << " peaks for sample " << measpeak.first->sample_number() << " in new SpecMeas" << endl;
    }
  }catch( std::exception &e )
  {
    cerr << "Caught excpetion in assembleCalFile(): " << e.what();
    return nullptr;
  }//try / catch
  
  if( answer->measurements().empty() )
    return nullptr;
  
  return answer;
}//std::shared_ptr<SpecMeas> assembleCalFile()


void MakeDrf::writeCsvSummary( std::ostream &out,
                               std::string drfname, std::string drfdescription )
{
  const char * const endline = "\r\n";
  
  //ToDo: implement proper reading in of escaped fields, including new lines in
  //      DrfSelect, then get rid of the below.
  SpecUtils::ireplace_all( drfname, "\r", " ");
  SpecUtils::ireplace_all( drfname, "\n", " ");
  SpecUtils::ireplace_all( drfdescription, "\r", " ");
  SpecUtils::ireplace_all( drfdescription, "\n", " ");
  
  const double effChi2 = m_effEqnChi2;
  const double fwhmChi2 = m_fwhmEqnChi2;
  const vector<float> fwhmCoefs = m_fwhmCoefs;
  const vector<float> fwhmCoefUncert = m_fwhmCoefUncerts;
  const vector<float> effEqnCoefs = m_effEqnCoefs;
  const vector<float> effEqnCoefsUncerts = m_effEqnCoefUncerts;
  
  const vector<MakeDrfChart::DataPoint> data = m_chart->currentDataPoints();
  
  const int effDof = static_cast<int>(data.size()) - effEqnCoefs.size();
  const int fwhmDof = static_cast<int>(data.size()) - fwhmCoefs.size();
  
  const bool effInMeV = isEffEqnInMeV();
  const bool fwhmIsGadrasEqn = isGadrasFwhmEqnType();
  const DetectorPeakResponse::ResolutionFnctForm resFcnForm = (fwhmIsGadrasEqn
                                                               ? DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn
                                                               : DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial);
  
  double diam;
  try
  {
    diam = detectorDiameter();
  }catch( std::exception &e )
  {
    out << "Invalid detector diameter: " << e.what() << endline;
    return;
  }
  
  if( effEqnCoefs.empty() )
  {
    out << "Detector response function not valid" << endline;
    
    if( m_errorMsg )
      out << m_errorMsg->text().toUTF8() << endline;
    return;
  }//if( effEqnCoefs.empty() )
  
  //Okay, we've got all the variables we need for this function I think, lets
  //  write stuff
  
  const float cs137Energy = (effInMeV ? 0.661657f : 661.657f);
  const float intrinsicEffAt661 = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( cs137Energy, effEqnCoefs );
  
  float releffuncert = -999.9f;
  if( effEqnCoefs.size() == effEqnCoefsUncerts.size() && effEqnCoefs.size() > 0 )
    releffuncert = effEqnUncert( cs137Energy, effEqnCoefs, effEqnCoefsUncerts );
  

  out << "# Detector Response Function generated by InterSpec "
  << SpecUtils::to_iso_string( boost::posix_time::second_clock::local_time() )
  << endline << endline;
  
  out << "# Name,";
  output_csv_field( out, drfname );
  out << endline;
  out << "# Description,";
  output_csv_field( out, drfdescription );

  out << endline << endline
  << "# Intrinsic Efficiency Coefficients for equation of form"
  " Eff(x) = exp( C_0 + C_1*log(x) + C_2*log(x)^2 + ...) where x is energy in "
  << (effInMeV ? "MeV" : "keV") << endline
  << "#  i.e. equation for probability of gamma that hits the face of the detector being detected in the full energy photopeak." << endline
  << "# Name,Relative Eff @ 661keV,eff.c name,c0,c1,c2,c3,c4,c5,c6,c7,p0,p1,p2,Calib Distance,Radius (cm),G factor" << endline
  << drfname << " Intrinsic," << (100.0*intrinsicEffAt661/ns_NaI3x3IntrinsicEff) << "%,";
  for( size_t i = 0; i < effEqnCoefs.size(); ++i )
    out << "," << effEqnCoefs[i];
  for( size_t i = effEqnCoefs.size(); i < 12; ++i )
    out << ",";
  out << "0.0," << (0.5*diam/PhysicalUnits::cm) << ",0.5" << endline;
  out << "# 1 sigma Uncertainties,";
  if( releffuncert >= 0.0 )
    out << 100*(releffuncert / ns_NaI3x3IntrinsicEff) << "%";
  out << ",";
  for( size_t i = 0; i < effEqnCoefsUncerts.size(); ++i )
    out << "," << effEqnCoefsUncerts[i];
  if( effChi2 > 0 )
    out << endline << "# Chi2 / DOF = " << effChi2 << " / " << (effDof-1)
        << " = " << (effDof >= 1 ? (effChi2/(effDof-1.0)) : 0.0);
  out << endline << endline;
  
  //Convert equation into an absolute efficiency at 25 cm, and output those equations
  const double solidAngleAt25cm = DetectorPeakResponse::fractionalSolidAngle( diam, 25*PhysicalUnits::cm );
  out << "# Absolute Efficiency Coefficients (i.e. probability of gamma emitted from source at 25cm being detected in the full energy photopeak) for equation of form"
  " Eff(x) = exp( C_0 + C_1*log(x) + C_2*log(x)^2 + ...) where x is energy in "
  << (effInMeV ? "MeV" : "keV") << " and at distance of 25 cm" << endline
  << "#  i.e. equation for probability of gamma emitted from source at 25cm being detected in the full energy photopeak." << endline
  << "# Name,Relative Eff @ 661keV,eff.c name,c0,c1,c2,c3,c4,c5,c6,c7,p0,p1,p2,Calib Distance,Radius (cm),G factor" << endline
  << drfname << " Absolute," << (100.0*intrinsicEffAt661/ns_NaI3x3IntrinsicEff) << "%,";
  for( size_t i = 0; i < effEqnCoefs.size(); ++i )
    out << "," << ( (i==0 ? log(solidAngleAt25cm) : 0.0) + effEqnCoefs[i]);  //todo: make sure its not
  for( size_t i = effEqnCoefs.size(); i < 12; ++i )
    out << ",";
  out << "25," << (0.5*diam/PhysicalUnits::cm) << "," << solidAngleAt25cm << endline;
  out << "# 1 sigma Uncertainties,";
  if( releffuncert >= 0.0 )
    out << 100*(releffuncert / ns_NaI3x3IntrinsicEff) << "%";
  out << ",";
  for( size_t i = 0; i < effEqnCoefsUncerts.size(); ++i )
    out << "," << effEqnCoefsUncerts[i];
  if( effChi2 > 0 )
    out << endline << "# Chi2 / DOF = " << effChi2 << " / " << (effDof-1)
        << " = " << (effDof >= 1 ? (effChi2/(effDof-1.0)) : 0.0);
  out << endline << endline;
  
  
  //Then need to give FWHM equation form and coefficient, if available.
  if( fwhmCoefs.size() )
  {
    out << "# Full width half maximum (FWHM) follows equation: ";
    if( fwhmIsGadrasEqn )
    {
      out << "GadrasEqn" << endline
          << "# P6 -> resolution @ E=0 (keV)" << endline
          << "# P7 -> % FWHM @ 661 keV" << endline
          << "# P8 -> resolution power" << endline
          << "# if( energy >= 661 or P6=0 ) FWHM = 6.61×P7×(energy/661)^P8" << endline
          << "# else if( P6 < 0.0 ) FWHM = 6.61 x P7 x (energy/661)^(P8^(1.0/log(1.0-P6))" << endline
          << "# else if( P6 > 6.61×P7 ) FWHM = P6" << endline
          << "# else FWHM = sqrt(P6^2 + (6.61 x (sqrt((6.61×P7)^2 - P6^2)/6.61) x (energy/661)^P8)^2)" << endline
          << "# Energy in keV"<< endline
          << "# ,P6,P7,P8" << endline;
    }else
    {
      out << "FWHM = sqrt( A0 + A1*(energy/1000) + A2*(energy/1000)^2 + A3*(energy/1000)^3 + ...)" << endline
           << "# Energy in keV" << endline
           << "# ";
      for( size_t i = 0; i < fwhmCoefs.size(); ++i )
        out << ",A" << i;
      out << endline;
    }
    
    out << "Values";
    for( size_t i = 0; i < fwhmCoefs.size(); ++i )
      out << "," << fwhmCoefs[i];
    out << endline;
    if( fwhmCoefUncert.size() == fwhmCoefs.size() )
    {
      out << "Uncertainties";
      for( size_t i = 0; i < fwhmCoefUncert.size(); ++i )
        out << "," << fwhmCoefUncert[i];
    }
    if( fwhmChi2 > 0 )
      out << endline << "# Chi2 / DOF = " << fwhmChi2 << " / " << (fwhmDof-1)
          << " = " << (fwhmDof >= 1 ? (fwhmChi2/(fwhmDof-1.0)) : 0.0);
  }//if( fwhmCoefs.size() )
  
  //Give upper/lower energy range
  float lowerEnergy = 0.0f, upperEnergy = 0.0f;
  if( !data.empty() )
  {
    lowerEnergy = 99999.9f;
    upperEnergy = -99999.9f;
    //Data is not sorted
    for( const auto &p : data )
    {
      lowerEnergy = std::min( lowerEnergy, p.energy );
      upperEnergy = std::max( upperEnergy, p.energy );
    }
  }//if( data.size() >= 2 )
  
  out << endline << endline
      << "Detector diameter = " << (diam/PhysicalUnits::cm) << " cm." << endline
      << "Valid energy range: " << lowerEnergy
      << " keV to " << upperEnergy << " keV." << endline
      << endline;
  
  //Give infomation on the peaks and sources used to do the fit
  out << "# Peaks used to create DRF" << endline;
  
  out << "# Energy (keV),LiveTime (s),PeakArea,PeakArea Uncert,FWHM (keV),FWHM Uncert (keV),"
         "Source CPS,Source CPS Uncert,Distance (cm),SourceInfo,BackgroundPeakCounts,BackgroundLiveTime"
         ",DetectionEfficiency,DetectionEfficiencyUncert,GeometryFactor" << endline;
  for( MakeDrfChart::DataPoint d: data )
  {
    SpecUtils::ireplace_all( d.source_information, ",", " ");
    
    const double deteff = d.peak_area / d.source_count_rate;
    const double deteffUncert = deteff * sqrt( pow(d.peak_area_uncertainty/d.peak_area,2)
                                               + pow(d.source_count_rate_uncertainty/d.source_count_rate,2) );
    const double geomFactor = DetectorPeakResponse::fractionalSolidAngle(diam, d.distance);
    out << d.energy << "," << d.livetime
        << "," << d.peak_area << "," << d.peak_area_uncertainty
        << "," << d.peak_fwhm << "," << d.peak_fwhm_uncertainty
        << "," << d.source_count_rate << "," << d.source_count_rate_uncertainty
        << "," << d.distance / PhysicalUnits::cm
        << "," << d.source_information
        << "," << (d.background_peak_area > 0.0 ? std::to_string(d.background_peak_area) : string("") )
        << "," << (d.background_peak_live_time > 0.0 ? std::to_string(d.background_peak_live_time) : string("") )
        << "," << deteff
        << "," << deteffUncert
        << "," << geomFactor
        << endline;
  }//for( MakeDrfChart::DataPoint d: data )
  
  out << endline << endline;
  
  //Then should give GADRAS Efficiency.csv style efficiecny and FWHM every 50 keV
  out << "# Energy (keV),IntrinsicEfficiency,IntrinsicEfficiency Uncert,FWHM" << endline;
  for( int i = 25; i <= 3000; i += 25 )
  {
    const float energy = i / (effInMeV ? 1000.0f : 1.0f);
    const float eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, effEqnCoefs );
    const float effUncert = effEqnUncert( energy, effEqnCoefs, effEqnCoefsUncerts );
    
    out << i << "," << 100.0*eff << ",";
    
    if( effEqnCoefs.size() == effEqnCoefsUncerts.size() && effEqnCoefs.size() > 0 )
      out << 100.0*effUncert;
    out << ",";
    
    if( !fwhmCoefs.empty() )
      out << DetectorPeakResponse::peakResolutionFWHM( i, resFcnForm, fwhmCoefs );
    out << endline;
  }//for( int i = 25; i <= 3000; i += 25 )
  
}//void MakeDrf::writeCsvSummary( std::ostream & )




void MakeDrf::writeRefSheet( std::ostream &output, std::string drfname, std::string drfdescrip )
{
  // Use a lambda to read in the template XML, mostly just to control scope and have result be const
  auto get_tmplt_txt = []() -> string {
    const string docroot = wApp->docRoot();
    const string tmpltpath = SpecUtils::append_path(docroot,"InterSpec_resources/static_text/drf_ref_card.xml");
    
    std::ifstream tmpltfile( tmpltpath.c_str(), ios::in | ios::binary );
    if( !tmpltfile )
      return "Error opening reference card template text '" + tmpltpath + "'\n";
    
    std::stringstream tmplttxtstrm;
    tmplttxtstrm << tmpltfile.rdbuf();
    
    return tmplttxtstrm.str();
  };//get_tmplt_txt lamda
  
  // Lets gather all our general information
  const string tmplttxt = get_tmplt_txt();
  
  const double diam = [this]() -> double {try{return detectorDiameter();}catch(...){return 0.0;}}();
  const WString diameter = m_detDiameter->text();
  
  const int offset = wApp->environment().timeZoneOffset();
  const WDateTime now = WDateTime::currentDateTime().addSecs(60*offset);
  const string date = now.date().toString("MMM d yyyy").toUTF8();
  
  // Energy range, m_chartLowerE, m_chartUpperE
  string eff_eqn = m_intrinsicEffAnswer->text().toUTF8();
  //SpecUtils::ireplace_all( eff_eqn, "*", "&times;" );
  
  // Lets keep the equation from spilling over to the next line without some proper formatting
  if( eff_eqn.size() > 110 )
  {
    const string::size_type plus_pos = eff_eqn.find( " + ", 100 );
    const string::size_type minus_pos = eff_eqn.find( " - ", 100 );
    const string::size_type break_pos = (plus_pos < minus_pos) ? plus_pos : minus_pos;
    if( break_pos != string::npos )
    {
      eff_eqn = "<div>" + eff_eqn.substr( 0, break_pos ) + "</div>"
                "<div style=\"padding-left: 7em\">" + eff_eqn.substr( break_pos ) + "</div>";
    }
  }//if( eff_eqn.size() > 110 )
  
  

  const bool effInMeV = isEffEqnInMeV();
  const float cs137Energy = (effInMeV ? 0.661657f : 661.657f);
  const float intrinsicEffAt661 = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( cs137Energy, m_effEqnCoefs );
  const float relEffPercent = 100.0 * intrinsicEffAt661 / ns_NaI3x3IntrinsicEff;
  char rel_eff_txt[256] = { '\0' };
  
  if( m_fwhmCoefs.empty() )
  {
    snprintf( rel_eff_txt, sizeof(rel_eff_txt), "%.1f%% eff. (rel. to 3x3 NaI) @661 keV", relEffPercent );
  }else
  {
    const auto form = isGadrasFwhmEqnType() ? DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn
                                            : DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
    const float fwhm661 = DetectorPeakResponse::peakResolutionSigma( 661.7, form, m_fwhmCoefs );
    const float relResolution = 100 * fwhm661 / 661.7;
    
    snprintf( rel_eff_txt, sizeof(rel_eff_txt),
              "@661 keV: %.1f%% eff (rel. to 3x3 NaI), and FWHM=%.2g%%",
             relEffPercent, relResolution );
  }//if( we have FWHM ) / else.
  
  string fwhm_eqn;
  if( m_fwhmCoefs.size() )
  {
    if( isGadrasFwhmEqnType() )
    {
      const float P6 = m_fwhmCoefs[0];
      const float P7 = (m_fwhmCoefs.size() > 1) ? m_fwhmCoefs[1] : 0.0f;
      const float P8 = (m_fwhmCoefs.size() > 2) ? m_fwhmCoefs[2] : 0.0f;
      
      char case1_buffer[256] = { '\0' }, case2_buffer[256] = { '\0' };
      char case3_buffer[256] = { '\0' }, case4_buffer[512] = { '\0' };
      snprintf( case1_buffer, sizeof(case1_buffer), "FWHM(keV) = 6.61 * %.5g * (x/661)<sup>%.5g</sup>", P7, P8 );
      snprintf( case2_buffer, sizeof(case2_buffer), "FWHM(keV) = 6.61 * %.5g * (x/661)^(%.5g^(1.0/log(1.0-%.5g))", P7, P8, P6 );
      snprintf( case3_buffer, sizeof(case3_buffer), "FWHM(keV) = %.5g", P6 );
      snprintf( case4_buffer, sizeof(case4_buffer), "FWHM(keV) = sqrt(%.5g<sup>2</sup>"
               " + (6.61 *"
               " (sqrt((6.61 * %.5g)<sup>2</sup> - %.5g<sup>2</sup>)/6.61)"
               " * (x/661)<sup>%.5g</sup>)<sup>2</sup>)", P6, P7, P6, P8 );
      
      if( fabs(P6) < 1.0E-8 )
      {
        fwhm_eqn = case1_buffer;
      }else if( P6 < 0.0 )
      {
        fwhm_eqn = case2_buffer;
      }else if( P6 > 6.61*P7 )
      {
        fwhm_eqn = case3_buffer;
      }else
      {
        fwhm_eqn = "<div class=\"gadras-fcn\">E &lt; 661: " + string(case1_buffer) +  "</div>"
                   "<div class=\"gadras-fcn\">E &gt; 661: " + string(case4_buffer) +  "</div>";
      }
    }else
    {
      fwhm_eqn = "FWHM(keV) = sqrt(";
      for( size_t i = 0; i < m_fwhmCoefs.size(); ++i )
      {
        if( i == 0 )
          fwhm_eqn += ((m_fwhmCoefs[i] < 0.0) ? " -" : "");
        else
          fwhm_eqn += ((m_fwhmCoefs[i] < 0.0) ? " - " : " + ");
        
        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%.4g", 0.001*fabs(m_fwhmCoefs[i]) );
        fwhm_eqn += buffer;
        
        if( i == 1 )
          fwhm_eqn += "*x";
        else if( i > 1 )
          fwhm_eqn += "*x<sup>" + std::to_string(i) + "</sup>";
      }
      fwhm_eqn += ")";
    }//if( isGadrasFwhmEqnType() ) / else
  }//if( m_fwhmCoefs.size() )
  
  stringstream efftable;
  const double eff_energies[] = { 59.5, 185.7, 413.7, 661.7, 1001.0, 1460.8, 2614.5 };
  const int eff_dist_cm[] = { 25, 50, 75, 100, 250 };
  
  efftable << "<table class=\"eff-table\">\n"
              "\t\t<caption>Full-Energy Detection Efficiency (% total gammas)</caption>\n"
              "\t<thead>\n"
              "\t\t<tr><th></th>";
  char buffer[128] = { '\0' };
  for( const double energy : eff_energies )
  {
    if( (energy < (m_effLowerEnergy - 10)) || (energy > (m_effUpperEnergy + 10)) )
      continue;
    
    snprintf( buffer, sizeof(buffer), "<th>%.0f&nbsp;keV</th>", std::round(energy) );
    efftable << buffer;
  }
  efftable << "</tr>\n\t</thead>\n"
              "\t<tbody>\n";
  
  efftable << "\t<tr><th>Intrinsic</th>";
  for( double raw_energy : eff_energies )
  {
    if( (raw_energy < (m_effLowerEnergy - 10)) || (raw_energy > (m_effUpperEnergy + 10)) )
      continue;
    
    const float energy = raw_energy / (effInMeV ? 1000.0f : 1.0f);
    const float intrinsic_eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, m_effEqnCoefs );
    snprintf( buffer, sizeof(buffer), "<td>%.3G</td>", 100.0*intrinsic_eff );
    efftable << buffer;
  }
  efftable << "\t</tr>\n";
  
  for( const int dist_cm : eff_dist_cm )
  {
    efftable << "\t<tr><th>" << dist_cm << " cm</th>";
    for( double raw_energy : eff_energies )
    {
      if( (raw_energy < (m_effLowerEnergy - 10)) || (raw_energy > (m_effUpperEnergy + 10)) )
        continue;
      
      const float energy = raw_energy / (effInMeV ? 1000.0f : 1.0f);
      const float intrinsic_eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, m_effEqnCoefs );
      const double frac_solid_angle = DetectorPeakResponse::fractionalSolidAngle( diam, dist_cm*PhysicalUnits::cm );
    
      snprintf( buffer, sizeof(buffer), "<td>%.2G</td>", 100.0*frac_solid_angle*intrinsic_eff );
      efftable << buffer;
    }
    efftable << "\t</tr>\n";
  }//for( const auto dist : eff_dist )
  
  efftable << "\t</tbody>\n"
              "</table>\n";
  
  WLength chart_width(5*96,WLength::Unit::Pixel), chart_height(3*96,WLength::Unit::Pixel);
  
  Wt::WSvgImage eff_chart( chart_width, chart_height, nullptr, false );

  {// begin make DRF chart
    MakeDrfChart chart;
    chart.resize( chart_width, chart_height );
    const std::vector<MakeDrfChart::DataPoint> &data = m_chart->currentDataPoints();
    const float detDiam = m_chart->currentDiameter();
    const float chartLower = 10.0f * std::round( 0.1f*(m_effLowerEnergy - 6.0f) );
    const float chartUpper = 10.0f * std::round( 0.1f*(m_effUpperEnergy + 6.0f) );
    
    chart.setDataPoints( data, detDiam, chartLower, chartUpper );
    
    const auto units = effInMeV ? MakeDrfChart::EqnEnergyUnits::MeV
                                : MakeDrfChart::EqnEnergyUnits::keV;
    
    MakeDrfChart::FwhmCoefType fwhmEqnType = MakeDrfChart::FwhmCoefType::Gadras;
    switch( m_fwhmEqnType->currentIndex() )
    {
      case 0: fwhmEqnType = MakeDrfChart::FwhmCoefType::Gadras; break;
      case 1: fwhmEqnType = MakeDrfChart::FwhmCoefType::SqrtEqn; break;
    }//switch( m_fwhmEqnType->currentIndex() )
    
    chart.setFwhmCoefficients( m_fwhmCoefs, m_fwhmCoefUncerts, fwhmEqnType, units );
    chart.setEfficiencyCoefficients( m_effEqnCoefs, m_effEqnCoefUncerts, units );
    
    WPainter p( &eff_chart );
    chart.paint( p );
    p.end();
  }// end make DRF chart
  
  stringstream eff_chart_svg;
  eff_chart.write( eff_chart_svg );
  
  
  WTemplate tmplt;
  // tmplttxt is unsafe HTML because of <head> tag and stuff, so we dont want Wt to filter it
  tmplt.setTemplateText( WString::fromUTF8(tmplttxt), TextFormat::XHTMLUnsafeText );
  tmplt.bindString("drf-name", drfname, TextFormat::XHTMLUnsafeText );
  tmplt.bindString("drf-desc", drfdescrip, TextFormat::XHTMLUnsafeText );
  tmplt.bindString("diameter", diameter, TextFormat::XHTMLUnsafeText );
  tmplt.bindString("date", date, TextFormat::XHTMLUnsafeText );
  tmplt.bindString("rel-eff", rel_eff_txt, TextFormat::XHTMLUnsafeText );
  tmplt.bindString("eff-eqn", eff_eqn, TextFormat::XHTMLUnsafeText );
  tmplt.bindString("fwhm-eqn", fwhm_eqn, TextFormat::XHTMLUnsafeText );
  tmplt.bindString("eff-table", efftable.str(), TextFormat::XHTMLUnsafeText );
  tmplt.bindString("eff-svg", eff_chart_svg.str(), TextFormat::XHTMLUnsafeText );
  
  tmplt.renderTemplate( output );
}//void writeRefSheet(...)


bool MakeDrf::isEffEqnInMeV() const
{
  return (m_effEqnUnits->currentIndex() == 1);
}//bool isEffEqnInMeV() const


bool MakeDrf::isGadrasFwhmEqnType() const
{
  switch( m_fwhmEqnType->currentIndex() )
  {
    case 0: return true;
    case 1: return false;
  }//switch( m_fwhmEqnType->currentIndex() )
  
  return false;
}//bool isGadrasFwhmEqnType() const


double MakeDrf::detectorDiameter() const
{
  const double diam = PhysicalUnits::stringToDistance( m_detDiameter->text().toUTF8() );
  if( diam <= 0.0 )
    throw runtime_error( "Detector diameter less than or equal to zero." );
  return diam;
}//double detectorDiameter() const
