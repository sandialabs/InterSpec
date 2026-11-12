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

#include <map>
#include <set>
#include <regex>
#include <deque>
#include <limits>
#include <fstream>
#include <algorithm>

#include <Wt/Utils.h>
#include <Wt/WText.h>
#include <Wt/WPoint.h>
#include <Wt/WImage.h>
#include <Wt/WPanel.h>
#include <Wt/WLabel.h>
#include <Wt/WTable.h>
#include <Wt/WAnchor.h>
#include <Wt/WServer.h>
#include <Wt/WPainter.h>
#include <Wt/WResource.h>
#include <Wt/WComboBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WCheckBox.h>
#include <Wt/WSvgImage.h>

#include <Wt/WTemplate.h>
#include <Wt/WTextArea.h>
#include <Wt/WTableCell.h>
#include <Wt/WFormWidget.h>
#include <Wt/WTabWidget.h>
#include <Wt/WIOService.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WEnvironment.h>
#include <Wt/WApplication.h>
#include <Wt/WPopupWidget.h>
#include <Wt/Http/Response.h>
#include <Wt/WDoubleSpinBox.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WSuggestionPopup.h>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/QrCode.h"
#include "InterSpec/MakeDrf.h"
#include "InterSpec/MakeMcResponseForDrf.h"

// CeeLo (external_libs/CeeLo/src): MC-parameterized response + grounding.
#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/MakeDrfChart.h"
#include "InterSpec/DrfChart.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MakeDrfCalc.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/DetectorGeometryInput.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/PeakFitDetPrefsGui.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/GroupBox.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhysicalUnitsLocalized.h"


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
  /** Sets a flag for the life of the scope; used to block a re-fit while the source widgets are
   only part-way through a bulk update.
   */
  struct BlockSourceUpdates
  {
    bool &m_flag;
    BlockSourceUpdates( bool &flag ) : m_flag( flag ){ m_flag = true; }
    ~BlockSourceUpdates(){ m_flag = false; }
  };//struct BlockSourceUpdates


  /** Width of the "Make Detector Response" window, for an app `app_width` pixels wide. */
  int make_drf_dialog_width( const int app_width )
  {
    return std::min( 1200, (app_width > 840) ? (9*app_width)/10 : (app_width - 20) );
  }


  /** The coarse kind of detector a crystal material is - for the ones that are clear-cut; Unknown
   for the rest.  `name` is one of `CeeLoUtils::crystalMaterialNames()`. */
  PeakFitUtils::CoarseResolutionType crystal_resolution_type( std::string name )
  {
    SpecUtils::to_lower_ascii( name );
    if( (name == "hpge") || (name == "ge") )
      return PeakFitUtils::CoarseResolutionType::High;
    if( (name == "nai") || (name == "csi") || (name == "bgo") || (name == "pvt") )
      return PeakFitUtils::CoarseResolutionType::Low;
    if( SpecUtils::istarts_with( name, "labr" ) || SpecUtils::istarts_with( name, "cebr" )
       || SpecUtils::istarts_with( name, "lacl" ) )
      return PeakFitUtils::CoarseResolutionType::LaBr;
    if( (name == "czt") || (name == "cdte") )
      return PeakFitUtils::CoarseResolutionType::CZT;
    return PeakFitUtils::CoarseResolutionType::Unknown;
  }//crystal_resolution_type(...)


  /** Whether two detector kinds plainly disagree: one is HPGe and the other is a scintillator or
   CZT.  The fuzzier kinds (medium resolution, low-or-medium, unknown) never disagree. */
  bool resolution_types_conflict( const PeakFitUtils::CoarseResolutionType a,
                                  const PeakFitUtils::CoarseResolutionType b )
  {
    auto is_hpge = []( PeakFitUtils::CoarseResolutionType t ){
      return (t == PeakFitUtils::CoarseResolutionType::High);
    };
    auto is_not_hpge = []( PeakFitUtils::CoarseResolutionType t ){
      return (t == PeakFitUtils::CoarseResolutionType::Low)
             || (t == PeakFitUtils::CoarseResolutionType::LaBr)
             || (t == PeakFitUtils::CoarseResolutionType::CZT);
    };
    return (is_hpge(a) && is_not_hpge(b)) || (is_not_hpge(a) && is_hpge(b));
  }//resolution_types_conflict(...)


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
      
      const SpecUtils::time_point_t datestr = SpecUtils::time_from_string( fields[2] );
      activityDate = to_ptime(datestr);
      
      comments = "";
      for( size_t i = 3; i < fields.size(); ++i )
        comments += ((i==3) ? "" : " ") + fields[i];
      
      return true;
    }//while( SpecUtils::safe_get_line( file, line) )
    
    return false;
  }//source_info_from_lib_file(...)


  /** Returns the distance given after an '@' in a spectrum title, or -1.0 if the title doesnt
   have one.

   The spacing around the '@', and between the value and its units, is arbitrary, and any text
   trailing the distance is ignored, so for example all of "U-232 @ 100 cm, H=100 cm",
   "60Co,58.599uCi @ 100cm H=83.5cm", and "Cs137@1m 50cm" give their distance.
   */
  double distance_from_spectrum_title( const std::string &title )
  {
    // PhysicalUnits::sm_distanceRegex is anchored to match an entire string; strip the anchors so
    //  we can search for the distance embedded in the title.
    string dist_regex = PhysicalUnits::sm_distanceRegex;
    SpecUtils::ireplace_all( dist_regex, "^", "" );
    SpecUtils::ireplace_all( dist_regex, "$", "" );

    try
    {
      const std::regex expr( "@\\s*(" + dist_regex + ")", std::regex::icase );

      std::smatch mtch;
      if( !std::regex_search( title, mtch, expr ) )
        return -1.0;

      const double distance = PhysicalUnits::stringToDistance( mtch[1].str() );
      if( (distance > 0.0) && !IsNan(distance) && !IsInf(distance) )
        return distance;
    }catch( std::exception & )
    {
      //stringToDistance(...) can throw, though shouldnt, given the regex matched
    }

    return -1.0;
  }//distance_from_spectrum_title(...)


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
    string m_description;
    Wt::WApplication *m_app;
    
  public:
    CalFileDownloadResource( const bool pcf, MakeDrf *parent )
     : Wt::WResource(), m_pcf(pcf), m_makedrf( parent ), m_filename( "" ),
       m_app( WApplication::instance() )
    {
      assert( m_app );

      // Have Wt hold the session lock across handleRequest instead of dropping it for us: by
      //  default `WResource::handle` unlocks before invoking us, so the session thread can be
      //  inside `~WResource`'s `beingDeleted()` - waiting on our use-count - while we block
      //  acquiring the lock it holds.  That is a deadlock.  Free: the `WApplication::UpdateLock`
      //  in handleRequest already spanned the whole body, so this only moves it earlier.
      setTakesUpdateLock( true );
    }
    
    virtual ~CalFileDownloadResource()
    {
      beingDeleted();
    }
    
    void setSuggestFileName( string filename )
    {
      m_filename = filename;
    }

    void setDescription( string description )
    {
      m_description = description;
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
      
      shared_ptr<SpecMeas> calfile = m_makedrf->assembleCalFile( m_filename, m_description );
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
      
      suggestFileName( filename, ContentDisposition::Attachment );
      response.setMimeType( "application/octet-stream" );
      
      if( m_pcf)
        calfile->write_pcf( response.out() );
      else
        calfile->write_2012_N42( response.out() );
    }
  };//class CalFileDownloadResource
  
  
  /** Uncertainty of the efficiency equation at `energy` (in the equation's units); negative if it
   can not be computed.  With the coefficient covariance (row-major) this is linear propagation;
   without it, only the individual uncertainties are combined - which greatly overstates it, since
   the coefficients are strongly correlated. */
  float effEqnUncert( const float energy, vector<float> coefs, const vector<float> &uncerts,
                      const vector<float> &cov )
  {
    assert( coefs.size() == uncerts.size() );
    
    const float eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, coefs );

    if( cov.size() == coefs.size()*coefs.size() )
    {
      const double uncert = DetectorPeakResponse::expOfLogPowerSeriesEfficiencyUncert( energy, coefs, cov );
      return (uncert >= 0.0) ? static_cast<float>( uncert ) : -999.9f;
    }
    
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
    : WResource(), m_filename(""), m_makedrf( parent ), m_app( WApplication::instance() )
    {
      assert( m_app );

      // Have Wt hold the session lock across handleRequest instead of dropping it for us: by
      //  default `WResource::handle` unlocks before invoking us, so the session thread can be
      //  inside `~WResource`'s `beingDeleted()` - waiting on our use-count - while we block
      //  acquiring the lock it holds.  That is a deadlock.  Free: the `WApplication::UpdateLock`
      //  in handleRequest already spanned the whole body, so this only moves it earlier.
      setTakesUpdateLock( true );
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
      filename += fileNameSeparator();
      filename += fileEnding();

      // The DRF name may hold characters a file name cant (as DrfLibraryWidgets' export replaces)
      for( char &c : filename )
      {
        if( string("\\/:?\"<>|*").find( c ) != string::npos )
          c = '_';
      }

      suggestFileName( filename, ContentDisposition::Attachment );
    }

    /** Between the DRF name and #fileEnding in the suggested file name. */
    virtual const char *fileNameSeparator() const { return "_"; }
    
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
      suggestFileName( fileEnding(), ContentDisposition::Attachment );
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
      suggestFileName( fileEnding(), ContentDisposition::Attachment );
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


  /** The DRF as InterSpec's own XML - what the DRF library stores, and what "Detector Response
   Select" imports - named "<name>.drf.xml" like the library's own export. */
  class DrfXmlDownloadResource : public DrfSummaryBase
  {
  public:
    DrfXmlDownloadResource( MakeDrf *parent )
     : DrfSummaryBase( parent )
    {
      suggestFileName( fileEnding(), ContentDisposition::Attachment );
    }

    virtual const char *fileEnding() const { return "drf.xml"; }
    virtual const char *fileNameSeparator() const { return "."; }

    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in DrfXmlDownloadResource.";
        response.out() << "Error grabbing application lock to form DrfXmlDownloadResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        return;
      }//if( !lock )

      try
      {
        const shared_ptr<DetectorPeakResponse> drf
             = m_makedrf->assembleDrf( m_filename.empty() ? string("Created DRF") : m_filename,
                                       m_description );
        if( !drf || !drf->isValid() )
          throw runtime_error( "the response is not valid" );

        response.setMimeType( "application/octet-stream" );
        rapidxml::xml_document<char> doc;
        drf->toXml( &doc, &doc );
        rapidxml::print( response.out(), doc, 0 );
      }catch( std::exception &e )
      {
        log("error") << "Failed to write DRF XML: " << e.what();
        response.setMimeType( "text/plain" );
        response.setStatus( 500 );
        response.out() << "Could not create the detector response: " << e.what();
      }//try / catch
    }//void handleRequest(...)
  };//class DrfXmlDownloadResource
  

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
             WContainerWidget * /*parent*/ = nullptr )
    : WContainerWidget(),
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
      m_useForEffCb = addNew<WCheckBox>( WString::tr("dp-use-cb") );
      m_useForEffCb->addStyleClass( "DrfPeakUseCb" );
      
      char buffer[256];
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
        m_userBr->setMargin( 5, Wt::Side::Left );
        m_userBr->setTextSize( 8 );
        
        m_useForEffCb->setUnChecked();
        m_userBr->hide();
        m_useForEffCb->checked().connect( this, [this](){
          m_userBr->setHidden( !m_useForEffCb->isChecked() || m_isBackground );
        } );
        m_useForEffCb->unChecked().connect( this, [this](){
          m_userBr->setHidden( !m_useForEffCb->isChecked() || m_isBackground );
        } );
      }
      
      m_descTxt = addNew<WText>( WString::fromUTF8(buffer) );
      // This is totally bizarre, but Safari (version 16.5 at least), and Chrome, will crash
      //  (the rendering process for the tab, by running up huge amounts of memory), the
      //  second time we create the "Make Detector Response" window, under some
      //  circumstances (seems to happen most reliably with multiple spectrum files open
      //  with lots of peaks).  If we either limit the text length of `m_descTxt`, or
      //  explicitly set its width to be less than 100%, then this crash doesnt happen.
      //  Again, totally bizarre, and took quiet a while to pin down.
      //  Perhaps there is something else going on???  Has to be, right?
      m_descTxt->setWidth( WLength(70,WLength::Unit::Percentage) );
      m_descTxt->addStyleClass( "DrfPeakInfoTxt" );
      m_backSubTxt = addNew<WText>( WString::fromUTF8("") );
      m_backSubTxt->addStyleClass( "DrfPeakBackSubTxt" );
      m_backSubTxt->hide();
      
      if( m_userBr )
      {
        WLabel *l = addNew<WLabel>( WString::fromUTF8("BR=") );
        l->addStyleClass( "DrfPeakNoNucBrLbl" );
        l->setBuddy( m_userBr );
        addWidget( std::unique_ptr<WWidget>(m_userBr) );
      }

      m_previewBtn = addNew<WPushButton>( WString::fromUTF8("") );
      m_previewBtn->setStyleClass( "DrfPeakShowChart" );
      m_previewBtn->setIcon( "InterSpec_resources/images/peak_small.png" );
      //m_previewBtn = new WImage( "InterSpec_resources/images/peak_small.png", this );
      
      m_previewBtn->clicked().connect( this, [this]( WMouseEvent a1 ){ togglePeakPreview( a1 ); } );
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
      if( !wasSubtractedFrom && !beingsubtracted )
        return;

      // The rate can change while it stays subtracted (e.g., a background peak is toggled)
      if( beingsubtracted )
      {
        char buffer[64];
        snprintf( buffer, sizeof(buffer)-1, "%.3f", background_cps );
        m_backSubTxt->setText( WString::tr("dp-back-sub-txt").arg(buffer) );
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
        m_previewBtn->removeStyleClass( "active", true );
        m_previewPopup->setHidden( true );
      }
    }//void hidePeakPreview()
    
    void togglePeakPreview( WMouseEvent event )
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
          m_previewBtn->removeStyleClass( "active", true );
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
        
        auto chartWidget = std::make_unique<WText>( strm.str(), Wt::TextFormat::UnsafeXHTML );
        chartWidget->addStyleClass( "DrfPeakChart" );
        m_previewPopup = addNew<WPopupWidget>( std::move(chartWidget) );

        //ToDo: add a "X" in the corner so users can close the showing window.
      }catch( std::exception & )
      {
        m_previewPopup = addNew<WPopupWidget>( std::make_unique<WText>( WString::tr("dp-no-peak-preview") ) );
      }//try / catch make a chart
      
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
      m_previewPopup->setAnchorWidget( m_previewBtn, Wt::Orientation::Vertical );
      m_previewPopup->setJavaScriptMember("wtNoReparent", "true");
      m_previewPopup->show();
#else
      // In Wt 3.7.1, it looks like there is no way to keep the JS from "re-parenting" the
      //  popup-widget when you call positionAtWidget (See positionAtWidget in Wt.js); Wt will
      //  re-parent the widget to the first ancestor of the anchor widget with scrollbars or
      //  something.
      //  This causes our window to then expand to scroll left/right, which is not what we want.
      //  So we'll hack around this by just positioning the popup at the mouse click.
      m_previewPopup->show();
      
      // Positioning taken from WPopupMenu.C
      m_previewPopup->setOffsets(42, Wt::Side::Left | Wt::Side::Top);
      m_previewPopup->setOffsets(-10000, Wt::Side::Left | Wt::Side::Top);
      doJavaScript(WT_CLASS ".positionXY('" + m_previewPopup->id() + "',"
                   + std::to_string(event.window().x) + ","
                   + std::to_string(event.window().y) + ");");
      
      m_previewPopup->setTransient( true, 0 );
      // connecting to WWebWidget::removeStyleClass(...) fails for some reason I cant quite tell why,
      //  and I dont want to use the std::bind(lamda) version for possible life-time issues, so we'll
      //  just remove the "active" style class in JS - has the downside if we then show the chart
      //  again, the "active" style class wont be added since then the C++ and JS are out of sync,
      //  so to fix this up, we'll also add the "active" style class in JS when its shown.  Of
      //  course now the roundtrips from JS -> C++ -> make notable delays, but this is a minor
      //  detail for the moment.
      //m_previewPopup->hidden().connect( boost::bind( &WPushButton::removeStyleClass, m_previewBtn, WString("active"), true ) );
      //m_previewPopup->hidden().connect( std::bind([this](){ m_previewBtn->removeStyleClass("active",true); }));
      m_previewPopup->hidden().connect( this, [btn = m_previewBtn, js = "document.getElementById('" + m_previewBtn->id() + "').classList.remove('active');"](){
        btn->doJavaScript( js );
      } );
      m_previewPopup->shown().connect( this, [btn = m_previewBtn, js = "document.getElementById('" + m_previewBtn->id() + "').classList.add('active');"](){
        btn->doJavaScript( js );
      } );
#endif
      
      m_previewBtn->addStyleClass( "active" );
      m_peakPreviewShow.emit( this );
    }//void togglePeakPreview()
    
  };//class DrfPeak
  
  
  class DrfSpecFileSample : public WPanel
  {
    std::shared_ptr<SpecMeas> m_meas;
    set<int> m_samples;
    double m_livetime;  ///< Summed over the samples and detectors
    double m_realtime;  ///< Summed over the samples; the detectors count at the same time

    WContainerWidget *m_peaks;
    WContainerWidget *m_sources;
    WCheckBox *m_background;
    WText *m_backgroundTxt;
    WCheckBox *m_allNoneSome;
    Wt::Signal<> m_srcInfoUpdated;
    DetectorPeakResponse::EffGeometryType m_geometry_type;
    
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
        title = filename + (filename.empty() ? "" : ": ") + WString::tr("dsfs-sample").toUTF8()
                  + " " + std::to_string(*m_samples.begin()) + ", &quot;" + title + "&quot;";
      }else if( m_samples.size() == 1 && title.empty() )
      {
        title = filename + (filename.empty() ? " " : ":  ") + WString::tr("dsfs-sample").toUTF8()
                + " " + std::to_string(*m_samples.begin());
      }else if( title.empty() )
      {
        title = filename + (filename.empty() ? " " : ":  ") + WString::tr("dsfs-samples").toUTF8() + " {";
        for( auto iter = begin(m_samples); iter != end(m_samples); ++iter )
          title += (iter==begin(m_samples) ? "" : ",") + std::to_string(*iter);
        title += "}";
      }else  //we have a title and multiple sample numbers
      {
        title = filename + ": &quot;" + title + "&quot; " + WString::tr("dsfs-samples").toUTF8() + " {";
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
      
      title += WString::tr("dsfs-using-num-peaks").arg(npeaksused).arg(npeaks).toUTF8();
      
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
                      WContainerWidget * /*parent*/ = nullptr )
    : WPanel(),
      m_meas( meas ),
      m_samples( samples ),
      m_livetime( 0.0 ),
      m_realtime( 0.0 ),
      m_peaks( nullptr ),
      m_sources( nullptr ),
      m_background( nullptr ),
      m_backgroundTxt( nullptr ),
      m_allNoneSome( nullptr ),
      m_geometry_type( DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      
      addStyleClass( "DrfSpecFileSample" );
      
      auto bodyOwned = std::make_unique<WContainerWidget>();
      WContainerWidget *body = bodyOwned.get();
      setCentralWidget( std::move(bodyOwned) );

      m_peaks = body->addNew<WContainerWidget>();

      WContainerWidget *bckgrndDiv = body->addNew<WContainerWidget>();
      bckgrndDiv->addStyleClass( "DrfSpecFileSampleBackground" );
      m_background = bckgrndDiv->addNew<WCheckBox>( WString::tr("dsfs-is-background-cb") );
      m_background->checked().connect( this, &DrfSpecFileSample::isBackgroundToggled );
      m_background->unChecked().connect( this, &DrfSpecFileSample::isBackgroundToggled );
      m_background->setChecked( false );
      m_backgroundTxt = bckgrndDiv->addNew<WText>( WString::tr("dsfs-back-peak-txt") );
      m_backgroundTxt->addStyleClass( "DrfSpecFileSampleBackgroundTxt" );
      m_backgroundTxt->hide();

      m_sources = body->addNew<WContainerWidget>();
      m_sources->addStyleClass( "DrfSampleSources" );
      
      setCollapsible( true );
      
      std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaksptr = meas->peaks( samples );
      
      //Make a widget representing sample number combination with peaks in file
      auto peaks = make_shared< deque<shared_ptr<const PeakDef>> >( *peaksptr );
      assert( peaks && !peaks->empty() );
      
      string title;
      //ToDo: The peaks may not belong to all the detectors (the user may have
      //      chosen to not show some detectors)
      const vector<string> &detnames = meas->detector_names();
      for( const int sample : samples )
      {
        double sample_realtime = 0.0;
        for( const string &detname : detnames )
        {
          shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, detname );
          if( m && !m->title().empty() )
            title = m->title();
          if( m && m->num_gamma_channels() > 5 )
          {
            m_livetime += m->live_time();
            sample_realtime = std::max( sample_realtime, static_cast<double>(m->real_time()) );
          }
        }//
        m_realtime += sample_realtime;
      }//for( const int sample : samples )
      
      std::shared_ptr<const SpecUtils::Measurement> summed_meas;
      if( samples.size()==1 && detnames.size()==1 )
        summed_meas = meas->measurement( *begin(samples), detnames[0] );
      if( !summed_meas )
        summed_meas = meas->sum_measurements( samples, meas->detector_names(), nullptr );
      
      //ToDo: in principle, summed_meas->live_time() should equal livetime; should check
      
      int npeaks = 0;
      for( auto peak : *peaks )
      {
        // Wt4: the DrfPeak ctor ignores its (commented-out) parent arg, so `new DrfPeak(...)`
        //  would be orphaned (leaked + never rendered).  addNew<> constructs and parents it to m_peaks.
        DrfPeak *p = m_peaks->addNew<DrfPeak>( peak, m_livetime, summed_meas );
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
        m_allNoneSome = titleBarWidget()->addNew<WCheckBox>( WString::tr("dsfs-all-peaks-cb") );
        m_allNoneSome->clicked().preventPropagation();
        m_allNoneSome->addStyleClass( "DrfSpecFileAllNoneSomeCb" );
        m_allNoneSome->checked().connect( this, &DrfSpecFileSample::handleUserToggleAllNoneSum );
        m_allNoneSome->unChecked().connect( this, &DrfSpecFileSample::handleUserToggleAllNoneSum );
      }
      
      
      refreshSourcesVisible();
      
      // Now loop over sources and check for information in the remarks about
      // the source - we assume the GADRAS convention for things.
      
      
      auto try_set_info = [meas,samples,db]( MakeDrfSrcDef *srcdef ){
        
        assert( srcdef );
        assert( srcdef->nuclide() );
        assert( !srcdef->nuclide()->isStable() );
        
        string shielding;
        boost::posix_time::ptime activityDate;
        double title_distance = -1.0;
        double distance = -1.0, activity = -1.0, age_at_meas = -1.0;
        
        bool have_set_nuc_info = false;
        
        string srcinfo;
        const vector<string> &dets = meas->detector_names();
        for( const string &det : dets )
        {
          if( have_set_nuc_info )
            break;
          
          for( const int sample : samples )
          {
            if( have_set_nuc_info )
              break;
            
            auto m = meas->measurement(sample, det);
            if( !m )
              continue;
            
            string spectitle = m->title();
            
            //Look for a distance after an '@' in the title. Examples:
            //  "U-232 @ 100 cm, H=100 cm"
            //  "60Co,58.599uCi @ 100cm H=83.5cm"   (the H=... measurement height is not used)
            //  "Background, H=100 cm"              (no distance)
            if( title_distance <= 0.0 )
              title_distance = distance_from_spectrum_title( spectitle );
            
            //ToDo: make this more robust for finding end of the nuclide, like requiring both numbers and letters.
            
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
                  age_at_meas = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( mtch[2].str(), srcdef->nuclide()->halfLife );
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
              
              
              
              /*
               //Untested regex equivalent (or maybe a bit stricter) of the above.
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
              
              
              // This next line wont give a nuclide if there is anything invalid trailing the nuclide (ex, 'Ba133,10uCi')
              const SandiaDecay::Nuclide *nuc = db->nuclide( remark );
              
              size_t commapos = remark.find(',');
              size_t underscorpos = remark.find('_'); //For a source database name like "133BA_<serial number>"
              
              if( !nuc && (commapos != string::npos) )
                nuc = db->nuclide( remark.substr(0,commapos) );
              
              if( !nuc && (underscorpos != string::npos) )
                nuc = db->nuclide( remark.substr(0,underscorpos) );
              
              if( !nuc || (nuc != srcdef->nuclide()) )
                continue;
              
              have_set_nuc_info = true;
              
              const size_t openCurlyPos = remark.find( '{' );
              if( shielding.empty() && (openCurlyPos != string::npos) )
              {
                const size_t closeCurlyPos = remark.find( '}', openCurlyPos );
                if( (closeCurlyPos != string::npos) && ((openCurlyPos+4)<closeCurlyPos) )
                {
                  shielding = remark.substr(openCurlyPos+1, closeCurlyPos-openCurlyPos-1);
                  // Keep what follows the shielding too - MakeDrfSrcDef::toGadrasLikeSourceString
                  //  writes " Age=... Distance=..." after it.
                  remark = remark.substr(0,openCurlyPos) + " " + remark.substr(closeCurlyPos+1);
                  SpecUtils::trim(remark);
                  
                  commapos = remark.find(',');
                  underscorpos = remark.find('_'); //For a source database name like "133BA_<serial number>"
                }
              }//if( source had shielding defined )
              
              
              distance = title_distance;
              
              std::smatch dist_mtch;
              
              // By default PhysicalUnits::sm_distanceRegex has a "^" character at begining, and "$"
              //  character at end - lets get rid of these
              string dist_regex = PhysicalUnits::sm_distanceRegex;
              SpecUtils::ireplace_all( dist_regex, "^", "" );
              SpecUtils::ireplace_all( dist_regex, "$", "" );
              
              std::regex dist_expr( string(".*\\b((?:distance|dist)\\s*\\=\\s*(") + dist_regex + ")).*?", std::regex::icase );
              
              if( std::regex_match( remark, dist_mtch, dist_expr ) )
              {
                try
                {
                  double trial_distance = PhysicalUnits::stringToDistance( dist_mtch[2].str() );
                  if( trial_distance >= 0.0 && !IsNan(trial_distance) && !IsInf(trial_distance) )
                    distance = trial_distance;
                }catch( std::exception & )
                {
                  cerr << "Failed to convert '" << dist_mtch[2].str() << "' to a distance." << endl;
                }
              }//if( std::regex_match( remark, dist_mtch, dist_expr ) )
              
              
              if( underscorpos != string::npos )
              {
                auto pos = remark.find_first_of( " \t", underscorpos );
                string src_name = remark.substr(0,pos); //fine if pos==string::npos
                SpecUtils::trim( src_name );
                
                for( const SrcLibLineInfo &src : srcdef->lib_srcs_for_nuc() )
                {
                  if( SpecUtils::icontains( src_name, src.m_source_name)  )
                  {
                    activityDate = src.m_activity_date;
                    activity = src.m_activity;
                    break;
                  }
                }//for(
              }//if( maybe a source name from a Source.lib )
              
              if( (commapos < underscorpos) && (activity <= 0.0) )
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
              }//if( commapos < underscorpos )
            }//for( string remark : m->remarks() )
          }//for( const int sample : samples )
        }//for( const string &det : dets )
        
        if( distance > 0.0 )
          srcdef->setDistance( distance );
        else if( title_distance > 0.0 )
          srcdef->setDistance( title_distance );
        
        if( (activity > 0.0) && activityDate.is_special() )
          srcdef->setActivity( activity );
        else if( activity > 0.0 )
          srcdef->setAssayInfo( activity, activityDate );
        
        if( age_at_meas >= 0.0 )
          srcdef->setAgeAtMeas( age_at_meas );
        
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
              srcdef->setShielding( an_ad[0], an_ad[1]*PhysicalUnits::gram/PhysicalUnits::cm2 );
            }
          }catch( std::exception &e )
          {
            cerr << "DrfSpecFileSample: Caught exception setting shielding from '" << shielding << "': " << e.what() << endl;
          }
        }//if( !shielding.empty() )
        
      };// try_set_info lamda
      
      
      bool is_background = false;
      if( samples.size() == 1 )
      {
        for( const string &det : meas->detector_names() )
        {
          for( const int sample : samples )
          {
            auto m = meas->measurement(sample, det);
            if( !m )
              continue;
            
            is_background |= (m->source_type()==SpecUtils::SourceType::Background);
            is_background |= SpecUtils::icontains( m->title(), "back" );
            is_background |= SpecUtils::icontains( m->title(), "bgr" );
          }//for( const int sample : samples )
        }//for( const string &det : dets )
      }//if( samples.size() == 1 )
      
      
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
      }else
      {
        for( auto w : m_sources->children() )
        {
          auto src = dynamic_cast<MakeDrfSrcDef *>(w);
          if( src && src->nuclide() )
            try_set_info( src );
        }//for( auto w : m_sources->children() )
      }//if( is_background ) / else
      
    }//DrfSpecFileSample constructor
    
    void handleUserToggleAllNoneSum()
    {
      assert( m_allNoneSome );
      
      bool useAll = true;
      switch( m_allNoneSome->checkState() )
      {
        case Wt::CheckState::Checked:          useAll = true;  break;
        case Wt::CheckState::Unchecked:        useAll = false; break;
        case Wt::CheckState::PartiallyChecked: useAll = false; break;
      }//switch( m_allNoneSome->checkState() )
      
      for( auto w : m_peaks->children() )
      {
        auto p = dynamic_cast<DrfPeak *>( w );
        if( p )
          p->setUse( useAll );
      }//for( auto w : m_peaks->children() )
      
      refreshSourcesVisible();
    }//void handleUserToggleAllNoneSum()
    
    void setIsEffGeometryType( const DetectorPeakResponse::EffGeometryType type )
    {
      m_geometry_type = type;
      for( auto w : m_sources->children() )
      {
        MakeDrfSrcDef *src = dynamic_cast<MakeDrfSrcDef *>( w );
        if( src )
          src->setIsEffGeometryType( static_cast<int>(type) );
      }//for( auto w : m_sources->children() )
      
      refreshSourcesVisible(); //I dont actually think this is necassary, but JIC
    }//void setIsFixedGeometry( const bool isFixedGeometry )
    
    
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
              if( m && !SpecUtils::is_special(m->start_time()) )
              {
                measDate = to_ptime(m->start_time());
                break;
              }
            }//for( const int sample : m_samples )
          }//for( size_t i = 0; i < detNames.size(); ++i )
          
          // Wt4: MakeDrfSrcDef has no parent ctor arg, so `new MakeDrfSrcDef(...)` would be orphaned
          //  (leaked + never rendered).  addNew<> constructs and parents it to m_sources.
          MakeDrfSrcDef *src = m_sources->addNew<MakeDrfSrcDef>( n, measDate );
          src->setIsEffGeometryType( static_cast<int>(m_geometry_type) );
          src->updated().connect( this, [this](){ m_srcInfoUpdated.emit(); } );
          
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
    
    
    void addSourceLibrary( const vector<shared_ptr<const SrcLibLineInfo>> &srcs, 
                          const bool auto_populate )
    {
      for( MakeDrfSrcDef *src : sources() )
        src->addSourceLibrary( srcs, auto_populate );
    }//void addSourceLibrary(...)
    
    
    Wt::Signal<> &srcInfoUpdated(){ return m_srcInfoUpdated; };
    
    std::shared_ptr<const SpecMeas> measurement() { return m_meas; }
    const set<int> &samples(){ return m_samples; }
    double liveTime() const { return m_livetime; }
    double realTime() const { return m_realtime; }

    /** Where this sample is, for messages: "file.n42", or "file.n42, Sample 3" when the file has
     more than one sample.  Plain text. */
    std::string locationLabel() const
    {
      string answer = m_meas ? SpecUtils::filename( m_meas->filename() ) : string();
      if( m_meas && (m_meas->sample_numbers().size() > 1) && !m_samples.empty() )
      {
        const string samples = (m_samples.size() == 1)
                                ? (WString::tr("dsfs-sample").toUTF8() + " " + std::to_string(*begin(m_samples)))
                                : (WString::tr("dsfs-samples").toUTF8() + " " + SpecUtils::sequencesToBriefString(m_samples));
        answer += (answer.empty() ? "" : ", ") + samples;
      }
      return answer;
    }//std::string locationLabel() const
    
  };//class DrfSpecFileSample
  
  
  class DrfSpecFile : public WPanel
  {
    std::shared_ptr<SpecMeas> m_meas;
    WContainerWidget *m_sampleWidgets;
    Wt::Signal<> m_updated;
    
  public:
    DrfSpecFile( std::shared_ptr<SpecMeas> meas,
                WContainerWidget */*parent*/ = nullptr )
      : WPanel(),
        m_meas( meas ),
        m_sampleWidgets( nullptr )
    {
      addStyleClass( "DrfSpecFile" );

      const set<set<int>> sampsWithPeaks = meas->sampleNumsWithPeaks();

      auto bodyOwned = std::make_unique<WContainerWidget>();
      m_sampleWidgets = bodyOwned.get();
      setCentralWidget( std::move(bodyOwned) );

      int nsamples = 0;
      for( const set<int> &peakSamps : sampsWithPeaks )
      {
        std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaksptr = meas->peaks( peakSamps );
        if( !peaksptr || peaksptr->empty() )
          continue;

        ++nsamples;
        DrfSpecFileSample *sample = m_sampleWidgets->addNew<DrfSpecFileSample>( meas, peakSamps, nullptr );
        sample->srcInfoUpdated().connect( this, [this](){ m_updated.emit(); } );
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
    
    void setIsEffGeometryType( const DetectorPeakResponse::EffGeometryType type )
    {
      for( DrfSpecFileSample *fs : fileSamples() )
        fs->setIsEffGeometryType( type );
    }//void setIsEffGeometryType( const bool is_fixed_geom )
    
    std::shared_ptr<const SpecMeas> measurement(){ return m_meas; }
  };//class DrfSpecFile
}//namespace to implement DrfSpecFile and DrfSpecFileSample



MakeDrfWindow::MakeDrfWindow( InterSpec *viewer )
: AuxWindow( WString::tr("window-title-create-drf"),
            (AuxWindowProperties::TabletNotFullScreen
             | AuxWindowProperties::SetCloseable
             | AuxWindowProperties::DisableCollapse
             | AuxWindowProperties::EnableResize
             | AuxWindowProperties::IsModal) ),
  m_tool( nullptr )
{
  const int ww = viewer->renderedWidth();
  const int wh = viewer->renderedHeight();
  if( ww > 100 && wh > 100 )
  {
    const int width = make_drf_dialog_width( ww );
    const int height = ((wh < 420) ? wh : (19*wh)/20 );
    
    resizeWindow( width, height );
    setMinimumSize( std::min(width,720), std::min(height,480) );
  }//if( ww > 100 && wh > 100 )
    
  auto toolOwned = std::make_unique<MakeDrf>( viewer );
  m_tool = toolOwned.get();

  stretcher()->addWidget( std::move(toolOwned), 0, 0 );
  stretcher()->setContentsMargins( 0, 0, 0, 0 );

  AuxWindow::addHelpInFooter( footer(), "make-drf" );

  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close") );
  closeButton->clicked().connect( this, &AuxWindow::hide );

  WPushButton *saveAs = footer()->addNew<WPushButton>( WString::tr("md-export-btn") );
  WidgetUtils::applyButtonRole( saveAs, WidgetUtils::ButtonRole::Neutral );
  saveAs->clicked().connect( m_tool, &MakeDrf::requestStoreExport );
  m_tool->farFieldIntrinsicEfficiencyIsValid().connect( saveAs, [saveAs]( bool a1 ){ saveAs->setEnabled( a1 ); } );
  saveAs->disable();

  // No characterization data: build an approximate response by hand in the Modify DRF tool
  WPushButton *noDataBtn = footer()->addNew<WPushButton>( WString::tr("md-make-without-data-btn") );
  WidgetUtils::applyButtonRole( noDataBtn, WidgetUtils::ButtonRole::Neutral );
  noDataBtn->clicked().connect( std::bind( [viewer,this](){
    hide();
    viewer->showDrfModifyWindow( nullptr, true );
  } ) );

  show();
    
  resizeToFitOnScreen();
  centerWindow();
  rejectWhenEscapePressed( false );
  
  m_tool->m_finished.connect( this, &AuxWindow::hide );
}//MakeDrfWindow constructor
  

MakeDrf *MakeDrfWindow::tool()
{
  return m_tool;
}
  


MakeDrf::MakeDrf( InterSpec *viewer )
: WContainerWidget(),
  m_interspec( viewer ),
  m_intrinsicEfficiencyIsValid(),
  m_finished(),
  m_chart( nullptr ),
  m_files( nullptr ),
  m_sidebar( nullptr ),
  m_detSection( nullptr ),
  m_detSectionBody( nullptr ),
  m_detNotUsedTxt( nullptr ),
  m_geomMode( nullptr ),
  m_geomGuidance( nullptr ),
  m_diameterDiv( nullptr ),
  m_mcTool( nullptr ),
  m_geomSeeded( false ),
  m_phoneTabs( nullptr ),
  m_blockSourceUpdates( false ),
  m_detDiameter( nullptr ),
  m_detSetback( nullptr ),
  m_geometry( nullptr ),
  m_showFwhmPoints( nullptr ),
  m_fwhmOptionGroup( nullptr ),
  m_fwhmEqnType( nullptr ),
  m_sqrtEqnOrder( nullptr ),
  m_sqrtEqnOrderLabel( nullptr ),
  m_effEqnOrder( nullptr ),
  m_effEqnUnits( nullptr ),
  m_effOptionGroup( nullptr ),
  m_airAttenuate( nullptr ),
  m_fitNeedsDataTxt( nullptr ),
  m_pfpPanel( nullptr ),
  m_chartRangeSet( false ),
  m_statusArea( nullptr ),
  m_statusSummary( nullptr ),
  m_problemList( nullptr ),
  m_intrinsicEffAnswer( nullptr ),
  m_inputProblems{},
  m_fitProblems{},
  m_numPeaksSelected( 0 ),
  m_numPeaksUsed( 0 ),
  m_peaksResolutionType( PeakFitUtils::CoarseResolutionType::Unknown ),
  m_storeAfterGenerationId( -1 ),
  m_fwhmFitId( 0 ),
  m_effEqnFitId( 0 ),
  m_fwhmEqnChi2( -999.9 ),
  m_effEqnChi2( -999.9 ),
  m_effEqnDof( 0 ),
  m_effEqnUserTerms( -1 ),
  m_sqrtEqnUserTerms( -1 ),
  m_effLowerEnergy( 0.0f ),
  m_effUpperEnergy( 0.0f ),
  m_effBirgeScale( 1.0 ),
  m_dataDiameter( 0.0 ),
  m_peakFitDetPrefsGui( nullptr )
{
  assert( m_interspec );
  assert( MaterialDB::instance() );

  wApp->useStyleSheet( "InterSpec_resources/MakeDrf.css" );
  if( m_interspec )
    m_interspec->useMessageResourceBundle( "MakeDrf" );
  
  addStyleClass( "MakeDrf" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>("ShowTooltips", m_interspec );
  const bool isPhone = m_interspec->isPhone();

  // ===== Left: the setup sidebar - what the detector is, how to fit, peak fit preferences =====
  //  The scrolling is on an inner div rather than on the layout item: Wt force-sets the size of a
  //  layout item that has a scroll-bar, after which it no longer follows updateSidebarWidth().
  auto sidebarOwned = std::make_unique<WContainerWidget>();
  m_sidebar = sidebarOwned.get();
  m_sidebar->addStyleClass( "MakeDrfSidebar" );
  WContainerWidget *sidebarScroll = m_sidebar->addNew<WContainerWidget>();
  sidebarScroll->addStyleClass( "MakeDrfSidebarScroll" );

  // Adds a "label: widget" row to a two-column form grid
  auto add_form_row = []( WContainerWidget *grid, const WString &label, WFormWidget *input ) -> WLabel * {
    WLabel *l = grid->insertWidget( grid->count() - 1, std::make_unique<WLabel>( label ) );
    l->setBuddy( input );
    return l;
  };

  // --- Detector: just its diameter, or its full geometry ---
  m_detSection = sidebarScroll->addNew<GroupBox>( WString::tr("md-section-detector") );
  m_detSection->addStyleClass( "MakeDrfSection" );
  m_detSectionBody = m_detSection->addNew<WContainerWidget>();

  WContainerWidget *modeRow = m_detSectionBody->addNew<WContainerWidget>();
  modeRow->addStyleClass( "MakeDrfGeomModeRow" );
  m_geomMode = modeRow->addNew<SwitchCheckbox>( WString::tr("md-geom-mode-diameter"),
                                                 WString::tr("md-geom-mode-geometry") );
  m_geomMode->setChecked( false );
  m_geomMode->checked().connect( this, &MakeDrf::handleGeometryModeChanged );
  m_geomMode->unChecked().connect( this, &MakeDrf::handleGeometryModeChanged );
  HelpSystem::attachToolTipOn( modeRow, WString::tr("md-tt-geom-mode"), showToolTips );

  // Tells the user what the "Detector geometry" mode buys them; only shown in diameter mode
  m_geomGuidance = m_detSectionBody->addNew<WText>( WString::tr("md-geom-guidance") );
  m_geomGuidance->addStyleClass( "MakeDrfGeomGuidance" );
  m_geomGuidance->setInline( false );

  m_diameterDiv = m_detSectionBody->addNew<WContainerWidget>();
  m_diameterDiv->addStyleClass( "MakeDrfFormGrid MakeDrfDiameterDiv" );

  auto distValidator = std::make_shared<WRegExpValidator>( PhysicalUnits::sm_distanceUnitOptionalRegex );
  distValidator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  distValidator->setInvalidNoMatchText( WString::tr("msd-tt-dist-format") );

  m_detDiameter = m_diameterDiv->addNew<WLineEdit>();
  add_form_row( m_diameterDiv, WString::tr("md-det-diam"), m_detDiameter );
  m_detDiameter->setValidator( distValidator );
  m_detDiameter->setText( "2.54 cm" );
  m_detDiameter->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  m_detDiameter->enterPressed().connect( this, &MakeDrf::handleSourcesUpdates );
  m_detDiameter->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_detDiameter->setAttributeValue( "autocorrect", "off" );
  m_detDiameter->setAttributeValue( "spellcheck", "off" );
#endif

  m_detSetback = m_diameterDiv->addNew<WLineEdit>();
  add_form_row( m_diameterDiv, WString::tr("md-det-setback"), m_detSetback );
  m_detSetback->setValidator( distValidator );
  m_detSetback->setPlaceholderText( WString::tr("md-det-setback-ph") );
  m_detSetback->changed().connect( this, &MakeDrf::handleSourcesUpdates );
  m_detSetback->enterPressed().connect( this, &MakeDrf::handleSourcesUpdates );
  m_detSetback->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_detSetback->setAttributeValue( "autocorrect", "off" );
  m_detSetback->setAttributeValue( "spellcheck", "off" );
#endif
  HelpSystem::attachToolTipOn( m_diameterDiv, WString::tr("md-tt-det-diam"), showToolTips );

  // Full geometry + location support: the same tool the Modify DRF dialog embeds.  Defaults to the
  //  instant measured-curve transfer (no Monte Carlo); the user can pick an MC method and generate.
  auto mcToolOwned = std::make_unique<MakeMcResponseForDrf>( m_interspec, nullptr );
  m_mcTool = mcToolOwned.get();
  m_mcTool->addStyleClass( "MakeDrfMcTool" );
  m_mcTool->setChartHidden( true );
  m_mcTool->setMethod( MakeMcResponseForDrf::Method::CurveTransfer );
  m_mcTool->geometryInput()->setCompact( true );
  m_mcTool->geometryInput()->setDiagramBelow( true );  //the sidebar has no room beside the form
  m_mcTool->setAutoBuildTransfer( false );  //assembleDrf builds it from the fit
  m_mcTool->setImportGeometryAlwaysOffered( true );  //the form starts as a guess; nothing to lose
  m_mcTool->setAnchorHidden( true );  //assembleDrf anchors the transfer itself
  m_mcTool->geometryChanged().connect( this, &MakeDrf::handleSourcesUpdates );
  m_mcTool->responseGenerated().connect( this, &MakeDrf::setGeneratedMcResponse );
  // The tool drops its result whenever the geometry, method or an option changes; a response held
  //  past that describes some other detector, so let it go too.
  // A run the user asked to store after that ended without a response - cancelled, failed, or
  //  superseded (a method change or restored state bumps the generation id) - is said so, rather
  //  than leave them waiting for a dialog that wont come.  A successful run has a result when
  //  validationChanged(true) comes, and is disarmed by setGeneratedMcResponse before userChanged().
  auto check_store_run = [this](){
    if( m_storeAfterGenerationId < 0 )
      return;
    const bool superseded = (m_storeAfterGenerationId != m_mcTool->generationId());
    if( !superseded && (m_mcTool->isGenerating() || m_mcTool->hasResult()) )
      return;
    m_storeAfterGenerationId = -1;
    passMessage( WString::tr("md-mc-not-stored"), WarningWidget::WarningMsgMedium );
  };//check_store_run

  m_mcTool->validationChanged().connect( this, [this,check_store_run]( const bool valid ){
    if( !valid )
      m_generatedMcResponse.reset();
    check_store_run();
  } );
  m_mcTool->userChanged().connect( this, check_store_run );
  m_detSectionBody->addWidget( std::move(mcToolOwned) );
  m_mcTool->hide();  //diameter mode to start

  // A detector file's geometry can be imported from "Diameter only" mode too; doing so switches to
  //  the geometry form, since a full geometry is the point.
  WPushButton *importGeomBtn = m_diameterDiv->addNew<WPushButton>( WString::tr("mmr-import-geom-btn") );
  importGeomBtn->addStyleClass( "LinkBtn MakeDrfImportGeom" );
  importGeomBtn->setToolTip( WString::tr("mmr-tt-import-geom") );
  importGeomBtn->clicked().connect( m_mcTool, &MakeMcResponseForDrf::showImportGeometryDialog );
  m_mcTool->geometryImported().connect( this, [this](){
    if( m_geomMode->isChecked() )
      return;
    m_geomSeeded = true;  //the imported geometry, not one guessed from the diameter
    m_geomMode->setChecked( true );
    handleGeometryModeChanged();
  } );

  // Fixed-geometry efficiencies dont use the detector description; say so rather than vanish
  m_detNotUsedTxt = m_detSection->addNew<WText>( WString::tr("md-det-not-used-fixed") );
  m_detNotUsedTxt->addStyleClass( "FainterTxt MakeDrfSectionNote" );
  m_detNotUsedTxt->setInline( false );
  m_detNotUsedTxt->hide();

  // --- Efficiency fit ---
  m_effOptionGroup = sidebarScroll->addNew<GroupBox>( WString::tr("md-section-eff") );
  m_effOptionGroup->addStyleClass( "MakeDrfSection" );
  WContainerWidget *effGrid = m_effOptionGroup->addNew<WContainerWidget>();
  effGrid->addStyleClass( "MakeDrfFormGrid" );

  m_effEqnOrder = effGrid->addNew<WComboBox>();
  add_form_row( effGrid, WString::tr("md-lbl-eqn-terms"), m_effEqnOrder );
  m_effEqnOrder->setNoSelectionEnabled( true );
  m_effEqnOrder->changed().connect( this, &MakeDrf::handleEffEqnOrderChange );

  m_effEqnUnits = effGrid->addNew<WComboBox>();
  add_form_row( effGrid, WString::tr("md-lbl-eqn-units"), m_effEqnUnits );
  m_effEqnUnits->addItem( WString::tr("md-eqn-in-kev") );
  m_effEqnUnits->addItem( WString::tr("md-eqn-in-mev") );
  m_effEqnUnits->setCurrentIndex( 1 );
  m_effEqnUnits->changed().connect( this, &MakeDrf::handleSourcesUpdates );

  m_geometry = effGrid->addNew<WComboBox>();
  add_form_row( effGrid, WString::tr("md-lbl-eff-type"), m_geometry );
  m_geometry->addItem( WString::tr("md-geom-intrinsic-eff") );
  m_geometry->addItem( WString::tr("md-geom-fixed-geom") );
  m_geometry->addItem( WString::tr("md-geom-fixed-cm2") );
  m_geometry->addItem( WString::tr("md-geom-fixed-m2") );
  m_geometry->addItem( WString::tr("md-geom-fixed-gram") );
  m_geometry->setCurrentIndex( 0 );
  HelpSystem::attachToolTipOn( m_geometry, WString::tr("md-tt-eff-interpret"), showToolTips,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_geometry->activated().connect( this, &MakeDrf::handleFixedGeometryChanged );

  m_airAttenuate = m_effOptionGroup->addNew<WCheckBox>( WString::tr("md-atten-for-air-cb") );
  m_airAttenuate->setChecked( true );
  m_airAttenuate->setInline( false );
  m_airAttenuate->addStyleClass( "AirAttenCb" );
  m_airAttenuate->checked().connect( this, &MakeDrf::handleSourcesUpdates );
  m_airAttenuate->unChecked().connect( this, &MakeDrf::handleSourcesUpdates );
  HelpSystem::attachToolTipOn( m_airAttenuate, WString::tr("md-tt-atten-for-air"), showToolTips,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );

  m_fitNeedsDataTxt = m_effOptionGroup->addNew<WText>( WString::tr("md-fit-needs-data") );
  m_fitNeedsDataTxt->addStyleClass( "FainterTxt MakeDrfSectionNote" );
  m_fitNeedsDataTxt->setInline( false );

  // --- FWHM fit ---
  m_fwhmOptionGroup = sidebarScroll->addNew<GroupBox>( WString::tr("md-section-fwhm") );
  m_fwhmOptionGroup->addStyleClass( "MakeDrfSection" );
  WContainerWidget *fwhmGrid = m_fwhmOptionGroup->addNew<WContainerWidget>();
  fwhmGrid->addStyleClass( "MakeDrfFormGrid" );

  m_fwhmEqnType = fwhmGrid->addNew<WComboBox>();
  add_form_row( fwhmGrid, WString::tr("md-lbl-fwhm-form"), m_fwhmEqnType );
  m_fwhmEqnType->setNoSelectionEnabled( true );
  for( auto i = DetectorPeakResponse::ResolutionFnctForm(0);
      i < DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm;
      i = DetectorPeakResponse::ResolutionFnctForm(i+1) )
  {
    switch( i )
    {
      case DetectorPeakResponse::kGadrasResolutionFcn:
        m_fwhmEqnType->addItem( WString::tr("md-fwhm-gadras") );
        break;
        
      case DetectorPeakResponse::kSqrtEnergyPlusInverse:
        m_fwhmEqnType->addItem( "sqrt(A + B*E + C/E)" );
        break;
        
      case DetectorPeakResponse::kConstantPlusSqrtEnergy:
        m_fwhmEqnType->addItem( "A + B*sqrt(E)" );
        break;
        
      case DetectorPeakResponse::kSqrtPolynomial:
        m_fwhmEqnType->addItem( WString::tr("md-sqrt-series") );
        break;
        
      case DetectorPeakResponse::kNumResolutionFnctForm:
        assert( 0 );
        break;
    }//switch( i )
  }//for( loop over DetectorPeakResponse::ResolutionFnctForm )
  m_fwhmEqnType->setCurrentIndex( 0 );
  m_fwhmEqnType->changed().connect( this, &MakeDrf::handleFwhmTypeChanged );

  m_sqrtEqnOrder = fwhmGrid->addNew<WComboBox>();
  m_sqrtEqnOrderLabel = add_form_row( fwhmGrid, WString::tr("md-lbl-eqn-terms"), m_sqrtEqnOrder );
  m_sqrtEqnOrder->setNoSelectionEnabled( true );
  m_sqrtEqnOrder->changed().connect( this, &MakeDrf::handleSqrtEqnOrderChange );
  m_sqrtEqnOrder->hide();
  m_sqrtEqnOrderLabel->hide();

  m_showFwhmPoints = m_fwhmOptionGroup->addNew<WCheckBox>( WString::tr("md-show-fwhm-markers") );
  m_showFwhmPoints->setChecked( true );
  m_showFwhmPoints->setInline( false );
  m_showFwhmPoints->checked().connect( this, &MakeDrf::handleShowFwhmPointsToggled );
  m_showFwhmPoints->unChecked().connect( this, &MakeDrf::handleShowFwhmPointsToggled );

  // --- Peak fit preferences (optional; saved into the DRF).  Collapsed to a one-line summary,
  //  since they are rarely changed, and expand to a tall block when a skew is chosen.
  m_pfpPanel = sidebarScroll->addNew<WPanel>();
  m_pfpPanel->addStyleClass( "MakeDrfPfpPanel" );
  m_peakFitDetPrefsGui = m_pfpPanel->setCentralWidget( std::make_unique<PeakFitDetPrefsGui>( m_interspec, false ) );
  m_pfpPanel->setCollapsible( true );
  m_pfpPanel->setCollapsed( true );
  m_interspec->peakFitDetPrefsChanged().connect( this, &MakeDrf::handlePeakFitPrefsChanged );

  // ===== Right: the chart, the status (problems and the fitted equation), the peaks and sources =====
  auto chartOwned = std::make_unique<DrfChart>();
  m_chart = chartOwned.get();
  m_chart->addStyleClass( "MakeDrfChart" );
  m_chart->setKeepZoomOnUpdate( true );
  m_chart->setMinimumSize( WLength::Auto, WLength(160, WLength::Unit::Pixel) );

  auto statusOwned = std::make_unique<WContainerWidget>();
  m_statusArea = statusOwned.get();
  m_statusArea->addStyleClass( "MakeDrfStatus" );
  m_statusSummary = m_statusArea->addNew<WText>( "", Wt::TextFormat::Plain );
  m_statusSummary->addStyleClass( "MakeDrfStatusSummary" );
  m_statusSummary->setInline( false );
  m_problemList = m_statusArea->addNew<WContainerWidget>();
  m_problemList->addStyleClass( "MakeDrfProblems" );
  m_problemList->hide();
  m_intrinsicEffAnswer = m_statusArea->addNew<WText>( WString::fromUTF8(""), Wt::TextFormat::XHTML );
  m_intrinsicEffAnswer->addStyleClass( "MakeDrfEqnAnswer" );
  m_intrinsicEffAnswer->setInline( false );
  m_intrinsicEffAnswer->hide();

  auto fileHolderOwned = std::make_unique<WContainerWidget>();
  WContainerWidget *fileHolder = fileHolderOwned.get();
  fileHolder->setOverflow( Overflow::Auto, Wt::Orientation::Vertical );
  fileHolder->addStyleClass( "MakeDrfFiles" );
  m_files = fileHolder->addNew<WContainerWidget>();

  WGridLayout *layout = setLayout( std::make_unique<WGridLayout>() );
  layout->setContentsMargins( 0, 0, 0, 0 );
  
  if( isPhone )
  {
    // One thing at a time: setup, then the chart, then the peaks and their sources.  The tabs load
    //  eagerly so a field the status list points at is always there to focus.
    m_phoneTabs = layout->addWidget( std::make_unique<WTabWidget>(), 0, 0 );
    m_phoneTabs->addTab( std::move(sidebarOwned), WString::tr("md-tab-setup"), ContentLoading::Eager );

    auto chartTabOwned = std::make_unique<WContainerWidget>();
    WGridLayout *chartLayout = chartTabOwned->setLayout( std::make_unique<WGridLayout>() );
    chartLayout->setContentsMargins( 0, 0, 0, 0 );
    chartLayout->addWidget( std::move(chartOwned), 0, 0 );
    chartLayout->addWidget( std::move(statusOwned), 1, 0 );
    chartLayout->setRowStretch( 0, 1 );
    m_phoneTabs->addTab( std::move(chartTabOwned), WString::tr("md-tab-chart"), ContentLoading::Eager );

    m_phoneTabs->addTab( std::move(fileHolderOwned), WString::tr("md-mi-peaks-phone"), ContentLoading::Eager );
  }else
  {
    layout->setHorizontalSpacing( 8 );
    layout->addWidget( std::move(sidebarOwned), 0, 0 );

    auto rightOwned = std::make_unique<WContainerWidget>();
    WContainerWidget *rightCol = rightOwned.get();
    rightCol->addStyleClass( "MakeDrfRight" );
    WGridLayout *rightLayout = rightCol->setLayout( std::make_unique<WGridLayout>() );
    rightLayout->setContentsMargins( 0, 0, 0, 0 );
    rightLayout->setVerticalSpacing( 8 );
    rightLayout->addWidget( std::move(chartOwned), 0, 0 );
    rightLayout->addWidget( std::move(statusOwned), 1, 0 );
    rightLayout->addWidget( std::move(fileHolderOwned), 2, 0 );
    rightLayout->setRowStretch( 2, 1 );
    fileHolder->setMinimumSize( WLength::Auto, WLength(120, WLength::Unit::Pixel) );

    // The divider between the chart and what is below it can be dragged.  Set once: setting it
    //  again would throw away the height the user dragged to.
    const int wh = m_interspec->renderedHeight();
    const double chart_height = (wh > 100) ? std::clamp( 0.45*(0.95*wh - 110.0), 200.0, 420.0 ) : 300.0;
    rightLayout->setRowResizable( 0, true, WLength(chart_height, WLength::Unit::Pixel) );
    WidgetUtils::compensateRowHandleLeak( rightCol );

    layout->addWidget( std::move(rightOwned), 0, 1 );
    layout->setColumnStretch( 1, 1 );
    layout->setRowStretch( 0, 1 );

    updateSidebarWidth();
  }//if( isPhone ) / else
  
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
    DrfSpecFile *fileWidget = m_files->addNew<DrfSpecFile>( meas );
    fileWidget->updated().connect( this, &MakeDrf::handleSourcesUpdates );
    added.push_back( fileWidget );
    

    for( DrfSpecFileSample *fs : fileWidget->fileSamples() )
    {
      for( DrfPeak *p : fs->peaks() )
        p->m_peakPreviewShow.connect( this, [this]( DrfPeak *a1 ){ peakPreviewShown( a1 ); } );
    }
    
  }//for( loop over opened files )
  
  if( added.size() == 1 )
  {
    added[0]->addStyleClass( "SingleFile" );
    if( !added[0]->title().empty() )
      added[0]->setTitle( "" );
  }
  
  handlePeakFitPrefsChanged();
  updateStatusDisplay();

  //If we directly call handleSourcesUpdates() now, getting the activity
  //  uncertainty will throw an exception because they wont validate.... whatever.
  //  Re-resolve by DOM id rather than capturing `this`: we live in a transient MakeDrfWindow and
  //  `WServer::post` only guarantees the session is alive (same pattern as fitFwhmEqn() below).
  const WidgetUtils::WidgetHandle self( this );
  WServer::instance()->post( wApp->sessionId(), [self](){
    MakeDrf * const drf = self.resolve_as<MakeDrf>();
    if( drf )
      drf->handleSourcesUpdates();
  } );
}//MakeDrf( constructor )


MakeDrf::~MakeDrf()
{
  if( m_saveAsWindow )
    AuxWindow::deleteAuxWindow( m_saveAsWindow.get() );
}//~MakeDrf()


void MakeDrf::startSaveAs()
{
  Wt::WFlags<AuxWindowProperties> windowprop = (AuxWindowProperties::IsModal
                                               | AuxWindowProperties::PhoneNotFullScreen
                                               | AuxWindowProperties::SetCloseable
                                               | AuxWindowProperties::DisableCollapse);
  if( m_effEqnCoefs.empty() )
  {
    AuxWindow *w = AuxWindow::make( WString::tr("md-err-window-title"), windowprop );
    WText *t = w->contents()->addNew<WText>( WString::tr("md-drf-not-valid"), Wt::TextFormat::XHTML );
    t->setInline( false );
    WPushButton *b = w->footer()->addNew<WPushButton>( WString::tr("Close") );
    WidgetUtils::applyButtonRole( b, WidgetUtils::ButtonRole::Dismiss );
    b->clicked().connect( w, &AuxWindow::hide );
    w->finished().connect( w, [w](Wt::DialogCode){ AuxWindow::deleteAuxWindow( w ); } );
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
  
  if( m_saveAsWindow )
    AuxWindow::deleteAuxWindow( m_saveAsWindow.get() );

  AuxWindow *w = AuxWindow::make( WString::tr("md-export-window-title"), windowprop );
  m_saveAsWindow = w;

  WTable *table = w->contents()->addNew<WTable>();
  table->addStyleClass( "MakeDrfSaveAsTable" );

  WTableCell *cell = table->elementAt(0, 0);
  WLabel *label = cell->addNew<WLabel>( WString::tr("Name") );

  cell = table->elementAt(0, 1);
  WLineEdit *name = cell->addNew<WLineEdit>();
  label->setBuddy( name );
  name->setTextSize( 32 );
  name->setMaxLength( 255 );
  name->setAutoComplete( false );

  name->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  name->setAttributeValue( "autocorrect", "off" );
  name->setAttributeValue( "spellcheck", "off" );
#endif

  //const char *valid_file_name_regex = "(^(?!\\.)(?!com[0-9]$)(?!con$)(?!lpt[0-9]$)(?!nul$)(?!prn$)[^\\|*\\?\\\\:<>/$\"]*[^\\.\\|*\\?\\\\:<>/$\"]+$";
  const char *valid_file_name_regex = R"MyRegexDelim(^[^\\\/\:\*\?\"\'\,\;\<\>\|]+$)MyRegexDelim";
  auto validatorOwned = std::make_shared<WRegExpValidator>( valid_file_name_regex );
  WRegExpValidator *validator = validatorOwned.get();
  name->setValidator( validatorOwned );

  cell = table->elementAt(0, 2);
  WImage *help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
  help->addStyleClass( "MakeDrfSaveHelp" );
  HelpSystem::attachToolTipOn( help, WString::tr("md-tt-name"), true,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );

  cell = table->elementAt(1, 0);
  label = cell->addNew<WLabel>( WString::tr("Description") );
  cell = table->elementAt(1, 1);
  cell->setRowSpan( 2 );
  WTextArea *description = cell->addNew<WTextArea>();
  label->setBuddy( description );
  description->setWidth( WLength(97,WLength::Unit::Percentage) );

  cell = table->elementAt(1, 2);
  help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
  help->addStyleClass( "MakeDrfSaveHelp" );
  HelpSystem::attachToolTipOn( help, WString::tr("md-tt-desc"), true,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  
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
    
      def_for_serial_cb = cell->addNew<WCheckBox>( WString::tr("md-make-default-for-serial-cb").arg(serial_number) );
      cell = table->elementAt(currentRow, 2);
      help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
      help->addStyleClass( "MakeDrfSaveHelp" );
      HelpSystem::attachToolTipOn( help, WString::tr("md-tt-make-default-serial"),
                                  true, HelpSystem::ToolTipPrefOverride::AlwaysShow );
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
      def_for_model_cb = cell->addNew<WCheckBox>( WString::tr("md-make-default-for-model-cb").arg(model) );
      cell = table->elementAt(currentRow, 2);
      help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
      help->addStyleClass( "MakeDrfSaveHelp" );
      HelpSystem::attachToolTipOn( help, WString::tr("md-tt-make-default-model"),
                                  true, HelpSystem::ToolTipPrefOverride::AlwaysShow );
    }//if( !model.empty() )
  }//if( representative_meas )
  
  
  currentRow = table->rowCount();
  cell = table->elementAt(currentRow, 0);
  cell->setColumnSpan( 2 );
  auto n42ResourceOwned = std::make_shared<CalFileDownloadResource>( false, this );
  CalFileDownloadResource *n42Resource = n42ResourceOwned.get();
  {
    Wt::WLink lnk( n42ResourceOwned );
    lnk.setTarget( Wt::LinkTarget::NewWindow );
    WAnchor *n42anchor = cell->addNew<WAnchor>( lnk, WString::tr("md-export-n42-anchor") );
    (void)n42anchor;
  }
  cell = table->elementAt(currentRow, 2);
  help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
  help->addStyleClass( "MakeDrfSaveHelp" );
  HelpSystem::attachToolTipOn( help, WString::tr("md-tt-export-n42"), true,
                  HelpSystem::ToolTipPrefOverride::AlwaysShow );

  currentRow = table->rowCount();
  cell = table->elementAt(currentRow, 0);
  cell->setColumnSpan( 2 );
  auto xmlResourceOwned = std::make_shared<DrfXmlDownloadResource>( this );
  DrfXmlDownloadResource *xmlResource = xmlResourceOwned.get();
  {
    Wt::WLink lnk( xmlResourceOwned );
    lnk.setTarget( Wt::LinkTarget::NewWindow );
    cell->addNew<WAnchor>( lnk, WString::tr("md-export-drf-xml") );
  }
  cell = table->elementAt(currentRow, 2);
  help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
  help->addStyleClass( "MakeDrfSaveHelp" );
  HelpSystem::attachToolTipOn( help, WString::tr("md-tt-export-drf-xml"), true,
                  HelpSystem::ToolTipPrefOverride::AlwaysShow );

  currentRow = table->rowCount();
  cell = table->elementAt(currentRow, 0);
  cell->setColumnSpan( 2 );
  auto csvResourceOwned = std::make_shared<CsvDrfDownloadResource>( this );
  CsvDrfDownloadResource *csvResource = csvResourceOwned.get();
  {
    Wt::WLink lnk( csvResourceOwned );
    lnk.setTarget( Wt::LinkTarget::NewWindow );
    WAnchor *csvanchor = cell->addNew<WAnchor>( lnk, WString::tr("md-export-as-csv") );
    (void)csvanchor;
  }

  cell = table->elementAt(currentRow, 2);
  help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
  help->addStyleClass( "MakeDrfSaveHelp" );
  HelpSystem::attachToolTipOn( help, WString::tr("md-tt-export-as-csv"), true,
                  HelpSystem::ToolTipPrefOverride::AlwaysShow );


  currentRow = table->rowCount();
  cell = table->elementAt(currentRow, 0);
  cell->setColumnSpan( 2 );
  auto refSheetResourceOwned = std::make_shared<RefSheetDownloadResource>( this );
  RefSheetDownloadResource *refSheetResource = refSheetResourceOwned.get();
  {
    Wt::WLink lnk( refSheetResourceOwned );
    lnk.setTarget( Wt::LinkTarget::NewWindow );
    WAnchor *refSheetAnchor = cell->addNew<WAnchor>( lnk, WString::tr("md-export-quick-ref") );
    (void)refSheetAnchor;
  }
  cell = table->elementAt(currentRow, 2);
  help = cell->addNew<WImage>( Wt::WLink("InterSpec_resources/images/help_mobile.svg") );
  help->addStyleClass( "MakeDrfSaveHelp" );
  HelpSystem::attachToolTipOn( help, WString::tr("md-tt-export-quick-ref"), true,
                  HelpSystem::ToolTipPrefOverride::AlwaysShow );


  auto updateName = [name,csvResource,n42Resource,refSheetResource,xmlResource](){
    const string drfname = (name->validate() == Wt::ValidationState::Valid) ? name->text().toUTF8()
                                                                            : string();
    n42Resource->setSuggestFileName( drfname );
    csvResource->setSuggestFileName( drfname );
    refSheetResource->setSuggestFileName( drfname );
    xmlResource->setSuggestFileName( drfname );
  };
  
  name->changed().connect( this, updateName );
  name->enterPressed().connect( this, updateName );
  name->blurred().connect( this, updateName );


  auto updateDesc = [description,n42Resource,csvResource,refSheetResource,xmlResource](){
    n42Resource->setDescription( description->text().toUTF8() );
    csvResource->setDescription( description->text().toUTF8() );
    refSheetResource->setDescription( description->text().toUTF8() );
    xmlResource->setDescription( description->text().toUTF8() );
  };

  description->changed().connect( this, updateDesc );
  description->enterPressed().connect( this, updateDesc );
  description->blurred().connect( this, updateDesc );

  
  
  WPushButton *b = w->addCloseButtonToFooter( WString::tr("Cancel") );
  b->clicked().connect( w, &AuxWindow::hide );
  
  // TODO: Need to display messsage to user about how it is saved, how to access it again, and figure out why menus go away for people (maybe only on Windows?) - also should consider not closing AuxWindow - also fix assay date being invalid error (seems to cause uncaught exception).
  auto doSave = [this, w, validator, name, description,
                 def_for_serial_cb, def_for_model_cb,
                 representative_meas ](){
    auto state = validator->validate(name->text()).state();
    if( name->text().empty() || state!=WValidator::State::Valid )
    {
      passMessage( WString::tr("md-det-name-not-valid"), 3 );
      return;
    }
    
    if( m_effEqnCoefs.empty() )
    {
      passMessage( WString::tr("md-eqn-not-valid"), 3 );
      return;
    }
    
    string drfname = name->text().toUTF8();
    string drfdescrip = description->text().toUTF8();
    SpecUtils::trim( drfname );
    SpecUtils::trim( drfdescrip );
    
    std::shared_ptr<DetectorPeakResponse> drf;
    try
    {
      drf = assembleDrf( drfname, drfdescrip );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("md-err-creating-drf").arg(e.what()), 3 );
      return;
    }
    
    try
    {
      auto sql = std::make_shared<DataBaseUtils::DbSession>( *m_interspec->sql() );
      
      //negladgeable chance we have another detector with same hash, so we'll
      // skip checking for it.
      
      DataBaseUtils::DbTransaction transaction( *sql );
      //Create a separate DetectorPeakResponse because shared_ptr and dbo::ptr don't work well together
      auto tempDetectorOwned = std::make_unique<DetectorPeakResponse>( *drf );
      tempDetectorOwned->m_user = static_cast<int>( m_interspec->user().id() );
      auto newDbDet = sql->session()->add( std::move(tempDetectorOwned) );
      
      transaction.commit();
    }catch( std::exception &e )
    {
      passMessage( WString::tr("md-err-saving-to-db").arg(e.what()), WarningWidget::WarningMsgHigh );
      return;
    }//try / catch
    
    m_interspec->detectorChanged().emit( drf );
    
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    const Wt::Dbo::ptr<InterSpecUser> &user = m_interspec->user();
    
    // On this (the session) thread: copying the `Dbo::ptr` user is not thread safe.
    if( def_for_serial_cb && def_for_serial_cb->isChecked() && representative_meas )
      DrfSelect::setUserPrefferedDetector( drf, sql, user, UseDrfPref::UseDrfType::UseDetectorSerialNumber,
                                           representative_meas );

    if( def_for_model_cb && def_for_model_cb->isChecked() && representative_meas )
      DrfSelect::setUserPrefferedDetector( drf, sql, user, UseDrfPref::UseDrfType::UseDetectorModelName,
                                           representative_meas );
    
    w->hide();
    
    //m_finished.emit();
    passMessage( WString::tr("md-saved-to-db").arg(drfname), 0 );
  };//auto doSave
  
  WPushButton *save = w->addCloseButtonToFooter( WString::tr("Save"), WidgetUtils::ButtonRole::Affirm );
  save->clicked().connect( this, doSave );
  save->disable();
  
  auto updateSaveBtn = [name,save,validator](){
    bool enable = false;
    switch( validator->validate(name->text()).state() )
    {
      case WValidator::State::Invalid:      enable = false; break;
      case WValidator::State::InvalidEmpty: enable = false; break;
      case WValidator::State::Valid:        enable = !name->text().empty(); break;
    }
    save->setEnabled( enable );
  };
  
  name->validated().connect( this, updateSaveBtn );  //doesnt seem to ever get called...
  //name->changed().connect( this, updateSaveBtn );
  //name->enterPressed().connect( this, updateSaveBtn );
  name->keyPressed().connect( this, updateSaveBtn ); //Uhg, I guess we'll just use this until I figure out whats wrong with above so save button will update in real time
  
  
  w->finished().connect( w, [w](Wt::DialogCode){ AuxWindow::deleteAuxWindow( w ); } );
  w->rejectWhenEscapePressed();
  w->show();
  w->centerWindow();
}//void startSaveAs();


Wt::Signal<bool> &MakeDrf::farFieldIntrinsicEfficiencyIsValid()
{
  return m_intrinsicEfficiencyIsValid;
}


bool MakeDrf::isIntrinsicEfficiencyValid()
{
  return !m_effEqnCoefs.empty();
}//bool MakeDrf::isDrfValid()


void MakeDrf::setGeneratedMcResponse( std::shared_ptr<ceelo::DetectorResponse> response )
{
  m_generatedMcResponse = response;

  // Only the run "Generate, then store" (or "Store when done") armed: a run the user started
  //  themselves, or the instant transfer rebuilding itself, must not pop up the dialog.
  const bool armed = (m_storeAfterGenerationId >= 0) && m_mcTool
                     && (m_storeAfterGenerationId == m_mcTool->generationId());
  m_storeAfterGenerationId = -1;
  if( armed && response )
    startSaveAs();
}//setGeneratedMcResponse(...)


void MakeDrf::requestStoreExport()
{
  m_storeAfterGenerationId = -1;  //this request replaces any earlier one

  // Only a Monte-Carlo method has to be generated first; the instant transfer from the measured
  //  curve (the default) is built as the DRF is assembled.
  const bool full_geom = m_mcTool && (m_geometry->currentIndex() == 0) && m_geomMode->isChecked();
  const bool mc_method = full_geom
                && (m_mcTool->selectedMethod() != MakeMcResponseForDrf::Method::CurveTransfer);
  const bool have_mc = m_generatedMcResponse
                && (m_generatedMcResponse->provenance.method != ceelo::ProductionMethod::CurveTransfer);

  if( !mc_method || have_mc || m_effEqnCoefs.empty() )
  {
    startSaveAs();
    return;
  }

  const bool running = m_mcTool->isGenerating();
  const bool can_generate = m_mcTool->generationReady();

  WString body;
  if( running )
    body = WString::tr("md-mc-running-body");
  else if( can_generate )
    body = WString::tr("md-mc-none-body");
  else  //the reason can quote typed text, so is escaped for the XHTML
    body = WString::tr("md-mc-cant-body")
             .arg( Wt::Utils::htmlEncode( m_mcTool->geometryProblem() ) );

  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>(
                  WString::tr( running ? "md-mc-running-title" : "md-mc-none-title" ), body );

  if( running )
  {
    WPushButton *wait = dialog->addButton( WString::tr("md-mc-store-when-done"),
                                           WidgetUtils::ButtonRole::Affirm );
    wait->clicked().connect( this, [this](){
      // Re-checked: the run may have finished while the dialog was up
      if( m_mcTool->isGenerating() )
      {
        m_storeAfterGenerationId = m_mcTool->generationId();
        passMessage( WString::tr("md-mc-will-store"), WarningWidget::WarningMsgInfo );
      }else
      {
        requestStoreExport();
      }
    } );
  }else if( can_generate )
  {
    WPushButton *generate = dialog->addButton( WString::tr("md-mc-generate-then-store"),
                                               WidgetUtils::ButtonRole::Affirm );
    generate->clicked().connect( this, [this](){
      // Armed to THIS run only: one that never starts must not leave a later one armed
      const bool started = m_mcTool->startGeneration();
      m_storeAfterGenerationId = started ? m_mcTool->generationId() : -1;
      if( started )
      {
        passMessage( WString::tr("md-mc-will-store"), WarningWidget::WarningMsgInfo );
      }else
      {
        passMessage( WString::tr("md-mc-not-started")
                       .arg( Wt::Utils::htmlEncode( m_mcTool->geometryProblem() ) ),
                     WarningWidget::WarningMsgHigh );
      }
    } );
  }//if( running ) / else if( can_generate )

  WPushButton *without = dialog->addButton( WString::tr("md-mc-store-without"),
                                            WidgetUtils::ButtonRole::Neutral );
  without->clicked().connect( this, &MakeDrf::startSaveAs );
  dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );
}//void requestStoreExport()


void MakeDrf::handleSourcesUpdates()
{
  // While the source widgets are being updated in bulk, they are not all in a consistent state
  //  (e.g. some still have their distance input hidden, while m_geometry already says far-field),
  //  so dont look at them until whoever is changing them is done and calls us itself.
  if( m_blockSourceUpdates )
    return;

  size_t numchan = 0;
  vector< std::shared_ptr<const PeakDef> > peaks, fwhm_peaks;  //the peaks used, and the ones the FWHM fit uses

  // Raw per-peak points (absolute efficiency at each source's own distance) and the sources they
  //  came from; per-source keys keep each source's certificate/distance error a common-mode block.
  vector<MeasuredEffPoint> measured_points;
  vector<pair<MakeDrfSrcDef *,MakeDrfCalc::CalSource>> used_sources;  //in the order first used, as read

  const bool is_fixed_geometry = (m_geometry->currentIndex() != 0);

  // What is wrong with the inputs, for the list under the chart: the detector first, then sources
  vector<StatusItem> problems = validateDetectorInputs();
  const bool detector_inputs_ok = std::none_of( begin(problems), end(problems),
                                                []( const StatusItem &i ){ return i.isError; } );

  // How the detector is described (a diameter, or its full geometry)
  string geomErrMsg;
  MakeDrfCalc::GeometryChoice geom_choice;
  try
  {
    geom_choice = geometryChoice();
  }catch( std::exception &e )
  {
    geomErrMsg = e.what();
    if( detector_inputs_ok )  //else validateDetectorInputs() has already said what is wrong
      problems.push_back( makeStatusItem( WString::fromUTF8(geomErrMsg), nullptr, true ) );
  }
  const double diameter = geom_choice.diameter;
  
  const bool doAirAtten = (!is_fixed_geometry && m_airAttenuate->isChecked());
  
  // The background peaks, with the background spectrum each is from; and how many source spectra
  //  there are - messages name the file and sample a source is in, unless there is only the one.
  vector<MakeDrfCalc::BackgroundPeak> backgroundpeaks;
  size_t num_source_samples = 0, num_background_samples = 0;
  for( auto w : m_files->children() )
  {
    auto fileWidget = dynamic_cast<DrfSpecFile *>( w );
    if( !fileWidget )
      continue;
    
    numchan = std::max( numchan, fileWidget->measurement()->num_gamma_channels() );
    
    for( DrfSpecFileSample *sample : fileWidget->fileSamples() )
    {
      // A peak shows it is background subtracted only while it is (set below)
      for( DrfPeak *peak : sample->peaks() )
        peak->setBackgroundBeingSubtractedInfo( false, 0.0 );

      if( !sample->isBackground() )
      {
        num_source_samples += 1;
        continue;
      }
      
      for( DrfPeak *peak : sample->peaks() )
      {
        if( !peak->useForEffFit() )
          continue;

        MakeDrfCalc::BackgroundPeak back;
        back.peak = peak->m_peak;
        back.liveTime = sample->liveTime();
        back.spectrum = num_background_samples;
        backgroundpeaks.push_back( back );
      }

      num_background_samples += 1;
    }//for( DrfSpecFileSample *sample : sampleWidgets )
  }//for( auto w : m_files->children()  )

  // A status item for a source: "<b>Na24</b> - file.n42, Sample 2: the problem"
  auto source_item = [num_source_samples]( MakeDrfSrcDef *src, DrfSpecFileSample *sample,
                                           const WString &msg, WFormWidget *field,
                                           const bool is_error ) -> StatusItem {
    const SandiaDecay::Nuclide * const nuc = src ? src->nuclide() : nullptr;
    const string nucname = nuc ? nuc->symbol : WString::tr("msd-non-specified-nuc").toUTF8();
    const string location = (num_source_samples > 1) ? sample->locationLabel() : string();
    const char * const key = location.empty() ? "md-prob-src-noloc" : "md-prob-src";

    StatusItem item;
    item.html = WString::tr(key).arg( "<b>" + Wt::Utils::htmlEncode(nucname) + "</b>" )
                                .arg( Wt::Utils::htmlEncode(location) )
                                .arg( Wt::Utils::htmlEncode(msg.toUTF8()) ).toUTF8();
    item.plain = WString::tr(key).arg( nucname ).arg( location ).arg( msg ).toUTF8();
    item.field = field;
    item.isError = is_error;
    return item;
  };//source_item

  //Go through and make MakeDrfChart::DataPoint points, and also make sure to
  //  indicate which peaks are being background subtracted.
  size_t num_peaks_selected = 0;
  float minenergy = 999999.0f, maxenergy = -9999999.0f;
  vector<MakeDrfChart::DataPoint> datapoints;
  vector<pair<string,float>> datapoint_fit_keys;  //source key and energy of a used point; parallel to `datapoints`
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
      

      const vector<pair<DrfPeak *,MakeDrfSrcDef *>> peaks_to_sources = sample->selected_peak_to_sources();
      num_peaks_selected += peaks_to_sources.size();

      // The sample's peaks and sources, as MakeDrfCalc::samplePoints takes them.  A source whose
      //  inputs have a problem is left out, but only it - the other sources in the sample are still
      //  used (save for peaks it would add to).
      const vector<DrfPeak *> sample_peaks = sample->peaks();
      vector<MakeDrfSrcDef *> src_widgets;  //parallel to `input.sources`
      vector<Wt::Core::observing_ptr<WFormWidget>> src_problem_field;  //an unusable source's first problem input
      MakeDrfCalc::SampleInput input;
      input.liveTime = sample->liveTime();
      input.realTime = sample->realTime();

      auto read_source = [&]( MakeDrfSrcDef * const src ) -> MakeDrfCalc::CalSource {
        MakeDrfCalc::CalSource cs;
        cs.nuclide = src->nuclide();

        Wt::Core::observing_ptr<WFormWidget> first_problem;
        for( const MakeDrfSrcDef::InputProblem &p : src->validate() )
        {
          problems.push_back( source_item( src, sample, p.message, p.field, p.blocking ) );
          if( p.blocking && cs.usable )
          {
            cs.usable = false;
            first_problem = p.field;
          }
        }//for( loop over problems with this source )
        src_problem_field.push_back( first_problem );

        if( !cs.usable )
          return cs;

        try
        {
          cs.activity = src->activityAtSpectrumTime();
          cs.fracActivityUncert = src->fractionalActivityUncertainty();
          if( cs.nuclide )
            cs.age = src->ageAtSpectrumTime();
          if( !is_fixed_geometry )
          {
            cs.distance = src->distance();
            cs.distanceUncert = src->distanceUncertainty();
          }

          ShieldingSelect * const shield = src->shielding();
          if( shield && shield->isGenericMaterial() )
          {
            cs.shieldAtomicNumber = shield->atomicNumber();
            cs.shieldArealDensity = shield->arealDensity();
          }else if( shield )
          {
            cs.shieldMaterial = shield->material();
            if( cs.shieldMaterial )
              cs.shieldThickness = shield->thickness();
          }
        }catch( std::exception &e )
        {
          problems.push_back( source_item( src, sample, WString::fromUTF8(e.what()), nullptr, true ) );
          cs.usable = false;
        }//try / catch

        return cs;
      };//read_source lambda

      for( DrfPeak * const drfPeak : sample_peaks )
      {
        MakeDrfCalc::SamplePeak sp;
        sp.peak = drfPeak->m_peak;

        const auto pos = std::find_if( begin(peaks_to_sources), end(peaks_to_sources),
                                [drfPeak]( const pair<DrfPeak *,MakeDrfSrcDef *> &pp ){ return pp.first == drfPeak; } );
        if( pos != end(peaks_to_sources) )
        {
          MakeDrfSrcDef * const src = pos->second;
          const auto src_pos = std::find( begin(src_widgets), end(src_widgets), src );
          sp.source = static_cast<int>( src_pos - begin(src_widgets) );
          if( src_pos == end(src_widgets) )
          {
            src_widgets.push_back( src );
            input.sources.push_back( read_source( src ) );
          }

          // A peak with no nuclide uses the branching ratio the user entered for it
          if( !src->nuclide() && drfPeak->m_userBr
             && (drfPeak->m_userBr->validate() == WValidator::State::Valid) )
            sp.userBr = drfPeak->m_userBr->value();
        }//if( the peak is used, for a source )

        input.peaks.push_back( sp );
      }//for( DrfPeak * const drfPeak : sample_peaks )

      const MakeDrfCalc::SampleResult result = MakeDrfCalc::samplePoints( input, backgroundpeaks, doAirAtten );

      for( const pair<size_t,string> &err : result.sourceErrors )
        problems.push_back( source_item( src_widgets[err.first], sample, WString::fromUTF8(err.second), nullptr, true ) );

      // Peaks left out for each (usable source, unusable source) pair, so one message lists them
      map<pair<size_t,size_t>,vector<double>> interfered_peaks;

      for( const MakeDrfCalc::PeakPoint &pt : result.points )
      {
        DrfPeak * const drfPeak = sample_peaks[pt.peak];
        const std::shared_ptr<const PeakDef> &peak = drfPeak->m_peak;
        const size_t src_index = static_cast<size_t>( input.peaks[pt.peak].source );
        MakeDrfSrcDef * const srcDef = src_widgets[src_index];
        const SandiaDecay::Nuclide * const nuc = srcDef->nuclide();

        char energy_str[32] = { '\0' };
        snprintf( energy_str, sizeof(energy_str), "%.1f", pt.energy );

        pair<string,float> fit_key;  //set for a peak that gives an efficiency point
        switch( pt.status )
        {
          case MakeDrfCalc::PeakPoint::Status::NoUserBranchRatio:
            snprintf( energy_str, sizeof(energy_str), "%.2f", peak->mean() );
            problems.push_back( source_item( srcDef, sample, WString::tr("md-err-peak-br").arg(energy_str),
                                             drfPeak->m_userBr, true ) );
            continue;

          case MakeDrfCalc::PeakPoint::Status::NeedsUnusableSource:
            interfered_peaks[{src_index, static_cast<size_t>(pt.interferer)}].push_back( pt.energy );
            continue;

          // The peak gives no efficiency; say why it is not being used rather than drop it silently
          //  (its width is still good for the FWHM fit)
          case MakeDrfCalc::PeakPoint::Status::NoSourceGammas:
            problems.push_back( source_item( srcDef, sample,
                                             WString::tr("md-warn-peak-no-src-gammas").arg(energy_str),
                                             nullptr, false ) );
            break;

          case MakeDrfCalc::PeakPoint::Status::NoNetArea:
            problems.push_back( source_item( srcDef, sample,
                                             WString::tr("md-warn-peak-no-net-area").arg(energy_str),
                                             nullptr, false ) );
            break;

          case MakeDrfCalc::PeakPoint::Status::Used:
          {
            // Per-source keys keep each source's certificate/distance error a common-mode block
            auto used_pos = std::find_if( begin(used_sources), end(used_sources),
                     [srcDef]( const pair<MakeDrfSrcDef *,MakeDrfCalc::CalSource> &u ){ return u.first == srcDef; } );
            if( used_pos == end(used_sources) )
            {
              used_sources.emplace_back( srcDef, input.sources[src_index] );
              used_pos = end(used_sources) - 1;
            }

            MeasuredEffPoint raw = pt.measured;
            raw.sourceKey = (nuc ? nuc->symbol : string("user"))
                            + "#" + std::to_string( used_pos - begin(used_sources) );
            raw.fileName = SpecUtils::filename( fileWidget->measurement()->filename() );
            raw.sampleNumbers = SpecUtils::sequencesToBriefString( sample->samples() );
            fit_key = { raw.sourceKey, raw.energy };
            measured_points.push_back( std::move(raw) );
            break;
          }
        }//switch( pt.status )

        const bool sub_back = (pt.bkgLiveTime > 0.0);
        drfPeak->setBackgroundBeingSubtractedInfo( sub_back, (sub_back ? (pt.bkgArea / pt.bkgLiveTime) : 0.0) );

        MakeDrfChart::DataPoint point;
        point.energy = static_cast<float>( pt.energy );
        point.livetime = static_cast<float>( input.liveTime );
        point.peak_area = static_cast<float>( pt.peakArea );
        point.peak_area_uncertainty = static_cast<float>( pt.peakAreaUncert );
        if( peak->type() == PeakDef::GaussianDefined )
        {
          point.peak_fwhm = peak->fwhm();
          point.peak_fwhm_uncertainty = 2.35482*peak->uncertainty(PeakDef::CoefficientType::Sigma);
        }else
        {
          point.peak_fwhm = point.peak_fwhm_uncertainty = 0.0f;
        }
        point.source_count_rate = static_cast<float>( pt.sourceRate );
        point.source_count_rate_uncertainty = static_cast<float>( pt.sourceRate * pt.fracSourceUncert );
        point.distance = static_cast<float>( input.sources[src_index].distance );
        point.peak_color = peak->lineColor();
        point.background_peak_area = static_cast<float>( pt.bkgArea );
        point.background_peak_live_time = static_cast<float>( pt.bkgLiveTime );
        point.use_for_fwhm_fit = pt.useForFwhm;

        char counts[32] = { '\0' };
        snprintf( counts, sizeof(counts), "%.1f", point.peak_area );
        if( nuc )
        {
          point.source_information = WString::tr("md-peak-counts")
                                      .arg(nuc->symbol)
                                      .arg(energy_str)
                                      .arg(counts)
                                      .toUTF8();
        }else
        {
          point.source_information = WString::tr("md-peak-no-nuc")
                                      .arg(energy_str)
                                      .arg(counts)
                                      .toUTF8();
        }

        if( pt.useForFwhm )
          fwhm_peaks.push_back( peak );
        peaks.push_back( peak );
        datapoints.push_back( point );
        datapoint_fit_keys.push_back( fit_key );
      }//for( const MakeDrfCalc::PeakPoint &pt : result.points )

      for( const auto &interfered : interfered_peaks )
      {
        string energies;
        for( const double energy : interfered.second )
        {
          char buffer[32] = { '\0' };
          snprintf( buffer, sizeof(buffer), "%.1f", energy );
          energies += (energies.empty() ? "" : ", ") + string(buffer);
        }

        const size_t bad_index = interfered.first.second;
        const SandiaDecay::Nuclide * const bad_nuc = src_widgets[bad_index]->nuclide();
        problems.push_back( source_item( src_widgets[interfered.first.first], sample,
                                         WString::tr("md-warn-peak-needs-bad-src")
                                           .arg( energies ).arg( bad_nuc ? bad_nuc->symbol : string() ),
                                         src_problem_field[bad_index].get(), false ) );
      }//for( const auto &interfered : interfered_peaks )
    }//for( DrfSpecFileSample *sample : sampleWidgets )
  }//for( auto w : m_files->children()  )
  
  
  if( maxenergy < minenergy )
  {
    minenergy = 0.0;
    maxenergy = 3000.0;
  }
  
  // The sources table (one entry per calibration source), for the DRFs provenance
  vector<MeasuredSourceInfo> sources;
  for( size_t i = 0; i < used_sources.size(); ++i )
  {
    const MakeDrfCalc::CalSource &cs = used_sources[i].second;
    MeasuredSourceInfo info;
    info.sourceKey = (cs.nuclide ? cs.nuclide->symbol : string("user")) + "#" + std::to_string( i );
    info.nuclide = cs.nuclide ? cs.nuclide->symbol : string();
    info.activity = cs.activity;
    info.fracActivityUncert = static_cast<float>( cs.fracActivityUncert );
    if( cs.nuclide )
      info.age = cs.age;
    info.distance = static_cast<float>( cs.distance );
    info.distanceUncert = static_cast<float>( cs.distanceUncert );
    if( cs.shieldMaterial )
    {
      info.shieldMaterial = cs.shieldMaterial->name;
      info.shieldAtomicNumber = static_cast<float>( cs.shieldMaterial->massWeightedAtomicNumber() );
      info.shieldArealDensity = static_cast<float>( cs.shieldMaterial->density * cs.shieldThickness );
    }else
    {
      info.shieldAtomicNumber = static_cast<float>( cs.shieldAtomicNumber );
      info.shieldArealDensity = static_cast<float>( cs.shieldArealDensity );
    }

    try
    {
      info.assayInfo = used_sources[i].first->assayInfoString();
    }catch( std::exception & )
    {
    }
    sources.push_back( std::move(info) );
  }//for( size_t i = 0; i < used_sources.size(); ++i )

  MeasuredDrfPoints measured;
  try
  {
    measured.setPoints( measured_points );
    measured.setSources( sources );
  }catch( std::exception &e )
  {
    cerr << "handleSourcesUpdates: failed to set measured points: " << e.what() << endl;
  }
  m_measuredPoints = measured;

  // The points the equation is fit to: far-field intrinsic, in the equations energy units
  const bool inMeV = isEffEqnInMeV();
  vector<MakeDrfFit::EffFitPoint> effpoints;
  if( geomErrMsg.empty() && !measured.empty() )
  {
    try
    {
      effpoints = MakeDrfCalc::intrinsicFitPoints( measured, geom_choice, inMeV );
    }catch( std::exception &e )
    {
      geomErrMsg = e.what();
      problems.push_back( makeStatusItem( WString::fromUTF8(geomErrMsg), nullptr, true ) );
    }
  }//if( geomErrMsg.empty() && !measured.empty() )

  m_inputProblems = std::move( problems );
  m_numPeaksSelected = num_peaks_selected;
  m_numPeaksUsed = effpoints.size();
  m_peaksResolutionType = peaks.empty() ? PeakFitUtils::CoarseResolutionType::Unknown
                                        : PeakFitUtils::coarse_resolution_from_peaks( peaks );
  
  
  // The equation term counts on offer follow the number of points.  By default two degrees of
  //  freedom are left (one, for three points), so chi2/dof can say how well the curve describes the
  //  points; a count the user picked is kept, as far as there are points for it.
  const auto update_terms_combo = []( WComboBox *combo, const int npoints, const int max_offered,
                                      const int max_default, const int user_terms ){
    const int num_offered = std::min( npoints, max_offered );
    if( combo->count() != num_offered )
    {
      combo->clear();
      for( int order = 0; order < num_offered; ++order )
        combo->addItem( WString::tr("md-n-fit-pars").arg(order+1) );
    }

    if( num_offered < 1 )
      return;

    const int default_terms = MakeDrfCalc::defaultEquationTerms( npoints, max_default );
    const int terms = (user_terms > 0) ? std::min( user_terms, num_offered ) : default_terms;
    combo->setCurrentIndex( terms - 1 );
  };//update_terms_combo

  update_terms_combo( m_effEqnOrder, static_cast<int>(effpoints.size()), 8, 7, m_effEqnUserTerms );
  update_terms_combo( m_sqrtEqnOrder, static_cast<int>(fwhm_peaks.size()), 6, 5, m_sqrtEqnUserTerms );
  
  m_dataPoints = datapoints;
  m_dataDiameter = diameter;

  // Chart points: the intrinsic efficiency each peak gave, with the total 1-sigma as its error bar,
  //  and the peaks FWHM.  A fit point is found by its source key and energy (converted exactly as
  //  intrinsicFitPoints does); points sharing both keep their order through the fit.
  {
    assert( datapoint_fit_keys.size() == datapoints.size() );
    map<pair<string,float>,deque<const MakeDrfFit::EffFitPoint *>> fit_points;
    for( const MakeDrfFit::EffFitPoint &fp : effpoints )
      fit_points[{fp.sourceKey, fp.energy}].push_back( &fp );

    vector<DrfChartPoint> chart_pts;
    for( size_t index = 0; index < datapoints.size(); ++index )
    {
      const MakeDrfChart::DataPoint &d = datapoints[index];
      DrfChartPoint cp;
      cp.energy = d.energy;
      cp.fwhm = d.use_for_fwhm_fit ? d.peak_fwhm : 0.0f;  //only the widths the FWHM fit uses
      cp.fwhmUncert = d.use_for_fwhm_fit ? d.peak_fwhm_uncertainty : 0.0f;
      cp.label = d.source_information;
      cp.color = d.peak_color.isDefault() ? string() : d.peak_color.cssText();

      const pair<string,float> &key = datapoint_fit_keys[index];
      const auto pos = key.first.empty() ? end(fit_points)
                         : fit_points.find( {key.first, inMeV ? (key.second / 1000.0f) : key.second} );
      if( (pos != end(fit_points)) && !pos->second.empty() )
      {
        const MakeDrfFit::EffFitPoint &fp = *pos->second.front();
        pos->second.pop_front();
        cp.efficiency = fp.efficiency;
        cp.efficiencyUncert = fp.efficiency * sqrt( fp.fracStatUncert*fp.fracStatUncert
                                                    + fp.fracCertUncert*fp.fracCertUncert
                                                    + fp.fracDistUncert*fp.fracDistUncert );
      }
      chart_pts.push_back( std::move(cp) );
    }//for( size_t index = 0; index < datapoints.size(); ++index )

    m_chart->setDataPoints( chart_pts );
  }
  
  // The fit options stay where they are; only the choices that need data wait for it
  m_effEqnOrder->setDisabled( datapoints.empty() );
  m_sqrtEqnOrder->setDisabled( datapoints.empty() );
  m_fitNeedsDataTxt->setHidden( !datapoints.empty() );
  
  fitEffEqn( effpoints );
  fitFwhmEqn( fwhm_peaks );

  updateStatusDisplay();
  
  wApp->triggerUpdate();
}//void handleSourcesUpdates()


void MakeDrf::handleFwhmTypeChanged()
{
  const auto fwhm_form = DetectorPeakResponse::ResolutionFnctForm( m_fwhmEqnType->currentIndex() );
  
  bool valid_fwhm_form = false;
  switch( fwhm_form )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
    case DetectorPeakResponse::kConstantPlusSqrtEnergy:
    case DetectorPeakResponse::kNumResolutionFnctForm:
      m_sqrtEqnOrder->hide();
      m_sqrtEqnOrderLabel->hide();
      valid_fwhm_form = true;
      break;
      
    case DetectorPeakResponse::kSqrtPolynomial:
      m_sqrtEqnOrder->show();
      m_sqrtEqnOrderLabel->show();
      valid_fwhm_form = true;
      break;
  }//switch( fwhm_form )
  
  assert( valid_fwhm_form );
  
  //We could update *just* the FWHM equation like:
  //size_t numchan = 0;
  //vector< std::shared_ptr<const PeakDef> > peaks;
  // ... get peaks and numchan ...
  //fitFwhmEqn( peaks, numchan );
  //But instead we'll be lazy and update everything
  
  handleSourcesUpdates();
}//void handleFwhmTypeChanged()


void MakeDrf::handleEffEqnOrderChange()
{
  m_effEqnUserTerms = m_effEqnOrder->currentIndex() + 1;
  handleSourcesUpdates();
}//void handleEffEqnOrderChange()


void MakeDrf::handleSqrtEqnOrderChange()
{
  m_sqrtEqnUserTerms = m_sqrtEqnOrder->currentIndex() + 1;

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
  m_chart->setShowFwhmPoints( m_showFwhmPoints->isChecked() );
}//void handleShowFwhmPointsToggled()


void MakeDrf::handleFixedGeometryChanged()
{
  const bool is_fixed_geom = (m_geometry->currentIndex() != 0);

  // Fixed-geometry efficiencies dont use the detectors size or shape, nor the air between it and
  //  the source; say so, rather than have those inputs vanish.
  m_detSectionBody->setHidden( is_fixed_geom );
  m_detNotUsedTxt->setHidden( !is_fixed_geom );
  m_airAttenuate->setDisabled( is_fixed_geom );
  updateSidebarWidth();
  
  DetectorPeakResponse::EffGeometryType geom_type = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
  switch( m_geometry->currentIndex() )
  {
    case 0: geom_type = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;   break;
    case 1: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct;   break;
    case 2: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2;  break;
    case 3: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2;   break;
    case 4: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram; break;
    default: assert( 0 ); break;
  }//switch( m_geometry->currentIndex() )
  
  // Go through all the sources and set distances visible/hidden.
  //  Each source emits "updated" as we change it, which would re-enter handleSourcesUpdates()
  //  while the sources we havent gotten to yet still show the previous geometrys inputs - so
  //  block those updates, and do the single re-fit ourselves below.
  {
    const BlockSourceUpdates block( m_blockSourceUpdates );

    for( auto w : m_files->children() )
    {
      DrfSpecFile *fileWidget = dynamic_cast<DrfSpecFile *>( w );
      if( fileWidget )
        fileWidget->setIsEffGeometryType( geom_type );
    }//for( auto w : m_files->children() )
  }

  handleSourcesUpdates();
}//void handleFixedGeometryChanged()


void MakeDrf::handleGeometryModeChanged()
{
  const bool full = m_geomMode->isChecked();

  if( full && !m_geomSeeded )
  {
    // Carry the diameter/setback over into the geometry form, once
    double diam = 0.0;
    try{ diam = PhysicalUnits::stringToDistance( m_detDiameter->text().toUTF8() ); }catch( std::exception & ){}
    m_mcTool->geometryInput()->seedFromDiameter( diam, enteredSetback(), detectorTypeHint() );
    m_geomSeeded = true;
  }else if( !full )
  {
    // ...and back-fill the diameter from the geometry
    const double d = m_mcTool->geometryInput()->enteredDiameter();
    if( d > 0.0 )
      m_detDiameter->setText( WString::fromUTF8( PhysicalUnits::printToBestLengthUnits(d) ) );
  }

  m_diameterDiv->setHidden( full );
  m_mcTool->setHidden( !full );
  m_geomGuidance->setHidden( full );

  // The geometry form needs a wider sidebar; the chart keeps whatever height the user gave it.
  updateSidebarWidth();

  handleSourcesUpdates();
}//handleGeometryModeChanged()


void MakeDrf::updateSidebarWidth()
{
  if( m_phoneTabs || !m_sidebar )
    return;  //on phones the sidebar is a tab, and fills it

  const bool fixed_geom = (m_geometry && (m_geometry->currentIndex() != 0));
  const bool full_geom = !fixed_geom && m_geomMode && m_geomMode->isChecked();

  int width = 340;
  if( full_geom )
  {
    // Enough for the geometry form's layer table; but leave the chart and peaks most of the window
    const int app_width = m_interspec ? m_interspec->renderedWidth() : 0;
    width = (app_width > 100) ? std::clamp( (45*make_drf_dialog_width(app_width))/100, 340, 450 ) : 450;
  }

  m_sidebar->setWidth( WLength(width, WLength::Unit::Pixel) );
}//void updateSidebarWidth()


std::string MakeDrf::detectorTypeHint() const
{
  // The crystal the spectra's instrument implies (its model text may name one, e.g. "LaBr3"),
  //  else the one the fitted peak widths imply.
  string hint;
  for( auto w : m_files->children() )
  {
    auto fileWidget = dynamic_cast<DrfSpecFile *>( w );
    if( !fileWidget || !fileWidget->measurement() )
      continue;
    const shared_ptr<const SpecMeas> meas = fileWidget->measurement();
    const PeakFitUtils::CoarseResolutionType type
                               = PeakFitUtils::coarse_type_for_detector_type( meas->detector_type() );
    hint += " " + meas->instrument_model() + " " + CeeLoUtils::crystalNameForResolutionType( type );
  }

  if( CeeLoUtils::crystalNameFromText( hint ).empty() )
    hint += " " + CeeLoUtils::crystalNameForResolutionType( m_peaksResolutionType );

  return hint;
}//detectorTypeHint()


MakeDrfCalc::GeometryChoice MakeDrf::geometryChoice() const
{
  MakeDrfCalc::GeometryChoice choice;
  choice.fixedGeometry = (m_geometry->currentIndex() != 0);
  if( choice.fixedGeometry )
    return choice;

  if( m_geomMode && m_geomMode->isChecked() && m_mcTool )
  {
    if( !m_mcTool->geometryValid() )
      throw runtime_error( WString::tr("md-err-geometry-invalid").toUTF8() );

    const ceelo::GeometryDescriptor gd = m_mcTool->geometryDescriptor();
    choice.diameter = 2.0 * gd.transverse_half_extent() * PhysicalUnits::cm;
    choice.setback = gd.endcap_front_offset_cm() * PhysicalUnits::cm;

    // The full geometry counts only once it is more than the fabricated diameter guess (see
    //  DetectorGeometryInput::generationReady); until then it is a flat disk like diameter mode.
    if( m_mcTool->generationReady() )
      choice.geometry = make_shared<const ceelo::GeometryDescriptor>( gd );

    return choice;
  }//if( geometry mode )

  try
  {
    choice.diameter = PhysicalUnits::stringToDistance( m_detDiameter->text().toUTF8() );
  }catch( std::exception & )
  {
    throw runtime_error( WString::tr("md-err-det-diam-invalid").toUTF8() );
  }
  if( (choice.diameter <= 0.0) || IsNan(choice.diameter) || IsInf(choice.diameter) )
    throw runtime_error( WString::tr("md-err-det-diam-invalid").toUTF8() );

  // Blank means no setback; anything else has to be a real one - a typo is not silently zero
  const string setback_txt = SpecUtils::trim_copy( m_detSetback->text().toUTF8() );
  if( !setback_txt.empty() )
  {
    double setback = -1.0;
    try
    {
      setback = PhysicalUnits::stringToDistance( setback_txt );
    }catch( std::exception & )
    {
    }

    if( !(setback >= 0.0) || IsInf(setback) )
      throw runtime_error( WString::tr("md-err-setback-invalid").toUTF8() );
    choice.setback = setback;
  }//if( !setback_txt.empty() )

  return choice;
}//geometryChoice()


MakeDrf::StatusItem MakeDrf::makeStatusItem( const Wt::WString &message, Wt::WFormWidget *field,
                                             const bool is_error )
{
  StatusItem item;
  item.plain = message.toUTF8();
  item.html = Wt::Utils::htmlEncode( item.plain );
  item.field = field;
  item.isError = is_error;
  return item;
}//makeStatusItem(...)


std::vector<MakeDrf::StatusItem> MakeDrf::validateDetectorInputs()
{
  vector<StatusItem> items;
  bool diameter_value_error = false, setback_value_error = false;

  const bool fixed_geom = (m_geometry->currentIndex() != 0);
  const bool full_geom = m_geomMode->isChecked();

  if( fixed_geom )
  {
    // The detector description isnt used
  }else if( full_geom )
  {
    if( !m_mcTool->geometryValid() )
    {
      const string problem = m_mcTool->geometryProblem();
      const WString msg = problem.empty() ? WString::tr("md-err-geometry-invalid")
                                          : WString::tr("md-err-geometry").arg( problem );
      items.push_back( makeStatusItem( msg, m_mcTool->geometryInput()->problemField(), true ) );
    }
  }else
  {
    const string diam_txt = SpecUtils::trim_copy( m_detDiameter->text().toUTF8() );
    if( diam_txt.empty() )
    {
      items.push_back( makeStatusItem( WString::tr("md-err-det-diam-missing"), m_detDiameter, true ) );
      diameter_value_error = true;
    }else
    {
      try
      {
        const double diam = PhysicalUnits::stringToDistance( diam_txt );
        if( !(diam > 0.0) || IsInf(diam) )
        {
          items.push_back( makeStatusItem( WString::tr("md-err-det-diam-nonpos"), m_detDiameter, true ) );
          diameter_value_error = true;
        }else if( (diam < 0.5*PhysicalUnits::cm) || (diam > 150.0*PhysicalUnits::cm) )
        {
          // Usable, but likely a unit slip (e.g. "3" meant as inches, but taken as cm, is fine;
          //  "76" meant as mm but taken as cm is not)
          items.push_back( makeStatusItem( WString::tr("md-warn-det-diam-unusual")
                                             .arg( PhysicalUnits::printToBestLengthUnits(diam) ),
                                           m_detDiameter, false ) );
        }
      }catch( std::exception & )
      {
        items.push_back( makeStatusItem( WString::tr("md-err-det-diam-invalid"), m_detDiameter, true ) );
      }
    }//if( diam_txt.empty() ) / else

    const string setback_txt = SpecUtils::trim_copy( m_detSetback->text().toUTF8() );
    if( !setback_txt.empty() )
    {
      try
      {
        const double setback = PhysicalUnits::stringToDistance( setback_txt );
        if( !(setback >= 0.0) || IsInf(setback) )
        {
          items.push_back( makeStatusItem( WString::tr("md-err-setback-invalid"), m_detSetback, true ) );
          setback_value_error = true;
        }
      }catch( std::exception & )
      {
        items.push_back( makeStatusItem( WString::tr("md-err-setback-invalid"), m_detSetback, true ) );
      }
    }//if( !setback_txt.empty() )
  }//if( fixed_geom ) / else if( full_geom ) / else

  WidgetUtils::setInputValueError( m_detDiameter, diameter_value_error );
  WidgetUtils::setInputValueError( m_detSetback, setback_value_error );

  return items;
}//validateDetectorInputs()


std::vector<MakeDrf::StatusItem> MakeDrf::detectorTypeMismatch() const
{
  vector<StatusItem> items;
  if( !m_peakFitDetPrefsGui || (m_geometry->currentIndex() != 0) )
    return items;

  const shared_ptr<const PeakFitDetPrefs> prefs = m_peakFitDetPrefsGui->currentPrefs();
  if( !prefs || (prefs->m_det_type == PeakFitUtils::CoarseResolutionType::Unknown) )
    return items;

  // What the detector looks like: its crystal, when the geometry is given, else its peak widths
  const bool full_geom = m_geomMode->isChecked();
  const string crystal = full_geom ? m_mcTool->geometryInput()->crystalMaterialName() : string();
  const PeakFitUtils::CoarseResolutionType seen = full_geom ? crystal_resolution_type( crystal )
                                                            : m_peaksResolutionType;
  if( !resolution_types_conflict( prefs->m_det_type, seen ) )
    return items;

  const WString prefs_type = PeakFitDetPrefsGui::detectorTypeLabel( prefs->m_det_type );
  const WString msg = full_geom
         ? WString::tr("md-warn-det-type-crystal").arg( prefs_type ).arg( crystal )
         : WString::tr("md-warn-det-type-peaks").arg( prefs_type )
                                                .arg( PeakFitDetPrefsGui::detectorTypeLabel(seen) );

  items.push_back( makeStatusItem( msg, m_peakFitDetPrefsGui->detectorTypeCombo(), false ) );
  return items;
}//detectorTypeMismatch()


void MakeDrf::updateStatusDisplay()
{
  m_problemList->clear();

  vector<StatusItem> items = m_inputProblems;
  for( const StatusItem &i : detectorTypeMismatch() )
    items.push_back( i );
  for( const StatusItem &i : m_fitProblems )
    items.push_back( i );

  for( const StatusItem &item : items )
  {
    WText *row = m_problemList->addNew<WText>( WString::fromUTF8(item.html), Wt::TextFormat::XHTML );
    row->setInline( false );
    row->addStyleClass( "MakeDrfProblem" );
    row->addStyleClass( item.isError ? "ErrorTxt" : "WarnTxt" );

    if( item.field )
    {
      row->addStyleClass( "HasField" );
      const Wt::Core::observing_ptr<WFormWidget> field = item.field;
      row->clicked().connect( this, [this,field](){ focusProblemField( field.get() ); } );
    }
  }//for( const StatusItem &item : items )

  m_problemList->setHidden( items.empty() );

  // A one-liner on how much of what was selected is going into the fit
  WString summary;
  if( m_numPeaksSelected == 0 )
    summary = WString::tr("md-status-no-peaks");
  else if( m_numPeaksUsed < m_numPeaksSelected )
    summary = WString::tr("md-status-using").arg( static_cast<int>(m_numPeaksUsed) )
                                             .arg( static_cast<int>(m_numPeaksSelected) );

  m_statusSummary->setText( summary );
  m_statusSummary->setHidden( summary.empty() );
}//void updateStatusDisplay()


void MakeDrf::focusProblemField( Wt::WFormWidget *field )
{
  if( !field )
    return;

  // On phones, first go to the tab the input is on
  if( m_phoneTabs )
  {
    for( int i = 0; i < m_phoneTabs->count(); ++i )
    {
      const WWidget * const content = m_phoneTabs->widget( i );
      for( const WWidget *p = field; p; p = p->parent() )
      {
        if( p == content )
        {
          m_phoneTabs->setCurrentIndex( i );
          break;
        }
      }
    }//for( loop over tabs )
  }//if( m_phoneTabs )

  WidgetUtils::revealAndFocus( field );
}//void focusProblemField( Wt::WFormWidget *field )


void MakeDrf::handlePeakFitPrefsChanged()
{
  if( !m_peakFitDetPrefsGui || !m_pfpPanel )
    return;

  m_peakFitDetPrefsGui->handlePrefsChanged();

  const string summary = m_peakFitDetPrefsGui->summaryText().toUTF8();
  m_pfpPanel->setTitle( WString::fromUTF8( WString::tr("md-peak-fit-prefs-title").toUTF8()
                                           + " <span class=\"MakeDrfPfpSummary\">"
                                           + Wt::Utils::htmlEncode( summary ) + "</span>" ) );

  updateStatusDisplay();
}//void handlePeakFitPrefsChanged()


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


void MakeDrf::fitFwhmEqn( std::vector< std::shared_ptr<const PeakDef> > peaks )
{
  ++m_fwhmFitId;
  
  const int fitid = m_fwhmFitId;
  const string sessionId = wApp->sessionId();
  
  m_fwhmEqnChi2 = -999.9;
  m_fwhmCoefs.clear();
  m_fwhmCoefUncerts.clear();
  updateChartDetector();

  int sqrtEqnOrder = -1;
  const auto fnctnlForm = DetectorPeakResponse::ResolutionFnctForm(m_fwhmEqnType->currentIndex());
  
  switch( fnctnlForm )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
    case DetectorPeakResponse::kNumResolutionFnctForm:
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
    case DetectorPeakResponse::kConstantPlusSqrtEnergy:
      break;
      
    case DetectorPeakResponse::kSqrtPolynomial:
      sqrtEqnOrder = m_sqrtEqnOrder->currentIndex() + 1;
      break;
      
    default:
      assert( 0 );
      return;
  }//switch( fnctnlForm )
  
  
  //ToDo: I'm not entirely sure the next line protects against updateFwhmEqn()
  //  not being called if this widget is deleted before fit is done.
  auto updater = [this, fnctnlForm = static_cast<int>(fnctnlForm), fitid](
      std::vector<float> coefs, std::vector<float> uncerts, double chi2 ){
    updateFwhmEqn( coefs, uncerts, chi2, fnctnlForm, fitid );
  };
  
  const string thisid = id();
  MakeDrf * const self = this;
  
  auto worker = [sessionId,fnctnlForm,peaks,sqrtEqnOrder,updater,thisid,self]() {
    try
    {
      auto peakdequ = std::make_shared<std::deque< std::shared_ptr<const PeakDef> > >( peaks.begin(), peaks.end() );
    
      //Takes between 5 and 500ms for a HPGe detector
      const double start_time = SpecUtils::get_wall_time();
      vector<float> fwhm_coefs, fwhm_coefs_uncert;
      const double chi2 = MakeDrfFit::performResolutionFit( peakdequ, fnctnlForm, sqrtEqnOrder, fwhm_coefs, fwhm_coefs_uncert );
    
      const double end_time = SpecUtils::get_wall_time();
    
      assert( fwhm_coefs.size() == fwhm_coefs_uncert.size() );
      (void)start_time;
      (void)end_time;
      
      WServer::instance()->post( sessionId, std::bind( [updater,fwhm_coefs,fwhm_coefs_uncert,thisid,chi2,self](){
        if( wApp->domRoot() && (dynamic_cast<MakeDrf *>(wApp->domRoot()->findById(thisid)) == self) )
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
  
  const auto fcntform = DetectorPeakResponse::ResolutionFnctForm( functionalForm );
  if( fcntform == DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm )
  {
    m_fwhmEqnType->setCurrentIndex( -1 );
    return;
  }

  m_fwhmEqnType->setCurrentIndex( static_cast<int>(fcntform) );
  m_fwhmEqnChi2 = chi2;
  m_fwhmCoefs = coefs;
  m_fwhmCoefUncerts = uncerts;
  updateChartDetector();
  
  wApp->triggerUpdate();
}//void updateFwhmEqn(...)


void MakeDrf::fitEffEqn( std::vector<MakeDrfFit::EffFitPoint> data )
{
  ++m_effEqnFitId;
  
  const int fitid = m_effEqnFitId;
  const string sessionId = wApp->sessionId();
  
  m_effEqnChi2 = -999.9;
  m_effEqnDof = 0;
  m_effEqnCoefs.clear();
  m_effEqnCoefUncerts.clear();
  m_effEqnCov.clear();
  m_effBirgeScale = 1.0;
  m_intrinsicEfficiencyIsValid.emit( false );
  m_intrinsicEffAnswer->setText( "" );
  m_intrinsicEffAnswer->hide();
  m_fitProblems.clear();

  // With nothing to fit, the data range is gone too (else the chart keeps shading the old one)
  if( data.empty() )
  {
    m_effLowerEnergy = m_effUpperEnergy = 0.0f;
    m_effEqnChi2 = -999.9;
  }

  updateChartDetector();
  
  const int eqnOrderIndex = m_effEqnOrder->currentIndex();
  if( data.empty() || (eqnOrderIndex < 0) )
    return;
  
  const int nfitpars = std::min( eqnOrderIndex+1, 8 );
  const bool inMeV = isEffEqnInMeV();
  
  // Re-resolve by DOM id rather than capturing `this`: the fit runs on a worker thread and the
  //  user may have closed the window before it finishes.
  auto updater = [this, fitid]( MakeDrfFit::EffFitResult result, float lowestEnergy,
                                float highestEnergy, std::string errmsg ){
    updateEffEqn( result, lowestEnergy, highestEnergy, fitid, errmsg );
  };
  const string thisid = id();
  MakeDrf * const self = this;
  
  auto worker = [sessionId,thisid,data,nfitpars,updater,inMeV,self]() {
    MakeDrfFit::EffFitResult result;
    string errmsg;
    try
    {
      result = MakeDrfFit::performEfficiencyFit( data, nfitpars );
    }catch( std::exception &e )
    {
      errmsg = e.what();
    }
    
    float lowestEnergy = 99999.0f, highestEnergy = -99999.0f;
    for( const MakeDrfFit::EffFitPoint &p : data )
    {
      lowestEnergy = std::min( lowestEnergy, (inMeV ? 1000.0f : 1.0f) * p.energy );
      highestEnergy = std::max( highestEnergy, (inMeV ? 1000.0f : 1.0f) * p.energy );
    }
    
    WServer::instance()->post( sessionId, std::bind( [updater,thisid,result,lowestEnergy,highestEnergy,errmsg,self](){
      //Make sure *this is still in the widget tree (incase user closed window while computation was being done)
      if( wApp->domRoot() && (dynamic_cast<MakeDrf *>(wApp->domRoot()->findById(thisid)) == self) )
        updater( result, lowestEnergy, highestEnergy, errmsg );
      else
        cerr << "MakeDrf widget was deleted while efficiency was being calculated" << endl;
    } ) );
  };
  
  WServer::instance()->ioService().boost::asio::io_service::post( worker );
}//void fitEffEqn(...)


void MakeDrf::updateEffEqn( MakeDrfFit::EffFitResult result,
                            const float lowestEnergy, const float highestEnergy,
                            const int fitid, const string errmsg )
{
  if( fitid != m_effEqnFitId )
    return;

  m_effEqnChi2 = result.chi2;
  m_effEqnDof = errmsg.empty() ? result.dof : 0;
  m_effLowerEnergy = lowestEnergy;
  m_effUpperEnergy = highestEnergy;
  m_effEqnCoefs = errmsg.empty() ? result.coefs : vector<float>{};
  m_effEqnCoefUncerts = errmsg.empty() ? result.uncerts : vector<float>{};
  m_effEqnCov = errmsg.empty() ? result.covRowMajor : vector<float>{};
  m_effBirgeScale = result.birgeScale;
  m_intrinsicEfficiencyIsValid.emit( !m_effEqnCoefs.empty() );
  
  // Fit problems go below the source/geometry ones handleSourcesUpdates() found
  m_fitProblems.clear();
  if( !errmsg.empty() )
    m_fitProblems.push_back( makeStatusItem( WString::fromUTF8(errmsg), nullptr, true ) );
  if( m_effBirgeScale > 1.5 )
  {
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%.1f", m_effBirgeScale );
    m_fitProblems.push_back( makeStatusItem( WString::tr("md-fit-warning-birge").arg(buffer), nullptr, false ) );
  }
  if( !result.warnings.empty() )
    m_fitProblems.push_back( makeStatusItem( WString::fromUTF8(result.warnings), nullptr, false ) );
  updateStatusDisplay();
  
  if( !m_effEqnCoefs.empty() )
  {
    const vector<float> &coefs = m_effEqnCoefs;
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
    eqn += ")&nbsp;&nbsp;" + WString::tr("md-eqn-x-units").arg( isEffEqnInMeV() ? "MeV" : "keV" ).toUTF8();
    
    m_intrinsicEffAnswer->setText( eqn );
    m_intrinsicEffAnswer->show();
  }else
  {
    m_intrinsicEffAnswer->setText( "" );
    m_intrinsicEffAnswer->hide();
  }//if( !coefs.empty() ) / else

  updateChartDetector();
  
  wApp->triggerUpdate();
}//void updateEffEqn(...)


void MakeDrf::updateChartDetector()
{
  shared_ptr<DetectorPeakResponse> preview;
  try
  {
    preview = assembleLegacyDrf();
  }catch( std::exception & )
  {
    //no valid fit yet
  }

  // The data range first: the chart scales its efficiency axis to the curve over that range
  m_chart->setDataRange( m_effLowerEnergy, m_effUpperEnergy );
  m_chart->updateChart( preview );

  if( !m_chartRangeSet && (m_effUpperEnergy > m_effLowerEnergy) )
  {
    const float span = m_effUpperEnergy - m_effLowerEnergy;
    m_chart->setXAxisRange( std::max( 0.0f, m_effLowerEnergy - 0.05f*span ), m_effUpperEnergy + 0.05f*span );
    m_chartRangeSet = true;
  }

  // The geometry/MC tool anchors and grounds to the DRF being made
  if( m_mcTool )
    m_mcTool->setSeedDrf( preview );
}//updateChartDetector()


MakeDrfCalc::FitResults MakeDrf::currentFitResults() const
{
  MakeDrfCalc::FitResults fit;
  fit.eff.coefs = m_effEqnCoefs;
  fit.eff.uncerts = m_effEqnCoefUncerts;
  fit.eff.covRowMajor = m_effEqnCov;
  fit.eff.chi2 = m_effEqnChi2;
  fit.eff.birgeScale = m_effBirgeScale;
  fit.effInMeV = isEffEqnInMeV();
  fit.lowerEnergy = m_effLowerEnergy;
  fit.upperEnergy = m_effUpperEnergy;
  fit.fwhm.coefs = m_fwhmCoefs;
  fit.fwhm.uncerts = m_fwhmCoefUncerts;
  fit.fwhm.chi2 = m_fwhmEqnChi2;
  fit.fwhmForm = m_fwhmCoefs.empty() ? DetectorPeakResponse::kNumResolutionFnctForm
                        : DetectorPeakResponse::ResolutionFnctForm( m_fwhmEqnType->currentIndex() );

  switch( m_geometry->currentIndex() )
  {
    case 0: fit.geometryType = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;   break;
    case 1: fit.geometryType = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct;   break;
    case 2: fit.geometryType = DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2;  break;
    case 3: fit.geometryType = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2;   break;
    case 4: fit.geometryType = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram; break;
    default: assert( 0 ); break;
  }//switch( m_geometry->currentIndex() )

  if( m_peakFitDetPrefsGui )
    fit.peakFitPrefs = m_peakFitDetPrefsGui->currentPrefs();

  return fit;
}//currentFitResults()


std::shared_ptr<DetectorPeakResponse> MakeDrf::assembleLegacyDrf() const
{
  MakeDrfCalc::GeometryChoice geom = geometryChoice();
  geom.geometry.reset();  //no response, no geometry: the plain fitted curve
  return MakeDrfCalc::assembleDrf( "preview", "", m_measuredPoints, geom, currentFitResults(), nullptr );
}//assembleLegacyDrf()


std::shared_ptr<SpecMeas> MakeDrf::assembleCalFile( const std::string &drfname,
                                                    const std::string &drfdescrip )
{
  auto answer = std::make_shared<SpecMeas>();
  
  const bool is_fixed_geometry = (m_geometry->currentIndex() != 0);
  
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

          // Everything about the source, in a form this tool can read back (the GADRAS form above
          //  cannot carry the uncertainties or assay details)
          try
          {
            const SandiaDecay::Nuclide * const nuc = source->nuclide();
            string remark = "InterSpec-Source: nuc=" + (nuc ? nuc->symbol : string("none"))
                            + "; act=" + PhysicalUnits::printToBestActivityUnits( source->activityAtSpectrumTime(), 6, true )
                            + "; actSig=" + std::to_string( 100.0*source->fractionalActivityUncertainty() ) + "%";
            if( nuc )
              remark += "; age=" + PhysicalUnits::printToBestTimeUnits( source->ageAtSpectrumTime(), 6 );
            if( !is_fixed_geometry )
            {
              remark += "; dist=" + PhysicalUnits::printToBestLengthUnits( source->distance(), 6 );
              if( source->distanceUncertainty() > 0.0 )
                remark += "; distSig=" + PhysicalUnits::printToBestLengthUnits( source->distanceUncertainty(), 6 );
            }
            ShieldingSelect * const shield = source->shielding();
            if( shield && shield->isGenericMaterial() )
              remark += "; shieldAN=" + std::to_string( shield->atomicNumber() )
                        + "; shieldAD=" + std::to_string( shield->arealDensity() / (PhysicalUnits::g/PhysicalUnits::cm2) ) + " g/cm2";
            else if( shield && shield->material() )
              remark += "; shield=" + shield->material()->name
                        + "; shieldThickness=" + PhysicalUnits::printToBestLengthUnits( shield->thickness(), 6 );
            const string assay = source->assayInfoString();
            if( !assay.empty() )
              remark += "; assay=" + assay;
            newremarks.push_back( remark );
          }catch( std::exception & )
          {
          }
          if( !is_fixed_geometry )
          {
            try
            {
              distance = source->distance();
            }catch(...){}
          }//if( !is_fixed_geometry )
        }//for( MakeDrfSrcDef *source : sample->sources() )
        
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
        }else if( is_fixed_geometry )
        {
          string title = SpecUtils::trim_copy(m->title());
          if( !SpecUtils::icontains(title, "Fixed") )
            title += " Fixed Geometry";
          answer->set_title( title, m );
        }//if( distance > 0.0 ) / else if( is_fixed_geometry )
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
    }

    // The DRF itself - equation, covariance, measured points with their sources, geometry and
    //  response - so the N42 file is the complete record of the characterization.
    try
    {
      answer->setDetector( assembleDrf( drfname.empty() ? string("Created DRF") : drfname, drfdescrip ) );
    }catch( std::exception &e )
    {
      cerr << "assembleCalFile: not embedding DRF: " << e.what() << endl;
    }
  }catch( std::exception &e )
  {
    cerr << "Caught exception in assembleCalFile(): " << e.what();
    return nullptr;
  }//try / catch
  
  if( answer->measurements().empty() )
    return nullptr;
  
  return answer;
}//std::shared_ptr<SpecMeas> assembleCalFile()


shared_ptr<DetectorPeakResponse> MakeDrf::assembleDrf( const string &name, const string &descrip ) const
{
  // A Monte-Carlo response generated in the embedded tool, or the one the tool currently holds
  //  (the measured-curve transfer is rebuilt from the fit by MakeDrfCalc::assembleDrf itself)
  shared_ptr<ceelo::DetectorResponse> mc_response = m_generatedMcResponse;
  if( mc_response && (mc_response->provenance.method == ceelo::ProductionMethod::CurveTransfer) )
    mc_response.reset();

  return MakeDrfCalc::assembleDrf( name, descrip, m_measuredPoints, geometryChoice(),
                                   currentFitResults(), mc_response );
}//std::make_shared<DetectorPeakResponse> assembleDrf() const;


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
  
  const vector<MakeDrfChart::DataPoint> data = m_dataPoints;
  
  const int effDof = m_effEqnDof;
  const int numFwhmPoints = static_cast<int>( std::count_if( begin(data), end(data),
                                      []( const MakeDrfChart::DataPoint &d ){ return d.use_for_fwhm_fit; } ) );
  const int fwhmDof = numFwhmPoints - static_cast<int>(fwhmCoefs.size());
  
  const bool effInMeV = isEffEqnInMeV();
  const auto resFcnForm = DetectorPeakResponse::ResolutionFnctForm( m_fwhmEqnType->currentIndex() );
  
  DetectorPeakResponse::EffGeometryType geom_type = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
  switch( m_geometry->currentIndex() )
  {
    case 0: geom_type = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;   break;
    case 1: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct;   break;
    case 2: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2;  break;
    case 3: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2;   break;
    case 4: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram; break;
    default: assert( 0 ); break;
  }//switch( m_geometry->currentIndex() )
  
  
  double diam = 0.0;
  const double setback = detectorSetback();
  try
  {
    diam = detectorDiameter();
  }catch( std::exception &e )
  {
    if( geom_type == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic )
    {
      out << "Invalid detector diameter: " << e.what() << endline;
      return;
    }
  }//try / catch

  if( effEqnCoefs.empty() )
  {
    out << "Detector response function not valid" << endline;
    
    for( const StatusItem &item : m_inputProblems )
      out << item.plain << endline;
    for( const StatusItem &item : m_fitProblems )
      out << item.plain << endline;
    return;
  }//if( effEqnCoefs.empty() )
  
  //Okay, we've got all the variables we need for this function I think, lets
  //  write stuff
  
  // The "Relative Eff" column: IEEE 325, of the DRF being saved (so a full geometry's response is
  //  used), with its uncertainty; empty where it has no meaning (fixed geometry, or a fit not
  //  reaching 1332 keV).
  string rel_eff_field, rel_eff_uncert_field;
  bool rel_eff_flat_disk = false;
  try
  {
    const shared_ptr<const DetectorPeakResponse> drf = assembleDrf( drfname, drfdescription );
    const double rel_eff = drf->relativeEfficiencyIeee325();
    if( rel_eff > 0.0 )
    {
      rel_eff_flat_disk = !drf->ceeloResponse();
      rel_eff_field = SpecUtils::printCompact( rel_eff, 4 ) + "%";
      const DetectorPeakResponse::EffEval eval = drf->efficiencyEval( 1332.49f, 25.0*PhysicalUnits::cm );
      if( (eval.value > 0.0) && (eval.sigma > 0.0) )
        rel_eff_uncert_field = SpecUtils::printCompact( rel_eff * eval.sigma / eval.value, 3 ) + "%";
    }
  }catch( std::exception & )
  {
  }

  std::chrono::minutes tzOffset( 0 );
  if( wApp )
    tzOffset = wApp->environment().timeZoneOffset();
  auto localtime = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
  localtime += tzOffset;

  out << "# Detector Response Function generated by InterSpec "
  << SpecUtils::to_iso_string( localtime ) << (tzOffset.count() ? "" : " UTC")
  << endline << endline;
  
  out << "# Name,";
  output_csv_field( out, drfname );
  out << endline;
  out << "# Description,";
  output_csv_field( out, drfdescription );

  // The coefficients are floats: nine significant digits read back as the same float, where six
  //  could move a many-term keV equation by tenths of a percent, and leave the covariance matrix
  //  not positive definite.
  const std::streamsize prev_precision = out.precision( std::numeric_limits<float>::max_digits10 );

  out << endline << endline
  << "# Intrinsic Efficiency Coefficients for equation of form"
  " Eff(x) = exp( C_0 + C_1*log(x) + C_2*log(x)^2 + ...) where x is energy in "
  << (effInMeV ? "MeV" : "keV") << endline
  << "#  i.e. equation for probability of gamma that hits the face of the detector being detected in the full energy photopeak." << endline
  << "#  Relative Eff is by IEEE 325: Co-60 1332 keV FEP efficiency 25 cm from the face, relative to 1.2E-3 for a 3x3 NaI"
  << (rel_eff_flat_disk ? " (estimated with the flat-disk model, which ignores the interaction depth)." : ".") << endline
  << "# Name,Relative Eff,eff.c name,c0,c1,c2,c3,c4,c5,c6,c7,p0,p1,p2,Calib Distance,Radius (cm),G factor,GeometryType" << endline
  << drfname << " Intrinsic," << rel_eff_field << ",";
  for( size_t i = 0; i < effEqnCoefs.size(); ++i )
    out << "," << effEqnCoefs[i];
  for( size_t i = effEqnCoefs.size(); i < 12; ++i )
    out << ",";
  out << "0.0," << (0.5*diam/PhysicalUnits::cm) << ",0.5,";
  switch( geom_type )
  {
    case DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic:
      out << "FarField";
      break;
      
    case DetectorPeakResponse::EffGeometryType::FarFieldAbsolute:
      out << "FarFieldAbsolute";
      break;
      
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      out << "FixedTotalAct";
      break;
      
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      out << "FixedPerCm2";
      break;
      
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      out << "FixedPerM2";
      break;
      
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      out << "FixedPerGram";
      break;
  }//switch( geom_type )
  out << endline;
  
  out << "# 1 sigma Uncertainties," << rel_eff_uncert_field << ",";
  for( size_t i = 0; i < effEqnCoefsUncerts.size(); ++i )
    out << "," << effEqnCoefsUncerts[i];
  if( effChi2 > 0 )
    out << endline << "# Chi2 / DOF = " << effChi2 << " / " << effDof
        << " = " << (effDof >= 1 ? (effChi2/effDof) : 0.0);
  out << endline;

  // The full coefficient covariance (the diagonal above is its square root)
  const size_t ncoef = effEqnCoefs.size();
  if( m_effEqnCov.size() == ncoef*ncoef )
  {
    out << "# Coefficient covariance matrix (row-major, c0..c" << (ncoef-1) << ")";
    if( m_effBirgeScale > 1.0 )
      out << ", inflated by chi2/dof = " << m_effBirgeScale;
    out << endline;
    for( size_t i = 0; i < ncoef; ++i )
    {
      out << "# Cov row " << i;
      for( size_t j = 0; j < ncoef; ++j )
        out << "," << m_effEqnCov[i*ncoef + j];
      out << endline;
    }
  }//if( have covariance )
  out << endline;
  
  if( (geom_type != DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic)
     && (geom_type != DetectorPeakResponse::EffGeometryType::FarFieldAbsolute) )
  {
    out << "# This is a DRF for a fixed geometry." << endline << endline;
  }else
  {
    //Convert equation into an absolute efficiency at 25 cm, and output those equations
    const double solidAngleAt25cm = DetectorPeakResponse::fractionalSolidAngle( diam, 25*PhysicalUnits::cm + setback );
    out << "# Absolute Efficiency Coefficients (i.e. probability of gamma emitted from source at 25cm being detected in the full energy photopeak) for equation of form"
    " Eff(x) = exp( C_0 + C_1*log(x) + C_2*log(x)^2 + ...) where x is energy in "
    << (effInMeV ? "MeV" : "keV") << " and at distance of 25 cm" << endline
    << "#  i.e. equation for probability of gamma emitted from source at 25cm being detected in the full energy photopeak." << endline
    << "# Name,Relative Eff,eff.c name,c0,c1,c2,c3,c4,c5,c6,c7,p0,p1,p2,Calib Distance,Radius (cm),G factor" << endline
    << drfname << " Absolute," << rel_eff_field << ",";
    for( size_t i = 0; i < effEqnCoefs.size(); ++i )
      out << "," << ( (i==0 ? log(solidAngleAt25cm) : 0.0) + effEqnCoefs[i]);  //todo: make sure its not
    for( size_t i = effEqnCoefs.size(); i < 12; ++i )
      out << ",";
    out << "25," << (0.5*diam/PhysicalUnits::cm) << "," << solidAngleAt25cm << endline;
    out << "# 1 sigma Uncertainties," << rel_eff_uncert_field << ",";
    for( size_t i = 0; i < effEqnCoefsUncerts.size(); ++i )
      out << "," << effEqnCoefsUncerts[i];
    if( effChi2 > 0 )
      out << endline << "# Chi2 / DOF = " << effChi2 << " / " << effDof
      << " = " << (effDof >= 1 ? (effChi2/effDof) : 0.0);
    out << endline << endline;
  }//if( !is_fixed_geometry )
  
  
  //Then need to give FWHM equation form and coefficient, if available.
  if( fwhmCoefs.size() )
  {
    out << "# Full width half maximum (FWHM) follows equation: ";
    switch( resFcnForm )
    {
      case DetectorPeakResponse::kGadrasResolutionFcn:
      {
        out << "GadrasEqn" << endline
            << "# P6 -> resolution offset: FWHM @ E=0 (keV) when positive; when negative, bends the low-energy power law" << endline
            << "# P7 -> % FWHM @ 661 keV" << endline
            << "# P8 -> resolution power" << endline
            << "# if( energy > 661 or P6 = 0 ) FWHM = 6.61*P7*(energy/661)^P8" << endline
            << "# else if( P6 > 0 ) FWHM = sqrt( (P6*(661 - energy)/661)^2 + (6.61*P7*(energy/661)^P8)^2 )" << endline
            << "# else FWHM = 6.61*P7*(max(energy,30)/661)^(P8^(1/ln(1 - P6)))" << endline
            << "# Energy in keV"<< endline
            << "# ,P6,P7,P8" << endline;
        break;
      }//case DetectorPeakResponse::kGadrasResolutionFcn:
        
      case DetectorPeakResponse::kSqrtPolynomial:
      {
        out << "FWHM = sqrt( A0 + A1*(energy/1000) + A2*(energy/1000)^2 + A3*(energy/1000)^3 + ...)" << endline
            << "# Energy in keV" << endline
            << "# ";
        for( size_t i = 0; i < fwhmCoefs.size(); ++i )
          out << ",A" << i;
        out << endline;
        break;
      }//case DetectorPeakResponse::kSqrtPolynomial:
        
      case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      {
        out << "FWHM = sqrt( A0 + A1*energy + A2/energy )" << endline
        << "# Energy in keV" << endline
        << "# ";
        for( size_t i = 0; i < fwhmCoefs.size(); ++i )
          out << ",A" << i;
        out << endline;
        
        break;
      }//case DetectorPeakResponse::kSqrtEnergyPlusInverse:
        
      case DetectorPeakResponse::kConstantPlusSqrtEnergy:
      {
        out << "FWHM = A0 + A1*sqrt( energy )" << endline
        << "# Energy in keV" << endline
        << "# ";
        for( size_t i = 0; i < fwhmCoefs.size(); ++i )
          out << ",A" << i;
        out << endline;
        
        break;
      }//case DetectorPeakResponse::kConstantPlusSqrtEnergy:
        
      case DetectorPeakResponse::kNumResolutionFnctForm:
        assert( 0 );
        break;
    }//switch( resFcnForm )
    
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
      out << endline << "# Chi2 / DOF = " << fwhmChi2 << " / " << fwhmDof
          << " = " << (fwhmDof >= 1 ? (fwhmChi2/fwhmDof) : 0.0);
  }//if( fwhmCoefs.size() )
  
  out.precision( prev_precision );

  // The energy range of the points the efficiency was fit to - what the DRF itself stores (`data`
  //  also holds peaks that gave no efficiency).
  const float lowerEnergy = m_effLowerEnergy, upperEnergy = m_effUpperEnergy;
  
  out << endline << endline;
  switch( geom_type )
  {
    case DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic:
    case DetectorPeakResponse::EffGeometryType::FarFieldAbsolute:
    {
      out << "Detector diameter = " << (diam/PhysicalUnits::cm) << " cm." << endline;
      if( setback > 0.0 )
        out << "Detector setback = " << (setback/PhysicalUnits::cm) << " cm." << endline;
      try
      {
        const MakeDrfCalc::GeometryChoice geom = geometryChoice();
        if( geom.geometry )
        {
          string xml = geom.geometry->to_xml_string();
          SpecUtils::ireplace_all( xml, "\r", "" );
          SpecUtils::ireplace_all( xml, "\n", "" );
          out << "# Detector geometry (CeeLo XML): " << xml << endline;
        }
      }catch( std::exception & )
      {
      }
      break;
    }
      
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      out << "Fixed Geometry DRF - total activity, # No detector dimensions or source distances recorded." << endline;
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      out << "Fixed Geometry DRF - activity per cm2, # No detector dimensions or source distances recorded." << endline;
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      out << "Fixed Geometry DRF - activity per m2, # No detector dimensions or source distances recorded." << endline;
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      out << "Fixed Geometry DRF - activity per gram, # No detector dimensions or source distances recorded." << endline;
      break;
  }//switch( geom_type )
  
  out << "Valid energy range: " << lowerEnergy
      << " keV to " << upperEnergy << " keV." << endline
      << endline;
  
  //Give infomation on the peaks and sources used to do the fit
  out << "# Peaks used to create DRF" << endline;
  
  out << "# Energy (keV),LiveTime (s),PeakArea,PeakArea Uncert,FWHM (keV),FWHM Uncert (keV),"
         "Source CPS,Source CPS Uncert,Distance (cm),SourceInfo,BackgroundPeakCounts,BackgroundLiveTime"
         ",DetectionEfficiency,DetectionEfficiencyUncert,FlatDiskGeometryFactor"
         ",DistanceUncert (cm),SourceKey,StatUncert (frac),CertUncert (frac),File,Samples" << endline;
  const vector<MeasuredEffPoint> &meas_pts = m_measuredPoints.points();
  for( MakeDrfChart::DataPoint d: data )
  {
    SpecUtils::ireplace_all( d.source_information, ",", " ");
    
    const double peakCps = d.peak_area / d.livetime;
    const double deteff = peakCps / d.source_count_rate;
    const double deteffUncert = deteff * sqrt( pow(d.peak_area_uncertainty/d.peak_area,2)
                                               + pow(d.source_count_rate_uncertainty/d.source_count_rate,2) );
    const double geomFactor = (d.distance < 0.0) ? 1.0 : DetectorPeakResponse::fractionalSolidAngle( diam, d.distance + setback );

    // The raw point this row came from (same energy and peak area)
    const MeasuredEffPoint *raw = nullptr;
    for( const MeasuredEffPoint &p : meas_pts )
    {
      if( (fabs(p.energy - d.energy) < 0.01f) && (fabs(p.peakArea - d.peak_area) < 0.01f*fabs(d.peak_area) + 1.0E-3f) )
      {
        raw = &p;
        break;
      }
    }

    out << d.energy << "," << d.livetime
        << "," << d.peak_area << "," << d.peak_area_uncertainty
        << "," << d.peak_fwhm << "," << d.peak_fwhm_uncertainty
        << "," << d.source_count_rate << "," << d.source_count_rate_uncertainty
        << "," << ((d.distance < 0.0) ? 0.0 : (d.distance / PhysicalUnits::cm))
        << "," << d.source_information
        << "," << (d.background_peak_area > 0.0 ? std::to_string(d.background_peak_area) : string("") )
        << "," << (d.background_peak_live_time > 0.0 ? std::to_string(d.background_peak_live_time) : string("") )
        << "," << deteff
        << "," << deteffUncert
        << "," << geomFactor
        << "," << (raw ? (raw->distanceUncert / PhysicalUnits::cm) : 0.0)
        << ",";
    if( raw )
      output_csv_field( out, raw->sourceKey );
    out << "," << (raw ? raw->fracStatUncert : 0.0f)
        << "," << (raw ? raw->fracCertUncert : 0.0f)
        << ",";
    if( raw )
      output_csv_field( out, raw->fileName );
    out << ",";
    if( raw )
      output_csv_field( out, raw->sampleNumbers );  //e.g. "1-3,5"
    out << endline;
  }//for( MakeDrfChart::DataPoint d: data )
  
  out << endline << endline;

  // The calibration sources the points came from
  if( !m_measuredPoints.sources().empty() )
  {
    out << "# Sources used to create DRF" << endline;
    out << "# SourceKey,Nuclide,Activity (Bq),ActivityUncert (%),Age (s),Distance (cm),DistanceUncert (cm)"
           ",ShieldAtomicNumber,ShieldArealDensity (g/cm2),ShieldMaterial,AssayInfo" << endline;
    for( const MeasuredSourceInfo &s : m_measuredPoints.sources() )
    {
      out << s.sourceKey << "," << s.nuclide
          << "," << (s.activity / PhysicalUnits::bq)
          << "," << (100.0f * s.fracActivityUncert)
          << "," << ((s.age >= 0.0) ? std::to_string( s.age / PhysicalUnits::second ) : string())
          << "," << ((s.distance >= 0.0f) ? std::to_string( s.distance / PhysicalUnits::cm ) : string())
          << "," << ((s.distanceUncert > 0.0f) ? std::to_string( s.distanceUncert / PhysicalUnits::cm ) : string())
          << "," << s.shieldAtomicNumber
          << "," << (s.shieldArealDensity / (PhysicalUnits::g/PhysicalUnits::cm2))
          << ",";
      output_csv_field( out, s.shieldMaterial );
      out << ",";
      output_csv_field( out, s.assayInfo );
      out << endline;
    }//for( const MeasuredSourceInfo &s : m_measuredPoints.sources() )
    out << endline << endline;
  }//if( !m_measuredPoints.sources().empty() )
  
  //Then should give GADRAS Efficiency.csv style efficiecny and FWHM every 50 keV
  out << "# Energy (keV),IntrinsicEfficiency,IntrinsicEfficiency Uncert,FWHM" << endline;
  for( int i = 25; i <= 3000; i += 25 )
  {
    const float energy = i / (effInMeV ? 1000.0f : 1.0f);
    const float eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, effEqnCoefs );
    const float effUncert = (effEqnCoefs.size() == effEqnCoefsUncerts.size())
                            ? effEqnUncert( energy, effEqnCoefs, effEqnCoefsUncerts, m_effEqnCov )
                            : -1.0f;
    
    out << i << "," << 100.0*eff << ",";
    
    if( effUncert >= 0.0f )
      out << 100.0*effUncert;
    out << ",";
    
    if( !fwhmCoefs.empty() )
      out << DetectorPeakResponse::peakResolutionFWHM( i, resFcnForm, fwhmCoefs );
    out << endline;
  }//for( int i = 25; i <= 3000; i += 25 )
  
}//void writeCsvSummary( std::ostream & )




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
  
  
  DetectorPeakResponse::EffGeometryType geom_type = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
  switch( m_geometry->currentIndex() )
  {
    case 0: geom_type = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;   break;
    case 1: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct;   break;
    case 2: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2;  break;
    case 3: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2;   break;
    case 4: geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram; break;
    default: assert( 0 ); break;
  }//switch( m_geometry->currentIndex() )
  
  const double diam = [this]() -> double {try{return detectorDiameter();}catch(...){return 0.0;}}();
  const double setback = detectorSetback();
  const WString diameter = ((geom_type == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic) || (diam > 0.0))
                              ? WString::fromUTF8( PhysicalUnits::printToBestLengthUnits(diam) ) : WString("N/A");

  // One line on how the detector was described: its full geometry, or just the flat disk
  string geometry_txt;
  try
  {
    const MakeDrfCalc::GeometryChoice geom = geometryChoice();
    if( geom.geometry )
    {
      const ceelo::GeometryDescriptor &gd = *geom.geometry;
      char buffer[128];
      if( gd.shape == ceelo::DetectorShape::Cylinder )
      {
        const ceelo::CylinderDims dims = gd.cylinder_dims();
        snprintf( buffer, sizeof(buffer), "%.3g cm dia. x %.3g cm long cylinder", 2.0*dims.radius_cm, dims.full_length_cm );
      }else
      {
        snprintf( buffer, sizeof(buffer), "%.3g x %.3g cm face rectangular crystal",
                  2.0*gd.dimensions_cm[0], 2.0*gd.dimensions_cm[1] );
      }
      geometry_txt = string("Geometry: ") + buffer;
      if( (gd.crystal_material_index >= 0) && (gd.crystal_material_index < static_cast<int>(gd.materials.size())) )
        geometry_txt += ", " + gd.materials[gd.crystal_material_index].name;
      if( !gd.layers.empty() )
        geometry_txt += ", " + std::to_string( gd.layers.size() ) + " housing layer" + (gd.layers.size() > 1 ? "s" : "");
    }else if( setback > 0.0 )
    {
      geometry_txt = "Flat disk model, crystal setback " + PhysicalUnits::printToBestLengthUnits(setback);
    }else
    {
      geometry_txt = "Flat disk model";
    }
  }catch( std::exception & )
  {
  }
  
  const int tzOffsetSecs = static_cast<int>( wApp->environment().timeZoneOffset().count() ) * 60;
  const WDateTime now = WDateTime::currentDateTime().addSecs( tzOffsetSecs );
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
  const auto fwhmForm = DetectorPeakResponse::ResolutionFnctForm( m_fwhmEqnType->currentIndex() );

  // The DRF being saved: its own efficiency (with a full geometry, its CeeLo response) gives the
  //  relative efficiency and the distance rows of the table, and it is what the QR code encodes.
  shared_ptr<const DetectorPeakResponse> drf;
  string drf_error;
  try
  {
    drf = assembleDrf( drfname, drfdescrip );
  }catch( std::exception &e )
  {
    drf_error = e.what();
  }

  // IEEE 325 relative efficiency (what vendor data sheets quote) - none for a fixed geometry, or a
  //  fit that does not reach 1332 keV - and the resolution at 661 keV.  Without a CeeLo response it
  //  is the flat-disk scaling of the intrinsic curve to 25 cm, which ignores the interaction depth,
  //  so is marked approximate.
  string rel_eff_txt;
  const double relEffPercent = drf ? drf->relativeEfficiencyIeee325() : -1.0;
  if( relEffPercent > 0.0 )
  {
    const char * const msg_id = drf->ceeloResponse() ? "md-html-rel-eff" : "md-html-rel-eff-approx";
    rel_eff_txt = WString::tr(msg_id).arg( SpecUtils::printCompact(relEffPercent, 3) ).toUTF8();
  }

  if( !m_fwhmCoefs.empty() )
  {
    const float fwhm661 = DetectorPeakResponse::peakResolutionFWHM( 661.7f, fwhmForm, m_fwhmCoefs );
    const string fwhm_txt = SpecUtils::printCompact( 100.0*fwhm661/661.7, 2 );
    rel_eff_txt += (rel_eff_txt.empty() ? "" : ", ") + WString::tr("md-html-fwhm-661").arg( fwhm_txt ).toUTF8();
  }//if( we have FWHM )

  string fwhm_eqn;
  if( m_fwhmCoefs.size() )
  {
    switch( fwhmForm )
    {
      case DetectorPeakResponse::kGadrasResolutionFcn:
      {
        // Mirrors DetectorPeakResponse::peakResolutionFWHM (GADRAS "form C"): a power law above
        //  661 keV; below it, the sign of P6 picks the branch.
        const float P6 = m_fwhmCoefs[0];
        const float P7 = (m_fwhmCoefs.size() > 1) ? m_fwhmCoefs[1] : 0.0f;
        const float P8 = (m_fwhmCoefs.size() > 2) ? m_fwhmCoefs[2] : 0.0f;

        char high_buffer[256] = { '\0' }, low_buffer[512] = { '\0' };
        snprintf( high_buffer, sizeof(high_buffer), "FWHM(keV) = 6.61 * %.5g * (x/661)<sup>%.5g</sup>", P7, P8 );

        if( P6 == 0.0f )
        {
          fwhm_eqn = high_buffer;
        }else
        {
          if( P6 > 0.0f )
          {
            snprintf( low_buffer, sizeof(low_buffer),
                      "FWHM(keV) = sqrt( (%.5g * (661 - x)/661)<sup>2</sup>"
                      " + (6.61 * %.5g * (x/661)<sup>%.5g</sup>)<sup>2</sup> )", P6, P7, P8 );
          }else
          {
            // The low-energy power is a constant, P8^(1/ln(1 - P6)); give its value
            const double power = std::pow( static_cast<double>(P8), 1.0/std::log(1.0 - P6) );
            snprintf( low_buffer, sizeof(low_buffer),
                      "FWHM(keV) = 6.61 * %.5g * (max(x,30)/661)<sup>%.5g</sup>", P7, power );
          }

          fwhm_eqn = "<div class=\"gadras-fcn\">x &#8804; 661: " + string(low_buffer) + "</div>"
                     "<div class=\"gadras-fcn\">x &gt; 661: " + string(high_buffer) + "</div>";
        }//if( P6 == 0 ) / else
        break;
      }//case DetectorPeakResponse::kGadrasResolutionFcn:
        
      case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      {
        fwhm_eqn = "FWHM(keV) = sqrt(";
        for( size_t i = 0; i < m_fwhmCoefs.size(); ++i )
        {
          if( i == 0 )
            fwhm_eqn += ((m_fwhmCoefs[i] < 0.0) ? " -" : "");
          else
            fwhm_eqn += ((m_fwhmCoefs[i] < 0.0) ? " - " : " + ");
          
          char buffer[64] = { '\0' };
          snprintf( buffer, sizeof(buffer), "%.4g", fabs(m_fwhmCoefs[i]) );
          fwhm_eqn += buffer;
          
          if( i == 1 )
            fwhm_eqn += "*x";
          else if( i == 2 )
            fwhm_eqn += "/x";
        }
        fwhm_eqn += ")";
        break;
      }//case DetectorPeakResponse::kSqrtEnergyPlusInverse:
        
      case DetectorPeakResponse::kConstantPlusSqrtEnergy:
      {
        assert( m_fwhmCoefs.size() == 2 );
        const float A0 = (m_fwhmCoefs.size()) > 0 ? m_fwhmCoefs[0] : 0.0f;
        const float A1 = (m_fwhmCoefs.size()) > 1 ? m_fwhmCoefs[1] : 0.0f;
        
        char buffer[128] = { '\0' };
        snprintf( buffer, sizeof(buffer), "FWHM(keV) = %.4g + %.4g*sqrt(x)", A0, A1 );
        fwhm_eqn = buffer;
        break;
      }//case DetectorPeakResponse::kConstantPlusSqrtEnergy:
        
      case DetectorPeakResponse::kSqrtPolynomial:
      {
        fwhm_eqn = "FWHM(keV) = sqrt(";
        for( size_t i = 0; i < m_fwhmCoefs.size(); ++i )
        {
          if( i == 0 )
            fwhm_eqn += ((m_fwhmCoefs[i] < 0.0) ? " -" : "");
          else
            fwhm_eqn += ((m_fwhmCoefs[i] < 0.0) ? " - " : " + ");
          
          // The polynomial is in MeV: FWHM = sqrt( A0 + A1*(x/1000) + A2*(x/1000)^2 + ... )
          char buffer[64] = { '\0' };
          snprintf( buffer, sizeof(buffer), "%.4g", fabs(m_fwhmCoefs[i]) );
          fwhm_eqn += buffer;
          
          if( i == 1 )
            fwhm_eqn += "*(x/1000)";
          else if( i > 1 )
            fwhm_eqn += "*(x/1000)<sup>" + std::to_string(i) + "</sup>";
        }
        fwhm_eqn += ")";
        
        break;
      }//case DetectorPeakResponse::kSqrtPolynomial:
        
        
      case DetectorPeakResponse::kNumResolutionFnctForm:
        assert( 0 );
        break;
    }//switch( fwhmForm )
  }//if( m_fwhmCoefs.size() )
  
  stringstream efftable;
  const double eff_energies[] = { 59.5, 185.7, 413.7, 661.7, 1001.0, 1460.8, 2614.5 };
  const int eff_dist_cm[] = { 25, 50, 75, 100, 250 };
  
  efftable << "<table class=\"eff-table\">\n"
              "\t\t<caption>" << WString::tr("md-html-eff-table-caption").toUTF8() << "</caption>\n"
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
  
  switch( geom_type )
  {
    case DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic:
      efftable << "\t<tr><th>" << WString::tr("md-html-Intrinsic").toUTF8() << "</th>";
      break;
    case DetectorPeakResponse::EffGeometryType::FarFieldAbsolute:
      efftable << "\t<tr><th>" << WString::tr("md-html-Intrinsic").toUTF8() << "</th>";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      efftable << "\t<tr><th>" << WString::tr("Efficiency").toUTF8() << "</th>";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      efftable << "\t<tr><th>" << WString::tr("Efficiency").toUTF8() << "/cm2</th>";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      efftable << "\t<tr><th>" << WString::tr("Efficiency").toUTF8() << "/m2</th>";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      efftable << "\t<tr><th>" << WString::tr("Efficiency").toUTF8() << "/gram</th>";
      break;
  }//switch( geom_type )
  
 
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
  
  if( (geom_type == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic)
     || (geom_type == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute) )
  {
    assert( geom_type != DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
    
    // The DRF's own answer at each distance - with a full geometry that includes the near field,
    //  which the flat-disk product below (the fallback) does not.
    for( const int dist_cm : eff_dist_cm )
    {
      efftable << "\t<tr><th>" << dist_cm << " cm</th>";
      for( double raw_energy : eff_energies )
      {
        if( (raw_energy < (m_effLowerEnergy - 10)) || (raw_energy > (m_effUpperEnergy + 10)) )
          continue;
        
        const float energy = raw_energy / (effInMeV ? 1000.0f : 1.0f);
        double abs_eff = -1.0;
        try
        {
          if( drf && drf->isValid() )
            abs_eff = drf->efficiency( static_cast<float>(raw_energy), dist_cm*PhysicalUnits::cm );
        }catch( std::exception & )
        {
        }

        if( abs_eff < 0.0 )
        {
          const float intrinsic_eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, m_effEqnCoefs );
          abs_eff = intrinsic_eff * DetectorPeakResponse::fractionalSolidAngle( diam, dist_cm*PhysicalUnits::cm + setback );
        }

        snprintf( buffer, sizeof(buffer), "<td>%.2G</td>", 100.0*abs_eff );
        efftable << buffer;
      }
      efftable << "\t</tr>\n";
    }//for( const auto dist : eff_dist )
  }//if( !is_fixed_geometry )
  
  efftable << "\t</tbody>\n"
              "</table>\n";
  
  WLength chart_width(5*96,WLength::Unit::Pixel), chart_height(3*96,WLength::Unit::Pixel);
  
  Wt::WSvgImage eff_chart( chart_width, chart_height, false );

  {// begin make DRF chart
    MakeDrfChart chart;
    chart.resize( chart_width, chart_height );
    std::vector<MakeDrfChart::DataPoint> data = m_dataPoints;

    // With a full geometry, the curve was fit to points divided by the geometry kernel's factor
    //  (MakeDrfCalc::intrinsicFitPoints), not the flat disk; draw them the same way, or near-field
    //  points sit off the curve they were fit to.
    try
    {
      const MakeDrfCalc::GeometryChoice choice = geometryChoice();
      if( !choice.fixedGeometry && choice.geometry )
      {
        CeeLoUtils::GeometryKernel kernel( *choice.geometry, 2048, choice.diameter / PhysicalUnits::cm );
        for( MakeDrfChart::DataPoint &d : data )
        {
          if( (d.distance >= 0.0f) && (d.energy > 0.0f) )
            d.geometric_factor = kernel.intrinsicFactor( d.energy, d.distance / PhysicalUnits::cm );
        }
      }//if( a full geometry )
    }catch( std::exception & )
    {
      //Fall back to the flat disk for any point not done
    }

    const float detDiam = ((geom_type != DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic)
                           && (geom_type != DetectorPeakResponse::EffGeometryType::FarFieldAbsolute))
                          ? -1.0f : static_cast<float>( m_dataDiameter );
    const float chartLower = 10.0f * std::round( 0.1f*(m_effLowerEnergy - 6.0f) );
    const float chartUpper = 10.0f * std::round( 0.1f*(m_effUpperEnergy + 6.0f) );
    
    chart.setDataPoints( data, detDiam, static_cast<float>( setback ), chartLower, chartUpper );
    
    const auto units = effInMeV ? MakeDrfChart::EqnEnergyUnits::MeV
                                : MakeDrfChart::EqnEnergyUnits::keV;
    
    MakeDrfChart::FwhmCoefType fwhmEqnType = MakeDrfChart::FwhmCoefType::Gadras;
    
    switch( DetectorPeakResponse::ResolutionFnctForm(m_fwhmEqnType->currentIndex()) )
    {
      case DetectorPeakResponse::kGadrasResolutionFcn:
        fwhmEqnType = MakeDrfChart::FwhmCoefType::Gadras;
        break;
        
      case DetectorPeakResponse::kSqrtEnergyPlusInverse:
        fwhmEqnType = MakeDrfChart::FwhmCoefType::SqrtEnergyPlusInverse;
        break;
      
      case DetectorPeakResponse::kConstantPlusSqrtEnergy:
        fwhmEqnType = MakeDrfChart::FwhmCoefType::ConstantPlusSqrtEnergy;
        break;
        
      case DetectorPeakResponse::kSqrtPolynomial:
        fwhmEqnType = MakeDrfChart::FwhmCoefType::SqrtEqn;
        break;
      
      case DetectorPeakResponse::kNumResolutionFnctForm:
        assert( 0 );
        break;
    }//switch( m_fwhmEqnType->currentIndex() )
    
    chart.setFwhmCoefficients( m_fwhmCoefs, m_fwhmCoefUncerts, fwhmEqnType, units );
    chart.setEfficiencyCoefficients( m_effEqnCoefs, m_effEqnCoefUncerts, m_effEqnCov, units );
    
    chart.update();
    
    WPainter p( &eff_chart );
    chart.paint( p );
    p.end();
  }// end make DRF chart
  
  stringstream eff_chart_svg;
  eff_chart.write( eff_chart_svg );
  
  string qr_code, qr_note;
  try
  {
    if( !drf )
      throw runtime_error( drf_error );
    
    // The URL carries the detector's shape, but not a Monte-Carlo response - that is orders of
    //  magnitude past what a QR code holds.  Say so on the sheet.
    if( drf->ceeloResponse() )
      qr_note = WString::tr("md-qr-no-mc-note").toUTF8();
    
    const string url = drf->toAppUrlQr();
    
    tuple<std::string,int,QrCode::ErrorCorrLevel> qr_and_size
              = QrCode::utf8_string_to_svg_qr( url, QrCode::ErrorCorrLevel::About30Percent, 5 );
    
    qr_code = get<0>(qr_and_size);
  }catch(std::exception &e )
  {
    qr_code = WString::tr("app-qr-err").arg(e.what()).toUTF8();
    cerr << qr_code << endl;
  }
  
  
  WTemplate tmplt;
  // tmplttxt is unsafe HTML because of <head> tag and stuff, so we dont want Wt to filter it
  tmplt.setTemplateText( WString::fromUTF8(tmplttxt), TextFormat::UnsafeXHTML );
  // User-entered text goes into the HTML as text, not markup
  const string drfname_html = Wt::Utils::htmlEncode( drfname );
  const string drfdescrip_html = Wt::Utils::htmlEncode( drfdescrip );
  tmplt.bindString("drf-name", drfname_html, TextFormat::UnsafeXHTML );
  tmplt.bindString("drf-desc", drfdescrip_html, TextFormat::UnsafeXHTML );
  tmplt.bindString("diameter", diameter, TextFormat::UnsafeXHTML );
  tmplt.bindString("geometry", Wt::Utils::htmlEncode(geometry_txt), TextFormat::UnsafeXHTML );
  tmplt.bindString("date", date, TextFormat::UnsafeXHTML );
  tmplt.bindString("rel-eff", rel_eff_txt, TextFormat::UnsafeXHTML );
  tmplt.bindString("eff-eqn", eff_eqn, TextFormat::UnsafeXHTML );
  tmplt.bindString("fwhm-eqn", fwhm_eqn, TextFormat::UnsafeXHTML );
  tmplt.bindString("eff-table", efftable.str(), TextFormat::UnsafeXHTML );
  tmplt.bindString("eff-svg", eff_chart_svg.str(), TextFormat::UnsafeXHTML );
  tmplt.bindString("qr-code-title", drfname_html, TextFormat::UnsafeXHTML );
  tmplt.bindString("qr-code", qr_code, TextFormat::UnsafeXHTML );
  tmplt.bindString("qr-code-summary", drfdescrip_html, TextFormat::UnsafeXHTML );
  tmplt.bindString("qr-code-note", qr_note, TextFormat::UnsafeXHTML );
  
  
  tmplt.renderTemplate( output );
}//void writeRefSheet(...)


bool MakeDrf::isEffEqnInMeV() const
{
  return (m_effEqnUnits->currentIndex() == 1);
}//bool isEffEqnInMeV() const


double MakeDrf::detectorDiameter() const
{
  const double diam = geometryChoice().diameter;
  if( diam <= 0.0 )
    throw runtime_error( "Detector diameter less than or equal to zero." );
  return diam;
}//double detectorDiameter() const


double MakeDrf::detectorSetback() const
{
  // In geometry mode the entered setback is a hidden, stale field; the geometry places the crystal.
  if( m_geomMode && m_geomMode->isChecked() )
  {
    try
    {
      return geometryChoice().setback;
    }catch( std::exception & )
    {
      return 0.0;
    }
  }//if( geometry mode )

  return enteredSetback();
}//double detectorSetback() const


double MakeDrf::enteredSetback() const
{
  try
  {
    const double sb = PhysicalUnits::stringToDistance( m_detSetback->text().toUTF8() );
    return (sb > 0.0) ? sb : 0.0;
  }catch(...)
  {
  }

  return 0.0;
}//double enteredSetback() const


void MakeDrf::useSourceLibrary( const vector<shared_ptr<const SrcLibLineInfo>> &srcs, 
                               const bool auto_populate )
{
  for( auto w : m_files->children() )
  {
    auto f = dynamic_cast<DrfSpecFile *>( w );
    if( f )
    {
      for( DrfSpecFileSample *sample : f->fileSamples() )
        sample->addSourceLibrary( srcs, auto_populate );
    }//if( this is a DrfSpecFile widget )
  }//for( auto w : m_files->children() )
  
}//void useSourceLibrary( const vector<SrcLibLineInfo> &srcs );
