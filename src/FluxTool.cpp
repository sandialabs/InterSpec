#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <vector>
#include <numeric>
#include <assert.h>
#include <algorithm>

#include <Wt/WText>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WResource>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>

#include "InterSpec/FluxTool.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/DetectorEdit.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "SpecUtils/UtilityFunctions.h"

using namespace std;
using namespace Wt;

namespace FluxToolImp
{
  class FluxCsvResource : public Wt::WResource
  {
  protected:
    FluxToolWidget *m_fluxtool;
    
  public:
    FluxCsvResource( FluxToolWidget *parent )
    : WResource( parent ),
      m_fluxtool( parent )
    {
    }
    
    virtual ~FluxCsvResource()
    {
      beingDeleted();
    }
    
    
  private:
    virtual void handleRequest( const Wt::Http::Request &, Wt::Http::Response &response )
    {
      const string eol_char = "\r\n"; //for windows - could potentially cosutomize this for the users operating system
      
      if( !m_fluxtool )
        return;
      
      string filename;
      auto meas = m_fluxtool->m_interspec->measurment(SpectrumType::kForeground);
      if( meas && !meas->filename().empty() )
      {
        filename = UtilityFunctions::filename( meas->filename() );
        const string extension = UtilityFunctions::file_extension( filename );
        filename = filename.substr( 0, filename.size() - extension.size() );
        if( filename.size() )
          filename += "_";
      }
      filename += "flux.csv";
      
      
      suggestFileName( filename, WResource::Attachment );
      
      for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
          col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
      {
        const WString &colname = m_fluxtool->m_colnamesCsv[col];
        
        response.out() << (col==0 ? "" : ",") << colname;
        
        //No uncertainty on energy.
        switch( col )
        {
          case FluxToolWidget::FluxEnergyCol:
          case FluxToolWidget::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxGeometricEffCol:
          case FluxToolWidget::FluxNumColumns:
            break;
            
          case FluxToolWidget::FluxPeakCpsCol:
          case FluxToolWidget::FluxFluxOnDetCol:
          case FluxToolWidget::FluxFluxPerCm2PerSCol:
          case FluxToolWidget::FluxGammasInto4PiCol:
            response.out() << "," << (colname + " Uncertainty");
            break;
        }//switch( col )
      }//for( loop over columns )
      
      response.out() << eol_char;
      
      for( size_t row = 0; row < m_fluxtool->m_data.size(); ++row )
      {
        for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
            col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
        {
          const double data = m_fluxtool->m_data[row][col];
          const double uncert = m_fluxtool->m_uncertainties[row][col];
          
          response.out() << (col==0 ? "" : ",") << std::to_string(data);
          switch( col )
          {
            case FluxToolWidget::FluxEnergyCol:
            case FluxToolWidget::FluxIntrinsicEffCol:
            case FluxToolWidget::FluxGeometricEffCol:
            case FluxToolWidget::FluxNumColumns:
              break;
              
            case FluxToolWidget::FluxPeakCpsCol:
            case FluxToolWidget::FluxFluxOnDetCol:
            case FluxToolWidget::FluxFluxPerCm2PerSCol:
            case FluxToolWidget::FluxGammasInto4PiCol:
              response.out() << ",";
              if( uncert > std::numeric_limits<double>::epsilon() )
                response.out() << std::to_string(uncert);
              break;
          }//switch( col )
        }//for( loop over columns )
        
        response.out() << eol_char;
      }//for( size_t row = 0; row < m_fluxtool->m_data.size(); ++row )
    }//handleRequest(...)
    
  };//class class FluxCsvResource
}//namespace


FluxToolWindow::FluxToolWindow( InterSpec *viewer )
: AuxWindow( "Flux Tool",
  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletModal)
   | AuxWindowProperties::SetCloseable
   | AuxWindowProperties::DisableCollapse) ),
  m_fluxTool( nullptr )
{
  assert( viewer );
  rejectWhenEscapePressed( true );
  
  m_fluxTool = new FluxToolWidget( viewer, contents() );
  m_fluxTool->setHeight( WLength(100, WLength::Percentage) );
  
  AuxWindow::addHelpInFooter( footer(), "flux-tool" );
  
  
  WContainerWidget *buttonDiv = footer();
  
  WResource *csv = new FluxToolImp::FluxCsvResource( m_fluxTool );
#if( BUILD_AS_OSX_APP )
  WAnchor *csvButton = new WAnchor( WLink(csv), buttonDiv );
  csvButton->setTarget( AnchorTarget::TargetNewWindow );
#else
  WPushButton *csvButton = new WPushButton( buttonDiv );
  csvButton->setIcon( "InterSpec_resources/images/download_small.png" );
  csvButton->setLink( WLink(csv) );
  csvButton->setLinkTarget( Wt::TargetNewWindow );
#endif
  
  csvButton->setText( "CSV" );
  csvButton->setStyleClass( "CsvLinkBtn" );
  
  auto enableDisableCsv = [csvButton,this](){
    csvButton->setEnabled( !m_fluxTool->m_data.empty() );
  };
  
  m_fluxTool->m_tableUpdated.connect( std::bind(enableDisableCsv) );
  csvButton->disable();
  
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  
  //const int screenW = viewer->renderedWidth();
  //const int screenH = viewer->renderedHeight();
  //const int width = ((screenW < 600) ? screenW : 600);
  //const int height = ((screenH < 420) ? screenH : 420);
  //resizeWindow( width, height );
  
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//FluxToolWindow(...) constrctor


FluxToolWindow::~FluxToolWindow()
{
}


FluxToolWidget::FluxToolWidget( InterSpec *viewer, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_interspec( viewer ),
    m_detector( nullptr ),
    m_msg( nullptr ),
    m_table( nullptr ),
    m_distance( nullptr ),
    m_needsTableRefresh( true ),
    m_tableUpdated( this )
{
  init();
}


FluxToolWidget::~FluxToolWidget()
{
}


void FluxToolWidget::init()
{
  assert( m_interspec );
  assert( !m_detector );
  
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
  {
    switch( col )
    {
      case FluxEnergyCol:
        m_colnames[col] = WString::fromUTF8("Energy (keV)");
        m_colnamesCsv[col] = WString::fromUTF8("Energy (keV)");
        break;
      case FluxPeakCpsCol:
        m_colnames[col] = WString::fromUTF8("Peak CPS");
        m_colnamesCsv[col] = WString::fromUTF8("Peak CPS");
        break;
      case FluxIntrinsicEffCol:
        m_colnames[col] = WString::fromUTF8("Intr. Eff.");
        m_colnamesCsv[col] = WString::fromUTF8("Intrinsic Efficiency");
        break;
      case FluxGeometricEffCol:
        m_colnames[col] = WString::fromUTF8("Geom. Eff.");
        m_colnamesCsv[col] = WString::fromUTF8("Geometric Efficiency");
        break;
      case FluxFluxOnDetCol:
        m_colnames[col] = WString::fromUTF8("Flux on Det. (&gamma;/s)");
        m_colnamesCsv[col] = WString::fromUTF8("Flux on Detector (gammas/s)");
        break;
      case FluxFluxPerCm2PerSCol:
        m_colnames[col] = WString::fromUTF8("Flux (&gamma;/cm&sup2;/s)");
        m_colnamesCsv[col] = WString::fromUTF8("Flux (gammas/cm2/s)");
        break;
      case FluxGammasInto4PiCol:
        m_colnames[col] = WString::fromUTF8("&gamma;/4&pi;/s");
        m_colnamesCsv[col] = WString::fromUTF8("gammas/4pi/s");
        break;
      case FluxNumColumns:        break;
    }//switch( col )
  }//for( loop over columns )
  
  
  wApp->useStyleSheet( "InterSpec_resources/FluxTool.css" );
  
  const bool showToolTipInstantly = m_interspec ? InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec ) : false;
  
  addStyleClass( "FluxToolWidget" );
  
  WTable *distDetRow = new WTable( this );
  distDetRow->addStyleClass( "FluxDistMsgDetTbl" );
  
  SpectraFileModel *specFileModel = m_interspec->fileManager()->model();
  m_detector = new DetectorDisplay( m_interspec, specFileModel );
  m_detector->addStyleClass( "FluxDet" );
  m_interspec->detectorChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  m_interspec->detectorModified().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  distDetRow->elementAt(0,1)->addWidget( m_detector );
  
  auto distCell = distDetRow->elementAt(0,0);
  distCell->addStyleClass( "FluxDistCell" );
  WLabel *label = new WLabel( "Distance:", distCell );
  label->addStyleClass( "FluxDistLabel" );
  
  m_distance = new WLineEdit( "100 cm", distCell );
  m_distance->addStyleClass( "FluxDistanceEnter" );
  label->setBuddy(m_distance);
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_distance->setValidator( validator );
  HelpSystem::attachToolTipOn( m_distance,
                              "Distance from center of source to face of detector. Number must be"
                              " followed by units; valid units are: meters, m, cm, mm, km, feet,"
                              " ft, ', in, inches, or \".  You may also add multiple distances,"
                              " such as '3ft 4in', or '3.6E-2 m 12 cm' which are equivalent to "
                              " 40inches and 15.6cm respectively.", showToolTipInstantly );
  m_distance->changed().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  m_distance->enterPressed().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  
  
  PeakModel *peakmodel = m_interspec->peakModel();
  assert( peakmodel );
  peakmodel->dataChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsRemoved().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsInserted().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->layoutChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );

  auto msgCell = distDetRow->elementAt(1,0);
  msgCell->setColumnSpan( 2 );
  msgCell->addStyleClass( "FluxMsgCell" );
  m_msg = new WText( "&nbsp;", Wt::XHTMLText, msgCell );
  //m_msg->setInline( false );
  m_msg->addStyleClass( "FluxMsg" );
  
  
  WContainerWidget *tableHolder = new WContainerWidget( this );
  tableHolder->addStyleClass( "FluxTableHolder" );
  
  m_table = new WTable( tableHolder );
  m_table->addStyleClass( "FluxTable" );
}//void init()


void FluxToolWidget::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_needsTableRefresh )
    refreshPeakTable();
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags );


void FluxToolWidget::setTableNeedsUpdating()
{
  m_needsTableRefresh = true;
  scheduleRender();
}//void setTableNeedsUpdating()


void FluxToolWidget::refreshPeakTable()
{
  PeakModel *peakmodel = m_interspec->peakModel();
  
  m_data.clear();
  m_uncertainties.clear();
  
  m_table->clear();
  m_msg->setText( "&nbsp;" );
  
  //ToDo: switch to RowStretchTreeView or Wt::WTreeView using a WAbstractItemDelegate to format the data... sounds like a lot of work
  
  //Set up table header
  m_table->setHeaderCount( 1, Wt::Horizontal );
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
  {
    new WText( m_colnames[col], m_table->elementAt(0,col) );
  }//for( loop over columns )
  
  
  float distance = static_cast<float>( 1.0*PhysicalUnits::meter );
  try
  {
    distance = static_cast<float>( PhysicalUnits::stringToDistance( m_distance->text().toUTF8() ) );
  }catch(...)
  {
    m_msg->setText( "Invalid Distance" );
    m_tableUpdated.emit();
    return;
  }
  
  auto det = m_detector->detector();
  if( !det || !det->isValid() )
  {
    m_msg->setText( "No Detector Response Function Chosen" );
    m_tableUpdated.emit();
    return;
  }
  
  auto spec = m_interspec->measurment(SpectrumType::kForeground);
  if( !spec )
  {
    m_msg->setText( "No foreground spectrum loaded" );
    m_tableUpdated.emit();
    return;
  }
  
  const float live_time = spec->gamma_live_time();
  if( live_time <= 0.0f )
  {
    m_msg->setText( "Invalid foregorund livetime" );
    m_tableUpdated.emit();
    return;
  }
  
  const vector<PeakDef> peaks = peakmodel->peakVec();
  
  const size_t npeaks = peaks.size();
  m_data.resize( npeaks );
  m_uncertainties.resize( npeaks );

  for( int i = 0; i < static_cast<int>(npeaks); ++i )
  {
    m_data[i].fill( 0.0 );
    m_uncertainties[i].fill( 0.0 );
    
    const PeakDef &peak = peaks[i];
    char buffer[128];
    
    const double energy = peak.mean();
    snprintf( buffer, sizeof(buffer), "%.2f", energy );
    new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,0) );
    
    const double amp = peak.peakArea();  //ToDO: make sure this works for non-Gaussian peaks
    const double ampUncert = peak.peakAreaUncert();
    const double cps = amp / live_time;
    const double cpsUncert = ampUncert / live_time;
    if( ampUncert > 0.0 )
      snprintf( buffer, sizeof(buffer), "%.4g &plusmn; %.1f", cps, cpsUncert );
    else
      snprintf( buffer, sizeof(buffer), "%.4g", cps );
    new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,1) );
    
    const double intrisic = det->intrinsicEfficiency(energy);
    snprintf( buffer, sizeof(buffer), "%.4g", intrisic );
    new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,2) );
    
    const double geomEff = det->fractionalSolidAngle( det->detectorDiameter(), distance );
    snprintf( buffer, sizeof(buffer), "%.4g", geomEff );
    new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,3) );
    
    const double totaleff = det->efficiency(energy, distance );
    
    m_data[i][FluxColumns::FluxEnergyCol] = energy;
    m_data[i][FluxColumns::FluxPeakCpsCol] = cps;
    m_uncertainties[i][FluxColumns::FluxPeakCpsCol] = cpsUncert;
    
    m_data[i][FluxColumns::FluxGeometricEffCol] = geomEff;
    m_data[i][FluxColumns::FluxIntrinsicEffCol] = intrisic;
    //ToDo: Check if there is an uncertainty on DRF, and if so include that.
    
    if( totaleff <= 0.0 || intrisic <= 0.0 )
    {
      new WText( "inf", m_table->elementAt(1+i,4) );
      new WText( "inf", m_table->elementAt(1+i,5) );
      new WText( "inf", m_table->elementAt(1+i,6) );
      m_data[i][FluxColumns::FluxFluxOnDetCol]      = std::numeric_limits<double>::infinity();
      m_data[i][FluxColumns::FluxFluxPerCm2PerSCol] = std::numeric_limits<double>::infinity();
      m_data[i][FluxColumns::FluxGammasInto4PiCol]  = std::numeric_limits<double>::infinity();
    }else
    {
      const double fluxOnDet = cps / intrisic;
      const double fluxOnDetUncert = cpsUncert / intrisic;
      if( ampUncert > 0.0 )
        snprintf( buffer, sizeof(buffer), "%.4g &plusmn; %.4g", fluxOnDet, fluxOnDetUncert );
      else
        snprintf( buffer, sizeof(buffer), "%.4g", fluxOnDet );
      new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,4) );
      
      //gammas into 4pi
      const double gammaInto4pi = cps / totaleff;
      const double gammaInto4piUncert = cpsUncert / totaleff;
      
      //Flux in g/cm2/s
      const double flux = gammaInto4pi / (4*M_PI*distance*distance);
      const double fluxUncert = gammaInto4piUncert / (4*M_PI*distance*distance);
      
      if( ampUncert > 0.0 )
        snprintf( buffer, sizeof(buffer), "%.4g &plusmn; %.4g", flux, fluxUncert );
      else
        snprintf( buffer, sizeof(buffer), "%.4g", flux );
      new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,5) );
      
      if( ampUncert > 0.0 )
        snprintf( buffer, sizeof(buffer), "%.4g &plusmn; %.4g", gammaInto4pi, gammaInto4piUncert );
      else
        snprintf( buffer, sizeof(buffer), "%.4g", gammaInto4pi );
      new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,6) );
      
      m_data[i][FluxColumns::FluxFluxOnDetCol] = fluxOnDet;
      m_uncertainties[i][FluxColumns::FluxFluxOnDetCol] = fluxOnDetUncert;
      
      m_data[i][FluxColumns::FluxFluxPerCm2PerSCol] = flux;
      m_uncertainties[i][FluxColumns::FluxFluxPerCm2PerSCol] = fluxUncert;
      
      m_data[i][FluxColumns::FluxGammasInto4PiCol] = gammaInto4pi;
      m_uncertainties[i][FluxColumns::FluxGammasInto4PiCol] = gammaInto4piUncert;
    }//if( eff > 0 ) / else
  }//for( const PeakDef &peak : peaks )
  
  m_tableUpdated.emit();
}//void refreshPeakTable()






