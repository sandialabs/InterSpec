#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <vector>
#include <assert.h>
#include <algorithm>

#include <Wt/WText>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>
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

using namespace std;
using namespace Wt;


FluxToolWindow::FluxToolWindow( InterSpec *viewer )
: AuxWindow( "Flux Tool",
(Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletModal)
 | AuxWindowProperties::SetCloseable
 | AuxWindowProperties::DisableCollapse) ),
  m_fluxTool( nullptr )
{
  rejectWhenEscapePressed( true );
  
  m_fluxTool = new FluxToolWidget( viewer, contents() );
  m_fluxTool->setHeight( WLength(100, WLength::Percentage) );
  
  AuxWindow::addHelpInFooter( footer(), "flux-tool" );
  
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
    m_needsTableRefresh( true )
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
  
  wApp->useStyleSheet( "InterSpec_resources/FluxTool.css" );
  
  const bool showToolTipInstantly = m_interspec ? InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec ) : false;
  
  addStyleClass( "FluxToolWidget" );
  
  WContainerWidget *distdiv = new WContainerWidget( this );
  distdiv->addStyleClass( "FluxDistRow" );
  
  WLabel *label = new WLabel( "Distance:", distdiv );
  label->addStyleClass( "FluxDistLabel" );
  
  m_distance = new WLineEdit( "100 cm", distdiv );
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
  
  
  SpectraFileModel *specFileModel = m_interspec->fileManager()->model();
  m_detector = new DetectorDisplay( m_interspec, specFileModel );
  m_detector->addStyleClass( "FluxDet" );
  m_interspec->detectorChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  m_interspec->detectorModified().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  
  PeakModel *peakmodel = m_interspec->peakModel();
  assert( peakmodel );
  peakmodel->dataChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsRemoved().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsInserted().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->layoutChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );

  addWidget( m_detector );
  
  m_msg = new WText( "&nbsp;", Wt::XHTMLText, this );
  m_msg->setInline( false );
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
  cout << "FluxToolWidget::refreshPeakTable()" << endl;
  PeakModel *peakmodel = m_interspec->peakModel();
  
  m_table->clear();
  m_msg->setText( "&nbsp;" );
  
  //ToDo: switch to RowStretchTreeView or Wt::WTreeView using a WAbstractItemDelegate to format the data... sounds like a lot of work
  
  //Set up table header
  m_table->setHeaderCount( 1, Wt::Horizontal );
  new WText( "Energy (keV)", Wt::XHTMLText, m_table->elementAt(0,0) );
  new WText( "Peak CPS", Wt::XHTMLText, m_table->elementAt(0,1) );
  new WText( "Intr. Eff.", Wt::XHTMLText, m_table->elementAt(0,2) );
  new WText( "Geom. Eff.", Wt::XHTMLText, m_table->elementAt(0,3) );
  new WText( "Flux on Det. (&gamma;/s)", Wt::XHTMLText, m_table->elementAt(0,4) );
  new WText( "&gamma;/4&pi;", Wt::XHTMLText, m_table->elementAt(0,5) );
  new WText( "Flux (&gamma;/cm&sup2;/s)", Wt::XHTMLText, m_table->elementAt(0,6) );
  
  
  
  float distance = static_cast<float>( 1.0*PhysicalUnits::meter );
  try
  {
    distance = static_cast<float>( PhysicalUnits::stringToDistance( m_distance->text().toUTF8() ) );
  }catch(...)
  {
    m_msg->setText( "Invalid Distance" );
    return;
  }
  
  auto det = m_detector->detector();
  if( !det || !det->isValid() )
  {
    m_msg->setText( "No Detector Response Function Chosen" );
    return;
  }
  
  auto spec = m_interspec->measurment(SpectrumType::kForeground);
  if( !spec )
  {
    m_msg->setText( "No foreground spectrum loaded" );
    return;
  }
  
  const float live_time = spec->gamma_live_time();
  if( live_time <= 0.0f )
  {
    m_msg->setText( "Invalid foregorund livetime" );
    return;
  }
  
  const vector<PeakDef> peaks = peakmodel->peakVec();
  for( int i = 0; i < static_cast<int>(peaks.size()); ++i )
  {
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
    
    if( totaleff <= 0.0 || intrisic <= 0.0 )
    {
      new WText( "inf", m_table->elementAt(1+i,4) );
      new WText( "inf", m_table->elementAt(1+i,5) );
      new WText( "inf", m_table->elementAt(1+i,6) );
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
      if( ampUncert > 0.0 )
        snprintf( buffer, sizeof(buffer), "%.4g &plusmn; %.4g", gammaInto4pi, gammaInto4piUncert );
      else
        snprintf( buffer, sizeof(buffer), "%.4g", gammaInto4pi );
      new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,5) );
      
      //Flux in g/cm2/s
      const double flux = gammaInto4pi / (4*M_PI*distance*distance);
      const double fluxUncert = gammaInto4piUncert / (4*M_PI*distance*distance);
      if( ampUncert > 0.0 )
        snprintf( buffer, sizeof(buffer), "%.4g &plusmn; %.4g", flux, fluxUncert );
      else
        snprintf( buffer, sizeof(buffer), "%.4g", flux );
      new WText( buffer, Wt::XHTMLText, m_table->elementAt(1+i,5) );
    }//if( eff > 0 ) / else
    
    
    
  }//for( const PeakDef &peak : peaks )
  
}//void refreshPeakTable()






