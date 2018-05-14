#include "InterSpec_config.h"

#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <boost/foreach.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include <Wt/WText>
#include <Wt/WPoint>
#include <Wt/WTable>
#include <Wt/WContainerWidget>


#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/SpectrumViewer.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/SimpleNuclideAssist.h"
//#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhotopeakLineDisplay.h"

using namespace std;
using namespace Wt;


#define foreach BOOST_FOREACH
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

const std::string SimpleNuclideAssist::sm_dataFileName = "data/gamma_rays_of_interest.xml";

const WColor s_main_nuc_highlight_color( 255, 255, 0, 100 );
const WColor s_contaminant_highlight_color(122, 122, 122, 75);
const WAnimation s_show_animation( WAnimation::SlideInFromTop, WAnimation::Linear, 250 );
const WAnimation s_appear_animation( WAnimation::Fade, WAnimation::Linear, 250 );

namespace
{
  template<size_t N>
  size_t lengthof(const char (&)[N])
  {
    return N - 1;
  }
  
  template<class Ch>
  std::string xml_value_str( const rapidxml::xml_base<Ch> *n )
  {
    return (n && n->value_size()) ? std::string(n->value(),n->value()+n->value_size()) : std::string();
  }
  
  template<class Ch>
  std::string xml_name_str( const rapidxml::xml_base<Ch> *n )
  {
    return (n && n->name_size()) ? std::string(n->name(), n->name() + n->name_size()) : std::string();
  }
}//anonomous namespace for XML utilities


//These macros are a evolution of what is in SpectrumDataStructs
#define XML_VALUE_COMPARE( node, cstr ) (rapidxml::internal::compare( node->value(), node->value_size(), cstr, lengthof(cstr), true) )
#define XML_VALUE_ICOMPARE( node, cstr )(rapidxml::internal::compare( node->value(), node->value_size(), cstr, lengthof(cstr), false))
#define XML_NAME_COMPARE( node, cstr )  (rapidxml::internal::compare( node->name(),  node->name_size(),  cstr, lengthof(cstr), true) )
#define XML_NAME_ICOMPARE( node, cstr ) (rapidxml::internal::compare( node->name(),  node->name_size(),  cstr, lengthof(cstr), false))

#define XML_FIRST_NODE(node,name)(node->first_node(name,lengthof(name),true))
#define XML_FIRST_ATTRIBUTE(node,name)(node->first_attribute(name,lengthof(name),true))
#define XML_NEXT_SAME_SIBLING(node)(node->next_sibling(node->name(),node->name_size(),true))
  
#define XML_FOREACH_DAUGHTER( nodename, parentnode, daughternamestr ) \
  for( const rapidxml::xml_node<char> *nodename = XML_FIRST_NODE(parentnode,daughternamestr); \
       nodename; \
       nodename = XML_NEXT_SAME_SIBLING(nodename) )


namespace
{
  
  template<typename T, typename Container>
  bool contains(const Container & c, const T & value)
  {
    return std::find(boost::begin(c), boost::end(c), value) != boost::end(c);
  }
  
struct GammaInfo
{
  string parent;
  string from;
  float energy;
  float branchingratio;
};//struct GammaInfo


struct SourceInfo
{
  string name;
  string halflife;
  string specificactivity;
  vector<string> comments;
  vector<string> contaminants;
  vector<string> equivalents;
  vector<GammaInfo> gammas;
};//struct SourceInfo

  
vector<SourceInfo> parseXmlFile()
{
  typedef const rapidxml::xml_node<char> xmlnode;
  typedef const rapidxml::xml_attribute<char> xmlattrib;
  
  vector<SourceInfo> infos;
  
  const char *filename = SimpleNuclideAssist::sm_dataFileName.c_str();
  ::rapidxml::file<char> input_file( filename );
  
  rapidxml::xml_document<char> doc;
  doc.parse<rapidxml::parse_non_destructive>( input_file.data() );
  xmlnode *catagories_node = XML_FIRST_NODE((&doc),"catagories");
  
  if( !catagories_node )
    throw runtime_error( SimpleNuclideAssist::sm_dataFileName
                         + " does not have a <catagories> node" );
  
  XML_FOREACH_DAUGHTER( catagory, catagories_node, "catagory" )
  {
    XML_FOREACH_DAUGHTER( source, catagory, "source")
    {
      SourceInfo info;
      info.name = xml_value_str( XML_FIRST_ATTRIBUTE(source, "name") );
      info.halflife = xml_value_str( XML_FIRST_NODE(source,"halflife") );
      info.specificactivity = xml_value_str( XML_FIRST_NODE(source,"specificactivity") );
      
      XML_FOREACH_DAUGHTER( comment, source, "comment" )
        info.comments.push_back( xml_value_str(comment) );
      
      XML_FOREACH_DAUGHTER( contaminant, source, "contaminant" )
        info.contaminants.push_back( xml_value_str(XML_FIRST_ATTRIBUTE(contaminant,"name")) );
      
      XML_FOREACH_DAUGHTER( equivalent, source, "equivalent" )
        info.equivalents.push_back( xml_value_str(XML_FIRST_ATTRIBUTE(equivalent,"name")) );
      
      XML_FOREACH_DAUGHTER( gamma, source, "gamma" )
      {
        GammaInfo thisgamma;
        thisgamma.parent = xml_value_str( XML_FIRST_ATTRIBUTE(gamma, "parent") );
        thisgamma.from = xml_value_str( XML_FIRST_ATTRIBUTE(gamma, "from") );
        const string energystr = xml_value_str( XML_FIRST_ATTRIBUTE(gamma, "energy") );
        const string brstr = xml_value_str( XML_FIRST_ATTRIBUTE(gamma, "branchingratio") );
        
        if( !(stringstream(energystr) >> thisgamma.energy) )
          throw runtime_error( "Could not convert " + energystr
                              + " to an energy for source " + info.name );
        
        if( !(stringstream(brstr) >> thisgamma.branchingratio) )
          thisgamma.branchingratio = -1.0f;
        
        info.gammas.push_back( thisgamma );
      }//foreach( gamma )
      
      infos.push_back( info );
    }//for( loop over source )
  }//for( loop over catagory )
  
  if( infos.empty() )
    throw runtime_error( SimpleNuclideAssist::sm_dataFileName
                         + " did not appear to be a gamma rays of interest file." );
  
  return infos;
}//vector<SourceInfo> parseDataFile()

} //namespace


class SimpleNuclide
  : public WContainerWidget
{
public:
  enum ShieldingLevel{ NoShielding, SomeShielding, MoreShielding };
  
public:
  const SourceInfo m_source;
  SimpleNuclideAssist *m_assist;
  SpectrumViewer *m_viewer;
  WContainerWidget *m_info;
  
  ShieldingLevel m_shieldLevel;
  WText *m_noshield, *m_someshield, *m_moreshield, *m_shielddesc;
public:
  SimpleNuclide( const SourceInfo &source, SimpleNuclideAssist *assist,
                 SpectrumViewer *viewer,
                 WContainerWidget *parent = 0 )
   : WContainerWidget( parent ),
     m_source( source ),
     m_assist( assist ),
     m_viewer( viewer ),
     m_info( 0 ),
     m_shieldLevel( NoShielding ),
     m_noshield( 0 ), m_someshield( 0 ), m_moreshield( 0 ), m_shielddesc( 0 )
  {
    addStyleClass( "SimpleNuclide" );
    
    string namestr = source.name;
    size_t close_paren_pos = namestr.find( ")" );
    if( close_paren_pos != string::npos )
      namestr = namestr.substr( 0, close_paren_pos + 1 );
    WText *name = new WText( namestr );
    name->addStyleClass( "SimpleNuclideNucName" );
    
    if( m_source.equivalents.size() )
    {
      WTable *table = new WTable( this );
      table->setWidth( WLength(100.0,WLength::Percentage) );
      
      WContainerWidget *cell = table->elementAt(0,0);
      table->columnAt(0)->setWidth( WLength(70.0,WLength::Percentage) );
      cell->addWidget( name );
      name->setAttributeValue( "style", "text-align: left; padding-left: 42.86%;" );
      
      string equivstr;
      foreach( string equiv, m_source.equivalents )
      {
        close_paren_pos = equiv.find( ")" );
        if( close_paren_pos != string::npos )
          equiv = equiv.substr( 0, close_paren_pos + 1 );
        equivstr += "<div>" + equiv + "</div>";
      }
      
      cell = table->elementAt(0,1);
      WText *equivtxt = new WText( equivstr, Wt::XHTMLText, cell );
      equivtxt->addStyleClass( "SimpleNuclideEquivs" );
    }else
    {
      addWidget( name );
      name->setInline( false );
      name->setAttributeValue( "style", "text-align: center;" );
    }
    
    
    
    m_info = new WContainerWidget( this );
    
    //TODO: add shieldings into XML to controll if the bellow even gets displayed
    WContainerWidget *shielding = new WContainerWidget( m_info );
    
    WContainerWidget *topline = new WContainerWidget( shielding );
    WText *shieldtxt = new WText( "Shielding:", topline );
    shieldtxt->setStyleClass( "SimpleNuclideShieldTxt" );
    
    m_shielddesc = new WText( "", topline );
    m_shielddesc->setStyleClass( "SimpleNuclideShieldDesc" );
    
    WContainerWidget *buttonholder = new WContainerWidget( shielding );
    buttonholder->setStyleClass( "SimpleNuclideShieldButHold" );
    m_noshield = new WText( "None", buttonholder );
    m_noshield->setStyleClass( "SimpleNuclideShieldButActive" );
    m_noshield->clicked().connect( boost::bind( &SimpleNuclide::changeShielding, this, NoShielding ) );
    
    m_someshield = new WText( "Some", buttonholder );
    m_someshield->setStyleClass( "SimpleNuclideShieldBut" );
    m_someshield->clicked().connect( boost::bind( &SimpleNuclide::changeShielding, this, SomeShielding ) );
    
    m_moreshield = new WText( "More", buttonholder );
    m_moreshield->setStyleClass( "SimpleNuclideShieldBut" );
    m_moreshield->clicked().connect( boost::bind( &SimpleNuclide::changeShielding, this, MoreShielding ) );
    
    const vector<string> &comments = m_source.comments;
    const vector<string> &contaminants = m_source.contaminants;
    
    
    string extraInfoHtml;
      
    if( m_source.halflife.size() )
    {
      string val = "&lambda;<sub>&#189;</sub>=" + m_source.halflife;
      if( m_source.equivalents.size() )
        val = m_source.name + " " + val;
      extraInfoHtml = "<div>" + val + "</div>";
    }//if( m_source.halflife.size() )
    
    
    if( contaminants.size() )
    {
      string contam;
      if( contaminants.size() > 2 )
      {
        contam += "Contaminants: <div style=\"padding-left: 10px;\">";
        for( size_t i = 0; i < contaminants.size(); ++i )
          contam += (i ? ", " : "") + contaminants[i];
        contam += "</div>";
      }else if( contaminants.size() == 2 )
      {
        contam = "Contaminants: " + contaminants[0] + ", " + contaminants[1];
      }else
      {
        contam = "Contaminant: " + m_source.contaminants[0];
      }
      
      extraInfoHtml += "<div class=\"SimpleNuclideContam\">" + contam + "</div>";
    }//if( m_source.contaminants.size() )
    
    foreach( const string &comment, comments )
      extraInfoHtml += "<div class=\"SimpleNuclideComment\">" + comment + "</div>";
  
    if( m_source.specificactivity.size() )
    {
      if( m_source.specificactivity.find( "SA" ) != string::npos )
        extraInfoHtml += "<div class=\"SimpleNuclideComment\">" + m_source.specificactivity + "</div>";
      else if( m_source.equivalents.empty() )
        extraInfoHtml += "<div class=\"SimpleNuclideComment\">SA of " + m_source.name + " is " + m_source.specificactivity + "</div>";
      else
        extraInfoHtml += "<div class=\"SimpleNuclideComment\">SA is " + m_source.specificactivity + "</div>";
    }
    
    if( extraInfoHtml.empty() )
      extraInfoHtml = "None Available";
      
    WText *moreInfoButton = new WText( "More Info", m_info );
    moreInfoButton->setInline( false );
    moreInfoButton->addStyleClass( "SimpleNuclideMoreInfoButton" );
    
    
    {//begin codeblock to make qtip
      const std::string showfunction
         = "function(event){ $('#" + moreInfoButton->id() + "').qtip({ "
          "overwrite: false,"
          "content:   { text: '" + extraInfoHtml + "'}, "
          "position:  { my: 'left center', at: 'right center', viewport: $(window), adjust: {method: 'flipinvert flipinvert', x:5} },"
          "show:      { ready: true, event: event.type },"
          "events:    { hide: function(){$(this).qtip('destroy');} }, "
          "style:     { classes: 'qtip-rounded qtip-shadow qtip-grey', tip: {corner: true, width: 18, height: 12} }"
        "}, event); "
      "}";
        
//      const std::string closefunction = "function(){$('#" + moreInfoButton->id() + "').qtip('destroy');}";
//      moreInfoButton->clicked().connect( showfunction );
//      moreInfoButton->mouseWentOver().connect( showfunction );
      
        const string showjs = "$(document).on('click mouseover', "
                              "'#" + moreInfoButton->id() + "', "
                              + showfunction + " );";
        moreInfoButton->doJavaScript( showjs );
    }//end codeblock to make qtip
  
    m_info->hide();
  }//SimpleNuclide constructor
  
  
  void changeShielding( ShieldingLevel level )
  {
    m_shieldLevel = level;
    
    switch( level )
    {
      case NoShielding:
        m_noshield->setStyleClass( "SimpleNuclideShieldButActive" );
        m_someshield->setStyleClass( "SimpleNuclideShieldBut" );
        m_moreshield->setStyleClass( "SimpleNuclideShieldBut" );
        m_shielddesc->setText( "" );
        m_assist->setShielding( "", "0 cm" );
      break;
        
      case SomeShielding:
        m_noshield->setStyleClass( "SimpleNuclideShieldBut" );
        m_someshield->setStyleClass( "SimpleNuclideShieldButActive" );
        m_moreshield->setStyleClass( "SimpleNuclideShieldBut" );
        m_shielddesc->setText( "1 cm Fe" );
        m_assist->setShielding( "Fe", "1 cm" );
      break;
        
      case MoreShielding:
        m_noshield->setStyleClass( "SimpleNuclideShieldBut" );
        m_someshield->setStyleClass( "SimpleNuclideShieldBut" );
        m_moreshield->setStyleClass( "SimpleNuclideShieldButActive" );
        m_shielddesc->setText( "1 in Pb" );
        m_assist->setShielding( "Pb", "2.54 cm" );
      break;
      
      default:
        break;
    }//switch( level )
    
  }//void changeShielding( int level )
  
  bool setActive()
  {
    if( !m_info->isHidden() )
      return false;
    
    addStyleClass( "ActiveSimpleNuclide" );
    m_info->setHidden( false, s_show_animation );
    changeShielding( m_shieldLevel );
    
    return true;
  }//void setActive()
  
  
  void removeActive()
  {
    if( m_info->isHidden() )
      return;
    
    removeStyleClass( "ActiveSimpleNuclide" );
    m_info->setHidden( true, s_show_animation );
  }//void removeActive()
  
  
  ~SimpleNuclide()
  {
  }
};//class SimpleNuclide



SimpleNuclideAssistPopup::SimpleNuclideAssistPopup( const float energy,
                           SpectrumViewer *viewer, int pagex, int pagey )
  : PopupDivMenu( 0, PopupDivMenu::TransientMenu ),
    m_item( 0 ),
    m_assist( 0 )
{
  setPositionScheme( Wt::Absolute );
  addStyleClass( "SimpleNuclideAssistPopup" );
  aboutToHide().connect( this, &SimpleNuclideAssistPopup::handleAboutToHide );
  setAttributeValue( "style", "background: none; border: none;" );
  
  
  WText *title = new WText( "From Characteristic &gamma;&#8217;s:" );
  title->setInline( false );
  title->addStyleClass( "SimpleNuclideAssistPopupTitle" );
  PopupDivMenuItem *titeleItem = addWidget( title );
  titeleItem->setAttributeValue( "style", "background: none; margin: 0; padding: 0;" );
  
  doJavaScript( "$('#" + title->id() + "').data('menuparent'," + jsRef() + ");" );
  
  const char *mouseDownJs = "function(o,e){$(o).data('hasmouse',true).data('ds',Wt.WT.pageCoordinates(e));}";
  const char *mouseUpJs = "function(o,e){$(o).data('hasmouse',false);}";
  const char *mouseMoveJs = INLINE_JAVASCRIPT(
    function(o,e){
      var self = $(o);
      if( !self.data('hasmouse') )
        return;
    
      var nowxy = Wt.WT.pageCoordinates(e);
      var el = self.data('menuparent');
      var ds = self.data('ds');
    
      if( el.style.right === 'auto' || el.style.right === "" ){
        el.style.left = (Wt.WT.px(el, 'left') + nowxy.x - ds.x) + 'px';
        el.style.right = "";
      }else{
        el.style.right = (Wt.WT.px(el, 'right') + ds.x - nowxy.x) + 'px';
        el.style.left = 'auto';
      }
    
      if( el.style.bottom === 'auto' || el.style.bottom === "" ){
        el.style.top = (Wt.WT.px(el, 'top') + nowxy.y - ds.y) + 'px';
        el.style.bottom = "";
      }else{
        el.style.bottom = (Wt.WT.px(el, 'bottom') + ds.y - nowxy.y) + 'px';
        el.style.top = 'auto';
      }
      self.data('ds',nowxy);
    }
  );
  
  title->mouseWentDown().connect( mouseDownJs );
  title->mouseWentUp().connect( mouseUpJs );
  title->mouseMoved().connect( mouseMoveJs );
  
  
  try
  {
    m_assist = new SimpleNuclideAssist( energy, viewer );
    m_item = addWidget( m_assist );
    m_item->setAttributeValue( "style", "background: none; margin: 0; padding: 2;" );
    if( viewer->isMobile() )
      showFromMouseOver();
    else
      popup( WPoint(pagex,pagey) );
  }catch( std::exception &e )
  {
    cerr << "Failed to initialize SimpleNuclideAssist: " << e.what() << endl;
  }//try / catch
}//SimpleNuclideAssistPopup constructor


SimpleNuclideAssistPopup::~SimpleNuclideAssistPopup()
{
  
}


bool SimpleNuclideAssistPopup::isValid() const
{
  return m_item;
}


void SimpleNuclideAssistPopup::handleAboutToHide()
{
  if( m_item )
    delete m_item;
  m_item = 0;
  m_assist = 0;
}


SimpleNuclideAssist::SimpleNuclideAssist( const float energy,
                                          SpectrumViewer *viewer,
                                          WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_viewer( viewer ),
    m_sources( 0 ),
    m_messages( 0 ),
    m_legend( 0 ),
    m_energyRangeTxt( 0 ),
    m_contaminantLegTxt( 0 ),
    m_contaminantLegImg( 0 ),
    m_initialized( false ),
    m_lowerEnergy( 0.0f ),
    m_upperEnergy( 0.0f )
{
  addStyleClass( "SimpleNuclideAssist" );
  
  m_sources = new WContainerWidget( this );
  m_messages = new WContainerWidget( this );
  m_messages->addStyleClass( "SimpleNuclideMsg" );
  
  m_legend = new WContainerWidget( this );
  m_legend->addStyleClass( "SimpleNuclideAssistLegend" );
  WGridLayout *leglayout = new WGridLayout( m_legend );
  m_energyRangeTxt = new WText();
  leglayout->addWidget( m_energyRangeTxt, 0, 0, 1, 2 );
  
  WText *txt = new WText( "Characteristic regions" );
  leglayout->addWidget( txt, 1, 0, AlignLeft );

  WContainerWidget *indicator = new WContainerWidget();
  indicator->decorationStyle().setBackgroundColor( s_main_nuc_highlight_color );
  indicator->setStyleClass( "SimpleNuclideAssistLegendROI" );
  leglayout->addWidget( indicator, 1, 1 );

  
  m_contaminantLegTxt = new WText( "Contaminants regions" );
  leglayout->addWidget( m_contaminantLegTxt, 2, 0, AlignLeft );
  
  m_contaminantLegImg = new WContainerWidget();
  m_contaminantLegImg->decorationStyle().setBackgroundColor( s_contaminant_highlight_color );
  m_contaminantLegImg->setStyleClass( "SimpleNuclideAssistLegendROI" );
  leglayout->addWidget( m_contaminantLegImg, 2, 1 );
  
  m_legend->hide();
  
  setEnergy( energy );
  
  //It would be nice to include the message somewhere contectually near the
  //  popup, bu in the popup it takes up to much space.
  const char *msg = "For more advanced searching, or to search on"
                    " non-characteristic peaks, use the <b>Nuclide Search</b>"
                    " tool.";
  passMessage( msg, "", WarningWidget::WarningMsgInfo );
}//SimpleNuclideAssist


SimpleNuclideAssist::~SimpleNuclideAssist()
{
  try
  {
    PhotopeakLineDisplay *photopeaks = m_viewer->isotopeLinesWidget();
    if( photopeaks )
    {
      if( photopeaks && m_initialized && m_refPhotopeakIntialStateXml.size() > 10
        /*&& UtilityFunctions::contains(m_refPhotopeakIntialStateXml, "PhotopeakLineDisplay")*/ )
      {
        photopeaks->clearAllLines();
        photopeaks->deSerialize( m_refPhotopeakIntialStateXml  );
      }
    }
  }catch( std::exception &e )
  {
    
  }
  
  for( size_t i = 0; i < m_regionids.size(); ++i )
    m_viewer->removeHighlightedEnergyRange( m_regionids[i] );
  m_regionids.clear();
  
  cerr << "~SimpleNuclideAssist()" << endl;
}//~SimpleNuclideAssist()


void SimpleNuclideAssist::setEnergy( const float energy )
{
  float fwhm = 10.0f;
  
  std::shared_ptr<const SpecMeas> meas = m_viewer->measurment( kForeground );
  //const set<int> detector_numbers = m_viewer->displayedDetectorNumbers();
  const set<int> samples = m_viewer->displayedSamples( kForeground );

  if( meas )
  {
    std::shared_ptr<const DetectorPeakResponse> det = meas->detector();
    std::shared_ptr< const deque< std::shared_ptr<const PeakDef> > > peaks;
    
    peaks = meas->automatedSearchPeaks( samples );
    if( !peaks )
      peaks = meas->peaks( samples );
    
    if( det && det->hasResolutionInfo() )
    {
      fwhm = det->peakResolutionSigma( energy );
    }else if( peaks && peaks->size() )
    {
      vector< std::shared_ptr<const PeakDef> > peakv( peaks->begin(), peaks->end() );
      std::stable_sort( peakv.begin(), peakv.end(), &PeakDef::lessThanByMeanShrdPtr );
      size_t index = 0;
      while( index < peakv.size() && energy > peakv[index]->mean() )
        ++index;
      if( index == 0 )
        fwhm = static_cast<float>( peakv[0]->fwhm() );
      else if( index >= peakv.size() )
        fwhm = static_cast<float>( peakv.back()->fwhm() );
      else
        fwhm = static_cast<float>( std::max( peakv[index]->fwhm(), peakv[index-1]->fwhm() ) );
    }
  }//if( meas )
  
  fwhm = std::max( fwhm, 5.0f );
  const float mindispl = 7.0f*m_viewer->kevPerPixel();
  if( !IsNan(mindispl) && !IsInf(mindispl) )
    fwhm = std::max( fwhm, mindispl );
  
  setEnergyRange( energy - fwhm, energy + fwhm );
}//void setEnergy()


void SimpleNuclideAssist::setEnergyRange( const float lowerEnergy,
                                          const float upperEnergy )
{
  m_lowerEnergy = lowerEnergy;
  m_upperEnergy = upperEnergy;
  
  char buffer[64];
  snprintf( buffer, sizeof(buffer), "Nucs in %.1f to %.1f keV",
            lowerEnergy, upperEnergy );
  m_energyRangeTxt->setText( buffer );
  
  m_sources->clear();
  m_messages->clear();
  
  vector<SourceInfo> relevantinfo;
  
  foreach( const SourceInfo &source, parseXmlFile() )
  {
    foreach( const GammaInfo &gamma, source.gammas )
    {
      if( (gamma.energy >= lowerEnergy) && (gamma.energy <= upperEnergy) )
      {
        relevantinfo.push_back( source );
        break;
      }
    }//foreach( const GammaInfo &gamma, source.gammas )
  }//foreach( const SourceInfo &source, parseXmlFile() )
  
  foreach( const SourceInfo &source, relevantinfo )
  {
    SimpleNuclide *widget = new SimpleNuclide( source, this, m_viewer, m_sources );
    widget->clicked().connect( boost::bind( &SimpleNuclideAssist::updateSourceDisplayed, this, widget ) );
    widget->mouseWentOver().connect( boost::bind( &SimpleNuclideAssist::updateSourceDisplayed, this, widget ) );
  }//foreach( const SourceInfo &source, relevantinfo )
  
  if( m_initialized )
  {
    PhotopeakLineDisplay *photopeaks = m_viewer->isotopeLinesWidget();
    if( photopeaks )
      photopeaks->clearAllLines();
    
    for( size_t i = 0; i < m_regionids.size(); ++i )
      m_viewer->removeHighlightedEnergyRange( m_regionids[i] );
    m_regionids.clear();
  }//if( m_initialized )
  
  if( relevantinfo.empty() )
  {
    m_messages->show();
    char buff[128];
    snprintf( buff, sizeof(buff), "No nuclides between %.1f and %.1f keV",
              lowerEnergy, upperEnergy );
    new WText( buff, m_messages );
  }else
  {
    m_messages->hide();
  }
  
}//void setEnergyRange( const double lowerEnergy, const double upperEnergy )


void SimpleNuclideAssist::setShielding( const std::string &name,
                                        const std::string &thickness )
{
  PhotopeakLineDisplay *photopeaks = m_viewer->isotopeLinesWidget();
  if( photopeaks )
  {
    try
    {
      photopeaks->setShieldingMaterialAndThickness( name, thickness );
    }catch( std::exception &e )
    {
      const string msg = "There apears to be an error in "
                         + SimpleNuclideAssist::sm_dataFileName
                         + string(": ") + e.what()
                         + string("<br />The shielding has not been set.");
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
    }
  }//if( photopeaks )
}//void setShielding(...)


void SimpleNuclideAssist::updateSourceDisplayed( SimpleNuclide *widget )
{
  if( !widget->setActive() )
    return;
  
  foreach( WWidget *kidw, m_sources->children() )
  {
    SimpleNuclide *kid = static_cast<SimpleNuclide *>(kidw);
    if( kid != widget )
      kid->removeActive();
  }
  
  const SourceInfo &source = widget->m_source;
  PhotopeakLineDisplay *photopeaks = m_viewer->isotopeLinesWidget();
  if( photopeaks )
  {
    if( !m_initialized )
      photopeaks->serialize( m_refPhotopeakIntialStateXml  );
    
    photopeaks->clearAllLines();
    
    for( size_t i = 0; i < m_regionids.size(); ++i )
      m_viewer->removeHighlightedEnergyRange( m_regionids[i] );
    m_regionids.clear();
    
    //TODO - make width vary with energy based on the detectors resolution if
    //       it is around.
    const float width = 0.5f*(m_upperEnergy - m_lowerEnergy);
    
    //TODO - need to add in a mechanism to make currently showing nuclide
    //  different than normal (set color, dashed lines, etc), to indicate a
    //  contaminant.
    foreach( const string &contam, source.contaminants )
    {
      photopeaks->setNuclideAndAge( contam, "" );
      photopeaks->persistCurentLines();
    }
    
    
    
    photopeaks->setNuclideAndAge( source.name, "" );
    
    const vector<string> &contaminants = source.contaminants;
    
    foreach( const GammaInfo &gamma, source.gammas )
    {
      const bool is_contaminant = contains(contaminants, gamma.parent);
      
      const WColor &color = (is_contaminant ? s_contaminant_highlight_color
                                            : s_main_nuc_highlight_color);
      
      const float lowerE = gamma.energy - width;
      const float upperE = gamma.energy + width;
      const size_t regionid
               = m_viewer->addHighlightedEnergyRange( lowerE, upperE, color );
      m_regionids.push_back( regionid );
    }//foreach( const GammaInfo &gamma, source.gammas )
  }//if( photopeaks )
  
  const bool has_contam = !source.contaminants.empty();
  
  if( m_legend->isHidden() )
  {
    if( has_contam == m_contaminantLegTxt->isHidden() )
    {
      m_contaminantLegTxt->setHidden( !has_contam );
      m_contaminantLegImg->setHidden( !has_contam );
    }
    
    m_legend->setHidden( false, s_show_animation );
  }//if( m_legend->isHidden() )
  
  
  
  if( has_contam == m_contaminantLegTxt->isHidden() )
  {
    m_contaminantLegTxt->setHidden( !has_contam, s_appear_animation );
    m_contaminantLegImg->setHidden( !has_contam, s_appear_animation );
  }
  
  
  m_initialized = true;
}//void updateSourceDisplayed( SimpleNuclide *widget )




