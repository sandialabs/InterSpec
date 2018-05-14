
#include <Wt/WConfig.h>
#include <Wt/WText>
#include <Wt/WLink>
#include <Wt/WImage>
#include <Wt/WCheckBox>
#ifdef WT_HAS_WRASTERIMAGE
#include <Wt/WRasterImage>
#else
#include <Wt/WSvgImage>
#endif
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WContainerWidget>

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "InterSpec/SrbActivityChart.h"
#include "InterSpec/SrbActivityApp.h"
#include "InterSpec/SrbActivityDiv.h"
//#include "InterSpec/SrbHeaderFooter.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ChartToImageResource.h"

using namespace Wt;
using namespace std;


#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


//use for served things like css sheet
std::string getSrbComonResourceUrl()
{
  const std::string &hostname = wApp->environment().hostName();
  const bool isOutrage = boost::algorithm::icontains( hostname, "outrage" )
                         || boost::algorithm::icontains( hostname, "hekili" );

  if( isOutrage )
  {
    const string url = wApp->url();

    //get the leading directory deployment
    string userDir = "";
    vector< string > splitVec;
    boost::algorithm::split( splitVec, url, boost::is_any_of("/"),
                            boost::algorithm::token_compress_on );

    if( splitVec.size() ) userDir = splitVec[0];
    //if url started with a '/', then splitVec[0]==""
    if( userDir.empty() && (splitVec.size() > 1) ) userDir = splitVec[1];

    return wApp->makeAbsoluteUrl( "/" + userDir + "/common" );
  }//if( on outrage )

  return "srb_common_for_not_on_outrage";
}//std::string getSrbComonResourceUrl()




//use for executing things in shells
std::string getSrbComonResourcePath()
{
  const std::string &hostname = wApp->environment().hostName();
  const bool isOutrage = boost::algorithm::icontains( hostname, "outrage" )
                        || boost::algorithm::icontains( hostname, "hekili" );


  if( isOutrage )
  {
    const string url = wApp->url();

    if( boost::starts_with(url, "/srb/" )     )
      return "/home/srb-prod/public_html/common";
    if( boost::starts_with(url, "/~srb/" )    )
      return "/home/srb/public_html/common";
    if( boost::starts_with(url, "/srb-dev/" )
       || boost::starts_with(url, "/~srb-dev/" ) )
      return "/home/srb-dev/public_html/common";
    if( boost::starts_with(url, "/srb-qc/" )
       || boost::starts_with(url, "/~srb-qc/" ) )
      return "/home/srb-qc/public_html/common";
    if( boost::starts_with(url, "/srb-prod/" )
       || boost::starts_with(url, "/~srb-prod/" ) )
      return "/home/srb-prod/public_html/common";

    string errmsg = "Unknown App url='" + wApp->url() + "'";
    cerr << errmsg << endl; //so it makes it into the apache log files
  }//if( on outrage )

  return "srb_common_for_not_on_outrage";
}//std::string getSrbComonResourcePath()



SrbActivityApp::SrbActivityApp( const Wt::WEnvironment& env , Wt::WContainerWidget *contentDiv)
  : WApplication( env ), m_activityDiv( NULL )
{
//  const string srb_common_url = getSrbComonResourceUrl();
//  useStyleSheet( srb_common_url + "/assets/css/styles.css" );
//  useStyleSheet( "local_resources/SrbActivityFuture.css" );
//  
//  useStyleSheet("resources/themes/default/wt.css");
//
//  if( env.agentIsIE() )
//    useStyleSheet( "local_resources/IESrbActivityFuture.css" );
//
//  wApp->doJavaScript( "if(typeof console==='undefined'){console={log:function(){}};}" );
//  
//  styleSheet().addRule( "input[type=\"text\"]", "font-size:0.95em;" );
//  styleSheet().addRule( "button[type=\"button\"]", "font-size:0.9em;" );
//  
//  if( isMobile() )
//  {
//#if( WT_VERSION < 0x3030000 )
//  doJavaScript(
//               "var t=document.createElement('meta');"
//               "t.name = \"viewport\";"
//               "t.content = \"user-scalable=no,width=device-width\";"
//               "document.getElementsByTagName('head')[0].appendChild(t);"
//               "document.ontouchmove = function(e){e.preventDefault();};"
//               "document.ongesturechange = function(e){e.preventDefault();};"
//               );
//#else
//  doJavaScript(
//               "var t=document.createElement('meta');"
//               "t.name = \"viewport\";"
//               "t.content = \"user-scalable=no,initial-scale=1.0\";"//,width=device-width,height=device-height\";"
//               "document.getElementsByTagName('head')[0].appendChild(t);" );
//#endif
//  }//if( isMobile() )
//
  init( env, contentDiv );
}//SrbActivityApp( constructor )


void SrbActivityApp::init( const Wt::WEnvironment &env, WContainerWidget *contentDiv )
{
  contentDiv->clear();
  contentDiv->setInline( false );

//  m_activityDiv = new SrbActivityDiv(contentDiv);

  //We will go through and see if the user wants to preload an isotope, and
  //  if the header and footer should be shown
//  typedef std::vector<std::string> ParameterValues;
//  typedef std::map<std::string, ParameterValues> ParameterMap;
  const Http::ParameterMap &parmap = env.getParameterMap();

  vector<string> isotopes;
  vector<string> activitystrs, agestrs;
  bool show_header = !isPhone();
  bool show_footer = show_header;
  bool show_photopeaks_first = false;
  bool image_only = false, logy = false;
  int pngwidth = 800, pngheight = 600;  
  string photopeaktimestr;
  double photopeaktime = -1.0;


  foreach( const Http::ParameterMap::value_type &p, parmap )
  {
    string parname = p.first;
    vector<std::string> parvalues = p.second;
    if( parvalues.empty() )
      continue;

//    boost::algorithm::to_lower( parname );
//    foreach( string &str, parvalues )
//      boost::algorithm::to_lower( str );

    if( parname=="header" )
    {
      if( parvalues.front()=="false" || parvalues.front()=="0" )
        show_header = false;
      else if( parvalues.front()=="true" || parvalues.front()=="1" )
        show_header = true;
    }else if( parname=="footer" )
    {
      if( parvalues.front()=="false" || parvalues.front()=="0" )
        show_footer = false;
      else if( parvalues.front()=="true" || parvalues.front()=="1" )
        show_footer = true;
    }else if( parname=="photopeaks" || parname=="photolines"
              || parname=="peaks" || parname=="lines" )
    {
      if( parvalues.front()=="false" || parvalues.front()=="0" )
        show_photopeaks_first = false;
      else if( parvalues.front()=="true" || parvalues.front()=="1" )
        show_photopeaks_first = true;
    }else if( parname=="iso" || parname=="nuc"
              || parname=="isotope" || parname=="nuclide" )
    {
      isotopes = parvalues;
    }else if( parname=="age" || parname=="initialage" )
    {
      agestrs = parvalues;
    }else if( parname=="activity" || parname=="act" )
    {
      activitystrs = parvalues;
    }else if( parname=="image" || parname=="imageonly"
              || parname=="pic" || parname=="piconly" )
    {
      if( parvalues.front()=="false" || parvalues.front()=="0" )
        image_only = false;
      else if( parvalues.front()=="true" || parvalues.front()=="1" )
        image_only = true;
    }else if( parname=="log" || parname=="logy" )
    {
      if( parvalues.front()=="false" || parvalues.front()=="0" )
        logy = false;
      else if( parvalues.front()=="true" || parvalues.front()=="1" )
        logy = true;
    }else if( parname=="width" || parname=="w" )
    {
      try { pngwidth = boost::lexical_cast<int>( parvalues.front() ); }catch(...){}
    }else if( parname=="height" || parname=="h" )
    {
      try { pngheight = boost::lexical_cast<int>( parvalues.front() ); }catch(...){}
    }else if( parname=="time" || parname=="t" )
    {
      photopeaktimestr = parvalues.front();
    }
  }//foreach( const Http::ParameterMap::value_type &p, parmap)

  size_t nNonEmtpyIso = 0;
  bool add_default_iso = true;

  for( size_t isoindex = 0; isoindex < isotopes.size(); ++isoindex )
  {
    const string isoname = isotopes[isoindex];
    string isoagestr, isoactivitystr;
    if( isoindex < activitystrs.size() )
      isoactivitystr = activitystrs[isoindex];
    if( isoindex < agestrs.size() )
      isoagestr = agestrs[isoindex];

//    std::cerr << "For isoname='" << isoname << "', isoactivitystr='"
//              << isoactivitystr << "', isoagestr='" << isoagestr << "'"
//              << std::endl;

    if( !isoname.empty() )
     nNonEmtpyIso++;

    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *nuc = db->nuclide( isoname );

    if( nuc && !isinf(nuc->halfLife) )
    {
      add_default_iso = false;
      double act = 1.0*PhysicalUnits::mCi, age = 0.0;
      try
      {
        act = PhysicalUnits::stringToActivity( isoactivitystr );
      }catch(...){}

      try
      {
        age = PhysicalUnits::stringToTimeDuration( isoagestr );
      }catch(...)
      {
        //lets also add the option to specify age in half-lives
        if( boost::algorithm::icontains(isoagestr,"hl")
            || boost::algorithm::icontains(isoagestr,"half") )
        {
          stringstream(isoagestr) >> age;
          age *= nuc->halfLife;
        }
      }//try / catch for age

      if( !photopeaktimestr.empty() )
      {
        try
        {
          photopeaktime = PhysicalUnits::stringToTimeDuration( photopeaktimestr );
        }catch(...)
        {
          //lets also add the option to specify age in half-lives
          if( boost::algorithm::icontains(photopeaktimestr,"hl")
              || boost::algorithm::icontains(photopeaktimestr,"half") )
          {
            if( (stringstream(photopeaktimestr) >> photopeaktime) )
              photopeaktime *= nuc->halfLife;
            else
              photopeaktime = -1.0;
          }
        }//try / catch for age
      }//if( !photopeaktimestr.empty() )

      age = max( 0.0, age );
      if( act <= 0.0 )
        act = 1.0*PhysicalUnits::mCi;

      const bool useCurries = !boost::algorithm::icontains( isoactivitystr, "b" );
      const PhysicalUnits::UnitNameValuePair &bestActUnit
                           = PhysicalUnits::bestActivityUnit( act, useCurries );
      m_activityDiv->addNuclide( nuc->atomicNumber, nuc->massNumber,
                               nuc->isomerNumber, act, bestActUnit.first, age );
      }//if( nuc )
  }//for( size_t isoindex = 0; isoindex < isotopes.size(); ++isoindex )

  //Must set photopeak time slider after adding all the nuclides, or else
  //  the photopeak time slider will just be changed anyway
  if( photopeaktime >= 0.0 )
    m_activityDiv->setPhopeakSliderTime( photopeaktime );

  if( add_default_iso && !nNonEmtpyIso )
  {
    m_activityDiv->addNuclide( 53, 135, false, 1.0*PhysicalUnits::MBq, "MBq", 0.0 );
    m_activityDiv->displayInitialLoadPopup( false );
  }//if(  add_default_iso )


  if( logy )
  {
    m_activityDiv->m_logYScale->setChecked();
    m_activityDiv->updateYScale();
    m_activityDiv->setPhotoPeakChartLogY( logy );
  }//if( logy )

  if( show_photopeaks_first )
    m_activityDiv->showPhotopeakTab();

  if( image_only )
  {
    m_activityDiv->m_photoPeakChart->setBackground( WBrush(transparent) );
    m_activityDiv->m_decayChart->setBackground( WBrush(transparent) );


    pngwidth = max( pngwidth, 100 );
    pngheight = max( pngheight, 100 );
assert(0);
/*
#ifdef WT_HAS_WRASTERIMAGE
    WRasterImage *pngImage = new WRasterImage( "png", pngwidth,
                                               pngheight, this );
#else
    WSvgImage *pngImage = new WSvgImage( pngwidth, pngheight, this );
#endif
    {//begin paint raster image code block
      Wt::WPainter p( pngImage );
      if( show_photopeaks_first )
        m_activityDiv->m_photoPeakChart->paint( p );
      else
        m_activityDiv->m_decayChart->paint( p );
    }//end paint raster image code block

    WImage *image = new WImage( WLink(pngImage), root() );
    image->setId( "image" );
    image->imageLoaded().connect( boost::bind( &WApplication::quit, this ) );

    delete m_activityDiv;
    m_activityDiv = NULL;
*/
  }else
  {
    WGridLayout *layout = new WGridLayout();
  //  layout->setSpacing( 0 );
    layout->setContentsMargins( 0, 0, 0, 0 );
    int row = 0;
//  #if(!DONT_SHOW_SRB_HEADER_FOOTER)
//    if( show_header )
//    {
//      string urlLocal = "hekili.ca.sandia.gov";  //url(); ?? wApp->url(); ??
//      SrbHeader *header = new SrbHeader( "", urlLocal );
//      layout->addWidget( header, row++, 0, 1, 1 );
//    }//if( show_header )
//  #endif
    layout->addWidget( m_activityDiv, row++, 0, 1, 1 );
    layout->setRowStretch( row - 1, 5 );
//  #if(!DONT_SHOW_SRB_HEADER_FOOTER)
//    if( show_footer )
//    {
//      SrbFooter *footer = new SrbFooter();
//      layout->addWidget( footer, row++, 0, 1, 1 );
//    }//if( show_footer )
//  #endif

contentDiv->setLayout( layout );
contentDiv->setStyleClass( "root" );
  }//if( image_only )
}//void init( const Wt::WEnvironment& env )


SrbActivityApp::~SrbActivityApp()
{
}//~SrbActivityApp()



bool SrbActivityApp::isMobile() const
{
  const WEnvironment &env = environment();
  const bool isMob = (   env.agentIsMobileWebKit()
                      || env.agentIsIEMobile()
                      || env.agentIsMobileWebKit()
                      || env.agent()==WEnvironment::MobileWebKitiPhone
                      || env.agent()==WEnvironment::MobileWebKitAndroid
                      || env.agent()==WEnvironment::MobileWebKit
                      || env.userAgent().find("Opera Mobi") != std::string::npos
                      || env.userAgent().find("Android") != std::string::npos
                      || env.userAgent().find("RIM ") != std::string::npos
                      || env.userAgent().find("iPad") != std::string::npos
                      || env.userAgent().find("iPhone") != std::string::npos
                      || env.userAgent().find("iPod") != std::string::npos
                      );
  return isMob;
}//bool isMobile() const


bool SrbActivityApp::isPhone() const
{
  //see: (iOS) http://www.enterpriseios.com/wiki/UserAgent
  //     (Android) http://www.gtrifonov.com/2011/04/15/google-android-user-agent-strings-2/
  
  const WEnvironment &env = environment();
  return ( env.userAgent().find("iPhone") != std::string::npos
          || env.userAgent().find("iPod") != std::string::npos
          || (env.userAgent().find("Android") != std::string::npos
              && env.userAgent().find("Mobile") != std::string::npos)
          || env.userAgent().find("RIM ") != std::string::npos);
}//bool isPhone()


