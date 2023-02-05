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

#include <list>
#include <tuple>
#include <math.h>
#include <ctype.h>
#include <algorithm>

#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WTable>
#include <Wt/WString>
#include <Wt/WDialog>
#include <Wt/WLength>
#include <Wt/WPainter>
#include <Wt/WRectArea>
#include <Wt/WTableCell>
#include <Wt/WGridLayout>
#include <Wt/Chart/WAxis>
#include <Wt/WPushButton>
#include <Wt/WPaintDevice>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesianChart>

#include <boost/ref.hpp>
#include <boost/any.hpp>
#include <boost/bind.hpp>

#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayChainChart.h"
#include "InterSpec/DecayActivityDiv.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

namespace
{
  const char *productname( const SandiaDecay::ProductType &type )
  {
    switch( type )
    {
      case SandiaDecay::AlphaParticle:           return "alpha";            break;
      case SandiaDecay::CaptureElectronParticle: return "electron capture"; break;
      case SandiaDecay::BetaParticle:            return "beta";             break;
      case SandiaDecay::GammaParticle:           return "gamma";            break;
      case SandiaDecay::PositronParticle:        return "positron";         break;
      case SandiaDecay::XrayParticle:            return "xray";             break;
      default:                                   return "unknown";          break;
    }//switch( par.type )
    
    return "";
  }//productname(...)
  
  
  vector<std::tuple<SandiaDecay::ProductType,double,double>> decay_particle_info( const SandiaDecay::Nuclide * const nuc )
  {
    vector<std::tuple<SandiaDecay::ProductType,double,double>> answer;
    if( !nuc )
      return answer;
    
    for( const SandiaDecay::Transition * const trans : nuc->decaysToChildren )
    {
      for( const SandiaDecay::RadParticle &par : trans->products )
      {
        answer.emplace_back( par.type, par.energy, trans->branchRatio * par.intensity );
      }//for(SandiaDecay::RadParticle par : trans->products)
    }//for (const SandiaDecay::Transition * trans : nuc->decaysToChildren)
    
    std::sort( begin(answer), end(answer),
      [](const tuple<SandiaDecay::ProductType,double,double> &lhs, const tuple<SandiaDecay::ProductType,double,double> &rhs) -> bool {
        //Define the order of particle types we would like in the table
        const SandiaDecay::ProductType partorder[] = { SandiaDecay::XrayParticle,
          SandiaDecay::GammaParticle, SandiaDecay::AlphaParticle,
          SandiaDecay::BetaParticle, SandiaDecay::PositronParticle,
          SandiaDecay::CaptureElectronParticle
        };
        
        auto typeindex = [&]( const SandiaDecay::ProductType &part ) -> int {
          return static_cast<int>( std::find( begin(partorder), end(partorder), part) - begin(partorder) );
        };
        
        const int lhstype = typeindex(std::get<0>(lhs));
        const int rhstype = typeindex(std::get<0>(rhs));
        if( lhstype == rhstype )
          return std::get<1>(lhs) < std::get<1>(rhs);
        return lhstype < rhstype;
    } );
    
    return answer;
  }//decay_particle_info(...)
}//namespace


DecayChainChart::DecayChainChart( WContainerWidget *parent  )
  : WContainerWidget( parent ),
  m_useCurrie( true ),
  m_jsLoaded( false ),
  m_nuclide( nullptr ),
  m_moreInfoDialog( nullptr ),
  m_showDecayParticles( this, "ShowDecayParticleInfo", true ),
  m_showDecaysThrough( this, "ShowDecaysThrough", true )
{
  //To keep Wt from loading d3.js again, if it loaded it fo D3SpectrumDisplay, then
  //  it looks like (for Wt 3.3.4 anyway) we need the same URL, not just the same
  //  second argument - so the loading here matches D3SpectrumDisplay.
#if( defined(NDEBUG) || IOS || ANDROID || BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP )
  //THe NDEBUG should be enough, but just making sure
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
#else
  wApp->require( "external_libs/SpecUtils/d3_resources/d3.v3.min.js", "d3.v3.js" );
#endif
  
  //wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/DecayChainChart.js" );
  wApp->useStyleSheet( "InterSpec_resources/DecayChainChart.css" );
  addStyleClass( "DecayChainChart" );
  
  m_showDecayParticles.connect( boost::bind( &DecayChainChart::showDecayParticleInfo, this,
                                            boost::placeholders::_1 ) );
  m_showDecaysThrough.connect( boost::bind( &DecayChainChart::showDecaysThrough, this,
                                           boost::placeholders::_1 ) );
}//DecayChainChart constructor


pair<AuxWindow *, DecayChainChart *>
  DecayChainChart::show_decay_chart_window( const SandiaDecay::Nuclide *const nuc,
         const DecayChainChart::DecayChainType type )
{
  if( !nuc )
    return pair<AuxWindow *, DecayChainChart *>( nullptr, nullptr );

  
  string title;
  switch( type )
  {
    case DecayChainChart::DecayChainType::DecayFrom:
      title = nuc->symbol + " Decay Chain";
      break;

    case DecayChainChart::DecayChainType::DecayThrough:
      title = "Decays through " + nuc->symbol;
      break;
  }//switch( type )
  
  DecayChainChart *chart = new DecayChainChart();

  AuxWindow *window = new AuxWindow( title,
      (Wt::WFlags<AuxWindowProperties>( AuxWindowProperties::DisableCollapse ) 
         | AuxWindowProperties::EnableResize) );
  
  WGridLayout *layout = window->stretcher();
  layout->addWidget( chart, 0, 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setRowStretch( 0, 1 );

  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", viewer );

  viewer->colorThemeChanged();
  viewer->colorThemeChanged().connect( chart, &DecayChainChart::colorThemeChanged );

  chart->setNuclide( nuc, !useBq, type );

  if( viewer && (viewer->renderedWidth() > 100) && (viewer->renderedHeight() > 100) )
  {
    const int w = std::min( 800, viewer->renderedWidth() - 20 );
    const int h = std::min( 600, viewer->renderedHeight() );
    window->resizeWindow( w, h );

    window->resizeToFitOnScreen();
    window->centerWindowHeavyHanded();
  }//if( we know the window size )

  return make_pair( window, chart );
}//show_decay_chart_window


void DecayChainChart::doJavaScript( const std::string &js )
{
  if( m_jsLoaded )
    WContainerWidget::doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//doJavaScript(...)


void DecayChainChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
}//render(...)
  
  
void DecayChainChart::setNuclide( const SandiaDecay::Nuclide * const nuc, const bool useCurrie, const DecayChainType decayType )
{
  const bool actuallyChanged = (m_nuclide != nuc);
  
  m_nuclide = nuc;
  m_useCurrie = useCurrie;
  
  switch( decayType )
  {
    case DecayChainType::DecayFrom:
      setJsonForDecaysFrom( nuc );
      break;
      
    case DecayChainType::DecayThrough:
      setJsonForDecaysThrough( nuc );
      break;
  }//switch( decayType )
  
  if( actuallyChanged )
    m_nuclideChanged.emit( m_nuclide );
}//void setNuclide(...)


const SandiaDecay::Nuclide *DecayChainChart::nuclide() const
{
  return m_nuclide;
}//const SandiaDecay::Nuclide *nuclide() const;


void DecayChainChart::colorThemeChanged()
{
  InterSpec *interspec = InterSpec::instance();
  if( !interspec )
    return;
  
  shared_ptr<const ColorTheme> theme = interspec->getColorTheme();
  if( !theme )
    return;
  
  auto csstxt = []( const Wt::WColor &c, const char * defcolor ) -> string {
    return "'" + (c.isDefault() ? string(defcolor) : c.cssText()) + "'";
  };
  
  
  string options = "{";
  options += " lineColor: " + csstxt(theme->spectrumAxisLines,"black");
  options += ", textColor: " + csstxt(theme->spectrumChartText,"black");
  options += ", backgroundColor: " + csstxt(theme->spectrumChartBackground,"rgba(0,0,0,0)");
  options += ", linkColor: " + csstxt(theme->defaultPeakLine,"blue");
  options += "}";
  
  doJavaScript( jsRef() + ".chart.setColorOptions(" + options + ");" );
}//void colorThemeChanged();


Wt::Signal<const SandiaDecay::Nuclide *> &DecayChainChart::nuclideChanged()
{
  return m_nuclideChanged;
}


void DecayChainChart::jsonInfoForNuclide( const SandiaDecay::Nuclide * const nuc, Wt::WStringStream &js ) const
{
  assert( nuc );
    
  vector<string> info = getTextInfoForNuclide( nuc, m_nuclide, m_useCurrie );
    
  const string hl = (IsInf(nuc->halfLife) ? std::string("stable") : PhysicalUnits::printToBestTimeUnits(nuc->halfLife, 2));
    
  js << "{ \"nuclide\": \"" << nuc->symbol << "\","
    << " \"massNumber\": " << static_cast<int>(nuc->massNumber) << ","
    << " \"atomicNumber\": " << static_cast<int>(nuc->atomicNumber) << ","
    << " \"iso\": " << static_cast<int>(nuc->isomerNumber) << ","
    << " \"halfLive\": \"" << hl << "\","
    << " \"isPrimaryNuc\": " << std::string((nuc==m_nuclide) ? "true," : "false,")
    << " \"txtinfo\": [";
    
  for( size_t i = 0; i < info.size(); ++i )
    js << std::string(i ? ", \"" : "\"") << "" << Wt::Utils::htmlEncode( info[i] ) << "\"";
  js << "],";
    
  js << " \"children\": [";
  const auto decays = nuc->decaysToChildren;
  for( size_t i = 0, kidnum = 0; i < decays.size(); ++i )
  {
    const SandiaDecay::Transition *trans = decays[i];
    
    if( !trans->child )
      continue;  //Spontaneous Fission, skip for now.
    
    if( trans->branchRatio <= 0.0 )
      continue;
    
    js << std::string(kidnum ? ", {" : " {")
    << " \"nuclide\": \"" << trans->child->symbol << "\","
    << " \"massNumber\": " << static_cast<int>(trans->child->massNumber) << ","
    << " \"atomicNumber\": " << static_cast<int>(trans->child->atomicNumber) << ","
    << " \"iso\": " << static_cast<int>(trans->child->isomerNumber) << ","
    << " \"br\": " << static_cast<double>(trans->branchRatio) << ","
    << " \"mode\": \"" << std::string(SandiaDecay::to_str( trans->mode )) << "\""
    << "}";
    //trans->products
    
    kidnum++;
  }
  js << "]}";
}//void jsonInfoForNuclide(...)
  
  
void DecayChainChart::setJsonForDecaysFrom( const SandiaDecay::Nuclide * const nuclide )
{
  if( !nuclide )
  {
    doJavaScript( jsRef() + ".chart.setDecayData('[]',true);" );
    return;
  }
  
  vector<const SandiaDecay::Nuclide *> descendants = nuclide->descendants();
  
  WStringStream js;
  js << "[";
  for( size_t i = 0; i < descendants.size(); ++i )
  {
    js << std::string(i ? "," : "");
    jsonInfoForNuclide( descendants[i], js );
  }
  js << "]";
  
  //cout << "JSON=" << js.str() << endl;
  
  doJavaScript( jsRef() + ".chart.setDecayData('" + js.str() + "',true);" );
}//void setJsonForDecaysFrom( const SandiaDecay::Nuclide * const nuclide )


void DecayChainChart::setJsonForDecaysThrough( const SandiaDecay::Nuclide * const nuclide )
{
  if( !nuclide )
  {
    doJavaScript( jsRef() + ".chart.setDecayData('[]',false);" );
    return;
  }
  
  vector<const SandiaDecay::Nuclide *> forebearers = nuclide->forebearers();
  
  //Get rid of elements above Californium since they arent typically encountered
  forebearers.erase( std::remove_if( begin(forebearers), end(forebearers),
                     [](const SandiaDecay::Nuclide *nuc) -> bool { return nuc->atomicNumber>98;} ),
                     end(forebearers) );
  
  WStringStream js;
  js << "[";
  for( size_t i = 0; i < forebearers.size(); ++i )
  {
    js << std::string(i ? "," : "");
    jsonInfoForNuclide( forebearers[i], js );
  }
  js << "]";
  
  //cout << "JSON=" << js.str() << endl;
  
  doJavaScript( jsRef() + ".chart.setDecayData('" + js.str() + "',false);" );
}//void setJsonForDecaysThrough( const SandiaDecay::Nuclide * const nuc )
  
  
std::vector<std::string> DecayChainChart::getTextInfoForNuclide( const SandiaDecay::Nuclide * const nuc,
                                                                  const SandiaDecay::Nuclide * const parentNuclide,
                                                                  const bool useCurrie )
{
  vector<string> information;
  
  if( !nuc )
    return information;
  
  char buffer[512];
  
  information.push_back( nuc->symbol );
  snprintf( buffer, sizeof(buffer), "Atomic Number: %i", nuc->atomicNumber );
  information.push_back( buffer );
  snprintf( buffer, sizeof(buffer), "Atomic Mass: %.2f", nuc->atomicMass );
  information.push_back( buffer );
  
  if( nuc == parentNuclide )
  {
    if( parentNuclide->canObtainSecularEquilibrium() )
      information.push_back( "Can reach secular equilibrium" );
    else
      information.push_back( "Cannot reach secular equilibrium" );
  }
  
  if( IsInf(nuc->halfLife) )
  {
    information.push_back("Half Life: stable" );
  }else
  {
    const string hl = PhysicalUnits::printToBestTimeUnits( nuc->halfLife );
    information.push_back("Half Life: " + hl );
  }
  
  if( parentNuclide && (nuc != parentNuclide) )
  {
    const float br_from = parentNuclide->branchRatioFromForebear(nuc);  //will be >0 when we are plotting decays through a nuclide
    const float br_to = parentNuclide->branchRatioToDecendant(nuc);     //will be >0 when plotting decays from a nuclide
    
    if( br_from <= 0.0 || br_to > 0.0 )
      snprintf( buffer, sizeof(buffer), "Branch Ratio from %s: %.5g", parentNuclide->symbol.c_str(), br_to );
    else
      snprintf( buffer, sizeof(buffer), "Branch Ratio through %s: %.5g", parentNuclide->symbol.c_str(), br_from );
      
    information.push_back( buffer );
  }//if( nuc == parentNuclide )
  
  if( !IsInf(nuc->halfLife) )
  {
    const double specificActivity = nuc->activityPerGram() / PhysicalUnits::gram;
    const string sa = PhysicalUnits::printToBestSpecificActivityUnits( specificActivity, 3, useCurrie );
    information.push_back("Specific Act: " + sa );
  }//if( not stable )
  
  
  for( const SandiaDecay::Transition * transition : nuc->decaysToChildren )
  {
    if( transition->child && (transition->branchRatio > 0.0) )
    {
      const char * const decay_mode = SandiaDecay::to_str( transition->mode );
      const string br = PhysicalUnits::printCompact(transition->branchRatio, 4);
      
      snprintf( buffer, sizeof(buffer),
               "Decays to %s by %s decay, BR %s",
               transition->child->symbol.c_str(),
               decay_mode, br.c_str() );
      
      information.push_back( buffer );
    }//if( transition->child )
  }//for( const SandiaDecay::Transition * transition : nuc->decaysToChildren)
  
  return information;
}//std::vector<std::string> getTextInfoForNuclide( const SandiaDecay::Nuclide *nuc )


void DecayChainChart::defineJavaScript()
{
  m_jsLoaded = true;
  
  InterSpec *interspec = InterSpec::instance();
  const bool isMobile = interspec && interspec->isMobile();
  
  shared_ptr<const ColorTheme> theme;
  if( interspec )  //should always be valid, but jic
    theme = interspec->getColorTheme();
  
  string options = string("{")
  + " isMobile: " + string( isMobile ? "true" : "false");
  
  if( theme )
  {
    auto csstxt = []( const Wt::WColor &c, const char * defcolor ) -> string {
      return "'" + (c.isDefault() ? string(defcolor) : c.cssText()) + "'";
    };
    
    options += ", lineColor: " + csstxt(theme->spectrumAxisLines,"black");
    options += ", textColor: " + csstxt(theme->spectrumChartText,"black");
    options += ", backgroundColor: " + csstxt(theme->spectrumChartBackground,"rgba(0,0,0,0)");
    options += ", linkColor: " + csstxt(theme->defaultPeakLine,"blue");
  }//if( theme )
  
  options += "}";
  
  
  setJavaScriptMember( "chart", "new DecayChainChart(" + jsRef() + "," + options + ");");
  setJavaScriptMember( WT_RESIZE_JS, "function(self, w, h, layout){ " + jsRef() + ".chart.handleResize();}" );
  
  for( const string &js : m_pendingJs )
    WContainerWidget::doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void DecayChainChart::showDecaysThrough( const std::string nuc )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nucptr = db->nuclide( nuc );
  
  //cout << "DecayChainChart::showDecaysThrough: '" << nuc << "'" << endl;
  
  showPossibleParents( nucptr );
}//void showDecaysThrough( const std::string nuc )


void DecayChainChart::showDecayParticleInfo( const std::string &csvIsotopeNames )
{
  vector<string> isotopenames;
  SpecUtils::split( isotopenames, csvIsotopeNames, "," );
  
  map<string,vector<tuple<SandiaDecay::ProductType,double,double>>> nucinfos;
  
  for( const string &isoname : isotopenames )
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *nuc = db->nuclide( isoname );
    if( nuc )
      nucinfos[isoname] = decay_particle_info( nuc );
  }//for( const string &isoname : isotopenames )
  
  deleteMoreInfoDialog();
  
  m_moreInfoDialog = new AuxWindow( "Decay Particle Energies" );
  m_moreInfoDialog->contents()->setMinimumSize(250, 80);
  m_moreInfoDialog->contents()->setMaximumSize(330, 400);
  m_moreInfoDialog->contents()->setOverflow(Wt::WContainerWidget::OverflowAuto, Wt::Vertical);
  m_moreInfoDialog->setClosable( true );
  m_moreInfoDialog->setModal( false );
  m_moreInfoDialog->rejectWhenEscapePressed();
  m_moreInfoDialog->disableCollapse();
  m_moreInfoDialog->centerWindow();
  m_moreInfoDialog->finished().connect( this, &DecayChainChart::deleteMoreInfoDialog );
  WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter("Close");
  ok->clicked().connect( this, &DecayChainChart::deleteMoreInfoDialog );
  
  int nucs = 0;
  for( const auto &nuc_infos : nucinfos )
  {
    auto header = new Wt::WText( "Particles from " + nuc_infos.first, m_moreInfoDialog->contents() );
    header->addStyleClass( "DecayPartInfoHeader" );
    header->setInline( false );
    
    if( nucs++ )
      header->setMargin( 20, Wt::Top );
    
    if( nuc_infos.second.empty() )
    {
      const char *msg = "<center>This nuclide has no particles from decay</center>";
      auto blank = new Wt::WText( msg, Wt::XHTMLText, m_moreInfoDialog->contents() );
      blank->setInline( false );
      continue;
    }
    
    WTable *table = new Wt::WTable( m_moreInfoDialog->contents() );
    table->addStyleClass( "DecayChainChartTable" );
    table->setHeaderCount( 1 );
    table->elementAt(0, 0)->addWidget( new Wt::WText("Particle") );
    table->elementAt(0, 1)->addWidget( new Wt::WText("Energy (keV)") );
    table->elementAt(0, 2)->addWidget( new Wt::WText("Intensity") );
    
    int row = 0;
    for( const auto &infos : nuc_infos.second )
    {
      ++row;
      
      char energyStr[32], intensityStr[32];
      snprintf( energyStr, sizeof(energyStr), "%.2f", std::get<1>(infos) );
      snprintf( intensityStr, sizeof(intensityStr), "%.5g", std::get<2>(infos) );
      
      table->elementAt(row, 0)->addWidget( new Wt::WText( productname(std::get<0>(infos)) ) );
      table->elementAt(row, 1)->addWidget( new Wt::WText( energyStr ) );
      table->elementAt(row, 2)->addWidget( new Wt::WText( intensityStr ) );
    }//for( loop over infos )
  }//for( loop over nucinfos )
  

  m_moreInfoDialog->show();
}//void showDecayParticleInfo( const std::string &csvIsotopeNames )


void DecayChainChart::showPossibleParents( const SandiaDecay::Nuclide *nuclide )
{
  if( !nuclide )
    return;
  
  InterSpec *interspec = InterSpec::instance();
  if( !interspec )  //shouldnt ever happen, but jic
    return;
  
  deleteMoreInfoDialog();
  
  const double ww = 0.8*interspec->renderedWidth();
  const double wh = 0.65*interspec->renderedHeight();
  
  Wt::WFlags<AuxWindowProperties> windowProp
  = Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
  | AuxWindowProperties::DisableCollapse
  | AuxWindowProperties::SetCloseable
  | AuxWindowProperties::EnableResize;
  
  m_moreInfoDialog = new AuxWindow( "Decays Through " + nuclide->symbol, windowProp );
  DecayChainChart *w = new DecayChainChart();
  w->setNuclide( nuclide, m_useCurrie, DecayChainChart::DecayChainType::DecayThrough );
  m_moreInfoDialog->stretcher()->addWidget( w, 0, 0  );
  
  m_moreInfoDialog->resizeWindow( std::max( 800, static_cast<int>(std::min(ww,1.3*wh) ) ), std::max( 600, static_cast<int>(wh) ) );
  
  
  m_moreInfoDialog->rejectWhenEscapePressed();
  m_moreInfoDialog->finished().connect( this, &DecayChainChart::deleteMoreInfoDialog );
  WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter("Close");
  ok->clicked().connect( this, &DecayChainChart::deleteMoreInfoDialog );
  
  
  m_moreInfoDialog->show();
  m_moreInfoDialog->centerWindowHeavyHanded();
}//void showPossibleParents(...)


void DecayChainChart::deleteMoreInfoDialog()
{
  if( m_moreInfoDialog )
    delete m_moreInfoDialog;
  m_moreInfoDialog = nullptr;
}//void deleteMoreInfoDialog()


