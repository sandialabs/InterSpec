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
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesianChart>

#include <boost/ref.hpp>
#include <boost/any.hpp>
#include <boost/bind.hpp>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
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
  
#if( !USE_D3_DECAY_CHAIN_IMP )
  //OnClickDownJs is not a full JS function; addiotinally the 'id' and 'iso'
  //  variables must be defined, an the function should accept a 's' and 'e'
  //  argument similar to OnClickUpJs
  const char *OnClickDownJs = INLINE_JAVASCRIPT
  (
   try{
     var a = $(s);
     var fcn = function()
     {
       try
       {
         a.data('touchdown',null);
         Wt.emit(id, {name: 'tapnhold', eventObject: s}, iso);
       }catch(a){}
     };
     var timeout = setTimeout(fcn,600);
     a.data('touchdown', timeout);
   }catch(a){}
   );
  
  
  const char *OnClickUpJs = INLINE_JAVASCRIPT
  (
   function(s,e){
     try{
       var a = $(s);
       var to = a.data('touchdown');
       if( to )
       {
         window.clearTimeout(to);
         a.data('touchdown',null);
       }
     }catch(a){}
   }
   );
  
  const char *particlename( const SandiaDecay::RadParticle &par )
  {
    return productname( par.type );
  }//const char *particlename( SandiaDecay::RadParticle par )
  
  
  void fill_information( const SandiaDecay::Nuclide *nuc,
                        std::vector<std::string> &types,
                        std::vector<double> &energies,
                        std::vector<double> &intensities,
                        std::vector<std::string> &names )
  {
    if( !nuc )
      return;
    
    const Wt::WString text = nuc->symbol;
    
    //the non metastable child (if any)
    const SandiaDecay::Nuclide * metaChild = nullptr;
    
    for(const SandiaDecay::Transition * transition : nuc->decaysToChildren)
    {
      if( transition->child
         && nuc->massNumber == transition->child->massNumber
         && nuc->atomicNumber == transition->child->atomicNumber )
      {
        metaChild = transition->child;
      }
    }//for( SandiaDecay::Transition * transition : nuc->decaysToChildren )
    
    names.push_back(nuc->symbol);
    for (const SandiaDecay::Transition * trans : nuc->decaysToChildren)
    {
      for(SandiaDecay::RadParticle par : trans->products)
      {
        string name = particlename( par );
        types.push_back( particlename( par ) );
        energies.push_back( par.energy );
        intensities.push_back( trans->branchRatio * par.intensity );
      }//for(SandiaDecay::RadParticle par : trans->products)
    }//for (const SandiaDecay::Transition * trans : nuc->decaysToChildren)
    
    //check if we need two tables
    if(  metaChild )
    {
      types.push_back("");
      energies.push_back(0);
      intensities.push_back(0);
      names.push_back(metaChild->symbol);
      for(const SandiaDecay::Transition * trans : metaChild->decaysToChildren)
      {
        for(SandiaDecay::RadParticle par : trans->products)
        {
          types.push_back( particlename( par ) );
          energies.push_back( par.energy );
          intensities.push_back( trans->branchRatio * par.intensity );
        }//for(SandiaDecay::RadParticle par : trans->products)
      }//for (const SandiaDecay::Transition * trans : metaChild->decaysToChildren)
    }//if( metaChild )
  }//void fill_information(...)
#endif //!USE_D3_DECAY_CHAIN_IMP
}//namespace


#if( USE_D3_DECAY_CHAIN_IMP )
DecayChainChartD3::DecayChainChartD3( WContainerWidget *parent  )
  : WContainerWidget( parent ),
  m_useCurrie( true ),
  m_jsLoaded( false ),
  m_nuclide( nullptr ),
  m_moreInfoDialog( nullptr ),
  m_showDecayParticles( this, "ShowDecayParticleInfo", true ),
  m_showDecaysThrough( this, "ShowDecaysThrough", true )
{
  wApp->require( "InterSpec_resources/d3.v3.min.js" );
  wApp->require( "InterSpec_resources/DecayChainChart.js" );
  wApp->useStyleSheet( "InterSpec_resources/DecayChainChart.css" );
  addStyleClass( "DecayChainChart" );
  
  m_showDecayParticles.connect( boost::bind( &DecayChainChartD3::showDecayParticleInfo, this, _1 ) );
  m_showDecaysThrough.connect( boost::bind( &DecayChainChartD3::showDecaysThrough, this, _1 ) );
}//DecayChainChartD3 constructor


void DecayChainChartD3::doJavaScript( const std::string &js )
{
  if( m_jsLoaded )
    WContainerWidget::doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//doJavaScript(...)


void DecayChainChartD3::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
}//render(...)
  
  
void DecayChainChartD3::setNuclide( const SandiaDecay::Nuclide * const nuc, const bool useCurrie, const DecayChainType decayType )
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


const SandiaDecay::Nuclide *DecayChainChartD3::nuclide() const
{
  return m_nuclide;
}//const SandiaDecay::Nuclide *nuclide() const;


void DecayChainChartD3::colorThemeChanged()
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  if( !app )
    return;
  
  InterSpec *interspec = app->viewer();
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


Wt::Signal<const SandiaDecay::Nuclide *> &DecayChainChartD3::nuclideChanged()
{
  return m_nuclideChanged;
}


void DecayChainChartD3::jsonInfoForNuclide( const SandiaDecay::Nuclide * const nuc, Wt::WStringStream &js ) const
{
  assert( nuc );
    
  vector<string> info = getTextInfoForNuclide( nuc, m_nuclide, false, m_useCurrie );
    
  const string hl = (IsInf(nuc->halfLife) ? std::string("stable") : PhysicalUnits::printToBestTimeUnits(nuc->halfLife));
    
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
  
  
void DecayChainChartD3::setJsonForDecaysFrom( const SandiaDecay::Nuclide * const nuclide )
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


void DecayChainChartD3::setJsonForDecaysThrough( const SandiaDecay::Nuclide * const nuclide )
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
  
  
std::vector<std::string> DecayChainChartD3::getTextInfoForNuclide( const SandiaDecay::Nuclide *nuc,
                                                                  const SandiaDecay::Nuclide *parentNuclide,
                                                                  const bool includeDaughterIsomerics,
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
      information.push_back("Can reach secular equilibrium");
    else
      information.push_back("Cannot reach secular equilibrium");
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
    snprintf( buffer, sizeof(buffer), "Branch Ratio from %s: %.5g",
             parentNuclide->symbol.c_str(),
             parentNuclide->branchRatioToDecendant(nuc) );
    information.push_back( buffer );
  }//if( nuc == parentNuclide )
  
  if( !IsInf(nuc->halfLife) )
  {
    const double specificActivity = nuc->activityPerGram() / PhysicalUnits::gram;
    const string sa = PhysicalUnits::printToBestSpecificActivityUnits( specificActivity, 3, useCurrie );
    information.push_back("Specific Act: " + sa );
  }//if( not stable )
  
  //the non metastable child (if any)
  const SandiaDecay::Nuclide * metaChild = NULL;
  
  for( const SandiaDecay::Transition * transition : nuc->decaysToChildren )
  {
    if( transition->child )
    {
      snprintf( buffer, sizeof(buffer),
               "Decays to %s by %s decay, BR %.5f",
               transition->child->symbol.c_str(),
               SandiaDecay::to_str( transition->mode ),
               transition->branchRatio );
      
      //strip off dangling zeros...
      string val = buffer;
      while( val.length() > 2 && val[val.length()-1]=='0' && val[val.length()-2]!='.' )
        val = val.substr(0,val.length()-1);
      
      information.push_back( val );
      
      if( nuc->massNumber == transition->child->massNumber
         && nuc->atomicNumber == transition->child->atomicNumber )
      {
        metaChild = transition->child;
      }
    }//if( transition->child )
  }//for( const SandiaDecay::Transition * transition : nuc->decaysToChildren)
  
  if( includeDaughterIsomerics && metaChild )
  {
    information.push_back("");
    information.push_back(metaChild->symbol);
    
    char anStr[64], amStr[64];
    snprintf( anStr, sizeof(anStr), "Atomic Number: %i", metaChild->atomicNumber );
    snprintf( amStr, sizeof(amStr), "Atomic Mass: %.f", metaChild->atomicMass );
    const string hl = PhysicalUnits::printToBestTimeUnits( metaChild->halfLife );
    
    information.push_back( "Half Life: " + hl );
    
    if( parentNuclide && (metaChild != parentNuclide) )
    {
      snprintf( buffer, sizeof(buffer), "Branch Ratio from %s: %.5g",
               parentNuclide->symbol.c_str(),
               parentNuclide->branchRatioToDecendant(metaChild) );
      information.push_back( buffer );
    }//if( nuc == parentNuclide )
    
    
    const double metaSpecificActivity = metaChild->activityPerGram() / PhysicalUnits::gram;
    const string metaSA = PhysicalUnits::printToBestSpecificActivityUnits( metaSpecificActivity, 3, useCurrie );
    information.push_back("Specific Act: " + metaSA );
    
    for(const SandiaDecay::Transition * transition : metaChild->decaysToChildren)
    {
      if( transition->child )
      {
        snprintf( buffer, sizeof(buffer),
                 "Decays to %s by %s decay, BR %.5f",
                 transition->child->symbol.c_str(),
                 SandiaDecay::to_str( transition->mode ),
                 transition->branchRatio );
        information.push_back( buffer );
      }//if( transition->child )
    }//for(const SandiaDecay::Transition * transition : metaChild->decaysToChildren)
    
    //add a blank line at the end to tell we have a lot of info to print
    information.push_back("");
  }//if( metaChild )
  
  return information;
}//std::vector<std::string> getTextInfoForNuclide( const SandiaDecay::Nuclide *nuc )


void DecayChainChartD3::defineJavaScript()
{
  m_jsLoaded = true;
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  const bool isMobile = app && app->isMobile();
  InterSpec *interspec = app ? app->viewer() : nullptr;
  
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
  setJavaScriptMember( "wtResize", "function(self, w, h, layout){ " + jsRef() + ".chart.handleResize();}" );
  
  for( const string &js : m_pendingJs )
    WContainerWidget::doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void DecayChainChartD3::showDecaysThrough( const std::string nuc )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nucptr = db->nuclide( nuc );
  
  //cout << "DecayChainChartD3::showDecaysThrough: '" << nuc << "'" << endl;
  
  showPossibleParents( nucptr );
}//void showDecaysThrough( const std::string nuc )


void DecayChainChartD3::showDecayParticleInfo( const std::string &csvIsotopeNames )
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
  m_moreInfoDialog->finished().connect( this, &DecayChainChartD3::deleteMoreInfoDialog );
  WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter("Close");
  ok->clicked().connect( this, &DecayChainChartD3::deleteMoreInfoDialog );
  
  int nucs = 0;
  for( const auto &nuc_infos : nucinfos )
  {
    auto header = new Wt::WText( "Particle from " + nuc_infos.first, m_moreInfoDialog->contents() );
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


void DecayChainChartD3::showPossibleParents( const SandiaDecay::Nuclide *nuclide )
{
  if( !nuclide )
    return;
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  InterSpec *interspec = app ? app->viewer() : nullptr;
  if( !interspec )  //shouldnt ever happen, but jic
    return;
  
  deleteMoreInfoDialog();
  
  const double ww = 0.8*interspec->renderedWidth();
  const double wh = 0.65*interspec->renderedHeight();
  
  Wt::WFlags<AuxWindowProperties> windowProp
  = Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
  | AuxWindowProperties::DisableCollapse
  | AuxWindowProperties::SetCloseable
  | AuxWindowProperties::EnableResize;
  
  m_moreInfoDialog = new AuxWindow( "Decays Through " + nuclide->symbol, windowProp );
  DecayChainChartD3 *w = new DecayChainChartD3();
  w->setNuclide( nuclide, m_useCurrie, DecayChainChartD3::DecayChainType::DecayThrough );
  m_moreInfoDialog->stretcher()->addWidget( w, 0, 0  );
  
  m_moreInfoDialog->resizeWindow( std::max( 800, static_cast<int>(std::min(ww,1.3*wh) ) ), std::max( 600, static_cast<int>(wh) ) );
  
  
  m_moreInfoDialog->rejectWhenEscapePressed();
  m_moreInfoDialog->finished().connect( this, &DecayChainChartD3::deleteMoreInfoDialog );
  WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter("Close");
  ok->clicked().connect( this, &DecayChainChartD3::deleteMoreInfoDialog );
  
  
  m_moreInfoDialog->show();
  m_moreInfoDialog->centerWindowHeavyHanded();
}//void showPossibleParents(...)


void DecayChainChartD3::deleteMoreInfoDialog()
{
  if( m_moreInfoDialog )
    delete m_moreInfoDialog;
  m_moreInfoDialog = nullptr;
}//void deleteMoreInfoDialog()


#else //USE_D3_DECAY_CHAIN_IMP

DecayChainChart::DecayChainChart(  Wt::WContainerWidget *parent )
  : Wt::WPaintedWidget( parent ),
    m_axisBuffer( 22 ),
    m_infoNuclide( nullptr ),
    m_nuclide( nullptr ),
    m_useCurrie( true ),
    m_moreInfoDialog( NULL ),
    m_tapnhold( this, "tapnhold" )
{
  wApp->useStyleSheet( "InterSpec_resources/DecayChainChart.css" );
  
  addStyleClass( "DecayChainChart" );
  m_tapnhold.connect( this, &DecayChainChart::makeDialog );
}//DecayChainChart( constructor )

DecayChainChart::~DecayChainChart() {}
  

/*gets the x coordinate to put the box
  params:
	x = number atomic mass
	minx = 0
	xRange = number of atomic masses
	windowWidth = width of the panel for the diagram (not the whole available panel)
*/
double DecayChainChart::getXCoordRect(int x, int minx, int xRange, int windowWidth)
{
  if( xRange == 1 )
    return windowWidth / 2.0 - 90;
  else
    return m_axisBuffer + windowWidth * (1 - (1.0 * x - minx) / xRange) - 90;
}

/*gets the y coordinate to put the box
  params:
	y = atomicNumber
	miny = min atomicNumber
	yRange = range of atomic numbers
	windowHeight = height of the panel for the diagram (not the whole available panel)
*/
double DecayChainChart::getYCoordRect(int y, int miny, int yRange, int windowHeight)
{
  if( yRange <= 0 )
    return 0.5 * windowHeight;
  return windowHeight * (1.0 - (static_cast<double>(y - miny) / yRange)) + 5.0;
}

/*gets the x coordinate for the arrow on the left side
  params:
	x = number atomic mass
	minx = 0
	xRange = number of atomic masses
	windowWidth = width of the panel for the diagram (not the whole available panel)
*/
double DecayChainChart::getXCoordArrowLeft(int x, int minx, int xRange, int windowWidth)
{
  return getXCoordRect(x, minx, xRange, windowWidth) + 10;
}

/*gets the y coordinate for the arrow
  params:
	y = atomicNumber
	miny = min atomicNumber
	yRange = range of atomic numbers
	windowHeight = height of the panel for the diagram (not the whole available panel)
*/
double DecayChainChart::getYCoordArrow(int y, int miny, int yRange, int windowHeight)
{
  return getYCoordRect(y, miny, yRange, windowHeight) + 7.5;
}

/*gets the x coordinate for the arrow on the right side
  params:
	x = number atomic mass
	minx = 0
	xRange = number of atomic masses
	windowWidth = width of the panel for the diagram (not the whole available panel)
*/
double DecayChainChart::getXCoordArrowRight(int x, int minx, int xRange, int windowWidth)
{
  return getXCoordRect(x, minx, xRange, windowWidth) + 65;
}

/*draws the triangle for the arrow
  params:
    parX = location of parent x coordinate
    parY = location of parent y coordinate
    childX = location of child x coordinate
    childY = location of child y coordinate
    ptr = pointer to array of Wpoints that has already been allocated
*/
void DecayChainChart::getTrianglePts( double parX, double parY,
                                         double childX, double childY,
                                         Wt::WPointF ptr[3] )
{
  ptr[2] = Wt::WPointF(childX, childY);

  //find the angle of the arrow
  const double deltaX = parX - childX;
  const double deltaY = parY - childY;
  const double angle = atan2(deltaY, deltaX);
  const double arrowLength = 7;
  const double x1 = childX + arrowLength * cos(angle + .5);
  const double y1 = childY + arrowLength * sin(angle + .5);
  const double x2 = childX + arrowLength * cos(angle - .5);
  const double y2 = childY + arrowLength * sin(angle - .5);

  ptr[0] = Wt::WPointF(x1, y1);
  ptr[1] = Wt::WPointF(x2, y2);
}

/*draws the mouseover text
  params:
    info = the informational text
*/
void DecayChainChart::displayTextInfo( const SandiaDecay::Nuclide *nuc )
{
  m_infoNuclide = nuc;
  update();
}



void DecayChainChart::showPossibleParents( const SandiaDecay::Nuclide *nuclide )
{
  if( !nuclide )
    return;

  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  InterSpec *interspec = app ? app->viewer() : nullptr;
  if( !interspec )  //shouldnt ever happen, but jic
    return;
  
  deleteMoreInfoDialog();
  
  const double ww = 0.8*interspec->renderedWidth();
  const double wh = 0.65*interspec->renderedHeight();
  
  Wt::WFlags<AuxWindowProperties> windowProp
    = Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
      | AuxWindowProperties::DisableCollapse
      | AuxWindowProperties::SetCloseable
      | AuxWindowProperties::EnableResize;
  
  m_moreInfoDialog = new AuxWindow( "Decays Through " + nuclide->symbol, windowProp );
  DecayChainChartD3 *w = new DecayChainChartD3();
  w->setNuclide( nuclide, m_useCurrie, DecayChainChartD3::DecayChainType::DecayThrough );
  m_moreInfoDialog->stretcher()->addWidget( w, 0, 0  );
  
  m_moreInfoDialog->resizeWindow( std::max( 800, static_cast<int>(std::min(ww,1.3*wh) ) ), std::max( 600, static_cast<int>(wh) ) );
  
  
  m_moreInfoDialog->rejectWhenEscapePressed();
  m_moreInfoDialog->finished().connect( this, &DecayChainChart::deleteMoreInfoDialog );
  WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter("Close");
  ok->clicked().connect( this, &DecayChainChart::deleteMoreInfoDialog );
  
  
  m_moreInfoDialog->show();
  m_moreInfoDialog->centerWindowHeavyHanded();
}//void showPossibleParents(...)


void DecayChainChart::setNuclide( const SandiaDecay::Nuclide *nuc, const bool useCurrie )
{
  m_useCurrie = useCurrie;
  
  if( nuc == m_nuclide )
    return;
  
  m_infoNuclide = nullptr;
  m_nuclide = nuc;
  update();
  
  m_nuclideChanged.emit( m_nuclide );
}//void setNuclide(...)


Wt::Signal<const SandiaDecay::Nuclide *> &DecayChainChart::nuclideChanged()
{
  return m_nuclideChanged;
}

const SandiaDecay::Nuclide *DecayChainChart::nuclide() const
{
  return m_nuclide;
}


void DecayChainChart::makeDialog( const std::string &csvIsotopeNames )
{
  vector<string> isotopenames;
  SpecUtils::split( isotopenames, csvIsotopeNames, "," );
  
  std::vector<std::string> types;
  std::vector<double> energies;
  std::vector<double> intensities;
  std::vector<std::string> names;
  
  for( const string &isoname : isotopenames )
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *nuc = db->nuclide( isoname );
    if( nuc )
      fill_information( nuc, types, energies, intensities, names );
  }//for( const string &isoname : isotopenames )
  
  if( m_moreInfoDialog )
  {
    m_moreInfoDialog->contents()->clear();
  }else
  {
    m_moreInfoDialog = new AuxWindow("Particle Energies");
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
  }//if( m_moreInfoDialog ) / else
    
  if( types.size() )
  {
    new Wt::WText( names[0], m_moreInfoDialog->contents() );
    WTable *table = new Wt::WTable( m_moreInfoDialog->contents() );
    table->addStyleClass( "DecayChainChartTable" );
    table->setHeaderCount( 1 );
    table->elementAt(0, 0)->addWidget( new Wt::WText("Particle") );
    table->elementAt(0, 1)->addWidget( new Wt::WText("Energy (keV)") );
    table->elementAt(0, 2)->addWidget( new Wt::WText("Intensity") );
    
    int offset = 0;
    for( int i = 0; i < static_cast<int>(types.size()); ++i )
    {
      //check to see if nuclide is metastable (will have two tables)
      if( types[i].compare("") == 0 )
      {
        //set up the new table
        new Wt::WText(names[1], m_moreInfoDialog->contents());
        table = new Wt::WTable(m_moreInfoDialog->contents());
        table->addStyleClass( "DecayChainChartTable" );
        table->setHeaderCount(1);
        table->elementAt(0, 0)->addWidget(new Wt::WText("Particle"));
        table->elementAt(0, 1)->addWidget(new Wt::WText("Energy (keV)"));
        table->elementAt(0, 2)->addWidget(new Wt::WText("num/decay"));
        offset = i + 1;
        
        continue;
      }//if (types[i].compare("") == 0)
      
      table->elementAt(i + 1 - offset, 0)->addWidget(new Wt::WText(types[i]));
      
      char energyStr[32], intensityStr[32];
      snprintf( energyStr, sizeof(energyStr), "%.2f", energies[i] );
      snprintf( intensityStr, sizeof(intensityStr), "%.5g", intensities[i] );
      table->elementAt(i + 1 - offset, 1)->addWidget( new Wt::WText(energyStr) );
      table->elementAt(i + 1 - offset, 2)->addWidget( new Wt::WText(intensityStr) );
    }//for( size_t i = 0; i < types.size(); ++i )
  }else
  {
    const char *msg = "<center>This nuclide has no particles from decay</center>";
    new Wt::WText( msg, Wt::XHTMLText, m_moreInfoDialog->contents() );
  }//if( types.size() ) / else

//  WPushButton *ok = new WPushButton( "OK", m_moreInfoDialog->contents() );
 
  m_moreInfoDialog->show();
}//void makeDialog( const std::string &csvIsotopeNames )


void DecayChainChart::deleteMoreInfoDialog()
{
  if( m_moreInfoDialog )
  {
    delete m_moreInfoDialog;
    m_moreInfoDialog = (AuxWindow *)0;
  }//if( m_moreInfoDialog )
}//void deleteMoreInfoDialog()





void DecayChainChart::paintEvent( WPaintDevice *paintDevice )
{
  WPainter painter( paintDevice );
  
  //clear the interactive areas
  for( Wt::WAbstractArea * area : this->areas() )
    delete area;

  //Might consider also calling InterSpec::isMobile() since its a little
  //  computationally more efficient (but small change)
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  const bool isMobile = app && app->isMobile();
  
    
  const int chart_width  = static_cast<int>( floor(0.5+painter.window().width()) );
  const int chart_height = static_cast<int>( floor(0.5+painter.window().height()) );
  const int workingWidth = chart_width - m_axisBuffer;
  const int workingHeight = chart_height - 15 - m_axisBuffer - 10;
  
  //draw some instructions!
  if( isMobile )
  {
    painter.drawText( chart_width - 300, 5, 300, 15, Wt::AlignLeft,
                      "Tap nuclide for more information.");
    painter.drawText( chart_width - 300, 20, 300, 15, Wt::AlignLeft,
                      "Touch nuclide for 1 second for particles it gives off.");
    painter.drawText( chart_width - 300, 35, 300, 15, Wt::AlignLeft,
                      "Click a different source nuclide above to");
    painter.drawText( chart_width - 300, 50, 300, 15, Wt::AlignLeft,
                      "change displayed decay chain.");
  }else
  {
    painter.drawText( chart_width - 300, 5, 300, 15, Wt::AlignLeft,
                      "Click nuclide for more information.");
    painter.drawText( chart_width - 300, 20, 300, 15, Wt::AlignLeft,
                      "Double-click nuclide for particles it gives off.");
    painter.drawText( chart_width - 300, 35, 300, 15, Wt::AlignLeft,
                      "Click a different source nuclide above to");
    painter.drawText( chart_width - 300, 50, 300, 15, Wt::AlignLeft,
                      "change displayed decay chain.");
  }//if( isMobile ) / else
  
  //check to see if we have a nuclide, or if we shouldn't draw anything
  if( m_nuclide == NULL )
    return;
  
  //A list to find the axes that should be used for the data
  std::list<const SandiaDecay::Nuclide *> axisFindingList;
  axisFindingList.push_back( m_nuclide );

  //A list of nuclides to draw once the axes are found
  std::list<const SandiaDecay::Nuclide *> toMakeList;

  //A list of the atomic masses to use
  std::list<int> atomicMasses;
  
  short int minAN = 999, maxAN = -999;
  const SandiaDecay::Nuclide * lastDescendant = NULL;
  while( !axisFindingList.empty() )
  {
    const SandiaDecay::Nuclide * parent = axisFindingList.front();
    axisFindingList.pop_front();
    for( const SandiaDecay::Nuclide *descendant : parent->descendants() )
    {
      //check if we have not already looked at this nuclide
      if( !std::count(toMakeList.begin(), toMakeList.end(), descendant) )
      {
        //add the atomic mass to our list
        atomicMasses.push_back(descendant->massNumber);	

        minAN = std::min( minAN, descendant->atomicNumber );
        maxAN = std::max( maxAN, descendant->atomicNumber );

        //then look at the nuclide's children
        axisFindingList.push_back( descendant );
        
        //also, make sure we keep the data on the nuclide
        toMakeList.push_back( descendant );
        lastDescendant = descendant;
      }//if( !std::count(toMakeList.begin(), toMakeList.end(), descendant) )
    }//for( const SandiaDecay::Nuclide *descendant : parent->descendants() )
  }//while( !axisFindingList.empty() )

  atomicMasses.sort();
  atomicMasses.unique();
  
  //You now know the range of atomic numbers and atomic masses you'll have to draw
  //  from minAN/maxAN and the number of masses
  const int xRange = static_cast<int>( atomicMasses.size() );
  const int yRange = maxAN - minAN;
  
  for( const SandiaDecay::Nuclide * nuc : toMakeList )
  {
    const Wt::WString text = nuc->symbol;
    
    const double text_height = 15, text_width = 75;
    
    WPen gridPen;
    gridPen.setStyle( SolidLine );
    gridPen.setColor(Wt::WColor("darkgray"));
    painter.setPen(gridPen);

    int i = 0;
    //xRange is the number of atomic masses
    for( int mass : atomicMasses )
    {
      if( nuc->massNumber == mass )
	      break;
      i++;
    }

    const double text_x = getXCoordRect(i, 0, xRange, workingWidth);
    const double text_y = getYCoordRect(nuc->atomicNumber, minAN, yRange, workingHeight);
    painter.setBrush( WBrush(Wt::white) );
    painter.drawRect( text_x, text_y, text_width, text_height );
    WRectArea *infoBox = new WRectArea( text_x, text_y, text_width, text_height );
    this->addArea( infoBox );
    
    infoBox->addStyleClass( "DecayChainChartIsotope" );

    infoBox->clicked().connect(
         boost::bind(&DecayChainChart::displayTextInfo, this, nuc ) );

    gridPen.setColor(Wt::WColor("crimson"));
    painter.setPen( gridPen );
    painter.drawText( text_x, text_y,
                      text_width, text_height, Wt::AlignCenter, text );
    
    infoBox->doubleClicked().connect(
              boost::bind(&DecayChainChart::makeDialog, this, nuc->symbol) );
    
    //we have to detect touch events, and there is no double tap event call,
    //  so instead we'll require the user to hold their finger down for ~600 ms
    //  (we cant emulate double tap events since the chart is redrawn after
    //   each tap event)
    //This currently doesnt seem to work on touch devices
    const string ondown = "function(s,e){"
                          "var id='" + id() + "';"
                          "var iso='" + nuc->symbol + "';"
                          + OnClickDownJs + "}";
    infoBox->mouseWentDown().connect( ondown );
    infoBox->mouseWentUp().connect( OnClickUpJs );

    painter.restore();
  }//for( const SandiaDecay::Nuclide * nuc : toMakeList )
  
  {//begin code block to draw axis
    painter.save();
    WPen axisPen( WColor( "black" ) );
    axisPen.setStyle( SolidLine );
    axisPen.setWidth(2);
    painter.setPen( axisPen );
    painter.setBrush( Wt::WBrush(WColor("black")) );
    
    const double triangleLen = 4.0;
    double linexstart = m_axisBuffer - 0.5*triangleLen;
    double linexend   = linexstart + std::min( 165, chart_width );
    double lineystart = chart_height - m_axisBuffer + 0.5*triangleLen;
    double lineyend   = chart_height - m_axisBuffer + 0.5*triangleLen;
    painter.drawLine( linexstart, lineystart, linexend, lineyend );
    
    WPointF trianglePts[3];
    trianglePts[0] = WPointF( linexend, lineyend + 0.5*triangleLen );
    trianglePts[1] = WPointF( linexend, lineyend - 0.5*triangleLen );
    trianglePts[2] = WPointF( linexend + triangleLen, lineyend );
    painter.drawPolygon( trianglePts, 3 );
    painter.drawText( linexstart + 5, lineystart+4, 22.0*8.0, 15.0,
                      Wt::AlignLeft, "Decreasing Atomic Mass" );
    
    linexstart = linexend = m_axisBuffer - 0.5*triangleLen;
    lineystart = chart_height - m_axisBuffer + 0.5*triangleLen;
    lineyend   = lineystart - std::min( 165.0, lineystart );
    painter.drawLine( linexstart, lineystart, linexend, lineyend );
    
    trianglePts[0] = WPointF( linexend + 0.5*triangleLen, lineyend );
    trianglePts[1] = WPointF( linexend - 0.5*triangleLen, lineyend );
    trianglePts[2] = WPointF( linexend, lineyend - triangleLen );
    painter.drawPolygon( trianglePts, 3 );
    painter.rotate( -90.0 );
    painter.drawText( WRectF(-lineystart + 5, 0, 24.0*8.0, 15.0),
                      Wt::AlignLeft, TextSingleLine,
                      "Increasing Atomic Number" );
    painter.restore();
  }//end code block to draw axis
  
  
  WPen gridPen;
  gridPen.setColor( WColor( "blue" ) );
  gridPen.setStyle( SolidLine );
  gridPen.setWidth(2);
  painter.setPen( gridPen );
  painter.setBrush(Wt::WBrush(WColor("blue")));
  
  
  const vector<string> information = DecayChainChartD3::getTextInfoForNuclide( m_infoNuclide, m_nuclide, true, m_useCurrie );
  
  if( m_infoNuclide && !information.empty() )
  {
    string info = "Decays through " + m_infoNuclide->symbol + "...";
    const auto oldpen = painter.pen();
    
    WPen pen(GlobalColor::blue);
    painter.setPen( pen );
    painter.drawText( chart_width - 300, 65 + 5, 8.0*info.size(),
                     15.0, Wt::AlignLeft, info );
    
    pen.setWidth(1);
    painter.setPen( pen );
    painter.drawLine( chart_width - 300.0, 83, chart_width - 300.0 + 0.9*(8.0*info.size() - 24), 83 );
    
    painter.setPen(oldpen);
    
    WRectArea *infoBox = new WRectArea( chart_width - 300.0, 65.0 + 3.0, 8.0*info.size(), 19.0 );
    this->addArea( infoBox );
    infoBox->clicked().connect( boost::bind(&DecayChainChart::showPossibleParents, this, m_infoNuclide) );
    infoBox->setCursor( Cursor::PointingHandCursor );
  }//if( !information.empty() )
  
  
  //add the box to display the mouseover text
  int line = 0;
  for( auto iter = information.rbegin(); iter != information.rend(); ++iter )
  {
    const string &info = *iter;
    
    if( !line && info.empty() )
      continue;
    double y = chart_height - m_axisBuffer - 15 - 9 - 15*(line++);
    painter.drawText( m_axisBuffer+9.0, y, 8.0*info.size(),
                      15.0, Wt::AlignLeft, info );
  }//backward loop over m_infoText
  
  for( const SandiaDecay::Nuclide *nuc : toMakeList )
  {
    int i = 0;
    //xRange is the number of atomic masses
    for( int mass : atomicMasses )
    {
      if (nuc->massNumber == mass)
	      break;
      i++;
    }//for( int mass : atomicMasses )

    for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
    {
      painter.save();
      if( trans->child )
      {
        //check to see if it's metastable/stable transition
        if( nuc->massNumber == trans->child->massNumber
             && nuc->atomicNumber == trans->child->atomicNumber )
        {
          Wt::WPointF trianglePts[3];
          int x = (int)getXCoordRect(i, 0, xRange, workingWidth);
          int y = (int)getYCoordRect(nuc->atomicNumber, minAN, yRange, workingHeight);

          getTrianglePts(x - 2, y + 13.5, x + 1, y + 11.5, trianglePts);
          for( int i = 0; i < 3; ++i )
            trianglePts[i].setY(trianglePts[i].y() + 1);
          
          painter.drawPolygon(trianglePts, 3);

          painter.setBrush(Wt::NoBrush);
          painter.drawArc(x - 17, y - 2.5, 20, 20, 0, 16*320);
          painter.setBrush(Wt::WBrush(Wt::black));
          continue;
        }//if( a matching isomer isotope )

        //first find which number mass this is
        int j = 0;
        for( const int mass : atomicMasses )
        {
          if( trans->child->massNumber == mass )
            break;
          j++;
        }

        //draw a line to the transition
        int parX = (int)getXCoordArrowRight(i, 0, xRange, workingWidth);
        int parY = (int)getYCoordArrow(nuc->atomicNumber, minAN, yRange, workingHeight);
        int childX = (int)getXCoordArrowLeft(j, 0, xRange, workingWidth);
        int childY = (int)getYCoordArrow(trans->child->atomicNumber, minAN, yRange, workingHeight);
        
        if (childX > parX)
          childX -= 1;
        else
          childX += 1;
        
        if (childY > parY)
          childY -= 1;
        else
          childY += 1;
        
        painter.drawLine(parX, parY, childX, childY);

        Wt::WPointF trianglePts[3];
        getTrianglePts(parX, parY, childX, childY, trianglePts);
        painter.drawPolygon(trianglePts, 3);

//        draw the type of transition - seems too clustered, will part for now
//        painter.drawText((parX + childX)/2.0, (parY + childY)/2.0 - 10, 10,
//                10, AlignCenter, Wt::WString::fromUTF8("\u0444"));
//        would also want to put BR
      }//if( trans->child )

      painter.restore();
    }//for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
  }//for( const SandiaDecay::Nuclide *nuc : toMakeList )
  
}//void paintEvent( WPaintDevice *paintDevice )

#endif //USE_D3_DECAY_CHAIN_IMP / else

