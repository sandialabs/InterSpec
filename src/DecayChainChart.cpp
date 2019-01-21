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
#include <math.h>
#include <ctype.h>
#include <algorithm>

#include <Wt/WText>
#include <Wt/WTable>
#include <Wt/WString>
#include <Wt/WDialog>
#include <Wt/WLength>
#include <Wt/WRectArea>
#include <Wt/WPainter>
#include <Wt/WTableCell>
#include <Wt/Chart/WAxis>
#include <Wt/WPushButton>
#include <Wt/WPaintDevice>
#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesianChart>

#include <boost/ref.hpp>
#include <boost/any.hpp>
#include <boost/bind.hpp>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayChainChart.h"
#include "InterSpec/DecayActivityDiv.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

namespace
{
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
  
  
  const char *particlename( SandiaDecay::RadParticle par )
  {
    switch( par.type )
    {
      case SandiaDecay::AlphaParticle:           return "alpha";            break;
      case SandiaDecay::CaptureElectronParticle: return "electron capture"; break;
      case SandiaDecay::BetaParticle:    return "beta";             break;
      case SandiaDecay::GammaParticle:   return "gamma";            break;
      case SandiaDecay::PositronParticle:        return "positron";         break;
      case SandiaDecay::XrayParticle:           return "xray";             break;
      default:                           return "unknown";          break;
    }//switch( par.type )
    
    return "";
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
    const SandiaDecay::Nuclide * metaChild = NULL;
    
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
    if( metaChild )
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
}//namespace


DecayChainChart::DecayChainChart(  Wt::WContainerWidget *parent )
  : Wt::WPaintedWidget( parent ),
    m_axisBuffer( 22 ),
    m_nuclide( NULL ),
    m_moreInfoDialog( NULL ),
    m_tapnhold( this, "tapnhold" )
{
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
void DecayChainChart::displayTextInfo(vector<string> info)
{
  m_infoText = info;
  update();
}
/*sets the starting nuclide for the decay chain
  params:
    nuc = the nuclide in question
*/
void DecayChainChart::setNuclide( const SandiaDecay::Nuclide *nuc )
{
  if( nuc == m_nuclide )
    return;
  
  m_infoText.clear();
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
  UtilityFunctions::split( isotopenames, csvIsotopeNames, "," );
  
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
  char buffer[256];
  WPainter painter( paintDevice );
  
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
  
  //clear the interactive areas
  for( Wt::WAbstractArea * area : this->areas() )
    delete area;

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

    vector<string> information;
    information.push_back( nuc->symbol );
    snprintf( buffer, sizeof(buffer), "Atomic Number: %i", nuc->atomicNumber );
    information.push_back( buffer );
    snprintf( buffer, sizeof(buffer), "Atomic Mass: %.2f", nuc->atomicMass );
    information.push_back( buffer );
    
    if( nuc == m_nuclide )
    {
      if( m_nuclide->canObtainSecularEquilibrium() )
        information.push_back("Can reach secular equilibrium");
      else
        information.push_back("Cannot reach secular equilibrium");

      snprintf( buffer, sizeof(buffer), "Branch Ratio to %s: %.5g",
                lastDescendant->symbol.c_str(),
                nuc->branchRatioToDecendant(lastDescendant) );
      information.push_back( buffer );
    }//if( nuc == m_nuclide )
    
    if( IsInf(nuc->halfLife) )
    {
      information.push_back("Half Life: stable" );
    }else
    {
      const string hl = PhysicalUnits::printToBestTimeUnits( nuc->halfLife );
      information.push_back("Half Life: " + hl );
    }
    
    //the non metastable child (if any)
    const SandiaDecay::Nuclide * metaChild = NULL;
    
    for( const SandiaDecay::Transition * transition : nuc->decaysToChildren )
    {
      if( transition->child )
      {
        snprintf( buffer, sizeof(buffer),
                  "Decays to %s by %s decay, with branching ratio %.5f",
                  transition->child->symbol.c_str(),
                  SandiaDecay::to_str( transition->mode ),
                  transition->branchRatio );
        information.push_back( buffer );
          
        if( nuc->massNumber == transition->child->massNumber
            && nuc->atomicNumber == transition->child->atomicNumber )
        {
          metaChild = transition->child;
        }
      }//if( transition->child )
    }//for( const SandiaDecay::Transition * transition : nuc->decaysToChildren)
    
    if( metaChild )
    {
      information.push_back("");
      information.push_back(metaChild->symbol);
      
      char anStr[64], amStr[64];
      snprintf( anStr, sizeof(anStr), "Atomic Number: %i", metaChild->atomicNumber );
      snprintf( amStr, sizeof(amStr), "Atomic Mass: %.f", metaChild->atomicMass );
      string hl = PhysicalUnits::printToBestTimeUnits( metaChild->halfLife );
      
      information.push_back( "Half Life: " + hl );

      for(const SandiaDecay::Transition * transition : metaChild->decaysToChildren)
      {
        if( transition->child )
        {
          snprintf( buffer, sizeof(buffer),
                    "Decays to %s by %s decay, with branching ratio %.5f",
                    transition->child->symbol.c_str(),
                    SandiaDecay::to_str( transition->mode ),
                    transition->branchRatio );
          information.push_back( buffer );
        }//if( transition->child )
      }//for(const SandiaDecay::Transition * transition : metaChild->decaysToChildren)
      
      //add a blank line at the end to tell we have a lot of info to print
      information.push_back("");
    }//if( metaChild )
    
    infoBox->addStyleClass( "DecayChainChartIsotope" );

    infoBox->clicked().connect(
         boost::bind(&DecayChainChart::displayTextInfo, this, information) );

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
  
  
  //add the box to display the mouseover text
  int line = 0;
  
  for( auto iter = m_infoText.rbegin(); iter != m_infoText.rend(); ++iter )
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



