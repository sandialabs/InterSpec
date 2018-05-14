#ifndef DecayChainChart_h
#define DecayChainChart_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <string>
#include <vector>

#include <Wt/WPainter>
#include <Wt/WPaintedWidget>

//Forward declarations
namespace Wt
{
  class WDialog;
  class WComboBox;
  class WPaintDevice;
}//namespace Wt

class AuxWindow;
namespace SandiaDecay{ struct Nuclide; }


class DecayChainChart : public Wt::WPaintedWidget
{
public:
  DecayChainChart(  Wt::WContainerWidget *parent = NULL  );
  virtual ~DecayChainChart();
  
  virtual void paintEvent( Wt::WPaintDevice *paintDevice );
  
  const SandiaDecay::Nuclide *nuclide() const;
  
  /*sets the starting nuclide for the decay chain
   params:
   nuc = the nuclide in question
   */
  void setNuclide( const SandiaDecay::Nuclide *nuc );
  
  //delete the more infor dialog dialog
  void deleteMoreInfoDialog();
  
  //nuclideChanged(): emmited in setNuclide(...), and only when nuclide actually
  //  changed
  Wt::Signal<const SandiaDecay::Nuclide *> &nuclideChanged();

  
protected:
  /*gets the x coordinate to put the box
    params:
	x = number atomic mass
	minx = 0
	xRange = number of atomic masses
	windowWidth = width of the panel for the diagram (not the whole available panel)
  */
  double getXCoordRect(int x, int minx, int xRange, int windowWidth);

  /*gets the y coordinate to put the box
    params:
	y = atomic number
	minx = minx atomic number
	xRange = range of atomic numbers
	windowHeight = height of the panel for the diagram (not the whole available panel)
  */
  double getYCoordRect(int y, int miny, int yRange, int windowHeight);

  /*gets the x coordinate to put the arrow on the left side
    params:
	x = number atomic mass
	minx = 0
	xRange = number of atomic masses
	windowWidth = width of the panel for the diagram (not the whole available panel)
  */
  double getXCoordArrowLeft(int x, int minx, int xRange, int windowWidth);

  /*gets the x coordinate to put the arrow on the right side
    params:
	x = number atomic mass
	minx = 0
	xRange = number of atomic masses
	windowWidth = width of the panel for the diagram (not the whole available panel)
  */
  double getXCoordArrowRight(int x, int minx, int xRange, int windowWidth);

  /*gets the y coordinate to put the box
    params:
	    y = atomic number
    	minx = minx atomic number
	    xRange = range of atomic numbers
	    windowHeight = height of the panel for the diagram (not the whole available panel)
  */
  double getYCoordArrow(int y, int miny, int yRange, int windowHeight);
 
  /*draws the triangle for the arrow
    params:
      parX = location of parent x coordinate
      parY = location of parent y coordinate
      childX = location of child x coordinate
      childY = location of child y coordinate
      ptr = pointer to array of Wpoints
  */
  void getTrianglePts( double parX, double parY,
                       double childX, double childY, Wt::WPointF ptr[3] );

  /*draws the mouseover text
    params:
      info = the informational text
  */
  void displayTextInfo(std::vector<std::string> info);
  
  //make a dialog
  void makeDialog( const std::string &csvIsotopeNames );

  int m_axisBuffer;
  std::vector<std::string> m_infoText;
  const SandiaDecay::Nuclide * m_nuclide;

  AuxWindow *m_moreInfoDialog;
  Wt::JSignal<std::string> m_tapnhold;

  Wt::Signal<const SandiaDecay::Nuclide *> m_nuclideChanged;
};//class DecayChainChart


#endif //DecayChainChart_h
