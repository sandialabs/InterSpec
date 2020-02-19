#ifndef DecayChainChart_h
#define DecayChainChart_h
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

//We will eventually remove the non-D3 based implementation, but for the moment
//  we will leave it in, but keep it seperate.
#define USE_D3_DECAY_CHAIN_IMP 1

#if( USE_D3_DECAY_CHAIN_IMP )
/** This class is intended to be the successor to DecayChainChart.
 This class does all the rending using D3, wheres the DecayChainChart
 uses Wt namtive functions.
 
 The hope is to make the charting a little more responsive and interactive.
 */
class DecayChainChartD3 : public Wt::WContainerWidget
{
public:
  enum class DecayChainType
  {
    DecayFrom,
    DecayThrough
  };
  
  DecayChainChartD3( WContainerWidget *parent = nullptr );
  
  virtual void doJavaScript( const std::string &js ) override;
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;
  void setNuclide( const SandiaDecay::Nuclide * const nuc, const bool useCurrie, const DecayChainType decayType );
  const SandiaDecay::Nuclide *nuclide() const;
  
  void colorThemeChanged();
  
  void jsonInfoForNuclide( const SandiaDecay::Nuclide * const nuc, Wt::WStringStream &js ) const;
  void setJsonForDecaysFrom( const SandiaDecay::Nuclide * const nuclide );
  void setJsonForDecaysThrough( const SandiaDecay::Nuclide * const nuclide );
  
  /** Emmited in setNuclide(...), when nuclide actually changes */
  Wt::Signal<const SandiaDecay::Nuclide *> &nuclideChanged();
  
  /** Returns lines of information for a nuclide to display to the user.
   
   @param infoNuclide The nuclide to produce information about.  Will also
   display information for all the isomer
   @param parentNuclide The ultimate parent nuclide the user was originally
   interested in (e.g. the nuclide the decay chain is drawn for); is
   used to give the branching ratio that the infoNuclide is decayed
   through.
   @param includeDaughterIsomerics Wether to also include information for
   isomeric daughters of this nuclide.
   @param useCurrie Wether to use Ci or Bq for displaying specific activities.
   */
  static std::vector<std::string> getTextInfoForNuclide( const SandiaDecay::Nuclide *infoNuclide,
                                                        const SandiaDecay::Nuclide *parentNuclide,
                                                        const bool includeDaughterIsomerics,
                                                        const bool useCurrie );
  
  void showDecayParticleInfo( const std::string &csvIsotopeNames );
  void showPossibleParents( const SandiaDecay::Nuclide *nuclide );
  void deleteMoreInfoDialog();
  
protected:
  void defineJavaScript();
  
  /** Shows the "decay through" decay chart in the m_moreInfoDialog.
   
   This function is what is called when the clientside emits the
   'ShowDecaysThrough' signal by double clicking or tap-and-holding.
   
   @arg nuc The name of the nuclide to show decays through.  Ex. "Th232",
        "U235m", etc.  If invalid nuclide name no action will be taken.
   */
  void showDecaysThrough( const std::string nuc );
  
private:
  /** Whether to use Curies or Bequerels for displaying specific activity. */
  bool m_useCurrie;
  
  /** Until the widget is fully rendered, the JavaScript wont have been fully
   defined; this variable keeps track of this so we know if the JS can be
   executated now, or stored in #m_pendingJs until things are rendered.
   */
  bool m_jsLoaded;
  
  /** Stores any JavaScript #DecayChainChartD3::doJavaScript would have done
   before widget has been fully rendered, so it can then be exucuted upon
   rendering.
   */
  std::vector<std::string> m_pendingJs;
  
  /** The nuclide whose decay chain is currently being shown for. */
  const SandiaDecay::Nuclide *m_nuclide;
  
  /** Dialog box used to show both the decay particle information (e.g., when
   user double clicks on a nuclide), or the "decay through" decay chain.
   Will be nullptr when dialog isnt currently showing.
   */
  AuxWindow *m_moreInfoDialog;
  
  /** Signal emitted by #DecayChainChartD3::setNuclide when nuclide is actually
   changed.
   */
  Wt::Signal<const SandiaDecay::Nuclide *> m_nuclideChanged;
  
  /** The 'ShowDecayParticleInfo' signal emitted from the client side when user
   double-clicks, or taps-and-holds a nuclide on the decay chart.
   The string argument is a CSV list of nuclides to show information for.
   Ex., for nuclides with meta-states wil be like 'Pa234,Pa234m'
   */
  Wt::JSignal<std::string> m_showDecayParticles;
  
  /** The 'ShowDecaysThrough' signal emitted from the client side when user
   wants to see what nuclides will decay through the given nuclide.
   */
  Wt::JSignal<std::string> m_showDecaysThrough;
};//DecayChainChartD3

#else //USE_D3_DECAY_CHAIN_IMP

class DecayChainChart : public Wt::WPaintedWidget
{
public:
  DecayChainChart(  Wt::WContainerWidget *parent = NULL  );
  virtual ~DecayChainChart();
  
  virtual void paintEvent( Wt::WPaintDevice *paintDevice );
  
  const SandiaDecay::Nuclide *nuclide() const;
  
  /** Set the starting nuclide for the decay chain.
   @param nuc The nuclide to show decay for
   @param useCurrie Wether to use curries or bequerels for displaying activity
          related quatities, like specific activity; this value should come from
          what the user has entered the nuclide activity as, or similar
   */
  void setNuclide( const SandiaDecay::Nuclide *nuc, const bool useCurrie );
  
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
  void displayTextInfo( const SandiaDecay::Nuclide *nuc );
  
  /** Shows a pop-up displaying all nuclides that decay through the specified
   nuclide.
   */
  void showPossibleParents( const SandiaDecay::Nuclide *nuc );
  
  //make a dialog
  void makeDialog( const std::string &csvIsotopeNames );

  int m_axisBuffer;
  
  /** The nuclide the user clicked on that we will display additional
   information about.
  */
  const SandiaDecay::Nuclide *m_infoNuclide;
  
  /** The nuclide to show decay information for. */
  const SandiaDecay::Nuclide *m_nuclide;

  bool m_useCurrie;
  
  AuxWindow *m_moreInfoDialog;
  Wt::JSignal<std::string> m_tapnhold;

  Wt::Signal<const SandiaDecay::Nuclide *> m_nuclideChanged;
};//class DecayChainChart
#endif //USE_D3_DECAY_CHAIN_IMP

#endif //DecayChainChart_h
