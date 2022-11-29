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


/**  Widget that renders the decay chain of a given nuclide.
 
 Uses d3.js to render in SVG.  The C++ code takes care of producing the JSON
 that gets sent to the client side to be rendered.
 
 @TODO Should send the raw nuclide information to client-side to generate
       the about text when a nuclide is clicked on, instead of doing it
       in the C++ function #DecayChainChart::getTextInfoForNuclide
 @TODO Should do all JSON organization in the C++ code, rather than reformatting
       the isomeric decays client side.
 @TODO Add branching ratio and decay type next to decay lines
 @TODO The instruction text, or nuclide info text can overlap with some of the
       nuclides, for some decay schemes - should fix.
 @TODO Make static HTML/JS example versions of this widget to include in non-Wt
       based webpages.
 @TODO For the "Show decays through ..." button should decide if the text should
       be underlined, or a box-border drawn around text.
 @TODO Add a SVG and/or PNG image download option
 */
class DecayChainChart : public Wt::WContainerWidget
{
public:
  enum class DecayChainType
  {
    DecayFrom,
    DecayThrough
  };//enum class DecayChainType
  
  DecayChainChart( WContainerWidget *parent = nullptr );
  
  /** Show the decay chain, or decay-through chart for a specified
  nuclide, as a AuxWindow.

   @param nuc Nuclide to show chart for. 
   @param type The type (decay-chain, or decay-through) of chart to show.

   @returns Pointers to the created AuxWindow and DecayChainChart; they could
            be nullptr.
  */
  static std::pair<AuxWindow *, DecayChainChart *>
    show_decay_chart_window( const SandiaDecay::Nuclide *const nuc,
      const DecayChainType type );


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
   @param useCurrie Wether to use Ci or Bq for displaying specific activities.
   */
  static std::vector<std::string>
    getTextInfoForNuclide( const SandiaDecay::Nuclide * const infoNuclide,
                           const SandiaDecay::Nuclide * const parentNuclide,
                           const bool useCurrie );
  
  void showDecayParticleInfo( const std::string &csvIsotopeNames );
  void showPossibleParents( const SandiaDecay::Nuclide *nuclide );
  void deleteMoreInfoDialog();
  
protected:
  /** Called during the first full render of the widget; loads all necassary
   javascript.
   */
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
  
  /** Stores any JavaScript #DecayChainChart::doJavaScript would have done
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
  
  /** Signal emitted by #DecayChainChart::setNuclide when nuclide is actually
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
};//DecayChainChart

#endif //DecayChainChart_h
