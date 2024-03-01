#ifndef MoreNuclideInfoDisplay_h
#define MoreNuclideInfoDisplay_h
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

#include <map>
#include <string>
#include <memory>
#include <vector>

#include <Wt/WTemplate>

#include "InterSpec/SimpleDialog.h"


// Forward declarations
namespace SandiaDecay
{
  struct Nuclide;
}


class MoreNuclideInfoDisplay : public Wt::WTemplate
{
public:
  /** Constructor for MoreNuclideInfoDisplay
   
   @par nuc The (initial) nuclide to display information for.
   @par displayTitle If the title (e.g., "More info on U238") should be displayed within the
        widget, or if you will maybe display your own title (e.x., in the title section of
        SimpleDialog)
   @par parent The standard Wt parent widget.
   */
  MoreNuclideInfoDisplay( const SandiaDecay::Nuclide * const nuc,
                         const bool displayTitle,
                          Wt::WContainerWidget *parent = nullptr);


  /** Renders WTemplate containing more information of the specified nuclide.

  Combines information from both MoreNucInfoDb (if available) and SandiaDecay.

  This is the content in the window shown when the "more info" button is
  clicked on the "Reference Photopeaks" tab.

  @param nuc The nuclide to display information for.
  @param history If non-emty, the option to go back to the previous nuclide
         will be given to the user.
  */
  void setNuclide( const SandiaDecay::Nuclide *const nuc,
                   std::vector<const SandiaDecay::Nuclide *> history );

  /** Signal emitted when the nuclide that information is being displayed for, is changed.
   
   Can be used, for example, to update the title, if you are displaying the title yourself.
   */
  Wt::Signal<const SandiaDecay::Nuclide *> &nuclideChanged();
  
protected:
  void setTemplateTxt();
  void showDecayChainChart();
  void showDecayThroughChart();

  const SandiaDecay::Nuclide *m_nuc;
  bool m_displayTitle;
  Wt::Signal<const SandiaDecay::Nuclide *> m_nuclideChanged;
};//class MoreNuclideInfoDisplay


class MoreNuclideInfoWindow : public SimpleDialog
{
public:
  MoreNuclideInfoWindow( const SandiaDecay::Nuclide *const nuc );

  void nuclideUpdated( const SandiaDecay::Nuclide *nuc );
protected:
  MoreNuclideInfoDisplay *m_display;
};//class MoreNuclideInfoWindow

#endif //MoreNuclideInfoDisplay_h
