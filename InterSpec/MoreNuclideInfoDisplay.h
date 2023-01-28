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


// Forward declerations
namespace SandiaDecay
{
  struct Nuclide;
}


class MoreNuclideInfoDisplay : public Wt::WTemplate
{
public:
  MoreNuclideInfoDisplay( const SandiaDecay::Nuclide *const nuc,
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

protected:
  void setTemplateTxt();
  void showDecayChainChart();
  void showDecayThroughChart();

  const SandiaDecay::Nuclide *m_nuc;
};//class MoreNuclideInfoDisplay


class MoreNuclideInfoWindow : public SimpleDialog
{
public:
  MoreNuclideInfoWindow( const SandiaDecay::Nuclide *const nuc );

protected:
  MoreNuclideInfoDisplay *m_display;
};//class MoreNuclideInfoWindow

#endif //MoreNuclideInfoDisplay_h
