#ifndef TerminalWidget_h
#define TerminalWidget_h
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

#include <tuple>
#include <memory>
#include <string>
#include <sstream>
#include <iostream>

#include <Wt/WSignal>
#include <Wt/WTextArea>
#include <Wt/WPopupMenu>
#include <Wt/WPushButton>
#include <Wt/WContainerWidget>

#include "InterSpec/PopupDiv.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/TerminalModel.h"

//Forward declarations
class InterSpec;

namespace Wt
{
  class WLineEdit;
}


class TerminalWidget : public Wt::WContainerWidget
{
  public:
  TerminalWidget( InterSpec *viewer, Wt::WContainerWidget *parent = 0 );
  void focusText();
  
  void chartClicked( float energy, float counts, int pageX, int pageY);
  
  protected:
  void handleEnterKey(); //an example function that gets called when user pressed the enter key
  void handleButtonClickedSignal( std::string argument );  //example function that gets called via m_buttonClickedSignal with the psosible JS filtered input text
  void commandMenuItemSelected();
  void commandMenuSearchInput();
  void handleTextInput();
  void textAdded();
  
  protected:
  typedef std::tuple<Wt::WMenuItem*,std::string,std::string> MenuItemTuple;
    
  InterSpec *m_viewer;
  Wt::WTextArea  *m_enteredtxt;
  Wt::WLineEdit  *m_edit, *m_commandsearch;
  PopupDivMenu   *m_commandmenu;
  std::unique_ptr<Wt::JSignal<std::string> > m_buttonClickedSignal;
    
  std::string darkenState();
  std::string lightenState();
    
  std::string searchToRegexLiteral( const std::string& search );
  void addHeaderToCommandSearchList( Wt::WMenuItem* header );
    
  private:
  TerminalModel *m_model;
  bool m_darkenState = false;

  std::vector< MenuItemTuple > m_commandMenuItems;
};//class TerminalWidget


#endif  //TerminalWidget_h
