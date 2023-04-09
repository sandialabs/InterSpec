//
//  TerminalWidget.cpp
//  InterSpec
//
//  Created by Christian Kenneth Morte on 7/15/16.
//
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

#include <regex>
#include <string>
#include <iostream>

#include <Wt/WText>
#include <Wt/WAnchor>
#include <Wt/WString>
#include <Wt/WLineEdit>
#include <Wt/WTextArea>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>
#include <Wt/WPopupMenu>
#include <Wt/WStackedWidget>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/TerminalWidget.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/HelpSystem.h"

#include "js/TerminalWidget.js"


// The regex in GCC 4.8.x does not have working regex
#if( defined(__GLIBCXX__) && (__cplusplus < 201402L) )
static_assert( defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) \
              || defined(_GLIBCXX_REGEX_STATE_LIMIT) \
              || (defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE > 4), "GCC 4.8 is not supported due to buggy regex implementation" );
#endif


using namespace std;
using namespace Wt;

/*
 * See also: http://www.webtoolkit.eu/wt/blog/2010/03/02/javascript_that_is_c__
 */
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

TerminalWidget::TerminalWidget( InterSpec *viewer, Wt::WContainerWidget *parent )
   : Wt::WContainerWidget( 0 ),
     m_viewer( viewer ),
     m_enteredtxt( 0 ),
     m_edit( 0 )
{
  m_model = new TerminalModel( m_viewer);

  addStyleClass( "TerminalWidget" );
  
  //wApp->require( "InterSpec_resources/assets/js/yourjs.js" );
    
  wApp->useStyleSheet( "InterSpec_resources/TerminalWidget.css" );
  //wApp->useStyleSheet( "InterSpec_resources/InterSpec.css" );

  WGridLayout *layout = new WGridLayout();
  this->setLayout( layout );
    
  m_enteredtxt = new WTextArea();
  m_enteredtxt->addStyleClass( "historytxtdiv" );
  m_enteredtxt->setReadOnly( true );
  layout->addWidget( m_enteredtxt, 0, 0, 1, 3 );
  
  m_edit = new WLineEdit();
  
  m_edit->setAutoComplete( false );
  m_edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_edit->setAttributeValue( "autocorrect", "off" );
  m_edit->setAttributeValue( "spellcheck", "off" );
#endif

  m_edit->setId( "m_edit" );
  m_edit->setObjectName( "m_edit" );
  m_edit->setFocus( true );
  m_edit->setEmptyText( "Enter your command/expression here." );
  layout->addWidget( m_edit, 1, 1, 0, 1 );
    
  WPushButton *commandButton = new WPushButton();
  commandButton->setFirstFocus();                                                               // button to show pop-up menu for list of commands
  commandButton->setIcon( "InterSpec_resources/images/bullet_arrow_down.png" );
  layout->addWidget( commandButton, 1, 0 );
    
  m_commandsearch = new Wt::WLineEdit;                                             // search bar inside pop-up menu for list of commands
  m_commandsearch->setPlaceholderText( "Search for a command/function here." );
  m_commandsearch->setStyleClass("TerminalSearchBox");
  m_commandsearch->setMinimumSize(Wt::WLength::Auto, WLength(1.5,Wt::WLength::FontEm));
  m_commandsearch->setMargin(WLength(5,WLength::Pixel));
    
  m_commandmenu = new PopupDivMenu( commandButton, PopupDivMenu::TransientMenu );   // pop-up menu for list of commands
  m_commandmenu->addStyleClass( "command-menu" );
  m_commandmenu->addStyleClass( "command-menu-content" );
  m_commandmenu->addStyleClass( "command-menu-content a" );
  m_commandmenu->addStyleClass( "command-menu-content a:hover" );
  m_commandmenu->addStyleClass( "command-menu:hover .command-menu-content" );
  m_commandmenu->setHeight( 125 );
  m_commandmenu->addWidget( m_commandsearch );
    
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
    
  for ( const TerminalModel::CommandHelperTuple s : m_model->commandsFunctionsList() ) {            // initialize the command menu with commands
      const std::string& name = std::get<0>(s),
                         tags = std::get<1>(s),
                         toolTip = std::get<2>(s);
      
      if ( name == "<separator>" ) continue;
      else if ( name.find("<header:") != std::string::npos ) {
        m_commandmenu->addSeparator();
        m_commandMenuItems.push_back( MenuItemTuple(m_commandmenu->addSectionHeader( name.substr( 8 ) ), tags, toolTip ) );
      }
      else {
        Wt::WMenuItem* item = m_commandmenu->addItem( name );
        if ( !toolTip.empty() )
            HelpSystem::attachToolTipOn( item, toolTip, showToolTips );
        m_commandMenuItems.push_back( MenuItemTuple(item, tags, toolTip ) );
      }
  }
  
  m_commandsearch->setWidth( 325 );       // fix style classes / css for these elements
  commandButton->addStyleClass( "command-button" );
  commandButton->addStyleClass( "command-button:hover" );
    
  m_commandmenu->itemSelected().connect( this, &TerminalWidget::commandMenuItemSelected );
  m_commandsearch->textInput().connect( this, &TerminalWidget::commandMenuSearchInput );
  
  WPushButton *button = new WPushButton( "Enter" );
  layout->addWidget( button, 1, 2 );

  layout->setRowStretch( 0, 1 );
  layout->setColumnStretch( 1, 1 );
    
  
  //An example of calling javascript from c++
  const char *myjs = INLINE_JAVASCRIPT(
    console.log( "created TerminalWidget" );
  );
  doJavaScript( myjs );
  
  //An example of creating a javascript function you can    call client side (no
  //  round trip to server). See js/TerminalWidget.js for function
  //  implementation.
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetInit );
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetDarken);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetLighten);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetHandleEnterKey);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetSaveInput);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetAccessNextInput);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetAccessPreviousInput);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetSelectFirstArgument);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetCommandMenuItemSelected);
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetChartClicked );
  
  //An example of calling the above function, from the c++
  doJavaScript( "Wt.WT.TerminalWidgetInit();" );

  //An example of calling a client side JavaScript function, based on a
  //  client side interaction, without the round-trip to the server.
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetExampleClientJsFcn );
  button->clicked().connect( "function(sender,event){Wt.WT.TerminalWidgetExampleClientJsFcn('" + id() + "'," + m_edit->jsRef() + ");}" );
    
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetKeyPressed );
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsTerminalWidgetKeyDown );

  m_edit->enterPressed().connect( this, &TerminalWidget::handleEnterKey );
  m_edit->keyPressed().connect( "function(sender,event){Wt.WT.TerminalWidgetKeyPressed(" + m_edit->jsRef() + ", event);}" );
  m_edit->keyWentDown().connect("function(sender,event){Wt.WT.TerminalWidgetKeyDown(" + m_edit->jsRef() + ", event);}" );
  
  
  //An example of using javascript signal slot so you can process the entered
  //  text in javascript before uploading to the server (c++)
  LOAD_JAVASCRIPT(wApp, "js/TerminalWidget.js", "TerminalWidget", wtjsButtonClickedJsSlot );
  const string slotjs = "function(sender,event){Wt.WT.ButtonClickedJsSlot('" + id() + "'," + m_edit->jsRef() + ");}";
  button->clicked().connect( slotjs );
  m_buttonClickedSignal.reset( new JSignal<std::string>( this, "lineentered", true ) );
  m_buttonClickedSignal->connect( boost::bind( &TerminalWidget::handleButtonClickedSignal, this,
                                              boost::placeholders::_1 ) );
}//TerminalWidget

void TerminalWidget::focusText()
{
    m_edit->setFocus( true );
}//focusText()

void TerminalWidget::chartClicked( float energy, float counts, int pageX, int pageY )
{
    doJavaScript("Wt.WT.TerminalWidgetChartClicked(" + m_edit->jsRef() + ", " + std::to_string(energy) + "," + std::to_string(counts) + "," + std::to_string(pageX) + "," + std::to_string(pageY) + ");");
}

std::string TerminalWidget::searchToRegexLiteral( const std::string& search )
{
    std::string result;
    for ( const char c : search ) {
        if ( c == '^' || c == '$' || c == '.' || c == '*' || c == '+' || c == '?' ||c == '(' || c == ')' ||
             c == '[' || c == ']' || c == '{' || c == '}' || c == '|')
             result.append( std::string("\\") + c );
        else result.push_back(c);
    }
    return result;
}

void TerminalWidget::addHeaderToCommandSearchList( Wt::WMenuItem* header )
{
    if ( header->isSectionHeader() ) {
        m_commandmenu->addSeparator();
        m_commandmenu->addItem( header );
    }
}

std::string TerminalWidget::darkenState()
{
    if ( !m_darkenState ) {
        doJavaScript("Wt.WT.TerminalWidgetDarken(" + m_enteredtxt->jsRef() + ", " + m_edit->jsRef() + ");");
        m_darkenState = true;
        return "Successfully darkened background.";
    }
    return "Background is already darkened.";
}

std::string TerminalWidget::lightenState()
{
    if ( m_darkenState ) {
        doJavaScript("Wt.WT.TerminalWidgetLighten(" + m_enteredtxt->jsRef() + ", " + m_edit->jsRef() + ");");
        m_darkenState = false;
        
        return "Successfully lightened background.";
    }
    return "Background is already lightened.";
}


void TerminalWidget::handleEnterKey()
{
    // If m_enteredtxt is a WContainerWidget
    const std::string& input = m_edit->text().narrow();
    
    doJavaScript(m_enteredtxt->jsRef() + ".value += " + "(" + m_enteredtxt->jsRef() + ".value.length == 0 ? '' : '\\n') + '" + ">>> ' + " + Wt::WString( input ).jsStringLiteral() + "+ '" + (input.empty() ? "" : "\\n") + "';");

    std::string output = m_model->evaluate(input);
    if ( output == "[:darken:]" )
        output = darkenState();
    else if ( output == "[:lighten:]" )
        output = lightenState();

    else if ( output.find_first_of("[:saveFile:") == 0 ) {
        // do something to save file here
    }
    
    doJavaScript(m_enteredtxt->jsRef() + ".value += " + "'  ' + " + Wt::WString( output ).jsStringLiteral() + ";");
  
    m_edit->setText( "" );
    m_edit->setFocus( true );
    
    // Scroll down to bottom of text
    wApp->doJavaScript(m_enteredtxt->jsRef() + ".scrollTop += " + m_enteredtxt->jsRef() + ".scrollHeight;");

    m_edit->setText( "" );
    m_edit->setFocus( true );
    
}//void handleEnterKey()


//Example function called whenever you emit the 'lineentered' signal client side
void TerminalWidget::handleButtonClickedSignal( std::string argument )
{
  cout << "Received signal txt of: " << argument << endl;
    
    // If m_enteredtxt is a WContainerWidget
    const std::string& input = m_edit->text().narrow();
    
    doJavaScript(m_enteredtxt->jsRef() + ".value += " + "(" + m_enteredtxt->jsRef() + ".value.length == 0 ? '' : '\\n') + '" + ">>> ' + " + Wt::WString( input ).jsStringLiteral() + "+ '" + (input.empty() ? "" : "\\n") + "';");

    std::string output = m_model->evaluate(argument);
    if ( output == "[:darken:]" )
        output = darkenState();
    else if ( output == "[:lighten:]" )
        output = lightenState();
        
    else if ( output.find_first_of("[:saveFile:") == 0 ) {
        // do something to save file here
    }
    
    doJavaScript(m_enteredtxt->jsRef() + ".value += " + "'  ' + " + Wt::WString( output ).jsStringLiteral() + ";");
    
    m_edit->setText( "" );
    m_edit->setFocus( true );
    
    // Scroll down to bottom of text
    wApp->doJavaScript(m_enteredtxt->jsRef() + ".scrollTop += " + m_enteredtxt->jsRef() + ".scrollHeight;");
}//void handleButtonClickedSignal( std::string argument )

void TerminalWidget::commandMenuSearchInput()
{
    const std::string& search = m_commandsearch->text().narrow();
    for ( Wt::WMenuItem* item : m_commandmenu->items() )
        if ( m_commandmenu->indexOf( item ) != 0 )
            m_commandmenu->removeItem( item );

    bool headerMatched = false;
    for ( int index = 0; index < m_commandMenuItems.size(); ++index ) {
        Wt::WMenuItem* menuItem = std::get<0>( m_commandMenuItems.at(index) );
        
        const std::string& menuItemText = menuItem->text().narrow(),
                           tags = std::get<1>( m_commandMenuItems.at(index) ),
                           toolTip = std::get<2>( m_commandMenuItems.at(index) ),
                           searchRegex = searchToRegexLiteral( search );
        
        const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
        
        if ( !toolTip.empty() )
            HelpSystem::attachToolTipOn( menuItem, toolTip, showToolTips );
        
        const bool itemMatched = std::regex_search( menuItemText, std::regex( searchRegex, std::regex_constants::icase ) )
                                 || std::regex_search( tags, std::regex( searchRegex, std::regex_constants::icase ) );
        
        if ( headerMatched ) {                                                      // if header matched, add all items from
            if ( menuItem->isSectionHeader() ) {                                    // that header to the list of items
                m_commandmenu->addSeparator();
                headerMatched = false;
            }
            m_commandmenu->addItem( menuItem );
            continue;
        }
        
        if ( menuItem->isSectionHeader() && itemMatched ) {                         // if header matched, add the header
            headerMatched = true;
            addHeaderToCommandSearchList( menuItem );
            continue;
        }
        
        if (  itemMatched ) {
            for (int headerIndex = index; headerIndex >= 0; --headerIndex) {        // find the section header for matched item
                Wt::WMenuItem* currentItem = std::get<0>( m_commandMenuItems.at(headerIndex) );
                
                if ( currentItem->isSectionHeader() ) {
                    if ( m_commandmenu->indexOf(currentItem) == -1 )                // if header not present in search list
                        addHeaderToCommandSearchList( currentItem );                // add the header
                    break;
                }
            }
            m_commandmenu->addItem( menuItem );
        }
        
        for ( Wt::WMenuItem* item : m_commandmenu->items() ) {                      // cleaning up the search list
            if ( m_commandmenu->indexOf( item ) != 0 && item->isSectionHeader() ) { // removes empty headers
                int index = m_commandmenu->indexOf( item );
                
                if ( index == m_commandmenu->items().size()-1 || m_commandmenu->itemAt( index+1 )->isSeparator() ) {
                    Wt::WMenuItem* itemBefore = m_commandmenu->itemAt( index-1 );
                    if ( itemBefore && itemBefore->isSeparator() )  m_commandmenu->removeItem( itemBefore );
                    m_commandmenu->removeItem( item );
                }
            }
        }
    }
}

void TerminalWidget::commandMenuItemSelected()
{
    const Wt::WString& selectedCommand = m_commandmenu->result()->text();
    doJavaScript("Wt.WT.TerminalWidgetCommandMenuItemSelected(" + m_edit->jsRef() + ", " + selectedCommand.jsStringLiteral() + ");");
//    m_edit->setFocus();
}//void m_m_commandMenuItemselected()

