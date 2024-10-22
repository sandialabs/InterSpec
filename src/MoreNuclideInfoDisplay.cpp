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
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>


#include <Wt/WText>
#include <Wt/WTable>
#include <Wt/WAnchor>
#include <Wt/WString>
#include <Wt/WTemplate>
#include <Wt/WPushButton>
#include <Wt/WApplication>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayChainChart.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MoreNuclideInfoDisplay.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"


using namespace std;
using namespace Wt;

namespace
{
  /** Forms the notes into HTML by doing things  like replacing '&lt;' with '<', or multiple line breaks with <br />
   as well as taking care of references.
   */
  std::string notes_html( const std::shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> &db,
                          const MoreNuclideInfo::NucInfo * const more_info )
  {
    if( !db || !more_info )
      return "";
    
    string notes = more_info->m_notes;
    SpecUtils::trim( notes ); //just to make sure (although I think the XML parser already did this)

    if( notes.empty() )
      return notes;
    
    SpecUtils::ireplace_all( notes, "\r", "" );
    vector<string> more_info_lines;

    // Replace things like &lt;li&gt;Some List Item&lt;/li&gt;
    //  to <li>Some List Item</li>
    std::regex tag_regex( "&lt;([a-zA-Z/ ]{1,4})&gt;" );
    notes = std::regex_replace( notes, tag_regex,"<$1>");
    
    boost::algorithm::split( more_info_lines, notes,
      boost::algorithm::is_any_of( "\n" ),
      boost::algorithm::token_compress_off );
    
    vector<string> references;
    
    notes = "<div class=\"MoreNucInfoSubSec\">";
    for( size_t i = 0; i < more_info_lines.size(); ++i )
    {
      string line = more_info_lines[i];
      SpecUtils::trim( line );

      std::regex ref_regex( "\\[(.+?)\\]" );
      auto ref_begin = std::sregex_iterator(begin(line), end(line), ref_regex);
      auto ref_end = std::sregex_iterator();
      for( std::sregex_iterator i = ref_begin; i != ref_end; ++i )
      {
        const string match_str = i->str();
        auto prev_pos = std::find( begin(references), end(references), match_str );
        if( prev_pos != end(references) )
        {
          const auto index = 1 + (prev_pos - begin(references));
          const string num_ref = "[" + std::to_string(index) + "]";
          SpecUtils::ireplace_all( line, match_str.c_str(), num_ref.c_str() );
        }else
        {
          references.push_back( match_str );
          const auto index = references.size();
          const string num_ref = "[" + std::to_string(index) + "]";
          SpecUtils::ireplace_all( line, match_str.c_str(), num_ref.c_str() );
        }
      }//for( std::sregex_iterator i = ref_begin; i != ref_end; ++i )
      
      notes += line;

      if( line.empty() )
      {
        notes += "</div><div class=\"MoreNucInfoSubSec\">";

        // Skip all the next empty lines
        for( ; (i + 1) < more_info_lines.size(); ++i )
        {
          string nextline = more_info_lines[i + 1];
          SpecUtils::trim( nextline );
          if( !nextline.empty() )
            break;
        }//for( skip consequitive empty lines )
      }//if( line.empty() )
    }//for( size_t i = 0; i < more_info_lines.size(); ++i )
    
    if( !references.empty() )
    {
      notes += "<div class=\"MoreNucReferences\">";
      for( size_t index = 0; index < references.size(); ++index )
      {
        string ref = references[index];
        
        if( !ref.empty() && (ref.front() == '[') )
          ref = ref.substr(1);
        
        if( !ref.empty() && (ref.back() == ']') )
          ref = ref.substr(0, ref.size() - 1);
        
        assert( !ref.empty() );
        if( ref.empty() )
          continue;
        
        notes += "<div class=\"MoreNucRef\">";
        
        notes += "<div>[" + std::to_string(index + 1) + "]</div><div>";
        
        const auto pos = db->m_references.find(ref);
        if( pos == end(db->m_references) )
        {
          notes += "Reference '" + ref + "' not found.";
        }else
        {
          const MoreNuclideInfo::RefInfo &ref_info = pos->second;
          
          notes += "<ul>";
          if( !ref_info.m_desc.empty() )
            notes += "<li>" + ref_info.m_desc + "</li>";
          
          if( !ref_info.m_url.empty() )
            notes += "<li><a href=\"" + ref_info.m_url + "\" target=\"_blank\">" + ref_info.m_url + "</a></li>";
          
          notes += "</ul>";
        }//if( we didnt find refernce ) / else
        
        notes += "</div></div>";
      }//for( string &ref : references )
      
      notes += "</div>";
    }//if( !references.empty() )
    
    notes += "</div>";
    
    return notes;
  }//string notes_html( MoreNucInfoDb &db, NucInfo * )
  
}//namespace


MoreNuclideInfoDisplay::MoreNuclideInfoDisplay( const SandiaDecay::Nuclide *const nuc, 
                                               const bool display_title,
                                               Wt::WContainerWidget *parent )
  : WTemplate( parent ),
    m_nuc( nullptr ),
    m_displayTitle( display_title ),
    m_nuclideChanged( this ),
    m_decayWindow( nullptr )
{
  try
  {
    setTemplateTxt();
  }catch( std::exception &e )
  {
    setTemplateText( "Error loading MoreNuclideInfoDisplay.tmplt.xml" );
    return;
  }

  wApp->useStyleSheet( "InterSpec_resources/MoreNuclideInfoDisplay.css" );

  setNuclide( nuc, {} );
}//MoreNuclideInfoDisplay constructor


void MoreNuclideInfoDisplay::handleDecayChartClose( AuxWindow *window )
{
  assert( window == m_decayWindow );
  if( window != m_decayWindow )
    cerr << "MoreNuclideInfoDisplay::handleDecayChartClose(...): passed in pointer ("
         << window << ") not same as internal pointer (" << m_decayWindow << ")." << endl;
  
  m_decayWindow = nullptr;
}//void handleDecayChartClose( AuxWindow *window )


void MoreNuclideInfoDisplay::programmaticallyCloseDecayChart()
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  if( m_decayWindow )
  {
    // Calling done seems to crash for some reason
    //m_decayWindow->done( Wt::WDialog::DialogCode::Accepted );
    AuxWindow::deleteAuxWindow( m_decayWindow );
  }
  m_decayWindow = nullptr;
}//void programmaticallyCloseDecayChart()


void MoreNuclideInfoDisplay::implementShowDecayCharts( const bool through )
{
  const auto type = through ? DecayChainChart::DecayChainType::DecayThrough
                            : DecayChainChart::DecayChainType::DecayFrom;
  const SandiaDecay::Nuclide *nuc = m_nuc;
  
  if( m_decayWindow )
  {
    m_decayWindow->done( Wt::WDialog::DialogCode::Accepted );
    assert( !m_decayWindow );
  }
  
  const pair<AuxWindow *, DecayChainChart *> results
                         = DecayChainChart::show_decay_chart_window( nuc, type );
  
  m_decayWindow = results.first;
  
  if( m_decayWindow )
  {
    m_decayWindow->finished().connect( boost::bind(&MoreNuclideInfoDisplay::handleDecayChartClose,
                                                   this, m_decayWindow) );
  }
  
  auto get_display = []() -> MoreNuclideInfoDisplay * {
    InterSpec *interspec = InterSpec::instance();
    ReferencePhotopeakDisplay *display = interspec ? interspec->referenceLinesWidget() : nullptr;
    MoreNuclideInfoWindow *window = display ? display->moreInfoWindow() : nullptr;
    MoreNuclideInfoDisplay *infoDisplay = window ? window->display() : nullptr;
    assert( infoDisplay );
    return infoDisplay;
  };
  
  m_decayWindow->finished().connect( std::bind( [through, get_display](){
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      auto undo = [through, get_display](){
        MoreNuclideInfoDisplay *disp = get_display();
        if( disp )
          disp->implementShowDecayCharts( through );
      };
      
      auto redo = [through, get_display](){
        MoreNuclideInfoDisplay *disp = get_display();
        if( disp )
          disp->programmaticallyCloseDecayChart();
      };
      
      undoRedo->addUndoRedoStep( undo, redo, "" );
    }//if( disp && undoRedo && undoRedo->canAddUndoRedoNow() )
  } ));
  

  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [get_display](){
      MoreNuclideInfoDisplay *infoDisplay = get_display();
      if( infoDisplay )
        infoDisplay->programmaticallyCloseDecayChart();
    };
    
    auto redo = [through, get_display](){
      MoreNuclideInfoDisplay *infoDisplay = get_display();
      if( infoDisplay )
        infoDisplay->implementShowDecayCharts( through );
    };
    
    const string desc = "Show " + (nuc ? nuc->symbol : string("null")) + " decay chart.";
    undoRedo->addUndoRedoStep( undo, redo, desc );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//void implementShowDecayCharts( const bool through )


void MoreNuclideInfoDisplay::showDecayChainChart()
{
  if( m_nuc )
    implementShowDecayCharts( false );
}//void showDecayChainChart()
  

void MoreNuclideInfoDisplay::showDecayThroughChart()
{
  if( m_nuc )
    implementShowDecayCharts( true );
}//void showDecayThroughChart()


void MoreNuclideInfoDisplay::setTemplateTxt()
{
  assert( wApp );
  const string filename = "InterSpec_resources/static_text/MoreNuclideInfoDisplay.tmplt.xml";
  const string docroot = wApp ? wApp->docRoot() : "";
  const string file_path = SpecUtils::append_path( docroot, filename );

#ifdef _WIN32
  ifstream infile( SpecUtils::convert_from_utf8_to_utf16( file_path ).c_str(), ios::in | ios::binary );
#else
  ifstream infile( file_path.c_str(), ios::in | ios::binary );
#endif
  assert( infile.is_open() );
  if( !infile.is_open() )
    throw runtime_error( "Failed to open more_nuc_info_tmplt.xml" );

  std::stringstream buffer;
  buffer << infile.rdbuf();
  WString tmplt_txt = WString::fromUTF8( buffer.str() );
  
  setTemplateText( tmplt_txt, Wt::TextFormat::XHTMLText );
}//void setTemplateTxt();


void MoreNuclideInfoDisplay::setNuclide( const SandiaDecay::Nuclide *const nuc, 
                                         vector<const SandiaDecay::Nuclide *> history )
{
  using namespace MoreNuclideInfo;
  
  WTemplate &tmplt = *this;
  
  const SandiaDecay::Nuclide * const orig_nuc = m_nuc;
  vector<const SandiaDecay::Nuclide *> orig_history = m_current_history;
  vector<const SandiaDecay::Nuclide *> input_history = history;
  
  try
  {
    if( !nuc )
      throw std::runtime_error( "Invalid Nuclide" );
    
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      throw runtime_error( "Error getting DecayDataBaseServer" );
    
    const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    tmplt.setCondition( "display-title", m_displayTitle );
    tmplt.setCondition( "if-have-prev", !history.empty() );
    if( !history.empty() )
    {
      const SandiaDecay::Nuclide *prev = history.back();
      assert( prev );
      if( !prev )
        throw runtime_error( "Unexpected nullptr" );
      vector<const SandiaDecay::Nuclide *> prev_history = history;
      prev_history.resize( prev_history.size() - 1 );
      
      WAnchor *anchor = new WAnchor();
      anchor->setText( "&larr;" + prev->symbol );
      anchor->clicked().connect( boost::bind( &MoreNuclideInfoDisplay::setNuclide, this, prev, prev_history ) );
      tmplt.bindWidget( "link-to-prev", anchor );
    }//if( !history.empty() )
    
    history.push_back( nuc );
    
    const SandiaDecay::Element *el = db->element( nuc->atomicNumber );
    assert( el );
    string elname = el ? el->name : string();
    if( !elname.empty() )
      elname[0] = static_cast<char>(std::tolower( static_cast<unsigned char>(elname[0]) ));
    
    string elsymbol = el ? el->symbol : string();
    if( !elsymbol.empty() )
      elsymbol[0] = static_cast<char>(std::toupper( static_cast<unsigned char>(elsymbol[0]) ));
    
    string meta_str = "";
    if( nuc->isomerNumber == 1 )
      meta_str = "m";
    else if( nuc->isomerNumber > 1 )
      meta_str = "m" + std::to_string( (int)nuc->isomerNumber );
    
    tmplt.bindString( "symbol", nuc->symbol, Wt::TextFormat::PlainText );
    tmplt.bindString( "element", elname, Wt::TextFormat::PlainText );
    tmplt.bindString( "element-symbol", elsymbol, Wt::TextFormat::PlainText );
    tmplt.bindInt( "mass-number", nuc->massNumber );
    tmplt.bindInt( "atomic-number", nuc->atomicNumber );
    tmplt.bindInt( "num-neutrons", nuc->massNumber - nuc->atomicNumber );
    
    tmplt.setCondition( "if-meta-stable", (nuc->isomerNumber != 0) );
    tmplt.bindInt( "meta-state-number", nuc->isomerNumber );
    tmplt.bindString( "meta-state", meta_str );
    
    
    if( !nuc->isStable() )
    {
      const string halflife = PhysicalUnits::printToBestTimeUnits( nuc->halfLife, 4, SandiaDecay::second );
      tmplt.bindString( "half-life", halflife, Wt::TextFormat::XHTMLText );
      
      const double specificActivity = nuc->activityPerGram() / PhysicalUnits::gram;
      const string sa = PhysicalUnits::printToBestSpecificActivityUnits( specificActivity, 3, !useBq );
      tmplt.bindString( "specific-activity", sa, Wt::TextFormat::XHTMLText );
    }else
    {
      tmplt.bindString( "half-life", "stable", Wt::TextFormat::PlainText );
      tmplt.bindString( "specific-activity", "", Wt::TextFormat::PlainText );
    }
    
    const string amu = SpecUtils::printCompact( nuc->atomicMass, 6 ) + " amu";
    tmplt.bindString( "atomic-mass", amu, Wt::TextFormat::PlainText );
    
    double natural_abundance = -1;
    if( el )
    {
      for( const auto &n : el->isotopes )
      {
        if( n.nuclide == nuc )
        {
          natural_abundance = n.abundance;
          break;
        }
      }
    }//if( el )
    
    tmplt.setCondition( "if-is-natural", (natural_abundance > 0.0) );
    const string natabun = SpecUtils::printCompact( 100.0*natural_abundance, 4 ) + "%";
    tmplt.bindString( "natural-abundance", natabun, Wt::TextFormat::PlainText );
    
    WPushButton *decayChainBtn = new WPushButton( "Show Decay Chain" );
    decayChainBtn->addStyleClass( "LightButton DecayBtn" );
    decayChainBtn->clicked().connect( this, &MoreNuclideInfoDisplay::showDecayChainChart );
    tmplt.bindWidget( "decay-chain-btn", decayChainBtn );
    
    tmplt.setCondition( "if-has-parents", !nuc->decaysFromParents.empty() );
    if( !nuc->decaysFromParents.empty() )
    {
      WPushButton *decayThroughBtn = new WPushButton( "Show Decays Through" );
      decayThroughBtn->addStyleClass( "LightButton DecayBtn" );
      decayThroughBtn->clicked().connect( this, &MoreNuclideInfoDisplay::showDecayThroughChart );
      tmplt.bindWidget( "decay-through-btn", decayThroughBtn );
    }
    
    //We'll overide the text of SandiaDecay::to_str( SandiaDecay::DecayMode ) for a few types
    auto trans_mode_txt = []( const SandiaDecay::DecayMode mode ){
      switch( mode )
      {
        case SandiaDecay::ElectronCaptureDecay:
          return "e- Capture";
        case SandiaDecay::BetaPlusAndAlphaDecay:
          return "Beta+ and Alpha";
        case SandiaDecay::DoubleBetaDecay:
          return "2 Beta+";
        case SandiaDecay::DoubleElectronCaptureDecay:
          return "2 e- Capture";
        case SandiaDecay::SpontaneousFissionDecay:
          return "Spont. Fission";
        default:
          break;
      };//switch( transition->mode )
      
      return SandiaDecay::to_str( mode );
    };//trans_mode_txt lamda
    
    
    // TODO: WTemplate has a Wt::WTemplate::Functions::while_f function, but from the documentation it isnt totally clear how to use it; it would be nice to use this to generate rows in a table or something to list all the individual decays
    WTable *decayTable = new WTable();
    decayTable->setHeaderCount( 1, Wt::Horizontal );
    WTableCell *cell = decayTable->elementAt( 0, 0 );
    new WText( "Child Iso.", cell );
    cell = decayTable->elementAt( 0, 1 );
    new WText( "Decay Mode", cell );
    cell = decayTable->elementAt( 0, 2 );
    new WText( "Branch Ratio", cell );
    
    
    for( const SandiaDecay::Transition *transition : nuc->decaysToChildren )
    {
      const int row = decayTable->rowCount();
      
      cell = decayTable->elementAt( row, 0 );
      if( transition->child )
      {
        WAnchor *anchor = new WAnchor( cell );
        anchor->setText( transition->child->symbol );
        anchor->clicked().connect( boost::bind( &MoreNuclideInfoDisplay::setNuclide, this, transition->child, history ) );
        
        string txt = "(t<sub>1/2</sub>=" + PhysicalUnits::printToBestTimeUnits( transition->child->halfLife, 2 ) + ")";
        WText *t = new WText( txt, cell );
        t->addStyleClass( "HalfLife" );
      }else
      {
        new WText( "<div style=\"text-align: center;\">--</div>", cell );
      }//if( transition->child )
      
      cell = decayTable->elementAt( row, 1 );
      new WText( trans_mode_txt( transition->mode ), cell );
      
      cell = decayTable->elementAt( row, 2 );
      new WText( SpecUtils::printCompact( transition->branchRatio, 5 ), cell );
    }//for( const SandiaDecay::Transition * transition : nuc->decaysToChildren)
    
    tmplt.bindWidget( "decays-to-table", decayTable );
    
    
    if( !nuc->decaysFromParents.empty() )
    {
      WTable *decayFromTable = new WTable();
      decayFromTable->setHeaderCount( 1, Wt::Horizontal );
      WTableCell *cell = decayFromTable->elementAt( 0, 0 );
      new WText( "Parent Isotope", cell );
      cell = decayFromTable->elementAt( 0, 1 );
      new WText( "Decay Mode", cell );
      cell = decayFromTable->elementAt( 0, 2 );
      new WText( "Branch Ratio", cell );
      
      string decaysFromHtml;
      for( const SandiaDecay::Transition *parentTrans : nuc->decaysFromParents )
      {
        const SandiaDecay::Nuclide *parentNuclide = parentTrans->parent;
        if( !parentNuclide || (nuc == parentNuclide) )
          continue;
        
        for( const SandiaDecay::Transition *trans : parentNuclide->decaysToChildren )
        {
          if( trans->child != nuc )
            continue;
          
          const int row = decayFromTable->rowCount();
          
          const float br_to = trans->branchRatio;
          
          cell = decayFromTable->elementAt( row, 0 );
          
          WAnchor *anchor = new WAnchor( cell );
          anchor->setText( parentNuclide->symbol );
          anchor->clicked().connect( boost::bind( &MoreNuclideInfoDisplay::setNuclide, this, parentNuclide, history ) );
          
          string txt = "(t<sub>1/2</sub>=" + PhysicalUnits::printToBestTimeUnits( parentNuclide->halfLife, 2 ) + ")";
          WText *t = new WText( txt, cell );
          t->addStyleClass( "HalfLife" );
          
          cell = decayFromTable->elementAt( row, 1 );
          new WText( trans_mode_txt( trans->mode ), cell );
          
          cell = decayFromTable->elementAt( row, 2 );
          new WText( SpecUtils::printCompact( br_to, 5 ), cell );
        }//for( const SandiaDecay::Transition *trans : parentNuclide->decaysToChildren )
      }//for( const SandiaDecay::Transition *parentTrans : nuc->decaysFromParents )
      
      
      tmplt.bindWidget( "decays-from-table", decayFromTable );
    }//if( !nuc->decaysFromParents.empty() )
    
    const auto more_info_db = MoreNucInfoDb::instance();
    if( more_info_db )
    {
      const NucInfo *const more_info = more_info_db->info( nuc );
      
      if( !more_info )
      {
        // We need to fill these conditions out, incase we are setting a nuclide that doesnt have
        //  them, but the previously displayed nuclide did.
        tmplt.setCondition( "if-has-more-info", false );
        tmplt.setCondition( "if-has-related", false );
        tmplt.bindWidget( "related-nucs", nullptr );
        tmplt.bindString( "more-info", "", Wt::TextFormat::XHTMLText );
      }else
      {
        tmplt.setCondition( "if-has-related", !more_info->m_associated.empty() );
        if( !more_info->m_associated.empty() )
        {
          WContainerWidget *associatedList = new WContainerWidget();
          associatedList->setList( true, false );
          associatedList->addStyleClass( "MoreNucInfoRelated" );
          
          for( size_t i = 0; i < more_info->m_associated.size(); ++i )
          {
            const SandiaDecay::Nuclide *associated = db->nuclide( more_info->m_associated[i] );
            WContainerWidget *li = new WContainerWidget( associatedList );
            if( associated )
            {
              WAnchor *anchor = new WAnchor( li );
              anchor->setText( more_info->m_associated[i] );
              anchor->clicked().connect( boost::bind( &MoreNuclideInfoDisplay::setNuclide, this, associated, history ) );
              
              string txt = " (t<sub>1/2</sub>=" + PhysicalUnits::printToBestTimeUnits( associated->halfLife, 2 ) + ")";
              WText *t = new WText( txt, li );
              t->addStyleClass( "HalfLife" );
            }else
            {
              new WText( more_info->m_associated[i], li );
            }
          }
          
          tmplt.bindWidget( "related-nucs", associatedList );
        }//if( !more_info->m_associated.empty() )
        
        
        
        string notes = notes_html( more_info_db, more_info );
        
        tmplt.setCondition( "if-has-more-info", !notes.empty() );
        tmplt.bindString( "more-info", notes, Wt::TextFormat::XHTMLText );
      }//if( more_info )
    }//if( more_info_db )
    
    //if( nuc->canObtainSecularEquilibrium() )
    //  information.push_back( "Can reach secular equilibrium" );
    //else
    //  information.push_back( "Cannot reach secular equilibrium" );
    
    // TODO: maybe compute gamma-dose information?
    
    tmplt.setCondition( "if-valid", true );
    tmplt.setCondition( "if-invalid", false );
    
    m_nuc = nuc;
    m_current_history = history;
    
    // Incase we are loading a nuclide, and the user had scrolled down - we'll scroll up
    doJavaScript( "try{" + jsRef() + ".parentElement.scrollTop=0;}catch{}" );
  }catch( std::exception &e )
  {
    m_nuc = nullptr;
    m_current_history = history;
    
    tmplt.setCondition( "if-valid", false );
    tmplt.setCondition( "if-invalid", true );
    tmplt.bindString( "error-message", e.what(), Wt::TextFormat::XHTMLText );
  }//try / catch
  
  m_nuclideChanged.emit( m_nuc );
  
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    // Pop the end nuclide off, to keep things consistent
    if( !orig_history.empty() && (orig_history.back() == orig_nuc) )
      orig_history.resize( orig_history.size() - 1 );
    
    auto setNuclideLambda = []( const SandiaDecay::Nuclide *newNuc, const vector<const SandiaDecay::Nuclide *> &newHist ){
      InterSpec *interspec = InterSpec::instance();
      ReferencePhotopeakDisplay *display = interspec ? interspec->referenceLinesWidget() : nullptr;
      MoreNuclideInfoWindow *window = display ? display->moreInfoWindow() : nullptr;
      MoreNuclideInfoDisplay *infoDisplay = window ? window->display() : nullptr;
      assert( infoDisplay );
      if( infoDisplay )
        infoDisplay->setNuclide( newNuc, newHist );
    };//setNuclideLambda
    
    
    auto undo = [orig_nuc, orig_history, setNuclideLambda](){
      setNuclideLambda( orig_nuc, orig_history );
    };
    
    auto redo = [nuc, input_history, setNuclideLambda](){
      setNuclideLambda( nuc, input_history );
    };
    
    undoRedo->addUndoRedoStep( undo, redo, "Change more info nuclide." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  
}//void setNuclide(...)


Wt::Signal<const SandiaDecay::Nuclide *> &MoreNuclideInfoDisplay::nuclideChanged()
{
  return m_nuclideChanged;
}

const SandiaDecay::Nuclide *MoreNuclideInfoDisplay::nuclide() const
{
  return m_nuc;
}//const SandiaDecay::Nuclide *nuclide() const;


MoreNuclideInfoWindow::MoreNuclideInfoWindow( const SandiaDecay::Nuclide *const nuc )
  : SimpleDialog( "More info dialog" ), //We need some title text so `SimpleDialog::m_title` is created
  m_display( nullptr ),
  m_orig_nuc( nuc )
{
  m_display = new MoreNuclideInfoDisplay( nuc, false, contents() );
  
  m_display->nuclideChanged().connect( boost::bind( &MoreNuclideInfoWindow::nuclideUpdated, this, boost::placeholders::_1 ) );
  
  nuclideUpdated( nuc );
  
  addButton( "Close" );
}//MoreNuclideInfoWindow


void MoreNuclideInfoWindow::nuclideUpdated( const SandiaDecay::Nuclide *nuc )
{
  if( m_title )
    m_title->setText( "More Info on " + (nuc ? nuc->symbol : string("---")) );
}//void nuclideUpdated( const SandiaDecay::Nuclide *nuc )


const SandiaDecay::Nuclide *MoreNuclideInfoWindow::currentNuclide() const
{
  return m_display ? m_display->nuclide() : nullptr;
}


const SandiaDecay::Nuclide *MoreNuclideInfoWindow::originalNuclide() const
{
  return m_orig_nuc;
}


MoreNuclideInfoDisplay *MoreNuclideInfoWindow::display()
{
  return m_display;
}//
