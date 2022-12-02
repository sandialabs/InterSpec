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
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>


#include <Wt/WText>
#include <Wt/WString>
#include <Wt/WTemplate>
#include <Wt/WPushButton>
#include <Wt/WApplication>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayChainChart.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MoreNuclideInfoDisplay.h"


using namespace std;
using namespace Wt;


MoreNuclideInfoDisplay::MoreNuclideInfoDisplay( const SandiaDecay::Nuclide *const nuc, Wt::WContainerWidget *parent )
  : WTemplate( parent ),
    m_nuc( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/MoreNuclideInfoDisplay.css" );
    
  setNuclide( nuc );
}//MoreNuclideInfoDisplay constructor


void MoreNuclideInfoDisplay::showDecayChainChart()
{
  if( m_nuc )
    DecayChainChart::show_decay_chart_window( m_nuc, DecayChainChart::DecayChainType::DecayFrom );
}//void showDecayChainChart()
  

void MoreNuclideInfoDisplay::showDecayThroughChart()
{
  if( m_nuc )
    DecayChainChart::show_decay_chart_window( m_nuc, DecayChainChart::DecayChainType::DecayThrough );
}//void showDecayThroughChart()


  void MoreNuclideInfoDisplay::setNuclide( const SandiaDecay::Nuclide *const nuc )
  {
    using namespace MoreNuclideInfo;

    WTemplate &tmplt = *this;

    try
    {
      if( !nuc )
        throw std::runtime_error( "Invalid Nuclide" );

      const SandiaDecay::SandiaDecayDataBase *const db = DecayDataBaseServer::database();
      assert( db );
      if( !db )
        throw runtime_error( "Error getting DecayDataBaseServer" );


      const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );

      {// begin reading in more_nuc_info_tmplt.xml
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
        tmplt.setTemplateText( tmplt_txt, Wt::TextFormat::XHTMLText );
      }// end reading in more_nuc_info_tmplt.xml

      const SandiaDecay::Element *el = db->element( nuc->atomicNumber );
      assert( el );
      string elname = el ? el->name : string();
      if( !elname.empty() )
        elname[0] = static_cast<char>(std::toupper( static_cast<unsigned char>(elname[0]) ));

      tmplt.bindString( "symbol", nuc->symbol, Wt::TextFormat::XHTMLText );
      tmplt.bindString( "element", elname, Wt::TextFormat::XHTMLText );
      tmplt.bindInt( "mass-number", nuc->massNumber );

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

      const string amu = PhysicalUnits::printCompact( nuc->atomicMass, 6 ) + " amu";
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
      const string natabun = PhysicalUnits::printCompact( natural_abundance, 4 );
      tmplt.bindString( "natural-abundance", natabun, Wt::TextFormat::XHTMLText );

      WPushButton *decayChainBtn = new WPushButton( "Show Decay Chain" );
      decayChainBtn->addStyleClass( "LinkBtn" );
      decayChainBtn->clicked().connect( this, &MoreNuclideInfoDisplay::showDecayChainChart );
      tmplt.bindWidget( "decay-chain-btn", decayChainBtn );

      tmplt.setCondition( "if-has-parents", !nuc->decaysFromParents.empty() );
      if( !nuc->decaysFromParents.empty() )
      {
        WPushButton *decayThroughBtn = new WPushButton( "Show Decays Through" );
        decayThroughBtn->addStyleClass( "LinkBtn" );
        decayThroughBtn->clicked().connect( this, &MoreNuclideInfoDisplay::showDecayThroughChart );
        tmplt.bindWidget( "decay-through-btn", decayThroughBtn );
      }

      // TODO: WTemplate has a Wt::WTemplate::Functions::while_f function, but from the documentation it isnt totally clear how to use it; it would be nice to use this to generate rows in a table or something to list all the individual decays
      string decaysToHtml;
      for( const SandiaDecay::Transition *transition : nuc->decaysToChildren )
      {
        if( transition->child )
        {
          decaysToHtml +=
            "<div>Decays to " + transition->child->symbol
            + " by "
            + SandiaDecay::to_str( transition->mode )
            + " decay, BR=" + PhysicalUnits::printCompact( transition->branchRatio, 5 )
            + "</div>";
        }//if( transition->child )
      }//for( const SandiaDecay::Transition * transition : nuc->decaysToChildren)

      tmplt.bindString( "decays-to-html", decaysToHtml, Wt::TextFormat::XHTMLText );


      string decaysFromHtml;
      for( const SandiaDecay::Transition *parentTrans : nuc->decaysFromParents )
      {
        const SandiaDecay::Nuclide *parentNuclide = parentTrans->parent;

        if( parentNuclide && (nuc != parentNuclide) )
        {
          const float br_from = parentNuclide->branchRatioFromForebear( nuc );
          const float br_to = parentNuclide->branchRatioToDecendant( nuc );

          if( br_from <= 0.0 || br_to > 0.0 )
            decaysFromHtml += "<div>Branch Ratio from " + parentNuclide->symbol + ": " + PhysicalUnits::printCompact( br_to, 5 ) + "</div>";
          else
            decaysFromHtml += "<div>Branch Ratio through " + parentNuclide->symbol + ": " + PhysicalUnits::printCompact( br_from, 5 ) + "</div>";
        }
      }//for( const SandiaDecay::Transition *parentTrans : nuc->decaysFromParents )

      tmplt.bindString( "decays-from-html", decaysFromHtml, Wt::TextFormat::XHTMLText );

      const auto more_info_db = MoreNucInfoDb::instance();
      if( more_info_db )
      {
        const NucInfo *const more_info = more_info_db->info( nuc );

        if( more_info )
        {
          string related;
          for( size_t i = 0; i < more_info->m_associated.size(); ++i )
          {
            related += (i ? ", " : "") + more_info->m_associated[i];
            const SandiaDecay::Nuclide *associated = db->nuclide( more_info->m_associated[i] );
            if( associated )
              related += " (hl=" + PhysicalUnits::printToBestTimeUnits( associated->halfLife, 2 ) + ")";
          }
          tmplt.setCondition( "if-has-related", !related.empty() );
          tmplt.bindString( "related-nucs", related, Wt::TextFormat::PlainText );

          string notes = more_info->m_notes;
          SpecUtils::trim( notes ); //just to make sure (although I think the XML parser already did this)

          tmplt.setCondition( "if-has-more-info", !notes.empty() );
          SpecUtils::ireplace_all( notes, "\r", "" );
          vector<string> more_info_lines;

          boost::algorithm::split( more_info_lines, notes,
            boost::algorithm::is_any_of( "\n" ),
            boost::algorithm::token_compress_off );

          notes = "<p>";
          for( size_t i = 0; i < more_info_lines.size(); ++i )
          {
            string line = more_info_lines[i];
            SpecUtils::trim( line );

            notes += line;

            if( line.empty() )
            {
              notes += "</p><p>";

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

          notes += "</p>";
          tmplt.bindString( "more-info", notes, Wt::TextFormat::XHTMLText );
        }//if( more_info )
      }//if( more_info_db )

      //if( nuc->canObtainSecularEquilibrium() )
      //  information.push_back( "Can reach secular equilibrium" );
      //else
      //  information.push_back( "Cannot reach secular equilibrium" );

      tmplt.setCondition( "if-valid", true );
      tmplt.setCondition( "if-invalid", false );

      m_nuc = nuc;
    } catch( std::exception &e )
    {
      m_nuc = nullptr;
      tmplt.setCondition( "if-valid", false );
      tmplt.setCondition( "if-invalid", true );
      tmplt.bindString( "error-message", e.what(), Wt::TextFormat::XHTMLText );
    }//try / catch
  }//void MoreNuclideInfoDisplay::setNuclide(...)



MoreNuclideInfoWindow::MoreNuclideInfoWindow( const SandiaDecay::Nuclide *const nuc )
  : SimpleDialog( "More Info on " + (nuc ? nuc->symbol : string("null")), "" ),
  m_display( nullptr )
{
  m_display = new MoreNuclideInfoDisplay( nuc, contents() );
  addButton( "Close" );
}//MoreNuclideInfoWindow