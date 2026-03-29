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

#include <Wt/Utils.h>
#include <Wt/WText.h>
#include <Wt/WGridLayout.h>
#include <Wt/WApplication.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActTxtResults.h"

using namespace Wt;
using namespace std;


RelActTxtResults::RelActTxtResults()
 : Wt::WContainerWidget(),
  m_txt( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/RelActTxtResults.css" );

  addStyleClass( "RelActTxtResults" );

  // Without using a layout the formatting messes up and the contents move over
  //  onto the tabs... not sure why, but using a layout fixes this.
  WGridLayout *lay = setLayout( std::make_unique<WGridLayout>() );
  m_txt = lay->addWidget( std::make_unique<WContainerWidget>(), 0, 0 );
  lay->setRowStretch( 0, 1 );
  lay->setColumnStretch( 0, 1 );

  m_txt->addNew<WText>( "No current results" );
}


void RelActTxtResults::setNoResults()
{
  m_txt->clear();
  m_txt->addNew<WText>( "No current results" );
}


void RelActTxtResults::updateResults( const RelActCalcAuto::RelActAutoSolution &solution )
{
  m_txt->clear();
  
  try
  {
    stringstream strm;
    solution.print_summary( strm );
    
    string res = strm.str();
    res = Wt::Utils::htmlEncode( res, Wt::Utils::HtmlEncodingFlag::EncodeNewLines );
    WText *txt = m_txt->addNew<WText>( res );
    txt->setInline( false );
  }catch( std::exception &e )
  {
    addNew<WText>( "Error displaying results" );
  }
}//void updateResults( const RelActCalcAuto::RelActAutoSolution &solution );
