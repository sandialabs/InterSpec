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
#include <sstream>
#include <iostream>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WAnchor>
#include <Wt/WApplication>
#include <Wt/WPushButton>
#include <Wt/WMemoryResource>
#include <Wt/WContainerWidget>
#include <Wt/WJavaScriptPreamble>

#include "InterSpec/InterSpec.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/LeafletRadMap.h"

using namespace Wt;
using namespace std;


SimpleDialog *LeafletRadMap::showForMeasurement( const std::shared_ptr<const SpecMeas> meas,
                                  std::function<void(LeafletRadMap *)> on_create )
{
  InterSpec *viewer = InterSpec::instance();
  
  const bool showWarning = InterSpecUser::preferenceValue<bool>( "ShowMapDataWarning", viewer );
  
  return nullptr;
}//static LeafletRadMap::showForMeasurement
  

SimpleDialog *LeafletRadMap::showForCoordinate( double latitude, double longitude,
                                 std::function<void(LeafletRadMap *)> on_create )
{
  return nullptr;
}//static SimpleDialog *showForCoordinate


LeafletRadMap::LeafletRadMap( Wt::WContainerWidget *parent )
 : Wt::WContainerWidget( parent ),
  m_meas(),
  m_jsmap( jsRef() + ".map" ),
  m_loadSelected( this )
{
}//LeafletRadMap


LeafletRadMap::~LeafletRadMap()
{
}//~LeafletRadMap()
  

void LeafletRadMap::displayMeasurementOnMap( const std::shared_ptr<const SpecMeas> meas )
{
  // blah blah blah
}//void displayMeasurementOnMap(...)


void LeafletRadMap::displayCoordinate( double latitude, double longitude )
{
  // blah blah blah
}//void displayCoordinate( double latitude, double longitude );



