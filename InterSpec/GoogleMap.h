#ifndef GoogleMap_h
#define GoogleMap_h
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

#include <set>
#include <vector>
#include <memory>

#include <Wt/WString>
#include <Wt/WContainerWidget>
#include <Wt/WResource>
#include <Wt/WGoogleMap>

#include "InterSpec/AuxWindow.h"

class SrbGoogleMap;
namespace SpecUtils { class SpecFile; }

class GoogleMap : public Wt::WContainerWidget
{
public:
  GoogleMap( const bool trackMapExtent, Wt::WContainerWidget *parent = 0 );
  virtual ~GoogleMap();
  
  void addMarker( double latitude, double longitude );
  void addInfoBox( double latitude, double longitude, const Wt::WString &html );
  void addMeasurment( std::shared_ptr<const SpecUtils::SpecFile> meas,
                      const Wt::WString &title,
                      const std::set<int> &displayed );
  std::string getStaticMapFromMeas( std::shared_ptr<const SpecUtils::SpecFile> meas,
                             const Wt::WString &title );
  std::string getStaticMap();
  void clearMeasurments();
  void adjustPanAndZoom();
  
 
  //getMapExtent(...): If 'trackMapExtent' wasnt true when this GoogleMap was
  //  created, then an exception will be thrown
  void getMapExtent( float &upperLatitude, float &leftLongitude,
                     float &lowerLatitude, float &rightLongitude );
  Wt::Signal<double /*Lat*/, double /*Lng*/> &mapClicked();
  
protected:
  void init();
  void updateMapGpsCoords( const double , const double ,
                           const double , const double );
  
  SrbGoogleMap *m_map;
  
  const bool m_trackMapExtent;
  float m_mapExtents[4];
  
  Wt::Signal<double,double> m_clicked;
};//class GoogleMap

#endif
