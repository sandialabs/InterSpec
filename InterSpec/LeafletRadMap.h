#ifndef LeafletRadMap_h
#define LeafletRadMap_h
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
#include <string>
#include <vector>
#include <functional>

#include <Wt/WContainerWidget>

// Forward declarations
class SpecMeas;
class InterSpec;
class SimpleDialog;
namespace SpecUtils{ enum class SpectrumType : int; }

class LeafletRadMap : public Wt::WContainerWidget
{
public:
  static SimpleDialog *showForMeasurement( const std::shared_ptr<const SpecMeas> meas,
                                  std::function<void(LeafletRadMap *)> on_create = nullptr );
  
  static SimpleDialog *showForCoordinate( double latitude, double longitude,
                                 std::function<void(LeafletRadMap *)> on_create = nullptr );

  LeafletRadMap( Wt::WContainerWidget *parent = nullptr );
  virtual ~LeafletRadMap();
  
  void displayMeasurementOnMap( const std::shared_ptr<const SpecMeas> meas );
  void displayCoordinate( double latitude, double longitude );
  
protected:
  std::shared_ptr<const SpecMeas> m_meas;
  
  
  /** The javascript variable name used to refer to the LeafletRadMap JS object.
      Currently is `jsRef() + ".map"`.
   */
  const std::string m_jsmap;
  
  
  Wt::Signal<SpecUtils::SpectrumType, std::shared_ptr<const SpecMeas>, std::set<int>> m_loadSelected;
};//class LeafletRadMap

#endif //LeafletRadMap_h
