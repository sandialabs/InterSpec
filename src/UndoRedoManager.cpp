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
#include <iostream>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/UndoRedoManager.h"


using namespace std;


UndoRedoManager::UndoRedoManager( InterSpec *parent )
 : Wt::WObject( parent ),
  m_interspec( parent )
{
  m_interspec->displayedSpectrumChanged().connect(
                boost::bind( &UndoRedoManager::handleSpectrumChange, this,
                             boost::placeholders::_1, boost::placeholders::_2,
                            boost::placeholders::_3, boost::placeholders::_4 ) );
}//UndoRedoManager


UndoRedoManager::~UndoRedoManager()
{
  
}//~UndoRedoManager()


void UndoRedoManager::handleSpectrumChange( const SpecUtils::SpectrumType type,
                                           const shared_ptr<SpecMeas> &meas,
                                           const set<int> &sample_nums,
                                           const vector<std::string> &detector_names )
{
  
}//void handleSpectrumChange( SpecUtils::SpectrumType type )


