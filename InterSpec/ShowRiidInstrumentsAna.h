#ifndef SHOW_RIID_INSTRUMENTS_ANA_h
#define SHOW_RIID_INSTRUMENTS_ANA_h
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
#include <memory>

class SpecMeas;
class SimpleDialog;

/** Gives a short HTML string summarizing results of the RIID analysis results in 'spec'.
 This is used to place RIID summary information into the application notifications when you open a file.
 
 Will return empty string if no (useful) information available.
 */
std::string riidAnaSummary( const std::shared_ptr<const SpecMeas> &spec );


/** Creates a popup containing the RIID analysis information.
 
 @returns The dialog that is created.
 */
SimpleDialog *showRiidInstrumentsAna( const std::shared_ptr<const SpecMeas> &spec );

#endif //SHOW_RIID_INSTRUMENTS_ANA_h


