#ifndef Daq_h
#define Daq_h
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

class Measurement;

namespace Daq
{
  /** Function to call to init connection to DAQ.
      Currently just parses example_spectra/passthrough.n42 then starts a thread
      that every few secons adds a sample number from passthrough.n42 to a
      spectrum being summed, and then calls #update_displayed_spectrum.
   
      Does not currently use argc/argv, but if you wanted to parse them to
      figure out which port to connect to for ZeroMQ or whatever, you could do
      this.
   */
  void init_daq( int argc, char **argv );
  
  
  /** Sets the displayed spectra for all connected clients to the ones passed in.
   Note: if you expect to have sources that change over time (e.g., rapid decay)
     you could instead change things to pass in std::shared_ptr<MeasurementInfo>
     objects that have many Measurements and that way get a time history
     displayed as well, much like when you show the example "passthrough.n42"
     file normally in InterSpec
   */
  void update_displayed_spectrum( std::shared_ptr<const Measurement> foreground,
                                  std::shared_ptr<const Measurement> background );
  
  /** Stops the daq thread, or ZeroMQ socket, or whatever. */
  void stop_daq();
  
}//namespace Daq

#endif //Daq_h
