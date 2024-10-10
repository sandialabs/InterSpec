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

#include <assert.h>

#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/CoincidenceDev.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace CoincidenceDev
{
  int run_coinc_dev_code()
  {
    // Example of printing out coincidence for a specific nuclide
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      throw runtime_error( "Could initialize decay database" );
    
    const SandiaDecay::Nuclide *nuc = db->nuclide( "U235" );
    assert( nuc );
    if( !nuc )
      throw runtime_error( "Uninitialized decay database" ); //shouldnt ever happen
    

    const vector<const SandiaDecay::Transition *> &decay_channels = nuc->decaysToChildren;
    for( const SandiaDecay::Transition *trans : decay_channels )
    {
      const vector<SandiaDecay::RadParticle> &products = trans->products;
      for( const SandiaDecay::RadParticle &particle : products )
      {
        if( particle.type != SandiaDecay::GammaParticle )
          continue;
        
        const vector<pair<unsigned short int,float>> &coincidences = particle.coincidences;
        for( const pair<unsigned short int,float> &coinc : coincidences )
        {
          const unsigned short int part_ind = coinc.first;
          const float fraction = coinc.second;
          assert( part_ind < products.size() );
          if( part_ind >= products.size() )
            throw runtime_error( "Invalid coincidence index" );
          
          const SandiaDecay::RadParticle &coinc_part = products[part_ind];
            
          if( coinc_part.type != SandiaDecay::ProductType::GammaParticle )
            continue;
          
          cout << "The " << particle.energy << " keV gamma is coincident with the "
                << coinc_part.energy << " keV gamma, " << fraction << " part of the time." << endl;
        }//for( loop over coincidences )
      }//for( loop over RadParticles )
    }//for( loop over transitions )
    
    DetectorPeakResponse det_eff;
  
    try
    {
      // To load a GADRAS efficiency function from a directory, the directory needs a Detector.dat,
      //  and Efficiency.csv file.
      det_eff.fromGadrasDirectory( "data/GenericGadrasDetectors/HPGe 40%/" );
    }catch( std::exception &e )
    {
      cerr << "Failed to open GADRAS detector: " << e.what() << endl;
      return 1;
    }
    
    // TODO: when loading a detector  peak response from a GADRAS source, we should also keep track of the total scatter column, and maybe other columns, so we can use that to help us compute coincidence effects
    
    vector<PeakDef> peaks;
    try
    {
      string peak_csv_filename = "some_peaks.CSV";
#ifdef _WIN32
      const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(peak_csv_filename);
      std::ifstream inpeaks( wfilename.c_str() );
#else
      std::ifstream inpeaks( peak_csv_filename.c_str() );
#endif
      
      if( !inpeaks )
        throw runtime_error( "Failed to open " + peak_csv_filename );
      shared_ptr<const SpecUtils::Measurement> meas; //We could load a spectrum file and get this, but we dont need to
      peaks = PeakModel::csv_to_candidate_fit_peaks( meas, inpeaks );
    }catch( std::exception &e )
    {
      cerr << "Failed to load peaks from CSV file: " << e.what() << endl;
    }//try / catch to load peaks from a CSV file
    
    return 0;
  }//int run_coinc_dev_code()
}//namespace CoincidenceDev
