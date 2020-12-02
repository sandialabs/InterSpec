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

#include <deque>
#include <string>
#include <vector>
#include <memory>
#include <cctype>
#include <algorithm>

#include <Wt/WResource>
#include <Wt/WModelIndex>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStringStream>
#include <Wt/WAbstractItemModel>


// Block out some warnings occurring in boost files.
#pragma warning(disable:4244)  // warning C4244: 'initializing' : conversion from 'std::streamoff' to 'size_t', possible loss of data

#include <boost/any.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayDataBaseServer.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/DateTime.h"
#include "InterSpec/SpecMeas.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PeakModel.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/PeakInfoDisplay.h"  //Only for ALLOW_PEAK_COLOR_DELEGATE
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/SpectrumDataModel.h"


using namespace std;
using namespace Wt;


bool PeakModel::recomendUseForFit( const SandiaDecay::Nuclide *nuc,
                                   const float energy )
{
  //These are the lines that should be recomended to users to use in a
  //  shielding source fit
  static const float Ag110m[] = {657.762f, 706.682f, 763.944f, 884.685f, 937.493f, 1384.29f, 1505.04f};
  static const float Al26[]   = {1129.67f, 1808.65f, 2938.0f};
  static const float Am241[]  = {26.35f, 59.54f, 125.3f, 208.01f, 335.37f, 376.65f, 662.4f, 722.01f};
  static const float Am243[]  = {74.66f, 103.032f, 106.123f, 209.753f, 228.183f, 277.599f, 315.88f, 334.31f};
  static const float Au198[]  = {411.802f, 675.884f, 1087.68f};
  static const float Ba133[]  = {80.9971f, 276.4f, 302.853f, 356.017f, 383.848f};
  static const float Ba140[]  = {487.021f, 537.3f, 815.8f, 925.2f, 1596.21f, 2521.7f};
  static const float Be7[]    = {477.595f};
  static const float Bi207[]  = {75.0f, 569.702f, 1063.66f, 1442.2f, 1770.24f};
  static const float Cd109[]  = {88.04f};
  static const float Ce139[]  = {165.85f};
  static const float Ce141[]  = {145.443f};
  static const float Cf249[]  = {109.2f, 252.8f, 333.37f, 388.16f};
  static const float Cm243[]  = {103.734f, 106.125f, 209.753f, 228.183f, 277.599f, 315.88f, 334.31f};
  static const float Cm244[]  = {98.86f, 104.23f, 152.63f, 554.5f, 817.8f};
  static const float Cm245[]  = {103.68f, 133.05f, 175.01f, 231.96f};
  static const float Co56[]   = {846.771f, 1238.28f, 1771.35f, 2598.46f, 3253.42f};
  static const float Co57[]   = {122.061f, 136.474f, 692.03f};
  static const float Co58[]   = {810.759f};
  static const float Co60[]   = {1173.23f, 1332.49f};
  static const float Cr51[]   = {320.082f};
  static const float Cs134[]  = {569.331f, 604.721f, 795.8f, 1038.5f, 1167.86f, 1365.13f};
  static const float Cs136[]  = {340.547f, 818.514f, 1048.07f, 1235.36f, 1537.9f};
  static const float Cs137[]  = {661.657f};
  static const float Cs138[]  = {462.796f, 871.8f, 1009.78f, 1435.86f, 2218.0f, 2639.59f};
  static const float Cu61[]   = {282.956f, 656.008f, 1185.23f};
  static const float Cu67[]   = {184.577f, 300.219f};
  static const float Eu152[]  = {121.781f, 244.697f, 344.278f, 411.116f, 443.965f, 778.904f, 867.378f, 964.079f, 1085.86f, 1112.07f, 1408.01f};
  static const float Eu154[]  = {123.14f, 723.305f, 873.19f, 1004.76f, 1274.39f, 1596.48f};
  static const float Eu155[]  = {86.54f, 105.308f};
  static const float Fe55[]   = {126.0f};
  static const float Fe59[]   = {192.35f, 334.8f, 1099.24f, 1291.59f, 1481.7f};
  static const float Ga67[]   = {93.311f, 184.577f, 300.219f, 393.529f, 887.693f};
  static const float Ge68[]   = {805.75f, 1077.35f, 1260.97f, 1883.09f};
  static const float Ge75[]   = {264.6f, 468.8f};
  static const float Hg194[]  = {293.58f, 328.5f, 948.29f, 1468.89f, 2043.67f};
  static const float Hg203[]  = {279.195f};
  static const float Ho166m[] = {184.41f, 280.459f, 410.944f, 711.683f, 810.276f, 1241.4f};
  static const float I123[]   = {158.97f, 528.96f, 624.57f, 1068.12f};
  static const float I125[]   = {27.472f};
  static const float I126[]   = {388.633f, 666.331f, 753.819f, 1420.17f, 2045.17f};
  static const float I131[]   = {284.305f, 364.489f, 503.004f, 636.989f, 722.911f};
  static const float In111[]  = {171.28f, 245.395f};
  static const float Ir192[]  = {295.957f, 316.506f, 468.069f, 604.411f, 884.537f, 1061.48f};
  static const float K40[]    = {1460.75f};
  static const float La140[]  = {487.0f, 815.8f, 925.2f, 1596.21f, 2521.7f};
  static const float Lu172[]  = {181.525f, 810.064f, 900.724f, 1093.63f, 1584.12f};
  static const float Lu177[]  = {112.94f, 208.36f, 321.31f};
  static const float Lu177m[] = {55.15f, 208.366f, 228.484f, 378.503f, 418.539f};
  static const float Mn54[]   = {834.848f};
  static const float Mn56[]   = {846.754f};
  static const float Mo99[]   = {140.511f, 181.063f, 366.421f, 739.5f, 777.921f};
  static const float Na22[]   = {1274.54f};
  static const float Na24[]   = {996.82f, 1368.63f, 2754.03f, 3866.19f};
  static const float Nb95[]   = {765.803f};
  static const float Np237[]  = {86.477f, 194.95f, 212.29f, 300.129f, 311.904f, 340.476f, 375.45f, 398.492f, 415.76f};
  static const float Np239[]  = {103.032f, 106.123f, 209.753f, 228.183f, 277.599f, 315.88f, 334.31f};
  static const float Pb202[]  = {439.56f};
  static const float Pd103[]  = {294.98f, 357.47f, 497.08f};
  static const float Po210[]  = {803.13f};
  static const float Pu238[]  = {43.498f, 99.853f, 152.72f, 742.81f, 766.39f, 1001.03f};
  static const float Pu239[]  = {129.297f, 203.55f, 255.38f, 345.01f, 375.05f, 413.71f, 451.48f, 639.97f, 645.896f, 769.37f};
  static const float Pu240[]  = {45.244f, 104.23f, 160.308f, 642.35f, 687.59f};
  static const float Pu241[]  = {59.54f, 98.97f, 103.68f, 148.567f, 208.0f};
  static const float Pu242[]  = {44.915f};
  static const float Ra226[]  = {74.815f, 77.107f, 186.211f, 241.997f, 295.24f, 351.932f, 609.312f, 768.356f, 1120.29f, 1238.11f, 1377.67f, 1407.98f, 1764.49f, 2204.21f, 2447.86f};
  static const float Th232[]  = {77.34f, 238.632f, 338.32f, 583.191f, 727.33f, 860.564f, 911.204f, 968.971f, 1588.2f, 2614.53f};
  static const float U233[]   = {146.345f, 164.522f, 291.354f, 317.16f, 1567.09f};
  static const float U235[]   = {106.587f, 143.76f, 163.38f, 185.71f, 205.309f};
  static const float U238[]   = {63.24f, 92.35f, 258.26f, 742.77f, 766.37f, 1000.99f, 1737.73f, 1831.36f};
  static const std::string symbols[]
         = { "Ag110m", "Al26", "Am241", "Am243", "Au198", "Ba133", "Ba140",
             "Be7", "Bi207", "Cd109", "Ce139", "Ce141", "Cf249", "Cm243",
             "Cm244", "Cm245", "Co56", "Co57", "Co58", "Co60", "Cr51", "Cs134",
             "Cs136", "Cs137", "Cs138", "Cu61", "Cu67", "Eu152", "Eu154",
             "Eu155", "Fe55", "Fe59", "Ga67", "Ge68", "Ge75", "Hg194", "Hg203",
             "Ho166m", "I123", "I125", "I126", "I131", "In111", "Ir192", "K40",
             "La140", "Lu172", "Lu177", "Lu177m", "Mn54", "Mn56", "Mo99",
             "Na22", "Na24", "Nb95", "Np237", "Np239", "Pb202", "Pd103",
             "Po210", "Pu238", "Pu239", "Pu240", "Pu241", "Pu242", "Ra226",
             "Th232", "U233", "U235", "U238"
           };
  static const size_t num_sumbols = sizeof(symbols) / sizeof(symbols[0]);
  static const size_t num_energies[]
           = { 7, 3, 8, 8, 3, 5, 6, 1, 5, 1, 1, 1, 4, 7, 5, 4, 5, 3, 1, 2, 1,
               6, 5, 1, 6, 3, 2, 11, 6, 2, 1, 5, 5, 4, 2, 5, 1, 6, 4, 1, 5, 5,
               2, 6, 1, 5, 5, 3, 5, 1, 1, 5, 1, 4, 1, 9, 7, 1, 3, 1, 6, 10, 5,
               5, 1, 15, 10, 5, 5, 8
             };
  static const float *isotope_energies[]
         = { &(Ag110m[0]), &(Al26[0]), &(Am241[0]), &(Am243[0]), &(Au198[0]),
             &(Ba133[0]), &(Ba140[0]), &(Be7[0]), &(Bi207[0]), &(Cd109[0]),
             &(Ce139[0]), &(Ce141[0]), &(Cf249[0]), &(Cm243[0]), &(Cm244[0]),
             &(Cm245[0]), &(Co56[0]), &(Co57[0]), &(Co58[0]), &(Co60[0]),
             &(Cr51[0]), &(Cs134[0]), &(Cs136[0]), &(Cs137[0]), &(Cs138[0]),
             &(Cu61[0]), &(Cu67[0]), &(Eu152[0]), &(Eu154[0]), &(Eu155[0]),
             &(Fe55[0]), &(Fe59[0]), &(Ga67[0]), &(Ge68[0]), &(Ge75[0]),
             &(Hg194[0]), &(Hg203[0]), &(Ho166m[0]), &(I123[0]), &(I125[0]),
             &(I126[0]), &(I131[0]), &(In111[0]), &(Ir192[0]), &(K40[0]),
             &(La140[0]), &(Lu172[0]), &(Lu177[0]), &(Lu177m[0]), &(Mn54[0]),
             &(Mn56[0]), &(Mo99[0]), &(Na22[0]), &(Na24[0]), &(Nb95[0]),
             &(Np237[0]), &(Np239[0]), &(Pb202[0]), &(Pd103[0]), &(Po210[0]),
             &(Pu238[0]), &(Pu239[0]), &(Pu240[0]), &(Pu241[0]), &(Pu242[0]),
             &(Ra226[0]), &(Th232[0]), &(U233[0]), &(U235[0]), &(U238[0])
           };

  if( !nuc )
    return false;
  
  const string *pos = lower_bound( symbols, symbols+num_sumbols, nuc->symbol );
  if( pos && (pos!=(symbols+num_sumbols)) && (*pos)==nuc->symbol )
  {
    const float *energies = isotope_energies[pos-symbols];
    const size_t nenergies = num_energies[pos-symbols];
    const float *epos = lower_bound( energies, energies+nenergies, energy );
    
    //*epos will be greater or equal to energy
    if( epos == (energies+nenergies) )
      return (fabs((*(epos-1))-energy)<0.1);  //assumes at least one element in array
    
    return (((*epos)==energy) || (fabs((*epos)-energy)<0.1)
             || ((epos!=energies) && (fabs((*(epos-1))-energy)<0.1)));
  }//if( we have info on this elemement )
  
  //todo - put in logic here as how to guess if a photopeak should be used or not
  
//  if( nuc && nuc->)
//  if( fabs(511.0-energy) < 1.0 )
//    return false;
  
  
  
  return true;
}//recomendUseForFit( const SandiaDecay::Nuclide *nuc )


std::vector<PeakDef> PeakModel::csv_to_candidate_fit_peaks(
                                      std::shared_ptr<const SpecUtils::Measurement> meas,
                                      std::istream &csv )
{
  //Info that will be parsed is based on PeakModel::PeakCsvResource::handleRequest(...)
  //  with a few small accommodations for how other programs may save peak CSVs.
  
  using SpecUtils::trim_copy;
  using SpecUtils::to_lower_ascii_copy;
  
  typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
  boost::escaped_list_separator<char> separator("\\",",\t", "\"");
  
  if( !meas || !meas->gamma_counts() || meas->gamma_counts()->size() < 7 )
    throw runtime_error( "input data invalid" );
  
  const float minenergy = meas->gamma_energy_min();
  const float maxenergy = meas->gamma_energy_max();
  
  
  string line;
  
  //Get first non-empty, non-comment ('#') line.
  while( SpecUtils::safe_get_line(csv, line, 2048) )
  {
    SpecUtils::trim(line);
    if( !line.empty() && line[0]!='#' )
      break;
  }
  
  if( line.empty() || !csv )
    throw runtime_error( "Failed to get first line" );
  
  //Columns garunteed to be in file, or we'll throw an exception.
  int mean_index = -1, area_index = -1, fwhm_index = -1;
  
  //Columns that may or not be in file, in whcih case will be >= 0.
  int roi_lower_index = -1, roi_upper_index = -1, nuc_index = -1, nuc_energy_index = -1;
  
  {//begin to get field_pos
    vector<string> headers;
    Tokeniser t( line, separator );
    for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
      headers.push_back( to_lower_ascii_copy( trim_copy(*it) ) );
    
    const auto centroid_pos = std::find( begin(headers), end(headers), "centroid");
    const auto net_area_pos = std::find( begin(headers), end(headers), "net_area");
    //Note First FWHM is in keV; there is a second FWHM that is %
    const auto fwhm_pos = std::find( begin(headers), end(headers), "fwhm");
    const auto nuc_pos = std::find( begin(headers), end(headers), "nuclide");
    const auto nuc_energy_pos = std::find( begin(headers), end(headers), "photopeak_energy");
    const auto roi_start_pos = std::find( begin(headers), end(headers), "roi_lower_energy");
    const auto roi_end_pos = std::find( begin(headers), end(headers), "roi_upper_energy");
    
    if( centroid_pos == end(headers) )
      throw runtime_error( "Header did not contain 'Centroid'" );
    mean_index = static_cast<int>( centroid_pos - begin(headers) );
    
    if( net_area_pos == end(headers) )
      throw runtime_error( "Header did not contain 'Net_Area'" );
    area_index = static_cast<int>( net_area_pos - begin(headers) );
    
    if( fwhm_pos == end(headers) )
      throw runtime_error( "Header did not contain 'FWHM'" );
    fwhm_index = static_cast<int>( fwhm_pos - begin(headers) );
    
    if( nuc_pos != end(headers) )
      nuc_index = static_cast<int>( nuc_pos - begin(headers) );
    
    if( nuc_energy_pos != end(headers) )
      nuc_energy_index = static_cast<int>( nuc_energy_pos - begin(headers) );
    
    if( roi_start_pos != end(headers) )
      roi_lower_index = static_cast<int>( roi_start_pos - begin(headers) );
    
    if( roi_end_pos != end(headers) )
      roi_upper_index = static_cast<int>( roi_end_pos - begin(headers) );
  }//end to get field_pos
  
  
  vector<PeakDef> answer;
  
  while( SpecUtils::safe_get_line(csv, line, 2048) )
  {
    SpecUtils::trim(line);
    if( line.empty() || line[0]=='#' || (!isdigit(line[0]) && line[0]!='+' && line[0]!='-') )
      continue;
    
    vector<string> fields;
    Tokeniser t( line, separator );
    for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
      fields.push_back( to_lower_ascii_copy( trim_copy(*it) ) );
    
    const int nfields = static_cast<int>( fields.size() );
    
    if( nfields < mean_index || nfields < area_index || nfields < fwhm_index )
      continue;
    
    try
    {
      const float centroid = std::stof( fields[mean_index] );
      const float fwhm = std::stof( fields[fwhm_index] );
      const float area = std::stof( fields[area_index] );
      
      if( centroid <= minenergy || centroid >= maxenergy || fwhm <= 0.0 || area <= 0.0 )
        continue;
      
      PeakDef peak( centroid, fwhm/2.35482, area );
      
      if( roi_lower_index >= 0 && roi_lower_index < nfields
         && roi_upper_index >= 0 && roi_upper_index < nfields )
      {
        const float roi_lower = std::max( minenergy, std::stof( fields[roi_lower_index] ) );
        const float roi_upper = std::min( maxenergy, std::stof( fields[roi_upper_index] ) );
        if( roi_lower >= roi_upper || centroid < roi_lower || centroid > roi_upper )
          throw runtime_error( "ROI range invalid." );
        
        peak.continuum()->setRange( roi_lower, roi_upper );
        peak.continuum()->calc_linear_continuum_eqn( meas, roi_lower, roi_upper, 1 );
      }else
      {
        double lowerEnengy, upperEnergy;
        findROIEnergyLimits( lowerEnengy, upperEnergy, peak, meas );
        
        peak.continuum()->setRange( lowerEnengy, upperEnergy );
        peak.continuum()->calc_linear_continuum_eqn( meas, lowerEnengy, upperEnergy, 1 );
      }//if( CSV five ROI extent ) / else( find from data )
      
      
      if( nuc_index >= 0 && nuc_energy_index >= 0
         && nuc_index < nfields && nuc_energy_index < nfields
         && !fields[nuc_index].empty() && !fields[nuc_energy_index].empty() )
      {
        const string nuctxt = fields[nuc_index] + " " + fields[nuc_energy_index] + " keV";
        const SetGammaSource result = setNuclideXrayReaction( peak, nuctxt, 4.0 );
        if( result == SetGammaSource::NoSourceChange )
          cerr << "csv_to_candidate_fit_peaks: could not assign src txt '"
          << nuctxt << "' as a nuc/xray/rctn" << endl;
      }//if( nuc_index >= 0 || nuc_energy_index >= 0 )
      
      //Go through existing peaks and if the new peak should share a ROI, do that here
      if( peak.continuum()->energyRangeDefined() )
      {
        for( auto &p : answer )
        {
          if( !p.continuum()->energyRangeDefined() )
            continue;
          if( fabs(p.continuum()->lowerEnergy() - peak.continuum()->lowerEnergy()) < 0.0001
             && fabs(p.continuum()->upperEnergy() - peak.continuum()->upperEnergy()) < 0.0001 )
          {
            peak.setContinuum( p.continuum() );
            //ToDo: determine if we need to make the continuum quadratic or higher...
            //peak.setType( OffsetType::Quadratic );
          }//if( ROI bounds are about same as another continuum )
        }//for( loop over peaks previously found to see if we should share continuum )
      }//if( new peak has a defined energy range )
      
      answer.push_back( peak );
    }catch( std::exception &e )
    {
      throw runtime_error( "Invalid value on line '" + line + "', " + string(e.what()) );
    }//try / catch to parse a line into a peak
  }//while( SpecUtils::safe_get_line(csv, line, 2048) )
  
  
  if( answer.empty() )
    throw runtime_error( "No peak rows found in file." );
  
  return answer;
}//csv_to_candidate_fit_peaks(...)


PeakModel::PeakCsvResource::PeakCsvResource( PeakModel *parent )
  : WResource( parent ),
    m_model( parent )
{
}


PeakModel::PeakCsvResource::~PeakCsvResource()
{
  beingDeleted();
}


void PeakModel::PeakCsvResource::handleRequest( const Wt::Http::Request &/*request*/,
                                                Wt::Http::Response& response )
{
  //If you update the fields or headers of this function, you should also update
  //  PeakModel::csv_to_candidate_fit_peaks(...)
  
  std::shared_ptr<SpecMeas> meas = m_model->m_measurment.lock();
  string filename = "peaks.CSV", specfilename = "Unknown";
  if( meas && !meas->filename().empty() )
  {
    specfilename = meas->filename();
    filename = "peaks_" + specfilename;
    const string::size_type pos = filename.find_last_of( "." );
    if( pos != string::npos )
      filename = filename.substr( 0, pos );
    filename += ".CSV";
  }//if( meas && !meas->filename().empty() )
    
  suggestFileName( filename, WResource::Attachment ); //WResource::NoDisposition

  response.setMimeType( "text/csv" );

  if( !m_model )
    return;

  const size_t npeaks = m_model->npeaks();
  const string eol_char = "\r\n"; //for windows - could potentially cosutomize this for the users operating system
  
  response.out() << "Centroid,  Net_Area,   Net_Area,      Peak, FWHM,   FWHM,Reduced, ROI_Total,ROI, "
                    "File,         ,     ,     , Nuclide, Photopeak_Energy, ROI_Lower_Energy, ROI_Upper_Energy"
                 << eol_char
                 << "     keV,    Counts,Uncertainty,       CPS,  keV,Percent,Chi_Sqr,    Counts,ID#, "
                    "Name, LiveTime, Date, Time,        ,              keV,              keV,              keV"
                 << eol_char;
  
  
  for( size_t peakn = 0; peakn < npeaks; ++peakn )
  {
    const PeakDef &peak = m_model->peak( peakn );

    std::shared_ptr<const SpecUtils::Measurement> data = m_model->m_dataModel->getData();
    float live_time = 0.0;
    boost::posix_time::ptime meastime;
    double region_area = 0.0, xlow = 0.0, xhigh = 0.0;
    if( data )
    {
      live_time = data->live_time();
      meastime = data->start_time();
      findROIEnergyLimits( xlow, xhigh, peak, data );
      region_area = gamma_integral( data, xlow, xhigh );
    }//if( data )

    string nuclide;
    char energy[32];
    energy[0] = '\0';
    try
    {
      const float gammaEnergy = peak.gammaParticleEnergy();
      snprintf( energy, sizeof(energy), "%.2f", gammaEnergy );
      
      if( peak.parentNuclide()
          && (peak.decayParticle() || peak.sourceGammaType()==PeakDef::AnnihilationGamma) )
      {
        nuclide = peak.parentNuclide()->symbol;
      }else if( peak.reaction() )
      {
        nuclide = peak.reaction()->name();
        SpecUtils::ireplace_all( nuclide, ",", " " );
      }else if( peak.xrayElement() )
      {
        nuclide = peak.xrayElement()->symbol + "-xray";
      }
    }catch( std::exception & )
    {
    }//try / catch
    
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.2f", peak.mean() );
    string meanstr = buffer;
    while( meanstr.size() < 8 )
      meanstr = " " + meanstr;
    
    snprintf( buffer, sizeof(buffer), "%.1f", peak.peakArea() );
    string areastr = buffer;
    while( areastr.size() < 10 )
      areastr = " " + areastr;
    
    snprintf( buffer, sizeof(buffer), "  %.1f", peak.peakAreaUncert() );
    string areauncertstr = buffer;
    while( areauncertstr.size() < 11 )
      areauncertstr = areauncertstr + " ";
    
    string cpststr;
    if( live_time > 0.0f )
    {
      snprintf( buffer, sizeof(buffer), "%1.4e", (peak.peakArea()/data->live_time()) );
      cpststr = buffer;
      size_t epos = cpststr.find( "e" );
      if( epos != string::npos )
      {
        if( cpststr[epos+1]!='-' && cpststr[epos+1]!='+' )
          cpststr.insert( cpststr.begin() + epos + 1, '+' );
        while( (cpststr.size()-epos) < 5 )
          cpststr.insert( cpststr.begin() + epos + 2, '0' );
      }
    }//if( data && data->live_time() > 0.0f )
    
    
    const double width = peak.gausPeak() ? (2.35482*peak.sigma()) : 0.5*peak.roiWidth();
    snprintf( buffer, sizeof(buffer), "%.2f", width );
    string widthstr = buffer;
    while( widthstr.size() < 5 )
      widthstr = " " + widthstr;
    
    snprintf( buffer, sizeof(buffer), "%.2f", (100.0*width/peak.mean()) );
    string widthprecentstr = buffer;
    widthprecentstr += "%";
    while( widthprecentstr.size() < 7 )
      widthprecentstr = " " + widthprecentstr;
    
    snprintf( buffer, sizeof(buffer), "  %.2f", peak.chi2dof() );
    string chi2str = buffer;
    while( chi2str.size() < 7 )
      chi2str = " " + chi2str;
    
    snprintf( buffer, sizeof(buffer), "  %.1f", region_area );
    string roiareastr = buffer;
    while( roiareastr.size() < 10 )
      roiareastr = " " + roiareastr;
    
    snprintf( buffer, sizeof(buffer), "  %i", int(peakn+1) );
    string numstr = buffer;
    while( numstr.size() < 3 )
      numstr = " " + numstr;
    
    specfilename = "  ";
    SpecUtils::ireplace_all( specfilename, ",", "-" );
    
    string live_time_str;
    if( live_time > 0.0f )
    {
      snprintf( buffer, sizeof(buffer), "%.3f", live_time );
      live_time_str = buffer;
    }
    
    string datestr, timestr;
    if( !meastime.is_special() )
    {
      const string tstr = SpecUtils::to_common_string( meastime, true );
      const auto pos = tstr.find(' ');
      if( pos != string::npos )
      {
        datestr = tstr.substr(0,pos);
        timestr = tstr.substr(pos+1);
      }
    }//if( !meastime.is_special() )
    
    
    response.out() << meanstr
            << ',' << areastr
            << ',' << areauncertstr
            << ',' << cpststr
            << ',' << widthstr
            << ',' << widthprecentstr
            << ',' << chi2str
            << ',' << roiareastr
            << ',' << numstr
            << ',' << specfilename
            << ',' << live_time_str
            << ',' << datestr
            << ',' << timestr
            << ',' << nuclide
            << ',' << energy
            << ',' << xlow
            << ',' << xhigh
            << eol_char;
  }//for( loop over peaks, peakn )
}//void handleRequest(...)






PeakModel::PeakModel( Wt::WObject *parent )
  : WAbstractItemModel( parent ),
    m_dataModel( NULL ),
    m_sortColumn( kMean ),
    m_sortOrder( Wt::AscendingOrder ),
    m_csvResource( NULL )
{
  m_csvResource = new PeakCsvResource( this );
}//PeakModel constructor

PeakModel::~PeakModel()
{
}//~PeakModel()



void PeakModel::setDataModel( SpectrumDataModel *dataModel )
{
  m_dataModel = dataModel;
}//void setDataModel( SpectrumDataModel *dataModel );



void PeakModel::setPeakFromSpecMeas( std::shared_ptr<SpecMeas> meas,
                                     const std::set<int> &samplenums )
{
  std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > > peaks;

  if( !!meas )
    peaks = meas->peaks( samplenums );
  
  m_measurment = meas;
  
  if( peaks == m_peaks )
    return;
  
  if( !!m_peaks && !m_peaks->empty() )
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks->size()-1) );
    m_peaks = std::make_shared< deque<std::shared_ptr<const PeakDef> > >();
    m_sortedPeaks.clear();
    endRemoveRows();
  }//if( !m_peaks->empty() )
  
  const size_t npeaksadd = !!peaks ? peaks->size() : size_t(0);

  if( !npeaksadd )
  {
    m_peaks = peaks;
  }else
  {
    boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn, meansort;
    
    shared_ptr<const SpecUtils::Measurement> data = m_dataModel ? m_dataModel->getData() : nullptr;
    sortfcn = boost::bind( &PeakModel::compare, _1, _2, m_sortColumn, m_sortOrder, data );
    meansort = boost::bind( &PeakModel::compare, _1, _2, kMean, Wt::AscendingOrder, data );

    //Make sure no null peaks
    peaks->erase( std::remove_if( peaks->begin(), peaks->end(), [](std::shared_ptr<const PeakDef> a){return !a;}), peaks->end() );
    
    
    std::sort( peaks->begin(), peaks->end(), meansort );
    
    beginInsertRows( WModelIndex(), 0, static_cast<int>(npeaksadd-1) );

    m_sortedPeaks = *peaks;
    std::stable_sort( m_sortedPeaks.begin(), m_sortedPeaks.end(), sortfcn );

    m_peaks = peaks;
    endInsertRows();
  }//if( peaks.size() )
  
  
  
}//void setPeakFromSpecMeas(...)



size_t PeakModel::npeaks() const
{
  if( !m_peaks )
    return 0;

  return m_peaks->size();
}//size_t npeaks() const


PeakModel::PeakShrdPtr PeakModel::peakPtr( const size_t peakn ) const
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before adding peak" );
  const PeakShrdPtr &peak = m_peaks->at( peakn );
  assert( peak );
  return peak;
}//PeakModel::PeakShrdPtr peakPtr( const size_t peakn ) const


PeakModel::PeakShrdPtr PeakModel::nearestPeak( double energy ) const
{
  if( !m_peaks || m_peaks->empty() )
    return PeakShrdPtr();
  
  boost::function<bool( const PeakShrdPtr &, const PeakShrdPtr &)> meansort;
  shared_ptr<const SpecUtils::Measurement> data = m_dataModel ? m_dataModel->getData() : nullptr;
  meansort = boost::bind( &PeakModel::compare, _1, _2, kMean, Wt::AscendingOrder, data );
  PeakDef *new_peak_ptr = new PeakDef();
  PeakShrdPtr peak_ptr( new_peak_ptr );
  new_peak_ptr->setMean( energy );
  deque< PeakShrdPtr >::const_iterator mean_pos
          = lower_bound( m_peaks->begin(), m_peaks->end(), peak_ptr, meansort );
  
  if( mean_pos == m_peaks->end() )
    return *(mean_pos-1);
  
  const double ahead_mean = ((mean_pos+1) != m_peaks->end())
                            ? (*(mean_pos+1))->mean()
                            : DBL_MAX;
  const double behind_mean = (mean_pos != m_peaks->begin())
                            ? (*(mean_pos-1))->mean()
                            : DBL_MAX;
  const double this_mean = (*mean_pos)->mean();
  
  const double ahead_diff = fabs( ahead_mean - energy );
  const double behind_diff = fabs( behind_mean - energy );
  const double this_diff = fabs( this_mean - energy );
  
  if( (ahead_diff < behind_diff) && (ahead_diff < this_diff) )
    return (*(mean_pos+1));
  if( behind_diff < this_diff )
    return (*(mean_pos-1));
  return (*mean_pos);
}//PeakShrdPtr nearestPeak( double energy ) const


Wt::WModelIndex PeakModel::indexOfPeak( const PeakShrdPtr &peak ) const
{
  int row = -1;  
  for( std::deque<PeakShrdPtr>::const_iterator i = m_sortedPeaks.begin();
      i != m_sortedPeaks.end(); ++i )
  {
    ++row;
    if( *i == peak )
      return PeakModel::index(row, 0);
  }
  
  return WModelIndex();
}//Wt::WModelIndex indexOfPeak( const PeakShrdPtr &peak ) const


const PeakDef &PeakModel::peak( const size_t peakn ) const
{
  return  *peakPtr(peakn);
}//const PeakDef &peak( const size_t peakn ) const


std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > PeakModel::peaks() const
{
  return m_peaks;
}//const deque<PeakDef> &peaks() const


std::vector<PeakDef> PeakModel::peakVec() const
{
  vector<PeakDef> answer;
  
  if( !m_peaks )
    return answer;
  
  deque<PeakModel::PeakShrdPtr>::const_iterator pos;
  for( pos = m_peaks->begin(); pos != m_peaks->end(); ++pos )
    answer.push_back( *(*pos) );
  
  return answer;
}//std::vector<PeakDef> peakVec() const


bool PeakModel::isWithinRange( const PeakDef &peak ) const
{
  if( !m_dataModel )
    return false;

  std::shared_ptr<const SpecUtils::Measurement> xaxis = m_dataModel->histUsedForXAxis();

  if( !xaxis )
  {
    cerr << "PeakModel::isWithinRange(...)\n\tThere is no xaxis!" << endl;
    return false;
  }

  const auto energycal = xaxis->energy_calibration();
  assert( energycal );
  
  const double lowerx = xaxis->gamma_energy_min();
  const double upperx = xaxis->gamma_energy_max();

  switch( peak.type() )
  {
    case PeakDef::GaussianDefined:
      if( peak.mean() < lowerx || peak.mean() > upperx )
        return false;
    break;
    
    case PeakDef::DataDefined:
      if( peak.lowerX() > upperx || peak.upperX() < lowerx )
        return false;
    break;
  }//switch( peak.type() )

  return true;
} // bool PeakModel::isWithinRange( const PeakDef &peak )


bool PeakModel::isOutOfRange( const PeakDef &peak ) const
{
  return !isWithinRange( peak );
}


WResource *PeakModel::peakCsvResource()
{
  return m_csvResource;
}


void PeakModel::notifySpecMeasOfPeakChange()
{
  std::shared_ptr<SpecMeas> meas = m_measurment.lock();
  if( meas )
    meas->setModified();
  else  //JIC, probably wont happen ever
    cerr << "\n\nnotifySpecMeasOfPeakChange: couldnt get lock" << endl;
}//void notifySpeakMeasOfPeakChange();


WModelIndex PeakModel::addNewPeak( const PeakDef &peak )
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before adding peak" );

  if( !isWithinRange( peak ) )
    return WModelIndex();

  notifySpecMeasOfPeakChange();
  
  PeakDef *new_peak_ptr = new PeakDef( peak );
  PeakShrdPtr peak_ptr( new_peak_ptr );

  definePeakXRange( *new_peak_ptr );
  //Need to go through and estimate the Chi2Dof here.
  
  boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn, meansort;
  shared_ptr<const SpecUtils::Measurement> data = m_dataModel ? m_dataModel->getData() : nullptr;
  sortfcn = boost::bind( &PeakModel::compare, _1, _2, m_sortColumn, m_sortOrder, data );
  meansort = boost::bind( &PeakModel::compare, _1, _2, kMean, Wt::AscendingOrder, data );

  deque< PeakShrdPtr >::iterator mean_pos
         = lower_bound( m_peaks->begin(), m_peaks->end(), peak_ptr, meansort );
  deque< PeakShrdPtr >::iterator sort_pos
       = lower_bound( m_sortedPeaks.begin(), m_sortedPeaks.end(), peak_ptr, sortfcn );
  
  const int indexpos = static_cast<int>( sort_pos - m_sortedPeaks.begin() );

  beginInsertRows( WModelIndex(), indexpos, indexpos );
  m_peaks->insert( mean_pos, peak_ptr );
  m_sortedPeaks.insert( sort_pos, peak_ptr );
  endInsertRows();
  
  return index( indexpos, 0 );
}//void addPeak( const PeakDef &peak )



void PeakModel::addPeaks( const vector<PeakDef> &peaks )
{
  for( size_t i = 0; i < peaks.size(); ++i )
    addNewPeak( peaks[i] );
}//void addPeaks( const vector<PeakDef> &peaks )


void PeakModel::definePeakXRange( PeakDef &peak )
{
  std::shared_ptr<const SpecUtils::Measurement> data, continuum;
  std::shared_ptr<PeakContinuum> peakcont = peak.continuum();
  
  if( m_dataModel )
  {
    data = m_dataModel->getData();
    
    if( peakcont->type() == PeakContinuum::External )
      continuum = peakcont->externalContinuum();
  }//if( m_dataModel )
  
  if( !peakcont->energyRangeDefined() && peakcont->isPolynomial() )
  {
    double lowerEnengy, upperEnergy;
    std::shared_ptr<const SpecUtils::Measurement> data;
    if( m_dataModel )
      data = m_dataModel->getData();
    findROIEnergyLimits( lowerEnengy, upperEnergy, peak, data );
    peakcont->setRange( lowerEnengy, upperEnergy );
  }//if( we should set the peak limits to save cpu (or rather memmorry-time) later
}//void definePeakXRange( PeakDef &peak )


void PeakModel::setPeaks( const std::vector<std::shared_ptr<const PeakDef> > &peaks )
{
  vector<PeakDef> peakcopies;
  peakcopies.reserve( peaks.size() );
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    if( peaks[i] )
      peakcopies.push_back( *peaks[i] );
  }
  
  setPeaks( peakcopies );
}//void PeakModel::setPeaks()


void PeakModel::setPeaks( vector<PeakDef> peaks )
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before adding peak" );

  if( !m_peaks->empty() )
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks->size()-1) );
    m_peaks->clear();
    m_sortedPeaks.clear();
    endRemoveRows();
  }//if( !m_peaks->empty() )

  // Remove any peaks that aren't within the range.
  const vector<PeakDef>::const_iterator peaks_end
              = remove_if( peaks.begin(), peaks.end(),
                           boost::bind( &PeakModel::isOutOfRange, this, _1 ) );
  
  
  const size_t npeaksadd = peaks_end - peaks.begin();
  
  if( npeaksadd )
  {
    for( size_t i = 0; i < npeaksadd; ++i )
      definePeakXRange( peaks[i] );
//    need to go through and estimate the chi2DOF here
    
    boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn, meansort;
    shared_ptr<const SpecUtils::Measurement> data = m_dataModel ? m_dataModel->getData() : nullptr;
    sortfcn = boost::bind( &PeakModel::compare, _1, _2, m_sortColumn, m_sortOrder, data );
    meansort = boost::bind( &PeakModel::compare, _1, _2, kMean, Wt::AscendingOrder, data );

    beginInsertRows( WModelIndex(), 0, int(npeaksadd - 1) );

    for( size_t i = 0; i < npeaksadd; ++i )
    {
      const PeakDef &inpeak = peaks[i];
      std::shared_ptr<PeakDef> peak = std::make_shared<PeakDef>( inpeak );
      m_peaks->push_back( peak );
      m_sortedPeaks.push_back( peak );
    }// for( loop over peaks )

    std::sort( m_peaks->begin(), m_peaks->end(), meansort );
    std::stable_sort( m_sortedPeaks.begin(), m_sortedPeaks.end(), sortfcn );
    
    endInsertRows();
    layoutChanged().emit();
  }//if( peaks.size() )
  
  notifySpecMeasOfPeakChange();
}//void setPeaks( const vector<PeakDef> &peaks )


void PeakModel::removePeak( const size_t peakn )  //throws if invalid peak numbers
{
  if( peakn >= m_peaks->size() )
    throw std::runtime_error( "PeakModel::removePeak(): invalid index" );

  PeakShrdPtr peak = (*m_peaks)[peakn];

  std::deque< PeakShrdPtr >::iterator sort_pos = find( m_sortedPeaks.begin(), m_sortedPeaks.end(), peak );
  assert( sort_pos != m_sortedPeaks.end() );
  const int index = static_cast<int>( sort_pos - m_sortedPeaks.begin() );

  beginRemoveRows( WModelIndex(), index, index );
  m_peaks->erase( m_peaks->begin() + peakn );
  m_sortedPeaks.erase( sort_pos );
  endRemoveRows();
  
  notifySpecMeasOfPeakChange();
}//void removePeak( const size_t peakn )


void PeakModel::removePeak( Wt::WModelIndex index )
{
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::removePeak(): non-valid index" );

  
  PeakShrdPtr peak = m_sortedPeaks.at( index.row() );
  
  std::deque< PeakShrdPtr >::iterator energy_pos = find( m_peaks->begin(), m_peaks->end(), peak );
  assert( energy_pos != m_peaks->end() );
  
  beginRemoveRows( WModelIndex(), index.row(), index.row() );
  m_peaks->erase( energy_pos );
  m_sortedPeaks.erase( m_sortedPeaks.begin() + index.row() );
  endRemoveRows();
  
  notifySpecMeasOfPeakChange();
}//void removePeak( Wt::WModelIndex index )


void PeakModel::removePeak( PeakModel::PeakShrdPtr peak )
{
  const size_t n = npeaks();
  for( size_t i = 0; i < n; ++ i )
  {
    if( peak == peakPtr(i) )
    {
      removePeak(i);
      return;
    }//if( peak == peakPtr(i) )
  }//for( size_t i = 0; i < n; ++ i )

  throw runtime_error( "PeakModel::removePeak(PeakShrdPtr):"
                       " input peak not in model" );
}//void PeakModel::removePeak( PeakModel::PeakShrdPtr peak )


void PeakModel::removeAllPeaks()
{
  setPeaks( vector<PeakDef>() );
}//void removeAllPeaks()


void PeakModel::setPeakFitFor( const Wt::WModelIndex index,
                               const PeakDef::CoefficientType coef,
                               const bool fitfor )
{
  if( !m_peaks )
    throw std::runtime_error( "PeakModel::setPeakFitFor(...): no foreground"
                              " spectrum loaded" );
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::setPeakFitFor(...): invalid index" );
  
  const int row = index.row();
  PeakShrdPtr old_peak = m_sortedPeaks.at( row );
  std::deque< PeakShrdPtr >::iterator energy_pos
                          = find( m_peaks->begin(), m_peaks->end(), old_peak );
  if( energy_pos == m_peaks->end() || (*energy_pos) != old_peak )
    throw runtime_error( "Serious logic error in setPeakFitFor(...)" );
  const size_t energy_index = energy_pos - m_peaks->begin();
  
  PeakDef newPeak = *old_peak;
  newPeak.setFitFor( coef, fitfor );
  
  m_sortedPeaks[row] = std::make_shared<PeakDef>( newPeak );
  (*m_peaks)[energy_index] = m_sortedPeaks[row];
  
  notifySpecMeasOfPeakChange();
}//setPeakFitFor(...)


void PeakModel::setContinuumPolynomialFitFor( const Wt::WModelIndex index,
                                              size_t polyCoefNum,
                                              const bool fitfor )
{
  //XXX - should refactor this function to share code with
  //      PeakModel::setPeakFitFor(...)
  if( !m_peaks )
    throw std::runtime_error( "PeakModel::setContinuumPolynomialFitFor(...): no"
                              " foreground spectrum loaded" );
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::setContinuumPolynomialFitFor(...):"
                              " invalid index" );
  
  const int row = index.row();
  PeakShrdPtr old_peak = m_sortedPeaks.at( row );
  std::deque< PeakShrdPtr >::iterator energy_pos
                           = find( m_peaks->begin(), m_peaks->end(), old_peak );
  if( energy_pos == m_peaks->end() || (*energy_pos) != old_peak )
    throw runtime_error( "Serious logic error in"
                         " setContinuumPolynomialFitFor(...)" );
  const size_t energy_index = energy_pos - m_peaks->begin();
  
  PeakDef newPeak = *old_peak;
//  newPeak.setContinuum( old_peak->continuum() );
  newPeak.continuum()->setPolynomialCoefFitFor( polyCoefNum, fitfor );
  
  m_sortedPeaks[row] = std::make_shared<PeakDef>( newPeak );
  (*m_peaks)[energy_index] = m_sortedPeaks[row];
  
  notifySpecMeasOfPeakChange();
}//void setContinuumPolynomialFitFor(...)


//Functions for the Wt::WAbstractItemModel interface
int PeakModel::rowCount( const WModelIndex & index) const
{
  //Returns # of children.  Root = size(), children = 0
  if (!index.isValid())
  { //root
    if( !m_peaks )
        return 0;
    else
      return static_cast<int>( m_peaks->size() );
  } //root
  else
  { //children
    return 0;
  } //children

}//int rowCount( const WModelIndex & ) const


int PeakModel::columnCount( const WModelIndex & ) const
{
  return PeakModel::kNumColumns;
}//int columnCount( const WModelIndex & ) const


WModelIndex PeakModel::parent( const Wt::WModelIndex & ) const
{
  return WModelIndex();
}//WModelIndex parent( const Wt::WModelIndex & ) const


const PeakModel::PeakShrdPtr &PeakModel::peak( const Wt::WModelIndex &index ) const
{
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::peak(WModelIndex): invalid input index" );

  return m_sortedPeaks.at( index.row() );
}//PeakShrdPtr peak( const Wt::WModelIndex &index ) const

boost::any PeakModel::data( const WModelIndex &index, int role ) const
{
  using boost::any;
  
  if( !m_peaks )
    return any();

  if( m_sortedPeaks.size() != m_peaks->size() )
  {
    stringstream msg;
    msg << "PeakModel::data(...)\n\tm_sortedPeaks.size()=" << m_sortedPeaks.size()
         << ", m_peaks->size()=" << m_peaks->size();
    cerr << endl << msg.str() << endl << endl;
    throw std::runtime_error( msg.str() );
  }//if( m_sortedPeaks.size() != m_peaks->size() )

  //should consider impementing ToolTipRole
  if( role != Wt::DisplayRole && role != Wt::EditRole
      && !((role==Wt::CheckStateRole) && ((index.column()==kUseForCalibration) || (index.column()==kUseForShieldingSourceFit)))
      )
    return any();


  assert( m_sortedPeaks.size() == m_peaks->size() );

  const int row = index.row();
  const int column = index.column();
  const int nrows = static_cast<int>( m_peaks->size() );

  if( row<0 || column<0 || column>=kNumColumns || row>=nrows )
    return any();


  const PeakShrdPtr &peak = m_sortedPeaks[row];

  
  auto getPeakArea = [&peak,this]() -> boost::any {
    switch( peak->type() )
    {
      case PeakDef::GaussianDefined:
        return peak->peakArea();
        
      case PeakDef::DataDefined:
      {
        if( !m_dataModel )
          return boost::any();
        
        double contArea = 0.0;
        std::shared_ptr<const SpecUtils::Measurement> dataH = m_dataModel->getData();
        
        if( !dataH )
          return boost::any();
        
        if( peak->continuum()->defined() )
        {
          double lowx(0.0), upperx(0.0);
          findROIEnergyLimits( lowx, upperx, *peak, dataH );
          contArea = peak->offset_integral( lowx, upperx );
        }else
        {
          std::shared_ptr<const SpecUtils::Measurement> continuum = m_dataModel->getBackground();
          if( continuum )
          {
            size_t lower_channel, upper_channel;
            estimatePeakFitRange( *peak, continuum, lower_channel, upper_channel );
            contArea = continuum->gamma_channels_sum(lower_channel, upper_channel);
          }
        }//if( peak->continuumDefined() ) / else
        
        size_t lower_channel, upper_channel;
        estimatePeakFitRange( *peak, dataH, lower_channel, upper_channel );
        const double dataArea = dataH->gamma_channels_sum(lower_channel, upper_channel);
        
        return (dataArea - contArea);
      }//case PeakDef::DataDefined:
    }//switch( peak->type() )
  };//auto getPeakArea lambda
  
  switch( column )
  {
    case kMean:
      return peak->mean();
      
    case kFwhm:
      switch( peak->type() )
      {
        case PeakDef::GaussianDefined:
          return 2.3548201*peak->sigma();
          
        case PeakDef::DataDefined:
          return boost::any();
      }//switch( peak->type() )
      
    case kAmplitude:
      return getPeakArea();
      
    case PeakModel::kCps:
    {
      boost::any areaAny = getPeakArea();
      
      if( areaAny.empty() )
        return areaAny;
      
      const double area = boost::any_cast<double>(areaAny);
      std::shared_ptr<const SpecUtils::Measurement> dataH = m_dataModel->getData();
      const float liveTime = dataH ? dataH->live_time() : -1.0f;
      if( liveTime <= 0.0f )
        return boost::any();
      
      const double cps = area / liveTime;
      const double uncert = peak->amplitudeUncert();
      
      char text[64];
      if( (peak->type() != PeakDef::GaussianDefined) || (uncert <= 0.0) )
      {
        snprintf( text, sizeof(text), "%.4g", cps );
        return WString(text);
      }
      const double cpsUncert = uncert / liveTime;
      const string txt = PhysicalUnits::printValueWithUncertainty( cps, cpsUncert, 4 );
      return WString::fromUTF8(txt);
    }//case PeakModel::kCps:
      
    case kIsotope:
    {
      const SandiaDecay::Nuclide *parentNuclide = peak->parentNuclide();

      if( parentNuclide )
        return WString( parentNuclide->symbol );
      else if( peak->xrayElement() )
        return WString( peak->xrayElement()->symbol + " xray" );
      else if( peak->reaction() )
        return WString( peak->reaction()->name() );
      
      return boost::any();
    }//case kIsotope:

    case kDifference:
    {
      try
      {
        const float energy = peak->gammaParticleEnergy() - peak->mean();
        char text[32];
        snprintf( text, sizeof(text), "%.2f keV", (floor(100.0*energy + 0.5)/100.0) );
        return WString( text );
      }catch( std::exception & )
      {
        return boost::any();
      }
    }
      
    case kPhotoPeakEnergy:
    {
      try
      {
        float energy = peak->gammaParticleEnergy();
        const char *prefix = "", *postfix = "";
        
        switch( peak->sourceGammaType() )
        {
          case PeakDef::NormalGamma:
          case PeakDef::AnnihilationGamma:
            break;
          case PeakDef::SingleEscapeGamma:
            prefix = "S.E. ";
            if( peak->decayParticle() )
              energy = peak->decayParticle()->energy;
          break;
          case PeakDef::DoubleEscapeGamma:
            prefix = "D.E. ";
            if( peak->decayParticle() )
              energy = peak->decayParticle()->energy;
          break;
          case PeakDef::XrayGamma:
            postfix = " xray";
          break;
        }//switch( peak->sourceGammaType() )
        
        char text[32];
        snprintf( text, sizeof(text), "%s%.2f keV%s", prefix, (floor(100.0*energy + 0.5)/100.0), postfix );
        return WString( text );
      }catch( std::exception & )
      {
        return boost::any();
      }
    }//case kPhotoPeakEnergy:

    case kUseForShieldingSourceFit:
      // Make so we only "use for shielding/source fit" show checkbox for decay gammas and x-rays,
      //  and not for peaks with no source associated, or with florescence x-rays, or for reactions
      if( !peak->nuclearTransition() || !peak->parentNuclide() || (peak->decayParticleIndex() < 0) )
        return boost::any();
      
      switch( peak->sourceGammaType() )
      {
        case PeakDef::NormalGamma:
        case PeakDef::XrayGamma:
        case PeakDef::AnnihilationGamma:
          break;
        
        case PeakDef::SingleEscapeGamma:
        case PeakDef::DoubleEscapeGamma:
          // Don't show "use for shielding/source fit" checkbox for single and double escape peaks
          return boost::any();
          break;
      }//switch( peak->sourceGammaType() )
      
      return peak->useForShieldingSourceFit();

    case kUseForCalibration:
      return peak->useForCalibration();
      
    case kPeakLineColor:
    {
      if( peak->lineColor().isDefault() )
        return boost::any();
      return boost::any( WString::fromUTF8(peak->lineColor().cssText(false)) );
    }//case kPeakLineColor:
      
    case kUserLabel:
    {
      if( peak->userLabel().empty() )
        return boost::any();
      return boost::any( WString::fromUTF8(peak->userLabel()) );
    }//case kUserLabel:
      
    case kHasSkew:
    {
      if( peak->skewType() == PeakDef::LandauSkew )
        return WString( "True" );
      else
        return WString( "False" );
    }//case kHasSkew:
      
    case kSkewAmount:
    {
      if( peak->skewType() == PeakDef::LandauSkew )
      {
        char text[64];
        snprintf( text, sizeof(text), "%.3f, W%.3f, %.3f",
                  peak->coefficient(PeakDef::LandauAmplitude),
                  peak->coefficient(PeakDef::LandauMode),
                  peak->coefficient(PeakDef::LandauSigma) );
        return WString( text );
      }else
      {
        return WString( "NA" );
      }
    }//case kSkewAmount:
      
    case kType:
    {
      if( peak->gausPeak() )
        return WString( "Gaussian" );
      else
        return WString( "Region" );
    }
      
    case kLowerX:
    case kUpperX:
    {
      if( !m_dataModel )
        return boost::any();
      std::shared_ptr<const SpecUtils::Measurement> dataH = m_dataModel->getData();
      double lowx(0.0), upperx(0.0);
      findROIEnergyLimits( lowx, upperx, *peak, dataH );
      if( column == kLowerX )
        return lowx;
      else
        return upperx;
    }//case kLowerX / case kUpperX:
    
      
    case kRoiCounts:
    {
      if( !m_dataModel )
        return boost::any();
      std::shared_ptr<const SpecUtils::Measurement> dataH = m_dataModel->getData();
      if( !dataH )
        return boost::any();
      
      double lowx(0.0), upperx(0.0);
      findROIEnergyLimits( lowx, upperx, *peak, dataH );
      return dataH->gamma_integral( lowx, upperx );
    }//case kRoiCounts:
      
    case kContinuumType:
    {
      switch( peak->continuum()->type() )
      {
        case PeakContinuum::NoOffset:   return WString( "None" );
        case PeakContinuum::Constant:   return WString( "Constant" );
        case PeakContinuum::Linear:     return WString( "Linear" );
        case PeakContinuum::Quadratic:  return WString( "Quadratic" );
        case PeakContinuum::Cubic:      return WString( "Cubic" );
        case PeakContinuum::External:   return WString( "Global Cont." );
      }//switch( peak->offsetType() )
      
      return boost::any();
    }//case kContinuumType:
      
    case kNumColumns:
      return any();
  }//switch( section )

  return any();
}//boost::any PeakModel::data( const WModelIndex &index, int role ) const


PeakModel::SetGammaSource PeakModel::setNuclideXrayReaction( PeakDef &peak,
                                                            std::string label,
                                                    const double nsigma_window )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  if( !db )
    return NoSourceChange;
  
  SpecUtils::to_lower_ascii( label );
  
  const bool hadSource = (peak.parentNuclide() || peak.xrayElement() || peak.reaction() );
  
  if( label.empty()
     || SpecUtils::contains( label, "none" )
     || SpecUtils::contains( label, "undef" )
     || label == "na" )
  {
    const bool use = peak.useForShieldingSourceFit();
    peak.clearSources();
    
    if( !hadSource )
      return NoSourceChange;
    
    if( use )
      return SourceAndUseChanged;
    
    return SourceChange;
  }//if( getting rid of source )
  
  double ref_energy = -1.0;
  const SandiaDecay::Nuclide *nuclide = NULL;
  
  PeakDef::SourceGammaType srcType;
  PeakDef::gammaTypeFromUserInput( label, srcType );
  

  //If there are 2 non-connected strings of numbers, and the second one is
  //  not followed by a m or meta, but rather nothingn or ev, kev, or MeV,
  //  then we will assume the second one is an energy
  //XXX - this assumes we will never deal with meta-stable states higher
  //      than 1 - e.g. if user types int 2m Co 60 55.1 keV, we will mess
  //      things up
  size_t number_start, number_stop = string::npos;
  number_start = label.find_first_of( "0123456789" );
  if( number_start != string::npos )
    number_stop = label.find_first_not_of( "0123456789", number_start );
  if( number_stop != string::npos )
  {
    number_start = label.find_first_of( "0123456789", number_stop );
    if( number_start != string::npos )
    {
      number_stop = label.find_first_not_of( ".0123456789", number_start );
      if( number_stop == string::npos )
        number_stop = label.size();
    }//if( number_start != string::npos )
  }//if( number_stop != string::npos )
  
  if( (number_start == string::npos) || (number_stop == string::npos) )
  {
    ref_energy = peak.mean();
    nuclide = db->nuclide( label );
  }else
  {
    const string nuclabel = label.substr( 0, number_start );
    nuclide = db->nuclide( nuclabel );
    
    if( number_stop == string::npos )
      number_stop = label.size();
    string energy( label.begin()+number_start, label.begin()+number_stop );
    ref_energy = std::stod( energy );
    
    if( number_stop < label.size() )
    {
      if( label.find( "mev", number_stop ) != string::npos )
        ref_energy *= 1000.0;
      else if( label.find( "ev", number_stop ) != string::npos
              && label.find( "kev", number_stop ) == string::npos )
        ref_energy /= 1000.0;
      // else assume keV
    }//if( number_stop < label.size() )
  }//if( number_start == string::npos )
  
  
  if( nuclide )
  {
    size_t transition_index = 0;
    const SandiaDecay::Transition *transition = NULL;
    
    const double width = peak.gausPeak() ? nsigma_window*peak.sigma() : (nsigma_window/8.0)*peak.roiWidth();
    
    double extraEnergy = 0.0;
    switch( srcType )
    {
      case PeakDef::NormalGamma:
      case PeakDef::AnnihilationGamma:
      case PeakDef::XrayGamma:
        break;
        
      case PeakDef::SingleEscapeGamma: extraEnergy = 510.99891;     break;
      case PeakDef::DoubleEscapeGamma: extraEnergy = 2.0*510.99891; break;
    }//switch( srcType )
    
    const bool xrayOnly = (srcType == PeakDef::SourceGammaType::XrayGamma);
    PeakDef::SourceGammaType sourceGammaType = srcType;
    
    PeakDef::findNearestPhotopeak( nuclide, ref_energy + extraEnergy, width,
                                    xrayOnly, transition, transition_index, sourceGammaType );
    
    //There wasnt any photopeaks within 4 sigma, so instead we'll just use
    //  the closest photpopeak
    if( !transition && (sourceGammaType!=PeakDef::AnnihilationGamma)
        && nsigma_window>=0.0 )
      PeakDef::findNearestPhotopeak( nuclide, ref_energy + extraEnergy, -1.0, xrayOnly,
                                    transition, transition_index, sourceGammaType );
    
    switch( srcType )
    {
      case PeakDef::NormalGamma:
      case PeakDef::AnnihilationGamma:
      case PeakDef::XrayGamma:
      break;
        
      case PeakDef::SingleEscapeGamma: case PeakDef::DoubleEscapeGamma:
        if( transition )
          sourceGammaType = srcType;
        else if( sourceGammaType == PeakDef::AnnihilationGamma )
          transition = NULL;  //dont expect this to ever happen
      break;
    }//switch( srcType )
    
    
    if( !transition && (sourceGammaType!=PeakDef::AnnihilationGamma) )
    {
      peak.clearSources();
      return (hadSource ? SourceChange : NoSourceChange);
    }else
    {
      bool changedFit = false, shouldFit = false;
      
      switch( sourceGammaType )
      {
        case PeakDef::NormalGamma:
          shouldFit = recomendUseForFit( nuclide, transition->products[transition_index].energy );
        break;
        
        case PeakDef::AnnihilationGamma:
          shouldFit = recomendUseForFit( nuclide, 510.99891f );
        break;
          
        case PeakDef::SingleEscapeGamma:
        case PeakDef::DoubleEscapeGamma:
        case PeakDef::XrayGamma:
          shouldFit = false;
        break;
      }//switch( srcType )

      
      changedFit |= (shouldFit==peak.useForShieldingSourceFit());
      changedFit |= (shouldFit && !peak.nuclearTransition());
      
      peak.useForShieldingSourceFit( shouldFit );
      peak.setNuclearTransition( nuclide, transition, int(transition_index), sourceGammaType );
      
      return (changedFit ? SourceAndUseChanged : SourceChange);
    }
  }else if( srcType == PeakDef::SourceGammaType::XrayGamma )
  {
      //Assume input like "fe xray 98.2 kev", which will have been transformed
      //  to "fe 98.2 kev" by PeakDef::gammaTypeFromUserInput
      size_t alpha_start = label.size(), alpha_end = label.size();
      for( size_t i = 0; i < label.size(); ++i )
        if( label[i]>='a' && label[i]<='z' ) { alpha_start = i; break; }
      
      for( size_t i = alpha_start; i < label.size(); ++i )
        if( label[i]<'a' || label[i]>'z' ){alpha_end = i; break;}
      
      string elname;
      if( alpha_start != alpha_end )
        elname = string(label.begin()+alpha_start, label.begin()+alpha_end);
      
      const SandiaDecay::Element *el = db->element( elname );
      if( !el )
      {
        peak.clearSources();
        return (hadSource ? SourceChange : NoSourceChange);
      }//if( !el )
      
      double energy = 0.0;
      const size_t numstart = label.find_first_of( "0123456789" );
      if( !(stringstream( string(label.c_str() + numstart) ) >> energy) )
        energy = peak.mean();
      
      const SandiaDecay::EnergyIntensityPair *nearXray
                                      = PeakDef::findNearestXray( el, energy );
      energy = (nearXray ? nearXray->energy : 0.0);
      peak.setXray( el, energy );
      
      return ((!nearXray && !hadSource) ? NoSourceChange : SourceChange);
    }else
    {
      //lets check for a reaction
      const ReactionGamma *rctndb = ReactionGammaServer::database();
      string rctstr = label;
      size_t pos = label.find_first_of( ')' );
      if( pos != string::npos )
        rctstr = label.substr( 0, pos+1 );
      vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
      try
      {
        rctndb->gammas( rctstr, possible_rctns );
        if( possible_rctns.empty() )
          return NoSourceChange;
      }catch(...)
      {
        return NoSourceChange;
      }
      
      //XXX - just taking first reaction, however there may be multiple
      const ReactionGamma::Reaction *rctn = possible_rctns[0].reaction;
      
      pos = label.find_first_of( "0123456789", pos );
      string energystr;
      if( pos != string::npos )
        energystr = label.substr( pos );
      double energy = 0.0;
      if( rctn )
      {
        if( !(stringstream(energystr) >> energy) )
          energy = peak.mean();
        double nearestE = -9999999.9;
        for( const ReactionGamma::EnergyAbundance &eip : rctn->gammas )
          if( fabs(eip.energy-energy) < fabs(nearestE-energy) )
            nearestE = eip.energy;
        energy = nearestE;
      }//if( (stringstream(energystr) >> energy) )
      
      peak.setReaction( rctn, static_cast<float>(energy), srcType );
      
      return (rctn ? SourceChange : NoSourceChange);
  }//if( nuclide ) / else if( xray ) / else
  
  return NoSourceChange;
}//bool PeakModel::setNuclideXrayReaction( PeakDef &peak, std::string )


bool PeakModel::setData( const WModelIndex &index,
                         const boost::any &value, int role )
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before setting peak data" );

  notifySpecMeasOfPeakChange();
  
  try
  {
    if( (role != Wt::EditRole)
        && (role != Wt::CheckStateRole) )
      return false;

    const int row = index.row();
    const int column = index.column();
    const int nrows = static_cast<int>( m_peaks->size() );

    bool changedFit = false;
    
    if( value.empty() || row < 0 || column < 0 || column>=kNumColumns || row >= nrows )
      return false;

    switch( column )
    {
      case kMean: case kFwhm: case kAmplitude:
      case kIsotope:
      case kPhotoPeakEnergy:
      case kCandidateIsotopes:
      case kUseForCalibration:
      case kUseForShieldingSourceFit:
      case kUserLabel:
      case kPeakLineColor:
      break;

      case kCps: case kHasSkew: case kSkewAmount: case kType: case kLowerX: case kUpperX:
      case kRoiCounts: case kContinuumType: case kNumColumns: case kDifference:
      default:
        cerr << "PeakModel::setData(...)\n\tUn Supported column" << endl;
        return false;
    }//switch( section )

    double dbl_val = 0.0;
    WString txt_val;

    try
    {
      txt_val = boost::any_cast<WString>( value );
    }catch(...)
    {}

    switch( column )
    {
      case kMean: case kFwhm: case kAmplitude:
        try
        {
          dbl_val = std::stod( txt_val.toUTF8() );
        }catch(...)
        {
          cerr << "PeakModel::setData(...)\n\tUnable to convert '" << txt_val
               << "' to a dbl" << endl;
          return false;
        }//try / catch
      break;
      default:
        break;
    }//switch( column )

    
    PeakShrdPtr old_peak = m_sortedPeaks[row];
    std::deque< PeakShrdPtr >::iterator energy_pos
                           = find( m_peaks->begin(), m_peaks->end(), old_peak );
    assert( energy_pos != m_peaks->end() );

    PeakDef new_peak( *old_peak );

    switch( column )
    {
      case kMean:
        switch( new_peak.type() )
        {
          case PeakDef::GaussianDefined:
            new_peak.setMean( dbl_val );
          break;
            
          case PeakDef::DataDefined:
            return false;
          break;
        }//switch( new_peak.type() )
      break;

      case kFwhm:
        switch( new_peak.type() )
        {
          case PeakDef::GaussianDefined:
            new_peak.setSigma( dbl_val / 2.3548201 );
          break;
          
          case PeakDef::DataDefined:
            return false;
          break;
        }//switch( new_peak.type() )
      break;

      case kAmplitude:
        switch( new_peak.type() )
        {
          case PeakDef::GaussianDefined:
            new_peak.setAmplitude( dbl_val );
          break;
          
          case PeakDef::DataDefined:
            return false;
          break;
        }//switch( new_peak.type() )
      break;

      case kIsotope:
      {
        const std::string srctxt = txt_val.toUTF8();
        const SetGammaSource result
                   = setNuclideXrayReaction( new_peak, srctxt, 4.0 );
        changedFit |= (result==SourceAndUseChanged);
        
        //Lambda to see if any peaks, other than the ignorePeak, that has the
        //  same nuc/xray/rctn assigned as srcPeak has a color, and if so
        //  return that color.  If the source has more than one color of peaks
        //  then dont return a color.
        auto sameSrcColor = [this]( const PeakDef &srcPeak, const PeakShrdPtr &ignorePeak ) -> WColor {
          if( !m_peaks || !ignorePeak )
            return WColor();
          
          vector<WColor> src_colors;
          for( const auto &p : *m_peaks )
          {
            if( p
               && p != ignorePeak
               && (srcPeak.parentNuclide() || srcPeak.xrayElement() || srcPeak.reaction())
               && srcPeak.parentNuclide()==p->parentNuclide()
               && srcPeak.xrayElement()==p->xrayElement()
               && srcPeak.reaction()==p->reaction()
               && !p->lineColor().isDefault() )
            {
              if( std::find( begin(src_colors), end(src_colors),p->lineColor()) == end(src_colors) )
                src_colors.push_back( p->lineColor() );
            }
          }//for( const auto &p : *m_peaks )
          
          return src_colors.size()==1 ? *begin(src_colors) : WColor();
        };//sameSrcColor lambda
        
        
        if( !new_peak.hasSourceGammaAssigned() && old_peak->hasSourceGammaAssigned()
            && (result != NoSourceChange)  )
        {
          if( old_peak->lineColor() == sameSrcColor(*old_peak,old_peak) )
            new_peak.setLineColor( WColor() );
        }else if( (result != NoSourceChange) )
        {
          const auto oldsrccolor = sameSrcColor(*old_peak,old_peak);
          const auto newcolor = sameSrcColor(new_peak,old_peak);
          
          //If the old peak had a different color than its source - then dont change the color.
          //If
          
          if( (!oldsrccolor.isDefault() && oldsrccolor==old_peak->lineColor()) || !newcolor.isDefault() )
          {
            if( oldsrccolor.isDefault() || old_peak->lineColor().isDefault() || oldsrccolor==old_peak->lineColor() )
              new_peak.setLineColor( newcolor );
          }
        }
        
        break;
      }//case kIsotope:

      case kPhotoPeakEnergy:
      {
        const SandiaDecay::Transition *trans = new_peak.nuclearTransition();
        const SandiaDecay::Element *el = new_peak.xrayElement();
        const ReactionGamma::Reaction *rctn = new_peak.reaction();
        
        string text = SpecUtils::to_lower_ascii_copy( txt_val.narrow() );

        PeakDef::SourceGammaType srcType;
        PeakDef::gammaTypeFromUserInput( text, srcType );
        
        double unit = 1.0;
        if( text.find("mev") != string::npos )
          unit = 1000.0;
        else if( (text.find("kev") == string::npos)
                 && (text.find("ev") != string::npos) )
          unit = 0.001;

        SpecUtils::trim( text );
        stringstream convertstr( text );

        double energy = -999.0;
        if( !(convertstr >> energy) )
          return false;
        
        energy *= unit;
        
        if( trans && trans->parent )
        {
          size_t trans_index = 0;
          const SandiaDecay::Nuclide *nuclide = new_peak.parentNuclide();
          const SandiaDecay::Transition *transition = NULL;
        
          //The 0.0 below means find the actual closest in energy, and not the
          //  most likely
          PeakDef::SourceGammaType sourceGammaType;
          const bool xrayOnly = (srcType == PeakDef::SourceGammaType::XrayGamma);
          PeakDef::findNearestPhotopeak( nuclide, energy, 0.0, xrayOnly,
                                transition, trans_index, sourceGammaType );
          
          if( !transition && (sourceGammaType!=PeakDef::AnnihilationGamma) )
            return false;
          
          switch( srcType )
          {
            case PeakDef::NormalGamma:
            case PeakDef::AnnihilationGamma:
            case PeakDef::XrayGamma:
              break;
            case PeakDef::SingleEscapeGamma:
            case PeakDef::DoubleEscapeGamma:
              if( energy > 1022.0 )
                sourceGammaType = srcType;
              break;
          }//switch( srcType )
          
          new_peak.setNuclearTransition( new_peak.parentNuclide(), transition,
                                         int(trans_index), sourceGammaType );
          
          bool shouldFit = false;
          switch( srcType )
          {
            case PeakDef::NormalGamma:
              shouldFit = recomendUseForFit( nuclide, transition->products[trans_index].energy );
              break;
            case PeakDef::AnnihilationGamma:
              shouldFit = recomendUseForFit( nuclide, 510.99891f );
              break;
            case PeakDef::SingleEscapeGamma:
            case PeakDef::DoubleEscapeGamma:
              break;
            case PeakDef::XrayGamma:
              shouldFit = false;
              break;
          }//switch( srcType )
          
          changedFit |= (shouldFit != old_peak->useForShieldingSourceFit());
          new_peak.useForShieldingSourceFit( shouldFit );
        }else if( el )
        {
          double nearestE = -999.9;
          for( const SandiaDecay::EnergyIntensityPair &eip : el->xrays )
            if( fabs(eip.energy-energy) < fabs(nearestE-energy) )
              nearestE = eip.energy;
          
          if( fabs(nearestE-energy) < 5.0 )
            energy = nearestE;
          else
            return false;
          
          new_peak.setXray( el, energy );
        }else if( rctn )
        {
          double nearestE = -999.9;
          for( const ReactionGamma::EnergyAbundance &eip : rctn->gammas )
            if( fabs(eip.energy-energy) < fabs(nearestE-energy) )
              nearestE = eip.energy;
          if( nearestE > 0.0 )
            new_peak.setReaction( rctn, static_cast<float>(nearestE), srcType );
        }else
        {
          return false;
        }
        
        break;
      }//case kPhotoPeakEnergy:

      case kCandidateIsotopes:
      {
        try
        {
          const vector<PeakDef::CandidateNuclide> candidates = boost::any_cast< vector<PeakDef::CandidateNuclide> >( value );
          new_peak.setCandidateNuclides( candidates );
        }catch(...)
        {
          return false;
        }

        break;
      }//case kCandidateIsotopes:

      case kUseForShieldingSourceFit:
      {
        try
        {
          const bool use = boost::any_cast<bool>( value );
          
//          cerr << "useForShieldingSourceFit=" << use << " type="
//               << value.type().name() << endl;
//          boost::any myfalsebool( bool(false) ), mytrueebool( bool(true) );
//          cerr << "myfalsebool value=" << boost::any_cast<bool>( myfalsebool)
//               << " mytrueebool value=" << boost::any_cast<bool>( mytrueebool)
//               << " and typename=" << myfalsebool.type().name() << endl;
          if( use && !new_peak.parentNuclide() )
            passMessage( "You must associate a nuclide with the"
                        " peak before using it for shielding/source-fitting." ,
                        "", WarningWidget::WarningMsgHigh );
          
          if( new_peak.useForShieldingSourceFit() == use )
            return false;
          
          const SandiaDecay::RadParticle *radpart = new_peak.decayParticle();
          if( use && radpart && (radpart->type == SandiaDecay::XrayParticle) )
          {
            passMessage( "Warning: using x-rays for fitting source nuclides " \
                         "is not usually a great idea, so please use caution", \
                         "", WarningWidget::WarningMsgLow );
          }
          
          new_peak.useForShieldingSourceFit( use );
        }catch(...)
        {
          cerr << "\n\tFailed in casting checkbox" << endl;
          return false;
        }

        break;
      }//case kUseForShieldingSourceFit:

      case kUseForCalibration:
      {
        try
        {
          const bool use = boost::any_cast<bool>( value );
          if( use && !new_peak.xrayElement() && !new_peak.parentNuclide()
              && !new_peak.reaction() )
            passMessage( "You must associate a nuclide, xray, or reaction with"
                         " the peak before using it for calibrarion" ,
                         "", WarningWidget::WarningMsgHigh );
          
          if( use == new_peak.useForCalibration() )
            return false;
          new_peak.useForCalibration( use );
        }catch(...)
        {
          cerr << "\n\tFailed in casting checkbox" << endl;
          return false;
        }

        break;
      }//case kUseForCalibration:

      case kPeakLineColor:
      {
        const string css_color = txt_val.toUTF8();
        if( SpecUtils::iequals_ascii(css_color, "none")
            || SpecUtils::iequals_ascii(css_color, "na")
            || css_color.empty() )
          new_peak.setLineColor( Wt::WColor() );
        else
        {
          try{ new_peak.setLineColor( Wt::WColor(css_color) ); }catch(...){ return false; }
        }
        break;
      }//case kPeakLineColor:
        
      case kUserLabel:
      {
        new_peak.setUserLabel( txt_val.toUTF8() );
        break;
      }
        
      case kHasSkew: case kSkewAmount: case kType:
      case kLowerX: case kUpperX: case kRoiCounts:
      case kCps: case kContinuumType:
      case kNumColumns:
        return false;
    }//switch( section )

    m_sortedPeaks[row] = std::make_shared<const PeakDef>( new_peak );
    (*energy_pos) = m_sortedPeaks[row];

    if( column != kIsotope
        && column != kPhotoPeakEnergy
        && column != kMean )
    {
      dataChanged().emit( index, index );
    }else if( column == kPhotoPeakEnergy )
    {
      if( changedFit )
        dataChanged().emit( index, PeakModel::index(row, kUseForShieldingSourceFit) );
      else
        dataChanged().emit( index, PeakModel::index(row, kDifference) );
    }else if( column == kMean )
    {
      if( new_peak.decayParticle() || new_peak.parentNuclide() || new_peak.xrayElement() )
        dataChanged().emit( index, PeakModel::index(row, kPhotoPeakEnergy) );
      else
        dataChanged().emit( index, index );
    }else
    {
      dataChanged().emit( index, PeakModel::index(row, kUserLabel) );
    }


    if( column == m_sortColumn )
      sort( m_sortColumn, m_sortOrder );

    return true;
  }catch( std::exception &e )
  {
    cerr << "PeakModel::setData(...)\n\tUnexpected exception: " << e.what() << endl;
  }
  return false;
}//bool setData(...)


WFlags<ItemFlag> PeakModel::flags( const WModelIndex &index ) const
{
  switch( index.column() )
  {
    case kMean: case kFwhm: case kAmplitude: case kUserLabel:
    case kIsotope: case kPhotoPeakEnergy:
#if( ALLOW_PEAK_COLOR_DELEGATE )
    case kPeakLineColor:
#endif
      return ItemIsEditable | ItemIsSelectable; //ItemIsSelectabl

    case kUseForShieldingSourceFit:
    case kUseForCalibration:
      return ItemIsUserCheckable | ItemIsSelectable;

    case kHasSkew: case kSkewAmount: case kType:
    case kLowerX: case kUpperX: case kRoiCounts:
    case kContinuumType: case kCps:
#if( !ALLOW_PEAK_COLOR_DELEGATE )
    case kPeakLineColor:
#endif
    case kNumColumns:
      return ItemIsSelectable;
  }//switch( section )

  return WFlags<ItemFlag>();
}//WFlags<ItemFlag> flags( const WModelIndex &index )


WModelIndex PeakModel::index( int row, int column, const WModelIndex & ) const
{
  return  WAbstractItemModel::createIndex( row, column, (void *)0 );
}//WModelIndex PeakModel::index( int row, int column, const WModelIndex & ) const



boost::any PeakModel::headerData( int section, Orientation orientation, int role ) const
{
  //When orientation is Horizontal, section is a column number,
  //  when orientation is Vertical, section is a row (peak) number.

  if( role == LevelRole )
  {
    return 0;
  } //LevelRole
  else if (role == DisplayRole)
  {
    //If we are here, we want the column title
    switch( section )
    {
      case kMean:           return boost::any( WString("Mean") );
      case kFwhm:           return boost::any( WString("FWHM") ); //\x03C3
      case kAmplitude:      return boost::any( WString("Area") );
      case kCps:            return boost::any( WString("CPS") );
      case kIsotope:        return boost::any( WString("Nuclide") );
      case kPhotoPeakEnergy:return boost::any( WString("Photopeak") );
      case kDifference:     return boost::any( WString("Difference") );
      case kUseForShieldingSourceFit: return boost::any( WString("Use") );
      case kCandidateIsotopes: return boost::any();
      case kUseForCalibration: return boost::any( WString("Calib. Peak") );
      case kUserLabel:      return boost::any( WString("Label") );
      case kPeakLineColor:  return boost::any( WString("Color") );
      case kHasSkew:        return boost::any( WString("Skew") );
      case kSkewAmount:     return boost::any( WString("Skew Amp") );
      case kType:           return boost::any( WString("Peak Type") );
      case kLowerX:         return boost::any( WString("Lower Energy") );
      case kUpperX:         return boost::any( WString("Upper Energy") );
      case kRoiCounts:      return boost::any( WString("ROI Counts") );
      case kContinuumType:  return boost::any( WString("Cont. Type") );
      case kNumColumns:     return boost::any();
    }//switch( section )
  } //DisplayRole
  else if (role == ToolTipRole)
  {
    switch( section )
    {
      case kMean:           return boost::any( WString("Mean Energy") );
      case kFwhm:           return boost::any( WString("Full Width at Half Maximum") ); //\x03C3
      case kAmplitude:      return boost::any( WString("Peak Area") );
      case kCps:            return boost::any( WString("Peak counts per second") );
      case kIsotope:        return boost::any();
      case kPhotoPeakEnergy:return boost::any( WString("Photopeak Energy") );
      case kDifference:     return boost::any( WString("Difference between photopeak and mean energy") );
      case kUseForShieldingSourceFit: return boost::any();
      case kCandidateIsotopes: return boost::any();
      case kUseForCalibration: return boost::any();
      case kPeakLineColor:     return boost::any( WString("Peak color") );
      case kUserLabel:         return boost::any( WString("User specified label") );
      case kRoiCounts:         return boost::any( WString("Integral of gamma counts over the region of interest") );
      case kHasSkew:
      case kSkewAmount:
      case kType:
      case kLowerX:
      case kUpperX:
      case kContinuumType:
      case kNumColumns:     return boost::any();
    }//switch( section )
    
  } //ToolTipRole
  else if( (orientation != Horizontal) || (role != DisplayRole) )
    return WAbstractItemModel::headerData( section, orientation, role );

  return boost::any();
}//any headerData( int section, Orientation orientation, int role ) const


bool PeakModel::removeRows( int row, int count, const WModelIndex & )
{
  if( !m_peaks )
    throw runtime_error( "Can not remove rows without a primary spectrum" );

  const int nrow = static_cast<int>( m_peaks->size() );

  if( !count )
    return false;

  if( (row >= nrow) || ((row+count) > nrow) || (count < 0) )
    throw runtime_error( "PeakModel::removeRows(...): invalid intput indexes" );

  beginRemoveRows( WModelIndex(), row, row+count-1 );

  deque< PeakShrdPtr >::iterator index, mean_index;
  const deque< PeakShrdPtr >::iterator start = m_sortedPeaks.begin() + row;
  const deque< PeakShrdPtr >::iterator end = start + count;

  for( index = start; index != end; ++index )
  {
    const PeakShrdPtr &peak = *index;
    mean_index = find( m_peaks->begin(), m_peaks->end(), peak );
    assert( mean_index != m_peaks->end() );
    m_peaks->erase( mean_index );
  }//for( index = start; index != end; ++index )

  m_sortedPeaks.erase( start, end );

  endRemoveRows();
  
  notifySpecMeasOfPeakChange();

  return true;
}//bool removeRows( int row, int count, const WModelIndex &)


//removeColumns(...) should not be used and will throw std::runtime_error
//  if called.
bool PeakModel::removeColumns( int, int, const Wt::WModelIndex & )
{
  throw std::runtime_error( "PeakModel::removeColumns(...): should not be called" );
  return false;
}//bool removeColumns( int, int, const Wt::WModelIndex & )


void PeakModel::sort( int col, Wt::SortOrder order )
{
  using boost::bind;

  if( m_sortedPeaks.empty() )
    return;

  m_sortOrder = order;
  m_sortColumn = Columns(col);

  if( m_sortColumn==kNumColumns )
    m_sortColumn = kMean;

  boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn;
  shared_ptr<const SpecUtils::Measurement> data = m_dataModel ? m_dataModel->getData() : nullptr;
  sortfcn = boost::bind( &PeakModel::compare, _1, _2, m_sortColumn, order, data );

  layoutAboutToBeChanged();
  stable_sort( m_sortedPeaks.begin(), m_sortedPeaks.end(), sortfcn );

  const int nrow = rowCount();
  const int ncol = columnCount();
  
  WModelIndex ul = index(0,0);
  WModelIndex lr = index(nrow-1,ncol-1);
  
  dataChanged().emit(ul,lr);
  
//XXX - I wouldnt think the next few lines should be necessary, but it appears
//      that it is
//  layoutChanged();
  
//  beginRemoveRows( WModelIndex(), 0, nrow-1 );
//  std::deque< PeakShrdPtr > peaks, sortedPeaks;
//  m_peaks->swap( peaks );
//  sortedPeaks.swap( m_sortedPeaks );
//  endRemoveRows();
//  beginInsertRows( WModelIndex(), 0, nrow-1 );
//  peaks.swap( *m_peaks );
//  m_sortedPeaks.swap( sortedPeaks );
//  endInsertRows();
}//void sort( int column, Wt::SortOrder order )


bool PeakModel::compare( const PeakShrdPtr &lhs, const PeakShrdPtr &rhs,
                         Columns column, Wt::SortOrder order,
                         const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  const bool asscend = (order==AscendingOrder);
  
  if( !lhs || !rhs )
    return (asscend ? lhs.get()<rhs.get() : lhs.get()>rhs.get());

  switch( column )
  {
    case kMean:
      return (asscend ? lhs->mean()<rhs->mean() : lhs->mean()>rhs->mean());
    case kFwhm:
    {
      const double lw = lhs->gausPeak() ? lhs->sigma() : 0.5*lhs->roiWidth();
      const double rw = rhs->gausPeak() ? rhs->sigma() : 0.5*rhs->roiWidth();
      return (asscend ? lw<rw : lw>rw);
    }//case kFwhm:
      
    case kCps:  //We dont have live time, but if peaks are for same spectrum, it doesnt matter
    case kAmplitude:
      return (asscend ? lhs->peakArea()<rhs->peakArea() : lhs->peakArea()>rhs->peakArea());
      
    case kUserLabel:
      return (asscend ? lhs->userLabel()<rhs->userLabel() : lhs->userLabel()>rhs->userLabel());
    case kPeakLineColor:
      return (asscend ? lhs->lineColor().cssText()<rhs->lineColor().cssText() : lhs->lineColor().cssText()>rhs->lineColor().cssText());
    case kHasSkew:
      return (asscend ? lhs->skewType()<rhs->skewType() : lhs->skewType()>rhs->skewType());
    case kSkewAmount:
      return true;
    case kType:
      return (asscend ? lhs->gausPeak()<rhs->gausPeak() : lhs->gausPeak()>rhs->gausPeak());
    case kLowerX:
      return (asscend ? lhs->lowerX()<rhs->lowerX() : lhs->lowerX()>rhs->lowerX());
    case kUpperX:
      return (asscend ? lhs->upperX()<rhs->upperX() : lhs->upperX()>rhs->upperX());
    case kRoiCounts:
    {
      double lhs_area( 0.0 ), rhs_area( 0.0 );
      if( data )
      {
        rhs_area = data->gamma_integral( rhs->lowerX(), rhs->upperX() );
        lhs_area = data->gamma_integral( lhs->lowerX(), lhs->upperX() );
      }else
      {
        rhs_area = rhs->offset_integral( rhs->lowerX(), rhs->upperX() );
        lhs_area = lhs->offset_integral( lhs->lowerX(), lhs->upperX() );
      }
      
      return (asscend ? (lhs_area<rhs_area) : (lhs_area>rhs_area));
    }//case kRoiCounts:

    case kContinuumType:
      return (asscend ? lhs->continuum()->type()<rhs->continuum()->type() : lhs->continuum()->type()>rhs->continuum()->type());
    case kIsotope:
    {
      const SandiaDecay::Nuclide *lhsParent = lhs->parentNuclide();
      const SandiaDecay::Nuclide *rhsParent = rhs->parentNuclide();

      bool thisorder = false;
      if( !lhsParent )
        thisorder = false;
      else if( !rhsParent )
        thisorder = true;
      else
        thisorder = (lhsParent->symbol < rhsParent->symbol);
      
      return (asscend ? thisorder : (!thisorder));
    }//case kIsotope

    case kPhotoPeakEnergy:
    {
      float lhsEnergy = 0.0f, rhsEnergy = 0.0f;
      try
      {
        lhsEnergy = lhs->gammaParticleEnergy();
      }catch(std::exception &){}
      try
      {
        rhsEnergy = rhs->gammaParticleEnergy();
      }catch(std::exception &){}
      
      if( lhs == rhs )
        return false;
      
      const bool thisorder = (lhsEnergy < rhsEnergy);
      return (asscend ? thisorder : (!thisorder));
    }//case kPhotoPeakEnergy:

    case kDifference:
    {
      double lhsdiff = DBL_MAX, rhsdiff = DBL_MAX;
      
      try
      {
        lhsdiff = lhs->gammaParticleEnergy() - lhs->mean();
      }catch(std::exception &){}
      try
      {
        rhsdiff = rhs->gammaParticleEnergy() - rhs->mean();
      }catch(std::exception &){}
      
      const bool thisorder = (lhsdiff < rhsdiff);
      return (asscend ? thisorder : (!thisorder));
    }//case kDifference:
      
    case kUseForShieldingSourceFit:
    {
      const bool thisorder = lhs->useForShieldingSourceFit() < rhs->useForShieldingSourceFit();
      return (asscend ? thisorder : (!thisorder));
    }//case kUseForShieldingSourceFit:

    case kUseForCalibration:
    {
      const bool thisorder = lhs->useForCalibration() < rhs->useForCalibration();
      return (asscend ? thisorder : (!thisorder));
    }//case kUseForCalibration:

    case kCandidateIsotopes:
    case kNumColumns:
    break;
  }//switch( section )

  cerr << "PeakModel::compare(...): invalid sorting column" << endl;

  return false;
}//bool compare(...)



