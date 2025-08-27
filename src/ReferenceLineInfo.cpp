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
#include <mutex>
#include <tuple>
#include <vector>
#include <string>
#include <stdexcept>

#include <Wt/WLocale>
#include <Wt/WWebWidget>  //For quoting strings only

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/DateTime.h" //only for debug timing
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/Integrate.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/ReferenceLinePredef.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/PhysicalUnitsLocalized.h"


using namespace std;

/**
 version 0 and 1: serialized the ReferenceLineInfo struct.
 version 2: implemented 20221220 serializes just RefLineInput - you can then use this to generate ReferenceLineInfo
 */
const int RefLineInput::sm_xmlSerializationVersion = 2;


namespace
{
std::string jsQuote( const std::string &str )
{
  return Wt::WWebWidget::jsStringLiteral(str,'\"');
}


std::mutex sm_nuc_mix_mutex;
bool sm_have_tried_init = false;
/** Protected by `sm_nuc_mix_mutex`; is only read in once, and never freed or changed until end of program. */
std::shared_ptr<const std::map<std::string,ReferenceLinePredef::NucMix>> sm_nuc_mixes;



/** Also protected by `sm_nuc_mix_mutex`; is only created once, and never freed or changed until program termination. */
std::shared_ptr<const std::map<std::string,ReferenceLinePredef::CustomSrcLines>> sm_custom_lines;

/** RIght now we will only hold info about fission files in memory - we wont hold fission data in memory.
 Maybe once things are working fully, and useful, could it be useful to do this.
 */
struct FissionDataSrcFile
{
  const SandiaDecay::Nuclide *nuclide;
  std::string filepath;
};//struct FissionDataSrcFile

/** Also protected by `sm_nuc_mix_mutex`; is only created once, and never freed or changed until program termination. */
std::shared_ptr<std::vector<FissionDataSrcFile>> sm_fission_products;

void sanitize_label_str( string &label )
{
  SpecUtils::trim( label );
  SpecUtils::ireplace_all( label, " ", "" );
  SpecUtils::to_lower_ascii( label );
}//void sanitize_label_str( string &label )



  
void load_custom_nuc_mixes()
{
  // Takes about 2ms to run this function in Debug mode on M1 mac
  //const double start_time = SpecUtils::get_wall_time();
  
  typedef rapidxml::xml_node<char>      XmlNode;
  typedef rapidxml::xml_attribute<char> XmlAttribute;
  
  sm_have_tried_init = true;
  
  
  const string data_dir = InterSpec::staticDataDirectory();
  const string add_lines_path = SpecUtils::append_path( data_dir, "add_ref_line.xml" );
  
  auto nuc_mixes = make_shared<map<string,ReferenceLinePredef::NucMix>>();
  auto custom_lines = make_shared<map<string,ReferenceLinePredef::CustomSrcLines>>();
  
  ReferenceLinePredef::load_ref_line_file( add_lines_path, *nuc_mixes, *custom_lines, nullptr );
  
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  
  string user_data_dir;
  
  // `InterSpec::writableDataDirectory()` will throw if it hasnt been set.
  //  - e.g., in testing or batch mode
  try
  {
    user_data_dir = InterSpec::writableDataDirectory();
  }catch( std::exception & )
  {
  }//try / catch
  
  if( !user_data_dir.empty() )
  {
    const string custom_lines_path = SpecUtils::append_path( user_data_dir, "add_ref_line.xml" );
    
    // Any duplicate names, will overwrite what comes with InterSpec
    if( SpecUtils::is_file( custom_lines_path ) )
      ReferenceLinePredef::load_ref_line_file( custom_lines_path, *nuc_mixes, *custom_lines, nullptr );
  }//if( !user_data_dir.empty() )
#endif
  
  // Look for fission data
  auto fission_products = make_shared<vector<FissionDataSrcFile>>();
  try
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      throw runtime_error( "couldnt open DecayDataBaseServer" );
    
    const std::string fission_dir = SpecUtils::append_path( data_dir, "fission_yields" );
    const vector<string> data_files = SpecUtils::recursive_ls( fission_dir, "_independent_fy.csv" );
    for( const string &name_path : data_files )
    {
      string nuclide = SpecUtils::filename( name_path );
      auto underscore_pos = nuclide.find("_");
      
      assert( underscore_pos != std::string::npos );
      if( underscore_pos == std::string::npos )
        continue;
      nuclide = nuclide.substr( 0, underscore_pos );
      
      const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide );
      assert( nuc ); //Not a coding problem necassarily, but perhaps an eroneously named file?
      if( nuc )
      {
        FissionDataSrcFile d;
        d.nuclide = nuc;
        d.filepath = name_path;
        fission_products->push_back( std::move(d) );
      }else
      {
        cerr << "Fission file '" << name_path << "' doesnt appear to be for a nuclide." << endl;
      }
    }//for( const string &name_path : data_files )
    
    std::sort( begin(*fission_products), end(*fission_products),
              []( const FissionDataSrcFile &lhs, const FissionDataSrcFile &rhs ) {
      return SandiaDecay::Nuclide::lessThanForOrdering(lhs.nuclide, rhs.nuclide);
    });
  }catch( std::exception &e )
  {
    cerr << "Failed to interpret fission data file name: " << e.what() << endl;
  }
  
  sm_nuc_mixes = nuc_mixes;
  sm_custom_lines = custom_lines;
  sm_fission_products = fission_products;
  
  //const double end_time = SpecUtils::get_wall_time();
  //cout << "load_custom_nuc_mixes(): took " << (end_time - start_time) << " s" << endl;
}//void load_custom_nuc_mixes()


const ReferenceLinePredef::NucMix *get_custom_nuc_mix( std::string label )
{
  sanitize_label_str( label );
  
  std::lock_guard<std::mutex> lock( sm_nuc_mix_mutex );
  
  if( !sm_have_tried_init )
    load_custom_nuc_mixes();
  
  if( !sm_nuc_mixes )
    return nullptr;
  
  const auto pos = sm_nuc_mixes->find( label );
  if( pos == end(*sm_nuc_mixes) )
    return nullptr;
  
  return &(pos->second);
}//const ReferenceLinePredef::NucMix *get_custom_nuc_mix( std::string label )


const ReferenceLinePredef::CustomSrcLines *get_custom_src_lines( std::string label )
{
  sanitize_label_str( label );
  
  std::lock_guard<std::mutex> lock( sm_nuc_mix_mutex );
  
  if( !sm_have_tried_init )
    load_custom_nuc_mixes();
  
  if( !sm_custom_lines )
    return nullptr;
  
  const auto pos = sm_custom_lines->find( label );
  if( pos == end(*sm_custom_lines) )
    return nullptr;
  
  return &(pos->second);
}//const ReferenceLinePredef::CustomSrcLines *get_custom_src_lines( std::string label )


//The efficiency of S.E. and D.E. peaks, relative to F.E. peak, for the 20% Generic GADRAS DRF
//  included in InterSpec, is given pretty well by the following (energy in keV):
double ns_single_escape_sf( const double x )
{
  return std::max( 0.0, (1.8768E-11 *x*x*x) - (9.1467E-08 *x*x) + (2.1565E-04 *x) - 0.16367 );
};

double ns_double_escape_sf( const double x )
{
  return std::max( 0.0, (1.8575E-11 *x*x*x) - (9.0329E-08 *x*x) + (2.1302E-04 *x) - 0.16176 );
};



/** Struct to hold info from each line of `data/fission_yields/u235_independent_fy.csv` and simial */
struct NuclideYield
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  double thermal_yield = 0.0;
  double fast_yield = 0.0;
  double fourteen_MeV_yield = 0.0;
};//struct NuclideYield

shared_ptr<const vector<NuclideYield>> fission_nuclide_info( const SandiaDecay::Nuclide *nuclide )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Couldnt open decay database" );
  
  if( !nuclide )
    throw runtime_error( "No nuclide specified to get fission product info for." );
  
  static std::mutex s_cached_yields_mutex;
  static std::map<const SandiaDecay::Nuclide *,shared_ptr<const vector<NuclideYield>>> s_cached_yields;
  
  {// Begin lock on `s_cached_yields_mutex`
    std::lock_guard<std::mutex> lock( s_cached_yields_mutex );
    
    auto pos = s_cached_yields.find( nuclide );
    if( (pos != end(s_cached_yields)) && pos->second )
      return pos->second;
  }// End lock on `s_cached_yields_mutex`
  
  
  // Now find the filename for the datafile, from where we had cached this in memory
  string filename;
  {// begin lock on sm_nuc_mix_mutex
    std::lock_guard<std::mutex> lock( sm_nuc_mix_mutex );
    
    if( !sm_have_tried_init ) // If we havent cahced to memory yet, do it now
      load_custom_nuc_mixes();
    
    assert( sm_fission_products );
    
    if( !sm_fission_products )
      throw runtime_error( "Fission products not initialized" ); //shouldnt happen
    
    for( const FissionDataSrcFile &src : *sm_fission_products )
    {
      if( src.nuclide == nuclide )
        filename = src.filepath;
    }
  }// end lock on sm_nuc_mix_mutex
  
  if( filename.empty() )
    throw runtime_error( "Could not find fission product data for " + nuclide->symbol + "." );
  
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream input( filename.c_str() );
#else
  ifstream input( filename.c_str() );
#endif
  
  if( !input.is_open() )
    throw runtime_error( "Unable to open fission yield file '" + filename + "'" );
  
  auto fission_yields = make_shared<vector<NuclideYield>>();
  
  try
  {
    string line;
    int line_num = 0;
    int progeny_el_index = -1, progeny_A_index = -1, progeny_level_index = -1;
    int fy_thermal_index = -1, fy_fast_index = -1, fy_14MeV_index = -1;
    
    
    const double min_halflife = 60*SandiaDecay::second;
    
    while( SpecUtils::safe_get_line(input, line, 16384) )
    {
      if( line.size() > 16380 )
      {
        throw runtime_error( "Line " + std::to_string(line_num)
                            + " is longer than max allowed length of 16380 characters;"
                            " not reading in file." );
      }
      
      ++line_num;
      
      SpecUtils::trim( line );
      if( line.empty() || line[0]=='#' )
        continue;
      
      vector<string> fields;
      SpecUtils::split_no_delim_compress( fields, line, "," );
      
      if( fields.size() < 2 )
        continue;
      
      for( string &field : fields )
      {
        SpecUtils::trim( field );
        SpecUtils::to_lower_ascii( field );
        SpecUtils::ireplace_all( field, "_", " " );
      }
      
      if( progeny_el_index < 0 )
      {
        auto index_of_field = [&fields,&line]( string val ) -> int {
          auto pos = std::find( begin(fields), end(fields), val );
          if( pos == end(fields) )
            throw runtime_error( "Failed to find header value '" + val + "', in line '" + line + "'" );
          
          return static_cast<int>( pos - begin(fields) );
        };//index_of_field(...)
        
        progeny_el_index = index_of_field( "element daughter" );
        progeny_A_index = index_of_field( "a daughter" );
        progeny_level_index = index_of_field( "daughter level idx" );
        fy_thermal_index = index_of_field( "independent thermal fy" );
        fy_fast_index = index_of_field( "independent fast fy" );
        fy_14MeV_index = index_of_field( "independent 14mev fy" );
        
        continue;
      }//if( progeny_el_index < 0 )
      
      
      if( (fields.size() < progeny_el_index) || (fields.size() < progeny_A_index)
         || (fields.size() < progeny_level_index) || (fields.size() < fy_thermal_index)
         || (fields.size() < fy_fast_index) || (fields.size() < fy_14MeV_index) )
      {
        throw runtime_error( "Fewer fields (" + std::to_string(fields.size()) + ")"
                            " than expected on line " + std::to_string(line_num) );
      }
      
      string label = fields[progeny_el_index] + fields[progeny_A_index];
      int level = 0;
      if( !(stringstream(fields[progeny_level_index]) >> level) )
      {
        cout << "Failed to convert level '" << fields[progeny_level_index] << "' to int." << endl;
      }
      
      // TODO: figure out exactly what level corresponded to, in terms of meta-stable state.
      if( level > 0 )
      {
        label += "m";
        if( level > 1 )
          cout << label << " has level=" << level << endl;
      }//if( level > 0 )
      
      const SandiaDecay::Nuclide * const nuc = db->nuclide(label);
      if( !nuc )
      {
        cerr << "Failed to get nuc '" << label << "' from decay database - skipping." << endl;
        continue;
      }
      
      if( nuc->isStable() )
        continue;
      
      const vector<const SandiaDecay::Nuclide *> kids = nuc->descendants();
      double max_hl = nuc->halfLife;
      for( const SandiaDecay::Nuclide *kid : kids )
      {
        if( !kid->isStable() )
          max_hl = std::max( max_hl, kid->halfLife );
      }
      
      if( max_hl < min_halflife )
        continue;
      
      NuclideYield yield;
      yield.nuclide = nuc;
      
      auto get_yield = [&line_num]( const string &strval, double &val ){
        if( !strval.empty() && !SpecUtils::parse_double(strval.c_str(), strval.size(), val ) )
          throw runtime_error( "Failed to convert '" + strval
                              + "' to fission yield; line " + std::to_string(line_num) );
      };//get_yield lambda
      
      get_yield( fields[fy_thermal_index], yield.thermal_yield );
      get_yield( fields[fy_fast_index], yield.fast_yield );
      get_yield( fields[fy_14MeV_index], yield.fourteen_MeV_yield );
      
      const double max_yield = std::max( yield.fourteen_MeV_yield,
                                        std::max(yield.thermal_yield, yield.fast_yield) );
      if( max_yield <= 0.0 )
        continue;
      
      fission_yields->emplace_back( std::move(yield) );
    }//while( SpecUtils::safe_get_line(input, line) )
  }catch( std::exception &e )
  {
    throw runtime_error( "Failed to parse fission file '" + filename + "': " + string(e.what()) );
  }//try / catch to parse file
  
  std::sort( begin(*fission_yields), end(*fission_yields),
            []( const NuclideYield &lhs, const NuclideYield &rhs ) -> bool {
    return lhs.fast_yield > rhs.fast_yield;
  } );
  
  
  // Update our cache with these results.
  {
    std::lock_guard<std::mutex> lock( s_cached_yields_mutex );
    s_cached_yields[nuclide] = fission_yields;
  }
  
  return fission_yields;
}//vector<NuclideYield> fission_nuclide_info( const SandiaDecay::Nuclide *nuclide )


/** Holds information for a sinlge gamma/x-ray line, resulting from the fission decay chain. */
struct FissionLine
{
  const SandiaDecay::Nuclide *parent = nullptr;
  const SandiaDecay::Transition * transition = nullptr;
  SandiaDecay::ProductType particle_type;
  double energy;
  double relative_amplitude;
};//struct FissionLine


enum class FissionType
{
  Thermal,
  FourteenMeV,
  Fast
};

std::shared_ptr<const vector<FissionLine>> fission_photons( const SandiaDecay::Nuclide *nuc, const FissionType type,
                                                             const double irradiation_time_seconds, const double cool_off_time )
{
  const size_t max_cache_photons = 5; //arbitrary - lets speed things a little for the most recent uses, but not waste too much memory
  static std::mutex s_cache_photons_mutex;
  static std::deque<pair<tuple<const SandiaDecay::Nuclide *,int,double,double>,std::shared_ptr<const vector<FissionLine>>>> s_cache_photons;
  
  const std::tuple<const SandiaDecay::Nuclide *,int,double,double> index{nuc,
    static_cast<int>(type), irradiation_time_seconds, cool_off_time
  };
  
  
  {// Begin lock on `s_cache_photons_mutex` to see if we've already computed this
    std::lock_guard<std::mutex> lock( s_cache_photons_mutex );
    
    for( const auto &index_value : s_cache_photons )
    {
      if( index_value.first == index )
        return index_value.second;
    }
  }// End lock on `s_cache_photons_mutex`
  
  
  shared_ptr<const vector<NuclideYield>> fission_yields = fission_nuclide_info( nuc );
  assert( fission_yields );
  if( !fission_yields )
    throw std::logic_error( "fission_nuclide_info returned nullptr." );
  
  // We will have the max-buildup nuclide, buildup by ~1 uCi per second (about 86mCi per day,
  //  or 31 Ci per year)
  double num_atoms_mult_cs137 = std::numeric_limits<double>::max();
  double num_atoms_mult_all = std::numeric_limits<double>::max();
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *cs137 = db->nuclide( "Cs137" );
  assert( cs137 );
  
  size_t num_nucs_with_yield = 0.0;
  for( const NuclideYield &n : *fission_yields )
  {
    double yield = 0.0;
    switch( type )
    {
      case FissionType::Thermal:     yield = n.thermal_yield;      break;
      case FissionType::FourteenMeV: yield = n.fourteen_MeV_yield; break;
      case FissionType::Fast:        yield = n.fast_yield;         break;
    }
    
    assert( yield >= 0.0 );
    assert( !IsInf(yield) && !IsNan(yield) );
    
    if( (yield <= 0) || IsInf(yield) || IsNan(yield) )
      continue;
    
    if( n.nuclide && (n.nuclide == cs137) )
    {
      const double num_atoms_in_uCi = n.nuclide->activityToNumAtoms( 1.0E-6*SandiaDecay::Ci );
      num_atoms_mult_cs137 = num_atoms_in_uCi / yield;
    }
    
    if( n.nuclide->halfLife > 60*SandiaDecay::second )
    {
      const double num_atoms_in_uCi = n.nuclide->activityToNumAtoms( 1.0E-6*SandiaDecay::Ci );
      num_atoms_mult_all = std::min( num_atoms_mult_all, num_atoms_in_uCi / yield );
    }
    
    num_nucs_with_yield += 1;
  }//for( const NuclideYield &n : fission_yields )
  
  
  if( !num_nucs_with_yield )
  {
    string error_msg = nuc->symbol + " does not have any ";
    switch( type )
    {
      case FissionType::Thermal:     error_msg += "thermal"; break;
      case FissionType::FourteenMeV: error_msg += "14 MeV";  break;
      case FissionType::Fast:        error_msg += "fast";    break;
    }
    error_msg += " fission yields.";
    throw runtime_error( error_msg );
  }//if( !num_nucs_with_yield )
  
  
  assert( num_atoms_mult_all != std::numeric_limits<double>::max() );
  assert( num_atoms_mult_cs137 != std::numeric_limits<double>::max() );
  assert( num_atoms_mult_all > 0.0 );
  assert( num_atoms_mult_cs137 > 0.0 );
  if( (num_atoms_mult_cs137 == std::numeric_limits<double>::max())
     || (num_atoms_mult_cs137 <= 0.0)
     || IsInf(num_atoms_mult_cs137)
     || IsNan(num_atoms_mult_cs137) )
  {
    num_atoms_mult_cs137 = num_atoms_mult_all;
  }
  
  SandiaDecay::NuclideMixture ager;
  for( const auto &n : *fission_yields )
  {
    double yield = 0.0;
    switch( type )
    {
      case FissionType::Thermal:     yield = n.thermal_yield;      break;
      case FissionType::FourteenMeV: yield = n.fourteen_MeV_yield; break;
      case FissionType::Fast:        yield = n.fast_yield;         break;
    }
    
    ager.addAgedNuclideByNumAtoms(n.nuclide, yield*num_atoms_mult_cs137, 0.0 );
  }//for( const auto &n : *fission_yields )
  
  // Now we will crate a mixture that will represent, for its t=0, the end of
  //  build up, and after the looping over all the time steps, it will have
  //  all the information we want.
  SandiaDecay::NuclideMixture mixture;
  
  // This integration is very niave, and could be greatly improved.
  //  However, informal checks show for irradiation time of months, and common
  //  neutron activation products in metals, accuracy and numeric error didn't
  //  become notable issues (checked smaller time deltas, as well as using 128 bit,
  //  instead of 64 bit internal precisions in SandiaDecay, as well as exact
  //  expectations).
  auto do_integrate = [&mixture, &ager, irradiation_time_seconds]( double start_time,
                                                                  const double end_integration_time, const double dt ) {
    for( ; start_time < end_integration_time; start_time += dt )
    {
      const double end_time = std::min( start_time + dt, end_integration_time );
      const double this_dt = end_time - start_time;
      const double mid_time = 0.5*(start_time + end_time);
      const double time_until_irad_end = irradiation_time_seconds - mid_time;
      
      // Get the number of atoms, for all activation products, and their progeny, we expect
      //  at the end of buildup time.
      const vector<SandiaDecay::NuclideNumAtomsPair> num_atoms = ager.numAtoms( time_until_irad_end );
      for( size_t index = 0; index < num_atoms.size(); ++index )
      {
        const SandiaDecay::NuclideNumAtomsPair &nuc_num = num_atoms[index];
        mixture.addNuclideByAbundance( nuc_num.nuclide, (this_dt / irradiation_time_seconds)*nuc_num.numAtoms );
      }
    }//for( loop over buildup time )
  };//do_integrate lambda
  
  // Define how many seconds each time-step should be.
  //  A smaller time step should be more accurate, but you should take into account the
  //  half-lives of the nuclides you care about are.
  const int num_irad_steps = 200; //chosen arbitrarily
  const double time_delta = irradiation_time_seconds / num_irad_steps;
  
  // We do a very niave integration at the half-poin of each time integral - but if that
  //  half-timestep is compible to the half-lives we care about, as indicated by cool off time,
  //  then we can do a finer-stepped integration the last little bit of the integral, so the
  //  answer will be more accurate, especially for short-lived isotope
  //  TODO: this is only a very niave decision of what interval to use for finer integration - should revist
  double initial_end_time = irradiation_time_seconds;
  if( cool_off_time > 0.0 && cool_off_time < time_delta )
    initial_end_time -= cool_off_time;
  else if( cool_off_time < 5*time_delta )
    initial_end_time -= 5*time_delta;
  
  do_integrate( 0.0, initial_end_time, time_delta );
  if( initial_end_time != irradiation_time_seconds )
  {
    const double fine_time_delta = (irradiation_time_seconds - initial_end_time) / num_irad_steps;
    do_integrate( initial_end_time, irradiation_time_seconds, fine_time_delta );
  }
  
  
  /*
   vector<SandiaDecay::NuclideActivityPair> irrad_end_activities = mixture.activity( cool_off_time );
   std::sort( begin(irrad_end_activities), end(irrad_end_activities),
   []( const SandiaDecay::NuclideActivityPair &lhs, const SandiaDecay::NuclideActivityPair &rhs ) -> bool {
   return lhs.activity > rhs.activity;
   });
   
   cout << "At the end of irradiation, the activities are:\n";
   for( size_t index = 0; index < irrad_end_activities.size() && (index < 100); ++index )
   cout << "\t" << irrad_end_activities[index].nuclide->symbol
   << ": " << irrad_end_activities[index].activity << " bq" << endl;
   
   vector<SandiaDecay::NuclideActivityPair> after_cool_off_activities = mixture.activity( cool_off_time );
   std::sort( begin(after_cool_off_activities), end(after_cool_off_activities),
   []( const SandiaDecay::NuclideActivityPair &lhs, const SandiaDecay::NuclideActivityPair &rhs ) -> bool {
   return lhs.activity > rhs.activity;
   });
   
   // Get the Cs137 activity
   double cs137_act = 0.0;
   const auto cs137 = db->nuclide( "Cs137" );
   for( const SandiaDecay::NuclideActivityPair &nap : after_cool_off_activities )
   {
   if( nap.nuclide == cs137 )
   {
   cs137_act = nap.activity;
   break;
   }
   }
   assert( cs137_act != 0.0 );
   
   cout << "\n\nAfter cooling off for " << cool_off_time << " seconds the activities, relative to Cs137 are:\n";
   for( size_t index = 0; index < after_cool_off_activities.size() && (index < 100); ++index )
   {
   cout << "\t" << after_cool_off_activities[index].nuclide->symbol
   << ": " << after_cool_off_activities[index].activity/cs137_act << endl;
   }
   */
  
  // We expect A = A_0 * (1 - exp(-lamda * t_activation), so lets check things, but
  //  please note that this is only a valid check if no other activation products
  //  decay through the activation nuclide of interest.
  /*
   cout << endl << endl;
   for( size_t index = 0; index < irrad_end_activities.size() && (index < 100); ++index )
   {
   const SandiaDecay::Nuclide * const output_nuc = irrad_end_activities[index].nuclide;
   const double out_act = irrad_end_activities[index].activity;
   for( size_t input_index = 0; input_index < nuclides_rates.size(); ++input_index )
   {
   const SandiaDecay::Nuclide * const input_nuc = nuclides_rates[input_index].first;
   const double input_rate = nuclides_rates[input_index].second;
   if( input_nuc != output_nuc )
   continue;
   
   const double lambda = input_nuc->decayConstant();
   const double expected_act = input_rate * (1.0 - exp( -lambda * irradiation_time_seconds) );
   cout << "For " << input_nuc->symbol << " analytically expected " << expected_act
   << " bq; our calculation is " << out_act << " bq" << endl;
   }//for( loop over input nuclides )
   }//for( loop over output nuclides )
   */
  
  const vector<SandiaDecay::NuclideActivityPair> after_cool_off_activities
  = mixture.activity( cool_off_time );
  
  double max_rel_amp = 0.0;
  vector<FissionLine> all_lines;
  
  for( const auto &nap : after_cool_off_activities )
  {
    const SandiaDecay::Nuclide * const nuc = nap.nuclide;
    assert( nuc );
    for( const SandiaDecay::Transition * const trans : nuc->decaysToChildren )
    {
      for( const SandiaDecay::RadParticle &part : trans->products )
      {
        switch( part.type )
        {
          case SandiaDecay::BetaParticle:
          case SandiaDecay::AlphaParticle:
          case SandiaDecay::CaptureElectronParticle:
            break;
            
          case SandiaDecay::XrayParticle:
          case SandiaDecay::GammaParticle:
          case SandiaDecay::PositronParticle:
          {
            FissionLine line;
            line.parent = nuc;
            line.transition = trans;
            line.particle_type = part.type;
            line.energy = part.energy;
            line.relative_amplitude = part.intensity * trans->branchRatio * nap.activity;
            if( part.type == SandiaDecay::PositronParticle )
              line.relative_amplitude *= 2.0;
            
            max_rel_amp = std::max( max_rel_amp, line.relative_amplitude );
            all_lines.push_back( std::move(line) );
            break;
          }
        }//switch( part.type )
      }//for( const SandiaDecay::RadParticle &part : trans->products )
    }//for( const SandiaDecay::Transition * const trans : nuc->decaysToChildren )
  }//for( const auto &nap : after_cool_off_activities )
  
  // If lines are all pretty much the same energy, we will sum thier amplitude, and keep them all
  vector<pair<double,size_t>> amp_to_index;
  amp_to_index.reserve( all_lines.size() );
  for( size_t i = 0; i < all_lines.size(); /* i incrememted in inner loop */ )
  {
    double rel_amp = 0.0;
    const FissionLine &line = all_lines[i];
    
    size_t end_run = i;
    for( ; (end_run < all_lines.size())
        && (fabs(all_lines[end_run].energy - line.energy) < 0.01);
        ++end_run )
    {
      rel_amp += all_lines[end_run].relative_amplitude;
    }//for( size_t j = i + 1; j < all_lines.size(); ++j, ++i )
    
    assert( end_run > i );
    
    for( size_t j = i; j < end_run; ++j, ++i )
      amp_to_index.emplace_back( rel_amp, j );
  }//for( size_t i = 0; i < all_lines.size(); ++i )
  
  assert( amp_to_index.size() == all_lines.size() );
  
  // Now sort from biggest to smallest amplitude.
  std::sort( begin(amp_to_index), end(amp_to_index),
            []( const pair<double,size_t> &lhs, const pair<double,size_t> &rhs ) -> bool {
    return lhs.first > rhs.first;
  } );
  
  
  const size_t max_lines = 2200; //arbitrary
  auto keeper_lines = make_shared<vector<FissionLine>>();
  for( size_t i = 0; i < amp_to_index.size() && i < max_lines; ++i )
  {
    keeper_lines->push_back(all_lines[amp_to_index[i].second] );
  }//for( size_t i = 0; i < answer.size(); ++i )
  
  
  // Sort all lines by energy
  std::sort( begin(*keeper_lines), end(*keeper_lines),
            []( const FissionLine &lhs, const FissionLine &rhs ) -> bool {
    if( lhs.energy == rhs.energy ) // If same energy, put larger amplitude to left
      return rhs.relative_amplitude < lhs.relative_amplitude;
    return lhs.energy < rhs.energy;
  } );
  
  
  {// Begin lock on `s_cache_photons_mutex` to see if we've already computed this
    std::lock_guard<std::mutex> lock( s_cache_photons_mutex );
    s_cache_photons.emplace_front( index, keeper_lines );
    if( s_cache_photons.size() > max_cache_photons )
      s_cache_photons.resize( max_cache_photons );
  }// End lock on `s_cache_photons_mutex`
  
  
  return keeper_lines;
}//shared_ptr<vector<FissionLine>> fission_photons(...)
}//namespace


/** Computes background lines by aging K40, Th232, U235, U238, and Ra226, and then transporting,
 as a trace source, through a 1m radius soil sphere.
 
 Since this computation takes ~0.15 seconds, the results are cached and returned on subsequent calls.
 
 
 TODO: Return #ReferenceLineInfo::RefLine instead, and completely get rid of #OtherRefLine, since it is no longer needed.
 */
const vector<OtherRefLine> &getBackgroundRefLines()
{
  using namespace GammaInteractionCalc;
  
  const float lower_photon_energy = 10.0f;
  
  static std::mutex answer_mutex;
  static bool have_computed = false;
  static vector<OtherRefLine> answer;
  
  std::lock_guard<std::mutex> lock( answer_mutex );
  if( have_computed )
    return answer;
  
  have_computed = true;
  
  try
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      return answer;
    
    // The threshold, relative to most intense line, to not include intensities below.
    //  1.0E-17 gives 2064 lines
    //  1.0E-16 gives 1927 lines
    //  1.0E-15 gives 1851 lines
    const double rel_threshold = 1.0E-17;
    
    const char *soil_chem_formula = "H0.022019C0.009009O0.593577Al0.066067Si0.272289K0.01001Fe0.027029 d=1.6";
    
    const Material soilobj = MaterialDB::materialFromChemicalFormula( soil_chem_formula, db );
    const Material * const soil = &soilobj;
    
    assert( soilobj.elements.size() == 7 );
    assert( fabs(soilobj.density - 1.6*PhysicalUnits::g/PhysicalUnits::cm3) < 0.001*soilobj.density );
    
    vector<OtherRefLine> prelim_answer;
    prelim_answer.resize( 3000 ); //we actually need 2890
    answer.reserve( 2100 ); //we actaully need 2062, when rel_threshold==1.0E-17;
    
    
// Pre-computing transport just prints out `integral_energies` and `integral_values` to stdout
//  Only used for development purposes to get the aformentioned arrays.
#define PRE_COMPUTE_BACKGROUND_TRANSPORT 0
    
#define USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT 1
    
    //const double start_wall = SpecUtils::get_wall_time();
    //const double start_cpu = SpecUtils::get_cpu_time();
    
#if( PRE_COMPUTE_BACKGROUND_TRANSPORT || !USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
    // Computation timings on M1 macbook, for integrating each background line:
    // - Single threaded, epsrel = 1e-5: wall=0.789006, cpu=0.787189
    // - Single threaded, epsrel = 1e-4: wall=0.593146, cpu=0.59266
    // - Single threaded, epsrel = 1e-3: wall=0.465307, cpu=0.464454
    // - Multithreaded (nthreads=10, batching in groups of 10), epsrel = 1e-4: wall=0.131257, cpu=0.669642
    // - Multithreaded (nthreads=10, batching in groups of 100), epsrel = 1e-4: wall=0.087524, cpu=0.72624
    // - Multithreaded (nthreads=10, no batching), epsrel = 1e-4: wall=0.09096, cpu=0.710744
    //
    // - Using precomputed transport, and single thread: wall=0.010582, cpu=0.009713
    
    DistributedSrcCalc soil_sphere;
    soil_sphere.m_geometry = GeometryType::Spherical;
    soil_sphere.m_materialIndex = 0;
    soil_sphere.m_attenuateForAir = false;
    soil_sphere.m_airTransLenCoef = 0.0;
    soil_sphere.m_isInSituExponential = false;
    soil_sphere.m_inSituRelaxationLength = 0.0;
    soil_sphere.m_detectorRadius  = 5.0 * PhysicalUnits::cm;
    soil_sphere.m_observationDist = 200.0 * PhysicalUnits::cm;
    
    const double sphereRad = 100.0 * PhysicalUnits::cm;
    soil_sphere.m_dimensionsTransLenAndType.push_back( {
      {sphereRad,0.0,0.0},
      0.0,
      DistributedSrcCalc::ShellType::Material
    } );
#endif // !USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT
  
#if( USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
    const array<double, 220> integral_energies{
           0,        5,        6,        7,        8,        9,       10,       11,       12,       13,       14,       15,       16,       17,
          18,       19,       20,       21,       22,       23,       24,       25,       26,       27,       28,       29,       30,       31,       32,
          33,       34,       35,       36,       37,       38,       39,       40,       41,       42,       43,       44,       45,       46,       47,
          48,       49,       50,       51,       52,       53,       54,       55,       56,       57,       58,       59,       60,       61,       62,
          63,       64,       65,       66,       67,       68,       69,       70,       71,       72,       73,       74,       75,       77,       79,
          81,       83,       85,       87,       89,       91,       93,       95,       97,       99,      102,      105,      108,      111,      114,
         117,      120,      124,      128,      132,      136,      140,      144,      149,      154,      159,      164,      169,      174,      180,
         186,      192,      198,      205,      215,      225,      235,      245,      255,      265,      275,      285,      295,      305,      315,
         325,      335,      345,      355,      365,      375,      385,      400,      410,      425,      440,      455,      470,      485,      500,
         515,      530,      545,      560,      575,      590,      605,      620,      635,      650,      665,      680,      700,      720,      740,
         760,      780,      800,      820,      840,      860,      880,      900,      920,      940,      965,      990,     1015,     1040,     1065,
        1090,     1115,     1140,     1165,     1190,     1215,     1240,     1270,     1295,     1325,     1355,     1385,     1415,     1445,     1475,
        1505,     1540,     1575,     1610,     1645,     1680,     1715,     1750,     1785,     1825,     1865,     1905,     1945,     1985,     2025,
        2065,     2110,     2155,     2200,     2245,     2290,     2340,     2390,     2440,     2490,     2540,     2595,     2650,     2705,     2760,
        2820,     2880,     2940,     3005,     3070,     3135,     3200,     3270,     3340,     3415,     3490};
    
    const array<double, 220> integral_values{
            0, 14.134587, 37.955691, 61.053456, 81.553097, 106.44044, 146.33149, 195.13829, 253.84184, 323.29028, 404.16482, 497.45935, 603.72919, 723.53085,
    857.62272, 1006.3214, 1169.7335, 1348.316, 1542.3132, 1751.454, 1975.4162, 2213.6208, 2466.3235, 2733.6567, 3014.2296, 3306.6471, 3610.0513, 3923.602, 4246.7122,
    4579.8723, 4921.6507, 5269.3928, 5621.8764, 5977.892, 6336.2679, 6696.0448, 7057.1923, 7420.2084, 7783.9253, 8146.2877, 8505.9111, 8861.7795, 9213.1772, 9559.4574,
    9900.0845, 10234.865, 10563.582, 10888.434, 11209.571, 11524.162, 11831.569, 12131.616, 12424.397, 12710.119, 12988.515, 13259.309, 13522.553, 13778.266, 14026.609,
    14267.881, 14505.645, 14740.374, 14968.793, 15190.695, 15406.222, 15615.472, 15818.601, 16015.768, 16207.127, 16392.842, 16573.098, 16748.05, 16999.027, 17319.145,
    17627.432, 17927.34, 18213.695,  18485.5, 18743.677, 18989.13, 19222.717, 19445.229, 19657.408, 19859.946, 20106.772, 20397.573, 20675.541, 20938.468, 21187.707,
    21424.455, 21649.762, 21898.576, 22175.379, 22450.301, 22717.998, 22973.563, 23218.067, 23480.578, 23759.766, 24027.788, 24297.555, 24568.681, 24830.316, 25107.653,
    25399.905,    25682, 25954.813, 26251.616, 26642.112, 27088.381, 27516.833, 27929.463, 28333.777, 28740.83, 29145.616, 29539.17, 29922.264, 30295.645, 30660.013,
    31027.158, 31399.069, 31764.893, 32123.481, 32475.234,  32820.5, 33159.572, 33576.174, 33995.696,  34418.2, 34915.801, 35402.684, 35879.466, 36346.686, 36804.84,
    37266.512, 37733.084, 38192.789, 38644.971, 39089.931, 39527.976, 39959.416, 40384.528, 40806.642, 41234.388, 41664.859, 42089.819, 42578.421, 43128.911, 43670.708,
    44204.202, 44729.725, 45251.585, 45780.138, 46311.632, 46836.248, 47354.257, 47865.884, 48371.341, 48870.845, 49425.509, 50034.01, 50643.808, 51261.626, 51877.974,
    52486.515, 53087.301, 53680.377, 54265.847, 54843.807, 55414.404, 55977.748, 56595.304, 57219.765, 57846.156, 58520.698, 59185.681, 59841.292, 60487.865, 61125.588,
    61754.475, 62425.916, 63138.406, 63852.858, 64575.719, 65293.349, 65999.372, 66695.524, 67382.016, 68106.916, 68868.824, 69617.361, 70354.821, 71080.753, 71808.509,
    72545.04, 73321.576, 74129.992, 74926.528, 75712.654, 76485.167, 77284.546, 78111.513, 78928.337, 79733.538, 80534.094, 81378.643, 82254.629, 83120.668, 83977.388,
    84854.37, 85750.145, 86625.844, 87517.263, 88433.815, 89339.424, 90237.042, 91169.853, 92123.699, 93086.383, 94066.595
    };

    /* const array<double, 220> integral_differences{
           0,        1,  0.40664,   0.3602,  0.15993,  0.28653,  0.26234,  0.24083,  0.22383,  0.20768,    0.194,  0.18226,  0.17085,  0.16116,
     0.15228,   0.1439,  0.13607,  0.12929,   0.1227,   0.1165,   0.1106,  0.10494,  0.10024, 0.095584, 0.090808, 0.086263, 0.082007, 0.077984, 0.074327,
    0.071275, 0.067737, 0.064361, 0.061138, 0.058065, 0.055137, 0.052396, 0.050014, 0.047884, 0.045622, 0.043391, 0.041213, 0.039143, 0.037176, 0.035306,
    0.033537, 0.031909, 0.030351, 0.029333, 0.027983, 0.026632, 0.025348, 0.024132, 0.023011,  0.02196, 0.020919, 0.019937, 0.019006,  0.01812, 0.017297,
     0.01653, 0.016255, 0.015599, 0.014926, 0.014294, 0.013689, 0.013115, 0.012571, 0.012054, 0.011563, 0.011098, 0.010657, 0.010237, 0.019204, 0.017776,
    0.017207, 0.016259, 0.015194, 0.014221, 0.013333, 0.012524, 0.011784, 0.011106, 0.010485, 0.0099147, 0.014602, 0.013916, 0.012979, 0.012141, 0.011391,
    0.010714, 0.010103, 0.012605, 0.012361, 0.012132, 0.011439, 0.010813, 0.010251, 0.012097, 0.011408, 0.010904, 0.011299, 0.010775, 0.010301, 0.011782,
    0.011233, 0.010738, 0.010287, 0.012313, 0.016961, 0.015995, 0.015152, 0.014401,  0.01414, 0.014186, 0.013595, 0.013054, 0.012555, 0.012097, 0.011673,
    0.011991,   0.0117, 0.011335, 0.010993, 0.010672, 0.010369, 0.010083, 0.014698, 0.010006, 0.014512, 0.013995, 0.013514, 0.013066, 0.012646, 0.012253,
    0.012522,  0.01221, 0.011865, 0.011539, 0.011229, 0.010936, 0.010659, 0.010396, 0.010293, 0.010453, 0.010212, 0.0099823, 0.012949, 0.012581, 0.012234,
    0.011906, 0.011594, 0.011472, 0.011618, 0.011336, 0.011067, 0.010812, 0.010567, 0.010333,  0.01011, 0.012321, 0.012004, 0.012077, 0.012028, 0.011736,
    0.011454, 0.011181, 0.010917, 0.010662, 0.010416, 0.010179, 0.0099493, 0.011863, 0.0099733, 0.011674, 0.011381, 0.011092, 0.010821, 0.010558, 0.010309,
     0.01006, 0.011444, 0.011127,  0.01125, 0.011138, 0.010845, 0.010551, 0.010325, 0.010052, 0.011228,   0.0109, 0.010606, 0.010359, 0.010068, 0.010201,
    0.010105, 0.011071, 0.010741, 0.010522, 0.010246, 0.009956, 0.010727, 0.010449,  0.01025, 0.0099488, 0.0099325, 0.010819, 0.010482, 0.010356, 0.010049,
    0.010619, 0.010276, 0.0099437, 0.010425, 0.010304, 0.009971, 0.0099239, 0.010536, 0.010174, 0.010508, 0.010334
    };
    */
#endif

    
#if( PRE_COMPUTE_BACKGROUND_TRANSPORT )
#warning "You probably dont mean to have PRE_COMPUTE_BACKGROUND_TRANSPORT enabled"
    
    // Pre-compute the transport through soil sphere, defining energy bounds so that the
    //  difference between lower and upper energy bounds, is 1% (arbitrarily chosen);
    //  For energies, below about 75 keV our error is larger due to stepping size
    const double allowed_error_fraction = 0.01;
    double prev_integral = 0.0, prev_energy = 0.0;
    vector<double> energies(1,0.0), integrals(1,0.0), errors(1,0.0);
    for( double energy = 5; energy < 3500; energy += (energy < 200 ? 1 : 5) )
    {
      DistributedSrcCalc sphere = soil_sphere;
      double transLenCoef = GammaInteractionCalc::transmition_length_coefficient( soil, energy );
      get<1>( sphere.m_dimensionsTransLenAndType[0] ) = transLenCoef;
      
      int nregions, neval, fail;
      double integral, error, prob;
      void *userdata = (void *)&sphere;
      
      const int ndim = 2;  //the number of dimensions of the integral.
      const double epsrel = 1e-8, epsabs = -1.0; //the requested relative and absolute accuracies
      const int mineval = 0, maxeval = 5000000;   //the min and (approx) max number of integrand evaluations allowed.
      
      
      Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_spherical, userdata, epsrel, epsabs,
                                Integrate::LastImportanceFcnt, mineval, maxeval, nregions, neval,
                                fail, integral, error, prob );
      
      if( fabs(integral - prev_integral) > allowed_error_fraction*prev_integral )
      {
        energies.push_back( energy );
        integrals.push_back( 0.5*(prev_integral + integral) );
        errors.push_back( fabs(integral - prev_integral)/integral );
        //cout << interval_num << ": energy=" << energy << "] keV" << ", diff="
        //     << (100*fabs(integral - prev_integral)/integral) << endl;
        prev_energy = energy;
        prev_integral = integral;
      }//if( integral at `energy` is more than `allowed_error_fraction` different from previous )
    }//for( loop over energies )
    
    cout << "const std::array<double, " << energies.size() << "> integral_energies{\n";
    for( size_t i = 0; i < energies.size(); ++i )
      cout << (i ? ", " : " ") << ((((i+1) % 15) == 0) ? "\n" : "")
           << std::setw(8) << std::setprecision(6) << energies[i];
    cout << "};" << endl;
    
    cout << "const std::array<double, " << energies.size() << "> integral_values{\n";
    for( size_t i = 0; i < integrals.size(); ++i )
      cout << (i ? ", " : " ") << ((((i+1) % 15) == 0) ? "\n" : "")
           << std::setw(8) << std::setprecision(8) << integrals[i];
    cout << "\n};\n" << endl;
    
    cout << "/* const std::array<double, " << energies.size() << "> integral_differences{\n";
    for( size_t i = 0; i < errors.size(); ++i )
      cout << (i ? ", " : " ") << ((((i+1) % 15) == 0) ? "\n" : "")
           << std::setw(8) << std::setprecision(5) << errors[i];
    cout << "\n};\n*/" << endl;
#endif //PRE_COMPUTE_BACKGROUND_TRANSPORT
    
    
    // The rel-activities below were just adjusted to match a representative background
    //  and are otherwise arbitrary.
    const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
    const SandiaDecay::Nuclide * const ra226 = db->nuclide( "Ra226" );
    const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
    const SandiaDecay::Nuclide * const th232 = db->nuclide( "Th232" );
    const SandiaDecay::Nuclide * const k40 = db->nuclide( "K40" );
    
    assert( u238 && ra226 && u235 && th232 && k40 );
    
    // All the denominators below were for observation distance of 200cm, 5cm detector radius, and
    //  a 1 m radius sphere
    // TODO: figure out the constants so components can be set as PPM, or % of soil, so it can be changed by looking up in a DB easily for a given location
    const vector<tuple<const SandiaDecay::Nuclide *,double,OtherRefLineType, double>> nuc_activity{
      //make the 1001 keV have amp 0.0004653, before norm; the denom is integral at 1001 keV times BR of 1001.
      { u238,   0.0004653/410.2892, OtherRefLineType::U238Series, 5.0*u238->promptEquilibriumHalfLife() },
      
      //make 609 keV have amp 0.02515, before norm
      { ra226,  0.02515/17990.5430, OtherRefLineType::Ra226Series, 5.0*ra226->promptEquilibriumHalfLife() },
      
      //make 185 keV have amp 0.001482, before norm
      { u235,   0.001482/14603.0156, OtherRefLineType::U235Series, 5.0*u235->promptEquilibriumHalfLife() },
      
      //make 2614 keV have amp 0.02038, before norm
      { th232,  0.02038/27897.2617, OtherRefLineType::Th232Series, 5.0*th232->secularEquilibriumHalfLife() },
      
      //make 1460 keV have amp 0.1066, before norm
      { k40,    0.1066/6523.8994, OtherRefLineType::K40Background, 0.0 }
    };//nuc_activity
    
    
    SandiaDecay::NuclideMixture mixture;
    for( const auto &src : nuc_activity )
    {
      const SandiaDecay::Nuclide * const nuc = get<0>(src);
      const double parent_activity = get<1>(src);
      const double age = get<3>(src);
      
      mixture.addAgedNuclideByActivity( nuc, parent_activity, age );
    }//for( const auto &src : nuc_activity )
    
    
    const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( 0.0 );
    
    size_t calc_index = 0;
#if( !USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
    SpecUtilsAsync::ThreadPool pool;
#endif
    
    for( const SandiaDecay::NuclideActivityPair &nap : activities )
    {
      const SandiaDecay::Nuclide *nuclide = nap.nuclide;
      const double activity = nap.activity;
      
      for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
      {
        for( const SandiaDecay::RadParticle &particle : transition->products )
        {
          // For the moment we'll keep gammas and x-rays only - but should
          switch( particle.type )
          {
            case SandiaDecay::BetaParticle:
            case SandiaDecay::AlphaParticle:
            case SandiaDecay::CaptureElectronParticle:
            case SandiaDecay::PositronParticle: //Posititrons are super-small, and not worth adding in here
              continue;
              break;
              
            case SandiaDecay::GammaParticle:
            case SandiaDecay::XrayParticle:
              if( particle.energy < lower_photon_energy )
                continue;
              break;
          };//switch( particle.type )
          
          
          const double br = activity * particle.intensity * transition->branchRatio;
          const double energy = particle.energy;
          const SandiaDecay::ProductType part_type = particle.type;
          
#if( !USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
          // Creating a `DistributedSrcCalc` here doesnt seem any slower than trying to share them
          //  between threads and such; so we'll just make a copy for each thread.
          DistributedSrcCalc sphere = soil_sphere;
          const double transLenCoef = GammaInteractionCalc::transmition_length_coefficient( soil, energy );
          get<1>( sphere.m_dimensionsTransLenAndType[0] ) = transLenCoef;
#endif
          
          auto do_calc = [soil, energy, calc_index, transition, part_type, br,
#if( USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
                          &integral_energies, &integral_values,
#else
                          sphere,
#endif
                          k40, ra226, th232, u238, u235, &prelim_answer]() {
            
#if( USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
            auto pos = std::lower_bound( begin(integral_energies), end(integral_energies), energy );
            const auto index = pos - begin(integral_energies);
            const double integral = (pos == end(integral_energies)) ? integral_values.back() : integral_values[index];

//            cout << "Energy=" << energy << " is in bucket ["
//            << (index ? integral_energies[index-1] : integral_energies[index]) << " to "
//            << (index >= integral_energies.size() ? integral_energies[index-1] : integral_energies[index])
//            << "] keV, with value integral=" << integral << endl;
#else
            int nregions, neval, fail;
            double integral, error, prob;
            void *userdata = (void *)&sphere;
            
            // Some constants for integration
            const int ndim = 2;  //the number of dimensions of the integral.
            const double epsrel = 1e-4, epsabs = -1.0; //the requested relative and absolute accuracies
            const int mineval = 0, maxeval = 500000;   //the min and (approx) max number of integrand evaluations allowed.
            
            Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_spherical, userdata, epsrel, epsabs,
                                      Integrate::LastImportanceFcnt, mineval, maxeval, nregions, neval,
                                      fail, integral, error, prob );
            //printf("%s: %.1f keV -> br=%.4f.\n", nuc->symbol.c_str(), energy, br*integral );
#endif // USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT / else
            
            string desc;
            if( transition->parent )
              desc = transition->parent->symbol;
            if( part_type == SandiaDecay::XrayParticle )
              desc += " x-ray";
            
            OtherRefLineType type = OtherRefLineType::BackgroundXRay;
            if( part_type == SandiaDecay::XrayParticle )
            {
              type = OtherRefLineType::BackgroundXRay;
            }else if( transition->parent )
            {
              if( k40->branchRatioToDecendant(transition->parent) > 0 )
                type = OtherRefLineType::K40Background;
              else if( ra226->branchRatioToDecendant(transition->parent) > 0 )
                type = OtherRefLineType::Ra226Series;
              else if( th232->branchRatioToDecendant(transition->parent) > 0 )
                type = OtherRefLineType::Th232Series;
              else if( u238->branchRatioToDecendant(transition->parent) > 0 )
                type = OtherRefLineType::U238Series;
              else if( u235->branchRatioToDecendant(transition->parent) > 0 )
                type = OtherRefLineType::U235Series;
              else{ assert( 0 ); }
            }
            
            const float amplitude = static_cast<float>( br * integral );
             
            prelim_answer[calc_index] = { static_cast<float>(energy), amplitude, desc, type, "" };
          };//do_calc
                    
#if( USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
          do_calc();
#else
          // We dont want to bog the thread pool down with ~2k jobs, so we'll batch things in
          //  groups of an arbitrarily selected number of 100
          const size_t batch_size = 100;
          if( (calc_index % batch_size) == 0 )
          {
            pool.join();
          
            assert( prelim_answer.size() > calc_index );
            if( prelim_answer.size() < (calc_index + batch_size) )
              prelim_answer.resize( calc_index + batch_size );
          }//if( (calc_index % 100) == 0 )
          
          pool.post( do_calc );
#endif //USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT / else
          
          calc_index += 1;
        }//for( const SandiaDecay::RadParticle &particle : transition->products )
      }//for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
    }//for( const SandiaDecay::NuclideActivityPair &nap : activities )
    
#if( !USE_PRE_COMPUTED_BACKGROUND_LINE_TRANSPORT )
    pool.join();
#endif
    
    assert( prelim_answer.size() >= calc_index );
    prelim_answer.resize( calc_index );
    
    //const double end_wall = SpecUtils::get_wall_time();
    //const double end_cpu = SpecUtils::get_cpu_time();
    
    //cout << "Background line computation took: wall=" << (end_wall - start_wall) << ", cpu=" << (end_cpu - start_cpu) << endl;
    //cout << "prelim_answer.size=" << prelim_answer.size() << endl;
    
    float max_br = 0.0f;
    for( const auto &i : prelim_answer )
    {
      //printf("prenorm %s: %.1f keV -> br=%.4f.\n", get<2>(i).c_str(), get<0>(i), get<1>(i) );
      max_br = std::max( max_br, get<1>(i) );
    }
    
    // Lets grab the amplitude of the 2614.5330 keV Th232 line
    double th232_2614_amp = 0;
    
    for( auto &i : prelim_answer )
    {
      get<1>(i) /= max_br;
      if( get<1>(i) > rel_threshold )
        answer.push_back( i );
      
      if( (get<0>(i) > 2614.52)
         && (get<0>(i) < 2614.55)
         && (get<3>(i) == OtherRefLineType::Th232Series) )
      {
        th232_2614_amp += get<1>(i);
      }
      
      // printf("%s: %.1f keV -> br=%.4f.\n", get<2>(i).c_str(), get<0>(i), get<1>(i) );
    }//for( auto &i : prelim_answer )
    
    // th232_2614_amp should be 0.19118185341358185 );
    assert( (th232_2614_amp > 0.05) && (th232_2614_amp < 0.5) );
    
    // Use S.E. and D.E. numbers from a 40% HPGe detector
    answer.emplace_back( 2614.533f - 510.9989f,      0.144*th232_2614_amp,
                         "Th232 S.E. 2614 keV", OtherRefLineType::OtherBackground, "" );
    answer.emplace_back( 2614.533f - 2.0f*510.9989f, 0.082*th232_2614_amp,
                        "Th232 D.E. 2614 keV", OtherRefLineType::OtherBackground, "" );
    
    // Now take care of annihilation; for a single background spectrum I looked at, it was about 5%
    //  the amplitude of the K40 line; didnt correct for DRF
    answer.emplace_back( 510.9989f, 0.052f, "",
                        OtherRefLineType::OtherBackground, "Annihilation radiation (beta+)" );
    
    // TODO: Do we want to add in extra x-rays or anything?
    // TODO: Do we group lines together later on, or should we do it now.
    
    std::sort( begin(answer), end(answer), []( const OtherRefLine &lhs, const OtherRefLine &rhs ) -> bool {
      return get<0>(lhs) < get<0>(rhs);
    });
    
    //cout << "answer.size=" << answer.size() << endl;
    //size_t nbelow = 0, nabove = 0;
    //for( auto &i : answer )
    //  (get<1>(i) > rel_threshold ? nabove : nbelow) += 1;
    //cout << "Got nabove=" << nabove << ", nbelow=" << nbelow << endl;
    //for 1.-E17, get: "Got nabove=2062, nbelow=829"
  }catch( std::exception &e )
  {
    assert( 0 );
    answer.clear();
  }//try / catc
  
  return answer;
}//vector<OtherRefLine> getBackgroundRefLines()


const char *to_str( const OtherRefLineType type )
{
  switch( type )
  {
    case OtherRefLineType::U238Series: return "U238Series";
    case OtherRefLineType::U235Series: return "U235Series";
    case OtherRefLineType::Th232Series: return "Th232Series";
    case OtherRefLineType::Ra226Series: return "Ra226Series";
    case OtherRefLineType::K40Background: return "K40Background";
    case OtherRefLineType::BackgroundXRay: return "BackgroundXRay";
    case OtherRefLineType::BackgroundReaction: return "BackgroundReaction";
    case OtherRefLineType::OtherBackground: return "Other";
  }//switch( type )

  return "InvalidOtherRefLineType";
}//to_str(...)


OtherRefLineType other_ref_line_type_from_str( const std::string &str )
{
  for( auto type = OtherRefLineType(0);
    type != OtherRefLineType::OtherBackground;
    type = OtherRefLineType( static_cast<int>(type) + 1) )
  {
    if( str == to_str( type ) )
      return type;
  }

  return OtherRefLineType::OtherBackground;
}



RefLineInput::RefLineInput()
  : m_input_txt(),
  m_age(),
  m_color(),
  m_lower_br_cutt_off( 0.0 ),
  m_promptLinesOnly( false ),
  m_showGammas( false ),
  m_showXrays( false ),
  m_showAlphas( false ),
  m_showBetas( false ),
  m_showCascades( false ),
  m_showEscapes( false ),
  m_detector_name(),
  m_det_intrinsic_eff(),
  m_shielding_name(),
  m_shielding_thickness(),
  m_shielding_an(),
  m_shielding_ad(),
  m_shielding_att()
{
}

void RefLineInput::reset()
{
  *this = RefLineInput{};
}


bool RefLineInput::operator==(const RefLineInput &rhs) const
{
  return ((m_input_txt == rhs.m_input_txt)
  && (m_age == rhs.m_age)
  && (m_color == rhs.m_color)
  && (m_lower_br_cutt_off == rhs.m_lower_br_cutt_off)
  && (m_promptLinesOnly == rhs.m_promptLinesOnly)
  && (m_showGammas == rhs.m_showGammas)
  && (m_showXrays == rhs.m_showXrays)
  && (m_showAlphas == rhs.m_showAlphas)
  && (m_showBetas == rhs.m_showBetas)
  && (m_showCascades == rhs.m_showCascades)
  && (m_showEscapes == rhs.m_showEscapes)
  && (m_detector_name == rhs.m_detector_name)
  //&& (m_det_intrinsic_eff == rhs.)
  && (m_shielding_name == rhs.m_shielding_name)
  && (m_shielding_thickness == rhs.m_shielding_thickness)
  && (m_shielding_an == rhs.m_shielding_an)
  && (m_shielding_ad == rhs.m_shielding_ad)
  //&& (m_shielding_att == rhs.m_shielding_att)
  );
}//bool operator==(const RefLineInput &rhs) const


ReferenceLineInfo::RefLine::RefLine()
  : m_energy( 0.0 ),
  m_normalized_intensity( 0.0 ),
  m_drf_factor( 1.0f ),
  m_shield_atten( 1.0f ),
  m_particle_sf_applied( 1.0f ),
  m_decaystr(),
  m_color(),
  m_decay_intensity( 0.0 ),
  m_particle_type( ReferenceLineInfo::RefLine::Particle::Gamma ),
  m_parent_nuclide( nullptr ),
  m_transition( nullptr ),
  m_source_type( ReferenceLineInfo::RefLine::RefGammaType::Normal ),
  m_attenuation_applies( true ),
  m_element( nullptr ),
  m_reaction( nullptr )
{
}


const std::string &ReferenceLineInfo::RefLine::particlestr() const
{
  static const string cascade_sum = "cascade-sum";
  static const string gamma_sum = "sum-gamma";
  static const string single_escape = "S.E.";
  static const string double_escape = "D.E.";
  
  static const string alpha = "alpha";
  static const string beta = "beta";
  static const string gamma = "gamma";
  static const string xray = "xray";
  static const string other = "invalid";
  
  
  switch( m_source_type )
  {
    case RefGammaType::Normal:
    case RefGammaType::Annihilation:
      break;
      
    case RefGammaType::SingleEscape:        return single_escape;
    case RefGammaType::DoubleEscape:        return double_escape;
    case RefGammaType::CoincidenceSumPeak:  return cascade_sum;
    case RefGammaType::SumGammaPeak:        return gamma_sum;
  }//switch( m_source_type )
  
  
  switch( m_particle_type )
  {
    case ReferenceLineInfo::RefLine::Particle::Alpha: return alpha;
    case ReferenceLineInfo::RefLine::Particle::Beta:  return beta;
    case ReferenceLineInfo::RefLine::Particle::Gamma: return gamma;
    case ReferenceLineInfo::RefLine::Particle::Xray:  return xray;
  }//switch( line.m_particle_type )
  
  assert( 0 );
  return other;
}//particlestr(...)

ReferenceLineInfo::ReferenceLineInfo()
{
  reset();
}

void ReferenceLineInfo::reset()
{
  m_ref_lines.clear();
  m_input_warnings.clear();
  m_validity = InputValidity::Blank;
  m_has_coincidences = false;
  m_input = RefLineInput();
  m_source_type = ReferenceLineInfo::SourceType::None;

  m_nuclide = nullptr;
  m_element = nullptr;
  m_reactions.clear();
}//void ReferenceLineInfo::reset()	  


bool ReferenceLineInfo::operator==(const ReferenceLineInfo &rhs) const
{
  return ((m_validity == rhs.m_validity)
          && (m_has_coincidences == rhs.m_has_coincidences)
          && (m_input == rhs.m_input)
          && (m_source_type == rhs.m_source_type)
          && (m_nuclide == rhs.m_nuclide)
          && (m_element == rhs.m_element)
          && (m_reactions == rhs.m_reactions));
}//bool operator==(const ReferenceLineInfo &rhs) const

void ReferenceLineInfo::toJson( string &json ) const
{
  // For reference, for Th232, the JSON returned is about 32 kb, and U238 is 70.5 kb (just gamma and xray).
  //  TODO: The "decay" for each line could be specified in a separate map, so like 'Ba133 to Cs133 via Electron Capture' isnt included in the JSON a bunch of times.  could also change "particle" to "p", "decay" to "d", and particle values from "gamma", "xray", etc, to "g", "x", etc.
  
  std::stringstream jsons;
  
  // We will put individual lines descriptions in a separate array, and give the index into
  //  the array that line should go; this is because decays like U238 have many gammas for
  //  many of the transitions, so we can save a decent about of space.
  //  For {U238,Th232,Ba133}, the JSON size goes from {71, 37, 1} kb, to {52,27,1} kb
  vector<string> decay_strs;
  map<string,size_t> decay_str_indexs;
  
  jsons << "{\"color\":\"" << (m_input.m_color.isDefault() ? "#0000FF" : m_input.m_color.cssText(false)) << "\","
  << "\"parent\":\"" << m_input.m_input_txt << "\",";
  
  if( m_input.m_promptLinesOnly )
    jsons << "\"prompt\":true,";
  if( m_source_type == ReferenceLineInfo::SourceType::Background )
    jsons << "\"age\":\"Primordial\",";
  else if( !m_input.m_age.empty() )
    jsons << "\"age\":" << jsQuote( m_input.m_age ) << ",";
  
  if( !m_input.m_detector_name.empty() )
    jsons << "\"detector\":" << jsQuote( m_input.m_detector_name ) << ",";
  
  if( !m_input.m_shielding_name.empty() && !m_input.m_shielding_thickness.empty() )
  {
    assert( m_input.m_shielding_an.empty() );
    assert( m_input.m_shielding_ad.empty() );
#ifndef NDEBUG
    PhysicalUnits::stringToDistance(m_input.m_shielding_thickness);
#endif

    jsons << "\"shielding\":" << jsQuote( m_input.m_shielding_name ) << ",";
    jsons << "\"shieldingThickness\":" << jsQuote( m_input.m_shielding_thickness ) << ",";
  }else if( !m_input.m_shielding_an.empty() && !m_input.m_shielding_ad.empty() )
  {
    assert( m_input.m_shielding_name.empty() );
    assert( m_input.m_shielding_thickness.empty() );
#ifndef NDEBUG
    double dummy;
    assert( (stringstream(m_input.m_shielding_an) >> dummy) );
    assert( (stringstream(m_input.m_shielding_ad) >> dummy) );
#endif
    
    const string shield = "AN=" + m_input.m_shielding_an + ", AD=" + m_input.m_shielding_ad + " g/cm2";
    jsons << "\"shielding\":" << jsQuote( shield ) << ",";
  }

  
  jsons << "\"lines\":[";
  
  bool printed = false;
  char intensity_buffer[32] = { '\0' };
  
  // Round to the nearest 10 eV; probably the extent to which any data useful, or even good to
  const auto round_energy = []( const double e ) -> double { return std::round(100.0*e)/100.0; };
  
  for( size_t index = 0; index < m_ref_lines.size(); ++index )
  {
    const RefLine &line = m_ref_lines[index];
    if( line.m_normalized_intensity <= std::numeric_limits<float>::min() )  //numeric_limits<float>::min()==1.17549e-38
      continue;
    
    const double energy = round_energy( line.m_energy );
    
    // There are situations where two lines have either the exact same energies, or super-close
    //  energies, and the spectrum chart is not particularly smart about this, so we effectively
    //  lose some amplitude, so we'll combine them here.
    // However, this additional munging has a overhead for a rare edge-case, so we'll split
    //  the code-paths, even though this adds code...
    // For U238, there are 39 pairs of lines that get combined, and 1 triplet of lines combined.
    
    // We will assume entries are sorted by energy, which is only guaranteed when
    //  ReferencePhotopeakDisplay::updateDisplayChange() sets the data - but also, I think this
    //  is the only place that sets the data!
    
    // TODO: This current way of doing things will sum a Cascade sum with a gamma (e.g., the 387.8
    //       keV of U235), which is probably not the right thing to do because it then shows up as
    //       giant gamma on the chart, which is deceptive
    const bool next_gamma_close = (((index+1) < m_ref_lines.size())
                                   && (round_energy(m_ref_lines[index+1].m_energy) == energy)
                                   //&& (m_ref_lines[index+1].m_source_type == line.m_source_type)
                                   //&& (m_ref_lines[index+1].m_particle_type == line.m_particle_type)
                                   //&& (m_ref_lines[index+1].m_particle_type == RefLine::Particle::Gamma)
                                   );
    
    if( next_gamma_close )
    {
      double intensity = 0.0;
      size_t num_combined = 0;
      map<const SandiaDecay::Nuclide *,double> nuc_to_frac;
      // TODO: be a little more efficient than allocating strings in these sets...
      set<string> particles, decays;
      vector<pair<string,double>> decay_fractions;
      
      for( size_t inner_index = index; inner_index < m_ref_lines.size(); ++inner_index )
      {
        const RefLine &inner_line = m_ref_lines[inner_index];
        
        const double this_energy = round_energy(inner_line.m_energy);
        if( this_energy != energy )
          break;
        
        if( inner_line.m_normalized_intensity <= 0.0 )
        {
          num_combined += 1;
          continue;
        }
        
        intensity += inner_line.m_normalized_intensity;
        particles.insert( inner_line.particlestr() );
        if( !inner_line.m_decaystr.empty() )
        {
          decays.insert( inner_line.m_decaystr );
          
          auto pos = std::lower_bound( begin(decay_fractions), end(decay_fractions),
                                      inner_line.m_decaystr,
                                      []( const pair<string,double> &el, const string &val) -> bool {
            return el.first < val;
          });
          
          if( (pos != end(decay_fractions)) && (inner_line.m_decaystr == pos->first) )
            pos->second += inner_line.m_normalized_intensity;
          else
            decay_fractions.insert( pos, {inner_line.m_decaystr, inner_line.m_normalized_intensity} );
          
          assert( decays.size() == decay_fractions.size() );
        }
        
        auto pos = nuc_to_frac.find(inner_line.m_parent_nuclide);
        if( pos == end(nuc_to_frac) )
          pos = nuc_to_frac.insert( {inner_line.m_parent_nuclide, 0.0} ).first;
        assert( pos != end(nuc_to_frac) );
        pos->second += inner_line.m_normalized_intensity;
        
        num_combined += 1;
      }//for( loop over inner_index to find all energies to cluster together )
      
      assert( num_combined != 0 ); //could tighten this up to (num_combined > 1)
      
      if( IsNan(intensity) || IsInf(intensity) )
        snprintf( intensity_buffer, sizeof(intensity_buffer), "0" );
      else
        snprintf( intensity_buffer, sizeof(intensity_buffer), "%.3g", intensity );
      
      auto combine_particle_strs = []( const set<string> &strs ) -> string {
        string answer;
        for( const auto &s : strs )
          answer += (answer.empty() ? "" : ", ") + s;
        return answer;
      };//combine_particle_strs lambda
      
      jsons << (printed ? "," : "") << "{\"e\":" << energy << ",\"h\":" << intensity_buffer;
      if( !particles.empty() )
        jsons << ",\"particle\":" << jsQuote(combine_particle_strs(particles));
      
      assert( decays.size() == decay_fractions.size() );
      
      if( !decay_fractions.empty() )
      {
        auto combine_decays = []( const vector<pair<string,double>> &fracs ) -> string {
          if( fracs.size() == 1 )
            return fracs[0].first;
          
          double sum = 0.0;
          for( const auto &i : fracs )
            sum += i.second;
          
          string answer;
          for( const auto &i : fracs )
          {
            answer += (answer.empty() ? "" : ", ") + i.first
                      + " (" + SpecUtils::printCompact( 100.0*i.second/sum, 4 ) + "%)";
          }
          return answer;
        };//combine_particle_strs lambda
        
        const string decay_str = combine_decays( decay_fractions );
        auto decay_str_iter = decay_str_indexs.find(decay_str);
        if( decay_str_iter == end(decay_str_indexs) )
        {
          decay_str_iter = decay_str_indexs.insert( {decay_str, decay_strs.size()} ).first;
          decay_strs.push_back(decay_str);
        }
        
        jsons << ",\"desc_ind\":" << decay_str_iter->second;
      }
      
      // For Nuclide Mixtures, need to add a special label to lines to indicate ultimate parent,
      //  to show when you mouse-over a line.
      if( (m_source_type == ReferenceLineInfo::SourceType::NuclideMixture) && line.m_parent_nuclide )
      {
        //Now need to make 'src_label'
        vector<pair<const SandiaDecay::Nuclide *,double>> srcs( begin(nuc_to_frac), end(nuc_to_frac) );
        std::sort( begin(srcs), end(srcs),
                  []( const pair<const SandiaDecay::Nuclide *,double> &lhs,
                      const pair<const SandiaDecay::Nuclide *,double> &rhs) -> bool
                  { return lhs.second > rhs.second; }
        );
        
        string src_label;
        if( srcs.size() <= 1 )
        {
          src_label = line.m_parent_nuclide->symbol;
        }else
        {
          // List the leading 3 sources, at most, putting their relative fraction in parenthesis
          int nsrcs = 0;
          for( const auto &src : srcs )
          {
            if( nsrcs++ >= 3 )
            {
              src_label += ",...";
              break;
            }
            
            src_label += (src_label.empty() ? "" : ", ")
                      + (src.first ? src.first->symbol : string("null"))
                      + " (" + SpecUtils::printCompact(src.second/intensity, 3) + ")";
          }//for( const auto &src : srcs )
        }//if( srcs.size() <= 1 ) / else
        
        jsons << ",\"src_label\":\"" << src_label << "\"";
      }//if( ReferenceLineInfo::SourceType::NuclideMixture )
      
      // Now increment 'i' so we'll skip over these lines we've already covered.
      index += (num_combined >= 1) ? (num_combined - 1) : size_t(0);
    }else  //if( next_gamma_close )
    {
      if( IsNan(line.m_normalized_intensity) || IsInf(line.m_normalized_intensity) )
        snprintf( intensity_buffer, sizeof(intensity_buffer), "0" );
      else
        snprintf( intensity_buffer, sizeof(intensity_buffer), "%.3g", line.m_normalized_intensity );
      
      jsons << (printed ? "," : "") << "{\"e\":" << energy << ",\"h\":" << intensity_buffer;
      jsons << ",\"particle\":" << jsQuote( line.particlestr() );
      
      if( !line.m_decaystr.empty() )
      {
        auto decay_str_iter = decay_str_indexs.find(line.m_decaystr);
        if( decay_str_iter == end(decay_str_indexs) )
        {
          decay_str_iter = decay_str_indexs.insert( {line.m_decaystr, decay_strs.size()} ).first;
          decay_strs.push_back(line.m_decaystr);
        }
        
        jsons << ",\"desc_ind\":" << decay_str_iter->second;
      }//if( !line.m_decaystr.empty() )
      
      // For Nuclide Mixtures, need to add a special label to lines to indicate ultimate parent,
      //  to show when you mouse-over a line.
      if( (m_source_type == ReferenceLineInfo::SourceType::NuclideMixture) && line.m_parent_nuclide )
        jsons << ",\"src_label\":\"" << line.m_parent_nuclide->symbol << "\"";
    }//if( next gamma line is close ) / else
    
    if( !line.m_color.isDefault() )
      jsons << ",\"color\":\"" << line.m_color.cssText(false) << "\"";
    
    jsons << "}";
    
    printed = true;
  }//for( size_t index = 0; index < m_ref_lines.size(); ++index )
  jsons <<"]";

  
  jsons <<", \"desc_strs\":[";
  for( size_t i = 0; i < decay_strs.size(); ++i )
    jsons << (i ? "," : "") << jsQuote( decay_strs[i] );
  jsons <<"]}";
  
  json += jsons.str();
}//std::string toJson( const ReferenceLineInfo &displnuc )


void ReferenceLineInfo::sortByEnergy()
{
  std::sort( begin( m_ref_lines ), end( m_ref_lines ), 
    []( const RefLine &lhs, const RefLine &rhs ) -> bool {
    if( lhs.m_energy == rhs.m_energy )
      return lhs.m_normalized_intensity < rhs.m_normalized_intensity;
    return lhs.m_energy < rhs.m_energy;
    } );
}//void sortByEnergy()


void RefLineInput::deSerialize( const rapidxml::xml_node<char> *base_node )
{
  // The schema used was for #ReferenceLineInfo, before 20221220, so its
  //  not exactly what you would write from scratch, but good-enough
  //  considering it gives us backward compatibility, mostly.
  // <DisplayedSource> is for pre-20221220, and <RefLineInput> is
  // post-20221220 (so that old-code wont even try to deserialize the stuff).
  
  if( !base_node
     || ((base_node->name() != string("DisplayedSource"))
         && (base_node->name() != string("RefLineInput"))) )
    throw runtime_error( "Invalid base node for DisplayedSource" );
  
  static_assert( sm_xmlSerializationVersion == 2, "You need to update xml serialization" );
  
  int version;
  const rapidxml::xml_attribute<char> * const version_attr = base_node->first_attribute( "version", 7 );
  if( !version_attr || !version_attr->value()
     || !(stringstream(version_attr->value()) >> version)
     || (version != 0 && version != 1 && version != 2) )
    throw runtime_error( "Missing or invalid DisplayedSource version" );
 
  reset();
  RefLineInput &input = *this;
  
  const rapidxml::xml_node<char> *node = base_node->first_node( "Nuclide", 7 );
  if( node && node->value_size() )
    input.m_input_txt = node->value();
  
  node = base_node->first_node( "Element", 7 );
  if( node && node->value_size() && input.m_input_txt.empty() )
    input.m_input_txt = node->value();
  
  node = base_node->first_node( "LineColor", 9 );
  if( node && node->value_size() )
    try{ input.m_color = Wt::WColor( node->value() ); }catch(...){ }
  
  //node = base_node->first_node( "DisplayLines", 12 );
  //if( node && node->value_size() )
  //  displayLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowGammas", 10 );
  if( node && node->value_size() )
    input.m_showGammas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowXrays", 9 );
  if( node && node->value_size() )
    input.m_showXrays = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowAlphas", 10 );
  if( node && node->value_size() )
    input.m_showAlphas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowBetas", 9 );
  if( node && node->value_size() )
    input.m_showBetas = (node->value()[0] == '1');
  
  node = base_node->first_node("ShowCascades", 12);
  if( node && node->value_size() )
    input.m_showCascades = (node->value()[0] == '1');
  
  node = base_node->first_node("ShowEscapes", 11);
  input.m_showEscapes = node && node->value_size() && (node->value()[0] == '1');
  
  //node = base_node->first_node( "ShowLines", 9 );
  //if( node && node->value_size() )
  //  showLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "PromptLinesOnly", 15 );
  if( node && node->value_size() )
    input.m_promptLinesOnly = (node->value()[0] == '1');
  
  //node = base_node->first_node( "IsBackground", 12 );
  //if( node && node->value_size() )
  //  isOtherRef = (node->value()[0] == '1');
  
  //node = base_node->first_node( "IsReaction", 10 );
  //if( node && node->value_size() )
  //  isReaction = (node->value()[0] == '1');
  
  node = base_node->first_node( "Age", 3 );
  if( node && node->value_size() )
    input.m_age = node->value();
  
  node = base_node->first_node( "LowestBranchRatio", 17 );
  if( node && node->value_size() && !(stringstream(node->value()) >> input.m_lower_br_cutt_off) )
    throw runtime_error( "Invalid LowerBrCuttoff parameter" );
  
  node = base_node->first_node( "ShieldingName", 13 );
  if( node && node->value_size() )
  {
    input.m_shielding_name = node->value();
    
    // Pre 20221220, generic shielding would look like "26, 12.1 g/cm2", so we'll convert
    //  this to AN and AD.
    if( (version < 2) && SpecUtils::icontains(input.m_shielding_name, "g/cm2") )
    {
      vector<string> parts;
      SpecUtils::split( parts, input.m_shielding_name, "," );
      if( parts.size() == 2 )
      {
        input.m_shielding_an = SpecUtils::trim_copy(parts[0]);
        input.m_shielding_ad = SpecUtils::trim_copy(parts[1]);
        double dummy;
        if( !(stringstream(input.m_shielding_an) >> dummy)
           || !(stringstream(input.m_shielding_ad) >> dummy) )
        {
          input.m_shielding_an = input.m_shielding_ad = "";
          cerr << "Old RefLineInput Shielding '" << input.m_shielding_name << "' was not valid AD, AN." << endl;
        }else
        {
          input.m_shielding_name = "";
        }
      }else
      {
        cerr << "\n\ninput.m_shielding_name_='" << input.m_shielding_name << "', isnt formed as expected" << endl;
      }
    }
    
    node = base_node->first_node( "ShieldingThickness", 18 );
    
    if( node && node->value_size() )
    {
      input.m_shielding_thickness = node->value();
      try
      {
        // Make sure a proper distance
        PhysicalUnits::stringToDistance(input.m_shielding_thickness);
      }catch( std::exception &e )
      {
        input.m_shielding_thickness = "";
        cerr << "Unexpected invalid shielding thickness when deserializing: " << e.what() << endl;
      }
    }//if( node && node->value_size() )
  }//if( node && node->value_size() && node->value_size() )
  
  node = base_node->first_node( "ShieldingAN", 11 );
  if( node && node->value_size() )
  {
    input.m_shielding_an = node->value();
  
    double dummy;
    if( !(stringstream(input.m_shielding_an) >> dummy) )
    {
      cerr << "Unexpected invalid AN (" << input.m_shielding_an << ")" << endl;
      input.m_shielding_an = input.m_shielding_ad = "";
    }else
    {
      assert( input.m_shielding_name.empty() );
      assert( input.m_shielding_thickness.empty() );
      
      input.m_shielding_name = input.m_shielding_thickness = "";
    }//
  }//if( node && node->value_size() )
  
  
  node = base_node->first_node( "ShieldingAD", 11 );
  if( node && node->value_size() )
  {
    input.m_shielding_ad = node->value();
    double dummy;
    if( !(stringstream(input.m_shielding_ad) >> dummy) )
    {
        cerr << "Unexpected invalid AD (" << input.m_shielding_ad << ")" << endl;
        input.m_shielding_an = input.m_shielding_ad = "";
    }else
    {
      assert( input.m_shielding_name.empty() );
      assert( input.m_shielding_thickness.empty() );
      input.m_shielding_name = input.m_shielding_thickness = "";
    }//
  }//if( node && node->value_size() )
  
  
  if( !input.m_shielding_an.empty() || !input.m_shielding_ad.empty() )
  {
    assert( input.m_shielding_name.empty() );
    assert( input.m_shielding_thickness.empty() );
    input.m_shielding_name.clear();
    input.m_shielding_thickness.clear();
  }
  
  if( !input.m_shielding_name.empty() || !input.m_shielding_thickness.empty() )
  {
    assert( input.m_shielding_an.empty() );
    assert( input.m_shielding_ad.empty() );
    input.m_shielding_an.clear();
    input.m_shielding_ad.clear();
  }
  
  node = base_node->first_node( "DetectorName", 12 );
  if( node && node->value_size() )
    input.m_detector_name = node->value();
}//void deSerialize( const rapidxml::xml_node<char> *node );


void RefLineInput::serialize( rapidxml::xml_node<char> *parent_node ) const
{
  // See notes at top of #RefLineInput::deSerialize for history of the function.
  
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  char buffer[64];
  
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node, *lines_node, *sf_node,
  *part_sf_node, *energy_node, *br_node, *str_node;
  rapidxml::xml_attribute<char> *attr;
  
  name = "RefLineInput";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  static_assert( sm_xmlSerializationVersion == 2, "You need to update xml serialization" );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  snprintf( buffer, sizeof(buffer), "%i", sm_xmlSerializationVersion );
  value = doc->allocate_string( buffer );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  name = "Nuclide";
  value = doc->allocate_string( m_input_txt.c_str() );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  if( !m_color.isDefault() )
  {
    name = "LineColor";
    value = doc->allocate_string( m_color.cssText(false).c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }
  
  if( !m_age.empty() )
  {
    name = "Age";

    string age_str = m_age;
    
    {// begin scope to convert age into English
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      const SandiaDecay::Nuclide *nuc = db->nuclide( m_input_txt );
      const double hl = nuc ? nuc->halfLife : -1.0;
      
      try
      {
        // If we can convert from text to duration, using non-localized function, we probably
        //  dont need to convert things
        PhysicalUnits::stringToTimeDurationPossibleHalfLife( age_str, hl );
      }catch( std::exception &e )
      {
        const Wt::WLocale &locale = Wt::WLocale::currentLocale();
        if( !locale.name().empty() && !SpecUtils::istarts_with(locale.name(), "en" ) )
        {
          try
          {
            const double duration = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( age_str, hl );
            age_str = PhysicalUnits::printToBestTimeUnits( duration );
          }catch( std::exception &e )
          {
            cerr << "Error converting age ('" << age_str << "') to English time-span: " << e.what() << endl;
          }
        }//if( we are in a non-English locale )
      }//try / catch - see if we can read it in
    }// End scope to convert age into English
    
    value = doc->allocate_string( age_str.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( age >= 0.0 )
  
  name = "LowestBranchRatio";
  snprintf( buffer, sizeof(buffer), "%g", m_lower_br_cutt_off );
  value = doc->allocate_string( buffer );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  //name = "DisplayLines";
  //value = (displayLines ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  name = "ShowGammas";
  value = (m_showGammas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowXrays";
  value = (m_showXrays ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowAlphas";
  value = (m_showAlphas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowBetas";
  value = (m_showBetas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowCascades";
  value = (m_showCascades ? "1" : "0");
  node = doc->allocate_node(rapidxml::node_element, name, value);
  base_node->append_node(node);

  name = "ShowEscapes";
  value = (m_showEscapes ? "1" : "0");
  node = doc->allocate_node(rapidxml::node_element, name, value);
  base_node->append_node(node);
  
  name = "PromptLinesOnly";
  value = (m_promptLinesOnly ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  //name = "ShowLines";
  //value = (showLines ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  //name = "IsBackground";
  //value = (isOtherRef ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  //name = "IsReaction";
  //value = (isReaction ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  if( !m_shielding_name.empty() )
  {
    assert( m_shielding_an.empty() );
    assert( m_shielding_ad.empty() );
    
    name = "ShieldingName";
    value = doc->allocate_string( m_shielding_name.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( shieldingName.size() )
  
  if( !m_shielding_thickness.empty() )
  {
    assert( m_shielding_an.empty() );
    assert( m_shielding_ad.empty() );
    
    name = "ShieldingThickness";
    value = doc->allocate_string( m_shielding_thickness.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }
  
  
  if( !m_shielding_an.empty() || !m_shielding_ad.empty() )
  {
    assert( m_shielding_name.empty() );
    assert( m_shielding_thickness.empty() );
    
#ifndef NDEBUG
    double dummy;
    assert( m_shielding_an.empty() || (stringstream(m_shielding_an) >> dummy) );
    assert( m_shielding_ad.empty() || (stringstream(m_shielding_ad) >> dummy) );
#endif
    
    if( !m_shielding_an.empty() )
    {
      name = "ShieldingAN";
      value = doc->allocate_string( m_shielding_an.c_str() );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      base_node->append_node( node );
    }
    
    if( !m_shielding_ad.empty() )
    {
      name = "ShieldingAD";
      value = doc->allocate_string( m_shielding_ad.c_str() );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      base_node->append_node( node );
    }
  }//if( !m_shielding_an.empty() && !m_shielding_ad.empty() )
  
  
  if( !m_detector_name.empty() )
  {
    name = "DetectorName";
    value = doc->allocate_string( m_detector_name.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( detectorName.size() )
}//void ReferenceLineInfo::serialize( rapidxml::xml_node<char> *parent_node )


void RefLineInput::setShieldingAttFcn( const MaterialDB *db )
{
  // Check for generic shielding
  if( !m_shielding_an.empty() && !m_shielding_ad.empty() )
  {
    assert( m_shielding_name.empty() );
    assert( m_shielding_thickness.empty() );
    
    double atomic_number, areal_density;
    if( !(stringstream(m_shielding_an) >> atomic_number) )
      throw runtime_error( "Couldnt convert atomic number string '"
                           + m_shielding_an + "' to a floating point." );
    
    if( !(stringstream(m_shielding_ad) >> areal_density) )
      throw runtime_error( "Couldnt convert areal density string '"
                          + m_shielding_ad + "' to a floating point." );
    
    if( (atomic_number < MassAttenuation::sm_min_xs_atomic_number)
       || (atomic_number > MassAttenuation::sm_max_xs_atomic_number) )
      throw runtime_error( "Invalid atomic value '" + m_shielding_an + "'." );
    
    if( areal_density < 0.0 )
      throw runtime_error( "Invalid areal density value '" + m_shielding_ad + "'." );
    
    areal_density *= (PhysicalUnits::gram / PhysicalUnits::cm2);
    
    m_shielding_att = [atomic_number, areal_density]( float energy ) -> double {
      const double att_coef = GammaInteractionCalc::transmition_coefficient_generic( atomic_number, areal_density, energy );
      return exp( -1.0 * att_coef );
    };
    
    return;
  }//if( !m_shielding_an.empty() && !m_shielding_ad.empty() )
  
  assert( m_shielding_an.empty() );
  assert( m_shielding_ad.empty() );
  
  if( m_shielding_name.empty() || m_shielding_thickness.empty() )
  {
    m_shielding_att = nullptr;
    return;
  }
  
  if( m_shielding_name.empty() )
  {
    assert( m_shielding_thickness.empty() );
    m_shielding_att = nullptr;
    return;
  }//if( m_shielding_name.empty() )
  
  if( !db )
    throw runtime_error( "No MaterialDB passed in." );
  
  const Material * const material_raw = db->material(m_shielding_name);
  if( !material_raw )
    throw runtime_error( "No material named '" + m_shielding_name + "' available." );
  
  const auto material = make_shared<Material>(*material_raw);
  const double thickness = PhysicalUnits::stringToDistance( m_shielding_thickness );
  if( thickness < 0.0 )
    throw runtime_error( "Distance '" + m_shielding_thickness + "' is negative." );
  
  m_shielding_att = [material, thickness]( float energy ) -> double {
    const double att_coef = GammaInteractionCalc::transmition_coefficient_material( material.get(), energy, thickness );
    return exp( -1.0 * att_coef );
  };
}//void ReferenceLineInfo::setShieldingAttFcn( const MaterialDB *db )


std::shared_ptr<ReferenceLineInfo> ReferenceLineInfo::generateRefLineInfo( RefLineInput input )
{
  // The gamma or xray energy below which we wont show lines for.
  //  x-rays for nuclides were limited at above 10 keV, so we'll just impose this
  //  as a lower limit to show to be consistent.
  const float lower_photon_energy = 10.0f;
  
  auto answer_ptr = make_shared<ReferenceLineInfo>();
  ReferenceLineInfo &answer = *answer_ptr;
  
  
  // We want to set the final _modified_ version of input to the answer, before
  // returning, so we'll just use a helper for this.
  //on_scope_exit on_exit( [&answer_ptr, &input](){
  //  answer_ptr->m_input = input;
  //  answer_ptr->lineColor = input.m_color;
  //  } );
  answer_ptr->m_input = input;
  
  if( input.m_input_txt.empty() )
  {
    answer.m_validity = ReferenceLineInfo::InputValidity::Blank; //Should already be this value, but being explicit
    
    return answer_ptr;
  }
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  //Inputs like "Fission U233 Fast, 1d buildup" will result in a valid nuclide of U233,
  //  so check for fission or reaction, and if so null-out nuclide
  const bool maybe_fission_rct = (SpecUtils::icontains(input.m_input_txt, "Fission")
                                  || SpecUtils::icontains(input.m_input_txt, "("));
  
  
  double age = 0.0;
  const SandiaDecay::Nuclide * const nuc = maybe_fission_rct ? nullptr
                                                             : db->nuclide( input.m_input_txt );
  
  if( nuc )
  {
    input.m_input_txt = nuc->symbol;
    input.m_promptLinesOnly = (input.m_promptLinesOnly && nuc->canObtainPromptEquilibrium());
    
    answer.m_nuclide = nuc;
    
    if( input.m_promptLinesOnly )
    {
      age = 5.0*nuc->promptEquilibriumHalfLife();
      input.m_age = PhysicalUnitsLocalized::printToBestTimeUnits(age, 2);
    }else if( input.m_age == "" )
    {
      age = PeakDef::defaultDecayTime( nuc, &input.m_age );
    }else if( nuc->isStable() )
    {
      // age = 0;
      input.m_age = "";
      answer.m_input_warnings.push_back( nuc->symbol + " is a stable isotope." );
      answer.m_validity = ReferenceLineInfo::InputValidity::InvalidSource;
      answer_ptr->m_input = input;
      
      return answer_ptr;
    }else
    {
      try
      {
        age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( input.m_age, nuc->halfLife );
      }catch( std::exception & )
      {
        answer.m_input_warnings.push_back( "Invalid nuclide age input." );
        answer.m_validity = ReferenceLineInfo::InputValidity::InvalidAge;
        
        answer_ptr->m_input = input;
        
        return answer_ptr;
      }//try /catch to get the age
      
      if( age > 100.0 * nuc->halfLife || age < 0.0 )
      {
        const string old_age_str = input.m_age;
        age = PeakDef::defaultDecayTime( nuc, &input.m_age );
        answer.m_input_warnings.push_back( "Changed age to a more reasonable value for "
                                          + nuc->symbol + " from '" + old_age_str + "' to " + input.m_age );
      }
    }//if( prompt only ) / else
  }//if( nuc )
  
  
  const bool check_element = (!nuc && (input.m_input_txt.find_first_of( "0123456789" ) == string::npos));
  const SandiaDecay::Element * el = check_element ? db->element( input.m_input_txt ) : nullptr;
  if( !el && check_element
     && (SpecUtils::icontains(input.m_input_txt, "xray") || SpecUtils::icontains(input.m_input_txt, "x-ray")) )
  {
    string elstr = input.m_input_txt;
    size_t pos = SpecUtils::ifind_substr_ascii( elstr, "xray");
    if( pos == string::npos )
      pos = SpecUtils::ifind_substr_ascii( elstr, "x-ray");
    if( pos != string::npos )
      elstr = elstr.substr(0,pos);
    SpecUtils::ireplace_all( elstr, "-_\t ,", "" );
    el = db->element( elstr );
  }//
  
  if( el )
  {
    input.m_input_txt = el->symbol;
    input.m_age = ""; //JIC
    input.m_showXrays = true;
    
    answer.m_element = el;
  }//if( !nuc )
  
  
  string reaction_txt; //CSV list of reactions - I think for our context, only ever a single reaction
  vector<ReactionGamma::ReactionPhotopeak> rctn_gammas;
  if( !nuc && !el )
  {
    const size_t open_paren = input.m_input_txt.find( "(" );
    const size_t close_paren = (open_paren == string::npos) ? string::npos
    : input.m_input_txt.find( ")", open_paren );
    
    if( close_paren != string::npos )
    {
      try
      {
        const ReactionGamma *rctnDb = ReactionGammaServer::database();
        if( rctnDb )
        {
          reaction_txt = rctnDb->gammas( input.m_input_txt, rctn_gammas );
          SpecUtils::ireplace_all( reaction_txt, "'", "" );
          
          // We will fill in answer.m_reactions later on.
          
          // Note: we are-not setting the input to the reaction_txt, as we dont need to bother
          //       the user that the underlying data is for the isotopics of the element, and also
          //       we dont want to change the user input, because we want to keep contributions
          //       normalized to natural abundance.
          //input.m_input_txt = reaction_txt;
          
          input.m_age = "";
          input.m_showGammas = true;
        }
      }catch( std::exception &e )
      {
        // Not a reaction
      }//try / catch
    }//if( (open_paren != string::npos) && (close_paren != string::npos) )
  }//if( !nuc && !el )
  
  
  vector<OtherRefLine> otherRefLinesToShow;
  const bool is_background = (nuc || el || !rctn_gammas.empty())
                               ? false
                               : SpecUtils::icontains( input.m_input_txt, "background" );
  if( is_background )
  {
    input.m_input_txt = "background";
    
    for( const OtherRefLine &bl : getBackgroundRefLines() )
    {
      const bool isXray = (std::get<3>(bl) == OtherRefLineType::BackgroundXRay);
      if( (isXray && input.m_showXrays) || (!isXray && input.m_showGammas) )
        otherRefLinesToShow.push_back( bl );
    }//for( const BackgroundLine &bl : getBackgroundRefLines() )
    
    input.m_age = "";
  }//if( is_background )
  
  bool is_custom_energy = false;
  if( !nuc && !el && rctn_gammas.empty() && !is_background )
  {
    try
    {
      const auto lineType = OtherRefLineType::OtherBackground;
      const float energy = static_cast<float>(PhysicalUnits::stringToEnergy( input.m_input_txt ));
      
      //BackgroundLine<Energy, RelBranchRatio, "Symbol", OtherRefLineType, "Description">
      OtherRefLine line{energy, 1.0f, "", lineType, input.m_input_txt};
      otherRefLinesToShow.push_back( line );
      
      if( input.m_showEscapes && (energy > 1022.0) )
      {
        const float se_energy = static_cast<float>( energy - 510.9989 );
        const float de_energy = static_cast<float>( energy - 2.0*510.9989 );
        const float se_rel_amp = static_cast<float>( ns_single_escape_sf(energy) );
        const float de_rel_amp = static_cast<float>( ns_double_escape_sf(energy) );
        
        OtherRefLine se_line{ se_energy, se_rel_amp, "", lineType, "Generic HPGe S.E. fraction" };
        OtherRefLine de_line{ de_energy, de_rel_amp, "", lineType, "Generic HPGe D.E. fraction" };
        
        otherRefLinesToShow.push_back( se_line );
        otherRefLinesToShow.push_back( de_line );
      }//if( input.m_showEscapes && (energy > 1022.0) )
      
      is_custom_energy = true;
      input.m_age = "";
    }catch( std::exception & )
    {
    }
  }//if( !nuc && !el && rctnGammas.empty() )
  
  const ReferenceLinePredef::NucMix *nuc_mix = nullptr;
  const ReferenceLinePredef::CustomSrcLines *src_lines = nullptr;
  if( !nuc && !el && rctn_gammas.empty() && !is_background && !is_custom_energy )
  {
    nuc_mix = get_custom_nuc_mix( input.m_input_txt );
    if( nuc_mix )
    {
      if( input.m_age == "" )
      {
        input.m_age = nuc_mix->m_default_age_str;
        age = nuc_mix->m_default_age;
      }else
      {
        try
        {
          age = PhysicalUnitsLocalized::stringToTimeDuration( input.m_age );
        }catch( std::exception & )
        {
          answer.m_input_warnings.push_back( "Invalid age input." );
          answer.m_validity = ReferenceLineInfo::InputValidity::InvalidAge;
          answer_ptr->m_input = input;
          return answer_ptr;
        }//try /catch to get the age
      }//if( input.m_age == "" ) / else
    }//if( nuc_mix )
    
    if( !nuc_mix )
    {
      src_lines = get_custom_src_lines( input.m_input_txt );
      if( src_lines )
      {
        age = 0.0;
        input.m_age = "";
      }
    }//if( !nuc_mix )
  }//if( not anything else so far )
  
  
  // Check if Fission product spectrum, and if so, parse out all the information.
  //  Should probably move all this information parsing to a seperate function (its
  //  length has grown quite a bit over initial imp).
  bool fission_src = false;
  FissionType fission_type = FissionType::Thermal;
  const SandiaDecay::Nuclide *fission_nuclide = nullptr;
  double fission_buildup_time = 0.0;
  if( !nuc && !el && rctn_gammas.empty() && !is_background && !is_custom_energy
     && SpecUtils::icontains( input.m_input_txt, "Fission" ) )
  {
    //Input is of the form "Fission U238 Thermal, 1d buildup"
    //  But we should make it a little more robust to formatting,
    //  so we'll kinda go through it searching for things
    string input_txt = input.m_input_txt;
    SpecUtils::ireplace_all( input_txt, "Fission", "" );
    SpecUtils::ireplace_all( input_txt, "  ", " " );
    SpecUtils::trim( input_txt );
    
    bool found_neut_type = true;
    if( SpecUtils::icontains(input_txt, "Thermal") )
    {
      fission_type = FissionType::Thermal;
      SpecUtils::ireplace_all( input_txt, "Thermal", "" );
    }else if( SpecUtils::icontains(input_txt, "Fast") )
    {
      fission_type = FissionType::Fast;
      SpecUtils::ireplace_all( input_txt, "Fast", "" );
    }else if( SpecUtils::icontains(input_txt, "14 MeV")
             || SpecUtils::icontains(input_txt, "14MeV")
             || SpecUtils::icontains(input_txt, "14-MeV"))
    {
      fission_type = FissionType::FourteenMeV;
      SpecUtils::ireplace_all( input_txt, "14MeV", "" );
      SpecUtils::ireplace_all( input_txt, "14 MeV", "" );
      SpecUtils::ireplace_all( input_txt, "14-MeV", "" );
    }else
    {
      found_neut_type = false;
      
      // Maybe we should just assume thermal? E.g.:
      //fission_type = FissionType::Thermal;
      //answer.m_input_warnings.push_back( "Count not interpret a neutron energy type"
      //                               " - assuming thermal (fast and 14MeV are other options)." );
    }
    
    if( found_neut_type )
    {
      fission_nuclide = db->nuclide( input_txt );
      if( !fission_nuclide )
      {
        vector<string> remaining_fields;
        SpecUtils::split( remaining_fields, input_txt, ", " );
        for( size_t i = 0; !fission_nuclide && (i < remaining_fields.size()); ++i )
        {
          fission_nuclide = db->nuclide( remaining_fields[i] );
          if( !fission_nuclide && ((i + 1) < remaining_fields.size()) )
            fission_nuclide = db->nuclide( remaining_fields[i] + remaining_fields[i+1] );
        }
      }//if( !fission_nuclide )
      
      //We want to remove the text that led to this nuclide from the string, but we dont know it
      //  here, so we'll do a bit of a work around, and delete the numbers from the nuclide
      //  and just assume a meta-stable state should be here
      if( fission_nuclide && (fission_nuclide->isomerNumber > 0) )
        fission_nuclide = nullptr;
      
      if( fission_nuclide )
      {
        const string numstr = std::to_string(static_cast<int>(fission_nuclide->massNumber));
        auto pos = input_txt.find( numstr );
        assert( pos != string::npos );
        if( pos != string::npos )
          input_txt.erase(pos, pos + numstr.size());
      }//if( fission_nuclide )
    }//if( found_neut_type )
    
    if( fission_nuclide )
    {
      bool found_buildup = false;
      vector<string> remaining_fields;
      SpecUtils::split( remaining_fields, input_txt, ", " );
      for( size_t i = 0; !found_buildup && (i < remaining_fields.size()); ++i )
      {
        const string &field = remaining_fields[i];
        const char c = field.empty() ? ' ' : field[0];
        if( field.empty() || (((c < '0') || (c > '9')) && (c != '.')) )
          continue;
        
        try
        {
          fission_buildup_time = PhysicalUnits::stringToTimeDuration( field );
          if( fission_buildup_time < 1*PhysicalUnits::second )
            throw runtime_error( "Buildup must be 1 second or greator" );
          found_buildup = true;
          break;
        }catch( std::exception &e )
        {
          
        }
        
        if( ((i + 1) < remaining_fields.size()) )
        {
          try
          {
            fission_buildup_time = PhysicalUnits::stringToTimeDuration( field + remaining_fields[i+1] );
            if( fission_buildup_time < 1*PhysicalUnits::second )
              throw runtime_error( "Buildup must be 1 second or greator" );
            found_buildup = true;
            break;
          }catch( std::exception &e )
          {
            
          }
        }//if( !found_buildup )
      }//for( size_t i = 0; i < remaining_fields.size(); ++i )
      
      if( !found_buildup )
      {
        cerr << "Failed to interpret buildup time in '" << input.m_input_txt
        << "', - using 1 day." << endl;
        
        fission_buildup_time = 24*2600*PhysicalUnits::second;
        input.m_input_txt = "Fission " + fission_nuclide->symbol + " ";
        switch( fission_type )
        {
          case FissionType::Thermal:     input.m_input_txt += "Thermal"; break;
          case FissionType::FourteenMeV: input.m_input_txt += "14 MeV";  break;
          case FissionType::Fast:        input.m_input_txt += "Fast";    break;
        }//switch( fission_type )
        
        input.m_input_txt += " 1d buildup";
        answer.m_input_warnings.push_back( "Could not interpret buildup time; using 1 day." );
      }//try / catch to
    }//if( fission_nuclide && (fields.size() > 3) )
    
    if( fission_nuclide )
    {
      // We'll do a quick check that we have fission data file for this nuclide - we could still run
      //  into problems later.
      std::lock_guard<std::mutex> lock( sm_nuc_mix_mutex );
        
      if( sm_fission_products )
      {
        for( size_t i = 0; !fission_src && (i < sm_fission_products->size()); ++i )
          fission_src = (fission_nuclide == (*sm_fission_products)[i].nuclide);
      }//if( sm_fission_products )
      
      if( input.m_age == "" )
      {
        input.m_age = "1 day";
        age = 24.0*3600.0*PhysicalUnits::second;
      }else
      {
        try
        {
          age = PhysicalUnitsLocalized::stringToTimeDuration( input.m_age );
        }catch( std::exception & )
        {
          answer.m_input_warnings.push_back( "Invalid age input." );
          answer.m_validity = ReferenceLineInfo::InputValidity::InvalidAge;
          answer_ptr->m_input = input;
          return answer_ptr;
        }//try /catch to get the age
      }//if( input.m_age == "" ) / else
    }//if( fission_nuclide )
  }//if( not anything else so far - so try fission products )
  
  
  // If we couldnt identify the source, set error and return
  if( !nuc && !el && rctn_gammas.empty() && !is_background && !is_custom_energy && !nuc_mix && !src_lines && !fission_src )
  {
    answer.m_validity = ReferenceLineInfo::InputValidity::InvalidSource;
    answer.m_input_warnings.push_back( "'" + input.m_input_txt + "' is not a valid isotope, element, reaction, or energy." );
    
    answer_ptr->m_input = input;
    
    return answer_ptr;
  }//if( we couldnt match input text to a source )
  
  answer_ptr->m_validity = ReferenceLineInfo::InputValidity::Valid;
  answer_ptr->m_input = input;
  
  
  if( nuc )
    answer.m_source_type = ReferenceLineInfo::SourceType::Nuclide;
  else if( el )
    answer.m_source_type = ReferenceLineInfo::SourceType::FluorescenceXray;
  else if( !rctn_gammas.empty() )
    answer.m_source_type = ReferenceLineInfo::SourceType::Reaction;
  else if( is_background )
    answer.m_source_type = ReferenceLineInfo::SourceType::Background;
  else if( is_custom_energy )
    answer.m_source_type = ReferenceLineInfo::SourceType::CustomEnergy;
  else if( nuc_mix )
    answer.m_source_type = ReferenceLineInfo::SourceType::NuclideMixture;
  else if( src_lines )
    answer.m_source_type = ReferenceLineInfo::SourceType::OneOffSrcLines;
  else if( fission_src )
    answer.m_source_type = ReferenceLineInfo::SourceType::FissionRefLines;
  else
    answer.m_source_type = ReferenceLineInfo::SourceType::None;
  
  bool use_particle[SandiaDecay::ProductType::XrayParticle + 1] = {false};
  
  // We'll loop + switch over SandiaDecay::ProductType so we'll at least get
  //  a compiler warning if SandiaDecay::ProductType changes.
  for( auto type = SandiaDecay::ProductType(0);
      type <= SandiaDecay::ProductType::XrayParticle;
      type = SandiaDecay::ProductType(type+1) )
  {
    switch( type )
    {
      case SandiaDecay::ProductType::BetaParticle:
        use_particle[type] = input.m_showBetas;
        break;
        
      case SandiaDecay::ProductType::GammaParticle:
        use_particle[type] = (input.m_showGammas || input.m_showCascades || input.m_showEscapes);
        break;
        
      case SandiaDecay::ProductType::AlphaParticle:
        use_particle[type] = input.m_showAlphas;
        break;
        
      case SandiaDecay::ProductType::PositronParticle:
        use_particle[type] = input.m_showGammas;
        break;
        
      case SandiaDecay::ProductType::CaptureElectronParticle:
        break;
        
      case SandiaDecay::ProductType::XrayParticle:
        use_particle[type] = (el || input.m_showXrays);
        break;
    }//switch(type)
  }//for( loop over SandiaDecay::ProductType )
  
  
  vector<ReferenceLineInfo::RefLine> lines;
  
  //transition, first gamma BR, first gamma energy, second gamma energy, coincidence fraction, second gamma BR (just for debug)
  typedef tuple<const SandiaDecay::Transition *, double, float, float, float, double> coincidence_info_t;
  vector<coincidence_info_t> gamma_coincidences;
  
  /** The `m_decay_intensity` value returned are normalized to 1 bq of the parent nuclide, at the specified age.
   */
  auto make_nuc_lines = [use_particle,lower_photon_energy]( const SandiaDecay::Nuclide *nuc, double age, const RefLineInput &input )
              -> tuple<vector<ReferenceLineInfo::RefLine>,vector<coincidence_info_t>,bool> {
    
                
    vector<ReferenceLineInfo::RefLine> nuc_lines;
    vector<coincidence_info_t> nuc_gamma_coincidences;
    bool has_coincidences = false;
                
    SandiaDecay::NuclideMixture mixture;
    
    if( input.m_promptLinesOnly )
    {
      age = 0.0;
      mixture.addNuclideInPromptEquilibrium( nuc, 1.0E-3 * SandiaDecay::curie );
    } else
    {
      mixture.addNuclideByActivity( nuc, 1.0E-3 * SandiaDecay::curie );
    }//if( we want prompt only ) / else
    
    
    const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( age );
    
    // x-rays are slightly problematic - we can *almost* treat them like gammas, but for
    //  some decays they essentially get duplicated - so instead we'll be a little inefficient
    //  and track them separately.
    vector<SandiaDecay::EnergyRatePair> xrays = use_particle[SandiaDecay::ProductType::XrayParticle]
                                                ? mixture.xrays( age )
                                                : vector<SandiaDecay::EnergyRatePair>{};
    
    const double parent_activity = nuc ? mixture.activity( age, nuc ) : 0.0;
    
    // We will accumulate positrons as a single line, and just assign the transition
    //  as the first one we run into
    //  TODO: should make sure the first transition we run into is the heaviest one.
    ReferenceLineInfo::RefLine positron_line;
    positron_line.m_energy = 510.9989 * PhysicalUnits::keV;
    positron_line.m_decay_intensity = 0.0;
    positron_line.m_parent_nuclide = nuc;
    positron_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
    positron_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Annihilation;
    
    
    for( const SandiaDecay::NuclideActivityPair &nap : activities )
    {
      const SandiaDecay::Nuclide *nuclide = nap.nuclide;
      const double activity = nap.activity;
      
      for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
      {
        for( const SandiaDecay::RadParticle &particle : transition->products )
        {
          assert( particle.type <= SandiaDecay::ProductType::XrayParticle );
          
          if( !use_particle[particle.type] )
            continue;
          
          if( ((particle.type == SandiaDecay::GammaParticle)
               || (particle.type == SandiaDecay::XrayParticle))
             && (particle.energy < lower_photon_energy) )
          {
            continue;
          }
          
          if( particle.type == SandiaDecay::PositronParticle )
          {
            if( !positron_line.m_transition )
              positron_line.m_transition = transition;
            
            const double br = activity * particle.intensity
            * transition->branchRatio / parent_activity;
            
            positron_line.m_decay_intensity += 2.0 * br;
            
            continue;
          }//if( particle.type == SandiaDecay::PositronParticle )
          
          
          ReferenceLineInfo::RefLine line;
          line.m_parent_nuclide = nuc;
          line.m_energy = particle.energy;
          line.m_transition = transition;
          
          if( particle.type == SandiaDecay::XrayParticle )
          {
            size_t index = 0;
            for( ; index < xrays.size(); ++index )
            {
              if( fabs( xrays[index].energy - particle.energy ) < 1.0E-6 )
                break;
            }
            
            if( index < xrays.size() )
            {
              line.m_decay_intensity = xrays[index].numPerSecond / parent_activity;
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
              line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
              
              nuc_lines.push_back( line );
              
              // Erase this x-ray so we dont double-count it
              xrays.erase( begin( xrays ) + index );
            } else
            {
              // We've already accounted for this energy.
            }
            
            continue;
          }//if( particle.type == SandiaDecay::XrayParticle )
          
          
          if( !has_coincidences && (particle.type == SandiaDecay::GammaParticle) )
            has_coincidences = !particle.coincidences.empty();
          
          const double br = activity * particle.intensity * transition->branchRatio / parent_activity;
          
          if( input.m_showCascades && (particle.type == SandiaDecay::GammaParticle) )
          {
            for( size_t coinc_index = 0; coinc_index < particle.coincidences.size(); ++coinc_index )
            {
              const unsigned short int part_ind = particle.coincidences[coinc_index].first;
              const float fraction = particle.coincidences[coinc_index].second;
              assert( part_ind < transition->products.size() );
              if( part_ind < transition->products.size() )
              {
                const SandiaDecay::RadParticle &coinc_part = transition->products[part_ind];
                
                // The BR of second gamma is just for debugging
                const double second_br = activity * coinc_part.intensity
                * transition->branchRatio / parent_activity;
                
                if( coinc_part.type == SandiaDecay::ProductType::GammaParticle )
                  nuc_gamma_coincidences.emplace_back( transition, br, particle.energy, coinc_part.energy, fraction, second_br );
              }//if (part_ind < transition->products.size())
            }//for( loop over coincidences )
          }//if( show cascade gammas )
          
          line.m_decay_intensity = br;
          line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
          
          switch( particle.type )
          {
            case SandiaDecay::ProductType::BetaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Beta;
              break;
              
            case SandiaDecay::ProductType::GammaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
              break;
              
            case SandiaDecay::ProductType::AlphaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Alpha;
              break;
              
            case SandiaDecay::ProductType::PositronParticle:
            case SandiaDecay::ProductType::CaptureElectronParticle:
            case SandiaDecay::ProductType::XrayParticle:
              assert( 0 );
              continue;
              break;
          }//switch( particle.type )
          
          if( input.m_showEscapes
             && (particle.type == SandiaDecay::GammaParticle)
             && (line.m_energy > 1022.0) && (br > 0.0) )
          {
            ReferenceLineInfo::RefLine se_line = line;
            ReferenceLineInfo::RefLine de_line = line;
            
            se_line.m_energy -= (510.9989 * PhysicalUnits::keV);
            de_line.m_energy -= (2.0 * 510.9989 * PhysicalUnits::keV);
            se_line.m_decay_intensity = br * ns_single_escape_sf( line.m_energy );
            de_line.m_decay_intensity = br * ns_double_escape_sf( line.m_energy );
            se_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::SingleEscape;
            de_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape;
            
            char descTxt[128] = { '\0' };
            snprintf( descTxt, sizeof(descTxt), "S.E. of %.1f keV (generic HPGe eff.)", line.m_energy );
            se_line.m_decaystr = descTxt + ((se_line.m_decaystr.empty() ? "" : " ") + se_line.m_decaystr );
            
            snprintf( descTxt, sizeof(descTxt), "D.E. of %.1f keV (generic HPGe eff.)", line.m_energy );
            de_line.m_decaystr = descTxt + ((de_line.m_decaystr.empty() ? "" : " ") + de_line.m_decaystr );
            
            nuc_lines.push_back( se_line );
            nuc_lines.push_back( de_line );
          }//if( input.m_showEscapes && (particle.type == SandiaDecay::GammaParticle) )
          
          
          // If type is GammaParticle, we could be here if user selected to show cascade, but
          //  not actual gammas, so we need to check for this, and if so not add this gamma
          //  in to be shown
          if( (particle.type == SandiaDecay::GammaParticle) && !input.m_showGammas )
            continue;
          
          nuc_lines.push_back( line );
        }//for( const SandiaDecay::RadParticle &particle : transition->products )
      }//for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
    }//for( const SandiaDecay::NuclideActivityPair &nap : activities )
    
    if( positron_line.m_decay_intensity > 0.0 )
      nuc_lines.push_back( positron_line );
                
    return { nuc_lines, nuc_gamma_coincidences, has_coincidences };
  };//make_nuc_lines lambda
  
  if( nuc )
  {
    tuple<vector<ReferenceLineInfo::RefLine>,vector<coincidence_info_t>,bool> nuc_line_info
                                                                = make_nuc_lines( nuc, age, input );
    
    lines = std::move( std::get<0>(nuc_line_info) );
    gamma_coincidences = std::move( std::get<1>(nuc_line_info) );
    answer.m_has_coincidences = std::get<2>(nuc_line_info);
  }//if( nuc )
  
  
  if( nuc_mix )
  {
    for( const ReferenceLinePredef::NucMixComp &comp : nuc_mix->m_components )
    {
      assert( comp.m_nuclide );
      const double nuc_age = std::max( 0.0, age - comp.m_age_offset ); //Make sure age is at least zero
      tuple<vector<ReferenceLineInfo::RefLine>,vector<coincidence_info_t>,bool> nuc_line_info
                                = make_nuc_lines( comp.m_nuclide, nuc_age, input );
      
      vector<ReferenceLineInfo::RefLine> &these_lines = get<0>(nuc_line_info);
      vector<coincidence_info_t> &these_coinc = get<1>(nuc_line_info);
      
      // Correct BR of line to account for fraction of parent nuclide
      //
      // The `m_decay_intensity` value returned are normalized to 1 bq of the parent nuclide, at
      //  the specified age, so if we are having fixed activity fraction, we just need to multiply
      //  by the activity fractions specified.  But otherwise, we need to figure out the parent
      //  nuclides activity at the time being used, and then multiply by that.
      double correction_factor = comp.m_rel_act;
      if( !nuc_mix->m_fixed_act_fractions )
      {
        // `NucMixComp::m_rel_act` is for age=0 in this case
        correction_factor *= exp( -age * comp.m_nuclide->decayConstant() );
      }//if( nuc_mix->m_fixed_act_fractions )
      
      //cout << "For " << nuc_mix->m_name << ", the rel act. of " << comp.m_nuclide->symbol
      //     << " is " << correction_factor << " at age " << PhysicalUnits::printToBestTimeUnits(age)
      //     << " (m_rel_act=" << comp.m_rel_act
      //     << ", m_fixed_act_fractions=" << (nuc_mix->m_fixed_act_fractions ? "yes" : "no") << ")"
      //     << endl;
      
      for( ReferenceLineInfo::RefLine &this_line : these_lines )
        this_line.m_decay_intensity *= correction_factor;
      
      for( coincidence_info_t &this_coinc : these_coinc )
      {
        get<1>(this_coinc) *= correction_factor; // first gamma BR
        get<5>(this_coinc) *= correction_factor; // second gamma BR (just for debug)
      }
      
      if( !comp.m_color.isDefault() )
      {
        for( ReferenceLineInfo::RefLine &this_line : these_lines )
          this_line.m_color = comp.m_color;
      }
      
      lines.insert( end(lines), begin(these_lines), end(these_lines) );
      gamma_coincidences.insert( end(gamma_coincidences), begin(these_coinc), end(these_coinc) );
      
      //combine 511 lines
      //positron_line.m_energy = 510.9989 * PhysicalUnits::keV;
      //positron_line.m_parent_nuclide = nuc;
      //positron_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      //positron_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Annihilation;
      
      answer.m_has_coincidences |= std::get<2>(nuc_line_info);
    }//for( const NucMixComp &comp : nuc_mix->m_components )
      
    // We will limit the total number of gammas/x-rays, according to intensity
    // 1250 is arbitrarily chosen, but its hard to imagine more than this being relevant.
    //  In principle, we would like to do this after combining similar energies, but maybe
    //  this is fine for this kinda rare use case.
    // Pu mixtures will have like 5000 lines.
    const size_t max_num_gamma_lines = 1250;
    if( lines.size() > max_num_gamma_lines )
    {
      std::sort( begin(lines), end(lines), []( const ReferenceLineInfo::RefLine &lhs,
                                              const ReferenceLineInfo::RefLine &rhs ) -> bool {
        return lhs.m_decay_intensity > rhs.m_decay_intensity;
      } );
      
      // cout << "Erasing " << (lines.size() - max_num_gamma_lines) << " lines, starting with amp: "
      //  << lines[max_num_gamma_lines].m_decay_intensity << endl;
      lines.erase( begin(lines) + max_num_gamma_lines, end(lines) );
    }//if( lines.size() > 1250 )
  }//if( nuc_mix )
  
  if( src_lines )
  {
    answer.m_has_coincidences = false;
    
    for( const ReferenceLinePredef::CustomLine &line : src_lines->m_lines )
    {
      ReferenceLineInfo::RefLine refline;
      refline.m_energy = line.m_energy;
      refline.m_normalized_intensity = line.m_branch_ratio / src_lines->m_max_branch_ratio;
      refline.m_particle_sf_applied = 1.0f / src_lines->m_max_branch_ratio;
      refline.m_decaystr = line.m_info;
      refline.m_decay_intensity = line.m_branch_ratio;
      refline.m_parent_nuclide = line.m_nuclide;
      refline.m_transition = line.m_transition;
      refline.m_color = line.m_color;
      refline.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      refline.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      refline.m_attenuation_applies = line.m_atten_applies;
      
      lines.push_back( std::move(refline) );
    }//for( const CustomLine &line : src_lines->m_lines )
  }//if( src_lines )
  
  // Update showing cascades based on if there are actually any present
  input.m_showCascades = (answer.m_has_coincidences && input.m_showCascades);
  answer.m_input.m_showCascades = (answer.m_has_coincidences && input.m_showCascades);
  
  
  
  if( el )
  {
    for( const SandiaDecay::EnergyIntensityPair &eip : el->xrays )
    {
      if( eip.energy < lower_photon_energy )
        continue;
      
      ReferenceLineInfo::RefLine line;
      line.m_element = el;
      line.m_energy = eip.energy;
      line.m_decay_intensity = eip.intensity;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      
      lines.push_back( line );
    }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
  }//if( m_showXrays->isChecked() )
  
  
  if( !rctn_gammas.empty() )
  {
    for( const ReactionGamma::ReactionPhotopeak &eip : rctn_gammas )
    {
      if( eip.reaction )
        answer.m_reactions.insert( eip.reaction );
      
      if( eip.energy < lower_photon_energy )
        continue;
      
      ReferenceLineInfo::RefLine line;
      line.m_reaction = eip.reaction;
      line.m_energy = eip.energy;
      line.m_decay_intensity = eip.abundance;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      line.m_decaystr = eip.remark;
      
      lines.push_back( line );
    }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
  }//if( !rctn_gammas.empty() )
  
  
  if( !otherRefLinesToShow.empty() )
  {
    assert( is_background || is_custom_energy );
    assert( !nuc && !el && rctn_gammas.empty() );
    
    const SandiaDecay::Nuclide *u238 = is_background ? db->nuclide( "U238" ) : nullptr;
    const SandiaDecay::Nuclide *u235 = is_background ? db->nuclide( "U235" ) : nullptr;
    const SandiaDecay::Nuclide *th232 = is_background ? db->nuclide( "Th232" ) : nullptr;
    const SandiaDecay::Nuclide *ra226 = is_background ? db->nuclide( "Ra226" ) : nullptr;
    const SandiaDecay::Nuclide *k40 = is_background ? db->nuclide( "K40" ) : nullptr;
    
    for( const OtherRefLine &bl : otherRefLinesToShow )
    {
      if( std::get<0>( bl ) < lower_photon_energy )
        continue;
      
      ReferenceLineInfo::RefLine line;
      line.m_energy = std::get<0>( bl );
      line.m_decay_intensity = std::get<1>( bl );
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      
      if( !std::get<2>(bl).empty() )
        line.m_decaystr = std::get<2>(bl) + ", ";
      
      switch( std::get<3>( bl ) )
      {
        case OtherRefLineType::U238Series:
          line.m_parent_nuclide = u238;
          line.m_decaystr += "U238 series";
          // TODO: set a color here
          break;
          
        case OtherRefLineType::U235Series:
          line.m_parent_nuclide = u235;
          line.m_decaystr += "U235 series";
          // TODO: set a color here
          break;
          
        case OtherRefLineType::Th232Series:
          line.m_parent_nuclide = th232;
          line.m_decaystr += "Th232 series";
          // TODO: set a color here
          break;
          
        case OtherRefLineType::Ra226Series:
          line.m_parent_nuclide = ra226;
          line.m_decaystr += "U238 (Ra226) series";
          // TODO: set a color here
          break;
          
        case OtherRefLineType::K40Background:
          line.m_parent_nuclide = k40;
          line.m_decaystr += "Primordial";
          // TODO: set a color here
          break;
          
        case OtherRefLineType::BackgroundXRay:
        {
          //std::get<2>( bl ) will be like "Pb xray"
          line.m_decaystr = std::get<2>( bl );
          if( !std::get<4>( bl ).empty() )
            line.m_decaystr += (line.m_decaystr.empty() ? "" : ", ") + std::get<4>( bl );
          line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
          vector<string> parts;
          SpecUtils::split( parts, std::get<2>( bl ), " " );
          if( !parts.empty() )
            line.m_element = db->element( parts[0] );
          // TODO: set a color here
          break;
        }//case OtherRefLineType::BackgroundXRay:
          
        case OtherRefLineType::OtherBackground:
        case OtherRefLineType::BackgroundReaction:
        {
          const string &nucstr = std::get<2>( bl );
          
          line.m_decaystr = nucstr;
          
          line.m_parent_nuclide = db->nuclide( nucstr );
          
          //A string like "Th232 S.E. 2614 keV" wont return a valid nuclide, so we'll fix this up
          const size_t single_escape_pos = nucstr.find("S.E.");
          const size_t double_escape_pos = nucstr.find("D.E.");
          
          if( !line.m_parent_nuclide )
          {
            const size_t pos = std::min( single_escape_pos, double_escape_pos );
            if( pos != string::npos )
              line.m_parent_nuclide = db->nuclide( nucstr.substr(0,pos) ); // Make like "Th232"
          }//if( !line.m_parent_nuclide )
          
          
          // TODO: try to get reaction if didnt get nuclide - also, nuclide list may be CSV, could split that
          //const ReactionGamma *rctnDb = ReactionGammaServer::database();
          //if( rctnDb )
          //{
          //  rctnDb->gammas( input.m_input_txt, rctn_gammas );
          //...
          //}
          
          if( !std::get<4>( bl ).empty() )
            line.m_decaystr += (line.m_decaystr.empty() ? "" : ", ") + std::get<4>( bl );
          
          // TODO: We should probably do a little better than string matching to define as a single/double escape peak.
          if( single_escape_pos != string::npos )
          {
            line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::SingleEscape;
          }else if( double_escape_pos != string::npos )
          {
            line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape;
          }
          
          // TODO: set a color here
          break;
        }
      }//switch( get<3>(*bl) )
      
      lines.push_back( line );
    }//for( otherRefLinesToShow )
  }//if( !otherRefLinesToShow.empty() )
  
  
  assert( !fission_src || ((!fission_src) == (!fission_nuclide)) );
  if( fission_src && fission_nuclide )
  {
    // At least for the moment, we'll avoid trying to deal with coincidences.
    answer.m_has_coincidences = false;
    
    try
    {
      const shared_ptr<const vector<FissionLine>> fission_lines
               = fission_photons( fission_nuclide, fission_type, fission_buildup_time, age );
      
      assert( fission_lines );
      if( !fission_lines )
        throw logic_error( "Unexpected nullptr from fission_photons." );
      
      double max_amp = 0.0;
      for( const FissionLine &info : *fission_lines )
        max_amp = std::max( max_amp, info.relative_amplitude );
      assert( !IsNan(max_amp) && !IsInf(max_amp) && (max_amp > 0.0) );
      if( IsNan(max_amp) || IsInf(max_amp) || (max_amp <= 0.0) )
        max_amp = 1.0; //good luck
      
      for( const FissionLine &info : *fission_lines )
      {
        assert( !IsNan(info.energy) && !IsInf(info.energy) );
        assert( !IsNan(info.relative_amplitude) && !IsInf(info.relative_amplitude) );
        if( IsNan(info.relative_amplitude) || IsInf(info.relative_amplitude)
           || (info.relative_amplitude <= 0.0) )
          continue;
        
        ReferenceLineInfo::RefLine line;
        line.m_energy = info.energy;
        line.m_decay_intensity = info.relative_amplitude / max_amp;
        line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
        line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
        line.m_parent_nuclide = info.parent;
        line.m_transition = info.transition;
        
        lines.push_back( std::move(line) );
      }//for( const FissionLine &info : *fission_lines )
    }catch( std::exception &e )
    {
      answer.m_input_warnings.push_back( "Error computing fission products: " + std::string(e.what()) );
      answer.m_validity = ReferenceLineInfo::InputValidity::InvalidSource;
      return answer_ptr;
    }//try / catch to
    
  }//if( fission_src )
  
  // Now calc detector response and shielding
  //  Up to now, we shouldnt have any escape or sum gammas in answer.m_ref_lines
  double max_alpha_br = 0.0, max_beta_br = 0.0, max_photon_br = 0.0;
  for( ReferenceLineInfo::RefLine &line : lines )
  {
    assert( (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal)
           || (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Annihilation)
           || (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::SingleEscape)
           || (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape) );
    
    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:
        max_alpha_br = std::max( max_alpha_br, line.m_decay_intensity );
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Beta:
        max_beta_br = std::max( max_beta_br, line.m_decay_intensity );
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Gamma:
      case ReferenceLineInfo::RefLine::Particle::Xray:
      {
        double energy = line.m_energy;
        switch( line.m_source_type )
        {
          case ReferenceLineInfo::RefLine::RefGammaType::Normal:
          case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
            // TODO: Need to put in S.E. DRF factor here
            energy += 510.998950;
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
            // TODO: Need to put in D.E. DRF factor here
            energy += 2.0*510.998950;
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
          case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
            assert( 0 );
            break;
        }//switch( line.m_source_type )
        
        if( input.m_det_intrinsic_eff && line.m_attenuation_applies )
          line.m_drf_factor = input.m_det_intrinsic_eff( energy );
        
        if( input.m_shielding_att && line.m_attenuation_applies )
          line.m_shield_atten = input.m_shielding_att( energy );
        
        max_photon_br = std::max( max_photon_br,
                                 line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten );
        break;
      }
    }//switch( line.m_particle_type )
  }//for( ReferenceLineInfo::RefLine &line : lines )
  
  
  const double alpha_sf = ((max_alpha_br > 0.0) && !IsNan(max_alpha_br)) ? (1.0 / max_alpha_br) : 1.0;
  const double beta_sf = ((max_beta_br > 0.0) && !IsNan(max_beta_br)) ? (1.0 / max_beta_br) : 1.0;
  const double photon_sf = ((max_photon_br > 0.0) && !IsNan(max_photon_br)) ? (1.0 / max_photon_br) : 1.0;
  
  for( ReferenceLineInfo::RefLine &line : lines )
  {
    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:
        line.m_particle_sf_applied = alpha_sf;
        line.m_normalized_intensity = line.m_decay_intensity * alpha_sf;
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Beta:
        line.m_particle_sf_applied = beta_sf;
        line.m_normalized_intensity = line.m_decay_intensity * beta_sf;
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Gamma:
      case ReferenceLineInfo::RefLine::Particle::Xray:
        line.m_particle_sf_applied = photon_sf;
        line.m_normalized_intensity = photon_sf * line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten;
        break;
    }//switch( line.m_particle_type )
    
    // We wont filter out lines smaller than wanted here
    if( (line.m_transition || is_background)
       && (line.m_normalized_intensity <= input.m_lower_br_cutt_off
           || IsInf( line.m_normalized_intensity )
           || IsNan( line.m_normalized_intensity )) )
    {
      continue;
    }
    
    
    // Now lets fill out line.m_decaystr
    if( line.m_reaction )
    {
      // TODO: we can probably come up with a better way to describe reactions
      if( line.m_reaction->targetNuclide && line.m_reaction->targetElement )
      {
        const string orig_deacystr = line.m_decaystr;
        switch( line.m_reaction->type )
        {
          case AlphaNeutron:   line.m_decaystr = "Alphas on "; break;
          case NeutronAlpha:   line.m_decaystr = "Neutrons on "; break;
          case AlphaProton:    line.m_decaystr = "Alphas on "; break;
          case NeutronCapture: line.m_decaystr = "Neutron capture by "; break;
          case NeutronInelasticScatter: line.m_decaystr = "Neutron inelastic scatter on "; break;
          case AlphaInelasticScatter:   line.m_decaystr = "Alpha inelastic scatter on "; break;
          case AnnihilationReaction:    line.m_decaystr = ""; break;
          case NumReactionType:         line.m_decaystr = ""; break;
        }//switch( line.m_reaction->type )
        
        line.m_decaystr += line.m_reaction->targetNuclide->symbol;
        
        if( !orig_deacystr.empty() )
          line.m_decaystr += ". " + orig_deacystr;
        
        //if( line.m_reaction->productNuclide )
        //  line.m_decaystr += " to give " + line.m_reaction->productNuclide->symbol;
      }//
    }else if( line.m_decaystr.empty() )
    {
      if( line.m_transition )
      {
        if( line.m_transition->parent )
          line.m_decaystr = line.m_transition->parent->symbol;
        if( line.m_transition->child )
          line.m_decaystr += " to " + line.m_transition->child->symbol;
        
        // TODO: for alphas and betas its pretty redundant to have this next line (I guess its redundant no matter what actually)
        line.m_decaystr += string( " via " ) + SandiaDecay::to_str( line.m_transition->mode );
        
        switch( line.m_transition->mode )
        {
          case SandiaDecay::AlphaDecay:
          case SandiaDecay::BetaDecay:
          case SandiaDecay::BetaPlusDecay:
          case SandiaDecay::ProtonDecay:
            line.m_decaystr += " decay";
            break;
            
          case SandiaDecay::IsometricTransitionDecay:
            line.m_decaystr += " transition";
            break;
          
          default:
            break;
        }//switch( line.m_transition->mode )
      }else if( line.m_element )
      {
        line.m_decaystr = line.m_element->name + " fluorescence";
      }
    }//if( line.m_decaystr.empty() )
    
    const double &amp = line.m_normalized_intensity;
    if( !IsNan( amp ) && !IsInf( amp )
       && (amp >= std::numeric_limits<float>::min()) // numeric_limits<float>::min()==1.17549e-38
       && (amp > input.m_lower_br_cutt_off)
       )
    {
      answer.m_ref_lines.push_back( line );
    }
  }//for( ReferenceLineInfo::RefLine &line : answer.m_ref_lines )
  
  
  // If we add in escape peaks - we could put them in here
  
  // Add in coincident gammas
  if( !gamma_coincidences.empty() )
  {
    double max_coincidence_br = 0.0;
    vector<ReferenceLineInfo::RefLine> coinc_ref_lines;
    for( const auto &casc : gamma_coincidences )
    {
      const SandiaDecay::Transition *const &trans = get<0>( casc );
      const double &first_br = get<1>( casc );
      const float &first_energy = get<2>( casc );
      const float &second_energy = get<3>( casc );
      const float &coinc_frac = get<4>( casc );
      const double &second_br = get<5>( casc );
      
      const float energy = first_energy + second_energy;
      
      ReferenceLineInfo::RefLine line;
      line.m_energy = energy;
      line.m_decay_intensity = first_br * coinc_frac;
      line.m_parent_nuclide = nuc;
      line.m_transition = trans;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak;
      
      if( input.m_det_intrinsic_eff )
        line.m_drf_factor = input.m_det_intrinsic_eff( first_energy ) * input.m_det_intrinsic_eff( second_energy );
      
      if( input.m_shielding_att )
        line.m_shield_atten = input.m_shielding_att( first_energy ) * input.m_shielding_att( second_energy );
      
      const double amp = line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten;
      assert( !IsNan( amp ) && !IsInf( amp ) );
      if( IsNan( amp ) || IsInf( amp ) )
      {
        cerr << "Unexpected NaN or Inf coincidence amp." << endl;
        continue;
      }
      
      line.m_decaystr = "Cascade sum";
      if( trans && trans->parent )
        line.m_decaystr += " " + trans->parent->symbol;
      if( trans && trans->child )
        line.m_decaystr += " to " + trans->child->symbol;
      
      char buffer[128];
      snprintf( buffer, sizeof( buffer ),
               " (%.1f + %.1f keV, coinc=%.3g)",
               first_energy, second_energy, coinc_frac );
      
      line.m_decaystr += buffer;
      
      coinc_ref_lines.push_back( std::move( line ) );
      
      max_coincidence_br = std::max( max_coincidence_br, amp );
    }//for( loop over cascades )
    
    assert( coinc_ref_lines.empty()
           || ((max_coincidence_br > 0.0) && !IsNan( max_coincidence_br )) );
    
    // Scale the coincidence line amplitudes to be between 0
    for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    {
      const double sf = 1.0 / max_coincidence_br;
      line.m_particle_sf_applied = sf;
      const double amp = line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten * sf;
      line.m_normalized_intensity = amp;
    }//for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    
    // There can be tons of cascade sums (4834 for U238), we'll limit the number
    //   we draw to an arbitrary 350, because this is even more than I expect to
    //   be relevant (although I didnt actually check this).
    //  TODO: limit based on importance, and not a flat limit, e.g., use something like
    //        yield(i)*sqrt(energy(i))/sum(yield*sqrt(energy))
    const size_t max_cascade_sums = 350;
    if( coinc_ref_lines.size() > max_cascade_sums )
    {
      std::sort( begin( coinc_ref_lines ), end( coinc_ref_lines ),
                []( const ReferenceLineInfo::RefLine &lhs, const ReferenceLineInfo::RefLine &rhs ) -> bool {
        if( lhs.m_normalized_intensity == rhs.m_normalized_intensity )
          return lhs.m_energy > rhs.m_energy;
        return lhs.m_normalized_intensity > rhs.m_normalized_intensity;
      } );
      
      cout << "Resizing cascade sums from " << coinc_ref_lines.size() << " to " << max_cascade_sums << endl;
      coinc_ref_lines.resize( max_cascade_sums );
    }//if( coinc_ref_lines.size() > 350 )
    
    answer.m_ref_lines.reserve( answer.m_ref_lines.size() + coinc_ref_lines.size() );
    
    for( const ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    {
      const double &amp = line.m_normalized_intensity;
      if( !IsNan( amp ) && !IsInf( amp )
         && (amp >= std::numeric_limits<float>::min()) // numeric_limits<float>::min()==1.17549e-38
         && (amp > input.m_lower_br_cutt_off)
         )
        answer.m_ref_lines.push_back( line );
    }//for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
  }//if( !gamma_coincidences.empty() )
  
  
  //Client-side javascript currently doesn't know about this guarantee that gamma
  //  lines will be sorted by energy.
  answer.sortByEnergy();
  
  return answer_ptr;
}//std::shared_ptr<ReferenceLineInfo> generateRefLineInfo()


vector<string> ReferenceLineInfo::additional_ref_line_sources()
{
  vector<string> answer;
  
  std::lock_guard<std::mutex> lock( sm_nuc_mix_mutex );
  
  if( !sm_have_tried_init )
    load_custom_nuc_mixes();
  
  if( sm_nuc_mixes )
  {
    for( const auto &n : *sm_nuc_mixes )
      answer.push_back( n.second.m_name );
  }//if( sm_nuc_mixes )
  
  if( sm_custom_lines )
  {
    for( const auto &n : *sm_custom_lines )
      answer.push_back( n.second.m_name );
  }//if( sm_custom_lines )
  
  if( sm_fission_products )
  {
    for( const FissionDataSrcFile &src : *sm_fission_products )
    {
      assert( src.nuclide );
      
      // U238 doesnt have any thermal fission - so lets not suggest it (probably other nuclides too)
      const bool has_thermal = (src.nuclide->symbol != "U238");
      
      if( has_thermal )
        answer.push_back( "Fission " + src.nuclide->symbol + " Thermal, 1d buildup" );
      answer.push_back( "Fission " + src.nuclide->symbol + " Fast, 1d buildup" );
      answer.push_back( "Fission " + src.nuclide->symbol + " 14 MeV, 1d buildup" );
      
      if( has_thermal )
        answer.push_back( "Fission " + src.nuclide->symbol + " Thermal, 30d buildup" );
      answer.push_back( "Fission " + src.nuclide->symbol + " Fast, 30d buildup" );
      answer.push_back( "Fission " + src.nuclide->symbol + " 14MeV, 30d buildup" );
      
      if( has_thermal )
        answer.push_back( "Fission " + src.nuclide->symbol + " Thermal, 1y buildup" );
      answer.push_back( "Fission " + src.nuclide->symbol + " Fast, 1y buildup" );
      answer.push_back( "Fission " + src.nuclide->symbol + " 14 MeV, 1y buildup" );
    }//for( const FissionDataSrcFile &src : *sm_fission_products )
  }//if( sm_fission_products )
  
  return answer;
}//vector<string> additional_ref_line_sources()


void ReferenceLineInfo::load_nuclide_mixtures()
{
  std::lock_guard<std::mutex> lock( sm_nuc_mix_mutex );
  if( !sm_have_tried_init )
    load_custom_nuc_mixes();
}//void load_nuclide_mixtures()
