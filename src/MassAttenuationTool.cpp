/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <map>
#include <mutex>
#include <cstdio>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <stdexcept>

#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"

using namespace std;

namespace
{
  //using phrase_parse inst of split_to_floats adds about 6.1 kb to binary
  //  size on 64 bit OSX MinSizeRel
  bool split_space_delim_flts( const string &str, vector<float> &res )
  {
    namespace qi = boost::spirit::qi;
    
    res.clear();
    
    //using phrase_parse inst of split_to_floats adds about 6.1 kb to binary
    //  size on 64 bit OSX MinSizeRel
    //    bool success = UtilityFunctions::split_to_floats( line.c_str(), line.size(), energies );
    return qi::phrase_parse( str.c_str(), str.c_str()+str.size(),
                            (*qi::float_) % qi::eol, qi::space, res );
  }//split_space_delim_flts(...)
  
  void exact_short_float( const float val, char *buffer, const size_t bufferlen )
  {
    snprintf( buffer, bufferlen, "%.9g", val );
  }//void exact_short_float( const float val, char *buffer, const size_t bufferlen )
  
  
  //Taken from UtilityFunctions, but copied to reduce cross-dependancies (and
  //  use on other projects that dont otherwise need UtilityFunctions)
  std::string append_path( const std::string &base, const std::string &name )
  {
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
    if( base.size() && (base[base.size()-1]=='\\'||base[base.size()-1]=='/') )
      return base + name;
    if( name.size() && (name[0]=='\\'||name[0]=='/') )
      return base + name;
    return base + '\\' + name;
#else
    if( base.size() && base[base.size()-1]=='/' )
      return base + name;
    if( name.size() && name[0]=='/' )
      return base + name;
    return base + '/' + name;
#endif
  }//std::string append_path( const std::string &base, const std::string &name )
  
  
  struct AttCoeffData
  {
    float ZA[100];           //atomic number / atomic weight
    float ELIM[100][13];     //energy limits photoelectric
    float ACP[100][13][4];   //photoelectric cross section coefficients
    float U[100][8];         //pair production coefficients
    char Element[100][2];
    float AtomicWeight[100];      //not used
    float EnergyKEdge[100];       //not used
    float FluorescenceYield[100]; //not used
  };//struct AttCoeffData
  
  struct ElementProccessCoeffients
  {
    MassAttenuation::GammaEmProcces m_proccess;
    std::vector<float> m_logEnergies;
    std::vector<float> m_logAttenuationCoeffs;
    
    size_t memsize() const;
  };//struct ElementProccessCoeffients
  
  struct ElementAttenuation
  {
    std::string m_symbol;
    float m_atomicMass;  //in PhysicalUnits grams per mole
    int m_atomicNumber;
    ElementProccessCoeffients m_proccesses[static_cast<int>(MassAttenuation::GammaEmProcces::NumGammaEmProcces)];
    
    size_t memsize() const;
    
    void saveTxt( std::string datapath );
    void loadTxt( std::string datapath, const int atomicNumber );
  };//struct ElementAttenuation


  /**
   Valid for for atomic numbers 1 through (and including) 98.
   Class should be thread safe... but not comprehensively tested yet.
   */
  class MassAttenuationTool
  {
  public:
    /** Constructor for the MassAttenuationTool.
     *
     * \param datapath Path to where the files similar to 1.xs.txt,... 98.xs.txt
     *        live at.  Typically will be "data/em_xs_data/".
     */
    MassAttenuationTool();
    
    /** Must be called before any cross-sections are loaded or an exception will
        be thrown
     */
    void set_data_dir( const std::string &dir );
    
    //Protect against copy or copy construction
    MassAttenuationTool( const MassAttenuationTool &tool ) = delete;
    MassAttenuationTool &operator=( const MassAttenuationTool &rhs ) = delete;
    
    ~MassAttenuationTool();
    
    /** Gives the mass attenuation coefficient in units of PhysicalUnits, that
     * is to print out to familiar units you would divide by
     * (PhysicalUnits::cm2 / PhysicalUnits::g) to put the result into units
     * of cm2/g.
     * To calculate the
     *
     *  Will throw ErrorLoadingDataException if the XS data cannot be loaded.
     *  Will throw std::runtime_error if atomic_number or energy is out of range.
     *
     * \param atomic_number Atomic number ranging from 1 to 98, inclusive
     * \param energy Energy (in keV)
     * \returns The mass attenuation coefficient for:
     *          compton + pair production + photo electric.
     */
    float massAttenuationCoeficient( const int atomic_number, const float energy );
    
    /** Similar to the other #massAttenuationCoeficient function, but instead
     * only for a specific sub-proccess.
     */
    float massAttenuationCoeficient( const int atomic_number, const float energy, MassAttenuation::GammaEmProcces process );
    
    /** Similar to #massAttenuationCoeficient, but interpolated the result between
     * floor(atomic_number) and ceil(atomic_number) linearly.
     */
    float massAttenuationCoeficientFracAN( const float atomic_number, const float energy );
    
    
    static float logLogInterpolate( const float energy,
                                   const std::vector<float> &logenergy,
                                   const std::vector<float> &logxs );
    
    /** Gives approximatly how much memorry is being taken up by this object.
     * Gives ~722 kb on my 64 bit mac.
     * \returns approximate memory this object is taking up, in bytes.
     */
    size_t memsize() const;
    
  protected:
    
    /**
     */
    const ElementAttenuation *attenuationData( const int atomic_number );
    
    std::string m_dataPath;
    
    std::atomic<const ElementAttenuation *> m_atten[98];
  };//class MassAttenuationTool
  
  inline float calcMassAttenuationCoeficient( const float energy,
                                             const MassAttenuation::GammaEmProcces process,
                                             const ElementAttenuation * const data )
  {
    if( data->m_proccesses[static_cast<int>(process)].m_logEnergies.empty() ) //wont happen unless catasrophy
      throw runtime_error( "Not-loaded data" );
    
    return MassAttenuationTool::logLogInterpolate( energy,
                                                  data->m_proccesses[static_cast<int>(process)].m_logEnergies,
                                                  data->m_proccesses[static_cast<int>(process)].m_logAttenuationCoeffs );
  }//calcMassAttenuationCoeficient(...)
}//namespace





namespace MassAttenuation
{
  static std::mutex sm_data_directory_mutex;  //probably a bit overkill
  static bool sm_have_loadded_data = false;
  static std::string sm_data_directory = "data";
  
  /** Have a singleton MassAttenuationTool tool object
   */
  static MassAttenuationTool sm_xs_tool;
  
  std::string data_directory()
  {
    std::lock_guard<std::mutex> lock( sm_data_directory_mutex );
    sm_have_loadded_data = true;
    return sm_data_directory;
  }
  
  void set_data_directory( const std::string &dir )
  {
    std::lock_guard<std::mutex> lock( sm_data_directory_mutex );
    if( sm_have_loadded_data )
      throw runtime_error( "You must not call MassAttenuation::set_data_directory(...) after using any cross sections" );
    sm_data_directory = dir;
    sm_xs_tool.set_data_dir( append_path( dir, "em_xs_data") );
  }
  
  float AttCoef( const float energy, const int atomic_number,
                float &scatter_mu, float &photoelectric_mu, float &pair_prod_mu )
  {
    scatter_mu = photoelectric_mu = pair_prod_mu = 0.0f;
    
    const float B[5] = { 1.148f, 0.06141f, 3.171f, 0.9328f, 0.02572f };
    
    //We will read in CrossSection.lib only once we need it, using a static atomic
    //  pointer to it so we can be both thread safe and not have to use a mutex.
    static std::atomic<const AttCoeffData *> s_data( 0 );
    
    const AttCoeffData *data = s_data.load();  //could probably do boost::memory_order_acquire instead of boost::memory_order_seq_cst
    
    if( !data )
    {
      const string path = append_path( MassAttenuation::data_directory(), "CrossSection.lib" );
      
      ifstream input( path.c_str() );
      
      if( !input.is_open() )
        throw runtime_error( "Failed reading " + path );
      
      AttCoeffData *newdata = new AttCoeffData();
      
      if( !input.read(   (char *)newdata->ZA,   sizeof(newdata->ZA) )
         || !input.read( (char *)newdata->ELIM, sizeof(newdata->ELIM) )
         || !input.read( (char *)newdata->ACP,  sizeof(newdata->ACP) )
         || !input.read( (char *)newdata->U,    sizeof(newdata->U) )
         || !input.read( (char *)newdata->Element,           sizeof(newdata->Element) )
         || !input.read( (char *)newdata->AtomicWeight,      sizeof(newdata->AtomicWeight) )
         || !input.read( (char *)newdata->EnergyKEdge,       sizeof(newdata->EnergyKEdge) )
         || !input.read( (char *)newdata->FluorescenceYield, sizeof(newdata->FluorescenceYield) ) )
        throw runtime_error( "Error reading CrossSection.lib" );
      
      //Should put in a check that we are at the end of the file
      
      const bool changed = s_data.compare_exchange_strong( data, newdata );
      
      if( changed ) //newdata was stored into s_data (was still prev. null)
        data = newdata;
      else //Another thread updated s_data; delete the copy we made, and other
        delete newdata;
    }//if( !data )
    
    if( energy < 0.0f || energy > 1.E6f )
      throw runtime_error( "Invalid energy to AttCoef" );
    
    if( atomic_number < 1 || atomic_number > 98 )
      throw runtime_error( "Invalid atomic number to AttCoef" );
    
    const int NA = min( atomic_number, 94 ) - 1;
    
    const float X = energy / 511.006f;
    
    scatter_mu = 0.4006f*data->ZA[NA]*(1.0f + X*(B[0]+X*B[1])) / (1.0f + X*(B[2]+X*(B[3]+X*B[4])));
    
    //Photoelectric
    int j = 0;
    while( energy > data->ELIM[NA][j] && j < 12 )
      ++j;
    
    float POWE = 1.0f;
    for( int k = 0; k < 4; ++k )
    {
      POWE = POWE / energy;
      photoelectric_mu += data->ACP[NA][j][k]*POWE;
    }
    
    //Pair-production
    if( energy < 1023.0f )
      pair_prod_mu = 0.0f;
    else if( energy < 1500.0f )
      pair_prod_mu = data->U[NA][7]*(energy-1022.012f)*(energy-1022.012f);
    else
    {
      float v = sqrt( energy - 1022.012f );
      float POW = 1.0;
      pair_prod_mu = 0.;
      for( int k = 0; k < 7; ++k )
      {
        POW *= v;
        pair_prod_mu += data->U[NA][k]*POW;
      }
      pair_prod_mu /= (1.0f + POW*1.7E-13f / v);
    }
    
    //C	Get total cross sections in units of cm-1.
    //cerr << "photoelectric_mu=" << photoelectric_mu << endl;
    //cerr << "scatter_mu=" << scatter_mu << endl;
    //cerr << "pair_prod_mu=" << pair_prod_mu << endl;
    
    return photoelectric_mu + scatter_mu + pair_prod_mu;
  }//float AttCoef( const float energy, const float atomic_number )
  
  
  float massAttenuationCoeficient( const int atomic_number, const float energy )
  {
    return sm_xs_tool.massAttenuationCoeficient( atomic_number, energy );
  }
  
  float massAttenuationCoeficient( const int atomic_number, const float energy, MassAttenuation::GammaEmProcces process )
  {
    return sm_xs_tool.massAttenuationCoeficient( atomic_number, energy, process );
  }
  
  float massAttenuationCoeficientFracAN( const float atomic_number, const float energy )
  {
    return sm_xs_tool.massAttenuationCoeficientFracAN( atomic_number, energy );
  }
}//namespace MassAttenuation




namespace
{

float MassAttenuationTool::logLogInterpolate( const float energy,
                                              const vector<float> &logenergy,
                                              const vector<float> &logxs )
{
  const float log_x = log10(energy);
  vector<float>::const_iterator ebegin = logenergy.begin();
  vector<float>::const_iterator eend = logenergy.end();
  const vector<float>::const_iterator iter = lower_bound( ebegin, eend, log_x );

  //Note: the (iter == (eend-1)) test below excludes values exactly equal
  //      to the highest energy value in the data file, but this is a detail
  if( iter == eend || (iter == (eend-1)) || iter == ebegin )
  {
    //Note that choosing 5keV to 10 MeV is arbitrary, and I didnt actually check
    //  the valid range of the cross-section files, but I think this should be
    //  fine
    if( energy > 5*PhysicalUnits::keV && energy < 10.0*PhysicalUnits::MeV )
      return 0;
    throw runtime_error( "logLogInterpolatedValue(...): Out of range" );
  }//if( iter == eend || (iter == (eend-1)) || iter == ebegin )

  const size_t bin = iter - ebegin - 1;
  assert( logenergy[bin] <= log_x );
  assert( logenergy[bin+1] >= log_x );
  const float f =(log_x - logenergy[bin])/(logenergy[bin+1] - logenergy[bin]);
  const float value = logxs[bin] + (logxs[bin+1] - logxs[bin])*f;
  const float answer = pow(float(10.0),value);

  if( IsNan(answer) )
  {
    cerr << "Found nan for input energy " << energy/PhysicalUnits::keV
         << " keV and bin=" << bin << ", f=" << f << ", log_x=" << log_x
         << ", logenergy[bin]=" << logenergy[bin] << ", logenergy[bin+1]="
         << logenergy[bin+1] << ", logxs[bin]=" << logxs[bin]
         << ", logxs[bin+1]=" << logxs[bin+1] << endl;
    return 0.0;
  }//if( IsNan(answer) )

  return answer;
}//float logLogInterpolate(...)


size_t ElementProccessCoeffients::memsize() const
{
  return sizeof(*this)
         + m_logEnergies.size()*sizeof(float)
         + m_logAttenuationCoeffs.size()*sizeof(float);
}//size_t memsize() const


size_t ElementAttenuation::memsize() const
{
  size_t size = sizeof(*this);
  size += m_symbol.capacity()*sizeof(std::string::value_type);
  
  for( const auto &p : m_proccesses )
    size += p.memsize();
  
  return size;
}//size_t memsize() const;


size_t MassAttenuationTool::memsize() const
{
  size_t size = sizeof(*this);
  size += m_dataPath.capacity();
  
  for( const auto &e : m_atten )
  {
    const ElementAttenuation *ptr = e.load();
    size += ptr ? ptr->memsize() : size_t(0);
  }
  
  return size;
}//size_t ElementAttenuation::memsize() const


void ElementAttenuation::saveTxt( std::string path )
{
  char filename[12];
  snprintf( filename, sizeof(filename), "%i.xs.txt", m_atomicNumber );
  path = append_path( path, filename );
  
  FILE *pFile = fopen( path.c_str(), "w" );
  
  if( !pFile )
    throw runtime_error( "Coulnt open file " + path );
  
  fprintf( pFile, "%s %1.8E %i\n", m_symbol.c_str(), m_atomicMass, m_atomicNumber );
  
  char buffer[128];
  
  for( const auto &proccess : m_proccesses )
  {
    const vector<float> &energies = proccess.m_logEnergies;
    for( size_t pos = 0; pos < energies.size(); ++pos )
    {
      exact_short_float( energies[pos], buffer, sizeof(buffer) );
      fprintf( pFile, (pos ? " %s" : "%s"), buffer );
    }
    fprintf( pFile, "\n" );
    
    const vector<float> &attcoefs = proccess.m_logAttenuationCoeffs;
    for( size_t pos = 0; pos < attcoefs.size(); ++pos )
    {
      exact_short_float( attcoefs[pos], buffer, sizeof(buffer) );
      fprintf( pFile, (pos ? " %s" : "%s"), buffer );
    }
    fprintf( pFile, "\n" );
  }//for( loop over proccesses )
  
  fclose( pFile );
}//void saveTxt( string datapath )



void ElementAttenuation::loadTxt( std::string datapath, const int atomicNumber )
{
  char filename[12];
  snprintf( filename, sizeof(filename), "%i.xs.txt", atomicNumber );
  datapath = append_path( datapath, filename );
  
  ifstream file( datapath.c_str(), ios_base::binary|ios_base::in );
  
  if( !file.is_open() || !file.good() )
    throw runtime_error( "Coulnt open file " + datapath );
  
  if( !(file >> m_symbol >> m_atomicMass >> m_atomicNumber) )
    throw runtime_error( "Error reading first line of: " + datapath );
  
  int nextchar = file.peek();
  while( (nextchar=='\n') && nextchar!=EOF )
  {
    file.get();
    nextchar = file.peek();
  }
  
  for( int i = 0; i < static_cast<int>(MassAttenuation::GammaEmProcces::NumGammaEmProcces); i += 1 )
  {
    string line;
    if( !std::getline( file, line, '\n' ) )
    {
      char msg[256];
      snprintf( msg, sizeof(msg),
                "Fail to read energy line for %s, proccess %i",
                 m_symbol.c_str(), int(i) );
      throw runtime_error( msg );
    }
    
    vector<float> &energies = m_proccesses[i].m_logEnergies;
    vector<float> &attcoefs = m_proccesses[i].m_logAttenuationCoeffs;
    
    bool success = split_space_delim_flts( line, energies );
    
    if( !success || energies.empty() )
    {
      char msg[256];
      snprintf( msg, sizeof(msg),
                "Failed to decode energy line for %s, proccess %i",
                m_symbol.c_str(), int(i) );
      throw runtime_error( msg );
    }
    
    if( !std::getline( file, line, '\n' ) )
    {
      char msg[256];
      snprintf( msg, sizeof(msg),
                "Fail to read attenuation line for %s, proccess %i",
                m_symbol.c_str(), int(i) );
      throw runtime_error( msg );
    }

    success = split_space_delim_flts( line, attcoefs );
    if( !success || attcoefs.empty() )
    {
      char msg[256];
      snprintf( msg, sizeof(msg), "Failed to decode attenuation line for %s, proccess %i", m_symbol.c_str(), int(i) );
      throw runtime_error( msg );
    }
    
    if( attcoefs.size() != energies.size() )
      throw runtime_error( "Attenuation coefficient size != energy size" );
  }//for( get proccess )
}//void loadTxt( std::string datapath, const int atomicNumber )



MassAttenuationTool::MassAttenuationTool()
  : m_dataPath( "data/em_xs_data" )
{
  static_assert( (sizeof(m_atten)/sizeof(m_atten[0])) == MassAttenuation::sm_max_xs_atomic_number, "" );
  
  for( auto &p : m_atten )
    p = nullptr;
}

void MassAttenuationTool::set_data_dir( const std::string &dir )
{
  for( auto &p : m_atten )
    if( p.load() )
      throw runtime_error( "MassAttenuationTool::set_data_dir(): you can not call this function after loading any cross-sections!" );
      
  m_dataPath = dir;
}

MassAttenuationTool::~MassAttenuationTool()
{
  for( std::atomic<const ElementAttenuation *> &el : m_atten )
  {
    const ElementAttenuation *expected = nullptr;
    const bool changed = el.compare_exchange_strong( expected, nullptr );
    if( !changed )
      delete expected;
  }
}//~MassAttenuationTool()


const ElementAttenuation *MassAttenuationTool::attenuationData(
                                                      const int atomic_number )
{
  if( atomic_number < 1 || atomic_number > 98 )
    throw runtime_error( "Invalid atomic number" );
    
  const ElementAttenuation *origptr = m_atten[atomic_number-1].load();
  if( origptr )
    return origptr;
  
  ElementAttenuation *thisData = nullptr;
  
  try
  {
    thisData = new ElementAttenuation();
    thisData->loadTxt( m_dataPath, atomic_number );

    const bool changed = m_atten[atomic_number-1].compare_exchange_strong( origptr, thisData );
      
    if( changed )
    {
      //thisData was stored into m_atten[atomic_number-1] since because it was
      //  still null
      origptr = thisData;
    }else
    {
      //Another thread updated m_atten[atomic_number-1], so delete the copy we
      //  made, and return the copy the other thread made (origptr)
      delete thisData;
    }
  }catch( std::exception &e )
  {
    if( thisData )
      delete thisData;
    thisData = nullptr;
    cerr << "Caught: " << e.what() << endl;
    throw MassAttenuation::ErrorLoadingDataException( string("Error Loading Data: ") + string(e.what()) );
  }//try / catch
    
  return origptr;
}//attenuationData(...)


float MassAttenuationTool::massAttenuationCoeficient( const int atomic_num,
                                                      const float energy )
{
#if( USE_SNL_GAMMA_ATTENUATION_VALUES )
  const float units = static_cast<float>( PhysicalUnits::cm2 / PhysicalUnits::gram );
  float s_mu, p_mu, pair_mu;
  return units * AttCoef( energy, atomic_num, s_mu, p_mu, pair_mu );
#else
  const ElementAttenuation * const data = attenuationData( atomic_num );

  float comptXs = 0.0, photoXs = 0.0, convXs = 0.0;
  try
  {
    comptXs = calcMassAttenuationCoeficient( energy, MassAttenuation::GammaEmProcces::ComptonScatter, data );
  }catch(...){}

  try
  {
    photoXs = calcMassAttenuationCoeficient( energy, MassAttenuation::GammaEmProcces::PhotoElectric, data );
  }catch(...){}

  try
  {
    if( energy > 1024.0*PhysicalUnits::keV )
    {
      convXs = calcMassAttenuationCoeficient( energy, MassAttenuation::GammaEmProcces::PairProduction, data );
    }
  }catch(...){}
  
  return comptXs + photoXs + convXs;
#endif
}//float massAttenuationCoeficient(...)


float MassAttenuationTool::massAttenuationCoeficientFracAN( const float atomic_number, const float energy )
{
  const int floor_an = static_cast<int>( std::floor( atomic_number ) );
  
  if( atomic_number >= MassAttenuation::sm_max_xs_atomic_number || (atomic_number == floor_an) )
    return massAttenuationCoeficient( static_cast<int>(atomic_number), energy );
  
  const int next_an = floor_an + 1;
  
  const float muf = massAttenuationCoeficient( floor_an, energy );
  const float mup1 = massAttenuationCoeficient( next_an, energy );
  
  const float anfrac = min( 1.0f, max( 0.0f, atomic_number-floor_an ) );  //the min/max can probably be removed, but leaving in JIC
  const float mu = (1.0f - anfrac)*muf + anfrac*mup1;

  return mu;
}//float MassAttenuationTool::massAttenuationCoeficientFracAN( const float atomic_number, const float energy )



float MassAttenuationTool::massAttenuationCoeficient( const int atomic_number,
                                                      const float energy,
                                                      MassAttenuation::GammaEmProcces process )
{
#if( USE_SNL_GAMMA_ATTENUATION_VALUES )
  float scatter_mu, photoelectric_mu, pair_prod_mu;
  const float units = static_cast<float>(PhysicalUnits::cm2 / PhysicalUnits::gram);
  
  AttCoef( energy, atomic_number, scatter_mu, photoelectric_mu, pair_prod_mu );
  
  switch( process )
  {
    case kComptonScatter:  return units*scatter_mu;
    case kPairProduction:  return units*pair_prod_mu;
    case kPhotoElectric:   return units*photoelectric_mu;
    case kRayleighScatter: return 0.0f;
    
    case kNumGammaEmProcces:
    default:
      throw runtime_error( "Invalis EM Proccess" );
      break;
  }
#endif
  
  if( static_cast<int>(process) >= static_cast<int>(MassAttenuation::GammaEmProcces::NumGammaEmProcces) )
    throw runtime_error( "Invalis EM Proccess" );

  const ElementAttenuation *data = attenuationData( atomic_number );

  if( data->m_proccesses[static_cast<int>(process)].m_logEnergies.empty() )
    throw runtime_error( "Not-loaded data" );

  return logLogInterpolate( energy,
                         data->m_proccesses[static_cast<int>(process)].m_logEnergies,
                         data->m_proccesses[static_cast<int>(process)].m_logAttenuationCoeffs );
}//float massAttenuationCoeficient(...)

}
