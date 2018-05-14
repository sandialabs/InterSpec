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

#include <mutex>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <condition_variable>

#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "sandia_decay/SandiaDecay.h"
#include "SpecUtils/UtilityFunctions.h"

using namespace std;


const Material MaterialDB::sm_voidMaterial( "void", 0.0 );

Material::Material()
{
}

Material::Material( const Material &rhs )
  : name( rhs.name ),
    description( rhs.description ),
    density( rhs.density ),
    source( rhs.source ),
    nuclides( rhs.nuclides ),
    elements( rhs.elements )
{
}

Material::Material( const std::string &_name, const float _density )
  : name( _name ), density( _density )
{
}

Material::Material( const std::string &_name,
          const std::string &_description, const float _density )
  : name( _name ), description( _description ), density( _density )
{
}

Material::~Material()
{

}


string Material::chemicalFormula() const
{
  stringstream formula;
  const double cm3PerG = PhysicalUnits::cm3 / PhysicalUnits::g;

  for( const NuclideFractionPair &n : nuclides )
  {
    if( n.second > 0.0f )
      formula << "'" <<  n.first->symbol << "'" << (n.second * density * cm3PerG);
  }

  for( const ElementFractionPair &e : elements )
  {
    if( e.second > 0.0f )
      formula << e.first->symbol << (e.second * density * cm3PerG);
  }

  return formula.str();
}//string chemicalFormula() const



void Material::parseChemicalFormula( string text,
                                    const SandiaDecay::SandiaDecayDataBase *db )
{
  density = 0.0;
  source = kUser;

  if( text.empty() )
  {
    name = description = "";
    throw runtime_error( "Material::parseChemicalFormula(...): I cant deal "
                         "with no input" );
  }//if( text.empty() )

  name = description = text;

  string::size_type densitypos = text.find( "d=" );

  string densitystr;
  if( densitypos != string::npos )
  {
    densitystr = text.substr( densitypos + 2 );
    text = text.substr( 0, densitypos );
    cerr << "densitystr=\"" << densitystr << "\", text=\"" << text << "\"" << endl;
  }//if( densitypos != string::npos )

  float totalfraction = 0.0;
  UtilityFunctions::erase_any_character( text, " \t\n_-,;" );

  const string numbers = "0123456789.";
  string::size_type pos = 0;

  //Could probably use a regex, but the bellow seems to work okay
  try
  {
    while( pos != string::npos && pos != text.size() )
    {
      float fraction = 0.0;
      string::size_type numend;

      if( text[pos] == '\'' || text[pos] == '\"' )
      { //we are reading a isotope here, e.g. 'U-235'
        string::size_type quotpos = text.find( text[pos], pos + 1 );
        if( quotpos == string::npos )
          throw runtime_error( "Didnt find matching \' in \"" + text + "\"");

        const string isostr = text.substr( pos+1, quotpos - pos - 1 );
        cerr << "isostr=\"" << isostr << "\"" << endl;
        string::size_type numstart = text.find_first_of( numbers, quotpos );
        if( pos == numstart || numstart==string::npos )
          throw runtime_error( "Error finding end of element string after pos "
                               + std::to_string(pos) );
        numend = text.find_first_not_of( numbers, numstart );
        if( numend == string::npos )
          numend = text.size();

        if( numend == numstart )
          throw runtime_error( "Error finding end of density string after pos "
                               + std::to_string(numstart) );

        const string numberstr = text.substr( numstart, numend - numstart );

        const SandiaDecay::Nuclide *nuc = db->nuclide( isostr );
        if( !nuc )
          throw runtime_error( isostr + " is not a nuclide" );

        if( !(stringstream(numberstr) >> fraction) )
          throw runtime_error( numberstr + " is not a valid decimal number" );

        const pair<const SandiaDecay::Nuclide *, float> thisone( nuc, fraction );
        nuclides.push_back( thisone );
      }else
      { //We can assume we are reading an element here
        string::size_type numstart = text.find_first_of( numbers, pos );
        if( pos == numstart || numstart==string::npos )
          throw runtime_error( "Error finding end of element string after pos "
                               + std::to_string(pos) );

        const string elementstr = text.substr( pos, numstart - pos );

        numend = text.find_first_not_of( numbers, numstart );
        if( numend == string::npos )
          numend = text.size();

        if( numend == numstart )
          throw runtime_error( "Error finding end of density string after pos "
                               + std::to_string(numstart) );

        const string numberstr = text.substr( numstart, numend - numstart );

        const SandiaDecay::Element *el = db->element( elementstr );
        if( !el )
          throw runtime_error( elementstr + " is not an element" );

        if( !(stringstream(numberstr) >> fraction) )
          throw runtime_error( numberstr + " is not a valid decimal number" );

        const pair<const SandiaDecay::Element *, float> thisone( el, fraction );
        elements.push_back( thisone );
      }//if( text[pos] == '\'' )

      totalfraction += fraction;
      pos = numend;
    }//while( pos != string::npos )

//    if( fabs(totalfraction-1.0) > 0.02 )
//      passMessage( "Total mass fraction does not add up to 1.0, rescaling all"
//                   " mass fractions to make this true",
//                   "GammaXsGui::parseMaterial", WarningWidget::WarningMsgInfo );

    density = totalfraction * static_cast<float>(PhysicalUnits::g / PhysicalUnits::cm3);

    for( Material::ElementFractionPair &nf : elements )
      nf.second /= totalfraction;

    for( Material::NuclideFractionPair &nf : nuclides )
      nf.second /= totalfraction;

    if( densitystr.size() )
    {
      try
      {
        density = static_cast<float>( std::stod( densitystr ) );
        density *= static_cast<float>(PhysicalUnits::g / PhysicalUnits::cm3);
      }catch(...)
      {
        passMessage( "Couldnt cast '" + densitystr + "' to a density, sorry,"
                     " not using is, becareful!",
                     "Material::parseChemicalFormula", WarningWidget::WarningMsgInfo );
      }
    }//if( densitystr.size() )
  }catch( std::exception &e )
  {
    density = 0.0;
    name = description = "";
    throw e;
  }//try / catch

}//void parseChemicalFormula(  )


void Material::writeGadrasStyleMaterialFile( ostream &file ) const
{
  const char *ending = "\r\n";

  size_t nstable = 0, nradioactive = 0;
  nstable += elements.size();

  for( const Material::NuclideFractionPair &nfp : nuclides )
    if( IsInf(nfp.first->halfLife) || IsNan(nfp.first->halfLife) )
      nstable++;
    else
      nradioactive++;

  file << "!" << ending
       << name << ending;

  if( description.empty() )
    file << name << ending;
  else
    file << description << ending;

  file << nstable << " " << nradioactive << " "
       << density * PhysicalUnits::cm3 / PhysicalUnits::g << ending;

  //Write elements out
  for( const Material::ElementFractionPair &elp : elements )
    file << elp.first->symbol << ending << 100.0*elp.second << ending;

  //Write stable nuclides out
  for( const Material::NuclideFractionPair &nfp : nuclides )
    if( IsInf(nfp.first->halfLife) || IsNan(nfp.first->halfLife) )
      file << nfp.first->symbol << ending << 100.0*nfp.second << ending;

  //Write "source" nuclides out
  for( const Material::NuclideFractionPair &nfp : nuclides )
    if( !(IsInf(nfp.first->halfLife) || IsNan(nfp.first->halfLife)) )
      file << nfp.first->symbol << ending << 100.0*nfp.second << ending;

}//void writeGadrasStyleMaterialFile( std::ostream file )


float Material::massWeightedAtomicNumber() const
{
  float totalFraction = 0.0f, massAnTotalFrac = 0.0f;
  
  for( const Material::ElementFractionPair &nf : elements )
  {
    totalFraction += nf.second;
    massAnTotalFrac += nf.second * nf.first->atomicNumber;
  }
  
  for( const Material::NuclideFractionPair &nf : nuclides )
  {
    totalFraction += nf.second;
    massAnTotalFrac += nf.second * nf.first->atomicNumber;
  }
  
  if( totalFraction == 0.0f )
    return 0.0f;
  
  return (massAnTotalFrac / totalFraction);
}//float massWeightedAtomicNumber() const


MaterialDB::MaterialDB()
{
  m_initFailure = false;
  m_materials.push_back( &sm_voidMaterial );
  refreshMaterialNames();
}



MaterialDB::~MaterialDB()
{
  for( const Material *m : m_materials )
    if( m != &sm_voidMaterial )
      delete m;
  m_materials.clear();
}


const std::vector<std::string> &MaterialDB::names() const
{
  if( !wait_material_ready() )
    throw std::runtime_error( "Material database is not initialized" );
  
  return m_materialNames;
}//std::vector<std::string> MaterialDB::names() const


bool MaterialDB::wait_material_ready() const
{
  std::unique_lock<std::mutex> lock( m_mutex );
  
  if( m_materials.size() > 1 )
    return true;
  
  if( m_initFailure )
    return false;
  
  if( m_materials.size() < 2 )
    m_condition.wait_for( lock, std::chrono::milliseconds(1000), [this](){return m_materials.size()>1;} );
  
  m_initFailure = (m_materials.size() < 2);
  
  return !m_initFailure;
}//bool wait_material_ready() const


const Material *MaterialDB::material( const std::string &name ) const
{
  if( name.empty() )
    throw std::runtime_error( "Empty material name" );

  if( !wait_material_ready() )
    throw std::runtime_error( "Material database is not initialized" );
  
  for( const Material *m : m_materials )
  {
    //we had added some comments in (...), so lets gets allow the user
    //  to either enter the whole thing eg "Fe (iron)" or just "Fe", or "iron"
    const string thisname = m->name;
    string thissymbol, thisdescrip;
    const string::size_type pos = thisname.find( " (" );
    if( pos != string::npos )
    {
      thissymbol = thisname.substr( 0, pos );
      const string::size_type closing = thisname.find( ")", pos+2 );
      if( closing != string::npos )
        thisdescrip = thisname.substr( pos+2, closing - pos - 2 );
    }//if( pos != string::npos )

    if( UtilityFunctions::iequals( thisname, name )
        || UtilityFunctions::iequals( thissymbol, name )
        || UtilityFunctions::iequals( thisdescrip, name ) )
      return m;
  }//for( const Material *m : m_materials )

  throw std::runtime_error( "Couldnt find material '" + name + "'" );
  return NULL;
}//const Material *material( &std::string &name ) const


const Material *MaterialDB::material( const size_t index ) const
{
  if( !wait_material_ready() )
    throw std::runtime_error( "Material database is not initialized" );
  
  if( index >= m_materials.size() )
    throw std::runtime_error( "MaterialDB::material( size_t index ): index out of bounds" );
  return m_materials[index];
}//const Material *material( &std::string &name ) const


void MaterialDB::refreshMaterialNames()
{
  //Try to get the lock for 10 milli-seconds
  //  Another (b=maybye better) approach would be to use boost::timed_mutex
  int ntry = 0;
  std::unique_lock<std::mutex> lock( m_mutex, std::defer_lock );
  while( (ntry < 1000) && !lock.try_lock() )
  {
    ++ntry;
    std::this_thread::sleep_for( std::chrono::microseconds(10) );
  }
  
  if( ntry >= 10 )
    throw runtime_error( "MaterialDB::refreshMaterialNames(): unable to get mutex lock" );
  
  m_materialNames.clear();
  m_materialNames.reserve( m_materials.size() );

  for( const Material *m : m_materials )
    m_materialNames.push_back( m->name );
}//void MaterialDB::refreshMaterialNames()


void MaterialDB::writeGadrasStyleMaterialFile( ostream &file ) const
{
  const Material *air = material( "Air" );
  if( air )
    air->writeGadrasStyleMaterialFile( file );

  std::unique_lock<std::mutex> lock( m_mutex );
  
  for( const Material *mat : m_materials )
    if( !UtilityFunctions::iequals(mat->name, "Air") )
      mat->writeGadrasStyleMaterialFile( file );
}//void writeGadrasStyleMaterialFile( std::ostream file )


void MaterialDB::parseGadrasMaterialFile( const std::string &file,
                                          const SandiaDecay::SandiaDecayDataBase *db,
                                          const bool swapNameDescription )
{
  //Gadras Material names are formatted like
  /*
    !
    Name
    Description
    [INT num stable elements/isotopes] [INT num radioactive elements/isotopes] [float density (g/cm3)]
    stable Element/Isotope 0 Symbol
    stable Element/Isotope 0 Mass percentage
    stable Element/Isotope 1 Symbol
    stable Element/Isotope 1 Mass percentage
    radioactive Element/Isotope 0 Symbol (ex LI6, Th, U, Co60)
    radioactive Element/Isotope 0 Mass percentage
    radioactive Element/Isotope 1 Symbol (ex LI6, Th, U, Co60)
    radioactive Element/Isotope 1 Mass percentage
    ...
  */
  
  int lineNumber = 0;
  Material *material = NULL;
  std::vector<const Material *> materials;

  try
  {
    ifstream is( file.c_str(), ios_base::binary|ios_base::in );
    if( !is.is_open() )
    {
      cerr << "Couldnt open file " << file << endl;
      throw runtime_error( "Couldnt open file " + file );
    }


    string line;
    bool has_hit_air = false;

    while( UtilityFunctions::safe_get_line( is, line ) )
    {
      ++lineNumber;

      UtilityFunctions::trim( line );

      if( !has_hit_air )
        has_hit_air = UtilityFunctions::starts_with(line,"Air");

      if( !has_hit_air )
        continue;

      if( UtilityFunctions::starts_with(line,"!") )
      {
        //We dont want to parse the "Notes" section at the bottom of the file
        if( line.find("Notes") != string::npos )
          break;
        continue;
      }//if( UtilityFunctions::starts_with(line,"!") )

      if( line.empty() || UtilityFunctions::starts_with(line,"#") )
        continue;

      //presumably, if we are here, we are at the _start_ of a new material
      material = new Material();
      material->source = Material::kGadras;
      material->name = line;
      bool gotit = !!UtilityFunctions::safe_get_line( is, line );  //should really test this
      if( !gotit )
        throw std::runtime_error( "Unexpected formating A" );
      ++lineNumber;

      material->description = line;
      const string::size_type exclamation_pos = line.find("!");
      if( exclamation_pos != string::npos )
        material->description = line.substr( 0, exclamation_pos );

      UtilityFunctions::safe_get_line( is, line );
      ++lineNumber;

      const string::size_type comment_pos = line.find( '!' );
      if( comment_pos != string::npos )
        line = line.substr( 0, comment_pos );

      UtilityFunctions::trim( line );

      vector<string> components_something_density;
      UtilityFunctions::split( components_something_density, line, " \t" );
      if( components_something_density.size() != 3 )
        throw runtime_error( "Expected 3 fields in line '" + line + "'" );

      if( !(stringstream(components_something_density[2]) >> (material->density)) )
        throw runtime_error( "What I expected to be density - '"
                             + components_something_density[2]
                             + " was not valid number" );

      material->density *= static_cast<float>(PhysicalUnits::g/PhysicalUnits::cm3);

      int nexpectedelement, nexpectednuclide;
      if( !(stringstream(components_something_density[0]) >> nexpectedelement) )
        throw runtime_error( "What I expected to be number of elements - '"
                             + components_something_density[0]
                             + " was not valid integer" );
      if( !(stringstream(components_something_density[1]) >> nexpectednuclide) )
        throw runtime_error( "What I expected to be number of nuclides - '"
                             + components_something_density[1]
                             + " was not valid integer" );


      int nelement = 0, nisotope = 0;
      while( UtilityFunctions::safe_get_line( is, line ) )
      {
        ++lineNumber;
        UtilityFunctions::trim( line );
        if( line.empty() || UtilityFunctions::starts_with(line,"#") )
          continue;

        if( UtilityFunctions::starts_with(line, "!") )
          break;

        const string symbol = line;
        gotit = !!UtilityFunctions::safe_get_line( is, line );  //should really test this
        if( !gotit )
          throw std::runtime_error( "Unexpected formating B" );
        ++lineNumber;

        const string mass_percentage = line;

        double percentage;
        if( !(stringstream(mass_percentage) >> percentage) )
          throw std::runtime_error( "Expected a numberical line following '"
                                    + symbol + "', but instead got '"
                                    + mass_percentage );
        const float mass_fraction = static_cast<float>(percentage / 100.0);

        const SandiaDecay::Element *element = NULL;
        const SandiaDecay::Nuclide *nuc = db->nuclide( symbol );
        if( !nuc )
          element = db->element( symbol );

        if( !nuc && !element )
          throw std::runtime_error( "'" + symbol
                                    + "' is not an element or isotope. "
                                    " On material " + material->name );

        if( nuc )
        {
          //We should check if isotope is 100% natural abundance of that element,
          //  and if so, just put the element in instead of the isotope
          ++nisotope;
          material->nuclides.push_back( Material::NuclideFractionPair(nuc,mass_fraction) );
        }

        if( element )
        {
          ++nelement;
          material->elements.push_back( Material::ElementFractionPair(element,mass_fraction) );
        }
      }//while( UtilityFunctions::safe_get_line( is, line ) )

      if( (nexpectedelement + nexpectednuclide) != (nelement + nisotope) )
        throw runtime_error( "Expected and found number of elements+nuclides"
                             " didnt match for material " + material->name );

//      if( nexpectedelement != nelement )
//        cerr << "Expected and found number of elements"
//                " didnt match for material " + material->name << endl;
//        throw runtime_error( "Expected and found number of elements"
//                             " didnt match for material " + material->name );

//      if( nexpectednuclide != nisotope )
//        cerr << "Expected and found number of isotopes"
//                " didnt match for material " + material->name << endl;
//        throw runtime_error( "Expected and found number of isotopes"
//                             " didnt match for material " + material->name );


      //I actually like the description in the GADRAS material library better
      //  as the names
      if( swapNameDescription )
        swap( material->name, material->description );


      if( material->elements.size() == 1 && material->nuclides.empty() )
      {
        const SandiaDecay::Element *el = db->element( material->name );
        if( el )
          material->name = el->symbol + " (" + el->name + ")";
      }//if( material->elements.size() == 1 && material->nuclides.empty() )

      //lets make sure we dont already have this material actually
      std::vector<const Material *>::iterator pos;
      pos = lower_bound( m_materials.begin(), m_materials.end(),
                         material, &MaterialDB::less_than_by_name );

      if( (pos != m_materials.end() &&
          (MaterialDB::less_than_by_name(material,*pos)
             == MaterialDB::less_than_by_name(*pos,material)))
          || UtilityFunctions::iequals(material->name,"void") )
      {
/*
        cerr << "Found duplicate material: " << material->name
             << " with density " << material->density << " (" << (*pos)->density
             << ") and " << material->elements.size() << " (" << (*pos)->elements.size()
             << ") elements and " << material->nuclides.size() << " (" << (*pos)->nuclides.size()
             << ")" << endl;
*/
        delete material;
        material = NULL;
      }//if( we already have this material )

      if( material )
        materials.push_back( material );
      material = NULL;
    }//while( safe_get_line( is, line ) )
  }catch( const std::exception &e )
  {
    if( material )
      delete material;

    for( const Material *m : materials )
      if( m )
        delete m;

    stringstream msg;
    msg << "MaterialDB::parseG4MaterialFile(...): Error parsing file  at or "
        << "near line " << lineNumber << ". error=" << e.what();

    cerr << "\n\n" << SRC_LOCATION << "\n\t" << msg.str() << endl << endl;
    
    m_condition.notify_all();
    
    throw std::runtime_error( msg.str() );
  }//try /catch

  {
    std::unique_lock<std::mutex> lock( m_mutex );
    m_materials.insert( m_materials.end(), materials.begin(), materials.end() );
    sort( m_materials.begin(), m_materials.end(), &MaterialDB::less_than_by_name );
  }

  refreshMaterialNames();
  
  m_condition.notify_all();
}//void parseGadrasMaterialFile(...)



void MaterialDB::parseG4MaterialFile( const std::string &file,
                            const SandiaDecay::SandiaDecayDataBase *db )
{
  int lineNumber = 0;
  Material *material = NULL;
  std::vector<const Material *> materials;

  try
  {
    ifstream is( file.c_str(), ios_base::binary|ios_base::in );
    if( !is.is_open() )
    {
      cerr << "Couldnt open file " << file << endl;
      throw runtime_error( "Couldnt open file " + file );
    }
/*
    is.seekg( 0, ios::end );
    ifstream::pos_type length = is.tellg();
    is.seekg( 0, ios::beg );
    std::unique_ptr<char> buffer( new char [length+1] );
    is.read( buffer.get(),length );
    buffer.get()[length] = '\0';
    is.close();
*/

    bool is_compound = true;
    string line;
    while( UtilityFunctions::safe_get_line( is, line ) )
//    while( getline( is, line ) )
    {
      ++lineNumber;

      if( line.empty() || UtilityFunctions::starts_with(line,"#") )
        continue;

      vector<string> fields;
      const bool is_material = (isdigit(line[0]) != 0);
      UtilityFunctions::trim( line );
      UtilityFunctions::split( fields, line, " \t" );

      if( fields.size() < 2 )
        continue;

      const bool is_component = !is_material && isdigit(fields[0][0]);  //we're garunteed fields[0].size()>0

      if( is_component )
        throw std::runtime_error( "Found a component line unexpectedly" );

      bool is_elements_header = false, is_material_header = false;

      if( UtilityFunctions::iequals(fields[1],"Name") )
      {
        is_elements_header = UtilityFunctions::iequals( fields[0], "Z" );
        is_material_header = UtilityFunctions::iequals( fields[0], "Ncomp" );
      }//if( UtilityFunctions::iequals(fields[0],"Name") )

      if( is_elements_header )
        is_compound = false;

      if( is_material_header )
        is_compound = true;

      if( is_elements_header || is_material_header )
        continue;

      if( !is_material  )
      {
        cerr << "Line: " << line << endl
             << "\tIs not a material like expected" << endl;
        continue;
      }//if( !is_material && !is_component )

      material = new Material();
      material->source = Material::kNist;

      if( is_compound )
      {
        assert( fields.size() >= 3 );
        if( fields.size() < 3 )
          throw std::runtime_error( "Compound line doesnt have "
                                    "at least 3 fields" );

        const int ncomp = std::stoi( fields[0] );
        material->name = fields[1];
        UtilityFunctions::to_lower( material->name );
        material->density = static_cast<float>( std::stod(fields[2]) * PhysicalUnits::g/PhysicalUnits::cm3 );

        for( int compn = 0; compn < ncomp; ++compn )
        {
//          if( !getline( is, line ) )
          if( !UtilityFunctions::safe_get_line( is, line ) )
            throw std::runtime_error( "Couldnt read expected number of "
                                      "components of compound "
                                      + material->name );

          UtilityFunctions::trim( line );
          stringstream linestrm( line );
          int z;
          float frac;
          linestrm >> z >> frac;

          if( !linestrm )
            throw std::runtime_error( "Error reading folowing line as a "
                                      "compound component '" + line + "'" );
          const SandiaDecay::Element *element = db->element( z );
          if( !element )
            throw std::runtime_error( "There is no element with Z="
                                      + std::to_string(z) );
          const Material::ElementFractionPair component( element, frac );
          material->elements.push_back( component );
        }//for( int compn = 0; compn < ncomp; ++compn )
      }else
      {
        assert( fields.size() >= 3 );
        if( fields.size() < 3 )
          throw std::runtime_error( "Element line doesnt have >= 3 fields" );

        const int z = std::stoi( fields[0] );
        material->name = fields[1];
        material->density = static_cast<float>( std::stod(fields[2]) * PhysicalUnits::g/PhysicalUnits::cm3 );

        const SandiaDecay::Element *element = db->element( z );
        if( !element )
          throw std::runtime_error( "There is no element with Z="
                                    + std::to_string(z) );
        const Material::ElementFractionPair component( element, 1.0f );
        material->elements.push_back( component );
      }//if( is_compound ) / else

      if( UtilityFunctions::istarts_with( material->name, "G4_" ) )
        material->name = material->name.substr( 3 );


      if( material->elements.size() == 1 && material->nuclides.empty() )
      {
        const SandiaDecay::Element *el = db->element( material->name );
        if( el )
          material->name = el->symbol + " (" + el->name + ")";
      }//if( material->elements.size() == 1 && material->nuclides.empty() )

      //lets make sure we dont already have this material actually
      std::vector<const Material *>::iterator pos;
      pos = lower_bound( m_materials.begin(), m_materials.end(),
                         material, &MaterialDB::less_than_by_name );

      if( pos != m_materials.end() &&
          (MaterialDB::less_than_by_name(material,*pos)
             == MaterialDB::less_than_by_name(*pos,material)) )
      {
        delete material;
        material = NULL;
      }//if( we already have this material )

      if( material )
        materials.push_back( material );
      material = NULL;
    }//while( safe_get_line( is, line ) )
  }catch( const std::exception &e )
  {
    if( material )
      delete material;

    for( const Material *m : materials )
      if( m )
        delete m;

    stringstream msg;
    msg << "MaterialDB::parseG4MaterialFile(...): Error parsing file  at or "
        << "near line " << lineNumber << ". error=" << e.what();

    cerr << "\n\n" << SRC_LOCATION << "\n\t" << msg.str() << endl << endl;

    m_condition.notify_all();
    
    throw std::runtime_error( msg.str() );
  }//try /catch

  {
    std::unique_lock<std::mutex> lock( m_mutex );
    m_materials.insert( m_materials.end(), materials.begin(), materials.end() );
    sort( m_materials.begin(), m_materials.end(), &MaterialDB::less_than_by_name );
  }

  refreshMaterialNames();
  
  m_condition.notify_all();
}//void parseG4MaterialFile(...)



const Material *MaterialDB::parseChemicalFormula( string text,
                                   const SandiaDecay::SandiaDecayDataBase *db )
{
  if( text.empty() )
    throw runtime_error( "GammaXsGui::parseChemicalFormula(...): I cant deal "
                         "with no input" );

  //first look to see if its in the database
  try
  {
    const Material *material = MaterialDB::material( text );
    assert( material );
    if( material )
      throw runtime_error( "Material '" + text + "' already exists in database" );
  }catch(...){}

  //If we're here, we have to parse the chemical formula, which might look like
  Material *material = new Material();

  try
  {
    material->parseChemicalFormula( text, db );
    vector<const Material *>::iterator iter = lower_bound( m_materials.begin(), m_materials.end(), material, &MaterialDB::less_than_by_name );
    
    {
      std::unique_lock<std::mutex> lock( m_mutex );
      m_materials.insert( iter, material );
    }
    
    refreshMaterialNames();

    cerr << "\n\nJust parsed and added: " << material->chemicalFormula() << endl << endl;
  }catch( std::exception &e )
  {
    delete material;
    throw e;
  }//try / catch

  return material;
}//const Material *parseChemicalFormula(  )




bool MaterialDB::less_than_by_name( const Material *lhs, const Material *rhs )
{
  string lhsname = UtilityFunctions::to_lower_copy(lhs->name);
  string rhsname = UtilityFunctions::to_lower_copy(rhs->name);
  UtilityFunctions::erase_any_character( lhsname, "_- \t " );
  UtilityFunctions::erase_any_character( rhsname, "_- \t " );
  return (lhsname < rhsname);
}

