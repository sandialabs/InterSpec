#ifndef MaterialDB_h
#define MaterialDB_h
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

#include <mutex>
#include <string>
#include <vector>
#include <utility>
#include <istream>
#include <condition_variable>

class MaterialDB;

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
  class SandiaDecayDataBase;
}//namespace SandiaDecay


struct Material
{
  typedef std::pair<const SandiaDecay::Element *,float> ElementFractionPair;
  typedef std::pair<const SandiaDecay::Nuclide *,float> NuclideFractionPair;

  enum MaterialDefintionsSrc
  {
    kGadras,
    kNist,
    kUser
  };//enum MaterialDefSrc

  //Material: this constructor is intended to be used for user editable
  //  materials which are not tracked by the MaterialDB class
  Material( const Material &rhs );

  ~Material();

  //chemicalFormula(): returns a string in the form of an elements followed
  //  by there density (in g/cm3).  Isotopic components have their names quoted
  //  in single quotes.  (fictional) Examples include: "C0.3O0.9Fe3.2",
  //  "'U235'0.1'U233'0.2'U238'1.1"
  std::string chemicalFormula() const;

  //parseChemicalFormula(...): parses chemical formulas such as
  //  "C0.5H0.2Ni0.3 d=2.2", where d=2.2 is the density in g/cm3 of the material
  //  or equivalently you could write "C1.1H0.44Ni0.66".
  //  If the "d=..." option is specified, than the elemental coefecients are
  //  treated as relaitive mass fractions (so sum of them is normalized to 1).
  //  Before parsing, name, description, density and compnents of the material
  //  are reset.  Isotopic symbols should be quoted (single or double), such
  //  as in "'U235'0.1'U233'0.2'U238'1.1Fe1.2"
  //  If the formula cant be parsed, than an exception is thrown.
  void parseChemicalFormula( std::string formula,
                             const SandiaDecay::SandiaDecayDataBase *db );


  //writeGadrasStyleMaterialFile(...):
  //  Writes the materials back out to a GADRAS compatible format, using Windows
  //  line endings.
  void writeGadrasStyleMaterialFile( std::ostream &file ) const;

  
  float massWeightedAtomicNumber() const;

  std::string name;
  std::string description;
  float density;  //in units of PhysicalUnits
  MaterialDefintionsSrc source;

  //nuclides and their fractions by weight
  std::vector< NuclideFractionPair >  nuclides;

  //elements and their fractions by weight
  std::vector< ElementFractionPair >  elements;


private:
  Material();
  Material( const std::string &_name, const float _density );
  Material( const std::string &_name,
            const std::string &_description, const float _density );

  friend class MaterialDB;
};//struct Material


class MaterialDB
{
/*
   MaterialDB: keeps an entry of materials in memorry.
   Text based files in either GEANT4 or GADRAS style formating may be used to
   populate database.  
   This class is thread safe, and may be populated by calling 
   parseGadrasMaterialFile() or parseG4MaterialFile() from a secondary client
   thread *assuming* only a single input file will be used to populate the 
   database.  If you are going to populate a second source file, then do not
   call any other member in parrelel during parsing.
*/
  
public:
  static const Material sm_voidMaterial;

public:
  MaterialDB();
  ~MaterialDB();

  //parseGadrasMaterialFile(...)
  //Will thow std::runtime_exception on failure or unexpected input formating.
  //See data/MaterialGadras.lib or GADRAS users manual for example formatting.
  //Ignores all lines of file before the "Air" material
  //text after a '!' character will be ignored for descriptions and density
  //  lines
  //If swapNameDescription==true, then the material name and descriptions will
  //  be swapped (I like this better from the GADRAS original material DB)
  void parseGadrasMaterialFile( const std::string &file,
                                const SandiaDecay::SandiaDecayDataBase *db,
                                const bool swapNameDescription = true );


  //parseG4MaterialFile(...)
  //Will thow std::runtime_exception on failure or unexpected input formating.
  //See data/material.txt for example formatting.
  void parseG4MaterialFile( const std::string &file,
                            const SandiaDecay::SandiaDecayDataBase *db );


  //writeGadrasStyleMaterialFile(...):
  //  Writes the materials back out to a GADRAS compatible format, using Windows
  //  line endings. Writes the material "Air" first in file, to be compatible
  //  with parseGadrasMaterialFile(...)
  //Has not been checked to work with GADRAS, just with this class
  void writeGadrasStyleMaterialFile( std::ostream &file ) const;

  //names(): will throw exception if the database has not been initialized, and
  //  the waite of 1 seconds elapesd for initialization in a background
  //  thread failed.
  const std::vector<std::string> &names() const;

  //material(...) is not case sensitive
  //  For the case of elements, the names have been altered to include both
  //  the symbol and the name, for instance "Fe (iron)" - in this case a
  //  material will be returned if it matches the string exactly
  //  e.g. "Fe (iron)", matches just the symbol, "Fe", or just the part in
  //  paranthesis, "iron".  If for some reason this ambiguous, the first
  //  material found in the sourted alphabetical listing, will be returned.
  //  throw std::runtime_error if it cant find material.
  //  Will throw exception if the database has not been initialized, and
  //  the waite of 1 seconds elapesd for initialization in a background
  //  thread failed.
  const Material *material( const std::string &name ) const;

  //material(...)
  //  throws std::runtime_error if index is out of bounds
  const Material *material( const size_t index ) const;


  //parseChemicalFormula(...): parses chemical formulas using
  //  Material::parseChemicalFormula(...).
  //  Once the furmula is parsed, the Material is created and added to the
  //  database and returned.
  //  If formula is the name of an existing Material, then exception is thrown.
  //  If the formula cant be parsed, than an exception is thrown.
  const Material *parseChemicalFormula( std::string formula,
                                   const SandiaDecay::SandiaDecayDataBase *db );


  //Performs a case-insensitive comparison (eg string::operator<) comparison
  //  of names
  static bool less_than_by_name( const Material *lhs, const Material *rhs );
  
  //Similar to #less_than_by_name
  static bool equal_by_name( const Material *lhs, const Material *rhs );

protected:
  //refreshMaterialNames(): should be called whenever m_materials is modified.
  //  Takes lock on m_mutex, so make sure m_mutex isnt already locked (will
  //  throw exception if cant imediately get it).
  void refreshMaterialNames();

  //wait_material_ready(): returns if the database (e.g. m_materials) has been
  //  initialized.  If it has not been initalized, and this funciton has not
  //  failed before, then will wait up to 1.0 second for the database to be
  //  initialized before returning false.
  bool wait_material_ready() const;
  
protected:

  std::vector<const Material *> m_materials;
  std::vector<std::string> m_materialNames;  //one-to-one correspondance to m_materials
  
  //m_initFailure: if a timeout while waiting for the database to be initialized
  //  has occured, then skip the waiting for subsequent calls.
  mutable bool m_initFailure;
  mutable std::mutex m_mutex;
  mutable std::condition_variable m_condition;
};//class MaterialDB


#endif//MaterialDB_h
