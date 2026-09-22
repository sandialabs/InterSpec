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

#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <istream>

class MaterialDB;

namespace rapidxml
{
  template<class Ch> class xml_node;
}//namespace rapidxml

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

  /** Returns the mass-fraction of the material a certain element is.
   */
  double massFractionOfElementInMaterial( const SandiaDecay::Element * const element ) const;

  float massWeightedAtomicNumber() const;

  /** Version of the XML written by #toXml.
   Change log:
   - 20260920, version 0: initial version.
   */
  static const int sm_xmlSerializationVersion;

  /** Serializes the full definition of this material (name, description, density, and
   elemental/nuclide composition) as a `<MaterialDefinition>` child of `parent`, so the material
   can later be reconstructed without the `MaterialDB` - e.g., after the user has changed its
   density, or if the material is removed from the database in a future version.

   Density is written in g/cm3; mass fractions are of the whole material.

   Returns the created node.  Throws std::exception on error.
   */
  rapidxml::xml_node<char> *toXml( rapidxml::xml_node<char> *parent ) const;

  /** Creates a material from a `<MaterialDefinition>` node written by #toXml.

   Element and nuclide symbols are resolved through `db`.

   Throws std::exception on any problem.
   */
  static std::shared_ptr<const Material> fromXml( const rapidxml::xml_node<char> *node,
                                                  const SandiaDecay::SandiaDecayDataBase *db );

  /** Exact comparison of name, description, density, and composition (#source is provenance
   metadata, and is not compared).
   */
  bool operator==( const Material &rhs ) const;
  bool operator!=( const Material &rhs ) const;

  /** Returns true if the two materials have the same name and identical element and nuclide
   components (same order, exact fractions) - i.e., they are the same material, differing at
   most in density or description.  Use to decide if a material change is only a density edit.
   */
  static bool sameComposition( const Material &lhs, const Material &rhs );

  /** Throws std::runtime_error describing the first difference, if the materials differ in
   name, description, composition, or density (relative tolerance 1E-6).
   */
  static void equalEnough( const Material &lhs, const Material &rhs );

  std::string name;
  std::string description;
  float density;  //in units of PhysicalUnits
  MaterialDefintionsSrc source;

  /** Nuclides and their fractions, by weight.

   The fraction is the fraction of the entire materials density, not just the element.

   The nuclides here are independent from the elements in #elements; normally if an
   isotope is specified here, then the other isotopes of that element will also be given,
   and the element will not be in #elements (but if it is, then you would add the
   mass-fractions).
   */
  std::vector< NuclideFractionPair >  nuclides;

  /** Elements and their fractions by weight.

   Isotopic compositions should be assumed to be natural.
   */
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
  /** MaterialDB: a singleton database of materials.

   The singleton is initialized at application startup via `initialize()`, and accessed
   via `instance()`.  All public query methods are const, making the singleton thread-safe
   after initialization.

   Materials are stored as `shared_ptr<const Material>`, so callers receive a ref-counted
   handle that remains valid regardless of the singleton's lifetime.
   */

public:
  /** Returns the singleton MaterialDB instance.

   If not yet initialized, will attempt lazy initialization.
   Throws std::runtime_error if initialization has not been done and fails.
   */
  static std::shared_ptr<const MaterialDB> instance();

  /** Explicitly initializes the singleton.

   Should be called at application startup, after DecayDataBaseServer::initialize().
   Parses `staticDataDirectory()/MaterialDataBase.txt`, then optionally
   `writableDataDirectory()/MaterialDataBase.txt` to allow user overrides.

   Thread-safe; subsequent calls after successful init are no-ops.
   */
  static void initialize();

  /** Returns whether the singleton has been successfully initialized. */
  static bool initialized();

  /** Returns the initialization error message, or empty string if no error. */
  static const std::string &init_error();

  ~MaterialDB();

  //writeGadrasStyleMaterialFile(...):
  //  Writes the materials back out to a GADRAS compatible format, using Windows
  //  line endings. Writes the material "Air" first in file, to be compatible
  //  with parseGadrasMaterialFile(...)
  //Has not been checked to work with GADRAS, just with this class
  void writeGadrasStyleMaterialFile( std::ostream &file ) const;

  /** Returns the names of all materials in the database.
   */
  const std::vector<std::string> &names() const;

  /** Returns all materials.
   */
  const std::vector<std::shared_ptr<const Material>> &materials() const;

  /** Returns the material with the given name (case-insensitive).

   For the case of elements, the names have been altered to include both
   the symbol and the name, for instance "Fe (iron)" - in this case a
   material will be returned if it matches the string exactly
   e.g. "Fe (iron)", matches just the symbol, "Fe", or just the part in
   parenthesis, "iron".  If for some reason this is ambiguous, the first
   material found in the sorted alphabetical listing will be returned.

   Throws std::runtime_error if it cant find material.
   */
  std::shared_ptr<const Material> material( const std::string &name ) const;

  /** Returns the material at the given index.

   Throws std::runtime_error if index is out of bounds.
   */
  std::shared_ptr<const Material> material( const size_t index ) const;

  /** Returns a Material from the given chemical formula.

   The new material is NOT added to the database.

   Useful for one-off uses in the code.

   Throws exception on error.
   */
  static std::shared_ptr<const Material> materialFromChemicalFormula(
    const std::string &formula,
    const SandiaDecay::SandiaDecayDataBase *db );

  /** Looks `text` up as a material name (see #material(const std::string &)), and if that
   fails, tries to parse it as a chemical formula (see #materialFromChemicalFormula).

   This is how user-entered material text should be resolved everywhere, since #material
   throws for unknown names.

   Returns nullptr if `text` is neither a material name nor a chemical formula, or if the
   singleton is not available.  Never throws.
   */
  static std::shared_ptr<const Material> materialFromNameOrFormula(
    const std::string &text,
    const SandiaDecay::SandiaDecayDataBase *db );

  /** Resolves a serialized material.

   If `definition_node` (a `<MaterialDefinition>` written by Material::toXml) is non-null and
   parses, that definition is used - unless the database has a material that is
   `Material::equalEnough` to it, in which case the databases shared instance is returned
   instead (stable address, no duplicated memory).  Otherwise, or if the definition fails to
   parse, `name` is resolved with #materialFromNameOrFormula.

   Returns nullptr if the material could not be resolved.  Never throws.
   */
  static std::shared_ptr<const Material> materialFromDefinitionOrName(
    const rapidxml::xml_node<char> *definition_node,
    const std::string &name,
    const SandiaDecay::SandiaDecayDataBase *db );

  //Performs a case-insensitive comparison (eg string::operator<) comparison
  //  of names
  static bool less_than_by_name( const Material *lhs, const Material *rhs );

  //Similar to #less_than_by_name
  static bool equal_by_name( const Material *lhs, const Material *rhs );

private:
  MaterialDB();

  //parseGadrasMaterialFile(...)
  //Will throw std::runtime_exception on failure or unexpected input formatting.
  //See data/MaterialGadras.lib or GADRAS users manual for example formatting.
  //Ignores all lines of file before the "Air" material.
  //text after a '!' character will be ignored for descriptions and density lines.
  //If swapNameDescription==true, then the material name and descriptions will
  //  be swapped (I like this better from the GADRAS original material DB).
  //If overwrite==true, materials with the same name as existing ones will replace them.
  void parseGadrasMaterialFile( const std::string &file,
                                const SandiaDecay::SandiaDecayDataBase *db,
                                const bool swapNameDescription = true,
                                const bool overwrite = false );


  //parseG4MaterialFile(...)
  //Will throw std::runtime_exception on failure or unexpected input formatting.
  //See data/material.txt for example formatting.
  void parseG4MaterialFile( const std::string &file,
                            const SandiaDecay::SandiaDecayDataBase *db );

  //refreshMaterialNames(): should be called whenever m_materials is modified.
  void refreshMaterialNames();

  std::vector<std::shared_ptr<const Material>> m_materials;
  std::vector<std::string> m_materialNames;  //one-to-one correspondence to m_materials
};//class MaterialDB


#endif//MaterialDB_h
