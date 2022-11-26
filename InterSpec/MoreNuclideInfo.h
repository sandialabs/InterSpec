#ifndef MoreNuclideInfo_h
#define MoreNuclideInfo_h
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

#include <map>
#include <string>
#include <memory>
#include <vector>

// Forward declerations
namespace SandiaDecay
{
  struct Nuclide;
}

namespace MoreNuclideInfo
{
  class NucInfo;
  class MoreNucInfoDb;
}//


/** Methods to provide additional nuclide information, supplied by an XML file.
* 
* The information isnt strictly limited to nuclide, but could include 
* reactions, x-rays, or whatever.
* 
* Similar to DecayDataBaseServer, as single global database of information is
* read in and used for all sessions.
* 
* As of 20221125, the additonal information XML file is under developmnet, and
* has not been released with the InterSpec code yet.
*/
namespace MoreNuclideInfo
{
  /** An enum to give the state of the database.  This is primarily useful
  for deciding if you want to use the information from the primary GUI 
  thread, or instead inualize things in a worker thread, and then update 
  the GUI (or not bother at all).
  */
  enum class InfoStatus
  {
    /** No attempt to initualize the database has been made. */
    NotInited,

    /** The database was attempted to be initualized, but failed
    because of missing or invalid XML file. 
    */
    FailedToInit,
    
    /** The database is inualized and ready for use. */
    Inited
  };//enum class InfoStatus


  /** Returns status of the database being inualized. */
  InfoStatus more_nuc_info_db_status();


  /** Class that holds the information in the XML file for each nuclide.
  */
  class NucInfo
  {
  public:
    std::string m_nuclide;
    std::vector<std::string> m_associated;
    std::string m_notes;

    size_t memsize() const;

    NucInfo(){};
    NucInfo( NucInfo && ) = default;
  private:

    friend class MoreNucInfoDb;
  };//struct NucInfo

  
  class RefInfo
  {
  public:
    std::string m_key;
    std::string m_url;
    std::string m_desc;

    /// Make sure move constructor is created
    RefInfo( RefInfo && ) = default;
    RefInfo(){};

    size_t memsize() const;

  private:
    RefInfo( std::string &&key, std::string &&url, std::string &&desc );

    friend class MoreNucInfoDb;
  };//class RefInfo


  class MoreNucInfoDb
  {
  public:
    const NucInfo *info( const std::string &nuc ) const;
    const NucInfo *info(const SandiaDecay::Nuclide *nuc) const;

    /** Returns the information database, inualizing it if needed.

      Returned pointer will be nullptr if XML file is missing, or invalid.
    */
    static std::shared_ptr<const MoreNucInfoDb> instance();
  
    std::map<const SandiaDecay::Nuclide *, NucInfo> m_nuc_infos;
    std::map<std::string, NucInfo> m_other_infos;
    std::map<std::string, RefInfo> m_references;

    /// Make sure move constructor is created
    MoreNucInfoDb( MoreNucInfoDb && ) = default;

    size_t memsize() const;
  
    MoreNucInfoDb();
  private:

    /** Throws exception on failure. */
    void init();
  };//class MoreNucInfoDb

}//namespace MoreNuclideInfo
#endif //MoreNuclideInfo_h
