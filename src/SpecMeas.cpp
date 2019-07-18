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
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

#include "external_libs/SpecUtils/3rdparty/rapidxml/rapidxml.hpp"
#include "external_libs/SpecUtils/3rdparty/rapidxml/rapidxml_utils.hpp"
#include "external_libs/SpecUtils/3rdparty/rapidxml/rapidxml_print.hpp"

#include <Wt/WServer>
#include <Wt/WIOService>

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "SpecUtils/UtilityFunctions.h"
#include "SpecUtils/SpectrumDataStructs.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace Wt;
using namespace std;

const int SpecMeas::sm_peakXmlSerializationVersion = 2;

namespace
{
  std::string xml_value( const rapidxml::xml_base<char> * const node )
  {
    return (node && node->value_size())
    ? string(node->value(),node->value()+node->value_size())
    : string("");
  }//xml_value(...)
  
  
  /** source is node to clone, and target is the node that will become
      equivalent of source (its name/attribs/children will become equivalent
      to source node).
   20180616: untested
   */
  template <class Ch>
  void clone_node_deep( const rapidxml::xml_node<Ch> *source, rapidxml::xml_node<Ch> *target )
  {
    using namespace rapidxml;
    
    assert( source );
    assert( target );
    
    //auto doc = ((target->type() == node_document) ? static_cast<xml_document<Ch> *>(target) : target->document());
    auto doc = target->document();
    
    target->remove_all_nodes();
    target->remove_all_attributes();
    target->type( source->type() );
    
    Ch *nametxt = doc->allocate_string( source->name(), source->name_size()+1 );
    nametxt[source->name_size()] = 0;    //prob not nec, but JIC
    
    Ch *valuetxt = doc->allocate_string( source->value(), source->value_size()+1 );
    valuetxt[source->value_size()] = 0;  //prob not nec, but JIC
    
    target->name( nametxt, source->name_size() );
    target->value( valuetxt, source->value_size() );
    
    // Clone child nodes and attributes
    for( const xml_attribute<Ch> *attr = source->first_attribute(); attr; attr = attr->next_attribute() )
    {
      Ch *attnametxt = doc->allocate_string( attr->name(), attr->name_size()+1 );
      attnametxt[attr->name_size()] = 0;  //jic
      
      Ch *attvaluetxt = doc->allocate_string( attr->value(), attr->value_size()+1 );
      attvaluetxt[attr->value_size()] = 0; //jic
      
      auto cloned_attrib = doc->allocate_attribute(attnametxt, attvaluetxt, attr->name_size(), attr->value_size());
      target->append_attribute( cloned_attrib );
    }//for( loop over and clone attributes )
    
    for( const xml_node<Ch> *child = source->first_node(); child; child = child->next_sibling() )
    {
      xml_node<Ch> *cloned_child = doc->allocate_node( child->type() );
      target->append_node( cloned_child );
      clone_node_deep( child, cloned_child );
    }
  }//clone_node_deep(...)

}//namespace

SpecMeas::SpecMeas()
  : MeasurementInfo()
{
  m_peaks.reset( new SampleNumsToPeakMap() );
  
//  m_detector.reset( new DetectorPeakResponse() );
  m_displayType.reset( new SpectrumType(kForeground) );
  m_displayedSampleNumbers.reset( new set<int>() );
  
  m_fileWasFromInterSpec = false;
}


SpecMeas::SpecMeas( const SpecMeas &rhs )
  : MeasurementInfo( rhs )
{
  assert(0);
}


const SpecMeas &SpecMeas::operator=( const SpecMeas &rhs )
{
  assert(0);
  return *this;
}//operator=

#if( PERFORM_DEVELOPER_CHECKS )

namespace
{
typedef std::deque< std::shared_ptr<const PeakDef> > PeakDeque;
typedef std::shared_ptr< PeakDeque >                 PeakDequeShrdPtr;
typedef std::map<std::set<int>, PeakDequeShrdPtr >     SampleNumsToPeakMap;


void testPeakMapEqual( const SampleNumsToPeakMap &orig_lhsmap,
                       const SampleNumsToPeakMap &orig_rhsmap )
{
  char buffer[1024];
  
  SampleNumsToPeakMap lhsmap = orig_lhsmap;
  SampleNumsToPeakMap rhsmap = orig_rhsmap;
  
  //Lets get rid of sets with no peaks
  vector< std::set<int> > lhstoremove, rhstoremove;
  for( auto i : lhsmap )
    if( !i.second || i.second->empty() )
      lhstoremove.push_back( i.first );
  for( auto i : rhsmap )
    if( !i.second || i.second->empty() )
      rhstoremove.push_back( i.first );
  for( auto i : lhstoremove )
    lhsmap.erase( lhsmap.find(i) );
  for( auto i : rhstoremove )
    rhsmap.erase( rhsmap.find(i) );
  
  if( lhsmap.size() != rhsmap.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "SpecMeas: Number of sets of peaks for LHS (%i) doesnt match RHS (%i)",
             int(lhsmap.size()), int(rhsmap.size()) );
    throw runtime_error( buffer );
  }

  std::map<std::set<int>, PeakDequeShrdPtr >::const_iterator lhsiter, rhsiter;
  int peakset = 0;
  for( lhsiter = lhsmap.begin(), rhsiter = rhsmap.begin();
      lhsiter != lhsmap.end(); ++lhsiter, ++rhsiter, ++peakset )
  {
    const set<int> &lhsset = lhsiter->first;
    const set<int> &rhsset = rhsiter->first;
    if( lhsset != rhsset )
    {
      stringstream msg;
      msg << "SpecMeas: Sample numbers for peak set " << peakset << " dont match; rhs={";
      for( int t : lhsset )
        msg << t << ", ";
      msg << "}, rhs={";
      for( int t : rhsset )
        msg << t << ", ";
      msg << "}";
      throw runtime_error( msg.str() );
    }
    
    std::shared_ptr< const deque< std::shared_ptr<const PeakDef> > > lhspeaks = lhsiter->second;
    std::shared_ptr< const deque< std::shared_ptr<const PeakDef> > > rhspeaks = rhsiter->second;
    
    if( (!lhspeaks) != (!rhspeaks) )
    {
      snprintf( buffer, sizeof(buffer), "SpecMeas: Peaks for set %i avaialblity for LHS (%s)"
               " doesnt match RHS (%s)",
               peakset,
               (!lhspeaks?"missing":"available"),
               (!rhspeaks?"missing":"available") );
      throw runtime_error(buffer);
    }
    
    if( !!lhspeaks )
    {
      if( lhspeaks->size() != rhspeaks->size() )
      {
        snprintf( buffer, sizeof(buffer),
                 "SpecMeas: Number of peaks for set %i for LHS (%i) doesnt match RHS (%i)",
                 peakset, int(lhspeaks->size()), int(rhspeaks->size()) );
        throw runtime_error( buffer );
      }
      
      deque< std::shared_ptr<const PeakDef> >::const_iterator lhspeakiter, rhspeakiter;
      
      int peakn = 0;
      for( lhspeakiter = lhspeaks->begin(), rhspeakiter = rhspeaks->begin();
          lhspeakiter != lhspeaks->end(); ++lhspeakiter, ++rhspeakiter, ++peakn )
      {
        std::shared_ptr<const PeakDef> lpeak = *lhspeakiter;
        std::shared_ptr<const PeakDef> rpeak = *rhspeakiter;
        
        try
        {
          if( !lpeak )
            throw runtime_error( "Invalid LHS peak pts" );
          if( !rpeak )
            throw runtime_error( "Invalid RHS peak pts" );
          
          PeakDef::equalEnough( *lpeak, *rpeak );
        }catch( std::exception &e )
        {
          stringstream msg;
          msg << "SpecMeas: For peak " << peakn << " of peak set " << peakset << ": "
          << e.what();
          throw runtime_error( msg.str() );
        }//try / catch
      }//for( loop over peaks )
    }//if( !!lhspeaks )
  }//for( loop over sets of peaks )
}//void testPeakMapEqual( const SampleNumsToPeakMap &lhs, const SampleNumsToPeakMap &rhs )
}


void SpecMeas::equalEnough( const SpecMeas &lhs, const SpecMeas &rhs )
{
  MeasurementInfo::equalEnough( lhs, rhs );
  
  char buffer[1024];
  if( (!lhs.m_peaks) != (!rhs.m_peaks) )
  {
    snprintf( buffer, sizeof(buffer), "SpecMeas: Peaks avaialblity for LHS (%s)"
             " doesnt match RHS (%s)",
             (!lhs.m_peaks?"missing":"available"),
             (!rhs.m_peaks?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  if( !!lhs.m_peaks )
  {
    const map<set<int>, PeakDequeShrdPtr > &lhsmap = *lhs.m_peaks;
    const map<set<int>, PeakDequeShrdPtr > &rhsmap = *rhs.m_peaks;
    try
    {
      testPeakMapEqual( lhsmap, rhsmap );
    }catch( std::exception &e )
    {
      throw runtime_error( string(e.what()) + string(", for m_peaks") );
    }
  }//if( !!lhs.m_peaks )
  
  try
  {
    testPeakMapEqual( lhs.m_autoSearchPeaks, rhs.m_autoSearchPeaks );
  }catch( std::exception &e )
  {
    throw runtime_error( string(e.what()) + string(", for m_autoSearchPeaks") );
  }
  
//  try
//  {
//    testPeakMapEqual( lhs.m_autoSearchInitialPeaks, rhs.m_autoSearchInitialPeaks );
//  }catch( std::exception &e )
//  {
//    throw runtime_error( string(e.what()) + string(", for m_autoSearchInitialPeaks") );
//  }
  
  
  if( !!lhs.m_detector != !!rhs.m_detector)
  {
    snprintf( buffer, sizeof(buffer), "SpecMeas: Detector avaialblity for LHS (%s)"
             " doesnt match RHS (%s)",
             (!lhs.m_detector?"missing":"available"),
             (!rhs.m_detector?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  if( !!lhs.m_detector )
    DetectorPeakResponse::equalEnough( *lhs.m_detector, *rhs.m_detector );
  
  if( !!lhs.m_displayType != !!rhs.m_displayType )
  {
    snprintf( buffer, sizeof(buffer), "SpecMeas: display type avaialblity for LHS (%s)"
             " doesnt match RHS (%s)",
             (!lhs.m_displayType?"missing":"available"),
             (!rhs.m_displayType?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  if( !!lhs.m_displayType && ( (*lhs.m_displayType) != (*rhs.m_displayType) ) )
  {
    snprintf( buffer, sizeof(buffer), "SpecMeas: display type for LHS (%i)"
             " doesnt match RHS (%i)",
             int(*lhs.m_displayType), int(*rhs.m_displayType) );
    throw runtime_error(buffer);
  }

  if( !!lhs.m_displayedSampleNumbers != !!rhs.m_displayedSampleNumbers )
  {
    snprintf( buffer, sizeof(buffer), "SpecMeas: displayed sample numbers avaialblity for LHS (%s)"
             " doesnt match RHS (%s)",
             (!lhs.m_displayedSampleNumbers?"missing":"available"),
             (!rhs.m_displayedSampleNumbers?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  
  if( !!lhs.m_displayedSampleNumbers && ( (*lhs.m_displayedSampleNumbers) != (*rhs.m_displayedSampleNumbers) ) )
  {
    stringstream msg;
    msg << "SpecMeas: diplayed sample numbers of LHS ({";
    for( int i : (*lhs.m_displayedSampleNumbers) )
      msg << i << ",";
    msg << "}) doesnt match RHS ({";
    for( int i : (*rhs.m_displayedSampleNumbers) )
      msg << i << ",";
    msg << "})";
    throw runtime_error( msg.str() );
  }
  
  
  if( (!lhs.m_shieldingSourceModel) != (!rhs.m_shieldingSourceModel) )
  {
    stringstream msg;
    msg << "SpecMeas: availability of shieldingSourceModel of LHS ("
    << (lhs.m_shieldingSourceModel ? "" : "not ") << "available)"
    << "doesnt match RHS (" << (rhs.m_shieldingSourceModel ? "" : "not ")
    << "available)";
    throw runtime_error( msg.str() );
  }
  
  if( lhs.m_shieldingSourceModel )
  {
    //ToDo: make a proper comparison by traversing nodes and comparing values.
    string lhsdata, rhsdata;
    rapidxml::print( std::back_inserter(lhsdata), *lhs.m_shieldingSourceModel, 0 );
    rapidxml::print( std::back_inserter(rhsdata), *rhs.m_shieldingSourceModel, 0 );
    if( lhsdata != rhsdata )
    {
      stringstream msg;
      msg << "The ShieldingSOurceModel of the LHS does not exactly match RHS;"
      " this could be a harmless error, or an actual issue - Will has not"
      " implemented a proper comparison yet:\n\tLHS=" << lhsdata << "\n\tRHS="
      << rhsdata << "\n";
      throw runtime_error( msg.str() );
    }
  }
}//void SpecMeas::equalEnough( const SpecMeas &lhs, const SpecMeas &rhs )
#endif //#if( PERFORM_DEVELOPER_CHECKS )



void SpecMeas::uniqueCopyContents( const SpecMeas &rhs )
{
  std::lock( mutex_, rhs.mutex_ );
  std::lock_guard<std::recursive_mutex> lhs_lock( mutex_, std::adopt_lock_t() );
  std::lock_guard<std::recursive_mutex> rhs_lock( rhs.mutex_, std::adopt_lock_t() );

  if( &rhs == this )
    return;

  MeasurementInfo::operator=( rhs );

  if( !m_peaks )
    m_peaks.reset( new SampleNumsToPeakMap );
  else
    m_peaks->clear();
  m_autoSearchPeaks.clear();
//  m_autoSearchInitialPeaks.clear();

  typedef std::shared_ptr<const PeakDef> PeakShrdPtr;
  if( rhs.m_peaks )
  {
    for( const SampleNumsToPeakMap::value_type &vt : *(rhs.m_peaks) )
    {
      (*m_peaks)[vt.first] = std::make_shared<PeakDeque>();
      for( PeakShrdPtr peak : *(vt.second) )
        (*m_peaks)[vt.first]->push_back( peak );
    }//for( const SampleNumsToPeakMap::value_type &vt : *(rhs.m_peaks) )
  }//if( rhs.m_peaks )
  
  for( const SampleNumsToPeakMap::value_type &vt : m_autoSearchPeaks )
  {
    m_autoSearchPeaks[vt.first] = std::make_shared<PeakDeque>();
    for( PeakShrdPtr peak : *(vt.second) )
      m_autoSearchPeaks[vt.first]->push_back( peak );
  }//for( const SampleNumsToPeakMap::value_type &vt : *(rhs.m_peaks) )
  
//  for( const SampleNumsToPeakMap::value_type &vt : m_autoSearchInitialPeaks )
//  {
//    m_autoSearchPeaks[vt.first] = std::make_shared<PeakDeque>();
//    for( PeakShrdPtr peak : *(vt.second) )
//      m_autoSearchPeaks[vt.first]->push_back( peak );
//  }//for( const SampleNumsToPeakMap::value_type &vt : *(rhs.m_peaks) )

  if( m_detector != rhs.m_detector )
  {
    if( !!rhs.m_detector )
    {
      m_detector.reset( new DetectorPeakResponse() );
      *m_detector = *rhs.m_detector;
    }else
    {
      m_detector.reset();
    }
  }//if( m_detector != rhs.m_detector )

  *m_displayType = *rhs.m_displayType;
  *m_displayedSampleNumbers = *rhs.m_displayedSampleNumbers;
  
  
  if( rhs.m_shieldingSourceModel && rhs.m_shieldingSourceModel->first_node())
  {
    m_shieldingSourceModel.reset( new rapidxml::xml_document<char>() );
    clone_node_deep( rhs.m_shieldingSourceModel.get(), m_shieldingSourceModel.get() );
  }else
    m_shieldingSourceModel.reset();
  
}//void uniqueCopyContents( const SpecMeas &rhs )


SpecMeas::~SpecMeas()
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );

//  cerr << "About to delete " << filename()
//       << ", connected=" << m_aboutToBeDeleted.isConnected()
//       << endl;

  m_aboutToBeDeleted.emit();

//  bool has_unique = false;
//  for( const MeasurementShrdPtr &p : measurements_ )
//  {
//    if( p && p->gamma_counts_ && p->gamma_counts_.unique() )
//    {
//      has_unique = true;
//      break;
//    }//if( the gamma_counts is unique )
//  }//if( measurements_.size() && measurements_[0] )

//  if( has_unique )
//    cerr << "\nSpecMeas destructing unique copy of '" << filename_ << "'" << endl;
}//SpecMeas destructor


void SpecMeas::save2012N42FileInSlaveThread( std::shared_ptr<SpecMeas> info,
                                         const std::string filename,
                                         boost::function<void()> error_callback )
{
  if( !info && error_callback )
  {
    error_callback();
    return;
  }//if( !info && error_callback )
  
  boost::function<void()> worker = boost::bind( &SpecMeas::save2012N42File, info, filename, error_callback );
  
  WServer *server = WServer::instance();  //can this ever be NULL?
  if( server )
  {
    server->ioService().post( worker );
  }else
  {
    //Probably wont ever get here - but just incase
    cerr << SRC_LOCATION << "\n\tWarning: couldnt get WServer - not good" << endl << endl;
    worker();
  }//if( server ) / else
}//save2012N42FileInSlaveThread(...)


bool SpecMeas::save2012N42File( const std::string &filename )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
#ifdef _WIN32
  const std::wstring wfilename = UtilityFunctions::convert_from_utf8_to_utf16(filename);
  ofstream ofs( wfilename.c_str(), ios::binary|ios::out );
#else
  ofstream ofs( filename.c_str(), ios::binary|ios::out );
#endif

  return ofs.is_open() && write_2012_N42( ofs );
}//bool save2012N42File( const std::string &filename );


void SpecMeas::save2012N42File( const std::string &filename,
                              boost::function<void()> error_callback )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  const bool status = save2012N42File( filename );
  if( !status && error_callback )
    error_callback();
}//void save2012N42File(...)


const char *SpecMeas::toString( const XmlPeakSource source )
{
  switch( source )
  {
    case UserPeaks:                   return "UserPeaks";
    case AutomatedSearchPeaks:        return "AutomatedSearchPeaks";
//    case AutomatedSearchInitialPeaks: return "AutomatedSearchInitialPeaks";
  }//switch( source )
  
  assert( 0 );
  return "";
}

SpecMeas::XmlPeakSource SpecMeas::xmlPeakSourceFromString( const string &val )
{
  if( val == "UserPeaks" )
    return UserPeaks;
//  if( val == "AutomatedSearchPeaks" )
//    return AutomatedSearchPeaks;
  if( val != "AutomatedSearchPeaks" )
    throw runtime_error( string(val) + " is not a valid SpecMeas::XmlPeakSource" );
  return AutomatedSearchPeaks;
  
//  if( val != "AutomatedSearchInitialPeaks" )
//    throw runtime_error( string(val) + " is not a valid SpecMeas::XmlPeakSource" );
  
//  return AutomatedSearchPeaks;
}//xmlPeakSourceFromString


void SpecMeas::addPeaksToXmlHelper( const SpecMeas::SampleNumsToPeakMap &inpeaks,
                    const SpecMeas::XmlPeakSource source,
                    rapidxml::xml_node<char> *peaksnode,
                    std::map<std::shared_ptr<PeakContinuum>,int> &continuums,
                    std::map<std::shared_ptr<const PeakDef>,int> &peakids )
{
  using namespace rapidxml;
  
  /* ToDo: use N42 MeasurementGroupReferences to match peaks up to their
           respective sample numbers (and also for analysis results and source
           shielding model), so we can a) be more N42 compliant, but b) use a
           a construct that isnt as fragile as sample number (which isnt really
           a thing in N42 files).
   */
  
  rapidxml::xml_document<char> *doc = peaksnode->document();
  
  const char *typeval = toString( source );
  
  
  for( SampleNumsToPeakMap::const_iterator iter = inpeaks.begin();
      iter != inpeaks.end(); ++iter )
  {
    const std::set<int> &nums = iter->first;
    const PeakDequeShrdPtr &peaks = iter->second;
    if( !peaks || peaks->empty() )
      continue;
    
    vector<int> thesepeakids;
    for( PeakDeque::const_iterator i = peaks->begin(); i != peaks->end(); ++i )
    {
      std::map<std::shared_ptr<const PeakDef>,int>::const_iterator pos;
      pos = peakids.find( *i );
      if( pos != peakids.end() )
      {
        thesepeakids.push_back( pos->second );
      }else
      {
        const int index = static_cast<int>( peakids.size() + 1 );
        peakids[*i] = index;
        thesepeakids.push_back( index );
        rapidxml::xml_node<char> *node = (*i)->toXml( peaksnode, peaksnode, continuums );
        
        if( node )
        {
          const string idstr = std::to_string( index );
          const char *value = doc->allocate_string( idstr.c_str() );
          xml_attribute<char> *attr = doc->allocate_attribute( "id", value );
          node->append_attribute( attr );
        }//if( node )
      }//
    }//for( PeakDeque::const_iterator i = peaks->begin(); i != peaks->end(); ++i )
    
    assert( thesepeakids.size() == peaks->size() );
    
    xml_node<char> *peaksset = doc->allocate_node( node_element, "PeakSet" );
    peaksnode->append_node( peaksset );
    
    xml_attribute<char> *typeatt = doc->allocate_attribute( "source", typeval );
    peaksset->append_attribute( typeatt );
    
    stringstream samples, thesepeakidstr;
    for( set<int>::const_iterator i = nums.begin(); i != nums.end(); ++i )
      samples << (i==nums.begin()?"": " ") << *i;
    
    for( vector<int>::const_iterator i = thesepeakids.begin(); i != thesepeakids.end(); ++i )
      thesepeakidstr << (i==thesepeakids.begin()?"": " ") << *i;
    
    const char *val = doc->allocate_string( samples.str().c_str() );
    xml_node<char> *samplenode = doc->allocate_node( node_element, "SampleNumbers", val );
    peaksset->append_node( samplenode );
    
    val = doc->allocate_string( thesepeakidstr.str().c_str() );
    xml_node<char> *peakidnode = doc->allocate_node( node_element, "PeakIds", val );
    peaksset->append_node( peakidnode );
  }//for( iterate over peaks
}//void addPeaksToXmlHelper(...)


void SpecMeas::addPeaksToXml( rapidxml::xml_node<char> *peaksnode ) const
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( !m_peaks )
    return;
  
  rapidxml::xml_document<char> *doc = peaksnode->document();
  
  rapidxml::xml_attribute<char> *versionatt = peaksnode->first_attribute( "version", 7 );
  char peakxmlverstr[12];
  snprintf( peakxmlverstr, sizeof(peakxmlverstr), "%i", sm_peakXmlSerializationVersion );
  
  if( !versionatt )
  {
    const char *value = doc->allocate_string( peakxmlverstr );
    versionatt = doc->allocate_attribute( "version", value );
    peaksnode->append_attribute( versionatt );
  }else
  {
    if( xml_value(versionatt) != peakxmlverstr )
      throw runtime_error( "SpecMeas::addPeaksToXml: Invalid peak XML version" );
  }
  
  std::map<std::shared_ptr<const PeakDef>,int> peakids;
  std::map<std::shared_ptr<PeakContinuum>,int> continuumids;
  addPeaksToXmlHelper( *m_peaks, UserPeaks, peaksnode, continuumids, peakids );
  addPeaksToXmlHelper( m_autoSearchPeaks, AutomatedSearchPeaks, peaksnode, continuumids, peakids );
//  addPeaksToXmlHelper( m_autoSearchInitialPeaks, AutomatedSearchInitialPeaks, peaksnode, continuumids, peakids );
}//void addPeaksToXml(...)



void SpecMeas::addPeaksFromXml( const ::rapidxml::xml_node<char> *peaksnode )
{
  using namespace rapidxml;
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( !m_peaks )
    m_peaks.reset( new SampleNumsToPeakMap() );
  
  int version = 1;
  const xml_attribute<char> *version_att = peaksnode->first_attribute("version",7);
  const std::string versionstr = xml_value( version_att );
  if( versionstr.size() )
  {
    if( !(stringstream(versionstr) >> version) )
      cerr << "Failed to extract integer version from Peak element; '"
           << versionstr << "' is not an integer.  Assuming version 1" << endl;
  }//if( versionstr.size() )
  
  std::map<int,std::shared_ptr<PeakContinuum> > continuums;
  for( const xml_node<char> *node = peaksnode->first_node("PeakContinuum",13);
      node; node = node->next_sibling("PeakContinuum",13) )
  {
    std::shared_ptr<PeakContinuum> continuum( new PeakContinuum() );
    int contId;
    continuum->fromXml( node, contId );
    continuums[contId] = continuum;
  }//for( loop over PeakSet nodes )
  
  if( version < 2 )
  {
    if( peaksnode->first_node("Peak") )
      throw runtime_error( "There can not be any Peak nodes under the Peaks"
                           " node for version 1 (or unversioned) XML." );
    
    for( const xml_node<char> *node = peaksnode->first_node("PeakSet",7);
         node; node = node->next_sibling("PeakSet",7) )
    {
      xml_node<char> *samples_node = node->first_node( "SampleNumbers", 13 );
      if( !samples_node || !samples_node->value() )
        throw runtime_error( "Did not find SampleNumbers under PeakSet" );
    
      //We'll assume a perfect convertion of int<-->float, which if the input
      //  data is well formed, should be true
      std::vector<float> contents;
      UtilityFunctions::split_to_floats( samples_node->value(), samples_node->value_size(), contents );
      set<int> samplenums;
      for( const float t : contents )
        samplenums.insert( samplenums.end(), static_cast<int>(t) );
    
      PeakDequeShrdPtr peaks = std::make_shared<PeakDeque>();
    
      for( const xml_node<char> *peak_node = node->first_node("Peak",4);
          peak_node; peak_node = peak_node->next_sibling("Peak",4) )
      {
        std::shared_ptr<PeakDef> peak( new PeakDef() );
        peak->fromXml( peak_node, continuums );
        peaks->push_back( peak );
      }//for( loop over Peak nodes )
    
      (*m_peaks)[samplenums] = peaks;
    }//for( loop over PeakSet nodes )
  }else if( version < 3 )
  {
    std::map<int,std::shared_ptr<const PeakDef> > id_to_peak;
    
    //Go through and read in all peaks and continuums
    for( const xml_node<char> *peak_node = peaksnode->first_node("Peak",4);
        peak_node; peak_node = peak_node->next_sibling("Peak",4) )
    {
      std::shared_ptr<PeakDef> peak( new PeakDef() );
      peak->fromXml( peak_node, continuums );
      
      const xml_attribute<char> *idatt = peak_node->first_attribute( "id" );
      const string peakidstr = xml_value( idatt );
      
      if( peakidstr.empty() )
        throw runtime_error( "Peak elements must have a \"id\" attribute." );
      
      int peakid;
      if( !(stringstream(peakidstr) >> peakid) )
        throw runtime_error( "Peak \"id\" attribute must be numerical: '"
                             + peakidstr + "' invalid" );
      if( id_to_peak.count(peakid) )
        throw runtime_error( "Peak elements \"id\" attribute must have a unique value: " + peakidstr );
      
      id_to_peak[peakid] = peak;
    }//for( loop over Peak nodes )
    
    
    for( const xml_node<char> *node = peaksnode->first_node("PeakSet",7);
        node; node = node->next_sibling("PeakSet",7) )
    {
      PeakDequeShrdPtr peaks = std::make_shared<PeakDeque>();
      const string srcstr = xml_value( node->first_attribute( "source", 6 ) );
      const XmlPeakSource source = xmlPeakSourceFromString( srcstr );
      
      xml_node<char> *samples_node = node->first_node( "SampleNumbers", 13 );
      if( !samples_node || !samples_node->value() )
        throw runtime_error( "Did not find SampleNumbers under PeakSet" );
      
      xml_node<char> *peakid_node = node->first_node( "PeakIds", 7 );
      if( !peakid_node || !peakid_node->value_size() )
        throw runtime_error( "Did not find PeakIds under PeakSet" );
      
      std::vector<int> samplenumvec, peakids;
      UtilityFunctions::split_to_ints( samples_node->value(), samples_node->value_size(), samplenumvec );
      UtilityFunctions::split_to_ints( peakid_node->value(), peakid_node->value_size(), peakids );
      set<int> samplenums;
      for( const float t : samplenumvec )
        samplenums.insert( samplenums.end(), static_cast<int>(t) );
      
      for( const int peakid : peakids )
      {
        const map<int,std::shared_ptr<const PeakDef> >::const_iterator iter
                                                    = id_to_peak.find( peakid );
        if( iter == id_to_peak.end() )
          throw runtime_error( "Could not find peak with id '"
                               + std::to_string(peakid)
                               + "' in the XML for sample numbers '"
                               + xml_value(samples_node) );
        peaks->push_back( iter->second );
      }//for( const int peakid : peakids )
      
      //The problem here is that the sample numbers can change; what we need to
      //  do is create a <RadMeasurementGroup> and use this to associate peaks
      //  with measurements.  We should also do the same for analysis, and also
      //  dispalyed sample numbers
      
      switch( source )
      {
        case UserPeaks:
          (*m_peaks)[samplenums] = peaks;
        break;
          
        case AutomatedSearchPeaks:
          m_autoSearchPeaks[samplenums] = peaks;
        break;
          
//        case AutomatedSearchInitialPeaks:
//          m_autoSearchInitialPeaks[samplenums] = peaks;
//        break;
      }//switch( source )
    }//for( loop over peaksets )
    
  }else
  {
    throw runtime_error( "Invalid version attribute for the Peak element" );
  }
}//void addPeaksFromXml( ::rapidxml::xml_node<char> *peaksnode );


rapidxml::xml_document<char> *SpecMeas::shieldingSourceModel()
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  return m_shieldingSourceModel.get();
}


void SpecMeas::setShieldingSourceModel( std::unique_ptr<rapidxml::xml_document<char>> &&model )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( !model && !m_shieldingSourceModel )
    return;
  
  bool is_diff = true;
  if( m_shieldingSourceModel && model && !modified_ )
  {
    //TODO: go through and compare nodes to see if they are actually different.
    //      for now, just do a string compare
    string lhsdata, rhsdata;
    rapidxml::print( std::back_inserter(lhsdata), *m_shieldingSourceModel, 0 );
    rapidxml::print( std::back_inserter(rhsdata), *model, 0 );
    is_diff = (lhsdata != rhsdata);
  }//
  
  m_shieldingSourceModel = std::move( model );
  
  if( is_diff )
    modified_ = modifiedSinceDecode_ = true;
}//void setShieldingSourceModel( std::unique_ptr<rapidxml::xml_document<char>> &&model )



bool SpecMeas::write_2006_N42( std::ostream &ostr ) const
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  {
    stringstream temporary;
    MeasurementInfo::write_2006_N42( temporary );
    string content = temporary.str();
    size_t pos = content.find( "</N42InstrumentData>" );
    if( pos == string::npos )
    {
      cerr << "Failed to find N42InstrumentData closing tag" << endl;
      ostr << content;
      return false;
    }
    
    ostr << content.substr( 0, pos );
  }
  
  
  {//begin codeblock to write SpecMeas specific stuff
    using namespace rapidxml;
    xml_document<char> doc;
    
    xml_node<char> *peaksnode = doc.allocate_node( node_element, "Peaks" );
    doc.append_node( peaksnode );
    
    xml_node<char> *interspecnode = appendSpecMeasStuffToXml( &doc );

    string data;
    rapidxml::print( std::back_inserter(data), *interspecnode, 0 );
    
    ostr << data << "\r\n";
  }//end codeblock to write SpecMeas specific stuff
  
  ostr << "</N42InstrumentData>" << "\r\n";
  
  return !ostr.bad();
}//bool write_2006_N42( std::ostream &ostr ) const




bool SpecMeas::write_iaea_spe( std::ostream &output,
                            std::set<int> sample_nums,
                            const std::set<int> &det_nums ) const
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  map<set<int>, PeakDequeShrdPtr>::const_iterator peakiter;
  
  if( m_peaks )
    peakiter = m_peaks->find(sample_nums);
  
  //If we dont have additional peak information to write to the file, we can
  //  just call MeasurementInfo's version of this function and return.
  if( !m_peaks || peakiter == m_peaks->end() || !peakiter->second || peakiter->second->empty() )
    return MeasurementInfo::write_iaea_spe( output, sample_nums, det_nums );
  
  
  //To convert from peak energy to channel number, we need the same energy
  //  calibration as the spectrum uses
  if( sample_nums.empty() )
    sample_nums = sample_numbers_;
  
  const size_t ndet = detector_numbers_.size();
  vector<bool> detectors( ndet, true );
  if( !det_nums.empty() )
  {
    for( size_t i = 0; i < ndet; ++i )
      detectors[i] = (det_nums.count(detector_numbers_[i]) != 0);
  }//if( det_nums.empty() )
  
  MeasurementShrdPtr summed = sum_measurements( sample_nums, detectors );
  if( !summed )
    return MeasurementInfo::write_iaea_spe( output, sample_nums, det_nums );
  
  
  //Call MeasurementInfo's version of this function, and get rid of last line,
  //  then append peak information.
  stringstream strm;
  if( !MeasurementInfo::write_iaea_spe( strm, sample_nums, det_nums ) )
    return false;
  
  //Lets get rid of the trailing: "$ENDRECORD:\r\n"
  //const stringstream::pos_type pos = strm.tellg();  //seems to give zero
  string str = strm.str();
  if( str.size() > 13 )
    str = str.substr( 0, str.size() - 13 );
  output << str;
  
  const deque< std::shared_ptr<const PeakDef> > &peaks = *peakiter->second;
  
  //I think PEAKLABELS is a PeakEasy specific addition
  output << "$PEAKLABELS:\r\n";
  
  for( deque< std::shared_ptr<const PeakDef> >::const_iterator iter = peaks.begin();
      iter != peaks.end(); ++iter )
  {
    const PeakDef &peak = **iter;
    
    if( peak.userLabel().size() || peak.gammaParticleEnergy() > 0.0f )
    {
      //channel should actually be a float (well, thats what Peak Easy uses),
      //  but for current ease of implementation will just truncate
      size_t channel = 0;
      try
      {
        channel = summed->find_gamma_channel( peak.gammaParticleEnergy() );
      }catch(...)
      {
        continue;
      }
      
      string label = peak.userLabel();
      if( peak.parentNuclide() )
        label += " " + peak.parentNuclide()->symbol;
      if( peak.xrayElement() )
        label += " " + peak.xrayElement()->symbol + " xray";
      if( peak.reaction() )
        label += " " + peak.reaction()->name();
      
      switch( peak.sourceGammaType() )
      {
        case PeakDef::NormalGamma:                            break;
        case PeakDef::AnnihilationGamma: label += " annih."; break;
        case PeakDef::SingleEscapeGamma: label += " s.e.";   break;
        case PeakDef::DoubleEscapeGamma: label += " d.e.";   break;
        case PeakDef::XrayGamma:         label += " xray";   break;
      }//switch( peak.sourceGammaType() )
      
      UtilityFunctions::ireplace_all( label, "\r\n", " " );
      UtilityFunctions::ireplace_all( label, "\r", " " );
      UtilityFunctions::ireplace_all( label, "\n", " " );
      UtilityFunctions::ireplace_all( label, "\"", "&quot;" );
      UtilityFunctions::ireplace_all( label, "  ", " " );
      
      if( label.size() && label[0]==' ' )
        label = label.substr( 1, label.size() - 1 );
      
      if( label.size() )
        output << channel << ", \"" << label << "\"\r\n";
    }//if( there is a label for the peak )
  }//for( loop over peaks )
  
   //$ROI: This group contains the regions of interest marked in the spectrum.
   //      The first line the number of regions, the following lines contain the
   //      start and stop channels for each region.
  
  output << "$ENDRECORD:\r\n";
  
  return true;
}//write_iaea_spe...



void SpecMeas::decodeSpecMeasStuffFromXml( const ::rapidxml::xml_node<char> *interspecnode )
{
  using namespace rapidxml;
  using ::rapidxml::internal::compare;
  
  if( !interspecnode || !interspecnode->value() 
     || !compare( interspecnode->name(), interspecnode->name_size(), "DHS:InterSpec", 13, false ) )
     throw std::logic_error( "SpecMeas::decodeSpecMeasStuffFromXml(...): invalid input" );
  
  m_displayedSampleNumbers.reset( new set<int>() );
  const xml_node<char> *node = node = interspecnode->first_node( "DisplayedSampleNumbers", 22 );
  if( node && node->value() )
  {
    std::vector<float> contents;
    UtilityFunctions::split_to_floats( node->value(), node->value_size(), contents );
    for( const float t : contents )
      m_displayedSampleNumbers->insert( static_cast<int>(t) );
  }//if( node )
  
  
  m_displayType.reset( new SpectrumType(kForeground) );
  node = interspecnode->first_node( "DisplayType", 11 );
  if( node && node->value() )
  {
    for( SpectrumType t = SpectrumType(0); 
         t <= kBackground; t = SpectrumType(t+1) )
    {
      const char *val = descriptionText( t );
      const size_t len = strlen( val );
      if( compare( node->value(), node->value_size(), val, len, false ) )
      {
        *m_displayType = t;
        break;
      }
    }//for( loop over SpectrumType );   
  }//if( node && node->value() )
  
  
  node = interspecnode->first_node( "Peaks", 5 );
  if( node )
    addPeaksFromXml( node );

  
  node = interspecnode->first_node( "DetectorPeakResponse", 20 );
  if( node )
  {
    m_detector.reset( new DetectorPeakResponse() );
    m_detector->fromXml( node );
  }else
    m_detector.reset();
  
  node = interspecnode->first_node( "ShieldingSourceFit", 18 );
  if( node )
  {
    m_shieldingSourceModel.reset( new rapidxml::xml_document<char>() );
    auto model_node = m_shieldingSourceModel->allocate_node(rapidxml::node_element);
    m_shieldingSourceModel->append_node( model_node );
    clone_node_deep( node, model_node );
  }else
    m_shieldingSourceModel.reset();
}//void decodeSpecMeasStuffFromXml( ::rapidxml::xml_node<char> *parent )

::rapidxml::xml_node<char> *SpecMeas::appendSampleNumbersToXml(
                                    ::rapidxml::xml_node<char> *interspec_node ) const
{
  using namespace rapidxml;
  
  rapidxml::xml_document<char> *doc = interspec_node->document();
  
  if( !!m_displayedSampleNumbers )
  {
    vector<int> samples( m_displayedSampleNumbers->begin(), m_displayedSampleNumbers->end() );
    stringstream dispSamples;
    for( size_t i = 0; i < samples.size(); ++i )
      dispSamples << (i?" ":"") << samples[i];
    
    const char *val = doc->allocate_string( dispSamples.str().c_str() );
    xml_node<char> *node = doc->allocate_node( node_element, "DisplayedSampleNumbers", val );
    interspec_node->append_node( node );
    
    return node;
  }//if( !!m_displayedSampleNumbers )
  
  return 0;
}//appendSampleNumbersToXml(...)


::rapidxml::xml_node<char> *SpecMeas::appendSpecMeasStuffToXml( 
                                    ::rapidxml::xml_node<char> *RadInstrumentData ) const
{
  using namespace rapidxml;
  
  rapidxml::xml_document<char> *doc = RadInstrumentData->document();
  xml_node<char> *interspec_node = doc->allocate_node( node_element, "DHS:InterSpec" );
  RadInstrumentData->append_node( interspec_node );
  
  appendSampleNumbersToXml( interspec_node );
  
  if( !!m_displayType )
  {
    const char *val = descriptionText( *m_displayType );
    xml_node<char> *node = doc->allocate_node( node_element, "DisplayType", val );
    interspec_node->append_node( node );
  }
  
  
  if( !!m_peaks && m_peaks->size() )
  {
    std::map<std::shared_ptr<PeakContinuum>,int> continuums;
    xml_node<char> *peaksnode = doc->allocate_node( node_element, "Peaks" );
    interspec_node->append_node( peaksnode );
    
    addPeaksToXml( peaksnode );
  }//if( !!m_peaks )
  
  if( !!m_detector && m_detector->isValid() )
    m_detector->toXml( interspec_node, doc );

  //Shielding source fit - if applicable.
  if( m_shieldingSourceModel )
  {
    xml_node<char> *modelnode = doc->allocate_node( node_element );
    interspec_node->append_node( modelnode );
    clone_node_deep( m_shieldingSourceModel->first_node(), modelnode );
  }//if( m_shieldingSourceModel && m_shieldingSourceModel->first_node() )
  
  
  return interspec_node;
}//appendSpecMeasStuffToXml(...);


std::shared_ptr< ::rapidxml::xml_document<char> > SpecMeas::create_2012_N42_xml() const
{
  using namespace rapidxml;
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  std::shared_ptr< xml_document<char> > doc = MeasurementInfo::create_2012_N42_xml();
  
  if( !doc )
    return doc;
  
  xml_node<char> *RadInstrumentData = doc->first_node( "RadInstrumentData", 17 );
  
  if( !RadInstrumentData )
    throw std::logic_error( "SpecMeas::create_2012_N42_xml failed to get RadInstrumentData node" );
  
  appendSpecMeasStuffToXml( RadInstrumentData );
  
  return doc;
}//create_2012_N42_xml() const


bool SpecMeas::load_from_N42( std::istream &input )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    ::rapidxml::file<char> input_file( input );
    
#if( RAPIDXML_USE_SIZED_INPUT_WCJOHNS )
    return SpecMeas::load_N42_from_data( input_file.data(), input_file.data()+input_file.size() );
#else
    input_file.check_for_premature_nulls( 2048, ' ' );
    return SpecMeas::load_N42_from_data( input_file.data() );
#endif
  }catch( std::exception & )
  {    
    reset();
    return false;
  }//try/catch
  
  return true;
}//bool load_from_N42( char *data )


bool SpecMeas::load_N42_file( const std::string &filename )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    std::vector<char> data;
    UtilityFunctions::load_file_data( filename.c_str(), data );
    
#if( RAPIDXML_USE_SIZED_INPUT_WCJOHNS )
    const bool loaded = SpecMeas::load_N42_from_data( &data.front(), (&data.front()) + data.size() );
#else
    input_file.check_for_premature_nulls( 2048, ' ' );
    const bool loaded = SpecMeas::load_N42_from_data( &data.front() );
#endif
    
    if( !loaded )
      throw runtime_error( "!loaded" );
    
    filename_ = filename;
  }catch( std::exception &e )
  {
    cout << endl << e.what() << endl;
    reset();
    return false;
  }//try/catch

  return true;
}//bool load_N42_file( const std::string &filename )


void SpecMeas::load_N42_from_doc( rapidxml::xml_document<char> &doc )
{
  const rapidxml::xml_node<char> *document_node = doc.first_node();
  const rapidxml::xml_node<char> *interspecnode = document_node ? document_node->first_node( "DHS:InterSpec", 13 ) : nullptr;
  
  //See notes for m_fileWasFromInterSpec.
  m_fileWasFromInterSpec = !!interspecnode;
  
  const bool parsed = load_from_N42_document( document_node );
  
  if( !parsed )
    throw runtime_error( "Couldnt Parse" );
  
  if( interspecnode )
    decodeSpecMeasStuffFromXml( interspecnode );
}//void load_N42_from_doc( rapidxml::xml_document<char> &doc )


bool SpecMeas::load_N42_from_data( char *data )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  reset();
  
  if( !is_candidate_n42_file(data) )
    return false;
  
  try
  {
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_trim_whitespace | rapidxml::allow_sloppy_parse; //rapidxml::parse_normalize_whitespace
    doc.parse<flags>( data );
    load_N42_from_doc( doc );
  }catch( rapidxml::parse_error &e )
  {
    ptrdiff_t predatalen = e.where<char>() - data;
    if( predatalen < 0 )
      predatalen = 0;
    if( predatalen > 256 )
      predatalen = 256;
    const size_t postlen = strnlen( e.where<char>(), 256 );
    
    m_fileWasFromInterSpec = false;
    
    const string predata( e.where<char>() - predatalen, e.where<char>() + 1 );
    const string postdata( e.where<char>(), e.where<char>() + postlen );
    
    if( (e.where<char>()-data) > 0 )
    {
      cerr << "SpecMeas::load_N42_from_data() caught parse_error: " << e.what()
      << " at location " << (e.where<char>()-data)
      << ", and pre-text: " << predata<< "\n\nAnd post text: " << postdata
      << endl << endl;
    }//if( (e.where<char>()-data) > 0 )
    
    reset();
    return false;
  }catch( std::exception & )
  {
    m_fileWasFromInterSpec = false;
    //cerr << "SpecMeas::load_N42_from_data() caught: " << e.what() << endl;
    
    reset();
    return false;
  }//try/catch
  
  return true;
}//bool load_N42_from_data( char *data )


#if( RAPIDXML_USE_SIZED_INPUT_WCJOHNS )
bool SpecMeas::load_N42_from_data( char *data, char *data_end )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  reset();
  
  data_end = convert_n42_utf16_xml_to_utf8( data, data_end );
  
  if( !is_candidate_n42_file(data,data_end) )
  {
    cout << "is_candidate_n42_file" << endl;
    return false;
  }
  
  rapidxml::xml_document<char> doc;
  try
  {
    const int flags = rapidxml::parse_trim_whitespace | rapidxml::allow_sloppy_parse; //rapidxml::parse_normalize_whitespace
    doc.parse<flags>( data, data_end );
    
    load_N42_from_doc( doc );
  }catch( rapidxml::parse_error &e )
  {
    ptrdiff_t predatalen = e.where<char>() - data;
    if( predatalen < 0 )
      predatalen = 0;
    if( predatalen > 256 )
      predatalen = 256;
    const size_t postlen = std::min( static_cast<size_t>(data_end - e.where<char>()), static_cast<size_t>(256) );
    
    m_fileWasFromInterSpec = false;
    
    const string predata( e.where<char>() - predatalen, e.where<char>() + 1 );
    const string postdata( e.where<char>(), e.where<char>() + postlen );
    
    if( (e.where<char>()-data) > 0 )
    {
      cerr << "SpecMeas::load_N42_from_data() caught parse_error: " << e.what()
      << " at location " << (e.where<char>()-data)
      << ", and pre-text: " << predata<< "\n\nAnd post text: " << postdata
      << endl << endl;
    }//if( (e.where<char>()-data) > 0 )
    
    reset();
    return false;
  }catch( std::exception &e )
  {
    cerr << "Caught: " << e.what() << endl;
    
   rapidxml::print( cout, doc, 0 );
    
    cerr << "Caught: " << e.what() << endl;
    
    m_fileWasFromInterSpec = false;
    
    reset();
    return false;
  }//try/catch
  
  return true;
}//bool load_N42_from_data( char *data )
#endif  // RAPIDXML_USE_SIZED_INPUT_WCJOHNS
      

std::shared_ptr<DetectorPeakResponse> SpecMeas::detector()
{
  return m_detector;
}

std::shared_ptr<const DetectorPeakResponse> SpecMeas::detector() const
{
  return m_detector;
}

void SpecMeas::setDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );

  if( det != m_detector )
    modified_ = modifiedSinceDecode_ = true;

  m_detector = det;
}//void setDetector( std::shared_ptr<const DetectorPeakResponse> );

SpectrumType SpecMeas::displayType() const
{
  if( !m_displayType )
  {
    cerr << SRC_LOCATION << "\n\tThere is a serious error here - definetly need"
         << " investigating!!!" << endl;
    return kForeground;
  }
  return *m_displayType;
}


const std::set<int> &SpecMeas::displayedSampleNumbers() const
{
  static const set<int> emptySet;
  if( !m_displayedSampleNumbers )
    return emptySet;
  return *m_displayedSampleNumbers;
}

void SpecMeas::detectorChangedCallback( std::shared_ptr<DetectorPeakResponse> det )
{
//  m_detector = det;
  setDetector( det );
}

std::set<std::set<int> > SpecMeas::sampleNumsWithPeaks() const
{
  std::set<std::set<int> > answer;
  if( !m_peaks )
    return answer;
  for( const SampleNumsToPeakMap::value_type &t : *m_peaks )
    answer.insert( t.first );
  return answer;
}//sampleNumsWithPeaks() const


void SpecMeas::removeAllPeaks()
{
  if( !!m_peaks )
    m_peaks->clear();
}//void removeAllPeaks()


void SpecMeas::setPeaks( const std::deque< std::shared_ptr<const PeakDef> > &peakdeque,
              const std::set<int> &samplenums )
{
  if( !m_peaks )
    m_peaks.reset( new SampleNumsToPeakMap() );
  (*m_peaks)[samplenums].reset( new PeakDeque( peakdeque ) );
}//void setPeaks(...)


std::shared_ptr< const SpecMeas::PeakDeque > SpecMeas::automatedSearchPeaks( const std::set<int> &samplenums ) const
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  SampleNumsToPeakMap::const_iterator pos = m_autoSearchPeaks.find( samplenums );
  if( pos != m_autoSearchPeaks.end() )
    return pos->second;
  return std::shared_ptr< const SpecMeas::PeakDeque >();
}


//std::shared_ptr< const SpecMeas::PeakDeque > SpecMeas::automatedSearchInitialPeaks( const std::set<int> &samplenums )
//{
//  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
//
//  SampleNumsToPeakMap::const_iterator pos = m_autoSearchInitialPeaks.find( samplenums );
//  if( pos != m_autoSearchPeaks.end() )
//    return pos->second;
//  return std::shared_ptr< const SpecMeas::PeakDeque >();
//}


void SpecMeas::setAutomatedSearchPeaks( const std::set<int> &samplenums,
                                       std::shared_ptr< PeakDeque > peaks
                                      /*, std::shared_ptr< PeakDeque > intitalPeaks*/ )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
#if( PERFORM_DEVELOPER_CHECKS )
  for( const int samplenum : samplenums )
  {
    if( !sample_numbers_.count(samplenum) )
    {
      char msg[512];
      snprintf( msg, sizeof(msg),
                "sample %i is not a valid sample number", samplenum );
      log_developer_error( BOOST_CURRENT_FUNCTION, msg );
    }//if( invalid sample number )
  }//for( const int samplenum : samplenums )
#endif
  
  m_autoSearchPeaks[samplenums] = peaks;
//  m_autoSearchInitialPeaks;[samplenums] = intitalPeaks;
  
  setModified();
}//setAutomatedSearchPeaks(...)


std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > >
                             SpecMeas::peaks( const std::set<int> &samplenums )
{
  if( !m_peaks )
    throw runtime_error( "SpecMeas::peaks(...): serious logic error" );
  
  if( m_peaks->find(samplenums) == m_peaks->end() )
    (*m_peaks)[samplenums].reset( new deque< std::shared_ptr<const PeakDef> >() );

//  cerr << "\nThere are " << (*m_peaks).size() << " sets of peaks; this set has "
//       << (*m_peaks)[samplenums]->size() << " peaks" << endl;
  
  return (*m_peaks)[samplenums];
}//PeakDequeShrdPtr peaks( const std::set<int> samplenums )


std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > >
                        SpecMeas::peaks( const std::set<int> &samplenums ) const
{
  if( !m_peaks )
    throw runtime_error( "SpecMeas::peaks(...): serious logic error" );
  
  SampleNumsToPeakMap::const_iterator pos = m_peaks->find(samplenums);
  if( pos == m_peaks->end() )
    return PeakDequeShrdPtr();
  return pos->second;
}//PeakDequeShrdPtr peaks( const std::set<int> samplenums )


void SpecMeas::displayedSpectrumChangedCallback( SpectrumType type,
                                       std::shared_ptr<SpecMeas> measurment,
                                       std::set<int> sample_numbers )
{
  if( measurment.get() != this )
    return;

  if( (*m_displayedSampleNumbers) != sample_numbers
      && ((*m_displayType) == type) )
  {
    modified_ = true;
    *m_displayedSampleNumbers = sample_numbers;
    *m_displayType = type;
  }
}

DetectorType SpecMeas::guessDetectorTypeFromFileName( std::string name )
{
  UtilityFunctions::to_lower( name );
  const bool detective = UtilityFunctions::contains( name, "detective" );
  const bool ex = (UtilityFunctions::contains( name, "detectiveex" )
                      || UtilityFunctions::contains( name, "detective-ex" )
                      || UtilityFunctions::contains( name, "detective ex" ) );
  const bool oneHundred = UtilityFunctions::contains( name, "100" );
  const bool ex100 = UtilityFunctions::contains( name, "ex-100" );
  const bool micro = UtilityFunctions::contains( name, "micro" );
  
  if( (detective && oneHundred) || ex100 )
    return kDetectiveEx100Detector;
  if( detective && ex )
    return kDetectiveExDetector;
  if( detective && micro )
    return kMicroDetectiveDetector;
  if( detective )
    return kDetectiveDetector;
  
  //Note: identiFINDER-NG are much more common that first generation
  //      identiFINDERs, so well just always assume the its an NG
  if( UtilityFunctions::contains( name, "identifinder" ) )
  {
    if( UtilityFunctions::contains( name, "labr" ) )
      return kIdentiFinderLaBr3Detector;
    return kIdentiFinderNGDetector; //kIdentiFinderDetector;
  }//if( UtilityFunctions::contains( name, "identifinder" ) )
  
//  if( UtilityFunctions::contains( name, "identifinder" )
//      && UtilityFunctions::contains( name, "ng" ) )
//    return kIdentiFinderNGDetector;
  
  if( UtilityFunctions::contains( name, "gr135" )
     || UtilityFunctions::contains( name, "gr-135" ) )
    return kGR135Detector;
  
  if( UtilityFunctions::contains( name, "falcon" ) )
     return kFalcon5000;
  
  
  return kUnknownDetector;
}

void SpecMeas::setModified()
{
  modified_ = modifiedSinceDecode_ = true;
}

void SpecMeas::cleanup_after_load( const unsigned int flags )
{
  if( m_fileWasFromInterSpec )
  {
    MeasurementInfo::cleanup_after_load( (flags | MeasurementInfo::DontChangeOrReorderSamples) );
  }else
  {
    MeasurementInfo::cleanup_after_load( flags );
  }

  //should detect if the detector was loaded, and if not, if we know the type,
  //  we could then load it.
  //   -instead for right now lets do this in InterSpec...
}//void SpecMeas::cleanup_after_load()


Wt::Signal<> &SpecMeas::aboutToBeDeleted()
{
  return m_aboutToBeDeleted;
}



void shiftPeaksHelper( std::map<std::set<int>, std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > > > &input,
                       map< std::shared_ptr<const PeakDef>, std::shared_ptr<const PeakDef> > &shiftedPeaks,
                       set< std::shared_ptr<const PeakContinuum> > &shiftedContinuums,
                       const std::vector<float> &old_pars,
                       const std::vector< std::pair<float,float> > &old_devpairs,
                       const Measurement::EquationType old_eqn_type,
                       const std::vector<float> &new_pars,
                       const std::vector< std::pair<float,float> > &new_devpairs,
                       const Measurement::EquationType new_eqn_type,
                       const size_t nbins )
{
  typedef std::shared_ptr<const PeakDef> PeakPtr;
  typedef std::deque< std::shared_ptr<const PeakDef> > PeakDeque;
  typedef std::shared_ptr< PeakDeque >                 PeakDequeShrdPtr;
  typedef std::map<std::set<int>, PeakDequeShrdPtr > SampleNumsToPeakMap;
  
  for( SampleNumsToPeakMap::value_type &vt : input )
  {
    PeakDequeShrdPtr peakdeque = vt.second;
    if( !peakdeque )
      continue;
    
    PeakDeque newpeaks;
    
    for( PeakPtr peak : *peakdeque )
    {
      if( !peak )
        continue;
      
      const map< PeakPtr, PeakPtr >::iterator pos = shiftedPeaks.find( peak );
      if( pos != shiftedPeaks.end() )
      {
        newpeaks.push_back( pos->second );
        continue;
      }
      
      std::shared_ptr<PeakDef> newpeak( new PeakDef(*peak) );
      
      const bool shiftCont = !shiftedContinuums.count( newpeak->continuum() );
      shiftedContinuums.insert( newpeak->continuum() );
      
      SpecMeas::translatePeakForCalibrationChange( *newpeak, old_pars, old_devpairs,
                                        old_eqn_type, new_pars, new_devpairs,
                                        new_eqn_type, nbins, shiftCont );
      
      shiftedPeaks[peak] = newpeak;
      newpeaks.push_back( newpeak );
    }//for( PeakPtr peak : *peakdeque )
    
    peakdeque->swap( newpeaks );
  }//for( SampleNumsToPeakMap::value_type : *m_peaks )
}//shiftPeaksHelper


void SpecMeas::shiftPeaksForRecalibration( std::vector<float> old_pars,
                                const std::vector< std::pair<float,float> > &old_devpairs,
                                Measurement::EquationType old_eqn_type,
                                std::vector<float> new_pars,
                                const std::vector< std::pair<float,float> > &new_devpairs,
                                Measurement::EquationType new_eqn_type )
{
  typedef std::shared_ptr<const PeakDef> PeakPtr;
  
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );

  const size_t nbins = num_gamma_channels();
  if( !nbins )
    return;
  
  map< PeakPtr, PeakPtr > shiftedPeaks;
  set< std::shared_ptr<const PeakContinuum> > shiftedContinuums;
  
  if( !!m_peaks )
    shiftPeaksHelper( *m_peaks,
                      shiftedPeaks, shiftedContinuums, old_pars, old_devpairs,
                     old_eqn_type, new_pars, new_devpairs, new_eqn_type, nbins );
  shiftPeaksHelper( m_autoSearchPeaks,
                   shiftedPeaks, shiftedContinuums, old_pars, old_devpairs,
                   old_eqn_type, new_pars, new_devpairs, new_eqn_type, nbins );
//  shiftPeaksHelper( m_autoSearchInitialPeaks,
//                   shiftedPeaks, shiftedContinuums, old_pars, old_devpairs,
//                   old_eqn_type, new_pars, new_devpairs, new_eqn_type, nbins );
}//SpecMeas::shiftPeaksForRecalibration(...)


void SpecMeas::translatePeakForCalibrationChange( PeakDef &peak,
                                          std::vector<float> old_pars,
                                          const std::vector< std::pair<float,float> > &old_devpairs,
                                          Measurement::EquationType old_eqn_type,
                                          std::vector<float> new_pars,
                                          const std::vector< std::pair<float,float> > &new_devpairs,
                                          Measurement::EquationType new_eqn_type,
                                          const size_t nbins,
                                          const bool translate_continuum )
{
//#if( PERFORM_DEVELOPER_CHECKS )
//  const double preGausArea = peak.gauss_integral(peak.lowerX(), peak.upperX() );
//  const double preContArea = peak.continuum()->offset_integral(peak.lowerX(), peak.upperX());
//#endif
  
  if( old_pars.size() < 2 || new_pars.size() < 2 )
    throw runtime_error( "translatePeakForCalibrationChange: invalid num coefs" );
  
  if( old_eqn_type==Measurement::LowerChannelEdge
      || old_eqn_type==Measurement::InvalidEquationType
      || new_eqn_type==Measurement::LowerChannelEdge
      || new_eqn_type==Measurement::InvalidEquationType )
    throw runtime_error( "translatePeakForCalibrationChange() can only handle"
                         " Polynomial or FullRangeFraction Calibrations" );
  
  if( old_eqn_type == Measurement::Polynomial || old_eqn_type == Measurement::UnspecifiedUsingDefaultPolynomial )
      old_pars = polynomial_coef_to_fullrangefraction( old_pars, nbins );
  if( new_eqn_type == Measurement::Polynomial || new_eqn_type == Measurement::UnspecifiedUsingDefaultPolynomial )
      new_pars = polynomial_coef_to_fullrangefraction( new_pars, nbins );
  
  
  if( !peak.gausPeak() )
  {
    if( !translate_continuum )
      return;
      
    const float oldMean = static_cast<float>(peak.mean());
    const float oldlow = static_cast<float>(peak.lowerX());
    const float oldhigh = static_cast<float>(peak.upperX());
    const float meanbin = find_bin_fullrangefraction( oldMean, old_pars, nbins, old_devpairs, 0.001f );
    const float lowbin = find_bin_fullrangefraction( oldlow, old_pars, nbins, old_devpairs, 0.001f );
    const float highbin = find_bin_fullrangefraction( oldhigh, old_pars, nbins, old_devpairs, 0.001f );
    const float newMean = fullrangefraction_energy( meanbin, new_pars, nbins, new_devpairs );
    const float newLower = fullrangefraction_energy( lowbin, new_pars, nbins, new_devpairs );
    const float newUpper = fullrangefraction_energy( highbin, new_pars, nbins, new_devpairs );
    
    peak.continuum()->setRange( newLower, newUpper );
    peak.set_coefficient( newMean, PeakDef::Mean );
    return;
  }//if( !peak.gausPeak() )
  
  
  const float oldMean = static_cast<float>(peak.mean());
  const float oldSigma = static_cast<float>(peak.sigma());
  const float oldbin = find_bin_fullrangefraction( oldMean, old_pars,
                                                  nbins, old_devpairs, 0.001f );
  const float newMean = fullrangefraction_energy( oldbin, new_pars, nbins, new_devpairs );
  
  const float oldneg2sigmabin = find_bin_fullrangefraction( oldMean - 2.0*oldSigma,
                                      old_pars, nbins, old_devpairs, 0.001f );
  const float oldpos2sigmabin = find_bin_fullrangefraction( oldMean + 2.0*oldSigma,
                                      old_pars, nbins, old_devpairs, 0.001f );
  const float newneg2sigma = fullrangefraction_energy( oldneg2sigmabin, new_pars, nbins, new_devpairs );
  const float newpos2sigma = fullrangefraction_energy( oldpos2sigmabin, new_pars, nbins, new_devpairs );
  
  float strech = 0.25f*(newpos2sigma - newneg2sigma) / oldSigma;
  
  if( IsNan(strech) || IsInf(strech) )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    const char * msg = "Found an invalid stretch value when claculated from the"
                       " mean";
    log_developer_error( BOOST_CURRENT_FUNCTION, msg );
#endif
    return;
  }//if( IsNan(strech) || IsInf(strech) )
  
#if( PERFORM_DEVELOPER_CHECKS )
  const float newbin = find_bin_fullrangefraction( newMean, new_pars,
                                                  nbins, new_devpairs, 0.001 );
  if( fabs(newbin - oldbin) > 0.025 )  //0.025 arbitrary
  {
    stringstream msg;
    msg.precision( 9 );
    msg << "When recalibrating from coefs={";
    for( size_t i = 0; i < old_pars.size(); ++i )
      msg << (i ? ", " : " ") << old_pars[i];
    msg << " }";
    if( old_devpairs.size() )
    {
      msg << ", devpairs={";
      for( size_t i = 0; i < old_devpairs.size(); ++i )
        msg << (i ? ", {" : " {") << old_devpairs[i].first
            << "," << old_devpairs[i].second << "}";
      msg << "}";
    }
    
    msg << " to {";
    for( size_t i = 0; i < new_pars.size(); ++i )
      msg << (i ? ", " : " ") << new_pars[i];
    msg << " }";
    if( new_devpairs.size() )
    {
      msg << ", devpairs={";
      for( size_t i = 0; i < new_devpairs.size(); ++i )
        msg << (i ? ", {" : " {") << new_devpairs[i].first << ","
            << new_devpairs[i].second << "}";
      msg << "}";
    }
    msg << ", peak at mean=" << oldMean << ", moved to mean=" << newMean
        << ", but some error caused it to move from bin " << oldbin << " to "
        << newbin << ", which shouldnt have happend (should have stayed same "
        << "bin number).";
    
    log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
  }//if( fabs(newbin - oldbin) > 0.025 )
#endif
  
  peak.setMean( newMean );
  peak.setSigma( strech*oldSigma );
  peak.setMeanUncert( strech*peak.meanUncert() );
  peak.setSigmaUncert( strech * peak.sigmaUncert() );
  
  
  if( !translate_continuum )
    return;
  
  std::shared_ptr<PeakContinuum> continuum = peak.continuum();
  
  if( continuum->energyRangeDefined() )
  {
    const float oldLowEnergy = static_cast<float>(continuum->lowerEnergy());
    const float oldlowbin = find_bin_fullrangefraction( oldLowEnergy, old_pars,
                                                       nbins, old_devpairs, 0.001 );
    const float new_lowenergy = fullrangefraction_energy( oldlowbin, new_pars, nbins, new_devpairs );
    
    
    const float oldHighEnergy = static_cast<float>(continuum->upperEnergy());
    const float oldhighbin = find_bin_fullrangefraction( oldHighEnergy, old_pars,
                                                        nbins, old_devpairs, 0.001 );
    const float new_highenergy = fullrangefraction_energy( oldhighbin, new_pars,
                                                          nbins, new_devpairs );
    
    strech = (new_highenergy - new_lowenergy) / (oldHighEnergy - oldLowEnergy);
    continuum->setRange( new_lowenergy, new_highenergy );
  }//if( peak.continuum().energyRangeDefined() )
  
  if( continuum->isPolynomial() )
  {
    const double oldref = continuum->referenceEnergy();
    
    const float oldrefbin = find_bin_fullrangefraction( oldref, old_pars,
                                                       nbins, old_devpairs, 0.001 );
    const float newref = fullrangefraction_energy( oldrefbin, new_pars, nbins, new_devpairs );
    
    if( !continuum->energyRangeDefined() )
      strech = static_cast<float>( newref / oldref );
    
    if( IsNan(strech) || IsInf(strech) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( BOOST_CURRENT_FUNCTION,
                          "Found an invalid stretch value when calculated "
                           "by the reference energy" );
#endif
      strech = newMean / oldMean;
    }//if( IsNan(strech) || IsInf(strech) )

    
    vector<double> vars = continuum->parameters();
    vector<double> uncerts = continuum->unertainties();

    for( size_t i = 0; i < vars.size(); ++i )
    {
      vars[i] = vars[i] / std::pow( strech, static_cast<float>(i+1.0f) );
      uncerts[i] = uncerts[i] / std::pow( strech, static_cast<float>(i+1.0f) );
    }//for( size_t i = 0; i < contvars.size(); ++i )
    
    continuum->setParameters( newref, vars, uncerts );
  }else if( !!continuum->externalContinuum() )
  {
    //XXX - each peak will now get their own continuum - a huge waste of memorry!
    std::shared_ptr<Measurement> newcont( new Measurement( *continuum->externalContinuum() ) );
    newcont->recalibrate_by_eqn( new_pars, new_devpairs, new_eqn_type );
    continuum->setExternalContinuum( newcont );
  }//if( peak.continuum().isPolynomial() ) / else if(
  
  
  
#if( PERFORM_DEVELOPER_CHECKS )
//  if( end == start )
//    log_developer_error( BOOST_CURRENT_FUNCTION, "data from database wasnt null terminated" );
#endif


  //    if( continuum->defined() )
  //    {
  //      const double delta_energy = new_energy-peak.mean();
  //      peak.offsetConstant -= peak.offsetSlope*delta_energy;
  //    }//if( continuum->defined() )
}//translatePeakForCalibrationChange(...)




