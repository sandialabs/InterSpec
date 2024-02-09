#include <deque>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include "SpecUtils/3rdparty/rapidxml/rapidxml.hpp"
#include "SpecUtils/3rdparty/rapidxml/rapidxml_utils.hpp"
#include "SpecUtils/3rdparty/rapidxml/rapidxml_print.hpp"

#include <Wt/WServer>
#include <Wt/WIOService>


#include "QLPeakDef.h"
#include "QLSpecMeas.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"


using namespace Wt;
using namespace std;

#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH

const int QLSpecMeas::sm_peakXmlSerializationVersion = 2;

namespace
{
  std::string xml_value( const rapidxml::xml_base<char> * const node )
  {
    return (node && node->value_size())
    ? string(node->value(),node->value()+node->value_size())
    : string("");
  }//xml_value(...)
}//namespace


void log_error_message( const std::string &message, const std::string &source, const int priority )
{
  std::cerr << source << ": " << message << std::endl;
}


QLSpecMeas::QLSpecMeas()
  : SpecUtils::SpecFile()
{
  m_peaks.reset( new SampleNumsToPeakMap() );
  
  m_displayedSampleNumbers.reset( new set<int>() );
}


QLSpecMeas::QLSpecMeas( const QLSpecMeas &rhs )
  : SpecUtils::SpecFile( rhs )
{
  assert(0);
}


const QLSpecMeas &QLSpecMeas::operator=( const QLSpecMeas &rhs )
{
  assert(0);
  return *this;
}//operator=



QLSpecMeas::~QLSpecMeas()
{
}//QLSpecMeas destructor



const char *QLSpecMeas::toString( const XmlPeakSource source )
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

QLSpecMeas::XmlPeakSource QLSpecMeas::xmlPeakSourceFromString( const string &val )
{
  if( val == "UserPeaks" )
    return UserPeaks;
//  if( val == "AutomatedSearchPeaks" )
//    return AutomatedSearchPeaks;
  if( val != "AutomatedSearchPeaks" )
    throw runtime_error( string(val) + " is not a valid QLSpecMeas::XmlPeakSource" );
  return AutomatedSearchPeaks;
  
//  if( val != "AutomatedSearchInitialPeaks" )
//    throw runtime_error( string(val) + " is not a valid QLSpecMeas::XmlPeakSource" );
  
//  return AutomatedSearchPeaks;
}//xmlPeakSourceFromString



void QLSpecMeas::addPeaksFromXml( const ::rapidxml::xml_node<char> *peaksnode )
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
  
  std::map<int,std::shared_ptr<QLPeakContinuum> > continuums;
  for( const xml_node<char> *node = peaksnode->first_node("PeakContinuum",13);
      node; node = node->next_sibling("PeakContinuum",13) )
  {
    std::shared_ptr<QLPeakContinuum> continuum( new QLPeakContinuum() );
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
      SpecUtils::split_to_floats( samples_node->value(), samples_node->value_size(), contents );
      set<int> samplenums;
      for( const float t : contents )
        samplenums.insert( samplenums.end(), static_cast<int>(t) );
    
      PeakDequeShrdPtr peaks = std::make_shared<PeakDeque>();
    
      for( const xml_node<char> *peak_node = node->first_node("Peak",4);
          peak_node; peak_node = peak_node->next_sibling("Peak",4) )
      {
        std::shared_ptr<QLPeakDef> peak( new QLPeakDef() );
        peak->fromXml( peak_node, continuums );
        peaks->push_back( peak );
      }//for( loop over Peak nodes )
    
      (*m_peaks)[samplenums] = peaks;
    }//for( loop over PeakSet nodes )
  }else if( version < 3 )
  {
    std::map<int,std::shared_ptr<const QLPeakDef> > id_to_peak;
    
    //Go through and read in all peaks and continuums
    for( const xml_node<char> *peak_node = peaksnode->first_node("Peak",4);
        peak_node; peak_node = peak_node->next_sibling("Peak",4) )
    {
      try
      {
      std::shared_ptr<QLPeakDef> peak( new QLPeakDef() );
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
      }catch( std::exception &e )
      {
        cerr << "Failed to parse peak: " << e.what() << endl;
      }
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
      SpecUtils::split_to_ints( samples_node->value(), samples_node->value_size(), samplenumvec );
      SpecUtils::split_to_ints( peakid_node->value(), peakid_node->value_size(), peakids );
      set<int> samplenums;
      for( const float t : samplenumvec )
        samplenums.insert( samplenums.end(), static_cast<int>(t) );
      
      for( const int peakid : peakids )
      {
        const map<int,std::shared_ptr<const QLPeakDef> >::const_iterator iter
                                                    = id_to_peak.find( peakid );
        if( iter == id_to_peak.end() )
          throw runtime_error( "Could not find peak with id '"
                               + boost::lexical_cast<string>(peakid)
                               + "' in the XML for sample numbers '"
                               + xml_value(samples_node) );
        peaks->push_back( iter->second );
      }//for( const int peakid, peakids )
      
      
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



void QLSpecMeas::decodeSpecMeasStuffFromXml( const ::rapidxml::xml_node<char> *interspecnode )
{
  using namespace rapidxml;
  using ::rapidxml::internal::compare;
  
  if( !interspecnode || !interspecnode->value() 
     || !compare( interspecnode->name(), interspecnode->name_size(), "DHS:InterSpec", 13, false ) )
     throw std::logic_error( "QLSpecMeas::decodeSpecMeasStuffFromXml(...): invalid input" );
  
  m_displayedSampleNumbers.reset( new set<int>() );
  const xml_node<char> *node = node = interspecnode->first_node( "DisplayedSampleNumbers", 22 );
  if( node && node->value() )
  {
    std::vector<float> contents;
    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    for( const float t : contents )
      m_displayedSampleNumbers->insert( static_cast<int>(t) );
  }//if( node )
  
  /*
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
   */
  
  
  node = interspecnode->first_node( "Peaks", 5 );
  if( node )
  {
    try
    {
      addPeaksFromXml( node );
    }
    catch(const std::exception& e)
    {
      std::cerr << "Failed to parse peaks: " << e.what() << '\n';
    }
  }
}//void decodeSpecMeasStuffFromXml( ::rapidxml::xml_node<char> *parent )



bool QLSpecMeas::load_from_N42( std::istream &input )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    ::rapidxml::file<char> input_file( input );
    return QLSpecMeas::load_N42_from_data( input_file.data() );
  }catch( std::exception & )
  {    
    reset();
    return false;
  }//try/catch
  
  return true;
}//bool load_from_N42( char *data )


bool QLSpecMeas::load_N42_file( const std::string &filename )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    rapidxml::file<char> input_file( filename.c_str() );  //throws runtime_error upon failure
    
    const bool loaded = QLSpecMeas::load_N42_from_data( input_file.data() );
    
    if( !loaded )
      throw runtime_error( "" );
    
    filename_ = filename;
  }catch( std::exception & )
  {
    reset();
    return false;
  }//try/catch

  return true;
}//bool load_N42_file( const std::string &filename )


bool QLSpecMeas::load_N42_from_data( char *data )
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    reset();
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_trim_whitespace; //rapidxml::parse_normalize_whitespace
    doc.parse<flags>( data );
    const rapidxml::xml_node<char> *document_node = doc.first_node();
    
    const bool parsed = load_from_N42_document( document_node );
    
    if( !parsed )
      throw runtime_error( "Couldnt Parse" );
    
    try
    {
      const rapidxml::xml_node<char> *interspecnode = document_node->first_node( "DHS:InterSpec", 13 );
      if( interspecnode )
        decodeSpecMeasStuffFromXml( interspecnode );
    }catch( std::exception &e )
    {
      cerr << "Failed to load InterSpec stuff form XML" << endl;
    }
  }catch( rapidxml::parse_error &e )
  {
    ptrdiff_t predatalen = e.where<char>() - data;
    if( predatalen < 0 )
      predatalen = 0;
    if( predatalen > 256 )
      predatalen = 256;
    const size_t postlen = strnlen( e.where<char>(), 256 );
    
    const string predata( e.where<char>() - predatalen, e.where<char>() + 1 );
    const string postdata( e.where<char>(), e.where<char>() + postlen );
    
    if( (e.where<char>()-data) > 0 )
    {
      cerr << "QLSpecMeas::load_N42_from_data() caught parse_error: " << e.what()
      << " at location " << (e.where<char>()-data)
      << ", and pre-text: " << predata<< "\n\nAnd post text: " << postdata
      << endl << endl;
    }//if( (e.where<char>()-data) > 0 )
    
    reset();
    return false;
  }catch( std::exception &e )
  {
    cerr << "QLSpecMeas::load_N42_from_data() caught: " << e.what() << endl;
    
    reset();
    return false;
  }//try/catch
  
  return true;
}//bool load_N42_from_data( char *data )


const std::set<int> &QLSpecMeas::displayedSampleNumbers() const
{
  static const set<int> emptySet;
  if( !m_displayedSampleNumbers )
    return emptySet;
  return *m_displayedSampleNumbers;
}

std::set<std::set<int> > QLSpecMeas::sampleNumsWithPeaks() const
{
  std::set<std::set<int> > answer;
  if( !m_peaks )
    return answer;
  for( const SampleNumsToPeakMap::value_type &t : *m_peaks )
    answer.insert( t.first );
  return answer;
}//sampleNumsWithPeaks() const


void QLSpecMeas::removeAllPeaks()
{
  if( !!m_peaks )
    m_peaks->clear();
}//void removeAllPeaks()


void QLSpecMeas::setPeaks( std::deque< std::shared_ptr<const QLPeakDef> > &peakdeque,
              const std::set<int> &samplenums )
{
  if( !m_peaks )
    m_peaks.reset( new SampleNumsToPeakMap() );
  (*m_peaks)[samplenums].reset( new PeakDeque( peakdeque ) );
}//void setPeaks(...)


std::shared_ptr< const QLSpecMeas::PeakDeque > QLSpecMeas::automatedSearchPeaks( const std::set<int> &samplenums ) const
{
  std::lock_guard<std::recursive_mutex> scoped_lock( mutex_ );
  
  SampleNumsToPeakMap::const_iterator pos = m_autoSearchPeaks.find( samplenums );
  if( pos != m_autoSearchPeaks.end() )
    return pos->second;
  return std::shared_ptr< const QLSpecMeas::PeakDeque >();
}


std::shared_ptr< std::deque< std::shared_ptr<const QLPeakDef> > >
                             QLSpecMeas::peaks( const std::set<int> &samplenums )
{
  if( !m_peaks )
    throw runtime_error( "QLSpecMeas::peaks(...): serious logic error" );
  
  if( m_peaks->find(samplenums) == m_peaks->end() )
    (*m_peaks)[samplenums].reset( new deque< std::shared_ptr<const QLPeakDef> >() );

//  cerr << "\nThere are " << (*m_peaks).size() << " sets of peaks; this set has "
//       << (*m_peaks)[samplenums]->size() << " peaks" << endl;
  
  return (*m_peaks)[samplenums];
}//PeakDequeShrdPtr peaks( const std::set<int> samplenums )


std::shared_ptr<const std::deque< std::shared_ptr<const QLPeakDef> > >
                        QLSpecMeas::peaks( const std::set<int> &samplenums ) const
{
  if( !m_peaks )
    throw runtime_error( "QLSpecMeas::peaks(...): serious logic error" );
  
  SampleNumsToPeakMap::const_iterator pos = m_peaks->find(samplenums);
  if( pos == m_peaks->end() )
    return PeakDequeShrdPtr();
  return pos->second;
}//PeakDequeShrdPtr peaks( const std::set<int> samplenums )



void QLSpecMeas::setModified()
{
  modified_ = modifiedSinceDecode_ = true;
}

void QLSpecMeas::cleanup_after_load( const unsigned int flags )
{
  SpecUtils::SpecFile::cleanup_after_load( flags );

  //should detect if the detector was loaded, and if not, if we know the type,
  //  we could then load it.
  //   -instead for right now lets do this in InterSpec...
}//void QLSpecMeas::cleanup_after_load()


