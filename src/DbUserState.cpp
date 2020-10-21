#include "InterSpec_config.h"

#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>


#include <Wt/Dbo/Dbo>
#include <Wt/WSpinBox>
#include <Wt/WCheckBox>
#include <Wt/WApplication>
#include <Wt/WRadioButton>
#include <Wt/WDoubleSpinBox>

#if( HAS_WTDBOSQLITE3 )
#include "Wt/Dbo/backend/Sqlite3"
#endif

#if( HAS_WTDBOMYSQL )
#include "Wt/Dbo/backend/MySQL"
#endif

#if( HAS_WTDBOPOSTGRES )
#include "Wt/Dbo/backend/Postgres"
#endif

#if( HAS_WTDBOFIREBIRD )
#include "Wt/Dbo/backend/Firebird"
#endif

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DbUserState.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "external_libs/diff-match-patch-cpp-stl/diff_match_patch.h"

using namespace Wt;
using namespace std;

namespace Wt
{
  namespace Dbo
  {
    const char *sql_value_traits<DbUserState::SpectrumData_t>
    ::type( SqlConnection *conn, int size)
    {
#if( HAS_WTDBOSQLITE3 )
      if( dynamic_cast<Wt::Dbo::backend::Sqlite3 *>(conn) )
        return "blob not null";
#endif
      
#if( HAS_WTDBOMYSQL )
      if( dynamic_cast<Wt::Dbo::backend::MySQL *>(conn) )
        return "MEDIUMBLOB";
#endif
      
#if( HAS_WTDBOPOSTGRES )
      if( dynamic_cast<Wt::Dbo::backend::Postgres *>(conn) )
        return "bytea not null";
#endif
      
#if( HAS_WTDBOFIREBIRD )
      if( dynamic_cast<Wt::Dbo::backend::Firebird *>(conn) )
        return "blob";
#endif
      cerr << "\n\n\nDbUserState:\n\tWarning, urogognized DB type\n";
      return conn->blobType();
    }
    
    void sql_value_traits< DbUserState::SpectrumData_t >
    ::bind( const DbUserState::SpectrumData_t &v, SqlStatement *statement,
            int column, int size)
    {
      statement->bind(column, v);
    }
    
    bool sql_value_traits<DbUserState::SpectrumData_t >
    ::read( DbUserState::SpectrumData_t &v, SqlStatement *statement,
            int column, int size )
    {
      return statement->getResult(column, &v, size);
    }
  }//namespace Dbo
}//namespace Wt


namespace DbUserState
{
  
Spectrum::Spectrum()
{
//  Wt::Dbo::ptr<InterSpecUser> m_user;
//  std::string m_uuid;
//  std::string m_filename;
//  std::string m_description;
  
  m_userHasModified = false;
//  Wt::WDateTime m_uploadTime;
//  Wt::WDateTime m_serializeTime;
//  std::string m_sessionID;
  
  m_numSamples = -1;
  m_isPassthrough = false;
  m_totalLiveTime = 0.0;
  m_totalRealTime = 0.0;
  m_totalGammaCounts = 0.0;
  m_totalNeutronCounts = 0.0;
  m_numDetectors = -1;
  m_hasNeutronDetector = false;
//  m_measurementsStartTime;
  
  m_writeprotected = false;
  m_isPartOfSaveState = false;
}//Spectrum()
  
  
  
void Spectrum::update( const SpecMeas &spec )
{
  diff_match_patch<string> dmp;
  dmp.Diff_Timeout = 0.01;
  
  
  SpectrumFile specfile;
  UserSpectrumStuff specstuff;
  specfile.setInformation( spec, DbUserState::SpectrumFile::SerializedFileFormat::k2011N42 );
  specstuff.setInformation( spec );
  
  Dbo::ptr<SpectrumFile> spectrumFile = m_spectrumFile.query();
  Dbo::ptr<UserSpectrumStuff> spectrumStuff = m_userSpectrumStuff.query();
  
  
  string str1 = "First string in dif that I dont care about";
  string str2 = "Second string in diff thatI care";
  
  std::list<diff_match_patch<string>::Patch> patches = dmp.patch_make(str1, str2);
  
  string strPatch = dmp.patch_toText( patches );
  
  pair<string, vector<bool> > out
  = dmp.patch_apply(dmp.patch_fromText(strPatch), str1);
  string strResult = out.first;

  
  
  throw runtime_error( "Spectrum::update not implementes" );
  

  
}//void update( std::shared_ptr<SpecMeas> spec )
  
  
std::shared_ptr<SpecMeas> Spectrum::undo()
{
  assert( 0 );
}//void undo()
  
  
std::shared_ptr<SpecMeas> Spectrum::redo()
{
  assert( 0 );
}//void redo()

  
std::shared_ptr<SpecMeas> Spectrum::assemble( Wt::Dbo::ptr<SpectrumFile> spectrumFile,
                                                 Wt::Dbo::ptr<UserSpectrumStuff> spectrumStuff )
{
  if( !spectrumFile || !spectrumStuff )
    throw runtime_error( "Spectrum::assemble: invalid input" );
  
  std::shared_ptr<SpecMeas> answer( new SpecMeas() );
  
  const SpectrumData_t &specdata = spectrumFile->m_spectrumData;
  const SpectrumData_t &peakdata = spectrumStuff->m_peaksXml;
  
  if( specdata.empty() )
    throw runtime_error( "Spucetrum data from DB SpectrumFIle is empty" );
  
  const bool loadedSpec = answer->SpecFile::load_N42_from_data( (char *)&specdata[0] );
  
  if( !loadedSpec )
    throw runtime_error( "Count load MeasurmentInformation from database data" );
  
  if( peakdata.size() )
  {
    ::rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_trim_whitespace; //rapidxml::parse_normalize_whitespace
    doc.parse<flags>( (char *)&peakdata[0] );
    const ::rapidxml::xml_node<char> *doc_node = doc.first_node();
    
    if( doc_node )
      answer->addPeaksFromXml( doc_node );
    set<set<int> > peakSamples = answer->sampleNumsWithPeaks();
    
    cerr << "Added peaks from XML: " << peakSamples.size() << endl;

#if( PERFORM_DEVELOPER_CHECKS )
    for( set<set<int> >::const_iterator i = peakSamples.begin(); i != peakSamples.end(); ++i )
    {
      std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > > peaks = answer->peaks( *i );
      
      if( !peaks )
      {
        log_developer_error( BOOST_CURRENT_FUNCTION, "Didnt get expected peaks after reading SpecMeas from database" );
        continue;
      }

//      for( const std::shared_ptr<const PeakDef> &peak : *peaks )
//      {
//      }//for( const std::shared_ptr<const PeakDef> &peak : *peaks )
    }//for( cosnt set<int> &i : peakSamples )
#endif
  }//if( peakdata.size() )
  
  //  int spectrumStuff->detectorID
  
  std::set<int> displayed_samples;
  
  const SpectrumData_t &sampleNumData = spectrumStuff->m_displayedSampleNumbers;
  if( sampleNumData.size() )
  {
    vector<int> samples;
    const char *horribleptr = (const char *)&sampleNumData[0];
    const bool ok = SpecUtils::split_to_ints( horribleptr, sampleNumData.size(), samples );
    if( !ok )
      throw runtime_error( "invalid displayed sample numbers" );
    
    displayed_samples.insert( samples.begin(), samples.end() );
  }//if( sampleNumData.size() )
  
  const vector<string> &detectors = answer->detector_names();
  
  answer->displayedSpectrumChangedCallback( spectrumStuff->m_spectrumType, answer,
                                            displayed_samples, detectors );
  
  return answer;
}//std::shared_ptr<SpecMeas> current()
  
bool Spectrum::isWriteProtected() const
{
  return m_writeprotected;
}//bool isWriteProtected() const

void Spectrum::makeWriteProtected( Wt::Dbo::ptr<Spectrum> ptr )
{
  if( !ptr || ptr->m_writeprotected )
    return;
  ptr.modify()->m_writeprotected = true;
}//void makeWriteProtected(...)
  
  
void Spectrum::removeWriteProtection( Wt::Dbo::ptr<Spectrum> ptr )
{
  if( !ptr || !ptr->m_writeprotected )
    return;
  ptr.modify()->m_writeprotected = false;
}//void removeWriteProtection(...)

  
SpectrumFile::SpectrumFile()
{
    
}//SpectrumFile()
  
  
void SpectrumFile::setInformation( const SpecUtils::SpecFile &spectrumFile,
                                const SerializedFileFormat format )
{
  m_fileFormat = format;
  m_compression = NoCompression;
  
  std::shared_ptr< ::rapidxml::xml_document<char> > xml = spectrumFile.SpecFile::create_2012_N42_xml();
  
  std::stringstream data;
  data << (*xml) << '\0';
  //  rapidxml::print( OutIt out, *xml );

  const string datastr = data.str();
  const unsigned char *horribleptr = (const unsigned char *)datastr.c_str();
  m_spectrumData.clear();
  m_spectrumData.insert( m_spectrumData.end(), horribleptr, horribleptr+datastr.size() );
}//void setFileData(...)
  
  
void UserSpectrumStuff::setInformation( const SpecMeas &spec )
{
//  detectorID
  
  {
    rapidxml::xml_document<char> doc;
    spec.addPeaksToXml( &doc );
    
    stringstream datastrm;
    datastrm << doc << '\0';
  
    const string datastr = datastrm.str();
    const unsigned char *horribleptr = (const unsigned char *)datastr.c_str();
    m_peaksXml.clear();
    m_peaksXml.insert( m_peaksXml.end(), horribleptr, horribleptr+datastr.size() );
  }
  
  
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *samples_node = spec.appendSampleNumbersToXml( &doc );
    
    stringstream datastrm;
    datastrm << (*samples_node) << '\0';
    
    const string datastr = datastrm.str();
    const unsigned char *horribleptr = (const unsigned char *)datastr.c_str();
    m_displayedSampleNumbers.clear();
    m_displayedSampleNumbers.insert( m_displayedSampleNumbers.end(), horribleptr, horribleptr+datastr.size() );
  }
  
  m_spectrumType = spec.displayType();

}//void setInformation( const SpecMeas &meas )

  
void SpectrumFile::setFileData( const std::string &path,
                                const SerializedFileFormat format )
{
    assert( 0 );
}//void setFileData(...)
  

void SpectrumFile::decodeSpectrum( std::shared_ptr<SpecUtils::SpecFile> &meas ) const
{
  if( !meas )
    throw runtime_error( "SpectrumFile::decodeSpectrum(): must have valid input pointer" );
  
  assert( 0 );
  
}//void decodeSpectrum(...)
  
}//namespace DbUserState





