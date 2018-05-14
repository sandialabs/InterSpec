#ifndef DbUserState_h
#define DbUserState_h

#include "InterSpec_config.h"

#include <string>

#include <Wt/Dbo/Dbo>
#include <Wt/WString>
#include <Wt/WDateTime>

#include "InterSpec/SpecMeas.h"

class MeasurementInfo;

namespace DbUserState
{
  class Spectrum;
  class SpectrumFile;
  class UserSpectrumStuff;
  
  class DataStr : public std::vector<unsigned char>{ public: };
  typedef DataStr SpectrumData_t;
}

namespace Wt
{
  namespace Dbo
  {
    //We have to specialize sql_value_traits<SpectrumData_t, void> so that for MySQL
    //  not just a BLOB will be used, but a MEDIUMBLOB (a pain!)
    template<>
    struct sql_value_traits<DbUserState::SpectrumData_t, void>
    {
      static const bool specialized = true;
      static const char *type(SqlConnection *conn, int size);
      static void bind(const DbUserState::SpectrumData_t& v,
                       SqlStatement *statement, int column, int size);
      static bool read(DbUserState::SpectrumData_t& v, SqlStatement *statement,
                       int column, int size);
    };
  }//namespace Dbo
}//namespace Wt


namespace DbUserState
{
class Spectrum;
class SpectrumFile;
class UserSpectrumStuff;
struct ShieldingSourceFitModel;
  
//InterSpecUserTmp: a urogate class for InterSpecUser for development purposes
class InterSpecUserTmp
{
public:
  InterSpecUserTmp()
  {
    m_userName = "DevelopingUser";
  }
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::field( a, m_userName, "UserName" );
    Wt::Dbo::hasMany( a, m_spectrums, Wt::Dbo::ManyToOne, "InterSpecUserTmp" );
    Wt::Dbo::hasMany( a, m_shieldSrcModels, Wt::Dbo::ManyToOne, "InterSpecUserTmp");
  }//void persist( Action &a )
  
protected:
  std::string m_userName;
  Wt::Dbo::collection< Wt::Dbo::ptr<Spectrum> >         m_spectrums;
  Wt::Dbo::collection< Wt::Dbo::ptr<ShieldingSourceFitModel> > m_shieldSrcModels;
//  Wt::Dbo::collection< Wt::Dbo::ptr<UserState> >            m_userStates;

};//class InterSpecUserTmp
  

class SpectrumDiff
{
public:
  enum FieldDiffed
  {
    SpectrumFileField,
    PeaksField,
    //DetectorField
    SpectrumTypeField,
    DisplayedSampleNumbersField
  };//enum FieldDiffed
  
  FieldDiffed m_field;
  
  //m_diff: applying this to m_spectrum, will bring this
  std::string m_diff;
  
  Wt::Dbo::ptr<Spectrum> m_spectrum;
  
  
  //m_reduDiff: the diff
  Wt::Dbo::ptr<SpectrumDiff> m_reduDiff;
  Wt::Dbo::weak_ptr<SpectrumDiff> m_undoDiff;
  
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, m_spectrum, "SpectrumDiff", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::belongsTo( a, m_reduDiff, "UndoDiff" );
    Wt::Dbo::hasOne( a, m_undoDiff, "UndoDiff" );
    
    Wt::Dbo::field( a, m_field, "FieldDiffed" );
    Wt::Dbo::field( a, m_diff, "DiffTxt" );
  }
  
};//class SpectrumDiff


class Spectrum
{
public:
  //sm_maxFileSizeBytes: maximum size of spectrum file we'll save to the
  //  database.  Dictated by MySQL medium blob size (a choice in
  //  sql_value_traits<FileData_t>).
  const static size_t sm_maxFileSizeBytes = 16777215;
  const static int sm_maxUuidLength = 40;
  const static int sm_maxSessionIdLength = 32;
  
public:
  Spectrum();
  
  void update( const SpecMeas &spec );
  std::shared_ptr<SpecMeas> undo();
  std::shared_ptr<SpecMeas> redo();
  
  static std::shared_ptr<SpecMeas> assemble( Wt::Dbo::ptr<SpectrumFile> spectrumFile,
                                  Wt::Dbo::ptr<UserSpectrumStuff> spectrumStuff );
  
  
  Wt::Dbo::ptr<InterSpecUserTmp> m_user;
  
  std::string m_uuid;
  std::string m_filename;
  std::string m_description;
  
  bool m_userHasModified;
  Wt::WDateTime m_uploadTime;
  Wt::WDateTime m_serializeTime;
  std::string m_sessionID;
  
  int m_numSamples;
  bool m_isPassthrough;
  double m_totalLiveTime;
  double m_totalRealTime;
  double m_totalGammaCounts;
  double m_totalNeutronCounts;
  int m_numDetectors;
  bool m_hasNeutronDetector;
  Wt::WDateTime m_measurementsStartTime;
  

  Wt::Dbo::weak_ptr<SpectrumFile> m_spectrumFile;
  Wt::Dbo::weak_ptr<UserSpectrumStuff> m_userSpectrumStuff;

  Wt::Dbo::ptr<Spectrum> m_snapshotParent;
  Wt::Dbo::collection< Wt::Dbo::ptr<Spectrum> > m_snapshots;
  
  Wt::Dbo::collection< Wt::Dbo::ptr<ShieldingSourceFitModel> > m_modelsUsedWith;
  
  
  Wt::Dbo::collection< Wt::Dbo::ptr<SpectrumDiff> > m_diffs;

  
  //isWriteProtected(): returns if this UserFileInDb has been made write
  //  protected
  bool isWriteProtected() const;
  
  //makeWriteProtected(...): sets the writeprotected flag so object wont
  //  be re-saved to the database in an altered state.  Should only be called
  //  from within an active transaction
  static void makeWriteProtected( Wt::Dbo::ptr<Spectrum> ptr );
  
  //removeWriteProtection(...): sets the writeprotected flag so object can
  //  be re-saved to the database.  Should only be called from within an active
  //  transaction
  static void removeWriteProtection( Wt::Dbo::ptr<Spectrum> ptr );
  
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, m_user, "InterSpecUserTmp", Wt::Dbo::OnDeleteCascade );
    
    Wt::Dbo::belongsTo( a, m_snapshotParent, "SnapshotParent", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::hasMany( a, m_snapshots, Wt::Dbo::ManyToOne, "SnapshotParent" );
    
    Wt::Dbo::hasMany( a, m_diffs, Wt::Dbo::ManyToOne, "SpectrumDiff" );
    
    Wt::Dbo::hasOne( a, m_spectrumFile, "SpectrumFile" );
    Wt::Dbo::hasOne( a, m_userSpectrumStuff, "UserSpectrumStuff" );

    Wt::Dbo::hasMany( a, m_modelsUsedWith,
                      Wt::Dbo::ManyToMany, "shielding_source_fit_model" );

    Wt::Dbo::field( a, m_uploadTime,      "UploadTime" );
    Wt::Dbo::field( a, m_uuid,            "UUID", sm_maxUuidLength );
    Wt::Dbo::field( a, m_filename,        "Filename" );
    Wt::Dbo::field( a, m_description,     "Description" );
    Wt::Dbo::field( a, m_writeprotected,  "WriteProtected" );
    Wt::Dbo::field( a, m_userHasModified, "UserHasModified" );
    Wt::Dbo::field( a, m_sessionID,       "SessionID", sm_maxSessionIdLength );
    
    Wt::Dbo::field( a, m_numSamples,            "NumSamples" );
    Wt::Dbo::field( a, m_isPassthrough,         "IsPassthrough" );
    Wt::Dbo::field( a, m_totalLiveTime,         "TotalLiveTime" );
    Wt::Dbo::field( a, m_totalRealTime,         "TotalRealTime" );
    Wt::Dbo::field( a, m_totalGammaCounts,      "TotalGammaCounts" );
    Wt::Dbo::field( a, m_totalNeutronCounts,    "TotalNeutronCounts" );
    Wt::Dbo::field( a, m_numDetectors,          "NumDetectors" );
    Wt::Dbo::field( a, m_hasNeutronDetector,    "HasNeutronDetector" );
    Wt::Dbo::field( a, m_measurementsStartTime, "MeasurementsStartTime" );
    Wt::Dbo::field( a, m_isPartOfSaveState,     "IsPartOfSaveState" );
    
    
    if( a.getsValue() )
    {
      m_serializeTime = Wt::WDateTime::currentDateTime();
      Wt::Dbo::field( a, m_serializeTime, "SerializeTime" );
    }//if( saving to DB )
    
    if( a.setsValue() || a.isSchema() )
    {
      Wt::Dbo::field( a, m_serializeTime, "SerializeTime" );
    }//if( reading from DB )
  }//void persist( Action &a )
  
private:
  //writeprotected: not yet fully implemented, but intended to allow a way to
  //  ensure database entry doesnt get changed in the future.  This is mostly
  //  relevant for saving the applications state, you dont want the user to
  //  change the individual components, thus corrupting the state
  bool m_writeprotected;
  
  bool m_isPartOfSaveState;
};//class Spectrum


//This class stores a MeasurementInfo object
class SpectrumFile
{
public:
  enum CompressionType
  {
    NoCompression
  };//enum CompressionType
  
  enum SerializedFileFormat
  {
    k2011N42
  };//enum SerializedFileFormat
  
  const static CompressionType sm_defaultCompressionType;
  const static SerializedFileFormat sm_defaultSerializationFormat;

  
public:
  SpectrumFile();
  Wt::Dbo::ptr<Spectrum> m_spectrum;
  
  //m_gzipCompressed: only set true when writing if ALLOW_SAVE_TO_DB_COMPRESSION
  //  has been set.
  CompressionType m_compression;
  
  //So we can choose to save files to the datbase in something besides native
  //  binary format
  SerializedFileFormat m_fileFormat;
  
  //fileData: the actual data of the serialized MeasurementInfo object, may be
  //  compressed.
  SpectrumData_t m_spectrumData;
  
  //setFileData(...): serializes spectrumFile to fileData as a binary native
  //  file format.
  //  Will throw FileToLargeForDbException if serialization is larger than
  //  UserFileInDb::sm_maxFileSizeBytes.
  //  Will throw runtime_error if any other issues.
  void setInformation( const MeasurementInfo &spectrumFile,
                       const SerializedFileFormat format );
  
  //setFileData(...): same as other setFileData(...) function, but instead
  //  sets the data from a file on the filesystem.  Will throw if the file
  //  reading fails for any reason.
  void setFileData( const std::string &path,
                    const SerializedFileFormat format );
  
  //decodeSpectrum(): de-serializes data currently in fileData.
  //  Will throw if de-serialization fails, otherwise will always return
  //  a valid SpecMeas object.
  void decodeSpectrum( std::shared_ptr<MeasurementInfo> &meas ) const;
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, m_spectrum, "SpectrumFile", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::field( a, m_compression, "CompressionType" );
    Wt::Dbo::field( a, m_fileFormat, "FileFormat" );
    Wt::Dbo::field( a, m_spectrumData, "SpectrumData" );
  }//void persist( Action &a )
};//class UserFileInDbData
  
  
//UserSpectrumStuff nominally hold the stuff the SpecMeas class holds, minus
//  the stuff the MeasurementInfo class holds
class UserSpectrumStuff
{
public:
  
  Wt::Dbo::ptr<Spectrum> m_spectrum;
  
  void setInformation( const SpecMeas &meas );
  
  SpectrumData_t m_peaksXml;
//  int detectorID
  
  SpectrumType m_spectrumType;
  SpectrumData_t m_displayedSampleNumbers;
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, m_spectrum, "UserSpectrumStuff", Wt::Dbo::OnDeleteCascade );
    
    Wt::Dbo::field( a, m_peaksXml, "Peaks" );
    Wt::Dbo::field( a, m_displayedSampleNumbers, "DisplayedSampleNumbers" );
    
    if( a.getsValue() )
    {
      std::string type = descriptionText( m_spectrumType );
      Wt::Dbo::field( a, type, "SpectrumType" );
    }//if( saving to DB )
    
    if( a.setsValue() )
    {
      std::string type;
      Wt::Dbo::field( a, type, "SpectrumType" );
      
      if( type == descriptionText(kForeground) )
        m_spectrumType = kForeground;
      else if( type == descriptionText(kSecondForeground) )
        m_spectrumType = kSecondForeground;
      else if( type == descriptionText(kBackground) )
        m_spectrumType = kBackground;
      else if( type == "" )
        m_spectrumType = kForeground;
      else throw std::runtime_error( "Invalid spectrum type in DB: " + type );
    }//if( reading from DB )
    
    if( a.isSchema() )
    {
      std::string type;
      Wt::Dbo::field( a, type, "SpectrumType" );
    }
  }//void persist( Action &a )
  
};//class UserSpectrumStuff
  
  
struct ShieldingSourceFitModel
{
public:
  ShieldingSourceFitModel()
    : writeprotected(false)
  {}
    
  Wt::Dbo::ptr<InterSpecUserTmp> user;
  Wt::WString name;
  Wt::WString description;
    
  Wt::WDateTime creationTime;
  Wt::WDateTime serializeTime;
  
  Wt::Dbo::collection< Wt::Dbo::ptr<DbUserState::Spectrum> > spectrumsUsedWith;
    
    
  std::string xmlData;
    
  bool isWriteProtected() const;
    
  //makeWriteProtected(...): an active transaction must exist when calling this
  //  function
  static void makeWriteProtected( Wt::Dbo::ptr<ShieldingSourceFitModel> ptr );
    
  //removeWriteProtection(): an active transaction must exist when calling this
  //  function
  static void removeWriteProtection( Wt::Dbo::ptr<ShieldingSourceFitModel> ptr );
    
    
  void shallowEquals( const ShieldingSourceFitModel &rhs );
    
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, user, "InterSpecUserTmp", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::hasMany( a, spectrumsUsedWith,
                      Wt::Dbo::ManyToMany, "shielding_source_fit_model" );
    Wt::Dbo::field( a, name, "Name", 255 );
    Wt::Dbo::field( a, description, "Description", 511 );
    Wt::Dbo::field( a, writeprotected, "WriteProtected" );
      
    Wt::Dbo::field( a, creationTime, "CreationTime" );
    Wt::Dbo::field( a, serializeTime, "SerializeTime" );
    Wt::Dbo::field( a, xmlData, "XmlData" );
  }//void persist( Action &a )
    
  private:
  //writeprotected: not yet implemented, but intended to allow a way to ensure
  //  database entry doesnt get changed in the future.  This is mostly relevant
  //  for saving the applications state, you dont want the user to change the
  //  individual components, thus corrupting the state
  bool writeprotected;
};//struct ShieldingSourceFitModel
  
} //namespace DbUserState

#endif //DbUserState_h
