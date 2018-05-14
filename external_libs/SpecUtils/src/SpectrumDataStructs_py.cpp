/**
 SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
 Copyright (C) 2016 William Johnson
 
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


#include "SpecUtils_config.h"

#include "SpecUtils/UtilityFunctions.h"
#include "SpecUtils/SpectrumDataStructs.h"


#include <set>
#include <iosfwd>
#include <vector>
#include <string>
#include <locale>
#include <limits>
#include <memory>
#include <iostream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <boost/python.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/concepts.hpp>


#include <datetime.h>  // compile with -I/path/to/python/include


using namespace std;


namespace
{
  //Begin structs to help convert between c++ and python
  
  class PythonInputDevice
  {
  public:
    typedef char char_type;
   
    // XXX - not working when tellg needs to be called
    
    //  struct category : boost::iostreams::source_tag, boost::iostreams::seekable_device_tag
    //  {};
    //  struct category : boost::iostreams::source_tag {};
    typedef boost::iostreams::seekable_device_tag category;
    
    explicit PythonInputDevice( boost::python::object object ) : object_(object)
    {}
    
    std::streamsize write( const char *s, std::streamsize n )
    {
      return 0;
    }
    
    std::streamsize read( char_type *buffer, std::streamsize buffer_size )
    {
       boost::python::object pyread = object_.attr( "read" );
      if( pyread.is_none() )
        throw std::runtime_error( "Python stream has no attibute 'read'" );
      
      boost::python::object py_data = pyread( buffer_size );
      const std::string data = boost::python::extract<std::string>( py_data );
      
      if( data.empty() && buffer_size != 0 )
        return -1;
      
      std::copy( data.begin(), data.end(), buffer );
      
      return data.size();
    }//read()
    
    
    boost::iostreams::stream_offset seek( boost::iostreams::stream_offset offset,
                                         std::ios_base::seekdir way )
    {
      int pway = 0;
      switch( way )
      {
        case std::ios_base::beg: pway = 0; break;
        case std::ios_base::cur: pway = 1; break;
        case std::ios_base::end: pway = 2; break;
      }//switch( way )
      
      boost::python::object pyseek = object_.attr( "seek" );
      boost::python::object pytell = object_.attr( "tell" );
      
      if( pyseek.is_none() )
        throw std::runtime_error( "Python stream has no attibute 'seek'" );
      
      if( pytell.is_none() )
        throw std::runtime_error( "Python stream has no attibute 'tell'" );
      
      const boost::python::object pynewpos = pyseek( offset, pway );
      const boost::python::object pyoffset = pytell();
      const boost::python::extract<std::streamoff> newpos( pyoffset );
    
      return newpos;
    }//seek()
    
  private:
    boost::python::object object_;
  };//class PythonInputDevice
  
  
  class PythonOutputDevice
  {
    //see http://stackoverflow.com/questions/26033781/converting-python-io-object-to-stdostream-when-using-boostpython
  public:
    typedef char char_type;
    
    struct category : boost::iostreams::sink_tag, boost::iostreams::flushable_tag
    {};
    
    explicit PythonOutputDevice( boost::python::object object ) : object_( object )
    {}
    
    std::streamsize write( const char *buffer, std::streamsize buffer_size )
    {
      boost::python::str data(buffer, buffer_size);
      boost::python::object pywrite = object_.attr( "write" );
      
      if( pywrite.is_none() )
        throw std::runtime_error( "Python stream has no attibute 'write'" );
      
      boost::python::object pynwrote = pywrite( data );
      const boost::python::extract<std::streamsize> bytes_written( pynwrote );
      return bytes_written.check() ? bytes_written : buffer_size;
    }
    
    bool flush()
    {
      boost::python::object flush = object_.attr( "flush" );
      
      if( flush.is_none() )
        throw std::runtime_error( "Python stream has no attibute 'flush'" );
        
      flush();
      return true;
    }
    
  private:
    boost::python::object object_;
  };//class PythonOutputDevice
  
  
  //wcjohns got the datetime convertion code 20150515 from
  //  http://en.sharejs.com/python/13125
  //There is also a time duration convertion code as wel, but removed it.
  /**
   * Convert boost::posix_ptime objects (ptime and time_duration)
   * to/from python datetime objects (datetime and timedelta).
   *
   * Credits:
   * http://libtorrent.svn.sourceforge.net/viewvc/libtorrent/trunk/bindings/python/src/datetime.cpp
   * http://www.nabble.com/boost::posix_time::ptime-conversion-td16857866.html
   */
  
  static long get_usecs(boost::posix_time::time_duration const& d)
  {
    static long resolution
    = boost::posix_time::time_duration::ticks_per_second();
    long fracsecs = d.fractional_seconds();
    if (resolution > 1000000)
      return fracsecs / (resolution / 1000000);
    else
      return fracsecs * (1000000 / resolution);
  }
  
  
  /* Convert ptime to/from python */
  struct ptime_to_python_datetime
  {
    static PyObject* convert(boost::posix_time::ptime const& pt)
    {
      boost::gregorian::date date = pt.date();
      boost::posix_time::time_duration td = pt.time_of_day();
      return PyDateTime_FromDateAndTime((int)date.year(),
                                        (int)date.month(),
                                        (int)date.day(),
                                        td.hours(),
                                        td.minutes(),
                                        td.seconds(),
                                        get_usecs(td));
    }
  };
  
  
  struct ptime_from_python_datetime
  {
    ptime_from_python_datetime()
    {
      boost::python::converter::registry::push_back(
                                                    &convertible,
                                                    &construct,
                                                    boost::python::type_id<boost::posix_time::ptime > ());
    }
    
    static void* convertible(PyObject * obj_ptr)
    {
      if ( ! PyDateTime_Check(obj_ptr))
        return 0;
      return obj_ptr;
    }
    
    static void construct(
                          PyObject* obj_ptr,
                          boost::python::converter::rvalue_from_python_stage1_data * data)
    {
      PyDateTime_DateTime const* pydate
      = reinterpret_cast<PyDateTime_DateTime*>(obj_ptr);
      
      // Create date object
      boost::gregorian::date _date(PyDateTime_GET_YEAR(pydate),
                                   PyDateTime_GET_MONTH(pydate),
                                   PyDateTime_GET_DAY(pydate));
      
      // Create time duration object
      boost::posix_time::time_duration
      _duration(PyDateTime_DATE_GET_HOUR(pydate),
                PyDateTime_DATE_GET_MINUTE(pydate),
                PyDateTime_DATE_GET_SECOND(pydate),
                0);
      // Set the usecs value
      _duration += boost::posix_time::microseconds(PyDateTime_DATE_GET_MICROSECOND(pydate));
      
      // Create posix time object
      void* storage = (
                       (boost::python::converter::rvalue_from_python_storage<boost::posix_time::ptime>*)
                       data)->storage.bytes;
      new (storage)
      boost::posix_time::ptime(_date, _duration);
      data->convertible = storage;
    }
  };
  
  
  
  
  //Begin wrapper functions
  void loadFile( MeasurementInfo *info,
                 const std::string &filename,
                 ParserType parser_type,
                 std::string file_ending_hint = "" )
  {
    const bool success = info->load_file( filename, parser_type, file_ending_hint );
    if( !success )
    {
      if( parser_type == kAutoParser )
        throw std::runtime_error( "Couldnt parse file " + filename );
      
      string type;
      switch( parser_type )
      {
        case k2006Icd1Parser: type = "2006 N42"; break;
        case K2011ICD1Parser: type = "2011 N42"; break;
        case kSpcParser: type = "SPC"; break;
        case kGR135Parser: type = "GR135"; break;
        case kPcfParser: type = "PCF"; break;
        case kChnParser: type = "CHN"; break;
        case kIaeaParser: type = "IAEA"; break;
        case kTxtOrCsvParser: type = "TXT or CSV"; break;
        case kCanberraCnfParser: type = "CNF"; break;
        case kTracsMpsParser: type = "MPS"; break;
        case kSPMDailyFile: type = "SpectroscopicPortalMonitor"; break;
        case kAmptekMca: type = "Amptek MCA"; break;
        case kMicroRaider: type = "Micro Raider"; break;
        case kAramParser: type = "Aram"; break;
        case kOrtecListMode: type = "Ortec Listmode"; break;
        case kAutoParser: type = ""; break;
      }//switch( parser_type )
      
      throw std::runtime_error( filename + " couldnt be parsed as a " + type + " file." );
    }//if( !success )
  }//loadFile(...)
  
  
  //I couldnt quite figure out how to get Python to play nicely with const
  //  references to smart pointers, so instead am using some thin wrapper
  //  functions to return a smart pointer by value.
  std::shared_ptr< const std::vector<float> > channel_energies_wrapper( Measurement *meas )
  {
    return meas->channel_energies();
  }
  
  std::shared_ptr< const std::vector<float> > gamma_counts_wrapper( Measurement *meas )
  {
    return meas->gamma_counts();
  }
  
  boost::posix_time::ptime start_time_wrapper( Measurement *meas )
  {
    return meas->start_time();
  }
  
  boost::python::list get_measurments_wrapper( MeasurementInfo *info )
  {
    //This function overcomes an issue where returning a vector of
    //  std::shared_ptr<const Measurment> objects to python, and then in python
    //  derefrencing an element and using it causes an issue in python to say
    //  that std::shared_ptr<const Measurment> is an unknown class.
    boost::python::list l;
    for( auto p : info->measurements() )
    l.append( p );
    return l;
  }
  
  boost::python::list Measurment_remarks_wrapper( Measurement *info )
  {
    boost::python::list l;
    for( const string &p : info->remarks() )
    l.append( p );
    return l;
  }
  
  
  boost::python::list MeasurmentInfo_remarks_wrapper( MeasurementInfo *info )
  {
    boost::python::list l;
    for( const string &p : info->remarks() )
    l.append( p );
    return l;
  }
  
  boost::python::list gamma_channel_counts_wrapper( MeasurementInfo *info )
  {
    boost::python::list l;
    for( size_t p : info->gamma_channel_counts() )
    l.append( p );
    return l;
  }
  
  boost::python::list sample_numbers_wrapper( MeasurementInfo *info )
  {
    boost::python::list l;
    for( size_t p : info->sample_numbers() )
    l.append( p );
    return l;
  }
  
  boost::python::list detector_numbers_wrapper( MeasurementInfo *info )
  {
    boost::python::list l;
    for( size_t p : info->detector_numbers() )
    l.append( p );
    return l;
  }
  
  boost::python::list neutron_detector_names_wrapper( MeasurementInfo *info )
  {
    boost::python::list l;
    for( const string &p : info->neutron_detector_names() )
      l.append( p );
    return l;
  }
  
  
  boost::python::list detector_names_wrapper( MeasurementInfo *info )
  {
    boost::python::list l;
    for( const string &p : info->detector_names() )
    l.append( p );
    return l;
  }
  
  
  std::shared_ptr<Measurement> sum_measurements_wrapper( MeasurementInfo *info,
                                                          boost::python::list py_samplenums,
                                                          boost::python::list py_detnums )
  {
    set<int> samplenums;
    vector<int> detnums;
    
    boost::python::ssize_t n = boost::python::len( py_samplenums );
    for( boost::python::ssize_t i = 0; i < n; ++i )
      samplenums.insert( boost::python::extract<int>( py_samplenums[i] ) );
    
    n = boost::python::len( py_detnums );
    for( boost::python::ssize_t i = 0; i < n; ++i )
      detnums.push_back( boost::python::extract<int>( py_detnums[i] ) );
    
    return info->sum_measurements(samplenums, detnums );
  }//sum_measurements_wrapper(...)
  
  
  
  void writePcf_wrapper( MeasurementInfo *info, boost::python::object pystream )
  {
    boost::iostreams::stream<PythonOutputDevice> output( pystream );
    if( !info->write_pcf( output ) )
      throw std::runtime_error( "Failed to write PCF file." );
  }
  
  void write2006N42_wrapper( MeasurementInfo *info, boost::python::object pystream )
  {
    boost::iostreams::stream<PythonOutputDevice> output( pystream );
    if( !info->write_2006_N42( output ) )
      throw std::runtime_error( "Failed to write 2006 N42 file." );
  }
  
  void write2011N42Xml_wrapper( MeasurementInfo *info, boost::python::object pystream )
  {
    boost::iostreams::stream<PythonOutputDevice> output( pystream );
    if( !info->write_2011_N42( output ) )
      throw std::runtime_error( "Failed to write 2011 N42 file." );
  }
  
  void writeCsv_wrapper( MeasurementInfo *info, boost::python::object pystream )
  {
    boost::iostreams::stream<PythonOutputDevice> output( pystream );
    if( !info->write_csv( output ) )
      throw std::runtime_error( "Failed to write CSV file." );
  }
  
  void writeTxt_wrapper( MeasurementInfo *info, boost::python::object pystream )
  {
    boost::iostreams::stream<PythonOutputDevice> output( pystream );
    if( !info->write_txt( output ) )
      throw std::runtime_error( "Failed to TXT file." );
  }
  
  void writeIntegerChn_wrapper( MeasurementInfo *info,
                               boost::python::object pystream,
                               boost::python::object py_sample_nums,
                               boost::python::object py_det_nums )
  {
    std::set<int> sample_nums, det_nums;
    
    boost::python::list sn_list = boost::python::extract<boost::python::list>(py_sample_nums);
    boost::python::list dn_list = boost::python::extract<boost::python::list>(py_det_nums);
    
    boost::python::ssize_t n = boost::python::len( sn_list );
    for( boost::python::ssize_t i = 0; i < n; ++i )
      sample_nums.insert( boost::python::extract<int>( sn_list[i] ) );
    
    n = boost::python::len( dn_list );
    for( boost::python::ssize_t i = 0; i < n; ++i )
      det_nums.insert( boost::python::extract<int>( dn_list[i] ) );
    
    boost::iostreams::stream<PythonOutputDevice> output( pystream );
    if( !info->write_integer_chn( output, sample_nums, det_nums ) )
      throw std::runtime_error( "Failed to write Integer CHN file." );
  }
  
  
  void setInfoFromN42File_wrapper( MeasurementInfo *info,
                                  boost::python::object pystream )
  {
    boost::iostreams::stream<PythonInputDevice> input( pystream );
    if( !info->load_from_N42( input ) )
      throw std::runtime_error( "Failed to decode input as a valid N42 file." );
  }//setInfoFromN42File_wrapper(...)
  
  
  void setInfoFromPcfFile_wrapper( MeasurementInfo *info,
                                  boost::python::object pystream )
  {
    boost::iostreams::stream<PythonInputDevice> input( pystream );
    if( !info->load_from_pcf( input ) )
      throw std::runtime_error( "Failed to decode input as a valid PCF file." );
  }//setInfoFromPcfFile_wrapper(...)
  
}//namespace



//Begin defintion of python functions

//make is so loadFile will take at least 3 arguments, but can take 4 (e.g. allow
//  default last argument)
BOOST_PYTHON_FUNCTION_OVERLOADS(loadFile_overloads, loadFile, 3, 4)

BOOST_PYTHON_MODULE(SpecUtils)
{
  using namespace boost::python;
  
  //enable python signature and user document string, but disable c++ signature
  //  for all docstrings created in the scope of local_docstring_options.
  docstring_options local_docstring_options( true, true, false );
  
  //Register our enums
  enum_<ParserType>( "ParserType" )
  .value( "k2006Icd1Parser", k2006Icd1Parser )
  .value( "K2011ICD1Parser", K2011ICD1Parser )
  .value( "kSpcParser", kSpcParser )
  .value( "kGR135Parser", kGR135Parser )
  .value( "kPcfParser", kPcfParser )
  .value( "kChnParser", kChnParser )
  .value( "kIaeaParser", kIaeaParser )
  .value( "kTxtOrCsvParser", kTxtOrCsvParser )
  .value( "kCanberraCnfParser", kCanberraCnfParser )
  .value( "kTracsMpsParser", kTracsMpsParser )
  .value( "kSPMDailyFile", kSPMDailyFile )
  .value( "kAmptekMca", kAmptekMca )
  .value( "kMicroRaider", kMicroRaider )
  .value( "kAutoParser", kAutoParser );
  
  
  enum_<DetectorType>( "DetectorType" )
  .value( "kGR135Detector", kGR135Detector )
  .value( "kIdentiFinderDetector", kIdentiFinderDetector )
  .value( "kIdentiFinderNGDetector", kIdentiFinderNGDetector )
  .value( "kIdentiFinderLaBr3Detector", kIdentiFinderLaBr3Detector )
  .value( "kDetectiveDetector", kDetectiveDetector )
  .value( "kDetectiveExDetector", kDetectiveExDetector )
  .value( "kDetectiveEx100Detector", kDetectiveEx100Detector )
  .value( "kOrtecIDMPortalDetector", kOrtecIDMPortalDetector )
  .value( "kSAIC8Detector", kSAIC8Detector )
  .value( "kFalcon5000", kFalcon5000 )
  .value( "kUnknownDetector", kUnknownDetector )
  .value( "kMicroDetectiveDetector", kMicroDetectiveDetector )
  .value( "kMicroRaiderDetector", kMicroRaiderDetector );
  
  
  {//begin Measurement class scope
    boost::python::scope Measurement_scope = class_<Measurement, boost::noncopyable>( "Measurement" )
    .def( "liveTime", &Measurement::live_time, "The live time help" )
    .def( "realTime", &Measurement::real_time )
    .def( "containedNeutron", &Measurement::contained_neutron )
    .def( "sampleNumber", &Measurement::sample_number )
    .def( "title", &Measurement::title, return_value_policy<copy_const_reference>() )
    .def( "occupied", &Measurement::occupied )
    .def( "gammaCountSum", &Measurement::gamma_count_sum )
    .def( "neutronCountsSum", &Measurement::neutron_counts_sum )
    .def( "speed", &Measurement::speed )
    .def( "latitude", &Measurement::latitude )
    .def( "longitude", &Measurement::longitude )
    .def( "hasGpsInfo", &Measurement::has_gps_info )
    //      .def( "positionTime", &Measurement::position_time, return_internal_reference<>() )
    .def( "detectorName", &Measurement::detector_name, return_value_policy<copy_const_reference>() )
    .def( "detectorNumber", &Measurement::detector_number )
    .def( "detectorType", &Measurement::detector_type, return_value_policy<copy_const_reference>() )
    .def( "qualityStatus", &Measurement::quality_status )
    .def( "sourceType", &Measurement::source_type )
    .def( "energyCalibrationModel", &Measurement::energy_calibration_model )
    //      .def( "remarks", &Measurement::remarks, return_internal_reference<>() )
    .def( "remarks", &Measurment_remarks_wrapper )
    //    .def( "startTime", &Measurement::start_time, return_internal_reference<>() )
    .def( "startTime", &start_time_wrapper )
    .def( "calibrationCoeffs", &Measurement::calibration_coeffs, return_internal_reference<>() )
    .def( "deviationPairs", &Measurement::deviation_pairs, return_internal_reference<>() )
    .def( "channelEnergies", &channel_energies_wrapper )
    .def( "gammaCounts", &gamma_counts_wrapper )
    .def( "neutronCounts", &Measurement::neutron_counts, return_internal_reference<>() )
    .def( "numGammaChannels", &Measurement::num_gamma_channels )
    .def( "findGammaChannel", &Measurement::find_gamma_channel )
    .def( "gammaChannelContent", &Measurement::gamma_channel_content )
    .def( "gammaChannelLower", &Measurement::gamma_channel_lower )
    .def( "gammaChannelCenter", &Measurement::gamma_channel_center )
    .def( "gammaChannelUpper", &Measurement::gamma_channel_upper )
    .def( "gammaChannelWidth", &Measurement::gamma_channel_width )
    .def( "gammaIntegral", &Measurement::gamma_integral )
    .def( "gammaChannelsSum", &Measurement::gamma_channels_sum )
    .def( "channelChannelEnergies", &channel_energies_wrapper )
    .def( "gammaChannelCounts", &gamma_counts_wrapper )
    .def( "gammaEnergyMin", &Measurement::gamma_energy_min )
    .def( "gammaEnergyMax", &Measurement::gamma_energy_max )
    //... setter functions here
    ;
    
    
    enum_<Measurement::SourceType>( "SourceType" )
    .value( "Background", Measurement::Background )
    .value( "Calibration", Measurement::Calibration )
    .value( "Foreground", Measurement::Foreground )
    .value( "IntrinsicActivity", Measurement::IntrinsicActivity )
    .value( "UnknownSourceType", Measurement::UnknownSourceType )
    .export_values();
    
    enum_<Measurement::EquationType>( "EquationType" )
    .value( "Polynomial", Measurement::Polynomial )
    .value( "FullRangeFraction", Measurement::FullRangeFraction )
    .value( "LowerChannelEdge", Measurement::LowerChannelEdge )
    .value( "UnknownEquationType", Measurement::UnknownEquationType )
    .export_values();
    
    enum_<Measurement::QualityStatus>( "QualityStatus" )
    .value( "Good", Measurement::Good )
    .value( "Suspect", Measurement::Suspect )
    .value( "Bad", Measurement::Bad )
    .value( "Missing", Measurement::Missing )
    .export_values();
    
    enum_<Measurement::OccupancyStatus>( "OccupancyStatus" )
    .value( "NotOccupied", Measurement::NotOccupied )
    .value( "Occupied", Measurement::Occupied )
    .value( "UnknownOccupancyStatus", Measurement::UnknownOccupancyStatus )
    .export_values();
  }//end Measurement class scope
  
  
  //Register smart pointers we will use with python.
  register_ptr_to_python< std::shared_ptr<Measurement> >();
  register_ptr_to_python< std::shared_ptr<const Measurement> >();
  register_ptr_to_python< std::shared_ptr< const std::vector<float> > >();
  
  implicitly_convertible< std::shared_ptr<Measurement>, std::shared_ptr<const Measurement> >();
  
  
  //Register vectors of C++ types we use to python
  class_< std::vector<float> >("FloatVec")
  .def( vector_indexing_suite<std::vector<float> >() );
  
  class_< std::vector<std::string> >("StringVec")
  .def( vector_indexing_suite<std::vector<std::string> >() );
  
  class_< std::vector<std::shared_ptr<const Measurement> > >("MeasurmentVec")
  .def( vector_indexing_suite<std::vector<std::shared_ptr<const Measurement> > >() );
  
  
  PyDateTime_IMPORT;
  ptime_from_python_datetime();
  boost::python::to_python_converter<const boost::posix_time::ptime, ptime_to_python_datetime>();
  
  //disambiguate a few functions that have overloads
  std::shared_ptr<const Measurement> (MeasurementInfo::*meas_fcn_ptr)(size_t) const = &MeasurementInfo::measurement;
  
  class_<MeasurementInfo>("MeasurementInfo")
  .def( "loadFile", &loadFile, loadFile_overloads(
        args( "file_name", "parser_type", "file_ending_hint" ),
        "Callling this function with parser_type==SpecUtils.ParserType.kAutoParser\n"
        "is the easiest way to load a spectrum file when you dont know the type of\n"
        "file.  The file_ending_hint is only used in the case of SpecUtils.ParserType.kAutoParser\n"
        "and uses the file ending to effect the order of parsers tried, example\n"
        "values for this might be: \"n24\", \"pcf\", \"chn\", etc. The entire filename\n"
        "can be passed in since only the letters after the last period are used.\n"
        "Throws RuntimeError if the file can not be opened or parsed." ) )
  .def( "modified", &MeasurementInfo::modified,
        "Indicates if object has been modified since last save." )
  .def( "numMeasurements", &MeasurementInfo::num_measurements,
        "Returns the number of measurments (sometimes called records) parsed." )
  .def( "measurement", meas_fcn_ptr, args("i"),
        "Returns the i'th measurment, where valid values are between 0 and\n"
        "MeasurementInfo.numMeasurements()-1.\n"
        "Throws RuntimeError if i is out of range." )
  .def( "measurements", &get_measurments_wrapper,
        "Returns a list of all Measurment's that were parsed." )
  .def( "gammaLiveTime", &MeasurementInfo::gamma_live_time,
        "Returns the sum of detector live times of the all the parsed Measurments." )
  .def( "gammaRealTime", &MeasurementInfo::gamma_real_time,
        "Returns the sum of detector real times (wall/clock time) of the all the\n"
        "parsed Measurments." )
  .def( "gammaCountSum", &MeasurementInfo::gamma_count_sum,
        "Returns the summed number of gamma counts from all parsed Measurments." )
  .def( "neutronCountsSum", &MeasurementInfo::neutron_counts_sum,
        "Returns the summed number of neutron counts from all parsed Measurments." )
  .def( "filename", &MeasurementInfo::filename, return_value_policy<copy_const_reference>(),
        "Returns the filename of parsed file; if the \"file\" was parsed from a\n"
        "stream, then may be empty unless user specifically set it using \n"
        "setFilename (not currently implemented for python)." )
  .def( "detectorNames", &detector_names_wrapper,
        "Returns a list of names for all detectors found within the parsed file.\n"
        "The list will be in the same order as (and correspond one-to-one with)\n"
        "the list MeasurementInfo.detectorNumbers() returns." )
  .def( "detectorNumbers", &detector_numbers_wrapper,
        "Returns a list of assigned detector numbers for all detectors found within\n"
        "the parsed file.  The list will be in the same order as (and correspond\n"
        "one-to-one with) the list MeasurementInfo.detectorNames() returns." )
  .def( "neutronDetectorNames", &neutron_detector_names_wrapper,
        "Returns list of names of detectors that contained neutron information." )
  .def( "uuid", &MeasurementInfo::uuid, return_value_policy<copy_const_reference>(),
        "Returns the unique ID string for this parsed spectrum file.  The UUID\n"
        "may have been specified in the input file itself, or if not, it is\n"
        "generated using the file contents.  This value will always be the same\n"
        "every time the file is parsed." )
  .def( "remarks", &MeasurmentInfo_remarks_wrapper,
        "Returns a list of remarks or comments found while parsing the spectrum file.\n"
        "May include parser generated warnings or notes." )
  .def( "laneNumber", &MeasurementInfo::lane_number,
        "Returns the lane number of the RPM if specified in the spectrum file, otherwise\n"
        "will have a value of -1." )
  .def( "measurementLocationName", &MeasurementInfo::measurement_location_name, return_value_policy<copy_const_reference>(),
        "Returns the location name specified in the spectrum file; will be an\n"
        "empty string if not specified." )
  .def( "inspection", &MeasurementInfo::inspection, return_value_policy<copy_const_reference>(),
        "Returns the inspection type (e.g. primary, secondary, etc.) specified\n"
        "in the spectrum file. If not specified an empty string will be returned." )
  .def( "measurmentOperator", &MeasurementInfo::measurment_operator, return_value_policy<copy_const_reference>(),
        "Returns the detector operators name if specified in the spectrum file.\n"
        "If not specified an empty string will be returned." )
  .def( "sampleNumbers", &sample_numbers_wrapper,
        "If a spectrum file contains multiple measurments (records) from multiple\n"
        "detectors, the the measurments for the same time intervals will be grouped\n"
        "into unique groupings of sample and detectors, with the sample number\n"
        "generally increasing for measurments taken later in time.\n"
        "This function returns a list of all sample numbers in the parsed file." )
  .def( "numMeasurements", &MeasurementInfo::num_measurements,
        "Returns the number of measurments (records) parsed from the spectrum file." )
  .def( "detectorType", &MeasurementInfo::detector_type,
        "Returns the detector type specified in the spectrum file, or an empty string\n"
        "if none was specified.  Example values could include: 'HPGe 50%' or 'NaI'.")
  .def( "instrumentType", &MeasurementInfo::instrument_type, return_value_policy<copy_const_reference>(),
        "Returns the instrument type if specified in (or infered from) the spectrum\n"
        "file, or an empty string otherwise. Example values could include: PortalMonitor,\n"
        "SpecPortal, RadionuclideIdentifier, etc." )
  .def( "manufacturer", &MeasurementInfo::manufacturer, return_value_policy<copy_const_reference>(),
        "Returns the detector manufacturer if specified (or infered), or an empty\n"
        "string otherwise." )
  .def( "instrumentModel", &MeasurementInfo::instrument_model, return_value_policy<copy_const_reference>(),
        "Returns the instrument model if specified, or infered from, the spectrum file.\n"
        "Returns empty string otherwise.  Examples include: 'Falcon 5000', 'ASP', \n"
        "'identiFINDER', etc." )
  .def( "instrumentId", &MeasurementInfo::instrument_id, return_value_policy<copy_const_reference>(),
        "Returns the instrument ID (typically the serial number) specified in the\n"
        "file, or an empty string otherwise." )
  //     inline std::shared_ptr<const DetectorAnalysis> detectors_analysis() const;
  .def( "hasGpsInfo", &MeasurementInfo::has_gps_info,
        "Returns True if any of the measurments contained valid GPS data." )
  .def( "meanLatitude", &MeasurementInfo::mean_latitude,
        "Returns the mean latitidue of all measurments with valid GPS data.  If no\n"
        "GPS data was availble, will return something close to -999.9." )
  .def( "meanLongitude", &MeasurementInfo::mean_longitude,
        "Returns the mean longitude of all measurments with valid GPS data.  If no\n"
        "GPS data was availble, will return something close to -999.9." )
  .def( "memmorysize", &MeasurementInfo::memmorysize,
        "Returns the approximate (lower bound) of bytes this object takes up in memory." )
  .def( "gammaChannelCounts", &gamma_channel_counts_wrapper,
        "Returns the set of number of channels the gamma data has. If all measurments\n"
        "in the file contained the same number of channels, then the resulting list\n"
        "will have one entry with the number of channels (so typically 1024 for Nai,\n"
        "16384 for HPGe, etc.).  If there are detectors with different numbers of bins,\n"
        "then the result returned will have multiple entries.")
  .def( "numGammaChannels", &MeasurementInfo::num_gamma_channels,
        "Returns the number of gamma channels of the first (gamma) detector found\n"
        " or 0 if there is no gamma data.")
  .def( "backgroundSampleNumber", &MeasurementInfo::background_sample_number,
        "Returns the first background sample number in the spectrum file, even if\n"
        "there is more than one background sample number." )
  .def( "reset", &MeasurementInfo::reset,
        "Resets the MeasurementInfo object to its initial (empty) state." )
  .def( "sumMeasurements", &sum_measurements_wrapper,
        args("SampleNumbers", "DetectorNumbers"),
        "Sums the measurments of the specified sample and detector numbers.\n"
        "SampleNumbers and DetectorNumbers are both lists of integers.\n"
        "If the measurments contain different energy binnings, one will be chosen\n"
        "and the other measurments rebinned before summing so that energies stay\n"
        "consistent (e.g. not just a bin-by-bin summing).\n"
        "Throws RuntimeError if SampleNumbers or DetectorNumbers contain invalid\n"
        "entries." )
  .def( "writePcf", &writePcf_wrapper, boost::python::arg("OutputStream"),
        "The PCF format is the binary native format of GADRAS.  Saving to this format\n"
        "will cause the loss of some information. However, Calibration,\n"
        "foreground/background, speed, sample, and spectrum title (up to 60 characters)\n"
        "will be preserved along with the spectral information and neutron counts.\n"
        "Throws RuntimeError on failure."
       )
  .def( "write2006N42", &write2006N42_wrapper, args("OutputStream"),
        "Writes a 2006 version of ICD1 N42 file to OutputStream; most information\n"
        "is preserved in the output.\n"
        "Throws RuntimeError on failure." )
  .def( "writeCsv", &writeCsv_wrapper, args("OutputStream"),
        "The spectra are written out in a two column format (seperated by a comma);\n"
        "the first column is gamma channel lower edge energy, the second column is\n"
        "channel counts.  Each spectrum in the file are written out contiguously and\n"
        "seperated by a header that reads \"Energy, Data\".  Windows style line endings\n"
        "are used (\\n\\r).  This format loses all non-spectral information, including\n"
        "live and real times, and is intended to be an easy way to import the spectral\n"
        "information into other programs like Excel.\n"
        "Throws RuntimeError on write failure." )
  .def( "writeTxt", &writeTxt_wrapper, args("OutputStream"),
        "Spectrum(s) will be written to an ascii text format.  At the beggining of the\n"
        "output the original file name, total live and real times, sum gamma counts,\n"
        "sum neutron counts, and any file level remarks will be written on seperate\n"
        "labeled lines. Then after two blank lines each spectrum in the current file\n"
        "will be written, seperated by two blank lines.  Each spectrum will contain\n"
        "all remarks, measurment start time (if valid), live and real times, sample\n"
        "number, detector name, detector type, GPS coordinates/time (if valid), \n"
        "serial number (if present), energy calibration type and coefficient values,\n"
        "and neutron counts (if valid); the channel number, channel lower energy,\n"
        "and channel counts is then provided with each channel being placed on a\n"
        "seperate line and each field being seperated by a space.\n"
        "Any detector provided analysis in the original program, as well\n"
        "manufacturer, UUID, deviation pairs, lane information, location name, or\n"
        "spectrum title is lost.\n"
        "Other programs may not be able to read back in all information written to\n"
        "the txt file.\n"
        "The Windows line ending convention is used (\\n\\r).\n"
        "This is not a standard format commonly read by other programs, and is\n"
        "intended as a easily human readable summary of the spectrum file information."
        "Throws RuntimeError on failure." )
  .def( "write2011N42Xml", &write2011N42Xml_wrapper, args("OutputStream"),
        "Saves to the 2011 N42 XML format.  Nearly all relevant"
        " information in most input spectrum files will also be saved in"
        " to the output stream."
        "Throws RuntimeError on failure." )
  .def( "writeIntegerChn", &writeIntegerChn_wrapper,
        args("OutputStream","SampleNumbers","DetectorNumbers"),
        "Writes an integer binary CHN file to OutputStream.  This format holds a\n"
        "single spectrum, so you must specify the sample and detector numbers you\n"
        "would like summed; if SampleNumbers or DetectorNumbers are empty, then all\n"
        "samples or detectors will be used.\n"
        "This format preserves the gamma spectrum, measurment start time, spectrum\n"
        "title (up to 63 characters), detector description, and energy calibration.\n"
        "Energy deviation pairs and neutron counts, as well as any other meta\n"
        "information is not preserved.\n"
        "SampleNumbers and DetectorNumbers are both lists of integers.\n"
        "If the measurments contain different energy binnings, one will be chosen\n"
        "and the other measurments rebinned before summing so that energies stay\n"
        "consistent (e.g. not just a bin-by-bin summing).\n"
        "Throws RuntimeError if SampleNumbers or DetectorNumbers contain invalid\n"
        "entries, or there is a error writing to OutputStream." )
  //    .def( "", &_wrapper )
  //probably write thin wrappers for the bellow
  //   enum SpcBinaryType{ IntegerSpcType, FloatSpcType };
  //   bool writeBinarySpcFile( std::ostream &ostr, const SpcBinaryType type, std::set<int> sample_nums, const std::set<int> &det_nums ) const;
  .def( "setInfoFromN42File", &setInfoFromN42File_wrapper, args("InputStream"),
        "Parses the InputStream as a N42 (2006, 2011 and HRPDS variants) spectrum\n"
        " file.\n"
        "Throws RuntimeError on parsing or data reading failure as well as reseting\n"
        "the input stream to its original position.\n"
        "InputStream must support random access seeking (one seek to end of the file\n"
        "is used to determine input size, then its reset to begining and read serially)." )
  .def( "setInfoFromPcfFile", &setInfoFromPcfFile_wrapper, args("InputStream"),
        "Parses the InputStream as a GADRAS PCF file."
        "InputStream must support random access seeking.\n"
        "Throws RuntimeError on parsing or data reading failure." )
  
  //... lots more functions here
  ;
  
}//BOOST_PYTHON_MODULE(SpecUtils)
