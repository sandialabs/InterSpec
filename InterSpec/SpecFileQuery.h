#ifndef SpecFileQuery_h
#define SpecFileQuery_h
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

#include <string>
#include <vector>
#include <memory>
#include <atomic>
#include <iostream>

#include <boost/any.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

struct SpecFileInfoToQuery;

/* Implimetation ideas
 -When user "hovers" over a row for abit, show a spectrum preview using the D3.js plotting (try to re-use logic from somewhere to what samples to plot).  Or similarish thing for when row is double-clicked.  (this all would be useful other places too...)
 -Help text! (could use this as an excuse to re-organize the help XML files).
 */

namespace SpecFileQuery
{
  //Stuff in this namespace is intended to be separate from the GUI and contain
  // most of the functional/logical code so it can eventually be used to form a
  // command line utility or something independent from the GUI.
  
  /** Enum representing which fields a user can search a spectrum file on.
   If you add an enum, you have to update both the C++ and the JavaScript
   */
  enum FileDataField
  {
    //String Fields to search
    ParentPath,
    Filename,
    DetectorName,
    SerialNumber,
    Manufacturer,
    Model,
    Uuid,
    Remark,
    LocationName,
    HasRIIDAnalysis,
    AnalysisResultText,
    AnalysisResultNuclide,
    
    //Discrete options
    DetectionSystemType,
    SearchMode,
    ContainedNuetronDetector,
    ContainedDeviationPairs,
    HasGps,
    EnergyCalibrationType,
    
    //More discrete types: ParserType,
    
    //Ranges of numeric values
    TotalLiveTime,
    TotalRealTime,
    IndividualSpectrumLiveTime,
    IndividualSpectrumRealTime,
    NumberOfSamples,
    NumberOfRecords,
    NumberOfGammaChannels,
    MaximumGammaEnergy,
    Latitude,
    Longitude,
    NeutronCountRate,
    GammaCountRate,
    
    //Time field
    StartTimeIoI,
    MeasurementsStartTimes,
    
    //More complex: peak location...
    
    NumFileDataFields
  };//enum FileDataField

  // Previous to 20231005 we printed live times as full duration strings, but to make data
  //  easier to work with in Excel, lets try just printing it out as seconds.
#define SpecFileQuery_TIME_AS_SECONDS 1
  
  const char *to_string( const FileDataField field );
  
  /** Returns true if field was succesfully found; if not found filed will be
      set to NumFileDataFields.  Does not throw.
   */
  bool from_string( const std::string &fieldName, FileDataField &field );
  
  
  enum TextFieldSearchType
  {
    TextIsExact,
    TextNotEqual,
    TextIsContained,
    TextDoesNotContain,
    TextStartsWith,
    TextDoesNotStartWith,
    TextEndsWith,
    TextDoesNotEndWith,
    TextRegex
  };//enum TextFieldSearchType
  
  const char *to_string( const TextFieldSearchType type );
  bool from_string( const std::string &val, TextFieldSearchType &type );
  
  enum NumericFieldMatchType
  {
    ValueIsExact,
    ValueIsNotEqual,
    ValueIsLessThan,
    ValueIsGreaterThan
  };//enum NumericFieldMatchType
  
  const char *to_string( const NumericFieldMatchType type );
  bool from_string( const std::string &val, NumericFieldMatchType &type );
  
  class SpecTest
  {
  public:
    SpecTest();
    
    //Will throw if field is not a string search filed
    void set_test( const FileDataField field, const std::string &searchstr, const TextFieldSearchType type );
    
    //Will throw if field is not a discrete field
    void set_discreet_test( const FileDataField field, const int val );
    
    //Will throw if field is not a numeric field
    void set_numeric_test( const FileDataField field, const double value, const NumericFieldMatchType type );
    
    //Will throw if not a time field
    void set_time_test( const FileDataField field, boost::posix_time::ptime comptime, const NumericFieldMatchType type );
    
    bool test( const SpecFileInfoToQuery &meas ) const;
    
    std::string summary() const;
    
    //Throw exception with explanation if not valid
    void isvalid();
    
    
    static bool test_string( const std::string &teststr, const TextFieldSearchType &t, const std::string &ss );
    /**
     \param teststr String representation of the date to compare to against the date user enetered in the GUI.
            If cant be converted to a valid date, returns false.
     \param test How to compare the dates
     \param ref_date Date user entered in the GUI that the file date should be tested against.
     */
    static bool test_date( const std::string &teststr, const NumericFieldMatchType &test, const boost::posix_time::ptime &ref_date );
    
  protected:
    FileDataField m_searchField;
    
    int m_discreteOption;
    double m_numeric;
    NumericFieldMatchType m_compareType;
    
    std::string m_searchString;
    TextFieldSearchType m_stringSearchType;
    
    boost::posix_time::ptime m_time;
  };//class SpecTest
  
  /** If you have an XML file in the same directory as spectrum files, the
   EventXmlTest allows you to test that XML against certain criteria to see if
   you should accept that spectrum file.  This XML file will be refered to as
   the Event XML file.
   
   Some conditions and limitations:
     -The Event XML file schema is dictated by the parameters passed into the
      #set_string_test_info or #set_date_test_info function.
     -The Event XML file size is restricted to be rather small (<64kb)
     -The same Event XML file applies to all spectrum files in a given directory
     -The Event XML file can be named anything, but must have the extension
      ".xml"
     -If there are multiple Event XML files in a given directory, only one of
      them will be used, with no specification of how that one is choosen.
     -An xpath is used to provide the string quantities to test against.
      (Note: If I dont get permission to add pugixml (which I want to start
      using instead of rapid XML anyway) then a very cripled subset of xpath
      will be used)
     -The Event XML is reparsed and tested for each spectrum file in a directory
     -All string comparisons are case insensitive.
   */
  class EventXmlTest
  {
  public:
    /** Default constructor (to play well with boost::any).
     You must call #set_string_test_info or set_date_test_info, or the #isvalid
     function will throw.
     */
    EventXmlTest();
    
    /** Function to set test conditions for.
     
       \param test_base_node The value of Event XML base node.  For example
              "event" for XML files that have a base-node of <event>.  This is
              used to help quickly discard invalid files without having to fully
              parse them.  If empty, all files will be tried.
       \param test_xpath The XPath expression to yeild the values the test
         string will compare against. May not be empty or #isvalid will throw.
       \param test_string The string to test the XPath yeiled values against.
              May not be empty or #isvalid will throw.
       \param comparison_type  The type of string comparison to do.
     */
    void set_string_test_info( const std::string &test_label,
                        const std::string &test_string,
                        const TextFieldSearchType &comparison_type );
    
    /** Same as #set_string_test_info, but for fields that are interpreted as a
        date.
     */
    void set_date_test_info( const std::string &test_label,
                              const boost::posix_time::ptime &test_time,
                              const NumericFieldMatchType &comparison_type );
    
    /** Function to perform the test.  Will search for a Event XML file in the
     spectrum files directory, and then perform the test.
     If no Event XML file is found, will return false.
     */
    bool test( const SpecFileInfoToQuery &meas ) const;
    
    /** Returns a summary string of the search. */
    std::string summary() const;
    
    //Throw exception with explanation if not valid
    void isvalid();
    
    /** Minimum file size an event XML file can be.  Defaults to 128 bytes. */
    static const size_t sm_min_event_xml_file_size;
    
    /** Maximum file size an event XML file can be.  Defaults to 64 kb. */
    static const size_t sm_max_event_xml_file_size;
    
  protected:
    
    enum class TestType{ String, Date, NotSet };
    
    TestType m_testType;
    
    std::string m_test_string;
    boost::posix_time::ptime m_test_time;
    std::string m_test_label;
    
    TextFieldSearchType m_fieldTestType;
    NumericFieldMatchType m_dateTestType;
  };//class EventXmlTest
  
  enum LogicType
  {
    LogicalOr,
    LogicalAnd,
    LogicalNot,
    LogicalOpenParan,
    LogicalCloseParan,
    NumLogicType
  };//enum LogicType
  
  
  const char *to_string( const LogicType type );
  
  
  //Holds all the contiditions
  class SpecLogicTest
  {
  public:
    SpecLogicTest();
    void addCondition( const SpecTest &test );
    void addLogic( const LogicType test );
    void addEventXmlTest( const EventXmlTest &test );
    

    bool test( const SpecFileInfoToQuery &meas ) const;
    
    //Throws exception if invalid
    void isvalid();
    
    std::string summary() const;
    
  protected:
    static bool evaluate( std::vector<boost::any> fields, const SpecFileInfoToQuery &meas );
    
    static std::ostream &print_equation( std::vector<boost::any> fields, std::ostream &strm );
    
    std::vector<boost::any> m_fields;  //Either LogicType or SpecTest
  };
}//namespace SpecFileQuery_h

#endif //
