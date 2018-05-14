#ifndef SpecFileQueryWidget_h
#define SpecFileQueryWidget_h
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

#include <string>
#include <vector>
#include <memory>
#include <atomic>
#include <iostream>

#include <boost/any.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <Wt/WContainerWidget>


class Measurement;
class InterSpec;
class MeasurementInfo;
class SpecMeasManager;
class ResultTableModel;
class RowStretchTreeView;

namespace Wt
{
  class WText;
  class WAnchor;
  class WSpinBox;
  class WCheckBox;
  class WLineEdit;
  class WPushButton;
  class WApplication;
}


namespace SpecQuery
{
  //Stuff in this namespace is intended to be seperate from the GUI and contain
  // most of the functional/logical code so it can eventually be used to form a
  // command line utilitiy or something independant from the GUI.
  
  
  enum FileDataField
  {
    //String Fields to search
    Filename,
    DetectorName,
    SerialNumber,
    Manufacturer,
    Model,
    Uuid,
    Remark,
    LocationName,
    AnalysisResultText,
    AnalysisResultNuclide,
  
    //Discrete options
    DetectionSystemType,
    SearchMode,
    //More discrete types: ContainedNuetronDetector, ContainedDeviationPairs, ParserType, InstrumentType, HasGps
 
    //Ranges of numeric values
    TotalLiveTime,
    TotalRealTime,
    IndividualSpectrumLiveTime,
    IndividualSpectrumRealTime,
    //More numeric types: NumberOfSamples, NumberOfRecords, NumberOfGammaChannels, MaximumGammaEnergy, Latitude, Longitude,
    
    //Time field
    StartTime,
    
    NumFileDataFields
  };//enum FileDataField
  
  const char *to_string( const FileDataField field );
  
  enum TextFieldSearchType
  {
    TextIsExact,
    TextIsContained,
    TextStartsWith,
    TextEndsWith,
    TextRegex
  };//enum TextFieldSearchType
  
  const char *to_string( const TextFieldSearchType type );
  
  enum NumericFieldMatchType
  {
    ValueIsExact,
    ValueIsLessThan,
    ValueIsGreaterThan
  };//enum NumericFieldMatchType
  
  const char *to_string( const NumericFieldMatchType type );
  
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
    
    bool test( const std::shared_ptr<const MeasurementInfo> meas ) const;
    
    std::string summary() const;
    
    //Throw exception with explanation if not valid
    void isvalid();
    
  protected:
    
    static bool test_string( const std::string &teststr, const TextFieldSearchType &t, const std::string &ss );
    
  protected:
    FileDataField m_searchField;
    
    int m_discreteOption;
    double m_numeric;
    NumericFieldMatchType m_compareType;
    
    std::string m_searchString;
    TextFieldSearchType m_stringSearchType;
    
    boost::posix_time::ptime m_time;
  };//class SpecTest

  
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
  class SpecFileQuery
  {
  public:
    SpecFileQuery();
    void addCondition( const SpecTest &test );
    void addLogic( const LogicType test );
    
    //
    bool test( const std::shared_ptr<const MeasurementInfo> meas ) const;
    
    //Throws exception if invalid
    void isvalid();
    
    std::string summary() const;
    
  protected:
    static bool evaluate( std::vector<boost::any> fields,
                          const std::shared_ptr<const MeasurementInfo> &meas );
    static std::ostream &print_equation( std::vector<boost::any> fields, std::ostream &strm );
    
    std::vector<boost::any> m_fields;  //Either LogicType or SpecTest
  };
}//namespace SpecQuery




class SpecFileQueryWidget : public Wt::WContainerWidget
{
public:
  SpecFileQueryWidget( InterSpec *viewer, Wt::WContainerWidget *parent = 0 );
  virtual ~SpecFileQueryWidget();
  
protected:
  void init();
  
  /** Counts the number of candidate files to query, in the given (optionally
      recursive) source directory, filtering them by file max size (MB) and file
      extension preferences.
      Results are are posted to WServer using the specified WApplication
      sessionid, calling updateNumberFilesInGui to actually update the GUI.
      This function can run in any thread.
      \param widgetdeleted should be a copy of #m_widgetDeleted
   */
  static void updateNumberFiles( const std::string srcdir, const bool recursive,
                          const bool extfilter, const size_t maxsize,
                          SpecFileQueryWidget *querywidget,
                          const std::string sessionid,
                          std::shared_ptr< std::atomic<bool> > widgetdeleted );
  
  /** Updates the GUI for the found number of candidate files, if the source
      directory, recursice and extension filter options are still all the same.
      This function must be run from the main WApplication thread for this
      instance.
      \param widgetdeleted should be a copy of #m_widgetDeleted to determine
      if this SpecFileQueryWidget had gotten deleted while the file search was
      happening.
   */
  static void updateNumberFilesInGui( const size_t nfiles,
                                      const std::string srcdir,
                                      const bool recursive,
                                      const bool extfilter,
                                      SpecFileQueryWidget *querywidget,
                               std::shared_ptr< std::atomic<bool> > widgetdeleted );
  
  void setResultFieldVisibility( const SpecQuery::FileDataField field, const bool visible );

  void basePathChanged();
  void setResultsStale();
  
  void addCondition();
  void addConditionAt( int index );
  void addConditionAfter( Wt::WWebWidget *ww );
  void removeCondition( Wt::WWebWidget *ww );
  
  void startUpdate();
  void finishUpdate( std::shared_ptr< std::vector< std::vector<std::string> > > result,
                     const std::string description,
                     const bool wasCanceled,
                     std::shared_ptr< std::atomic<bool> > widgetDeleted );
  void cancelUpdate();
  void updateSearchStatus( const size_t nfilestotal, const size_t nfileschecked,
                           const size_t nfilesaccepted,
                           const std::string specialmsg,
                           std::shared_ptr< std::atomic<bool> > widgetDeleted );
  
  void doSearch( const std::string basedir, unsigned long options,
                 const size_t maxsize,
                 const SpecQuery::SpecFileQuery query,
                 const std::string sesssionID,
                 std::shared_ptr< std::atomic<bool> > stopUpdate,
                 std::shared_ptr< std::atomic<bool> > widgetDeleted );
  
  void selectionChanged();
  void loadSelected();
  
  SpecQuery::SpecFileQuery query();
  
protected:
  
  enum PrefilterOptions
  {
    FilterByFilename = 0x1,
    FilterDuplicates = 0x2,
    SearchRecursive  = 0x4
  };//enum PrefilterOptions
  
  std::shared_ptr< std::atomic<bool> > m_stopUpdate;
  std::shared_ptr< std::atomic<bool> > m_widgetDeleted;
  
  Wt::WApplication * const m_app;
  InterSpec * const m_viewer;
  
  Wt::WContainerWidget *m_conditions;
  Wt::WText *m_addCondition;
  Wt::WPushButton *m_update;
  Wt::WPushButton *m_cancelUpdate;
  Wt::WPushButton *m_loadSelectedFile;
  
  ResultTableModel *m_resultmodel;
  RowStretchTreeView *m_resultview;
  
  Wt::WLineEdit *m_baseLocation;
  Wt::WCheckBox *m_recursive;
  Wt::WCheckBox *m_filterByExtension;
  Wt::WSpinBox *m_maxFileSize;
  Wt::WCheckBox *m_filterUnique;
  
  Wt::WText *m_numberFiles;
  Wt::WText *m_numberResults;
  Wt::WAnchor *m_csv;
};//class SpecFileQueryWidget


#endif
