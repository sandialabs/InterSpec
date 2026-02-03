#ifndef SpecFileQueryDbCache_h
#define SpecFileQueryDbCache_h
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

#include <set>
#include <mutex>
#include <ctime>
#include <memory>
#include <string>
#include <vector>
#include <condition_variable>

#include <Wt/Dbo/Dbo>
#include <Wt/Dbo/WtSqlTraits>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/FarmOptions.h"

//Forward declarations and Wt::Dbo::overhead ish.
namespace Wt {
  namespace Dbo {
    class Session;
    namespace backend {
      class Sqlite3;
    }
    
    template<>
    struct sql_value_traits<size_t, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind(size_t v, SqlStatement *statement, int column, int size);
      static bool read(size_t& v, SqlStatement *statement, int column, int size);
    };

    template<>
    struct sql_value_traits<std::set<float>, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind(const std::set<float> &v, SqlStatement *statement, int column, int size);
      static bool read(std::set<float> &v, SqlStatement *statement, int column, int size);
    };
        
    template<>
    struct sql_value_traits<std::set<size_t>, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind(const std::set<size_t> &v, SqlStatement *statement, int column, int size);
      static bool read(std::set<size_t> &v, SqlStatement *statement, int column, int size);
    };
    
    template<>
    struct sql_value_traits<std::set<SpecUtils::EnergyCalType>, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind(const std::set<SpecUtils::EnergyCalType> &v, SqlStatement *statement, int column, int size);
      static bool read(std::set<SpecUtils::EnergyCalType> &v, SqlStatement *statement, int column, int size);
    };
    
    template<>
    struct sql_value_traits<std::set<std::time_t>, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind(const std::set<std::time_t> &v, SqlStatement *statement, int column, int size);
      static bool read(std::set<std::time_t> &v, SqlStatement *statement, int column, int size);
    };
    
    template<>
    struct sql_value_traits<std::set<std::string>, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind(const std::set<std::string> &v, SqlStatement *statement, int column, int size);
      static bool read(std::set<std::string> &v, SqlStatement *statement, int column, int size);
    };
        
    template<>
    struct sql_value_traits<SpecUtils::DetectorAnalysis, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind( const SpecUtils::DetectorAnalysis &v, SqlStatement *statement, int column, int size);
      static bool read( SpecUtils::DetectorAnalysis &v, SqlStatement *statement, int column, int size );
    };
    
    template<>
    struct sql_value_traits<std::map<std::string,std::vector<std::string>>, void>
    {
      static const bool specialized = true;
      static std::string type(SqlConnection *conn, int size);
      static void bind( const std::map<std::string,std::vector<std::string>> &v, SqlStatement *statement, int column, int size);
      static bool read( std::map<std::string,std::vector<std::string>> &v, SqlStatement *statement, int column, int size );
    };
  }//namespace Dbo
}//namespace Wt

/** Information parsed from an external XML file to construct the GUI element
 that will be used to construct a #EventXmlTest test, if the user selects
 that option.
 */
struct EventXmlFilterInfo
{
  /** The type of input from the user*/
  enum class InputType
  {
    /** A drop down box where user can select one from many choices.  You
     must specify at least two #m_discreet_options.
     */
    Select,
    
    /** A user string input. */
    Text,
    
    /** A Date, e.g., 'yyyy/mm/dd' */
    Date
  };//enum class InputType
  
  std::string m_base_node_test; //String that should be in the first 128 bytes of the file. ex., "event"
  std::string m_label; //must be between 1 and 64 characters.  Must be unique wrt the widget.
  bool show_in_gui_table;
  std::string m_placeholder; //only for m_type==InputType::Text
  std::string m_xpath; //Must be valid xpath
  InputType m_type;
  std::vector<std::string> m_operators; //can only be from
  std::vector<std::string> m_discreet_options;  //m_type==InputType::Select must have at least two entries
  
  
  /** Parses a JSON string into EventXmlFilterInfo.
   Throws exception with message appropriate for displaying to user upon
   parsing error.  Also enforces relevant constraints and ensures text is
   escaped.
   
   ToDo: Document JSON format.  Untill then see implementation of function.
   */
  static std::vector<EventXmlFilterInfo> parseJsonString( const std::string &str );
};//struct EventXmlFilterInfo

/** Instead of determining query criteria directly from a SpecUtils::SpecFile, we
 will instead copy all relevant info to a SpecFileInfoToQuery struct, and make
 a decision off of that.  This allows us to cache the relevant information from
 a spectrum file, and save it to a sqlite3 database to greatly speed up future
 queries by avoiding having to re-parse the spectrum file.
 */
struct SpecFileInfoToQuery
{
  SpecFileInfoToQuery();
  void reset();
  void fill_info_from_file( const std::string filepath, const Farm::FarmOptions &farm_options );
  void fill_event_xml_filter_values( const std::string &filepath,
                                const std::vector<EventXmlFilterInfo> &xmlfilters );
  
  std::string file_path;
  long long int file_size;
  long long int file_path_hash;
  bool is_file;
  bool is_spectrum_file;
  bool is_event_xml_file;
  
  std::string filename;
  std::set<std::string> detector_names;
  std::string serial_number;
  std::string manufacturer;
  std::string model;
  std::string uuid;
  std::set<std::string> file_remarks;
  std::set<std::string> record_remarks;
  std::string location_name;
  
  bool has_riid_analysis;
  SpecUtils::DetectorAnalysis riid_ana;
  SpecUtils::DetectorType detector_type;
  bool passthrough;
  /** Total live-time of foreground spectra,  Background, intrinsic, calibration filtered out, unless no foreground identified.  */
  float total_livetime;
  /** Total real-time of foreground spectra,  Background, intrinsic, calibration filtered out, unless no foreground identified.  */
  float total_realtime;
  bool contained_neutron;
  bool contained_dev_pairs;
  bool contained_gps;
  bool contains_background;
  bool contains_calibration;
  bool contains_intrinsic;
  std::set<SpecUtils::EnergyCalType> energy_cal_types;
  std::set<float> individual_spectrum_live_time;
  std::set<float> individual_spectrum_real_time;
  size_t number_of_samples;
  size_t number_of_records;
  std::set<size_t> number_of_gamma_channels;
  std::set<float> max_gamma_energy;
  float mean_latitude;
  float mean_longitude;
  
  std::set<float> neutron_count_rate;
  std::set<float> gamma_count_rate;
  std::time_t start_time_ioi; //time of first non-background or intrinsic spectrum
  std::set<std::time_t> start_times;  //long; good to the second - good enough.
  
  /** For a particularly large dataset (~40k specfiles, 14k events),
      and maybe a dozen Event XML search criteria, caching the Event XML
      information increases initial search time from somethign like 3:40 (for
      no Event XML search criterial defined) to somewhere between 4:00 and 4:40.
      Also increases subsequent search times from 9 seconds, to 11 seconds, or
      18 seconds if a blanket xpath is defined for all event tag values.
   */
  std::map<std::string,std::vector<std::string>> event_xml_filter_values;

  // ============= FARM Analysis Fields =============

  // Peaks found by automated peak search (simplified summary)
  // JSON array: [{"mean":662.5,"fwhm":2.3,"amplitude":1000,"area":2500,"chi2dof":1.2}, ...]
  std::string farm_peaks_json;

  // Full peak JSON from PeakDef::peak_json() — ROI-grouped, includes continuum
  // coefficients, skew parameters, colors, etc.  Same format as D3 chart receives.
  std::string farm_peaks_full_json;

  // GADRAS Full Spectrum Isotope ID results as JSON (ExternalRidResults format)
  std::string gadras_rid_json;

  // Enrichment/isotopics results — JSON array of EnrichmentResults objects
  std::string isotopics_result_json;

  // Statistical moments of foreground spectrum
  double spectrum_mean = 0.0;      // Weighted mean energy (keV)
  double spectrum_variance = 0.0;  // Weighted variance
  double spectrum_skewness = 0.0;  // Third standardized moment
  double spectrum_kurtosis = 0.0;  // Fourth standardized moment (excess kurtosis)

  // Foreground channel/count summary
  int farm_min_channel_with_data = -1;   // First channel index with nonzero counts (-1 if none)
  int farm_max_channel_with_data = -1;   // Last channel index with nonzero counts  (-1 if none)
  double farm_foreground_total_gamma_counts = 0.0;
  int farm_foreground_num_gamma_channels = 0;
  bool farm_foreground_has_neutrons = false;
  double farm_foreground_neutron_count = 0.0;
  float farm_foreground_neutron_live_time = 0.0f;
  float farm_foreground_min_gamma_count = 0.0f;  // Minimum nonzero channel count
  float farm_foreground_max_gamma_count = 0.0f;

  // Energy calibration info as JSON object; empty when cal is invalid/unavailable.
  // Contains: polynomial_energy_coeffs, average_amp_gain (LCE only),
  //   nonlinear_deviation_pairs, min_energy, max_energy
  std::string farm_energy_cal_json;

  //bool peaks_fit;
  //std::vector<std::tuple<float,float,float>> peak_info; //mean, fwhm, rate

  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::id( a, file_path_hash, "file_path_hash" );
    Wt::Dbo::field( a, file_path, "file_path" );
    Wt::Dbo::field( a, file_size, "file_size" );
    Wt::Dbo::field( a, is_file, "is_file" );
    Wt::Dbo::field( a, is_spectrum_file, "is_spectrum_file" );
    Wt::Dbo::field( a, is_event_xml_file, "is_event_xml_file" );
    Wt::Dbo::field( a, filename, "filename" );
    Wt::Dbo::field( a, detector_names, "detector_names" );
    Wt::Dbo::field( a, serial_number, "serial_number" );
    Wt::Dbo::field( a, manufacturer, "manufacturer" );
    Wt::Dbo::field( a, model, "model" );
    Wt::Dbo::field( a, uuid, "uuid" );
    Wt::Dbo::field( a, file_remarks, "file_remarks" );
    Wt::Dbo::field( a, record_remarks, "record_remarks" );
    Wt::Dbo::field( a, location_name, "location_name" );
    Wt::Dbo::field( a, has_riid_analysis, "has_riid_analysis" );
    Wt::Dbo::field( a, riid_ana, "riid_ana" );
    Wt::Dbo::field( a, detector_type, "detector_type" );
    Wt::Dbo::field( a, passthrough, "passthrough" );
    Wt::Dbo::field( a, total_livetime, "total_livetime" );
    Wt::Dbo::field( a, total_realtime, "total_realtime" );
    Wt::Dbo::field( a, contained_neutron, "contained_neutron" );
    Wt::Dbo::field( a, contained_dev_pairs, "contained_dev_pairs" );
    Wt::Dbo::field( a, contained_gps, "contained_gps" );
    Wt::Dbo::field( a, contains_background, "contains_background" );
    Wt::Dbo::field( a, contains_calibration, "contains_calibration" );
    Wt::Dbo::field( a, contains_intrinsic, "contains_intrinsic" );
    Wt::Dbo::field( a, energy_cal_types, "energy_cal_types" );
    Wt::Dbo::field( a, individual_spectrum_live_time, "individual_spectrum_live_time" );
    Wt::Dbo::field( a, individual_spectrum_real_time, "individual_spectrum_real_time" );
    Wt::Dbo::field( a, number_of_samples, "number_of_samples" );
    Wt::Dbo::field( a, number_of_records, "number_of_records" );
    Wt::Dbo::field( a, number_of_gamma_channels, "number_of_gamma_channels" );
    Wt::Dbo::field( a, max_gamma_energy, "max_gamma_energy" );
    Wt::Dbo::field( a, mean_latitude, "mean_latitude" );
    Wt::Dbo::field( a, mean_longitude, "mean_longitude" );
    Wt::Dbo::field( a, neutron_count_rate, "neutron_count_rate" );
    Wt::Dbo::field( a, gamma_count_rate, "gamma_count_rate" );
    Wt::Dbo::field( a, start_time_ioi, "start_time_ioi" );
    Wt::Dbo::field( a, start_times, "start_times" );
    Wt::Dbo::field( a, event_xml_filter_values, "event_xml_filter_values" );

    // FARM fields
    Wt::Dbo::field( a, farm_peaks_json, "farm_peaks_json" );
    Wt::Dbo::field( a, farm_peaks_full_json, "farm_peaks_full_json" );
    Wt::Dbo::field( a, gadras_rid_json, "gadras_rid_json" );
    Wt::Dbo::field( a, isotopics_result_json, "isotopics_result_json" );
    Wt::Dbo::field( a, spectrum_mean, "spectrum_mean" );
    Wt::Dbo::field( a, spectrum_variance, "spectrum_variance" );
    Wt::Dbo::field( a, spectrum_skewness, "spectrum_skewness" );
    Wt::Dbo::field( a, spectrum_kurtosis, "spectrum_kurtosis" );
    Wt::Dbo::field( a, farm_min_channel_with_data, "farm_min_channel_with_data" );
    Wt::Dbo::field( a, farm_max_channel_with_data, "farm_max_channel_with_data" );
    Wt::Dbo::field( a, farm_foreground_total_gamma_counts, "farm_foreground_total_gamma_counts" );
    Wt::Dbo::field( a, farm_foreground_num_gamma_channels, "farm_foreground_num_gamma_channels" );
    Wt::Dbo::field( a, farm_foreground_has_neutrons, "farm_foreground_has_neutrons" );
    Wt::Dbo::field( a, farm_foreground_neutron_count, "farm_foreground_neutron_count" );
    Wt::Dbo::field( a, farm_foreground_neutron_live_time, "farm_foreground_neutron_live_time" );
    Wt::Dbo::field( a, farm_foreground_min_gamma_count, "farm_foreground_min_gamma_count" );
    Wt::Dbo::field( a, farm_foreground_max_gamma_count, "farm_foreground_max_gamma_count" );
    Wt::Dbo::field( a, farm_energy_cal_json, "farm_energy_cal_json" );

    //bool peaks_fit;
    //std::vector<std::tuple<float,float,float>> peak_info; //mean, fwhm, rate
  }//void persist( Action &a )
};//struct SpecFileInfoToQuery



/** This class parses a spectrum file and returns a
 #SpecFileInfoToQuery struct that can then be tested against
 a #SpecFileQuery::SpecTest.
 
 This class can be used either to parse and extract the information on demand
 when #spec_file_info is called, or it can cache the spectrum file information
 in a temporary database either at the first call of #spec_file_info or ahead
 of time via #cache_results, greatly speading up searches.
 */
class SpecFileQueryDbCache
{
public:
  /**
   @param use_db_caching Determines if a temporary database will be used to
          store the results of parsing spectrum files.
   */
  SpecFileQueryDbCache( const bool use_db_caching,
                        const std::string &base_path,
                        const std::vector<EventXmlFilterInfo> &xmlfilters );
  
  /** Class detructor makes sure that if #cache_results is running in another
      thread it will exit before the constructor closes.
   */
  ~SpecFileQueryDbCache();
  
  /** Function meant to be called from an axuilary thread to start parsing
   spectrum files and storing their #SpecFileInfoToQuery before the user clicks
   "Search".
   */
  void cache_results( const std::vector<std::string> &&files );
  
  /** Stops #cache_results if executing in another thread.
   Subsequent calls to #cache_results will immediately return until
   #allow_start_caching is called.
   */
  void stop_caching();
  
  /** Called after #stop_caching to re-enable caching. */
  void allow_start_caching();
  
  /** Returns whether cahcing to the database is enabled.
   */
  bool caching_enabled();
  
  /** Throws exception on error. */
  void set_persist( const bool persist );
  
  bool is_persistand() const;
  
  /** Sets the FARM analysis options used when populating SpecFileInfoToQuery.
   Must be called before caching/querying if FARM analysis is desired.
   */
  void set_farm_options( const Farm::FarmOptions &opts );

  /** Returns the XML filters that can be queried on.
   */
  const std::vector<EventXmlFilterInfo> &xml_filters() const;
  
  /** Returns the #SpecFileInfoToQuery information for a given
      file on the filesystem.  If database caching is enabled, will first check
      the database for the information, and if found, return that (assuming
      file size hasnt changed; this may switch to files modification date in the
      future).  If not from the database then the spectrum file will be parsed
      and information filled out from that; if DB caching is enabled the info
      will also be stored to the databsae.
   */
  std::unique_ptr<SpecFileInfoToQuery> spec_file_info( const std::string &filepath );

protected:
  bool open_db( const std::string &path, const bool create_tables );
  
  /** Returns true if a persisted, valid, cache DB was found for m_fs_path, and
      if if was, m_db and m_db_session are set.
      You should have a lock on m_db_mutex while calling this function.
   */
  bool init_existing_persisted_db();

  /** Constructs the persisted DB filename.  Appends "_FARM" when FARM analysis
      is enabled so that changing FARM options naturally invalidates the cache.
   */
  std::string construct_persisted_db_filename() const;

  /** Closes the current DB and opens a fresh one at the correct (possibly
      renamed) path.  Call after m_farm_options changes while caching is active.
   */
  void resetCache();
  
protected:
  bool m_use_db_caching;
  bool m_using_persist_caching;
  
  std::condition_variable m_cv;
  std::mutex m_cv_mutex;
  bool m_stop_caching;
  bool m_doing_caching;
  
  std::string m_fs_path;
  std::string m_db_location;
  
  std::mutex m_db_mutex;
  std::unique_ptr<Wt::Dbo::backend::Sqlite3> m_db;
  std::unique_ptr<Wt::Dbo::Session> m_db_session;
  
  const std::vector<EventXmlFilterInfo> m_xmlfilters;
  Farm::FarmOptions m_farm_options;
};//class SpecFileQueryDbCache


#endif //SpecFileQueryDbCache_h
