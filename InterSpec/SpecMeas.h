#ifndef SpecMeas_h
#define SpecMeas_h
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
#include <deque>
#include <memory>

#include <Wt/WSignal>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"


//Temporary include - only used to get RAPIDXML_USE_SIZED_INPUT_WCJOHNS macro
#include "external_libs/SpecUtils/3rdparty/rapidxml/rapidxml.hpp"

class PeakDef;
class PeakModel;
class Recalibrator;
struct PeakContinuum;
class DetectorPeakResponse;

namespace rapidxml
{
  template<class Ch> class xml_document;
}

class SpecMeas : public SpecUtils::SpecFile
{
/*
  This SpecMeas class is the class which hold not just the SpecUtils::SpecFile
  information via inheritance (the SpecUtils::SpecFile holds all info that comes
  from a measurment file, sucha as a PCF or N42 file), but all the user
  'value added' information, such as peak information, highlighted regions,
  etc.

  I do have a bit intrepidation about integrating this class so tightly with
  all the other classes in this project.
*/
  
public:
  typedef std::deque< std::shared_ptr<const PeakDef> > PeakDeque;
  typedef std::shared_ptr< PeakDeque >                 PeakDequeShrdPtr;
  typedef std::map<std::set<int>, PeakDequeShrdPtr >     SampleNumsToPeakMap;
  
  //
  //typedef std::map<std::pair<std::set<int>,std::set<std::string>>, PeakDequeShrdPtr >     SampleNumsToPeakMap;
  
public:
  //Default constructor for this class
  SpecMeas();

  //uniqueCopyContents(...):  makes *this a copy copy of 'rhs'.  this->m_peaks,
  //  this->m_detector, etc. will not point to the objects owned by 'rhs', but
  //  rather unique copies of those objects.
  //If *this has any "UserSpaceObservers", then they _will_ "see" the changes,
  //  and future changes to *this.
  //Also, all the measurment objects of *this and rhs will point to different
  //  locations in memorry (although channel data and bin energies will have
  //  same shared const pointers).
  //
  //  Note 20180426: This function is only used when INCLUDE_ANALYSIS_TEST_SUITE
  //    And can probably be eliminated (it is a vestigial of a previous (poor)
  //    design)
  void uniqueCopyContents( const SpecMeas &rhs );

  //~SpecMeas(): will emit the aboutToBeDeleted() signal
  virtual ~SpecMeas();
  
  //load_N42_file(...): similar to MeasurmentInfo::load_from_N42(...),
  //  but also tries to read peak info (see also SpecMeas::write_2006_N42(...))
  virtual bool load_N42_file( const std::string &filename );
  
  //load_from_N42: loads spectrum from the raw XML data. The data must be
  //  null terminated.
  virtual bool load_from_N42( std::istream &input );
  
  //load_N42_from_data(...): raw xml file data - must be 0 terminated
  virtual bool load_N42_from_data( char *data );
  
  /** Loads N42 file from raw XML file data specified by begin and end of data
      (does not need to be null terminated)
   */
  virtual bool load_N42_from_data( char *data, char *data_end );
  
  virtual void load_N42_from_doc( rapidxml::xml_document<char> &doc );
  
  virtual bool save2012N42File( const std::string &filename );
  virtual void save2012N42File( const std::string &filename,
                              boost::function<void()> error_callback ); //usefull when saving in a detached slave thread - be carful of race conditions (eg make error_callback will remain valid)

  
  //write_2006_N42(...): writes a 2006 N42 simple spectrometer file (similar to
  //  Cambio, see SpectrumInfo::write_2006_N42(...)), but with the addition of
  //  m_peaks information.
  //  Should also consider adding: m_detector, m_continuum, m_continuumVisible,
  //                               m_displayType, m_displayedSampleNumbers;
  virtual bool write_2006_N42( std::ostream &ostr ) const;

  //write_iaea_spe(...) Writes ascii IAEA SPE file; the specialization in
  //  SpecMeas also adds in peak information.
  virtual bool write_iaea_spe( std::ostream &output,
                               std::set<int> sample_nums,
                               const std::set<int> &det_nums ) const;
  
  virtual std::shared_ptr< ::rapidxml::xml_document<char> > create_2012_N42_xml() const;
  
  
  //guessDetectorTypeFromFileName(...): not called by default
  static SpecUtils::DetectorType guessDetectorTypeFromFileName( std::string name );
  

  std::shared_ptr<SpecUtils::Measurement> continuum();
  bool continuumVisible() const;
  std::shared_ptr<DetectorPeakResponse> detector();
  std::shared_ptr<const DetectorPeakResponse> detector() const;

  /**
  \deprecated Not consistently used if SpecMeas is used for foreground and background or whatever.
  */
  SpecUtils::SpectrumType displayType() const;
  
  /**
   \deprecated Not consistently used if SpecMeas is used for foreground and background or whatever.
   */
  const std::set<int> &displayedSampleNumbers() const;

  //aboutToBeDeleted(): signal emited right before object destructions - useful
  //  if you want to serialize changes to disk
  Wt::Signal<> &aboutToBeDeleted();

  
  static void save2012N42FileInSlaveThread( std::shared_ptr<SpecMeas> info,
                                            const std::string filename,
                                            boost::function<void()> error_callback );

  //setDetector(): set not only the detector of *this, but also of all of its
  //  observers, so they all point to the same object in memmory
  void setDetector( std::shared_ptr<DetectorPeakResponse> det );

  //detectorChangedCallback(): right now just calls setDetector(), but may
  //  change in the future or be removed.
  void detectorChangedCallback( std::shared_ptr<DetectorPeakResponse> det );


  void displayedSpectrumChangedCallback( SpecUtils::SpectrumType type,
                                         std::shared_ptr<SpecMeas> measurment,
                                         std::set<int> sample_numbers,
                                         std::vector<std::string> detectors  );
  virtual void cleanup_after_load( const unsigned int flags
                                         = SpecFile::StandardCleanup );
  
  std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > >
                                 peaks( const std::set<int> &samplenums );
  /** CAUTION: the returned deque may be modified elsewhere in the app, if for
   example the user adds or removes a peak.  However, the individual peaks will
   not be modified.  So if you plan to use the returned deque from outside of
   the main thread, you should make a new deque and fill it with pointers to the
   peaks.
   ToDo: I fixed one place I ran into the above issue, however this function
         should be rewritten to return a unique deque, instead of a pointer to
         a shared one.
   */
  std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > >
                                 peaks( const std::set<int> &samplenums ) const;
  
  std::set<std::set<int> > sampleNumsWithPeaks() const;
  
  //removeAllPeaks(): removes all peaks for all samplenumbers.
  //  Note: does not notify PeakModel, or anywhere else.
  void removeAllPeaks();
  
  //setPeaks(): sets peaks for the specified sample numbers, getting rid of
  //  any peaks that may have existed for those sample numbers previously.
  //  The peaks stored by *this will be same pointers to peaks as passed in,
  //  but the deque will be different than passed in.
  //  Note: does not notify PeakModel, or anywhere else.
  void setPeaks( const std::deque< std::shared_ptr<const PeakDef> > &peakdeque,
                 const std::set<int> &samplenums );
  
  std::shared_ptr< const PeakDeque > automatedSearchPeaks(
                                        const std::set<int> &samplenums ) const;
//  std::shared_ptr< const PeakDeque > automatedSearchInitialPeaks(
//                                            const std::set<int> &samplenums );
  
  void setAutomatedSearchPeaks( const std::set<int> &samplenums,
                                std::shared_ptr< PeakDeque > foundPeaks
                                /*, std::shared_ptr< PeakDeque > intitalPeaks*/ );
  
  //peaksHaveBeenAdded(): marks this SpecUtils::SpecFile object as
  void setModified();

  
  //shiftPeaksForRecalibration: shift the peaks for when you apply a
  //  recalibration to the spectrum, for instance after calling
  //  SpecFile::recalibrate_by_eqn(...).  Note that the PeakModel is not
  //  notified of the changes, and the shift is applied to all peaks of the
  //  SpecMeas object.
  //void shiftPeaksForRecalibration( std::vector<float> old_pars,
  //                  const std::vector< std::pair<float,float> > &old_devpairs,
  //                  SpecUtils::EnergyCalType old_eqn_type,
  //                  std::vector<float> new_pars,
  //                  const std::vector< std::pair<float,float> > &new_devpairs,
  //                  SpecUtils::EnergyCalType new_eqn_type );
  
  //translatePeakForCalibrationChange(...): translates a peaks definition from
  //  one calibration to another.
  //  XXX - Currently does not translate skew values (also, uncertainties may
  //        not be handled correctly either)!
  //
  //ToDo: Only handles polynomial or full range fraction calibration type;
  //      should handle lower channel energy.
  //static void translatePeakForCalibrationChange( PeakDef &peak,
  //                  std::vector<float> old_pars,
  //                  const std::vector< std::pair<float,float> > &old_devpairs,
  //                  SpecUtils::EnergyCalType old_eqn_type,
  //                  std::vector<float> new_pars,
  //                  const std::vector< std::pair<float,float> > &new_devpairs,
  //                  SpecUtils::EnergyCalType new_eqn_type,
  //                  const size_t nbins,
  //                  const bool translate_continuum );
  
  //appenSpecMeasStuffToXml(...): Adds the SpecMeas specific stuff (detector, 
  //  peaks, etc) to a node named <DHS:InterSpec> under 'parent'.
  //  Returns the node appended to parent.
  ::rapidxml::xml_node<char> *appendSpecMeasStuffToXml( 
                                ::rapidxml::xml_node<char> *parent ) const;
  
  
  //appendSampleNumbersToXml(...): appends sample number XML to parent, returns
  //  appended node (or null if !m_displayedSampleNumbers )
  ::rapidxml::xml_node<char> *appendSampleNumbersToXml(
                                    ::rapidxml::xml_node<char> *parent ) const;
  
  //appendDisplayedDetectorsToXml(...): appends displayed detector names XML to parent, returns
  //  appended node (or null if !m_displayedSampleNumbers )
  ::rapidxml::xml_node<char> *appendDisplayedDetectorsToXml(
                                    ::rapidxml::xml_node<char> *parent ) const;
  
  
  //decodeSpecMeasStuffFromXml(...): Parent node should be named "DHS:InterSpec"
  //Throws on error.
  void decodeSpecMeasStuffFromXml( const ::rapidxml::xml_node<char> *parent );
  
  
  void addPeaksToXml( ::rapidxml::xml_node<char> *peaksnode ) const;
  void addPeaksFromXml( const ::rapidxml::xml_node<char> *peaksnode );
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  //equalEnough(...): tests whether the passed in Measurement objects are
  //  equal, for most intents and purposes.  Allows some small numerical
  //  rounding to occur.
  //Throws an std::exception with a brief explanaition when an issue is found.
  static void equalEnough( const SpecMeas &lhs, const SpecMeas &rhs );
#endif
  
  
  /** Returns the shielding source model associated with thsi SpecMeas
      Be carefule, as you should take a lock on this SpecMeass mutex to
      make sure another thread doesnt change the model on you.
   */
  rapidxml::xml_document<char> *shieldingSourceModel();
  
  /** Sets the shielding source model that was serialized by the gui. */
  void setShieldingSourceModel( std::unique_ptr<rapidxml::xml_document<char>> &&model );
  
  
private:
  //Do not call operator= or copy constructor, will (purposely) crash program.
  //  See makeUserSpaceObserver(...) and uniqueCopyContents(...).
  //Design rational: in various places we make copies of a SpecMeas object
  //  so we can 'undo' an action, or split a file; however, I dont want to
  //  carelessly call a sequence of operator= and then end up with a SpecMeas
  //  not hooked up to the notification signals that I want hooked up, therefore
  //  we have to explicitly call makeUserSpaceObserver(...) or
  //  uniqueCopyContents(...) to make a copy, to really understand the sequence
  //  of object copies
  SpecMeas( const SpecMeas &rhs );
  const SpecMeas &operator=( const SpecMeas &meas );

protected:
  //m_peaks: is only accessed by the PeakModel class - it shouldnt be set or
  //         modified anywhere else. When a new primary spectrum is loaded in
  //         InterSpec::setSpectrum(...), the m_peaks variable is passed
  //         to the PeakModel class, so it then controls inserting or removing
  //         peaks.
  
  enum XmlPeakSource{ UserPeaks, AutomatedSearchPeaks /*, AutomatedSearchInitialPeaks */ };
  static const char *toString( const XmlPeakSource source );
  static XmlPeakSource xmlPeakSourceFromString( const std::string &val );
  
  static void addPeaksToXmlHelper(
                  const SpecMeas::SampleNumsToPeakMap &peaks,
                  const XmlPeakSource source,
                  rapidxml::xml_node<char> *peaksnode,
                  std::map<std::shared_ptr<PeakContinuum>,int> &continuumids,
                  std::map<std::shared_ptr<const PeakDef>,int> &peakids );

  
  //Note that if the user select many permutaions of display sample numbers for
  //  a single spetrum file, the memmorry usage could kinda grow (but I expect
  //  a PeakDeque to be relatively small in size, in general)
  std::shared_ptr< SampleNumsToPeakMap > m_peaks;
  std::shared_ptr<DetectorPeakResponse> m_detector;

  /** \TODO: Currently m_displayType, m_displayedSampleNumbers, and m_displayedDetectors track
   foreground, or background, or secondary, but this doesnt cover case where this SpecMeas is
   being used for two of these catagories... Should upgrade this mechanism to track all of these
   cases, or get rid of it totally.
   */
  std::shared_ptr<SpecUtils::SpectrumType> m_displayType;
  std::shared_ptr<std::set<int> > m_displayedSampleNumbers;
  std::shared_ptr<std::vector<std::string> > m_displayedDetectors;

  
  SampleNumsToPeakMap m_autoSearchPeaks;
//  SampleNumsToPeakMap m_autoSearchInitialPeaks;
  
  std::unique_ptr<rapidxml::xml_document<char>> m_shieldingSourceModel;
  
  Wt::Signal<> m_aboutToBeDeleted;
  
  /** ToDo/hack: we are currently using sample numbers to match peaks fit by
   InterSpec up to specific Measurement's.  This variable tells us if we need
   to strickly enforce the sample numbers found in the N42 file attributes
   (which is a non-standad thing we are doing).  The real solution is to use
   MeasurementGroupReferences, for both InterSpec specific stuff, but also for
   analysis results.
   */
  bool m_fileWasFromInterSpec;
  
  /** Version of XML serialization of the <DHS:InterSpec> node.
   Changes:
   - Added version field to xml 20200807, with initial value 1.  Added <DisplayedDetectors> field.
   
   */
  static const int sm_specMeasSerializationVersion;
  
  static const int sm_peakXmlSerializationVersion;
};//class SpecMeas


#endif //SpecMeas_h
