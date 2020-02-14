#ifndef QLSpecMeas_h
#define QLSpecMeas_h

/* QLSpecMeas is a version of InterSpecs SpecMeas class that was
 branched 20171227 to create a version slightly more acceptable for QuickLook
 functionality.
 */


#include <set>
#include <deque>
#include <memory>

#include <Wt/WSignal>

#include "SpecUtils/SpecFile.h"

class QLPeakDef;
class QLPeakModel;
class Recalibrator;
struct QLPeakContinuum;

void log_error_message( const std::string &message, const std::string &source, const int priority );


class QLSpecMeas : public MeasurementInfo
{
/* Adds storing peaks to the MeasurementInfo class.  Used more extensively in
   InterSpec.
*/
  
public:
  typedef std::deque< std::shared_ptr<const QLPeakDef> > PeakDeque;
  typedef std::shared_ptr< PeakDeque >                 PeakDequeShrdPtr;
  typedef std::map<std::set<int>, PeakDequeShrdPtr >     SampleNumsToPeakMap;
  
  
public:
  //Default constructor for this class
  QLSpecMeas();

  //~QLSpecMeas():
  virtual ~QLSpecMeas();
  
  //load_N42_file(...): similar to MeasurmentInfo::load_from_N42(...),
  //  but also tries to read peak info (see also QLSpecMeas::write_2006_N42(...))
  virtual bool load_N42_file( const std::string &filename );
  
  //load_from_N42: loads spectrum from the raw XML data. The data must be
  //  null terminated.
  virtual bool load_from_N42( std::istream &input );

  //load_N42_from_data(...): raw xml file data - must be 0 terminated
  virtual bool load_N42_from_data( char *data );

  std::shared_ptr<Measurement> continuum();
  bool continuumVisible() const;

  const std::set<int> &displayedSampleNumbers() const;

  virtual void cleanup_after_load( const unsigned int flags
                                         = MeasurementInfo::StandardCleanup );
  
  std::shared_ptr< std::deque< std::shared_ptr<const QLPeakDef> > >
                                 peaks( const std::set<int> &samplenums );
  std::shared_ptr<const std::deque< std::shared_ptr<const QLPeakDef> > >
                                 peaks( const std::set<int> &samplenums ) const;
  
  std::set<std::set<int> > sampleNumsWithPeaks() const;
  
  //removeAllPeaks(): removes all peaks for all samplenumbers.
  //  Note: does not notify QLPeakModel, or anywhere else.
  void removeAllPeaks();
  
  //setPeaks(): sets peaks for the specified sample numbers, getting rid of
  //  any peaks that may have existed for those sample numbers previously.
  //  The peaks stored by *this will be same pointers to peaks as passed in,
  //  but the deque will be different than passed in.
  //  Note: does not notify QLPeakModel, or anywhere else.
  void setPeaks( std::deque< std::shared_ptr<const QLPeakDef> > &peakdeque,
                 const std::set<int> &samplenums );
  
  std::shared_ptr< const PeakDeque > automatedSearchPeaks(
                                        const std::set<int> &samplenums ) const;

  //peaksHaveBeenAdded(): marks this MeasurementInfo object as
  void setModified();

  
  //decodeSpecMeasStuffFromXml(...): Parent node should be named "DHS:InterSpec"
  //Throws on error.
  void decodeSpecMeasStuffFromXml( const ::rapidxml::xml_node<char> *parent );
  
  void addPeaksFromXml( const ::rapidxml::xml_node<char> *peaksnode );

private:
  //Do not call operator= or copy constructor, will (purposely) crash program.
  //  See makeUserSpaceObserver(...) and uniqueCopyContents(...).
  //Design rational: in various places we make copies of a QLSpecMeas object
  //  so we can 'undo' an action, or split a file; however, I dont want to
  //  carelessly call a sequence of operator= and then end up with a QLSpecMeas
  //  not hooked up to the notification signals that I want hooked up, therefore
  //  we have to explicitly call makeUserSpaceObserver(...) or
  //  uniqueCopyContents(...) to make a copy, to really understand the sequence
  //  of object copies
  QLSpecMeas( const QLSpecMeas &rhs );
  const QLSpecMeas &operator=( const QLSpecMeas &meas );

protected:
  //m_peaks: is only accessed by the QLPeakModel class - it shouldnt be set or
  //         modified anywhere else. When a new primary spectrum is loaded in
  //         InterSpec::setSpectrum(...), the m_peaks variable is passed
  //         to the QLPeakModel class, so it then controls inserting or removing
  //         peaks.
  
  enum XmlPeakSource{ UserPeaks, AutomatedSearchPeaks /*, AutomatedSearchInitialPeaks */ };
  static const char *toString( const XmlPeakSource source );
  static XmlPeakSource xmlPeakSourceFromString( const std::string &val );
  
  //Note that if the user select many permutaions of display sample numbers for
  //  a single spetrum file, the memmorry usage could kinda grow (but I expect
  //  a PeakDeque to be relatively small in size, in general)
  std::shared_ptr< SampleNumsToPeakMap > m_peaks;

  std::shared_ptr<std::set<int> > m_displayedSampleNumbers;

  
  SampleNumsToPeakMap m_autoSearchPeaks;
//  SampleNumsToPeakMap m_autoSearchInitialPeaks;
  
  static const int sm_peakXmlSerializationVersion;
};//class QLSpecMeas


#endif //SpecMeas_h
