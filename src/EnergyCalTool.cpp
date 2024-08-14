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

#include <map>
#include <string>
#include <vector>
#include <ctime>
#include <memory>
#include <iostream>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WResource>
#include <Wt/WCheckBox>
#include <Wt/WFileUpload>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/EnergyCalGraphical.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/EnergyCalAddActions.h"
#include "InterSpec/IsotopeSelectionAids.h"


using namespace std;
using namespace Wt;

#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif


namespace
{
template<class T> struct index_compare_assend
{
  index_compare_assend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
  bool operator()(const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }
  const T arr;
};//struct index_compare
  
}// namepsace


namespace
{
  // For undo/redo, we will use a `EnergyCalUndoRedoSentry` struct to combine mutliple operations
  //  (like multiple calls to #EnergyCalTool::setEnergyCal or #EnergyCalTool::addDeviationPair)
  //  as part of a single user operation.  #EnergyCalUndoRedoSentry is a thread local object that
  //  allows you to construct as many of them as you would like, and any changes from any of them,
  //  will be combined together, and when the last #EnergyCalUndoRedoSentry destructs, an undo/redo
  //  step will be inserted into the history.
  //
  
  // For undo/redo, we will need to store the mappings, for each measurement, from old to new
  //  energy calibration.  Using weak ptrs for Measurement, although I dont think it is strickly
  //  necassary.
  //  TODO: store energy calibration a little more compactly; could store just the coefficients, and not the channel lower-energies
  typedef vector< tuple<weak_ptr<const SpecUtils::Measurement>, \
            shared_ptr<const SpecUtils::EnergyCalibration>, \
            shared_ptr<const SpecUtils::EnergyCalibration>> > \
          meas_old_new_cal_t;
  
  // For undo/redo, keep track of new and old peaks
  //  TODO: translate peaks on-the-fly, to reduce memory use, but will need to track pre-post energy cal to enable this
  typedef vector< tuple< set<int>, \
          deque< std::shared_ptr<const PeakDef> >, \
          deque< std::shared_ptr<const PeakDef> > > >
        meas_old_new_peaks_t;
  
  /** Function that actually does the work of undo/redo. */
  void do_undo_or_redo( const bool is_undo,
                        const SpecUtils::SpectrumType type,
                        const meas_old_new_peaks_t &meas_old_new_peaks,
                        const meas_old_new_peaks_t &meas_old_new_hint_peaks,
                        const meas_old_new_cal_t &meas_old_new_cal,
                        const std::weak_ptr<SpecMeas> &specfile_weak )
  {
    using namespace SpecUtils;
    
    const shared_ptr<SpecMeas> specfile = specfile_weak.lock();
    assert( specfile );
    InterSpec * const viewer = InterSpec::instance();
    assert( viewer );
    if( !viewer )
      return;

    const shared_ptr<SpecMeas> specfile_now = viewer ? viewer->measurment( type ) : nullptr;
    assert( specfile == specfile_now );
    if( !specfile || (specfile != specfile_now) )
    {
      Wt::log("error") << "SpecFile not same, as was expected during undo/redo.";
      return;
    }
    
    const shared_ptr<SpecMeas> foreground = viewer->measurment( SpecUtils::SpectrumType::Foreground );
    const set<int> &foresamples = viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );
    
  
    EnergyCalTool *tool = viewer->energyCalTool();
    PeakModel *peakModel = viewer->peakModel();
    assert( tool && peakModel );
    if( !tool || !peakModel )
    {
      Wt::log("error") << "Failed to get EnergyCalTool or PeakModel during undo/redo of dev. pairs.";
      return;
    }
    
    for( const auto &m_o_n : meas_old_new_cal )
    {
      const shared_ptr<const Measurement> lm = get<0>(m_o_n).lock();
      assert( lm );
      if( !lm )
      {
        Wt::log("error") << "Failed to get Measurement during undo/redo of dev. pairs.";
        continue;
      }
      
      const auto &from_cal = is_undo ? get<2>(m_o_n) : get<1>(m_o_n);
      const auto &to_cal = is_undo ? get<1>(m_o_n) : get<2>(m_o_n);
      assert( to_cal );
      
      shared_ptr<const Measurement> m = specfile->measurement( lm->sample_number(), lm->detector_name() );
      assert( m && (lm == m) );
      
      if( lm != m )
      {
        Wt::log("error") << "Failed to update a Measurement during undo/redo of dev. pairs.";
        continue;
      }
      
      assert( lm && (lm->energy_calibration() == from_cal) );
      specfile->set_energy_calibration( to_cal, m );
    }//for( const auto &m_o_n : meas_old_new_cal )
      
    
    for( const auto &m_o_n : meas_old_new_peaks )
    {
      const set<int> &samples = get<0>(m_o_n);
      //const deque<shared_ptr<const PeakDef>> &from_peaks = is_undo ? get<2>(m_o_n) : get<1>(m_o_n);
      const deque<shared_ptr<const PeakDef>> &to_peaks = is_undo ? get<1>(m_o_n) : get<2>(m_o_n);
      
      specfile->setPeaks( to_peaks, samples );
      if( peakModel && (specfile == foreground) && (samples == foresamples) )
        peakModel->setPeakFromSpecMeas(foreground, foresamples);
    }//for( loop over changed peaks )
    
    for( const auto &m_o_n : meas_old_new_hint_peaks )
    {
      const set<int> &samples = get<0>(m_o_n);
      const deque<shared_ptr<const PeakDef>> &to_peaks = is_undo ? get<1>(m_o_n) : get<2>(m_o_n);
      
      auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( to_peaks );
      specfile->setAutomatedSearchPeaks( samples, peaks );
    }//for( loop over changed peaks )
    
    viewer->refreshDisplayedCharts();
    tool->refreshGuiFromFiles();
  }//void do_undo_or_redo(...)
  
  
  /** We may make multiple calls to #EnergyCalTool::setEnergyCal and or #EnergyCalTool::addDeviationPair, or other function, for
   a single user-instigated change.  We want to combine multiple of these calls, so its a single user undo/redo operation, so we'll use
   thread-local storage, and this `EnergyCalUndoRedoSentry` object to aggregate all the changes, for a single user operation,
   and have the destructor of #EnergyCalUndoRedoSentry actually insert the undo/redo step.
  */
  
  typedef map<weak_ptr<SpecMeas>,meas_old_new_cal_t,owner_less<weak_ptr<SpecMeas>>> SpecMeasToCalHistoryMap;
  typedef map<weak_ptr<SpecMeas>,meas_old_new_peaks_t,owner_less<weak_ptr<SpecMeas>>> SpecMeasToPeakHistoryMap;
  
  // If we were using c++17, we could declare the following as `inline` variables of
  //  EnergyCalUndoRedoSentry, but for the moment we will just declare outside the class.
  //  (or instead we could make `cal_info` and `peak_info` static functions, but this defaeats
  //   thier purpose a bit)
  thread_local static int sm_undo_redo_level;
  thread_local static unique_ptr<SpecMeasToCalHistoryMap> sm_meas_old_new_cal_map;
  thread_local static unique_ptr<SpecMeasToPeakHistoryMap> sm_meas_old_new_peaks_map;
  thread_local static unique_ptr<SpecMeasToPeakHistoryMap> sm_meas_old_new_hint_peaks_map;
  
  struct EnergyCalUndoRedoSentry
  {
    EnergyCalUndoRedoSentry()
    {
      if( !sm_meas_old_new_cal_map )
      {
        sm_undo_redo_level = 1;
        sm_meas_old_new_cal_map = make_unique<SpecMeasToCalHistoryMap>();
        sm_meas_old_new_peaks_map = make_unique<SpecMeasToPeakHistoryMap>();
        sm_meas_old_new_hint_peaks_map = make_unique<SpecMeasToPeakHistoryMap>();
      }else
      {
        sm_undo_redo_level += 1;
        assert( sm_meas_old_new_peaks_map && sm_meas_old_new_hint_peaks_map );
      }
    }//EnergyCalUndoRedoSentry()
      
    meas_old_new_cal_t &cal_info( const shared_ptr<SpecMeas> &meas )
    {
      assert( sm_meas_old_new_cal_map );
      if( !sm_meas_old_new_cal_map )
        throw logic_error( "EnergyCalUndoRedoSentry: cal info not initied?" );
      
      SpecMeasToCalHistoryMap &m = *sm_meas_old_new_cal_map;
      return m[meas];
    }
    
    meas_old_new_peaks_t &peak_info( const shared_ptr<SpecMeas> &meas )
    {
      assert( sm_meas_old_new_peaks_map );
      if( !sm_meas_old_new_peaks_map )
        throw logic_error( "EnergyCalUndoRedoSentry: peak info not inited?" );
      
      SpecMeasToPeakHistoryMap &m = *sm_meas_old_new_peaks_map;
      return m[meas];
    }
    
    meas_old_new_peaks_t &hint_peak_info( const shared_ptr<SpecMeas> &meas )
    {
      assert( sm_meas_old_new_hint_peaks_map );
      if( !sm_meas_old_new_hint_peaks_map )
        throw logic_error( "EnergyCalUndoRedoSentry: hint peak info not inited?" );
      
      SpecMeasToPeakHistoryMap &m = *sm_meas_old_new_hint_peaks_map;
      return m[meas];
    }
    
    
    ~EnergyCalUndoRedoSentry()
    {
      using namespace SpecUtils;
      
      sm_undo_redo_level -= 1;
      assert( sm_meas_old_new_cal_map && sm_meas_old_new_peaks_map && sm_meas_old_new_hint_peaks_map );
      if( !sm_meas_old_new_cal_map || !sm_meas_old_new_peaks_map || !sm_meas_old_new_hint_peaks_map )
        return;
      
      if( sm_undo_redo_level > 0 )
        return;
      
      assert( sm_undo_redo_level == 0 );
      
      // Note, with C++14, we could capture unique_ptr into the lambdas, but we'll worry about that later
      const SpecMeasToCalHistoryMap meas_old_new_cal_map( std::move(*sm_meas_old_new_cal_map) );
      const SpecMeasToPeakHistoryMap meas_old_new_peaks_map( std::move(*sm_meas_old_new_peaks_map) );
      const SpecMeasToPeakHistoryMap meas_old_new_hint_peaks_map( std::move(*sm_meas_old_new_hint_peaks_map) );
      
      sm_meas_old_new_cal_map.reset();
      sm_meas_old_new_peaks_map.reset();
      sm_meas_old_new_hint_peaks_map.reset();
      
      // No changes registered.
      if( meas_old_new_cal_map.empty() 
         && meas_old_new_peaks_map.empty()
         && meas_old_new_hint_peaks_map.empty() )
      {
        return;
      }
      
      InterSpec *viewer = InterSpec::instance();
      assert( viewer );
      UndoRedoManager *undoManager = viewer ? viewer->undoRedoManager() : nullptr;
      if( !undoManager )
        return;
      
      
      // Create a map from SpecMeas to SpectrumType; although this is only used as a check
      //  in do_undo_or_redo.
      map<weak_ptr<SpecMeas>,SpecUtils::SpectrumType,std::owner_less<std::weak_ptr<SpecMeas>>> meas_to_type;
      const SpectrumType types[] = {
        SpectrumType::SecondForeground,
        SpectrumType::Background,
        SpectrumType::Foreground
      };
      for( const SpectrumType type : types )
      {
        const shared_ptr<SpecMeas> m = viewer->measurment(type);
        if( m && (meas_old_new_cal_map.count(m) 
                  || meas_old_new_peaks_map.count(m)
                  || meas_old_new_hint_peaks_map.count(m)) )
        {
          meas_to_type[m] = type;
        }
      }
      
      if( meas_to_type.empty() )
        return;
      
      // Lets avoid creating two copies of everything, and create a ptr to the undo/redo fcn
      auto doUndoOrRedo = make_shared<function<void(bool)>>(
        [meas_to_type, meas_old_new_peaks_map, 
         meas_old_new_hint_peaks_map, meas_old_new_cal_map]( const bool is_undo ){
          
          size_t num_peak_sets_used = 0;
          for( const auto &key : meas_old_new_cal_map )
          {
            const weak_ptr<SpecMeas> &specfile_weak = key.first;
            const meas_old_new_cal_t &meas_old_new_cal = key.second;
            
            meas_old_new_peaks_t meas_old_new_peaks;
            const auto peak_iter = meas_old_new_peaks_map.find(specfile_weak);
            if( peak_iter != end(meas_old_new_peaks_map) )
            {
              num_peak_sets_used += 1;
              meas_old_new_peaks = peak_iter->second;
            }
            
            meas_old_new_peaks_t meas_old_new_hint_peaks;
            const auto hint_peak_iter = meas_old_new_hint_peaks_map.find(specfile_weak);
            if( hint_peak_iter != end(meas_old_new_hint_peaks_map) )
            {
              //num_peak_sets_used += 1;
              meas_old_new_hint_peaks = hint_peak_iter->second;
            }
            
            SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
            const auto type_iter = meas_to_type.find(specfile_weak);
            assert( type_iter != end(meas_to_type) );
            if( type_iter != end(meas_to_type) )
              type = type_iter->second;
            
            do_undo_or_redo( is_undo, type, meas_old_new_peaks, 
                            meas_old_new_hint_peaks, meas_old_new_cal, specfile_weak );
          }//for( const meas_old_new_cal_t &meas_old_new_cal : meas_old_new_cal_map )
          
          assert( num_peak_sets_used == meas_old_new_peaks_map.size() );
      } );//create doUndoOrRedo function
      
      auto undo = [doUndoOrRedo](){ doUndoOrRedo->operator()( true ); };
      auto redo = [doUndoOrRedo](){ doUndoOrRedo->operator()( false ); };
      
      undoManager->addUndoRedoStep( undo, redo, "Edit energy cal" );
    }//~EnergyCalUndoRedoSentry()
  };//struct EnergyCalUndoRedoSentry
  
}//namespace


namespace EnergyCalImp
{

  class DeviationPairDisplay;
  class DevPair : public Wt::WContainerWidget
  {
  protected:
    DevPair( Wt::WContainerWidget *parent = 0 );
    void setDevPair( const std::pair<float,float> &d );
    std::pair<float,float> devPair() const;
    void visuallyIndicateChanged();
    //Wt::WDoubleSpinBox *m_energy, *m_offset;
    NativeFloatSpinBox *m_energy, *m_offset;
    Wt::WContainerWidget *m_delete;
    friend class DeviationPairDisplay;
  };//class DevPair
    
  class DeviationPairDisplay : public Wt::WContainerWidget
  {
  public:
    DeviationPairDisplay( Wt::WContainerWidget *parent = 0 );
    void setDeviationPairs( std::vector< std::pair<float,float> > d );
    std::vector< std::pair<float,float> > deviationPairs() const;
    void removeDevPair( DevPair *devpair );
    DevPair *newDevPair( const bool emitChangedNow );
    void setInvalidValues();
    void setValidValues();
    
    //void setMsg( const string &msg );
    
    enum class UserFieldChanged : int
    {
      AddedDeviationPair,
      RemovedDeviationPair,
      EnergyChanged,
      OffsetChanged
    };//enum UserFieldChanged
    
    /** Signal that gets emited when dev pair is added, deleted, or changed.
     Int argument is of type UserFieldChanged endum
     */
    Wt::Signal<int> &changed();
    
  protected:
    void sortDisplayOrder( const bool indicateVisually );
    
    void emitChanged( const UserFieldChanged whatChanged );
      
    Wt::Signal<int> m_changed;
    Wt::WContainerWidget *m_pairs;
    //Wt::WText *m_msg;
  };//class DeviationPairDisplay


class CALpDownloadResource : public Wt::WResource
{
  Wt::WApplication *m_app;
  InterSpec *m_interspec;
  EnergyCalTool *m_tool;
  
public:
  CALpDownloadResource( EnergyCalTool *tool, InterSpec *viewer, WObject* parent = nullptr )
  : WResource( parent ), m_app( WApplication::instance() ), m_interspec( viewer ), m_tool( tool )
  {
    assert( m_app );
    assert( m_tool );
    assert( m_interspec );
  }
  
  virtual ~CALpDownloadResource()
  {
    beingDeleted();
  }
  
  virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response )
  {
    assert( m_app );
    assert( m_interspec );
    
    try
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
        throw std::runtime_error( "Error grabbing application lock to from CALpDownloadResource resource." );
  
      const SpecUtils::SpectrumType type = m_tool->typeOfCurrentlyShowingCoefficients();
      shared_ptr<SpecMeas> meas = m_interspec->measurment( type );
      if( !meas )
        throw std::runtime_error( "Error getting spectrum file currently being shown." );
      
      string filename = meas->filename();
      if( filename.empty() )
        filename = "energy_calibration";
      const string orig_extension = SpecUtils::file_extension(filename);
      if( orig_extension.size() && (orig_extension.size() < filename.size()) )
        filename = filename.substr(0,filename.size() - orig_extension.size());
      filename += ".CALp";
      
      //Remove bad filename characters
      const string notallowed = "\\/:?\"<>|*";
      for( auto it = begin(filename) ; it < end(filename) ; ++it )
      {
        if( notallowed.find(*it) != string::npos )
          *it = ' ';
      }
      
      suggestFileName( filename, WResource::Attachment );
      response.setMimeType( "application/octet-stream" );
      
      // First loop over visible measurements, adding calibrations for each new detector name,
      //  then loop over all measurements to pick up the rest.  This is because there may be multiple
      //  calibrations for a single detector, but we want to prefer the ones in the currently
      //  displayed spectra, and not deal with the complexity of handling multiple calibrations per
      //  detector
      set<string> dets_so_far;
      
      //Note that disp_detectors has both gamma and neutron detectors, but we only care about gamma
      const vector<string> &gamma_detectors = meas->gamma_detector_names();
      const vector<string> &detectors = meas->detector_names();
      const vector<string> disp_detectors = m_interspec->detectorsToDisplay(type);
      const vector<string> &neut_dets = meas->neutron_detector_names();
      
      for( const int sample : m_interspec->displayedSamples(type) )
      {
        for( const string det : disp_detectors )
        {
          if( dets_so_far.count(det) )
            continue;
          
          const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det );
          shared_ptr<const SpecUtils::EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
          
          if( !cal || !cal->valid() || (cal->num_channels() < 3) )
          {
            // We'll assume that a detector that only has neutrons, will always only have neutrons...
            //  TODO: this may not actually be that good of an assumptions; re-evaluate later.
            if( std::find(begin(neut_dets), end(neut_dets), det) != end(neut_dets) )
              dets_so_far.insert( det );
            
            continue;
          }//if( energy cal is not valid )
          
          // Dont write the detector name if its unambiguous
          const string detname = (gamma_detectors.size() == 1) ? string() : det;
          
          if( SpecUtils::write_CALp_file(response.out(), cal, detname) )
          {
            dets_so_far.insert( det );
          }else
          {
            log("error") << "Error writing CALp file to WResource output stream.";
          }
        }//
        
        if( disp_detectors.size() == dets_so_far.size() )
          break;
      }//for( const int sample : m_interspec->displayedSamples(type) )
      
      if( disp_detectors.size() != dets_so_far.size() )
      {
        // Note that we are going over all samples and detectors here - being super thorough - we
        //  could probably tighten this up a lot.
        
        for( const int sample : meas->sample_numbers() )
        {
          for( const string &det : detectors )
          {
            if( dets_so_far.count(det) )
              continue;
            
            const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det );
            shared_ptr<const SpecUtils::EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
            
            if( !cal || !cal->valid() || (cal->num_channels() < 3) )
            {
              continue;
            }//if( energy cal is not valid )
            
            if( SpecUtils::write_CALp_file(response.out(), cal, det) )
            {
              dets_so_far.insert( det );
            }else
            {
              log("error") << "Error writing CALp file to WResource output stream (2).";
            }
          }//for( const string &det : gamma_dets )
        }//for( const int sample : meas->sample_numbers() )
      }//if( disp_detectors.size() != dets_so_far.size() )
    }catch( std::exception &e )
    {
      log("error") << "Error handling request for CalFileDownloadResource: " << e.what();
      response.out() << "Error creating CALp file: " << e.what()
                     << "\n\nPlease report to InterSpec@sandia.gov.";
      
      //passMessage( "Error getting spectrum file currently being shown", WarningWidget::WarningMsgHigh );
      
      response.setStatus(500);
      assert( 0 );
    }//try / catch
  }//void handleRequest(...)
};//class CalFileDownloadResource


DevPair::DevPair( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    //m_energy( new WDoubleSpinBox() ),
    //m_offset( new WDoubleSpinBox() ),
    m_energy( new NativeFloatSpinBox() ),
    m_offset( new NativeFloatSpinBox() ),
    m_delete( new WContainerWidget() )
{
  WGridLayout* layout = new WGridLayout();
  layout->setContentsMargins(0, 0, 0, 0);
  layout->setVerticalSpacing( 0 );
  
  setLayout(layout);
  setStyleClass( "DevPair" );
  
  layout->addWidget(m_energy,0,0);
  layout->addWidget(m_offset,0,1);
  layout->addWidget(m_delete,0,2);
  layout->setColumnStretch(0, 1);
  layout->setColumnStretch(1, 1);
  layout->setColumnStretch(2, 0);
          
  m_energy->setStyleClass( "DevEnergy" );
  m_offset->setStyleClass( "DevOffset" );
  m_energy->setPlaceholderText( "Energy" );
  m_offset->setPlaceholderText( "Offset" );
  m_energy->setValueText( "" );
  m_offset->setValueText( "" );
          
  m_delete->addStyleClass( "Wt-icon DeleteDevPair" );
}//DevPair constructor


void DevPair::setDevPair( const std::pair<float,float> &d )
{
  m_energy->setValue( d.first );
  m_offset->setValue( d.second );
  
  auto printval = []( float val ) -> std::string {
    char buffer[64];
    const float fraction = val - std::floor(val);
    if( fraction == 0.0 )
      snprintf( buffer, sizeof(buffer), "%.0f", val );
    else if( fabs(fraction - 0.1f) < 1.0E-4f )
      snprintf( buffer, sizeof(buffer), "%.1f", val );
    else
      snprintf( buffer, sizeof(buffer), "%.2f", val );
    return buffer;
  };
 
  m_energy->setText( printval(d.first) );
  m_offset->setText( printval(d.second) );
}


std::pair<float,float> DevPair::devPair() const
{
  const float energy = m_energy->value();
  const float offset = m_offset->value();
  return std::make_pair( energy, offset );
}


void DevPair::visuallyIndicateChanged()
{
  doJavaScript( "$('#" + id() + "').fadeIn(100).fadeOut(100).fadeIn(100).fadeOut(100).fadeIn(100);" );
}//void visuallyIndicateChanged();


DeviationPairDisplay::DeviationPairDisplay( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_pairs( NULL )
    //, m_msg( nullptr )
{
  addStyleClass( "DevPairDisplay" );
  WLabel *title = new WLabel( "Deviation Pairs", this );
  title->setStyleClass( "Wt-itemview Wt-header Wt-label DevPairTitle" );
  title->setInline( false );

  m_pairs = new WContainerWidget(this);
  m_pairs->setStyleClass( "DevPairsContainer" );
          
  //m_msg = new WText( "&nbsp;", this );
  //m_msg->addStyleClass( "DevPairMsg" );
  //m_msg->setHidden( true );
  
  auto footer = new WContainerWidget(this);
  footer->setStyleClass( "DevPairsFooter" );
  
  auto addBtn = new WContainerWidget( footer );
  addBtn->addStyleClass( "Wt-icon AddDevPair" );
  addBtn->clicked().connect( boost::bind(&DeviationPairDisplay::newDevPair, this, true) );
  addBtn->setToolTip( "Add another deviation pair" );
}//DeviationPairDisplay constructor

/*
void DeviationPairDisplay::setMsg( const string &msg )
{
  if( msg == m_msg->text().toUTF8() )
    return;
  
  const bool wasHidden = m_msg->isHidden();
  const bool shouldBeHidden = msg.empty();
  if( wasHidden != shouldBeHidden )
  {
    //WAnimation anim (WAnimation::AnimationEffect::Fade | WAnimation::AnimationEffect::SlideInFromBottom, WAnimation::Linear, 200);
    m_msg->setHidden( shouldBeHidden, anim );
  }
  m_msg->setText( WString::fromUTF8(msg) );
}
 */

void DeviationPairDisplay::setDeviationPairs( vector< pair<float,float> > d )
{
  m_pairs->clear();
  std::sort( d.begin(), d.end() );
  
  for( size_t i = 0; i < d.size(); ++i )
  {
    DevPair *dev = newDevPair( false );
    dev->setDevPair( d[i] );
  }//for( size_t i = 0; i < d.size(); ++i )
  
  sortDisplayOrder(false);
}//setDeviationPairs(...)


void DeviationPairDisplay::sortDisplayOrder( const bool indicateVisually )
{
  vector<size_t> sort_indices;
  vector<DevPair *> displays;
  
  vector<float> offsets;
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    DevPair *p = dynamic_cast<DevPair *>( t );
    if( p )
    {
      const size_t index = sort_indices.size();
      offsets.push_back( p->devPair().first );
      displays.push_back( p );
      sort_indices.push_back( index );
    }
  }//for( WWidget *t : childs )
  
  std::stable_sort( sort_indices.begin(), sort_indices.end(),
                    index_compare_assend<vector<float>&>(offsets) );

  bool order_changed = false;
  for( size_t i = 0; i < sort_indices.size(); ++i )
    order_changed |= (sort_indices[i] != i );
  
  if( !order_changed )
    return;
          
  for( size_t i = 0; i < displays.size(); ++i )
    m_pairs->removeWidget( displays[i] );
    
  for( size_t i = 0; i < displays.size(); ++i )
  {
    DevPair *p = displays[ sort_indices[i] ];
    m_pairs->addWidget( p );
    if( indicateVisually && (i != sort_indices[i]) )
      p->visuallyIndicateChanged();
  }
}//void sortDisplayOrder()


void DeviationPairDisplay::emitChanged( const UserFieldChanged whatChanged )
{
  m_changed.emit( static_cast<int>(whatChanged) );
}


vector< pair<float,float> > DeviationPairDisplay::deviationPairs() const
{
  vector< pair<float,float> > answer;
  const vector<WWidget *> childs = m_pairs->children();
  for( const WWidget *t : childs )
  {
    const DevPair *p = dynamic_cast<const DevPair *>( t );
    
    if( !p )
      continue;
    
    // We will only insert a deviation pair if both fields have text entered.
    if( p->m_energy->text().empty() || p->m_offset->text().empty() )
      continue;
    
    answer.push_back( p->devPair() );
  }//for( WWidget *t : childs )
  
  std::sort( answer.begin(), answer.end() );
  
  // The EnergyCalibration code implicitly will add in a {0,0} deviation pair if its needed, so
  //  might as well just add it in now so it isnt implict to the user
  //if( (answer.size() == 1) && (answer[0].first > 0.01f) )
  //  answer.insert( begin(answer), {0.0f,0.0f} );
  
  // Remove duplicates, usually {0,0}
  //const double epsilon = 1.0;
  //for( size_t i = 1; i < answer.size(); ++i )
  //{
  //  const auto &prev = answer[i-1];
  //  if( (fabs(answer[i].first - prev.first) < epsilon)
  //      && (fabs(answer[i].second - prev.second) < epsilon) )
  //    answer.erase( begin(answer) + i );
  //}
  
  return answer;
}//deviationPairs()


void DeviationPairDisplay::removeDevPair( DevPair *devpair )
{
  if( !devpair )
    return;
  
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    auto tt = dynamic_cast<DevPair *>( t ); //dynamic_cast prob not necessary
    if( devpair == tt )
    {
      delete devpair;
      sortDisplayOrder(false);
      emitChanged( UserFieldChanged::RemovedDeviationPair );
      return;
    }
  }//for( WWidget *t : childs )
}//removeDevPair(...)


void DeviationPairDisplay::setInvalidValues()
{
  addStyleClass( "InvalidDevPairs" );
}


void DeviationPairDisplay::setValidValues()
{
  removeStyleClass( "InvalidDevPairs" );
}



DevPair *DeviationPairDisplay::newDevPair( const bool emitChangedNow )
{
  DevPair *dev = new DevPair( m_pairs );
  dev->m_delete->clicked().connect( boost::bind( &DeviationPairDisplay::removeDevPair, this, dev ) );
  //dev->m_energy->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  //dev->m_offset->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_energy->valueChanged().connect( boost::bind( &DeviationPairDisplay::emitChanged, this, UserFieldChanged::EnergyChanged ) );
  dev->m_offset->valueChanged().connect( boost::bind( &DeviationPairDisplay::emitChanged, this, UserFieldChanged::OffsetChanged ) );
  dev->m_energy->blurred().connect( boost::bind(&DeviationPairDisplay::sortDisplayOrder, this, true) );
  
  if( emitChangedNow )
    emitChanged( UserFieldChanged::AddedDeviationPair );
  
  return dev;
}//newDevPair()

Wt::Signal<int> &DeviationPairDisplay::changed()
{
  return m_changed;
}


class CoefDisplay : public WContainerWidget
{
public:
  const size_t m_order;
  WLabel *m_label;
  WCheckBox *m_fit;
  NativeFloatSpinBox *m_value;
  
  CoefDisplay( const size_t order, WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
    m_order( order ),
    m_label( nullptr ),
    m_fit( nullptr ),
    m_value( nullptr )
  {
    addStyleClass( "CoefDisplay" );
    
    switch( order )
    {
      case 0:  m_label = new WLabel( "Offset", this );    break;
      case 1:  m_label = new WLabel( "Linear", this );    break;
      case 2:  m_label = new WLabel( "Quadratic", this ); break;
      case 3:  m_label = new WLabel( "Cubic", this );   break;
      default: m_label = new WLabel( std::to_string(order) + "'th order", this ); break;
    }//switch( order )
    
    m_label->addStyleClass( "CoefLabel" );
    
    m_value = new NativeFloatSpinBox( this );
    m_value->setSpinnerHidden( true );
    
    m_fit = new WCheckBox( "Fit", this );
    m_fit->addStyleClass( "CoefFit" );
  }//CoefDisplay
};//class CoefDisplay

class CalDisplay : public WContainerWidget
{
  static const size_t sm_min_coef_display ;
  
  EnergyCalTool *m_tool;
  const SpecUtils::SpectrumType m_cal_type;
  const std::string m_det_name;
  
  WText *m_type;
  WText *m_convertMsg;
  WContainerWidget *m_coefficients;
  DeviationPairDisplay *m_devPairs;
  
  //I cant decide if I like hiding deviation pair widget when there is none, or not
#define HIDE_EMPTY_DEV_PAIRS 0
#if( HIDE_EMPTY_DEV_PAIRS )
  WPushButton *m_addPairs;
#endif
  
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  WPushButton *m_fitCoeffs;
#endif
  
#if( IMP_CALp_BTN_NEAR_COEFS )
  WPushButton *m_downloadCALp;
  WPushButton *m_uploadCALp;
#endif
  
  shared_ptr<const SpecUtils::EnergyCalibration> m_cal;
  
public:
  CalDisplay( EnergyCalTool *tool,
             const SpecUtils::SpectrumType type,
             const std::string &detname,
             const bool isWideLayout,
             WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
   m_tool( tool ),
   m_cal_type( type ),
   m_det_name( detname ),
   m_type( nullptr ),
   m_convertMsg( nullptr ),
   m_coefficients( nullptr ),
   m_devPairs( nullptr )
#if( HIDE_EMPTY_DEV_PAIRS )
   , m_addPairs( nullptr )
#endif
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  , m_fitCoeffs( nullptr )
#endif
#if( IMP_CALp_BTN_NEAR_COEFS )
  , m_downloadCALp( nullptr )
  , m_uploadCALp( nullptr )
#endif
  {
    addStyleClass( "CalDisplay" );
    
    WGridLayout *layout = new WGridLayout( this );
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    
    WContainerWidget *coefDiv = new WContainerWidget();
    coefDiv->addStyleClass( "CoefCol" );
    layout->addWidget( coefDiv, 0, 0 );
    
    m_type = new WText( "&nbsp;", coefDiv );
    m_type->setInline( false );
    m_type->addStyleClass( "Wt-itemview Wt-header Wt-label CalType" );
    
    m_coefficients = new WContainerWidget( coefDiv );
    m_coefficients->addStyleClass( "CoefContent" );
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    WContainerWidget *btndiv = new WContainerWidget();
    btndiv->addStyleClass( "CalCoefsBtnDiv" );
  
#if( IMP_CALp_BTN_NEAR_COEFS )
    
    WResource *csv = m_model->peakCsvResource();
#if( BUILD_AS_OSX_APP || IOS )
    m_downloadCALp = new WAnchor( WLink(m_tool->calpResources()), btndiv );
    m_downloadCALp->setTarget( AnchorTarget::TargetNewWindow );
    m_downloadCALp->setStyleClass( "LinkBtn DownloadLink" );
#else
    m_downloadCALp = new WPushButton( btndiv );
    m_downloadCALp->setIcon( "InterSpec_resources/images/download_small.svg" );
    m_downloadCALp->setLink( WLink( m_tool->calpResources() ) );
    m_downloadCALp->setLinkTarget( Wt::TargetNewWindow );
    m_downloadCALp->setStyleClass( "LinkBtn DownloadBtn CALp" );
    
#if( ANDROID )
    // Using hacked saving to temporary file in Android, instead of via network download of file.
    m_downloadCALp->clicked().connect( std::bind([this](){
      android_download_workaround(m_tool->calpResources(), "energy_cal.CALp");
    }) );
#endif //ANDROID
    
#endif //#if( BUILD_AS_OSX_APP || IOS ) / #else

    m_downloadCALp->setText( "CALp" );
    
    m_uploadCALp = new WPushButton( btndiv );
    m_uploadCALp->setIcon( "InterSpec_resources/images/upload_small.svg" );
    m_uploadCALp->setStyleClass( "LinkBtn UploadBtn CALp" );
    m_uploadCALp->clicked().connect( m_tool, &EnergyCalTool::handleRequestToUploadCALp );
#endif //#if( IMP_CALp_BTN_NEAR_COEFS )
    
    WContainerWidget *spacer = new WContainerWidget( btndiv );
    spacer->addStyleClass( "Spacer" );
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
    m_fitCoeffs = new WPushButton( "Fit Coeffs", btndiv );
    m_fitCoeffs->addStyleClass( "CalCoefFitBtn" );
#endif
    
    layout->addWidget( btndiv, 1, 0 );
    layout->setRowStretch( 0, 1 );
#endif //#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    
    m_devPairs = new DeviationPairDisplay();
    
#if( HIDE_EMPTY_DEV_PAIRS )
    //For files with multiple detectors, the "Add dev. pairs" buttons doesnt show up right
    // for the detectors not currently showing - I guess should toggle dev pairs for all detectors.
    m_devPairs->setHidden( true );
    m_addPairs = new WPushButton( "Add dev. pairs" );
    m_addPairs->addStyleClass( "LinkBtn" );
    //m_addPairs->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
    m_addPairs->setHidden( true );
    m_addPairs->clicked().connect( this, &CalDisplay::showDevPairs );
#endif
    
    if( isWideLayout )
    {
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
      layout->addWidget( m_devPairs, 0, 1, 2, 1 );
#else
      layout->addWidget( m_devPairs, 0, 1 );
#endif
      
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( m_addPairs, 1, 0, AlignmentFlag::AlignRight );
      layout->setRowStretch( 0, 1 );
#endif
    }else
    {
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( m_addPairs, layout->rowCount(), 0, AlignmentFlag::AlignCenter );
      layout->addWidget( m_devPairs, layout->rowCount(), 0 );
#else
      layout->addWidget( m_devPairs, layout->rowCount(), 0 );
#endif
      m_devPairs->setHeight( 100 );
    }
    
    m_devPairs->changed().connect( boost::bind( &EnergyCalTool::userChangedDeviationPair, m_tool, this,
                                               boost::placeholders::_1 ) );
  }//CalDisplay( constructor )
  
  SpecUtils::SpectrumType spectrumType() const { return m_cal_type; }
  const std::string &detectorName() const { return m_det_name; }
  
    
  //void setDeviationPairMsg( const std::string &msg )
  //{
  //  m_devPairs->setMsg( msg );
  //}//void setDeviationPairMsg( const std::string &msg )
  
  
  shared_ptr<const SpecUtils::EnergyCalibration> lastSetCalibration()
  {
    return m_cal;
  }
  
  /// @param fitfor The order coefficients that should be set checked to fit for
  void setFitFor( const set<size_t> &fitfor )
  {
    for( auto w : m_coefficients->children() )
    {
      auto ww = dynamic_cast<const CoefDisplay *>( w );
      assert( ww );
      if( ww )
        ww->m_fit->setChecked( fitfor.count(ww->m_order) );
    }//for( auto w : m_coefficients->children() )
  }//void setFitFor( const set<size_t> &fitfor )
  
#if( HIDE_EMPTY_DEV_PAIRS )
  void showDevPairs()
  {
    m_addPairs->setHidden( true );
    m_devPairs->setHidden( false, WAnimation(WAnimation::Fade, WAnimation::Linear, 200) );
  }
#endif
  
  
  /// @returns The order coefficents that are checked to be fit for
  set<size_t> fitForCoefficents() const
  {
    set<size_t> coeffs;
    
    for( auto w : m_coefficients->children() )
    {
      auto ww = dynamic_cast<const CoefDisplay *>( w );
      assert( ww );
      if( !ww )
        continue;
      
      if( ww->m_fit->isChecked() )
        coeffs.insert( ww->m_order );
    }//for( auto w : m_coefficients->children() )
    
    return coeffs;
  }//vector<bool> fitForCoefficents() const
  
  
  vector<float> displayedCoefficents()
  {
    vector<float> coeffs;
  
    const auto existing = m_coefficients->children();
    for( size_t i = 0; i < existing.size(); ++i )
    {
      if( !existing[i] )
        continue;  //shouldnt happen, but JIC
      
      auto ww = dynamic_cast<const CoefDisplay *>( existing[i] );
      assert( ww );
      if( !ww )
        continue;
    
      //NativeFloatSpinBox avoids round-off errors, and we will rely on this here.
      const float dispvalue = ww->m_value->value();
      coeffs.push_back( dispvalue );
    }//for( size_t i = 0; i < existing.size(); ++i )
    
    //remove trailing zeros
    while( !coeffs.empty() && coeffs.back()==0.0f )
      coeffs.resize( coeffs.size() - 1 );
    
    return coeffs;
  }//vector<float> displayedCoefficents() const
  
  
  std::vector< std::pair<float,float> > displayedDeviationPairs() const
  {
    return m_devPairs->deviationPairs();
  }
  
  void setDeviationPairsInvalid()
  {
    m_devPairs->setInvalidValues();
  }
  
  void setDeviationPairsValid()
  {
    m_devPairs->setValidValues();
  }
  
  void updateToGui( const shared_ptr<const SpecUtils::EnergyCalibration> &cal )
  {
#if( HIDE_EMPTY_DEV_PAIRS )
    const bool hadCal = !!m_cal;
#endif
    
    m_cal = cal;
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    bool fitCoeffsVisible = false;
    const auto type = m_cal ? m_cal->type() : SpecUtils::EnergyCalType::InvalidEquationType;
    
    switch( type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        fitCoeffsVisible = true;
        break;
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        fitCoeffsVisible = false;
        break;
    }//switch( m_cal->type() )
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
    if( m_fitCoeffs )
      m_fitCoeffs->setHidden( !fitCoeffsVisible );
#endif
    
#if( IMP_CALp_BTN_NEAR_COEFS )
    if( m_downloadCALp )
      m_downloadCALp->setHidden( (type == SpecUtils::EnergyCalType::InvalidEquationType) );
#endif
#endif //IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS
    
    
    if( !m_cal )
    {
      m_type->setText( "No Calibration" );
      m_coefficients->clear();
      m_devPairs->setDeviationPairs( {} );
#if( HIDE_EMPTY_DEV_PAIRS )
      m_devPairs->setHidden( true );
      m_addPairs->setHidden( true );
#endif
      return;
    }//if( !m_cal )
    
    const char *typetxt = "";
    switch( m_cal->type() )
    {
      case SpecUtils::EnergyCalType::LowerChannelEdge:    typetxt = "Lower Channel Energy"; break;
      case SpecUtils::EnergyCalType::InvalidEquationType: typetxt = "Not Defined";          break;
      case SpecUtils::EnergyCalType::Polynomial:          typetxt = "Polynomial";           break;
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                                                          typetxt = "Default Polynomial";   break;
      case SpecUtils::EnergyCalType::FullRangeFraction:   typetxt = "Full Range Fraction";  break;
    }//switch( m_cal->type() )
    
    m_type->setText( typetxt );
    
    switch( m_cal->type() )
    {
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        m_coefficients->clear();
        m_devPairs->setDeviationPairs( {} );
        if( !m_devPairs->isHidden() )
          m_devPairs->hide();

        if( !m_convertMsg )
        {
          if( auto p = dynamic_cast<WContainerWidget *>( m_type->parent() ) )
          {
            m_convertMsg = new WText( "Please convert to Polynomial calibration to edit", p );
            m_convertMsg->addStyleClass( "ConvertToPolyMsg" );
            m_convertMsg->setInline( false );
          }
        }//if( !m_convertMsg )
        
#if( HIDE_EMPTY_DEV_PAIRS )
        m_addPairs->setHidden( true );
        m_devPairs->setHidden( true );
#endif
        return;
              
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
#if( !HIDE_EMPTY_DEV_PAIRS )
        if( m_devPairs->isHidden() )
          m_devPairs->show();
#endif
        if( m_convertMsg )
        {
          delete m_convertMsg;
          m_convertMsg = nullptr;
        }
        break;
    }//switch( m_cal->type() )
    
    const auto &devpairs = m_cal->deviation_pairs();
    const vector<float> &coeffs = m_cal->coefficients();
    
#if( HIDE_EMPTY_DEV_PAIRS )
    // Once deviation pairs are showing, we will leave them showing, even if there are no deviation
    //  pairs; this is because if you click to add deviation pairs, dont add any, then adjust the
    //  the gain, the deviation pair display getting hidden, causing a layout update, is really
    //  jarring.
    const auto anim = hadCal ? WAnimation(WAnimation::Fade, WAnimation::Linear, 200) : WAnimation{};
    if( m_devPairs->isHidden() && !devpairs.empty() )
      m_devPairs->setHidden( false, anim );
    m_addPairs->setHidden( !m_devPairs->isHidden(), anim );
#endif
    
    m_devPairs->setDeviationPairs( devpairs );
    //m_devPairs->changed().connect( std::bind( &EnergyCalTool::userChangedDeviationPair, m_tool, this) );
    
    const size_t num_coef_disp = std::max( coeffs.size(), sm_min_coef_display );
    vector<CoefDisplay *> coef_disps( num_coef_disp, nullptr );
    
    size_t coefnum = 0;
    const auto existing = m_coefficients->children();
    for( size_t i = 0; i < existing.size(); ++i )
    {
      auto ww = dynamic_cast<CoefDisplay *>( existing[i] );
      assert( ww || !existing[i] );
      
      if( ww )
      {
        if( coefnum >= num_coef_disp )
        {
          m_coefficients->removeWidget( existing[i] );
          delete existing[i];
        }else
        {
          assert( coefnum < coef_disps.size() );
          coef_disps[coefnum] = ww;
          const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
          ww->m_value->setValue( value );
        }
        
        ++coefnum;
      }//if( ww )
    }//for( size_t i = 0; i < existing.size(); ++i )
    
    for( ; coefnum < num_coef_disp; ++coefnum )
    {
      CoefDisplay *disp = new CoefDisplay( coefnum, m_coefficients );
      coef_disps[coefnum] = disp;
      const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
      disp->m_value->setValue( value );
      disp->m_fit->setChecked( (coefnum < 2) );
      
      disp->m_fit->changed().connect( m_tool, &EnergyCalTool::updateFitButtonStatus );
      
      /* Note: if the user uses the up.down arrows in a NativeFloatSpinBox to change values, things
               get all messed up (new values get set via c++ messing  up current values, or the
               valueChanged() callback gets called like 10 times per second, causing changes faster
               than everything can keep up, and just generally poor working), so for now I have
               disabled these spinners via #NativeFloatSpinBox::setSpinnerHidden()
       */
      disp->m_value->valueChanged().connect( boost::bind(&EnergyCalTool::userChangedCoefficient, m_tool, coefnum, this) );
    }
    
    
    //Set the step size to move the upper range of energy by about 1 keV per step
    // Set up the little tick/spin/whatever boxes
    /*
     //The NativeFloatSpinBox::setSpinnerHidden() call is currently removing the spin-box up/down
     //  arrow, so we wont set the step-size, as on Firefox if we fit for a value, then it will turn
     //  red if the new value doesnt hit on the step size
    for( size_t i = 0; i < coef_disps.size(); ++i )
    {
      CoefDisplay *disp = coef_disps[i];
      assert( disp );
      
      float stepsize = 1.0f;
      if( m_cal && m_cal->num_channels() > 4 )
      {
        const size_t nchannel = m_cal->num_channels();
        switch( m_cal->type() )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            stepsize = 1.0f / std::pow(nchannel,i);
            break;
            
          case SpecUtils::EnergyCalType::FullRangeFraction:
            stepsize = 1.0;
            break;
            
          case SpecUtils::EnergyCalType::InvalidEquationType:
          case SpecUtils::EnergyCalType::LowerChannelEdge:
            stepsize = 0.0f;
            break;
        }//switch( m_coeffEquationType )
      }//if( valid calibration )
      
      disp->m_value->setSingleStep( stepsize );
    }//for( int i = 0; i < sm_numCoefs; ++i )
     */
    
  }//updateToGui(...)

#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  void setFitButtonEnabled( const bool canFitCeofs )
  {
    if( !m_cal )
      return;
      
    if( m_fitCoeffs->isEnabled() != canFitCeofs )
      m_fitCoeffs->setEnabled( canFitCeofs );
  }
  
  Wt::EventSignal<Wt::WMouseEvent> &doFitCoeffs()
  {
    return m_fitCoeffs->clicked();
  }
#endif
  
};//class CalDisplay

const size_t CalDisplay::sm_min_coef_display = 4;


}//namespace


EnergyCalTool::EnergyCalTool( InterSpec *viewer, PeakModel *peakModel, WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_peakModel( peakModel ),
  m_calpResource( new EnergyCalImp::CALpDownloadResource(this, viewer, this) ),
  m_tallLayoutContent( nullptr ),
  m_peakTable( nullptr ),
  m_specTypeMenu( nullptr ),
  m_specTypeMenuStack( nullptr ),
  m_detectorMenu{ nullptr },
  m_calInfoDisplayStack( nullptr ),
  m_noCalTxt( nullptr ),
  m_moreActionsColumn( nullptr ),
  m_applyToColumn( nullptr ),
  m_detColumn( nullptr ),
  m_detColLayout( nullptr ),
  m_calColumn( nullptr ),
  m_peakTableColumn( nullptr ),
  m_layout( nullptr ),
  m_applyToCbs{ nullptr },
  m_moreActions{ nullptr },
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  m_fitCalBtn( nullptr ),
#endif
#if( !IMP_CALp_BTN_NEAR_COEFS )
  m_downloadCALp( nullptr ),
  m_uploadCALp( nullptr ),
#endif
  m_lastGraphicalRecal( 0 ),
  m_lastGraphicalRecalType( EnergyCalGraphicalConfirm::NumRecalTypes ),
  m_lastGraphicalRecalEnergy( -999.0f ),
  m_graphicalRecal( nullptr ),
  m_addActionWindow( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/EnergyCalTool.css" );
  
  assert( viewer );
  viewer->useMessageResourceBundle( "EnergyCalTool" );
    
  addStyleClass( "EnergyCalTool" );
  
  initWidgets( EnergyCalTool::LayoutType::Wide );
}


void EnergyCalTool::initWidgets( EnergyCalTool::LayoutType layoutType )
{
  const bool wide = (layoutType == LayoutType::Wide);
  
  if( (wide && !m_tallLayoutContent && m_layout) || (!wide && m_tallLayoutContent) )
    return;
  
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal );
    m_graphicalRecal = nullptr;
  }//if( m_graphicalRecal )
  
  
  
  if( wide )
  {
    removeStyleClass( "TallEnergyCal" );
    if( m_tallLayoutContent )
      delete m_tallLayoutContent;
    m_tallLayoutContent = nullptr;
    
    m_layout = new WGridLayout( this );
  }else
  {
    addStyleClass( "TallEnergyCal" );
    if( m_layout )
      delete m_layout;
    
    m_tallLayoutContent = new WContainerWidget( this );
    m_layout = new WGridLayout( m_tallLayoutContent );
  }//if( wide ) / else
  
  
  // \TODO: null out all the other memener variables
  m_peakTable = nullptr;
  m_specTypeMenu = nullptr;
  m_specTypeMenuStack = nullptr;
  for( auto &m : m_detectorMenu )
    m = nullptr;
  m_calInfoDisplayStack = nullptr;
  m_noCalTxt = nullptr;
  m_moreActionsColumn = nullptr;
  m_applyToColumn = nullptr;
  m_detColumn = nullptr;
  m_detColLayout = nullptr;
  m_calColumn = nullptr;
  m_peakTableColumn = nullptr;
  for( auto &m : m_applyToCbs )
    m = nullptr;
  for( auto &m : m_moreActions )
    m = nullptr;
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  m_fitCalBtn = nullptr;
#endif
#if( !IMP_CALp_BTN_NEAR_COEFS )
  m_downloadCALp = nullptr;
  m_uploadCALp = nullptr;
#endif
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setVerticalSpacing( 0 );
  m_layout->setHorizontalSpacing( 0 );
  
  m_noCalTxt = new WText( WString::tr("ect-no-spec") );
  m_noCalTxt->addStyleClass( "NoCalContentTxt" );
  if( wide )
    m_layout->addWidget( m_noCalTxt, 0, 0, AlignmentFlag::AlignCenter | AlignmentFlag::AlignMiddle );
  else
    m_layout->addWidget( m_noCalTxt, 0, 0, 1, 2, AlignmentFlag::AlignCenter | AlignmentFlag::AlignMiddle );
  
  //Create the more actions column...
  m_moreActionsColumn = new WContainerWidget();
  m_moreActionsColumn->addStyleClass( "ToolTabTitledColumn MoreActionCol" );
  if( wide )
    m_layout->addWidget( m_moreActionsColumn, 0, 1 );
  else
    m_layout->addWidget( m_moreActionsColumn, 3, 0 );
  
  WGridLayout *collayout = new WGridLayout( m_moreActionsColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  WText *header = new WText( WString::tr("ect-more-act") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *moreActionsDiv = new WContainerWidget();
  moreActionsDiv->addStyleClass( "ToolTabTitledColumnContent MoreActionsMenuContent" );
  collayout->addWidget( moreActionsDiv, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  WContainerWidget *moreActionsList = new WContainerWidget( moreActionsDiv );
  moreActionsList->addStyleClass( "MoreActionsMenuList" );
  moreActionsList->setList( true );
  
  
  for( MoreActionsIndex index = static_cast<MoreActionsIndex>(0);
      index < MoreActionsIndex::NumMoreActionsIndex;
      index = MoreActionsIndex(static_cast<int>(index) + 1) )
  {
    const char *label = "", *tooltip = nullptr;
    switch( index )
    {
      case MoreActionsIndex::Linearize:
        label = "ect-linearize";
        tooltip = "ect-tt-linearize";
        break;
        
      case MoreActionsIndex::Truncate:
        label = "ect-truncate";
        tooltip = "ect-tt-truncate";
        break;
        
      case MoreActionsIndex::CombineChannels:
        label = "ect-combine";
        tooltip = "ect-tt-combine";
        break;
        
      case MoreActionsIndex::ConvertToFrf:
        label = "ect-to-frf";
        tooltip = "ect-tt-to-frf";
        break;
        
      case MoreActionsIndex::ConvertToPoly:
        label = "ect-to-poly";
        tooltip = "ect-tt-to-poly";
        break;
        
      case MoreActionsIndex::MultipleFilesCal:
        label = "ect-multi-file";
        tooltip = "ect-tt-multi-file";
        break;
        
      case MoreActionsIndex::NumMoreActionsIndex:
        assert(0);
        break;
    }//switch( index )
    
    WContainerWidget *holder = new WContainerWidget( moreActionsList );
    m_moreActions[static_cast<int>(index)] = new WAnchor( WLink(), WString::tr(label), holder );
    m_moreActions[static_cast<int>(index)]->clicked().connect( boost::bind(&EnergyCalTool::moreActionBtnClicked, this, index) );
    
    assert( tooltip );
    if( tooltip )
      HelpSystem::attachToolTipOn( holder, WString::tr(tooltip), showToolTips );
  }//for( loop over more actions )
  
  WContainerWidget *btndiv = new WContainerWidget();
  btndiv->addStyleClass( "BtmBtnDiv" );
  collayout->addWidget( btndiv, 2, 0 );
  
  auto helpBtn = new WContainerWidget( btndiv );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "energy-calibration" ) );

  
#if( !IMP_CALp_BTN_NEAR_COEFS )
  m_uploadCALp = new WPushButton( btndiv );
  m_uploadCALp->setIcon( "InterSpec_resources/images/upload_small.svg" );
  m_uploadCALp->setStyleClass( "LinkBtn UploadBtn CALp" );
  m_uploadCALp->clicked().connect( this, &EnergyCalTool::handleRequestToUploadCALp );
  
#if( BUILD_AS_OSX_APP || IOS )
  m_downloadCALp = new WAnchor( WLink(m_calpResource), btndiv );
  m_downloadCALp->setTarget( AnchorTarget::TargetNewWindow );
  m_downloadCALp->setStyleClass( "LinkBtn DownloadLink CALp" );
#else
  m_downloadCALp = new WPushButton( btndiv );
  m_downloadCALp->setIcon( "InterSpec_resources/images/download_small.svg" );
  m_downloadCALp->setLink( WLink( m_calpResource ) );
  m_downloadCALp->setLinkTarget( Wt::TargetNewWindow );
  m_downloadCALp->setStyleClass( "LinkBtn DownloadBtn CALp" );
  
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  m_downloadCALp->clicked().connect( std::bind([this](){
    android_download_workaround( m_calpResource, "energy_cal.CALp");
  }) );
#endif //ANDROID
  
#endif
  m_downloadCALp->setText( "CALp" );
  
  m_downloadCALp->clicked().connect( std::bind([this](){
    m_interspec->logMessage( WString::tr("ect-export-CALp-msg"), WarningWidget::WarningMsgInfo );
  }) );
  HelpSystem::attachToolTipOn( m_downloadCALp, WString::tr("ect-tt-CALp"), showToolTips );
  
  m_downloadCALp->setHidden( true );
  m_uploadCALp->setHidden( true );
#endif // !IMP_CALp_BTN_NEAR_COEFS
  
  
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  m_fitCalBtn = new WPushButton( WString::tr("ect-fit-coeff-btn"), btndiv );
  m_fitCalBtn->addStyleClass( "FitCoefBtn" );
  m_fitCalBtn->clicked().connect( this, &EnergyCalTool::fitCoefficients );
  m_fitCalBtn->setDisabled( true );
  HelpSystem::attachToolTipOn( m_fitCalBtn, WString::tr("ect-tt-fit-coeff-btn"), showToolTips );
#endif // !IMP_COEF_FIT_BTN_NEAR_COEFS
  
  // Create the "Apply To" column that determines what to apply changes to
  m_applyToColumn = new WContainerWidget();
  m_applyToColumn->addStyleClass( "ToolTabTitledColumn ApplyToCol" );
  if( wide )
    m_layout->addWidget( m_applyToColumn, 0, 2 );
  else
    m_layout->addWidget( m_applyToColumn, 3, 1 );
  
  collayout = new WGridLayout( m_applyToColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( WString::tr("ect-apply-changes-to") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *applyToDiv = new WContainerWidget();
  applyToDiv->addStyleClass( "ToolTabTitledColumnContent ApplyToMenuContent" );
  collayout->addWidget( applyToDiv, 1, 0 );
  
  WContainerWidget *applyToList = new WContainerWidget( applyToDiv );
  applyToList->addStyleClass( "ApplyToMenuList" );
  applyToList->setList( true );
  
  
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    const char *label = "";
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:         label = "Foreground";       break;
      case ApplyToCbIndex::ApplyToBackground:         label = "Background";       break;
      case ApplyToCbIndex::ApplyToSecondary:          label = "Secondary";        break;
      case ApplyToCbIndex::ApplyToDisplayedDetectors: label = "ect-disp-dets";    break;
      case ApplyToCbIndex::ApplyToAllDetectors:       label = "ect-all-dets";     break;
      case ApplyToCbIndex::ApplyToDisplayedSamples:   label = "ect-disp-samples"; break;
      case ApplyToCbIndex::ApplyToAllSamples:         label = "ect-all-samples";  break;
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    WContainerWidget *item = new WContainerWidget( applyToList );
    item->addStyleClass( "ApplyToItem" );
    auto cb = new WCheckBox( WString::tr(label), item );
    cb->setWordWrap( false );
    cb->addStyleClass( "ApplyToItem" );
    cb->setInline( false );
    
    
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:
      case ApplyToCbIndex::ApplyToBackground:
      case ApplyToCbIndex::ApplyToSecondary:
      case ApplyToCbIndex::ApplyToAllDetectors:
      case ApplyToCbIndex::ApplyToAllSamples:
        cb->setChecked( true );
        break;
        
      case ApplyToCbIndex::ApplyToDisplayedDetectors:
      case ApplyToCbIndex::ApplyToDisplayedSamples:
      case ApplyToCbIndex::NumApplyToCbIndex:
        cb->setChecked( false );
        break;
    }//switch( index )
    
    cb->checked().connect( boost::bind( &EnergyCalTool::applyToCbChanged, this, index ) );
    cb->unChecked().connect( boost::bind( &EnergyCalTool::applyToCbChanged, this, index ) );
    
    m_applyToCbs[index] = cb;
  }//for( loop over ApplyToCbIndex )
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);

  
  // Create the "Coefficients" column that show the polynomial/FRF coefficents.
  m_calColumn = new WContainerWidget();
  m_calColumn->addStyleClass( "ToolTabTitledColumn CoefColumn" );
  if( wide )
    m_layout->addWidget( m_calColumn, 0, 3 );
  else
    m_layout->addWidget( m_calColumn, 1, 0, 1, 2 );
  
  
  collayout = new WGridLayout( m_calColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  if( wide )
    collayout->setColumnStretch( 1, 1 );
  
  header = new WText( WString::tr("ect-calib-coeffs") );
  header->addStyleClass( "ToolTabColumnTitle" );
  
  collayout->addWidget( header, 0, 0, 1, 2 );
  //collayout->addWidget( m_calInfoDisplayStack, 1, 0 );
  
  
  // Create the "Detector" column that determines which coefficients to show
  m_detColumn = new WContainerWidget();
  m_detColumn->addStyleClass( "DetCol" );
  collayout->addWidget( m_detColumn, 1, 0 );
  
  m_detColLayout = new WGridLayout( m_detColumn );
  m_detColLayout->setContentsMargins( 0, 0, 0, 0 );
  m_detColLayout->setVerticalSpacing( 0 );
  m_detColLayout->setHorizontalSpacing( 0 );
  
  auto detheader = new WText( WString::tr("Detector") );
  detheader->setInline( false );
  detheader->addStyleClass( "DetHdr Wt-itemview Wt-header Wt-label" );
  //detheader->resize( WLength::Auto, WLength(20,WLength::Unit::Pixel) );
  //collayout->addWidget( detheader, 0, 0 );
  m_detColLayout->addWidget( detheader, 0, 0  );
  m_detColLayout->setRowStretch( 2, 1 );
  
  // Create the "Cal Peaks" table
  m_peakTableColumn = new WContainerWidget();
  m_peakTableColumn->addStyleClass( "ToolTabTitledColumn PeakTableCol" );
  if( wide )
  {
    m_layout->addWidget( m_peakTableColumn, 0, 4 );
    m_layout->setColumnStretch( 4, 1 );
  }else
  {
    m_layout->addWidget( m_peakTableColumn, 2, 0, 1, 2 );
  }

  collayout = new WGridLayout( m_peakTableColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  if( wide )
    collayout->setRowStretch( 1, 1 );
  
  header = new WText( WString::tr("ect-cal-peaks") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  m_peakTable = new RowStretchTreeView();
  m_peakTable->addStyleClass( "ToolTabTitledColumnContent PeakTable" );
  collayout->addWidget( m_peakTable, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  m_peakTable->setRootIsDecorated( false ); //makes the tree look like a table! :)
  m_peakTable->setModel( m_peakModel );
  const int numModelCol = m_peakModel->columnCount();
  for( int col = 0; col < numModelCol; ++col )
    m_peakTable->setColumnHidden( col, true );
  
  m_peakTable->setSortingEnabled( true );
  m_peakTable->setAlternatingRowColors( true );
  m_peakTable->setSelectable( true );
  m_peakTable->setSelectionMode( SingleSelection );
  m_peakTable->setEditTriggers( WAbstractItemView::SingleClicked
                               | WAbstractItemView::DoubleClicked );
  
  m_peakTable->setColumnHidden( PeakModel::kUseForCalibration, false );
  m_peakTable->setColumnHidden( PeakModel::kMean, false );
  m_peakTable->setColumnHidden( PeakModel::kIsotope, false );
  m_peakTable->setColumnHidden( PeakModel::kPhotoPeakEnergy, false );
  m_peakTable->setColumnHidden( PeakModel::kDifference, false );
  
  
  m_peakTable->setColumnWidth( PeakModel::kUseForCalibration, WLength(3.7, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kMean, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kIsotope, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(6.25, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kDifference, WLength(5, WLength::FontEm) );
  
  
  
  WItemDelegate *dblDelagate = new WItemDelegate( m_peakTable );
  dblDelagate->setTextFormat( "%.2f" );
  m_peakTable->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  
  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );
  
  PhotopeakDelegate *photopeakDelegate = new PhotopeakDelegate( PhotopeakDelegate::GammaEnergyDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );
  
  m_peakModel->dataChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsRemoved().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsInserted().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->layoutChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  
  m_interspec->displayedSpectrumChanged().connect(
              boost::bind( &EnergyCalTool::displayedSpectrumChanged,
                           this, boost::placeholders::_1, boost::placeholders::_2,
                          boost::placeholders::_3, boost::placeholders::_4 ) );
  
  m_renderFlags |= EnergyCalToolRenderFlags::FullGuiUpdate;
  scheduleRender();
}//void initWidgets( EnergyCalTool::LayoutType layout )


void EnergyCalTool::setWideLayout()
{
  initWidgets( EnergyCalTool::LayoutType::Wide );
}//void setWideLayout()


void EnergyCalTool::setTallLayout()
{
  initWidgets( EnergyCalTool::LayoutType::Tall );
}//void setTallLayout()


EnergyCalTool::~EnergyCalTool()
{
}
  

set<string> EnergyCalTool::gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type )
{
  //We want the names of just the detectors that have gamma calibration information, of the
  //  currently displayed samples.
  //  We will assume the first Measurement for a given named detector will have gamma data if
  //  any of the Measurements from that detector will.
  //  \TODO: evaluate if that is true.
  
  auto meas = m_interspec->measurment( type );
  
  if( !meas )
    return {};
  
  const vector<string> &detnames = meas->gamma_detector_names();
  const vector<string> displayedDets = m_interspec->detectorsToDisplay(type);
  const set<int> &samples = m_interspec->displayedSamples(type);
  
  set<string> detectors, nongammadets;
  
  for( const int sample : samples )
  {
    for( const string &name : detnames )
    {
      if( detectors.count(name) || nongammadets.count(name) )
        continue;
      
      auto m = meas->measurement( sample, name );
      if( m && (m->num_gamma_channels() > 4) )
        detectors.insert( name );
      else if( m && !detectors.count(name) )
        nongammadets.insert(name);
    }//for( const string &name : detnames )
    
    if( (detectors.size() + nongammadets.size()) == detnames.size() )
      break;
  }//for( const int sample : samples )
  
  return detectors;
}//set<string> gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type )


#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
vector<EnergyCalImp::CalDisplay *> EnergyCalTool::calDisplays()
{
  vector<EnergyCalImp::CalDisplay *> answer;
  for( int i = 0; i < 3; ++i )
  {
    WMenu *detMenu = m_detectorMenu[i];
    
    if( !detMenu )
      continue;
    
    WStackedWidget *stack = detMenu->contentsStack();
    
    if( !stack )
      continue;
    
    for( WWidget *w : stack->children() )
    {
      EnergyCalImp::CalDisplay *caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( w );
      assert( caldisp );
      if( caldisp )
        answer.push_back( caldisp );
    }
  }//for( int i = 0; i < 3; ++i )
  
  return answer;
}//vector<EnergyCalImp::CalDisplay *> calDisplays()
#endif





vector<MeasToApplyCoefChangeTo> EnergyCalTool::measurementsToApplyCoeffChangeTo()
{
  std::vector<MeasToApplyCoefChangeTo> answer;
  
  //Lets loop over spectrum types (Foreground, Background, Secondary), and decide if we should apply
  //  changes to that file, and if so, decide which sample numbers/detectors
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  for( const auto spectype : spectypes )
  {
    auto meas = m_interspec->measurment( spectype );
    if( !meas )
      continue;
    
    switch( spectype )
    {
      case SpecUtils::SpectrumType::Foreground:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked() )
          continue;
        break;
        
      case SpecUtils::SpectrumType::SecondForeground:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToSecondary]->isChecked() )
          continue;
        break;
        
      case SpecUtils::SpectrumType::Background:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToBackground]->isChecked() )
          continue;
        break;
    }//switch( spectype )
    
    //If we're here, we should apply changes to this spectrum type.
    
    //It could be that the background SpecFile is the same as the foreground, so check if we
    //  already have an entry for this file in answer, and if so, use it.  We dont want duplicate
    //  entries for SpecFiles since this could cause us to maybe move peaks multiple times or
    //  something
    //  \TODO: if the background and foreground use different detectors (no overlap) or different
    //         sample numbers (no overlap) should return multiple entries for the SpecFile to handle
    //         the edge-case correctly
    MeasToApplyCoefChangeTo *changes = nullptr;
    for( size_t i = 0; !changes && i < answer.size(); ++i )
      changes = (answer[i].meas == meas) ? &(answer[i]) : changes;
    
    if( !changes )
    {
      //We havent seen this SpecFile yet, create an entry in answer for it
      answer.emplace_back();
      changes = &(answer.back()); //C++14 returns a reference to the emplaced object.
      changes->meas = meas;
    }
    
    assert( changes );
    
    //m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors]->parent()->isHidden();
    const set<string> detectors = gammaDetectorsForDisplayedSamples(spectype);
    const vector<string> displayed_dets = m_interspec->detectorsToDisplay(spectype);
    
    bool displayingAllDets = true;
    for( const auto &det : detectors )
    {
      if( std::find(begin(displayed_dets), end(displayed_dets), det) == end(displayed_dets) )
        displayingAllDets = false;
    }
    
    if( displayingAllDets )
    {
      //We will insert detector names since different SpecType from the same file could have
      //  different detector names available.
      // \TODO: if InterSpec class is upgraded to select detector by SpecType, we will have to
      //        upgrade this part of the code.
      for( const auto &det : detectors )
        changes->detectors.insert( det );
    }else
    {
      const auto applyToAllCb = m_applyToCbs[ApplyToCbIndex::ApplyToAllDetectors];
      const auto applyToDisplayedCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors];
      
      const bool toAll = (applyToAllCb->parent()
                          && !applyToAllCb->parent()->isHidden()
                          && applyToAllCb->isChecked());
      const bool toDisplayed = (applyToDisplayedCb->parent()
                                && !applyToDisplayedCb->parent()->isHidden()
                                && applyToDisplayedCb->isChecked());
      if( toAll == toDisplayed )
      {
        cerr << "EnergyCalTool::measurementsToApplyCoeffChangeTo:"
                " got (toAll == toDisplayed) which shouldnt have happended" << endl;
      }
        
      if( toAll || (toAll == toDisplayed) )
      {
        for( const auto &det : detectors )
          changes->detectors.insert( det );
      }else
      {
        for( const auto &dispdet : displayed_dets )
        {
          if( detectors.count(dispdet) )
            changes->detectors.insert( dispdet );
        }
      }//if( apply to all detectors ) / else ( only displayed detectors )
    }//if( displayingAllDets ) / else
    
    
    const auto toAllSamplesCb = m_applyToCbs[ApplyToCbIndex::ApplyToAllSamples];
    const bool onlyDispSamples = (toAllSamplesCb->parent()
                                  && !toAllSamplesCb->parent()->isHidden()
                                  && !toAllSamplesCb->isChecked());
    
    if( onlyDispSamples )
    {
      const set<int> &displayed_samples = m_interspec->displayedSamples(spectype);
      for( const int sample : displayed_samples )
        changes->sample_numbers.insert( sample );
    }else
    {
      changes->sample_numbers = meas->sample_numbers();
    }
  }//for( const auto spectype : spectypes )
  
  
  return answer;
}//std::vector<MeasToApplyCoefChangeTo> measurementsToApplyCoeffChangeTo()


void EnergyCalTool::applyCALpEnergyCal( std::map<std::string,std::shared_ptr<const SpecUtils::EnergyCalibration>> det_to_cal,
                                         const SpecUtils::SpectrumType specfile,
                                         const bool all_detectors, const bool all_samples )
{
  // TODO: add option to not set deviation pairs
  EnergyCalUndoRedoSentry undo_sentry;
  
  set<string> fore_gamma_dets;
  const shared_ptr<SpecMeas> measurment = m_interspec->measurment( specfile );
  const shared_ptr<const SpecUtils::Measurement> disp_spec = m_interspec->displayedHistogram( specfile );
  const shared_ptr<const SpecUtils::EnergyCalibration> old_disp_cal
                                           = disp_spec ? disp_spec->energy_calibration() : nullptr;
  const set<int> &disp_samples = m_interspec->displayedSamples( specfile );
  const vector<string> disp_detectors = m_interspec->detectorsToDisplay( specfile );
  
  const char * const desc = SpecUtils::descriptionText(specfile);
  
  if( !measurment || !disp_spec || !old_disp_cal || !old_disp_cal->valid() )
    throw runtime_error( WString::tr("ect-CALp-no-meas").arg(desc).toUTF8() );
  
  MeasToApplyCoefChangeTo tochange;
  tochange.meas = measurment;
  if( all_samples )
    tochange.sample_numbers = measurment->sample_numbers();
  else
    tochange.sample_numbers = disp_samples;
  
  const vector<string> &gamma_det_names = measurment->gamma_detector_names();
  if( all_detectors )
  {
    tochange.detectors.insert( begin(gamma_det_names), end(gamma_det_names) );
  }else
  {
    for( const string &det : disp_detectors )
    {
      const auto pos = std::find( begin(gamma_det_names), end(gamma_det_names), det );
      if( pos != end(gamma_det_names) )
        tochange.detectors.insert( det );
    }
  }//if( all_detectors ) / else
  
  
  if( tochange.detectors.empty() )
    throw runtime_error( WString::tr("ect-CALp-not-applic-det").arg(desc).toUTF8() );
  
  if( tochange.sample_numbers.empty() )
    throw runtime_error( WString::tr("ect-CALp-not-applic-sample").arg(desc).toUTF8() );
  
  
  if( det_to_cal.size() == 1 )
  {
    const shared_ptr<const SpecUtils::EnergyCalibration> new_cal = det_to_cal.begin()->second;
    
    if( !new_cal || !new_cal->valid() )
      throw runtime_error( WString::tr("ect-CALp-invalid-input").toUTF8() );
    
    // Its possible things could work out if there is a different number of channels, but this
    //  doesnt make much sense, so we'll require it, at least for now.
    if( new_cal->num_channels() != old_disp_cal->num_channels() )
    {
      throw runtime_error( WString::tr("ect-CALp-num-channel-mismatch")
                          .arg( static_cast<int>(new_cal->num_channels()) )
                          .arg( static_cast<int>(old_disp_cal->num_channels()) )
                          .toUTF8() );
    }
    
    if( tochange.detectors.size() == 1 )
    {
      setEnergyCal( new_cal, tochange, true );
    }else
    {
      // Check to see if all detectors are using the same energy calibration, and if so, set the
      //  calibration, otherwise do the applyCalChange(...)
      set<shared_ptr<const SpecUtils::EnergyCalibration>> old_energy_cals;
      for( const string &det : tochange.detectors )
      {
        for( const int sample : tochange.sample_numbers )
        {
          const auto m = tochange.meas->measurement( sample, det );
          const shared_ptr<const SpecUtils::EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
          if( cal && cal->valid() && (cal->num_channels() >= 3) )
            old_energy_cals.insert( m->energy_calibration() );
        }//for( loop over samples )
      }//for( loop over detector names )
      
      if( old_energy_cals.size() > 1 )
      {
        // Note: we ARE NOT updating the deviation pairs here... I'm not totally sure how to
        //       make everything consistent from the users perspective... perhaps just ignoring
        //       this corner case of using a CALp file with a single calibration, with data that
        //       currently uses multiple calibrations is fine - after all all that is lost is the
        //       deviation pairs, which its not clear that you do actually want the deviation pairs.
        //       The work around from the user perspective would be to display a single detector at
        //       a time, and then apply the CALp file.  I doubt anyone will ever actually hit this
        //       corner case, and even if they did, the behaviour is probably what they want anyway.
        applyCalChange( old_disp_cal, new_cal, { tochange }, false );
      }else
      {
        setEnergyCal( new_cal, tochange, true );
      }
    }//if( we are changing a single detector ) / else
  }else
  {
    // We will apply the input calibration on a detector-by-detector basis, requiring the
    //  input to have info for every relevant detector.
    //
    //  We will only adjust peak positions if the detectors energy calibration is the
    //   `suggested_sum_energy_calibration`, and even then, only once.
    //   (note: not trusting disp_spec->energy_calibration() to be same as
    //          `suggested_sum_energy_calibration` - but I think we could).
    bool adjusted_peaks = false;
    shared_ptr<const SpecUtils::EnergyCalibration> display_cal;
    try
    {
      display_cal = tochange.meas->suggested_sum_energy_calibration( disp_samples, disp_detectors );
    }catch( std::exception & )
    {
    }
    
    if( !display_cal )
      display_cal = old_disp_cal;
    
    string missing_dets, extra_dets;
    for( const auto &det : tochange.detectors )
    {
      const auto pos = det_to_cal.find(det);
      
      if( (pos == end(det_to_cal)) || !pos->second || !pos->second->valid() )
        missing_dets += (missing_dets.empty() ? "" : ", ") + det;
    }
    
    for( const auto &det_cal : det_to_cal )
    {
      if( !tochange.detectors.count(det_cal.first) )
        extra_dets += (extra_dets.empty() ? "" : ", ") + det_cal.first;
    }
    
    // If we have extra detectors, no problem - if we dont have calibration for some detectors,
    //  its a bigger problem.
    if( missing_dets.size() )
    {
      WString msg = WString::tr("ect-CALp-missing-det").arg( missing_dets );
      if( extra_dets.size() )
        msg.arg( WString::tr("ect-CALp-extra-det").arg( extra_dets ) );
      else
        msg.arg( "" );
      
      throw runtime_error( msg.toUTF8() );
    }//if( missing_dets.size() )
    
    const set<string> detectors_to_apply = tochange.detectors;
    for( const string &det : detectors_to_apply )
    {
      const auto pos = det_to_cal.find(det);
      
      // TODO: we can probably do a better job of _trying_ something...
      if( pos == end(det_to_cal) )
      {
        if( det.empty() )
          throw runtime_error( WString::tr("ect-CALp-named-vs-unamed").toUTF8() );
        throw runtime_error( WString::tr("ect-CALp-no-cal-for-det").arg(det).toUTF8() );
      }//if( pos == end(det_to_cal) )
        
      const shared_ptr<const SpecUtils::EnergyCalibration> new_cal = pos->second;
      assert( new_cal && new_cal->valid() );
      if( !new_cal || !new_cal->valid() )
        throw runtime_error( "unexpected logic error" );
      
      shared_ptr<const SpecUtils::EnergyCalibration> old_cal;
      for( const int sample : tochange.sample_numbers )
      {
        const auto m = measurment->measurement( sample, det );
        const auto cal = m ? m->energy_calibration() : nullptr;
        if( cal && cal->valid() )
        {
          old_cal = cal;
          break;
        }
      }
      
      if( !old_cal )
        throw runtime_error( WString::tr("ect-CALp-no-prev-for-det").arg(det).toUTF8() );
      
      if( old_cal->num_channels() != new_cal->num_channels() )
        throw runtime_error( WString::tr("ect-CALp-nchannel-mismatch-det")
                               .arg( static_cast<int>(new_cal->num_channels()) )
                               .arg( static_cast<int>(old_cal->num_channels()) )
                               .arg(det).toUTF8() );
      
      MeasToApplyCoefChangeTo det_specific_tochange = tochange;
      
      det_specific_tochange.detectors.clear();
      det_specific_tochange.detectors.insert( det );
      
      // We will adjust the peak energies, only if the calibrations being changed contain the
      //  display energy calibration.
      const bool adjust_peaks = !adjusted_peaks && (display_cal == old_cal);
      adjusted_peaks |= adjust_peaks;
      
      setEnergyCal( new_cal, det_specific_tochange, adjust_peaks );
    }//for( const string &det : detectors_to_apply )
  }//if( det_to_cal.size() == 1 ) / else
  
  
  m_interspec->refreshDisplayedCharts();
  doRefreshFromFiles();
}//void applyCALpEnergyCal(...)


SpecUtils::SpectrumType EnergyCalTool::typeOfCurrentlyShowingCoefficients() const
{
  assert( m_specTypeMenu );
  
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  const int selectedType = m_specTypeMenu->currentIndex();
  if( selectedType < 0 || selectedType > 2 )
  {
    assert( 0 );  //shouldnt ever happen, right?
    throw runtime_error( "EnergyCalTool::typeOfCurrentlyShowingCoefficients(): invalid spec type" );
  }
  
  return spectypes[selectedType];
}//SpecUtils::SpectrumType typeOfCurrentlyShowingCoefficients() const

/*
 //Commented out because it is unused/untested
std::string EnergyCalTool::detectorNameOfCurrentlyShowingCoefficients() const
{
  assert( m_specTypeMenu );
  
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  const int selectedType = m_specTypeMenu->currentIndex();
  if( selectedType < 0 || selectedType > 2 )
  {
    assert( 0 );  //shouldnt ever happen, right?
    return;
  }
  
  WMenu *detMenu = m_detectorMenu[selectedType];
  assert( detMenu );
  
  WMenuItem *detitem = detMenu->currentItem();
  if( !detitem || !detMenu->count() )
  {
    assert( 0 );
    return "";
  }
  
  return detitem->text().toUTF8();
}//std::string detectorNameOfCurrentlyShowingCoefficients() const
*/


EnergyCalImp::CALpDownloadResource *EnergyCalTool::calpResources()
{
  return m_calpResource;
}



void EnergyCalTool::handleRequestToUploadCALp()
{
  SimpleDialog *dialog = new SimpleDialog();
  WPushButton *closeButton = dialog->addButton( "Cancel" );
  WGridLayout *stretcher = new WGridLayout();
  stretcher->setContentsMargins( 0, 0, 0, 0 );
  dialog->contents()->setLayout( stretcher );
  dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowVisible,
                                  Wt::Horizontal | Wt::Vertical );
  WText *title = new WText( WString::tr("ect-import-CALp") );
  title->addStyleClass( "title" );
  stretcher->addWidget( title, 0, 0 );
  
  WText *t = new WText( WString::tr("ect-select-CALp") );
  stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  t->setTextAlignment( Wt::AlignCenter );
  
  
  WFileUpload *upload = new WFileUpload();
  upload->fileTooLarge().connect( std::bind( [=](){
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    WPushButton *closeButton = dialog->addButton( WString::tr("Close") );
    WGridLayout *stretcher = new WGridLayout();
    stretcher->setContentsMargins( 0, 0, 0, 0 );
    dialog->contents()->setLayout( stretcher );
    WText *title = new WText( WString::tr("ect-upload-CALp-to-large") );
    title->addStyleClass( "title" );
    stretcher->addWidget( title, 0, 0 );
  }) );
  
  upload->changed().connect( upload, &WFileUpload::upload );
  upload->uploaded().connect( std::bind( [dialog,upload](){
    InterSpec *interspec = InterSpec::instance();
    SpecMeasManager *measmn = interspec ? interspec->fileManager() : nullptr;
    
    assert( measmn );
    if( !measmn )
      return;
    
    ifstream input( upload->spoolFileName().c_str(), ios::in | ios::binary );
    
    if( !measmn->handleCALpFile( input, dialog, true ) )
    {
      dialog->contents()->clear();
      dialog->footer()->clear();
      
      WPushButton *closeButton = dialog->addButton( WString::tr("Close") );
      WGridLayout *stretcher = new WGridLayout();
      stretcher->setContentsMargins( 0, 0, 0, 0 );
      dialog->contents()->setLayout( stretcher );
      WText *title = new WText( WString::tr("ect-invalid-CALp") );
      title->addStyleClass( "title" );
      stretcher->addWidget( title, 0, 0 );
      
      return;
    }//if( was not a valid CALp file )
    
    //wApp->doJavaScript( "$('.Wt-dialogcover').hide();" ); // JIC
    //dialog->done( Wt::WDialog::DialogCode::Accepted );
  } ) );
  
  stretcher->addWidget( upload, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  
  InterSpec *interspec = InterSpec::instance();
  if( interspec && !interspec->isPhone() )
  {
    t = new WText( WString::tr("ect-CALp-drag-n-drop-note") );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
  }
  
  /*
   //In case we want to use AuxWindow instead of SimpleDialog
   AuxWindow *window = new AuxWindow( "Import CALp file",
   (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
   | AuxWindowProperties::PhoneNotFullScreen
   | AuxWindowProperties::DisableCollapse
   | AuxWindowProperties::SetCloseable) );
   
   //...
   
   window->rejectWhenEscapePressed();
   window->show();
   window->resizeToFitOnScreen();
   window->centerWindow();
   
   WPushButton *close = window->addCloseButtonToFooter( "Cancel" );
   close->clicked().connect( boost::bind( &AuxWindow::hide, window ) );
   
   window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
   
   // TODO: add link to relevant section of documentation
   //AuxWindow::addHelpInFooter( window->footer(), "energy-cal-CALp" );
   */
}//void handleRequestToUploadCALp();



void EnergyCalTool::applyCalChange( std::shared_ptr<const SpecUtils::EnergyCalibration> disp_prev_cal,
                                    std::shared_ptr<const SpecUtils::EnergyCalibration> new_disp_cal,
                                    const vector<MeasToApplyCoefChangeTo> &changemeas,
                                    const bool isOffsetOnly )
{
  // We get here when:
  //  - we manually change a coefficiecnt on the GUI
  //  - we fit calibration coefficients
  //  - sometimes when we load in a CALp file
  
  using namespace SpecUtils;
  
  assert( disp_prev_cal && disp_prev_cal->valid() );
  assert( new_disp_cal && new_disp_cal->valid() );
  
  EnergyCalUndoRedoSentry undo_sentry;
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  const set<int> &foresamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  
  
  // Create a cache of modified calibration both to save time/memory, but also keep it so previous
  //  samples that share a energy calibration will continue to do so (if possible based on what user
  //  wanted calibration applied to).  Also, we wont set any new calibrations until we know all
  //  updated calibrations and peaks are valid
  //Note: we could take this oppritunity to share calibration across SpecFile objects by not just
  //      comparing pointers, but also the actual EnergyCalibration object.  But for now we'll
  //      skip this to avoid trouble, and it isnt clear that it would actually be overall beneficial
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  
  // We will store updated peaks and not set any of them until we know all the energy calibrations
  //  and peak shifts were successfully done.
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;
  
  // const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  //We will loop over the changes to apply twice.  Once to calculate new calibrations, and make sure
  //  they are valid, then a second time to actually set them.  If a new calibration is invalid,
  //  an exception will be thrown so we will catch that.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    //string dbgmsg = "For '" + change.meas->filename() + "' will apply changes to Detectors: {";
    //for( auto iter = begin(change.detectors); iter != end(change.detectors); ++iter )
    //  dbgmsg += (iter==begin(change.detectors) ? "" : ",") + (*iter);
    //dbgmsg += "} and Samples: {";
    //for( auto iter = begin(change.sample_numbers); iter != end(change.sample_numbers); ++iter )
    //  dbgmsg += (iter==begin(change.sample_numbers) ? "" : ",") + std::to_string(*iter);
    //dbgmsg += "}";
    //cout << dbgmsg << endl;
    //wApp->log("app:debug") << dbgmsg;
    
    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          auto m = change.meas->measurement( sample, detname );
          if( !m || m->num_gamma_channels() <= 4 )
            continue;
          
          const auto meas_old_cal = m->energy_calibration();
          assert( meas_old_cal );
          
          if( !meas_old_cal || !meas_old_cal->valid() )
            continue;
          
          //If we have already computed the new calibration for a EnergyCalibration object, lets not
          //  re-due it.
          if( old_to_new_cals.count(meas_old_cal) )
            continue;
          
          shared_ptr<const SpecUtils::EnergyCalibration> new_meas_cal;
          if( meas_old_cal == disp_prev_cal )
          {
            new_meas_cal = new_disp_cal;
          }else if( isOffsetOnly )
          {
            const vector<float> &new_disp_coefs = new_disp_cal->coefficients();
            const vector<float> &prev_disp_coefs = disp_prev_cal->coefficients();
            const vector<pair<float,float>> &dev_pairs = meas_old_cal->deviation_pairs();
            
            vector<float> new_coefs = meas_old_cal->coefficients();
            new_coefs[0] += (new_disp_coefs[0] - prev_disp_coefs[0]);
            
            auto cal = make_shared<EnergyCalibration>();
            switch( meas_old_cal->type() )
            {
              case SpecUtils::EnergyCalType::Polynomial:
              case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                cal->set_polynomial( meas_old_cal->num_channels(), new_coefs, dev_pairs );
                break;
                
              case SpecUtils::EnergyCalType::FullRangeFraction:
                cal->set_full_range_fraction( meas_old_cal->num_channels(), new_coefs, dev_pairs );
                break;
                
              case SpecUtils::EnergyCalType::LowerChannelEdge:
                cal->set_lower_channel_energy( meas_old_cal->num_channels(), new_coefs ); //eh, whatever
                break;
                
              case SpecUtils::EnergyCalType::InvalidEquationType:
                assert( 0 );
                break;
            }//switch( meas_old_cal->type() )
            
            new_meas_cal = cal;
          }else
          {
            new_meas_cal = EnergyCal::propogate_energy_cal_change( disp_prev_cal, new_disp_cal, meas_old_cal );
          }
          assert( new_meas_cal && new_meas_cal->valid() );
          old_to_new_cals[meas_old_cal] = new_meas_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      WString msg = WString::tr("ect-change-made-invalid");
      if( (backgrnd && (backgrnd != change.meas)) || (secgrnd && (secgrnd != change.meas)) )
      {
        if( change.meas == forgrnd )
          msg.arg( WString::tr("ect-for-the-fore") );
        else if( change.meas == backgrnd )
          msg.arg( WString::tr("ect-for-the-back") );
        else if( change.meas == secgrnd )
          msg.arg( WString::tr("ect-for-the-sec") );
        else
          msg.arg( "" );
      }else
      {
        msg.arg( "" );
      }//if( it is necessary to say which spectrum had the error )
      
      msg.arg( e.what() );
      
      throw runtime_error( msg.toUTF8() );
    }//try catch
    
    
    // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
    //  until we know we can update all the peaks
    
    const set<set<int>> samples_with_peaks = change.meas->sampleNumsWithPeaks();
    const set<set<int>> samples_with_hint_peaks = change.meas->sampleNumsWithAutomatedSearchPeaks();
    
    // We may not be updating all samples, so we will only update peaks who are owned
    //  by sample numbers that are all in the samples being updated.
    set<set<int>> peaksamples, hintPeakSamples;
    
    for( const set<int> &samples : samples_with_peaks )
    {
      bool all_samples = true;
      
      // Check if the peaks sample numbers are all getting re-calibrated
      for( auto sample_num_iter = begin(samples);
          all_samples && (sample_num_iter != end(samples));
          ++sample_num_iter )
      {
        all_samples = (change.sample_numbers.count(*sample_num_iter) != 0u);
      }
      
      if( all_samples )
        peaksamples.insert( samples );
    }//for( const set<int> &samples : samples_with_peaks )
    
    for( const set<int> &samples : samples_with_hint_peaks )
    {
      bool all_samples = true;

      // Check if the peaks sample numbers are all getting re-calibrated
      for( auto sample_num_iter = begin(samples);
          all_samples && (sample_num_iter != end(samples));
          ++sample_num_iter )
      {
        all_samples = (change.sample_numbers.count(*sample_num_iter) != 0u);
      }
      
      if( all_samples )
        hintPeakSamples.insert( samples );
    }//for( const set<int> &samples : samples_with_peaks )
    
    
    // The peaks position (i.e., mean channel number) is determined by
    //  #SpecFile::suggested_sum_energy_calibration, however we may be applying an energy
    //  calibration change to one or two detectors, that arent that detector, in which case we
    //  dont want to move the peaks.
    //  This is why we dont just use change.detectors.
    vector<string> display_detectors;
    if( change.meas == forgrnd )
    {
      display_detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
    }else if( change.meas == backgrnd )
    {
      display_detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Background);
    }else if( change.meas == secgrnd )
    {
      display_detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::SecondForeground);
    }
    
    if( display_detectors.empty() )
    {
      // We probably shouldnt ever get here, but if we do, it will probably be fine to just use
      //  change.detectors (
      cerr << __func__ <<  ": Apparently no changing For/Back/Sec spectrum - doesnt seem right!" << endl;
      assert(0);
      display_detectors.insert( end(display_detectors), begin(change.detectors), end(change.detectors) );
    }//if( foreground being change ) / else background / else
    
    
    for( const set<int> &samples : peaksamples )
    {
      //If there is any overlap between 'samples' and 'change.sample_numbers', then apply the change
      //  Note: this isnt correct, but I cant think of a better solution at the moment.
      auto oldpeaks = change.meas->peaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, display_detectors );
      
      if( !oldpeaks || oldpeaks->empty() || !oldcal || !oldcal->valid() )
      {
        if( !oldpeaks || !oldcal || !oldcal->valid() )
        {
          //development sanity check, shouldnt normally get here I dont think
          cerr << __func__ << "Failed to get peaks or oldcal!" << endl; //just for sanity check
          assert( 0 );
        }
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        //development sanity check, shouldnt normally happen I dont think
        cerr << __func__ << "Failed to get newcal for peaks shift!" << endl;
        assert( 0 );
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      
      if( oldcal == newcal )
      {
        cerr << __func__ <<  ": oldcal == newcal - skipping shifting peak" << endl;
        continue;
      }
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
        updated_peaks[oldpeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
        throw runtime_error( msg );
      }//try / catch
    }//for( const set<int> &samples : peaksampels )
    
    
    // Do similar loop to the above, but for hint peaks
    for( const set<int> &samples : hintPeakSamples )
    {
      auto oldHintPeaks = change.meas->automatedSearchPeaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, display_detectors );
      
      if( !oldHintPeaks || oldHintPeaks->empty() || !oldcal || !oldcal->valid() )
      {
        if( !oldHintPeaks || !oldcal || !oldcal->valid() )
        {
          assert( 0 );  //shouldnt get here
        }
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        assert( 0 ); //shouldnt get here
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      if( oldcal == newcal )
        continue;
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks, oldcal, newcal );
        updated_hint_peaks[oldHintPeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating hint peaks for this energy change,"
        " still applying, Error: " + string(e.what());
  #if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
  #endif
      }//try / catch
    }//for( const set<int> &samples : hintPeakSamples )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  if( old_to_new_cals.find(disp_prev_cal) == end(old_to_new_cals) )
  {
    //Shouldnt ever happen; check is for development
    string msg = "There was an internal error updating energy calibration - energy cal"
    " associated with GUI wasn't updated - energy calibration state is suspect";
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.c_str() );
#endif
    
    m_interspec->logMessage( msg, 3 );
    
#if( PERFORM_DEVELOPER_CHECKS )
    assert( 0 );
#endif
  }//if( old_to_new_cals.find(disp_prev_cal) == end(old_to_new_cals) )
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info( change.meas );
    meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info( change.meas );
    meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info( change.meas );
    
    
    for( const int sample : change.sample_numbers )
    {
      for( const string &detname : change.detectors )
      {
        auto m = change.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
        
        const auto measoldcal = m->energy_calibration();
        assert( measoldcal );
        
        auto iter = old_to_new_cals.find( measoldcal );
        if( iter == end(old_to_new_cals) )
        {
          //Shouldnt ever happen
          string msg = "There was an internal error updating energy calibration - precomputed"
          " calibration couldnt be found - energy calibration will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
          
          m_interspec->logMessage( msg, 3 );
          assert( 0 );
          continue;
        }//if( we havent already computed a new energy cal )
        
        assert( iter->second );
        assert( iter->second->num_channels() == m->num_gamma_channels() );
        
        meas_old_new_cal.emplace_back( m, measoldcal, iter->second );
        
        change.meas->set_energy_calibration( iter->second, m );
      }//for( loop over detector names )
    }//for( loop over sample numbers )
    
    
    //Now actually set the updated peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    
    for( const set<int> &samples : peaksamples )
    {
      auto oldpeaks = change.meas->peaks(samples);
      if( oldpeaks )
      {
        const auto pos = updated_peaks.find(oldpeaks);
        if( pos == end(updated_peaks) )
        {
          if( oldpeaks && !oldpeaks->empty() )
            cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        }else
        {
          meas_old_new_peaks.emplace_back( samples, *oldpeaks, pos->second );
          
          change.meas->setPeaks( pos->second, samples );
          if( m_peakModel && (change.meas == forgrnd) && (samples == foresamples) )
            m_peakModel->setPeakFromSpecMeas(forgrnd, foresamples);
        }//if( pos == end(updated_peaks) ) / else
      }//if( oldpeaks )
    }//for( const set<int> &samples : peaksampels )
    
    
    // Do similar thing for the hint peaks
    const set<set<int>> hintPeakSamples = change.meas->sampleNumsWithAutomatedSearchPeaks();
    for( const set<int> &samples : hintPeakSamples )
    {
      shared_ptr<const SpecMeas::PeakDeque> hintPeaks = change.meas->automatedSearchPeaks(samples);
      assert( hintPeaks && !hintPeaks->empty() );
      
      if( hintPeaks && !hintPeaks->empty() )
      {
        const auto pos = updated_hint_peaks.find( hintPeaks );
        if( pos == end(updated_hint_peaks) )
        {
          if( hintPeaks && !hintPeaks->empty() )
            cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        }else if( !pos->second.empty() )
        {
          meas_old_new_hint_peaks.emplace_back( samples, *hintPeaks, pos->second );
          
          auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( pos->second );
          change.meas->setAutomatedSearchPeaks( samples, peaks );
        }//if( pos == end(updated_hint_peaks) ) / else
      }//if( hintPeaks )
    }//for( const set<int> &samples : hintPeakSamples )
  }//for( loop over SpecFiles for change )
  
  m_interspec->refreshDisplayedCharts();
  doRefreshFromFiles();
}//applyCalChange(...)


void EnergyCalTool::setEnergyCal( shared_ptr<const SpecUtils::EnergyCalibration> new_cal,
                   const MeasToApplyCoefChangeTo &changemeas, const bool adjust_peaks )
{
  // We get here from applying CALp to files, sometimes.
  using namespace SpecUtils;
  
  assert( new_cal && new_cal->valid() );
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  // Currently its pointless to map old energy calibrations to the new_cal being set, but in the
  //  future we might be a little more lenient about matching number of channels or whatever, so
  //  we'll just toss this mechanism in for now.
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  
  // We will store updated peaks and not set any of them until we know all the energy calibrations
  //  and peak shifts were successfully done.
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;
  
  
  // const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  //First calculate new calibrations, and make sure they are valid, then actually set them.
  assert( changemeas.meas );
  if( !changemeas.meas )
    throw logic_error( "EnergyCalTool::setEnergyCal: nullptr measurement to change passed in." );
  
  SpecUtils::SpectrumType spectype;
  if( changemeas.meas == forgrnd )
    spectype = SpecUtils::SpectrumType::Foreground;
  else if( changemeas.meas == backgrnd )
    spectype = SpecUtils::SpectrumType::Background;
  else if( changemeas.meas == secgrnd )
    spectype = SpecUtils::SpectrumType::SecondForeground;
  else
  {
    assert( 0 );
    throw runtime_error( "EnergyCalTool::setEnergyCal: measurement to change must be foreground,"
                         " background, or secondary spec displayed" );
  }
  
  
  const set<int> &dispsamples = m_interspec->displayedSamples( spectype );
  const vector<string> dispdets = m_interspec->detectorsToDisplay( spectype );
  const shared_ptr<const SpecUtils::EnergyCalibration> display_cal
                       = changemeas.meas->suggested_sum_energy_calibration( dispsamples, dispdets );
  
  try
  {
    for( const int sample : changemeas.sample_numbers )
    {
      for( const string &detname : changemeas.detectors )
      {
        auto m = changemeas.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
          
        const auto meas_old_cal = m->energy_calibration();
        assert( meas_old_cal );
          
        if( !meas_old_cal || !meas_old_cal->valid() )
          continue;
          
        // TODO: maybe be a little more lenient on matching number of channels, by possible making
        //       new energy calibrations when possible.
        
        if( meas_old_cal->num_channels() != new_cal->num_channels() )
        {
          throw runtime_error( WString::tr("ect-set-cal-invalid-nchan")
                              .arg(sample)
                              .arg(detname)
                              .arg( static_cast<int>(meas_old_cal->num_channels()) )
                              .arg( static_cast<int>(new_cal->num_channels()) )
                              .toUTF8() );
        }
        //if( display_cal == meas_old_cal )
        //  adjust_peaks = true;
        
        old_to_new_cals[meas_old_cal] = new_cal;
      }//for( const string &detname : change.detectors )
    }//for( loop over sample numbers )
  }catch( std::exception &e )
  {
    throw runtime_error( WString::tr("ect-error-setting-cal").arg(e.what()).toUTF8() );
  }//try catch
  
  // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
  //  until we know we can update all the peaks
  const set<set<int>> peaksamples = adjust_peaks ? changemeas.meas->sampleNumsWithPeaks()
                                                 : set<set<int>>();
    
  // The peaks position (i.e., mean channel number) is determined by
  //  #SpecFile::suggested_sum_energy_calibration, however we may be applying an energy
  //  calibration change to one or two detectors, that arent that detector, in which case we
  //  dont want to move the peaks.
  for( const set<int> &samples : peaksamples )
  {
    //If there is any overlap between 'samples' and 'change.sample_numbers', then apply the change
    //  Note: this isnt correct, but I cant think of a better solution at the moment.
    auto oldpeaks = changemeas.meas->peaks(samples);
    auto oldcal = changemeas.meas->suggested_sum_energy_calibration( samples, dispdets );
      
    if( !oldpeaks || oldpeaks->empty() || !oldcal || !oldcal->valid() )
    {
      if( !oldpeaks || !oldcal || !oldcal->valid() )
      {
        // we can get here if we dont have any display detectors
        cerr << __func__ << "Failed to get peaks or oldcal!" << endl; //just for sanity check
      }
      continue;
    }
    
    const auto newcal_pos = old_to_new_cals.find(oldcal);
    if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
    {
      cout << __func__ << ": not applying energy cal change to peaks." << endl;
      continue;
    }
      
    const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
    assert( newcal && newcal->valid() );
      
    if( oldcal == newcal )
    {
      assert( 0 );
      continue;
    }
      
    try
    {
      auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
      updated_peaks[oldpeaks] = newpeaks;
    }catch( std::exception &e )
    {
      string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
      throw runtime_error( msg );
    }//try / catch
  }//for( const set<int> &samples : peaksampels )
  
  
  // Now do same thing as for peaks, but for the hint peaks
  const set<set<int>> hintPeakSamples = adjust_peaks
                                            ? changemeas.meas->sampleNumsWithAutomatedSearchPeaks()
                                            : set<set<int>>();
  for( const set<int> &samples : hintPeakSamples )
  {
    auto oldHintPeaks = changemeas.meas->automatedSearchPeaks( samples );
    auto oldcal = changemeas.meas->suggested_sum_energy_calibration( samples, dispdets );
      
    if( !oldHintPeaks || oldHintPeaks->empty() || !oldcal || !oldcal->valid() )
      continue;
    
    const auto newcal_pos = old_to_new_cals.find(oldcal);
    if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      continue;
      
    const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
    assert( newcal && newcal->valid() );
      
    if( !newcal || !newcal->valid() || (oldcal == newcal) )
    {
      assert( 0 );
      continue;
    }
    
    try
    {
      auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks, oldcal, newcal );
      updated_hint_peaks[oldHintPeaks] = newpeaks;
    }catch( std::exception &e )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      string msg = "There was an issue translating hint peaks for this energy change.  Error: "
                   + string(e.what());
      log_developer_error( __func__, msg.c_str() );
#endif
    }//try / catch
  }//for( const set<int> &samples : hintPeakSamples )
  
    
  // For undo/redo, store changes to energy cal, and peaks.
  EnergyCalUndoRedoSentry undo_sentry;
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks
  for( const int sample : changemeas.sample_numbers )
  {
    meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info(changemeas.meas);
    
    for( const string &detname : changemeas.detectors )
    {
      auto m = changemeas.meas->measurement( sample, detname );
      if( !m || m->num_gamma_channels() <= 4 )
        continue;
        
      const auto measoldcal = m->energy_calibration();
      assert( measoldcal );
        
      auto iter = old_to_new_cals.find( measoldcal );
      if( iter == end(old_to_new_cals) )
      {
        //Shouldnt ever happen
        string msg = "There was an internal error updating energy calibration - precomputed"
        " calibration couldnt be found - energy calibration will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
          
        m_interspec->logMessage( msg, 3 );
        assert( 0 );
        continue;
      }//if( we havent already computed a new energy cal )
        
      assert( iter->second );
      assert( iter->second->num_channels() == m->num_gamma_channels() );
        
      meas_old_new_cal.emplace_back( m, measoldcal, iter->second );
      
      changemeas.meas->set_energy_calibration( iter->second, m );
    }//for( loop over detector names )
  }//for( loop over sample numbers )
    
    
  //Now actually set the updated peaks
  const set<int> &foresamples = m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );
  for( const set<int> &samples : peaksamples )
  {
    auto oldpeaks = changemeas.meas->peaks(samples);
    if( oldpeaks )
    {
      const auto pos = updated_peaks.find(oldpeaks);
      if( pos == end(updated_peaks) )
      {
        if( oldpeaks && !oldpeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_peaks" << endl;
      }else
      {
        meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info( changemeas.meas );
        meas_old_new_peaks.emplace_back( samples, *oldpeaks, pos->second );
        
        changemeas.meas->setPeaks( pos->second, samples );
        if( m_peakModel && (changemeas.meas == forgrnd) && (samples == foresamples) )
          m_peakModel->setPeakFromSpecMeas(forgrnd, foresamples);
      }//if( pos == end(updated_peaks) ) / else
    }//if( oldpeaks )
  }//for( const set<int> &samples : peaksampels )
  
  // And set the updated hint peaks
  for( const set<int> &samples : hintPeakSamples )
  {
    auto hintPeaks = changemeas.meas->automatedSearchPeaks(samples);
    if( hintPeaks && !hintPeaks->empty() )
    {
      const auto pos = updated_hint_peaks.find(hintPeaks);
      if( pos == end(updated_hint_peaks) )
      {
        if( hintPeaks && !hintPeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_hint_peaks" << endl;
      }else
      {
        meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info( changemeas.meas );
        meas_old_new_hint_peaks.emplace_back( samples, *hintPeaks, pos->second );
        auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( pos->second );
        changemeas.meas->setAutomatedSearchPeaks( samples, peaks );
      }//if( pos == end(updated_hint_peaks) ) / else
    }//if( hintPeaks && !hintPeaks->empty() )
  }//for( const set<int> &samples : peaksampels )
}//void setEnergyCal( new_cal, changemeas )


void EnergyCalTool::addDeviationPair( const std::pair<float,float> &new_pair )
{
  // We get here when we graphically (i.e., cntrl+option+drag) add in a deviation pair
  
  using namespace SpecUtils;
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  const set<int> &foresamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  
  // We will calculate all new energy calibrations, to make sure we actually can, befor setting them
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  
  // We will also pre-calculate updated peaks
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;
  
  const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  // For undo/redo, store changes to energy cal, and peaks.
  EnergyCalUndoRedoSentry undo_sentry;
  
  // Do first loop to calculate new calibrations
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info(change.meas);
    meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info(change.meas);
    
    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          auto m = change.meas->measurement( sample, detname );
          if( !m || m->num_gamma_channels() <= 4 )
            continue;
          
          const auto old_cal = m->energy_calibration();
          assert( old_cal );
          
          if( !m || (m->num_gamma_channels() <= 4) || !old_cal || !old_cal->valid()
             || (old_cal->type() == EnergyCalType::LowerChannelEdge) )
            continue;
          
          //If we have already computed the new calibration for a EnergyCalibration object, lets not
          //  re-due it.
          if( old_to_new_cals.count(old_cal) )
            continue;
          
          const size_t nchannel = old_cal->num_channels();
          const vector<float> &coefs = old_cal->coefficients();
          vector<pair<float,float>> dev_pairs = old_cal->deviation_pairs();
          if( dev_pairs.empty() )
            dev_pairs.push_back( {0.0f, 0.0f} );
          dev_pairs.push_back( new_pair );
          
          std::sort( begin(dev_pairs), end(dev_pairs) );
          
          shared_ptr<EnergyCalibration> new_cal = make_shared<EnergyCalibration>();
          
          switch( old_cal->type() )
          {
            case SpecUtils::EnergyCalType::Polynomial:
            case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
              new_cal->set_polynomial( nchannel, coefs, dev_pairs );
              break;
              
            case SpecUtils::EnergyCalType::FullRangeFraction:
              new_cal->set_full_range_fraction( nchannel, coefs, dev_pairs );
              break;
              
            case SpecUtils::EnergyCalType::LowerChannelEdge:
            case SpecUtils::EnergyCalType::InvalidEquationType:
              assert( 0 );
              break;
          }//switch( meas_old_cal->type() )
          
          assert( new_cal && new_cal->valid() );
          old_to_new_cals[old_cal] = new_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      WString msg = WString::tr("ect-add-dev-pair-made-invalid");
      if( (backgrnd && (backgrnd != change.meas)) || (secgrnd && (secgrnd != change.meas)) )
      {
        if( change.meas == forgrnd )
          msg.arg( WString::tr("ect-for-the-fore") );
        else if( change.meas == backgrnd )
          msg.arg( WString::tr("ect-for-the-back") );
        else if( change.meas == secgrnd )
          msg.arg( WString::tr("ect-for-the-sec") );
        else
          msg.arg( "" );
      }else
      {
        msg.arg( "" );
      }//if( it is necessary to say which spectrum had the error ) / else
      
      throw runtime_error( msg.arg(e.what()).toUTF8() );
    }//try catch
    
    
    // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
    //  until we know we can update all the peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    const vector<string> detnamesv( begin(change.detectors), end(change.detectors) );
    
    for( const set<int> &samples : peaksamples )
    {
      //If there is any overlap between 'samples' and 'change.sample_numbers', then apply the change
      //  Note: this isnt correct, but I cant think of a better solution at the moment.
      auto oldpeaks = change.meas->peaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );
      
      if( !oldpeaks || oldpeaks->empty()
         || !oldcal || (oldcal->type() == EnergyCalType::InvalidEquationType) )
      {
        if( !oldpeaks || !oldcal || !oldcal->valid() )
          cerr << "Failed to get peaks or oldcal!" << endl; //just for development
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        cerr << "Failed to get newcal for peaks shift!" << endl; //just for development, shouldnt happen I dont think
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      
      if( oldcal == newcal )
      {
        cerr << __func__ <<  ": oldcal == newcal - skipping shifting peak" << endl;
        continue;
      }
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
        updated_peaks[oldpeaks] = newpeaks;
        
        meas_old_new_peaks.push_back( {samples, *oldpeaks, newpeaks} );
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
        throw runtime_error( msg );
      }//try / catch
    }//for( const set<int> &samples : peaksampels )
    
    
    //Now get new hint peaks
    const set<set<int>> hintPeakSamples = change.meas->sampleNumsWithAutomatedSearchPeaks();
    
    for( const set<int> &samples : hintPeakSamples )
    {
      auto oldHintPeaks = change.meas->automatedSearchPeaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );
      
      if( !oldHintPeaks || oldHintPeaks->empty()
         || !oldcal || (oldcal->type() == EnergyCalType::InvalidEquationType) )
      {
        if( !oldHintPeaks || !oldcal || !oldcal->valid() )
          cerr << "Failed to get peaks or oldcal!" << endl; //just for development
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        cerr << "Failed to get newcal for peaks shift!" << endl; //just for development, shouldnt happen I dont think
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      if( !newcal || !newcal->valid() || (oldcal == newcal) )
        continue;
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks, oldcal, newcal );
        updated_hint_peaks[oldHintPeaks] = newpeaks;
        
        meas_old_new_hint_peaks.push_back( {samples, *oldHintPeaks, newpeaks} );
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating hint peaks for this energy change/"
        "  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
      }//try / catch
    }//for( const set<int> &samples : peaksampels )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info( change.meas );
    
    for( const int sample : change.sample_numbers )
    {
      for( const string &detname : change.detectors )
      {
        auto m = change.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
        
        const auto measoldcal = m->energy_calibration();
        assert( measoldcal );
        
        auto iter = old_to_new_cals.find( measoldcal );
        if( iter == end(old_to_new_cals) )
        {
          //Shouldnt ever happen
          string msg = "There was an internal error updating energy calibration - precomputed"
          " calibration couldnt be found - energy calibation will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
          
          m_interspec->logMessage( msg, 3 );
          assert( 0 );
          continue;
        }//if( we havent already computed a new energy cal )
        
        assert( iter->second );
        assert( iter->second->num_channels() == m->num_gamma_channels() );
        
        change.meas->set_energy_calibration( iter->second, m );
        
        meas_old_new_cal.emplace_back( m, measoldcal, iter->second );
      }//for( loop over detector names )
    }//for( loop over sample numbers )
    
    
    //Now actually set the updated peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    
    for( const set<int> &samples : peaksamples )
    {
      auto oldpeaks = change.meas->peaks(samples);
      if( !oldpeaks || oldpeaks->empty() )
        continue;
      
      const auto pos = updated_peaks.find(oldpeaks);
      if( pos == end(updated_peaks) )
      {
        if( oldpeaks && !oldpeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        continue;
      }
      
      change.meas->setPeaks( pos->second, samples );
      if( m_peakModel && (change.meas == forgrnd) && (samples == foresamples) )
        m_peakModel->setPeakFromSpecMeas(forgrnd, foresamples);
    }//for( const set<int> &samples : peaksampels )
    
    // Also grab the updated hint peaks
    const set<set<int>> hintPeakSamples = change.meas->sampleNumsWithAutomatedSearchPeaks();
    for( const set<int> &samples : hintPeakSamples )
    {
      auto oldHintPeaks = change.meas->automatedSearchPeaks( samples );
      if( !oldHintPeaks || oldHintPeaks->empty() )
        continue;
      
      const auto pos = updated_hint_peaks.find(oldHintPeaks);
      if( pos == end(updated_hint_peaks) )
      {
        if( oldHintPeaks && !oldHintPeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_hint_peaks" << endl;
      }else
      {
        auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( pos->second );
        change.meas->setAutomatedSearchPeaks( samples, peaks );
      }
    }//for( const set<int> &samples : peaksampels )
  }//for( loop over SpecFiles for change )
  

  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//void addDeviationPair( const std::pair<float,float> &new_pair );



void EnergyCalTool::userChangedCoefficient( const size_t coefnum, EnergyCalImp::CalDisplay *display )
{
  //cout << "EnergyCalTool::userChangedCoefficient" << endl;
  using namespace SpecUtils;
  assert( coefnum < 10 );  //If we ever allow lower channel energy adjustment this will need to be removed
  
  shared_ptr<const EnergyCalibration> disp_prev_cal = display->lastSetCalibration();
  if( !disp_prev_cal )
  {
    cerr << "unexpected error getting updated energy calibration coefficients" << endl;
    m_interspec->logMessage( WString::tr("ect-unexpected-error-prev-coefs"), 2 );
    doRefreshFromFiles();
    return;
  }//if( !disp_prev_cal )
  
  
  {// Begin check to make sure the changed energy cal is actually checked for it to be applied to...
    bool willBeAppliedToDisplay = false;
    const std::string &detname = display->detectorName();
    const SpecUtils::SpectrumType type = display->spectrumType();
    std::shared_ptr<SpecMeas> cal_disp_meas = m_interspec->measurment(type);
    assert( cal_disp_meas );
    
    const vector<MeasToApplyCoefChangeTo> applyTo = measurementsToApplyCoeffChangeTo();
    
    for( const MeasToApplyCoefChangeTo &delta : applyTo )
    {
      const shared_ptr<SpecMeas> &meas = delta.meas;
      if( meas != cal_disp_meas )
        continue;
      
      const set<string> &detectors = delta.detectors;
      if( !detectors.count(detname) )
        continue;
      
      // Actually, I think if we're here, we're probably good, but we'll check a little deeper to
      //   make sure check in EnergyCalTool::applyCalChange will be satisfied
      const set<int> &samples = delta.sample_numbers;
      for( const int sample : samples )
      {
        auto m = meas->measurement( sample, detname);
        willBeAppliedToDisplay = (m && (m->energy_calibration() == disp_prev_cal));
        if( willBeAppliedToDisplay )
          break;
      }//for( const int sample : samples )
    
      if( willBeAppliedToDisplay )
        break;
    }//for( loop over changes )
    
    if( !willBeAppliedToDisplay )
    {
      m_interspec->logMessage( WString::tr("ect-changed-cal-not-selected"), 2 );
      doRefreshFromFiles();
      return;
    }
  }// End check to make sure the changed energy cal is actually checked for it to be applied to...
  
  
  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  vector<float> dispcoefs = display->displayedCoefficents();
  if( dispcoefs.size() <= coefnum )
    dispcoefs.resize( coefnum+1, 0.0f );
  
  vector<float> prev_disp_coefs = disp_prev_cal->coefficients();
  if( prev_disp_coefs.size() <= coefnum )
    prev_disp_coefs.resize( coefnum+1, 0.0f );
  
  vector<float> new_disp_coefs = prev_disp_coefs;
  new_disp_coefs[coefnum] = dispcoefs[coefnum];
  
  const size_t dispnchannel = disp_prev_cal->num_channels();
  const auto &disp_dev_pairs = disp_prev_cal->deviation_pairs();
  
  shared_ptr<const EnergyCalibration> new_disp_cal;
  try
  {
    auto cal = make_shared<EnergyCalibration>();
    switch( disp_prev_cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        cal->set_polynomial( dispnchannel, new_disp_coefs, disp_dev_pairs );
        break;
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
        cal->set_full_range_fraction( dispnchannel, new_disp_coefs, disp_dev_pairs );
        break;
      
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "Invalid calibration type changed?  Something is way wack." );
        break;
    }//switch( disp_prev_cal->type() )
    
    new_disp_cal = cal;
  }catch( std::exception &e )
  {
    display->updateToGui( disp_prev_cal );
    m_interspec->logMessage( WString::tr("ect-change-made-invalid").arg("").arg(e.what()), 2 );
    
    return;
  }//try / catch to create new_disp_cal
  
  assert( new_disp_cal && new_disp_cal->valid() );
  
  try
  {
    const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
    applyCalChange( disp_prev_cal, new_disp_cal, changemeas, coefnum==0 );
  }catch( std::exception &e )
  {
    display->updateToGui( disp_prev_cal );
    m_interspec->logMessage( WString::tr("ect-change-made-invalid").arg("").arg(e.what()), 2 );
  }//try / catch
}//userChangedCoefficient(...)


void EnergyCalTool::userChangedDeviationPair( EnergyCalImp::CalDisplay *display, const int fieldTypeChanged )
{
  using namespace SpecUtils;
  
  const auto field = EnergyCalImp::DeviationPairDisplay::UserFieldChanged( fieldTypeChanged );
  
  switch( field )
  {
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::AddedDeviationPair:
      return;  //hasnt been filled out yet; no need to do anything
      
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::RemovedDeviationPair:
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::EnergyChanged:
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::OffsetChanged:
      break;
  };//enum UserFieldChanged

  
  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const set<int> &foresamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  
  assert( forgrnd );
  if( !forgrnd )
    return;
  
  const shared_ptr<const EnergyCalibration> old_cal = display->lastSetCalibration();
  const vector<pair<float,float>> old_dev_pairs = old_cal
                                        ? old_cal->deviation_pairs() : vector<pair<float,float>>{};
  const vector<pair<float,float>> new_dev_pairs = display->displayedDeviationPairs();
  
  // After the user clicks to add a new deviation pair, and then fills in one of the fields, if the
  //  other field is not yet filled in, then the GUI wont insert that deviation pair; so for this
  //  case where the deviation pairs havent yet changed, lets skip doing anything yet
  if( new_dev_pairs.size() == old_dev_pairs.size() )
  {
    bool equal = true;
    for( size_t i = 0; equal && (i < new_dev_pairs.size()); ++i )
      equal = ( (fabs(new_dev_pairs[i].first - old_dev_pairs[i].first) < 1.0E-4)
                && (fabs(new_dev_pairs[i].second - old_dev_pairs[i].second) < 1.0E-4) );
    if( equal )
      return;
  }//if( new_dev_pairs.size() == old_dev_pairs.size() )
  
  
  const SpectrumType type = display->spectrumType();
  const std::string &detname = display->detectorName();
  
  auto specfile = m_interspec->measurment( type );
  if( !specfile )  //Shouldnt ever happen
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( "Internal error retrieving correct measurement", 2 );
    return;
  }

  // We will store updated calibrations and not set any of them until we know they can all be
  //  successfully altered
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cal;
  
  try
  {
    for( auto &m : specfile->measurements() )
    {
      if( (m->detector_name() != detname) || (m->num_gamma_channels() < 5) )
        continue;
      
      const auto cal = m->energy_calibration();
      if( !cal || !cal->valid() || (cal->type() == EnergyCalType::LowerChannelEdge) )
        continue;
      
      if( old_to_new_cal.find(cal) != end(old_to_new_cal) )
        continue;
      
      const size_t nchannel = cal->num_channels();
      const vector<float> &coefficients = cal->coefficients();
      
      auto new_cal = make_shared<EnergyCalibration>();
      switch( cal->type() )
      {
        case EnergyCalType::Polynomial:
        case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          new_cal->set_polynomial( nchannel, coefficients, new_dev_pairs );
          break;
          
        case EnergyCalType::FullRangeFraction:
          new_cal->set_full_range_fraction( nchannel, coefficients, new_dev_pairs );
          break;
          
        case EnergyCalType::LowerChannelEdge:
        case EnergyCalType::InvalidEquationType:
          assert( 0 );
          break;
      }//switch( old_cal->type() )
      
      assert( new_cal->valid() );
      old_to_new_cal[cal] = new_cal;
    }//for( auto &m : specfile->measurements() )
  }catch( std::exception &e )
  {
    display->updateToGui( old_cal );
    //display->setDeviationPairsInvalid();
    
    m_interspec->logMessage( WString::tr("ect-dev-pair-change-invalid").arg(e.what()), 2 );
    
    return;
  }//try / catch
  
  // We will store updated peaks and not set any of them until we know all the energy calibrations
  //  and peak shifts were sucessfully done.
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;
  
  //const vector<string> &detnames = specfile->gamma_detector_names();
  const vector<string> &detnames = m_interspec->detectorsToDisplay( type );
  const set<set<int>> samplesWithPeaks = specfile->sampleNumsWithPeaks();
  for( const set<int> &samples : samplesWithPeaks )
  {
    try
    {
      auto dispcal = specfile->suggested_sum_energy_calibration( samples, detnames );
      
      const auto dispcaliter = old_to_new_cal.find(dispcal);
      if( dispcaliter == end(old_to_new_cal) )
        continue;
      
      bool detInSample = false;
      for( auto iter = begin(samples); !detInSample && (iter != end(samples)); ++iter )
      {
        auto m = specfile->measurement( *iter, detname );
        detInSample = !!m;
      }
      
      if( !detInSample )
        continue;
      
      // I *think* that applying the update to deviation pairs shouldnt compound when you then apply
      //  them to another detector since the dispaly energy cal should update, but could there be
      //  any edge-cases?
      auto oldpeaks = specfile->peaks(samples);
      if( !oldpeaks || oldpeaks->empty() )
        continue;
      
    
      auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, dispcaliter->first,
                                                                     dispcaliter->second );
      updated_peaks[oldpeaks] = newpeaks;
    }catch( std::exception &e )
    {
      //display->updateToGui( old_cal );
      display->setDeviationPairsInvalid();
      
      string msg = "There was an issue translating peaks for this deviation pair change;"
                   " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, msg.c_str() );
#endif
      m_interspec->logMessage( msg, 2 );
      
      return;
    }//try / catch
  }//for( const set<int> &samples : samplesWithPeaks )
  
  
  const set<set<int>> samplesWithHintPeaks = specfile->sampleNumsWithAutomatedSearchPeaks();
  for( const set<int> &samples : samplesWithHintPeaks )
  {
    try
    {
      auto dispcal = specfile->suggested_sum_energy_calibration( samples, detnames );
      
      const auto dispcaliter = old_to_new_cal.find(dispcal);
      if( dispcaliter == end(old_to_new_cal) )
        continue;
      
      bool detInSample = false;
      for( auto iter = begin(samples); !detInSample && (iter != end(samples)); ++iter )
      {
        auto m = specfile->measurement( *iter, detname );
        detInSample = !!m;
      }
      
      if( !detInSample )
        continue;
      
      auto oldHintPeaks = specfile->automatedSearchPeaks(samples);
      if( !oldHintPeaks || oldHintPeaks->empty() )
        continue;
      
      auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks, dispcaliter->first,
                                                                    dispcaliter->second );
      updated_hint_peaks[oldHintPeaks] = newpeaks;
    }catch( std::exception &e )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      string msg = "There was an issue translating hint peaks for this deviation pair change."
      " Error: " + string(e.what());
      log_developer_error( __func__, msg.c_str() );
#endif
    }//try / catch
  }//for( const set<int> &samples : samplesWithHintPeaks )
  
  display->setDeviationPairsValid();
  
  
  size_t num_updated = 0;
  
  // Track some info for Undo/Redo
  EnergyCalUndoRedoSentry undo_sentry;
  meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info(forgrnd);
  meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info(forgrnd);
  meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info(forgrnd);
  
  for( auto &m : specfile->measurements() )
  {
    // I'm a little torn if we should update just the one energy calibration, or all occurances of
    //  the detectors energy calibration.
    //  Maybe we should update according to the checked GUI, but also restrict on name as well?
    //auto cal = m->energy_calibration();
    //if( cal == old_cal )
    if( (m->detector_name() != detname) || (m->num_gamma_channels() < 5) )
      continue;
    
    const auto cal = m->energy_calibration();
    if( !cal || !cal->valid() || (cal->type() == EnergyCalType::LowerChannelEdge) )
      continue;
    
    auto calpos = old_to_new_cal.find(cal);
    if( (calpos == end(old_to_new_cal)) || !calpos->second || !calpos->second->valid() )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Unexpectedly found invalid calibration in old_to_new_cal" );
#endif
      continue;
    }//if( sanity check that new calibration is valid - should always be )
    
    meas_old_new_cal.emplace_back( m, cal, calpos->second );
    
    specfile->set_energy_calibration( calpos->second, m );
    ++num_updated;
  }//for( loop over measurements )
  
  if( num_updated == 0 )
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( WString::tr("ect-set-dev-pair-err"), 2 );
    return;
  }
  
  // Now actually set the new peaks
  for( const set<int> &samples : samplesWithPeaks )
  {
    auto oldpeaks = specfile->peaks(samples);
    if( !oldpeaks || oldpeaks->empty() )
      continue;
    
    const auto peakpos = updated_peaks.find(oldpeaks);
    if( peakpos == end(updated_peaks) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Unexpectedly couldn't find peaks in updated_peaks" );
#endif
      continue;
    }//if( sanity check that shouldnt ever happen )
      
    meas_old_new_peaks.emplace_back( samples, *oldpeaks, peakpos->second );
    
    specfile->setPeaks( peakpos->second, samples );
    if( m_peakModel && (specfile == forgrnd) && (samples == foresamples) )
      m_peakModel->setPeakFromSpecMeas(forgrnd, foresamples);
  }//for( const set<int> &samples : samplesWithPeaks )
  
  // And set the automated hint peaks
  for( const set<int> &samples : samplesWithHintPeaks )
  {
    auto oldHintPeaks = specfile->automatedSearchPeaks(samples);
    if( !oldHintPeaks || oldHintPeaks->empty() )
      continue;
    
    const auto peakpos = updated_hint_peaks.find(oldHintPeaks);
    if( peakpos == end(updated_hint_peaks) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Unexpectedly coudlnt find peaks in updated_peaks" );
#endif
      continue;
    }//if( sanity check that shouldnt ever happen )
    
    meas_old_new_hint_peaks.emplace_back( samples, *oldHintPeaks, peakpos->second );
    
    auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( peakpos->second );
    specfile->setAutomatedSearchPeaks( samples, peaks );
  }//for( const set<int> &samples : samplesWithPeaks )
  
  const size_t ndets = specfile->gamma_detector_names().size();
  const size_t nsamples = specfile->sample_numbers().size();
  
  if( (ndets > 1) || (nsamples > 1) )
  {
    WString msg;
    if( (ndets > 1) && (nsamples > 1) )
      msg = WString::tr("ect-dev-applied-dets-samples").arg(detname);
    else if( nsamples > 1 )
      msg = WString::tr("ect-dev-applied-to-samples");
    
    int nfiles = 0;
    for( auto t : {0,1,2} )
      nfiles += (m_interspec->measurment(static_cast<SpectrumType>(t)) != specfile);
    
    if( nfiles && !msg.empty() )
    {
      switch( type )
      {
        case SpectrumType::Foreground:       msg.arg( WString::tr("ect-of-the-fore") ); break;
        case SpectrumType::SecondForeground: msg.arg( WString::tr("ect-of-the-sec") ); break;
        case SpectrumType::Background:       msg.arg( WString::tr("ect-of-the-back") ); break;
      }//switch( type )
    }else if( !msg.empty() )
    {
      msg.arg( "" );
    }//if( nfiles )
  
    /// \TODO: keep from issuing this message for ever single change!
    if( !msg.empty() )
      m_interspec->logMessage( msg, 1 );
    
    //Calling setDeviationPairMsg(...) is useless because we currently completely replace 
    //display->setDeviationPairMsg( msg );
  }else
  {
    //display->setDeviationPairMsg( "" );
  }//if( more than one gamma detector and more than one sample number ) / else
  
  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//void userChangedDeviationPair( CalDisplay *display )


void EnergyCalTool::displayedSpectrumChanged( const SpecUtils::SpectrumType type,
                                              const std::shared_ptr<SpecMeas> &meas,
                                              const std::set<int> &samples,
                                              const std::vector<std::string> &detectors )
{
  static_assert( static_cast<int>(SpecUtils::SpectrumType::Foreground) == 0, "" );
  static_assert( static_cast<int>(SpecUtils::SpectrumType::SecondForeground) == 1, "" );
  static_assert( static_cast<int>(SpecUtils::SpectrumType::Background) == 2, "" );

  const int index = static_cast<int>( type );
  assert( index >= 0 && index < 3 );
  
  if( meas != m_currentSpecMeas[index] )
  {
    //whole new file
    cout << "EnergyCalTool::displayedSpectrumChanged: new file" << endl;
    
    //We want to cache original energy calibration, if we havent already
  }else if( samples != m_currentSampleNumbers[index] )
  {
    // Just changed what was displayed
    cout << "EnergyCalTool::displayedSpectrumChanged: changed sample numbers" << endl;
  }else
  {
    //no change...
    cout << "EnergyCalTool::displayedSpectrumChanged: same file and sample numbers" << endl;
  }

  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  m_currentSpecMeas[index] = meas;
  m_currentSampleNumbers[index] = samples;
  
  refreshGuiFromFiles();
}//void displayedSpectrumChanged(...)


void EnergyCalTool::setShowNoCalInfo( const bool nocal )
{
  m_noCalTxt->setHidden( !nocal );
  m_moreActionsColumn->setHidden( nocal );
  m_applyToColumn->setHidden( nocal );
  m_detColumn->setHidden( nocal );
  m_calColumn->setHidden( nocal );
  m_peakTableColumn->setHidden( nocal );
}//void setShowNoCalInfo( const bool nocal )


void EnergyCalTool::setWasGraphicalRecal( int type, float energy )
{
  time( &m_lastGraphicalRecal );
  m_lastGraphicalRecalType = type;
  m_lastGraphicalRecalEnergy = energy;
}//void setWasGraphicalRecal( int type, double energy )


void EnergyCalTool::specTypeToDisplayForChanged()
{
  const int selectedType = m_specTypeMenu->currentIndex();
  if( selectedType < 0 || selectedType > 2 )
  {
    assert( 0 );  //shouldnt ever happen, right?
    return;
  }
  
  WMenu *detMenu = m_detectorMenu[selectedType];
  WMenuItem *detitem = detMenu->currentItem();
  if( !detitem && detMenu->count() )
    detMenu->select( 0 );
  if( detitem )
    detMenu->select( detitem );
  
  updateFitButtonStatus();
}//void specTypeToDisplayForChanged();


bool EnergyCalTool::canDoEnergyFit()
{
  // Check that "Apply Changes To" for the "Foreground" is checked, otherwise fitting makes no sense
  if( !m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked() )
    return false;
  
  // Check if there are any peaks currently showing.
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
  if( !peaks )
    return false;
  
  size_t nPeaksToUse = 0;
  for( const PeakModel::PeakShrdPtr &p : *peaks )
    nPeaksToUse += (p && p->useForEnergyCalibration());
  
  if( nPeaksToUse < 1 )
    return false;
  
  if( !m_calInfoDisplayStack )
    return false;
  
  // We are actually going to fit the coefficients for the currently showing CalDisplay, so only
  //  consult the checkboxes on that one display.
  auto caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
  if( !caldisp )
    return false;
  
  auto cal = caldisp->lastSetCalibration();
  switch( cal->type() )
  {
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      return false;
      
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      break;
  }//switch( cal->type() )
  
  const set<size_t> ordersToFit = caldisp->fitForCoefficents();
  if( ordersToFit.empty() || ordersToFit.size() > nPeaksToUse )
    return false;
  
  return true;
}//bool canDoEnergyFit()


void EnergyCalTool::fitCoefficients()
{
  try
  {
    if( !canDoEnergyFit() )
    {
      m_interspec->logMessage( WString::tr("ect-err-not-enough-peaks"), 2 );
      return;
    }//if( double check we can actually do the fit )
    
    
    // Check if there are any peaks currently showing.
    shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
    if( !peaks || peaks->empty() )  //shouldnt ever happen.
      throw runtime_error( WString::tr("ect-no-peaks").toUTF8() );
    
    // We are actually going to fit the coefficients for the currently showing CalDisplay, so only
    //  consult the checkboxes on that one display.
    auto caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
    if( !caldisp )  //shouldnt ever happen, but JIC
      throw runtime_error( "Unexpected error determining current calibration" );
    
    // The spectrum may not be displaying the detector we are currently seeing the calibration for
    //  lets make sure of this, and if so, switch to one we are displaying
    string detname = caldisp->detectorName();
    int previousSpecTypeInd = m_specTypeMenu ? m_specTypeMenu->currentIndex() : 0;
    const auto type = (previousSpecTypeInd == 1)
                          ? SpecUtils::SpectrumType::Background
                          : (previousSpecTypeInd==2
                            ? SpecUtils::SpectrumType::SecondForeground
                            : SpecUtils::SpectrumType::Foreground);
    const vector<string> displayed = m_interspec->detectorsToDisplay( type );
    
    if( std::find(begin(displayed), end(displayed), detname) == end(displayed) )
    {
      if( previousSpecTypeInd >= 0
         && previousSpecTypeInd < 3
         && m_detectorMenu[previousSpecTypeInd] )
      {
        for( WMenuItem *item : m_detectorMenu[previousSpecTypeInd]->items() )
        {
          const string thisdetname = item->text().toUTF8();
          if( std::find(begin(displayed), end(displayed), thisdetname) != end(displayed) )
          {
            item->select();
            cout << "Starting from " << caldisp << endl;
            caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
            cout << "  we moved to " << caldisp << endl;
            if( !caldisp )  //shouldnt ever happen, but JIC
              throw runtime_error( "Unexpected error determining current calibration" );
            
            detname = caldisp->detectorName();
          }//if( we found a displayed detector )
        }//for( loop over menu items to find a displayed detector )
      }else
      {
        throw runtime_error( WString::tr("ect-select-cal-of-disp-det").toUTF8() );
      }
    }//if( std::find(begin(displayed), end(displayed), detname) == end(displayed) )
    
        
    auto original_cal = caldisp->lastSetCalibration();
    
    
    
    switch( original_cal->type() )
    {
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "Unexpected calibration type from display" ); //shouldnt ever happen, but JIC
        
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        break;
    }//switch( cal->type() )
    
    if( original_cal->num_channels() < 5 )
      throw runtime_error( WString::tr("ect-not-enough-channel").toUTF8() );
    
    
    const set<size_t> orders_to_fit = caldisp->fitForCoefficents();
    if( orders_to_fit.empty() )  //shouldnt ever happen
      throw runtime_error( WString::tr("ect-no-coeff-selected").toUTF8() );
    
    
    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contain the error of the fit mean for that peak
    vector<EnergyCal::RecalPeakInfo> peakInfos;
    
    for( const auto &peakptr : *peaks )
    {
      if( !peakptr )  //shouldnt be necassary, but JIC
        continue;
      const PeakDef &peak = *peakptr;
      
      if( !peak.useForEnergyCalibration() )
        continue;
      
      const double wantedEnergy = peak.gammaParticleEnergy();
      
      EnergyCal::RecalPeakInfo peakInfo;
      peakInfo.peakMean = peak.mean();
      // Clamp to peak mean uncertainty to be at least 0.25 keV.  This is an arbitrary decision, but
      //  motivated by not wanting a single peak to way, way, dominate the other peaks, when it is
      //  likely non-linearities in the detector may actually dominate the effects
      peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
      if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
        peakInfo.peakMeanUncert = 0.5;
      
      peakInfo.photopeakEnergy = wantedEnergy;
      peakInfo.peakMeanBinNumber = original_cal->channel_for_energy( peak.mean() );
      
      peakInfos.push_back( peakInfo );
    }//for( int col = 0; col < numModelCol; ++col )
    
    if( orders_to_fit.size() > peakInfos.size() )
      throw runtime_error( WString::tr("ect-err-not-enough-peaks").toUTF8() );
    
    auto answer = make_shared<SpecUtils::EnergyCalibration>();
    
    const size_t eqn_order = std::max( original_cal->coefficients().size(), (*orders_to_fit.rbegin()) + 1 );
    const size_t nchannel = original_cal->num_channels();
    const auto &devpairs = original_cal->deviation_pairs();
    
    double chi2 = -999;
    
    try
    {
      vector<bool> fitfor( eqn_order, false );
      
      for( auto order : orders_to_fit )
      {
        assert( order < fitfor.size() );
        fitfor[order] = true;
      }
      
      vector<float> coefficent_uncerts;
      vector<float> coefficents = original_cal->coefficients();
      if( coefficents.size() < eqn_order )
        coefficents.resize( eqn_order, 0.0f );
      
      switch ( original_cal->type() )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          chi2 = EnergyCal::fit_energy_cal_poly( peakInfos, fitfor, nchannel, devpairs,
                                                coefficents, coefficent_uncerts );
          answer->set_polynomial( nchannel, coefficents, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          chi2 = EnergyCal::fit_energy_cal_frf( peakInfos, fitfor, nchannel, devpairs,
                                               coefficents, coefficent_uncerts );
          answer->set_full_range_fraction( nchannel, coefficents, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
        case SpecUtils::EnergyCalType::InvalidEquationType:
          throw runtime_error( "Didnt expect lower channel or invalid eqn type." );
          break;
      }//switch ( original_cal->type() )
      
      //Print some developer info to terminal
      stringstream msg;
      msg << "\nfit_energy_cal_poly gave chi2=" << chi2 << " with coefs={";
      for( size_t i = 0; i < coefficents.size(); ++i )
        msg << coefficents[i] << "+-" << coefficent_uncerts[i] << ", ";
      msg << "}\n";
      cout << msg.str() << endl;
      
    }catch( std::exception &e )
    {
      cerr << "fit_energy_cal_poly threw: " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      char buffer[512] = { '\0' };
      snprintf( buffer, sizeof(buffer)-1, "fit_energy_cal_poly threw: %s", e.what() );
      log_developer_error( __func__, buffer );
#endif
    }//try / catch fit for coefficents using least linear squares
    
    
    if( !answer->valid() )
    {
      vector<bool> fitfor( eqn_order, false );
      
      for( auto order : orders_to_fit )
        fitfor[order] = true;
      
      vector<float> calib_coefs = original_cal->coefficients();
      if( calib_coefs.size() < eqn_order )
        calib_coefs.resize( eqn_order, 0.0f );
      
      std::string warning_msg;
      std::vector<float> coefs, coefs_uncert;
      chi2 = EnergyCal::fit_energy_cal_iterative( peakInfos, nchannel, original_cal->type(), fitfor,
                                          calib_coefs, devpairs, coefs, coefs_uncert, warning_msg );
      
      if( warning_msg.size() )
        m_interspec->logMessage( warning_msg, 3 );
      
      switch ( original_cal->type() )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          answer->set_polynomial( nchannel, coefs, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          answer->set_full_range_fraction( nchannel, coefs, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
        case SpecUtils::EnergyCalType::InvalidEquationType:
          assert( 0 );
          break;
      }//switch ( original_cal->type() )
    }//if( !fit_coefs )
    
    if( !answer->valid() )
      throw runtime_error( WString::tr("ect-fail-min").toUTF8() );
    
  
    if( !answer || !answer->valid() )
      return;
    
    const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
    applyCalChange( original_cal, answer, changemeas, false );
    
    
    //To show Chi2 in the message, uncomment out this next section
    /*
    if( peakInfos.size() > orders_to_fit.size() )
    {
      double dof = peakInfos.size() - 1;
      dof -= orders_to_fit.size();
      dof = (dof < 1) ? 1.0 : dof;
      
      char buffer[64];
      snprintf( buffer, sizeof(buffer), " &chi;&sup2;/dof=%.2g", (chi2/dof) );
      msg += buffer;
    }
    */
    
    m_interspec->logMessage( WString::tr("ect-fit-successful"), 1 );
  }catch( std::exception &e )
  {
    WString msg = WString::tr("ect-fail-fit").arg( e.what() );
    cerr << "EnergyCalTool::fitCoefficients():\n\tCaught: " << msg.toUTF8() << endl;
    m_interspec->logMessage( msg, 3 );
  }//try / catch
}//void fitCoefficients()


void EnergyCalTool::updateFitButtonStatus()
{
  const bool canFit = canDoEnergyFit();
  
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  for( EnergyCalImp::CalDisplay *disp : calDisplays() )
    disp->setFitButtonEnabled( canFit );
#else
  if( canFit != m_fitCalBtn->isEnabled() )
    m_fitCalBtn->setDisabled( !canFit );
#endif
}//void updateFitButtonStatus()


#if( !IMP_CALp_BTN_NEAR_COEFS )
void EnergyCalTool::updateCALpButtonsStatus()
{
  shared_ptr<const SpecUtils::Measurement> foregrnd
                           = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  const bool showDownload = (foregrnd && foregrnd->energy_calibration()
                              && foregrnd->energy_calibration()->valid());
  m_downloadCALp->setHidden( !showDownload );
  
  const bool showUpload = (foregrnd && foregrnd->num_gamma_channels());
  m_uploadCALp->setHidden( !showUpload );
}//void updateCALpButtonsStatus()
#endif


void EnergyCalTool::displayedSpecChangedCallback( const SpecUtils::SpectrumType,
                                                  const std::shared_ptr<SpecMeas>,
                                                  const std::set<int>,
                                                  const std::vector<std::string> )
{
  /// \TODO: set the various m_applyToCbs if it is a new spectrum being shown.
  /// \TODO: if this is the first time seeing a SpecMeas, cache all of its energy calibration
  ///        information
  
  // \TODO: we could maybe save a little time by inspecting what was changed, but the added
  //        complexity probably isnt worth it, so we'll skip this.
  refreshGuiFromFiles();
}//void displayedSpecChangedCallback(...)


void EnergyCalTool::refreshGuiFromFiles()
{
  m_renderFlags |= EnergyCalToolRenderFlags::FullGuiUpdate;
  scheduleRender();
}//void refreshGuiFromFiles()


void EnergyCalTool::handleGraphicalRecalRequest( double xstart, double xfinish )
{
  try
  {
    auto foreground = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::EnergyCalibration> energycal = foreground
                                                               ? foreground->energy_calibration()
                                                               : nullptr;
    if( !energycal || !energycal->valid()
        || (energycal->type() == SpecUtils::EnergyCalType::LowerChannelEdge) )
      return;
  
    if( m_graphicalRecal )
    {
      m_graphicalRecal->setEnergies( xstart, xfinish );
    }else
    {
      m_graphicalRecal = new EnergyCalGraphicalConfirm( xstart, xfinish, this,
                                    m_lastGraphicalRecal,
                                    EnergyCalGraphicalConfirm::RecalTypes(m_lastGraphicalRecalType),
                                    m_lastGraphicalRecalEnergy );
    
      m_graphicalRecal->finished().connect( this, &EnergyCalTool::deleteGraphicalRecalConfirmWindow );
    }
  }catch( std::runtime_error & )
  {
    m_interspec->logMessage( "Internal error doing graphical recal; sorry :(", 3 );
  }
}//void handleGraphicalRecalRequest( double xstart, double xfinish )


void EnergyCalTool::deleteGraphicalRecalConfirmWindow()
{
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal );
    m_graphicalRecal = nullptr;
  }//if( m_graphicalRecal )
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  if( showToolTips )
  {
    m_interspec->logMessage( WString::tr("ect-del-graphical-msg"), 1 );
  }//if( showToolTips )
}//void deleteGraphicalRecalConfirmWindow()


string EnergyCalTool::applyToSummaryTxt() const
{
  string answer;
  
  if( m_applyToColumn->isHidden() )
    return answer;
  
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    Wt::WCheckBox *cb = m_applyToCbs[index];
    const auto cbparent = cb->parent();
    assert( cbparent );
    assert( dynamic_cast<WContainerWidget *>(cbparent) );
    
    if( cbparent->isHidden() || !cb->isChecked() )
      continue;
    
    if( !answer.empty() )
      answer += ", ";
    
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:         answer += WString::tr("foreground").toUTF8(); break;
      case ApplyToCbIndex::ApplyToBackground:         answer += WString::tr("background").toUTF8(); break;
      case ApplyToCbIndex::ApplyToSecondary:          answer += WString::tr("secondary").toUTF8();  break;
      case ApplyToCbIndex::ApplyToDisplayedDetectors: answer += WString::tr("ect-lc-disp-dets").toUTF8(); break;
      case ApplyToCbIndex::ApplyToAllDetectors:       answer += WString::tr("ect-lc-all-dets").toUTF8(); break;
      case ApplyToCbIndex::ApplyToDisplayedSamples:   answer += WString::tr("ect-lc-disp-samples").toUTF8(); break;
      case ApplyToCbIndex::ApplyToAllSamples:         answer += WString::tr("ect-lc-all-samples").toUTF8(); break;
      case ApplyToCbIndex::NumApplyToCbIndex:
        break;
    }//switch( index )
  }//for( loop over ApplyToCbIndex )
  
  return answer;
}//string applyToSummaryTxt() const


void EnergyCalTool::doRefreshFromFiles()
{
  //Labels for horizontal labels when you have multiple spectra shown, and at least one of them has
  //  more than one detectors
  const char * const spec_type_labels[3] = {"ect-short-fore","ect-short-back","ect-short-secondary"};
  const char * const spec_type_labels_vert[3] = {"Foreground", "Background", "Secondary"};
  
  string prevdet[3];
  int previousSpecInd = m_specTypeMenu ? m_specTypeMenu->currentIndex() : 0;
  
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  shared_ptr<const SpecMeas> specfiles[3];
  for( int i = 0; i < 3; ++i )
    specfiles[i] = m_interspec->measurment( spectypes[i] );
  
  set<string> specdetnames[3]; //Just the names of gamma detectors with at least 4 channels
  
/*
  //Delete calibration contents and menu items - we will add in all the current ones below.
  //  - this appears to not work well - the stacks need replacing...
  for( int i = 0; i < 3; ++i)
  {
    Wt::WMenu *menu = m_detectorMenu[i];
    
    if( menu->currentItem() )
      prevdet[i] = menu->currentItem()->text().toUTF8();
    
    for( WMenuItem *item : menu->items() )
    {
      WWidget *content = item->contents();
      assert( content );
      menu->removeItem( item );
      
      delete item;
      delete content;
    }//
  }//for( Wt::WMenu *menu : m_detectorMenu )
*/
  
  const bool isWide = (m_tallLayoutContent ? false : true);
  
  //If each spectrum file has only a single detector, then instead of having vertical menu display
  //  detector name, we will have "For.", "Back", "Sec."
  bool specTypeInForgrndMenu = true;
  
  //If we try to re-use the menu and stacks, for some reason the calibration coefficents wont show
  //  up if we alter which background/secondary spectra are showing... not sure why, but for the
  //  moment we'll just re-create the menus and stacks... not great, but works, for the moment.
  bool needStackRefresh = false;
  
  // Get the detector names, for the displayed sample numbers, for each spectrum type
  set<string> disp_det_names[3];
  
  int nFilesWithCalInfo = 0;
  
  // We want to preserve wich "Fit" check boxes are set.
  map< pair<SpecUtils::SpectrumType,string>, set<size_t> > set_fit_for_cbs;
  
  
  {//begin code-block to see if we need to refresh stack
    
    //TODO: need to check that isWide hasnt changed, and if it has set needStackRefresh=true
    
    //Check if we can put "For.", "Back", "Sec." on the vertical menu instead of detector names
    for( int i = 0; i < 3; ++i )
    {
      disp_det_names[i] = gammaDetectorsForDisplayedSamples( spectypes[i] );
      specTypeInForgrndMenu = (specTypeInForgrndMenu && (disp_det_names[i].size() <= 1));
      
      if( !disp_det_names[i].empty() )
        nFilesWithCalInfo += 1;
    }//for( int i = 0; i < 3; ++i )
    
    if( !m_specTypeMenu )
    {
      //cout << "needStackRefresh: m_specTypeMenu == nullptr" << endl;
      needStackRefresh = true;
    }
    
    for( int i = 0; !needStackRefresh && (i < 3); ++i )
    {
      WMenuItem *typeItem = m_specTypeMenu->itemAt(i);
      if( !typeItem )
      {
        //cout << "needStackRefresh: !typeItem" << endl;
        needStackRefresh = true;
        continue;
      }
      
      WMenu *detMenu = m_detectorMenu[(specTypeInForgrndMenu ? 0 : i)];
      
      if( !detMenu )
      {
        //cout << "needStackRefresh: !detMenu" << endl;
        needStackRefresh = true;
        continue;
      }
      
      if( !specTypeInForgrndMenu && ((!detMenu) != (!specfiles[i])) )
      {
        //cout << "needStackRefresh: (!detMenu != !specfiles[i])" << endl;
        needStackRefresh = true;
        continue;
      }
      
      if( specTypeInForgrndMenu )
      {
        if( specfiles[i] )
        {
          //Make sure one of the widgets have spec_type_labels_vert[i] in it
          needStackRefresh = true;
          for( int w = 0; needStackRefresh && (w < detMenu->count()); ++w )
          {
            auto item = detMenu->itemAt(w);
            needStackRefresh = !(item && (item->text() == WString::tr(spec_type_labels_vert[i])));
          }
          
          //if( needStackRefresh )
          //  cout << "needStackRefresh: Did not have a menu entry for specfile[" << i << "]" << endl;
        }else
        {
          //Make sure none of the widgets have spec_type_labels_vert[i] in them
          for( int w = 0; !needStackRefresh && (w < detMenu->count()); ++w )
          {
            auto item = detMenu->itemAt(w);
            needStackRefresh = (!item || (item->text() == WString::tr(spec_type_labels_vert[i])));
          }
          
          //if( needStackRefresh )
          //  cout << "needStackRefresh: For empty place " << i << " had a menu entry, but no spectrum file" << endl;
        }//if( specfiles[i] ) / else
      }else
      {
        assert( m_specTypeMenu );
        WMenuItem *item = m_specTypeMenu->itemAt(i);
        assert( item );
     
        if( !item )
        {
          //Shouldnt ever get here I think
          needStackRefresh = true;
          continue;
        }//if( !item )
        
        if( !specfiles[i] )
        {
          //If there is no spectrum file for index i, item should be hidden, and if not need a refresh
          needStackRefresh = !item->isHidden();
          //if( needStackRefresh )
          //  cout << "needStackRefresh: item is not hidden for missing spectrum " << i << endl;
          continue;
        }//if( !specfiles[i] )
        
        if( item->isHidden() )
        {
          //If item is hidden, but here we do have a spectrum file, we need a refresh
          //cout << "needStackRefresh: item is hidden for " << i << endl;
          needStackRefresh = true;
          continue;
        }//if( item->isHidden() )
        
        //if( needStackRefresh )
        //  cout << "needStackRefresh: New det names dont equal old for type=" << i << endl;
      }//if( specTypeInForgrndMenu ) / else
      
      
      //Need to check all the detectors are the same names
      set<string> detsInMenu;
      for( int j = 0; j < detMenu->count(); ++j )
      {
        auto detitem = detMenu->itemAt(j);
        if( detitem )
          detsInMenu.insert( detitem->text().toUTF8() );
      }//for( int j = 0; j < detMenu->count(); ++j )
      
      needStackRefresh = (detsInMenu != disp_det_names[i]);
    }//for( int i = 0; i < 3; ++i )
  }//end code-block to see if we need to refresh stack
  
  //Dont show spectype menu (the vertical "For.", "Back", "Sec." menu), if we dont need to
  const bool hideSpecType = ( specTypeInForgrndMenu
                             || (nFilesWithCalInfo < 2)
                             || ( (!specfiles[1] || (specfiles[0]==specfiles[1]))
                                 && (!specfiles[2] || (specfiles[0]==specfiles[2]))) );
  
  if( !m_specTypeMenu || (m_specTypeMenu->isHidden() != hideSpecType) )
  {
    //cout << "needStackRefresh: m_specTypeMenu->isHidden() != hideSpecType" << endl;
    needStackRefresh = true;
  }
  
  //cout << "needStackRefresh=" << needStackRefresh << endl;
  
  // TODO: instead of the above logic to catch when we need to refresh (\e.g., create all new)
  //       widgets, should combine it with the logic to create new widgets, but only delte or create
  
  if( needStackRefresh )
  {
    for( int i = 0; i < 3; ++i )
    {
      if( m_detectorMenu[i] )
      {
        for( WMenuItem *item : m_detectorMenu[i]->items() )
        {
          /// \TODO: the item text isnt necassarily the detector name - when specTypeInForgrndMenu
          ///        is true, this breaks down - need to fix this
          const string detname = item->text().toUTF8();
          auto display = dynamic_cast<EnergyCalImp::CalDisplay *>( item->contents() );
          if( display )
            set_fit_for_cbs[{spectypes[i],detname}] = display->fitForCoefficents();
          else
            cerr << "Unexpected widget type as a sub-menu!" << endl;
        }//for( loop over menu types )
        
        if( m_detectorMenu[i]->currentItem() )
          prevdet[i] = m_detectorMenu[i]->currentItem()->text().toUTF8();
        delete m_detectorMenu[i];
      }
      m_detectorMenu[i] = nullptr;
    }//for( int i = 0; i < 3; ++i )
    
    delete m_specTypeMenuStack;
    delete m_specTypeMenu;
    m_specTypeMenuStack = nullptr;
    m_specTypeMenu = nullptr;
    
    WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
    
    m_specTypeMenuStack = new WStackedWidget();
    m_specTypeMenuStack->addStyleClass( "CalSpecStack" );
    m_specTypeMenuStack->setTransitionAnimation( animation );
    
    auto callayout = dynamic_cast<WGridLayout *>( m_calColumn->layout() );
    assert( callayout );
    
    m_specTypeMenu = new WMenu( m_specTypeMenuStack );
    m_specTypeMenu->addStyleClass( "CalSpecMenu" );
    m_specTypeMenu->itemSelected().connect( this, &EnergyCalTool::specTypeToDisplayForChanged );
    m_detColLayout->addWidget( m_specTypeMenu, 1, 0 );  
    m_detColLayout->addWidget( m_specTypeMenuStack, 2, 0 );
    
    if( m_calInfoDisplayStack )
      delete m_calInfoDisplayStack;
    m_calInfoDisplayStack = new WStackedWidget();
    m_calInfoDisplayStack->addStyleClass( "ToolTabTitledColumnContent CalStack" );
    m_calInfoDisplayStack->setTransitionAnimation( animation );
    callayout->addWidget( m_calInfoDisplayStack, 1, 1 );
    
    /// \TODO: only create these menus when actually needed, so we wont need to
    for( int i = 0; i < 3; ++i )
    {
      if( !specfiles[i] )
      {
        //Add a dummy entry into the menu or else the 'm_specTypeMenu->itemAt(i)' call below will
        //  segfault or not necassarily give the wanted answer.
        WContainerWidget *detMenuDiv = new WContainerWidget();
        WMenuItem *item = m_specTypeMenu->addItem( "", detMenuDiv );
        item->setHidden( true );
        continue;
      }
      
      WContainerWidget *detMenuDiv = new WContainerWidget();  //this holds the WMenu for this SpecFile
      detMenuDiv->addStyleClass( "DetMenuDiv" );
      
      WMenuItem *item = m_specTypeMenu->addItem( WString::tr(spec_type_labels[i]), detMenuDiv, WMenuItem::LoadPolicy::PreLoading );
      //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
      item->clicked().connect( boost::bind(&WMenuItem::select, item) );
      
      m_detectorMenu[i] = new WMenu( m_calInfoDisplayStack, detMenuDiv );
      m_detectorMenu[i]->addStyleClass( "VerticalNavMenu HeavyNavMenu DetCalMenu" );
      
      m_detectorMenu[i]->itemSelected().connect( this, &EnergyCalTool::updateFitButtonStatus );
    }//for( int i = 0; i < 3; ++i )
  }//if( needStackRefresh )

#if( !IMP_CALp_BTN_NEAR_COEFS )
  updateCALpButtonsStatus();
#endif
  
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  const bool canFitCeofs = canDoEnergyFit();
#endif
  
  if( !specfiles[0] )
  {
    setShowNoCalInfo( true );
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
    for( EnergyCalImp::CalDisplay *disp : calDisplays() )
      disp->setFitButtonEnabled( false );
#else
    m_fitCalBtn->disable();
#endif
    return;
  }
  
  if( previousSpecInd < 0 ||  previousSpecInd > 2
     || (previousSpecInd == 1 && !specfiles[1])
     || (previousSpecInd == 2 && !specfiles[2])  )
  {
    previousSpecInd = 0;
  }
  
  
  bool selectedDetToShowCalFor = false;
  bool hasFRFCal = false, hasPolyCal = false, hasLowerChanCal = false;
  
  
  for( int i = 0; i < 3; ++i )
  {
    if( !specfiles[i] )
      continue;
    
    Wt::WMenu *detMenu = m_detectorMenu[(specTypeInForgrndMenu ? 0 : i)];
    if( !detMenu )
      continue;
    
    WMenuItem *specItem = m_specTypeMenu->itemAt(i);
    assert( specItem );
    
    const SpecUtils::SpectrumType type = spectypes[i];
    shared_ptr<const SpecMeas> meas = specfiles[i];
    assert( meas );
    
    const set<int> &samples = m_interspec->displayedSamples(type);
    const set<string> &detectors = disp_det_names[i];
    
    specdetnames[i] = detectors;
    
    if( detectors.empty() )
    {
      if( m_specTypeMenu->currentIndex() == i )
        m_specTypeMenu->select(0);
      specItem->setHidden( true );
      continue;
    }
    

    assert( specItem );
    specItem->setHidden( false );
    
    for( const string &detname : detectors )
    {
      for( const int sample : samples )
      {
        auto m = meas->measurement( sample, detname );
        if( !m || (m->num_gamma_channels() <= 4) )
          continue;
        
        shared_ptr<const SpecUtils::EnergyCalibration> energycal = m->energy_calibration();
        
        WString displayname = WString::fromUTF8( detname );
        
        if( specTypeInForgrndMenu )
        {
          /// \TODO: when specTypeInForgrndMenu is true, we may not match detector name because
          ///        we are getting the menu item text above and assuming the detector name.
          ///        should fix this.
          displayname = WString::tr(spec_type_labels_vert[i]);
        }//if( specTypeInForgrndMenu )
        
        WMenuItem *item = nullptr;
        if( !needStackRefresh )
        {
          for( int j = 0; !item && (j < detMenu->count()); ++j )
          {
            auto jitem = detMenu->itemAt(j);
            if( jitem && (jitem->text() == WString(displayname)) )
              item = jitem;
          }
        }//if( !needStackRefresh )
        
        if( item )
        {
          WWidget *ww = item->contents();
          EnergyCalImp::CalDisplay *calcontent = dynamic_cast<EnergyCalImp::CalDisplay *>(ww);
          assert( calcontent );
          if( calcontent )
            calcontent->updateToGui( energycal );
          else
            item = nullptr;
        }//if( item )
        
        if( !item )
        {
          auto calcontent = new EnergyCalImp::CalDisplay( this, type, detname, isWide );
          item = detMenu->addItem( displayname, calcontent, WMenuItem::LoadPolicy::PreLoading );
          //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
          item->clicked().connect( boost::bind(&WMenuItem::select, item) );
          
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
          calcontent->setFitButtonEnabled( canFitCeofs );
          calcontent->doFitCoeffs().connect( this, &EnergyCalTool::fitCoefficients );
#endif
          
          calcontent->updateToGui( energycal );
          
          const auto fitfor_iter = set_fit_for_cbs.find( {type,displayname.toUTF8()} );
          if( fitfor_iter != end(set_fit_for_cbs) )
            calcontent->setFitFor( fitfor_iter->second );
          
          if( (displayname == prevdet[i]) && (i == previousSpecInd) )
          {
            m_specTypeMenu->select( previousSpecInd );
            detMenu->select( item );
            selectedDetToShowCalFor = true;
          }
        }//if( !item )
        
        if( energycal )
        {
          switch( energycal->type() )
          {
            case SpecUtils::EnergyCalType::Polynomial:
            case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
              hasPolyCal = true;
              break;
              
            case SpecUtils::EnergyCalType::FullRangeFraction:
              hasFRFCal = true;
              break;
              
            case SpecUtils::EnergyCalType::LowerChannelEdge:
              hasLowerChanCal = true;
              break;
            
            case SpecUtils::EnergyCalType::InvalidEquationType:
              break;
          }//switch( energycal->type() )
        }//if( energycal )
      
        break;
      }
    }//for( const string &detname : detectors )
  }//for( int i = 0; i < 3; ++i )
  
  if( needStackRefresh && !selectedDetToShowCalFor && m_detectorMenu[0] && m_detectorMenu[0]->count() )
  {
    m_specTypeMenu->select( 0 );
    m_detectorMenu[0]->select( 0 );
  }
  
  setShowNoCalInfo( !nFilesWithCalInfo );
  m_specTypeMenu->setHidden( hideSpecType );
  
  bool anyApplyToCbShown = false;
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    Wt::WCheckBox *cb = m_applyToCbs[index];
    const auto cbparent = cb->parent();
    assert( cbparent );
    assert( dynamic_cast<WContainerWidget *>(cbparent) );
    
    bool hideRow = false;
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:
        hideRow = (!specfiles[0] || (!specfiles[1] && !specfiles[2]));
        break;
        
      case ApplyToCbIndex::ApplyToBackground:
        // We'll ignore case where foreground and background is same {SpecFile, Samples, Detectors}
        hideRow = !specfiles[1];
        break;
        
      case ApplyToCbIndex::ApplyToSecondary:
        // We'll ignore case where secondary and background is same {SpecFile, Samples, Detectors}
        hideRow = !specfiles[2];
        break;
        
      case ApplyToCbIndex::ApplyToDisplayedDetectors:
      {
        bool displayingAll = true;
        for( int i = 0; displayingAll && (i < 3); ++i )
        {
          const auto meas = specfiles[i];
          if( !meas )
            continue;
          
          const auto type = spectypes[i];
          const vector<string> displayed = m_interspec->detectorsToDisplay( type );
          for( const auto &name : meas->gamma_detector_names() )
          {
            if( std::find( begin(displayed), end(displayed), name ) == end(displayed) )
              displayingAll = false;
          }//for( const auto &name : meas->gamma_detector_names() )
        }//for( loop over the types of spectrum files )
        
        hideRow = displayingAll;
        if( displayingAll )
          cb->setChecked( false );
        
        break;
      }//case ApplyToCbIndex::ApplyToDisplayedDetectors:
        
      case ApplyToCbIndex::ApplyToAllDetectors:
      {
        static_assert( ApplyToCbIndex::ApplyToDisplayedDetectors
                       < ApplyToCbIndex::ApplyToAllDetectors, "" );
        
        const WCheckBox *dispDetsCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors];
        assert( dispDetsCb->parent() );
        
        const bool dispDetsHid = dispDetsCb->parent()->isHidden();
        
        hideRow = dispDetsHid;
        if( dispDetsHid )
          cb->setChecked( true );
        
        if( dispDetsCb->isChecked() )
          cb->setChecked( false );
        
        break;
      }//case ApplyToAllDetectors:
        
      case ApplyToCbIndex::ApplyToDisplayedSamples:
      {
        // We want to check that for each displayed unique SpecMeas we are displaying all samples.
        /// \TODO: avoid allocating all the set<int>'s below, and then the looping threw every value
        ///        Can probably easily use a hyristic to avoid most of the time, or change how
        ///        things are tracked to avoid probably all the time (at the cost of adding
        ///        complexity to the code)
        map<shared_ptr<const SpecMeas>,set<int>> undisplayed;
        for( int i = 0; i < 3; ++i )
        {
          const auto meas = specfiles[i];
          if( !meas )
            continue;
          
          if( !undisplayed.count(meas) )
            undisplayed[meas] = meas->sample_numbers();
          
          set<int> &undispsamples = undisplayed[meas];
          for( const int sample : m_interspec->displayedSamples(spectypes[i]) )
            undispsamples.erase( sample );
        }//for( loop over the types of spectrum files )
        
        bool displayingAll = true;
        for( const auto &p : undisplayed )
          displayingAll = (displayingAll && p.second.empty());
        
        hideRow = displayingAll;
        if( displayingAll )
          cb->setChecked( false );
        
        break;
      }//case ApplyToDisplayedSamples:
        
      case ApplyToCbIndex::ApplyToAllSamples:
      {
        static_assert( ApplyToCbIndex::ApplyToDisplayedSamples
                         < ApplyToCbIndex::ApplyToAllSamples, "" );
        
        const WCheckBox *dispSamplesCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedSamples];
        assert( dispSamplesCb->parent() );
        const bool dispSamplesHid = dispSamplesCb->parent()->isHidden();
        
        hideRow = dispSamplesHid;
        if( dispSamplesHid )
          cb->setChecked( true );
        
        if( dispSamplesCb->isChecked() )
          cb->setChecked( false );
        
        break;
      }//case ApplyToCbIndex::ApplyToAllSamples:
        
        
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    cbparent->setHidden( hideRow );
    
    if( !hideRow )
      anyApplyToCbShown = true;
  }//for( loop over ApplyToCbIndex )
  
  m_applyToColumn->setHidden( !anyApplyToCbShown );
  
  bool hideDetCol = true;
  if( specfiles[0] && specdetnames[0].size() > 1 )
    hideDetCol = false;
  if( specfiles[1] && (specfiles[0] != specfiles[1]) )
    hideDetCol = false;
  if( specfiles[2] && (specfiles[0] != specfiles[2]) )
    hideDetCol = false;
  
  m_detColumn->setHidden( hideDetCol );
  
  for( MoreActionsIndex index = MoreActionsIndex(0);
      index < MoreActionsIndex::NumMoreActionsIndex;
      index = MoreActionsIndex(static_cast<int>(index) + 1) )
  {
    Wt::WAnchor *anchor = m_moreActions[static_cast<int>(index)];
    assert( anchor );
    auto aparent = anchor->parent();
    assert( dynamic_cast<WContainerWidget *>(aparent) );
    
    switch( index )
    {
      case MoreActionsIndex::Linearize:
      case MoreActionsIndex::Truncate:
      case MoreActionsIndex::CombineChannels:
        aparent->setHidden( !specfiles[0] );
        break;
        
      case MoreActionsIndex::ConvertToFrf:
        aparent->setHidden( !hasPolyCal );
        break;
        
      case MoreActionsIndex::ConvertToPoly:
        aparent->setHidden( !(hasFRFCal || hasLowerChanCal) );
        break;
        
      case MoreActionsIndex::MultipleFilesCal:
      {
        SpecMeasManager *manager = m_interspec->fileManager();
        SpectraFileModel *fmodel = manager ? manager->model() : nullptr;
        const int nfiles = fmodel ? fmodel->rowCount() : 0;
        int nRecordsWithPeaks = 0;
        for( int row = 0; row < nfiles && (nRecordsWithPeaks < 2); ++row )
        {
          shared_ptr<SpectraFileHeader> header = fmodel->fileHeader( row );
          if( !header )
            continue;
          
          int nsamples = header->numSamples();
          shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
          if( meas )
            nRecordsWithPeaks += meas->sampleNumsWithPeaks().size();
          else
            nRecordsWithPeaks += nsamples; //not worth readin file from disk, so we'll be hopeful
        }
        
        aparent->setHidden( nRecordsWithPeaks < 2 );
        break;
      }//case MoreActionsIndex::MultipleFilesCal:
        
      case MoreActionsIndex::NumMoreActionsIndex:
        break;
    }//switch( index )
  }//for( loop over
  
  // Update the "Fit Coeffs" button to be enabled/disabled
  updateFitButtonStatus();
  
  //const int currentwidget = m_detectorMenu[0]->contentsStack()->currentIndex();
  //cout << "currentwidget=" << currentwidget << endl;
}//void doRefreshFromFiles()



void EnergyCalTool::moreActionBtnClicked( const MoreActionsIndex index )
{
  const vector<MeasToApplyCoefChangeTo> measToChange = measurementsToApplyCoeffChangeTo();

  if( m_addActionWindow )
  {
    AuxWindow::deleteAuxWindow( m_addActionWindow );
    m_addActionWindow = nullptr;
  }
  
  m_addActionWindow = new EnergyCalAddActionsWindow( index, measToChange, this );
  m_addActionWindow->finished().connect( this, &EnergyCalTool::cancelMoreActionWindow );
  
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( !undoManager )
    return;
  
  auto undo = [](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( tool )
      tool->cancelMoreActionWindow();
  };
  
  auto redo = [index](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( tool )
      tool->moreActionBtnClicked(index);
  };
  
  undoManager->addUndoRedoStep( undo, redo, "Show additional energy cal tool.");
}//void moreActionBtnClicked( const MoreActionsIndex index )


void EnergyCalTool::cancelMoreActionWindow()
{
  if( m_addActionWindow )
  {
    AuxWindow::deleteAuxWindow( m_addActionWindow );
    m_addActionWindow = nullptr;
  }
}//void cancelMoreActionWindow()


void EnergyCalTool::render( Wt::WFlags<Wt::RenderFlag> flags)
{
  //flags.testFlag(RenderFlag::RenderFull) will only be true on initial rending of widget, and
  //  after that only the RenderFlag::RenderUpdate flag will be set
  
  if( flags.testFlag(Wt::RenderFlag::RenderFull)
      || m_renderFlags.testFlag(EnergyCalToolRenderFlags::FullGuiUpdate) )
  {
    doRefreshFromFiles();
    m_renderFlags.clear( EnergyCalToolRenderFlags::FullGuiUpdate );
  }
  
  WContainerWidget::render(flags);
}//void render( Wt::WFlags<Wt::RenderFlag> flags)


void EnergyCalTool::applyToCbChanged( const EnergyCalTool::ApplyToCbIndex index )
{
  // We only get here if the user checked/unchecked a checkbox, or in a undo/redo step.
  assert( index >= 0 && index <= EnergyCalTool::NumApplyToCbIndex );
  WCheckBox *cb = m_applyToCbs[index];
  const bool isChecked = cb->isChecked();
  
  // Grab the starting state of all checkboxed
  bool startingState[ApplyToCbIndex::NumApplyToCbIndex];
  for( ApplyToCbIndex i = ApplyToCbIndex(0);
      i < ApplyToCbIndex::NumApplyToCbIndex;
      i = ApplyToCbIndex(i + 1) )
  {
    startingState[i] = m_applyToCbs[i]->isChecked();
  }
  
  // Assume `index` was actually not what it is now.
  startingState[index] = !startingState[index];
  
  
  switch( index )
  {
    case EnergyCalTool::ApplyToForeground:
      updateFitButtonStatus();
      break;
      
    case EnergyCalTool::ApplyToBackground:
    case EnergyCalTool::ApplyToSecondary:
      break;
    
    case EnergyCalTool::ApplyToDisplayedDetectors:
      m_applyToCbs[EnergyCalTool::ApplyToAllDetectors]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToAllDetectors:
      m_applyToCbs[EnergyCalTool::ApplyToDisplayedDetectors]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToDisplayedSamples:
      m_applyToCbs[EnergyCalTool::ApplyToAllSamples]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToAllSamples:
      m_applyToCbs[EnergyCalTool::ApplyToDisplayedSamples]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::NumApplyToCbIndex:
      break;
  }//switch( index )
  
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( !undoManager || undoManager->isInUndoOrRedo() )
    return;
  
  
  bool finalState[ApplyToCbIndex::NumApplyToCbIndex];
  for( ApplyToCbIndex i = ApplyToCbIndex(0);
      i < ApplyToCbIndex::NumApplyToCbIndex;
      i = ApplyToCbIndex(i + 1) )
  {
    finalState[i] = m_applyToCbs[i]->isChecked();
  }
  
  
  auto undo = [index,isChecked,startingState](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( !tool )
      return;
    
    for( ApplyToCbIndex i = ApplyToCbIndex(0);
        i < ApplyToCbIndex::NumApplyToCbIndex;
        i = ApplyToCbIndex(i + 1) )
    {
      tool->m_applyToCbs[i]->setChecked( startingState[i] );
    }
    
    tool->applyToCbChanged( index );
  };
  
  auto redo = [index,isChecked,finalState](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( !tool )
      return;
    
    for( ApplyToCbIndex i = ApplyToCbIndex(0);
        i < ApplyToCbIndex::NumApplyToCbIndex;
        i = ApplyToCbIndex(i + 1) )
    {
      tool->m_applyToCbs[i]->setChecked( finalState[i] );
    }
    
    tool->applyToCbChanged( index );
  };
  
  const string label = cb->text().toUTF8();
  undoManager->addUndoRedoStep( undo, redo, "Change " + label );
}//void applyToCbChanged( const ApplyToCbIndex index )
