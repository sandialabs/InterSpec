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

#include <atomic>
#include <string>
#include <iostream>

#include <Wt/WServer>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/UndoRedoManager.h"


using namespace std;

namespace
{
  /** The number of undo/redo steps to keep in memory; There may be `ns_max_steps + ns_nsteps_histerious`
   steps in memory before cleanup is triggered.
   
   This value applies to all sessions - in the future we could add a member variable to #UndoRedoManager to allow more
   fine-grained control.
   */
  std::atomic<int> ns_max_steps( 250 );
  
  /** The amount of steps we can go over #ns_max_steps, before triggering a cleanup.
   
   Current value of 10 is arbitrarily chosen, just so the #UndoRedoManager::limitTotalStepsInMemory
   isnt constantly being called, but also so we dont go to much over #ns_max_steps.
   */
  const int ns_nsteps_histerious = 10;
  
  
  struct ClearStateOnDestruct
  {
  protected:
    UndoRedoManager::State &m_state;
    
  public:
    ClearStateOnDestruct( UndoRedoManager::State &state, const UndoRedoManager::State new_state )
      : m_state( state )
    {
      m_state = new_state;
    }
    
    ~ClearStateOnDestruct()
    {
      m_state = UndoRedoManager::State::Neither;
    }
  };//struct DoWorkOnDestruct
}//namespace


bool UndoRedoManager::spec_key_equal( const spec_key_t &lhs, const spec_key_t &rhs )
{
  const std::weak_ptr<SpecMeas> &lptr = std::get<0>(lhs);
  const std::weak_ptr<SpecMeas> &rptr = std::get<0>(rhs);
  std::owner_less<std::weak_ptr<SpecMeas>> comp;
  const bool llr = comp(lptr,rptr);
  const bool rlr = comp(rptr,lptr);
  
  return ((llr == rlr) && (std::get<1>(lhs) == std::get<1>(rhs)));
}//bool spec_key_equal( const spec_key_t &lhs, const spec_key_t &rhs )



UndoRedoManager::UndoRedoManager( InterSpec *parent )
 : Wt::WObject( parent ),
  m_state( UndoRedoManager::State::Neither ),
  m_steps{},
  m_step_offset{ 0 },
  m_current_specs{ nullptr },
  m_current_samples{},
  m_current_detectors{},
  m_num_steps_in_mem( 0 ),
  m_prev{},
  m_interspec( parent ),
  m_undoMenuDisableUpdate( this ),
  m_redoMenuDisableUpdate( this ),
  m_undoMenuToolTipUpdate( this ),
  m_redoMenuToolTipUpdate( this ),
  m_PeakModelChange_counter( 0 ),
  m_PeakModelChange_starting_peaks{},
  m_BlockUndoRedoInserts_counter( 0 ),
  m_BlockGuiUndoRedo_counter( 0 )
{
  assert( m_interspec );
  if( !m_interspec )
    throw logic_error( "UndoRedoManager must be initializes with valid InterSpec pointer." );
  
  const SpecUtils::SpectrumType types[] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  for( const SpecUtils::SpectrumType type : types )
  {
    const int index = static_cast<int>( type );
    m_current_specs[index] = m_interspec->measurment( type );
    m_current_samples[index] = m_interspec->displayedSamples( type );
    m_current_detectors[index] = m_interspec->detectorsToDisplay( type );
  }
    
  m_interspec->displayedSpectrumChanged().connect(
                boost::bind( &UndoRedoManager::handleSpectrumChange, this,
                             boost::placeholders::_1, boost::placeholders::_2,
                            boost::placeholders::_3, boost::placeholders::_4 ) );
}//UndoRedoManager


UndoRedoManager::~UndoRedoManager()
{
  
}//~UndoRedoManager()

void UndoRedoManager::PeakModelChange::setToCurrentPeaks()
{
  UndoRedoManager *manager = UndoRedoManager::instance();
  if( !manager || !manager->m_interspec )
    return;
  
  auto &start_peaks = manager->m_PeakModelChange_starting_peaks;
  assert( start_peaks.empty() );
  start_peaks.clear();
  
  PeakModel *pmodel = manager->m_interspec->peakModel();
  assert( pmodel );
  if( !pmodel )
    return;
  
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = pmodel->peaks();
  if( peaks )
    start_peaks.insert( end(start_peaks), begin(*peaks), end(*peaks) );
}//void UndoRedoManager::PeakModelChange::setToCurrentPeaks()


UndoRedoManager::PeakModelChange::PeakModelChange()
{
  UndoRedoManager *manager = UndoRedoManager::instance();
  if( !manager || !manager->m_interspec )
    return;
  
  if( manager->m_PeakModelChange_counter == 0 )
  {
    setToCurrentPeaks();
  }//if( manager->m_PeakModelChange_counter == 0 )
  
  manager->m_PeakModelChange_counter += 1;
}//PeakModelChange constructor


UndoRedoManager::PeakModelChange::~PeakModelChange()
{
  UndoRedoManager *manager = UndoRedoManager::instance();
  if( !manager || !manager->m_interspec )
    return;
  
  assert( manager->m_PeakModelChange_counter > 0 );
  if( manager->m_PeakModelChange_counter == 0 )
    return;
  
  manager->m_PeakModelChange_counter -= 1;
  if( manager->m_PeakModelChange_counter != 0 )
    return;
  
  const vector<shared_ptr<const PeakDef>> starting_peaks = manager->m_PeakModelChange_starting_peaks;
  manager->m_PeakModelChange_starting_peaks.clear();
  
  PeakModel *pmodel = manager->m_interspec->peakModel();
  assert( pmodel );
  if( !pmodel )
    return;
  
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks_now = pmodel->peaks();
  vector<shared_ptr<const PeakDef>> final_peaks;
  if( peaks_now )
    final_peaks.insert( end(final_peaks), begin(*peaks_now), end(*peaks_now) );
  
  if( starting_peaks == final_peaks )
    return;
  
  function<void(bool)> undo_redo = [starting_peaks,final_peaks]( const bool is_undo ){
    InterSpec *viewer = InterSpec::instance();
    assert( viewer );
    if( !viewer )
      return;
    
    PeakModel *pmodel = viewer->peakModel();
    assert( pmodel );
    if( !pmodel )
      return;
    
    pmodel->setPeaks( is_undo ? starting_peaks : final_peaks );
  };//undo_redo
  
  function<void()> undo = [undo_redo](){ undo_redo(true); };
  function<void()> redo = [undo_redo](){ undo_redo(false); };
  
  const int dpeaks = static_cast<int>(final_peaks.size()) - static_cast<int>(starting_peaks.size());
  string desc = (dpeaks == 0) ? "edit peak" : ((dpeaks < 0) ? "remove peak" : "add peak");
  if( abs(dpeaks) > 1 )
    desc += "s";
  
  manager->addUndoRedoStep( undo, redo, desc );
}//~PeakModelChange()


UndoRedoManager::BlockUndoRedoInserts::BlockUndoRedoInserts()
{
  UndoRedoManager *manager = UndoRedoManager::instance();
  if( manager )
    manager->m_BlockUndoRedoInserts_counter += 1;
}//BlockUndoRedoInserts constructor


UndoRedoManager::BlockUndoRedoInserts::~BlockUndoRedoInserts()
{
  UndoRedoManager *manager = UndoRedoManager::instance();
  if( manager )
  {
    assert( manager->m_BlockUndoRedoInserts_counter > 0 );
    if( manager->m_BlockUndoRedoInserts_counter > 0 )
      manager->m_BlockUndoRedoInserts_counter -= 1;
  }
}//~BlockUndoRedoInserts()


UndoRedoManager::BlockGuiUndoRedo::BlockGuiUndoRedo( Wt::WObject *parent )
  : Wt::WObject( parent ),
    m_valid( false ),
    m_peak_change( nullptr )
{
  UndoRedoManager *manager = UndoRedoManager::instance();
  if( !manager || !manager->canAddUndoRedoNow() )
    return;
  
  m_valid = true;
  m_peak_change = make_unique<UndoRedoManager::PeakModelChange>();
  
  if( manager->m_BlockGuiUndoRedo_counter == 0 )
  {
    manager->m_undoMenuDisableUpdate.emit( true );
    manager->m_redoMenuDisableUpdate.emit( true );
    manager->m_undoMenuToolTipUpdate.emit( Wt::WString() );
    manager->m_redoMenuToolTipUpdate.emit( Wt::WString() );
  }//if( manager->m_BlockGuiUndoRedo_counter == 0 )
  
  manager->m_BlockGuiUndoRedo_counter += 1;
}//BlockGuiUndoRedo constructor


UndoRedoManager::BlockGuiUndoRedo::~BlockGuiUndoRedo()
{
  if( !m_valid )
    return;
  
  UndoRedoManager *manager = UndoRedoManager::instance();
  if( !manager )
    return;
  
  assert( manager->m_BlockGuiUndoRedo_counter > 0 );
  if( manager->m_BlockGuiUndoRedo_counter > 0 )
    manager->m_BlockGuiUndoRedo_counter -= 1;
  
  if( manager->m_BlockGuiUndoRedo_counter == 0 )
    manager->updateMenuItemStates();
}//BlockGuiUndoRedo destructor


UndoRedoManager *UndoRedoManager::instance()
{
  InterSpec *viewer = InterSpec::instance();
  
  //assert( viewer );  //this will hit sometimes in the destruction of the app
  if( !viewer )
    return nullptr;
  
  UndoRedoManager *manager = viewer->undoRedoManager();
  return manager;
}//UndoRedoManager::instance()


int UndoRedoManager::maxUndoRedoSteps()
{
  return ns_max_steps.load();
}


void UndoRedoManager::setMaxUndoRedoSteps( const int steps )
{
  ns_max_steps = steps;
}


void UndoRedoManager::addUndoRedoStep( std::function<void()> undo,
                                       std::function<void()> redo,
                                       const std::string &description )
{
  assert( undo || redo );
  
  if( !m_steps )
  {
    Wt::log("debug") << "No spectrum file set for undo/redo.";
    return;
  }
  
  assert( m_step_offset <= m_steps->size() );
  
  if( m_state != State::Neither )
  {
    Wt::log("debug") << "UndoRedoManager::addUndoRedoStep(): not adding undo/redo step"
                        " - currently executing undo or redo.";
    return;
  }//if( m_state != State::Neither )
  
  if( m_BlockUndoRedoInserts_counter > 0 )
  {
    Wt::log("debug") << "UndoRedoManager::addUndoRedoStep(): not adding undo/redo step"
                        " - currently an active BlockUndoRedoInserts.";
    return;
  }
  
  if( m_BlockGuiUndoRedo_counter > 0 )
  {
    Wt::log("debug") << "UndoRedoManager::addUndoRedoStep(): not adding undo/redo step"
                        " - currently an active BlockGuiUndoRedo.";
    return;
  }
  
  const int max_steps = ns_max_steps;
  if( max_steps < 0 )
  {
    Wt::log("debug") << "UndoRedoManager::addUndoRedoStep(): undo/redo disabled app-wide.";
    return;
  }
  
  const size_t num_steps = m_steps->size();
  
  if( (m_step_offset != 0) && (m_step_offset <= num_steps) )
  {
    // If we are here, we have hit 'undo' one or more times, and now we are making a new edit
    //  that we want an undo step for.  Instead of just discarding `m_step_offset` undo steps,
    //  like some programs, we'll instead make it so if the user hits undo, after this new step
    //  they are adding, the undo's will unwind the previous undo steps, and then re-execute
    //  the undo steps.
    for( size_t index = 0; index < m_step_offset; ++index )
    {
      UndoRedoStep step = (*m_steps)[num_steps - 1 - index];
      std::swap( step.m_redo, step.m_undo );
      m_num_steps_in_mem += 1;
      m_steps->push_back( step );
    }//
  }//if( m_step_offset != 0 )
  
  m_step_offset = 0;
  m_num_steps_in_mem += 1;
  m_steps->push_back( {undo, redo, description, std::chrono::system_clock::now()} );
  
  m_undoMenuDisableUpdate.emit( false );
  m_undoMenuToolTipUpdate.emit( Wt::WString::fromUTF8(description) );
  m_redoMenuDisableUpdate.emit( true );
  m_redoMenuToolTipUpdate.emit( Wt::WString() );
  
  
  // Cleaning up the history isnt a super-cheap operation, so we'll We'll wait until we
  //  are #ns_nsteps_histerious over the max limit, to bother to clean things up.
  //  Also, we'll clean things up outside of the main event loop
  if( (max_steps != 0) && (m_num_steps_in_mem > (max_steps + ns_nsteps_histerious)) )
    Wt::WServer::instance()->schedule( 100, wApp->sessionId(),
                                  boost::bind( &UndoRedoManager::limitTotalStepsInMemory, this ) );
}//void addUndoRedoStep(...)


bool UndoRedoManager::canUndo() const
{
  return m_steps && (m_step_offset < m_steps->size());
}


bool UndoRedoManager::canRedo() const
{
  return (m_step_offset > 0);
}

void UndoRedoManager::executeUndo()
{
  if( m_BlockGuiUndoRedo_counter > 0 )
  {
    Wt::log("debug") << "Currently blocking GUI undo.";
    return;
  }
  
  if( !m_steps || (m_step_offset >= m_steps->size()) )
  {
    updateMenuItemStates();
    Wt::log("debug") << "No more undo steps to execute.";
    return;
  }
  
  ClearStateOnDestruct state_guard( m_state, UndoRedoManager::State::InUndo );
  
  UndoRedoStep *step = nullptr;
  while( m_step_offset < m_steps->size() )
  {
    UndoRedoStep *this_step = &((*m_steps)[m_steps->size() - 1 - m_step_offset]);
    m_step_offset += 1;
    
    if( this_step->m_undo )
    {
      step = this_step;
      break;
    }
  }//while( m_step_offset < m_steps->size() )
  
  assert( !step || step->m_undo );
  if( !step )
  {
    assert( m_step_offset == m_steps->size() );
    m_step_offset = m_steps->size();
    updateMenuItemStates();
    Wt::log("debug") << "No non-empty undo steps to execute.";
    return;
  }//if( !step )
  
  
  try
  {
    step->m_undo();
  }catch( std::exception &e )
  {
    const string message = "Error executing undo step: " + std::string(e.what());
    Wt::log("error") << message;
    passMessage( message, WarningWidget::WarningMsgHigh );
  }// try / catch
  
  updateMenuItemStates();
}//UndoRedoManager::executeUndo()


void UndoRedoManager::executeRedo()
{
  if( m_BlockGuiUndoRedo_counter > 0 )
  {
    Wt::log("debug") << "Currently blocking GUI redo.";
    return;
  }
  
  if( !m_steps || (m_step_offset == 0) || m_steps->empty() )
  {
    updateMenuItemStates();
    Wt::log("debug") << "No redo steps to execute.";
    return;
  }
  
  ClearStateOnDestruct state_guard( m_state, UndoRedoManager::State::InRedo );
  
  assert( m_step_offset <= m_steps->size() );
  if( m_step_offset > m_steps->size() )
    m_step_offset = m_steps->size();
  
  UndoRedoStep *step = nullptr;
  while( m_step_offset > 0 )
  {
    UndoRedoStep *this_step = &((*m_steps)[m_steps->size() - m_step_offset]);
    m_step_offset -= 1;
    
    if( this_step->m_redo )
    {
      step = this_step;
      break;
    }
  }//while( m_step_offset > 0 )
  
  assert( !step || step->m_redo );
  if( !step )
  {
    assert( m_step_offset == 0 );
    m_step_offset = 0;
    updateMenuItemStates();
    Wt::log("debug") << "No non-empty redo steps to execute.";
    return;
  }//if( !step )
  
  try
  {
    step->m_redo();
  }catch( std::exception &e )
  {
    const string message = "Error executing redo step: " + std::string(e.what());
    Wt::log("error") << message;
    passMessage( message, WarningWidget::WarningMsgHigh );
  }// try / catch
  
  updateMenuItemStates();
}//void UndoRedoManager::executeRedo()


bool UndoRedoManager::isInUndo() const
{
  return (m_state == State::InUndo);
}


bool UndoRedoManager::isInRedu() const
{
  return (m_state == State::InRedo);
}


bool UndoRedoManager::isInUndoOrRedo() const
{
  return (m_state != State::Neither);
}


bool UndoRedoManager::canAddUndoRedoNow() const
{
  if( !m_steps )
    return false;
  
  if( m_state != State::Neither )
    return false;
  
  if( m_BlockUndoRedoInserts_counter > 0 )
    return false;
  
  if( m_BlockGuiUndoRedo_counter > 0 )
    return false;
  
  return true;
}//bool canAddUndoRedoNow() const


Wt::Signal<bool> &UndoRedoManager::undoMenuDisableUpdate()
{
  return m_undoMenuDisableUpdate;
}


Wt::Signal<bool> &UndoRedoManager::redoMenuDisableUpdate()
{
  return m_redoMenuDisableUpdate;
}


Wt::Signal<Wt::WString> &UndoRedoManager::undoMenuToolTipUpdate()
{
  return m_undoMenuToolTipUpdate;
}


Wt::Signal<Wt::WString> &UndoRedoManager::redoMenuToolTipUpdate()
{
  return m_redoMenuToolTipUpdate;
}


void UndoRedoManager::clearUndoRedu()
{
  if( m_steps )
    m_steps->clear();
  
  updateMenuItemStates();
}//void clearUndoRedu()


void UndoRedoManager::updateMenuItemStates()
{
  if( !m_steps || m_steps->empty() )
  {
    m_undoMenuDisableUpdate.emit( true );
    m_redoMenuDisableUpdate.emit( true );
    m_undoMenuToolTipUpdate.emit( Wt::WString() );
    m_redoMenuToolTipUpdate.emit( Wt::WString() );
    
    return;
  }//if( !m_steps )
  
  
  // Find next undo step
  UndoRedoStep *undo_step = nullptr;
  for( size_t offset = m_step_offset; (offset < m_steps->size()) && !undo_step; ++offset )
  {
    UndoRedoStep *this_step = &((*m_steps)[m_steps->size() - 1 - m_step_offset]);
    
    if( this_step->m_undo )
      undo_step = this_step;
  }//while( m_step_offset < m_steps->size() )
  
  // Find next redo step
  UndoRedoStep *redo_step = nullptr;
  for( size_t offset = min(m_step_offset, m_steps->size()); (offset > 0) && !redo_step; --offset )
  {
    UndoRedoStep *this_step = &((*m_steps)[m_steps->size() - offset]);
    
    if( this_step->m_redo )
      redo_step = this_step;
  }//while( m_step_offset > 0 )
  
  
  m_undoMenuDisableUpdate.emit( !undo_step );
  auto undo_tooltip = undo_step ? Wt::WString::fromUTF8(undo_step->m_description) : Wt::WString();
  m_undoMenuToolTipUpdate.emit( undo_tooltip );
  
  m_redoMenuDisableUpdate.emit( !redo_step );
  auto redo_tooltip = redo_step ? Wt::WString::fromUTF8(redo_step->m_description) : Wt::WString();
  m_redoMenuToolTipUpdate.emit( redo_tooltip );
}//void updateMenuItemStates()


void UndoRedoManager::handleSpectrumChange( const SpecUtils::SpectrumType type,
                                           const shared_ptr<SpecMeas> &meas,
                                           const set<int> &sample_nums,
                                           const vector<string> &detector_names )
{
  const int type_index = static_cast<int>( type );
  if( (meas == m_current_specs[type_index]) && (sample_nums == m_current_samples[type_index]) )
    return;
 
  const shared_ptr<SpecMeas> prev_meas = std::move( m_current_specs[type_index] );
  const set<int> prev_samples = std::move( m_current_samples[type_index] );
  const vector<string> prev_dets = std::move( m_current_detectors[type_index] );
  const weak_ptr<SpecMeas> prev_weak = prev_meas;
  
  
  m_current_specs[type_index] = meas;
  m_current_samples[type_index] = sample_nums;
  m_current_detectors[type_index] = detector_names;
  
  const weak_ptr<SpecMeas> current_weak = meas;
  
  const bool was_valid = !!prev_meas;
  const bool is_valid = !!meas;
  
  
  if( type == SpecUtils::SpectrumType::Foreground )
  {
    // TODO: go through and cleanup files we are holding undo/redo for
    
    
    // We dont want to set the old peaks to the new spectrum, so we'll update it
    //  to the current spectrums peaks... not ideal, but something.
    if( m_PeakModelChange_counter > 0 )
    {
      UndoRedoManager::PeakModelChange::setToCurrentPeaks();
    }
    
    if( m_steps && !m_steps->empty() )
    {
      const spec_key_t key{ prev_weak, prev_samples };
      m_prev.push_front( {key,m_steps} );
    }//if( !m_steps.empty() )
    
    m_step_offset = 0;
    m_steps.reset();
    
    if( meas && sample_nums.size() )
    {
      const spec_key_t key{ current_weak, sample_nums };
      
      for( auto iter = begin(m_prev); iter != end(m_prev); ++iter )
      {
        const auto &prev_key = std::get<0>(*iter);
        if( spec_key_equal(key, prev_key) )
        {
          m_steps = std::get<1>( *iter );
          m_prev.erase( iter );
          break;
        }//if( we found a previous instance of this file/sample )
      }//for( loop over prev spec-file/sample-nums to see if we have undo/redo history )
      
      if( !m_steps )
        m_steps = make_shared<deque<UndoRedoStep>>();
    }//if( meas && sample_nums.size() )
    
    // TODO: we could insert an undo/redo step here to change back to the previous spectrum
  }else
  {
    auto undo = [prev_weak,prev_samples,prev_dets,type,was_valid](){
      InterSpec *viewer = InterSpec::instance();
      const shared_ptr<SpecMeas> prev_meas = prev_weak.lock();
      
      if( !viewer || (was_valid && !prev_meas) )
      {
        Wt::log("error") << "Previous back/second meas no longer in memory for undo.";
        return;
      }
      
      Wt::WFlags<InterSpec::SetSpectrumOptions> options;
      options |= InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal;
      options |= InterSpec::SetSpectrumOptions::CheckForRiidResults;
      options |= InterSpec::SetSpectrumOptions::SkipParseWarnings;
      
      viewer->setSpectrum( prev_meas, prev_samples, type, options );
    };//undo lamda
    
    auto redo = [current_weak,sample_nums,detector_names,type,is_valid](){
      InterSpec *viewer = InterSpec::instance();
      const shared_ptr<SpecMeas> current = current_weak.lock();
      
      if( !viewer || (is_valid && !current) )
      {
        Wt::log("error") << "Previous back/second meas no longer in memory for redo.";
        return;
      }
      
      Wt::WFlags<InterSpec::SetSpectrumOptions> options;
      options |= InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal;
      options |= InterSpec::SetSpectrumOptions::CheckForRiidResults;
      options |= InterSpec::SetSpectrumOptions::SkipParseWarnings;
      
      viewer->setSpectrum( current, sample_nums, type, options );
    };//redo lamda
    
    const string spec_type = (type == SpecUtils::SpectrumType::Background)
                              ? "background" : "secondary";
    
    addUndoRedoStep( undo, redo, "Change " + spec_type + " spectrum"  );
  }//if( type != SpecUtils::SpectrumType::Foreground )
  
  updateMenuItemStates();
}//void handleSpectrumChange( SpecUtils::SpectrumType type )


void UndoRedoManager::limitTotalStepsInMemory()
{
  // We shouldnt be in undo or redo state, but just in case, well check,
  //  and if so, try again later.
  assert( m_state == State::Neither );
  if( m_state != State::Neither )
  {
    Wt::WServer::instance()->schedule( 1000, wApp->sessionId(),
                                  boost::bind( &UndoRedoManager::limitTotalStepsInMemory, this ) );
    return;
  }
    
  const int max_steps = ns_max_steps;
  
  if( max_steps == 0 )
    return;
  
  if( max_steps < 0 )
  {
    if( m_steps || !m_prev.empty() )
    {
      m_steps = nullptr;
      m_prev.clear();
      m_step_offset = 0;
      m_num_steps_in_mem = 0;
      m_undoMenuDisableUpdate.emit(true);
      m_redoMenuDisableUpdate.emit(true);
      m_undoMenuToolTipUpdate.emit(Wt::WString());
      m_redoMenuToolTipUpdate.emit(Wt::WString());
      wApp->triggerUpdate();
    }//if( m_steps || !m_prev.empty() )
    
    return;
  }//if( max_steps < 0 )
  
  const size_t umax_steps = static_cast<size_t>( max_steps );
  
  m_num_steps_in_mem = 0;
  if( m_steps )
  {
    if( m_steps->size() > umax_steps )
    {
      // The most recent undo/redo steps are in the back of m_steps (e.g., we always
      //  do m_steps->push_back(...) )
      m_steps->erase( std::begin(*m_steps), std::end(*m_steps) - umax_steps );
      assert( m_steps->size() == umax_steps );
    }
    
    assert( m_steps->size() <= umax_steps );
    m_num_steps_in_mem = std::min(m_steps->size(), umax_steps);
  }//if( m_steps )
  
  assert( m_num_steps_in_mem <= umax_steps );
  
  for( auto iter = begin(m_prev); iter != end(m_prev); ++iter )
  {
    const auto &deque_ptr = std::get<1>( *iter );
    assert( deque_ptr );
    if( !deque_ptr )
      continue;
    
    if( (deque_ptr->size() + m_num_steps_in_mem) <= umax_steps )
    {
      // We can keep all items for this spec-file/sample-num
      m_num_steps_in_mem += deque_ptr->size();
      
      assert( m_num_steps_in_mem <= umax_steps );
    }else
    {
      // We have to limit the number of items for this spec-file/sample-num, and then
      //  remove all subsequent spec-file/sample-nums
      const size_t num_in_deque = deque_ptr->size();
      assert( m_num_steps_in_mem <= umax_steps );
      const size_t num_keep = umax_steps - m_num_steps_in_mem;
      
      if( num_keep == 0 )
      {
        m_prev.erase( iter, end(m_prev) );
      }else
      {
        assert( num_keep <= num_in_deque );
        const size_t num_erase = num_in_deque - num_keep;
        assert( num_erase <= num_in_deque );
        deque_ptr->erase( std::begin(*deque_ptr), std::begin(*deque_ptr) + num_erase );
        
        m_num_steps_in_mem += deque_ptr->size();
        assert( m_num_steps_in_mem == umax_steps );
        
        // Now erase all less-recent spec-file/sample-nums
        m_prev.erase( iter + 1, end(m_prev) );
      }//if( num_keep == 0 ) / else
      
      // We've deleted all the rest of m_prev, so we are done in this loop
      break;
    }//if( we can keep all these items ) / else
  }//for( auto iter = begin(m_prev); iter != end(m_prev); ++iter )
  
  
  m_step_offset = std::min(m_step_offset, (m_steps ? m_steps->size() : 0) );
}//void limitTotalStepsInMemory();
