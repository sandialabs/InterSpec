#ifndef UndoRedoManager_h
#define UndoRedoManager_h
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
#include <tuple>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <functional>

#include <Wt/WObject>
#include <Wt/WSignal>
#include <Wt/WString>

class PeakDef;
class SpecMeas;
class InterSpec;

namespace SpecUtils
{
  enum class SpectrumType : int;
}

/**
 TODO items:
 - Have the SpecMeas hold the history for itself
 - Make it so this class collects all the undo/redo steps for a given event loop, and mark this object as needing update, so
   during render, all steps are collected up, and made into a single step - but this would require making this class inherit
   from Wt::WWidget (or more likely Wt::WCompositeWidget), instead of just Wt::WObject - or hooking into this class from
   `InterSpecApp::notify(...)` (but I am assuming all updating/rendering happens from within `notify`, which I'm 100% on)
 - `canAddUndoRedoNow()` and similar could also instead reference a variable that is not reset until the end of
   `InterSpecApp::notify(...)`, so objects that set undo/redo steps during `render(...)`, will properly be able
   to detect if they should add a step or not (particularly from within `addUndoRedoStep(...)`)
 */
class UndoRedoManager : public Wt::WObject
{
public:
  UndoRedoManager( InterSpec *parent );
  
  virtual ~UndoRedoManager();

  /** Get the UndoRedoManager for this wApp instance.
   
   Useful if you do not want to include InterSpec.h.
   */
  static UndoRedoManager *instance();
  
  /** The approximate maximum number of undo/redo steps that should be kept in memory.
   
   A value of zero indicates unlimited, a negative value indicates disabled.
   Default value is 250, but may be set be the `DesktopAppConfig` mechanism
   */
  static InterSpec_API int maxUndoRedoSteps();
  
  /** Sets the approximate maximum number of undo/redo steps to be kept in memory.
   
   Set to `0` for unlimited, or a negative value to disable.
   */
  static InterSpec_API void setMaxUndoRedoSteps( const int steps );
  
  /** Adds an undo/redo step.
   
   If the user has "undo" one or more times, then the undo/redo functions of the undid steps will swapped and new
   steps added to the stack, and the new undo/redo step is then placed on the stack after these.
   
   You can pass in nullptr function for undo or redo; if you do this, then this step will be skipped over for the relevant
   undo/redo step.  This is so you can use "undo" to close an one-off dialog, and not specify a redo, and then later on
   there wont be a seeming empty redo step later on (although, depending if use use WApplication::bind(...) to close a
   specific window, the undo slot will seem empty to the user).
   
   If you are currently executing an undo or redo step, the new step passed in will not be added.
   
   Be careful of:
   - If you capture any shared pointers, particularly to SpecMeas objects, keep in mind you could create dependency
      cycles that could cause the object to never destruct, even if all other shared pointers are gone.  Using std::weak_ptr<SpecMeas>
      can help with this.
   */
  void addUndoRedoStep( std::function<void()> undo,
                        std::function<void()> redo,
                        const std::string &description );
  
  bool canUndo() const;
  bool canRedo() const;
  
  void executeUndo();
  void executeRedo();
  
  bool isInUndo() const;
  bool isInRedu() const;
  bool isInUndoOrRedo() const;
  
  /** Returns if a call to #addUndoRedoStep will actually add an undo/redo step.
   i.e., if not currently executing a undo/redo, and there are no currently active #BlockUndoRedoInserts.
   */
  bool canAddUndoRedoNow() const;
  
  /** Signal to update the "undo" menu item.  Emitted argument is wether item should be disabled. */
  Wt::Signal<bool> &undoMenuDisableUpdate();
  
  /** Signal to update the "redo" menu item.  Emitted argument is wether item should be disabled. */
  Wt::Signal<bool> &redoMenuDisableUpdate();
  
  /** Updates the tool tip for the undo menu item, based on `description` passed into #addUndoRedoStep. */
  Wt::Signal<Wt::WString> &undoMenuToolTipUpdate();
  
  /** Updates the tool tip for the redo menu item, based on `description` passed into #addUndoRedoStep. */
  Wt::Signal<Wt::WString> &redoMenuToolTipUpdate();
  
  /** Clears all undo/redo history. */
  void clearUndoRedu();
  
  /** A struct that will insert a single peak-change undo/redo step, for all changes between the construction of
   the first `PeakModelChange` and the  last `PeakModelChange` destructed, for an InterSpec session.
   */
  struct PeakModelChange
  {
    PeakModelChange();
    ~PeakModelChange();
    static void setToCurrentPeaks();
    
    PeakModelChange( const PeakModelChange& ) = delete; // non construction-copyable
    PeakModelChange &operator=( const PeakModelChange & ) = delete; // non copyable
  };//struct PeakModelChange
  
  
  /** A struct that as long as it is in scope, all calls to #addUndoRedoStep will result in the undo/redo step NOT being
   added.
   
   \sa m_BlockUndoRedoInserts_counter
   */
  struct BlockUndoRedoInserts
  {
    BlockUndoRedoInserts();
    ~BlockUndoRedoInserts();
    
    BlockUndoRedoInserts( const BlockUndoRedoInserts& ) = delete; // non construction-copyable
    BlockUndoRedoInserts &operator=( const BlockUndoRedoInserts & ) = delete; // non copyable
  };//BlockUndoRedoInserts
  
  /** A struct that blocks the GUI from allowing undo/redo (disables menu items), as well as blocks undo/redo*/
  struct BlockGuiUndoRedo : public Wt::WObject
  {
    BlockGuiUndoRedo( Wt::WObject *parent );
    ~BlockGuiUndoRedo();
    
    BlockGuiUndoRedo( const BlockGuiUndoRedo & ) = delete; // non construction-copyable
    BlockGuiUndoRedo &operator=( const BlockGuiUndoRedo & ) = delete; // non copyable
    
  private:
    bool m_valid;
    
    /** We will save any peak changes between construction and destruction of this object - but that is currently the
     only items that will be saved.
     */
    std::unique_ptr<PeakModelChange> m_peak_change;
    
    // TODO: should also save energy calibration...
    //typedef std::map<std::weak_ptr<SpecMeas>,std::shared_ptr<const SpecUtils::EnergyCalibration>,std::owner_less<std::weak_ptr<SpecMeas>>> MeasToEnergyCal_t;
    //MeasToEnergyCal_t m_energy_cals[3];
  };//BlockGuiUndoRedo
  
  
  enum class State : int
  {
    Neither,
    InUndo,
    InRedo
  };//enum class State
  
protected:
  /** Emits `m_[undo|redo]MenuDisableUpdate` and `m_[undo|redo|MenuToolTipUpdate` signals
   to set the state that the menu items should be in.
   
   Currently we have things hooked up so values will be updated to the WWidget, even if they dont need to be.
   */
  void updateMenuItemStates();
  
  void handleSpectrumChange( const SpecUtils::SpectrumType type,
                            const std::shared_ptr<SpecMeas> &meas,
                            const std::set<int> &sample_nums,
                            const std::vector<std::string> &detector_names );
  
  /** Limits total number of steps held in #m_steps plus #m_prev to be less than that retunred by #maxUndoRedoSteps */
  void limitTotalStepsInMemory();
  
protected:
  
  State m_state;
  
  struct UndoRedoStep
  {
    std::function<void()> m_undo;
    std::function<void()> m_redo;
    std::string m_description;
    std::chrono::time_point<std::chrono::system_clock> m_time;
  };//struct UndoRedoStep
  
  /** Each new Undo/Redo step, we will `push` onto #m_steps.  But if we have executed an undo step,
   we will use #m_step_offset to track were we are at.
   */
  std::shared_ptr<std::deque<UndoRedoStep>> m_steps;
  
  /** Where we are at, relative to end of #m_steps.
   So if we havent done any undoes, this value will be zero. */
  size_t m_step_offset;
  
  /** Track current SpecMeas and sample numbers; indexed by SpecUtils::SpectrumType enum. */
  std::shared_ptr<SpecMeas> m_current_specs[3];
  std::set<int> m_current_samples[3];
  std::vector<std::string> m_current_detectors[3];
  
  /** The total number of undo/redo steps in memory, for all spectrum files (i.e., for all entries in #m_prev).*/
  size_t m_num_steps_in_mem;
  
  /** Eventually we want to have `SpecMeas` itself track its history, but for the moment we'll just track it here. */
  typedef std::tuple<std::weak_ptr<SpecMeas>,std::set<int> > spec_key_t;
  static bool spec_key_equal( const spec_key_t &lhs, const spec_key_t &rhs );
    
  /** The most recent spectrum-file/sample-numbers undo/redo items will be in the front.
   
   The current spectrum-file/sample-numbers will not be in this deque.
   */
  std::deque< std::tuple<spec_key_t, std::shared_ptr<std::deque<UndoRedoStep>>> > m_prev;
  
  InterSpec *m_interspec;
  
  Wt::Signal<bool> m_undoMenuDisableUpdate;
  Wt::Signal<bool> m_redoMenuDisableUpdate;
  Wt::Signal<Wt::WString> m_undoMenuToolTipUpdate;
  Wt::Signal<Wt::WString> m_redoMenuToolTipUpdate;
  
  size_t m_PeakModelChange_counter;
  std::vector<std::shared_ptr<const PeakDef>> m_PeakModelChange_starting_peaks;
  
  /** A counter, only ever changed by #BlockUndoRedoInserts, that if it is not zero, then no undo/redo steps
   will be added by #addUndoRedoStep.
   
   \sa BlockUndoRedoInserts
   */
  size_t m_BlockUndoRedoInserts_counter;
  
  /** A counter, only changed by #BlockGuiUndoRedo, that if is not zero, then no undo/redo steps will
   be added or executed; also menu items should be disabled.
   */
  size_t m_BlockGuiUndoRedo_counter;
  
  friend class PeakModelChange;
};//class UndoRedoManager

#endif //UndoRedoManager_h
