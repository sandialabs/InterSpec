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
 - Be able to insert, or maybe just modify, an undo/redo step from within a undo/redo step
 - Add undo/redo menu items to the app
 - Add in support a ton more places
 - Add ability to remove a step later on (e.g., the close a one-off dialog via "undo", we should remove this step once
    the dialog is closed otherwise).  We could add a helper function to this class to close one-off dialogs.
 */
class UndoRedoManager : public Wt::WObject
{
public:
  UndoRedoManager( InterSpec *parent );
  
  virtual ~UndoRedoManager();

  static UndoRedoManager *instance();
  
  /** Adds an undo/redo step.
   
   If the user has "undo" one or more times, then the undo/redo functions of the undid steps will swapped and new
   steps added to the stack, and the new undo/redo step is then placed on the stack after these.
   
   You can pass in nullptr function for undo or redo; if you do this, then this step will be skipped over for the relevant
   undo/redo step.  This is so you can use "undo" to close an one-off dialog, and not specify a redo, and then later on
   there wont be a seeming empty redo step later on (although, depending if use use WApplication::bind(...) to close a
   specific window, the undo slot will seem empty to the user).
   
   Be careful of:
   - If you capture any shared pointers, particularly to SpecMeas objects, keep in mind you could create dependency
      cycles that could cause the object to never destruct, even if all other shared pointers are gone.  Using std::weak_ptr<SpecMeas>
      can help with this.
   - Do not create new undo/redo steps, from within the undo/redo functions, even incidentally - you'll create a undo/redo cycle
      that will make previous history inaccessible.
   */
  void addUndoRedoStep( std::function<void()> undo,
                        std::function<void()> redo,
                        const std::string &description = "" );
  
  bool canUndo() const;
  bool canRedo() const;
  
  void executeUndo();
  void executeRedo();
  
  /**
   
   TODO: add constructor with a function call for before applying changes, and after applying changes
   */
  struct PeakModelChange
  {
    PeakModelChange();
    ~PeakModelChange();
    
    PeakModelChange( const PeakModelChange& ) = delete; // non construction-copyable
    PeakModelChange &operator=( const PeakModelChange & ) = delete; // non copyable
  };//struct PeakModelChange
  
protected:
  void handleSpectrumChange( const SpecUtils::SpectrumType type,
                            const std::shared_ptr<SpecMeas> &meas,
                            const std::set<int> &sample_nums,
                            const std::vector<std::string> &detector_names );
  
protected:
  
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
  
  std::shared_ptr<SpecMeas> m_current_spec;
  std::set<int> m_current_samples;
  
  /** Eventually we want to have `SpecMeas` itself track its history, but for the moment we'll just track it here. */
  typedef std::tuple<std::weak_ptr<SpecMeas>,std::set<int> > spec_key_t;
  struct SpecKeyLess
  {
    bool operator()( const spec_key_t &lhs, const spec_key_t &rhs ) const noexcept;
  };//struct SpecKeyLess
  
  std::map< spec_key_t, std::shared_ptr<std::deque<UndoRedoStep>>, SpecKeyLess > m_prev;
  
  InterSpec *m_interspec;
  
  
  size_t m_PeakModelChange_counter;
  std::vector<std::shared_ptr<const PeakDef>> m_PeakModelChange_starting_peaks;
  
  friend class PeakModelChange;
};//class UndoRedoManager

#endif //UndoRedoManager_h
