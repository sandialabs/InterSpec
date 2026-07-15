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
#include <memory>
#include <cassert>
#include <stdexcept>
#include <functional>

#include <Wt/WLogger>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/EnergyCalUndoRedo.h"

using namespace std;

using EnergyCalImp::meas_old_new_cal_t;
using EnergyCalImp::meas_old_new_peaks_t;

namespace
{
  typedef map<weak_ptr<SpecMeas>,meas_old_new_cal_t,owner_less<weak_ptr<SpecMeas>>> SpecMeasToCalHistoryMap;
  typedef map<weak_ptr<SpecMeas>,meas_old_new_peaks_t,owner_less<weak_ptr<SpecMeas>>> SpecMeasToPeakHistoryMap;

  thread_local static int sm_undo_redo_level;
  thread_local static unique_ptr<SpecMeasToCalHistoryMap> sm_meas_old_new_cal_map;
  thread_local static unique_ptr<SpecMeasToPeakHistoryMap> sm_meas_old_new_peaks_map;
  thread_local static unique_ptr<SpecMeasToPeakHistoryMap> sm_meas_old_new_hint_peaks_map;

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
    const set<int> &foreSamples = viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );

    const shared_ptr<SpecMeas> background = viewer->measurment( SpecUtils::SpectrumType::Background );
    const set<int> &backSamples = viewer->displayedSamples( SpecUtils::SpectrumType::Background );

    const shared_ptr<SpecMeas> secondary = viewer->measurment( SpecUtils::SpectrumType::SecondForeground );
    const set<int> &secoSamples = viewer->displayedSamples( SpecUtils::SpectrumType::SecondForeground );


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
      if( peakModel && (specfile == foreground) && (samples == foreSamples) )
        peakModel->setPeakFromSpecMeas( foreground, foreSamples, SpecUtils::SpectrumType::Foreground );
      else if( peakModel && (specfile == background) && (samples == backSamples) )
        peakModel->setPeakFromSpecMeas( background, backSamples, SpecUtils::SpectrumType::Background );
      else if( peakModel && (specfile == secondary) && (samples == secoSamples) )
        peakModel->setPeakFromSpecMeas( secondary, secoSamples, SpecUtils::SpectrumType::SecondForeground );
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
}//namespace


namespace EnergyCalImp
{

EnergyCalUndoRedoSentry::EnergyCalUndoRedoSentry()
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


meas_old_new_cal_t &EnergyCalUndoRedoSentry::cal_info( const shared_ptr<SpecMeas> &meas )
{
  assert( sm_meas_old_new_cal_map );
  if( !sm_meas_old_new_cal_map )
    throw logic_error( "EnergyCalUndoRedoSentry: cal info not initied?" );

  SpecMeasToCalHistoryMap &m = *sm_meas_old_new_cal_map;
  return m[meas];
}


meas_old_new_peaks_t &EnergyCalUndoRedoSentry::peak_info( const shared_ptr<SpecMeas> &meas )
{
  assert( sm_meas_old_new_peaks_map );
  if( !sm_meas_old_new_peaks_map )
    throw logic_error( "EnergyCalUndoRedoSentry: peak info not inited?" );

  SpecMeasToPeakHistoryMap &m = *sm_meas_old_new_peaks_map;
  return m[meas];
}


meas_old_new_peaks_t &EnergyCalUndoRedoSentry::hint_peak_info( const shared_ptr<SpecMeas> &meas )
{
  assert( sm_meas_old_new_hint_peaks_map );
  if( !sm_meas_old_new_hint_peaks_map )
    throw logic_error( "EnergyCalUndoRedoSentry: hint peak info not inited?" );

  SpecMeasToPeakHistoryMap &m = *sm_meas_old_new_hint_peaks_map;
  return m[meas];
}


EnergyCalUndoRedoSentry::~EnergyCalUndoRedoSentry()
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

  // No changes registered.  Note: the per-SpecMeas accessors create map entries on lookup, so a
  //  map entry whose change vectors are all empty (e.g., an operation that bailed out after it
  //  started recording) is also "no changes" - dont insert a no-op undo/redo step for those.
  bool any_changes = false;
  for( const auto &m : meas_old_new_cal_map )
    any_changes = (any_changes || !m.second.empty());
  for( const auto &m : meas_old_new_peaks_map )
    any_changes = (any_changes || !m.second.empty());
  for( const auto &m : meas_old_new_hint_peaks_map )
    any_changes = (any_changes || !m.second.empty());

  if( !any_changes )
    return;

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

      // Iterate the union of all the SpecMeas that had anything recorded; callers should always
      //  record calibration and peak changes together per SpecMeas, but a peaks-only record
      //  shouldnt get silently dropped if a future caller doesnt.
      map<weak_ptr<SpecMeas>,const meas_old_new_cal_t *,owner_less<weak_ptr<SpecMeas>>> per_meas_cals;
      for( const auto &key : meas_old_new_cal_map )
        per_meas_cals[key.first] = &key.second;
      for( const auto &key : meas_old_new_peaks_map )
      {
        if( !per_meas_cals.count(key.first) )
          per_meas_cals[key.first] = nullptr;
      }
      for( const auto &key : meas_old_new_hint_peaks_map )
      {
        if( !per_meas_cals.count(key.first) )
          per_meas_cals[key.first] = nullptr;
      }

      for( const auto &key : per_meas_cals )
      {
        const weak_ptr<SpecMeas> &specfile_weak = key.first;

        static const meas_old_new_cal_t no_cal_changes;
        const meas_old_new_cal_t &meas_old_new_cal = key.second ? *key.second : no_cal_changes;

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

}//namespace EnergyCalImp
