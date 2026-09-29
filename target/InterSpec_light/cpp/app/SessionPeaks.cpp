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

#include <cmath>
#include <limits>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"

#include "PeakCsv.h"
#include "Session.h"

using namespace std;
using json = nlohmann::json;

typedef vector<shared_ptr<const PeakDef>> PeakShrdVec;

namespace
{
  const char * const sm_default_peak_color = "rgb(0,51,255)";

  /** Copies peaks, giving each ROI its own new PeakContinuum, so modifying the copies' continua
   cannot affect the originals (PeakDef copies otherwise share their continuum).
   */
  vector<PeakDef> deep_copy_peaks( const Session::PeakDeque &peaks )
  {
    map<shared_ptr<const PeakContinuum>,shared_ptr<PeakContinuum>> new_conts;
    vector<PeakDef> answer;
    for( const shared_ptr<const PeakDef> &p : peaks )
    {
      PeakDef copy = *p;
      shared_ptr<PeakContinuum> &cont = new_conts[p->continuum()];
      if( !cont )
        cont = make_shared<PeakContinuum>( *p->continuum() );
      copy.setContinuum( cont );
      answer.push_back( copy );
    }
    return answer;
  }//deep_copy_peaks(...)


  void sort_peaks( Session::PeakDeque &peaks )
  {
    std::sort( begin(peaks), end(peaks), &PeakDef::lessThanByMeanShrdPtr );
  }

  void remove_peaks( Session::PeakDeque &peaks, const PeakShrdVec &to_remove )
  {
    peaks.erase( std::remove_if( begin(peaks), end(peaks), [&to_remove]( const shared_ptr<const PeakDef> &p ){
      return std::find( begin(to_remove), end(to_remove), p ) != end(to_remove);
    } ), end(peaks) );
  }

  PeakShrdVec peaks_sharing_roi( const Session::PeakDeque &peaks, const shared_ptr<const PeakDef> &peak )
  {
    PeakShrdVec answer;
    for( const auto &p : peaks )
    {
      if( p->continuum() == peak->continuum() )
        answer.push_back( p );
    }
    return answer;
  }

  /** The color for `peak` given its new source `src`, from the library source `from`: the color of
   the shown reference lines it is from (as for a peak just fit), else of shown lines of its nuclide,
   element, or reaction; else its current color if its nuclide is unchanged; else the color of the
   other peaks of the nuclide, if they all have the same one (as InterSpec's peak table); else the
   default.
   */
  Wt::WColor color_for_source( const shared_ptr<const PeakDef> &peak, const PeakDef::Source &src,
                               const RefLib::Source *from, const vector<RefLib::Shown> &shown,
                               const Session::PeakDeque &peaks )
  {
    for( const RefLib::Shown &sh : shown )
    {
      if( sh.source && (sh.source == from) && !sh.color.empty() )
        return Wt::WColor( sh.color );
    }

    for( const RefLib::Shown &sh : shown )
    {
      const bool has_nuclide = sh.source && std::any_of( begin(sh.source->lines), end(sh.source->lines),
                                                         [&src]( const RefLib::Line &line ){
        return RefLib::same_nuclide( line.source, src );
      } );
      if( has_nuclide && !sh.color.empty() )
        return Wt::WColor( sh.color );
    }

    if( RefLib::same_nuclide( peak->source(), src ) )
      return peak->lineColor();

    vector<Wt::WColor> colors;
    for( const shared_ptr<const PeakDef> &p : peaks )
    {
      if( (p != peak) && RefLib::same_nuclide( p->source(), src ) && !p->lineColor().isDefault()
         && (std::find( begin(colors), end(colors), p->lineColor() ) == end(colors)) )
        colors.push_back( p->lineColor() );
    }
    return (colors.size() == 1) ? colors.front() : Wt::WColor();
  }//color_for_source(...)
}//namespace


Session::PeakDeque &Session::foregroundPeaks()
{
  Slot &fore = m_slots[Foreground];
  if( !fore.file || !fore.summed )
    throw runtime_error( "No foreground spectrum." );
  return fore.file->peaks[fore.samples];
}


const Session::PeakDeque *Session::foregroundPeaksConst() const
{
  const Slot &fore = m_slots[Foreground];
  if( !fore.file )
    return nullptr;
  const auto pos = fore.file->peaks.find( fore.samples );
  return (pos == end(fore.file->peaks)) ? nullptr : &pos->second;
}


int Session::effectiveDetType() const
{
  const Slot &fore = m_slots[Foreground];
  if( !fore.file || !fore.summed )
    return static_cast<int>( PeakFitUtils::CoarseResolutionType::Unknown );
  return static_cast<int>( PeakFitUtils::effective_det_type( fore.file->prefs, fore.summed, fore.file->spec ) );
}


bool Session::isHighRes() const
{
  return (effectiveDetType() == static_cast<int>(PeakFitUtils::CoarseResolutionType::High));
}


void Session::resetDragCaches()
{
  m_drag_continuum.reset();
  m_drag_last_peaks.clear();
  m_drag_time = {};
  m_create_roi_peaks.clear();
}


string Session::peaksToChartJson( const PeakShrdVec &peaks ) const
{
  const Slot &fore = m_slots[Foreground];
  return PeakDef::peak_json( peaks, fore.summed, Wt::WColor(sm_default_peak_color), 255 );
}


json Session::peaksJson() const
{
  const PeakDeque *peaks = foregroundPeaksConst();
  if( !peaks || peaks->empty() || !m_slots[Foreground].summed )
    return json::array();
  return json::parse( peaksToChartJson( PeakShrdVec( begin(*peaks), end(*peaks) ) ) );
}


json Session::peakListJson() const
{
  json answer = json::array();
  const PeakDeque *peaks = foregroundPeaksConst();
  if( !peaks )
    return answer;

  for( const auto &p : *peaks )
  {
    json pj;
    pj["mean"] = p->mean();
    pj["meanUnc"] = p->meanUncert();
    pj["fwhm"] = p->gausPeak() ? p->fwhm() : 0.0;
    pj["area"] = p->peakArea();
    pj["areaUnc"] = p->peakAreaUncert();
    pj["lower"] = p->lowerX();
    pj["upper"] = p->upperX();
    pj["gaussian"] = p->gausPeak();
    pj["continuumType"] = static_cast<int>( p->continuum()->type() );
    pj["skewType"] = static_cast<int>( p->skewType() );
    pj["chi2dof"] = p->chi2dof();
    pj["source"] = RefLib::source_label( p->source() );
    pj["useForCal"] = p->useForEnergyCalibration();
    pj["wantUseForCal"] = p->m_useForEnergyCal;
    pj["color"] = p->lineColor().cssText();
    answer.push_back( pj );
  }//for( loop over peaks )

  return answer;
}//json peakListJson() const


shared_ptr<const PeakDef> Session::peakContaining( const double energy ) const
{
  // Port of InterSpec::nearestPeak: nearest mean, among peaks whose ROI contains the energy
  const PeakDeque *peaks = foregroundPeaksConst();
  shared_ptr<const PeakDef> answer;
  if( !peaks )
    return answer;

  double min_de = std::numeric_limits<double>::infinity();
  for( const auto &p : *peaks )
  {
    const double de = fabs( p->mean() - energy );
    if( (de < min_de) && (energy > p->lowerX()) && (energy < p->upperX()) )
    {
      min_de = de;
      answer = p;
    }
  }
  return answer;
}//peakContaining(...)


void Session::assignSourceToNewPeak( PeakDef &peak, const string &refLineParent, PeakDeque &peaks )
{
  if( peak.hasSourceGammaAssigned() || m_shown_ref.empty() )
    return;

  const PeakShrdVec previous( begin(peaks), end(peaks) );
  const bool hpge = isHighRes();

  const auto apply_swap = [&peaks]( unique_ptr<pair<shared_ptr<const PeakDef>,PeakDef::Source>> swap ){
    if( !swap || !swap->first )
      return;
    const auto pos = std::find( begin(peaks), end(peaks), swap->first );
    if( pos == end(peaks) )
      return;
    auto newpeak = make_shared<PeakDef>( **pos );
    newpeak->setSource( swap->second );
    *pos = newpeak;
  };

  // A line under the mouse takes precedence (the chart passes e.g., "Th232;S.E. of 2614.5 keV")
  if( !refLineParent.empty() )
  {
    string parent = refLineParent.substr( 0, refLineParent.find( ';' ) );
    SpecUtils::trim( parent );
    for( const RefLib::Shown &sh : m_shown_ref )
    {
      if( sh.source && SpecUtils::iequals_ascii( sh.source->parent, parent ) )
        apply_swap( RefLib::assign_source( peak, previous, { sh }, hpge ) );
    }
  }//if( !refLineParent.empty() )

  if( !peak.hasSourceGammaAssigned() )
    apply_swap( RefLib::assign_source( peak, previous, m_shown_ref, hpge ) );
}//void assignSourceToNewPeak(...)


json Session::setRefLibrary( const json &p )
{
  m_ref_lib.load( p.at( "lib" ) );
  m_shown_ref.clear();
  json answer;
  answer["numSources"] = m_ref_lib.size();
  return answer;
}


json Session::setShownRefLines( const json &p )
{
  vector<RefLib::Shown> shown;
  for( const json &s : p.at( "sources" ) )
  {
    RefLib::Shown sh;
    sh.source = m_ref_lib.find( s.at( "parent" ).get<string>() );
    sh.color = s.value( "color", string() );
    if( sh.source )
      shown.push_back( sh );
  }
  m_shown_ref = shown;
  return json::object();
}//json setShownRefLines( const json &p )


json Session::fitPeakAt( const json &p )
{
  // Port of PeakSearchGuiUtils::fit_peak_from_double_click and InterSpec::addPeak
  const double x = p.at( "energy" ).get<double>();
  const double pixPerKeV = p.value( "pixPerKeV", 2.0 );
  const string refLineParent = p.value( "refLineParent", string() );

  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();

  const PeakShrdVec orig_peaks( begin(peaks), end(peaks) );
  for( const auto &pk : orig_peaks )
  {
    if( !pk->gausPeak() && (x >= pk->lowerX()) && (x <= pk->upperX()) )
      return { {"message", "Can not fit a peak inside a data-defined ROI."} };
  }

  const pair<PeakShrdVec,PeakShrdVec> found
             = searchForPeakFromUser( x, pixPerKeV, fore.summed, orig_peaks, nullptr, nullptr, fore.file->prefs );

  if( found.first.empty() || (found.second.size() >= found.first.size()) )
  {
    char msg[128];
    snprintf( msg, sizeof(msg), "Couldn't find a peak near %.1f keV", x );
    return { {"message", msg} };
  }

  PeakShrdVec peakstoadd = found.first;
  PeakShrdVec existingpeaks = found.second;
  PeakShrdVec replacement_peaks;

  // Refit versions of existing peaks keep their source: match those with a source first...
  for( const auto &pk : found.first )
  {
    if( !pk->hasSourceGammaAssigned() )
      continue;

    int nearest = -1;
    double smallest = std::numeric_limits<double>::max();
    for( size_t i = 0; i < existingpeaks.size(); ++i )
    {
      const PeakDef::Source &a = existingpeaks[i]->source(), &b = pk->source();
      if( (a.kind != b.kind) || (a.name != b.name) )
        continue;
      const double dist = fabs( pk->mean() - existingpeaks[i]->mean() );
      if( dist < smallest )
      {
        nearest = static_cast<int>( i );
        smallest = dist;
      }
    }

    if( nearest >= 0 )
    {
      replacement_peaks.push_back( pk );
      existingpeaks.erase( begin(existingpeaks) + nearest );
      peakstoadd.erase( std::find( begin(peakstoadd), end(peakstoadd), pk ) );
    }
  }//for( loop over found peaks )

  // ...then the remaining existing peaks, by energy.
  for( const auto &prev : existingpeaks )
  {
    if( peakstoadd.empty() )
      break;
    size_t nearest = 0;
    double smallest = std::numeric_limits<double>::max();
    for( size_t i = 0; i < peakstoadd.size(); ++i )
    {
      const double dist = fabs( prev->mean() - peakstoadd[i]->mean() );
      if( dist < smallest )
      {
        nearest = i;
        smallest = dist;
      }
    }
    replacement_peaks.push_back( peakstoadd[nearest] );
    peakstoadd.erase( begin(peakstoadd) + nearest );
  }//for( loop over remaining existing peaks )

  remove_peaks( peaks, found.second );
  peaks.insert( end(peaks), begin(replacement_peaks), end(replacement_peaks) );
  sort_peaks( peaks );

  for( const auto &newpk : peakstoadd )
  {
    PeakDef peak = *newpk;
    assignSourceToNewPeak( peak, refLineParent, peaks );
    peaks.push_back( make_shared<const PeakDef>( peak ) );
    sort_peaks( peaks );
  }

  return state( PartPeaks );
}//json fitPeakAt( const json &p )


json Session::deletePeakAt( const json &p )
{
  const double energy = p.at( "energy" ).get<double>();
  const shared_ptr<const PeakDef> peak = peakContaining( energy );
  if( !peak )
    return { {"message", "No peak to delete there."} };

  remove_peaks( foregroundPeaks(), { peak } );
  resetDragCaches();
  return state( PartPeaks );
}//json deletePeakAt( const json &p )


json Session::clearPeaks( const json & )
{
  foregroundPeaks().clear();
  resetDragCaches();
  return state( PartPeaks );
}


json Session::erasePeaks( const json &p )
{
  // Port of InterSpec::excludePeaksFromRange
  double x0 = p.at( "e0" ).get<double>(), x1 = p.at( "e1" ).get<double>();
  if( x0 > x1 )
    std::swap( x0, x1 );

  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();
  const shared_ptr<const SpecUtils::Measurement> &data = fore.summed;

  vector<PeakDef> all_peaks = deep_copy_peaks( peaks );
  const vector<PeakDef> peaks_in_range = peaksTouchingRange( x0, x1, all_peaks );
  if( peaks_in_range.empty() )
    return json::object();

  for( const PeakDef &peak : peaks_in_range )
  {
    const auto iter = std::find( begin(all_peaks), end(all_peaks), peak );
    if( iter != end(all_peaks) )
      all_peaks.erase( iter );
  }

  const bool isHPGe = isHighRes();
  vector<PeakDef> peaks_to_keep;
  for( PeakDef peak : peaks_in_range )
  {
    if( (peak.mean() >= x0) && (peak.mean() <= x1) )
      continue;

    double lowx = 0.0, upperx = 0.0;
    findROIEnergyLimits( lowx, upperx, peak, data, isHPGe );

    const bool x0InPeak = ((x0 >= lowx) && (x0 <= upperx));
    const bool x1InPeak = ((x1 >= lowx) && (x1 <= upperx));
    if( (x0 >= upperx) || (x1 <= lowx) || (!x0InPeak && !x1InPeak) )
    {
      peaks_to_keep.push_back( peak );
      continue;
    }

    if( x0InPeak && x1InPeak )
    {
      if( x0 > peak.mean() )
        upperx = x0;
      else
        lowx = x1;
    }else if( x0InPeak )
    {
      upperx = x0;
    }else
    {
      lowx = x1;
    }

    peak.continuum()->setRange( lowx, upperx );
    peaks_to_keep.push_back( peak );
  }//for( loop over peaks in range )

  map<shared_ptr<const PeakContinuum>,PeakShrdVec> peaksinroi;
  for( const PeakDef &peak : peaks_to_keep )
    peaksinroi[peak.continuum()].push_back( make_shared<PeakDef>( peak ) );

  const auto det_type = static_cast<PeakFitUtils::CoarseResolutionType>( effectiveDetType() );
  for( const auto &roi : peaksinroi )
  {
    const PeakShrdVec newpeaks = refitPeaksThatShareROI( data, nullptr, roi.second, det_type,
                                                          PeakFitLM::PeakFitLMOptions::SmallRefinementOnly );
    const PeakShrdVec &use = (newpeaks.size() == roi.second.size()) ? newpeaks : roi.second;
    for( const auto &pk : use )
      all_peaks.push_back( *pk );
  }

  peaks.clear();
  for( const PeakDef &pk : all_peaks )
    peaks.push_back( make_shared<const PeakDef>( pk ) );
  sort_peaks( peaks );

  resetDragCaches();
  return state( PartPeaks );
}//json erasePeaks( const json &p )


json Session::roiDrag( const json &p )
{
  // Port of D3SpectrumDisplayDiv::performExistingRoiEdgeDragWork (foreground only)
  double new_lower = p.at( "newLower" ).get<double>();
  double new_upper = p.at( "newUpper" ).get<double>();
  const double new_roi_px = p.value( "px", 100.0 );
  const double orig_lower = p.at( "origLower" ).get<double>();
  const bool isfinal = p.value( "isFinal", false );

  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();

  const auto drag_reply = [this]( const PeakShrdVec &newpeaks ) -> json {
    return { {"roiDragPeaks", newpeaks.empty() ? json(nullptr)
                : json::parse( PeakDef::gaus_peaks_to_json( newpeaks, m_slots[Foreground].summed,
                                                            Wt::WColor(sm_default_peak_color), 255 ) )} };
  };

  double min_de = std::numeric_limits<double>::max();
  shared_ptr<const PeakContinuum> continuum;
  for( const auto &pk : peaks )
  {
    const double de = fabs( pk->continuum()->lowerEnergy() - orig_lower );
    if( de < min_de )
    {
      min_de = de;
      continuum = pk->continuum();
    }
  }

  if( !continuum || (min_de > 1.0) )
  {
    json answer = drag_reply( {} );
    answer["message"] = "Could not find the ROI being dragged.";
    return answer;
  }

  // Keep the C++ value of the edge not being dragged, to avoid rounding drift
  const bool dragging_upper = (fabs(new_lower - continuum->lowerEnergy()) < fabs(new_upper - continuum->upperEnergy()));
  if( dragging_upper )
    new_lower = continuum->lowerEnergy();
  else
    new_upper = continuum->upperEnergy();

  auto new_continuum = make_shared<PeakContinuum>( *continuum );
  new_continuum->setRange( new_lower, new_upper );

  bool all_gaus = true;
  PeakShrdVec new_roi_peaks, orig_roi_peaks;
  for( const auto &pk : peaks )
  {
    if( pk->continuum() != continuum )
      continue;
    orig_roi_peaks.push_back( pk );
    auto newpeak = make_shared<PeakDef>( *pk );
    newpeak->setContinuum( new_continuum );
    all_gaus = (all_gaus && newpeak->gausPeak());
    new_roi_peaks.push_back( newpeak );
  }

  // Drop peaks more than a sigma outside the new ROI, furthest first, keeping at least one
  if( all_gaus && (new_roi_peaks.size() > 1) )
  {
    map<double,shared_ptr<const PeakDef>> distance_to_roi;
    for( const auto &pk : new_roi_peaks )
    {
      const double m = pk->mean();
      if( ((m + pk->sigma()) < new_lower) || ((m - pk->sigma()) > new_upper) )
        distance_to_roi[-std::min(m - new_lower, m - new_upper)] = pk;
    }
    for( const auto &dp : distance_to_roi )
    {
      if( new_roi_peaks.size() > 1 )
        new_roi_peaks.erase( std::find( begin(new_roi_peaks), end(new_roi_peaks), dp.second ) );
    }
  }//if( multiple Gaussian peaks )

  if( !all_gaus )
  {
    if( !isfinal )
      return drag_reply( new_roi_peaks );

    remove_peaks( peaks, orig_roi_peaks );
    peaks.insert( end(peaks), begin(new_roi_peaks), end(new_roi_peaks) );
    sort_peaks( peaks );
    return state( PartPeaks );
  }//if( !all_gaus )

  if( new_roi_px <= 10.0 )
  {
    // User dragged the ROI (nearly) closed: delete its peaks
    if( !isfinal )
      return drag_reply( {} );

    remove_peaks( peaks, orig_roi_peaks );
    resetDragCaches();
    return state( PartPeaks );
  }//if( ROI dragged closed )

  // Start from the previous drag's fit, if recent, so parameters dont wander
  const auto now = std::chrono::steady_clock::now();
  if( (m_drag_continuum == continuum)
     && (new_roi_peaks.size() == m_drag_last_peaks.size())
     && !m_drag_last_peaks.empty()
     && ((now - m_drag_time) < std::chrono::seconds(15)) )
  {
    auto cont = make_shared<PeakContinuum>( *m_drag_last_peaks.front()->continuum() );
    cont->setRange( new_lower, new_upper );
    new_roi_peaks.clear();
    for( const auto &pk : m_drag_last_peaks )
    {
      auto newpeak = make_shared<PeakDef>( *pk );
      newpeak->setContinuum( cont );
      new_roi_peaks.push_back( newpeak );
    }
  }//if( reuse previous drag fit )

  Wt::WFlags<PeakFitLM::PeakFitLMOptions> refit_options( PeakFitLM::PeakFitLMOptions::MediumAmplitudeRefinementOnly );
  {
    vector<shared_ptr<PeakDef>> mutable_peaks;
    for( const auto &pk : new_roi_peaks )
      mutable_peaks.push_back( make_shared<PeakDef>( *pk ) );
    apply_fwhm_method_to_peaks( mutable_peaks, nullptr, *fore.file->prefs, refit_options );
    new_roi_peaks.assign( begin(mutable_peaks), end(mutable_peaks) );
    if( fore.file->prefs->m_fwhm_method == PeakFitDetPrefs::FwhmMethod::Normal )
      refit_options |= PeakFitLM::PeakFitLMOptions::MediumFwhmRefinementOnly;
  }

  const auto det_type = static_cast<PeakFitUtils::CoarseResolutionType>( effectiveDetType() );
  const PeakShrdVec refit = refitPeaksThatShareROI( fore.summed, nullptr, new_roi_peaks, det_type, refit_options );

  m_drag_continuum = continuum;
  m_drag_last_peaks = refit;
  m_drag_time = std::chrono::steady_clock::now();

  const PeakShrdVec &newpeaks = refit.empty() ? new_roi_peaks : refit;
  if( !isfinal )
    return drag_reply( newpeaks );

  remove_peaks( peaks, orig_roi_peaks );
  peaks.insert( end(peaks), begin(newpeaks), end(newpeaks) );
  sort_peaks( peaks );
  resetDragCaches();
  return state( PartPeaks );
}//json roiDrag( const json &p )


json Session::fitRoiDrag( const json &p )
{
  // Port of D3SpectrumDisplayDiv::performDragCreateRoiWork (synchronous, foreground only)
  double lower = p.at( "lower" ).get<double>(), upper = p.at( "upper" ).get<double>();
  const int nForcedPeaks = p.value( "npeaks", -1 );
  const bool isfinal = p.value( "isFinal", false );
  if( upper < lower )
    std::swap( lower, upper );

  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();
  const shared_ptr<const SpecUtils::Measurement> &data = fore.summed;
  const shared_ptr<const PeakFitDetPrefs> &prefs = fore.file->prefs;
  const auto det_type = static_cast<PeakFitUtils::CoarseResolutionType>( effectiveDetType() );

  vector<shared_ptr<PeakDef>> best;
  const bool use_cached = isfinal && !m_create_roi_peaks.empty()
                          && ((nForcedPeaks <= 0) || (static_cast<size_t>(nForcedPeaks) == m_create_roi_peaks.size()));
  if( use_cached )
  {
    for( const auto &pk : m_create_roi_peaks )
      best.push_back( make_shared<PeakDef>( *pk ) );
  }else
  {
    const float erange = static_cast<float>( upper - lower );
    const float midenergy = static_cast<float>( 0.5*(lower + upper) );

    float min_sigma, max_sigma;
    expected_peak_width_limits( midenergy, det_type, data, min_sigma, max_sigma );
    if( erange < min_sigma )
    {
      m_create_roi_peaks.clear();
      return { {"roiDragPeaks", nullptr} };
    }

    const size_t start_channel = data->find_gamma_channel( static_cast<float>(lower) );
    const size_t end_channel = data->find_gamma_channel( static_cast<float>(upper) );
    vector<tuple<float,float,float>> candidates;
    secondDerivativePeakCanidates( data, prefs, start_channel, end_channel, candidates );
    const int ncandidates = static_cast<int>( candidates.size() );

    int npeaks = std::max( ncandidates, 1 );
    npeaks = std::min( npeaks, static_cast<int>(2*erange/min_sigma) );
    npeaks = std::max( npeaks, 1 );

    vector<int> npeakstry;
    if( (nForcedPeaks > 0) && (nForcedPeaks < 10) )
    {
      npeakstry.push_back( nForcedPeaks );
    }else
    {
      if( npeaks > 1 )
        npeakstry.push_back( npeaks - 1 );
      npeakstry.push_back( npeaks );
      if( (erange/(npeaks + 1)) > min_sigma )
        npeakstry.push_back( npeaks + 1 );
    }

    int best_choice = -1;
    vector<double> chi2s( npeakstry.size(), std::numeric_limits<double>::quiet_NaN() );
    vector<vector<shared_ptr<PeakDef>>> results( npeakstry.size() );
    for( size_t i = 0; i < npeakstry.size(); ++i )
    {
      const auto method = (npeakstry[i] == ncandidates) ? MultiPeakInitialGuessMethod::FromDataInitialGuess
                                                        : MultiPeakInitialGuessMethod::UniformInitialGuess;
      findPeaksInUserRange( lower, upper, npeakstry[i], method, data, nullptr, prefs, results[i], chi2s[i] );
      if( results[i].empty() || !std::isfinite(chi2s[i]) )
        continue;

      // Require ~10% better chi2 for each extra peak
      const int npeaksdiff = npeakstry[i] - ((best_choice >= 0) ? npeakstry[best_choice] : 0);
      const double weight = std::max( 0.25, 1.0 - 0.10*npeaksdiff );
      if( (best_choice < 0) || (chi2s[i] < weight*chi2s[best_choice]) )
        best_choice = static_cast<int>( i );
    }//for( loop over number of peaks to try )

    if( best_choice < 0 )
    {
      m_create_roi_peaks.clear();
      json answer = { {"roiDragPeaks", nullptr} };
      if( isfinal )
        answer["message"] = "Failed to fit any peaks in the selected range.";
      return answer;
    }

    best = results[best_choice];

    // Apply the fit preferences (FWHM method and skew), and refit
    Wt::WFlags<PeakFitLM::PeakFitLMOptions> prefs_options;
    apply_fit_prefs_to_peaks( best, data, nullptr, *prefs, prefs_options );
    prefs_options |= PeakFitLM::SmallAmplitudeRefinementOnly;
    const PeakShrdVec refit = refitPeaksThatShareROI( data, nullptr, PeakShrdVec( begin(best), end(best) ),
                                                      det_type, prefs_options );
    if( !refit.empty() )
    {
      best.clear();
      for( const auto &pk : refit )
        best.push_back( make_shared<PeakDef>( *pk ) );
    }

    // Assign sources, treating the already-assigned new peaks as existing peaks
    PeakDeque scratch = peaks;
    for( auto &pk : best )
    {
      assignSourceToNewPeak( *pk, "", scratch );
      scratch.push_back( pk );
    }
  }//if( use_cached ) / else

  if( !isfinal )
  {
    m_create_roi_peaks.assign( begin(best), end(best) );
    return { {"roiDragPeaks", json::parse( PeakDef::gaus_peaks_to_json( PeakShrdVec( begin(best), end(best) ),
                                                data, Wt::WColor(sm_default_peak_color), 255 ) )} };
  }

  for( const auto &pk : best )
    peaks.push_back( pk );
  sort_peaks( peaks );
  resetDragCaches();
  return state( PartPeaks );
}//json fitRoiDrag( const json &p )


json Session::peakInfoAt( const json &p )
{
  const double energy = p.at( "energy" ).get<double>();
  const shared_ptr<const PeakDef> peak = peakContaining( energy );
  if( !peak )
    return { {"peak", nullptr} };

  const PeakDeque *peaks = foregroundPeaksConst();
  json pj;
  pj["mean"] = peak->mean();
  pj["gaussian"] = peak->gausPeak();
  pj["continuumType"] = static_cast<int>( peak->continuum()->type() );
  pj["skewType"] = static_cast<int>( peak->skewType() );
  pj["numPeaksInRoi"] = peaks ? peaks_sharing_roi( *peaks, peak ).size() : 1;
  pj["source"] = RefLib::source_label( peak->source() );

  // Sources to offer from the shown reference lines
  json candidates = json::array();
  for( const PeakDef::Source &src : RefLib::suggest_sources( *peak, m_shown_ref ) )
  {
    if( candidates.size() >= 5 )
      break;
    candidates.push_back( RefLib::source_label( src ) );
  }
  pj["candidates"] = candidates;

  return { {"peak", pj} };
}//json peakInfoAt( const json &p )


json Session::setContinuumType( const json &p )
{
  // Port of PeakSearchGuiUtils::change_continuum_type_from_right_click
  const double energy = p.at( "energy" ).get<double>();
  const int type_int = p.at( "type" ).get<int>();
  if( (type_int < 0) || (type_int > static_cast<int>(PeakContinuum::External)) )
    throw runtime_error( "Invalid continuum type" );
  const auto type = static_cast<PeakContinuum::OffsetType>( type_int );

  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();
  const shared_ptr<const PeakDef> peak = peakContaining( energy );
  if( !peak )
    return { {"message", "No ROI there."} };

  const PeakShrdVec old_roi_peaks = peaks_sharing_roi( peaks, peak );
  if( peak->continuum()->type() == type )
    return json::object();

  auto new_continuum = make_shared<PeakContinuum>( *peak->continuum() );
  new_continuum->setType( type );
  if( type == PeakContinuum::External )
    new_continuum->setExternalContinuum( estimateContinuum( fore.summed ) );

  PeakShrdVec new_peaks;
  for( const auto &pk : old_roi_peaks )
  {
    auto newpeak = make_shared<PeakDef>( *pk );
    newpeak->setContinuum( new_continuum );
    new_peaks.push_back( newpeak );
  }

  const auto det_type = static_cast<PeakFitUtils::CoarseResolutionType>( effectiveDetType() );

  if( !peak->gausPeak() )
  {
    // Data-defined peaks: just change the type
  }else if( new_peaks.size() > 1 )
  {
    const PeakShrdVec result = refitPeaksThatShareROI( fore.summed, nullptr, new_peaks, det_type, {} );
    if( result.size() != new_peaks.size() )
      return { {"message", "Changing the continuum type made the peaks insignificant; not changed."} };
    new_peaks = result;
  }else
  {
    vector<PeakDef> input{ *new_peaks.front() };
    const vector<PeakDef> output = fitPeaksInRange( peak->mean() - 0.1, peak->mean() + 0.1, 0.0, 0.0, 0.0,
                                                    input, fore.summed, {}, det_type );
    if( output.size() != 1 )
      return { {"message", "Changing the continuum type made the peak insignificant; not changed."} };
    new_peaks = { make_shared<const PeakDef>( output.front() ) };
  }//if( data defined ) / else if( multiple peaks ) / else

  remove_peaks( peaks, old_roi_peaks );
  peaks.insert( end(peaks), begin(new_peaks), end(new_peaks) );
  sort_peaks( peaks );
  resetDragCaches();
  return state( PartPeaks );
}//json setContinuumType( const json &p )


json Session::setSkewType( const json &p )
{
  // Port of PeakSearchGuiUtils::change_skew_type_from_right_click
  const double energy = p.at( "energy" ).get<double>();
  const int type_int = p.at( "type" ).get<int>();
  if( (type_int < 0) || (type_int >= static_cast<int>(PeakDef::NumSkewType)) )
    throw runtime_error( "Invalid skew type" );
  const auto type = static_cast<PeakDef::SkewType>( type_int );

  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();
  const shared_ptr<const PeakDef> peak = peakContaining( energy );
  if( !peak )
    return { {"message", "No ROI there."} };
  if( peak->skewType() == type )
    return json::object();

  const PeakShrdVec roi_peaks = peaks_sharing_roi( peaks, peak );
  for( const auto &pk : roi_peaks )
  {
    if( !pk->gausPeak() )
      return { {"message", "Can not change skew of a data-defined peak."} };
  }

  // Start from the skew values of the nearest peak with this skew type, if any
  shared_ptr<const PeakDef> near_skew_peak;
  double near_dist = std::numeric_limits<double>::max();
  for( const auto &pk : peaks )
  {
    const double dist = fabs( pk->mean() - peak->mean() );
    if( (pk->skewType() == type) && (dist < near_dist) )
    {
      near_dist = dist;
      near_skew_peak = pk;
    }
  }

  const size_t num_skew_pars = PeakDef::num_skew_parameters( type );
  vector<double> skew_pars( num_skew_pars, 0.0 );
  vector<bool> fit_for( num_skew_pars, true );
  for( size_t i = 0; i < num_skew_pars; ++i )
  {
    const auto ct = PeakDef::CoefficientType( PeakDef::SkewPar0 + i );
    double lower, upper, starting, step;
    if( !PeakDef::skew_parameter_range( type, ct, lower, upper, starting, step ) )
      throw logic_error( "inconsistent skew par def" );

    skew_pars[i] = starting;
    fit_for[i] = PeakDef::skew_parameter_fit_by_default( type, ct );
    if( near_skew_peak )
    {
      const double val = near_skew_peak->coefficient( ct );
      if( (val >= lower) && (val <= upper) )
      {
        skew_pars[i] = val;
        fit_for[i] = near_skew_peak->fitFor( ct );
      }
    }
  }//for( loop over skew parameters )

  PeakShrdVec new_peaks;
  for( const auto &pk : roi_peaks )
  {
    auto newpeak = make_shared<PeakDef>( *pk );
    newpeak->setSkewType( type );
    for( size_t i = 0; i < num_skew_pars; ++i )
    {
      const auto ct = PeakDef::CoefficientType( PeakDef::SkewPar0 + i );
      newpeak->set_coefficient( skew_pars[i], ct );
      newpeak->set_uncertainty( 0.0, ct );
      newpeak->setFitFor( ct, fit_for[i] );
    }
    new_peaks.push_back( newpeak );
  }

  const auto det_type = static_cast<PeakFitUtils::CoarseResolutionType>( effectiveDetType() );
  if( new_peaks.size() > 1 )
  {
    const PeakShrdVec result = refitPeaksThatShareROI( fore.summed, nullptr, new_peaks, det_type, {} );
    if( result.size() != new_peaks.size() )
      return { {"message", "Changing the skew type made the peaks insignificant; not changed."} };
    new_peaks = result;
  }else
  {
    vector<PeakDef> input{ *new_peaks.front() };
    const vector<PeakDef> output = fitPeaksInRange( peak->mean() - 0.1, peak->mean() + 0.1, 0.0, 0.0, 0.0,
                                                    input, fore.summed, {}, det_type );
    if( output.size() != 1 )
      return { {"message", "Changing the skew type made the peak insignificant; not changed."} };
    new_peaks = { make_shared<const PeakDef>( output.front() ) };
  }

  remove_peaks( peaks, roi_peaks );
  peaks.insert( end(peaks), begin(new_peaks), end(new_peaks) );
  sort_peaks( peaks );
  resetDragCaches();
  return state( PartPeaks );
}//json setSkewType( const json &p )


json Session::refitRoi( const json &p )
{
  // Port of PeakSearchGuiUtils::refit_peaks_from_right_click (standard refit)
  const double energy = p.at( "energy" ).get<double>();
  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();
  const shared_ptr<const PeakDef> peak = peakContaining( energy );
  if( !peak )
    return { {"message", "No peak to refit there."} };

  PeakShrdVec roi_peaks = peaks_sharing_roi( peaks, peak );
  std::sort( begin(roi_peaks), end(roi_peaks), &PeakDef::lessThanByMeanShrdPtr );
  const auto det_type = static_cast<PeakFitUtils::CoarseResolutionType>( effectiveDetType() );

  PeakShrdVec new_peaks;
  if( roi_peaks.size() > 1 )
    new_peaks = refitPeaksThatShareROI( fore.summed, nullptr, roi_peaks, det_type, {} );

  if( new_peaks.size() != roi_peaks.size() )
  {
    vector<PeakDef> input;
    for( const auto &pk : roi_peaks )
      input.push_back( *pk );
    const vector<PeakDef> output = fitPeaksInRange( input.front().mean() - 0.1, input.back().mean() + 0.1,
                                                    0.0, 0.0, 0.0, input, fore.summed, {}, det_type );
    if( output.size() != input.size() )
      return { {"message", "Failed to refit (peak became insignificant)."} };

    new_peaks.clear();
    for( const PeakDef &pk : output )
      new_peaks.push_back( make_shared<const PeakDef>( pk ) );
  }

  remove_peaks( peaks, roi_peaks );
  peaks.insert( end(peaks), begin(new_peaks), end(new_peaks) );
  sort_peaks( peaks );
  resetDragCaches();
  return state( PartPeaks );
}//json refitRoi( const json &p )


json Session::setPeakProperty( const json &p )
{
  // The peak with this mean (from the peak list), or else the peak at this energy (as the chart gives)
  PeakDeque &peaks = foregroundPeaks();
  PeakDeque::iterator pos = end(peaks);
  if( p.contains( "mean" ) )
  {
    const double mean = p["mean"].get<double>();
    pos = std::min_element( begin(peaks), end(peaks), [mean]( const auto &a, const auto &b ){
      return fabs(a->mean() - mean) < fabs(b->mean() - mean);
    } );
    if( (pos != end(peaks)) && (fabs((*pos)->mean() - mean) > 0.01) )
      pos = end(peaks);
  }else
  {
    pos = std::find( begin(peaks), end(peaks), peakContaining( p.at( "energy" ).get<double>() ) );
  }

  if( pos == end(peaks) )
    throw runtime_error( "Could not find peak" );

  auto newpeak = make_shared<PeakDef>( **pos );
  if( p.contains( "useForCal" ) )
    newpeak->useForEnergyCalibration( p["useForCal"].get<bool>() );
  if( p.value( "clearSource", false ) )
  {
    newpeak->clearSources();
    newpeak->setLineColor( Wt::WColor() );
  }

  // Source text the user typed, or picked from the suggestions
  string note;
  if( p.contains( "source" ) )
  {
    const RefLib::Source *from = nullptr;
    const optional<PeakDef::Source> src
                = RefLib::source_from_text( m_ref_lib, *newpeak, p["source"].get<string>(), from, note );
    if( !src )
    {
      newpeak->clearSources();
      newpeak->setLineColor( Wt::WColor() );
    }else if( !(*src == newpeak->source()) )
    {
      newpeak->setLineColor( color_for_source( *pos, *src, from, m_shown_ref, peaks ) );
      newpeak->setSource( *src );
    }
  }//if( p.contains( "source" ) )

  *pos = newpeak;

  json answer = state( PartPeaks );
  if( !note.empty() )
    answer["message"] = note;
  return answer;
}//json setPeakProperty( const json &p )


json Session::exportPeakCsv( const json &p )
{
  const Slot &fore = m_slots[Foreground];
  const PeakDeque *peaks = foregroundPeaksConst();
  if( !fore.file || !peaks || peaks->empty() )
    throw runtime_error( "There are no peaks to export." );

  const string path = p.value( "path", string("/tmp/interspec_light_peaks.csv") );
  std::ofstream out( path.c_str(), ios::out | ios::binary | ios::trunc );
  if( !out )
    throw runtime_error( "Could not open peak CSV file." );

  PeakCsv::write( out, *peaks, fore.summed );

  if( !out )
    throw runtime_error( "Failed writing peak CSV." );
  out.close();

  string base = fore.file->name;
  const size_t dot = base.find_last_of( '.' );
  if( (dot != string::npos) && (dot > 0) )
    base = base.substr( 0, dot );

  return { {"path", path}, {"filename", "peaks_" + base + ".CSV"} };
}//json exportPeakCsv( const json &p )


json Session::searchPeaks( const json & )
{
  // Port of PeakSearchGuiUtils::automated_search_for_peaks / search_for_peaks_worker: existing peaks
  //  are kept, and new peaks get sources from the shown reference lines (assign_srcs_from_ref_lines).
  const Slot &fore = m_slots[Foreground];
  PeakDeque &peaks = foregroundPeaks();
  const size_t nbefore = peaks.size();

  const auto existing = make_shared<const PeakDeque>( peaks );
  const PeakShrdVec found = ExperimentalAutomatedPeakSearch::search_for_peaks( fore.summed, existing, fore.file->prefs );

  // Assign sources to peaks without one, largest first, so smaller peaks cant take their lines
  PeakDeque answer;
  vector<shared_ptr<PeakDef>> unassigned;
  for( const auto &pk : found )
  {
    if( pk->hasSourceGammaAssigned() )
      answer.push_back( pk );
    else
      unassigned.push_back( make_shared<PeakDef>( *pk ) );
  }
  std::sort( begin(unassigned), end(unassigned), []( const shared_ptr<PeakDef> &a, const shared_ptr<PeakDef> &b ){
    return a->amplitude() > b->amplitude();
  } );

  for( const shared_ptr<PeakDef> &pk : unassigned )
  {
    if( !m_shown_ref.empty() )
    {
      auto swap = RefLib::assign_source( *pk, PeakShrdVec( begin(answer), end(answer) ), m_shown_ref, isHighRes() );
      const auto pos = swap ? std::find( begin(answer), end(answer), swap->first ) : end(answer);
      if( pos != end(answer) )
      {
        auto changed = make_shared<PeakDef>( **pos );
        changed->setSource( swap->second );
        *pos = changed;
      }
    }//if( there are reference lines to assign from )
    answer.push_back( pk );
  }//for( loop over peaks without a source )

  peaks = std::move( answer );
  sort_peaks( peaks );
  resetDragCaches();

  const size_t nfound = (peaks.size() > nbefore) ? (peaks.size() - nbefore) : 0;
  char msg[128];
  snprintf( msg, sizeof(msg), "Found %i new peak%s.", static_cast<int>(nfound), (nfound == 1) ? "" : "s" );
  json result = state( PartPeaks );
  result["message"] = msg;
  return result;
}//json searchPeaks( const json & )
