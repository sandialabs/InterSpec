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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef RelActCalcAuto_PreFit_hpp
#define RelActCalcAuto_PreFit_hpp

#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcManual.h"

/** How RelActCalcAuto starts a fit when it models sum peaks (cascade summing or random pile-up): which peaks
 its initial estimate uses, and the order of the pre-fit stages.  Kept apart from the solver so the two
 decisions can be tested on their own.
 */
namespace RelActCalcAutoImp
{
namespace PreFit
{

/** Which of the peaks found in the spectrum to leave out of the initial estimate of the relative efficiency
 and activities, as they may be as much sum peak as line.

 Before anything is fit, each peak is taken to be the lines at its energy.  With the source close to the
 detector, or at a high count rate, a sum peak on top of a weak line can be several times the line (e.g.,
 U-235's 345.9 and 387.8 keV); an estimate made from such peaks can be far enough off that the fit does not
 recover.  So a peak is left out if the sums within its line-matching window could be more than
 `max_sum_share` of it - unless that would leave one of its sources with no peak at all, in which case that
 source keeps the left-out peak with its largest yield.

 @param peaks The peaks, with the lines (`m_source_gammas`) they were matched to.
 @param possible_sums For each modeled sum peak: its energy (keV), and the most counts it could have in the
        spectrum.
 @param cluster_num_sigma A line is matched to a peak within this many of the peak's widths (sigma).
 @param max_sum_share See above.
 @returns For each peak, whether to leave it out.
 */
inline std::vector<bool> peaks_mostly_sums( const std::vector<RelActCalcManual::GenericPeakInfo> &peaks,
                                            const std::vector<std::pair<double,double>> &possible_sums,
                                            const double cluster_num_sigma,
                                            const double max_sum_share )
{
  std::vector<bool> leave_out( peaks.size(), false );
  if( possible_sums.empty() )
    return leave_out;

  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const RelActCalcManual::GenericPeakInfo &peak = peaks[i];
    const double window = cluster_num_sigma * peak.m_fwhm / PhysicalUnits::fwhm_nsigma;
    double sum_counts = 0.0;
    for( const std::pair<double,double> &sum : possible_sums )
    {
      if( std::fabs( sum.first - peak.m_energy ) <= window )
        sum_counts += sum.second;
    }
    leave_out[i] = (sum_counts > (max_sum_share * peak.m_counts));
  }//for( loop over peaks )

  // A source left without a peak keeps the one its largest yield is in.
  std::map<std::string,bool> source_has_peak;
  std::map<std::string,std::pair<double,size_t>> source_best_left_out;  //{yield, peak index}
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    for( const RelActCalcManual::GenericLineInfo &line : peaks[i].m_source_gammas )
    {
      bool &has_peak = source_has_peak[line.m_isotope];
      has_peak = (has_peak || !leave_out[i]);
      if( leave_out[i] )
      {
        std::pair<double,size_t> &best
                    = source_best_left_out.emplace( line.m_isotope, std::make_pair( 0.0, i ) ).first->second;
        if( line.m_yield > best.first )
          best = { line.m_yield, i };
      }
    }//for( loop over the peak's lines )
  }//for( loop over peaks )

  for( const std::pair<const std::string,bool> &entry : source_has_peak )
  {
    if( !entry.second )
      leave_out[source_best_left_out.at( entry.first ).second] = false;
  }

  return leave_out;
}//peaks_mostly_sums(...)


/** The parameters of a fit, in the groups the pre-fit stages hold fixed or free (indices into the parameter
 vector; only parameters the fit varies - none that are held for the whole fit).
 */
struct StageGroups
{
  /** Every parameter the fit varies. */
  std::vector<size_t> free;

  /** All the energy-calibration parameters; of them, the offset and gain; and the deviation pairs. */
  std::vector<size_t> energy_cal, linear_cal, dev_pairs;

  std::vector<size_t> fwhm;

  /** The Physical model's correction function (e.g., Hoerl), self-attenuator, and external attenuators;
   and the other relative efficiency curves' coefficients. */
  std::vector<size_t> correction, self_atten, ext_atten, other_rel_eff;

  /** The peak skew, and with it the sum-peak parameters: these are held through every stage below except
   the one for the sums alone. */
  std::vector<size_t> skew;

  /** The sum-peak parameters - the cascade-summing solid angle and the random-sum scales (also in `skew`). */
  std::vector<size_t> sums;
};//struct StageGroups


/** The pre-fit stages: fits with some parameters held fixed, run in this order before the main fit to get
 near the solution.  Each entry is the parameters held in that stage.

 1. Everything held but the activities, ages and free peaks.
 2. With a correction function, everything freed but it and the skew.
 3. With the energy calibration fit: its offset and gain freed (deviation pairs still held); then its
    deviation pairs too.
 4. If no stage above freed the widths, they are fit once with the skew still held.
 5. The sum-peak parameters fit by themselves, with everything else where the stages above left it -
    before the deviation pairs are freed, if they are fit.  A ROI whose peaks are mostly sum peaks is, while
    the sums are still held at next to nothing, a ROI of unmodeled peaks; a deviation there runs to its
    limit to put some weak line under them, which moves the sum peaks off the data and leaves their scale
    at none.  So until this stage the deviation pairs are held too.
 */
inline std::vector<std::vector<size_t>> stages( const StageGroups &g )
{
  const auto join = []( const std::vector<std::vector<size_t>> &lists ) -> std::vector<size_t> {
    std::vector<size_t> answer;
    for( const std::vector<size_t> &list : lists )
    {
      for( const size_t index : list )
      {
        if( std::find( begin(answer), end(answer), index ) == end(answer) )
          answer.push_back( index );
      }
    }
    return answer;
  };//join

  const bool has_sums = !g.sums.empty();
  std::vector<size_t> all_but_sums;
  for( const size_t index : g.free )
  {
    if( std::find( begin(g.sums), end(g.sums), index ) == end(g.sums) )
      all_but_sums.push_back( index );
  }

  std::vector<std::vector<size_t>> answer;
  bool sums_done = !has_sums;
  const std::vector<size_t> dev_pairs_held_for_sums = has_sums ? g.dev_pairs : std::vector<size_t>{};

  answer.push_back( join( {g.energy_cal, g.fwhm, g.correction, g.self_atten, g.ext_atten, g.other_rel_eff, g.skew} ) );

  if( !g.correction.empty() )
    answer.push_back( join( {g.correction, g.skew, dev_pairs_held_for_sums} ) );

  if( !g.linear_cal.empty() )
  {
    answer.push_back( join( {g.dev_pairs, g.correction, g.skew} ) );

    if( !sums_done )
    {
      answer.push_back( all_but_sums );
      sums_done = true;
    }

    if( !g.dev_pairs.empty() )
      answer.push_back( join( {g.correction, g.skew} ) );
  }//if( the energy calibration is fit )

  if( !g.skew.empty() && !g.fwhm.empty() && g.correction.empty() && g.linear_cal.empty() )
    answer.push_back( g.skew );

  if( !sums_done )
    answer.push_back( all_but_sums );

  return answer;
}//stages(...)

}//namespace PreFit
}//namespace RelActCalcAutoImp

#endif //RelActCalcAuto_PreFit_hpp
