/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the GNU
 Free Software Foundation; either version 2.1 of the License, or (at your option)
 any later version.
 */

/* Local robustness sweep of RelActCalcAuto over a corpus of spectra x presets.

 This is a developer tool, not a unit test: the IDB/JRC corpora are ~1 GB of licensed scratch data
 that cannot live in the repository, so it is built but never registered with ctest.  Its job is to
 answer two questions that a single hand-picked spectrum cannot:

   1. Does every (spectrum, preset) pair either succeed or fail *gracefully*?  A preset may legitimately
      be unable to analyze a spectrum - the 610-775 keV window on low-statistics data, or a low-energy
      preset under heavy random/cascade summing - but it must never emit NaN/Inf, an out-of-range mass
      fraction, an unusable status with no explanation, a crash, or a hang.
   2. How do the profile-likelihood scans behave across a population?  Per-sample delta-chi2 is emitted
      so that non-monotonicity - a conditional solve which failed to minimize and therefore reports a
      value that is too *high* - can be measured over many spectra instead of anecdotally.

 Output is two TSVs: one row per (spectrum, preset) case, and one row per profile sample.  Shard with
 `--shard=i/n` to run several processes in parallel.
 */

#include "InterSpec_config.h"

#include <cmath>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/BatchRelActAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"

using namespace std;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

namespace
{
using Solution = RelActCalcAuto::RelActAutoSolution;

/** Every way a case can be non-graceful.  These are bugs regardless of whether the fit "worked";
    a preset failing to analyze a spectrum is expected and is NOT counted here. */
struct Anomalies
{
  vector<string> reasons;

  void add( const string &reason ){ reasons.push_back(reason); }

  string joined() const
  {
    string answer;
    for( size_t i = 0; i < reasons.size(); ++i )
      answer += (i ? "; " : "") + reasons[i];
    return answer.empty() ? string("-") : answer;
  }
};


bool finite_or_flag( const double value, const char *what, Anomalies &anomalies )
{
  if( std::isfinite(value) )
    return true;
  anomalies.add( string("non-finite ") + what );
  return false;
}


shared_ptr<const SpecUtils::Measurement> load_foreground( const string &path, string &error )
{
  try
  {
    auto file = make_shared<SpecUtils::SpecFile>();
    if( !file->load_file(path,SpecUtils::ParserType::Auto) )
    {
      error = "could not parse";
      return nullptr;
    }
    for( const shared_ptr<const SpecUtils::Measurement> &m : file->measurements() )
      if( m && (m->num_gamma_channels() > 64) )
        return m;
    error = "no gamma measurement";
  }catch( const std::exception &e )
  {
    error = string("parse threw: ") + e.what();
  }
  return nullptr;
}


/** Check everything a caller could read off a solution.  Deliberately independent of whether the
    status is usable: an unusable solution is still read by the GUI/report layer, so its fields must
    be finite and self-consistent too. */
void audit_solution( const Solution &solution, Anomalies &anomalies )
{
  const bool usable = Solution::is_usable_status(solution.m_status);
  switch( solution.m_status )
  {
    case Solution::Status::Success:
    case Solution::Status::UsableWithWarnings:
    case Solution::Status::NotInitiated:
    case Solution::Status::FailedToSetupProblem:
    case Solution::Status::FailToSolveProblem:
    case Solution::Status::UserCanceled:
      break;
    default:
      anomalies.add( "undefined status enum value" );
      break;
  }

  if( !usable && solution.m_error_message.empty() )
    anomalies.add( "unusable status with empty error message" );

  if( usable )
  {
    finite_or_flag( solution.m_chi2, "chi2", anomalies );
    finite_or_flag( solution.m_cov_scale, "cov_scale", anomalies );
    if( std::isfinite(solution.m_cov_scale) && (solution.m_cov_scale < 0.0) )
      anomalies.add( "negative cov_scale" );
  }

  for( size_t curve = 0; curve < solution.m_rel_activities.size(); ++curve )
  {
    for( const RelActCalcAuto::NuclideRelAct &act : solution.m_rel_activities[curve] )
    {
      finite_or_flag( act.rel_activity, "rel_activity", anomalies );
      finite_or_flag( act.rel_activity_uncertainty, "rel_activity_uncertainty", anomalies );
      if( std::isfinite(act.rel_activity) && (act.rel_activity < 0.0) )
        anomalies.add( "negative rel_activity" );
      if( std::isfinite(act.rel_activity_uncertainty) && (act.rel_activity_uncertainty < 0.0) )
        anomalies.add( "negative rel_activity_uncertainty" );
      if( act.age_was_fit )
      {
        finite_or_flag( act.age, "age", anomalies );
        finite_or_flag( act.age_uncertainty, "age_uncertainty", anomalies );
      }
    }
  }

  for( const vector<double> &coefficients : solution.m_rel_eff_coefficients )
    for( const double coefficient : coefficients )
      finite_or_flag( coefficient, "rel_eff_coefficient", anomalies );

  for( const vector<double> &row : solution.m_covariance )
    for( const double value : row )
      finite_or_flag( value, "covariance entry", anomalies );
}


/** Read every mass fraction through the structured accessor and check the physical invariants a
    renderer relies on.  Returns the per-curve sums so the caller can report normalization drift. */
void audit_mass_fractions( const Solution &solution, Anomalies &anomalies,
                           ostream &sample_out, const string &case_key )
{
  if( !Solution::is_usable_status(solution.m_status) )
    return;

  const SandiaDecay::Nuclide * const pu242 = DecayDataBaseServer::database()->nuclide("Pu242");

  for( size_t curve = 0; curve < solution.m_options.rel_eff_curves.size(); ++curve )
  {
    // The reported set is not the input set: the Pu-242 correlation deliberately generates a
    // reported isotope which is never a fitted NucInputInfo, and it participates in the
    // renormalization.  Summing only the inputs would under-count by the Pu-242 fraction.
    vector<const SandiaDecay::Nuclide *> reported;
    bool any_plutonium = false, has_pu242 = false;
    for( const RelActCalcAuto::NucInputInfo &input : solution.m_options.rel_eff_curves[curve].nuclides )
    {
      const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(input.source);
      if( !nuclide )
        continue;
      any_plutonium = any_plutonium || (nuclide->atomicNumber == 94);
      has_pu242 = has_pu242 || (nuclide == pu242);
      reported.push_back( nuclide );
    }
    if( any_plutonium && !has_pu242 && pu242 )
      reported.push_back( pu242 );

    double sum = 0.0;
    size_t counted = 0;
    for( const SandiaDecay::Nuclide * const nuclide : reported )
    {

      Solution::MassFractionResult result;
      try
      {
        result = solution.mass_enrichment_result( nuclide, curve );
      }catch( const std::exception &e )
      {
        anomalies.add( string("mass_enrichment_result threw for ") + nuclide->symbol
                       + ": " + e.what() );
        continue;
      }

      if( !std::isfinite(result.fraction) )
      {
        anomalies.add( "non-finite mass fraction for " + nuclide->symbol );
        continue;
      }
      if( (result.fraction < 0.0) || (result.fraction > 1.0) )
        anomalies.add( "mass fraction outside [0,1] for " + nuclide->symbol );
      sum += result.fraction;
      ++counted;

      if( result.covariance_one_sigma
          && (!std::isfinite(*result.covariance_one_sigma) || (*result.covariance_one_sigma < 0.0)) )
        anomalies.add( "invalid local sigma for " + nuclide->symbol );
      if( result.uncorrected_fraction && !std::isfinite(*result.uncorrected_fraction) )
        anomalies.add( "non-finite uncorrected fraction for " + nuclide->symbol );

      if( !result.profile )
        continue;

      for( const Solution::MassFractionProfileInterval &interval : result.profile->intervals )
      {
        if( !std::isfinite(interval.lower) || !std::isfinite(interval.upper) )
          anomalies.add( "non-finite profile interval for " + nuclide->symbol );
        else if( interval.lower > interval.upper )
          anomalies.add( "inverted profile interval for " + nuclide->symbol );
        else if( (interval.lower < 0.0) || (interval.upper > 1.0) )
          anomalies.add( "profile interval outside [0,1] for " + nuclide->symbol );
      }

      // Emit every profile sample.  A profile likelihood is a minimum over nuisance parameters, so a
      // sample whose delta-chi2 exceeds one further from the baseline is a failed minimization, not a
      // likelihood value.  Recording them lets that be measured over a population.
      for( const pair<double,double> &sample : result.profile->sampled_delta_chi2 )
        sample_out << case_key << '\t' << nuclide->symbol << '\t'
                   << static_cast<int>(result.profile->status) << '\t'
                   << setprecision(17) << sample.first << '\t'
                   << setprecision(17) << sample.second << '\n';
    }

    if( counted && (std::fabs(sum - 1.0) > 1.0e-6) )
    {
      ostringstream msg;
      msg << "mass fractions sum to " << setprecision(12) << sum << " on curve " << curve;
      anomalies.add( msg.str() );
    }
  }
}


string status_name( const Solution::Status status )
{
  switch( status )
  {
    case Solution::Status::Success:              return "Success";
    case Solution::Status::UsableWithWarnings:   return "UsableWithWarnings";
    case Solution::Status::NotInitiated:         return "NotInitiated";
    case Solution::Status::FailedToSetupProblem: return "FailedToSetupProblem";
    case Solution::Status::FailToSolveProblem:   return "FailToSolveProblem";
    case Solution::Status::UserCanceled:         return "UserCanceled";
  }
  return "Unknown";
}


string argument_value( const string &argument, const string &name )
{
  const string prefix = "--" + name + "=";
  return SpecUtils::istarts_with(argument,prefix) ? argument.substr(prefix.size()) : string();
}
}//anonymous namespace


int main( int argc, char **argv )
{
  string data_dir, out_path = "rel_act_sweep.tsv", sample_path, list_path, background_path;
  vector<string> preset_paths, spectrum_paths;
  bool auto_profile = false, resume = false;
  size_t shard_index = 0, shard_count = 1;

  for( int i = 1; i < argc; ++i )
  {
    const string argument = argv[i];
    string value;
    if( !(value = argument_value(argument,"datadir")).empty() )        data_dir = value;
    else if( !(value = argument_value(argument,"preset")).empty() )    preset_paths.push_back(value);
    else if( !(value = argument_value(argument,"spectrum")).empty() )  spectrum_paths.push_back(value);
    else if( !(value = argument_value(argument,"list")).empty() )      list_path = value;
    else if( !(value = argument_value(argument,"background")).empty() ) background_path = value;
    else if( !(value = argument_value(argument,"out")).empty() )       out_path = value;
    else if( !(value = argument_value(argument,"samples")).empty() )   sample_path = value;
    else if( argument == "--auto-profile" )                            auto_profile = true;
    else if( argument == "--resume" )                                  resume = true;
    else if( !(value = argument_value(argument,"shard")).empty() )
    {
      const size_t slash = value.find('/');
      if( slash == string::npos )
      {
        cerr << "Bad --shard (expected i/n): " << value << endl;
        return 1;
      }
      shard_index = static_cast<size_t>(std::stoul(value.substr(0,slash)));
      shard_count = static_cast<size_t>(std::stoul(value.substr(slash+1)));
    }
    else
    {
      cerr << "Unknown argument: " << argument << "\n"
           << "Usage: rel_act_corpus_sweep --datadir=DIR --preset=XML [--preset=XML ...]\n"
           << "         (--spectrum=FILE ... | --list=FILE) [--background=FILE]\n"
           << "         [--auto-profile] [--resume] [--out=TSV] [--samples=TSV] [--shard=i/n]" << endl;
      return 1;
    }
  }

  if( data_dir.empty() || preset_paths.empty() )
  {
    cerr << "--datadir and at least one --preset are required." << endl;
    return 1;
  }

  if( !list_path.empty() )
  {
    ifstream input( list_path.c_str() );
    if( !input )
    {
      cerr << "Could not open --list file: " << list_path << endl;
      return 1;
    }
    string line;
    while( std::getline(input,line) )
    {
      SpecUtils::trim(line);
      if( !line.empty() && (line[0] != '#') )
        spectrum_paths.push_back(line);
    }
  }

  if( spectrum_paths.empty() )
  {
    cerr << "No spectra supplied." << endl;
    return 1;
  }

  if( !shard_count || (shard_index >= shard_count) )
  {
    cerr << "Bad --shard values." << endl;
    return 1;
  }

  try
  {
    InterSpec::setStaticDataDirectory( data_dir );
  }catch( const std::exception &e )
  {
    cerr << "Could not set data directory: " << e.what() << endl;
    return 1;
  }
  if( !DecayDataBaseServer::database() )
  {
    cerr << "Decay database unavailable." << endl;
    return 1;
  }

  // Load every preset once; a malformed preset should stop the sweep rather than be silently skipped.
  vector<pair<string,RelActCalcAuto::Options>> presets;
  for( const string &path : preset_paths )
  {
    try
    {
      const shared_ptr<RelActCalcAuto::RelActAutoGuiState> state
          = BatchRelActAuto::load_state_from_xml_file( path );
      if( !state )
        throw runtime_error( "no state" );
      RelActCalcAuto::Options options = state->options;
      options.auto_profile_weak_mass_fractions = auto_profile;
      options.robust_solve = auto_profile;
      presets.emplace_back( SpecUtils::filename(path), std::move(options) );
    }catch( const std::exception &e )
    {
      cerr << "Could not load preset '" << path << "': " << e.what() << endl;
      return 1;
    }
  }

  shared_ptr<const SpecUtils::Measurement> background;
  if( !background_path.empty() )
  {
    string error;
    background = load_foreground( background_path, error );
    if( !background )
    {
      cerr << "Could not load --background '" << background_path << "': " << error << endl;
      return 1;
    }
  }

  // A sweep exists to find crashes, so it must survive one.  Recording which cases already have a
  // row lets a supervising shell restart a shard after a hard crash without redoing hours of work,
  // and without the crashed case silently disappearing from the results.
  const auto read_case_keys = []( const string &path, const bool skip_header ){
    std::set<string> keys;
    ifstream input( path.c_str() );
    string line;
    while( std::getline(input,line) )
    {
      if( skip_header && SpecUtils::istarts_with(line,"spectrum\t") )
        continue;
      const size_t tab = line.find('\t');
      if( tab == string::npos )
        continue;
      const size_t second_tab = line.find('\t',tab+1);
      keys.insert( (second_tab == string::npos) ? line : line.substr(0,second_tab) );
    }
    return keys;
  };

  // A case that was started but produced no row crashed the process outright.  Recording it as a
  // `Crashed` result is the whole point: a hard crash is the least graceful failure there is, and
  // without this a resumed shard would simply run into it again forever.
  const string attempted_path = out_path + ".attempted.tsv";
  std::set<string> completed_cases, crashed_cases;
  if( resume )
  {
    completed_cases = read_case_keys( out_path, true );
    for( const string &attempted : read_case_keys(attempted_path,false) )
      if( !completed_cases.count(attempted) )
        crashed_cases.insert( attempted );
    if( !completed_cases.empty() || !crashed_cases.empty() )
      cerr << "Resuming: " << completed_cases.size() << " cases already recorded, "
           << crashed_cases.size() << " previously crashed and will be recorded as such." << endl;
  }

  ofstream attempted_out( attempted_path.c_str(), resume ? ios::app : ios::trunc );
  if( !attempted_out )
  {
    cerr << "Could not open '" << attempted_path << "'" << endl;
    return 1;
  }

  ofstream out( out_path.c_str(), resume ? ios::app : ios::trunc );
  if( !out )
  {
    cerr << "Could not open --out '" << out_path << "'" << endl;
    return 1;
  }
  if( completed_cases.empty() )
    out << "spectrum\tpreset\tstatus\tseconds\tchi2\tcov_scale\tkappa\trank_deficient"
           "\tnum_warnings\tprofiles_complete\tprofiles_failed\tmax_profile_fits\tanomalies\n";

  ofstream sample_out;
  if( sample_path.empty() )
    sample_path = out_path + ".samples.tsv";
  sample_out.open( sample_path.c_str() );
  if( !sample_out )
  {
    cerr << "Could not open --samples '" << sample_path << "'" << endl;
    return 1;
  }
  sample_out << "case\tnuclide\tprofile_status\tfraction\tdelta_chi2\n";

  size_t case_index = 0, anomaly_cases = 0, unusable_cases = 0, thrown_cases = 0;
  for( size_t spectrum_index = 0; spectrum_index < spectrum_paths.size(); ++spectrum_index )
  {
    // Shard on the spectrum, not on the flat case index: with `case_index % shard_count` and a
    // preset count sharing a factor with the shard count, each shard would draw the same preset
    // every time (6 presets over 6 shards gives each shard exactly one preset).  That both
    // unbalances the run badly and lets one crashed shard silently lose an entire preset.
    if( (spectrum_index % shard_count) != shard_index )
      continue;

    const string &spectrum_path = spectrum_paths[spectrum_index];
    string load_error;
    const shared_ptr<const SpecUtils::Measurement> foreground
        = load_foreground( spectrum_path, load_error );
    const string spectrum_name = SpecUtils::filename( spectrum_path );

    for( const pair<string,RelActCalcAuto::Options> &preset : presets )
    {
      ++case_index;

      const string case_key = spectrum_name + "|" + preset.first;
      const string resume_key = spectrum_name + "\t" + preset.first;
      if( completed_cases.count(resume_key) )
        continue;
      if( crashed_cases.count(resume_key) )
      {
        out << resume_key << "\tCrashed\t\t\t\t\t\t0\t0\t0\t0\t"
            << "the solve crashed the process outright on a previous attempt\n" << std::flush;
        ++anomaly_cases;
        continue;
      }
      attempted_out << resume_key << '\n' << std::flush;
      if( !foreground )
      {
        out << spectrum_name << '\t' << preset.first << "\tLoadFailed\t0\t\t\t\t\t0\t0\t0\t0\t"
            << load_error << '\n' << std::flush;
        ++anomaly_cases;
        continue;
      }

      Anomalies anomalies;
      Solution solution;
      const auto started = std::chrono::steady_clock::now();
      bool threw = false;
      try
      {
        solution = RelActCalcAuto::solve( preset.second, foreground, background, nullptr, {},
                                          PeakFitUtils::CoarseResolutionType::High, nullptr );
      }catch( const std::exception &e )
      {
        // A thrown solve is never graceful: callers are documented to receive a status.
        threw = true;
        anomalies.add( string("solve threw: ") + e.what() );
      }
      const double seconds = std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - started ).count();

      size_t profiles_complete = 0, profiles_failed = 0, max_profile_fits = 0;
      if( !threw )
      {
        audit_solution( solution, anomalies );
        audit_mass_fractions( solution, anomalies, sample_out, case_key );

        for( const auto &curve_profiles : solution.m_mass_fraction_profiles )
        {
          for( const auto &entry : curve_profiles )
          {
            max_profile_fits = (std::max)(max_profile_fits,entry.second.num_fits);
            if( entry.second.status == Solution::MassFractionProfileStatus::Failed )
              ++profiles_failed;
            else
              ++profiles_complete;
          }
        }
      }

      const bool usable = !threw && Solution::is_usable_status(solution.m_status);
      unusable_cases += !usable;
      thrown_cases += threw;
      anomaly_cases += !anomalies.reasons.empty();

      out << spectrum_name << '\t' << preset.first << '\t'
          << (threw ? string("Threw") : status_name(solution.m_status)) << '\t'
          << setprecision(4) << seconds << '\t'
          << setprecision(12) << solution.m_chi2 << '\t'
          << solution.m_cov_scale << '\t'
          << solution.m_jacobian_condition_number << '\t'
          << solution.m_num_rank_deficient_dirs << '\t'
          << solution.m_warnings.size() << '\t'
          << profiles_complete << '\t' << profiles_failed << '\t' << max_profile_fits << '\t'
          << anomalies.joined() << '\n' << std::flush;
    }
  }

  cerr << "Sweep finished: " << case_index << " cases considered, "
       << unusable_cases << " unusable, " << thrown_cases << " threw, "
       << anomaly_cases << " with anomalies." << endl;
  return anomaly_cases ? 2 : 0;
}
