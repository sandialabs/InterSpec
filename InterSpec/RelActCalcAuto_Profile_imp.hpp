/*
  Post-fit mass-fraction profile likelihood implementation.

  This file is included by RelActCalcAuto.cpp after RelActAutoCostFcn is fully defined and while
  namespace RelActCalcAuto is open.  Keeping the conditional-solve orchestration here prevents the
  already-large solver translation unit from obscuring the public solve flow.
*/
#ifndef InterSpec_RelActCalcAuto_Profile_imp_hpp
#define InterSpec_RelActCalcAuto_Profile_imp_hpp

namespace
{
/** Add bounded profiles after the final ROI/model/basin solve.

 A conditional fixes one coordinate and reoptimizes every nuisance parameter under it.  Which
 coordinate is fixed is the central design trade of this file: the interval is quoted in the
 REPORTED coordinate (the mass fraction, after Pu-242 correlation and renormalization), so a scan
 that fixes anything else is biased narrow by the misalignment between the two - measured (2026-08
 coverage study) as 95% coverage tracking the squared pin/reported correlation in lockstep.

 Conditional points run one way: PINNED, in place, on the CARRIER chart.  The winning solve hands
 its converged `ceres::Problem` forward on the solution; `install_carrier_reparam` reparameterizes
 the target's element (through a temporary wide, non-binding [0,1] mass-fraction constraint, or a
 user constraint that already has that shape) so the target's own slot IS its pre-Pu242-correction
 mass fraction; that slot is removed from the tangent space with a `ceres::SubsetManifold`, and
 each point loads the previous point's parameter vector, writes the pinned coordinate, and
 re-solves.  The pinned statistic is therefore `min over {q_pre == q0}` - the exact profile of the
 reported coordinate, aligned by construction (for Pu the correlation renormalizes on top, a
 second-order effect measured negligible on the validation corpus).  No residual is
 added, so the Jacobian sparsity is untouched and a point costs about a second rather than minutes.

 There is deliberately NO inexact engine below this: a target the carrier route refuses (no slot
 scans the reported fraction exactly), or a solve that retained no problem to pin within, reports
 a structured `Failed` profile carrying the reason.

 Before any scan runs, a PROBE pass solves just the two opening probes of every would-be target
 (measured: ~93% of better-baseline discoveries happen there), so a baseline sitting in the wrong
 basin is re-seeded through a cheap matrix-free restart before scans spend their budget on it;
 `needs_matrix_polish` carries the obligation to run the one full candidate-matrix polish at the
 end of such a probe-hop chain.  Scans keep their own deferred-rebaseline flow as the backstop.
 */
/** Everything one profile target contributes to its own scan, so the scan itself can be written
 once for every `ProfileTarget::Kind`.

 Each lane's preamble does the kind-specific work - the mass-fraction lane builds its nominal
 through `mass_enrichment_result` and projects a control window through the sibling simplex; the
 activity/ratio/age lane picks an identity or ratio chart and caps an open box - and then fills
 this.  Everything after that (the evaluator, the probe pair, the scan, and all reporting) is
 shared, so a fix to the conditional-point handling cannot land in one kind and miss the others.

 Two pairs of bounds here are easy to confuse and must not be:
  - `parameter_box` is the RAW Ceres box.  Only a genuine edge of it is a feasibility statement,
    which is why `reached_feasible_bound` tests this and never the scan limits: an open activity
    or ratio box gets a synthetic cap in `scan_lower/scan_upper`, and reporting THAT as a bound is
    exactly what `EndpointKind::ScanRangeLimit` exists to prevent.
  - `domain_*` is the read-back REJECTION window (a value outside it means the conditional solve
    produced nonsense, so the sample is discarded), while `feasible_*` is the window whose edges
    `lower_bound_kind`/`upper_bound_kind` describe.  They coincide for activities and ages, but
    for a mass fraction the domain is `[0,1]` while the feasible window is the projected control
    window, which is usually much narrower.
 */
struct ScanSetup
{
  RelActCalcAuto::Options::ProfileTarget target;
  RelActAutoSolution::MassFractionProfileReason reason
      = RelActAutoSolution::MassFractionProfileReason::AutomaticWeak;

  size_t slot = 0;                       //!< the pinned parameter
  double chart_scale = 0.0;              //!< dq/dx of the pinned slot
  double baseline_parameter = 0.0;
  double nominal_reported = 0.0;         //!< the scan's nominal, and the level the refine holds
  std::pair<double,double> parameter_box{ 0.0, 0.0 };   //!< RAW Ceres box; see above
  double scan_lower = 0.0, scan_upper = 0.0;            //!< scan limits, already capped
  double parameter_one_sigma = 0.0;      //!< in the SCANNED coordinate

  double domain_lower = 0.0, domain_upper = 0.0;        //!< read-back rejection window
  double feasible_lower = 0.0, feasible_upper = 0.0;    //!< what the endpoint kinds describe
  /** Precomputed per lane: the mass-fraction lane scales it to the control window's span, the
   others to the reported magnitude.  Deriving one formula for both would give an infinite
   tolerance on an open window and flag every sample as at-a-bound. */
  double feasible_tolerance = 0.0;
  RelActAutoSolution::MassFractionProfileEndpointKind lower_bound_kind
      = RelActAutoSolution::MassFractionProfileEndpointKind::PhysicalLimit;
  RelActAutoSolution::MassFractionProfileEndpointKind upper_bound_kind
      = RelActAutoSolution::MassFractionProfileEndpointKind::PhysicalLimit;

  std::array<double,2> thresholds{{ 0.0, 0.0 }};

  /** Level-set refinement is MASS-FRACTION ONLY, and this must stay false for every other kind:
   `refine_on_level_set` takes its direction from the gradient of the reported mass fraction
   regardless of the target's kind, so refining another quantity would move along the wrong level
   set and quote the sample at a value the pin never held - a silently NARROW interval.  The lane
   folds `have_covariance` in here, so the evaluator does not recompute it. */
  bool use_level_set_refine = false;

  size_t curve_index = 0;                //!< for stats lines and the deferred discovery
  std::string semantic_key;              //!< full deterministic ranking key for a discovery

  // Per-lane user-visible text and stats keys, kept distinct so this consolidation is provably
  // behavior-preserving; unifying the wording is a separate, reviewable change.
  std::string stats_name;                //!< how stats lines name this target
  std::string dq_key = "dq_dx";          //!< the stats key for `chart_scale`
  std::string cancel_message;
  std::string domain_message;
  std::string failed_warning_suffix;     //!< e.g. " mass-fraction profile failed: "
  std::string plain_warning_suffix;      //!< e.g. " mass-fraction profile: "
  bool emit_interval_stats = false;
};//struct ScanSetup


/** What the caller must do after one target's scan. */
enum class TargetOutcome
{
  Continue,   //!< move to the next target
  AbortPass   //!< the user canceled; abandon the whole pass immediately
};


void add_mass_fraction_profiles( RelActAutoSolution &solution,
                                 std::shared_ptr<const SpecUtils::Measurement> foreground,
                                 std::shared_ptr<const SpecUtils::Measurement> background,
                                 std::shared_ptr<const DetectorPeakResponse> input_drf,
                                 const std::vector<std::shared_ptr<const PeakDef>> &all_peaks,
                                 const PeakFitUtils::CoarseResolutionType det_type,
                                 const std::shared_ptr<std::atomic_bool> &cancel_calc,
                                 const unsigned baseline_restart_count = 0,
                                 const bool needs_matrix_polish = false )
{
  using CovQuality = RelActAutoSolution::MassFractionCovarianceQuality;
  using EndpointKind = RelActAutoSolution::MassFractionProfileEndpointKind;
  using Interval = RelActAutoSolution::MassFractionProfileInterval;
  using Profile = RelActAutoSolution::MassFractionProfileResult;
  using ProfileReason = RelActAutoSolution::MassFractionProfileReason;
  using ProfileStatus = RelActAutoSolution::MassFractionProfileStatus;

  // A returned solution never carries a live optimizer problem, and the shared cost functor is
  // always left describing exactly the point the solution reports (invariant 7).  Both hold however
  // this function exits - including the early returns just below and any exception - so this guard
  // is the first thing constructed.
  const struct ProfileHostRelease
  {
    RelActAutoSolution &sol;
    std::shared_ptr<RelActCalcAutoImp::ProfileConditionalHost> host;
    ~ProfileHostRelease()
    {
      if( host )
        host->restore_optimum();
      sol.m_profile_host.reset();
    }
  } profile_host_release{ solution, solution.m_profile_host };

  if( !RelActAutoSolution::is_usable_status(solution.m_status) )
    return;

  const size_t num_curves = solution.m_options.rel_eff_curves.size();
  solution.m_profile_results.clear();

  // Automatic profiling is the single most expensive thing a solve can do - up to
  // `sm_profile_max_points_per_quantity` conditional optimizations per weak quantity - so it
  // belongs to the opt-in robust budget.  An *explicit* per-nuclide request is a different thing:
  // the user asked for that specific quantity, and it is honored whatever the budget.
  //
  // `solve_may_profile` is the same predicate `solve_ceres` used to decide whether to retain its
  // converged problem, and the two MUST agree: if they diverge the failure is silent either way -
  // a profile with no problem to pin within, or a needlessly retained problem after every ordinary
  // solve.
  const bool automatic_profiling_enabled = solution.m_options.robust_solve
                                        && solution.m_options.auto_profile_weak_mass_fractions;
  bool any_forced_profile = !solution.m_options.profile_targets.empty();
  for( const RelEffCurveInput &curve : solution.m_options.rel_eff_curves )
    for( const NucInputInfo &nuc : curve.nuclides )
      any_forced_profile = any_forced_profile || nuc.force_profile_mass_fraction;
  if( !automatic_profiling_enabled && !any_forced_profile )
    return;
  assert( RelActCalcAutoImp::solve_may_profile(solution.m_options) );

  // A profile is a conditional optimization of the selected physical problem.  Snapshot every
  // BR-nuisance interval (including an intentionally empty interval set) from the main solution;
  // conditional solves and any better-baseline reselection both consume these exact values.
  if( !solution.m_cost_functor )
  {
    solution.m_warnings.push_back( "Mass-fraction profiles were skipped because the usable solution"
                                   " did not retain its frozen objective." );
    return;
  }
  (void)input_drf;
  (void)all_peaks;
  // Public callers commonly leave peaks/DRF empty and let the primary solve derive both.  The
  // fitted output peaks are not that input: substituting them changes FWHM initialization and
  // nonlinear-calibration anchors.  Freeze the exact objects retained by the primary cost functor,
  // including an intentionally empty peak-search result.
  const std::vector<std::shared_ptr<const PeakDef>> frozen_profile_peaks
      = solution.m_cost_functor->m_all_peaks;
  const std::shared_ptr<const DetectorPeakResponse> frozen_profile_drf
      = solution.m_cost_functor->m_drf;
  if( !frozen_profile_drf )
  {
    solution.m_warnings.push_back( "Mass-fraction profiles were skipped because the usable solution"
                                   " did not retain its detector-response input." );
    return;
  }
  const std::vector<std::pair<double,double>> frozen_profile_peak_ranges
                                      = solution.m_cost_functor->m_peak_ranges_with_uncert;
  const std::uint64_t frozen_profile_gamma_hash = solution.m_frozen_gamma_membership_hash;

  // The winning frame's live `ceres::Problem`, already sitting at the unconstrained optimum.  With
  // it, a conditional point is one warm `ceres::Solve` on the same problem and the same parameter
  // buffer.  It is also the ONLY conditional mechanism: with no host (or an untrusted one, nulled
  // just below), every target reports a structured `Failed` profile rather than silently running
  // some slower engine - so a null host is handled here, never treated as an error.
  //
  // Two things are checked before it is trusted, because every delta-chi2 below is measured against
  // `solution.m_chi2` and a host built around a DIFFERENT point would shift all of them silently:
  // the problem must belong to the reported solution's functor, and the objective at its stored
  // optimum must be the objective the solution reports.
  std::shared_ptr<RelActCalcAutoImp::ProfileConditionalHost> profile_host = solution.m_profile_host;
  if( profile_host
      && (profile_host->cost_functor() != solution.m_cost_functor.get()) )
  {
    profile_host = nullptr;
    solution.m_warnings.push_back( "Profiles are unavailable because the retained optimizer"
                                   " problem did not belong to the reported solution." );
  }
  if( profile_host )
  {
    const double host_objective = profile_host->optimum_objective();
    const double tolerance = 1.0e-8*(1.0 + std::fabs(solution.m_chi2));
    if( !std::isfinite(host_objective) || !std::isfinite(solution.m_chi2)
        || (std::fabs(host_objective - solution.m_chi2) > tolerance) )
    {
      profile_host = nullptr;
      solution.m_warnings.push_back( "Profiles are unavailable because the retained optimizer"
          " problem was not at the reported objective ("
          + SpecUtils::printCompact(host_objective,12) + " vs "
          + SpecUtils::printCompact(solution.m_chi2,12) + ")." );
    }
  }//if( profile_host )

  // Backward elimination is part of the selected physical model, but its fixed mask is not
  // encoded by Options.  Snapshot it by canonical name and exact value so every conditional solve
  // (and any better-baseline reselection) uses precisely the same nuisance-parameter manifold.
  if( (solution.m_parameter_fixed_by_model_selection.size()
                                                != solution.m_final_parameters.size())
      || (solution.m_parameter_names.size() != solution.m_final_parameters.size()) )
  {
    solution.m_warnings.push_back( "Mass-fraction profiles were skipped because the selected-model"
                                   " parameter policy was incomplete." );
    return;
  }
  RelActCalcAutoImp::FrozenModelPolicy frozen_profile_model_policy;
  for( size_t index = 0; index < solution.m_final_parameters.size(); ++index )
    if( solution.m_parameter_fixed_by_model_selection[index] )
      frozen_profile_model_policy.push_back( {
          solution.m_parameter_names[index],solution.m_final_parameters[index] } );
  try
  {
    frozen_profile_model_policy = RelActCalcAutoImp::canonical_frozen_model_policy(
                                             std::move(frozen_profile_model_policy));
  }catch( const std::exception &e )
  {
    solution.m_warnings.push_back( std::string("Mass-fraction profiles were skipped because the")
        + " selected-model policy was invalid: " + e.what() );
    return;
  }
  const std::uint64_t frozen_profile_model_hash
      = RelActCalcAutoImp::frozen_model_policy_hash(frozen_profile_model_policy);
  if( solution.m_frozen_model_policy_hash != frozen_profile_model_hash )
  {
    solution.m_warnings.push_back( "Mass-fraction profiles were skipped because the retained"
                                   " selected-model policy hash did not match its parameters." );
    return;
  }

  const auto matches_frozen_profile_objective
      = [&]( const RelActAutoSolution &candidate ) -> bool {
    if( !candidate.m_cost_functor
        || (candidate.m_cost_functor->m_drf.get() != frozen_profile_drf.get())
        || (candidate.m_cost_functor->m_all_peaks.size() != frozen_profile_peaks.size()) )
      return false;
    for( size_t peak = 0; peak < frozen_profile_peaks.size(); ++peak )
      if( candidate.m_cost_functor->m_all_peaks[peak].get()
                                                  != frozen_profile_peaks[peak].get() )
        return false;
    return candidate.m_cost_functor
        && (candidate.m_cost_functor->m_peak_ranges_with_uncert
                                                == frozen_profile_peak_ranges)
        && (candidate.m_frozen_gamma_membership_hash == frozen_profile_gamma_hash)
        && (candidate.m_frozen_model_policy_hash == frozen_profile_model_hash);
  };

  // Do not let the first weak quantity visited choose a baseline restart on its own (each pass
  // performs at most one reselection, from the best deferred discovery).  A
  // different quantity can enter a still lower basin (the Pu free-age problem is a concrete
  // example), and source/caller order must not decide which of those seeds reaches final candidate
  // selection.  Collect one independently evaluated discovery per interrupted profile, rank all
  // of them by the frozen physical objective and a semantic key, and perform exactly one
  // unconstrained reselection after the complete first pass.
  struct DeferredBaselineDiscovery
  {
    RelActCalcAutoImp::ProfileLikelihood::PendingBaselineDiscovery rank;
    std::unique_ptr<RelActAutoSolution> conditional_solution;
    /** The discovering quantity's full descriptor, so a failed reselection can mark its slot and
     name the quantity correctly (never "mass-fraction" for an activity or ratio target). */
    Options::ProfileTarget target;
    /** Why the discovering quantity was being profiled, so a failure entry created here does not
     default a forced request to `AutomaticWeak`. */
    ProfileReason reason = ProfileReason::AutomaticWeak;
    size_t curve_index = 0;
    std::string target_symbol;
  };
  // The reselection trail: a restarted solution is a fresh solve carrying only its own warnings,
  // so the hop history has to be carried across every adoption or a multi-hop chain reports only
  // its final hop.  CONTRACT: the trail is selected by this exact prefix, so it exists ONCE here
  // and every hop message is composed through `push_reselection_warning` - there is no second
  // place to keep in sync.
  static const char * const sm_reselection_prefix = "Baseline reselection ";
  const auto adopt_restarted_solution = [&solution]( RelActAutoSolution &&restarted ){
    std::vector<std::string> reselection_trail;
    for( const std::string &warning : solution.m_warnings )
      if( warning.rfind(sm_reselection_prefix,0) == 0 )
        reselection_trail.push_back( warning );
    solution = std::move(restarted);
    for( const std::string &warning : reselection_trail )
      solution.m_warnings.push_back( warning );
  };
  const auto push_reselection_warning = [&solution]( const std::string &text ){
    solution.m_warnings.push_back( std::string(sm_reselection_prefix) + text );
  };

  std::vector<DeferredBaselineDiscovery> deferred_baseline_discoveries;

  // Profile results are keyed by their full target descriptor (kind + source + curve).  The
  // returned reference is valid only until the next call - the vector may reallocate.
  const auto profile_slot = [&solution]( const Options::ProfileTarget &slot_target )
                                                                          -> Profile & {
    for( RelActAutoSolution::ProfileResultEntry &entry : solution.m_profile_results )
      if( entry.target == slot_target )
        return entry.profile;
    solution.m_profile_results.push_back( RelActAutoSolution::ProfileResultEntry{} );
    solution.m_profile_results.back().target = slot_target;
    return solution.m_profile_results.back().profile;
  };
  const auto mass_fraction_target = []( const size_t curve, const SrcVariant &src ) {
    Options::ProfileTarget slot_target;
    slot_target.kind = Options::ProfileTarget::Kind::MassFraction;
    slot_target.source = src;
    slot_target.rel_eff_curve_index = curve;
    return slot_target;
  };
  // Warn only once per pass when the reselection budget was exhausted (see the RejectAfterRestart
  // branch): the taint applies to the whole solution, not per discovering target.
  bool budget_exhaustion_warned = false;

  // Fixed for the whole pass: it depends only on the reported objective and covariance scale.
  const double baseline_improvement_tolerance
      = RelActCalcAutoImp::ProfileLikelihood::baseline_improvement_tolerance(
                                          solution.m_chi2,solution.m_cov_scale);

  /** Records a better-baseline discovery for the deferred reselection; false when the discovery
   could not be turned into a usable seed (the caller then reports that as a structured Failed). */
  const auto record_discovery = [&]( const ScanSetup &setup,
                                     std::unique_ptr<RelActAutoSolution> &better ) -> bool {
    if( !better || !std::isfinite(better->m_chi2) )
      return false;
    DeferredBaselineDiscovery deferred;
    deferred.rank.full_objective = better->m_chi2;
    deferred.rank.semantic_key = setup.semantic_key;
    deferred.conditional_solution = std::move(better);
    deferred.target = setup.target;
    deferred.reason = setup.reason;
    deferred.curve_index = setup.curve_index;
    deferred.target_symbol = setup.stats_name;
    deferred_baseline_discoveries.push_back( std::move(deferred) );
    return true;
  };//record_discovery

  /** One target's conditional work and all of its reporting: the pin, the evaluator, either the
   two pre-scan probes or the full scan, and the result entry.  Shared by every
   `ProfileTarget::Kind` - see `ScanSetup` for what each lane must decide first.

   Assumes the caller has already installed the chart and holds a guard that unwinds it.

   @returns `AbortPass` when the user canceled, which the caller must propagate by returning from
            the pass immediately - a cancel that is merely `continue`d keeps spending solves. */
  const auto run_one_target = [&]( const ScanSetup &setup, const bool probe_only ) -> TargetOutcome
  {
    using Eval = RelActCalcAutoImp::ProfileLikelihood::Evaluation;
    using EvalStatus = RelActCalcAutoImp::ProfileLikelihood::EvaluationStatus;
    using PLDirection = RelActCalcAutoImp::ProfileLikelihood::Direction;
    using ScanStatus = RelActCalcAutoImp::ProfileLikelihood::ScanStatus;
    using RelActCalcAutoImp::ProfileLikelihood::has_real_bound;

    assert( profile_host );
    assert( !setup.use_level_set_refine
            || (setup.target.kind == Options::ProfileTarget::Kind::MassFraction) );

    Profile profile;
    profile.reason = setup.reason;

    const size_t host_solves_before = profile_host->num_conditional_solves();
    const std::chrono::steady_clock::time_point profile_start
                                                  = std::chrono::steady_clock::now();
    // Conditional solves run the SAME cost functor as the main fit, so its own counters give the
    // profile phase's share directly: `evals` is full-problem residual evaluations (double and
    // Jet alike) and `us_in_eval` the time inside them, which together answer what a conditional
    // evaluation costs relative to a main-fit one - the chart decodes installed for a scan sit on
    // that path, so the two are not assumed equal.
    const RelActCalcAutoImp::RelActAutoCostFcn * const eval_counter = profile_host->cost_functor();
    const size_t evals_before = eval_counter ? eval_counter->m_ncalls.load() : 0;
    const size_t eval_ns_before
        = eval_counter ? eval_counter->m_nanoseconds_spent_in_eval.load() : 0;

    // The warm chains must be copied AFTER the chart is installed: `optimum()` returns a reference
    // into the live chart's re-encoded vector while a reparameterization is engaged.
    std::vector<double> lower_warm = profile_host->optimum(), upper_warm = lower_warm;
    std::unique_ptr<RelActAutoSolution> better_conditional;

    if( RelActCalcAutoImp::profile_stats_enabled() )
      std::cerr << "profile-scan " << setup.stats_name << " curve=" << setup.curve_index
                << " index=" << setup.slot << " a0=" << setup.baseline_parameter
                << " box=[" << setup.parameter_box.first << "," << setup.parameter_box.second << "]"
                << " " << setup.dq_key << "=" << setup.chart_scale
                << " q0=" << setup.nominal_reported
                << (probe_only ? " mode=probe" : "") << std::endl;

    const auto evaluate = [&]( const double requested_parameter,
                               const PLDirection direction,
                               const size_t ) -> Eval {
      Eval evaluation;
      if( cancel_calc && cancel_calc->load() )
      {
        evaluation.status = EvalStatus::Canceled;
        evaluation.diagnostic = setup.cancel_message;
        return evaluation;
      }

      // Points are visited outward from the optimum, so the previous point of this direction is
      // always the nearest one and is chronologically available - which is what makes the
      // scanner's proximity-based warm-start policy unnecessary here.
      std::vector<double> &warm = (direction == PLDirection::Lower) ? lower_warm : upper_warm;

      const RelActCalcAutoImp::ProfileConditionalHost::PinnedResult point
          = profile_host->solve_pinned( setup.slot, requested_parameter, warm );
      if( RelActCalcAutoImp::profile_stats_enabled() )
        std::cerr << "profile-point " << setup.stats_name
                  << " dir=" << ((direction == PLDirection::Lower) ? "lo" : "hi")
                  << " req=" << requested_parameter
                  << " got=" << point.achieved_parameter
                  << " conv=" << point.converged
                  << " chi2=" << point.physical_chi2
                  << " " << point.diagnostic << std::endl;

      if( point.canceled )
      {
        evaluation.status = EvalStatus::Canceled;
        evaluation.diagnostic = setup.cancel_message;
        return evaluation;
      }
      if( !point.converged )
      {
        evaluation.status = EvalStatus::Failed;
        evaluation.diagnostic = point.diagnostic;
        return evaluation;
      }

      // The pinned point's own reported value: the sample coordinate, and the level the
      // refinement below holds fixed.
      double pinned_value = std::numeric_limits<double>::quiet_NaN();
      try
      {
        pinned_value = profile_host->reported_quantity( setup.target, point.parameters );
      }catch( const std::exception &e )
      {
        evaluation.status = EvalStatus::Failed;
        evaluation.diagnostic = std::string("A conditional point could not be read back: ")
                                + e.what();
        return evaluation;
      }

      // Correct the one known bias of pinning, by EVALUATION rather than by subtracting an
      // estimate.  See `refine_on_level_set`: the refined point is still a rigorous upper bound on
      // the exact profile at its own achieved `q`, so this can only move a sample toward the truth
      // from the narrow side and can never widen the interval.
      //
      // Skipped near the vertex: there the gap is a small fraction of an already-small delta chi2,
      // the recoverable amount sits near the objective's own reproducibility floor, and those
      // samples do not move the crossing anyway.
      std::vector<double> sample_parameters = point.parameters;
      double sample_chi2 = point.physical_chi2;
      const double refine_floor = 0.3*setup.thresholds[0];
      if( setup.use_level_set_refine && std::isfinite(pinned_value)
          && ((point.physical_chi2 - profile_host->optimum_objective()) >= refine_floor) )
      {
        const RelActCalcAutoImp::ProfileConditionalHost::LevelSetStep refined
            = profile_host->refine_on_level_set( point.parameters, point.physical_chi2,
                                                 pinned_value, setup.nominal_reported,
                                                 solution.m_covariance,
                                                 setup.target.source, setup.curve_index );
        if( refined.improved )
        {
          // Recovery is ~0 at 95-98% of corpus points; the rare 1e5-scale values mark conditional
          // solves that stalled far above their true minima at far-out requests - the stall
          // population this line exists to measure.
          if( RelActCalcAutoImp::profile_stats_enabled() )
            std::cerr << "profile-refine " << setup.stats_name
                      << " recovered=" << (point.physical_chi2 - refined.physical_chi2)
                      << " at_delta="
                      << (point.physical_chi2 - profile_host->optimum_objective())
                      << std::endl;
          sample_parameters = refined.parameters;
          sample_chi2 = refined.physical_chi2;
        }
      }//if( refinement is worth attempting )

      // Only a refinement can have moved the point, and re-reading an unmoved vector would just
      // pay for the same answer twice.
      double achieved_value = pinned_value;
      if( sample_parameters != point.parameters )
      {
        try
        {
          achieved_value = profile_host->reported_quantity( setup.target, sample_parameters );
        }catch( const std::exception &e )
        {
          evaluation.status = EvalStatus::Failed;
          evaluation.diagnostic = std::string("A conditional point could not be read back: ")
                                  + e.what();
          return evaluation;
        }

        // The refinement holds `q` fixed only to first order, so a refined point lands at its own
        // slightly different reported value.  That is legitimate - the sample is quoted at
        // whatever it achieved - but a LARGE drift could reorder samples and make the hygiene pass
        // see a fold where the map never folded, truncating the side.  Keep the original point
        // rather than risk that.
        const double drift = std::fabs( achieved_value - pinned_value );
        if( drift > 0.1*std::fabs(pinned_value - setup.nominal_reported) )
        {
          sample_parameters = point.parameters;
          sample_chi2 = point.physical_chi2;
          achieved_value = pinned_value;
        }
      }//if( the refinement moved the point )

      // Checked on whatever value is actually going to be quoted, so a reverted refinement cannot
      // ship a value that was never itself checked.
      if( !std::isfinite(achieved_value)
          || (achieved_value < setup.domain_lower) || (achieved_value > setup.domain_upper) )
      {
        evaluation.status = EvalStatus::Failed;
        evaluation.diagnostic = setup.domain_message;
        return evaluation;
      }

      const double delta = sample_chi2 - profile_host->optimum_objective();

      // A conditional is a restricted minimization, so it cannot beat the unconstrained optimum.
      // One that does means the reported baseline was not the optimum; hand it to the deferred
      // rebaseline flow rather than reporting a negative delta-chi2.
      if( delta < -baseline_improvement_tolerance )
      {
        try
        {
          better_conditional = std::make_unique<RelActAutoSolution>(
                  profile_host->warm_start_solution(solution,sample_parameters,sample_chi2) );
        }catch( const std::exception & )
        {
          better_conditional.reset();
        }
        evaluation.status = EvalStatus::BetterBaseline;
        evaluation.diagnostic = "A conditional profile point reached a lower physical objective"
                                " than the reported solution.";
        return evaluation;
      }

      // Deliberately the PINNED point, never the refined one: a refined point is not a conditional
      // minimum of anything, and the outward sweep's warm-start chain - and the exact
      // `x_optimum - x` chord the refinement itself relies on - are only meaningful between
      // genuine conditional minima.
      warm = point.parameters;

      evaluation.status = EvalStatus::Success;
      evaluation.reported_fraction = achieved_value;
      evaluation.delta_chi2 = (std::max)( 0.0, delta );

      // The honest bound test: a genuine Ceres box edge, or a reported value that has arrived at
      // its feasible limit.  Deliberately the RAW box - a synthetic scan cap is not a bound (see
      // `ScanSetup`), and quoting one as if it were would assert a limit that does not exist.
      const bool at_box_edge
          = (has_real_bound(setup.parameter_box.first)
               && (point.achieved_parameter <= setup.parameter_box.first))
            || (has_real_bound(setup.parameter_box.second)
               && (point.achieved_parameter >= setup.parameter_box.second));
      const bool at_reported_limit
          = (achieved_value <= (setup.feasible_lower + setup.feasible_tolerance))
            || (achieved_value >= (setup.feasible_upper - setup.feasible_tolerance));
      evaluation.reached_feasible_bound = at_box_edge || at_reported_limit;
      return evaluation;
    };//evaluate lambda

    // --- The pre-scan basin check ---------------------------------------------------------
    //
    // One conditional solve per side, at the scan's own opening displacement, so a probe lands
    // exactly where the scan's first point would - the premise the measured hit rate rests on.
    // A probe landing below the optimum is the basin-miss signature; the evaluator routes it into
    // `better_conditional` exactly as a scan point would.  Probes write no per-target result: all
    // reporting belongs to the scan pass (the cancel status below is the one exception).
    if( probe_only )
    {
      using RelActCalcAutoImp::ProfileLikelihood::opening_probe_distance;
      const double steps[2] = {
        opening_probe_distance( setup.parameter_one_sigma,
                                setup.baseline_parameter - setup.scan_lower ),
        opening_probe_distance( setup.parameter_one_sigma,
                                setup.scan_upper - setup.baseline_parameter ) };
      for( int side = 0; side < 2; ++side )
      {
        if( !(steps[side] > 0.0) )
          continue;
        const double request = setup.baseline_parameter + (side ? steps[side] : -steps[side]);
        const Eval eval = evaluate( request,
                                    side ? PLDirection::Upper : PLDirection::Lower, 0 );
        if( eval.status == EvalStatus::Canceled )
        {
          solution.m_status = RelActAutoSolution::Status::UserCanceled;
          solution.m_error_message = setup.cancel_message;
          return TargetOutcome::AbortPass;
        }
        if( eval.status == EvalStatus::BetterBaseline )
        {
          if( RelActCalcAutoImp::profile_stats_enabled() )
            std::cerr << "profile-probe " << setup.stats_name << " curve=" << setup.curve_index
                      << " side=" << (side ? "hi" : "lo") << " discovery=1" << std::endl;
          record_discovery( setup, better_conditional );
          break;
        }
      }//for( each probe side )
      return TargetOutcome::Continue;
    }//if( probe_only )

    // --- The scan, and everything it reports ----------------------------------------------
    const RelActCalcAutoImp::ProfileLikelihood::ScanResult scan
        = RelActCalcAutoImp::ProfileLikelihood::fit_profile(
              setup.baseline_parameter, setup.nominal_reported,
              setup.scan_lower, setup.scan_upper, setup.thresholds,
              RelActCalcAutoImp::ProfileLikelihood::sm_profile_max_points_per_quantity,
              evaluate, setup.parameter_one_sigma );

    profile.num_fits = scan.num_fits;
    for( const RelActCalcAutoImp::ProfileLikelihood::Sample &sample : scan.samples )
      profile.sampled_delta_chi2.emplace_back(sample.reported_fraction,sample.delta_chi2);

    const auto report_failure = [&]( const std::string &message, const bool plain_prefix ){
      profile.status = ProfileStatus::Failed;
      profile.message = message;
      profile_slot( setup.target ) = profile;
      solution.m_warnings.push_back( setup.stats_name
          + (plain_prefix ? setup.plain_warning_suffix : setup.failed_warning_suffix) + message );
    };

    if( scan.status == ScanStatus::Canceled )
    {
      solution.m_status = RelActAutoSolution::Status::UserCanceled;
      solution.m_error_message = scan.diagnostic.empty() ? setup.cancel_message : scan.diagnostic;
      return TargetOutcome::AbortPass;
    }

    if( scan.status == ScanStatus::BetterBaseline )
    {
      const auto disposition
          = RelActCalcAutoImp::ProfileLikelihood::baseline_discovery_disposition(
                                                      baseline_restart_count);
      if( disposition
          == RelActCalcAutoImp::ProfileLikelihood::BaselineDiscoveryDisposition::RejectAfterRestart )
      {
        report_failure( "A profile found a still better point after the baseline-reselection"
                        " budget was exhausted; the selected baseline is still not optimal."
                        + (scan.diagnostic.empty() ? std::string()
                                                   : std::string(" ") + scan.diagnostic), true );
        // The proof of non-optimality taints every interval of this pass, not just the discovering
        // target's - their delta-chi2 baselines are the same nominal.  Say so once.
        if( !budget_exhaustion_warned )
        {
          budget_exhaustion_warned = true;
          solution.m_warnings.push_back( "The reported baseline is known not to be the global"
              " optimum (the baseline-reselection budget was exhausted while profiles kept"
              " finding better points), so every profile interval of this solution is measured"
              " against a non-optimal nominal and may be unreliable." );
        }
        if( RelActCalcAutoImp::profile_stats_enabled() )
          std::cerr << "profile-rebaseline-exhausted target=" << setup.stats_name
                    << " hops=" << baseline_restart_count << std::endl;
        return TargetOutcome::Continue;
      }//if( the reselection budget is exhausted )

      if( !record_discovery( setup, better_conditional ) )
      {
        report_failure( "A better profile point was reported but could not seed baseline"
                        " reselection.", true );
        return TargetOutcome::Continue;
      }

      // This entry is temporary: a successful reselection replaces the solution and restarts all
      // profiles.  If reselection itself fails, retaining a structured Failed result is more
      // useful than silently dropping the quantity which discovered the better point.
      profile.status = ProfileStatus::Failed;
      profile.message = "A better physical point is pending final candidate selection.";
      profile_slot( setup.target ) = profile;
      return TargetOutcome::Continue;
    }//if( a scan point beat the reported optimum )

    if( (scan.status == ScanStatus::Failed) || (scan.status == ScanStatus::FitCapReached) )
    {
      report_failure( scan.diagnostic.empty()
                      ? std::string("The bounded profile could not be classified.")
                      : scan.diagnostic, false );
      return TargetOutcome::Continue;
    }

    const double confidences[2] = { 0.6827, 0.95 };
    for( size_t level = 0; level < 2; ++level )
    {
      Interval interval;
      interval.confidence_level = confidences[level];
      interval.delta_chi2 = setup.thresholds[level];
      interval.lower = scan.intervals[level].lower.reported_fraction;
      interval.upper = scan.intervals[level].upper.reported_fraction;
      interval.lower_kind = scan.intervals[level].lower.likelihood_crossing
                            ? EndpointKind::LikelihoodCrossing : setup.lower_bound_kind;
      interval.upper_kind = scan.intervals[level].upper.likelihood_crossing
                            ? EndpointKind::LikelihoodCrossing : setup.upper_bound_kind;
      if( interval.lower > interval.upper )
        std::swap(interval.lower, interval.upper);
      profile.intervals.push_back(interval);
    }

    if( scan.status == ScanStatus::NonIdentifiable )
      profile.status = ProfileStatus::NonIdentifiable;
    else if( scan.status == ScanStatus::BoundaryLimited )
      profile.status = ProfileStatus::BoundaryLimited;
    else
      profile.status = ProfileStatus::Complete;
    profile.message = scan.diagnostic;

    // GAUGE CAVEAT (same text as `Options::ProfileTarget::Kind::RelativeActivity`): an individual
    // relative activity is defined only against the curve's normalization, so the interval must be
    // presented per gauge, never as gauge-free.  Mass fractions, same-curve ratios, and ages are
    // gauge invariant and need no note.
    if( setup.target.kind == Options::ProfileTarget::Kind::RelativeActivity )
      profile.message += (profile.message.empty() ? "" : "  ")
          + std::string("Interval is quoted in this fit's relative-activity gauge (the curve"
                        " normalization convention): differences and relative widths are"
                        " meaningful; the absolute value is not gauge-free.");

    profile_slot( setup.target ) = profile;

    if( RelActCalcAutoImp::profile_stats_enabled() )
    {
      std::cerr << "profile-stats " << setup.stats_name << " curve=" << setup.curve_index
                << " points=" << profile.num_fits
                << " ceres_solves=" << (profile_host->num_conditional_solves() - host_solves_before)
                << " seconds=" << std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - profile_start).count()
                << " evals=" << (eval_counter ? (eval_counter->m_ncalls.load() - evals_before) : 0)
                << " us_in_eval="
                << (eval_counter ? ((eval_counter->m_nanoseconds_spent_in_eval.load()
                                     - eval_ns_before)/1000) : 0)
                << " status=" << static_cast<int>(profile.status) << std::endl;
      if( setup.emit_interval_stats && (profile.intervals.size() >= 2) )
        std::cerr << "profile-intervals " << setup.stats_name
                  << " q0=" << setup.nominal_reported
                  << " l68=" << profile.intervals[0].lower
                  << " u68=" << profile.intervals[0].upper
                  << " l95=" << profile.intervals[1].lower
                  << " u95=" << profile.intervals[1].upper << std::endl;
    }

    return TargetOutcome::Continue;
  };//run_one_target

  // One pass over every would-be profile target.  `probe_only` is the pre-scan basin check: the
  // same eligibility, control window, carrier install, and evaluator as a scan, but only the two
  // opening probes are solved and NOTHING is written to the solution - discoveries are collected
  // into `deferred_baseline_discoveries` and every other side effect belongs to the scan pass.
  const auto run_target_pass = [&]( const bool probe_only ) -> void
  {
    for( size_t curve_index = 0; curve_index < num_curves; ++curve_index )
    {
      const RelEffCurveInput &base_curve = solution.m_options.rel_eff_curves[curve_index];

      // Every profile target is a fitted source with a parameter of its own.  Correlation-generated
      // Pu-242 is deliberately NOT among them: it is derived from the other Pu isotopes rather than
      // fitted, so it has no likelihood direction and a profile interval over it would not describe
      // anything this data constrains.  Its reporting is unaffected - the correlation still supplies
      // its value and still renormalizes its siblings - it simply is never a target here, and
      // `ProfileTarget::why_not_usable()` rejects it with that reasoning.
      struct TargetNuclide
      {
        const SandiaDecay::Nuclide *nuclide = nullptr;
        const NucInputInfo *input = nullptr;   //!< never null
      };
      // INVARIANT: only nuclides EXPLICITLY in the input list are ever profiled - in every
      // conditional engine.  In particular, correlation-generated Pu-242 is never a target: users
      // accept that it is known only by correlation, and care about it only as a small correction
      // to the enrichment of the modeled isotopes.
      std::vector<TargetNuclide> profile_targets;
      profile_targets.reserve( base_curve.nuclides.size() );
      for( const NucInputInfo &input : base_curve.nuclides )
      {
        const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(input.source);
        if( nuclide )
          profile_targets.push_back( {nuclide,&input} );
      }

      for( const TargetNuclide &profile_target : profile_targets )
      {
        const SandiaDecay::Nuclide * const target = profile_target.nuclide;
        if( !target )
          continue;
        const SrcVariant target_source(target);
        assert( profile_target.input );
        // An explicit `Options::profile_targets` mass-fraction entry is the same request as the
        // per-nuclide checkbox sugar; either forces this quantity.
        const bool forced_by_target = std::any_of(
            begin(solution.m_options.profile_targets), end(solution.m_options.profile_targets),
            [&target_source,curve_index]( const Options::ProfileTarget &request ){
              return (request.kind == Options::ProfileTarget::Kind::MassFraction)
                     && (request.source == target_source)
                     && (request.rel_eff_curve_index == curve_index);
            } );
        const bool forced = profile_target.input->force_profile_mass_fraction || forced_by_target;

        // Which quantities carry information about this source (see
        // `Options::ProfileTarget::Kind`).  Only two shapes short-circuit here; everything else -
        // in particular a CONSTRAINT-BLOCKED mass fraction - deliberately falls through to the
        // ordinary flow below, where the weak gate decides whether it is reported at all and
        // `install_carrier_reparam` refuses it with the reason.  That refusal is the structured
        // Failed both forced and automatically-weak targets have always received; short-circuiting
        // them here on the (newer, options-level) eligibility check would silently drop the
        // automatic ones, leaving only the covariance band that automatic profiling exists to
        // replace, with nothing saying a profile was refused.
        const std::vector<Options::ProfileTarget::Kind> profilable_kinds
            = RelActCalcAuto::profilable_quantity_kinds( solution.m_options,target_source,curve_index );

        // Pu-242 under an active correlation is derived rather than fitted, so it is never a
        // target and never reported as a failed one (see the INVARIANT above).
        if( (target->atomicNumber == 94) && (target->massNumber == 242)
            && (base_curve.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable) )
          continue;

        // A SOLE-ISOTOPE element has no mass fraction to profile at all - the normalized fraction
        // is identically one - so a forced request is honored through the quantity that does carry
        // information, its relative activity (synthesized before the kinds loop below), or through
        // the boundary placeholder when there is no activity to scan either.
        if( RelActCalcAuto::same_element_input_count(base_curve,target) < 2 )
        {
          const bool activity_would_carry_it
              = std::find( begin(profilable_kinds),end(profilable_kinds),
                           Options::ProfileTarget::Kind::RelativeActivity )
                != end(profilable_kinds);
          if( forced && !probe_only && !activity_would_carry_it )
          {
            Profile profile;
            profile.reason = ProfileReason::Forced;
            profile.status = ProfileStatus::BoundaryLimited;
            profile.num_fits = 0;
            profile.message = "This is the only modeled isotope of its element; its normalized mass"
                              " fraction is fixed at the physical limit.";
            Interval p68;
            p68.confidence_level = 0.6827;
            p68.delta_chi2 = solution.m_cov_scale;
            p68.lower = p68.upper = 1.0;
            p68.lower_kind = p68.upper_kind = EndpointKind::PhysicalLimit;
            Interval p95 = p68;
            p95.confidence_level = 0.95;
            p95.delta_chi2 = solution.m_cov_scale*3.841458820694124;
            profile.intervals = {p68,p95};
            profile_slot( mass_fraction_target(curve_index,target_source) ) = std::move(profile);
          }
          continue;
        }//if( sole isotope of its element )

        RelActAutoSolution::MassFractionResult nominal;
        try
        {
          nominal = solution.mass_enrichment_result( target, curve_index );
        }catch( const std::exception &e )
        {
          // Do not silently lose an explicit request when even the structured nominal result cannot
          // be formed: retain the otherwise usable main fit and expose a failed profile.
          if( forced && !probe_only )
          {
            Profile failed;
            failed.status = ProfileStatus::Failed;
            failed.reason = ProfileReason::Forced;
            failed.message = std::string("Could not initialize profile: ") + e.what();
            profile_slot( mass_fraction_target(curve_index,target_source) ) = failed;
            solution.m_warnings.push_back( target->symbol + " mass-fraction profile failed: " + e.what() );
          }
          continue;
        }

        const bool has_pu_reporting_transform
            = (target->atomicNumber == 94)
              && (base_curve.pu242_correlation_method
                  != RelActCalc::PuCorrMethod::NotApplicable);
        double control_baseline = has_pu_reporting_transform
                                ? nominal.fraction
                                : nominal.uncorrected_fraction.value_or(nominal.fraction);
        double control_lower = 0.0;
        double control_upper = 1.0;
        EndpointKind lower_bound_kind = EndpointKind::PhysicalLimit;
        EndpointKind upper_bound_kind = EndpointKind::PhysicalLimit;
        if( !has_pu_reporting_transform )
        {
          // Project every same-element input window through the unit simplex.  A sibling lower bound
          // limits this target from above; sibling upper bounds can limit it from below when all
          // remaining isotopes are constrained.  Without this projection an infeasible conditional
          // solve was previously reported as a profile failure instead of an input-bound endpoint.
          std::vector<std::pair<double,double>> sibling_bounds;
          for( const NucInputInfo &sibling_info : base_curve.nuclides )
          {
            const SandiaDecay::Nuclide * const sibling
                = RelActCalcAuto::nuclide(sibling_info.source);
            if( !sibling || (sibling == target)
                || (sibling->atomicNumber != target->atomicNumber) )
              continue;

            const auto sibling_constraint = std::find_if(
                begin(base_curve.mass_fraction_constraints),
                end(base_curve.mass_fraction_constraints),
                [sibling]( const RelEffCurveInput::MassFractionConstraint &constraint ){
                  return constraint.nuclide == sibling;
                } );
            if( sibling_constraint == end(base_curve.mass_fraction_constraints) )
            {
              // An unconstrained sibling may carry any positive remainder.
              sibling_bounds.emplace_back(0.0,1.0);
            }else
            {
              sibling_bounds.emplace_back(sibling_constraint->lower_mass_fraction,
                                          sibling_constraint->upper_mass_fraction);
            }
          }

          for( const RelEffCurveInput::MassFractionConstraint &constraint
               : base_curve.mass_fraction_constraints )
          {
            if( constraint.nuclide != target )
              continue;
            control_lower = constraint.lower_mass_fraction;
            control_upper = constraint.upper_mass_fraction;
            lower_bound_kind = EndpointKind::InputConstraintLimit;
            upper_bound_kind = EndpointKind::InputConstraintLimit;
            break;
          }

          const std::pair<double,double> projected
              = RelActCalcAutoImp::ProfileLikelihood::project_simplex_component_bounds(
                  control_lower,control_upper,sibling_bounds );
          if( projected.first > control_lower + 1.0e-14 )
          {
            control_lower = projected.first;
            lower_bound_kind = EndpointKind::InputConstraintLimit;
          }
          if( projected.second < control_upper - 1.0e-14 )
          {
            control_upper = projected.second;
            upper_bound_kind = EndpointKind::InputConstraintLimit;
          }

          // Collapse every ordinary same-element activity (including fixed activity-ratio chains)
          // into its independent root coordinate.  Mass-fraction-parameterized siblings consume a
          // separately bounded share of the simplex; the ordinary groups share the remainder.  The
          // selected ordinary targets fraction is monotone in each root activity, so activity boxes
          // have an exact corner projection and conditional scans never ask Ceres to cross them.
          struct RootMassGroup
          {
            SrcVariant root;
            RelActCalcAutoImp::ProfileLikelihood::ActivityBoxMassGroup mass;
          };
          std::vector<RootMassGroup> root_groups;
          double constrained_fraction_lower_sum = 0.0;
          double constrained_fraction_upper_sum = 0.0;
          for( const NucInputInfo &member_info : base_curve.nuclides )
          {
            const SandiaDecay::Nuclide * const member
                = RelActCalcAuto::nuclide(member_info.source);
            if( !member || (member->atomicNumber != target->atomicNumber) )
              continue;

            const auto mass_constraint = std::find_if(
                begin(base_curve.mass_fraction_constraints),
                end(base_curve.mass_fraction_constraints),
                [member]( const RelEffCurveInput::MassFractionConstraint &constraint ){
                  return constraint.nuclide == member;
                } );
            if( mass_constraint != end(base_curve.mass_fraction_constraints) )
            {
              constrained_fraction_lower_sum += mass_constraint->lower_mass_fraction;
              constrained_fraction_upper_sum += mass_constraint->upper_mass_fraction;
              continue;
            }

            SrcVariant root = member_info.source;
            double activity_multiple = 1.0;
            for( size_t step = 0; step < base_curve.act_ratio_constraints.size(); ++step )
            {
              const auto controlling = std::find_if(
                  begin(base_curve.act_ratio_constraints),end(base_curve.act_ratio_constraints),
                  [&root]( const RelEffCurveInput::ActRatioConstraint &constraint ){
                    return constraint.constrained_source == root;
                  } );
              if( controlling == end(base_curve.act_ratio_constraints) )
                break;
              activity_multiple *= controlling->constrained_to_controlled_activity_ratio;
              root = controlling->controlling_source;
            }

            const double mass_coefficient = activity_multiple / member->activityPerGram();
            if( !(std::isfinite(mass_coefficient) && (mass_coefficient > 0.0)) )
              continue;
            auto group = std::find_if( begin(root_groups),end(root_groups),
                [&root]( const RootMassGroup &candidate ){ return candidate.root == root; } );
            if( group == end(root_groups) )
            {
              RootMassGroup added;
              added.root = root;
              const auto root_input = std::find_if(
                  begin(base_curve.nuclides),end(base_curve.nuclides),
                  [&root]( const NucInputInfo &candidate ){ return candidate.source == root; } );
              if( root_input != end(base_curve.nuclides) )
              {
                added.mass.lower_root_activity
                    = root_input->min_rel_act.value_or(0.0);
                added.mass.upper_root_activity = root_input->max_rel_act.value_or(
                    std::numeric_limits<double>::infinity() );
              }
              root_groups.push_back(std::move(added));
              group = std::prev(end(root_groups));
            }
            group->mass.total_mass_per_root_activity += mass_coefficient;
            if( member == target )
              group->mass.target_mass_per_root_activity += mass_coefficient;
          }

          std::vector<RelActCalcAutoImp::ProfileLikelihood::ActivityBoxMassGroup>
              activity_mass_groups;
          activity_mass_groups.reserve(root_groups.size());
          for( const RootMassGroup &group : root_groups )
            activity_mass_groups.push_back(group.mass);
          const auto ordinary_composition_bounds
              = RelActCalcAutoImp::ProfileLikelihood::activity_box_fraction_bounds(
                                                           activity_mass_groups);
          if( ordinary_composition_bounds )
          {
            // The sigma-block can independently vary its constrained total over the sum of the
            // member windows.  The ordinary activity groups occupy the complementary mass.
            const double ordinary_fraction_lower = (std::max)(0.0,
                1.0-(std::min)(1.0,constrained_fraction_upper_sum));
            const double ordinary_fraction_upper = (std::min)(1.0,
                (std::max)(0.0,1.0-constrained_fraction_lower_sum));
            const double activity_lower
                = ordinary_fraction_lower*ordinary_composition_bounds->first;
            const double activity_upper
                = ordinary_fraction_upper*ordinary_composition_bounds->second;
            if( activity_lower > control_lower + 1.0e-14 )
            {
              control_lower = activity_lower;
              lower_bound_kind = EndpointKind::InputConstraintLimit;
            }
            if( activity_upper < control_upper - 1.0e-14 )
            {
              control_upper = activity_upper;
              upper_bound_kind = EndpointKind::InputConstraintLimit;
            }
          }
        }

        // Weakness is quantity- and input-domain-specific.  The structured nominal accessor can
        // classify coverage of the physical [0,1] range, but a target may have a much narrower
        // interval after its own window, sibling-simplex projection, and ratio-chain structure are
        // applied.  A local 95% band covering that actual interval is just as non-identifying and
        // must trigger the bounded profile.
        const bool spans_actual_feasible_range = nominal.covariance_one_sigma
            && RelActCalcAutoImp::ProfileLikelihood::local_gaussian_band_spans_feasible_range(
                 control_baseline,*nominal.covariance_one_sigma,control_lower,control_upper );
        const bool weak = (nominal.covariance_quality == CovQuality::Unavailable)
                          || (nominal.covariance_quality == CovQuality::LocallyUnreliable)
                          || (nominal.covariance_quality == CovQuality::SpansFeasibleRange)
                          || spans_actual_feasible_range;
        if( !forced && (!automatic_profiling_enabled || !weak) )
          continue;

        Profile profile;
        profile.reason = forced ? ProfileReason::Forced : ProfileReason::AutomaticWeak;

        control_baseline = (std::max)(control_lower, (std::min)(control_upper, control_baseline));
        if( !(std::isfinite(control_baseline) && std::isfinite(solution.m_chi2)) )
        {
          profile.status = ProfileStatus::Failed;
          profile.message = "The nominal fraction or objective was non-finite.";
          if( !probe_only )
          {
            profile_slot( mass_fraction_target(curve_index,target_source) ) = profile;
            solution.m_warnings.push_back( target->symbol + " mass-fraction profile failed: "
                                           + profile.message );
          }
          continue;
        }

        if( std::fabs(control_upper - control_lower) <= 1.0e-14 )
        {
          Interval p68;
          p68.confidence_level = 0.6827;
          p68.delta_chi2 = solution.m_cov_scale;
          p68.lower = p68.upper = nominal.fraction;
          p68.lower_kind = p68.upper_kind = EndpointKind::InputConstraintLimit;
          Interval p95 = p68;
          p95.confidence_level = 0.95;
          p95.delta_chi2 = solution.m_cov_scale * 3.841458820694124;
          profile.intervals = {p68, p95};
          profile.status = ProfileStatus::BoundaryLimited;
          profile.message = "The input fixes this mass fraction exactly.";
          if( !probe_only )
            profile_slot( mass_fraction_target(curve_index,target_source) ) = profile;
          continue;
        }

        // A conditional point that lands BELOW the reported optimum is impossible for a restricted
        // minimization; it means the reported baseline was not the optimum, and it is routed to the
        // deferred-rebaseline flow rather than reported as a negative delta chi2.
        std::unique_ptr<RelActAutoSolution> better_conditional;
        const double baseline_improvement_tolerance
            = RelActCalcAutoImp::ProfileLikelihood::baseline_improvement_tolerance(
                                                solution.m_chi2,solution.m_cov_scale);

        const size_t host_solves_before = profile_host ? profile_host->num_conditional_solves() : 0;
        const std::chrono::steady_clock::time_point profile_start
                                                      = std::chrono::steady_clock::now();

        const std::array<double,2> thresholds{{
          solution.m_cov_scale,
          solution.m_cov_scale * 3.841458820694124
        }};

        // --- The carrier engine: pin the reported coordinate itself and re-solve in place ---------
        //
        // `install_carrier_reparam` makes the target the sole range nuclide of its element (through
        // a temporary wide, non-binding [0,1] mass-fraction constraint, or a user constraint that
        // already has that shape), so the target's own slot IS its pre-Pu242-correction mass
        // fraction.  Pinning that slot scans the REPORTED coordinate directly - the textbook exact
        // profile, aligned by construction (a dev-check below verifies the gradient structure).
        // Measured (2026-08 coverage study): pinning an ACTIVITY instead biased the interval narrow
        // by ~|r|, with 95% coverage tracking the squared pin/reported correlation r^2 in lockstep
        // (U-235 r^2 ~ 0.09 covered 0.17 of a nominal 0.95), while this construction
        // restored coverage to the exact reference at the same per-point cost - so an inexact scan
        // is refused rather than run.  No residual is added and the layout is untouched, so a
        // conditional point stays one warm in-place `ceres::Solve`.
        RelActCalcAutoImp::ProfileConditionalHost::CarrierReparam carrier;
        std::string pin_refusal;
        if( profile_host )
        {
          carrier = profile_host->install_carrier_reparam( target_source, curve_index );
          pin_refusal = carrier.why_not;
        }
        const std::optional<size_t> pin_index = carrier.index;

        bool used_pinned_engine = false;

        if( pin_index )
        {
          // One manifold for the whole of this target's scan; the restore also unwinds any installed
          // carrier reparameterization, on every exit path.
          struct PinGuard
          {
            RelActCalcAutoImp::ProfileConditionalHost *host;
            ~PinGuard(){ if( host ) host->restore_optimum(); }
          } pin_guard{ profile_host.get() };

          used_pinned_engine = profile_host->set_pin( *pin_index );
          if( !used_pinned_engine )
            pin_refusal = "The profile manifold could not be installed.";

          if( used_pinned_engine )
          {
            const std::vector<double> &optimum = profile_host->optimum();
            const double baseline_parameter = profile_host->optimum_parameter( *pin_index );
            const std::pair<double,double> parameter_box
                = profile_host->parameter_bounds( *pin_index );

            // The slope of the reported quantity with respect to the pinned coordinate is the
            // chart's own scale, exactly: `q_pre = sig_lo + (x - offset)*chart_scale` (for Pu the
            // correlation renormalizes on top, a small correction irrelevant to probe placement).
            const double reported_per_parameter = carrier.chart_scale;

            const bool have_covariance
                = (solution.m_covariance.size() == optimum.size())
                  && (*pin_index < optimum.size());

  #if( PERFORM_DEVELOPER_CHECKS )
            // The derivative-structure counterpart of the install's nominal-invariance check: on the
            // carrier chart the reported (pre-correction) fraction is a function of the pinned slot
            // ALONE, so its gradient must have support only there.  Off-slot support would mean some
            // eval-path consumer routes the fraction through another parameter - the scan would then
            // silently be of the wrong coordinate.  Pu is exempt only because the Pu-242 correlation
            // legitimately couples the REPORTED coordinate to the siblings downstream of the carrier.
            if( !has_pu_reporting_transform )
            {
              try
              {
                const std::vector<double> gradient
                    = profile_host->reported_gradient( target_source, curve_index, optimum );
                double off_slot = 0.0;
                for( size_t par = 0; par < gradient.size(); ++par )
                  if( par != *pin_index )
                    off_slot = (std::max)( off_slot, std::fabs(gradient[par]) );
                assert( (*pin_index < gradient.size())
                        && (std::fabs(gradient[*pin_index]) > 0.0)
                        && (off_slot <= 1.0e-9*std::fabs(gradient[*pin_index])) );
              }catch( const std::exception & )
              {
                //an unevaluable gradient is not a structure violation
              }
            }
  #endif

            // Scan limits: a carrier slot always has a genuine finite box (the chart's [offset,
            // 1+offset], possibly narrowed by the user's window), so the box edges ARE the feasible
            // limits - unlike the open activity boxes the pre-carrier engine had to cap.  The honest
            // bound detection is still the read-back in the evaluator, never these edges.
            const auto has_real_bound = []( const double edge ){
              return std::isfinite(edge)
                     && (std::fabs(edge) < 0.5*std::numeric_limits<double>::max());
            };
            assert( has_real_bound(parameter_box.first) && has_real_bound(parameter_box.second) );

            // Everything kind-specific is decided; hand the rest to the shared runner.
            ScanSetup setup;
            setup.target = mass_fraction_target( curve_index, target_source );
            setup.reason = profile.reason;
            setup.slot = *pin_index;
            setup.chart_scale = reported_per_parameter;
            setup.baseline_parameter = baseline_parameter;
            setup.nominal_reported = nominal.fraction;
            setup.parameter_box = parameter_box;
            setup.scan_lower = (std::min)( parameter_box.first, baseline_parameter );
            setup.scan_upper = (std::max)( parameter_box.second, baseline_parameter );
            setup.domain_lower = 0.0;      //a mass fraction outside [0,1] is a failed solve
            setup.domain_upper = 1.0;
            setup.feasible_lower = control_lower;   //the projected control window, whose edges
            setup.feasible_upper = control_upper;   // `lower/upper_bound_kind` describe
            setup.feasible_tolerance
                = (std::max)( 1.0e-12, 1.0e-9*(control_upper - control_lower) );
            setup.lower_bound_kind = lower_bound_kind;
            setup.upper_bound_kind = upper_bound_kind;
            setup.thresholds = thresholds;
            setup.use_level_set_refine = have_covariance;
            setup.curve_index = curve_index;
            setup.semantic_key = semantic_curve_key(solution.m_options,curve_index)
                                 + "|target=" + semantic_source_key(target_source)
                                 + "|kind=" + Options::ProfileTarget::to_str(
                                                  Options::ProfileTarget::Kind::MassFraction);
            // The local one-sigma, carried into the pinned coordinate so the probe starts
            // somewhere sensible.  Absent for a weak quantity, which is the common case; the
            // placement is self-calibrating and falls back to a fraction of the span.
            double parameter_one_sigma = 0.0;
            if( nominal.covariance_one_sigma && (reported_per_parameter > 0.0) )
              parameter_one_sigma
                  = std::fabs( *nominal.covariance_one_sigma / reported_per_parameter );

            setup.stats_name = target->symbol;
            setup.dq_key = "dq_da";
            setup.cancel_message = "Mass-fraction profiling was canceled.";
            setup.domain_message = "A conditional solve produced a mass fraction outside [0,1].";
            setup.failed_warning_suffix = " mass-fraction profile failed: ";
            setup.plain_warning_suffix = " mass-fraction profile: ";
            setup.parameter_one_sigma = parameter_one_sigma;

            if( run_one_target( setup, probe_only ) == TargetOutcome::AbortPass )
              return;
            continue;
          }//if( used_pinned_engine )
        }//if( pin_index )

        // Probes have no per-target result: every write - including the pin-refusal Failed just
        // below - belongs to the scan pass.
        if( probe_only )
          continue;

        // The carrier scan is the only conditional mechanism there is.  When no slot scans the
        // reported fraction exactly - a mass fraction fixed by input, an element whose constraints
        // leave the target sharing its block, a ratio-involved or activity-bounded target, a slot
        // held constant by the selected model, or no retained problem to pin within - there is
        // deliberately nothing to fall back to, and saying so plainly is the whole point: a
        // structured `Failed` carrying the reason is honest, where an inexact scan would quote a
        // confidently narrow interval.
        assert( !used_pinned_engine );
        profile.status = ProfileStatus::Failed;
        profile.message = pin_refusal.empty()
            ? std::string("No retained optimizer problem was available to profile against.")
            : pin_refusal;
        profile_slot( mass_fraction_target(curve_index,target_source) ) = profile;
        solution.m_warnings.push_back( target->symbol + " mass-fraction profile failed: "
                                       + profile.message );
      }
    }

    // --- Explicitly requested non-mass-fraction targets ---------------------------------------
    //
    // `RelativeActivity` and `Age` scan an IDENTITY chart: the reported quantity is exactly
    // linear in one existing slot (`activity_identity_chart` / `age_identity_chart`), so there
    // is nothing to install or unwind and the pinned statistic is the exact profile.
    // `ActivityRatio` installs the slot-driven ratio reparameterization
    // (`install_ratio_reparam`): pinning a bare numerator activity would NOT fix the ratio (the
    // denominator re-optimizes - the same misalignment the carrier eliminated for fractions), so
    // the numerator's slot is temporarily reinterpreted as the ratio coordinate instead.  All
    // three reuse the same pin/solve/scan machinery, probe pass, and deferred-rebaseline flow as
    // the mass-fraction targets above.  Automatic selection never chooses these kinds: every
    // entry is an explicit request, honored whatever the budget.
    //
    // No `refine_on_level_set` here: on these charts the pin IS the reported coordinate (the
    // squared pin/reported correlation is 1 by construction), so the level-set correction the
    // mass-fraction path applies has nothing systematic to recover.
    // A forced per-nuclide request on a sole-isotope element is honored through the quantity
    // that carries information there - its RELATIVE ACTIVITY (see the placeholder branch above).
    // Synthesize those requests here so both passes see the identical deterministic list.
    std::vector<Options::ProfileTarget> kind_requests = solution.m_options.profile_targets;
    for( size_t curve_index = 0; curve_index < num_curves; ++curve_index )
    {
      for( const NucInputInfo &input : solution.m_options.rel_eff_curves[curve_index].nuclides )
      {
        if( !input.force_profile_mass_fraction )
          continue;
        const std::vector<Options::ProfileTarget::Kind> kinds
            = RelActCalcAuto::profilable_quantity_kinds( solution.m_options, input.source,
                                                         curve_index );
        const bool mass_fraction_meaningful
            = std::find( begin(kinds), end(kinds), Options::ProfileTarget::Kind::MassFraction )
              != end(kinds);
        const bool activity_offered
            = std::find( begin(kinds), end(kinds),
                         Options::ProfileTarget::Kind::RelativeActivity ) != end(kinds);
        // Substitution is ONLY for sole-isotope elements (fraction identically one); a
        // constraint-blocked mass fraction keeps its structured Failed instead.
        const SandiaDecay::Nuclide * const input_nuc = RelActCalcAuto::nuclide(input.source);
        size_t same_element_count = 0;
        for( const NucInputInfo &other
                  : solution.m_options.rel_eff_curves[curve_index].nuclides )
        {
          const SandiaDecay::Nuclide * const other_nuc = RelActCalcAuto::nuclide(other.source);
          same_element_count += (input_nuc && other_nuc
                                 && (other_nuc->atomicNumber == input_nuc->atomicNumber));
        }
        if( mass_fraction_meaningful || !activity_offered || (same_element_count >= 2) )
          continue;
        Options::ProfileTarget request;
        request.kind = Options::ProfileTarget::Kind::RelativeActivity;
        request.source = input.source;
        request.rel_eff_curve_index = curve_index;
        if( std::find( begin(kind_requests), end(kind_requests), request )
            == end(kind_requests) )
          kind_requests.push_back( request );
      }
    }//for( synthesize sole-isotope RelativeActivity requests )

    for( const Options::ProfileTarget &request : kind_requests )
    {
      using Kind = Options::ProfileTarget::Kind;
      using EvalStatus = RelActCalcAutoImp::ProfileLikelihood::EvaluationStatus;
      using PLDirection = RelActCalcAutoImp::ProfileLikelihood::Direction;
      if( request.kind == Kind::MassFraction )
        continue;   //folded into the per-nuclide mass-fraction loop above

      Profile profile;
      profile.reason = ProfileReason::Forced;
      const std::string request_name = request.display_name();

      const auto fail_request = [&]( const std::string &message ){
        if( probe_only )
          return;
        profile.status = ProfileStatus::Failed;
        profile.message = message;
        profile_slot( request ) = profile;
        solution.m_warnings.push_back( request_name + " profile failed: " + message );
      };

      if( !profile_host )
      {
        fail_request( "No retained optimizer problem was available to profile against." );
        continue;
      }

      RelActCalcAutoImp::ProfileConditionalHost::CarrierReparam chart;
      switch( request.kind )
      {
        case Kind::MassFraction:
          assert( 0 );   //handled by the per-nuclide loop
          break;
        case Kind::RelativeActivity:
          chart = profile_host->activity_identity_chart( request.source,
                                                         request.rel_eff_curve_index );
          break;
        case Kind::ActivityRatio:
          chart = profile_host->install_ratio_reparam( request );
          break;
        case Kind::Age:
          chart = profile_host->age_identity_chart( request.source, request.rel_eff_curve_index );
          break;
      }//switch( request.kind )

      // Unwinds any installed ratio reparameterization (and the pin below) on every exit path.
      struct PinGuard
      {
        RelActCalcAutoImp::ProfileConditionalHost *host;
        ~PinGuard(){ if( host ) host->restore_optimum(); }
      } pin_guard{ profile_host.get() };

      if( !chart.index )
      {
        fail_request( chart.why_not.empty() ? std::string("No slot scans this quantity exactly.")
                                            : chart.why_not );
        continue;
      }
      const size_t slot = *chart.index;

      if( !profile_host->set_pin( slot ) )
      {
        fail_request( "The profile manifold could not be installed." );
        continue;
      }

      const std::vector<double> &optimum = profile_host->optimum();
      const double baseline_parameter = profile_host->optimum_parameter( slot );
      const std::pair<double,double> parameter_box = profile_host->parameter_bounds( slot );

      double nominal_reported = std::numeric_limits<double>::quiet_NaN();
      try
      {
        nominal_reported = profile_host->reported_quantity( request, optimum );
      }catch( const std::exception &e )
      {
        fail_request( std::string("The nominal value could not be evaluated: ") + e.what() );
        continue;
      }
      if( !std::isfinite(nominal_reported) || !std::isfinite(solution.m_chi2) )
      {
        fail_request( "The nominal value or objective was non-finite." );
        continue;
      }

      std::unique_ptr<RelActAutoSolution> better_conditional;
      const double baseline_improvement_tolerance
          = RelActCalcAutoImp::ProfileLikelihood::baseline_improvement_tolerance(
                                              solution.m_chi2,solution.m_cov_scale);
      const size_t host_solves_before = profile_host->num_conditional_solves();
      const std::chrono::steady_clock::time_point profile_start
                                                    = std::chrono::steady_clock::now();
      const std::array<double,2> thresholds{{
        solution.m_cov_scale,
        solution.m_cov_scale * 3.841458820694124
      }};

      const auto has_real_bound = []( const double edge ){
        return std::isfinite(edge)
               && (std::fabs(edge) < 0.5*std::numeric_limits<double>::max());
      };

      // Scan limits.  A real box edge is a feasible limit (age boxes are finite on both sides);
      // an open activity/ratio box gets the pre-carrier engine's deliberately MODEST synthetic
      // cap - ten times the parameter's own natural scale, far outside any plausible crossing
      // while still somewhere the optimizer can work.  The honest bound detection is the
      // read-back in the evaluator, never this cap.  (Ceres writes "no bound" as
      // +-numeric_limits<double>::max(), which isfinite() accepts - hence `has_real_bound`.)
      const double parameter_scale = (std::max)( 1.0e-9,
          has_real_bound(parameter_box.first) ? (baseline_parameter - parameter_box.first)
                                              : std::fabs(baseline_parameter) );
      const double lower_parameter = has_real_bound(parameter_box.first)
          ? (std::min)( parameter_box.first, baseline_parameter )
          : (baseline_parameter - 10.0*parameter_scale);
      const double upper_parameter = has_real_bound(parameter_box.second)
          ? (std::max)( parameter_box.second, baseline_parameter )
          : (baseline_parameter + 10.0*parameter_scale);

      // What a non-crossing endpoint on each side actually MEANS, which is three different
      // statements the report must not conflate:
      //  - no real box edge: the side stopped at the synthetic cap above, which is a scan-range
      //    choice and NOT a feasibility bound (activities and ratios have no upper limit);
      //  - a real edge at the coordinate's zero: the physical floor (activity, ratio, or age
      //    cannot be negative);
      //  - any other real edge: it came from the input (min/max_rel_act, a fit_age window, or the
      //    model's age cap), so it is an input constraint.
      const double zero_coordinate
          = (request.kind == Kind::Age) ? 0.0
                        : RelActCalcAutoImp::RelActAutoCostFcn::sm_activity_par_offset;
      const auto endpoint_kind_for = [&]( const double edge, const bool real_bound ){
        if( !real_bound )
          return EndpointKind::ScanRangeLimit;
        return (std::fabs(edge - zero_coordinate) <= 1.0e-9*(std::max)(1.0,std::fabs(edge)))
               ? EndpointKind::PhysicalLimit : EndpointKind::InputConstraintLimit;
      };
      const EndpointKind lower_bound_kind
          = endpoint_kind_for( parameter_box.first, has_real_bound(parameter_box.first) );
      const EndpointKind upper_bound_kind
          = endpoint_kind_for( parameter_box.second, has_real_bound(parameter_box.second) );

      // Every one of these reported quantities is bounded below by zero and (activities, ratios)
      // open above.
      const double reported_lower_limit = 0.0;
      const double reported_tolerance
          = (std::max)( 1.0e-12, 1.0e-9*(std::max)( 1.0, std::fabs(nominal_reported) ) );

      // The local one sigma IN THE SCANNED COORDINATE, which is what places the opening probe.
      //
      // For an identity chart the pinned slot IS a production parameter, so its own covariance
      // diagonal is exactly that.  For the RATIO chart it is NOT: `install_ratio_reparam` has
      // re-encoded the slot (nominal at `offset + 1`, `dq/dx == chart_scale == r0`), so the
      // production diagonal describes a coordinate the scan is no longer moving - using it
      // mis-scales the probe by the fitted activity's distance from the offset and ignores the
      // denominator's contribution entirely.  Take the reported ratio's own uncertainty and carry
      // it into the chart with the chart scale; an unavailable covariance leaves 0, which the
      // scanner reads as "no usable sigma" and falls back to a fraction of the span.
      double parameter_one_sigma = 0.0;
      if( request.kind == Kind::ActivityRatio )
      {
        try
        {
          const double ratio_sigma = solution.activity_ratio_uncertainty(
                          request.source, request.rel_eff_curve_index,
                          request.denominator, request.denominator_curve_index );
          if( std::isfinite(ratio_sigma) && (ratio_sigma > 0.0) && (chart.chart_scale > 0.0) )
            parameter_one_sigma = ratio_sigma/chart.chart_scale;
        }catch( const std::exception & )
        {
          //no covariance for the ratio; the span fallback places the probe
        }
      }else if( (solution.m_covariance.size() == optimum.size())
                && (slot < solution.m_covariance.size())
                && (slot < solution.m_covariance[slot].size())
                && (solution.m_covariance[slot][slot] > 0.0) )
      {
        parameter_one_sigma = std::sqrt( solution.m_covariance[slot][slot] );
      }

      ScanSetup setup;
      setup.target = request;
      setup.reason = ProfileReason::Forced;   //every kind target is an explicit request
      setup.slot = slot;
      setup.chart_scale = chart.chart_scale;
      setup.baseline_parameter = baseline_parameter;
      setup.nominal_reported = nominal_reported;
      setup.parameter_box = parameter_box;
      setup.scan_lower = lower_parameter;
      setup.scan_upper = upper_parameter;
      // Activities, ratios and ages are bounded below by zero and open above, so the rejection
      // window and the feasible window coincide here (unlike a mass fraction).
      setup.domain_lower = reported_lower_limit;
      setup.domain_upper = std::numeric_limits<double>::infinity();
      setup.feasible_lower = reported_lower_limit;
      setup.feasible_upper = std::numeric_limits<double>::infinity();
      setup.feasible_tolerance = reported_tolerance;
      setup.lower_bound_kind = lower_bound_kind;
      setup.upper_bound_kind = upper_bound_kind;
      setup.thresholds = thresholds;
      setup.use_level_set_refine = false;   //identity/ratio charts pin the reported coordinate
      setup.curve_index = request.rel_eff_curve_index;
      setup.semantic_key = semantic_curve_key(solution.m_options,request.rel_eff_curve_index)
                           + "|target=" + semantic_source_key(request.source)
                           + "|kind=" + Options::ProfileTarget::to_str(request.kind);
      setup.stats_name = request_name;
      setup.dq_key = "dq_dx";
      setup.cancel_message = "Profiling was canceled.";
      setup.domain_message = "A conditional solve produced a value outside its domain.";
      setup.failed_warning_suffix = " profile failed: ";
      setup.plain_warning_suffix = " profile: ";
      setup.parameter_one_sigma = parameter_one_sigma;
      setup.emit_interval_stats = true;

      if( run_one_target( setup, probe_only ) == TargetOutcome::AbortPass )
        return;
    }//for( non-mass-fraction profile_targets )
  };//run_target_pass( probe_only )

  // Consume the best deferred discovery: restart the solve on the frozen problem from the
  // discovered basin's exact seed, and on success recurse so covariance, thresholds, and every
  // profile use the new baseline.  `from_probe` selects the cheap restart shape (no candidate
  // matrix; the one full polish is deferred to the end of the probe-hop chain); a scan-phase
  // discovery keeps the historical full-matrix restart.  Returns true when the caller must
  // return (the reselection recursed, was canceled, or - for a scan-phase discovery - failed
  // terminally); false only for a failed probe reselection, where falling back to ordinary
  // scans against the current baseline is safe and loses nothing.
  const auto attempt_reselection = [&]( const bool from_probe ) -> bool
  {
    std::vector<RelActCalcAutoImp::ProfileLikelihood::PendingBaselineDiscovery> ranks;
    ranks.reserve(deferred_baseline_discoveries.size());
    for( const DeferredBaselineDiscovery &discovery : deferred_baseline_discoveries )
      ranks.push_back(discovery.rank);
    const std::optional<size_t> best_index
        = RelActCalcAutoImp::ProfileLikelihood::best_pending_baseline_discovery_index(ranks);
    assert( best_index && (*best_index < deferred_baseline_discoveries.size()) );

    // A PROBE-phase reselection that fails costs the caller nothing: it clears the discoveries and
    // falls back to the ordinary scan pass, which re-discovers and reports on its own terms.
    // Writing failures here would leave permanent warnings claiming a profile failed beside the
    // interval that same profile then produced successfully - so probe-phase failures are silent,
    // and only the scan phase (whose budget really is spent) reports.
    const auto fail_deferred_profiles = [&]( const std::string &message ) {
      if( from_probe )
        return;
      for( const DeferredBaselineDiscovery &discovery : deferred_baseline_discoveries )
      {
        Profile &profile = profile_slot( discovery.target );
        profile.status = ProfileStatus::Failed;
        profile.reason = discovery.reason;
        profile.message = message;
        solution.m_warnings.push_back(
            discovery.target.display_name() + " profile: " + message );
      }
    };

    if( !best_index )
    {
      fail_deferred_profiles( "No finite profile-discovered point was available for baseline"
                              " reselection." );
      return !from_probe;
    }

    DeferredBaselineDiscovery &selected = deferred_baseline_discoveries[*best_index];
    assert( selected.conditional_solution );
    if( !selected.conditional_solution )
    {
      fail_deferred_profiles( "The selected profile-discovered point could not seed baseline"
                              " reselection." );
      return !from_probe;
    }

    // A scan-phase discovery is a deterministic-search trigger: polish the complete named
    // candidate matrix on the same frozen ROIs.  A probe-phase discovery instead re-minimizes
    // matrix-free from the exact seed - hops are cheap, and the one full polish is deferred to the
    // end of the probe-hop chain (see `needs_matrix_polish`).  Either way every requested profile
    // restarts so covariance and thresholds use the selected baseline.
    Options restart_options = solution.m_options;
    for( RoiRange &roi : restart_options.rois )
      roi.range_limits_type = RoiRange::RangeLimitsType::Fixed;
    const auto run_restart = [&]( const bool with_matrix ) -> RelActAutoSolution {
      return RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres(
          restart_options,foreground,background,frozen_profile_drf,frozen_profile_peaks,det_type,
          cancel_calc,RelActCalcAutoImp::SearchSeedVariant::Default,
          /*force_candidate_search=*/with_matrix,selected.conditional_solution.get(),
          /*allow_candidate_search=*/with_matrix,&frozen_profile_peak_ranges,
          &frozen_profile_model_policy,true,/*may_host_profile=*/true );
    };
    RelActAutoSolution restarted;
    bool used_matrix = !from_probe;
    try
    {
      restarted = run_restart( !from_probe );

      // MEASURED escalation (2026-08-31, Pu 4h arm): a matrix-free hop occasionally stalls short
      // of its seed - 6 of 71 hops - because the semantic warm start's re-expression carries a
      // small loss the candidate matrix historically absorbed (fidelity was zero with the matrix,
      // never without).  An exact parameter-vector copy was tried instead and REVERTED: the
      // QR/gauge frame shifts between functor instances, so an identical vector does not
      // evaluate to the seed's objective (see the note in `solve_ceres`'s warm-start block).
      // A cheap hop that misses its seed therefore re-runs WITH the matrix rather than shipping
      // the stall; the fidelity warning below then only marks post-matrix misses.
      const double seed_gate = RelActCalcAutoImp::ProfileLikelihood::baseline_improvement_tolerance(
                                   selected.rank.full_objective, solution.m_cov_scale );
      if( from_probe
          && (restarted.m_status != RelActAutoSolution::Status::UserCanceled)
          && (!cancel_calc || !cancel_calc->load())
          && (!RelActAutoSolution::is_usable_status(restarted.m_status)
              || !std::isfinite(restarted.m_chi2)
              || (restarted.m_chi2 > (selected.rank.full_objective + seed_gate))) )
      {
        if( RelActCalcAutoImp::profile_stats_enabled() )
          std::cerr << "profile-rebaseline-escalate hop=" << (baseline_restart_count + 1)
                    << " cheap=" << restarted.m_chi2
                    << " seed=" << selected.rank.full_objective << std::endl;
        restarted = run_restart( true );
        used_matrix = true;
      }
    }catch( const std::exception &e )
    {
      fail_deferred_profiles(
          std::string("Baseline reselection after profile-discovered points failed: ") + e.what() );
      return !from_probe;
    }

    if( restarted.m_status == RelActAutoSolution::Status::UserCanceled
        || (cancel_calc && cancel_calc->load()) )
    {
      solution.m_status = RelActAutoSolution::Status::UserCanceled;
      solution.m_error_message = "Mass-fraction baseline reselection was canceled.";
      return true;
    }

    if( RelActAutoSolution::is_usable_status(restarted.m_status) )
    {
      const bool same_frozen_objective = matches_frozen_profile_objective(restarted);
      assert( same_frozen_objective );
      if( !same_frozen_objective )
      {
        fail_deferred_profiles( "Baseline reselection changed the frozen branching-ratio,"
                                " gamma-membership, or selected-model policy." );
        return !from_probe;
      }
    }

    // The same scale that decided this reselection was worth attempting decides whether it succeeded;
    // a restart which only shaved round-off off the objective has not established a better baseline.
    const double baseline_improvement_tolerance
        = RelActCalcAutoImp::ProfileLikelihood::baseline_improvement_tolerance(
                                            solution.m_chi2,solution.m_cov_scale);
    if( !RelActAutoSolution::is_usable_status(restarted.m_status)
        || !std::isfinite(restarted.m_chi2)
        || !(restarted.m_chi2 < solution.m_chi2-baseline_improvement_tolerance) )
    {
      if( RelActCalcAutoImp::profile_stats_enabled() )
        std::cerr << "profile-rebaseline-rejected hop=" << (baseline_restart_count + 1)
                  << " from=" << solution.m_chi2 << " restarted=" << restarted.m_chi2
                  << " seed=" << selected.rank.full_objective << std::endl;
      fail_deferred_profiles( "Profiles found better physical points, but final candidate selection"
                              " did not establish a better unconstrained baseline." );
      return !from_probe;
    }

    const double old_objective = solution.m_chi2;
    const double seed_objective = selected.rank.full_objective;
    const size_t num_discoveries = deferred_baseline_discoveries.size();
    const std::string selected_target = selected.target_symbol;
    std::vector<size_t> discovery_order;
    discovery_order.reserve(num_discoveries);
    for( size_t index = 0; index < num_discoveries; ++index )
      discovery_order.push_back(index);
    std::sort( begin(discovery_order),end(discovery_order),
        [&deferred_baseline_discoveries]( const size_t lhs, const size_t rhs ) {
          return deferred_baseline_discoveries[lhs].rank.semantic_key
               < deferred_baseline_discoveries[rhs].rank.semantic_key;
        } );
    std::string discovery_summary;
    for( const size_t index : discovery_order )
    {
      const DeferredBaselineDiscovery &discovery = deferred_baseline_discoveries[index];
      if( !discovery_summary.empty() )
        discovery_summary += ", ";
      discovery_summary += discovery.target_symbol + "="
          + SpecUtils::printCompact(discovery.rank.full_objective,12);
    }

    adopt_restarted_solution( std::move(restarted) );
    const std::string discovery_phrase = from_probe
        ? std::string(" pre-scan probe(s) found better points; the discovered basin's exact vector"
                      " re-seeded a matrix-free re-minimization, from the best frozen-objective"
                      " seed (")
        : std::string(" mass-fraction profile(s) found better points; deterministic final candidate"
                      " selection used the best frozen-objective seed (");
    push_reselection_warning( std::to_string(baseline_restart_count + 1) + " of "
        + std::to_string(RelActCalcAutoImp::ProfileLikelihood::sm_max_baseline_restarts) + ": "
        + std::to_string(num_discoveries)
        + discovery_phrase + selected_target + "; discoveries: "
        + discovery_summary + ") and lowered the objective from "
        + std::to_string(old_objective) + " to " + std::to_string(solution.m_chi2)
        + ", so covariance and all profiles were restarted against the new baseline." );
    if( RelActCalcAutoImp::profile_stats_enabled() )
      std::cerr << "profile-rebaseline hop=" << (baseline_restart_count + 1)
                << " of=" << RelActCalcAutoImp::ProfileLikelihood::sm_max_baseline_restarts
                << " from=" << old_objective << " to=" << solution.m_chi2
                << " seed=" << seed_objective << " target=" << selected_target
                << " discoveries=" << num_discoveries
                << " matrix=" << used_matrix << std::endl;

    // The seed is a physically realizable point of the SAME unconstrained problem, so a re-solve
    // warm-started from it finishing ABOVE it is proof the semantic warm start did not reach the
    // discovered basin (or the solve terminated on a plateau).  The improved baseline is still kept
    // - it beat the old one - but say so loudly: this is the signature that separates genuine
    // basin chains from a restart crawling toward one basin it keeps rediscovering.
    if( std::isfinite(seed_objective)
        && (solution.m_chi2 > seed_objective + baseline_improvement_tolerance) )
    {
      push_reselection_warning( std::to_string(baseline_restart_count + 1)
          + " did not reach its seed's known objective ("
          + std::to_string(solution.m_chi2) + " vs " + std::to_string(seed_objective)
          + "); the warm start likely failed to enter the discovered basin, which may be"
            " rediscovered on the next profile pass." );
      if( RelActCalcAutoImp::profile_stats_enabled() )
        std::cerr << "profile-rebaseline-fidelity hop=" << (baseline_restart_count + 1)
                  << " final=" << solution.m_chi2 << " seed=" << seed_objective << std::endl;
    }
    add_merged_single_curve_comparison( solution,restart_options,foreground,background,
                                        solution.m_drf,solution.m_spectrum_peaks,
                                        det_type,cancel_calc );
    add_mass_fraction_profiles( solution,foreground,background,solution.m_drf,
                                solution.m_spectrum_peaks,
                                det_type,cancel_calc,baseline_restart_count+1,
                                /*needs_matrix_polish=*/from_probe );

    return true;
  };//attempt_reselection( from_probe )

  // --- Pass orchestration -------------------------------------------------------------------
  //
  // Probe first.  MEASURED (pu4h_multirestart, 242 discovering scans): ~93% of better-baseline
  // discoveries happen at a side's opening probe or first ladder rung, so two conditional solves
  // per target vet the baseline at a small fraction of a scan pass's cost, and a discovered basin
  // re-seeds the solve BEFORE any scan spends its budget against the doomed baseline.  The scan
  // pass keeps its own deferred-rebaseline flow as the backstop for the far-displacement tail
  // (measured: Pu-240 upper-side discoveries at 3-16x the opening step, and degenerate baselines
  // whose collapsed sigma leaves probes near-nominal).
  const bool may_defer_discoveries
      = (RelActCalcAutoImp::ProfileLikelihood::baseline_discovery_disposition(baseline_restart_count)
           == RelActCalcAutoImp::ProfileLikelihood::BaselineDiscoveryDisposition::DeferUntilPassComplete);
  if( may_defer_discoveries )
  {
    run_target_pass( /*probe_only=*/true );
    if( solution.m_status == RelActAutoSolution::Status::UserCanceled )
      return;
    if( !deferred_baseline_discoveries.empty() )
    {
      if( attempt_reselection( /*from_probe=*/true ) )
        return;
      deferred_baseline_discoveries.clear();
    }
  }//if( may_defer_discoveries )

  // A probe-hop restart deliberately skipped the candidate-matrix polish; the chain has now
  // settled (this level's probes found nothing better), so apply the original design intent - "a
  // profile-discovered point is a deterministic-search trigger" - ONCE, to the final baseline,
  // before any scan trusts it.  An accepted polish is a hop like any other: everything restarts
  // against it.  A rejected polish CONFIRMS the probe-selected baseline is matrix-stable, which
  // is the insurance this one extra solve buys.
  // Deliberately NOT gated on `may_defer_discoveries`: the polish is an obligation carried from
  // earlier cheap hops, not a new discovery, so a chain that exhausted its reselection budget is
  // exactly the case that must not ship intervals measured against a never-matrix-polished
  // baseline.  (An accepted polish recurses with the obligation cleared, so it runs at most once.)
  if( needs_matrix_polish )
  {
    Options polish_options = solution.m_options;
    for( RoiRange &roi : polish_options.rois )
      roi.range_limits_type = RoiRange::RangeLimitsType::Fixed;
    RelActAutoSolution polished;
    bool polish_usable = false;
    try
    {
      polished = RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres(
          polish_options,foreground,background,frozen_profile_drf,frozen_profile_peaks,det_type,
          cancel_calc,RelActCalcAutoImp::SearchSeedVariant::Default,
          /*force_candidate_search=*/true,&solution,/*allow_candidate_search=*/true,
          &frozen_profile_peak_ranges,&frozen_profile_model_policy,true,
          /*may_host_profile=*/true );
      polish_usable = RelActAutoSolution::is_usable_status( polished.m_status );
    }catch( const std::exception & )
    {
      polish_usable = false;  //the un-polished baseline remains vetted and trusted
    }

    if( (polished.m_status == RelActAutoSolution::Status::UserCanceled)
        || (cancel_calc && cancel_calc->load()) )
    {
      solution.m_status = RelActAutoSolution::Status::UserCanceled;
      solution.m_error_message = "Mass-fraction baseline reselection was canceled.";
      return;
    }

    const double polish_tolerance
        = RelActCalcAutoImp::ProfileLikelihood::baseline_improvement_tolerance(
                                            solution.m_chi2,solution.m_cov_scale);
    const bool same_frozen = !polish_usable || matches_frozen_profile_objective(polished);
    assert( same_frozen );
    const bool polish_accepted = polish_usable && same_frozen
        && std::isfinite(polished.m_chi2)
        && (polished.m_chi2 < solution.m_chi2 - polish_tolerance);
    if( RelActCalcAutoImp::profile_stats_enabled() )
      std::cerr << "profile-rebaseline-polish hop=" << (baseline_restart_count + 1)
                << " from=" << solution.m_chi2 << " to=" << polished.m_chi2
                << " accepted=" << polish_accepted << std::endl;

    if( polish_accepted )
    {
      const double old_objective = solution.m_chi2;
      adopt_restarted_solution( std::move(polished) );
      push_reselection_warning( std::to_string(baseline_restart_count + 1) + " of "
          + std::to_string(RelActCalcAutoImp::ProfileLikelihood::sm_max_baseline_restarts)
          + ": the final candidate-matrix polish of the probe-selected baseline lowered the"
            " objective from " + std::to_string(old_objective) + " to "
          + std::to_string(solution.m_chi2)
          + ", so covariance and all profiles were restarted against the polished baseline." );
      add_merged_single_curve_comparison( solution,polish_options,foreground,background,
                                          solution.m_drf,solution.m_spectrum_peaks,
                                          det_type,cancel_calc );
      add_mass_fraction_profiles( solution,foreground,background,solution.m_drf,
                                  solution.m_spectrum_peaks,
                                  det_type,cancel_calc,baseline_restart_count+1,
                                  /*needs_matrix_polish=*/false );
      return;
    }//if( polish_accepted )
  }//if( needs_matrix_polish )

  run_target_pass( /*probe_only=*/false );
  if( solution.m_status == RelActAutoSolution::Status::UserCanceled )
    return;
  if( !deferred_baseline_discoveries.empty() )
    attempt_reselection( /*from_probe=*/false );

}
}//anonymous namespace

#endif //InterSpec_RelActCalcAuto_Profile_imp_hpp
