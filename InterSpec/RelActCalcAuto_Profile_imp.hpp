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
 */
void add_mass_fraction_profiles( RelActAutoSolution &solution,
                                 std::shared_ptr<const SpecUtils::Measurement> foreground,
                                 std::shared_ptr<const SpecUtils::Measurement> background,
                                 std::shared_ptr<const DetectorPeakResponse> input_drf,
                                 const std::vector<std::shared_ptr<const PeakDef>> &all_peaks,
                                 const PeakFitUtils::CoarseResolutionType det_type,
                                 const std::shared_ptr<std::atomic_bool> &cancel_calc,
                                 const unsigned baseline_restart_count = 0 )
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
  solution.m_mass_fraction_profiles.clear();
  solution.m_mass_fraction_profiles.resize( num_curves );

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
  bool any_forced_profile = false;
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
    size_t curve_index = 0;
    std::string target_symbol;
  };
  std::vector<DeferredBaselineDiscovery> deferred_baseline_discoveries;
  // Warn only once per pass when the reselection budget was exhausted (see the RejectAfterRestart
  // branch): the taint applies to the whole solution, not per discovering target.
  bool budget_exhaustion_warned = false;

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
      const bool forced = profile_target.input->force_profile_mass_fraction;

      // Which quantities carry information about this source.  Automatic selection and every
      // reporting surface below are mass-fraction typed today, so a source whose only profilable
      // quantity is its activity or its age is skipped here rather than silently reported as a mass
      // fraction; the constraint machinery itself handles all four kinds (see
      // `Options::ProfileTarget::Kind`).
      const std::vector<Options::ProfileTarget::Kind> profilable_kinds
          = RelActCalcAuto::profilable_quantity_kinds( solution.m_options,target_source,curve_index );
      const bool mass_fraction_is_meaningful
          = std::find( begin(profilable_kinds),end(profilable_kinds),
                       Options::ProfileTarget::Kind::MassFraction )
            != end(profilable_kinds);
      if( !mass_fraction_is_meaningful )
      {
        // An explicit request always produces a structured result.  With no second isotope of the
        // element in the model, the normalized fraction is identically one and there is no
        // conditional optimization to perform.
        if( forced )
        {
          const bool activity_would_carry_it
              = std::find( begin(profilable_kinds),end(profilable_kinds),
                           Options::ProfileTarget::Kind::RelativeActivity )
                != end(profilable_kinds);
          Profile profile;
          profile.reason = ProfileReason::Forced;
          profile.status = ProfileStatus::BoundaryLimited;
          profile.num_fits = 0;
          profile.message = "This is the only modeled isotope of its element; its normalized mass"
                            " fraction is fixed at the physical limit.";
          if( activity_would_carry_it )
            profile.message += "  Its relative activity is the quantity that carries information"
                               " here.";
          Interval p68;
          p68.confidence_level = 0.6827;
          p68.delta_chi2 = solution.m_cov_scale;
          p68.lower = p68.upper = 1.0;
          p68.lower_kind = p68.upper_kind = EndpointKind::PhysicalLimit;
          Interval p95 = p68;
          p95.confidence_level = 0.95;
          p95.delta_chi2 = solution.m_cov_scale*3.841458820694124;
          profile.intervals = {p68,p95};
          solution.m_mass_fraction_profiles[curve_index][target->symbol] = std::move(profile);
        }
        continue;
      }

      RelActAutoSolution::MassFractionResult nominal;
      try
      {
        nominal = solution.mass_enrichment_result( target, curve_index );
      }catch( const std::exception &e )
      {
        // Do not silently lose an explicit request when even the structured nominal result cannot
        // be formed: retain the otherwise usable main fit and expose a failed profile.
        if( forced )
        {
          Profile failed;
          failed.status = ProfileStatus::Failed;
          failed.reason = ProfileReason::Forced;
          failed.message = std::string("Could not initialize profile: ") + e.what();
          solution.m_mass_fraction_profiles[curve_index][target->symbol] = failed;
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
        solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
        solution.m_warnings.push_back( target->symbol + " mass-fraction profile failed: " + profile.message );
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
        solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
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

      RelActCalcAutoImp::ProfileLikelihood::ScanResult scan;
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

          if( RelActCalcAutoImp::profile_stats_enabled() )
            std::cerr << "profile-scan " << target->symbol << " curve=" << curve_index
                      << " index=" << *pin_index
                      << " a0=" << baseline_parameter
                      << " box=[" << parameter_box.first << "," << parameter_box.second << "]"
                      << " dq_da=" << reported_per_parameter
                      << " q0=" << nominal.fraction << std::endl;

          const double lower_parameter = (std::min)( parameter_box.first, baseline_parameter );
          const double upper_parameter = (std::max)( parameter_box.second, baseline_parameter );

          // The local one-sigma, carried into the pinned coordinate so the probe starts somewhere
          // sensible.  Absent for a weak quantity, which is the common case; the placement is
          // self-calibrating and falls back to a fraction of the span.
          double parameter_one_sigma = 0.0;
          if( nominal.covariance_one_sigma && (reported_per_parameter > 0.0) )
            parameter_one_sigma
                = std::fabs( *nominal.covariance_one_sigma / reported_per_parameter );

          std::vector<double> lower_warm = optimum, upper_warm = optimum;

          const auto evaluate_pinned = [&]( const double requested_parameter,
                             const RelActCalcAutoImp::ProfileLikelihood::Direction direction,
                             const size_t )
                             -> RelActCalcAutoImp::ProfileLikelihood::Evaluation {
            using Eval = RelActCalcAutoImp::ProfileLikelihood::Evaluation;
            using EvalStatus = RelActCalcAutoImp::ProfileLikelihood::EvaluationStatus;
            using PLDirection = RelActCalcAutoImp::ProfileLikelihood::Direction;

            Eval evaluation;
            if( cancel_calc && cancel_calc->load() )
            {
              evaluation.status = EvalStatus::Canceled;
              evaluation.diagnostic = "Mass-fraction profiling was canceled.";
              return evaluation;
            }

            // Points are visited outward from the optimum, so the previous point of this direction
            // is always the nearest one and is chronologically available - which is what makes the
            // scanner's proximity-based warm-start policy unnecessary here.
            std::vector<double> &warm = (direction == PLDirection::Lower) ? lower_warm : upper_warm;

            const RelActCalcAutoImp::ProfileConditionalHost::PinnedResult point
                = profile_host->solve_pinned( *pin_index, requested_parameter, warm );
            if( RelActCalcAutoImp::profile_stats_enabled() )
              std::cerr << "profile-point " << target->symbol
                        << " dir=" << ((direction == PLDirection::Lower) ? "lo" : "hi")
                        << " req=" << requested_parameter
                        << " got=" << point.achieved_parameter
                        << " conv=" << point.converged
                        << " chi2=" << point.physical_chi2
                        << " " << point.diagnostic << std::endl;

            if( point.canceled )
            {
              evaluation.status = EvalStatus::Canceled;
              evaluation.diagnostic = "Mass-fraction profiling was canceled.";
              return evaluation;
            }
            if( !point.converged )
            {
              evaluation.status = EvalStatus::Failed;
              evaluation.diagnostic = point.diagnostic;
              return evaluation;
            }

            // The pinned point's own reported value, needed both as the sample coordinate and as the
            // level the refinement below holds fixed.
            double pinned_fraction = std::numeric_limits<double>::quiet_NaN();
            try
            {
              pinned_fraction = profile_host->reported_mass_fraction(
                                        target_source, curve_index, point.parameters );
            }catch( const std::exception &e )
            {
              evaluation.status = EvalStatus::Failed;
              evaluation.diagnostic = std::string("A conditional point could not be read back: ")
                                      + e.what();
              return evaluation;
            }

            // Correct the one known bias of pinning, by EVALUATION rather than by subtracting an
            // estimate.  See `refine_on_level_set`: the refined point is still a rigorous upper
            // bound on the exact profile at its own achieved `q`, so this can only move a sample
            // toward the truth from the narrow side and can never widen the interval.
            //
            // Skipped near the vertex: there the gap is a small fraction of an already-small delta
            // chi2, the recoverable amount sits near the objective's own reproducibility floor, and
            // those samples do not move the crossing anyway.
            std::vector<double> sample_parameters = point.parameters;
            double sample_chi2 = point.physical_chi2;
            const double refine_floor = 0.3*thresholds[0];
            if( have_covariance && std::isfinite(pinned_fraction)
                && ((point.physical_chi2 - profile_host->optimum_objective()) >= refine_floor) )
            {
              const RelActCalcAutoImp::ProfileConditionalHost::LevelSetStep refined
                  = profile_host->refine_on_level_set( point.parameters, point.physical_chi2,
                                                       pinned_fraction, nominal.fraction,
                                                       solution.m_covariance,
                                                       target_source, curve_index );
              if( refined.improved )
              {
                sample_parameters = refined.parameters;
                sample_chi2 = refined.physical_chi2;
              }
            }

            double achieved_fraction = std::numeric_limits<double>::quiet_NaN();
            try
            {
              achieved_fraction = profile_host->reported_mass_fraction(
                                        target_source, curve_index, sample_parameters );
            }catch( const std::exception &e )
            {
              evaluation.status = EvalStatus::Failed;
              evaluation.diagnostic = std::string("A conditional point could not be read back: ")
                                      + e.what();
              return evaluation;
            }

            if( !std::isfinite(achieved_fraction) || (achieved_fraction < 0.0)
                || (achieved_fraction > 1.0) )
            {
              evaluation.status = EvalStatus::Failed;
              evaluation.diagnostic = "A conditional solve produced a mass fraction outside [0,1].";
              return evaluation;
            }

            // The refinement holds `q` fixed only to first order, so a refined point lands at its
            // own slightly different reported value.  That is legitimate - the sample is quoted at
            // whatever it achieved - but a LARGE drift could reorder samples and make the hygiene
            // pass see a fold where the map never folded, truncating the side.  Keep the original
            // point rather than risk that.
            if( sample_parameters != point.parameters )
            {
              const double drift = std::fabs( achieved_fraction - pinned_fraction );
              if( drift > 0.1*std::fabs(pinned_fraction - nominal.fraction) )
              {
                sample_parameters = point.parameters;
                sample_chi2 = point.physical_chi2;
                achieved_fraction = pinned_fraction;
              }
            }

            const double delta = sample_chi2 - profile_host->optimum_objective();

            // A conditional is a restricted minimization, so it cannot beat the unconstrained
            // optimum.  One that does means the reported baseline was not the optimum; hand it to
            // the deferred-rebaseline flow rather than reporting a negative delta-chi2.
            if( delta < -baseline_improvement_tolerance )
            {
              try
              {
                better_conditional = std::make_unique<RelActAutoSolution>(
                        profile_host->warm_start_solution(solution,sample_parameters) );
                better_conditional->m_chi2 = sample_chi2;
              }catch( const std::exception & )
              {
                better_conditional.reset();
              }
              evaluation.status = EvalStatus::BetterBaseline;
              evaluation.diagnostic = "A conditional profile point reached a lower physical"
                                      " objective than the reported solution.";
              return evaluation;
            }

            // Deliberately the PINNED point, never the refined one: a refined point is not a
            // conditional minimum of anything, and the outward sweep's warm-start chain - and the
            // exact `x_optimum - x` chord the refinement itself relies on - are only meaningful
            // between genuine conditional minima.
            warm = point.parameters;

            evaluation.status = EvalStatus::Success;
            evaluation.reported_fraction = achieved_fraction;
            evaluation.delta_chi2 = (std::max)( 0.0, delta );

            // The honest bound test, and the primary one: a genuine Ceres box edge, or a reported
            // fraction that has arrived at its feasible limit.  The synthetic scan cap above is
            // deliberately NOT what decides this.
            const bool at_box_edge
                = (has_real_bound(parameter_box.first)
                     && (point.achieved_parameter <= parameter_box.first))
                  || (has_real_bound(parameter_box.second)
                     && (point.achieved_parameter >= parameter_box.second));
            const double fraction_tolerance
                = (std::max)( 1.0e-12, 1.0e-9*(control_upper - control_lower) );
            const bool at_fraction_limit
                = (achieved_fraction <= (control_lower + fraction_tolerance))
                  || (achieved_fraction >= (control_upper - fraction_tolerance));
            evaluation.reached_feasible_bound = at_box_edge || at_fraction_limit;
            return evaluation;
          };//evaluate_pinned

          scan = RelActCalcAutoImp::ProfileLikelihood::fit_profile(
                     baseline_parameter, nominal.fraction, lower_parameter, upper_parameter,
                     thresholds,
                     RelActCalcAutoImp::ProfileLikelihood::sm_profile_max_points_per_quantity,
                     evaluate_pinned, parameter_one_sigma );
        }//if( used_pinned_engine )
      }//if( pin_index )

      // The carrier scan is the only conditional mechanism there is.  When no slot scans the
      // reported fraction exactly - a mass fraction fixed by input, an element whose constraints
      // leave the target sharing its block, a ratio-involved or activity-bounded target, a slot
      // held constant by the selected model, or no retained problem to pin within - there is
      // deliberately nothing to fall back to, and saying so plainly is the whole point: a
      // structured `Failed` carrying the reason is honest, where an inexact scan would quote a
      // confidently narrow interval.
      if( !used_pinned_engine )
      {
        profile.status = ProfileStatus::Failed;
        profile.message = pin_refusal.empty()
            ? std::string("No retained optimizer problem was available to profile against.")
            : pin_refusal;
        solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
        solution.m_warnings.push_back( target->symbol + " mass-fraction profile failed: "
                                       + profile.message );
        continue;
      }
      profile.num_fits = scan.num_fits;
      for( const RelActCalcAutoImp::ProfileLikelihood::Sample &sample : scan.samples )
        profile.sampled_delta_chi2.emplace_back(sample.reported_fraction,sample.delta_chi2);

      using ScanStatus = RelActCalcAutoImp::ProfileLikelihood::ScanStatus;
      if( scan.status == ScanStatus::Canceled )
      {
        solution.m_status = RelActAutoSolution::Status::UserCanceled;
        solution.m_error_message = scan.diagnostic.empty()
                                 ? "Mass-fraction profiling was canceled." : scan.diagnostic;
        return;
      }

      if( scan.status == ScanStatus::BetterBaseline )
      {
        const auto discovery_disposition
            = RelActCalcAutoImp::ProfileLikelihood::baseline_discovery_disposition(
                                                        baseline_restart_count);
        if( discovery_disposition
            == RelActCalcAutoImp::ProfileLikelihood::BaselineDiscoveryDisposition::RejectAfterRestart )
        {
          profile.status = ProfileStatus::Failed;
          profile.message = "A profile found a still better point after the baseline-reselection"
                            " budget was exhausted; the selected baseline is still not optimal."
                            + (scan.diagnostic.empty() ? std::string()
                                                       : std::string(" ") + scan.diagnostic);
          solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
          solution.m_warnings.push_back( target->symbol + " mass-fraction profile: " + profile.message );
          // The proof of non-optimality taints every interval of this pass, not just the
          // discovering target's - their delta-chi2 baselines are the same nominal.  Say so once.
          if( !budget_exhaustion_warned )
          {
            budget_exhaustion_warned = true;
            solution.m_warnings.push_back( "The reported baseline is known not to be the global"
                " optimum (the baseline-reselection budget was exhausted while profiles kept"
                " finding better points), so every profile interval of this solution is measured"
                " against a non-optimal nominal and may be unreliable." );
          }
          continue;
        }

        if( !better_conditional )
        {
          profile.status = ProfileStatus::Failed;
          profile.message = "A better profile point was reported but could not seed baseline reselection.";
          solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
          solution.m_warnings.push_back( target->symbol + " mass-fraction profile: " + profile.message );
          continue;
        }

        if( !std::isfinite(better_conditional->m_chi2) )
        {
          profile.status = ProfileStatus::Failed;
          profile.message = "A better profile point had a non-finite frozen physical objective.";
          solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
          solution.m_warnings.push_back( target->symbol + " mass-fraction profile: " + profile.message );
          continue;
        }

        DeferredBaselineDiscovery deferred;
        deferred.rank.full_objective = better_conditional->m_chi2;
        deferred.rank.semantic_key = semantic_curve_key(solution.m_options,curve_index)
                                   + "|target=" + semantic_source_key(target_source);
        deferred.conditional_solution = std::move(better_conditional);
        deferred.curve_index = curve_index;
        deferred.target_symbol = target->symbol;
        deferred_baseline_discoveries.push_back(std::move(deferred));

        // This entry is temporary: a successful reselection replaces the solution and restarts all
        // profiles.  If reselection itself fails, retaining a structured Failed result is more
        // useful than silently dropping the quantity which discovered the better point.
        profile.status = ProfileStatus::Failed;
        profile.message = "A better physical point is pending final candidate selection.";
        solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
        continue;
      }

      if( scan.status == ScanStatus::Failed || scan.status == ScanStatus::FitCapReached )
      {
        profile.status = ProfileStatus::Failed;
        profile.message = scan.diagnostic.empty()
                        ? "The bounded profile could not be classified." : scan.diagnostic;
        solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
        solution.m_warnings.push_back( target->symbol + " mass-fraction profile failed: " + profile.message );
        continue;
      }

      const double confidences[2] = { 0.6827, 0.95 };
      for( size_t level = 0; level < 2; ++level )
      {
        Interval interval;
        interval.confidence_level = confidences[level];
        interval.delta_chi2 = thresholds[level];
        interval.lower = scan.intervals[level].lower.reported_fraction;
        interval.upper = scan.intervals[level].upper.reported_fraction;
        interval.lower_kind = scan.intervals[level].lower.likelihood_crossing
                              ? EndpointKind::LikelihoodCrossing : lower_bound_kind;
        interval.upper_kind = scan.intervals[level].upper.likelihood_crossing
                              ? EndpointKind::LikelihoodCrossing : upper_bound_kind;
        if( interval.lower > interval.upper )
          std::swap(interval.lower, interval.upper);
        profile.intervals.push_back(interval);
      }

      if( scan.status == ScanStatus::NonIdentifiable )
      {
        profile.status = ProfileStatus::NonIdentifiable;
      }else if( scan.status == ScanStatus::BoundaryLimited )
      {
        profile.status = ProfileStatus::BoundaryLimited;
      }else
      {
        profile.status = ProfileStatus::Complete;
      }
      profile.message = scan.diagnostic;

      solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;

      if( RelActCalcAutoImp::profile_stats_enabled() )
      {
        const size_t ceres_solves = profile_host
            ? (profile_host->num_conditional_solves() - host_solves_before) : 0;
        std::cerr << "profile-stats " << target->symbol << " curve=" << curve_index
                  << " points=" << profile.num_fits
                  << " ceres_solves=" << ceres_solves
                  << " seconds=" << std::chrono::duration<double>(
                                       std::chrono::steady_clock::now() - profile_start).count()
                  << " status=" << static_cast<int>(profile.status) << std::endl;
      }
    }
  }

  if( deferred_baseline_discoveries.empty() )
    return;

  std::vector<RelActCalcAutoImp::ProfileLikelihood::PendingBaselineDiscovery> ranks;
  ranks.reserve(deferred_baseline_discoveries.size());
  for( const DeferredBaselineDiscovery &discovery : deferred_baseline_discoveries )
    ranks.push_back(discovery.rank);
  const std::optional<size_t> best_index
      = RelActCalcAutoImp::ProfileLikelihood::best_pending_baseline_discovery_index(ranks);
  assert( best_index && (*best_index < deferred_baseline_discoveries.size()) );

  const auto fail_deferred_profiles = [&]( const std::string &message ) {
    for( const DeferredBaselineDiscovery &discovery : deferred_baseline_discoveries )
    {
      Profile &profile
          = solution.m_mass_fraction_profiles[discovery.curve_index][discovery.target_symbol];
      profile.status = ProfileStatus::Failed;
      profile.message = message;
      solution.m_warnings.push_back(
          discovery.target_symbol + " mass-fraction profile: " + message );
    }
  };

  if( !best_index )
  {
    fail_deferred_profiles( "No finite profile-discovered point was available for baseline"
                            " reselection." );
    return;
  }

  DeferredBaselineDiscovery &selected = deferred_baseline_discoveries[*best_index];
  assert( selected.conditional_solution );
  if( !selected.conditional_solution )
  {
    fail_deferred_profiles( "The selected profile-discovered point could not seed baseline"
                            " reselection." );
    return;
  }

  // A profile-discovered point is a deterministic-search trigger.  Polish the complete named
  // candidate matrix on the same frozen ROIs, then restart every requested profile so covariance
  // and thresholds use the selected baseline.
  Options restart_options = solution.m_options;
  for( RoiRange &roi : restart_options.rois )
    roi.range_limits_type = RoiRange::RangeLimitsType::Fixed;
  RelActAutoSolution restarted;
  try
  {
    restarted = RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres(
        restart_options,foreground,background,frozen_profile_drf,frozen_profile_peaks,det_type,
        cancel_calc,RelActCalcAutoImp::SearchSeedVariant::Default,
        true,selected.conditional_solution.get(),true,&frozen_profile_peak_ranges,
        &frozen_profile_model_policy,true,/*may_host_profile=*/true );
  }catch( const std::exception &e )
  {
    fail_deferred_profiles(
        std::string("Baseline reselection after profile-discovered points failed: ") + e.what() );
    return;
  }

  if( restarted.m_status == RelActAutoSolution::Status::UserCanceled
      || (cancel_calc && cancel_calc->load()) )
  {
    solution.m_status = RelActAutoSolution::Status::UserCanceled;
    solution.m_error_message = "Mass-fraction baseline reselection was canceled.";
    return;
  }

  if( RelActAutoSolution::is_usable_status(restarted.m_status) )
  {
    const bool same_frozen_objective = matches_frozen_profile_objective(restarted);
    assert( same_frozen_objective );
    if( !same_frozen_objective )
    {
      fail_deferred_profiles( "Baseline reselection changed the frozen branching-ratio,"
                              " gamma-membership, or selected-model policy." );
      return;
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
    fail_deferred_profiles( "Profiles found better physical points, but final candidate selection"
                            " did not establish a better unconstrained baseline." );
    return;
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

  // The reselection trail must survive the move below: `restarted` is a fresh solve carrying only
  // its own warnings, so without this a multi-hop chain would report only its final hop.
  // CONTRACT: every warning this flow wants carried across restarts - the hop summaries and the
  // seed-fidelity note pushed below - must begin with exactly "Baseline reselection ", because
  // that prefix is the only thing selecting them here.  Reword the prefix in ALL places or hops
  // silently vanish from the trail.
  std::vector<std::string> reselection_trail;
  for( const std::string &warning : solution.m_warnings )
    if( warning.rfind("Baseline reselection ",0) == 0 )
      reselection_trail.push_back( warning );

  solution = std::move(restarted);
  for( const std::string &warning : reselection_trail )
    solution.m_warnings.push_back( warning );
  solution.m_warnings.push_back( "Baseline reselection "
      + std::to_string(baseline_restart_count + 1) + " of "
      + std::to_string(RelActCalcAutoImp::ProfileLikelihood::sm_max_baseline_restarts) + ": "
      + std::to_string(num_discoveries)
      + " mass-fraction profile(s) found better points; deterministic final candidate selection"
        " used the best frozen-objective seed (" + selected_target + "; discoveries: "
      + discovery_summary + ") and lowered the objective from "
      + std::to_string(old_objective) + " to " + std::to_string(solution.m_chi2)
      + ", so covariance and all profiles were restarted against the new baseline." );

  // The seed is a physically realizable point of the SAME unconstrained problem, so a re-solve
  // warm-started from it finishing ABOVE it is proof the semantic warm start did not reach the
  // discovered basin (or the solve terminated on a plateau).  The improved baseline is still kept
  // - it beat the old one - but say so loudly: this is the signature that separates genuine
  // basin chains from a restart crawling toward one basin it keeps rediscovering.
  if( std::isfinite(seed_objective)
      && (solution.m_chi2 > seed_objective + baseline_improvement_tolerance) )
    solution.m_warnings.push_back( "Baseline reselection "
        + std::to_string(baseline_restart_count + 1) + " did not reach its seed's known objective ("
        + std::to_string(solution.m_chi2) + " vs " + std::to_string(seed_objective)
        + "); the warm start likely failed to enter the discovered basin, which may be"
          " rediscovered on the next profile pass." );
  add_merged_single_curve_comparison( solution,restart_options,foreground,background,
                                      solution.m_drf,solution.m_spectrum_peaks,
                                      det_type,cancel_calc );
  add_mass_fraction_profiles( solution,foreground,background,solution.m_drf,
                              solution.m_spectrum_peaks,
                              det_type,cancel_calc,baseline_restart_count+1 );
}
}//anonymous namespace

#endif //InterSpec_RelActCalcAuto_Profile_imp_hpp
