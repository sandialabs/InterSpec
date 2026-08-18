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
 coordinate is fixed is the central design trade of this file: constraining the REPORTED coordinate
 - the one an interval is quoted in, after Pu-242 correlation, age correction and renormalization -
 keeps the statistic exact, but needs a residual row whose Jacobian couples the whole problem.
 Pinning the target's own fit PARAMETER instead is roughly two orders of magnitude cheaper and is
 what ships; it reports intervals that are somewhat too narrow, by an amount the endpoint-refit
 harness measures rather than assumes.

 Conditional points run one way: PINNED, in place.  The winning solve hands its converged
 `ceres::Problem` forward on the solution; the target's own parameter is removed from the tangent
 space with a `ceres::SubsetManifold`, and each point loads the previous point's parameter vector,
 writes the pinned coordinate, and re-solves.  No residual is added, so the Jacobian sparsity is
 untouched and a point costs about a second rather than minutes.  The statistic is
 `min over {a == a0}` rather than `min over {q == q0}` - an upper bound, so intervals come out
 somewhat narrow, which is why the endpoint-refit harness measures the shortfall directly and the
 level-set refinement recovers most of it.  There is no slower fallback engine: a target
 `pin_index_for` refuses, or a solve that retained no problem to pin within, reports a structured
 `Failed` profile carrying the reason.
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

  // A profile is a conditional optimization of the selected physical problem, not a request to
  // classify a nearby piecewise model again.  Snapshot every continuum rank/active-set decision
  // and every BR-nuisance interval (including an intentionally empty interval set) from the main
  // solution.  Conditional solves and a one-time better-baseline search both consume these exact
  // values and defer their ordinary continuum outer check.
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
  std::vector<PeakFit::ContinuumFitPolicy> frozen_profile_continuum_policies;
  frozen_profile_continuum_policies.reserve(solution.m_cost_functor->m_energy_ranges.size());
  for( const RelActCalcAutoImp::RoiRangeChannels &roi
       : solution.m_cost_functor->m_energy_ranges )
    frozen_profile_continuum_policies.push_back(roi.frozen_continuum_policy);
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
  // (and a one-time better-baseline search) uses precisely the same nuisance-parameter manifold.
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
        && candidate.m_cost_functor->continuum_policies_match(
                                                frozen_profile_continuum_policies)
        && (candidate.m_cost_functor->m_peak_ranges_with_uncert
                                                == frozen_profile_peak_ranges)
        && (candidate.m_frozen_gamma_membership_hash == frozen_profile_gamma_hash)
        && (candidate.m_frozen_model_policy_hash == frozen_profile_model_hash);
  };

  // Do not let the first weak quantity visited choose the one permitted baseline restart.  A
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

      bool ratio_constrained = false;
      for( const RelEffCurveInput::ActRatioConstraint &constraint : base_curve.act_ratio_constraints )
        ratio_constrained = ratio_constrained || (constraint.controlling_source == target_source)
                            || (constraint.constrained_source == target_source);
      const bool has_pu_reporting_transform
          = (target->atomicNumber == 94)
            && (base_curve.pu242_correlation_method
                != RelActCalc::PuCorrMethod::NotApplicable);
      (void)ratio_constrained;   //the pin walks ratio chains to their root; see `pin_index_for`
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

      // --- The production engine: pin a fit parameter and re-solve in place ---------------------
      //
      // The alternative - the deleted augmented-Lagrangian engine - constrained the REPORTED
      // coordinate exactly, but did so with a residual row whose Jacobian touched every activity
      // and every age of the element, so one dense row coupled the whole problem and each
      // conditional point cost seconds.  Restricting a parameter instead adds no residual and
      // leaves the sparsity untouched.  Measured on JRC Pu70: 1.09 s/point against a
      // whole-scan-set cost of 6,532 s.
      //
      // The statistic is not identical.  Pinning gives `min over {a == a0}`, which is a feasible
      // point of the reported-coordinate constraint set and therefore an UPPER bound on
      // `min over {q == q0}`: delta-chi2 is overstated at a given q, the threshold is crossed too
      // early, and the interval comes out somewhat NARROW.  That is accepted only because it is
      // measured - see the endpoint-refit harness - and because the one case where the two
      // genuinely diverge (correlation-generated Pu-242, whose reported value depends on Pu-239 as
      // well) is no longer a profile target at all.
      std::optional<size_t> pin_index;
      std::string pin_refusal;
      if( profile_host )
      {
        const RelActCalcAutoImp::ProfileConditionalHost::PinSelection selection
            = profile_host->pin_index_for( target_source, curve_index );
        pin_index = selection.index;
        pin_refusal = selection.why_not;
      }

      RelActCalcAutoImp::ProfileLikelihood::ScanResult scan;
      bool used_pinned_engine = false;

      if( pin_index )
      {
        // One manifold for the whole of this target's scan, restored on every exit path.
        struct PinGuard
        {
          RelActCalcAutoImp::ProfileConditionalHost *host;
          ~PinGuard(){ if( host ) host->restore_optimum(); }
        } pin_guard{ profile_host.get() };

        used_pinned_engine = profile_host->set_pin( *pin_index );
        if( !used_pinned_engine )
        {
          pin_guard.host = nullptr;
          pin_refusal = "The profile manifold could not be installed.";
        }

        if( used_pinned_engine )
        {
          const std::vector<double> &optimum = profile_host->optimum();
          const double baseline_parameter = profile_host->optimum_parameter( *pin_index );
          const std::pair<double,double> parameter_box
              = profile_host->parameter_bounds( *pin_index );

          // The slope of the ACHIEVED reported quantity with respect to the pinned coordinate.
          //
          // This must be the CONDITIONAL response `(Cov u)_k / Cov_kk`, not the partial derivative
          // `dq/da` at fixed everything-else.  A conditional solve reoptimizes the siblings, and the
          // covariance cross-terms are the entire difference between the two: a frozen finite
          // difference under-measured the true response by about sevenfold on Pu-240, which put the
          // opening probe in the wrong place.
          //
          // Alongside it, `r^2 = (Cov u)_k^2 / (Cov_kk * u.Cov.u)` - the squared correlation between
          // the pinned parameter and the delta-method linearization of `q`, a cosine in the
          // covariance metric.  In the quadratic regime the pinned profile overstates delta chi2 by
          // exactly `1/r^2`, so the interval is narrow by `|r|`: `r^2 == 1` means the pin IS the
          // reported coordinate and there is no bias at all (the sole-range-nuclide carrier case),
          // while a smaller value predicts, for free, how much the level-set refinement below has to
          // recover.
          double reported_per_parameter = 0.0;
          double pin_alignment_r2 = std::numeric_limits<double>::quiet_NaN();
          std::vector<double> reported_gradient;
          const bool have_covariance
              = (solution.m_covariance.size() == optimum.size())
                && (*pin_index < optimum.size());
          if( have_covariance )
          {
            try
            {
              reported_gradient = profile_host->reported_gradient( target_source,curve_index,optimum );
              const size_t num_pars = optimum.size();
              double cov_u_k = 0.0, u_cov_u = 0.0;
              for( size_t row = 0; row < num_pars; ++row )
              {
                double cov_u_row = 0.0;
                for( size_t col = 0; col < num_pars; ++col )
                  cov_u_row += solution.m_covariance[row][col]*reported_gradient[col];
                u_cov_u += reported_gradient[row]*cov_u_row;
                if( row == *pin_index )
                  cov_u_k = cov_u_row;
              }
              const double cov_kk = solution.m_covariance[*pin_index][*pin_index];
              if( std::isfinite(cov_u_k) && std::isfinite(cov_kk) && (cov_kk > 0.0) )
                reported_per_parameter = cov_u_k/cov_kk;
              if( std::isfinite(u_cov_u) && (u_cov_u > 0.0) && (cov_kk > 0.0) )
                pin_alignment_r2 = (cov_u_k*cov_u_k)/(cov_kk*u_cov_u);
            }catch( const std::exception & )
            {
              reported_per_parameter = 0.0;
            }
          }

          // Fall back to the frozen partial derivative only when the covariance is unavailable; it
          // is a poor estimate of the response, but it is better than no scale at all.
          if( !(std::fabs(reported_per_parameter) > 0.0) )
          {
            const double step = (std::max)( 1.0e-8, 1.0e-5*std::fabs(baseline_parameter) );
            std::vector<double> probe = optimum;
            probe[*pin_index] = (std::min)( (std::max)(baseline_parameter + step,
                                                       parameter_box.first), parameter_box.second );
            const double actual_step = probe[*pin_index] - baseline_parameter;
            try
            {
              if( actual_step != 0.0 )
              {
                const double probe_fraction = profile_host->reported_mass_fraction(
                                                    target_source, curve_index, probe );
                reported_per_parameter = (probe_fraction - nominal.fraction)/actual_step;
              }
            }catch( const std::exception & )
            {
              reported_per_parameter = 0.0;
            }
          }

          // Scan limits in the pinned coordinate.  A genuine Ceres box edge is a real feasible
          // limit and is used directly.  Where the box is open there is no such limit at all - a
          // mass fraction approaches 1 asymptotically as its activity grows, so the "distance to
          // the feasible bound" is genuinely infinite and any cap is a choice rather than a fact.
          //
          // The choice is deliberately MODEST: four times the local linear estimate of the distance
          // to the reported limit.  For a target at fraction 0.7 that is roughly an activity of
          // 5.7x nominal, i.e. a fraction near 0.93 - far outside any plausible 95% crossing, while
          // still a place the optimizer can actually solve.  A larger cap (an earlier draft used
          // 20x) buys nothing and costs a great deal: it puts the outermost ladder point at ~30x
          // nominal activity, a badly conditioned corner where each conditional solve is slow and
          // the sample it produces is far outside the region the fit cares about.
          //
          // The honest bound detection is the read-back below, never this cap.
          //
          // Whatever the estimate says, the excursion is capped at ten times the parameter's own
          // natural scale.  Without that cap a target already near saturation - a mass fraction of
          // 0.94, where the remaining 0.06 costs a THOUSANDFOLD activity increase - produces a span
          // of thousands, and then even the modest `0.02*span` opening probe lands at forty times
          // nominal activity, in a corner where the model cannot be evaluated at all and every
          // sample on the side is discarded.  Ten times the parameter scale is far outside any
          // plausible crossing while remaining somewhere the optimizer can work; a side that
          // genuinely does not reach the threshold inside it is correctly `BoundaryLimited`.
          // Ceres represents "no bound" with `numeric_limits<double>::max()` and `lowest()`, NOT
          // with an infinity - so `std::isfinite` accepts the sentinel as a genuine bound, and the
          // scan then marches the pin out toward 1e308.  Every sample there is unevaluable, which
          // is how an entire side comes back empty on a problem whose other side is perfectly
          // clean: five converged points below the optimum, nothing at all above it.
          const auto has_real_bound = []( const double edge ){
            return std::isfinite(edge)
                   && (std::fabs(edge) < 0.5*std::numeric_limits<double>::max());
          };

          const double parameter_scale = (std::max)( 1.0e-9,
              has_real_bound(parameter_box.first) ? (baseline_parameter - parameter_box.first)
                                                  : std::fabs(baseline_parameter) );

          const auto scan_limit = [&]( const double box_edge, const double fraction_limit ){
            if( has_real_bound(box_edge) )
              return box_edge;

            const double direction = (fraction_limit < nominal.fraction) ? -1.0 : 1.0;
            const double ceiling = baseline_parameter + direction*10.0*parameter_scale;
            if( (reported_per_parameter != 0.0) && std::isfinite(reported_per_parameter) )
            {
              const double linear = (fraction_limit - nominal.fraction)/reported_per_parameter;
              if( std::isfinite(linear) && (linear != 0.0) )
              {
                const double estimate = baseline_parameter + 4.0*linear;
                return (direction > 0.0) ? (std::min)(estimate,ceiling)
                                         : (std::max)(estimate,ceiling);
              }
            }
            return ceiling;
          };

          if( RelActCalcAutoImp::profile_stats_enabled() )
            std::cerr << "profile-pin " << target->symbol << " curve=" << curve_index
                      << " index=" << *pin_index
                      << " a0=" << baseline_parameter
                      << " box=[" << parameter_box.first << "," << parameter_box.second << "]"
                      << " dq_da=" << reported_per_parameter
                      << " r2=" << pin_alignment_r2
                      << " q0=" << nominal.fraction << std::endl;

          double lower_parameter = scan_limit( parameter_box.first, control_lower );
          double upper_parameter = scan_limit( parameter_box.second, control_upper );
          if( lower_parameter > baseline_parameter )
            lower_parameter = baseline_parameter;
          if( upper_parameter < baseline_parameter )
            upper_parameter = baseline_parameter;

          // The local one-sigma, carried into the pinned coordinate so the probe starts somewhere
          // sensible.  Absent for a weak quantity, which is the common case; the placement is
          // self-calibrating and falls back to a fraction of the span.
          double parameter_one_sigma = 0.0;
          if( nominal.covariance_one_sigma && (reported_per_parameter != 0.0)
              && std::isfinite(reported_per_parameter) )
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
            if( have_covariance && !reported_gradient.empty() && std::isfinite(pinned_fraction)
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

      // Pinning is the only conditional mechanism there is.  When a target has no coordinate to pin
      // - a mass fraction fixed by input, a slot held constant by the selected model, or no retained
      // problem to pin within - there is nothing to fall back to, and saying so plainly is the whole
      // point: a structured `Failed` carrying the reason is honest, where silence would read as a
      // quantity nobody thought to profile.
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
          profile.message = "A profile found a second better point after baseline reselection;"
                            " the selected baseline is still not optimal."
                            + (scan.diagnostic.empty() ? std::string()
                                                       : std::string(" ") + scan.diagnostic);
          solution.m_mass_fraction_profiles[curve_index][target->symbol] = profile;
          solution.m_warnings.push_back( target->symbol + " mass-fraction profile: " + profile.message );
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
        true,selected.conditional_solution.get(),true,
        &frozen_profile_continuum_policies,0,true,&frozen_profile_peak_ranges,
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
      fail_deferred_profiles( "Baseline reselection changed the frozen continuum, branching-ratio,"
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
  solution = std::move(restarted);
  solution.m_warnings.push_back( std::to_string(num_discoveries)
      + " mass-fraction profile(s) found better points; deterministic final candidate selection"
        " used the best frozen-objective seed (" + selected_target + "; first-pass discoveries: "
      + discovery_summary + ") and lowered the objective from "
      + std::to_string(old_objective) + " to " + std::to_string(solution.m_chi2)
      + ", so covariance and all profiles were restarted once." );
  add_merged_single_curve_comparison( solution,restart_options,foreground,background,
                                      solution.m_drf,solution.m_spectrum_peaks,
                                      det_type,cancel_calc );
  add_mass_fraction_profiles( solution,foreground,background,solution.m_drf,
                              solution.m_spectrum_peaks,
                              det_type,cancel_calc,baseline_restart_count+1 );
}
}//anonymous namespace

#endif //InterSpec_RelActCalcAuto_Profile_imp_hpp
