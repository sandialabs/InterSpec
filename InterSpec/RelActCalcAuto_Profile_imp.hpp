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

 The conditional equality is imposed in the reported post-correction coordinate, and all nuisance
 parameters are reoptimized.  Where no reporting transform is present, the production exact
 mass-fraction parameterization is used.  Pu-correlation profiles use a transient augmented-
 Lagrangian row evaluated after correlation and renormalization.  Conditional solves call
 solve_ceres directly with fixed ROIs and profiling disabled, avoiding recursion; reported
 likelihood differences always exclude the augmented-Lagrangian row.
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

  if( !RelActAutoSolution::is_usable_status(solution.m_status) )
    return;

  const size_t num_curves = solution.m_options.rel_eff_curves.size();
  solution.m_mass_fraction_profiles.clear();
  solution.m_mass_fraction_profiles.resize( num_curves );

  // Automatic profiling is the single most expensive thing a solve can do - up to
  // `max_conditional_fits` conditional optimizations per weak quantity - so it belongs to the opt-in
  // robust budget.  An *explicit* per-nuclide request is a different thing: the user asked for that
  // specific quantity, and it is honored whatever the budget.
  const bool automatic_profiling_enabled = solution.m_options.robust_solve
                                        && solution.m_options.auto_profile_weak_mass_fractions;
  bool any_forced_profile = false;
  for( const RelEffCurveInput &curve : solution.m_options.rel_eff_curves )
    for( const NucInputInfo &nuc : curve.nuclides )
      any_forced_profile = any_forced_profile || nuc.force_profile_mass_fraction;
  if( !automatic_profiling_enabled && !any_forced_profile )
    return;

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

    // The Pu-242 correlation creates a reported quantity which deliberately is not a fitted
    // NucInputInfo.  Treat it as a first-class automatic profile candidate: its covariance and
    // weakness classification are available through mass_enrichment_result(), and the conditional
    // solver below constrains the true post-correlation coordinate.  Keep the optional input
    // pointer solely for the persisted per-input force flag.
    struct ProfileTarget
    {
      const SandiaDecay::Nuclide *nuclide = nullptr;
      const NucInputInfo *input = nullptr;
    };
    std::vector<ProfileTarget> profile_targets;
    profile_targets.reserve( base_curve.nuclides.size() + 1 );
    for( const NucInputInfo &input : base_curve.nuclides )
    {
      const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(input.source);
      if( nuclide )
        profile_targets.push_back( {nuclide,&input} );
    }

    if( base_curve.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      const SandiaDecay::Nuclide * const pu242 = db ? db->nuclide("Pu242") : nullptr;
      const bool already_present = pu242 && std::any_of(
          begin(profile_targets),end(profile_targets),
          [pu242]( const ProfileTarget &candidate ){ return candidate.nuclide == pu242; } );
      if( pu242 && !already_present )
        profile_targets.push_back( {pu242,nullptr} );
    }

    for( const ProfileTarget &profile_target : profile_targets )
    {
      const SandiaDecay::Nuclide * const target = profile_target.nuclide;
      if( !target )
        continue;
      const SrcVariant target_source(target);
      const bool forced = profile_target.input
                       && profile_target.input->force_profile_mass_fraction;

      size_t same_element_count = 0;
      for( const NucInputInfo &other : base_curve.nuclides )
      {
        const SandiaDecay::Nuclide * const other_nuc = RelActCalcAuto::nuclide(other.source);
        same_element_count += (other_nuc && (other_nuc->atomicNumber == target->atomicNumber));
      }
      if( same_element_count < 2 )
      {
        // An explicit request always produces a structured result.  With no second isotope of the
        // element in the model, the normalized fraction is identically one and there is no
        // conditional optimization to perform.
        if( forced )
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
        // A generated correlation quantity has no persisted force flag.  If automatic profiling
        // was requested, do not silently lose it when even its structured nominal result cannot be
        // formed: retain the otherwise usable main fit and expose a failed automatic profile.
        const bool generated_automatic = !profile_target.input
                                      && automatic_profiling_enabled;
        if( forced || generated_automatic )
        {
          Profile failed;
          failed.status = ProfileStatus::Failed;
          failed.reason = forced ? ProfileReason::Forced : ProfileReason::AutomaticWeak;
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
      const bool target_has_activity_coordinate_input
          = profile_target.input
            && (profile_target.input->min_rel_act.has_value()
                || profile_target.input->max_rel_act.has_value()
                || profile_target.input->starting_rel_act.has_value());
      const bool use_reported_fraction_constraint
          = has_pu_reporting_transform || ratio_constrained
            || target_has_activity_coordinate_input;
      // An unconstrained fitted source can use the exact production mass-fraction
      // reparameterization.  Pu correlation must target the post-correction result, while an
      // activity-ratio target must vary its controller rather than trying to free the constrained
      // activity.  A target carrying an activity box/start must likewise remain an ordinary
      // activity coordinate: replacing it with a sigma-block mass-fraction coordinate would
      // silently discard the callers physical activity feasible set.  These cases therefore use
      // the reported-coordinate augmented Lagrangian below.
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

      Options base_profile_options = solution.m_options;
      base_profile_options.auto_profile_weak_mass_fractions = false;
      base_profile_options.auto_simplify_model = false;
      for( RoiRange &roi : base_profile_options.rois )
        roi.range_limits_type = RoiRange::RangeLimitsType::Fixed;
      for( RelEffCurveInput &curve : base_profile_options.rel_eff_curves )
        for( NucInputInfo &nuc : curve.nuclides )
          nuc.force_profile_mass_fraction = false;

      RelEffCurveInput &base_profile_curve = base_profile_options.rel_eff_curves[curve_index];
      if( !use_reported_fraction_constraint )
      {
        base_profile_curve.mass_fraction_constraints.erase(
          std::remove_if( begin(base_profile_curve.mass_fraction_constraints),
                          end(base_profile_curve.mass_fraction_constraints),
                          [target]( const RelEffCurveInput::MassFractionConstraint &constraint ){
                            return constraint.nuclide == target;
                          } ),
          end(base_profile_curve.mass_fraction_constraints) );
      }

      for( NucInputInfo &nuc_info : base_profile_curve.nuclides )
      {
        const bool controlled_by_activity_ratio = std::any_of(
            begin(base_profile_curve.act_ratio_constraints),
            end(base_profile_curve.act_ratio_constraints),
            [&nuc_info]( const RelEffCurveInput::ActRatioConstraint &constraint ){
              return constraint.constrained_source == nuc_info.source;
            } );
        if( curve_index < solution.m_rel_activities.size() )
        {
          const auto pos = std::find_if( begin(solution.m_rel_activities[curve_index]),
                                         end(solution.m_rel_activities[curve_index]),
                       [&nuc_info]( const NuclideRelAct &fit ){ return fit.source == nuc_info.source; } );
          if( pos != end(solution.m_rel_activities[curve_index]) )
          {
            nuc_info.age = pos->age;
            if( use_reported_fraction_constraint && !nuc_info.starting_rel_act
                && !controlled_by_activity_ratio
                && std::isfinite(pos->rel_activity)
                && (pos->rel_activity > 0.0) )
              nuc_info.starting_rel_act = pos->rel_activity;
          }
        }
      }

      static constexpr size_t max_conditional_fits = 32;
      struct ConditionalWarmPoint
      {
        Options options;
        RelActAutoSolution solution;
      };
      std::vector<double> lower_warm_controls, upper_warm_controls;
      std::vector<ConditionalWarmPoint> lower_warm_points, upper_warm_points;
      std::unique_ptr<RelActAutoSolution> better_conditional;
      const double baseline_improvement_tolerance
          = RelActCalcAutoImp::ProfileLikelihood::baseline_improvement_tolerance(
                                              solution.m_chi2,solution.m_cov_scale);

      auto evaluate = [&]( const double requested_control,
                           const RelActCalcAutoImp::ProfileLikelihood::Direction direction,
                           const size_t remaining_fit_budget )
                           -> RelActCalcAutoImp::ProfileLikelihood::Evaluation {
        using Eval = RelActCalcAutoImp::ProfileLikelihood::Evaluation;
        using EvalStatus = RelActCalcAutoImp::ProfileLikelihood::EvaluationStatus;
        Eval evaluation;
        evaluation.num_fits = 0;
        if( cancel_calc && cancel_calc->load() )
        {
          evaluation.status = EvalStatus::Canceled;
          evaluation.diagnostic = "Mass-fraction profiling was canceled.";
          return evaluation;
        }

        std::vector<double> &warm_controls = (direction
                   == RelActCalcAutoImp::ProfileLikelihood::Direction::Lower)
                   ? lower_warm_controls : upper_warm_controls;
        std::vector<ConditionalWarmPoint> &warm_points = (direction
                   == RelActCalcAutoImp::ProfileLikelihood::Direction::Lower)
                   ? lower_warm_points : upper_warm_points;
        const size_t warm_index
            = RelActCalcAutoImp::ProfileLikelihood::nearest_control_index(
                warm_controls,requested_control );
        Options trial_options = (warm_index < warm_points.size())
                              ? warm_points[warm_index].options : base_profile_options;
        // The unconstrained baseline is the nearest available point until this direction has a
        // successful conditional fit.  Brent and the nested 68/95 passes are not chronological,
        // so every later request selects the closest cached coordinate rather than whichever fit
        // happened to run last.
        std::unique_ptr<RelActAutoSolution> semantic_warm
            = std::make_unique<RelActAutoSolution>(
                (warm_index < warm_points.size())
                  ? warm_points[warm_index].solution : solution );

        RelActAutoSolution conditional;
        double actual_reported = std::numeric_limits<double>::quiet_NaN();
        double previous_reported = std::numeric_limits<double>::quiet_NaN();
        double last_reported_change = std::numeric_limits<double>::infinity();
        // The nearest successful conditional point retains the exact multiplier and penalty used
        // by its final solve in `trial_options`.  Reuse that AL state as well as its fitted
        // parameters: resetting every nearby Brent request to (lambda,mu)=(0,1e12) made each
        // trace-isotope point repeat all three continuation solves and could exhaust the 32-fit
        // profile cap on a single crossing.  The previous multiplier is a natural warm estimate of
        // the equality KKT multiplier; the existing convergence check and at-most-three-pass
        // continuation remain unchanged.
        const Options::ProfileOnlyMassFractionConstraint *warm_profile_constraint = nullptr;
        if( trial_options.profile_only_mass_fraction_constraint )
        {
          const Options::ProfileOnlyMassFractionConstraint &candidate
              = *trial_options.profile_only_mass_fraction_constraint;
          if( (candidate.rel_eff_curve_index == curve_index)
              && (candidate.nuclide == target)
              && std::isfinite(candidate.lagrange_multiplier)
              && std::isfinite(candidate.penalty) && (candidate.penalty > 0.0) )
            warm_profile_constraint = &candidate;
        }
        double lagrange_multiplier = warm_profile_constraint
            ? warm_profile_constraint->lagrange_multiplier : 0.0;
        // A reported fraction is O(1), while the spectrum objective commonly has thousands of
        // rows.  A cold start at 1e12 normally satisfies the equality in one conditional solve;
        // multiplier updates remain available for genuinely difficult points.
        double penalty = warm_profile_constraint ? warm_profile_constraint->penalty : 1.0e12;
        // Preserve every valid interior input coordinate exactly, even for trace-isotope windows
        // far below one ppm.  Only mathematical zero/one need an inward representative for the
        // open production parameterization, and both that displacement and the AL convergence
        // tolerance scale with this targets actual feasible span.
        const double requested_for_solve
            = RelActCalcAutoImp::ProfileLikelihood::conditional_solve_control(
                       requested_control,control_baseline,control_lower,control_upper,
                       2.0e-6,use_reported_fraction_constraint ? 0.0 : 2.0e-6);
        const double equality_tolerance
            = RelActCalcAutoImp::ProfileLikelihood::conditional_equality_tolerance(
                       requested_for_solve,control_baseline,control_lower,control_upper);
        const size_t max_augmented_lagrangian_passes
            = use_reported_fraction_constraint ? (std::min)(size_t(3),remaining_fit_budget)
                                               : (std::min)(size_t(1),remaining_fit_budget);
        bool equality_satisfied = false;

        for( size_t al_pass = 0; al_pass < max_augmented_lagrangian_passes; ++al_pass )
        {
          RelEffCurveInput &trial_curve = trial_options.rel_eff_curves[curve_index];
          bool appended_exact_constraint = false;
          if( use_reported_fraction_constraint )
          {
            Options::ProfileOnlyMassFractionConstraint constraint;
            constraint.rel_eff_curve_index = curve_index;
            constraint.nuclide = target;
            constraint.reported_fraction = requested_for_solve;
            constraint.lagrange_multiplier = lagrange_multiplier;
            constraint.penalty = penalty;
            trial_options.profile_only_mass_fraction_constraint = constraint;
          }else
          {
            trial_options.profile_only_mass_fraction_constraint.reset();
            // Exact one leaves no positive remainder for an unconstrained sibling.  Its inward
            // limit is numerically equivalent while the returned endpoint remains physical.
            RelEffCurveInput::MassFractionConstraint fixed;
            fixed.nuclide = target;
            fixed.lower_mass_fraction = requested_for_solve;
            fixed.upper_mass_fraction = requested_for_solve;
            trial_curve.mass_fraction_constraints.push_back( fixed );
            appended_exact_constraint = true;
          }

          try
          {
            ++evaluation.num_fits;
            conditional = RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres(
                              trial_options, foreground, background, frozen_profile_drf,
                              frozen_profile_peaks, det_type, cancel_calc,
                              RelActCalcAutoImp::SearchSeedVariant::Default,
                              false, semantic_warm.get(), false,
                              &frozen_profile_continuum_policies,0,true,
                              &frozen_profile_peak_ranges,
                              &frozen_profile_model_policy,true );
          }catch( const std::exception &e )
          {
            if( appended_exact_constraint )
              trial_curve.mass_fraction_constraints.pop_back();
            evaluation.status = EvalStatus::Failed;
            evaluation.diagnostic = std::string("A conditional profile solve threw: ") + e.what();
            return evaluation;
          }
          if( appended_exact_constraint )
            trial_curve.mass_fraction_constraints.pop_back();

          if( conditional.m_status == RelActAutoSolution::Status::UserCanceled )
          {
            evaluation.status = EvalStatus::Canceled;
            evaluation.diagnostic = "Mass-fraction profiling was canceled.";
            return evaluation;
          }
          if( !RelActAutoSolution::is_usable_status(conditional.m_status)
              || !std::isfinite(conditional.m_chi2) )
          {
            evaluation.status = EvalStatus::Failed;
            evaluation.diagnostic = conditional.m_error_message.empty()
                ? "A conditional profile solve did not produce a usable finite point."
                : conditional.m_error_message;
            return evaluation;
          }
          const bool same_frozen_objective = matches_frozen_profile_objective(conditional);
          assert( same_frozen_objective );
          if( !same_frozen_objective )
          {
            evaluation.status = EvalStatus::Failed;
            evaluation.diagnostic = "A conditional profile solve changed the frozen continuum,"
                                    " branching-ratio, gamma-membership, or selected-model policy.";
            return evaluation;
          }

          try
          {
            // Conditional solves disable profiling, so the structured nominal accessor cannot
            // recurse.  Avoid the compatibility wrapper here: it also computes/floors a legacy
            // Gaussian uncertainty and can mutate its historical reliability flag even though
            // this equality check needs only the reported fraction.
            actual_reported = conditional.mass_enrichment_result(target,curve_index).fraction;
          }catch( const std::exception &e )
          {
            evaluation.status = EvalStatus::Failed;
            evaluation.diagnostic = std::string("Could not read a conditional mass fraction: ")
                                    + e.what();
            return evaluation;
          }
          if( !std::isfinite(actual_reported) )
          {
            evaluation.status = EvalStatus::Failed;
            evaluation.diagnostic = "A conditional solve returned a non-finite reported fraction.";
            return evaluation;
          }

          last_reported_change = std::isfinite(previous_reported)
                               ? std::fabs(actual_reported-previous_reported)
                               : std::numeric_limits<double>::infinity();
          previous_reported = actual_reported;

          if( !use_reported_fraction_constraint
              || (std::fabs(actual_reported-requested_for_solve) <= equality_tolerance) )
          {
            equality_satisfied = true;
            break;
          }

          // Standard augmented-Lagrangian multiplier update, retaining at most three passes for a
          // difficult point while charging every pass against the global conditional-fit cap.
          const double violation = actual_reported - requested_for_solve;
          lagrange_multiplier += penalty*violation;
          penalty *= 100.0;
          semantic_warm = std::make_unique<RelActAutoSolution>(conditional);
        }//for( augmented-Lagrangian outer passes )

        bool reached_structural_bound = false;
        if( !equality_satisfied && use_reported_fraction_constraint )
        {
          // A post-correlation fraction need not span the full mathematical [0,1] interval.  For
          // example, ByPu239Only gives some corrected isotopes a structural maximum well inside
          // that interval.  Successive augmented-Lagrangian penalties which converge to the same
          // independently evaluated point between the baseline and request identify that limit,
          // even when geometric bracketing encounters it before requesting zero or one.  Requiring
          // both a genuine unresolved gap and a stable repeated point distinguishes this from an
          // ordinary equality solve which merely needed another iteration.
          const double directional_span
              = (direction == RelActCalcAutoImp::ProfileLikelihood::Direction::Lower)
                ? (control_baseline-control_lower) : (control_upper-control_baseline);
          reached_structural_bound
              = RelActCalcAutoImp::ProfileLikelihood::stable_reported_structural_bound(
                    requested_for_solve,actual_reported,control_baseline,last_reported_change,
                    equality_tolerance,(std::max)(0.0,directional_span),
                    evaluation.num_fits,direction );
        }

        if( !equality_satisfied && !reached_structural_bound )
        {
          evaluation.status = (evaluation.num_fits >= remaining_fit_budget)
              ? EvalStatus::FitCapReached : EvalStatus::Failed;
          evaluation.diagnostic = "The reported-coordinate equality did not converge within the"
              " augmented-Lagrangian/profile fit budget (requested="
              + SpecUtils::printCompact(requested_for_solve,12) + ", actual="
              + SpecUtils::printCompact(actual_reported,12) + ", tolerance="
              + SpecUtils::printCompact(equality_tolerance,12) + ", last change="
              + SpecUtils::printCompact(last_reported_change,12) + ").";
          return evaluation;
        }

        try
        {
          // Always expose the independently evaluated reported coordinate.  In particular, a
          // mathematical zero/one endpoint uses a scale-aware inward limit and must not label that
          // objective as though it had been evaluated at another coordinate.
          evaluation.reported_fraction = actual_reported;
          evaluation.reached_feasible_bound = reached_structural_bound;
          evaluation.control_tolerance = use_reported_fraction_constraint
                                       ? equality_tolerance : 0.0;
          evaluation.delta_chi2 = conditional.m_chi2 - solution.m_chi2;

          // Preserve this successful conditional point and its source-level starts.  A later
          // request in this direction selects the cached point nearest in profile coordinate.
          // The target itself remains controlled by the equality; sibling activities and all
          // fitted ages are safely transferable.
          for( size_t ci = 0; ci < trial_options.rel_eff_curves.size()
                              && ci < conditional.m_rel_activities.size(); ++ci )
          {
            for( NucInputInfo &input_nuc : trial_options.rel_eff_curves[ci].nuclides )
            {
              const RelEffCurveInput &input_curve = trial_options.rel_eff_curves[ci];
              const bool controlled_by_activity_ratio = std::any_of(
                  begin(input_curve.act_ratio_constraints),end(input_curve.act_ratio_constraints),
                  [&input_nuc]( const RelEffCurveInput::ActRatioConstraint &constraint ){
                    return constraint.constrained_source == input_nuc.source;
                  } );
              const auto fitted = std::find_if( begin(conditional.m_rel_activities[ci]),
                                                end(conditional.m_rel_activities[ci]),
                    [&input_nuc]( const NuclideRelAct &act ){
                      return act.source == input_nuc.source;
                    } );
              if( fitted == end(conditional.m_rel_activities[ci]) )
                continue;
              if( std::isfinite(fitted->age) && (fitted->age >= 0.0) )
                input_nuc.age = fitted->age;
              const bool exact_profile_target = !use_reported_fraction_constraint
                                             && (ci == curve_index)
                                             && (input_nuc.source == target_source);
              if( !exact_profile_target && !input_nuc.starting_rel_act
                  && !controlled_by_activity_ratio
                  && std::isfinite(fitted->rel_activity)
                  && (fitted->rel_activity > 0.0) )
                  input_nuc.starting_rel_act = fitted->rel_activity;
            }
          }

          if( evaluation.delta_chi2 < -baseline_improvement_tolerance )
          {
            better_conditional = std::make_unique<RelActAutoSolution>(std::move(conditional));
            evaluation.status = EvalStatus::BetterBaseline;
            evaluation.diagnostic = "A conditional fit improved the frozen physical objective by "
                + std::to_string(-evaluation.delta_chi2) + ".";
            return evaluation;
          }

          evaluation.delta_chi2 = (std::max)(0.0,evaluation.delta_chi2);
          evaluation.status = EvalStatus::Success;
          warm_controls.push_back(evaluation.reported_fraction);
          warm_points.push_back( {std::move(trial_options),std::move(conditional)} );
          return evaluation;
        }catch( const std::exception &e )
        {
          evaluation.status = EvalStatus::Failed;
          evaluation.diagnostic = std::string("Could not read a conditional mass fraction: ")
                                  + e.what();
          return evaluation;
        }
      };

      const std::array<double,2> thresholds{{
        solution.m_cov_scale,
        solution.m_cov_scale * 3.841458820694124
      }};
      const RelActCalcAutoImp::ProfileLikelihood::ScanResult scan
          = RelActCalcAutoImp::ProfileLikelihood::scan(
              control_baseline,nominal.fraction,control_lower,control_upper,
              thresholds,max_conditional_fits,evaluate );
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

  // A profile-discovered point is a deterministic-search trigger.  Remove the conditional
  // equality, polish the complete named candidate matrix on the same frozen ROIs, and then
  // restart every requested profile so covariance and thresholds use the selected baseline.
  Options restart_options = solution.m_options;
  restart_options.profile_only_mass_fraction_constraint.reset();
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
        &frozen_profile_model_policy,true );
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
