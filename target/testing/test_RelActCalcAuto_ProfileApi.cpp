/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; either version 2.1 of the License, or (at your option)
 any later version.
*/

#include "InterSpec_config.h"

#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <type_traits>

#define BOOST_TEST_MODULE RelActCalcAuto_ProfileApi_suite
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

namespace
{
void initialize_decay_database()
{
  static bool initialized = false;
  if( initialized )
    return;

  string data_dir;
  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      data_dir = arg.substr( 10 );
  }
  SpecUtils::ireplace_all( data_dir, "%20", " " );

  if( data_dir.empty() )
  {
    for( const char *candidate : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(candidate, "sandia.decay.xml") ) )
      {
        data_dir = candidate;
        break;
      }
    }
  }

  BOOST_REQUIRE_MESSAGE( !data_dir.empty(), "Could not locate sandia.decay.xml" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory(data_dir) );
  BOOST_REQUIRE( DecayDataBaseServer::database() );
  initialized = true;
}

int xml_version( const rapidxml::xml_node<char> *node )
{
  BOOST_REQUIRE( node );
  const rapidxml::xml_attribute<char> * const attr = node->first_attribute( "version" );
  BOOST_REQUIRE( attr );
  return std::stoi( string(attr->value(), attr->value_size()) );
}

rapidxml::xml_node<char> *new_root( rapidxml::xml_document<char> &doc )
{
  rapidxml::xml_node<char> * const root
      = doc.allocate_node( rapidxml::node_element, "Root" );
  doc.append_node( root );
  return root;
}

RelActCalcAuto::NucInputInfo pu239_input( const bool force_profile )
{
  initialize_decay_database();
  const SandiaDecay::Nuclide * const pu239 = DecayDataBaseServer::database()->nuclide( "Pu239" );
  BOOST_REQUIRE( pu239 );

  RelActCalcAuto::NucInputInfo input;
  input.source = pu239;
  input.age = 0.0;
  input.fit_age = false;
  input.force_profile_mass_fraction = force_profile;
  return input;
}
}


BOOST_AUTO_TEST_CASE( status_values_are_backward_compatible )
{
  using Solution = RelActCalcAuto::RelActAutoSolution;
  using Status = Solution::Status;

  BOOST_CHECK_EQUAL( static_cast<int>(Status::Success), 0 );
  BOOST_CHECK_EQUAL( static_cast<int>(Status::NotInitiated), 1 );
  BOOST_CHECK_EQUAL( static_cast<int>(Status::FailedToSetupProblem), 2 );
  BOOST_CHECK_EQUAL( static_cast<int>(Status::FailToSolveProblem), 3 );
  BOOST_CHECK_EQUAL( static_cast<int>(Status::UserCanceled), 4 );
  BOOST_CHECK_EQUAL( static_cast<int>(Status::UsableWithWarnings), 5 );

  BOOST_CHECK( Solution::is_usable_status(Status::Success) );
  BOOST_CHECK( Solution::is_usable_status(Status::UsableWithWarnings) );
  BOOST_CHECK( !Solution::is_usable_status(Status::NotInitiated) );
  BOOST_CHECK( !Solution::is_usable_status(Status::FailedToSetupProblem) );
  BOOST_CHECK( !Solution::is_usable_status(Status::FailToSolveProblem) );
  BOOST_CHECK( !Solution::is_usable_status(Status::UserCanceled) );
}


BOOST_AUTO_TEST_CASE( pu242_correlation_rejects_unrepresentable_pu_isotopes_at_validation )
{
  initialize_decay_database();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const pu236 = db->nuclide("Pu236");
  const SandiaDecay::Nuclide * const pu239 = db->nuclide("Pu239");
  const SandiaDecay::Nuclide * const pu240 = db->nuclide("Pu240");
  BOOST_REQUIRE( pu236 && pu239 && pu240 );

  const auto input_for = []( const SandiaDecay::Nuclide * const nuc ) {
    RelActCalcAuto::NucInputInfo input;
    input.source = nuc;
    input.age = 0.0;
    input.fit_age = false;
    return input;
  };

  RelActCalcAuto::RelEffCurveInput supported;
  supported.pu242_correlation_method = RelActCalc::PuCorrMethod::ByPu239Only;
  supported.nuclides = { input_for(pu239),input_for(pu240) };
  BOOST_CHECK_NO_THROW( supported.check_nuclide_constraints() );

  RelActCalcAuto::RelEffCurveInput unsupported = supported;
  unsupported.nuclides.push_back( input_for(pu236) );
  BOOST_CHECK_EXCEPTION( unsupported.check_nuclide_constraints(),std::logic_error,
    []( const std::logic_error &error ) {
      return SpecUtils::icontains(error.what(),"Pu236")
          && SpecUtils::icontains(error.what(),"cannot be used");
    } );

  // Without a correlation the exact same isotope list remains a valid model input.
  unsupported.pu242_correlation_method = RelActCalc::PuCorrMethod::NotApplicable;
  BOOST_CHECK_NO_THROW( unsupported.check_nuclide_constraints() );
}


BOOST_AUTO_TEST_CASE( nuc_profile_request_uses_lowest_compatible_xml_version )
{
  using RelActCalcAuto::NucInputInfo;

  // The default request remains byte/schema compatible with legacy v0 readers.
  NucInputInfo legacy_compatible = pu239_input( false );
  rapidxml::xml_document<char> default_doc;
  rapidxml::xml_node<char> * const default_root = new_root( default_doc );
  legacy_compatible.toXml( default_root );
  const rapidxml::xml_node<char> * const default_node
      = default_root->first_node( "NucInputInfo" );
  BOOST_REQUIRE( default_node );
  BOOST_CHECK_EQUAL( xml_version(default_node), 0 );
  BOOST_CHECK( !default_node->first_node("ForceMassFractionProfile") );

  // Explicit profiling is the only NucInputInfo feature that requires v1.
  NucInputInfo forced = pu239_input( true );
  rapidxml::xml_document<char> forced_doc;
  rapidxml::xml_node<char> * const forced_root = new_root( forced_doc );
  forced.toXml( forced_root );
  const rapidxml::xml_node<char> * const forced_node
      = forced_root->first_node( "NucInputInfo" );
  BOOST_REQUIRE( forced_node );
  BOOST_CHECK_EQUAL( xml_version(forced_node), 1 );
  BOOST_REQUIRE( forced_node->first_node("ForceMassFractionProfile") );
  BOOST_CHECK_EQUAL( string(forced_node->first_node("ForceMassFractionProfile")->value()), "true" );

  NucInputInfo restored;
  BOOST_REQUIRE_NO_THROW( restored.fromXml(forced_node) );
  BOOST_CHECK( restored.force_profile_mass_fraction );
  BOOST_CHECK( restored == forced );

  // Reading a legacy record must actively restore the safe default, even when the
  // destination object previously held a forced request.
  restored.force_profile_mass_fraction = true;
  BOOST_REQUIRE_NO_THROW( restored.fromXml(default_node) );
  BOOST_CHECK( !restored.force_profile_mass_fraction );
  BOOST_CHECK( restored == legacy_compatible );
}


BOOST_AUTO_TEST_CASE( options_profile_behavior_round_trips_and_preserves_old_versions )
{
  using RelActCalcAuto::Options;

  // Default-on behavior requires no new element and therefore retains Options v2.
  Options defaults;
  BOOST_CHECK( defaults.auto_profile_weak_mass_fractions );
  rapidxml::xml_document<char> default_doc;
  rapidxml::xml_node<char> * const default_root = new_root( default_doc );
  const rapidxml::xml_node<char> * const default_node = defaults.toXml( default_root );
  BOOST_REQUIRE( default_node );
  BOOST_CHECK_EQUAL( xml_version(default_node), 2 );
  BOOST_CHECK( !default_node->first_node("AutoProfileWeakMassFractions") );

  Options restored_defaults;
  restored_defaults.auto_profile_weak_mass_fractions = false;
  BOOST_REQUIRE_NO_THROW( restored_defaults.fromXml(default_node) );
  BOOST_CHECK( restored_defaults.auto_profile_weak_mass_fractions );

  // Disabling automatic profiling is non-default behavior and is persisted as v4.
  Options disabled;
  disabled.auto_profile_weak_mass_fractions = false;
  rapidxml::xml_document<char> disabled_doc;
  rapidxml::xml_node<char> * const disabled_root = new_root( disabled_doc );
  const rapidxml::xml_node<char> * const disabled_node = disabled.toXml( disabled_root );
  BOOST_REQUIRE( disabled_node );
  BOOST_CHECK_EQUAL( xml_version(disabled_node), 4 );
  BOOST_REQUIRE( disabled_node->first_node("AutoProfileWeakMassFractions") );
  BOOST_CHECK_EQUAL( string(disabled_node->first_node("AutoProfileWeakMassFractions")->value()), "false" );

  Options restored_disabled;
  BOOST_REQUIRE_NO_THROW( restored_disabled.fromXml(disabled_node) );
  BOOST_CHECK( !restored_disabled.auto_profile_weak_mass_fractions );

  // The preceding v3 metadata feature must still select v3 when no v4 feature is used.
  Options with_filename;
  with_filename.foreground_filename = "example.n42";
  rapidxml::xml_document<char> v3_doc;
  const rapidxml::xml_node<char> * const v3_node = with_filename.toXml( new_root(v3_doc) );
  BOOST_REQUIRE( v3_node );
  BOOST_CHECK_EQUAL( xml_version(v3_node), 3 );
  BOOST_CHECK( !v3_node->first_node("AutoProfileWeakMassFractions") );
}


BOOST_AUTO_TEST_CASE( robust_solve_round_trips_and_does_not_disturb_older_versions )
{
  using RelActCalcAuto::Options;

  // Off is the historical behavior, so it must add no element and must not bump the version - an
  // existing config that never asks for a robust solve has to stay readable by older builds.
  Options defaults;
  BOOST_CHECK( !defaults.robust_solve );
  rapidxml::xml_document<char> default_doc;
  const rapidxml::xml_node<char> * const default_node = defaults.toXml( new_root(default_doc) );
  BOOST_REQUIRE( default_node );
  BOOST_CHECK_EQUAL( xml_version(default_node), 2 );
  BOOST_CHECK( !default_node->first_node("RobustSolve") );

  Options restored_defaults;
  restored_defaults.robust_solve = true;
  BOOST_REQUIRE_NO_THROW( restored_defaults.fromXml(default_node) );
  BOOST_CHECK( !restored_defaults.robust_solve );

  // Asking for a robust solve is the only thing that selects v5.
  Options robust;
  robust.robust_solve = true;
  rapidxml::xml_document<char> robust_doc;
  const rapidxml::xml_node<char> * const robust_node = robust.toXml( new_root(robust_doc) );
  BOOST_REQUIRE( robust_node );
  BOOST_CHECK_EQUAL( xml_version(robust_node), 5 );
  BOOST_REQUIRE( robust_node->first_node("RobustSolve") );
  BOOST_CHECK_EQUAL( string(robust_node->first_node("RobustSolve")->value()), "true" );

  Options restored_robust;
  BOOST_REQUIRE_NO_THROW( restored_robust.fromXml(robust_node) );
  BOOST_CHECK( restored_robust.robust_solve );
  BOOST_CHECK( restored_robust == robust );

  // v5 dominates the older conditional fields rather than being ignored by them.
  Options robust_and_disabled;
  robust_and_disabled.robust_solve = true;
  robust_and_disabled.auto_profile_weak_mass_fractions = false;
  rapidxml::xml_document<char> both_doc;
  const rapidxml::xml_node<char> * const both_node
      = robust_and_disabled.toXml( new_root(both_doc) );
  BOOST_REQUIRE( both_node );
  BOOST_CHECK_EQUAL( xml_version(both_node), 5 );
  Options restored_both;
  BOOST_REQUIRE_NO_THROW( restored_both.fromXml(both_node) );
  BOOST_CHECK( restored_both.robust_solve );
  BOOST_CHECK( !restored_both.auto_profile_weak_mass_fractions );

  // The two flags are independent: automatic profiling still persists at v4 on its own.
  Options disabled_only;
  disabled_only.auto_profile_weak_mass_fractions = false;
  rapidxml::xml_document<char> disabled_doc;
  const rapidxml::xml_node<char> * const disabled_node
      = disabled_only.toXml( new_root(disabled_doc) );
  BOOST_REQUIRE( disabled_node );
  BOOST_CHECK_EQUAL( xml_version(disabled_node), 4 );
  BOOST_CHECK( !disabled_node->first_node("RobustSolve") );

  // Differing only by this flag must not compare equal, or a state change could be lost.
  BOOST_CHECK( !(robust == defaults) );
}


BOOST_AUTO_TEST_CASE( forced_nuclide_promotes_containing_options_to_v4 )
{
  RelActCalcAuto::Options options;
  RelActCalcAuto::RelEffCurveInput curve;
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnY;
  curve.rel_eff_eqn_order = 2;
  curve.nuclides.push_back( pu239_input(true) );
  options.rel_eff_curves.push_back( curve );

  rapidxml::xml_document<char> doc;
  const rapidxml::xml_node<char> * const options_node = options.toXml( new_root(doc) );
  BOOST_REQUIRE( options_node );
  BOOST_CHECK_EQUAL( xml_version(options_node), 4 );
  BOOST_CHECK( !options_node->first_node("AutoProfileWeakMassFractions") );

  const rapidxml::xml_node<char> * const curves
      = options_node->first_node( "RelEffCurveInputs" );
  BOOST_REQUIRE( curves );
  const rapidxml::xml_node<char> * const curve_node
      = curves->first_node( "RelEffCurveInput" );
  BOOST_REQUIRE( curve_node );
  const rapidxml::xml_node<char> * const nucs
      = curve_node->first_node( "NucInputInfoList" );
  BOOST_REQUIRE( nucs );
  const rapidxml::xml_node<char> * const nuc_node
      = nucs->first_node( "NucInputInfo" );
  BOOST_REQUIRE( nuc_node );
  BOOST_CHECK_EQUAL( xml_version(nuc_node), 1 );

  RelActCalcAuto::Options restored;
  BOOST_REQUIRE_NO_THROW( restored.fromXml(options_node) );
  BOOST_REQUIRE_EQUAL( restored.rel_eff_curves.size(), 1u );
  BOOST_REQUIRE_EQUAL( restored.rel_eff_curves[0].nuclides.size(), 1u );
  BOOST_CHECK( restored.auto_profile_weak_mass_fractions );
  BOOST_CHECK( restored.rel_eff_curves[0].nuclides[0].force_profile_mass_fraction );
}


BOOST_AUTO_TEST_CASE( structured_mass_fraction_result_types_are_public_and_bounded )
{
  using Solution = RelActCalcAuto::RelActAutoSolution;

  // This compile-time check guards the public accessor signature independently of a
  // particular fitted spectrum.
  using Accessor = Solution::MassFractionResult (Solution::*)(
      const SandiaDecay::Nuclide *, size_t) const;
  Accessor accessor = &Solution::mass_enrichment_result;
  BOOST_CHECK( accessor != nullptr );

  Solution::MassFractionProfileInterval interval;
  interval.confidence_level = 0.6827;
  interval.delta_chi2 = 1.0;
  interval.lower = 0.12;
  interval.upper = 0.31;
  interval.lower_kind = Solution::MassFractionProfileEndpointKind::LikelihoodCrossing;
  interval.upper_kind = Solution::MassFractionProfileEndpointKind::InputConstraintLimit;

  Solution::MassFractionProfileResult profile;
  profile.status = Solution::MassFractionProfileStatus::BoundaryLimited;
  profile.reason = Solution::MassFractionProfileReason::Forced;
  profile.intervals.push_back( interval );
  profile.num_fits = 7;

  Solution::MassFractionResult result;
  result.fraction = 0.2;
  result.covariance_one_sigma = 0.04;
  result.covariance_quality = Solution::MassFractionCovarianceQuality::LocallyUnreliable;
  result.profile = profile;
  result.uncorrected_fraction = 0.19;

  BOOST_REQUIRE( result.profile );
  BOOST_REQUIRE_EQUAL( result.profile->intervals.size(), 1u );
  BOOST_CHECK( result.profile->status == Solution::MassFractionProfileStatus::BoundaryLimited );
  BOOST_CHECK( result.profile->reason == Solution::MassFractionProfileReason::Forced );
  BOOST_CHECK_GE( result.profile->intervals[0].lower, 0.0 );
  BOOST_CHECK_LE( result.profile->intervals[0].upper, 1.0 );
  BOOST_CHECK_LT( result.profile->intervals[0].lower, result.profile->intervals[0].upper );
}


BOOST_AUTO_TEST_CASE( forced_single_isotope_profile_contract_is_zero_fit_physical_boundary )
{
  using Solution = RelActCalcAuto::RelActAutoSolution;
  Solution::MassFractionProfileInterval p68;
  p68.confidence_level = 0.6827;
  p68.delta_chi2 = 1.0;
  p68.lower = p68.upper = 1.0;
  p68.lower_kind = p68.upper_kind
      = Solution::MassFractionProfileEndpointKind::PhysicalLimit;
  Solution::MassFractionProfileInterval p95 = p68;
  p95.confidence_level = 0.95;
  p95.delta_chi2 = 3.841458820694124;

  Solution::MassFractionProfileResult profile;
  profile.status = Solution::MassFractionProfileStatus::BoundaryLimited;
  profile.reason = Solution::MassFractionProfileReason::Forced;
  profile.num_fits = 0;
  profile.intervals = {p68,p95};

  BOOST_REQUIRE_EQUAL( profile.intervals.size(),2U );
  BOOST_CHECK( profile.status == Solution::MassFractionProfileStatus::BoundaryLimited );
  BOOST_CHECK( profile.reason == Solution::MassFractionProfileReason::Forced );
  BOOST_CHECK_EQUAL( profile.num_fits,0U );
  for( const Solution::MassFractionProfileInterval &interval : profile.intervals )
  {
    BOOST_CHECK_EQUAL( interval.lower,1.0 );
    BOOST_CHECK_EQUAL( interval.upper,1.0 );
    BOOST_CHECK( interval.lower_kind
                 == Solution::MassFractionProfileEndpointKind::PhysicalLimit );
    BOOST_CHECK( interval.upper_kind
                 == Solution::MassFractionProfileEndpointKind::PhysicalLimit );
  }
}


/** The §3.2 rule: "profile the mass fraction" is only meaningful where there ARE isotopics.

 A fit of Cs-137 plus I-131 has none - each element's normalized fraction is identically one - so a
 mass-fraction interval there would be an artifact of the reporting code rather than a statement
 about the data, and the quantity carrying the information is the activity.
 */
BOOST_AUTO_TEST_CASE( profilable_kind_selection_follows_what_carries_information )
{
  using Kind = RelActCalcAuto::Options::ProfileTarget::Kind;

  initialize_decay_database();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const cs137 = db->nuclide("Cs137");
  const SandiaDecay::Nuclide * const i131  = db->nuclide("I131");
  const SandiaDecay::Nuclide * const u235  = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238  = db->nuclide("U238");
  const SandiaDecay::Nuclide * const pu239 = db->nuclide("Pu239");
  const SandiaDecay::Nuclide * const pu240 = db->nuclide("Pu240");
  const SandiaDecay::Nuclide * const pu242 = db->nuclide("Pu242");
  const SandiaDecay::Nuclide * const co60  = db->nuclide("Co60");
  BOOST_REQUIRE( cs137 && i131 && u235 && u238 && pu239 && pu240 && pu242 && co60 );

  const auto input_for = []( const SandiaDecay::Nuclide * const nuc, const bool fit_age ) {
    RelActCalcAuto::NucInputInfo input;
    input.source = nuc;
    input.age = 0.0;
    input.fit_age = fit_age;
    return input;
  };
  const auto options_with = []( const RelActCalcAuto::RelEffCurveInput &curve ) {
    RelActCalcAuto::Options options;
    options.rel_eff_curves = { curve };
    return options;
  };
  const auto kinds_for = []( const RelActCalcAuto::Options &options,
                             const SandiaDecay::Nuclide * const nuc ) {
    return RelActCalcAuto::profilable_quantity_kinds( options,
                                                      RelActCalcAuto::SrcVariant(nuc),0 );
  };

  // Single-isotope elements: the activity is the only thing that can be profiled, and a mass
  // fraction must NEVER be offered.
  RelActCalcAuto::RelEffCurveInput singles;
  singles.nuclides = { input_for(cs137,false),input_for(i131,false) };
  const RelActCalcAuto::Options single_isotope_options = options_with( singles );
  const vector<Kind> cs137_kinds = kinds_for( single_isotope_options,cs137 );
  const vector<Kind> i131_kinds  = kinds_for( single_isotope_options,i131 );
  BOOST_REQUIRE_EQUAL( cs137_kinds.size(),1U );
  BOOST_CHECK( cs137_kinds[0] == Kind::RelativeActivity );
  BOOST_REQUIRE_EQUAL( i131_kinds.size(),1U );
  BOOST_CHECK( i131_kinds[0] == Kind::RelativeActivity );
  BOOST_CHECK( std::find(begin(cs137_kinds),end(cs137_kinds),Kind::MassFraction)
               == end(cs137_kinds) );

  // Two isotopes of one element: the isotopic fraction is the informative quantity.
  RelActCalcAuto::RelEffCurveInput uranium;
  uranium.nuclides = { input_for(u235,false),input_for(u238,false) };
  const RelActCalcAuto::Options uranium_options = options_with( uranium );
  const vector<Kind> u235_kinds = kinds_for( uranium_options,u235 );
  BOOST_REQUIRE_EQUAL( u235_kinds.size(),1U );
  BOOST_CHECK( u235_kinds[0] == Kind::MassFraction );

  // Pu-242 by correlation is a reported quantity that is deliberately NOT a fitted input, so it is
  // never offered as a profile target: it has no free parameter and hence no likelihood direction
  // of its own, and an interval over it would describe the correlation's systematic rather than
  // anything this data constrains.  Its *reporting* is untouched - the correlation still supplies
  // its value and still renormalizes its siblings' fractions - which is why the sibling Pu-239
  // below is still profilable.
  RelActCalcAuto::RelEffCurveInput plutonium;
  plutonium.pu242_correlation_method = RelActCalc::PuCorrMethod::ByPu239Only;
  plutonium.nuclides = { input_for(pu239,false),input_for(pu240,false) };
  const RelActCalcAuto::Options plutonium_options = options_with( plutonium );
  const vector<Kind> pu242_kinds = kinds_for( plutonium_options,pu242 );
  BOOST_CHECK( pu242_kinds.empty() );
  BOOST_CHECK( kinds_for(plutonium_options,pu239).size() == 1U );

  // Without the correlation the same Pu-242 is not part of the model at all - which is now simply
  // the general rule for a source that is not a fitted input, rather than a special case.
  RelActCalcAuto::RelEffCurveInput uncorrelated = plutonium;
  uncorrelated.pu242_correlation_method = RelActCalc::PuCorrMethod::NotApplicable;
  BOOST_CHECK( kinds_for(options_with(uncorrelated),pu242).empty() );

  // A fitted age is an independent second quantity, not an alternative to the first.
  RelActCalcAuto::RelEffCurveInput aged_uranium;
  aged_uranium.nuclides = { input_for(u235,true),input_for(u238,false) };
  const vector<Kind> aged_kinds = kinds_for( options_with(aged_uranium),u235 );
  BOOST_REQUIRE_EQUAL( aged_kinds.size(),2U );
  BOOST_CHECK( aged_kinds[0] == Kind::MassFraction );
  BOOST_CHECK( aged_kinds[1] == Kind::Age );

  // A lone isotope with a fitted age offers both of the quantities that mean something.
  RelActCalcAuto::RelEffCurveInput aged_single;
  aged_single.nuclides = { input_for(cs137,true),input_for(i131,false) };
  const vector<Kind> aged_single_kinds = kinds_for( options_with(aged_single),cs137 );
  BOOST_REQUIRE_EQUAL( aged_single_kinds.size(),2U );
  BOOST_CHECK( aged_single_kinds[0] == Kind::RelativeActivity );
  BOOST_CHECK( aged_single_kinds[1] == Kind::Age );

  // A source that is not on the curve has nothing to profile.
  BOOST_CHECK( kinds_for(uranium_options,co60).empty() );

  // A ratio needs a denominator only the caller can pick, so it is never auto-offered even though
  // it is the one gauge-invariant kind.
  for( const vector<Kind> &kinds : { cs137_kinds,u235_kinds,aged_kinds } )
    BOOST_CHECK( std::find(begin(kinds),end(kinds),Kind::ActivityRatio) == end(kinds) );
}


/** Invariant 9 of the profile design: a profile target's kind and source must be coherent,
 checked before any solve rather than producing a silently meaningless interval. */
BOOST_AUTO_TEST_CASE( profile_target_coherence_is_checked_up_front )
{
  using Target = RelActCalcAuto::Options::ProfileTarget;
  using Kind = Target::Kind;

  initialize_decay_database();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  const SandiaDecay::Nuclide * const co60 = db->nuclide("Co60");
  BOOST_REQUIRE( u235 && u238 && co60 );

  // Every kind is nameable, and the names are distinct - reports and diagnostics quote them.
  const Kind all_kinds[] = { Kind::MassFraction,Kind::RelativeActivity,
                             Kind::ActivityRatio,Kind::Age };
  set<string> names;
  for( const Kind kind : all_kinds )
  {
    const string name = Target::to_str(kind);
    BOOST_CHECK( !name.empty() );
    names.insert( name );
  }
  BOOST_CHECK_EQUAL( names.size(),4U );

  RelActCalcAuto::NucInputInfo u235_input;
  u235_input.source = u235;
  u235_input.age = 0.0;
  u235_input.fit_age = false;
  RelActCalcAuto::NucInputInfo u238_input = u235_input;
  u238_input.source = u238;

  RelActCalcAuto::RelEffCurveInput curve;
  curve.nuclides = { u235_input,u238_input };
  RelActCalcAuto::Options options;
  options.rel_eff_curves = { curve };

  Target good;
  good.kind = Kind::MassFraction;
  good.source = u235;
  good.rel_eff_curve_index = 0;
  BOOST_CHECK_EQUAL( good.why_not_usable(options),string() );

  const auto rejects = [&options]( const Target &target,const char * const because ) {
    const string why = target.why_not_usable( options );
    BOOST_CHECK_MESSAGE( !why.empty(),
        string("expected rejection because ") + because + ", but the target was accepted" );
  };

  Target absent = good;       absent.source = co60;
  rejects( absent,"Co60 is not a source on the curve" );

  Target bad_curve = good;    bad_curve.rel_eff_curve_index = 7;
  rejects( bad_curve,"curve 7 does not exist" );

  Target fixed_age = good;    fixed_age.kind = Kind::Age;
  rejects( fixed_age,"U235's age is not being fitted" );

  Target ratio = good;
  ratio.kind = Kind::ActivityRatio;
  ratio.denominator = u235;      //same as the numerator
  rejects( ratio,"a source divided by itself carries no information" );

  ratio.denominator = co60;
  rejects( ratio,"the denominator is not a source on its curve" );

  ratio.denominator = u238;
  BOOST_CHECK_EQUAL( ratio.why_not_usable(options),string() );

  // An age target becomes usable exactly when the age becomes a fitted parameter.  `RelEffCurveInput`
  // shares one age per element by default, so ask for independent ages to isolate the per-source flag.
  RelActCalcAuto::Options fitted_age_options = options;
  fitted_age_options.rel_eff_curves[0].nucs_of_el_same_age = false;
  fitted_age_options.rel_eff_curves[0].nuclides[0].fit_age = true;
  Target age = good;
  age.kind = Kind::Age;
  BOOST_CHECK_EQUAL( age.why_not_usable(fitted_age_options),string() );

  // A single fitted age shared across the element may be named through any of its isotopes: the
  // functor resolves to the one parameter, so constraining U238's age constrains U235's too.
  RelActCalcAuto::Options shared_age_options = fitted_age_options;
  shared_age_options.rel_eff_curves[0].nucs_of_el_same_age = true;
  Target follower_age = age;
  follower_age.source = u238;                   //fit_age false on this input
  BOOST_CHECK_EQUAL( follower_age.why_not_usable(shared_age_options),string() );
  // ...but with independent ages, U238's own flag is what decides.
  BOOST_CHECK( !follower_age.why_not_usable(fitted_age_options).empty() );
}


/** Constraint shapes that leave no slot scanning the reported quantity exactly are rejected
 UP FRONT by `why_not_usable`, mirroring the scan-time carrier-install refusals, so callers (and
 the GUI) can warn before paying for a solve. */
BOOST_AUTO_TEST_CASE( constraint_shapes_are_rejected_up_front )
{
  using Target = RelActCalcAuto::Options::ProfileTarget;
  using Kind = Target::Kind;

  initialize_decay_database();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const u234 = db->nuclide("U234");
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  const SandiaDecay::Nuclide * const cs137 = db->nuclide("Cs137");
  BOOST_REQUIRE( u234 && u235 && u238 && cs137 );

  const auto input_for = []( const SandiaDecay::Nuclide * const nuc ) {
    RelActCalcAuto::NucInputInfo input;
    input.source = nuc;
    input.age = 0.0;
    input.fit_age = false;
    return input;
  };

  RelActCalcAuto::RelEffCurveInput curve;
  curve.nuclides = { input_for(u234),input_for(u235),input_for(u238),input_for(cs137) };
  RelActCalcAuto::Options base_options;
  base_options.rel_eff_curves = { curve };

  Target u238_frac;
  u238_frac.kind = Kind::MassFraction;
  u238_frac.source = u238;
  BOOST_REQUIRE_EQUAL( u238_frac.why_not_usable(base_options),string() );

  const auto with_window = [&base_options,&input_for]( const SandiaDecay::Nuclide * const nuc,
                                                       const double lower,const double upper ) {
    RelActCalcAuto::Options options = base_options;
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = nuc;
    constraint.lower_mass_fraction = lower;
    constraint.upper_mass_fraction = upper;
    options.rel_eff_curves[0].mass_fraction_constraints.push_back( constraint );
    return options;
  };

  // A windowed SIBLING makes the target unprofilable (the shared block's carrier pins the
  // constrained total, not the target's fraction)...
  BOOST_CHECK( !u238_frac.why_not_usable( with_window(u235,0.10,0.50) ).empty() );
  // ...while the target's own window is the exactly-scannable sole-range chart.
  Target u235_frac = u238_frac;
  u235_frac.source = u235;
  BOOST_CHECK_EQUAL( u235_frac.why_not_usable( with_window(u235,0.10,0.50) ),string() );
  // A fixed target has no direction at all.
  BOOST_CHECK( !u235_frac.why_not_usable( with_window(u235,0.30,0.30) ).empty() );
  // An unconstrained target beside a FIXED sibling still has no exactly-scannable slot.
  BOOST_CHECK( !u238_frac.why_not_usable( with_window(u235,0.30,0.30) ).empty() );

  // Activity bounds / start values sit on the very slot the reparameterization would swap.
  {
    RelActCalcAuto::Options options = base_options;
    options.rel_eff_curves[0].nuclides[2].min_rel_act = 0.0;
    BOOST_CHECK( !u238_frac.why_not_usable(options).empty() );
  }

  // A ratio-constrained target's slot is inert for a mass fraction...
  RelActCalcAuto::Options ratio_options = base_options;
  {
    RelActCalcAuto::RelEffCurveInput::ActRatioConstraint link;
    link.controlling_source = u235;
    link.constrained_source = u238;
    link.constrained_to_controlled_activity_ratio = 2.0;
    ratio_options.rel_eff_curves[0].act_ratio_constraints.push_back( link );
  }
  BOOST_CHECK( !u238_frac.why_not_usable(ratio_options).empty() );
  // ...and the controlling same-element isotope is circular for a mass fraction too...
  BOOST_CHECK( !u235_frac.why_not_usable(ratio_options).empty() );
  // ...but a RELATIVE ACTIVITY of the constrained source stays profilable: the scan pins the
  // chain root's slot, in which the reported activity is exactly linear.
  Target u238_act = u238_frac;
  u238_act.kind = Kind::RelativeActivity;
  BOOST_CHECK_EQUAL( u238_act.why_not_usable(ratio_options),string() );

  // A mass-fraction-constrained source's activity is not slot-linear.
  Target u235_act = u238_act;
  u235_act.source = u235;
  BOOST_CHECK( !u235_act.why_not_usable( with_window(u235,0.10,0.50) ).empty() );

  // Ratios: cross-curve is not gauge invariant; rigid links make the ratio a model constant.
  Target cs_over_u235;
  cs_over_u235.kind = Kind::ActivityRatio;
  cs_over_u235.source = cs137;
  cs_over_u235.denominator = u235;
  BOOST_CHECK_EQUAL( cs_over_u235.why_not_usable(base_options),string() );
  Target cross_curve = cs_over_u235;
  cross_curve.denominator_curve_index = 1;
  {
    RelActCalcAuto::Options two_curves = base_options;
    two_curves.rel_eff_curves.push_back( base_options.rel_eff_curves[0] );
    BOOST_CHECK( !cross_curve.why_not_usable(two_curves).empty() );
  }
  Target linked_ratio;
  linked_ratio.kind = Kind::ActivityRatio;
  linked_ratio.source = u238;
  linked_ratio.denominator = u235;
  BOOST_CHECK( !linked_ratio.why_not_usable(ratio_options).empty() );
}


/** `Options::profile_targets` round-trips through XML, gates the v6 version bump, and feeds the
 options-level usability check. */
BOOST_AUTO_TEST_CASE( profile_targets_round_trip_and_gate_versioning )
{
  using RelActCalcAuto::Options;
  using Target = Options::ProfileTarget;
  using Kind = Target::Kind;

  initialize_decay_database();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE( u235 && u238 );

  // An empty target list adds no element and keeps the default (v2) serialization.
  Options defaults;
  rapidxml::xml_document<char> default_doc;
  const rapidxml::xml_node<char> * const default_node = defaults.toXml( new_root(default_doc) );
  BOOST_REQUIRE( default_node );
  BOOST_CHECK_EQUAL( xml_version(default_node), 2 );
  BOOST_CHECK( !default_node->first_node("ProfileTargets") );

  Options with_targets;
  Target frac;
  frac.kind = Kind::MassFraction;
  frac.source = u235;
  frac.rel_eff_curve_index = 0;
  Target ratio;
  ratio.kind = Kind::ActivityRatio;
  ratio.source = u238;
  ratio.denominator = u235;
  with_targets.profile_targets = { frac,ratio };

  rapidxml::xml_document<char> doc;
  const rapidxml::xml_node<char> * const node = with_targets.toXml( new_root(doc) );
  BOOST_REQUIRE( node );
  BOOST_CHECK_EQUAL( xml_version(node), 6 );
  BOOST_REQUIRE( node->first_node("ProfileTargets") );

  Options restored;
  BOOST_REQUIRE_NO_THROW( restored.fromXml(node) );
  BOOST_REQUIRE_EQUAL( restored.profile_targets.size(), 2U );
  BOOST_CHECK( restored.profile_targets[0] == frac );
  BOOST_CHECK( restored.profile_targets[1] == ratio );
  BOOST_CHECK( restored.profile_targets == with_targets.profile_targets );

  // An explicitly requested target that cannot be scanned fails the options check up front, with
  // the target's own reason.
  RelActCalcAuto::NucInputInfo u235_input;
  u235_input.source = u235;
  u235_input.age = 0.0;
  u235_input.fit_age = false;
  RelActCalcAuto::NucInputInfo u238_input = u235_input;
  u238_input.source = u238;
  RelActCalcAuto::RelEffCurveInput curve;
  curve.nuclides = { u235_input,u238_input };
  RelActCalcAuto::RoiRange roi;
  roi.lower_energy = 100.0;
  roi.upper_energy = 210.0;
  Options solve_options;
  solve_options.rel_eff_curves = { curve };
  solve_options.rois = { roi };
  solve_options.profile_targets = { frac };
  BOOST_CHECK_EQUAL( solve_options.why_not_usable(),string() );
  Target bad = frac;
  bad.kind = Kind::Age;                       //no age is being fitted
  solve_options.profile_targets = { bad };
  BOOST_CHECK( !solve_options.why_not_usable().empty() );
}


/** The shape that used to be an uncatchable crash: an activity ratio whose DENOMINATOR is a
 mass-fraction-constrained isotope of the numerator's own element.  With the link installed the
 numerator evaluates as `r * activity(denominator)`, and the denominator's sigma-block decode sums
 the element's unconstrained members - the numerator among them - so the evaluation re-enters the
 link and recurses until the stack dies.  It must be refused up front, before any solve. */
BOOST_AUTO_TEST_CASE( circular_ratio_denominator_is_refused_before_solving )
{
  using Target = RelActCalcAuto::Options::ProfileTarget;
  using Kind = Target::Kind;

  initialize_decay_database();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const u234 = db->nuclide("U234");
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  const SandiaDecay::Nuclide * const cs137 = db->nuclide("Cs137");
  BOOST_REQUIRE( u234 && u235 && u238 && cs137 );

  const auto input_for = []( const SandiaDecay::Nuclide * const nuc ) {
    RelActCalcAuto::NucInputInfo input;
    input.source = nuc;
    input.age = 0.0;
    input.fit_age = false;
    return input;
  };

  RelActCalcAuto::RelEffCurveInput curve;
  curve.nuclides = { input_for(u234), input_for(u235), input_for(u238), input_for(cs137) };
  RelActCalcAuto::Options options;
  options.rel_eff_curves = { curve };

  Target ratio;
  ratio.kind = Kind::ActivityRatio;
  ratio.source = u238;          //numerator: unconstrained, so its slot can carry the ratio
  ratio.denominator = u234;
  BOOST_CHECK_EQUAL( ratio.why_not_usable(options), string() );

  // Constraining the DENOMINATOR, an isotope of the numerator's element, closes the cycle.
  RelActCalcAuto::Options circular = options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint window;
  window.nuclide = u234;
  window.lower_mass_fraction = 0.0001;
  window.upper_mass_fraction = 0.01;
  circular.rel_eff_curves[0].mass_fraction_constraints.push_back( window );
  BOOST_CHECK_MESSAGE( !ratio.why_not_usable(circular).empty(),
      "a mass-fraction-constrained same-element denominator must be refused" );

  // The same constraint on a DIFFERENT element's source cannot recurse, so it stays usable.
  RelActCalcAuto::Options other_element = options;
  Target cs_ratio;
  cs_ratio.kind = Kind::ActivityRatio;
  cs_ratio.source = cs137;
  cs_ratio.denominator = u235;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint u235_window;
  u235_window.nuclide = u235;
  u235_window.lower_mass_fraction = 0.10;
  u235_window.upper_mass_fraction = 0.50;
  other_element.rel_eff_curves[0].mass_fraction_constraints.push_back( u235_window );
  BOOST_CHECK_EQUAL( cs_ratio.why_not_usable(other_element), string() );

  // An explicit request for the circular shape fails the whole solve's options check rather than
  // reaching the engine.
  circular.profile_targets = { ratio };
  RelActCalcAuto::RoiRange roi;
  roi.lower_energy = 100.0;
  roi.upper_energy = 210.0;
  circular.rois = { roi };
  BOOST_CHECK( !circular.why_not_usable().empty() );
}


/** Every profile enum spells itself through one shared `to_str`, and a target names itself one
 way, so the report, the LLM tool, warnings, and the text summary cannot drift apart. */
BOOST_AUTO_TEST_CASE( profile_enum_and_target_naming_is_single_sourced )
{
  using Solution = RelActCalcAuto::RelActAutoSolution;
  using Target = RelActCalcAuto::Options::ProfileTarget;

  initialize_decay_database();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE( u235 && u238 );

  const Solution::MassFractionProfileStatus statuses[] = {
    Solution::MassFractionProfileStatus::NotRequested,
    Solution::MassFractionProfileStatus::Complete,
    Solution::MassFractionProfileStatus::BoundaryLimited,
    Solution::MassFractionProfileStatus::NonIdentifiable,
    Solution::MassFractionProfileStatus::Failed };
  set<string> status_names;
  for( const Solution::MassFractionProfileStatus status : statuses )
  {
    const string name = Solution::to_str(status);
    BOOST_CHECK( !name.empty() );
    status_names.insert( name );
  }
  BOOST_CHECK_EQUAL( status_names.size(), 5U );

  // A synthetic scan cap is NOT a feasibility bound, so it must not spell itself as one.
  const Solution::MassFractionProfileEndpointKind kinds[] = {
    Solution::MassFractionProfileEndpointKind::LikelihoodCrossing,
    Solution::MassFractionProfileEndpointKind::PhysicalLimit,
    Solution::MassFractionProfileEndpointKind::InputConstraintLimit,
    Solution::MassFractionProfileEndpointKind::ScanRangeLimit };
  set<string> kind_names;
  for( const Solution::MassFractionProfileEndpointKind kind : kinds )
    kind_names.insert( Solution::to_str(kind) );
  BOOST_CHECK_EQUAL( kind_names.size(), 4U );
  BOOST_CHECK( string(Solution::to_str(Solution::MassFractionProfileEndpointKind::ScanRangeLimit))
               != string(Solution::to_str(Solution::MassFractionProfileEndpointKind::PhysicalLimit)) );

  Target frac;
  frac.kind = Target::Kind::MassFraction;
  frac.source = u235;
  BOOST_CHECK_EQUAL( frac.display_name(), "U235 MassFraction" );

  Target ratio;
  ratio.kind = Target::Kind::ActivityRatio;
  ratio.source = u235;
  ratio.denominator = u238;
  BOOST_CHECK_EQUAL( ratio.display_name(), "U235/U238 ActivityRatio" );

  Target activity;
  activity.kind = Target::Kind::RelativeActivity;
  activity.source = u238;
  BOOST_CHECK_EQUAL( activity.display_name(), "U238 RelativeActivity" );
}
