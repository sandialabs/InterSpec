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

#include <string>
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
