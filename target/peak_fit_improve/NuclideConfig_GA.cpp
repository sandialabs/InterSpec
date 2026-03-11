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

#include <deque>
#include <mutex>
#include <atomic>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <functional>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PeakFit.h"

#include "NuclideConfig_GA.h"
#include "InitialFit_GA.h"
#include "FinalFit_GA.h"
#include "CandidatePeak_GA.h"
#include "ClassifyDetType_GA.h"
#include "FitPeaksForNuclideDev.h"

using namespace std;


namespace NuclideConfig_GA
{

// Module-level state for the GA (following InitialFit_GA pattern)
static std::function<double( const PeakFitForNuclideConfig & )> ns_ga_eval_fcn;
static std::atomic<size_t> ns_num_evals_this_generation{0};
static bool sm_has_been_called = false;
static bool sm_set_best_genes = false;
static NuclideConfigSolution sm_best_genes;
static double sm_best_total_cost = std::numeric_limits<double>::max();
static std::ofstream sm_output_file;

// Store precomputed data pointer for use by SO_report_generation
static const std::vector<PrecomputedNuclideData> *sm_precomputed_ptr = nullptr;


std::vector<RelActCalcAuto::SrcVariant> resolve_sources( const std::string &src_name )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    throw runtime_error( "Failed to open SandiaDecayDataBase" );

  // Strip shielding/config suffix from name
  string name = src_name;
  const auto underscore_pos = name.find( '_' );
  if( underscore_pos != string::npos )
    name = name.substr( 0, underscore_pos );

  vector<RelActCalcAuto::SrcVariant> sources;

  if( name == "Tl201woTl202" || name == "Tl201wTl202" || name == "Tl201" )
  {
    sources.push_back( db->nuclide( "Tl201" ) );
    sources.push_back( db->nuclide( "Tl202" ) );
  }
  else if( name == "I125" )
  {
    sources.push_back( db->nuclide( "I125" ) );
    sources.push_back( db->nuclide( "I126" ) );
  }
  else if( name == "U233" )
  {
    sources.push_back( db->element( "U" ) );
    sources.push_back( db->nuclide( "U232" ) );
    sources.push_back( db->nuclide( "U233" ) );
  }
  else if( name == "Pu238" )
  {
    sources.push_back( db->element( "Pu" ) );
    sources.push_back( db->nuclide( "Pu238" ) );
    sources.push_back( db->nuclide( "Pu239" ) );
    sources.push_back( db->nuclide( "Pu241" ) );
  }
  else if( name == "Pu239" )
  {
    sources.push_back( db->element( "Pu" ) );
    sources.push_back( db->nuclide( "Pu239" ) );
    sources.push_back( db->nuclide( "Pu241" ) );
  }
  else if( name == "Uore" )
  {
    sources.push_back( db->element( "U" ) );
    sources.push_back( db->nuclide( "U235" ) );
    sources.push_back( db->nuclide( "U238" ) );
    sources.push_back( db->nuclide( "Ra226" ) );
  }
  else if( name == "U235" )
  {
    sources.push_back( db->element( "U" ) );
    sources.push_back( db->nuclide( "U235" ) );
    sources.push_back( db->nuclide( "U238" ) );
  }
  else if( name == "Np237" )
  {
    sources.push_back( db->element( "Np" ) );
    sources.push_back( db->nuclide( "Np237" ) );
  }
  else if( name == "Xe133" )
  {
    sources.push_back( db->nuclide( "Xe133" ) );
    sources.push_back( db->nuclide( "Xe133m" ) );
  }
  else if( name == "Cf252" || name == "Am241Li" )
  {
    // Skip these sources
    return sources;
  }
  else
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( name );
    if( !nuc )
    {
      cout << "Unable to get nuclide for '" << src_name << "' - so will skip" << endl;
      return sources;
    }
    sources.push_back( nuc );
  }

  return sources;
}//resolve_sources


std::vector<PrecomputedNuclideData> precompute_nuclide_data(
  const std::vector<DataSrcInfo> &srcs_info,
  const BackgroundMode bg_mode )
{
  std::vector<PrecomputedNuclideData> result;
  result.reserve( srcs_info.size() );

  cout << "Precomputing auto-search peaks for " << srcs_info.size() << " spectra..." << endl;

  for( size_t i = 0; i < srcs_info.size(); ++i )
  {
    const DataSrcInfo &info = srcs_info[i];
    const InjectSourceInfo &src = info.src_info;

    if( src.src_spectra.empty() )
      continue;

    // Resolve sources
    std::vector<RelActCalcAuto::SrcVariant> sources = resolve_sources( src.src_name );
    if( sources.empty() || RelActCalcAuto::is_null( sources.front() ) )
      continue;

    PrecomputedNuclideData pd;
    pd.src_info = &info;
    pd.sources = std::move( sources );
    pd.foreground = src.src_spectra.front();
    pd.drf = nullptr;

    // Set background based on mode
    switch( bg_mode )
    {
      case BackgroundMode::BackgroundSubtracted:
        pd.background = src.long_background;
        break;
      case BackgroundMode::NoBackground:
      case BackgroundMode::NoBackgroundFitNorm:
        pd.background = nullptr;
        break;
    }//switch( bg_mode )

    // Use det_type and skew_type from the DataSrcInfo (set during data loading)
    pd.det_type = info.det_type;

    // Create PeakFitDetPrefs with correct type
    pd.peak_fit_prefs = std::make_shared<PeakFitDetPrefs>();
    pd.peak_fit_prefs->m_det_type = pd.det_type;

    const bool singleThreaded = false;
    std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> dummy_origpeaks;

    // The expensive call - done once
    pd.auto_search_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks(
      pd.foreground, pd.drf, dummy_origpeaks, singleThreaded, pd.peak_fit_prefs );

    result.push_back( std::move( pd ) );

    if( ((i + 1) % 10) == 0 || (i + 1) == srcs_info.size() )
      cout << "  Precomputed " << (i + 1) << " of " << srcs_info.size() << " spectra" << endl;
  }//for( size_t i = 0; i < srcs_info.size(); ++i )

  cout << "Precomputation complete: " << result.size() << " spectra ready for GA evaluation." << endl;

  return result;
}//precompute_nuclide_data


std::string NuclideConfigSolution::to_string( const std::string &separator ) const
{
  std::ostringstream oss;
  oss << std::setprecision(6);
  oss << "fwhm_functional_form=" << fwhm_functional_form << separator
      << "rel_eff_manual_base_rel_eff_uncert=" << rel_eff_manual_base_rel_eff_uncert << separator
      << "initial_nuc_match_cluster_num_sigma=" << initial_nuc_match_cluster_num_sigma << separator
      << "manual_eff_cluster_num_sigma=" << manual_eff_cluster_num_sigma << separator
      << "initial_manual_relEff_1peak_eqn_order=" << initial_manual_relEff_1peak_eqn_order << separator
      << "initial_manual_relEff_1peak_form=" << initial_manual_relEff_1peak_form << separator
      << "initial_manual_relEff_2peak_eqn_order=" << initial_manual_relEff_2peak_eqn_order << separator
      << "initial_manual_relEff_2peak_form=" << initial_manual_relEff_2peak_form << separator
      << "initial_manual_relEff_3peak_eqn_order=" << initial_manual_relEff_3peak_eqn_order << separator
      << "initial_manual_relEff_3peak_form=" << initial_manual_relEff_3peak_form << separator
      << "initial_manual_relEff_4peak_physical_use_hoerl=" << initial_manual_relEff_4peak_physical_use_hoerl << separator
      << "initial_manual_relEff_4peak_eqn_order=" << initial_manual_relEff_4peak_eqn_order << separator
      << "initial_manual_relEff_4peak_form=" << initial_manual_relEff_4peak_form << separator
      << "initial_manual_relEff_many_peak_physical_use_hoerl=" << initial_manual_relEff_many_peak_physical_use_hoerl << separator
      << "initial_manual_relEff_many_peak_eqn_order=" << initial_manual_relEff_many_peak_eqn_order << separator
      << "initial_manual_relEff_manypeak_form=" << initial_manual_relEff_manypeak_form << separator
      << "manual_rel_eff_sol_min_data_area_keep=" << manual_rel_eff_sol_min_data_area_keep << separator
      << "manual_rel_eff_sol_min_est_peak_area_keep=" << manual_rel_eff_sol_min_est_peak_area_keep << separator
      << "manual_rel_eff_sol_min_est_significance_keep=" << manual_rel_eff_sol_min_est_significance_keep << separator
      << "manual_rel_eff_sol_min_fwhm_roi=" << manual_rel_eff_sol_min_fwhm_roi << separator
      << "manual_rel_eff_sol_min_fwhm_quad_cont=" << manual_rel_eff_sol_min_fwhm_quad_cont << separator
      << "manual_rel_eff_sol_max_fwhm=" << manual_rel_eff_sol_max_fwhm << separator
      << "manual_rel_eff_roi_width_num_fwhm_lower=" << manual_rel_eff_roi_width_num_fwhm_lower << separator
      << "manual_rel_eff_roi_width_num_fwhm_upper=" << manual_rel_eff_roi_width_num_fwhm_upper << separator
      << "fwhm_form=" << fwhm_form << separator
      << "rel_eff_auto_base_rel_eff_uncert=" << rel_eff_auto_base_rel_eff_uncert << separator
      << "auto_rel_eff_cluster_num_sigma=" << auto_rel_eff_cluster_num_sigma << separator
      << "auto_rel_eff_sol_min_data_area_keep=" << auto_rel_eff_sol_min_data_area_keep << separator
      << "auto_rel_eff_sol_min_est_peak_area_keep=" << auto_rel_eff_sol_min_est_peak_area_keep << separator
      << "auto_rel_eff_sol_min_est_significance_keep=" << auto_rel_eff_sol_min_est_significance_keep << separator
      << "auto_rel_eff_roi_width_num_fwhm_lower=" << auto_rel_eff_roi_width_num_fwhm_lower << separator
      << "auto_rel_eff_roi_width_num_fwhm_upper=" << auto_rel_eff_roi_width_num_fwhm_upper << separator
      << "auto_rel_eff_sol_max_fwhm=" << auto_rel_eff_sol_max_fwhm << separator
      << "auto_rel_eff_sol_min_fwhm_roi=" << auto_rel_eff_sol_min_fwhm_roi << separator
      << "auto_rel_eff_sol_min_fwhm_quad_cont=" << auto_rel_eff_sol_min_fwhm_quad_cont << separator
      << "rel_eff_eqn_type=" << rel_eff_eqn_type << separator
      << "rel_eff_eqn_order=" << rel_eff_eqn_order << separator
      << "desperation_phys_model_atomic_number=" << desperation_phys_model_atomic_number << separator
      << "desperation_phys_model_areal_density_g_per_cm2=" << desperation_phys_model_areal_density_g_per_cm2 << separator
      << "nucs_of_el_same_age=" << nucs_of_el_same_age << separator
      << "phys_model_use_hoerl=" << phys_model_use_hoerl << separator
      << "fit_energy_cal=" << fit_energy_cal << separator
      << "roi_significance_min_chi2_reduction=" << roi_significance_min_chi2_reduction << separator
      << "roi_significance_min_peak_sig=" << roi_significance_min_peak_sig << separator
      << "observable_peak_initial_significance_threshold=" << observable_peak_initial_significance_threshold << separator
      << "observable_peak_final_significance_threshold=" << observable_peak_final_significance_threshold << separator
      << "step_cont_min_peak_area=" << step_cont_min_peak_area << separator
      << "step_cont_min_peak_significance=" << step_cont_min_peak_significance << separator
      << "step_cont_left_right_nsigma=" << step_cont_left_right_nsigma;
  return oss.str();
}//NuclideConfigSolution::to_string


PeakFitForNuclideConfig genes_to_settings( const NuclideConfigSolution &p )
{
  PeakFitForNuclideConfig config;

  config.fwhm_functional_form = static_cast<DetectorPeakResponse::ResolutionFnctForm>(
    std::clamp( p.fwhm_functional_form, 0, static_cast<int>(DetectorPeakResponse::kNumResolutionFnctForm) - 1 ) );

  config.rel_eff_manual_base_rel_eff_uncert = p.rel_eff_manual_base_rel_eff_uncert;
  config.initial_nuc_match_cluster_num_sigma = p.initial_nuc_match_cluster_num_sigma;
  config.manual_eff_cluster_num_sigma = p.manual_eff_cluster_num_sigma;

  config.initial_manual_relEff_1peak_eqn_order = static_cast<size_t>( std::clamp( p.initial_manual_relEff_1peak_eqn_order, 0, 2 ) );
  config.initial_manual_relEff_1peak_form = static_cast<RelActCalc::RelEffEqnForm>( std::clamp( p.initial_manual_relEff_1peak_form, 0, 3 ) );
  config.initial_manual_relEff_2peak_eqn_order = static_cast<size_t>( std::clamp( p.initial_manual_relEff_2peak_eqn_order, 0, 3 ) );
  config.initial_manual_relEff_2peak_form = static_cast<RelActCalc::RelEffEqnForm>( std::clamp( p.initial_manual_relEff_2peak_form, 0, 3 ) );
  config.initial_manual_relEff_3peak_eqn_order = static_cast<size_t>( std::clamp( p.initial_manual_relEff_3peak_eqn_order, 0, 4 ) );
  config.initial_manual_relEff_3peak_form = static_cast<RelActCalc::RelEffEqnForm>( std::clamp( p.initial_manual_relEff_3peak_form, 0, 3 ) );
  config.initial_manual_relEff_4peak_physical_use_hoerl = (p.initial_manual_relEff_4peak_physical_use_hoerl != 0);
  config.initial_manual_relEff_4peak_eqn_order = static_cast<size_t>( std::clamp( p.initial_manual_relEff_4peak_eqn_order, 0, 5 ) );
  config.initial_manual_relEff_4peak_form = static_cast<RelActCalc::RelEffEqnForm>( std::clamp( p.initial_manual_relEff_4peak_form, 0, 3 ) );
  config.initial_manual_relEff_many_peak_physical_use_hoerl = (p.initial_manual_relEff_many_peak_physical_use_hoerl != 0);
  config.initial_manual_relEff_many_peak_eqn_order = static_cast<size_t>( std::clamp( p.initial_manual_relEff_many_peak_eqn_order, 0, 6 ) );
  config.initial_manual_relEff_manypeak_form = static_cast<RelActCalc::RelEffEqnForm>( std::clamp( p.initial_manual_relEff_manypeak_form, 0, 3 ) );

  config.manual_rel_eff_sol_min_data_area_keep = p.manual_rel_eff_sol_min_data_area_keep;
  config.manual_rel_eff_sol_min_est_peak_area_keep = p.manual_rel_eff_sol_min_est_peak_area_keep;
  config.manual_rel_eff_sol_min_est_significance_keep = p.manual_rel_eff_sol_min_est_significance_keep;
  config.manual_rel_eff_sol_min_fwhm_roi = p.manual_rel_eff_sol_min_fwhm_roi;
  config.manual_rel_eff_sol_min_fwhm_quad_cont = p.manual_rel_eff_sol_min_fwhm_quad_cont;
  config.manual_rel_eff_sol_max_fwhm = p.manual_rel_eff_sol_max_fwhm;
  config.manual_rel_eff_roi_width_num_fwhm_lower = p.manual_rel_eff_roi_width_num_fwhm_lower;
  config.manual_rel_eff_roi_width_num_fwhm_upper = p.manual_rel_eff_roi_width_num_fwhm_upper;

  // FwhmForm: limited to Berstein_2 through Berstein_5
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2) == 8, "Berstein_2 enum value changed" );
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) == 11, "Berstein_5 enum value changed" );
  config.fwhm_form = static_cast<RelActCalcAuto::FwhmForm>(
    std::clamp( p.fwhm_form,
                static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2),
                static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) ) );

  config.rel_eff_auto_base_rel_eff_uncert = p.rel_eff_auto_base_rel_eff_uncert;
  config.auto_rel_eff_cluster_num_sigma = p.auto_rel_eff_cluster_num_sigma;
  config.auto_rel_eff_sol_min_data_area_keep = p.auto_rel_eff_sol_min_data_area_keep;
  config.auto_rel_eff_sol_min_est_peak_area_keep = p.auto_rel_eff_sol_min_est_peak_area_keep;
  config.auto_rel_eff_sol_min_est_significance_keep = p.auto_rel_eff_sol_min_est_significance_keep;
  config.auto_rel_eff_roi_width_num_fwhm_lower = p.auto_rel_eff_roi_width_num_fwhm_lower;
  config.auto_rel_eff_roi_width_num_fwhm_upper = p.auto_rel_eff_roi_width_num_fwhm_upper;
  config.auto_rel_eff_sol_max_fwhm = p.auto_rel_eff_sol_max_fwhm;
  config.auto_rel_eff_sol_min_fwhm_roi = p.auto_rel_eff_sol_min_fwhm_roi;
  config.auto_rel_eff_sol_min_fwhm_quad_cont = p.auto_rel_eff_sol_min_fwhm_quad_cont;

  // RelEffEqnForm: 0=LnX, 1=LnY, 2=LnXLnY, 3=FramEmpirical, 4=FramPhysicalModel
  config.rel_eff_eqn_type = static_cast<RelActCalc::RelEffEqnForm>(
    std::clamp( p.rel_eff_eqn_type, 0, static_cast<int>(RelActCalc::RelEffEqnForm::FramPhysicalModel) ) );
  config.rel_eff_eqn_order = static_cast<size_t>( std::clamp( p.rel_eff_eqn_order, 0, 6 ) );

  config.desperation_phys_model_atomic_number = p.desperation_phys_model_atomic_number;
  config.desperation_phys_model_areal_density_g_per_cm2 = p.desperation_phys_model_areal_density_g_per_cm2;

  config.nucs_of_el_same_age = (p.nucs_of_el_same_age != 0);
  config.phys_model_use_hoerl = (p.phys_model_use_hoerl != 0);
  config.fit_energy_cal = (p.fit_energy_cal != 0);

  config.roi_significance_min_chi2_reduction = p.roi_significance_min_chi2_reduction;
  config.roi_significance_min_peak_sig = p.roi_significance_min_peak_sig;
  config.observable_peak_initial_significance_threshold = p.observable_peak_initial_significance_threshold;
  config.observable_peak_final_significance_threshold = p.observable_peak_final_significance_threshold;

  config.step_cont_min_peak_area = p.step_cont_min_peak_area;
  config.step_cont_min_peak_significance = p.step_cont_min_peak_significance;
  config.step_cont_left_right_nsigma = p.step_cont_left_right_nsigma;

  // skew_type is not GA-optimized - leave at default (NoSkew)
  // Caller should override per detector type (e.g., DoubleSidedCrystalBall for CZT)

  return config;
}//genes_to_settings


void init_genes( NuclideConfigSolution &p, const std::function<double(void)> &rnd01 )
{
  // FWHM functional form [0, 3]
  p.fwhm_functional_form = static_cast<int>( 4 * rnd01() );
  if( p.fwhm_functional_form > 3 ) p.fwhm_functional_form = 3;

  // Manual RelEff parameters
  p.rel_eff_manual_base_rel_eff_uncert = 0.0 + 0.5 * rnd01();
  p.initial_nuc_match_cluster_num_sigma = 0.5 + 3.5 * rnd01();
  p.manual_eff_cluster_num_sigma = 0.5 + 3.5 * rnd01();

  // Manual RelEff equation forms/orders
  p.initial_manual_relEff_1peak_eqn_order = static_cast<int>( 3 * rnd01() );
  if( p.initial_manual_relEff_1peak_eqn_order > 2 ) p.initial_manual_relEff_1peak_eqn_order = 2;
  p.initial_manual_relEff_1peak_form = static_cast<int>( 4 * rnd01() );
  if( p.initial_manual_relEff_1peak_form > 3 ) p.initial_manual_relEff_1peak_form = 3;

  p.initial_manual_relEff_2peak_eqn_order = static_cast<int>( 4 * rnd01() );
  if( p.initial_manual_relEff_2peak_eqn_order > 3 ) p.initial_manual_relEff_2peak_eqn_order = 3;
  p.initial_manual_relEff_2peak_form = static_cast<int>( 4 * rnd01() );
  if( p.initial_manual_relEff_2peak_form > 3 ) p.initial_manual_relEff_2peak_form = 3;

  p.initial_manual_relEff_3peak_eqn_order = static_cast<int>( 5 * rnd01() );
  if( p.initial_manual_relEff_3peak_eqn_order > 4 ) p.initial_manual_relEff_3peak_eqn_order = 4;
  p.initial_manual_relEff_3peak_form = static_cast<int>( 4 * rnd01() );
  if( p.initial_manual_relEff_3peak_form > 3 ) p.initial_manual_relEff_3peak_form = 3;

  p.initial_manual_relEff_4peak_physical_use_hoerl = rnd01() < 0.5 ? 0 : 1;
  p.initial_manual_relEff_4peak_eqn_order = static_cast<int>( 6 * rnd01() );
  if( p.initial_manual_relEff_4peak_eqn_order > 5 ) p.initial_manual_relEff_4peak_eqn_order = 5;
  p.initial_manual_relEff_4peak_form = static_cast<int>( 4 * rnd01() );
  if( p.initial_manual_relEff_4peak_form > 3 ) p.initial_manual_relEff_4peak_form = 3;

  p.initial_manual_relEff_many_peak_physical_use_hoerl = rnd01() < 0.5 ? 0 : 1;
  p.initial_manual_relEff_many_peak_eqn_order = static_cast<int>( 7 * rnd01() );
  if( p.initial_manual_relEff_many_peak_eqn_order > 6 ) p.initial_manual_relEff_many_peak_eqn_order = 6;
  p.initial_manual_relEff_manypeak_form = static_cast<int>( 4 * rnd01() );
  if( p.initial_manual_relEff_manypeak_form > 3 ) p.initial_manual_relEff_manypeak_form = 3;

  // Manual clustering thresholds
  p.manual_rel_eff_sol_min_data_area_keep = 1.0 + 99.0 * rnd01();
  p.manual_rel_eff_sol_min_est_peak_area_keep = 1.0 + 49.0 * rnd01();
  p.manual_rel_eff_sol_min_est_significance_keep = 0.5 + 5.5 * rnd01();
  p.manual_rel_eff_sol_min_fwhm_roi = 0.5 + 2.5 * rnd01();
  p.manual_rel_eff_sol_min_fwhm_quad_cont = 3.0 + 12.0 * rnd01();
  p.manual_rel_eff_sol_max_fwhm = 5.0 + 25.0 * rnd01();
  p.manual_rel_eff_roi_width_num_fwhm_lower = 1.0 + 5.0 * rnd01();
  p.manual_rel_eff_roi_width_num_fwhm_upper = 1.0 + 5.0 * rnd01();

  // Auto RelEff parameters
  // FwhmForm: limited to Berstein_2(8) through Berstein_5(11)
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2) == 8, "Berstein_2 enum value changed" );
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) == 11, "Berstein_5 enum value changed" );
  p.fwhm_form = 8 + static_cast<int>( 4 * rnd01() );
  if( p.fwhm_form > 11 ) p.fwhm_form = 11;

  p.rel_eff_auto_base_rel_eff_uncert = 0.0 + 0.5 * rnd01();
  p.auto_rel_eff_cluster_num_sigma = 0.5 + 4.5 * rnd01();
  p.auto_rel_eff_sol_min_data_area_keep = 1.0 + 99.0 * rnd01();
  p.auto_rel_eff_sol_min_est_peak_area_keep = 1.0 + 49.0 * rnd01();
  p.auto_rel_eff_sol_min_est_significance_keep = 0.5 + 7.5 * rnd01();
  p.auto_rel_eff_roi_width_num_fwhm_lower = 1.0 + 6.0 * rnd01();
  p.auto_rel_eff_roi_width_num_fwhm_upper = 1.0 + 6.0 * rnd01();
  p.auto_rel_eff_sol_max_fwhm = 4.0 + 21.0 * rnd01();
  p.auto_rel_eff_sol_min_fwhm_roi = 0.5 + 2.5 * rnd01();
  p.auto_rel_eff_sol_min_fwhm_quad_cont = 3.0 + 12.0 * rnd01();

  // RelActAuto model
  // RelEffEqnForm: 0..4 (LnX, LnY, LnXLnY, FramEmpirical, FramPhysicalModel)
  p.rel_eff_eqn_type = static_cast<int>( 5 * rnd01() );
  if( p.rel_eff_eqn_type > 4 ) p.rel_eff_eqn_type = 4;
  p.rel_eff_eqn_order = static_cast<int>( 7 * rnd01() );
  if( p.rel_eff_eqn_order > 6 ) p.rel_eff_eqn_order = 6;

  // Desperation physical model
  p.desperation_phys_model_atomic_number = 6.0 + 76.0 * rnd01();
  p.desperation_phys_model_areal_density_g_per_cm2 = 0.1 + 19.9 * rnd01();

  // Booleans
  p.nucs_of_el_same_age = rnd01() < 0.5 ? 0 : 1;
  p.phys_model_use_hoerl = rnd01() < 0.5 ? 0 : 1;
  p.fit_energy_cal = rnd01() < 0.5 ? 0 : 1;

  // ROI significance and observable peak thresholds
  p.roi_significance_min_chi2_reduction = 2.0 + 28.0 * rnd01();
  p.roi_significance_min_peak_sig = 1.0 + 7.0 * rnd01();
  p.observable_peak_initial_significance_threshold = 1.0 + 4.0 * rnd01();
  p.observable_peak_final_significance_threshold = 0.5 + 4.5 * rnd01();

  // Step continuum
  p.step_cont_min_peak_area = 100.0 + 4900.0 * rnd01();
  p.step_cont_min_peak_significance = 5.0 + 95.0 * rnd01();
  p.step_cont_left_right_nsigma = 1.0 + 7.0 * rnd01();

  // skew_type is not GA-optimized
}//init_genes


bool eval_solution( const NuclideConfigSolution &p, NuclideConfigCost &c )
{
  const PeakFitForNuclideConfig config = genes_to_settings( p );

  c.objective1 = ns_ga_eval_fcn( config );

  ns_num_evals_this_generation += 1;
  if( (ns_num_evals_this_generation % 10) == 0 )
    cout << "Have evaluated " << ns_num_evals_this_generation.load() << " individuals this generation." << endl;

  if( std::isnan( c.objective1 ) || std::isinf( c.objective1 ) )
  {
    cerr << "Got an objective of " << c.objective1 << " for " << p.to_string( ", " ) << endl;
    return false;
  }

  return true;
}//eval_solution


// Helper macro to reduce repetition in mutate for double genes
// Each gene: if random < threshold, apply mutation, check range
#define MUTATE_DOUBLE( field, lo, hi ) \
  if( rnd01() < mutate_threshold ) { \
    X_new.field += mu * (rnd01() - rnd01()); \
    in_range = in_range && (X_new.field >= (lo) && X_new.field < (hi)); \
  }

// For integer genes: randomly reassign with small probability
#define MUTATE_INT( field, lo, hi ) \
  if( rnd01() < mutate_threshold ) { \
    X_new.field = static_cast<int>( (lo) + ((hi) - (lo) + 1) * rnd01() ); \
    if( X_new.field > (hi) ) X_new.field = (hi); \
    in_range = in_range && (X_new.field >= (lo) && X_new.field <= (hi)); \
  }

// For boolean genes: toggle with small probability
#define MUTATE_BOOL( field ) \
  if( rnd01() < mutate_threshold ) { \
    X_new.field = (X_new.field == 0) ? 1 : 0; \
  }


NuclideConfigSolution mutate( const NuclideConfigSolution &X_base,
                              const std::function<double(void)> &rnd01,
                              double shrink_scale )
{
  NuclideConfigSolution X_new;
  const double mu = 0.2 * shrink_scale;
  const double mutate_threshold = PeakFitImprove::sm_ga_mutate_threshold;
  bool in_range;

  size_t num_tries = 0;

  do
  {
    num_tries += 1;
    if( num_tries > 1000 )
    {
      cerr << "NuclideConfig_GA::mutate: took over " << num_tries << " tries - reinitializing." << endl;
      std::random_device rng;
      std::uniform_real_distribution<double> unif_dist( 0.0, 1.0 );
      init_genes( X_new, [&](){ return unif_dist( rng ); } );
      return X_new;
    }

    in_range = true;
    X_new = X_base;

    // FWHM functional form
    MUTATE_INT( fwhm_functional_form, 0, 3 )

    // Manual RelEff doubles
    MUTATE_DOUBLE( rel_eff_manual_base_rel_eff_uncert, 0.0, 0.5 )
    MUTATE_DOUBLE( initial_nuc_match_cluster_num_sigma, 0.5, 4.0 )
    MUTATE_DOUBLE( manual_eff_cluster_num_sigma, 0.5, 4.0 )

    // Manual equation forms/orders (ints)
    MUTATE_INT( initial_manual_relEff_1peak_eqn_order, 0, 2 )
    MUTATE_INT( initial_manual_relEff_1peak_form, 0, 3 )
    MUTATE_INT( initial_manual_relEff_2peak_eqn_order, 0, 3 )
    MUTATE_INT( initial_manual_relEff_2peak_form, 0, 3 )
    MUTATE_INT( initial_manual_relEff_3peak_eqn_order, 0, 4 )
    MUTATE_INT( initial_manual_relEff_3peak_form, 0, 3 )
    MUTATE_BOOL( initial_manual_relEff_4peak_physical_use_hoerl )
    MUTATE_INT( initial_manual_relEff_4peak_eqn_order, 0, 5 )
    MUTATE_INT( initial_manual_relEff_4peak_form, 0, 3 )
    MUTATE_BOOL( initial_manual_relEff_many_peak_physical_use_hoerl )
    MUTATE_INT( initial_manual_relEff_many_peak_eqn_order, 0, 6 )
    MUTATE_INT( initial_manual_relEff_manypeak_form, 0, 3 )

    // Manual clustering doubles
    MUTATE_DOUBLE( manual_rel_eff_sol_min_data_area_keep, 1.0, 100.0 )
    MUTATE_DOUBLE( manual_rel_eff_sol_min_est_peak_area_keep, 1.0, 50.0 )
    MUTATE_DOUBLE( manual_rel_eff_sol_min_est_significance_keep, 0.5, 6.0 )
    MUTATE_DOUBLE( manual_rel_eff_sol_min_fwhm_roi, 0.5, 3.0 )
    MUTATE_DOUBLE( manual_rel_eff_sol_min_fwhm_quad_cont, 3.0, 15.0 )
    MUTATE_DOUBLE( manual_rel_eff_sol_max_fwhm, 5.0, 30.0 )
    MUTATE_DOUBLE( manual_rel_eff_roi_width_num_fwhm_lower, 1.0, 6.0 )
    MUTATE_DOUBLE( manual_rel_eff_roi_width_num_fwhm_upper, 1.0, 6.0 )

    // Auto RelEff - fwhm_form limited to Berstein_2 through Berstein_5
    static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2) == 8, "Berstein_2 enum value changed" );
    static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) == 11, "Berstein_5 enum value changed" );
    MUTATE_INT( fwhm_form, 8, 11 )
    MUTATE_DOUBLE( rel_eff_auto_base_rel_eff_uncert, 0.0, 0.5 )
    MUTATE_DOUBLE( auto_rel_eff_cluster_num_sigma, 0.5, 5.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_min_data_area_keep, 1.0, 100.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_min_est_peak_area_keep, 1.0, 50.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_min_est_significance_keep, 0.5, 8.0 )
    MUTATE_DOUBLE( auto_rel_eff_roi_width_num_fwhm_lower, 1.0, 7.0 )
    MUTATE_DOUBLE( auto_rel_eff_roi_width_num_fwhm_upper, 1.0, 7.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_max_fwhm, 4.0, 25.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_min_fwhm_roi, 0.5, 3.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_min_fwhm_quad_cont, 3.0, 15.0 )

    // RelActAuto model
    MUTATE_INT( rel_eff_eqn_type, 0, 4 )
    MUTATE_INT( rel_eff_eqn_order, 0, 6 )

    // Desperation physical model
    MUTATE_DOUBLE( desperation_phys_model_atomic_number, 6.0, 82.0 )
    MUTATE_DOUBLE( desperation_phys_model_areal_density_g_per_cm2, 0.1, 20.0 )

    // Booleans
    MUTATE_BOOL( nucs_of_el_same_age )
    MUTATE_BOOL( phys_model_use_hoerl )
    MUTATE_BOOL( fit_energy_cal )

    // ROI significance
    MUTATE_DOUBLE( roi_significance_min_chi2_reduction, 2.0, 30.0 )
    MUTATE_DOUBLE( roi_significance_min_peak_sig, 1.0, 8.0 )
    MUTATE_DOUBLE( observable_peak_initial_significance_threshold, 1.0, 5.0 )
    MUTATE_DOUBLE( observable_peak_final_significance_threshold, 0.5, 5.0 )

    // Step continuum
    MUTATE_DOUBLE( step_cont_min_peak_area, 100.0, 5000.0 )
    MUTATE_DOUBLE( step_cont_min_peak_significance, 5.0, 100.0 )
    MUTATE_DOUBLE( step_cont_left_right_nsigma, 1.0, 8.0 )

  } while( !in_range );

  return X_new;
}//mutate

#undef MUTATE_DOUBLE
#undef MUTATE_INT
#undef MUTATE_BOOL


// Helper for crossover of double genes
#define CROSSOVER_DOUBLE( field ) \
  if( rnd01() < crossover_threshold ) { \
    const double r = rnd01(); \
    X_new.field = r * X1.field + (1.0 - r) * X2.field; \
  } else if( rnd01() < 0.5 ) { \
    X_new.field = X2.field; \
  }

// For int genes: pick from one parent
#define CROSSOVER_INT( field ) \
  if( rnd01() < 0.5 ) { \
    X_new.field = X2.field; \
  }


NuclideConfigSolution crossover( const NuclideConfigSolution &X1,
                                 const NuclideConfigSolution &X2,
                                 const std::function<double(void)> &rnd01 )
{
  const double crossover_threshold = PeakFitImprove::sm_ga_crossover_threshold;

  NuclideConfigSolution X_new = X1;

  // FWHM functional form
  CROSSOVER_INT( fwhm_functional_form )

  // Manual RelEff doubles
  CROSSOVER_DOUBLE( rel_eff_manual_base_rel_eff_uncert )
  CROSSOVER_DOUBLE( initial_nuc_match_cluster_num_sigma )
  CROSSOVER_DOUBLE( manual_eff_cluster_num_sigma )

  // Manual equation forms/orders
  CROSSOVER_INT( initial_manual_relEff_1peak_eqn_order )
  CROSSOVER_INT( initial_manual_relEff_1peak_form )
  CROSSOVER_INT( initial_manual_relEff_2peak_eqn_order )
  CROSSOVER_INT( initial_manual_relEff_2peak_form )
  CROSSOVER_INT( initial_manual_relEff_3peak_eqn_order )
  CROSSOVER_INT( initial_manual_relEff_3peak_form )
  CROSSOVER_INT( initial_manual_relEff_4peak_physical_use_hoerl )
  CROSSOVER_INT( initial_manual_relEff_4peak_eqn_order )
  CROSSOVER_INT( initial_manual_relEff_4peak_form )
  CROSSOVER_INT( initial_manual_relEff_many_peak_physical_use_hoerl )
  CROSSOVER_INT( initial_manual_relEff_many_peak_eqn_order )
  CROSSOVER_INT( initial_manual_relEff_manypeak_form )

  // Manual clustering doubles
  CROSSOVER_DOUBLE( manual_rel_eff_sol_min_data_area_keep )
  CROSSOVER_DOUBLE( manual_rel_eff_sol_min_est_peak_area_keep )
  CROSSOVER_DOUBLE( manual_rel_eff_sol_min_est_significance_keep )
  CROSSOVER_DOUBLE( manual_rel_eff_sol_min_fwhm_roi )
  CROSSOVER_DOUBLE( manual_rel_eff_sol_min_fwhm_quad_cont )
  CROSSOVER_DOUBLE( manual_rel_eff_sol_max_fwhm )
  CROSSOVER_DOUBLE( manual_rel_eff_roi_width_num_fwhm_lower )
  CROSSOVER_DOUBLE( manual_rel_eff_roi_width_num_fwhm_upper )

  // Auto RelEff
  CROSSOVER_INT( fwhm_form )
  CROSSOVER_DOUBLE( rel_eff_auto_base_rel_eff_uncert )
  CROSSOVER_DOUBLE( auto_rel_eff_cluster_num_sigma )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_min_data_area_keep )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_min_est_peak_area_keep )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_min_est_significance_keep )
  CROSSOVER_DOUBLE( auto_rel_eff_roi_width_num_fwhm_lower )
  CROSSOVER_DOUBLE( auto_rel_eff_roi_width_num_fwhm_upper )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_max_fwhm )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_min_fwhm_roi )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_min_fwhm_quad_cont )

  // RelActAuto model
  CROSSOVER_INT( rel_eff_eqn_type )
  CROSSOVER_INT( rel_eff_eqn_order )

  // Desperation physical model
  CROSSOVER_DOUBLE( desperation_phys_model_atomic_number )
  CROSSOVER_DOUBLE( desperation_phys_model_areal_density_g_per_cm2 )

  // Booleans
  CROSSOVER_INT( nucs_of_el_same_age )
  CROSSOVER_INT( phys_model_use_hoerl )
  CROSSOVER_INT( fit_energy_cal )

  // ROI significance
  CROSSOVER_DOUBLE( roi_significance_min_chi2_reduction )
  CROSSOVER_DOUBLE( roi_significance_min_peak_sig )
  CROSSOVER_DOUBLE( observable_peak_initial_significance_threshold )
  CROSSOVER_DOUBLE( observable_peak_final_significance_threshold )

  // Step continuum
  CROSSOVER_DOUBLE( step_cont_min_peak_area )
  CROSSOVER_DOUBLE( step_cont_min_peak_significance )
  CROSSOVER_DOUBLE( step_cont_left_right_nsigma )

  return X_new;
}//crossover

#undef CROSSOVER_DOUBLE
#undef CROSSOVER_INT


double calculate_SO_total_fitness( const GA_Type::thisChromosomeType &X )
{
  double final_cost = 0.0;
  final_cost += X.middle_costs.objective1;
  return final_cost;
}//calculate_SO_total_fitness


void SO_report_generation( int generation_number,
                           const EA::GenerationType<NuclideConfigSolution,NuclideConfigCost> &last_generation,
                           const NuclideConfigSolution &best_genes )
{
  bool best_yet = false;
  if( !sm_set_best_genes || (last_generation.best_total_cost < sm_best_total_cost) )
  {
    best_yet = true;
    sm_set_best_genes = true;
    sm_best_genes = best_genes;
    sm_best_total_cost = last_generation.best_total_cost;
  }

  sm_output_file
    << generation_number << "\t"
    << last_generation.average_cost << "\t"
    << last_generation.best_total_cost << "\t"
    << "{" << best_genes.to_string( ", " ) << "}"
    << "\t" << (best_yet ? "BestYet" : "NoImprovement")
    << "\n\n";

  cout
    << "Generation [" << generation_number << "], "
    << "Best=" << last_generation.best_total_cost << ", "
    << "Average=" << last_generation.average_cost << ", "
    << (best_yet ? "Best generation yet" : "no improvement")
    << "\n  Best genes: {\n\t" << best_genes.to_string( "\n\t" ) << "\n}\n"
    << "Exe_time=" << last_generation.exe_time
    << endl;

  // Write HTML output for the best solution of this generation
  if( best_yet && sm_precomputed_ptr && !sm_precomputed_ptr->empty() )
  {
    try
    {
      const PeakFitForNuclideConfig config = genes_to_settings( sm_best_genes );
      const string html_filename = "nuclide_config_ga_best.html";
      const string n42_dir = "output_n42_nuclide_config_ga";
      write_results_html_and_n42( *sm_precomputed_ptr, config, sm_background_mode, html_filename, n42_dir );
      cout << "Wrote HTML results to " << html_filename << endl;
    }
    catch( const std::exception &e )
    {
      cerr << "Warning: failed to write HTML results for generation " << generation_number
           << ": " << e.what() << endl;
    }
  }//if( best_yet && sm_precomputed_ptr )

  ns_num_evals_this_generation = 0;
}//SO_report_generation


void write_results_html_and_n42(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const PeakFitForNuclideConfig &config,
  const BackgroundMode bg_mode,
  const std::string &html_filename,
  const std::string &n42_output_dir )
{
  // Open HTML file
  std::ofstream html_output( html_filename );

  try
  {
    D3SpectrumExport::write_html_page_header( html_output, "NuclideConfig GA Results", "InterSpec_resources" );
  }
  catch( std::exception &e )
  {
    cerr << "Error writing HTML header: " << e.what() << endl;
    cerr << "You probably need to symbolically link InterSpec_resources to your CWD." << endl;
    return;
  }

  html_output << "<body>" << endl;
  html_output << "<style>"
    << ".TopLinesTable{ margin-top: 25px; margin-left: auto; margin-right: auto; border-collapse: collapse; border: 1px solid black; }" << endl
    << "table, th, td{ border: 1px solid black; }" << endl
    << "fieldset{width: 90vw; margin-left: auto; margin-right: auto; margin-top: 20px;}" << endl
    << "</style>" << endl;

  size_t chart_counter = 0;
  double total_score = 0.0;

  const double num_sigma_contribution = 1.5;

  for( const PrecomputedNuclideData &pd : precomputed )
  {
    const InjectSourceInfo &src = pd.src_info->src_info;
    const string src_name = src.src_name;

    // Build options based on background mode
    Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options;
    if( bg_mode == BackgroundMode::NoBackgroundFitNorm )
      options |= FitPeaksForNuclides::FitNormBkgrndPeaks;

    try
    {
      const std::vector<std::shared_ptr<const PeakDef>> user_peaks;

      const FitPeaksForNuclides::PeakFitResult result = FitPeaksForNuclides::fit_peaks_for_nuclides(
        pd.auto_search_peaks, pd.foreground, pd.sources, user_peaks,
        pd.background, pd.drf, options, config, pd.peak_fit_prefs );

      if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        continue;

      const RelActCalcAuto::RelActAutoSolution &solution = result.solution;
      const std::vector<PeakDef> &fit_peaks = result.observable_peaks;

      // Score the results
      CombinedPeakFitScore combined_score;
      combined_score.final_fit_score = FinalFit_GA::calculate_final_fit_score(
        fit_peaks, pd.src_info->expected_signal_photopeaks, num_sigma_contribution );
      combined_score.initial_fit_weights = InitialFit_GA::calculate_peak_find_weights(
        fit_peaks, pd.src_info->expected_signal_photopeaks, num_sigma_contribution );
      combined_score.candidate_peak_score = CandidatePeak_GA::calculate_candidate_peak_score_for_source(
        fit_peaks, pd.src_info->expected_signal_photopeaks );
      CandidatePeak_GA::correct_score_for_escape_peaks(
        combined_score.candidate_peak_score, pd.src_info->expected_signal_photopeaks );

      combined_score.final_weight = combined_score.initial_fit_weights.find_weight
                                  + combined_score.final_fit_score.total_weight
                                  + combined_score.candidate_peak_score.score;
      total_score += combined_score.final_weight;

      // Write N42 file
      if( !n42_output_dir.empty() )
      {
        string outdir = n42_output_dir;
        if( !SpecUtils::is_directory( outdir ) )
          SpecUtils::create_directory( outdir );

        outdir = SpecUtils::append_path( outdir, pd.src_info->detector_name );
        if( !SpecUtils::is_directory( outdir ) )
          SpecUtils::create_directory( outdir );

        const string out_n42 = SpecUtils::append_path( outdir, src_name ) + "_nuclide_config_ga.n42";

        SpecMeas output;
        output.remove_measurements( output.measurements() );

        std::shared_ptr<SpecUtils::Measurement> out_fg = std::make_shared<SpecUtils::Measurement>( *pd.foreground );
        out_fg->set_sample_number( 1 );
        out_fg->set_title( src_name + " - Foreground" );
        output.add_measurement( out_fg, false );

        if( pd.background )
        {
          std::shared_ptr<SpecUtils::Measurement> out_bg = std::make_shared<SpecUtils::Measurement>( *pd.background );
          out_bg->set_sample_number( 2 );
          out_bg->set_title( src_name + " - Background" );
          output.add_measurement( out_bg, false );
        }

        std::deque<std::shared_ptr<const PeakDef>> peaks;
        for( const PeakDef &p : fit_peaks )
          peaks.push_back( std::make_shared<PeakDef>( p ) );
        output.setPeaks( peaks, {1} );

        output.save2012N42File( out_n42, [&](){ cerr << "Failed to write '" << out_n42 << "'" << endl; } );
      }//if( !n42_output_dir.empty() )

      // Write HTML chart
      {
        // Generate reference lines for nuclides
        std::map<std::string, std::string> reference_lines_json;

        const std::vector<Wt::WColor> colors = {
          Wt::WColor( 0, 0, 139 ),
          Wt::WColor( 0, 139, 139 ),
          Wt::WColor( 0, 100, 0 ),
          Wt::WColor( 139, 0, 139 ),
          Wt::WColor( 255, 140, 0 ),
          Wt::WColor( 139, 0, 0 ),
        };

        size_t color_index = 0;
        for( size_t rel_eff_index = 0; rel_eff_index < solution.m_rel_activities.size(); ++rel_eff_index )
        {
          for( const RelActCalcAuto::NuclideRelAct &nuc_info : solution.m_rel_activities[rel_eff_index] )
          {
            RefLineInput input;
            input.m_input_txt = nuc_info.name();
            if( nuc_info.age_was_fit && nuc_info.age > 0 )
              input.m_age = std::to_string( nuc_info.age ) + " s";
            input.m_color = colors[color_index % colors.size()];
            color_index++;
            input.m_showGammas = true;
            input.m_showXrays = true;
            input.m_showAlphas = false;
            input.m_showBetas = false;
            input.m_promptLinesOnly = false;
            input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
            if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string ref_json;
              ref_info->toJson( ref_json );
              reference_lines_json[nuc_info.name()] = ref_json;
            }
          }//for( nuc_info )
        }//for( rel_eff_index )

        // Add red reference lines for missing peaks
        const std::vector<ExpectedPhotopeakInfo> &not_detected = combined_score.candidate_peak_score.def_expected_but_not_detected;
        for( size_t idx = 0; idx < not_detected.size(); ++idx )
        {
          RefLineInput red_input;
          red_input.m_input_txt = std::to_string( not_detected[idx].effective_energy ) + " keV";
          red_input.m_color = Wt::WColor( 255, 0, 0 );
          red_input.m_showGammas = true;
          red_input.m_showXrays = false;
          red_input.m_showAlphas = false;
          red_input.m_showBetas = false;
          red_input.m_promptLinesOnly = false;
          red_input.m_lower_br_cutt_off = 0.0;

          std::shared_ptr<ReferenceLineInfo> red_ref = ReferenceLineInfo::generateRefLineInfo( red_input );
          if( red_ref && red_ref->m_validity == ReferenceLineInfo::InputValidity::Valid )
          {
            std::string json;
            red_ref->toJson( json );
            reference_lines_json["Missing Peak " + std::to_string( idx + 1 )] = json;
          }
        }//for( not_detected )

        // Chart options
        float xMin = 0.0f, xMax = 3000.0f;
        if( !solution.m_final_roi_ranges.empty() )
        {
          xMin = static_cast<float>( solution.m_final_roi_ranges.front().lower_energy );
          xMax = static_cast<float>( solution.m_final_roi_ranges.back().upper_energy );
        }
        else if( pd.foreground )
        {
          xMin = static_cast<float>( pd.foreground->gamma_energy_min() );
          xMax = static_cast<float>( pd.foreground->gamma_energy_max() );
        }

        const string title = src_name + " - NuclideConfig GA Fit";
        D3SpectrumExport::D3SpectrumChartOptions chart_options(
          title, "Energy (keV)", "Counts/Channel",
          "", true, false, false, true, true,
          false, false, false, false, false, false,
          false, false, false, xMin, xMax, reference_lines_json );

        D3SpectrumExport::D3SpectrumOptions fg_opts;
        fg_opts.line_color = "black";
        fg_opts.title = src_name;
        fg_opts.display_scale_factor = 1.0;
        fg_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;

        vector<shared_ptr<const PeakDef>> fit_peaks_ptrs;
        for( const PeakDef &p : fit_peaks )
          fit_peaks_ptrs.push_back( make_shared<PeakDef>( p ) );
        fg_opts.peaks_json = PeakDef::peak_json( fit_peaks_ptrs, pd.foreground, Wt::WColor(), false );

        const string div_id = "chart_" + std::to_string( chart_counter );
        chart_counter++;

        html_output << "<fieldset>" << endl
          << "<legend>" << src_name << " (chi2/dof=" << solution.m_chi2 << "/" << solution.m_dof << ")</legend>" << endl;
        html_output << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>" << endl;
        html_output << "<script>" << endl;

        D3SpectrumExport::write_js_for_chart( html_output, div_id, chart_options.m_dataTitle, chart_options.m_xAxisTitle, chart_options.m_yAxisTitle );

        std::vector<std::pair<const SpecUtils::Measurement *, D3SpectrumExport::D3SpectrumOptions>> measurements;
        measurements.emplace_back( pd.foreground.get(), fg_opts );

        if( pd.background )
        {
          D3SpectrumExport::D3SpectrumOptions bg_opts;
          bg_opts.line_color = "steelblue";
          bg_opts.title = "Background";
          bg_opts.display_scale_factor = pd.foreground->live_time() / pd.background->live_time();
          bg_opts.spectrum_type = SpecUtils::SpectrumType::Background;
          measurements.emplace_back( pd.background.get(), bg_opts );
        }

        D3SpectrumExport::write_and_set_data_for_chart( html_output, div_id, measurements );

        html_output << R"delim(
const resizeChart)delim" << chart_counter - 1 << R"delim( = function(){
  let height = window.innerHeight;
  let width = document.documentElement.clientWidth;
  let el = spec_chart_)delim" << div_id << R"delim(.chart;
  el.style.width = 0.8*width + "px";
  el.style.height = Math.min(500,Math.max(250, Math.min(0.4*width,height-175))) + "px";
  el.style.marginLeft = 0.05*width + "px";
  el.style.marginRight = 0.05*width + "px";
  )delim"
          << "  spec_chart_" << div_id << R"delim(.handleResize();
};

window.addEventListener('resize', resizeChart)delim" << chart_counter - 1 << R"delim();
)delim" << endl;

        D3SpectrumExport::write_set_options_for_chart( html_output, div_id, chart_options );
        html_output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );" << endl;
        html_output << "resizeChart" << chart_counter - 1 << "();" << endl;
        html_output << "</script>" << endl;

        // Score summary
        html_output << "<div style=\"margin: 10px auto; max-width: 800px;\">" << endl;
        html_output << "<h4>Score: " << combined_score.final_weight << "</h4>" << endl;
        html_output << "<p>Find weight: " << combined_score.initial_fit_weights.find_weight
                   << ", Final fit: " << combined_score.final_fit_score.total_weight
                   << ", Candidate: " << combined_score.candidate_peak_score.score << "</p>" << endl;
        html_output << "<p>Missing peaks: " << combined_score.candidate_peak_score.num_def_wanted_not_found
                   << ", Extra peaks: " << combined_score.candidate_peak_score.num_extra_peaks << "</p>" << endl;
        html_output << "</div>" << endl;

        html_output << "</fieldset>" << endl;
      }// HTML chart block

    }
    catch( const std::exception &e )
    {
      cerr << "Error processing " << src_name << ": " << e.what() << endl;
    }
  }//for( const PrecomputedNuclideData &pd : precomputed )

  html_output << "<h2 style=\"text-align: center;\">Total Score: " << total_score << "</h2>" << endl;
  html_output << "</body>" << endl;
  html_output << "</html>" << endl;
  html_output.close();
}//write_results_html_and_n42


PeakFitForNuclideConfig do_nuclide_config_ga(
  const std::vector<PrecomputedNuclideData> &precomputed,
  std::function<double( const PeakFitForNuclideConfig & )> ga_eval_fcn )
{
  assert( !sm_has_been_called );
  if( sm_has_been_called )
  {
    cerr << "You may only call NuclideConfig_GA::do_nuclide_config_ga(...) once per program execution." << endl;
    exit( 1 );
  }

  sm_has_been_called = true;

  assert( !!ga_eval_fcn );
  if( !ga_eval_fcn )
    throw runtime_error( "Invalid eval function passed in." );

  ns_ga_eval_fcn = ga_eval_fcn;
  sm_precomputed_ptr = &precomputed;

  sm_output_file.open( "nuclide_config_ga_results.txt" );
  sm_output_file << "step" << "\t" << "cost_avg" << "\t" << "cost_best" << "\t" << "solution_best" << "\n";

  EA::Chronometer timer;
  timer.tic();

  GA_Type ga_obj;
  ga_obj.problem_mode = EA::GA_MODE::SOGA;
  ga_obj.multi_threading = true;
  ga_obj.idle_delay_us = 1;
  ga_obj.dynamic_threading = (PeakFitImprove::sm_ga_population > PeakFitImprove::sm_num_optimization_threads);
  ga_obj.verbose = false;
  ga_obj.population = static_cast<unsigned int>( PeakFitImprove::sm_ga_population );
  ga_obj.N_threads = static_cast<int>( PeakFitImprove::sm_num_optimization_threads );
  ga_obj.generation_max = static_cast<int>( PeakFitImprove::sm_ga_generation_max );
  ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;
  ga_obj.init_genes = init_genes;
  ga_obj.eval_solution = eval_solution;
  ga_obj.mutate = mutate;
  ga_obj.crossover = crossover;
  ga_obj.SO_report_generation = SO_report_generation;
  ga_obj.crossover_fraction = PeakFitImprove::sm_ga_crossover_fraction;
  ga_obj.mutation_rate = PeakFitImprove::sm_ga_mutation_rate;
  ga_obj.best_stall_max = static_cast<int>( PeakFitImprove::sm_ga_best_stall_max );
  ga_obj.elite_count = static_cast<int>( PeakFitImprove::sm_ga_elite_count );

  const EA::StopReason stop_reason = ga_obj.solve();

  cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;
  cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason ) << endl;

  sm_output_file.close();
  sm_precomputed_ptr = nullptr;

  // Print the best config in copy-pasteable format
  const PeakFitForNuclideConfig best_config = genes_to_settings( sm_best_genes );
  cout << "\n========================================" << endl;
  cout << "Best NuclideConfig GA result (cost=" << sm_best_total_cost << "):" << endl;
  cout << "Genes:\n\t" << sm_best_genes.to_string( "\n\t" ) << endl;
  cout << "========================================\n" << endl;

  // Write final HTML with best solution
  try
  {
    write_results_html_and_n42( precomputed, best_config, sm_background_mode,
                                "nuclide_config_ga_final.html", "output_n42_nuclide_config_ga_final" );
    cout << "Wrote final HTML results to nuclide_config_ga_final.html" << endl;
  }
  catch( const std::exception &e )
  {
    cerr << "Warning: failed to write final HTML results: " << e.what() << endl;
  }

  return best_config;
}//do_nuclide_config_ga

}//namespace NuclideConfig_GA
