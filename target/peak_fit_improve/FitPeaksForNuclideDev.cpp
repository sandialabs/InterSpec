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

#include "InterSpec_config.h"

#include <chrono>
#include <string>
#include <fstream>
#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp> // For boost::asio::post

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "PeakFitImprove.h"
#include "InitialFit_GA.h"
#include "CandidatePeak_GA.h"

using namespace std;
namespace FitPeaksForNuclideDev
{

// This is development playground for `AnalystChecks::fit_peaks_for_nuclides(...)`
void eval_peaks_for_nuclide( const std::vector<DataSrcInfo> &srcs_info )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    throw runtime_error( "Failed to open SandiaDecayDataBase" );

  for( const DataSrcInfo &src_info : srcs_info )
  {
    const InjectSourceInfo &src = src_info.src_info;
    assert( src.src_spectra.size() >= 12 );
    if( src.src_spectra.size() < 12 )
      throw runtime_error( "Unexpected number of measurements." );

    string src_name = src.src_name;
    const auto underscore_pos = src_name.find( '_' );
    if( underscore_pos != string::npos )
    {
      const string last_part = src_name.substr(underscore_pos);
      assert( (last_part == "_Sh")  || (last_part == "_Sh-Point") || (last_part == "_Sh_BPE") || (last_part == "_Sh_PE")
             || (last_part == "_Sh_100g") || (last_part == "_ShHeavy") || (last_part == "_Unsh")
             || (last_part == "_Unsh-Point") || (last_part == "_Phantom")
             || (last_part == "_Unsh_0200") || (last_part == "_Unsh_0300") || (last_part == "_Unsh_0400")
             || (last_part == "_Unsh_0500") || (last_part == "_Unsh_1000") || (last_part == "_Unsh_2000")
             || (last_part == "_Unsh_4000") || (last_part == "_Unsh_5000") || (last_part == "_Unsh_6000")
             || (last_part == "_Unsh_8000") || (last_part == "_Unsh_9000") || (last_part == "_Unsh_9330")
      );
      src_name = src_name.substr(0,underscore_pos);
    }//if( underscore_pos != string::npos )

    vector<const SandiaDecay::Nuclide *> source_nuclides;
    if( src_name == "Tl201woTl202" )
    {
      source_nuclides.push_back( db->nuclide("Tl201"));
    }else if( src_name == "Tl201wTl202" )
    {
      source_nuclides.push_back( db->nuclide("Tl201"));
      source_nuclides.push_back( db->nuclide("Tl202"));
    }else if( src_name == "Uore" )
    {
      source_nuclides.push_back( db->nuclide("U235") );
      source_nuclides.push_back( db->nuclide("U238") );
      source_nuclides.push_back( db->nuclide("Ra226") );
    }else if( src_name == "Am241Li" )
    {
      continue;
    }else
    {
      const SandiaDecay::Nuclide * const nuc = db->nuclide( src_name );
      if( !nuc )
      {
        cout << "Unable to get nuclide for '" << src_name << "' - so will skip" << endl;
        continue;
      }
      source_nuclides.push_back( nuc );
      assert( nuc->symbol == src_name );
    }

    assert( !source_nuclides.empty() && source_nuclides.front() );
    if( source_nuclides.empty() || !source_nuclides.front() )
      throw runtime_error( "Failed to get a source" );

    const shared_ptr<SpecUtils::SpecFile> &spec_file = src.spec_file;
    const shared_ptr<const SpecUtils::Measurement> &foreground = src.src_spectra.front(); // TODO: we cold loop over all 16 of these histograms
    const shared_ptr<const SpecUtils::Measurement> &long_background = src.long_background;

    const bool isHPGe = true;
    const bool singleThreaded = false;
    shared_ptr<const DetectorPeakResponse> drf = nullptr; //Probably

    try
    {
      shared_ptr<const deque<shared_ptr<const PeakDef>>> dummy_origpeaks;

      const vector<shared_ptr<const PeakDef> > auto_search_peaks
          = ExperimentalAutomatedPeakSearch::search_for_peaks( foreground, drf, dummy_origpeaks, singleThreaded, isHPGe );


      // Create base nuclide info that will be reused across all attempts
      std::vector<RelActCalcAuto::NucInputInfo> base_nuclides;
      SandiaDecay::NuclideMixture mix;

      for( const SandiaDecay::Nuclide *nuc : source_nuclides)
      {
        RelActCalcAuto::NucInputInfo nuc_info;
        nuc_info.source = nuc;
        nuc_info.age = PeakDef::defaultDecayTime(nuc);
        nuc_info.fit_age = false; // Don't fit age by default

        base_nuclides.push_back(nuc_info);
        mix.addAgedNuclideByActivity( nuc, 1.0E-6 * SandiaDecay::curie, nuc_info.age );
      }

      double highest_energy_gamma = 0.0, lowest_energy_gamma = std::numeric_limits<double>::max();
      const vector<SandiaDecay::EnergyRatePair> photons = mix.photons( 0.0 );
      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        highest_energy_gamma = (std::max)( highest_energy_gamma, photon.energy );
        lowest_energy_gamma = (std::min)( lowest_energy_gamma, photon.energy );
      }

      lowest_energy_gamma = (std::max)( lowest_energy_gamma - (isHPGe ? 5 : 25), (double)foreground->gamma_energy_min() );
      highest_energy_gamma = (std::min)( highest_energy_gamma + (isHPGe ? 5 : 25), (double)foreground->gamma_energy_max() );



      // TODO: Implement the ROI breaking up logic here that `RelActCalcAuto::solve` and `RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres` will eventually be udpated to
      throw runtime_error( "TODO: Implement the ROI breaking up logic here that `RelActCalcAuto::solve` and `RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres` will eventually be udpated to" );
      // I think first round should be based on peaks in the spectrum...
      // Next round use RelEff fit




      // Set up base ROI to cover the full spectrum energy range
      RelActCalcAuto::RoiRange base_roi;
      base_roi.lower_energy = lowest_energy_gamma;
      base_roi.upper_energy = highest_energy_gamma;
      base_roi.continuum_type = PeakContinuum::OffsetType::Linear;  //TODO: perhaps `PeakContinuum::OffsetType::External` should be "auto"
      base_roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::CanBeBrokenUp;

      // Create all trial options combinations to run in parallel
      std::vector<RelActCalcAuto::Options> trial_options;

      // Define skew types to try
      std::vector<PeakDef::SkewType> skew_types = {
        PeakDef::SkewType::NoSkew,
        PeakDef::SkewType::GaussExp
      };

      for (PeakDef::SkewType skew_type : skew_types) {
        // LnX equation form
        {
          RelActCalcAuto::Options solve_options;

          // Create relative efficiency curve with LnX
          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "LnX Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
          rel_eff_curve.rel_eff_eqn_order = 3;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;

          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois.push_back(base_roi);

          // Set other options
          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = RelActCalcAuto::FwhmForm::Gadras;
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;

          trial_options.push_back(std::move(solve_options));
        }

        // FramPhysicalModel with Aluminum shielding
        {
          RelActCalcAuto::Options solve_options;

          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "Physical Model Al Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;
          rel_eff_curve.phys_model_use_hoerl = true;

          // Create aluminum external shielding
          auto al_shield = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
          al_shield->atomic_number = 13.0; // Aluminum
          al_shield->fit_atomic_number = false;
          al_shield->areal_density = 0.1; // Starting value in g/cm2
          al_shield->fit_areal_density = true;
          al_shield->lower_fit_areal_density = 0.001;
          al_shield->upper_fit_areal_density = 10.0;

          rel_eff_curve.phys_model_external_atten.push_back(al_shield);

          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois.push_back(base_roi);

          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = RelActCalcAuto::FwhmForm::Gadras; //RelActCalcAuto::FwhmForm::Polynomial_3
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;

          trial_options.push_back(std::move(solve_options));
        }

        // FramPhysicalModel with Lead shielding
        {
          RelActCalcAuto::Options solve_options;

          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "Physical Model Pb Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;
          rel_eff_curve.phys_model_use_hoerl = true;

          // Create lead external shielding
          auto pb_shield = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
          pb_shield->atomic_number = 82.0; // Lead
          pb_shield->fit_atomic_number = false;
          pb_shield->areal_density = 0.1; // Starting value in g/cm2
          pb_shield->fit_areal_density = true;
          pb_shield->lower_fit_areal_density = 0.0;
          pb_shield->upper_fit_areal_density = 10.0;

          rel_eff_curve.phys_model_external_atten.push_back(pb_shield);

          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois.push_back(base_roi);

          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = RelActCalcAuto::FwhmForm::Gadras;
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;

          trial_options.push_back(std::move(solve_options));
        }
      }

      // Execute all trials in parallel using ThreadPool
      std::vector<RelActCalcAuto::RelActAutoSolution> solutions(trial_options.size());
      std::mutex solutions_mutex;

      SpecUtilsAsync::ThreadPool pool;
      for (size_t i = 0; i < trial_options.size(); ++i) {
        pool.post([&, i]() {
          try {
            RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
              trial_options[i], foreground, long_background, drf, auto_search_peaks
            );

            std::lock_guard<std::mutex> lock(solutions_mutex);
            solutions[i] = std::move(solution);
          } catch (const std::exception &e) {
            std::lock_guard<std::mutex> lock(solutions_mutex);
            solutions[i].m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
            solutions[i].m_error_message = e.what();
          }
        });
      }

      pool.join();

      // Find the best solution from all parallel trials
      RelActCalcAuto::RelActAutoSolution best_solution;
      double best_chi2 = std::numeric_limits<double>::max();
      bool found_valid_solution = false;
      string last_error_message;

      for (size_t i = 0; i < solutions.size(); ++i) {
        const RelActCalcAuto::RelActAutoSolution &solution = solutions[i];

        if (solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success &&
            solution.m_chi2 < best_chi2) {
          best_chi2 = solution.m_chi2;
          best_solution = solution;
          found_valid_solution = true;
        } else if (solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success) {
          last_error_message = solution.m_error_message;
          cerr << "RelActCalcAuto::solve(trial " << i << "): failed: " << solution.m_error_message << endl;
        }
      }

      if (!found_valid_solution) {
        throw std::runtime_error("RelActCalcAuto::solve failed for all attempted configurations ('" + last_error_message + "')");
      }

      // Use the best solution
      RelActCalcAuto::RelActAutoSolution solution = std::move(best_solution);

      //solution.m_status != RelActAutoSolution::Status::Success
      //solution.m_chi2 / solution.m_dof
    } catch (const std::exception &e) {
      //throw std::runtime_error("Error in fit_peaks_for_nuclides: " + string(e.what()));
      cerr << e.what() << endl;
    }
  }//for( const DataSrcInfo &src_info : srcs_info )
}//void eval_peaks_for_nuclide( const DataSrcInfo &src_info )


}//namespace FitPeaksForNuclideDev
