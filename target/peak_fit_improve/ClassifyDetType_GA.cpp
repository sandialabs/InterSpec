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

#include <map>
#include <cmath>
#include <mutex>
#include <random>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#define OPENGA_EXTERN_LOCAL_VARS 1

#include "openGA.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitSpecImp.h"

#include "PeakFitImprove.h"
#include "CandidatePeak_GA.h"
#include "ClassifyDetType_GA.h"

using namespace std;


// ============================================================================
// Anonymous namespace GA state variables
// ============================================================================
namespace
{
  // Stage 1 (HighRes: High vs NotHigh) GA state
  std::function<double( const HighResClassifySettings & )> ns_highres_eval_fcn;
  std::ofstream sm_highres_output_file;
  bool sm_highres_set_best_genes = false;
  ClassifyDetType_GA::HighResSolution sm_highres_best_genes;
  double sm_highres_best_total_cost = 1.0E99;
  bool sm_highres_has_been_called = false;

}//anonymous namespace


// ============================================================================
// HighResClassifySettings print/to_json
// ============================================================================

string HighResClassifySettings::print( const string &var_name ) const
{
  string result;
  result += var_name + ".hpge_fwhm_bias_mult = " + to_string( hpge_fwhm_bias_mult ) + ";\n";
  result += var_name + ".nai_fwhm_bias_mult = " + to_string( nai_fwhm_bias_mult ) + ";\n";
  result += var_name + ".hpge_distance_weight = " + to_string( hpge_distance_weight ) + ";\n";
  result += var_name + ".nai_distance_weight = " + to_string( nai_distance_weight ) + ";\n";
  result += var_name + ".amp_clamp_denom = " + to_string( amp_clamp_denom ) + ";\n";
  result += var_name + ".weight_clamp_max = " + to_string( weight_clamp_max ) + ";\n";
  result += var_name + ".min_peak_energy = " + to_string( min_peak_energy ) + ";\n";
  result += var_name + ".max_peak_energy = " + to_string( max_peak_energy ) + ";\n";
  result += var_name + ".unknown_threshold = " + to_string( unknown_threshold ) + ";\n";
  result += var_name + ".min_peaks_for_classify = " + to_string( min_peaks_for_classify ) + ";\n";
  result += var_name + ".min_peak_significance = " + to_string( min_peak_significance ) + ";\n";
  result += var_name + ".sig_fraction_of_max = " + to_string( sig_fraction_of_max ) + ";\n";
  result += var_name + ".narrow_penalty_mult = " + to_string( narrow_penalty_mult ) + ";\n";
  result += var_name + ".narrow_sig_threshold = " + to_string( narrow_sig_threshold ) + ";\n";
  result += var_name + ".best_peak_czt_penalty = " + to_string( best_peak_czt_penalty ) + ";\n";
  result += var_name + ".best_peak_conf_threshold = " + to_string( best_peak_conf_threshold ) + ";\n";
  return result;
}//HighResClassifySettings::print


string HighResClassifySettings::to_json() const
{
  string result = "{\n";
  result += "  \"hpge_fwhm_bias_mult\": " + to_string( hpge_fwhm_bias_mult ) + ",\n";
  result += "  \"nai_fwhm_bias_mult\": " + to_string( nai_fwhm_bias_mult ) + ",\n";
  result += "  \"hpge_distance_weight\": " + to_string( hpge_distance_weight ) + ",\n";
  result += "  \"nai_distance_weight\": " + to_string( nai_distance_weight ) + ",\n";
  result += "  \"amp_clamp_denom\": " + to_string( amp_clamp_denom ) + ",\n";
  result += "  \"weight_clamp_max\": " + to_string( weight_clamp_max ) + ",\n";
  result += "  \"min_peak_energy\": " + to_string( min_peak_energy ) + ",\n";
  result += "  \"max_peak_energy\": " + to_string( max_peak_energy ) + ",\n";
  result += "  \"unknown_threshold\": " + to_string( unknown_threshold ) + ",\n";
  result += "  \"min_peaks_for_classify\": " + to_string( min_peaks_for_classify ) + ",\n";
  result += "  \"min_peak_significance\": " + to_string( min_peak_significance ) + ",\n";
  result += "  \"sig_fraction_of_max\": " + to_string( sig_fraction_of_max ) + ",\n";
  result += "  \"narrow_penalty_mult\": " + to_string( narrow_penalty_mult ) + ",\n";
  result += "  \"narrow_sig_threshold\": " + to_string( narrow_sig_threshold ) + ",\n";
  result += "  \"best_peak_czt_penalty\": " + to_string( best_peak_czt_penalty ) + ",\n";
  result += "  \"best_peak_conf_threshold\": " + to_string( best_peak_conf_threshold ) + "\n";
  result += "}";
  return result;
}//HighResClassifySettings::to_json


// ============================================================================
// Detector truth mapping and utilities
// ============================================================================

namespace ClassifyDetType_GA
{

PeakFitUtils::CoarseResolutionType true_det_type_for_name( const string &detector_name )
{
  using PeakFitUtils::CoarseResolutionType;

  // HPGe detectors
  if( SpecUtils::icontains( detector_name, "Falcon" )
     || SpecUtils::icontains( detector_name, "Fulcrum" )
     || SpecUtils::icontains( detector_name, "LANL_X" )
     || SpecUtils::icontains( detector_name, "Detective" )
     || SpecUtils::icontains( detector_name, "HPGe_Planar" ) )
  {
    return CoarseResolutionType::High;
  }

  // CZT detectors
  if( SpecUtils::icontains( detector_name, "1.5cm-2cm" )
     || SpecUtils::icontains( detector_name, "1cm-1cm" )
     || SpecUtils::icontains( detector_name, "CZT_H3D" )
     || SpecUtils::icontains( detector_name, "Kromek-GR1" )
     || SpecUtils::icontains( detector_name, "nanoRaider" )
     || SpecUtils::icontains( detector_name, "Interceptor" )
     || SpecUtils::icontains( detector_name, "Raider" ) )
  {
    // "Raider" must come after more specific "nanoRaider" check to avoid matching both,
    //  but since we use icontains and nanoRaider contains "Raider", both will match anyway.
    //  The result is the same for both (MedRes), so this is fine.
    return CoarseResolutionType::MedRes;
  }

  // LaBr detectors
  if( SpecUtils::icontains( detector_name, "InSpector 1000" )
     || SpecUtils::icontains( detector_name, "SAM-Eagle-LaBr" )
     || SpecUtils::icontains( detector_name, "LaBr3" )
     || SpecUtils::icontains( detector_name, "IdentiFINDER-LaBr" )
     || SpecUtils::icontains( detector_name, "Radseeker-LaBr" ) )
  {
    return CoarseResolutionType::MedRes;
  }

  // Everything else is NaI/CsI (low resolution)
  return CoarseResolutionType::Low;
}//true_det_type_for_name


const char *det_type_str( PeakFitUtils::CoarseResolutionType type )
{
  switch( type )
  {
    case PeakFitUtils::CoarseResolutionType::Low:        return "NaI/CsI";
    case PeakFitUtils::CoarseResolutionType::LaBr:       return "LaBr";
    case PeakFitUtils::CoarseResolutionType::CZT:        return "CZT";
    case PeakFitUtils::CoarseResolutionType::MedRes:     return "MedRes";
    case PeakFitUtils::CoarseResolutionType::LowOrMedRes: return "LowOrMedRes";
    case PeakFitUtils::CoarseResolutionType::High:       return "HPGe";
    case PeakFitUtils::CoarseResolutionType::Unknown:    return "Unknown";
  }
  return "???";
}//det_type_str


// ============================================================================
// Stage 1: delegates to PeakFitSpec::classify_lowres_highres_from_peak_candidates
// ============================================================================

PeakFitUtils::CoarseResolutionType classify_highres_from_peaks(
  const vector<PeakDef> &candidate_peaks,
  const HighResClassifySettings &settings,
  double *confidence_out )
{
  const PeakFitSpec::SpecClassType result
    = PeakFitSpec::classify_lowres_highres_from_peak_candidates( candidate_peaks, settings, confidence_out );
  switch( result )
  {
    case PeakFitSpec::SpecClassType::High:        return PeakFitUtils::CoarseResolutionType::High;
    case PeakFitSpec::SpecClassType::LowOrMedRes: return PeakFitUtils::CoarseResolutionType::LowOrMedRes;
    case PeakFitSpec::SpecClassType::Unknown:     return PeakFitUtils::CoarseResolutionType::Unknown;
  }
  return PeakFitUtils::CoarseResolutionType::Unknown;
}//classify_highres_from_peaks


// ============================================================================
// Precompute candidate peaks
// ============================================================================

vector<PrecomputedCandidates> precompute_candidates(
  const vector<DataSrcInfo> &input_srcs,
  const FindCandidateSettings &candidate_settings )
{
  vector<PrecomputedCandidates> result( input_srcs.size() );

  for( size_t i = 0; i < input_srcs.size(); ++i )
  {
    const DataSrcInfo &info = input_srcs[i];
    PrecomputedCandidates &pc = result[i];

    pc.truth = info.det_type;

    const size_t nchannel = info.src_info.src_spectra.empty() ? 0
                              : info.src_info.src_spectra[0]->num_gamma_channels();

    if( !info.src_info.src_spectra.empty() && info.src_info.src_spectra[0] )
    {
      pc.src_candidates = PeakFitSpec::find_candidate_peaks(
        info.src_info.src_spectra[0], 0, nchannel, candidate_settings );
    }

    if( info.src_info.short_background )
    {
      pc.bg_candidates = PeakFitSpec::find_candidate_peaks(
        info.src_info.short_background, 0, nchannel, candidate_settings );
      pc.has_background = true;
    }
  }//for( input_srcs )

  return result;
}//precompute_candidates


// ============================================================================
// Stage 1 evaluation: High vs NotHigh scoring
// ============================================================================

double eval_highres_precomputed( const HighResClassifySettings &settings,
                                 const vector<PrecomputedCandidates> &precomputed )
{
  using PeakFitUtils::CoarseResolutionType;

  double total_score = 0.0;
  size_t num_spectra = 0;

  for( const PrecomputedCandidates &pc : precomputed )
  {
    // For Stage 1: truth is either High (HPGe) or NotHigh (everything else)
    const bool is_hpge = (pc.truth == CoarseResolutionType::High);

    // Source spectrum
    {
      const CoarseResolutionType classified
        = classify_highres_from_peaks( pc.src_candidates, settings, nullptr );

      if( is_hpge )
      {
        if( classified == CoarseResolutionType::High )
          total_score += 1.0;
        else if( classified == CoarseResolutionType::Unknown )
          total_score -= 1.0;
        else
          total_score -= 3.0;
      }
      else
      {
        if( classified == CoarseResolutionType::LowOrMedRes )
          total_score += 1.0;
        else if( classified == CoarseResolutionType::Unknown )
          total_score -= 1.0;
        else
          total_score -= 3.0;
      }

      num_spectra += 1;
    }

    // Background spectrum (softer penalties)
    if( pc.has_background )
    {
      const CoarseResolutionType bg_classified
        = classify_highres_from_peaks( pc.bg_candidates, settings, nullptr );

      if( is_hpge )
      {
        if( bg_classified == CoarseResolutionType::High )
          total_score += 0.5;
        else if( bg_classified == CoarseResolutionType::Unknown )
          total_score -= 0.0;
        else
          total_score -= 1.5;
      }
      else
      {
        if( bg_classified == CoarseResolutionType::LowOrMedRes )
          total_score += 0.5;
        else if( bg_classified == CoarseResolutionType::Unknown )
          total_score -= 0.0;
        else
          total_score -= 1.5;
      }

      num_spectra += 1;
    }
  }//for( precomputed )

  if( num_spectra == 0 )
    return 0.0;

  return total_score / static_cast<double>( num_spectra );
}//eval_highres_precomputed


// ============================================================================
// Stage 1 confusion matrix and per-detector results
// ============================================================================

void print_stage1_results_precomputed(
  const vector<DataSrcInfo> &input_srcs,
  const vector<PrecomputedCandidates> &precomputed,
  const HighResClassifySettings &highres_settings )
{
  using PeakFitUtils::CoarseResolutionType;

  assert( input_srcs.size() == precomputed.size() );

  cout << "\n=== Stage 1 (HPGe vs NotHPGe) Results ===" << endl;

  // Confusion matrix: [true][predicted]
  // Rows: 0=NotHPGe, 1=HPGe
  // Cols: 0=NotHPGe, 1=HPGe, 2=Unknown
  size_t confusion[2][3] = {};

  struct DetectorStats
  {
    string name;
    bool is_hpge = false;
    size_t num_correct = 0;
    size_t num_unknown = 0;
    size_t num_incorrect = 0;
  };

  map<string, DetectorStats> detector_stats;

  for( size_t idx = 0; idx < input_srcs.size(); ++idx )
  {
    const DataSrcInfo &info = input_srcs[idx];
    const PrecomputedCandidates &pc = precomputed[idx];

    const bool is_hpge = (pc.truth == CoarseResolutionType::High);
    const size_t true_row = is_hpge ? 1 : 0;

    DetectorStats &stats = detector_stats[info.detector_name];
    stats.name = info.detector_name;
    stats.is_hpge = is_hpge;

    // Source spectrum
    {
      const CoarseResolutionType classified
        = classify_highres_from_peaks( pc.src_candidates, highres_settings, nullptr );

      size_t pred_col;
      if( classified == CoarseResolutionType::High )
        pred_col = 1;
      else if( classified == CoarseResolutionType::Unknown )
        pred_col = 2;
      else
        pred_col = 0;

      confusion[true_row][pred_col] += 1;

      const bool correct = (is_hpge && pred_col == 1) || (!is_hpge && pred_col == 0);
      if( correct )
        stats.num_correct += 1;
      else if( pred_col == 2 )
        stats.num_unknown += 1;
      else
        stats.num_incorrect += 1;
    }

    // Background
    if( pc.has_background )
    {
      const CoarseResolutionType bg_classified
        = classify_highres_from_peaks( pc.bg_candidates, highres_settings, nullptr );

      size_t pred_col;
      if( bg_classified == CoarseResolutionType::High )
        pred_col = 1;
      else if( bg_classified == CoarseResolutionType::Unknown )
        pred_col = 2;
      else
        pred_col = 0;

      confusion[true_row][pred_col] += 1;

      const bool correct = (is_hpge && pred_col == 1) || (!is_hpge && pred_col == 0);
      if( correct )
        stats.num_correct += 1;
      else if( pred_col == 2 )
        stats.num_unknown += 1;
      else
        stats.num_incorrect += 1;
    }
  }//for( input_srcs )

  // Print per-detector summary
  cout << "\n--- Per-Detector Stage 1 Summary ---" << endl;
  for( const auto &kv : detector_stats )
  {
    const DetectorStats &ds = kv.second;
    const size_t total = ds.num_correct + ds.num_unknown + ds.num_incorrect;
    cout << "  " << setw( 30 ) << left << ds.name
         << " (truth=" << (ds.is_hpge ? "HPGe" : "NotHPGe") << "): "
         << ds.num_correct << "/" << total << " correct, "
         << ds.num_unknown << " unknown, "
         << ds.num_incorrect << " incorrect" << endl;
  }

  // Print confusion matrix
  cout << "\n--- Stage 1 Confusion Matrix ---" << endl;
  const char *pred_labels[] = { "NotHPGe", "HPGe", "Unknown" };
  cout << setw( 16 ) << left << "True\\Pred";
  for( int j = 0; j < 3; ++j )
    cout << setw( 10 ) << right << pred_labels[j];
  cout << endl;

  const char *true_labels[] = { "NotHPGe", "HPGe" };
  for( int i = 0; i < 2; ++i )
  {
    cout << setw( 16 ) << left << true_labels[i];
    for( int j = 0; j < 3; ++j )
      cout << setw( 10 ) << right << confusion[i][j];
    cout << endl;
  }

  // Overall counts
  const size_t total = confusion[0][0] + confusion[0][1] + confusion[0][2]
                     + confusion[1][0] + confusion[1][1] + confusion[1][2];
  const size_t correct = confusion[0][0] + confusion[1][1];
  const size_t incorrect = confusion[0][1] + confusion[1][0];
  const size_t unknown = confusion[0][2] + confusion[1][2];

  cout << "\nStage 1 overall (" << total << " spectra): "
       << correct << " correct (" << fixed << setprecision( 1 )
       << (total > 0 ? 100.0 * correct / total : 0.0) << "%), "
       << incorrect << " incorrect (" << (total > 0 ? 100.0 * incorrect / total : 0.0) << "%), "
       << unknown << " unknown (" << (total > 0 ? 100.0 * unknown / total : 0.0) << "%)" << endl;

  // Specifically call out HPGe FN and FP
  cout << "  HPGe false negatives (HPGe->NotHPGe): " << confusion[1][0] << endl;
  cout << "  HPGe false positives (NotHPGe->HPGe): " << confusion[0][1] << endl;
  cout << "  HPGe unknown: " << confusion[1][2] << endl;
}//print_stage1_results_precomputed


// ============================================================================
// Diagnostic output for combined two-stage classifier
// ============================================================================

void print_diagnostic_output_precomputed(
  const vector<DataSrcInfo> &input_srcs,
  const vector<PrecomputedCandidates> &precomputed,
  const HighResClassifySettings &highres_settings )
{
  using PeakFitUtils::CoarseResolutionType;

  assert( input_srcs.size() == precomputed.size() );

  cout << "\n=== Detector Type Classification Diagnostics (precomputed) ===" << endl;

  // Confusion matrix: [true_type][predicted_type]
  // Index: 0=Low/LowOrMedRes, 1=MedRes, 2=High, 3=Unknown
  auto type_to_idx = []( CoarseResolutionType t ) -> size_t {
    switch( t )
    {
      case CoarseResolutionType::Low:        return 0;
      case CoarseResolutionType::LowOrMedRes: return 0;
      case CoarseResolutionType::LaBr:       return 1;  // mapped to MedRes
      case CoarseResolutionType::CZT:        return 1;  // mapped to MedRes
      case CoarseResolutionType::MedRes:     return 1;
      case CoarseResolutionType::High:       return 2;
      case CoarseResolutionType::Unknown:    return 3;
    }
    return 3;
  };

  size_t confusion[4][4] = {};

  struct DetectorStats
  {
    string name;
    CoarseResolutionType truth;
    size_t num_correct = 0;
    size_t num_unknown = 0;
    size_t num_incorrect = 0;
    double total_score = 0.0;
  };

  map<string, DetectorStats> detector_stats;
  map<string, size_t> incorrect_examples_printed;
  const size_t max_incorrect_examples_per_detector = 3;

  // Collect all HPGe-involved misclassifications for summary at end
  struct HpgeMisclass
  {
    string detector_name;
    string location_name;
    string live_time_name;
    string src_name;
    bool is_background;
    CoarseResolutionType truth;
    CoarseResolutionType predicted;
    double confidence;
    vector<PeakDef> peaks;  // the candidate peaks used for classification
  };
  vector<HpgeMisclass> hpge_misclassifications;

  for( size_t idx = 0; idx < input_srcs.size(); ++idx )
  {
    const DataSrcInfo &info = input_srcs[idx];
    const PrecomputedCandidates &pc = precomputed[idx];

    DetectorStats &stats = detector_stats[info.detector_name];
    stats.name = info.detector_name;
    stats.truth = pc.truth;

    // Classify source spectrum
    {
      double confidence = 0.0;
      const CoarseResolutionType classified
        = classify_highres_from_peaks( pc.src_candidates, highres_settings, &confidence );

      const bool is_correct = (type_to_idx( classified ) == type_to_idx( pc.truth ));
      const bool is_unknown = (classified == CoarseResolutionType::Unknown);

      // Print per-peak detail for incorrect/unknown results (limited per detector)
      if( !is_correct )
      {
        size_t &count = incorrect_examples_printed[info.detector_name];
        if( count < max_incorrect_examples_per_detector )
        {
          count += 1;

          cout << "  [" << info.detector_name << "] Source: " << info.src_info.src_name
               << " (truth=" << det_type_str( pc.truth ) << "):" << endl;

          // Show Stage 1 per-peak diagnostics (HPGe vs NaI distances)
          cout << "    --- Stage 1 (HPGe vs NaI) per-peak ---" << endl;
          size_t num_passing = 0;
          for( const PeakDef &pk : pc.src_candidates )
          {
            if( !pk.gausPeak() )
              continue;
            if( (pk.mean() < highres_settings.min_peak_energy) || (pk.mean() > highres_settings.max_peak_energy) )
              continue;
            if( pk.amplitude() < highres_settings.min_peak_significance )
              continue;
            num_passing += 1;

            const double ref_hpge = highres_settings.hpge_fwhm_bias_mult
              * PeakFitUtils::hpge_fwhm_fcn( static_cast<float>( pk.mean() ) );
            const double ref_nai = highres_settings.nai_fwhm_bias_mult
              * PeakFitUtils::nai_fwhm_fcn( static_cast<float>( pk.mean() ) );

            const double mfwhm = pk.fwhm();
            const double d_hpge = fabs( log( mfwhm / ref_hpge ) ) * highres_settings.hpge_distance_weight;
            const double d_nai  = fabs( log( mfwhm / ref_nai ) )  * highres_settings.nai_distance_weight;

            const char *vote_str = (d_hpge <= d_nai) ? "HPGe" : "NaI";

            cout << "    Peak " << fixed << setprecision( 1 ) << pk.mean() << " keV:"
                 << " FWHM=" << setprecision( 2 ) << mfwhm
                 << ", amp=" << setprecision( 1 ) << pk.amplitude()
                 << " | dists: HPGe=" << setprecision( 3 ) << d_hpge
                 << " NaI=" << d_nai
                 << " -> " << vote_str << endl;
          }//for( pc.src_candidates )

          // Also show Stage 1 result
          double s1_conf = 0.0;
          const CoarseResolutionType s1_result
            = classify_highres_from_peaks( pc.src_candidates, highres_settings, &s1_conf );
          cout << "    Stage 1 peaks passing: " << num_passing
               << ", Result: " << det_type_str( s1_result )
               << " (conf=" << setprecision( 2 ) << s1_conf << ")" << endl;

          cout << "    Final result: " << det_type_str( classified )
               << " (confidence=" << setprecision( 2 ) << confidence << ")"
               << " -> " << (is_unknown ? "UNKNOWN" : "INCORRECT") << endl;
        }
      }//if( !is_correct )

      confusion[type_to_idx( pc.truth )][type_to_idx( classified )] += 1;

      if( is_correct )
      {
        stats.num_correct += 1;
        stats.total_score += 1.0;
      }
      else if( is_unknown )
      {
        stats.num_unknown += 1;
        stats.total_score -= 0.5;
      }
      else
      {
        stats.num_incorrect += 1;
        stats.total_score -= 3.0;

        const bool involves_hpge = (pc.truth == CoarseResolutionType::High)
                                || (classified == CoarseResolutionType::High);
        if( involves_hpge )
          hpge_misclassifications.push_back( { info.detector_name, info.location_name,
            info.live_time_name, info.src_info.src_name,
            false, pc.truth, classified, confidence, pc.src_candidates } );
      }
    }

    // Background spectrum
    if( pc.has_background )
    {
      double bg_confidence = 0.0;
      const CoarseResolutionType bg_classified
        = classify_highres_from_peaks( pc.bg_candidates, highres_settings, &bg_confidence );

      const bool bg_correct = (type_to_idx( bg_classified ) == type_to_idx( pc.truth ));
      const bool bg_unknown = (bg_classified == CoarseResolutionType::Unknown);

      confusion[type_to_idx( pc.truth )][type_to_idx( bg_classified )] += 1;

      if( bg_correct )
      {
        stats.num_correct += 1;
        stats.total_score += 0.5;
      }
      else if( bg_unknown )
      {
        stats.num_unknown += 1;
        stats.total_score -= 0.0;
      }
      else
      {
        stats.num_incorrect += 1;
        stats.total_score -= 1.5;

        const bool involves_hpge = (pc.truth == CoarseResolutionType::High)
                                || (bg_classified == CoarseResolutionType::High);
        if( involves_hpge )
          hpge_misclassifications.push_back( { info.detector_name, info.location_name,
            info.live_time_name, info.src_info.src_name,
            true, pc.truth, bg_classified, bg_confidence, pc.bg_candidates } );
      }
    }
  }//for( input_srcs )

  // Print per-detector summary
  cout << "\n=== Per-Detector Summary ===" << endl;
  for( const auto &kv : detector_stats )
  {
    const DetectorStats &ds = kv.second;
    const size_t total = ds.num_correct + ds.num_unknown + ds.num_incorrect;
    cout << "  " << setw( 30 ) << left << ds.name
         << " (" << det_type_str( ds.truth ) << "): "
         << ds.num_correct << "/" << total << " correct, "
         << ds.num_unknown << " unknown, "
         << ds.num_incorrect << " incorrect, "
         << "score=" << fixed << setprecision( 1 ) << ds.total_score << endl;
  }

  // Print confusion matrix (final 3-way result: Low/MedRes/High/Unknown)
  cout << "\n=== Confusion Matrix ===" << endl;
  const char *type_labels[] = { "NaI/CsI", "MedRes", "HPGe", "Unknown" };
  cout << setw( 16 ) << left << "True\\Pred";
  for( int j = 0; j < 4; ++j )
    cout << setw( 10 ) << right << type_labels[j];
  cout << endl;

  for( int i = 0; i < 3; ++i )  // only 3 truth types (Low, MedRes, High)
  {
    cout << setw( 16 ) << left << type_labels[i];
    for( int j = 0; j < 4; ++j )
      cout << setw( 10 ) << right << confusion[i][j];
    cout << endl;
  }

  // Overall score and classification counts
  double total_score = 0.0;
  size_t total_spectra = 0;
  size_t total_correct = 0, total_incorrect = 0, total_unknown = 0;
  for( const auto &kv : detector_stats )
  {
    total_score += kv.second.total_score;
    total_correct += kv.second.num_correct;
    total_incorrect += kv.second.num_incorrect;
    total_unknown += kv.second.num_unknown;
    total_spectra += kv.second.num_correct + kv.second.num_unknown + kv.second.num_incorrect;
  }

  const double pct_correct   = total_spectra > 0 ? 100.0 * total_correct   / total_spectra : 0.0;
  const double pct_incorrect = total_spectra > 0 ? 100.0 * total_incorrect / total_spectra : 0.0;
  const double pct_unknown   = total_spectra > 0 ? 100.0 * total_unknown   / total_spectra : 0.0;

  cout << "\nOverall score: " << fixed << setprecision( 3 )
       << (total_spectra > 0 ? total_score / total_spectra : 0.0)
       << " (" << total_spectra << " spectra)"
       << " | Correct: " << total_correct << " (" << setprecision( 1 ) << pct_correct << "%)"
       << ", Incorrect: " << total_incorrect << " (" << pct_incorrect << "%)"
       << ", Unknown: " << total_unknown << " (" << pct_unknown << "%)" << endl;

  // Print all HPGe-involved misclassifications with peak details
  if( !hpge_misclassifications.empty() )
  {
    cout << "\n=== HPGe-Involved Misclassifications (" << hpge_misclassifications.size() << " total) ===" << endl;
    for( const HpgeMisclass &mc : hpge_misclassifications )
    {
      cout << "  " << setw( 30 ) << left << mc.detector_name
           << " " << (mc.is_background ? "[BG] " : "[Src]")
           << " truth=" << setw( 8 ) << det_type_str( mc.truth )
           << " predicted=" << setw( 8 ) << det_type_str( mc.predicted )
           << " conf=" << fixed << setprecision( 3 ) << mc.confidence
           << " " << mc.location_name << "/" << mc.live_time_name
           << "/" << mc.src_name << endl;

      // Print candidate peaks with reference FWHM values
      for( const PeakDef &pk : mc.peaks )
      {
        if( !pk.gausPeak() )
          continue;
        const float e = static_cast<float>( pk.mean() );
        const double pk_sig = std::sqrt( std::max( pk.peakArea(), 0.0 ) );
        cout << "      peak: mean=" << fixed << setprecision( 1 ) << pk.mean()
             << " keV, area=" << setprecision( 1 ) << pk.peakArea()
             << ", sig=" << setprecision( 1 ) << pk_sig
             << ", fwhm=" << setprecision( 2 ) << pk.fwhm()
             << "  (ref: hpge=" << setprecision( 2 ) << PeakFitUtils::hpge_fwhm_fcn( e )
             << ", czt=" << PeakFitUtils::czt_fwhm_fcn( e )
             << ", labr=" << PeakFitUtils::labr_fwhm_fcn( e )
             << ", nai=" << PeakFitUtils::nai_fwhm_fcn( e ) << ")"
             << endl;
      }
    }
  }
}//print_diagnostic_output_precomputed


// ============================================================================
// Stage 1 GA functions (HighRes: High vs NotHigh)
// ============================================================================

string HighResSolution::to_string( const string &separator ) const
{
  return
    string( "hpge_fwhm_bias_mult: " ) + std::to_string( hpge_fwhm_bias_mult )
    + separator + "nai_fwhm_bias_mult: " + std::to_string( nai_fwhm_bias_mult )
    + separator + "hpge_distance_weight: " + std::to_string( hpge_distance_weight )
    + separator + "nai_distance_weight: " + std::to_string( nai_distance_weight )
    + separator + "amp_clamp_denom: " + std::to_string( amp_clamp_denom )
    + separator + "weight_clamp_max: " + std::to_string( weight_clamp_max )
    + separator + "min_peak_energy: " + std::to_string( min_peak_energy )
    + separator + "max_peak_energy: " + std::to_string( max_peak_energy )
    + separator + "unknown_threshold: " + std::to_string( unknown_threshold )
    + separator + "min_peaks_for_classify: " + std::to_string( min_peaks_for_classify )
    + separator + "min_peak_significance: " + std::to_string( min_peak_significance )
    + separator + "sig_fraction_of_max: " + std::to_string( sig_fraction_of_max )
    + separator + "narrow_penalty_mult: " + std::to_string( narrow_penalty_mult )
    + separator + "narrow_sig_threshold: " + std::to_string( narrow_sig_threshold )
    + separator + "best_peak_czt_penalty: " + std::to_string( best_peak_czt_penalty )
    + separator + "best_peak_conf_threshold: " + std::to_string( best_peak_conf_threshold );
}//HighResSolution::to_string


void highres_init_genes( HighResSolution &p, const std::function<double(void)> &rnd01 )
{
  p.hpge_fwhm_bias_mult = 0.4 + 0.9 * rnd01();     // [0.4, 1.3]
  p.nai_fwhm_bias_mult  = 0.7 + 0.8 * rnd01();     // [0.7, 1.5]

  p.hpge_distance_weight = 0.1 + 4.9 * rnd01();    // [0.1, 5.0]
  p.nai_distance_weight  = 0.1 + 4.9 * rnd01();    // [0.1, 5.0]

  p.amp_clamp_denom  = 1.0 + 49.0 * rnd01();       // [1, 50]
  p.weight_clamp_max = 1.0 + 29.0 * rnd01();       // [1, 30]

  p.min_peak_energy = 30.0 + 70.0 * rnd01();       // [30, 100]
  p.max_peak_energy = 2000.0 + 3000.0 * rnd01();   // [2000, 5000]

  p.unknown_threshold = 0.15 + 0.65 * rnd01();     // [0.15, 0.8]

  p.min_peaks_for_classify = 1 + static_cast<int>( 4 * rnd01() ); // [1, 4]

  p.min_peak_significance = 2.25 + 2.75 * rnd01(); // [2.25, 5.0]

  p.sig_fraction_of_max   = 0.005 + 0.295 * rnd01();  // [0.005, 0.3]
  p.narrow_penalty_mult   = 0.5 + 9.5 * rnd01();      // [0.5, 10.0]
  p.narrow_sig_threshold  = 5.0 + 95.0 * rnd01();     // [5, 100]
  p.best_peak_czt_penalty      = 5.0 * rnd01();       // [0, 5]
  p.best_peak_conf_threshold   = 0.5 + 0.5 * rnd01(); // [0.5, 1.0]
}//highres_init_genes


HighResClassifySettings highres_genes_to_settings( const HighResSolution &p )
{
  HighResClassifySettings s;
  s.hpge_fwhm_bias_mult  = p.hpge_fwhm_bias_mult;
  s.nai_fwhm_bias_mult   = p.nai_fwhm_bias_mult;
  s.hpge_distance_weight = p.hpge_distance_weight;
  s.nai_distance_weight  = p.nai_distance_weight;
  s.amp_clamp_denom        = p.amp_clamp_denom;
  s.weight_clamp_max       = p.weight_clamp_max;
  s.min_peak_energy        = p.min_peak_energy;
  s.max_peak_energy        = p.max_peak_energy;
  s.unknown_threshold      = p.unknown_threshold;
  s.min_peaks_for_classify = p.min_peaks_for_classify;
  s.min_peak_significance  = p.min_peak_significance;
  s.sig_fraction_of_max    = p.sig_fraction_of_max;
  s.narrow_penalty_mult    = p.narrow_penalty_mult;
  s.narrow_sig_threshold   = p.narrow_sig_threshold;
  s.best_peak_czt_penalty      = p.best_peak_czt_penalty;
  s.best_peak_conf_threshold   = p.best_peak_conf_threshold;
  return s;
}//highres_genes_to_settings


bool highres_eval_solution( const HighResSolution &p, HighResCost &c )
{
  const HighResClassifySettings settings = highres_genes_to_settings( p );

  assert( ns_highres_eval_fcn );

  try
  {
    c.objective1 = -1.0 * ns_highres_eval_fcn( settings );
  }
  catch( std::exception & )
  {
    return false;
  }

  return true;
}//highres_eval_solution


HighResSolution highres_mutate( const HighResSolution &X_base,
                                const std::function<double(void)> &rnd01,
                                double shrink_scale )
{
  HighResSolution X_new;
  const double mu = 0.2 * shrink_scale;
  const double mutate_threshold = PeakFitImprove::sm_ga_mutate_threshold;

  bool in_range;
  size_t num_tries = 0;

  do
  {
    num_tries += 1;
    if( num_tries > 1000 )
    {
      std::random_device rng;
      std::uniform_real_distribution<double> unif_dist( 0.0, 1.0 );
      highres_init_genes( X_new, [&](){ return unif_dist( rng ); } );
      return X_new;
    }

    in_range = true;
    X_new = X_base;

    if( rnd01() > mutate_threshold )
    {
      X_new.hpge_fwhm_bias_mult += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.hpge_fwhm_bias_mult >= 0.4) && (X_new.hpge_fwhm_bias_mult <= 1.3);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.nai_fwhm_bias_mult += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.nai_fwhm_bias_mult >= 0.7) && (X_new.nai_fwhm_bias_mult <= 1.5);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.hpge_distance_weight += mu * 5.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.hpge_distance_weight >= 0.1) && (X_new.hpge_distance_weight <= 5.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.nai_distance_weight += mu * 5.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.nai_distance_weight >= 0.1) && (X_new.nai_distance_weight <= 5.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.amp_clamp_denom += mu * 50.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.amp_clamp_denom >= 1.0) && (X_new.amp_clamp_denom <= 50.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.weight_clamp_max += mu * 30.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.weight_clamp_max >= 1.0) && (X_new.weight_clamp_max <= 30.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.min_peak_energy += mu * 70.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.min_peak_energy >= 30.0) && (X_new.min_peak_energy <= 100.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.max_peak_energy += mu * 3000.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.max_peak_energy >= 2000.0) && (X_new.max_peak_energy <= 5000.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.unknown_threshold += mu * 0.65 * (rnd01() - rnd01());
      in_range = in_range && (X_new.unknown_threshold >= 0.15) && (X_new.unknown_threshold <= 0.8);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.min_peaks_for_classify += static_cast<int>( mu * 4.0 * (rnd01() - rnd01()) );
      in_range = in_range && (X_new.min_peaks_for_classify >= 1) && (X_new.min_peaks_for_classify <= 5);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.min_peak_significance += mu * 5.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.min_peak_significance >= 2.25) && (X_new.min_peak_significance <= 5.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.sig_fraction_of_max += mu * 0.3 * (rnd01() - rnd01());
      in_range = in_range && (X_new.sig_fraction_of_max >= 0.005) && (X_new.sig_fraction_of_max <= 0.3);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.narrow_penalty_mult += mu * 10.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.narrow_penalty_mult >= 0.5) && (X_new.narrow_penalty_mult <= 10.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.narrow_sig_threshold += mu * 95.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.narrow_sig_threshold >= 5.0) && (X_new.narrow_sig_threshold <= 100.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.best_peak_czt_penalty += mu * 5.0 * (rnd01() - rnd01());
      in_range = in_range && (X_new.best_peak_czt_penalty >= 0.0) && (X_new.best_peak_czt_penalty <= 5.0);
    }
    if( rnd01() > mutate_threshold )
    {
      X_new.best_peak_conf_threshold += mu * 0.5 * (rnd01() - rnd01());
      in_range = in_range && (X_new.best_peak_conf_threshold >= 0.5) && (X_new.best_peak_conf_threshold <= 1.0);
    }
  } while( !in_range );

  return X_new;
}//highres_mutate


HighResSolution highres_crossover( const HighResSolution &X1,
                                   const HighResSolution &X2,
                                   const std::function<double(void)> &rnd01 )
{
  const double crossover_threshold = PeakFitImprove::sm_ga_crossover_threshold;
  HighResSolution X_new = X1;
  double r;

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.hpge_fwhm_bias_mult = r * X1.hpge_fwhm_bias_mult + (1.0 - r) * X2.hpge_fwhm_bias_mult;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.hpge_fwhm_bias_mult = X2.hpge_fwhm_bias_mult;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.nai_fwhm_bias_mult = r * X1.nai_fwhm_bias_mult + (1.0 - r) * X2.nai_fwhm_bias_mult;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.nai_fwhm_bias_mult = X2.nai_fwhm_bias_mult;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.hpge_distance_weight = r * X1.hpge_distance_weight + (1.0 - r) * X2.hpge_distance_weight;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.hpge_distance_weight = X2.hpge_distance_weight;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.nai_distance_weight = r * X1.nai_distance_weight + (1.0 - r) * X2.nai_distance_weight;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.nai_distance_weight = X2.nai_distance_weight;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.amp_clamp_denom = r * X1.amp_clamp_denom + (1.0 - r) * X2.amp_clamp_denom;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.amp_clamp_denom = X2.amp_clamp_denom;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.weight_clamp_max = r * X1.weight_clamp_max + (1.0 - r) * X2.weight_clamp_max;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.weight_clamp_max = X2.weight_clamp_max;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.min_peak_energy = r * X1.min_peak_energy + (1.0 - r) * X2.min_peak_energy;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.min_peak_energy = X2.min_peak_energy;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.max_peak_energy = r * X1.max_peak_energy + (1.0 - r) * X2.max_peak_energy;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.max_peak_energy = X2.max_peak_energy;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.unknown_threshold = r * X1.unknown_threshold + (1.0 - r) * X2.unknown_threshold;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.unknown_threshold = X2.unknown_threshold;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.min_peaks_for_classify = static_cast<int>(
      r * X1.min_peaks_for_classify + (1.0 - r) * X2.min_peaks_for_classify );
  }
  else if( rnd01() < 0.5 )
  {
    X_new.min_peaks_for_classify = X2.min_peaks_for_classify;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.min_peak_significance = r * X1.min_peak_significance + (1.0 - r) * X2.min_peak_significance;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.min_peak_significance = X2.min_peak_significance;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.sig_fraction_of_max = r * X1.sig_fraction_of_max + (1.0 - r) * X2.sig_fraction_of_max;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.sig_fraction_of_max = X2.sig_fraction_of_max;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.narrow_penalty_mult = r * X1.narrow_penalty_mult + (1.0 - r) * X2.narrow_penalty_mult;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.narrow_penalty_mult = X2.narrow_penalty_mult;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.narrow_sig_threshold = r * X1.narrow_sig_threshold + (1.0 - r) * X2.narrow_sig_threshold;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.narrow_sig_threshold = X2.narrow_sig_threshold;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.best_peak_czt_penalty = r * X1.best_peak_czt_penalty + (1.0 - r) * X2.best_peak_czt_penalty;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.best_peak_czt_penalty = X2.best_peak_czt_penalty;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.best_peak_conf_threshold = r * X1.best_peak_conf_threshold + (1.0 - r) * X2.best_peak_conf_threshold;
  }
  else if( rnd01() < 0.5 )
  {
    X_new.best_peak_conf_threshold = X2.best_peak_conf_threshold;
  }

  return X_new;
}//highres_crossover


namespace
{
  typedef EA::Genetic<ClassifyDetType_GA::HighResSolution, ClassifyDetType_GA::HighResCost> HighRes_GA_Type;

  double highres_calculate_SO_total_fitness( const HighRes_GA_Type::thisChromosomeType &X )
  {
    return X.middle_costs.objective1;
  }


  void highres_SO_report_generation(
    int generation_number,
    const EA::GenerationType<ClassifyDetType_GA::HighResSolution, ClassifyDetType_GA::HighResCost> &last_generation,
    const ClassifyDetType_GA::HighResSolution &best_genes )
  {
    if( !sm_highres_set_best_genes || (last_generation.best_total_cost < sm_highres_best_total_cost) )
    {
      sm_highres_set_best_genes = true;
      sm_highres_best_genes = best_genes;
      sm_highres_best_total_cost = last_generation.best_total_cost;
    }

    cout << "HighRes Generation [" << generation_number << "], "
         << "Best=" << last_generation.best_total_cost << ", "
         << "Average=" << last_generation.average_cost << ", "
         << "Best genes: {\n\t" << best_genes.to_string( "\n\t" ) << "\n}\n"
         << "Exe_time=" << last_generation.exe_time
         << endl << endl;

    sm_highres_output_file
      << generation_number << "\t"
      << last_generation.average_cost << "\t"
      << last_generation.best_total_cost << "\t"
      << "{" << best_genes.to_string( ", " ) << "}\n\n";
  }//highres_SO_report_generation


  // Stage 1 GA wrapper functions
  void highres_init_genes_wrapper( ClassifyDetType_GA::HighResSolution &p,
                                   const std::function<double(void)> &rnd01 )
  {
    ClassifyDetType_GA::highres_init_genes( p, rnd01 );
  }

  bool highres_eval_solution_wrapper( const ClassifyDetType_GA::HighResSolution &p,
                                      ClassifyDetType_GA::HighResCost &c )
  {
    return ClassifyDetType_GA::highres_eval_solution( p, c );
  }

  ClassifyDetType_GA::HighResSolution highres_mutate_wrapper(
    const ClassifyDetType_GA::HighResSolution &X_base,
    const std::function<double(void)> &rnd01,
    double shrink_scale )
  {
    return ClassifyDetType_GA::highres_mutate( X_base, rnd01, shrink_scale );
  }

  ClassifyDetType_GA::HighResSolution highres_crossover_wrapper(
    const ClassifyDetType_GA::HighResSolution &X1,
    const ClassifyDetType_GA::HighResSolution &X2,
    const std::function<double(void)> &rnd01 )
  {
    return ClassifyDetType_GA::highres_crossover( X1, X2, rnd01 );
  }
}//anonymous namespace


HighResClassifySettings do_highres_ga_eval(
  std::function<double( const HighResClassifySettings & )> eval_fcn )
{
  assert( !sm_highres_has_been_called );
  if( sm_highres_has_been_called )
  {
    cerr << "ClassifyDetType_GA::do_highres_ga_eval should only be called once!" << endl;
    exit( 1 );
  }

  sm_highres_has_been_called = true;

  assert( !!eval_fcn );
  if( !eval_fcn )
    throw runtime_error( "Invalid eval function passed to do_highres_ga_eval." );

  ns_highres_eval_fcn = eval_fcn;

  sm_highres_output_file.open( "det_classify_highres_results.txt" );
  sm_highres_output_file << "step\tcost_avg\tcost_best\tsolution_best\n";

  EA::Chronometer timer;
  timer.tic();

  HighRes_GA_Type ga_obj;
  ga_obj.problem_mode = EA::GA_MODE::SOGA;
  ga_obj.multi_threading = true;
  ga_obj.idle_delay_us = 1;
  ga_obj.dynamic_threading = true;
  ga_obj.verbose = false;
  ga_obj.population = static_cast<unsigned int>( PeakFitImprove::sm_ga_population );
  ga_obj.generation_max = static_cast<int>( PeakFitImprove::sm_ga_generation_max );
  ga_obj.calculate_SO_total_fitness = highres_calculate_SO_total_fitness;
  ga_obj.init_genes = highres_init_genes_wrapper;
  ga_obj.eval_solution = highres_eval_solution_wrapper;
  ga_obj.mutate = highres_mutate_wrapper;
  ga_obj.crossover = highres_crossover_wrapper;
  ga_obj.SO_report_generation = highres_SO_report_generation;
  ga_obj.crossover_fraction = PeakFitImprove::sm_ga_crossover_fraction;
  ga_obj.mutation_rate = PeakFitImprove::sm_ga_mutation_rate;
  ga_obj.best_stall_max = static_cast<int>( PeakFitImprove::sm_ga_best_stall_max );
  ga_obj.elite_count = static_cast<int>( PeakFitImprove::sm_ga_elite_count );
  ga_obj.N_threads = static_cast<int>( PeakFitImprove::sm_num_optimization_threads );

  const EA::StopReason stop_reason = ga_obj.solve();

  cout << "HighRes stage stop reason: " << ga_obj.stop_reason_to_string( stop_reason ) << endl;
  cout << "HighRes stage optimized in " << timer.toc() << " seconds." << endl;

  sm_highres_output_file.close();

  return highres_genes_to_settings( sm_highres_best_genes );
}//do_highres_ga_eval


}//namespace ClassifyDetType_GA
