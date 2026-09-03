#pragma once
/** Recorded extended-source Monte-Carlo truth for the near-field volumetric-efficiency validation.

 Absolute full-energy-peak efficiency per emitted gamma, averaged over the source volume, for the
 scenarios in VolumetricNearFieldScenarios.h against the ANGLE GEM35-70 in test_data/det_eff.

 *** GENERATED - do not hand-edit. ***
 Regenerate with:
   ./test_VolumetricNearField --run_test=VolumetricNearFieldTruth -- --datadir=... --testfiledir=...
 and paste the printed table below.  Note the FEP window it was generated at (printed in the header
 line): a model that credits in-window Compton must use the same one.
 */

#include <string>
#include <vector>

/** One point of the Monte-Carlo-derived anchor curve: absolute FEP efficiency for a bare POINT
 source on axis at kMcAnchorDistanceCm. */
struct AnchorRow
{
  double energy_keV = 0.0;
  double eff = 0.0;
  double frac_sigma = 0.0;
};

/** Transfer-vs-MC for a bare POINT source on axis, using the file's MEASURED anchor: the shared
 baseline every volumetric row would otherwise inherit.  Recorded so the check is free. */
struct PointBaselineRow
{
  double distance_cm = 0.0;
  double energy_keV = 0.0;
  double mc_fep_eff = 0.0;
};

struct TruthRow
{
  std::string scenario;
  double energy_keV = 0.0;
  double fep_eff = 0.0;      ///< absolute FEP efficiency per emitted gamma
  double fep_uncert = 0.0;   ///< 1-sigma MC statistical uncertainty
};

// Regenerated 2026-09-02 by test_VolumetricLadder --run_test=Rung4_ScenarioMatrixTruth at a 0.25%
//  target FEP precision (the previous bank was 0.5%), through the ladder's on-disk Monte-Carlo cache
//  in scratch/20260902_volumetric_ladder.  FEP window: CeeLo's kDefaultFepWindowKeV = 0.75 keV, a
//  HALF-width (the MC scores |E_dep - E_src| < window), so a model that credits in-window Compton
//  must use the same one.  No row stopped on the event or wall cap.
static const std::vector<AnchorRow> sm_mc_anchor = {
  { 40, 0.0002166922291, 0.002468 },
  { 50, 0.0007000196867, 0.002393 },
  { 60, 0.001205237029, 0.002303 },
  { 75, 0.001739125546, 0.002187 },
  { 88, 0.002007897455, 0.002142 },
  { 105, 0.002178171418, 0.002093 },
  { 122, 0.00223576164, 0.002079 },
  { 150, 0.002199547186, 0.002109 },
  { 200, 0.00197907852, 0.002137 },
  { 280, 0.001604614818, 0.002217 },
  { 344, 0.001370186157, 0.002264 },
  { 500, 0.00102519318, 0.002314 },
  { 661.7, 0.0008390627995, 0.00235 },
  { 900, 0.0006702607004, 0.002383 },
  { 1332.5, 0.0005048321619, 0.002405 },
  { 1800, 0.0004036696217, 0.00242 },
  { 2614, 0.0002924228087, 0.00244 },
};

static const std::vector<PointBaselineRow> sm_point_baseline = {
  { 10, 60, 0.006348818681 },
  { 10, 88, 0.01060806327 },
  { 10, 122, 0.01163132708 },
  { 10, 344, 0.006748015781 },
  { 10, 661.7, 0.003932524577 },
  { 10, 1332.5, 0.002346853963 },
  { 25, 60, 0.001213150487 },
  { 25, 88, 0.002015041987 },
  { 25, 122, 0.002241559521 },
  { 25, 344, 0.00138201105 },
  { 25, 661.7, 0.0008307100209 },
  { 25, 1332.5, 0.0005041898495 },
  { 50, 60, 0.0003186948267 },
  { 50, 88, 0.0005319992951 },
  { 50, 122, 0.000593328362 },
  { 50, 344, 0.0003775475808 },
  { 50, 661.7, 0.0002307336055 },
  { 50, 1332.5, 0.0001409865908 },
};

static const std::vector<TruthRow> sm_truth = {
  { "small-near-light", 60, 0.05832350746, 0.0001432 },
  { "small-near-light", 88, 0.1044469388, 0.0002523 },
  { "small-near-light", 122, 0.1161603053, 0.0002799 },
  { "small-near-light", 344, 0.0595221374, 0.0001462 },
  { "small-near-light", 661.7, 0.0337739151, 8.334e-05 },
  { "small-near-light", 1332.5, 0.01952172903, 4.84e-05 },
  { "large-near-light", 60, 0.01812808219, 4.508e-05 },
  { "large-near-light", 88, 0.03440501089, 8.508e-05 },
  { "large-near-light", 122, 0.03986379747, 9.844e-05 },
  { "large-near-light", 344, 0.0236441791, 5.87e-05 },
  { "large-near-light", 661.7, 0.01446608771, 3.597e-05 },
  { "large-near-light", 1332.5, 0.008823442557, 2.199e-05 },
  { "small-far-light", 60, 0.0002852179473, 7.089e-07 },
  { "small-far-light", 88, 0.0004798720633, 1.183e-06 },
  { "small-far-light", 122, 0.0005380795956, 1.247e-06 },
  { "small-far-light", 344, 0.0003499377921, 8.364e-07 },
  { "small-far-light", 661.7, 0.0002177296272, 5.275e-07 },
  { "small-far-light", 1332.5, 0.0001352776223, 3.314e-07 },
  { "large-far-light", 60, 0.0002104192649, 5.261e-07 },
  { "large-far-light", 88, 0.0003617282186, 9.102e-07 },
  { "large-far-light", 122, 0.0004087901935, 9.467e-07 },
  { "large-far-light", 344, 0.0002846832096, 6.792e-07 },
  { "large-far-light", 661.7, 0.0001834396183, 4.455e-07 },
  { "large-far-light", 1332.5, 0.0001173904691, 2.88e-07 },
  { "wide-angle-far-light", 60, 0.001926823862, 4.752e-06 },
  { "wide-angle-far-light", 88, 0.003536504548, 8.631e-06 },
  { "wide-angle-far-light", 122, 0.004124745375, 9.933e-06 },
  { "wide-angle-far-light", 344, 0.002550400673, 6.3e-06 },
  { "wide-angle-far-light", 661.7, 0.001560152469, 3.84e-06 },
  { "wide-angle-far-light", 1332.5, 0.0009473556593, 2.343e-06 },
  { "shielded-near-light", 60, 0.0002436049699, 6.09e-07 },
  { "shielded-near-light", 88, 0.0092079145, 2.296e-05 },
  { "shielded-near-light", 122, 0.02640951586, 6.552e-05 },
  { "shielded-near-light", 344, 0.02759982548, 6.844e-05 },
  { "shielded-near-light", 661.7, 0.01819712647, 4.518e-05 },
  { "shielded-near-light", 1332.5, 0.01176795326, 2.927e-05 },
  { "box-large-near-light", 60, 0.01595145436, 3.968e-05 },
  { "box-large-near-light", 88, 0.03086523438, 7.643e-05 },
  { "box-large-near-light", 122, 0.03607025172, 8.92e-05 },
  { "box-large-near-light", 344, 0.02198030513, 5.46e-05 },
  { "box-large-near-light", 661.7, 0.01339211196, 3.332e-05 },
  { "box-large-near-light", 1332.5, 0.008283014277, 2.064e-05 },
  { "box-shielded-near-light", 60, 0.0002301288875, 5.753e-07 },
  { "box-shielded-near-light", 88, 0.008726710454, 2.176e-05 },
  { "box-shielded-near-light", 122, 0.02514, 6.237e-05 },
  { "box-shielded-near-light", 344, 0.02643405676, 6.555e-05 },
  { "box-shielded-near-light", 661.7, 0.01748083389, 4.343e-05 },
  { "box-shielded-near-light", 1332.5, 0.01139150747, 2.834e-05 },
  { "box-slab-near-light", 60, 0.03933175, 9.719e-05 },
  { "box-slab-near-light", 88, 0.0739952381, 0.0001806 },
  { "box-slab-near-light", 122, 0.08419891304, 0.0002047 },
  { "box-slab-near-light", 344, 0.04591695906, 0.0001132 },
  { "box-slab-near-light", 661.7, 0.02677909582, 6.629e-05 },
  { "box-slab-near-light", 1332.5, 0.01574277971, 3.911e-05 },
  { "small-near-dense", 60, 0.00829958355, 2.07e-05 },
  { "small-near-dense", 88, 0.03444148472, 8.521e-05 },
  { "small-near-dense", 122, 0.05872030075, 0.0001441 },
  { "small-near-dense", 344, 0.04420985915, 0.0001091 },
  { "small-near-dense", 661.7, 0.02705819629, 6.691e-05 },
  { "small-near-dense", 1332.5, 0.0166302547, 4.128e-05 },
  { "large-near-dense", 60, 0.0009057817543, 2.264e-06 },
  { "large-near-dense", 88, 0.004298869752, 1.073e-05 },
  { "large-near-dense", 122, 0.008760164835, 2.184e-05 },
  { "large-near-dense", 344, 0.01026374759, 2.558e-05 },
  { "large-near-dense", 661.7, 0.007371897842, 1.838e-05 },
  { "large-near-dense", 1332.5, 0.005273964911, 1.316e-05 },
  { "small-far-dense", 60, 3.782969763e-05, 9.456e-08 },
  { "small-far-dense", 88, 0.0001536347874, 3.84e-07 },
  { "small-far-dense", 122, 0.0002729087088, 6.819e-07 },
  { "small-far-dense", 344, 0.000261982745, 6.589e-07 },
  { "small-far-dense", 661.7, 0.0001745058739, 4.235e-07 },
  { "small-far-dense", 1332.5, 0.0001146357118, 2.81e-07 },
  { "large-far-dense", 60, 9.358319116e-06, 2.34e-08 },
  { "large-far-dense", 88, 3.970300351e-05, 9.925e-08 },
  { "large-far-dense", 122, 8.085086955e-05, 2.021e-07 },
  { "large-far-dense", 344, 0.0001142096628, 2.852e-07 },
  { "large-far-dense", 661.7, 8.848907387e-05, 2.207e-07 },
  { "large-far-dense", 1332.5, 6.737428208e-05, 1.677e-07 },
  { "wide-angle-far-dense", 60, 0.0002363582579, 5.909e-07 },
  { "wide-angle-far-dense", 88, 0.001071480309, 2.675e-06 },
  { "wide-angle-far-dense", 122, 0.00199417362, 4.964e-06 },
  { "wide-angle-far-dense", 344, 0.001865156432, 4.611e-06 },
  { "wide-angle-far-dense", 661.7, 0.00123176773, 3.037e-06 },
  { "wide-angle-far-dense", 1332.5, 0.0007989212389, 1.978e-06 },
  { "shielded-near-dense", 60, 1.842988553e-05, 4.607e-08 },
  { "shielded-near-dense", 88, 0.001741509639, 4.352e-06 },
  { "shielded-near-dense", 122, 0.008644685466, 2.156e-05 },
  { "shielded-near-dense", 344, 0.01618584521, 4.027e-05 },
  { "shielded-near-dense", 661.7, 0.01209552346, 3.01e-05 },
  { "shielded-near-dense", 1332.5, 0.008780976886, 2.187e-05 },
  { "box-large-near-dense", 60, 0.0007703337989, 1.926e-06 },
  { "box-large-near-dense", 88, 0.003735109864, 9.327e-06 },
  { "box-large-near-dense", 122, 0.007741893204, 1.931e-05 },
  { "box-large-near-dense", 344, 0.009287296037, 2.316e-05 },
  { "box-large-near-dense", 661.7, 0.006748194208, 1.683e-05 },
  { "box-large-near-dense", 1332.5, 0.004892197465, 1.221e-05 },
  { "box-shielded-near-dense", 60, 1.732921037e-05, 4.332e-08 },
  { "box-shielded-near-dense", 88, 0.001636789845, 4.09e-06 },
  { "box-shielded-near-dense", 122, 0.008156982097, 2.034e-05 },
  { "box-shielded-near-dense", 344, 0.0154786965, 3.85e-05 },
  { "box-shielded-near-dense", 661.7, 0.01164496133, 2.898e-05 },
  { "box-shielded-near-dense", 1332.5, 0.008429220168, 2.099e-05 },
  { "box-slab-near-dense", 60, 0.004952233251, 1.236e-05 },
  { "box-slab-near-dense", 88, 0.02205368567, 5.477e-05 },
  { "box-slab-near-dense", 122, 0.03985367089, 9.842e-05 },
  { "box-slab-near-dense", 344, 0.03323242105, 8.224e-05 },
  { "box-slab-near-dense", 661.7, 0.02083458517, 5.167e-05 },
  { "box-slab-near-dense", 1332.5, 0.01324086649, 3.291e-05 },
};
