// Baseline known-good GADRAS Detector.dat / Efficiency.csv content embedded
// as raw string literals.  Used by the GADRAS fuzz harnesses so each fuzzer
// can mutate one of the two files while the other is held constant.  Content
// taken from data/GenericGadrasDetectors/HPGe 40%/.

#pragma once

namespace fuzz_gadras_baseline {

inline constexpr const char *k_detector_dat = R"GADRAS(  1      0.00000      0   A1: e-cal ord. 1
  2   3000.00000      0   A2: e-cal ord. 2
  3      0.00000      0   A3: e-cal ord. 3
  4      0.00000      0   A4: e-cal ord. 4
  5      0.00000      0   A5: e-cal, low-E
  6      1.55000      0   res. @ E=0 (keV)
  7      0.25000      0   % FWHM @ 661 keV
  8      0.35000      0   resolution Power
  9      0.02400      0   solid angle (%)
 10      0.00000      0   det. length (cm)
 11      5.85000      0   det. width (cm)
 12      1.00000      0   height/width
 13     20.00000      0   shape factor
 14     13.00000      0   attenuator Z
 15      2.00000      0   attenuator g/cm2
 16      0.00000      0   porosity (%)
 17    100.00000      0   distance (cm)
 18      1.00000      0   eff. scalar
 19      4.00000      0   outer atten. Z
 20      0.00000      0   outer atten AD
 21      0.00000      1   attenuate scat
 22      0.00000      1   clutter
 23      0.00000      0   scat prob @ 0
 24      0.00000      0   scat prob @ 45
 25      0.00000      0   scat prob @ 90
 26      0.00000      0   scat prob @ 135
 27      0.00000      0   scat prob @ 180
 28      0.00000      0   scat @ E-> Edge
 29      0.00000      0   high-E skew
 30      0.00000      0   scat @ E -> 0
 31      0.00000      0   scat = f(E)
 32      0.00000      0   outer porosity
 33      0.00000      0   ext annihilation
 34      0.00000      0   shaping time, us
 35      0.00000      0   LLD(keV)
 36      0.00000      0   low-E noise
 37    100.00000      0   height (cm)
 38      0.00000      0   Frisch grid (%)
 39      0.00000      0   sheild xray 1 scalar
 40      0.50000      0   det setback (cm)
 41      0.00000      0   shield angle (%)
 42      0.00000      0   shield thick, cm
 43      0.00000      0   dead layer Z
 44      0.00000      0   dead layer g/cm2
 45      0.00000      0   shield LLD (keV)
 46      0.00000      0   Bremsstrahlung
 47      0.00000      0   low-E skew
 48      0.00000      0   hole mu-tau (cm)
 49      0.00000      0   % neutron shield
 50      0.00000      0   neutron reflect
 51      0.00000      0   dead layer (mm)
 52      0.00000      0   bad pole zero       / INBIN
 53      0.00000      0   efficiency holder   / COINC
 54      0.00000      1   # coincidence array / PILEUP
 55      0.00000      0   enhance Compton escape (coincidence arrays)
 56      0.00000      0   alt collimator dia. / REBIN
 57      0.00000      0   skew low-E power    / adjustable if value=1
 58      0.00000      0   neutron environment / local xray 1 Z
 59      0.00000      0   scatter environment / gamma detector
 60      0.00000      0   unused              / anticoincidence detector
 61      1.00000      0   air pressure, atm   / 0) AC shield; 1) CC or 2) CA imagers
 62      0.00000      0   template error      / min energy, keV
 63      0.00000      0   chi-square          / max energy, keV
 64      0.00000      0   local xray 2 scalar / # of channels
Version 18.7.6
 65      0.00000      0   side shield (%)     / local xray 2 Z
 66      0.00000      0   AN side shield
 67      0.00000      0   AD side shield
 68      0.00000      0   back shield (%)
 69      0.00000      0   AN back shield
 70      0.00000      0   AD back shield
 71      0.00000      0   skew extent +/-
 72      0.00000      0   LLD sharpness
 73      0.00000      0   skew high-E power   / adjustable if value=1
 74      0.00000      0   peak scatter angle  / adjustable if value=1
 75      0.00000      0   peak scatter width  / adjustable if value=1
 76      0.00000      0   peak scatter amount / adjustable if value=1
 77      0.00000      0   reserved
 78      0.00000      0   reserved
 79      0.00000      0   reserved
 80      0.00000      0   reserved
)GADRAS";

inline constexpr const char *k_efficiency_csv = R"GADRAS(Energy,Peak,PCOM,PCOM1,SE,DE,PTOT,Ei*Rij
(keV),(%),(%),(%),(%),(%),(%)
10,0.0000,0.000,0.000,0.000,0.000,0.000,0.000
50,4.9709,0.001,0.001,0.000,0.000,7.661,3.266
100,49.899,1.161,1.324,0.000,0.000,108.204,80.715
500,25.746,41.521,28.947,0.000,0.000,111.589,286.673
1000,15.589,44.900,31.528,0.000,0.000,91.908,454.125
2000,9.0851,38.169,27.823,0.455,0.365,73.159,719.086
3000,6.3107,32.875,24.108,1.033,0.827,64.916,952.074
)GADRAS";

}  // namespace fuzz_gadras_baseline
