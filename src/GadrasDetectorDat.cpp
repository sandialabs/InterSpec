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

#include <cmath>
#include <limits>
#include <array>
#include <cstdio>
#include <string>
#include <cctype>
#include <memory>
#include <vector>
#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/GadrasDetectorDat.h"

using namespace std;

namespace
{
  const double ns_pi = 3.14159265358979323846;

  /** Canonical GADRAS parameter labels (1-based; index 0 unused), transcribed from
   `GADRAS_DRF_parameter_documentation.md` §2.  Used to write a text file when a parsed file did
   not supply a verbatim label (e.g. when converting an XML-sourced detector to text). */
  const char * const sm_canonical_labels[GadrasDetectorDat::sm_num_params + 1] = {
    /* 0 unused */                    "",
    /* 1  */ "A1: e-cal ord. 1",      /* 2  */ "A2: e-cal ord. 2",
    /* 3  */ "A3: e-cal ord. 3",      /* 4  */ "A4: e-cal ord. 4",
    /* 5  */ "A5: e-cal, low-E",      /* 6  */ "res. @ E=0 (keV)",
    /* 7  */ "% FWHM @ 661 keV",      /* 8  */ "resolution Power",
    /* 9  */ "solid angle (%)",       /* 10 */ "det. length (cm)",
    /* 11 */ "det. width (cm)",       /* 12 */ "height/width",
    /* 13 */ "shape factor",          /* 14 */ "attenuator Z",
    /* 15 */ "attenuator g/cm2",      /* 16 */ "porosity (%)",
    /* 17 */ "distance (cm)",         /* 18 */ "eff. scalar",
    /* 19 */ "outer atten. Z",        /* 20 */ "outer atten AD",
    /* 21 */ "attenuate scat",        /* 22 */ "clutter",
    /* 23 */ "scat prob @ 0",         /* 24 */ "scat prob @ 45",
    /* 25 */ "scat prob @ 90",        /* 26 */ "scat prob @ 135",
    /* 27 */ "scat prob @ 180",       /* 28 */ "scat @ E-> Edge",
    /* 29 */ "high-E skew",           /* 30 */ "scat @ E -> 0",
    /* 31 */ "scat = f(E)",           /* 32 */ "outer porosity",
    /* 33 */ "ext annihilation",      /* 34 */ "shaping time, us",
    /* 35 */ "LLD(keV)",              /* 36 */ "low-E noise",
    /* 37 */ "height (cm)",           /* 38 */ "Frisch grid (%)",
    /* 39 */ "shield xray 1 scalar",  /* 40 */ "det setback (cm)",
    /* 41 */ "shield angle (%)",      /* 42 */ "shield thick, cm",
    /* 43 */ "dead layer Z",          /* 44 */ "dead layer g/cm2",
    /* 45 */ "shield LLD (keV)",      /* 46 */ "Bremsstrahlung",
    /* 47 */ "low-E skew",            /* 48 */ "hole mu-tau (cm)",
    /* 49 */ "% neutron shield",      /* 50 */ "neutron reflect",
    /* 51 */ "dead layer (mm)",       /* 52 */ "bad pole zero / INBIN",
    /* 53 */ "efficiency holder / COINC",
    /* 54 */ "# coincidence array / PILEUP",
    /* 55 */ "enhance Compton escape",
    /* 56 */ "alt collimator dia. / REBIN",
    /* 57 */ "skew low-E power",
    /* 58 */ "neutron environment / local xray 1 Z",
    /* 59 */ "scatter environment / gamma detector",
    /* 60 */ "unused / anticoincidence detector",
    /* 61 */ "air pressure / imager type",
    /* 62 */ "template error / min energy keV",
    /* 63 */ "chi-square / max energy keV",
    /* 64 */ "local xray 2 scalar / # channels",
    /* 65 */ "side shield +X (%) / local xray 2 Z",
    /* 66 */ "AN side shield",        /* 67 */ "AD side shield",
    /* 68 */ "back shield (%)",       /* 69 */ "AN back shield",
    /* 70 */ "AD back shield",        /* 71 */ "skew low-e extent",
    /* 72 */ "LLD sharpness",         /* 73 */ "skew high-E power",
    /* 74 */ "peak scatter angle",    /* 75 */ "peak scatter width",
    /* 76 */ "peak scatter amount",   /* 77 */ "skew high-e extent",
    /* 78 */ "default dead time per pulse (us)",
    /* 79 */ "aerial start (cm)",     /* 80 */ "aerial stop (cm)",
    /* 81 */ "side shield -X (%)",    /* 82 */ "side shield +Y (%)",
    /* 83 */ "side shield -Y (%)",    /* 84 */ "source delta height (cm)"
  };


  /** Default GADRAS material table (`DetectorData.gadras`), transcribed from
   `GADRAS_DRF_parameter_documentation.md` §3.1.2 (density g/cm3, resolution % FWHM @ 661 keV).
   `specialCased` marks the names GADRAS keys material-specific physics branches off of (§3.1.3).

     #  Name         Formula                              rho     res%
     1  BGO          Bi4Ge3O12                            7.125   16
     2  CsI          CsI                                  4.51     8
     3  HPGe         Ge                                   5.33     0.25
     4  HgI          HgI2                                 6.4      1
     5  NaI          NaI                                  3.667    7
     6  PVT          C1H1.104                             1.03    20
     7  Si           Si                                   2.33     0.3
     8  CZT          Cd0.96Zn0.04Te                       5.78     2
     9  CdTe         CdTe                                 5.85     3
    10  LaCl3        LaCl3                                3.84     4
    11  LaBr3        LaBr3                                5.06     2.5
    12  CdWO4        CdWO4                                7.90    14
    13  Xe           Xe                                   0.1      4
    14  CaF2         CaF2                                 3.18    10
    15  CLYC(nat)    Cs2LiYCl6                            3.31     4
    16  CLYC(95%)    Cs2LiYCl6                            3.31     4
    17  SrI2         SrI2                                 4.59     3
    18  YAP          YAlO3                                5.37     6
    19  PFCBBr       C8H3OBr4                             3.6     20
    20  BaFBr:Eu     BaFBr                                4.96    20
    21  CeBr3        CeBr3                                5.06     4.5
    22  LiI(nat)     LiI                                  4.076   10
    23  LiI(95%)     LiI                                  4.076   10
    24  Gd2O2S       Gd2O2S                               7.3     20
    25  CLLBC(95%)   Cs2LiLaBr4.8Cl1.2                    4.06     3.5
    26  LSI(95%)     LiSr2I2                              4.4      3.5
    27  BiPlastic40  C616H657N11O44Bi10                   1.4     15
    28  CLLB(95%)    Cs2LiLaBr6                           4.2      3
    29  TlBr         TlBr                                 7.56     3
    30  GAGG:Ce      Gd1.49Y1.49Ga2.2Al2.8O12Ce0.02      6.7      5.5
    31  GYGAG        Gd1.49Y1.49Ga2.5Al2.5O12Ce0.01      5.8     10.0
    32  Stilbene     C4H12                                1.16     8.0
    33  Deut.Stilb.  C14[H2]12                            1.24     8.0
    34  EJ301        H1C0.826                             0.874   11.1
    35  EJ301D       H0.0246[H2]1C0.845                   0.948   11.1
    36  BiPiv17      H239C205O6Bi1                        1.4     28.0
   */
  const std::array<GadrasDetectorDat::MaterialInfo, 36> sm_material_table = {{
    {  1, "BGO",                 "Bi4Ge3O12",                        7.125f, 16.0f,  true  },
    {  2, "CsI",                 "CsI",                              4.51f,   8.0f,  true  },
    {  3, "HPGe",                "Ge",                               5.33f,   0.25f, true  },
    {  4, "HgI",                 "HgI2",                             6.4f,    1.0f,  false },
    {  5, "NaI",                 "NaI",                              3.667f,  7.0f,  true  },
    {  6, "PVT",                 "C1H1.104",                         1.03f,  20.0f,  true  },
    {  7, "Si",                  "Si",                               2.33f,   0.3f,  false },
    {  8, "CZT",                 "Cd0.96Zn0.04Te",                   5.78f,   2.0f,  true  },
    {  9, "CdTe",                "CdTe",                             5.85f,   3.0f,  true  },
    { 10, "LaCl3",               "LaCl3",                            3.84f,   4.0f,  true  },
    { 11, "LaBr3",               "LaBr3",                            5.06f,   2.5f,  true  },
    { 12, "CdWO4",               "CdWO4",                            7.90f,  14.0f,  true  },
    { 13, "Xe",                  "Xe",                               0.1f,    4.0f,  false },
    { 14, "CaF2",                "CaF2",                             3.18f,  10.0f,  true  },
    { 15, "CLYC(nat)",           "Cs2LiYCl6",                        3.31f,   4.0f,  true  },
    { 16, "CLYC(95%)",           "Cs2LiYCl6",                        3.31f,   4.0f,  true  },
    { 17, "SrI2",                "SrI2",                             4.59f,   3.0f,  true  },
    { 18, "YAP",                 "YAlO3",                            5.37f,   6.0f,  false },
    { 19, "PFCBBr",              "C8H3OBr4",                         3.6f,   20.0f,  false },
    { 20, "BaFBr:Eu",            "BaFBr",                            4.96f,  20.0f,  false },
    { 21, "CeBr3",               "CeBr3",                            5.06f,   4.5f,  false },
    { 22, "LiI(nat)",            "LiI",                              4.076f, 10.0f,  true  },
    { 23, "LiI(95%)",            "LiI",                              4.076f, 10.0f,  true  },
    { 24, "Gd2O2S",              "Gd2O2S",                           7.3f,   20.0f,  false },
    { 25, "CLLBC(95%)",          "Cs2LiLaBr4.8Cl1.2",                4.06f,   3.5f,  true  },
    { 26, "LSI(95%)",            "LiSr2I2",                          4.4f,    3.5f,  false },
    { 27, "BiPlastic40",         "C616H657N11O44Bi10",               1.4f,   15.0f,  false },
    { 28, "CLLB(95%)",           "Cs2LiLaBr6",                       4.2f,    3.0f,  true  },
    { 29, "TlBr",                "TlBr",                             7.56f,   3.0f,  false },
    { 30, "GAGG:Ce",             "Gd1.49Y1.49Ga2.2Al2.8O12Ce0.02",  6.7f,    5.5f,  false },
    { 31, "GYGAG",               "Gd1.49Y1.49Ga2.5Al2.5O12Ce0.01",  5.8f,   10.0f,  false },
    { 32, "Stilbene",            "C4H12",                            1.16f,   8.0f,  false },
    { 33, "Deuterated Stilbene", "C14[H2]12",                        1.24f,   8.0f,  false },
    { 34, "EJ301",               "H1C0.826",                         0.874f, 11.1f,  false },
    { 35, "EJ301D",              "H0.0246[H2]1C0.845",               0.948f, 11.1f,  false },
    { 36, "BiPiv17",             "H239C205O6Bi1",                    1.4f,   28.0f,  false }
  }};


  /** Descriptor for a numeric XML leaf and the GADRAS parameter it maps to.  Shared by the XML
   reader and writer so the two stay in sync.  @p path is slash-delimited relative to the root
   `<gamma_detector>`.  @p use_int_col selects which column of the parameter the leaf feeds (the
   value column, or the overloaded integer column - doc §2).  @p bare is true for leaves stored
   as node text directly (no nested `<value>`/`<varying>`). */
  struct XmlNumMap
  {
    int         param;
    bool        use_int_col;
    bool        bare;
    const char *path;
  };

  const XmlNumMap sm_xml_num_map[] = {
    // top-level wrapped
    { 17, false, false, "distance" },
    { 37, false, false, "height" },
    { 61, false, false, "air_pressure" },
    { 62, true,  false, "weight_range_lower" },
    { 63, true,  false, "weight_range_upper" },
    { 54, false, false, "number_coincidence_array" },
    { 55, false, false, "enhance_compton_escape" },
    { 33, false, false, "external_annihilation" },
    { 46, false, false, "external_bremsstrahlung" },
    { 38, false, false, "frisch_grid_percent" },
    { 36, false, false, "low_e_noise" },
    // top-level bare
    { 62, false, true,  "template_error" },
    { 64, true,  true,  "default_channel_count" },
    // dimensions
    { 40, false, false, "dimensions/setback" },
    { 10, false, false, "dimensions/length" },
    { 11, false, false, "dimensions/width" },
    { 12, false, false, "dimensions/height_to_width_ratio" },
    { 13, false, false, "dimensions/shape_factor" },
    { 51, false, false, "dimensions/dead_layer" },
    { 18, false, false, "dimensions/scalar" },
    // peak_shape
    {  6, false, false, "peak_shape/fwhm_offset" },
    {  7, false, false, "peak_shape/fwhm_at_661keV" },
    {  8, false, false, "peak_shape/fwhm_power" },
    { 47, false, false, "peak_shape/low_energy_skew" },
    { 29, false, false, "peak_shape/high_energy_skew" },
    { 57, false, false, "peak_shape/low_energy_skew_power" },
    { 73, false, false, "peak_shape/high_energy_skew_power" },
    { 71, false, false, "peak_shape/low_energy_skew_extent" },
    { 77, false, false, "peak_shape/high_energy_skew_extent" },
    { 48, false, false, "peak_shape/percent_holes_trapped" },
    // lower_level_discrimination
    { 35, false, false, "lower_level_discrimination/cutoff_energy" },
    { 72, false, false, "lower_level_discrimination/sharpness" },
    // photon_scatter
    { 22, false, false, "photon_scatter/clutter" },
    { 23, false, false, "photon_scatter/degree0" },
    { 24, false, false, "photon_scatter/degree45" },
    { 25, false, false, "photon_scatter/degree90" },
    { 26, false, false, "photon_scatter/degree135" },
    { 27, false, false, "photon_scatter/degree180" },
    { 28, false, false, "photon_scatter/compton_edge" },
    { 30, false, false, "photon_scatter/low_energy" },
    { 31, false, false, "photon_scatter/increase_with_energy" },
    { 21, false, false, "photon_scatter/attenuation" },
    { 74, false, false, "photon_scatter/pref_angle" },
    { 75, false, false, "photon_scatter/pref_angle_width" },
    { 76, false, false, "photon_scatter/pref_angle_magnitude" },
    // neutron_scatter
    { 49, false, false, "neutron_scatter/thermals_stopped_percent" },
    { 50, false, false, "neutron_scatter/reflection_scalar" },
    // energy_calibration
    {  1, false, false, "energy_calibration/offset" },
    {  2, false, false, "energy_calibration/range" },
    {  3, false, false, "energy_calibration/quadratic" },
    {  4, false, false, "energy_calibration/cubic" },
    {  5, false, false, "energy_calibration/low_energy" },
    // inner_attenuator
    { 16, false, false, "inner_attenuator/porosity" },
    { 14, false, false, "inner_attenuator/material/atomic_number" },
    { 15, false, false, "inner_attenuator/material/areal_density" },
    // outer_attenuator
    { 32, false, false, "outer_attenuator/porosity" },
    { 19, false, false, "outer_attenuator/material/atomic_number" },
    { 20, false, false, "outer_attenuator/material/areal_density" },
    // timing
    { 34, false, false, "timing/shaping_time" },
    { 78, false, false, "timing/dead_time_per_pulse" },
    // side_shield
    { 65, false, false, "side_shield/coverage_percent_plus_x" },
    { 81, false, false, "side_shield/coverage_percent_minus_x" },
    { 82, false, false, "side_shield/coverage_percent_plus_y" },
    { 83, false, false, "side_shield/coverage_percent_minus_y" },
    { 66, false, false, "side_shield/material/atomic_number" },
    { 67, false, false, "side_shield/material/areal_density" },
    { 56, false, false, "side_shield/alt_collimator_diameter" },
    // back_shield
    { 68, false, false, "back_shield/coverage_percent" },
    { 69, false, false, "back_shield/material/atomic_number" },
    { 70, false, false, "back_shield/material/areal_density" },
    // fluorescent x-rays
    { 39, false, false, "fluorescent_xray1/magnitude" },
    { 58, true,  true,  "fluorescent_xray1/atomic_number" },
    { 64, false, false, "fluorescent_xray2/magnitude" },
    { 65, true,  true,  "fluorescent_xray2/atomic_number" },
    // anti_coincidence_shield
    { 41, false, false, "anti_coincidence_shield/solid_angle_percent" },
    { 42, false, false, "anti_coincidence_shield/thickness" },
    { 45, false, false, "anti_coincidence_shield/lower_level_discriminator" },
    { 43, false, false, "anti_coincidence_shield/dead_layer_material/atomic_number" },
    { 44, false, false, "anti_coincidence_shield/dead_layer_material/areal_density" }
  };


  /** Format a value as GADRAS writes it: fixed with 6 decimals. */
  string fmt_f( const double v )
  {
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.6f", v );
    return string( buffer );
  }

  /** Parse a float from XML text, tolerating scientific notation ("8.060000e-01"). */
  float parse_float( const string &str )
  {
    float val = 0.0f;
    if( !(stringstream(str) >> val) )
      throw runtime_error( "Non-numeric value '" + str + "'" );
    return val;
  }

  /** Walk a slash-delimited @p path from @p parent, returning the leaf node or nullptr. */
  const rapidxml::xml_node<char> *node_by_path( const rapidxml::xml_node<char> *parent,
                                                const string &path )
  {
    vector<string> segs;
    SpecUtils::split( segs, path, "/" );
    const rapidxml::xml_node<char> *node = parent;
    for( const string &seg : segs )
    {
      if( !node )
        return nullptr;
      node = node->first_node( seg.c_str(), seg.size() );
    }
    return node;
  }
}//anonymous namespace


float GadrasDetectorDat::value( const int idx ) const
{
  if( (idx < 1) || (idx > sm_num_params) )
    return 0.0f;
  return m_params[idx].value;
}


int GadrasDetectorDat::intCol( const int idx ) const
{
  if( (idx < 1) || (idx > sm_num_params) )
    return 0;
  return m_params[idx].int_col;
}


void GadrasDetectorDat::setValue( const int idx, const float v )
{
  if( (idx < 1) || (idx > sm_num_params) )
    return;
  m_params[idx].value = v;
  m_highest_param = std::max( m_highest_param, idx );
}


const char *GadrasDetectorDat::canonicalLabel( const int idx )
{
  if( (idx < 1) || (idx > sm_num_params) )
    return "";
  return sm_canonical_labels[idx];
}


float GadrasDetectorDat::equivalentCircularDiameterCm() const
{
  const float surface_area = width() * width() * heightToWidth();
  if( surface_area <= 0.0f )
    return 0.0f;
  return 2.0f * std::sqrt( surface_area / static_cast<float>(ns_pi) );
}


GadrasDetectorDat::Attenuator GadrasDetectorDat::innerAttenuator() const
{
  Attenuator a;
  a.atomicNumber = value(14);
  a.arealDensity = value(15);
  a.porosityPercent = value(16);
  return a;
}


GadrasDetectorDat::Attenuator GadrasDetectorDat::outerAttenuator() const
{
  Attenuator a;
  a.atomicNumber = value(19);
  a.arealDensity = value(20);
  a.porosityPercent = value(32);
  return a;
}


GadrasDetectorDat::SideShield GadrasDetectorDat::sideShield() const
{
  SideShield s;
  s.atomicNumber = value(66);
  s.arealDensity = value(67);
  s.covPlusX = value(65);
  s.covMinusX = value(81);
  s.covPlusY = value(82);
  s.covMinusY = value(83);
  return s;
}


GadrasDetectorDat::BackShield GadrasDetectorDat::backShield() const
{
  BackShield b;
  b.atomicNumber = value(69);
  b.arealDensity = value(70);
  b.coverage = value(68);
  return b;
}


int GadrasDetectorDat::materialIndex() const
{
  return intCol(59);
}


std::string GadrasDetectorDat::materialName() const
{
  if( !m_material_name.empty() )
    return m_material_name;

  const MaterialInfo * const info = materialByIndex( materialIndex() );
  return info ? string(info->name) : string();
}


const std::array<GadrasDetectorDat::MaterialInfo, 36> &GadrasDetectorDat::materialTable()
{
  return sm_material_table;
}


const GadrasDetectorDat::MaterialInfo *GadrasDetectorDat::materialByIndex( const int idx )
{
  if( (idx < 1) || (idx > static_cast<int>(sm_material_table.size())) )
    return nullptr;
  return &sm_material_table[idx - 1];
}


const GadrasDetectorDat::MaterialInfo *GadrasDetectorDat::materialByName( const std::string &name )
{
  if( name.empty() )
    return nullptr;

  string wanted = name;
  SpecUtils::trim( wanted );
  SpecUtils::to_lower_ascii( wanted );

  for( const MaterialInfo &info : sm_material_table )
  {
    string candidate = info.name;
    SpecUtils::to_lower_ascii( candidate );
    if( candidate == wanted )
      return &info;
  }

  return nullptr;
}



GadrasDetectorDat::InferredShape
GadrasDetectorDat::inferShape( const std::string &materialOverride ) const
{
  // Implements the decision recipe of GADRAS_DRF_parameter_documentation.md §5.1.
  InferredShape result;

  const float h2w = heightToWidth();
  const float len = length();
  const float wid = width();

  // Signal 1: HeightToWidth != 1  =>  rectangular parallelepiped (literal dimensions).
  if( std::fabs(h2w - 1.0f) > 0.01f )
  {
    result.shape = Shape::Rectangular;
    result.dimA = len;              // axial depth
    result.dimB = wid;              // width
    result.dimC = h2w * wid;        // height
    return result;
  }

  // Square cross-section: disambiguate by material family (Signal 2), shapeFactor as fallback.
  string mat = materialOverride.empty() ? materialName() : materialOverride;
  SpecUtils::to_lower_ascii( mat );

  auto contains = [&mat]( const char * const s ) -> bool {
    return mat.find(s) != string::npos;
  };

  if( contains("hpge") || contains("germanium") )
  {
    result.shape = Shape::CoaxialCylinder;   // bore not recoverable from Detector.dat
    result.dimA = 0.5f * wid;                // radius
    result.dimB = len;                       // height
  }else if( contains("czt") || contains("cdte") || contains("cadmium")
            || contains("tlbr") || (mat == "si") )
  {
    result.shape = Shape::Box;               // solid-state block
    result.dimA = len;
    result.dimB = wid;
    result.dimC = wid;
  }else if( contains("nai") || contains("csi") || contains("labr") || contains("lacl")
            || contains("cebr") || contains("sri") || contains("clyc") || contains("bgo")
            || contains("cllb") || contains("cdwo") || contains("caf") )
  {
    result.shape = Shape::Cylinder;
    result.dimA = 0.5f * wid;
    result.dimB = len;
  }else
  {
    // Weak fallback: low shape factor => rounded (cylinder), high => box.
    if( shapeFactor() <= 30.0f )
    {
      result.shape = Shape::Cylinder;
      result.dimA = 0.5f * wid;
      result.dimB = len;
    }else
    {
      result.shape = Shape::Box;
      result.dimA = len;
      result.dimB = wid;
      result.dimC = wid;
    }
  }

  return result;
}


GadrasDetectorDat::ChordDistribution
GadrasDetectorDat::computeChordDistribution( const double length, const double width,
                                             const double heightToWidthRatio,
                                             const double shapeFactor, const double distance )
{
  // Reproduces GADRAS PathLengthComputer.f90:Compute (doc §6.2).
  ChordDistribution out;

  const double height = heightToWidthRatio * width;
  const double L = 0.5 * length;                        // axial half-depth
  double H, W;                                          // H = larger, W = smaller (taper axis)
  if( width > height ) { H = 0.5*width;  W = 0.5*height; }
  else                { W = 0.5*width;  H = 0.5*height; }

  if( (W <= 0.0) || (L <= 0.0) || (shapeFactor <= 0.0) )
    return out;

  const double D = std::min( distance, 10.0*(H + W) );  // clamp far source
  const double P = 1.0 / shapeFactor;                   // profile exponent
  const double thetaMax = std::atan( H / D );           // aperture half-angle (set by H)
  const double phiMax   = 0.5 * ns_pi;                  // one transverse quadrant

  const int    KMAX = 256;                              // polar samples
  const int    JMAX = std::max( 4, KMAX/4 );            // azimuth samples (= 64)
  const double DT   = 0.01 * std::min( L, std::min(W, H) );  // ray-march step

  auto halfThickness = [&]( const double y ) { return L * std::pow( 1.0 - y/W, P ); };

  double SAlast = 0.0;
  for( int k = 1; k <= KMAX; ++k )
  {
    const double theta = thetaMax * (double(k) / KMAX);
    const double SA    = 2.0 * ns_pi * std::sin(theta);
    const double solidAngle = SA - SAlast;
    SAlast = SA;

    const double ct = std::cos(theta), st = std::sin(theta);

    for( int j = 1; j <= JMAX; ++j )
    {
      const double phi = phiMax * ((j - 0.5) / JMAX);
      const double dX = ct;
      const double dY = st * std::cos(phi);
      const double dZ = st * std::sin(phi);

      const double T = 0.98 * D / dX;
      double X = T*dX - D;
      double Y = T*dY;
      double Z = T*dZ;

      if( (Y > W) || (Z > H) )
        continue;

      double XL1 = L - halfThickness(Y);
      while( (Z < H) && (Y < W) && (X < XL1) )
      {
        X += dX*DT; Y += dY*DT; Z += dZ*DT;
        XL1 = L - halfThickness(Y);
      }
      if( !((Z < H) && (Y < W)) )
        continue;

      const double X1 = X, Y1 = Y, Z1 = Z;

      double XL2 = L + halfThickness(Y);
      while( (Z < H) && (Y < W) && (X < XL2) )
      {
        X += dX*DT; Y += dY*DT; Z += dZ*DT;
        XL2 = L + halfThickness(Y);
      }

      const double dx = X - X1, dy = Y - Y1, dz = Z - Z1;
      const double pathLength = std::sqrt( dx*dx + dy*dy + dz*dz );

      const double FL = pathLength / (2.0 * L);
      const int i = int( std::max( 1.0, std::floor(32.0*FL + 0.5) ) );
      if( i < 128 )
        out.lengths[i] += solidAngle;
    }//for( j )
  }//for( k )

  double sum = 0.0;
  for( double v : out.lengths )
    sum += v;
  if( sum > 0.0 )
  {
    const double sc = 1.0 / sum;
    for( int i = 0; i < 128; ++i )
    {
      if( out.lengths[i] > 0.0 )
      {
        out.lengths[i] *= sc;
        out.count = i;
      }
    }
  }

  return out;
}


GadrasDetectorDat GadrasDetectorDat::fromStream( std::istream &input )
{
  GadrasDetectorDat result;

  const istream::pos_type orig_pos = input.tellg();

  string first_line;
  if( !SpecUtils::safe_get_line( input, first_line, 64 * 1024 ) )
    throw runtime_error( "Could not read first line of Detector.dat" );

  const bool is_xml = (first_line.find("xml") != string::npos);

  input.clear();
  input.seekg( orig_pos );

  if( is_xml )
    result.parseXml( input );
  else
    result.parseText( input );

  return result;
}


GadrasDetectorDat GadrasDetectorDat::fromFile( const std::string &path )
{
#ifdef _WIN32
  const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16( path );
  ifstream input( wpath.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream input( path.c_str(), ios_base::binary | ios_base::in );
#endif
  if( !input.is_open() )
    throw runtime_error( "Could not open Detector.dat file '" + path + "'." );

  return fromStream( input );
}


bool GadrasDetectorDat::isCandidateDetectorDat( std::istream &input )
{
  const std::istream::pos_type start = input.tellg();

  bool answer = false;
  try
  {
    const GadrasDetectorDat dat = GadrasDetectorDat::fromStream( input );

    if( dat.m_format == SourceFormat::Xml )
    {
      // fromStream only reaches the XML branch on a <gamma_detector> root.
      answer = true;
    }else
    {
      // A text file has to look like a real parameter table, not just contain a
      //  line that happens to start with an integer.
      int populated = 0;
      for( int i = 1; i <= sm_num_params; ++i )
      {
        if( (dat.m_params[i].value != 0.0f) || (dat.m_params[i].int_col != 0) )
          populated += 1;
      }

      answer = (dat.m_highest_param >= 40) && (populated >= 10)
               && (dat.resFWHM661() > 0.0f) && (dat.width() > 0.0f);
    }
  }catch( std::exception & )
  {
    answer = false;
  }

  input.clear();
  input.seekg( start );

  return answer;
}//isCandidateDetectorDat(...)


void GadrasDetectorDat::parseText( std::istream &input )
{
  m_format = SourceFormat::LegacyText;

  string line;
  while( SpecUtils::safe_get_line( input, line, 64 * 1024 ) )
  {
    // The literal text "NaN" is replaced with 0 before parsing (doc §1.1).
    SpecUtils::ireplace_all( line, "NaN", "0" );
    SpecUtils::trim( line );

    if( line.empty() )
      continue;

    // The newer text files carry a "Version <x.y.z>" line between params 64 and 65.
    if( SpecUtils::istarts_with( line, "Version" ) )
    {
      string ver = line.substr( 7 );
      SpecUtils::trim( ver );
      m_response_version = ver;
      continue;
    }

    if( !std::isdigit( static_cast<unsigned char>(line[0]) ) )
      continue;

    vector<string> parts;
    SpecUtils::split( parts, line, " \t" );
    if( parts.size() < 2 )
      continue;

    try
    {
      const int parnum = std::stoi( parts[0] );
      if( (parnum < 1) || (parnum > sm_num_params) )
        continue;

      m_params[parnum].value = static_cast<float>( std::stod( parts[1] ) );
      if( parts.size() > 2 )
      {
        try{ m_params[parnum].int_col = std::stoi( parts[2] ); }catch( ... ){}
      }

      // Capture the verbatim label (everything after the first three tokens) for lossless
      //  text round-trip.
      if( parts.size() > 3 )
      {
        string label;
        for( size_t i = 3; i < parts.size(); ++i )
          label += (i == 3 ? "" : " ") + parts[i];
        m_labels[parnum] = label;
      }

      m_highest_param = std::max( m_highest_param, parnum );
    }catch( ... )
    {
      // Skip malformed lines rather than abort the whole parse.
      continue;
    }
  }//while( read a line )

  if( m_highest_param < 1 )
    throw runtime_error( "Detector.dat contained no valid parameter lines." );

  // Back-compat (doc §2): files with <=80 params default side-shield -X/+Y/-Y coverage (81-83)
  //  to the +X value (65), and source-delta-height (84) to 0.
  if( m_highest_param < 81 )
  {
    m_params[81].value = m_params[82].value = m_params[83].value = m_params[65].value;
  }
}//parseText


void GadrasDetectorDat::parseXml( std::istream &input )
{
  m_format = SourceFormat::Xml;

  try
  {
    rapidxml::file<char> input_file( input );
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace>( input_file.data() );

    const rapidxml::xml_node<char> * const root = doc.first_node( "gamma_detector" );
    if( !root )
      throw runtime_error( "Missing <gamma_detector> node." );

    // Numeric leaves from the shared mapping table.
    for( const XmlNumMap &m : sm_xml_num_map )
    {
      const rapidxml::xml_node<char> * const node = node_by_path( root, m.path );
      if( !node )
        continue;

      string val_str;
      if( m.bare )
      {
        val_str = SpecUtils::xml_value_str( node );
      }else
      {
        const rapidxml::xml_node<char> * const value_node
                                          = SpecUtils::xml_first_node( node, "value" );
        if( !value_node )
          continue;
        val_str = SpecUtils::xml_value_str( value_node );
      }

      if( val_str.empty() )
        continue;

      const float val = parse_float( val_str );
      if( m.use_int_col )
        m_params[m.param].int_col = static_cast<int>( std::lround(val) );
      else
        m_params[m.param].value = val;
    }//for( each numeric leaf )

    // String / meta leaves, preserved for round-trip.
    auto capture_str = [&]( const string &path ) {
      const rapidxml::xml_node<char> * const node = node_by_path( root, path );
      if( node )
      {
        const string txt = SpecUtils::xml_value_str( node );
        if( !txt.empty() )
          m_xml_extras[path] = txt;
      }
    };

    capture_str( "photon_scatter/environment" );
    capture_str( "neutron_scatter/environment" );
    capture_str( "timing/compute_pileup" );
    capture_str( "custom_input_channel_widths" );
    capture_str( "rebin_all_spectra" );
    capture_str( "anti_coincidence_shield/material" );

    if( const rapidxml::xml_node<char> * const n = node_by_path( root, "file_version" ) )
      m_file_version = SpecUtils::xml_value_str( n );
    if( const rapidxml::xml_node<char> * const n = node_by_path( root, "response_version" ) )
      m_response_version = SpecUtils::xml_value_str( n );

    // The crystal material is a name in XML; resolve to the material-table index so that
    //  materialIndex() is uniform across formats.
    if( const rapidxml::xml_node<char> * const n = root->first_node( "material", 8 ) )
    {
      m_material_name = SpecUtils::xml_value_str( n );
      const MaterialInfo * const info = materialByName( m_material_name );
      m_params[59].int_col = info ? info->index : 0;
    }

    m_highest_param = sm_num_params;
  }catch( std::exception &e )
  {
    throw runtime_error( "Failed to read XML Detector.dat: " + string(e.what()) );
  }
}//parseXml


void GadrasDetectorDat::toText( std::ostream &out ) const
{
  const int highest = (m_highest_param > 0) ? m_highest_param : sm_num_params;

  for( int i = 1; i <= highest; ++i )
  {
    const string label = m_labels[i].empty() ? string(canonicalLabel(i)) : m_labels[i];

    char buffer[128] = { '\0' };
    // GADRAS format: (I3,F13.5,I7,3X,A)
    snprintf( buffer, sizeof(buffer), "%3d%13.5f%7d   %s",
              i, m_params[i].value, m_params[i].int_col, label.c_str() );
    out << buffer << "\n";

    // The "Version" line is wedged in between params 64 and 65 in newer files.
    if( (i == 64) && !m_response_version.empty() )
      out << "Version " << m_response_version << "\n";
  }
}//toText


void GadrasDetectorDat::toXml( std::ostream &out ) const
{
  // Emit the real GADRAS XML schema (verified against sample Detector.dat files).  Written
  //  directly rather than via rapidxml so the fixed-decimal formatting and node order match
  //  what GADRAS produces.  Node -> parameter mapping mirrors sm_xml_num_map above.
  const string t1 = "\t", t2 = "\t\t", t3 = "\t\t\t", t4 = "\t\t\t\t";

  auto wrapped = [&]( const string &indent, const char * const name, const int param,
                      const bool use_int_col ) {
    const double v = use_int_col ? static_cast<double>( intCol(param) )
                                 : static_cast<double>( value(param) );
    out << indent << "<" << name << ">\n"
        << indent << t1 << "<value>" << fmt_f(v) << "</value>\n"
        << indent << t1 << "<varying>false</varying>\n"
        << indent << "</" << name << ">\n";
  };

  auto bare_num = [&]( const string &indent, const char * const name, const int param,
                       const bool use_int_col ) {
    if( use_int_col )
      out << indent << "<" << name << ">" << intCol(param) << "</" << name << ">\n";
    else
      out << indent << "<" << name << ">" << fmt_f( value(param) ) << "</" << name << ">\n";
  };

  auto extra = [&]( const string &path, const char * const dflt ) -> string {
    const auto pos = m_xml_extras.find( path );
    return (pos != m_xml_extras.end()) ? pos->second : string(dflt);
  };

  const string mat_name = materialName();

  out << "<?xml version=\"1.0\"?>\n";
  out << "<gamma_detector>\n";

  out << t1 << "<file_version>" << (m_file_version.empty() ? "1.2.0" : m_file_version)
      << "</file_version>\n";
  out << t1 << "<response_version>" << m_response_version << "</response_version>\n";

  wrapped( t1, "distance", 17, false );
  wrapped( t1, "height", 37, false );

  if( !mat_name.empty() )
    out << t1 << "<material>" << mat_name << "</material>\n";

  out << t1 << "<dimensions>\n";
  wrapped( t2, "setback", 40, false );
  wrapped( t2, "length", 10, false );
  wrapped( t2, "width", 11, false );
  wrapped( t2, "height_to_width_ratio", 12, false );
  wrapped( t2, "shape_factor", 13, false );
  wrapped( t2, "dead_layer", 51, false );
  wrapped( t2, "scalar", 18, false );
  out << t1 << "</dimensions>\n";

  out << t1 << "<peak_shape>\n";
  wrapped( t2, "fwhm_offset", 6, false );
  wrapped( t2, "fwhm_at_661keV", 7, false );
  wrapped( t2, "fwhm_power", 8, false );
  wrapped( t2, "low_energy_skew", 47, false );
  wrapped( t2, "high_energy_skew", 29, false );
  wrapped( t2, "low_energy_skew_power", 57, false );
  wrapped( t2, "high_energy_skew_power", 73, false );
  wrapped( t2, "low_energy_skew_extent", 71, false );
  wrapped( t2, "high_energy_skew_extent", 77, false );
  wrapped( t2, "percent_holes_trapped", 48, false );
  out << t1 << "</peak_shape>\n";

  out << t1 << "<lower_level_discrimination>\n";
  wrapped( t2, "cutoff_energy", 35, false );
  wrapped( t2, "sharpness", 72, false );
  out << t1 << "</lower_level_discrimination>\n";

  out << t1 << "<photon_scatter>\n";
  out << t2 << "<environment>" << extra("photon_scatter/environment", "Unknown")
      << "</environment>\n";
  wrapped( t2, "clutter", 22, false );
  wrapped( t2, "degree0", 23, false );
  wrapped( t2, "degree45", 24, false );
  wrapped( t2, "degree90", 25, false );
  wrapped( t2, "degree135", 26, false );
  wrapped( t2, "degree180", 27, false );
  wrapped( t2, "compton_edge", 28, false );
  wrapped( t2, "low_energy", 30, false );
  wrapped( t2, "increase_with_energy", 31, false );
  wrapped( t2, "attenuation", 21, false );
  wrapped( t2, "pref_angle", 74, false );
  wrapped( t2, "pref_angle_width", 75, false );
  wrapped( t2, "pref_angle_magnitude", 76, false );
  out << t1 << "</photon_scatter>\n";

  out << t1 << "<neutron_scatter>\n";
  out << t2 << "<environment>" << extra("neutron_scatter/environment", "Unknown")
      << "</environment>\n";
  wrapped( t2, "thermals_stopped_percent", 49, false );
  wrapped( t2, "reflection_scalar", 50, false );
  out << t1 << "</neutron_scatter>\n";

  out << t1 << "<energy_calibration>\n";
  wrapped( t2, "offset", 1, false );
  wrapped( t2, "range", 2, false );
  wrapped( t2, "quadratic", 3, false );
  wrapped( t2, "cubic", 4, false );
  wrapped( t2, "low_energy", 5, false );
  out << t1 << "</energy_calibration>\n";

  out << t1 << "<inner_attenuator>\n";
  wrapped( t2, "porosity", 16, false );
  out << t2 << "<material>\n";
  wrapped( t3, "atomic_number", 14, false );
  wrapped( t3, "areal_density", 15, false );
  out << t2 << "</material>\n";
  out << t1 << "</inner_attenuator>\n";

  out << t1 << "<outer_attenuator>\n";
  wrapped( t2, "porosity", 32, false );
  out << t2 << "<material>\n";
  wrapped( t3, "atomic_number", 19, false );
  wrapped( t3, "areal_density", 20, false );
  out << t2 << "</material>\n";
  out << t1 << "</outer_attenuator>\n";

  out << t1 << "<timing>\n";
  out << t2 << "<compute_pileup>" << extra("timing/compute_pileup", "false")
      << "</compute_pileup>\n";
  wrapped( t2, "shaping_time", 34, false );
  wrapped( t2, "dead_time_per_pulse", 78, false );
  out << t1 << "</timing>\n";

  out << t1 << "<side_shield>\n";
  wrapped( t2, "coverage_percent_plus_x", 65, false );
  wrapped( t2, "coverage_percent_minus_x", 81, false );
  wrapped( t2, "coverage_percent_plus_y", 82, false );
  wrapped( t2, "coverage_percent_minus_y", 83, false );
  out << t2 << "<material>\n";
  wrapped( t3, "atomic_number", 66, false );
  wrapped( t3, "areal_density", 67, false );
  out << t2 << "</material>\n";
  wrapped( t2, "alt_collimator_diameter", 56, false );
  out << t1 << "</side_shield>\n";

  out << t1 << "<back_shield>\n";
  wrapped( t2, "coverage_percent", 68, false );
  out << t2 << "<material>\n";
  wrapped( t3, "atomic_number", 69, false );
  wrapped( t3, "areal_density", 70, false );
  out << t2 << "</material>\n";
  out << t1 << "</back_shield>\n";

  out << t1 << "<fluorescent_xray1>\n";
  bare_num( t2, "atomic_number", 58, true );
  wrapped( t2, "magnitude", 39, false );
  out << t1 << "</fluorescent_xray1>\n";

  out << t1 << "<fluorescent_xray2>\n";
  bare_num( t2, "atomic_number", 65, true );
  wrapped( t2, "magnitude", 64, false );
  out << t1 << "</fluorescent_xray2>\n";

  wrapped( t1, "air_pressure", 61, false );
  wrapped( t1, "weight_range_lower", 62, true );
  wrapped( t1, "weight_range_upper", 63, true );
  bare_num( t1, "default_channel_count", 64, true );
  bare_num( t1, "template_error", 62, false );
  out << t1 << "<custom_input_channel_widths>" << extra("custom_input_channel_widths", "false")
      << "</custom_input_channel_widths>\n";
  out << t1 << "<rebin_all_spectra>" << extra("rebin_all_spectra", "false")
      << "</rebin_all_spectra>\n";
  wrapped( t1, "number_coincidence_array", 54, false );
  wrapped( t1, "enhance_compton_escape", 55, false );
  wrapped( t1, "external_annihilation", 33, false );
  wrapped( t1, "external_bremsstrahlung", 46, false );
  wrapped( t1, "frisch_grid_percent", 38, false );

  out << t1 << "<anti_coincidence_shield>\n";
  out << t2 << "<material>" << extra("anti_coincidence_shield/material", "None") << "</material>\n";
  wrapped( t2, "solid_angle_percent", 41, false );
  wrapped( t2, "thickness", 42, false );
  wrapped( t2, "lower_level_discriminator", 45, false );
  out << t2 << "<dead_layer_material>\n";
  wrapped( t3, "atomic_number", 43, false );
  wrapped( t3, "areal_density", 44, false );
  out << t2 << "</dead_layer_material>\n";
  out << t1 << "</anti_coincidence_shield>\n";

  wrapped( t1, "low_e_noise", 36, false );

  out << "</gamma_detector>\n";
}//toXml
