#ifndef VolumetricIntegrationFixtures_h
#define VolumetricIntegrationFixtures_h
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

/** Shared volumetric-source geometry fixtures for the self-attenuation /
 trace-source integration tests.

 These fixtures define a matrix of {geometry, shell-stack, energy, distance}
 scenarios, and a builder that constructs a `DistributedSrcCalc` the same way
 `ShieldingSourceChi2Fcn::energy_chi_contributions` does.  They are used to:
  - record Cuhre baseline integral values (regression net),
  - validate any replacement integration backend against those baselines,
  - exercise off-axis and effective-AN/AD extensions on identical geometry.
 */

#include <array>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>

#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"

namespace VolumetricFixture
{

/** One layer of shielding; either a material from the MaterialDB, or a
 zero-thickness generic layer defined by atomic number and areal density.
 */
struct ShellSpec
{
  /** Material name in MaterialDB (e.g. "Fe (iron)"); empty means generic. */
  std::string material;

  /** Only used when `material` is empty. */
  double atomic_number = 0.0;
  /** Only used when `material` is empty; in PhysicalUnits (e.g. g/cm2). */
  double areal_density = 0.0;

  /** Thicknesses this layer adds to the cumulative outer dimensions; meaning
   depends on geometry: {radius} for Spherical, {radius, half-length} for
   cylinders, {half-width, half-height, half-depth} for Rectangular.
   Ignored (must be zero) for generic layers, which add no physical extent.
   */
  std::array<double,3> thicknesses = { 0.0, 0.0, 0.0 };
};//struct ShellSpec


/** A complete volumetric-source scenario. */
struct VolumetricSrcSpec
{
  std::string name;
  GammaInteractionCalc::GeometryType geometry;
  std::vector<ShellSpec> shells;

  /** Index into `shells` of the layer the source is distributed in;
   must be a material (non-generic) layer.
   */
  size_t source_shell_index = 0;

  double energy = 661.657 * PhysicalUnits::keV;

  /** Distance from assembly center to detector face. */
  double distance = 100.0 * PhysicalUnits::cm;

  double detector_radius = 2.54 * PhysicalUnits::cm;
  double detector_setback = 0.0;

  bool attenuate_for_air = false;

  /** Lateral offsets of the source assembly from the detector axis. */
  double offset_x = 0.0;
  double offset_y = 0.0;

  /** In-situ exponentially-distributed (trace) source options. */
  bool in_situ_exponential = false;
  double relaxation_length = -1.0;
};//struct VolumetricSrcSpec


/** Builds a `DistributedSrcCalc` from a scenario, mirroring how
 `ShieldingSourceChi2Fcn::energy_chi_contributions` populates one.
 `m_srcVolumetricActivity` is set to 1.0 so `integral` is directly comparable.

 Throws std::runtime_error on invalid specs or unknown materials.
 */
inline GammaInteractionCalc::DistributedSrcCalc
   make_distributed_src_calc( const VolumetricSrcSpec &spec, const MaterialDB &matdb )
{
  GammaInteractionCalc::DistributedSrcCalc calc;

  calc.m_geometry = spec.geometry;
  calc.m_materialIndex = spec.source_shell_index;
  calc.m_detectorRadius = spec.detector_radius;
  calc.m_detectorSetback = spec.detector_setback;
  calc.m_observationDist = spec.distance;
  calc.m_srcOffsetX = spec.offset_x;
  calc.m_srcOffsetY = spec.offset_y;
  calc.m_attenuateForAir = spec.attenuate_for_air;
  calc.m_airTransLenCoef = spec.attenuate_for_air
              ? GammaInteractionCalc::transmission_length_coefficient_air( static_cast<float>(spec.energy) )
              : 0.0;
  calc.m_isInSituExponential = spec.in_situ_exponential;
  calc.m_inSituRelaxationLength = spec.relaxation_length;
  calc.m_srcVolumetricActivity = 1.0;
  calc.m_energy = spec.energy;
  calc.m_nuclide = nullptr;
  calc.integral = 0.0;

  if( spec.source_shell_index >= spec.shells.size() )
    throw std::runtime_error( "make_distributed_src_calc: source_shell_index out of range" );

  std::array<double,3> outer_dims = { 0.0, 0.0, 0.0 };

  for( size_t shell_index = 0; shell_index < spec.shells.size(); ++shell_index )
  {
    const ShellSpec &shell = spec.shells[shell_index];

    if( shell.material.empty() )
    {
      // Generic layers have no physical extent - same convention as production
      //  code, which skips a generic layer at the very center (zero volume).
      if( shell_index == 0 )
        throw std::runtime_error( "make_distributed_src_calc: generic layer cant be innermost" );
      if( shell_index == spec.source_shell_index )
        throw std::runtime_error( "make_distributed_src_calc: source layer cant be generic" );

      const double trans_coef = GammaInteractionCalc::transmition_coefficient_generic(
                                          static_cast<float>(shell.atomic_number),
                                          static_cast<float>(shell.areal_density),
                                          static_cast<float>(spec.energy) );

      calc.m_dimensionsTransLenAndType.push_back( {outer_dims, trans_coef,
                    GammaInteractionCalc::DistributedSrcCalc::ShellType::Generic} );
      continue;
    }//if( a generic layer )

    const std::shared_ptr<const Material> material = matdb.material( shell.material );
    if( !material )
      throw std::runtime_error( "make_distributed_src_calc: no material '" + shell.material + "'" );

    switch( spec.geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        outer_dims[0] += shell.thicknesses[0];
        break;

      case GammaInteractionCalc::GeometryType::CylinderEndOn:
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        outer_dims[0] += shell.thicknesses[0];
        outer_dims[1] += shell.thicknesses[1];
        break;

      case GammaInteractionCalc::GeometryType::Rectangular:
        outer_dims[0] += shell.thicknesses[0];
        outer_dims[1] += shell.thicknesses[1];
        outer_dims[2] += shell.thicknesses[2];
        break;

      case GammaInteractionCalc::GeometryType::NumGeometryType:
        throw std::runtime_error( "make_distributed_src_calc: invalid geometry" );
    }//switch( spec.geometry )

    const double trans_len_coef = GammaInteractionCalc::transmition_length_coefficient(
                                          material.get(), static_cast<float>(spec.energy) );

    calc.m_dimensionsTransLenAndType.push_back( {outer_dims, trans_len_coef,
                  GammaInteractionCalc::DistributedSrcCalc::ShellType::Material} );
  }//for( loop over shells )

  return calc;
}//make_distributed_src_calc(...)


/** Builds the templated `DistributedSrcCalcT<double>` equivalent of
 #make_distributed_src_calc, for validating the templated ray-tracing and
 Gauss-Legendre integration backend against the Cuhre/double path.
 */
inline GammaInteractionCalc::DistributedSrcCalcT<double>
   make_distributed_src_calc_t( const VolumetricSrcSpec &spec, const MaterialDB &matdb )
{
  const GammaInteractionCalc::DistributedSrcCalc legacy = make_distributed_src_calc( spec, matdb );

  GammaInteractionCalc::DistributedSrcCalcT<double> calc;
  calc.m_geometry = legacy.m_geometry;
  calc.m_materialIndex = legacy.m_materialIndex;
  calc.m_detector = GammaInteractionCalc::detector_geom_from_config<double>(
                          spec.geometry, spec.distance, spec.detector_radius, spec.detector_setback,
                          spec.offset_x, spec.offset_y );
  calc.m_attenuateForAir = legacy.m_attenuateForAir;
  calc.m_airTransLenCoef = legacy.m_airTransLenCoef;
  calc.m_isInSituExponential = legacy.m_isInSituExponential;
  calc.m_inSituRelaxationLength = legacy.m_inSituRelaxationLength;
  calc.m_srcVolumetricActivity = legacy.m_srcVolumetricActivity;
  calc.m_energy = legacy.m_energy;
  calc.m_nuclide = legacy.m_nuclide;
  calc.integral = 0.0;

  size_t legacy_index = 0;
  for( const std::tuple<std::array<double,3>,double,GammaInteractionCalc::DistributedSrcCalc::ShellType> &shell
                                                          : legacy.m_dimensionsTransLenAndType )
  {
    GammaInteractionCalc::DistributedSrcCalcT<double>::ShellInfo info;
    info.dims = std::get<0>( shell );
    info.trans_len_coef = std::get<1>( shell );
    info.type = std::get<2>( shell );

    // Metadata for the effective-AN/AD/H accumulation
    const ShellSpec &shell_spec = spec.shells[legacy_index];
    if( shell_spec.material.empty() )
    {
      info.areal_density = shell_spec.areal_density;
      info.effective_an = shell_spec.atomic_number;
    }else
    {
      const std::shared_ptr<const Material> mat = matdb.material( shell_spec.material );
      if( mat )
      {
        info.density = mat->density;
        info.effective_an = GammaInteractionCalc::material_mass_weighted_atomic_number( *mat );
        info.hydrogen_mass_frac = GammaInteractionCalc::material_hydrogen_mass_fraction( *mat );
      }
    }//if( generic ) / else

    calc.m_shells.push_back( info );
    ++legacy_index;
  }//for( loop over legacy shells )

  return calc;
}//make_distributed_src_calc_t(...)


/** The standard fixture matrix used to record Cuhre baselines, and to compare
 any new integration backend / capability against.
 */
inline std::vector<VolumetricSrcSpec> standard_volumetric_fixtures()
{
  using GammaInteractionCalc::GeometryType;
  const double cm = PhysicalUnits::cm;
  const double keV = PhysicalUnits::keV;
  const double g_per_cm2 = PhysicalUnits::g / PhysicalUnits::cm2;

  std::vector<VolumetricSrcSpec> fixtures;

  {// Spherical: bare self-attenuating U sphere, strongly attenuated low-energy line
    VolumetricSrcSpec spec;
    spec.name = "sph-U-self-atten-185keV";
    spec.geometry = GeometryType::Spherical;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 0.0, 0.0}} };
    spec.energy = 185.72*keV;
    fixtures.push_back( spec );
  }

  {// Spherical: same sphere, higher-energy line
    VolumetricSrcSpec spec;
    spec.name = "sph-U-self-atten-1001keV";
    spec.geometry = GeometryType::Spherical;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 0.0, 0.0}} };
    spec.energy = 1001.03*keV;
    fixtures.push_back( spec );
  }

  {// Spherical: U source sphere inside an iron shell
    VolumetricSrcSpec spec;
    spec.name = "sph-U-plus-Fe-661keV";
    spec.geometry = GeometryType::Spherical;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 0.0, 0.0}},
                    {"Fe (iron)", 0.0, 0.0, {1.0*cm, 0.0, 0.0}} };
    fixtures.push_back( spec );
  }

  {// Spherical: U source sphere inside a generic (AN/AD of iron) layer
    VolumetricSrcSpec spec;
    spec.name = "sph-U-plus-generic-661keV";
    spec.geometry = GeometryType::Spherical;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 0.0, 0.0}},
                    {"", 26.0, 7.874*g_per_cm2, {0.0, 0.0, 0.0}} };
    fixtures.push_back( spec );
  }

  {// Spherical: source in an outer shell, with a non-source inner core
   //  (exercises the ray-through-inner-shells path)
    VolumetricSrcSpec spec;
    spec.name = "sph-innervoid-Fe-then-U-661keV";
    spec.geometry = GeometryType::Spherical;
    spec.shells = { {"Fe (iron)", 0.0, 0.0, {1.5*cm, 0.0, 0.0}},
                    {"U (uranium)", 0.0, 0.0, {1.0*cm, 0.0, 0.0}} };
    spec.source_shell_index = 1;
    fixtures.push_back( spec );
  }

  {// Spherical: in-situ exponentially-distributed trace source in soil
    VolumetricSrcSpec spec;
    spec.name = "sph-trace-exp-soil-661keV";
    spec.geometry = GeometryType::Spherical;
    spec.shells = { {"Soil", 0.0, 0.0, {25.0*cm, 0.0, 0.0}} };
    spec.in_situ_exponential = true;
    spec.relaxation_length = 3.0*cm;
    fixtures.push_back( spec );
  }

  {// Spherical: bare U sphere with air attenuation, at 2 m
    VolumetricSrcSpec spec;
    spec.name = "sph-U-air-atten-661keV";
    spec.geometry = GeometryType::Spherical;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 0.0, 0.0}} };
    spec.distance = 200.0*cm;
    spec.attenuate_for_air = true;
    fixtures.push_back( spec );
  }

  {// CylinderEndOn: single iron cylinder (fast single-shell code path, 2D)
    VolumetricSrcSpec spec;
    spec.name = "cyl-end-Fe-661keV";
    spec.geometry = GeometryType::CylinderEndOn;
    spec.shells = { {"Fe (iron)", 0.0, 0.0, {5.0*cm, 5.0*cm, 0.0}} };
    fixtures.push_back( spec );
  }

  {// CylinderEndOn: U cylinder inside iron (multi-shell path, 2D)
    VolumetricSrcSpec spec;
    spec.name = "cyl-end-U-plus-Fe-661keV";
    spec.geometry = GeometryType::CylinderEndOn;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 2.0*cm, 0.0}},
                    {"Fe (iron)", 0.0, 0.0, {1.0*cm, 1.0*cm, 0.0}} };
    fixtures.push_back( spec );
  }

  {// CylinderEndOn: source in outer layer, non-source inner core
    VolumetricSrcSpec spec;
    spec.name = "cyl-end-innervoid-661keV";
    spec.geometry = GeometryType::CylinderEndOn;
    spec.shells = { {"Fe (iron)", 0.0, 0.0, {1.5*cm, 1.5*cm, 0.0}},
                    {"U (uranium)", 0.0, 0.0, {1.0*cm, 1.0*cm, 0.0}} };
    spec.source_shell_index = 1;
    fixtures.push_back( spec );
  }

  {// CylinderSideOn: single iron cylinder (3D)
    VolumetricSrcSpec spec;
    spec.name = "cyl-side-Fe-661keV";
    spec.geometry = GeometryType::CylinderSideOn;
    spec.shells = { {"Fe (iron)", 0.0, 0.0, {5.0*cm, 10.0*cm, 0.0}} };
    fixtures.push_back( spec );
  }

  {// CylinderSideOn: U cylinder inside iron (3D)
    VolumetricSrcSpec spec;
    spec.name = "cyl-side-U-plus-Fe-661keV";
    spec.geometry = GeometryType::CylinderSideOn;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 5.0*cm, 0.0}},
                    {"Fe (iron)", 0.0, 0.0, {1.0*cm, 1.0*cm, 0.0}} };
    fixtures.push_back( spec );
  }

  {// Rectangular: single iron box (3D)
    VolumetricSrcSpec spec;
    spec.name = "rect-Fe-661keV";
    spec.geometry = GeometryType::Rectangular;
    spec.shells = { {"Fe (iron)", 0.0, 0.0, {5.0*cm, 5.0*cm, 5.0*cm}} };
    fixtures.push_back( spec );
  }

  {// Rectangular: U box inside iron (3D)
    VolumetricSrcSpec spec;
    spec.name = "rect-U-plus-Fe-661keV";
    spec.geometry = GeometryType::Rectangular;
    spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 2.0*cm, 2.0*cm}},
                    {"Fe (iron)", 0.0, 0.0, {1.0*cm, 1.0*cm, 1.0*cm}} };
    fixtures.push_back( spec );
  }

  {// Rectangular: source in outer layer, non-source inner core
    VolumetricSrcSpec spec;
    spec.name = "rect-innervoid-661keV";
    spec.geometry = GeometryType::Rectangular;
    spec.shells = { {"Fe (iron)", 0.0, 0.0, {1.5*cm, 1.5*cm, 1.5*cm}},
                    {"U (uranium)", 0.0, 0.0, {1.0*cm, 1.0*cm, 1.0*cm}} };
    spec.source_shell_index = 1;
    fixtures.push_back( spec );
  }

  {// Rectangular: polyethylene slab at 3 m - the analytic effective-AN/AD/H
   //  sanity case for the (future) effective-shielding accumulation
    VolumetricSrcSpec spec;
    spec.name = "rect-poly-slab-661keV";
    spec.geometry = GeometryType::Rectangular;
    spec.shells = { {"Polyethylene", 0.0, 0.0, {10.0*cm, 10.0*cm, 2.0*cm}} };
    spec.distance = 300.0*cm;
    fixtures.push_back( spec );
  }

  {// Rectangular: in-situ exponentially-distributed trace source (soil slab)
    VolumetricSrcSpec spec;
    spec.name = "rect-trace-exp-soil-661keV";
    spec.geometry = GeometryType::Rectangular;
    spec.shells = { {"Soil", 0.0, 0.0, {50.0*cm, 50.0*cm, 15.0*cm}} };
    spec.in_situ_exponential = true;
    spec.relaxation_length = 3.0*cm;
    fixtures.push_back( spec );
  }

  return fixtures;
}//standard_volumetric_fixtures()

}//namespace VolumetricFixture

#endif //VolumetricIntegrationFixtures_h
