/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
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

// Export GDMLs for extended/shielded cascade-summing GEANT4 cross-checks.
// Both geometries also absorb the Co-60 betas (shield / source self-absorption),
// removing the beta-gamma-summing confound of a bare source.
//
//   A) shielded point source : NaI 3x3, Co-60 point at -2 cm, 0.2 cm Fe shell
//   B) extended source       : NaI 3x3, Co-60 in a r=1,hz=1 cm water cylinder at -4 cm
//
// Vacuum world (no air) to match compute_cascade, which does not transport the
// source-to-detector air gap.

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"

#include <cstdio>

using namespace ceelo;

int main() {
    static Material nai = make_NaI();
    static Material fe = make_Iron();
    static Material water = make_Water();

    {   // A) shielded point source
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -2.0));
        calc.add_source_shield(&fe, 0.2);  // 0.2 cm Fe shell (stops betas, ~8% gamma atten)
        calc.export_geant4_gdml("detector_casc_shield.gdml", /*vacuum_world=*/true);
        printf("A) detector_casc_shield.gdml  : NaI 3x3, point -2 cm, 0.2 cm Fe shell\n");
    }
    {   // B) extended water cylinder source
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -4.0), 1.0, 1.0);
        calc.set_source_material(&water);
        calc.export_geant4_gdml("detector_casc_extended.gdml", /*vacuum_world=*/true);
        printf("B) detector_casc_extended.gdml: NaI 3x3, water cyl r=1 hz=1 at -4 cm\n");
    }
    return 0;
}
