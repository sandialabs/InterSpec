#pragma once
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

/// @file SpectralMeasure.h
/// @brief A mixed discrete/continuous representation of transported photons.

#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace ceelo::full_spectrum {

/// Physics label retained while training and validating an approximation.
enum class SpectralComponent : std::uint8_t {
    Uncollided,
    ElasticLine,
    SingleCompton,
    MultipleCompton,
    DoubleCompton,
    HigherCompton,
    ComptonRayleigh,
    SingleComptonRayleigh,
    DoubleComptonRayleigh,
    HigherComptonRayleigh,
    MultipleComptonRayleigh,
    PairAnnihilation,
    PairAnnihilationContinuum,
    PairAnnihilationInFlight,
    SourceElectronBrems,
    PhotoelectronBrems,
    ComptonElectronBrems,
    PairElectronBrems,
    PrimaryElectronBrems,
    PrimaryPositronBrems,
    PrimaryPositronAnnihilation,
    PrimaryPositronAnnihilationContinuum,
    PrimaryPositronAnnihilationInFlight,
    PrimaryPositronAnnihilationInFlightContinuum,
    PhotonuclearSecondary,
    DetectorZeroDeposit,
    DetectorFullEnergy,
    DetectorPairFullEnergy,
    DetectorSingleEscape,
    DetectorDoubleEscape,
    DetectorAnnihilation511,
    DetectorKAlphaEscape,
    DetectorKBetaEscape,
    DetectorKAlphaAbsorbed,
    DetectorKBetaAbsorbed,
    DetectorPhotoelectricEscape,
    DetectorPairNoCompton,
    DetectorPairSingleCompton,
    DetectorPairMultipleCompton,
    DetectorContinuum,
    /// Shield K and L x-rays escaping unscattered (source escape). The
    /// umbrella SourceFluorescence collects the three families.
    SourceFluorescence,
    SourceFluorescenceKAlpha,
    SourceFluorescenceKBeta,
    SourceFluorescenceL
};

template <typename Scalar>
struct SpectralAtom {
    double energy_keV = 0.0;
    Scalar yield = Scalar(0);
    SpectralComponent component = SpectralComponent::Uncollided;
};

/// An integral of a continuous density over [lower_keV, upper_keV).
///
/// The approximation evaluates a continuous function and integrates it over
/// caller-requested bins. These are not samples from a stored histogram.
template <typename Scalar>
struct SpectralBin {
    double lower_keV = 0.0;
    double upper_keV = 0.0;
    Scalar integral = Scalar(0);
    SpectralComponent component = SpectralComponent::DetectorContinuum;
};

/// A finite measure containing line-like atoms and a continuous density.
///
/// A source measure is an expected photon yield and can integrate above one
/// (for example, pair production can produce two annihilation photons). A
/// detector measure is normalized per photon incident on the detector envelope
/// and includes a zero-deposit atom.
template <typename Scalar>
struct SpectralMeasure {
    std::vector<SpectralAtom<Scalar>> atoms;
    std::vector<SpectralBin<Scalar>> continuous;

    Scalar integrated_yield() const {
        Scalar answer = Scalar(0);
        for (const auto& atom : atoms)
            answer += atom.yield;
        for (const auto& bin : continuous)
            answer += bin.integral;
        return answer;
    }

    Scalar component_yield(SpectralComponent component) const {
        const auto matches = [component](SpectralComponent candidate) {
            if (candidate == component)
                return true;
            if (component == SpectralComponent::SingleCompton) {
                return candidate == SpectralComponent::ComptonRayleigh
                    || candidate
                        == SpectralComponent::SingleComptonRayleigh;
            }
            if (component == SpectralComponent::DoubleCompton)
                return candidate
                    == SpectralComponent::DoubleComptonRayleigh;
            if (component == SpectralComponent::HigherCompton)
                return candidate
                    == SpectralComponent::HigherComptonRayleigh;
            if (component == SpectralComponent::PairAnnihilation)
                return candidate
                        == SpectralComponent::PairAnnihilationContinuum
                    || candidate
                        == SpectralComponent::PairAnnihilationInFlight;
            if (component == SpectralComponent::SourceFluorescence) {
                return candidate == SpectralComponent::SourceFluorescenceKAlpha
                    || candidate == SpectralComponent::SourceFluorescenceKBeta
                    || candidate == SpectralComponent::SourceFluorescenceL;
            }
            if (component == SpectralComponent::SourceElectronBrems) {
                return candidate == SpectralComponent::PhotoelectronBrems
                    || candidate
                        == SpectralComponent::ComptonElectronBrems
                    || candidate
                        == SpectralComponent::PairElectronBrems;
            }
            return component == SpectralComponent::MultipleCompton
                && (candidate == SpectralComponent::DoubleCompton
                    || candidate == SpectralComponent::HigherCompton
                    || candidate
                        == SpectralComponent::DoubleComptonRayleigh
                    || candidate
                        == SpectralComponent::HigherComptonRayleigh
                    || candidate
                        == SpectralComponent::MultipleComptonRayleigh);
        };
        Scalar answer = Scalar(0);
        for (const auto& atom : atoms)
            if (matches(atom.component))
                answer += atom.yield;
        for (const auto& bin : continuous)
            if (matches(bin.component))
                answer += bin.integral;
        return answer;
    }
};

inline void validate_bin_edges(const std::vector<double>& edges) {
    if (edges.size() < 2)
        throw std::invalid_argument("A spectral grid needs at least two edges");
    if (!std::isfinite(edges.front()))
        throw std::invalid_argument("Spectral bin edges must be finite");
    for (std::size_t i = 1; i < edges.size(); ++i) {
        if (!std::isfinite(edges[i]) || !(edges[i] > edges[i - 1]))
            throw std::invalid_argument(
                "Spectral bin edges must be finite and strictly increasing");
    }
}

} // namespace ceelo::full_spectrum
