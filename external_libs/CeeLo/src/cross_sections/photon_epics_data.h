/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC.
 This library is distributed under the GNU Lesser General Public License,
 version 2.1 or (at your option) any later version.
 */

#pragma once

#include <cstdint>

namespace ceelo {

/// One independently sampled photon cross-section curve. data_offset selects
/// parallel slices in g_photon_process_grid_index and
/// g_photon_process_log_value_q. Grid indices address the shared, bit-exact
/// float32 log10(MeV) coordinate pool. Values use a per-curve affine uint16
/// encoding; zero is the exact-zero/log-floor sentinel.
struct PhotonProcessCurve {
    uint16_t size;
    uint16_t data_offset;
    float value_offset;
    float value_scale;
};

/// Per-element affine value decoding for a Rayleigh total-cross-section row.
/// Four adjacent elements share each adaptive energy grid; group and row
/// offsets are derived from Z at runtime rather than repeated per element.
struct PhotonRayleighValues {
    float value_offset;
    float value_scale;
};

/// One independently sampled angular-factor curve. x is in inverse angstroms.
/// data_offset selects parallel slices in g_photon_angular_log_x and
/// g_photon_angular_log_value_q. Values are packed in log10 space. F is
/// interpolated in log-value; S is decoded at the bracket nodes and
/// interpolated in linear value.
struct PhotonAngularCurve {
    uint16_t size;
    uint16_t data_offset;
    float value_offset;
    float value_scale;
};

static_assert(sizeof(PhotonProcessCurve) == 12,
              "Photon process descriptors must remain compact");
static_assert(sizeof(PhotonRayleighValues) == 8,
              "Rayleigh value descriptors must remain compact");
static_assert(sizeof(PhotonAngularCurve) == 12,
              "Photon angular descriptors must remain compact");

struct PhotonEpicsElementData {
    PhotonRayleighValues rayleigh;        // EPDL MF=23 MT=502
    PhotonProcessCurve compton;           // EPDL MF=23 MT=504
    PhotonProcessCurve pair_production;   // EPDL MF=23 MT=516
    PhotonProcessCurve photoelectric;     // EPDL MF=23 MT=522
    PhotonProcessCurve k_photoelectric;   // EPDL MF=23 MT=534
    PhotonAngularCurve scattering_function; // EPDL MF=27 MT=504
};

static_assert(sizeof(PhotonEpicsElementData) == 68,
              "Per-element photon descriptors must remain compact");

extern const PhotonEpicsElementData g_photon_epics_data[92];
inline constexpr uint16_t kRayleighXsElementsPerGroup = 4;
inline constexpr uint16_t kRayleighXsGroups = 23;
extern const float g_rayleigh_log_energy[];
extern const uint16_t g_rayleigh_log_value_q[];
extern const uint16_t g_rayleigh_group_grid_offset[kRayleighXsGroups + 1];
extern const uint16_t g_rayleigh_group_value_offset[kRayleighXsGroups + 1];
extern const float g_photon_energy_pool[];
extern const uint16_t g_photon_process_grid_index[];
extern const uint16_t g_photon_process_log_value_q[];
extern const float g_photon_angular_log_x[];
extern const uint16_t g_photon_angular_log_value_q[];

/// Offline-generated normalized inverse-CDF data for coherent scattering.
/// Each row represents G(x)/G(x_max) on a shared logarithmic x grid; 65535 is
/// exactly one. Runtime interpolation is performed directly in the codes.
inline constexpr uint16_t kRayleighSamplingNodes = 256;
inline constexpr float kRayleighSamplingLogXMin = -3.0f;
inline constexpr float kRayleighSamplingLogXMax = 4.0f;
extern const float g_rayleigh_sampling_x[kRayleighSamplingNodes];
extern const uint16_t
    g_rayleigh_sampling_cdf_q[92][kRayleighSamplingNodes];

} // namespace ceelo
