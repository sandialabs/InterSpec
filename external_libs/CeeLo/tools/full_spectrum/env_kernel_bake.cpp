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

/// @file env_kernel_bake.cpp
/// @brief Offline generator for environment-scatter kernels, and the driver used
///        to inspect a single geometry.
///
/// This is the OFFLINE truth/bake tool for the air + environment scatter tier.
/// It is deliberately NOT wired into EfficiencyCalculator::compute(): running a
/// full analog environment inside the per-event efficiency loop would add cost
/// and an FEP-coupling hazard to the main path, for a peak effect already
/// covered by the analytic AirAttenuation modes. The production path is
/// analytic transmission plus the kernel this tool bakes.
///
/// WHAT IT PRODUCES. The photon flux arriving at the detector, resolved by
/// energy and incidence angle, and split by which surface of the detector
/// envelope it enters (front face / curved side / rear face). NO detector
/// efficiency is applied; detector response is folded in downstream.
///
/// With --fit-block it instead writes checksummed CEEFIT-BLOCK-v1 worker
/// records for the fit-first campaign driven by
/// tools/full_spectrum/env_scatter_fit.
///
/// Examples:
///   # One geometry, printed: 662 keV, 5 m, source and detector 1 m up, soil
///   env_kernel_bake --mode field --energies 662 --distances 500 \
///                   --h-src 100 --h-det 100 --histories 2000000
///
///   # Open air with a bounded halo, several distances, written to a table
///   env_kernel_bake --mode air --energies 59.5,662,2614 \
///                   --distances 100,500,2500 --extent 1000 \
///                   --histories 2000000 --out air_kernel.txt
///
///   # Cross-check the estimator against the analog truth on one node
///   env_kernel_bake --mode field --energies 662 --distances 40 \
///                   --h-src 25 --h-det 25 --histories 200000 --brute 40000000

#include "geometry/EnvironmentGeometry.h"
#include "io/EnvFitBlockBuilder.h"
#include "io/EnvScatterKernel.h"
#include "materials/Material.h"
#include "transport/NextEventEstimator.h"

#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace ceelo;

namespace {

std::vector<double> parse_list(const std::string& s)
{
    std::vector<double> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(std::atof(tok.c_str()));
    }
    return out;
}

std::string json_string(const std::string& value)
{
    std::ostringstream out;
    out << '"';
    static constexpr char hex[] = "0123456789abcdef";
    for (unsigned char c : value) {
        switch (c) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\b': out << "\\b"; break;
        case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (c < 0x20) {
                out << "\\u00" << hex[c >> 4] << hex[c & 0x0f];
            } else {
                out << static_cast<char>(c);
            }
        }
    }
    out << '"';
    return out.str();
}

std::string json_number(double value)
{
    if (!std::isfinite(value))
        throw std::invalid_argument("non-finite value in canonical configuration");
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return out.str();
}

std::string material_json(const Material& material)
{
    std::ostringstream out;
    out << "{\"composition\":[";
    for (std::size_t i = 0; i < material.composition().size(); ++i) {
        if (i != 0) out << ',';
        const auto& component = material.composition()[i];
        out << "{\"Z\":" << static_cast<unsigned>(component.Z)
            << ",\"mass_fraction\":" << json_number(component.mass_fraction)
            << '}';
    }
    out << "],\"density_g_cm3\":" << json_number(material.density())
        << ",\"name\":" << json_string(material.name()) << '}';
    return out.str();
}

std::string canonical_fit_config_json(const std::string& mode,
                                      const Material& air,
                                      const Material& slab,
                                      const EnvScatterConfig& cfg,
                                      int registered_bins)
{
    // Keys are emitted in lexical order. This exact byte string is hashed and
    // therefore must not depend on locale, indentation, or unordered maps.
    std::ostringstream out;
    out << "{\"air\":" << material_json(air)
        << ",\"axis\":{\"energy_bins\":" << registered_bins
        << ",\"energy_cut_keV\":" << json_number(cfg.axes.e_min_keV)
        << ",\"mu_bins\":" << cfg.axes.n_mu
        << ",\"transform\":\"registered-component-v1\"}"
        << ",\"block_histories\":262144"
        << ",\"geometry_tier\":" << json_string(mode)
        << ",\"material\":" << (mode == "air" ? "null" : material_json(slab))
        << ",\"physics\":{\"brems_electron_floor_keV\":"
        << json_number(cfg.brems_electron_floor_keV)
        << ",\"doppler\":" << (cfg.doppler ? "true" : "false")
        << ",\"free_electron_compton\":"
        << (cfg.free_electron_compton ? "true" : "false")
        << ",\"include_annihilation\":"
        << (cfg.include_annihilation ? "true" : "false")
        << ",\"include_brems\":" << (cfg.include_brems ? "true" : "false")
        << ",\"include_rayleigh\":" << (cfg.include_rayleigh ? "true" : "false")
        << ",\"max_vertices\":" << cfg.max_vertices
        << ",\"positron_range\":" << (cfg.positron_range ? "true" : "false")
        << ",\"roulette_tau\":" << json_number(cfg.roulette_tau)
        << ",\"walk_annihilation\":"
        << (cfg.walk_annihilation ? "true" : "false")
        << ",\"walk_brems\":" << (cfg.walk_brems ? "true" : "false")
        << ",\"weight_floor\":" << json_number(cfg.weight_floor) << '}'
        << ",\"scoring_radius_cm\":" << json_number(cfg.scoring_radius_cm)
        << ",\"schema\":\"ceelo-env-fit-config-v1\"}";
    return out.str();
}

std::string fit_block_output_path(const std::string& pattern,
                                  std::uint64_t block_id,
                                  std::uint64_t block_count)
{
    static const std::string token = "{block_id}";
    if (block_count == 1 && pattern.find(token) == std::string::npos)
        return pattern;
    if (pattern.find(token) == std::string::npos)
        throw std::invalid_argument(
            "--fit-block must contain {block_id} when --block-count > 1");
    std::string path = pattern;
    std::size_t position = 0;
    while ((position = path.find(token, position)) != std::string::npos) {
        const std::string value = std::to_string(block_id);
        path.replace(position, token.size(), value);
        position += value.size();
    }
    return path;
}

void usage()
{
    std::cout <<
"Usage: env_kernel_bake [options]\n"
"  --mode air|field|room   Geometry family (default field)\n"
"  --material soil|concrete  Slab material for --mode field (default soil)\n"
"  --energies E1,E2,...    Source energies in keV (default 662)\n"
"  --distances D1,D2,...   Source-detector distances in cm (default 500)\n"
"  --h-src H1,H2,...       Source heights above the slab, cm (default 100)\n"
"  --h-det H1,H2,...       Detector heights above the slab, cm (default 100)\n"
"  --extent X1,X2,...      Transverse air extent, cm; 0 = unbounded (air mode)\n"
"  --room H,W,BACK,FRONT   Room height, half-width, back, front extents in cm\n"
"                          (room mode; default 300,250,150,150)\n"
"  --histories N           Primary photons per node, or per BATCH when\n"
"                          --precision is given (default 1000000)\n"
"  --precision P           Per-node target relative sigma (e.g. 0.003); the\n"
"                          node runs batches until it is met. Off by default.\n"
"  --max-histories N       Cap for --precision (default 100x --histories)\n"
"  --landmark [N]          Per-node landmark-blocked energy axis with ~N bins\n"
"                          (default 192): edges pinned at the backscatter\n"
"                          edge, 511 keV and E0, so no bin straddles them\n"
"  --registered [N]        Fit-ready registered energy axis with N bins\n"
"                          (default 384): equal log blocks below/above the\n"
"                          backscatter edge; bin indices align across E0\n"
"  --shard I --nshards N   Bake only nodes with (index %% N) == I. Shards are\n"
"                          disjoint and seeds depend on the GLOBAL node index,\n"
"                          so a sharded bake equals the single-machine one.\n"
"  --binary-out FILE       Write the binary float32 kernel (the production\n"
"                          container; merge shards with env_kernel_merge)\n"
"  --fit-block FILE        Write checksummed CEEFIT-BLOCK-v1 scientific\n"
"                          records (requires exactly one geometry node)\n"
"  --block-id N            Nonzero continuation block id for --fit-block\n"
"  --block-count N         Run N consecutive blocks in this worker process;\n"
"                          FILE must then contain {block_id} (default 1)\n"
"  --role ROLE             pilot|train|validation|final|diagnostic\n"
"  --brute N              Also run N analog-truth histories as a cross-check\n"
"  --scoring-radius R      Scoring sphere radius, cm (default 4)\n"
"  --det-radius R          Detector envelope radius, cm (default 3.81)\n"
"  --det-length L          Detector envelope length, cm (default 7.62)\n"
"  --threads N             Worker threads (default: hardware concurrency)\n"
"  --no-rayleigh           Drop Rayleigh scattering from the estimator\n"
"  --no-positron-range     Annihilate the e+ at the pair-production vertex\n"
"  --no-walk-annihilation  Score the 511 pair but do not walk it (A/B)\n"
"  --no-brems              Drop the bremsstrahlung channel (A/B)\n"
"  --no-walk-brems         Score brems emissions but do not walk them (A/B)\n"
"  --brems-floor KE        Electron KE floor for brems, keV (default 150)\n"
"  --e-floor KE            Transport/tally energy floor, keV (default 30;\n"
"                          fit campaigns use 10)\n"
"  --mu-bins N             Uniform incidence-angle bins (default 32; fit\n"
"                          pilot uses 64)\n"
"  --max-vertices N        Cap scatter order (1 = single scatter only)\n"
"  --out FILE              Write the baked kernel table to FILE\n"
"  --spectrum FILE         Write the per-energy-bin flux spectrum to FILE\n";
}

/// Bake one node, multi-threaded, and return the merged result.
///
/// With `precision > 0` the node runs in batches of `histories` until the relative sigma of its total meets that target (or
/// `max_histories` is exhausted); otherwise exactly one batch runs.
///
/// Precision targeting exists because a FIXED history count does not buy a
/// fixed accuracy: at 4M histories the same field geometry lands at 0.41%
/// relative sigma at 50 cm but 1.27% at 50 m, so a uniform-history bake
/// delivers a product whose quality silently degrades along its own grid --
/// and it is the far nodes, where interpolation is hardest, that come out
/// worst. Time per history is flat with distance (measured 21.8/21.0/21.3 s
/// per 4M at 0.5/5/50 m), so equalizing precision simply moves histories to
/// where they are needed.
EnvScatterResult run_node(const EnvironmentGeometry& env,
                          const EnvScatterConfig& cfg,
                          const Eigen::Vector3d& src, double e0_keV,
                          uint64_t histories, int threads, bool brute,
                          uint64_t node_index, double precision = 0.0,
                          uint64_t max_histories = 0)
{
    EnvScatterResult total(cfg.axes, e0_keV);
    if (threads < 1) threads = 1;
    if (histories == 0) return total;

    // Construct ONCE on this thread. The estimator is const and keeps no
    // mutable state across run() (each call owns its RNG), so one instance is
    // safely shared. Constructing it here also means its validation (e.g. the
    // scoring sphere overlapping a slab) surfaces as a catchable error instead
    // of an uncaught exception inside a worker thread, which would terminate.
    const NextEventEstimator est(env, cfg, src);

    for (uint64_t round = 0; ; ++round) {
        const uint64_t done = total.tally.num_histories();
        uint64_t round_histories = histories;
        if (max_histories > 0) {
            if (done >= max_histories) break;
            round_histories = std::min(round_histories, max_histories - done);
        }
        const int active_threads = static_cast<int>(std::min<uint64_t>(
            static_cast<uint64_t>(threads), round_histories));
        std::vector<EnvScatterResult> parts(threads,
                                            EnvScatterResult(cfg.axes, e0_keV));
        std::vector<std::thread> pool;
        pool.reserve(active_threads);
        const uint64_t base = round_histories /
                              static_cast<uint64_t>(active_threads);
        const uint64_t remainder = round_histories %
                                   static_cast<uint64_t>(active_threads);
        for (int t = 0; t < active_threads; ++t) {
            const uint64_t thread_histories = base +
                (static_cast<uint64_t>(t) < remainder ? 1u : 0u);
            pool.emplace_back([&, t, thread_histories]() {
                // Salt the seed with the run MODE, the node index AND the
                // batch round. Without the mode salt the estimator and the
                // analog cross-check would consume byte-identical RNG
                // streams, so their early histories would share emission
                // directions and flight draws -- and the printed z assumes an
                // independence it would not have. Without the round salt the
                // precision loop would re-run the SAME histories and drive
                // sigma down without adding information.
                const uint64_t seed = 0x9E3779B97F4A7C15ULL * (t + 1) +
                                      (brute ? 0xB5297A4DU : 0x68E31DA4U) +
                                      0x9E3779B1ULL * node_index +
                                      0x2545F4914F6CDD1DULL * round;
                if (brute) est.run_brute(e0_keV, thread_histories, seed, parts[t]);
                else       est.run(e0_keV, thread_histories, seed, parts[t]);
            });
        }
        for (auto& th : pool) th.join();

        for (const auto& p : parts) {
            total.tally.merge(p.tally);
            total.histories += p.histories;
            total.direct_transmission = p.direct_transmission;
            total.direct_fluence = p.direct_fluence;
            total.source_distance_cm = p.source_distance_cm;
        }

        if (precision <= 0.0) break;
        const double m = total.tally.total_mean();
        const double s = total.tally.total_sigma();
        if (m > 0.0 && s / m <= precision) break;
        if (max_histories > 0 && total.tally.num_histories() >= max_histories)
            break;
    }
    return total;
}

void print_flux(const char* label, const SurfaceFlux& f, const SurfaceFlux& s)
{
    std::cout << "    " << std::left << std::setw(10) << label << std::right
              << "  face " << std::scientific << std::setprecision(4) << f.face
              << " +/- " << s.face
              << "   side " << f.side << " +/- " << s.side
              << "   back " << f.back << " +/- " << s.back
              << "   total " << f.total() << "\n";
}

} // namespace

static int run_main(int argc, char** argv);

int main(int argc, char** argv)
{
    try {
        return run_main(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}

static int run_main(int argc, char** argv)
{
    std::string mode = "field", material = "soil", out_path, spectrum_path;
    std::string binary_out, fit_block_path;
    std::string fit_role = "train";
    uint64_t fit_block_id = 1, fit_block_count = 1;
    double precision = 0.0;
    uint64_t max_histories = 0;
    int landmark_bins = 0;      // 0 = log-uniform table-level axis
    int registered_bins = 0;    // 0 = disabled; fit-ready E0-registered axis
    int shard = 0, nshards = 1;
    std::vector<double> energies{662.0}, distances{500.0};
    std::vector<double> h_src{100.0}, h_det{100.0}, extents{0.0};
    uint64_t histories = 1000000, brute_histories = 0;
    bool histories_explicit = false, e_floor_explicit = false;
    bool mu_bins_explicit = false;
    double scoring_radius = 4.0;
    DetectorEnvelope det;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    bool rayleigh = true;
    bool positron_range = true, walk_ann = true, free_e = false;
    bool include_brems = true, walk_brems = true;
    double brems_floor = -1.0;   // <0 = library default
    double e_floor = 30.0;
    int mu_bins = 32;
    int max_vtx = 0;   // 0 = library default
    std::vector<double> room_dims{300.0, 250.0, 150.0, 150.0};

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            return (i + 1 < argc) ? argv[++i] : std::string();
        };
        if (a == "--mode") mode = next();
        else if (a == "--material") material = next();
        else if (a == "--energies") energies = parse_list(next());
        else if (a == "--distances") distances = parse_list(next());
        else if (a == "--h-src") h_src = parse_list(next());
        else if (a == "--h-det") h_det = parse_list(next());
        else if (a == "--extent") extents = parse_list(next());
        else if (a == "--histories") {
            histories = std::strtoull(next().c_str(), nullptr, 10);
            histories_explicit = true;
        }
        else if (a == "--brute") brute_histories = std::strtoull(next().c_str(), nullptr, 10);
        else if (a == "--scoring-radius") scoring_radius = std::atof(next().c_str());
        else if (a == "--det-radius") det.radius_cm = std::atof(next().c_str());
        else if (a == "--det-length") det.length_cm = std::atof(next().c_str());
        else if (a == "--threads") threads = std::atoi(next().c_str());
        else if (a == "--room") room_dims = parse_list(next());
        else if (a == "--no-rayleigh") rayleigh = false;
        else if (a == "--free-electron") free_e = true;
        else if (a == "--max-vertices") max_vtx = std::atoi(next().c_str());
        else if (a == "--no-positron-range") positron_range = false;
        else if (a == "--no-walk-annihilation") walk_ann = false;
        else if (a == "--no-brems") include_brems = false;
        else if (a == "--no-walk-brems") walk_brems = false;
        else if (a == "--brems-floor") brems_floor = std::atof(next().c_str());
        else if (a == "--e-floor") {
            e_floor = std::atof(next().c_str());
            e_floor_explicit = true;
        }
        else if (a == "--mu-bins") {
            mu_bins = std::atoi(next().c_str());
            mu_bins_explicit = true;
        }
        else if (a == "--out") out_path = next();
        else if (a == "--binary-out") binary_out = next();
        else if (a == "--fit-block") fit_block_path = next();
        else if (a == "--block-id")
            fit_block_id = std::strtoull(next().c_str(), nullptr, 10);
        else if (a == "--block-count")
            fit_block_count = std::strtoull(next().c_str(), nullptr, 10);
        else if (a == "--role") fit_role = next();
        else if (a == "--precision") precision = std::atof(next().c_str());
        else if (a == "--max-histories")
            max_histories = std::strtoull(next().c_str(), nullptr, 10);
        else if (a == "--landmark") {
            // Optional numeric argument.
            if (i + 1 < argc && argv[i + 1][0] != '-')
                landmark_bins = std::atoi(next().c_str());
            else
                landmark_bins = 192;
        }
        else if (a == "--registered") {
            // Optional numeric argument.
            if (i + 1 < argc && argv[i + 1][0] != '-')
                registered_bins = std::atoi(next().c_str());
            else
                registered_bins = 384;
        }
        else if (a == "--shard") shard = std::atoi(next().c_str());
        else if (a == "--nshards") nshards = std::atoi(next().c_str());
        else if (a == "--spectrum") spectrum_path = next();
        else if (a == "--help" || a == "-h") { usage(); return 0; }
        else { std::cerr << "Unknown option: " << a << "\n"; usage(); return 1; }
    }

    const std::vector<std::string> valid_modes{"air", "field", "room"};
    if (std::find(valid_modes.begin(), valid_modes.end(), mode) ==
        valid_modes.end()) {
        std::cerr << "--mode must be air, field, or room\n";
        return 1;
    }
    if (material != "soil" && material != "concrete") {
        std::cerr << "--material must be soil or concrete\n";
        return 1;
    }

    const bool fit_mode = !fit_block_path.empty();
    if (fit_mode) {
        if (!histories_explicit) histories = 262144;
        if (!e_floor_explicit) e_floor = 10.0;
        if (!mu_bins_explicit) mu_bins = 64;
        if (registered_bins == 0 && landmark_bins == 0) registered_bins = 384;
        if (fit_block_id == 0 || fit_block_count == 0 ||
            fit_block_count - 1 >
                std::numeric_limits<std::uint64_t>::max() - fit_block_id) {
            std::cerr << "--block-id/count must define a nonzero, non-overflowing "
                         "range\n";
            return 1;
        }
        if (fit_block_count > 1 &&
            fit_block_path.find("{block_id}") == std::string::npos) {
            std::cerr << "--fit-block must contain {block_id} when "
                         "--block-count > 1\n";
            return 1;
        }
        if (histories != 262144) {
            std::cerr << "--fit-block uses exactly 262144 histories; continuation "
                         "is expressed with --block-id\n";
            return 1;
        }
        if (registered_bins != 384) {
            std::cerr << "--fit-block uses exactly 384 registered energy bins "
                         "(192 below and 192 above the component landmark)\n";
            return 1;
        }
        const std::vector<std::string> roles{
            "pilot", "train", "validation", "final", "diagnostic"};
        if (std::find(roles.begin(), roles.end(), fit_role) == roles.end()) {
            std::cerr << "invalid --role for --fit-block\n";
            return 1;
        }
        const std::size_t node_count = energies.size() * distances.size() *
            h_src.size() * h_det.size() * extents.size();
        if (node_count != 1) {
            std::cerr << "--fit-block requires exactly one requested node\n";
            return 1;
        }
        if (mode == "room") {
            std::cerr << "--fit-block initial scope is air or field, not room\n";
            return 1;
        }
        if (!out_path.empty() || !binary_out.empty() ||
            brute_histories != 0 || max_histories != 0) {
            std::cerr << "--fit-block cannot be combined with runtime-table, "
                         "brute-force, or max-history outputs\n";
            return 1;
        }
        if (landmark_bins > 0) {
            std::cerr << "--fit-block requires the registered axis, not --landmark\n";
            return 1;
        }
    }

    const Material air = make_Air();
    const Material slab = (material == "concrete") ? make_Concrete() : make_Soil();
    const bool is_air_mode = (mode == "air");
    const bool is_room_mode = (mode == "room");
    if (is_room_mode && room_dims.size() != 4) {
        std::cerr << "--room needs exactly 4 values: H,HALF_WIDTH,BACK,FRONT\n";
        return 1;
    }

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = scoring_radius;
    if (!(e_floor >= 10.0)) {
        std::cerr << "--e-floor must be at least 10 keV\n";
        return 1;
    }
    if (mu_bins < 1) {
        std::cerr << "--mu-bins must be positive\n";
        return 1;
    }
    cfg.e_floor_keV = e_floor;
    cfg.axes.e_min_keV = e_floor;
    cfg.axes.n_mu = mu_bins;
    cfg.include_rayleigh = rayleigh;
    cfg.positron_range = positron_range;
    cfg.walk_annihilation = walk_ann;
    cfg.include_brems = include_brems;
    cfg.walk_brems = walk_brems;
    if (brems_floor >= 0.0) cfg.brems_electron_floor_keV = brems_floor;
    cfg.free_electron_compton = free_e;
    if (max_vtx > 0) cfg.max_vertices = max_vtx;

    EnvKernelTable table;
    table.kind = is_air_mode ? "air" : (is_room_mode ? ("room_" + material) : material);
    table.air_name = air.name();
    table.air_density = air.density();
    if (!is_air_mode) { table.slab_name = slab.name(); table.slab_density = slab.density(); }
    table.scoring_radius_cm = scoring_radius;
    table.axes = cfg.axes;
    table.e0_grid_keV = energies;
    table.distance_grid_cm = distances;
    table.h_src_grid_cm = is_air_mode ? std::vector<double>{} : h_src;
    table.h_det_grid_cm = is_air_mode ? std::vector<double>{} : h_det;
    table.extent_grid_cm = is_air_mode ? extents : std::vector<double>{};
    table.provenance = "env_kernel_bake; incident flux only, no detector efficiency";

    std::cout << "Environment-scatter bake: mode=" << mode
              << (is_air_mode ? "" : (" material=" + material))
              << " histories/node=" << histories
              << " threads=" << threads << "\n"
              << "Detector envelope: R=" << det.radius_cm
              << " cm, L=" << det.length_cm << " cm; scoring sphere R="
              << scoring_radius << " cm\n"
              << "NOTE: fluxes below are INCIDENT photons per source photon.\n"
              << "      No detector efficiency is applied.\n\n";

    const std::vector<double>& hs = is_air_mode ? std::vector<double>{0.0} : h_src;
    const std::vector<double>& hd = is_air_mode ? std::vector<double>{0.0} : h_det;
    const std::vector<double>& ex = is_air_mode ? extents : std::vector<double>{0.0};

    if (nshards < 1 || shard < 0 || shard >= nshards) {
        std::cerr << "--shard must be in [0, --nshards)\n";
        return 1;
    }
    if (landmark_bins > 0 && registered_bins > 0) {
        std::cerr << "--landmark and --registered are mutually exclusive\n";
        return 1;
    }
    if (precision > 0.0 && max_histories == 0) max_histories = 100 * histories;

    // GLOBAL node index: increments over every grid point whether or not this
    // shard bakes it, so seeds (and therefore results) are independent of how
    // the work was split -- a 1-shard bake and an 8-shard bake produce the
    // same numbers.
    uint64_t global_index = 0;
    int n_degenerate = 0, n_baked = 0;

    for (double e0 : energies)
    for (double d : distances)
    for (double xs : ex)
    for (double hsrc : hs)
    for (double hdet : hd) {
        const uint64_t node_index = global_index++;
        if (static_cast<int>(node_index % nshards) != shard) continue;

        const Eigen::Vector3d src(0.0, 0.0, -d);
        EnvironmentGeometry env =
            is_air_mode
                ? (xs > 0.0 ? EnvironmentGeometry::open_air(air, src,
                                  Eigen::Vector3d::Zero(), xs)
                            : EnvironmentGeometry(air))
            : is_room_mode
                ? EnvironmentGeometry::room(air, slab, src, hsrc, hdet,
                                            room_dims[0], room_dims[1],
                                            room_dims[2], room_dims[3])
                : EnvironmentGeometry::field(air, slab, src, hsrc, hdet);

        if (!is_air_mode && env.num_slabs() == 0) {
            // EnvironmentGeometry::field returns a slab-free environment on
            // degenerate input (zero separation, or a height difference larger
            // than the separation). Baking that would silently label an open-air
            // node as a field node. On a full grid sweep these corners are
            // EXPECTED (the geometry axes are coupled: |h_src - h_det| < d is
            // a validity constraint, not a rectangle), so skip and count them
            // -- env_kernel_is_complete() applies the same rule when checking
            // a merged table, so a skipped corner is not reported as a hole.
            ++n_degenerate;
            std::cout << "  SKIP degenerate (d=" << d << " h_src=" << hsrc
                      << " h_det=" << hdet << "): |h_src-h_det| >= d\n";
            continue;
        }

        // Per-node energy axis in landmark mode: the backscatter edge and the
        // upper limit are properties of THIS node's E0.
        EnvScatterConfig ncfg = cfg;
        if (registered_bins > 0)
            ncfg.axes = make_registered_landmark_axes(
                e0, registered_bins, ncfg.axes.n_mu, ncfg.axes.e_min_keV);
        else if (landmark_bins > 0)
            ncfg.axes = make_landmark_axes(e0, landmark_bins, ncfg.axes.n_mu,
                                           ncfg.axes.e_min_keV);

        if (fit_mode) {
            EnvFitManifest manifest;
            manifest.canonical_config_json = canonical_fit_config_json(
                mode, air, slab, ncfg, registered_bins);
            manifest.config_sha256 =
                env_fit_sha256_hex(manifest.canonical_config_json);
            manifest.material_name = is_air_mode ? air.name() : slab.name();
            manifest.material_definition_json =
                is_air_mode ? material_json(air) : material_json(slab);
            manifest.geometry_tier = mode;
            manifest.rng_algorithm = "mt19937_64/sha256-tuple-stream-v1";
            manifest.estimator_flags =
                "next-event;rayleigh=" + std::string(rayleigh ? "on" : "off") +
                ";brems=" + (include_brems ? "on" : "off") +
                ";annihilation=" + (ncfg.include_annihilation ? "on" : "off");
            manifest.units = "energy=keV,length=cm,fluence=cm^-2";
            manifest.observable =
                "bin-integrated scalar incident fluence per source photon";
            manifest.normalization =
                "raw S1/S2; s_hat divisor=T(E0)/(4*pi*d^2)";

            EnvFitBlockRequest request;
            request.manifest = std::move(manifest);
            request.role = fit_role;
            request.e0_keV = e0;
            request.distance_cm = d;
            request.h_src_cm = hsrc;
            request.h_det_cm = hdet;
            request.transverse_extent_cm = xs;
            request.histories = histories;
            request.threads = threads;
            request.target_total_rel_se =
                precision > 0.0 ? precision : 0.0025;

            const auto t0 = std::chrono::steady_clock::now();
            for (std::uint64_t offset = 0; offset < fit_block_count; ++offset) {
                request.block_id = fit_block_id + offset;
                EnvFitBlock block = run_env_fit_block(env, ncfg, src, request);
                const std::string output_path = fit_block_output_path(
                    fit_block_path, request.block_id, fit_block_count);
                std::string write_error;
                if (!save_env_fit_block_atomic(block, output_path, &write_error))
                    throw std::runtime_error(
                        "failed to write fit block: " + write_error);
                std::cout << "CEEFIT block node=" << std::hex << block.node_id
                          << std::dec << " block=" << block.block_id
                          << " histories=" << block.node.histories
                          << " rel_se=" << std::scientific
                          << block.node.achieved_total_rel_se
                          << " worker_wall=" << std::fixed
                          << std::setprecision(1) << block.node.wall_seconds
                          << " s\n"
                          << "Wrote fit block to " << output_path << "\n";
            }
            ++n_baked;
            const double secs = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            std::cout << "Completed " << fit_block_count
                      << " fit block(s) in " << std::fixed
                      << std::setprecision(1) << secs << " s\n";
            continue;
        }

        const auto t0 = std::chrono::steady_clock::now();
        EnvScatterResult r = run_node(env, ncfg, src, e0, histories, threads,
                                      false, node_index, precision,
                                      max_histories);
        ++n_baked;
        const double secs = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();

        std::cout << "E0=" << e0 << " keV  d=" << d << " cm";
        if (is_air_mode) { if (xs > 0.0) std::cout << "  extent=" << xs << " cm"; }
        else std::cout << "  h_src=" << hsrc << " h_det=" << hdet << " cm";
        std::cout << "   [" << std::fixed << std::setprecision(1) << secs
                  << " s, " << r.tally.num_histories() << " histories]\n";

        std::cout << std::scientific << std::setprecision(4)
                  << "    T(E0) = " << r.direct_transmission
                  << "   direct fluence = " << r.direct_fluence << " /cm2\n";

        // The direct beam arrives along +z (source on the -z axis), so mu_in = 1.
        const SurfaceFlux direct = direct_surface_flux(r.direct_fluence, 1.0, det);
        const SurfaceFlux scat = surface_flux(r.tally, det);
        const SurfaceFlux scat_s = surface_flux_sigma(r.tally, det);
        print_flux("direct", direct, SurfaceFlux{});
        print_flux("scattered", scat, scat_s);
        // The SCALAR fluence is what the kernel stores and what a GEANT4
        // track-length scorer (phi = L/V) measures, so print it explicitly --
        // it is the apples-to-apples validation observable. The surface fluxes
        // above are this field projected onto the envelope.
        std::cout << "    scalar fluence: scattered "
                  << r.tally.total_mean() << " +/- " << r.tally.total_sigma()
                  << " /cm2   (direct " << r.direct_fluence << ")\n";
        // Line channels: discrete-energy arrivals kept out of the continuum.
        for (int l = 0; l < r.tally.n_lines(); ++l) {
            const double lm = r.tally.line_total_mean(l);
            if (lm <= 0.0) continue;
            std::cout << "    line " << r.tally.line(l).energy_keV << " keV ("
                      << (r.tally.line(l).origin == LineOrigin::ElasticE0
                              ? "elastic-E0" : "annihilation")
                      << "): " << lm << " +/- " << r.tally.line_total_sigma(l)
                      << " /cm2   = " << 100.0 * lm / std::max(1e-300, r.tally.total_mean())
                      << "% of scattered\n";
        }
        std::cout << "    scatter/direct (total surfaces) = "
                  << (direct.total() > 0.0 ? scat.total() / direct.total() : 0.0)
                  << "\n";

        // Origin decomposition: reads directly against the G4 harness's
        // ORIGIN DECOMPOSITION block (same creator-with-chain-inheritance
        // semantics on both sides).
        {
            const double m_all = std::max(1e-300, r.tally.total_mean());
            const double m_pri =
                r.tally.total_mean_origin(FluxOrigin::Primary);
            const double m_ann =
                r.tally.total_mean_origin(FluxOrigin::Annihilation);
            const double m_brm = r.tally.total_mean_origin(FluxOrigin::Brems);
            std::cout << "    origins: primary " << m_pri << " ("
                      << 100.0 * m_pri / m_all << "%)  annihil " << m_ann
                      << " +/- " << r.tally.total_sigma_origin(
                                        FluxOrigin::Annihilation)
                      << " (" << 100.0 * m_ann / m_all << "%)  brems "
                      << m_brm << " +/- "
                      << r.tally.total_sigma_origin(FluxOrigin::Brems)
                      << " (" << 100.0 * m_brm / m_all << "%)\n";
        }

        const double m_air = r.tally.total_mean_class(ChainClass::AirOnly);
        const double m_one = r.tally.total_mean_class(ChainClass::OneBounce);
        const double m_same = r.tally.total_mean_class(ChainClass::MultiSameSlab);
        const double m_dist = r.tally.total_mean_class(ChainClass::MultiDistinct);
        const double m_mul = m_same + m_dist;
        const double m_all = std::max(1e-300, r.tally.total_mean());
        std::cout << "    chain classes: air-only " << (m_air / m_all)
                  << "  1-bounce " << (m_one / m_all)
                  << "  multi-same " << (m_same / m_all)
                  << "  multi-distinct " << (m_dist / m_all) << "\n";
        if (r.tally.dropped_weight() > 0.0) {
            std::cout << "    NOTE: weight below the "
                      << ncfg.axes.e_min_keV << " keV tally floor: "
                      << r.tally.dropped_weight() / r.tally.num_histories()
                      << " per history\n";
        }

        if (brute_histories > 0) {
            EnvScatterConfig bcfg = ncfg;
            // The analog truth tallies the flight LEAVING a vertex, so it needs
            // one more vertex of budget to cover the same physics.
            bcfg.max_vertices = ncfg.max_vertices + 1;
            EnvScatterResult b = run_node(env, bcfg, src, e0, brute_histories,
                                          threads, true, table.nodes.size());
            const double mn = r.tally.total_mean(), sn = r.tally.total_sigma();
            const double mb = b.tally.total_mean(), sb = b.tally.total_sigma();
            const double z = (mn - mb) / std::sqrt(sn * sn + sb * sb);
            std::cout << "    ANALOG CROSS-CHECK: next-event " << mn << " +/- " << sn
                      << ", brute " << mb << " +/- " << sb
                      << "  (ratio " << std::fixed << std::setprecision(4)
                      << (mb > 0.0 ? mn / mb : 0.0)
                      << ", z = " << std::setprecision(2) << z << ")\n"
                      << std::scientific;
        }

        // Record the node, normalized to the direct fluence.
        EnvKernelNode node;
        node.e0_keV = e0;
        node.distance_cm = d;
        node.h_src_cm = hsrc;
        node.h_det_cm = hdet;
        node.transverse_extent_cm = xs;
        node.direct_transmission = r.direct_transmission;
        node.direct_fluence = r.direct_fluence;
        node.histories = r.tally.num_histories();
        node.multi_slab_fraction = m_dist / m_all;
        // Landmark mode: the axis belongs to the node, so it travels with it.
        node.energy_edges = ncfg.axes.edges;
        const double scale = (r.direct_fluence > 0.0) ? 1.0 / r.direct_fluence : 0.0;
        node.s_hat.reserve(ncfg.axes.n_bins());
        node.s_hat_slab.reserve(ncfg.axes.n_bins());
        double worst = 0.0;
        for (int ie = 0; ie < ncfg.axes.n_energy; ++ie) {
            for (int j = 0; j < ncfg.axes.n_mu; ++j) {
                const double m = r.tally.mean(ie, j);
                node.s_hat.push_back(m * scale);
                node.s_hat_slab.push_back(
                    (r.tally.mean_class(ie, j, ChainClass::OneBounce) +
                     r.tally.mean_class(ie, j, ChainClass::MultiSameSlab) +
                     r.tally.mean_class(ie, j, ChainClass::MultiDistinct)) * scale);
                if (m > 0.0) worst = std::max(worst, r.tally.sigma(ie, j) / m);
            }
        }
        // Persist the line channels alongside the continuum grid.
        for (int l = 0; l < r.tally.n_lines(); ++l) {
            EnvKernelNode::Line L;
            L.energy_keV = r.tally.line(l).energy_keV;
            L.origin = static_cast<uint8_t>(r.tally.line(l).origin);
            for (int j = 0; j < ncfg.axes.n_mu; ++j) {
                L.s_hat.push_back(r.tally.line_mean(l, j) * scale);
                L.s_hat_slab.push_back(
                    (r.tally.line_mean_class(l, j, ChainClass::OneBounce) +
                     r.tally.line_mean_class(l, j, ChainClass::MultiSameSlab) +
                     r.tally.line_mean_class(l, j, ChainClass::MultiDistinct)) * scale);
            }
            node.lines.push_back(std::move(L));
        }
        node.max_rel_sigma = worst;
        table.nodes.push_back(std::move(node));

        if (!spectrum_path.empty()) {
            // One file per node: writing them all to `spectrum_path` would leave
            // only the last node's spectrum behind on a multi-node sweep.
            std::string sp_name = spectrum_path;
            if (table.nodes.size() > 0 || energies.size() * distances.size() > 1) {
                const std::size_t dot = spectrum_path.rfind('.');
                const std::string stem = (dot == std::string::npos)
                                           ? spectrum_path
                                           : spectrum_path.substr(0, dot);
                const std::string ext = (dot == std::string::npos)
                                          ? std::string(".csv")
                                          : spectrum_path.substr(dot);
                std::ostringstream tag;
                tag << stem << "_node" << table.nodes.size() << ext;
                sp_name = tag.str();
            }
            std::ofstream sp(sp_name);
            sp << "# Incident scattered flux spectrum (photons per source photon)\n"
               << "# E0: " << e0 << " keV, distance " << d << " cm\n"
               << "# No detector efficiency applied.\n"
               << "Energy (keV),Face,Side,Back\n";
            for (int ie = 0; ie < ncfg.axes.n_energy; ++ie) {
                const SurfaceFlux f = surface_flux_in_bin(r.tally, det, ie);
                sp << ncfg.axes.energy_low(ie) << "," << f.face << "," << f.side
                   << "," << f.back << "\n";
            }
            std::cout << "    wrote spectrum to " << sp_name << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "Baked " << n_baked << " node(s)";
    if (nshards > 1)
        std::cout << " (shard " << shard << " of " << nshards << ", "
                  << global_index << " grid points total)";
    if (n_degenerate > 0)
        std::cout << "; skipped " << n_degenerate
                  << " degenerate |h_src-h_det| >= d corner(s)";
    std::cout << "\n";
    if (n_baked == 0 && n_degenerate > 0) {
        std::cerr << "ERROR: every requested node was degenerate; check that "
                     "|h_src - h_det| < d\n";
        return 1;
    }

    if (!out_path.empty()) {
        if (save_env_kernel(table, out_path))
            std::cout << "Wrote kernel table (" << table.nodes.size()
                      << " nodes) to " << out_path << "\n";
        else
            std::cerr << "FAILED to write " << out_path << "\n";
    }
    if (!binary_out.empty()) {
        if (save_env_kernel_binary(table, binary_out)) {
            std::ifstream fb(binary_out, std::ios::binary | std::ios::ate);
            std::cout << "Wrote BINARY kernel (" << table.nodes.size()
                      << " nodes, " << (fb ? fb.tellg() / 1024 : 0)
                      << " kB) to " << binary_out << "\n";
        } else {
            std::cerr << "FAILED to write " << binary_out << "\n";
        }
    }
    return 0;
}
