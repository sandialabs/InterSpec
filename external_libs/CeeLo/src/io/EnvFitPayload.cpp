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

#include "io/EnvFitPayload.h"

#include "io/EnvFitArchive.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

#if defined(_WIN32)
#include <fcntl.h>
#include <io.h>
#include <process.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

namespace ceelo {
namespace {

constexpr std::array<std::uint8_t, 16> kMagic{{
    'C', 'E', 'E', 'F', 'I', 'T', '-', 'P',
    'A', 'Y', 'L', 'O', 'A', 'D', '1', '\0'
}};
constexpr std::uint32_t kEndianTag = 0x01020304U;
constexpr std::uint64_t kHeaderBytes = 64;
constexpr std::uint64_t kMaximumPayloadBytes = 1024ULL * 1024ULL * 1024ULL;
constexpr std::uint64_t kMaximumStringBytes = 16ULL * 1024ULL * 1024ULL;
constexpr std::uint64_t kMaximumVectorElements = 1ULL << 28;
constexpr unsigned kSplineDegree = 3;

bool fail(std::string* error, const std::string& message)
{
    if (error) *error = message;
    return false;
}

bool durable_sync_file(const std::string& path)
{
#if defined(_WIN32)
    const int descriptor = _open(path.c_str(), _O_RDWR | _O_BINARY);
    if (descriptor < 0) return false;
    const bool ok = _commit(descriptor) == 0;
    return _close(descriptor) == 0 && ok;
#else
    const int descriptor = open(path.c_str(), O_RDONLY);
    if (descriptor < 0) return false;
    const bool ok = fsync(descriptor) == 0;
    return close(descriptor) == 0 && ok;
#endif
}

bool durable_sync_parent_directory(const std::string& path)
{
#if defined(_WIN32)
    (void)path;
    return true;
#else
    const std::size_t separator = path.find_last_of('/');
    const std::string parent = separator == std::string::npos
        ? "." : separator == 0 ? "/" : path.substr(0, separator);
    const int descriptor = open(parent.c_str(), O_RDONLY);
    if (descriptor < 0) return false;
    const bool ok = fsync(descriptor) == 0;
    return close(descriptor) == 0 && ok;
#endif
}

bool finite(double value)
{
    return std::isfinite(value);
}

bool checked_product(std::size_t a, std::size_t b, std::size_t& result)
{
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) return false;
    result = a * b;
    return true;
}

std::size_t feature_count(std::uint32_t level)
{
    if (level < 1U || level > 4U) return 0U;
    const std::size_t n = (std::size_t{1} << (level - 1U)) + 3U;
    return n * n * n * n;
}

bool is_continuum_kind(EnvFitAxisKind kind)
{
    return kind == EnvFitAxisKind::Primary ||
           kind == EnvFitAxisKind::Annihilation ||
           kind == EnvFitAxisKind::Brems;
}

bool is_line_kind(EnvFitAxisKind kind)
{
    return kind == EnvFitAxisKind::LineE0 || kind == EnvFitAxisKind::Line511;
}

bool canonical_model_identity(const std::string& name, EnvFitAxisKind kind,
                              bool line)
{
    const auto match = [&](const char* expected_name,
                           EnvFitAxisKind expected_kind,
                           bool expected_line) {
        return name == expected_name && kind == expected_kind &&
               line == expected_line;
    };
    return match("primary_order1", EnvFitAxisKind::Primary, false) ||
           match("primary_order2plus", EnvFitAxisKind::Primary, false) ||
           match("annihilation_order1", EnvFitAxisKind::Annihilation, false) ||
           match("annihilation_order2plus", EnvFitAxisKind::Annihilation, false) ||
           match("brems_order0", EnvFitAxisKind::Brems, false) ||
           match("brems_order1", EnvFitAxisKind::Brems, false) ||
           match("brems_order2plus", EnvFitAxisKind::Brems, false) ||
           match("elastic_E0", EnvFitAxisKind::LineE0, true) ||
           match("annihilation_511", EnvFitAxisKind::Line511, true);
}

bool valid_sha256(const std::string& value)
{
    if (value.size() != 64) return false;
    return std::all_of(value.begin(), value.end(), [](char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    });
}

bool valid_knots(const std::vector<double>& knots)
{
    if (knots.size() < 2U * (kSplineDegree + 1U)) return false;
    if (!std::all_of(knots.begin(), knots.end(), finite)) return false;
    for (std::size_t i = 1; i < knots.size(); ++i)
        if (knots[i] < knots[i - 1]) return false;
    if (knots.front() != 0.0 || knots.back() != 1.0) return false;
    for (unsigned i = 0; i <= kSplineDegree; ++i) {
        if (knots[i] != 0.0 || knots[knots.size() - 1U - i] != 1.0)
            return false;
    }
    const std::size_t n_basis = knots.size() - kSplineDegree - 1U;
    for (std::size_t i = 0; i < n_basis; ++i)
        if (!(knots[i + kSplineDegree + 1U] > knots[i])) return false;
    return true;
}

bool valid_edges(const std::vector<double>& edges)
{
    if (edges.size() < 2 || !std::all_of(edges.begin(), edges.end(), finite))
        return false;
    for (std::size_t i = 1; i < edges.size(); ++i)
        if (!(edges[i] > edges[i - 1])) return false;
    return true;
}

bool validate_model(const EnvFitRuntimeModel& model, bool line,
                    double brems_activation_threshold_keV,
                    std::string* error)
{
    const std::string prefix = "model '" + model.name + "': ";
    if (model.name.empty()) return fail(error, "runtime model has an empty name");
    if (line ? !is_line_kind(model.axis_kind) : !is_continuum_kind(model.axis_kind))
        return fail(error, prefix + "axis kind does not match model collection");
    if (!canonical_model_identity(model.name, model.axis_kind, line))
        return fail(error, prefix + "name/axis is not in the canonical runtime taxonomy");
    if (!finite(model.e_cut_keV) || model.e_cut_keV <= 0.0)
        return fail(error, prefix + "Ecut must be finite and positive");
    if (model.has_energy_threshold) {
        if (!finite(model.energy_threshold_keV) || model.energy_threshold_keV <= 0.0)
            return fail(error, prefix + "activation threshold must be positive");
    } else if (model.energy_threshold_keV != 0.0) {
        return fail(error, prefix + "disabled activation threshold must be zero");
    }
    if (model.rank == 0 || model.level < 1 || model.level > 4)
        return fail(error, prefix + "rank/geometry spline level is invalid");
    const bool annihilation = model.axis_kind == EnvFitAxisKind::Annihilation ||
                              model.axis_kind == EnvFitAxisKind::Line511;
    if (annihilation &&
        (!model.has_energy_threshold ||
         std::abs(model.energy_threshold_keV -
                  2.0 * kEnvFitElectronRestEnergyKeV) > 1.0e-9))
        return fail(error, prefix +
            "annihilation support must begin at the exact pair threshold");
    if ((model.axis_kind == EnvFitAxisKind::Primary ||
         model.axis_kind == EnvFitAxisKind::LineE0) &&
        model.has_energy_threshold)
        return fail(error, prefix + "primary support must not have an activation threshold");
    if (model.axis_kind == EnvFitAxisKind::Brems) {
        const bool expected_flag = brems_activation_threshold_keV > 0.0;
        if (model.has_energy_threshold != expected_flag ||
            (expected_flag &&
             std::abs(model.energy_threshold_keV -
                      brems_activation_threshold_keV) > 1.0e-9))
            return fail(error, prefix +
                "brems support disagrees with the configured activation threshold");
    }
    if (model.axis_kind == EnvFitAxisKind::Annihilation &&
        model.e_cut_keV >= kEnvFitElectronRestEnergyKeV)
        return fail(error, prefix + "annihilation Ecut must be below 511 keV");
    if (!valid_knots(model.mu_knots))
        return fail(error, prefix + "mu knots are not clamped cubic knots on [0,1]");
    if (line) {
        if (!model.energy_knots.empty())
            return fail(error, prefix + "line model must not carry energy knots");
    } else if (!valid_knots(model.energy_knots)) {
        return fail(error, prefix + "energy knots are not clamped cubic knots on [0,1]");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        if (!finite(model.q_min[i]) || !finite(model.q_span[i]) ||
            model.q_span[i] <= 0.0)
            return fail(error, prefix + "q domain is non-finite or empty");
    }

    const std::size_t features = feature_count(model.level);
    std::size_t logits = 0;
    if (!checked_product(features, model.rank, logits) ||
        model.log_yield_coefficients.size() != features ||
        model.logit_coefficients.size() != logits)
        return fail(error, prefix + "geometry coefficient dimensions do not match rank/level");
    if (!std::all_of(model.log_yield_coefficients.begin(),
                     model.log_yield_coefficients.end(), finite) ||
        !std::all_of(model.logit_coefficients.begin(),
                     model.logit_coefficients.end(), finite))
        return fail(error, prefix + "geometry coefficients contain non-finite values");

    const std::size_t n_mu_basis = model.mu_knots.size() - kSplineDegree - 1U;
    const std::size_t n_energy_basis = line
        ? 1U : model.energy_knots.size() - kSplineDegree - 1U;
    std::size_t archetype_size = 0;
    if (!checked_product(n_energy_basis, n_mu_basis, archetype_size) ||
        !checked_product(archetype_size, model.rank, archetype_size) ||
        model.archetype_coefficients.size() != archetype_size)
        return fail(error, prefix + "archetype dimensions do not match knots/rank");
    const std::size_t per_rank = n_energy_basis * n_mu_basis;
    for (std::size_t r = 0; r < model.rank; ++r) {
        long double sum = 0.0L;
        for (std::size_t i = 0; i < per_rank; ++i) {
            const double value = model.archetype_coefficients[r * per_rank + i];
            if (!finite(value) || value < 0.0)
                return fail(error, prefix + "archetype coefficients must be finite/nonnegative");
            sum += value;
        }
        if (std::abs(sum - 1.0L) > 1.0e-10L)
            return fail(error, prefix + "each archetype coefficient row must sum to one");
    }
    return true;
}

class Writer {
public:
    const std::vector<std::uint8_t>& bytes() const { return bytes_; }

    void u8(std::uint8_t value) { bytes_.push_back(value); }
    void u32(std::uint32_t value)
    {
        for (unsigned i = 0; i < 4; ++i)
            bytes_.push_back(static_cast<std::uint8_t>(value >> (8U * i)));
    }
    void u64(std::uint64_t value)
    {
        for (unsigned i = 0; i < 8; ++i)
            bytes_.push_back(static_cast<std::uint8_t>(value >> (8U * i)));
    }
    void f64(double value)
    {
        std::uint64_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "double must be binary64");
        std::memcpy(&bits, &value, sizeof(bits));
        u64(bits);
    }
    void raw(const void* data, std::size_t size)
    {
        const auto* p = static_cast<const std::uint8_t*>(data);
        bytes_.insert(bytes_.end(), p, p + size);
    }
    void string(const std::string& value)
    {
        u64(static_cast<std::uint64_t>(value.size()));
        raw(value.data(), value.size());
    }
    void doubles(const std::vector<double>& values)
    {
        u64(static_cast<std::uint64_t>(values.size()));
        for (double value : values) f64(value);
    }

private:
    std::vector<std::uint8_t> bytes_;
};

class Reader {
public:
    Reader(const std::uint8_t* data, std::size_t size)
        : data_(data), size_(size) {}

    bool done() const { return position_ == size_; }
    std::uint8_t u8()
    {
        require(1);
        return data_[position_++];
    }
    std::uint32_t u32()
    {
        require(4);
        std::uint32_t value = 0;
        for (unsigned i = 0; i < 4; ++i)
            value |= static_cast<std::uint32_t>(data_[position_++]) << (8U * i);
        return value;
    }
    std::uint64_t u64()
    {
        require(8);
        std::uint64_t value = 0;
        for (unsigned i = 0; i < 8; ++i)
            value |= static_cast<std::uint64_t>(data_[position_++]) << (8U * i);
        return value;
    }
    double f64()
    {
        const std::uint64_t bits = u64();
        double value = 0.0;
        std::memcpy(&value, &bits, sizeof(value));
        return value;
    }
    void raw(void* output, std::size_t size)
    {
        require(size);
        std::memcpy(output, data_ + position_, size);
        position_ += size;
    }
    std::string string()
    {
        const std::uint64_t count = u64();
        if (count > kMaximumStringBytes || count > size_ - position_)
            throw std::runtime_error("invalid payload string length");
        std::string value(reinterpret_cast<const char*>(data_ + position_),
                          static_cast<std::size_t>(count));
        position_ += static_cast<std::size_t>(count);
        return value;
    }
    std::vector<double> doubles()
    {
        const std::uint64_t count = u64();
        if (count > kMaximumVectorElements ||
            count > (size_ - position_) / sizeof(double))
            throw std::runtime_error("invalid payload coefficient-vector length");
        std::vector<double> result(static_cast<std::size_t>(count));
        for (double& value : result) value = f64();
        return result;
    }

private:
    void require(std::size_t count)
    {
        if (count > size_ - position_)
            throw std::runtime_error("truncated fitted payload");
    }

    const std::uint8_t* data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t position_ = 0;
};

std::array<std::uint8_t, 32> digest_bytes(const void* data, std::size_t size)
{
    const std::string hex = env_fit_sha256_hex(data, size);
    std::array<std::uint8_t, 32> result{};
    const auto nibble = [](char c) -> std::uint8_t {
        if (c >= '0' && c <= '9') return static_cast<std::uint8_t>(c - '0');
        return static_cast<std::uint8_t>(10 + c - 'a');
    };
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = static_cast<std::uint8_t>((nibble(hex[2 * i]) << 4U) |
                                              nibble(hex[2 * i + 1]));
    return result;
}

void serialize_model(const EnvFitRuntimeModel& model, Writer& writer)
{
    writer.string(model.name);
    writer.u32(static_cast<std::uint32_t>(model.axis_kind));
    writer.f64(model.e_cut_keV);
    writer.u8(model.has_energy_threshold ? 1U : 0U);
    writer.f64(model.energy_threshold_keV);
    writer.u32(model.rank);
    writer.u32(model.level);
    writer.doubles(model.energy_knots);
    writer.doubles(model.mu_knots);
    for (double value : model.q_min) writer.f64(value);
    for (double value : model.q_span) writer.f64(value);
    writer.doubles(model.log_yield_coefficients);
    writer.doubles(model.logit_coefficients);
    writer.doubles(model.archetype_coefficients);
}

EnvFitRuntimeModel deserialize_model(Reader& reader)
{
    EnvFitRuntimeModel model;
    model.name = reader.string();
    model.axis_kind = static_cast<EnvFitAxisKind>(reader.u32());
    model.e_cut_keV = reader.f64();
    const std::uint8_t threshold = reader.u8();
    if (threshold > 1U) throw std::runtime_error("invalid threshold flag");
    model.has_energy_threshold = threshold != 0U;
    model.energy_threshold_keV = reader.f64();
    model.rank = reader.u32();
    model.level = reader.u32();
    model.energy_knots = reader.doubles();
    model.mu_knots = reader.doubles();
    for (double& value : model.q_min) value = reader.f64();
    for (double& value : model.q_span) value = reader.f64();
    model.log_yield_coefficients = reader.doubles();
    model.logit_coefficients = reader.doubles();
    model.archetype_coefficients = reader.doubles();
    return model;
}

std::vector<const EnvFitRuntimeModel*> sorted_models(
    const std::vector<EnvFitRuntimeModel>& models)
{
    std::vector<const EnvFitRuntimeModel*> result;
    result.reserve(models.size());
    for (const auto& model : models) result.push_back(&model);
    std::sort(result.begin(), result.end(),
              [](const auto* a, const auto* b) { return a->name < b->name; });
    return result;
}

std::vector<const EnvFitModelInventoryEntry*> sorted_inventory(
    const std::vector<EnvFitModelInventoryEntry>& inventory)
{
    std::vector<const EnvFitModelInventoryEntry*> result;
    result.reserve(inventory.size());
    for (const auto& entry : inventory) result.push_back(&entry);
    std::sort(result.begin(), result.end(), [](const auto* a, const auto* b) {
        return std::tie(a->line, a->name) < std::tie(b->line, b->name);
    });
    return result;
}

void serialize_inventory(const EnvFitPayload& payload, Writer& writer)
{
    writer.u64(static_cast<std::uint64_t>(payload.model_inventory.size()));
    for (const auto* entry : sorted_inventory(payload.model_inventory)) {
        writer.u8(entry->line ? 1U : 0U);
        writer.string(entry->name);
        writer.u32(static_cast<std::uint32_t>(entry->axis_kind));
    }
}

void serialize_model_content(const EnvFitPayload& payload, Writer& writer)
{
    writer.string(payload.config_sha256);
    writer.string(payload.material_name);
    writer.string(payload.geometry_tier);
    writer.u32(static_cast<std::uint32_t>(payload.normalization));
    writer.f64(payload.brems_activation_threshold_keV);
    for (double value : payload.physical_domain_min) writer.f64(value);
    for (double value : payload.physical_domain_max) writer.f64(value);
    serialize_inventory(payload, writer);
    writer.u64(static_cast<std::uint64_t>(payload.continuum_models.size()));
    for (const auto* model : sorted_models(payload.continuum_models))
        serialize_model(*model, writer);
    writer.u64(static_cast<std::uint64_t>(payload.line_models.size()));
    for (const auto* model : sorted_models(payload.line_models))
        serialize_model(*model, writer);
}

std::string model_content_sha256(const EnvFitPayload& payload)
{
    Writer writer;
    serialize_model_content(payload, writer);
    return env_fit_sha256_hex(writer.bytes().data(), writer.bytes().size());
}

void serialize_payload(const EnvFitPayload& payload, Writer& writer)
{
    writer.string(payload.config_sha256);
    writer.string(payload.material_name);
    writer.string(payload.geometry_tier);
    writer.u32(static_cast<std::uint32_t>(payload.normalization));
    writer.u8(payload.accepted_for_publication ? 1U : 0U);
    writer.f64(payload.brems_activation_threshold_keV);
    writer.string(model_content_sha256(payload));
    for (double value : payload.physical_domain_min) writer.f64(value);
    for (double value : payload.physical_domain_max) writer.f64(value);
    serialize_inventory(payload, writer);
    writer.string(payload.validation_metadata_json);
    writer.u64(static_cast<std::uint64_t>(payload.continuum_models.size()));
    for (const auto* model : sorted_models(payload.continuum_models))
        serialize_model(*model, writer);
    writer.u64(static_cast<std::uint64_t>(payload.line_models.size()));
    for (const auto* model : sorted_models(payload.line_models))
        serialize_model(*model, writer);
}

EnvFitPayload deserialize_payload(const std::uint8_t* data, std::size_t size)
{
    Reader reader(data, size);
    EnvFitPayload payload;
    payload.config_sha256 = reader.string();
    payload.material_name = reader.string();
    payload.geometry_tier = reader.string();
    payload.normalization = static_cast<EnvFitPayloadNormalization>(reader.u32());
    const std::uint8_t publication = reader.u8();
    if (publication > 1U) throw std::runtime_error("invalid publication flag");
    payload.accepted_for_publication = publication != 0U;
    payload.brems_activation_threshold_keV = reader.f64();
    payload.model_content_sha256 = reader.string();
    for (double& value : payload.physical_domain_min) value = reader.f64();
    for (double& value : payload.physical_domain_max) value = reader.f64();
    const std::uint64_t n_inventory = reader.u64();
    if (n_inventory > 4096U) throw std::runtime_error("too many inventory entries");
    payload.model_inventory.reserve(static_cast<std::size_t>(n_inventory));
    for (std::uint64_t i = 0; i < n_inventory; ++i) {
        const std::uint8_t line = reader.u8();
        if (line > 1U) throw std::runtime_error("invalid inventory collection flag");
        EnvFitModelInventoryEntry entry;
        entry.line = line != 0U;
        entry.name = reader.string();
        entry.axis_kind = static_cast<EnvFitAxisKind>(reader.u32());
        payload.model_inventory.push_back(std::move(entry));
    }
    payload.validation_metadata_json = reader.string();
    const std::uint64_t n_continuum = reader.u64();
    if (n_continuum > 4096U) throw std::runtime_error("too many continuum models");
    payload.continuum_models.reserve(static_cast<std::size_t>(n_continuum));
    for (std::uint64_t i = 0; i < n_continuum; ++i)
        payload.continuum_models.push_back(deserialize_model(reader));
    const std::uint64_t n_lines = reader.u64();
    if (n_lines > 4096U) throw std::runtime_error("too many line models");
    payload.line_models.reserve(static_cast<std::size_t>(n_lines));
    for (std::uint64_t i = 0; i < n_lines; ++i)
        payload.line_models.push_back(deserialize_model(reader));
    if (!reader.done()) throw std::runtime_error("trailing bytes in fitted payload");
    return payload;
}

std::vector<double> bspline_values(const std::vector<double>& knots, double x)
{
    const std::size_t n_basis = knots.size() - kSplineDegree - 1U;
    x = std::max(0.0, std::min(1.0, x));
    std::vector<double> basis(knots.size() - 1U, 0.0);
    for (std::size_t i = 0; i + 1U < knots.size(); ++i)
        if (x >= knots[i] && x < knots[i + 1U]) basis[i] = 1.0;
    if (x == 1.0) basis[n_basis - 1U] = 1.0;
    for (unsigned order = 1; order <= kSplineDegree; ++order) {
        std::vector<double> next(basis.size() - 1U, 0.0);
        for (std::size_t i = 0; i < next.size(); ++i) {
            const double left_denom = knots[i + order] - knots[i];
            const double right_denom = knots[i + order + 1U] - knots[i + 1U];
            if (left_denom > 0.0)
                next[i] += (x - knots[i]) / left_denom * basis[i];
            if (right_denom > 0.0)
                next[i] += (knots[i + order + 1U] - x) / right_denom *
                           basis[i + 1U];
        }
        basis = std::move(next);
    }
    basis.resize(n_basis);
    return basis;
}

/// Exact (to roundoff) interval masses of individually normalized cubic
/// B-spline PDFs.  Two-point Gauss-Legendre is exact for each cubic piece; the
/// interval is split at every knot before quadrature.
std::vector<double> basis_interval_masses(const std::vector<double>& knots,
                                          double low, double high)
{
    const std::size_t n_basis = knots.size() - kSplineDegree - 1U;
    std::vector<double> result(n_basis, 0.0);
    low = std::max(0.0, std::min(1.0, low));
    high = std::max(0.0, std::min(1.0, high));
    if (high <= low) return result;

    std::vector<double> cuts{low};
    for (double knot : knots)
        if (knot > low && knot < high && knot != cuts.back()) cuts.push_back(knot);
    cuts.push_back(high);
    constexpr double inv_sqrt_three = 0.577350269189625764509148780501957;
    for (std::size_t segment = 0; segment + 1U < cuts.size(); ++segment) {
        const double midpoint = 0.5 * (cuts[segment] + cuts[segment + 1U]);
        const double half = 0.5 * (cuts[segment + 1U] - cuts[segment]);
        const auto left = bspline_values(knots, midpoint - half * inv_sqrt_three);
        const auto right = bspline_values(knots, midpoint + half * inv_sqrt_three);
        for (std::size_t i = 0; i < n_basis; ++i)
            result[i] += half * (left[i] + right[i]);
    }
    for (std::size_t i = 0; i < n_basis; ++i) {
        const double norm = (knots[i + kSplineDegree + 1U] - knots[i]) /
                            static_cast<double>(kSplineDegree + 1U);
        result[i] = norm > 0.0 ? result[i] / norm : 0.0;
        if (std::abs(result[i]) < 1.0e-15) result[i] = 0.0;
    }
    return result;
}

std::vector<double> geometry_features(const EnvFitRuntimeModel& model,
                                      const std::array<double, 4>& q,
                                      bool& inside)
{
    const std::size_t intervals = std::size_t{1} << (model.level - 1U);
    const std::size_t n_basis = intervals + 3U;
    std::array<std::vector<double>, 4> one_dim;
    inside = true;
    for (std::size_t dim = 0; dim < 4; ++dim) {
        double value = (q[dim] - model.q_min[dim]) / model.q_span[dim];
        constexpr double tolerance = 1.0e-12;
        if (value < -tolerance || value > 1.0 + tolerance) inside = false;
        value = std::max(0.0, std::min(1.0, value));
        std::vector<double> knots(kSplineDegree + 1U, 0.0);
        for (std::size_t i = 1; i < intervals; ++i)
            knots.push_back(static_cast<double>(i) /
                            static_cast<double>(intervals));
        knots.insert(knots.end(), kSplineDegree + 1U, 1.0);
        one_dim[dim] = bspline_values(knots, value);
    }
    std::vector<double> features;
    features.reserve(n_basis * n_basis * n_basis * n_basis);
    for (std::size_t i0 = 0; i0 < n_basis; ++i0)
        for (std::size_t i1 = 0; i1 < n_basis; ++i1)
            for (std::size_t i2 = 0; i2 < n_basis; ++i2)
                for (std::size_t i3 = 0; i3 < n_basis; ++i3)
                    features.push_back(one_dim[0][i0] * one_dim[1][i1] *
                                       one_dim[2][i2] * one_dim[3][i3]);
    return features;
}

bool model_yield_mixture(const EnvFitRuntimeModel& model,
                         const std::array<double, 4>& q,
                         double& yield,
                         std::vector<double>& mixture,
                         std::string* error)
{
    bool inside = false;
    const std::vector<double> features = geometry_features(model, q, inside);
    if (!inside)
        return fail(error, "query lies outside q domain of model '" + model.name + "'");
    double log_yield = 0.0;
    for (std::size_t f = 0; f < features.size(); ++f)
        log_yield += features[f] * model.log_yield_coefficients[f];
    yield = std::exp(log_yield);
    if (!finite(yield) || yield <= 0.0)
        return fail(error, "non-finite fitted yield in model '" + model.name + "'");

    mixture.assign(model.rank, 0.0);
    for (std::size_t r = 0; r < model.rank; ++r)
        for (std::size_t f = 0; f < features.size(); ++f)
            mixture[r] += features[f] * model.logit_coefficients[f * model.rank + r];
    const double maximum = *std::max_element(mixture.begin(), mixture.end());
    double sum = 0.0;
    for (double& value : mixture) {
        value = std::exp(value - maximum);
        sum += value;
    }
    if (!finite(sum) || sum <= 0.0)
        return fail(error, "non-finite fitted softmax in model '" + model.name + "'");
    for (double& value : mixture) value /= sum;
    return true;
}

double component_endpoint(const EnvFitRuntimeModel& model, double e0_keV)
{
    switch (model.axis_kind) {
    case EnvFitAxisKind::Primary:
    case EnvFitAxisKind::Brems:
        return e0_keV;
    case EnvFitAxisKind::Annihilation:
        return kEnvFitElectronRestEnergyKeV;
    default:
        return 0.0;
    }
}

double physical_to_registered(const EnvFitRuntimeModel& model,
                              double energy_keV, double e0_keV)
{
    const double endpoint = component_endpoint(model, e0_keV);
    if (energy_keV <= model.e_cut_keV) return 0.0;
    if (energy_keV >= endpoint) return 1.0;
    if (model.axis_kind == EnvFitAxisKind::Brems)
        return std::log(energy_keV / model.e_cut_keV) /
               std::log(endpoint / model.e_cut_keV);
    const double landmark = model.axis_kind == EnvFitAxisKind::Primary
        ? e0_keV / (1.0 + 2.0 * e0_keV / kEnvFitElectronRestEnergyKeV)
        : kEnvFitElectronRestEnergyKeV / 3.0;
    if (energy_keV <= landmark)
        return 0.5 * std::log(energy_keV / model.e_cut_keV) /
               std::log(landmark / model.e_cut_keV);
    return 0.5 + 0.5 * std::log(energy_keV / landmark) /
                       std::log(endpoint / landmark);
}

std::vector<std::vector<double>> requested_energy_integrals(
    const EnvFitRuntimeModel& model, double e0_keV,
    const std::vector<double>& requested_edges)
{
    const std::size_t n_basis = model.energy_knots.size() - kSplineDegree - 1U;
    std::vector<std::vector<double>> result(requested_edges.size() - 1U,
                                            std::vector<double>(n_basis, 0.0));
    const double endpoint = component_endpoint(model, e0_keV);
    if (!(endpoint > model.e_cut_keV)) return result;
    for (std::size_t bin = 0; bin + 1U < requested_edges.size(); ++bin) {
        const double low = std::max(requested_edges[bin], model.e_cut_keV);
        const double high = std::min(requested_edges[bin + 1U], endpoint);
        if (high > low)
            result[bin] = basis_interval_masses(
                model.energy_knots,
                physical_to_registered(model, low, e0_keV),
                physical_to_registered(model, high, e0_keV));
    }
    return result;
}

std::vector<std::vector<double>> requested_mu_integrals(
    const std::vector<double>& knots, const std::vector<double>& requested_edges)
{
    const std::size_t n_basis = knots.size() - kSplineDegree - 1U;
    std::vector<std::vector<double>> result(requested_edges.size() - 1U,
                                            std::vector<double>(n_basis, 0.0));
    for (std::size_t bin = 0; bin + 1U < requested_edges.size(); ++bin) {
        const double low = std::max(-1.0, requested_edges[bin]);
        const double high = std::min(1.0, requested_edges[bin + 1U]);
        if (high > low)
            result[bin] = basis_interval_masses(knots, 0.5 * (low + 1.0),
                                                0.5 * (high + 1.0));
    }
    return result;
}

bool inactive_at_energy(const EnvFitRuntimeModel& model, double e0_keV)
{
    return model.has_energy_threshold && e0_keV < model.energy_threshold_keV;
}

bool evaluate_continuum_model(const EnvFitRuntimeModel& model,
                              const std::array<double, 4>& q,
                              double e0_keV,
                              const std::vector<double>& energy_edges,
                              const std::vector<double>& mu_edges,
                              std::vector<double>& output,
                              std::string* error)
{
    const std::size_t n_e_out = energy_edges.size() - 1U;
    const std::size_t n_mu_out = mu_edges.size() - 1U;
    output.assign(n_e_out * n_mu_out, 0.0);
    if (inactive_at_energy(model, e0_keV)) return true;

    double yield = 0.0;
    std::vector<double> mixture;
    if (!model_yield_mixture(model, q, yield, mixture, error)) return false;
    const std::size_t n_e_basis = model.energy_knots.size() - kSplineDegree - 1U;
    const std::size_t n_mu_basis = model.mu_knots.size() - kSplineDegree - 1U;
    const std::size_t per_rank = n_e_basis * n_mu_basis;
    std::vector<double> coefficients(per_rank, 0.0);
    for (std::size_t r = 0; r < model.rank; ++r)
        for (std::size_t cell = 0; cell < per_rank; ++cell)
            coefficients[cell] += mixture[r] *
                model.archetype_coefficients[r * per_rank + cell];

    const auto e_integrals = requested_energy_integrals(model, e0_keV, energy_edges);
    const auto mu_integrals = requested_mu_integrals(model.mu_knots, mu_edges);
    for (std::size_t de = 0; de < n_e_out; ++de) {
        for (std::size_t dm = 0; dm < n_mu_out; ++dm) {
            long double mass = 0.0L;
            for (std::size_t se = 0; se < n_e_basis; ++se) {
                if (e_integrals[de][se] == 0.0) continue;
                for (std::size_t sm = 0; sm < n_mu_basis; ++sm)
                    mass += static_cast<long double>(coefficients[se * n_mu_basis + sm]) *
                            e_integrals[de][se] * mu_integrals[dm][sm];
            }
            output[de * n_mu_out + dm] = yield * static_cast<double>(mass);
        }
    }
    return true;
}

bool evaluate_line_model(const EnvFitRuntimeModel& model,
                         const std::array<double, 4>& q,
                         double e0_keV,
                         const std::vector<double>& mu_edges,
                         EnvFitEvaluatedLine& output,
                         std::string* error)
{
    double yield = 0.0;
    std::vector<double> mixture;
    if (!model_yield_mixture(model, q, yield, mixture, error)) return false;
    const std::size_t n_mu_basis = model.mu_knots.size() - kSplineDegree - 1U;
    std::vector<double> coefficients(n_mu_basis, 0.0);
    for (std::size_t r = 0; r < model.rank; ++r)
        for (std::size_t m = 0; m < n_mu_basis; ++m)
            coefficients[m] += mixture[r] *
                model.archetype_coefficients[r * n_mu_basis + m];
    const auto mu_integrals = requested_mu_integrals(model.mu_knots, mu_edges);
    output.name = model.name;
    output.energy_keV = model.axis_kind == EnvFitAxisKind::LineE0
        ? e0_keV : kEnvFitElectronRestEnergyKeV;
    output.s_hat_mu_bin_masses.assign(mu_edges.size() - 1U, 0.0);
    for (std::size_t dm = 0; dm < output.s_hat_mu_bin_masses.size(); ++dm) {
        long double mass = 0.0L;
        for (std::size_t sm = 0; sm < n_mu_basis; ++sm)
            mass += static_cast<long double>(coefficients[sm]) * mu_integrals[dm][sm];
        output.s_hat_mu_bin_masses[dm] = yield * static_cast<double>(mass);
    }
    return true;
}

bool evaluate_valid_payload(const EnvFitPayload& payload,
                            const EnvFitRuntimeQuery& query,
                            const std::vector<double>& energy_edges_keV,
                            const std::vector<double>& mu_edges,
                            EnvFitEvaluation& result,
                            std::string* error)
{
    if (!valid_edges(energy_edges_keV) || energy_edges_keV.front() < 0.0)
        return fail(error, "requested energy edges must be finite, increasing, and nonnegative");
    if (!valid_edges(mu_edges) || mu_edges.front() < -1.0 || mu_edges.back() > 1.0)
        return fail(error, "requested mu edges must be finite, increasing, and within [-1,1]");
    if (!finite(query.direct_fluence_per_cm2) || query.direct_fluence_per_cm2 <= 0.0)
        return fail(error, "runtime query requires positive analytic direct fluence");
    const std::array<double, 4> physical{{
        query.e0_keV, query.distance_m,
        query.source_height_m, query.detector_height_m
    }};
    for (std::size_t i = 0; i < physical.size(); ++i) {
        const double scale = std::max(
            {1.0, std::abs(payload.physical_domain_min[i]),
             std::abs(payload.physical_domain_max[i])});
        const double tolerance = 1.0e-12 * scale;
        if (physical[i] < payload.physical_domain_min[i] - tolerance ||
            physical[i] > payload.physical_domain_max[i] + tolerance)
            return fail(error, "runtime query lies outside the payload physical domain");
    }
    std::array<double, 4> q{};
    if (!env_fit_transform_query(query, q, error)) return false;

    EnvFitEvaluation evaluated;
    evaluated.n_energy = static_cast<std::uint32_t>(energy_edges_keV.size() - 1U);
    evaluated.n_mu = static_cast<std::uint32_t>(mu_edges.size() - 1U);
    evaluated.direct_fluence_per_cm2 = query.direct_fluence_per_cm2;
    evaluated.continuum_s_hat_cell_masses.assign(
        static_cast<std::size_t>(evaluated.n_energy) * evaluated.n_mu, 0.0);
    evaluated.incident_continuum_cell_masses_per_cm2.assign(
        evaluated.continuum_s_hat_cell_masses.size(), 0.0);
    evaluated.components.reserve(payload.continuum_models.size());
    for (const auto& model : payload.continuum_models) {
        EnvFitEvaluatedComponent component;
        component.name = model.name;
        if (!evaluate_continuum_model(model, q, query.e0_keV,
                                      energy_edges_keV, mu_edges,
                                      component.s_hat_cell_masses, error))
            return false;
        component.incident_cell_masses_per_cm2.resize(component.s_hat_cell_masses.size());
        for (std::size_t i = 0; i < component.s_hat_cell_masses.size(); ++i) {
            evaluated.continuum_s_hat_cell_masses[i] +=
                component.s_hat_cell_masses[i];
            component.incident_cell_masses_per_cm2[i] =
                component.s_hat_cell_masses[i] * query.direct_fluence_per_cm2;
            evaluated.incident_continuum_cell_masses_per_cm2[i] +=
                component.incident_cell_masses_per_cm2[i];
        }
        evaluated.components.push_back(std::move(component));
    }
    evaluated.lines.reserve(payload.line_models.size());
    for (const auto& model : payload.line_models) {
        if (inactive_at_energy(model, query.e0_keV)) continue;
        EnvFitEvaluatedLine line;
        if (!evaluate_line_model(model, q, query.e0_keV, mu_edges, line, error))
            return false;
        line.incident_mu_bin_masses_per_cm2.resize(line.s_hat_mu_bin_masses.size());
        for (std::size_t i = 0; i < line.s_hat_mu_bin_masses.size(); ++i)
            line.incident_mu_bin_masses_per_cm2[i] =
                line.s_hat_mu_bin_masses[i] * query.direct_fluence_per_cm2;
        evaluated.lines.push_back(std::move(line));
    }
    result = std::move(evaluated);
    return true;
}

}  // namespace

std::string env_fit_payload_model_content_sha256(const EnvFitPayload& payload)
{
    return model_content_sha256(payload);
}

bool env_fit_transform_query(const EnvFitRuntimeQuery& query,
                             std::array<double, 4>& q,
                             std::string* error)
{
    if (!finite(query.e0_keV) || !finite(query.distance_m) ||
        !finite(query.source_height_m) || !finite(query.detector_height_m) ||
        query.e0_keV <= 0.0 || query.distance_m <= 0.0 ||
        query.source_height_m < 0.0 || query.detector_height_m < 0.0)
        return fail(error, "runtime query contains invalid physical coordinates");
    const double dz = query.source_height_m - query.detector_height_m;
    if (std::abs(dz) >= query.distance_m)
        return fail(error, "distance must exceed |source height - detector height|");
    const double rho2 = query.distance_m * query.distance_m - dz * dz;
    const double rho = std::sqrt(std::max(rho2,
        std::numeric_limits<double>::epsilon() * query.distance_m * query.distance_m));
    q = {{std::log(query.e0_keV), std::log(rho),
          std::atan2(query.source_height_m, rho),
          std::atan2(query.detector_height_m, rho)}};
    return true;
}

bool validate_env_fit_payload(const EnvFitPayload& payload, std::string* error)
{
    if (!valid_sha256(payload.config_sha256))
        return fail(error, "payload configuration SHA-256 is not lowercase hexadecimal");
    if (payload.material_name.empty() || payload.geometry_tier.empty())
        return fail(error, "payload material and geometry tier are required");
    if (payload.normalization != EnvFitPayloadNormalization::SHatToDirectFluence)
        return fail(error, "unsupported fitted-payload normalization");
    if (!finite(payload.brems_activation_threshold_keV) ||
        payload.brems_activation_threshold_keV < 0.0)
        return fail(error, "payload brems activation threshold is invalid");
    for (std::size_t i = 0; i < payload.physical_domain_min.size(); ++i) {
        if (!finite(payload.physical_domain_min[i]) ||
            !finite(payload.physical_domain_max[i]) ||
            payload.physical_domain_max[i] < payload.physical_domain_min[i])
            return fail(error, "payload physical domain bounds are invalid");
    }
    if (payload.physical_domain_min[0] <= 0.0 ||
        payload.physical_domain_min[1] <= 0.0 ||
        payload.physical_domain_min[2] < 0.0 ||
        payload.physical_domain_min[3] < 0.0)
        return fail(error, "payload physical domain has nonphysical lower bounds");
    const double height_gap = std::max({
        0.0,
        payload.physical_domain_min[2] - payload.physical_domain_max[3],
        payload.physical_domain_min[3] - payload.physical_domain_max[2]
    });
    if (height_gap >= payload.physical_domain_max[1])
        return fail(error, "payload physical domain contains no valid geometry");
    if (payload.validation_metadata_json.size() > kMaximumStringBytes)
        return fail(error, "validation metadata exceeds the format safety limit");
    if (payload.validation_metadata_json.size() < 2 ||
        payload.validation_metadata_json.front() != '{' ||
        payload.validation_metadata_json.back() != '}')
        return fail(error, "validation metadata must be a canonical JSON object");
    if (payload.continuum_models.empty() && payload.line_models.empty())
        return fail(error, "payload has no fitted models");
    if (payload.model_inventory.empty())
        return fail(error, "payload requires an explicit runtime model inventory");
    using InventoryKey = std::tuple<bool, std::string, std::uint32_t>;
    std::set<InventoryKey> expected_inventory;
    for (const auto& entry : payload.model_inventory) {
        if (!canonical_model_identity(entry.name, entry.axis_kind, entry.line))
            return fail(error, "model inventory contains an unknown name/axis entry");
        const InventoryKey key{
            entry.line, entry.name, static_cast<std::uint32_t>(entry.axis_kind)};
        if (!expected_inventory.insert(key).second)
            return fail(error, "model inventory contains a duplicate entry");
    }
    std::set<InventoryKey> actual_inventory;
    std::set<std::string> names;
    const EnvFitRuntimeModel* domain_reference = nullptr;
    const auto check_common_domain = [&](const EnvFitRuntimeModel& model) {
        if (!domain_reference) {
            domain_reference = &model;
            return true;
        }
        for (std::size_t i = 0; i < 4; ++i) {
            const auto close = [](double a, double b) {
                return std::abs(a - b) <= 1.0e-13 *
                    std::max({1.0, std::abs(a), std::abs(b)});
            };
            if (!close(model.q_min[i], domain_reference->q_min[i]) ||
                !close(model.q_span[i], domain_reference->q_span[i]))
                return false;
        }
        return true;
    };
    for (const auto& model : payload.continuum_models) {
        if (!validate_model(model, false,
                            payload.brems_activation_threshold_keV, error))
            return false;
        if (!check_common_domain(model))
            return fail(error, "all models must use the campaign-wide q domain");
        if (!names.insert(model.name).second)
            return fail(error, "duplicate fitted model name '" + model.name + "'");
        actual_inventory.emplace(false, model.name,
                                 static_cast<std::uint32_t>(model.axis_kind));
    }
    for (const auto& model : payload.line_models) {
        if (!validate_model(model, true,
                            payload.brems_activation_threshold_keV, error))
            return false;
        if (!check_common_domain(model))
            return fail(error, "all models must use the campaign-wide q domain");
        if (!names.insert(model.name).second)
            return fail(error, "duplicate fitted model name '" + model.name + "'");
        actual_inventory.emplace(true, model.name,
                                 static_cast<std::uint32_t>(model.axis_kind));
    }
    if (actual_inventory != expected_inventory)
        return fail(error, "explicit runtime inventory does not exactly match fitted models");
    const double q_energy_min = std::log(payload.physical_domain_min[0]);
    const double q_energy_max = std::log(payload.physical_domain_max[0]);
    if (domain_reference->q_min[0] > q_energy_min + 1.0e-12 ||
        domain_reference->q_min[0] + domain_reference->q_span[0] <
            q_energy_max - 1.0e-12)
        return fail(error, "campaign q domain does not cover the physical energy domain");
    if (!payload.model_content_sha256.empty()) {
        if (!valid_sha256(payload.model_content_sha256))
            return fail(error, "payload model-content SHA-256 is invalid");
        if (payload.model_content_sha256 != model_content_sha256(payload))
            return fail(error, "payload model-content SHA-256 mismatch");
    }
    return true;
}

bool evaluate_env_fit_payload(const EnvFitPayload& payload,
                              const EnvFitRuntimeQuery& query,
                              const std::vector<double>& energy_edges_keV,
                              const std::vector<double>& mu_edges,
                              EnvFitEvaluation& result,
                              std::string* error)
{
    if (!validate_env_fit_payload(payload, error)) return false;
    return evaluate_valid_payload(payload, query, energy_edges_keV, mu_edges,
                                  result, error);
}

bool save_env_fit_payload(const EnvFitPayload& payload,
                          const std::string& path,
                          std::string* error)
{
    if (path.empty()) return fail(error, "output path is empty");
    if (!validate_env_fit_payload(payload, error)) return false;
    Writer body;
    serialize_payload(payload, body);
    if (body.bytes().size() > kMaximumPayloadBytes)
        return fail(error, "fitted payload exceeds the format safety limit");
    const auto checksum = digest_bytes(body.bytes().data(), body.bytes().size());

    static std::atomic<std::uint64_t> counter{0};
#if defined(_WIN32)
    const unsigned long process_id = static_cast<unsigned long>(_getpid());
#else
    const unsigned long process_id = static_cast<unsigned long>(getpid());
#endif
    std::ostringstream suffix;
    suffix << ".tmp." << process_id << "."
           << counter.fetch_add(1, std::memory_order_relaxed);
    const std::string temporary = path + suffix.str();
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) return fail(error, "cannot open adjacent fitted-payload temporary file");

    Writer header;
    header.raw(kMagic.data(), kMagic.size());
    header.u32(kEnvFitPayloadFormatVersion);
    header.u32(kEndianTag);
    header.u64(static_cast<std::uint64_t>(body.bytes().size()));
    header.raw(checksum.data(), checksum.size());
    output.write(reinterpret_cast<const char*>(header.bytes().data()),
                 static_cast<std::streamsize>(header.bytes().size()));
    output.write(reinterpret_cast<const char*>(body.bytes().data()),
                 static_cast<std::streamsize>(body.bytes().size()));
    output.flush();
    output.close();
    if (!output) {
        std::remove(temporary.c_str());
        return fail(error, "failed while writing fitted-payload temporary file");
    }
    if (!durable_sync_file(temporary)) {
        std::remove(temporary.c_str());
        return fail(error, "cannot durably flush fitted-payload temporary file");
    }
    if (std::rename(temporary.c_str(), path.c_str()) != 0) {
        std::remove(temporary.c_str());
        return fail(error, "cannot atomically rename fitted-payload temporary file");
    }
    if (!durable_sync_parent_directory(path))
        return fail(error, "cannot durably flush fitted-payload output directory");
    return true;
}

bool load_env_fit_payload(EnvFitPayload& payload,
                          const std::string& path,
                          std::string* error,
                          bool allow_candidate)
{
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) return fail(error, "cannot open fitted-payload file");
    const std::streamoff end = input.tellg();
    if (end < static_cast<std::streamoff>(kHeaderBytes))
        return fail(error, "fitted-payload file is truncated");
    if (static_cast<std::uint64_t>(end) > kHeaderBytes + kMaximumPayloadBytes)
        return fail(error, "fitted-payload file exceeds the format safety limit");
    input.seekg(0);
    std::vector<std::uint8_t> bytes(static_cast<std::size_t>(end));
    input.read(reinterpret_cast<char*>(bytes.data()), end);
    if (!input) return fail(error, "cannot read complete fitted-payload file");

    try {
        Reader header(bytes.data(), static_cast<std::size_t>(kHeaderBytes));
        std::array<std::uint8_t, 16> magic{};
        std::array<std::uint8_t, 32> expected{};
        header.raw(magic.data(), magic.size());
        const std::uint32_t version = header.u32();
        const std::uint32_t endian = header.u32();
        const std::uint64_t body_size = header.u64();
        header.raw(expected.data(), expected.size());
        if (!header.done() || magic != kMagic)
            return fail(error, "not a CEEFIT-PAYLOAD-v1 file");
        if (version != kEnvFitPayloadFormatVersion || endian != kEndianTag)
            return fail(error, "unsupported fitted-payload version or byte order");
        if (body_size != bytes.size() - kHeaderBytes)
            return fail(error, "fitted-payload length does not match file size");
        const std::uint8_t* body = bytes.data() + kHeaderBytes;
        if (digest_bytes(body, static_cast<std::size_t>(body_size)) != expected)
            return fail(error, "fitted-payload SHA-256 checksum mismatch");
        EnvFitPayload decoded = deserialize_payload(body,
            static_cast<std::size_t>(body_size));
        if (!validate_env_fit_payload(decoded, error)) return false;
        if (!decoded.accepted_for_publication && !allow_candidate)
            return fail(error,
                "fitted payload is a candidate; validation tooling must opt in explicitly");
        payload = std::move(decoded);
        return true;
    } catch (const std::exception& exception) {
        return fail(error, std::string("cannot decode fitted payload: ") +
                           exception.what());
    }
}

bool EnvFitPayloadEvaluator::set_payload(EnvFitPayload payload,
                                         std::string* error,
                                         bool allow_candidate)
{
    if (!validate_env_fit_payload(payload, error)) return false;
    if (!payload.accepted_for_publication && !allow_candidate)
        return fail(error,
            "fitted payload is a candidate; validation tooling must opt in explicitly");
    payload_ = std::move(payload);
    ready_ = true;
    return true;
}

bool EnvFitPayloadEvaluator::load(const std::string& path, std::string* error,
                                  bool allow_candidate)
{
    EnvFitPayload decoded;
    if (!load_env_fit_payload(decoded, path, error, allow_candidate)) return false;
    payload_ = std::move(decoded);
    ready_ = true;
    return true;
}

bool EnvFitPayloadEvaluator::evaluate(
    const EnvFitRuntimeQuery& query,
    const std::vector<double>& energy_edges_keV,
    const std::vector<double>& mu_edges,
    EnvFitEvaluation& result,
    std::string* error) const
{
    if (!ready_) return fail(error, "fitted-payload evaluator is not initialized");
    return evaluate_valid_payload(payload_, query, energy_edges_keV, mu_edges,
                                  result, error);
}

}  // namespace ceelo
