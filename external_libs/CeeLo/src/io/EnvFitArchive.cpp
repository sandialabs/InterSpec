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

#include "io/EnvFitArchive.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>

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
    'C', 'E', 'E', 'F', 'I', 'T', '-', 'B',
    'L', 'O', 'C', 'K', '-', 'v', '1', '\0'
}};
constexpr std::uint64_t kFileHeaderBytes = 64;
constexpr std::uint64_t kMaximumPayloadBytes = 8ULL * 1024ULL * 1024ULL * 1024ULL;
constexpr std::uint64_t kMaximumStringBytes = 64ULL * 1024ULL * 1024ULL;
constexpr std::uint64_t kMaximumVectorElements = 1ULL << 30;

bool fail(std::string* error, const std::string& message)
{
    if (error) *error = message;
    return false;
}

bool durable_sync_file(const std::string& path)
{
#if defined(_WIN32)
    const int descriptor = _open(path.c_str(), _O_RDONLY | _O_BINARY);
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
    // ReplaceFile/MoveFileEx durability is platform-specific; the file itself
    // has already been committed above. There is no portable directory fsync.
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

// Minimal dependency-free SHA-256 (FIPS 180-4).  It is used both for the file
// payload checksum and for canonical campaign-configuration identities.
class Sha256 {
public:
    Sha256()
        : state_{{0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
                  0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U}}
    {
    }

    void update(const void* input, std::size_t size)
    {
        const auto* data = static_cast<const std::uint8_t*>(input);
        while (size != 0) {
            const std::size_t take = std::min(size, block_.size() - block_size_);
            std::memcpy(block_.data() + block_size_, data, take);
            block_size_ += take;
            total_bytes_ += static_cast<std::uint64_t>(take);
            data += take;
            size -= take;
            if (block_size_ == block_.size()) {
                transform(block_.data());
                block_size_ = 0;
            }
        }
    }

    std::array<std::uint8_t, 32> finish()
    {
        const std::uint64_t bit_count = total_bytes_ * 8ULL;
        block_[block_size_++] = 0x80U;
        if (block_size_ > 56) {
            std::fill(block_.begin() + static_cast<std::ptrdiff_t>(block_size_),
                      block_.end(), 0U);
            transform(block_.data());
            block_size_ = 0;
        }
        std::fill(block_.begin() + static_cast<std::ptrdiff_t>(block_size_),
                  block_.begin() + 56, 0U);
        for (std::size_t i = 0; i < 8; ++i)
            block_[63U - i] = static_cast<std::uint8_t>(bit_count >> (8U * i));
        transform(block_.data());

        std::array<std::uint8_t, 32> digest{};
        for (std::size_t i = 0; i < state_.size(); ++i) {
            digest[4 * i] = static_cast<std::uint8_t>(state_[i] >> 24);
            digest[4 * i + 1] = static_cast<std::uint8_t>(state_[i] >> 16);
            digest[4 * i + 2] = static_cast<std::uint8_t>(state_[i] >> 8);
            digest[4 * i + 3] = static_cast<std::uint8_t>(state_[i]);
        }
        return digest;
    }

private:
    static std::uint32_t rotr(std::uint32_t x, unsigned n)
    {
        return (x >> n) | (x << (32U - n));
    }

    void transform(const std::uint8_t* data)
    {
        static constexpr std::array<std::uint32_t, 64> k{{
            0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
            0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
            0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
            0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
            0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
            0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
            0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
            0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
            0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
            0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
            0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
            0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
            0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
            0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
            0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
            0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U
        }};
        std::array<std::uint32_t, 64> w{};
        for (std::size_t i = 0; i < 16; ++i) {
            w[i] = (static_cast<std::uint32_t>(data[4 * i]) << 24) |
                   (static_cast<std::uint32_t>(data[4 * i + 1]) << 16) |
                   (static_cast<std::uint32_t>(data[4 * i + 2]) << 8) |
                   static_cast<std::uint32_t>(data[4 * i + 3]);
        }
        for (std::size_t i = 16; i < w.size(); ++i) {
            const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^
                                     (w[i - 15] >> 3);
            const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^
                                     (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }

        std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (std::size_t i = 0; i < 64; ++i) {
            const std::uint32_t s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const std::uint32_t choose = (e & f) ^ (~e & g);
            const std::uint32_t t1 = h + s1 + choose + k[i] + w[i];
            const std::uint32_t s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t t2 = s0 + majority;
            h = g; g = f; f = e; e = d + t1;
            d = c; c = b; b = a; a = t1 + t2;
        }
        state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
        state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
    }

    std::array<std::uint32_t, 8> state_;
    std::array<std::uint8_t, 64> block_{};
    std::size_t block_size_ = 0;
    std::uint64_t total_bytes_ = 0;
};

std::array<std::uint8_t, 32> sha256(const void* data, std::size_t size)
{
    Sha256 hash;
    if (size != 0) hash.update(data, size);
    return hash.finish();
}

std::string hex_digest(const std::array<std::uint8_t, 32>& digest)
{
    static constexpr char hex[] = "0123456789abcdef";
    std::string result(64, '0');
    for (std::size_t i = 0; i < digest.size(); ++i) {
        result[2 * i] = hex[digest[i] >> 4];
        result[2 * i + 1] = hex[digest[i] & 0x0fU];
    }
    return result;
}

class Writer {
public:
    const std::vector<std::uint8_t>& bytes() const { return bytes_; }

    void u8(std::uint8_t value) { bytes_.push_back(value); }

    void u32(std::uint32_t value)
    {
        for (unsigned i = 0; i < 4; ++i)
            bytes_.push_back(static_cast<std::uint8_t>(value >> (8 * i)));
    }

    void u64(std::uint64_t value)
    {
        for (unsigned i = 0; i < 8; ++i)
            bytes_.push_back(static_cast<std::uint8_t>(value >> (8 * i)));
    }

    void f32(float value)
    {
        std::uint32_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "float must be binary32");
        std::memcpy(&bits, &value, sizeof(bits));
        u32(bits);
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

    template<class T, class Put>
    void vector(const std::vector<T>& values, Put put)
    {
        u64(static_cast<std::uint64_t>(values.size()));
        for (const T& value : values) put(value);
    }

private:
    std::vector<std::uint8_t> bytes_;
};

class Reader {
public:
    Reader(const std::uint8_t* data, std::size_t size)
        : data_(data), size_(size)
    {
    }

    std::size_t remaining() const { return size_ - position_; }
    bool done() const { return ok_ && position_ == size_; }
    bool ok() const { return ok_; }
    const std::string& error() const { return error_; }

    std::uint8_t u8()
    {
        if (!need(1)) return 0;
        return data_[position_++];
    }

    std::uint32_t u32()
    {
        if (!need(4)) return 0;
        std::uint32_t value = 0;
        for (unsigned i = 0; i < 4; ++i)
            value |= static_cast<std::uint32_t>(data_[position_++]) << (8 * i);
        return value;
    }

    std::uint64_t u64()
    {
        if (!need(8)) return 0;
        std::uint64_t value = 0;
        for (unsigned i = 0; i < 8; ++i)
            value |= static_cast<std::uint64_t>(data_[position_++]) << (8 * i);
        return value;
    }

    float f32()
    {
        const std::uint32_t bits = u32();
        float value = 0.0f;
        std::memcpy(&value, &bits, sizeof(value));
        return value;
    }

    double f64()
    {
        const std::uint64_t bits = u64();
        double value = 0.0;
        std::memcpy(&value, &bits, sizeof(value));
        return value;
    }

    bool raw(void* destination, std::size_t size)
    {
        if (!need(size)) return false;
        std::memcpy(destination, data_ + position_, size);
        position_ += size;
        return true;
    }

    std::string string()
    {
        const std::uint64_t count = u64();
        if (!ok_) return {};
        if (count > kMaximumStringBytes || count > remaining()) {
            set_error("invalid or excessive string length");
            return {};
        }
        std::string value(reinterpret_cast<const char*>(data_ + position_),
                          static_cast<std::size_t>(count));
        position_ += static_cast<std::size_t>(count);
        return value;
    }

    std::uint64_t count(std::size_t minimum_element_size = 1)
    {
        const std::uint64_t value = u64();
        if (!ok_) return 0;
        if (value > kMaximumVectorElements ||
            (minimum_element_size != 0 &&
             value > remaining() / minimum_element_size)) {
            set_error("invalid or excessive vector length");
            return 0;
        }
        return value;
    }

    void set_error(const std::string& message)
    {
        if (ok_) error_ = message;
        ok_ = false;
    }

private:
    bool need(std::size_t size)
    {
        if (!ok_ || size > remaining()) {
            set_error("truncated payload");
            return false;
        }
        return true;
    }

    const std::uint8_t* data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t position_ = 0;
    bool ok_ = true;
    std::string error_;
};

template<class T, class Put>
void put_vector(Writer& writer, const std::vector<T>& values, Put put)
{
    writer.u64(static_cast<std::uint64_t>(values.size()));
    for (const T& value : values) put(value);
}

template<class T, class Get>
std::vector<T> get_vector(Reader& reader, std::size_t minimum_element_size,
                          Get get)
{
    const std::uint64_t count = reader.count(minimum_element_size);
    if (!reader.ok()) return {};
    std::vector<T> values;
    values.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t i = 0; i < count && reader.ok(); ++i)
        values.push_back(get());
    return values;
}

void put_moments(Writer& writer, const EnvFitMomentArray& value)
{
    put_vector(writer, value.n, [&](std::uint64_t x) { writer.u64(x); });
    put_vector(writer, value.s1, [&](double x) { writer.f64(x); });
    put_vector(writer, value.s2, [&](double x) { writer.f64(x); });
}

EnvFitMomentArray get_moments(Reader& reader)
{
    EnvFitMomentArray value;
    value.n = get_vector<std::uint64_t>(reader, 8, [&] { return reader.u64(); });
    value.s1 = get_vector<double>(reader, 8, [&] { return reader.f64(); });
    value.s2 = get_vector<double>(reader, 8, [&] { return reader.f64(); });
    return value;
}

bool finite(double value) { return std::isfinite(value); }
bool finite(float value) { return std::isfinite(value); }

bool valid_sha256(const std::string& value)
{
    if (value.size() != 64) return false;
    return std::all_of(value.begin(), value.end(), [](char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    });
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
            if (c < 0x20)
                out << "\\u00" << hex[c >> 4] << hex[c & 0x0f];
            else
                out << static_cast<char>(c);
        }
    }
    out << '"';
    return out.str();
}

bool add_size(std::size_t& total, std::size_t add)
{
    if (add > std::numeric_limits<std::size_t>::max() - total) return false;
    total += add;
    return true;
}

bool multiply_size(std::size_t a, std::size_t b, std::size_t& result)
{
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) return false;
    result = a * b;
    return true;
}

void serialize_payload(const EnvFitBlock& block, Writer& writer);
bool deserialize_payload(const std::uint8_t* data, std::size_t size,
                         EnvFitBlock& block, std::string* error);
bool ascending_finite(const std::vector<double>& values);
bool validate_moments(const EnvFitMomentArray& moments, std::size_t expected,
                      std::uint64_t histories, const std::string& name,
                      std::string* error);
bool validate_replica_sums(const std::vector<float>& values,
                           std::size_t cells, const std::string& name,
                           std::string* error);
bool validate_replica_closure(const std::vector<float>& values,
                              std::size_t cells,
                              const EnvFitMomentArray& moments,
                              const std::string& name,
                              std::string* error);
bool validate_line_replica_additivity(const EnvFitBlock& block,
                                      std::string* error);
bool same_manifest(const EnvFitManifest& a, const EnvFitManifest& b);
bool same_axes(const EnvFitAxes& a, const EnvFitAxes& b);
bool same_components(const std::vector<EnvFitComponent>& a,
                     const std::vector<EnvFitComponent>& b);
bool same_node_identity(const EnvFitNodeMetadata& a,
                        const EnvFitNodeMetadata& b);
bool same_diagnostic_layout(const EnvFitDiagnostics& a,
                            const EnvFitDiagnostics& b);
bool ranges_overlap(const EnvFitReplicaRecord& a,
                    const EnvFitReplicaRecord& b);
void add_moments(EnvFitMomentArray& destination,
                 const EnvFitMomentArray& source);
template<class T>
void add_vector(std::vector<T>& destination, const std::vector<T>& source);

} // namespace

std::string env_fit_stopping_definition_json(const EnvFitAxes& axes)
{
    if (axes.total_energy_axis >= axes.energy_axes.size()) return {};
    const EnvFitEnergyAxis& axis = axes.energy_axes[axes.total_energy_axis];
    if (axis.name.empty() || axis.physical_edges_keV.size() < 2) return {};
    const std::size_t energy_bins = axis.physical_edges_keV.size() - 1;
    constexpr std::size_t bands = 12;
    std::ostringstream out;
    out << "{\"angular_basis\":\"standard Legendre P_l(mu), "
           "mu=cos(direction_of_travel,+detector_axis)\",\"dimension\":16,"
           "\"entries\":[";
    for (std::size_t band = 0; band < bands; ++band) {
        if (band != 0) out << ',';
        const std::size_t begin = (band * energy_bins + bands - 1) / bands;
        const std::size_t end = ((band + 1) * energy_bins + bands - 1) / bands;
        out << "{\"energy_axis_index\":" << axes.total_energy_axis
            << ",\"energy_axis_name\":" << json_string(axis.name)
            << ",\"energy_bin_begin\":" << begin
            << ",\"energy_bin_end\":" << end
            << ",\"kind\":\"registered_energy_index_band\",\"name\":\"registered_band_"
            << std::setw(2) << std::setfill('0') << band << "\"}";
    }
    for (std::size_t order = 1; order <= 4; ++order) {
        out << ",{\"kind\":\"legendre_angular_moment\",\"legendre_order\":"
            << order << ",\"name\":\"legendre_l" << order << "\"}";
    }
    out << "],\"normalization\":\"raw per-history scattered continuum-plus-lines "
           "fluence moments\",\"version\":\"CEEFIT-STOPPING-v1\"}";
    return out.str();
}

bool validate_env_fit_manifest(const EnvFitManifest& manifest,
                               std::string* error)
{
    if (error) error->clear();
    if (manifest.canonical_config_json.empty())
        return fail(error, "canonical configuration JSON is empty");
    if (!valid_sha256(manifest.config_sha256) ||
        manifest.config_sha256 != env_fit_sha256_hex(manifest.canonical_config_json))
        return fail(error, "config_sha256 does not match canonical configuration");
    if (manifest.material_name.empty() ||
        manifest.material_definition_json.empty() || manifest.geometry_tier.empty() ||
        manifest.rng_algorithm.empty() || manifest.estimator_flags.empty() ||
        manifest.units.empty() || manifest.observable.empty() ||
        manifest.normalization.empty())
        return fail(error,
                    "manifest is missing a required material, estimator, or "
                    "observable field");
    return true;
}

bool validate_env_fit_block(const EnvFitBlock& block, std::string* error)
{
    if (error) error->clear();
    if (block.node_id == 0 || block.block_id == 0)
        return fail(error, "node_id and block_id must be nonzero");

    const EnvFitManifest& manifest = block.manifest;
    if (!validate_env_fit_manifest(manifest, error)) return false;

    if (block.axes.energy_axes.empty() ||
        block.axes.total_energy_axis >= block.axes.energy_axes.size())
        return fail(error, "missing or invalid energy axes");
    std::set<std::string> axis_names;
    for (const EnvFitEnergyAxis& axis : block.axes.energy_axes) {
        if (axis.name.empty() || !axis_names.insert(axis.name).second ||
            axis.transform_version == 0 ||
            axis.registered_edges.size() != axis.physical_edges_keV.size() ||
            !ascending_finite(axis.registered_edges) ||
            !ascending_finite(axis.physical_edges_keV) ||
            axis.physical_edges_keV.front() < 0.0)
            return fail(error, "invalid or duplicate energy axis");
    }
    if (!ascending_finite(block.axes.mu_edges) ||
        block.axes.mu_edges.front() < -1.0 || block.axes.mu_edges.back() > 1.0)
        return fail(error, "invalid incidence-cosine axis");

    std::set<std::string> line_names;
    for (const EnvFitLineDefinition& line : block.axes.lines) {
        if (line.name.empty() || !line_names.insert(line.name).second ||
            !finite(line.energy_keV) || line.energy_keV <= 0.0)
            return fail(error, "invalid or duplicate line definition");
    }

    if (block.components.empty())
        return fail(error, "worker block contains no continuum components");
    std::set<std::string> component_names;
    std::set<std::pair<std::uint32_t, std::uint32_t>> component_tags;
    for (const EnvFitComponent& component : block.components) {
        if (component.name.empty() || !component_names.insert(component.name).second ||
            !component_tags.insert({component.origin,
                                    component.post_creation_order}).second ||
            component.energy_axis >= block.axes.energy_axes.size())
            return fail(error, "invalid or duplicate continuum component");
        const std::size_t bins =
            block.axes.energy_axes[component.energy_axis].physical_edges_keV.size() - 1;
        if (component.support_mask.size() != bins ||
            std::any_of(component.support_mask.begin(), component.support_mask.end(),
                        [](std::uint8_t value) { return value > 1U; }))
            return fail(error, "component support mask has the wrong dimensions");
    }

    const EnvFitNodeMetadata& node = block.node;
    const double node_scalars[] = {
        node.e0_keV, node.distance_cm, node.h_src_cm, node.h_det_cm,
        node.transverse_extent_cm, node.direct_transmission,
        node.direct_fluence, node.wall_seconds, node.cpu_seconds,
        node.target_total_rel_se, node.achieved_total_rel_se
    };
    if (std::any_of(std::begin(node_scalars), std::end(node_scalars),
                    [](double value) { return !finite(value); }) ||
        std::any_of(node.transformed_coordinates.begin(),
                    node.transformed_coordinates.end(),
                    [](double value) { return !finite(value); }))
        return fail(error, "node metadata contains a non-finite value");
    if (node.e0_keV <= 0.0 || node.distance_cm <= 0.0 || node.h_src_cm < 0.0 ||
        node.h_det_cm < 0.0 || node.transverse_extent_cm < 0.0 ||
        node.direct_transmission < 0.0 || node.direct_fluence < 0.0 ||
        node.wall_seconds < 0.0 || node.cpu_seconds < 0.0 ||
        node.target_total_rel_se < 0.0 || node.achieved_total_rel_se < 0.0 ||
        node.histories == 0 || node.role.empty() || node.status.empty())
        return fail(error, "node metadata is outside its physical domain");
    const double vertical = node.h_src_cm - node.h_det_cm;
    if (!(node.distance_cm > std::abs(vertical)))
        return fail(error, "node distance is shorter than its vertical separation");
    const double distance_m = 0.01 * node.distance_cm;
    const double h_src_m = 0.01 * node.h_src_cm;
    const double h_det_m = 0.01 * node.h_det_cm;
    const double dz_m = h_src_m - h_det_m;
    const double rho_m = std::sqrt(distance_m * distance_m - dz_m * dz_m);
    const std::array<double, 4> expected_coordinates{{
        std::log(node.e0_keV), std::log(rho_m),
        std::atan2(h_src_m, rho_m), std::atan2(h_det_m, rho_m)
    }};
    for (std::size_t i = 0; i < expected_coordinates.size(); ++i) {
        const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(expected_coordinates[i]));
        if (std::abs(node.transformed_coordinates[i] - expected_coordinates[i]) >
            tolerance)
            return fail(error, "stored transformed coordinates do not match node");
    }
    if (block.node_id != env_fit_node_id(manifest.config_sha256, node))
        return fail(error, "node_id does not match configuration and coordinates");

    std::uint64_t history_sum = 0;
    std::uint64_t expected_replica_count = 0;
    std::set<std::uint64_t> seeds;
    for (std::size_t i = 0; i < block.replicas.size(); ++i) {
        const EnvFitReplicaRecord& replica = block.replicas[i];
        if (replica.replica_id != i || replica.history_count == 0 ||
            replica.history_begin >
                std::numeric_limits<std::uint64_t>::max() - replica.history_count ||
            !seeds.insert(replica.rng_seed).second)
            return fail(error, "invalid replica id, range, count, or duplicate seed");
        if (i == 0) expected_replica_count = replica.history_count;
        if (replica.history_count != expected_replica_count)
            return fail(error, "replica history counts are not equal");
        if (history_sum > std::numeric_limits<std::uint64_t>::max() -
                              replica.history_count)
            return fail(error, "block history count overflows uint64");
        history_sum += replica.history_count;
        for (std::size_t j = 0; j < i; ++j) {
            if (ranges_overlap(replica, block.replicas[j]))
                return fail(error, "replica history ranges overlap");
        }
    }
    if (history_sum != node.histories)
        return fail(error, "replica counts do not sum to node.histories");

    const std::size_t component_cells =
        env_fit_continuum_component_cells(block);
    const std::size_t total_cells = env_fit_continuum_total_cells(block);
    if (component_cells == 0 || total_cells == 0)
        return fail(error, "continuum tensor dimensions overflow or are empty");
    const std::size_t line_component_cells =
        env_fit_line_component_cells(block);
    const std::size_t line_total_cells = env_fit_line_total_cells(block);
    if (!validate_moments(block.continuum_components, component_cells,
                          node.histories, "component continuum", error) ||
        !validate_moments(block.continuum_total, total_cells, node.histories,
                          "total continuum", error) ||
        !validate_moments(block.line_components, line_component_cells,
                          node.histories, "component line", error) ||
        !validate_moments(block.line_total, line_total_cells,
                          node.histories, "total line", error) ||
        !validate_replica_sums(block.replica_sums.continuum_components,
                               component_cells, "component continuum", error) ||
        !validate_replica_sums(block.replica_sums.continuum_total,
                               total_cells, "total continuum", error) ||
        !validate_replica_sums(block.replica_sums.line_components,
                               line_component_cells, "component line", error) ||
        !validate_replica_sums(block.replica_sums.line_total,
                               line_total_cells, "total line", error))
        return false;
    if (!validate_replica_closure(
            block.replica_sums.continuum_components, component_cells,
            block.continuum_components, "component continuum", error) ||
        !validate_replica_closure(
            block.replica_sums.continuum_total, total_cells,
            block.continuum_total, "total continuum", error) ||
        !validate_replica_closure(
            block.replica_sums.line_components, line_component_cells,
            block.line_components, "component line", error) ||
        !validate_replica_closure(
            block.replica_sums.line_total, line_total_cells,
            block.line_total, "total line", error) ||
        !validate_line_replica_additivity(block, error))
        return false;

    const std::string expected_stopping_definition =
        env_fit_stopping_definition_json(block.axes);
    if (expected_stopping_definition.empty() ||
        block.stopping.definition_json != expected_stopping_definition)
        return fail(error,
                    "stopping-vector definition is absent or noncanonical");
    if (block.stopping.n != node.histories)
        return fail(error, "stopping-moment N does not equal block histories");
    for (double value : block.stopping.s1)
        if (!finite(value)) return fail(error, "stopping S1 is non-finite");
    for (std::size_t i = 0; i < kEnvFitStoppingDimension; ++i) {
        const double diagonal =
            block.stopping.outer_s2[i * kEnvFitStoppingDimension + i];
        if (!finite(diagonal) || diagonal < 0.0)
            return fail(error, "stopping outer moment has an invalid diagonal");
        const long double minimum =
            static_cast<long double>(block.stopping.s1[i]) *
            block.stopping.s1[i] / static_cast<long double>(node.histories);
        const long double tolerance = 128.0L *
            std::numeric_limits<double>::epsilon() *
            std::max(1.0L, std::fabs(minimum));
        if (static_cast<long double>(diagonal) + tolerance < minimum)
            return fail(error, "stopping diagonal has S2 < S1^2/N");
        for (std::size_t j = 0; j < kEnvFitStoppingDimension; ++j) {
            const double a =
                block.stopping.outer_s2[i * kEnvFitStoppingDimension + j];
            const double b =
                block.stopping.outer_s2[j * kEnvFitStoppingDimension + i];
            if (!finite(a) ||
                std::abs(a - b) > 64.0 * std::numeric_limits<double>::epsilon() *
                                      std::max({1.0, std::abs(a), std::abs(b)}))
                return fail(error, "stopping outer moment is non-finite or asymmetric");
        }
    }

    const EnvFitDiagnostics& diagnostics = block.diagnostics;
    std::set<std::string> diagnostic_names;
    for (const EnvFitDiagnosticChannel& channel : diagnostics.path_channels) {
        if (channel.name.empty() || !diagnostic_names.insert(channel.name).second)
            return fail(error, "invalid or duplicate diagnostic channel");
        if (!validate_moments(channel.continuum, total_cells, node.histories,
                              "diagnostic continuum", error) ||
            !validate_moments(channel.lines, line_total_cells, node.histories,
                              "diagnostic lines", error))
            return false;
    }
    if (diagnostics.interaction_names.size() !=
        diagnostics.interaction_counts.size())
        return fail(error, "interaction name/count dimensions differ");
    std::set<std::string> interaction_names;
    for (const std::string& name : diagnostics.interaction_names)
        if (name.empty() || !interaction_names.insert(name).second)
            return fail(error, "invalid or duplicate interaction counter");
    if (diagnostics.scattered_total.n != node.histories ||
        !finite(diagnostics.scattered_total.s1) ||
        !finite(diagnostics.scattered_total.s2) ||
        diagnostics.scattered_total.s1 < 0.0 ||
        diagnostics.scattered_total.s2 < 0.0 ||
        diagnostics.dropped_weight.n != node.histories ||
        !finite(diagnostics.dropped_weight.s1) ||
        !finite(diagnostics.dropped_weight.s2) ||
        diagnostics.dropped_weight.s1 < 0.0 ||
        diagnostics.dropped_weight.s2 < 0.0 ||
        !finite(diagnostics.maximum_history_contribution) ||
        diagnostics.maximum_history_contribution < 0.0)
        return fail(error, "invalid total, dropped-weight, or maximum-history diagnostic");
    const long double total_minimum =
        static_cast<long double>(diagnostics.scattered_total.s1) *
        diagnostics.scattered_total.s1 / static_cast<long double>(node.histories);
    if (static_cast<long double>(diagnostics.scattered_total.s2) +
            128.0L * std::numeric_limits<double>::epsilon() *
                std::max(1.0L, std::fabs(total_minimum)) <
        total_minimum)
        return fail(error, "scattered-total S2 is less than S1^2/N");
    double band_total_s1 = 0.0;
    double band_total_s2 = 0.0;
    for (std::size_t i = 0; i < 12; ++i) {
        band_total_s1 += block.stopping.s1[i];
        for (std::size_t j = 0; j < 12; ++j)
            band_total_s2 += block.stopping.outer_s2[
                i * kEnvFitStoppingDimension + j];
    }
    const auto closes = [](double observed, double expected) {
        return std::abs(observed - expected) <=
            256.0 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::abs(observed), std::abs(expected)});
    };
    if (!closes(band_total_s1, diagnostics.scattered_total.s1) ||
        !closes(band_total_s2, diagnostics.scattered_total.s2))
        return fail(error,
                    "registered energy-band moments do not close to scattered total");
    const long double dropped_minimum =
        static_cast<long double>(diagnostics.dropped_weight.s1) *
        diagnostics.dropped_weight.s1 / static_cast<long double>(node.histories);
    if (static_cast<long double>(diagnostics.dropped_weight.s2) +
            128.0L * std::numeric_limits<double>::epsilon() *
                std::max(1.0L, std::fabs(dropped_minimum)) <
        dropped_minimum)
        return fail(error, "dropped-weight S2 is less than S1^2/N");

    return true;
}

bool save_env_fit_block_atomic(const EnvFitBlock& block,
                               const std::string& path,
                               std::string* error)
{
    if (path.empty()) return fail(error, "output path is empty");
    if (!validate_env_fit_block(block, error)) return false;
    if (block.node.histories != kEnvFitFormalBlockHistories)
        return fail(error,
                    "worker blocks must contain exactly 262144 histories");

    Writer payload;
    try {
        serialize_payload(block, payload);
    } catch (const std::exception& exception) {
        return fail(error, std::string("cannot serialize worker block: ") +
                           exception.what());
    }
    if (payload.bytes().size() > kMaximumPayloadBytes)
        return fail(error, "worker-block payload exceeds the format safety limit");
    const auto checksum = sha256(payload.bytes().data(), payload.bytes().size());

    static std::atomic<std::uint64_t> counter{0};
#if defined(_WIN32)
    const unsigned long process_id = static_cast<unsigned long>(_getpid());
#else
    const unsigned long process_id = static_cast<unsigned long>(getpid());
#endif
    const std::uint64_t serial = counter.fetch_add(1, std::memory_order_relaxed);
    std::ostringstream suffix;
    suffix << ".tmp." << process_id << "." << serial;
    const std::string temporary = path + suffix.str();

    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) return fail(error, "cannot open adjacent temporary output file");

    Writer header;
    header.raw(kMagic.data(), kMagic.size());
    header.u32(kEnvFitBlockFormatVersion);
    header.u32(0U);  // reserved flags
    header.u64(static_cast<std::uint64_t>(payload.bytes().size()));
    header.raw(checksum.data(), checksum.size());
    output.write(reinterpret_cast<const char*>(header.bytes().data()),
                 static_cast<std::streamsize>(header.bytes().size()));
    output.write(reinterpret_cast<const char*>(payload.bytes().data()),
                 static_cast<std::streamsize>(payload.bytes().size()));
    output.flush();
    output.close();
    if (!output) {
        std::remove(temporary.c_str());
        return fail(error, "failed while writing worker-block temporary file");
    }
    if (!durable_sync_file(temporary)) {
        std::remove(temporary.c_str());
        return fail(error, "cannot durably flush worker-block temporary file");
    }
    if (std::rename(temporary.c_str(), path.c_str()) != 0) {
        std::remove(temporary.c_str());
        return fail(error, "cannot atomically rename worker-block temporary file");
    }
    if (!durable_sync_parent_directory(path))
        return fail(error, "cannot durably flush worker-block output directory");
    return true;
}

bool load_env_fit_block(EnvFitBlock& out, const std::string& path,
                        std::string* error)
{
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) return fail(error, "cannot open worker-block file");
    const std::streamoff end = input.tellg();
    if (end < static_cast<std::streamoff>(kFileHeaderBytes))
        return fail(error, "worker-block file is truncated");
    if (static_cast<std::uint64_t>(end) >
        kFileHeaderBytes + kMaximumPayloadBytes)
        return fail(error, "worker-block file exceeds the format safety limit");
    input.seekg(0);
    std::vector<std::uint8_t> bytes(static_cast<std::size_t>(end));
    input.read(reinterpret_cast<char*>(bytes.data()), end);
    if (!input) return fail(error, "cannot read complete worker-block file");

    Reader header(bytes.data(), static_cast<std::size_t>(kFileHeaderBytes));
    std::array<std::uint8_t, 16> magic{};
    std::array<std::uint8_t, 32> stored_checksum{};
    header.raw(magic.data(), magic.size());
    const std::uint32_t version = header.u32();
    const std::uint32_t flags = header.u32();
    const std::uint64_t payload_size = header.u64();
    header.raw(stored_checksum.data(), stored_checksum.size());
    if (!header.done()) return fail(error, "malformed worker-block header");
    if (magic != kMagic) return fail(error, "not a CEEFIT-BLOCK-v1 file");
    if (version != kEnvFitBlockFormatVersion || flags != 0U)
        return fail(error, "unsupported worker-block version or flags");
    if (payload_size != bytes.size() - kFileHeaderBytes)
        return fail(error, "worker-block payload length does not match file size");

    const std::uint8_t* payload = bytes.data() + kFileHeaderBytes;
    const auto actual_checksum = sha256(payload, static_cast<std::size_t>(payload_size));
    if (actual_checksum != stored_checksum)
        return fail(error, "worker-block SHA-256 checksum mismatch");

    EnvFitBlock decoded;
    if (!deserialize_payload(payload, static_cast<std::size_t>(payload_size),
                             decoded, error) ||
        !validate_env_fit_block(decoded, error))
        return false;
    out = std::move(decoded);
    return true;
}

EnvFitEstimate env_fit_estimate(std::uint64_t n, double s1, double s2)
{
    EnvFitEstimate result;
    if (n == 0) {
        result.mean = std::numeric_limits<double>::quiet_NaN();
        result.standard_error = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    result.mean = s1 / static_cast<double>(n);
    if (n < 2 || !finite(s1) || !finite(s2)) {
        result.standard_error = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    const long double numerator =
        static_cast<long double>(s2) -
        static_cast<long double>(s1) * s1 / static_cast<long double>(n);
    const long double scale = std::max(
        1.0L, std::max(std::fabs(static_cast<long double>(s2)),
                       std::fabs(static_cast<long double>(s1) * s1 /
                                 static_cast<long double>(n))));
    const long double tolerance = 128.0L *
        std::numeric_limits<double>::epsilon() * scale;
    if (numerator < -tolerance) {
        result.standard_error = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    const long double variance_of_mean =
        std::max(0.0L, numerator) /
        (static_cast<long double>(n) * static_cast<long double>(n - 1));
    result.standard_error = std::sqrt(static_cast<double>(variance_of_mean));
    return result;
}

EnvFitEstimate env_fit_estimate(const EnvFitMomentArray& moments,
                                std::size_t index)
{
    if (index >= moments.n.size() || index >= moments.s1.size() ||
        index >= moments.s2.size())
        return env_fit_estimate(0, 0.0, 0.0);
    return env_fit_estimate(moments.n[index], moments.s1[index],
                            moments.s2[index]);
}

bool EnvFitBlockReducer::add(const EnvFitBlock& block, std::string* error)
{
    if (!validate_env_fit_block(block, error)) return false;

    if (empty()) {
        aggregate_ = block;
        aggregate_.block_id = 0;  // Explicitly a synthetic, non-serializable sum.
        for (EnvFitReplicaRecord& replica : aggregate_.replicas) {
            replica.history_begin = 0;
            replica.rng_seed = 0;
        }
        block_ids_.push_back(block.block_id);
        history_ranges_.reserve(kEnvFitReplicaCount);
        for (const EnvFitReplicaRecord& replica : block.replicas)
            history_ranges_.push_back({block.block_id, replica});
        return true;
    }

    if (std::find(block_ids_.begin(), block_ids_.end(), block.block_id) !=
        block_ids_.end())
        return fail(error, "duplicate block_id");
    if (!same_manifest(block.manifest, aggregate_.manifest))
        return fail(error, "block manifest/configuration mismatch");
    if (block.node_id != aggregate_.node_id)
        return fail(error, "block belongs to a different node_id");
    if (!same_axes(block.axes, aggregate_.axes) ||
        !same_components(block.components, aggregate_.components))
        return fail(error, "block axes/component layout mismatch");
    if (!same_node_identity(block.node, aggregate_.node))
        return fail(error, "block node coordinates or deterministic values mismatch");
    if (!same_diagnostic_layout(block.diagnostics, aggregate_.diagnostics))
        return fail(error, "block diagnostic layout mismatch");
    for (const EnvFitReplicaRecord& incoming : block.replicas) {
        for (const EnvFitSourceRange& existing : history_ranges_) {
            if (ranges_overlap(incoming, existing.replica))
                return fail(error, "history range overlaps an already merged block");
            if (incoming.rng_seed == existing.replica.rng_seed)
                return fail(error, "RNG seed is reused by an already merged block");
        }
    }
    if (aggregate_.node.histories >
        std::numeric_limits<std::uint64_t>::max() - block.node.histories)
        return fail(error, "merged history count overflows uint64");
    for (std::size_t i = 0; i < kEnvFitReplicaCount; ++i) {
        if (aggregate_.replicas[i].history_count >
            std::numeric_limits<std::uint64_t>::max() -
                block.replicas[i].history_count)
            return fail(error, "merged replica history count overflows uint64");
    }
    for (std::size_t i = 0; i < aggregate_.diagnostics.interaction_counts.size(); ++i) {
        if (aggregate_.diagnostics.interaction_counts[i] >
            std::numeric_limits<std::uint64_t>::max() -
                block.diagnostics.interaction_counts[i])
            return fail(error, "merged interaction count overflows uint64");
    }

    add_moments(aggregate_.continuum_components, block.continuum_components);
    add_moments(aggregate_.continuum_total, block.continuum_total);
    add_moments(aggregate_.line_components, block.line_components);
    add_moments(aggregate_.line_total, block.line_total);
    add_vector(aggregate_.replica_sums.continuum_components,
               block.replica_sums.continuum_components);
    add_vector(aggregate_.replica_sums.continuum_total,
               block.replica_sums.continuum_total);
    add_vector(aggregate_.replica_sums.line_components,
               block.replica_sums.line_components);
    add_vector(aggregate_.replica_sums.line_total,
               block.replica_sums.line_total);

    aggregate_.stopping.n += block.stopping.n;
    for (std::size_t i = 0; i < aggregate_.stopping.s1.size(); ++i)
        aggregate_.stopping.s1[i] += block.stopping.s1[i];
    for (std::size_t i = 0; i < aggregate_.stopping.outer_s2.size(); ++i)
        aggregate_.stopping.outer_s2[i] += block.stopping.outer_s2[i];

    for (std::size_t i = 0; i < aggregate_.diagnostics.path_channels.size(); ++i) {
        add_moments(aggregate_.diagnostics.path_channels[i].continuum,
                    block.diagnostics.path_channels[i].continuum);
        add_moments(aggregate_.diagnostics.path_channels[i].lines,
                    block.diagnostics.path_channels[i].lines);
    }
    add_vector(aggregate_.diagnostics.interaction_counts,
               block.diagnostics.interaction_counts);
    aggregate_.diagnostics.scattered_total.n +=
        block.diagnostics.scattered_total.n;
    aggregate_.diagnostics.scattered_total.s1 +=
        block.diagnostics.scattered_total.s1;
    aggregate_.diagnostics.scattered_total.s2 +=
        block.diagnostics.scattered_total.s2;
    aggregate_.diagnostics.dropped_weight.n +=
        block.diagnostics.dropped_weight.n;
    aggregate_.diagnostics.dropped_weight.s1 +=
        block.diagnostics.dropped_weight.s1;
    aggregate_.diagnostics.dropped_weight.s2 +=
        block.diagnostics.dropped_weight.s2;
    aggregate_.diagnostics.transport_failures +=
        block.diagnostics.transport_failures;
    aggregate_.diagnostics.maximum_history_contribution = std::max(
        aggregate_.diagnostics.maximum_history_contribution,
        block.diagnostics.maximum_history_contribution);
    aggregate_.diagnostics.capped =
        aggregate_.diagnostics.capped || block.diagnostics.capped;
    if (!block.diagnostics.failure_reason.empty()) {
        if (!aggregate_.diagnostics.failure_reason.empty())
            aggregate_.diagnostics.failure_reason += "; ";
        aggregate_.diagnostics.failure_reason += block.diagnostics.failure_reason;
    }

    aggregate_.node.histories += block.node.histories;
    aggregate_.node.wall_seconds += block.node.wall_seconds;
    aggregate_.node.cpu_seconds += block.node.cpu_seconds;
    aggregate_.node.status = aggregate_.diagnostics.transport_failures != 0
                                 ? "failed"
                                 : (aggregate_.diagnostics.capped ? "capped" : "merged");
    // The first twelve stopping entries are the registered-energy bands. Their
    // exact outer moment gives the uncertainty of their summed total.
    double total_s1 = 0.0;
    double total_s2 = 0.0;
    for (std::size_t i = 0; i < 12; ++i) {
        total_s1 += aggregate_.stopping.s1[i];
        for (std::size_t j = 0; j < 12; ++j)
            total_s2 += aggregate_.stopping.outer_s2[
                i * kEnvFitStoppingDimension + j];
    }
    const EnvFitEstimate total =
        env_fit_estimate(aggregate_.stopping.n, total_s1, total_s2);
    aggregate_.node.achieved_total_rel_se =
        (total.mean > 0.0 && finite(total.standard_error))
            ? total.standard_error / total.mean
            : 0.0;
    for (std::size_t i = 0; i < kEnvFitReplicaCount; ++i)
        aggregate_.replicas[i].history_count += block.replicas[i].history_count;

    block_ids_.push_back(block.block_id);
    for (const EnvFitReplicaRecord& replica : block.replicas)
        history_ranges_.push_back({block.block_id, replica});
    return true;
}

std::string env_fit_sha256_hex(const void* data, std::size_t size)
{
    if (size != 0 && data == nullptr) return {};
    return hex_digest(sha256(data, size));
}

std::string env_fit_sha256_hex(const std::string& text)
{
    return env_fit_sha256_hex(text.data(), text.size());
}

std::uint64_t env_fit_stable_id(const std::string& canonical_key)
{
    const auto digest = sha256(canonical_key.data(), canonical_key.size());
    std::uint64_t value = 0;
    // Interpret the leading digest bytes in network order so the printed id is
    // stable independently of host endianness.
    for (std::size_t i = 0; i < 8; ++i)
        value = (value << 8) | digest[i];
    return value;
}

std::uint64_t env_fit_node_id(const std::string& config_sha256,
                              const EnvFitNodeMetadata& node)
{
    static constexpr char domain[] = "CEEFIT-NODE-v1";
    Writer identity;
    identity.raw(domain, sizeof(domain)); // Include the terminating domain NUL.
    identity.raw(config_sha256.data(), config_sha256.size());
    identity.f64(node.e0_keV);
    identity.f64(node.distance_cm);
    identity.f64(node.h_src_cm);
    identity.f64(node.h_det_cm);
    identity.f64(node.transverse_extent_cm);
    const auto digest = sha256(identity.bytes().data(), identity.bytes().size());
    std::uint64_t value = 0;
    for (std::size_t i = 0; i < 8; ++i)
        value = (value << 8) | digest[i];
    return value;
}

std::size_t env_fit_continuum_component_cells(const EnvFitBlock& block)
{
    if (block.axes.mu_edges.size() < 2) return 0;
    const std::size_t n_mu = block.axes.mu_edges.size() - 1;
    std::size_t total = 0;
    for (const EnvFitComponent& component : block.components) {
        if (component.energy_axis >= block.axes.energy_axes.size()) return 0;
        const auto& edges = block.axes.energy_axes[component.energy_axis].physical_edges_keV;
        if (edges.size() < 2) return 0;
        std::size_t cells = 0;
        if (!multiply_size(edges.size() - 1, n_mu, cells) ||
            !add_size(total, cells))
            return 0;
    }
    return total;
}

std::size_t env_fit_continuum_total_cells(const EnvFitBlock& block)
{
    if (block.axes.mu_edges.size() < 2 ||
        block.axes.total_energy_axis >= block.axes.energy_axes.size())
        return 0;
    const auto& edges =
        block.axes.energy_axes[block.axes.total_energy_axis].physical_edges_keV;
    if (edges.size() < 2) return 0;
    std::size_t result = 0;
    return multiply_size(edges.size() - 1, block.axes.mu_edges.size() - 1,
                         result) ? result : 0;
}

std::size_t env_fit_line_component_cells(const EnvFitBlock& block)
{
    if (block.axes.mu_edges.size() < 2) return 0;
    std::size_t component_lines = 0;
    if (!multiply_size(block.components.size(), block.axes.lines.size(),
                       component_lines))
        return 0;
    std::size_t result = 0;
    return multiply_size(component_lines, block.axes.mu_edges.size() - 1,
                         result) ? result : 0;
}

std::size_t env_fit_line_total_cells(const EnvFitBlock& block)
{
    if (block.axes.mu_edges.size() < 2) return 0;
    std::size_t result = 0;
    return multiply_size(block.axes.lines.size(), block.axes.mu_edges.size() - 1,
                         result) ? result : 0;
}

namespace {

// Payload layout (all integers and IEEE-754 values little-endian; every string
// and variable-length vector starts with a uint64 byte/element count):
// ids, manifest, axes, components, node metadata, replica count+records,
// three moment tensors, three replica-S1 tensors, stopping moments, diagnostics.
// The 64-byte outer header is magic[16], version u32, flags u32, payload_size
// u64, and SHA-256(payload)[32].
void serialize_payload(const EnvFitBlock& block, Writer& writer)
{
    writer.u64(block.node_id);
    writer.u64(block.block_id);

    const EnvFitManifest& manifest = block.manifest;
    writer.string(manifest.canonical_config_json);
    writer.string(manifest.config_sha256);
    writer.string(manifest.material_name);
    writer.string(manifest.material_definition_json);
    writer.string(manifest.geometry_tier);
    writer.string(manifest.rng_algorithm);
    writer.string(manifest.estimator_flags);
    writer.string(manifest.units);
    writer.string(manifest.observable);
    writer.string(manifest.normalization);

    const EnvFitAxes& axes = block.axes;
    writer.u64(static_cast<std::uint64_t>(axes.energy_axes.size()));
    for (const EnvFitEnergyAxis& axis : axes.energy_axes) {
        writer.string(axis.name);
        writer.u32(axis.transform_version);
        put_vector(writer, axis.registered_edges,
                   [&](double value) { writer.f64(value); });
        put_vector(writer, axis.physical_edges_keV,
                   [&](double value) { writer.f64(value); });
    }
    writer.u32(axes.total_energy_axis);
    put_vector(writer, axes.mu_edges,
               [&](double value) { writer.f64(value); });
    writer.u64(static_cast<std::uint64_t>(axes.lines.size()));
    for (const EnvFitLineDefinition& line : axes.lines) {
        writer.string(line.name);
        writer.f64(line.energy_keV);
        writer.u32(line.origin);
        writer.u32(line.post_creation_order);
    }

    writer.u64(static_cast<std::uint64_t>(block.components.size()));
    for (const EnvFitComponent& component : block.components) {
        writer.string(component.name);
        writer.u32(component.origin);
        writer.u32(component.post_creation_order);
        writer.u32(component.energy_axis);
        put_vector(writer, component.support_mask,
                   [&](std::uint8_t value) { writer.u8(value); });
    }

    const EnvFitNodeMetadata& node = block.node;
    writer.f64(node.e0_keV);
    writer.f64(node.distance_cm);
    writer.f64(node.h_src_cm);
    writer.f64(node.h_det_cm);
    writer.f64(node.transverse_extent_cm);
    for (double value : node.transformed_coordinates) writer.f64(value);
    writer.string(node.role);
    writer.string(node.status);
    writer.f64(node.direct_transmission);
    writer.f64(node.direct_fluence);
    writer.u64(node.histories);
    writer.f64(node.wall_seconds);
    writer.f64(node.cpu_seconds);
    writer.f64(node.target_total_rel_se);
    writer.f64(node.achieved_total_rel_se);

    writer.u32(static_cast<std::uint32_t>(kEnvFitReplicaCount));
    for (const EnvFitReplicaRecord& replica : block.replicas) {
        writer.u32(replica.replica_id);
        writer.u64(replica.history_begin);
        writer.u64(replica.history_count);
        writer.u64(replica.rng_seed);
    }

    put_moments(writer, block.continuum_components);
    put_moments(writer, block.continuum_total);
    put_moments(writer, block.line_components);
    put_moments(writer, block.line_total);
    put_vector(writer, block.replica_sums.continuum_components,
               [&](float value) { writer.f32(value); });
    put_vector(writer, block.replica_sums.continuum_total,
               [&](float value) { writer.f32(value); });
    put_vector(writer, block.replica_sums.line_components,
               [&](float value) { writer.f32(value); });
    put_vector(writer, block.replica_sums.line_total,
               [&](float value) { writer.f32(value); });

    writer.string(block.stopping.definition_json);
    writer.u64(block.stopping.n);
    for (double value : block.stopping.s1) writer.f64(value);
    for (double value : block.stopping.outer_s2) writer.f64(value);

    const EnvFitDiagnostics& diagnostics = block.diagnostics;
    writer.u64(static_cast<std::uint64_t>(diagnostics.path_channels.size()));
    for (const EnvFitDiagnosticChannel& channel : diagnostics.path_channels) {
        writer.string(channel.name);
        put_moments(writer, channel.continuum);
        put_moments(writer, channel.lines);
    }
    writer.u64(static_cast<std::uint64_t>(diagnostics.interaction_names.size()));
    for (const std::string& name : diagnostics.interaction_names)
        writer.string(name);
    put_vector(writer, diagnostics.interaction_counts,
               [&](std::uint64_t value) { writer.u64(value); });
    writer.u64(diagnostics.scattered_total.n);
    writer.f64(diagnostics.scattered_total.s1);
    writer.f64(diagnostics.scattered_total.s2);
    writer.u64(diagnostics.dropped_weight.n);
    writer.f64(diagnostics.dropped_weight.s1);
    writer.f64(diagnostics.dropped_weight.s2);
    writer.u64(diagnostics.transport_failures);
    writer.f64(diagnostics.maximum_history_contribution);
    writer.u8(diagnostics.capped ? 1U : 0U);
    writer.string(diagnostics.failure_reason);
}

bool deserialize_payload(const std::uint8_t* data, std::size_t size,
                         EnvFitBlock& block, std::string* error)
{
    Reader reader(data, size);
    try {
        block.node_id = reader.u64();
        block.block_id = reader.u64();

        EnvFitManifest& manifest = block.manifest;
        manifest.canonical_config_json = reader.string();
        manifest.config_sha256 = reader.string();
        manifest.material_name = reader.string();
        manifest.material_definition_json = reader.string();
        manifest.geometry_tier = reader.string();
        manifest.rng_algorithm = reader.string();
        manifest.estimator_flags = reader.string();
        manifest.units = reader.string();
        manifest.observable = reader.string();
        manifest.normalization = reader.string();

        EnvFitAxes& axes = block.axes;
        const std::uint64_t n_axes = reader.count(28);
        axes.energy_axes.reserve(static_cast<std::size_t>(n_axes));
        for (std::uint64_t i = 0; i < n_axes && reader.ok(); ++i) {
            EnvFitEnergyAxis axis;
            axis.name = reader.string();
            axis.transform_version = reader.u32();
            axis.registered_edges = get_vector<double>(
                reader, 8, [&] { return reader.f64(); });
            axis.physical_edges_keV = get_vector<double>(
                reader, 8, [&] { return reader.f64(); });
            axes.energy_axes.push_back(std::move(axis));
        }
        axes.total_energy_axis = reader.u32();
        axes.mu_edges = get_vector<double>(reader, 8, [&] { return reader.f64(); });
        const std::uint64_t n_lines = reader.count(24);
        axes.lines.reserve(static_cast<std::size_t>(n_lines));
        for (std::uint64_t i = 0; i < n_lines && reader.ok(); ++i) {
            EnvFitLineDefinition line;
            line.name = reader.string();
            line.energy_keV = reader.f64();
            line.origin = reader.u32();
            line.post_creation_order = reader.u32();
            axes.lines.push_back(std::move(line));
        }

        const std::uint64_t n_components = reader.count(28);
        block.components.reserve(static_cast<std::size_t>(n_components));
        for (std::uint64_t i = 0; i < n_components && reader.ok(); ++i) {
            EnvFitComponent component;
            component.name = reader.string();
            component.origin = reader.u32();
            component.post_creation_order = reader.u32();
            component.energy_axis = reader.u32();
            component.support_mask = get_vector<std::uint8_t>(
                reader, 1, [&] { return reader.u8(); });
            block.components.push_back(std::move(component));
        }

        EnvFitNodeMetadata& node = block.node;
        node.e0_keV = reader.f64();
        node.distance_cm = reader.f64();
        node.h_src_cm = reader.f64();
        node.h_det_cm = reader.f64();
        node.transverse_extent_cm = reader.f64();
        for (double& value : node.transformed_coordinates) value = reader.f64();
        node.role = reader.string();
        node.status = reader.string();
        node.direct_transmission = reader.f64();
        node.direct_fluence = reader.f64();
        node.histories = reader.u64();
        node.wall_seconds = reader.f64();
        node.cpu_seconds = reader.f64();
        node.target_total_rel_se = reader.f64();
        node.achieved_total_rel_se = reader.f64();

        const std::uint32_t replica_count = reader.u32();
        if (replica_count != kEnvFitReplicaCount)
            reader.set_error("worker block does not contain exactly 32 replicas");
        for (EnvFitReplicaRecord& replica : block.replicas) {
            replica.replica_id = reader.u32();
            replica.history_begin = reader.u64();
            replica.history_count = reader.u64();
            replica.rng_seed = reader.u64();
        }

        block.continuum_components = get_moments(reader);
        block.continuum_total = get_moments(reader);
        block.line_components = get_moments(reader);
        block.line_total = get_moments(reader);
        block.replica_sums.continuum_components = get_vector<float>(
            reader, 4, [&] { return reader.f32(); });
        block.replica_sums.continuum_total = get_vector<float>(
            reader, 4, [&] { return reader.f32(); });
        block.replica_sums.line_components = get_vector<float>(
            reader, 4, [&] { return reader.f32(); });
        block.replica_sums.line_total = get_vector<float>(
            reader, 4, [&] { return reader.f32(); });

        block.stopping.definition_json = reader.string();
        block.stopping.n = reader.u64();
        for (double& value : block.stopping.s1) value = reader.f64();
        for (double& value : block.stopping.outer_s2) value = reader.f64();

        EnvFitDiagnostics& diagnostics = block.diagnostics;
        const std::uint64_t n_channels = reader.count(24);
        diagnostics.path_channels.reserve(static_cast<std::size_t>(n_channels));
        for (std::uint64_t i = 0; i < n_channels && reader.ok(); ++i) {
            EnvFitDiagnosticChannel channel;
            channel.name = reader.string();
            channel.continuum = get_moments(reader);
            channel.lines = get_moments(reader);
            diagnostics.path_channels.push_back(std::move(channel));
        }
        const std::uint64_t n_interactions = reader.count(8);
        diagnostics.interaction_names.reserve(static_cast<std::size_t>(n_interactions));
        for (std::uint64_t i = 0; i < n_interactions && reader.ok(); ++i)
            diagnostics.interaction_names.push_back(reader.string());
        diagnostics.interaction_counts = get_vector<std::uint64_t>(
            reader, 8, [&] { return reader.u64(); });
        diagnostics.scattered_total.n = reader.u64();
        diagnostics.scattered_total.s1 = reader.f64();
        diagnostics.scattered_total.s2 = reader.f64();
        diagnostics.dropped_weight.n = reader.u64();
        diagnostics.dropped_weight.s1 = reader.f64();
        diagnostics.dropped_weight.s2 = reader.f64();
        diagnostics.transport_failures = reader.u64();
        diagnostics.maximum_history_contribution = reader.f64();
        const std::uint8_t capped = reader.u8();
        if (capped > 1U) reader.set_error("invalid capped flag");
        diagnostics.capped = capped != 0;
        diagnostics.failure_reason = reader.string();
    } catch (const std::exception& exception) {
        return fail(error, std::string("cannot allocate/parse payload: ") + exception.what());
    }

    if (!reader.ok()) return fail(error, reader.error());
    if (!reader.done()) return fail(error, "payload contains trailing bytes");
    return true;
}

bool ascending_finite(const std::vector<double>& values)
{
    if (values.size() < 2) return false;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!finite(values[i])) return false;
        if (i != 0 && !(values[i] > values[i - 1])) return false;
    }
    return true;
}

bool validate_moments(const EnvFitMomentArray& moments, std::size_t expected,
                      std::uint64_t histories, const std::string& name,
                      std::string* error)
{
    if (moments.n.size() != expected || moments.s1.size() != expected ||
        moments.s2.size() != expected)
        return fail(error, name + " tensor has the wrong dimensions");
    for (std::size_t i = 0; i < expected; ++i) {
        if (moments.n[i] != histories)
            return fail(error, name + " N does not equal block history count");
        if (!finite(moments.s1[i]) || !finite(moments.s2[i]) ||
            moments.s1[i] < 0.0 || moments.s2[i] < 0.0)
            return fail(error, name + " contains a negative or non-finite moment");
        if (histories != 0) {
            const long double minimum =
                static_cast<long double>(moments.s1[i]) * moments.s1[i] /
                static_cast<long double>(histories);
            const long double tolerance = 128.0L *
                std::numeric_limits<double>::epsilon() *
                std::max(1.0L, std::fabs(minimum));
            if (static_cast<long double>(moments.s2[i]) + tolerance < minimum)
                return fail(error, name + " has S2 < S1^2/N");
        }
    }
    return true;
}

bool validate_replica_sums(const std::vector<float>& values,
                           std::size_t cells, const std::string& name,
                           std::string* error)
{
    std::size_t expected = 0;
    if (!multiply_size(kEnvFitReplicaCount, cells, expected) ||
        values.size() != expected)
        return fail(error, name + " replica tensor has the wrong dimensions");
    for (float value : values) {
        if (!finite(value) || value < 0.0f)
            return fail(error, name + " replica tensor contains invalid values");
    }
    return true;
}

bool validate_replica_closure(const std::vector<float>& values,
                              std::size_t cells,
                              const EnvFitMomentArray& moments,
                              const std::string& name,
                              std::string* error)
{
    for (std::size_t cell = 0; cell < cells; ++cell) {
        long double replica_sum = 0.0L;
        long double absolute_sum = 0.0L;
        for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
            const long double value = values[replica * cells + cell];
            replica_sum += value;
            absolute_sum += std::fabs(value);
        }
        const long double expected = moments.s1[cell];
        const long double scale =
            std::max(std::fabs(expected), absolute_sum);
        const long double tolerance =
            128.0L * std::numeric_limits<float>::epsilon() * scale +
            static_cast<long double>(kEnvFitReplicaCount) *
                std::numeric_limits<float>::denorm_min();
        if (std::fabs(replica_sum - expected) > tolerance)
            return fail(error, name + " replica S1 maps do not close to raw S1");
    }
    return true;
}

bool validate_line_replica_additivity(const EnvFitBlock& block,
                                      std::string* error)
{
    const std::size_t total_cells = env_fit_line_total_cells(block);
    const std::size_t component_cells = env_fit_line_component_cells(block);
    if (total_cells == 0) return true;
    for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
        for (std::size_t cell = 0; cell < total_cells; ++cell) {
            long double component_sum = 0.0L;
            long double absolute_sum = 0.0L;
            for (std::size_t component = 0;
                 component < block.components.size(); ++component) {
                const long double value =
                    block.replica_sums.line_components[
                        replica * component_cells + component * total_cells + cell];
                component_sum += value;
                absolute_sum += std::fabs(value);
            }
            const long double total = block.replica_sums.line_total[
                replica * total_cells + cell];
            const long double scale =
                std::max({std::fabs(total), std::fabs(component_sum), absolute_sum});
            const long double tolerance =
                128.0L * std::numeric_limits<float>::epsilon() * scale +
                static_cast<long double>(block.components.size() + 1) *
                    std::numeric_limits<float>::denorm_min();
            if (std::fabs(total - component_sum) > tolerance)
                return fail(error,
                            "replica total-line S1 is not the sum of components");
        }
    }
    return true;
}

bool same_manifest(const EnvFitManifest& a, const EnvFitManifest& b)
{
    return a.canonical_config_json == b.canonical_config_json &&
           a.config_sha256 == b.config_sha256 &&
           a.material_name == b.material_name &&
           a.material_definition_json == b.material_definition_json &&
           a.geometry_tier == b.geometry_tier &&
           a.rng_algorithm == b.rng_algorithm &&
           a.estimator_flags == b.estimator_flags && a.units == b.units &&
           a.observable == b.observable && a.normalization == b.normalization;
}

bool same_axes(const EnvFitAxes& a, const EnvFitAxes& b)
{
    if (a.total_energy_axis != b.total_energy_axis ||
        a.mu_edges != b.mu_edges || a.energy_axes.size() != b.energy_axes.size() ||
        a.lines.size() != b.lines.size())
        return false;
    for (std::size_t i = 0; i < a.energy_axes.size(); ++i) {
        const auto& x = a.energy_axes[i];
        const auto& y = b.energy_axes[i];
        if (x.name != y.name || x.transform_version != y.transform_version ||
            x.registered_edges != y.registered_edges ||
            x.physical_edges_keV != y.physical_edges_keV)
            return false;
    }
    for (std::size_t i = 0; i < a.lines.size(); ++i) {
        const auto& x = a.lines[i];
        const auto& y = b.lines[i];
        if (x.name != y.name || x.energy_keV != y.energy_keV ||
            x.origin != y.origin ||
            x.post_creation_order != y.post_creation_order)
            return false;
    }
    return true;
}

bool same_components(const std::vector<EnvFitComponent>& a,
                     const std::vector<EnvFitComponent>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].name != b[i].name || a[i].origin != b[i].origin ||
            a[i].post_creation_order != b[i].post_creation_order ||
            a[i].energy_axis != b[i].energy_axis ||
            a[i].support_mask != b[i].support_mask)
            return false;
    }
    return true;
}

bool same_node_identity(const EnvFitNodeMetadata& a,
                        const EnvFitNodeMetadata& b)
{
    return a.e0_keV == b.e0_keV && a.distance_cm == b.distance_cm &&
           a.h_src_cm == b.h_src_cm && a.h_det_cm == b.h_det_cm &&
           a.transverse_extent_cm == b.transverse_extent_cm &&
           a.transformed_coordinates == b.transformed_coordinates &&
           a.role == b.role &&
           a.direct_transmission == b.direct_transmission &&
           a.direct_fluence == b.direct_fluence &&
           a.target_total_rel_se == b.target_total_rel_se;
}

bool same_diagnostic_layout(const EnvFitDiagnostics& a,
                            const EnvFitDiagnostics& b)
{
    if (a.path_channels.size() != b.path_channels.size() ||
        a.interaction_names != b.interaction_names)
        return false;
    for (std::size_t i = 0; i < a.path_channels.size(); ++i)
        if (a.path_channels[i].name != b.path_channels[i].name) return false;
    return true;
}

bool ranges_overlap(const EnvFitReplicaRecord& a,
                    const EnvFitReplicaRecord& b)
{
    const std::uint64_t a_end = a.history_begin + a.history_count;
    const std::uint64_t b_end = b.history_begin + b.history_count;
    return a.history_begin < b_end && b.history_begin < a_end;
}

void add_moments(EnvFitMomentArray& destination,
                 const EnvFitMomentArray& source)
{
    for (std::size_t i = 0; i < destination.n.size(); ++i) {
        destination.n[i] += source.n[i];
        destination.s1[i] += source.s1[i];
        destination.s2[i] += source.s2[i];
    }
}

template<class T>
void add_vector(std::vector<T>& destination, const std::vector<T>& source)
{
    for (std::size_t i = 0; i < destination.size(); ++i)
        destination[i] += source[i];
}

} // namespace

} // namespace ceelo
