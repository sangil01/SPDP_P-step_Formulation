#include "TypeTravelTimeCoverCuts.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <utility>

#include "KMinComputation.h"

namespace spdp {
namespace {

class NullBuffer final : public std::streambuf {
protected:
    int overflow(int ch) override { return ch; }
};

void fnv_mix(std::uint64_t& hash, std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
}

std::uint64_t double_bits(double value) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

std::uint64_t instance_fingerprint(
    const SPDPData& data,
    const ShortestTravelTimes& shortest,
    const TypeTravelTimeCoverOptions& options
) {
    std::uint64_t hash = 1469598103934665603ULL;
    fnv_mix(hash, 0x54544331ULL);  // "TTC1": certificate format version
    fnv_mix(hash, static_cast<std::uint64_t>(data.requests.size()));
    for (const Request& request : data.requests) {
        fnv_mix(hash, static_cast<std::uint64_t>(request.from_id));
        fnv_mix(hash, static_cast<std::uint64_t>(request.to_id));
        fnv_mix(hash, static_cast<std::uint64_t>(request.container_type));
    }
    fnv_mix(hash, double_bits(data.time_limit));
    fnv_mix(hash, double_bits(data.time_pickup));
    fnv_mix(hash, double_bits(data.time_empty));
    fnv_mix(hash, double_bits(data.time_delivery));
    fnv_mix(hash, static_cast<std::uint64_t>(data.locations));
    fnv_mix(hash, shortest.fingerprint());
    // The formulation used to certify the bound is part of the key so that a
    // certificate produced by one duration model is never reused by another.
    fnv_mix(hash, options.add_time_flow_formulation ? 1ULL : 0ULL);
    fnv_mix(hash, options.time_flow_state_disaggregated ? 1ULL : 0ULL);
    fnv_mix(hash, options.add_time_constraints ? 1ULL : 0ULL);
    return hash;
}

std::string hex(std::uint64_t value) {
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

std::vector<int> mask_types(const std::vector<int>& types, std::uint64_t mask) {
    std::vector<int> result;
    for (std::size_t i = 0; i < types.size(); ++i) {
        if ((mask >> i) & 1ULL) {
            result.push_back(types[i]);
        }
    }
    return result;
}

struct CachedCertificate {
    int rho = 0;
    double lower_bound = 0.0;
    double upper_bound = 0.0;
};

std::filesystem::path cache_path(const std::string& cache_dir, std::uint64_t fingerprint) {
    return std::filesystem::path(cache_dir) / ("ttcover_" + hex(fingerprint) + ".txt");
}

std::map<std::uint64_t, CachedCertificate> load_cache(
    const std::string& cache_dir,
    std::uint64_t fingerprint
) {
    std::map<std::uint64_t, CachedCertificate> cache;
    if (cache_dir.empty()) {
        return cache;
    }
    std::ifstream in(cache_path(cache_dir, fingerprint));
    if (!in) {
        return cache;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream fields(line);
        std::uint64_t mask = 0;
        CachedCertificate certificate;
        if (fields >> mask >> certificate.rho >> certificate.lower_bound >> certificate.upper_bound) {
            cache[mask] = certificate;
        }
    }
    return cache;
}

void append_cache(
    const std::string& cache_dir,
    std::uint64_t fingerprint,
    const TypeTravelTimeCoverEntry& entry
) {
    if (cache_dir.empty() || !entry.certified) {
        return;
    }
    std::filesystem::create_directories(cache_dir);
    const std::filesystem::path path = cache_path(cache_dir, fingerprint);
    const bool fresh = !std::filesystem::exists(path);
    std::ofstream out(path, std::ios::app);
    if (!out) {
        return;
    }
    if (fresh) {
        out << "# type travel-time cover certificates; fingerprint " << hex(fingerprint) << '\n';
        out << "# mask rho duration_lb duration_ub\n";
    }
    out << entry.type_mask << ' ' << entry.rho << ' '
        << std::setprecision(17) << entry.duration_lower_bound << ' '
        << entry.duration_upper_bound << '\n';
}

}  // namespace

const char* type_travel_time_cover_mode_name(TypeTravelTimeCoverMode mode) {
    switch (mode) {
        case TypeTravelTimeCoverMode::Static: return "static";
        case TypeTravelTimeCoverMode::Screened: return "screened";
    }
    return "unknown";
}

std::vector<int> instance_container_types(const SPDPData& data) {
    std::set<int> types;
    for (const Request& request : data.requests) {
        types.insert(request.container_type);
    }
    return std::vector<int>(types.begin(), types.end());
}

SPDPData build_type_restricted_instance(
    const SPDPData& data,
    const std::vector<std::vector<double>>& closure_travel_time,
    const std::vector<int>& types,
    std::uint64_t type_mask
) {
    SPDPData restricted = data;
    restricted.time = closure_travel_time;
    restricted.requests.clear();
    for (const Request& request : data.requests) {
        const auto found = std::find(types.begin(), types.end(), request.container_type);
        if (found == types.end()) {
            throw std::runtime_error("Request type is missing from the type list.");
        }
        const auto bit = static_cast<std::size_t>(found - types.begin());
        if ((type_mask >> bit) & 1ULL) {
            restricted.requests.push_back(request);
        }
    }
    return restricted;
}

TypeTravelTimeCoverResult compute_type_travel_time_cover_bounds(
    const SPDPData& data,
    const TypeTravelTimeCoverOptions& options
) {
    const auto overall_start = std::chrono::steady_clock::now();
    TypeTravelTimeCoverResult result;
    result.types = instance_container_types(data);
    if (result.types.size() >= 63U) {
        throw std::runtime_error("Too many container types for a 64-bit type mask.");
    }

    const ShortestTravelTimes shortest = compute_physical_shortest_travel_times(data);
    const std::vector<std::vector<double>> closure =
        metric_closure_travel_time_matrix(data, shortest);
    result.shortest_time_fingerprint = shortest.fingerprint();
    result.non_metric_pair_count = shortest.non_metric_pair_count;
    result.max_travel_time_reduction = shortest.max_reduction;
    result.instance_fingerprint = instance_fingerprint(data, shortest, options);

    const std::map<std::uint64_t, CachedCertificate> cache =
        load_cache(options.cache_dir, result.instance_fingerprint);

    const std::size_t type_count = result.types.size();
    const std::uint64_t full_mask = (std::uint64_t{1} << type_count) - 1ULL;
    std::vector<std::uint64_t> masks;
    for (std::uint64_t mask = 1; mask <= full_mask; ++mask) {
        const auto size = static_cast<std::size_t>(__builtin_popcountll(mask));
        if (options.max_type_set_size > 0U && size > options.max_type_set_size) {
            continue;
        }
        masks.push_back(mask);
    }
    std::sort(masks.begin(), masks.end(), [](std::uint64_t lhs, std::uint64_t rhs) {
        const int lp = __builtin_popcountll(lhs);
        const int rp = __builtin_popcountll(rhs);
        if (lp != rp) {
            return lp < rp;
        }
        return lhs < rhs;
    });

    NullBuffer null_buffer;
    std::ostream null_stream(&null_buffer);

    for (std::uint64_t mask : masks) {
        TypeTravelTimeCoverEntry entry;
        entry.type_mask = mask;
        entry.types = mask_types(result.types, mask);
        for (const Request& request : data.requests) {
            if (std::find(entry.types.begin(), entry.types.end(), request.container_type) !=
                entry.types.end()) {
                ++entry.request_count;
            }
        }

        if (mask == full_mask) {
            entry.skipped = true;
            entry.skip_reason = "full type set (departure bound covered by VI44)";
            result.entries.push_back(std::move(entry));
            continue;
        }

        const auto cached = cache.find(mask);
        if (cached != cache.end()) {
            entry.rho = cached->second.rho;
            entry.duration_lower_bound = cached->second.lower_bound;
            entry.duration_upper_bound = cached->second.upper_bound;
            entry.has_lower_bound = true;
            entry.certified = true;
            entry.from_cache = true;
            ++result.cache_hits;
            result.entries.push_back(std::move(entry));
            continue;
        }

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - overall_start).count();
        double remaining = options.subproblem_time_limit;
        if (options.total_time_limit > 0.0) {
            remaining = std::min(remaining, options.total_time_limit - elapsed);
            if (remaining <= 1e-3) {
                result.total_time_limit_hit = true;
                entry.skipped = true;
                entry.skip_reason = "preprocessing time budget exhausted";
                result.entries.push_back(std::move(entry));
                continue;
            }
        }

        const SPDPData restricted =
            build_type_restricted_instance(data, closure, result.types, mask);
        const MultiDiGraph restricted_graph =
            build_multigraph(restricted, GraphBuildOptions{}, &null_stream);
        const DerivedMultiGraph duration_graph = derive_multigraph(
            restricted_graph,
            GraphPurpose::DurationBound,
            DerivedGraphOptions{true, true}
        );
        entry.restricted_edge_count = duration_graph.graph.number_of_edges();
        entry.restricted_graph_fingerprint = duration_graph.stats.fingerprint;

        DurationRoundedBoundOptions bound_options;
        bound_options.solver_time_limit = remaining;
        bound_options.add_time_constraints = options.add_time_constraints;
        bound_options.add_time_flow_formulation = options.add_time_flow_formulation;
        bound_options.time_flow_state_disaggregated = options.time_flow_state_disaggregated;
        bound_options.rounded_bound_stop = true;
        bound_options.gurobi_threads = options.gurobi_threads;
        bound_options.name_prefix = "ttcover_" + std::to_string(mask);
        const DurationRoundedBoundResult bound =
            solve_duration_rounded_bound(restricted, duration_graph.graph, bound_options);

        entry.status = bound.status;
        entry.runtime_seconds = bound.runtime_seconds;
        entry.hit_time_limit = bound.hit_time_limit;
        entry.stopped_by_rounded_bound = bound.stopped_by_rounded_bound;
        entry.duration_upper_bound = bound.exact_upper_bound;
        if (bound.has_lower_bound) {
            entry.has_lower_bound = true;
            entry.duration_lower_bound = bound.safe_lower_bound;
            entry.rho = bound.rounded_lower_bound;
            entry.certified = bound.certified;
        }
        append_cache(options.cache_dir, result.instance_fingerprint, entry);
        result.entries.push_back(std::move(entry));
    }

    result.total_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - overall_start).count();
    if (options.log_stream != nullptr) {
        write_type_travel_time_cover_log(*options.log_stream, result);
    }
    return result;
}

std::vector<TypeTravelTimeCoverRow> build_type_travel_time_cover_rows(
    const SPDPData& data,
    const MultiDiGraph& main_graph,
    const TypeTravelTimeCoverResult& result,
    int min_rho
) {
    std::vector<TypeTravelTimeCoverRow> rows;
    const std::size_t request_count = data.requests.size();
    for (const TypeTravelTimeCoverEntry& entry : result.entries) {
        if (entry.skipped || !entry.certified || entry.rho < std::max(1, min_rho)) {
            continue;
        }
        std::vector<bool> in_set(main_graph.number_of_nodes(), false);
        TypeTravelTimeCoverRow row;
        row.type_mask = entry.type_mask;
        row.rho = entry.rho;
        for (std::size_t request_index = 0; request_index < request_count; ++request_index) {
            const int type = data.requests[request_index].container_type;
            if (std::find(entry.types.begin(), entry.types.end(), type) == entry.types.end()) {
                continue;
            }
            const NodeId pickup = static_cast<NodeId>(request_index + 1U);
            const NodeId delivery = static_cast<NodeId>(request_count + request_index + 1U);
            row.nodes.push_back(pickup);
            row.nodes.push_back(delivery);
            in_set[static_cast<std::size_t>(pickup)] = true;
            in_set[static_cast<std::size_t>(delivery)] = true;
        }
        row.row.name = "ttcover_mask" + std::to_string(entry.type_mask) + "_rho" +
            std::to_string(entry.rho);
        row.row.sense = PstepValidInequalitySense::GreaterEqual;
        row.row.rhs = static_cast<double>(entry.rho);
        for (std::size_t edge_id = 0; edge_id < main_graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = main_graph.edges()[edge_id];
            if (in_set[static_cast<std::size_t>(edge.v)] &&
                !in_set[static_cast<std::size_t>(edge.u)]) {
                row.row.edge_terms.emplace_back(edge_id, 1.0);
            }
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

std::vector<std::size_t> screen_type_travel_time_cover_rows(
    const std::vector<TypeTravelTimeCoverRow>& rows,
    const std::vector<double>& y_values,
    const std::vector<bool>& already_added,
    double violation_tolerance,
    std::size_t max_rows
) {
    std::vector<std::pair<double, std::size_t>> violated;
    for (std::size_t index = 0; index < rows.size(); ++index) {
        if (index < already_added.size() && already_added[index]) {
            continue;
        }
        double lhs = 0.0;
        for (const auto& term : rows[index].row.edge_terms) {
            lhs += term.second * y_values[term.first];
        }
        const double violation = rows[index].row.rhs - lhs;
        if (violation > violation_tolerance) {
            violated.emplace_back(violation, index);
        }
    }
    std::sort(violated.begin(), violated.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first > rhs.first;
    });
    std::vector<std::size_t> indices;
    for (const auto& [violation, index] : violated) {
        if (indices.size() >= max_rows) {
            break;
        }
        indices.push_back(index);
    }
    return indices;
}

void write_type_travel_time_cover_log(
    std::ostream& out,
    const TypeTravelTimeCoverResult& result
) {
    out << "[ttcover] types=" << result.types.size()
        << " entries=" << result.entries.size()
        << " cache_hits=" << result.cache_hits
        << " non_metric_pairs=" << result.non_metric_pair_count
        << " max_travel_time_reduction=" << result.max_travel_time_reduction
        << " shortest_time_fingerprint=" << hex(result.shortest_time_fingerprint)
        << " instance_fingerprint=" << hex(result.instance_fingerprint)
        << " total_seconds=" << result.total_seconds
        << " budget_hit=" << (result.total_time_limit_hit ? 1 : 0) << '\n';
    for (const TypeTravelTimeCoverEntry& entry : result.entries) {
        out << "[ttcover] mask=" << entry.type_mask << " types=";
        for (std::size_t i = 0; i < entry.types.size(); ++i) {
            out << (i == 0 ? "" : ",") << entry.types[i];
        }
        out << " requests=" << entry.request_count;
        if (entry.skipped) {
            out << " skipped=1 reason=\"" << entry.skip_reason << "\"\n";
            continue;
        }
        out << " edges=" << entry.restricted_edge_count
            << " duration_lb=" << entry.duration_lower_bound
            << " duration_ub=" << entry.duration_upper_bound
            << " rho=" << entry.rho
            << " certified=" << (entry.certified ? 1 : 0)
            << " rounded_stop=" << (entry.stopped_by_rounded_bound ? 1 : 0)
            << " time_limit=" << (entry.hit_time_limit ? 1 : 0)
            << " cache=" << (entry.from_cache ? 1 : 0)
            << " status=" << entry.status
            << " seconds=" << entry.runtime_seconds << '\n';
    }
}

}  // namespace spdp

