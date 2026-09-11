#include "ShortestTravelTimes.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>

namespace spdp {
namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

void fnv_mix(std::uint64_t& hash, std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
}

std::uint64_t double_bits(double value) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

}  // namespace

double ShortestTravelTimes::between(int from, int to) const {
    return distance.at(static_cast<std::size_t>(from)).at(static_cast<std::size_t>(to));
}

bool ShortestTravelTimes::reachable(int from, int to) const {
    return std::isfinite(between(from, to));
}

std::vector<int> ShortestTravelTimes::path(int from, int to) const {
    if (!reachable(from, to)) {
        return {};
    }
    std::vector<int> reversed;
    int current = to;
    reversed.push_back(current);
    while (current != from) {
        current = predecessor[static_cast<std::size_t>(from)][static_cast<std::size_t>(current)];
        if (current < 0) {
            throw std::runtime_error("Shortest-path predecessor chain is broken.");
        }
        reversed.push_back(current);
        if (reversed.size() > location_count + 1U) {
            throw std::runtime_error("Shortest-path predecessor chain contains a cycle.");
        }
    }
    std::reverse(reversed.begin(), reversed.end());
    return reversed;
}

std::uint64_t ShortestTravelTimes::fingerprint() const {
    std::uint64_t hash = 1469598103934665603ULL;
    fnv_mix(hash, static_cast<std::uint64_t>(location_count));
    for (const std::vector<double>& row : distance) {
        for (double value : row) {
            fnv_mix(hash, double_bits(std::isfinite(value) ? value : -1.0));
        }
    }
    return hash;
}

ShortestTravelTimes compute_shortest_travel_times(
    const std::vector<std::vector<double>>& travel_time
) {
    const std::size_t n = travel_time.size();
    for (const std::vector<double>& row : travel_time) {
        if (row.size() != n) {
            throw std::runtime_error("Travel-time matrix must be square.");
        }
        for (double value : row) {
            if (!std::isfinite(value) || value < 0.0) {
                throw std::runtime_error(
                    "Travel-time matrix entries must be finite and nonnegative."
                );
            }
        }
    }

    ShortestTravelTimes result;
    result.location_count = n;
    result.distance.assign(n, std::vector<double>(n, kInfinity));
    result.predecessor.assign(n, std::vector<int>(n, -1));

    using Item = std::pair<double, int>;
    for (std::size_t source = 0; source < n; ++source) {
        std::vector<double>& dist = result.distance[source];
        std::vector<int>& pred = result.predecessor[source];
        std::vector<bool> settled(n, false);
        dist[source] = 0.0;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
        heap.emplace(0.0, static_cast<int>(source));
        while (!heap.empty()) {
            const auto [d, u_int] = heap.top();
            heap.pop();
            const auto u = static_cast<std::size_t>(u_int);
            if (settled[u] || d > dist[u]) {
                continue;
            }
            settled[u] = true;
            for (std::size_t v = 0; v < n; ++v) {
                if (v == u) {
                    continue;
                }
                const double candidate = d + travel_time[u][v];
                if (candidate < dist[v] - 1e-12) {
                    dist[v] = candidate;
                    pred[v] = u_int;
                    heap.emplace(candidate, static_cast<int>(v));
                }
            }
        }
    }

    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = 0; b < n; ++b) {
            if (a == b) {
                continue;
            }
            const double reduction = travel_time[a][b] - result.distance[a][b];
            if (reduction > 1e-9) {
                ++result.non_metric_pair_count;
                result.max_reduction = std::max(result.max_reduction, reduction);
            }
        }
    }
    return result;
}

ShortestTravelTimes compute_physical_shortest_travel_times(const SPDPData& data) {
    const auto physical = static_cast<std::size_t>(std::max(0, data.locations));
    if (data.time.size() < physical) {
        throw std::runtime_error(
            "The instance travel-time matrix is smaller than its location count."
        );
    }
    std::vector<std::vector<double>> block(physical, std::vector<double>(physical, 0.0));
    for (std::size_t a = 0; a < physical; ++a) {
        if (data.time[a].size() < physical) {
            throw std::runtime_error(
                "The instance travel-time matrix row is shorter than its location count."
            );
        }
        for (std::size_t b = 0; b < physical; ++b) {
            block[a][b] = data.time[a][b];
        }
    }
    return compute_shortest_travel_times(block);
}

std::vector<std::vector<double>> metric_closure_travel_time_matrix(
    const SPDPData& data,
    const ShortestTravelTimes& shortest
) {
    const auto physical = static_cast<std::size_t>(std::max(0, data.locations));
    if (shortest.location_count != physical) {
        throw std::runtime_error(
            "Shortest travel times do not match the instance location count."
        );
    }
    std::vector<std::vector<double>> matrix = data.time;
    for (std::size_t a = 0; a < physical; ++a) {
        for (std::size_t b = 0; b < physical; ++b) {
            if (a == b) {
                continue;
            }
            const double value = shortest.distance[a][b];
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    "Metric closure requires every physical location pair to be reachable."
                );
            }
            matrix[a][b] = value;
        }
    }
    return matrix;
}

}  // namespace spdp

