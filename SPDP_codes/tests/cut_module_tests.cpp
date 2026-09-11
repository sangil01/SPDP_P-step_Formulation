// Module tests for the capacity blossom / type travel-time cover cut families
// and the shared shortest-travel-time layer. Run through CTest
// (spdp_cut_module_tests); requires the SPDP_data directory (instance A8).
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "CapacityBlossomCuts.h"
#include "GenMultiGraph.h"
#include "ReadData.h"
#include "ShortestTravelTimes.h"
#include "TwoIndexSolver.h"
#include "TypeTravelTimeCoverCuts.h"

namespace {

int g_failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++g_failures;
        std::cerr << "FAIL: " << message << '\n';
    } else {
        std::cout << "ok: " << message << '\n';
    }
}

void test_shortest_travel_times_small() {
    // Asymmetric, non-metric: 0->2 direct 100, via 1 costs 20 + 5 = 25.
    const std::vector<std::vector<double>> t = {
        {0, 20, 100, 7},
        {30, 0, 5, 50},
        {60, 9, 0, 3},
        {1, 80, 2, 0},
    };
    const spdp::ShortestTravelTimes s = spdp::compute_shortest_travel_times(t);
    check(std::abs(s.between(0, 2) - 9.0) < 1e-12, "shortest 0->2 uses 0->3->2 (7 + 2 = 9)");
    check(std::abs(s.between(1, 0) - 9.0) < 1e-12, "shortest 1->0 uses 1->2->3->0 (5 + 3 + 1 = 9)");
    check(std::abs(s.between(0, 0)) < 1e-12, "diagonal is zero");
    const std::vector<int> path = s.path(1, 0);
    check(path == std::vector<int>{1, 2, 3, 0}, "path reconstruction 1->2->3->0");
    double path_length = 0.0;
    for (std::size_t i = 1; i < path.size(); ++i) {
        path_length += t[static_cast<std::size_t>(path[i - 1])][static_cast<std::size_t>(path[i])];
    }
    check(std::abs(path_length - s.between(1, 0)) < 1e-12, "reconstructed path length equals distance");
    check(s.non_metric_pair_count > 0, "non-metric pairs detected");
    // Floyd--Warshall reference.
    std::vector<std::vector<double>> fw = t;
    for (std::size_t k = 0; k < 4; ++k) {
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 4; ++j) {
                fw[i][j] = std::min(fw[i][j], fw[i][k] + fw[k][j]);
            }
        }
    }
    bool same = true;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            same = same && std::abs(fw[i][j] - s.between(static_cast<int>(i), static_cast<int>(j))) < 1e-12;
        }
    }
    check(same, "Dijkstra matches Floyd--Warshall");
}

void test_metric_closure_instance(const spdp::SPDPData& data) {
    const spdp::ShortestTravelTimes s = spdp::compute_physical_shortest_travel_times(data);
    const std::vector<std::vector<double>> closure = spdp::metric_closure_travel_time_matrix(data, s);
    const auto L = static_cast<std::size_t>(data.locations);
    check(closure.size() == data.time.size(), "closure keeps the augmented matrix size");
    bool not_longer = true;
    bool triangle = true;
    for (std::size_t a = 0; a < L; ++a) {
        for (std::size_t b = 0; b < L; ++b) {
            not_longer = not_longer && closure[a][b] <= data.time[a][b] + 1e-9;
            for (std::size_t c = 0; c < L; ++c) {
                triangle = triangle && closure[a][b] <= closure[a][c] + closure[c][b] + 1e-9;
            }
        }
    }
    check(not_longer, "closure never exceeds the raw travel time");
    check(triangle, "closure satisfies the triangle inequality");
    bool depot_zero = true;
    for (std::size_t a = 0; a <= L; ++a) {
        depot_zero = depot_zero && closure[a][L] == data.time[a][L] && closure[L][a] == data.time[L][a];
    }
    check(depot_zero, "virtual depot row/column untouched");
    std::cout << "info: non-metric pairs = " << s.non_metric_pair_count
              << ", max reduction = " << s.max_reduction << '\n';
}

// Random fractional pickup/delivery adjacency with z-degree <= 1 on the real
// multigraph edge set; exact separator versus exhaustive odd-set enumeration.
void test_blossom_exactness(const spdp::MultiDiGraph& graph) {
    const spdp::CapacityBlossomSeparator separator = spdp::build_capacity_blossom_separator(graph);
    std::size_t treatment_pair_edges = 0;
    for (const auto& entry : separator.pickup.pair_edges) {
        for (std::size_t edge_id : entry.second) {
            if (!graph.edges()[edge_id].data.sequence_pi.empty()) {
                ++treatment_pair_edges;
            }
        }
    }
    check(treatment_pair_edges > 0, "pickup-pickup edges with embedded treatment are counted as internal");

    std::mt19937 rng(20260912);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (int family_index = 0; family_index < 2; ++family_index) {
        const spdp::CapacityFamily family = family_index == 0
            ? spdp::CapacityFamily::Pickup : spdp::CapacityFamily::Delivery;
        const spdp::CapacityFamilyData& fam = family_index == 0 ? separator.pickup : separator.delivery;
        const std::size_t n = fam.nodes.size();
        int matched = 0;
        int violated_cases = 0;
        for (int trial = 0; trial < 25; ++trial) {
            std::vector<double> y(graph.number_of_edges(), 0.0);
            const auto pair_key = [n](std::size_t i, std::size_t j) {
                if (i > j) std::swap(i, j);
                return i * n + j;
            };
            // Half of the trials plant an odd cycle with weights near 1/2 (a
            // violated blossom), the rest are sparse random points.
            if (trial % 2 == 0) {
                const std::size_t k = 3 + 2 * static_cast<std::size_t>(rng() % 3);  // 3, 5, 7
                std::vector<std::size_t> order(n);
                for (std::size_t i = 0; i < n; ++i) order[i] = i;
                std::shuffle(order.begin(), order.end(), rng);
                for (std::size_t c = 0; c < k && k <= n; ++c) {
                    const auto found = fam.pair_edges.find(pair_key(order[c], order[(c + 1) % k]));
                    if (found == fam.pair_edges.end()) continue;
                    y[found->second[static_cast<std::size_t>(rng() % found->second.size())]] =
                        0.5 - 0.08 * uniform(rng);
                }
                for (const auto& entry : fam.pair_edges) {
                    if (uniform(rng) < 0.05) {
                        y[entry.second[static_cast<std::size_t>(rng() % entry.second.size())]] += 0.1 * uniform(rng);
                    }
                }
            } else {
                for (const auto& entry : fam.pair_edges) {
                    if (uniform(rng) < 0.35) {
                        y[entry.second[static_cast<std::size_t>(rng() % entry.second.size())]] = uniform(rng);
                    }
                }
            }
            // Scale down until every node has z-degree <= 1.
            for (int pass = 0; pass < 50; ++pass) {
                std::vector<double> degree(n, 0.0);
                for (const auto& entry : fam.pair_edges) {
                    double w = 0.0;
                    for (std::size_t e : entry.second) w += y[e];
                    degree[entry.first / n] += w;
                    degree[entry.first % n] += w;
                }
                bool ok = true;
                for (const auto& entry : fam.pair_edges) {
                    const double d = std::max(degree[entry.first / n], degree[entry.first % n]);
                    if (d > 1.0 + 1e-12) {
                        ok = false;
                        for (std::size_t e : entry.second) y[e] /= d;
                    }
                }
                if (ok) break;
            }
            const double reference = spdp::max_capacity_blossom_violation_by_enumeration(
                separator, family, y, 20);
            spdp::CapacityBlossomOptions options;
            options.enabled = true;
            options.min_violation = 1e-9;
            options.max_overlap_jaccard = 1.1;  // keep every candidate
            spdp::CapacityBlossomPool pool;
            spdp::CapacityBlossomStats stats;
            const std::vector<spdp::CapacityBlossomCut> cuts = spdp::separate_capacity_blossom_cuts(
                separator, y, options, 1000, pool, stats);
            double found = -1e300;
            bool rows_consistent = true;
            for (const spdp::CapacityBlossomCut& cut : cuts) {
                if (cut.family != family) continue;
                found = std::max(found, cut.violation);
                double lhs = 0.0;
                for (const auto& term : cut.row.edge_terms) lhs += term.second * y[term.first];
                rows_consistent = rows_consistent &&
                    std::abs((lhs - cut.row.rhs) - cut.violation) < 1e-9 &&
                    cut.compact_set.size() % 2 == 1 && cut.compact_set.size() >= 3;
            }
            if (!rows_consistent) {
                check(false, "capacity blossom row violation equals reported violation");
            }
            if (reference > 1e-9) {
                ++violated_cases;
                if (std::abs(found - reference) < 1e-7) ++matched;
                else std::cerr << "mismatch: reference " << reference << " found " << found << '\n';
            } else if (found > 1e-7) {
                std::cerr << "spurious cut: reference " << reference << " found " << found << '\n';
            } else {
                ++matched;
            }
        }
        std::ostringstream label;
        label << spdp::capacity_family_name(family) << " family: separator matches enumeration on 25 random points ("
              << violated_cases << " violated)";
        check(matched == 25, label.str());
    }
}

// Validity on an integer solution: no blossom cut and no certified type cover
// row may be violated by a feasible integer point.
void test_validity_on_integer_solution(const spdp::SPDPData& data, const spdp::MultiDiGraph& graph) {
    spdp::DirectTwoIndexOptions options;
    options.model_type = spdp::DirectTwoIndexModelType::IP;
    options.objective = spdp::DirectTwoIndexObjective::OriginalCost;
    options.add_time_flow_formulation = true;
    options.time_flow_state_disaggregated = true;
    options.solver_time_limit = 60.0;
    options.gurobi_threads = 4;
    const spdp::DirectTwoIndexResult result = spdp::solve_direct_two_index_model(data, graph, options);
    check(result.has_feasible_solution, "reference IP produced an integer solution");
    if (!result.has_feasible_solution) return;
    std::vector<double> y = result.edge_values;
    for (double& v : y) v = v > 0.5 ? 1.0 : 0.0;

    const spdp::CapacityBlossomSeparator separator = spdp::build_capacity_blossom_separator(graph);
    spdp::CapacityBlossomOptions blossom;
    blossom.enabled = true;
    blossom.min_violation = 1e-9;
    spdp::CapacityBlossomPool pool;
    spdp::CapacityBlossomStats stats;
    const auto cuts = spdp::separate_capacity_blossom_cuts(separator, y, blossom, 1000, pool, stats);
    check(cuts.empty(), "integer solution violates no capacity blossom cut");
    const double p_ref = spdp::max_capacity_blossom_violation_by_enumeration(separator, spdp::CapacityFamily::Pickup, y, 20);
    const double d_ref = spdp::max_capacity_blossom_violation_by_enumeration(separator, spdp::CapacityFamily::Delivery, y, 20);
    check(p_ref <= 1e-9 && d_ref <= 1e-9, "enumeration confirms no violated odd set on the integer solution");

    spdp::TypeTravelTimeCoverOptions cover;
    cover.enabled = true;
    cover.subproblem_time_limit = 20.0;
    cover.total_time_limit = 90.0;
    cover.gurobi_threads = 4;
    const spdp::TypeTravelTimeCoverResult bounds = spdp::compute_type_travel_time_cover_bounds(data, cover);
    const auto rows = spdp::build_type_travel_time_cover_rows(data, graph, bounds, 1);
    std::size_t certified = 0;
    for (const auto& e : bounds.entries) if (e.certified) ++certified;
    std::cout << "info: type cover entries " << bounds.entries.size() << ", certified " << certified
              << ", rows " << rows.size() << ", seconds " << bounds.total_seconds << '\n';
    spdp::write_type_travel_time_cover_log(std::cout, bounds);
    bool all_satisfied = true;
    for (const auto& row : rows) {
        double lhs = 0.0;
        for (const auto& term : row.row.edge_terms) lhs += term.second * y[term.first];
        if (lhs < row.row.rhs - 1e-9) {
            all_satisfied = false;
            std::cerr << "violated cover row mask " << row.type_mask << " rho " << row.rho << " lhs " << lhs << '\n';
        }
    }
    check(all_satisfied, "integer solution satisfies every certified type travel-time cover row");
    // Restricted instance sanity: single-type mask keeps exactly that type.
    const std::vector<int> types = spdp::instance_container_types(data);
    const spdp::ShortestTravelTimes s = spdp::compute_physical_shortest_travel_times(data);
    const auto closure = spdp::metric_closure_travel_time_matrix(data, s);
    const spdp::SPDPData restricted = spdp::build_type_restricted_instance(data, closure, types, 1ULL);
    bool single_type = !restricted.requests.empty();
    for (const auto& r : restricted.requests) single_type = single_type && r.container_type == types[0];
    check(single_type, "restricted instance for mask 1 contains only the first type");
}

}  // namespace

int main(int argc, char** argv) {
    const std::string instance = argc > 1 ? argv[1] : "RecDep_day_A8.dat";
    const bool quick = argc > 2 && std::string(argv[2]) == "quick";
    test_shortest_travel_times_small();
    const spdp::SPDPData data = spdp::read_spdp_data(instance);
    test_metric_closure_instance(data);
    std::ostringstream sink;
    const spdp::MultiDiGraph graph = spdp::build_multigraph(data, spdp::GraphBuildOptions{}, &sink);
    test_blossom_exactness(graph);
    if (!quick) {
        test_validity_on_integer_solution(data, graph);
    }
    if (g_failures != 0) {
        std::cerr << g_failures << " check(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "all checks passed\n";
    return EXIT_SUCCESS;
}

