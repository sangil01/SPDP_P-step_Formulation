#include "CpIncumbentAdapter.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "PstepFormulation.h"
#include "RouteState.h"

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;

bool same_state(const State& left, const State& right) {
    return left == right;
}

int service_node_id(CpActionKind kind, int request, int n) {
    if (kind == CpActionKind::Pickup) {
        return request + 1;
    }
    if (kind == CpActionKind::Delivery) {
        return n + request + 1;
    }
    return -1;
}

std::vector<CpActionVisit> normalize_treatments(
    const SPDPData& data,
    const CpActionRoute& route
) {
    const int n = static_cast<int>(data.requests.size());
    std::vector<bool> handled(static_cast<std::size_t>(n), false);
    std::set<int> expected_treatments;
    for (const CpActionVisit& visit : route.actions) {
        if (visit.kind == CpActionKind::Treatment) {
            if (!expected_treatments.insert(visit.request_index).second) {
                throw std::runtime_error(
                    "CP route repeats a treatment action."
                );
            }
        }
    }
    std::vector<CpActionVisit> normalized;
    normalized.reserve(route.actions.size());
    RouteState state;

    for (const CpActionVisit& visit : route.actions) {
        if (visit.request_index < 0 || visit.request_index >= n) {
            throw std::runtime_error("CP route contains an invalid request index.");
        }
        const Request& request =
            data.requests[static_cast<std::size_t>(visit.request_index)];
        if (visit.kind == CpActionKind::Pickup) {
            if (!state.try_pickup(OnboardSkip{
                    visit.request_index,
                    request.container_type,
                    request.to_id,
                    true
                })) {
                throw std::runtime_error(
                    "CP route violates capacity during treatment normalization."
                );
            }
            normalized.push_back(visit);
            continue;
        }
        if (visit.kind == CpActionKind::Delivery) {
            if (!state.try_delivery(request.container_type)) {
                throw std::runtime_error(
                    "CP route delivers without a compatible empty skip."
                );
            }
            normalized.push_back(visit);
            continue;
        }
        if (handled[static_cast<std::size_t>(visit.request_index)]) {
            continue;
        }

        std::vector<int> block;
        for (const OnboardSkip& skip : state.onboard()) {
            if (skip.is_full && skip.treatment_location == request.to_id &&
                !handled[static_cast<std::size_t>(skip.request_index)]) {
                block.push_back(skip.request_index);
            }
        }
        if (std::find(block.begin(), block.end(), visit.request_index) == block.end()) {
            throw std::runtime_error(
                "CP treatment action does not match a full onboard skip."
            );
        }
        std::sort(block.begin(), block.end());
        for (int request_index : block) {
            handled[static_cast<std::size_t>(request_index)] = true;
            normalized.push_back(CpActionVisit{
                CpActionKind::Treatment,
                request_index,
                -1,
                -1,
            });
        }
        if (state.empty_at_treatment(request.to_id) !=
            static_cast<int>(block.size())) {
            throw std::runtime_error(
                "Treatment normalization failed to batch all eligible skips."
            );
        }
    }

    const bool all_expected_handled = std::all_of(
        expected_treatments.begin(),
        expected_treatments.end(),
        [&handled](int request_index) {
            return request_index >= 0 &&
                static_cast<std::size_t>(request_index) < handled.size() &&
                handled[static_cast<std::size_t>(request_index)];
        }
    );
    if (!state.is_empty() || !all_expected_handled ||
        std::count(handled.begin(), handled.end(), true) !=
            static_cast<int>(expected_treatments.size())) {
        throw std::runtime_error(
            "Normalized CP route does not end empty or omits a treatment action."
        );
    }
    return normalized;
}

double normalized_transition_time(
    const SPDPData& data,
    int from_location,
    const std::vector<int>& treatments,
    int treatment_action_count,
    int to_location,
    double destination_service
) {
    double travel = 0.0;
    int current = from_location;
    for (int treatment : treatments) {
        if (current >= 0) {
            travel += data.time[static_cast<std::size_t>(current)]
                               [static_cast<std::size_t>(treatment)];
        }
        current = treatment;
    }
    if (current >= 0 && to_location >= 0) {
        travel += data.time[static_cast<std::size_t>(current)]
                           [static_cast<std::size_t>(to_location)];
    }
    return travel + treatment_action_count * data.time_empty +
        destination_service;
}

int find_edge(
    const MultiDiGraph& graph,
    NodeId u,
    NodeId v,
    const State& start_state,
    const State& end_state,
    const std::vector<int>& sequence,
    double normalized_time
) {
    int exact = -1;
    int fallback = -1;
    double fallback_time = std::numeric_limits<double>::infinity();
    for (std::size_t edge_id : graph.outgoing_edge_indices(u)) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.v != v || !same_state(edge.data.start_state, start_state) ||
            !same_state(edge.data.end_state, end_state)) {
            continue;
        }
        if (edge.data.sequence_pi == sequence) {
            if (exact < 0 || static_cast<int>(edge_id) < exact) {
                exact = static_cast<int>(edge_id);
            }
            continue;
        }
        if (edge.data.time <= normalized_time + kTolerance &&
            (edge.data.time < fallback_time - kTolerance ||
             (std::fabs(edge.data.time - fallback_time) <= kTolerance &&
              (fallback < 0 || static_cast<int>(edge_id) < fallback)))) {
            fallback = static_cast<int>(edge_id);
            fallback_time = edge.data.time;
        }
    }
    return exact >= 0 ? exact : fallback;
}

}  // namespace

CpMappedIncumbent map_cp_incumbent_to_multigraph(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<CpActionRoute>& routes
) {
    CpMappedIncumbent mapped;
    try {
        const int n = static_cast<int>(data.requests.size());
        std::set<int> globally_seen_service_nodes;
        for (const CpActionRoute& route : routes) {
            const std::vector<CpActionVisit> normalized =
                normalize_treatments(data, route);
            RouteState state;
            NodeId current_node = 0;
            int current_location = -1;
            State transition_start_state = state.canonical_state();
            std::vector<int> treatments;
            int treatment_action_count = 0;

            auto flush_to = [&](NodeId next_node,
                                int next_location,
                                double destination_service,
                                const CpActionVisit* destination) {
                const State start_state = transition_start_state;
                for (int treatment : treatments) {
                    state.empty_at_treatment(treatment);
                }
                if (destination != nullptr) {
                    const Request& request = data.requests[static_cast<std::size_t>(
                        destination->request_index)];
                    if (destination->kind == CpActionKind::Pickup) {
                        if (!state.try_pickup(OnboardSkip{
                                destination->request_index,
                                request.container_type,
                                request.to_id,
                                true
                            })) {
                            throw std::runtime_error(
                                "Normalized route exceeds capacity at pickup."
                            );
                        }
                    } else if (!state.try_delivery(request.container_type)) {
                        throw std::runtime_error(
                            "Normalized route has an incompatible delivery."
                        );
                    }
                }
                const State end_state = state.canonical_state();
                const double transition_time = normalized_transition_time(
                    data,
                    current_location,
                    treatments,
                    treatment_action_count,
                    next_location,
                    destination_service
                );
                const int edge_id = find_edge(
                    graph, current_node, next_node, start_state, end_state,
                    treatments, transition_time
                );
                if (edge_id < 0) {
                    std::ostringstream error;
                    error << "No surviving multigraph edge maps normalized transition "
                          << current_node << " -> " << next_node << '.';
                    throw std::runtime_error(error.str());
                }
                mapped.active_edge_ids.push_back(edge_id);
                mapped.total_duration += graph.edges()[static_cast<std::size_t>(edge_id)].data.time;
                mapped.total_original_cost += graph.edges()[static_cast<std::size_t>(edge_id)].data.cost;
                current_node = next_node;
                current_location = next_location;
                transition_start_state = end_state;
                treatments.clear();
                treatment_action_count = 0;
            };

            int last_treatment = -1;
            for (const CpActionVisit& visit : normalized) {
                const Request& request =
                    data.requests[static_cast<std::size_t>(visit.request_index)];
                if (visit.kind == CpActionKind::Treatment) {
                    if (request.to_id != last_treatment) {
                        treatments.push_back(request.to_id);
                        last_treatment = request.to_id;
                    }
                    ++treatment_action_count;
                    continue;
                }
                last_treatment = -1;
                const NodeId next_node = service_node_id(
                    visit.kind, visit.request_index, n
                );
                if (!globally_seen_service_nodes.insert(next_node).second) {
                    throw std::runtime_error(
                        "CP incumbent repeats a pickup or delivery service node."
                    );
                }
                flush_to(
                    next_node,
                    request.from_id,
                    visit.kind == CpActionKind::Pickup
                        ? data.time_pickup : data.time_delivery,
                    &visit
                );
            }
            if (!treatments.empty() || treatment_action_count != 0) {
                throw std::runtime_error(
                    "Normalized CP route ends with treatment actions."
                );
            }
            flush_to(graph.end_node_id(), -1, 0.0, nullptr);
            if (!state.is_empty()) {
                throw std::runtime_error("Mapped CP route does not end empty.");
            }
        }
        if (globally_seen_service_nodes.size() != 2U * data.requests.size()) {
            throw std::runtime_error(
                "CP incumbent does not cover every pickup and delivery exactly once."
            );
        }

        const RecoveredSolution recovered = recover_selected_edge_solution(
            data,
            graph,
            mapped.active_edge_ids,
            mapped.total_original_cost,
            0.0
        );
        if (recovered.routes.size() != routes.size()) {
            throw std::runtime_error(
                "Existing route validator recovered a different vehicle count."
            );
        }
        mapped.success = true;
    } catch (const std::exception& error) {
        mapped.success = false;
        mapped.error_message = error.what();
        mapped.active_edge_ids.clear();
        mapped.total_duration = 0.0;
        mapped.total_original_cost = 0.0;
    }
    return mapped;
}

}  // namespace spdp
