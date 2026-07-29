#include "PstepValidInequality.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gurobi_c++.h"

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;
constexpr double kAuxiliaryLPSolverTolerance = 1e-9;

enum class LocationFamily {
    Pickup,
    Treatment,
    Delivery,
};

struct RequestServiceNodes {
    NodeId pickup_node = -1;
    NodeId delivery_node = -1;
};

void append_rows(
    std::vector<PstepValidInequalityRow>& destination,
    std::vector<PstepValidInequalityRow> source
) {
    destination.insert(
        destination.end(),
        std::make_move_iterator(source.begin()),
        std::make_move_iterator(source.end())
    );
}

void add_unit_edge_term(PstepValidInequalityRow& row, std::size_t edge_id) {
    row.edge_terms.emplace_back(edge_id, 1.0);
}

int cover_rhs(std::size_t request_count) {
    return static_cast<int>((request_count + 1U) / 2U);
}

bool edge_has_treatment(const EdgeRecord& edge) {
    return !edge.data.sequence_pi.empty();
}

std::vector<PstepValidInequalityRow> build_vi35_rows(
    const SPDPData& data,
    const MultiDiGraph& graph
) {
    const double rhs = static_cast<double>(cover_rhs(data.requests.size()));
    std::array<PstepValidInequalityRow, 6> rows;
    const std::array<std::string, 6> names{
        "vi35_pickup_in",
        "vi35_pickup_out",
        "vi35_treatment_in",
        "vi35_treatment_out",
        "vi35_delivery_in",
        "vi35_delivery_out",
    };
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        rows[row_index].name = names[row_index];
        rows[row_index].rhs = rhs;
    }

    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        const NodeSpec& u = graph.node(edge.u);
        const NodeSpec& v = graph.node(edge.v);
        const bool has_treatment = edge_has_treatment(edge);

        if (v.kind == NodeSpec::Kind::Pickup &&
            (has_treatment || u.kind == NodeSpec::Kind::Start ||
             u.kind == NodeSpec::Kind::Delivery)) {
            add_unit_edge_term(rows[0], edge_id);
        }
        if (u.kind == NodeSpec::Kind::Pickup &&
            (has_treatment || v.kind == NodeSpec::Kind::Delivery ||
             v.kind == NodeSpec::Kind::End)) {
            add_unit_edge_term(rows[1], edge_id);
        }
        if (has_treatment) {
            add_unit_edge_term(rows[2], edge_id);
            add_unit_edge_term(rows[3], edge_id);
        }
        if (v.kind == NodeSpec::Kind::Delivery &&
            (has_treatment || u.kind == NodeSpec::Kind::Pickup ||
             u.kind == NodeSpec::Kind::Start)) {
            add_unit_edge_term(rows[4], edge_id);
        }
        if (u.kind == NodeSpec::Kind::Delivery &&
            (has_treatment || v.kind == NodeSpec::Kind::Pickup ||
             v.kind == NodeSpec::Kind::End)) {
            add_unit_edge_term(rows[5], edge_id);
        }
    }

    return std::vector<PstepValidInequalityRow>(
        std::make_move_iterator(rows.begin()),
        std::make_move_iterator(rows.end())
    );
}

bool contains_location(const std::vector<int>& subset, int location) {
    return std::binary_search(subset.begin(), subset.end(), location);
}

void enumerate_fixed_size_location_subsets(
    const std::vector<int>& locations,
    std::size_t target_size,
    std::size_t next_index,
    std::vector<int>& current,
    std::vector<std::vector<int>>& subsets
) {
    if (current.size() == target_size) {
        subsets.push_back(current);
        return;
    }

    const std::size_t remaining = target_size - current.size();
    for (std::size_t index = next_index; index + remaining <= locations.size(); ++index) {
        current.push_back(locations[index]);
        enumerate_fixed_size_location_subsets(
            locations,
            target_size,
            index + 1U,
            current,
            subsets
        );
        current.pop_back();
    }
}

std::vector<std::vector<int>> enumerate_location_subsets(
    const std::map<int, int>& count_by_location,
    std::size_t max_subset_size,
    bool skip_full_location_set
) {
    std::vector<int> locations;
    locations.reserve(count_by_location.size());
    for (const auto& entry : count_by_location) {
        locations.push_back(entry.first);
    }

    const std::size_t capped_max_size = std::min(max_subset_size, locations.size());
    std::vector<std::vector<int>> subsets;
    std::vector<int> current;
    for (std::size_t subset_size = 1U; subset_size <= capped_max_size; ++subset_size) {
        if (skip_full_location_set && subset_size == locations.size()) {
            continue;
        }
        enumerate_fixed_size_location_subsets(
            locations,
            subset_size,
            0U,
            current,
            subsets
        );
    }
    return subsets;
}

std::string location_family_name(LocationFamily family) {
    switch (family) {
        case LocationFamily::Pickup:
            return "pickup";
        case LocationFamily::Treatment:
            return "treatment";
        case LocationFamily::Delivery:
            return "delivery";
    }
    throw std::runtime_error("Unknown VI-36 location family.");
}

std::string location_subset_name_suffix(const std::vector<int>& subset) {
    std::string suffix = subset.size() == 1U ? "loc" : "locs";
    for (int location : subset) {
        suffix += "_" + std::to_string(location);
    }
    return suffix;
}

int action_count_in_subset(
    const std::map<int, int>& count_by_location,
    const std::vector<int>& subset
) {
    int action_count = 0;
    for (int location : subset) {
        action_count += count_by_location.at(location);
    }
    return action_count;
}

std::pair<int, int> pickup_or_delivery_coefficients(
    const NodeSpec& u,
    const NodeSpec& v,
    bool has_treatment,
    const std::vector<int>& subset,
    NodeSpec::Kind service_kind
) {
    const bool u_in_subset =
        u.kind == service_kind && contains_location(subset, u.location);
    const bool v_in_subset =
        v.kind == service_kind && contains_location(subset, v.location);
    const int inbound = v_in_subset && (has_treatment || !u_in_subset) ? 1 : 0;
    const int outbound = u_in_subset && (has_treatment || !v_in_subset) ? 1 : 0;
    return {inbound, outbound};
}

std::pair<int, int> treatment_coefficients(
    const std::vector<int>& sequence_pi,
    const std::vector<int>& subset
) {
    int inbound = 0;
    int outbound = 0;
    for (std::size_t index = 0; index < sequence_pi.size(); ++index) {
        if (!contains_location(subset, sequence_pi[index])) {
            continue;
        }
        const bool previous_in_subset =
            index > 0U && contains_location(subset, sequence_pi[index - 1U]);
        const bool next_in_subset =
            index + 1U < sequence_pi.size() &&
            contains_location(subset, sequence_pi[index + 1U]);
        if (!previous_in_subset) {
            ++inbound;
        }
        if (!next_in_subset) {
            ++outbound;
        }
    }
    return {inbound, outbound};
}

std::vector<PstepValidInequalityRow> build_vi36_family_rows(
    const MultiDiGraph& graph,
    const std::map<int, int>& count_by_location,
    LocationFamily family,
    std::size_t max_subset_size,
    bool skip_full_location_set
) {
    const std::vector<std::vector<int>> subsets = enumerate_location_subsets(
        count_by_location,
        max_subset_size,
        skip_full_location_set
    );
    std::vector<PstepValidInequalityRow> rows;
    rows.reserve(subsets.size());

    for (const std::vector<int>& subset : subsets) {
        PstepValidInequalityRow row;
        row.name = "vi36_combined_" + location_family_name(family) + "_" +
                   location_subset_name_suffix(subset);
        const int action_count = action_count_in_subset(count_by_location, subset);
        row.rhs = 2.0 * static_cast<double>((action_count + 1) / 2);

        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            std::pair<int, int> coefficients;
            if (family == LocationFamily::Treatment) {
                coefficients = treatment_coefficients(edge.data.sequence_pi, subset);
            } else {
                const NodeSpec::Kind service_kind =
                    family == LocationFamily::Pickup
                        ? NodeSpec::Kind::Pickup
                        : NodeSpec::Kind::Delivery;
                coefficients = pickup_or_delivery_coefficients(
                    graph.node(edge.u),
                    graph.node(edge.v),
                    edge_has_treatment(edge),
                    subset,
                    service_kind
                );
            }

            const int combined_coefficient = coefficients.first + coefficients.second;
            if (combined_coefficient > 0) {
                row.edge_terms.emplace_back(
                    edge_id,
                    static_cast<double>(combined_coefficient)
                );
            }
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

std::vector<PstepValidInequalityRow> build_vi36_rows(
    const MultiDiGraph& graph,
    std::size_t max_subset_size,
    bool skip_full_location_sets
) {
    if (max_subset_size == 0U) {
        return {};
    }

    std::map<int, int> pickup_count_by_location;
    std::map<int, int> treatment_count_by_location;
    std::map<int, int> delivery_count_by_location;
    for (NodeId node_id = 1; node_id < graph.end_node_id(); ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const NodeSpec& node = graph.node(node_id);
        if (node.kind == NodeSpec::Kind::Pickup) {
            ++pickup_count_by_location[node.location];
            if (!node.landfill_location.has_value()) {
                throw std::runtime_error("Pickup node must have treatment location for VI-36.");
            }
            ++treatment_count_by_location[node.landfill_location.value()];
        } else if (node.kind == NodeSpec::Kind::Delivery) {
            ++delivery_count_by_location[node.location];
        }
    }

    std::vector<PstepValidInequalityRow> rows;
    append_rows(
        rows,
        build_vi36_family_rows(
            graph,
            pickup_count_by_location,
            LocationFamily::Pickup,
            max_subset_size,
            skip_full_location_sets
        )
    );
    append_rows(
        rows,
        build_vi36_family_rows(
            graph,
            treatment_count_by_location,
            LocationFamily::Treatment,
            max_subset_size,
            skip_full_location_sets
        )
    );
    append_rows(
        rows,
        build_vi36_family_rows(
            graph,
            delivery_count_by_location,
            LocationFamily::Delivery,
            max_subset_size,
            skip_full_location_sets
        )
    );
    return rows;
}

void enumerate_fixed_size_request_subsets(
    std::size_t request_count,
    std::size_t target_size,
    std::size_t next_index,
    std::vector<std::size_t>& current,
    std::vector<std::vector<std::size_t>>& subsets
) {
    if (current.size() == target_size) {
        subsets.push_back(current);
        return;
    }

    const std::size_t remaining = target_size - current.size();
    for (std::size_t index = next_index; index + remaining <= request_count; ++index) {
        current.push_back(index);
        enumerate_fixed_size_request_subsets(
            request_count,
            target_size,
            index + 1U,
            current,
            subsets
        );
        current.pop_back();
    }
}

std::string request_subset_name_suffix(
    const std::vector<std::size_t>& request_indices
) {
    std::string suffix;
    for (std::size_t request_index : request_indices) {
        suffix += "_" + std::to_string(request_index + 1U);
    }
    return suffix;
}

std::vector<PstepValidInequalityRow> build_request_block_sec_rows(
    const SPDPData& data,
    const MultiDiGraph& graph,
    std::size_t max_subset_size
) {
    if (max_subset_size == 0U) {
        return {};
    }
    if (max_subset_size < 2U) {
        throw std::runtime_error(
            "Request-block SEC max size must be zero or at least 2."
        );
    }

    const std::size_t request_count = data.requests.size();
    std::vector<RequestServiceNodes> service_nodes_by_request(request_count);
    for (NodeId node_id = 1; node_id < graph.end_node_id(); ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const NodeSpec& node = graph.node(node_id);
        if (!node.request_idx.has_value()) {
            throw std::runtime_error("Physical service node is missing its request index.");
        }
        const int request_index = node.request_idx.value();
        if (request_index < 0 || static_cast<std::size_t>(request_index) >= request_count) {
            throw std::runtime_error("Invalid request index on physical service node.");
        }
        RequestServiceNodes& service_nodes =
            service_nodes_by_request[static_cast<std::size_t>(request_index)];
        if (node.kind == NodeSpec::Kind::Pickup) {
            service_nodes.pickup_node = node_id;
        } else if (node.kind == NodeSpec::Kind::Delivery) {
            service_nodes.delivery_node = node_id;
        }
    }

    for (const RequestServiceNodes& service_nodes : service_nodes_by_request) {
        if (service_nodes.pickup_node < 0 || service_nodes.delivery_node < 0) {
            throw std::runtime_error("Each request must have one pickup and one delivery node.");
        }
    }

    const std::size_t capped_max_size = std::min(max_subset_size, request_count);
    std::vector<std::vector<std::size_t>> request_subsets;
    std::vector<std::size_t> current;
    for (std::size_t subset_size = 2U; subset_size <= capped_max_size; ++subset_size) {
        enumerate_fixed_size_request_subsets(
            request_count,
            subset_size,
            0U,
            current,
            request_subsets
        );
    }

    std::vector<PstepValidInequalityRow> rows;
    rows.reserve(request_subsets.size());
    for (const std::vector<std::size_t>& request_indices : request_subsets) {
        std::vector<std::uint8_t> in_block(graph.number_of_nodes(), 0U);
        for (std::size_t request_index : request_indices) {
            const RequestServiceNodes& service_nodes = service_nodes_by_request[request_index];
            in_block[static_cast<std::size_t>(service_nodes.pickup_node)] = 1U;
            in_block[static_cast<std::size_t>(service_nodes.delivery_node)] = 1U;
        }

        PstepValidInequalityRow row;
        row.name = "vi_request_block_sec_req" +
                   request_subset_name_suffix(request_indices);
        row.sense = PstepValidInequalitySense::LessEqual;
        row.rhs = static_cast<double>(2U * request_indices.size() - 1U);
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            if (in_block[static_cast<std::size_t>(edge.u)] != 0U &&
                in_block[static_cast<std::size_t>(edge.v)] != 0U) {
                add_unit_edge_term(row, edge_id);
            }
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

int cor_route_lower_bound_kmin(const SPDPData& data) {
    if (data.time_limit <= 0.0) {
        throw std::runtime_error("Time limit must be positive for VI-44.");
    }

    double pickup_to_treatment_time_sum = 0.0;
    double treatment_to_best_delivery_time_sum = 0.0;
    for (const Request& request : data.requests) {
        pickup_to_treatment_time_sum +=
            data.time[static_cast<std::size_t>(request.from_id)]
                     [static_cast<std::size_t>(request.to_id)];

        double best_compatible_delivery_time = std::numeric_limits<double>::infinity();
        for (const Request& candidate : data.requests) {
            if (candidate.container_type != request.container_type) {
                continue;
            }
            best_compatible_delivery_time = std::min(
                best_compatible_delivery_time,
                data.time[static_cast<std::size_t>(request.to_id)]
                         [static_cast<std::size_t>(candidate.from_id)]
            );
        }
        if (!std::isfinite(best_compatible_delivery_time)) {
            throw std::runtime_error(
                "Failed to compute VI-44 k_min because no compatible delivery location exists."
            );
        }
        treatment_to_best_delivery_time_sum += best_compatible_delivery_time;
    }

    const double total_service_time =
        static_cast<double>(data.requests.size()) *
        (data.time_pickup + data.time_empty + data.time_delivery);
    const double lower_bound =
        0.5 * (pickup_to_treatment_time_sum + treatment_to_best_delivery_time_sum) +
        total_service_time;
    return std::max(
        0,
        static_cast<int>(std::ceil(lower_bound / data.time_limit - kTolerance))
    );
}

State canonical_auxiliary_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

struct AuxiliaryDurationSubproblemResult {
    int status = 0;
    bool hit_time_limit = false;
    bool has_certified_bound = false;
    double objective_value = -1.0;
    double objective_bound = -1.0;
    double runtime_seconds = 0.0;
};

AuxiliaryDurationSubproblemResult solve_auxiliary_duration_subproblem(
    const SPDPData& data,
    const MultiDiGraph& graph,
    double solver_time_limit,
    bool binary_y,
    bool add_time_constraints
) {
    if (data.time_limit <= 0.0) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem requires a positive route-duration limit."
        );
    }
    if (!std::isfinite(solver_time_limit) || solver_time_limit < 0.0) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem time limit must be finite and nonnegative."
        );
    }
    if (graph.number_of_nodes() != 2U * data.requests.size() + 2U) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem requires exactly two action nodes per request."
        );
    }

    const NodeId start_node_id = 0;
    const NodeId end_node_id = graph.end_node_id();
    if (graph.node(start_node_id).kind != NodeSpec::Kind::Start ||
        graph.node(end_node_id).kind != NodeSpec::Kind::End) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem requires node 0 and the last node to be the depots."
        );
    }

    std::size_t dummy_edge_count = 0U;
    for (const EdgeRecord& edge : graph.edges()) {
        if (!std::isfinite(edge.data.time) || edge.data.time < -kTolerance) {
            throw std::runtime_error(
                "The auxiliary VI-44 subproblem requires finite, nonnegative edge durations."
            );
        }
        if (edge.u == start_node_id && edge.v == end_node_id) {
            ++dummy_edge_count;
            if (std::abs(edge.data.time) > kTolerance) {
                throw std::runtime_error(
                    "The auxiliary VI-44 subproblem requires zero-duration dummy depot edges."
                );
            }
        }
        if (edge.u == start_node_id && edge.v != end_node_id &&
            graph.node(edge.v).kind != NodeSpec::Kind::Pickup) {
            throw std::runtime_error(
                "An auxiliary VI-44 subproblem departure edge must enter a pickup node."
            );
        }
        if (edge.v == end_node_id && edge.u != start_node_id &&
            graph.node(edge.u).kind != NodeSpec::Kind::Delivery) {
            throw std::runtime_error(
                "An auxiliary VI-44 subproblem return edge must leave a delivery node."
            );
        }
    }
    if (dummy_edge_count == 0U) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem graph has no start-to-end dummy edge."
        );
    }

    GRBEnv environment(true);
    environment.set(GRB_IntParam_OutputFlag, 0);
    environment.start();
    GRBModel model(environment);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_Threads, 1);
    model.set(GRB_DoubleParam_FeasibilityTol, kAuxiliaryLPSolverTolerance);
    model.set(GRB_DoubleParam_OptimalityTol, kAuxiliaryLPSolverTolerance);
    if (solver_time_limit > 0.0) {
        model.set(GRB_DoubleParam_TimeLimit, solver_time_limit);
    }

    std::vector<GRBVar> y_vars;
    y_vars.reserve(graph.number_of_edges());
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        const bool is_dummy = edge.u == start_node_id && edge.v == end_node_id;
        y_vars.push_back(model.addVar(
            0.0,
            is_dummy ? 0.0 : 1.0,
            edge.data.time,
            binary_y ? GRB_BINARY : GRB_CONTINUOUS,
            "vi44_aux_y_" + std::to_string(edge_id)
        ));
    }
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.update();

    std::vector<std::set<State>> states_by_node(graph.number_of_nodes());
    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const NodeSpec::Kind node_kind = graph.node(node_id).kind;
        if (node_kind != NodeSpec::Kind::Pickup && node_kind != NodeSpec::Kind::Delivery) {
            throw std::runtime_error(
                "The auxiliary VI-44 subproblem found a non-service physical node."
            );
        }

        GRBLinExpr incoming_degree = 0.0;
        for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
            incoming_degree += y_vars[edge_id];
        }
        model.addConstr(
            incoming_degree == 1.0,
            "vi44_aux_in_degree_" + std::to_string(node_id)
        );

        GRBLinExpr outgoing_degree = 0.0;
        for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
            outgoing_degree += y_vars[edge_id];
        }
        model.addConstr(
            outgoing_degree == 1.0,
            "vi44_aux_out_degree_" + std::to_string(node_id)
        );

        for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
            states_by_node[static_cast<std::size_t>(node_id)].insert(
                canonical_auxiliary_state(graph.edges()[edge_id].data.end_state)
            );
        }
        for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
            states_by_node[static_cast<std::size_t>(node_id)].insert(
                canonical_auxiliary_state(graph.edges()[edge_id].data.start_state)
            );
        }

        std::size_t state_index = 0U;
        for (const State& state : states_by_node[static_cast<std::size_t>(node_id)]) {
            GRBLinExpr state_balance = 0.0;
            for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
                if (canonical_auxiliary_state(graph.edges()[edge_id].data.end_state) == state) {
                    state_balance += y_vars[edge_id];
                }
            }
            for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
                if (canonical_auxiliary_state(graph.edges()[edge_id].data.start_state) == state) {
                    state_balance -= y_vars[edge_id];
                }
            }
            model.addConstr(
                state_balance == 0.0,
                "vi44_aux_state_" + std::to_string(node_id) + "_" +
                    std::to_string(state_index++)
            );
        }
    }

    if (add_time_constraints) {
        std::vector<std::map<State, GRBVar>> time_vars_by_node(
            graph.number_of_nodes()
        );
        for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
            std::size_t state_index = 0U;
            for (const State& state :
                 states_by_node[static_cast<std::size_t>(node_id)]) {
                time_vars_by_node[static_cast<std::size_t>(node_id)].emplace(
                    state,
                    model.addVar(
                        0.0,
                        data.time_limit,
                        0.0,
                        GRB_CONTINUOUS,
                        "vi44_aux_B_" + std::to_string(node_id) + "_" +
                            std::to_string(state_index++)
                    )
                );
            }
        }
        model.update();

        const auto time_var = [&](NodeId node_id, const State& state) -> GRBVar {
            const State canonical_state = canonical_auxiliary_state(state);
            const auto& variables =
                time_vars_by_node.at(static_cast<std::size_t>(node_id));
            const auto found = variables.find(canonical_state);
            if (found == variables.end()) {
                throw std::runtime_error(
                    "The auxiliary VI-44 subproblem found no time variable for an "
                    "edge endpoint state."
                );
            }
            return found->second;
        };

        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            if (edge.u == start_node_id && edge.v == end_node_id) {
                continue;
            }

            const std::string constraint_name =
                "vi44_aux_time_" + std::to_string(edge_id);
            if (edge.u == start_node_id) {
                const GRBVar end_time =
                    time_var(edge.v, edge.data.end_state);
                model.addConstr(
                    end_time >= edge.data.time * y_vars[edge_id],
                    constraint_name
                );
                continue;
            }
            if (edge.v == end_node_id) {
                const GRBVar start_time =
                    time_var(edge.u, edge.data.start_state);
                model.addConstr(
                    start_time + edge.data.time * y_vars[edge_id] <=
                        data.time_limit,
                    constraint_name
                );
                continue;
            }

            const GRBVar start_time =
                time_var(edge.u, edge.data.start_state);
            const GRBVar end_time =
                time_var(edge.v, edge.data.end_state);
            const double big_m = data.time_limit + edge.data.time;
            model.addConstr(
                end_time >=
                    start_time + edge.data.time -
                    big_m * (1.0 - y_vars[edge_id]),
                constraint_name
            );
        }
    }

    GRBLinExpr actual_departures = 0.0;
    GRBLinExpr actual_returns = 0.0;
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == start_node_id &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            actual_departures += y_vars[edge_id];
        }
        if (graph.node(edge.u).kind == NodeSpec::Kind::Delivery &&
            edge.v == end_node_id) {
            actual_returns += y_vars[edge_id];
        }
    }
    model.addConstr(actual_departures == actual_returns, "vi44_aux_depot_balance");

    const auto solve_start = std::chrono::steady_clock::now();
    model.optimize();
    const auto solve_end = std::chrono::steady_clock::now();

    AuxiliaryDurationSubproblemResult result;
    result.status = model.get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.runtime_seconds =
        std::chrono::duration<double>(solve_end - solve_start).count();

    if (model.get(GRB_IntAttr_SolCount) > 0) {
        try {
            const double objective_value = model.get(GRB_DoubleAttr_ObjVal);
            if (std::isfinite(objective_value) &&
                std::abs(objective_value) < 0.5 * GRB_INFINITY) {
                result.objective_value = objective_value;
            }
        } catch (const GRBException& error) {
            if (error.getErrorCode() != GRB_ERROR_DATA_NOT_AVAILABLE) {
                throw;
            }
        }
    }

    try {
        const double objective_bound = model.get(GRB_DoubleAttr_ObjBound);
        if (std::isfinite(objective_bound) &&
            std::abs(objective_bound) < 0.5 * GRB_INFINITY) {
            result.objective_bound = objective_bound;
            result.has_certified_bound = true;
        }
    } catch (const GRBException& error) {
        if (error.getErrorCode() != GRB_ERROR_DATA_NOT_AVAILABLE) {
            throw;
        }
    }

    if (result.status != GRB_OPTIMAL && result.status != GRB_TIME_LIMIT) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem stopped with an unsupported status (" +
            std::to_string(result.status) + ")."
        );
    }
    if (result.status == GRB_OPTIMAL && !result.has_certified_bound) {
        throw std::runtime_error(
            "The optimal auxiliary VI-44 subproblem returned no finite certified bound."
        );
    }
    return result;
}

const char* vi44_k_min_mode_name(VI44KMinMode mode) {
    switch (mode) {
        case VI44KMinMode::COR:
            return "cor";
        case VI44KMinMode::SubLP:
            return "sub-lp";
        case VI44KMinMode::SubIP:
            return "sub-ip";
    }
    throw std::runtime_error("Unknown VI-44 k_min mode.");
}

PstepValidInequalityRow build_vi44_row(
    const MultiDiGraph& graph,
    int k_min
) {
    PstepValidInequalityRow row;
    row.name = "vi44_route_lower_bound";
    row.rhs = static_cast<double>(k_min);
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == 0 && graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            add_unit_edge_term(row, edge_id);
        }
    }
    return row;
}

}  // namespace

VI44KMinResult compute_vi44_k_min(
    const SPDPData& data,
    const MultiDiGraph& graph,
    VI44KMinMode mode,
    double sub_lp_time_limit,
    bool add_time_constraints
) {
    if (!std::isfinite(sub_lp_time_limit) || sub_lp_time_limit < 0.0) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem time limit must be finite and nonnegative."
        );
    }
    VI44KMinResult result;
    result.cor_k_min = cor_route_lower_bound_kmin(data);
    result.selected_k_min = result.cor_k_min;
    if (mode == VI44KMinMode::COR) {
        return result;
    }
    if (mode != VI44KMinMode::SubLP && mode != VI44KMinMode::SubIP) {
        throw std::runtime_error("Unknown VI-44 k_min mode.");
    }

    const AuxiliaryDurationSubproblemResult auxiliary =
        solve_auxiliary_duration_subproblem(
            data,
            graph,
            sub_lp_time_limit,
            mode == VI44KMinMode::SubIP,
            add_time_constraints
        );
    result.sub_lp_status = auxiliary.status;
    result.sub_lp_hit_time_limit = auxiliary.hit_time_limit;
    result.sub_lp_has_certified_bound = auxiliary.has_certified_bound;
    result.sub_lp_objective_value = auxiliary.objective_value;
    result.sub_lp_objective_bound = auxiliary.objective_bound;
    result.sub_lp_runtime_seconds = auxiliary.runtime_seconds;

    if (!auxiliary.has_certified_bound) {
        return result;
    }

    const double scale = std::max(
        {1.0, data.time_limit, std::abs(auxiliary.objective_bound)}
    );
    result.sub_lp_numerical_tolerance =
        10.0 * kAuxiliaryLPSolverTolerance * scale;
    result.sub_lp_safe_lower_bound = std::max(
        0.0,
        auxiliary.objective_bound - result.sub_lp_numerical_tolerance
    );
    result.sub_lp_k_min = std::max(
        0,
        static_cast<int>(std::ceil(result.sub_lp_safe_lower_bound / data.time_limit))
    );
    result.selected_k_min = std::max(result.cor_k_min, result.sub_lp_k_min);
    return result;
}

std::vector<PstepValidInequalityRow> build_pstep_valid_inequality_rows(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const PstepValidInequalityOptions& options
) {
    std::vector<PstepValidInequalityRow> rows;
    if (options.add_vi_35) {
        append_rows(rows, build_vi35_rows(data, graph));
    }
    if (options.add_vi_36_combined && options.vi_36_subset_max_size > 0U) {
        append_rows(
            rows,
            build_vi36_rows(
                graph,
                options.vi_36_subset_max_size,
                options.add_vi_35
            )
        );
    }
    if (options.add_vi_request_block_sec &&
        options.vi_request_block_sec_max_size > 0U) {
        append_rows(
            rows,
            build_request_block_sec_rows(
                data,
                graph,
                options.vi_request_block_sec_max_size
            )
        );
    }
    if (options.add_vi_44) {
        const VI44KMinResult result =
            compute_vi44_k_min(
                data,
                graph,
                options.vi_44_k_min_mode,
                options.vi_44_sub_lp_time_limit,
                options.vi_44_subproblem_add_time_constraints
            );
        rows.push_back(build_vi44_row(graph, result.selected_k_min));
        if (options.log_stream != nullptr) {
            std::ostringstream message;
            message << std::setprecision(15)
                    << "[pstep-vi] VI44 k_min mode="
                    << vi44_k_min_mode_name(options.vi_44_k_min_mode)
                    << " cor=" << result.cor_k_min;
            if (options.vi_44_k_min_mode != VI44KMinMode::COR) {
                message << " sub_lp_status=" << result.sub_lp_status
                        << " sub_lp_time_limit=" << options.vi_44_sub_lp_time_limit
                        << " subproblem_time_constraints="
                        << (options.vi_44_subproblem_add_time_constraints ? 1 : 0)
                        << " sub_lp_hit_time_limit="
                        << (result.sub_lp_hit_time_limit ? 1 : 0)
                        << " sub_lp_has_certified_bound="
                        << (result.sub_lp_has_certified_bound ? 1 : 0)
                        << " sub_lp_objective=" << result.sub_lp_objective_value
                        << " sub_lp_bound=" << result.sub_lp_objective_bound
                        << " sub_lp_safe_bound=" << result.sub_lp_safe_lower_bound
                        << " sub_lp_tolerance=" << result.sub_lp_numerical_tolerance
                        << " sub_lp_k_min=" << result.sub_lp_k_min
                        << " sub_lp_runtime_seconds=" << result.sub_lp_runtime_seconds;
            }
            message << " selected=" << result.selected_k_min;
            *options.log_stream << message.str() << '\n';
        }
    }
    return rows;
}

}  // namespace spdp
