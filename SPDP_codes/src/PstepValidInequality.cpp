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
#include <memory>
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

enum class TwoIndexCoreObjective {
    OriginalCost,
    Duration,
};

struct TwoIndexCoreOptions {
    TwoIndexCoreObjective objective = TwoIndexCoreObjective::Duration;
    bool binary_y = false;
    bool add_time_constraints = false;
    double solver_time_limit = 0.0;
    int gurobi_threads = 1;
    bool output_enabled = false;
    std::string gurobi_log_path;
    std::string name_prefix = "two_index";
};

struct TwoIndexCoreModel {
    std::unique_ptr<GRBEnv> environment;
    std::unique_ptr<GRBModel> model;
    std::vector<GRBVar> y_vars;
};

void validate_two_index_inputs(
    const SPDPData& data,
    const MultiDiGraph& graph,
    double solver_time_limit
) {
    if (data.time_limit <= 0.0) {
        throw std::runtime_error(
            "The two-index formulation requires a positive route-duration limit."
        );
    }
    if (!std::isfinite(solver_time_limit) || solver_time_limit < 0.0) {
        throw std::runtime_error(
            "The two-index formulation time limit must be finite and nonnegative."
        );
    }
    if (graph.number_of_nodes() != 2U * data.requests.size() + 2U) {
        throw std::runtime_error(
            "The two-index formulation requires exactly two action nodes per request."
        );
    }

    const NodeId start_node_id = 0;
    const NodeId end_node_id = graph.end_node_id();
    if (graph.node(start_node_id).kind != NodeSpec::Kind::Start ||
        graph.node(end_node_id).kind != NodeSpec::Kind::End) {
        throw std::runtime_error(
            "The two-index formulation requires node 0 and the last node to be the depots."
        );
    }

    std::size_t dummy_edge_count = 0U;
    for (const EdgeRecord& edge : graph.edges()) {
        if (!std::isfinite(edge.data.time) || edge.data.time < -kTolerance ||
            !std::isfinite(edge.data.cost)) {
            throw std::runtime_error(
                "The two-index formulation requires finite costs and nonnegative durations."
            );
        }
        if (edge.u == start_node_id && edge.v == end_node_id) {
            ++dummy_edge_count;
            if (std::abs(edge.data.time) > kTolerance) {
                throw std::runtime_error(
                    "The two-index formulation requires zero-duration dummy depot edges."
                );
            }
        }
        if (edge.u == start_node_id && edge.v != end_node_id &&
            graph.node(edge.v).kind != NodeSpec::Kind::Pickup) {
            throw std::runtime_error(
                "A two-index departure edge must enter a pickup node."
            );
        }
        if (edge.v == end_node_id && edge.u != start_node_id &&
            graph.node(edge.u).kind != NodeSpec::Kind::Delivery) {
            throw std::runtime_error(
                "A two-index return edge must leave a delivery node."
            );
        }
    }
    if (dummy_edge_count == 0U) {
        throw std::runtime_error(
            "The two-index formulation graph has no start-to-end dummy edge."
        );
    }
}

TwoIndexCoreModel build_two_index_model_core(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const TwoIndexCoreOptions& options
) {
    validate_two_index_inputs(data, graph, options.solver_time_limit);

    TwoIndexCoreModel core;
    core.environment = std::make_unique<GRBEnv>(true);
    core.environment->set(
        GRB_IntParam_OutputFlag,
        options.output_enabled ? 1 : 0
    );
    if (!options.gurobi_log_path.empty()) {
        core.environment->set(GRB_StringParam_LogFile, options.gurobi_log_path);
    }
    core.environment->start();
    core.model = std::make_unique<GRBModel>(*core.environment);
    core.model->set(
        GRB_IntParam_OutputFlag,
        options.output_enabled ? 1 : 0
    );
    if (options.gurobi_threads >= 0) {
        core.model->set(GRB_IntParam_Threads, options.gurobi_threads);
    }
    core.model->set(GRB_DoubleParam_FeasibilityTol, kAuxiliaryLPSolverTolerance);
    core.model->set(GRB_DoubleParam_OptimalityTol, kAuxiliaryLPSolverTolerance);
    if (options.solver_time_limit > 0.0) {
        core.model->set(GRB_DoubleParam_TimeLimit, options.solver_time_limit);
    }

    const NodeId start_node_id = 0;
    const NodeId end_node_id = graph.end_node_id();
    core.y_vars.reserve(graph.number_of_edges());
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        const bool is_dummy = edge.u == start_node_id && edge.v == end_node_id;
        const double objective_coefficient =
            options.objective == TwoIndexCoreObjective::Duration
                ? edge.data.time
                : edge.data.cost;
        core.y_vars.push_back(core.model->addVar(
            0.0,
            is_dummy ? 0.0 : 1.0,
            objective_coefficient,
            options.binary_y ? GRB_BINARY : GRB_CONTINUOUS,
            options.name_prefix + "_y_" + std::to_string(edge_id)
        ));
    }
    core.model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    core.model->update();

    std::vector<std::set<State>> states_by_node(graph.number_of_nodes());
    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const NodeSpec::Kind node_kind = graph.node(node_id).kind;
        if (node_kind != NodeSpec::Kind::Pickup &&
            node_kind != NodeSpec::Kind::Delivery) {
            throw std::runtime_error(
                "The two-index formulation found a non-service physical node."
            );
        }

        GRBLinExpr incoming_degree = 0.0;
        for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
            incoming_degree += core.y_vars[edge_id];
        }
        core.model->addConstr(
            incoming_degree == 1.0,
            options.name_prefix + "_in_degree_" + std::to_string(node_id)
        );

        GRBLinExpr outgoing_degree = 0.0;
        for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
            outgoing_degree += core.y_vars[edge_id];
        }
        core.model->addConstr(
            outgoing_degree == 1.0,
            options.name_prefix + "_out_degree_" + std::to_string(node_id)
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
                if (canonical_auxiliary_state(graph.edges()[edge_id].data.end_state) ==
                    state) {
                    state_balance += core.y_vars[edge_id];
                }
            }
            for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
                if (canonical_auxiliary_state(graph.edges()[edge_id].data.start_state) ==
                    state) {
                    state_balance -= core.y_vars[edge_id];
                }
            }
            core.model->addConstr(
                state_balance == 0.0,
                options.name_prefix + "_state_" + std::to_string(node_id) + "_" +
                    std::to_string(state_index++)
            );
        }
    }

    if (options.add_time_constraints) {
        std::vector<std::map<State, GRBVar>> time_vars_by_node(
            graph.number_of_nodes()
        );
        for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
            std::size_t state_index = 0U;
            for (const State& state :
                 states_by_node[static_cast<std::size_t>(node_id)]) {
                time_vars_by_node[static_cast<std::size_t>(node_id)].emplace(
                    state,
                    core.model->addVar(
                        0.0,
                        data.time_limit,
                        0.0,
                        GRB_CONTINUOUS,
                        options.name_prefix + "_B_" + std::to_string(node_id) + "_" +
                            std::to_string(state_index++)
                    )
                );
            }
        }
        core.model->update();

        const auto time_var = [&](NodeId node_id, const State& state) -> GRBVar {
            const State canonical_state = canonical_auxiliary_state(state);
            const auto& variables =
                time_vars_by_node.at(static_cast<std::size_t>(node_id));
            const auto found = variables.find(canonical_state);
            if (found == variables.end()) {
                throw std::runtime_error(
                    "The two-index formulation found no time variable for an edge "
                    "endpoint state."
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
                options.name_prefix + "_time_" + std::to_string(edge_id);
            if (edge.u == start_node_id) {
                const GRBVar end_time = time_var(edge.v, edge.data.end_state);
                core.model->addConstr(
                    end_time >= edge.data.time * core.y_vars[edge_id],
                    constraint_name
                );
                continue;
            }
            if (edge.v == end_node_id) {
                const GRBVar start_time = time_var(edge.u, edge.data.start_state);
                core.model->addConstr(
                    start_time + edge.data.time * core.y_vars[edge_id] <=
                        data.time_limit,
                    constraint_name
                );
                continue;
            }

            const GRBVar start_time = time_var(edge.u, edge.data.start_state);
            const GRBVar end_time = time_var(edge.v, edge.data.end_state);
            const double big_m = data.time_limit + edge.data.time;
            core.model->addConstr(
                end_time >= start_time + edge.data.time -
                    big_m * (1.0 - core.y_vars[edge_id]),
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
            actual_departures += core.y_vars[edge_id];
        }
        if (graph.node(edge.u).kind == NodeSpec::Kind::Delivery &&
            edge.v == end_node_id) {
            actual_returns += core.y_vars[edge_id];
        }
    }
    core.model->addConstr(
        actual_departures == actual_returns,
        options.name_prefix + "_depot_balance"
    );
    core.model->update();
    return core;
}

AuxiliaryDurationSubproblemResult solve_auxiliary_duration_subproblem(
    const SPDPData& data,
    const MultiDiGraph& graph,
    double solver_time_limit,
    bool binary_y,
    bool add_time_constraints
) {
    TwoIndexCoreOptions core_options;
    core_options.objective = TwoIndexCoreObjective::Duration;
    core_options.binary_y = binary_y;
    core_options.add_time_constraints = add_time_constraints;
    core_options.solver_time_limit = solver_time_limit;
    core_options.gurobi_threads = 1;
    core_options.output_enabled = false;
    core_options.name_prefix = "vi44_aux";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    const auto solve_start = std::chrono::steady_clock::now();
    core.model->optimize();
    const auto solve_end = std::chrono::steady_clock::now();

    AuxiliaryDurationSubproblemResult result;
    result.status = core.model->get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.runtime_seconds =
        std::chrono::duration<double>(solve_end - solve_start).count();

    if (core.model->get(GRB_IntAttr_SolCount) > 0) {
        try {
            const double objective_value = core.model->get(GRB_DoubleAttr_ObjVal);
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
        const double objective_bound = core.model->get(GRB_DoubleAttr_ObjBound);
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

const char* vi44_subproblem_type_name(VI44SubproblemType type) {
    switch (type) {
        case VI44SubproblemType::LP:
            return "lp";
        case VI44SubproblemType::IP:
            return "ip";
    }
    throw std::runtime_error("Unknown VI-44 subproblem type.");
}

struct VehicleAssignmentMatchingArc {
    std::size_t pickup_index = 0U;
    std::size_t delivery_index = 0U;
    std::vector<GRBVar> variables;
};

struct VehicleAssignmentPhysicalArc {
    int tail = -1;
    int head = -1;
    GRBVar selected;
    GRBVar flow;
};

struct VehicleAssignmentFeasibilityResult {
    int status = 0;
    bool hit_time_limit = false;
    bool feasible = false;
    bool proven_infeasible = false;
    double runtime_seconds = 0.0;
    std::vector<VI44VehicleAssignmentVehicleResult> vehicles;
};

VehicleAssignmentFeasibilityResult solve_vehicle_assignment_feasibility(
    const SPDPData& data,
    int vehicle_count,
    const VI44VehicleAssignmentOptions& options,
    GRBEnv& environment,
    double solver_time_limit
) {
    const std::size_t request_count = data.requests.size();
    const int location_count = data.locations;
    if (vehicle_count <= 0 ||
        static_cast<std::size_t>(vehicle_count) > request_count) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment vehicle count must be in [1, |I|]."
        );
    }
    if (location_count <= 0) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment problem requires physical locations."
        );
    }

    GRBModel model(environment);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_Threads, 1);
    model.set(GRB_IntParam_DualReductions, 0);
    model.set(GRB_IntParam_SolutionLimit, 1);
    model.set(GRB_DoubleParam_FeasibilityTol, kAuxiliaryLPSolverTolerance);
    model.set(GRB_DoubleParam_OptimalityTol, kAuxiliaryLPSolverTolerance);
    if (solver_time_limit > 0.0) {
        model.set(GRB_DoubleParam_TimeLimit, solver_time_limit);
    }

    std::vector<std::vector<GRBVar>> pickup_assignment(request_count);
    std::vector<std::vector<GRBVar>> delivery_assignment(request_count);
    for (std::size_t request_index = 0U; request_index < request_count; ++request_index) {
        pickup_assignment[request_index].reserve(static_cast<std::size_t>(vehicle_count));
        delivery_assignment[request_index].reserve(static_cast<std::size_t>(vehicle_count));
        for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
            pickup_assignment[request_index].push_back(model.addVar(
                0.0,
                1.0,
                0.0,
                GRB_BINARY,
                "vi44_va_u_" + std::to_string(request_index) + "_" +
                    std::to_string(vehicle)
            ));
            delivery_assignment[request_index].push_back(model.addVar(
                0.0,
                1.0,
                0.0,
                GRB_BINARY,
                "vi44_va_v_" + std::to_string(request_index) + "_" +
                    std::to_string(vehicle)
            ));
        }
    }

    std::vector<VehicleAssignmentMatchingArc> matching_arcs;
    for (std::size_t pickup_index = 0U;
         pickup_index < request_count;
         ++pickup_index) {
        for (std::size_t delivery_index = 0U;
             delivery_index < request_count;
             ++delivery_index) {
            if (data.requests[pickup_index].container_type !=
                data.requests[delivery_index].container_type) {
                continue;
            }
            VehicleAssignmentMatchingArc arc;
            arc.pickup_index = pickup_index;
            arc.delivery_index = delivery_index;
            arc.variables.reserve(static_cast<std::size_t>(vehicle_count));
            for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
                arc.variables.push_back(model.addVar(
                    0.0,
                    1.0,
                    0.0,
                    GRB_BINARY,
                    "vi44_va_mu_" + std::to_string(pickup_index) + "_" +
                        std::to_string(delivery_index) + "_" +
                        std::to_string(vehicle)
                ));
            }
            matching_arcs.push_back(std::move(arc));
        }
    }

    std::vector<std::vector<GRBVar>> required_location;
    std::vector<std::vector<GRBVar>> path_start;
    std::vector<std::vector<GRBVar>> path_end;
    std::vector<std::vector<GRBVar>> source_flow;
    std::vector<std::vector<VehicleAssignmentPhysicalArc>> physical_arcs;
    if (options.add_tsp_bound) {
        required_location.resize(static_cast<std::size_t>(location_count));
        path_start.resize(static_cast<std::size_t>(location_count));
        path_end.resize(static_cast<std::size_t>(location_count));
        source_flow.resize(static_cast<std::size_t>(location_count));
        for (int location = 0; location < location_count; ++location) {
            for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
                required_location[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_b_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
                path_start[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_rho_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
                path_end[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_eta_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
                source_flow[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        static_cast<double>(location_count),
                        0.0,
                        GRB_CONTINUOUS,
                        "vi44_va_fs_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
            }
        }

        physical_arcs.resize(static_cast<std::size_t>(vehicle_count));
        for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
            auto& vehicle_arcs = physical_arcs[static_cast<std::size_t>(vehicle)];
            vehicle_arcs.reserve(
                static_cast<std::size_t>(location_count) *
                static_cast<std::size_t>(std::max(0, location_count - 1))
            );
            for (int tail = 0; tail < location_count; ++tail) {
                for (int head = 0; head < location_count; ++head) {
                    if (tail == head) {
                        continue;
                    }
                    VehicleAssignmentPhysicalArc arc;
                    arc.tail = tail;
                    arc.head = head;
                    arc.selected = model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_x_" + std::to_string(tail) + "_" +
                            std::to_string(head) + "_" + std::to_string(vehicle)
                    );
                    arc.flow = model.addVar(
                        0.0,
                        static_cast<double>(location_count),
                        0.0,
                        GRB_CONTINUOUS,
                        "vi44_va_f_" + std::to_string(tail) + "_" +
                            std::to_string(head) + "_" + std::to_string(vehicle)
                    );
                    vehicle_arcs.push_back(std::move(arc));
                }
            }
        }
    }

    model.update();
    GRBLinExpr feasibility_objective = 0.0;
    model.setObjective(feasibility_objective, GRB_MINIMIZE);

    for (std::size_t request_index = 0U; request_index < request_count; ++request_index) {
        GRBLinExpr pickup_cover = 0.0;
        GRBLinExpr delivery_cover = 0.0;
        for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
            pickup_cover +=
                pickup_assignment[request_index][static_cast<std::size_t>(vehicle)];
            delivery_cover +=
                delivery_assignment[request_index][static_cast<std::size_t>(vehicle)];
        }
        model.addConstr(
            pickup_cover == 1.0,
            "vi44_va_pickup_cover_" + std::to_string(request_index)
        );
        model.addConstr(
            delivery_cover == 1.0,
            "vi44_va_delivery_cover_" + std::to_string(request_index)
        );
    }

    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        GRBLinExpr pickup_count = 0.0;
        for (std::size_t request_index = 0U;
             request_index < request_count;
             ++request_index) {
            pickup_count +=
                pickup_assignment[request_index][static_cast<std::size_t>(vehicle)];
        }
        model.addConstr(
            pickup_count >= 1.0,
            "vi44_va_nonempty_" + std::to_string(vehicle)
        );
        if (vehicle + 1 < vehicle_count) {
            GRBLinExpr next_pickup_count = 0.0;
            for (std::size_t request_index = 0U;
                 request_index < request_count;
                 ++request_index) {
                next_pickup_count += pickup_assignment[request_index]
                    [static_cast<std::size_t>(vehicle + 1)];
            }
            model.addConstr(
                pickup_count >= next_pickup_count,
                "vi44_va_vehicle_size_order_" + std::to_string(vehicle)
            );
        }
    }
    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        for (std::size_t pickup_index = 0U;
             pickup_index < request_count;
             ++pickup_index) {
            GRBLinExpr matching_out = 0.0;
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                if (arc.pickup_index == pickup_index) {
                    matching_out += arc.variables[static_cast<std::size_t>(vehicle)];
                }
            }
            model.addConstr(
                matching_out ==
                    pickup_assignment[pickup_index][static_cast<std::size_t>(vehicle)],
                "vi44_va_match_out_" + std::to_string(pickup_index) + "_" +
                    std::to_string(vehicle)
            );
        }
        for (std::size_t delivery_index = 0U;
             delivery_index < request_count;
             ++delivery_index) {
            GRBLinExpr matching_in = 0.0;
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                if (arc.delivery_index == delivery_index) {
                    matching_in += arc.variables[static_cast<std::size_t>(vehicle)];
                }
            }
            model.addConstr(
                matching_in ==
                    delivery_assignment[delivery_index][static_cast<std::size_t>(vehicle)],
                "vi44_va_match_in_" + std::to_string(delivery_index) + "_" +
                    std::to_string(vehicle)
            );
        }
    }

    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        GRBLinExpr service_time = 0.0;
        for (std::size_t request_index = 0U;
             request_index < request_count;
             ++request_index) {
            service_time +=
                (data.time_pickup + data.time_empty) *
                pickup_assignment[request_index][static_cast<std::size_t>(vehicle)];
            service_time +=
                data.time_delivery *
                delivery_assignment[request_index][static_cast<std::size_t>(vehicle)];
        }

        if (options.add_container_bound) {
            GRBLinExpr container_workload = 0.0;
            for (std::size_t pickup_index = 0U;
                 pickup_index < request_count;
                 ++pickup_index) {
                const Request& pickup = data.requests[pickup_index];
                container_workload +=
                    data.time[static_cast<std::size_t>(pickup.from_id)]
                             [static_cast<std::size_t>(pickup.to_id)] *
                    pickup_assignment[pickup_index][static_cast<std::size_t>(vehicle)];
            }
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                const Request& pickup = data.requests[arc.pickup_index];
                const Request& delivery = data.requests[arc.delivery_index];
                container_workload +=
                    data.time[static_cast<std::size_t>(pickup.to_id)]
                             [static_cast<std::size_t>(delivery.from_id)] *
                    arc.variables[static_cast<std::size_t>(vehicle)];
            }
            model.addConstr(
                service_time + 0.5 * container_workload <= data.time_limit,
                "vi44_va_container_time_" + std::to_string(vehicle)
            );
        }

        if (!options.add_tsp_bound) {
            continue;
        }

        for (int location = 0; location < location_count; ++location) {
            GRBLinExpr activity = 0.0;
            GRBLinExpr pickup_root_activity = 0.0;
            for (std::size_t request_index = 0U;
                 request_index < request_count;
                 ++request_index) {
                const Request& request = data.requests[request_index];
                if (request.from_id == location) {
                    activity += pickup_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                    activity += delivery_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                    pickup_root_activity += pickup_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                }
                if (request.to_id == location) {
                    activity += pickup_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                }
            }

            const GRBVar b =
                required_location[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                b <= activity,
                "vi44_va_location_upper_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_start[static_cast<std::size_t>(location)]
                          [static_cast<std::size_t>(vehicle)] <= pickup_root_activity,
                "vi44_va_root_eligible_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_start[static_cast<std::size_t>(location)]
                          [static_cast<std::size_t>(vehicle)] <= b,
                "vi44_va_root_selected_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_end[static_cast<std::size_t>(location)]
                        [static_cast<std::size_t>(vehicle)] <= b,
                "vi44_va_end_selected_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );

            for (std::size_t request_index = 0U;
                 request_index < request_count;
                 ++request_index) {
                const Request& request = data.requests[request_index];
                if (request.from_id == location) {
                    model.addConstr(
                        b >= pickup_assignment[request_index]
                            [static_cast<std::size_t>(vehicle)],
                        "vi44_va_location_pickup_" + std::to_string(location) + "_" +
                            std::to_string(request_index) + "_" +
                            std::to_string(vehicle)
                    );
                    model.addConstr(
                        b >= delivery_assignment[request_index]
                            [static_cast<std::size_t>(vehicle)],
                        "vi44_va_location_delivery_" + std::to_string(location) + "_" +
                            std::to_string(request_index) + "_" +
                            std::to_string(vehicle)
                    );
                }
                if (request.to_id == location) {
                    model.addConstr(
                        b >= pickup_assignment[request_index]
                            [static_cast<std::size_t>(vehicle)],
                        "vi44_va_location_landfill_" + std::to_string(location) + "_" +
                            std::to_string(request_index) + "_" +
                            std::to_string(vehicle)
                    );
                }
            }
        }

        GRBLinExpr start_count = 0.0;
        GRBLinExpr end_count = 0.0;
        GRBLinExpr source_outflow = 0.0;
        GRBLinExpr selected_location_count = 0.0;
        for (int location = 0; location < location_count; ++location) {
            start_count += path_start[static_cast<std::size_t>(location)]
                                     [static_cast<std::size_t>(vehicle)];
            end_count += path_end[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            source_outflow += source_flow[static_cast<std::size_t>(location)]
                                         [static_cast<std::size_t>(vehicle)];
            selected_location_count +=
                required_location[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                source_flow[static_cast<std::size_t>(location)]
                           [static_cast<std::size_t>(vehicle)] <=
                    static_cast<double>(location_count) *
                    path_start[static_cast<std::size_t>(location)]
                              [static_cast<std::size_t>(vehicle)],
                "vi44_va_source_flow_link_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
        }
        model.addConstr(start_count == 1.0, "vi44_va_one_start_" + std::to_string(vehicle));
        model.addConstr(end_count == 1.0, "vi44_va_one_end_" + std::to_string(vehicle));
        model.addConstr(
            source_outflow == selected_location_count,
            "vi44_va_source_balance_" + std::to_string(vehicle)
        );

        GRBLinExpr tsp_time = 0.0;
        auto& vehicle_arcs = physical_arcs[static_cast<std::size_t>(vehicle)];
        for (VehicleAssignmentPhysicalArc& arc : vehicle_arcs) {
            const GRBVar tail_selected =
                required_location[static_cast<std::size_t>(arc.tail)]
                                 [static_cast<std::size_t>(vehicle)];
            const GRBVar head_selected =
                required_location[static_cast<std::size_t>(arc.head)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                arc.selected <= tail_selected,
                "vi44_va_arc_tail_" + std::to_string(arc.tail) + "_" +
                    std::to_string(arc.head) + "_" + std::to_string(vehicle)
            );
            model.addConstr(
                arc.selected <= head_selected,
                "vi44_va_arc_head_" + std::to_string(arc.tail) + "_" +
                    std::to_string(arc.head) + "_" + std::to_string(vehicle)
            );
            model.addConstr(
                arc.flow <= static_cast<double>(location_count) * arc.selected,
                "vi44_va_flow_arc_" + std::to_string(arc.tail) + "_" +
                    std::to_string(arc.head) + "_" + std::to_string(vehicle)
            );
            tsp_time +=
                data.time[static_cast<std::size_t>(arc.tail)]
                         [static_cast<std::size_t>(arc.head)] *
                arc.selected;
        }

        for (int location = 0; location < location_count; ++location) {
            GRBLinExpr incoming_arc = 0.0;
            GRBLinExpr outgoing_arc = 0.0;
            GRBLinExpr incoming_flow =
                source_flow[static_cast<std::size_t>(location)]
                           [static_cast<std::size_t>(vehicle)];
            GRBLinExpr outgoing_flow = 0.0;
            for (const VehicleAssignmentPhysicalArc& arc : vehicle_arcs) {
                if (arc.head == location) {
                    incoming_arc += arc.selected;
                    incoming_flow += arc.flow;
                }
                if (arc.tail == location) {
                    outgoing_arc += arc.selected;
                    outgoing_flow += arc.flow;
                }
            }
            const GRBVar b =
                required_location[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                path_start[static_cast<std::size_t>(location)]
                          [static_cast<std::size_t>(vehicle)] +
                    incoming_arc == b,
                "vi44_va_path_in_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_end[static_cast<std::size_t>(location)]
                        [static_cast<std::size_t>(vehicle)] +
                    outgoing_arc == b,
                "vi44_va_path_out_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                incoming_flow - outgoing_flow == b,
                "vi44_va_flow_balance_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
        }
        model.addConstr(
            service_time + tsp_time <= data.time_limit,
            "vi44_va_tsp_time_" + std::to_string(vehicle)
        );
    }

    const auto solve_start = std::chrono::steady_clock::now();
    model.optimize();
    const auto solve_end = std::chrono::steady_clock::now();

    VehicleAssignmentFeasibilityResult result;
    result.status = model.get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.runtime_seconds =
        std::chrono::duration<double>(solve_end - solve_start).count();
    result.feasible = model.get(GRB_IntAttr_SolCount) > 0;
    result.proven_infeasible = result.status == GRB_INFEASIBLE;

    if (!result.feasible && !result.proven_infeasible &&
        result.status != GRB_TIME_LIMIT) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment problem stopped with unsupported status " +
            std::to_string(result.status) + "."
        );
    }

    if (!result.feasible) {
        return result;
    }

    result.vehicles.reserve(static_cast<std::size_t>(vehicle_count));
    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        VI44VehicleAssignmentVehicleResult vehicle_result;
        vehicle_result.vehicle_index = vehicle;
        for (std::size_t request_index = 0U;
             request_index < request_count;
             ++request_index) {
            if (pickup_assignment[request_index][static_cast<std::size_t>(vehicle)]
                    .get(GRB_DoubleAttr_X) > 0.5) {
                ++vehicle_result.pickup_count;
            }
            if (delivery_assignment[request_index][static_cast<std::size_t>(vehicle)]
                    .get(GRB_DoubleAttr_X) > 0.5) {
                ++vehicle_result.delivery_count;
            }
        }
        vehicle_result.service_time =
            (data.time_pickup + data.time_empty) *
                static_cast<double>(vehicle_result.pickup_count) +
            data.time_delivery * static_cast<double>(vehicle_result.delivery_count);

        double active_travel_lb = 0.0;
        if (options.add_container_bound) {
            double workload = 0.0;
            for (std::size_t pickup_index = 0U;
                 pickup_index < request_count;
                 ++pickup_index) {
                if (pickup_assignment[pickup_index][static_cast<std::size_t>(vehicle)]
                        .get(GRB_DoubleAttr_X) <= 0.5) {
                    continue;
                }
                const Request& pickup = data.requests[pickup_index];
                workload +=
                    data.time[static_cast<std::size_t>(pickup.from_id)]
                             [static_cast<std::size_t>(pickup.to_id)];
            }
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                if (arc.variables[static_cast<std::size_t>(vehicle)]
                        .get(GRB_DoubleAttr_X) <= 0.5) {
                    continue;
                }
                const Request& pickup = data.requests[arc.pickup_index];
                const Request& delivery = data.requests[arc.delivery_index];
                workload +=
                    data.time[static_cast<std::size_t>(pickup.to_id)]
                             [static_cast<std::size_t>(delivery.from_id)];
            }
            vehicle_result.container_lb = 0.5 * workload;
            active_travel_lb = vehicle_result.container_lb;
        }
        if (options.add_tsp_bound) {
            double tsp_lb = 0.0;
            for (const VehicleAssignmentPhysicalArc& arc :
                 physical_arcs[static_cast<std::size_t>(vehicle)]) {
                if (arc.selected.get(GRB_DoubleAttr_X) > 0.5) {
                    tsp_lb +=
                        data.time[static_cast<std::size_t>(arc.tail)]
                                 [static_cast<std::size_t>(arc.head)];
                }
            }
            vehicle_result.tsp_lb = tsp_lb;
            active_travel_lb = std::max(active_travel_lb, tsp_lb);
        }
        vehicle_result.active_duration_lb =
            vehicle_result.service_time + active_travel_lb;
        result.vehicles.push_back(vehicle_result);
    }
    return result;
}

VI44KMinResult compute_vehicle_assignment_k_min(
    const SPDPData& data,
    const VI44VehicleAssignmentOptions& options
) {
    if (!options.add_tsp_bound && !options.add_container_bound) {
        throw std::runtime_error(
            "VI-44 vehicle assignment requires the TSP bound, the container bound, "
            "or both."
        );
    }
    if (!std::isfinite(options.time_limit) || options.time_limit < 0.0) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment time limit must be finite and nonnegative."
        );
    }

    VI44KMinResult result;
    result.vehicle_assignment_enabled = true;
    if (data.requests.empty()) {
        result.vehicle_assignment_has_certified_bound = true;
        result.vehicle_assignment_optimal = true;
        return result;
    }

    GRBEnv environment(true);
    environment.set(GRB_IntParam_OutputFlag, 0);
    environment.start();

    const auto overall_start = std::chrono::steady_clock::now();
    int certified_k_min = 1;
    for (int vehicle_count = 1;
         vehicle_count <= static_cast<int>(data.requests.size());
         ++vehicle_count) {
        double remaining_time = 0.0;
        if (options.time_limit > 0.0) {
            const double elapsed =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - overall_start
                ).count();
            remaining_time = std::max(0.0, options.time_limit - elapsed);
            if (remaining_time <= kTolerance) {
                result.vehicle_assignment_hit_time_limit = true;
                result.vehicle_assignment_k_min = certified_k_min;
                result.vehicle_assignment_has_certified_bound = true;
                result.vehicle_assignment_runtime_seconds = elapsed;
                return result;
            }
        }

        const VehicleAssignmentFeasibilityResult feasibility =
            solve_vehicle_assignment_feasibility(
                data,
                vehicle_count,
                options,
                environment,
                remaining_time
            );
        result.vehicle_assignment_status = feasibility.status;
        result.vehicle_assignment_tested_vehicle_count = vehicle_count;

        if (feasibility.feasible) {
            result.vehicle_assignment_k_min = vehicle_count;
            result.vehicle_assignment_has_certified_bound = true;
            result.vehicle_assignment_has_feasible_solution = true;
            result.vehicle_assignment_optimal = true;
            result.vehicle_assignment_vehicles = feasibility.vehicles;
            break;
        }
        if (feasibility.proven_infeasible) {
            certified_k_min = vehicle_count + 1;
            continue;
        }

        result.vehicle_assignment_hit_time_limit = feasibility.hit_time_limit;
        result.vehicle_assignment_k_min = certified_k_min;
        result.vehicle_assignment_has_certified_bound = true;
        break;
    }

    result.vehicle_assignment_runtime_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - overall_start
        ).count();
    if (!result.vehicle_assignment_has_certified_bound) {
        result.vehicle_assignment_k_min = certified_k_min;
        result.vehicle_assignment_has_certified_bound = true;
    }
    return result;
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
    const VI44KMinOptions& options
) {
    if (options.precomputed_result != nullptr) {
        return *options.precomputed_result;
    }
    if (!options.use_cor && !options.use_subproblem &&
        !options.use_vehicle_assignment) {
        throw std::runtime_error(
            "VI-44 is enabled, but every k_min method is disabled."
        );
    }

    VI44KMinResult result;
    if (options.use_cor) {
        result.cor_enabled = true;
        result.cor_k_min = cor_route_lower_bound_kmin(data);
        result.selected_k_min = std::max(result.selected_k_min, result.cor_k_min);
    }

    if (options.use_subproblem) {
        if (!std::isfinite(options.subproblem.time_limit) ||
            options.subproblem.time_limit < 0.0) {
            throw std::runtime_error(
                "The auxiliary VI-44 subproblem time limit must be finite and "
                "nonnegative."
            );
        }
        result.subproblem_enabled = true;
        const AuxiliaryDurationSubproblemResult auxiliary =
            solve_auxiliary_duration_subproblem(
                data,
                graph,
                options.subproblem.time_limit,
                options.subproblem.type == VI44SubproblemType::IP,
                options.subproblem.add_time_constraints
            );
        result.subproblem_status = auxiliary.status;
        result.subproblem_hit_time_limit = auxiliary.hit_time_limit;
        result.subproblem_has_certified_bound = auxiliary.has_certified_bound;
        result.subproblem_objective_value = auxiliary.objective_value;
        result.subproblem_objective_bound = auxiliary.objective_bound;
        result.subproblem_runtime_seconds = auxiliary.runtime_seconds;

        if (auxiliary.has_certified_bound) {
            const double scale = std::max(
                {1.0, data.time_limit, std::abs(auxiliary.objective_bound)}
            );
            result.subproblem_numerical_tolerance =
                10.0 * kAuxiliaryLPSolverTolerance * scale;
            result.subproblem_safe_lower_bound = std::max(
                0.0,
                auxiliary.objective_bound -
                    result.subproblem_numerical_tolerance
            );
            result.subproblem_k_min = std::max(
                0,
                static_cast<int>(
                    std::ceil(
                        result.subproblem_safe_lower_bound / data.time_limit
                    )
                )
            );
            result.selected_k_min =
                std::max(result.selected_k_min, result.subproblem_k_min);
        }
    }

    if (options.use_vehicle_assignment) {
        VI44KMinResult vehicle_assignment_result =
            compute_vehicle_assignment_k_min(data, options.vehicle_assignment);
        result.vehicle_assignment_enabled = true;
        result.vehicle_assignment_k_min =
            vehicle_assignment_result.vehicle_assignment_k_min;
        result.vehicle_assignment_status =
            vehicle_assignment_result.vehicle_assignment_status;
        result.vehicle_assignment_hit_time_limit =
            vehicle_assignment_result.vehicle_assignment_hit_time_limit;
        result.vehicle_assignment_has_certified_bound =
            vehicle_assignment_result.vehicle_assignment_has_certified_bound;
        result.vehicle_assignment_has_feasible_solution =
            vehicle_assignment_result.vehicle_assignment_has_feasible_solution;
        result.vehicle_assignment_optimal =
            vehicle_assignment_result.vehicle_assignment_optimal;
        result.vehicle_assignment_tested_vehicle_count =
            vehicle_assignment_result.vehicle_assignment_tested_vehicle_count;
        result.vehicle_assignment_runtime_seconds =
            vehicle_assignment_result.vehicle_assignment_runtime_seconds;
        result.vehicle_assignment_vehicles =
            std::move(vehicle_assignment_result.vehicle_assignment_vehicles);
        if (result.vehicle_assignment_has_certified_bound) {
            result.selected_k_min =
                std::max(result.selected_k_min, result.vehicle_assignment_k_min);
        }
    }

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
                options.vi_44_k_min_options
            );
        rows.push_back(build_vi44_row(graph, result.selected_k_min));
        if (options.log_stream != nullptr) {
            std::ostringstream message;
            message << std::setprecision(15)
                    << "[pstep-vi] VI44 k_min"
                    << " cor_enabled=" << (result.cor_enabled ? 1 : 0)
                    << " cor_k_min=" << result.cor_k_min
                    << " subproblem_enabled="
                    << (result.subproblem_enabled ? 1 : 0);
            if (result.subproblem_enabled) {
                message << " subproblem_type="
                        << vi44_subproblem_type_name(
                               options.vi_44_k_min_options.subproblem.type
                           )
                        << " subproblem_status=" << result.subproblem_status
                        << " subproblem_time_limit="
                        << options.vi_44_k_min_options.subproblem.time_limit
                        << " subproblem_time_constraints="
                        << (options.vi_44_k_min_options.subproblem
                                    .add_time_constraints
                                ? 1
                                : 0)
                        << " subproblem_hit_time_limit="
                        << (result.subproblem_hit_time_limit ? 1 : 0)
                        << " subproblem_has_certified_bound="
                        << (result.subproblem_has_certified_bound ? 1 : 0)
                        << " subproblem_objective="
                        << result.subproblem_objective_value
                        << " subproblem_bound="
                        << result.subproblem_objective_bound
                        << " subproblem_safe_bound="
                        << result.subproblem_safe_lower_bound
                        << " subproblem_tolerance="
                        << result.subproblem_numerical_tolerance
                        << " subproblem_k_min=" << result.subproblem_k_min
                        << " subproblem_runtime_seconds="
                        << result.subproblem_runtime_seconds
                        // Backward-compatible token consumed by the summary script.
                        << " sub_lp_runtime_seconds="
                        << result.subproblem_runtime_seconds;
            }
            message << " vehicle_assignment_enabled="
                    << (result.vehicle_assignment_enabled ? 1 : 0);
            if (result.vehicle_assignment_enabled) {
                message << " vehicle_assignment_tsp_bound="
                        << (options.vi_44_k_min_options.vehicle_assignment
                                    .add_tsp_bound
                                ? 1
                                : 0)
                        << " vehicle_assignment_container_bound="
                        << (options.vi_44_k_min_options.vehicle_assignment
                                    .add_container_bound
                                ? 1
                                : 0)
                        << " vehicle_assignment_status="
                        << result.vehicle_assignment_status
                        << " vehicle_assignment_hit_time_limit="
                        << (result.vehicle_assignment_hit_time_limit ? 1 : 0)
                        << " vehicle_assignment_has_certified_bound="
                        << (result.vehicle_assignment_has_certified_bound ? 1 : 0)
                        << " vehicle_assignment_has_feasible_solution="
                        << (result.vehicle_assignment_has_feasible_solution ? 1 : 0)
                        << " vehicle_assignment_optimal="
                        << (result.vehicle_assignment_optimal ? 1 : 0)
                        << " vehicle_assignment_tested_vehicle_count="
                        << result.vehicle_assignment_tested_vehicle_count
                        << " vehicle_assignment_k_min="
                        << result.vehicle_assignment_k_min
                        << " vehicle_assignment_runtime_seconds="
                        << result.vehicle_assignment_runtime_seconds;
            }
            message << " selected=" << result.selected_k_min;
            *options.log_stream << message.str() << '\n';

            for (const VI44VehicleAssignmentVehicleResult& vehicle :
                 result.vehicle_assignment_vehicles) {
                *options.log_stream
                    << std::setprecision(15)
                    << "[pstep-vi] VI44 vehicle_assignment_vehicle="
                    << vehicle.vehicle_index
                    << " pickup_count=" << vehicle.pickup_count
                    << " delivery_count=" << vehicle.delivery_count
                    << " service_time=" << vehicle.service_time
                    << " container_lb=" << vehicle.container_lb
                    << " tsp_lb=" << vehicle.tsp_lb
                    << " active_duration_lb=" << vehicle.active_duration_lb
                    << '\n';
            }
        }
    }
    return rows;
}

DirectTwoIndexResult solve_direct_two_index_ip(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& options
) {
    TwoIndexCoreOptions core_options;
    core_options.objective =
        options.objective == DirectTwoIndexObjective::Duration
            ? TwoIndexCoreObjective::Duration
            : TwoIndexCoreObjective::OriginalCost;
    core_options.binary_y = true;
    core_options.add_time_constraints = options.add_time_constraints;
    core_options.solver_time_limit = options.solver_time_limit;
    core_options.gurobi_threads = options.gurobi_threads;
    core_options.output_enabled = true;
    core_options.gurobi_log_path = options.gurobi_log_path;
    core_options.name_prefix = "direct_two_index";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    if (!options.initial_edge_start.empty()) {
        if (options.initial_edge_start.size() != core.y_vars.size()) {
            throw std::runtime_error(
                "The direct two-index MIP start size does not match the edge count."
            );
        }
        for (std::size_t edge_id = 0; edge_id < core.y_vars.size(); ++edge_id) {
            core.y_vars[edge_id].set(
                GRB_DoubleAttr_Start,
                options.initial_edge_start[edge_id]
            );
        }
    }

    const std::vector<PstepValidInequalityRow> vi_rows =
        build_pstep_valid_inequality_rows(
            data,
            graph,
            options.valid_inequalities
        );
    for (const PstepValidInequalityRow& row : vi_rows) {
        GRBLinExpr expression = 0.0;
        for (const auto& term : row.edge_terms) {
            if (term.first >= core.y_vars.size()) {
                throw std::runtime_error(
                    "A direct two-index valid inequality references an invalid edge."
                );
            }
            expression += term.second * core.y_vars[term.first];
        }
        if (row.sense == PstepValidInequalitySense::GreaterEqual) {
            core.model->addConstr(expression >= row.rhs, row.name);
        } else {
            core.model->addConstr(expression <= row.rhs, row.name);
        }
    }
    core.model->update();

    DirectTwoIndexResult result;
    result.variable_count = core.model->get(GRB_IntAttr_NumVars);
    result.constraint_count = core.model->get(GRB_IntAttr_NumConstrs);
    result.valid_inequality_count = static_cast<int>(vi_rows.size());

    core.model->optimize();
    result.status = core.model->get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.solved_to_optimality = result.status == GRB_OPTIMAL;
    result.runtime_seconds = core.model->get(GRB_DoubleAttr_Runtime);
    result.has_feasible_solution = core.model->get(GRB_IntAttr_SolCount) > 0;

    try {
        const double objective_bound = core.model->get(GRB_DoubleAttr_ObjBound);
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

    if (!result.has_feasible_solution) {
        return result;
    }

    result.objective_value = core.model->get(GRB_DoubleAttr_ObjVal);
    if (result.has_certified_bound && std::abs(result.objective_value) > kTolerance) {
        result.gap_percent =
            100.0 * std::abs(result.objective_value - result.objective_bound) /
            std::abs(result.objective_value);
    } else if (result.has_certified_bound) {
        result.gap_percent = 0.0;
    }

    result.total_duration = 0.0;
    result.total_original_cost = 0.0;
    result.edge_values.resize(graph.number_of_edges(), 0.0);
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const double value = core.y_vars[edge_id].get(GRB_DoubleAttr_X);
        result.edge_values[edge_id] = value;
        result.total_duration += graph.edges()[edge_id].data.time * value;
        result.total_original_cost += graph.edges()[edge_id].data.cost * value;
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (value > 0.5 && edge.u == 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            ++result.vehicle_count;
        }
    }
    return result;
}

InitialIncumbentSolveResult solve_fixed_k_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
) {
    if (options.vehicle_count <= 0) {
        throw std::runtime_error(
            "Initial-incumbent generation requires a positive fixed vehicle count."
        );
    }
    if (!std::isfinite(options.solver_time_limit) ||
        options.solver_time_limit < 0.0) {
        throw std::runtime_error(
            "Initial-incumbent time limit must be finite and nonnegative."
        );
    }

    TwoIndexCoreOptions core_options;
    core_options.objective = TwoIndexCoreObjective::Duration;
    core_options.binary_y = true;
    core_options.add_time_constraints = true;
    core_options.solver_time_limit = options.solver_time_limit;
    core_options.gurobi_threads = options.gurobi_threads;
    core_options.output_enabled = true;
    core_options.gurobi_log_path = options.gurobi_log_path;
    core_options.name_prefix = "initial_incumbent";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    GRBLinExpr actual_departures = 0.0;
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            actual_departures += core.y_vars[edge_id];
        }
    }
    core.model->addConstr(
        actual_departures == static_cast<double>(options.vehicle_count),
        "initial_incumbent_fixed_vehicle_count"
    );
    core.model->set(GRB_IntParam_SolutionLimit, 1);
    core.model->update();
    core.model->optimize();

    InitialIncumbentSolveResult result;
    result.status = core.model->get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.hit_solution_limit = result.status == GRB_SOLUTION_LIMIT;
    result.infeasible = result.status == GRB_INFEASIBLE;
    result.runtime_seconds = core.model->get(GRB_DoubleAttr_Runtime);
    result.has_feasible_solution = core.model->get(GRB_IntAttr_SolCount) > 0;
    if (!result.has_feasible_solution) {
        return result;
    }

    result.total_duration = 0.0;
    result.total_original_cost = 0.0;
    result.edge_values.resize(graph.number_of_edges(), 0.0);
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const double value = core.y_vars[edge_id].get(GRB_DoubleAttr_X);
        result.edge_values[edge_id] = value;
        result.total_duration += graph.edges()[edge_id].data.time * value;
        result.total_original_cost += graph.edges()[edge_id].data.cost * value;
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (value > 0.5 && edge.u == 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            ++result.vehicle_count;
        }
    }
    return result;
}

}  // namespace spdp
