#include "PstepValidInequality.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iterator>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace spdp {
namespace {

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

const char* vi44_subproblem_type_name(VI44SubproblemType type) {
    switch (type) {
        case VI44SubproblemType::LP:
            return "lp";
        case VI44SubproblemType::IP:
            return "ip";
    }
    throw std::runtime_error("Unknown VI-44 subproblem type.");
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

}  // namespace spdp
