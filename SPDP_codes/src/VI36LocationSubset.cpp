#include "VI36LocationSubset.h"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace spdp {
namespace {

enum class LocationFamily {
    Pickup,
    Treatment,
    Delivery,
};

bool contains_location(const std::vector<int>& subset, int location) {
    return std::binary_search(subset.begin(), subset.end(), location);
}

void enumerate_fixed_size_subsets(
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
        enumerate_fixed_size_subsets(locations, target_size, index + 1U, current, subsets);
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
        enumerate_fixed_size_subsets(locations, subset_size, 0U, current, subsets);
    }
    return subsets;
}

std::string family_name(LocationFamily family) {
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

std::string subset_name_suffix(const std::vector<int>& subset) {
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

std::vector<VI36CombinedRow> build_family_rows(
    const MultiDiGraph& graph,
    const std::map<int, int>& count_by_location,
    LocationFamily family,
    std::size_t max_subset_size,
    bool skip_full_location_set
) {
    std::vector<VI36CombinedRow> rows;
    const std::vector<std::vector<int>> subsets = enumerate_location_subsets(
        count_by_location,
        max_subset_size,
        skip_full_location_set
    );
    rows.reserve(subsets.size());

    for (const std::vector<int>& subset : subsets) {
        VI36CombinedRow row;
        row.name = "vi36_combined_" + family_name(family) + "_" +
                   subset_name_suffix(subset);
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
                    !edge.data.sequence_pi.empty(),
                    subset,
                    service_kind
                );
            }

            const int combined_coefficient = coefficients.first + coefficients.second;
            if (combined_coefficient > 0) {
                row.edge_terms.emplace_back(edge_id, combined_coefficient);
            }
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

}  // namespace

std::vector<VI36CombinedRow> build_vi36_combined_location_subset_rows(
    const MultiDiGraph& graph,
    std::size_t max_subset_size,
    bool skip_full_location_sets
) {
    if (max_subset_size < 1U) {
        throw std::runtime_error("VI-36 subset max size must be at least 1.");
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

    std::vector<VI36CombinedRow> rows;
    const auto append_family = [&](const std::map<int, int>& counts, LocationFamily family) {
        std::vector<VI36CombinedRow> family_rows = build_family_rows(
            graph,
            counts,
            family,
            max_subset_size,
            skip_full_location_sets
        );
        rows.insert(
            rows.end(),
            std::make_move_iterator(family_rows.begin()),
            std::make_move_iterator(family_rows.end())
        );
    };
    append_family(pickup_count_by_location, LocationFamily::Pickup);
    append_family(treatment_count_by_location, LocationFamily::Treatment);
    append_family(delivery_count_by_location, LocationFamily::Delivery);
    return rows;
}

}  // namespace spdp
