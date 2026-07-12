#include "RequestBlockSEC.h"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace spdp {
namespace {

struct RequestServiceNodes {
    NodeId pickup_node = -1;
    NodeId delivery_node = -1;
};

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

std::string subset_name_suffix(const std::vector<std::size_t>& request_indices) {
    std::string suffix;
    for (std::size_t request_index : request_indices) {
        suffix += "_" + std::to_string(request_index + 1U);
    }
    return suffix;
}

}  // namespace

std::vector<RequestBlockSECRow> build_request_block_sec_rows(
    const SPDPData& data,
    const MultiDiGraph& graph,
    std::size_t max_subset_size
) {
    if (max_subset_size < 2U) {
        throw std::runtime_error("Request-block SEC max size must be at least 2.");
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

    std::vector<RequestBlockSECRow> rows;
    rows.reserve(request_subsets.size());
    for (const std::vector<std::size_t>& request_indices : request_subsets) {
        std::vector<std::uint8_t> in_block(graph.number_of_nodes(), 0U);
        for (std::size_t request_index : request_indices) {
            const RequestServiceNodes& service_nodes = service_nodes_by_request[request_index];
            in_block[static_cast<std::size_t>(service_nodes.pickup_node)] = 1U;
            in_block[static_cast<std::size_t>(service_nodes.delivery_node)] = 1U;
        }

        RequestBlockSECRow row;
        row.name = "vi_request_block_sec_req" + subset_name_suffix(request_indices);
        row.rhs = static_cast<double>(2U * request_indices.size() - 1U);
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            if (in_block[static_cast<std::size_t>(edge.u)] != 0U &&
                in_block[static_cast<std::size_t>(edge.v)] != 0U) {
                row.edge_ids.push_back(edge_id);
            }
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

}  // namespace spdp
