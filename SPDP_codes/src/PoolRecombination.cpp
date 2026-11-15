#include "PoolRecombination.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "PstepFormulation.h"
#include "TwoIndexModelCore.h"
#include "gurobi_c++.h"

namespace spdp {
namespace {

using Clock = std::chrono::steady_clock;

// Kind of an edge between two service nodes: tail kind, head kind, tail
// location, head location, head container type, facility sequence in between.
using EdgeKind = std::tuple<char, char, int, int, int, std::vector<int>>;

struct ServiceOp {
    char kind = 'P';
    int location = -1;
    int container_type = -1;
    std::vector<int> facilities;  // visited between the previous service and this one
};

// The multigraph keeps one order of the operations made at one location without
// travelling (pruning rules 41a, 41b, 41c): a delivery before a pickup, a pickup
// before an emptying, an emptying before a delivery. PALNS applies the same swaps
// before it writes a solution. After the swaps the vehicle may become empty
// inside the block, which cuts it into several blocks.
std::vector<std::vector<PoolOp>> canonical_blocks(std::vector<PoolOp> ops) {
    const auto swaps = [](char a, char b) {
        return (a == 'P' && b == 'D') || (a == 'E' && b == 'P') || (a == 'D' && b == 'E');
    };
    for (int pass = 0; pass < 10000; ++pass) {
        bool changed = false;
        for (std::size_t i = 0; i + 1 < ops.size(); ++i) {
            if (swaps(ops[i].kind, ops[i + 1].kind) && ops[i].location == ops[i + 1].location) {
                std::swap(ops[i], ops[i + 1]);
                changed = true;
            }
        }
        if (!changed) break;
    }
    std::vector<std::vector<PoolOp>> out;
    std::vector<PoolOp> current;
    int load = 0;
    for (const PoolOp& op : ops) {
        current.push_back(op);
        if (op.kind == 'P') ++load;
        else if (op.kind == 'D') --load;
        if (load == 0) {
            out.push_back(current);
            current.clear();
        }
    }
    if (!current.empty()) {
        throw std::runtime_error("A pool block does not end with an empty vehicle.");
    }
    return out;
}

void add_block_edge_kinds(const std::vector<PoolOp>& block, std::set<EdgeKind>& kinds) {
    std::vector<ServiceOp> services;
    std::vector<int> facilities;
    for (const PoolOp& op : block) {
        if (op.kind == 'E') {
            if (facilities.empty() || facilities.back() != op.location) {
                facilities.push_back(op.location);
            }
            continue;
        }
        services.push_back({op.kind, op.location, op.container_type, facilities});
        facilities.clear();
    }
    for (std::size_t i = 1; i < services.size(); ++i) {
        const ServiceOp& u = services[i - 1];
        const ServiceOp& v = services[i];
        kinds.emplace(u.kind, v.kind, u.location, v.location, v.container_type, v.facilities);
        // Two pickups (deliveries) in a row at one location: the multigraph keeps
        // them in increasing request order (rules 40a, 40c), so the classes of the
        // two may appear in either order.
        if (v.facilities.empty() && u.kind == v.kind && u.location == v.location) {
            kinds.emplace(v.kind, u.kind, v.location, u.location, u.container_type, std::vector<int>{});
            if (i >= 2) {
                const ServiceOp& w = services[i - 2];
                kinds.emplace(w.kind, v.kind, w.location, v.location, v.container_type, u.facilities);
            }
            if (i + 1 < services.size()) {
                const ServiceOp& w = services[i + 1];
                kinds.emplace(u.kind, w.kind, u.location, w.location, w.container_type, w.facilities);
            }
        }
    }
}

char kind_char(NodeSpec::Kind kind) {
    switch (kind) {
        case NodeSpec::Kind::Start: return 'S';
        case NodeSpec::Kind::End: return 'T';
        case NodeSpec::Kind::Pickup: return 'P';
        case NodeSpec::Kind::Delivery: return 'D';
    }
    return '?';
}

// Records the time at which each new incumbent of the restricted model was found.
class TraceCallback : public GRBCallback {
public:
    double best_value = GRB_INFINITY;
    double best_seconds = -1.0;
    int solutions = 0;

protected:
    void callback() override {
        if (where != GRB_CB_MIPSOL) return;
        const double value = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
        ++solutions;
        if (value < best_value - 1e-6) {
            best_value = value;
            best_seconds = getDoubleInfo(GRB_CB_RUNTIME);
        }
    }
};

}  // namespace

PoolRecombinationResult recombine_pool_edges(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<PoolBlock>& blocks,
    const std::vector<int>& start_edge_ids,
    int vehicle_count,
    const PoolRecombinationOptions& options
) {
    PoolRecombinationResult result;
    const auto wall_start = Clock::now();
    const auto seconds = [&]() {
        return std::chrono::duration<double>(Clock::now() - wall_start).count();
    };
    std::ostringstream log;
    log.precision(10);
    log << "[palns-recomb] k=" << vehicle_count
        << " delta=" << options.delta
        << " time_limit=" << options.time_limit
        << " threads=" << options.threads
        << " blocks=" << blocks.size();

    double start_cost = 0.0;
    for (int edge_id : start_edge_ids) {
        start_cost += graph.edges()[static_cast<std::size_t>(edge_id)].data.cost;
    }
    log << " start_cost=" << start_cost;

    const auto finish = [&](const std::string& status) {
        result.seconds = seconds();
        log << " status=" << status << " improved=" << (result.improved ? 1 : 0)
            << " seconds=" << result.seconds;
        result.log = log.str();
        return result;
    };

    if (blocks.empty() || start_edge_ids.empty()) {
        return finish("skipped-empty-pool");
    }

    try {
        // 1. Edge kinds of the good blocks.
        double best = blocks.front().min_cost;
        for (const PoolBlock& block : blocks) best = std::min(best, block.min_cost);
        std::set<EdgeKind> kinds;
        std::size_t used_blocks = 0;
        for (const PoolBlock& block : blocks) {
            if (block.min_cost > best + options.delta + 1e-9) continue;
            ++used_blocks;
            for (const std::vector<PoolOp>& part : canonical_blocks(block.ops)) {
                add_block_edge_kinds(part, kinds);
            }
        }

        // 2. Edge subset.
        const State empty_state{};
        const NodeId end_node_id = graph.end_node_id();
        std::vector<char> keep(graph.number_of_edges(), 0);
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            const NodeSpec& u = graph.node(edge.u);
            const NodeSpec& v = graph.node(edge.v);
            if (edge.data.start_state == empty_state &&
                (u.kind == NodeSpec::Kind::Start || u.kind == NodeSpec::Kind::Delivery)) {
                keep[edge_id] = 1;
                continue;
            }
            if (edge.u == 0 || edge.v == end_node_id) continue;
            const EdgeKind kind(
                kind_char(u.kind), kind_char(v.kind), u.location, v.location,
                v.container_type ? *v.container_type : -1, edge.data.sequence_pi);
            if (kinds.count(kind) > 0) keep[edge_id] = 1;
        }
        std::size_t start_added = 0;
        for (int edge_id : start_edge_ids) {
            if (!keep[static_cast<std::size_t>(edge_id)]) {
                keep[static_cast<std::size_t>(edge_id)] = 1;
                ++start_added;
            }
        }

        MultiDiGraph sub;
        for (const NodeSpec& node : graph.nodes()) sub.add_node(node);
        std::vector<int> local_to_main;
        std::vector<int> main_to_local(graph.number_of_edges(), -1);
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            if (!keep[edge_id]) continue;
            const EdgeRecord& edge = graph.edges()[edge_id];
            main_to_local[edge_id] = static_cast<int>(local_to_main.size());
            local_to_main.push_back(static_cast<int>(edge_id));
            sub.add_edge(edge.u, edge.v, edge.data);
        }
        log << " best=" << best << " blocks_used=" << used_blocks
            << " edge_kinds=" << kinds.size()
            << " edges=" << local_to_main.size() << " edges_all=" << graph.number_of_edges()
            << " start_edges_added=" << start_added;

        // 3. Real model on the subset: two-index model with the time flow.
        detail::TwoIndexCoreOptions core_options;
        core_options.objective = detail::TwoIndexCoreObjective::OriginalCost;
        core_options.binary_y = true;
        core_options.add_time_flow_formulation = true;
        core_options.time_flow_state_disaggregated = true;
        core_options.gurobi_threads = options.threads;
        core_options.output_enabled = false;
        core_options.name_prefix = "pool_recomb";
        detail::TwoIndexCoreModel core =
            detail::build_two_index_model_core(data, sub, core_options);

        GRBLinExpr departures = 0.0;
        for (std::size_t local = 0; local < local_to_main.size(); ++local) {
            const EdgeRecord& edge = sub.edges()[local];
            if (edge.u == 0 && sub.node(edge.v).kind == NodeSpec::Kind::Pickup) {
                departures += core.y_vars[local];
            }
        }
        core.model->addConstr(
            departures == static_cast<double>(vehicle_count), "pool_recomb_fixed_vehicle_count");
        for (std::size_t local = 0; local < local_to_main.size(); ++local) {
            core.y_vars[local].set(GRB_DoubleAttr_Start, 0.0);
        }
        for (int edge_id : start_edge_ids) {
            core.y_vars[static_cast<std::size_t>(main_to_local[static_cast<std::size_t>(edge_id)])]
                .set(GRB_DoubleAttr_Start, 1.0);
        }
        core.model->set(GRB_IntParam_MIPFocus, options.mip_focus);
        const double build_seconds = seconds();
        const double remaining = options.time_limit - build_seconds;
        log << " build_seconds=" << build_seconds;
        if (remaining <= 0.0) {
            return finish("skipped-no-time");
        }
        core.model->set(GRB_DoubleParam_TimeLimit, remaining);

        TraceCallback trace;
        core.model->setCallback(&trace);
        core.model->optimize();
        result.ran = true;

        const int status = core.model->get(GRB_IntAttr_Status);
        const int solution_count = core.model->get(GRB_IntAttr_SolCount);
        const std::string status_name =
            status == GRB_OPTIMAL ? "optimal"
            : (status == GRB_TIME_LIMIT ? "time-limit" : "gurobi-" + std::to_string(status));
        log << " solve_seconds=" << core.model->get(GRB_DoubleAttr_Runtime)
            << " solutions=" << solution_count
            << " nodes=" << core.model->get(GRB_DoubleAttr_NodeCount);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            log << " bound=" << core.model->get(GRB_DoubleAttr_ObjBound);
        }
        if (solution_count <= 0) {
            return finish(status_name + "-no-solution");
        }
        const double value = core.model->get(GRB_DoubleAttr_ObjVal);
        log << " value=" << value << " best_found_at=" << trace.best_seconds;
        if (value >= start_cost - 1e-6) {
            return finish(status_name);
        }

        // 4. A cheaper solution: validate it as routes of the full graph.
        std::vector<int> active;
        double duration = 0.0;
        double cost = 0.0;
        for (std::size_t local = 0; local < local_to_main.size(); ++local) {
            if (core.y_vars[local].get(GRB_DoubleAttr_X) <= 0.5) continue;
            const int main_id = local_to_main[local];
            active.push_back(main_id);
            duration += graph.edges()[static_cast<std::size_t>(main_id)].data.time;
            cost += graph.edges()[static_cast<std::size_t>(main_id)].data.cost;
        }
        const RecoveredSolution recovered =
            recover_selected_edge_solution(data, graph, active, cost, 0.0);
        if (static_cast<int>(recovered.routes.size()) != vehicle_count) {
            throw std::runtime_error("The recombined solution has a different number of routes.");
        }
        std::size_t service_nodes = 0;
        for (const RecoveredRouteSolution& route : recovered.routes) {
            if (route.total_time > data.time_limit + 1e-6) {
                throw std::runtime_error("A recombined route violates the duration limit.");
            }
        }
        for (int main_id : active) {
            if (graph.edges()[static_cast<std::size_t>(main_id)].v != end_node_id) ++service_nodes;
        }
        if (service_nodes != 2U * data.requests.size()) {
            throw std::runtime_error("The recombined solution does not visit every service node once.");
        }
        if (std::fabs(cost - value) > 1e-4) {
            throw std::runtime_error("The recombined solution cost differs from the model value.");
        }
        result.improved = true;
        result.active_edge_ids = std::move(active);
        result.total_duration = duration;
        result.total_original_cost = cost;
        return finish(status_name);
    } catch (const GRBException& error) {
        result.improved = false;
        result.active_edge_ids.clear();
        log << " error=\"gurobi " << error.getErrorCode() << ": " << error.getMessage() << "\"";
        return finish("error");
    } catch (const std::exception& error) {
        result.improved = false;
        result.active_edge_ids.clear();
        log << " error=\"" << error.what() << "\"";
        return finish("error");
    }
}

}  // namespace spdp
