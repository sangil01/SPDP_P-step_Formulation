#include "ArcElimination.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "PstepValidInequality.h"
#include "gurobi_c++.h"

namespace spdp {
namespace {

using Clock = std::chrono::steady_clock;

// Costs are integers, so "cost <= U - 1" is tested as "bound > U - 1 + margin".
constexpr double kBoundMargin = 1e-4;


double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

bool costs_are_integral(const MultiDiGraph& graph) {
    for (const EdgeRecord& edge : graph.edges()) {
        if (std::abs(edge.data.cost - std::round(edge.data.cost)) > 1e-9) {
            return false;
        }
    }
    return true;
}

State canonical_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

bool no_solution_status(int status) {
    // Used for models with binary edge variables and bounded time variables,
    // which are never unbounded.
    return status == GRB_CUTOFF || status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD;
}

struct SubGraph {
    MultiDiGraph graph;
    std::vector<int> local_to_main;
};

SubGraph build_subgraph(const MultiDiGraph& graph, const std::vector<char>& keep) {
    SubGraph sub;
    for (const NodeSpec& node : graph.nodes()) {
        sub.graph.add_node(node);
    }
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        if (!keep[edge_id]) {
            continue;
        }
        const EdgeRecord& edge = graph.edges()[edge_id];
        sub.local_to_main.push_back(static_cast<int>(edge_id));
        sub.graph.add_edge(edge.u, edge.v, edge.data);
    }
    return sub;
}

std::vector<double> to_main_values(
    const SubGraph& sub,
    const std::vector<double>& local_values,
    std::size_t edge_count
) {
    std::vector<double> values(edge_count, 0.0);
    for (std::size_t local = 0; local < sub.local_to_main.size() && local < local_values.size();
         ++local) {
        values[static_cast<std::size_t>(sub.local_to_main[local])] = local_values[local];
    }
    return values;
}

// Merged time-free model (Step 2). Service nodes with the
// same kind, location and container type (pickups: also the same treatment
// facility, which is part of the full skip F(h, l) they load; a delivery only
// unloads an empty skip E(h)) form a group g of n_g nodes.
// The E1 edges with the same class and the same groups at both ends become one
// integer variable x_a in [0, number of such edges] with the smallest of their
// costs. Rows: in- and out-degree n_g per group, state flow per (group, state),
// the depot balance, and every valid inequality row whose coefficient is the
// same on all edges of each merged variable (other rows and the capacity cuts are
// left out). Every M0 solution maps to a merged solution with the same classes
// and at most the same cost, so the merged model is a relaxation of M0: "no
// solution of cost <= U - 1" still proves U optimal, and the classes of the
// merged solutions of cost <= U - 1 cover those of the M0 solutions.
class MergedM0 {
public:
    struct Outcome {
        int status = 0;
        bool has_solution = false;
        double value = -1.0;
        double bound = -1.0;
        double nodes = -1.0;
        double seconds = 0.0;
        std::vector<int> classes;  // classes of the solution returned (x_a > 0.5)
    };
    // Called for every solution with its classes. Returns 0 (accept), 1 (add
    // the lazy row "at least one class outside E2") or 2 (add it and stop).
    using SolutionHandler = std::function<int(const std::vector<int>&)>;

    MergedM0(
        const SPDPData& data,
        const MultiDiGraph& graph,
        const std::vector<char>& keep,
        const std::vector<int>& class_of,
        const DirectTwoIndexOptions& options
    )
        : graph_(graph), class_of_(class_of), options_(options) {
        const NodeId end = graph.end_node_id();
        std::map<std::tuple<int, int, int, int>, int> group_ids;
        group_of_.assign(static_cast<std::size_t>(end) + 1, -1);
        for (NodeId v = 1; v < end; ++v) {
            if (!graph.is_physical_service_node(v)) {
                continue;
            }
            const NodeSpec& node = graph.node(v);
            const int landfill = node.kind == NodeSpec::Kind::Pickup
                                     ? node.landfill_location.value_or(-1)
                                     : -1;
            const auto key = std::make_tuple(
                static_cast<int>(node.kind), node.location, node.container_type.value_or(-1),
                landfill);
            const auto inserted = group_ids.emplace(key, static_cast<int>(group_ids.size()));
            const int g = inserted.first->second;
            group_of_[static_cast<std::size_t>(v)] = g;
            if (inserted.second) {
                group_size_.push_back(0);
                group_kind_.push_back(node.kind);
            }
            ++group_size_[static_cast<std::size_t>(g)];
            ++services_;
        }
        std::map<std::tuple<int, int, int>, int> var_ids;
        member_of_.assign(graph.number_of_edges(), -1);
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            if (!keep[edge_id] || (edge.u == 0 && edge.v == end)) {
                continue;
            }
            ++e1_edges_;
            const int tail = edge.u == 0 ? -1 : group_of_[static_cast<std::size_t>(edge.u)];
            const int head = edge.v == end ? -2 : group_of_[static_cast<std::size_t>(edge.v)];
            const auto inserted = var_ids.emplace(
                std::make_tuple(class_of[edge_id], tail, head), static_cast<int>(vars_.size()));
            if (inserted.second) {
                Var var;
                var.tail = tail;
                var.head = head;
                var.cls = class_of[edge_id];
                var.cost = edge.data.cost;
                var.start = canonical_state(edge.data.start_state);
                var.end = canonical_state(edge.data.end_state);
                vars_.push_back(var);
            }
            Var& var = vars_[static_cast<std::size_t>(inserted.first->second)];
            if (edge.data.cost != var.cost) {
                ++mixed_cost_;
                var.cost = std::min(var.cost, edge.data.cost);
            }
            ++var.count;
            member_of_[edge_id] = inserted.first->second;
        }
        // Valid inequality rows with one coefficient per merged variable.
        DirectTwoIndexOptions vi_options = options;
        const std::vector<PstepValidInequalityRow> rows =
            build_pstep_valid_inequality_rows(data, graph, vi_options.valid_inequalities);
        for (const PstepValidInequalityRow& row : rows) {
            std::map<int, std::pair<double, int>> seen;  // var -> (coefficient, terms)
            bool uniform = true;
            for (const auto& [edge_id, coefficient] : row.edge_terms) {
                const int a = member_of_[edge_id];
                if (a < 0) {
                    continue;  // edge not in E1 (fixed to 0) or the dummy edge
                }
                auto found = seen.find(a);
                if (found == seen.end()) {
                    seen.emplace(a, std::make_pair(coefficient, 1));
                } else if (found->second.first != coefficient) {
                    uniform = false;
                    break;
                } else {
                    ++found->second.second;
                }
            }
            for (const auto& [a, item] : seen) {
                if (item.second != vars_[static_cast<std::size_t>(a)].count) {
                    uniform = false;  // some edges of the variable are not in the row
                }
            }
            if (!uniform) {
                ++rows_dropped_;
                continue;
            }
            MergedRow merged;
            merged.sense = row.sense;
            merged.rhs = row.rhs;
            for (const auto& [a, item] : seen) {
                merged.terms.emplace_back(a, item.first);
            }
            rows_.push_back(std::move(merged));
        }
    }

    std::string describe() const {
        std::ostringstream text;
        text << "groups=" << group_size_.size() << " services=" << services_
             << " vars=" << vars_.size() << " e1_edges=" << e1_edges_
             << " mixed_cost_edges=" << mixed_cost_ << " vi_rows_kept=" << rows_.size()
             << " vi_rows_dropped=" << rows_dropped_;
        return text.str();
    }

    // Root LP value of the merged model (diagnostic).
    double root_lp(double time_limit) {
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, 0);
        env.start();
        Model built = build(env, -1.0, time_limit);
        GRBModel relaxed = built.model->relax();
        relaxed.optimize();
        return relaxed.get(GRB_IntAttr_Status) == GRB_OPTIMAL
            ? relaxed.get(GRB_DoubleAttr_ObjVal) : -1.0;
    }

    // first_only: stop at the first solution. Otherwise every solution goes to
    // handler (lazy constraints) until the search ends or handler returns 2.
    Outcome solve(
        double cutoff,
        double time_limit,
        bool first_only,
        const std::vector<char>* in_e2,
        const SolutionHandler& handler,
        const std::string& log_path
    ) {
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, log_path.empty() ? 0 : 1);
        env.set(GRB_IntParam_LogToConsole, 0);
        if (!log_path.empty()) {
            env.set(GRB_StringParam_LogFile, log_path);
        }
        env.start();
        Model built = build(env, cutoff, time_limit);
        GRBModel& model = *built.model;
        if (!first_only && in_e2 != nullptr) {
            GRBLinExpr outside = 0.0;
            for (std::size_t a = 0; a < vars_.size(); ++a) {
                if (vars_[a].count > 0 && !(*in_e2)[static_cast<std::size_t>(vars_[a].cls)]) {
                    outside += built.x[a];
                }
            }
            model.addConstr(outside >= 1.0);
        }
        Callback callback(*this, built.x, first_only, in_e2, handler);
        if (!first_only) {
            model.set(GRB_IntParam_LazyConstraints, 1);
        }
        model.setCallback(&callback);
        model.optimize();
        Outcome out;
        out.status = model.get(GRB_IntAttr_Status);
        out.nodes = model.get(GRB_DoubleAttr_NodeCount);
        out.seconds = model.get(GRB_DoubleAttr_Runtime);
        out.has_solution = model.get(GRB_IntAttr_SolCount) > 0;
        try {
            out.bound = model.get(GRB_DoubleAttr_ObjBound);
        } catch (const GRBException&) {
            out.bound = -1.0;
        }
        if (out.has_solution) {
            out.value = model.get(GRB_DoubleAttr_ObjVal);
            std::set<int> classes;
            for (std::size_t a = 0; a < vars_.size(); ++a) {
                if (built.x[a].get(GRB_DoubleAttr_X) > 0.5) {
                    classes.insert(vars_[a].cls);
                }
            }
            out.classes.assign(classes.begin(), classes.end());
        }
        return out;
    }

private:
    struct Var {
        int tail = -1;  // group, -1 the start node
        int head = -2;  // group, -2 the end node
        int cls = -1;
        double cost = 0.0;
        int count = 0;
        State start{};
        State end{};
    };
    struct MergedRow {
        PstepValidInequalitySense sense = PstepValidInequalitySense::GreaterEqual;
        double rhs = 0.0;
        std::vector<std::pair<int, double>> terms;
    };
    struct Model {
        std::unique_ptr<GRBModel> model;
        std::vector<GRBVar> x;
    };

    Model build(GRBEnv& env, double cutoff, double time_limit) {
        Model built;
        built.model = std::make_unique<GRBModel>(env);
        GRBModel& model = *built.model;
        if (options_.gurobi_threads >= 0) {
            model.set(GRB_IntParam_Threads, options_.gurobi_threads);
        }
        model.set(GRB_DoubleParam_FeasibilityTol, detail_tolerance());
        model.set(GRB_DoubleParam_OptimalityTol, detail_tolerance());
        if (time_limit > 0.0) {
            model.set(GRB_DoubleParam_TimeLimit, time_limit);
        }
        if (cutoff >= 0.0) {
            model.set(GRB_DoubleParam_Cutoff, cutoff);
        }
        built.x.reserve(vars_.size());
        for (const Var& var : vars_) {
            built.x.push_back(model.addVar(
                0.0, static_cast<double>(var.count), var.cost, GRB_INTEGER));
        }
        const std::size_t groups = group_size_.size();
        std::vector<GRBLinExpr> in(groups);
        std::vector<GRBLinExpr> out(groups);
        std::map<std::pair<int, State>, GRBLinExpr> flow;
        std::map<std::pair<int, State>, GRBLinExpr> by_state;
        GRBLinExpr departures = 0.0;
        GRBLinExpr returns = 0.0;
        for (std::size_t a = 0; a < vars_.size(); ++a) {
            const Var& var = vars_[a];
            const GRBVar& x = built.x[a];
            if (var.head >= 0) {
                in[static_cast<std::size_t>(var.head)] += x;
                flow[{var.head, var.end}] += x;
                by_state[{var.head, var.end}] += x;
            }
            if (var.tail >= 0) {
                out[static_cast<std::size_t>(var.tail)] += x;
                flow[{var.tail, var.start}] -= x;
            }
            if (var.tail == -1 && var.head >= 0 &&
                group_kind_[static_cast<std::size_t>(var.head)] == NodeSpec::Kind::Pickup) {
                departures += x;
            }
            if (var.head == -2 && var.tail >= 0 &&
                group_kind_[static_cast<std::size_t>(var.tail)] == NodeSpec::Kind::Delivery) {
                returns += x;
            }
        }
        for (std::size_t g = 0; g < groups; ++g) {
            model.addConstr(in[g] == static_cast<double>(group_size_[g]));
            model.addConstr(out[g] == static_cast<double>(group_size_[g]));
        }
        for (auto& [key, expression] : flow) {
            model.addConstr(expression == 0.0);
        }
        model.addConstr(departures == returns);
        for (const MergedRow& row : rows_) {
            GRBLinExpr expression = 0.0;
            for (const auto& [a, coefficient] : row.terms) {
                expression += coefficient * built.x[static_cast<std::size_t>(a)];
            }
            if (row.sense == PstepValidInequalitySense::GreaterEqual) {
                model.addConstr(expression >= row.rhs);
            } else if (row.sense == PstepValidInequalitySense::LessEqual) {
                model.addConstr(expression <= row.rhs);
            } else {
                model.addConstr(expression == row.rhs);
            }
        }
        // Number of nodes of each group in each state, branched on first.
        for (auto& [key, expression] : by_state) {
            GRBVar w = model.addVar(
                0.0, static_cast<double>(group_size_[static_cast<std::size_t>(key.first)]),
                0.0, GRB_INTEGER);
            w.set(GRB_IntAttr_BranchPriority, 10);
            model.addConstr(w == expression);
        }
        model.update();
        return built;
    }

    static double detail_tolerance() { return 1e-6; }

    class Callback : public GRBCallback {
    public:
        Callback(
            MergedM0& owner,
            std::vector<GRBVar>& x,
            bool first_only,
            const std::vector<char>* in_e2,
            const SolutionHandler& handler
        )
            : owner_(owner), x_(x), first_only_(first_only), in_e2_(in_e2), handler_(handler) {}

    protected:
        void callback() override {
            if (where != GRB_CB_MIPSOL) {
                return;
            }
            if (first_only_) {
                abort();
                return;
            }
            std::unique_ptr<double[]> values(getSolution(x_.data(), static_cast<int>(x_.size())));
            std::set<int> classes;
            for (std::size_t a = 0; a < x_.size(); ++a) {
                if (values[a] > 0.5) {
                    classes.insert(owner_.vars_[a].cls);
                }
            }
            const std::vector<int> list(classes.begin(), classes.end());
            const int action = handler_(list);
            if (action == 0) {
                return;
            }
            GRBLinExpr outside = 0.0;
            for (std::size_t a = 0; a < x_.size(); ++a) {
                const auto& var = owner_.vars_[a];
                if (var.count > 0 && !(*in_e2_)[static_cast<std::size_t>(var.cls)]) {
                    outside += x_[a];
                }
            }
            addLazy(outside >= 1.0);
            if (action == 2) {
                abort();
            }
        }

    private:
        MergedM0& owner_;
        std::vector<GRBVar>& x_;
        bool first_only_;
        const std::vector<char>* in_e2_;
        const SolutionHandler& handler_;
    };

    const MultiDiGraph& graph_;
    const std::vector<int>& class_of_;
    DirectTwoIndexOptions options_;
    std::vector<int> group_of_;
    std::vector<int> group_size_;
    std::vector<NodeSpec::Kind> group_kind_;
    std::vector<int> member_of_;
    std::vector<Var> vars_;
    std::vector<MergedRow> rows_;
    int services_ = 0;
    int e1_edges_ = 0;
    int mixed_cost_ = 0;
    int rows_dropped_ = 0;
};

// Class of an edge: the same move for different requests.
std::vector<int> edge_classes(const MultiDiGraph& graph, int& class_count) {
    using Key = std::tuple<int, int, int, int, int, State, State, std::vector<int>>;
    std::map<Key, int> ids;
    std::vector<int> class_of(graph.number_of_edges(), -1);
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        const NodeSpec& tail = graph.node(edge.u);
        const NodeSpec& head = graph.node(edge.v);
        Key key(
            static_cast<int>(tail.kind), static_cast<int>(head.kind), tail.location,
            head.location, head.container_type.value_or(-1), edge.data.start_state,
            edge.data.end_state, edge.data.sequence_pi);
        const auto inserted = ids.emplace(std::move(key), static_cast<int>(ids.size()));
        class_of[edge_id] = inserted.first->second;
    }
    class_count = static_cast<int>(ids.size());
    return class_of;
}

double cost_of(const MultiDiGraph& graph, const std::vector<double>& edge_values) {
    double cost = 0.0;
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        if (edge_values[edge_id] > 0.5) {
            cost += graph.edges()[edge_id].data.cost;
        }
    }
    return cost;
}

// The best known solution reported as the answer of the main model.
DirectTwoIndexResult solution_as_result(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<double>& edge_values,
    bool proved_optimal,
    double lower_bound,
    double runtime_seconds
) {
    DirectTwoIndexResult result;
    result.model_type = DirectTwoIndexModelType::IP;
    result.status = proved_optimal ? GRB_OPTIMAL : GRB_TIME_LIMIT;
    result.solved_to_optimality = proved_optimal;
    result.hit_time_limit = !proved_optimal;
    result.has_feasible_solution = true;
    result.runtime_seconds = runtime_seconds;
    result.variable_count = static_cast<int>(graph.number_of_edges());
    result.edge_values.assign(graph.number_of_edges(), 0.0);
    result.total_duration = 0.0;
    result.total_original_cost = 0.0;
    result.departure_flow = 0.0;
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        if (edge_values[edge_id] <= 0.5) {
            continue;
        }
        result.edge_values[edge_id] = 1.0;
        const EdgeRecord& edge = graph.edges()[edge_id];
        result.total_duration += edge.data.time;
        result.total_original_cost += edge.data.cost;
        if (edge.u == 0 && graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            result.departure_flow += 1.0;
            ++result.vehicle_count;
        }
    }
    result.total_duration_plus_fixed =
        result.total_duration + data.fixed_vehicle_cost * result.departure_flow;
    result.total_travel_cost =
        result.total_original_cost - data.fixed_vehicle_cost * result.departure_flow;
    result.objective_value = result.total_original_cost;
    if (proved_optimal) {
        result.objective_bound = result.objective_value;
        result.has_certified_bound = true;
        result.gap_percent = 0.0;
    } else if (lower_bound >= 0.0) {
        result.objective_bound = std::min(lower_bound, result.objective_value);
        result.has_certified_bound = true;
        result.gap_percent = result.objective_value > 0.0
            ? 100.0 * (result.objective_value - result.objective_bound) / result.objective_value
            : 0.0;
    }
    return result;
}

// Options of a solve on a subgraph: the main options without the parts that
// refer to edge ids of the full graph.
DirectTwoIndexOptions subgraph_options(
    const DirectTwoIndexOptions& main_options,
    const std::string& root_lp_log_path,
    const std::string& mip_log_path
) {
    DirectTwoIndexOptions options = main_options;
    options.model_type = DirectTwoIndexModelType::IP;
    options.initial_edge_start.clear();
    options.type_cover = TypeTravelTimeCoverOptions{};
    options.type_cover_rows.clear();
    options.mip_start_check_only = false;
    options.valid_inequalities.log_stream = nullptr;
    options.gurobi_log_path = mip_log_path;
    options.pre_mip_lp_gurobi_log_path = root_lp_log_path;
    return options;
}

// The merged model writes its own logs (MergedM0::solve), so none here.
DirectTwoIndexOptions time_free_options(const DirectTwoIndexOptions& main_options) {
    DirectTwoIndexOptions options = subgraph_options(main_options, "", "");
    options.add_time_constraints = false;
    options.add_time_flow_formulation = false;
    options.time_flow_state_disaggregated = false;
    return options;
}

bool parse_number(const std::string& value, double& out) {
    try {
        std::size_t used = 0;
        out = std::stod(value, &used);
        return used == value.size();
    } catch (const std::exception&) {
        return false;
    }
}

}  // namespace


bool parse_arc_elimination_config(
    const std::string& text,
    ArcEliminationOptions& options,
    std::string& error
) {
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty()) {
            continue;
        }
        const std::size_t eq = item.find('=');
        if (eq == std::string::npos) {
            error = "'" + item + "' is not key=value";
            return false;
        }
        const std::string key = item.substr(0, eq);
        double number = 0.0;
        if (!parse_number(item.substr(eq + 1), number)) {
            error = "bad value for '" + key + "'";
            return false;
        }
        const bool whole = std::abs(number - std::round(number)) < 1e-9;
        const int integer = static_cast<int>(std::lround(number));
        bool ok = true;
        if (key == "rc") {
            ok = whole && integer >= 0 && integer <= 2;
            options.rc_mode = integer;
        } else if (key == "rc_tl") {
            ok = number >= 0.0;
            options.rc_time_limit = number;
        } else if (key == "total_tl") {
            ok = number >= 0.0;
            options.total_time_limit = number;
        } else {
            error = "unknown key '" + key + "'";
            return false;
        }
        if (!ok) {
            error = "bad value for '" + key + "'";
            return false;
        }
    }
    return true;
}

ArcEliminationResult run_arc_elimination(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& main_options,
    const std::vector<double>& incumbent_edge_values,
    double incumbent_cost,
    const ArcEliminationOptions& options,
    std::ostream* log
) {
    ArcEliminationResult out;
    const Clock::time_point start = Clock::now();
    const auto say = [&](const std::ostringstream& line) {
        if (log != nullptr) {
            *log << line.str() << '\n';
            log->flush();
        }
    };
    const auto skip = [&](const char* reason) {
        std::ostringstream line;
        line << "[arc-elim] skipped reason=" << reason;
        say(line);
        return out;
    };
    if (main_options.model_type != DirectTwoIndexModelType::IP ||
        main_options.objective != DirectTwoIndexObjective::OriginalCost) {
        return skip("needs-ip-with-original-cost");
    }
    const std::size_t edge_count = graph.number_of_edges();
    if (incumbent_edge_values.size() != edge_count) {
        return skip("no-incumbent");
    }
    if (!costs_are_integral(graph) ||
        std::abs(incumbent_cost - std::round(incumbent_cost)) > 1e-6 ||
        std::abs(cost_of(graph, incumbent_edge_values) - incumbent_cost) > 1e-6) {
        return skip("non-integer-costs-or-inconsistent-incumbent");
    }
    out.ran = true;

    // Budget of the whole procedure.
    const double total_limit = main_options.solver_time_limit;
    const auto remaining = [&]() {
        return total_limit > 0.0
            ? std::max(0.0, total_limit - seconds_since(start))
            : 1e18;
    };
    const auto limited = [&](double configured) {
        const double left = remaining();
        const double value = configured > 0.0 ? std::min(configured, left) : left;
        return std::max(value, 1e-3);
    };
    // Time left for Steps 1 and 2 under total_time_limit.
    const auto total_left = [&]() {
        return options.total_time_limit > 0.0
            ? options.total_time_limit - options.elapsed_before - seconds_since(start)
            : 1e18;
    };

    std::vector<double> incumbent = incumbent_edge_values;
    double upper = std::round(incumbent_cost);
    double lower = -1.0;
    {
        std::ostringstream line;
        line << "[arc-elim] start incumbent=" << upper << " edges=" << edge_count
             << " rc=" << options.rc_mode << " rc_tl=" << options.rc_time_limit
             << " total_tl=" << options.total_time_limit
             << " elapsed_before=" << options.elapsed_before
             << " budget=" << remaining();
        say(line);
    }
    const auto finish = [&](bool proved, const char* reason) {
        out.finished = true;
        out.proved_optimal = proved;
        out.incumbent_cost = upper;
        out.lower_bound = proved ? upper : lower;
        out.seconds = seconds_since(start);
        out.result = solution_as_result(data, graph, incumbent, proved, lower, out.seconds);
        if (std::abs(out.result.total_original_cost - upper) > 1e-6) {
            throw std::runtime_error(
                "Arc elimination: the reported solution does not have the reported cost.");
        }
        std::ostringstream line;
        line << "[arc-elim] result=" << reason << " proved=" << (proved ? 1 : 0)
             << " value=" << upper << " lower_bound=" << out.lower_bound
             << " seconds=" << out.seconds;
        say(line);
        return out;
    };

    // The start-to-end dummy edge is part of every model and is never deleted.
    std::vector<char> structural(edge_count, 0);
    for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        structural[edge_id] = edge.u == 0 && edge.v == graph.end_node_id() ? 1 : 0;
    }

    // ---- Step 1: root LP of the main model, bound check and (rc > 0) deletion by the dual bound.
    std::vector<char> keep1(edge_count, 1);
    bool have_lp = false;
    double lp_bound = 0.0;
    std::vector<double> gain;
    // E1: the edges kept by the dual bound.
    const auto refilter = [&]() {
        if (!have_lp) {
            return;
        }
        const double limit = upper - 1.0 + kBoundMargin;
        std::size_t kept = 0;
        std::size_t incumbent_edges = 0;
        std::size_t incumbent_kept = 0;
        for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
            keep1[edge_id] = options.rc_mode == 0 || structural[edge_id] ||
                                     !(lp_bound + gain[edge_id] > limit)
                                 ? 1
                                 : 0;
            kept += keep1[edge_id] ? 1 : 0;
            if (incumbent[edge_id] > 0.5) {
                ++incumbent_edges;
                incumbent_kept += keep1[edge_id] ? 1 : 0;
            }
        }
        std::ostringstream line;
        line << "[arc-elim-rc] incumbent=" << upper << " kept=" << kept
             << " of=" << edge_count
             << " incumbent_edges_kept=" << incumbent_kept << "/" << incumbent_edges;
        say(line);
    };

    bool limit_hit = false;
    if (options.total_time_limit > 0.0 && total_left() < 120.0) {
        limit_hit = true;
        std::ostringstream line;
        line << "[arc-elim-limit] before=step1 total_left=" << total_left();
        say(line);
    }
    if (!limit_hit) {
        const Clock::time_point lp_start = Clock::now();
        DirectTwoIndexOptions lp = main_options;
        lp.model_type = DirectTwoIndexModelType::LP;
        lp.initial_edge_start.clear();
        lp.mip_start_check_only = false;
        lp.valid_inequalities.log_stream = nullptr;
        lp.mip_method = 2;
        lp.mip_node_method = -1;
        lp.mip_crossover = options.rc_mode == 1 ? -1 : 0;
        lp.collect_lp_dual_bound = true;
        lp.solver_time_limit = limited(std::min(
            options.rc_time_limit > 0.0 ? options.rc_time_limit : 1e18, total_left()));
        lp.gurobi_log_path = p0_gurobi_log_path(options.log_prefix, kLogE0RootLp);
        const DirectTwoIndexResult lp_result = solve_direct_two_index_model(data, graph, lp);
        std::size_t deleted_by_gurobi_values = 0;
        if (lp_result.lp_dual_bound_valid) {
            have_lp = true;
            lp_bound = lp_result.lp_dual_bound;
            gain = lp_result.edge_fix_gain;
            lower = std::max(lower, std::ceil(lp_bound - 1e-6));
            for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
                if (lp_result.objective_value + lp_result.edge_reduced_costs[edge_id] >
                    upper - 1.0 + kBoundMargin) {
                    ++deleted_by_gurobi_values;
                }
            }
        }
        std::ostringstream line;
        line << "[arc-elim-lp] status=" << lp_result.status
             << " valid=" << (have_lp ? 1 : 0)
             << " lp_value=" << (lp_result.has_feasible_solution ? lp_result.objective_value : -1.0)
             << " dual_bound=" << (have_lp ? lp_bound : -1.0)
             << " clipped_duals=" << lp_result.lp_clipped_dual_count
             << " cut_rounds=" << lp_result.lp_cut_rounds
             << " blossom_cuts=" << lp_result.capacity_blossom_stats.management.cuts_added
             << " deleted_with_gurobi_values=" << deleted_by_gurobi_values
             << " seconds=" << seconds_since(lp_start);
        say(line);
    }
    refilter();
    if (have_lp && lower >= upper) {
        return finish(true, "incumbent-optimal-by-lp-bound");
    }

    // ---- Step 2: classes of the merged time-free solutions of cost <= U - 1.
    // E2 starts with the classes of the incumbent edges.
    int class_count = 0;
    const std::vector<int> class_of = edge_classes(graph, class_count);
    std::vector<char> in_e2(static_cast<std::size_t>(class_count), 0);
    for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
        if (structural[edge_id] || incumbent[edge_id] > 0.5) {
            in_e2[static_cast<std::size_t>(class_of[edge_id])] = 1;
        }
    }
    const auto e2_keep = [&]() {
        std::vector<char> keep(edge_count, 0);
        for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
            keep[edge_id] =
                keep1[edge_id] && in_e2[static_cast<std::size_t>(class_of[edge_id])] ? 1 : 0;
        }
        return keep;
    };
    const auto count = [](const std::vector<char>& keep) {
        return static_cast<std::size_t>(std::count(keep.begin(), keep.end(), 1));
    };

    const Clock::time_point collect_start = Clock::now();
    if (limit_hit) {
        std::ostringstream line;
        line << "[arc-elim-skip] incumbent=" << upper << " reason=total_tl";
        say(line);
    }
    // Step 2 ends at the total_tl deadline (or with the budget when total_tl = 0).
    const auto collect_left = [&]() {
        return limit_hit ? 0.0 : std::min(remaining(), total_left());
    };
    bool verified = false;
    int m0_solutions = 0;
    double last_solution_seconds = -1.0;  // Step 2b: time of the last solution found

    // Merged time-free model: built once on E1.
    std::unique_ptr<MergedM0> merged;
    if (collect_left() > 1.0) {
        const Clock::time_point merged_start = Clock::now();
        const DirectTwoIndexOptions merged_options =
            time_free_options(main_options);
        merged = std::make_unique<MergedM0>(data, graph, keep1, class_of, merged_options);
        const double merged_lp = merged->root_lp(std::max(1.0, std::min(60.0, collect_left())));
        std::ostringstream line;
        line << "[arc-elim-merged] " << merged->describe() << " root_lp=" << merged_lp
             << " incumbent=" << upper << " seconds=" << seconds_since(merged_start);
        say(line);
    }
    const auto merged_none = [](int status) {
        return status == GRB_CUTOFF || status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD;
    };

    // Step 2a: merged model with "cost <= U - 1" only, stopped at the first solution.
    if (merged && collect_left() > 1.0) {
        const Clock::time_point first_start = Clock::now();
        const MergedM0::Outcome first = merged->solve(
            upper - 0.5, std::max(collect_left(), 1e-3), true, nullptr, {},
            p0_gurobi_log_path(options.log_prefix, kLogE1InfeasibleCheck));
        const bool none = merged_none(first.status);
        {
            std::ostringstream line;
            line << "[arc-elim-first] kind=merged status=" << first.status
                 << " no_solution_below_incumbent=" << (none ? 1 : 0)
                 << " value=" << first.value << " incumbent=" << upper
                 << " nodes=" << first.nodes
                 << " seconds=" << seconds_since(first_start);
            say(line);
        }
        if (none) {
            return finish(true, "incumbent-optimal");
        }
        if (first.has_solution) {
            ++m0_solutions;
            for (const int cls : first.classes) {
                if (!in_e2[static_cast<std::size_t>(cls)]) {
                    in_e2[static_cast<std::size_t>(cls)] = 1;
                }
            }
        }
    }

    // Step 2b: one search that collects the classes.
    const auto e1_outside_e2 = [&]() {
        for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
            if (keep1[edge_id] && !structural[edge_id] &&
                !in_e2[static_cast<std::size_t>(class_of[edge_id])]) {
                return true;
            }
        }
        return false;
    };
    if (collect_left() > 1.0 && merged) {
        if (!e1_outside_e2()) {
            verified = true;  // E2 = E1
        } else {
            bool everything = false;
            const Clock::time_point search_start = Clock::now();
            const MergedM0::SolutionHandler handler = [&](const std::vector<int>& classes) {
                ++m0_solutions;
                last_solution_seconds = seconds_since(search_start);
                for (const int cls : classes) {
                    in_e2[static_cast<std::size_t>(cls)] = 1;
                }
                if (!e1_outside_e2()) {
                    everything = true;
                    return 2;
                }
                return 1;
            };
            const Clock::time_point round_start = Clock::now();
            const MergedM0::Outcome search = merged->solve(
                upper - 0.5, std::max(collect_left(), 1e-3), false, &in_e2, handler,
                p0_gurobi_log_path(options.log_prefix, kLogE1FeasibleArcCollect));
            const bool none = merged_none(search.status);
            {
                std::ostringstream line;
                line << "[arc-elim-round] kind=merged status=" << search.status
                     << " complete=" << ((none || everything) ? 1 : 0)
                     << " incumbent=" << upper
                     << " e2_edges=" << count(e2_keep())
                     << " m0_solutions=" << m0_solutions
                     << " last_solution_seconds=" << last_solution_seconds
                     << " mip_seconds=" << search.seconds
                     << " nodes=" << search.nodes
                     << " seconds=" << seconds_since(round_start);
                say(line);
            }
            if (none || everything) {
                verified = true;
            }
        }
    }
    {
        std::ostringstream line;
        line << "[arc-elim-collect] complete=" << (verified ? 1 : 0) << " rounds=1"
             << " m0_solutions=" << m0_solutions
             << " e1_edges=" << count(keep1) << " e2_edges=" << count(e2_keep())
             << " edges=" << edge_count
             << " seconds=" << seconds_since(collect_start);
        say(line);
    }
    // ---- Step 3 (E2 complete): main model on E2 with "cost <= U - 1".
    // ---- Fallback (E2 not complete): main model on E1, started from the known
    //      solution when all its edges are in E1, with "cost <= U - 1" otherwise.
    const std::vector<char> keep_final = verified ? e2_keep() : keep1;
    const SubGraph sub = build_subgraph(graph, keep_final);
    DirectTwoIndexOptions final_options = subgraph_options(
        main_options,
        p0_gurobi_log_path(options.log_prefix, verified ? kLogE2RootLp : kLogE1RootLp),
        p0_gurobi_log_path(options.log_prefix, verified ? kLogE2Milp : kLogE1Milp));
    bool started = false;
    if (!verified) {
        std::vector<double> start_values(sub.local_to_main.size(), 0.0);
        std::size_t present = 0;
        std::size_t wanted = 0;
        for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
            wanted += incumbent[edge_id] > 0.5 ? 1 : 0;
        }
        for (std::size_t local = 0; local < sub.local_to_main.size(); ++local) {
            if (incumbent[static_cast<std::size_t>(sub.local_to_main[local])] > 0.5) {
                start_values[local] = 1.0;
                ++present;
            }
        }
        if (present == wanted) {
            final_options.initial_edge_start = start_values;
            started = true;
        }
    }
    if (!started) {
        final_options.cutoff = upper - 0.5;
    }
    final_options.solver_time_limit = limited(0.0);
    final_options.pre_mip_root_lp_time_limit = limited(main_options.pre_mip_root_lp_time_limit);
    const Clock::time_point final_start = Clock::now();
    const DirectTwoIndexResult final_result =
        solve_direct_two_index_model(data, sub.graph, final_options);
    const bool none = no_solution_status(final_result.status);
    const double before = upper;
    if (final_result.has_feasible_solution) {
        const std::vector<double> values = to_main_values(sub, final_result.edge_values, edge_count);
        const double value = std::round(cost_of(graph, values));
        if (value < upper - 0.5) {
            upper = value;
            incumbent = values;
        }
    }
    const bool proved = none || final_result.solved_to_optimality;
    if (!proved && final_result.has_certified_bound) {
        // The optimum is min(U, optimum on the final edge set).
        lower = std::max(
            lower, std::min(before, std::ceil(final_result.objective_bound - 1e-6)));
    }
    {
        std::ostringstream line;
        line << "[arc-elim-final] path=" << (verified ? "e2" : "fallback-e1")
             << " start=" << (started ? 1 : 0)
             << " status=" << final_result.status
             << " edges=" << sub.local_to_main.size()
             << " incumbent_before=" << before << " incumbent_after=" << upper
             << " bound=" << (final_result.has_certified_bound ? final_result.objective_bound : -1.0)
             << " pre_mip_seconds=" << final_result.pre_mip_seconds
             << " blossom_cuts=" << final_result.capacity_blossom_stats.management.cuts_added
             << " mip_seconds=" << final_result.runtime_seconds
             << " seconds=" << seconds_since(final_start);
        say(line);
    }
    ArcEliminationResult done = finish(
        proved,
        proved ? (upper < before ? "better-solution-optimal" : "incumbent-optimal")
               : "not-proved");
    done.result.capacity_blossom_stats = final_result.capacity_blossom_stats;
    done.result.constraint_count = final_result.constraint_count;
    done.result.valid_inequality_count = final_result.valid_inequality_count;
    done.result.pre_mip_rounds = final_result.pre_mip_rounds;
    done.result.pre_mip_rows = final_result.pre_mip_rows;
    done.result.pre_mip_seconds = final_result.pre_mip_seconds;
    return done;
}

}  // namespace spdp
